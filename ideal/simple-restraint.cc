/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// #define ANALYSE_REFINEMENT_TIMING

#include <string.h> // for strcmp

#ifdef ANALYSE_REFINEMENT_TIMING
#include <sys/time.h>
#endif // ANALYSE_REFINEMENT_TIMING

// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL


#include <fstream>
#include <algorithm> // for sort
#include <stdexcept>
#include <iomanip>

#ifdef HAVE_CXX_THREAD
#include <thread>
#include <chrono>
#endif // HAVE_CXX_THREAD

#include "geometry/main-chain.hh"
#include "simple-restraint.hh"

//
#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

// #include "mmdb.h" // for printing of mmdb::Atom pointers as info not raw
                     // pointers.  Removed. Too much (linking issues in)
                     // Makefile pain.

#include "compat/coot-sysdep.h"



// iend_res is inclusive, so that 17,17 selects just residue 17.
//   have_disulfide_residues: other residues are included in the
//				residues_mol for disphide restraints.
// 
coot::restraints_container_t::restraints_container_t(int istart_res_in, int iend_res_in,
						     short int have_flanking_residue_at_start,
						     short int have_flanking_residue_at_end,
						     short int have_disulfide_residues,
						     const std::string &altloc,
						     const std::string &chain_id,
						     mmdb::Manager *mol_in, 
						     const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   init(true);
   are_all_one_atom_residues = false;
   init_from_mol(istart_res_in, iend_res_in, 
		 have_flanking_residue_at_start, 
		 have_flanking_residue_at_end,
		 have_disulfide_residues,
		 altloc,
		 chain_id, mol_in, fixed_atom_specs);
}

// Used in omega distortion graph
// 
coot::restraints_container_t::restraints_container_t(atom_selection_container_t asc_in,
						     const std::string &chain_id) {

   init(true);
   mol = asc_in.mol;
   are_all_one_atom_residues = false;

   istart_res = 999999;
   iend_res = -9999999;

   mmdb::PResidue *SelResidues = NULL;
   int nSelResidues;

   // -------- Find the max and min res no -----------------------------
   int selHnd = mol->NewSelection();
   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 1,
	       chain_id.c_str(),
	       mmdb::ANY_RES, "*",
	       mmdb::ANY_RES, "*",
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       mmdb::SKEY_NEW // selection key
	       );
   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   int resno;
   for (int ires=0; ires<nSelResidues; ires++) {
      resno = SelResidues[ires]->GetSeqNum();
      if (resno < istart_res)
	 istart_res = resno;
      if (resno > iend_res)
	 iend_res = resno;
   }
   mol->DeleteSelection(selHnd);
   // 
   // -------- Found the max and min res no -----------------------------

   // -------------------------------------------------------------------
   // Set class variables atom to the selection that includes the
   // chain (which is not the same as the input atom selection)
   //
   int SelHnd = mol->NewSelection();
   atom = NULL;
   mol->SelectAtoms(SelHnd, 0,
		    chain_id.c_str(),
		    mmdb::ANY_RES, // starting resno, an int
		    "*", // any insertion code
		    mmdb::ANY_RES, // ending resno
		    "*", // ending insertion code
		    "*", // any residue name
		    "*", // atom name
		    "*", // elements
		    "*"  // alt loc.
		    );
   mol->GetSelIndex(SelHnd, atom, n_atoms);

   // -------------------------------------------------------------------

   initial_position_params_vec.resize(3*n_atoms);
   for (int i=0; i<n_atoms; i++) {
      initial_position_params_vec[3*i  ] = atom[i]->x; 
      initial_position_params_vec[3*i+1] = atom[i]->y; 
      initial_position_params_vec[3*i+2] = atom[i]->z;
      // std::cout << "    " << i << "  " << coot::atom_spec_t(atom[i]) << "\n";
   }
}

coot::restraints_container_t::restraints_container_t(mmdb::PResidue *SelResidues, int nSelResidues,
						     const std::string &chain_id,
						     mmdb::Manager *mol_in) { 
   
   init(true);
   are_all_one_atom_residues = false;

   std::vector<coot::atom_spec_t> fixed_atoms_dummy;
   int istart_res = 999999;
   int iend_res = -9999999;
   int resno;
   
   for (int i=0; i<nSelResidues; i++) { 
      resno = SelResidues[i]->seqNum;
      if (resno < istart_res)
	 istart_res = resno;
      if (resno > iend_res)
	 iend_res = resno;
   }
   
   short int have_flanking_residue_at_start = 0;
   short int have_flanking_residue_at_end = 0;
   short int have_disulfide_residues = 0;
   const char *chn = chain_id.c_str();

   // std::cout << "DEBUG:  ==== istart_res iend_res " << istart_res << " "
   // << iend_res << std::endl; 

   init_from_mol(istart_res, iend_res, 
		 have_flanking_residue_at_start,
		 have_flanking_residue_at_end,
		 have_disulfide_residues, 
		 std::string(""), chn, mol_in, fixed_atoms_dummy);

}

coot::restraints_container_t::restraints_container_t(int istart_res_in, int iend_res_in,
						     short int have_flanking_residue_at_start,
						     short int have_flanking_residue_at_end,
						     short int have_disulfide_residues,
						     const std::string &altloc,
						     const std::string &chain_id,
						     mmdb::Manager *mol,
						     const std::vector<coot::atom_spec_t> &fixed_atom_specs,
						     const clipper::Xmap<float> &map_in,
						     float map_weight_in) {

   init(true);
   init_from_mol(istart_res_in, iend_res_in, 		 
		 have_flanking_residue_at_start, 
		 have_flanking_residue_at_end,
		 have_disulfide_residues,
		 altloc,
		 chain_id, mol, fixed_atom_specs);
   are_all_one_atom_residues = false;
   map = map_in;
   map_weight = map_weight_in;
   include_map_terms_flag = 1;

}

// 20081106 construct from a vector of residues, each of which
// has a flag attached that denotes whether or not it is a fixed
// residue (it would be set, for example in the case of flanking
// residues).
coot::restraints_container_t::restraints_container_t(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						     const std::vector<mmdb::Link> &links,
						     const coot::protein_geometry &geom,
						     mmdb::Manager *mol,
						     const std::vector<atom_spec_t> &fixed_atom_specs) {

   init(true);
   from_residue_vector = 1;
   are_all_one_atom_residues = false;
   init_from_residue_vec(residues, geom, mol, fixed_atom_specs);
}



// What are the rules for dealing with alt conf in flanking residues?
// 
// First we want to try to select only atom that have the same alt
// conf, if that fails to select atoms in flanking residues (and there
// are atoms in flanking residues (which we know from
// have_flanking_residue_at_end/start), then we should try an mmdb construction [""|"A"]
// (i.e. blank or "A").  If that fails, give up - it's a badly formed pdb file (it 
// seems to me).
// 
void
coot::restraints_container_t::init_from_mol(int istart_res_in, int iend_res_in,
					    short int have_flanking_residue_at_start,
					    short int have_flanking_residue_at_end,
					    short int have_disulfide_residues,
					    const std::string &altloc,
					    const std::string &chain_id,
					    mmdb::Manager *mol_in, 
					    const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   init_shared_pre(mol_in);

   istart_res = istart_res_in;
   iend_res   = iend_res_in;
   chain_id_save = chain_id;

   // internal flags that mirror passed have_flanking_residue_at_* variables
   istart_minus_flag = have_flanking_residue_at_start;
   iend_plus_flag    = have_flanking_residue_at_end;

   int iselection_start_res = istart_res;
   int iselection_end_res   = iend_res;
   // std::cout << "start res range: " << istart_res << " " << iend_res << " " << chain_id << "\n";

   // Are the flanking atoms available in mol_in?  (mol_in was
   // constrcted outside so the mol_in constructing routine know if
   // they were there or not.
   // 
   if (have_flanking_residue_at_start) iselection_start_res--;
   if (have_flanking_residue_at_end)   iselection_end_res++;

   SelHnd_atom = mol->NewSelection();
   mol->SelectAtoms(SelHnd_atom,
		    0,
		    chain_id.c_str(),
		    iselection_start_res, "*",
		    iselection_end_res,   "*",
		    "*", // rnames
		    "*", // anames
		    "*", // elements
		    "*"  // altLocs 
		    );

   // set the mmdb::PPAtom atom (class variable) and n_atoms:
   // 
   mol->GetSelIndex(SelHnd_atom, atom, n_atoms);

   if (0) // debugging
      for (int iat=0; iat<n_atoms; iat++)
	 std::cout << "   " << iat << "  "  << coot::atom_spec_t(atom[iat]) << "  with altloc :"
		   << altloc << ":" << std::endl;

   bool debug = 0;
   if (debug) { 
      std::cout << "DEBUG:: Selecting residues in chain \"" << chain_id << "\" gives "
		<< n_atoms << " atoms " << std::endl;
      for (int iat=0; iat<n_atoms; iat++) {
	 std::cout << "   " << iat << " " << atom[iat]->name << " "  << atom[iat]->GetSeqNum()
		   << " " << atom[iat]->GetChainID() << std::endl;
      }
   }

   // debugging, you need to uncomment mmdb.h at the top and add coords to the linking
   // atom_selection_container_t tmp_res_asc = make_asc(mol_in);
   //    std::cout << "There are " << tmp_res_asc.n_selected_atoms
   //              << " atoms in tmp_res_asc\n";
//    for (int kk=0; kk<tmp_res_asc.n_selected_atoms; kk++) { 
//       std::cout << "In simple rest " << kk << " "
//                 << tmp_res_asc.atom_selection[kk] << "\n";
//    } 

//    std::cout << "INFO::" << n_atoms << " atoms selected from molecule for refinement" ;
//    std::cout << " (this includes fixed and flanking atoms)." << std::endl;

   if (n_atoms == 0) { 
      std::cout << "ERROR:: atom selection disaster:" << std::endl;
      std::cout << "   This should not happen" << std::endl;
      std::cout << "   residue range: " << iselection_start_res << " " 
		<< iselection_end_res << " chain-id \"" << chain_id << "\" " 
		<< "flanking flags: " << have_flanking_residue_at_start 
		<< " " << have_flanking_residue_at_end << std::endl;
   }

   init_shared_post(fixed_atom_specs);
}

void
coot::restraints_container_t::init_shared_pre(mmdb::Manager *mol_in) {

   do_numerical_gradients_flag = 0;
   verbose_geometry_reporting = NORMAL;
   have_oxt_flag = false; // set in mark_OXT()
   // the smaller the alpha, the more like least squares
   geman_mcclure_alpha = 0.2; // Is this a good value? Talk to Rob.
   mol = mol_in;
   cryo_em_mode = false;
}

void
coot::restraints_container_t::init_shared_post(const std::vector<atom_spec_t> &fixed_atom_specs) {


   bonded_atom_indices.resize(n_atoms);

   initial_position_params_vec.resize(3*n_atoms); 
   for (int i=0; i<n_atoms; i++) {
      initial_position_params_vec[3*i  ] = atom[i]->x; 
      initial_position_params_vec[3*i+1] = atom[i]->y; 
      initial_position_params_vec[3*i+2] = atom[i]->z; 
   }

   // Set the UDD have_bond_or_angle to initally all "not".  They get
   // set to "have" (1) in make_restraints (and functions thereof).
   // 
   // udd_handle becomes member data so that it can be used in
   // make_restraints() without passing it back (this function is part
   // of a constructor, don't forget).
   //
   // 20131213-PE: I dont see the point of udd_bond_angle.
   // 
   if (mol) { 
      udd_bond_angle = mol->RegisterUDInteger (mmdb::UDR_ATOM, "bond or angle");
      if (udd_bond_angle < 0) { 
	 std::cout << "ERROR:: can't make udd_handle in init_from_mol\n";
      } else { 
	 for (int i=0; i<n_atoms; i++) {
	    atom[i]->PutUDData(udd_bond_angle,0);
	 }
      }
   }

   // Set the UDD of the indices in the atom array (i.e. the thing
   // that get_asc_index returns)
   // 
   if (mol) {
      udd_atom_index_handle = mol->RegisterUDInteger ( mmdb::UDR_ATOM, "atom_array_index");
      if (udd_atom_index_handle < 0) {
	 std::cout << "ERROR:: can't make udd_handle in init_from_mol\n";
      } else {
	 for (int i=0; i<n_atoms; i++) {
	    atom[i]->PutUDData(udd_atom_index_handle,i);
	    // std::cout << "init_shared_post() atom " << atom_spec_t(atom[i])
	    // << " gets udd_atom_index_handle value " << i << std::endl;
	 }
      }
   }

   use_map_gradient_for_atom.resize(n_atoms, false);
   if (! from_residue_vector) {
      // convential way
      for (int i=0; i<n_atoms; i++) {
	 if (atom[i]->residue->seqNum >= istart_res &&
	     atom[i]->residue->seqNum <= iend_res) {
	    if (! is_hydrogen(atom[i]))
	       use_map_gradient_for_atom[i] = true;
	 } else {
	    use_map_gradient_for_atom[i] = false;
	 }
      }
   } else {
      // blank out the non moving atoms (i.e. flanking residues)
      for (int i=0; i<n_atoms; i++) {
	 mmdb::Residue *res_p = atom[i]->residue;
	 if (is_a_moving_residue_p(res_p)) {
	    if (! is_hydrogen(atom[i]))
	       use_map_gradient_for_atom[i] = true;
	 } else {
	    // std::cout << "blanking out density for atom " << i << std::endl;
	    use_map_gradient_for_atom[i] = false;
	 }
      }
   }

   // z weights:
   //
   atom_z_weight.resize(n_atoms);
   std::vector<std::pair<std::string, int> > atom_list = coot::util::atomic_number_atom_list();
   for (int i=0; i<n_atoms; i++) {
      double z = coot::util::atomic_number(atom[i]->element, atom_list);
      double weight = 1.0;
      if (cryo_em_mode) {
	 // is-side-chain? would be a better test
	 if (! is_main_chain_or_cb_p(atom[i]))
	    weight = 0.3;
      }

      if (z < 0.0) {
	 std::cout << "Unknown element :" << atom[i]->element << ": " << std::endl;
	 z = weight * 6.0; // as for carbon
      } 
      atom_z_weight[i] = weight * z;
   }
   
   // the fixed atoms:   
   // 
   assign_fixed_atom_indices(fixed_atom_specs); // convert from std::vector<atom_spec_t>
   				                // to std::vector<int> fixed_atom_indices;

   // blank out those atoms from seeing electron density map gradients
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      use_map_gradient_for_atom[fixed_atom_indices[ifixed]] = false;
   } 
   
   if (verbose_geometry_reporting == VERBOSE)
      for (int i=0; i<n_atoms; i++)
	 std::cout << atom[i]->name << " " << atom[i]->residue->seqNum << " "
		   << use_map_gradient_for_atom[i] << std::endl;

} 

void
coot::restraints_container_t::init_from_residue_vec(const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
						    const coot::protein_geometry &geom,
						    mmdb::Manager *mol,
						    const std::vector<atom_spec_t> &fixed_atom_specs) {


   init_shared_pre(mol);
   residues_vec = residues;

   // Need to set class members mmdb::PPAtom atom and int n_atoms.
   // ...
   // 20090620: or do we?

   // debug:
   bool debug = false;
   if (debug) { 
      for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
	 mmdb::PAtom *res_atom_selection = NULL;
	 int n_res_atoms;
	 residues_vec[ir].second->GetAtomTable(res_atom_selection, n_res_atoms);
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " has : "
		   << n_res_atoms << " atom " << std::endl;
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " " << residue_spec_t(residues_vec[ir].second)
		   << std::endl;
	 if (false)
	    for (int iat=0; iat<n_res_atoms; iat++) {
	       mmdb::Atom *at =  res_atom_selection[iat];
	       std::cout << "DEBUG:: in init_from_residue_vec: atom "
			 << iat << " of " << n_res_atoms << " \"" 
			 << at->name << "\" \"" << at->altLoc << "\" " 
			 << at->GetSeqNum() << " \"" << at->GetInsCode() << "\" \"" 
			 << at->GetChainID() << "\"" << std::endl;
	 }
      }
   }

   // what about adding the flanking residues?  How does the atom
   // indexing of that work when (say) adding a bond?
   bonded_pair_container_t bpc = bonded_flanking_residues_by_residue_vector(geom);

   // internal variable non_bonded_neighbour_residues is set by this
   // function:
   set_non_bonded_neighbour_residues_by_residue_vector(bpc, geom);

   // std::cout << "   DEBUG:: made " << bpc.size() << " bonded flanking pairs " << std::endl;

   // passed and flanking
   // 
   std::vector<mmdb::Residue *> all_residues;
   std::vector<mmdb::Residue *>::const_iterator it;
   for (unsigned int i=0; i<residues.size(); i++) {
      all_residues.push_back(residues[i].second);
   }

   // Include only the fixed residues, because they are the flankers,
   // the other residues are the ones in the passed residues vector.
   // We don't have members of bonded_pair_container_t that are both
   // fixed.
   //
   // 20151128 only include the residues once (test that they are not there first)
   //
   int n_bonded_flankers_in_total = 0; // debug/info counter
   for (unsigned int i=0; i<bpc.size(); i++) {
      if (bpc[i].is_fixed_first) {
	 it = std::find(all_residues.begin(),
			all_residues.end(),
			bpc[i].res_1);
	 if (it == all_residues.end()) { 
	    all_residues.push_back(bpc[i].res_1);
	    n_bonded_flankers_in_total++;
	 }
      } 
      if (bpc[i].is_fixed_second) {
	 it = std::find(all_residues.begin(),
			all_residues.end(),
			bpc[i].res_2);
	 if (it == all_residues.end()) {
	    all_residues.push_back(bpc[i].res_2);
	    n_bonded_flankers_in_total++;
	 }
      } 
   }

   // Finally add the neighbour residues that are not bonded:
   for (unsigned int ires=0; ires<non_bonded_neighbour_residues.size(); ires++) {

      it = std::find(all_residues.begin(),
		     all_residues.end(),
		     non_bonded_neighbour_residues[ires]);
      if (it == all_residues.end())
	 all_residues.push_back(non_bonded_neighbour_residues[ires]);
   }

   if (0) {
      std::cout << "   DEBUG:: There are " << residues.size() << " passed residues and "
		<< all_residues.size() << " residues total (including flankers)"
		<< " with " << non_bonded_neighbour_residues.size()
		<< " non-bonded neighbours" << std::endl;
      std::cout << "These are the non-bonded neighbour residues " << std::endl;
      for (unsigned int i=0; i<non_bonded_neighbour_residues.size(); i++) { 
	 std::cout << "   " << residue_spec_t(non_bonded_neighbour_residues[i]) << std::endl;
      }
   }

   n_atoms = 0;

   // usualy this is reset in the following loop
   are_all_one_atom_residues = true;  // all refining residue (in residues) that is - do
                                      // not consider flanking residues here.
   // Now, is everything to be refined a water or MG or some such?  If
   // not (as is most often the case) are_all_one_atom_residues is
   // false.
   //
   // Something of a kludge because this will fail for residues with
   // alt-conf waters (but perhaps not, because before we get here,
   // we have done an atom selection so that there is only one alt
   // conf in each residue).
   // 
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      if (residues[ires].second->GetNumberOfAtoms() > 1) {
	 are_all_one_atom_residues = false;
	 break;
      } 
   }
   
   
   for (unsigned int i=0; i<all_residues.size(); i++)
      n_atoms += all_residues[i]->GetNumberOfAtoms();
   atom = new mmdb::PAtom[n_atoms];
   int atom_index = 0;
   for (unsigned int i=0; i<all_residues.size(); i++) {
      mmdb::PPAtom residue_atoms = 0;
      int n_res_atoms;
      all_residues[i]->GetAtomTable(residue_atoms, n_res_atoms);
      for (int iat=0; iat<n_res_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 atom[atom_index] = at;
	 atom_index++;
      }
   }

   init_shared_post(fixed_atom_specs); // use n_atoms, fills fixed_atom_indices

   if (0) { 
      std::cout << "---- after init_shared_post(): here are the "<< fixed_atom_indices.size()
		<< " fixed atoms " << std::endl;
      for (unsigned int i=0; i<fixed_atom_indices.size(); i++)
	 std::cout << "    " << i << " " << atom_spec_t(atom[fixed_atom_indices[i]]) << std::endl;
   }
   
   add_fixed_atoms_from_flanking_residues(bpc);

   if (debug) {
      std::cout << "DEBUG:: Selecting residues gives " << n_atoms << " atoms " << std::endl;
      for (int iat=0; iat<n_atoms; iat++) {
	 bool fixed_flag = false;
	 if (std::find(fixed_atom_indices.begin(),
		       fixed_atom_indices.end(), iat) != fixed_atom_indices.end())
	    fixed_flag = true;
	 std::cout << "   " << std::setw(3) << iat << " " << atom[iat]->name << " "
		   << atom[iat]->GetSeqNum() << " " << atom[iat]->GetChainID() << " "
		   << atom[iat]->GetResName() << " fixed: " << fixed_flag << std::endl;
      }
   }
   
   
}



void
coot::restraints_container_t::assign_fixed_atom_indices(const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   fixed_atom_indices.clear();
//    std::cout << "Finding atom indices for " << fixed_atom_specs.size()
// 	     << " fixed atoms " << std::endl;
   for (unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      for (int iat=0; iat<n_atoms; iat++) {
	 if (fixed_atom_specs[i].matches_spec(atom[iat])) {
	    fixed_atom_indices.push_back(iat);
	 }
      }
   }
   //    std::cout << "Found indices for " << fixed_atom_indices.size()
   // << " fixed atoms" << std::endl;
}


void
coot::restraints_container_t::debug_atoms() const {

   std::cout << "---- " << n_atoms << " atoms" << std::endl;
   for (int iat=0; iat<n_atoms; iat++) {
      std::cout << iat << " " << atom_spec_t(atom[iat]) << "  "
		<< atom[iat]->x << " "
		<< atom[iat]->y << " "
		<< atom[iat]->z << std::endl;
   }
}


// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
// 
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags) {

   short int print_chi_sq_flag = 1;
   return minimize(usage_flags, 1000, print_chi_sq_flag);

}


#include <gsl/gsl_blas.h> // for debugging norm of gradient
 
// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
//
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags, 
				       int nsteps_max,
				       short int print_initial_chi_sq_flag) {

   restraints_usage_flag = usage_flags;
   // restraints_usage_flag = BONDS_AND_ANGLES;
   // restraints_usage_flag = GEMAN_MCCLURE_DISTANCE_RESTRAINTS;
   // restraints_usage_flag = NO_GEOMETRY_RESTRAINTS;
   
   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;


   // check that we have restraints before we start to minimize:
   if (restraints_vec.size() == 0) {
      if (restraints_usage_flag != NO_GEOMETRY_RESTRAINTS) {
      std::cout << "SPECIFICATION ERROR:  There are no restraints. ";
      std::cout << "No minimization will happen" << std::endl;
      return coot::refinement_results_t(0, 0, "No Restraints!");
      }
   } 
   
   setup_gsl_vector_variables();  //initial positions

   setup_multimin_func(); // provide functions for f, df, fdf

   // T = gsl_multimin_fdfminimizer_conjugate_fr; // not as good as pr
   // T = gsl_multimin_fdfminimizer_steepest_descent; // pathetic
   // T = gsl_multimin_fminimizer_nmsimplex; // you can't just drop this in, 
                                             // because simplex is a 
                                             // non-gradient method. 
   // T = gsl_multimin_fdfminimizer_vector_bfgs;
   T = gsl_multimin_fdfminimizer_conjugate_pr;

   //   This is the Polak-Ribiere conjugate gradient algorithm. It is
   //   similar to the Fletcher-Reeves method, differing only in the
   //   choice of the coefficient \beta. Both methods work well when
   //   the evaluation point is close enough to the minimum of the
   //   objective function that it is well approximated by a quadratic
   //   hypersurface.

   s = gsl_multimin_fdfminimizer_alloc (T, n_variables());

   // restraints_usage_flag = BONDS; 
   // restraints_usage_flag = BONDS_AND_ANGLES; 
   // restraints_usage_flag = BONDS_ANGLES_AND_TORSIONS;
   // restraints_usage_flag = BONDS_ANGLES_TORSIONS_AND_PLANES;
   // restraints_usage_flag = BONDS_ANGLES_AND_PLANES;

   // restraints_usage_flag = BONDS_MASK; 
   // restraints_usage_flag = ANGLES_MASK; 
   // restraints_usage_flag = TORSIONS_MASK;
   // restraints_usage_flag = PLANES_MASK;
   
   // We get ~1ms/residue with bond and angle terms and no density terms.
   // 
   // for (int i=0; i<100; i++) { // time testing

   float tolerance = 0.06; // was 0.035
   if (! include_map_terms())
      tolerance = 0.18;

   double step_size = 0.1 * gsl_blas_dnrm2(x);

   // std::cout << ":::: starting with step_size " << step_size << std::endl;

   gsl_multimin_fdfminimizer_set (s, &multimin_func, x, step_size, tolerance);

   if (print_initial_chi_sq_flag) { 
      double d = coot::distortion_score(x,(double *)this);
      if (verbose_geometry_reporting != QUIET)
	 std::cout << "initial distortion_score: " << d << std::endl; 
      chi_squareds("Initial RMS Z values", s->x);
   }

//      std::cout << "pre minimization atom positions\n";
//      for (int i=0; i<n_atoms*3; i+=3)
//         std::cout << i/3 << " ("
// 		  << gsl_vector_get(x, i) << ", "
// 		  << gsl_vector_get(x, i+1) << ", "
// 		  << gsl_vector_get(x, i+2) << ")"
// 		  << std::endl;


   size_t iter = 0; 
   int status;
   std::vector<coot::refinement_lights_info_t> lights_vec;
   bool done_final_chi_squares = false;
   do
      {
	 iter++;
	 status = gsl_multimin_fdfminimizer_iterate(s);
// 	 std::cout << "debug:: iteration number " << iter << " of " << nsteps_max
// 		   << " status from gsl_multimin_fdfminimizer_iterate() " << status << std::endl;

	 if (status) {
	    std::cout << "unexpected error from gsl_multimin_fdfminimizer_iterate" << endl;
	    if (status == GSL_ENOPROG) {
	       cout << "Error in gsl_multimin_fdfminimizer_iterate was GSL_ENOPROG"
		    << endl; 
	       lights_vec = chi_squareds("Final Estimated RMS Z Scores", s->x);
	    }
	    break;
	 }

	 // back of envelope calculation suggests g_crit = 0.1 for
	 // coordinate shift of 0.001:  So let's choose 0.05
	 // when we have 100 restraints, 0.5 is OK.
	 // wehen we have 200,000 restraints, 0.5 is not OK.
	 // grad_lim = sqrt(n_restraints) * 0.05
	 double grad_lim = sqrt(size()) * 0.15;
	 if (grad_lim < 0.3)
	    grad_lim = 0.3;
	 status = gsl_multimin_test_gradient (s->gradient, grad_lim);

	 if (false) { // debug
	    double norm = gsl_blas_dnrm2(s->gradient); 
	    std::cout << "debug:: iteration number " << iter << " of " << nsteps_max
		      << " status from gsl_multimin_test_gradient() " << status << " for norm "
		      << norm << std::endl;
	 }
	 
	 if (status == GSL_SUCCESS) {
	    if (verbose_geometry_reporting != QUIET) { 
	       std::cout << "Minimum found (iteration number " << iter << ") at ";
	       std::cout << s->f << "\n";
	    }
	 }
	 
	 if (status == GSL_SUCCESS || status == GSL_ENOPROG) {
	    std::string title = "Final Estimated RMS Z Scores:";
	    if (status == GSL_ENOPROG)
	       title = "(No Progress) Final Estimated RMS Z Scores:";
	    std::vector<coot::refinement_lights_info_t> results = chi_squareds(title, s->x);
	    lights_vec = results;
	    done_final_chi_squares = true;
	 }

	 if (verbose_geometry_reporting == VERBOSE)
	    cout << "iteration number " << iter << " " << s->f << endl;

      }
   while ((status == GSL_CONTINUE) && (int(iter) < nsteps_max));

   if (false)
      std::cout << " in minimize() done_final_chi_squares: " << done_final_chi_squares
		<< " and status " << status
		<< " c.f GSL_SUCCESS " << GSL_SUCCESS
		<< " c.f GSL_CONTINUE " << GSL_CONTINUE
		<< " c.f GSL_ENOPROG " << GSL_ENOPROG
		<< std::endl;
   
   if (! done_final_chi_squares)
      if (status != GSL_CONTINUE)
	 chi_squareds("Final Estimated RMS Z Scores:", s->x);

   // } time testing
   update_atoms(s->x); // do OXT here
   
   // if there were bad Hs at the end of refinement
   if (status != GSL_ENOPROG) {
      if (check_pushable_chiral_hydrogens(s->x)) {
	 update_atoms(s->x);
      }
      // check and correct them if needed. 
      if (check_through_ring_bonds(s->x)) {
	 update_atoms(s->x);
      }
   }

   gsl_multimin_fdfminimizer_free(s);
   gsl_vector_free(x);
   
   // (we don't get here unless restraints were found)
   coot::refinement_results_t rr(1, status, lights_vec);

   if (false)
      std::cout << "DEBUG:: returning from minimize() with progress/status " << rr.progress
		<< ":" << std::endl;
   
   return rr;
}


std::ostream &
coot::operator<<(std::ostream &s, const simple_restraint &r) {

   s << "{restraint: ";
   if (r.restraint_type == coot::BOND_RESTRAINT)
      s << "Bond   ";
   if (r.restraint_type == coot::ANGLE_RESTRAINT)
      s << "Angle  ";
   if (r.restraint_type == coot::TORSION_RESTRAINT)
      s << "Torsion";
   if (r.restraint_type == coot::PLANE_RESTRAINT)
      s << "Plane  ";
   if (r.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT)
      s << "NBC    ";
   if (r.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
      s << "Chiral ";
      s << r.atom_index_centre;
   }
   if (r.restraint_type == coot::RAMACHANDRAN_RESTRAINT)
      s << "Rama   ";
   s << "}";
   return s;
}


void
coot::restraints_container_t::adjust_variables(const atom_selection_container_t &asc) { 

}

// 
double
starting_structure_diff_score(const gsl_vector *v, void *params) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   double d;
   double dist = 0; 

   // if (v->size != restraints->init_positions_size() ) {

   for (int i=0; i<restraints->init_positions_size(); i++) { 
      d = restraints->initial_position(i) - gsl_vector_get(v, i);
      dist += 0.01*d*d;
   }
   cout << "starting_structure_diff_score: " << dist << endl; 
   return dist; 
}


std::vector<coot::refinement_lights_info_t>
coot::restraints_container_t::chi_squareds(std::string title, const gsl_vector *v) const {

   bool print_summary = true;
   if (verbose_geometry_reporting == QUIET) print_summary = false;
   
   std::vector<refinement_lights_info_t> lights_vec;
   int n_bond_restraints = 0; 
   int n_angle_restraints = 0; 
   int n_torsion_restraints = 0; 
   int n_plane_restraints = 0; 
   int n_non_bonded_restraints = 0;
   int n_chiral_volumes = 0;
   int n_rama_restraints = 0;
   int n_start_pos_restraints = 0;
   int n_geman_mcclure_distance = 0;

   double bond_distortion = 0; 
   double gm_distortion = 0; 
   double angle_distortion = 0; 
   double torsion_distortion = 0; 
   double plane_distortion = 0; 
   double non_bonded_distortion = 0;
   double chiral_vol_distortion = 0;
   double rama_distortion = 0;
   double start_pos_distortion = 0;
   
   void *params = (double *)this;
   std::pair<int, double> dist_max_bonds(0,0);
   std::pair<int, double> dist_max_angles(0,0);
   std::pair<int, double> dist_max_planes(0,0);
   std::pair<int, double> dist_max_nbc(0,0);

   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & BONDS_MASK) {
	 if ( restraints_vec[i].restraint_type == BOND_RESTRAINT) {
	    n_bond_restraints++;
	    // 	    bond_distortion += distortion_score_bond(restraints_vec[i], v);
	    double dist = distortion_score_bond(restraints_vec[i], v);
	    const simple_restraint &rest = restraints_vec[i];
	    bond_distortion += dist;
	    if (dist > dist_max_bonds.second) {
	       dist_max_bonds.first = i;
	       dist_max_bonds.second = dist;
	    }
	 }
      }
      
      if (restraints_usage_flag & GEMAN_MCCLURE_DISTANCE_MASK) {
	 if ( restraints_vec[i].restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    n_geman_mcclure_distance++;
	    double d = distortion_score_geman_mcclure_distance(restraints_vec[i], v, geman_mcclure_alpha);
	    gm_distortion += d;
	    // std::cout << "distortion_score_geman_mcclure_distance " << i << " " << d << "\n";
	 }
      }

      if (restraints_usage_flag & ANGLES_MASK) { // 2: angles
	 if ( restraints_vec[i].restraint_type == coot::ANGLE_RESTRAINT) {
	    n_angle_restraints++;
	    double dist = coot::distortion_score_angle(restraints_vec[i], v);
	    angle_distortion += dist;
	    if (dist > dist_max_angles.second) {
	       dist_max_angles.first = i;
	       dist_max_angles.second = dist;
	    }
	 }
      }

      if (restraints_usage_flag & TORSIONS_MASK) { // 4: torsions
	 if ( restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) {
	    try { 
	       torsion_distortion += coot::distortion_score_torsion(restraints_vec[i], v); 
	       n_torsion_restraints++;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "WARNING:: caught runtime_error " << rte.what() << std::endl;
	    } 
	 }
      }

      if (restraints_usage_flag & PLANES_MASK) { // 8: planes
	 if ( restraints_vec[i].restraint_type == coot::PLANE_RESTRAINT) {
	    n_plane_restraints++;
	    double dist = coot::distortion_score_plane(restraints_vec[i], v); 
	    plane_distortion += dist;
	    if (dist > dist_max_planes.second) {
	       dist_max_planes.first = i;
	       dist_max_planes.second = dist;
	    }
            if (false) {  // debugging plane restraints.
                std::cout << " plane distortion " << i << " " 
                          << coot::distortion_score_plane(restraints_vec[i], v) << " " 
                          << restraints_vec[i];
                for (unsigned int jj = 0; jj<restraints_vec[i].plane_atom_index.size(); jj+=3) { 
                    std::cout << "\n                                ";
		    unsigned int idx = restraints_vec[i].plane_atom_index[jj].first;
                    std::cout << idx << " " << coot::atom_spec_t(atom[idx]);
                    if ((jj+1) < restraints_vec[i].plane_atom_index.size()) { 
		       unsigned int idx_1 = restraints_vec[i].plane_atom_index[jj+1].first;
                       std::cout << " " << idx_1 << " " << coot::atom_spec_t(atom[idx_1]);
                    }
                    if ((jj+2) < restraints_vec[i].plane_atom_index.size()) { 
		       unsigned int idx_2 = restraints_vec[i].plane_atom_index[jj+1].first;
                       std::cout << " " << idx_2 << " " << coot::atom_spec_t(atom[idx_2]);
                    }
                }
                std::cout << std::endl;
            }
	 }
      }

      if (restraints_usage_flag & coot::NON_BONDED_MASK) { 
	 if ( restraints_vec[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 
	    n_non_bonded_restraints++;
	    double dist = coot::distortion_score_non_bonded_contact(restraints_vec[i], v);
	    non_bonded_distortion += dist;
	    if (dist > dist_max_nbc.second) {
	       dist_max_nbc.first = i;
	       dist_max_nbc.second = dist;
	    }
	 }
      }

      if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) { 
  	 if ( restraints_vec[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
  	    n_chiral_volumes++;
  	    chiral_vol_distortion += coot::distortion_score_chiral_volume(restraints_vec[i], v);
  	 }
      }

      if (restraints_usage_flag & coot::RAMA_PLOT_MASK) {
  	 if ( restraints_vec[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) { 
  	    n_rama_restraints++;
  	    rama_distortion += coot::distortion_score_rama(restraints_vec[i], v, lograma);
  	 }
      }
      
      if ( restraints_vec[i].restraint_type == coot::START_POS_RESTRAINT) {
         n_start_pos_restraints++;
         start_pos_distortion += coot::distortion_score_start_pos(restraints_vec[i], params, v);
      }
   }

   std::string r = "";

   r += title;
   r += "\n";
   if (print_summary)
      std::cout << "    " << title << std::endl;
   if (n_bond_restraints == 0) {
      if (print_summary)
	 std::cout << "bonds:      N/A " << std::endl;
   } else {
      double bd = bond_distortion/double(n_bond_restraints);
      double sbd = 0.0;
      if (bd > 0)
	 sbd = sqrt(bd);
      if (print_summary)
	 std::cout << "bonds:      " << sbd << std::endl;
      r += "   bonds:  ";
      r += coot::util::float_to_string_using_dec_pl(sbd, 3);
      r += "\n";
      std::string s = "Bonds:  ";
      s += coot::util::float_to_string_using_dec_pl(sbd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Bonds", s, sbd));
   } 
   if (n_angle_restraints == 0) {
      if (print_summary)
	 std::cout << "angles:     N/A " << std::endl;
   } else {
      double ad = angle_distortion/double(n_angle_restraints);
      double sad = 0.0;
      if (ad > 0.0)
	 sad = sqrt(ad);
      if (print_summary)
	 std::cout << "angles:     " << sad << std::endl;
      r += "   angles: ";
      r += coot::util::float_to_string_using_dec_pl(sad, 3);
      r += "\n";
      std::string s = "Angles: ";
      s += coot::util::float_to_string_using_dec_pl(sad, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Angles", s, sad));
   } 
   if (n_torsion_restraints == 0) {
      if (print_summary)
	 std::cout << "torsions:   N/A " << std::endl;
   } else {
      double td = torsion_distortion/double(n_torsion_restraints);
      double std = 0.0;
      if (td > 0.0)
	 std = sqrt(td);
      if (print_summary)
	 std::cout << "torsions:   " << std << std::endl;
      r += "   torsions: ";
      r += coot::util::float_to_string_using_dec_pl(std, 3);
      r += "\n";
      std::string s = "Torsions: ";
      s += coot::util::float_to_string_using_dec_pl(std, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Torsions", s, std));
   } 
   if (n_plane_restraints == 0) {
      if (print_summary)
	 std::cout << "planes:     N/A " << std::endl;
   } else {
      double pd = plane_distortion/double(n_plane_restraints);
      double spd = 0.0;
      if (pd > 0.0)
	 spd = sqrt(pd);
      if (print_summary)
	 std::cout << "planes:     " << spd << std::endl;
      r += "   planes: ";
      r += coot::util::float_to_string_using_dec_pl(spd, 3);
      r += "\n";
      std::string s = "Planes: ";
      s += coot::util::float_to_string_using_dec_pl(spd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Planes", s, spd));
   }
   if (n_non_bonded_restraints == 0) {
      if (print_summary)
	 std::cout << "non-bonded: N/A " << std::endl;
   } else {
      double nbd = non_bonded_distortion/double(n_non_bonded_restraints);
      double snbd = 0.0;
      if (nbd > 0.0)
	 snbd = sqrt(nbd);
      if (print_summary)
	 std::cout << "non-bonded: " << nbd << std::endl;
      r += "   non-bonded: ";
      r += coot::util::float_to_string_using_dec_pl(snbd, 3);
      r += "\n";
      std::string s = "Non-bonded: ";
      s += coot::util::float_to_string_using_dec_pl(snbd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Non-bonded", s, snbd));
   }
   if (n_chiral_volumes == 0) { 
      if (print_summary)
	 std::cout << "chiral vol: N/A " << std::endl;
   } else {
      double cd = chiral_vol_distortion/double(n_chiral_volumes);
      double scd = 0.0;
      if (cd > 0.0)
	 scd = sqrt(cd);
      if (print_summary)
	 std::cout << "chiral vol: " << scd << std::endl;
      r += "   chirals: ";
      r += coot::util::float_to_string_using_dec_pl(scd, 3);
      std::string s = "Chirals: ";
      s += coot::util::float_to_string_using_dec_pl(scd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Chirals", s, scd));
   }
   if (n_rama_restraints == 0) { 
      if (print_summary)
	 std::cout << "rama plot:  N/A " << std::endl;
   } else {
      double rd = rama_distortion/double(n_rama_restraints);
      if (print_summary)
	 std::cout << "rama plot:  " << rd << std::endl;
      r += "   rama plot: ";
      r += coot::util::float_to_string_using_dec_pl(rd, 3);
      std::string s = "Rama Plot: ";
      s += coot::util::float_to_string_using_dec_pl(rd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Rama", s, rd));
   }
   if (n_start_pos_restraints == 0) {
      if (print_summary)
	 std::cout << "start_pos:  N/A " << std::endl;
   } else {
      double spd = start_pos_distortion/double(n_start_pos_restraints);
      double sspd = 0.0;
      if (spd > 0.0)
	 sspd = sqrt(spd);
      if (print_summary)
	 std::cout << "start_pos:  " << sspd << std::endl;
      r += "startpos:  ";
      r += coot::util::float_to_string_using_dec_pl(sspd, 3);
      r += "\n";
      std::string s = "Start pos: ";
      s += coot::util::float_to_string_using_dec_pl(sspd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Start_pos", s, sspd));
   } 
   if (n_geman_mcclure_distance == 0) {
      if (print_summary)
	 std::cout << "GemanMcCl:  N/A " << std::endl;
   } else {
      double spd = gm_distortion/double(n_geman_mcclure_distance);
      double sspd = 0.0;
      if (spd > 0.0)
	 sspd = sqrt(spd);
      if (print_summary)
	 std::cout << "GemanMcCl:  " << sspd << " from " << n_geman_mcclure_distance << " distances" << std::endl;
      r += "GemanMcCl:  ";
      r += coot::util::float_to_string_using_dec_pl(sspd, 3);
      r += "\n";
      std::string s = "GemanMcCl: ";
      s += coot::util::float_to_string_using_dec_pl(sspd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("GemanMcCl", s, sspd));
   } 
   return lights_vec;
} 


// public
coot::model_bond_deltas
coot::restraints_container_t::resolve_bonds() {

   setup_gsl_vector_variables();
   return resolve_bonds(x);

}

coot::model_bond_deltas
coot::restraints_container_t::resolve_bonds(const gsl_vector *v) const {

   model_bond_deltas resultant;

   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & BONDS_MASK) {
	 const simple_restraint &rest = restraints_vec[i];
	 if (rest.restraint_type == BOND_RESTRAINT) {

	    int idx = 3*(rest.atom_index_1);
	    clipper::Coord_orth a1(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));
	    idx = 3*(rest.atom_index_2);
	    clipper::Coord_orth a2(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));
	    double ideal = rest.target_value;
	    double bl = clipper::Coord_orth::length(a1,a2);
	    double delta = bl - ideal;
	    clipper::Coord_orth b_uv((a2-a1).unit());
	    clipper::Coord_orth b_uv_abs(std::fabs(b_uv.x()),
					 std::fabs(b_uv.y()),
					 std::fabs(b_uv.z()));
	    clipper::Coord_orth frag(b_uv_abs * delta);
	    resultant.add(clipper::Coord_orth(b_uv_abs * delta));
	 }
      }
   }
   return resultant;
}


// Ah, but (c.f. distortion) we want to return a low value for a good
// fit and a high one for a bad.
double
coot::electron_density_score(const gsl_vector *v, void *params) { 

   // We weight and sum to get the score and negate.  That will do?
   // 
   double score = 0; 
   // double e = 2.718281828; 
   
   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (restraints->include_map_terms() == 1) { 
      
      // convert from variables to coord_orths of where the atoms are
      
      for (unsigned int i=0; i< v->size; i += 3) { 
	 int iat = i/3;
	 if (restraints->use_map_gradient_for_atom[iat]) {
	    bool use_it = 1;
// 	    for (unsigned int ifixed=0; ifixed<restraints->fixed_atom_indices.size(); ifixed++) {
// 	       if (restraints->fixed_atom_indices[ifixed] == iat) { 
// 		  std::cout << "ignoring density term for atom " << iat << std::endl;
// 		  use_it = 0;
// 		  break;
// 	       }
// 	    }
	    if (use_it) { 
	       clipper::Coord_orth ao(gsl_vector_get(v,i), 
				      gsl_vector_get(v,i+1), 
				      gsl_vector_get(v,i+2));
	       
	       score += restraints->Map_weight() *
		  restraints->atom_z_weight[iat] *
		  restraints->electron_density_score_at_point(ao);
	    }
	 }
      }
   }
   
   // return pow(e,-score*0.01);
   return -score;

}

// Note that the gradient for the electron density is opposite to that
// of the gradient for the geometry (consider a short bond on the edge
// of a peak - in that case the geometry gradient will be negative as
// the bond is lengthened and the electron density gradient will be
// positive).
//
// So we want to change that positive gradient for a low score when
// the atoms cooinside with the density - hence the contributions that
// we add are negated.
// 
void coot::my_df_electron_density(const gsl_vector *v, 
				  void *params, 
				  gsl_vector *df) {

#ifdef ANALYSE_REFINEMENT_TIMING
   timeval start_time;
   timeval current_time;
   gettimeofday(&start_time, NULL);
#endif // ANALYSE_REFINEMENT_TIMING

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints_p = static_cast<restraints_container_t *> (params); 

   if (restraints_p->include_map_terms() == 1) { 
      
#ifdef HAVE_CXX_THREAD

      std::atomic<unsigned int> done_count_for_threads(0);

      if (restraints_p->thread_pool_p) {
	 int idx_max = (v->size)/3;
	 unsigned int n_per_thread = idx_max/restraints_p->n_threads;

	 for (unsigned int i_thread=0; i_thread<restraints_p->n_threads; i_thread++) {
	    int idx_start = i_thread * n_per_thread;
	    int idx_end   = idx_start + n_per_thread;
	    // for the last thread, set the end atom index
	    if (i_thread == (restraints_p->n_threads - 1))
	       idx_end = idx_max; // for loop uses iat_start and tests for < iat_end

	    restraints_p->thread_pool_p->push(my_df_electron_density_threaded_single,
					      v, restraints_p, df, idx_start, idx_end,
					      std::ref(done_count_for_threads));

	 }

	 while (done_count_for_threads != restraints_p->n_threads) {
// 	    std::cout << "comparing " << restraints_p->done_count_for_threads
// 		      << " "  << restraints_p->n_threads << std::endl;
	    std::this_thread::sleep_for(std::chrono::microseconds(1));
	 }

      } else {
	 my_df_electron_density_single(v, restraints_p, df, 0, v->size/3);
      }
#else
      my_df_electron_density_single(v, restraints_p, df, 0, v->size/3);
#endif // HAVE_CXX_THREAD      

   }
#ifdef ANALYSE_REFINEMENT_TIMING
   gettimeofday(&current_time, NULL);
   double td = current_time.tv_sec - start_time.tv_sec;
   td *= 1000.0;
   td += double(current_time.tv_usec - start_time.tv_usec)/1000.0;
   std::cout << "------------- mark my_df_electron_density: " << td << std::endl;
#endif // ANALYSE_REFINEMENT_TIMING
}

// Note that the gradient for the electron density is opposite to that
// of the gradient for the geometry (consider a short bond on the edge
// of a peak - in that case the geometry gradient will be negative as
// the bond is lengthened and the electron density gradient will be
// positive).
//
// So we want to change that positive gradient for a low score when
// the atoms cooinside with the density - hence the contributions that
// we add are negated.
// 
void coot::my_df_electron_density_old_2017(const gsl_vector *v, 
					   void *params, 
					   gsl_vector *df) {

#ifdef ANALYSE_REFINEMENT_TIMING
   timeval start_time;
   timeval current_time;
   gettimeofday(&start_time, NULL);
#endif // ANALYSE_REFINEMENT_TIMING

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints = static_cast<restraints_container_t *> (params); 

   if (restraints->include_map_terms() == 1) { 

      clipper::Grad_orth<double> grad_orth;
      float scale = restraints->Map_weight();
      float zs;
      
      for (unsigned int i=0; i<v->size; i+=3) {

	 int iat = i/3;
	 
	 if (restraints->use_map_gradient_for_atom[iat]) {

	    clipper::Coord_orth ao(gsl_vector_get(v,i), 
				   gsl_vector_get(v,i+1), 
				   gsl_vector_get(v,i+2));
	    
	    // std::cout << "gradients: " << grad_orth.format() << std::endl;
	    // std::cout << "adding to " << gsl_vector_get(df, i  ) << std::endl;
	    // add this density term to the gradient
	    //
	    // 
	    grad_orth = restraints->electron_density_gradient_at_point(ao);
	    zs = scale * restraints->atom_z_weight[iat];

	    if (0) {
	       std::cout << "electron density df: adding "
			 <<  - zs * grad_orth.dx() << " "
			 <<  - zs * grad_orth.dy() << " "
			 <<  - zs * grad_orth.dz() << " to "
			 <<  gsl_vector_get(df, i  ) << " "
			 <<  gsl_vector_get(df, i+1) << " "
			 <<  gsl_vector_get(df, i+2) << "\n";
	    }
	    
// 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
// 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
// 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());
	    
	    *gsl_vector_ptr(df, i  ) -= zs * grad_orth.dx();
	    *gsl_vector_ptr(df, i+1) -= zs * grad_orth.dy();
	    *gsl_vector_ptr(df, i+2) -= zs * grad_orth.dz();

	 } else {
	    // atom is private	    
// 	    std::cout << "  Not adding elecron density for atom "
// 		      << restraints->atom[iat]->GetChainID() << " "
// 		      << restraints->atom[iat]->GetSeqNum() << " "
// 		      << restraints->atom[iat]->GetAtom() << std::endl;
	 } 
      }
   }
#ifdef ANALYSE_REFINEMENT_TIMING
   gettimeofday(&current_time, NULL);
   double td = current_time.tv_sec - start_time.tv_sec;
   td *= 1000.0;
   td += double(current_time.tv_usec - start_time.tv_usec)/1000.0;
   std::cout << "------------- mark my_df_electron_density: " << td << std::endl;
#endif // ANALYSE_REFINEMENT_TIMING
}


#ifdef HAVE_CXX_THREAD

// restraints are modified by atomic done_count_for_threads changing.
//
void coot::my_df_electron_density_threaded_single(int thread_idx, const gsl_vector *v,
						  coot::restraints_container_t *restraints,
						  gsl_vector *df,
						  int atom_idx_start, int atom_idx_end,
						  std::atomic<unsigned int> &done_count_for_threads) {

   for (int iat=atom_idx_start; iat<atom_idx_end; ++iat) {
      if (restraints->use_map_gradient_for_atom[iat]) {

	 int idx = 3 * iat;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx), 
				gsl_vector_get(v,idx+1), 
				gsl_vector_get(v,idx+2));
	    
	 clipper::Grad_orth<double> grad_orth = restraints->electron_density_gradient_at_point(ao);
	 float zs = restraints->Map_weight() * restraints->atom_z_weight[iat];

	 if (0) { 
	    std::cout << "electron density df: adding "
		      <<  - zs * grad_orth.dx() << " "
		      <<  - zs * grad_orth.dy() << " "
		      <<  - zs * grad_orth.dz() << " to "
		      <<  gsl_vector_get(df, idx  ) << " "
		      <<  gsl_vector_get(df, idx+1) << " "
		      <<  gsl_vector_get(df, idx+2) << "\n";
	 }
	    
	 // 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	 // 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	 // 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());

	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
      }
   }
   ++done_count_for_threads; // atomic
}
#endif	 


//
void coot::my_df_electron_density_single(const gsl_vector *v,
					 coot::restraints_container_t *restraints,
					 gsl_vector *df,
					 int atom_idx_start, int atom_idx_end) {

   for (int iat=atom_idx_start; iat<atom_idx_end; ++iat) {
      if (restraints->use_map_gradient_for_atom[iat]) {

	 int idx = 3 * iat;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx), 
				gsl_vector_get(v,idx+1), 
				gsl_vector_get(v,idx+2));
	    
	 clipper::Grad_orth<double> grad_orth = restraints->electron_density_gradient_at_point(ao);
	 float zs = restraints->Map_weight() * restraints->atom_z_weight[iat];

	 if (0) { 
	    std::cout << "electron density df: adding "
		      <<  - zs * grad_orth.dx() << " "
		      <<  - zs * grad_orth.dy() << " "
		      <<  - zs * grad_orth.dz() << " to "
		      <<  gsl_vector_get(df, idx  ) << " "
		      <<  gsl_vector_get(df, idx+1) << " "
		      <<  gsl_vector_get(df, idx+2) << "\n";
	 }
	    
	 // 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	 // 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	 // 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());

	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
      }
   }
}
void coot::my_df_electron_density_old (gsl_vector *v, 
				       void *params, 
				       gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (restraints->include_map_terms() == 1) { 

      double new_S_minu, new_S_plus, tmp, val; 

      cout << "density_gradients" << endl; 
      for (unsigned int i=0; i<v->size; i++) { 
      
	 tmp = gsl_vector_get(v, i); 
	 gsl_vector_set(v, i, tmp+0.01); 
	 new_S_plus = coot::electron_density_score(v, params); 
	 gsl_vector_set(v, i, tmp-0.01); 
	 new_S_minu = coot::electron_density_score(v, params);
	 // new_S_minu = 2*tmp - new_S_plus; 

	 // restore the initial value: 
	 gsl_vector_set(v, i, tmp);

	 val = (new_S_plus - new_S_minu)/(2*0.01); 
	 cout << "density gradient: " << i << " " << val << endl;

	 // add this density term to the gradient
	 gsl_vector_set(df, i, gsl_vector_get(df, i) + val);
      } 
   }
}

// Compute both f and df together.
void coot::my_fdf(const gsl_vector *x, void *params, 
		  double *f, gsl_vector *df) { 

   // 20170423 these can be done in parallel? ... check the timings at least.
   *f = coot::distortion_score(x, params); 
    coot::my_df(x, params, df); 
}


unsigned int
coot::restraints_container_t::test_function(const coot::protein_geometry &geom) {
   std::cout << "----- test_function() with geom of size : " << geom.size() << std::endl;
   std::cout << "    geom ref pointer " << &geom << std::endl;
   return geom.size();
} 

unsigned int
coot::restraints_container_t::const_test_function(const coot::protein_geometry &geom) const {
   std::cout << "----- const_test_function() with geom of size : " << geom.size() << std::endl;
   std::cout << "    geom ref pointer " << &geom << std::endl;
   return geom.size();
} 


// We need to fill restraints_vec (which is a vector of
// simple_restraint) using the coordinates () and the dictionary of
// restraints, protein_geometry geom.
//
// The plan is to get a list of residues, and for each of those
// residues, look in geom for a residue of that type and if found,
// fill restraints_vec appropriately.
int
coot::restraints_container_t::make_restraints(int imol,
					      const coot::protein_geometry &geom,
					      coot::restraint_usage_Flags flags_in, 
					      bool do_residue_internal_torsions,
					      bool do_trans_peptide_restraints,
					      float rama_plot_target_weight,
					      bool do_rama_plot_restraints, 
					      coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds,
					      bool do_link_restraints,
					      bool do_flank_restraints) {

   // if a peptider is trans, add a restraint to penalize non-trans configuration
   // (currently a torsion restraint on peptide w of 180)
   // 

   // debugging SRS inclusion.
   if (false) {
      std::cout << "----- make restraints() called with geom of size : " << geom.size() << std::endl;
      std::cout << "    geom ref pointer " << &geom << std::endl;
   }
   
   restraints_usage_flag = flags_in; // also set in minimize() and geometric_distortions()
   // restraints_usage_flag = BONDS_AND_ANGLES;
   // restraints_usage_flag = GEMAN_MCCLURE_DISTANCE_RESTRAINTS;
   // restraints_usage_flag = NO_GEOMETRY_RESTRAINTS;

   if (n_atoms) {

      mark_OXT(geom);
      make_monomer_restraints(imol, geom, do_residue_internal_torsions);

      bool do_link_restraints_internal = true;
      bool do_flank_restraints_internal = true;

      if (! from_residue_vector) {
	 if (istart_res == iend_res)
	    do_link_restraints_internal = false;
	 if (! istart_minus_flag && !iend_plus_flag)
	    do_flank_restraints_internal = false;
      }

      if (! do_link_restraints)
	 do_link_restraints_internal = false;
      if (! do_flank_restraints)
	 do_flank_restraints_internal = false;

      if (do_link_restraints_internal)
	 make_link_restraints(geom, do_rama_plot_restraints, do_trans_peptide_restraints);

      // don't do torsions, ramas maybe.   
      coot::bonded_pair_container_t bpc;

      if (do_flank_restraints_internal)
	 bpc = make_flanking_atoms_restraints(geom,
					      do_rama_plot_restraints,
					      do_trans_peptide_restraints);
      int iret_prev = restraints_vec.size();

      if (sec_struct_pseudo_bonds == coot::HELIX_PSEUDO_BONDS) {
	 make_helix_pseudo_bond_restraints();
      } 
      if (sec_struct_pseudo_bonds == coot::STRAND_PSEUDO_BONDS) {
	 make_strand_pseudo_bond_restraints();
      }
      if (restraints_usage_flag & coot::NON_BONDED_MASK) {
	 if ((iret_prev > 0) || are_all_one_atom_residues) {
	    reduced_angle_info_container_t ai(restraints_vec);
	    int n_nbcr = make_non_bonded_contact_restraints(imol, bpc, ai, geom);
	    if (verbose_geometry_reporting != QUIET)
	       std::cout << "INFO:: make_restraints(): made " << n_nbcr << " non-bonded restraints\n";
	 }
      }
      make_restraint_types_index_limits();
   }
   return restraints_vec.size();
}


void
coot::restraints_container_t::make_restraint_types_index_limits() {

   unsigned int unset = 9999999;
   restraints_limits_bonds = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_angles = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_torsions = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_chirals = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_planes =  std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_parallel_planes =  std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_non_bonded_contacts = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_geman_mclure = std::pair<unsigned int, unsigned int> (unset,0);
   restraints_limits_start_pos = std::pair<unsigned int, unsigned int> (unset,0);

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      const simple_restraint &restraint = restraints_vec[i];
      if (restraint.restraint_type == coot::BOND_RESTRAINT) {
	 if (restraints_limits_bonds.first == unset)
	    restraints_limits_bonds.first = i;
	 if (i > restraints_limits_bonds.second)
	    restraints_limits_bonds.second = i;
      }
      if (restraint.restraint_type == coot::ANGLE_RESTRAINT) {
	 if (restraints_limits_angles.first == unset)
	    restraints_limits_angles.first = i;
	 if (i > restraints_limits_angles.second)
	    restraints_limits_angles.second = i;
      }
      if (restraint.restraint_type == coot::TORSION_RESTRAINT) {
	 if (restraints_limits_torsions.first == unset)
	    restraints_limits_torsions.first = i;
	 if (i > restraints_limits_torsions.second)
	    restraints_limits_torsions.second = i;
      }
      if (restraint.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	 if (restraints_limits_chirals.first == unset)
	    restraints_limits_chirals.first = i;
	 if (i > restraints_limits_chirals.second)
	    restraints_limits_chirals.second = i;
      }
      if (restraint.restraint_type == coot::PLANE_RESTRAINT) {
	 if (restraints_limits_planes.first == unset)
	    restraints_limits_planes.first = i;
	 if (i > restraints_limits_planes.second)
	    restraints_limits_planes.second = i;
      }
      if (restraint.restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
	 if (restraints_limits_parallel_planes.first == unset)
	    restraints_limits_parallel_planes.first = i;
	 if (i > restraints_limits_parallel_planes.second)
	    restraints_limits_parallel_planes.second = i;
      }
      if (restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	 if (restraints_limits_non_bonded_contacts.first == unset)
	    restraints_limits_non_bonded_contacts.first = i;
	 if (i > restraints_limits_non_bonded_contacts.second)
	    restraints_limits_non_bonded_contacts.second = i;
      }
      if (restraint.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	 if (restraints_limits_geman_mclure.first == unset)
	    restraints_limits_geman_mclure.first = i;
	 if (i > restraints_limits_geman_mclure.second)
	    restraints_limits_geman_mclure.second = i;
      }
      if (restraint.restraint_type == coot::START_POS_RESTRAINT) {
	 if (restraints_limits_start_pos.first == unset)
	    restraints_limits_start_pos.first = i;
	 if (i > restraints_limits_start_pos.second)
	    restraints_limits_start_pos.second = i;
      }
   }

   // now check for unsets
   if (restraints_limits_bonds.first  == unset)   restraints_limits_bonds.first = 0;
   if (restraints_limits_angles.first == unset)   restraints_limits_angles.first = 0;
   if (restraints_limits_torsions.first == unset) restraints_limits_torsions.first = 0;
   if (restraints_limits_chirals.first == unset)  restraints_limits_chirals.first = 0;
   if (restraints_limits_planes.first == unset)   restraints_limits_planes.first = 0;
   if (restraints_limits_parallel_planes.first == unset) restraints_limits_parallel_planes.first = 0;
   if (restraints_limits_non_bonded_contacts.first == unset) restraints_limits_non_bonded_contacts.first = 0;
   if (restraints_limits_geman_mclure.first == unset) restraints_limits_geman_mclure.first = 0;
   if (restraints_limits_start_pos.first == unset) restraints_limits_start_pos.first = 0;

   if (false) {
      std::cout << "restraints limits bonds "
	        << restraints_limits_bonds.first << " " << restraints_limits_bonds.second << std::endl;
      std::cout << "restraints limits angles "
	        << restraints_limits_angles.first << " " << restraints_limits_angles.second << std::endl;
      std::cout << "restraints limits torsions "
	        << restraints_limits_torsions.first << " " << restraints_limits_torsions.second << std::endl;
      std::cout << "restraints limits chirals "
	        << restraints_limits_chirals.first << " " << restraints_limits_chirals.second << std::endl;
      std::cout << "restraints limits planes "
	        << restraints_limits_planes.first << " " << restraints_limits_planes.second << std::endl;
      std::cout << "restraints limits nbc "
	        << restraints_limits_non_bonded_contacts.first << " " << restraints_limits_non_bonded_contacts.second
	        << std::endl;
   }

}



// This only marks the first OXT atom we find that has all its
// reference atoms (which is reasonable (I hope)).
// 
void 
coot::restraints_container_t::mark_OXT(const coot::protein_geometry &geom) {

   std::string oxt(" OXT");
   for (int i=0; i<n_atoms; i++) { 
      if (std::string(atom[i]->name) == oxt) {

	 mmdb::Residue *residue = atom[i]->residue;
	 mmdb::Atom *res_atom = NULL;
	 
	 std::string res_name = residue->GetResName();
	 if (coot::util::is_standard_residue_name(res_name)) {
	    // add it if it has not been added before.
	    if (std::find(residues_with_OXTs.begin(),
			  residues_with_OXTs.end(),
			  residue) == residues_with_OXTs.end()) {
	       residues_with_OXTs.push_back(residue);
	       have_oxt_flag = true; // added 20160612 for Miguel's flying OXT bug.
	    }
	 }
      }
   }
}

bool
coot::restraints_container_t::fixed_check(int index_1) const {

   bool r = 0;
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index_1 == fixed_atom_indices[ifixed]) {
	 r = 1;
	 break;
      }
   }
   return r;
} 

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2) const {

   std::vector<bool> r(2,0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed])
	 r[0] = 1;
      if (index2 == fixed_atom_indices[ifixed])
	 r[1] = 1;
   }

   return r;
}

std::vector<bool>
coot::restraints_container_t::make_non_bonded_fixed_flags(int index1, int index2) const {

   std::vector<bool> r(2,0);
   bool set_0 = 0;
   bool set_1 = 0;
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed]) { 
	 r[0] =  true;
	 set_0 = true;
	 break;
      }
   } 
      
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index2 == fixed_atom_indices[ifixed]) { 
	 r[1] = true;
	 set_1 = true;
	 break;
      }
   }

   if (set_0 && set_1 ) {
      return r;  // yay, fast.
   }

   if (! set_0) {
      mmdb::Residue *res = atom[index1]->residue;
      if (std::find(non_bonded_neighbour_residues.begin(),
		    non_bonded_neighbour_residues.end(),
		    res) != non_bonded_neighbour_residues.end())
	 r[0] = 1; // if we found the residue in non_bonded_neighbour_residues
                   // then that atom of that residue is fixed
   }
   if (! set_1) {
      mmdb::Residue *res = atom[index2]->residue;
      if (std::find(non_bonded_neighbour_residues.begin(),
		    non_bonded_neighbour_residues.end(),
		    res) != non_bonded_neighbour_residues.end())
	 r[1] = 1; 
   }
   return r;
} 


std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2, int index3) const {

   std::vector<bool> r(3,0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed])
	 r[0] = 1;
      if (index2 == fixed_atom_indices[ifixed])
	 r[1] = 1;
      if (index3 == fixed_atom_indices[ifixed])
	 r[2] = 1;
   }
   return r;
} 

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2, int index3, int index4) const {

   std::vector<bool> r(4,0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed])
	 r[0] = 1;
      if (index2 == fixed_atom_indices[ifixed])
	 r[1] = 1;
      if (index3 == fixed_atom_indices[ifixed])
	 r[2] = 1;
      if (index4 == fixed_atom_indices[ifixed])
	 r[3] = 1;
   }
   return r;
}

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(const std::vector<int> &indices) const {

   std::vector<bool> r(indices.size(), 0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      for (unsigned int i_index=0; i_index<indices.size(); i_index++) {
	 if (indices[i_index] == fixed_atom_indices[ifixed])
	    r[i_index] = 1;
      }
   }
   return r;
} 


void
coot::restraints_container_t::make_helix_pseudo_bond_restraints() {

   // This method of making pseudo bonds relies on the residue range
   // being continuous in sequence number (seqNum) and no insertion
   // codes messing up the number scheme.  If these are not the case
   // then we will make bonds between the wrongs atoms (of the wrong
   // residues)... a more sophisticated algorithm would be needed.
   // 
   // The method ignores residues with alt confs.
   //
   float pseudo_bond_esd = 0.04; // just a guess
   int selHnd = mol->NewSelection();
   int nSelResidues;
   mmdb::PPResidue SelResidue;
   mmdb::PPAtom res_1_atoms = NULL;
   mmdb::PPAtom res_2_atoms = NULL;
   int n_res_1_atoms;
   int n_res_2_atoms;
   int index1 = -1; 
   int index2 = -1; 
   mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"",  // altLocs
		mmdb::SKEY_NEW // selection key
		);
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 0) {
      for (int i=4; i<nSelResidues; i++) {
	 // nN -> (n-4)O 2.91
         // nN -> (n-3)O 3.18
         // nO -> (n+3)N 3.18  // symmetric.  No need to specify both forwards
         // nO -> (n+4)N 2.91  // and backwards directions.
	 SelResidue[i]->GetAtomTable(res_1_atoms, n_res_1_atoms);
	 for (int iat1=0; iat1<n_res_1_atoms; iat1++) {
	    std::string at_1_name(res_1_atoms[iat1]->name);
	    
	    if (at_1_name == " N  ") {
	       mmdb::Residue *contact_res = SelResidue[i-4];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 4)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    2.91, pseudo_bond_esd, 1.2);
			std::cout << "Helix Bond restraint (" << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 2.91" << std::endl;
		     }
		  }
	       }

	       contact_res = SelResidue[i-3];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 3)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    3.18, pseudo_bond_esd, 1.2);
			std::cout << "Helix Bond restraint (" << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 3.18" << std::endl;
		     }
		  }
	       }
	       
	    } 
	 } 
      }
   }
   mol->DeleteSelection(selHnd);
}


void
coot::restraints_container_t::make_strand_pseudo_bond_restraints() { 

   // This method of making pseudo bonds relies on the residue range
   // being continuous in sequence number (seqNum) and no insertion
   // codes messing up the number scheme.  If these are not the case
   // then we will make bonds between the wrongs atoms (of the wrong
   // residues)... a more sophisticated algorithm would be needed.
   // 
   // The method ignores residues with alt confs.
   //
   float pseudo_bond_esd = 0.08; // just a guess
   int selHnd = mol->NewSelection();
   int nSelResidues;
   mmdb::PPResidue SelResidue;
   mmdb::PPAtom res_1_atoms = NULL;
   mmdb::PPAtom res_2_atoms = NULL;
   mmdb::PPAtom res_3_atoms = NULL;
   int n_res_1_atoms;
   int n_res_2_atoms;
   int n_res_3_atoms;
   int index1 = -1; 
   int index2 = -1; 
   int index3 = -1; 
   mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"",  // altLocs
		mmdb::SKEY_NEW // selection key
		);
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 0) {
      for (int i=1; i<nSelResidues; i++) {
         // nO -> (n-1)O 4.64  // symmetric.  No need to specify both forwards
	 // Angle nO-(n+1)O-(n+2)O: 
	 SelResidue[i]->GetAtomTable(res_1_atoms, n_res_1_atoms);
	 if (res_1_atoms) { 
	    for (int iat1=0; iat1<n_res_1_atoms; iat1++) {
	       std::string at_1_name(res_1_atoms[iat1]->name);
	       // O Pseudo bonds and angles
	       if (at_1_name == " O  ") {
		  mmdb::Residue *contact_res = SelResidue[i-1];
		  if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 1)) {
		     contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		     if (res_2_atoms) { 
			for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
			   std::string at_2_name(res_2_atoms[iat2]->name);
			   if (at_2_name == " O  ") {
			      std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			      res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			      res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			      add(BOND_RESTRAINT, index1, index2, fixed_flags,
				  4.64, pseudo_bond_esd, 1.2);
			      std::cout << "Strand Bond restraint ("
					<< res_1_atoms[iat1]->name << " "
					<< res_1_atoms[iat1]->GetSeqNum() << ") to ("
					<< res_2_atoms[iat2]->name << " "
					<< res_2_atoms[iat2]->GetSeqNum() << ") 4.64" << std::endl;

			      // now the pseudo angle
			      if (i<(nSelResidues-1)) { 
				 mmdb::Residue *contact_res_2 = SelResidue[i+1];
				 if (SelResidue[i]->GetSeqNum() == (contact_res_2->GetSeqNum() - 1)) {
				    contact_res_2->GetAtomTable(res_3_atoms, n_res_3_atoms);
				    for (int iat3=0; iat3<n_res_3_atoms; iat3++) {
				       std::string at_3_name(res_3_atoms[iat3]->name);
				       if (at_3_name == " O  ") {
					  std::vector<bool> fixed_flag =
					     make_fixed_flags(index2, index1, index3);
					  res_3_atoms[iat3]->GetUDData(udd_atom_index_handle, index3);
					  // 98.0 degrees
					  add(ANGLE_RESTRAINT, index2, index1, index3,
					      fixed_flag, 98.0, 0.5, 1.2);
					  std::cout << "Strand Angle restraint ("
						    << res_1_atoms[iat1]->name << " "
						    << res_1_atoms[iat1]->GetSeqNum() << ") to ("
						    << res_2_atoms[iat2]->name << " "
						    << res_2_atoms[iat2]->GetSeqNum()
						    << ") to ("
						    << res_3_atoms[iat3]->name << " "
						    << res_3_atoms[iat3]->GetSeqNum()
						    << ") 98.0 " << std::endl;
					  break;
				       }
				    }
				 }
			      }
			      break;
			   }
			}
		     }
		  }
	       } // end of O

	       // Now make a CA-CA-CA pseudo angle of 120 degrees
	       if (at_1_name == " CA ") {
		  mmdb::Residue *contact_res_2 = SelResidue[i-1];
		  if (SelResidue[i]->GetSeqNum() == (contact_res_2->GetSeqNum() + 1)) {
		     contact_res_2->GetAtomTable(res_2_atoms, n_res_2_atoms);
		     if (res_2_atoms) { 
			for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
			   std::string at_2_name(res_2_atoms[iat2]->name);
			   if (at_2_name == " CA ") {
			      if (i<(nSelResidues-1)) {
				 mmdb::Residue *contact_res_3 = SelResidue[i+1];
				 if (SelResidue[i]->GetSeqNum() == (contact_res_3->GetSeqNum() - 1)) {
				    contact_res_3->GetAtomTable(res_3_atoms, n_res_3_atoms);
				    for (int iat3=0; iat3<n_res_3_atoms; iat3++) {
				       std::string at_3_name(res_3_atoms[iat3]->name);
				       if (at_3_name == " CA ") {
					  std::vector<bool> fixed_flag =
					     make_fixed_flags(index1, index2, index3);
					  res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
					  res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
					  res_3_atoms[iat3]->GetUDData(udd_atom_index_handle, index3);
					  add(ANGLE_RESTRAINT, index2, index1, index3,
					      fixed_flag, 120.0, 0.5, 1.2);
					  std::cout << "Strand Angle restraint ("
						    << res_1_atoms[iat1]->name << " "
						    << res_1_atoms[iat1]->GetSeqNum() << ") to ("
						    << res_2_atoms[iat2]->name << " "
						    << res_2_atoms[iat2]->GetSeqNum()
						    << ") to ("
						    << res_3_atoms[iat3]->name << " "
						    << res_3_atoms[iat3]->GetSeqNum()
						    << ") 120.0 " << std::endl;
					  break;
				       }
				    }
				 }
			      }
			      break;
			   }
			}
		     }
		  }
	       }
	    }
	 } 
      }
   }
   mol->DeleteSelection(selHnd);
}


int
coot::restraints_container_t::make_monomer_restraints(int imol,
						      const coot::protein_geometry &geom,
						      short int do_residue_internal_torsions) {

   // std::cout << "------------------------ in make_monomer_restraints() "
   // << from_residue_vector << std::endl;
   
   if (from_residue_vector)
      return make_monomer_restraints_from_res_vec(imol, geom, do_residue_internal_torsions);
   else
      return make_monomer_restraints_by_linear(imol, geom, do_residue_internal_torsions);

}

int
coot::restraints_container_t::make_monomer_restraints_by_linear(int imol,
								const coot::protein_geometry &geom,
								bool do_residue_internal_torsions) {

   // note: mini-rsr uses only the residue vector method
   
   int iret = 0;
   
   int selHnd = mol->NewSelection();
   int nSelResidues;
   coot::restraints_container_t::restraint_counts_t sum;

   mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		mmdb::SKEY_NEW // selection key
		);
   SelResidue_active = NULL;
   mol->GetSelIndex (selHnd, SelResidue_active, nSelResidues);
//    std::cout << "INFO:: GetSelIndex returned " << nSelResidues
// 	     << " residues (monomer restraints) " << std::endl;
   // save the (new (7Nov2003)) class variables (used in non_bonded
   // stuff) that keep the "active" (as opposed to "flanking") residues:
   nSelResidues_active = nSelResidues;
   // std::cout << "------------------------ in make_monomer_restraints_by_linear() nSelResidues "
   // << nSelResidues << std::endl;
   if (nSelResidues > 0) { 
      for (int i=0; i<nSelResidues; i++) {
	 if (SelResidue_active[i]) {
	    // std::cout << "------- calling make_monomer_restraints_by_residue() " << std::endl;
	    coot::restraints_container_t::restraint_counts_t local = 
	       make_monomer_restraints_by_residue(imol, SelResidue_active[i], geom,
						  do_residue_internal_torsions);
	    sum += local;
	 }
      }
   } else {
      std::cout << "get_monomer_restraints: There were no residues selected!? "
		<< std::endl;
   }
   // mol->DeleteSelection(selHnd); // -> makes crash! 6-feb-2004.
   //
   // This is because SelResidue_active is used elsewhere.
   // 
   // This should go into the destructor, I guess.

   sum.report(do_residue_internal_torsions);
   if (verbose_geometry_reporting != QUIET) {
      std::cout << "created " << restraints_vec.size() << " restraints" << std::endl;
      std::cout << std::endl;
   }
   return iret; // return 1 on success.  Hmm... how is this set? (and subsequently used?)
}

int
coot::restraints_container_t::make_monomer_restraints_from_res_vec(int imol,
								   const coot::protein_geometry &geom,
								   bool do_residue_internal_torsions) {

   bool print_summary = true;
   int iret = 0;

   coot::restraints_container_t::restraint_counts_t sum;

   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
      coot::restraints_container_t::restraint_counts_t local = 
	 make_monomer_restraints_by_residue(imol, residues_vec[ir].second, geom,
					    do_residue_internal_torsions);
      sum += local;
   } 

   if (verbose_geometry_reporting != QUIET) {
      std::cout << "INFO:: created " << restraints_vec.size() << " restraints" << std::endl;
      std::cout << std::endl;
      if (print_summary)
	 sum.report(do_residue_internal_torsions);
   }
   return iret;
}



coot::restraints_container_t::restraint_counts_t
coot::restraints_container_t::make_monomer_restraints_by_residue(int imol, mmdb::Residue *residue_p,
								 const protein_geometry &geom,
								 bool do_residue_internal_torsions) {

   coot::restraints_container_t::restraint_counts_t local;
   int i_no_res_atoms;
   mmdb::PPAtom res_selection = NULL;
   std::string pdb_resname(residue_p->name);
   if (pdb_resname == "UNK") pdb_resname = "ALA";

   if (false)
      std::cout << "--------------- make_monomer_restraints_by_residue() called "
		<< residue_spec_t(residue_p)
		<<  " and using type :" << pdb_resname << ": and imol "
		<< imol << std::endl;

   // idr: index dictionary residue
   int idr = geom.get_monomer_restraints_index(pdb_resname, imol, false);
   if (idr >= 0) {

      // if (geom[idr].comp_id == pdb_resname) {
      // old style comp_id usage
      // if (dictionary_name_matches_coords_resname(geom[idr].comp_id,pdb_resname)) {

      // OK, we need the 3 letter code for carbohydrates, the
      // comp_id for nucleotides:
      //
      // comp_id 3-letter-code name group
      // Ar         A        'Adenosine                    ' RNA                33  22 .
      // GAL-b-D    GAL      'beta_D_galactose             ' D-pyranose         24  12 .


      // now get a list of atoms in that residue
      // (SelResidue[i]) and compare them to the atoms in
      // geom[idr].bond_restraint[ib].

      residue_p->GetAtomTable(res_selection, i_no_res_atoms);
		  
      if (i_no_res_atoms > 0) {

	 if (restraints_usage_flag & BONDS_MASK)
	    local.n_bond_restraints += add_bonds(idr, res_selection, i_no_res_atoms,
						 residue_p, geom);
	    
	 if (restraints_usage_flag & ANGLES_MASK)
	    local.n_angle_restraints += add_angles(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);

	 if (restraints_usage_flag & TORSIONS_MASK) {
	    if (do_residue_internal_torsions) {
	       std::string residue_type = residue_p->GetResName();
	       if (residue_type != "PRO")
		  local.n_torsion_restr += add_torsions(idr, res_selection, i_no_res_atoms,
							residue_p, geom);
	    }
	 }

	 if (restraints_usage_flag & PLANES_MASK)
	    local.n_plane_restraints += add_planes(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);

	 if (restraints_usage_flag & CHIRAL_VOLUME_MASK)
	    local.n_chiral_restr += add_chirals(idr, res_selection, i_no_res_atoms, 
						residue_p, geom);

	 coot::restraints_container_t::restraint_counts_t mod_counts =
	    apply_mods(idr, res_selection, i_no_res_atoms, residue_p, geom);
	 // now combine mod_counts with local
      }
   }
   return local;
}


// return in millisecs
// double 
// coot::restraints_container_t::time_diff(const timeval &current, const timeval &start) const {

//    double d = current.tv_sec - start.tv_sec;
//    d *= 1000.0;
//    d += double(current.tv_usec - start.tv_usec)/1000.0;
//    return d;


// Need to add test that residues are linked with trans
//
std::vector<coot::rama_triple_t>
coot::restraints_container_t::make_rama_triples(int SelResHnd,
						const coot::protein_geometry &geom) const {
   std::vector<coot::rama_triple_t> v;
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);

   for (int i=0; i<(nSelResidues-2); i++) {
      if (SelResidue[i] && SelResidue[i+1] && SelResidue[i+2]) {
	 coot::rama_triple_t t(SelResidue[i], SelResidue[i+1], SelResidue[i+2]);
	 v.push_back(t);
      }
   }
   return v;
}


coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_conventional(int selHnd,
							   const coot::protein_geometry &geom) const {

   float dist_crit = 3.0; // if atoms in different residues are closer
			  // than this, then they are considered
			  // bonded (potentially).
   
   // First add "linear" links
   // 
   coot::bonded_pair_container_t c = bonded_residues_by_linear(selHnd, geom);

   // Now, are there any other links?
   //
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {
      for (int ii=0; ii<nSelResidues; ii++) {
	 for (int jj=0; jj<nSelResidues; jj++) {
	    if (jj>ii) {
	       if (! c.linked_already_p(SelResidue[ii], SelResidue[jj])) {
		  std::pair<bool, float> d = closest_approach(SelResidue[ii], SelResidue[jj]);
		  if (d.first) {
		     if (d.second < dist_crit) {
			std::pair<std::string, bool> l =
			   find_link_type_complicado(SelResidue[ii], SelResidue[jj], geom);
			if (l.first != "") {

			   // Eeek!  Fill me? [for 0.7]
			   
			} 
		     } 
		  }
	       }
	    }
	 }
      }
   }
   return c;
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_by_linear(int SelResHnd,
							const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t c;
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {
   
      std::string link_type("TRANS");
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage

      for (int i=0; i<(nSelResidues-1); i++) {
	 if (SelResidue[i] && SelResidue[i+1]) {


	    // There are a couple of ways we allow links.  First, that
	    // the residue numbers are consecutive.
	    //
	    // If there is a gap in residue numbering, the C and N
	    // have to be within 3A.
	    
	    // Are these residues neighbours?  We can add some sort of
	    // ins code test here in the future.
	    if ( abs(SelResidue[i]->GetSeqNum() - SelResidue[i+1]->GetSeqNum()) <= 1) { 
	       link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       // std::cout << "DEBUG ------------ in bonded_residues_by_linear() link_type is :"
	       // << link_type << ":" << std::endl;
	       if (link_type != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					whole_first_residue_is_fixed,
					whole_second_residue_is_fixed, link_type);
		  c.try_add(p);
	       }
	    }

	    // distance check this one.... it could be the opposite of
	    // an insertion code, or simply a gap - and we don't want
	    // to make a bond for a gap.
	    // 
	    if (abs(SelResidue[i]->index - SelResidue[i+1]->index) <= 1) {
	       // link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       std::pair<std::string, bool> link_info =
		  find_link_type_complicado(SelResidue[i], SelResidue[i+1], geom);
	       if (false)
		  std::cout << "DEBUG:: ---------- in bonded_residues_by_linear() link_info is :"
			    << link_info.first << " " << link_info.second << ":" << std::endl;
	       
	       if (link_info.first != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  if (link_info.second == 0) { 
		     coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_info.first);
		     c.try_add(p);
		  } else {
		     // order switch
		     coot::bonded_pair_t p(SelResidue[i+1], SelResidue[i],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_info.first);
		     c.try_add(p);
		  }
	       } 
	    }
	 }
      }
   }
   return c;
} 

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_from_res_vec(const coot::protein_geometry &geom) const {

   bool debug = false;

   coot::bonded_pair_container_t bpc;
   float dist_crit = 3.0;

   if (verbose_geometry_reporting == VERBOSE)
      debug = true;
   
   if (debug) {
      std::cout << "  debug:: bonded_residues_from_res_vec() residues_vec.size() " << residues_vec.size() << std::endl;
      for (unsigned int i=0; i<residues_vec.size(); i++) {
	 std::cout << "   " << residues_vec[i].first << " " << residue_spec_t(residues_vec[i].second) << std::endl;
      }
   }

   int nres = residues_vec.size();
   for (unsigned int ii=0; ii<residues_vec.size(); ii++) {
      mmdb::Residue *res_f = residues_vec[ii].second;
      for (unsigned int jj=ii+1; jj<residues_vec.size(); jj++) {
	 mmdb::Residue *res_s = residues_vec[jj].second;
	 std::pair<bool, float> d = closest_approach(res_f, res_s);

	 if (debug) { 
	    std::cout << " closest approach given " << coot::residue_spec_t(res_f)
		      << " and " << coot::residue_spec_t(res_s) << std::endl;
	    std::cout << " closest approach d " << d.first << " " << d.second << std::endl;
	 }
	 if (d.first) {
	    if (d.second < dist_crit) {
	       std::pair<std::string, bool> l  = find_link_type_complicado(res_f, res_s, geom);
	       std::string link_type = l.first;
	       if (link_type != "") {

		  // too verbose?
		  if (debug)
		     std::cout << "   INFO:: find_link_type_complicado(): "
			       << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
			       << " link_type -> :" << link_type << ":" << std::endl;
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  bool order_switch_flag = l.second;

		  if (!order_switch_flag) {
		     coot::bonded_pair_t p(res_f, res_s,
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_type);
		     bool added_flag = bpc.try_add(p);
		  } else {
		     coot::bonded_pair_t p(res_s, res_f,
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed,
					   link_type);
		     bool added_flag = bpc.try_add(p);
		  }
	       } else {
		  if (debug)
		     std::cout << "DEBUG:: blank link_type find_link_type_complicado() returns \""
			       << l.first << "\" " << l.second << std::endl;
	       } 
	    }
	 }
      }
   }
   bpc.filter();
   return bpc;
}



// a pair, first is if C and N are close and second if and order
// switch is needed to make it so.
std::pair<bool, bool>
coot::restraints_container_t::peptide_C_and_N_are_close_p(mmdb::Residue *r1, mmdb::Residue *r2) const {

   // needs PDBv3 fixup.
   float dist_crit = 2.8; // 20150714: this used to be 2.0.  That is too long, I think.
                          //           Try 2.8.
                          // 
                          // 2.0 A for a peptide link - so that we
			  // don't find unintentional peptides - which
			  // would be a disaster.

   std::string C_atom_name = " C  ";
   std::string N_atom_name = " N  ";
   
   mmdb::Atom *at_c_1 = NULL;
   mmdb::Atom *at_n_1 = NULL;
   mmdb::Atom *at_c_2 = NULL;
   mmdb::Atom *at_n_2 = NULL;

   mmdb::PPAtom residue_atoms_1 = NULL;
   mmdb::PPAtom residue_atoms_2 = NULL;
   int n_residue_atoms_1;
   int n_residue_atoms_2;
   r1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   r2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   for (int iat=0; iat<n_residue_atoms_1; iat++) {
      std::string atom_name(residue_atoms_1[iat]->name);
      if (atom_name == C_atom_name) {
	 at_c_1 = residue_atoms_1[iat];
      } 
      if (atom_name == N_atom_name) {
	 at_n_1 = residue_atoms_1[iat];
      } 
   }

   for (int iat=0; iat<n_residue_atoms_2; iat++) {
      std::string atom_name(residue_atoms_2[iat]->name);
      if (atom_name == C_atom_name) {
	 at_c_2 = residue_atoms_2[iat];
      } 
      if (atom_name == N_atom_name) {
	 at_n_2 = residue_atoms_2[iat];
      } 
   }

   if (at_c_1 && at_n_2) {
      clipper::Coord_orth c1(at_c_1->x, at_c_1->y, at_c_1->z);
      clipper::Coord_orth n2(at_n_2->x, at_n_2->y, at_n_2->z);
      float d = clipper::Coord_orth::length(c1, n2);
      // std::cout << "   C1->N2 dist " << d << std::endl;
      if (d < dist_crit)
	 return std::pair<bool, bool> (1,0);
   } 

   if (at_n_1 && at_c_2) {
      clipper::Coord_orth n1(at_n_1->x, at_n_1->y, at_n_1->z);
      clipper::Coord_orth c2(at_c_2->x, at_c_2->y, at_c_2->z);
      float d = clipper::Coord_orth::length(n1, c2);
      // std::cout << "   N1->C2 dist " << d << std::endl;
      if (d < dist_crit)
	 return std::pair<bool, bool> (1,1);
   } 

   return std::pair<bool, bool> (0, 0);

}



std::pair<bool,float>
coot::restraints_container_t::closest_approach(mmdb::Residue *r1, mmdb::Residue *r2) const {

   return coot::closest_approach(mol, r1, r2);
} 



// find residues in the neighbourhood that are not in the refining set
// and are not already marked as bonded flankers.
// 
void
coot::restraints_container_t::set_non_bonded_neighbour_residues_by_residue_vector(const coot::bonded_pair_container_t &bonded_flanking_pairs, const coot::protein_geometry &geom) {

   std::vector<mmdb::Residue *> nbr; // non-bonded residues 
   float dist_crit = 3.0;

   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
      std::vector<mmdb::Residue *> neighbours =
	 coot::residues_near_residue(residues_vec[ir].second, mol, dist_crit);
      for (unsigned int ineighb=0; ineighb<neighbours.size(); ineighb++) {
	 mmdb::Residue *test_res = neighbours[ineighb];
	 if (std::find(nbr.begin(), nbr.end(), test_res) == nbr.end()) {
	    // not already there...
	    bool found = 0;

	    if (0)
	       std::cout << ".... about to compare " << residue_spec_t(test_res) << " to "
			 << residues_vec.size() << " refining residues " << std::endl;
	    for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	       if (test_res == residues_vec[ires].second) {
		  found = 1;
		  break;
	       }
	    }

	    if (! found) {
	       // OK, so this neighbour was not in the passed set of
	       // moving residues (and not already in nbr)... it can
	       // be a flanking residue then...

	       // check that it is not a bonded flanking residue...
	       for (unsigned int iflank=0; iflank<bonded_flanking_pairs.size(); iflank++) { 
		  if (bonded_flanking_pairs[iflank].res_1 == test_res) {
		     found = 1;
		     // std::cout << "      oops bonded flanking residue res1 " << std::endl;
		     break;
		  } 
		  if (bonded_flanking_pairs[iflank].res_2 == test_res) {
		     found = 1;
		     // std::cout << "   oops bonded flanking residue res2 " << std::endl;
		     break;
		  }
	       }

	       if (! found) {
		  // std::cout << ".... adding non-bonded neighbour " << residue_spec_t(test_res) << std::endl;
		  nbr.push_back(test_res);
	       } 
	    }
	 }
      }
   }
   non_bonded_neighbour_residues = nbr;
} 


int 
coot::restraints_container_t::make_non_bonded_contact_restraints(int imol, const coot::bonded_pair_container_t &bpc,
								 const coot::protein_geometry &geom) {

   // is this function used any more?
   //
   coot::restraints_container_t::reduced_angle_info_container_t ai(restraints_vec);
   ai.write_angles_map("angles_map.tab");
   return make_non_bonded_contact_restraints(imol, bpc, ai, geom);
   
} 

// Atoms that are not involved in bonds or angles, but are in the
// residue selection should be at least 2.7A away from each other.
// 
// Here are my anti-bumping notes:
//
//
//     Anti-bumping restraints in regularization:
//
//     Considering totally screwed-up geometry: We should add a strong
//     repulsion for atoms that are not bonded so that they go away from
//     each other.  
//
//     Something like a triangle function between 0->2A and 0 beyond that.
//
//     Each atom has to check the distance to each other atom in the
//     selection: if they are not bonded, get a repulsion score for that
//     distance.
//
//     Derivative of that should be not too tricky, similar to bonds, but
//     not the same.
//
//     Instead of 500, use 10*matrix?  Doesn't really matter, I think.
//
//     Instead of using 2.0 as the critical distance, let's instead use
//     d_crit:
//
//     Infact,
//
//     f = 1000-d*(1000/d_crit)   for d<d_crit
//     df/dd = -1000/d_crit
//     df/dx = df/dd dd/dx
//           = -1000/d_crit
//
//     It's like bonds:
//     d = sqrt[ (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 ]
//     => dd/dx = (xi-xj)/d
//
//     So df/dx = -1000/d_crit * (xi-xj)/d
//
//     Need to keep a list of repulsing atom pairs so that we don't have
//     to calculate them once each for distortion_score and derivates..?
//
// Note that if residue-2 is not moving then it will not have angle restraints.  If
// it doesn't have angle resraints then the is_1_4_related test will fail.
// e.g (if n-1 is fixed residue): C(n-1)-N(n)-Ca(n)-C(n) or C(n-1)-N(n)-Ca(n)-CB(n)
// will not be seen as 1-4 related. So that's where strange_exception comes in.
//
int
coot::restraints_container_t::make_non_bonded_contact_restraints(int imol, const coot::bonded_pair_container_t &bpc,
								 const coot::restraints_container_t::reduced_angle_info_container_t &ai,
								 const coot::protein_geometry &geom) {

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;
   construct_non_bonded_contact_list(bpc, geom);

   // so now filtered_non_bonded_atom_indices is filled.
   // but it is not necessarily symmetric - so we can't do a j > 1 test (yet).
   // 
   // Write a debug/test for symmetry of filtered_non_bonded_atom_indices.
   symmetry_non_bonded_contacts(0);

   // We need to find if atom pairs are in the same ring.
   // We do that by finding the restraints of each residue and putting them in the map.
   // To make things faster in cases where the restraints look-up fails, we add a flag 
   // to the value of the map which let's us know that we have searched this dictionary 
   // type before.
   std::map<mmdb::Residue *, std::pair<bool, dictionary_residue_restraints_t> > restraints_map;

   if (false) {
      std::cout << "--------- make_non_bonded_contact_restraints() the atom array: " << std::endl;
      for (int iat=0; iat<n_atoms; iat++)
	 std::cout << "------- " << iat << " " << atom_spec_t(atom[iat]) << std::endl;
   }

   // Thinking of setting this to true? is the (link) angle in the dictionary? Is one of the
   // residues non-moving? (see above notes).
   if (false)
      ai.write_angles_map("angles-map.tab");
   
   if (false) {
      std::cout << "--------------------------------------------------\n";
      std::cout << "   non-bonded list:" << std::endl;
      std::cout << "--------------------------------------------------\n";
      for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
	 std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
	 for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) { 
	    std::cout << filtered_non_bonded_atom_indices[i][j] << " ";
	 } 
	 std::cout << std::endl;
      } 
      std::cout << "--------------------------------------------------\n";
   }

   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      mmdb::Atom *at = atom[i];
      std::string res_type = at->GetResName();
      std::map<mmdb::Residue *, std::pair<bool, dictionary_residue_restraints_t> >::const_iterator it;
      it = restraints_map.find(at->residue);
      if (it == restraints_map.end()) {
	 // have_restraints_for() is faster?
	 std::pair<bool, dictionary_residue_restraints_t> p = geom.get_monomer_restraints(res_type, imol);
	 // p.first is false if this is not a filled dictionary
	 restraints_map[at->residue] = p;
      }
   }

   // cache the energy types:
   std::map<mmdb::Atom *, std::string> energy_type_cache;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) {
      mmdb::Atom *at = atom[i];
      energy_type_cache[at] = get_type_energy(imol, at, geom);
   }

   int n_nbc_r = 0;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) {

	 mmdb::Atom *at_1 = atom[i];
	 mmdb::Atom *at_2 = atom[filtered_non_bonded_atom_indices[i][j]];

	 std::vector<bool> fixed_atom_flags =
	    make_non_bonded_fixed_flags(i, filtered_non_bonded_atom_indices[i][j]);

	 timeval start_time;
	 timeval current_time;
	 double d;
         if (at_1 && at_2) {

	    std::string type_1 = energy_type_cache[at_1];
	    std::string type_2 = energy_type_cache[at_2];

	    bool add_it = true;

	    // no H-H bumps in the same residue
	    //
	    // [20131212: Why not?  I suppose that there was a reason,
	    // it is not clear to me what it is now].  This needs to
	    // be investigated/fixed.
	    //
	    if (at_2->residue == at_1->residue)
	       if (is_hydrogen(at_1))
	          if (is_hydrogen(at_2))
		     add_it = false;

   	    if (filtered_non_bonded_atom_indices[i][j] < int(i))
  	       add_it = false;

	    int res_no_1 = at_1->GetSeqNum();
	    int res_no_2 = at_2->GetSeqNum();
	    
	    if (add_it) { 

	       // Don't make a bump between the CD of a PRO at residue(n) and the atoms of n-1
	    
	       std::string res_name_1 = at_1->GetResName();
	       std::string res_name_2 = at_2->GetResName();
	    
	       if (res_name_1 == "PRO" || res_name_1 == "HYP") {
		  int res_no_pro   = res_no_1;
		  int res_no_other = res_no_2;
		  if (res_no_pro == (res_no_other + 1)) {
		     std::string atom_name = at_1->name;
		     if (atom_name == " CD ") {  // PDBv3 FIXME
			add_it = false;
		     }
		  }
	       }
	       if (res_name_2 == "PRO" || res_name_2 == "HYP") {
		  int res_no_pro   = res_no_2;
		  int res_no_other = res_no_1;
		  if (res_no_pro == (res_no_other + 1)) {
		     std::string atom_name = at_2->name;
		     if (atom_name == " CD ") {  // PDBv3 FIXME
			add_it = false;
		     } 
		  }
	       }
	       // hack to remove C1-OD1 NBC on N-linked glycosylation
	       //
	       if (res_name_1 == "ASN" || res_name_2 == "NAG") {
		  std::string atom_name_1(at_1->name);
		  std::string atom_name_2(at_2->name);
		  if (atom_name_1 == " OD1")
		     if (atom_name_2 == " C1 ")
			add_it = false;
	       }
	       if (res_name_1 == "NAG" || res_name_2 == "ASN") {
		  std::string atom_name_1(at_1->name);
		  std::string atom_name_2(at_2->name);
		  if (atom_name_1 == " C1 ")
		     if (atom_name_2 == " OD1")
			add_it = false;
	       }
	    }

	    // -------------- OK add_it was set -----
	    
	    if (add_it) {

	       double dist_min = 3.4;

	       bool in_same_ring_flag    = true;
	       bool in_same_residue_flag = true;
	       
	       if (at_2->residue != at_1->residue) {
		  in_same_ring_flag    = false;
		  in_same_residue_flag = false;
	       }
	       
	       if (in_same_ring_flag) {
		  std::string atom_name_1 = at_1->GetAtomName();
		  std::string atom_name_2 = at_2->GetAtomName();

		  // in_same_ring_flag = restraints_map[at_2->residue].second.in_same_ring(atom_name_1,
		  //                                                                       atom_name_2);

		  in_same_ring_flag = is_in_same_ring(imol, at_2->residue,
						      residue_ring_map_cache,
						      atom_name_1, atom_name_2, geom);
	       }
	       
	       // this doesn't check 1-4 over a moving->non-moving peptide link (see comment above function)
	       // because the non-moving atom doesn't have angle restraints.
	       //
	       bool is_1_4_related = ai.is_1_4(i, filtered_non_bonded_atom_indices[i][j]);

	       if (false)
		  std::cout << "here C with at_1 " << atom_spec_t(at_1) << " at_2 " << atom_spec_t(at_2)
			    << " is_1_4_related " << is_1_4_related << std::endl;

	       if (is_1_4_related) {
		  dist_min = 2.64; // was 2.7 but c.f. guanine ring distances
		  if (is_hydrogen(at_1))
		      dist_min -= 0.7;
		  if (is_hydrogen(at_2))
		      dist_min -= 0.7;
	       } else {

		  std::pair<bool, double> nbc_dist = geom.get_nbc_dist(type_1, type_2,
								       in_same_residue_flag,
								       in_same_ring_flag);

		  if (nbc_dist.first) {

		     // In a helix O(n) is close to C(n+1), we should allow it.
		     // 
		     bool is_O_C_1_5_related = check_for_O_C_1_5_relation(at_1, at_2);

		     if (is_O_C_1_5_related) {
			dist_min = 2.84;
		     } else {

			// Perhaps we don't have angle restraints to both atoms because one
			// of the atoms is fixed (and thus miss that these have a 1-4 relationship).
			// e.g. O(n) [moving] -> CA(n+1) [fixed]
			// 
			// (this test will fail on insertion codes)
			//

			bool strange_exception = false;
			int rn_diff = abs(res_no_2 - res_no_1);
			if (rn_diff == 1) {
			   std::string atom_name_1 = at_1->GetAtomName();
			   std::string atom_name_2 = at_2->GetAtomName();
			   if (fixed_atom_flags.size()) {
			      if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
				 if (atom_name_1 == " O  ")
				    if (atom_name_2 == " CA ") 
				       strange_exception = true;
				 if (atom_name_1 == " CA ")
				    if (atom_name_2 == " O  ")
				       strange_exception = true;
				 if (atom_name_1 == " N  ")
				    if (atom_name_2 == " CB ")
				       strange_exception = true;
				 if (atom_name_1 == " CB ")
				    if (atom_name_2 == " N  ")
				       strange_exception = true;
				 if (atom_name_1 == " C  ")
				    if (atom_name_2 == " CB ")
				       strange_exception = true;
			      }
			   }
			   if (strange_exception)
			      dist_min = 2.7;

			   // Strange that these are not marked as 1-4 related.  Fix here...
			   // HA-CA-N-C can be down to ~2.4A.
			   // HA-CA-C-N can be down to ~2.41A.
			   if (res_no_2 > res_no_1) {
			      if (atom_name_1 == " C  ") {
				 if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			      if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
				 if (atom_name_2 == " N  ") {
				    strange_exception = true;
				    dist_min = 2.41;
				 }
			      }
			      if (atom_name_1 == " N  ") {
				 if (atom_name_2 == " H  ") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			   } else {
			      if (atom_name_1 == " HA " || atom_name_1 == "HA2" || atom_name_1 == " HA3") {
				 if (atom_name_2 == " C  ") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			      if (atom_name_1 == " N  ") {
				 if (atom_name_2 == " HA " || atom_name_2 == "HA2" || atom_name_2 == " HA3") {
				    strange_exception = true;
				    dist_min = 2.41;
				 }
			      }
			      if (atom_name_2 == " N  ") {
				 if (atom_name_1 == " H  ") {
				    strange_exception = true;
				    dist_min = 2.4;
				 }
			      }
			   }
			}
			if (rn_diff == 2) { 
			   if (fixed_atom_flags.size()) {
			      if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
				 std::string atom_name_1 = at_1->GetAtomName();
				 std::string atom_name_2 = at_2->GetAtomName();
				 if (atom_name_1 == " C  ")
				    if (atom_name_2 == " N  ")
				       strange_exception = true;
				 if (atom_name_1 == " N  ")
				    if (atom_name_2 == " C  ")
				       strange_exception = true; // 3.1 would be enough

				 if (strange_exception)
				    dist_min = 2.7;
			      }
			   }
			}

			if (! strange_exception)
			   dist_min = nbc_dist.second;
		     }
		  } else {
		     // short/standard value
		     dist_min = 2.8;
		  }
	       }

	 
	       if (false) { // debug.
	          clipper::Coord_orth pt1(atom[i]->x, atom[i]->y, atom[i]->z);
	          clipper::Coord_orth pt2(at_2->x,    at_2->y,    at_2->z);
	          double d = sqrt((pt1-pt2).lengthsq());

	          std::cout << "adding non-bonded contact restraint index " 
			    << i << " to index " << filtered_non_bonded_atom_indices[i][j]
			    << " "
			    << atom_spec_t(atom[i]) << " to " 
			    << atom_spec_t(atom[filtered_non_bonded_atom_indices[i][j]])
			    << "  types: " << type_1 <<  " " << type_2 <<  " fixed: "
			    << fixed_atom_flags[0] << " " << fixed_atom_flags[1] << "   current: " << d
			    << " dist_min: " << dist_min << std::endl;
	       }

	       if (is_hydrogen(at_1)) // should check from donor
		  if (is_acceptor(type_2, geom))
		     dist_min -= 0.7;
	       if (is_hydrogen(at_2)) // should check from donor
		  if (is_acceptor(type_1, geom))
		      dist_min -= 0.7;

	       simple_restraint r(NON_BONDED_CONTACT_RESTRAINT,
				  i, filtered_non_bonded_atom_indices[i][j],
				  type_1, type_2, 
				  fixed_atom_flags, dist_min);

	       restraints_vec.push_back(r);

	       n_nbc_r++;
	    }
	 }
      }
   }
   return n_nbc_r;
}

bool
coot::restraints_container_t::is_acceptor(const std::string &energy_type,
					  const coot::protein_geometry &geom) const {

   // get_energy_lib_atom() returns a blank atom on failure to look up energy_type
   energy_lib_atom ela = geom.get_energy_lib_atom(energy_type);
   bool acceptor_flag = ((ela.hb_type == HB_ACCEPTOR) || (ela.hb_type == HB_BOTH));
   
   return acceptor_flag;
}

// the bool in the residue_ring_map_cache is a flag that means "I've
// tried before to look this residue up and failed".
// 
bool
coot::restraints_container_t::is_in_same_ring(int imol, mmdb::Residue *residue_p,
					      std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > &residue_ring_map_cache,
					      const std::string &atom_name_1,
					      const std::string &atom_name_2,
					      const coot::protein_geometry &geom) const {
   bool r = false;

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map;
   std::list<std::string> r1;
   std::list<std::string> r2;
   std::list<std::string> r3;
   std::list<std::string> r4;

   // HIS
   r1.push_back(" CG ");
   r1.push_back(" CD2");
   r1.push_back(" ND1");
   r1.push_back(" CE1");
   r1.push_back(" NE2");

   // PHE/TYR
   r2.push_back(" CG ");
   r2.push_back(" CD1");
   r2.push_back(" CD2");
   r2.push_back(" CE1");
   r2.push_back(" CE2");
   r2.push_back(" CZ ");

   // TRP
   r3.push_back(" CG ");
   r3.push_back(" CD1");
   r3.push_back(" CD2");
   r3.push_back(" CE2");
   r3.push_back(" NE1");
   
   r4.push_back(" CD2");
   r4.push_back(" CE2");
   r4.push_back(" CE3");
   r4.push_back(" CZ2");
   r4.push_back(" CZ3");
   r4.push_back(" CH2");

   if (residue_ring_map_cache.size() == 0) {
      r1.sort();
      r2.sort();
      r3.sort();
      r4.sort();
      residue_ring_map["HIS"].second.push_back(r1);
      residue_ring_map["PHE"].second.push_back(r2);
      residue_ring_map["TYR"].second.push_back(r2);
      residue_ring_map["TRP"].second.push_back(r3);
      residue_ring_map["TRP"].second.push_back(r4);
      residue_ring_map["HIS"].first = false;
      residue_ring_map["PHE"].first = false;
      residue_ring_map["TYR"].first = false;
      residue_ring_map["TRP"].first = false;
   }

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > >::const_iterator it;
   std::string res_name = residue_p->GetResName();

   it = residue_ring_map_cache.find(res_name);
   if (it != residue_ring_map_cache.end()) {

      if (it->second.first == 0) { // not looked up before and failed 
	 for (unsigned int i=0; i<it->second.second.size(); i++) {
	    std::list<std::string>::const_iterator it_1 = std::find(it->second.second[i].begin(), it->second.second[i].end(), atom_name_1);
	    std::list<std::string>::const_iterator it_2 = std::find(it->second.second[i].begin(), it->second.second[i].end(), atom_name_2);
	    if (it_1 != it->second.second[i].end()) {
	       if (it_2 != it->second.second[i].end()) {
		  r = true;
		  break;
	       }
	    }
	 }
      } else {
	 // We tried to look it up before and failed
      }
   } else {

      // add it then
      std::pair<bool, dictionary_residue_restraints_t> rest =
	 geom.get_monomer_restraints(res_name, imol);
      if (rest.first) {
	 std::vector<std::vector<std::string> > ri = rest.second.get_ligand_ring_list();
	 residue_ring_map_cache[res_name].first = false; // not looked up before and failed	 
	 for (unsigned int ii=0; ii<ri.size(); ii++) {
	    std::list<std::string> l;
	    for (unsigned int jj=0; jj<ri[ii].size(); jj++)
	       l.push_back(ri[ii][jj]);
	    l.sort();
	    residue_ring_map_cache[res_name].second.push_back(l);
	 }

	 std::vector<std::list<std::string> > &vl = residue_ring_map_cache[res_name].second;

	 for (unsigned int ii=0; ii<vl.size(); ii++) {
	    std::list<std::string>::const_iterator it_1 = std::find(vl[ii].begin(), vl[ii].end(), atom_name_1);
	    std::list<std::string>::const_iterator it_2 = std::find(vl[ii].begin(), vl[ii].end(), atom_name_2);
	    if (it_1 != vl[ii].end()) {
	       if (it_2 != vl[ii].end()) {
		  r = true;
		  break;
	       }
	    }
	 }
      } else {
	 // OK, the lookup failed
	 std::vector<std::list<std::string> > fv;
	 std::pair<bool, std::vector<std::list<std::string> > > failed_data(true, fv);
	 residue_ring_map_cache[res_name] = failed_data;
      }
   }
   return r;
}


bool
coot::restraints_container_t::check_for_1_4_relation(int idx_1, int idx_2,
						     const reduced_angle_info_container_t &ai) const {

   bool is_1_4 = false;
   is_1_4 = ai.is_1_4(idx_1, idx_2);
   // std::cout << "debug:: check_for_1_4_relation(ai) " << idx_1 << " " << idx_2 << " is " << is_1_4
   // << std::endl;
   return is_1_4;
}

coot::restraints_container_t::reduced_angle_info_container_t::reduced_angle_info_container_t(const std::vector<coot::simple_restraint> &r) {

   // this map is constructed correctly.  If you are here it's because
   // you expect an angle restraint that not there.
   // 
   for (unsigned int ii=0; ii<r.size(); ii++) {
      if (r[ii].restraint_type == coot::ANGLE_RESTRAINT) {
	 std::pair<int, int> p_1(r[ii].atom_index_2, r[ii].atom_index_3);
	 std::pair<int, int> p_2(r[ii].atom_index_2, r[ii].atom_index_1);
	 angles[r[ii].atom_index_1].push_back(p_1);
	 angles[r[ii].atom_index_3].push_back(p_2);
      }
   }
}

void
coot::restraints_container_t::reduced_angle_info_container_t::write_angles_map(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   if (f) {
      std::map<int, std::vector<std::pair<int, int> > >::const_iterator it;
      for (it=angles.begin(); it!=angles.end(); it++) {
	 const std::vector<std::pair<int, int> > &v = it->second;
	 for (unsigned int i=0; i<v.size(); i++) {
	    f << "key: ";
	    f << it->first;
	    f << " value ";
	    f << " " << v[i].first <<  " " << v[i].second << "\n";
	 }
      }
      f.close();
   }

} 

bool
coot::restraints_container_t::reduced_angle_info_container_t::is_1_4(int indx_1, int indx_2) const {

   // this function can be const because we don't use [] operator on the angles map.
   
   bool f = false;

   std::map<int, std::vector<std::pair<int, int> > >::const_iterator it_1, it_2;
   it_1 = angles.find(indx_1);
   if (it_1 != angles.end()) {
      const std::vector<std::pair<int, int> > &v = it_1->second;
      for (unsigned int ii=0; ii<v.size(); ii++) {
	 
	 // what are the angles that have atom_mid as atom_1?  We can ask this because angles
	 // go into this object both way rounds: A-B-C, C-B-A.
	 
	 int idx_mid = v[ii].first;

	 it_2 = angles.find(idx_mid);
	 if (it_2 != angles.end()) {
	    const std::vector<std::pair<int, int> > &v_2 = it_2->second;
	    // are any of these indx_2?
	    for (unsigned int jj=0; jj<v_2.size(); jj++) { 
	       if (v_2[jj].second == indx_2) {
		  f = true;
		  break;
	       }
	    }
	 }

	 if (f)
	    break;

      }
   }

   return f;
} 

bool
coot::restraints_container_t::check_for_1_4_relation(int idx_1, int idx_2) const {

   bool is_1_4 = false;

   for (unsigned int ii=0; ii<restraints_vec.size(); ii++) { 
      if (restraints_vec[ii].restraint_type == coot::ANGLE_RESTRAINT) {

	 if (idx_1 == restraints_vec[ii].atom_index_1 ||
	     idx_1 == restraints_vec[ii].atom_index_3 ||
	     idx_2 == restraints_vec[ii].atom_index_1 ||
	     idx_2 == restraints_vec[ii].atom_index_3) { 

	    for (unsigned int jj=ii; jj<restraints_vec.size(); jj++) {
	       if (jj != ii) { 
		  if (restraints_vec[jj].restraint_type == coot::ANGLE_RESTRAINT) {

		     if (idx_2 == restraints_vec[jj].atom_index_1 ||
			 idx_2 == restraints_vec[jj].atom_index_3 ||
			 idx_1 == restraints_vec[jj].atom_index_1 ||
			 idx_1 == restraints_vec[jj].atom_index_3) {

			if (false)
			   std::cout << "check_for_1_4_relation() indices "
				     << idx_1 << " " << idx_2
				     << " examining angle restraint pair "
				     << restraints_vec[ii].atom_index_1 << " "
				     << restraints_vec[ii].atom_index_2 << " "
				     << restraints_vec[ii].atom_index_3 << " and "
				     << restraints_vec[jj].atom_index_1 << " "
				     << restraints_vec[jj].atom_index_2 << " "
				     << restraints_vec[jj].atom_index_3 << std::endl;

			if ((restraints_vec[ii].atom_index_2 == restraints_vec[jj].atom_index_1) ||
			    (restraints_vec[ii].atom_index_2 == restraints_vec[jj].atom_index_3)) {
			   
			   if ((restraints_vec[jj].atom_index_2 == restraints_vec[ii].atom_index_1) ||
			       (restraints_vec[jj].atom_index_2 == restraints_vec[ii].atom_index_3)) {
			      
			      is_1_4 = true;
			      break;
			   }
			} 
		     }
		  }
	       }
	    }
	 }
      }
      if (is_1_4)
	 break;
   }
   // std::cout << "debug:: check_for_1_4_relation() " << idx_1 << " " << idx_2 << " is " << is_1_4 << std::endl;
   return is_1_4;
}

// check either way round
bool
coot::restraints_container_t::check_for_O_C_1_5_relation(mmdb::Atom *at_1, mmdb::Atom *at_2) const {

   // PDBv3 FIXME.
   
   bool match = false;
   if (at_2->residue != at_1->residue) {

      // std::cout << "debug check_for_O_C_1_5_relation " << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << std::endl;

      // Check first at_1 is O(n) and at_2 is C(n+1)
      // 
      if ((at_1->GetSeqNum() + 1) == at_2->GetSeqNum()) {
	 std::string atom_name_1 = at_1->GetAtomName();
	 std::string atom_name_2 = at_2->GetAtomName();

	 if (atom_name_1 == " O  ") { 
	    if (atom_name_2 == " C  ") { 
	 
	       std::string chain_id_1 = at_1->GetChainID();
	       std::string chain_id_2 = at_2->GetChainID();
	       
	       if (chain_id_2 == chain_id_1) {
		  match = true;
	       } 
	    }
	 }
      }

      if (match) return match;

      // Check now that at_1 is C(n+1) and at_2 is O(n)
      // 
      if ((at_2->GetSeqNum() + 1) == at_1->GetSeqNum()) {
	 std::string atom_name_1 = at_1->GetAtomName();
	 std::string atom_name_2 = at_2->GetAtomName();

	 if (atom_name_1 == " C  ") { 
	    if (atom_name_2 == " O  ") { 
	 
	       std::string chain_id_1 = at_1->GetChainID();
	       std::string chain_id_2 = at_2->GetChainID();
	       
	       if (chain_id_2 == chain_id_1) {
		  match = true;
	       } 
	    }
	 }
      }
   }
   return match;
}


void
coot::restraints_container_t::symmetry_non_bonded_contacts(bool print_table) {

   int n_non_symmetric = 0;
   int n_ele = 0;
   int idx;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      n_ele += filtered_non_bonded_atom_indices[i].size();
      for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) {
	 idx = filtered_non_bonded_atom_indices[i][j];
	 // is i in the idx set?
	 if (std::find(filtered_non_bonded_atom_indices[idx].begin(),
		       filtered_non_bonded_atom_indices[idx].end(),
		       i) == filtered_non_bonded_atom_indices[idx].end()) {
	    // it wasn't - i.e. non-symmetry
	    if (0) {
	       std::cout << "   " << atom_spec_t(atom[idx]) << " was an unreciprocated neighbour of "
			 << atom_spec_t(atom[i]) << std::endl;
	       std::cout << "  to  " << idx << " added " << i << std::endl;
	    }
	    int prev_size = filtered_non_bonded_atom_indices[idx].size();
	    filtered_non_bonded_atom_indices[idx].push_back(i);
	    if (0)
	       std::cout << "  filtered_non_bonded_atom_indices[" << idx << "] was of size "
			 << prev_size << " and now " << filtered_non_bonded_atom_indices[idx].size()
			 << std::endl;
	    n_non_symmetric++;
	 }
      }
   }

   if (print_table) { 
      for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
	 std::cout << "  " << i << " : ";
	 for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++)
	    std::cout << " " << filtered_non_bonded_atom_indices[i][j];
	 std::cout << "\n";
      }
   }
} 


// fill the member data filtered_non_bonded_atom_indices
// 
void
coot::restraints_container_t::construct_non_bonded_contact_list(const coot::bonded_pair_container_t &bpc,
								const coot::protein_geometry &geom) {

   if (from_residue_vector)
      construct_non_bonded_contact_list_by_res_vec(bpc, geom);
   else 
      construct_non_bonded_contact_list_conventional();

}

std::pair<bool, double>
coot::simple_restraint::get_nbc_dist(const std::string &atom_1_type,
				     const std::string &atom_2_type, 
				     const protein_geometry &geom) {

   // This function reduces NBC distance for ring systems types and
   // H-B donor&acceptor combinations.  There is no provision of 1-4
   // distances and I dont know about 1-3 distances (which is a bit of
   // a worry).
   // 
   return geom.get_nbc_dist(atom_1_type, atom_2_type);
}

double
coot::simple_restraint::torsion_distortion(double model_theta) const {

   if ((restraint_type != TORSION_RESTRAINT) && (restraint_type != TRANS_PEPTIDE_RESTRAINT)) return 0;

   // this functions needs to mirror distortion_score_torsion()
   double diff = 99999.9; 
   double tdiff; 
   double trial_target; 
   int per = periodicity;
   for(int i=0; i<per; i++) { 
      // trial_target = torsion_restraint.target_value + double(i)*360.0/double(per);  ??
      trial_target = target_value + double(i)*360.0/double(per); 
      if (trial_target >= 360.0) trial_target -= 360.0; 
      tdiff = model_theta - trial_target;
      if (tdiff < -180) tdiff += 360;
      if (tdiff >  180) tdiff -= 360;
      if (abs(tdiff) < abs(diff)) { 
	 diff = tdiff;
      }
   }
   if (diff < -180.0) { 
      diff += 360.; 
   } else { 
      if (diff > 180.0) { 
	 diff -= 360.0; 
      }
   }
   if (false)
      std::cout << "in torsion_distortion() per: " << per << " diff: " << diff << " sigma: " << sigma
		<< "   penalty: " << diff*diff/(sigma * sigma) << std::endl;
   return diff*diff/(sigma * sigma);
} 


   
void
coot::restraints_container_t::construct_non_bonded_contact_list_conventional() {

   // So first, I need a method/list to determine what is not-bonded
   // to what.
   // 

//    std::cout << "bonded list:" << std::endl;
//    std::cout << "--------------------------------------------------\n";
//    for (int i=0; i<bonded_atom_indices.size(); i++) { 
//       std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
//       for (int j=0; j<bonded_atom_indices[i].size(); j++) { 
// 	 std::cout << bonded_atom_indices[i][j] << " ";
//       } 
//       std::cout << std::endl;
//    } 
//    std::cout << "--------------------------------------------------\n";

  // Now we need to know which indices into the mmdb::PPAtom atoms are in
  // the moving set (rather than the flanking atoms).
  // 
  std::vector<std::vector<int> > non_bonded_atom_indices;
  //
  short int was_bonded_flag;
  // Note: bonded_atom_indices is sized to n_atoms in init_shared_post().
  non_bonded_atom_indices.resize(bonded_atom_indices.size());

  // Set up some mmdb::PPAtom things needed in the loop:
  mmdb::PPAtom res_selection_local;
  int n_res_atoms;
  mmdb::PPAtom res_selection_local_inner;
  int n_res_atoms_inner;
  int atom_index, atom_index_inner;
  int ierr;

  // iar = i active residue, nSelResidues_active is class variable

  // std::cout << "INFO:: There are " << nSelResidues_active << " active residues\n";
  for (int iar=0; iar<nSelResidues_active; iar++) { 

     SelResidue_active[iar]->GetAtomTable(res_selection_local, n_res_atoms);
     // std::cout << "There are " << n_res_atoms << " active atoms in this active residue\n";

     for (int iat=0; iat<n_res_atoms; iat++) { 

	// set atom_index
	ierr = res_selection_local[iat]->GetUDData(udd_atom_index_handle, atom_index);
	if (ierr != mmdb::UDDATA_Ok) { 
	   std::cout << "ERROR:: in getting UDDATA res_selection_local, ierr=" 
		     << ierr << " "
		     << res_selection_local[iat]->GetSeqNum() << " " 
		     << res_selection_local[iat]->GetAtomName() << " \n";
	}
	
	bool matched_oxt = false;
	if (have_oxt_flag) {
	   if (std::string(res_selection_local[iat]->name) == " OXT") {  // PDBv3 FIXME
	      matched_oxt = true;
	   } else { 
	      matched_oxt = false;
	   }
	}

	if (! matched_oxt) { 

	   // For each of the bonds of atom with atom_index index
	   // we need to check if bonded_atom_indices[atom_index][j]
	   // matches any atom index of the active atoms:

	   for (int jar=0; jar<nSelResidues_active; jar++) { 
	   
	      SelResidue_active[jar]->GetAtomTable(res_selection_local_inner, 
						   n_res_atoms_inner);
	   
	      for (int jat=0; jat<n_res_atoms_inner; jat++) { 
	      
		 // set atom_index_inner
		 ierr =  res_selection_local_inner[jat]->GetUDData(udd_atom_index_handle, 
								   atom_index_inner);

		 if (atom_index == atom_index_inner) { 
		    // std::cout << "skipping same index " << std::endl;
		 } else {
		    
// 		    std::cout << "DEBUG:: checking bond pair " << atom_index << " " 
// 			      << atom_index_inner << " " 
// 			      << atom[atom_index]->name << " " << atom[atom_index]->GetSeqNum() << "    " 
// 			      << atom[atom_index_inner]->name << " " << atom[atom_index_inner]->GetSeqNum()
//          		      << std::endl;
	      
		    was_bonded_flag = 0;

		    // PDBv3 FIXME
		    if (have_oxt_flag) 
		       if (! strcmp(res_selection_local_inner[jat]->name, " OXT")) // matched
			  matched_oxt = true;

		    if (! matched_oxt) { 

		       for (unsigned int j=0; j<bonded_atom_indices[atom_index].size(); j++) { 
		 
			  if (bonded_atom_indices[atom_index][j] == atom_index_inner) { 
			     was_bonded_flag = 1;
			     break;
			  } 
		       }
		 
		       if (was_bonded_flag == 0) { 
			  non_bonded_atom_indices[atom_index].push_back(atom_index_inner);
		       }
		    }
		 }
	      }
	   }
	}
     }
  }

  if (false) {
     std::cout << "--------------------------------------------------\n";
     std::cout << "   conventional non-bonded list (unfiltered by distance):" << std::endl;
     std::cout << "--------------------------------------------------\n";
     for (unsigned int i=0; i<non_bonded_atom_indices.size(); i++) { 
	std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
	for (unsigned int j=0; j<non_bonded_atom_indices[i].size(); j++) { 
	   std::cout << non_bonded_atom_indices[i][j] << " ";
	} 
	std::cout << std::endl;
     } 
     std::cout << "--------------------------------------------------\n";
  }

  filter_non_bonded_by_distance(non_bonded_atom_indices, 8.0);
}

// This function is called by make_non_bonded_contact_restraints()
//
void
coot::restraints_container_t::construct_non_bonded_contact_list_by_res_vec(const coot::bonded_pair_container_t &bpc,
									   const coot::protein_geometry &geom) {

#ifdef HAVE_CXX_THREAD
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
#endif

   // How frequently does this function get called? - needs optimizing

   //  on a whole chain:
   //  8 -> 2.9 s
   // 11 -> 3.1 s
   //
   const double dist_crit = 11.0;
   
   filtered_non_bonded_atom_indices.resize(bonded_atom_indices.size());

   if (false) { // debug
      std::cout << "DEBUG:: construct_non_bonded_contact_list_by_res_vec ::::::::::::::::" << std::endl;
      for (unsigned int i=0; i<bpc.size(); i++)
	 std::cout << "   "
		   << coot::residue_spec_t(bpc[i].res_1) << " "
		   << coot::residue_spec_t(bpc[i].res_2) << " "
		   << bpc[i].is_fixed_first << " " 
		   << bpc[i].is_fixed_second << " " 
		   << std::endl;

      std::cout << "--------------- debug:: bonded_atom_indices size "
		<< bonded_atom_indices.size() << std::endl;
      std::cout << "--------------- debug:: n_atoms " << n_atoms << std::endl;

      std::cout << "Bonded atom indices:" << std::endl;
      for (unsigned int i=0; i<bonded_atom_indices.size(); i++) {
	 std::cout << "  " << i << " " << atom_spec_t(atom[i]) << " " << bonded_atom_indices[i].size()
		   << " |  ";
	 for (unsigned int j=0; j<bonded_atom_indices[i].size(); j++)
	    std::cout << " " << bonded_atom_indices[i][j];
	 std::cout << "\n";
      }
   }

   // Yes, this adds symmetry to filtered_non_bonded_atom_indices, i.e.
   // filtered_non_bonded_atom_indices[0] contains 1
   // filtered_non_bonded_atom_indices[1] contains 0
   // 
   for (unsigned int i=0; i<bonded_atom_indices.size(); i++) {

      // This is a hack.  It removes OXT from all NBCs.  Not the Right Way
      // 
      bool matched_oxt = false;
      if (have_oxt_flag) {
	 if (std::string(atom[i]->name) == " OXT") {  // PDBv3 FIXME
	    matched_oxt = true;
	 }
      }

      for (unsigned int j=0; j<bonded_atom_indices.size(); j++) {

	 if (i != j) {

	    if (have_oxt_flag) {
	       if (std::string(atom[j]->name) == " OXT") {  // PDBv3 FIXME
		  matched_oxt = true;
	       }
	    }

	    if (false)
	       std::cout << "moving->moving: here with atoms "
			 << atom_spec_t(atom[i]) <<  " " << atom_spec_t(atom[j])
			 << " have_oxt_flag: " << have_oxt_flag
			 << " matched_oxt: " << matched_oxt << std::endl;

	    if (! matched_oxt) {
	    
	       // In this section, we don't want NCBs within or to fixed
	       // residues (including the flanking residues), so if both
	       // atoms are in residues that are not in residue_vec, then
	       // we don't add a NCB for that atom pair.

	       // bonded_atom_indices contains indices of atoms that
	       // are angle-related (not just directly bonded)
	       // 
	       if (is_member_p(bonded_atom_indices[i], j)) {

		  // debug bonded atoms
		  
	       } else {
		  
		  // atom j is not bonded to atom i, is it close? (i.e. within dist_crit?)
		  clipper::Coord_orth pt1(atom[i]->x, atom[i]->y, atom[i]->z);
		  clipper::Coord_orth pt2(atom[j]->x, atom[j]->y, atom[j]->z);
		  double d = clipper::Coord_orth::length(pt1, pt2);
		  if (d < dist_crit) {
		     mmdb::Residue *r1 = atom[i]->residue;
		     mmdb::Residue *r2 = atom[j]->residue;

		     std::string alt_conf_1 = atom[i]->altLoc;
		     std::string alt_conf_2 = atom[j]->altLoc;

 		     if ((alt_conf_1 == alt_conf_2) ||
 			 (alt_conf_1.length() == 0) ||
 			 (alt_conf_2.length() == 0)) {

			if (is_a_moving_residue_p(r1) && is_a_moving_residue_p(r2)) {
			   filtered_non_bonded_atom_indices[i].push_back(j);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   // now add NBC restraints between atoms that are moving and atoms
   // of the neighbour residues.
   // 
   for (int iat=0; iat<n_atoms; iat++) {
      if (bonded_atom_indices[iat].size()) { 
	 mmdb::Residue *bonded_atom_residue = atom[iat]->residue;
	 for (int jat=0; jat<n_atoms; jat++) {

	    if (iat != jat) {
		     
	       mmdb::Residue *other_atom_residue = atom[jat]->residue;
	       if (bonded_atom_residue != other_atom_residue) {

		  if (is_a_moving_residue_p(bonded_atom_residue) &&
		      ! is_a_moving_residue_p(other_atom_residue)) {

		     bool matched_oxt = false;
		     if (have_oxt_flag) {
			if (std::string(atom[jat]->name) == " OXT") {  // PDBv3 FIXME
			   matched_oxt = true;
			} else { 
			   matched_oxt = false;
			}
		     }

		     if (false)
			std::cout << "moving->non-moving: here with atom "
				  << atom_spec_t(atom[jat]) << " have_oxt_flag: "
				  << have_oxt_flag << " matched_oxt: " << matched_oxt
				  << std::endl;

		     if (! matched_oxt) {

			bonded_pair_match_info_t mi = 
			   bpc.match_info(bonded_atom_residue, other_atom_residue);

			if (! mi.state) {
		      
			   // Simple part, the residues were not bonded to each other.
		     
			   if (! is_member_p(bonded_atom_indices[iat], jat)) {
			
			      // atom j is not bonded to atom i, is it close? (i.e. within dist_crit?)
			      clipper::Coord_orth pt1(atom[iat]->x, atom[iat]->y, atom[iat]->z);
			      clipper::Coord_orth pt2(atom[jat]->x, atom[jat]->y, atom[jat]->z);
			      double d = clipper::Coord_orth::length(pt1, pt2);
			      if (d < dist_crit) {
				 if (false)
				    std::cout << " ///////////////////////////// NBC   here "
					      << iat << " " << jat << " "
					      << coot::atom_spec_t(atom[iat]) << " "
					      << coot::atom_spec_t(atom[jat]) << " "
					      << std::endl;
				 filtered_non_bonded_atom_indices[iat].push_back(jat);
			      }
			   }
			   
			} else {

			   // they were bonded to each other.
			   
			   // add to filtered_non_bonded_atom_indices (which is a class variable)
			
			   if (mi.swap_needed)
			      construct_nbc_for_moving_non_moving_bonded(jat, iat, mi.link_type, geom);
			   else
			      construct_nbc_for_moving_non_moving_bonded(iat, jat, mi.link_type, geom);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (false) {
      std::cout << "--------------------------------------------------\n";
      std::cout << "  res-vec non-bonded list:" << std::endl;
      std::cout << "--------------------------------------------------\n";
      for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
	 std::cout << i << "  " << atom_spec_t(atom[i]) << " ";
	 for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) { 
	    std::cout << filtered_non_bonded_atom_indices[i][j] << " "
	       //  << atom_spec_t(atom[filtered_non_bonded_atom_indices[i][j]])
		      << " ";
	    if (j%20==0)
	       if (j > 0)
		  if (j != (filtered_non_bonded_atom_indices[i].size()-1))
		     std::cout << "\n          ";
	 } 
	 std::cout << std::endl;
      } 
      std::cout << "--------------------------------------------------\n";
   }
   
#ifdef HAVE_CXX_THREAD
   end = std::chrono::system_clock::now();

   std::chrono::duration<double> elapsed_seconds = end-start;
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);

   std::cout << "INFO:: nbc computation " // std::ctime(&end_time)
	     << "elapsed time: " << elapsed_seconds.count() << "s\n";
#endif // HAVE_CXX_THREAD

}

// Add non-bonded contacts for atoms that are in residues that are
// bonded to each other.  Atom iat is in a moving residue and atom jat
// is in a residue that is not moving.
// 
// Add to filtered_non_bonded_atom_indices (which is a class variable).
// 
void
coot::restraints_container_t::construct_nbc_for_moving_non_moving_bonded(unsigned int iat, unsigned int jat,
									 const std::string &link_type,
									 const coot::protein_geometry &geom) {

   // We dont know if res_1 is in the moving residue or not.  We do
   // know that res_1 and res_2 are in the correct order for the given
   // link_type link.
   // 
   mmdb::Residue *res_1 = atom[iat]->residue;
   mmdb::Residue *res_2 = atom[jat]->residue;

   dictionary_residue_link_restraints_t link = geom.link(link_type);
   // std::cout << "link: " << link.link_id << " " << link.link_bond_restraint.size() << std::endl;
   if (! link.empty()) {
      std::string atom_name_1 = atom[iat]->name;
      std::string atom_name_2 = atom[jat]->name;
      bool add_it = true;
      for (unsigned int i=0; i<link.link_bond_restraint.size(); i++) {
	 if (atom_name_1 == link.link_bond_restraint[i].atom_id_1_4c() && 
	     atom_name_2 == link.link_bond_restraint[i].atom_id_2_4c()) {
	    add_it = false;
	    break;
	 }
      }
      for (unsigned int i=0; i<link.link_angle_restraint.size(); i++) { 
	 if (atom_name_1 == link.link_angle_restraint[i].atom_id_1_4c() && 
	     atom_name_2 == link.link_angle_restraint[i].atom_id_3_4c()) {
	    add_it = false;
	    break;
	 }
      }
      for (unsigned int i=0; i<link.link_torsion_restraint.size(); i++) { 
	 if (atom_name_1 == link.link_torsion_restraint[i].atom_id_1_4c() && 
	     atom_name_2 == link.link_torsion_restraint[i].atom_id_4_4c()) {
	    add_it = false;
	    break;
	 }
      }
      if (add_it) {

	 filtered_non_bonded_atom_indices[iat].push_back(jat);
	 
	 if (0) { // debug.
	    clipper::Coord_orth pt1(atom[iat]->x, atom[iat]->y, atom[iat]->z);
	    clipper::Coord_orth pt2(atom[jat]->x, atom[jat]->y, atom[jat]->z);
	    double d = sqrt((pt1-pt2).lengthsq());
		     
	    std::cout << "moving-non-moving: adding filtered non-bonded atom indices: " 
		      << atom_spec_t(atom[iat]) << " to  " 
		      << atom_spec_t(atom[jat]) << " dist: " << d
		      << std::endl;
	 }
	 
      } else {

	 if (0) 
	    std::cout << "moving-non-moving: REJECT filtered non-bonded atom indices: " 
		      << atom_spec_t(atom[iat]) << " to  " 
		      << atom_spec_t(atom[jat]) 
		      << std::endl;
      } 
   }
}

void 
coot::restraints_container_t::filter_non_bonded_by_distance(const std::vector<std::vector<int> > &non_bonded_atom_indices, double dist) { 

   filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());

   mmdb::Atom *atom_1;
   mmdb::Atom *atom_2;
   double dist2;
   double dist_lim2 = dist*dist;
   int i_at_ind;
   
   for (unsigned int i=0; i<non_bonded_atom_indices.size(); i++) { 
      for (unsigned int j=0; j<non_bonded_atom_indices[i].size(); j++) {
	 
	 atom_1 = atom[i];
	 atom_2 = atom[non_bonded_atom_indices[i][j]];

// 	 dist2 = clipper::Coord_orth::lengthsq(clipper::Coord_orth(atom_1->x, atom_1->y, atom_1->z), 
// 					       clipper::Coord_orth(atom_2->x, atom_2->y, atom_2->z));

	 dist2 = (clipper::Coord_orth(atom_1->x, atom_1->y, atom_1->z) -
		  clipper::Coord_orth(atom_2->x, atom_2->y, atom_2->z)).lengthsq();

	 if (dist2 < dist_lim2) { 
// 	    std::cout << "accepting non-bonded contact between " << atom_1->GetSeqNum() 
// 		      << " " << atom_1->name << " and " << atom_2->GetSeqNum() 
// 		      << " " << atom_2->name  << "\n";
	    atom_2->GetUDData(udd_atom_index_handle, i_at_ind); // sets i_at_ind.
	    filtered_non_bonded_atom_indices[i].push_back(i_at_ind);
	 } else { 
// 	    std::cout << "          reject non-bonded contact between " << atom_1->GetSeqNum() 
// 		      << " " << atom_1->name  << " and " << atom_2->GetSeqNum() 
// 		      << " " << atom_2->name << " rejected by distance\n";
	 } 
      }
   }
}

bool
coot::restraints_container_t::is_a_moving_residue_p(mmdb::Residue *r) const {

   bool ret = 0;
   for (unsigned int i=0; i<residues_vec.size(); i++) {
      if (residues_vec[i].second == r) {
	 ret = 1;
	 break;
      } 
   }
   return ret;
}


int
coot::restraints_container_t::add_bonds(int idr, mmdb::PPAtom res_selection,
					int i_no_res_atoms,
					mmdb::PResidue SelRes,
					const coot::protein_geometry &geom) {

   int n_bond_restr = 0;
   int index1, index2;
   bool debug = false;

   if (debug)
      std::cout << "in add_bonds() for " << residue_spec_t(SelRes) << std::endl;

   for (unsigned int ib=0; ib<geom[idr].second.bond_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name1(res_selection[iat]->name);

	 if (debug)
	    std::cout << "comparing first (pdb) :" << pdb_atom_name1
		      << ": with (dict) :"
		      << geom[idr].second.bond_restraint[ib].atom_id_1_4c()
		      << ":" << std::endl; 

	 if (pdb_atom_name1 == geom[idr].second.bond_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       std::string pdb_atom_name2(res_selection[iat2]->name);

	       if (debug)
		  std::cout << "comparing second (pdb) :" << pdb_atom_name2
			    << ": with (dict) :"
			    << geom[idr].second.bond_restraint[ib].atom_id_2_4c()
			    << ":" << std::endl;
	       
	       if (pdb_atom_name2 == geom[idr].second.bond_restraint[ib].atom_id_2_4c()) {

		  // check that the alt confs aren't different
		  std::string alt_1(res_selection[iat ]->altLoc);
		  std::string alt_2(res_selection[iat2]->altLoc);
		  if (alt_1 == "" || alt_2 == "" || alt_1 == alt_2) { 

		     if (debug) { 
			std::cout << "atom match 1 " << pdb_atom_name1;
			std::cout << " atom match 2 " << pdb_atom_name2
				  << std::endl;
		     }

		     // now we need the indices of
		     // pdb_atom_name1 and
		     // pdb_atom_name2 in asc.atom_selection:

		     //  		  int index1_old = get_asc_index(pdb_atom_name1,
		     //  					     SelRes->seqNum,
		     //  					     SelRes->GetChainID());
		     //  		  int index2_old = get_asc_index(pdb_atom_name2,
		     //  					     SelRes->seqNum,
		     //  					     SelRes->GetChainID());

		     int udd_get_data_status_1 = res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
		     int udd_get_data_status_2 = res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);

		     // set the UDD flag for this residue being bonded/angle with 
		     // the other

		     if (udd_get_data_status_1 == mmdb::UDDATA_Ok &&
			 udd_get_data_status_2 == mmdb::UDDATA_Ok) { 
		  
			bonded_atom_indices[index1].push_back(index2);
			bonded_atom_indices[index2].push_back(index1);

			// this needs to be fixed for fixed atom (rather
			// than just knowing that these are not flanking
			// atoms).
			// 
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);

			if (false)
			   std::cout << "creating (monomer) bond restraint, idr " << idr
				     << " with fixed flags "
				     << fixed_flags[0] << " " << fixed_flags[1] << " "
				     << atom[index1]->GetSeqNum() << " "
				     << atom[index1]->name << " to "
				     << atom[index2]->GetSeqNum() << " "
				     << atom[index2]->name
				     << " restraint index " << n_bond_restr << "\n";
			try { 
			   add(BOND_RESTRAINT, index1, index2,
			       fixed_flags,
			       geom[idr].second.bond_restraint[ib].value_dist(),
			       geom[idr].second.bond_restraint[ib].value_esd(),
			       1.2);  // junk value
			   n_bond_restr++;
			}

			catch (const std::runtime_error &rte) {
			   
			   // do nothing, it's not really an error if the dictionary
			   // doesn't have target geometry (the bonding description came
			   // from a Chemical Component Dictionary entry for example).
			   std::cout << "trapped a runtime_error on adding bond restraint " << std::endl;
			} 
		     } else {
			std::cout << "ERROR:: Caught Enrico Stura bug.  How did it happen?" << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_bond_restr;
}

int
coot::restraints_container_t::add_angles(int idr, mmdb::PPAtom res_selection,
					 int i_no_res_atoms,
					 mmdb::PResidue SelRes,
					 const coot::protein_geometry &geom) {

   int n_angle_restr = 0;
   int index1, index2, index3;

//    std::cout << "There are " << geom[idr].angle_restraint.size()
// 	     << " angle restraints for this residue type" << std::endl; 

   for (unsigned int ib=0; ib<geom[idr].second.angle_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name1(res_selection[iat]->name);

//  	 std::cout << "angle:  comparing :" << pdb_atom_name1 << ": with :"
//  		   << geom[idr].angle_restraint[ib].atom_id_1_4c()
//  		   << ":" << std::endl;
	 
	 if (pdb_atom_name1 == geom[idr].second.angle_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       std::string pdb_atom_name2(res_selection[iat2]->name);
	       if (pdb_atom_name2 == geom[idr].second.angle_restraint[ib].atom_id_2_4c()) {
				    
// 		  std::cout << "angle: atom match 1 " << pdb_atom_name1;
// 		  std::cout << " atom match 2 " << pdb_atom_name2
// 			    << std::endl;

		  for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		     
		     std::string pdb_atom_name3(res_selection[iat3]->name);
		     if (pdb_atom_name3 == geom[idr].second.angle_restraint[ib].atom_id_3_4c()) {
		  

			// now we need the indices of
			// pdb_atom_name1 and
			// pdb_atom_name2 in asc.atom_selection:

//  			int index1_old = get_asc_index(pdb_atom_name1,
//  						   SelRes->seqNum,
// 						   SelRes->GetChainID());
//  			int index2_old = get_asc_index(pdb_atom_name2,
//  						   SelRes->seqNum,
//  						   SelRes->GetChainID());
//  			int index3_old = get_asc_index(pdb_atom_name3,
//  						   SelRes->seqNum,
//  						   SelRes->GetChainID());

			std::string alt_1(res_selection[iat ]->altLoc);
			std::string alt_2(res_selection[iat2]->altLoc);
			std::string alt_3(res_selection[iat3]->altLoc);

			if (((alt_1 == alt_2) && (alt_1 == alt_3)) ||
			    ((alt_1 == ""   ) && (alt_2 == alt_3)) ||
			    ((alt_2 == ""   ) && (alt_1 == alt_3)) ||
			    ((alt_3 == ""   ) && (alt_1 == alt_2)))
			   {
			
			   res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
			   res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
			   res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);

			   // std::cout << "add_angles: " << index1_old << " " << index1 << std::endl;
			   // std::cout << "add_angles: " << index2_old << " " << index2 << std::endl;
			   // std::cout << "add_angles: " << index3_old << " " << index3 << std::endl;
			
			   // set the UDD flag for this residue being bonded/angle with 
			   // the other
			
			   bonded_atom_indices[index1].push_back(index3);
			   bonded_atom_indices[index3].push_back(index1);
		  
			   // this needs to be fixed for fixed atom (rather
			   // than just knowing that these are not flanking
			   // atoms).
			   // 
			   std::vector<bool> fixed_flag = make_fixed_flags(index1, index2, index3);

			   add(ANGLE_RESTRAINT, index1, index2, index3,
			       fixed_flag,
			       geom[idr].second.angle_restraint[ib].angle(),
			       geom[idr].second.angle_restraint[ib].esd(),
			       1.2);  // junk value
			   n_angle_restr++;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_angle_restr;
}

int
coot::restraints_container_t::add_torsions(int idr, mmdb::PPAtom res_selection,
					   int i_no_res_atoms,
					   mmdb::PResidue SelRes,
					   const coot::protein_geometry &geom) {

   int n_torsion_restr = 0; 

   for (unsigned int ib=0; ib<geom[idr].second.torsion_restraint.size(); ib++) {

      // Joel Bard fix: Don't add torsion restraints for torsion that
      // have either s.d. or period 0

      if (geom[idr].second.torsion_restraint[ib].periodicity() > 0) { // we had this test most inner
	 if (geom[idr].second.torsion_restraint[ib].esd() > 0.000001) { // new test
	 
	    // now find the atoms
	    for (int iat=0; iat<i_no_res_atoms; iat++) {
	       std::string pdb_atom_name1(res_selection[iat]->name);

	       if (pdb_atom_name1 == geom[idr].second.torsion_restraint[ib].atom_id_1_4c()) {
		  for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

		     std::string pdb_atom_name2(res_selection[iat2]->name);
		     if (pdb_atom_name2 == geom[idr].second.torsion_restraint[ib].atom_id_2_4c()) {
				    
			// 		  std::cout << "atom match 1 " << pdb_atom_name1;
			// 		  std::cout << " atom match 2 " << pdb_atom_name2
			// 			    << std::endl;

			for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		     
			   std::string pdb_atom_name3(res_selection[iat3]->name);
			   if (pdb_atom_name3 == geom[idr].second.torsion_restraint[ib].atom_id_3_4c()) {
		  
			      for (int iat4=0; iat4<i_no_res_atoms; iat4++) {
		     
				 std::string pdb_atom_name4(res_selection[iat4]->name);
				 if (pdb_atom_name4 == geom[idr].second.torsion_restraint[ib].atom_id_4_4c()) {
		  
				    // now we need the indices of
				    // pdb_atom_name1 and
				    // pdb_atom_name2 in asc.atom_selection:

				    int index1;
				    int index2;
				    int index3;
				    int index4;

				    res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iat4]->GetUDData(udd_atom_index_handle, index4);

				    double torsion_angle = geom[idr].second.torsion_restraint[ib].angle();
				    if (torsion_angle < 0)
				       torsion_angle += 360;
				    if (torsion_angle > 360)
				       torsion_angle -= 360;

				    std::vector<bool> fixed_flags =
				       make_fixed_flags(index1, index2, index3, index4);
				    add(TORSION_RESTRAINT, index1, index2, index3, index4,
					fixed_flags,
					torsion_angle,
					geom[idr].second.torsion_restraint[ib].esd(),
					1.2,  // junk value
					geom[idr].second.torsion_restraint[ib].periodicity());
				    if (0) // debug
				       std::cout << "Adding monomer torsion restraint: "
						 << index1 << " "
						 << index2 << " "
						 << index3 << " "
						 << index4 << " angle "
						 << geom[idr].second.torsion_restraint[ib].angle() << " esd " 
						 << geom[idr].second.torsion_restraint[ib].esd() << " period " 
						 << geom[idr].second.torsion_restraint[ib].periodicity()
						 << std::endl;
				    n_torsion_restr++;
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   return n_torsion_restr;
}


int
coot::restraints_container_t::add_chirals(int idr, mmdb::PPAtom res_selection,
					  int i_no_res_atoms,
					  mmdb::PResidue SelRes,
					  const coot::protein_geometry &geom) { 

   int n_chiral_restr = 0;
   int index1, index2, index3, indexc;
   
   //   std::cout << "DEBUG:: trying to add chirals for this residue..." << std::endl;
   
   for (unsigned int ic=0; ic<geom[idr].second.chiral_restraint.size(); ic++) {
      // for now, let's just reject restraints that are a "both",
      // better would be to check the geometry and refine to the one
      // that is closest.
      if (!geom[idr].second.chiral_restraint[ic].is_a_both_restraint()) { 
	 for (int iat1=0; iat1<i_no_res_atoms; iat1++) {
	    std::string pdb_atom_name1(res_selection[iat1]->name);
	    if (pdb_atom_name1 == geom[idr].second.chiral_restraint[ic].atom_id_1_4c()) {
	       
	       for (int iat2=0; iat2<i_no_res_atoms; iat2++) {
		  std::string pdb_atom_name2(res_selection[iat2]->name);
		  if (pdb_atom_name2 == geom[idr].second.chiral_restraint[ic].atom_id_2_4c()) {
		     
		     for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
			std::string pdb_atom_name3(res_selection[iat3]->name);
			if (pdb_atom_name3 == geom[idr].second.chiral_restraint[ic].atom_id_3_4c()) {
			   
			   for (int iatc=0; iatc<i_no_res_atoms; iatc++) {
			      std::string pdb_atom_namec(res_selection[iatc]->name);
			      if (pdb_atom_namec == geom[idr].second.chiral_restraint[ic].atom_id_c_4c()) {
				 
//   			      std::cout << "DEBUG:: adding chiral number " << ic << " for " 
//   					<< res_selection[iatc]->GetSeqNum() << " "
//   					<< res_selection[0]->GetResName()
//   					<< pdb_atom_namec << " bonds to "
//   					<< pdb_atom_name1 << " "
//   					<< pdb_atom_name2 << " "
//   					<< pdb_atom_name3 << " "
//   					<< std::endl;

				 std::string alt_conf_c = res_selection[iatc]->altLoc;
				 std::string alt_conf_1 = res_selection[iat1]->altLoc;
				 std::string alt_conf_2 = res_selection[iat2]->altLoc;
				 std::string alt_conf_3 = res_selection[iat3]->altLoc;

				 if (((alt_conf_1 == alt_conf_c) || (alt_conf_1 == "")) &&
				     ((alt_conf_2 == alt_conf_c) || (alt_conf_2 == "")) && 
				     ((alt_conf_3 == alt_conf_c) || (alt_conf_3 == ""))) { 

				    res_selection[iat1]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iatc]->GetUDData(udd_atom_index_handle, indexc);


				    // does this chiral centre have exactly one hydrogen?
				    // If so, set chiral_hydrogen_index to the atom index.
				    // If not set to -1.
				    //
				    int chiral_hydrogen_index = get_chiral_hydrogen_index(indexc, index1, index2, index3);
				    
				    if (fabs(geom[idr].second.chiral_restraint[ic].target_volume()) < 1000.0 &&
					fabs(geom[idr].second.chiral_restraint[ic].target_volume()) > 0.00001) {

				       if (false) // debug
					  std::cout << "   Adding chiral restraint for "
						    << res_selection[iatc]->name
						    << " " << res_selection[iatc]->GetSeqNum() <<  " "
						    << res_selection[iatc]->GetChainID()
						    << " with target volume "
						    << geom[idr].second.chiral_restraint[ic].target_volume()
						    << " with volume sigma "
						    << geom[idr].second.chiral_restraint[ic].volume_sigma()
						    << " with volume sign "
						    << geom[idr].second.chiral_restraint[ic].volume_sign
						    << " idr index: " << idr << " ic index: " << ic
						    << " chiral_hydrogen_index: " << chiral_hydrogen_index
						    << std::endl;

				       std::vector<bool> fixed_flags =
					  make_fixed_flags(indexc, index1, index2, index3);
				       restraints_vec.push_back(simple_restraint(CHIRAL_VOLUME_RESTRAINT, indexc,
										 index1, index2, index3,
										 geom[idr].second.chiral_restraint[ic].volume_sign,
										 geom[idr].second.chiral_restraint[ic].target_volume(),
										 geom[idr].second.chiral_restraint[ic].volume_sigma(),
										 fixed_flags, chiral_hydrogen_index));
				       n_chiral_restr++;
				    } else {
				       std::cout << "WARNING:: Reject chiral restraint for "
						 << res_selection[iatc]->name
						 << " " << res_selection[iatc]->GetSeqNum() <<  " "
						 << res_selection[iatc]->GetChainID()
						 << " with target volume "
						 << geom[idr].second.chiral_restraint[ic].target_volume()
						 << " with volume sigma "
						 << geom[idr].second.chiral_restraint[ic].volume_sigma()
						 << " with volume sign "
						 << geom[idr].second.chiral_restraint[ic].volume_sign
						 << " idr index: " << idr << " ic index: " << ic
						 << " chiral_hydrogen_index: " << chiral_hydrogen_index
						 << std::endl;
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_chiral_restr;
}


// is there a single hydrogen connected to this chiral centre?
// If so, return the index, if not return -1
// 
int
coot::restraints_container_t::get_chiral_hydrogen_index(int indexc, int index1, int index2, int index3) const {

   int r = -1;
   int n_H = 0;
   int H_atom_index = -1;

   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & coot::BONDS_MASK) {
	 if ( (*this)[i].restraint_type == coot::BOND_RESTRAINT) {
	    mmdb::Atom *at_1 = atom[(*this)[i].atom_index_1]; 
	    mmdb::Atom *at_2 = atom[(*this)[i].atom_index_2];
	    if ((*this)[i].atom_index_1 == indexc) {
	       if (is_hydrogen(at_2)) {
		  H_atom_index = (*this)[i].atom_index_2;
		  n_H++;
	       } 
	    }
	    if ((*this)[i].atom_index_2 == indexc) {
	       if (is_hydrogen(at_1)) {
		  H_atom_index = (*this)[i].atom_index_1;
		  n_H++;
	       }
	    }
	 }
      }
   }
   if (n_H == 1)
      return H_atom_index;
   else
      return -1;
} 

// Creates any number of simple_restraints for this monomer and adds
// them to restraints_vec.
// 
// idr provides the index of the comp_id (e.g. "ALA") match in geom.
// 
int
coot::restraints_container_t::add_planes(int idr, mmdb::PPAtom res_selection,
					 int i_no_res_atoms,
					 mmdb::PResidue SelRes,
					 const coot::protein_geometry &geom) {

   bool debug = false;

   if (debug)
      std::cout << "There are " << geom[idr].second.plane_restraint.size()
		<< " dictionary plane restraints for " << SelRes->seqNum << " type: "
		<< geom[idr].second.residue_info.comp_id << std::endl;

   int n_plane_restr = 0;
   // either altconfs are all the same,
   // or they are different, in which case, add the atoms with blank altconfs to
   // (only) each of the non-blank ones
   //
   std::vector<std::string> altconfs = util::get_residue_alt_confs(SelRes);
   bool all_altconfs_the_same = true;
   if (altconfs.size() > 1)
      all_altconfs_the_same = false;

   for (unsigned int ip=0; ip<geom[idr].second.plane_restraint.size(); ip++) {
      std::map<std::string, std::vector <std::pair<int, double> > > idx_and_sigmas;
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name(res_selection[iat]->name);
	 std::string alt_conf(res_selection[iat]->altLoc);
	 for (int irest_at=0; irest_at<geom[idr].second.plane_restraint[ip].n_atoms(); irest_at++) {
	    if (pdb_atom_name == geom[idr].second.plane_restraint[ip].atom_id(irest_at)) {
	       // is this slow?
// 	       int idx = get_asc_index(res_selection[iat]->name,
// 				       res_selection[iat]->altLoc,
// 				       SelRes->seqNum,
// 				       SelRes->GetInsCode(),
// 				       SelRes->GetChainID());

	       int idx = get_asc_index(res_selection[iat]);

	       if (idx >= 0) {
		  double sigma = geom[idr].second.plane_restraint[ip].dist_esd(irest_at);
		  if (sigma > 0) {
		     std::pair<int, double> idx_sigma_pair(idx, sigma);
		     if (alt_conf.empty()) {
			if (all_altconfs_the_same) {
			   idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
			} else {
			   idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
			   for (unsigned int ialtconf=0; ialtconf<altconfs.size(); ialtconf++) {
			      if (! altconfs[ialtconf].empty()) {
				 idx_and_sigmas[altconfs[ialtconf]].push_back(idx_sigma_pair);
			      }
			   }
			}
		     } else {
			idx_and_sigmas[alt_conf].push_back(idx_sigma_pair);
		     }
		  }
	       }
	    }
	 }
      }

      std::map<std::string, std::vector <std::pair<int, double> > >::const_iterator it;
      for (it=idx_and_sigmas.begin(); it != idx_and_sigmas.end(); it++) {
	 if (it->second.size() > 3 ) {
	    std::vector<int> pos(it->second.size());
	    for (unsigned int i=0; i<it->second.size(); i++) pos[i] = it->second[i].first;
	    std::vector<bool> fixed_flags = make_fixed_flags(pos);
	    add_plane(it->second, fixed_flags);
	    n_plane_restr++;
	 }
      }

   }
   return n_plane_restr; 
}


// make RAMACHANDRAN_RESTRAINTs, not TORSION_RESTRAINTs these days.
int
coot::restraints_container_t::add_rama(std::string link_type,
				       mmdb::PResidue prev_res,
				       mmdb::PResidue this_res,
				       mmdb::PResidue post_res,
				       bool is_fixed_first,
				       bool is_fixed_second,
				       bool is_fixed_third,
				       const coot::protein_geometry &geom) {

   // Old notes:
   // TRANS    psi      1 N      1 CA     1 C      2 N   
   // TRANS    phi      1 C      2 N      2 CA     2 C   
   // TRANS    omega    1 CA     1 C      2 N      2 CA
   //
   // New assignements:
   // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
   // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C) 
   // 
   // So Rama_atoms in this order:
   //   0       1        2      3         4
   // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

   
   // std::cout << "DEBUG:: --------- :: Adding RAMA phi_psi_restraints_type" << std::endl;
   
   int n_rama = 0;
      
   mmdb::PPAtom prev_sel;
   mmdb::PPAtom this_sel;
   mmdb::PPAtom post_sel;
   int n_first_res_atoms, n_second_res_atoms, n_third_res_atoms;

   prev_res->GetAtomTable(prev_sel,  n_first_res_atoms); 
   this_res->GetAtomTable(this_sel, n_second_res_atoms);
   post_res->GetAtomTable(post_sel,  n_third_res_atoms);

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
      // throw 
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
      // throw
   }
   if (n_third_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
      // throw
   }

   std::vector<bool> fixed_flag(5, 0);
   if (is_fixed_first) {
      fixed_flag[0] = 1;
   }
   if (is_fixed_second) {
      fixed_flag[1] = 1;
      fixed_flag[2] = 1;
      fixed_flag[3] = 1;
   }
   if (is_fixed_third) {
      fixed_flag[4] = 1;
   }
   std::vector<mmdb::Atom *> rama_atoms(5);
   for (int ir=0; ir<5; ir++)
      rama_atoms[ir] = 0;
   
   for (int i=0; i<n_first_res_atoms; i++) {
      std::string atom_name(prev_sel[i]->name);
      if (atom_name == " C  ")
	 rama_atoms[0] = prev_sel[i];
   }
   for (int i=0; i<n_second_res_atoms; i++) {
      std::string atom_name(this_sel[i]->name);
      if (atom_name == " N  ")
	 rama_atoms[1] = this_sel[i];
      if (atom_name == " CA ")
	 rama_atoms[2] = this_sel[i];
      if (atom_name == " C  ")
	 rama_atoms[3] = this_sel[i];
   }
   for (int i=0; i<n_third_res_atoms; i++) {
      std::string atom_name(post_sel[i]->name);
      if (atom_name == " N  ")
	 rama_atoms[4] = post_sel[i];
   }

   if (rama_atoms[0] && rama_atoms[1] && rama_atoms[2] && 
       rama_atoms[3] && rama_atoms[4]) {

      std::vector<int> atom_indices(5, -1);
      for (int i=0; i<5; i++) {
// 	 atom_indices[i] = get_asc_index(rama_atoms[i]->name,
// 					 rama_atoms[i]->altLoc,
// 					 rama_atoms[i]->residue->seqNum,
// 					 rama_atoms[i]->GetInsCode(),
// 					 rama_atoms[i]->GetChainID());
	 atom_indices[i] = get_asc_index(rama_atoms[i]);
      }

      if ( (atom_indices[0] != -1) && (atom_indices[1] != -1) && (atom_indices[2] != -1) && 
	   (atom_indices[3] != -1) && (atom_indices[4] != -1)) { 

// 	 std::cout << "in add_rama() Adding RAMACHANDRAN_RESTRAINT\n       "
// 		   << coot::atom_spec_t(atom[atom_indices[0]]) << " " 
// 		   << coot::atom_spec_t(atom[atom_indices[1]]) << " " 
// 		   << coot::atom_spec_t(atom[atom_indices[2]]) << " " 
// 		   << coot::atom_spec_t(atom[atom_indices[3]]) << " " 
// 		   << coot::atom_spec_t(atom[atom_indices[4]]) << " fixed: "
// 		   << fixed_flag[0] << " " << fixed_flag[1] << " " 
// 		   << fixed_flag[2] << " " << fixed_flag[3] << " " 
// 		   << fixed_flag[4]
// 		   << std::endl;
	 
	 add(RAMACHANDRAN_RESTRAINT,
	     atom_indices[0], atom_indices[1], atom_indices[2],
	     atom_indices[3], atom_indices[4], fixed_flag);
	 n_rama++;
      }
   }
   // std::cout << "returning..." << n_rama << std::endl;
   return n_rama; 
}

coot::atom_spec_t
coot::restraints_container_t::get_atom_spec(int atom_index) const {

   if (atom)
      return atom_spec_t(atom[atom_index]);
   else
      return atom_spec_t();
}


int
coot::restraints_container_t::get_asc_index(const coot::atom_spec_t &spec) const {

   return get_asc_index_new(spec.atom_name.c_str(), spec.alt_conf.c_str(), spec.res_no,
			    spec.ins_code.c_str(), spec.chain_id.c_str());
}


int
coot::restraints_container_t::get_asc_index(const char *at_name,
					    const char *alt_loc,
					    int resno,
					    const char *ins_code,
					    const char *chain_id) const {

   return get_asc_index_new(at_name, alt_loc, resno, ins_code, chain_id);
   
}

int
coot::restraints_container_t::get_asc_index(mmdb::Atom *at) {

   int idx = -1;
   at->GetUDData(udd_atom_index_handle, idx);
   return idx;
}

int
coot::restraints_container_t::get_asc_index_new(const char *at_name,
						const char *alt_loc,
						int resno,
						const char *ins_code,
						const char *chain_id) const {

   int index = -1;

   if (mol) { 
      int SelHnd = mol->NewSelection(); // d
      mol->SelectAtoms(SelHnd,
		       0,
		       chain_id,
		       resno, ins_code,
		       resno, ins_code,
		       "*",      // resnames
		       at_name,  // anames
		       "*",      // elements
		       alt_loc  // altLocs 
		       );

      int nSelAtoms;
      mmdb::PPAtom SelAtom = NULL;
      mol->GetSelIndex(SelHnd, SelAtom, nSelAtoms);

      if (nSelAtoms > 0) {
	 if (udd_atom_index_handle >= 0) { 
	    SelAtom[0]->GetUDData(udd_atom_index_handle, index); // sets index
	 } else { 
	    index = get_asc_index_old(at_name, resno, chain_id);
	 } 
      }
      mol->DeleteSelection(SelHnd);
   }
   return index;
}

int
coot::restraints_container_t::get_asc_index_old(const std::string &at_name,
						int resno,
						const char *chain_id) const {

   int index = -1;
   int SelHnd = mol->NewSelection();
   
   mol->SelectAtoms(SelHnd,
			0,
			chain_id,
			resno, "*",
			resno, "*",
			"*", // rnames
			at_name.c_str(), // anames
			"*", // elements
			"*" // altLocs 
			);

   int nSelAtoms;
   mmdb::PPAtom SelAtom;
   mol->GetSelIndex(SelHnd, SelAtom, nSelAtoms);

   if (nSelAtoms > 0) {
      // now check indices.
      // 
      // Sigh. Is mmdb really this shit or am I missing something?
      //
      // It's not shit, you are not using it as it is supposed to be
      // used, I think. Instead of passing around the index to an
      // atom selection, you should simply be passing a pointer to an
      // atom.
      for (int i=0; i<n_atoms; i++) {
	 if (atom[i] == SelAtom[0]) {
	    index = i;
	    break;
	 }
      }
   }
   mol->DeleteSelection(SelHnd);

   if (index == -1 ) { 
      std::cout << "ERROR:: failed to find atom index for "
		<< at_name << " " << resno << " " << chain_id
		<< std::endl;
   } 
   return index; 
}
					    

// now setup the gsl_vector with initial values
//
// We presume that the atoms in mmdb::PPAtom are exactly the same 
// order as they are in the pdb file that refmac/libcheck uses
// to generate the restraints. 
//  
void 
coot::restraints_container_t::setup_gsl_vector_variables() {

   // recall that x is a class variable, 
   // (so are n_atoms and atom, which were set in the constructor)
   //  

   x = gsl_vector_alloc(3*n_atoms);

   // If atom is going out of date, check how atom is handled when the
   // restraints (this object) goes out of scope.  the destructor.
   
//    std::cout << "DEBUG:: using atom array pointer " << atom << std::endl;
//    std::cout << "DEBUG:: Top few atoms" << std::endl;
//    for (int i=0; i<10; i++) {
//       std::cout << "Top atom " << i << " "
// 		<< atom[i] << "   "
// 		<< atom[i]->x << " "
// 		<< atom[i]->y << " "
// 		<< atom[i]->z << " "
// 		<< std::endl;
//    } 

   for (int i=0; i<n_atoms; i++) {
      int idx = 3*i;
      gsl_vector_set(x, idx,   atom[i]->x);
      gsl_vector_set(x, idx+1, atom[i]->y);
      gsl_vector_set(x, idx+2, atom[i]->z);
   }

   setup_gsl_vector_atom_pos_deriv_locks();
}

void
coot::restraints_container_t::setup_gsl_vector_atom_pos_deriv_locks() {

   // setup gsl vector atom pos deriv locks.
   // We don't lock every derivative, we lock every atom (which corresponds to 3 derivs)
   //
#ifdef HAVE_CXX_THREAD

   // we need only do this once per instance gsl_vector_atom_pos_deriv_locks is set to 0 in init().
   //
   if (! gsl_vector_atom_pos_deriv_locks) {

      gsl_vector_atom_pos_deriv_locks = std::shared_ptr<std::atomic<unsigned int> > (new std::atomic<unsigned int>[n_atoms]);
      for (int ii=0; ii<n_atoms; ii++)
	 gsl_vector_atom_pos_deriv_locks.get()[ii] = 0; // unlocked
   }
#endif
}



void 
coot::restraints_container_t::update_atoms(gsl_vector *s) { 

   int idx;

   if (false) { 
      std::cout << "update_atom(0): from " << atom[0]->x  << " " << atom[0]->y << " " << atom[0]->z
		<< std::endl;
      double x = gsl_vector_get(s, 0);
      double y = gsl_vector_get(s, 1);
      double z = gsl_vector_get(s, 2);
      std::cout << "                  to " << x  << " " << y << " " << z << std::endl;
   }
   
   for (int i=0; i<n_atoms; i++) { 
      idx = 3*i; 
      atom[i]->x = gsl_vector_get(s,idx);
      atom[i]->y = gsl_vector_get(s,idx+1);
      atom[i]->z = gsl_vector_get(s,idx+2);
   }
} 


void
coot::restraints_container_t::position_OXT() { 

   if (oxt_reference_atom_pos.size()== 4) { 
      // std::cout << "DEBUG:: Positioning OXT by dictionary" << std::endl;
      double tors_o = 
	 clipper::Coord_orth::torsion(oxt_reference_atom_pos[0], 
				      oxt_reference_atom_pos[1], 
				      oxt_reference_atom_pos[2], 
				      oxt_reference_atom_pos[3]);
      double angl_o = clipper::Util::d2rad(120.8);
      clipper::Coord_orth oxt_pos(oxt_reference_atom_pos[0], 
				  oxt_reference_atom_pos[1], 
				  oxt_reference_atom_pos[2], 
				  1.231, angl_o, tors_o + M_PI);
      atom[oxt_index]->x = oxt_pos.x();
      atom[oxt_index]->y = oxt_pos.y();
      atom[oxt_index]->z = oxt_pos.z();
   }
} 

int
coot::restraints_container_t::write_new_atoms(std::string pdb_file_name) { 

   //
   int status = -1;
   if (mol != NULL) {
      // return 0 on success, non-zero on failure.
      status = mol->WritePDBASCII(pdb_file_name.c_str());
      if (status == 0)
	 std::cout << "INFO:: output file: " << pdb_file_name
		   << " written." << std::endl;
      else
	 std::cout << "WARNING:: output file: " << pdb_file_name
		   << " not written." << std::endl;
   } else { 
      cout << "not constructed from asc, not writing coords" << endl; 
   }
   return status;
}

void
coot::restraints_container_t::info() const {

   std::cout << "INFO:: There are " << restraints_vec.size() << " restraints" << std::endl;

   for (unsigned int i=0; i< restraints_vec.size(); i++) {
      if (restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) {
	 std::cout << "INFO:: restraint " << i << " is of type "
		   << restraints_vec[i].restraint_type << std::endl;

	 std::cout << restraints_vec[i].atom_index_1 << " "
		   << restraints_vec[i].atom_index_2 << " "
		   << restraints_vec[i].atom_index_3 << " "
		   << restraints_vec[i].atom_index_4 << " "
		   << restraints_vec[i].target_value << " "
		   << restraints_vec[i].sigma << " " << std::endl
		   << " with "
  		   << restraints_vec[i].plane_atom_index.size() << " vector atoms " << std::endl
		   << " with periodicity "
  		   << restraints_vec[i].periodicity << std::endl;
      }

      std::cout << "restraint number " << i << " is restraint_type " <<
	 restraints_vec[i].restraint_type << std::endl;
   }
} 

void
coot::simple_refine(mmdb::Residue *residue_p,
		    mmdb::Manager *mol,
		    const coot::dictionary_residue_restraints_t &dict_restraints) {

   if (residue_p) {
      if (mol) {

	 int imol = 0; // shouldn't matter
	 protein_geometry geom;
	 geom.replace_monomer_restraints(residue_p->GetResName(), imol, dict_restraints);
   
	 short int have_flanking_residue_at_start = 0;
	 short int have_flanking_residue_at_end = 0;
	 short int have_disulfide_residues = 0;
	 std::string altloc("");
	 std::vector<coot::atom_spec_t> fixed_atom_specs;

	 char *chain_id = residue_p->GetChainID();
	 int istart_res = residue_p->GetSeqNum();
	 int iend_res   = istart_res;

	 coot::restraints_container_t restraints(istart_res,
						 iend_res,
						 have_flanking_residue_at_start,
						 have_flanking_residue_at_end,
						 have_disulfide_residues,
						 altloc,
						 chain_id,
						 mol,
						 fixed_atom_specs);
   
	 // restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	 restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
	 pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	 bool do_internal_torsions = true;
	 bool do_trans_peptide_restraints = true;
	 restraints.make_restraints(imol, geom, flags, do_internal_torsions,
				    do_trans_peptide_restraints, 0, 0, pseudos);
	 restraints.minimize(flags, 3000, 1);
      }
   }
}



#endif // HAVE_GSL
