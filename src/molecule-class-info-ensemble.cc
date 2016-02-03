/* src/molecule-class-info-ensemble.cc
 * 
 * Copyright 2016
 * Author: Bernhard Lohkamp, Paul Emsley
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA.
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#ifdef _MSC_VER
#include <windows.h>
#endif

#include <stdexcept>
#include <algorithm>
#include <string.h>  // strncpy


#include <mmdb2/mmdb_manager.h>
#include "coords/Cartesian.h"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h"
#include "molecule-class-info.h"

#include "coot-utils/coot-map-utils.hh"
#include "xmap-utils.h"
#include "graphics-info.h"

// This is called by make_bonds_type_checked(), which is called by
// update_molecule_after_additions().
// 
void
molecule_class_info_t::update_ensemble_ghosts() {

   if (show_ensemble_ghosts_flag == 1) {
      if (ensemble_ghosts.size() > 0) {
         // FIXME (maybe) rather than updating each model ghost,
         // we just make model 0 (i.e. all) a ghost display
         ensemble_ghosts[0].update_ensemble_bonds(atom_sel.mol, 0);
      }
   }
}

// FIXME (maybe) not needed I think?!
//void
//molecule_class_info_t::delete_ensemble_ghost_selections() {

//   // 20060923 Bill Scott was reporting crashes in DeleteSelection
//   // here.  He was not using ghosts.  But had been manipulating his
//   // molecule a lot.
//   //
//   // So, I suspect that the selection handle for the ghosts was going
//   // out of date, which means that it's bad to delete it (of course).
//   //
//   // So, let's not delete selection if ghosts are not being used.  I
//   // use the is_empty() test (which seems to return 0 even if ghosts
//   // where not turned on! - oh dear?) and so also use the displayed?
//   // flag.
//   //
//   // Which means of course that the SelectionHandles of the NCS
//   // ghosts can be bogus because we are not updating them properly
//   // (and not only in this function - *any* function that changes the
//   // atom selection at all is suspect (moving atom coords/bfacs/occ
//   // is fine of course)).
//   //
//   // fill_ghost_info has a ncs_ghosts.resize(0), which is a potential
//   // memory leak, but that's not as bad as a crash (which would
//   // happen if we tried to DeleteSelection on SelectionHandles of out
//   // of date atom selections) is it?  For the safety checks to be
//   // removed all atom manipulation functions must be checked for
//   // proper operation with NCS ghosts atom selection.
//   //
//   // Question: why is is_empty() 0 for a not-turned-on ghost?
//   // (e.g. RNASA).
//   //
//   // Hmmm... reflection: NCS code is complex and crash-prone.

////    std::cout << "::::::::::::::::::::;; wwwwwoooooo!  ghosts! ::::::"
////  	     << std::endl;
   
//   if (ncs_ghosts.size() > 0) {
//      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
//// 	 std::cout << "Ghost " << ighost << " state "
//// 		   << ncs_ghosts[ighost].is_empty() << std::endl;
//	 if (! ncs_ghosts[ighost].is_empty()) {
//	    if (ncs_ghosts[ighost].display_it_flag) {
//	       atom_sel.mol->DeleteSelection(ncs_ghosts[ighost].SelectionHandle);
//	    }
//	 }
//      }
//   }
//}

// This is a ghost_molecule_display_t member function
void
coot::ghost_molecule_display_t::update_ensemble_bonds(mmdb::Manager *mol, int imod) {

   atom_selection_container_t asc;
   asc.mol = mol;

   // We should update the atom selection here: Yes, this needs to
   // happen.  Otherwise: a modification is made and when we come to
   // Bond_lines_container constructor below, which is given an atom
   // selection, it may well point to atoms that have been removed.
   // And then distaster [Bush was re-elected today].
   // 
   // (I guess that we don't need to regenerate the transformation
   // matrix)
   // 
   // (Note: trash the current selection first)
   //

   asc.atom_selection = NULL;
   SelectionHandle = mol->NewSelection();
   // std::cout << "ghost:: update_bonds new SelectionHandle: " << SelectionHandle << std::endl;

   mol->SelectAtoms(SelectionHandle, imod,
                    "*",
		    mmdb::ANY_RES, "*",
		    mmdb::ANY_RES, "*",
		    "*", "*", "*", "*");
   
   asc.mol->GetSelIndex(SelectionHandle, asc.atom_selection,
			asc.n_selected_atoms);
   asc.SelectionHandle = SelectionHandle;

   //    std::cout << "update_bonds ghost selection selected "
   // << asc.n_selected_atoms
   // << " atoms from chain " << chain_id << "\n";

   float min_dist = 0.1;
   float max_dist = 1.85;

   //    std::cout << "ghost molecule bonds molecule has " << asc.n_selected_atoms
   // 	     << " selected atoms" << std::endl;

   Bond_lines_container bonds(asc, min_dist, max_dist);
   bonds_box = bonds.make_graphical_bonds();

   // now run through bonds_box and change the coordinates therein by
   // rtop.  This is not apparently not a sensible place to put this
   // functionality because it should be part of the
   // Bond_lines_container, perhaps.  But Bond_lines_container does
   // not use clipper at all and I don't want to start now.

}

// public interface
int
molecule_class_info_t::update_ensemble_ghosts_size() {

   update_ensemble_ghosts();
   return ensemble_ghosts.size();
}

// Called on read pdb:
//
// fill ensemble ghosts:  detect models automatically.
int
molecule_class_info_t::fill_ensemble_info() {

   //nstd::cout << "DEBUG::   --------------- in fill_ghost_info ------- with homology_lev "
   // << homology_lev << std::endl;

   std::vector<std::string> chain_ids;
   std::vector<std::vector<std::pair<std::string, int> > > residue_types;
   std::vector<int> chain_atom_selection_handles;
   std::vector<short int> first_chain_of_this_type;

   bool allow_offset_flag = 0;
   if (is_from_shelx_ins_flag)
      allow_offset_flag = 1;

   // start from a blank slate:
   ensemble_ghosts.resize(0); // potential memory leak.

   if (atom_sel.n_selected_atoms > 0) {

      int n_models = atom_sel.mol->GetNumberOfModels();
      if (n_models > 0) {
         //cycle thru models now
         for (int imod=1; imod < n_models + 1; imod+=1){
            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
            mmdb::Chain *chain_p;
            int iselhnd = atom_sel.mol->NewSelection();
            mmdb::PAtom *atom_selection = NULL;
            int nSelAtoms;
            atom_sel.mol->SelectAtoms(iselhnd, imod,
                                      "*",
                                      mmdb::ANY_RES, "*",
                                      mmdb::ANY_RES, "*",
                                      "*", "*", "*", "*");
            atom_sel.mol->GetSelIndex(iselhnd, atom_selection, nSelAtoms);
            coot::ghost_molecule_display_t ghost;
            ghost.display_it_flag = 1;
            ghost.name = "Ensemble model ";
            ghost.name += coot::util::int_to_string(imod);
            // ghost.bonds_box filled by update_ghosts().
            ensemble_ghosts.push_back(ghost);
         }
      }

      if (ensemble_ghosts.size() > 0) {
         update_ensemble_ghosts();
         std::cout << "  INFO:: fill_ghost_info Constructed " << ensemble_ghosts.size() << " ghosts\n";
         for (unsigned int ighost=0; ighost<ensemble_ghosts.size(); ighost++) {
            std::cout << "      Ghost " << ighost << " name: \"" << ensemble_ghosts[ighost].name << "\""
                      << std::endl;
         }
      }
   }
   return ensemble_ghosts.size();
}

void
molecule_class_info_t::set_show_ensemble_ghosts(short int state) {

   show_ensemble_ghosts_flag = state;
   // caller redraws
}


void
molecule_class_info_t::set_ensemble_ghost_bond_thickness(float f) {

   ensemble_ghost_bond_width = f;
}


// FIXME (maybe) should have a test? Always good. (for what - without graphics?!)
//int
//molecule_class_info_t::test_function() {
//   int imol;
//   graphics_info_t g;

//   if (ncs_ghosts.size() > 0) {
//      if (ncs_ghosts_have_rtops_flag == 0) {
//	 float homology_lev =0.7;
//	 fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
//      }
//   }
//   std::cout <<  "make_dynamically_transformed_maps on " << ncs_ghosts.size()
//	     << " maps\n";

//   std::vector<coot::ghost_molecule_display_t> local_ncs_ghosts = ncs_ghosts;
//   int imol_base = graphics_info_t::n_molecules();
//   for(unsigned int ighost=0; ighost<10; ighost++) {
//      std::cout << "DEBUG:: pre-create molecule " << ighost << "/"
//		<< local_ncs_ghosts.size() << std::endl;
//      std::cout << "DEBUG:: This is imol=" << imol_no << std::endl;
//      imol = graphics_info_t::create_molecule();
//   }
   
//   imol = imol_base;
//   std::cout << "DEBUG:: pre-second-loop: This is imol=" << imol_no << std::endl;
//   for(unsigned int ighost=0; ighost<local_ncs_ghosts.size(); ighost++) {

//      std::cout << "DEBUG:: This is imol=" << imol_no << std::endl;
//      for (int itmp=0; itmp<=imol; itmp++)
//	 std::cout << "DEBUG:: molecule names: " << itmp << " :"
//		   << graphics_info_t::molecules[itmp].name_ << ":" << std::endl;
      
//      std::cout << "DEBUG:: NCS Copy to map number " << imol << std::endl;
//      std::cout << "DEBUG:: pre-install of ghost map " << ighost << "/"
//		<< local_ncs_ghosts.size() << std::endl;
//      std::cout << "DEBUG:: Post install of ghost map " << ighost << "/"
//	<< local_ncs_ghosts.size() << std::endl;
//   }

//   return imol;
//}


std::vector<coot::ghost_molecule_display_t>
molecule_class_info_t::get_ensemble_ghosts() const {

   return ensemble_ghosts;

}

// FIXME:: maybe somethgin to have!?

// not const because we can update (fill) ncs_ghosts.
//coot::ncs_differences_t
//molecule_class_info_t::ncs_chain_differences(std::string master_chain_id,
//					     float main_chain_weight) {

//   std::vector<coot::ncs_chain_difference_t> diffs;

//   // Note to self: recall:
//   //
//   // class ncs_chain_difference_t {
//   //    std::string peer_chain_id;
//   //    std::vector<ncs_residue_info_t> residue_info;

   
//   if (ncs_ghosts.size() > 0) {
//      if (ncs_ghosts_have_rtops_flag == 0) {
//	 float homology_lev =0.7;
//	 fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
//      }
//   }

//   if (ncs_ghosts.size() > 0) {
//      if (!ncs_ghosts_have_rtops_flag) {
//	 float homology_lev =0.7;
//	 fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
//      }
      
//      for (unsigned int ighost = 0; ighost<ncs_ghosts.size(); ighost++) {
//	 int imod = 1;
//	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
//	 mmdb::Chain *chain_p;
//	 // run over chains of the existing mol
//	 int nchains = model_p->GetNumberOfChains();
//	 mmdb::Chain *this_chain_p = 0;
//	 mmdb::Chain *master_chain_p = 0;
//	 for (int ichain=0; ichain<nchains; ichain++) {
//	    chain_p = model_p->GetChain(ichain);
//	    if (std::string(chain_p->GetChainID()) == ncs_ghosts[ighost].chain_id) {
//	       this_chain_p = chain_p;
//	    }
//	    if (std::string(chain_p->GetChainID()) == ncs_ghosts[ighost].target_chain_id) {
//	       master_chain_p = chain_p;
//	    }
//	 }
//	 if (this_chain_p && master_chain_p) {
	       
//	    int nres_this     = this_chain_p->GetNumberOfResidues();
//	    int nres_master = master_chain_p->GetNumberOfResidues();
//	    mmdb::PResidue this_residue_p;
//	    mmdb::PResidue master_residue_p;
//	    std::vector<coot::ncs_residue_info_t> residue_info;

//	    // return first == 0 if residues not found in chain
//	    std::pair<short int, int> mm_master =
//	       coot::util::min_resno_in_chain(master_chain_p);
//	    std::pair<short int, int> mm_this =
//	       coot::util::min_resno_in_chain(this_chain_p);

//	    if (mm_master.first && mm_this.first) {
//	       int resno_offset = mm_this.second - mm_master.second;

//	       for (int ires=0; ires<nres_this && ires<nres_master; ires++) {
//		  master_residue_p = master_chain_p->GetResidue(ires);
//		  this_residue_p = 0;
//		  if ( ((ires-resno_offset)<nres_this) && ((ires-resno_offset)>=0)) {
//		     this_residue_p = this_chain_p->GetResidue(master_residue_p->GetSeqNum(),
//							       master_residue_p->GetInsCode());
//		  }
//		  if (this_residue_p && master_residue_p) {
//		     if (this_residue_p->GetSeqNum() == master_residue_p->GetSeqNum()) {
//			coot::ncs_residue_info_t ds =
//			   ncs_ghosts[ighost].get_differences(this_residue_p, master_residue_p,
//							      main_chain_weight);
//			if (ds.filled) {
//			   if (0)
//			      std::cout << "     pushing back a residue_info with resno "
//					<< ds.resno << std::endl;
//			   residue_info.push_back(ds);
//			}
//		     }
//		  }
//	       }
//	    }
//	    coot::ncs_chain_difference_t d(ncs_ghosts[ighost].chain_id, residue_info);
//	    if (residue_info.size() > 0)
//	       diffs.push_back(d);
//	 }
//      }
//   }
//   return coot::ncs_differences_t(master_chain_id, diffs);
//}


//// Return a ncs_residue_info_t, note we cannot use that
//// ncs_residue_info_t if filled is 0.
////
//coot::ncs_residue_info_t
//coot::ghost_molecule_display_t::get_differences(mmdb::Residue *this_residue_p,
//						mmdb::Residue *master_residue_p,
//						float main_chain_weight) const {
//   // has access to rtop
//   coot::ncs_residue_info_t r;

//   if (std::string(this_residue_p->GetResName()) == std::string(master_residue_p->GetResName())) {
//      std::vector<std::pair<int, int> > index_pairs =
//	 coot::util::pair_residue_atoms(this_residue_p, master_residue_p);
//      mmdb::PPAtom residue_atoms_1 = NULL;
//      mmdb::PPAtom residue_atoms_2 = NULL;
//      int n_residue_atoms_1, n_residue_atoms_2;
//        this_residue_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
//      master_residue_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
//      float n_weighted_atoms = 0.0;
//      double sum_dist = 0.0;
//      for (unsigned int i=0; i<index_pairs.size(); i++) {
//	 float atom_weight = 1.0;
//	 mmdb::Atom *at1 = residue_atoms_1[index_pairs[i].first];
//	 mmdb::Atom *at2 = residue_atoms_2[index_pairs[i].second];
//	 std::string resname_1 = at1->GetResName(); // PDB puts waters in protein chains, bleugh.
//	 std::string resname_2 = at2->GetResName();
//	 if ((resname_1 != "HOH") && (resname_2 != "HOH")) {
//	    if (!at1->isTer() && !at2->isTer()) {
//	       if (coot::is_main_chain_p(at1))
//		  atom_weight = main_chain_weight;
//	       clipper::Coord_orth pt1(at1->x, at1->y, at1->z);
//	       clipper::Coord_orth pt2(at2->x, at2->y, at2->z);
//	       double len = clipper::Coord_orth::length(pt1.transform(rtop), pt2);
//	       sum_dist +=len;
//	       n_weighted_atoms += atom_weight;
//	    }
//	 }
//      }
//      if (n_weighted_atoms > 0.0) {
//	 r = coot::ncs_residue_info_t(this_residue_p->GetSeqNum(),
//				      this_residue_p->GetInsCode(),
//				      this_residue_p->index,
//				      master_residue_p->GetSeqNum(),
//				      master_residue_p->GetInsCode(),
//				      master_residue_p->index); // sets r.filled
//	 r.mean_diff = sum_dist/float(n_weighted_atoms);
//	 r.n_weighted_atoms = n_weighted_atoms;
//      }
//   } else {
//      std::cout << "different residue types " << this_residue_p->GetSeqNum()
//		<< " " << this_residue_p->GetInsCode() << " "
//		<< this_residue_p->GetResName() << "   vs  "
//		<< master_residue_p->GetSeqNum() << " " << master_residue_p->GetInsCode()
//		<< " " << master_residue_p->GetResName() << std::endl;
//   }
//   return r;
//}

