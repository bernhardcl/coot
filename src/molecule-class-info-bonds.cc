/* src/molecule-class-info-bonds.cc
 * 
 * Copyright 2016 by Medical Research Council
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

#include "molecule-class-info.h"


void
molecule_class_info_t::set_user_defined_colour_indices_by_residues(const std::vector<std::pair<coot::residue_spec_t, int> > &cis) {

   if (atom_sel.mol) {
      int udd_handle = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
      // std::cout << "got udd_handle: " << udd_handle << " " << cis.size() << " specs " << std::endl;
      if (udd_handle == 0)
	 udd_handle = atom_sel.mol->RegisterUDInteger(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
      
      for (unsigned int i=0; i<cis.size(); i++) {
	 const coot::residue_spec_t &spec = cis[i].first;
	 mmdb::Residue *residue_p = get_residue(spec);
	 if (residue_p) {
	    mmdb::Atom **residue_atoms = 0;
	    int n_residue_atoms;
	    residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       int ierr = at->PutUDData(udd_handle, cis[i].second);
	       if (ierr != mmdb::UDDATA_Ok) {
		  std::cout << "WARNING:: problem setting udd on atom " << coot::atom_spec_t(at) << std::endl;
	       }
	    }
	 } else {
	    std::cout << "WARNING:: No residue for " << spec << std::endl;
	 }
      }
   }
}


void
molecule_class_info_t::set_user_defined_colour_indices(const std::vector<std::pair<coot::atom_spec_t, int> > &cis) {

   if (atom_sel.mol) {
      int udd_handle = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
      if (udd_handle == 0)
	 udd_handle = atom_sel.mol->RegisterUDInteger(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
      for (unsigned int i=0; i<cis.size(); i++) {
	 const coot::atom_spec_t &spec = cis[i].first;
	 mmdb::Atom *at = get_atom(spec);
	 if (at) {
	    int ierr = at->PutUDData(udd_handle, cis[i].second);
	    if (ierr != mmdb::UDDATA_Ok) {
	       std::cout << "problem setting udd on atom " << coot::atom_spec_t(at) << std::endl;
	    }
	 }
      }
   }
}


void
molecule_class_info_t::user_defined_colours_representation(coot::protein_geometry *geom_p,
							   bool all_atoms_mode) { 

   bonds_box.clear_up();
   if (all_atoms_mode) {

      Bond_lines_container bonds(atom_sel, Bond_lines_container::COLOUR_BY_USER_DEFINED_COLOURS);
      bonds_box = bonds.make_graphical_bonds_no_thinning();
      bonds_box_type = coot::COLOUR_BY_USER_DEFINED_COLOURS_BONDS;
      
   } else {

      Bond_lines_container bonds(geom_p);
      bonds.do_Ca_plus_ligands_bonds(atom_sel, imol_no, geom_p, 2.4, 4.7, coot::COLOUR_BY_USER_DEFINED_COLOURS, false);
      bonds_box = bonds.make_graphical_bonds_no_thinning();
      bonds_box_type = coot::COLOUR_BY_USER_DEFINED_COLOURS_BONDS;
   }
}


void
molecule_class_info_t::clear_user_defined_atom_colours() {

   if (atom_sel.mol) {
      int udd_handle = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
      if (udd_handle != 0) {
	 udd_handle = 0; // reset
      }
   }
}
