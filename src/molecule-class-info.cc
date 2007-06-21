
/* src/molecule-class-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by Paul Emsley, The
 * University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#include <stdlib.h>

#ifndef _MSC_VER
#include <unistd.h>
#else
#include <windows.h>
#include <direct.h>
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#define snprintf _snprintf
#undef AddAtom
#undef GetAtomName
#endif

#include <iostream>
#include <string>
#include <vector>

#include "CIsoSurface.h"

#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/core/xmap.h"
#include "clipper/cns/cns_map_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats
#include "clipper/core/resol_basisfn.h"
#include "clipper/core/resol_targetfn.h"
#include "clipper/mmdb/clipper_mmdb.h"
#include "clipper/clipper-phs.h"
#include "clipper/contrib/sfcalc_obs.h"
#include "clipper/contrib/sfscale.h"
#include "clipper/contrib/sfweight.h"


#define HAVE_CIF

#ifdef HAVE_CIF
#include "clipper/clipper-cif.h"
#include "clipper/contrib/sfcalc.h"
#endif // HAVE_CIF

// using namespace clipper::float_data; old clipper version (before datatypes)
// 

using namespace std; // change me in a spare half hour when in York.

#include "mmdb_manager.h"
#include "mmdb_tables.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

// For stat, mkdir:
#include <sys/types.h>
#include <sys/stat.h>

#include "globjects.h"
#include "Bond_lines.h"

#include "gl-matrix.h"
#include "graphics-info.h"

#include "Bond_lines_ext.h"  
#include "graphical_skel.h"

#include "xmap-utils.h"


#include "coot-coord-utils.hh"
#include "coot-utils.hh"

#include <GL/glut.h> // needed (only?) for wirecube

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
// window magic jiggery pokery.
#define AddAtomA AddAtom
#define GetAtomNameA GetAtomName
#endif

//
// Return the molecule number of the molecule that we just filled.
// Return -1 if there was a failure.
// 
int
molecule_class_info_t::handle_read_draw_molecule(std::string filename,
						 short int reset_rotation_centre, 
						 short int is_undo_or_redo) {

   //
   graphics_info_t g;
   imol_no = g.n_molecules; // g.n_molecules gets updated outside afterwards

   if (! is_undo_or_redo) 
      bond_width = g.default_bond_width; // bleugh, perhaps this should
                                         // be a passed parameter?

   // std::cout << "DEBUG:: ---- imol_no is now " << imol_no << std::endl;

   // need to check that filename exists and is a file.
   //
   struct stat s;
   int status = stat(filename.c_str(), &s);

   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   }

				    
   // Read in pdb, [shelx files use the read_shelx_ins_file method]
   //
   atom_sel = get_atom_selection(filename);
   
   if ( atom_sel.read_success == 1 ) {

      //
      // and move mol_class_info to indexed molecule[n_molecules];
      // note indexing difficulties/considerations.
      
      // save it in the static molecule_class_info_t
      //
      
      // CMMDBCryst *cryst_p =  (atom_sel.mol)->get_cell_p();
      
      mat44 my_matt;

      // 
      // 
      int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (err != SYMOP_Ok) {
	 cout << "!! Warning:: No symmetry available for this molecule"
	      << endl;
      } else { 
	 cout << "Symmetry available for this molecule" << endl;
      }
      
      set_have_unit_cell_flag_maybe();

      if (molecule_is_all_c_alphas()) {
	 ca_representation();
      } else {

	 if (! is_undo_or_redo) { 
	    // 	 std::cout << "DEBUG:: filling ghost info in
	    // 	 handle_read_draw_molecule" << std::endl;
	    short int do_rtops_flag = 0;
	    // 0.7 is not used (I think) if do_rtops_flag is 0.
	    // hack to fix Mac bug/strangeness

	    // It only makes sense to to NCS chain searching on
	    // crystallographic models.  Which generally only have one
	    // model per manager.  This may change in future...
	    int nmodels = atom_sel.mol->GetNumberOfModels();
	    if (nmodels == 1) { 
	       int nghosts = fill_ghost_info(do_rtops_flag, 0.7);
	       // std::cout << "INFO:: found " << nghosts << " ghosts\n";
	    }
	 }
	    
	 // Generate bonds and save them in the graphical_bonds_container
	 // which has static data members.
	 //
	 if (bonds_box_type == coot::UNSET_TYPE)
	    bonds_box_type = coot::NORMAL_BONDS;
	 make_bonds_type_checked();
	 
      }

      //
      if (g.recentre_on_read_pdb || g.n_molecules == 0)  // n_molecules
							 // is updated
							 // in
							 // c-interface.cc
	 if (reset_rotation_centre) 
	    g.setRotationCentre(::centre_of_molecule(atom_sel)); 

      drawit = 1;
      if (g.show_symmetry == 1) {
	 if (show_symmetry) {  // internal
	    update_symmetry();
	 }
      }

      // Now, we have no map assocaited with this molecule, 
      // so we set the draw_vects to zero.
      // 
      // However, in future, we will have maps associated with
      // coordinates, so we should put a test here first before
      // we:
      //

      // initialize some things.
      //
      initialize_coordinate_things_on_read_molecule_internal(filename, is_undo_or_redo); 
      
      // update the maps so that they appear around the new centre. 
      // 
      for (int ii=0; ii<g.n_molecules; ii++) { 
	 g.molecules[ii].update_map(); 
      }

      // save state strings
      save_state_command_strings_.push_back("handle-read-draw-molecule");
      save_state_command_strings_.push_back(single_quote(coot::util::intelligent_debackslash(filename)));

      // We don't need any of this!  We have try_dynamic_add already (I'd forgotten!)
      // 
//       // 
//       std::vector<std::string> residues_types 
//               = coot::util::residues_in_molecule(atom_sel.mol);
//       for(int i=0; i<residues_types.size(); i++) { 
// 	 std::cout << i << " " << residues_types[i] << std::endl;
//       }
//       g.load_needed_monomers(residues_types);

      // we force a draw after we have updated the number of
      // molecules, n_molecules.

      return 1;

   } else {
      std::cout << "There was a coordinates read error\n";
      return -1;
   }
}

// cleaner interface to molecule's attributes:
std::pair<bool, clipper::Spacegroup>
molecule_class_info_t::space_group() const {

   clipper::Spacegroup sg;
   std::pair<bool, clipper::Spacegroup> p(0, sg);
   return p;
}

std::pair<bool, clipper::Cell>
molecule_class_info_t::cell() const {

   clipper::Cell cell;
   std::pair<bool, clipper::Cell> p(0, cell);
   return p;
}


coot::Cartesian
molecule_class_info_t::centre_of_molecule() const {

   double xs=0, ys=0, zs=0;
   if (atom_sel.n_selected_atoms > 0) {
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 xs += atom_sel.atom_selection[i]->x;
	 ys += atom_sel.atom_selection[i]->y;
	 zs += atom_sel.atom_selection[i]->z;
      }
      xs /= double(atom_sel.n_selected_atoms);
      ys /= double(atom_sel.n_selected_atoms);
      zs /= double(atom_sel.n_selected_atoms);
   }
   return coot::Cartesian(xs, ys, zs);
}


std::string
molecule_class_info_t::show_spacegroup() const { 

   std::string s("No spacegroup");
   
   if (has_model()) {
      char *st = atom_sel.mol->GetSpaceGroup();
      if (st)
	 s = st;
   }
   
   if (has_map()) 
      s = xmap_list[0].spacegroup().symbol_hm();

   return s;
}


coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt) const {

   coot::at_dist_info_t at_info = closest_atom(pt, 1);
   if (at_info.atom)
      return at_info;
   else
      return closest_atom(pt, 0);

}

coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt, bool ca_check_flag) const {

   coot::at_dist_info_t at_info(0,0,0);
   CAtom *at_best = 0;
   float dist_best = 99999999999.9;

   for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
      CAtom *at = atom_sel.atom_selection[iat];
      float d2 = (at->x - pt.x()) * (at->x - pt.x());
      d2 += (at->y - pt.y()) * (at->y - pt.y());
      d2 += (at->z - pt.z()) * (at->z - pt.z());
      if (d2 < dist_best) {
	 dist_best = d2;
	 at_best = at;
	 // Now, does this at belong to a residue that has a CA?  If
	 // it does, reset at_best to be the CA of the residue, but
	 // keep dist_best as it was, of course.
	 if (ca_check_flag == 1) {
	    CResidue *res = at->residue;
	    int natoms;
	    PPCAtom residue_atoms;
	    res->GetAtomTable(residue_atoms, natoms);
	    for (int iatom=0; iatom<natoms; iatom++) {
	       if (! strcmp(residue_atoms[iatom]->name, " CA ")) {
		  at_best = residue_atoms[iatom];
	       }
	    }
	 }
      }
   }
   if (at_best) { 
      at_info.dist = sqrt(dist_best);
      at_info.atom = at_best;
      at_info.imol = imol_no;
   }
   return at_info;
}




std::string
molecule_class_info_t::single_quote(const std::string &s) const {
   std::string r("\"");
   r += s;
   r += "\"";
   return r;
}

void
molecule_class_info_t::install_model(atom_selection_container_t asc,
				     const std::string &name, 
				     short int display_in_display_control_widget_status) {

   graphics_info_t g;
   imol_no = g.n_molecules;

   bond_width = g.default_bond_width; // bleugh, perhaps this should
				      // be a passed parameter?

   atom_sel = asc;

   CMMDBCryst *cryst_p =  (atom_sel.mol)->get_cell_p();
   mat44 my_matt;
   
   int err = cryst_p->GetTMatrix(my_matt, 0, 0, 0, 0);
   if (err != 0) {
      std::cout << "!! Warning:: No symmetry available for this molecule"
		<< std::endl;
   } else { 
      std::cout << "Symmetry available for this molecule" << std::endl;
   }
   set_have_unit_cell_flag_maybe();
   
   makebonds();
   if (g.show_symmetry == 1)
      if (show_symmetry) 
	 update_symmetry();

   have_unsaved_changes_flag = 1; 

   short int is_undo_or_redo = 0;

   if (display_in_display_control_widget_status == 0) {
      // treat as undo then (e.g. "terminal residue") made by add_cb_to_terminal_res().
      is_undo_or_redo = 1;
   } else {
      pickable_atom_selection = 1;
   }
   initialize_coordinate_things_on_read_molecule_internal(name, is_undo_or_redo);
}


void
molecule_class_info_t::update_map() {

   if (xmap_is_filled[0]) {
      coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
			 graphics_info_t::RotationCentre_y(),
			 graphics_info_t::RotationCentre_z());
      
      update_map_triangles(graphics_info_t::box_radius, rc); 
      if (graphics_info_t::display_lists_for_maps_flag) {
	 compile_density_map_display_list();
      }
   }
}

void
molecule_class_info_t::label_atoms(int brief_atom_labels_flag) {

   if (drawit) {
      
      if (has_model()) { 

	 // keep a list of atoms that have labels (either in graphics_info
	 // or mol_class_info) and loop over them calling label_atom(i)
	 // which labels the i'th atom of the atom selection in
	 // mol_class_info.
	 //
	 //
	 //
	 int max = max_labelled_atom();

	 if (max > 0) { 
	    // also remove labels from atom indexes list of over the end.
	    for (int ii=0; ii< max ; ii++)
	       label_atom(labelled_atom(ii), brief_atom_labels_flag);
	 }

	 max = max_labelled_symm_atom();

	 if (max > 0) { 
      
	    for (int ii=0; ii< max ; ii++) {

	       // symm_trans_t st = g.labelled_symm_atom_symm_trans(ii);

	       // cout << "symm label for atom: " <<  ii << " ppc index: "
	       //     << g.labelled_symm_atom(ii) << " " << st << endl;

	       //label_symm_atom(g.labelled_symm_atom(ii), st);
	       test_label_symm_atom(ii);
	    }
	 }
      }
   }
}

void
molecule_class_info_t::trim_atom_label_table() {

   int new_max = atom_sel.n_selected_atoms;
   int current_max = max_labelled_atom();
   std::vector<int> running_tab;
   std::vector<int> running_symm_tab;
   
   for (int i=0; i<current_max; i++) {
      if (i < new_max) { 
	 running_tab.push_back(labelled_atom_index_list[i]);
      }
   }
   n_labelled_atoms = running_tab.size();
   for (int i=0; i<n_labelled_atoms; i++) {
      labelled_atom_index_list[i] = running_tab[i];
   }

   // and now for symmetry index
   //
   current_max = max_labelled_symm_atom();
   for (int i=0; i<current_max; i++) {
      if (i < new_max) {
	 running_symm_tab.push_back(labelled_symm_atom_index_list[i]);
      }
   }
   n_labelled_symm_atoms = running_symm_tab.size();
   for (int i=0; i<n_labelled_symm_atoms; i++) {
      labelled_symm_atom_index_list[i] = running_symm_tab[i];
   }
}


void
molecule_class_info_t::anisotropic_atoms() {
   
   int c; // atom colour

   if (has_model()) { 
      graphics_info_t g;
      if (drawit) {
	 if (g.show_aniso_atoms_flag == 1 ) {
	    glPushMatrix();

	    float rx = g.X();
	    float ry = g.Y();
	    float rz = g.Z();

	    float x1, y1, z1;
	    float x_diff, y_diff, z_diff;
	    float d2, mc_r2 = g.show_aniso_atoms_radius*g.show_aniso_atoms_radius;
	    float rad_50, r;

	    for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      
	       // put a wiresphere at the atom positions
	 
	       if (atom_sel.atom_selection[i]->u11 > 0) {
	 
		  glLineWidth(1.0);
		  glPushMatrix();
	 
		  x1 = atom_sel.atom_selection[i]->x;
		  y1 = atom_sel.atom_selection[i]->y;
		  z1 = atom_sel.atom_selection[i]->z;

		  x_diff = x1 - rx;
		  y_diff = y1 - ry;
		  z_diff = z1 - rz;

		  d2 = x_diff*x_diff + y_diff*y_diff + z_diff*z_diff;

		  // are we either inside the distance or there is no distance set?
		  //
		  if ( (d2 <= mc_r2) || (g.show_aniso_atoms_radius_flag == 0) ) { 
	    
		     c = atom_colour(atom_sel.atom_selection[i]->element);
		     set_bond_colour_by_mol_no(c);
	       
		     GL_matrix mat(atom_sel.atom_selection[i]->u11, 
				   atom_sel.atom_selection[i]->u12,
				   atom_sel.atom_selection[i]->u13,
				   atom_sel.atom_selection[i]->u12, 
				   atom_sel.atom_selection[i]->u22,
				   atom_sel.atom_selection[i]->u23,
				   atom_sel.atom_selection[i]->u13, 
				   atom_sel.atom_selection[i]->u23,
				   atom_sel.atom_selection[i]->u33);

		     glTranslatef(x1, y1, z1);
		     // glMultMatrixf(mat.get());
		     glMultMatrixf(mat.cholesky().get());
		     rad_50 = r_50(atom_sel.atom_selection[i]->element);
		     r = rad_50_and_prob_to_radius(rad_50,
						   g.show_aniso_atoms_probability);
		     // note: g.show_aniso_atoms_probability is in the range
		     // 0.0 -> 100.0
		     glutWireSphere(r,10,10);

		  } 
		  glPopMatrix();
	       }
	    }
	    glPopMatrix();
	 }
      }
   }
}

// fix the name to something involving rotation perhaps?
//
void 
molecule_class_info_t::set_bond_colour_by_mol_no(int i) {

   if (bonds_rotate_colour_map_flag == 0) {
      set_bond_colour(i);
   } else {
      std::vector<float> rgb(3);
      // float rotation_size = float(imol_no + 1) * bonds_colour_map_rotation/360.0;
      float rotation_size = bonds_colour_map_rotation/360.0;

//       std::cout << " ::::::::::: in set_bond_colour bonds_colour_map_rotation is "
// 		<< bonds_colour_map_rotation << " for imol " << imol_no << std::endl;

      // rotation_size typically then: 2*32/360 = 0.178
      
      while (rotation_size > 1.0) { // no more black bonds?
	 rotation_size -= 1.0;
      } 
      switch (i) {
      case yellow: 
	  rgb[0] = 0.8; rgb[1] =  0.8; rgb[2] =  0.3;
	 break;
      case blue: 
	  rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] =  1.0;
	 break;
      case red: 
	  rgb[0] = 1.0; rgb[1] =  0.3; rgb[2] =  0.3;
	 break;
      case green:
	 rgb[0] = 0.1; rgb[1] =  0.99; rgb[2] =  0.1;
	 break;
      case grey: 
	 rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] =  0.7;
	 break;
// replaced in mmdb-extras.h
//       case white:   
// 	 rgb[0] = 0.99; rgb[1] =  0.99; rgb[2] = 0.99;
// 	 break;
      case magenta:
	 rgb[0] = 0.99; rgb[1] =  0.2; rgb[2] = 0.99;
	 break;
      case orange:
	 rgb[0] = 0.89; rgb[1] =  0.89; rgb[2] = 0.1;
	 break;
      case cyan:
	 rgb[0] = 0.1; rgb[1] =  0.89; rgb[2] = 0.89;
	 break;
      
      default:
	 rgb[0] = 0.8; rgb[1] =  0.2; rgb[2] =  0.2;
	 rgb = rotate_rgb(rgb, float(i*26.0/360.0));

      }

      // "correct" for the +1 added in the calculation of the rotation
      // size.
      // 21. is the default colour map rotation
      rgb = rotate_rgb(rgb, float(1.0 - 21.0/360.0));

      if (graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag) {
	 if (i == yellow) { 
	    std::vector<float> rgb_new = rotate_rgb(rgb, rotation_size);
	    bond_colour_internal = rgb_new;
	    glColor3f(rgb_new[0],rgb_new[1], rgb_new[2]);
	 } else {
	    bond_colour_internal = rgb;
	    glColor3f(rgb[0],rgb[1], rgb[2]);
	 }
      } else {
//  	 std::cout << "DEBUG: rotating coordinates colour map by "
//  		   << rotation_size * 360.0 << " degrees " << std::endl;
	 std::vector<float> rgb_new = rotate_rgb(rgb, rotation_size);
	 bond_colour_internal = rgb_new;
	 glColor3f(rgb_new[0], rgb_new[1], rgb_new[2]);
      }
   }
}


// aka rainbow - or maybe b factor, occupancy
void
molecule_class_info_t::set_bond_colour_by_colour_wheel_position(int i, int bond_type) {

   float max_colour = 30;
   std::vector<float> rgb(3);
   rgb[0] = 0.2; rgb[1] =  0.2; rgb[2] =  0.8; // blue
   
   // 30 is the size of rainbow colours, 0 -> 1.0 is the range of rainbow colours
   
   float rotation_size = 1.0 - float(i) * 0.7/max_colour;
   rgb = rotate_rgb(rgb, rotation_size);
   bond_colour_internal = rgb;
   glColor3f(rgb[0], rgb[1], rgb[2]);
}



// We find a box (symm: 2 trans: 0 0 0), but we don't find any atoms
// in it, but we still get though to bonds,
// bonds.make_graphical_symmetry_bonds
// 
// But there are no bonds generated.
// 
void
molecule_class_info_t::update_symmetry() {

   graphics_info_t g;

   // a bit of a hack...
   int shift_search_size = graphics_info_t::symmetry_shift_search_size;
      
   if ((graphics_info_t::show_symmetry == 1) &&
       (show_symmetry == 1)) {

      // don't do stuff until we have read in a molecule.
      //
      if (drawit == 1) {
      
	 molecule_extents_t extents(atom_sel, g.symmetry_search_radius);
	 graphics_info_t g;
	 coot::Cartesian point = g.RotationCentre();

      
	 // cout << "extents " << extents << endl;
	 // cout << "point:  " << point << endl;
	 std::vector<std::pair<symm_trans_t, Cell_Translation> > symm_trans_boxes =
	    extents.which_boxes(point, atom_sel, shift_search_size);

//   	 std::cout << "DEBUG:: symm_trans_boxes.size() is " 
//   		   << symm_trans_boxes.size() << std::endl;
//   	 std::cout << "Here are the symms we should check:" << std::endl;
//    	 for(int ii=0; ii<symm_trans_boxes.size(); ii++)
//   	    std::cout << ii << " " << symm_trans_boxes[ii].first << " "
//  		      << symm_trans_boxes[ii].second << std::endl;
	 

	 if (symm_trans_boxes.size() > 0) {

	    // when bonds goes out of scope (i.e. immediate after
	    // this) then the class data member vector "bonds" of the
	    // Bond_lines_container gets given back.
	    //
	    // It is with the "new"ly allocated graphical_symmetry_bonds
	    // that we need to concern ourselves.
	    //
	    Bond_lines_container bonds;

	    //
	    // delete the old symmetry_bonds_box
	    //
	    // symmetry_bonds_box.clear_up();
	    clear_up_all_symmetry();
 	    symmetry_bonds_box.resize(0);

// 	    symmetry_bonds_box.resize(0);
// 	    symmetry_bonds_box.push_back(bonds.addSymmetry(atom_sel,
// 						   point,
// 						   graphics_info_t::symmetry_search_radius,
// 						   symm_trans_boxes,
// 						   g.symmetry_as_calphas,
// 						   g.symmetry_whole_chain_flag));

// 	    for (unsigned int ibox=0; ibox<symm_trans_boxes.size(); ibox++)
// 	       std::cout << "box " << ibox << "/" << symm_trans_boxes.size()
// 			 << " " << symm_trans_boxes[ibox] << "\n";

	    symmetry_bonds_box = 
	       bonds.addSymmetry_vector_symms(atom_sel,
					      point,
					      graphics_info_t::symmetry_search_radius,
					      symm_trans_boxes,
					      symmetry_as_calphas,
					      symmetry_whole_chain_flag);
	    
	    //Bond_lines_container bonds(atom_sel,
	    //			       point,
	    //			       graphics_info_t::symmetry_search_radius,
	    //			       symm_trans);
	    
	 } else {
	    Bond_lines_container bonds(NO_SYMMETRY_BONDS);
	 }

	 if (show_strict_ncs_flag == 1) {
	    if (strict_ncs_matrices.size() > 0) {
	       update_strict_ncs_symmetry(point, extents);
	    }
	 }

      } else { 
	 // cout << "update_symmetry: no molecule yet" << endl;
      }
   }
}

//
void
molecule_class_info_t::draw_coord_unit_cell(const coot::colour_holder &cell_colour) {

   // Don't display if we have closed this molecule
   // (perhaps use (atom_sel.mol==NULL) instead?) (no).

   // (same test as has_model()):
   if (atom_sel.n_selected_atoms > 0) { 

      if (show_unit_cell_flag == 1) {

	 if (drawit) { 
	 
	    if (have_unit_cell == 1) {

	       glLineWidth(2.0);
	       glColor3f(cell_colour.red, cell_colour.green, cell_colour.blue);
	 
	       float corners[8][3] = {
		  {0,0,0}, //0 
		  {0,0,1}, //1 
		  {0,1,0}, //2 
		  {0,1,1}, //3 
		  {1,0,0}, //4 
		  {1,0,1}, //5 
		  {1,1,0}, //6 
		  {1,1,1}};//7 

	       realtype x_orth, y_orth, z_orth;
	       // rsc = real_space_corners
	       float rsc[8][3];

	       for (int ii=0; ii<8; ii++) {
	    
		  atom_sel.mol->Frac2Orth(corners[ii][0], corners[ii][1], corners[ii][2],
					  x_orth, y_orth, z_orth);
	    
		  rsc[ii][0] = x_orth;
		  rsc[ii][1] = y_orth;
		  rsc[ii][2] = z_orth;
	       }
	 
	       draw_unit_cell_internal(rsc);
	    }
	 }
      }
   }
}

//
void
molecule_class_info_t::draw_map_unit_cell(const coot::colour_holder &cell_colour) {

   if (has_map()) { 
      if (show_unit_cell_flag == 1) {
      
	 if ( max_xmaps > 0 ) {
	    if ( xmap_is_filled[0] ) {

	       if (drawit_for_map) { 
	    
		  // rsc = real_space_corners
		  float rsc[8][3];
	    
		  glLineWidth(2.0);
		  glColor3f(cell_colour.red, cell_colour.green, cell_colour.blue);
	    
		  float corners[8][3] = {
		     {0,0,0}, //0 
		     {0,0,1}, //1 
		     {0,1,0}, //2 
		     {0,1,1}, //3 
		     {1,0,0}, //4 
		     {1,0,1}, //5 
		     {1,1,0}, //6 
		     {1,1,1}};//7 
	    
		  for (int ii=0; ii<8; ii++) {
	       
		     clipper::Coord_frac c_f(corners[ii][0],corners[ii][1],corners[ii][2]);
	       
		     clipper::Coord_orth c_o = c_f.coord_orth( xmap_list[0].cell());
	       
		     rsc[ii][0] = c_o.x();
		     rsc[ii][1] = c_o.y();
		     rsc[ii][2] = c_o.z();
		  }
		  draw_unit_cell_internal(rsc); 
	       }
	    }
	 }
      }
   }
}

//
void
molecule_class_info_t::draw_unit_cell_internal(float rsc[8][3]) {

   std::vector<coot::CartesianPair> p;

   // bottom left connections
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[0][0], rsc[0][1], rsc[0][2]),
				   coot::Cartesian(rsc[1][0], rsc[1][1], rsc[1][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[0][0], rsc[0][1], rsc[0][2]),
				   coot::Cartesian(rsc[2][0], rsc[2][1], rsc[2][2])));
	 
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[0][0], rsc[0][1], rsc[0][2]),
				   coot::Cartesian(rsc[4][0], rsc[4][1], rsc[4][2])));

   // top right front connections
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[6][0], rsc[6][1], rsc[6][2]),
				   coot::Cartesian(rsc[4][0], rsc[4][1], rsc[4][2])));
	 
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[6][0], rsc[6][1], rsc[6][2]),
				   coot::Cartesian(rsc[2][0], rsc[2][1], rsc[2][2])));
	 
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[6][0], rsc[6][1], rsc[6][2]),
				   coot::Cartesian(rsc[7][0], rsc[7][1], rsc[7][2])));
	 
   // from 5
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[5][0], rsc[5][1], rsc[5][2]), 
				   coot::Cartesian(rsc[7][0], rsc[7][1], rsc[7][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[5][0], rsc[5][1], rsc[5][2]), 
				   coot::Cartesian(rsc[4][0], rsc[4][1], rsc[4][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[5][0], rsc[5][1], rsc[5][2]), 
				   coot::Cartesian(rsc[1][0], rsc[1][1], rsc[1][2])));

   // from 3
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[3][0], rsc[3][1], rsc[3][2]), 
				   coot::Cartesian(rsc[1][0], rsc[1][1], rsc[1][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[3][0], rsc[3][1], rsc[3][2]),
				   coot::Cartesian(rsc[7][0], rsc[7][1], rsc[7][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[3][0], rsc[3][1], rsc[3][2]), 
				   coot::Cartesian(rsc[2][0], rsc[2][1], rsc[2][2])));

   float x1, y1, z1;
   float x2, y2, z2;
   glBegin(GL_LINES);
   for (unsigned i=0; i<p.size(); i++) {
//       glVertex3f(p[i].getStart().x(),  p[i].getStart().y(),  p[i].getStart().z());
//       glVertex3f(p[i].getFinish().x(), p[i].getFinish().y(), p[i].getFinish().z());
      coot::Cartesian diff = p[i].getFinish() - p[i].getStart();
      for (float j=0.0; j<0.999; j+=0.1) {
	 x1 = p[i].getStart().x() + (j)*diff.x();
	 y1 = p[i].getStart().y() + (j)*diff.y();
	 z1 = p[i].getStart().z() + (j)*diff.z();
	 x2 = p[i].getStart().x() + (j+0.1)*diff.x();
	 y2 = p[i].getStart().y() + (j+0.1)*diff.y();
	 z2 = p[i].getStart().z() + (j+0.1)*diff.z();
	 
	 glVertex3f(x1, y1, z1);
	 glVertex3f(x2, y2, z2);
      }
   }
   glEnd();

	   
	 // add a label
// 	 glColor3f(1.0, 1.0, 1.0);
// 	 glColor3f(1.0, 0.2, 1.0);
	 glRasterPos3f(-1.6, -1.6,-1.6);      
	 printString("0");
	 glRasterPos3f(rsc[1][0]-1, rsc[1][1], rsc[1][2]+1);
	 printString("C");
	 glRasterPos3f(rsc[2][0]+1, rsc[2][1], rsc[2][2]+1);      
	 printString("B");
	 glRasterPos3f(rsc[4][0]+1, rsc[4][1]+1, rsc[4][2]-1);      
	 printString("A");


}

// --------------------------------------------------------------------
//   Conversion functions
// --------------------------------------------------------------------
//
void
molecule_class_info_t::initialize_coordinate_things_on_read_molecule(std::string molecule_name) {
   
   // presume not an undo/redo by default:
   initialize_coordinate_things_on_read_molecule_internal(molecule_name, 0);
}

// If we are a redo/undo, then we don't want to update (add a) mol in
// display control widget
// 
// Or a non graphics_info_t::molecules[] usage of this class.
// 
void
molecule_class_info_t::initialize_coordinate_things_on_read_molecule_internal(std::string molecule_name,
									      short int is_undo_or_redo) {

   // Atom lable initialization:

   labelled_atom_index_list = new int[coot::MAX_LABELLED_ATOMS];
   n_labelled_atoms = 0;

   labelled_symm_atom_index_list = new int[coot::MAX_LABELLED_ATOMS]; 
   n_labelled_symm_atoms = 0;

   // Dear oh dear, how anachronistic...
   labelled_symm_atom_symm_trans_ = new std::pair<symm_trans_t, Cell_Translation>[coot::MAX_LABELLED_ATOMS];

   // we use xmap_is_filled[0] to see if this molecule is a map
   //
   // FIXME.  Delete these lines, max_xmaps should be used instead.
   //
   xmap_is_filled = new int[10];
   xmap_is_filled[0] = 0;

   //
   name_ = molecule_name;

   // 
   drawit = 1; // by default, display it, we change change this later, if we want.

   //
   if (! is_undo_or_redo) { 
      bonds_colour_map_rotation = (imol_no + 1) * graphics_info_t::rotate_colour_map_on_read_pdb;
      while (bonds_colour_map_rotation > 360.0)
	 bonds_colour_map_rotation -= 360.0;
      bonds_rotate_colour_map_flag = graphics_info_t::rotate_colour_map_on_read_pdb_flag;
//       std::cout << "::::::: in initialization setting bonds_colour_map_rotation "
// 		<< bonds_colour_map_rotation << " for imol no " << imol_no << std::endl;
   }

   if (! is_undo_or_redo) { 
      // std::cout << "DEBUG:: not an undo/redo!\n";
      new_mol_in_display_control_widget(); // uses drawit
   }
}

void
molecule_class_info_t::set_symm_bond_colour_mol(int icol) {

   switch (icol) {
      case green:
	 glColor3f (combine_colour(0.1,0),
		    combine_colour(0.8,1),
		    combine_colour(0.1,2));
	 break;
      case blue: 
	 glColor3f (combine_colour(0.2,0),
		    combine_colour(0.2,1),
		    combine_colour(0.8,2));
	 break;
      case red: 
	 glColor3f (combine_colour(0.8,0),
		    combine_colour(0.1,1),
		    combine_colour(0.1,2));
	 break;
      case yellow: 
	 glColor3f (combine_colour(0.7,0),
		    combine_colour(0.7,1),
		    combine_colour(0.0,2));
	 break;
      
      default:
	 glColor3f (combine_colour(0.7, 0),
		    combine_colour(0.8, 1),
		    combine_colour(0.8, 2));
   }
}

void
molecule_class_info_t::set_symm_bond_colour_mol_and_symop(int icol, int isymop) {

//    std::cout << "in set_symm_bond_colour_mol_and_symop " << imol_no << " " << icol << " "
// 	     << isymop << " symmetry_rotate_colour_map_flag: "
// 	     << symmetry_rotate_colour_map_flag << "\n";
   
   if (symmetry_rotate_colour_map_flag) { 
      if (symmetry_colour_by_symop_flag) { 
	 set_symm_bond_colour_mol_rotate_colour_map(icol, isymop);
      } else { 
	 set_symm_bond_colour_mol_rotate_colour_map(icol, 0);
      }
   } else {
      set_symm_bond_colour_mol(icol);
   } 

}

void
molecule_class_info_t::set_symm_bond_colour_mol_rotate_colour_map(int icol, int isymop) {

   float step = graphics_info_t::rotate_colour_map_on_read_pdb/360.0;
   float rotation_size = float(icol + isymop) * step;
   std::vector<float> orig_colours(3);
   std::vector<float> rgb_new(3);
   std::vector<float> t_colours(3);

   switch (icol) {
   case green:
      t_colours[0] = combine_colour(0.1, 0);
      t_colours[1] = combine_colour(0.8, 1);
      t_colours[2] = combine_colour(0.1, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case blue: 
      t_colours[0] = combine_colour(0.2, 0);
      t_colours[1] = combine_colour(0.2, 1);
      t_colours[2] = combine_colour(0.8, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case red: 
      t_colours[0] = combine_colour(0.8, 0);
      t_colours[1] = combine_colour(0.1, 1);
      t_colours[2] = combine_colour(0.1, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case yellow: 
      t_colours[0] = combine_colour(0.7, 0);
      t_colours[1] = combine_colour(0.7, 1);
      t_colours[2] = combine_colour(0.0, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
      
   default:
      t_colours[0] = combine_colour(0.6, 0);
      t_colours[1] = combine_colour(0.7, 1);
      t_colours[2] = combine_colour(0.7, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
   }
}



float
molecule_class_info_t::combine_colour(float v, int col_part_index) {

   // col_part_index is 0,1,2 for red gree blue components of the colour
   double w = graphics_info_t::symm_colour_merge_weight[0];
   return w*graphics_info_t::symm_colour[0][col_part_index] + v*(1.0-w);
}

// amount is not in degrees, it is in fractions of a circle, e.g. 10/360.
// 
void
molecule_class_info_t::rotate_rgb_in_place(float *rgb, const float &amount) const {

   float *hsv;
   hsv = new float[3];
   convert_rgb_to_hsv_in_place(rgb, hsv);
   hsv[0] += amount;
   if (hsv[0] > 1.0) hsv[0] -= 1.0;
   convert_hsv_to_rgb_in_place(hsv, rgb);
   delete [] hsv;
}

// This allocated memory for xmap_is_diff_map, xmap_is_filled and
// contour_level, but *does not filll them!*
// So they need to be filled after calling this function.
void 
molecule_class_info_t::initialize_map_things_on_read_molecule(std::string molecule_name,
							      int is_diff_map,
							      short int swap_difference_map_colours) {

   // unset coordinates, this is not a set of coordinates:
   atom_sel.n_selected_atoms = 0;
   atom_sel.mol = 0;  // tested (in set_undo_molecule()) to see if this
                      // was a coordinates molecule.  So maps have to
                      // set this to NULL.

   // Map initialization:
   n_draw_vectors = 0;
   n_diff_map_draw_vectors = 0;
   xmap_is_filled = new int[10];

   // give it some memory:
   xmap_is_diff_map = new int[10];
   xmap_is_diff_map[0] = is_diff_map;  // and set the first (only) one. 

   contour_level = new float[10];

   draw_vectors = NULL;
   diff_map_draw_vectors = NULL;

   show_unit_cell_flag = 0;
   have_unit_cell      = 0; // hmmm - CHECKME.

   map_colour = new double*[10];

   map_colour[0] = new double[4];
   if (is_diff_map == 1) { 
      if (! swap_difference_map_colours) { 
	 map_colour[0][0] = 0.2; 
	 map_colour[0][1] = 0.6; 
	 map_colour[0][2] = 0.2;
      } else { 
	 map_colour[0][0] = 0.6; 
	 map_colour[0][1] = 0.2; 
	 map_colour[0][2] = 0.2; 
      } 
   } else {
      std::vector<float> orig_colours(3);
      orig_colours[0] =  0.2;
      orig_colours[1] =  0.5;
      orig_colours[2] =  0.7;
      float rotation_size = float(imol_no) * graphics_info_t::rotate_colour_map_for_map/360.0;
      // std::cout << "rotating map colour by " << rotation_size * 360.0 << std::endl;
      std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
      map_colour[0][0] = rgb_new[0];
      map_colour[0][1] = rgb_new[1];
      map_colour[0][2] = rgb_new[2];
   } 
      
   // negative contour level
   // 
   map_colour[1] = new double[4];
   if (! swap_difference_map_colours) { 
      map_colour[1][0] = 0.6; 
      map_colour[1][1] = 0.2; 
      map_colour[1][2] = 0.2; 
   } else { 
      map_colour[1][0] = 0.2; 
      map_colour[1][1] = 0.6; 
      map_colour[1][2] = 0.2;
   } 
   name_ = molecule_name;

   drawit_for_map = 1; // display the map initially, by default

   update_map_in_display_control_widget();
   
}

// Create a new combo box for this newly created map.
void
molecule_class_info_t::update_map_in_display_control_widget() const { 

   graphics_info_t g; 

   std::string dmn = name_for_display_manager();
   if (g.display_control_window())
      display_control_map_combo_box(g.display_control_window(), 
				    dmn.c_str(),
				    imol_no_ptr);

} 

void
molecule_class_info_t::update_mol_in_display_control_widget() const { 

   graphics_info_t g;

   // we don't want to add a display control hbox if we are simply
   // doing an undo: This is now deal with by the calling function.
   // 
//    std::cout << "update_mol_in_display_control_widget() now" << std::endl;
//    std::cout << "update_mol_in_display_control_widget() passed derefrerence imol_no_ptr: "
// 	     << *imol_no_ptr << std::endl;
   std::string dmn = name_for_display_manager();
   if (g.display_control_window()) 
      update_name_in_display_control_molecule_combo_box(g.display_control_window(), 
							dmn.c_str(), 
							imol_no_ptr);
}

void
molecule_class_info_t::new_mol_in_display_control_widget() const { 

   graphics_info_t g;

   // we don't want to add a display control hbox if we are simply
   // doing an undo: This is now deal with by the calling function.
   // 
   std::string dmn = name_for_display_manager();
   if (g.display_control_window()) 
      display_control_molecule_combo_box(g.display_control_window(), 
							dmn.c_str(), 
							imol_no_ptr);
}

std::string
molecule_class_info_t::name_for_display_manager() const { 

   std::string s("");
   if (graphics_info_t::show_paths_in_display_manager_flag) { 
      s = name_;
   } else {
      if (has_model()) {
	 std::string::size_type islash = name_.find_last_of("/");
	 if (islash == std::string::npos) { 
	    s = name_;
	 } else {
	    s = name_.substr(islash+1, name_.length());
	 }
      } else {
	 // This is a map, so we want to strip of xxx/ from each of
	 // the (space separated) strings.
	 // e.g.:
	 // thing/other.mtz -> other.mtz
	 // but
	 // Averged -> Averaged

	 std::vector<std::string> v = coot::util::split_string(name_, " ");
	 for (unsigned int i=0; i<v.size(); i++) {
	    if (i > 0) 
	       s += " ";
	    std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(v[i]);
	    if (p.second == "")
	       s += v[i];
	    else 
	       s += p.second;
	 }
      } 
   }
   return s;
}

std::string
molecule_class_info_t::dotted_chopped_name() const {

   std::string ss = coot::util::int_to_string(imol_no);
   ss += " " ;
   int ilen = name_.length();
   int left_size = ilen-graphics_info_t::go_to_atom_menu_label_n_chars_max;
   if (left_size <= 0) {
      // no chop
      left_size = 0;
   } else {
      // chop
      ss += "...";
   } 
   ss += name_.substr(left_size, ilen);
   return ss;
}



void
molecule_class_info_t::add_to_labelled_atom_list(int atom_index) {

   // note initialization n_labelled_atoms is 0;
   // 
   if (n_labelled_atoms < coot::MAX_LABELLED_ATOMS) {
      if ( is_in_labelled_list(atom_index) == 1 ) {
	 unlabel_atom(atom_index);
      } else { 
	 if (! is_in_labelled_list(atom_index)) {
	    // cout << "adding atom index " << atom_index << " to  "
	    // << labelled_atom_index_list << endl;
	    
	    labelled_atom_index_list[n_labelled_atoms] = atom_index;
	    n_labelled_atoms++;
	 }
      }
   } else {
      cout << "Too many labelled atoms" << endl;
   } 
}

int
molecule_class_info_t::labelled_atom(int i) {

   return labelled_atom_index_list[i];

} 


int
molecule_class_info_t::max_labelled_atom() {

   return n_labelled_atoms;
} 

// or as we would say in lisp: rember
void
molecule_class_info_t::unlabel_atom(int i) {

   // This is a classic recursive function.  I really should learn how
   // to do recussion in c++, because currently this is terrible and
   // ugly and makes my fingers itch to type it.
   // 
   // Remove i from the list of atoms to be labelled.
   //
   int offset = 0;
   //   cout << "trying to unlabel atom " << i << endl; 
   
   for (int ii=0; ii<n_labelled_atoms; ii++) {
      if (labelled_atom_index_list[ii] == i) {
	 // cout << "found " << i << " at label_atom index " << ii << endl; 
	 offset++; 
      }
      if (offset > 0) {
	 labelled_atom_index_list[ii] =
	    labelled_atom_index_list[ii+offset];
      }
   }
   if (offset > 0)
      n_labelled_atoms--;

}

void
molecule_class_info_t::unlabel_last_atom() {
   // remove the last atom from the list (if
   // there *are* atoms in the list, else do
   // nothing).
   if (n_labelled_atoms > 0) {
      unlabel_atom(labelled_atom_index_list[n_labelled_atoms-1]);
   }
}

// or as we would say in lisp: member?
bool
molecule_class_info_t::is_in_labelled_list(int i) {

   // is the i'th atom in the list of atoms to be labelled?

   for (int ii=0; ii<n_labelled_atoms; ii++) {

      if (labelled_atom_index_list[ii] == i) {
	 return 1;
      }
   }

   return 0;
}

// ------------------------- residue exists? -------------------------

int
molecule_class_info_t::does_residue_exist_p(const std::string &chain_id,
					    int resno,
					    const std::string &inscode) const {
   int state = 0;
   if (atom_sel.n_selected_atoms > 0) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
	 
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 CChain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "ERROR:: bad nchains in molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in ... " << std::endl;
	       } else {
		  PCResidue residue_p;
		  if (chain_id == chain_p->GetChainID()) { 
		     int nres = chain_p->GetNumberOfResidues();
		     for (int ires=0; ires<nres; ires++) { 
			residue_p = chain_p->GetResidue(ires);
			if (resno == residue_p->seqNum) {
			   if (inscode == residue_p->GetInsCode()) {
			      state = 1;
			      break;
			   }
			}
		     }
		  }
	       }
	    }
	 }
	 if (state)
	    break;
      }
   }
   return state;
}



// ------------------------- symmmetry atom labels -------------------------

void
molecule_class_info_t::add_atom_to_labelled_symm_atom_list(int atom_index,
							   const symm_trans_t &symm_trans,
							   const Cell_Translation &pre_shift_cell_trans) {

   // note initialization n_labelled_atoms is 0;
   // 
   if (n_labelled_symm_atoms < coot::MAX_LABELLED_ATOMS) {
      if ( is_in_labelled_symm_list(atom_index) == 1 ) {
	 unlabel_symm_atom(atom_index);
      } else { 
	 labelled_symm_atom_index_list[n_labelled_symm_atoms] = atom_index;

// 	 cout << "add_atom_to_labelled_symm_atom_list(..) adding "
// 	      << symm_trans << " at position " << n_labelled_symm_atoms
// 	      << endl;
	    
	 labelled_symm_atom_symm_trans_[n_labelled_symm_atoms] =
	    std::pair<symm_trans_t, Cell_Translation> (symm_trans, pre_shift_cell_trans);
	 n_labelled_symm_atoms++;

      }
   } else {
      cout << "Too many symmetry labelled atoms" << endl;
   } 
}


int
molecule_class_info_t::labelled_symm_atom(int i) {

   return labelled_symm_atom_index_list[i];

}

std::pair<symm_trans_t, Cell_Translation> 
molecule_class_info_t::labelled_symm_atom_symm_trans(int i) {

   return labelled_symm_atom_symm_trans_[i];
}



int
molecule_class_info_t::max_labelled_symm_atom() {

   return n_labelled_symm_atoms;
}

void
molecule_class_info_t::unlabel_symm_atom(int i) {
   
   int offset = 0;   
   
   for (int ii=0; ii<n_labelled_symm_atoms - 1; ii++) {
      if (labelled_symm_atom_index_list[ii] == i) {
	 offset++; 
      }
      if (offset > 0) {
	 labelled_symm_atom_index_list[ii] =
	    labelled_symm_atom_index_list[ii+offset];
      }
   }
   n_labelled_symm_atoms--;

}

// shall we pass the symm_trans too?  Ideally we should, I think.
//
bool
molecule_class_info_t::is_in_labelled_symm_list(int i) {

   // is the i'th atom in the list of atoms to be labelled?

   for (int ii=0; ii<n_labelled_symm_atoms; ii++) {

      if (labelled_symm_atom_index_list[ii] == i) {
	 return 1;
      }
   }

   return 0;
}


int molecule_class_info_t::add_atom_label(char *chain_id, int iresno, char *atom_id) { 

   // int i = atom_index(chain_id, iresno, atom_id);
   int i = atom_spec_to_atom_index(std::string(chain_id),
				   iresno,
				   std::string(atom_id));
   if (i > 0) 
      add_to_labelled_atom_list(i);
   else
      std::cout << atom_id << "/" << iresno << "/" << chain_id
		<< " is not found in this molecule: (" <<  imol_no << ") "
		<< name_ << std::endl; 

   return i; 
}


int molecule_class_info_t::remove_atom_label(char *chain_id, int iresno, char *atom_id) {

   int i = atom_index(chain_id, iresno, atom_id);
   if (i > 0) 
      unlabel_atom(i);
   return i; 
}

void
molecule_class_info_t::compile_density_map_display_list() {

   // std::cout << "Deleting theMapContours " << theMapContours << std::endl;
   glDeleteLists(theMapContours, 1);
   theMapContours = glGenLists(1);
   glNewList(theMapContours, GL_COMPILE);

   draw_density_map(0); // don't use theMapContours (make them!)
   glEndList();
}


// old and apparently faster (both at rotating (surprisingly) and
// generation of the lines) immediate mode
// 
// This is for electron density of course, not a surface as molecular
// modellers would think of it.
// 
void 
molecule_class_info_t::draw_density_map(short int display_lists_for_maps_flag) {

   // std::cout << "   draw_density_map() called for " << imol_no << std::endl;
   int nvecs = n_draw_vectors;
   if (drawit_for_map) {

      // same test as has_map():
      if (xmap_is_filled[0]) {

	 if (display_lists_for_maps_flag) {

	    // std::cout << " debug  Call list" << std::endl;
	    glCallList(theMapContours);
	    
	 } else { 

	    // std::cout << "   debug draw immediate mode " << std::endl;
	    if ( nvecs > 0 ) {
	       // std::cout << " debug some vectors " << nvecs << std::endl;

	       coot::Cartesian start, finish;
	       int linesdrawn = 0;
	       //

	       glColor3dv (map_colour[0]);
	       glLineWidth(graphics_info_t::map_line_width);
      
	       glBegin(GL_LINES);
	       for (int i=0; i< nvecs; i++) { 
		  glVertex3f(draw_vectors[i].getStart().x(),
			     draw_vectors[i].getStart().y(),
			     draw_vectors[i].getStart().z());
		  glVertex3f(draw_vectors[i].getFinish().x(),
			     draw_vectors[i].getFinish().y(),
			     draw_vectors[i].getFinish().z());
		  if((++linesdrawn & 1023) == 0){
		     linesdrawn = 0;
		     glEnd();
		     glBegin(GL_LINES);
		  }
	       }
	       glEnd();
	    }

	    if (xmap_is_diff_map[0] == 1) {

	       if (n_diff_map_draw_vectors > 0) { 
	       
		  glColor3dv (map_colour[1]);
		  // we only need to do this if it wasn't done above.
		  if (n_draw_vectors == 0)
		     glLineWidth(graphics_info_t::map_line_width);
	       
		  int linesdrawn_dm = 0;
		  glBegin(GL_LINES);
		  for (int i=0; i< n_diff_map_draw_vectors; i++) { 
		  
		     glVertex3f(diff_map_draw_vectors[i].getStart().get_x(),
				diff_map_draw_vectors[i].getStart().get_y(),
				diff_map_draw_vectors[i].getStart().get_z());
		     glVertex3f(diff_map_draw_vectors[i].getFinish().get_x(),
				diff_map_draw_vectors[i].getFinish().get_y(),
				diff_map_draw_vectors[i].getFinish().get_z());
		     if((++linesdrawn_dm & 1023) == 0) {
			linesdrawn_dm = 0;
			glEnd();
			glBegin(GL_LINES);
		     }
		  }
		  glEnd();
	       }
	    }
	 }
      }
   }
}

 
void
molecule_class_info_t::draw_molecule(short int do_zero_occ_spots) {


   //
   //
   //  We also need a c++ object to store molecular information
   //

   if (has_model()) { 
      if (drawit == 1) {
	 if (!cootsurface) { 
	    if (do_zero_occ_spots)
	       zero_occupancy_spots();
	    display_bonds();

	    // ghosts
// 	    std::cout << "debug ghosts: " << show_ghosts_flag << " " << ncs_ghosts.size()
// 		      << std::endl;
	    if (show_ghosts_flag) {
	       if (ncs_ghosts.size() > 0) {
		  for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
		     display_ghost_bonds(ighost);
		  }
	       }
	    }
	 }
      }
   }
}


void
molecule_class_info_t::zero_occupancy_spots() const {

   if (bonds_box.n_zero_occ_spot > 0) { 

      glColor3f(0.8, 0.7, 0.7);
      glPointSize(6.5);
      glBegin(GL_POINTS); 
      for (int i=0; i<bonds_box.n_zero_occ_spot; i++) { 
	 glVertex3f(bonds_box.zero_occ_spot[i].x(),
		    bonds_box.zero_occ_spot[i].y(),
		    bonds_box.zero_occ_spot[i].z());
      }
      glEnd();
   }
}

void
molecule_class_info_t::display_ghost_bonds(int ighost) {


   // debug
//     int ilines = 0;
//     for (int i=0; i<ncs_ghosts[ighost].bonds_box.num_colours; i++)
//        ilines += ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines;
//     std::cout << " ghost " << ighost << " has "
// 	      << ncs_ghosts[ighost].bonds_box.num_colours
// 	      << " colours and " << ilines << " lines\n";

//    std::cout << "debug ighost: " << ighost << " ncs_ghosts.size(): " << ncs_ghosts.size()
// 	     << std::endl;
   if (ighost<int(ncs_ghosts.size())) {
//       std::cout << "debug ncs_ghosts[" << ighost << "].display_it_flag "
// 		<< ncs_ghosts[ighost].display_it_flag << std::endl;
      if (ncs_ghosts[ighost].display_it_flag) {
	 glLineWidth(ghost_bond_width);
	 int c;
	 for (int i=0; i<ncs_ghosts[ighost].bonds_box.num_colours; i++) {
	    c = atom_colour(atom_sel.atom_selection[i]->element);
	    if (ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines > 0) 
	       set_bond_colour_by_mol_no(ighost);
	    glBegin(GL_LINES);
	    for (int j=0; j< ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines; j++) {
	       glVertex3f(ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getStart().get_x(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getStart().get_y(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getStart().get_z());
	       glVertex3f(ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getFinish().get_x(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getFinish().get_y(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getFinish().get_z());
	    }
	    glEnd();
	 }
      }
   }
}


// This used to use an int_grid.  Goodbye.  Clipper skeletonization is better.
// 
// void
// molecule_class_info_t::update_skeleton() {

// }


void
molecule_class_info_t::display_bonds() {

   //

   coot::CartesianPair pair;
   Lines_list ll;
   
   glLineWidth(bond_width);
   //

   for (int i=0; i<bonds_box.num_colours; i++) {

      ll = bonds_box.bonds_[i];
      //cout << "j range: for i = " << i << " is "
      //	   << bonds_box.bonds_[i].num_lines << endl;

      if ( bonds_box.bonds_[i].num_lines > 256000) {
	 std::cout << "Fencepost heuristic failure bonds_box.bonds_[i].num_lines "
	      << bonds_box.bonds_[i].num_lines << std::endl;
      }
      // std::cout << "debug:: bonds_box_type " << bonds_box_type << std::endl;
      if (bonds_box_type != coot::COLOUR_BY_RAINBOW_BONDS) {
	 // if test suggested by Ezra Peisach.
	 if (bonds_box.bonds_[i].num_lines > 0)
	    set_bond_colour_by_mol_no(i); // outside inner loop
      } else {
	 set_bond_colour_by_colour_wheel_position(i, coot::COLOUR_BY_RAINBOW);
      }
      int linesdrawn = 0;
      glBegin(GL_LINES); 
      for (int j=0; j< bonds_box.bonds_[i].num_lines; j++) {

// 	 if ( j > 200000) {
// 	    cout << "Heuristics fencepost failure j " << j << endl;
// 	    exit(1);
// 	 }
	 
	 glVertex3f(ll.pair_list[j].getStart().get_x(),
		    ll.pair_list[j].getStart().get_y(),
		    ll.pair_list[j].getStart().get_z());
	 glVertex3f(ll.pair_list[j].getFinish().get_x(),
		    ll.pair_list[j].getFinish().get_y(),
		    ll.pair_list[j].getFinish().get_z());
	 if ( (++linesdrawn & 1023) == 0) {
	    glEnd();
	    glBegin(GL_LINES);
	    linesdrawn = 0;
	 }
      }
      glEnd();
   }

   if ((show_symmetry == 1) && (graphics_info_t::show_symmetry == 1)) {
      int isymop;

      for (unsigned int isym=0; isym<symmetry_bonds_box.size(); isym++) { 
	 // isymop = isym;
	 isymop = symmetry_bonds_box[isym].second.first.isym();

	 if (symmetry_bonds_box[isym].first.symmetry_has_been_created == 1) {
	 
	    for (int icol=0; icol<symmetry_bonds_box[isym].first.num_colours; icol++) {
	 
	       set_symm_bond_colour_mol_and_symop(icol, isymop);
	       int linesdrawn = 0;
	    
	       ll = symmetry_bonds_box[isym].first.symmetry_bonds_[icol];
	 
	       glBegin(GL_LINES); 
	       for (int j=0; j< symmetry_bonds_box[isym].first.symmetry_bonds_[icol].num_lines; j++) {

		  // pair = ll.pair_list[j];
	    
		  glVertex3f(ll.pair_list[j].getStart().get_x(),
			     ll.pair_list[j].getStart().get_y(),
			     ll.pair_list[j].getStart().get_z());
		  glVertex3f(ll.pair_list[j].getFinish().get_x(),
			     ll.pair_list[j].getFinish().get_y(),
			     ll.pair_list[j].getFinish().get_z());
		  if ( (++linesdrawn & 1023) == 0) {
		     glEnd();
		     glBegin(GL_LINES);
		     linesdrawn = 0;
		  }
	       }
	       glEnd();
	    }
	 }
      }

      if (show_strict_ncs_flag == 1) {
	 // isn -> i_strict_ncs
	 for (unsigned int isn=0; isn<strict_ncs_bonds_box.size(); isn++) {

// 	    std::cout << "here 1 "
// 		      << isn << " " 
// 		      << strict_ncs_bonds_box[isn].first.symmetry_has_been_created
// 		      << " \n" ;

	    if (strict_ncs_bonds_box[isn].first.symmetry_has_been_created == 1) {
	 
	       for (int icol=0; icol<strict_ncs_bonds_box[isn].first.num_colours; icol++) {
	 
// 		  std::cout << "here 3 - isn: " << isn << " "
// 			    << "icol: " << icol << " num lines: "
// 			    << strict_ncs_bonds_box[isn].first.symmetry_bonds_[icol].num_lines
// 			    << "\n" ;
	       
		  set_symm_bond_colour_mol_and_symop(icol, isn);
		  int linesdrawn = 0;
	    
		  ll = strict_ncs_bonds_box[isn].first.symmetry_bonds_[icol];
	 
		  glBegin(GL_LINES); 
		  for (int j=0; j< strict_ncs_bonds_box[isn].first.symmetry_bonds_[icol].num_lines; j++) {

		     // pair = ll.pair_list[j];
	    
		     glVertex3f(ll.pair_list[j].getStart().get_x(),
				ll.pair_list[j].getStart().get_y(),
				ll.pair_list[j].getStart().get_z());
		     glVertex3f(ll.pair_list[j].getFinish().get_x(),
				ll.pair_list[j].getFinish().get_y(),
				ll.pair_list[j].getFinish().get_z());
		     if ( (++linesdrawn & 1023) == 0) {
			glEnd();
			glBegin(GL_LINES);
			linesdrawn = 0;
		     }
		  }
		  glEnd();
	       }
	    }
	 }
      }
   }
}

// 
void
molecule_class_info_t::update_map_triangles(float radius, coot::Cartesian centre) {
   
   CIsoSurface<float> my_isosurface;
   coot::CartesianPairInfo v;
   int isample_step = 1;
   graphics_info_t g;

   // std::cout   << "DEBUG:: g.zoom: " << g.zoom << std::endl;

   if (g.dynamic_map_resampling == 1)
      // isample_step = 1 + int (0.009*g.zoom);
      isample_step = 1 + int (0.009*(g.zoom + g.dynamic_map_zoom_offset));

   if (isample_step > 15) 
      isample_step = 15;

   // for critical points of size display and resampling being different:
   // 
   float dy_radius = radius;
   if (g.dynamic_map_size_display == 1)
      if (isample_step <= 15 )
	 dy_radius *= float(isample_step);
      else
	 dy_radius *= 15.0;

   // 
   if (isample_step <= 0) { 
      std::cout << "WARNING:: Bad zoom   ("<< g.zoom 
		<< "):  setting isample_step to 1" << std::endl;
      isample_step = 1;
   }
   if (dy_radius <= 0.0) { 
      std::cout << "WARNING:: Bad radius (" << dy_radius 
		<< ") setting to 10" << std::endl;
      dy_radius = 10.0;
   }

   // dynamically transformed maps get their vectors from molecule B
   // (we are looking at molecule at atoms in molecule A) which then
   // have a the inverse of that transformation applied to them.
   //
   // But note that to get from the centre in the A molecule to the
   // corresponding centre in the B molecule we need to apply the
   // *inverse* of the transformation in the map_ghost_info.

   if (is_dynamically_transformed_map_flag) {
      clipper::Coord_orth c(centre.x(), centre.y(), centre.z());
      clipper::Coord_orth ct = c.transform(map_ghost_info.rtop.inverse());
      centre = coot::Cartesian(ct.x(), ct.y(), ct.z());
   }

   for (int i=0; i< max_xmaps; i++) {
      if (xmap_is_filled[i] == 1) {
	 v = my_isosurface.GenerateSurface_from_Xmap(xmap_list[0],
						     contour_level[i],
						     dy_radius, centre,
						     isample_step);
	 if (is_dynamically_transformed_map_flag)
	    dynamically_transform(v);
	 set_draw_vecs(v.data, v.size);
      }
      if (xmap_is_diff_map[i] == 1) {
	 v = my_isosurface.GenerateSurface_from_Xmap(xmap_list[0],
						     -contour_level[i],
						     dy_radius, centre,
						     isample_step);
	 if (is_dynamically_transformed_map_flag)
	    dynamically_transform(v);
	 set_diff_map_draw_vecs(v.data, v.size);
      }
   }
}

// modify v
void
molecule_class_info_t::dynamically_transform(coot::CartesianPairInfo v) {

   int s = v.size;
   for (int i=0; i<s; i++) {
      clipper::Coord_orth c1(v.data[i].getStart().x(),
			     v.data[i].getStart().y(),
			     v.data[i].getStart().z());
      clipper::Coord_orth c2(v.data[i].getFinish().x(),
			     v.data[i].getFinish().y(),
			     v.data[i].getFinish().z());
      clipper::Coord_orth ct1 = c1.transform(map_ghost_info.rtop);
      clipper::Coord_orth ct2 = c2.transform(map_ghost_info.rtop);
      v.data[i] = coot::CartesianPair(coot::Cartesian(ct1.x(), ct1.y(), ct1.z()),
				coot::Cartesian(ct2.x(), ct2.y(), ct2.z()));
   }
   
}


// 
void
molecule_class_info_t::map_fill_from_mtz(std::string mtz_file_name,
					 std::string f_col,
					 std::string phi_col,
					 std::string weight_col,
					 int use_weights,
					 int is_diff_map) {

   short int use_reso_flag = 0;
   short int is_anomalous_flag = 0;
   map_fill_from_mtz_with_reso_limits(mtz_file_name,
				      f_col,
				      phi_col,
				      weight_col,
				      use_weights,
				      is_anomalous_flag,
				      is_diff_map,
				      use_reso_flag, 0.0, 0.0); // don't use these reso limits.

}


// 
void
molecule_class_info_t::map_fill_from_mtz_with_reso_limits(std::string mtz_file_name,
							  std::string f_col,
							  std::string phi_col,
							  std::string weight_col,
							  int use_weights,
							  short int is_anomalous_flag,
							  int is_diff_map,
							  short int use_reso_limits,
							  float low_reso_limit,
							  float high_reso_limit) {

   graphics_info_t g;

   // save for potential phase recombination in refmac later
   if (use_weights) { 
      fourier_f_label = f_col; 
      fourier_phi_label = phi_col;
      fourier_weight_label = weight_col; // magic label, we can go
					 // combining if this is not
					 // "";
//       std::cout << "DEBUG:: saving fourier_weight_label: " <<
// 	 fourier_weight_label << std::endl;
   }


   // std::cout << "DEBUG:: reso tinkering " << use_reso_limits << std::endl;
   clipper::Resolution user_resolution(high_reso_limit);
   clipper::Resolution fft_reso; // filled later

   clipper::HKL_info myhkl; 
   clipper::MTZdataset myset; 
   clipper::MTZcrystal myxtl; 

   std::cout << "reading mtz file..." << std::endl; 

   long T0 = 0; // timer
   T0 = glutGet(GLUT_ELAPSED_TIME);

   clipper::CCP4MTZfile mtzin; 
   mtzin.open_read( mtz_file_name );       // open new file 
   mtzin.import_hkl_info( myhkl );         // read sg, cell, reso, hkls
   clipper::HKL_data< clipper::datatypes::F_sigF<float> >   f_sigf_data(myhkl, myxtl);
   clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data(myhkl, myxtl);
   clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphidata(myhkl, myxtl); 
   
   std::string mol_name = mtz_file_name + " "; 
   mol_name += f_col; 
   mol_name += " ";
   mol_name += phi_col; 
   
   if (use_weights) { 
      mol_name += " ";
      mol_name += weight_col; 
   }

   if (use_reso_limits) {
      mol_name += " ";
      mol_name += g.float_to_string(low_reso_limit);
      mol_name += " ";
      mol_name += g.float_to_string(high_reso_limit);
   }
   
   initialize_map_things_on_read_molecule(mol_name,
					  is_diff_map,
					  g.swap_difference_map_colours);
   xmap_list = new clipper::Xmap<float>[5];  // 18Kb each (empty) Xmap
   max_xmaps++;

   // If use weights, use both strings, else just use the first
   std::pair<std::string, std::string> p = make_import_datanames(f_col, phi_col, weight_col, use_weights);

   if (p.first.length() == 0) { // mechanism to signal an error
      std::cout << "ERROR:: fill_map.. - There was a column label error.\n";
   } else { 
      
      if (use_weights) {
	 // 	 std::cout << "DEBUG:: Importing f_sigf_data: " << p.first << std::endl;
	 mtzin.import_hkl_data( f_sigf_data, myset, myxtl, p.first );
	 // std::cout << "DEBUG:: Importing phi_fom_data: " << p.second << std::endl;
	 mtzin.import_hkl_data(phi_fom_data, myset, myxtl, p.second);
	 mtzin.close_read();
	 fphidata.compute(f_sigf_data, phi_fom_data,
			  clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());
      } else {
	 // std::cout << "DEBUG:: Importing f_phi_data: " << p.first << std::endl;
	 mtzin.import_hkl_data(fphidata, myset, myxtl, p.first);
	 mtzin.close_read();
      }
   
      long T1 = glutGet(GLUT_ELAPSED_TIME);

      int n_reflections = myhkl.num_reflections();
      std::cout << "Number of reflections: " << n_reflections << "\n";
      if (n_reflections <= 0) {
	 std::cout << "WARNING:: No reflections in mtz file!?" << std::endl;
      } else { 
	 if (use_reso_limits) {
	    fft_reso = user_resolution;
	    filter_by_resolution(&fphidata, low_reso_limit, high_reso_limit);
	 } else {
	    // fft_reso = myhkl.resolution();
	    // Kevin says do this instead:
	    fft_reso = clipper::Resolution(1.0/sqrt(fphidata.invresolsq_range().max()));
	 }
      
	 if (is_anomalous_flag) {
	    fix_anomalous_phases(&fphidata);
	 } 
   
   
	 cout << "finding ASU unique map points with sampling rate "
   	      << graphics_info_t::map_sampling_rate	<< endl;
         clipper::Grid_sampling gs(myhkl.spacegroup(),
				    myhkl.cell(),
				    fft_reso,
				    graphics_info_t::map_sampling_rate);
	 cout << "done grid sampling..." << gs.format() << endl; 
	 xmap_list[0].init( myhkl.spacegroup(), myhkl.cell(), gs); // 1.5 default
	 cout << "Grid..." << xmap_list[0].grid_sampling().format() << "\n";
   
	 long T2 = glutGet(GLUT_ELAPSED_TIME);
// 	 std::cout << "MTZ:: debug:: " << myhkl.spacegroup().symbol_hm() << " " 
// 		   << myhkl.cell().descr().a() << " " 
// 		   << myhkl.cell().descr().b() << " " 
// 		   << myhkl.cell().descr().c() << " " 
// 		   << clipper::Util::rad2d(myhkl.cell().descr().alpha()) << " " 
// 		   << clipper::Util::rad2d(myhkl.cell().descr().beta ()) << " " 
// 		   << clipper::Util::rad2d(myhkl.cell().descr().gamma()) << std::endl;
// 	 std::cout << "MTZ:: debug:: n_reflections: " << myhkl.num_reflections()
// 		   << std::endl;
// 	 int ncount = 0;
// 	 clipper::HKL_info::HKL_reference_index hri;
// 	 for (hri=fphidata.first(); !hri.last(); hri.next()) {
// 	    if (ncount < 500) 
// 	       std::cout << " MTZ fphi: " << hri.hkl().h() << " "
// 			 << hri.hkl().k() << " " << hri.hkl().l() << " "
// 			 << fphidata[hri].f() << " "
// 			 << clipper::Util::rad2d(fphidata[hri].phi()) << std::endl;
// 	    ncount++;
// 	 } 
   
	 cout << "doing fft..." << endl;
	 xmap_list[0].fft_from( fphidata );                  // generate map
	 cout << "done fft..." << endl;
   
	 long T3 = glutGet(GLUT_ELAPSED_TIME);
	 std::cout << "INFO:: " << float(T1-T0)/1000.0 << " seconds to read MTZ file\n";
	 std::cout << "INFO:: " << float(T2-T1)/1000.0 << " seconds to initialize map\n";
	 std::cout << "INFO:: " << float(T3-T2)/1000.0 << " seconds for FFT\n";
	 xmap_is_filled[0] = 1;  // set the map-is-filled? flag
   
  
	 // Fill the class variables:
	 //   clipper::Map_stats stats(xmap_list[0]);
	 //   map_mean_ = stats.mean();
	 //   map_sigma_ = stats.std_dev();

	 mean_and_variance<float> mv = map_density_distribution(xmap_list[0], 0);

	 save_mtz_file_name = mtz_file_name;
	 save_f_col = f_col;
	 save_phi_col = phi_col;
	 save_weight_col = weight_col;
	 save_use_weights = use_weights;
	 save_is_anomalous_map_flag = is_anomalous_flag;
	 save_is_diff_map_flag = is_diff_map;
	 save_high_reso_limit = high_reso_limit;
	 save_low_reso_limit = low_reso_limit;
	 save_use_reso_limits = use_reso_limits;

	 // 
	 map_mean_  = mv.mean; 
	 map_sigma_ = sqrt(mv.variance);
	 map_max_   = mv.max_density;
	 map_min_   = mv.min_density;

	 long T4 = glutGet(GLUT_ELAPSED_TIME);
	 std::cout << "INFO:: " << float(T4-T3)/1000.0 << " seconds for statistics\n";

	 std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
	 std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
	 std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
	 std::cout << "      Map minimum: ..... " << map_min_ << std::endl;

	 set_initial_contour_level();

	 // update_map_colour_menu_manual(g.n_molecules, name_.c_str()); 
	 // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str()); 

	 imol_no = graphics_info_t::n_molecules; // gets updated outside afterwards
	 update_map();
	 long T5 = glutGet(GLUT_ELAPSED_TIME);
	 std::cout << "INFO:: " << float(T5-T4)/1000.0 << " seconds for contour map\n";
	 std::cout << "INFO:: " << float(T5-T0)/1000.0 << " seconds in total\n";

	 // save state strings

	 //
  
	 if (have_sensible_refmac_params) { 
	    save_state_command_strings_.push_back("make-and-draw-map-with-refmac-params");
	    save_state_command_strings_.push_back(single_quote(coot::util::intelligent_debackslash(mtz_file_name)));
	    save_state_command_strings_.push_back(single_quote(f_col));
	    save_state_command_strings_.push_back(single_quote(phi_col));
	    save_state_command_strings_.push_back(single_quote(weight_col));
	    save_state_command_strings_.push_back(g.int_to_string(use_weights));
	    save_state_command_strings_.push_back(g.int_to_string(is_diff_map));
	    save_state_command_strings_.push_back(g.int_to_string(1)); // have refmac params
	    save_state_command_strings_.push_back(single_quote(refmac_fobs_col));
	    save_state_command_strings_.push_back(single_quote(refmac_sigfobs_col));
	    save_state_command_strings_.push_back(single_quote(refmac_r_free_col));
	    save_state_command_strings_.push_back(g.int_to_string(refmac_r_free_flag_sensible));
	 } else { 
	    save_state_command_strings_.push_back("make-and-draw-map");
	    save_state_command_strings_.push_back(single_quote(coot::util::intelligent_debackslash(mtz_file_name)));
	    save_state_command_strings_.push_back(single_quote(f_col));
	    save_state_command_strings_.push_back(single_quote(phi_col));
	    save_state_command_strings_.push_back(single_quote(weight_col));
	    save_state_command_strings_.push_back(g.int_to_string(use_weights));
	    save_state_command_strings_.push_back(g.int_to_string(is_diff_map));
	 }
      }
   }
}

// Return a pair.first string of length 0 on error to construct dataname(s).
std::pair<std::string, std::string>
molecule_class_info_t::make_import_datanames(const std::string &f_col_in,
					     const std::string &phi_col_in,
					     const std::string &weight_col_in,
					     int use_weights) const {

   // If use_weights return 2 strings, else set something useful only for pair.first

   std::string f_col = f_col_in;
   std::string phi_col = phi_col_in;
   std::string weight_col = weight_col_in;

   std::string::size_type islash_f   =      f_col.find_last_of("/");
   std::string::size_type islash_phi =    phi_col.find_last_of("/");
   short int label_error = 0; 

   if (islash_f != std::string::npos) {
      // f_col is of form e.g. xxx/yyy/FWT
      if (f_col.length() > islash_f)
	 f_col = f_col.substr(islash_f+1);
      else
	 label_error = 1;
   }

   if (islash_phi != std::string::npos) {
      // phi_col is of form e.g. xxx/yyy/PHWT
      if (phi_col.length() > islash_phi)
	 phi_col = phi_col.substr(islash_phi+1);
      else
	 label_error = 1;
   }

   if (use_weights) { 
      std::string::size_type islash_fom = weight_col.find_last_of("/");
      if (islash_fom != std::string::npos) {
	 // weight_col is of form e.g. xxx/yyy/WT
	 if (weight_col.length() > islash_fom)
	    weight_col = weight_col.substr(islash_fom+1);
	 else
	    label_error = 1;
      }
   }

  
   std::pair<std::string, std::string> p("", "");

   if (!label_error) {
      std::string no_xtal_dataset_prefix= "/*/*/";
      if (use_weights) { 
	 p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " +      f_col + "]";
	 p.second = no_xtal_dataset_prefix + "[" + phi_col + " " + weight_col + "]";
      } else {
	 p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " + phi_col + "]";
      }
   }
   return p;
}



void
molecule_class_info_t::filter_by_resolution(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata,
					    const float &reso_low,
					    const float &reso_high) const {

   float inv_low  = 1.0/(reso_low*reso_low);
   float inv_high = 1.0/(reso_high*reso_high);
   int n_data = 0;
   int n_reset = 0;

      
   for (clipper::HKL_info::HKL_reference_index hri = fphidata->first(); !hri.last(); hri.next()) {
//        std::cout << "high: " << inv_high << " low: " << inv_low
//  		<< " data: " << hri.invresolsq() << std::endl;
      n_data++;

      if ( hri.invresolsq() > inv_low &&
	   hri.invresolsq() < inv_high) {
      } else {
	 (*fphidata)[hri].f() = 0.0;
	 n_reset++;
      } 
   }
   std::cout << "Chopped " << n_reset << " data out of " << n_data << std::endl;
}

void
molecule_class_info_t::set_initial_contour_level() {

   float level = 1.0;
   if (xmap_is_diff_map[0]) {
      if (map_sigma_ > 0.05) {
	 level = nearest_step(map_mean_ +
			      graphics_info_t::default_sigma_level_for_fofc_map*map_sigma_, 0.01);
      } else {
	 level = 3.0*map_sigma_;
      }
   } else { 
      if (map_sigma_ > 0.05) {
	 level = nearest_step(map_mean_ + graphics_info_t::default_sigma_level_for_map*map_sigma_, 0.01);
      } else {
	 level = graphics_info_t::default_sigma_level_for_map * map_sigma_;
      }
   }
   contour_level[0] = level;

}

void
molecule_class_info_t::fix_anomalous_phases(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata) const {

   for (clipper::HKL_info::HKL_reference_index hri = fphidata->first(); !hri.last(); hri.next()) {
      (*fphidata)[hri].shift_phase(-M_PI_2);
   }
} 



void
molecule_class_info_t::test_label_symm_atom(int i) {
   //
 
   // same test as has_model():
   if (has_model()) {
      
      int max = max_labelled_symm_atom();
   
      if (max > 0) {       

	 if (i < atom_sel.n_selected_atoms) { 

	    int iatom_index = labelled_symm_atom(i);
	    std::pair <symm_trans_t, Cell_Translation> st = labelled_symm_atom_symm_trans(i);
	 
	    std::string label =
	       make_symm_atom_label_string(atom_sel.atom_selection[iatom_index], st.first);

	    GLfloat blueish[3] = { 0.7, 0.7, 1.0 };
	 
	    glColor3fv(blueish);
	 
	    coot::Cartesian symm_point = translate_atom_with_pre_shift(atom_sel, iatom_index, st);

	    glRasterPos3f(symm_point.get_x(),
			  symm_point.get_y()+0.02,
			  symm_point.get_z()+0.02);
	 
	    printString(label);
	 }
      }
   }
}

void
molecule_class_info_t::label_symm_atom(int i, symm_trans_t symm_trans) {
   //

   // same test as has_model():
   if (atom_sel.n_selected_atoms > 0 ) { 

      if (atom_sel.n_selected_atoms > 0) {

	 std::string label =
	    make_symm_atom_label_string(atom_sel.atom_selection[i], symm_trans);

	 GLfloat blueish[3] = { 0.7, 0.8, 1.0 };

	 glColor3fv(blueish);

	 std::cout << "label_symm_atom :" << symm_trans << std::endl;
	 coot::Cartesian symm_point =
	    translate_atom(atom_sel, i, symm_trans);
            
	 glRasterPos3f(symm_point.get_x(),
		       symm_point.get_y()+0.02,
		       symm_point.get_z()+0.02);

	 std::cout << "adding symm label: " << label << std::endl;
	 printString(label);
	 std::cout << "done  symm label: in  label_symm_atom(..): " << label << std::endl;

      }
   }   
}

// Put a label at the ith atom of mol_class_info::atom_selection. 
//
void
molecule_class_info_t::label_atom(int i, int brief_atom_labels_flag) {

   // same test as has_model():
   if (atom_sel.n_selected_atoms > 0 ) { 

      if (i < atom_sel.n_selected_atoms) { 

	 PCAtom atom = (atom_sel.atom_selection)[i];

	 if (atom) { 

	    std::string label = make_atom_label_string(atom, brief_atom_labels_flag);
      
	    // GLfloat white[3] = { 1.0, 1.0, 1.0 };
	    GLfloat pink[3] =  { graphics_info_t::font_colour.red,
				 graphics_info_t::font_colour.green,
				 graphics_info_t::font_colour.blue };
      
	    // glClear(GL_COLOR_BUFFER_BIT);
	    glColor3fv(pink);
	    // glShadeModel (GL_FLAT);
      
	    glRasterPos3f((atom)->x, (atom)->y+0.02, (atom)->z +0.02);
	    printString(label);
      
	 }
      } else { 
	 std::cout << "INFO:: trying to label atom out of range: " 
		   << i << " " << atom_sel.n_selected_atoms 
		   << " Removing label\n";
	 unlabel_atom(i);
      }
   }
}


void
molecule_class_info_t::set_have_unit_cell_flag_maybe() {
   
   CMMDBCryst *cryst_p = atom_sel.mol->get_cell_p();

   mat44 my_matt;

   int err = cryst_p->GetTMatrix(my_matt, 0, 0, 0, 0);

   if (err != 0) {
      have_unit_cell = 0;
      std::cout << "No Symmetry for this model" << std::endl;
   } else { 
      have_unit_cell = 1;
   }
}

void
molecule_class_info_t::save_previous_map_colour() {

   if (has_map()) { 
      previous_map_colour.resize(3);
      for (int i=0; i<3; i++) 
	 previous_map_colour[i] = map_colour[0][i];
   }
}


void
molecule_class_info_t::restore_previous_map_colour() {

   if (has_map()) { 
      if (previous_map_colour.size() == 3) { 
	 for (int i=0; i<3; i++) 
	    map_colour[0][i] = previous_map_colour[i];
      }
   }
   compile_density_map_display_list();
}



// Where is this used?
int
molecule_class_info_t::next_free_map() {

   return 0; // placeholder FIXME.
}


// This function is not used in anger, right?
//
// (a debugging function).
// 
void
check_static_vecs_extents() {

   //
   int imol = 0;
   
   graphics_info_t g;

   cout << "checking extents of the " << g.molecules[imol].n_draw_vectors
	<< " in the graphics static" << endl;

   coot::Cartesian first, second;
   float max_x = -9999, min_x = 9999;
   float max_y = -9999, min_y = 9999;
   float max_z = -9999, min_z = 9999;
   
   for (int i=0; i<g.molecules[imol].n_draw_vectors; i++) {
      first  = g.molecules[imol].draw_vectors[i].getStart();
      second = g.molecules[imol].draw_vectors[i].getFinish();

      if (first.get_x() < min_x) min_x = first.get_x();
      if (first.get_y() < min_y) min_y = first.get_y();
      if (first.get_z() < min_z) min_z = first.get_z();

      if (second.get_x() < min_x) min_x = second.get_x();
      if (second.get_y() < min_y) min_y = second.get_y();
      if (second.get_z() < min_z) min_z = second.get_z();

      if (first.get_x() > max_x) max_x = first.get_x();
      if (first.get_y() > max_y) max_y = first.get_y();
      if (first.get_z() > max_z) max_z = first.get_z();

      if (second.get_x() > max_x) max_x = second.get_x();
      if (second.get_y() > max_y) max_y = second.get_y();
      if (second.get_z() > max_z) max_z = second.get_z();
      
   }
   cout << min_x << " " << max_x << endl
	<< min_y << " " << max_y << endl
	<< min_z << " " << max_z << endl;
}


void
molecule_class_info_t::makebonds(float min_dist, float max_dist) {

   Bond_lines_container bonds(atom_sel, min_dist, max_dist);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::NORMAL_BONDS;
}

void
molecule_class_info_t::makebonds(float max_dist) {
   
   Bond_lines_container bonds(atom_sel, max_dist);

   // bonds.check(); 

   bonds_box = bonds.make_graphical_bonds();

   // cout << "makebonds bonds_box.num_colours "
   // << bonds_box.num_colours << endl;
}

void
molecule_class_info_t::makebonds() {

   int do_disulphide_flag = 1;
   Bond_lines_container bonds(atom_sel, do_disulphide_flag, draw_hydrogens_flag);
   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::NORMAL_BONDS;

}

void
molecule_class_info_t::make_ca_bonds(float min_dist, float max_dist) {

   Bond_lines_container bonds;
   bonds.do_Ca_bonds(atom_sel, min_dist, max_dist);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::CA_BONDS;
   // std::cout << "ca: bonds_box_type is now " << bonds_box_type << std::endl;

}

void
molecule_class_info_t::make_ca_bonds() { 
   make_ca_bonds(2.4, 4.7); 
}

void
molecule_class_info_t::make_ca_plus_ligands_bonds() { 

   Bond_lines_container bonds;
   bonds.do_Ca_plus_ligands_bonds(atom_sel, 2.4, 4.7);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::CA_BONDS;
   bonds_box_type = coot::COLOUR_BY_RAINBOW_BONDS; // FIXME
   // std::cout << "ca: bonds_box_type is now " << bonds_box_type << std::endl;
}

void
molecule_class_info_t::make_colour_by_chain_bonds(short int change_c_only_flag) {
   // 
   Bond_lines_container bonds;
   bonds.do_colour_by_chain_bonds(atom_sel, draw_hydrogens_flag, change_c_only_flag);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;

   if (graphics_info_t::glarea) 
      graphics_info_t::graphics_draw();
} 

void
molecule_class_info_t::make_colour_by_molecule_bonds() { 

   // 
   Bond_lines_container bonds;
   bonds.do_colour_by_molecule_bonds(atom_sel, draw_hydrogens_flag);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_MOLECULE_BONDS;

   if (graphics_info_t::glarea) 
      graphics_info_t::graphics_draw();

} 



void 
molecule_class_info_t::make_bonds_type_checked() { 

   if (bonds_box_type == coot::NORMAL_BONDS)
      makebonds();
   if (bonds_box_type == coot::CA_BONDS)
      make_ca_bonds();
   if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS)
      // Baah, we have to use the static in graphics_info_t here as it
      // is not a per-molecule property.
      make_colour_by_chain_bonds(graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag);
   if (bonds_box_type == coot::COLOUR_BY_MOLECULE_BONDS)
      make_colour_by_molecule_bonds();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS)
      make_ca_plus_ligands_bonds();
   if (bonds_box_type == coot::BONDS_NO_WATERS)
      bonds_no_waters_representation();
   if (bonds_box_type == coot::BONDS_SEC_STRUCT_COLOUR)
      bonds_sec_struct_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR)
      ca_plus_ligands_sec_struct_representation();
   if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS)
      ca_plus_ligands_rainbow_representation();
   if (bonds_box_type == coot::COLOUR_BY_OCCUPANCY_BONDS)
      occupancy_representation();
   if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS)
      b_factor_representation();
   
   update_ghosts();
}


// 
void
molecule_class_info_t::draw_skeleton() {


   if (has_map()) { 

      coot::CartesianPair pair;

      set_bond_colour(grey);
      glLineWidth(2.0);

      if (greer_skeleton_draw_on == 1) {
      
	 //cout << "greer_skeleton_draw_on: "
	 //	   << greer_skel_box.bonds_[0].num_lines<< endl;

	 glBegin(GL_LINES);
	 for (int j=0; j<greer_skel_box.bonds_[0].num_lines; j++) {

            glVertex3f(greer_skel_box.bonds_[0].pair_list[j].getStart().get_x(),
		       greer_skel_box.bonds_[0].pair_list[j].getStart().get_y(),
		       greer_skel_box.bonds_[0].pair_list[j].getStart().get_z());
            glVertex3f(greer_skel_box.bonds_[0].pair_list[j].getFinish().get_x(),
		       greer_skel_box.bonds_[0].pair_list[j].getFinish().get_y(),
		       greer_skel_box.bonds_[0].pair_list[j].getFinish().get_z());
	 }
	 glEnd();
      }

      if (fc_skeleton_draw_on == 1) { 

	 for (int l=0; l<fc_skel_box.num_colours; l++) {
 	    if (colour_skeleton_by_random) {
	       //  	       set_skeleton_bond_colour_random(l, colour_table);
	       set_skeleton_bond_colour(0.96);
 	    } else {
// 	       std::cout << "skel: " << l
// 			 << " of  " <<  fc_skel_box.num_colours <<  " " 
// 			 << (float(l)/float(fc_skel_box.num_colours)+0.01)/1.011
// 			 << std::endl;
	       set_skeleton_bond_colour( (float(l)/float(fc_skel_box.num_colours)+0.01)/1.011 );
	    }


	    glBegin(GL_LINES);
	    for (int j=0; j<fc_skel_box.bonds_[l].num_lines; j++) {

	       glVertex3f(fc_skel_box.bonds_[l].pair_list[j].getStart().get_x(),
			  fc_skel_box.bonds_[l].pair_list[j].getStart().get_y(),
			  fc_skel_box.bonds_[l].pair_list[j].getStart().get_z());
	       glVertex3f(fc_skel_box.bonds_[l].pair_list[j].getFinish().get_x(),
			  fc_skel_box.bonds_[l].pair_list[j].getFinish().get_y(),
			  fc_skel_box.bonds_[l].pair_list[j].getFinish().get_z());
	    }
	    glEnd();
	 }

	 // debugging:
	 //       cout << "Molecule name: " << name_ << endl; 
	 //       for (int l=0; l<fc_skel_box.num_colours; l++) {
	 // 	 cout << "skeleton levels in draw_skeleton: level: " 
	 // 	      << l << " " << skel_levels[l] << endl; 
	 //       }
      }
   }
}

// Added rotate colour_map for EJD 5/5/2004.
void
molecule_class_info_t::set_skeleton_bond_colour(float f) {

   float rotation_size = float(imol_no) * 2.0*graphics_info_t::rotate_colour_map_on_read_pdb/360.0;
   while (rotation_size > 1.0) {
      rotation_size -= 1.0;
   }

   std::vector<float> c(3);
   c[0] = 0.1+0.6*f*graphics_info_t::skeleton_colour[0];
   c[1] = 0.1+0.9*f*graphics_info_t::skeleton_colour[1];
   c[2] = 0.1+0.2*f*graphics_info_t::skeleton_colour[2];
   std::vector<float> rgb_new = rotate_rgb(c, rotation_size);

   glColor3f(rgb_new[0], rgb_new[1], rgb_new[2]);
}



void
molecule_class_info_t::set_colour_skeleton_by_segment() { // use random colouring
   
   colour_skeleton_by_random = 1;
} 

void
molecule_class_info_t::set_colour_skeleton_by_level() { // use random colouring
   
   colour_skeleton_by_random = 0;
} 


//
void
molecule_class_info_t::draw_fc_skeleton() {
   
}

//
void
molecule_class_info_t::update_clipper_skeleton() {

   // Create map extents (the extents of the skeletonization)
   // from the current centre.

   if (xskel_is_filled == 1) { 

      graphics_info_t g;

      if (xmap_is_filled[0] == 1 && xmap_is_diff_map[0] != 1) { 
	 //
	 float skeleton_box_radius = g.skeleton_box_radius; 
	 graphics_info_t g;

	 GraphicalSkel cowtan; 

	 // fc_skel_box: class object type graphical_bonds_container
	 //
	 cout<< "making graphical bonds..." << endl; 
	 fc_skel_box = cowtan.make_graphical_bonds(xmap_list[0],xskel_cowtan,
						   g.RotationCentre(),
						   skeleton_box_radius,
						   g.skeleton_level);

// 	 cout << "DEBUG: " << "fc_skel_box.num_lines = "
// 	      << fc_skel_box.num_colours << endl;
// 	 cout << "DEBUG: " << "fc_skel_box.bonds_ = "
// 	      << fc_skel_box.bonds_ << endl;
// 	 cout << "DEBUG: " << "fc_skel_box.bonds_[0].num_lines = "
// 	      << fc_skel_box.bonds_[0].num_lines << endl;
      }
   }
}

void
molecule_class_info_t::unskeletonize_map() { 

   //    std::cout << "DEBUG:: unskeletonize_map" << std::endl;
   fc_skeleton_draw_on = 0;
   xskel_is_filled = 0;
   clipper::Xmap<int> empty; 
   xskel_cowtan = empty;
   // std::cout << "DEBUG:: done unskeletonize_map" << std::endl;
} 

// //
// void
// molecule_class_info_t::update_fc_skeleton_old() {

//    // Create map extents (the extents of the skeletonization)
//    // from the current centre.

//    if (xmap_is_filled[0] == 1 && xmap_is_diff_map[0] != 1) { 
//       //
//       float skeleton_box_radius = 20.0;
//       graphics_info_t g;

//       fc_skeleton foadi_cowtan; 

//       // fc_skel_box: class object type graphical_bonds_container
//       //
//       fc_skel_box = foadi_cowtan.itskel(xmap_list[0], 0.50);

//       cout << "DEBUG: " << "fc_skel_box.num_lines = "
// 	   << fc_skel_box.num_colours << endl;
//       cout << "DEBUG: " << "fc_skel_box.bonds_ = "
// 	   << fc_skel_box.bonds_ << endl;
//       cout << "DEBUG: " << "fc_skel_box.bonds_[0].num_lines = "
// 	   << fc_skel_box.bonds_[0].num_lines << endl;
//    }
// }

// Return -1 on error
int
molecule_class_info_t::read_ccp4_map(std::string filename, int is_diff_map_flag,
				     const std::vector<std::string> &acceptable_extensions) {

   // For now, where we try to read in a map and we crash in importing
   // a file that is not a ccp4 map, lets do some checking: first that
   // the file exists and is not a directory, then that the file has
   // an extension of ".map" or ".ext".  If not, then complain and
   // return having done nothing.

   // stat filename 
   struct stat s;
   int status = stat(filename.c_str(), &s);
   if (status != 0 ||  !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // an error
   } else { 
      
      // was a regular file, let's check the extension:
      // 
      std::string::size_type islash = filename.find_last_of("/");
      std::string tstring;
      if (islash == std::string::npos) { 
	 // no slash found
	 tstring = filename;
      } else { 
	 tstring = filename.substr(islash + 1);
      }

      bool good_extension_flag = 0;
      for (unsigned int iextension=0; iextension<acceptable_extensions.size(); iextension++) {
	 std::string::size_type imap = tstring.rfind(acceptable_extensions[iextension]);
	 if (imap != std::string::npos) {
	    good_extension_flag = 1;
	    break;
	 }
      }
      
      // not really extension checking, just that it has it in the
      // filename:
      if (good_extension_flag == 0) { 
	 
	 std::cout << "Filename for a CCP4 map must end in .map or .ext "
		   << "or some other approved extension - sorry\n";
	 return -1;
	 std::string ws = "The filename for a CCP4 map must\n";
	 ws += "currently end in .map or .ext - sorry.\n\n";
	 ws += "The map must be a CCP4 map or Badness Will Happen! :-)\n";
	 GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(ws);
	 gtk_widget_show(w);
      }
   }

   // KDC: check map type
   enum MAP_FILE_TYPE { CCP4, CNS };
   MAP_FILE_TYPE map_file_type;
   {
     FILE* file = fopen( filename.c_str(), "r" );
     int c1, c2;
     c1 = c2 = 0;
     for ( int i = 0; i < 16; i++ ) {
       int c = getc(file);
       if ( c == EOF ) break;
       if ( c == 0 )                                 c1++;
       if ( isalpha(c) || isdigit(c) || isspace(c) ) c2++;
     }
     if ( c1 > c2 ) map_file_type = CCP4;
     else           map_file_type = CNS;
   }

   if (max_xmaps == 0) {
      std::cout << "allocating space in read_ccp4_map" << std::endl;
      xmap_list = new clipper::Xmap<float>[10];
      xmap_is_filled   = new int[10];
      xmap_is_diff_map = new int[10];
      contour_level    = new float[10];
   }

   short int bad_read = 0;
   max_xmaps++; 

   if ( map_file_type == CCP4 ) {
     std::cout << "attempting to read CCP4 map: " << filename << std::endl;
     clipper::CCP4MAPfile file;
     file.open_read(filename);
     try {
       file.import_xmap( xmap_list[0] );
     }
     catch (clipper::Message_base exc) {
       std::cout << "WARNING:: failed to read " << filename
		 << " Bad ASU (inconsistant gridding?)." << std::endl;
       bad_read = 1;
     }
     file.close_read();
   } else {
     std::cout << "attempting to read CNS map: " << filename << std::endl;
     clipper::CNSMAPfile file;
     file.open_read(filename);
     try {
       file.import_xmap( xmap_list[0] );
     }
     catch (clipper::Message_base exc) {
       std::cout << "WARNING:: failed to read " << filename << std::endl;
       bad_read = 1;
     }
     file.close_read();
   }

   if (bad_read == 0) { 
      initialize_map_things_on_read_molecule(filename, is_diff_map_flag, 
					     graphics_info_t::swap_difference_map_colours);

      mean_and_variance<float> mv = map_density_distribution(xmap_list[0], 1); 

      float mean = mv.mean; 
      float var = mv.variance;
      contour_level[0]    = nearest_step(mean + 1.5*sqrt(var), 0.05);
      xmap_is_filled[0] = 1;  // set the map-is-filled? flag
      xmap_is_diff_map[0] = is_diff_map_flag; // but it may be...
      // fill class variables
      map_mean_ = mv.mean;
      map_sigma_ = sqrt(mv.variance);
      map_max_   = mv.max_density;
      map_min_   = mv.min_density;

      set_initial_contour_level();

      std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
      std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
      std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
      std::cout << "      Map minimum: ..... " << map_min_ << std::endl;

      // save state strings
      save_state_command_strings_.push_back("handle-read-ccp4-map");
      save_state_command_strings_.push_back(single_quote(filename));
      save_state_command_strings_.push_back(graphics_info_t::int_to_string(is_diff_map_flag));
   
      update_map();
   }

   int stat = 0;
   if (bad_read)
      stat = -1;
   return stat;
}

void
molecule_class_info_t::new_map(const clipper::Xmap<float> &map_in, std::string name_in) {

   if (max_xmaps == 0) {
      xmap_list        = new clipper::Xmap<float>[1];
      xmap_is_filled   = new int[1];
      xmap_is_diff_map = new int[1];
      contour_level    = new float[1];
   }
   xmap_list[0] = map_in; 
   max_xmaps++;
   xmap_is_filled[0] = 1;
   // the map name is filled by using set_name(std::string)
   // sets name_ to name_in:
   initialize_map_things_on_read_molecule(name_in, 0, 0); // not a diff_map

   mean_and_variance<float> mv = map_density_distribution(xmap_list[0], 1); 

   float mean = mv.mean; 
   float var = mv.variance; 
   contour_level[0]  = nearest_step(mean + 1.5*sqrt(var), 0.05);
   xmap_is_filled[0] = 1;  // set the map-is-filled? flag

   // fill class variables
   map_mean_ = mv.mean;
   map_sigma_ = sqrt(mv.variance);

   update_map(); 
   
}

void
molecule_class_info_t::set_name(std::string name) {
   name_ = name;
}


int
molecule_class_info_t::make_map_from_phs(std::string pdb_filename,
                                         std::string phs_filename) {

   int iret = -1; // default error return status
   //
   std::cout << "Make a map from " << phs_filename << " using "
	     << pdb_filename << " for the cell and symmetry " << std::endl; 

   atom_selection_container_t SelAtom = get_atom_selection(pdb_filename);

   if (SelAtom.read_success == 1) { // success
      mat44 my_matt;

      int err = SelAtom.mol->GetTMatrix(my_matt, 0,0,0,0);
      if (err != 0) {
         cout << "!! Cryst::GetTMatrix() fails in make_map_from_phs"
              << endl;
      } else { // success
         cout << "Cryst::GetTMatrix() is good in make_map_from_phs"
              << endl;

         // we get to the cell and symmetry thusly:
         
         cout << "PHS Cell:  "
              << SelAtom.mol->get_cell().a << " "
              << SelAtom.mol->get_cell().b << " "
              << SelAtom.mol->get_cell().c << " "
              << SelAtom.mol->get_cell().alpha << " "
              << SelAtom.mol->get_cell().beta << " "
              << SelAtom.mol->get_cell().gamma << " "
              << endl;

         cout << "PHS Spacegroup: " << SelAtom.mol->GetSpaceGroup() << endl; 

	 clipper::Spacegroup spacegroup(clipper::Spgr_descr(SelAtom.mol->GetSpaceGroup()));
	 clipper::Cell cell(clipper::Cell_descr(SelAtom.mol->get_cell().a,
						SelAtom.mol->get_cell().b,
						SelAtom.mol->get_cell().c,
						clipper::Util::d2rad(SelAtom.mol->get_cell().alpha),
						clipper::Util::d2rad(SelAtom.mol->get_cell().beta),
						clipper::Util::d2rad(SelAtom.mol->get_cell().gamma))); 

	 iret = make_map_from_phs(spacegroup, cell, phs_filename);
      }
   }
   return iret;
}


int
molecule_class_info_t::make_map_from_phs_using_reso(std::string phs_filename,
						    const clipper::Spacegroup &sg,
						    const clipper::Cell &cell,
						    float reso_limit_low,
						    float reso_limit_high) {

   clipper::PHSfile phs;

   phs.open_read(phs_filename);

   // std::cout << "creating resolution" << std::endl;
   clipper::Resolution resolution(reso_limit_high);

   clipper::HKL_info mydata(sg, cell, resolution);
   clipper::HKL_data<clipper::datatypes::F_sigF<float>  >  myfsig(mydata);
   clipper::HKL_data<clipper::datatypes::Phi_fom<float> >  myphwt(mydata);
   clipper::HKL_data<clipper::datatypes::F_phi<float>   >  fphidata(mydata); 

   std::cout << "importing info" << std::endl;
   phs.import_hkl_info(mydata);
   std::cout << "importing data" << std::endl;
   phs.import_hkl_data(myfsig); 
   phs.import_hkl_data(myphwt);

   phs.close_read(); 

   std::cout << "PHS file: Number of reflections: " << mydata.num_reflections() << "\n";

   fphidata.update();

   fphidata.compute(myfsig, myphwt, 
 		    clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

//    for (int i=0; i<10; i++) {
//       std::cout << "checking phi weight: " << i << " " << myphwt[i].phi() << "  " << myphwt[i].fom()
// 		<< std::endl;
//        std::cout << "checking f    sigf: " << i << " " << myfsig[i].f() << "   "
// 		 << myfsig[i].sigf() << std::endl;
//        std::cout << "checking missing: " << i << " " << myfsig[i].missing() << " "
// 		 << myphwt[i].missing() << " " << fphidata[i].missing() << std::endl;
//        // << " " << fphidata[i].phi() <<
//    }
   
  if (max_xmaps == 0) {
     
     xmap_list = new clipper::Xmap<float>[10];
     xmap_is_filled   = new int[10];
     xmap_is_diff_map = new int[10];
     contour_level    = new float[10];
   }

   max_xmaps++; 

  std::string mol_name = phs_filename; 

  initialize_map_things_on_read_molecule(mol_name, 0, 0); // not diff map

  cout << "initializing map..."; 
  xmap_list[0].init(mydata.spacegroup(), 
		    mydata.cell(), 
		    clipper::Grid_sampling(mydata.spacegroup(),
					   mydata.cell(), 
					   mydata.resolution()));
  cout << "done."<< endl; 

//   cout << "Map Grid (from phs file)..." 
//        << xmap_list[0].grid_sampling().format()
//        << endl;  

  cout << "doing fft..." ; 
  xmap_list[0].fft_from( fphidata );                  // generate map
  cout << "done." << endl;

  mean_and_variance<float> mv = map_density_distribution(xmap_list[0],0);

  cout << "Mean and sigma of map from PHS file: " << mv.mean 
       << " and " << sqrt(mv.variance) << endl;

  // fill class variables
  map_mean_ = mv.mean;
  map_sigma_ = sqrt(mv.variance);

  xmap_is_diff_map[0] = 0; 
  xmap_is_filled[0] = 1; 
  contour_level[0] = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);

  graphics_info_t g; 
  // update_map_colour_menu_manual(g.n_molecules, name_.c_str()); 
  // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str());

  g.scroll_wheel_map = g.n_molecules; // change the current scrollable map.

  std::cout << "updating map..." << std::endl;
  update_map();
  std::cout << "done updating map..." << std::endl;

  // how do we restore this map?
  save_state_command_strings_.push_back("read-phs-and-make-map-using-cell-symm");
  save_state_command_strings_.push_back(single_quote(phs_filename));
  save_state_command_strings_.push_back(single_quote(sg.symbol_hm()));
  save_state_command_strings_.push_back(g.float_to_string(cell.descr().a()));
  save_state_command_strings_.push_back(g.float_to_string(cell.descr().b()));
  save_state_command_strings_.push_back(g.float_to_string(cell.descr().c()));
  save_state_command_strings_.push_back(g.float_to_string(clipper::Util::rad2d(cell.descr().alpha())));
  save_state_command_strings_.push_back(g.float_to_string(clipper::Util::rad2d(cell.descr().beta())));
  save_state_command_strings_.push_back(g.float_to_string(clipper::Util::rad2d(cell.descr().gamma())));

  return g.n_molecules;
}



// and the molecule number imol_coords where the coordinates are.
int
molecule_class_info_t::make_map_from_cif(std::string cif_file_name,
					 int imol_coords) {

   return make_map_from_cif(cif_file_name,
			    graphics_info_t::molecules[imol_coords].atom_sel); 
}

// and the molecule number imol_coords where the coordinates are.
int
molecule_class_info_t::make_map_from_cif_2fofc(std::string cif_file_name,
					       int imol_coords) {

   return make_map_from_cif_2fofc(cif_file_name,
				  graphics_info_t::molecules[imol_coords].atom_sel); 
}

// and the molecule number imol_coords where the coordinates are.
int
molecule_class_info_t::make_map_from_cif_fofc(std::string cif_file_name,
					      int imol_coords) {

   return make_map_from_cif_generic(cif_file_name,
				    graphics_info_t::molecules[imol_coords].atom_sel,
				    2);  // 2 -> is Fo-Fc map
}


int
molecule_class_info_t::make_map_from_cif(std::string cif_file_name,
					 atom_selection_container_t SelAtom) {

   // 0 is not is_2fofc_type map (is sigmaa)
   return make_map_from_cif_generic(cif_file_name, SelAtom, 0);

}

int
molecule_class_info_t::make_map_from_cif_2fofc(std::string cif_file_name,
					 atom_selection_container_t SelAtom) {

   // 1 is is_2fofc_type map (not sigmaa)
   return make_map_from_cif_generic(cif_file_name, SelAtom, 1); 

}


int
molecule_class_info_t::make_map_from_cif_generic(std::string cif_file_name,
						 atom_selection_container_t SelAtom,
						 short int is_2fofc_type) {

#ifdef HAVE_CIF

   clipper::HKL_info mydata;
   clipper::CIFfile cif; 
      
   cif.open_read ( cif_file_name ); 
   cif.import_hkl_info(mydata); // set spacegroup, cell and get hkl list. 
   clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf(mydata);
   
   cif.import_hkl_data(myfsigf);
   cif.close_read();

   std::cout << "make_map_from_cif: Read " << mydata.num_reflections()
	     << " from CIF file." << std::endl;

   if (mydata.num_reflections() == 0) {
      std::cout << "Error reading cif file, can't make a map" << std::endl;
      return -1; // Error
   }
   return calculate_sfs_and_make_map(cif_file_name, mydata, myfsigf,
				     SelAtom, is_2fofc_type);
}
   
int
molecule_class_info_t::calculate_sfs_and_make_map(const std::string &mol_name,
						  const clipper::HKL_info &mydata,
						  const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf,
						  atom_selection_container_t SelAtom,
						  short int is_2fofc_type) {

   if (max_xmaps == 0) {
      xmap_list = new clipper::Xmap<float>[10];
      xmap_is_filled   = new int[10];
      xmap_is_diff_map = new int[10];
      contour_level    = new float[10];
   }
   max_xmaps++; 
   initialize_map_things_on_read_molecule(mol_name, 0, 0); // not diff map
   
   std::cout << "calculating structure factors..." << std::endl;

   // Fix up fphidata to contain the calculated structure factors

   // Calculated structure factors go here:
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphidata(mydata);
   // map coefficients ((combined Fo and scaled Fc) and calc phi) go here:

   clipper::HKL_data< clipper::datatypes::F_phi<float> > map_fphidata(mydata);
  
   // get a list of all the atoms
   clipper::MMDBAtom_list atoms(SelAtom.atom_selection, SelAtom.n_selected_atoms);
  
   std::cout << "isotropic fft of " << SelAtom.n_selected_atoms
	     << " atoms..." << std::endl;
   clipper::SFcalc_iso_fft<float>(fphidata, atoms);
   std::cout << "done iso fft..." << std::endl;

   // debug:: examine fphidata and myfsigf:
   std::cout << "DEBUG:: myfsigf  has " <<  myfsigf.data_size() << " data"
	     << std::endl;
   std::cout << "DEBUG:: fphidata has " << fphidata.data_size() << " data"
	     << std::endl;

   if (0) { // debug
      float sum_fo = 0;
      float sum_fc = 0;
      int n_fo = 0;
      int n_fc = 0;
      for (clipper::HKL_info::HKL_reference_index ih=mydata.first();
	   !ih.last(); ih.next()) {
	 if (!myfsigf[ih].missing()) {
	    n_fo++;
	    n_fc++;
	    sum_fo += myfsigf[ih].f();
	    sum_fc = fphidata[ih].f();
	 }
      }
      
      std::cout << "DEBUG:: fo: sum average: " << sum_fo << " " << sum_fo/float(n_fo)
		<< std::endl; 
      std::cout << "DEBUG:: fc: sum average: " << sum_fc << " " << sum_fc/float(n_fc)
		<< std::endl;
      for (clipper::HKL_info::HKL_reference_index ih=mydata.first();
	   !ih.last(); ih.next())
	      std::cout << "DEBUG::  myfsigf " <<  " " <<  myfsigf[ih].f() << " "
		   << myfsigf[ih].sigf() << " " << myfsigf[ih].missing() << std::endl;
      for (int i=0; i<10; i++)
	 std::cout << "DEBUG:: fphidata " << i << " " << fphidata[i].f()
		   << " " << fphidata[i].phi() << std::endl;
   }
   
   int nprm = 10;
   std::vector<clipper::ftype> params_init( nprm, 1.0 );
   // clipper::BasisFn_spline basis_f1f2( mydata, nprm, 2.0 );
   //  target_f1f2( fc, fo );
   //clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<float>,
   //                          clipper::datatypes::F_sigF<float> >
   // target_f1f2( fphidata, myfsigf );
   //clipper::ResolutionFn fscale( mydata, basis_f1f2,
   //                              target_f1f2, params_init );

   int i=0;
   float r_top = 0.0, r_bot = 0.0;
   float sum_fo = 0.0, sum_fc = 0.0, sum_scale = 0.0;
   int n_data = 0;

   if (is_2fofc_type == molecule_map_type::TYPE_2FO_FC ||
       is_2fofc_type == molecule_map_type::TYPE_FO_FC) {

      if (is_2fofc_type == molecule_map_type::TYPE_2FO_FC)
	 std::cout << "INFO:: calculating 2fofc map..." << std::endl; 
      if (is_2fofc_type == molecule_map_type::TYPE_FO_FC)
	 std::cout << "INFO:: calculating fofc map..." << std::endl;
      
      clipper::BasisFn_spline basis_f1f2( mydata, nprm, 2.0 );
      //  target_f1f2( fc, fo );
      clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<float>,
	 clipper::datatypes::F_sigF<float> >
	 target_f1f2( fphidata, myfsigf );
      clipper::ResolutionFn fscale( mydata, basis_f1f2, target_f1f2, params_init );

      float multiplier = 2.0;
      if (is_2fofc_type == molecule_map_type::TYPE_FO_FC)
	 multiplier = 1.0;
     
      for ( clipper::HKL_info::HKL_reference_index ih=mydata.first();
	    !ih.last(); ih.next() ) { 
	 map_fphidata[ih].phi() = fphidata[ih].phi(); 
	 if (!myfsigf[ih].missing()) {
	    map_fphidata[ih].f() = multiplier*myfsigf[ih].f() -
	       fphidata[ih].f()*sqrt(fscale.f(ih));
	    float top_tmp = fabs(myfsigf[ih].f() - fphidata[ih].f()*sqrt(fscale.f(ih)));
	    if (0) { 
	       std::cout << "debug:: fobs: " << myfsigf[ih].f() << " fscale: "
			 << fphidata[ih].f() << " scale: " << fscale.f(ih)
			 << std::endl;
	    }
	    
	    r_top += top_tmp;
	    r_bot += fabs(myfsigf[ih].f());
// 	    std::cout << "debug:: adding to top: " << top_tmp << " bot: "
// 		      << fabs(myfsigf[ih].f()) << std::endl;
	    sum_fo += myfsigf[ih].f();
	    sum_fc += fphidata[ih].f();
	    sum_scale += sqrt(fscale.f(ih));
	    n_data++; 
	 } else {
	    map_fphidata[ih].f() = 0.0; 
	 }
      }

   } else { // not 2fofc-style, i.e. is sigmaa style
      
      if (is_2fofc_type == molecule_map_type::TYPE_SIGMAA) {
     
	 std::cout << "sigmaa and scaling..." << std::endl;

	 // need an mmdb
	 CMMDBManager *mmdb = SelAtom.mol;
	 // clipper::MTZcrystal cxtl; // and a cxtl, whatever that is...
	 clipper::Cell cxtl = myfsigf.hkl_info().cell();
	 // hkls
	 // fo
	 // 
	    
	 // get a list of all the atoms
	 clipper::mmdb::PPCAtom psel;
	 int hndl, nsel;
	 hndl = mmdb->NewSelection();
	 mmdb->SelectAtoms( hndl, 0, 0, SKEY_NEW );
	 mmdb->GetSelIndex( hndl, psel, nsel );
	 clipper::MMDBAtom_list atoms( psel, nsel );
	 mmdb->DeleteSelection( hndl );

	 // calculate structure factors
	 const clipper::HKL_info& hkls = mydata;
	 const clipper::HKL_data<clipper::datatypes::F_sigF<float> >& fo = myfsigf;
	 clipper::HKL_data<clipper::datatypes::F_phi<float> > fc( hkls, cxtl );
	 clipper::SFcalc_obs_bulk<float> sfcb;
	 sfcb( fc, fo, atoms );

	 // do anisotropic scaling
	 clipper::SFscale_aniso<float> sfscl;
	 // sfscl( fo, fc );  // scale Fobs
	 sfscl( fc, fo );  // scale Fcal

	 // now do sigmaa calc
	 clipper::HKL_data<clipper::datatypes::F_phi<float> >   fb( hkls, cxtl ), fd( hkls, cxtl );
	 clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phiw( hkls, cxtl );
	 clipper::HKL_data<clipper::datatypes::Flag>    flag( hkls, cxtl );
	 typedef clipper::HKL_data_base::HKL_reference_index HRI;
	 // If no free flag is available, then use all reflections..
	 for (HRI ih = flag.first(); !ih.last(); ih.next() )
	    flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	 /* This code uses free reflections only for sigmaa and scaling...
	 for (HRI ih = flag.first(); !ih.last(); ih.next() )
	    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==0) )
	       flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	    else
	       flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	 */

	 // do sigmaa calc
	 int n_refln = mydata.num_reflections();
	 int n_param = 10;
	 clipper::SFweight_spline<float> sfw( n_refln, n_param );
	 sfw( fb, fd, phiw, fo, fc, flag );

	 // OK, so now fb and fd contain F_phis, one for "best"
	 // sigmaa, one for difference map.  Let's just use the "best"
	 // map for now.
	 map_fphidata = fb;
      }
   } // is 2fofc else sigmaa style check

   std::cout << "DEBUG:: rdiffsum/rsum: " << r_top << "/" << r_bot << std::endl;
   if (r_bot>0.0) { 
      std::cout << "Isotropic R-factor: " << 100.0*r_top/r_bot << "%"
		<< " for " << n_data  << " reflections" <<  std::endl;
      std::cout << "DEBUG:: sums: fo: " << sum_fo/float(n_data) << " fc: "
		<< sum_fc/n_data << " scale: " << sum_scale/n_data << " with "
		<< n_data << " data" << std::endl;
   } else {
      std::cout << "Problem with structure factors, no structure factor sum!?"
		<< std::endl;
   } 
   cout << "initializing map..."; 
   xmap_list[0].init(mydata.spacegroup(), 
		     mydata.cell(), 
		     clipper::Grid_sampling(mydata.spacegroup(),
					    mydata.cell(), 
					    mydata.resolution()));
   cout << "done."<< endl; 


   cout << "doing fft..." ; 
   xmap_list[0].fft_from( map_fphidata ); // generate map
   cout << "done." << endl;

   mean_and_variance<float> mv = map_density_distribution(xmap_list[0],0);

   cout << "Mean and sigma of map " << mol_name << " " << mv.mean 
	<< " and " << sqrt(mv.variance) << endl; 

   // fill class variables
   map_mean_ = mv.mean;
   map_sigma_ = sqrt(mv.variance);
   map_max_   = mv.max_density;
   map_min_   = mv.min_density;
  
   xmap_is_diff_map[0] = 0; 
   xmap_is_filled[0] = 1; 

   std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
   std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
   std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
   std::cout << "      Map minimum: ..... " << map_min_ << std::endl;
   
  set_initial_contour_level();

   graphics_info_t g; 
   // update_map_colour_menu_manual(g.n_molecules, name_.c_str()); 
   // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str()); 

   g.scroll_wheel_map = g.n_molecules; // change the current scrollable map.

   int imol = g.n_molecules;
   update_map(); 
   return imol;

#else // HAVE_CIF
   return -1; 
#endif   
}

// This needs to be rationalized with the version that *does* pass the
// coordinates.
//
// This is the version that gets called when we use a file selector to
// get the file (i.e. it doesn't specify the coordinates (molecule))
// because there are calculated structure factors in the file.
//
// We make a Fc alpha-c map.  Which is not usually what we want.
// 
int
molecule_class_info_t::make_map_from_cif(std::string cif_file_name) {
   return make_map_from_cif_sigmaa(cif_file_name, molecule_map_type::TYPE_SIGMAA);
}

int
molecule_class_info_t::make_map_from_cif_diff_sigmaa(std::string cif_file_name) {
   return make_map_from_cif_sigmaa(cif_file_name,
				   molecule_map_type::TYPE_DIFF_SIGMAA);
}

// SigmaA map type, either molecule_map_type::TYPE_SIGMAA or TYPE_DIFF_SIGMAA.
// 
int
molecule_class_info_t::make_map_from_cif_sigmaa(std::string cif_file_name,
					       int sigmaa_map_type) {

#ifdef HAVE_CIF

   clipper::HKL_info mydata;
   clipper::CIFfile cif; 
      
   cif.open_read ( cif_file_name );
   cif.import_hkl_info(mydata); // set spacegroup, cell and get hkl list. 
   clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf(mydata); // Fobs
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fc(mydata); // FC PHIC

   cif.import_hkl_data(myfsigf);
   cif.import_hkl_data(fc); 

   cif.close_read(); 
      
   std::cout << "Read " << mydata.num_reflections() << " from CIF file." 
	     << std::endl; 

   if (mydata.num_reflections() == 0) {
      return -1;
   } else {

      // Are all the calculated sfs missing/zero?
      // 
      int non_zero = 0;
      for(int i=0; i< mydata.num_reflections(); i++) {
	 if (! fc[i].missing()) {
	    if (fc[i].f() > 0.0) {
	       non_zero = 1;
	       break;
	    }
	 }
      }

      if (non_zero == 0) {
	 std::cout << "WARNING:: Ooops - all the structure factor amplitudes "
		   << " appear to be zero - or missing.  " << std::endl;
	 std::cout << "WARNING:: Are you sure this file (" << cif_file_name
		   << ") contains calculated structure factors?" << std::endl;
	 std::cout << "WARNING:: No map calculated." << std::endl;
	 std::cout << "INFO:: if you want to calculate structure factors from a"
		   << " set of coordinates, " << std::endl
		   << "       consider the function (read-cif-data cif-file imol)"
		   << std::endl;
      } else {

	 if (max_xmaps == 0) {
     
	    xmap_list = new clipper::Xmap<float>[10];
	    // these are allocated in initialize_map_things_on_read_molecule()
// 	    xmap_is_filled   = new int[10];
// 	    xmap_is_diff_map = new int[10];
// 	    contour_level    = new float[10];
	 }
	 max_xmaps++; 
	 std::string mol_name = cif_file_name;
	 if (sigmaa_map_type == molecule_map_type::TYPE_SIGMAA)
	    mol_name += " SigmaA";
	 if (sigmaa_map_type == molecule_map_type::TYPE_DIFF_SIGMAA)
	    mol_name += " Difference SigmaA";

	 // new sigmaA code... needs to be updated to new Kevin
	 // code... but that is slightly tricky because here we have
	 // sfs, whereas KC code calculates them.
	 
	 std::cout << "sigmaa and scaling..." << std::endl; 
	 
	 clipper::HKL_data< clipper::datatypes::F_phi<float> > map_fphidata(mydata);
	 clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phifom(mydata);

	 clipper::Cell cxtl = myfsigf.hkl_info().cell();
	 // "Aliases" to fix Kevin's sigmaA code into mine
	 const clipper::HKL_data<clipper::datatypes::F_sigF<float> >& fo = myfsigf;
	 const clipper::HKL_info& hkls = mydata;

	 // now do sigmaa calc
	 clipper::HKL_data<clipper::datatypes::F_phi<float> >   fb( hkls, cxtl ), fd( hkls, cxtl );
	 clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phiw( hkls, cxtl );
	 clipper::HKL_data<clipper::datatypes::Flag>    flag( hkls, cxtl );
	 typedef clipper::HKL_data_base::HKL_reference_index HRI;
	 // If no free flag is available, then use all reflections..
	 for (HRI ih = flag.first(); !ih.last(); ih.next() )
	    flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;

	 /* This code uses free reflections only for sigmaa and scaling...
	 for (HRI ih = flag.first(); !ih.last(); ih.next() )
	    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==0) )
	       flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	    else
	       flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	 */

	 // do sigmaa calc
	 int n_refln = mydata.num_reflections();
	 int n_param = 10;
	 clipper::SFweight_spline<float> sfw( n_refln, n_param );
	 sfw( fb, fd, phiw, fo, fc, flag );
	 // fb is F+phi for "Best"
	 // fd is F+phi for difference map
	 short int is_diff = 0;
	 if (sigmaa_map_type == molecule_map_type::TYPE_DIFF_SIGMAA) {
	    map_fphidata = fd;
	    is_diff = 1;
	 } else { 
	    map_fphidata = fb;
	 }
	 
	 // back to old code 
	 //
	 cout << "initializing map..."; 
	 xmap_list[0].init(mydata.spacegroup(), 
			   mydata.cell(), 
			   clipper::Grid_sampling(mydata.spacegroup(),
						  mydata.cell(), 
						  mydata.resolution()));
	 cout << "done."<< endl; 

	 cout << "doing fft..." ; 
	 // xmap_list[0].fft_from( fphidata );       // generate Fc alpha-c map
	 xmap_list[0].fft_from( map_fphidata );       // generate sigmaA map 20050804
	 cout << "done." << endl;
	 initialize_map_things_on_read_molecule(mol_name, is_diff, 0);
	 // now need to fill contour_level, xmap_is_diff_map xmap_is_filled
	 if (is_diff)
	    xmap_is_diff_map[0] = 1;
	 else 
	    xmap_is_diff_map[0] = 0;

	 mean_and_variance<float> mv = map_density_distribution(xmap_list[0],0);

	 cout << "Mean and sigma of map from CIF file (make_map_from_cif): "
	      << mv.mean << " and " << sqrt(mv.variance) << endl; 

	 xmap_is_filled[0] = 1; 

	 map_mean_  = mv.mean; 
	 map_sigma_ = sqrt(mv.variance);
	 map_max_   = mv.max_density;
	 map_min_   = mv.min_density;
	 
	 set_initial_contour_level();

	 int imol = graphics_info_t::n_molecules;
	 // update_map_colour_menu_manual(graphics_info_t::n_molecules, name_.c_str()); 
	 // update_map_scroll_wheel_menu_manual(graphics_info_t::n_molecules,
	 // name_.c_str());
	 update_map(); 

	 if (sigmaa_map_type != molecule_map_type::TYPE_DIFF_SIGMAA) {
	    save_state_command_strings_.push_back("read-cif-data-with-phases-sigmaa");
	    save_state_command_strings_.push_back(single_quote(cif_file_name));
	 } else {
	    save_state_command_strings_.push_back("read-cif-data-with-phases-diff-sigmaa");
	    save_state_command_strings_.push_back(single_quote(cif_file_name));
	 } 
	 return imol;
      } 
	 
   }
   return -1; 
   
#else // HAVE_CIF

   return -1;

#endif
   
}


//
// This is the version that gets called when we use a file selector to
// get the file (i.e. it doesn't specify the coordinates (molecule)) because
// this cif file has (or it is hoped that it has) calculated structure factors.
int
molecule_class_info_t::make_map_from_cif_nfofc(std::string cif_file_name,
					       int map_type,
					       short int swap_difference_map_colours) {

   int ir = -1;
   
#ifdef HAVE_CIF

   clipper::HKL_info mydata;
   clipper::CIFfile cif; 
      
   cif.open_read ( cif_file_name );
   cif.import_hkl_info(mydata); // set spacegroup, cell and get hkl list. 
   clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf(mydata);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphidata(mydata);

   cif.import_hkl_data(myfsigf);
   cif.import_hkl_data(fphidata); 

   cif.close_read(); 
      
   std::cout << "Read " << mydata.num_reflections() << " from CIF file." 
	     << std::endl; 

   if (mydata.num_reflections() == 0) {
      return -1;
   } else {

      int non_zero = 0;
      for(int i=0; i< mydata.num_reflections(); i++) {
	 if (! fphidata[i].missing() ) {
	    if (fphidata[i].f() > 0.0) {
	       non_zero++;
	       break;
	    }
	 }
      }

      if (! non_zero) {
	 std::cout << "WARNING:: Ooops - all the calculated structure factor "
		   << "amplitudes appear"
		   << " to be zero - or missing.  " << std::endl;
	 std::cout << "WARNING:: Are you sure this file (" << cif_file_name
		   << ") contains calculated structure factors?" << std::endl;
	 std::cout << "WARNING:: No map calculated." << std::endl;
      } else {

	 if (max_xmaps == 0) {
     
	    xmap_list = new clipper::Xmap<float>[10];
	    xmap_is_filled   = new int[10];
	    xmap_is_diff_map = new int[10];
	    contour_level    = new float[10];
	 }
	 max_xmaps++; 
	 std::string mol_name = cif_file_name;

	 int is_diff_map_flag = 0;
	 if (map_type == molecule_map_type::TYPE_FO_FC) { 
	    is_diff_map_flag = 1;
	    mol_name += " Fo-Fc";
	 }
	 if (map_type == molecule_map_type::TYPE_2FO_FC) {
	    mol_name += " 2Fo-Fc";
	 }
	 if (map_type == molecule_map_type::TYPE_FO_ALPHA_CALC) {
	    mol_name += " Fo ac";
	 }
	 
	 initialize_map_things_on_read_molecule(mol_name, is_diff_map_flag,
						swap_difference_map_colours);
	
	 cout << "initializing map..."; 
	 xmap_list[0].init(mydata.spacegroup(), 
			   mydata.cell(), 
			   clipper::Grid_sampling(mydata.spacegroup(),
						  mydata.cell(), 
						  mydata.resolution()));
	 std::cout << "done."<< std::endl;

	 // Here we need to fix up fphidata to be a combination
	 // of fsigf data and fphidata.
	 //
	 float fo_multiplier = 2.0;
	 float fc_multiplier = 1.0;
	 if (map_type == molecule_map_type::TYPE_FO_FC)
	    fo_multiplier = 1.0;
	 if (map_type == molecule_map_type::TYPE_FO_ALPHA_CALC) {
	    fo_multiplier = 1.0;
	    fc_multiplier = 0.0;
	 }

	 int nprm = 10;
	 std::vector<clipper::ftype> params_init(nprm, 1.0);
	 clipper::BasisFn_spline basis_f1f2( mydata, nprm, 2.0 );
	 clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<float>,
	    clipper::datatypes::F_sigF<float> >
	    target_f1f2( fphidata, myfsigf );
	 clipper::ResolutionFn fscale( mydata, basis_f1f2, target_f1f2, params_init );

	 int nrefl = 0;
	 int nmissing = 0;
	 for (clipper::HKL_info::HKL_reference_index ih=mydata.first();
	      !ih.last(); ih.next()) {
	    nrefl++;
	    if (!myfsigf[ih].missing()) {
	       fphidata[ih].f() = fo_multiplier * myfsigf[ih].f() -
		  fc_multiplier * fphidata[ih].f() * sqrt(fscale.f(ih));
	       // std::cout << "scale: " << sqrt(fscale.f(ih)) << std::endl;
	    } else {
	       nmissing++;
	       // std::cout << "missing reflection: " << ih << std::endl;
	       fphidata[ih].f() = 0;
	    }
	 }
	 std::cout << "There were " << nrefl << " reflections of which "
		   << nmissing << " were missing\n";
	 

	 std::cout << "doing fft..." ; 
	 xmap_list[0].fft_from( fphidata );                  // generate map
	 std::cout << "done." << std::endl;

	 mean_and_variance<float> mv = map_density_distribution(xmap_list[0],0);

	 std::cout << "Mean and sigma of map from CIF file (make_map_from_cif_nfofc): "
		   << mv.mean << " and " << sqrt(mv.variance) << std::endl; 

	 if (is_diff_map_flag == 1) {
	    xmap_is_diff_map[0] = 1; 
	    contour_level[0] = nearest_step(mv.mean + 2.5*sqrt(mv.variance), 0.01);
	 } else { 
	    xmap_is_diff_map[0] = 0; 
	    contour_level[0] = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);
	 }

	 // fill class variables
	 map_mean_ = mv.mean;
	 map_sigma_ = sqrt(mv.variance);
	 xmap_is_diff_map[0] = is_diff_map_flag; 
	 xmap_is_filled[0] = 1; 

	 graphics_info_t g; 
	 // update_map_colour_menu_manual(g.n_molecules, name_.c_str()); 
	 // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str());

	 g.scroll_wheel_map = g.n_molecules; // change the current scrollable map.

	 int imol = g.n_molecules;
   
	 update_map();

	 have_unsaved_changes_flag = 0;
	 std::vector<std::string> strings;
	 if (map_type == molecule_map_type::TYPE_FO_FC)
	    strings.push_back("read-cif-data-with-phases-fo-fc");
	 else 
	    strings.push_back("read-cif-data-with-phases-2fo-fc");
	 strings.push_back(single_quote(cif_file_name));
	 save_state_command_strings_ = strings;

	 return imol;
      }
   }
#endif // HAVE_CIF
   return ir;
}

int
molecule_class_info_t::make_map_from_mtz_by_calc_phases(const std::string &mtz_file_name,
							const std::string &f_col,
							const std::string &sigf_col,
							atom_selection_container_t SelAtom,
							short int is_2fofc_type) {

   clipper::HKL_info mydata;
   clipper::CCP4MTZfile mtz;
   clipper::MTZdataset myset; 
   clipper::MTZcrystal myxtal; 

   std::cout << "reading mtz file..." << mtz_file_name << std::endl; 
   mtz.open_read(mtz_file_name);

   // make the data names for import:
   std::pair<std::string, std::string> p = make_import_datanames(f_col, sigf_col, "", 0);
   mtz.import_hkl_info(mydata); // set spacegroup, cell and get hkl list. 
   clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf(mydata, myxtal);
   mtz.import_hkl_data(myfsigf, myset, myxtal, p.first);
   mtz.close_read();
   
   return calculate_sfs_and_make_map(mtz_file_name, mydata, myfsigf,
				     SelAtom, is_2fofc_type);
}



// The rest was all interface fluff.  Here is where we do the real work
// (or get clipper to do it :).
// 
int
molecule_class_info_t::make_map_from_phs(const clipper::Spacegroup &sg,
					 const clipper::Cell &cell,
                                         std::string phs_filename) {

   // clipper::Resolution resolution(reso);  // no.

   // clipper::HKL_info mydata(sg, cell, resolution);

   clipper::PHSfile phs;

   std::cout << "PHS:: reading " << phs_filename << std::endl;
   phs.open_read(phs_filename);

   std::cout << "PHS:: creating resolution" << std::endl;
   clipper::Resolution resolution = phs.resolution(cell);
   // mydata.init(sg, cell, resolution);

   std::cout << "PHS:: creating mydata" << std::endl;
   clipper::HKL_info mydata(sg, cell, resolution);
   clipper::HKL_data<clipper::datatypes::F_sigF<float>  >  myfsig(mydata);
   clipper::HKL_data<clipper::datatypes::Phi_fom<float> >  myphwt(mydata);
   clipper::HKL_data<clipper::datatypes::F_phi<float>   >  fphidata(mydata); 

   std::cout << "PHS:: importing info" << std::endl;
   phs.import_hkl_info(mydata);
   std::cout << "PHS:: importing data" << std::endl;
   phs.import_hkl_data(myfsig); 
   phs.import_hkl_data(myphwt);

   phs.close_read();

   std::cout << "PHS:: using cell and symmetry: "
	     << cell.descr().a() << " "
	     << cell.descr().b() << " "
	     << cell.descr().c() << " "
	     << clipper::Util::rad2d(cell.descr().alpha()) << " "
	     << clipper::Util::rad2d(cell.descr().beta())  << " "
	     << clipper::Util::rad2d(cell.descr().gamma()) << " "
	     << single_quote(sg.symbol_hm()) << std::endl;

   std::cout << "PHS:: number of reflections: " << mydata.num_reflections()
	     << "\n";

   fphidata.update();

   int ncount = 0; 
   clipper::HKL_info::HKL_reference_index hri;
//    for (hri=myfsig.first(); !hri.last(); hri.next()) {
//       if (ncount < 300) 
// 	 std::cout << " PHS fsigf: " << hri.hkl().h() << " "
// 		   << hri.hkl().k() << " "
// 		   << hri.hkl().l() << " " << myfsig[hri].f() << " "
// 		   << (myfsig[hri].sigf()) << std::endl;
//       ncount++;
//    }

   ncount = 0; 
//    for (hri=myphwt.first(); !hri.last(); hri.next()) {
//       if (ncount < 300) 
// 	 std::cout << " PHS myphwt: " << hri.hkl().h() << " " << hri.hkl().k() << " "
// 		   << hri.hkl().l() << " " << myphwt[hri].fom() << " "
// 		   << clipper::Util::rad2d(myphwt[hri].phi()) << std::endl;
//       ncount++;
//    }

  fphidata.compute(myfsig, myphwt, 
 		    clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

//    for (int i=0; i<10; i++) {
//       std::cout << "checking phi weight: " << i << " " << myphwt[i].phi() << "  "
//              << myphwt[i].fom() << std::endl;
//        std::cout << "checking f    sigf: " << i << " " << myfsig[i].f() << "   "
// 		 << myfsig[i].sigf() << std::endl;
//        std::cout << "checking missing: " << i << " " << myfsig[i].missing() << " "
// 		 << myphwt[i].missing() << " " << fphidata[i].missing() << std::endl;
//        // << " " << fphidata[i].phi() <<
//    }
   
  if (max_xmaps == 0) {
     
     xmap_list = new clipper::Xmap<float>[10];
     xmap_is_filled   = new int[10];
     xmap_is_diff_map = new int[10];
     contour_level    = new float[10];
   }

   max_xmaps++; 

  std::string mol_name = phs_filename; 

  initialize_map_things_on_read_molecule(mol_name, 0, 0); // not diff map

  cout << "initializing map..."; 
  xmap_list[0].init(mydata.spacegroup(), 
		    mydata.cell(), 
		    clipper::Grid_sampling(mydata.spacegroup(),
					   mydata.cell(), 
					   mydata.resolution()));
  cout << "done."<< endl;

  std::cout << "PHS:: debug:: " << mydata.spacegroup().symbol_hm() << " " 
	    << mydata.cell().descr().a() << " " 
	    << mydata.cell().descr().b() << " " 
	    << mydata.cell().descr().c() << " " 
	    << clipper::Util::rad2d(mydata.cell().descr().alpha()) << " " 
	    << clipper::Util::rad2d(mydata.cell().descr().beta ()) << " " 
	    << clipper::Util::rad2d(mydata.cell().descr().gamma()) << std::endl;

  std::cout << "PHS:: debug:: n_reflections: " << mydata.num_reflections()
		   << std::endl;

  ncount = 0; 
  // clipper::HKL_info::HKL_reference_index hri;
//   for (hri=fphidata.first(); !hri.last(); hri.next()) {
//      if (ncount < 300) 
// 	std::cout << " PHS fphi: " << hri.hkl().h() << " " << hri.hkl().k() << " "
// 		  << hri.hkl().l() << " " << fphidata[hri].f() << " "
// 		  << clipper::Util::rad2d(fphidata[hri].phi()) << std::endl;
//      ncount++;
//   } 


//   cout << "Map Grid (from phs file)..." 
//        << xmap_list[0].grid_sampling().format()
//        << endl;  

  cout << "doing fft..." ; 
  xmap_list[0].fft_from( fphidata );                  // generate map
  cout << "done." << endl;

  mean_and_variance<float> mv = map_density_distribution(xmap_list[0],0);

  cout << "Mean and sigma of map from PHS file: " << mv.mean 
       << " and " << sqrt(mv.variance) << endl;

  // fill class variables
  map_mean_ = mv.mean;
  map_sigma_ = sqrt(mv.variance);

  xmap_is_diff_map[0] = 0; 
  xmap_is_filled[0] = 1; 
  contour_level[0] = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);

  graphics_info_t g; 
  // update_map_colour_menu_manual(g.n_molecules, name_.c_str()); 
  // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str());

  g.scroll_wheel_map = g.n_molecules; // change the current scrollable map.

  std::cout << "updating map..." << std::endl;
  update_map();
  std::cout << "done updating map..." << std::endl;

  // how do we restore this map?
  save_state_command_strings_.push_back("read-phs-and-make-map-using-cell-symm");
  save_state_command_strings_.push_back(single_quote(phs_filename));
  save_state_command_strings_.push_back(single_quote(sg.symbol_hm()));
  save_state_command_strings_.push_back(g.float_to_string(cell.descr().a()));
  save_state_command_strings_.push_back(g.float_to_string(cell.descr().b()));
  save_state_command_strings_.push_back(g.float_to_string(cell.descr().c()));
  save_state_command_strings_.push_back(g.float_to_string(clipper::Util::rad2d(cell.descr().alpha())));
  save_state_command_strings_.push_back(g.float_to_string(clipper::Util::rad2d(cell.descr().beta())));
  save_state_command_strings_.push_back(g.float_to_string(clipper::Util::rad2d(cell.descr().gamma())));

  return g.n_molecules;
}
 

// export the molecule in atom_selection_container_t atom_sel;
// 
int
molecule_class_info_t::export_coordinates(std::string filename) const { 

   //
   int err = atom_sel.mol->WritePDBASCII((char *)filename.c_str()); 
   
   if (err) { 
      std::cout << "WARNING:: export coords: There was an error in writing "
		<< filename << std::endl; 
      std::cout << GetErrorDescription(err) << std::endl;
      graphics_info_t g;
      std::string s = "ERROR:: writing coordinates file ";
      s += filename;
      g.statusbar_text(s);
   } else {
      std::string s = "INFO:: coordinates file ";
      s += filename;
      s += " saved successfully";
      graphics_info_t g;
      g.statusbar_text(s);
   } 
   return err;
}

// Perhaps this should be a util function?
CMMDBManager
*molecule_class_info_t::get_residue_range_as_mol(const std::string &chain_id,
						 int resno_start,
						 int resno_end) const {

   int imod = 1;

   CMMDBManager *mol_new = new CMMDBManager;
   CModel *model_new = new CModel;
   CChain *chain_new = new CChain;
   
   realtype cell[6];
   realtype vol;
   int orthcode;
   char *spacegroup_str = atom_sel.mol->GetSpaceGroup();
   atom_sel.mol->GetCell(cell[0], cell[1], cell[2],
			 cell[3], cell[4], cell[5],
			 vol, orthcode); 
   mol_new->SetCell(cell[0], cell[1], cell[2],
		    cell[3], cell[4], cell[5], orthcode);
   mol_new->SetSpaceGroup(spacegroup_str);

   CModel *model_p = atom_sel.mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (std::string(chain_p->GetChainID()) == chain_id) { 
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p->GetSeqNum() >= resno_start) {
	       if (residue_p->GetSeqNum() <= resno_end) {
		  CResidue *res_new =
		     coot::util::deep_copy_this_residue(residue_p, "", 1);
		  chain_new->AddResidue(res_new);
	       }
	    }
	 }
      }
   }

   model_new->AddChain(chain_new);
   mol_new->AddModel(model_new);
   mol_new->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   mol_new->FinishStructEdit();
   return mol_new;
} 


std::string
molecule_class_info_t::make_symm_atom_label_string(PCAtom atom, symm_trans_t symm_trans) const {

   std::string s = make_atom_label_string(atom, 0);
   s += " ";
   s += symm_trans.str(graphics_info_t::symmetry_atom_labels_expanded_flag);
   return s;
}

// This took one minute to write.
//
// And then another 2 to add the alt location code (and test it).
// 
std::string
molecule_class_info_t::make_atom_label_string(PCAtom atom, int brief_atom_labels_flag) const {

   char *chain_id  = atom->GetChainID();
   char *res_name  = atom->GetResName();
   int   res_no    = atom->GetSeqNum();
   char *atom_name = atom->name;
   char *ins_code  = atom->GetInsCode();

   // format: atom_name/res_no res_name/chain_id
   // new format: atom_name,alt_conf/res_no res_name/chain_id if altconf != ""
   graphics_info_t g;

   std::string s(atom_name);
   std::string alt_loc(atom->altLoc);
   if (alt_loc != "") {
      int slen = s.length();
      if (slen > 0) {
	 if (s[slen-1] == ' ') {
	    s = s.substr(0,slen-1) + ",";
	 } else {
	    s += ",";
	 }
      } else {
	 s += ",";
      }
      s += alt_loc;
   }

   if (brief_atom_labels_flag) {
      s += g.int_to_string(res_no);
      if (strlen(ins_code) > 0) {
	 s += ins_code;
	 s += " ";
      }
      s += chain_id;
   } else { 
      s += "/";
      s += g.int_to_string(res_no);
      s += ins_code;
      s += " ";
      s += res_name;
      s += "/";
      s += chain_id;
   }

   return s;
}

// Don't use this function.
// 
// For labelling from guile and replacing coordinates and others.
// 
// Return -1 on failure to find match
int 
molecule_class_info_t::atom_spec_to_atom_index(std::string chain, int resno, 
					       std::string atom_name) const { 

   int iatom_index = -1; 
   int selHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(selHnd, 0, (char *) chain.c_str(), 
			    resno, "*", // start, insertion code
			    resno, "*", // end, insertion code
			    "*", // residue name
			    (char *) atom_name.c_str(),
			    "*", // elements
			    "*"); // alt locs

   int nSelAtoms;
   PPCAtom local_SelAtom; 
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "Sorry (atom_spec_to_atom_index): Could not find " << atom_name << "/"
		<< resno << "/" << chain << " in this molecule: ("
		<<  imol_no << ") " << name_ << std::endl; 

      // debug:
      selHnd = atom_sel.mol->NewSelection();
      
      atom_sel.mol->SelectAtoms(selHnd, 0, "*",
				ANY_RES, "*", // start, insertion code
				ANY_RES, "*", // end, insertion code
				"*", // residue name
				(char *) atom_name.c_str(),
				"*", // elements
				"*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

      std::cout << "There were " << nSelAtoms << " atoms with resno "
		<< resno << std::endl;

      
   } else {
      // compare pointers
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
	    iatom_index = i;
	    break; 
	 }
      }
   } 
   return iatom_index; 
}			    
			    

// return -1 on no atom found.
int
molecule_class_info_t::full_atom_spec_to_atom_index(const std::string &chain,
						    int resno,
						    const std::string &insertion_code,
						    const std::string &atom_name,
						    const std::string &alt_conf) const {

   int iatom_index = -1; 
   int selHnd = atom_sel.mol->NewSelection();
   int idx = 0;

   atom_sel.mol->SelectAtoms(selHnd, 0, (char *) chain.c_str(), 
			    resno, (char *) insertion_code.c_str(), // start, insertion code
			    resno, (char *) insertion_code.c_str(), // end, insertion code
			    "*", // residue name
			    (char *) atom_name.c_str(),
			    "*", // elements
			    (char *) alt_conf.c_str()); // alt locs

   int nSelAtoms;
   PPCAtom local_SelAtom; 
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

//    std::cout << chain << " " << resno << " " << insertion_code << " " 
// 	     << atom_name << " " << alt_conf << " finds " << nSelAtoms <<  " atoms\n";

   if (nSelAtoms == 0) { 

      std::cout << "Sorry (full_atom_spec_to_atom_index) Could not find "
		<< atom_name << "," << "\"" << alt_conf  << "\"" << "/"
		<< resno << insertion_code << "/" << chain << " in this molecule: ("
		<<  imol_no << ") " << name_ << std::endl; 
      // debug:
      int selHnd2 = atom_sel.mol->NewSelection();
      
      atom_sel.mol->SelectAtoms(selHnd2, 0, 
				(char *) chain.c_str(),
				resno, "*", // start, insertion code
				resno, "*", // end, insertion code
				"*", // residue name
				"*", // atom name
				"*", // elements
				"*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd2, local_SelAtom, nSelAtoms);

      std::cout << "There were " << nSelAtoms << " atoms with in that residue:\n";
      for (int i=0; i<nSelAtoms; i++) { 
	 std::cout << "      " << local_SelAtom[i] << "\n";
      }

      atom_sel.mol->DeleteSelection(selHnd2);

   } else { 

      if (nSelAtoms != 1) {
	 // the wildcard atom selection case "*HO2"
	 short int found = 0;
	 for (int i=0; i<nSelAtoms; i++) { 
	    if (std::string(local_SelAtom[i]->GetChainID()) == chain) { 
	       if (local_SelAtom[i]->residue->seqNum == resno) { 
		  if (std::string(local_SelAtom[i]->GetInsCode()) == insertion_code) { 
		     if (std::string(local_SelAtom[i]->name) == atom_name) { 
			if (std::string(local_SelAtom[i]->altLoc) == alt_conf) { 
			   found = 0;
			   idx = i;
			   break;
			} 
		     }
		  }
	       }
	    }
	 } 
      }
      
      int iatom_index_udd = -1;
      int ic;
      if (local_SelAtom[idx]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
	 iatom_index_udd = ic;
      }
      iatom_index = iatom_index_udd;
   }
   return iatom_index; 
}
//       // compare pointers
//       for (int i=0; i<atom_sel.n_selected_atoms; i++) {
// 	 if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
// 	    iatom_index = i;
// 	    break; 
// 	 }
//       }
	 // 	 if (iatom_index != iatom_index_udd) { 
	 // 	    std::cout << "ERROR: atom indexes (UDD test) dont match "
	 // 		      << iatom_index << " " << iatom_index_udd << std::endl;
	 // 	 }

// Another attempt at that:
// 
// find "by hand" the atom with the given characteristics in
// the atom selection.
//
// return -1 if atom not found.
// Note we have to search for " CA " etc
// 
int molecule_class_info_t::atom_index(const char *chain_id, int iresno, const char *atom_id) {

   int n = atom_sel.n_selected_atoms;
   for (int i=0; i<n; i++) {
      if ( ( ! strcmp(atom_id,atom_sel.atom_selection[i]->name) ) &&
	   (atom_sel.atom_selection[i]->residue->seqNum == iresno)  &&
	   ( ! strcmp(chain_id,atom_sel.atom_selection[i]->residue->GetChainID()) )
	   ) {
	 return i;
      }
   }

   return -1; 
}

int 
molecule_class_info_t::atom_index_first_atom_in_residue(const std::string &chain_id,
		                  		        int iresno,
				                        const std::string &ins_code) {

   int index = -1; // failure
   int selHnd = atom_sel.mol->NewSelection();
   int nSelResidues;
   PPCResidue SelResidues;
   atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 1,
			(char *) chain_id.c_str(), 
			iresno, (char *) ins_code.c_str(),
			iresno, (char *) ins_code.c_str(),
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			SKEY_NEW // selection key
			);
   atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
   if (nSelResidues > 0) {
      int ic = -1;
      for (int ires=0; ires<nSelResidues; ires++) {
	 int natoms;
	 PPCAtom residue_atoms;
	 SelResidues[ires]->GetAtomTable(residue_atoms, natoms);
	 for (int iatom=0; iatom<natoms; iatom++) {
	    if (residue_atoms[iatom]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
	       index = ic;
	       break;
	    }
	 }
	 if (index > -1)
	    break;
      }
   }
   atom_sel.mol->DeleteSelection(selHnd);
   // std::cout << "DEBUG:: atom_index_first_atom_in_residue returns " << index << std::endl;
   return index;

} 

// Put the regularization results back into the molecule:
//
//// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
molecule_class_info_t::replace_coords(const atom_selection_container_t &asc) {

   int idx = -1;
   CAtom *atom;
   CAtom *mol_atom;
   int n_atom = 0;
   int tmp_index;

   make_backup();

   
//    std::cout << "DEBUG:: --------------- replace_coords replacing " << asc.n_selected_atoms
//  	     << " atoms " << std::endl;

   for (int i=0; i<asc.n_selected_atoms; i++) {
      atom = asc.atom_selection[i];
//       idx = atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
// 				    atom->residue->seqNum,
// 				    std::string(atom->name));
      if (asc.UDDOldAtomIndexHandle >= 0) { // OK for fast atom indexing
	 if (atom->GetUDData(asc.UDDOldAtomIndexHandle, tmp_index) == UDDATA_Ok) {
	    if (tmp_index >= 0) { 
	       // std::cout << "successfully found old atom index" << std::endl;
	       idx = tmp_index;
	    } else {
	       // This shouldn't happen.
	       std::cout << "Good Handle, bad index found for old atom: specing" << std::endl;
	       idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
						  atom->residue->seqNum,
						  std::string(atom->GetInsCode()),
						  std::string(atom->name),
						  std::string(atom->altLoc));
	    }
	 } else { 
	    std::cout << "ERROR:: non-bad handle (" << asc.UDDOldAtomIndexHandle 
		      <<  "), bad GetUDData for this atom " << std::endl;
	 } 
      } else {
//  	 std::cout << "DEBUG:: asc.UDDOldAtomIndexHandle is " 
//  		   << asc.UDDOldAtomIndexHandle << " using full atom spec to atom index...";
	    
	 idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
					    atom->residue->seqNum,
					    std::string(atom->GetInsCode()),
					    std::string(atom->name),
					    std::string(atom->altLoc));
	 // std::cout << "DEBUG:: idx: " << idx << "\n";
      }
      if (idx >= 0) {
	 n_atom++;
	 mol_atom = atom_sel.atom_selection[idx];
	 float atom_occ = atom->occupancy;
	 // if this is a shelx molecule, then we don't change
	 // occupancies this way.  We do it by changing the FVAR
	 if (is_from_shelx_ins_flag) { 
	    atom_occ = mol_atom->occupancy;

	    // OK, one more go.  We have an occupancy of 31 or -31
	    // say.  Now, the alt conf atoms has been immmediately
	    // added with the old occupancy for the actual FVAR number
	    // - this happens before we get to twiddle the occupancy
	    // slider.  So here we have to find out the index of the
	    // replaced atom and set it's fvar to whatever the slider
	    // value had been set to.

	    int fvar_number = coot::ShelxIns::shelx_occ_to_fvar(atom_occ);
	    if (fvar_number > 1) { 
// 	       std::cout << "DEBUG:: replace_coords: setting fvar number "
// 			 <<  fvar_number << " (generated from occ " << atom_occ << ") to "
// 			 << graphics_info_t::add_alt_conf_new_atoms_occupancy << std::endl;
	       shelxins.set_fvar(fvar_number, graphics_info_t::add_alt_conf_new_atoms_occupancy);
	    }
	    
	    mol_atom->SetCoordinates(atom->x,
				     atom->y,
				     atom->z,
				     atom_occ,
				     mol_atom->tempFactor);
	 } else { 
	    mol_atom->SetCoordinates(atom->x,
				     atom->y,
				     atom->z,
				     atom_occ,
				     mol_atom->tempFactor);
	 }

	 // similarly we adjust occupancy if this is not a shelx molecule
	 if (! is_from_shelx_ins_flag) 
	    adjust_occupancy_other_residue_atoms(mol_atom, mol_atom->residue, 0);
	 // std::cout << atom << " coords replace " << idx << " " << mol_atom << std::endl;
      } else {
	 std::cout << "ERROR:: bad atom index in replace_coords replacing atom: "
		   << atom << std::endl;
      } 
   }
   std::cout << n_atom << " atoms updated." << std::endl;
   have_unsaved_changes_flag = 1; 

   make_bonds_type_checked();
}


// This relies on the mol of the asc being different to the mol of the
// atom_sel.
//
// If it is the same, then error and do nothing.
//
// Typically, this is called by fit terminal residue, which has its
// mol created from a pcmmdbmanager() from the molecule of the
// residue, so this is fine in this case.
// 
void
molecule_class_info_t::insert_coords(const atom_selection_container_t &asc) {

   // for each residue in the asc, do a InsResidue into its chain:

   if (! (atom_sel.n_selected_atoms > 0) ) {
      std::cout << "ERROR: Can't insert_coords this asc  - no atoms in molecule!\n"; 
   } else {

      // pointer comparison
      if (asc.mol == atom_sel.mol) {

	 std::cout << "ERROR:: matching asc.mol and atom_sel.mol in insert_coords\n";
	 std::cout << "ERROR:: new algorithm required\n";

      } else {

         make_backup();
	 insert_coords_internal(asc);

      }
   }
} 

void 
molecule_class_info_t::insert_coords_internal(const atom_selection_container_t &asc) {

   // run over each chain, residue of the asc (if terminal residue
   // fit only one chain, one residue, of course).

   short int inserted = 0; // not inserted yet
   CChain *asc_chain;
   int asc_n_chains = asc.mol->GetNumberOfChains(1);
   for (int i_asc_chain=0; i_asc_chain<asc_n_chains; i_asc_chain++) {
      asc_chain = asc.mol->GetChain(1,i_asc_chain);
      int nres_asc = asc_chain->GetNumberOfResidues();
//       std::cout << "DEBUG:: There are " << nres_asc << " residues in "
// 		<< "asc_chain (chain id: " << asc_chain->GetChainID()
// 		<< ")." << std::endl;

      int udd_atom_index = asc.UDDAtomIndexHandle; 

      for (int ires_asc=0; ires_asc<nres_asc; ires_asc++) {
	 CResidue *asc_residue = asc_chain->GetResidue(ires_asc);

	 // Now find the corresponding chain in our atom_sel.mol:

	 CChain *chain;
	 int n_chains = atom_sel.mol->GetNumberOfChains(1);
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {

	    chain = atom_sel.mol->GetChain(1,i_chain);

	    // test chains
	    std::string asc_chain_str(asc_chain->GetChainID());
	    std::string mol_chain_str(    chain->GetChainID());
	    if (asc_chain_str == mol_chain_str) {

	       // insert that residue!
	       CResidue *res = 
		  coot::deep_copy_this_residue(asc_residue, "", 1, udd_atom_index);
// 	       std::cout  << "DEBUG:: inserting residue in chain "
// 			  << mol_chain << " residue number "
// 			  << asc_residue->GetSeqNum()
// 			  << " i_chain = " << i_chain
// 			  << " ires_asc = " << ires_asc
// 			  << std::endl;

	       int serial_number =
		  find_serial_number_for_insert(asc_residue->GetSeqNum(),
						mol_chain_str);

// 	       std::cout << "DEBUG:: returned serial_number: " << serial_number
// 			 << std::endl;

	       if (serial_number != -1) { 
		  chain->InsResidue(res, serial_number);
		  inserted = 1;
	       } else { 
		  std::cout << "insert coord add residue\n";
		  chain->AddResidue(res);
	       }
	    }
	    //if (inserted) break;
	 }
	 //if (inserted) break;
      }
      //if (inserted) break;
   }
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
}



void 
molecule_class_info_t::insert_coords_change_altconf(const atom_selection_container_t &asc) {

   // There are 2 things we want to do here.
   //
   // For matching atoms, if there is only one atom that matches the
   // spec (appart from the altconf), we move the original atoms
   // altconf to "A".  If they are not "" we leave don't change the
   // altconf.
   // (see table in // molecule_class_info_t::make_new_alt_conf()
   // [molecule-class-info-other.cc]
   // 

   // The second thing is the change the occ of the existing atoms:
   // 1 atom:  new_occ_existing_atom = 1 - occ_new_atom;
   // general: new_occ_existing_atom = current_occ_existing_atom - occ_new_atom/n_similar_atoms
   // where n_similar_atoms is the number of alt confs we have for that atom.
   std::cout << "DEBUG:: ----------------- in insert_coords_change_altconf ------ " << std::endl;
   std::cout << "DEBUG:: IN insert_coords_change_altconf" << std::endl;
   make_backup();
   
   // OK if we were from a shelx ins file, then we have to create a
   // new FVAR for this new alt conf thingy.
   int shelx_occ_fvar_number = -1; 
   if (is_from_shelx_ins_flag) {
      // OK, what was the occupancy?
      if (asc.n_selected_atoms > 0) {
	 float occ = asc.atom_selection[0]->occupancy;
// 	 std::cout << "DEBUG:: IN insert_coords_change_altconf adding fvar "
// 		   << occ << std::endl;
	 shelx_occ_fvar_number = 10 * shelxins.add_fvar(occ); // FVAR 1 is not written
	 // to SHELX file so if shelx ins file has 1 FVAR value, then we've just
	 // created shelx FVAR 3.
	 shelx_occ_fvar_number += 1;  // so thats (e.g.) 1 x the 20th FVAR
      } 
   }
   
   char *chain_id;
   char *atom_name;
   int  resno;
   float occ; 
   CAtom *at;
   for(int i=0; i<asc.n_selected_atoms; i++) { 
      at = asc.atom_selection[i];
      chain_id  = at->GetChainID();
      atom_name = at->GetAtomName();
      resno     = at->GetSeqNum();
      occ       = at->occupancy;
      char *inscode = at->GetInsCode();

      // Now find that atom the corresponding atom (with altconf "" in
      // the original atoms).  We skip over atoms that don't have
      // altconf "").

      int selHnd = atom_sel.mol->NewSelection();
      atom_sel.mol->SelectAtoms(selHnd, 0, chain_id,
				resno, inscode,
				resno, inscode,
				"*", // residue name
				atom_name,
				"*", 
				"*"); // alt-loc
      int nSelAtoms;
      PPCAtom local_SelAtom;
      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

      if (nSelAtoms == 0) { 
	 // debugging
// 	       std::cout << "add alt conf: skipping " << atom_name << "/"
// 		<< resno << "/" << chain_id << " in this molecule: ("
// 		<<  imol_no << ") " << name_ << std::endl; 
      } else {

	 // Let's deal with every atom in this residue that that has
	 // the same atom name.

	 // first check for atoms with alt conf "" and move them over if needed
	 for (int iat=0; iat<nSelAtoms; iat++) {
	    std::string current_alt_conf = local_SelAtom[iat]->altLoc;
	    if (current_alt_conf == "") { 
	       std::string new_alt_conf("A");
	       // force it down the atom's throat :)
	       strncpy(local_SelAtom[0]->altLoc, new_alt_conf.c_str(), 2);
	    }
	 }

	 if (shelx_occ_fvar_number == -1) {
	    // i.e. was not from a shelx ins file (normal case):
	    // 
	    // now stuff occupancies in...
	    for (int iat=0; iat<nSelAtoms; iat++) {
	       // local_SelAtom[0]->occupancy = 1.0 - occ; // complemetary (1+1 atom case)
	       local_SelAtom[iat]->occupancy -= occ/float(nSelAtoms); // complemetary (general case)
	    }

	 } else {

	    // This is the SHELX case:
	    if (nSelAtoms > 1) {
	       // so we added an alt conf to a residue that already
	       // has an alt conf.  This involves messing with SUMP.
	       // Let's not handle that for now.
	       std::cout << "WARNING:: SHELX occupancy handler under-resourced on handling "
			 << at << std::endl;
	    }  else {
	       local_SelAtom[0]->occupancy = -shelx_occ_fvar_number;
	    }
	 } 
      }
      atom_sel.mol->DeleteSelection(selHnd);
   } 
   insert_coords_atoms_into_residue_internal(asc, shelx_occ_fvar_number);

}

// In this instance, we don't want to install a whole residue, we want
// to install atoms in this residue (alt conf B) into a a atom_sel mol
// residue that contains (say) "" and "A".
//
// -1 is passed as shelx_occ_fvar_number if this atom was not from a
// SHELX ins file.  If shelx_occ_fvar_number > 1, then use this as the
// new atom's occupancy.
// 
void 
molecule_class_info_t::insert_coords_atoms_into_residue_internal(const atom_selection_container_t &asc,
								 int shelx_occ_fvar_number) {

   char *chain_id;
   char *atom_name;
   int  resno;
   CAtom *at;
   int afix_handle_this_mol = -1;
   int afix_handle_intermediate_mol = -1;

   afix_handle_this_mol    = atom_sel.mol->GetUDDHandle(UDR_ATOM, "shelx afix");
   afix_handle_intermediate_mol = asc.mol->GetUDDHandle(UDR_ATOM, "shelx afix");

   std::cout << "DEBUG in insert_coords_atoms_into_residue_internal afix handles:"
	     << afix_handle_this_mol << " " << afix_handle_intermediate_mol << std::endl;
   
   for(int i=0; i<asc.n_selected_atoms; i++) { 
      at = asc.atom_selection[i];
      chain_id  = at->GetChainID();
      atom_name = at->GetAtomName();
      resno     = at->GetSeqNum();

      // Now find the corresponding residue in atom_sel.mol;

      int selHnd = atom_sel.mol->NewSelection();
      int nSelResidues;
      PPCResidue SelResidues;
      atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 1,
			   chain_id, 
			   resno, "*",
			   resno, "*",
			   "*",  // residue name
			   "*",  // Residue must contain this atom name?
			   "*",  // Residue must contain this Element?
			   "*",  // altLocs
			   SKEY_NEW // selection key
			   );
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      
      if (nSelResidues != 1) { 
	 std::cout << "ERROR:: something broken in residue selection in ";
	 std::cout << "insert_coords_atoms_into_residue_internal: got " << nSelResidues
		   << " residues." << std::endl;
      } else { 
	 CAtom *t = new CAtom;
	 t->Copy(at);
	 SelResidues[0]->AddAtom(t);
	 // if these coords were from a shelx ins file, then we are
	 // passed the free varible number (FVAR) for this atom's
	 // occupancy.
	 if (shelx_occ_fvar_number > 1)
	    t->occupancy = shelx_occ_fvar_number;

	 // If these coords were from a shelx ins file, then we need
	 // to copy the AFIX numbers too:
	 int afix_number; // set by getUDD
	 if (afix_handle_intermediate_mol > -1) {
	    int ierr = at->GetUDData(afix_handle_intermediate_mol, afix_number);
	    if (ierr == UDDATA_Ok) {
	       if (afix_handle_this_mol > -1) {
		  t->PutUDData(afix_handle_this_mol, afix_number);
	       } else {
		  std::cout << "ERROR:: bad afix handle for this molecule in "
			    << "insert_coords_atoms_into_residue_internal"
			    << afix_handle_this_mol << " " << at << std::endl;
	       }

 	    } else {
	       if (is_from_shelx_ins_flag) 
		  std::cout << "ERROR:: attempt to get UDD afix number from "
			    << "intermediate molecule failed " << at << std::endl;
 	    } 
	 } else {
	    std::cout << "ERROR:: bad afix handle for intermediate molecule in "
		      << "insert_coords_atoms_into_residue_internal"
		      << afix_handle_intermediate_mol << " " << at << std::endl;
	 } 
      } 
      atom_sel.mol->DeleteSelection(selHnd);
   } 
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
} 

// We need to find the serial number of the residue after the residue
// we want to insert (i.e. the new residue will be inserted just
// before the residue whose serial number we return).
//
// return -1 on error.
int
molecule_class_info_t::find_serial_number_for_insert(int seqnum_new,
						     const std::string &chain_id) const {

   int iserial_no = -1;
   int current_diff = 999999;
   CChain *chain;
   int n_chains = atom_sel.mol->GetNumberOfChains(1);
   for (int i_chain=0; i_chain<n_chains; i_chain++) {
      
      chain = atom_sel.mol->GetChain(1,i_chain);
      
      // test chains
      std::string mol_chain(    chain->GetChainID());
      if (chain_id == mol_chain) {

	 // find the 
	 int nres = chain->GetNumberOfResidues();
	 
	 for (int ires=0; ires<nres; ires++) { // ires is a serial number
	    CResidue *residue = chain->GetResidue(ires);

	    // we are looking for the smallest negative diff:
	    // 
	    int diff = residue->GetSeqNum() - seqnum_new;

	    if ( (diff > 0) && (diff < current_diff) ) {
	       iserial_no = ires;
	       current_diff = diff;
	    }
	 }
      }
   }
   return iserial_no;
} 


// Put the regularization results back into the molecule:
//
// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
molecule_class_info_t::add_coords(const atom_selection_container_t &asc) {

   CAtom *atom;
   CAtom *mol_atom; // an atom already existing in mol
   int n_atom = 0;
   CChain *chain;

   // std::cout << "DEBUG:: ----------------- in add_coords ----------- " << std::endl;

   make_backup();
   
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idone = 0;
      atom = asc.atom_selection[i];
      // chain = atom->GetChain();
      
      // run over chains of the existing mol
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) { 
	 
	 chain = atom_sel.mol->GetChain(1,ichain);
	 std::string atom_chain_id(atom->GetChainID());
	 std::string mol_chain_id(chain->GetChainID());
	 
	 if (atom_chain_id == mol_chain_id) {

	    // int iseqno_at = atom->GetSeqNum();
	    int nres = chain->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) { 
	       PCResidue res = chain->GetResidue(ires);
	       if (res) { // fixes bug 030813, but best way?
		          // 030814, ah, we discover FinishStructEdit().
		  if (res->GetSeqNum() == atom->GetSeqNum()) { 
		     int natom = res->GetNumberOfAtoms();
		     for (int iat=0; iat<natom; iat++) {
			
			mol_atom = res->GetAtom(atom->GetAtomName());
			if (mol_atom) { // we got a match
			   // replace the coordinates then
			   
			   // This should not happen very often
			   std::cout << "add_coords: replacing " << mol_atom 
				     << " with new atom " << atom << std::endl;
			   mol_atom->SetCoordinates(atom->x,
						    atom->y,
						    atom->z,
						    mol_atom->occupancy,
						    mol_atom->tempFactor);
			   idone = 1;
			   break;
			   
			} else { 
			   
			   std::cout << "adding atom to existing residue " 
				     << atom << " (already has " 
				     << res->GetNumberOfAtoms() << " atoms)" 
				     << std::endl;
			   CAtom *new_atom = new CAtom;
			   new_atom->Copy(atom);
			   res->AddAtom(new_atom);
			   new_atom->occupancy = 1.0;
			   new_atom->tempFactor = 10.0;
			   // chain id:
			   std::cout << "setting chainid of this new atom from "
				     << new_atom->GetChainID() << " to : "
				     << atom->GetChainID() << std::endl;
			   new_atom->residue->chain->SetChainID(atom->GetChainID());
			   idone = 1;
			   n_atom++;
			   break;
			}
		     }
		  }
	       }
	       if (idone == 1) break;
	    } // residue loop
	 }
      }

      if (idone == 0) { 

	 std::cout << "adding whole residue triggered by atom " 
		   << atom << std::endl;
	 std::cout << "     with element " << atom->element << std::endl;

	 // in this bit of code, atom is an atom from the asc and
	 // atom_p is a new atom that we are adding to a new residue
	 // (that we are adding to an existing chain (that we do a
	 // lookup to find)).
	 // 
	 CResidue *res_p = new CResidue;
	 CAtom *atom_p = new CAtom;
	 // CChain *chain_p = atom_sel.mol->GetChain(1,0);
	 PCChain chain_p = atom_sel.mol->GetChain(1,atom->GetChainID());
	 chain_p->AddResidue(res_p);
	 atom_p->SetAtomName(atom->name);
	 atom_p->SetCoordinates(atom->x, atom->y, atom->z,
				atom->occupancy, atom->tempFactor);
	 atom_p->SetElementName(atom->element);
	 res_p->AddAtom(atom_p);
	 res_p->seqNum = atom->GetSeqNum();
	 res_p->SetResID(atom->residue->name,
			 atom->GetSeqNum(),
			 atom->GetInsCode());
	 
	 // add to end:
	 atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	 atom_sel.mol->FinishStructEdit();
      }
   }

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;

   // now regenerate the atom_selection
   //

   // Uncomment when we have fixed the bug, this seems not to be it.
   // 
   // clear out the old
   // if (atom_sel.atom_selection != NULL)
   // delete [] atom_sel.atom_selection;
   
   // and in with the new:
   int selHnd = atom_sel.mol->NewSelection(); 
   atom_sel.mol->SelectAtoms(selHnd, 0,"*",ANY_RES,"*",ANY_RES,
			    "*","*",  // EndInsertionCode, RNames
			    "*","*",  // ANames, Elements
			    "*" );    // Alternate locations.

   int old_n_atoms = atom_sel.n_selected_atoms;
   atom_sel.mol->GetSelIndex(selHnd, 
			     atom_sel.atom_selection, 
			     atom_sel.n_selected_atoms);

   std::cout << "INFO:: old n_atoms: " << old_n_atoms << " new: " 
	     << atom_sel.n_selected_atoms << std::endl;

   debug_selection();
   have_unsaved_changes_flag = 1; 

   make_bonds_type_checked();
   // std::cout << "DEBUG:: ---------------- done add_coords ----------- " << std::endl;
} 



void
molecule_class_info_t::close_yourself() {

   // Deletion causing problems on application closure

   short int was_map = 0;
   short int was_coords = 0;

   name_ = ""; // not "Baton Atoms" or anything now.
   
   if (atom_sel.n_selected_atoms > 0)
      was_coords = 1;

   if (xmap_is_filled[0])
      was_map = 1;

   // delete from display manager combo box
   // 
   graphics_info_t g;
   GtkWidget *display_control_window = g.display_control_window();
   // 
   if (display_control_window) { // is being displayed
      std::string display_frame_name = "display_mol_frame_";
      if (was_map)
	 display_frame_name = "display_map_frame_";
      display_frame_name += g.int_to_string(imol_no);
      // std::cout << "DEBUG:: looking up " << display_frame_name << std::endl;
      GtkWidget *display_frame = lookup_widget(display_control_window,
					       display_frame_name.c_str());
      if (display_frame)
	 gtk_widget_destroy(display_frame);
   } else {
      // std::cout << "close: display_control_window is not active" << std::endl;
   }

   if (was_coords) { 
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
      delete atom_sel.mol;
      // atom_sel.mol = 0; done later
   }

   if (was_map) {
      fc_skeleton_draw_on = 0; // turn off the skeleton
      delete [] xmap_list;
      xmap_list = NULL;
      xmap_is_filled[0] = 0;
      max_xmaps = 0;
   }

   bonds_box.clear_up();
   // symmetry_bonds_box?  (It is a vector of pairs)
   drawit = 0;
   drawit_for_map = 0;

   // Do these whatever the molecule type:
   atom_sel.n_selected_atoms = 0;
   atom_sel.atom_selection = NULL;
   atom_sel.mol = NULL;

   // set the number of labels to 0.
   n_labelled_symm_atoms = 0;
   n_labelled_atoms = 0;
   // 20060104
   delete [] labelled_symm_atom_symm_trans_;
   delete [] labelled_symm_atom_index_list;   
   delete [] labelled_atom_index_list;
   labelled_symm_atom_symm_trans_ = 0;
   labelled_symm_atom_index_list = 0;
   labelled_atom_index_list = 0;
   //
   // gl widget redraw is done in close_molecule
}


// Return the atom index of the "next" atom
// -1 on failure.
int
molecule_class_info_t::intelligent_next_atom(const std::string &chain,
					     int resno,
					     const std::string &atom_name,
					     const std::string &ins_code) {

   // First we try adding 1 to the residue number
   // If that fails,
   // 
   // Get the atom index of the current atom and step through
   // atom_selection until the residue number that atom is not the
   // same as the starting atom

   int i_atom_index = -1; // failure initially.

   // std::cout << "intelligent_next_atom: start" << std::endl;
   
   if (atom_sel.n_selected_atoms <= 0 || atom_sel.mol == NULL) { 

      std::cout << "ERROR:: trying to move to (next) atom of a closed molecule!\n";

   } else {
      int selHnd1 = atom_sel.mol->NewSelection();

      atom_sel.mol->SelectAtoms(selHnd1, 0, (char *) chain.c_str(), 
				resno+1, "*", // start, insertion code
				resno+1, "*", // end, insertion code
				"*", // residue name
				(char *) atom_name.c_str(),
				"*", // elements
				"*"); // alt locs

      int nSelAtoms;
      PPCAtom local_SelAtom; 
      atom_sel.mol->GetSelIndex(selHnd1, local_SelAtom, nSelAtoms);
      // std::cout << "intelligent_next_atom: nSelAtoms: " << nSelAtoms << std::endl;

      // Deleted at end.
      if (nSelAtoms > 0) {

	 // Ye Olde Compare Pointers senario:
	 for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	    if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
	       i_atom_index = i;
	       break; 
	    }
	 }
      } else { 
	 //       std::cout << "Intelligent:: Atom " << chain << "/" << resno << "/"
	 // 		<< atom_name << " not found " << std::endl;
	 //       std::cout << "Intelligent:: trying something else..." << std::endl;

	 int selHnd2 = atom_sel.mol->NewSelection();

	 // get the current atom index:
	 // 
	 atom_sel.mol->SelectAtoms(selHnd2, 0, (char *) chain.c_str(), 
				   resno, "*", // start, insertion code
				   resno, "*", // end, insertion code
				   "*", // residue name
				   (char *) atom_name.c_str(),
				   "*", // elements
				   "*"); // alt locs

	 atom_sel.mol->GetSelIndex(selHnd2, local_SelAtom, nSelAtoms);
	 if (nSelAtoms > 0) {

	    int icurrent = -1;
	    // Ye Olde Compare Pointers senario again:
	    for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	       if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
		  icurrent = i;
		  break; 
	       }
	    }

	    if (icurrent != -1) {
	    
	       // 	    std::cout << "Intelligent:: current index is " << icurrent << std::endl;
	       for (int i=icurrent; i<atom_sel.n_selected_atoms; i++) {
		  if (atom_sel.atom_selection[i]->GetSeqNum() != resno ||
		      std::string(atom_sel.atom_selection[i]->GetChainID()) != chain) {
		     // 		  std::cout << "Intelligent:: found new res at atom index "
		     // 			    << i << std::endl;
		     i_atom_index = i;

		     // So far, we have an atom in a residue, often this
		     // wil be the "N" atom of a peptide residue.
		     //
		     // Now, we want to move to the CA in that residue if
		     // we can.  So lets search the next 20 atoms in the
		     // list (should be only 1 required, actually to see
		     // if we can find the CA of that residue)

		     std::string chain_id = atom_sel.atom_selection[i]->GetChainID();
		     int ires_no =  atom_sel.atom_selection[i]->GetSeqNum();
		     std::string ca_name = " CA ";
		  
		     for (int in=i; in<atom_sel.n_selected_atoms && in<i+20; in++) {
			if (std::string(atom_sel.atom_selection[in]->GetChainID()) == chain_id &&
			    atom_sel.atom_selection[in]->GetSeqNum() == ires_no &&
			    std::string(atom_sel.atom_selection[in]->name) == ca_name) {
			   i_atom_index = in;
			   break;
			}
		     }
		     break; 
		  }
	       }
	    } else {
	       std::cout << "IMPOSSIBLE: Cannot happen in intelligent_next_atom\n";
	    }
	    // now give back the (selection) memory
	    atom_sel.mol->DeleteSelection(selHnd2);
	 } else {
	    // The initial atom was not found
	    // 
	    std::cout << "WARNING intelligent_next_atom: "
		      << "inital atom not found - giving up\n";
	 }
      }
      atom_sel.mol->DeleteSelection(selHnd1);
   }

   // std::cout << "Intelligent:: returning index " << i_atom_index << std::endl;
   return i_atom_index;
} 


// Return the atom index of the "next" atom
// -1 on failure.
int
molecule_class_info_t::intelligent_previous_atom(const std::string &chain,
						 int resno,
						 const std::string &atom_name,
						 const std::string &ins_code) {

   // First we try adding 1 to the residue number
   // If that fails,
   // 
   // Get the atom index of the current atom and step through
   // atom_selection until the residue number that atom is not the
   // same as the starting atom

   int i_atom_index = -1; // failure initially.

   if (! atom_sel.mol) {

      std::cout << "ERROR:: trying to move to (prev) atom of a closed molecule!\n";

   } else { 

   
      int selHnd = atom_sel.mol->NewSelection();

      atom_sel.mol->SelectAtoms(selHnd, 0, (char *) chain.c_str(), 
				resno-1, "*", // start, insertion code
				resno-1, "*", // end, insertion code
				"*", // residue name
				(char *) atom_name.c_str(),
				"*", // elements
				"*"); // alt locs

      int nSelAtoms;
      PPCAtom local_SelAtom; 
      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);
      if (nSelAtoms > 0) {

	 // Ye Olde Compare Pointers senario:
	 for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	    if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
	       i_atom_index = i;
	       break; 
	    }
	 }
      } else { 
	 std::cout << "Atom " << chain << "/" << resno << "/"
		   << atom_name << " not found " << std::endl;
	 std::cout << "trying something else..." << std::endl;

	 selHnd = atom_sel.mol->NewSelection();

	 // get the atom index:
	 // 
	 atom_sel.mol->SelectAtoms(selHnd, 0, (char *) chain.c_str(), 
				   resno, "*", // start, insertion code
				   resno, "*", // end, insertion code
				   "*", // residue name
				   (char *) atom_name.c_str(),
				   "*", // elements
				   "*"); // alt locs

	 atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);
	 if (nSelAtoms > 0) {

	    int icurrent = -1;
	    // Ye Olde Compare Pointers senario again:
	    for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	       if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
		  icurrent = i;
		  break; 
	       }
	    }

	    if (icurrent != -1) {
	       for (int i=icurrent; i>=0; i--) {
		  if (atom_sel.atom_selection[i]->GetSeqNum() != resno ||
		      std::string(atom_sel.atom_selection[i]->GetChainID()) != chain) {
		     i_atom_index = i;

		     // So far, we have an atom in a residue, often this
		     // wil be the last atom of a peptide residue.
		     //
		     // Now, we want to move to the CA in that residue if
		     // we can.  So lets search the back for 20 atoms in the
		     // list.
		     // 
		     std::string chain_id = atom_sel.atom_selection[i]->GetChainID();
		     int ires_no =  atom_sel.atom_selection[i]->GetSeqNum();
		     std::string ca_name = " CA ";
		  
		     for (int in=i; in>=0 && in<(i-20); in--) {
			if (std::string(atom_sel.atom_selection[in]->GetChainID()) == chain_id &&
			    atom_sel.atom_selection[in]->GetSeqNum() == ires_no &&
			    std::string(atom_sel.atom_selection[in]->name) == ca_name) {
			   std::cout << "Found a CA for that residue\n";
			   i_atom_index = in;
			   break;
			}
		     }
		     break; 
		  }
	       }
	    } else {
	       std::cout << "IMPOSSIBLE: Cannot happen in intelligent_previous_atom\n";
	    }
	    // now give back the (selection) memory
	    atom_sel.mol->DeleteSelection(selHnd);
	 } else {
	    // The initial atom was not found
	    // 
	    std::cout << "WARNING intelligent_previous_atom: "
		      << "inital atom not found - giving up\n";
	 }
      }
   }   
   return i_atom_index;
}

// If there is a CA in this residue then return the index of that
// atom, if not, then return the index of the first atom in the
// residue.
// 
// Return -1 on no atoms in residue.
//
int
molecule_class_info_t::intelligent_this_residue_atom(CResidue *res_p) const {

   PPCAtom residue_atoms;
   int nResidueAtoms;
   
   res_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " CA ") {
	 return atom_spec_to_atom_index(residue_atoms[i]->GetChainID(),
					residue_atoms[i]->GetSeqNum(),
					residue_atoms[i]->name);
      }
   }

   if (nResidueAtoms > 0) {
      return atom_spec_to_atom_index(residue_atoms[0]->GetChainID(),
				     residue_atoms[0]->GetSeqNum(),
				     residue_atoms[0]->name);
   }

   return -1;
}

// If there is a CA in this residue then return that atom (pointer)
// atom, if not, then return the index of the first atom in the
// residue.
// 
// Return NULL on no atoms in residue.
//
CAtom *
molecule_class_info_t::intelligent_this_residue_mmdb_atom(CResidue *res_p) const {
   PPCAtom residue_atoms;
   int nResidueAtoms;
   
   res_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " CA ") {
	 return residue_atoms[i];
      }
   }

   if (nResidueAtoms > 0) {
      return residue_atoms[0];
   }

   // failure
   return NULL;

}

// Return pointer to atom " CA ", or the first atom in the residue, or
// null (no residue or atoms error):
// 
CAtom *
molecule_class_info_t::atom_intelligent(const std::string &chain_id, int resno,
					const std::string &ins_code) const { 

   CAtom *at = NULL; 

   if (atom_sel.n_selected_atoms > 0) { 
      int selHnd = atom_sel.mol->NewSelection();
      PPCResidue SelResidue;
      int nSelResidues;

      char *ins_code_search = (char *) ins_code.c_str(); // bleugh (as usual)
      
      atom_sel.mol->Select (selHnd, STYPE_RESIDUE, 0,
			    (char *)chain_id.c_str(), 
			    resno, ins_code_search,
			    resno, ins_code_search,
			    "*",  // residue name
			    "*",  // Residue must contain this atom name?
			    "*",  // Residue must contain this Element?
			    "*",  // altLocs
			    SKEY_NEW // selection key
			    );

      atom_sel.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
      
      if (nSelResidues == 0) {
	 std::cout << "INFO:: No selected residues" << std::endl;
      } else {
	 
	 PPCAtom residue_atoms;
	 int nResidueAtoms;
	 SelResidue[0]->GetAtomTable(residue_atoms, nResidueAtoms);
	 if (nResidueAtoms == 0) {
	    std::cout << "INFO:: No atoms in residue" << std::endl;
	 } else {
	    short int found_it = 0;
	    for (int i=0; i<nResidueAtoms; i++) { 
	       if (std::string(residue_atoms[i]->name) == std::string(" CA ")) { 
		  at = residue_atoms[i];
		  found_it = 1;
		  break;
	       }
	    } 
	    if (! found_it) 
	       at = residue_atoms[0];
	 }
      }
      atom_sel.mol->DeleteSelection(selHnd); // Safe to put it here?
   }
   return at;
} 


// ----------------------------------------------------------------------
//               Pointer Atoms
// ----------------------------------------------------------------------
void
molecule_class_info_t::add_pointer_atom(coot::Cartesian pos) {

   int nchains = atom_sel.mol->GetNumberOfChains(1);

   if (nchains != 1) {
      std::cout << "failed to add (untyped) pointer atom to molecule " 
		<< name_ << std::endl;
      return;
   }

   make_backup();
      
   CChain *chain_p = atom_sel.mol->GetChain(1,0);
   
   std::string mol_chain_id(chain_p->GetChainID());
   // int ires_prev = chain_p->GetNumberOfResidues();
   int ires_prev = coot::util::max_resno_in_chain(chain_p).second;

   CResidue *res_p = new CResidue;
   CAtom *atom_p = new CAtom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" O  ");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" O");
   res_p->AddAtom(atom_p);
   res_p->seqNum = ires_prev + 1;
   res_p->SetResName("HOH");

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1; 
   make_bonds_type_checked();

}

// This is a bit messy, I'm afraid - we test single atom twice. If you use this for 
// a multiatom other than SO4 and P04, you will need to add it to the type test
// 
void 
molecule_class_info_t::add_typed_pointer_atom(coot::Cartesian pos, const std::string &type) { 

   short int single_atom = 1; // true

   std::cout << "INFO:: adding atom of type " << type << " at " << pos << std::endl;
   make_backup();

   // we get a chain pointer or NULL, if there is not yet a chain only
   // of the given type:
   CChain *single_type = coot::util::chain_only_of_type(atom_sel.mol, type);

   // We do different things (e.g adding the chain) if this is a new
   // chain or a pre-existing one, let's set a flag.
   short int pre_existing_chain_flag;
   CChain *chain_p;
   if (single_type) {
      chain_p = single_type;
      pre_existing_chain_flag = 1;
   } else {
      chain_p = new CChain;
      pre_existing_chain_flag = 0;
   }
   
   std::pair<short int, std::string> mol_chain_id = unused_chain_id();
   CResidue *res_p = new CResidue;

   // type test
   if (type == "PO4") single_atom = 0;
   if (type == "SO4") single_atom = 0;

   if (single_atom) { 
      CAtom *atom_p = new CAtom;
      float occ;
      if (is_from_shelx_ins_flag)
	 occ = 11.0;
      else
	 occ = 1.0;
      atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), occ, graphics_info_t::default_new_atoms_b_factor);

      if (type == "Water") { 

	 // special rule for water: we add a water to a water chain if
	 // possible

	 atom_p->SetAtomName(" O  ");
	 atom_p->SetElementName(" O");
	 res_p->SetResName("HOH");

	 CChain *w = water_chain();
	 int wresno = 1;
	 
	 if (w) { 
	    // add atom to chain w
	    std::pair<short int, int> wresno_pair = next_residue_in_chain(w);
	    if (wresno_pair.first) { 
	       wresno = wresno_pair.second;
	    } else { 
	       wresno = 1;
	    } 
	    res_p->seqNum = wresno;
	    res_p->AddAtom(atom_p);
	    w->AddResidue(res_p);
	    std::cout << atom_p << " added to molecule" << std::endl;
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	    
	 } else { 
	    // There was no water chain
	    res_p->AddAtom(atom_p);
	    std::cout << atom_p << " added to molecule" << std::endl;
	    if (!pre_existing_chain_flag) { 
	       chain_p->SetChainID(mol_chain_id.second.c_str());
	       atom_sel.mol->GetModel(1)->AddChain(chain_p);
	    }
	    res_p->seqNum = 1; // start of a new chain.
	    chain_p->AddResidue(res_p);
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	 } 
      } else { 
	
  	 // Not water

	 if (mol_chain_id.first || pre_existing_chain_flag) { 

	    if (type == "Br") { 
	       atom_p->SetAtomName("BR  ");
	       atom_p->SetElementName("Br");
	       res_p->SetResName("BR ");
	    } else { 
	      if (type == "Ca") { 
		atom_p->SetAtomName("CA  ");
		atom_p->SetElementName("Ca");
		res_p->SetResName("CA ");
	      } else { 
		if (type == "Na") { 
		  atom_p->SetAtomName("NA  ");
		  atom_p->SetElementName("Na");
		  res_p->SetResName("NA ");
		} else { 
		  if (type == "Cl") { 
		    atom_p->SetAtomName("CL  ");
		    atom_p->SetElementName("Cl");
		    res_p->SetResName("CL ");
		  } else { 
		    if (type == "Mg") { 
		      atom_p->SetAtomName("MG  ");
		      atom_p->SetElementName("Mg");
		      res_p->SetResName("MG ");
		    } else { 

		      // User Typed atom:

		      // make up (guess) the residue type and element
		      std::string at_name = type;
		      std::string ele = type;
		      std::string resname = type;
		      if (type.length() > 4)
			at_name = at_name.substr(0,4);
		      if (type.length() > 3)
			resname = at_name.substr(0,3);
		      if (type.length() > 2)
			ele = at_name.substr(0,2);

		      res_p->seqNum = 1; // start of a new chain.
		      atom_p->SetAtomName(at_name.c_str());
		      atom_p->SetElementName(ele.c_str());
		      res_p->SetResName(resname.c_str());
		    } 
		  }
		}
	      }
	    }

	    res_p->AddAtom(atom_p);
	    std::cout << atom_p << " added to molecule" << std::endl;
	    if (! pre_existing_chain_flag) { 
	       chain_p->SetChainID(mol_chain_id.second.c_str());
	       atom_sel.mol->GetModel(1)->AddChain(chain_p);
	    }
	    std::pair<short int, int> ires_prev_pair = coot::util::max_resno_in_chain(chain_p);
	    int previous_max = 0;
	    if (ires_prev_pair.first) { // was not an empty chain
	       previous_max =  ires_prev_pair.second;
	    }
	    chain_p->AddResidue(res_p);
	    res_p->seqNum = previous_max + 1;
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	 } else { 
	    std::cout << "WARNING:: Can't find new chain for new atom\n";
	 } 
      } // type was water, or not

   } else { 
      // multi atom:

      if (mol_chain_id.first || pre_existing_chain_flag) { 
	 add_pointer_multiatom(res_p, pos, type);
	 if (! pre_existing_chain_flag) { 
	    chain_p->SetChainID(mol_chain_id.second.c_str());
	    atom_sel.mol->GetModel(1)->AddChain(chain_p);
	 }
	 std::pair<short int, int> ires_prev_pair = coot::util::max_resno_in_chain(chain_p);
	 int previous_max = 0;
	 if (ires_prev_pair.first) { // was not an empty chain
	    previous_max =  ires_prev_pair.second;
	 } 
	 res_p->seqNum = previous_max + 1;
	 
	 chain_p->AddResidue(res_p);
	 atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	 atom_sel.mol->FinishStructEdit();
	 atom_sel = make_asc(atom_sel.mol);
	 have_unsaved_changes_flag = 1;
	 make_bonds_type_checked();
      } else { 
	 std::cout << "WARNING:: Can't find new chain for new atom\n";
      } 
   }
   // or we could just use update_molecule_after_additions() there.
}

// return status [1 means "usable"] and a chain id [status = 0 when
// there are 2*26 chains...]
// 
std::pair<short int, std::string>
molecule_class_info_t::unused_chain_id() const { 
   
   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
   std::pair<short int, std::string> s(0,""); 
   CChain *chain_p;
   if (atom_sel.n_selected_atoms > 0) { 
      CModel *model_p = atom_sel.mol->GetModel(1);
      int nchains = model_p->GetNumberOfChains();
      int idx;
   
      for (int ich=0; ich<nchains; ich++) {
	 chain_p = model_p->GetChain(ich);
	 idx = r.find(std::string(chain_p->GetChainID()));
	 // 	 while (idx != std::string::npos) { 
	    r = r.substr(0, idx) + r.substr(idx+1);
	    idx = r.find(std::string(chain_p->GetChainID()));
	    s.first = 1;
	    // 	 } // Take out the while, as per Ezra's suggestion.
      }
      std::string tstring = r.substr(0,1);
      s.second = tstring;
   } else {
      s.first = 1;
      s.second = "A";
   } 
   return s;
} 

void 
molecule_class_info_t::add_pointer_multiatom(CResidue *res_p, 
					     const coot::Cartesian &pos, const std::string &type) {
   
   coot::Cartesian p;
   float bf = graphics_info_t::default_new_atoms_b_factor;
   res_p->SetResName(type.c_str());
   if (type == "SO4") { 
      CAtom *atom_p;

      atom_p = new CAtom;
      p = pos + coot::Cartesian(0.000, 0.000, 0.088);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" S  ");
      atom_p->SetElementName(" S");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 1.227, 0.000, -0.813);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O1 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian(-1.227, 0.000, -0.813);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O2 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian(  0.000, -1.263, 0.740);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O3 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 0.000, 1.263, 0.740);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O4 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
   } 

   if (type == "PO4") { 
      CAtom *atom_p;

      atom_p = new CAtom;
      p = pos + coot::Cartesian(0.000, 0.021, 0.036);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" P  ");
      atom_p->SetElementName(" P");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 1.315,   0.599,  -0.691 );
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O1 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( -1.315,   0.599,  -0.691);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O2 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 0.000,  -1.587,  -0.055);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O3 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 0.000,   0.434,   1.457);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O4 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
   } 

} 


// ----------------------------------------------------------------------
//                 Save yourself
// ----------------------------------------------------------------------
//
// return 0 on success.
int
molecule_class_info_t::save_coordinates(const std::string filename) {

   int ierr = 0;
   std::string ext = coot::util::file_name_extension(filename);
   if (coot::util::extension_is_for_shelx_coords(ext)) {
      write_shelx_ins_file(filename);
   } else {
      ierr = write_atom_selection_file(atom_sel, filename);
   }

   if (ierr) {
      std::cout << "WARNING!! Coordinates write to " << filename
		<< " failed!" << std::endl;
      std::string ws = "WARNING:: export coords: There was an error ";
      ws += "in writing ";
      ws += filename;
      GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(ws);
      gtk_widget_show(w);
   } else {
      have_unsaved_changes_flag = 0;

      // Now we have updated the molecule name, how shall we restore
      // this from the state file?
      std::vector<std::string> strings;
      strings.push_back("handle-read-draw-molecule");
      strings.push_back(single_quote(coot::util::intelligent_debackslash(filename)));
      save_state_command_strings_ = strings;

      name_ = filename;  // hmm... // update go to atom widget now? FIXME.
      std::string::size_type icoot = filename.rfind("-coot-");
      if (icoot != std::string::npos) { 
	 coot_save_index++;
      }
      update_mol_in_display_control_widget();  // FIXME. 
   }
   return ierr;
} 



// Return 1 on yes, unsaved changes present,
//        0 on no
int
molecule_class_info_t::Have_unsaved_changes_p() const {
   if (has_model())
      return have_unsaved_changes_flag;
   else
      return 0;
}


// ----------------------------------------------------------------------
//               Baton Atoms
// ----------------------------------------------------------------------


// Recall that the chain is set by the creation of the empty molecule in 
// graphics_info_t::baton_build_atoms_molecule();
// direction_flag is +1 for forward building, -1 for backwards direction.
//
// In the case of direction_flag being negative, I think that residues
// (atoms) will go into the chain with their seqNums in decreasing
// order.  The may need to be sorted later, I think (maybe not).
// 
CAtom *
molecule_class_info_t::add_baton_atom(coot::Cartesian pos, 
				      int istart_resno,
				      short int iresno_active,
				      short int direction_flag) {

   int nchains = atom_sel.mol->GetNumberOfChains(1);

   if (nchains != 1) {
      std::cout << "failed to add baton atom" << std::endl;
      return NULL;
   }
      
   make_backup();
   CChain *chain_p = atom_sel.mol->GetChain(1,0);
   
   std::string mol_chain_id(chain_p->GetChainID());
   int n_res = chain_p->GetNumberOfResidues();


   // if this is the first atom to be added in the chain, we get the
   // seqnum from the passed istart_resno otherwise it is the seqnum
   // of the previous residue, plus or minus one.
   // 
   int this_res_seqnum; 
   if (n_res == 0) {
      this_res_seqnum = istart_resno;
   } else { 

      if (iresno_active == 0) { 
	 int ires_prev = chain_p->GetResidue(n_res-1)->seqNum; // seqnum of the last
                                                               // residue in chain.
	 this_res_seqnum = ires_prev + 1*direction_flag;
      } else { 
	 this_res_seqnum = istart_resno;
      } 
   }

   CResidue *res_p = new CResidue;
   CAtom *atom_p = new CAtom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" CA ");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0,
			  graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" C");
   res_p->AddAtom(atom_p);
   res_p->seqNum = this_res_seqnum;
   res_p->SetResName("ALA");

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1; 
   make_ca_bonds(2.4, 4.7);

   return atom_p;
}

// Return a vector of upto 3 positions of the most latestly added
// atoms with the most lastest atom addition (that is the passed atom)
// in the back() slot of the vector.
// 
std::vector<clipper::Coord_orth>
molecule_class_info_t::previous_baton_atom(const CAtom* latest_atom_addition,
					   short int direction) const {
   std::vector<clipper::Coord_orth> positions;
   int direction_sign = +1; 

   if (direction == 1) { // building forward, look in negative
			 // direction for previously build atoms.
      direction_sign = +1;
   } else { 
      direction_sign = -1; // building backward, look in positive
			   // direction for previously build atoms.
   }
   int ires_last_atom = ((CAtom *) latest_atom_addition)->GetSeqNum();

   char *chain = ((CAtom *) latest_atom_addition)->GetChainID();
   // does the CA for the (ires_last_atom-2) exist?
   int selHnd = atom_sel.mol->NewSelection();
	 
   atom_sel.mol->SelectAtoms(selHnd, 0, chain,
			     ires_last_atom-2*direction_sign, "*", 
			     ires_last_atom-2*direction_sign, "*",
			     "*", " CA ", "*", "*"); 
      
   int nSelAtoms;
   PPCAtom local_SelAtom; 
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);
      
   if (nSelAtoms == 0) {
      std::cout << "residue with sequence number " << ires_last_atom - 2*direction_sign
		<< " not found for ires_last_atom = " << ires_last_atom 
		<< " with direction_sign = " << direction_sign << "\n";
   } else {
      positions.push_back(clipper::Coord_orth(local_SelAtom[0]->x,
					      local_SelAtom[0]->y,
					      local_SelAtom[0]->z));
   }
   atom_sel.mol->DeleteSelection(selHnd);

   // rinse, lather, repeat...

   // does the CA for the (ires_last_atom-1) exist?
   selHnd = atom_sel.mol->NewSelection();
      
   atom_sel.mol->SelectAtoms(selHnd, 0, chain,
			     ires_last_atom-direction_sign, "*", 
			     ires_last_atom-direction_sign, "*",
			     "*", " CA ", "*", "*"); 

   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);
      
   if (nSelAtoms == 0) {
      std::cout << "residue with sequence number " << ires_last_atom - direction_sign
		<< " not found\n";
   } else {
      positions.push_back(clipper::Coord_orth(local_SelAtom[0]->x,
					      local_SelAtom[0]->y,
					      local_SelAtom[0]->z));
   }
   atom_sel.mol->DeleteSelection(selHnd);

   // And finally this one is guaranteed to exist:
   // 
   positions.push_back(clipper::Coord_orth(latest_atom_addition->x,
					   latest_atom_addition->y,
					   latest_atom_addition->z));

   return positions;
   
} 

#include "CalphaBuild.hh"

std::vector<coot::scored_skel_coord>
molecule_class_info_t::next_ca_by_skel(const std::vector<clipper::Coord_orth> &previous_ca_positions,
				       const clipper::Coord_grid &coord_grid_start,
				       short int use_coord_grid_start_flag,
				       float ca_ca_bond_length,
				       float map_cut_off,
				       int max_skeleton_search_depth) const {

   std::vector<coot::scored_skel_coord> t; 
   coot::CalphaBuild buildca(max_skeleton_search_depth);

   if (skeleton_treenodemap_is_filled) { 
      t = buildca.next_ca_by_skel(previous_ca_positions,
				  coord_grid_start,
				  use_coord_grid_start_flag,
				  ca_ca_bond_length,
				  xskel_cowtan, xmap_list[0],
				  map_cut_off,
				  skeleton_treenodemap);
   } else {
      std::cout << "treenodemap is not filled" << std::endl;
   }
   return t;
} 

#include <time.h>

// ----------------------------------------------------------------------
//               Dummy Atoms (not bonded)
// ----------------------------------------------------------------------
void
molecule_class_info_t::add_dummy_atom(coot::Cartesian pos) {

   int nchains = atom_sel.mol->GetNumberOfChains(1);

   if (nchains != 1) {
      std::cout << "failed to add dummy atom" << std::endl;
      return;
   }

   make_backup();

   CChain *chain_p = atom_sel.mol->GetChain(1,0);
   
   std::string mol_chain_id(chain_p->GetChainID());
   int ires_prev = chain_p->GetNumberOfResidues();

   CResidue *res_p = new CResidue;
   CAtom *atom_p = new CAtom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" DUM");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0,
			  graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" O");
   res_p->AddAtom(atom_p);
   res_p->seqNum = ires_prev + 1;
   res_p->SetResName("DUM");

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   // std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1; 
   makebonds(0.0, 0.0);

}


void
molecule_class_info_t::fill_skeleton_treenodemap() {

   // if we have a skeleton map but not treenodemap:
   // 
   if (xskel_is_filled && !skeleton_treenodemap_is_filled) {

      // Chomp up the lovely memory! Yum!
      // 
      skeleton_treenodemap.init(xskel_cowtan.spacegroup(), 
				xskel_cowtan.cell(),
				xskel_cowtan.grid_sampling()); 
      clipper::Coord_grid c_g; 
      clipper::Skeleton_basic::Neighbours skel_neighbs(xskel_cowtan);
      
//       std::cout << "Build tree: there are " << skel_neighbs.size() << " skel_neighbs"
// 		<< std::endl;  18, actually.
      
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xskel_cowtan.first(); !ix.last(); ix.next() ) {
	 if (xskel_cowtan[ix] > 0) { 

	    coot::SkeletonTreeNode stn; 
	    
	    for(int i=0; i< skel_neighbs.size(); i++) {
	       c_g = ix.coord() + skel_neighbs[i]; 
	    
	       if (xskel_cowtan.get_data(c_g) > 0 ) {
	       
		  // OK, so this node has a neighbour:
		  // 
		  stn.neighbs.push_back(c_g);
	       }
	    }
	    stn.near_grid_point = ix.coord();  // Strange but true!
	    // 
	    // We do this because "out of cell" reference
	    // (e.g.  uvw = (  -1, -12, -19)) will get wrapped 
	    // to some (hidden) value.  To get the wrapped
	    // value (i.e the grid), we look it up here. 
	    // Cunning (if it works). 
	    skeleton_treenodemap[ix] = stn; 
	 }
      }
      // set the flag
      skeleton_treenodemap_is_filled = 1;
   }
}


float
molecule_class_info_t::density_at_point(const clipper::Coord_orth &co) const {

   if (!xmap_is_filled[0]) {
      std::cout << " returning bogus value from density_at_point: " << std::endl;
      return -1000.0;
   } else {
      clipper::Coord_frac cf = co.coord_frac(xmap_list[0].cell());
      clipper::Coord_grid cg = cf.coord_grid(xmap_list[0].grid_sampling());
      return xmap_list[0].get_data(cg);
   }
}


// backups:

// Backup filename: return a stub.
// 
std::string
molecule_class_info_t::save_molecule_filename(const std::string &dir) { 

   std::string time_string = save_time_string;
   graphics_info_t g;
   
   if ((history_index == 0) ||
       history_index != max_history_index) { 

      time_string = dir;

      // unix dependent logic here:  Don't know how to do this on other systems...
      // We want a filename proceeded by a directory name:
      // i.e. we end up with something like
      // "coot-backup/a.pdb_Tues_Aug_19_20:16:00_2003_modification_0.mmdbbin"

      time_string += "/";

      std::string clean_name = name_;
      if (g.unpathed_backup_file_names_flag) {
	 clean_name = name_for_display_manager();
      }
      // convert "/" to "_"
      int slen = clean_name.length();
      for (int i=0; i<slen; i++)
	 if (clean_name[i] == '/')
	    clean_name[i] = '_';

      time_string += clean_name;
      time_string += "_";
      
      // add in the time component:

#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)

      // but not if we are in windows:

#else      
      time_t t;
      time(&t);
      char *chars_time = ctime(&t);
      time_string += chars_time;
#endif

      // strip off the trailing newline:
      slen = time_string.length();
      if (slen > 2) 
	 time_string = time_string.substr(0,slen-1);
      
      // convert spaces to underscores
      // 
      for (unsigned int i=0; i<time_string.length(); i++)
	 if (time_string[i] == ' ')
	    time_string[i] = '_';
	    
#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)

      // convert : to underscores in windows
      // 
      for (int i=0; i<time_string.length(); i++)
	 if (time_string[i] == ':')
	    time_string[i] = '_';
#endif
      
      time_string += "_modification_";

      save_time_string = time_string; // why do we do this?  Ah, because we want the
                                      // time to calculated at the start:
                                      // and use that as a stub.

      time_string += g.int_to_string(history_index);
      //time_string += ".mmdbbin";
      if (! is_from_shelx_ins_flag) 
	 time_string += ".pdb";
      else 
	 time_string += ".res";

#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)

#else
      if (! is_from_shelx_ins_flag)
	 time_string += ".gz"; // 'cos we can do compression.  Groovy baby!
#endif
      
   } else {
      // (this is not the first save molecule that we have done)
      
      // add to the stub that we have previously generated.
      //
      time_string += g.int_to_string(history_index);
      if (! is_from_shelx_ins_flag) 
	 time_string += ".pdb";
      else 
	 time_string += ".res";
#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)
#else
      if (! is_from_shelx_ins_flag) 
	 time_string += ".gz";
#endif
   } 
   return time_string;
}

// Return like mkdir: mkdir returns zero on success, or -1 if an  error  occurred
//
// if it already exists as a dir, return 0 of course.
// 
int
molecule_class_info_t::make_maybe_backup_dir(const std::string &backup_dir) const {

   return coot::util::create_directory(backup_dir);
}

// Ignore return value.
//
// If successful, increase history_index and if not in a backup
// increase max_history_index too.
// 
int
molecule_class_info_t::make_backup() { // changes history details

   if (backup_this_molecule) { 
      std::string backup_dir("coot-backup");

      //shall we use the environment variable instead?
      char *env_var = getenv("COOT_BACKUP_DIR");
      if (env_var) { 
	 struct stat buf;
	 int err = stat(env_var, &buf);
	 if (!err) {
	    if (! S_ISDIR(buf.st_mode)) {
	       env_var = NULL;
	    }
	 } else {
	    env_var = NULL;
	 }
      }
      if (env_var)
	 backup_dir = env_var;

      if (atom_sel.mol) {
	 int dirstat = make_maybe_backup_dir(backup_dir);

	 // all is hunkey-dorey.  Directory exists.
	 if (dirstat == 0) { 
	    
	    std::string backup_file_name = save_molecule_filename(backup_dir);
 	    std::cout << "INFO:: backup file " << backup_file_name << std::endl;

#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)
	    byte gz = GZM_NONE;
#else
	    byte gz = GZM_ENFORCE;
#endif
	    // Writing out a modified binary mmdb like this results in the
	    // file being unreadable (crash in mmdb read).
	    // 
	    // int istat = atom_sel.mol->WriteMMDBF((char *)backup_file_name.c_str(), gz);
	    int istat;
	    if (! is_from_shelx_ins_flag) {
	       istat = atom_sel.mol->WritePDBASCII((char *)backup_file_name.c_str(), gz);
	       // WriteMMDBF returns 0 on success, else mmdb:Error_CantOpenFile (15)
	       if (istat) { 
		  std::cout<< "WARNING:: WritePDBASCII failed! Return status " << istat << std::endl;
	       }
	    } else { 
	       std::pair<int, std::string> p = write_shelx_ins_file(backup_file_name);
	       istat = p.first;
	    }
		  
	    save_history_file_name(backup_file_name);
	    if (history_index == max_history_index)
	       max_history_index++;
	    history_index++;
	 } else {
	    std::cout << "BACKUP:: directory "<< backup_dir << " failure" << std::endl;
	 } 
      } else {
	 std::cout << "BACKUP:: Ooops - no atoms to backup for this empty molecule"
		   << std::endl;
      }
   } else {
      // Occasionally useful but mostly tedious...
      // std::cout << "INFO:: backups turned off on this molecule"
      // << std::endl;
   }
   return 0;
}


void
molecule_class_info_t::save_history_file_name(const std::string &file) {

   // First, history_index is zero and the vec is zero,
   // normal service, then another backup: history_index is 1 and vec is 1.
   // 
   if (history_index == int(history_filename_vec.size())) {
      history_filename_vec.push_back(file);
   } else {
      // we have gone back in history.
      // 
      if (history_index < int(history_filename_vec.size())) {
	 history_filename_vec[history_index] = file;
      }
   }
} 

// restore from (previous) backup
void
molecule_class_info_t::restore_from_backup(int history_offset) {

   int hist_vec_index = history_index + history_offset;
   if (int(history_filename_vec.size()) > hist_vec_index) {
      std::cout << "restoring from backup " << history_filename_vec.size()
		<< " " << history_index << std::endl;
      std::string save_name = name_;
      std::string filename = history_filename_vec[hist_vec_index];
      //      history_index = hist_index;
      short int reset_rotation_centre = 0;
      // handle_read_draw_molecule uses graphics_info_t::n_molecules
      // to determine its molecule number.  We don't want it to
      // change.
      int save_imol = imol_no;
      // similarly, it messes with the save_state_command_strings_, we
      // don't want that either:
      std::vector<std::string> save_save_state = save_state_command_strings_;
      short int is_undo_or_redo = 1;
      handle_read_draw_molecule(filename, reset_rotation_centre, is_undo_or_redo);
      save_state_command_strings_ = save_save_state;
      imol_no = save_imol; 
      name_ = save_name;
   } else {
      std::cout << "not restoring from backup because "
		<< history_filename_vec.size()
		<< " " << history_index << std::endl;
   }
}


// I need to write an essay on how that backup system works.
//
// Insight: if we are at hist_index = max_hist_index (i.e. not in a
// backup) then when an undo is requested, we should make a backup.
// This makes the indexing a bit tricky,
//
// So imagine this situation:
//
//             hist_index max_hist_index                filenames filled
// pepflip         1         1                              [0]
// rotate          2         2                              [0,1]
// undo            3         3 [first step is a backup]     [0,1,2]
//
// filenames get routinely pushed back onto history_filename_vec
// (i.e. not in an undo situation)
//
// So imagine we make one mod, then the history_filename_vec size() is 1.
// on undo:
//    we restore from backup using history_filename_vec index 0.
//    we have added to history_filename_vec [now has size 2] in this proceedure
//    history_index was 1 on starting the undo
//    at end of undo it is 0.
//
// how do we redo that?
//    (obviously) we restore from backup using history_filename_vec index 1.
//    we have not added to history_filename_vec [now has size 2]
//    history_index was 0 on starting the redo.
//    at end of redo, history_index is 1
//
// So having done 2 mods:
//
// on undo:
//    restore from backup using history_filename_vec index 1
//    we have added to history_filename_vec [now has size 3] in this proceedure
//    history_index was 2 on starting undo
//    at end of undo it is 1.
//
// on redo:
//    restore from backup using history_filename_vec index 2
//    we have not added to history_filename_vec [now has size 3]
//    history_index was 1 on starting the redo
//    at end of redo, history_index was 2


//
// [It would be cool to have the Redo button greyed out when there are
// no redos availabile (set its state when either it or undo is
// pressed)]
// 
// initially it should be greyed out (insensitive).

// restore from (next) backup
void
molecule_class_info_t::apply_undo() {

//    std::cout << std::endl << "DEBUG:: in apply undo start hist_index: "
// 	     << history_index
// 	     << " max_history_index: " << max_history_index << std::endl;

   if (history_index > 0) {
      int offset = -1;
      if (history_index == max_history_index) { 
	 make_backup(); // increments history_index
	 offset--;
      }
      restore_from_backup(offset);
      history_index += offset;

      // So that we don't get asked to save the molecule on exist when
      // we have reverted all our modifications:
      // 
      if (history_index == 0) { 
	 have_unsaved_changes_flag = 0;
      }
   }

   std::cout << "DEBUG:: apply_undo: (end) history_index: " <<
      history_index << " max_history_index: " << max_history_index << std::endl;

}

void
molecule_class_info_t::apply_redo() {

   if (history_index < max_history_index) {
      std::cout << "DEBUG:: molecule applying redo " << history_index << std::endl;

      // When there are 3 backups made and we are viewing molecule 2,
      // we don't want to restore from history_filename_vec[3]:
      // 
      if (int(history_filename_vec.size()) > (history_index + 1)) { 
	 restore_from_backup(+1); 
	 history_index++; 
	 have_unsaved_changes_flag = 1;
      } else {
	 std::cout << "Not redoing history file vec: " << history_filename_vec.size()
		   << " " << history_index << std::endl;
      } 
   } else {
      std::cout << "Not redoing history: " << max_history_index
		<< " " << history_index << std::endl;
   }
} 



// For model view (go to atom)
//
std::vector<coot::model_view_residue_button_info_t>
molecule_class_info_t::model_view_residue_button_labels() const {

   std::vector<coot::model_view_residue_button_info_t> v;

   if (atom_sel.n_selected_atoms > 0) { 

      int nchains = atom_sel.mol->GetNumberOfChains(1);

      if (nchains < 1) {
	 std::cout << "failed to find chains for atom in "
		   << " model_view_residue_button_info_t" << std::endl;
      } else {

	 graphics_info_t g;
	 CModel *model_p = atom_sel.mol->GetModel(1);

	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in model_view_residue_button_info_t: "
		      << nchains << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in model_view_residue_button_info_t: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) { 
		     PCResidue residue_p = chain->GetResidue(ires);
		     std::string button_label = 
			g.int_to_string(residue_p->GetSeqNum());
		     button_label += " ";
		     button_label += residue_p->GetChainID();
		     button_label += " ";
		     button_label += residue_p->name;

		     v.push_back(coot::model_view_residue_button_info_t(button_label,
									residue_p));
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

// return vector of atom list (aka button) info for this residue
// 
std::vector<coot::model_view_atom_button_info_t>
molecule_class_info_t::model_view_atom_button_labels(char *chain_id, int seqno) const {

   graphics_info_t g;
   std::vector<coot::model_view_atom_button_info_t> v;

   // protection against the molecule having been deleted after the
   // gtklist widget was created:
   if (atom_sel.n_selected_atoms > 0) { 

      CChain *chain;
      
      // first we have to find the residue res_p (from which we wil get the atoms)
      //
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain = atom_sel.mol->GetChain(1, ichain);
	 if (chain == NULL) { 
	    // This should not be necessary. It seem to be a result of
	    // mmdb corruption elsewhere - possibly DeleteChain in
	    // update_molecule_to().
	    std::cout << "ERROR getting chain in model_view_atom_button_info_t\n";
	 } else { 
	    std::string residue_chain_id(chain_id); // passed from residue list
	    std::string mol_chain_id(chain->GetChainID());
	    if (residue_chain_id == mol_chain_id) {
	       int nres = chain->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) { 
		  PCResidue res_p = chain->GetResidue(ires);
		  if (res_p->GetSeqNum() == seqno) {
      
		     PPCAtom residue_atoms;
		     int nResidueAtoms;

		     res_p->GetAtomTable(residue_atoms, nResidueAtoms);
		     for (int i=0; i<nResidueAtoms; i++) {
			std::string button_label = residue_atoms[i]->name;
			std::string altConf = residue_atoms[i]->altLoc;
			if (altConf != "") { 
			   button_label += ",";
			   button_label += altConf;
			}
			button_label += " occ=";
			button_label += g.float_to_string(residue_atoms[i]->occupancy);
			button_label += " bf=";
			button_label += g.float_to_string(residue_atoms[i]->tempFactor);
	    
			v.push_back(coot::model_view_atom_button_info_t(button_label, residue_atoms[i]));
		  
		     }
		  }
	       }
	    } 
	 } 
      }
   }
   return v;
}


std::vector<coot::model_view_atom_tree_chain_t>
molecule_class_info_t::model_view_residue_tree_labels() const {

   std::vector<coot::model_view_atom_tree_chain_t> v;

   if (atom_sel.n_selected_atoms > 0) {

      CChain *chain_p;
      int im = 1;
      int nchains = atom_sel.mol->GetNumberOfChains(im);
      for (int ichain=0; ichain<nchains; ichain++) {

	 chain_p = atom_sel.mol->GetChain(im, ichain);
	 std::string chain_label("Chain ");
	 chain_label += chain_p->GetChainID();
	 v.push_back(coot::model_view_atom_tree_chain_t(chain_label));
	 
	 if (! chain_p) {
	    std::cout << "ERROR getting chain in model_view_residue_tree_labels\n";
	 } else {
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       std::string label = residue_p->GetChainID();
	       label += " ";
	       label += coot::util::int_to_string(residue_p->GetSeqNum());
	       label += residue_p->GetInsCode();
	       label += " ";
	       label += residue_p->name;
	       coot::model_view_atom_tree_item_info_t res(label, residue_p);
	       v.back().add_residue(res);
	    }
	 }
      }
   }
   return v;
}



// return success on residue type match
// success: 1, failure: 0.
int 
molecule_class_info_t::mutate_single_multipart(int ires_serial, const char *chain_id, const std::string &target_res_type) {

   int istat = 0;
   if (atom_sel.n_selected_atoms > 0) {
      CChain *chain_p;
      int nres;
      int nchains = atom_sel.mol->GetNumberOfChains(1) ;
      for (int ichain =0; ichain<nchains; ichain++) {
	 chain_p = atom_sel.mol->GetChain(1,ichain);
	 if (std::string(chain_id) == std::string(chain_p->GetChainID())) { 
	    nres = chain_p->GetNumberOfResidues();
	    if (ires_serial < nres) {
	       CResidue *res_p = chain_p->GetResidue(ires_serial);
	       if (res_p) {

		  if (std::string(res_p->name) == target_res_type) {

		     std::cout << "residue type match for ires = " << ires_serial << std::endl;
		     istat = 1; // success
		  
		  } else {
		  
		     // OK, do the mutation:
	       
		     // get an instance of a standard residue of type target_res_type
		     CResidue *std_res = get_standard_residue_instance(target_res_type); // a deep copy
		     // move the standard res to position of res_p
		     // move_std_residue(moving_residue, (const) reference_residue);
		     if (std_res) { 
			istat = move_std_residue(std_res, (const CResidue *)res_p);
			
			if (istat) { 
			   mutate_internal(res_p, (const CResidue *)std_res);
			} else { 
			   std::cout << "WARNING:  Not mutating residue due to missing atoms!\n";
			} 
			// atom_selection and bonds regenerated in mutate_internal
		     } else {
			std::cout << "ERROR failed to get residue of type :" << target_res_type
				  << ":" << std::endl;
		     }
		  }
	       } else {
		  std::cout << "ERROR:: in mutate_single_multipart oops - can't get residue"
			    << " with ires_serial: " << ires_serial << std::endl;
	       } 
	    } else {
	       std::cout << "PROGRAMMER ERROR: out of range residue indexing" << std::endl;
	    } 
	 }
      }
   }
   return 0 + istat;
} 

// Return 0 on failure.
short int 
molecule_class_info_t::move_std_residue(CResidue *moving_residue,
					const CResidue *reference_residue) const {

   std::pair<clipper::RTop_orth, short int> pair = 
      coot::util::get_ori_to_this_res((CResidue *)reference_residue); 

   short int istat = 0; // success

   if (!reference_residue) { 
      std::cout << "This should not happen!" << std::endl;
      std::cout << "null reference residue in move_std_residue" << std::endl;
   } else { 

      if (pair.second == 1) { // successful attempt to get the matrix
	 PPCAtom residue_atoms = NULL;
	 int nResidueAtoms;
	 moving_residue->GetAtomTable(residue_atoms, nResidueAtoms);
	 if (nResidueAtoms == 0) {
	    std::cout << " something broken in atom residue selection in ";
	    std::cout << "mutate, got 0 atoms" << std::endl;
	    istat = 0;
	 } else {
	    istat = 1;
// 	    std::cout << "DEBUG:: move_std_residue: " << nResidueAtoms
// 		      << " atoms in residue " 
// 		      << moving_residue << " " << moving_residue->seqNum << " " 
// 		      << moving_residue->GetChainID() << std::endl;
	    for(int iat=0; iat<nResidueAtoms; iat++) {
	       if (residue_atoms[iat]) { 
// 		  std::cout  << "residue atom " << iat << " coords: " 
// 			     << residue_atoms[iat]->x << " "
// 			     << residue_atoms[iat]->y << " "
// 			     << residue_atoms[iat]->z << std::endl;
		  clipper::Coord_orth co(residue_atoms[iat]->x,
					 residue_atoms[iat]->y,
					 residue_atoms[iat]->z);
		  clipper::Coord_orth rotted = co.transform(pair.first); // an rtop
		  residue_atoms[iat]->x = rotted.x();
		  residue_atoms[iat]->y = rotted.y();
		  residue_atoms[iat]->z = rotted.z();
	       } else { 
		  istat = 0;
		  std::cout << "ERROR:: bad residue atom in move_std_residue: iat: "
			    << iat << std::endl;
	       }
	    }
	 }
      } else { 
	 istat = 0; // failure
	 std::cout << "DISASTER - failed to generate RTop for move_std_residue\n";
	 if (reference_residue) { 
	    // 	 molecule-class-info.cc:4184: passing `const CResidue' as `this' 
	    // argument of `int CResidue::GetSeqNum ()' discards qualifiers
	    CResidue *tmp = (CResidue *) reference_residue;
	    std::cout << "mainchain atoms missing from residue " 
		      << tmp->GetSeqNum() 
		      << tmp->GetChainID() << std::endl;
	 } else { 
	    std::cout << "This should not happen!" << std::endl;
	    std::cout << "null residue in move_std_residue" << std::endl;
	 }
      }
   }
   return istat;
}


void
molecule_class_info_t::make_backup_from_outside() {  // when we have a multi mutate, we
				    // want the wrapper to make a
				    // backup when we start and set
				    // changes when when finish.
				    // Rather crap that this needs to
				    // be done externally, I think.

   std::cout << " calling make_backup from make_backup_from_outside" << std::endl;
   make_backup();
}


void
molecule_class_info_t::set_have_unsaved_changes_from_outside() {

   have_unsaved_changes_flag = 1;

}


//
// Get a deep copy:
// return NULL on failure
// 
CResidue *
molecule_class_info_t::get_standard_residue_instance(const std::string &residue_type) {

   graphics_info_t g;
   CResidue *std_residue = NULL;
   
   if (g.standard_residues_asc.read_success) { 
//      std::cout << "DEBUG:: There are " << g.standard_residues_asc.n_selected_atoms
// 	       << " atoms in standard_residues_asc" << std::endl;
     int selHnd = g.standard_residues_asc.mol->NewSelection();
     g.standard_residues_asc.mol->Select ( selHnd,STYPE_RESIDUE, 1, // .. TYPE, iModel
					   "*", // Chain(s) it's "A" in this case.
					   ANY_RES,"*",  // starting res
					   ANY_RES,"*",  // ending res
					   (char *) residue_type.c_str(),  // residue name
					   "*",  // Residue must contain this atom name?
					   "*",  // Residue must contain this Element?
					   "*",  // altLocs
					   SKEY_NEW // selection key
					   );
     // get the standard orientation residue for this residue type
     PPCResidue SelResidue;
     int nSelResidues;

     g.standard_residues_asc.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   
     if (nSelResidues != 1) {
       std::cout << "This should never happen - ";
       std::cout << "badness in get_standard_residue_instance, we selected " << nSelResidues
		 << " residues looking for residues of type :" << residue_type << ":\n";
     } else {
       std_residue = coot::deep_copy_this_residue(SelResidue[0], "", 1, 
						  g.standard_residues_asc.UDDAtomIndexHandle);
     }
     g.standard_residues_asc.mol->DeleteSelection(selHnd);
   }
   return std_residue;
}

// 1: success
// 0: failure
// 
short int
molecule_class_info_t::progressive_residues_in_chain_check_by_chain(const char *chain_id) const {

   short int r = 0;
   
   if (atom_sel.n_selected_atoms > 0) { 
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 const CChain *chain_p = (const CChain *)atom_sel.mol->GetChain(1,ichain);
	 std::string mol_chain_id(((CChain*)chain_p)->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) { 
	    r = coot::progressive_residues_in_chain_check(chain_p);
	    break;
	 }
      }
   }

   return r; 
} 


// return the number of residues in chain with chain_id, return -1 on error
// 
int
molecule_class_info_t::chain_n_residues(const char *chain_id) const {

   int r = -1;
   
   if (atom_sel.n_selected_atoms > 0) {
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 const CChain *chain_p = (const CChain *)atom_sel.mol->GetChain(1,ichain);
	 std::string mol_chain_id(((CChain*)chain_p)->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    r = ((CChain *)chain_p)->GetNumberOfResidues();
	 }
      }
   }
   return r; 
} 


void
molecule_class_info_t::store_refmac_params(const std::string &mtz_filename,
					   const std::string &fobs_col,
					   const std::string &sigfobs_col,
					   const std::string &r_free_col, 
					   int r_free_flag) { 

   have_sensible_refmac_params = 1; // true
   refmac_mtz_filename = mtz_filename; 
   refmac_fobs_col = fobs_col;
   refmac_sigfobs_col = sigfobs_col;
   refmac_r_free_col = r_free_col;
   refmac_r_free_flag_sensible = r_free_flag;

   std::cout << "INFO:: Stored refmac parameters: " 
	     << refmac_fobs_col << " "
	     << refmac_sigfobs_col;
   if (r_free_flag)
      std::cout << " " << refmac_r_free_col << " is sensible." << std::endl;
   else
      std::cout << " the r-free-flag is not sensible" << std::endl;
} 


// return 0 on success
// 
int
molecule_class_info_t::write_pdb_file(const std::string &filename) {

   int err = 1; // fail
   if (atom_sel.n_selected_atoms > 0) { 
      std::string ext = coot::util::file_name_extension(filename);
      if (coot::util::extension_is_for_shelx_coords(ext)) {
	 write_shelx_ins_file(filename);
      } else {
	 err = write_atom_selection_file(atom_sel, filename);
      }
   }
   return err; 
} 


// Add this molecule (typically of waters to this molecule by trying
// to put them into an already-existing solvent chain).  If a solvent
// chain does not already exist, put create a new chain id for the
// water_mol atoms.
//
// All the atoms of water_mol need to be in a
// chain that has a different chain id to all the chains in this
// molecule.  Else fail (return status 0).
// 
int
molecule_class_info_t::insert_waters_into_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0;  // set to failure initially

   // So run over the the chains of the existing molecule looking for
   // a solvent chain.  If there isn't one we simply use
   // append_to_molecule()
   //
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   CChain *chain_p = NULL;
   CChain *solvent_chain_p = NULL; 
   short int i_have_solvent_chain_flag = 0;
   for (int ichain=0; ichain<nchains; ichain++) { 
      
      chain_p = atom_sel.mol->GetChain(1,ichain);
      if (chain_p->isSolventChain()) {
	 solvent_chain_p = chain_p;
	 std::string mol_chain_id(chain_p->GetChainID());
	 i_have_solvent_chain_flag = 1;
      }
   }


   // For every atom in water_mol, create a new atom and a new residue
   // for it. Add the residue to our model's solvent chain and the
   // atom the the residue (of course).
   //
   if (i_have_solvent_chain_flag == 0) {
      
      // We didn't manage to find a solvent chain.
      // We need to create a new chain.
      chain_p = new CChain;
      atom_sel.mol->GetModel(1)->AddChain(chain_p);
      std::pair<short int, std::string> u = unused_chain_id();
      if (u.first)
	 chain_p->SetChainID(u.second.c_str());
      else 
	 chain_p->SetChainID("Z");
   } else {
      chain_p = solvent_chain_p; // put it back, (kludgey, should use
				 // solvent_chain_p from here, not chain_p).
   } 

//    std::cout << "Debug:: choose chain " << chain_p->GetChainID()
// 	     << " with have_solvent flag: " << i_have_solvent_chain_flag
// 	     << std::endl;
//    std::cout << "Debug:: isSolvent for each residue of chain: " << std::endl;
//    for (int tmp_r=0; tmp_r<chain_p->GetNumberOfResidues(); tmp_r++) {
//       CResidue *rtmp = chain_p->GetResidue(tmp_r);
//       short int flag = isSolvent(rtmp->name);
// 	 std::cout << rtmp->name << " is solvent? " << flag << std::endl;
//    }
   
   std::pair<short int, int> p = coot::util::max_resno_in_chain(chain_p);
   float bf = graphics_info_t::default_new_atoms_b_factor; // 20.0 by default
   int max_resno;
   if (p.first) { 
      max_resno = p.second;
   } else {
      max_resno = 0;
   }
   if (p.first || (i_have_solvent_chain_flag == 0)) {
      make_backup();
      std::cout << "INFO:: Adding to solvent chain: " << chain_p->GetChainID()
		<< std::endl;
      int prev_max_resno = max_resno;
      CResidue *new_residue_p = NULL;
      CAtom    *new_atom_p = NULL;
      int water_count = 0;
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {
	 for (int ires=water_mol[ifrag].min_res_no();
	      ires<=water_mol[ifrag].max_residue_number();
	      ires++) {
	    for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {
	       new_residue_p = new CResidue;
	       new_residue_p->SetResName("HOH");
	       new_residue_p->seqNum = prev_max_resno + 1 + water_count; 
	       water_count++; 
	       new_atom_p = new CAtom;
	       new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
					  water_mol[ifrag][ires][iatom].pos.y(),
					  water_mol[ifrag][ires][iatom].pos.z(), 1.0, bf);
	       new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
	       strncpy(new_atom_p->element, water_mol[ifrag][ires][iatom].element.c_str(), 3);
	       strncpy(new_atom_p->altLoc, water_mol[ifrag][ires][iatom].altLoc.c_str(), 2);

	       // residue number, atom name, occ, coords, b factor

	       // add the atom to the residue and the residue to the chain
	       new_residue_p->AddAtom(new_atom_p);
	       chain_p->AddResidue(new_residue_p);
	    }
	 }
      }
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions(); // sets unsaved changes flag
   }
   
   return istat;
}

// Add this molecule (typically of waters to this
// molecule... somehow).  All the atoms of water_mol need to be in a
// chain that has a different chain id to all the chains in this
// molecule.  Else fail (return status 0).
// 
int
molecule_class_info_t::append_to_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0; // fail status initially.
   int n_atom = 0;  // 0 new atoms added initially.
   
   if (atom_sel.n_selected_atoms > 0) {

      make_backup();

      // run over the chains in water_mol (there is only one for waters)
      //
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {

// 	 std::cout << "DEBUG:: append_to_molecule: fragment id id for frag " << ifrag
// 		   << " is " << water_mol[ifrag].fragment_id << std::endl;

	 short int imatch = 0;
	 
	 // Run over chains of the existing mol, to see if there
	 // already exists a chain with the same chain id as the
	 // waters we want to add.  Only if imatch is 0 does this
	 // function do anything.
	 // 
	 int nchains = atom_sel.mol->GetNumberOfChains(1);
	 CChain *chain;
	 for (int ichain=0; ichain<nchains; ichain++) { 
	 
	    chain = atom_sel.mol->GetChain(1,ichain);
	    std::string mol_chain_id(chain->GetChainID());
	 
	    if (water_mol.fragments[ifrag].fragment_id == mol_chain_id) {
	       //
	       imatch = 1;
	       istat = 1;
	       std::cout << "INFO:: Can't add waters from additional molecule "
			 << "chain id = " << mol_chain_id << std::endl
			 << "INFO:: That chain id already exists in this molecule"
			 << std::endl;
	       break;
	    }
	 }

	 CModel *model_p = atom_sel.mol->GetModel(1);
	 if (imatch == 0) {
	    // There was not already a chain in this molecule of that name.

	    CChain *new_chain_p;
	    CAtom *new_atom_p;
	    CResidue *new_residue_p;

	    new_chain_p = new CChain;
	    std::cout << "DEBUG INFO:: chain id of new chain :"
		      << water_mol[ifrag].fragment_id << ":" << std::endl;
	    new_chain_p->SetChainID(water_mol[ifrag].fragment_id.c_str());
	    model_p->AddChain(new_chain_p);

	    for (int ires=water_mol[ifrag].min_res_no();
		 ires<=water_mol[ifrag].max_residue_number();
		 ires++) {

	       if (water_mol[ifrag][ires].atoms.size() > 0) {
		  new_residue_p = new CResidue;
		  new_residue_p->seqNum = ires;
		  strcpy(new_residue_p->name, water_mol[ifrag][ires].name.c_str());
		  new_chain_p->AddResidue(new_residue_p);
		  for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {
		     
		     new_atom_p = new CAtom;
		     new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
		     new_atom_p->SetElementName(water_mol[ifrag][ires][iatom].element.c_str());
		     new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
						water_mol[ifrag][ires][iatom].pos.y(),
						water_mol[ifrag][ires][iatom].pos.z(),
						1.0, graphics_info_t::default_new_atoms_b_factor);
		     new_residue_p->AddAtom(new_atom_p);
		     n_atom++; 
		  }
	       }
	    }
	 }
      }

      std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;
      if (n_atom > 0) { 
	 atom_sel.mol->FinishStructEdit();
	 update_molecule_after_additions(); // sets unsaved changes flag
      }
   }

   return istat;
} 


void
molecule_class_info_t::update_molecule_after_additions() {

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);

   atom_sel = make_asc(atom_sel.mol); // does the udd stuff too.

//    std::cout << "old_n_atoms: " << old_n_atoms << " new: " 
// 	     << atom_sel.n_selected_atoms << std::endl;

   have_unsaved_changes_flag = 1;
   make_bonds_type_checked();
} 

std::string
molecule_class_info_t::Refmac_in_name() const {

   return Refmac_name_stub() + "-pre.pdb";

}

std::string
molecule_class_info_t::Refmac_out_name() const {

   return Refmac_name_stub() + ".pdb";

}

std::string
molecule_class_info_t::Refmac_mtz_out_name() const {

   return Refmac_name_stub() + ".mtz";

}

// combine and strip molecule and refmac count to come up with a pdb filename
// for refmac
std::string
molecule_class_info_t::Refmac_name_stub() const {

   // shall we try to take into account ccp4i refmac naming?
   // OK:
   // Here is an example: 
   // demo.pdb         -> demo_refmac1.pdb 
   // demo_refmac1.pdb -> demo_refmac2.pdb 

   std::string refmac_name = "pre-refmac.pdb"; // default

   // First strip off the path of name_:
   std::string stripped_name; 
   // /a/b.mtz -> b.mtz
   std::string::size_type islash = name_.find_last_of("/");
   if (islash == string::npos) {
      // std::cout << "DEBUG:: slash not found in " << name_ << std::endl;
      stripped_name = name_;
   } else {
      // std::cout << "DEBUG:: slash found at " << islash << std::endl;
      // stripped_name = name_.substr(islash+1, name_.length());
      stripped_name = name_.substr(islash+1);
   } 
   // std::cout << "DEBUG:: stripped_name: " << stripped_name << std::endl;      

   
   std::string::size_type irefmac = stripped_name.rfind("-refmac");
   std::string::size_type irefmac_ccp4i = stripped_name.rfind("_refmac");
   
   if (irefmac == string::npos) { // not found

      // so was it a ccp4i refmac pdb file?

      if ( ! (irefmac_ccp4i == string::npos) ) { 
	 // it *was* a ccp4i pdb file:
	 //
	 refmac_name = stripped_name.substr(0,irefmac_ccp4i) + "_refmac";
	 refmac_name += graphics_info_t::int_to_string(refmac_count);
      }
      // std::cout << "DEBUG:: irefmac not found in " << stripped_name  << std::endl;
      // lets strip off ".pdb", ".pdb.gz"
      std::string::size_type ipdb = stripped_name.rfind(".pdb");

      if (ipdb == string::npos) { // not a pdb

	 // std::cout << "DEBUG:: ipdb not found" << std::endl;
	 // just tack "refmac-2.pdb" on to the name then
	 refmac_name = stripped_name + "_refmac"; 
	 refmac_name += graphics_info_t::int_to_string(refmac_count); 

      } else {
	 // is a pdb:

	 // std::cout << "DEBUG:: ipdb *was* found" << std::endl;
	 refmac_name = stripped_name.substr(0,ipdb) + "_refmac"; 
	 refmac_name += graphics_info_t::int_to_string(refmac_count); 
      }
   } else {

      // refmac *was* found as part of the name
      // std::cout << "DEBUG:: irefmac *was* found in " << stripped_name << std::endl;
      refmac_name = stripped_name.substr(0,irefmac) + "_refmac";
      refmac_name += graphics_info_t::int_to_string(refmac_count);
   }

   // std::cout << "DEBUG:: returning refmac_name: " << refmac_name << std::endl;
   return refmac_name;

}


std::string
molecule_class_info_t::name_sans_extension(short int include_path_flag) const {

   std::string outstring = name_;
   
   std::string::size_type ipdb = name_.rfind(".pdb");
   if (ipdb != std::string::npos)
      outstring = name_.substr(0, ipdb);

   std::string::size_type islash = outstring.rfind("/");
   if (islash != std::string::npos)
      outstring = outstring.substr(islash+1);
      
   return outstring;
}

void
molecule_class_info_t::update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position) {

   if (has_model()) {
      CModel *model_p = atom_sel.mol->GetModel(1); 
      
      if (! model_p) { 
	 std::cout << "ERROR:: Disaster in finding model_p in update_molecule_to" 
		   << std::endl;
      } else {
	 CChain *chain_p;
	 int n_chains = atom_sel.mol->GetNumberOfChains(1);
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    model_p->DeleteChain(i_chain);
	 }

	 // Now add new chain to model:
	 chain_p = new CChain;
	 if (chain_p) {
	    model_p->AddChain(chain_p);
	    add_multiple_dummies(chain_p, pos_position);
	 } else { 
	    std::cout << "ERROR:: creating chain in mol::update_molecule_to" << std::endl;
	 }
      }
   }
}

// function callable from graphics_info_t::create_molecule_and_display()
// 
void
molecule_class_info_t::add_multiple_dummies(const std::vector<coot::scored_skel_coord> &pos_position) {

   if (has_model()) {
      CModel *model_p = atom_sel.mol->GetModel(1);
      int n_chains = atom_sel.mol->GetNumberOfChains(1);
      if (n_chains > 0) {
	 CChain *chain_p = model_p->GetChain(0);
	 add_multiple_dummies(chain_p, pos_position);
      }
   }
}


// we presume that the chain exists.  This exists so that we dont do a
// backup every time we add a dummy atom (as is done using
// add_dummy_atom().
// 
void
molecule_class_info_t::add_multiple_dummies(CChain *chain_p,
					    const std::vector<coot::scored_skel_coord> &pos_position) {


   if (pos_position.size() > 0) {
      make_backup(); // maybe
   }
   
   for (unsigned int i=0; i<pos_position.size(); i++) {
      CResidue *res_p = new CResidue;
      CAtom *atom_p = new CAtom;
      chain_p->AddResidue(res_p);
      atom_p->SetAtomName(" DUM");
      atom_p->SetCoordinates(pos_position[i].position.x(),
			     pos_position[i].position.y(),
			     pos_position[i].position.z(), 1.0,
			     graphics_info_t::default_new_atoms_b_factor);

      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      res_p->seqNum = i + 1;
      res_p->SetResName("DUM");

      // std::cout << atom_p << " added to molecule" << std::endl;
   }

   // std::cout << "DEBUG:: add_multiple_dummies finishing.. "
   // << pos_position.size() << std::endl;
   // if (pos_position.size() > 0) {
   
   // Actually, we want to run this code when there are no new guide
   // points too.  This sets atom_sel.SelectionHandle properly, which
   // is needed in close_yourself, where a DeleteSelection() is done
   // to give back the memory.
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1; 
      makebonds(0.0, 0.0);
      // }
} 

void
molecule_class_info_t::add_multiple_dummies(const std::vector<coot::Cartesian> &pos_position) {

   if (atom_sel.mol) {
      CModel *model_p = atom_sel.mol->GetModel(1);
      int n_chains = atom_sel.mol->GetNumberOfChains(1);
      if (n_chains > 0) {
	 CChain *chain_p = model_p->GetChain(0);
	 if (pos_position.size() > 0) {
	    make_backup(); // maybe
	    
	    for (unsigned int i=0; i< pos_position.size(); i++) { 
	       CResidue *res_p = new CResidue;
	       CAtom *atom_p = new CAtom;
	       chain_p->AddResidue(res_p);
	       atom_p->SetAtomName(" DUM");
	       atom_p->SetCoordinates(pos_position[i].x(),
				      pos_position[i].y(),
				      pos_position[i].z(), 1.0,
				      graphics_info_t::default_new_atoms_b_factor);
	    
	       atom_p->SetElementName(" O");
	       res_p->AddAtom(atom_p);
	       res_p->seqNum = i + 1;
	       res_p->SetResName("DUM");
	    }
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1; 
	    makebonds(0.0, 0.0);
	 }
      }
   }
}


// return an empty vector on failure, a vector of size 6 on success:
// 
std::pair<std::vector<float>, std::string> 
molecule_class_info_t::get_cell_and_symm() const { 

   std::pair<std::vector<float>, std::string> cell_spgr;
   
   mat44 my_matt;
   if (atom_sel.mol) { 
      int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (err != 0) {
	 std::cout << "!! Warning:: No symmetry available for this template molecule"
		   << std::endl;
      } else {
	 realtype a[6];
	 realtype vol;
	 int orthcode;
	 atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
	 for (int i=0; i<6; i++) cell_spgr.first.push_back(a[i]);
	 cell_spgr.second = std::string(atom_sel.mol->GetSpaceGroup());
      }
   }
   return cell_spgr;
} 

void 
molecule_class_info_t::set_mmdb_cell_and_symm(std::pair<std::vector<float>, std::string> cell_spgr) {

   if (cell_spgr.first.size() == 6) { 
      std::vector<float> a = cell_spgr.first; // short name
      atom_sel.mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
      atom_sel.mol->SetSpaceGroup((char *)cell_spgr.second.c_str());
      std::cout << "successfully set cell and symmetry" << std::endl;
   } else { 
      std::cout << "WARNING:: failure to set cell on this molecule" << std::endl;
   } 
}

void
molecule_class_info_t::set_mmdb_symm(const std::string &spg) {

   atom_sel.mol->SetSpaceGroup((char *)spg.c_str());

} 

// Return atom_index of -1 when no nearest atom.
// 
std::pair<float, int>
molecule_class_info_t::nearest_atom(const coot::Cartesian &pos) const { 

   float min_dist = 999999999; 
   float d;
   int atom_index = -1;

   for(int i=0; i<atom_sel.n_selected_atoms; i++) { 
      coot::Cartesian a(atom_sel.atom_selection[i]->x, atom_sel.atom_selection[i]->y, atom_sel.atom_selection[i]->z);
      d =  fabs((pos-a).length());
      if (d < min_dist) { 
	 min_dist = d;
	 atom_index = i;
      }
   }

   std::pair<float, int> r;
   r.first = min_dist;
   r.second = atom_index;
   return r;
}

// return an empty
// vector
// if closed or is a coords
// mol, 3 elements and 6
// elements for a
// difference map.
std::vector<float>
molecule_class_info_t::map_colours() const {

   std::vector<float> v;
   if (has_map()) {
      if (is_difference_map_p()) {
	 v.resize(6,0.3);
	 v[0] = map_colour[0][0];
	 v[1] = map_colour[0][1];
	 v[2] = map_colour[0][2];
	 v[3] = map_colour[1][0];
	 v[4] = map_colour[1][1];
	 v[5] = map_colour[1][2];
      } else {
	 v.resize(3,0.3);
	 v[0] = map_colour[0][0];
	 v[1] = map_colour[0][1];
	 v[2] = map_colour[0][2];
      }
   }
   return v;
}

// perhaps there is a better place for this?
// 
std::vector<std::string> 
molecule_class_info_t::set_map_colour_strings() const { 

   // return something like
   // (list "set_last_map_colour" "0.2" "0.3" "0.4")


   std::vector<std::string> r;

   r.push_back("set-last-map-colour");
   r.push_back(graphics_info_t::float_to_string(map_colour[0][0]));
   r.push_back(graphics_info_t::float_to_string(map_colour[0][1]));
   r.push_back(graphics_info_t::float_to_string(map_colour[0][2]));

   return r;   
} 


// and symm labels.
void
molecule_class_info_t::remove_atom_labels() { 

//    std::cout << "DEBUG:: n_labelled_atoms " << n_labelled_atoms 
// 	    << " " << n_labelled_symm_atoms << std::endl;

   int init_count = n_labelled_atoms;
   for (int i=0; i<init_count; i++) { 
      unlabel_atom(labelled_atom_index_list[0]);
   }

   init_count = n_labelled_symm_atoms;
   for (int i=0; i<init_count; i++) {
      unlabel_symm_atom(labelled_symm_atom_index_list[0]);
   }
} 


// So that we can move around all the atoms of a ligand (typically)
// 
void 
molecule_class_info_t::translate_by(float x, float y, float z) { 

   if (has_model()) {
      make_backup();
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 atom_sel.atom_selection[i]->x += x;
	 atom_sel.atom_selection[i]->y += y;
	 atom_sel.atom_selection[i]->z += z;
      }
      make_bonds_type_checked();
      have_unsaved_changes_flag = 1;
   }
} 

// Sets coot_save_index maybe (using set_coot_save_index()).
// 
std::string 
molecule_class_info_t::stripped_save_name_suggestion() { 

   std::string s;
   
   std::string stripped_name1;
   std::string::size_type islash = name_.find_last_of("/");
   if (islash == std::string::npos) { 
      stripped_name1 = name_;
   } else {
      stripped_name1 = name_.substr(islash+1, name_.length());
   }
   // so we have got rid of the pathname.
   // now lets get rid of the extension
   // 
   std::string::size_type ibrk   = stripped_name1.rfind(".brk");
   std::string::size_type ibrkgz = stripped_name1.rfind(".brk.gz");
   std::string::size_type ipdb   = stripped_name1.rfind(".pdb");
   std::string::size_type ires   = stripped_name1.rfind(".res");
   std::string::size_type ipdbgz = stripped_name1.rfind(".pdb.gz");
   std::string::size_type icoot  = stripped_name1.rfind("-coot-");

   // std::cout << "DEBUG:: icoot: " << icoot << " ipdb " << ipdb << std::endl;
   
   std::string stripped_name2;
   if (icoot == std::string::npos) { 
      if (ibrk == std::string::npos) { 
	 if (ibrkgz == std::string::npos) { 
	    if (ipdb == std::string::npos) { 
	       if (ires == std::string::npos) { 
		  if (ipdbgz == std::string::npos) { 
		     stripped_name2 = stripped_name1;
		  } else {
		     stripped_name2 = stripped_name1.substr(0, ipdbgz);
		  }
	       } else {
		  stripped_name2 = stripped_name1.substr(0, ires);
	       }
	    } else {
	       stripped_name2 = stripped_name1.substr(0, ipdb);
	    }
	 } else {
	    stripped_name2 = stripped_name1.substr(0, ibrkgz);
	 }
      } else { 
	 stripped_name2 = stripped_name1.substr(0, ibrk);
      }
   } else {
      set_coot_save_index(stripped_name1.substr(icoot));
      stripped_name2 = stripped_name1.substr(0,icoot);
   }

   stripped_name2 += "-coot-";
   stripped_name2 += graphics_info_t::int_to_string(coot_save_index);
   // As per George Sheldrick's suggestion, if this was from shelx,
   // suggest a .ins extension, not .pdb
   if (!is_from_shelx_ins_flag) {
      stripped_name2 += ".pdb";
   } else { 
      stripped_name2 += ".ins";
   }

//    std::cout << "DEBUG:: stripped_save_name_suggestion: " 
// 	     << stripped_name2 << std::endl;

   return stripped_name2;
} 

int 
molecule_class_info_t::set_coot_save_index(const std::string &filename) { 

   std::cout << "extracting from :" << filename << std::endl;

   // filename is something like: "-coot-12.pdb".
   // 
   // We want to find 12 and set coot_save_index to 12 + 1, which gets
   // used to suggest the next saved filename.
   // 

   std::string twelve_pdb = filename.substr(6);
   // std::cout << "twelve_pdb:"<< twelve_pdb << std::endl;
   
   std::string::size_type ipdb   = twelve_pdb.rfind(".pdb");
   if (ipdb != std::string::npos) { 
      // .pdb was found
      std::string twelve = twelve_pdb.substr(0,ipdb);
      int i = atoi(twelve.c_str());
      // std::cout << "found i: " << i << std::endl;
      if (i >= 0 && i<100000)
	 coot_save_index = i+1;
   } 
   return coot_save_index; 
} 


void
molecule_class_info_t::transform_by(mat44 mat) { 

   if (has_model()) { 
      clipper::Coord_orth co;
      clipper::Coord_orth trans_pos; 
      make_backup();
      clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
					 mat[1][0], mat[1][1], mat[1][2],
					 mat[2][0], mat[2][1], mat[2][2]);
      clipper::Coord_orth cco(mat[0][3], mat[1][3], mat[2][3]);
      clipper::RTop_orth rtop(clipper_mat, cco);
      std::cout << "INFO:: coordinates transformed by orthonal matrix: \n"
		<< rtop.format() << std::endl;
      for (int i=0; i<atom_sel.n_selected_atoms; i++) { 
	 // atom_sel.atom_selection[i]->Transform(mat); // doesn't compile!
	 // Argh.  sigh.  Use clipper. c.f. graphics_info_t::fill_hybrid_atoms()
	 co = clipper::Coord_orth(atom_sel.atom_selection[i]->x, 
				  atom_sel.atom_selection[i]->y, 
				  atom_sel.atom_selection[i]->z);
	 trans_pos = co.transform(rtop);
	 atom_sel.atom_selection[i]->x = trans_pos.x();
	 atom_sel.atom_selection[i]->y = trans_pos.y();
	 atom_sel.atom_selection[i]->z = trans_pos.z();
      } 
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
}


void
molecule_class_info_t::transform_by(const clipper::RTop_orth &rtop) {

   make_backup();
   std::cout << "INFO:: coordinates transformed by orthonal maxix: \n"
	     << rtop.format() << std::endl;
   if (have_unit_cell) {
      clipper::Cell cell(clipper::Cell_descr(atom_sel.mol->get_cell().a,
					     atom_sel.mol->get_cell().b,
					     atom_sel.mol->get_cell().c,
					     clipper::Util::d2rad(atom_sel.mol->get_cell().alpha),
					     clipper::Util::d2rad(atom_sel.mol->get_cell().beta),
					     clipper::Util::d2rad(atom_sel.mol->get_cell().gamma))); 
      std::cout << "INFO:: fractional coordinates matrix:" << std::endl;
      std::cout << rtop.rtop_frac(cell).format() << std::endl;
   } else {
      std::cout << "No unit cell for this molecule, hence no fractional matrix." << std::endl;
   }
   clipper::Coord_orth co;
   clipper::Coord_orth trans_pos; 
   if (has_model()) { 
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 co = clipper::Coord_orth(atom_sel.atom_selection[i]->x, 
				  atom_sel.atom_selection[i]->y, 
				  atom_sel.atom_selection[i]->z);
	 trans_pos = co.transform(rtop);
	 atom_sel.atom_selection[i]->x = trans_pos.x();
	 atom_sel.atom_selection[i]->y = trans_pos.y();
	 atom_sel.atom_selection[i]->z = trans_pos.z();
      }
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
}

void
molecule_class_info_t::transform_by(const clipper::RTop_orth &rtop, CResidue *residue_moving) {

   make_backup();
   std::cout << "INFO:: coordinates transformed_by: \n"
	     << rtop.format() << std::endl;
   clipper::Coord_orth co;
   clipper::Coord_orth trans_pos; 
   if (has_model()) {

      PPCAtom residue_atoms;
      int n_residue_atoms;
      residue_moving->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iatom=0; iatom<n_residue_atoms; iatom++) { 
	 clipper::Coord_orth p(residue_atoms[iatom]->x,
			       residue_atoms[iatom]->y,
			       residue_atoms[iatom]->z);
	 clipper::Coord_orth p2 = p.transform(rtop);
	 residue_atoms[iatom]->x = p2.x();
	 residue_atoms[iatom]->y = p2.y();
	 residue_atoms[iatom]->z = p2.z();
      }
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
}




void
molecule_class_info_t::set_contour_level(float f) {
   contour_level[0] = f;
}

void
molecule_class_info_t::set_contour_level_by_sigma(float f) {

   contour_level[0] = f * map_sigma_;
}

std::vector <std::string>
molecule_class_info_t::get_map_contour_strings() const {

   std::vector <std::string> s; 
   s.push_back("set-last-map-contour-level");
   char cs[100];
   snprintf(cs, 99, "%e", contour_level[0]);
   s.push_back(cs);

   return s;
} 

std::vector <std::string>
molecule_class_info_t::get_map_contour_sigma_step_strings() const {

   std::vector <std::string> s; 
   s.push_back("set-last-map-sigma-step");
   s.push_back(graphics_info_t::float_to_string(contour_sigma_step));
   
   return s;
}

short int
molecule_class_info_t::contoured_by_sigma_p() const { 
   return contour_by_sigma_flag;
} 


// Return status, was the contour level changed?  In that way, we
// don't try to recontour (which is a slow process) when the contour
// level has not been changed.
// 
// We don't change the contour level if the contour level goes too
// low (typically below 0).
//
// We don't change the contour level if the contour level goes too
// high (above the maximum level of the map).
// 
short int
molecule_class_info_t::change_contour(int direction) {

   short int istat = 0;
   // std::cout << "DEBUG:: contour_by_sigma_flag " << contour_by_sigma_flag << std::endl;
   // std::cout << "DEBUG:: adding " << contour_sigma_step << " * " << map_sigma_
   // << " to  " << contour_level[0] << std::endl;
   if (has_map()) {

      float shift = graphics_info_t::diff_map_iso_level_increment;
      if (contour_by_sigma_flag) { 
	 shift = contour_sigma_step * map_sigma_;
      } else { 
	 if (xmap_is_diff_map[0]) { 
	    shift = graphics_info_t::diff_map_iso_level_increment;
	 } else { 
	    shift = graphics_info_t::iso_level_increment;
	 }
      }

      if (xmap_is_diff_map[0]) {
	 if (direction == -1) {
	    if (graphics_info_t::stop_scroll_diff_map_flag) {
	       if ((contour_level[0] - shift) > 
		   graphics_info_t::stop_scroll_diff_map_level) { 
		  contour_level[0] -= shift;
		  istat = 1;
	       }
	    } else {
	       contour_level[0] -= shift;
	       istat = 1;
	    }
	 } else {
	    // add, but don't go past the top of the map or the bottom of the map
	    // 
	    if (contour_level[0] <= map_max_ || contour_level[0] <= -map_min_) {
	       contour_level[0] += shift;
	       istat = 1;
	    }
	 }
      } else {
	 // iso map

	 if (direction == -1) {
	    if (graphics_info_t::stop_scroll_iso_map_flag) {
	       if ((contour_level[0] - shift) >
		   graphics_info_t::stop_scroll_iso_map_level) {
		  contour_level[0] -= shift;
		  istat = 1;
	       }
	    } else {
	       contour_level[0] -= shift;
	       istat = 1;
	    }
	 } else {
	    if (contour_level[0] <= map_max_) {
	       contour_level[0] += shift;
	       istat = 1;
	    }
	 } 
      }
   }
   return istat;
}

CChain *
molecule_class_info_t::water_chain() const { 

   CChain *water_chain = 0;

   if (has_model()) { 

      CModel *model_p = atom_sel.mol->GetModel(1);
      CResidue *residue_p;
      CChain *chain_p;
      
      if (model_p) {

	 if (is_from_shelx_ins_flag) {
	    water_chain = water_chain_from_shelx_ins();
	 } else { 
	    int nchains = model_p->GetNumberOfChains();
	    for (int ich=0; ich<nchains; ich++) {
	       chain_p = model_p->GetChain(ich);
	       int nres = chain_p->GetNumberOfResidues();
	       short int all_water_flag = 1; 
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  std::string resname(residue_p->name);
		  if (! ( (resname == "WAT") || (resname == "HOH"))) {
		     all_water_flag = 0;
		     break;
		  }
	       }
	       if (all_water_flag) { 
		  water_chain = chain_p;
		  break;
	       }
	    }
	 }
      } 
   }
   return water_chain;
}


// there is only one chain from a shelxl ins file.
CChain *
molecule_class_info_t::water_chain_from_shelx_ins() const {

   CChain *water_chain = 0;
   CModel *model_p = atom_sel.mol->GetModel(1);

   if (has_model()) { 
      int nchains = model_p->GetNumberOfChains();
      for (int ich=0; ich<nchains; ich++) {
	 water_chain = model_p->GetChain(ich);
      }
   }
   return water_chain;
}


// return state, max_resno + 1, or 0, 1 of no residues in chain.
// 
std::pair<short int, int> 
molecule_class_info_t::next_residue_in_chain(CChain *w) const { 

   std::pair<short int, int> p(0,1);
   int max_res_no = -9999;

   if (w) { 
      int nres = w->GetNumberOfResidues();
      CResidue *residue_p;
      if (nres > 0) { 
	 for (int ires=nres-1; ires>=0; ires--) { 
	    residue_p = w->GetResidue(ires);
	    if (residue_p->seqNum > max_res_no) {
	       max_res_no = residue_p->seqNum;
	       p = std::pair<short int, int>(1, max_res_no+1);
	    }
	 }
      }
   }
   return p;
}


// 
void
molecule_class_info_t::set_map_is_difference_map() { 

   if (has_map()) { 
      xmap_is_diff_map[0] = 1;
      update_map();
   }
}

short int
molecule_class_info_t::is_difference_map_p() const {

   short int istat = 0;
   if (has_map())
      if (xmap_is_diff_map[0])
	 istat = 1;
   return istat;
}


void
molecule_class_info_t::set_contour_by_sigma_step(float v, short int state) { 
   contour_by_sigma_flag = state;
   if (state)
      contour_sigma_step = v;
}




// add a factor to scale the colours in b factor representation:.
// It goes into the atom_sel.mol
void
molecule_class_info_t::set_b_factor_bonds_scale_factor(float f) {

   std::cout << "Here Adding b-factor scale " << f << std::endl;
   if (atom_sel.mol) {
      // bleugh, casting.
      int udd_handle =
	 atom_sel.mol->RegisterUDReal(UDR_HIERARCHY,
				      (char *) coot::b_factor_bonds_scale_handle_name.c_str());
      if (udd_handle > 0) {
// 	 std::cout << "Adding b-factor scale " << f << " with handle "
// 		   << udd_handle << std::endl;
	 atom_sel.mol->PutUDData(udd_handle, f);

	 // test getting the uddata:
	 int udd_b_factor_handle =
	    atom_sel.mol->GetUDDHandle(UDR_HIERARCHY, (char *) coot::b_factor_bonds_scale_handle_name.c_str());
// 	 std::cout << "debug:: test Got b factor udd handle: "
// 		   << udd_b_factor_handle << std::endl;
	 if (udd_b_factor_handle > 0) {
	    realtype scale;
	    if (atom_sel.mol->GetUDData(udd_b_factor_handle, scale) == UDDATA_Ok) {
// 	       std::cout << " test got b factor scale: " << scale << std::endl;
	    } else {
 	       std::cout << "ERROR:: bad get b factor scale " << std::endl;
	    }
	 }
      }
   }
   make_bonds_type_checked();
}
