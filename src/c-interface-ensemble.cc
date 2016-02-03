/* src/c-interface-ensemble.cc
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#include <stdlib.h>
#include <iostream>

#if defined _MSC_VER
#include <windows.h>
#endif

 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // functions are built by glade.

#include <vector>
#include <string>
#include <algorithm>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h"

#include "graphics-info.h"

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"

#include "guile-fixups.h"


/*  ----------------------------------------------------------------------- */
/*                  Ensemble view                                           */
/*  ----------------------------------------------------------------------- */
// BL says:: wonder if this should be streamlined to have one ghost function
// which serves both NCS and ensembles as to avoid some code dublication
// FIXME:: not all headers are required here...

void set_draw_ensemble_ghosts(int imol, int istate) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_ensemble_ghosts(istate);
      graphics_draw();
   }
}

/*! \brief return the drawing state of ensemble ghosts for molecule number imol   */
int draw_ensemble_ghosts_state(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].draw_ensemble_ghosts_p();
   }
   return r;
} 


void set_ensemble_ghost_bond_thickness(int imol, float f) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_ensemble_ghost_bond_thickness(f);
      graphics_draw();
   }
}


void ensemble_update_ghosts(int imol) {

   if (is_valid_model_molecule(imol)) {
      int incs = graphics_info_t::molecules[imol].update_ensemble_ghosts_size();
      if (incs)
	 graphics_draw();
   }
}

// BL says:: I don tthink we need this certainly FIXME!
// The NCS "Yes" was pressed, do we need to make the rtops for this
// molecule?  We do if there ncs_ghosts.size() > 0.
//
// If ncs_ghosts.size() == 0, then this widget should not be active.
// 
// 
void
make_ensemble_ghosts_maybe(int imol) {

   if (is_valid_model_molecule(imol)) {  // it should be!
      if (graphics_info_t::molecules[imol].has_models_p()) {
         graphics_info_t::molecules[imol].fill_ensemble_info();
      }
   }
}


// BL says:: could have some ensemble model difference functions?!

//#ifdef USE_GUILE
//SCM ncs_chain_differences_scm(int imol, const char *master_chain_id) {

//   float mc_weight = 1.0;
//   SCM r = SCM_BOOL_F;
//   if (is_valid_model_molecule(imol)) {
//      coot::ncs_differences_t diffs =
//	 graphics_info_t::molecules[imol].ncs_chain_differences(master_chain_id,
//								mc_weight);
//      if (diffs.size() == 0) {
//	 std::cout << "no diffs" << std::endl;
//      } else {
//	 r = SCM_EOL;
//	 for (int idiff=int(diffs.size()-1); idiff>=0; idiff--) {
//	    SCM l_residue_data = SCM_EOL;
//	    coot::ncs_chain_difference_t cd = diffs.diffs[idiff];
//	    if (cd.residue_info.size() > 0) {
//	       std::cout << "NCS target chain has " << cd.residue_info.size()
//			 << " peers." << std::endl;
//	       //	       for (int iresinf=0; iresinf<cd.residue_info.size(); iresinf++) {
//	       for (int iresinf=(int(cd.residue_info.size())-1); iresinf>=0; iresinf--) {
//		  if (0)
//		     std::cout << "resinfo: "
//			       << cd.residue_info[iresinf].resno << " "
//			       << cd.residue_info[iresinf].inscode << " "
//			       << cd.residue_info[iresinf].serial_number << " to "
//			       << cd.residue_info[iresinf].target_resno << " "
//			       << cd.residue_info[iresinf].target_inscode << " "
//			       << cd.residue_info[iresinf].target_serial_number << " diff: "
//			       << cd.residue_info[iresinf].mean_diff
//			       << std::endl;
//		  coot::residue_spec_t this_res(cd.peer_chain_id,
//						cd.residue_info[iresinf].resno,
//						cd.residue_info[iresinf].inscode);
//		  coot::residue_spec_t target_res(diffs.target_chain_id,
//						  cd.residue_info[iresinf].target_resno,
//						  cd.residue_info[iresinf].target_inscode);
//		  SCM res_l = SCM_EOL;
//		  res_l = scm_cons(scm_double2num(cd.residue_info[iresinf].mean_diff), res_l);
////		  res_l = scm_cons(scm_cdr(scm_residue(target_res)), res_l);
////		  res_l = scm_cons(scm_cdr(scm_residue(this_res)), res_l);
//		  l_residue_data = scm_cons(res_l, l_residue_data);
//	       }
//	       r = scm_cons(l_residue_data, SCM_EOL);
//	       r = scm_cons(scm_makfrom0str(diffs.target_chain_id.c_str()), r);
//	       r = scm_cons(scm_makfrom0str(cd.peer_chain_id.c_str()), r);
//	    }
//	 }
//      }
//   }
//   return r;
//}
//#endif	/* USE_GUILE */
//#ifdef USE_PYTHON
//PyObject *ncs_chain_differences_py(int imol, const char *master_chain_id) {

//   PyObject *r = Py_False;
//   if (is_valid_model_molecule(imol)) {
//      coot::ncs_differences_t diffs =
//	 graphics_info_t::molecules[imol].ncs_chain_differences(master_chain_id, 1.0);
//      if (diffs.size() == 0) {
//	 std::cout << "no diffs" << std::endl;
//      } else {
//	 r = PyList_New(0);
//	 for (unsigned int idiff=0; idiff<diffs.size(); idiff++) {
//	    PyObject *l_residue_data = PyList_New(0);
//	    coot::ncs_chain_difference_t cd = diffs.diffs[idiff];
//	    if (cd.residue_info.size() > 0) {
//	       std::cout << "NCS target chain has " << cd.residue_info.size()
//			 << " peers." << std::endl;
//	       //	       for (int iresinf=0; iresinf<cd.residue_info.size(); iresinf++) {
//	       for (unsigned int iresinf=0; iresinf<cd.residue_info.size(); iresinf++) {
//		  if (0)
//		     std::cout << "resinfo: "
//			       << cd.residue_info[iresinf].resno << " "
//			       << cd.residue_info[iresinf].inscode << " "
//			       << cd.residue_info[iresinf].serial_number << " to "
//			       << cd.residue_info[iresinf].target_resno << " "
//			       << cd.residue_info[iresinf].target_inscode << " "
//			       << cd.residue_info[iresinf].target_serial_number << " diff: "
//			       << cd.residue_info[iresinf].mean_diff
//			       << std::endl;
//		  coot::residue_spec_t this_res(cd.peer_chain_id,
//						cd.residue_info[iresinf].resno,
//						cd.residue_info[iresinf].inscode);
//		  coot::residue_spec_t target_res(diffs.target_chain_id,
//						  cd.residue_info[iresinf].target_resno,
//						  cd.residue_info[iresinf].target_inscode);
//		  // according to Paul's documentation we should have resno and inscode
//		  // for both residues here too
//		  PyObject *thisr = PyList_GetSlice(py_residue(this_res), 2, 4);
		  
//		  PyObject *masta = PyList_GetSlice(py_residue(target_res), 2, 4);

//		  // res_l list seems only to have one element in Paul's scm code
//		  // currently?! Correct?!
//		  PyObject *res_l = PyList_New(3);
//		  PyList_SetItem(res_l, 0, thisr);
//		  PyList_SetItem(res_l, 1, masta);
//		  PyList_SetItem(res_l, 2, PyFloat_FromDouble(cd.residue_info[iresinf].mean_diff));
//		  PyList_Append(l_residue_data, res_l);
//		  Py_XDECREF(res_l);
//	       }
//	       PyList_Append(r, PyString_FromString(cd.peer_chain_id.c_str()));
//	       PyList_Append(r, PyString_FromString(diffs.target_chain_id.c_str()));
//	       PyList_Append(r, l_residue_data);
//	       Py_XDECREF(l_residue_data);
//	    }
//	 }
//      }
//   }
//   if (PyBool_Check(r)) {
//     Py_INCREF(r);
//   }
//   return r;
//}
//#endif	/* USE_PYTHON */


//// This should be  in c-interface-ncs-gui.cc
//void validation_graph_ncs_diffs_mol_selector_activate (GtkMenuItem     *menuitem,
//						      gpointer         user_data) {
   
//   int imol = GPOINTER_TO_INT(user_data);
//#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
//#ifdef HAVE_GSL
//   graphics_info_t g;
//   g.ncs_diffs_from_mol(imol);
//#else
//   printf("not compiled with GSL - remake\n");
//#endif /* HAVE_GSL */
//#else
//   printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n");
//#endif /* HAVE_GTK_CANVAS */

//}

