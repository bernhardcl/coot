

#ifdef USE_MOLECULES_TO_TRIANGLES

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <cstddef>

// #include "globjects.h" //includes gtk/gtk.h
#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "cc-interface-molecular-representation.hh"

#ifdef USE_PYTHON

// Martin's Triangles

// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_py(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
#ifdef USE_MOLECULES_TO_TRIANGLES
      // check that these are strings
      std::string atom_selection = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_selection_py));
      std::string ColorScheme    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ColorScheme_py));
      std::string style          = PyBytes_AS_STRING(PyUnicode_AsUTF8String(style_py));
      // status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
      graphics_info_t g;
      status = g.add_molecular_representation(imol, atom_selection, ColorScheme, style);
      graphics_draw();
#endif
   }
   return status;

}
#endif // USE_PYTHON

#ifdef USE_GUILE
// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_scm(int imol, SCM atom_selection_scm, SCM ColorScheme_scm, SCM style_scm) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
#ifdef USE_MOLECULES_TO_TRIANGLES
      // check that these are strings
      std::string atom_selection = scm_to_locale_string(atom_selection_scm);
      std::string ColorScheme    = scm_to_locale_string(ColorScheme_scm);
      std::string style          = scm_to_locale_string(style_scm);
      std::cout << "Calling add_molecular_representation with " << atom_selection << " " << ColorScheme << "  " << style << std::endl;
      // status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
      graphics_info_t g;
      status = g.add_molecular_representation(imol, atom_selection, ColorScheme, style);
      graphics_draw();
#endif
   }
   return status;
}
#endif // USE_GUILE

void remove_molecular_representation(int imol, int rep_no) {

   if (is_valid_model_molecule(imol)) {
      // graphics_info_t::molecules[imol].remove_molecular_representation(rep_no);
      graphics_info_t g;
      g.remove_molecular_representation(imol, rep_no);
      graphics_draw();
   }
}

extern "C" void add_molecular_representation_test() {
#ifdef USE_MOLECULES_TO_TRIANGLES
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
   if (active_atom.first) {
      int imol = active_atom.second.first;
      std::cout << "Ribbons on molecule " << imol << std::endl;
      if (is_valid_model_molecule(imol)) {
         std::string atom_selection = "//A";
         std::string ColorScheme = "colorRampChainsScheme";
         std::string style = "Ribbon";
         // status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
         graphics_info_t g;
         g.add_molecular_representation(imol, atom_selection, ColorScheme, style);
         graphics_info_t::graphics_draw();
      }
   }
#endif
}

#else

#endif // USE_MOLECULES_TO_TRIANGLES

