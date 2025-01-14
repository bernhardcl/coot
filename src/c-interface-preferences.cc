/* src/c-interface-preferences.cc
 * 
 * Copyright 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007 The University of Oxford
 * Copyright 2014, 2015 by Medical Research Council
 * Author: Paul Emsley, Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include <algorithm>

#include <string.h> // strlen, strcpy
#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#include <dirent.h>   // for extra scheme dir
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#include <windows.h>
#include <direct.h>
#include "compat/dirent.h"   // for extra scheme dir
#endif // _MSC_VER


#include "guile-fixups.h"

#include "graphics-info.h"
#include "interface.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "c-interface-preferences.h"
#include "c-interface-widgets.hh"
#include "cc-interface.hh"
#include "c-interface-gui.hh" // for set_transient_for_main_window()
#include "coot-preferences.h"

#include "widget-from-builder.hh"

/* insert some new stuff with the new method */

// Define preference value types
using preferences_value = std::variant<int, double, std::string, bool, std::vector<float>>;

preferences_manager coot_preferences;

// Register a preference with setter/getter functions
void preferences_manager::register_preference(
      const std::string& key,
      const std::function<void(const preferences_value&)>& setter,
      const std::function<preferences_value()>& getter,
      const preferences_value& default_value) {

   if (preferences_registry.find(key) != preferences_registry.end()) {
      throw std::runtime_error("Preference already registered: " + key);
   }

   // Store the custom getter and setter
   preferences_registry[key] = {setter, getter};

   // Set the default value
   preferences_defaults[key] = default_value;
   try {
      setter(default_value); // Initialize with the default value
   } catch (std::bad_variant_access) {
      std::cout<<"BL INFO:: try to set the wrong value for key " << key <<
                 " should be return of setter but is index " << default_value.index()<<std::endl;
   }
}

// Set a preference value
void preferences_manager::set_preference(const std::string& key, const preferences_value& value) {
   auto it = preferences_registry.find(key);
   if (it != preferences_registry.end()) {
      it->second.preference_set_function(value); // Call the setter function
   } else {
      throw std::runtime_error("Preference key not registered: " + key);
   }
}

// Get a preference value
preferences_value preferences_manager::get_preference(const std::string& key) const {
   auto it = preferences_registry.find(key);
   if (it != preferences_registry.end()) {
      return it->second.preference_get_function(); // Call the getter function
   } else if (preferences_defaults.find(key) != preferences_defaults.end()) {
      return preferences_defaults.at(key); // Return default if not explicitly set
   }
   throw std::runtime_error("Preference key not registered: " + key);
}

// Reset a preference to its default value
void preferences_manager::reset_preference_to_default(const std::string& key) {
   auto def_it = preferences_defaults.find(key);
   if (def_it != preferences_defaults.end()) {
      set_preference(key, def_it->second); // Use the default value
   } else {
      throw std::runtime_error("No default value for preference: " + key);
   }
}

// List all registered preferences
// do we need this function?
void preferences_manager::list_preferences() const {
   for (const auto& [key, callbacks] : preferences_registry) {
      std::cout << key << " = ";
      try {
      auto value = get_preference(key);

      // other way?
      // Serialize the value to Python syntax
      if (std::holds_alternative<bool>(value)) {
         std::cout << (std::get<bool>(value) ? "True" : "False");
      } else if (std::holds_alternative<int>(value)) {
          std::cout << std::get<int>(value);
      } else if (std::holds_alternative<double>(value)) {
         std::cout << std::get<double>(value); // could set with some precision but not needed here!?
      } else if (std::holds_alternative<std::string>(value)) {
          std::cout << "\"" << std::get<std::string>(value) << "\"";
      } else if (std::holds_alternative<std::vector<float>>(value)) {
          const auto& vec = std::get<std::vector<float>>(value);
          std::cout << "[";
          for (size_t i = 0; i < vec.size(); ++i) {
              std::cout << vec[i];
              if (i < vec.size() - 1) std::cout << ", ";
          }
          std::cout << "]";
      }
      std::cout << ")\n";
      } catch (...) {
      // Skip if the preference cannot be retrieved
      }
   }

}

// reset all preferences:
void preferences_manager::reset_all_preferences_to_defaults() {
   for (const auto& [key, default_value] : preferences_defaults) {
      set_preference(key, default_value); // Reset each preference to its default
   }
}

/* save/load preferences */

// maybe should use another/previous function!?
void preferences_manager::savePreferencesToScript(const std::string& filename) {
   std::ofstream file(filename);
   if (!file.is_open()) {
      throw std::runtime_error("Unable to open file for saving preferences: " + filename);
      }

   int precision = 2;
   double double_value;
   // Write preferences as Python code
   file << "# Auto-generated coot preferences script\n";
   file << "# Modify this file to customize preferences\n\n";
   file << "import coot\n\n";

   for (const auto& [key, callbacks] : preferences_registry) {
      try {
      auto value = get_preference(key);
      file << "coot.set_preference(\"" << key << "\", ";

      // Serialize the value to Python syntax
      if (std::holds_alternative<bool>(value)) {
          file << (std::get<bool>(value) ? "True" : "False");
      } else if (std::holds_alternative<int>(value)) {
          file << std::get<int>(value);
      } else if (std::holds_alternative<double>(value)) {
         double_value = std::get<double>(value);
         double_value < 1. ? precision = 6 : precision = 2;
         file << std::fixed << std::setprecision(precision) << double_value;
      } else if (std::holds_alternative<std::string>(value)) {
          file << "\"" << std::get<std::string>(value) << "\"";
      } else if (std::holds_alternative<std::vector<float>>(value)) {
          const auto& vec = std::get<std::vector<float>>(value);
          file << "[";
          for (size_t i = 0; i < vec.size(); ++i) {
             vec[i] < 1. ? precision = 6 : precision = 2;
             file << std::setprecision(precision) << vec[i];
             if (i < vec.size() - 1) file << ", ";
             }
          file << "]";
      }
      file << ")\n";
      } catch (...) {
      // Skip if the preference cannot be retrieved
      }
   }

   file.close();
}

// probbaly dont need, should use existing function.
void preferences_manager::loadPreferencesFromScript(const std::string& filename) {
   FILE* file = fopen(filename.c_str(), "r");
   if (!file) {
      throw std::runtime_error("Unable to open file for loading preferences: " + filename);
      }

   // Execute the Python script
   PyRun_SimpleFile(file, filename.c_str());
   fclose(file);
}




/* register all preferences */


/* bind all the preference */
// Variables get from existing functions...

void initialize_preferences() {
   graphics_info_t g;
   std::vector<float> def_vec;

   // Mouse rotation button
   coot_preferences.register_preference("use_trackpad",[](const preferences_value& value) {
      set_use_primary_mouse_button_for_rotation(std::get<bool>(value));},
         []() -> preferences_value { graphics_info_t gg; return gg.using_trackpad;},
   g.using_trackpad);
   std::cout<<"BL DEBUG:: done with 1" <<std::endl;

   // Virtual trackball
   coot_preferences.register_preference("virtual_trackball",[](const preferences_value& value) {
      vt_surface(std::get<int>(value));},
         []() -> preferences_value { return vt_surface_status();},
   g.vt_surface_status());
   std::cout<<"BL DEBUG:: done with 2" <<std::endl;

   // Noughty refinement physics
   coot_preferences.register_preference("noughty_refinement_physics",[](const preferences_value& value) {
      set_refine_use_noughties_physics(std::get<int>(value));},
         []() -> preferences_value { return get_refine_use_noughties_physics_state();},
   int(g.noughties_physics));
   std::cout<<"BL DEBUG:: done with 3" <<std::endl;

   // recentre coordinates
   coot_preferences.register_preference("recentre_coordinates",[](const preferences_value& value) {
      set_recentre_on_read_pdb(std::get<int>(value));},
         []() -> preferences_value { return recentre_on_read_pdb();},
   g.recentre_on_read_pdb);
   std::cout<<"BL DEBUG:: done with 4" <<std::endl;

   // Recentre smooth scrolling
   coot_preferences.register_preference("smooth_scroll",[](const preferences_value& value) {
      set_smooth_scroll_flag(std::get<int>(value));},
         []() -> preferences_value { return get_smooth_scroll();},
   g.smooth_scroll);
   std::cout<<"BL DEBUG:: done with 5" <<std::endl;

   // Recentre smooth scrolling steps
   coot_preferences.register_preference("smooth_scroll_steps",[](const preferences_value& value) {
      set_smooth_scroll_steps(std::get<int>(value));},
         []() -> preferences_value { graphics_info_t gg; return gg.smooth_scroll_n_steps;},
   g.smooth_scroll_n_steps);
   std::cout<<"BL DEBUG:: done with 6" <<std::endl;

   // Recentre smooth scrolling limit
   coot_preferences.register_preference("smooth_scroll_limit",[](const preferences_value& value) {
      set_smooth_scroll_limit(std::get<double>(value));},
         []() -> preferences_value { graphics_info_t gg; return gg.smooth_scroll_limit;},
   g.smooth_scroll_limit);
   std::cout<<"BL DEBUG:: done with 7" <<std::endl;


   // ask read state
   coot_preferences.register_preference("map_radius",
                                        [](const preferences_value& value) {
       set_map_radius(std::get<double>(value));},
          []() -> preferences_value { return get_map_radius();},
    g.box_radius_xray);
    def_vec = {g.font_colour.red, g.font_colour.green, g.font_colour.blue};

    coot_preferences.register_preference("font_colour", [](const preferences_value& value) {
       set_font_colour(std::get<std::vector<float>>(value)[0],std::get<std::vector<float>>(value)[1],std::get<std::vector<float>>(value)[2]);},
          []() -> preferences_value { graphics_info_t gg; std::vector<float> ret = {gg.font_colour.red,
                                      gg.font_colour.green,gg.font_colour.blue}; return ret;},
    def_vec);
}

/* FIXME: check ref counting as well, could this just be void!? or we could/should return something
 useful*/
#ifdef USE_PYTHON
PyObject* set_preference(const char *key, PyObject* value) {

   try {
   // bool needs to come first otherwise it will be int...
   if (PyBool_Check(value)) {
      if (PyObject_IsTrue(value)) {
         coot_preferences.set_preference(key, true);
      } else {
         coot_preferences.set_preference(key, false);
      }
   } else if (PyFloat_Check(value)) {
      coot_preferences.set_preference(key, PyFloat_AsDouble(value));
   } else if (PyUnicode_Check(value)) {
      coot_preferences.set_preference(key, PyUnicode_AsUTF8(value));
   } else if (PyLong_Check(value)) {
      coot_preferences.set_preference(key, PyLong_AsLong(value));
   } else if (PyList_Check(value)) {
      // Handle Python list -> std::vector<float>
      std::vector<float> vec;
      for (Py_ssize_t i = 0; i < PyList_Size(value); ++i) {
         PyObject* item = PyList_GetItem(value, i);
         if (!PyFloat_Check(item)) {
            PyErr_SetString(PyExc_ValueError, "List must contain only floats");
            return nullptr;
         }
         vec.push_back(static_cast<float>(PyFloat_AsDouble(item)));
      }
      coot_preferences.set_preference(key, vec);
   } else {
      PyErr_SetString(PyExc_TypeError, "Unsupported preference value type");
      return nullptr;
   }
   } catch (const std::exception& e) {
      PyErr_SetString(PyExc_KeyError, e.what());
      return nullptr;
   }

   Py_RETURN_NONE;
}

PyObject* get_preference(const char *key) {

   try {
   PyObject *ret = Py_None;
   auto value = coot_preferences.get_preference(key);
   if (std::holds_alternative<bool>(value)) {
      if (std::get<bool>(value) ? ret = Py_True : ret = Py_False);
   } else if (std::holds_alternative<double>(value)) {
      ret = PyFloat_FromDouble(std::get<double>(value));
   } else if (std::holds_alternative<std::string>(value)) {
      ret = PyUnicode_FromString(std::get<std::string>(value).c_str());
   } else if (std::holds_alternative<int>(value)) {
      ret = PyLong_FromLong(std::get<int>(value));
   } else if (std::holds_alternative<std::vector<float>>(value)) {
      const auto& vec = std::get<std::vector<float>>(value);
      ret = PyList_New(0);
      for (size_t i = 0; i < vec.size(); ++i) {
         PyList_Append(ret, PyFloat_FromDouble(static_cast<float>(vec[i])));
      }
   };
   return ret;

   } catch (const std::exception& e) {
       PyErr_SetString(PyExc_KeyError, e.what());
       return nullptr;
   }
}


/* reset preferences py */

void reset_all_preferences() {
   try {
   coot_preferences.reset_all_preferences_to_defaults();
   //Py_RETURN_NONE;
   } catch (const std::exception& e) {
      PyErr_SetString(PyExc_RuntimeError, e.what());
      //return nullptr;
   }
}
#endif // PYTHON

/* this is old cruft */

void preferences() {

  show_preferences();
  update_preference_gui();

}

void show_preferences() {

   GtkWidget *w = widget_from_preferences_builder("preferences_dialog");

   graphics_info_t::preferences_widget = w;

   GtkWidget *scrolled_win_model_toolbar = widget_from_preferences_builder("preferences_model_toolbar_icons_scrolledwindow");
   fill_preferences_model_toolbar_icons(w, scrolled_win_model_toolbar);
   GtkWidget *scrolled_win_main_toolbar = widget_from_preferences_builder("preferences_main_toolbar_icons_scrolledwindow");
   fill_preferences_main_toolbar_icons(w, scrolled_win_main_toolbar);

   // we don't want to see the non-General tabs when we first start
   GtkWidget *togglebutton = widget_from_preferences_builder("preferences_general_radiotoolbutton");
   show_hide_preferences_tabs(GTK_TOGGLE_BUTTON(togglebutton), COOT_GENERAL_PREFERENCES);

   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}

void clear_preferences() {

   graphics_info_t::preferences_widget = NULL;

}

void set_mark_cis_peptides_as_bad(int istate) {

   graphics_info_t::mark_cis_peptides_as_bad_flag = istate;
}

int show_mark_cis_peptides_as_bad_state() {

   return graphics_info_t::mark_cis_peptides_as_bad_flag;

}

void show_hide_preferences_tabs(GtkToggleButton *toggletoolbutton, int preference_type) {

   std::vector<std::string> preferences_tabs;

   if (preference_type == COOT_GENERAL_PREFERENCES) {
      preferences_tabs = graphics_info_t::preferences_general_tabs;
   }
   if (preference_type == COOT_BOND_PREFERENCES) {
      preferences_tabs = graphics_info_t::preferences_bond_tabs;
   }
   if (preference_type == COOT_GEOMETRY_PREFERENCES) {
      preferences_tabs = graphics_info_t::preferences_geometry_tabs;
   }
   if (preference_type == COOT_COLOUR_PREFERENCES) {
      preferences_tabs = graphics_info_t::preferences_colour_tabs;
   }
   if (preference_type == COOT_MAP_PREFERENCES) {
      preferences_tabs = graphics_info_t::preferences_map_tabs;
   }
   if (preference_type == COOT_OTHER_PREFERENCES) {
      preferences_tabs = graphics_info_t::preferences_other_tabs;
   }

   // BL says:: shouldnt this be a static vector rather than making it again for each signal?!
   // FIXME BL maybe at some point...
   auto append_tabs = [] (std::vector<std::string> &all_tabs,
                          const std::vector<std::string> &other_tabs) {
      all_tabs.insert(all_tabs.end(), other_tabs.begin(), other_tabs.end());
   };

   std::vector<std::string> all_tabs;
   append_tabs(all_tabs, graphics_info_t::preferences_general_tabs);
   append_tabs(all_tabs, graphics_info_t::preferences_bond_tabs);
   append_tabs(all_tabs, graphics_info_t::preferences_geometry_tabs);
   append_tabs(all_tabs, graphics_info_t::preferences_colour_tabs);
   append_tabs(all_tabs, graphics_info_t::preferences_map_tabs);
   append_tabs(all_tabs, graphics_info_t::preferences_other_tabs);

   GtkWidget *notebook = widget_from_preferences_builder("preferences_notebook");
   int left_tab = 999;
   guint pos = 999;

   for (const auto &tab : all_tabs) {

      GtkWidget *frame = widget_from_preferences_builder(tab);
      if (frame) {
         if (std::find(preferences_tabs.begin(), preferences_tabs.end(), tab) != preferences_tabs.end()) {
            gtk_widget_set_visible(frame, TRUE);
            pos= gtk_notebook_page_num(GTK_NOTEBOOK(notebook), frame);
            // find the tab with the lowest index, i.e. the leftmost
            if (pos < left_tab) {
               left_tab = pos;
               }
         } else {
            gtk_widget_set_visible(frame, FALSE);
         }
      } else {
         // 20240913-PE One day clear this up, don't just remove the error message.
         // speed, antialias bond_parameters and tips have gone.
         // std::cout << "No frame " << preference_type << " " << tab << std::endl;
      }
   }
   // focus on the first tab (alt: save the position, otherwise "random"?!)
   if (left_tab < 999) {
      gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), left_tab);
   }

}

#include "c-interface-preferences.h"

void make_preferences_internal() {

  graphics_info_t g;
  g.make_preferences_internal();
}

void make_preferences_internal_default() {

  make_preferences_internal();
  graphics_info_t g;
  g.preferences_internal_default = g.preferences_internal;

}

void reset_preferences() {

  graphics_info_t g;
  //std::vector<coot::preference_info_t> *ret = g.preferences_internal_default;
  g.preferences_internal = g.preferences_internal_default;
  update_preference_gui();

}


// update and populate Preferences GUI according to (preference file) setting
void update_preference_gui() {

  GtkWidget *w;
  GtkWidget *colour_button = nullptr;
  GtkAdjustment *adjustment;
  GtkWidget *entry;
  std::string text;
  int preference_type;
  int ivalue;
  int ivalue2;
  float fval1;
  float fval2;
  float fval3;
  unsigned short int v = 4;
  graphics_info_t g;

  bool debug = false;

  if (debug)
     std::cout << "--------------------------- update_preference_gui() " << std::endl;

  GtkWidget *dialog = widget_from_preferences_builder("preferences");

  if (debug)
     std::cout << "--------------------------- update_preference_gui() preferences internal size "
               << g.preferences_internal.size() << std::endl;

  for (unsigned int i=0; i<g.preferences_internal.size(); i++) {
     auto preference_type = g.preferences_internal[i].preference_type;
     if (debug)
        std::cout << " -------------- update_preference_gui() "
                  << preference_type << " " << g.preferences_internal[i].ivalue1 << " "
                  << g.preferences_internal[i].fvalue1 << std::endl;
  }

  preferences_value value;
  //case PREFERENCES_VIEW_ROTATION_MOUSE_BUTTON:
  w = widget_from_preferences_builder("preferences_view_rotation_left_mouse_checkbutton");
  value = coot_preferences.get_preference("use_trackpad");
  gtk_check_button_set_active(GTK_CHECK_BUTTON(w), std::get<bool>(value));

  for (unsigned int i=0; i<g.preferences_internal.size(); i++) {
     preference_type = g.preferences_internal[i].preference_type;

     if (debug)
        std::cout << "----------------------------------------- update_preference_gui() "
                  << preference_type << " " << g.preferences_internal[i].ivalue1 << " "
                  << g.preferences_internal[i].fvalue1 << std::endl;

     switch (preference_type) {
      
//     case PREFERENCES_VIEW_ROTATION_MOUSE_BUTTON:
//        w = widget_from_preferences_builder("preferences_view_rotation_left_mouse_checkbutton");
//        ivalue = g.preferences_internal[i].ivalue1;
//        if (ivalue == 1)
//           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
//        else
//           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), FALSE);
//        break;

        case PREFERENCES_VT_SURFACE:
           w = widget_from_preferences_builder("preferences_hid_spherical_radiobutton");
           ivalue = g.preferences_internal[i].ivalue1;
           if (ivalue == 2) {
              gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           } else {
              w = widget_from_preferences_builder("preferences_hid_flat_radiobutton");
              gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           }
           break;

     case PREFERENCES_RECENTRE_PDB:
        w = widget_from_preferences_builder("preferences_recentre_pdb_on_radiobutton");
        if (g.preferences_internal[i].ivalue1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           w = widget_from_preferences_builder("preferences_recentre_pdb_off_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        }
        break;

     // 20240916-PE this has gone
     // case PREFERENCES_BONDS_THICKNESS:
     //    w = widget_from_preferences_builder("preferences_bond_width_combobox");
     //    ivalue = g.preferences_internal[i].ivalue1;
     //    ivalue -= 1;      // offset
     //    gtk_combo_box_set_active(GTK_COMBO_BOX(w), ivalue);
     //    break;

     case PREFERENCES_BOND_COLOURS_MAP_ROTATION:
        w = widget_from_preferences_builder("preferences_bond_colours_hscale");
        fval1 = g.preferences_internal[i].fvalue1;
        adjustment = gtk_range_get_adjustment(GTK_RANGE(w));
        gtk_adjustment_set_value(adjustment, fval1);
        break;

     case PREFERENCES_BOND_COLOUR_ROTATION_C_ONLY:
        w = widget_from_preferences_builder("preferences_bond_colours_checkbutton");
        if (g.preferences_internal[i].ivalue1 == 1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), FALSE);
        }
        break;

     case PREFERENCES_MAP_RADIUS:
        w = widget_from_preferences_builder("preferences_map_radius_entry");
        text = graphics_info_t::float_to_string(g.preferences_internal[i].fvalue1);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_MAP_ISOLEVEL_INCREMENT:
        w = widget_from_preferences_builder("preferences_map_increment_size_entry");
        text = graphics_info_t::float_to_string_using_dec_pl(g.preferences_internal[i].fvalue1, v);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_DIFF_MAP_ISOLEVEL_INCREMENT:
        w = widget_from_preferences_builder("preferences_map_diff_increment_entry");
        text = graphics_info_t::float_to_string_using_dec_pl(g.preferences_internal[i].fvalue1, v);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_MAP_SAMPLING_RATE:
        w = widget_from_preferences_builder("preferences_map_sampling_entry");
        text = graphics_info_t::float_to_string_using_dec_pl(g.preferences_internal[i].fvalue1, v);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_DYNAMIC_MAP_SAMPLING:
        w = widget_from_preferences_builder("preferences_map_dynamic_sampling_checkbutton");
        if (g.preferences_internal[i].ivalue1 == 1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), FALSE);
        }
        break;

     case PREFERENCES_DYNAMIC_MAP_SIZE_DISPLAY:
        w = widget_from_preferences_builder("preferences_map_dynamic_size_checkbutton");
        if (g.preferences_internal[i].ivalue1 == 1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), FALSE);
        }
        break;

     case PREFERENCES_SWAP_DIFF_MAP_COLOURS:
        w = widget_from_preferences_builder("preferences_diff_map_colours_o_radiobutton");
        if (g.preferences_internal[i].ivalue1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           w = widget_from_preferences_builder("preferences_diff_map_colours_coot_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        }
        break;

     case PREFERENCES_MAP_COLOURS_MAP_ROTATION:
        w = widget_from_preferences_builder("preferences_map_colours_hscale");
        fval1 = g.preferences_internal[i].fvalue1;
        adjustment = gtk_range_get_adjustment(GTK_RANGE(w));
        gtk_adjustment_set_value(adjustment, fval1);
        break;

     case PREFERENCES_SMOOTH_SCROLL:
        w = widget_from_preferences_builder("preferences_smooth_scroll_on_radiobutton");
        if (g.preferences_internal[i].ivalue1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           w = widget_from_preferences_builder("preferences_smooth_scroll_off_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        }
        break;

     case PREFERENCES_SMOOTH_SCROLL_STEPS:
        w = widget_from_preferences_builder("preferences_smooth_scroll_steps_entry");
        text = graphics_info_t::int_to_string(g.preferences_internal[i].ivalue1);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_SMOOTH_SCROLL_LIMIT:
        w = widget_from_preferences_builder("preferences_smooth_scroll_limit_entry");
        text = graphics_info_t::float_to_string(g.preferences_internal[i].fvalue1);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_MAP_DRAG:
        w = widget_from_preferences_builder("preferences_map_drag_on_radiobutton");
        if (g.preferences_internal[i].ivalue1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           w = widget_from_preferences_builder("preferences_map_drag_off_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        }
        break;

     // case PREFERENCES_MARK_CIS_BAD:
     //    w = widget_from_preferences_builder("preferences_geometry_cis_peptide_bad_yes_radiobutton");
     //    if (g.preferences_internal[i].ivalue1) {
     //       gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
     //    } else {
     //       w = widget_from_preferences_builder("preferences_geometry_cis_peptide_bad_no_radiobutton");
     //       gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
     //    }
     //    break;

     case PREFERENCES_DEFAULT_B_FACTOR:
        w = widget_from_preferences_builder("preferences_defaults_b_factor_entry");
        {
           std::string s = coot::util::float_to_string(graphics_info_t::default_new_atoms_b_factor);
           gtk_editable_set_text(GTK_EDITABLE(w), s.c_str());
        }
        break;

     case PREFERENCES_BG_COLOUR:

        fval1 = g.preferences_internal[i].fvalue1;  // red
        fval2 = g.preferences_internal[i].fvalue2;  // green
        fval3 = g.preferences_internal[i].fvalue3;  // blue

        GdkRGBA bg_colour;

        if (fval1 < 0.01 && fval2 < 0.01 && fval3 < 0.01) {
           // black
           w = widget_from_preferences_builder("preferences_bg_colour_black_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           bg_colour.red = 0;
           bg_colour.green = 0;
           bg_colour.blue = 0;
        } else if (fval1 > 0.99 && fval2 > 0.99 && fval3 > 0.99) {
           // white
           w = widget_from_preferences_builder("preferences_bg_colour_white_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           bg_colour.red = 65535;
           bg_colour.green = 65535;
           bg_colour.blue = 65535;
        } else {
           // other colour
           w = widget_from_preferences_builder("preferences_bg_colour_own_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           bg_colour.red = (guint)(fval1 * 65535);
           bg_colour.green = (guint)(fval2 * 65535);
           bg_colour.blue = (guint)(fval3 * 65535);
        }
        {
           GtkWidget *colour_button_box = widget_from_preferences_builder("preferences_bg_colour_vbox");
           if (colour_button_box) {
              std::cout << "about to gtk_color_button_set_color() colour_button: " << colour_button
                        << " bg_colour " << bg_colour.red << " " << bg_colour.green << " " << bg_colour.blue << std::endl;

              if (colour_button_box) {
                 GtkWidget *child_item = gtk_widget_get_first_child(colour_button_box);
                 if (child_item) {
                    // the colour button has already been added
                 } else {

                    // c.f. colour button in wrapped_create_show_symmetry_window()

                    // 20230513-PE color dialog is not in GTK 4.4.0 (it is in 4.10+)
#if GTK_MAJOR_VERSION == 5 && GTK_MINOR_VERSION >= 10
                    GtkWidget *col_dialog = gtk_color_dialog_new();
                    // this will need a callback
                    GtkWidget *colour_button_dialog = gtk_color_dialog_button_new(col_dialog);
                    gtk_box_append(GTK_BOX(box_for_colour_button), colour_button_dialog);
#else

                    auto on_color_set_func = +[] (GtkColorButton *self, gpointer user_data) {
                       GdkRGBA rgba;
                       gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(self), &rgba);
                       // std::cout << "Selected color: " << gdk_rgba_to_string(&rgba) << std::endl;
                       float fval1 = static_cast<float>(rgba.red);
                       float fval2 = static_cast<float>(rgba.green);
                       float fval3 = static_cast<float>(rgba.blue);
                       preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, fval1, fval2, fval3);
                       // std::cout << "........  " << fval1 << " " << fval2 << " " << fval3 << std::endl;
                       set_background_colour(fval1, fval2, fval3);
                       graphics_info_t::graphics_draw();
                    };

                    GtkWidget *colour_button = gtk_color_button_new_with_rgba(&bg_colour);
                    gtk_box_append(GTK_BOX(colour_button_box), colour_button);
                    g_signal_connect(G_OBJECT(colour_button), "color-set", G_CALLBACK(on_color_set_func), nullptr);
#endif
                 }
              }
           }
        }
        break;

        // case PREFERENCES_ANTIALIAS:
        //    w = widget_from_preferences_builder("preferences_antialias_on_radiobutton");
        //    if (g.preferences_internal[i].ivalue1) {
        //       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
        //    } else {
        //       w = widget_from_preferences_builder("preferences_antialias_off_radiobutton");
        //       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
        //    }
        //    break;

     case PREFERENCES_CONSOLE_COMMANDS:
        w = widget_from_preferences_builder("preferences_console_info_on_radiobutton");
        if (g.preferences_internal[i].ivalue1) {
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        } else {
           w = widget_from_preferences_builder("preferences_console_info_off_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        }
        break;

     case PREFERENCES_FONT_COLOUR:

        fval1 = g.preferences_internal[i].fvalue1;  // red
        fval2 = g.preferences_internal[i].fvalue2;  // green
        fval3 = g.preferences_internal[i].fvalue3;  // blue
        colour_button = widget_from_preferences_builder("preferences_font_colorbutton");

        GdkRGBA font_colour;

        if (fval1 >= 0.999 &&
            fval2 >= 0.799 && fval2 <= 0.801 &&
            fval3 >= 0.799 && fval3 <= 0.801) {
           // default
           w = widget_from_preferences_builder("preferences_font_colour_default_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           font_colour.red   = (guint)(1.0 * 65535);
           font_colour.green = (guint)(0.8 * 65535);
           font_colour.blue  = (guint)(0.8 * 65535);

        } else {
           // other colour

           w = widget_from_preferences_builder("preferences_font_colour_own_radiobutton");
           gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
           font_colour.red   = (guint)(fval1 * 65535);
           font_colour.green = (guint)(fval2 * 65535);
           font_colour.blue  = (guint)(fval3 * 65535);
        }

        if (colour_button)
           std::cout << "about to gtk_color_button_set_color() colour_button: " << colour_button
                     << " font_colour " << font_colour.red << " " << font_colour.green << " " << font_colour.blue << std::endl;
        else
           std::cout << "about to gtk_color_button_set_color() null colour_button: "
                     << " font_colour " << font_colour.red << " " << font_colour.green << " " << font_colour.blue << std::endl;


        break;

     case PREFERENCES_PINK_POINTER:
        w = widget_from_preferences_builder("preferences_pink_pointer_entry");
        text = graphics_info_t::float_to_string(g.preferences_internal[i].fvalue1);
        gtk_editable_set_text(GTK_EDITABLE(w), text.c_str());
        break;

     case PREFERENCES_PHYSICS:
        w = widget_from_preferences_builder("noughties_physics_on_checkbutton");
        {
           int state = get_refine_use_noughties_physics_state();
           if (state)
              gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
        }
     }
  }
}

#include "utils/xdg-base.hh"

void save_preferences() {

   graphics_info_t g;
   short int istat = 1;
   short int il;
   std::string preferences_name;
   std::filesystem::path full_file_name_path;
   xdg_t xdg;

#ifdef USE_GUILE
   preferences_name = "coot-preferences.scm";
   il = 1;
   full_file_name_path = xdg.get_config_home().append(preferences_name);
   istat = g.save_preference_file(full_file_name_path.string(), il);
   if (istat == 0) {
      std::cout << "WARNING:: failed to save preferences " << full_file_name_path.string() << std::endl;
   }
#endif // USE_GUILE

   preferences_name = "coot_preferences.py";
   full_file_name_path = xdg.get_config_home().append(preferences_name);
   il = 2;
   istat = g.save_preference_file(full_file_name_path.string(), il);
   if (istat == 0) {
      std::cout << "WARNING:: failed to save preferences " << full_file_name_path.string() << std::endl;
   }

}
 

void preferences_internal_change_value_int(int preference_type, int ivalue) {
  graphics_info_t g;
  g.preferences_internal_change_value(preference_type, ivalue);
}

void preferences_internal_change_value_int2(int preference_type, int ivalue1, int ivalue2) {
  graphics_info_t g;
  g.preferences_internal_change_value(preference_type, ivalue1, ivalue2);
}

void preferences_internal_change_value_float(int preference_type, float fvalue) {
  graphics_info_t g;
  g.preferences_internal_change_value(preference_type, fvalue);
}

void preferences_internal_change_value_float3(int preference_type,
                                              float fvalue1, float fvalue2, float fvalue3) {
  graphics_info_t g;
  g.preferences_internal_change_value(preference_type, fvalue1, fvalue2, fvalue3);
}

// FIXME:: make this generic with TOOLBAR enums
void
show_model_toolbar_icon(int pos) {
  graphics_info_t g;
  g.show_hide_toolbar_icon_pos(pos, 1, MODEL_TOOLBAR);
}

void
hide_model_toolbar_icon(int pos) {
  graphics_info_t g;
  g.show_hide_toolbar_icon_pos(pos, 0, MODEL_TOOLBAR);
}

void
show_main_toolbar_icon(int pos) {
  graphics_info_t g;
  g.show_hide_toolbar_icon_pos(pos, 1, MAIN_TOOLBAR);
}

void
hide_main_toolbar_icon(int pos) {
  graphics_info_t g;
  g.show_hide_toolbar_icon_pos(pos, 0, MAIN_TOOLBAR);
}

void
fill_preferences_model_toolbar_icons(GtkWidget *preferences,
				     GtkWidget *scrolled_window) {

  graphics_info_t g;
  g.fill_preferences_model_toolbar_icons(preferences, scrolled_window);
}

void
fill_preferences_main_toolbar_icons(GtkWidget *preferences,
				     GtkWidget *scrolled_window) {

  graphics_info_t g;
  g.fill_preferences_main_toolbar_icons(preferences, scrolled_window);
}


GtkWidget *popup_window(const char *str) {

   // GtkWidget *w = create_popup_info_window();
   GtkWidget *w = widget_from_preferences_builder("popup_info_window");
   GtkWidget *label = widget_from_preferences_builder("info_label");
   gtk_label_set_text(GTK_LABEL(label), str);
   return w;
}

void add_status_bar_text(const std::string &s) {

   graphics_info_t g;
   g.add_status_bar_text(std::string(s));
} 



/*  ----------------------------------------------------------------------- */
/*                  Other interface preferences                            */
/*  ----------------------------------------------------------------------- */

void set_model_fit_refine_dialog_stays_on_top(int istate) { 
   graphics_info_t::model_fit_refine_dialog_stays_on_top_flag = istate;
}

int model_fit_refine_dialog_stays_on_top_state() {

   return graphics_info_t::model_fit_refine_dialog_stays_on_top_flag;
} 

void save_accept_reject_dialog_window_position(GtkWidget *acc_rej_dialog) {
   graphics_info_t g;
   g.save_accept_reject_dialog_window_position(acc_rej_dialog);
}


void set_accept_reject_dialog(GtkWidget *w) { /* used by callbacks to unset the widget.
						 (errr... it wasn't but it is now (as it
						 should be)). */

   graphics_info_t::accept_reject_dialog = w;
}

/* \brief set position of Model/Fit/Refine dialog */
void set_model_fit_refine_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::model_fit_refine_x_position = x_pos;
   graphics_info_t::model_fit_refine_y_position = y_pos;
}

/* \brief set position of Display Control dialog */
void set_display_control_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::display_manager_x_position = x_pos;
   graphics_info_t::display_manager_y_position = y_pos;
}

/* \brief set position of Go To Atom dialog */
void set_go_to_atom_window_position(int x_pos, int y_pos) {

   graphics_info_t::go_to_atom_window_x_position = x_pos;
   graphics_info_t::go_to_atom_window_y_position = y_pos;
}

/* \brief set position of Delete dialog */
void set_delete_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::delete_item_widget_x_position = x_pos;
   graphics_info_t::delete_item_widget_y_position = y_pos;
}

void set_rotate_translate_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::rotate_translate_x_position = x_pos;
   graphics_info_t::rotate_translate_y_position = y_pos;
}

/*! \brief set position of the Accept/Reject dialog */
void set_accept_reject_dialog_position(int x_pos, int y_pos) {
   graphics_info_t::accept_reject_dialog_x_position = x_pos;
   graphics_info_t::accept_reject_dialog_y_position = y_pos;
}

/*! \brief set position of the Ramachadran Plot dialog */
void set_ramachandran_plot_dialog_position(int x_pos, int y_pos) {
   graphics_info_t::ramachandran_plot_x_position = x_pos;
   graphics_info_t::ramachandran_plot_y_position = y_pos;
}

/*! \brief set edit chi angles dialog position */
void set_edit_chi_angles_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::edit_chi_angles_dialog_x_position = x_pos;
   graphics_info_t::edit_chi_angles_dialog_y_position = y_pos;
}


/*! \brief set rotamer selection dialog position */
void set_rotamer_selection_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::rotamer_selection_dialog_x_position = x_pos;
   graphics_info_t::rotamer_selection_dialog_y_position = y_pos;
} 



/*  ------------------------------------------------------------------------ */
/*                     user define clicks                                    */
/*  ------------------------------------------------------------------------ */
#ifdef USE_GUILE
void user_defined_click_scm(int n_clicks, SCM func) {
  if (n_clicks > 0) {
    graphics_info_t g;
    g.user_defined_atom_pick_specs.clear();
    g.in_user_defined_define = n_clicks;
    SCM dest = SCM_BOOL_F;
    SCM mess = scm_from_locale_string("~s");
    SCM v = scm_simple_format(dest, mess, scm_list_1(func));
    std::string func_string = scm_to_locale_string(v);
    g.user_defined_click_scm_func = func;
    g.pick_cursor_maybe();
  } else {
    std::cout<<"INFO:: number of clicks less than 1, cannot define user click"<<std::endl;
  } 
} 
#endif // USE_GUILE

#ifdef USE_PYTHON
void user_defined_click_py(int n_clicks, PyObject *func) {
  if (n_clicks > 0) {
    graphics_info_t g;
    g.user_defined_atom_pick_specs.clear();
    g.in_user_defined_define = n_clicks;
    g.user_defined_click_py_func = func;
    Py_XINCREF(g.user_defined_click_py_func);
    g.pick_cursor_maybe();
  } else {
    std::cout<<"INFO:: number of clicks less than 1, cannot define user click"<<std::endl;
  } 
} 
#endif // USE_PYTHON
/*  ------------------------------------------------------------------------ */
/*                     state (a graphics_info thing)                         */
/*  ------------------------------------------------------------------------ */
void set_save_state_file_name(const char *filename) {
   graphics_info_t::save_state_file_name = filename; 
}

const char *save_state_file_name_raw() {
   return graphics_info_t::save_state_file_name.c_str();
}


#ifdef USE_GUILE
SCM save_state_file_name_scm() {

//    char *f = (char *) malloc(graphics_info_t::save_state_file_name.length() +1);
//    strcpy(f, graphics_info_t::save_state_file_name.c_str());
//    return f;

   std::string f = graphics_info_t::save_state_file_name;
   return scm_from_locale_string(f.c_str());
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *save_state_file_name_py() {
   std::string f = graphics_info_t::save_state_file_name;
   return myPyString_FromString(f.c_str());
}
#endif // USE_PYTHON




char *unmangle_hydrogen_name(const char *pdb_hydrogen_name) {

   std::string atom_name(pdb_hydrogen_name);
   std::string new_atom_name = atom_name;

   if (atom_name.length() == 4) { 
      if (atom_name[0] == '1' ||
	  atom_name[0] == '2' ||
	  atom_name[0] == '3' ||
	  atom_name[0] == '4' ||
	  atom_name[0] == '*') {
	 // switch it.
	 if (atom_name[3] == ' ') {
	    new_atom_name = " "; // lead with a space, testing.
	    new_atom_name += atom_name.substr(1,2) + atom_name[0];
	 } else { 
	    new_atom_name = atom_name.substr(1,3) + atom_name[0];
	 }
      }
   } else { 
      if (atom_name[3] != ' ') { 
	 if (atom_name[3] == ' ') {
	    new_atom_name = atom_name.substr(1,2) + atom_name[0];
	    new_atom_name += ' ';
	 }
	 if (atom_name[2] == ' ') {
	    new_atom_name = atom_name.substr(1,1) + atom_name[0];
	 new_atom_name += ' ';
	 new_atom_name += ' ';
	 }
      } else {
	 // atom_name length is 3 presumably
	 new_atom_name = ' ';
	 new_atom_name += atom_name.substr(1,2) + atom_name[0];
      }
   }

   int new_length = strlen(pdb_hydrogen_name) + 1;
   char *s = new char[new_length];
   for (int i=0; i<new_length; i++) s[i] = 0;
   strncpy(s, new_atom_name.c_str(), new_length);

//    std::cout << "mangle debug:: :" << pdb_hydrogen_name << ": to :" << s << ":" << std::endl;
   return s;
}



short int do_probe_dots_on_rotamers_and_chis_state() {
   return graphics_info_t::do_probe_dots_on_rotamers_and_chis_flag;
} 

void set_do_probe_dots_on_rotamers_and_chis(short int state) {
   graphics_info_t::do_probe_dots_on_rotamers_and_chis_flag = state;
}

void set_do_probe_dots_post_refine(short int state) {
   graphics_info_t::do_probe_dots_post_refine_flag = state;
} 

short int do_probe_dots_post_refine_state() {
   return graphics_info_t::do_probe_dots_post_refine_flag;
}

/* state is 1 for on and 0 for off */
void set_do_coot_probe_dots_during_refine(short int state) {
   graphics_info_t::do_coot_probe_dots_during_refine_flag = state;
}

/* state is 1 for on and 0 for off */
short int get_do_coot_probe_dots_during_refine() {
   return graphics_info_t::do_coot_probe_dots_during_refine_flag;
}



// This is tedious and irritating to parse in C++.
// 
// a const when a member function
// 
// Note that the filenames have a trailing "/".
// 
std::vector<std::pair<std::string, std::string> > 
parse_ccp4i_defs(const std::string &filename) {

   std::vector<std::pair<std::string, std::string> > v;

   // put the current directory in, whether or not we can find the
   // ccp4 project dir
   // on Windows (without mingw or cygwin) there is no PWD,
   // so we set it to "", should work
   char *pwd = getenv("PWD");
   if (pwd) {
      v.push_back(std::pair<std::string, std::string> (std::string(" - Current Dir - "),
						       std::string(pwd) + "/"));
   } else {
#ifdef WINDOWS_MINGW
      v.push_back(std::pair<std::string, std::string> (std::string(" - Current Dir - "),
						       std::string("")));  
#endif // MINGW
   }
   
   struct stat buf;
   int stat_status = stat(filename.c_str(), &buf);
   if (stat_status != 0) {
     // silently return nothing if we can't find the file.
     return v;
   } 

   std::ifstream c_in(filename.c_str());

   // Let's also add ccp4_scratch to the list if the environment
   // variable is declared and if directory exists
   char *scratch = getenv("CCP4_SCR");
   if (scratch) {
      // struct stat buf; no shadow
      // in Windows stat needs to have a last / or \ removed, if existent
#ifdef WINDOWS_MINGW
      if (scratch[strlen(scratch) - 1] == '/') {
	scratch[strlen(scratch) - 1] = '\0';
      }
      if (scratch[strlen(scratch) - 1] == '\\') {
	scratch[strlen(scratch) - 1] = '\0';
      }
#endif // MINGW
      int istat_scratch = stat(scratch, &buf);
      if (istat_scratch == 0) {
	 if (S_ISDIR(buf.st_mode)) {
	    v.push_back(std::pair<std::string, std::string>(std::string("CCP4_SCR"),
							    std::string(scratch) + "/"));
	 }
      }
   }

   if (! c_in) {
      std::cout << "WARNING:: failed to open " << filename << std::endl;
   } else {
      // std::string s;
      char s[1000];
      std::vector <coot::alias_path_t> alias;
      std::vector <coot::alias_path_t> path;
      std::string::size_type ipath;
      std::string::size_type ialias;
      int index = -1;
      int icomma;
      short int path_coming = 0;
      short int alias_coming = 0;
      bool alias_flag = 0;
      while (! c_in.eof()) {
	 c_in >> s;
	 std::string ss(s);
	 // std::cout << "parsing:" << ss << std::endl;
	 if (path_coming == 2) {
	    path.push_back(coot::alias_path_t(index, ss, alias_flag));
	    path_coming = 0;
	 }
	 if (alias_coming == 2) {
	    alias.push_back(coot::alias_path_t(index, ss, alias_flag));
	    alias_coming = 0;
	 }
	 if ( path_coming == 1)  path_coming++;
	 if (alias_coming == 1) alias_coming++;
	 ipath  = ss.find("PROJECT_PATH,");
	 ialias = ss.find("PROJECT_ALIAS,");
	 if (ipath != std::string::npos) {
	    // std::cout << "DEBUG::  found a project path..." << std::endl;
	    path_coming = 1;
	    alias_flag = 0;
	    icomma = ss.find_last_of(",");
	    // std::cout << icomma << " " << ss.length() << std::endl;
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	       // std::cout << "index: " << index << std::endl;
	    }
	 }
	 if (ialias != std::string::npos) {
	    alias_coming = 1;
	    alias_flag = 0;
	    icomma = ss.find_last_of(",");
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	    }
	 }

	 // Things called ALIASES at at the CCP4 top level are
	 // actually speciified by DEF_DIR_PATH and DEF_DIR_ALIAS 
	 // in the same way that PROJECT_ALIAS and PROJECT_PATH work.
	 //

	 ipath  = ss.find("DEF_DIR_PATH,");
	 ialias = ss.find("DEF_DIR_ALIAS,");
	 if (ipath != std::string::npos) {
	    // std::cout << "DEBUG::  found an ALIAS path..." << ss << std::endl;
	    path_coming = 1;
	    alias_flag = 1;
	    icomma = ss.find_last_of(",");
	    // std::cout << icomma << " " << ss.length() << std::endl;
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	       // std::cout << "index: " << index << std::endl;
	    }
	 }
	 if (ialias != std::string::npos) {
	    alias_coming = 1;
	    // std::cout << "DEBUG::  found an ALIAS name..." << ss << std::endl;
	    icomma = ss.find_last_of(",");
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	    }
	 }
	 
      }

//       std::cout << "----------- path pairs: ------------" << std::endl;
//       for (int i=0; i<path.size(); i++) {
//   	 std::cout << path[i].index << "  " << path[i].s << " " << path[i].flag << std::endl;
//       }
//       std::cout << "----------- alias pairs: ------------" << std::endl;
//       for (int i=0; i<alias.size(); i++)
//   	 std::cout << alias[i].index << "  " << alias[i].s << " " << alias[i].flag << std::endl;
//       std::cout << "-------------------------------------" << std::endl;
      

      std::string alias_str;
      std::string path_str;
      for (unsigned int j=0; j<alias.size(); j++) {
	 for (unsigned int i=0; i<path.size(); i++) {
	    if (path[i].index == alias[j].index) {
	       if (path[i].flag == alias[j].flag) {
		  // check for "" "" pair here.
		  alias_str = alias[j].s;
		  path_str  = path[i].s;
		  // if the file is a directory, we need to put a "/" at the
		  // end so that went we set that filename in the fileselection
		  // widget, we go into the directory, rather than being in the
		  // directory above with the tail as the selected file.
		  //
		  struct stat buf_l;
		  int status = stat(path_str.c_str(), &buf_l);
	       
		  // valgrind says that buf.st_mode is uninitialised here
		  // strangely.  Perhaps we should first test for status?
		  // Yes - that was it.  I was using S_ISDIR() on a file
		  // that didn't exist.  Now we skip if the file does not
		  // exist or is not a directory.

		  // std::cout << "stating "<< path_str << std::endl;

		  if (status == 0) { 
		     if (S_ISDIR(buf_l.st_mode)) {
			path_str += "/";

			if (alias_str == "\"\"") {
			   alias_str = "";
			   path_str  = "";
			}
			v.push_back(std::pair<std::string, std::string> (alias_str, path_str));
		     }
		     // } else { 
		     // // This is too boring to see every time we open a file selection
		     // std::cout << "INFO:: directory for a CCP4i project: " 
		     // << path_str << " was not found\n";
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

std::string
ccp4_project_directory(const std::string &ccp4_project_name) {

   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > v = 
      parse_ccp4i_defs(ccp4_defs_file_name);
   std::string r = "";
   for (unsigned int i=0; i<v.size(); i++) {
      if (v[i].first == ccp4_project_name) {
	 r = v[i].second;
	 break;
      }
   }
   return r;
}


/* movies */
void set_movie_file_name_prefix(const char *file_name) {
   graphics_info_t::movie_file_prefix = file_name;
}

void set_movie_frame_number(int frame_number) {
   graphics_info_t::movie_frame_number = frame_number;
}

#ifdef USE_GUILE
SCM movie_file_name_prefix() {
   SCM r = scm_from_locale_string(graphics_info_t::movie_file_prefix.c_str());
   return r;
}
#endif
#ifdef USE_PYTHON
PyObject * movie_file_name_prefix_py() {
   PyObject *r;
   r = myPyString_FromString(graphics_info_t::movie_file_prefix.c_str());
   return r;
}
#endif // PYTHON

int movie_frame_number() {
   return graphics_info_t::movie_frame_number;
}

void set_make_movie_mode(int make_movie_flag) {
   graphics_info_t::make_movie_flag = make_movie_flag;
}


#ifdef USE_GUILE
void try_load_scheme_extras_dir() {

   // Function no longer used? Check and delete.

   char *s = getenv("COOT_SCHEME_EXTRAS_DIR");
   if (s) {

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
      std::vector<std::string> dirs = coot::util::split_string(s, ";");
#else
      std::vector<std::string> dirs = coot::util::split_string(s, ":");
#endif
      for (unsigned int i=0; i<dirs.size(); i++) { 
	 struct stat buf;
	 int status = stat(dirs[i].c_str(), &buf);
	 if (status != 0) {
	    std::cout << "WARNING:: no directory \"" << dirs[i] << "\""
		      << " in COOT_SCHEME_EXTRAS_DIR " << s
		      << std::endl;
	 } else {
	    if (S_ISDIR(buf.st_mode)) {

	       DIR *lib_dir = opendir(dirs[i].c_str());
	       if (lib_dir == NULL) {
		  std::cout << "An ERROR occured on opening the directory "
			    << dirs[i] << std::endl;
	       } else {

		  struct dirent *dir_ent;

		  // loop until the end of the filelist (readdir returns NULL)
		  // 
		  while (1) {
		     dir_ent = readdir(lib_dir);
		     if (dir_ent == NULL) {
			break;
		     } else {
			std::string sub_part(std::string(dir_ent->d_name));
			struct stat buf2;
			std::string fp = s;
			fp += "/";
			fp += sub_part;
			int status2 = stat(fp.c_str(), &buf2);
			if (status2 != 0) {
			   // std::cout << "WARNING:: no file " << sub_part << " in directory "
			   // << dirs[i] << std::endl;
			} else {
			   if (S_ISREG(buf2.st_mode)) {
			      if (coot::util::file_name_extension(sub_part) == ".scm") {
				 std::cout << "loading extra: " << fp << std::endl;
				 scm_c_primitive_load(fp.c_str());
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
#endif // USE_GUILE


#ifdef USE_PYTHON
void try_load_python_extras_dir() {

   char *s = getenv("COOT_PYTHON_EXTRAS_DIR");
   if (s) {
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
      std::vector<std::string> dirs = coot::util::split_string(s, ";");
#else
      std::vector<std::string> dirs = coot::util::split_string(s, ":");
#endif
      for (unsigned int i=0; i<dirs.size(); i++) {
	 struct stat buf;
	 int status = stat(dirs[i].c_str(), &buf);
	 if (status != 0) {
	    std::cout << "WARNING:: no directory \"" << dirs[i] << "\""
		      << " in COOT_PYTHON_EXTRAS_DIR " << s
		      << std::endl;
	 } else {
	    if (S_ISDIR(buf.st_mode)) {

	       DIR *lib_dir = opendir(dirs[i].c_str());
	       if (lib_dir == NULL) {
		  std::cout << "An ERROR occured on opening the directory "
			    << dirs[i] << std::endl;
	       } else {

		  struct dirent *dir_ent;

		  // loop until the end of the filelist (readdir returns NULL)
		  // 
		  while (1) {
		     dir_ent = readdir(lib_dir);
		     if (dir_ent == NULL) {
			break;
		     } else {
			std::string sub_part(std::string(dir_ent->d_name));
			struct stat buf2;
			std::string fp = s;
			fp += "/";
			fp += sub_part;
			int status2 = stat(fp.c_str(), &buf2);
			if (status2 != 0) {
			   // std::cout << "WARNING:: no file " << sub_part << std::endl;
			} else {
			   if (S_ISREG(buf2.st_mode)) {
			      if (coot::util::file_name_extension(sub_part) == ".py") {
				 std::cout << "loading python extra: " << fp << std::endl;
				 run_python_script(fp.c_str()); 
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
#endif // USE_PYTHON


void set_button_label_for_external_refinement(const char *button_label) {
   graphics_info_t::external_refinement_program_button_label = button_label;
}


int preferences_internal_font_own_colour_flag() {

   int r = -1; 
   graphics_info_t g;
   for (unsigned int i=0; i<g.preferences_internal.size(); i++) {
      if (g.preferences_internal[i].preference_type == PREFERENCES_FONT_OWN_COLOUR_FLAG) {
	 r = g.preferences_internal[i].ivalue1;
	 break;
      }
   }

   return r;
} 
