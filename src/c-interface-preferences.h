/* src/c-interface-preferences.h
 * 
 * Copyright 2011 by the University of Oxford
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifndef C_INTERFACE_PREFERENCES_H
#define C_INTERFACE_PREFERENCES_H

#include <gtk/gtk.h>

#include <functional>
//#include <unordered_map>
#include <map>
#include <variant>
#include <iostream>
#include <stdexcept>
#include <fstream>


// New stuff
using preferences_value = std::variant<int, double, std::string, bool, std::vector<float>>;

class preferences_manager {

   struct preferences_callbacks {
       std::function<void(const preferences_value&)> preference_set_function;
       std::function<preferences_value()> preference_get_function;
   };

   // can be unordered_map too... this will then be random
   // using map it will be alphabetically sorted by key...
   std::map<std::string, preferences_callbacks> preferences_registry;
   std::map<std::string, preferences_value> preferences_defaults;
   std::map<std::string, preferences_value> preferences_values;

public:
   void register_preference(
       const std::string& key,
       const std::function<void(const preferences_value&)>& set_function,
       const std::function<preferences_value()>& get_function,
       const preferences_value& default_value);
   void set_preference(const std::string& key, const preferences_value& value);
   preferences_value get_preference(const std::string& key) const;
   void reset_preference_to_default(const std::string& key);
   preferences_value get_preference_default(const std::string& key) const;
   void list_preferences() const;
   void reset_all_preferences_to_defaults();
   // these may be replaced by old/existing functions.
   int save_preferences_to_file(const std::string& filename);
   void load_preferences_from_file(const std::string& filename);
};

extern preferences_manager coot_preferences;
void initialize_preferences();

#ifdef __cplusplus
#ifdef USE_PYTHON
// python binding for preferences
PyObject* set_preference(const char *key, PyObject* args);
PyObject* get_preference(const char *key);
void reset_all_preferences();
#endif // PYTHON
#endif

// This is old and for the gtk interface

#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

/*  ----------------------------------------------------------------------- */
/*                  Preferences Notebook                                    */
/*  ----------------------------------------------------------------------- */
/* section Preferences */
void preferences();
void show_preferences();
void setup_preferences_gui();
void clear_preferences();
void set_mark_cis_peptides_as_bad(int istate); /* in geometry graph */
int show_mark_cis_peptides_as_bad_state();
#ifndef SWIG
#if (GTK_MAJOR_VERSION >= 4)
void show_hide_preferences_tabs(GtkToggleButton *toggletoolbutton, int preference_type);
#else
void show_hide_preferences_tabs(GtkToggleToolButton *toggletoolbutton, int preference_type);
#endif
#endif

void update_preference_gui();
void make_preferences_internal();
void make_preferences_internal_default();
void reset_preferences();
void save_preferences();
void preferences_internal_change_value_int(int preference_type, int ivalue);
void preferences_internal_change_value_int2(int preference_type, int ivalue1, int ivalue2);
void preferences_internal_change_value_float(int preference_type, float fvalue);
void preferences_internal_change_value_float3(int preference_type, 
					float fvalue1, float fvalue2, float fvalue3);
void show_model_toolbar_icon(int pos);
void hide_model_toolbar_icon(int pos);

void show_main_toolbar_icon(int pos);
void hide_main_toolbar_icon(int pos);

int preferences_internal_font_own_colour_flag();

#ifndef SWIG
void fill_preferences_model_toolbar_icons(GtkWidget *preferences,
				     	  GtkWidget *scrolled_window);
void fill_preferences_main_toolbar_icons(GtkWidget *preferences,
				     	  GtkWidget *scrolled_window);
#endif

END_C_DECLS

#endif /* C_INTERFACE_PREFERENCES_H */

