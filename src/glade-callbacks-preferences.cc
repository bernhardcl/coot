/* src/glade-callbacks-preferences.cc
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2008 The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
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


#include "Python.h"

#include <iostream>
#include <gtk/gtk.h>

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"
#include "positioned-widgets.h"
#include "interface.h"
#include "coot-references.h"

// put preferences functions into their own file, not here.
#include "coot-preferences.h"
#include "c-interface-preferences.h"
#include "rotate-translate-modes.hh"
#include "restraints-editor-c.h"
#include "generic-display-objects-c.h"
#include "c-interface-refmac.h"
#include "gtk-widget-conversion-utils.h"
#include "curlew.h"
#include "read-phs.h"
#include "gtk-manual.h"
#include "c-interface-refine.h"
#include "utils/coot-utils.hh"

#include "widget-from-builder.hh"

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;


extern "C" G_MODULE_EXPORT
void
on_preferences_general_radiotoolbutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                         gpointer         user_data) {
  show_hide_preferences_tabs(toggletoolbutton, COOT_GENERAL_PREFERENCES);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bond_radiotoolbutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_BOND_PREFERENCES);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_radiotoolbutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_MAP_PREFERENCES);
}



extern "C" G_MODULE_EXPORT
void
on_preferences_geometry_radiotoolbutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_GEOMETRY_PREFERENCES);
}




extern "C" G_MODULE_EXPORT
void
on_preferences_colour_radiotoolbutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_COLOUR_PREFERENCES);
}




extern "C" G_MODULE_EXPORT
void
on_preferences_other_radiotoolbutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_OTHER_PREFERENCES);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *w = lookup_widget(GTK_WIDGET(button), "preferences");
  GtkWidget *w = widget_from_preferences_builder("preferences_dialog");
  save_preferences();
  gtk_widget_set_visible(w, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_preferences_reset_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_preferences_builder("preferences_dialog");
  reset_all_preferences();
  update_preference_gui();
  // hide or not after reset?
  //gtk_widget_set_visible(w, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_preferences_destroy                 (GtkWidget       *object,
                                        gpointer         user_data)
{
  clear_preferences();
}

void set_use_trackpad(short int state); // or #include cc-interface.hh

extern "C" G_MODULE_EXPORT
void
on_preferences_view_rotation_left_mouse_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                            gpointer         user_data) {
   coot_preferences.set_preference("use_trackpad",
                                   static_cast<bool>(gtk_check_button_get_active(checkbutton)));
}


extern "C" G_MODULE_EXPORT
void
on_preferences_hid_spherical_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                 gpointer         user_data) {

   coot_preferences.set_preference("virtual_trackball",
                                   static_cast<int>(gtk_check_button_get_active(checkbutton)) + 1);

}


extern "C" G_MODULE_EXPORT
void
on_preferences_bond_colours_hscale_value_changed(GtkRange        *range,
                                                 gpointer         user_data) {

   GtkAdjustment *adjustment;
   float fvalue;
   adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
   fvalue = gtk_adjustment_get_value(adjustment);
   coot_preferences.set_preference("bond_colour_map_rotation", fvalue);
}

extern "C" G_MODULE_EXPORT
void
on_preferences_bond_colours_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                gpointer         user_data) {

   coot_preferences.set_preference("bond_colour_map_rotation_c_only",
                                   static_cast<int>(gtk_check_button_get_active(checkbutton)));

}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_black_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {

   std::vector<float> bg_colour(3, 0.);
   if (gtk_check_button_get_active(checkbutton)) {
      coot_preferences.set_preference("background_colour", bg_colour);
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_nearlyblack_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {

   std::vector<float> bg_colour(3, 0.035);
   if (gtk_check_button_get_active(checkbutton)) {
      coot_preferences.set_preference("background_colour", bg_colour);
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_darkgrey_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {
   std::vector<float> bg_colour(3, 0.07);
   if (gtk_check_button_get_active(checkbutton)) {
      coot_preferences.set_preference("background_colour", bg_colour);
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_semidarkgrey_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {
   std::vector<float> bg_colour(3, 0.207);
   if (gtk_check_button_get_active(checkbutton)) {
      coot_preferences.set_preference("background_colour", bg_colour);
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_lightgrey_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {
   std::vector<float> bg_colour(3, 0.83);
   if (gtk_check_button_get_active(checkbutton)) {
      coot_preferences.set_preference("background_colour", bg_colour);
   }
}

extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_white_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {
   std::vector<float> bg_colour = {1., 1., 1.};
   if (gtk_check_button_get_active(checkbutton)) {
      coot_preferences.set_preference("background_colour", bg_colour);
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_own_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                 gpointer         user_data)
{

   // BL says;: only do something if active
   if (gtk_check_button_get_active(checkbutton)) {
      const GdkRGBA *bg_colour;
      // 20220528-PE more here
      GtkWidget *w = widget_from_preferences_builder("preferences_background_color_button");
      bg_colour = gtk_color_dialog_button_get_rgba(GTK_COLOR_DIALOG_BUTTON(w));
      std::vector<float> bg_col = {
                                   (float) bg_colour->red,
                                   (float) bg_colour->green,
                                   (float) bg_colour->blue};
      coot_preferences.set_preference("background_colour", bg_col);
   }

}



// Generic function to deal with changing entries
// requires a struct to pass multiple arguments to a timeout function
typedef struct {
   GtkWidget *widget;
   std::string preference_key;
   guint timeout_id;
   bool is_int = false;
} EntryTimeoutData;

// a map to keep track of the timeout_ids
static std::unordered_map<std::string, guint> timeout_id_map;

extern "C" G_MODULE_EXPORT
void
on_preferences_map_radius_entry_activate(GtkEntry        *entry,
                                         gpointer         user_data) {

   // not strictly needed but probably good to have the activate here to ensure
   // instant application of the map radius change.
   // Note: not done for other entries which have changed signals since actually
   // redundant

   // when we handle activate explicitly we may remove the corresponding timeout
   std::string key = "map_radius";
   // Remove existing timeout if present
   auto it = timeout_id_map.find(key);
   if (it != timeout_id_map.end()) {
      g_source_remove(it->second);
      timeout_id_map.erase(it);
   }

   std::cout << "debug:: on_preferences_map_radius_entry_activate() entry " << entry << std::endl;
   const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   float fval = 0;
   fval = atof(text);
   if ((fval > 0) && (fval <1000)) {
      coot_preferences.set_preference("map_radius", fval);
   }

}

gboolean
process_preferences_entry_text(gpointer user_data) {

   if (!user_data) return G_SOURCE_REMOVE;

   auto* data = static_cast<EntryTimeoutData*>(user_data);
   const gchar* text = gtk_editable_get_text(GTK_EDITABLE(data->widget));
   float fval = 0.;
   fval = atof(text);
   if ((fval > 0.) && (fval <1999.9)) {
      std::string key = data->preference_key;
      if (data->is_int) {
         // have int entry
         int ival = atoi(text);
         coot_preferences.set_preference(key, ival);
      } else {
         // have float
         coot_preferences.set_preference(key, fval);
      }
   }

   timeout_id_map.erase(data->preference_key); // Remove timeout ID entry
   delete data; // clean up
   return G_SOURCE_REMOVE; // Remove timeout
}

// use a timeout function as not to process the entry change whilst editing (at least not too much)
extern "C" G_MODULE_EXPORT
void
on_preferences_map_radius_entry_changed(GtkEditable     *editable,
                                        gpointer         user_data) {

   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_map_radius_entry"));
   if (entry) {
      std::string key = "map_radius";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

      // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
    }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_increment_size_entry_changed(GtkEditable     *editable,
                                                gpointer         user_data) {

   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_map_increment_size_entry"));
   if (entry) {
      std::string key = "map_iso_level_increment";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

      // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }
}


/* BL says:: FIXME again, not used so remove!?
extern "C" G_MODULE_EXPORT
void
on_preferences_map_diff_increment_entry_activate(GtkEntry        *entry,
                                                 gpointer         user_data) {

   // not used any more FIXME
   const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   float fval = 0;
   fval = atof(text);
   if (fval > 0) {
      coot_preferences.set_preference("diff_map_iso_level_increment", fval);
   }

}
*/


extern "C" G_MODULE_EXPORT
void
on_preferences_map_diff_increment_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_map_diff_increment_entry"));
   if (entry) {
      std::string key = "diff_map_iso_level_increment";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

              // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }

}


// BL says:: maybe not needed but good for instance response...
extern "C" G_MODULE_EXPORT
void
on_preferences_map_sampling_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  std::string key = "map_sampling_rate";
  // Remove existing timeout if present
  auto it = timeout_id_map.find(key);
  if (it != timeout_id_map.end()) {
     g_source_remove(it->second);
     timeout_id_map.erase(it);
  }
  const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  float fval = 0;
  fval = atof(text);
  if ((fval < 100) && (fval > 1)) {
    coot_preferences.set_preference("map_sampling_rate", fval);
  }

}



extern "C" G_MODULE_EXPORT
void
on_preferences_map_sampling_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_map_sampling_entry"));
   if (entry) {
      std::string key = "map_sampling_rate";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

      // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_dynamic_sampling_checkbutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data)
{
    coot_preferences.set_preference("dynamic_map_sampling",
                                    static_cast<int>(gtk_check_button_get_active(checkbutton)));

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_dynamic_size_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                    gpointer         user_data) {
      coot_preferences.set_preference("dynamic_map_display_size",
                                      static_cast<int>(gtk_check_button_get_active(checkbutton)));
}


extern "C" G_MODULE_EXPORT
void
on_preferences_diff_map_colours_coot_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                         gpointer        user_data) {
   coot_preferences.set_preference("swap_diff_map_colours",
                                   static_cast<int>(!gtk_check_button_get_active(checkbutton)));
}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_colours_hscale_value_changed(GtkRange        *range,
                                                gpointer         user_data) {
   GtkAdjustment *adjustment;
   float fvalue;
   adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
   fvalue = gtk_adjustment_get_value(adjustment);
   coot_preferences.set_preference("map_colour_map_rotation", fvalue);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_on_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                    gpointer        user_data) {

   coot_preferences.set_preference("smooth_scroll",
                                   static_cast<int>(gtk_check_button_get_active(checkbutton)));

}


// BL note:: at this point there is only "changed" signal for smooth_scroll_step
// dont need activate!?
extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_steps_entry_changed(GtkEditable     *editable,
                                                 gpointer         user_data) {

   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_smooth_scroll_steps_entry"));
   if (entry) {
      std::string key = "smooth_scroll_steps";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

      // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0, true};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_limit_entry_activate(GtkEntry        *entry,
                                                  gpointer         user_data) {

   // BL says:: currently not used - changed signal should be enough!?

}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_limit_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_smooth_scroll_limit_entry"));
   if (entry) {
      std::string key = "smooth_scroll_limit";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

              // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_drag_on_radiobutton_toggled(GtkCheckButton *checkbutton,
                                               gpointer         user_data) {

   coot_preferences.set_preference("map_drag",
                                   static_cast<int>(gtk_check_button_get_active(checkbutton)));

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_drag_off_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                gpointer         user_data) {

   // not used any more FIXME
   if (gtk_check_button_get_active(checkbutton)) {
      //preferences_internal_change_value_int(PREFERENCES_MAP_DRAG, 0);
      set_active_map_drag_flag(0);
   }
}


extern "C" G_MODULE_EXPORT
void on_preferences_default_b_factor_entry_activate(GtkEntry        *entry,
                                                    gpointer         user_data) {

   /* BL says:: FIXME maybe not needed any more since replaced by changed signal
   const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      float f = coot::util::string_to_float(std::string(text));
      coot_preferences.set_preference("default_b_factor", f);
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING:: in on_preferences_default_b_factor_entry_activate(): " << e.what() << std::endl;
   }
*/
}

extern "C" G_MODULE_EXPORT
void on_preferences_default_b_factor_entry_changed(GtkEditable     *editable,
                                                   gpointer         user_data)
{
   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_default_b_factor_entry"));
   if (entry) {
      std::string key = "default_b_factor";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

      // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }

}



extern "C" G_MODULE_EXPORT
void
on_preferences_recentre_pdb_on_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {

   coot_preferences.set_preference("recentre_coordinates",
                                   static_cast<int>(gtk_check_button_get_active(checkbutton)));

}



extern "C" G_MODULE_EXPORT
void
on_preferences_console_info_on_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {

   coot_preferences.set_preference("console_display_commands",
                                   static_cast<bool>(gtk_check_button_get_active(checkbutton)));

}


extern "C" G_MODULE_EXPORT
void
on_preferences_console_info_off_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                    gpointer         user_data) {

   // just do this in the above radio button callback - it's 2-state
   // FIXME not used any more/currently

   if (gtk_check_button_get_active(checkbutton)) {
      //preferences_internal_change_value_int(PREFERENCES_CONSOLE_COMMANDS, 0);
      set_console_display_commands_state(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_font_colour_default_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                       gpointer         user_data) {

   if (gtk_check_button_get_active(checkbutton)) {
      // get from default value
      preferences_value font_colour = coot_preferences.get_preference_default("font_colour");
      coot_preferences.set_preference("font_colour", font_colour);
   }
}

extern "C" G_MODULE_EXPORT
void
on_preferences_font_colour_own_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                   gpointer         user_data) {

   // 20230716-PE There is no function to set the colour of the colorbutton if
   // the button is in the .ui file (which it is in this case (we look up the
   // widget from the name)).

   if (gtk_check_button_get_active(checkbutton)) {
      const GdkRGBA *font_colour;
      GtkWidget *colorbutton = widget_from_preferences_builder("preferences_font_color_button");
      font_colour = gtk_color_dialog_button_get_rgba(GTK_COLOR_DIALOG_BUTTON(colorbutton));
      std::vector<float> font_col = {
                                   (float) font_colour->red,
                                   (float) font_colour->green,
                                   (float) font_colour->blue};
      coot_preferences.set_preference("font_colour", font_col);
   }

}


// BL says:: not sure if we need the activate signal but probably good for instance change
extern "C" G_MODULE_EXPORT
void
on_preferences_pink_pointer_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  // when we handle activate explicitly we may remove the corresponding timeout
  std::string key = "rotation_centre_cube_size";
  // Remove existing timeout if present
  auto it = timeout_id_map.find(key);
  if (it != timeout_id_map.end()) {
     g_source_remove(it->second);
     timeout_id_map.erase(it);
  }

  float fval;
  const gchar* entry_text = gtk_editable_get_text(GTK_EDITABLE(entry));
  try {
     fval = coot::util::string_to_float(std::string(entry_text));
     if ((fval > 1000) || (fval < 0)) {
        printf("Invalid cube size: %s Assuming default 0.1 A \n", entry_text);
        fval  = 0.1;
     }
     coot_preferences.set_preference("rotation_centre_cube_size", fval);
  }
  catch (const std::runtime_error &e) {
     std::cout << "WARNING:: in on_preferences_pink_pointer_entry_activate(): " << e.what() << std::endl;
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_pink_pointer_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
   GtkEntry *entry = GTK_ENTRY(widget_from_preferences_builder("preferences_pink_pointer_entry"));
   if (entry) {
      std::string key = "rotation_centre_cube_size";
      // Remove existing timeout if present
      auto it = timeout_id_map.find(key);
      if (it != timeout_id_map.end()) {
         g_source_remove(it->second);
         timeout_id_map.erase(it);
      }

              // Set a timeout and pass the entry as user_data
      auto* data = new EntryTimeoutData{GTK_WIDGET(editable), key, 0};
      guint id = g_timeout_add(400, process_preferences_entry_text, data);
      data->timeout_id = id; // store the ids
      timeout_id_map[key] = id;
   }
}


extern "C" G_MODULE_EXPORT
void
noughties_physics_switch_state_set(GtkSwitch *switch_widget,
                                   gboolean   state,
                                   gpointer   user_data) {

   if (state)
      set_show_unit_cells_all(1);
  else
     set_show_unit_cells_all(0);
}

extern "C" G_MODULE_EXPORT
void
on_noughties_physics_checkbutton_toggled(GtkCheckButton *toggletoolbutton,
                                         gpointer         user_data) {

   coot_preferences.set_preference("noughty_refinement_physics",
                                   static_cast<int>(gtk_check_button_get_active(toggletoolbutton)));

}

extern "C" G_MODULE_EXPORT
void on_preferences_background_color_selected(GtkColorDialogButton *button,
                                              GParamSpec *pspec, gpointer user_data) {

   const GdkRGBA *color = gtk_color_dialog_button_get_rgba(button);
   std::vector<float> bg_col = {
      (float) color->red,
      (float) color->green,
      (float) color->blue};

   coot_preferences.set_preference("background_colour", bg_col);

}

extern "C" G_MODULE_EXPORT
void on_preferences_font_color_selected(GtkColorDialogButton *button,
                                        GParamSpec *pspec, gpointer user_data) {

   const GdkRGBA *color = gtk_color_dialog_button_get_rgba(button);
   std::vector<float> f_col = {
                               (float) color->red,
                               (float) color->green,
                               (float) color->blue};

   coot_preferences.set_preference("font_colour", f_col);

}

