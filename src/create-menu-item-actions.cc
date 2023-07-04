
#include <iostream>
#include <gtk/gtk.h>

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"
#include "widget-from-builder.hh"
#include "c-interface-gui.hh" // set transient for main window
#include "cc-interface-scripting.hh"  // move this up
#include "cc-interface.hh"
#include "fit-loop-gui.hh"
#include "read-molecule.hh" // 20230621-PE now with std::string args

#include "support.h" // for internationalizations.
#include "c-interface-preferences.h"
#include "c-interface-ligands-swig.hh"
#include "curlew-gtk4.hh"

void fill_and_show_shader_preferences(); // should this be in a header?
                                         // it's in glade-callbacks.cc (that doesn't seem right today)

extern "C" { void load_tutorial_model_and_data(); }

extern "C" G_MODULE_EXPORT
void on_coords_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file   = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);

#if 0
      GSList *files_list = gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(dialog));
      while (files_list) {

         const char *fnc = static_cast<const char *>(files_list->data);
         if (fnc) {
            std::string fn(fnc);
            handle_read_draw_molecule_with_recentre(fn, 0);
         }
         files_list = g_slist_next(files_list);
      }
#endif

      bool move_molecule_here_flag = false;
      bool recentre_on_read_pdb_flag = true; // was false;

      const char *r = gtk_file_chooser_get_choice(GTK_FILE_CHOOSER(dialog), "recentering");
      if (r) {
         std::string sr(r);
         if (sr == "No Recentre")
            recentre_on_read_pdb_flag = false;
         if (sr == "Move Molecule Here")
            move_molecule_here_flag = true;
      }

      if (file_name) {
         std::cout << "INFO: " << file_name << " move_molecule_here_flag: " << move_molecule_here_flag
                   << " recentre_on_read_pdb_flag " << recentre_on_read_pdb_flag << std::endl;
         if (move_molecule_here_flag) {
            handle_read_draw_molecule_and_move_molecule_here(file_name);
         } else {
            if (recentre_on_read_pdb_flag)
               handle_read_draw_molecule_with_recentre(file_name, 1);
            else
               handle_read_draw_molecule_with_recentre(file_name, 0); // no recentre
         }
      }
   }
   gtk_window_close(GTK_WINDOW(dialog));
}


void on_dataset_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                 int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);
      bool auto_read_flag = false, ismtz = false, ismtzauto = false, iscnsauto = false;
      auto_read_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(chooser), "auto_read_flag"));
      ismtz = is_mtz_file_p(file_name);

      if (ismtz)
         ismtzauto = mtz_file_has_phases_p(file_name);
      else
         iscnsauto = cns_file_has_phases_p(file_name);

      if ( ismtzauto || iscnsauto ) {
         if (auto_read_flag)
            wrapped_auto_read_make_and_draw_maps(file_name);
         else
            manage_column_selector(file_name); /* calls create_column_label_window(), fills and displays.*/
      } else {
         manage_column_selector(file_name); // try read a cif (strangely)
      }

   }
   gtk_window_close(GTK_WINDOW(dialog));
}


void on_map_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                             int response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
      GFile *file = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);
      int is_diff_map_flag = 0; // needs fixing obviously... FIXME
      handle_read_ccp4_map(file_name, is_diff_map_flag);
   }
   gtk_window_close(GTK_WINDOW(dialog));
}

void open_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   // Ancient GTK3
   // open_coords_dialog();

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File",
                                                   parent_window,
                                                   action,
                                                   _("_Cancel"),
                                                  GTK_RESPONSE_CANCEL,
                                                   _("_Open"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);

   // void gtk_file_chooser_add_choice (GtkFileChooser* chooser,
   //                                   const char* id,
   //                                   const char* label,
   //                                   const char** options,
   //                                   const char** option_labels)


   // gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "recentre",      "Centre on New Molecule", NULL, NULL);
   // gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "no-recentre",   "No Recentre",            NULL, NULL);
   // gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "move-mol-here", "Move Molecule Here",     NULL, NULL);

   const gchar *labels[]  = {"Centre on New Molecule", "No Recentre",      "Move Molecule Here", NULL};
   const gchar *options[] = {"recentre-view",          "no-recentre-view", "move-mol-here",      NULL};

   // I don't follow the options and labels, but this works strangely.
   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "recentering", "Recentre", options, labels);

   add_filename_filter_button(dialog, COOT_COORDS_FILE_SELECTION);

   g_signal_connect(dialog, "response", G_CALLBACK(on_coords_filechooser_dialog_response_gtk4), NULL);
   gtk_widget_set_visible(dialog, TRUE);

}

void open_dataset_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {
#if 0
   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   GtkWidget *main_window = graphics_info_t::get_main_window();
   gtk_window_set_transient_for(GTK_WINDOW(dataset_chooser), GTK_WINDOW(main_window));

   set_directory_for_filechooser(dataset_chooser);
   set_file_selection_dialog_size(dataset_chooser);
   add_filechooser_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_set_visible(dataset_chooser, TRUE);
#endif

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", parent_window, action,
                                                   _("_Cancel"), GTK_RESPONSE_CANCEL,
                                                   _("_Open"), GTK_RESPONSE_ACCEPT,
                                                   NULL);
   g_signal_connect(dialog, "response", G_CALLBACK(on_dataset_filechooser_dialog_response_gtk4), NULL);
   g_object_set_data(G_OBJECT(dialog), "auto_read_flag", GINT_TO_POINTER(FALSE));
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.mtz");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_set_visible(dialog, TRUE);
}

void auto_open_mtz_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {

#if 0
   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   int is_auto_read_fileselection = 1;
   set_directory_for_filechooser(dataset_chooser);
   add_filename_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   g_object_set_data(G_OBJECT(dataset_chooser), "imol", GINT_TO_POINTER(-1)); // 20220627-PE do I need this?
   g_object_set_data(G_OBJECT(dataset_chooser), "is_auto", GINT_TO_POINTER(is_auto_read_fileselection));
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_set_visible(dataset_chooser, TRUE);
#endif

   // How were is the user-data set?
   GtkWindow *parent_window = GTK_WINDOW(user_data);
   if (user_data) {
      parent_window = GTK_WINDOW(user_data);
   }
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", parent_window, action,
                                                   _("_Cancel"), GTK_RESPONSE_CANCEL,
                                                   _("_Open"), GTK_RESPONSE_ACCEPT,
                                                   NULL);
   g_signal_connect(dialog, "response", G_CALLBACK(on_dataset_filechooser_dialog_response_gtk4), NULL);
   g_object_set_data(G_OBJECT(dialog), "auto_read_flag", GINT_TO_POINTER(TRUE));
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.mtz");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_set_visible(dialog, TRUE);
}

void open_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {

#if 0 // still problems with the the Builder FileChooser dialogs
   GtkWidget *dataset_chooser = widget_from_builder("map_name_filechooser_dialog");
   set_directory_for_filechooser(dataset_chooser);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_set_visible(dataset_chooser, TRUE);
#endif


   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File",
                                                   parent_window,
                                                   action,
                                                   _("_Cancel"),
                                                   GTK_RESPONSE_CANCEL,
                                                   _("_Open"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);

   g_signal_connect(dialog, "response", G_CALLBACK(on_map_filechooser_dialog_response_gtk4), NULL);

   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.map");
   gtk_file_filter_add_pattern(filterselect, "*.mrc");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_set_visible(dialog, TRUE);
}


void load_tutorial_model_and_data_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                         G_GNUC_UNUSED GVariant *parameter,
                                         G_GNUC_UNUSED gpointer user_data) {
   load_tutorial_model_and_data();
}

void exit_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   coot_checked_exit(0);
}

void calculate_hydrogen_bonds_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                     G_GNUC_UNUSED GVariant *parameter,
                                     G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      calculate_hydrogen_bonds(imol);
   }
}


void curlew_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                   G_GNUC_UNUSED GVariant *parameter,
                   G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *dialog = curlew_dialog();
   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

}

void get_monomer_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *frame    = widget_from_builder("get_monomer_frame");
   GtkWidget *no_entry_frame = widget_from_builder("get_monomer_no_entry_frame");
   if (no_entry_frame)
      gtk_widget_set_visible(no_entry_frame, FALSE); // each time "get_monomer" is shown

   GtkWidget *entry = widget_from_builder("get_monomer_entry");
   gtk_widget_grab_focus(entry);
   gtk_widget_set_visible(frame, TRUE);
}


void on_cif_dictionary_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                        int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);
      int imol_enc = -999997;

      bool create_ligand = false;
      const char *r = gtk_file_chooser_get_choice(GTK_FILE_CHOOSER(dialog), "create-molecule");
      if (r) {
         std::string sr(r);
         if (sr == "Create New Instance")
            create_ligand = true;
      }

      if (create_ligand) {
         std::cout << "create ligand! " << std::endl;
      }
      int monomer_index = handle_cif_dictionary_for_molecule(file_name, imol_enc, create_ligand);
   }
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

void import_cif_dictionary_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                  G_GNUC_UNUSED GVariant *parameter,
                                  G_GNUC_UNUSED gpointer user_data) {

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", parent_window, action,
                                                   _("_Cancel"), GTK_RESPONSE_CANCEL,
                                                   _("_Open"), GTK_RESPONSE_ACCEPT,
                                                   NULL);

   const gchar *labels[]  = {"No Instance", "Create New Instance", NULL};
   const gchar *options[] = {"no-instance", "create-new-instance", NULL};
   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "create-molecule", "Create Molecule", options, labels);
   g_signal_connect(dialog, "response", G_CALLBACK(on_cif_dictionary_filechooser_dialog_response_gtk4), NULL);
   g_object_set_data(G_OBJECT(dialog), "auto_read_flag", GINT_TO_POINTER(FALSE));
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.cif");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

}

void toggle_display_frames_per_second_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                             G_GNUC_UNUSED GVariant *parameter,
                                             G_GNUC_UNUSED gpointer user_data) {

int state = get_fps_flag();
if (state == 1)
   set_show_fps(0);
else
   set_show_fps(1);
}

void
search_monomer_library_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = widget_from_builder("monomer_search_dialog");
   gtk_widget_set_visible(w, TRUE);
}

void show_accession_code_fetch_frame(G_GNUC_UNUSED GSimpleAction *simple_action,
                                     G_GNUC_UNUSED GVariant *parameter,
                                     G_GNUC_UNUSED gpointer user_data) {
   gchar* mode_name_cstr;
   g_variant_get(parameter,"s",&mode_name_cstr);
   std::string mode_name(mode_name_cstr);
   auto mode_num_from_name = [](const std::string& mode_name){
      if(mode_name == "oca") {
         return COOT_ACCESSION_CODE_WINDOW_OCA;
      } else if(mode_name == "eds") {
         return COOT_ACCESSION_CODE_WINDOW_EDS;
      } else if (mode_name == "pdb-redo"){
         return COOT_ACCESSION_CODE_WINDOW_PDB_REDO;
      } else if(mode_name == "uniprot-id") {
         return COOT_UNIPROT_ID;
      } else {
         g_error("Unrecognized mode name for the accession code frame: %s",mode_name.c_str());
         // fallback
         return COOT_ACCESSION_CODE_WINDOW_OCA;
      }
   };
   int mode_num = mode_num_from_name(mode_name);
   g_debug("Accession code fetch frame mode number: %i",mode_num);
   GtkWidget *frame = widget_from_builder("accession_code_frame");
   g_object_set_data(G_OBJECT(frame), "mode", GINT_TO_POINTER(mode_num));
   GtkWidget *label = widget_from_builder("accession_code_label");
   switch(mode_num) {
      case COOT_ACCESSION_CODE_WINDOW_EDS:
      case COOT_ACCESSION_CODE_WINDOW_OCA:
      case COOT_ACCESSION_CODE_WINDOW_PDB_REDO:
      {
         gtk_label_set_text(GTK_LABEL(label), "PDB Accession Code: ");
         break;
      }
      case COOT_UNIPROT_ID:{
         gtk_label_set_text(GTK_LABEL(label), "UniProt ID: ");
         break;
      }
      default: {
         g_error("Unrecognized mode number for the accession code frame: %i",mode_num);
         // fallback
         gtk_label_set_text(GTK_LABEL(label), "PDB Accession Code: ");
         break;
      }
   }

   GtkWidget* entry = widget_from_builder("accession_code_entry");
   gtk_widget_grab_focus(entry);
   // this is probably equivalent
   //gtk_widget_set_visible(frame,TRUE);
   gtk_widget_set_visible(frame, TRUE);
}

void
fetch_and_superpose_alphafold_models_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                     G_GNUC_UNUSED GVariant *parameter,
                                     G_GNUC_UNUSED gpointer user_data) {
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      fetch_and_superpose_alphafold_models(imol);
   }
}

void
fetch_map_from_emdb_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *frame = widget_from_builder("emdb_map_code_frame");
   gtk_widget_set_visible(frame, TRUE);
   
}


void
fetch_pdbe_ligand_description_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                     G_GNUC_UNUSED GVariant *parameter,
                                     G_GNUC_UNUSED gpointer user_data) {
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::residue_spec_t res_spec(pp.second.second);
      const auto &m = graphics_info_t::molecules[imol];
      std::string comp_id = m.get_residue_name(res_spec);
      // python-function: coot_utils.get_SMILES_for_comp_id_from_pdbe arg: comp_id
      std::cout << "run python function coot_utils.get_SMILES_for_comp_id_from_pdbe " << comp_id << std::endl;
      short int lang = coot::STATE_PYTHON;
      std::vector<coot::command_arg_t> args = { coot::command_arg_t(comp_id) };
      std::string sc = g.state_command("coot_utils", "get_SMILES_for_comp_id_from_pdbe", args, lang);
      std::cout << ":::::::::::::::::::::: python command: " << sc << std::endl;
      safe_python_command("import coot_utils"); // Hack. This has already happened, but python has forgotten.
      safe_python_command(sc);
   }
}

void
save_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   GCallback callback_func = G_CALLBACK(save_molecule_coords_combobox_changed);
   int imol = first_coords_imol();
   int imol_unsaved = first_unsaved_coords_imol();
   if (imol_unsaved != -1)
      imol = imol_unsaved;
   std::cout << "DEBUG:: in on_save_coordinates1_activate() with imol_unsaved "
             << imol_unsaved << std::endl;
   set_save_molecule_number(imol); /* set *save* molecule number */

   // this is the molecule chooser, not the file chooser
   //
   GtkWidget *widget = widget_from_builder("save_coords_dialog");
   GtkWidget *combobox = widget_from_builder("save_coordinates_combobox");

   if (combobox) {
      fill_combobox_with_coordinates_options(combobox, callback_func, imol);
      set_transient_and_position(COOT_UNDEFINED_WINDOW, widget);
      gtk_widget_set_visible(widget, TRUE);
      gtk_window_present(GTK_WINDOW(widget));
   } else {
      std::cout << "ERROR:: in on_save_coordinates1_activate() bad combobox!\n";
   }
}

void
save_symmetry_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {

   setup_save_symmetry_coords();

}

void
on_save_state_dialog_response(GtkDialog *dialog,
                              int response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      char *file_name = g_file_get_path(file);
      std::cout << "Now save state script to file " << file_name << std::endl;
      short int il = coot::PYTHON_SCRIPT;
      graphics_info_t g;
      g.save_state_file(file_name, il);
   }

   // maybe save the dialog in graphics_info_t, and just hide it?
   // gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

   gtk_window_destroy(GTK_WINDOW(dialog));
}

void
save_state_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   GtkWindow *parent_window = GTK_WINDOW(graphics_info_t::get_main_window());
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Save State",
                                                   parent_window,
                                                   action,
                                                   _("Cancel"),
                                                   GTK_RESPONSE_CANCEL,
                                                   _("Save"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.py");
   gtk_file_filter_add_pattern(filterselect, "*.scm");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_set_size_request(dialog, 800, 700);
   g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_save_state_dialog_response), NULL);
   gtk_widget_set_visible(dialog, TRUE);
}


void
on_save_views_dialog_response(GtkDialog *dialog,
                              int response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      char *file_name = g_file_get_path(file);
      std::cout << "Now save views script to file " << file_name << std::endl;
      save_views(file_name);
   }
   gtk_window_destroy(GTK_WINDOW(dialog));
}

void
on_save_views_clicked(GtkButton *button, gpointer user_data) {

   GtkWindow *parent_window = GTK_WINDOW(graphics_info_t::get_main_window());
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Save Views",
                                                   parent_window,
                                                   action,
                                                   _("Cancel"),
                                                   GTK_RESPONSE_CANCEL,
                                                   _("Save"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.py");
   gtk_file_filter_add_pattern(filterselect, "*.scm");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_set_size_request(dialog, 800, 700);
   g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_save_views_dialog_response), NULL);
   gtk_widget_set_visible(dialog, TRUE);
}

void
recover_session_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   recover_session();
}

void
file_export_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   short int is_fragment = false;
   export_map_gui(is_fragment);

}

void
file_export_map_fragment_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {
   short int is_fragment = true;
   export_map_gui(is_fragment);

}

void
close_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *widget = wrapped_create_new_close_molecules_dialog(); // uses builder
   set_transient_for_main_window(widget);
   gtk_widget_set_visible(widget, TRUE);
}

void
change_chain_ids_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_change_chain_id_dialog(); // uses builder
   gtk_widget_set_visible(w, TRUE);
}

void
make_link_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {

   std::cout << "make_link_action(): coot user_defined click 2" << std::endl;
   // 20220828-PE needs check_if_in_range_define
}

#include "cc-interface.hh" // for fullscreen()

void
fix_nomenclature_errors_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      fix_nomenclature_errors(imol);
   }
}


void
invert_this_chiral_centre_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      const coot::atom_spec_t &atom_spec = pp.second.second;
      invert_chiral_centre(imol, atom_spec.chain_id, atom_spec.res_no, atom_spec.ins_code, atom_spec.atom_name);
   }
}

void
merge_solvent_chains_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::cout << "merge solvent chains for imol " << imol << std::endl;
   }
}

void
copy_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
   do_edit_copy_molecule();
}

void
copy_molecule_fragment_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
   do_edit_copy_fragment();
}

void
merge_molecules_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}

void
move_molecule_here_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = widget_from_builder("move_molecule_here_dialog");
   fill_move_molecule_here_dialog(w);
   gtk_widget_set_visible(w, TRUE);
}

void
mutate_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_mutate_sequence_dialog();
   gtk_widget_set_visible(w, TRUE);
}

void
edit_replace_fragment_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {
     do_edit_replace_fragment();
}

void
edit_replace_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {
     do_edit_replace_residue();
}

void
renumber_residues_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_renumber_residue_range_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}

void
renumber_waters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      renumber_waters(imol);
   }
}

void
undo_molecule_chooser_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {
   show_set_undo_molecule_chooser();
}

void
residue_info_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {
   // do_residue_info_dialog(); // this waits for a click - the old way
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::residue_spec_t res_spec(pp.second.second);
      graphics_info_t g;
      g.output_residue_info_dialog(imol, res_spec);
   }
}

void
edit_restraints_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {

   // fills residue types in the combobox
   GtkWidget *w =  wrapped_create_residue_editor_select_monomer_type_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}

void
exchange_chain_ids_for_seg_ids_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                      G_GNUC_UNUSED GVariant *parameter,
                                      G_GNUC_UNUSED gpointer user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      exchange_chain_ids_for_seg_ids(imol);
   }
}


void
show_preferences_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   show_preferences();
   update_preference_gui();
}


void
show_shader_preferences_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {
   fill_and_show_shader_preferences();
}


void
align_and_mutate_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_align_and_mutate_dialog();
   gtk_widget_set_visible(w, TRUE);
}

void
fit_loop_by_database_search(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

  wrapped_fit_loop_db_loop_dialog();

}

void
fit_loop_by_ramachandran_search(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = create_fit_loop_rama_search_dialog();
   gtk_widget_set_visible(w, TRUE);

}


void
ligand_builder_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {

#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "FIXME:: start_ligand_builder_gui_internal() " << std::endl;
#else
   start_ligand_builder_gui_internal(menuitem, user_data);
#endif
}


void
on_run_script_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                               int response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER(dialog);
      GFile *file = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);

      std::cout << "Run this script file: " << file_name << std::endl;
      run_script(file_name);
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

   }
}

void
run_script_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   // reading from the ui file means that there are no buttons
   // std::cout << "lllllllllllllllllllllllll run script action" << std::endl;
   // GtkWidget *dialog = widget_from_builder("run_script_filechooser_dialog");
   // add_filename_filter_button(dialog, COOT_SCRIPTS_FILE_SELECTION);
   // gtk_widget_set_visible(dialog, TRUE);

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Run Script File",
                                                   parent_window,
                                                   action,
                                                   _("_Cancel"),
                                                   GTK_RESPONSE_CANCEL,
                                                   _("_Open"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);

   g_signal_connect(dialog, "response", G_CALLBACK(on_run_script_filechooser_dialog_response_gtk4), NULL);
   add_filename_filter_button(dialog, COOT_SCRIPTS_FILE_SELECTION);
   gtk_widget_set_visible(dialog, TRUE);

}

void
scripting_python_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   reveal_python_scripting_entry();
}

void
scripting_scheme_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   std::cout << "launch the scheme dialog here" << std::endl;

}

void
use_clustalw_for_alignment_then_mutate_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                              G_GNUC_UNUSED GVariant *parameter,
                                              G_GNUC_UNUSED gpointer user_data) {

   std::cout << "launch a python gui for clustalw" << std::endl;

   // generic_chooser_entry_and_file_selector("Align Sequence to Model: ",
   //                                         coot_utils.valid_model_molecule_qm,
   //                                         "Chain ID",
   //                                         "",
   //                                         "Select PIR Alignment file",
   //                                         lambda imol, chain_id, target_sequence_pif_file:
   //                                         coot.run_clustalw_alignment(imol, chain_id, target_sequence_pif_file)))

}


void
mask_map_by_atom_selection_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                  G_GNUC_UNUSED GVariant *parameter,
                                  G_GNUC_UNUSED gpointer user_data) {

   // 20230427-PE surprisingly similar to wrapped_create_unmodelled_blobs_dialog()

   graphics_info_t g;
   GtkWidget *dialog = widget_from_builder("mask_map_by_atom_selection_dialog");
   GtkWidget *model_combobox = widget_from_builder("mask_map_by_atom_selection_model_combobox");
   GtkWidget *map_combobox   = widget_from_builder("mask_map_by_atom_selection_map_combobox");
   
   int imol_mol_active = -1;
   int imol_map_active = -1;
   GCallback func = G_CALLBACK(nullptr); // we don't care until this dialog is read

   auto get_model_molecule_vector = [] () {
                                       graphics_info_t g;
                                       std::vector<int> vec;
                                       int n_mol = g.n_molecules();
                                       for (int i=0; i<n_mol; i++)
                                          if (g.is_valid_model_molecule(i))
                                             vec.push_back(i);
                                       return vec;
                                    };

   auto get_map_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_map_molecule(i))
                                           vec.push_back(i);
                                     return vec;
                                  };

   auto model_list = get_model_molecule_vector();
   auto   map_list = get_map_molecule_vector();
   if (! model_list.empty()) imol_mol_active = model_list[0];
   if (!   map_list.empty()) imol_map_active =   map_list[0];

   g.fill_combobox_with_molecule_options(model_combobox, func, imol_mol_active, model_list);
   g.fill_combobox_with_molecule_options(  map_combobox, func, imol_map_active,   map_list);

   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

}

void
copy_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                G_GNUC_UNUSED GVariant *parameter,
                G_GNUC_UNUSED gpointer user_data) {

   // not modern but works
   std::string cmd = "import coot; import coot_gui; coot_gui.map_molecule_chooser_gui(\"Molecule to Copy...\", lambda imol: coot.copy_molecule(imol))";
   safe_python_command(cmd);

}

void
make_a_smoother_copy_of_a_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                     G_GNUC_UNUSED GVariant *parameter,
                                     G_GNUC_UNUSED gpointer user_data) {

   safe_python_command("import coot_gui");
   std::string sc = "coot_gui.map_molecule_chooser_gui(\"Map Molecule to Smoothenize...\", lambda imol: coot.smooth_map(imol, 1.25))";
   safe_python_command(sc);

}

void
make_a_very_smooth_copy_of_a_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

   safe_python_command("import coot_gui");
   std::string sc = "coot_gui.map_molecule_chooser_gui(\"Map Molecule to Smoothenize...\", lambda imol: coot.smooth_map(imol, 2.0))";
   safe_python_command(sc);
}

void
make_a_difference_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {
}

void
transform_map_by_lsq_model_fit_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                      G_GNUC_UNUSED GVariant *parameter,
                                      G_GNUC_UNUSED gpointer user_data) {
}

void
average_maps_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {
}

void
brighten_maps_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
}

void
set_map_is_difference_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {
}

void
another_contour_level_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {
}

void
multichicken_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {
}

void
copy_ncs_residue_range_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
}

void
copy_ncs_chain_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {
}

void
ncs_jumping_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                   G_GNUC_UNUSED GVariant *parameter,
                   G_GNUC_UNUSED gpointer user_data) {
}

void
ncs_ligands_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                   G_GNUC_UNUSED GVariant *parameter,
                   G_GNUC_UNUSED gpointer user_data) {
}

void add_HOLE_module_action(GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   short int lang = coot::STATE_PYTHON;
   std::vector<coot::command_arg_t> args = {};
   std::string sc = g.state_command("coot_hole", "hole_ify", args, lang);
   safe_python_command("import coot_gui");
   safe_python_command("import coot_hole");
   std::cout << "calling this: " << sc << std::endl;
   safe_python_command(sc);

   // needed?
   // g_simple_action_set_enabled(simple_action,FALSE);
}

void add_ccp4_module_action(GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {
   safe_python_command("import coot_gui");
   safe_python_command("coot_gui.add_module_ccp4()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_carbohydrate_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                    G_GNUC_UNUSED GVariant *parameter,
                                    G_GNUC_UNUSED gpointer user_data) {
   safe_python_command("import gui_add_linked_cho");
   safe_python_command("gui_add_linked_cho.add_module_carbohydrate_gui()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_cryo_em_module_action(GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   safe_python_command("import coot_gui");
   safe_python_command("coot_gui.add_module_cryo_em()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_ligand_module_action(GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   safe_python_command("import coot_utils");
   safe_python_command("import coot_gui");
   safe_python_command("coot_gui.add_module_ligand()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_prosmart_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {
   safe_python_command("import gui_prosmart");
   safe_python_command("gui_prosmart.add_module_prosmart()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_rcrane_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
   std::cout << "INFO:: no RCrane" << std::endl;
   info_dialog("INFO:: No RCrane interface yet");
}

void add_restraints_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                  G_GNUC_UNUSED GVariant *parameter,
                                  G_GNUC_UNUSED gpointer user_data) {
   safe_python_command("import gui_prosmart");
   safe_python_command("gui_prosmart.add_module_restraints()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_refine_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
   safe_python_command("import coot_gui");
   safe_python_command("coot_gui.add_module_refine()");
   g_simple_action_set_enabled(simple_action,FALSE);
}

void add_shelx_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   std::cout << "No SHELXL yet! " << std::endl;
   info_dialog("INFO:: No SHELXL interface yet! - sorry");
}

void add_views_module_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   // Add a button in the toolbar where we can add a ligand popupdialog
   GtkWidget *toolbar_hbox = graphics_info_t::get_widget_from_builder("main_window_toolbar_hbox"); // toolbar style hbox

   GtkApplication *app = graphics_info_t::application;
   GMenuModel *menubar = gtk_application_get_menubar(app);

   // This is how the menubar is used in python:

   // menu_refine = Gtk.Menu()
   // menuitem = Gtk.MenuItem(s)
   // menuitem.set_submenu(menu_refine)
   // main_menubar = coot_gui_api.main_menubar()
   // main_menubar.append(menuitem)
   //
   // then use it:
   // coot_gui.add_simple_coot_menu_menuitem(menu_refine, "All-atom Refine", lambda arg: all_atom_refine_func())

   GtkWidget *view_menubutton = gtk_menu_button_new();
   gtk_menu_button_set_label(GTK_MENU_BUTTON(view_menubutton), "Views");
   GtkWidget *popover = gtk_popover_new();
   gtk_popover_set_position(GTK_POPOVER(popover), GTK_POS_BOTTOM);
   gtk_menu_button_set_popover(GTK_MENU_BUTTON(view_menubutton), popover);

   // Create the content for the popover
   GtkWidget *outer_box  = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 6);
   GtkWidget *views_hbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
   GtkWidget *content_box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
   gtk_widget_set_margin_start(content_box, 6);
   gtk_widget_set_margin_end(content_box, 6);
   gtk_widget_set_margin_top(content_box, 6);
   gtk_widget_set_margin_bottom(content_box, 6);
   gtk_popover_set_child(GTK_POPOVER(popover), outer_box);

   GtkWidget *add_view_button   = gtk_button_new_with_label("Add View");
   GtkWidget *save_views_button = gtk_button_new_with_label("Save View");
   GtkWidget *play_views_button = gtk_button_new_with_label("Play Views");

   auto add_view_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {

      auto view_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {
         int view_idx = GPOINTER_TO_INT(user_data);
         int snap_mode = 0;
         go_to_view_number(view_idx, snap_mode);
      };

      std::string label = "A View";
      int view_idx = add_view_here(label.c_str());
      GtkWidget *views_box = GTK_WIDGET(user_data);
      GtkWidget *new_button = gtk_button_new_with_label(label.c_str());
      g_signal_connect(G_OBJECT(new_button), "clicked", G_CALLBACK(view_button_clicked_callback), GINT_TO_POINTER(view_idx));
      gtk_box_append(GTK_BOX(views_box), new_button);
   };

   auto save_views_button_clicked_callback = +[] (G_GNUC_UNUSED GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      on_save_views_clicked(button, user_data); // arguments are not needed.
   };

   auto play_views_button_clicked_callback = +[] (G_GNUC_UNUSED GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      play_views();
   };

   g_signal_connect(G_OBJECT(  add_view_button), "clicked", G_CALLBACK(add_view_button_clicked_callback),   views_hbox);
   g_signal_connect(G_OBJECT(save_views_button), "clicked", G_CALLBACK(save_views_button_clicked_callback), nullptr);
   g_signal_connect(G_OBJECT(play_views_button), "clicked", G_CALLBACK(play_views_button_clicked_callback), nullptr);

   GtkWidget *h_sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
   gtk_box_append(GTK_BOX(content_box), add_view_button);
   gtk_box_append(GTK_BOX(content_box), play_views_button);
   gtk_box_append(GTK_BOX(content_box), h_sep);
   gtk_box_append(GTK_BOX(content_box), save_views_button);
   gtk_box_append(GTK_BOX(outer_box), content_box);
   gtk_box_append(GTK_BOX(outer_box), views_hbox);

   gtk_box_append(GTK_BOX(toolbar_hbox), view_menubutton);
 
}



void
associate_sequence_file_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {
}

void
assign_sequence_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
}

void
lsq_superpose_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_least_squares_dialog(); // uses builder
   gtk_widget_set_visible(w, TRUE);
}

// ---------------- where is LSQ Plane? --------------------------


void
other_modelling_tools_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_other_model_tools_dialog();
   gtk_widget_set_visible(w, TRUE);
}

void
whats_this_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      graphics_info_t g;
      auto &m = g.molecules[imol];
      coot::residue_spec_t rs(pp.second.second);
      std::string rn = m.get_residue_name(rs);
      std::string s = rs.format() + std::string(" ") + rn;
      add_status_bar_text(s);
   }
}





void
ssm_superposition_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_superpose_dialog(); // uses builder

   /* we get returned w = 0 when there is no MMDBSSM. (We are doing it
      this way because we don't have to introduce HAVE_MMDBSSM into the
      *c* compiler arguments (this is simpler)).  */
  if (w)
     gtk_widget_set_visible(w, TRUE);
}


void
sharpen_blur_for_xray_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   std::cout << "sharpen/blur for x-ray" << std::endl;

#if 0
   GtkWidget *dialog   = widget_from_builder("map_sharpening_dialog");
   GtkWidget *combobox = widget_from_builder("map_sharpening_molecule_combobox");
   GtkWidget *scale    = widget_from_builder("map_sharpening_hscale");

   auto get_map_molecule_vector = [] () {
      graphics_info_t g;
      std::vector<int> vec;
      int n_mol = g.n_molecules();
      for (int i=0; i<n_mol; i++)
         if (g.is_valid_map_molecule(i))
            vec.push_back(i);
      return vec;
   };

   graphics_info_t g;
   int imol_active = -1;
   GCallback func = G_CALLBACK(nullptr); // we don't care until this dialog is read
   auto model_list = get_map_molecule_vector();
   g.fill_combobox_with_molecule_options(combobox, func, imol_active, model_list);

   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

#endif

   GtkWidget *dialog = wrapped_create_map_sharpening_dialog();
   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

}

void
calculate_updating_maps_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   // show_calculate_updating_maps_pythonic_gui();

   auto get_model_molecule_vector = [] () {
                                       graphics_info_t g;
                                       std::vector<int> vec;
                                       int n_mol = g.n_molecules();
                                       for (int i=0; i<n_mol; i++)
                                          if (g.is_valid_model_molecule(i))
                                             vec.push_back(i);
                                       return vec;
                                    };

   auto get_map_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_map_molecule(i))
                                           vec.push_back(i);
                                     return vec;
                                  };

   auto get_diff_map_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_map_molecule(i))
                                           if (g.is_difference_map(i))
                                              vec.push_back(i);
                                     return vec;
                                  };

   graphics_info_t g;
   GtkWidget *dialog            = widget_from_builder("updating_maps_dialog");
   GtkWidget *model_combobox    = widget_from_builder("updating_maps_model_combobox");
   GtkWidget *map_combobox      = widget_from_builder("updating_maps_map_combobox");
   GtkWidget *diff_map_combobox = widget_from_builder("updating_maps_diff_map_combobox");

   int imol_mol_active = -1;
   int imol_map_active = -1;

   auto    model_list =    get_model_molecule_vector();
   auto      map_list =      get_map_molecule_vector();
   auto diff_map_list = get_diff_map_molecule_vector();

   std::cout << "::::::::::::::::::::::: diff_map_list size " << diff_map_list.size() << std::endl;

   GCallback func = G_CALLBACK(nullptr); // we don't care until this dialog is read

   g.fill_combobox_with_molecule_options(   model_combobox, func, imol_mol_active,    model_list);
   g.fill_combobox_with_molecule_options(     map_combobox, func, imol_map_active,      map_list);
   g.fill_combobox_with_molecule_options(diff_map_combobox, func, imol_map_active, diff_map_list);

   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);
}



void
background_black_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t::background_colour = glm::vec3(0,0,0);
   graphics_info_t::graphics_draw();
   
}

void
background_dark_grey_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t::background_colour = glm::vec3(0.07f,0.07f,0.07f);
   graphics_info_t::graphics_draw();
}

void
background_light_grey_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t::background_colour = glm::vec3(0.83f, 0.83f, 0.83f);
   graphics_info_t::graphics_draw();
}

void
background_white_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t::background_colour = glm::vec3(1,1,1);
   graphics_info_t::graphics_draw();
}

void
bond_colours_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = widget_from_builder("coords_colour_control_dialog");
   graphics_info_t g;
   g.fill_bond_colours_dialog_internal(w);
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);

}


void
bond_parameters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_bond_parameters_dialog(); // uses builder
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}

void add_hydrogen_atoms_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot_add_hydrogen_atoms(imol);
   }
}

void add_hydrogen_atoms_using_refmac_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                            G_GNUC_UNUSED GVariant *parameter,
                                            G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      graphics_info_t g;
      int imol = pp.second.first;
      short int lang = coot::STATE_PYTHON;
      std::vector<coot::command_arg_t> args = { coot::command_arg_t(imol) };
      std::string sc = g.state_command("coot_utils", "add_hydrogens_using_refmac", args, lang);
      safe_python_command("import coot_utils");
      safe_python_command(sc);
   }
}

void add_an_atom_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   // GtkWidget *dialog = widget_from_builder("add_an_atom_dialog");
   // set_transient_for_main_window(dialog);
   // gtk_widget_set_visible(dialog, TRUE);

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   gtk_widget_set_visible(box, TRUE);

}

void add_other_solvent_molecules_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

   safe_python_command("import coot_gui");
   safe_python_command("coot_gui.solvent_ligands_gui()");

}


// this should be in a header, I suppose.
GtkWidget *wrapped_create_find_waters_dialog();

void find_waters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_find_waters_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);

}

void dna_rna_models_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = widget_from_builder("nucleotide_builder_dialog");
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}


void place_helix_here_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   place_helix_here();
}



void find_ligands_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {

   do_find_ligands_dialog();

}

void cis_trans_convert_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *at =  aa.second;
      std::string atom_name = at->name;
      cis_trans_convert(imol, at->GetChainID(), at->GetSeqNum(), at->GetInsCode());
   }
}

void add_OXT_to_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_add_OXT_dialog(); // uses builder
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);

}

void reverse_chain_direction_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                    G_GNUC_UNUSED GVariant *parameter,
                                    G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      if (is_valid_model_molecule(imol)) {
         // 20230520-PE pass a std::string here.
	 reverse_direction_of_fragment(imol, pp.second.second.chain_id.c_str(), pp.second.second.res_no);
      }
   }

}

void arrange_waters_around_protein_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                          G_GNUC_UNUSED GVariant *parameter,
                                          G_GNUC_UNUSED gpointer user_data) {

}

void assign_hetatms_for_this_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                            G_GNUC_UNUSED GVariant *parameter,
                                            G_GNUC_UNUSED gpointer user_data) {

}

void assign_hetatoms_to_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

}

void backrub_rotamers_for_chain_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                       G_GNUC_UNUSED GVariant *parameter,
                                       G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      graphics_info_t g;
      int imol = pp.second.first;
      std::string chain_id = pp.second.second.chain_id;
      short int lang = coot::STATE_PYTHON;
      std::vector<coot::command_arg_t> args = { coot::command_arg_t(imol), coot::command_arg_t(chain_id) };
      std::string sc = g.state_command("fitting", "backrub_rotamers_for_chain", args, lang);
      safe_python_command("import fitting");
      safe_python_command(sc);
   }
}

void find_helices_in_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {

}

void find_strands_in_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {

}

void
fill_partial_residues_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

}

void phosphorylate_this_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                       G_GNUC_UNUSED GVariant *parameter,
                                       G_GNUC_UNUSED gpointer user_data) {

}

void rebuild_fragment_using_dbloop_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                          GVariant *parameter,
                                          G_GNUC_UNUSED gpointer user_data) {

   if (parameter) {
      gchar *result;
      g_variant_get(parameter, "s", &result);
      std::string ss(result);
      std::cout << "db-loop size " << ss << std::endl;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
         int imol = pp.second.first;
         std::cout << "do something here with db loop fit for " << imol << " " << pp.second.second << std::endl;
      }
   }

}

void replace_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

}

void rigid_body_fit_residue_ranges_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                          G_GNUC_UNUSED GVariant *parameter,
                                          G_GNUC_UNUSED gpointer user_data) {

}

void rigid_body_fit_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                    G_GNUC_UNUSED GVariant *parameter,
                                    G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      rigid_body_refine_by_atom_selection(imol, "/");
      graphics_draw();
   }
}

void superpose_ligands_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   safe_python_command("import coot_gui");
   safe_python_command("coot_gui.superpose_ligand_gui()");
}

void symm_shift_reference_chain_here_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                            G_GNUC_UNUSED GVariant *parameter,
                                            G_GNUC_UNUSED gpointer user_data) {
   move_reference_chain_to_symm_chain_position();
}



void
draw_cell_and_symmetry_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *show_symm_window = wrapped_create_show_symmetry_window();
   if (show_symm_window) {
      set_transient_for_main_window(show_symm_window);
      gtk_widget_set_visible(show_symm_window, TRUE);
   }
}


void
display_only_active_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   display_only_active();
}

void
fullscreen_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {
   fullscreen();
}


#include "generic-display-objects-c.h"
void
generic_objects_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   generic_objects_gui_wrapper();
}



void
go_to_atom_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   // wrapped_create_show_symmetry_window() fills the window also
   GtkWidget *widget = wrapped_create_goto_atom_window(); // uses gtkbuilder

				/* now we need to fill the entry boxes
				   with default vaules and the option
				   menu according to molecules that
				   have coordinates. */

   gtk_widget_set_visible(widget, TRUE);
   graphics_info_t g;
}


void
label_atoms_in_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   label_atoms_in_residue();
}


void
label_CA_atoms_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {

   auto label_all_CAs = [] (int imol) {
      graphics_info_t::molecules[imol].add_labels_for_all_CAs();
   };

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      label_all_CAs(imol);
      graphics_draw();
   }
}


void
label_neighbours_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   label_neighbours();
}


void show_map_parameters_dialog(); // in glade-callbacks

void
map_parameters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {
   show_map_parameters_dialog();
}


void
ghost_control_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_ncs_control_dialog(); // uses builder
   gtk_widget_set_visible(w, TRUE);
}


void
spin_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   toggle_idle_spin_function();
}


void
rock_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   toggle_idle_rock_function();
}


void
ribbons_colour_by_chain_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "Chain";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}


void
ribbons_colour_rainbow_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "colorRampChainsScheme";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}


void
ribbons_colour_by_secondary_structure_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                             G_GNUC_UNUSED GVariant *parameter,
                                             G_GNUC_UNUSED gpointer user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "colorBySecondaryScheme";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}

void
screenshot_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {
    GtkWidget *file_chooser = coot_screendump_chooser();

    set_transient_and_position(COOT_UNDEFINED_WINDOW, file_chooser);
    g_object_set_data(G_OBJECT(file_chooser), "image_type", GINT_TO_POINTER(COOT_SCREENDUMP_SIMPLE));
    gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(file_chooser), "coot-screendump.tga");
    gtk_widget_set_visible(file_chooser, TRUE);
    check_for_dark_blue_density(); /* give a dialog if density it too dark (blue) */
}



void
sequence_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {

   g_warning("todo: after the removal of dynamic menus, this has to be reworked");
}

void
undo_last_navigation_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

   undo_last_move();
}

void
undo_symmetry_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {
   undo_symmetry_view();
}

void
clear_atom_labels_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {
   remove_all_atom_labels();
}


void
distances_and_angles_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

  GtkWidget* widget = widget_from_builder("geometry_frame");
  store_geometry_dialog(widget); /* needed to deactivate the distance
				    togglebutton after 2nd atoms
				    clicked in graphics */
  gtk_widget_set_visible(widget, TRUE);

}

void
pointer_distances_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = widget_from_builder("pointer_distances_dialog");
   fill_pointer_distances_widget(w);
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}

void
environment_distances_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *widget = widget_from_builder("environment_distance_dialog");
   fill_environment_widget(widget);
   gtk_widget_set_visible(widget, TRUE);

}


void
check_delete_waters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_check_waters_dialog();
   int imol_map = imol_refinement_map();
   gtk_widget_set_visible(w, TRUE);
   if (imol_map < 0)
      show_select_map_dialog();

}

void
difference_map_peaks_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *dialog = wrapped_create_generate_diff_map_peaks_dialog();
   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

}


void show_validation_graphs_dialog(G_GNUC_UNUSED GSimpleAction *simple_action, G_GNUC_UNUSED GVariant *parameter, G_GNUC_UNUSED gpointer user_data) {

   // 20230415-PE a common motif
   GtkWidget *di = widget_from_builder("validation_graph_dialog");
   GtkWidget *main_window = graphics_info_t::get_main_window();
   gtk_window_set_transient_for(GTK_WINDOW(di), GTK_WINDOW(main_window));

   auto get_first_model_molecule = [] () {
      graphics_info_t g;
      int imol = -1;
      int n_mol = g.molecules.size();
      for (int ii=0; ii<n_mol; ii++) {
         if (is_valid_model_molecule(ii))
            return ii;
      }
      return imol;
   };

   graphics_info_t g;
   GtkWidget *model_combobox = widget_from_builder("validation_graph_model_combobox");

   int imol = g.active_validation_graph_model_idx;
   if (! g.is_valid_model_molecule(imol))
      imol = get_first_model_molecule();

   // I don't think that it's imol that I want to use for the index.
   std::cout << "--------- in show_validation_graphs_dialog() " << model_combobox << " " << imol << std::endl;
   gtk_combo_box_set_active(GTK_COMBO_BOX(model_combobox), imol);

   gtk_widget_set_visible(di, TRUE);
}

void ramachandran_plot_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *di = widget_from_builder("ramachandran_plot_molecule_chooser_dialog");
   GtkWidget *main_window = graphics_info_t::get_main_window();
   gtk_window_set_transient_for(GTK_WINDOW(di), GTK_WINDOW(main_window));

   GtkWidget *model_combobox = widget_from_builder("ramachandran_plot_molecule_chooser_model_combobox");
   int imol_idx = 0; // not imol.
   gtk_combo_box_set_active(GTK_COMBO_BOX(model_combobox), imol_idx);
   gtk_widget_set_visible(di, TRUE);

}


void alignment_vs_pir_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

}

void atoms_with_zero_occupancies_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
   }
}

void atoms_overlaps_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
      coot_all_atom_contact_dots(imol);
   }

}

void all_atom_contact_dots_molprobity_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                             G_GNUC_UNUSED GVariant *parameter,
                                             G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
     int imol = pp.second.first;

      short int lang = coot::STATE_PYTHON;
      std::string module = "generic_objects";
      std::string function = "probe";
      std::vector<coot::command_arg_t> args = { coot::command_arg_t(imol)};
      std::string sc = g.state_command(module, function, args, lang);
      safe_python_command("import generic_objects");
      safe_python_command(sc);
   }
}

void highly_coordinates_waters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                          G_GNUC_UNUSED GVariant *parameter,
                                          G_GNUC_UNUSED gpointer user_data) {

   
   graphics_info_t g;
   short int lang = coot::STATE_PYTHON;
   std::string module = "coot_gui";
   std::string function = "water_coordination_gui";
   std::vector<coot::command_arg_t> args;
   std::string sc = g.state_command(module, function, args, lang);
   safe_python_command("import coot_gui");
   safe_python_command(sc);

}

#include "dynamic-validation.hh"

void overlaps_peptides_cbeta_ramas_and_rotas_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                                    G_GNUC_UNUSED GVariant *parameter,
                                                    G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   std::cout << "in overlaps_peptides_cbeta_ramas_and_rotas_action() with pp " << pp.first << " " << pp.second.first << std::endl;
   if (pp.first) {

      int imol = pp.second.first;
      overlaps_peptides_cbeta_ramas_and_rotas_internal(imol);
   }
   
   // graphics_info_t g;
   // std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   // if (pp.first) {
   //    int imol = pp.second.first;
   //    short int lang = coot::STATE_PYTHON;
   //    std::string module = "dynamic_atom_overlaps_and_other_outliers";
   //    std::string function = "quick_test_validation_outliers_dialog";
   //    std::vector<coot::command_arg_t> args = { coot::command_arg_t(imol)};
   //    std::string sc = g.state_command(module, function, args, lang);
   //    safe_python_command("import dynamic_atom_overlaps_and_other_outliers");
   //    safe_python_command(sc);
   // }

}

void refmac_log_validation_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                  G_GNUC_UNUSED GVariant *parameter,
                                  G_GNUC_UNUSED gpointer user_data) {
   info_dialog("Oops! No Refmac Log Validation yet");
}

void pepflips_from_difference_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                         G_GNUC_UNUSED GVariant *parameter,
                                         G_GNUC_UNUSED gpointer user_data) {
   pepflips_by_difference_map_dialog();
}



void validation_outliers_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {

      int imol = pp.second.first;
      int imol_map = imol_refinement_map();
      dynamic_validation_internal(imol, imol_map);

      // Goodbye wretched python
      //
      // short int lang = coot::STATE_PYTHON;
      // std::string module = "find_baddies";
      // std::string function = "validation_outliers_dialog";
      // std::vector<coot::command_arg_t> args = { coot::command_arg_t(imol), coot::command_arg_t(imol_map)};
      // std::string sc = g.state_command(module, function, args, lang);
      // safe_python_command("import find_baddies");
      // safe_python_command(sc);
   }

}

void
unmodelled_blobs_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_unmodelled_blobs_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}


void
remarks_browser_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_remarks_browser_molecule_chooser_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}


void
about_coot_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   std::cout << "About Coot" << std::endl;

   GtkWidget *dialog = widget_from_builder("about_dialog");
   if (dialog) {
      set_transient_for_main_window(dialog);
      gtk_widget_set_visible(dialog, TRUE);
   }
}

void
orthographic_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {
   set_use_perspective_projection(0);
}


void
perspective_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   set_use_perspective_projection(1);
}

void residue_type_selection_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                   G_GNUC_UNUSED GVariant *parameter,
                                   G_GNUC_UNUSED gpointer user_data) {

   std::cout << "residue_type_selection action" << std::endl;

}

void residues_with_alt_confs_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                   G_GNUC_UNUSED GVariant *parameter,
                                   G_GNUC_UNUSED gpointer user_data) {

   std::cout << "residues with alt confs action" << std::endl;
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
   }
}

void residues_with_cis_peptides_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                       G_GNUC_UNUSED GVariant *parameter,
                                       G_GNUC_UNUSED gpointer user_data) {

   std::cout << "residues with cis peptides action" << std::endl;
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
   }
}

void residues_with_missing_atoms_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

   std::cout << "residues with missing atoms action" << std::endl;
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
   }
}

#include "rsr-functions.hh"

void
refine_sphere(G_GNUC_UNUSED GSimpleAction *simple_action,
              G_GNUC_UNUSED GVariant *parameter,
              G_GNUC_UNUSED gpointer user_data) {

   rsr_sphere_refine();

}

void
refine_sphere_big(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   rsr_sphere_refine_plus();
}

void
refine_tandem_3(G_GNUC_UNUSED GSimpleAction *simple_action,
                G_GNUC_UNUSED GVariant *parameter,
                G_GNUC_UNUSED gpointer user_data) {

   rsr_refine_tandem_3();
}

void
refine_tandem_5(G_GNUC_UNUSED GSimpleAction *simple_action,
                G_GNUC_UNUSED GVariant *parameter,
                G_GNUC_UNUSED gpointer user_data) {

   rsr_refine_tandem_5();

}

void
refine_single_residue(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {

   rsr_refine_residue();
}

void
refine_chain(G_GNUC_UNUSED GSimpleAction *simple_action,
             G_GNUC_UNUSED GVariant *parameter,
             G_GNUC_UNUSED gpointer user_data) {

   rsr_refine_chain();
}

void
refine_all_atoms(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {

   rsr_refine_all_atoms();
}

void
refine_with_range_picked_atoms() {

   graphics_info_t g;
   std::string alt_conf; // needs to be set correctly
   short int is_water_flag = false; // needs to be set correctly

   const std::string &alt_conf_1 = g.in_range_first_picked_atom.alt_conf;
   const std::string &alt_conf_2 = g.in_range_second_picked_atom.alt_conf;

   if (alt_conf_1 == alt_conf_2)
      if (! alt_conf_1.empty())
         alt_conf = alt_conf_1;

   int imol_1 = g.in_range_first_picked_atom.int_user_data;
   int imol_2 = g.in_range_second_picked_atom.int_user_data;

   const std::string &chain_id_1 = g.in_range_first_picked_atom.chain_id;
   const std::string &chain_id_2 = g.in_range_second_picked_atom.chain_id;

   int res_no_1 = g.in_range_first_picked_atom.res_no;
   int res_no_2 = g.in_range_second_picked_atom.res_no;

   const std::string &ins_code_1 = g.in_range_first_picked_atom.ins_code;
   const std::string &ins_code_2 = g.in_range_second_picked_atom.ins_code;

   if (g.is_valid_model_molecule(imol_1)) {

      if (imol_1 == imol_2) {
         g.refine_residue_range(imol_1, chain_id_1, chain_id_2, res_no_1, ins_code_1, res_no_2, ins_code_2,
                                alt_conf, is_water_flag);
      }
   }

}

void
refine_range(G_GNUC_UNUSED GSimpleAction *simple_action,
             G_GNUC_UNUSED GVariant *parameter,
             G_GNUC_UNUSED gpointer user_data) {

   graphics_info_t g;
   std::cout << "in refine_range with in_range_define " << g.in_range_define << std::endl;
   if (g.in_range_define == 2) {
      // so what were the two atoms?

      refine_with_range_picked_atoms();

   } else {
      std::string m = "Use the Range button to define a residue range (pick 2 atoms)";
      g.add_status_bar_text(m);
   }
}

void
repeat_refine_range(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {

   refine_with_range_picked_atoms();
}

void
refine_regularize_sphere(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {

   regularize_sphere();
}

void
refine_regularize_tandem_3(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   regularize_tandem_3();
}


void
refine_regularize_single_residue(G_GNUC_UNUSED GSimpleAction *simple_action,
                G_GNUC_UNUSED GVariant *parameter,
                G_GNUC_UNUSED gpointer user_data) {

   regularize_residue();
}


void
fix_atom(GSimpleAction *simple_action,
         GVariant *parameter,
         gpointer user_data) {
   
   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
      g.attach_buffers(); // 20220823-PE needed?
      g.mark_atom_as_fixed(imol, pp.second.second, true);
      g.graphics_draw(); // maybe not needed here
   }
}


void
unfix_atom(GSimpleAction *simple_action,
           GVariant *parameter,
           gpointer user_data) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
      g.attach_buffers(); // 20220823-PE needed?
      g.mark_atom_as_fixed(imol, pp.second.second, true);
      g.graphics_draw(); // maybe not needed here
   }
}


void
unfix_all_atoms(GSimpleAction *simple_action,
                GVariant *parameter,
                gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      graphics_info_t g;
      g.molecules[imol].clear_all_fixed_atoms();
   }
}


void
mutate_to_type(GSimpleAction *simple_action,
               GVariant *parameter,
               gpointer user_data) {

   if (parameter) {
      gchar *result;
      g_variant_get(parameter, "s", &result);
      std::string ss(result);
      std::cout << "mutate_to type parameter " << ss << std::endl;
      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
         int imol = pp.second.first;
         g.mutate_residue_imol = imol;
         g.mutate_auto_fit_residue_imol = imol;
         coot::residue_spec_t res_spec(pp.second.second);
         g.do_mutation(imol, res_spec, ss, false); // not stub
      }
   }
}

void
delete_item(GSimpleAction *simple_action,
            GVariant *parameter,
            gpointer user_data) {

   auto delete_residue_range = [] () {

      graphics_info_t g;
      int imol_1 = g.in_range_first_picked_atom.int_user_data;
      int imol_2 = g.in_range_second_picked_atom.int_user_data;
      if (g.is_valid_model_molecule(imol_1)) {
         if (imol_1 == imol_2) {
            coot::residue_spec_t rs1(g.in_range_first_picked_atom);
            coot::residue_spec_t rs2(g.in_range_second_picked_atom);
            g.delete_residue_range(imol_1, rs1, rs2);
         }
      }
   };

   if (parameter) {
      gchar *result;
      g_variant_get(parameter, "s", &result);
      std::string par(result);
      std::cout << "debug:: delete_item parameter " << par << std::endl;
      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
      if (pp.first) {
         auto atom_spec = pp.second.second;
         coot::residue_spec_t res_spec(atom_spec);
         int imol = pp.second.first;
         if (par == "atom") {
            auto &m = g.molecules[imol];
            // change this signature to use an atom spec.
            m.delete_atom(atom_spec);
            g.graphics_draw();
         }
         if (par == "residue") {
            g.delete_active_residue(); // does a redraw
         }
         if (par == "chain") {
            auto &m = g.molecules[imol];
            m.delete_chain(atom_spec.chain_id);
            g.graphics_draw();
         }
         if (par == "hydrogen-atoms-in-residue") {
            auto &m = g.molecules[imol];
            // change this signature to use an residue spec.
            std::cout << "callign delete_residue_hydrogens()!!!! " << res_spec << std::endl;
            m.delete_residue_hydrogens(res_spec.chain_id, res_spec.res_no, res_spec.ins_code, atom_spec.alt_conf);
            graphics_draw();
         }
         if (par == "residue-range") {
            // use old-style "setup"
            // Needs "check_if_in_range_defines" to be working.
            // Here we need to turn on the expecting the delet residue range "start" flag
            // and unset the others c.f. set_delete_residue_zone_mode()

            delete_residue_range();
            g.graphics_draw();

         }
         if (par == "side-chain") {
            auto &m = g.molecules[imol];
            // change this signature to use an residue spec.
            m.delete_residue_sidechain(res_spec.chain_id, res_spec.res_no, res_spec.ins_code);
            g.graphics_draw();
         }
         if (par == "side-chain-residue-range") {
            // use old-style "setup"
            std::cout << "delete side-chain-residue-range needs fixing" << std::endl;
         }
         if (par == "side-chains-in-chain") {
            delete_sidechains_for_chain(imol, atom_spec.chain_id);
         }
         if (par == "water") {

            std::cout << "....................... delete water! " << atom_spec << std::endl;
            auto &m = g.molecules[imol];
            m.delete_water(atom_spec);
            graphics_draw();
         }
      }
   }
}

void
create_actions(GtkApplication *application) {

   auto add_action = [application] (const std::string &action_name,
                                    void (*action_function) (GSimpleAction *simple_action,
                                                             GVariant *parameter,
                                                             gpointer user_data)) {
      GSimpleAction *simple_action = g_simple_action_new(action_name.c_str(), NULL);
      g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
      g_signal_connect(simple_action, "activate", G_CALLBACK(action_function), NULL);
   };

   auto add_action_with_param = [application] (const std::string &action_name,
                                    void (*action_function) (GSimpleAction *simple_action,
                                                             GVariant *parameter,
                                                             gpointer user_data)) {
      GSimpleAction *simple_action = g_simple_action_new(action_name.c_str(), G_VARIANT_TYPE_STRING);
      g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
      g_signal_connect(simple_action, "activate", G_CALLBACK(action_function), NULL);
   };


   // File

   add_action(     "open_coordinates_action",      open_coordinates_action);
   add_action(        "auto_open_mtz_action",         auto_open_mtz_action);
   add_action(         "open_dataset_action",          open_dataset_action);
   add_action(             "open_map_action",              open_map_action);
   add_action("import_cif_dictionary_action", import_cif_dictionary_action);

   add_action("get_monomer_action", get_monomer_action);
   add_action(     "curlew_action",      curlew_action);
   add_action(       "exit_action",        exit_action);

   add_action(             "search_monomer_library_action",           search_monomer_library_action);
   add_action_with_param("show_accession_code_fetch_frame",         show_accession_code_fetch_frame);
   add_action(    "fetch_pdbe_ligand_description_action",    fetch_pdbe_ligand_description_action);
   add_action( "fetch_and_superpose_alphafold_models_action", fetch_and_superpose_alphafold_models_action);
   add_action(              "fetch_map_from_emdb_action",              fetch_map_from_emdb_action);
   add_action(                 "save_coordinates_action",                 save_coordinates_action);
   add_action(        "save_symmetry_coordinates_action",        save_symmetry_coordinates_action);
   add_action(                       "save_state_action",                       save_state_action);
   add_action(                  "recover_session_action",                  recover_session_action);
   add_action(                  "file_export_map_action",                  file_export_map_action);
   add_action(         "file_export_map_fragment_action",                  file_export_map_action);
   add_action(                   "close_molecule_action",                   close_molecule_action);

   // Edit

   add_action(                  "make_link_action",                   make_link_action); // add header link
   add_action(           "change_chain_ids_action",            change_chain_ids_action);
   add_action(              "copy_molecule_action",               copy_molecule_action);
   add_action(     "copy_molecule_fragment_action",      copy_molecule_fragment_action);
   add_action(            "merge_molecules_action",             merge_molecules_action);
   add_action(         "move_molecule_here_action",          move_molecule_here_action);
   add_action(            "mutate_molecule_action",             mutate_molecule_action);
   add_action(      "edit_replace_fragment_action",       edit_replace_fragment_action);
   add_action(       "edit_replace_residue_action",        edit_replace_residue_action);
   add_action(          "renumber_residues_action",           renumber_residues_action);
   add_action(            "renumber_waters_action",             renumber_waters_action);
   add_action(               "residue_info_action",                residue_info_action);
   add_action(            "edit_restraints_action",             edit_restraints_action);
   add_action(           "show_preferences_action",            show_preferences_action);
   add_action(       "merge_solvent_chains_action",        merge_solvent_chains_action);
   add_action(      "undo_molecule_chooser_action",       undo_molecule_chooser_action);
   add_action(    "show_shader_preferences_action",     show_shader_preferences_action);
   add_action(    "fix_nomenclature_errors_action",     fix_nomenclature_errors_action);
   add_action(  "invert_this_chiral_centre_action",    invert_this_chiral_centre_action);
   add_action("exchange_chain_ids_for_seg_ids_action", exchange_chain_ids_for_seg_ids_action);

   // Calculate

   add_action(       "align_and_mutate_action",        align_and_mutate_action);
   add_action(   "fit_loop_by_database_search",     fit_loop_by_database_search);
   add_action("fit_loop_by_ramachandran_search",fit_loop_by_ramachandran_search);
   add_action(         "ligand_builder_action",          ligand_builder_action);
   add_action(          "lsq_superpose_action",           lsq_superpose_action);
   add_action(             "run_script_action",              run_script_action);
   add_action(      "ssm_superposition_action",       ssm_superposition_action);
   add_action("calculate_updating_maps_action", calculate_updating_maps_action);
   add_action(  "sharpen_blur_for_xray_action",   sharpen_blur_for_xray_action);
   add_action(       "scripting_python_action",        scripting_python_action);
   add_action(       "scripting_scheme_action",        scripting_scheme_action);

   add_action("calculate_hydrogen_bonds_action", calculate_hydrogen_bonds_action);
   add_action(          "load_tutorial_model_and_data_action",           load_tutorial_model_and_data_action);
   add_action("use_clustalw_for_alignment_then_mutate_action", use_clustalw_for_alignment_then_mutate_action);

   // Calculate -> Map Tools

   add_action(                        "copy_map_action",                         copy_map_action);
   add_action(                    "average_maps_action",                     average_maps_action);
   add_action(                    "multichicken_action",                     multichicken_action);
   add_action(                   "brighten_maps_action",                    brighten_maps_action);
   add_action(           "another_contour_level_action",            another_contour_level_action);
   add_action(           "make_a_difference_map_action",            make_a_difference_map_action);
   add_action(       "set_map_is_difference_map_action",        set_map_is_difference_map_action);
   add_action(      "mask_map_by_atom_selection_action",       mask_map_by_atom_selection_action);
   add_action(   "make_a_smoother_copy_of_a_map_action",    make_a_smoother_copy_of_a_map_action);
   add_action(  "transform_map_by_lsq_model_fit_action",   transform_map_by_lsq_model_fit_action);
   add_action("make_a_very_smooth_copy_of_a_map_action", make_a_very_smooth_copy_of_a_map_action);

   // Calculate -> Modules

   add_action(        "add_HOLE_module_action",         add_HOLE_module_action);
   add_action(        "add_ccp4_module_action",         add_ccp4_module_action);
   add_action("add_carbohydrate_module_action", add_carbohydrate_module_action);
   add_action(     "add_cryo_em_module_action",      add_cryo_em_module_action);
   add_action(    "add_prosmart_module_action",     add_prosmart_module_action);
   add_action(      "add_ligand_module_action",       add_ligand_module_action);
   add_action(      "add_rcrane_module_action",       add_rcrane_module_action);
   add_action(  "add_restraints_module_action",   add_restraints_module_action);
   add_action(      "add_refine_module_action",       add_refine_module_action);
   add_action(       "add_shelx_module_action",        add_shelx_module_action);
   add_action(       "add_views_module_action",        add_views_module_action);

   // Calculate -> NCS

   add_action("copy_ncs_residue_range_action", copy_ncs_residue_range_action);
   add_action(        "copy_ncs_chain_action",         copy_ncs_chain_action);
   add_action(           "ncs_jumping_action",            ncs_jumping_action);
   add_action(           "ncs_ligands_action",            ncs_ligands_action);

   // Calculate -> Assign Sequnence

   add_action("associate_sequence_file_action", associate_sequence_file_action);
   add_action(        "assign_sequence_action",         assign_sequence_action);

   // Modelling

   add_action(                    "add_an_atom_action",                     add_an_atom_action);
   add_action(             "add_hydrogen_atoms_action",              add_hydrogen_atoms_action);
   add_action("add_hydrogen_atoms_using_refmac_action", add_hydrogen_atoms_using_refmac_action);
   add_action(    "add_other_solvent_molecules_action",     add_other_solvent_molecules_action);
   add_action(                   "find_ligands_action",                    find_ligands_action);
   add_action(                    "find_waters_action",                     find_waters_action);
   add_action(                 "dna_rna_models_action",                  dna_rna_models_action);
   add_action(               "place_helix_here_action",                place_helix_here_action);
   add_action(              "cis_trans_convert_action",               cis_trans_convert_action);
   add_action(             "add_OXT_to_residue_action",              add_OXT_to_residue_action);
   add_action(        "reverse_chain_direction_action",         reverse_chain_direction_action);
   add_action(  "arrange_waters_around_protein_action",   arrange_waters_around_protein_action);
   add_action("assign_hetatms_for_this_residue_action", assign_hetatms_for_this_residue_action);
   add_action(    "assign_hetatoms_to_molecule_action",     assign_hetatoms_to_molecule_action);
   add_action(     "backrub_rotamers_for_chain_action",      backrub_rotamers_for_chain_action);
   add_action(            "find_helices_in_map_action",             find_helices_in_map_action);
   add_action(            "find_strands_in_map_action",             find_strands_in_map_action);
   add_action(          "fill_partial_residues_action",           fill_partial_residues_action);
   add_action(     "phosphorylate_this_residue_action",      phosphorylate_this_residue_action);
   add_action(                "replace_residue_action",                 replace_residue_action);
   add_action(  "rigid_body_fit_residue_ranges_action",   rigid_body_fit_residue_ranges_action);
   add_action(        "rigid_body_fit_molecule_action",         rigid_body_fit_molecule_action);
   add_action(              "superpose_ligands_action",               superpose_ligands_action);
   add_action("symm_shift_reference_chain_here_action", symm_shift_reference_chain_here_action);
   add_action(          "other_modelling_tools_action",           other_modelling_tools_action);
   add_action(                     "whats_this_action",                      whats_this_action);

   add_action_with_param("rebuild_fragment_using_dbloop_action", rebuild_fragment_using_dbloop_action);

   // Draw

   add_action(       "background_black_action",        background_black_action);
   add_action(  "background_light_grey_action",   background_light_grey_action);
   add_action(   "background_dark_grey_action",    background_dark_grey_action);
   add_action(       "background_white_action",        background_white_action);
   add_action(    "display_only_active_action",     display_only_active_action);
   add_action(        "bond_parameters_action",         bond_parameters_action);
   add_action(           "bond_colours_action",            bond_colours_action);
   add_action(             "fullscreen_action",              fullscreen_action);
   add_action(             "go_to_atom_action",              go_to_atom_action);
   add_action(         "label_CA_atoms_action",          label_CA_atoms_action);
   add_action(         "map_parameters_action",          map_parameters_action);
   add_action(        "generic_objects_action",         generic_objects_action);
   add_action(       "label_neighbours_action",        label_neighbours_action);
   add_action( "label_atoms_in_residue_action",  label_atoms_in_residue_action);
   add_action( "draw_cell_and_symmetry_action",  draw_cell_and_symmetry_action);

   add_action(   "ghost_control_action",      ghost_control_action);
   add_action(        "spin_view_action",         spin_view_action);
   add_action(        "rock_view_action",         rock_view_action);
   add_action( "perspective_view_action",  perspective_view_action);
   add_action("orthographic_view_action", orthographic_view_action);

   add_action(     "residue_type_selection_action",      residue_type_selection_action);
   add_action(    "residues_with_alt_confs_action",     residues_with_alt_confs_action);
   add_action( "residues_with_cis_peptides_action",  residues_with_cis_peptides_action);
   add_action("residues_with_missing_atoms_action", residues_with_missing_atoms_action);

   add_action("screenshot_action",                           screenshot_action);
   add_action("sequence_view_action",                     sequence_view_action);
   add_action("undo_symmetry_view_action",           undo_symmetry_view_action);
   add_action("undo_last_navigation_action",       undo_last_navigation_action);
   add_action("ribbons_colour_by_chain_action", ribbons_colour_by_chain_action);
   add_action("ribbons_colour_rainbow_action",   ribbons_colour_rainbow_action);
   add_action("ribbons_colour_by_secondary_structure_action", ribbons_colour_by_secondary_structure_action);

   add_action( "toggle_display_frames_per_second_action", toggle_display_frames_per_second_action);

   // Measures

   add_action(    "clear_atom_labels_action",     clear_atom_labels_action);
   add_action(    "pointer_distances_action",     pointer_distances_action);
   add_action( "distances_and_angles_action",  distances_and_angles_action);
   add_action("environment_distances_action", environment_distances_action);

   // Validate

   add_action(                "unmodelled_blobs_action",                 unmodelled_blobs_action);
   add_action(             "check_delete_waters_action",              check_delete_waters_action);
   add_action(            "difference_map_peaks_action",             difference_map_peaks_action);
   add_action(          "show_validation_graphs_dialog",           show_validation_graphs_dialog);
   add_action(               "ramachandran_plot_action",                ramachandran_plot_action);
   add_action(                "alignment_vs_pir_action",                 alignment_vs_pir_action);
   add_action(                  "atoms_overlaps_action",                   atoms_overlaps_action);
   add_action(             "validation_outliers_action",              validation_outliers_action);
   add_action(           "refmac_log_validation_action",            refmac_log_validation_action);
   add_action(       "highly_coordinates_waters_action",        highly_coordinates_waters_action);
   add_action(     "atoms_with_zero_occupancies_action",      atoms_with_zero_occupancies_action);
   add_action(    "pepflips_from_difference_map_action",     pepflips_from_difference_map_action);
   add_action("all_atom_contact_dots_molprobity_action", all_atom_contact_dots_molprobity_action);
   add_action("overlaps_peptides_cbeta_ramas_and_rotas_action", overlaps_peptides_cbeta_ramas_and_rotas_action);

   // About

   add_action("remarks_browser_action", remarks_browser_action);
   add_action("about_coot_action", about_coot_action);

   // Refine menu

   add_action("refine_sphere",                    refine_sphere);
   add_action("refine_sphere_big",                refine_sphere_big);
   add_action("refine_tandem_3",                  refine_tandem_3);
   add_action("refine_tandem_5",                  refine_tandem_5);
   add_action("refine_single_residue",            refine_single_residue);
   add_action("refine_chain",                     refine_chain);
   add_action("refine_all_atoms",                 refine_all_atoms);
   add_action("refine_range",                     refine_range);
   add_action("repeat_refine_range",              repeat_refine_range);
   add_action("refine_regularize_sphere",         refine_regularize_sphere);
   add_action("refine_regularize_tandem_3",       refine_regularize_tandem_3);
   add_action("refine_regularize_single_residue", refine_regularize_single_residue);

   // Fix Atoms

   add_action(  "fix_atom",        fix_atom);
   add_action("unfix_atom",      unfix_atom);
   add_action("unfix_all_atoms", unfix_all_atoms);

   // Mutate menu
   add_action_with_param("mutate_to_type", mutate_to_type);

   // Delete menu
   add_action_with_param("delete_item", delete_item);
}
