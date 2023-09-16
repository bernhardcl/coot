
#include <iostream>
#include <gtk/gtk.h>
#include "utils/coot-utils.hh"
#include "layla-interface.hh"

void
layla() {

   std::string dir = coot::package_data_dir();
   // all ui files should live here:
   std::string dir_ui = coot::util::append_dir_dir(dir, "ui");
   std::string ui_file_name = "layla.ui";
   std::string ui_file_full = coot::util::append_dir_file(dir_ui, ui_file_name);
   // allow local override
   if (coot::file_exists(ui_file_name)) ui_file_full = ui_file_name;
   GError* error = NULL;
   GtkBuilder* builder = gtk_builder_new();
   gboolean status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);
   if (status == FALSE) {
      std::cout << "ERROR:: Failure to read or parse " << ui_file_full << std::endl;
      std::cout << error->message << std::endl;
      return;
   }

   std::cout << "Now do something with Laya" << std::endl;
   
}
