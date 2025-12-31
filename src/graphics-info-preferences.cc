/* src/graphics-info-preferences.cc
 * 
 * Copyright 2008 The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#if defined _MSC_VER
#include <windows.h>
#endif

#include <fstream>

#include <gtk/gtk.h>
#include "interface.h"
#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-preferences.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh"
#include "coot-preferences.h"
#include "utils/coot-utils.hh"
#include "widget-from-builder.hh"


std::string
graphics_info_t::get_preferences_directory() const {

   xdg_t xdg;
   std::string preferences_dir = xdg.get_config_home().u8string();
   std::string pkgdatadir = coot::package_data_dir();

   std::string fn;

   if (preferences_dir.empty()) {
      fn = coot::util::append_dir_dir(pkgdatadir, ".coot");
   } else {
      fn = preferences_dir;
   }

   return fn;
}

void
graphics_info_t::add_to_preferences(const std::string &file_name, const std::string &contents) const {

   std::string pref_dir = get_preferences_directory();
   std::string pref_subdir = coot::util::append_dir_dir(pref_dir, "preferences");
   std::string fn = coot::util::append_dir_file(pref_subdir, file_name);

   std::ofstream f(fn.c_str());
   if (f) {
      f << contents << std::endl;
   }
   f.close();
}


