/* src/main.cc
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007, 2009, 2011, 2012 by The University of Oxford
 * Copyright 2014, 2015 by Medical Research Council
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

#include <gtk/gtk.h>


#include "coot-setup-python.hh"
#include "python-classes.hh"

#include "compat/coot-sysdep.h"

/*
 * Initial main.c file generated by Glade. Edit as required.
 * Glade will not overwrite this file.
 */

#include <sys/time.h>
#include <string.h> // strcmp

#include <iostream>

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif


// We are not using NLS yet.
// #ifndef WINDOWS_MINGW
// #define ENABLE_NLS
// #endif
// #ifdef DATADIR
// #endif // DATADIR

// #include <GL/glut.h> // for glutInit()

// #include "lbg/lbg.hh"

#include "interface.h"
#ifndef HAVE_SUPPORT_H
#define HAVE_SUPPORT_H
#include "support.h"
#endif /* HAVE_SUPPORT_H */


#include <sys/types.h> // for stating
#include <sys/stat.h>

#ifndef _MSC_VER
#include <unistd.h>
#else
#define PKGDATADIR "C:/coot/share"
#endif

// not in the world of GTK4
// #include "globjects.h"

#include <vector>
#include <string>

// #include <mmdb2/mmdb_manager.h>
// #include "coords/mmdb-extras.h"
// #include "coords/mmdb.h"
// #include "coords/mmdb-crystal.h"

#include "clipper/core/test_core.h"
#include "clipper/contrib/test_contrib.h"

#include "utils/coot-utils.hh"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "command-line.hh"

#include "graphics-info.h"
// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
// BL says:: and (2.3 - dewinter), i.e. is a Mac - Python issue
// since the follwing two include python graphics-info.h is moved up
//
#if defined (WINDOWS_MINGW)
#include <locale.h>
#ifdef DATADIR
#undef DATADIR
#endif // DATADIR
#endif
#include "compat/sleep-fixups.h"

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "c-interface-preferences.h"

// #include "coot-surface/rgbreps.h"

#include "coot-database.hh"

#include <glob.h>

#ifdef USE_GUILE
#include <libguile.h>
#endif


// #include "c-inner-main.h"
#include "coot-glue.hh"

#include "rotate-translate-modes.hh"

#include "change-dir.hh"
#include "curlew.hh"

void show_citation_request();
void load_gtk_resources();
void setup_splash_screen();
void desensitive_scripting_menu_item_maybe(GtkWidget *window);
int setup_screen_size_settings();
void setup_application_icon(GtkWindow *window);
void setup_symm_lib();
void check_reference_structures_dir();
#include "boot-python.hh"

#ifdef USE_MYSQL_DATABASE
#include "mysql/mysql.h"
int setup_database();
#endif

#include "testing.hh" // for test_internal();

#include "scm-boot-guile.hh"

#include "widget-headers.hh" // put these somewhere else? better name? -------- GTK-FIME
#include "sound.hh"

#include "draw.hh" // for test_gtk3_adjustment_changed() - maybe that should go elsewhere?
#include "draw-2.hh"

#include "widget-from-builder.hh"

void
windows_set_error_mode() {

#ifdef WINDOWS_MINGW
      // in Windows we don't want a crash dialog if no-graphics
      SetErrorMode(SetErrorMode(SEM_NOGPFAULTERRORBOX) | SEM_NOGPFAULTERRORBOX);
#endif // MINGW
}




int
do_self_tests() {

   std::cout << "INFO:: Running internal self tests" << std::endl;
   // return true on success
   clipper::Test_core test_core;       bool result_core    = test_core();
   clipper::Test_contrib test_contrib; bool result_contrib = test_contrib();
   std::cout<<" INFO:: Test Clipper core   : "<<(result_core   ?"OK":"FAIL")<<std::endl;
   std::cout<<" INFO:: Test Clipper contrib: "<<(result_contrib?"OK":"FAIL")<<std::endl;
   // return 1 on success
   int gis = test_internal();
   int shell_exit_code = 1;
   if (result_core)
      if (result_contrib)
         if (gis == 1)
            shell_exit_code = 0;
   return shell_exit_code;

}


#include "widget-from-builder.hh"
#include "coot-application.hh"

int new_startup(int argc, char **argv);

// This main is used for both python/guile useage and unscripted.
int
main(int argc, char *argv[]) {

   new_startup(argc, argv);
   return 0;

}



void
show_citation_request() {

   std::cout << "\n   If you have found this software to be useful, you are requested to cite:\n"
	     << "   Coot: model-building tools for molecular graphics" << std::endl;
   std::cout << "   Emsley P, Cowtan K" << std::endl;
   std::cout << "   ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY\n";
   std::cout << "   60: 2126-2132 Part 12 Sp. Iss. 1 DEC 2004\n\n";

   std::cout << "   The reference for the REFMAC5 Dictionary is:\n";
   std::cout << "   REFMAC5 dictionary: organization of prior chemical knowledge and\n"
	     << "   guidelines for its use" << std::endl;
   std::cout << "   Vagin AA, Steiner RA, Lebedev AA, Potterton L, McNicholas S,\n"
	     << "   Long F, Murshudov GN" << std::endl;
   std::cout << "   ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY " << std::endl;
   std::cout << "   60: 2184-2195 Part 12 Sp. Iss. 1 DEC 2004" << std::endl;

#ifdef HAVE_SSMLIB
    std::cout << "\n   If using \"SSM Superposition\", please cite:\n";

    std::cout << "   Secondary-structure matching (SSM), a new tool for fast\n"
	      << "   protein structure alignment in three dimensions" << std::endl;
    std::cout << "   Krissinel E, Henrick K" << std::endl;
    std::cout << "   ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY" << std::endl;
    std::cout << "   60: 2256-2268 Part 12 Sp. Iss. 1 DEC 2004\n" << std::endl;
#endif // HAVE_SSMLIB

}

