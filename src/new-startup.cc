
#include <iostream>
#include <string>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "graphics-info.h"
#include "create-menu-item-actions.hh"

#include "coot-setup-python.hh"

void print_opengl_info();

void init_framebuffers(GtkWidget *glarea) {

   // put this into graphics-info I suppose.

   std::cout << "DEBUG:: use_framebuffers: " << graphics_info_t::use_framebuffers << std::endl;

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

      if (graphics_info_t::use_framebuffers) {
         unsigned int index_offset = 0;
         GLenum err;
         graphics_info_t::screen_framebuffer.init(w, h, index_offset, "screen/occlusion");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 1;
         graphics_info_t::blur_y_framebuffer.init(w, h, index_offset, "blur-y");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_y_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 2;
         graphics_info_t::blur_x_framebuffer.init(w, h, index_offset, "blur-x");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_x_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 3;
         graphics_info_t::combine_textures_using_depth_framebuffer.init(w, h, index_offset, "new-blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_combine framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 4;
         graphics_info_t::blur_framebuffer.init(w, h, index_offset, "blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                                << err << std::endl;
      }

}


#include "text-rendering-utils.hh"

void
new_startup_realize(GtkWidget *gl_area) {

   std::cout << "new_startup_realize() ------------------- start ------------------"
             << std::endl;

   gtk_gl_area_make_current(GTK_GL_AREA (gl_area));

   if (gtk_gl_area_get_error(GTK_GL_AREA (gl_area)) != NULL)
      return;

   GdkGLContext *context = gtk_gl_area_get_context(GTK_GL_AREA(gl_area));

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(gl_area), TRUE);

   print_opengl_info();

   graphics_info_t g;
   init_framebuffers(gl_area); // Hmm - I don't know what this does compared to below.
   g.init_framebuffers();
   g.init_buffers();
   g.init_joey_ssao_stuff();
   g.init_shaders();
   g.setup_lights();

   float x_scale = 4.4;  // what are these numbers!?
   float y_scale = 1.2;
   x_scale = 1.002;
   y_scale = 1.002;
   g.tmesh_for_labels.setup_camera_facing_quad(x_scale, y_scale, 0.0, 0.0);
   g.setup_hud_geometry_bars();
   g.setup_hud_buttons();
   g.setup_rama_balls();
   g.setup_key_bindings();
   float double_rama_size = 0.8; // scaled by 0.5 in the gl-rama draw call.
   g.gl_rama_plot.setup_buffers(double_rama_size); // rama relative size, put it into graphics_info_t
   // and allow it to be set in the API
   g.setup_draw_for_happy_face_residue_markers_init();
   g.setup_draw_for_bad_nbc_atom_pair_markers();
   g.setup_draw_for_anchored_atom_markers_init();
   g.lines_mesh_for_hud_lines.set_name("lines mesh for fps graph");
   unsigned int frame_time_history_list_max_n_elements = 500;
   // +40 for base and grid lines
   std::vector<s_generic_vertex> empty_vertices(frame_time_history_list_max_n_elements + 40);
   std::vector<unsigned int> empty_indices(1500, 0); // or some number
   g.lines_mesh_for_hud_lines.setup_vertices_and_indices(empty_vertices, empty_indices);

   int w = 500;
   int h = 500;

   setup_hud_text(w, h, graphics_info_t::shader_for_hud_text, false);
   setup_hud_text(w, h, graphics_info_t::shader_for_atom_labels, true);
   g.tmesh_for_hud_refinement_dialog_arrow = HUDTextureMesh("HUD tmesh for refinement dialog arrow");
   g.tmesh_for_hud_refinement_dialog_arrow.setup_quad();
   g.texture_for_hud_refinement_dialog_arrow             = Texture("refinement-dialog-arrrow.png", Texture::DIFFUSE);
   g.texture_for_hud_refinement_dialog_arrow_highlighted = Texture("refinement-dialog-arrrow-highlighted.png", Texture::DIFFUSE);

   g.tmesh_for_shadow_map.setup_quad();

   g.setup_key_bindings();

   GLenum err = glGetError();
   if (err)
      std::cout << "ERROR:: new_startup_realize() --start-- err is " << err << std::endl;

   // Hmm! - causes weird graphics problems
   // setup_python(0, NULL); // needs to called after GTK has started - because it depends on gtk.
                             // 20220629-PE not at the moment though - I removed the gobject parts from the code path
}


void
new_startup_unrealize(GtkWidget *widget) {

   gtk_gl_area_make_current (GTK_GL_AREA (widget));
   if (gtk_gl_area_get_error (GTK_GL_AREA (widget)) != NULL)
      return;

}


gboolean
new_startup_on_glarea_render(GtkGLArea *glarea) {

   // std::cout << "DEBUG: new_startup_on_glarea_render()!" << std::endl;
   bool screen_dump_frame_buffer = false;
   return graphics_info_t::render(screen_dump_frame_buffer);
}


void
new_startup_on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   std::cout << "DEBUG: new_startup_on_glarea_resize() " <<  width << " " << height << std::endl;
   graphics_info_t g;
   // for the GL widget, not the window.
   g.graphics_x_size = width;
   g.graphics_y_size = height;
   // g.reset_frame_buffers(width, height); // currently makes the widget blank (not drawn)
   g.init_shaders();

}

// void on_glarea_realize(GtkWidget *widget); // using this give linking problems.

GtkWidget *new_startup_create_glarea_widget() {

   GtkWidget *gl_area = gtk_gl_area_new();
   g_signal_connect(gl_area, "realize",   G_CALLBACK(new_startup_realize),   NULL);
   g_signal_connect(gl_area, "unrealize", G_CALLBACK(new_startup_unrealize), NULL);
   g_signal_connect(gl_area, "render",    G_CALLBACK(new_startup_on_glarea_render),  NULL);
   g_signal_connect(gl_area, "resize",    G_CALLBACK(new_startup_on_glarea_resize),  NULL);

   gtk_widget_set_can_focus(gl_area, TRUE);
   gtk_widget_set_focusable(gl_area, TRUE);

   gtk_widget_set_hexpand(gl_area, TRUE);
   gtk_widget_set_vexpand(gl_area, TRUE);

   return gl_area;

}

void
on_open_clicked(GSimpleAction *action,
                GVariant *parameter,
                gpointer data) {

   std::cout << "open clicked" << std::endl;

}

void
on_close_clicked(GSimpleAction *action,
                 GVariant *parameter,
                 gpointer data) {

   std::cout << "close clicked" << std::endl;

}


GMenu *create_menu_by_hand(const GtkApplication *application) {
   const GActionEntry entries[] = {
      { "open",  on_open_clicked,  NULL, NULL, NULL, { 0, 0, 0 } },
      { "close", on_close_clicked, NULL, NULL, NULL, { 0, 0, 0 } }
   };

   g_action_map_add_action_entries(G_ACTION_MAP(application), entries, G_N_ELEMENTS(entries), NULL);

   GMenu *menu      = g_menu_new();
   GMenu *file_menu = g_menu_new();

   GMenuItem *item;
   item = g_menu_item_new("Open", "app.open");
   g_menu_append_item(file_menu, item);

   item = g_menu_item_new("Close", "app.close");
   g_menu_append_item(file_menu, item);

   g_menu_append_submenu(menu, "File", G_MENU_MODEL(file_menu));

   return menu;
}


void on_glarea_drag_begin_primary(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_begin_primary(gesture, x, y, area);
}

void on_glarea_drag_update_primary(GtkGestureDrag *gesture,
                           double          delta_x,
                           double          delta_y,
                           GtkWidget      *area) {

   graphics_info_t g;
   // Hack for mac. Needs more thought.
   // g.on_glarea_drag_update_primary(gesture, delta_x, delta_y, area);
   g.on_glarea_drag_update_secondary(gesture, delta_x, delta_y, area);

}

void on_glarea_drag_end_primary(GtkGestureDrag *gesture,
                                double          x,
                                double          y,
                                GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_end_primary(gesture, x, y, area);
}


void on_glarea_drag_begin_secondary(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_begin_secondary(gesture, x, y, area);
}

void on_glarea_drag_update_secondary(GtkGestureDrag *gesture,
                                     double          delta_x,
                                     double          delta_y,
                                     GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_update_secondary(gesture, delta_x, delta_y, area);
}

void on_glarea_drag_end_secondary(GtkGestureDrag *gesture,
                                  double          x,
                                  double          y,
                                  GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_end_secondary(gesture, x, y, area);
}



void on_glarea_drag_begin_middle(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_begin_middle(gesture, x, y, area);
}

void on_glarea_drag_update_middle(GtkGestureDrag *gesture,
                                  double          delta_x,
                                  double          delta_y,
                                  GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_update_middle(gesture, delta_x, delta_y, area);
}

void on_glarea_drag_end_middle(GtkGestureDrag *gesture,
                               double          x,
                               double          y,
                               GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_end_middle(gesture, x, y, area);
}


gboolean
on_glarea_key_controller_key_pressed(GtkEventControllerKey *controller,
                                     guint                  keyval,
                                     guint                  keycode,
                                     guint                  modifiers,
                                     GtkButton             *button) {

   graphics_info_t g;
   // allow other controllers to act (say TAB has been pressed)
   gboolean handled = g.on_glarea_key_controller_key_pressed(controller, keyval, keycode, modifiers);
   return gboolean(handled);
}


void
on_glarea_key_controller_key_released(GtkEventControllerKey *controller,
                                      guint                  keyval,
                                      guint                  keycode,
                                      guint                  modifiers,
                                      GtkButton             *button) {

   graphics_info_t g;
   g.on_glarea_key_controller_key_released(controller, keyval, keycode, modifiers);

}


void
on_glarea_click(GtkGestureClick* click_gesture,
                gint n_press,
                gdouble x,
                gdouble y,
                gpointer user_data) {

   // std::cout << "click gesture " << std::endl;
   graphics_info_t g;
   g.on_glarea_click(click_gesture, n_press, x, y, user_data);

}

void
on_glarea_scrolled(GtkEventControllerScroll *controller,
                   double                    dx,
                   double                    dy,
                   gpointer                  user_data) {

   graphics_info_t g;
   g.on_glarea_scrolled(controller, dx, dy, user_data);

}


void setup_gestures(GtkWidget *glarea) {

      std::cout << "================= setting up GTK4 style event controlllers ====================" << std::endl;

      GtkEventController *key_controller = gtk_event_controller_key_new();

      g_signal_connect(key_controller, "key-pressed",  G_CALLBACK(on_glarea_key_controller_key_pressed),  glarea);
      g_signal_connect(key_controller, "key-released", G_CALLBACK(on_glarea_key_controller_key_released), glarea);
      gtk_widget_add_controller(GTK_WIDGET(glarea), key_controller);

      GtkGesture *drag_controller_secondary = gtk_gesture_drag_new();
      GtkGesture *drag_controller_primary   = gtk_gesture_drag_new();
      GtkGesture *drag_controller_middle    = gtk_gesture_drag_new();
      GtkGesture *click_controller          = gtk_gesture_click_new();

      GtkEventControllerScrollFlags scroll_flags = GTK_EVENT_CONTROLLER_SCROLL_VERTICAL;
      GtkEventController *scroll_controller = gtk_event_controller_scroll_new(scroll_flags);

      // #ifdef __APPLE__
      //    mouse_view_rotate_button_mask = GDK_BUTTON1_MASK; // GDK_BUTTON_PRIMARY
      //    mouse_pick_button_mask        = GDK_BUTTON1_MASK; // GDK_BUTTON_PRIMARY
      // #endif

      gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_primary), GDK_BUTTON_PRIMARY);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER (drag_controller_primary));
      g_signal_connect(drag_controller_primary, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_primary),  glarea);
      g_signal_connect(drag_controller_primary, "drag-update", G_CALLBACK(on_glarea_drag_update_primary), glarea);
      g_signal_connect(drag_controller_primary, "drag-end",    G_CALLBACK(on_glarea_drag_end_primary),    glarea);

      gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_secondary), GDK_BUTTON_SECONDARY);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER (drag_controller_secondary));
      g_signal_connect(drag_controller_secondary, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_secondary),  glarea);
      g_signal_connect(drag_controller_secondary, "drag-update", G_CALLBACK(on_glarea_drag_update_secondary), glarea);
      g_signal_connect(drag_controller_secondary, "drag-end",    G_CALLBACK(on_glarea_drag_end_secondary),    glarea);

      gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_middle), GDK_BUTTON_MIDDLE);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER (drag_controller_middle));
      g_signal_connect(drag_controller_middle, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_middle),  glarea);
      g_signal_connect(drag_controller_middle, "drag-update", G_CALLBACK(on_glarea_drag_update_middle), glarea);
      g_signal_connect(drag_controller_middle, "drag-end",    G_CALLBACK(on_glarea_drag_end_middle),    glarea);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(click_controller));
      g_signal_connect(click_controller, "pressed",  G_CALLBACK(on_glarea_click),  glarea);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(scroll_controller));
      g_signal_connect(scroll_controller, "scroll",  G_CALLBACK(on_glarea_scrolled),  glarea);

}

void
new_startup_application_activate(GtkApplication *application,
                                 G_GNUC_UNUSED gpointer user_data) {

   GtkBuilder *builder = gtk_builder_new();
   if (GTK_IS_BUILDER(builder)) {
   } else {
      std::cout << "in new_startup_application_activate() builder was NOT a builder" << std::endl;
   }

   std::string dir = coot::package_data_dir();
   std::string dir_glade = coot::util::append_dir_dir(dir, "glade");
   std::string glade_file_name = "coot-gtk4.ui";
   std::string glade_file_full = coot::util::append_dir_file(dir_glade, glade_file_name);
   if (coot::file_exists(glade_file_name))
      glade_file_full = glade_file_name;

   GError* error = NULL;
   gboolean status = gtk_builder_add_from_file(builder, glade_file_full.c_str(), &error);
   if (status == FALSE) {
      std::cout << "ERROR:: Failure to read or parse " << glade_file_full << std::endl;
      std::cout << error->message << std::endl;
      exit(0);
   }

   GtkWidget *sb = GTK_WIDGET(gtk_builder_get_object(builder, "main_window_statusbar"));
   graphics_info_t::statusbar = sb;

   GtkWidget *app_window = gtk_application_window_new(application);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_title(GTK_WINDOW(app_window), "Coot App Main Window");
   graphics_info_t::set_main_window(app_window);

   guint id = gtk_application_window_get_id(GTK_APPLICATION_WINDOW(app_window));
   std::cout << "debug:: new_startup_application_activate(): Window id: " << id << std::endl;

   graphics_info_t g;
   g.set_gtkbuilder(builder);

   //GMenu *menu = create_menu_by_hand(application);
   GMenu *menubar = G_MENU(g.get_gobject_from_builder("menubar"));
   gtk_application_set_menubar(application, G_MENU_MODEL(menubar));
   gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(app_window), TRUE);

   // GtkWidget *graphics_hbox = widget_from_builder("crows_graphics_hbox", builder);
   // GtkWidget *main_window   = widget_from_builder("crows_main_window",   builder);
   GtkWidget *graphics_hbox = widget_from_builder("main_window_hbox", builder);
   GtkWidget *graphics_vbox = widget_from_builder("main_window_vbox", builder);
   // GObject *menubar  = g.get_gobject_from_builder("main_window_menubar");

   gtk_window_set_child(GTK_WINDOW(app_window), graphics_vbox);

   gtk_window_present(GTK_WINDOW(app_window));
   // gtk_widget_show(window);


   GtkWidget *gl_area = new_startup_create_glarea_widget();
   graphics_info_t::glareas.push_back(gl_area);
   gtk_widget_show(gl_area);
   // gtk_box_prepend(GTK_BOX(graphics_hbox), gl_area); // crows
   gtk_box_prepend(GTK_BOX(graphics_hbox), gl_area);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_default_size(GTK_WINDOW(app_window), 300, 300);
   gtk_widget_set_size_request(gl_area, 700, 400); // bigger than the window size - for testing.
   gtk_widget_show(app_window);

   setup_gestures(gl_area);

   create_actions(application);

   // load_tutorial_model_and_data();
}

// move these to the top.
void setup_symm_lib();
void check_reference_structures_dir();

int new_startup(int argc, char **argv) {

   graphics_info_t graphics_info;
   setup_symm_lib();
   check_reference_structures_dir();
   graphics_info.init();
   gtk_init();

   // set this by parsing the command line arguments
   graphics_info.use_graphics_interface_flag = true;

   // setup_python(argc, argv) needs to called after GTK has started - because it depends on gtk

   g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);

   GError *error = NULL;
   GtkApplication *app = gtk_application_new ("org.coot.crows", G_APPLICATION_FLAGS_NONE);
   g_signal_connect(app, "activate", G_CALLBACK(new_startup_application_activate), NULL);
   g_application_register(G_APPLICATION(app), NULL, &error);
   int status = g_application_run (G_APPLICATION (app), argc, argv);
   std::cout << "--- g_application_run() returns with status " << status << std::endl;
   g_object_unref (app);
   return status;
}