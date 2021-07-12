#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#include "compat/coot-sysdep.h"

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "c-interface.h" // for update_go_to_atom_from_current_position()
#include "globjects.h"
#include "graphics-info.h"

#include "draw.hh"
#include "draw-2.hh"
#include "framebuffer.hh"

#include "text-rendering-utils.hh"
#include "cc-interface-scripting.hh"
#include "cylinder-with-rotation-translation.hh"

#include "screendump-tga.hh"

enum {VIEW_CENTRAL_CUBE, ORIGIN_CUBE};


float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
      // positions   // texCoords
      -1.0f,  1.0f,  0.0f, 1.0f,
      -1.0f, -1.0f,  0.0f, 0.0f,
       1.0f, -1.0f,  1.0f, 0.0f,

      -1.0f,  1.0f,  0.0f, 1.0f,
       1.0f, -1.0f,  1.0f, 0.0f,
       1.0f,  1.0f,  1.0f, 1.0f
};


void
graphics_info_t::init_screen_quads() {

   graphics_info_t::shader_for_screen.Use();
   // screen quad VAO
   unsigned int quadVBO;
   glGenVertexArrays(1, &graphics_info_t::screen_quad_vertex_array_id);
   glBindVertexArray(graphics_info_t::screen_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   GLenum err = glGetError();
   if (err) std::cout << "init_screen_quads() err is " << err << std::endl;

}
void
graphics_info_t::init_blur_quads() {

   graphics_info_t::shader_for_blur.Use();
   // screen quad VAO
   unsigned int quadVBO;
   glGenVertexArrays(1, &graphics_info_t::blur_quad_vertex_array_id);
   glBindVertexArray(graphics_info_t::blur_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   GLenum err = glGetError();
   if (err) std::cout << "init_blur_quads() err is " << err << std::endl;

}

// void graphics_info_t::init_central_cube();

void
graphics_info_t::init_buffers() {
   init_central_cube();
   init_screen_quads();
   init_blur_quads();
}

void
graphics_info_t::init_central_cube() {

   float positions[24] = {
                          -0.5,  -0.5, -0.5,
                          -0.5,  -0.5,  0.5,
                          -0.5,   0.5, -0.5,
                          -0.5,   0.5,  0.5,
                           0.5,  -0.5, -0.5,
                           0.5,  -0.5,  0.5,
                           0.5,   0.5, -0.5,
                           0.5,   0.5,  0.5
   };

   glUseProgram(graphics_info_t::shader_for_central_cube.get_program_id());
   GLenum err = glGetError();
   if (err) std::cout << "init_central_cube() glUseProgram() err is " << err << std::endl;

   // number of lines * 2:
   unsigned int indices[24] { 0,1, 1,5, 5,4, 4,0, 2,3, 3,7, 7,6, 6,2, 0,2, 1,3, 5,7, 4,6 };

   // GLuint VertexArrayID;
   glGenVertexArrays(1, &graphics_info_t::central_cube_vertexarray_id);
   glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);

   // GLuint vertexbuffer;
   glGenBuffers(1, &graphics_info_t::central_cube_array_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::central_cube_array_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 24, &positions[0], GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

   // unsigned int ibo;
   glGenBuffers(1, &graphics_info_t::central_cube_index_buffer_id);
   err = glGetError();
   if (err) std::cout << "init_central_cube() index glGenBuffers() err is " << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::central_cube_index_buffer_id);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 24, &indices[0], GL_STATIC_DRAW);
   err = glGetError();
   if (err) std::cout << "init_central_cube() glBufferData() err is " << err << std::endl;
   glBindVertexArray(0);

}

void
graphics_info_t::init_hud_text() {

   graphics_info_t g;
   g.load_freetype_font_textures();
   glUseProgram(g.shader_for_hud_text.get_program_id());
   GLenum err = glGetError();
   if (err) std::cout << "init_hud_text() glUseProgram() err is " << err << std::endl;
   glGenVertexArrays(1, &graphics_info_t::hud_text_vertexarray_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glGenVertexArrays() err is " << err << std::endl;
   glBindVertexArray(graphics_info_t::hud_text_vertexarray_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glBindVertexArray() err is " << err << std::endl;
   glGenBuffers(1, &graphics_info_t::hud_text_array_buffer_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glGenBuffers() err is " << err << std::endl;
   glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::hud_text_array_buffer_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glBindBuffer() err is " << err << std::endl;
   glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
   err = glGetError(); if (err) std::cout << "init_hud_text() glBufferData() err is " << err << std::endl;
   glEnableVertexAttribArray(0);
   err = glGetError(); if (err) std::cout << "init_hud_text() glEnableVertexAttribArray() err is " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindVertexArray(0);
}


void
graphics_info_t::handle_delete_item_curor_change(GtkWidget *widget) {

   if (delete_item_widget) {
      if (delete_item_water) {
         graphics_info_t g;
	 pick_info naii = g.atom_pick_gtk3(false);
         GdkDisplay *display = gdk_display_get_default();
         GdkWindow *window = 0;
#if (GTK_MAJOR_VERSION < 4)
         window = gtk_widget_get_window(GTK_WIDGET(widget));
#endif
         if (window) {
            GdkCursor *current_cursor = gdk_window_get_cursor(window);
            // std::cout << "current cursor " << gdk_cursor_get_cursor_type(current_cursor) << std::endl;
            if (naii.success == GL_TRUE) {
               int imol = naii.imol;
               molecule_class_info_t &m = graphics_info_t::molecules[imol];
               std::string res_name = m.atom_sel.atom_selection[naii.atom_index]->GetResName();
               if (res_name == "HOH") {
                  GdkCursor *c = gdk_cursor_new_from_name (display, "crosshair");
                  // std::cout << "crosshair type " << gdk_cursor_get_cursor_type(c) << std::endl;
                  gdk_window_set_cursor(window, c);
               } else {
                  GdkCursor *c = gdk_cursor_new_from_name (display, "not-allowed");
                  // std::cout << "not-allowed type " << gdk_cursor_get_cursor_type(c) << std::endl;
                  gdk_window_set_cursor(window, c);
               }
            } else {
               GdkCursor *c = gdk_cursor_new_from_name (display, "not-allowed");
               // std::cout << "not-allowed type " << gdk_cursor_get_cursor_type(c) << std::endl;
               gdk_window_set_cursor(window, c);
            }
         }
      }
   }
}


// static
void
graphics_info_t::mouse_zoom(double delta_x, double delta_y) {

   // Zooming
   double fx = 1.0 + delta_x/300.0;
   double fy = 1.0 + delta_y/300.0;
   if (fx > 0.0) graphics_info_t::zoom /= fx;
   if (fy > 0.0) graphics_info_t::zoom /= fy;
   if (false)
      std::cout << "zooming with perspective_projection_flag "
                << graphics_info_t::perspective_projection_flag
                << " " << graphics_info_t::zoom << std::endl;
   if (! graphics_info_t::perspective_projection_flag) {
      // std::cout << "now zoom: " << g.zoom << std::endl;
   } else {
      // Move the eye towards the rotation centre (don't move the rotation centre)
      if (fabs(delta_y) > fabs(delta_x))
         delta_x = delta_y;
      float sf = 1.0 - delta_x * 0.003;
      graphics_info_t::eye_position.z *= sf;

      { // own graphics_info_t function - c.f. adjust clipping
         double  l = graphics_info_t::eye_position.z;
         double zf = graphics_info_t::screen_z_far_perspective;
         double zn = graphics_info_t::screen_z_near_perspective;

         graphics_info_t::screen_z_near_perspective *= sf;
         graphics_info_t::screen_z_far_perspective  *= sf;

         float screen_z_near_perspective_limit = l * 0.95;
         float screen_z_far_perspective_limit  = l * 1.05;
         if (graphics_info_t::screen_z_near_perspective < 2.0)
            graphics_info_t::screen_z_near_perspective = 2.0;
         if (graphics_info_t::screen_z_far_perspective > 1000.0)
            graphics_info_t::screen_z_far_perspective = 1000.0;

         if (graphics_info_t::screen_z_near_perspective > screen_z_near_perspective_limit)
            graphics_info_t::screen_z_near_perspective = screen_z_near_perspective_limit;
         if (graphics_info_t::screen_z_far_perspective < screen_z_far_perspective_limit)
            graphics_info_t::screen_z_far_perspective = screen_z_far_perspective_limit;
         if (false)
            std::cout << "on_glarea_motion_notify(): debug l: " << l << " post-manip: "
                      << graphics_info_t::screen_z_near_perspective << " "
                      << graphics_info_t::screen_z_far_perspective << std::endl;
      }
   }
   graphics_draw(); // or should this be called by the function that calls this function?
}

// static
void
graphics_info_t::scroll_zoom(int direction) {

   // c.f. mouse_zoom() it was copied from there. I am not sure that I like this yet.
   // and don't want to refactor if I'm going to dump either this or that later.

   // scroll up mean direction -1

   // Zooming
   double delta_x = 15.0;
   if (direction == 1) delta_x = -delta_x;
   double fx = 1.0 + delta_x/300.0;
   if (fx > 0.0) graphics_info_t::zoom /= fx;
   if (false)
      std::cout << "zooming with perspective_projection_flag "
                << graphics_info_t::perspective_projection_flag
                << " " << graphics_info_t::zoom << std::endl;
   if (! graphics_info_t::perspective_projection_flag) {
      // std::cout << "now zoom: " << g.zoom << std::endl;
   } else {
      // Move the eye towards the rotation centre (don't move the rotation centre)
      float sf = 1.0 - delta_x * 0.003;
      graphics_info_t::eye_position.z *= sf;

      { // own graphics_info_t function - c.f. adjust clipping
         double  l = graphics_info_t::eye_position.z;

         graphics_info_t::screen_z_near_perspective *= sf;
         graphics_info_t::screen_z_far_perspective  *= sf;

         float screen_z_near_perspective_limit = l * 0.95;
         float screen_z_far_perspective_limit  = l * 1.05;
         if (graphics_info_t::screen_z_near_perspective < 2.0)
            graphics_info_t::screen_z_near_perspective = 2.0;
         if (graphics_info_t::screen_z_far_perspective > 1000.0)
            graphics_info_t::screen_z_far_perspective = 1000.0;

         if (graphics_info_t::screen_z_near_perspective > screen_z_near_perspective_limit)
            graphics_info_t::screen_z_near_perspective = screen_z_near_perspective_limit;
         if (graphics_info_t::screen_z_far_perspective < screen_z_far_perspective_limit)
            graphics_info_t::screen_z_far_perspective = screen_z_far_perspective_limit;
         if (false)
            std::cout << "on_glarea_motion_notify(): debug l: " << l << " post-manip: "
                      << graphics_info_t::screen_z_near_perspective << " "
                      << graphics_info_t::screen_z_far_perspective << std::endl;
      }
   }
   graphics_draw(); // or should this be called by the function that calls this function?
}




// static
glm::mat4
graphics_info_t::get_model_view_matrix() {

   // where is this used?

   if (! perspective_projection_flag) {
      return glm::toMat4(glm_quat);
   } else {
      return glm::mat4(1.0);
   }
}

// static
glm::mat4
graphics_info_t::get_molecule_mvp(bool debug_matrices) {

   // presumes that we are in the correct programID

   float w = static_cast<float>(graphics_info_t::graphics_x_size);
   float h = static_cast<float>(graphics_info_t::graphics_y_size);
   float screen_ratio = static_cast<float>(w)/static_cast<float>(h);

   // I don't think that the quaternion belongs to the model matrix, it should be
   // part of the view matrix I think.
   // Yes. That's right.
   glm::mat4 model_matrix = glm::mat4(1.0);

   float z = graphics_info_t::zoom * 0.04;
   glm::vec3 sc(z,z,z);

   GLfloat near = -0.1 * zoom * clipping_front;
   GLfloat far  =  0.3 * zoom * clipping_back;

   if (false)
      std::cout << "near " << near << " far " << far << " clipping front "
                << clipping_front << " back " << clipping_back << std::endl;

   float sr = screen_ratio;
   glm::mat4 projection_matrix = glm::ortho(-0.3f*zoom*sr, 0.3f*zoom*sr,
                                            -0.3f*zoom,    0.3f*zoom,
                                            near, far);


   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   // std::cout << "rotation centre " << glm::to_string(rc) << std::endl;
   glm::mat4 view_matrix = glm::toMat4(graphics_info_t::glm_quat);

   view_matrix = glm::translate(view_matrix, -rc);
   // view_matrix = glm::scale(view_matrix, reverse_z); causes weirdness - not sure about handedness
   glm::mat4 mvp = projection_matrix * view_matrix * model_matrix;

   if (graphics_info_t::perspective_projection_flag) {

      // for fun/testing
      // turn off view scaling when tinkering with this?
      // there should not be a concept of "zoom" with perspective view, just translation
      // along screen-Z.

      glm::mat4 trackball_matrix = glm::toMat4(graphics_info_t::glm_quat);

      glm::vec3 ep = eye_position; // in view space i.e. (0,0,z) (z = 40, say)
      glm::vec3 up(0,1,0);
      glm::vec3 origin(0,0,0);

      model_matrix = glm::mat4(1.0);
      model_matrix = glm::translate(model_matrix, -rc);
      model_matrix = trackball_matrix * model_matrix;

      view_matrix = glm::lookAt(ep, origin, up);

      float fov = 50.0/zoom; // degrees
      if (fov > 50.0) fov = 50.0;
      fov = 35.0;
      fov = 30.0;

      glm::mat4 projection_matrix_persp = glm::perspective(glm::radians(fov),
                                                           screen_ratio,
                                                           screen_z_near_perspective,
                                                           screen_z_far_perspective);
      projection_matrix = projection_matrix_persp; // for debugging below
      mvp = projection_matrix_persp * view_matrix * model_matrix;
   }

   // bool debug_matrices = false;
   if (debug_matrices) {
      std::cout << "model, view, projection, mvp" << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(model_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(view_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(projection_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(mvp) << std::endl;
   }
   return mvp;
}

// can we work out the eye position without needing to unproject? (because that depends
// on get_molecule_mvp()...
//
glm::vec3
graphics_info_t::get_world_space_eye_position() {

   if (! graphics_info_t::perspective_projection_flag) {

      // orthograph eye position is inferred from centre position
      // and zoom and view rotation.

      // does this work? How can I tell?

      glm::vec3 test_vector_1(0.0, 0.0, 1.0);
      glm::vec3 test_vector_2(1.0, 1.0, 0.0);

      glm::mat4 vr = get_view_rotation();
      glm::vec4 rot_test_vector_1 = glm::vec4(test_vector_1, 1.0) * vr;
      glm::vec4 rot_test_vector_2 = glm::vec4(test_vector_2, 1.0) * vr;

      glm::vec3 ep = graphics_info_t::zoom * glm::vec3(rot_test_vector_1);
      glm::vec3 rc = graphics_info_t::get_rotation_centre();
      ep += rc;

      return ep;

   } else {

      // the eye_position is in view-coordinates is stored directly and is
      // by default and often (always?) (0,0,40).

      // I need to convert that to world coordinates and then rotate
      // and translate the world according to rotation centre and mouse-based
      // quaternion

      glm::vec3 ep = eye_position;
      glm::vec4 ep_4(ep, 1.0);
      glm::vec3 up(0,1,0);
      glm::vec3 origin(0,0,0);

      glm::mat4 trackball_matrix = glm::toMat4(graphics_info_t::glm_quat);
      glm::vec3 rc = graphics_info_t::get_rotation_centre();
      glm::mat4 model_matrix = glm::mat4(1.0);
      model_matrix = glm::translate(model_matrix, -rc);
      model_matrix = trackball_matrix * model_matrix;
      glm::mat4 model_inv = glm::inverse(model_matrix);

      glm::mat4 view_matrix = glm::lookAt(ep, origin, up);
      glm::mat4 view_inv = glm::inverse(view_matrix);
      glm::vec4 ep_world_4 = model_inv * ep_4; // * view_inv;

      if (false) {
         std::cout << "model_inv " << glm::to_string(model_inv) << std::endl;
         std::cout << "view_inv  " << glm::to_string( view_inv) << std::endl;
      }

      glm::vec3 ep_world(ep_world_4);

      return ep_world;
   }

}

glm::vec4
graphics_info_t::unproject(float z) {

   // z is 1 and -1 for front and back (or vice verse).

   GtkAllocation allocation;
   gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
   float w = allocation.width;
   float h = allocation.height;
   graphics_info_t g;
   float mouseX = 2.0 *    g.GetMouseBeginX()/w  - 1.0f;
   float mouseY = 2.0 * (h-g.GetMouseBeginY())/h - 1.0f;
   std::cout << "debug in new_unproject widget w and h " << w << " " << h << std::endl;
   std::cout << "debug in new_unproject mouse x and y widget  "
             << g.GetMouseBeginX() << " "
             << g.GetMouseBeginY() << std::endl;
   std::cout << "debug in new_unproject mouse x and y GL      " << mouseX << " " << mouseY << std::endl;
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   float real_y = - mouseY; // in range -1 -> 1
   glm::vec4 screenPos_f = glm::vec4(mouseX, real_y, z, 1.0f);
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   std::cout << "debug in new_unproject() screen_pos " << glm::to_string(screenPos_f) << std::endl;
   std::cout << "debug in new_unproject() world_pos " << glm::to_string(worldPos_f) << std::endl;
   return worldPos_f;

}


glm::mat4
graphics_info_t::get_view_rotation() {

   // need to be in the correct program (well, the model-drawing part)

   return glm::toMat4(graphics_info_t::glm_quat);

}


void
graphics_info_t::setup_map_uniforms(Shader *shader_p,
                                    const glm::mat4 &mvp,
                                    const glm::mat4 &view_rotation,
                                    float density_surface_opacity) {

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   GLenum err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniformMatrix4fv() mvp " << err << std::endl;
   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniformMatrix4fv() vr  " << err << std::endl;

   GLuint background_colour_uniform_location = shader_p->background_colour_uniform_location;
   glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
   glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
   err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniform4fv() for bg  " << err << std::endl;

   // GLuint opacity_uniform_location = shader.map_opacity_uniform_location;
   shader_p->set_float_for_uniform("map_opacity", density_surface_opacity);
   err = glGetError(); if (err) std::cout << "   setup_map_uniforms() glUniformf() for opacity "
                                          << err << std::endl;

   GLuint eye_position_uniform_location = shader_p->eye_position_uniform_location;
   glm::vec4 ep = glm::vec4(get_world_space_eye_position(), 1.0);
   glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
   err = glGetError(); if (err) std::cout << "   setup_map_uniforms() glUniform4fv() for eye position "
                                          << err << std::endl;

}


void
graphics_info_t::draw_map_molecules(bool draw_transparent_maps) {

   // run through this molecule loop twice - for opaque then transparent maps
   // first, a block that decides if we need to do anything.

   bool needs_blend_reset = false;

   //

   bool cosine_dependent_map_opacity = true;

   unsigned int n_transparent_maps = 0;
   unsigned int n_maps_to_draw = 0;
   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
      const molecule_class_info_t &m = graphics_info_t::molecules[ii];
      if (draw_transparent_maps) {
         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         if (! m.draw_it_for_map) continue;
         if (! m.is_an_opaque_map()) {
            n_transparent_maps++;
            n_maps_to_draw += 1;
         }
      } else {
         if (m.is_an_opaque_map()) {
            if (m.draw_it_for_map)
               n_maps_to_draw += 1;
         }
      }
   }

   if (n_maps_to_draw == 0) return;

   if (n_transparent_maps > 0) {
      needs_blend_reset = true;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   if (cosine_dependent_map_opacity) {
      needs_blend_reset = true;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   GLenum err = glGetError();
   if (err) std::cout << "gtk3_draw_map_molecules() A " << err << std::endl;

   if (!draw_transparent_maps || n_transparent_maps > 0) {

      glLineWidth(map_line_width * framebuffer_scale);
      err = glGetError();
      if (err) std::cout << "gtk3_draw_map_molecules() glLineWidth " << err << std::endl;

      shader_for_maps.Use();

      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation = get_view_rotation();

      glEnable(GL_DEPTH_TEST); // this needs to be in the draw loop!?
      glDepthFunc(GL_LESS);

      Shader &shader = shader_for_maps;

      glm::vec4 ep(get_world_space_eye_position(), 1.0);
      glm::vec3 ep3 = ep/ep.w;

      for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {

         molecule_class_info_t &m = graphics_info_t::molecules[ii]; // not const because shader changes
         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         if (! m.draw_it_for_map) continue;
         if (draw_transparent_maps) {
            if (m.is_an_opaque_map())
               continue; // not this round
         } else {
            // only draw (completely) opaque (that's what the question means)
            if (! m.is_an_opaque_map())
               continue;
         }

         if (m.n_vertices_for_map_VertexArray > 0) {

            err = glGetError();
            if (err) std::cout << "draw_map_molecules() --- draw map loop start --- error "
                               << std::endl;

            bool draw_with_lines = true;
            if (!m.draw_it_for_map_standard_lines) draw_with_lines = false;

            //glUniform1i(shader.is_perspective_projection_uniform_location,
            // graphics_info_t::perspective_projection_flag);
            shader.set_bool_for_uniform("is_perspective_projection", perspective_projection_flag);
            err = glGetError(); if (err) std::cout << "   draw_map_molecules() error B " << std::endl;

            shader.set_bool_for_uniform("do_depth_fog", graphics_info_t::shader_do_depth_fog_flag);
            shader.set_bool_for_uniform("do_diffuse_lighting", true);
            shader.set_float_for_uniform("shininess", m.shader_shininess);
            shader.set_float_for_uniform("specular_strength", m.shader_specular_strength);

            // --- lights ----

            std::map<unsigned int, lights_info_t>::const_iterator it; // iterate over the lights map
            for (it=lights.begin(); it!=lights.end(); it++) {
               unsigned int light_idx = it->first;
               const lights_info_t &light = it->second;
               shader.setup_light(light_idx, light, view_rotation);
            }

            // --- material ---

            Material &material = m.material_for_maps;
            shader.set_vec4_for_uniform( "material.ambient",   material.ambient);
            shader.set_vec4_for_uniform( "material.diffuse",   material.diffuse);
            shader.set_vec4_for_uniform( "material.specular",  material.specular);
            shader.set_float_for_uniform("material.shininess", material.shininess);
            shader.set_float_for_uniform("material.specular_strength", material.specular_strength);

            // --- background ---

            GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
            glm::vec4 bgc(background_colour, 1.0);
            glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniform4fv() for bg  " << err << std::endl;

            // --- fresnel ---

            shader.set_bool_for_uniform("do_fresnel",     m.fresnel_settings.state);
            shader.set_float_for_uniform("fresnel_bias",  m.fresnel_settings.bias);
            shader.set_float_for_uniform("fresnel_scale", m.fresnel_settings.scale);
            shader.set_float_for_uniform("fresnel_power", m.fresnel_settings.power);
            shader.set_vec4_for_uniform("fresnel_colour", m.fresnel_settings.colour);

            // --- draw ---

            if (draw_with_lines) {
               // I don't see why this is needed - but it is.
               if (! m.is_an_opaque_map())
                  glEnable(GL_BLEND);

               glBindVertexArray(m.m_VertexArrayID_for_map);
               err = glGetError();
               if (err) std::cout << "ERROR:: draw_map_molecules() glBindVertexArray() "
                                  << m.m_VertexArrayID_for_map
                                  << " with GL err " << err << std::endl;
               if (err) {
                  // no point in continuing
                  std::cout << "### Catastrophic failure in draw_map_molecules() returning now " << std::endl;
                  return;
               }

               glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_map_lines_ID);

               setup_map_uniforms(&shader, mvp, view_rotation, m.density_surface_opacity);
               glDrawElements(GL_LINES, m.n_vertices_for_map_VertexArray,
                              GL_UNSIGNED_INT, nullptr);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glDrawElements() n_vertices: "
                                  << m.n_vertices_for_map_VertexArray
                                  << " with GL err " << err << std::endl;
            }

            if (!draw_with_lines) { // draw as a solid object
               if (false)
                  std::cout << "   draw_map_molecules(): imol " << ii
                            << " array_id and n_vertices_for_VertexArray: "
                            << m.m_VertexArrayID_for_map << " "
                            << m.n_indices_for_triangles
                            << std::endl;

               if (! m.is_an_opaque_map()) {
                  // sort the triangles
                  clipper::Coord_orth eye_pos_co(ep.x, ep.y, ep.z);
                  m.sort_map_triangles(eye_pos_co);
               }

               glBindVertexArray(m.m_VertexArrayID_for_map);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glBindVertexArray() "
                                  << m.m_VertexArrayID_for_map << " with GL err " << err << std::endl;
               glEnable(GL_BLEND);
               glBindBuffer(GL_ARRAY_BUFFER,         m.m_VertexBufferID);
               glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_map_triangles_ID);

               glUniformMatrix4fv(graphics_info_t::shader_for_maps.mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() " << err << std::endl;
               glUniformMatrix4fv(shader.view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() " << err << std::endl;

               // opacity:
               // GLuint opacity_uniform_location = graphics_info_t::shader_for_maps.map_opacity_uniform_location;
               // float opacity = m.density_surface_opacity;
               // glUniform1f(opacity_uniform_location, opacity);
               shader.set_float_for_uniform("map_opacity", m.density_surface_opacity);
               err = glGetError(); if (err) std::cout << "   draw_map_molecules() glUniformf() for opacity "
                                                      << err << std::endl;
               // cosine_dependent_map_opacity
               shader.set_bool_for_uniform("cosine_dependent_map_opacity", cosine_dependent_map_opacity);

               GLuint eye_position_uniform_location = graphics_info_t::shader_for_maps.eye_position_uniform_location;
               glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));

               glDrawElements(GL_TRIANGLES, m.n_indices_for_triangles, GL_UNSIGNED_INT, nullptr);

               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glDrawElements() n_indices_for_triangles "
                                  << graphics_info_t::molecules[ii].n_indices_for_triangles
                                  << " with GL err " << err << std::endl;
            }
         }
      }


   }

   // to be clean we should use
   // glDisableVertexAttribArray(0);
   // here.
   // that would mean adding glEnableVertexAttribArray() for the attributes (position, normal, colour).
   // in the above block.

   if (needs_blend_reset) {
      glDisable(GL_BLEND);
   }
}

void
graphics_info_t::draw_model_molecules() {

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();

   glm::vec4 bgc(background_colour, 1.0);
   for (int ii=n_molecules()-1; ii>=0; ii--) {
      if (! is_valid_model_molecule(ii)) continue;
      molecule_class_info_t &m = molecules[ii];
      if (! m.draw_it) continue;
      Shader &shader_p = shader_for_model_as_meshes;
      m.draw_molecule_as_meshes(&shader_p, mvp, view_rotation, lights, eye_position, bgc, shader_do_depth_fog_flag);
      if (show_symmetry) {
         Shader &symm_shader_p = shader_for_symmetry_atoms_bond_lines;
         m.draw_symmetry(&symm_shader_p, mvp, view_rotation, lights, eye_position, bgc, shader_do_depth_fog_flag);
      }
   }

   // this block of code should be a member function of molecule_class_info_t

   Shader &shader = graphics_info_t::shader_for_models;
   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {

      molecule_class_info_t &m = graphics_info_t::molecules[ii];
      if (! graphics_info_t::is_valid_model_molecule(ii)) continue;
      if (! m.draw_it) continue;

      if (false)
         std::cout << "imol " << ii << " n_vertices_for_model_VertexArray "
                   << m.n_vertices_for_model_VertexArray << std::endl;
      if (m.n_vertices_for_model_VertexArray > 0) {

         glDisable(GL_BLEND); // stop semi-transparent bonds - but why do we have them?
         gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));

         shader.Use();
         GLuint err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glUseProgram() "
                                                       << err << std::endl;

         // glUniform1i(shader.is_perspective_projection_uniform_location,
         // graphics_info_t::perspective_projection_flag);
         shader.set_bool_for_uniform("is_perspective_projection", perspective_projection_flag);

         glBindVertexArray(m.m_VertexArray_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glBindVertexArray() "
                            << m.m_VertexArray_for_model_ID
                            << " with GL err " << err << std::endl;

         glBindBuffer(GL_ARRAY_BUFFER, m.m_VertexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() v " << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() i " << err << std::endl;

         GLuint mvp_location           = shader.mvp_uniform_location;
         GLuint view_rotation_location = shader.view_rotation_uniform_location;

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() pre mvp " << err << std::endl;
         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for mvp " << err << std::endl;
         glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for view_rotation "
                            << err << std::endl;

         // std::cout << glm::to_string(mvp) << std::endl;
         // std::cout << glm::to_string(view_rotation) << std::endl;

         GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
         glm::vec4 bgc(background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec3 ep(get_world_space_eye_position());
         glUniform3fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform3fv() for eye position " << err << std::endl;

         shader.set_bool_for_uniform("do_depth_fog", shader_do_depth_fog_flag);
         shader.set_bool_for_uniform("do_diffuse_lighting", true); // false for demo c.f. old style graphics

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() pre-lights glDrawElements() "
                            << shader.name << " with err " << err << std::endl;

         // --- lights ----

         std::map<unsigned int, lights_info_t>::const_iterator it; // iterate over the lights map
         for (it=lights.begin(); it!=lights.end(); it++) {
            unsigned int light_idx = it->first;
            const lights_info_t &light = it->second;
            shader.setup_light(light_idx, light, view_rotation);
         }
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() post-lights "
                            << shader.name << " with err " << err << std::endl;


         // --- material ---

         const Material &material = m.material_for_models;
         shader.set_vec4_for_uniform( "material.ambient",   material.ambient);
         shader.set_vec4_for_uniform( "material.diffuse",   material.diffuse);
         shader.set_vec4_for_uniform( "material.specular",  material.specular);
         shader.set_float_for_uniform("material.shininess", material.shininess);
         shader.set_float_for_uniform("material.specular_strength", material.specular_strength);

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() post-material "
                            << shader.name << " with err " << err << std::endl;

         // draw with the vertex count, not the index count.
         GLuint n_verts = graphics_info_t::molecules[ii].n_indices_for_model_triangles;

         // std::cout << "   Drawing " << n_verts << " model vertices" << std::endl;
         err = glGetError();
         if (err) std::cout << "   error pre draw_model_molecules() glDrawElements() " << shader.name
                            << " n_vertices " << n_verts << " with GL err " << err << std::endl;
         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glDrawElements() " << shader.name
                            << " n_vertices " << n_verts << " with GL err " << err << std::endl;

         draw_molecule_atom_labels(m, mvp, view_rotation);

      }
   }
}

void
graphics_info_t::draw_molecule_atom_labels(molecule_class_info_t &m,
                                           const glm::mat4 &mvp,
                                           const glm::mat4 &view_rotation) {

   // pass the glarea widget width and height.

   // put a triangle or square where the atom label should be, facing the camera
   // "billboarding"

   glm::vec4 label_colour(font_colour.red, font_colour.green, font_colour.blue, 1.0);

   if (false) { // test label
      // Put atom label test at 42, 9, 13
      glm::vec3 point(42, 9, 13);
      // point = glm::vec3(0,0,0);

      glm::vec4 projected_point_2 = mvp * glm::vec4(point, 1.0);
      if (true)
         std::cout << "projected point " << glm::to_string(projected_point_2) << std::endl;

      projected_point_2.x = 0.5 * (projected_point_2.x + 1.0f);
      projected_point_2.y = 0.5 * (projected_point_2.y + 1.0f);

      projected_point_2.x *= 900.0;
      projected_point_2.y *= 900.0;

      glm::vec3 pp(projected_point_2);

      glEnable(GL_DEPTH_TEST); // or we don't see the label. Either a blurred label or nothing.
      render_atom_label(shader_for_atom_labels, ". Test Label", pp, 1.0, label_colour);
   }

   int n_atoms_to_label = m.labelled_atom_index_list.size();
   if (n_atoms_to_label == 0) return;

   // maybe pass these?
   GtkAllocation allocation;
   GtkWidget *widget = graphics_info_t::glareas[0];
   if (! widget) return;
   gtk_widget_get_allocation(widget, &allocation);

   // this doesn't seem sensibly arranged.
   glm::vec3 unused(0,0,0);
   m.draw_atom_labels(brief_atom_labels_flag,
                      seg_ids_in_atom_labels_flag,
                      label_colour,
                      mvp, view_rotation, unused);

   glDisable(GL_BLEND);

}

void
graphics_info_t::draw_intermediate_atoms() { // draw_moving_atoms()

   // all these draw functions should be moved int graphics_info_t.

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();

   Shader &shader = graphics_info_t::shader_for_models;
   if (true) {
      const molecule_class_info_t &m = graphics_info_t::moving_atoms_molecule;

      if (m.n_vertices_for_model_VertexArray > 0) {

         glDisable(GL_BLEND); // stop semi-transparent bonds - but why do we have them?
         gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));

         glBindVertexArray(m.m_VertexArray_for_model_ID);
         GLenum err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glBindVertexArray() "
                            << m.m_VertexArray_for_model_ID
                            << " with GL err " << err << std::endl;

         shader.Use();
         err = glGetError();
         if (err) std::cout << "   error draw_intermediate_atoms() glUseProgram() "
                            << err << std::endl;

#if 0
         // should not be needed? - the VAO contains this information.  Needs testing.

         glBindBuffer(GL_ARRAY_BUFFER, m.m_VertexBuffer_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_intermediate_atoms() glBindBuffer() v "
                            << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_intermediate_atoms() glBindBuffer() i "
                            << err << std::endl;

#endif

         GLuint mvp_location           = shader.mvp_uniform_location;
         GLuint view_rotation_location = shader.view_rotation_uniform_location;

         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniformMatrix4fv() pre mvp "
                            << err << std::endl;
         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniformMatrix4fv() for mvp "
                            << err << std::endl;
         glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniformMatrix4fv() "
                            << "for view_rotation " << err << std::endl;

         GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
         glm::vec4 bgc(background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "error draw_model_molecules() glUniform4fv() for background "
                            << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec3 ep = get_world_space_eye_position();
         glUniform3fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniform4fv() for eye position "
                            << err << std::endl;

         shader.set_bool_for_uniform("do_depth_fog", shader_do_depth_fog_flag);

         // draw with the vertex count, not the index count.
         GLuint n_verts = m.n_indices_for_model_triangles;
         // std::cout << "   Drawing " << n_verts << " model vertices" << std::endl;
         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glDrawElements() "
                            << n_verts << " with GL err " << err << std::endl;

      }
   }
}


#include "Instanced-Markup-Mesh.hh"

void
graphics_info_t::setup_rama_balls() {

   rama_balls_mesh.setup_octasphere(2);
   rama_balls_mesh.setup_instancing_buffers(3000);
}

void
graphics_info_t::update_rama_balls(std::vector<Instanced_Markup_Mesh_attrib_t> *balls) {

   // so this function should be called set_rama_balls_new positions_and_colours

   // the calling function calls
   // rama_balls_mesh.update_instancing_buffers(balls) after this function

   auto rr = saved_dragged_refinement_results;

   balls->clear();

  glm::vec3 screen_up_dir(0.2, 0.3, 0.3);

   // std::cout << "update rama ball for " << rr.all_ramas.size() << " balls " << std::endl;
   for (unsigned int i=0; i<rr.all_ramas.size(); i++) {

      float rama_score = rr.all_ramas[i].distortion;

      const coot::atom_spec_t &spec_CA = rr.all_ramas[i].atom_spec_CA;
      mmdb::Atom *at = spec_CA.get_atom(moving_atoms_asc->mol);
      if (at) {

         float d = rr.all_ramas[i].distortion;
         glm::vec3 atom_position(at->x, at->y, at->z);
         glm::vec3 ball_position(rr.all_ramas[i].ball_pos_x,
                                 rr.all_ramas[i].ball_pos_y,
                                 rr.all_ramas[i].ball_pos_z);
         float size = 0.38;
         // std::cout << "debug d " << d << std::endl;
         float ra = hud_geometry_distortion_to_rotation_amount_rama(d);
         coot::colour_t cc(0.1, 0.9, 0.2);
         cc.rotate(ra);
         glm::vec4 col = cc.to_glm();
         Instanced_Markup_Mesh_attrib_t ball(col, ball_position, size);
         float d1 = d + 85.0;
         float d2 = - d1 * 0.016;
         if (d2 < 0.0) d2 = 0.0;
         if (d2 > 1.0) d2 = 1.0;
         ball.specular_strength = 0.01 + d2;
         ball.shininess = 0.9 + 155.0 * d2;
         balls->push_back(ball);
      }
   }
}



void
graphics_info_t::draw_intermediate_atoms_rama_balls() {

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   Shader &shader = graphics_info_t::shader_for_rama_balls;

   glm::mat4 mvp = get_molecule_mvp();
   glm::vec3 eye_position = get_world_space_eye_position();
   glm::mat4 view_rotation = get_view_rotation();
   glm::vec4 bg_col(background_colour, 1.0);
   bool do_depth_fog = true;
   // this is a bit ugly

   // the balls are updated after a refinement cycle has finished -
   // no need to do it here
   // graphics_info_t g;
   // std::vector<Instanced_Markup_Mesh_attrib_t> balls;
   // update_rama_balls(&balls);
   // rama_balls_mesh.update_instancing_buffers(balls);
   rama_balls_mesh.draw(&shader, mvp, view_rotation, lights, eye_position, bg_col, do_depth_fog);

}

void
graphics_info_t::setup_atom_pull_restraints_glsl() {

   // build the triangles for the cylinder and
   // set m_VertexArray_for_pull_restraints_ID
   //     m_VertexBuffer_for_pull_restraints_ID
   //     n_indices_for_atom_pull_triangles

   n_atom_pulls = 0; // now a class variable.
   for (std::size_t i=0; i<atom_pulls.size(); i++) {
      const atom_pull_info_t &atom_pull = atom_pulls[i];
      if (atom_pull.get_status()) {
         std::pair<bool, int> spec = atom_pull.find_spec(moving_atoms_asc->atom_selection,
                                                         moving_atoms_asc->n_selected_atoms);
         if (spec.first)
            n_atom_pulls++;
      }
   }

   if (n_atom_pulls > 0) {
      unsigned int n_slices = 10;
      unsigned int n_stacks = 2;
      n_vertices_for_atom_pull_restraints  = n_slices * (n_stacks +1) * n_atom_pulls;
      n_triangles_for_atom_pull_restraints = n_slices * n_stacks * 2  * n_atom_pulls;
      // add in the triangles for the arrow-head (lots of vertices at the arrow tip because from cylinder)
      unsigned int n_stacks_for_arrow_tip = 6;
      n_vertices_for_atom_pull_restraints  += n_stacks_for_arrow_tip * n_slices * 2 * n_atom_pulls;
      n_triangles_for_atom_pull_restraints += n_stacks_for_arrow_tip * n_slices * 2 * n_atom_pulls;
      // the indices of the vertices in the triangles (3 indices per triangle)
      unsigned int *flat_indices = new unsigned int[n_triangles_for_atom_pull_restraints * 3];
      unsigned int *flat_indices_start = flat_indices;
      unsigned int ifi = 0; // index into flat indices - running
      vertex_with_rotation_translation *vertices =
         new vertex_with_rotation_translation[n_vertices_for_atom_pull_restraints];
      vertex_with_rotation_translation *vertices_start = vertices;
      unsigned int iv = 0; // index into vertices - running

      for (std::size_t i=0; i<atom_pulls.size(); i++) {
         const atom_pull_info_t &atom_pull = atom_pulls[i];
         if (atom_pull.get_status()) {
            std::pair<bool, int> spec = atom_pull.find_spec(moving_atoms_asc->atom_selection,
                                                            moving_atoms_asc->n_selected_atoms);
            if (spec.first) {
               float arrow_head_length = 0.2;

               mmdb::Atom *at = moving_atoms_asc->atom_selection[spec.second];

               // coot::Cartesian pt_start_c(at->x, at->y, at->z);
               // coot::Cartesian pt_end_c(atom_pull.pos.x(), atom_pull.pos.y(), atom_pull.pos.z());
               // coot::Cartesian b = pt_end_c - pt_start_c;

               glm::vec3 pt_start_g(at->x, at->y, at->z);
               glm::vec3 pt_end_g(atom_pull.pos.x(), atom_pull.pos.y(), atom_pull.pos.z());
               glm::vec3 b = pt_end_g - pt_start_g;

               float bl_pull = glm::distance(b, glm::vec3(0,0,0));
               float bl = bl_pull - arrow_head_length;
               if (arrow_head_length > bl_pull)
                  arrow_head_length = bl_pull;
               glm::vec3 b_uv = glm::normalize(b);
               float bl_stick = bl;
               if (bl_stick < 0.0) bl_stick = 0.0;

               // coot::Cartesian meeting_point = pt_start_c + b_uv * bl_stick;
               // coot::CartesianPair pos_pair(pt_start_c, meeting_point);
               glm::vec3 meeting_point = pt_start_g + b_uv * bl_stick;
               std::pair<glm::vec3, glm::vec3> pos_pair(pt_start_g, meeting_point);
               float radius = 0.1;

               cylinder_with_rotation_translation c(pos_pair, radius, radius, bl, n_slices, n_stacks);
               for (std::size_t j=0; j<c.triangle_indices_vec.size(); j++) {
                  flat_indices[ifi  ] = c.triangle_indices_vec[j].point_id[0]+iv;
                  flat_indices[ifi+1] = c.triangle_indices_vec[j].point_id[1]+iv;
                  flat_indices[ifi+2] = c.triangle_indices_vec[j].point_id[2]+iv;
                  ifi += 3;
               }
               for (std::size_t j=0; j<c.vertices.size(); j++) {
                  // Use a constructor here when hmt code has been correctly integrated
                  vertices[iv] = c.vertices[j];
                  vertices[iv].colour = glm::vec4(0.8, 0.5, 0.3, 1.0);
                  iv++;
               }

               // coot::CartesianPair pp(pt_end_c, meeting_point);
               std::pair<glm::vec3, glm::vec3> pp(pt_end_g, meeting_point);
               // std::cout << "arrow-head: " << radius << " " << arrow_head_length << std::endl;
               cylinder_with_rotation_translation c_arrow_head(pp, 2.0 * radius, 0.0, arrow_head_length,
                                                               n_slices, n_stacks_for_arrow_tip);
               for (std::size_t j=0; j<c_arrow_head.triangle_indices_vec.size(); j++) {
                  if (ifi < n_triangles_for_atom_pull_restraints * 3) {
                     flat_indices[ifi  ] = c_arrow_head.triangle_indices_vec[j].point_id[0]+iv;
                     flat_indices[ifi+1] = c_arrow_head.triangle_indices_vec[j].point_id[1]+iv;
                     flat_indices[ifi+2] = c_arrow_head.triangle_indices_vec[j].point_id[2]+iv;
                  } else {
                     std::cout << "ERROR:: indexing for c_arrow_head "
                               << ifi << " " << n_triangles_for_atom_pull_restraints << std::endl;
                  }
                  ifi += 3;
               }
               for (std::size_t j=0; j<c_arrow_head.vertices.size(); j++) {
                  vertices[iv] = c_arrow_head.vertices[j];
                  vertices[iv].colour = glm::vec4(0.8,0.5,0.3,1.0);
                  iv++;
               }
            }
         }
      }

      // does this need to be done every time? I doubt it. Needs check.
      //
      glGenVertexArrays(1, &m_VertexArray_for_pull_restraints_ID);
      GLenum err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() A"
                          << " with GL err " << err << std::endl;
      glBindVertexArray(m_VertexArray_for_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() B"
                          << " with GL err " << err << std::endl;
      glGenBuffers(1, &m_VertexBuffer_for_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() C"
                          << " with GL err " << err << std::endl;
      glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer_for_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() D"
                          << " with GL err " << err << std::endl;
      GLuint n_bytes = sizeof(vertex_with_rotation_translation) * n_vertices_for_atom_pull_restraints;
      // maybe STATIC_DRAW, maybe not
      glBufferData(GL_ARRAY_BUFFER, n_bytes, vertices, GL_DYNAMIC_DRAW);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() E"
                          << " with GL err " << err << std::endl;


      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);
      glEnableVertexAttribArray(2);
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17c\n";
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(0 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17d\n";
      glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(1 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17e\n";
      glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17f\n";

      // translate position, 3, size 3 floats
      glEnableVertexAttribArray(3);

      // surely this (annd below) has been set-up already? -- CheckMe.
      glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(3 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17aa\n";

      // positions, 4, size 3 floats
      glEnableVertexAttribArray(4);
      err = glGetError(); if (err) std::cout << "GL error bonds 6\n";
      glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(4 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 7\n";

      //  normals, 5, size 3 floats
      glEnableVertexAttribArray(5);
      err = glGetError(); if (err) std::cout << "GL error bonds 11\n";
      glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(5 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 12\n";

      //  colours, 6, size 4 floats
      glEnableVertexAttribArray(6);
      err = glGetError(); if (err) std::cout << "GL error bonds 16\n";
      glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(6 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17\n";

      glGenBuffers(1, &m_IndexBuffer_for_atom_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() G"
                         << " with GL err " << err << std::endl;
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_atom_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() H"
                         << " with GL err " << err << std::endl;
      n_bytes = n_triangles_for_atom_pull_restraints * 3 * sizeof(unsigned int);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, flat_indices, GL_STATIC_DRAW);
      err = glGetError();
      if (err) std::cout << "setup_atom_pull_restraints_glsl() --end-- err " << err << std::endl;

      delete [] flat_indices;
      delete [] vertices;

   }
}


// static
void
graphics_info_t::draw_atom_pull_restraints() {

   // Note to self: do this first with standard (modern) OpenGL.
   //
   // Then do it again with instances. It will be faster to draw bonds and atoms that way.
   // Maybe density lines too.

   // don't draw this if intermediate atoms are not shown
   //
   if (! regularize_object_bonds_box.empty()) {
      if (!moving_atoms_asc) return;
      if (moving_atoms_asc->n_selected_atoms > 0) {

         if (n_atom_pulls > 0) { // class variable now.

            Shader &shader = shader_for_models;
            shader.Use();
            GLuint err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glUseProgram() "
                               << err << std::endl;

            glBindVertexArray(m_VertexArray_for_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindVertexArray()"
                               << " with GL err " << err << std::endl;
            glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer_for_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindBuffer()"
                               << " with GL err " << err << std::endl;

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_atom_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindBuffer() for index"
                               << " with GL err " << err << std::endl;

            // Uniforms - how are these not used!?
            glm::mat4 mvp = get_molecule_mvp();
            glm::mat4 view_rotation = get_view_rotation();
            GLuint mvp_location = shader.mvp_uniform_location;

            glEnableVertexAttribArray(0);
            glEnableVertexAttribArray(1);
            glEnableVertexAttribArray(2);
            glEnableVertexAttribArray(3);
            glEnableVertexAttribArray(4);
            glEnableVertexAttribArray(5);
            glEnableVertexAttribArray(6);

            GLuint n_verts = 3 * n_triangles_for_atom_pull_restraints;
            err = glGetError();
            if (err) std::cout << "      error draw_atom_pull_restraints() pre-glDrawElements() "
                               << n_verts << " with GL err " << err << std::endl;
            glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
            err = glGetError();
            if (err) std::cout << "   error in draw_atom_pull_restraints() glDrawElements() n_verts: "
                               << n_verts << " with GL err " << err << std::endl;
         }
      }
   }
}

void
graphics_info_t::draw_molecular_triangles() {
#ifdef USE_MOLECULES_TO_TRIANGLES
   // Martin's triangular molecules
   //
   // centre of the screen
   FCXXCoord pos(graphics_info_t::RotationCentre_x(),
                 graphics_info_t::RotationCentre_y(),
                 graphics_info_t::RotationCentre_z());

   glm::vec3 eye_position = get_world_space_eye_position();
   FCXXCoord eye_pos(eye_position.x, eye_position.y, eye_position.z);

   // std::cout << "eye_pos: " << eye_pos << "\n";
   // coot::Cartesian eye_cart = pos + 20 * diff;
   // FCXXCoord eye_pos(eye_cart.x(), eye_cart.y(), eye_cart.z());
   if (graphics_info_t::mol_tri_scene_setup) {
      if (graphics_info_t::mol_tri_renderer) {
         //Can retrieve reference to the light if so preferred
         FCXXCoord light_pos = pos;
         FCXXCoord neg_light_pos = pos;

         graphics_info_t::mol_tri_scene_setup->getLight(0)->setTranslation(light_pos);
         graphics_info_t::mol_tri_scene_setup->getLight(1)->setTranslation(neg_light_pos);

         for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
            if (graphics_info_t::is_valid_model_molecule(ii)) {
               if (graphics_info_t::molecules[ii].draw_it) {
                  if (graphics_info_t::molecules[ii].molrepinsts.size()) {
                     std::cout << "----------------------- in draw_molecular_triangles() calling Martin code now... \n";
                     // molrepinsts get added to mol_tri_scene_setup when then are made
                     GLenum err = glGetError();
                     if (err) std::cout << "gl error pre-renderer in draw_molecular_triangles() " << err << std::endl;
                     // turns on glLighting.
                     graphics_info_t::mol_tri_scene_setup->renderWithRendererFromViewpoint(graphics_info_t::mol_tri_renderer,
                                                                                           eye_pos);
                     err = glGetError();
                     if (err) std::cout << "gl error in draw_molecular_triangles() " << err << std::endl;
                  }
               }
            }
         }
      }
   }
#endif
}

// static
void
graphics_info_t::draw_particles() {

   if (! particles.empty()) {
      if (mesh_for_particles.have_instances()) {
         glm::mat4 mvp = get_molecule_mvp();
         glm::mat4 view_rotation = get_view_rotation();
         mesh_for_particles.draw_particles(&shader_for_particles, mvp, view_rotation);
      }
   }
}

// static
void
graphics_info_t::draw_happy_face_residue_markers() {

   // make it work (somewhat) like particles, but we are using a screen-facing
   // texture, not a bespoke screen-facing n-triangle polygon.
   // So it's somewhat like the atom labels (a texture in 3D), in perspective
   // happy faces at the back are smaller.
   // But unlike labels, happy faces should be scaled by distance eye to rotation
   // centre, so that they always appear at a constant size (constant n-pixels wide
   // on the screen).

   if (tmesh_for_happy_face_residues_markers.draw_this_mesh) {
      
      if (tmesh_for_happy_face_residues_markers.have_instances()) {

         // the update of the instanced positions is done in the tick function

         graphics_info_t g; // needed for draw_count_max_for_happy_face_residue_markers. Use a better way?
         glm::mat4 mvp = get_molecule_mvp();
         glm::mat4 view_rotation = get_view_rotation();
         texture_for_happy_face_residue_marker.Bind(0);
         unsigned int draw_count = draw_count_for_happy_face_residue_markers;
         unsigned int draw_count_max = g.draw_count_max_for_happy_face_residue_markers;
         tmesh_for_happy_face_residues_markers.draw_instances(&shader_for_happy_face_residue_markers,
                                                              mvp, view_rotation,
                                                              draw_count, draw_count_max);
      }
   }
}


void
graphics_info_t::draw_molecules() {

   // opaque things

   draw_intermediate_atoms();

   draw_intermediate_atoms_rama_balls();

   draw_atom_pull_restraints();

   draw_meshed_generic_display_object_meshes();

   draw_instanced_meshes();

   draw_map_molecules(false); // transparency

   draw_unit_cells();

   draw_environment_graphics_object();

   draw_generic_objects();

   draw_hydrogen_bonds_mesh(); // like boids

   draw_boids();

   draw_particles();

   draw_happy_face_residue_markers();

   // this is the last opaque thing to be drawn because the atom labels are blended.
   // It should be easy to break out the atom label code into its own function. That
   // might be better.
   //
   draw_model_molecules();

   // transparent things...

   draw_map_molecules(true);

}


// This does (draws) symmetry too.
//
// static
void
graphics_info_t::draw_environment_graphics_object() {

#if 0   // old... keep for reference (for a while)
   graphics_info_t g;
   if (is_valid_model_molecule(mol_no_for_environment_distances)) {
      if (g.molecules[mol_no_for_environment_distances].is_displayed_p()) {
      g.environment_graphics_object_internal(environment_object_bonds_box);
      if (g.show_symmetry)
         g.environment_graphics_object_internal(symmetry_environment_object_bonds_box);
      }
   }
#endif

   if (is_valid_model_molecule(mol_no_for_environment_distances)) {
      molecule_class_info_t &m = molecules[mol_no_for_environment_distances];
      if (m.is_displayed_p()) {
         if (environment_show_distances) {
            glm::mat4 mvp = get_molecule_mvp();
            glm::vec3 eye_position = get_world_space_eye_position();
            glm::mat4 view_rotation = get_view_rotation();
            glm::vec4 bg_col(background_colour, 1.0);

            bool do_depth_fog = shader_do_depth_fog_flag;
            mesh_for_environment_distances.mesh.draw(&shader_for_moleculestotriangles,
                                                     mvp, view_rotation,
                                                     lights, eye_position, bg_col,
                                                     do_depth_fog);

            Shader *shader_p = &shader_for_atom_labels;

            GLenum err = glGetError();
            if (err) std::cout << "error draw_environment_graphics_object() before labela err "
                               << err << std::endl;

            if (! labels.empty()) {
               for (unsigned int i=0; i<labels.size(); i++) {
                  const std::string &label  = labels[i].label;
                  const glm::vec3 &position = labels[i].position;
                  const glm::vec4 &colour   = labels[i].colour;
                  // caches these textures in a map std::map<std::string, thing> where
                  // the key is the label.
                  tmesh_for_labels.draw_atom_label(label, position, colour, shader_p,
                                                   mvp, view_rotation, lights, eye_position,
                                                   bg_col, do_depth_fog,
                                                   perspective_projection_flag);
               }
            }

            if (show_symmetry) {

               // Fill me.

            }
         }
      }
   }
}


void
graphics_info_t::draw_unit_cells() {

   glm::mat4 mvp = get_molecule_mvp();
   for (int ii=n_molecules()-1; ii>=0; ii--) {
      molecule_class_info_t &m = molecules[ii];
      m.draw_unit_cell(&shader_for_lines, mvp);
   }

}

void
graphics_info_t::draw_meshed_generic_display_object_meshes() {

   bool draw_meshes = true;
   bool draw_mesh_normals = false;

   glm::vec3 eye_position = get_world_space_eye_position();
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();
   glm::vec4 bg_col(background_colour, 1.0);

   bool do_depth_fog = true;

   //std::cout << "mvp diag "
   // << mvp[0][0] << " " << mvp[1][1] << " " << mvp[2][2] << std::endl;

   glm::mat3 vrm(glm::toMat4(graphics_info_t::glm_quat));
   glm::mat3 vrmt = glm::transpose(vrm);
   glm::mat3 p = vrmt * vrm;

   // Yes, identity matrix
   // std::cout << "p: " << glm::to_string(p) << std::endl;

   if (draw_meshes) { //local, debugging
      bool have_meshes_to_draw = false;
      for (int i=n_molecules()-1; i>=0; i--) {
         if (! molecules[i].meshes.empty()) {
            have_meshes_to_draw = true;
            break;
         }
      }

      if (have_meshes_to_draw) {
         glDisable(GL_BLEND);
         for (int ii=n_molecules()-1; ii>=0; ii--) {
            molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
            for (unsigned int jj=0; jj<m.meshes.size(); jj++) {
               m.meshes[jj].draw(&shader_for_moleculestotriangles, mvp,
                                 view_rotation, lights, eye_position, bg_col, do_depth_fog);
            }
            glUseProgram(0);
         }

         if (draw_mesh_normals) {
            if (draw_normals_flag) {
               for (int ii=n_molecules()-1; ii>=0; ii--) {
                  molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
                  m.mesh_draw_normals(mvp);
               }
            }
         }
      }
   }
}

void
graphics_info_t::draw_instanced_meshes() {

   // presumes opaque-only

   bool have_meshes_to_draw = false;
   for (int i=n_molecules()-1; i>=0; i--) {
      if (! molecules[i].instanced_meshes.empty()) {
         if (molecules[i].draw_it) {
            have_meshes_to_draw = true;
            break;
         }
      }
   }
   if (! instanced_meshes.empty())
      have_meshes_to_draw = true;

   if (have_meshes_to_draw) {
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation = get_view_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool do_depth_fog = shader_do_depth_fog_flag;
      glDisable(GL_BLEND);
      for (int ii=n_molecules()-1; ii>=0; ii--) {
         molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
         if (molecules[ii].draw_it) {
            for (unsigned int jj=0; jj<m.instanced_meshes.size(); jj++) {
               m.instanced_meshes[jj].draw(&shader_for_rama_balls, mvp,
                                           view_rotation, lights, eye_position, bg_col, do_depth_fog);
            }
         }
      }

      // And draw our own

      if (! instanced_meshes.empty()) {
         for (unsigned int jj=0; jj<instanced_meshes.size(); jj++) {
            instanced_meshes[jj].draw(&shader_for_rama_balls, mvp,
                                      view_rotation, lights, eye_position, bg_col, do_depth_fog);
         }
      }
   }
}

void
graphics_info_t::draw_meshes() {

   // presumes opaque-only

   draw_meshed_generic_display_object_meshes();
   draw_instanced_meshes();
}

void
graphics_info_t::draw_cube(GtkGLArea *glarea, unsigned int cube_type) {

   gtk_gl_area_make_current(glarea);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_central_cube() A0 err " << err << std::endl;

// wrap this if you use it again. myglLineWidth()
#ifdef __APPLE__
   // Shut up Mac OS X. This should not give an error
#else
   glLineWidth(2.0);  // GLv4 antialiasing - OpenGL implementations are not required to support this
                      // but they should not create an error. Mac does.
   err = glGetError();
   if (err) std::cout << "error draw_central_cube() A1 glLineWidth() err " << err << std::endl;
#endif


   // To see the possible values of the line width in aliased mode:
   // GLfloat line_width_max_min[2] = {0.0f, 0.0f};
   // glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, lineWidthRange);
   // This may not be possible in GL_LINE_SMOOTH mode.

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation(); // hhmm... naming

   glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);
   err = glGetError(); if (err) std::cout << "   error::draw_central_cube() B err " << err << std::endl;
   glUseProgram(graphics_info_t::shader_for_central_cube.get_program_id());
   err = glGetError(); if (err) std::cout << "   error::draw_central_cube() C err " << err << std::endl;
   glm::mat4 view_orientation = glm::toMat4(graphics_info_t::glm_quat);
   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   if (cube_type == VIEW_CENTRAL_CUBE) {
      mvp = glm::translate(mvp, rc);
      float s = rotation_centre_cube_size;
      glm::vec3 sc(s,s,s);
      mvp = glm::scale(mvp, sc);
   }
   if (cube_type == ORIGIN_CUBE) {
      glm::vec3 sc(0.3f, 0.3f, 0.3f);
      mvp = glm::scale(mvp, sc);
   }

   // we don't diverge here on the cube tye. Maybe change the name of the shader
   // because it does both
   Shader &shader = graphics_info_t::shader_for_central_cube;

   // we do this for all the shaders - Hmm.
   {
      GLuint mvp_location           = shader.mvp_uniform_location;
      GLuint view_rotation_location = shader.view_rotation_uniform_location;

      glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniformMatrix4fv() for mvp " << err << std::endl;
      glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniformMatrix4fv() for view_rotation " << err
                         << std::endl;

      GLuint line_colour_uniform_location = shader.line_colour_uniform_location;
      glm::vec4 lc(0.5, 0.4, 0.4, 1.0);
      if (cube_type == ORIGIN_CUBE)
         lc = glm::vec4(0.6, 0.6, 0.4, 1.0);
      glUniform4fv(line_colour_uniform_location, 1, glm::value_ptr(lc));
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniform4fv() for line colour " << err << std::endl;

      GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
      glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
      glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniform4fv() for background " << err << std::endl;

   }

   glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "error::draw_central_cube() F glDrawElements() err " << err << std::endl;

   glBindVertexArray(0); // unbind
   glUseProgram(0);
}


void
graphics_info_t::draw_central_cube(GtkGLArea *glarea) {
   draw_cube(glarea, VIEW_CENTRAL_CUBE);
}

void
graphics_info_t::draw_origin_cube(GtkGLArea *glarea) {
   draw_cube(glarea, ORIGIN_CUBE);
}

GtkWidget *create_and_pack_gtkglarea(GtkWidget *vbox, bool use_gtk_builder) {

   // the use_gtk_builder flag really means "was invoked from the path that..."

   GtkWidget *w = gtk_gl_area_new();

   auto get_gl_widget_dimension_scale_factor  = [] () {
                                                   int sf = 1;
                                                   char *e = getenv("COOT_OPENGL_WIDGET_SCALE_FACTOR");
                                                   if (e) {
                                                      std::string ee(e);
                                                      sf = std::stoi(ee);
                                                   }
                                                   return sf;
                                                };

   // allow the user to set the major and minor version (for debugging)

   int opengl_major_version = 3;
   int opengl_minor_version = 3;
   char *e1 = getenv("COOT_OPENGL_MAJOR_VERSION");
   char *e2 = getenv("COOT_OPENGL_MINOR_VERSION");
   if (e1) {
      std::string e1s(e1);
      opengl_major_version = std::stoi(e1s);
   }
   if (e2) {
      std::string e2s(e2);
      opengl_minor_version = std::stoi(e2s);
   }
   std::cout << "DEBUG:: seting OpenGL required version to "
             << opengl_major_version << " " << opengl_minor_version << std::endl;

   gtk_gl_area_set_required_version(GTK_GL_AREA(w), opengl_major_version, opengl_minor_version);

   unsigned int dimensions = 700;
   if (! use_gtk_builder) dimensions = 900;
   int gl_widget_dimension_scale_factor = get_gl_widget_dimension_scale_factor();
   gtk_widget_set_size_request(w,
                               gl_widget_dimension_scale_factor * dimensions,
                               gl_widget_dimension_scale_factor * dimensions);
   gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 0);
   return w;
}

void
graphics_info_t::setup_lights() {

   lights_info_t light;
   light.position = glm::vec4(-2.0f, 2.0f, 5.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(0.5, 0, 1.0));
   graphics_info_t::lights[0] = light;

   light.position = glm::vec4(3.0f, -2.0f, 4.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(-1.0, 0.5, 1.0));
   graphics_info_t::lights[1] = light;
}

// this is called from realize()
void
graphics_info_t::setup_hud_geometry_bars() {

   if (! glareas[0]) return;

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   float aspect_ratio = static_cast<float>(w)/static_cast<float>(h);

   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   shader_for_hud_geometry_bars.Use();

   mesh_for_hud_geometry.setup_camera_facing_quad_for_bar();
   mesh_for_hud_geometry.setup_instancing_buffer(500);

   // If not found in this directory, then try default directory.
   texture_for_hud_geometry_labels.set_default_directory(coot::package_data_dir());
   texture_for_hud_geometry_labels.init("hud-label-nbc-rama.png");

   // above?
   texture_for_hud_tooltip_background.set_default_directory(coot::package_data_dir());
   // texture_for_hud_tooltip_background.init("hud-label-nbc-rama.png");
   texture_for_hud_tooltip_background.init("hud-tooltip.png"); // 94x47
   float sc_x = 0.001 * static_cast<float>(103) / aspect_ratio;
   float sc_y = 0.001 * static_cast<float>(50);

   // Do I need to Use() the shader_for_hud_geometry_labels here?
   shader_for_hud_geometry_labels.Use();
   mesh_for_hud_geometry_labels.setup_quad();
   glm::vec2 position(-0.98, 0.903);
   glm::vec2 scales(0.04/aspect_ratio, 0.06);
   mesh_for_hud_geometry_labels.set_position_and_scales(position, scales); // ""NBC, Pull"" texture

   mesh_for_hud_tooltip_background.setup_quad(); // does setup_buffers()
   mesh_for_hud_tooltip_background.set_scales(glm::vec2(sc_x, sc_y));
   // mesh_for_hud_tooltip_background.setup_texture_coords_for_nbcs_and_rama();  // surely not needed?

   tmesh_for_hud_geometry_tooltip_label.setup_quad();
   glm::vec2 label_scale(0.00015/aspect_ratio, 0.00015);
   tmesh_for_hud_geometry_tooltip_label.set_scales(label_scale);

}

float
graphics_info_t::hud_geometry_distortion_to_bar_size_rama(float distortion) {
   float d1 = distortion + 200.0;
   float d2 = d1 * 0.0003;
   if (d2 < 0.0) d2 = 0.0;
   float d3 = 100.0 * d2 * d2;
   return d3;
}

float
graphics_info_t::hud_geometry_distortion_to_bar_size_nbc(float distortion) {
   return distortion * 0.002;
}


float
graphics_info_t::hud_geometry_distortion_to_bar_size_atom_pull(float distortion) {
   return distortion * 0.0008;
}


float
graphics_info_t::hud_geometry_distortion_to_rotation_amount_rama(float distortion) {
   distortion += 200.0;
   // float rotation_amount = 1.0 - 0.0022 * distortion;
   float rotation_amount = 1.0 - 0.0028 * distortion;
   if (rotation_amount < 0.68) rotation_amount = 0.68; // red cap
   if (rotation_amount > 1.0) rotation_amount = 1.0;
   return rotation_amount;
}


void
graphics_info_t::draw_hud_geometry_bars() {

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glEnable(GL_DEPTH_TEST); // needed?
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   // first draw the text (labels) texture

   coot::refinement_results_t &rr = saved_dragged_refinement_results;

   if (! rr.refinement_results_contain_overall_rama_plot_score) {
      // std::cout << "chopping the texture " << std::endl;
      mesh_for_hud_geometry_labels.setup_texture_coords_for_nbcs_only();
   } else {
      mesh_for_hud_geometry_labels.setup_texture_coords_for_nbcs_and_rama();
   }

   // Restore me when fixed                                                       
   // texture_for_hud_geometry_labels.Bind(0);
   // mesh_for_hud_geometry_labels.draw(&shader_for_hud_geometry_labels);

   // now draw the bars

   auto probability_to_rotation_amount = [] (float probability) {
                                            // high probability should have low rotation
                                            float q_1 = 0.01 * (100.0 - probability);
                                            // if (q < 0) q = 0;
                                            // if (q > 1) q = 1;
                                            float q_2 = 0.68 * q_1;
                                            return q_2;
                                         };

   auto distortion_to_rotation_amount_nbc = [] (float distortion) {
                                               // we want to rotate to red (which is negative direction) but
                                               // rotate() doesn't work with negative rotations, so make it
                                               // 1.0 - amount (1.0 being a full rotation).
                                               float rotation_amount = 1.0 - 0.012 * distortion;
                                               if (rotation_amount < 0.68) rotation_amount = 0.68;
                                               return rotation_amount;
                                            };

   auto add_bars = [] (const std::vector<std::pair<coot::atom_spec_t, float> > &baddies,
                       unsigned int bar_index,
                       std::vector<HUD_bar_attribs_t> *new_bars_p,
                       float (*distortion_to_rotation_amount)(float),
                       float (*distortion_to_bar_size)(float)) {

                         glm::vec2 to_top_left(-0.90, 0.943 - 0.05 * static_cast<float>(bar_index));
                         float sum_l = 0;
                         int n = baddies.size();
                         glm::vec4 col_white(0.8,0.8, 0.8, 0.7);
                         for (int i=(n-1); i>=0; i--) {
                            coot::colour_t cc(0.1, 0.9, 0.2);
                            float d = baddies[i].second;
                            float rotation_amount = distortion_to_rotation_amount(d);
                            cc.rotate(rotation_amount);
                            glm::vec4 col = cc.to_glm();
                            col.w = 0.7;
                            glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                            float bar_length = distortion_to_bar_size(d);
                            bool this_atom_is_in_a_moving_atoms_residue = baddies[i].first.int_user_data;

                            // the vector of HUD_bar_attribs_t is fed directly to a opengl buffer.
                            // So I need "expand" to 2 bars right here - one with a "thin bar" attribute
                            //
                            if (! this_atom_is_in_a_moving_atoms_residue) {
                               float bar_height = 0.03; // universal
                               float bar_slither_y_scale = 0.3; // looks nice
                               float y_offset_main = bar_height * bar_slither_y_scale;
                               glm::vec2 position_offset_for_main = position_offset + glm::vec2(0, y_offset_main);
                               HUD_bar_attribs_t bar_main(col, position_offset_for_main, bar_length);
                               bar_main.scale_y = 1.0 - bar_slither_y_scale;
                               new_bars_p->push_back(bar_main);
                               // slither bar
                               glm::vec2 position_offset_for_slither = position_offset + glm::vec2(0,0);
                               HUD_bar_attribs_t bar_slither(col_white, position_offset_for_slither, bar_length);
                               bar_slither.scale_y = bar_slither_y_scale;
                               new_bars_p->push_back(bar_slither);
                            } else {
                               HUD_bar_attribs_t bar(col, position_offset, bar_length);
                               new_bars_p->push_back(bar);
                            }
                            sum_l += bar_length + 0.005; // with a gap between bars
                          }
                       };

   auto rota_sorter = [] (const rotamer_markup_container_t &rmc_1,
                          const rotamer_markup_container_t &rmc_2) {
                         if (rmc_2.rpi.probability < rmc_1.rpi.probability)
                            return true;
                         else
                            return false;
                      };

   auto add_rotamer_bars = [rota_sorter] (std::vector<HUD_bar_attribs_t> *new_bars_p,
                                          unsigned int bar_index,
                                          rotamer_markup_container_t *rotamer_markups,
                                          int n_rotamer_markups) {

                              // this code has to be the same as the check_if_hud_bar_clicked code

                              glm::vec4 col_white(0.8,0.8, 0.8, 0.7);
                              glm::vec2 to_top_left(-0.90, 0.943 - 0.05 * static_cast<float>(bar_index));
                              std::vector<rotamer_markup_container_t> v;
                              // filter out the goodies
                              for (int i=0; i<n_rotamer_markups; i++)
                                 if (rotamer_markups[i].rpi.probability < 40) // 40 %
                                    if (rotamer_markups[i].rpi.probability >= 0)
                                       v.push_back(rotamer_markups[i]);
                              // sort the baddies
                              std::sort(v.begin(), v.end(), rota_sorter);
                              unsigned int n_rota_max = 20;
                              if (v.size() > n_rota_max) {
                                 unsigned int n_for_deletion = v.size() - n_rota_max;
                                 std::vector<rotamer_markup_container_t>::iterator v_begin = v.begin();
                                 std::vector<rotamer_markup_container_t>::iterator v_last;
                                 v_last = v_begin + n_for_deletion; // (line length)
                                 v.erase(v_begin, v_last);
                              }

                              float sum_l = 0;
                              for (unsigned int i=0; i<v.size(); i++) {

                                 bool this_atom_is_in_a_moving_atoms_residue = v[i].spec.int_user_data;

                                 float pr = v[i].rpi.probability;
                                 float q = 0.01 * (48.0f - v[i].rpi.probability);
                                 if (q > 1.0) q = 1.0;
                                 if (q < 0.0) q = 0.0;
                                 float bar_length = std::pow(q, 6.0) * 4.0;

                                 if (false)
                                    std::cout << "bar i " << i << " " << v[i].spec << " " << v[i].col
                                              << " pr " << pr << " length " << bar_length << std::endl;

                                 const coot::colour_holder &ch = v[i].col;
                                 glm::vec4 col(ch.red, ch.green, ch.blue, 0.7);

                                 if (! this_atom_is_in_a_moving_atoms_residue) {
                                    float bar_height = 0.03; // universal
                                    float bar_slither_y_scale = 0.3; // looks nice
                                    float y_offset_main = bar_height * bar_slither_y_scale;
                                    glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                                    glm::vec2 position_offset_for_main = position_offset + glm::vec2(0, y_offset_main);
                                    HUD_bar_attribs_t bar_main(col, position_offset_for_main, bar_length);
                                    bar_main.scale_y = 1.0 - bar_slither_y_scale;
                                    new_bars_p->push_back(bar_main);
                                    // slither bar
                                    glm::vec2 position_offset_for_slither = position_offset + glm::vec2(0,0);
                                    HUD_bar_attribs_t bar_slither(col_white, position_offset_for_slither, bar_length);
                                    bar_slither.scale_y = bar_slither_y_scale;
                                    new_bars_p->push_back(bar_slither);


                                 } else {
                                    glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                                    HUD_bar_attribs_t bar(col, position_offset, bar_length);
                                    new_bars_p->push_back(bar);
                                 }
                                 sum_l += bar_length + 0.005; // with a gap between bars
                              }
                           };

   std::vector<HUD_bar_attribs_t> new_bars;

   // add to new_bars
   add_bars(rr.sorted_atom_pulls, 0, &new_bars,
            distortion_to_rotation_amount_nbc, hud_geometry_distortion_to_bar_size_atom_pull);

   if (rr.refinement_results_contain_overall_nbc_score)
      add_bars(rr.sorted_nbc_baddies, 1, &new_bars,
               distortion_to_rotation_amount_nbc, hud_geometry_distortion_to_bar_size_nbc);

   if (rr.refinement_results_contain_overall_rama_plot_score)
      add_bars(rr.sorted_rama_baddies, 2, &new_bars,
               hud_geometry_distortion_to_rotation_amount_rama, hud_geometry_distortion_to_bar_size_rama);

   // add rotas to new_bars
   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         int nrms = moving_atoms_molecule.bonds_box.n_rotamer_markups;
         if (nrms > 0) {
            add_rotamer_bars(&new_bars, 3, moving_atoms_molecule.bonds_box.rotamer_markups, nrms);
         }
      }
   }

   if (! new_bars.empty()) {
      // std::cout << "new bar size " << new_bars.size() << std::endl;
      mesh_for_hud_geometry.update_instancing_buffer_data(new_bars);
      Shader &shader = shader_for_hud_geometry_bars;
      mesh_for_hud_geometry.draw(&shader);
   }
   glDisable(GL_BLEND);

}

std::pair<bool, mmdb::Atom *>
graphics_info_t::check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(double mouse_x, double mouse_y, bool act_on_hit) {

   std::pair<bool, mmdb::Atom *> status_pair(false, 0);
   if (! moving_atoms_asc) return std::pair<bool, mmdb::Atom *>(false, 0);
   if (! moving_atoms_asc->mol) return std::pair<bool, mmdb::Atom *>(false, 0);

   coot::refinement_results_t &rr = saved_dragged_refinement_results;

   // this values in this loop must match those in the loop above
   // (draw_hud_geometry_bars())

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   double frac_x = mouse_x/static_cast<double>(w);
   double frac_y = 1.0 - mouse_y/static_cast<double>(h);
   glm::vec2 mouse_in_opengl_coords(2.0 * frac_x - 1.0, 2.0 * frac_y - 1.0);

   // these functions are copies of those in the above function. If you edit them again, make them
   // member functions.

   auto rota_sorter = [] (const rotamer_markup_container_t &rmc_1,
                          const rotamer_markup_container_t &rmc_2) {
                         if (rmc_2.rpi.probability < rmc_1.rpi.probability)
                            return true;
                         else
                            return false;
                      };

   // the act_on_hit flag is check to see if the move should be made, or that we merely return
   // a success status (we want to act when clicked, but return a status when moused-over)
   //
   auto check_blocks = [mouse_in_opengl_coords] (const std::vector<std::pair<coot::atom_spec_t, float> > &baddies,
                                                 unsigned int bar_index,
                                                 float (*distortion_to_bar_size)(float),
                                                 bool act_on_hit) {

                          bool status = false;
                          mmdb::Atom *at_out = 0;
                          glm::vec2 to_top_left(-0.90, 0.943 - 0.05 * static_cast<float>(bar_index));
                          float sum_l = 0;
                          int n = baddies.size();
                          for (int i=(n-1); i>=0; i--) {
                             float d = baddies[i].second;
                             glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                             float bar_length = distortion_to_bar_size(d);
                             sum_l += bar_length + 0.005; // with a gap between bars

                             if (false) {
                                glm::vec2 position_offset_far_point = position_offset;
                                position_offset_far_point.x += bar_length;
                                std::cout << "checking "
                                          << glm::to_string(position_offset) << " "
                                          << glm::to_string(position_offset_far_point) << " "
                                          << " " << bar_length
                                          << " vs mouse " << glm::to_string(mouse_in_opengl_coords)
                                          << std::endl;
                             }

                             if (mouse_in_opengl_coords.x >= position_offset.x) {
                                if (mouse_in_opengl_coords.x <= (position_offset.x + bar_length)) {

                                   // std::cout << ":::::::::: x hit bar_index " << bar_index
                                   // << " i " << i << " " << baddies[i].first << std::endl;

                                   float tiny_y_offset = -0.01; // not sure why I need this
                                   if (mouse_in_opengl_coords.y >= (to_top_left.y + tiny_y_offset)) {
                                      // 0.03 is the bar height in setup_camera_facing_quad()
                                      float bar_height = 0.03;
                                      if (mouse_in_opengl_coords.y <= (to_top_left.y+tiny_y_offset+bar_height)) {
                                         coot::atom_spec_t spec(baddies[i].first);
                                         if (moving_atoms_asc->mol) {
                                            mmdb::Atom *at = spec.get_atom(moving_atoms_asc->mol);
                                            if (at) {
                                               at_out = at;
                                               status = true;
                                               if (act_on_hit) {
                                                  clipper::Coord_orth pt = coot::co(at);
                                                  std::cout << "INFO: geom bar atom: " << coot::atom_spec_t(at)
                                                            << std::endl;
                                                  set_rotation_centre(pt);
                                               }
                                            }
                                         } else {
                                            std::cout << "ERROR:: no moving atoms mol" << std::endl;
                                         }
                                      }
                                   }
                                }
                             }
                          }
                          return std::pair<bool, mmdb::Atom *>(status, at_out);
                       };

   auto check_rota_blocks = [mouse_in_opengl_coords,
                             rota_sorter] (unsigned int bar_index,
                                           rotamer_markup_container_t *rotamer_markups,
                                           int n_rotamer_markups,
                                           bool act_on_hit) {

                               bool status = false;
                               mmdb::Atom *at_out = 0;
                               coot::residue_spec_t spec_for_at_out;
                               glm::vec2 to_top_left(-0.90, 0.943 - 0.05 * static_cast<float>(bar_index));
                               std::vector<rotamer_markup_container_t> v;
                               // filter out the goodies
                               for (int i=0; i<n_rotamer_markups; i++)
                                  if (rotamer_markups[i].rpi.probability < 40) // 40 %
                                    if (rotamer_markups[i].rpi.probability >= 0)
                                       v.push_back(rotamer_markups[i]);
                               // sort the baddies
                               std::sort(v.begin(), v.end(), rota_sorter);
                               unsigned int n_rota_max = 20;
                               if (v.size() > n_rota_max) {
                                  unsigned int n_for_deletion = v.size() - n_rota_max;
                                  std::vector<rotamer_markup_container_t>::iterator v_begin = v.begin();
                                  std::vector<rotamer_markup_container_t>::iterator v_last  = v_begin + n_for_deletion;
                                  v.erase(v_begin, v_last);
                               }

                               float sum_l = 0;
                               for (unsigned int i=0; i<v.size(); i++) {
                                  // float pr = v[i].rpi.probability;
                                  float q = 0.01 * (48.0f - v[i].rpi.probability);
                                  if (q > 1.0) q = 1.0;
                                  if (q < 0.0) q = 0.0;
                                  float bar_length = std::pow(q, 6.0) * 4.0;
                                  glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);

                                  if (mouse_in_opengl_coords.x >= position_offset.x) {
                                     if (mouse_in_opengl_coords.x <= (position_offset.x + bar_length)) {
                                        // std::cout << ":::::::::: x hit bar_index " << bar_index
                                        //           << " i " << i << " " << baddies[i].first << std::endl;
                                        float tiny_y_offset = -0.01; // not sure why I need this
                                        if (mouse_in_opengl_coords.y >= (to_top_left.y + tiny_y_offset)) {
                                           // 0.03 is the bar height in setup_camera_facing_quad()
                                           float bar_height = 0.03;
                                           if (mouse_in_opengl_coords.y <= (to_top_left.y+tiny_y_offset+bar_height)) {

                                              if (false)
                                                 std::cout << "rama bar hit! " << i << " "
                                                           << v[i].spec << " "
                                                           << v[i].col << " "
                                                           << "probability " << v[i].rpi.probability << std::endl;

                                              if (moving_atoms_asc->mol) {
                                                 status = true;
                                                 spec_for_at_out = v[i].spec;
                                                 if (act_on_hit) {
                                                    clipper::Coord_orth pos = v[i].pos;
                                                    graphics_info_t::set_rotation_centre(pos);
                                                 }
                                              }
                                           }
                                        }
                                     }
                                  }
                                  sum_l += bar_length + 0.005; // with a gap between bars
                               }
                               // I need this function to return an atom (so that it's like the other geometry bars)
                               // - because that atom spec gets turned into an tooltip label.
                               if (status) {
                                  // I first need to find the first residue
                                  mmdb::Residue *residue_p = spec_for_at_out.get_residue(moving_atoms_asc->mol);
                                  if (residue_p) {
                                     int n_atoms = residue_p->GetNumberOfAtoms();
                                     if (n_atoms > 0) at_out = residue_p->GetAtom(0);
                                     if (n_atoms > 1) at_out = residue_p->GetAtom(1); // CA, usually
                                  }
                               }
                               return std::pair<bool, mmdb::Atom *>(status, at_out);
                            };

   status_pair = check_blocks(rr.sorted_atom_pulls, 0, hud_geometry_distortion_to_bar_size_atom_pull, act_on_hit);

   if (!status_pair.first)
      if (rr.refinement_results_contain_overall_nbc_score)
         status_pair = check_blocks(rr.sorted_nbc_baddies, 1, hud_geometry_distortion_to_bar_size_nbc, act_on_hit);

   if (!status_pair.first)
      if (rr.refinement_results_contain_overall_rama_plot_score)
         status_pair = check_blocks(rr.sorted_rama_baddies, 2, hud_geometry_distortion_to_bar_size_rama, act_on_hit);

   if (!status_pair.first) {
      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            int nrms = moving_atoms_molecule.bonds_box.n_rotamer_markups;
            if (nrms > 0) {
               status_pair = check_rota_blocks(3, moving_atoms_molecule.bonds_box.rotamer_markups, nrms, act_on_hit);
            }
         }
      }
   }


   if (act_on_hit) {
      if (status_pair.first) {
         mmdb::Atom *at = status_pair.second;
         if (at) {
            mmdb::Residue *residue_p = at->residue;
            moving_atoms_visited_residues.insert(residue_p);
            active_atom_for_hud_geometry_bar = at;
         }
      }
   }

   return status_pair;
}

std::pair<bool, mmdb::Atom *>
graphics_info_t::check_if_moused_over_hud_bar(double mouse_x, double mouse_y) {

   // copied from check_if_hud_bar_clicked().
   // Now that we have act_on_hit, we can extract most of this code to a common function
   // called check_if_hud_bar_mouse_over_or_act_on_hurd_bar_click()

   bool act_on_hit = false;
   return check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(mouse_x, mouse_y, act_on_hit);
}

bool
graphics_info_t::check_if_hud_bar_clicked(double mouse_x, double mouse_y) {

   if (! moving_atoms_asc) return false;
   if (! moving_atoms_asc->mol) return false;
   bool act_on_hit = true;
   std::pair<bool, mmdb::Atom *> r = check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(mouse_x, mouse_y, act_on_hit);
   return  r.first;
}

void
graphics_info_t::draw_hud_geometry_tooltip() {

   // this flag is set when the user mouses over a HUD bar
   // and removed when they move from a hud bar.

   glEnable(GL_DEPTH_TEST);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   if (draw_hud_tooltip_flag) {

      texture_for_hud_tooltip_background.Bind(0);
      mesh_for_hud_tooltip_background.draw(&shader_for_hud_geometry_labels);

      // now the text that goes into (on top of) the background

      std::string label = "W 356 CA"; // checking the HUD bars (on mouse-over) should
                                     // return this - and the position, which means
                                     // that they need to return
                                     // something other than a bool. Store it in
                                     // graphics_info_t somewhere.
                                     // HUD_geometry_tooltip_text_position
                                     // HUD_geometry_tooltip_text_label
      label = label_for_hud_geometry_tooltip;

      // do this elsewhere (on_glarea_motion_notify())
      // glm::vec2 label_position(-0.64, 0.72);
      // tmesh_for_hud_geometry_tooltip_label.set_position(label_position);

      // this is now done in setup_hud_geometry_bars() which is called by resize()
      // glm::vec2 label_scale(0.00015, 0.00015); // fixed.
      // tmesh_for_hud_geometry_tooltip_label.set_scales(label_scale);


      bool use_label_highlight = true;
      mmdb::Residue *residue_p = 0;
      if (active_atom_for_hud_geometry_bar)
         residue_p = active_atom_for_hud_geometry_bar->residue;
      if (residue_p)
         if (moving_atoms_visited_residues.find(residue_p) != moving_atoms_visited_residues.end())
            use_label_highlight = false;

      tmesh_for_hud_geometry_tooltip_label.draw_label(label, use_label_highlight,
                                                      &shader_for_hud_geometry_tooltip_text,
                                                      ft_characters);
   }
}


gboolean
graphics_info_t::render(bool to_screendump_framebuffer, const std::string &output_file_name) {

   // auto tp_0 = std::chrono::high_resolution_clock::now();

   auto render_scene = [] (GtkGLArea *gl_area) {
                          
                          //  ------------------- render scene ----------------------------
                          const glm::vec3 &bg = graphics_info_t::background_colour;
                          glClearColor (bg[0], bg[1], bg[2], 1.0); // what difference does this make?
                          GLenum err = glGetError(); if (err) std::cout << "render_scene lambda B err " << err << std::endl;
                          glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                          err = glGetError(); if (err) std::cout << "render_scene lambda C err " << err << std::endl;

                          draw_central_cube(gl_area);
                          draw_origin_cube(gl_area);
                          err = glGetError(); if (err) std::cout << "render scene lambda post cubes err " << err << std::endl;

                          draw_molecules(); // includes particles, happy-faces and boids (should they be there (maybe not))

                          draw_hud_geometry_bars();

                          draw_hud_geometry_tooltip(); // background and text

                          draw_identification_pulse();

                          draw_delete_item_pulse();

                          draw_invalid_residue_pulse();

                          draw_ligand_view();

                          glBindVertexArray(0); // here is not the place to call this.
                       };

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   if (use_framebuffers) { // static class variable

      GLenum err = glGetError();
      if (err) std::cout << "render() start " << err << std::endl;

      // is this needed? - does the context ever change?
      gtk_gl_area_make_current(gl_area);
      err = glGetError();
      if (err) std::cout << "render() post gtk_gl_area_make_current() err " << err << std::endl;

      glViewport(0, 0, framebuffer_scale * w, framebuffer_scale * h);
      err = glGetError();
      if (err) std::cout << "render() post glViewport() err " << err << std::endl;
      screen_framebuffer.bind();
      err = glGetError();
      if (err) std::cout << "render() post screen_framebuffer bind() err " << err << std::endl;

      glEnable(GL_DEPTH_TEST);

      render_scene(gl_area);

      if (to_screendump_framebuffer) {

         glDisable(GL_DEPTH_TEST);
         unsigned int sf = framebuffer_scale;
         glViewport(0, 0, sf * w, sf * h);
         framebuffer screendump_framebuffer;
         unsigned int index_offset = 0;
         screendump_framebuffer.init(sf * w, sf * h, index_offset, "screendump");
         screendump_framebuffer.bind();
         render_scene_to_base_framebuffer();
         gtk_gl_area_attach_buffers(gl_area);
         screendump_tga_internal(output_file_name, w, h, sf, screendump_framebuffer.get_fbo());

      } else {
         glViewport(0, 0, w, h);
         // use this, rather than glBindFramebuffer(GL_FRAMEBUFFER, 0); ... just Gtk things.
         gtk_gl_area_attach_buffers(gl_area);
         render_scene_to_base_framebuffer(); // render current framebuffer to base framebuffer
      }
   } else {
      //  simple/direct - for debugging framebuffers
      gtk_gl_area_attach_buffers(gl_area);
      glEnable(GL_DEPTH_TEST); // needed?
      render_scene(gl_area);
   }

   glEnable(GL_DEPTH_TEST); // what's this for?

   // auto tp_1 = std::chrono::high_resolution_clock::now();
   // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: render() " << d10 << " microseconds" << std::endl;

   return FALSE;
}

void
graphics_info_t::render_scene_to_base_framebuffer() {

   glEnable(GL_DEPTH_TEST);
   shader_for_screen.Use();
   glBindVertexArray(screen_quad_vertex_array_id);

   const glm::vec3 &bg = background_colour;
   glClearColor(bg[0], bg[1], bg[2], 1.0); // this can be seen
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glActiveTexture(GL_TEXTURE0 + 0);
   glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_depth());
   shader_for_screen.set_int_for_uniform("screenTexture", 0);
   shader_for_screen.set_int_for_uniform("screenDepth", 1);
   GLenum err = glGetError(); if (err) std::cout << "render() D err " << err << std::endl;
   shader_for_screen.set_bool_for_uniform("do_ambient_occlusion", shader_do_ambient_occlusion_flag);
   shader_for_screen.set_bool_for_uniform("do_outline", shader_do_outline_flag);

   glDrawArrays(GL_TRIANGLES, 0, 6);
   err = glGetError(); if (err) std::cout << "render() E err " << err << std::endl;

}


void
graphics_info_t::reset_frame_buffers(int width, int height) {

   if (use_framebuffers) {
      unsigned int sf = framebuffer_scale;
      unsigned int index_offset = 0;
      // std::cout << "debug:: reset_frame_buffers() with sf " << sf << " "
      // << width << " x " << height << std::endl;
      screen_framebuffer.init(sf * width, sf * height, index_offset, "screen");
      GLenum err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      // index_offset = 0;
      // g.blur_framebuffer.init(width, height, index_offset, "blur");
   }

}

void
graphics_info_t::try_label_unlabel_active_atom() {

   std::pair<int, mmdb::Atom *> aa = get_active_atom();
   int im = aa.first;
   if (im >= 0) {
      mmdb::Atom *at = aa.second;
      if (at) {
         int atom_index;
         // this is a bit convoluted :-)
         int ierr = at->GetUDData(molecules[im].atom_sel.UDDAtomIndexHandle, atom_index);
	 if (ierr == mmdb::UDDATA_Ok) {
            molecules[im].add_to_labelled_atom_list(atom_index);
            graphics_draw();
         } else {
            std::cout << "WARNING:: Bad UDData for atom_index for atom " << std::endl;
         }
      }
   }
}


// static
glm::vec3
graphics_info_t::get_screen_y_uv() {

   glm::vec3 minus_y = graphics_info_t::unproject_to_world_coordinates(glm::vec3(0.0f, -1.0f, 0.0f));
   glm::vec3  plus_y = graphics_info_t::unproject_to_world_coordinates(glm::vec3(0.0f,  1.0f, 0.0f));
   glm::vec3 delta = plus_y - minus_y;
   glm::vec3 d_uv = glm::normalize(delta);
   return d_uv;
}


void
graphics_info_t::translate_in_screen_z(float step_size) {

   // this looks a bit weird without perspective view

   glm::vec3 ep = get_world_space_eye_position();
   glm::vec3 rc = get_rotation_centre();
   glm::vec3 delta = rc - ep;
   glm::vec3 delta_uv = normalize(delta);

   // more zoomed in has smaller zoom than zoomed out. Zoomed out is ~100. Zoomed in is ~25
   glm::vec3 step = 0.005 * step_size * zoom * delta_uv;

   if (false) // debug
      std::cout << "ep " << glm::to_string(ep) << " rc " << glm::to_string(rc)
                << " zoom " << zoom << " step " << glm::to_string(step) << std::endl;

   add_to_rotation_centre(step);

}

// static
std::vector<glm::vec3>
graphics_info_t::get_particle_centre_positions() {

   auto mmdb_to_glm = [] (mmdb::Atom *at) { return glm::vec3(at->x, at->y, at->z); };

   get_moving_atoms_lock(__FUNCTION__);

   std::vector<glm::vec3> v;
   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
            mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
            if (! at->isTer()) {
               std::string atom_name(at->GetAtomName());
               if (atom_name == " CA " || atom_name == " N1 " || atom_name == " N9 ") {
                  glm::vec3 p = mmdb_to_glm(at);
                  v.push_back(p);
               }
            }
         }
      }
   }
   release_moving_atoms_lock(__FUNCTION__);

   if (v.empty()) {
      glm::vec3 rc = get_rotation_centre();
      v.push_back(rc);
   }

   return v;
}

void
graphics_info_t::setup_draw_for_particles() {

   if (false) // from the days when particle drawing was a problem!
      std::cout << "setup_draw_for_particles(): -- start -- n_particles " << particles.size()
                <<  std::endl;

   if (particles.empty()) {
      std::cout << "setup_draw_for_particles(): let's make new particles " << std::endl;

      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
      GLenum err = glGetError();
      if (err) std::cout << "Error:: setup_draw_for_particles() Post attach buffers err is "
                         << err << std::endl;

      shader_for_particles.Use();

      err = glGetError();
      if (err) std::cout << "Error::- setup_draw_for_particles() Post Use() err is "
                         << err << std::endl;

      std::vector<glm::vec3> positions = get_particle_centre_positions();
      particles.make_particles(n_particles, positions);
      std::cout << "setup_draw_for_particles(): done making " << n_particles << " particles "
                << " for " << positions.size() << " positions" << std::endl;

      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
      mesh_for_particles.setup_vertex_and_instancing_buffers_for_particles(particles.size());
      mesh_for_particles.update_instancing_buffer_data_for_particles(particles);
      glUseProgram(0);
   }
   // passing user_data and Notify function at the end
   if (! do_tick_particles) {
      do_tick_particles = true;
      int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
      graphics_info_t::idle_function_spin_rock_token = new_tick_id;
   }

   // std::cout << "setup_draw_for_particles(): -- done -- " << std::endl;
}

void
graphics_info_t::setup_draw_for_happy_face_residue_markers_init() {

   // run this once - call from realize()

   const unsigned int max_happy_faces = 200; // surely enough?
   
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   GLenum err = glGetError();
   if (err) std::cout << "Error::- setup_draw_for_happy_face_residue_markers_init() "
                      << "Post attach buffers err is " << err << std::endl;

   // If not found in this directory, then try default directory.
   texture_for_happy_face_residue_marker.set_default_directory(coot::package_data_dir());
   texture_for_happy_face_residue_marker.init("happy-face-marker.png");

   shader_for_happy_face_residue_markers.Use();
   tmesh_for_happy_face_residues_markers.setup_camera_facing_quad(&shader_for_happy_face_residue_markers, 1.0, 1.0);
   tmesh_for_happy_face_residues_markers.setup_instancing_buffers(max_happy_faces);
   tmesh_for_happy_face_residues_markers.draw_this_mesh = false;


}

void
graphics_info_t::setup_draw_for_happy_face_residue_markers() {

   // run this at the start of a "show animated happy faces"

   std::vector<glm::vec3> positions = get_happy_face_residue_marker_positions();
   happy_face_residue_marker_starting_positions = positions;
   
   glm::vec3 up_uv = get_screen_y_uv();
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?

   if (false)
      std::cout << "setup_draw_for_happy_face_residue_markers() calling update_instancing_buffer_data()"
                << " with draw_count_for_happy_face_residue_markers "
                << draw_count_for_happy_face_residue_markers << std::endl;
   unsigned int n_max = draw_count_max_for_happy_face_residue_markers;
   tmesh_for_happy_face_residues_markers.update_instancing_buffer_data_for_happy_faces(positions, 0, n_max, up_uv);
   tmesh_for_happy_face_residues_markers.draw_this_mesh = true;
   do_tick_happy_face_residue_markers = true;
   draw_count_for_happy_face_residue_markers = 0;
   gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);

}

#include "coot-utils/fib-sphere.hh"

std::vector<glm::vec3>
graphics_info_t::get_happy_face_residue_marker_positions() {

   const unsigned int max_happy_faces = 200; // surely enough? - If changed, change above
   std::vector<glm::vec3> v;
   bool make_fake_points = false;

   auto clipper_to_glm = [] (const clipper::Coord_orth &co) {
                            return glm::vec3(co.x(), co.y(), co.z());
                         };

   if (make_fake_points) {
      glm::vec3 sc = get_rotation_centre();
      std::vector<clipper::Coord_orth> cv = coot::fibonacci_sphere(80);
      for (auto p : cv)
         v.push_back(sc + 5.0 * clipper_to_glm(p));
   } else {

      // This is just a bit of fun... actually, I will need to ask something like
      // last_restraints->get_improved_residues();
      // last_restraints->get_damaged_residues();

      // How about the set-fixed-during-refinement udd?
      // see set_fixed_during_refinement_udd()

      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            int uddHnd = moving_atoms_asc->mol->GetUDDHandle(mmdb::UDR_ATOM , "FixedDuringRefinement");
            std::vector<mmdb::Residue *> residues;
            int imod = 1;
            mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        // I need to not add this if it's a fixed residue. How can I know that?
                        bool is_fixed = false;
                        mmdb::Atom **residue_atoms = 0;
                        int n_residue_atoms = 0;
                        residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                        for(int iat=0; iat<n_residue_atoms; iat++) {
                           mmdb::Atom *at = residue_atoms[iat];
                           int is_fixed_status = 0;
                           int ierr = at->GetUDData(uddHnd, is_fixed_status);
                           if (ierr == mmdb::Error_Ok) {
                              if (is_fixed_status == 1) {
                                 is_fixed = true;
                                 break;
                              }
                           }
                        }
                        if (! is_fixed)
                           residues.push_back(residue_p);
                     }
                  }
               }
            }

            for (auto r : residues) {
               std::pair<bool, clipper::Coord_orth> rc = coot::util::get_residue_centre(r);
               if (rc.first) {
                  glm::vec3 p = clipper_to_glm(rc.second);
                  v.push_back(p);
               }
            }
         }
      }
   }

   
   if (v.size() > max_happy_faces)
      std::cout << "error:: ------------------ too many happy faces" << std::endl;

   return v;
}


//static
gboolean
graphics_info_t::wait_for_hooray_refinement_tick_func(GtkWidget *widget,
                                                      GdkFrameClock *frame_clock,
                                                      gpointer data) {
   gboolean continue_status = 1;

   if (setup_draw_for_particles_semaphore) {
      if (! particles_have_been_shown_already_for_this_round_flag) {
         graphics_info_t g;
         g.setup_draw_for_particles();
         setup_draw_for_particles_semaphore = false; // it's done it's job
         particles_have_been_shown_already_for_this_round_flag = true; // only once per round
         continue_status = 0; // job done.
      }
   }
   return continue_status;
}


void
graphics_info_t::setup_draw_for_boids() {

   if (boids.size() == 0) {
      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?

      unsigned int n_boids = 30;
      boids.make_boids(n_boids);

      meshed_generic_display_object m;
      coot::colour_holder col(0.4, 0.5, 0.6);
      std::pair<glm::vec3, glm::vec3> start_end(glm::vec3(0.95,0,0), glm::vec3(-0.95,0,0));
      m.add_cone(start_end, col, 1.0, 0.0, 12, false, true,
                 meshed_generic_display_object::FLAT_CAP,
                 meshed_generic_display_object::FLAT_CAP);
      mesh_for_boids = m.mesh;

      std::vector<glm::mat4>    mats(n_boids);
      std::vector<glm::vec4> colours(n_boids);
      for (unsigned int i=0; i<n_boids; i++) {
         const fun::boid boid = boids[i];
         mats[i] = glm::mat4(1.0f);
         colours[i] = glm::vec4(0.2, 0.6, 0.4, 1.0);
      }
      Material material;
      mesh_for_boids.setup_rtsc_instancing(&shader_for_instanced_objects,
                                           mats, colours, n_boids, material);

      do_tick_boids = true;
      int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);

      // boids box

      std::vector<s_generic_vertex> vertices;
      std::vector<unsigned int> indices;
      glm::vec3 n(0,0,1);
      glm::vec4 c(0.3f, 0.3f, 0.3f, 1.0f);
      float boids_box_lim = boids.boids_box_limit;

      float corners[8][3] = {
                          {0,0,0}, //0
                          {0,0,1}, //1
                          {0,1,0}, //2
                          {0,1,1}, //3
                          {1,0,0}, //4
                          {1,0,1}, //5
                          {1,1,0}, //6
                          {1,1,1}};//7
      for (unsigned int i=0; i<8; i++) {
         for (unsigned int j=0; j<3; j++) {
            corners[i][j] *= 2.0f;
            corners[i][j] -= 1.0;
            corners[i][j] *= boids_box_lim;
         }
      }
      for (unsigned int ii=0; ii<8; ii++)
         vertices.push_back(s_generic_vertex(glm::vec3(corners[ii][0],corners[ii][1],corners[ii][2]), n, c));

      indices.push_back(0); indices.push_back(1);
      indices.push_back(1); indices.push_back(3);
      indices.push_back(3); indices.push_back(2);
      indices.push_back(2); indices.push_back(0);

      indices.push_back(4); indices.push_back(5);
      indices.push_back(5); indices.push_back(7);
      indices.push_back(7); indices.push_back(6);
      indices.push_back(6); indices.push_back(4);

      indices.push_back(0); indices.push_back(4);
      indices.push_back(1); indices.push_back(5);
      indices.push_back(2); indices.push_back(6);
      indices.push_back(3); indices.push_back(7);

      lines_mesh_for_boids_box = LinesMesh(vertices, indices);
      lines_mesh_for_boids_box.setup(&shader_for_lines);
   }
}

void
graphics_info_t::draw_ligand_view() {

   GtkAllocation allocation;
   gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
   float w = allocation.width;
   float h = allocation.height;
   graphics_ligand_mesh_molecule.draw(&shader_for_ligand_view,
                                      &shader_for_hud_geometry_tooltip_text,
                                      w, h, ft_characters);
}



void
graphics_info_t::draw_boids() {

   if (boids.size() > 0) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      mesh_for_boids.draw(&shader_for_instanced_objects,
                          mvp, view_rotation_matrix, lights, eye_position, bg_col,
                          shader_do_depth_fog_flag);

      lines_mesh_for_boids_box.draw(&shader_for_lines, mvp, view_rotation_matrix);
   }
}

void
graphics_info_t::draw_hydrogen_bonds_mesh() {

   if (mesh_for_hydrogen_bonds.draw_this_mesh) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      glm::vec4 bg_col(background_colour, 1.0);

      mesh_for_hydrogen_bonds.draw(&shader_for_instanced_objects,
                                   mvp, view_rotation_matrix, lights, eye_position, bg_col,
                                   shader_do_depth_fog_flag);
   }
}


std::vector<glm::vec3>
graphics_info_t::residue_to_positions(mmdb::Residue *residue_p) const {
   std::vector<glm::vec3> v;
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for(int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         glm::vec3 p(at->x, at->y, at->z);
         v.push_back(p);
      }
   }
   return v;
};

#include "pulse-data.hh"

void
graphics_info_t::setup_delete_item_pulse(mmdb::Residue *residue_p) {

   // next you use this functionn make it a member of graphics_info_t
   // gboolean delete_item_pulse_func(GtkWidget *widget,
   //                                 GdkFrameClock *frame_clock,
   //                                 gpointer data)
   // 
   auto delete_item_pulse_func = [] (GtkWidget *widget,
                                     GdkFrameClock *frame_clock,
                                     gpointer data) {

                                    gboolean continue_status = 1;
                                    pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
                                    pulse_data->n_pulse_steps += 1;
                                    if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
                                       continue_status = 0;
                                       lines_mesh_for_delete_item_pulse.clear();
                                       delete_item_pulse_centres.clear();
                                    } else {
                                       float ns = pulse_data->n_pulse_steps;
                                       lines_mesh_for_delete_item_pulse.update_buffers_for_pulse(ns, -1);
                                    }
                                    graphics_draw();
                                    return gboolean(continue_status);
                                 };

   pulse_data_t *pulse_data = new pulse_data_t(0, 20); // 20 matches the number in update_buffers_for_pulse()
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> positions = residue_to_positions(residue_p);
   delete_item_pulse_centres = positions;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = true;
   lines_mesh_for_delete_item_pulse.setup_pulse(&shader_for_lines_pulse, broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], delete_item_pulse_func, user_data, NULL);

};

void
graphics_info_t::setup_delete_residues_pulse(const std::vector<mmdb::Residue *> &residues) {

   // next you use this functionn make it a member of graphics_info_t
   // gboolean delete_item_pulse_func(GtkWidget *widget,
   //                                 GdkFrameClock *frame_clock,
   //                                 gpointer data)
   // 
   auto delete_item_pulse_func = [] (GtkWidget *widget,
                                     GdkFrameClock *frame_clock,
                                     gpointer data) {

                                    gboolean continue_status = 1;
                                    pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
                                    pulse_data->n_pulse_steps += 1;
                                    if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
                                       continue_status = 0;
                                       lines_mesh_for_delete_item_pulse.clear();
                                       delete_item_pulse_centres.clear();
                                    } else {
                                       float ns = pulse_data->n_pulse_steps;
                                       lines_mesh_for_delete_item_pulse.update_buffers_for_pulse(ns, -1);
                                    }
                                    graphics_draw();
                                    return gboolean(continue_status);
                                 };

   pulse_data_t *pulse_data = new pulse_data_t(0, 20); // 20 matches the number in update_buffers_for_pulse()
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> all_positions;
   for (unsigned int i=0; i<residues.size(); i++) {
      mmdb::Residue *residue_p = residues[i];
      std::vector<glm::vec3> residue_positions = residue_to_positions(residue_p);
      all_positions.insert(all_positions.end(), residue_positions.begin(), residue_positions.end());
   }
   delete_item_pulse_centres = all_positions;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = true;
   lines_mesh_for_delete_item_pulse.setup_pulse(&shader_for_lines_pulse, broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], delete_item_pulse_func, user_data, NULL);

};

// static
gboolean
graphics_info_t::invalid_residue_pulse_function(GtkWidget *widget,
                                                GdkFrameClock *frame_clock,
                                                gpointer data) {

   gboolean continue_status = 1;
   pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
   pulse_data->n_pulse_steps += 1;
   if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
      continue_status = 0;
      lines_mesh_for_identification_pulse.clear();
      delete_item_pulse_centres.clear(); // we sneakily use this vector (but no longer)
   } else {
      float ns = pulse_data->n_pulse_steps;
      lines_mesh_for_identification_pulse.update_buffers_for_invalid_residue_pulse(ns);
   }
   graphics_draw();
   return gboolean(continue_status);
}

void
graphics_info_t::setup_invalid_residue_pulse(mmdb::Residue *residue_p) {

   pulse_data_t *pulse_data = new pulse_data_t(0, 24);
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> residue_positions = residue_to_positions(residue_p);
   delete_item_pulse_centres = residue_positions; // sneakily use a wrongly named function
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = false;
   lines_mesh_for_identification_pulse.setup_pulse(&shader_for_lines_pulse, broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], invalid_residue_pulse_function, user_data, NULL);

}


void
graphics_info_t::draw_invalid_residue_pulse() {

   if (! lines_mesh_for_identification_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      glLineWidth(3.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_invalid_residue_pulse() glLineWidth " << err << std::endl;
      for (auto pulse_centre : delete_item_pulse_centres)
         lines_mesh_for_identification_pulse.draw(&shader_for_lines_pulse,
                                                  pulse_centre, mvp,
                                                  view_rotation_matrix, true);
   }
}


void
graphics_info_t::draw_identification_pulse() {

   if (! lines_mesh_for_identification_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      glLineWidth(2.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_identification_pulse() glLineWidth " << err << std::endl;
      lines_mesh_for_identification_pulse.draw(&shader_for_lines_pulse,
                                               identification_pulse_centre,
                                               mvp, view_rotation_matrix, true);
   }
}

void
graphics_info_t::draw_delete_item_pulse() {

   if (! lines_mesh_for_delete_item_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      glLineWidth(2.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_delete_item_pulse() glLineWidth " << err << std::endl;
      for (unsigned int i=0; i<delete_item_pulse_centres.size(); i++) {
         lines_mesh_for_delete_item_pulse.draw(&shader_for_lines_pulse,
                                               delete_item_pulse_centres[i],
                                               mvp, view_rotation_matrix, true);
      }
   }
}


void
graphics_info_t::move_forwards() {
   // these are the other way round in perspective - that's interesting.
   translate_in_screen_z(3.0);
}

void
graphics_info_t::move_backwards() {
   translate_in_screen_z(-3.0);
}

#include <glm/gtx/rotate_vector.hpp>
#include "matrix-utils.hh"

void
graphics_info_t::setup_key_bindings() {

   graphics_info_t g;

   // if we are serious about user-defined key-bindings all of these functions should be thunks in the user API
   // (and returning gboolean).

   auto l1 = []() { graphics_info_t g; g.adjust_clipping(-0.1); return gboolean(TRUE); };
   auto l2 = []() { graphics_info_t g; g.adjust_clipping( 0.1); return gboolean(TRUE); };
   auto l5 = []() { graphics_info_t g; g.blob_under_pointer_to_screen_centre(); return gboolean(TRUE); };

   auto l6 = []() {
                if (idle_function_spin_rock_token != -1) {
                   std::cout << "Removing the idle function\n";
                   gtk_widget_remove_tick_callback(glareas[0], idle_function_spin_rock_token);
                   idle_function_spin_rock_token = -1;
                } else {

                   do_tick_spin = true;
                   int spin_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);

                   // this is not a good name if we are storing a generic tick function id.
                   idle_function_spin_rock_token = spin_tick_id;
                }
                return gboolean(TRUE);
             };

   auto l7 = []() {
                int imol_scroll = graphics_info_t::scroll_wheel_map;
                if (graphics_info_t::is_valid_map_molecule(imol_scroll))
                   graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
                if (graphics_info_t::glareas.size() > 0)
                   int contour_idle_token = g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
                graphics_info_t g;
                g.set_density_level_string(imol_scroll, graphics_info_t::molecules[imol_scroll].contour_level);
                graphics_info_t::display_density_level_this_image = 1;
                return gboolean(TRUE);
             };

   auto l8 = []() {
                int imol_scroll = graphics_info_t::scroll_wheel_map;
                if (graphics_info_t::is_valid_map_molecule(imol_scroll))
                   graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;
                if (graphics_info_t::glareas.size() > 0)
                   int contour_idle_token = g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
                graphics_info_t g;
                g.set_density_level_string(imol_scroll, graphics_info_t::molecules[imol_scroll].contour_level);
                graphics_info_t::display_density_level_this_image = 1;
                return gboolean(TRUE);
             };

   auto l9 = []() { update_go_to_atom_from_current_position(); return gboolean(TRUE); };

   auto l10 = []() { graphics_info_t::zoom *= 0.9; return gboolean(TRUE); };

   auto l11 = []() { graphics_info_t::zoom *= 1.1; return gboolean(TRUE); };

   auto l12 = []() { graphics_info_t g; g.move_forwards(); return gboolean(TRUE); };

   auto l13 = []() { graphics_info_t g; g.move_backwards(); return gboolean(TRUE); };

   auto l14 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('forward')"); return gboolean(TRUE); };

   auto l15 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('backward')"); return gboolean(TRUE); };

   auto l16 = []() { graphics_info_t g; g.undo_last_move(); return gboolean(TRUE); };

   auto l18 = []() { graphics_info_t g; g.accept_moving_atoms(); return gboolean(TRUE); };

   auto l19 = []() { graphics_info_t g; g.clear_up_moving_atoms_wrapper(); return gboolean(TRUE); };

   auto l20 = []() { graphics_info_t g; g.eigen_flip_active_residue(); return gboolean(TRUE); };

   auto l21 = []() { graphics_info_t g; g.try_label_unlabel_active_atom(); return gboolean(TRUE); };

   auto l22 = []() {
                  graphics_info_t g;
                  g.setup_draw_for_particles();
                  return gboolean(TRUE);
             };

   // boids
   auto l23 = [] () {
                 graphics_info_t g;
                 if (! graphics_info_t::do_tick_boids)
                    graphics_info_t::do_tick_boids = true;
                 else
                    graphics_info_t::do_tick_boids = false;

                 g.setup_draw_for_boids();

                 if (! graphics_info_t::do_tick_boids)
                    std::cout << "--------- key press ----------- do_tick_boids "
                           << graphics_info_t::do_tick_boids << std::endl;
                 return gboolean(TRUE);
              };

   auto l24 = [] () {
                 // using the C API
                 // do_add_terminal_residue(1); // waits for user click :-)
                 graphics_info_t g;
                 g.add_terminal_residue_using_active_atom();
                 return gboolean(TRUE);
      };

   auto l25 = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    int imol_map = g.imol_refinement_map;
                    if (residue_p) {
                       mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].fill_partial_residue(residue_spec, g.Geom_p(), imol_map);

                       // now refine that
                       int saved_state = g.refinement_immediate_replacement_flag;
                       g.refinement_immediate_replacement_flag = 1;
                       std::string alt_conf("");
                       std::vector<mmdb::Residue *> rs = { residue_p };
                       g.refine_residues_vec(imol, rs, alt_conf, mol);
                       g.conditionally_wait_for_refinement_to_finish();
                       g.accept_moving_atoms();
                       g.refinement_immediate_replacement_flag = saved_state;
                    }
                 }
                 return gboolean(TRUE);
              };


   auto l26 = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].delete_residue_sidechain(residue_spec);
                    }
                 }
                 return gboolean(TRUE);
              };

   auto l28 = [] () {

                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       std::string this_chain_id = residue_p->GetChainID();
                       coot::residue_spec_t residue_spec(residue_p);
                       std::vector<std::vector<std::string> > ghost_chains_sets = molecules[imol].ncs_ghost_chains();
                       unsigned int n_ghost_chain_sets = ghost_chains_sets.size();
                       for (unsigned int i=0; i<n_ghost_chain_sets; i++) {
                          const std::vector<std::string> &chain_ids = ghost_chains_sets[i];
                          if (std::find(chain_ids.begin(), chain_ids.end(), this_chain_id) != chain_ids.end()) {
                             unsigned int idx_next = 0;
                             for (unsigned int j=0; j<chain_ids.size(); j++) {
                                if (chain_ids[j] == this_chain_id) {
                                   idx_next = j + 1;
                                   if (idx_next == chain_ids.size())
                                      idx_next = 0;
                                   break;
                                }
                             }
                             std::string chain_id_next = chain_ids[idx_next];
                             clipper::Coord_orth current_position = coot::co(at);
                             bool forward_flag = true;
                             glm::mat4 quat_mat = glm::toMat4(glm_quat);
                             clipper::Mat33<double> current_view_mat = glm_to_mat33(quat_mat);

                             if (molecules[imol].ncs_ghosts_have_rtops_p() == 0)
                                molecules[imol].fill_ghost_info(1, ncs_homology_level);

                             std::pair<bool, clipper::RTop_orth> new_ori =
                                molecules[imol].apply_ncs_to_view_orientation(current_view_mat,
                                                                              current_position,
                                                                              this_chain_id, chain_id_next,
                                                                              forward_flag);
                             if (new_ori.first) {
                                coot::util::quaternion q(new_ori.second.rot());
                                glm::quat q_ncs = coot_quaternion_to_glm(q);
                                glm_quat = glm::normalize(glm_quat * q_ncs); // wrong
                                clipper::Coord_orth t(new_ori.second.trn());
                                set_rotation_centre(t);

                                glm::quat q_ncs_1 = glm::rotate(q_ncs,   3.1415926f, glm::vec3(1,0,0));
                                glm::quat q_ncs_2 = glm::rotate(q_ncs_1, 3.1415926f, glm::vec3(1,0,0));
                                glm::quat q_ncs_3 = glm::inverse(q_ncs_2);

                                coot::util::quaternion cq = glm_to_coot_quaternion(q_ncs_3);

                                std::cout << "debug q_ncs  : " << glm::to_string(q_ncs)   << std::endl;
                                std::cout << "debug q_ncs_1: " << glm::to_string(q_ncs_1) << std::endl;
                                std::cout << "debug q_ncs_2: " << glm::to_string(q_ncs_2) << std::endl;
                                std::cout << "debug q_ncs_3: " << glm::to_string(q_ncs_3) << std::endl;
                                std::cout << "before: " << q << " after " << cq << std::endl;

                                graphics_info_t g;
                                g.update_things_on_move();

                             }
                             break;
                          }
                       }
                    } else {
                       std::cout << "ERROR:: no residue" << std::endl;
                    }
                 }
                 graphics_draw();
                 return gboolean(TRUE);
              };

   // Note to self, Space and Shift Space are key *Release* functions

   std::vector<std::pair<keyboard_key_t, key_bindings_t> > kb_vec;
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_d,      key_bindings_t(l1, "increase clipping")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_f,      key_bindings_t(l2, "decrease clipping")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_g,      key_bindings_t(l5, "go to blob")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_i,      key_bindings_t(l6, "spin")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_plus,   key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_equal,  key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_minus,  key_bindings_t(l7, "decrease contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_p,      key_bindings_t(l9, "update go-to atom by position")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_n,      key_bindings_t(l10, "Zoom in")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_m,      key_bindings_t(l11, "Zoom out")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_w,      key_bindings_t(l12, "Move forward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_s,      key_bindings_t(l13, "Move backward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l14, "NCS Skip forward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_O,      key_bindings_t(l15, "NCS Skip backward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_u,      key_bindings_t(l16, "Undo Move")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Return, key_bindings_t(l18, "Accept Moving Atoms")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Escape, key_bindings_t(l19, "Reject Moving Atoms")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_l,      key_bindings_t(l21, "Label/Unlabel Active Atom")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_q,      key_bindings_t(l22, "Particles")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_b,      key_bindings_t(l23, "Murmuration")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_y,      key_bindings_t(l24, "Add Terminal Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_k,      key_bindings_t(l25, "Fill Partial Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_K,      key_bindings_t(l26, "Delete Sidechain")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l28, "NCS Other Chain")));

   // control keys

   auto lc1 = []() { show_go_to_residue_keyboarding_mode_window(); return gboolean(TRUE); };
   key_bindings_t go_to_residue_key_binding(lc1, "Show Go To Residue Keyboarding Window");
   std::pair<keyboard_key_t, key_bindings_t> p1(keyboard_key_t(GDK_KEY_g, true), go_to_residue_key_binding);
   kb_vec.push_back(p1);

   auto lc2 = []() { graphics_info_t g; g.apply_undo(); return gboolean(TRUE); };
   key_bindings_t undo_key_binding(lc2, "Undo");
   std::pair<keyboard_key_t, key_bindings_t> p2(keyboard_key_t(GDK_KEY_z, true), undo_key_binding);
   kb_vec.push_back(p2);

   auto lc3 = []() { graphics_info_t g; g.apply_redo(); return gboolean(TRUE);};
   key_bindings_t redo_key_binding(lc3, "Redo");
   std::pair<keyboard_key_t, key_bindings_t> p3(keyboard_key_t(GDK_KEY_y, true), redo_key_binding);
   kb_vec.push_back(p3);

   auto ldr = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       // for this to work I need to move setup_delete_item_pulse() into
                       // graphics_info_t. Not today.
                       g.setup_delete_item_pulse(residue_p);
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].delete_residue(residue_spec);
                    }
                 }
                 return gboolean(TRUE);
              };
   key_bindings_t delete_residue_key_binding(ldr, "Delete Residue");
   std::pair<keyboard_key_t, key_bindings_t> pdel(keyboard_key_t(GDK_KEY_d, true), delete_residue_key_binding);
   kb_vec.push_back(pdel);


   // ctrl left
   auto lc4 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed. No need to test it again here.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Left);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Left);
                 } else {
                    keypad_translate_xyz(1, 1);
                 }
                 return gboolean(TRUE);
              };

   // ctrl right
   auto lc5 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Right);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Right);
                 } else {
                    keypad_translate_xyz(1, -1);
                 }
                 return gboolean(TRUE);
              };

   // ctrl up
   auto lc6 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Up);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Up);
                 } else {
                    keypad_translate_xyz(2, 1);
                 }
                 return gboolean(TRUE);
              };
   // ctrl down
   auto lc7 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Down);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Down);
                 } else {
                    keypad_translate_xyz(2, -1);
                 }
                 return gboolean(TRUE);
              };

   auto lc_qsa = [] () {
                    graphics_info_t g;
                    g.quick_save();
                    return gboolean(TRUE);
                 };

   key_bindings_t ctrl_arrow_left_key_binding(lc4, "R/T Left");
   key_bindings_t ctrl_arrow_right_key_binding(lc5, "R/T Right");
   key_bindings_t ctrl_arrow_up_key_binding(lc6, "R/T Up");
   key_bindings_t ctrl_arrow_down_key_binding(lc7, "R/T Down");
   key_bindings_t ctrl_eigen_flip(l20, "Eigen-Flip");
   key_bindings_t ctrl_quick_save(lc_qsa, "Quick Save");

   std::pair<keyboard_key_t, key_bindings_t> p4(keyboard_key_t(GDK_KEY_Left,  true), ctrl_arrow_left_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p5(keyboard_key_t(GDK_KEY_Right, true), ctrl_arrow_right_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p6(keyboard_key_t(GDK_KEY_Up,    true), ctrl_arrow_up_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p7(keyboard_key_t(GDK_KEY_Down,  true), ctrl_arrow_down_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p8(keyboard_key_t(GDK_KEY_e,     true), ctrl_eigen_flip);
   std::pair<keyboard_key_t, key_bindings_t> p9(keyboard_key_t(GDK_KEY_s,     true), ctrl_quick_save);

   kb_vec.push_back(p4);
   kb_vec.push_back(p5);
   kb_vec.push_back(p6);
   kb_vec.push_back(p7);
   kb_vec.push_back(p8);
   kb_vec.push_back(p9);

   std::vector<std::pair<keyboard_key_t, key_bindings_t> >::const_iterator it;
   for (it=kb_vec.begin(); it!=kb_vec.end(); it++)
     g.key_bindings_map[it->first] = it->second;

}


void
graphics_info_t::contour_level_scroll_scrollable_map(int direction) {

   int imol_scroll = scroll_wheel_map;
   if (! is_valid_map_molecule(imol_scroll)) {

      std::vector<int> dm = displayed_map_imols();
      if (std::find(dm.begin(), dm.end(), imol_scroll) == dm.end()) {
         if (dm.size() > 0)
            imol_scroll = dm[0];
      }
   }

   if (is_valid_map_molecule(imol_scroll)) {
      if (! molecules[imol_scroll].is_displayed_p()) {
         // don't scroll the map if the map is not displayed. Scroll the
         // map that *is* displayed
         std::vector<int> dm = displayed_map_imols();
         if (dm.size() > 0)
            imol_scroll = dm[0];
      }
   }

   if (is_valid_map_molecule(imol_scroll)) {
      // use direction
      if (direction == 1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
      if (direction == -1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;

      // std::cout << "INFO:: contour level for map " << imol_scroll << " is "
      // << molecules[imol_scroll].contour_level << std::endl;
      set_density_level_string(imol_scroll, molecules[imol_scroll].contour_level);
      display_density_level_this_image = 1;

      graphics_draw(); // queue
   }
}
