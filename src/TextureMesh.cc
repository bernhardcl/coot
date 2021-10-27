
#ifdef USE_PYTHON
#include "Python.h"
#endif

#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>

#include "ft-character.hh"
#include "TextureMesh.hh"

#ifdef THIS_IS_HMT
#include "display-info.hh"
#else
#include "graphics-info.h"
#endif

#include <Texture.hh> // experimental

// for the moment make the scales explict, when fixed make the scales default
void
TextureMesh::setup_camera_facing_quad(Shader *shader_p, float scale_x, float scale_y) {

   shader_p->Use();

   draw_this_mesh = true;

   glm::vec3 n(0,0,1);
   glm::vec4 col(1.0, 1.0, 1.0, 1.0);

   vertices.clear();
   triangles.clear();

   // the indexing might well be wrong here - I'm sort of guessing
   vertices.push_back(TextureMeshVertex(glm::vec3(-scale_x,  scale_y, 0.0f), n, col, glm::vec2(0,0)));
   vertices.push_back(TextureMeshVertex(glm::vec3( scale_x,  scale_y, 0.0f), n, col, glm::vec2(1,0)));
   vertices.push_back(TextureMeshVertex(glm::vec3( scale_x, -scale_y, 0.0f), n, col, glm::vec2(1,1)));
   vertices.push_back(TextureMeshVertex(glm::vec3(-scale_x, -scale_y, 0.0f), n, col, glm::vec2(0,1)));

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));

   setup_buffers();

}

void
TextureMesh::set_colour(const glm::vec4 &col_in) {

   for (unsigned int i=0; i<vertices.size(); i++) {
      vertices[i].color = col_in;
   }
}

void
TextureMesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   glGenVertexArrays(1, &vao);
   glBindVertexArray(vao);

   glGenBuffers(1, &buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   unsigned int n_vertices = vertices.size();
   glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_STATIC_DRAW);

   // position
   glEnableVertexAttribArray(0);
   glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex), 0);

   // normal
   glEnableVertexAttribArray (1);
   glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                          reinterpret_cast<void *>(sizeof(glm::vec3)));

   // colour
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                         reinterpret_cast<void *>(2 * sizeof(glm::vec3)));

   // texture coordinates
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                         reinterpret_cast<void *>(2 * sizeof(glm::vec3) + sizeof(glm::vec4)));


   glGenBuffers(1, &index_buffer_id);
   GLenum err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";
   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);
   if (false)
      std::cout << "debug:: glBufferData for index buffer_id " << index_buffer_id
                << " n_triangles: " << n_triangles
                << " allocating with size: " << n_bytes << " bytes" << std::endl;
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_STATIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error setup_simple_triangles()\n";

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

}


void
TextureMesh::draw_atom_label(const std::string &atom_label,
                             const glm::vec3 &atom_label_position,
                             const glm::vec4 &text_colour, // set using glBufferSubData
                             Shader *shader_p,
                             const glm::mat4 &mvp,
                             const glm::mat4 &view_rotation_matrix,
                             const glm::vec4 &background_colour,
                             bool do_depth_fog,
                             bool is_perspective_projection) {

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "   error draw_atom_label() " << shader_p->name << " -- start -- error "
                      << err << std::endl;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   shader_p->set_vec3_for_uniform("label_position", atom_label_position);

   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post set label_position " << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);

   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post set mvp " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);

   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post set view rotation " << err << std::endl;

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post background_colour " << err << std::endl;
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_atom_label() " << name << " " << shader_p->name
                      << " post do_depth_fog " << err << std::endl;
   shader_p->set_bool_for_uniform("is_perspective_projection", is_perspective_projection);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_atom_label() " << name << " " << shader_p->name
                      << " post is_perspective_projection " << err << std::endl;

   if (vao == 99999999)
      std::cout << "You forget to setup this TextureMesh " << name << " " << shader_p->name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error TextureMesh::draw_atom_label()) " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw_atom_label() A3 " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "error TextureMesh::draw_atom_label() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "error TextureMesh::draw_atom_label() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_atom_label() " << name << " pre-draw " << err << std::endl;

   // scale of 20 make letter separation good
   GLfloat scale = 100.1;

   // ------------------------------- text texture code here -----------------------

   // consider refactoring when working

   // display_info_t di;

   float x = 0;
   float y = 0;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   err = glGetError(); if (err) std::cout << "GL ERROR:: draw_atom_label() A0 " << err << std::endl;

#ifdef THIS_IS_HMT
   std::map<GLchar, FT_character> &ft_characters = display_info_t::ft_characters;
#else
   std::map<GLchar, FT_character> &ft_characters = graphics_info_t::ft_characters;
#endif

   for (std::string::const_iterator it_c = atom_label.begin(); it_c != atom_label.end(); ++it_c) {
      // bitsy spider...
      std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*it_c);
      if (it == ft_characters.end()) {
         std::cout << "ERROR:: TextureMesh::draw_atom_label():: Failed to lookup glyph for " << *it_c
                   << " in " << atom_label << std::endl;
         continue;
      };
      const FT_character &ch = it->second;
      GLfloat xpos = x + ch.Bearing.x * scale;
      GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;
      GLfloat w = ch.Size.x * scale;
      GLfloat h = ch.Size.y * scale;

      // std::cout << "---------- char " << *c << " x " << x << std::endl;

      // Update the vertices for each character
      //
      std::vector<TextureMeshVertex> texture_mesh_vertices = vertices;

      // 0,0 -> add h
      // 1,0 -> add w,h
      // 1,1 -> add w

      bool debug = false;

      if (debug) {
         std::cout << "texture_mesh_vertices 0 " << glm::to_string(texture_mesh_vertices[0].position) << std::endl;
         std::cout << "texture_mesh_vertices 1 " << glm::to_string(texture_mesh_vertices[1].position) << std::endl;
         std::cout << "texture_mesh_vertices 2 " << glm::to_string(texture_mesh_vertices[2].position) << std::endl;
         std::cout << "texture_mesh_vertices 3 " << glm::to_string(texture_mesh_vertices[3].position) << std::endl;
         std::cout << "here with w " << w << " and h " << h << std::endl;
      }

      for (unsigned int i=0; i<4; i++)
         texture_mesh_vertices[i].color = text_colour;

      for (unsigned int i=0; i<4; i++)
         texture_mesh_vertices[i].position += glm::vec3(xpos, ypos, 0.0f);
      
      texture_mesh_vertices[0].position.y += h;
      texture_mesh_vertices[1].position.x += w;
      texture_mesh_vertices[1].position.y += h;
      texture_mesh_vertices[2].position.x += w;

      if (debug) {
         std::cout << "post texture_mesh_vertices 0 " << glm::to_string(texture_mesh_vertices[0].position) << std::endl;
         std::cout << "post texture_mesh_vertices 1 " << glm::to_string(texture_mesh_vertices[1].position) << std::endl;
         std::cout << "post texture_mesh_vertices 2 " << glm::to_string(texture_mesh_vertices[2].position) << std::endl;
         std::cout << "post texture_mesh_vertices 3 " << glm::to_string(texture_mesh_vertices[3].position) << std::endl;
      }

      glBindTexture(GL_TEXTURE_2D, ch.TextureID);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      unsigned int n_vertices = vertices.size();
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(TextureMeshVertex), &texture_mesh_vertices[0]);

      unsigned int n_draw_verts = 6;
      glDrawElements(GL_TRIANGLES, n_draw_verts, GL_UNSIGNED_INT, nullptr);

      err = glGetError(); if (err) std::cout << "draw_atom_label() glDrawArrays() " << err << std::endl;

       // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
      x += (ch.Advance >> 6) * scale * 1.0;

   }


   // ------------------------------- done text texture code  ----------------------

   err = glGetError();
   if (err) std::cout << "   error TextureMesh::draw() glDrawElements()"
                      << " of Mesh \"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);

   glUseProgram (0);
}

void
TextureMesh::draw(Shader *shader_p,
                  const glm::mat4 &mvp,
                  const glm::mat4 &view_rotation_matrix,
                  const std::map<unsigned int, lights_info_t> &lights,
                  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                  const glm::vec4 &background_colour,
                  bool do_depth_fog) {

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "   error draw() " << shader_p->name << " -- start -- " << err << std::endl;

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "   error:: " << shader_p->name << " draw() post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "   error:: " << shader_p->name << " draw() post view rotation uniform "
                      << err << std::endl;

   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " pre-set eye position "
                      << " with GL err " << err << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);
   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " post-set eye position "
                      << " with GL err " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   // this lights block can be in it's own function (same as Mesh)
   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   if (vao == 99999999)
      std::cout << "You forget to setup this mesh " << name << " "
                << shader_p->name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw() A3 " << err << std::endl;
   glActiveTexture(GL_TEXTURE1);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw() A4 " << err << std::endl;

   shader_p->set_int_for_uniform("base_texture", 0);
   shader_p->set_int_for_uniform("normal_map",   1);

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id); // 20211017-PE  needed?
   err = glGetError();
   if (err) std::cout << "   error draw() glBindBuffer() v " << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError();
   if (err) std::cout << "   error draw() glBindBuffer() i " << err << std::endl;

   glEnableVertexAttribArray(0);  // position
   glEnableVertexAttribArray(1);  // normal
   glEnableVertexAttribArray(2);  // colour (not used)
   glEnableVertexAttribArray(3);  // texture coordinates

   err = glGetError();
   if (err) std::cout << "   error draw() " << name << " pre-draw " << err << std::endl;

   // if (use_blending) {
   // glEnable(GL_BLEND);
   // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   //}

   // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
   // before making a new VAO?

   if (false)
      std::cout << "debug:: TextureMesh::draw() " << name << " shader " << shader_p->name
                << " vao " << vao
                << " drawing " << n_verts << " triangle vertices"  << std::endl;

   glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "   error TextureMesh::draw() glDrawElements()"
                      << " of Mesh \"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   // if (use_blending) {
   // glDisable(GL_BLEND);
   // }

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);

   glUseProgram (0);

}


void
TextureMesh::import(const std::vector<TextureMeshVertex> &verts_in, const std::vector<g_triangle> &triangles_in) {

   vertices = verts_in;
   triangles = triangles_in;
   draw_this_mesh = true;
}


void
TextureMesh::import(const IndexedModel &ind_model, float scale) {

   // std::cout << "TextureMesh::import(const IndexedModel &ind_model)" << std::endl;

   bool sane_input = false;
   if (ind_model.positions.size() == ind_model.texCoords.size())
      if (ind_model.positions.size() == ind_model.normals.size())
         sane_input = true;

   if (ind_model.positions.size() == 0)
      sane_input = false;

   std::cout << "TextureMesh::import() indices.size() " << ind_model.indices.size() << std::endl;

   if (sane_input) {
      for (unsigned int i=0; i<ind_model.positions.size(); i++) {
         glm::vec4 col(0.5, 0.5, 0.5, 1.0);
         // std::cout << "debug normal " << i << " " << glm::to_string(ind_model.normals[i]) << std::endl;
         // std::cout << "debug texCoords " << i << " " << glm::to_string(ind_model.texCoords[i]) << std::endl;
         TextureMeshVertex v(scale * ind_model.positions[i],
                             ind_model.normals[i],
                             col,
                             ind_model.texCoords[i]);
         vertices.push_back(v);
      }

      for (unsigned int i=0; i<ind_model.indices.size(); i += 3) {
         g_triangle gt(ind_model.indices[i],
                       ind_model.indices[i+1],
                       ind_model.indices[i+2]);
         triangles.push_back(gt);
      }
   }

   setup_buffers();

}

// for happy faces that drift up the screen
void
TextureMesh::update_instancing_buffer_data_for_happy_faces(const std::vector<glm::vec3> &positions_in, // original positions
                                           unsigned int draw_count_in,
                                           unsigned int draw_count_max,
                                           const glm::vec3 &screen_y_uv) {

   glBindVertexArray(vao);
   draw_count = draw_count_in; // do I need this to be a class data item?
   std::vector<glm::vec3> positions(positions_in);
   int n_positions = positions.size();
   if (n_positions > n_instances_allocated) {
      std::cout << "Too many TextureMesh instances " << n_positions << " " << n_instances_allocated
                << std::endl;
   } else {

      // Use the screen centre to generate tp.
      // Use the index of the position to change the phaase of the wiggle.

      auto get_position_delta = [draw_count_max] (const glm::vec3 &screen_y_uv,
                                                  unsigned int draw_count_in,
                                                  unsigned int index) {
                                   float f1 = static_cast<float>(draw_count_in)/static_cast<float>(draw_count_max);
                                   float f2 = f1 * f1 * 2.5f;
                                   glm::vec3 f_uv = f2 * screen_y_uv;
                                   glm::vec3 tp = glm::normalize(glm::vec3(0.1, 0.2, 0.3));
                                   glm::vec3 cp_1 = glm::cross(screen_y_uv, tp);
                                   glm::vec3 cp_2 = glm::cross(screen_y_uv, cp_1);
                                   float phase = 0.1 * static_cast<float>(index);
                                   f_uv += 0.9 * sinf(9.0 * f1 + phase) * cp_2;
                                   return f_uv;
                                };

      n_instances = positions.size();

      // now update the positions
      for (unsigned int i=0; i<positions.size(); i++)
         positions[i] += get_position_delta(screen_y_uv, draw_count_in, i);

      glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_positions * sizeof(glm::vec3), &(positions[0]));
   }
}

void
TextureMesh::setup_instancing_buffers(unsigned int n_happy_faces_max) {

   n_instances = 0;
   n_instances_allocated = n_happy_faces_max;
   is_instanced = true;

   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "GL error ####"
                      << " TextureMesh::setup_instancing_buffers() A " << err << std::endl;
   unsigned int n_bytes = n_happy_faces_max * sizeof(glm::vec3);

   glGenBuffers(1, &inst_positions_id);
   // std::cout << "setup_instancing_buffers: glGenBuffers inst_positions_id " << inst_positions_id << std::endl;
   glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
   glBufferData(GL_ARRAY_BUFFER, n_bytes, nullptr, GL_DYNAMIC_DRAW);
   // prevous attributes are position, normal, colour, texCoords
   glEnableVertexAttribArray(4);

   void *step_over_previous_attrib_bytes = 0; // for instanced attributes, we don't need to step
                                              // over the standard vertex attributes.

   glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), step_over_previous_attrib_bytes);
   glVertexAttribDivisor(4, 1);
   err = glGetError();
   if (err) std::cout << "GL error #####"
                      << " TextureMesh::setup_instancing_buffers() B " << err << std::endl;
}


void
TextureMesh::draw_instances(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation,
                            unsigned int draw_count, unsigned int draw_count_max) {

   // std::cout << "TextureMesh::draw_instances() A " << n_instances << " " << triangles.size()
   // <<std::endl;

   if (! draw_this_mesh) return;
   // this can happen when all the particles have life 0 - and have been removed.
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   const float pi = 3.1415926f;
   float draw_count_frac = static_cast<float>(draw_count)/static_cast<float>(draw_count_max);
   const float &f = draw_count_frac;  // shorthand
   float opacity = sinf(sqrt(f) * pi); // maybe this goes negative?
   // std::cout << "opacity: f " << f << " opacity " << opacity << std::endl;

   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_instances() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // vertex colours
   glEnableVertexAttribArray(3); // texture coords
   glEnableVertexAttribArray(4); // instanced position

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_instances() " << shader_p->name
                      << " draw_instances() post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_instances() " << shader_p->name
                      << " draw_instances() post view_rotation uniform " << err << std::endl;

   shader_p->set_float_for_uniform("opacity", opacity);

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw_instances() activetexture "
                                          << err << std::endl;

   glEnable(GL_DEPTH_TEST); // not the problem
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   unsigned int n_verts = 6;
   // std::cout << "TextureMesh::draw_instances() C " << n_verts << " " << n_instances << std::endl;
   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);

   err = glGetError();
   if (err) std::cout << "error draw_instances() on glDrawElementsInstanced() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   
}

// Define these only in *one* .cc file.
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "tiny_gltf.h"

bool
TextureMesh::load_from_glTF(const std::string &file_name_in, bool include_call_to_setup_buffers) {

   // is this what's in a node?
   class extracted_buffer_info_t {
   public:
      int buffer_view_index;
      std::vector<glm::vec3> positions;
      std::vector<glm::vec3> normals;
      std::vector<glm::vec2> texture_coords;
      std::vector<g_triangle> triangles;
      glm::vec4 base_colour;
      int normal_map_texture_image_index; // NormalTextureInfo is part of a Material. It tells us
                                          // which image has the normal info (if any)
   };


   auto indent = [] (unsigned int il) {
                   std::string s;
                   for (unsigned int i=0; i<il; i++) s += "   ";
                   return s;
                   };

   auto proc_indices = [indent] (tinygltf::Model &model, const tinygltf::Accessor &indices_accessor) {

                          std::vector<g_triangle> triangles;
                          std::cout << indent(3) << "proc_indices()" << std::endl;

                          std::vector<unsigned int> mesh_indices;
                          const tinygltf::BufferView &buffer_view = model.bufferViews[indices_accessor.bufferView];
                          const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                          const uint8_t *base = &buffer.data.at(buffer_view.byteOffset + indices_accessor.byteOffset);

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                             const uint32_t *p = reinterpret_cast<const uint32_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                             const uint16_t *p = reinterpret_cast<const uint16_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_BYTE) {
                             const uint8_t *p = reinterpret_cast<const uint8_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          unsigned int n_triangles = mesh_indices.size()/3;
                          if (n_triangles * 3 == mesh_indices.size()) {
                             for (unsigned int ii=0; ii<mesh_indices.size(); ii+=3) {
                                g_triangle t(mesh_indices[ii], mesh_indices[ii+1], mesh_indices[ii+2]);
                                triangles.push_back(t);
                             }
                          }
                          return triangles;
                       };

   auto proc_material = [] (const tinygltf::Material &material) {

                           std::cout << "debug:: proc_material(): material.pbrMetallicRoughness() basecolour size"
                                     << material.pbrMetallicRoughness.baseColorFactor.size() << std::endl;
                           std::cout << "debug:: proc_material(): material.normalTexture:"
                                     << " index    " << material.normalTexture.index
                                     << " texCoord " << material.normalTexture.texCoord
                                     << " scale    " << material.normalTexture.scale
                                     << std::endl;

                           if (false)
                              std::cout << "Material name: \"" << material.name << "\""
                                        << " alphaMode" << material.alphaMode
                                        << " pbr-metallicroughness colour "
                                        << material.pbrMetallicRoughness.baseColorFactor[0] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[1] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[2] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[3] << " "
                                        << " metalicFactor " << material.pbrMetallicRoughness.metallicFactor
                                        << " roughnessFactor " << material.pbrMetallicRoughness.roughnessFactor
                                        << std::endl;
                           glm::vec4 colour(material.pbrMetallicRoughness.baseColorFactor[0],
                                            material.pbrMetallicRoughness.baseColorFactor[1],
                                            material.pbrMetallicRoughness.baseColorFactor[2],
                                            material.pbrMetallicRoughness.baseColorFactor[3]);
                           return colour;
                   };

   auto proc_primitive = [proc_indices, proc_material, indent] (tinygltf::Model &model, const tinygltf::Primitive &primitive) {

                            std::cout << indent(2) << "proc_primitive()" << std::endl;

                            extracted_buffer_info_t ebi; // returned object

                            // ------------------- primitive material -------------------

                            std::cout << indent(3) << "primitive material " << primitive.material << std::endl;
                            if (primitive.material >= 0) {
                               int material_index = primitive.material;
                               std::cout << "debug:: proc_primitive(): material_index " << material_index << " of "
                                         << " materials.size() " << model.materials.size() << std::endl;
                               const tinygltf::Material &material = model.materials[material_index];
                               ebi.base_colour = proc_material(material);
                               ebi.normal_map_texture_image_index = material.normalTexture.index;
                               // other maps are available, emissive, roughness, occlussion
                            } else {
                               std::cout << "Ooops skipping proc_material()" << std::endl;
                            }

                            // ------------------- primitive attributes -------------------

                            std::map<std::string, int>::const_iterator it;
                            unsigned int icount = 0;
                            for (it=primitive.attributes.begin(); it!=primitive.attributes.end(); ++it) {
                               const std::string &key = it->first;
                               std::cout << "::: key " << std::setw(8) << key << " for attribute index "
                                         << icount << " of " << primitive.attributes.size() << std::endl;
                               icount++;
                            }

                            for (it=primitive.attributes.begin(); it!=primitive.attributes.end(); ++it) {
                               const std::string &key = it->first;
                               int index = it->second;
                               std::cout << indent(3) << "primvitive attribute " << key << " " << index << std::endl;

                               const tinygltf::Accessor &accessor = model.accessors[index];
                               const tinygltf::BufferView &buffer_view = model.bufferViews[accessor.bufferView];
                               const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER)
                                  std::cout << indent(3) << "an array buffer" << std::endl;
                               if (buffer_view.target == TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER)
                                  std::cout << indent(3) << "an index buffer" << std::endl;

                               if (false)
                                  std::cout << indent(2) << "accessor componentType " << accessor.componentType
                                            << " accessor type " << accessor.type << std::endl;

                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_BYTE)
                                  std::cout << indent(3) <<"component type byte " << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE)
                                  std::cout << indent(3) <<"component type unsigned byte " << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_SHORT)
                                  std::cout << indent(3) <<"component type short" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT)
                                  std::cout << indent(3) <<"component type unsigned short" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_INT)
                                  std::cout << indent(3) <<"component type int" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT)
                                  std::cout << indent(3) <<"component type unsigned int" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT)
                                  std::cout << indent(3) <<"component type float" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE)
                                  std::cout << indent(3) <<"component type double" << std::endl;

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "POSITION") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                           // std::cout << "      ---- float array for positions " << std::endl;
                                           const float *positions = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = positions[ii*3 + 0];
                                              float y = positions[ii*3 + 1];
                                              float z = positions[ii*3 + 2];
                                              ebi.positions.push_back(glm::vec3(x,y,z));
                                           }
                                        }
                                     }
                                  }
                               }

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "NORMAL") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                           // std::cout << "      ---- float array for normals " << std::endl;
                                           const float *normals = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = normals[ii*3 + 0];
                                              float y = normals[ii*3 + 1];
                                              float z = normals[ii*3 + 2];
                                              ebi.normals.push_back(glm::vec3(x,y,z));
                                           }
                                        }
                                     }
                                  }
                               }

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "TEXCOORD_0") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC2) {
                                           std::cout << "      ---- float array for texture coordinates " << std::endl;
                                           const float *tcs = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = tcs[ii*2 + 0];
                                              float y = tcs[ii*2 + 1];
                                              ebi.texture_coords.push_back(glm::vec2(x,y));
                                           }
                                        }
                                     }
                                  }
                               }
                            }

                            // ------------------- primitive indices -------------------

                            const tinygltf::Accessor &indices_accessor = model.accessors[primitive.indices];
                            std::vector<g_triangle> triangles = proc_indices(model, indices_accessor);
                            ebi.triangles = triangles; // std::move here?

                            std::cout << indent(2) << "debug:: proc_primitive() returns ebi with " << ebi.positions.size() << " positions "
                                      << ebi.normals.size() << " normals " << ebi.texture_coords.size() << " texture-coords and "
                                      << triangles.size() << " triangles" << std::endl;

                            return ebi;
                         };

   auto proc_mesh = [proc_primitive, indent] (tinygltf::Model &model, const tinygltf::Mesh &mesh) {

                       std::cout << indent(1) << "proc_mesh(): " << mesh.name << std::endl;
                       std::vector<extracted_buffer_info_t> ebi_vec;
                       for (unsigned int i=0; i<mesh.primitives.size(); i++) {
                          std::cout << indent(1) << "proc_mesh(): primitive number " << i << " of " << mesh.primitives.size() << std::endl;
                          const tinygltf::Primitive &primitive = mesh.primitives[i];
                          if (primitive.indices >= 0) {
                             auto new_mesh_primitive = proc_primitive(model, primitive);
                             ebi_vec.push_back(new_mesh_primitive);
                          }
                       }
                       return ebi_vec;
                   };


   auto proc_node = [proc_mesh] (tinygltf::Model &model, size_t i_node_index, const tinygltf::Node &node) {

                       std::vector<extracted_buffer_info_t> r;

                       std::cout << "Node " << i_node_index << " info:" << std::endl;
                       std::cout << "   Node name: " << node.name << std::endl;
                       std::cout << "   Node mesh index: " << node.mesh << std::endl;
                       std::cout << "   Node rotation vec elements count: " << node.rotation.size() << std::endl;
                       std::cout << "   Node scale vec elements count:    " << node.scale.size() << std::endl;
                       std::cout << "   Node scale vec translation count: " << node.translation.size() << std::endl;
                       if (node.mesh > -1) { // not a light or a camera
                          r = proc_mesh(model, model.meshes[node.mesh]);
                       } else {
                       }
                       return r;
                   };

   auto proc_image = [indent] (const tinygltf::Image &image) {

                        std::cout << "Image:" << std::endl;
                        std::cout << indent(1) << "Image name: "   << image.name   << std::endl;
                        std::cout << indent(1) << "Image bits: "   << image.bits   << std::endl;
                        std::cout << indent(1) << "Image width: "  << image.width  << std::endl;
                        std::cout << indent(1) << "Image height: " << image.height << std::endl;
                        std::cout << indent(1) << "Image component: " << image.component << std::endl;
                        std::cout << indent(1) << "Image data buffer size:" << image.image.size() << std::endl;
                        std::cout << indent(1) << "Image.bufferView:" << image.bufferView << std::endl;

                        Texture texture;
                        texture.handle_raw_image_data(image.name, image.image, image.width, image.height);
                        std::cout << "debug:: made a texture with m_texture_handle " << texture.m_texture_handle << std::endl;

                        return texture;
                   };

   auto proc_model_v2 = [proc_node, proc_image] (tinygltf::Model &model) {

                           std::vector<extracted_buffer_info_t> r;
                           std::vector<Texture> textures;

                           unsigned int n_accessors = model.accessors.size();
                           std::cout << "debug:: model has " << n_accessors << " accessors" << std::endl;

                           // --- Materials ---

                           std::cout << "INFO:: this model contains " << model.materials.size() << " materials" << std::endl;
                           for (unsigned int imat=0; imat<model.materials.size(); imat++) {
                              const tinygltf::Material &material = model.materials[imat];
                              std::cout << "Material " << imat << " name: " << material.name << " alphaMode" << material.alphaMode
                                        << " pbr-metallicroughness colour "
                                        << material.pbrMetallicRoughness.baseColorFactor[0] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[1] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[2] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[3] << " "
                                        << " metalicFactor " << material.pbrMetallicRoughness.metallicFactor
                                        << " roughnessFactor " << material.pbrMetallicRoughness.roughnessFactor
                                        << std::endl;
                           }

                           // --- Vertices and Indices ---
#if 0
                           // scenes

                           for (unsigned int iscene=0; iscene<model.scenes.size(); iscene++) {
                              const tinygltf::Scene &scene = model.scenes[iscene];
                              std::cout << ":::: This scene contains " << scene.nodes.size() << " nodes "<< std::endl;
                              for (size_t i = 0; i < scene.nodes.size(); i++) {
                                 auto nodes = proc_node(model, i, model.nodes[scene.nodes[i]]);
                                 r.insert(r.end(), nodes.begin(), nodes.end());
                              }
                           }
#endif

                           // model

                           std::cout << "This model contains " << model.nodes.size() << " nodes" << std::endl;

                           for (size_t i = 0; i < model.nodes.size(); i++) {
                              auto nodes = proc_node(model, i, model.nodes[i]);
                              r.insert(r.end(), nodes.begin(), nodes.end());
                           }

                           // images

                           std::cout << "This model contains " << model.images.size() << " images" << std::endl;
                           for (unsigned int i=0; i<model.images.size(); i++) {
                              Texture texture = proc_image(model.images[i]);
                              textures.push_back(texture);
                           }

                           return std::make_pair(r, textures);
                        };

   bool status = true;
   tinygltf::Model model;
   tinygltf::TinyGLTF loader;
   std::string err;
   std::string warn;

   std::string file_name = file_name_in;
   if (coot::file_exists(file_name)) {
      // do nothing
   } else {
      std::string dir = coot::package_data_dir();
      std::string dir_2 = coot::util::append_dir_dir(dir, "glTF");
      file_name = coot::util::append_dir_file(dir_2, file_name);
   }

   bool use_binary = false;
   std::string ext = coot::util::file_name_extension(file_name_in);
   std::cout << "debug:: ext is " << ext << std::endl;
   if (ext == ".glb") use_binary = true;

   // std::cout << "debug:: TextureMesh::load_from_glTF(): " << file_name_in << " use_binary: " << use_binary << std::endl;

   // use the extension to check which function to use
   //
   // bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
   bool read_success = false;

   if (use_binary) {
      read_success = loader.LoadBinaryFromFile(&model, &err, &warn, file_name); // for binary glTF(.glb)
   } else {
      read_success = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
   }

   if (!warn.empty()) {
      std::cout << "WARNING:: load_from_glTF():" << warn << std::endl;
   }

   if (!err.empty()) {
      std::cout << "ERROR:: load_from_glTF(): " << err << std::endl;
   }

   if (read_success == false) {
      std::cout << "WARNING:: failed to parse glTF from " << file_name << std::endl;
   }

   if (err.empty() && read_success) { // OK, good, let's go!

      auto r_pair = proc_model_v2(model);
      std::vector<extracted_buffer_info_t> &r = r_pair.first;
      std::vector<Texture> &extracted_textures = r_pair.second;

      if (! extracted_textures.empty()) {
         for (unsigned int idx_texture=0; idx_texture<extracted_textures.size(); idx_texture++) {
            const auto &texture = extracted_textures[idx_texture];
            int idx_texture_as_int = idx_texture; // sigh
            TextureInfoType ti(texture, "something", "base_texture", 0);
            // now, was this marked as a normal map by any of the buffer infos?
            for (unsigned int i=0; i<r.size(); i++) {
               const extracted_buffer_info_t &ebi = r[i];
               if (ebi.normal_map_texture_image_index == idx_texture_as_int) {
                  ti.texture_type = TextureInfoType::NORMAL_MAP;
                  ti.unit = 1; // for rendering using textures-meshes.shader
                  ti.sampler_name = "normal_map";
               }
            }
            std::cout << "load_from_glTF() pushing back a textureinfo " << ti.name << " " << ti.sampler_name << " unit "
                      << ti.unit << " texture_type_t " << ti.texture_type << " filename " << ti.texture.file_name << std::endl;
            textures.push_back(ti);
         }
      }

      // std::cout << "load_from_glTF(): found " << r.size() << " mesh primitives" << std::endl;
      for (unsigned int i=0; i<r.size(); i++) {
         const extracted_buffer_info_t &ebi = r[i];
         if (ebi.normals.size() > 0) {
            if (ebi.normals.size() == ebi.positions.size()) {
               unsigned int max_vertex_index = 0; // set this
               if (max_vertex_index < ebi.normals.size()) {
                  unsigned int idx_vert_base = vertices.size();
                  unsigned int idx_tri_base = triangles.size();
                  std::cout << "debug ebi sizes " << ebi.positions.size() << " " << ebi.normals.size() << " " << ebi.texture_coords.size()
                            << std::endl;
                  if (ebi.texture_coords.size() == ebi.positions.size()) {
                     for (unsigned int j=0; j<ebi.normals.size(); j++) {
                        // s_generic_vertex g(ebi.positions[j], ebi.normals[j], ebi.base_colour);
                        TextureMeshVertex tmv(ebi.positions[j], ebi.normals[j], ebi.base_colour, ebi.texture_coords[j]);
                        vertices.push_back(tmv);
                     }
                  }

                  // 20211011-PE If ebi.texture_coords.size() is 0 - as we would expect for a untextured model, then perhaps we can
                  // just add a dummy texture_coords and set a flag to say to use the colours and not the textures.
                  // Another time.

                  triangles.insert(triangles.end(), ebi.triangles.begin(), ebi.triangles.end());
                  if (idx_vert_base != 0) {
                     for (unsigned int i=idx_tri_base; i<triangles.size(); i++) {
                        triangles[i].rebase(idx_vert_base);
                     }
                  }
               }
            }
         }
      }

      if (include_call_to_setup_buffers) {
         std::cout << "pre-setup_buffer() " << vertices.size() << " vertices "  << triangles.size() << " triangles "  << std::endl;
         setup_buffers();
      }
   } else {
      std::cout << "::::::::::::: non-success path" << std::endl;
      status = false; // boo
   }

   std::cout << "load_from_glTF() returns status " << status << std::endl;
   return status;
}

