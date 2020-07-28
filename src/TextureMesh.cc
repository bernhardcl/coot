
#ifdef USE_PYTHON
#include "Python.h"
#endif

#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>

#include "ft-character.hh"
#include "TextureMesh.hh"

#ifdef THIS_IS_HMT
#include "display-info.hh"
#else
#include "graphics-info.h"
#endif

void
TextureMesh::setup_camera_facing_quad(const Shader &sh) {

   float scale_x = 4.4; // pass?
   float scale_y = 1.2;

   // I don't understand why these don't affect the size of the atom label;

   // scale_x = 0.00001;
   // scale_y = 0.00001;

   shader = sh;
   shader.Use();

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
   if (true)
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

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   shader_p->set_vec3_for_uniform("label_position", atom_label_position);

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);
   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);
   err = glGetError();
   if (err) std::cout << "   error:: " << name << " " << shader_p->name << " draw() post uniforms "
                      << err << std::endl;

   if (vao == 99999999)
      std::cout << "You forget to setup this mesh " << name << " " << shader_p->name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: render_atom_label A3 " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);

   err = glGetError();
   if (err) std::cout << "   error draw() " << name << " pre-draw " << err << std::endl;

   // scale of 20 make letter separation good
   GLfloat scale = 100.1;

   // ------------------------------- text texture code here -----------------------

   // consider refactoring when working

   // display_info_t di;

   float x = 0;
   float y = 0;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   shader.Use();
   err = glGetError(); if (err) std::cout << "error:: render_atom_label A0 " << err << std::endl;

#if THIS_IS_HMT
   std::map<GLchar, FT_character> &ft_characters = display_info_t::ft_characters;
#else
   std::map<GLchar, FT_character> &ft_characters = graphics_info_t::ft_characters;
#endif

   std::string::const_iterator c; // call this it_c
   for (c = atom_label.begin(); c != atom_label.end(); c++) {
      std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*c);
      if (it == ft_characters.end()) {
         std::cout << "Failed to lookup glyph for " << *c << std::endl;
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
   err = glGetError(); if (err) std::cout << "   error:: " << shader_p->name << " draw() post mvp uniform "
                                          << err << std::endl;
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

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);

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
