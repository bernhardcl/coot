
#ifdef USE_PYTHON
#include <Python.h>
#endif

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>
#include "graphics-info.h"


bool
graphics_info_t::init_shaders() {

   std::vector<std::reference_wrapper<Shader> > shaders = {shader_for_maps,
                                                           shader_for_map_caps,
                                                           shader_for_models,
                                                           shader_for_outline_of_active_residue,
                                                           shader_for_model_as_meshes,
                                                           shader_for_symmetry_atoms_bond_lines,
                                                           shader_for_central_cube,
                                                           shader_for_origin_cube,
                                                           shader_for_hud_text,
                                                           shader_for_hud_geometry_bars,
                                                           shader_for_hud_geometry_labels,
                                                           shader_for_hud_geometry_tooltip_text,
                                                           shader_for_hud_buttons,
                                                           shader_for_hud_image_texture,
                                                           shader_for_atom_labels,
                                                           shader_for_moleculestotriangles,
                                                           shader_for_hud_lines,
                                                           shader_for_lines,
                                                           shader_for_lines_pulse,
                                                           shader_for_rama_balls,
                                                           shader_for_particles,
                                                           shader_for_instanced_objects,
                                                           shader_for_happy_face_residue_markers,
                                                           shader_for_rama_plot_phi_phis_markers,
                                                           shader_for_rama_plot_axes_and_ticks,
                                                           shader_for_ligand_view,
                                                           shader_for_screen,
                                                           shader_for_x_blur,
                                                           shader_for_y_blur,
                                                           shader_for_dof_blur_by_texture_combination,
                                                           shader_for_blur, // old, 2020 version
                                                           shader_for_texture_meshes,
                                                           camera_facing_quad_shader
   };

   bool status = true;  // success

   std::string p = coot::package_data_dir();
   std::string d = coot::util::append_dir_dir(p, "shaders");
   if (false)
      std::cout << "INFO:: shader default dir: " << d << std::endl;
   std::vector<std::reference_wrapper<Shader> >::iterator it;
   for (it=shaders.begin(); it!=shaders.end(); ++it)
      it->get().set_default_directory(d);

   shader_for_outline_of_active_residue.init("outline-of-active-residue.shader", Shader::Entity_t::MODEL);
   shader_for_maps.init("map.shader", Shader::Entity_t::MAP);
   shader_for_map_caps.init("draw-map-cap.shader", Shader::Entity_t::MAP);
   shader_for_models.init("model.shader", Shader::Entity_t::MODEL);
   shader_for_central_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   shader_for_origin_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   shader_for_hud_text.init("hud-text.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_geometry_bars.init("hud-bars.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_geometry_labels.init("hud-labels.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_image_texture.init("hud-image-texture.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_atom_labels.init("atom-label.shader", Shader::Entity_t::MODEL);
   shader_for_moleculestotriangles.init("moleculestotriangles.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_lines.init("lines.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_lines_pulse.init("lines-pulse.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_rama_balls.init("rama-balls.shader", Shader::Entity_t::MODEL);
   shader_for_particles.init("particles.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_instanced_objects.init("instanced-objects.shader", Shader::Entity_t::INSTANCED_DISPLAY_OBJECT);
   shader_for_hud_geometry_tooltip_text.init("hud-geometry-tooltip-text.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_happy_face_residue_markers.init("happy-face-residue-markers.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_ligand_view.init("ligand-view.shader", Shader::Entity_t::NONE);
   shader_for_model_as_meshes.init("model-as-mesh.shader", Shader::Entity_t::MODEL);
   shader_for_symmetry_atoms_bond_lines.init("symmetry-atoms-lines.shader", Shader::Entity_t::MAP);
   shader_for_hud_buttons.init("hud-bars.shader", Shader::Entity_t::HUD_TEXT); // ! needs a better name
   shader_for_rama_plot_axes_and_ticks.init("rama-plot-axes-and-ticks.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_rama_plot_phi_phis_markers.init("rama-plot-phi-psi-markers.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_lines.init("hud-lines.shader", Shader::Entity_t::MODEL);
   shader_for_texture_meshes.init("texture-meshes.shader", Shader::Entity_t::MAP);

   // testing image textures
   camera_facing_quad_shader.init("camera-facing-quad-shader-for-testing.shader", Shader::Entity_t::MODEL);

   // we use the above to make an image/texture in the framebuffer and use then
   // shader_for_screen to convert that framebuffer to the screen buffer.
   shader_for_screen.init("screen.shader", Shader::Entity_t::SCREEN);
   shader_for_x_blur.init("blur-x.shader", Shader::Entity_t::SCREEN);
   shader_for_y_blur.init("blur-y.shader", Shader::Entity_t::SCREEN);
   shader_for_dof_blur_by_texture_combination.init("depth-of-field.shader", Shader::Entity_t::SCREEN);
   shader_for_blur.init("blur.shader", Shader::Entity_t::SCREEN); // old

   for (it=shaders.begin(); it!=shaders.end(); ++it) {
      if (! it->get().get_success_status()) {
         std::cout << "ERROR:: shader \"" <<it->get().name << "\" failed" << std::endl;
         status = false;
      }
   }
   return status;
}



glm::vec4
graphics_info_t::unproject(float x, float y, float z) {

   // z is 1 and -1 for front and back (or vice verse).

   if (! glareas[0]) return glm::vec4(0,0,0,0);

   GtkAllocation allocation;
   gtk_widget_get_allocation(glareas[0], &allocation);
   int w = allocation.width;
   int h = allocation.height;

   float mouseX = x / (w * 0.5f) - 1.0f;
   float mouseY = (h - y) / (h * 0.5f) - 1.0f;

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos_f = glm::vec4(mouseX, mouseY, z, 1.0f); // maybe +1
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   if (false) {
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(mvp) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(vp_inv) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(screenPos_f) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(worldPos_f) << std::endl;
   }

   // to turn these into points in the world space, don't forget to divide by w.
   return worldPos_f;

}

// projected_coords are in clip space, so mouse position will have to be converted.
//
// static
glm::vec3
graphics_info_t::unproject_to_world_coordinates(const glm::vec3 &projected_coords) {

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos = glm::vec4(projected_coords, 1.0f);
   glm::vec4 c = vp_inv * screenPos;

   // std::cout << "debug:: unproject_to_world_coordinates() c " << glm::to_string(c) << std::endl;
   double oow = 1.0/c.w;
   return glm::vec3(c.x * oow, c.y * oow, c.z * oow);

}


int
graphics_info_t::blob_under_pointer_to_screen_centre() {

   // a quick slide would be better than a jump

   graphics_info_t g; // needed?
   int r = 0;
   if (use_graphics_interface_flag) {
      int imol_map = Imol_Refinement_Map();
      if (imol_map != -1) {
	 // OK we have a map to search.
	 // coot::Cartesian front = unproject(0.0);
	 // coot::Cartesian back  = unproject(1.0);
         // glm::vec4 glm_front = new_unproject(-0.3);
         // glm::vec4 glm_back  = new_unproject( 1.0);

         GtkAllocation allocation = graphics_info_t::get_glarea_allocation();
         int w = allocation.width;
         int h = allocation.height;

         glm::mat4 mvp = graphics_info_t::get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion
         glm::mat4 vp_inv = glm::inverse(mvp);

         float mouseX_2 = g.mouse_current_x  / (w * 0.5f) - 1.0f;
         float mouseY_2 = g.mouse_current_y  / (h * 0.5f) - 1.0f;
         // I revered the sign here - it does the right thing now.
         glm::vec4 screenPos_1 = glm::vec4(mouseX_2, -mouseY_2, -1.0f, 1.0f);
         glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2,  1.0f, 1.0f);
         glm::vec4 worldPos_1 = vp_inv * screenPos_1;
         glm::vec4 worldPos_2 = vp_inv * screenPos_2;

         double oowp_1 = 1.0/worldPos_1.w;
         double oowp_2 = 1.0/worldPos_2.w;
         coot::Cartesian front(worldPos_1.x * oowp_1, worldPos_1.y * oowp_1, worldPos_1.z * oowp_1);
	 coot::Cartesian  back(worldPos_2.x * oowp_2, worldPos_2.y * oowp_2, worldPos_2.z * oowp_2);

	 clipper::Coord_orth p1(front.x(), front.y(), front.z());
	 clipper::Coord_orth p2( back.x(),  back.y(),  back.z());
         if (false) {
            std::cout << "blob_under_pointer_to_screen_centre() " << glm::to_string(screenPos_1) << " "
                      << glm::to_string(screenPos_2) << std::endl;
            std::cout << "blob_under_pointer_to_screen_centre() " << front << " " << back << std::endl;
            // std::cout << "blob_under_pointer_to_screen_centre() " << p1.format() << " "
            // << p2.format() << std::endl;
         }
         coot::Cartesian rc = g.RotationCentre();

	 try {
	    clipper::Coord_orth blob =
	       molecules[imol_refinement_map].find_peak_along_line_favour_front(p1, p2);
	    coot::Cartesian cc(blob.x(), blob.y(), blob.z());
            // coot::Cartesian cc = front.mid_point(back);
            coot::Cartesian delta = rc - cc;
            // std::cout << "Delta: " << delta << std::endl;
	    g.setRotationCentre(cc);
	    for(int ii=0; ii<n_molecules(); ii++) {
	       molecules[ii].update_map(auto_recontour_map_flag);
	       molecules[ii].update_symmetry();
	    }
	    g.make_pointer_distance_objects();
	    graphics_draw();
	 }
	 catch (const std::runtime_error &mess) {
            std::cout << "debug:: given front " << front << " and back " << back << std::endl;
	    std::cout << mess.what() << std::endl;
	 }
      } else {
	 std::string s = "WARNING:: Refinement map not selected - no action";
	 std::cout << s << std::endl;
	 // add_status_bar_text(s.c_str());
	 info_dialog(s.c_str());
      }
   }
   return r;
}


void
graphics_info_t::set_clipping_front(float v) {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_near_perspective_limit = l * 0.99;
      if (v < screen_z_near_perspective_limit)
         if (v > 2.0)
            screen_z_near_perspective = v;
   } else {
      clipping_front = v;
   }
   graphics_draw();
}


void
graphics_info_t::set_clipping_back(float v) {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_far_perspective_limit = l * 1.01;
      if (v > screen_z_far_perspective_limit)
         if (v < 1000.0)
            screen_z_far_perspective = v;
   } else {
      clipping_back = v;
   }
   graphics_draw();
}


void
graphics_info_t::adjust_clipping(float d) {

   if (! graphics_info_t::perspective_projection_flag) {

      clipping_front = clipping_front * (1.0 + d);
      clipping_back  = clipping_back  * (1.0 + d);

   } else {

      // --- perspective ---

      double l = eye_position.z;
      double zf = screen_z_far_perspective;
      double zn = screen_z_near_perspective;

      // we should (and now do) concern ourselves with the distance to
      // the rotation centre so that the clipping planes are not
      // changed so that rotation centre is clipped.

      if (d < 0) {

         // close down (narrow)

         screen_z_near_perspective = l - (l-zn) * 0.97;
         screen_z_far_perspective  = l + (zf-l) * 0.95;

      } else {

         // expand

         screen_z_far_perspective  = l + (zf-l) * 1.05;
         screen_z_near_perspective = l - (l-zn) * 1.03;

      }

      float screen_z_near_perspective_limit = l * 0.99;
      float screen_z_far_perspective_limit  = l * 1.01;
      if (screen_z_near_perspective > screen_z_near_perspective_limit)
         screen_z_near_perspective = screen_z_near_perspective_limit;
      if (screen_z_far_perspective < screen_z_far_perspective_limit)
         screen_z_far_perspective = screen_z_far_perspective_limit;

      if (screen_z_near_perspective <    2.0) screen_z_near_perspective =    2.0;
      if (screen_z_far_perspective  > 1000.0) screen_z_far_perspective  = 1000.0;

      std::cout << "adjust_clipping(): debug l " << l
                << "    post-manip: " << screen_z_near_perspective << " "
                << screen_z_far_perspective << std::endl;
   }
}

void
graphics_info_t::set_view_quaternion(float i, float j, float k, float l) {

   // currently sets glm_quat (that's not a good name)
   // change it it view_quaternion

   // maybe the order will be wrong
   glm::quat q(i, j, k, l);
   glm_quat = q;
}


//static
void
graphics_info_t::update_view_quaternion(int area_width, int area_height) {

   graphics_info_t g;
   float tbs = g.get_trackball_size();

   glm::quat tb_quat =
      g.trackball_to_quaternion((2.0*g.GetMouseBeginX() - area_width)/area_width,
                                (area_height - 2.0*g.GetMouseBeginY())/area_height,
                                (2.0*g.mouse_current_x - area_width)/area_width,
                                (area_height - 2.0*g.mouse_current_y)/area_height,
                                tbs);

   tb_quat = glm::conjugate(tb_quat); // hooray, no more "backwards" mouse motion
   glm::quat product = tb_quat * glm_quat;
   glm_quat = glm::normalize(product);

}



#include "glarea_tick_function.hh"

// static
// extra_annotation is a default argument, default false;
void
graphics_info_t::setup_cylinder_clashes(const coot::atom_overlaps_dots_container_t &c,
                                        int imol, float tube_radius, bool extra_annotation) {

   // clashes - 20211008-PE this should be in a lambda for clarity
   //
   //             We can't do cylinders with this shader! So make a ball instead.
   // 20210910-PE Let's try to make another separate instancing mesh for clashes.
   //
   if (c.clashes.size() > 0) {
      graphics_info_t g;
      std::string clashes_name = "  Molecule " + coot::util::int_to_string(imol) + ":";
      clashes_name += " clashes insta-mesh";
      int clashes_obj_index = generic_object_index(clashes_name);
      if (clashes_obj_index == -1)
         clashes_obj_index = g.new_generic_object_number_for_molecule(clashes_name, imol); // make static?
      else
         g.generic_display_objects[clashes_obj_index].clear();

      float dimmer = 0.66;
      float z_scale = 0.37;
      float unstubby_cap_factor = 1.1/z_scale; // see below
      if (extra_annotation) {
         dimmer = 0.95;
         unstubby_cap_factor = 2.1;
         z_scale = 0.7;
      }
      coot::colour_holder clash_col = colour_values_from_colour_name("#ff59c9");
      clash_col.scale_intensity(dimmer);
      glm::vec4 clash_col_glm(clash_col.red, clash_col.green, clash_col.blue, 1.0);

      // instancing for capped cylinders
      meshed_generic_display_object &obj = g.generic_display_objects[clashes_obj_index];
      float line_radius = 0.062f;
      line_radius = tube_radius;
      const unsigned int n_slices = 16;
      std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(0,0,0), glm::vec3(0,0,1));
      obj.add_cylinder(pos_pair, clash_col, line_radius, n_slices, true, true,
                       meshed_generic_display_object::ROUNDED_CAP,
                       meshed_generic_display_object::ROUNDED_CAP, extra_annotation, unstubby_cap_factor); // does obj.mesh.import()
      // now I need to init the buffers of obj.mesh.
      Material material;
      if (extra_annotation) {
         material.shininess = 199.9;
         material.specular_strength = 0.9;
      }
      obj.mesh.setup(material); // calls setup_buffers()
      //
      // now accumulate the instancing matrices (the colours will stay the same)
      std::vector<glm::mat4> mats;
      for (unsigned int i=0; i<c.clashes.size(); i++) {
         std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                  glm::vec3(coord_orth_to_glm(c.clashes[i].second)));
         const glm::vec3 &start  = pos_pair.first;
         const glm::vec3 &finish = pos_pair.second;
         glm::vec3 b = finish - start;
         glm::vec3 normalized_bond_orientation(glm::normalize(b));
         glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0));
         glm::vec3 sc(1.1, 1.1, z_scale);
         glm::mat4 unit(1.0);
         glm::mat4 mt_1 = glm::translate(unit, start);
         glm::mat4 mt_2 = mt_1 * ori;
         glm::mat4 mt_3 = glm::scale(mt_2, sc);
         mats.push_back(mt_3);
      }

      auto set_display_generic_object_simple = [] (int object_number, short int istate) {
                                                  if (object_number >=0  && object_number < int(generic_display_objects.size())) {
                                                     generic_display_objects[object_number].mesh.set_draw_this_mesh(istate);
                                                  } else {
                                                     std::cout << "ERROR:: object_number in to_generic_object_add_point: "
                                                               << object_number << std::endl;
                                                  }
                                               };

      // mats.resize(1);
      std::vector<glm::vec4> cols(c.clashes.size(), clash_col_glm);
      unsigned int n_instances = mats.size();
      obj.mesh.setup_rtsc_instancing(nullptr, mats, cols, n_instances, material); // also does setup_buffers()
      obj.mesh.update_instancing_buffer_data(mats, cols); // is this needed?
      set_display_generic_object_simple(clashes_obj_index, 1);

      // add continuous updating
      if (! tick_function_is_active()) {
         tick_function_id = gtk_widget_add_tick_callback(g.glareas[0], glarea_tick_func, 0, 0); // turn off by turn off FPS monitor
      }
      g.do_tick_constant_draw = true;
   }
};


//static
bool
graphics_info_t::get_exta_annotation_state() {

   bool extra_annotation = false;
   time_t times = time(NULL);
   struct tm result;
   localtime_r(&times, &result);

   if (result.tm_mday == 1)
      if (result.tm_mon == 3)
         if (result.tm_sec%5==0)
            extra_annotation = true;
   if (result.tm_mon == 9)
      if (result.tm_mday > 15)
         if (result.tm_sec%5==0)
            extra_annotation = true;
   return extra_annotation;
}



// probably not the right place for this function
//
void
graphics_info_t::coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol) {

   unsigned int octasphere_subdivisions = 1; // make a member of graphics_info_t with an API
   octasphere_subdivisions = contact_dot_sphere_subdivisions; // 20211129-PE done

   if (true) {
      bool ignore_waters = true;
      // I don't like this part
      std::map<std::string, coot::colour_holder> colour_map;
      colour_map["blue"      ] = colour_values_from_colour_name("blue");
      colour_map["sky"       ] = colour_values_from_colour_name("sky");
      colour_map["sea"       ] = colour_values_from_colour_name("sea");
      colour_map["greentint" ] = colour_values_from_colour_name("greentint");
      colour_map["darkpurple"] = colour_values_from_colour_name("darkpurple");
      colour_map["green"     ] = colour_values_from_colour_name("green");
      colour_map["orange"    ] = colour_values_from_colour_name("orange");
      colour_map["orangered" ] = colour_values_from_colour_name("orangered");
      colour_map["yellow"    ] = colour_values_from_colour_name("yellow");
      colour_map["yellowtint"] = colour_values_from_colour_name("yellowtint");
      colour_map["red"       ] = colour_values_from_colour_name("red");
      colour_map["#55dd55"   ] = colour_values_from_colour_name("#55dd55");
      colour_map["hotpink"   ] = colour_values_from_colour_name("hotpink");
      colour_map["grey"      ] = colour_values_from_colour_name("grey");
      colour_map["magenta"   ] = colour_values_from_colour_name("magenta");
      colour_map["royalblue" ] = colour_values_from_colour_name("royalblue");

      auto colour_string_to_colour_holder = [colour_map] (const std::string &c) {
                                               return colour_map.find(c)->second;
                                            };

      coot::atom_overlaps_dots_container_t c;

#if 0 // why did I want to add contact dots for a (static) molecule when moving atoms were being displayed?
      // Weird.

      // get_moving_atoms_lock(__FUNCTION__);
      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            if (moving_atoms_asc->n_selected_atoms > 0) {
               coot::atom_overlaps_container_t overlaps(mol, graphics_info_t::Geom_p(), ignore_waters, 0.5, 0.25);
               c = overlaps.all_atom_contact_dots(contact_dots_density, true);
            }
         }
      }
      // release_moving_atoms_lock(__FUNCTION__);
#endif

      // more sensible
      coot::atom_overlaps_container_t overlaps(mol, graphics_info_t::Geom_p(), ignore_waters, 0.5, 0.25);
      bool do_vdw_surface = all_atom_contact_dots_do_vdw_surface; // static class variable
      c = overlaps.all_atom_contact_dots(contact_dots_density, do_vdw_surface);

      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
      std::string molecule_name_stub = "Contact Dots for Molecule ";
      molecule_name_stub += coot::util::int_to_string(imol);
      molecule_name_stub += ": ";

      float point_size = 0.05f;
      std::unordered_map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
      for (it=c.dots.begin(); it!=c.dots.end(); ++it) {
         float point_size_inner = point_size;
	 const std::string &type = it->first;
	 const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
         // std::cout << "dots type " << type << " count " << v.size() << std::endl;
         float specular_strength = 0.5; // default
         if (type == "vdw-surface") specular_strength= 0.1; // dull, reduces zoomed out speckles
         if (type == "vdw-surface") point_size_inner = 0.03;
         std::string mesh_name = molecule_name_stub + type;

         Instanced_Markup_Mesh &im = graphics_info_t::molecules[imol].find_or_make_new(mesh_name);
         im.clear();
         im.setup_octasphere(octasphere_subdivisions);
         im.setup_instancing_buffers(v.size());
         std::vector<Instanced_Markup_Mesh_attrib_t> balls;
         balls.resize(v.size());
         // I now do simple-minded colour caching here
         std::string previous_colour_string;
         glm::vec4 previous_colour(0,0,0,1);
         for (unsigned int i=0; i<v.size(); i++) {
            glm::vec3 position(v[i].pos.x(), v[i].pos.y(), v[i].pos.z());
            const std::string &colour_string = v[i].col;
            if (false) // debugging
               if (type != "vdw-surface")
                  std::cout << " type " << type << " " << i << " colour_string: " << colour_string << std::endl;
            if (colour_string == previous_colour_string) {
               Instanced_Markup_Mesh_attrib_t attribs(previous_colour, position, point_size_inner);
               attribs.specular_strength = specular_strength;
               balls[i] = attribs;
            } else {
               coot::colour_holder ch = colour_string_to_colour_holder(colour_string);
               glm::vec4 colour(ch.red, ch.green, ch.blue, 1.0);
               Instanced_Markup_Mesh_attrib_t attribs(colour, position, point_size_inner);
               attribs.specular_strength = specular_strength;
               balls[i] = attribs;
               previous_colour_string = colour_string;
               previous_colour = colour;
            }
         }
         im.update_instancing_buffers(balls);
      }

      float tube_size = point_size * 1.1f;
      setup_cylinder_clashes(c, imol, tube_size, get_exta_annotation_state());

   }

}
