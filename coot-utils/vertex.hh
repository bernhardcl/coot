#ifndef VERTEX_HH
#define VERTEX_HH

#include <glm/glm.hpp>

// for standard objects at the origin - typically used in instancing
namespace coot {

   namespace api {

      //! a vertex with (only) postion and normal. Useful for instancing, perhaps
      class vn_vertex {
      public:
         glm::vec3 pos;
         glm::vec3 normal; // normalized on input
         vn_vertex(const glm::vec3 &pos_in,
                   const glm::vec3 &norm_in) :
            pos(pos_in), normal(norm_in) {}
         vn_vertex() {}
      };

      //! a vertex with postion, normal and colour
      class vnc_vertex {
      public:
         glm::vec3 pos;
         glm::vec3 normal; // normalized on input
         glm::vec4 color;  // or colour?
         vnc_vertex(const glm::vec3 &pos_in,
                    const glm::vec3 &norm_in,
                    const glm::vec4 &col_in) : pos(pos_in), normal(norm_in), color(col_in) {}
         explicit vnc_vertex(const vn_vertex &vn) :
            pos(vn.pos), normal(vn.normal), color(glm::vec4(0.5, 0.5, 0.5, 1.0)) {}
         vnc_vertex(const vn_vertex &vn, const glm::vec4 &c) : pos(vn.pos), normal(vn.normal), color(c) {}
         vnc_vertex() {}
      };

      //! 20230108-PE copied from generic-vertex.hh (then edited).
      //!
      //! vertex with rotation and translation (e.g. for oriented bonds)
      class vertex_with_rotation_translation {
      public:
         glm::mat3 model_rotation_matrix; // orientation
         glm::vec3 model_translation; // the coordinates of the first atom of the bond
         glm::vec3 pos;
         glm::vec3 normal; // normalized when set
         glm::vec4 colour;
         //! constructor
         vertex_with_rotation_translation(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &c) : pos(p), normal(n), colour(c) {}
         //! constructor
         vertex_with_rotation_translation(const vnc_vertex &v, const glm::vec3 &atom_position, float scale) :
            model_rotation_matrix(glm::mat3(1.0f)), model_translation(atom_position),
            pos(v.pos * scale), normal(v.normal), colour(v.color) {}
         //! constructor
         vertex_with_rotation_translation() {}
      };
   }
}


#endif // VERTEX_HH
