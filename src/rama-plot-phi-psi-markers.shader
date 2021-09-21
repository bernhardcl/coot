
#shader vertex

#version 330 core

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 texCoord;
layout (location = 2) in vec2 positions; // instanced (e.g. phi, psi positions)

out vec2 texCoord_transfer;

uniform vec2 position;
uniform vec2 scales;

uniform vec2 window_resize_position_correction;
uniform vec2 window_resize_scales_correction;

// This shader is for textures, uses HUDTextureMesh with the
// draw_instances() draw method.
//
// It is used for Rama-plot points

void main() {

   vec2 scaled_vertices = (vertex + positions) * scales;
   vec2 p1 = scaled_vertices + position;
   vec2 p2 = p1 * window_resize_scales_correction;
   vec2 p3 = p2 + window_resize_position_correction;
   // vec2 p3 = p1;
   gl_Position = vec4(p3, -1.0, 1.0);
   texCoord_transfer = texCoord;
}


#shader fragment

#version 330 core

uniform sampler2D text; // change this name to "image_texture"

in vec2 texCoord_transfer;

out vec4 outputColor;

void main() {

   vec4 sampled = texture(text, texCoord_transfer);
   outputColor = sampled;
   if (outputColor.a < 0.6)
      discard;
}
