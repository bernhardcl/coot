
#shader vertex

#version 330 core

// this shader uses instancing: The basic quad has the vertices
// and each block has colour and position/displacement (which places the
// bar on the screen)
//
// This shader is also used for HUD buttons

layout (location = 0) in vec2 vertex;
layout (location = 1) in float shade;
layout (location = 2) in vec4 colour;  // this and below are instanced
layout (location = 3) in vec2 position_offset;
layout (location = 4) in float scale_x;
layout (location = 5) in float scale_y;

out vec4 colour_transfer;

// culture clash! should the bars/buttons each know their offset position (and scales)
// (and they will stretch if the window is widened)
uniform vec2 position_offset_as_uniform;
uniform vec2 scales_as_uniform;

void main()
{
   // we want the tooltip to appear "over" this bar. Maybe
   // I could turn off depth test for that? Hmm. Not done at the moment
   gl_Position = vec4(scale_x * vertex.x + position_offset.x,
                      scale_y * vertex.y + position_offset.y,
                      -0.999, 1.0);
   // we adds colour adjust so that the buttons look smooth shaded
   float c = 0.24;
   float s3 = shade * shade * shade * shade * shade;
   vec4 colour_adjust = vec4(s3 * c, s3 * c, s3 * c, 0.1);
   colour_transfer = colour + colour_adjust;
}

#shader fragment

#version 330 core

in vec4 colour_transfer;
out vec4 colour;

void main()
{
   colour = colour_transfer;
}
