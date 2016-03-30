#ifndef __TEXTURE_GLOBALS_HELP__
#define __TEXTURE_GLOBALS_HELP__

#define __GL_HELP__
#ifdef __APPLE_CC__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <utility>
#include "ppmutil.h"
enum {NEAREST, LINEAR, MIPMAP};
#if defined (_WIN32) && not defined (WINDOWS_MINGW)
#define TEXTURE_DLL __declspec(dllimport)
unsigned int __stdcall TEXTURE_DLL load_texture(const image_info &iinfo, int style);
std::pair<unsigned,unsigned> __stdcall TEXTURE_DLL GetCompatibleTextureSize(const unsigned width_in, const unsigned height_in);
image_info __stdcall TEXTURE_DLL ResizeWithEmptySpace(const image_info &iinfo, const unsigned width, const unsigned height);
void __stdcall TEXTURE_DLL set_texture_coord(GLfloat x, GLfloat y, int textured);
#else
unsigned int load_texture(const image_info &iinfo, int style);
std::pair<unsigned,unsigned> GetCompatibleTextureSize(const unsigned width_in, const unsigned height_in);
image_info ResizeWithEmptySpace(const image_info &iinfo, const unsigned width, const unsigned height);
void set_texture_coord(GLfloat x, GLfloat y, int textured);
#endif

#endif
