#pragma once

#include "os_gl.h"

#ifdef PURE_OPENGL_ES_2
# define GLEW_EXT_gpu_shader4 false
# define GLEW_EXT_geometry_shader4 false
# ifdef _WEBGL
#  include <emscripten/val.h>
#  define GLEW_EXT_draw_buffers2 !emscripten::val::module_property("ONEBUFFER").as<bool>()
# else
#  define GLEW_EXT_draw_buffers2 false
# endif
#endif

#define TM3_IS_ONEBUF !GLEW_EXT_draw_buffers2

