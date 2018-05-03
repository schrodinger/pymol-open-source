#include <stdlib.h>
#if !defined(_WIN32) && !defined(_WEBGL)
#include <execinfo.h>
#endif
#ifdef _WEBGL
#endif
#include <iostream>
#include "GraphicsUtil.h"
// -----------------------------------------------------------------------------
// UTIL
namespace {
  const int stack_frames = 12;

  // Prints a backtrace during runtime of the last ^ stack frames
  void print_trace() {
#if !defined(_WIN32) && !defined(_WEBGL)
    void *array[stack_frames];
    size_t size;
    char **strings;

    size = backtrace (array, stack_frames);
    strings = backtrace_symbols (array, size);

    printf ("Obtained %zd stack frames.\n", size);

    for (size_t i = 1; i < size; i++)
      printf ("%s\n", strings[i]);

    free (strings);
#endif
#ifdef _WEBGL
#endif
  }
};

bool glCheckOkay() {
  int err = 0;
  if ((err = glGetError()) != 0) {
    printf("GL_ERROR : %d\n", err);
    print_trace();
    return false;
  }
  return true;
}

/*
 * GL debugging callback - enable with "pymol --gldebug"
 *
 * glDebugMessageCallback(gl_debug_proc, NULL);
 * glEnable(GL_DEBUG_OUTPUT);
 */
void gl_debug_proc(
    GLenum source,
    GLenum type,
    GLuint id,
    GLenum severity,
    GLsizei length,
    const GLchar *msg,
    const void *)
{
#ifdef GL_DEBUG_TYPE_ERROR
  if (type == GL_DEBUG_TYPE_ERROR) {
    printf("glDebug: %s\n", msg);
    print_trace();
  }
#endif
}

// Returns the size in bytes of the opengl type
size_t gl_sizeof(GLenum type){
  size_t size = 1;
  switch (type) {
  case GL_UNSIGNED_BYTE:
  case GL_BYTE:
    size = 1;
    break;
  case GL_UNSIGNED_SHORT:
  case GL_SHORT:
    size = 2;
    break;
  case GL_UNSIGNED_INT:
  case GL_INT:
  case GL_FLOAT:
    size = 4;
    break;
  default :
    printf("Unsupported GL Type!");
    break;
  }
  return size;
}
