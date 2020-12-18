#include <stdlib.h>
#ifdef _WEBGL
#endif
#include <iostream>
#include "GraphicsUtil.h"
// -----------------------------------------------------------------------------
// UTIL
  // Prints a backtrace during runtime of the last ^ stack frames
  void print_trace();
  void print_trace() {
#ifdef _WEBGL
#else
    printf("Use debugger with `b %s` to get a backtrace\n", __func__);
#endif
  }

bool glCheckOkay() {
  int err = 0;
  if ((err = glGetError()) != 0) {
    printf("GL_ERROR : 0x%04x\n", err);
#if 0
    print_trace();
#endif
    return false;
  }
  return true;
}

/**
 * GL debugging callback - enable with "pymol --gldebug"
 *
 * glDebugMessageCallback(gl_debug_proc, NULL);
 * glEnable(GL_DEBUG_OUTPUT);
 */
void GLAPIENTRY gl_debug_proc(
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
void GLAPIENTRY gl_debug_proc(
    GLenum source,
    GLenum type,
    GLuint id,
    GLenum severity,
    GLsizei length,
    const GLchar *msg,
    void* userParams)
{
    gl_debug_proc(source, type, id, severity, length, msg,
                  const_cast<const void*>(userParams));
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
