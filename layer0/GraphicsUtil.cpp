#include "GraphicsUtil.h"

#include "PyMOLGlobals.h"
#include "Feedback.h"

#include <stdlib.h>
#ifdef _WEBGL
#include <emscripten.h>
#endif

#include <iostream>
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

bool CheckGLErrorOK(PyMOLGlobals* G, std::string_view errString)
{
  GLenum err;
  if ((err = glGetError()) != 0) {
#ifdef _WEBGL
    print_trace();
#else
    if (G) {
      PRINTFB(G, FB_CGO, FB_Errors)
      "GL_Error: 0x%04x @ %s\n", err, errString.data() ENDFB(G);
    } else {
      printf("GL_ERROR : 0x%04x %s\n", err, errString.data());
    }
    std::terminate();
#endif
  }
  return err == 0;
}

/**
 * GL debugging callback - enable with "pymol --gldebug"
 *
 * glDebugMessageCallback(gl_debug_proc, nullptr);
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
