#ifndef _H_os_gl
#define _H_os_gl

#include"os_predef.h"
#include"os_proprietary.h"

// hardcode either true, or (x)
#define ALWAYS_IMMEDIATE_OR(x) true

#if 1
  #define _PYMOL_NO_AA_SHADERS
#endif

#if !defined(GL_GLEXT_PROTOTYPES) && !defined(_WIN32)
#define GL_GLEXT_PROTOTYPES
#endif

#ifndef GLEW_NO_GLU
#define GLEW_NO_GLU
#endif

#ifndef PURE_OPENGL_ES_2
#include <GL/glew.h>
#endif

#ifdef PURE_OPENGL_ES_2
#include "os_gl_es.h"
#elif defined(_PYMOL_OSX)
#import <OpenGL/gl.h>
#import <OpenGL/glext.h>
#else
#include <GL/gl.h>
#endif

#include "os_gl_glut.h"

void PyMOLReadPixels(GLint x,
                     GLint y,
                     GLsizei width,
                     GLsizei height, GLenum format, GLenum type, GLvoid * pixels);

void PyMOLDrawPixels(GLsizei width,
                     GLsizei height, GLenum format, GLenum type, const GLvoid * pixels);

int PyMOLCheckOpenGLErr(const char *pos);

#if defined(_PYMOL_IOS) && !defined(_WEBGL)
#define GL_C_INT_TYPE ushort
#define GL_C_INT_ENUM GL_UNSIGNED_SHORT
#define SceneGLClearColor(red,green,blue,alpha) if (!SceneGetBackgroundColorAlreadySet(G)) glClearColor(red,green,blue,alpha);
#else
#define GL_C_INT_TYPE uint
#define GL_C_INT_ENUM GL_UNSIGNED_INT
#define SceneGLClearColor(red,green,blue,alpha) glClearColor(red,green,blue,alpha);
#endif

#ifndef GLAPIENTRY
#define GLAPIENTRY
#endif

#define hasFrameBufferBinding() false

#ifndef PURE_OPENGL_ES_2
#define GL_DEBUG_PUSH(title) \
  GLEW_KHR_debug ? glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, 0, -1, title) : (void)0

#define GL_DEBUG_POP() \
  GLEW_KHR_debug ? glPopDebugGroup() : (void)0

#ifdef __cplusplus

class glDebugBlock {
public:
  explicit glDebugBlock(char const* title) {
    GL_DEBUG_PUSH(title);
  }
  ~glDebugBlock() {
    GL_DEBUG_POP();
  }
};

#define GL_DEBUG_FUN() \
  glDebugBlock glDebugBlockVariable(__FUNCTION__)

#endif /* __cplusplus */
#else
#define GL_DEBUG_FUN()
#endif

#endif
