#ifndef _H_os_gl_pre
#define _H_os_gl_pre

#include"os_predef.h"
#include"os_proprietary.h"

#if defined(_PYMOL_IOS) || defined(ANDROID)
#define PURE_OPENGL_ES
#define _PYMOL_PURE_OPENGL_ES
#define _PYMOL_CGO_DRAWARRAYS
#define _PYMOL_GL_DRAWARRAYS
#define _PYMOL_CGO_DRAWBUFFERS
#if !defined(OPENGL_ES_1) && !defined(OPENGL_ES_2)
#define OPENGL_ES_1
#endif
#endif

#if !defined(OPENGL_ES_1) && !defined(OPENGL_ES_2)
#define OPENGL_ES_1 1
#endif

#if defined(_PYMOL_PURE_OPENGL_ES) && defined(OPENGL_ES_2)
#define PURE_OPENGL_ES_2 1
#endif

#if defined(_PYMOL_PURE_OPENGL_ES) && !defined(_PYMOL_GL_DRAWARRAYS)
#define _PYMOL_GL_DRAWARRAYS
#endif

#endif

