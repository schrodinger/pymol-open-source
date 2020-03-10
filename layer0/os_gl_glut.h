/**
 * GLUT (real or pretended)
 */

#pragma once

#ifdef _PYMOL_NO_MAIN
#define _PYMOL_NO_GLUT
#endif

#ifdef _PYMOL_NO_GLUT
#define _PYMOL_PRETEND_GLUT
#endif

#ifndef _PYMOL_PRETEND_GLUT
#include "os_gl_glut_real.h"
#else
#include "os_gl_glut_pretend.h"
#endif

#define P_GLUT_BUTTON_SCROLL_FORWARD  3
#define P_GLUT_BUTTON_SCROLL_BACKWARD 4
#define P_GLUT_SINGLE_LEFT      100
#define P_GLUT_SINGLE_MIDDLE    101
#define P_GLUT_SINGLE_RIGHT     102
#define P_GLUT_DOUBLE_LEFT      200
#define P_GLUT_DOUBLE_MIDDLE    201
#define P_GLUT_DOUBLE_RIGHT     202
