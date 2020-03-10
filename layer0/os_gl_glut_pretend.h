/**
 * GLUT emulation mode
 */

#pragma once

#ifndef _PYMOL_PRETEND_GLUT
#error "do not include directly, include os_gl.h"
#endif

#define P_GLUT_DOWN             0
#define P_GLUT_UP               1

#define P_GLUT_KEY_LEFT         100
#define P_GLUT_KEY_UP           101
#define P_GLUT_KEY_RIGHT        102
#define P_GLUT_KEY_DOWN         103
#define P_GLUT_KEY_PAGE_UP      104
#define P_GLUT_KEY_PAGE_DOWN    105
#define P_GLUT_KEY_HOME         106
#define P_GLUT_KEY_END          107
#define P_GLUT_KEY_INSERT       108

#define P_GLUT_LEFT_BUTTON      0
#define P_GLUT_MIDDLE_BUTTON    1
#define P_GLUT_RIGHT_BUTTON     2
