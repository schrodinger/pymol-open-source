/**
 * GLUT key and button constants (emulated)
 */

#pragma once

// There is no GLUT; these are the key/button codes PyMOL uses internally.
#define _PYMOL_PRETEND_GLUT
#include "os_gl_glut_pretend.h"


#define P_GLUT_BUTTON_SCROLL_FORWARD  3
#define P_GLUT_BUTTON_SCROLL_BACKWARD 4
#define P_GLUT_SINGLE_LEFT      100
#define P_GLUT_SINGLE_MIDDLE    101
#define P_GLUT_SINGLE_RIGHT     102
#define P_GLUT_DOUBLE_LEFT      200
#define P_GLUT_DOUBLE_MIDDLE    201
#define P_GLUT_DOUBLE_RIGHT     202
