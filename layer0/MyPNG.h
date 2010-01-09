

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_MyPNG
#define _H_MyPNG

#include"PyMOLGlobals.h"

#define cMyPNG_FormatPNG 0
#define cMyPNG_FormatPPM 1

int MyPNGWrite(PyMOLGlobals * G, char *file_name, unsigned char *p,
               unsigned int width, unsigned int height, float dpi, int format, int quiet);
int MyPNGRead(char *file_name, unsigned char **p_ptr, unsigned int *width_ptr,
              unsigned int *height_ptr);

#endif
