/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

#include"os_std.h"
#include"os_gl.h"

#include"Grap.h"

void GrapDrawStr(char *c,int x,int y)
{
  glRasterPos4d((double)(x),(double)(y),0.0,1.0);
  while(*c) 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
}

void GrapContStr(char *c)
{
  while(*c) 
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
}

