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
-* Mark J. Kilgard
-* 
-*
Z* -------------------------------------------------------------------
*/

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Grap.h"

void GrapMoveTo(int x,int y)
{
  glRasterPos4d((double)(x),(double)(y),0.0,1.0);
}

void GrapDrawStr(char *c,int x,int y)
{
  glRasterPos4d((double)(x),(double)(y),0.0,1.0);
  while(*c) 
    p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c++));
}


void GrapDrawSubStrSafe(char *c,int x,int y,int start,int n)
{
  glRasterPos4d((double)(x),(double)(y),0.0,1.0);
  while(*c) {
    if((start--)<1) {
      p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c));
      n--;
      if(n<=0) break;
    }
    c++;
  }
}

void GrapDrawSubStrFast(char *c,int x,int y,int start,int n)
{
  c+=start;
  glRasterPos4d((double)(x),(double)(y),0.0,1.0);
  if(n)
    while(*c) {
      p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c));
      n--;
      if(n<=0) break;
      c++;
    }
}

void GrapContStr(char *c)
{
  while(*c) 
    p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c++));
}

int GrapMeasureStr(char *c)
{
  return strlen(c)*8;
}
