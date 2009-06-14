
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
#ifndef _H_Menu
#define _H_Menu

#include"os_python.h"
#include"PyMOLGlobals.h"

void MenuActivate(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                  int passive, char *name, char *sele);
void MenuActivate0Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                      int passive, char *name);
void MenuActivate1Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                      int passive, char *name, char *arg1);
void MenuActivate2Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                      int passive, char *name, char *sele1, char *sele2);
void MenuActivate3fv(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                     char *name, float *xyz);

#endif
