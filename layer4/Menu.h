
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

struct PyMOLGlobals;
struct Block;

void MenuActivate(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                  int passive, const char *name, const char *sele);
void MenuActivate0Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                      int passive, const char *name);
Block *MenuActivate1Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                      int passive, const char *name, const char *arg1);
void MenuActivate2Arg(PyMOLGlobals * G, int x, int y, int last_x, int last_y,
                      int passive, const char *name, const char *sele1, const char *sele2);
void MenuActivate3fv(PyMOLGlobals * G, int x, int y, int last_x, int last_y, int passive,
                     const char *name, const float *xyz);

#endif
