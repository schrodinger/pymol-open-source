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
#ifndef _H_TestPyMOL
#define _H_TestPyMOL

#include"PyMOLGlobals.h"

typedef struct {
  PyMOLGlobals *G;
  int i;
} CTestPyMOL;

int TestPyMOLRun(PyMOLGlobals *G,CTestPyMOL *I,int group,int test);

#endif
