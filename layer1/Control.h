
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
#ifndef _H_Control
#define _H_Control

#include"Ortho.h"
#include "Result.h"

int ControlInit(PyMOLGlobals * G);
void ControlFree(PyMOLGlobals * G);
Block *ControlGetBlock(PyMOLGlobals * G);
int ControlIdling(PyMOLGlobals * G);
void ControlInterrupt(PyMOLGlobals * G);
pymol::Result<bool> ControlRock(PyMOLGlobals * G, int mode);
int ControlRocking(PyMOLGlobals * G);
int ControlSdofUpdate(PyMOLGlobals * G, float tx, float ty, float tz, float rx, float ry,
                      float rz);
int ControlSdofIterate(PyMOLGlobals * G);
int ControlSdofButton(PyMOLGlobals * G, int button);

#endif
