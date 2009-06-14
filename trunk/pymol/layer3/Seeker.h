
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
#ifndef _H_Seeker
#define _H_Seeker

#include"Ortho.h"
#include"ObjectMolecule.h"

int SeekerInit(PyMOLGlobals * G);
void SeekerFree(PyMOLGlobals * G);
void SeekerUpdate(PyMOLGlobals * G);
char SeekerGetAbbr(PyMOLGlobals * G, char *abbr, char water, char unknown);

PyObject *SeekerGetRawAlignment(PyMOLGlobals * G, int align_sele, int active_only);

#endif
