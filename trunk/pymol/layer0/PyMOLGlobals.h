/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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
#ifndef _H_PyMOLGlobals
#define _H_PyMOLGlobals

typedef struct _CIsosurf CIsosurf;
typedef struct _CTetsurf CTetsurf;
typedef struct _CSphere CSphere;
typedef struct _CUtil CUtil;

typedef struct PyMOLGlobals {

  CIsosurf *Isosurf;
  CTetsurf *Tetsurf;
  CSphere  *Sphere;
  CUtil    *Util;
  
} PyMOLGlobals;

extern PyMOLGlobals *TempPyMOLGlobals;

#endif

