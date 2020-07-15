
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
#ifndef _H_Crystal
#define _H_Crystal

#include"Vector.h"
#include"Base.h"
#include"os_python.h"

struct CCrystal {
  PyMOLGlobals *G;
  float Dim[3] = {1.0, 1.0, 1.0};
  float Angle[3] = {90.0, 90.0, 90.0}; /* stored in degrees for convenience */
  float RealToFrac[9];
  float FracToReal[9];
  float UnitCellVolume = 1.0;
  float Norm[3];
  float RecipDim[3];
  CCrystal(PyMOLGlobals* GParam);
};

void CrystalUpdate(CCrystal* I);
void CrystalDump(const CCrystal* I);
CGO* CrystalGetUnitCellCGO(const CCrystal* I);
int CrystalFromPyList(CCrystal* I, PyObject* list);
PyObject* CrystalAsPyList(const CCrystal* I);

#endif
