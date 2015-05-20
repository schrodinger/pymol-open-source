
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
#ifndef _H_Symmetry
#define _H_Symmetry

#include <vector>
#include <string>

#include"Crystal.h"
#include"Word.h"
#include"os_python.h"

struct CSymmetry {
  PyMOLGlobals *G;
  CCrystal *Crystal;
  int PDBZValue;
  WordType SpaceGroup;
  float *SymMatVLA;

  // get the number of symmetry matrices
  int getNSymMat() const;

  // get the i'th symmetry matrix (pointer to float[16])
  const float * getSymMat(int i) const {
    return SymMatVLA + i * 16;
  }
};

int SymmetryAttemptGeneration(CSymmetry * I, int quiet=false);
void SymmetryFree(CSymmetry * I);
void SymmetryClear(CSymmetry * I);
CSymmetry *SymmetryNew(PyMOLGlobals * G);
void SymmetryUpdate(CSymmetry * I);
void SymmetryDump(CSymmetry * I);
CSymmetry *SymmetryCopy(const CSymmetry * other);
PyObject *SymmetryAsPyList(CSymmetry * I);
CSymmetry *SymmetryNewFromPyList(PyMOLGlobals * G, PyObject * list);

void SymmetrySpaceGroupRegister(PyMOLGlobals * G, const char* sg, const std::vector<std::string>& sym_op);

#endif
