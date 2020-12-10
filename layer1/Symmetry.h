
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
#include"vla.h"

struct CSymmetry {
  PyMOLGlobals *G;
  CCrystal Crystal;
  int PDBZValue = 0;

private:
  WordType SpaceGroup{};
  pymol::vla<float> SymMatVLA;

  bool updateSymMatVLA() const;

public:
  /// Space group name
  const char* spaceGroup() const { return SpaceGroup; }
  void setSpaceGroup(const char*);

  /// Get the number of symmetry matrices
  int getNSymMat() const;

  /// Get the i'th symmetry matrix (pointer to float[16]).
  /// Only valid after getNSymMat() has been called.
  const float * getSymMat(int i) const {
    return SymMatVLA + i * 16;
  }
  CSymmetry(PyMOLGlobals* G) : G(G), Crystal(G){};
};

#ifdef SYM_TO_MAT_LIST_IN_C
#endif
PyObject *SymmetryAsPyList(CSymmetry * I);
CSymmetry *SymmetryNewFromPyList(PyMOLGlobals * G, PyObject * list);
void SymmetrySpaceGroupRegister(PyMOLGlobals * G, const char* sg, const std::vector<std::string>& sym_op);

#endif
