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

#include"Crystal.h"
#include"Word.h"

typedef struct { 
  CCrystal *Crystal;
  WordType PDBSpaceGroup;
  int PDBZValue;
  WordType SpaceGroup;
  int NSymMat;
  float *SymMatVLA;
  int NSymOp;
  WordType *SymOpVLA;
} CSymmetry;

void SymmetryAttemptGeneration(CSymmetry *I);
void SymmetryFree(CSymmetry *I);
CSymmetry *SymmetryNew(void);
void SymmetryUpdate(CSymmetry *I);
void SymmetryDump(CSymmetry *I);

#endif









