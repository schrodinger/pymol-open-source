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
#ifndef _H_CoordSet
#define _H_CoordSet

#include<Python.h>
#include"Rep.h"
#include"Symmetry.h"
#include"Word.h"

typedef struct CoordSet {
  void (*fUpdate)(struct CoordSet *I);
  void (*fRender)(struct CoordSet *I,CRay *ray,Pickable **pick);
  void (*fFree)(struct CoordSet *I);
  void (*fEnumIndices)(struct CoordSet *I);
  void (*fAppendIndices)(struct CoordSet *I,int existingAtoms);
  void (*fExtendIndices)(struct CoordSet *I,int nAtom);
  void (*fInvalidateRep)(struct CoordSet *I,int type,int level);
  struct ObjectMolecule *Obj;
  float *Coord;
  int *Color;
  int *IdxToAtm;
  int *AtmToIdx;
  int NIndex,NAtIndex;
  Rep *Rep[cRepCnt]; /* an array of pointers to representations */
  int Active[cRepCnt]; /* active flags */
  int NRep;
  int NTmpBond; /* optional + temporary (for coord set transfers) */
  int *TmpBond; /* actual bond info is stored in ObjectMolecule */
  CSymmetry *TmpSymmetry;
  WordType Name;
} CoordSet;

#include"ObjectMolecule.h"

CoordSet *CoordSetNew(void);
void CoordSetAtomToPDBStrVLA(char **charVLA,int *c,AtomInfoType *ai,float *v,int cnt);
void CoordSetAtomToTERStrVLA(char **charVLA,int *c,AtomInfoType *ai,int cnt);
CoordSet *CoordSetCopy(CoordSet *cs);

void CoordSetTransform44f(CoordSet *I,float *mat);
void CoordSetRealToFrac(CoordSet *I,CCrystal *cryst);
void CoordSetFracToReal(CoordSet *I,CCrystal *cryst);
void CoordSetGetAverage(CoordSet *I,float *v0);
PyObject *CoordSetAtomToChemPyAtom(AtomInfoType *ai,float *v);

#endif

