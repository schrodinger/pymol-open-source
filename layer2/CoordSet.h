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

#include"os_python.h"
#include"Rep.h"
#include"Symmetry.h"
#include"Word.h"
#include"Setting.h"
#include"ObjectMolecule.h"

typedef struct CoordSet {
  void (*fUpdate)(struct CoordSet *I);
  void (*fRender)(struct CoordSet *I,CRay *ray,Pickable **pick,int pass);
  void (*fFree)(struct CoordSet *I);
  void (*fEnumIndices)(struct CoordSet *I);
  void (*fAppendIndices)(struct CoordSet *I,int existingAtoms);
  void (*fExtendIndices)(struct CoordSet *I,int nAtom);
  void (*fInvalidateRep)(struct CoordSet *I,int type,int level);
  ObjectMolecule *Obj;
  float *Coord;
  int *Color;
  int *IdxToAtm;
  int *AtmToIdx;
  int NIndex,NAtIndex;
  Rep *Rep[cRepCnt]; /* an array of pointers to representations */
  int Active[cRepCnt]; /* active flags */
  int NRep;
  int NTmpBond; /* optional, temporary (for coord set transfers) */
  BondType *TmpBond; /* actual bond info is stored in ObjectMolecule */
  int NTmpLinkBond; /* optional, temporary storage of linkage  info. */
  BondType *TmpLinkBond; /* first atom is in obj, second is in cset */
  CSymmetry *Symmetry;
  WordType Name;
  float *Spheroid;
  float *SpheroidNormal;
  int NSpheroid;
  int SpheroidSphereSize;
  CSetting *Setting;
} CoordSet;


CoordSet *CoordSetNew(void);
void CoordSetAtomToPDBStrVLA(char **charVLA,int *c,AtomInfoType *ai,float *v,int cnt);
void CoordSetAtomToTERStrVLA(char **charVLA,int *c,AtomInfoType *ai,int cnt);
CoordSet *CoordSetCopy(CoordSet *cs);

void CoordSetTransform44f(CoordSet *I,float *mat);
void CoordSetRealToFrac(CoordSet *I,CCrystal *cryst);
void CoordSetFracToReal(CoordSet *I,CCrystal *cryst);
void CoordSetGetAverage(CoordSet *I,float *v0);
PyObject *CoordSetAtomToChemPyAtom(AtomInfoType *ai,float *v,int index);
int CoordSetGetAtomVertex(CoordSet *I,int at,float *v);
int CoordSetSetAtomVertex(CoordSet *I,int at,float *v);
int CoordSetMoveAtom(CoordSet *I,int at,float *v,int mode);
int CoordSetTransformAtom(CoordSet *I,int at,float *TTT);
void CoordSetPurge(CoordSet *I);
void CoordSetAdjustAtmIdx(CoordSet *I,int *lookup,int nAtom);
void CoordSetMerge(CoordSet *I,CoordSet *cs); /* must be non-overlapping */

#endif


