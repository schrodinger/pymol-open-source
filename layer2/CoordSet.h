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

#include"Rep.h"

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
  Rep **Rep; /* an array of pointers to representations */
  int NRep;
  int NTmpBond; /* optional + temporary (for coord set transfers) */
  int *TmpBond; /* actual bond info is stored in ObjectMolecule */
} CoordSet;

#include"ObjectMolecule.h"

CoordSet *CoordSetNew(void);
void CoordSetAtomToPDBStrVLA(char **charVLA,unsigned int *c,AtomInfoType *ai,float *v,int cnt);

#endif

