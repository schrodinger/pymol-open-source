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

#define COORD_SET_HAS_ANISOU 0x01

typedef struct CoordSet {
  void (*fUpdate)(struct CoordSet *I,int state);
  void (*fRender)(struct CoordSet *I,RenderInfo *info);
  void (*fFree)(struct CoordSet *I);
  void (*fEnumIndices)(struct CoordSet *I);
  void (*fAppendIndices)(struct CoordSet *I,int existingAtoms);
  void (*fExtendIndices)(struct CoordSet *I,int nAtom);
  void (*fInvalidateRep)(struct CoordSet *I,int type,int level);
  CObjectState State;
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
  /* for periodic MD boxes -- may be merge into symmetry lattice later... */
  CCrystal *PeriodicBox;
  int PeriodicBoxType;
  int tmp_index; /* for saving */

  int NMatrix; /* number of matrices for this coordinate set */
  double *MatrixVLA; /* end-to-end array of 16x16 matrices */
  LabPosType *LabPos;

  /* not saved in state */

  /* idea:  
  int start_atix, stop_atix <-- for discrete objects, we need
  something like this that would enable pymol to skip atoms not in the
  discrete state...question is: are these atoms sorted together right
  now or not? probably not, and if not then we need to change sorting
  for discrete objects to be state-dependent, but this could screw up
  byres/bychain actions which assume such atoms to be adjancent...
  */

  CGO *SculptCGO;
  MapType *Coord2Idx;
  float Coord2IdxReq,Coord2IdxDiv;

  /* temporary / optimization */

  int objMolOpInvalidated;
} CoordSet;

typedef void (*fUpdateFn)(CoordSet *,int);

#define cCSet_NoPeriodicity 0
#define cCSet_Orthogonal 1
#define cCSet_Octahedral 2

int BondInOrder(BondType *a,int b1,int b2);
int BondCompare(BondType *a,BondType *b);

PyObject *CoordSetAsPyList(CoordSet *I);
int CoordSetFromPyList(PyMOLGlobals *G,PyObject *list,CoordSet **cs);

CoordSet *CoordSetNew(PyMOLGlobals *G);
void CoordSetAtomToPDBStrVLA(PyMOLGlobals *G,char **charVLA,int *c,AtomInfoType *ai,
                             float *v,int cnt,PDBInfoRec *pdb_info);
void CoordSetAtomToTERStrVLA(PyMOLGlobals *G,char **charVLA,int *c,AtomInfoType *ai,int cnt);
CoordSet *CoordSetCopy(CoordSet *cs);

void CoordSetTransform44f(CoordSet *I,float *mat);
void CoordSetTransform33f(CoordSet *I,float *mat);
void CoordSetRealToFrac(CoordSet *I,CCrystal *cryst);
void CoordSetFracToReal(CoordSet *I,CCrystal *cryst);
void CoordSetGetAverage(CoordSet *I,float *v0);
PyObject *CoordSetAtomToChemPyAtom(PyMOLGlobals *G,AtomInfoType *ai,float *v,int index);
int CoordSetGetAtomVertex(CoordSet *I,int at,float *v);
int CoordSetGetAtomTxfVertex(CoordSet *I,int at,float *v);
int CoordSetSetAtomVertex(CoordSet *I,int at,float *v);
int CoordSetMoveAtom(CoordSet *I,int at,float *v,int mode);
int CoordSetMoveAtomLabel(CoordSet *I,int at,float *v,int mode);

int CoordSetTransformAtomTTTf(CoordSet *I,int at,float *TTT);
int CoordSetTransformAtomR44f(CoordSet *I,int at,float *matrix);

void CoordSetPurge(CoordSet *I);
void CoordSetAdjustAtmIdx(CoordSet *I,int *lookup,int nAtom);
void CoordSetMerge(CoordSet *I,CoordSet *cs); /* must be non-overlapping */
void CoordSetRecordTxfApplied(CoordSet *I,float *TTT, int homogenous);
void CoordSetUpdateCoord2IdxMap(CoordSet *I, float cutoff);

typedef struct _CCoordSetUpdateThreadInfo CCoordSetUpdateThreadInfo;

void CoordSetUpdateThread(CCoordSetUpdateThreadInfo *T);

#endif


