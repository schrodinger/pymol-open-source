
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
  // methods (not fully refactored yet)
  void fFree();

  // methods
  void update(int state);
  void render(RenderInfo * info);
  void enumIndices();
  void appendIndices(int offset);
  int extendIndices(int nAtom);
  void invalidateRep(int type, int level);
  int atmToIdx(int atm) const;

  // read/write pointer to coordinate
  float * coordPtr(int idx) {
    return Coord + idx * 3;
  }

  // read pointer to coordinate
  const float * coordPtr(int idx) const {
    return Coord + idx * 3;
  }

  AtomInfoType * getAtomInfo(int idx) {
    return Obj->AtomInfo + IdxToAtm[idx];
  }

  // true if any atom in this coord set has any of the reps in "bitmask" shown
  bool hasRep(int bitmask) {
    if (Obj->RepVisCache & bitmask)
      for (int idx = 0; idx < NIndex; idx++)
        if (getAtomInfo(idx)->visRep & bitmask)
          return true;
    return false;
  }

  CObjectState State;
  ObjectMolecule *Obj;
  float *Coord;
  int *IdxToAtm;
  int *AtmToIdx;
  int NIndex, NAtIndex, prevNIndex, prevNAtIndex;
  ::Rep *Rep[cRepCnt];            /* an array of pointers to representations */
  int Active[cRepCnt];          /* active flags */
  int NTmpBond;                 /* optional, temporary (for coord set transfers) */
  BondType *TmpBond;            /* actual bond info is stored in ObjectMolecule */
  int NTmpLinkBond;             /* optional, temporary storage of linkage  info. */
  BondType *TmpLinkBond;        /* first atom is in obj, second is in cset */
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
  int tmp_index;                /* for saving */

  LabPosType *LabPos;

  /* not saved in state */

  RefPosType *RefPos;

  /* idea:  
     int start_atix, stop_atix <-- for discrete objects, we need
     something like this that would enable pymol to skip atoms not in the
     discrete state...question is: are these atoms sorted together right
     now or not? probably not, and if not then we need to change sorting
     for discrete objects to be state-dependent, but this could screw up
     byres/bychain actions which assume such atoms to be adjancent...
   */

  CGO *SculptCGO, *SculptShaderCGO;
  MapType *Coord2Idx;
  float Coord2IdxReq, Coord2IdxDiv;

  /* temporary / optimization */

  int objMolOpInvalidated;
#ifndef NO_MMLIBS
  bool validMMStereo;
  bool validTextType;
#endif

#ifdef _PYMOL_IP_EXTRAS
#endif
} CoordSet;

typedef void (*fUpdateFn) (CoordSet *, int);

#define cCSet_NoPeriodicity 0
#define cCSet_Orthogonal 1
#define cCSet_Octahedral 2

int BondInOrder(BondType * a, int b1, int b2);
int BondCompare(BondType * a, BondType * b);

PyObject *CoordSetAsNumPyArray(CoordSet * cs, short copy);
PyObject *CoordSetAsPyList(CoordSet * I);
int CoordSetFromPyList(PyMOLGlobals * G, PyObject * list, CoordSet ** cs);

CoordSet *CoordSetNew(PyMOLGlobals * G);
void CoordSetAtomToPDBStrVLA(PyMOLGlobals * G, char **charVLA, int *c,
                             const AtomInfoType * ai,
                             const float *v, int cnt,
                             const PDBInfoRec * pdb_info,
                             const double *matrix);
void CoordSetAtomToTERStrVLA(PyMOLGlobals * G, char **charVLA, int *c, AtomInfoType * ai,
                             int cnt);
CoordSet *CoordSetCopy(const CoordSet * cs);

void CoordSetTransform44f(CoordSet * I, const float *mat);
void CoordSetTransform33f(CoordSet * I, const float *mat);
void CoordSetRealToFrac(CoordSet * I, const CCrystal * cryst);
void CoordSetFracToReal(CoordSet * I, const CCrystal * cryst);

bool CoordSetInsureOrthogonal(PyMOLGlobals * G,
    CoordSet * cset,
    const float * sca,
    const CCrystal *cryst=NULL,
    bool quiet=true);

void CoordSetGetAverage(CoordSet * I, float *v0);
PyObject *CoordSetAtomToChemPyAtom(PyMOLGlobals * G, AtomInfoType * ai, const float *v,
                                   const float *ref, int index, const double *matrix);
int CoordSetGetAtomVertex(CoordSet * I, int at, float *v);
int CoordSetGetAtomTxfVertex(CoordSet * I, int at, float *v);
int CoordSetSetAtomVertex(CoordSet * I, int at, const float *v);
int CoordSetMoveAtom(CoordSet * I, int at, const float *v, int mode);
int CoordSetMoveAtomLabel(CoordSet * I, int at, const float *v, int mode);

int CoordSetTransformAtomTTTf(CoordSet * I, int at, const float *TTT);
int CoordSetTransformAtomR44f(CoordSet * I, int at, const float *matrix);

int CoordSetValidateRefPos(CoordSet * I);

void CoordSetPurge(CoordSet * I);
void CoordSetAdjustAtmIdx(CoordSet * I, int *lookup, int nAtom);
int CoordSetMerge(ObjectMolecule *OM, CoordSet * I, CoordSet * cs);        /* must be non-overlapping */
void CoordSetRecordTxfApplied(CoordSet * I, const float *TTT, int homogenous);
void CoordSetUpdateCoord2IdxMap(CoordSet * I, float cutoff);

typedef struct _CCoordSetUpdateThreadInfo CCoordSetUpdateThreadInfo;

void CoordSetUpdateThread(CCoordSetUpdateThreadInfo * T);

void LabPosTypeCopy(const LabPosType * src, LabPosType * dst);
void RefPosTypeCopy(const RefPosType * src, RefPosType * dst);

// object-state level setting
template <typename V> void SettingSet(int index, V value, CoordSet *cs) {
  SettingSet(cs->State.G, &cs->Setting, index, value);
}

#endif
