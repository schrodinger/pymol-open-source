
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
#ifndef _H_DistSet
#define _H_DistSet

#include"Base.h"
#include"Rep.h"
#include"Setting.h"
#include"ObjectDist.h"

typedef struct CMeasureInfo {
  /* AtomInfoType.unique_id */
  int id[4];
  /* offset into this distance set's Coord list */
  int offset;
  /* save object state */
  int state[4];
  /* distance, angle, or dihedral */
  int measureType;
  struct CMeasureInfo* next;
} CMeasureInfo;

typedef struct DistSet {
  // methods (not fully refactored yet)
  void fFree();

  // methods
  void update(int state);
  void render(RenderInfo *);
  void invalidateRep(int type, int level);

  CObjectState State;
  struct ObjectDist *Obj;
  float *Coord;
  int NIndex;
  ::Rep **Rep;                    /* an array of pointers to representations */
  int NRep;
  CSetting *Setting;
  /* extended for mobile distance labels */
  float *LabCoord;
  LabPosType *LabPos;
  int NLabel;
  /* extended for angles and torsions, with labels embedded with coordinates */
  float *AngleCoord;
  int NAngleIndex;
  float *DihedralCoord;
  int NDihedralIndex;
  /* -- JV */
  CMeasureInfo* MeasureInfo;
  /* -- JV end */
} DistSet;

DistSet *DistSetNew(PyMOLGlobals * G);
PyObject *DistSetAsPyList(DistSet * I);
int DistSetFromPyList(PyMOLGlobals * G, PyObject * list, DistSet ** cs);
int DistSetGetExtent(DistSet * I, float *mn, float *mx);
int DistSetMoveLabel(DistSet * I, int at, float *v, int mode);
int DistSetGetLabelVertex(DistSet * I, int at, float *v);
/* -- JV */
int DistSetMoveWithObject(DistSet* I, struct ObjectMolecule * O);
/* -- JV end */

#endif
