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

#include "Base.h"
#include "pymol/math_defines.h"
#include "PyMOLObject.h"
#include "Rep.h"
#include "Result.h"
#include "vla.h"

#include <array>
#include <forward_list>

struct ObjectDist;

struct CMeasureInfo {
  /* AtomInfoType.unique_id */
  int id[4];
  /* offset into this distance set's Coord list */
  int offset;
  /* save object state */
  int state[4];
  /* distance, angle, or dihedral */
  int measureType;
};

struct DistSet : CObjectState {
  DistSet(PyMOLGlobals *);

  // methods
  void update(int state);
  void render(RenderInfo *);
  void invalidateRep(cRep_t type, cRepInv_t level);

  ObjectDist *Obj = nullptr;

  pymol::vla<float> Coord;
  int NIndex = 0;

  /* an array of pointers to representations */
  std::array<pymol::cache_ptr<::Rep>, cRepCnt> Rep;
  static int getNRep() { return cRepCnt; }

  std::vector<pymol::Vec3> LabCoord;
  std::vector<LabPosType> LabPos;

  pymol::vla<float> AngleCoord;
  int NAngleIndex = 0;

  pymol::vla<float> DihedralCoord;
  int NDihedralIndex = 0;

  std::forward_list<CMeasureInfo> MeasureInfo;

  pymol::Result<pymol::Vec3> getLabelOffset(int atm) const;
  pymol::Result<> setLabelOffset(int atm, const float* pos);
};

#define DistSetNew(G) (new DistSet(G))
PyObject *DistSetAsPyList(DistSet * I);
DistSet* DistSetFromPyList(PyMOLGlobals * G, PyObject * list);
int DistSetGetExtent(DistSet * I, float *mn, float *mx);
int DistSetMoveLabel(DistSet * I, int at, float *v, int mode);
int DistSetGetLabelVertex(DistSet * I, int at, float *v);
/* -- JV */
int DistSetMoveWithObject(DistSet* I, struct ObjectMolecule * O);
/* -- JV end */

#endif
