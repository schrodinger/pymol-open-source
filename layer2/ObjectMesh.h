
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
#ifndef _H_ObjectMesh
#define _H_ObjectMesh

#include"ObjectMap.h"
#include "CGO.h"
#include"Word.h"
#include"Symmetry.h"
#include"Result.h"

struct ObjectMeshState : public CObjectState {
  ObjectNameType MapName{};
  int MapState;
  CCrystal Crystal;
  int Active = true;
  pymol::vla<int> N;
  std::vector<int> RC;
  int VCsize, base_n_V;
  int OneColor;
  pymol::vla<float> V;
  std::vector<float> VC;
  int Range[6]{};
  float ExtentMin[3]{}, ExtentMax[3]{};
  int ExtentFlag = false;
  float Level, Radius;
  int RefreshFlag;
  int ResurfaceFlag = true;
  int quiet = true;
  int RecolorFlag = false;
  pymol::vla<float> AtomVertex;
  int CarveFlag = false;
  float CarveBuffer = 0.0f;
  cIsomeshMode MeshMode;
  pymol::cache_ptr<CGO> UnitCellCGO;
  WordType caption{};
  float AltLevel;
  pymol::copyable_ptr<Isofield> Field;
  /* not stored */
  pymol::cache_ptr<CGO> shaderCGO;
  pymol::cache_ptr<CGO> shaderUnitCellCGO;
  ObjectMeshState(PyMOLGlobals* G);
};

struct ObjectMesh : public pymol::CObject {
  std::vector<ObjectMeshState> State;
  int NState = 0;
  ObjectMesh(PyMOLGlobals* G);

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  int getNFrame() const override;
  pymol::CObject* clone() const override;
};

ObjectMesh *ObjectMeshFromBox(PyMOLGlobals * G, ObjectMesh * obj, ObjectMap * map,
                              int map_state,
                              int state, float *mn, float *mx,
                              float level, cIsomeshMode,
                              float carve, float *vert_vla, float alt_level, int quiet);
ObjectMesh *ObjectMeshFromXtalSym(PyMOLGlobals * G, ObjectMesh * obj, ObjectMap * map,
                                  CSymmetry * sym,
                                  int map_state,
                                  int state, float *mn, float *mx,
                                  float level, cIsomeshMode,
                                  float carve, float *vert_vla,
                                  float alt_level, int quiet);
void ObjectMeshDump(ObjectMesh * I, const char *fname, int state, int quiet);

PyObject *ObjectMeshAsPyList(ObjectMesh * I);
int ObjectMeshNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectMesh ** result);
int ObjectMeshSetLevel(ObjectMesh * I, float level, int state, int quiet);
pymol::Result<float> ObjectMeshGetLevel(ObjectMesh * I, int state);
int ObjectMeshInvalidateMapName(ObjectMesh * I, const char *name, const char * new_name);
int ObjectMeshAllMapsInStatesExist(ObjectMesh * I);

#endif
