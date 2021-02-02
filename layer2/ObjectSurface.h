
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
#ifndef _H_ObjectSurface
#define _H_ObjectSurface

#include"os_gl.h"
#include"ObjectMap.h"
#include"Result.h"
#include"CGO.h"
#include"PyMOLEnums.h"

struct ObjectSurfaceState : public CObjectState
{
  ObjectNameType MapName;
  int MapState;
  CCrystal Crystal;
  int Active = false;
  pymol::vla<int> N;
  int nT = 0;
  int base_n_V;
  pymol::vla<float> V;
  std::vector<float> VC;
  std::vector<int> RC;
  int OneColor;
  int VCsize() { return VC.size() / 3; }
  int Range[6];
  float ExtentMin[3], ExtentMax[3];
  int ExtentFlag = false;
  float Level, Radius;
  int RefreshFlag;
  int ResurfaceFlag = true;
  int RecolorFlag = false;
  int quiet = true;
  pymol::vla<float> AtomVertex;
  int CarveFlag = false;
  float CarveBuffer;
  cIsosurfaceMode Mode;
  int DotFlag;
  pymol::cache_ptr<CGO> UnitCellCGO;
  cIsosurfaceSide Side = cIsosurfaceSide::front;
  pymol::cache_ptr<CGO> shaderCGO;
  ObjectSurfaceState(PyMOLGlobals* G);
};

struct ObjectSurface : public pymol::CObject {
  std::vector<ObjectSurfaceState> State;
  ObjectSurface(PyMOLGlobals* G);

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  int getNFrame() const override;
  pymol::CObject* clone() const override;
};

ObjectSurface *ObjectSurfaceFromBox(PyMOLGlobals * G, ObjectSurface * obj,
                                    ObjectMap * map, int map_state, int state, float *mn,
                                    float *mx, float level, cIsosurfaceMode, float carve,
                                    pymol::vla<float>&& vert_vla, cIsosurfaceSide, int quiet);
void ObjectSurfaceDump(ObjectSurface * I, const char *fname, int state, int quiet);

int ObjectSurfaceNewFromPyList(PyMOLGlobals * G, PyObject * list,
				   ObjectSurface ** result);
PyObject *ObjectSurfaceAsPyList(ObjectSurface * I);
int ObjectSurfaceSetLevel(ObjectSurface * I, float level, int state, int quiet);
pymol::Result<float> ObjectSurfaceGetLevel(ObjectSurface * I, int state);
int ObjectSurfaceInvalidateMapName(ObjectSurface * I, const char *name, const char * new_name);

#endif
