
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
#ifndef _H_ObjectSlice
#define _H_ObjectSlice

#include"ObjectMap.h"

struct ObjectSliceState {
  PyMOLGlobals *G;
  /* stored in a session */

  int Active = true;
  ObjectNameType MapName{};
  int MapState = 0;
  float MapMean = 0.0f;
  float MapStdev = 0.0f;

  float ExtentMin[3]{};
  float ExtentMax[3]{};
  int ExtentFlag = false;

  float origin[3]{};              /* the origin of the plane */
  float system[9]{              /* x, y, and z of the system */
    1.f, 0.f, 0.f,
    0.f, 1.f, 0.f,
    0.f, 0.f, 1.0f};

  /* not stored in session */

  int RefreshFlag = true;
  int min[2]{}, max[2]{};           /* extents of the arrays */
  float last_scale = 0.0f;

  /* the data is normalized for easier ploting */
  int n_points = 0;

  pymol::vla<float> values;
  pymol::vla<float> points;
  pymol::vla<int> flags;
  pymol::vla<float> colors;
  pymol::vla<float> normals;

  int n_strips = 0;
  pymol::vla<int> strips;

  pymol::cache_ptr<CGO> shaderCGO;
  float Corner[24];

  float outline_points[36];
  int outline_n_points = 0;
  float outline_zaxis[3];
  ObjectSliceState(PyMOLGlobals* G)
      : G(G){};
};

struct ObjectSlice : public pymol::CObject {
  std::vector<ObjectSliceState> State;
  PickContext context{};
  ObjectSlice(PyMOLGlobals* G);

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  int getNFrame() const override;
  pymol::CObject* clone() const override;
};

ObjectSlice *ObjectSliceFromMap(PyMOLGlobals * G, ObjectSlice * obj, ObjectMap * map,
                                int state, int map_state);

PyObject *ObjectSliceAsPyList(ObjectSlice * I);
int ObjectSliceNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectSlice ** result);

ObjectSliceState *ObjectSliceStateGetActive(ObjectSlice * I, int state);
void ObjectSliceDrag(ObjectSlice * I, int state, int mode, float *pt, float *mov,
                     float *z_dir);
int ObjectSliceGetVertex(ObjectSlice * I, int index, int base, float *v);

#endif
