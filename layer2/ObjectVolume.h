/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright (C) by Schrodinger.
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
#ifndef _H_ObjectVolume
#define _H_ObjectVolume

#include"CGO.h"
#include"ObjectMap.h"
#include"Word.h"
#include"Symmetry.h"
#include"Util2.h"
#include"Result.h"

struct ObjectVolumeState : public CObjectState {
  ObjectNameType MapName;
  int MapState;
  CCrystal Crystal;
  int Active = false;
  int Range[6];
  float ExtentMin[3], ExtentMax[3];
  int ExtentFlag = false;
  // TODO difference between Resurface, Recolor, Refresh???
  int RefreshFlag;
  int ResurfaceFlag = true;
  int RecolorFlag = true;
  pymol::vla<float> AtomVertex;
  float CarveBuffer = 0.0f;
  WordType caption{};
  float Corner[24];
  /* not stored */
  pymol::cache_array<std::size_t, 3> textures{}; // 3D volume (map), 1D/2D color table, 3D carvemask
  pymol::copyable_ptr<CField> carvemask;
  unsigned int dim[3]{};
  pymol::copyable_ptr<Isofield> Field;
  float min_max_mean_stdev[4];
  float ramp_min, ramp_range;
  int RampSize() const { return Ramp.size() / 5; };
  std::vector<float> Ramp;
  int isUpdated = false;

  ObjectVolumeState(PyMOLGlobals* G);
  ~ObjectVolumeState();
};

struct ObjectVolume : public pymol::CObject {
  std::vector<ObjectVolumeState> State;
  ObjectVolume(PyMOLGlobals* G);

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  int getNFrame() const override;
  pymol::CObject* clone() const override;
};

ObjectVolume *ObjectVolumeFromBox(PyMOLGlobals * G, ObjectVolume * obj, ObjectMap * map,
                              int map_state,
                              int state, float *mn, float *mx,
                              float level, int meshMode,
                              float carve, float *vert_vla, int quiet);
ObjectVolume *ObjectVolumeFromXtalSym(PyMOLGlobals * G, ObjectVolume * obj, ObjectMap * map,
                                  CSymmetry * sym,
                                  int map_state,
                                  int state, float *mn, float *mx,
                                  float level, int meshMode,
                                  float carve, float *vert_vla,
                                  int quiet);

PyObject *ObjectVolumeAsPyList(ObjectVolume * I);
int ObjectVolumeNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectVolume ** result);
int ObjectVolumeInvalidateMapName(ObjectVolume * I, const char *name, const char * new_name);

CField   * ObjectVolumeGetField(ObjectVolume* I);
PyObject * ObjectVolumeGetRamp(ObjectVolume* I);
pymol::Result<>  ObjectVolumeSetRamp(ObjectVolume* I, std::vector<float>&& ramp_list);

ObjectMapState * ObjectVolumeGetMapState(ObjectVolume * I);

#endif
