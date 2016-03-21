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

#include"ObjectMap.h"
#include"Word.h"
#include"Symmetry.h"

typedef struct {
  CObjectState State;
  ObjectNameType MapName;
  int MapState;
  CCrystal Crystal;
  int Active;
  int Range[6];
  float ExtentMin[3], ExtentMax[3];
  int ExtentFlag;
  // TODO difference between Resurface, Recolor, Refresh???
  int RefreshFlag;
  int ResurfaceFlag;
  int RecolorFlag;
  float *AtomVertex;
  float CarveBuffer;
  CGO *UnitCellCGO;
  WordType caption;
  float Corner[24];
  /* not stored */
  int textures[3];
  CField *carvemask;
  unsigned int dim[3];
  Isofield *Field;
  float min_max_mean_stdev[4];
  float ramp_min, ramp_range;
  int RampSize;
  float *Ramp;
  int isUpdated; 
} ObjectVolumeState;

typedef struct ObjectVolume {
  CObject Obj;
  ObjectVolumeState *State;
  int NState;
} ObjectVolume;

ObjectVolume *ObjectVolumeFromBox(PyMOLGlobals * G, ObjectVolume * obj, ObjectMap * map,
                              int map_state,
                              int state, float *mn, float *mx,
                              float level, int meshMode,
                              float carve, float *vert_vla, float alt_level, int quiet);
ObjectVolume *ObjectVolumeFromXtalSym(PyMOLGlobals * G, ObjectVolume * obj, ObjectMap * map,
                                  CSymmetry * sym,
                                  int map_state,
                                  int state, float *mn, float *mx,
                                  float level, int meshMode,
                                  float carve, float *vert_vla,
                                  float alt_level, int quiet);

PyObject *ObjectVolumeAsPyList(ObjectVolume * I);
int ObjectVolumeNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectVolume ** result);
int ObjectVolumeInvalidateMapName(ObjectVolume * I, const char *name, const char * new_name);

int ObjectVolumeColor(ObjectVolume * I, float * colors, int ncolors);

CField   * ObjectVolumeGetField(ObjectVolume* I);
PyObject * ObjectVolumeGetRamp(ObjectVolume* I);
int        ObjectVolumeSetRamp(ObjectVolume* I, float *ramp_list, int list_size);

ObjectMapState * ObjectVolumeGetMapState(ObjectVolume * I);

#endif
