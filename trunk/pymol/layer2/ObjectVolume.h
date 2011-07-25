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
  int *N, *RC, VCsize, base_n_V;
  int OneColor;
  float *V, *VC;
  int Range[6];
  float ExtentMin[3], ExtentMax[3];
  int ExtentFlag;
  float Level, Radius;
  int RefreshFlag;
  int ResurfaceFlag;
  int quiet;
  int RecolorFlag;
  float *AtomVertex;
  int CarveFlag;
  float CarveBuffer;
  int VolumeMode;
  CGO *UnitCellCGO;
  int displayList;
  int displayListInvalid;
  WordType caption;
  float AltLevel;
  float Corner[24];
  int textures[3];
  CField *volume;
  float *Histogram;
  /* not stored */
  float *colors;
  Isofield *Field; 
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
void ObjectVolumeDump(ObjectVolume * I, char *fname, int state);

PyObject *ObjectVolumeAsPyList(ObjectVolume * I);
int ObjectVolumeNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectVolume ** result);
int ObjectVolumeSetLevel(ObjectVolume * I, float level, int state, int quiet);
int ObjectVolumeGetLevel(ObjectVolume * I, int state, float *result);
int ObjectVolumeInvalidateMapName(ObjectVolume * I, char *name);

int ObjectVolumeColor(ObjectVolume * I, float * colors, int ncolors);

PyObject * ObjectVolumeGetField(ObjectVolume* I);
PyObject * ObjectVolumeGetHistogram(ObjectVolume* I);
PyObject * ObjectVolumeGetRamp(ObjectVolume* I);
PyObject * ObjectVolumeSetRamp(ObjectVolume* I, float *ramp_list, int list_size);
PyObject * ObjectVolumeGetIsUpdated(ObjectVolume *I);

#endif
