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
#ifndef _H_ObjectMap
#define _H_ObjectMap

#include"os_python.h"

#include"PyMOLObject.h"
#include"Crystal.h"
#include"Isosurf.h"
#include"CGO.h"

#define cMapSourceUndefined 0
#define cMapSourceXPLOR 1
#define cMapSourceCCP4 2
#define cMapSourcePHI 3
#define cMapSourceDesc 4
#define cMapSourceFLD 5
#define cMapSourceBRIX 6
#define cMapSourceGRD 7
#define cMapSourceChempyBrick 8
#define cMapSourceVMDPlugin 9

typedef struct ObjectMapState {
  CObjectState State;
  int Active;
  CCrystal *Crystal;
  int Div[3]; /* NOTE: Div is only meaningful for maps defined relative to a unit cell */
  int Min[3],Max[3]; /* valid min and max indices, required. */
  int FDim[4]; /* Array dimensions with 3 in fourth slot, required */
  int MapSource;
  Isofield *Field;
  float Corner[24];
  int *Dim; /* this field is redundant and should be eliminated -- if exists, must match FDim */
  float *Origin; /* Origin for non-xtal maps */
  float *Range; /* Range for non-xtal maps */
  float *Grid; /* Spacing for non-xtal maps */
  float ExtentMin[3],ExtentMax[3];

  /* below not stored */

  int have_range;
  float high_cutoff,low_cutoff;
} ObjectMapState;

typedef struct ObjectMap {
  CObject Obj;
  ObjectMapState *State;
  int NState;
} ObjectMap;

#define cObjectMap_OrthoMinMaxGrid 0

typedef struct ObjectMapDesc { /* information for creating a new map */
  int mode; 
  float Grid[3];
  int Dim[3];
  float MinCorner[3],MaxCorner[3];  
  int init_mode; /* -1 = nothing
                     0 = zeros
                     1 = ones */
} ObjectMapDesc;

ObjectMap *ObjectMapNew(PyMOLGlobals *G);
int ObjectMapNewCopy(PyMOLGlobals *G,ObjectMap *src,ObjectMap **result,int source_state, int target_state);
ObjectMapState *ObjectMapNewStateFromDesc(PyMOLGlobals *G,ObjectMap *I,ObjectMapDesc *md,int state,int quiet);
int ObjectMapStateGetExcludedStats(PyMOLGlobals *G,ObjectMapState *ms,float *vert_vla,
                                   float beyond, float within, float *level);

int ObjectMapValidXtal(ObjectMap *I, int state);

int ObjectMapStateGetRange(PyMOLGlobals *G,ObjectMapState *ms,float *range);

void ObjectMapStateClamp(ObjectMapState *I,float clamp_floor, float clamp_ceiling);

ObjectMap *ObjectMapLoadXPLOR(PyMOLGlobals *G,ObjectMap *obj,char *fname,
                                  int state,int is_file,int quiet);

ObjectMap *ObjectMapLoadCCP4(PyMOLGlobals *G,ObjectMap *obj,char *fname,
                             int state,int is_string,int bytes,int quiet);

ObjectMap *ObjectMapLoadPHI(PyMOLGlobals *G,ObjectMap *obj,char *fname,int state,
                            int is_string,int bytes, int quiet);

ObjectMap *ObjectMapLoadDXFile(PyMOLGlobals *G,ObjectMap *obj,char *fname,int state,int quiet);
ObjectMap *ObjectMapLoadFLDFile(PyMOLGlobals *G,ObjectMap *obj,char *fname,int state,int quiet);
ObjectMap *ObjectMapLoadBRIXFile(PyMOLGlobals *G,ObjectMap *obj,char *fname,int state,int quiet);
ObjectMap *ObjectMapLoadGRDFile(PyMOLGlobals *G,ObjectMap *obj,char *fname,int state,int quiet);

ObjectMap *ObjectMapLoad(PyMOLGlobals *G,ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapLoadChemPyBrick(PyMOLGlobals *G,ObjectMap *I,PyObject *Map,
                                    int state,int discrete,int quiet);
ObjectMap *ObjectMapLoadCObject(PyMOLGlobals *G,ObjectMap *obj,int state);
ObjectMap *ObjectMapLoadChemPyMap(PyMOLGlobals *G,ObjectMap *I,PyObject *Map,
                                  int state,int discrete,int quiet);
int ObjectMapDouble(ObjectMap *I,int state);
int ObjectMapHalve(ObjectMap *I,int state,int smooth);
int ObjectMapTrim(ObjectMap *I,int state, float *mn, float *mx,int quiet);
int ObjectMapSetBorder(ObjectMap *I,float level,int state);
int ObjectMapStateSetBorder(ObjectMapState *I,float level);
void ObjectMapStateInit(PyMOLGlobals *G,ObjectMapState *I);
void ObjectMapStatePurge(PyMOLGlobals *G,ObjectMapState *I);
int ObjectMapStateInterpolate(ObjectMapState *ms,float *array,float *result,int *flag, int n);
int ObjectMapStateContainsPoint(ObjectMapState *ms,float *point);
ObjectMapState *ObjectMapStatePrime(ObjectMap *I,int state);
ObjectMapState *ObjectMapStateGetActive(ObjectMap *I,int state);
int ObjectMapGetNStates(ObjectMap *I);
void ObjectMapUpdateExtents(ObjectMap *I);
ObjectMapState *ObjectMapGetState(ObjectMap *I,int state);

PyObject *ObjectMapAsPyList(ObjectMap *I);
int ObjectMapNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectMap **result);

int ObjectMapInterpolate(ObjectMap *I,int state,float *array,float *result,int *flag,int n);

void ObjectMapTransformMatrix(ObjectMap *I, int state, double *matrix);
void ObjectMapResetMatrix(ObjectMap *I, int state);
int ObjectMapGetMatrix(ObjectMap *I,int state,double **matrix);
int ObjectMapSetMatrix(ObjectMap *I,int state,double *matrix);
void ObjectMapStateRegeneratePoints(ObjectMapState *ms);
int ObjectMapStateGetDataRange(PyMOLGlobals *G,ObjectMapState *ms, float *min, float *max);

#endif











