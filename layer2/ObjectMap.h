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

typedef struct ObjectMapState {
  int Active;
  CCrystal *Crystal;
  int Div[3],Min[3],Max[3],FDim[4];
  int MapSource;
  Isofield *Field;
  float Corner[8][3];
  int *Dim;
  float *Origin;
  float *Range;
  float *Grid;
  float ExtentMin[3],ExtentMax[3];
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

ObjectMap *ObjectMapNew(void);
ObjectMapState *ObjectMapNewStateFromDesc(ObjectMap *I,ObjectMapDesc *md,int state);
int ObjectMapStateGetExcludedStats(ObjectMapState *ms,float *vert_vla,
                                   float beyond, float within, float *level);

ObjectMap *ObjectMapLoadXPLORFile(ObjectMap *obj,char *fname,int state,int is_file);
ObjectMap *ObjectMapReadXPLORStr(ObjectMap *I,char *XPLORStr,int state);
int ObjectMapXPLORStrToMap(ObjectMap *I,char *XPLORStr,int state);

ObjectMap *ObjectMapLoadCCP4File(ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapReadCCP4Str(ObjectMap *I,char *XPLORStr,int bytes,int state);
int ObjectMapCCP4StrToMap(ObjectMap *I,char *XPLORStr,int bytes,int state);

ObjectMap *ObjectMapLoadDXFile(ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapLoadPHIFile(ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapLoadFLDFile(ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapLoadBRIXFile(ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapLoadGRDFile(ObjectMap *obj,char *fname,int state);

ObjectMap *ObjectMapLoad(ObjectMap *obj,char *fname,int state);
ObjectMap *ObjectMapLoadChemPyBrick(ObjectMap *I,PyObject *Map,
                                    int state,int discrete);
ObjectMap *ObjectMapLoadCObject(ObjectMap *obj,int state);
ObjectMap *ObjectMapLoadChemPyMap(ObjectMap *I,PyObject *Map,
                                  int state,int discrete);
int ObjectMapDouble(ObjectMap *I,int state);
int ObjectMapSetBorder(ObjectMap *I,float level);
int ObjectMapStateSetBorder(ObjectMapState *I,float level);
void ObjectMapStateInit(ObjectMapState *I);
void ObjectMapStatePurge(ObjectMapState *I);
int ObjectMapStateInterpolate(ObjectMapState *ms,float *array,float *result,int *flag, int n);
int ObjectMapStateContainsPoint(ObjectMapState *ms,float *point);
ObjectMapState *ObjectMapStatePrime(ObjectMap *I,int state);
ObjectMapState *ObjectMapStateGetActive(ObjectMap *I,int state);
int ObjectMapGetNStates(ObjectMap *I);
void ObjectMapUpdateExtents(ObjectMap *I);
ObjectMapState *ObjectMapGetState(ObjectMap *I,int state);

PyObject *ObjectMapAsPyList(ObjectMap *I);
int ObjectMapNewFromPyList(PyObject *list,ObjectMap **result);

int ObjectMapInterpolate(ObjectMap *I,int state,float *array,float *result,int *flag,int n);

#endif











