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

#include"ObjectMap.h"

typedef struct {
  CObjectState State;
  ObjectNameType MapName;
  int MapState;
  CCrystal Crystal;
  int Active;
  int *N,nT,base_n_V;
  float *V;
  float *VC;
  int *RC;
  int OneColor;
  int VCsize;
  int Range[6];
  float ExtentMin[3],ExtentMax[3];
  int ExtentFlag;
  float Level,Radius;
  int RefreshFlag;
  int ResurfaceFlag;
  int RecolorFlag;
  float *AtomVertex;
  int CarveFlag;
  float CarveBuffer;
  int Mode; /* 0 dots, 1 lines, 2 triangles */
  int DotFlag;
  CGO *UnitCellCGO;
  int Side;
  int displayList;
  int displayListInvalid;
} ObjectSurfaceState;

typedef struct ObjectSurface {
  CObject Obj;
  ObjectSurfaceState *State;
  int NState;
} ObjectSurface;

ObjectSurface *ObjectSurfaceFromBox(PyMOLGlobals *G,ObjectSurface *obj,ObjectMap *map,
                                    int map_state,
                              int state,float *mn,float *mx,
                              float level,int mode,
                              float carve,float *vert_vla,int side);
void ObjectSurfaceDump(ObjectSurface *I,char *fname,int state);

int ObjectSurfaceNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectSurface **result);
PyObject *ObjectSurfaceAsPyList(ObjectSurface *I);
int ObjectSurfaceSetLevel(ObjectSurface *I,float level,int state);
int ObjectSurfaceInvalidateMapName(ObjectSurface *I,char *name);

#endif

