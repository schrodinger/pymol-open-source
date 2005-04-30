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

typedef struct {
  CObjectState State;
  char MapName[ObjNameMax];
  int MapState;
  CCrystal Crystal;
  int Active;
  int *N;
  float *V;
  int Range[6];
  float ExtentMin[3],ExtentMax[3];
  int ExtentFlag;
  float Level,Radius;
  int RefreshFlag;
  int ResurfaceFlag;
  float *AtomVertex;
  int CarveFlag;
  float CarveBuffer;
  int DotFlag;
  CGO *UnitCellCGO;
  int displayList;
} ObjectMeshState;

typedef struct ObjectMesh {
  CObject Obj;
  ObjectMeshState *State;
  int NState;
} ObjectMesh;

ObjectMesh *ObjectMeshFromBox(PyMOLGlobals *G,ObjectMesh *obj,ObjectMap* map,
                              int map_state,
                              int state,float *mn,float *mx,
                              float level,int dotFlag,
                              float carve,float *vert_vla);
void ObjectMeshDump(ObjectMesh *I,char *fname,int state);

PyObject *ObjectMeshAsPyList(ObjectMesh *I);
int ObjectMeshNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectMesh **result);
int ObjectMeshSetLevel(ObjectMesh *I,float level,int state);

#endif

