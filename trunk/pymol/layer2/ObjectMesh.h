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

typedef struct ObjectMesh {
  Object Obj;
  ObjectMap *Map;
  int *N;
  float *V;
  int Range[6];
  float Level;
  int ResurfaceFlag;
} ObjectMesh;

ObjectMesh *ObjectMeshFromBox(ObjectMap *map,float *mn,float *mx,float level);

#endif











