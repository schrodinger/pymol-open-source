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
#ifndef _H_ObjectCGO
#define _H_ObjectCGO

#include"os_python.h"

#include"PyMOLObject.h"
#include"CGO.h"

typedef struct ObjectCGOState {
  CGO *std;
  CGO *ray;
} ObjectCGOState;

typedef struct ObjectCGO {
  CObject Obj;
  ObjectCGOState *State;
  int NState;
} ObjectCGO;

ObjectCGO *ObjectCGONew(PyMOLGlobals *G);
ObjectCGO *ObjectCGODefine(PyMOLGlobals *G,ObjectCGO *obj,PyObject *pycgo,int state);
ObjectCGO *ObjectCGOFromFloatArray(PyMOLGlobals *G,ObjectCGO *obj,float *array, int size, int state, int quiet);
ObjectCGO *ObjectCGOFromCGO(PyMOLGlobals *G,ObjectCGO *obj,CGO *cgo,int state);
void ObjectCGORecomputeExtent(ObjectCGO *I);

PyObject *ObjectCGOAsPyList(ObjectCGO *I);
int ObjectCGONewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectCGO **result,int version);
ObjectCGO *ObjectCGONewVFontTest(PyMOLGlobals *G,char *text,float *pos);

#endif











