
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
#ifndef _H_ObjectCallback
#define _H_ObjectCallback

#include"os_python.h"

#include"PyMOLObject.h"

typedef struct {
  PyObject *PObj;
} ObjectCallbackState;

typedef struct ObjectCallback {
  CObject Obj;
  ObjectCallbackState *State;
  int NState;
} ObjectCallback;

ObjectCallback *ObjectCallbackNew(PyMOLGlobals * G);
ObjectCallback *ObjectCallbackDefine(PyMOLGlobals * G, ObjectCallback * obj,
                                     PyObject * PObj, int state);
void ObjectCallbackRecomputeExtent(ObjectCallback * I);

#endif
