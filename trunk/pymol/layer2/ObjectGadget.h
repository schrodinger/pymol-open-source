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
#ifndef _H_ObjectGadget
#define _H_ObjectGadget

#include"os_python.h"

#include"PyMOLObject.h"
#include"CGO.h"

#include"GadgetSet.h"


typedef struct ObjectGadgetState {
} ObjectGadgetState;

typedef struct ObjectGadget {
  CObject Obj;
  struct GadgetSet **GSet;
  int NGSet;
  int CurGSet;
  int GadgetType;
} ObjectGadget;


#define cGadget     0
#define cGadgetRamp 1

ObjectGadget *ObjectGadgetNew(void);
void ObjectGadgetInit(ObjectGadget *I);
void ObjectGadgetPurge(ObjectGadget *I);
void ObjectGadgetFree(ObjectGadget *I);
ObjectGadget *ObjectGadgetDefine(ObjectGadget *obj,PyObject *pycgo,int state);
ObjectGadget *ObjectGadgetFromCGO(ObjectGadget *obj,CGO *cgo,int state);
void ObjectGadgetRecomputeExtent(ObjectGadget *I);

PyObject *ObjectGadgetAsPyList(ObjectGadget *I);
int ObjectGadgetNewFromPyList(PyObject *list,ObjectGadget **result);
ObjectGadget *ObjectGadgetTest(void);
int ObjectGadgetGetVertex(ObjectGadget *I,int index,int base, float *v); /* in current state */
int ObjectGadgetSetVertex(ObjectGadget *I,int index,int base, float *v); /* in current state */
void ObjectGadgetUpdate(ObjectGadget *I);
void ObjectGadgetUpdateExtents(ObjectGadget *I);
void ObjectGadgetUpdateStates(ObjectGadget *I);

#endif











