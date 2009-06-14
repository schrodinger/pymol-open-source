
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

typedef struct ObjectGadget {
  CObject Obj;
  struct GadgetSet **GSet;
  int NGSet;
  int CurGSet;
  int GadgetType;
  int Changed;
} ObjectGadget;

#define cGadgetPlain 0
#define cGadgetRamp 1

ObjectGadget *ObjectGadgetNew(PyMOLGlobals * G);
void ObjectGadgetInit(PyMOLGlobals * G, ObjectGadget * I);
void ObjectGadgetPurge(ObjectGadget * I);
void ObjectGadgetFree(ObjectGadget * I);
ObjectGadget *ObjectGadgetDefine(PyMOLGlobals * G, ObjectGadget * obj, PyObject * pycgo,
                                 int state);
ObjectGadget *ObjectGadgetFromCGO(PyMOLGlobals * G, ObjectGadget * obj, CGO * cgo,
                                  int state);
void ObjectGadgetRecomputeExtent(ObjectGadget * I);

PyObject *ObjectGadgetAsPyList(ObjectGadget * I);
PyObject *ObjectGadgetAsPyList(ObjectGadget * I);
PyObject *ObjectGadgetPlainAsPyList(ObjectGadget * I);

int ObjectGadgetNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGadget ** result,
                              int version);
int ObjectGadgetInitFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGadget * I,
                               int version);

ObjectGadget *ObjectGadgetTest(PyMOLGlobals * G);
int ObjectGadgetGetVertex(ObjectGadget * I, int index, int base, float *v);     /* in current state */
int ObjectGadgetSetVertex(ObjectGadget * I, int index, int base, float *v);     /* in current state */
void ObjectGadgetUpdate(ObjectGadget * I);
void ObjectGadgetUpdateExtents(ObjectGadget * I);
void ObjectGadgetUpdateStates(ObjectGadget * I);

#endif
