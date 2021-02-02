
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
  bool is_callable;
} ObjectCallbackState;

struct ObjectCallback : public pymol::CObject {
  ObjectCallbackState *State = nullptr;
  int NState = 0;
  ObjectCallback(PyMOLGlobals* G);
  ~ObjectCallback();

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  int getNFrame() const override;
};

ObjectCallback *ObjectCallbackDefine(PyMOLGlobals * G, ObjectCallback * obj,
                                     PyObject * PObj, int state);
void ObjectCallbackRecomputeExtent(ObjectCallback * I);

#ifndef _PYMOL_NOPY
PyObject *ObjectCallbackAsPyList(ObjectCallback * I);
int ObjectCallbackNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectCallback ** result);
#endif

#endif
