
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
#include"Base.h"

#include"GadgetSet.h"

struct ObjectGadget : public pymol::CObject {
  pymol::vla<GadgetSet*> GSet;
  int NGSet = 0;
  int CurGSet = 0;
  int GadgetType;
  bool Changed = true;
  ObjectGadget(PyMOLGlobals* G);
  virtual ~ObjectGadget();

  // virtual methods
  virtual void update() override;
  void render(RenderInfo* info) override;
  int getNFrame() const override;
  pymol::RenderContext getRenderContext() const override;
};

#define cGadgetPlain 0
#define cGadgetRamp 1

PyObject *ObjectGadgetAsPyList(ObjectGadget * I);
PyObject *ObjectGadgetPlainAsPyList(ObjectGadget * I, bool incl_cgos=true);

int ObjectGadgetNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGadget ** result,
                              int version);
int ObjectGadgetInitFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGadget * I,
                               int version);

ObjectGadget *ObjectGadgetTest(PyMOLGlobals * G);
int ObjectGadgetGetVertex(const ObjectGadget * I, int index, int base, float *v);     /* in current state */
int ObjectGadgetSetVertex(ObjectGadget * I, int index, int base, const float *v);     /* in current state */
void ObjectGadgetUpdateExtents(ObjectGadget * I);
void ObjectGadgetUpdateStates(ObjectGadget * I);

#endif
