
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

struct ObjectCGOState {
  pymol::cache_ptr<CGO> origCGO;
  pymol::cache_ptr<CGO> renderCGO;
  PyMOLGlobals* G;
  bool renderWithShaders, hasTransparency, cgo_lighting, hasOpaque;
  ObjectCGOState(PyMOLGlobals* G);
  ObjectCGOState(const ObjectCGOState& other);
};

struct ObjectCGO : public pymol::CObject {
  std::vector<ObjectCGOState> State;
  ObjectCGO(PyMOLGlobals* G);
  ObjectCGO(const ObjectCGO& other);

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  pymol::CObject* clone() const override;
  int getNFrame() const override;
};

ObjectCGO *ObjectCGODefine(PyMOLGlobals * G, ObjectCGO * obj, PyObject * pycgo,
                           int state);
ObjectCGO *ObjectCGOFromFloatArray(PyMOLGlobals * G, ObjectCGO * obj, float *array,
                                   int size, int state, int quiet);
ObjectCGO *ObjectCGOFromCGO(PyMOLGlobals * G, ObjectCGO * obj, CGO * cgo, int state);
void ObjectCGORecomputeExtent(ObjectCGO * I);

PyObject *ObjectCGOAsPyList(ObjectCGO * I);
int ObjectCGONewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectCGO ** result,
			       int version);
ObjectCGO *ObjectCGONewVFontTest(PyMOLGlobals * G, const char *text, float *pos);

#endif
