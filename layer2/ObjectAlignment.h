
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2006 by Warren Lyford Delano of DeLano Scientific. 
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
#ifndef _H_ObjectAlignment
#define _H_ObjectAlignment

#include"os_python.h"

#include"PyMOLObject.h"
#include"CGO.h"
#include"ObjectMolecule.h"
#include"pymol/memory.h"

struct ObjectAlignmentState {
  pymol::vla<int> alignVLA;
  WordType guide;
  /* not stored */
  int valid;
  std::unordered_map<int, int> id2tag;
  pymol::cache_ptr<CGO> primitiveCGO;
  pymol::cache_ptr<CGO> renderCGO;
  bool renderCGO_has_cylinders;
  bool renderCGO_has_trilines;
};

struct ObjectAlignment : public pymol::CObject {
  std::vector<ObjectAlignmentState> State;
  int SelectionState = -1;
  int ForceState = -1;
  ObjectAlignment(PyMOLGlobals* G);

  // virtual methods
  void update() override;
  void render(RenderInfo* info) override;
  void invalidate(cRep_t rep, cRepInv_t level, int state) override;
  int getNFrame() const override;
  pymol::CObject* clone() const override;
};

ObjectAlignment *ObjectAlignmentDefine(PyMOLGlobals * G,
                                       ObjectAlignment * obj,
                                       const pymol::vla<int>& align_vla,
                                       int state,
                                       int merge,
                                       ObjectMolecule * guide, ObjectMolecule * flush);

void ObjectAlignmentRecomputeExtent(ObjectAlignment * I);

PyObject *ObjectAlignmentAsPyList(ObjectAlignment * I);

int ObjectAlignmentNewFromPyList(PyMOLGlobals * G, PyObject * list,
				     ObjectAlignment ** result, int version);

int ObjectAlignmentAsStrVLA(PyMOLGlobals * G, ObjectAlignment * I, int state, int format,
                            char **str_vla);

#endif
