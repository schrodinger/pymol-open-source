
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
#include"OVOneToAny.h"

typedef struct ObjectAlignmentState {
  CObjectState state;
  int *alignVLA;
  WordType guide;
  /* not stored */
  int valid;
  OVOneToAny *id2tag;
  CGO *std;
  CGO *ray;
} ObjectAlignmentState;

typedef struct ObjectAlignment {
  CObject Obj;
  ObjectAlignmentState *State;
  int NState;
  int SelectionState, ForceState;
} ObjectAlignment;

void ObjectAlignmentUpdate(ObjectAlignment * I);

ObjectAlignment *ObjectAlignmentDefine(PyMOLGlobals * G,
                                       ObjectAlignment * obj,
                                       int *align_vla,
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
