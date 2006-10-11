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
#ifndef _H_ObjectGroup
#define _H_ObjectGroup

#include"os_python.h"

#include"PyMOLObject.h"


typedef struct ObjectGroup {
  CObject Obj;
} ObjectGroup;


ObjectGroup *ObjectGroupNew(PyMOLGlobals *G);

void ObjectGroupRecomputeExtent(ObjectGroup *I);

PyObject *ObjectGroupAsPyList(ObjectGroup *I);

int ObjectGroupNewFromPyList(PyMOLGlobals *G,PyObject *list,
                                 ObjectGroup **result,int version);


#endif











