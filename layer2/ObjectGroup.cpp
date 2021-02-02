
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectGroup.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"PConv.h"

int ObjectGroupNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectGroup ** result,
                             int version)
{
  int ok = true, ll = 0;
  ObjectGroup *I = NULL;
  (*result) = NULL;
  if(ok)
    ok = (list != Py_None);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  I = new ObjectGroup(G);
  if(ok)
    ok = (I != NULL);
  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if(ok)
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 1, &I->OpenOrClosed);
  if(ok && (ll > 2)){
    // State removed because it was unused
  }
  if(ok) {
    *result = I;
  } else {
    /* to do: cleanup */
  }
  return (ok);
}

PyObject *ObjectGroupAsPyList(ObjectGroup * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->OpenOrClosed));

  // State removed in PyMOL 2.4.0
  CObjectState tmp_state(I->G);
  PyList_SetItem(result, 2, ObjectStateAsPyList(&tmp_state));

  return (PConvAutoNone(result));
}


/*========================================================================*/

ObjectGroup::~ObjectGroup()
{
}


/*========================================================================*/
ObjectGroup::ObjectGroup(PyMOLGlobals * G) : pymol::CObject(G)
{
  auto I = this;
  I->type = cObjectGroup;
}
