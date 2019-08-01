
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
  I = ObjectGroupNew(G);
  if(ok)
    ok = (I != NULL);
  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->OpenOrClosed);
  if(ok && (ll > 2))
    ok = ObjectStateFromPyList(G, PyList_GetItem(list, 2), &I->State);
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
  PyList_SetItem(result, 2, ObjectStateAsPyList(&I->State));
  return (PConvAutoNone(result));
}


/*========================================================================*/

static void ObjectGroupFree(ObjectGroup * I)
{
  ObjectStatePurge(&I->State);
  ObjectPurge(I);
  OOFreeP(I);
}


/*========================================================================*/
static CObjectState *ObjectGroupGetObjectState(ObjectGroup * I, int state)
{
  return &I->State;
}


/*========================================================================*/
ObjectGroup *ObjectGroupNew(PyMOLGlobals * G)
{
  OOCalloc(G, ObjectGroup);

  ObjectInit(G, (CObject *) I);

  I->type = cObjectGroup;
  I->fFree = (void (*)(CObject *)) ObjectGroupFree;
  I->fRender = NULL;
  I->OpenOrClosed = false;
  I->fGetObjectState = (CObjectState * (*)(CObject *, int state))
    ObjectGroupGetObjectState;

  ObjectStateInit(G, &I->State);
  return (I);
}

void ObjectGroupResetMatrix(ObjectGroup * I, int state)
{
  ObjectStateResetMatrix(&I->State);
}

int ObjectGroupGetMatrix(ObjectGroup * I, int state, double **matrix)
{
  *matrix = ObjectStateGetMatrix(&I->State);
  return true;
}

int ObjectGroupSetMatrix(ObjectGroup * I, int state, double *matrix)
{
  ObjectStateSetMatrix(&I->State, matrix);
  return true;
}

void ObjectGroupTransformMatrix(ObjectGroup * I, int state, double *matrix)
{
  ObjectStateTransformMatrix(&I->State, matrix);
}
