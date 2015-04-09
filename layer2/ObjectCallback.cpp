
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectCallback.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"P.h"
#include"Scene.h"
#include"PConv.h"
#include"main.h"
#include"Setting.h"

static void ObjectCallbackFree(ObjectCallback * I);


/*========================================================================*/

static void ObjectCallbackFree(ObjectCallback * I)
{
#ifndef _PYMOL_NOPY
  int a;
  PyMOLGlobals *G = I->Obj.G;
  int blocked = PAutoBlock(G);
  for(a = 0; a < I->NState; a++) {
    if(I->State[a].PObj) {
      Py_DECREF(I->State[a].PObj);
      I->State[a].PObj = NULL;
    }
  }
  PAutoUnblock(G, blocked);
#endif
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}


/*========================================================================*/

static void ObjectCallbackUpdate(ObjectCallback * I)
{
  SceneInvalidate(I->Obj.G);
}


/*========================================================================*/

static void ObjectCallbackRender(ObjectCallback * I, RenderInfo * info)
{
#ifndef _PYMOL_NOPY
  int state = info->state;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int pass = info->pass;
  PyMOLGlobals *G = I->Obj.G;
  ObjectCallbackState *sobj = NULL;

  if(pass != 1) /* for now, the callback should be called during the first pass (opaque), so 
		  that it is possible to set positions for any object that is rendered in the 
		  opaque pass.  This is still a kludge, since the callback should probably 
		  happen in a pass before this (should we add a new pass 2?  this will probably
		  need changes all over.. and will introduce some minimal overhead (another pass
		  on the objects).  For now, these callbacks need to be in the object list before
		  any objects it tries to set (this came up on the crosshair.py screen stabilized example)
		  We also might want to have arguments on the callback for whether it gets called
		  at the beginning (pre-first pass) or end (post-last pass).
	       */
    return;

  if(ray || pick)
    return;

  if(!(G->HaveGUI && G->ValidContext))
    return;

  if(!I->State || I->NState == 0)
    return;

  ObjectPrepareContext(&I->Obj, ray);

  if((I->Obj.visRep & cRepCallbackBit)) {
    int blocked = PAutoBlock(G);
    for(StateIterator iter(G, I->Obj.Setting, state, I->NState); iter.next();) {
      sobj = I->State + iter.state;
      if(!sobj->is_callable)
        continue;

      Py_DecRef(PyObject_CallObject(sobj->PObj, NULL));
      if(PyErr_Occurred())
        PyErr_Print();
    }
    PAutoUnblock(G, blocked);
  }
#endif
}


/*========================================================================*/
static int ObjectCallbackGetNStates(ObjectCallback * I)
{
  return (I->NState);
}


/*========================================================================*/
ObjectCallback *ObjectCallbackNew(PyMOLGlobals * G)
{
  OOAlloc(G, ObjectCallback);

  ObjectInit(G, (CObject *) I);

  I->State = VLACalloc(ObjectCallbackState, 10);       /* autozero */
  I->NState = 0;

  I->Obj.type = cObjectCallback;
  I->Obj.fFree = (void (*)(CObject *)) ObjectCallbackFree;
  I->Obj.fUpdate = (void (*)(CObject *)) ObjectCallbackUpdate;
  I->Obj.fRender = (void (*)(CObject *, RenderInfo *))
    ObjectCallbackRender;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectCallbackGetNStates;

  return (I);
}


/*========================================================================*/
ObjectCallback *ObjectCallbackDefine(PyMOLGlobals * G, ObjectCallback * obj,
                                     PyObject * pobj, int state)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  ObjectCallback *I = NULL;

  if(!obj) {
    I = ObjectCallbackNew(G);
  } else {
    I = obj;
  }

  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCallbackState, state);
    I->NState = state + 1;
  }

  if(I->State[state].PObj) {
    Py_DECREF(I->State[state].PObj);
  }
  I->State[state].is_callable = PyCallable_Check(pobj);
  I->State[state].PObj = pobj;
  Py_INCREF(pobj);
  if(I->NState <= state)
    I->NState = state + 1;

  if(I) {
    ObjectCallbackRecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
#endif
}


/*========================================================================*/

void ObjectCallbackRecomputeExtent(ObjectCallback * I)
{
  int extent_flag = false;

#ifndef _PYMOL_NOPY
  float mx[3], mn[3];
  int a;
  PyObject *py_ext;

  for(a = 0; a < I->NState; a++)
    if(I->State[a].PObj) {
      if(PyObject_HasAttrString(I->State[a].PObj, "get_extent")) {
        py_ext = PYOBJECT_CALLMETHOD(I->State[a].PObj, "get_extent", "");
        if(PyErr_Occurred())
          PyErr_Print();
        if(py_ext) {
          if(PConvPyListToExtent(py_ext, mn, mx)) {
            if(!extent_flag) {
              extent_flag = true;
              copy3f(mx, I->Obj.ExtentMax);
              copy3f(mn, I->Obj.ExtentMin);
            } else {
              max3f(mx, I->Obj.ExtentMax, I->Obj.ExtentMax);
              min3f(mn, I->Obj.ExtentMin, I->Obj.ExtentMin);
            }
          }
          Py_DECREF(py_ext);
        }
      }
    }
#endif
  I->Obj.ExtentFlag = extent_flag;

}

/*========================================================================*/

#ifndef _PYMOL_NOPY

static int ObjectCallbackStateFromPyObject(PyMOLGlobals * G, ObjectCallbackState * I,
    PyObject * pobj)
{
  Py_XINCREF(pobj);
  I->PObj = pobj;
  I->is_callable = PyCallable_Check(pobj);

  return true;
}

static int ObjectCallbackAllStatesFromPyObject(ObjectCallback * I, PyObject * obj)
{
  int result = false;
  PyObject *list = NULL;

  if(PyList_Check(obj)) {
    list = obj;
    Py_INCREF(list);
  } else {
    // unpickle list from string
    ok_assert(1, list = PConvPickleLoads(obj));
    ok_assert(1, PyList_Check(list));
  }

  I->NState = PyList_Size(list);
  VLACheck(I->State, ObjectCallbackState, I->NState);

  for(int a = 0; a < I->NState; a++) {
    PyObject *val = PyList_GetItem(list, a);
    ObjectCallbackStateFromPyObject(I->Obj.G, I->State + a, val);
  }

  result = true;
ok_except1:
  if(PyErr_Occurred()) {
    PyErr_Print();

    PRINTFB(I->Obj.G, FB_ObjectCallback, FB_Warnings)
      " Warning: could not load callback object\n"
      ENDFB(I->Obj.G);
  }

  Py_XDECREF(list);
  return result;
}

int ObjectCallbackNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectCallback ** result)
{
  int ll;
  ObjectCallback *I;
  PyObject *val;

  ok_assert(1, list != NULL);
  ok_assert(1, PyList_Check(list));

  ll = PyList_Size(list);

  ok_assert(1, I = ObjectCallbackNew(G));

  val = PyList_GetItem(list, 0);
  ok_assert(2, ObjectFromPyList(G, val, &I->Obj));

  val = PyList_GetItem(list, 1);
  ok_assert(2, ObjectCallbackAllStatesFromPyObject(I, val));

  ObjectCallbackRecomputeExtent(I);

  *result = I;
  return true;
ok_except2:
  ObjectCallbackFree(I);
ok_except1:
  *result = NULL;
  return false;
}

/*========================================================================*/

static PyObject *ObjectCallbackStateAsPyObject(ObjectCallbackState * I)
{
  Py_XINCREF(I->PObj);
  return I->PObj;
}

static PyObject *ObjectCallbackAllStatesAsPyObject(ObjectCallback * I)
{
  int a;
  PyObject *result = NULL;
  PyObject *list = PyList_New(I->NState);

  for(a = 0; a < I->NState; a++) {
    PyList_SetItem(list, a, ObjectCallbackStateAsPyObject(I->State + a));
  }

  // pickle the list to a string
  result = PConvPickleDumps(list);

  Py_XDECREF(list);

  if(PyErr_Occurred()) {
    PyErr_Print();

    PRINTFB(I->Obj.G, FB_ObjectCallback, FB_Warnings)
      " Warning: callable needs to be picklable for session storage\n"
      ENDFB(I->Obj.G);
  }

  return result;
}

PyObject *ObjectCallbackAsPyList(ObjectCallback * I)
{
  PyObject *result = NULL, *states;

  ok_assert(1, states = ObjectCallbackAllStatesAsPyObject(I));

  result = PyList_New(2);
  PyList_SetItem(result, 0, ObjectAsPyList(&I->Obj));
  PyList_SetItem(result, 1, states);

ok_except1:
  return PConvAutoNone(result);
}

#endif
