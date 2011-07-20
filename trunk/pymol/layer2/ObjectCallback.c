
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
  PBlock(G);
  for(a = 0; a < I->NState; a++) {
    if(I->State[a].PObj) {
      Py_DECREF(I->State[a].PObj);
      I->State[a].PObj = NULL;
    }
  }
  PUnblock(G);
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
  int a;

  PyObject *pobj = NULL;

  if(!pass) {

    ObjectPrepareContext(&I->Obj, ray);
    if(I->Obj.RepVis[cRepCallback]) {
      if(state < I->NState) {
        sobj = I->State + state;
      }
      if(state < 0) {
        if(I->State) {
          PBlock(G);
          for(a = 0; a < I->NState; a++) {
            sobj = I->State + a;
            pobj = sobj->PObj;
            if(ray) {
            } else if(G->HaveGUI && G->ValidContext) {
              if(pick) {
              } else {
                if(PyObject_HasAttrString(pobj, "__call__")) {
                  Py_XDECREF(PyObject_CallMethod(pobj, "__call__", ""));
                }
                if(PyErr_Occurred())
                  PyErr_Print();
              }
            }
          }
          PUnblock(G);
        }
      } else {
        if(!sobj) {
          if(I->NState && SettingGet(G, cSetting_static_singletons))
            sobj = I->State;
        }
        if(ray) {
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
          } else {
            if(sobj) {
              pobj = sobj->PObj;
              PBlock(G);
              if(PyObject_HasAttrString(pobj, "__call__")) {
                Py_XDECREF(PyObject_CallMethod(pobj, "__call__", ""));
              }
              if(PyErr_Occurred())
                PyErr_Print();
              PUnblock(G);
            }
          }
        }
      }
    }
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

  I->State = VLAMalloc(10, sizeof(ObjectCallbackState), 5, true);       /* autozero */
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
        py_ext = PyObject_CallMethod(I->State[a].PObj, "get_extent", "");
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
