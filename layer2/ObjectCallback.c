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


static void ObjectCallbackFree(ObjectCallback *I);

/*========================================================================*/

static void ObjectCallbackFree(ObjectCallback *I) {
  if(I->PObj) {
    Py_DECREF(I->PObj);
  }
  OOFreeP(I);
}

/*========================================================================*/

static void ObjectCallbackUpdate(ObjectCallback *I) {
  SceneDirty();
}

/*========================================================================*/

static void ObjectCallbackRender(ObjectCallback *I,int frame,CRay *ray,Pickable **pick)
{
  if(ray) {    
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
    if(I->PObj) {
      PBlock();
      if(PyObject_HasAttrString(I->PObj,"__call__"))
         PyObject_CallMethod(I->PObj,"__call__","");
      if(PyErr_Occurred())
        PyErr_Print();
      PUnblock();
    }
  }
}

/*========================================================================*/
ObjectCallback *ObjectCallbackNew(void)
{
  OOAlloc(ObjectCallback);

  ObjectInit((Object*)I);

  I->PObj = NULL;

  I->Obj.type = cObjectCallback;
  I->Obj.fFree = (void (*)(struct Object *))ObjectCallbackFree;
  I->Obj.fUpdate =  (void (*)(struct Object *)) ObjectCallbackUpdate;
  I->Obj.fRender =(void (*)(struct Object *, int, CRay *, Pickable **))ObjectCallbackRender;

#ifdef _NOT_YET_NEEDED
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectCallbackGetNFrames;
#endif

  return(I);
}
/*========================================================================*/
ObjectCallback *ObjectCallbackDefine(ObjectCallback *obj,PyObject *pobj)
{
  ObjectCallback *I = NULL;
  I=ObjectCallbackNew();
  I->PObj = pobj; /* need to set extents */
  Py_INCREF(I->PObj);
  SceneChanged();
  return(I);
}
