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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectCGO.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"VFont.h"

static void ObjectCGOFree(ObjectCGO *I);
CGO *ObjectCGOPyListFloatToCGO(PyObject *list);



static PyObject *ObjectCGOStateAsPyList(ObjectCGOState *I)
{
  PyObject *result = NULL;

  result = PyList_New(2);
  if(I->std) 
    PyList_SetItem(result,0,CGOAsPyList(I->std));
  else
    PyList_SetItem(result,0,PConvAutoNone(NULL));
  if(I->ray) 
    PyList_SetItem(result,1,CGOAsPyList(I->ray));
  else
    PyList_SetItem(result,1,PConvAutoNone(NULL));
  return(PConvAutoNone(result));  
}


static PyObject *ObjectCGOAllStatesAsPyList(ObjectCGO *I)
{
  
  PyObject *result=NULL;
  int a;
  result = PyList_New(I->NState);
  for(a=0;a<I->NState;a++) {
    PyList_SetItem(result,a,ObjectCGOStateAsPyList(I->State+a));
  }
  return(PConvAutoNone(result));  

}

static int ObjectCGOStateFromPyList(ObjectCGOState *I,PyObject *list,int version)
{
  int ok=true;
  int ll;
  PyObject *tmp;
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) {
    tmp = PyList_GetItem(list,0);
    if(tmp == Py_None)
      I->std = NULL;
    else 
      ok = ((I->std=CGONewFromPyList(PyList_GetItem(list,0),version))!=NULL);
  }
  if(ok) {
    tmp = PyList_GetItem(list,1);
    if(tmp == Py_None)
      I->ray = NULL;
    else 
      ok = ((I->ray=CGONewFromPyList(PyList_GetItem(list,1),version))!=NULL);
  }
  return(ok);
}

static int ObjectCGOAllStatesFromPyList(ObjectCGO *I,PyObject *list,int version)
{
  int ok=true;
  int a;
  VLACheck(I->State,ObjectCGOState,I->NState);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    for(a=0;a<I->NState;a++) {
      ok = ObjectCGOStateFromPyList(I->State+a,PyList_GetItem(list,a),version);
      if(!ok) break;
    }
  }
  return(ok);
}

int ObjectCGONewFromPyList(PyObject *list,ObjectCGO **result,int version)
{
  int ok = true;
  ObjectCGO *I=NULL;
  (*result) = NULL;
  if(ok) ok=(list!=Py_None);
  if(ok) ok=PyList_Check(list);

  I=ObjectCGONew();
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NState);
  if(ok) ok = ObjectCGOAllStatesFromPyList(I,PyList_GetItem(list,2),version);
  if(ok) {
    (*result) = I;
    ObjectCGORecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return(ok);
}




PyObject *ObjectCGOAsPyList(ObjectCGO *I)
{
  
  PyObject *result=NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NState));
  PyList_SetItem(result,2,ObjectCGOAllStatesAsPyList(I));

  return(PConvAutoNone(result));  
}


/*========================================================================*/

static void ObjectCGOFree(ObjectCGO *I) {
  int a;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].std)
      CGOFree(I->State[a].std);
    if(I->State[a].ray)
      CGOFree(I->State[a].ray);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}

/*========================================================================*/

void ObjectCGORecomputeExtent(ObjectCGO *I)
{
  float mx[3],mn[3];
  int extent_flag = false;
  int a;
  for(a=0;a<I->NState;a++) 
    if(I->State[a].std) {
      if(CGOGetExtent(I->State[a].std,mn,mx)) {
        if(!extent_flag) {
          extent_flag=true;
          copy3f(mx,I->Obj.ExtentMax);
          copy3f(mn,I->Obj.ExtentMin);
        } else {
          max3f(mx,I->Obj.ExtentMax,I->Obj.ExtentMax);
          min3f(mn,I->Obj.ExtentMin,I->Obj.ExtentMin);
        }
      }
    }
  I->Obj.ExtentFlag=extent_flag;
}
/*========================================================================*/
static void ObjectCGOUpdate(ObjectCGO *I)
{
  SceneDirty();/* needed ?*/
}

/*========================================================================*/

static int ObjectCGOGetNState(ObjectCGO *I) {
  return(I->NState);
}

/*========================================================================*/

static void ObjectCGORender(ObjectCGO *I,int state,CRay *ray,Pickable **pick,int pass)
{
  ObjectCGOState *sobj = NULL;
  int a;
  float *color;

  ObjectPrepareContext(&I->Obj,ray);

  color = ColorGet(I->Obj.Color);

  if(!pass) {
    if(I->Obj.RepVis[cRepCGO]) {
      
      if(state<I->NState) {
        sobj = I->State+state;
      }
      if(state<0) {
        if(I->State) {
          for(a=0;a<I->NState;a++) {
            sobj = I->State+a;
            if(ray) {    
              if(sobj->ray)
                CGORenderRay(sobj->ray,ray,color,I->Obj.Setting,NULL);
              else
                CGORenderRay(sobj->std,ray,color,I->Obj.Setting,NULL);
            } else if(pick&&PMGUI) {
            } else if(PMGUI) {
              if(sobj->std)
                CGORenderGL(sobj->std,color,I->Obj.Setting,NULL);
            }
          }
        }
      } else {
        if(!sobj) {
          if(I->NState&&SettingGet(cSetting_static_singletons)) 
            sobj = I->State;
        }
        if(ray) {    
          if(sobj)
            {
              if(sobj->ray)
                CGORenderRay(sobj->ray,ray,color,I->Obj.Setting,NULL);
              else
                CGORenderRay(sobj->std,ray,color,I->Obj.Setting,NULL);
            }
        } else if(pick&&PMGUI) {
        } else if(PMGUI) {
          if(sobj)
            if(sobj->std)
              CGORenderGL(sobj->std,color,I->Obj.Setting,NULL);
        }
      }
    }
  }
}

/*========================================================================*/
ObjectCGO *ObjectCGONew(void)
{
  OOAlloc(ObjectCGO);

  ObjectInit((CObject*)I);

  I->State=VLAMalloc(10,sizeof(ObjectCGOState),5,true);
  I->NState=0;

  I->Obj.type = cObjectCGO;
  I->Obj.fFree = (void (*)(struct CObject *))ObjectCGOFree;
  I->Obj.fUpdate =(void (*)(struct CObject *)) ObjectCGOUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, int, CRay *, Pickable **,int))ObjectCGORender;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectCGOGetNState;

  return(I);
}

/*========================================================================*/
CGO *ObjectCGOPyListFloatToCGO(PyObject *list)
{
  CGO *cgo=NULL;
  int len;
  int ok = true;
  int result;
  float *raw=NULL;
  if(PyList_Check(list)) {
    len = PConvPyListToFloatArray(list,&raw);
    if(len<0) len = 0;
    if(raw) {
      if(ok) {
        cgo=CGONewSized(len);
        if(cgo) {
          result = CGOFromFloatArray(cgo,raw,len);
          if(result) {
            PRINTF " FloatToCGO: error encountered on element %d\n", result ENDF
          }
          CGOStop(cgo);
        }
      }
      FreeP(raw);
    }
  }
  return(cgo);
}
/*========================================================================*/
ObjectCGO *ObjectCGOFromCGO(ObjectCGO *obj,CGO *cgo,int state)
{
  ObjectCGO *I = NULL;
  int est;

  if(obj) {
    if(obj->Obj.type!=cObjectCGO) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectCGONew();
  } else {
    I=obj;
  }
  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectCGOState,state);
    I->NState=state+1;
  }

  if(I->State[state].std) {
    CGOFree(I->State[state].std);
  }
  if(I->State[state].ray) {
    CGOFree(I->State[state].ray);
  }
  est=CGOCheckComplex(cgo);
  if(est) {
    I->State[state].ray=cgo;
    I->State[state].std=CGOSimplify(cgo,est);
  } else 
    I->State[state].std=cgo;
  if(I) {
    ObjectCGORecomputeExtent(I);
  }
  SceneChanged();
  SceneCountFrames();
  return(I);
}
/*========================================================================*/

ObjectCGO *ObjectCGONewVFontTest(char *text,float *pos)
{
  ObjectCGO *I = NULL;
  int font_id;
  CGO *cgo = NULL;
  float scale[2] = {1.0,1.0};

  font_id = VFontLoad(1,1,1,true);
  cgo = CGONew();
  VFontWriteToCGO(font_id,cgo,text,pos,scale,NULL);
  I = ObjectCGOFromCGO(NULL,cgo,0);

  return(I);
}


/*========================================================================*/
ObjectCGO *ObjectCGODefine(ObjectCGO *obj,PyObject *pycgo,int state)
{ /* assumes blocked interpreter */
  ObjectCGO *I = NULL;

  CGO *cgo,*font_cgo;
  int est;

  if(obj) {
    if(obj->Obj.type!=cObjectCGO) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectCGONew();
  } else {
    I=obj;
  }
  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectCGOState,state);
    I->NState=state+1;
  }

  if(I->State[state].std) {
    CGOFree(I->State[state].std);
  }
  if(I->State[state].ray) {
    CGOFree(I->State[state].ray);
  }
  if(PyList_Check(pycgo)) {
    if(PyList_Size(pycgo)) {
      if(PyFloat_Check(PyList_GetItem(pycgo,0))) {
        cgo=ObjectCGOPyListFloatToCGO(pycgo);
        if(cgo) {
          est=CGOCheckForText(cgo);
          if(est) {
            CGOPreloadFonts(cgo);
            font_cgo = CGODrawText(cgo,est,NULL);
            CGOFree(cgo);
            cgo=font_cgo;
          }
          est=CGOCheckComplex(cgo);
          if(est) {
            I->State[state].ray=cgo;
            I->State[state].std=CGOSimplify(cgo,est);
          } else 
            I->State[state].std=cgo;
          
        } else {
          ErrMessage("ObjectCGO","could not parse CGO List.");
        }
      }
    }
  }
  if(I) {
    ObjectCGORecomputeExtent(I);
  }
  SceneChanged();
  SceneCountFrames();
  return(I);
}
