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
#include"ObjectGadget.h"
#include"ObjectGadgetRamp.h"
#include"GadgetSet.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"VFont.h"

CGO *ObjectGadgetPyListFloatToCGO(PyObject *list);



int ObjectGadgetGetVertex(ObjectGadget *I,int index,int base,float *v)
{
  GadgetSet *gs;
  int ok=false;
  if(I->CurGSet<I->NGSet) {
    gs=I->GSet[I->CurGSet];
    if(gs) {
      ok = GadgetSetGetVertex(gs,index,base,v);
    }
  }
  return(ok);
}

int ObjectGadgetSetVertex(ObjectGadget *I,int index,int base, float *v)
{
  GadgetSet *gs;
  int ok=false;
  if(I->CurGSet<I->NGSet) {
    gs=I->GSet[I->CurGSet];
    if(gs) {
      ok = GadgetSetSetVertex(gs,index,base,v);
    }
  }
  I->Changed=true;
  return(ok);
}

 /* in current state */
ObjectGadget *ObjectGadgetTest(PyMOLGlobals *G)
{
  ObjectGadget *I = NULL;
  GadgetSet *gs = NULL;
  CGO *cgo = NULL;
  int a;

  float coord[] = {
    0.5F ,  0.5F , 0.0F ,
    0.0F ,  0.0F , 0.0F ,
    0.3F ,  0.0F , 0.0F ,
    0.0F , -0.3F , 0.0F ,
    0.3F , -0.3F , 0.0F ,
    0.03F, -0.03F, 0.03F,
    0.27F, -0.03F, 0.03F,
    0.03F, -0.27F, 0.03F,
    0.27F, -0.27F, 0.03F,
    0.02F, -0.02F, 0.01F,
    0.28F, -0.02F, 0.01F,
    0.02F, -0.28F, 0.01F,
    0.28F, -0.28F, 0.01F,
  };

  float normal[] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
   -1.0, 0.0, 0.0,
    0.0,-1.0, 0.0,
  };

  I=ObjectGadgetNew(G);
  gs = GadgetSetNew(G);

  gs->NCoord = 13;
  gs->Coord = VLAlloc(float,gs->NCoord*3);
  for(a=0;a<gs->NCoord*3;a++) {
    gs->Coord[a]=coord[a];
  }

  gs->NNormal = 5;
  gs->Normal = VLAlloc(float,gs->NNormal*3);
  for(a=0;a<gs->NNormal*3;a++) {
    gs->Normal[a]=normal[a];
  }

  cgo = CGONewSized(G,100);
  CGOColor(cgo,1.0,1.0,1.0);

  /* top */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGONormal(cgo,2.0,2.0,0.0);
  CGOVertex(cgo,1.0,5.0,0.0);
  CGOVertex(cgo,1.0,6.0,0.0);

  CGONormal(cgo,2.0,1.0,0.0);
  CGOVertex(cgo,1.0,1.0,0.0);
  CGOVertex(cgo,1.0,2.0,0.0);
  CGOEnd(cgo);

  /* bottom */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGONormal(cgo,2.0,4.0,0.0);
  CGOVertex(cgo,1.0,3.0,0.0);
  CGOVertex(cgo,1.0,4.0,0.0);

  CGONormal(cgo,2.0,2.0,0.0);
  CGOVertex(cgo,1.0,7.0,0.0);
  CGOVertex(cgo,1.0,8.0,0.0);
  CGOEnd(cgo);

  /* left */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGONormal(cgo,2.0,3.0,0.0);
  CGOVertex(cgo,1.0,1.0,0.0);
  CGOVertex(cgo,1.0,3.0,0.0);

  CGONormal(cgo,2.0,2.0,0.0);
  CGOVertex(cgo,1.0,5.0,0.0);
  CGOVertex(cgo,1.0,7.0,0.0);
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGONormal(cgo,2.0,2.0,0.0);
  CGOVertex(cgo,1.0,6.0,0.0);
  CGOVertex(cgo,1.0,8.0,0.0);

  CGONormal(cgo,2.0,0.0,0.0);
  CGOVertex(cgo,1.0,2.0,0.0);
  CGOVertex(cgo,1.0,4.0,0.0);
  CGOEnd(cgo);

  CGOColor(cgo,1.0,0.0,0.0);

  /* center */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGONormal(cgo,2.0,2.0,0.0);
  CGOVertex(cgo,1.0,5.0,0.0);
  CGOVertex(cgo,1.0,7.0,0.0);
  CGOVertex(cgo,1.0,6.0,0.0);
  CGOVertex(cgo,1.0,8.0,0.0);
  CGOEnd(cgo);

  CGOColor(cgo,0.0,1.0,0.0);
  /* backr */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGONormal(cgo,2.0,2.0,0.0);
  CGOVertex(cgo,1.0,9.0,0.0);
  CGOVertex(cgo,1.0,10.0,0.0);
  CGOVertex(cgo,1.0,11.0,0.0);
  CGOVertex(cgo,1.0,12.0,0.0);
  CGOEnd(cgo);



  CGOStop(cgo);

  gs->ShapeCGO = cgo;
  
  cgo = CGONewSized(G,100);
  CGODotwidth(cgo,5);

  CGOPickColor(cgo,0,cPickableGadget);

  /* top */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGOVertex(cgo,1.0,1.0,0.0);
  CGOVertex(cgo,1.0,2.0,0.0);
  CGOVertex(cgo,1.0,5.0,0.0);
  CGOVertex(cgo,1.0,6.0,0.0);
  CGOEnd(cgo);

  /* bottom */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGOVertex(cgo,1.0,3.0,0.0);
  CGOVertex(cgo,1.0,4.0,0.0);
  CGOVertex(cgo,1.0,7.0,0.0);
  CGOVertex(cgo,1.0,8.0,0.0);
  CGOEnd(cgo);

  /* left */

  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGOVertex(cgo,1.0,1.0,0.0);
  CGOVertex(cgo,1.0,3.0,0.0);
  CGOVertex(cgo,1.0,5.0,0.0);
  CGOVertex(cgo,1.0,7.0,0.0);
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo,GL_TRIANGLE_STRIP);
  CGOVertex(cgo,1.0,6.0,0.0);
  CGOVertex(cgo,1.0,8.0,0.0);
  CGOVertex(cgo,1.0,2.0,0.0);
  CGOVertex(cgo,1.0,4.0,0.0);
  CGOEnd(cgo);
  
  CGOEnd(cgo);
  CGOStop(cgo);
  gs->PickShapeCGO = cgo;

  gs->Obj = I;
  gs->State = 0;

  I->GSet[0] = gs;
  I->NGSet = 1;
  I->Obj.Context=1;
  gs->fUpdate(gs);
  ObjectGadgetUpdateExtents(I);
  return(I);
  
}

void ObjectGadgetUpdateExtents(ObjectGadget *I)
{
  float maxv[3] = {FLT_MAX,FLT_MAX,FLT_MAX};
  float minv[3] = {-FLT_MAX,-FLT_MAX,-FLT_MAX};
  int a;
  GadgetSet *ds;

  /* update extents */
  copy3f(maxv,I->Obj.ExtentMin);
  copy3f(minv,I->Obj.ExtentMax);
  I->Obj.ExtentFlag=false;
  for(a=0;a<I->NGSet;a++) {
    ds = I->GSet[a];
    if(ds) {
      if(GadgetSetGetExtent(ds,I->Obj.ExtentMin,I->Obj.ExtentMax))
        I->Obj.ExtentFlag=true;
    }
  }
}

#ifndef _PYMOL_NOPY
static PyObject *ObjectGadgetGSetAsPyList(ObjectGadget *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NGSet);
  for(a=0;a<I->NGSet;a++) {
    if(I->GSet[a]) {
      PyList_SetItem(result,a,GadgetSetAsPyList(I->GSet[a]));
    } else {
      PyList_SetItem(result,a,PConvAutoNone(Py_None));
    }
  }
  return(PConvAutoNone(result));

}
#endif
#ifndef _PYMOL_NOPY
static int ObjectGadgetGSetFromPyList(ObjectGadget *I,PyObject *list,int version)
{

  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->GSet,GadgetSet*,I->NGSet);
    for(a=0;a<I->NGSet;a++) {
      if(ok) ok = GadgetSetFromPyList(I->Obj.G,PyList_GetItem(list,a),&I->GSet[a],version);
      if(ok&&I->GSet[a]) {
        I->GSet[a]->Obj = I;
        I->GSet[a]->State = a;
      }
    }
  }
  return(ok);
}
#endif
int ObjectGadgetInitFromPyList(PyMOLGlobals *G,PyObject *list,ObjectGadget *I,int version)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  int ll;
  if(ok) ok = (I!=NULL)&&(list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) ok = ObjectFromPyList(G,PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->GadgetType);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->NGSet);
  if(ok) ok = ObjectGadgetGSetFromPyList(I,PyList_GetItem(list,3),version);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,4),&I->CurGSet);
  
  /*  ObjectGadgetInvalidateRep(I,cRepAll);*/
  if(ok) {
    ObjectGadgetUpdateExtents(I);
  }
  else {
    /* cleanup? */
  }
  return(ok);
#endif
}

int ObjectGadgetNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectGadget **result,int version)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  ObjectGadget *I=NULL;
  int gadget_type = -1;
  PyObject *plain;
  (*result) = NULL;

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);

  /* NOTE there is a serious screw-up here...ramp gadgets aren't saved right, but
     we've got to maintain backward compat...ugh */

  if(ok) ok=((plain=PyList_GetItem(list,0))!=NULL); 
  if(ok) ok=PyList_Check(plain);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(plain,1),&gadget_type);
  if(ok) switch(gadget_type) { /* call the right routine to restore the gadget! */
  case cGadgetRamp:
    ok = ObjectGadgetRampNewFromPyList(G,list,(ObjectGadgetRamp**)result,version);
    break;
  case cGadgetPlain:
    I=ObjectGadgetNew(G);
    if(ok) ok = (I!=NULL);
    if(ok) ok = ObjectGadgetInitFromPyList(G,list,I,version);
    if(ok) (*result) = I;
    break;
  default:
    ok=false;
    break;
  } 
  return(ok);
#endif
}

PyObject *ObjectGadgetPlainAsPyList(ObjectGadget *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  /* first, dump the atoms */

  result = PyList_New(5);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->GadgetType));
  PyList_SetItem(result,2,PyInt_FromLong(I->NGSet));
  PyList_SetItem(result,3,ObjectGadgetGSetAsPyList(I));
  PyList_SetItem(result,4,PyInt_FromLong(I->CurGSet));
  return(PConvAutoNone(result));  
#endif
}

PyObject *ObjectGadgetAsPyList(ObjectGadget *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  /* first, dump the atoms */

  switch(I->GadgetType) { 
  case cGadgetRamp:
    result = ObjectGadgetRampAsPyList((ObjectGadgetRamp*)I);
    break;
  case cGadgetPlain:
    result = ObjectGadgetPlainAsPyList(I);
    break;
  } 
  return(PConvAutoNone(result));  
#endif
}

void ObjectGadgetPurge(ObjectGadget *I)
{
  int a;

  SceneObjectDel(I->Obj.G,(CObject*)I);
  for(a=0;a<I->NGSet;a++)
	 if(I->GSet[a]) {
      if(I->GSet[a]->fFree)
        I->GSet[a]->fFree(I->GSet[a]);
		I->GSet[a]=NULL;
	 }
  VLAFreeP(I->GSet);
  ObjectPurge(&I->Obj);
}

void ObjectGadgetFree(ObjectGadget *I) {
  ObjectGadgetPurge(I);
  OOFreeP(I);
}

void ObjectGadgetUpdateStates(ObjectGadget *I)
{
  int a;
  OrthoBusyPrime(I->Obj.G);
  for(a=0;a<I->NGSet;a++)
    if(I->GSet[a]) {	
      OrthoBusySlow(I->Obj.G,a,I->NGSet);
      /*	   printf(" ObjectGadget: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);*/
      if(I->GSet[a]->fUpdate)
        I->GSet[a]->fUpdate(I->GSet[a]);
    }
}

/*========================================================================*/
void ObjectGadgetUpdate(ObjectGadget *I)
{
  if(I->Changed) {
    ObjectGadgetUpdateStates(I);
    ObjectGadgetUpdateExtents(I);
    I->Changed=false;
  }
}

/*========================================================================*/

static int ObjectGadgetGetNState(ObjectGadget *I) {
  return(I->NGSet);
}

/*========================================================================*/
static void ObjectGadgetRender(ObjectGadget *I,RenderInfo *info)
{
  int state = info->state;
  CRay *ray = info->ray;
  int pass = info->pass;
  int a;
  if(!pass) {

    ObjectPrepareContext(&I->Obj,ray);
    if(state<0) {
      for(a=0;a<I->NGSet;a++)
        if(I->GSet[a])
          if(I->GSet[a]->fRender)
            I->GSet[a]->fRender(I->GSet[a],info);
    } else if(state<I->NGSet) {
      I->CurGSet=state;
      if(I->GSet[I->CurGSet]) {
        if(I->GSet[I->CurGSet]->fRender)
          I->GSet[I->CurGSet]->fRender(I->GSet[I->CurGSet],info);
      }
    } else if(I->NGSet==1) { /* if only one coordinate set, assume static */
      if(I->GSet[0]->fRender)
        I->GSet[0]->fRender(I->GSet[0],info);
      I->CurGSet=0;
    }
  }
}

/*========================================================================*/
void ObjectGadgetInit(PyMOLGlobals *G,ObjectGadget *I)
{
  ObjectInit(G,(CObject*)I);

  I->Obj.type=cObjectGadget;
  I->GSet=VLAMalloc(10,sizeof(GadgetSet*),5,true); /* auto-zero */
  I->NGSet=0;
  I->Changed=true;

  I->Obj.fFree = (void (*)(CObject *))ObjectGadgetFree;
  I->Obj.fUpdate =(void (*)(CObject *)) ObjectGadgetUpdate;
  I->Obj.fRender =(void (*)(CObject *, RenderInfo *info))ObjectGadgetRender;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectGadgetGetNState;
  I->Obj.fDescribeElement = NULL;
  I->CurGSet=0;
}

/*========================================================================*/
ObjectGadget *ObjectGadgetNew(PyMOLGlobals *G)
{
  OOAlloc(G,ObjectGadget);

  ObjectGadgetInit(G,I);
  return(I);
}

#if 0

/*========================================================================*/
CGO *ObjectGadgetPyListFloatToCGO(PyObject *list)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

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
        cgo=CGONewSized(G,len);
        if(cgo) {
          result = CGOFromFloatArray(cgo,raw,len);
          if(result) {
            PRINTF " FloatToCGO: error encountered on element %d\n", result ENDF(G);
          }
          CGOStop(cgo);
        }
      }
      FreeP(raw);
    }
  }
  return(cgo);
#endif
}


ObjectGadget *ObjectGadgetFromCGO(PyMOLGlobals *G,ObjectGadget *obj,CGO *cgo,int state)
{
  ObjectGadget *I = NULL;
  int est;

  if(obj) {
    if(obj->Obj.type!=cObjectGadget) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectGadgetNew(G);
  } else {
    I=obj;
  }
  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectGadgetState,state);
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
    ObjectGadgetRecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return(I);
}


/*========================================================================*/
ObjectGadget *ObjectGadgetDefine(PyMOLGlobals *G,ObjectGadget *obj,PyObject *pycgo,int state)
{ /* assumes blocked interpreter */
#ifdef _PYMOL_NOPY
  return NULL;
#else
  ObjectGadget *I = NULL;

  CGO *cgo,*font_cgo;
  int est;

  if(obj) {
    if(obj->Obj.type!=cObjectGadget) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectGadgetNew(G);
  } else {
    I=obj;
  }
  I->GadgetType = cGadgetPlain;

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectGadgetState,state);
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
        cgo=ObjectGadgetPyListFloatToCGO(pycgo);
        if(cgo) {
          est=CGOCheckForText(cgo);
          if(est) {
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
          ErrMessage(G,"ObjectGadget","could not parse CGO List.");
        }
      }
    }
  }
  if(I) {
    ObjectGadgetRecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return(I);
#endif

}
#endif
