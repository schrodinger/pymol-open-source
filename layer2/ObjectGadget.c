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

void ObjectGadgetRender(ObjectGadget *I,int state,CRay *ray,Pickable **pick,int pass);

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
  return(ok);
}

 /* in current state */
ObjectGadget *ObjectGadgetTest(void)
{
  ObjectGadget *I = NULL;
  GadgetSet *gs = NULL;
  CGO *cgo = NULL;
  int a;

  float coord[] = {
    0.5 ,  0.5 , 0.0 ,
    0.0 ,  0.0 , 0.0 ,
    0.3 ,  0.0 , 0.0 ,
    0.0 , -0.3 , 0.0 ,
    0.3 , -0.3 , 0.0 ,
    0.03, -0.03, 0.03,
    0.27, -0.03, 0.03,
    0.03, -0.27, 0.03,
    0.27, -0.27, 0.03,
    0.02, -0.02, 0.01,
    0.28, -0.02, 0.01,
    0.02, -0.28, 0.01,
    0.28, -0.28, 0.01,
  };

  float normal[] = {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
   -1.0, 0.0, 0.0,
    0.0,-1.0, 0.0,
  };

  I=ObjectGadgetNew();
  gs = GadgetSetNew();

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

  cgo = CGONewSized(100);
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
  
  cgo = CGONewSized(100);
  CGODotwidth(cgo,5);

  CGOPickColor(cgo,0,0);

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

static PyObject *ObjectGadgetGSetAsPyList(ObjectGadget *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NGSet);
  for(a=0;a<I->NGSet;a++) {
    if(I->GSet[a]) {
      PyList_SetItem(result,a,GadgetSetAsPyList(I->GSet[a]));
    } else {
      PyList_SetItem(result,a,Py_None);
      Py_INCREF(Py_None);
    }
  }
  return(PConvAutoNone(result));
}

static int ObjectGadgetGSetFromPyList(ObjectGadget *I,PyObject *list)
{
  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->GSet,GadgetSet*,I->NGSet);
    for(a=0;a<I->NGSet;a++) {
      if(ok) ok = GadgetSetFromPyList(PyList_GetItem(list,a),&I->GSet[a]);
      /*      if(ok) I->GSet[a]->Obj = I;*/
    }
  }
  return(ok);
}

int ObjectGadgetNewFromPyList(PyObject *list,ObjectGadget **result)
{
  int ok = true;
  ObjectGadget *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=PyList_Check(list);

  I=ObjectGadgetNew();
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NGSet);
  if(ok) ok = ObjectGadgetGSetFromPyList(I,PyList_GetItem(list,2));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,3),&I->CurGSet);
  
  /*  ObjectGadgetInvalidateRep(I,cRepAll);*/
  if(ok) {
    (*result) = I;
    ObjectGadgetUpdateExtents(I);
  }
  else {
    /* cleanup? */
  }

  return(ok);
}




PyObject *ObjectGadgetAsPyList(ObjectGadget *I)
{
  PyObject *result = NULL;

  /* first, dump the atoms */

  result = PyList_New(4);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NGSet));
  PyList_SetItem(result,2,ObjectGadgetGSetAsPyList(I));
  PyList_SetItem(result,3,PyInt_FromLong(I->CurGSet));

#if 0

  CObject Obj;
  struct GadgetSet **GSet;
  int NGSet;
  int CurGSet;

#endif

  return(PConvAutoNone(result));  
}

void ObjectGadgetPurge(ObjectGadget *I)
{
  int a;

  SceneObjectDel((CObject*)I);
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
  OrthoBusyPrime();
  for(a=0;a<I->NGSet;a++)
    if(I->GSet[a]) {	
      OrthoBusySlow(a,I->NGSet);
      /*	   printf(" ObjectGadget: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);*/
      if(I->GSet[a]->fUpdate)
        I->GSet[a]->fUpdate(I->GSet[a]);
    }
}

/*========================================================================*/
void ObjectGadgetUpdate(ObjectGadget *I)
{
  ObjectGadgetUpdateStates(I);
  ObjectGadgetUpdateExtents(I);
}

/*========================================================================*/

static int ObjectGadgetGetNState(ObjectGadget *I) {
  return(I->NGSet);
}

/*========================================================================*/
void ObjectGadgetRender(ObjectGadget *I,int state,CRay *ray,Pickable **pick,int pass)
{
  int a;
  if(!pass) {

    ObjectPrepareContext(&I->Obj,ray);
    if(state<0) {
      for(a=0;a<I->NGSet;a++)
        if(I->GSet[a])
          if(I->GSet[a]->fRender)
            I->GSet[a]->fRender(I->GSet[a],ray,pick,pass);        
    } else if(state<I->NGSet) {
      I->CurGSet=state;
      if(I->GSet[I->CurGSet]) {
        if(I->GSet[I->CurGSet]->fRender)
          I->GSet[I->CurGSet]->fRender(I->GSet[I->CurGSet],ray,pick,pass);
      }
    } else if(I->NGSet==1) { /* if only one coordinate set, assume static */
      if(I->GSet[0]->fRender)
        I->GSet[0]->fRender(I->GSet[0],ray,pick,pass);    
      I->CurGSet=0;
    }
  }
}

/*========================================================================*/
void ObjectGadgetInit(ObjectGadget *I)
{
  ObjectInit((CObject*)I);

  I->Obj.type=cObjectGadget;
  I->GSet=VLAMalloc(10,sizeof(GadgetSet*),5,true); /* auto-zero */
  I->NGSet=0;

  I->Obj.type = cObjectGadget;
  I->Obj.fFree = (void (*)(struct CObject *))ObjectGadgetFree;
  I->Obj.fUpdate =(void (*)(struct CObject *)) ObjectGadgetUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, int, CRay *, Pickable **,int))ObjectGadgetRender;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectGadgetGetNState;
  I->Obj.fDescribeElement = NULL;
  I->CurGSet=0;
}

/*========================================================================*/
ObjectGadget *ObjectGadgetNew(void)
{
  OOAlloc(ObjectGadget);

  ObjectGadgetInit(I);
  return(I);
}

#if 0

/*========================================================================*/
CGO *ObjectGadgetPyListFloatToCGO(PyObject *list)
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


ObjectGadget *ObjectGadgetFromCGO(ObjectGadget *obj,CGO *cgo,int state)
{
  ObjectGadget *I = NULL;
  int est;

  if(obj) {
    if(obj->Obj.type!=cObjectGadget) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectGadgetNew();
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
  SceneChanged();
  SceneCountFrames();
  return(I);
}


/*========================================================================*/
ObjectGadget *ObjectGadgetDefine(ObjectGadget *obj,PyObject *pycgo,int state)
{ /* assumes blocked interpreter */
  ObjectGadget *I = NULL;

  CGO *cgo,*font_cgo;
  int est;

  if(obj) {
    if(obj->Obj.type!=cObjectGadget) /* TODO: handle this */
      obj=NULL;
  }
  if(!obj) {
    I=ObjectGadgetNew();
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
          ErrMessage("ObjectGadget","could not parse CGO List.");
        }
      }
    }
  }
  if(I) {
    ObjectGadgetRecomputeExtent(I);
  }
  SceneChanged();
  SceneCountFrames();
  return(I);
}
#endif
