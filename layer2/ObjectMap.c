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
#include"ObjectMap.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Parse.h"
#include"Isosurf.h"
#include"Vector.h"
#include"Color.h"
#include"main.h"
#include"Scene.h"
#include"PConv.h"
#include"Word.h"
#include"Vector.h"

#define cMapSourceUndefined 0
#define cMapSourceXPLOR 1
#define cMapSourceCCP4 2
#define cMapSourcePHI 3
#define cMapSourceDesc 4

#ifdef _PYMOL_NUMPY
typedef struct {
  PyObject_HEAD
  char *data;
  int nd;
  int *dimensions, *strides;
  PyObject *base;
  void *descr;
  int flags;
} MyArrayObject;
#endif


int ObjectMapInterpolate(ObjectMap *I,int state,float *array,float *result,int n)
{
  int ok=false;
  if(state<0) state=0;
  if(state<I->NState)
    if(I->State[state].Active)
      ok = ObjectMapStateInterpolate(&I->State[state],array,result,n);
  return(ok);
}

int ObjectMapStateInterpolate(ObjectMapState *ms,float *array,float *result,int n)
{
  int ok=true;
  float *inp;
  float *out;
  int a,b,c;
  float x,y,z;
  inp = array;
  out = result;
  
  switch(ms->MapSource) {
  case cMapSourcePHI:
  case cMapSourceDesc:
    while(n--) {

      x = (inp[0] - ms->Origin[0])/ms->Grid[0];
      y = (inp[1] - ms->Origin[1])/ms->Grid[1];
      z = (inp[2] - ms->Origin[2])/ms->Grid[2];
      inp+=3;
      
      a=(int)floor(x);
      b=(int)floor(y);
      c=(int)floor(z);
      x-=a;
      y-=b;
      z-=c;

      if(a<ms->Min[0]) {
        x=0.0F;
        a=ms->Min[0];
        ok=false;
      } else if(a>=ms->Max[0]) {
        x=1.0F;
        a=ms->Max[0]-1;
        ok=false;
      }

      if(b<ms->Min[1]) {
        y=0.0F;
        b=ms->Min[1];
        ok=false;
      } else if(b>=ms->Max[1]) {
        y=1.0F;
        b=ms->Min[1];
        ok=false;
      }

      if(c<ms->Min[2]) {
        z=0.0F;
        c=ms->Min[2];
        ok=false;
      } else if(c>=ms->Max[2]) {
        z=1.0F;
        c=ms->Max[2]-1;
        ok=false;
      }
      /*      printf("%d %d %d %8.3f %8.3f %8.3f\n",a,b,c,x,y,z);*/
      *(result++)=FieldInterpolatef(ms->Field->data,
                                    a-ms->Min[0],
                                    b-ms->Min[1],
                                    c-ms->Min[2],x,y,z);
    }
    break;
  default:
    ok=false;
    break;
  }
  return(ok);
}

int ObjectMapNumPyArrayToMapState(ObjectMapState *I,PyObject *ary);

static void ObjectMapStateRegeneratePoints(ObjectMapState *ms)
{
  int a,b,c,e;
  float v[3],vr[3];
  switch(ms->MapSource) {
  case cMapSourceXPLOR:
  case cMapSourceCCP4:
    for(c=0;c<ms->FDim[2];c++)
      {
        v[2]=(c+ms->Min[2])/((float)ms->Div[2]);
        for(b=0;b<ms->FDim[1];b++) {
          v[1]=(b+ms->Min[1])/((float)ms->Div[1]);
          for(a=0;a<ms->FDim[0];a++) {
            v[0]=(a+ms->Min[0])/((float)ms->Div[0]);
            transform33f3f(ms->Crystal->FracToReal,v,vr);
            for(e=0;e<3;e++) 
              F4(ms->Field->points,a,b,c,e) = vr[e];
          }
        }
      }
    break;
  case cMapSourcePHI:
    for(c=0;c<ms->FDim[2];c++) {
      v[2]=ms->Origin[2]+ms->Grid[2]*(c+ms->Min[2]);
      for(b=0;b<ms->FDim[1];b++) {
        v[1]=ms->Origin[1]+ms->Grid[1]*(b+ms->Min[1]);
        for(a=0;a<ms->FDim[0];a++) {
          v[0]=ms->Origin[0]+ms->Grid[0]*(a+ms->Min[0]);
          for(e=0;e<3;e++) {
            F4(ms->Field->points,a,b,c,e) = v[e];
          }
        }
      }
  }
  default: 
    break;
  }
}

static PyObject *ObjectMapStateAsPyList(ObjectMapState *I)
{
  PyObject *result = NULL;

  result = PyList_New(15);
  PyList_SetItem(result,0,PyInt_FromLong(I->Active));
  if(I->Crystal) {
    PyList_SetItem(result,1,CrystalAsPyList(I->Crystal));
  } else {
    PyList_SetItem(result,1,PConvAutoNone(Py_None));
  }
  if(I->Origin) {
    PyList_SetItem(result,2,PConvFloatArrayToPyList(I->Origin,3));
  } else {
    PyList_SetItem(result,2,PConvAutoNone(Py_None));
  }
  if(I->Range) {
    PyList_SetItem(result,3,PConvFloatArrayToPyList(I->Range,3));
  } else {
    PyList_SetItem(result,3,PConvAutoNone(Py_None));
  }
  if(I->Dim) {
    PyList_SetItem(result,4,PConvIntArrayToPyList(I->Dim,3));
  } else {
    PyList_SetItem(result,4,PConvAutoNone(Py_None));
  }
  if(I->Grid) {
    PyList_SetItem(result,5,PConvFloatArrayToPyList(I->Grid,3));
  } else {
    PyList_SetItem(result,5,PConvAutoNone(Py_None));
  }
  PyList_SetItem(result,6,PConvFloatArrayToPyList(&I->Corner[0][0],24));
  PyList_SetItem(result,7,PConvFloatArrayToPyList(I->ExtentMin,3));
  PyList_SetItem(result,8,PConvFloatArrayToPyList(I->ExtentMax,3));
  PyList_SetItem(result,9,PyInt_FromLong(I->MapSource));

  PyList_SetItem(result,10,PConvIntArrayToPyList(I->Div,3));
  PyList_SetItem(result,11,PConvIntArrayToPyList(I->Min,3));
  PyList_SetItem(result,12,PConvIntArrayToPyList(I->Max,3));
  PyList_SetItem(result,13,PConvIntArrayToPyList(I->FDim,4));
  
  PyList_SetItem(result,14,IsosurfAsPyList(I->Field));
#if 0
  int Active;
  CCrystal *Crystal;
  int Div[3],Min[3],Max[3],FDim[4];
  Isofield *Field;
  float Corner[8][3];
  int *Dim;
  float *Origin;
  float *Range;
  float *Grid;
  float ExtentMin[3],ExtentMax[3];
#endif

  return(PConvAutoNone(result));  
}

static PyObject *ObjectMapAllStatesAsPyList(ObjectMap *I)
{
  PyObject *result=NULL;
  int a;
  result = PyList_New(I->NState);
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result,a,ObjectMapStateAsPyList(I->State+a));
    } else {
      PyList_SetItem(result,a,PConvAutoNone(NULL));
    }
  }
  return(PConvAutoNone(result));  
}

static int ObjectMapStateFromPyList(ObjectMapState *I,PyObject *list)
{
  int ok=true;
  int ll;
  PyObject *tmp;
  if(ok) ok=(list!=NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active=false;
    else {
      if(ok) ok=PyList_Check(list);
      if(ok) ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */

      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->Active);
      if(ok) {
        tmp = PyList_GetItem(list,1);
        if(tmp == Py_None)
          I->Crystal = NULL;
        else 
          ok = ((I->Crystal=CrystalNewFromPyList(tmp))!=NULL);
      }
      if(ok) {
        tmp = PyList_GetItem(list,2);
        if(tmp == Py_None)
          I->Origin = NULL;
        else 
          ok = PConvPyListToFloatArray(tmp,&I->Origin);
      }
      if(ok) {
        tmp = PyList_GetItem(list,3);
        if(tmp == Py_None)
          I->Range = NULL;
        else 
          ok = PConvPyListToFloatArray(tmp,&I->Range);
      }
      if(ok) {
        tmp = PyList_GetItem(list,4);
        if(tmp == Py_None)
          I->Dim = NULL;
        else 
          ok = PConvPyListToIntArray(tmp,&I->Dim);
      }
      if(ok) {
        tmp = PyList_GetItem(list,5);
        if(tmp == Py_None)
          I->Grid = NULL;
        else 
          ok = PConvPyListToFloatArray(tmp,&I->Grid);
      }
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,6),&I->Corner[0][0],24);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,7),I->ExtentMin,3);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,8),I->ExtentMax,3);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,9),&I->MapSource);
      if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,10),I->Div,3);
      if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,11),I->Min,3);
      if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,12),I->Max,3);
      if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,13),I->FDim,4);
      if(ok) ok = ((I->Field=IsosurfNewFromPyList(PyList_GetItem(list,14)))!=NULL);
      if(ok) ObjectMapStateRegeneratePoints(I);
    }
  }
  return(ok);
}

static int ObjectMapAllStatesFromPyList(ObjectMap *I,PyObject *list)
{
  int ok=true;
  int a;
  VLACheck(I->State,ObjectMapState,I->NState);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    for(a=0;a<I->NState;a++) {
      ok = ObjectMapStateFromPyList(I->State+a,PyList_GetItem(list,a));
      if(!ok) break;
    }
  }
  return(ok);
}


PyObject *ObjectMapAsPyList(ObjectMap *I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NState));
  PyList_SetItem(result,2,ObjectMapAllStatesAsPyList(I));

  return(PConvAutoNone(result));  
}

int ObjectMapNewFromPyList(PyObject *list,ObjectMap **result)
{
  int ok = true;
  int ll;
  ObjectMap *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  I=ObjectMapNew();
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NState);
  if(ok) ok = ObjectMapAllStatesFromPyList(I,PyList_GetItem(list,2));
  if(ok) {
    (*result) = I;
    ObjectMapUpdateExtents(I);
  } else {
    /* cleanup? */
  }

  return(ok);
}


ObjectMapState *ObjectMapGetState(ObjectMap *I,int state)
{
  ObjectMapState *result = NULL;

  if(state<I->NState)
    result = &I->State[state];
  return(result);
}

ObjectMapState *ObjectMapStatePrime(ObjectMap *I,int state)
{
  ObjectMapState *ms = NULL;
  if(state<0)
    state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectMapState,state);
    I->NState=state+1;
  }
  ms=&I->State[state];
  ObjectMapStateInit(ms);
  return(ms);
}

ObjectMapState *ObjectMapStateGetActive(ObjectMap *I,int state)
{
  ObjectMapState *ms = NULL;
  if(state>=0) {
    if(state<I->NState) {
      ms=&I->State[state];
      if(!ms->Active)
        ms = NULL;
    }
  }
  return(ms);
}


void ObjectMapUpdateExtents(ObjectMap *I)
{
  int a;
  I->Obj.ExtentFlag=false;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active)
      {
        if(!I->Obj.ExtentFlag) {
          copy3f(I->State[a].ExtentMin,I->Obj.ExtentMin);
          copy3f(I->State[a].ExtentMax,I->Obj.ExtentMax);
          I->Obj.ExtentFlag=true;
        } else {
          min3f(I->State[a].ExtentMin,I->Obj.ExtentMin,I->Obj.ExtentMin);
          max3f(I->State[a].ExtentMax,I->Obj.ExtentMax,I->Obj.ExtentMax);
        }
      }
  }
  PRINTFD(FB_ObjectMap)
    " ObjectMapUpdateExtents-DEBUG: ExtentFlag %d\n",I->Obj.ExtentFlag
    ENDFD;
}

int ObjectMapStateSetBorder(ObjectMapState *I,float level)
{
  int result = false;
  int a,b,c;

  c=I->FDim[2]-1;
  for(a=0;a<I->FDim[0];a++) 
    for(b=0;b<I->FDim[1];b++)
      {
        F3(I->Field->data,a,b,0) = level;
        F3(I->Field->data,a,b,c) = level;
      }

  a=I->FDim[0]-1;
  for(b=0;b<I->FDim[1];b++) 
    for(c=0;c<I->FDim[2];c++)
      {
        F3(I->Field->data,0,b,c) = level;
        F3(I->Field->data,a,b,c) = level;
      }

  b=I->FDim[1]-1;
  for(a=0;a<I->FDim[0];a++) 
    for(c=0;c<I->FDim[2];c++)
      {
        F3(I->Field->data,a,0,c) = level;
        F3(I->Field->data,a,b,c) = level;
      }
  return(result);
}

void ObjectMapStatePurge(ObjectMapState *I)
{
  if(I->Field) {
    IsosurfFieldFree(I->Field);
    I->Field=NULL;
  }
  FreeP(I->Origin);
  FreeP(I->Dim);
  FreeP(I->Range);
  FreeP(I->Grid);
  OOFreeP(I->Crystal);
  I->Active=false;
}

static void ObjectMapFree(ObjectMap *I) {

  int a;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active)
      ObjectMapStatePurge(I->State+a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
}

static void ObjectMapUpdate(ObjectMap *I) {
  SceneDirty();
}

static void ObjectMapRender(ObjectMap *I,int state,CRay *ray,Pickable **pick,int pass)
{
  ObjectMapState *ms = NULL;
  if(!pass) {
    if(state<I->NState)
      if(I->State[state].Active)
        ms=&I->State[state];
    
    if(ms) {
      ObjectPrepareContext(&I->Obj,ray);

      if(I->Obj.RepVis[cRepExtent]) {
        if(ray) {
          float *vc;
          vc = ColorGet(I->Obj.Color);
          ray->fColor3fv(ray,vc);
          ray->fSausage3fv(ray,ms->Corner[0],ms->Corner[1],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[0],ms->Corner[2],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[2],ms->Corner[3],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[1],ms->Corner[3],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[0],ms->Corner[4],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[1],ms->Corner[5],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[2],ms->Corner[6],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[3],ms->Corner[7],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[4],ms->Corner[5],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[4],ms->Corner[6],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[6],ms->Corner[7],0.20F,vc,vc);
          ray->fSausage3fv(ray,ms->Corner[5],ms->Corner[7],0.20F,vc,vc);
        } else if(pick&&PMGUI) {
        } else if(PMGUI) {
          ObjectUseColor(&I->Obj);
          glDisable(GL_LIGHTING); 
          glBegin(GL_LINES);
          glVertex3fv(ms->Corner[0]);
          glVertex3fv(ms->Corner[1]);
          
          glVertex3fv(ms->Corner[0]);
          glVertex3fv(ms->Corner[2]);
          
          glVertex3fv(ms->Corner[2]);
          glVertex3fv(ms->Corner[3]);
          
          glVertex3fv(ms->Corner[1]);
          glVertex3fv(ms->Corner[3]);
          
          glVertex3fv(ms->Corner[0]);
          glVertex3fv(ms->Corner[4]);
          
          glVertex3fv(ms->Corner[1]);
          glVertex3fv(ms->Corner[5]);
          
          glVertex3fv(ms->Corner[2]);
          glVertex3fv(ms->Corner[6]);
          
          glVertex3fv(ms->Corner[3]);
          glVertex3fv(ms->Corner[7]);
          
          glVertex3fv(ms->Corner[4]);
          glVertex3fv(ms->Corner[5]);
          
          glVertex3fv(ms->Corner[4]);
          glVertex3fv(ms->Corner[6]);
          
          glVertex3fv(ms->Corner[6]);
          glVertex3fv(ms->Corner[7]);
          
          glVertex3fv(ms->Corner[5]);
          glVertex3fv(ms->Corner[7]);
          
          glEnd();
          glEnable(GL_LIGHTING);
        }
      }
    }
  }
}

void ObjectMapStateInit(ObjectMapState *I) 
{
  ObjectMapStatePurge(I);
  I->Crystal = CrystalNew();
  I->Field = NULL;
  I->Origin = NULL;
  I->Dim = NULL;
  I->Range = NULL;
  I->Grid = NULL;
}
int ObjectMapGetNStates(ObjectMap *I)     
{
  return(I->NState);
}
/*========================================================================*/
ObjectMap *ObjectMapNew(void)
{
  OOAlloc(ObjectMap);

  ObjectInit((CObject*)I);
  I->Obj.type = cObjectMap;

  
  I->NState = 0;
  I->State=VLAMalloc(1,sizeof(ObjectMapState),5,true); /* autozero important */

  I->Obj.RepVis[cRepExtent]=true; 
  I->Obj.fFree = (void (*)(struct CObject *))ObjectMapFree;
  I->Obj.fUpdate =  (void (*)(struct CObject *)) ObjectMapUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, int, CRay *, Pickable **,int))ObjectMapRender;

  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectMapGetNStates;

  return(I);
}
/*========================================================================*/
ObjectMapState *ObjectMapNewStateFromDesc(ObjectMap *I,ObjectMapDesc *md,int state)
{
  int ok=true;
  float v[3];
  int a,b,c,d;
  float *fp;
  ObjectMapState *ms = NULL;

  ms=ObjectMapStatePrime(I,state);
  
  if(I) {
    ms->Origin=Alloc(float,3);
    ms->Range=Alloc(float,3);
    ms->Dim=Alloc(int,3);
    ms->Grid=Alloc(float,3);
    ms->MapSource=cMapSourceDesc;
  }
  switch(md->mode) {
  case cObjectMap_OrthoMinMaxGrid: /* Orthorhombic: min, max, spacing, centered over range  */

    subtract3f(md->MaxCorner,md->MinCorner,v);
    for(a=0;a<3;a++) { if(v[a]<0.0) swap1f(md->MaxCorner+a,md->MinCorner+a); };
    subtract3f(md->MaxCorner,md->MinCorner,v);
    for(a=0;a<3;a++) {
      md->Dim[a] = (int)(v[a]/md->Grid[a]);
      if(md->Dim[a]<1) md->Dim[a]=1;
      if((md->Dim[a]*md->Grid[a])<v[a]) md->Dim[a]++;
    }

    PRINTFB(FB_ObjectMap,FB_Blather)
      " ObjectMap: Dim %d %d %d\n",md->Dim[0],md->Dim[1],md->Dim[2]
      ENDFB;

    average3f(md->MaxCorner,md->MinCorner,v);
    for(a=0;a<3;a++) { md->MinCorner[a] = v[a]-0.5F*md->Dim[a]*md->Grid[a]; }

    if(Feedback(FB_ObjectMap,FB_Blather)) {
      dump3f(md->MinCorner," ObjectMap: MinCorner:");
      dump3f(md->MaxCorner," ObjectMap: MaxCorner:");
      dump3f(md->Grid," ObjectMap: Grid:");
    }
    
    /* now populate the map data structure */

    copy3f(md->MinCorner,ms->Origin);
    copy3f(md->Grid,ms->Grid);
    for(a=0;a<3;a++) ms->Range[a] = md->Grid[a] * md->Dim[a];

    /* these maps start at zero */
    for(a=0;a<3;a++) ms->Min[a]=0; 
    copy3f(md->Dim,ms->Max);

    /* define corners */

    for(a=0;a<8;a++) copy3f(ms->Origin,ms->Corner[a]);

    d = 0;
    for(c=0;c<2;c++) {
      {
        v[2] = (c ? ms->Range[2] : 0.0F);
        for(b=0;b<2;b++) {
          v[1]= (b ? ms->Range[1] : 0.0F);
          for(a=0;a<2;a++) {
            v[0]= (a ? ms->Range[0] : 0.0F);
            add3f(v,ms->Corner[d],ms->Corner[d]);
            d++;
          }
        }
      }
    }
    for(a=0;a<3;a++) ms->FDim[a] = ms->Max[a];
    ms->FDim[3] = 3; 

    ms->Field=IsosurfFieldAlloc(ms->FDim);
    if(!ms->Field) 
      ok=false;
    else {
      for(a=0;a<md->Dim[0];a++) {
        v[0] = md->MinCorner[0] + a * md->Grid[0];
        for(b=0;b<md->Dim[1];b++) {
          v[1] = md->MinCorner[1] + b * md->Grid[1];
          for(c=0;c<md->Dim[2];c++) {
            v[2] = md->MinCorner[2] + c * md->Grid[2];
            fp = F4Ptr(ms->Field->points,a,b,c,0);
            copy3f(v,fp);
          }
        }
      }
    }
    break;
  default:
    ok = false;
  }
  if(ok) {
    switch(md->init_mode) {
    case 0:
      for(a=0;a<md->Dim[0];a++) {      
        for(b=0;b<md->Dim[1];b++) {
          for(c=0;c<md->Dim[2];c++) {
            F3(ms->Field->data,a,b,c)=0.0F;
          }
        }
      }
      break;
    case 1:
      for(a=0;a<md->Dim[0];a++) {      
        for(b=0;b<md->Dim[1];b++) {
          for(c=0;c<md->Dim[2];c++) {
            F3(ms->Field->data,a,b,c)=1.0F;
          }
        }
      }
      break;
    case -2: /* for testing... */
      for(a=0;a<md->Dim[0];a++) {      
        for(b=0;b<md->Dim[1];b++) {
          for(c=0;c<md->Dim[2];c++) {
            F3(ms->Field->data,a,b,c)=(float)sqrt1d(a*a+b*b+c*c);
          }
        }
      }
      break;
    }
  }
  
  if(ok) {
    copy3f(ms->Origin,ms->ExtentMin);
    copy3f(ms->Origin,ms->ExtentMax);
    add3f(ms->Range,ms->ExtentMax,ms->ExtentMax);
    ObjectMapUpdateExtents(I);
  }
  if(!ok) {
    ErrMessage("ObjectMap","Unable to create map");
    ObjectMapFree(I);
    I=NULL;
  } else {
    PRINTFB(FB_ObjectMap,FB_Actions) 
      " ObjectMap: Map created.\n"
      ENDFB;
  }
  
  return(ms);
}
/*========================================================================*/
int ObjectMapCCP4StrToMap(ObjectMap *I,char *CCP4Str,int bytes,int state) {
  
  char *p;
  int *i;
  unsigned int *u;
  unsigned char *uc,c0,c1,c2,c3;
  float *f;
  float dens;
  int a,b,c,d,e;
  float v[3],vr[3],maxd,mind;
  int ok = true;
  int little_endian = 1,map_endian;
  /* CCP4 named from their docs */
  int nc,nr,ns;
  int map_mode;
  int ncstart,nrstart,nsstart;
  int nx,ny,nz;
  float xlen,ylen,zlen,alpha,beta,gamma;
  int n_skip;
  int mapc,mapr,maps;
  int cc[3],xref[3];
  int n_pts;
  double sum,sumsq;
  float mean,stdev;
  int normalize;
  ObjectMapState *ms;

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectMapState,state);
    I->NState=state+1;
  }
  ms=&I->State[state];
  ObjectMapStateInit(ms);

  normalize=(int)SettingGet(cSetting_normalize_ccp4_maps);
  maxd = FLT_MIN;
  mind = FLT_MAX;
  p=CCP4Str;
  little_endian = *((char*)&little_endian);
  map_endian = (*p||*(p+1));

  if(bytes<256*sizeof(int)) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      " ObjectMapCCP4: Map appears to be truncated -- aborting."
      ENDFB;
    return(0);
  }
  if(little_endian!=map_endian) {
    PRINTFB(FB_ObjectMap,FB_Blather)
      " ObjectMapCCP4: Map appears to be reverse endian, swapping...\n"
      ENDFB;
    c = bytes;
    u = (unsigned int*)p;
    uc = (unsigned char *)u;
    while(c>3) {
      c0 = *(uc++);
      c1 = *(uc++);
      c2 = *(uc++);
      c3 = *(uc++);
      uc = (unsigned char *)u;
      *(uc++) = c3;
      *(uc++) = c2;
      *(uc++) = c1;
      *(uc++) = c0;
      u++;
      c-=4;
    }
  }
  i = (int*)p;
  nc=*(i++); /* columns */
  nr=*(i++); /* rows */
  ns=*(i++); /* sections */
  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: NC %d   NR %d   NS %d\n",
    nc,nr,ns
    ENDFB;
  map_mode = *(i++); /* mode */

  
  if(map_mode!=2) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      "ObjectMapCCP4-ERR: Only map mode 2 currently supported (this map is mode %d)",map_mode
      ENDFB;
    return(0);
  }

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: Map is mode %d.\n",map_mode
    ENDFB;

  ncstart = *(i++);
  nrstart = *(i++);
  nsstart = *(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: NCSTART %d   NRSTART %d   NSSTART  %d\n",
    ncstart,nrstart,nsstart
    ENDFB;

  nx = *(i++);
  ny = *(i++);
  nz = *(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: NX %d   NY %d   NZ  %d \n",
    nx,ny,nz
    ENDFB;

  xlen = *(float*)(i++);
  ylen = *(float*)(i++);
  zlen = *(float*)(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: X %8.3f   Y %8.3f  Z  %8.3f \n",
    xlen,ylen,zlen
    ENDFB;

  alpha = *(float*)(i++);
  beta = *(float*)(i++);
  gamma = *(float*)(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: alpha %8.3f   beta %8.3f  gamma %8.3f \n",
    alpha,beta,gamma
    ENDFB;

  mapc = *(i++);
  mapr = *(i++);
  maps = *(i++);

  PRINTFB(FB_ObjectMap,FB_Blather)
    " ObjectMapCCP4: MAPC %d   MAPR %d  MAPS  %d \n",
    mapc,mapr,maps
    ENDFB;

  i+=4;
  n_skip = *(i++);
  
  if(*(i++)) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      "ObjectMapCCP4-ERR: PyMOL doesn't know how to handle skewed maps. Sorry!\n"
      ENDFB;
    return(0);
  }

  n_pts = nc*ns*nr;
  if((unsigned)bytes<(n_skip + sizeof(int)*(256+n_pts))) {
    PRINTFB(FB_ObjectMap,FB_Errors)
      " ObjectMapCCP4: Map appears to be truncated -- aborting.\n"
      ENDFB;
    return(0);
  }

  if(n_pts>1) {
    f = (float*)(p+(sizeof(int)*256)+n_skip);
    c = n_pts;
    sum = 0.0;
    sumsq = 0.0;
    while(c--) {
      sumsq+=(*f)*(*f);
      sum+=*f++;
    }
    mean = (float)(sum/n_pts);
    stdev = (float)sqrt1d((sumsq - (sum*sum/n_pts))/(n_pts-1));

    if(normalize) {
      PRINTFB(FB_ObjectMap,FB_Details)
        " ObjectMapCCP4: Normalizing with mean = %8.6f and stdev = %8.6f.\n",
        mean,stdev
        ENDFB;
    } else {
      PRINTFB(FB_ObjectMap,FB_Details)
        " ObjectMapCCP4: Map will not be normalized.\n ObjectMapCCP4: Current mean = %8.6f and stdev = %8.6f.\n",
        mean,stdev
        ENDFB;
    }

    if(stdev<0.000001)
      stdev = 1.0;
    
  } else {
    mean = 1.0;
    stdev = 1.0;
  }
  
  f = (float*)(p+(sizeof(int)*256)+n_skip);
  mapc--; /* convert to C indexing... */
  mapr--;
  maps--;

  xref[maps]=0;
  xref[mapr]=1;
  xref[maps]=2;

  xref[maps]=0;
  xref[mapr]=1;
  xref[mapc]=2;

  ms->Div[0] = nx;
  ms->Div[1] = ny;
  ms->Div[2] = nz;

  ms->FDim[mapc] = nc;
  ms->Min[mapc] = ncstart;
  ms->Max[mapc] = nc+ncstart;

  ms->FDim[mapr] = nr;
  ms->Min[mapr] = nrstart;
  ms->Max[mapr] = nr+nrstart;

  ms->FDim[maps] = ns;
  ms->Min[maps] = nsstart;
  ms->Max[maps] = ns+nsstart;
  
  ms->Crystal->Dim[0] = xlen;
  ms->Crystal->Dim[1] = ylen;
  ms->Crystal->Dim[2] = zlen;

  ms->Crystal->Angle[0] = alpha;
  ms->Crystal->Angle[1] = beta;
  ms->Crystal->Angle[2] = gamma;


  ms->FDim[3]=3;
  if(!(ms->FDim[0]&&ms->FDim[1]&&ms->FDim[2])) 
    ok=false;
  else {
    CrystalUpdate(ms->Crystal);
    CrystalDump(ms->Crystal);
    fflush(stdout);
    ms->Field=IsosurfFieldAlloc(ms->FDim);
    ms->MapSource = cMapSourceCCP4;
    ms->Field->save_points=false;

    for(cc[maps]=0;cc[maps]<ms->FDim[maps];cc[maps]++)
      {
        v[maps]=(cc[maps]+ms->Min[maps])/((float)ms->Div[maps]);

        for(cc[mapr]=0;cc[mapr]<ms->FDim[mapr];cc[mapr]++) {
          v[mapr]=(cc[mapr]+ms->Min[mapr])/((float)ms->Div[mapr]);

          for(cc[mapc]=0;cc[mapc]<ms->FDim[mapc];cc[mapc]++) {
            v[mapc]=(cc[mapc]+ms->Min[mapc])/((float)ms->Div[mapc]);

            if(normalize) 
              dens = (*f-mean)/stdev;
            else 
              dens = *f;
            F3(ms->Field->data,cc[0],cc[1],cc[2]) = dens;
            if(maxd<*f) maxd = dens;
            if(mind>*f) mind = dens;
            f++;
            transform33f3f(ms->Crystal->FracToReal,v,vr);
            for(e=0;e<3;e++) 
              F4(ms->Field->points,cc[0],cc[1],cc[2],e) = vr[e];
          }
        }
      }
  }
  if(ok) {
    d = 0;
    for(c=0;c<ms->FDim[2];c+=(ms->FDim[2]-1))
      {
        v[2]=(c+ms->Min[2])/((float)ms->Div[2]);
        for(b=0;b<ms->FDim[1];b+=(ms->FDim[1]-1)) {
          v[1]=(b+ms->Min[1])/((float)ms->Div[1]);
          for(a=0;a<ms->FDim[0];a+=(ms->FDim[0]-1)) {
            v[0]=(a+ms->Min[0])/((float)ms->Div[0]);
            transform33f3f(ms->Crystal->FracToReal,v,vr);
            copy3f(vr,ms->Corner[d]);
            d++;
          }
        }
      }
  }
  
  if(ok) {
    v[2]=(ms->Min[2])/((float)ms->Div[2]);
    v[1]=(ms->Min[1])/((float)ms->Div[1]);
    v[0]=(ms->Min[0])/((float)ms->Div[0]);
    
    transform33f3f(ms->Crystal->FracToReal,v,ms->ExtentMin);
    
    v[2]=((ms->FDim[2]-1)+ms->Min[2])/((float)ms->Div[2]);
    v[1]=((ms->FDim[1]-1)+ms->Min[1])/((float)ms->Div[1]);
    v[0]=((ms->FDim[0]-1)+ms->Min[0])/((float)ms->Div[0]);
    
    transform33f3f(ms->Crystal->FracToReal,v,ms->ExtentMax);
  }
#ifdef _UNDEFINED
  printf("%d %d %d %d %d %d %d %d %d\n",
         ms->Div[0],
         ms->Min[0],
         ms->Max[0],
         ms->Div[1],
         ms->Min[1],
         ms->Max[1],
         ms->Div[2],
         ms->Min[2],
         ms->Max[2]);
  printf("Okay? %d\n",ok);
  fflush(stdout);
#endif
  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    ms->Active=true;
    ObjectMapUpdateExtents(I);
    printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
  }
  return(ok);
}
/*========================================================================*/
static int ObjectMapPHIStrToMap(ObjectMap *I,char *PHIStr,int bytes,int state) {
  
  char *p;
  float dens,dens_rev;
  int a,b,c,d,e;
  float v[3],maxd,mind;
  int ok = true;
  int little_endian = 1;
  /* PHI named from their docs */
  int map_endian = 0;
  int map_dim;

  ObjectMapState *ms;

  char cc[MAXLINELEN];
  char *rev;

  little_endian = *((char*)&little_endian);


  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectMapState,state);
    I->NState=state+1;
  }
  ms=&I->State[state];
  ObjectMapStateInit(ms);

  maxd = FLT_MIN;
  mind = FLT_MAX;
  p=PHIStr;

  p+=4;
  ParseNCopy(cc,p,20);
  PRINTFB(FB_ObjectMap,FB_Details)
    " PHIMapToStr: %s\n",cc
    ENDFB;
  p+=20;
  p+=4;

  p+=4;
  ParseNCopy(cc,p,10);
  PRINTFB(FB_ObjectMap,FB_Details)
    " PHIMapToStr: %s\n",cc
    ENDFB;
  p+=10;
  ParseNCopy(cc,p,60);
  PRINTFB(FB_ObjectMap,FB_Details)
    " PHIMapToStr: %s\n",cc
    ENDFB;
  p+=60;
  p+=4;

  p+=4;

  map_dim = 65;

  ms->FDim[0] = map_dim;
  ms->FDim[1] = map_dim;
  ms->FDim[2] = map_dim;
  ms->FDim[3] = 3;

  ms->Div[0] = (map_dim-1)/2;
  ms->Div[1] = (map_dim-1)/2;
  ms->Div[2] = (map_dim-1)/2;
  ms->Min[0] = -ms->Div[0];
  ms->Min[1] = -ms->Div[1];
  ms->Min[2] = -ms->Div[2];
  ms->Max[0] = ms->Div[0];
  ms->Max[1] = ms->Div[1];
  ms->Max[2] = ms->Div[2];

  ms->Field=IsosurfFieldAlloc(ms->FDim);
  ms->MapSource = cMapSourcePHI;
  ms->Field->save_points=false;

  rev = (char*)&dens_rev;

  for(c=0;c<ms->FDim[2];c++) { /* z y x ordering into c b a  so that x = a, etc. */
    for(b=0;b<ms->FDim[1];b++) {
      for(a=0;a<ms->FDim[0];a++) {

        if(little_endian!=map_endian) {
          rev[0]=p[3];
          rev[1]=p[2];
          rev[2]=p[1];
          rev[3]=p[0];
        } else {
          rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
          rev[1]=p[1];
          rev[2]=p[2];
          rev[3]=p[3];
        }
        dens = *((float*)rev);
        F3(ms->Field->data,a,b,c) = dens;
        if(maxd<dens) maxd = dens;
        if(mind>dens) mind = dens;
        p+=4;
      }
    }
  }
  p+=4;


  p+=4;
  ParseNCopy(cc,p,16);
  PRINTFB(FB_ObjectMap,FB_Details)
    " PHIMapToStr: %s\n",cc
    ENDFB;
  p+=16;
  p+=4;

  ms->Grid = Alloc(float,3);
  p+=4;
  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];
  }
  ms->Grid[0] = 1.0F/(*((float*)rev));
  ms->Grid[1] = ms->Grid[0];
  ms->Grid[2] = ms->Grid[0];
  p+=4;

  ms->Origin = Alloc(float,3);
  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];;
  }
  ms->Origin[0] = *((float*)rev);
  p+=4;

  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];;
  }
  ms->Origin[1] = *((float*)rev);
  p+=4;
  if(little_endian!=map_endian) {
    rev[0]=p[3];
    rev[1]=p[2];
    rev[2]=p[1];
    rev[3]=p[0];
  } else {
    rev[0]=p[0]; /* gotta go char by char because of memory alignment issues ... */
    rev[1]=p[1];
    rev[2]=p[2];
    rev[3]=p[3];;
  }
  ms->Origin[2] = *((float*)rev);
  p+=4;

  p+=4;

  if(ok) {
    for(e=0;e<3;e++) {
      ms->ExtentMin[e] = ms->Origin[e]+ms->Grid[e]*ms->Min[e];
      ms->ExtentMax[e] = ms->Origin[e]+ms->Grid[e]*ms->Max[e];
    }
  }

  for(c=0;c<ms->FDim[2];c++) {
    v[2]=ms->Origin[2]+ms->Grid[2]*(c+ms->Min[2]);
    for(b=0;b<ms->FDim[1];b++) {
      v[1]=ms->Origin[1]+ms->Grid[1]*(b+ms->Min[1]);
      for(a=0;a<ms->FDim[0];a++) {
        v[0]=ms->Origin[0]+ms->Grid[0]*(a+ms->Min[0]);
        for(e=0;e<3;e++) {
          F4(ms->Field->points,a,b,c,e) = v[e];
        }
      }
    }
  }

  d=0;
  for(c=0;c<ms->FDim[2];c+=ms->FDim[2]-1) {
    v[2]=ms->Origin[2]+ms->Grid[2]*(c+ms->Min[2]);

    for(b=0;b<ms->FDim[1];b+=ms->FDim[1]-1) {
      v[1]=ms->Origin[1]+ms->Grid[1]*(b+ms->Min[1]);

      for(a=0;a<ms->FDim[2];a+=ms->FDim[2]-1) {
        v[0]=ms->Origin[0]+ms->Grid[0]*(a+ms->Min[0]);
        copy3f(v,ms->Corner[d]);
        d++;
      }
    }
  }

  /* interpolation test code 
  { 
    float test[3];
    float result;
    float cmp;

    for(c=0;c<ms->FDim[0]-1;c++) {
      for(b=0;b<ms->FDim[1]-1;b++) {
        for(a=0;a<ms->FDim[2]-1;a++) {
          for(e=0;e<3;e++) {
            test[e] = (F4(ms->Field->points,a,b,c,e)+
                       F4(ms->Field->points,a,b,c,e))/2.0;
          }
          ObjectMapStateInterpolate(ms,test,&result,1);
          cmp = (F3(ms->Field->data,a,b,c)+
                 F3(ms->Field->data,a,b,c))/2.0;
          if(fabs(cmp-result)>0.001) {
            printf("%d %d %d\n",a,b,c);
            printf("%8.5f %8.5f\n",
                   cmp,
                   result);
          }
        }
      }
    }
  }
  */

  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    ms->Active=true;
    ObjectMapUpdateExtents(I);
    printf(" ObjectMap: Map Read.  Range = %5.6f to %5.6f\n",mind,maxd);
  }
  return(ok);
}
/*========================================================================*/
int ObjectMapXPLORStrToMap(ObjectMap *I,char *XPLORStr,int state) {
  
  char *p;
  int a,b,c,d,e;
  float v[3],vr[3],dens,maxd,mind;
  char cc[MAXLINELEN];
  int n;
  int ok = true;
  ObjectMapState *ms;

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectMapState,state);
    I->NState=state+1;
  }

  ms=&I->State[state];
  ObjectMapStateInit(ms);

  maxd = FLT_MIN;
  mind = FLT_MAX;
  p=XPLORStr;

  while(*p) {
    p = ParseNCopy(cc,p,8);
    if(!*cc) 
      p = ParseNextLine(p);
    else if(sscanf(cc,"%i",&n)==1) {
      p=ParseWordCopy(cc,p,MAXLINELEN);
      if(strstr(cc,"!NTITLE")||(!*cc)) {
        p=ParseNextLine(p);
        while(n--) {
          p=ParseNextLine(p);          
        } 
      } else if(strstr(cc,"REMARKS")) {
        p=ParseNextLine(p);          
      } else {
        break;
      }
    }
  }
  if(*p) { /* n contains first dimension */
    ms->Div[0]=n;
    if(sscanf(cc,"%i",&ms->Min[0])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Max[0])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Div[1])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Min[1])!=1) ok=false;    
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Max[1])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Div[2])!=1) ok=false;
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Min[2])!=1) ok=false;    
    p = ParseNCopy(cc,p,8); if(sscanf(cc,"%i",&ms->Max[2])!=1) ok=false;
    p=ParseNextLine(p);
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&ms->Crystal->Dim[0])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&ms->Crystal->Dim[1])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&ms->Crystal->Dim[2])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&ms->Crystal->Angle[0])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&ms->Crystal->Angle[1])!=1) ok=false;
    p = ParseNCopy(cc,p,12); if(sscanf(cc,"%f",&ms->Crystal->Angle[2])!=1) ok=false;
    p=ParseNextLine(p);
    p = ParseNCopy(cc,p,3);
    if(strcmp(cc,"ZYX")) ok=false;
    p=ParseNextLine(p);
    
  } else {
    ok=false;
  }
  if(ok) {
    
    ms->FDim[0]=ms->Max[0]-ms->Min[0]+1;
    ms->FDim[1]=ms->Max[1]-ms->Min[1]+1;
    ms->FDim[2]=ms->Max[2]-ms->Min[2]+1;
    ms->FDim[3]=3;
    if(!(ms->FDim[0]&&ms->FDim[1]&&ms->FDim[2])) 
      ok=false;
    else {
      CrystalUpdate(ms->Crystal);
      ms->Field=IsosurfFieldAlloc(ms->FDim);
      ms->MapSource = cMapSourceXPLOR;
      ms->Field->save_points=false;
      for(c=0;c<ms->FDim[2];c++)
        {
          v[2]=(c+ms->Min[2])/((float)ms->Div[2]);
          p=ParseNextLine(p);
          for(b=0;b<ms->FDim[1];b++) {
            v[1]=(b+ms->Min[1])/((float)ms->Div[1]);
            for(a=0;a<ms->FDim[0];a++) {
              v[0]=(a+ms->Min[0])/((float)ms->Div[0]);
              p=ParseNCopy(cc,p,12);
              if(!cc[0]) {
                p=ParseNextLine(p);
                p=ParseNCopy(cc,p,12);                
              }
              if(sscanf(cc,"%f",&dens)!=1) {
                ok=false;
              } else {
                F3(ms->Field->data,a,b,c) = dens;
                if(maxd<dens) maxd = dens;
                if(mind>dens) mind = dens;
              }
              transform33f3f(ms->Crystal->FracToReal,v,vr);
              for(e=0;e<3;e++) {
                F4(ms->Field->points,a,b,c,e) = vr[e];
              }
            }
          }
          p=ParseNextLine(p);
        }
      if(ok) {
        d = 0;
        for(c=0;c<ms->FDim[2];c+=(ms->FDim[2]-1))
          {
            v[2]=(c+ms->Min[2])/((float)ms->Div[2]);
            for(b=0;b<ms->FDim[1];b+=(ms->FDim[1]-1)) {
              v[1]=(b+ms->Min[1])/((float)ms->Div[1]);
              for(a=0;a<ms->FDim[0];a+=(ms->FDim[0]-1)) {
                v[0]=(a+ms->Min[0])/((float)ms->Div[0]);
                transform33f3f(ms->Crystal->FracToReal,v,vr);
                copy3f(vr,ms->Corner[d]);
                d++;
              }
            }
          }
      }
    }
  }

  if(ok) {
    v[2]=(ms->Min[2])/((float)ms->Div[2]);
    v[1]=(ms->Min[1])/((float)ms->Div[1]);
    v[0]=(ms->Min[0])/((float)ms->Div[0]);

    transform33f3f(ms->Crystal->FracToReal,v,ms->ExtentMin);
    
    v[2]=((ms->FDim[2]-1)+ms->Min[2])/((float)ms->Div[2]);
    v[1]=((ms->FDim[1]-1)+ms->Min[1])/((float)ms->Div[1]);
    v[0]=((ms->FDim[0]-1)+ms->Min[0])/((float)ms->Div[0]);

    transform33f3f(ms->Crystal->FracToReal,v,ms->ExtentMax);

  }
#ifdef _UNDEFINED
    printf("%d %d %d %d %d %d %d %d %d\n",
           ms->Div[0],
           ms->Min[0],
           ms->Max[0],
           ms->Div[1],
           ms->Min[1],
           ms->Max[1],
           ms->Div[2],
           ms->Min[2],
           ms->Max[2]);
    printf("Okay? %d\n",ok);
    fflush(stdout);
#endif
  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    ms->Active=true;
    ObjectMapUpdateExtents(I);
    printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
  }
    
  return(ok);
}
/*========================================================================*/
ObjectMap *ObjectMapReadXPLORStr(ObjectMap *I,char *XPLORStr,int state)
{
  int ok=true;
  int isNew = true;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;
  if(ok) {
	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    ObjectMapXPLORStrToMap(I,XPLORStr,state);
    SceneChanged();
    SceneCountFrames();
  }
  return(I);
}
/*========================================================================*/
ObjectMap *ObjectMapReadCCP4Str(ObjectMap *I,char *XPLORStr,int bytes,int state)
{
  int ok=true;
  int isNew = true;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;
  if(ok) {
	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    ObjectMapCCP4StrToMap(I,XPLORStr,bytes,state);
    SceneChanged();
    SceneCountFrames();
  }
  return(I);
}
/*========================================================================*/
ObjectMap *ObjectMapLoadCCP4File(ObjectMap *obj,char *fname,int state)
{
  ObjectMap *I = NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;
  float mat[9];

  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMapLoadCCP4File","Unable to open file!");
  else
	 {
		if(Feedback(FB_ObjectMap,FB_Actions))
		  {
			printf(" ObjectMapLoadCCP4File: Loading from '%s'.\n",fname);
		  }
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		fclose(f);

		I=ObjectMapReadCCP4Str(obj,buffer,size,state);

		mfree(buffer);
      if(state<0)
        state=I->NState-1;
      if(state<I->NState) {
        ObjectMapState *ms;
        ms = &I->State[state];
        if(ms->Active) {
          CrystalDump(ms->Crystal);
          multiply33f33f(ms->Crystal->FracToReal,ms->Crystal->RealToFrac,mat);
        }
      }
    }
#ifdef _UNDEFINED
  {int a;
  for(a=0;a<9;a++)
    printf("%10.5f\n",mat[a]);
  }
#endif
  return(I);

}

/*========================================================================*/
static ObjectMap *ObjectMapReadPHIStr(ObjectMap *I,char *MapStr,int bytes,int state)
{
  int ok=true;
  int isNew = true;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;
  if(ok) {
	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    ObjectMapPHIStrToMap(I,MapStr,bytes,state);
    SceneChanged();
    SceneCountFrames();
  }
  return(I);
}
/*========================================================================*/
ObjectMap *ObjectMapLoadPHIFile(ObjectMap *obj,char *fname,int state)
{
  ObjectMap *I = NULL;
  int ok=true;
  FILE *f;
  long size;
  char *buffer,*p;
  float mat[9];

  f=fopen(fname,"rb");
  if(!f)
	 ok=ErrMessage("ObjectMapLoadPHIFile","Unable to open file!");
  else
	 {
		if(Feedback(FB_ObjectMap,FB_Actions))
		  {
			printf(" ObjectMapLoadPHIFile: Loading from '%s'.\n",fname);
		  }
		
		fseek(f,0,SEEK_END);
      size=ftell(f);
		fseek(f,0,SEEK_SET);

		buffer=(char*)mmalloc(size);
		ErrChkPtr(buffer);
		p=buffer;
		fseek(f,0,SEEK_SET);
		fread(p,size,1,f);
		fclose(f);

		I=ObjectMapReadPHIStr(obj,buffer,size,state);

		mfree(buffer);
      if(state<0)
        state=I->NState-1;
      if(state<I->NState) {
        ObjectMapState *ms;
        ms = &I->State[state];
        if(ms->Active) {
          CrystalDump(ms->Crystal);
          multiply33f33f(ms->Crystal->FracToReal,ms->Crystal->RealToFrac,mat);
        }
      }
    }
#ifdef _UNDEFINED
  {int a;
  for(a=0;a<9;a++)
    printf("%10.5f\n",mat[a]);
  }
#endif
  return(I);

}

/*========================================================================*/
ObjectMap *ObjectMapLoadXPLORFile(ObjectMap *obj,char *fname,int state,int is_file)
{
  ObjectMap *I = NULL;
  int ok=true;
  FILE *f = NULL;
  long size;
  char *buffer,*p;
  float mat[9];

  if(is_file) {
    f=fopen(fname,"rb");
    if(!f)
      ok=ErrMessage("ObjectMapLoadXPLORFile","Unable to open file!");
  }
  if(ok)
	 {
      if(Feedback(FB_ObjectMap,FB_Actions))
        {
          if(is_file) {
            printf(" ObjectMapLoadXPLORFile: Loading from '%s'.\n",fname);
          } else {
            printf(" ObjectMapLoadXPLORFile: Loading...\n");
          }
          
        }
      
      if(is_file) {
        fseek(f,0,SEEK_END);
        size=ftell(f);
        fseek(f,0,SEEK_SET);
        
        buffer=(char*)mmalloc(size+255);
        ErrChkPtr(buffer);
        p=buffer;
        fseek(f,0,SEEK_SET);
        fread(p,size,1,f);
        p[size]=0;
        fclose(f);
      } else {
        buffer = fname;
      }
      
		I=ObjectMapReadXPLORStr(obj,buffer,state);

      if(is_file) 
        mfree(buffer);
      if(state<0)
        state=I->NState-1;
      if(state<I->NState) {
        ObjectMapState *ms;
        ms = &I->State[state];
        if(ms->Active) {
          CrystalDump(ms->Crystal);
          multiply33f33f(ms->Crystal->FracToReal,ms->Crystal->RealToFrac,mat);
        }
      }
    }
#ifdef _UNDEFINED
  {int a;
  for(a=0;a<9;a++)
    printf("%10.5f\n",mat[a]);
  }
#endif
  return(I);

}
/*========================================================================*/
int ObjectMapSetBorder(ObjectMap *I,float level)
{
  int a;
  int result=false;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active)
      result = result && ObjectMapStateSetBorder(&I->State[a],level);
  }
  return(result);
}
/*========================================================================*/
int ObjectMapNumPyArrayToMapState(ObjectMapState *ms,PyObject *ary) {
  
  int a,b,c,d,e;
  float v[3],dens,maxd,mind;
  int ok = true;
#ifdef _PYMOL_NUMPY
  MyArrayObject *pao;
#endif


#ifdef _PYMOL_NUMPY
  pao = (MyArrayObject*)ary;
#endif
  maxd = FLT_MIN;
  mind = FLT_MAX;
  if(ok) {
    ms->FDim[0]=ms->Dim[0];
    ms->FDim[1]=ms->Dim[1];
    ms->FDim[2]=ms->Dim[2];
    ms->FDim[3]=3;

    if(!(ms->FDim[0]&&ms->FDim[1]&&ms->FDim[2])) 
      ok=false;
    else {
      ms->Field=IsosurfFieldAlloc(ms->FDim);
      for(c=0;c<ms->FDim[2];c++)
        {
          v[2]=ms->Origin[2]+ms->Grid[2]*c;
          for(b=0;b<ms->FDim[1];b++) {
            v[1]=ms->Origin[1]+ms->Grid[1]*b;
            for(a=0;a<ms->FDim[0];a++) {
              v[0]=ms->Origin[0]+ms->Grid[0]*a;
#ifdef _PYMOL_NUMPY
              dens = (float)(*((double*)
                (pao->data+
                 (pao->strides[0]*a)+
                 (pao->strides[1]*b)+
                 (pao->strides[2]*c))));
#else
              dens = 0.0;
#endif
              F3(ms->Field->data,a,b,c) = dens;
              if(maxd<dens) maxd = dens;
              if(mind>dens) mind = dens;
              for(e=0;e<3;e++) 
                F4(ms->Field->points,a,b,c,e) = v[e];
            }
          }
        }
      d = 0;
      for(c=0;c<ms->FDim[2];c+=(ms->FDim[2]-1))
        {
          v[2]=ms->Origin[2]+ms->Grid[2]*c;
          for(b=0;b<ms->FDim[1];b+=(ms->FDim[1]-1)) {
            v[1]=ms->Origin[1]+ms->Grid[1]*b;
            for(a=0;a<ms->FDim[0];a+=(ms->FDim[0]-1)) {
              v[0]=ms->Origin[0]+ms->Grid[0]*a;
              copy3f(v,ms->Corner[d]);
              d++;
            }
          }
        }
    }
  }
  if(ok) {
    copy3f(ms->Origin,ms->ExtentMin);
    copy3f(ms->Origin,ms->ExtentMax);
    add3f(ms->Range,ms->ExtentMax,ms->ExtentMax);
  }
  if(!ok) {
    ErrMessage("ObjectMap","Error reading map");
  } else {
    ms->Active=true;
    if(Feedback(FB_ObjectMap,FB_Actions)) {
      printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
    }
  }
  return(ok);
}
/*========================================================================*/
ObjectMap *ObjectMapLoadChemPyBrick(ObjectMap *I,PyObject *Map,
                                           int state,int discrete)
{
  int ok=true;
  int isNew = true;
  PyObject *tmp;
  ObjectMapState *ms;


  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }

    if(state<0) state=I->NState;
    if(I->NState<=state) {
      VLACheck(I->State,ObjectMapState,state);
      I->NState=state+1;
    }
    ms=&I->State[state];
    ObjectMapStateInit(ms);
    
    if(PyObject_HasAttrString(Map,"origin")&&
       PyObject_HasAttrString(Map,"dim")&&
       PyObject_HasAttrString(Map,"range")&&
       PyObject_HasAttrString(Map,"grid")&&
       PyObject_HasAttrString(Map,"lvl"))
      {
        tmp = PyObject_GetAttrString(Map,"origin");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&ms->Origin);
          Py_DECREF(tmp);
        } else 
          ok=ErrMessage("ObjectMap","missing brick origin.");
        tmp = PyObject_GetAttrString(Map,"dim");
        if(tmp) {
          PConvPyListToIntArray(tmp,&ms->Dim);
          Py_DECREF(tmp);
        } else 
          ok=ErrMessage("ObjectMap","missing brick dimension.");
        tmp = PyObject_GetAttrString(Map,"range");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&ms->Range);
          Py_DECREF(tmp);
        } else 
          ok=ErrMessage("ObjectMap","missing brick range.");
        tmp = PyObject_GetAttrString(Map,"grid");
        if(tmp) {
          PConvPyListToFloatArray(tmp,&ms->Grid);
          Py_DECREF(tmp);
        } else
          ok=ErrMessage("ObjectMap","missing brick grid.");
        tmp = PyObject_GetAttrString(Map,"lvl");
        if(tmp) {

          ObjectMapNumPyArrayToMapState(ms,tmp);	 

          Py_DECREF(tmp);
        } else
          ok=ErrMessage("ObjectMap","missing brick density.");

      }
    SceneChanged();
    SceneCountFrames();
    if(ok) {
      ms->Active=true;
      ObjectMapUpdateExtents(I);
          
    }
  }
  return(I);
}

/*========================================================================*/
ObjectMap *ObjectMapLoadChemPyMap(ObjectMap *I,PyObject *Map,
                                  int state,int discrete)
{

  int ok=true;
  int isNew = true;
  float *cobj;
  WordType format;
  float v[3],vr[3],dens,maxd,mind;
  int a,b,c,d,e;
  ObjectMapState *ms;



  /*  
  double test[1000];
  for(a=0;a<1000;a++) {
      test[a]=rand()/(1.0+INT_MAX);
      }
      PyObject_SetAttrString(Map,"c_object",
      PyCObject_FromVoidPtr(test,NULL));
  */
  maxd = FLT_MIN;
  mind = FLT_MAX;

  if(!I) 
	 isNew=true;
  else 
	 isNew=false;

  if(ok) {

	 if(isNew) {
		I=(ObjectMap*)ObjectMapNew();
		isNew = true;
	 } else {
		isNew = false;
	 }
    
    if(state<0) state=I->NState;
    if(I->NState<=state) {
      VLACheck(I->State,ObjectMapState,state);
      I->NState=state+1;
    }
    ms=&I->State[state];
    ObjectMapStateInit(ms);

    if(!PConvAttrToStrMaxLen(Map,"format",format,sizeof(WordType)-1))
      ok=ErrMessage("LoadChemPyMap","bad 'format' parameter.");
    else if(!PConvAttrToFloatArrayInPlace(Map,"cell_dim",ms->Crystal->Dim,3))
      ok=ErrMessage("LoadChemPyMap","bad 'cell_dim' parameter.");
    else if(!PConvAttrToFloatArrayInPlace(Map,"cell_ang",ms->Crystal->Angle,3))
      ok=ErrMessage("LoadChemPyMap","bad 'cell_ang' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map,"cell_div",ms->Div,3))
      ok=ErrMessage("LoadChemPyMap","bad 'cell_div' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map,"first",ms->Min,3))
      ok=ErrMessage("LoadChemPyMap","bad 'first' parameter.");
    else if(!PConvAttrToIntArrayInPlace(Map,"last",ms->Max,3))
      ok=ErrMessage("LoadChemPyMap","bad 'last' parameter.");

    if(ok) {
      if (strcmp(format,"CObjectZYXfloat")==0) {
        ok = PConvAttrToPtr(Map,"c_object",(void**)&cobj);
        if(!ok)
          ErrMessage("LoadChemPyMap","CObject unreadable.");        
      } else {
        ok=ErrMessage("LoadChemPyMap","unsupported format.");        
      }
    }
    /* good to go */

    if(ok) {
      if (strcmp(format,"CObjectZYXfloat")==0) {

        ms->FDim[0]=ms->Max[0]-ms->Min[0]+1;
        ms->FDim[1]=ms->Max[1]-ms->Min[1]+1;
        ms->FDim[2]=ms->Max[2]-ms->Min[2]+1;
        if(Feedback(FB_ObjectMap,FB_Actions)) {
          printf(" LoadChemPyMap: CObjectZYXdouble %dx%dx%d\n",ms->FDim[0],ms->FDim[1],ms->FDim[2]);        
        }
        ms->FDim[3]=3;
        if(!(ms->FDim[0]&&ms->FDim[1]&&ms->FDim[2])) 
          ok=false;
        else {
          CrystalUpdate(ms->Crystal);
          ms->Field=IsosurfFieldAlloc(ms->FDim);
          for(c=0;c<ms->FDim[2];c++)
            {
              v[2]=(c+ms->Min[2])/((float)ms->Div[2]);
              for(b=0;b<ms->FDim[1];b++) {
                v[1]=(b+ms->Min[1])/((float)ms->Div[1]);
                for(a=0;a<ms->FDim[0];a++) {
                  v[0]=(a+ms->Min[0])/((float)ms->Div[0]);
                  
                  dens = *(cobj++);

                  F3(ms->Field->data,a,b,c) = dens;
                  if(maxd<dens) maxd = dens;
                  if(mind>dens) mind = dens;
                  transform33f3f(ms->Crystal->FracToReal,v,vr);
                  for(e=0;e<3;e++) 
                    F4(ms->Field->points,a,b,c,e) = vr[e];
                }
              }
            }

          if(ok) {
            d = 0;
            for(c=0;c<ms->FDim[2];c+=(ms->FDim[2]-1))
              {
                v[2]=(c+ms->Min[2])/((float)ms->Div[2]);
                for(b=0;b<ms->FDim[1];b+=(ms->FDim[1]-1)) {
                  v[1]=(b+ms->Min[1])/((float)ms->Div[1]);
                  for(a=0;a<ms->FDim[0];a+=(ms->FDim[0]-1)) {
                    v[0]=(a+ms->Min[0])/((float)ms->Div[0]);
                    transform33f3f(ms->Crystal->FracToReal,v,vr);
                    copy3f(vr,ms->Corner[d]);
                    d++;
                  }
                }
              }
          }
        }
      }
    }
    
    if(ok) {
      CrystalDump(ms->Crystal);
      
      v[2]=(ms->Min[2])/((float)ms->Div[2]);
      v[1]=(ms->Min[1])/((float)ms->Div[1]);
      v[0]=(ms->Min[0])/((float)ms->Div[0]);
      
      transform33f3f(ms->Crystal->FracToReal,v,ms->ExtentMin);
      
      v[2]=((ms->FDim[2]-1)+ms->Min[2])/((float)ms->Div[2]);
      v[1]=((ms->FDim[1]-1)+ms->Min[1])/((float)ms->Div[1]);
      v[0]=((ms->FDim[0]-1)+ms->Min[0])/((float)ms->Div[0]);
      
      transform33f3f(ms->Crystal->FracToReal,v,ms->ExtentMax);

    }

    if(!ok) {
      ErrMessage("ObjectMap","Error reading map");
    } else {
      ms->Active=true;
      ObjectMapUpdateExtents(I);
		if(Feedback(FB_ObjectMap,FB_Actions)) {
        printf(" ObjectMap: Map Read.  Range = %5.3f to %5.3f\n",mind,maxd);
      }
    }

    if(ok) {
      SceneChanged();
      SceneCountFrames();
    }
  }
  return(I);
}

