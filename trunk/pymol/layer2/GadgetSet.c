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

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Scene.h"
#include"GadgetSet.h"
#include"Color.h"
#include"PConv.h"
#include"main.h"

void GadgetSetUpdate(GadgetSet *I);
void GadgetSetFree(GadgetSet *I);
void GadgetSetRender(GadgetSet *I,CRay *ray,Pickable **pick,int pass);
void GadgetSetStrip(GadgetSet *I);
void GadgetSetInvalidateRep(GadgetSet *I,int type,int level);

int GadgetSetGetVertex(GadgetSet *I,int index,int base,float *v)
{
  int ok=true;
  float *v0,*v1;
  if(index<I->NCoord) {
    v0=I->Coord+3*index;
    if(base<0) {
      copy3f(v0,v);
      if(index) add3f(I->Coord,v,v);
    } else if(base<I->NCoord) {
      v1=I->Coord+3*base;
      add3f(v1,v0,v);
      if(index) add3f(I->Coord,v,v);
    } else {
      ok=false;
    }
  } else
    ok=false;
  return(ok);
}

int GadgetSetSetVertex(GadgetSet *I,int index,int base, float *v)
{
  int ok=true;
  float *v0,*v1;
  if(index<I->NCoord) {
    v0=I->Coord+3*index;
    if(base<0) {
      copy3f(v,v0);
      if(index) subtract3f(v0,I->Coord,v0);
    } else if(base<I->NCoord) {
      v1=I->Coord+3*base;
      subtract3f(v,v1,v0);
      if(index) subtract3f(v0,I->Coord,v0);
    } else {
      ok=false;
    }
  } else
    ok=false;
  return(ok);
}

int GadgetSetFetch(GadgetSet *I,float *inp,float *out)
{
  int ok=true;
  int idx = (int)inp[1],base;
  float *v,*b;
  switch((int)inp[0]) {
  case 0: /* absolute global vertex  */
    if(idx<I->NCoord) {
      v=I->Coord+3*idx;
      copy3f(v,out);
    } else 
      ok=false;
    break;
  case 1: /* relative vertex in gadget (relative to gadget origin)  */
    if(idx<I->NCoord) {
      v=I->Coord+3*idx;
      copy3f(v,out);
      if(idx) add3f(I->Coord,out,out);
    } else 
      ok=false;
    break;
  case 2: /* offset vertex in gadget (relative to vertex relative to gadget origin) */
    base = (int)inp[2];
    if((idx<I->NCoord)&&
       (base<I->NCoord)) {
      v = I->Coord+3*idx;
      b = I->Coord+3*base;
      add3f(b,v,out);
      if(idx) add3f(I->Coord,out,out);
    } else 
      ok=false;
    break;
  case 3: /* normal */
    if(idx<I->NNormal) {
      v=I->Normal+3*idx;
      copy3f(v,out);
    } else 
      ok=false;
    break;
  case 4: /* color */
    if(idx<I->NColor) {
      v=I->Color+3*idx;
      copy3f(v,out);
    } else 
      ok=false;
    break;
  default:
    ok=false;
    break;
  }
  return(ok);
}

int GadgetSetFetchColor(GadgetSet *I,float *inp,float *out)
{
  int ok=true;
  float *v;
  int idx;
  if(inp[0]<1.1) {/* explicit color */
    copy3f(inp,out);
  } else {
    idx = (int)inp[1];
    /* lookup color */
    if(idx<I->NColor) {
      v=I->Color+3*idx;
      copy3f(v,out);
    } else 
      ok=false;
  }
  return(ok);
}

int GadgetSetFetchNormal(GadgetSet *I,float *inp,float *out)
{
  int ok=true;
  int idx;
  float *v;
  if(inp[0]<1.1) {/* explicit normal */
    copy3f(inp,out);
  } else {
    idx = (int)inp[1];
    /* lookup normal */
    if(idx<I->NNormal) {
      v=I->Normal+3*idx;
      copy3f(v,out);
    } else 
      ok=false;
  }
  return(ok);
}


int GadgetSetFromPyList(PyObject *list,GadgetSet **cs)
{
  int ok = true;
  *cs=NULL;
#if TO_DO
  GadgetSet *I = NULL;
  if(*cs) {
    GadgetSetFree(*cs);
    *cs=NULL;
  }

  if(list==Py_None) { /* allow None for CSet */
    *cs = NULL;
  } else {
  
    if(ok) I=GadgetSetNew();
    if(ok) ok = (I!=NULL);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->NCoord);
    if(ok) ok = PConvPyListToFloatVLA(PyList_GetItem(list,1),&I->Coord);
    if(!ok) {
      if(I)
        GadgetSetFree(I);
    } else {
      *cs = I;
    }
  }
#endif
  return(ok);
}

PyObject *GadgetSetAsPyList(GadgetSet *I)
{
  PyObject *result = NULL;
#if TO_DO

  if(I) {
    result = PyList_New(2);
    
    PyList_SetItem(result,0,PyInt_FromLong(I->NCoord));
    PyList_SetItem(result,1,PConvFloatArrayToPyList(I->Coord,I->NCoord*3));
    /* TODO setting ... */

  }
  return(PConvAutoNone(result));
#endif
  return(PConvAutoNone(result));
}

/*========================================================================*/
int GadgetSetGetExtent(GadgetSet *I,float *mn,float *mx)
{
  float *v;
  int a;
  v = I->Coord;
  for(a=0;a<I->NCoord;a++) {
    min3f(v,mn,mn);
    max3f(v,mx,mx);
    v+=3;
  }
  return(I->NCoord);
}
/*========================================================================*/
void GadgetSetInvalidateRep(GadgetSet *I,int type,int level)
{
#if TO_DO
  int a;
  PRINTFD(FB_GadgetSet)
    " GadgetSetInvalidateRep: entered.\n"
    ENDFD;
  if(type>=0) {
	 if(type<I->NRep)	{
		SceneChanged();		
		if(I->Rep[type]) {
		  I->Rep[type]->fFree(I->Rep[type]);
		  I->Rep[type] = NULL;
		}
	 }
  } else {
	 for(a=0;a<I->NRep;a++)	{
		SceneChanged();
		if(I->Rep[a]) {
		  switch(level) {
		  case cRepInvColor:
			 if(I->Rep[a]->fRecolor) {
				I->Rep[a]->fInvalidate(I->Rep[a],(struct CoordSet*)I,level);
			 } else {
				I->Rep[a]->fFree(I->Rep[a]);
				I->Rep[a] = NULL;
			 }
			 break;
		  default:
			 I->Rep[a]->fFree(I->Rep[a]);
			 I->Rep[a] = NULL;
			 break;
		}
	 }
  }
  }
#endif
}
/*========================================================================*/
void GadgetSetUpdate(GadgetSet *I)
{
  CGO *cgo = NULL,*font_cgo = NULL;
  
  int est;

  if(I->StdCGO) {
    CGOFree(I->StdCGO);
    I->StdCGO=NULL;
  }

  if(I->RayCGO) {
    CGOFree(I->RayCGO);  
    I->RayCGO = NULL;
  }

  if(I->PickShapeCGO) {
    I->PickCGO = CGOProcessShape(I->PickShapeCGO,I,I->PickCGO);
  }

  if(I->ShapeCGO) {
    font_cgo = CGOProcessShape(I->ShapeCGO,I,NULL);
    est=CGOCheckForText(font_cgo);
    if(est) {
      /* assume we've already preloaded fonts... */
      cgo = CGODrawText(font_cgo,est,NULL);
      CGOFree(font_cgo);
    } else {
      cgo = font_cgo;
      font_cgo = NULL;
    }
  }
  if(cgo) {
    est=CGOCheckComplex(cgo);
    if(est) {
      I->RayCGO=cgo;
      I->StdCGO=CGOSimplify(cgo,est);
    } else 
      I->StdCGO=cgo;
  }
}
/*========================================================================*/
void GadgetSetRender(GadgetSet *I,CRay *ray,Pickable **pick,int pass)
{

  float *color;

  color = ColorGet(I->Obj->Obj.Color);

  if(!pass) {
    if(ray) {    
      if(I->RayCGO)
        CGORenderRay(I->RayCGO,ray,color,I->Obj->Obj.Setting,NULL);
      else
        CGORenderRay(I->StdCGO,ray,color,I->Obj->Obj.Setting,NULL);
    } else if(pick&&PMGUI) {
      if(I->PickCGO) {
        CGORenderGLPickable(I->PickCGO,pick,(void*)I->Obj,
                            I->Obj->Obj.Setting,NULL);
      }
    } else if(PMGUI) {
      if(I->StdCGO) {
        /*CGORenderGL(I->PickCGO,color,I->Obj->Obj.Setting,NULL);*/
        CGORenderGL(I->StdCGO,color,I->Obj->Obj.Setting,NULL);
      }
    }
  }
  }
/*========================================================================*/
GadgetSet *GadgetSetNew(void)
{
  OOAlloc(GadgetSet);

  I->fFree=GadgetSetFree;
  I->fRender=GadgetSetRender;
  I->fUpdate=GadgetSetUpdate;
  I->fInvalidateRep=GadgetSetInvalidateRep;
  I->NCoord=0;
  I->NColor=0;
  I->NNormal = 0;
  I->Coord = NULL;
  I->Normal = NULL;
  I->Color = NULL;
  I->Setting=NULL;
  I->PickCGO=NULL;
  I->StdCGO=NULL;
  I->RayCGO=NULL;
  I->ShapeCGO=NULL;
  I->PickShapeCGO=NULL;

  return(I);
}
/*========================================================================*/
void GadgetSetStrip(GadgetSet *I)
{
}
/*========================================================================*/
void GadgetSetFree(GadgetSet *I)
{
  if(I) 
	 {
      CGOFree(I->PickCGO);
      CGOFree(I->PickShapeCGO);
      CGOFree(I->StdCGO);
      CGOFree(I->RayCGO);
      CGOFree(I->ShapeCGO);
      VLAFreeP(I->Coord);
      VLAFreeP(I->Normal);
      VLAFreeP(I->Color);
      OOFreeP(I);
	 }
}
