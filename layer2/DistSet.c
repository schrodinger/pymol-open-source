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
#include"DistSet.h"
#include"Color.h"
#include"RepDistDash.h"
#include"RepDistLabel.h"
#include"PConv.h"

void DistSetUpdate(DistSet *I);
void DistSetFree(DistSet *I);
void DistSetRender(DistSet *I,CRay *ray,Pickable **pick,int pass);
void DistSetStrip(DistSet *I);
void DistSetInvalidateRep(DistSet *I,int type,int level);

int DistSetFromPyList(PyMOLGlobals *G,PyObject *list,DistSet **cs)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  DistSet *I = NULL;
  int ok = true;
  int ll;
  if(*cs) {
    DistSetFree(*cs);
    *cs=NULL;
  }

  if(list==Py_None) { /* allow None for CSet */
    *cs = NULL;
  } else {
  
    if(ok) I=DistSetNew(G);
    if(ok) ok = (I!=NULL);
    if(ok) ok = (list!=NULL);
    if(ok) ok = PyList_Check(list);
    if(ok) ll = PyList_Size(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->NIndex);
    if(ok) ok = PConvPyListToFloatVLA(PyList_GetItem(list,1),&I->Coord);
    if(!ok) {
      if(I)
        DistSetFree(I);
    } else {
      *cs = I;
    }
  }
  return(ok);
#endif
}

PyObject *DistSetAsPyList(DistSet *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(2);
    
    PyList_SetItem(result,0,PyInt_FromLong(I->NIndex));
    PyList_SetItem(result,1,PConvFloatArrayToPyList(I->Coord,I->NIndex*3));
    /* TODO setting ... */

  }
  return(PConvAutoNone(result));
#endif
}

/*========================================================================*/
int DistSetGetExtent(DistSet *I,float *mn,float *mx)
{
  float *v;
  int a;
  v = I->Coord;
  for(a=0;a<I->NIndex;a++) {
    min3f(v,mn,mn);
    max3f(v,mx,mx);
    v+=3;
  }
  return(I->NIndex);
}
/*========================================================================*/
void DistSetInvalidateRep(DistSet *I,int type,int level)
{
  int a;
  PRINTFD(I->G,FB_DistSet)
    " DistSetInvalidateRep: entered.\n"
    ENDFD;
  if(type>=0) {
	 if(type<I->NRep)	{
		SceneChanged(I->G);		
		if(I->Rep[type]) {
		  I->Rep[type]->fFree(I->Rep[type]);
		  I->Rep[type] = NULL;
		}
	 }
  } else {
	 for(a=0;a<I->NRep;a++)	{
		SceneChanged(I->G);
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
}
/*========================================================================*/
void DistSetUpdate(DistSet *I)
{

  OrthoBusyFast(I->G,0,I->NRep);
  if(!I->Rep[cRepDash]) {
    I->Rep[cRepDash]=RepDistDashNew(I);
    SceneDirty(I->G);
  }
  if(!I->Rep[cRepLabel]) {
    I->Rep[cRepLabel]=RepDistLabelNew(I);
    SceneDirty(I->G);
  }
  OrthoBusyFast(I->G,1,1);
}
/*========================================================================*/
void DistSetRender(DistSet *I,CRay *ray,Pickable **pick,int pass)
{
  int a;
  if(!pass) { /* only render on zero/default pass */
    for(a=0;a<I->NRep;a++)
      if(I->Rep[a]) 
        if(I->Obj->Obj.RepVis[a])
          {
            if(!ray) {
              ObjectUseColor((CObject*)I->Obj);
            } else {
              ray->fColor3fv(ray,ColorGet(I->G,I->Obj->Obj.Color));
            }			 
            I->Rep[a]->fRender(I->Rep[a],ray,pick);
        }
  }
}
/*========================================================================*/
DistSet *DistSetNew(PyMOLGlobals *G)
{
  int a;
  OOAlloc(G,DistSet);
  I->G=G;
  I->fFree=DistSetFree;
  I->fRender=DistSetRender;
  I->fUpdate=DistSetUpdate;
  I->fInvalidateRep=DistSetInvalidateRep;
  I->NIndex=0;
  I->Coord = NULL;
  I->Rep=VLAlloc(Rep*,cRepCnt);
  I->NRep=cRepCnt;
  I->Setting=NULL;
  for(a=0;a<I->NRep;a++)
	 I->Rep[a] = NULL;
  return(I);
}
/*========================================================================*/
void DistSetStrip(DistSet *I)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a])
		I->Rep[a]->fFree(I->Rep[a]);
  I->NRep=0;
}
/*========================================================================*/
void DistSetFree(DistSet *I)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a])
		I->Rep[a]->fFree(I->Rep[a]);
  if(I) 
	 {
	 VLAFreeP(I->Coord);
	 VLAFreeP(I->Rep);
	 OOFreeP(I);
	 }
}


