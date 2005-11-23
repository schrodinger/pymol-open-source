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
#include"RepAngle.h"
#include"RepDihedral.h"
#include"PConv.h"

void DistSetUpdate(DistSet *I);
void DistSetFree(DistSet *I);
void DistSetStrip(DistSet *I);
void DistSetInvalidateRep(DistSet *I,int type,int level);

int DistSetFromPyList(PyMOLGlobals *G,PyObject *list,DistSet **cs)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  DistSet *I = NULL;
  int ok = true;
  int ll = 0;
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
    if(ok) ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list,1),&I->Coord);

    if(ok && (ll>2)) { /* have angles and torsions too... */
      if(ok) ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list,2),&I->LabelCoord);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,3),&I->NAngleIndex);
      if(ok) ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list,4),&I->AngleCoord);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,5),&I->NDihedralIndex);
      if(ok) ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list,6),&I->DihedralCoord);
    }
    if(ok&&(ll>7)) I->Setting = SettingNewFromPyList(G,PyList_GetItem(list,7)); /* state settings */

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
    result = PyList_New(8);
    
    PyList_SetItem(result,0,PyInt_FromLong(I->NIndex));
    PyList_SetItem(result,1,PConvFloatArrayToPyListNullOkay(I->Coord,I->NIndex*3));
    PyList_SetItem(result,2,PConvFloatArrayToPyListNullOkay(I->LabelCoord,I->NIndex*3));
    PyList_SetItem(result,3,PyInt_FromLong(I->NAngleIndex));
    PyList_SetItem(result,4,PConvFloatArrayToPyListNullOkay(I->AngleCoord,I->NAngleIndex*3));
    PyList_SetItem(result,5,PyInt_FromLong(I->NDihedralIndex));
    PyList_SetItem(result,6,PConvFloatArrayToPyListNullOkay(I->DihedralCoord,I->NDihedralIndex*3));
    PyList_SetItem(result,7,SettingAsPyList(I->Setting));
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

  v = I->AngleCoord;
  for(a=0;a<I->NAngleIndex;a++) {
    min3f(v,mn,mn);
    max3f(v,mx,mx);
    v+=3;
  }

  v = I->DihedralCoord;
  for(a=0;a<I->NDihedralIndex;a++) {
    min3f(v,mn,mn);
    max3f(v,mx,mx);
    v+=3;
  }
  return(I->NIndex+I->NAngleIndex+I->NDihedralIndex);
}
/*========================================================================*/
void DistSetInvalidateRep(DistSet *I,int type,int level)
{
  int a;
  PRINTFD(I->State.G,FB_DistSet)
    " DistSetInvalidateRep: entered.\n"
    ENDFD;
  if(type>=0) {
	 if(type<I->NRep)	{
		SceneChanged(I->State.G);		
		if(I->Rep[type]) {
		  I->Rep[type]->fFree(I->Rep[type]);
		  I->Rep[type] = NULL;
		}
	 }
  } else {
	 for(a=0;a<I->NRep;a++)	{
		SceneChanged(I->State.G);
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

  OrthoBusyFast(I->State.G,0,I->NRep);
  if(!I->Rep[cRepDash]) {
    I->Rep[cRepDash]=RepDistDashNew(I);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepLabel]) {
    I->Rep[cRepLabel]=RepDistLabelNew(I);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepAngle]) {
    I->Rep[cRepAngle]=RepAngleNew(I);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepDihedral]) {
    I->Rep[cRepDihedral]=RepDihedralNew(I);
    SceneInvalidate(I->State.G);
  }
  OrthoBusyFast(I->State.G,1,1);
}
/*========================================================================*/
static void DistSetRender(DistSet *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  int pass = info->pass;
  Picking **pick = info->pick;
  int float_labels = SettingGet_i(I->State.G,
                                  I->Setting,
                                  I->Obj->Obj.Setting,
                                  cSetting_float_labels);
  int a;
  Rep *r;
  for(a=0;a<I->NRep;a++)
    if(I->Rep[a])
      if(I->Obj->Obj.RepVis[a]) {
        r = I->Rep[a];
        if(ray||pick) {
          if(ray) 
            ray->fColor3fv(ray,ColorGet(I->State.G,I->Obj->Obj.Color));
          r->fRender(r,info);
        } else {
          ObjectUseColor((CObject*)I->Obj);
          switch(a) {
          case cRepLabel:
            if(float_labels) {
              if(pass==-1) 
                r->fRender(r,info);
            } else if(pass==0)
              r->fRender(r,info);                  
            break;
          default:
            if(pass==0) {
              r->fRender(r,info);
            }
            break;
          }
        }
      }
}
/*========================================================================*/
DistSet *DistSetNew(PyMOLGlobals *G)
{
  int a;
  OOAlloc(G,DistSet);
  I->State.G=G;
  I->State.Matrix=NULL;
  I->fFree=DistSetFree;
  I->fRender=DistSetRender;
  I->fUpdate=DistSetUpdate;
  I->fInvalidateRep=DistSetInvalidateRep;
  I->NIndex=0;
  I->Coord = NULL;
  I->Rep=VLAlloc(Rep*,cRepCnt);
  I->NRep=cRepCnt;
  I->Setting=NULL;
  I->LabelCoord = NULL;
  I->AngleCoord = NULL;
  I->NAngleIndex = 0;
  I->DihedralCoord = NULL;
  I->NDihedralIndex = 0;
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
	 VLAFreeP(I->AngleCoord);
	 VLAFreeP(I->DihedralCoord);
     VLAFreeP(I->LabelCoord);
	 VLAFreeP(I->Coord);
	 VLAFreeP(I->Rep);
     SettingFreeP(I->Setting);
	 OOFreeP(I);
	 }
}


