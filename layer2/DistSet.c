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
#include<stdlib.h>
#include<math.h>
#include<stdio.h>

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Scene.h"
#include"DistSet.h"
#include"Color.h"
#include"RepDistDash.h"
#include"RepDistLabel.h"

void DistSetUpdate(DistSet *I);
void DistSetFree(DistSet *I);
void DistSetRender(DistSet *I,CRay *ray,Pickable **pick);
void DistSetStrip(DistSet *I);
void DistSetInvalidateRep(DistSet *I,int type,int level);
/*========================================================================*/
void DistSetInvalidateRep(DistSet *I,int type,int level)
{
  int a;
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
				I->Rep[a]->fInvalidate(I->Rep[a],level);
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

  OrthoBusyFast(0,I->NRep);
  if(!I->Rep[cRepDash]) {
    I->Rep[cRepDash]=RepDistDashNew(I);
    SceneDirty();
  }
  if(!I->Rep[cRepLabel]) {
    I->Rep[cRepLabel]=RepDistLabelNew(I);
    SceneDirty();
  }
  OrthoBusyFast(1,1);
}
/*========================================================================*/
void DistSetRender(DistSet *I,CRay *ray,Pickable **pick)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a]) 
      if(I->Obj->Obj.RepVis[a])
        {
          if(!ray) {
            ObjectUseColor((Object*)I->Obj);
          } else {
            ray->fColor3fv(ray,ColorGet(I->Obj->Obj.Color));
          }			 
          I->Rep[a]->fRender(I->Rep[a],ray,pick);
        }
}
/*========================================================================*/
DistSet *DistSetNew(void)
{
  int a;
  OOAlloc(DistSet);

  I->fFree=DistSetFree;
  I->fRender=DistSetRender;
  I->fUpdate=DistSetUpdate;
  I->fInvalidateRep=DistSetInvalidateRep;
  I->NIndex=0;
  I->Coord = NULL;
  I->Rep=VLAlloc(Rep*,10);
  I->NRep=cRepCnt;

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


