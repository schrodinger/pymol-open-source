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
#include<string.h>

#include"Base.h"
#include"Debug.h"
#include"OOMac.h"
#include"Vector.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Setting.h"
#include"Scene.h"
#include"Ray.h"
#include"ObjectDist.h"
#include"Selector.h"

void ObjectDistRender(ObjectDist *I,int frame,CRay *ray,Pickable **pick);
void ObjectDistFree(ObjectDist *I);
void ObjectDistUpdate(ObjectDist *I);
int ObjectDistGetNFrames(ObjectDist *I);

/*========================================================================*/
int ObjectDistGetNFrames(ObjectDist *I)
{
  return I->NDSet;
}
/*========================================================================*/
void ObjectDistUpdate(ObjectDist *I)
{
  int a;
  OrthoBusyPrime();
  for(a=0;a<I->NDSet;a++)
	 if(I->DSet[a]) {	
	   OrthoBusySlow(a,I->NDSet);
	   printf(" ObjectDist: updating state %d of \"%s\".\n" , a+1, I->Obj.Name);
      if(I->DSet[a]->fUpdate)
        I->DSet[a]->fUpdate(I->DSet[a]);
	 }
}
/*========================================================================*/
void ObjectDistInvalidateRep(ObjectDist *I,int rep)
{
  int a;
  for(a=0;a<I->NDSet;a++) 
	 if(I->DSet[a]) {	 
      if(I->DSet[a]->fInvalidateRep)
        I->DSet[a]->fInvalidateRep(I->DSet[a],rep,0);
	 }
}
/*========================================================================*/
void ObjectDistRender(ObjectDist *I,int frame,CRay *ray,Pickable **pick)
{
  int a;
  if(frame<0) {
    for(a=0;a<I->NDSet;a++)
      if(I->DSet[a])
        if(I->DSet[a]->fRender)
          I->DSet[a]->fRender(I->DSet[a],ray,pick);        
  } else if(frame<I->NDSet) {
	 I->CurDSet=frame % I->NDSet;
	 if(I->DSet[I->CurDSet]) {
      if(I->DSet[I->CurDSet]->fRender)
        I->DSet[I->CurDSet]->fRender(I->DSet[I->CurDSet],ray,pick);
	 }
  } else if(I->NDSet==1) { /* if only one coordinate set, assume static */
    if(I->DSet[0]->fRender)
      I->DSet[0]->fRender(I->DSet[0],ray,pick);    
  }
}
/*========================================================================*/
ObjectDist *ObjectDistNew(int sele1,int sele2,int mode,float cutoff)
{
  int a,mn,n;
  OOAlloc(ObjectDist);
  ObjectInit((Object*)I);
  I->Obj.type=cObjectDist;
  I->NAtom=0;
  I->DSet=VLAMalloc(10,sizeof(DistSet*),5,true); /* auto-zero */
  I->NDSet=0;
  I->Obj.fRender=(void (*)(struct Object *, int, CRay *, Pickable **))ObjectDistRender;
  I->Obj.fFree= (void (*)(struct Object *))ObjectDistFree;
  I->Obj.fUpdate= (void (*)(struct Object *)) ObjectDistUpdate;
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectDistGetNFrames;
  I->Obj.fDescribeElement = NULL;
  I->CurDSet=0;
  I->Obj.Color=ColorGetIndex("dash");
  mn = 0;
  SelectorUpdateTable();
  mn = SelectorGetSeleNCSet(sele1);
  n = SelectorGetSeleNCSet(sele2);
  if(mn>n) mn = n;
  if(mn) {
    for(a=0;a<mn;a++)
      {
        VLACheck(I->DSet,DistSet*,a);
        I->DSet[a] = SelectorGetDistSet(sele1,a,sele2,a,mode,cutoff);
        if(I->DSet[a]) {
          I->DSet[a]->Obj = I;
          I->NDSet=a+1;
        }
      }  
  } else 
    OOFreeP(I);

  SceneChanged();
  return(I);
}
/*========================================================================*/
void ObjectDistFree(ObjectDist *I)
{
  int a;
  SceneObjectDel((Object*)I);
  for(a=0;a<I->NDSet;a++)
	 if(I->DSet[a]) {
      if(I->DSet[a]->fFree)
        I->DSet[a]->fFree(I->DSet[a]);
		I->DSet[a]=NULL;
	 }
  VLAFreeP(I->DSet);
  OOFreeP(I);
}

