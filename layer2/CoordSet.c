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
#include"CoordSet.h"
#include"Color.h"

#include"RepWireBond.h"
#include"RepCylBond.h"
#include"RepDot.h"
#include"RepMesh.h"
#include"RepSphere.h"
#include"RepRibbon.h"
#include"RepSurface.h"

void CoordSetUpdate(CoordSet *I);

void CoordSetFree(CoordSet *I);
void CoordSetRender(CoordSet *I,CRay *ray,Pickable **pick);
void CoordSetEnumIndices(CoordSet *I);
void CoordSetStrip(CoordSet *I);
void CoordSetInvalidateRep(CoordSet *I,int type,int level);
void CoordSetExtendIndices(CoordSet *I,int nAtom);
void CoordSetAppendIndices(CoordSet *I,int offset);


/*========================================================================*/
static  char sATOM[]="ATOM  ";
static  char sHETATM[]="HETATM";
/*========================================================================*/
void CoordSetRealToFrac(CoordSet *I,CCrystal *cryst)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform33f3f(cryst->RealToFrac,v,v);
    v+=3;
  }
}
/*========================================================================*/
void CoordSetTransform44f(CoordSet *I,float *mat)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform44f3f(mat,v,v);
    v+=3;
  }  
}
/*========================================================================*/
void CoordSetGetAverage(CoordSet *I,float *v0)
{
  int a;
  float *v;
  double accum[3];
  if(I->NIndex) {
    accum[0]=*(v++);
    accum[1]=*(v++);
    accum[2]=*(v++);
    v=I->Coord;
    for(a=1;a<I->NIndex;a++) {
      accum[0]+=*(v++);
      accum[1]+=*(v++);
      accum[2]+=*(v++);
    }
    v0[0]=accum[0]/I->NIndex;
    v0[1]=accum[1]/I->NIndex;
    v0[2]=accum[2]/I->NIndex;
  }
}
/*========================================================================*/
void CoordSetFracToReal(CoordSet *I,CCrystal *cryst)
{
  int a;
  float *v;
  v=I->Coord;
  for(a=0;a<I->NIndex;a++) {
    transform33f3f(cryst->FracToReal,v,v);
    v+=3;
  }
}
/*========================================================================*/
void CoordSetAtomToPDBStrVLA(char **charVLA,int *c,AtomInfoType *ai,float *v,int cnt)
{
  char *aType;
  AtomName name;
  ResIdent resi;
  int rl;

  if(ai->hetatm)
	aType=sHETATM;
  else
	aType=sATOM;

  strcpy(resi,ai->resi);
  rl = strlen(resi)-1;
  if(rl>=0)
    if((resi[rl]>='0')&&(resi[rl]<='9')) {
        resi[rl+1]=' ';
        resi[rl+2]=0;
    }
  VLACheck(*charVLA,char,(*c)+1000);  

  if(strlen(ai->name)<4)
	{
	  name[0]=' ';	
	  strcpy(name+1,ai->name);
	} else {
	  strcpy(name,ai->name);
	}
  (*c)+=sprintf((*charVLA)+(*c),"%6s%5i %-4s %3s %1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
				aType,cnt+1,name,ai->resn,
				ai->chain,resi,*v,*(v+1),*(v+2),ai->q,ai->b,ai->segi);
  
}

/*========================================================================*/
void CoordSetInvalidateRep(CoordSet *I,int type,int level)
{
  int a;
  if(level>=cRepInvColor) 
	 VLAFreeP(I->Color);
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
void CoordSetUpdate(CoordSet *I)
{
  int a;
  int i;
  if(!I->Color) /* colors invalidated */
	 {
		I->Color=VLAlloc(int,I->NIndex);
		if(I->Color) {
		  for(a=0;a<I->Obj->NAtom;a++)
			 {
				i=I->AtmToIdx[a];
				if(i>=0) 
				  I->Color[i]=I->Obj->AtomInfo[a].color;
			 }
		}
	 }
  OrthoBusyFast(0,I->NRep);
  if(!I->Rep[cRepLine]) {
    I->Rep[cRepLine]=RepWireBondNew(I);
    SceneDirty();
  }
  OrthoBusyFast(1,I->NRep);
  if(!I->Rep[cRepCyl]) {
    I->Rep[cRepCyl]=RepCylBondNew(I);
    SceneDirty();
  }
  OrthoBusyFast(2,I->NRep);
  if(!I->Rep[cRepDot]) {
    I->Rep[cRepDot]=RepDotNew(I,0);
    SceneDirty();
  }
  OrthoBusyFast(3,I->NRep);
  if(!I->Rep[cRepMesh]) {
    I->Rep[cRepMesh]=RepMeshNew(I);
    SceneDirty();
  } else {
    I->Rep[cRepMesh]->fUpdate(I->Rep[cRepMesh],I);
    SceneDirty();
  }
  OrthoBusyFast(4,I->NRep);
  if(!I->Rep[cRepSphere]) {
    I->Rep[cRepSphere]=RepSphereNew(I);
    SceneDirty();
  }
  OrthoBusyFast(5,I->NRep);
  if(!I->Rep[cRepRibbon]) {
    I->Rep[cRepRibbon]=RepRibbonNew(I);
    SceneDirty();
  }
  OrthoBusyFast(6,I->NRep);
  if(!I->Rep[cRepSurface]) {
    I->Rep[cRepSurface]=RepSurfaceNew(I);
    SceneDirty();
  } else {
    I->Rep[cRepSurface]->fUpdate(I->Rep[cRepSurface],I);
    SceneDirty();
  }
  OrthoBusyFast(1,1);
}
/*========================================================================*/
void CoordSetRender(CoordSet *I,CRay *ray,Pickable **pick)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a]) 
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
CoordSet *CoordSetNew(void)
{
  int a;
  OOAlloc(CoordSet);

  I->fFree=CoordSetFree;
  I->fRender=CoordSetRender;
  I->fUpdate=CoordSetUpdate;
  I->fEnumIndices=CoordSetEnumIndices;
  I->fExtendIndices=CoordSetExtendIndices;
  I->fAppendIndices=CoordSetAppendIndices;
  I->fInvalidateRep=CoordSetInvalidateRep;
  I->NIndex=0;
  I->Coord = NULL;
  I->Color = NULL;
  I->AtmToIdx = NULL;
  I->IdxToAtm = NULL;
  I->TmpBond = NULL;
  I->Rep=VLAlloc(Rep*,10);
  I->NRep=cRepCnt;
  I->TmpSymmetry = NULL;
  for(a=0;a<I->NRep;a++)
	 I->Rep[a] = NULL;
  return(I);
}
/*========================================================================*/
CoordSet *CoordSetCopy(CoordSet *cs)
{
  int a;
  int nAtom;
  float *v0,*v1;
  int *i0,*i1;
  OOAlloc(CoordSet);

  (*I)=(*cs);
  I->Coord = VLAlloc(float,I->NIndex*3);
  v0=I->Coord;
  v1=cs->Coord;
  for(a=0;a<I->NIndex;a++) {
    *(v0++)=*(v1++);
    *(v0++)=*(v1++);
    *(v0++)=*(v1++);
  }

  nAtom = cs->Obj->NAtom;
  I->AtmToIdx = Alloc(int,nAtom);
  i0=I->AtmToIdx;
  i1=cs->AtmToIdx;
  for(a=0;a<nAtom;a++)
    *(i0++)=*(i1++);

  I->IdxToAtm = Alloc(int,I->NIndex);
  i0=I->IdxToAtm;
  i1=cs->IdxToAtm;
  for(a=0;a<I->NIndex;a++)
    *(i0++)=*(i1++);
  
  I->Rep=VLAlloc(Rep*,I->NRep);
  for(a=0;a<I->NRep;a++)	
    I->Rep[a] = NULL;

  I->TmpBond=NULL;
  I->Color=NULL;

  return(I);
}
/*========================================================================*/
void CoordSetExtendIndices(CoordSet *I,int nAtom)
{
  int a;
  if(I->NAtIndex<nAtom)
	 {
		if(I->AtmToIdx) {
		  I->AtmToIdx = Realloc(I->AtmToIdx,int,nAtom);
		  ErrChkPtr(I->AtmToIdx);
		  for(a=I->NAtIndex;a<nAtom;a++)
			 I->AtmToIdx[a]=-1;
		  I->NAtIndex = nAtom;
		} else {
		  I->AtmToIdx = Alloc(int,nAtom);
		  for(a=0;a<nAtom;a++)
			 I->AtmToIdx[a]=-1;
		  I->NAtIndex = nAtom;
		}
	 }
}
/*========================================================================*/
void CoordSetAppendIndices(CoordSet *I,int offset) 
	  /* this is going to get impractical down the road...
		  we need to revise the index system */
{
  int a;
  I->AtmToIdx = Alloc(int,I->NIndex+offset);
  I->IdxToAtm = Alloc(int,I->NIndex);
  ErrChkPtr(I->AtmToIdx);
  ErrChkPtr(I->IdxToAtm);
  for(a=0;a<offset;a++)
	 I->AtmToIdx[a]=-1;
  for(a=0;a<I->NIndex;a++) {
	 I->AtmToIdx[a+offset]=a;
	 I->IdxToAtm[a]=a+offset;
  }
  I->NAtIndex = I->NIndex + offset;
}
/*========================================================================*/
void CoordSetEnumIndices(CoordSet *I)
{
  /* set up for simple case where 1 = 1, etc. */
  int a;
  I->AtmToIdx = Alloc(int,I->NIndex);
  I->IdxToAtm = Alloc(int,I->NIndex);
  ErrChkPtr(I->AtmToIdx);
  ErrChkPtr(I->IdxToAtm);
  for(a=0;a<I->NIndex;a++)
	 {
	 I->AtmToIdx[a]=a;
	 I->IdxToAtm[a]=a;
	 }
  I->NAtIndex = I->NIndex;
}
/*========================================================================*/
void CoordSetStrip(CoordSet *I)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a])
		I->Rep[a]->fFree(I->Rep[a]);
  I->NRep=0;
}
/*========================================================================*/
void CoordSetFree(CoordSet *I)
{
  int a;
  for(a=0;a<I->NRep;a++)
	 if(I->Rep[a])
		I->Rep[a]->fFree(I->Rep[a]);
  if(I) 
	 {
      
	 FreeP(I->AtmToIdx);
	 FreeP(I->IdxToAtm);
	 VLAFreeP(I->Color);
	 VLAFreeP(I->Coord);
	 VLAFreeP(I->Rep);
	 VLAFreeP(I->TmpBond);
    if(I->TmpSymmetry) SymmetryFree(I->TmpSymmetry);
	 OOFreeP(I);
	 }
}


