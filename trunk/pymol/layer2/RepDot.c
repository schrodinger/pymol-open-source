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


#include<GL/gl.h>
#include"Base.h"
#include"OOMac.h"
#include"RepDot.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"ObjectMolecule.h"

void RepDotRender(RepDot *I,CRay *ray,Pickable **pick);
void RepDotFree(RepDot *I);

void RepDotInit(void)
{
}

void RepDotFree(RepDot *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->T);
  FreeP(I->F);
  FreeP(I->VN);
  FreeP(I->A);
  OOFreeP(I);
}

void RepDotRender(RepDot *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  int cc=0;

  if(ray) {
	 glBegin(GL_POINTS);
	 while(c--)
		{
		  if(!cc) /* load up the current vertex color */
			 {
				cc=(*(v++));
				ray->fColor3fv(ray,v);
				v+=3;
			 }
		  v+=3;
		  ray->fSphere3fv(ray,v,I->dotSize);
		  glVertex3fv(v);
		  v+=3;
		  cc--;
		}
	 glEnd();

	 /*	 v=I->VC;
	 c=I->NC;
	 while(c--) {
		ray->fColor3fv(ray,v);
		v+=3;
		ray->fSphere3fv(ray,v,*(v+3));
		v+=4;
		}*/

  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 glBegin(GL_POINTS);
	 while(c--)
		{
		  if(!cc) /* load up the current vertex color */
			 {
				cc=(*(v++));
				glColor3fv(v);
				v+=3;
			 }
		  glNormal3fv(v);
		  v+=3;
		  glVertex3fv(v);
		  v+=3;
		  cc--;
		}
	 glEnd();
  }
}

Rep *RepDotNew(CoordSet *cs,int mode)
{
  ObjectMolecule *obj;
  int a,b,a1,a2,flag,h,k,l,i,j,c1;
  float *v,*v0,*vc,vdw,*vn;
  float *aa=NULL;
  int *tp=NULL;
  int *tf=NULL;
  float *countPtr = NULL;
  int colorCnt,lastColor;
  Vector3f v1;
  MapType *map;
  SphereRec *sp = Sphere0; 
  int ds;
  float max_vdw = MAX_VDW;
  float solv_rad=0.0;
  int inclH = true;
  int cullByFlag = false;
  int visFlag;
  OOAlloc(RepDot);

  obj = cs->Obj;
  visFlag=false;
  for(a=0;a<cs->NIndex;a++) {
	 if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepDot])
		{
		  visFlag=true;
		  break;
		}
  }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  RepInit(&I->R);
  I->dotSize = SettingGet(cSetting_dot_size);
  cullByFlag = SettingGet(cSetting_trim_dots);
  inclH = SettingGet(cSetting_dot_hydrogens);

  I->A=NULL;
  I->T=NULL;
  I->F=NULL;
  I->V=NULL;
  I->VC=NULL;
  I->VN=NULL;
  I->R.fRecolor=NULL;

  if(SettingGet(cSetting_dot_surface)>0.0) {
	 solv_rad = SettingGet(cSetting_solvent_radius);
	 max_vdw+=solv_rad;
  }
 
 /* get current dot sampling */
  ds = (int)SettingGet(cSetting_dot_density);
  if(ds<0) ds=0;
  switch(ds) {
  case 0: sp=Sphere0; break;
  case 1: sp=Sphere1; break;
  case 2: sp=Sphere2; break;
  default: sp=Sphere3; break;
  }

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepDotRender;
  I->R.fFree=(void (*)(struct Rep *))RepDotFree;
  
  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*sp->nDot*10);
  ErrChkPtr(I->V);

  if(mode==cRepDotAreaType) {
	 I->A=Alloc(float,cs->NIndex*sp->nDot);
	 I->T=Alloc(int,cs->NIndex*sp->nDot);
	 I->F=Alloc(int,cs->NIndex*sp->nDot);
	 I->VN=Alloc(float,cs->NIndex*sp->nDot*3);
	 aa=I->A;
	 tp=I->T;
	 tf=I->F;
	 inclH=true;
	 cullByFlag=true;
  }
  vn=I->VN;

  I->N=0;
  lastColor=-1;
  colorCnt=0;
  map=MapNew(max_vdw,cs->Coord,cs->NIndex,NULL);
  v=I->V;
  if(map)
	 {
		MapSetupExpress(map);
		for(a=0;a<cs->NIndex;a++)
		  {
			 a1 = cs->IdxToAtm[a];
			 if(obj->AtomInfo[a1].visRep[cRepDot])
				if((inclH||(obj->AtomInfo[a1].name[0]!='H'))&&
					((!cullByFlag)||(!(obj->AtomInfo[a1].customFlag&0x2)))) { /* ignore if the "2" bit is set */
				  c1=*(cs->Color+a);
				  v0 = cs->Coord+3*a;
				  vdw = cs->Obj->AtomInfo[a1].vdw+solv_rad;
				  for(b=0;b<sp->nDot;b++)
					 {
						v1[0]=v0[0]+vdw*sp->dot[b].v[0];
						v1[1]=v0[1]+vdw*sp->dot[b].v[1];
						v1[2]=v0[2]+vdw*sp->dot[b].v[2];
						
						MapLocus(map,v1,&h,&k,&l);

						flag=true;

						i=*(MapEStart(map,h,k,l));
						if(i) {
						  j=map->EList[i++];
						  while(j>=0) {
							 a2 = cs->IdxToAtm[j];
							 if((inclH||(obj->AtomInfo[a2].name[0]!='H'))&&
								 ((!cullByFlag)||(!(obj->AtomInfo[a2].customFlag&0x2))))  /* ignore if the "2" bit is set */
								if(j!=a)
								  if(within3f(cs->Coord+3*j,v1,
												  cs->Obj->AtomInfo[a2].vdw+solv_rad))
									 {
										flag=false;
										break;
									 }
							 j=map->EList[i++];
						  }
						}
						if(flag)
						  {
							 switch(mode) {
							 case cRepDotNormal:
								if(lastColor!=c1) /* new color */
								  {
									 if(countPtr) /* after first pass */
										*countPtr=(float)colorCnt; /* save count */
									 colorCnt=1;
									 countPtr=v++;
									 vc = ColorGet(c1); /* save new color */
									 lastColor=c1;
									 *(v++)=*(vc++);
									 *(v++)=*(vc++);
									 *(v++)=*(vc++);
								  }
								else 
								  colorCnt++;
								*(v++)=sp->dot[b].v[0];
								*(v++)=sp->dot[b].v[1];
								*(v++)=sp->dot[b].v[2];
								*(v++)=v1[0];
								*(v++)=v1[1];
								*(v++)=v1[2];
								I->N++;
								break;
							 case cRepDotAreaType:
								*(v++)=v1[0];
								*(v++)=v1[1];
								*(v++)=v1[2];
								*(aa++)=vdw*vdw*sp->dot[b].area;
								*(tp++)=cs->Obj->AtomInfo[a1].customType;
								*(tf++)=cs->Obj->AtomInfo[a1].customFlag;
								*(vn++)=sp->dot[b].v[0];
								*(vn++)=sp->dot[b].v[1];
								*(vn++)=sp->dot[b].v[2];
								I->N++;
								break;
							 }
						  }
					 }
				}
		  }
		if(countPtr) *countPtr=(float)colorCnt; /* save count */
		MapFree(map);
	 }
  
  I->V = Realloc(I->V,float,(v-I->V));
  
  if(mode==cRepDotAreaType) {
	 I->A = Realloc(I->A,float,(aa-I->A));
	 I->T= Realloc(I->T,int,(tp-I->T));
	 I->F= Realloc(I->F,int,(tf-I->F));
	 I->VN= Realloc(I->VN,float,(vn-I->VN));
  }
  return((void*)(struct Rep*)I);
}




