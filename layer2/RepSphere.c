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
#include"RepSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"

typedef struct RepSphere {
  Rep R;
  float *V;
  float *VC;
  SphereRec *SP;
  int N,NC;

} RepSphere;

#include"ObjectMolecule.h"

void RepSphereRender(RepSphere *I,CRay *ray,Pickable **pick);
void RepSphereFree(RepSphere *I);

void RepSphereInit(void)
{
}

void RepSphereFree(RepSphere *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  OOFreeP(I);
}

void RepSphereRender(RepSphere *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  int cc=0;
  int a;
  SphereRec *sp;

  if(ray) {
	 v=I->VC;
	 c=I->NC;
	 while(c--) {
		ray->fColor3fv(ray,v);
		v+=3;
		ray->fSphere3fv(ray,v,*(v+3));
		v+=4;
	 }
  } else if(pick) {
  } else {
	 sp=I->SP;

	 while(c--)
		{
		  glColor3fv(v);
		  v+=3;
		  for(a=0;a<sp->NStrip;a++) {
			 glBegin(GL_TRIANGLE_STRIP);
			 cc=sp->StripLen[a];
			 while(cc--) {
				glNormal3fv(v);
				v+=3;
				glVertex3fv(v);
				v+=3;
			 }
			 glEnd();
		  }
		}
  }
}

Rep *RepSphereNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,b,c,a1,c1;
  float *v,*v0,*vc,vdw;
  int *q, *s;
  SphereRec *sp = Sphere0; 
  int ds;

  OOAlloc(RepSphere);

 /* get current dot sampling */
  ds = (int)SettingGet(cSetting_dot_density);
  ds=1;
  if(ds<0) ds=0;
  switch(ds) {
  case 0: sp=Sphere0; break;
  case 1: sp=Sphere1; break;
  case 2: sp=Sphere2; break;
  default: sp=Sphere3; break;
  }

  obj = cs->Obj;
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepSphereRender;
  I->R.fFree=(void (*)(struct Rep *))RepSphereFree;
  
  I->VC=(float*)mmalloc(sizeof(float)*cs->NIndex*7);
  ErrChkPtr(I->VC);
  I->NC=0;

  v=I->VC; 
  for(a=0;a<cs->NIndex;a++)
	 {
		a1 = cs->IdxToAtm[a];
		if(obj->AtomInfo[a1].visRep[cRepSphere])
		  {
			 I->NC++;
			 c1=*(cs->Color+a);
			 vc = ColorGet(c1); /* save new color */
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 v0 = cs->Coord+3*a;			 
			 *(v++)=*(v0++);
			 *(v++)=*(v0++);
			 *(v++)=*(v0++);
			 *(v++)=obj->AtomInfo[a1].vdw;
		  }
	 }

  if(I->NC) 
	 I->VC=(float*)mrealloc(I->VC,sizeof(float)*(v-I->VC));
  else
	 I->VC=(float*)mrealloc(I->VC,1);

  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*(3+sp->NVertTot*6));
  ErrChkPtr(I->V);

  I->N=0;
  I->SP=sp;
  v=I->V;

  for(a=0;a<cs->NIndex;a++)
	 {
		a1 = cs->IdxToAtm[a];
		if(obj->AtomInfo[a1].visRep[cRepSphere])
		  {
			 I->N++;
			 c1=*(cs->Color+a);
			 v0 = cs->Coord+3*a;
			 vdw = cs->Obj->AtomInfo[a1].vdw;
			 q=sp->Sequence;
			 s=sp->StripLen;
			 vc = ColorGet(c1);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 for(b=0;b<sp->NStrip;b++)
				{
				  for(c=0;c<(*s);c++)
					 {
						*(v++)=sp->dot[*q].v[0]; /* normal */
						*(v++)=sp->dot[*q].v[1];
						*(v++)=sp->dot[*q].v[2];
						*(v++)=v0[0]+vdw*sp->dot[*q].v[0]; /* point */
						*(v++)=v0[1]+vdw*sp->dot[*q].v[1];
						*(v++)=v0[2]+vdw*sp->dot[*q].v[2];
						q++;
					 }
				  s++;
				}
		  }
	 }

  /*TODO: NEED TO SHRINK POINTERS HERE*/

  if(I->N) 
	 I->V=(float*)mrealloc(I->V,sizeof(float)*(v-I->V));
  else
	 I->V=(float*)mrealloc(I->V,1);
  
  return((void*)(struct Rep*)I);
}


