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

#include"os_gl.h"

#include"Base.h"
#include"OOMac.h"
#include"RepNonbondedSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Setting.h"
#include"main.h"

typedef struct RepNonbondedSphere {
  Rep R;
  float *V;
  float *VC;
  SphereRec *SP;
  int N,NC;
} RepNonbondedSphere;

#include"ObjectMolecule.h"

void RepNonbondedSphereRender(RepNonbondedSphere *I,CRay *ray,Pickable **pick);
void RepNonbondedSphereFree(RepNonbondedSphere *I);

void RepNonbondedSphereInit(void)
{
}

void RepNonbondedSphereFree(RepNonbondedSphere *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  OOFreeP(I);
}

void RepNonbondedSphereRender(RepNonbondedSphere *I,CRay *ray,Pickable **pick)
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
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
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

Rep *RepNonbondedSphereNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,c,d,c1;
  float *v,*v0,*vc;
  float nonbonded_size;
  int *q, *s;
  SphereRec *sp = Sphere0; 
  int ds;
  int *active=NULL;
  AtomInfoType *ai;
  int nSphere = 0;
  OOAlloc(RepNonbondedSphere);
  obj = cs->Obj;

  active = Alloc(int,cs->NIndex);
  
  for(a=0;a<cs->NIndex;a++) {
    ai = obj->AtomInfo+cs->IdxToAtm[a];
    active[a] =(!ai->bonded) && (ai->visRep[ cRepNonbondedSphere ]);
    if(active[a]) nSphere++;
  }
  if(!nSphere) {
    OOFreeP(I);
    FreeP(active);
    return(NULL); /* skip if no dots are visible */
  }
  
  nonbonded_size = SettingGet(cSetting_nonbonded_size);
  
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

  RepInit(&I->R);
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepNonbondedSphereRender;
  I->R.fFree=(void (*)(struct Rep *))RepNonbondedSphereFree;
  I->R.fRecolor=NULL;
  I->N=0;
  I->NC=0;
  I->V=NULL;
  I->VC=NULL;
  I->SP=NULL;

  /* raytracing primitives */

  I->VC=(float*)mmalloc(sizeof(float)*nSphere*7);
  ErrChkPtr(I->VC);
  I->NC=0;

  v=I->VC; 

  for(a=0;a<cs->NIndex;a++)
	 {
      if(active[a])
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
			 *(v++)=nonbonded_size;
		  }
	 }

  if(I->NC) 
	 I->VC=(float*)mrealloc(I->VC,sizeof(float)*(v-I->VC));
  else
	 I->VC=(float*)mrealloc(I->VC,1);

  I->V=(float*)mmalloc(sizeof(float)*nSphere*(3+sp->NVertTot*6));
  ErrChkPtr(I->V);

  /* rendering primitives */

  I->N=0;
  I->SP=sp;
  v=I->V;

  for(a=0;a<cs->NIndex;a++)
	 {
		if(active[a])
		  {
			 c1=*(cs->Color+a);
			 v0 = cs->Coord+3*a;
			 vc = ColorGet(c1);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);

          q=sp->Sequence;
          s=sp->StripLen;
          
          for(d=0;d<sp->NStrip;d++)
            {
              for(c=0;c<(*s);c++)
                {
                  *(v++)=sp->dot[*q].v[0]; /* normal */
                  *(v++)=sp->dot[*q].v[1];
                  *(v++)=sp->dot[*q].v[2];
                  *(v++)=v0[0]+nonbonded_size*sp->dot[*q].v[0]; /* point */
                  *(v++)=v0[1]+nonbonded_size*sp->dot[*q].v[1];
                  *(v++)=v0[2]+nonbonded_size*sp->dot[*q].v[2];
                  q++;
                }
              s++;
            }
			 I->N++;
		  }
	 }
  
  if(I->N) 
    I->V=(float*)mrealloc(I->V,sizeof(float)*(v-I->V));
  else
    I->V=(float*)mrealloc(I->V,1);

  FreeP(active);
  return((void*)(struct Rep*)I);
}


