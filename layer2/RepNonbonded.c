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

#include"OOMac.h"
#include"RepNonbonded.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Setting.h"

typedef struct RepNonbonded {
  Rep R;
  float *V,*VP;
  Pickable *P;
  int N,NP;
} RepNonbonded;

#include"ObjectMolecule.h"

void RepNonbondedRender(RepNonbonded *I,CRay *ray,Pickable **pick);
void RepNonbondedFree(RepNonbonded *I);

void RepNonbondedFree(RepNonbonded *I)
{
  FreeP(I->VP);
  FreeP(I->V);
  RepFree(&I->R);
  OOFreeP(I);
}

void RepNonbondedRender(RepNonbonded *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  unsigned int i,j;
  Pickable *p;

  if(ray) {

	 v=I->V;
	 c=I->N;
	 
	 while(c--) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]);*/
      ray->fCylinder3fv(ray,v+3,v+6,0.15,v,v);
      ray->fCylinder3fv(ray,v+9,v+12,0.15,v,v);
      ray->fCylinder3fv(ray,v+15,v+18,0.15,v,v);
		v+=21;
	 }

  } else if(pick&&PMGUI) {
	 
	 i=(*pick)->index;

	 v=I->VP;
	 c=I->NP;
	 p=I->R.P;
	 
	 glBegin(GL_LINES);
	 
	 while(c--) {

		i++;

		if(!(*pick)[0].ptr) {
		  /* pass 1 - low order bits */

		  glColor3ub((i&0xF)<<4,(i&0xF0)|0x8,(i&0xF00)>>4); 
		  VLACheck((*pick),Pickable,i);
		  p++;
		  (*pick)[i] = *p; /* copy object and atom info */
		} else { 
		  /* pass 2 - high order bits */

		  j=i>>12;

		  glColor3ub((j&0xF)<<4,(j&0xF0)|0x8,(j&0xF00)>>4); 

		}			 

		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;

	 }
	 glEnd();

	 (*pick)[0].index = i;
  } else if(PMGUI) {

	 v=I->V;
	 c=I->N;
    if(c) {
      glDisable(GL_LIGHTING);
      glBegin(GL_LINES);	 
      SceneResetNormal(true);
      while(c--) {
        
        glColor3fv(v);
        v+=3;
        
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        
      }
      glEnd();
      glEnable(GL_LIGHTING);
    }
  }
}

Rep *RepNonbondedNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,c1;
  float nonbonded_size;
  float *v,*v0,*v1;
  int *active;
  AtomInfoType *ai;
  int nAtom = 0;
  OOAlloc(RepNonbonded);
  obj = cs->Obj;

  active = Alloc(int,cs->NIndex);
  
  for(a=0;a<cs->NIndex;a++) {
    ai = obj->AtomInfo+cs->IdxToAtm[a];
    active[a] =(!ai->bonded) && (ai->visRep[ cRepNonbonded]);
    if(active[a]) nAtom++;
  }
  if(!nAtom) {
    OOFreeP(I);
    FreeP(active);
    return(NULL); /* skip if no dots are visible */
  }

  nonbonded_size = SettingGet(cSetting_nonbonded_size);
  RepInit(&I->R);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepNonbondedRender;
  I->R.fFree=(void (*)(struct Rep *))RepNonbondedFree;

  I->N=0;
  I->NP=0;
  I->V=NULL;
  I->VP=NULL;
  I->R.P=NULL;
  I->R.fRecolor=NULL;

  I->V=(float*)mmalloc(sizeof(float)*nAtom*21);
  ErrChkPtr(I->V);
  v=I->V;
  for(a=0;a<cs->NIndex;a++) 
    if(active[a]) {
      c1=*(cs->Color+a);
      v0 = ColorGet(c1);
      v1 = cs->Coord+3*a;
      *(v++)=*(v0++);
      *(v++)=*(v0++);
      *(v++)=*(v0++);
      *(v++)=v1[0]-nonbonded_size;
      *(v++)=v1[1];
      *(v++)=v1[2];
      *(v++)=v1[0]+nonbonded_size;
      *(v++)=v1[1];
      *(v++)=v1[2];
      *(v++)=v1[0];
      *(v++)=v1[1]-nonbonded_size;
      *(v++)=v1[2];
      *(v++)=v1[0];
      *(v++)=v1[1]+nonbonded_size;
      *(v++)=v1[2];
      *(v++)=v1[0];
      *(v++)=v1[1];
      *(v++)=v1[2]-nonbonded_size;
      *(v++)=v1[0];
      *(v++)=v1[1];
      *(v++)=v1[2]+nonbonded_size;
      I->N++;
    }
  I->V = Realloc(I->V,float,I->N*21);

  /* now create pickable verson */
  
  if(SettingGet(cSetting_pickable)) {
    I->VP=(float*)mmalloc(sizeof(float)*nAtom*18);
    ErrChkPtr(I->VP);
    
    I->R.P=Alloc(Pickable,cs->NIndex+1);
    ErrChkPtr(I->R.P);
    
    v=I->VP;
    
    for(a=0;a<cs->NIndex;a++) 
      if(active[a]) {
        
        I->NP++;

        a1=cs->IdxToAtm[a];
        
        I->R.P[I->NP].ptr = (void*)obj;
        I->R.P[I->NP].index = a1;
        
        v1 = cs->Coord+3*a;
        
        *(v++)=v1[0]-nonbonded_size;
        *(v++)=v1[1];
        *(v++)=v1[2];
        *(v++)=v1[0]+nonbonded_size;
        *(v++)=v1[1];
        *(v++)=v1[2];
        *(v++)=v1[0];
        *(v++)=v1[1]-nonbonded_size;
        *(v++)=v1[2];
        *(v++)=v1[0];
        *(v++)=v1[1]+nonbonded_size;
        *(v++)=v1[2];
        *(v++)=v1[0];
        *(v++)=v1[1];
        *(v++)=v1[2]-nonbonded_size;
        *(v++)=v1[0];
        *(v++)=v1[1];
        *(v++)=v1[2]+nonbonded_size;
      }
    I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
    I->R.P[0].index = I->NP;
    I->VP = Realloc(I->VP,float,I->NP*21);
  }
  FreeP(active);
  return((void*)(struct Rep*)I);
}

