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

#include"OOMac.h"
#include"RepWireBond.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Setting.h"

typedef struct RepWireBond {
  Rep R;
  float *V,*VP;
  Pickable *P;
  int N,NP;
} RepWireBond;

#include"ObjectMolecule.h"

void RepWireBondRender(RepWireBond *I,CRay *ray,Pickable **pick);
void RepWireBondFree(RepWireBond *I);

void RepWireBondFree(RepWireBond *I)
{
  FreeP(I->VP);
  FreeP(I->V);
  RepFree(&I->R);
  OOFreeP(I);
}

int flip;

void RepWireBondRender(RepWireBond *I,CRay *ray,Pickable **pick)
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
		v+=9;
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

	 }
	 glEnd();

	 (*pick)[0].index = i;
  } else if(PMGUI) {

	 
	 v=I->V;
	 c=I->N;
	 
	 glBegin(GL_LINES);	 
	 SceneResetNormal(true);
	 while(c--) {
		glColor3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
	 }
	 glEnd();
  }
}

Rep *RepWireBondNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,a2,*b,c1,c2,s1,s2,b1,b2;
  int half_bonds;
  float *v,*v0,*v1,*v2,h[3];
  int visFlag;
  OOAlloc(RepWireBond);
  
  obj = cs->Obj;

  visFlag=false;
  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = *(b++);
      b2 = *(b++);
      b++;
      if(obj->AtomInfo[b1].visRep[cRepLine]||
         obj->AtomInfo[b2].visRep[cRepLine]) {
        visFlag=true;
        break;
      }
    }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  RepInit(&I->R);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepWireBondRender;
  I->R.fFree=(void (*)(struct Rep *))RepWireBondFree;

  half_bonds = SettingGet(cSetting_half_bonds);

  I->N=0;
  I->NP=0;
  I->V=NULL;
  I->VP=NULL;
  I->R.P=NULL;
  I->R.fRecolor=NULL;

  if(obj->NBond) {
	 I->V=(float*)mmalloc(sizeof(float)*obj->NBond*9*3);
	 ErrChkPtr(I->V);
	 
	 
	 v=I->V;
	 b=obj->Bond;
	 for(a=0;a<obj->NBond;a++)
		{
		  b1 = *(b++);
		  b2 = *(b++);
        b++;
		  a1=cs->AtmToIdx[b1];
		  a2=cs->AtmToIdx[b2];
		  
		  if((a1>=0)&&(a2>=00))
			 {
				s1=obj->AtomInfo[b1].visRep[cRepLine];
				s2=obj->AtomInfo[b2].visRep[cRepLine];

				if(!(s1&&s2))
              if(!half_bonds) {
                s1 = 0;
                s2 = 0;
              }

				if(s1||s2)
				  {	
					 c1=*(cs->Color+a1);
					 c2=*(cs->Color+a2);
					 
					 v1 = cs->Coord+3*a1;
					 v2 = cs->Coord+3*a2;
					 
					 if((c1==c2)&&s1&&s2) {
						
						I->N++;
						
						v0 = ColorGet(c1);
						
						*(v++)=*(v0++);
						*(v++)=*(v0++);
						*(v++)=*(v0++);
						
						*(v++)=*(v1++);
						*(v++)=*(v1++);
						*(v++)=*(v1++);
						
						*(v++)=*(v2++);
						*(v++)=*(v2++);
						*(v++)=*(v2++);
						
					 } else {
						
						h[0]=(v1[0]+v2[0])/2;
						h[1]=(v1[1]+v2[1])/2;
						h[2]=(v1[2]+v2[2])/2;
						
						if(s1)
						  {
							 I->N++;
							 
							 v0 = ColorGet(c1);
							 
							 *(v++)=*(v0++);
							 *(v++)=*(v0++);
							 *(v++)=*(v0++);
							 
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
						  }
						if(s2)
						  {
							 I->N++;
							 v0 = ColorGet(c2);
							 
							 *(v++)=*(v0++);
							 *(v++)=*(v0++);
							 *(v++)=*(v0++);
							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
							 
							 *(v++)=*(v2++);
							 *(v++)=*(v2++);
							 *(v++)=*(v2++);
						  }
					 }
				  }
			 }
		}

	 I->V = Realloc(I->V,float,I->N*9);

	 /* now create pickable verson */

	 if(SettingGet(cSetting_pickable)) {
		I->VP=(float*)mmalloc(sizeof(float)*obj->NBond*6*2);
		ErrChkPtr(I->VP);
		
		I->R.P=Alloc(Pickable,2*obj->NBond+1);
		ErrChkPtr(I->R.P);
		
		v=I->VP;
		b=obj->Bond;
		for(a=0;a<obj->NBond;a++)
		  {
			 b1 = *(b++);
			 b2 = *(b++);
			 b++;
			 a1=cs->AtmToIdx[b1];
			 a2=cs->AtmToIdx[b2];
			 
			 if((a1>=0)&&(a2>=0))
				{
				  s1=obj->AtomInfo[b1].visRep[cRepLine];
				  s2=obj->AtomInfo[b2].visRep[cRepLine];
				  
				  if(s1||s2)
					 {	
						v1 = cs->Coord+3*a1;
						v2 = cs->Coord+3*a2;
						
						h[0]=(v1[0]+v2[0])/2;
						h[1]=(v1[1]+v2[1])/2;
						h[2]=(v1[2]+v2[2])/2;
						
						if(s1)
						  {
							 I->NP++;
							 
							 I->R.P[I->NP].ptr = (void*)obj;
							 I->R.P[I->NP].index = *(b-3);
							 
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 *(v++)=*(v1++);
							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
						  }
						if(s2)
						  {
							 I->NP++;
							 I->R.P[I->NP].ptr = (void*)obj;
							 I->R.P[I->NP].index = *(b-2);
							 
							 
							 *(v++)=h[0];
							 *(v++)=h[1];
							 *(v++)=h[2];
							 
							 *(v++)=*(v2++);
							 *(v++)=*(v2++);
							 *(v++)=*(v2++);
						  }
					 }
				}
		  }
		I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
		I->R.P[0].index = I->NP;
		I->VP = Realloc(I->VP,float,I->NP*9);
	 }
  }
  return((void*)(struct Rep*)I);
}


