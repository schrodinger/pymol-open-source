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
#include<math.h>
#include"Base.h"
#include"OOMac.h"
#include"RepRibbon.h"
#include"Color.h"
#include"Setting.h"
#include"Word.h"

typedef struct RepRibbon {
  Rep R;
  float *V;
  float *VC;
  int N,NC;

} RepRibbon;

#include"ObjectMolecule.h"

void RepRibbonRender(RepRibbon *I,CRay *ray,Pickable **pick);
void RepRibbonFree(RepRibbon *I);

void RepRibbonInit(void)
{
}

void RepRibbonFree(RepRibbon *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  OOFreeP(I);
}

void RepRibbonRender(RepRibbon *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;

  if(ray) {
	 v=I->VC;
	 c=I->NC-1;
	 if(c>0)
		while(c--) {
		  ray->fCylinder3fv(ray,v+4,v+7,*(v+3),v,v);
		  v+=10;
		}
  } else if(pick) {
  } else {
	 glBegin(GL_LINES);
	 while(c--)
		{
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
Rep *RepRibbonNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,b,a1,a2,c1,c2,*i,*s,*at,*seg,nAt,*atp;
  float *v,*v0,*v1,*v2,*v3;
  float *pv=NULL;
  float *dv=NULL;
  float *nv=NULL;
  float *tv=NULL;
  float *vc=NULL;
  float f0,f1,f2,f3,f4;
  float *d;
  float *dl=NULL;
  int nSeg;
  int sampling;
  float  power_a = 5;
  float power_b = 5;
  float radius;
  OOAlloc(RepRibbon);
  RepInit(&I->R);
  power_a=SettingGet(cSetting_ribbon_power);
  power_b=SettingGet(cSetting_ribbon_power_b);

  sampling=SettingGet(cSetting_ribbon_sampling);
  radius=SettingGet(cSetting_ribbon_radius);

  obj = cs->Obj;
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepRibbonRender;
  I->R.fFree=(void (*)(struct Rep *))RepRibbonFree;
  I->R.fRecolor=NULL;

  /* find all of the CA points */

  at = Alloc(int,cs->NIndex);
  pv = Alloc(float,cs->NIndex*3);
  seg = Alloc(int,cs->NIndex);
  
  i=at;
  v=pv;
  s=seg;

  nAt = 0;
  nSeg = 0;
  a2=-1;
  for(a1=0;a1<cs->NAtIndex;a1++)
	 {
		a=cs->AtmToIdx[a1];
		if(a>=0)
		  if(obj->AtomInfo[a1].visRep[cRepRibbon])
			 if(!obj->AtomInfo[a1].hetatm)
				if(WordMatch("CA",obj->AtomInfo[a1].name,1)<0)
				  {
					 if(a2>=0) {
						if((abs(obj->AtomInfo[a1].resv-obj->AtomInfo[a2].resv)>1)||
							(obj->AtomInfo[a1].chain[0]!=obj->AtomInfo[a2].chain[0])||
							(!WordMatch(obj->AtomInfo[a1].segi,obj->AtomInfo[a2].segi,1)))
						  {
							 a2=-1;
						  }
					 }
					 if(a2<=0)
						nSeg++;
					 *(s++) = nSeg;
					 nAt++;
					 *(i++)=a;
					 v1 = cs->Coord+3*a;		
					 *(v++)=*(v1++);
					 *(v++)=*(v1++);
					 *(v++)=*(v1++);
					 
					 a2=a1;
				  }
	 }
  if(nAt)
	 {
		/* compute differences and normals */

		s=seg;
		v=pv;
		
		dv = Alloc(float,nAt*3);
		nv = Alloc(float,nAt*3);
		dl = Alloc(float,nAt);
		v1=dv;
		v2=nv;
		d=dl;
		
		for(a=0;a<(nAt-1);a++)
		  {
			 if(*s==*(s+1))
				{
				  subtract3f(v+3,v,v1);
				  *d = length3f(v1);
				  scale3f(v1,1.0/(*d),v2);
				}
			 d++;
			 v+=3;
			 v1+=3;
			 v2+=3;
			 s++;
		  }
		
		/* compute tangents */
		
		s=seg;
		v=nv;
		
		tv = Alloc(float,nAt*3);
		v1=tv;
		
		*(v1++)=*(v++); /* first segment */
		*(v1++)=*(v++);
		*(v1++)=*(v++);
		s++;
		
		for(a=1;a<(nAt-1);a++)
		  {
			 if((*s==*(s-1))&&(*s==*(s+1)))
				{
				  add3f(v,(v-3),v1);
				  normalize3f(v1);			 
				}
			 else if(*s==*(s-1))
				{
				  *(v1)=*(v-3);  /* end a segment */
				  *(v1+1)=*(v-2); 
				  *(v1+2)=*(v-1); 
				}
			 else if(*s==*(s+1))
				{
				  *(v1)=*(v);   /* new segment */
				  *(v1+1)=*(v+1); 
				  *(v1+2)=*(v+2); 
				}
			 v+=3;
			 v1+=3;
			 s++;
		  }
		
		*(v1++)=*(v-3); /* last segment */
		*(v1++)=*(v-2);
		*(v1++)=*(v-1);
		
	 }

  /* okay, we now have enough info to generate smooth interpolations */
  
  I->VC=(float*)mmalloc(sizeof(float)*cs->NIndex*10*sampling);
  ErrChkPtr(I->VC);
  I->NC=0;

  /*  */

  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*9*sampling);
  ErrChkPtr(I->V);

  I->N=0;
  v=I->V;
  vc=I->VC;
  if(nAt) {
	 v1=pv;
	 v2=tv;
	 v3=dv;
	 s=seg;
	 atp=at;
	 for(a=0;a<(nAt-1);a++)
		{
		  if(*s==*(s+1))
			 {
				c1=*(cs->Color+*atp);
				c2=*(cs->Color+*(atp+1));

				for(b=0;b<sampling;b++) /* needs optimization */
				  {

					 f0=((float)b)/sampling;
					 if(f0<0.5) {
						v0 = ColorGet(c1);
					 } else {
						v0 = ColorGet(c2);
					 }

					 *(v++)=*(v0);
					 *(vc++)=*(v0++);
					 *(v++)=*(v0);
					 *(vc++)=*(v0++);
					 *(v++)=*(v0);
					 *(vc++)=*(v0);

					 f1=2*(f0-0.5);
					 if(f1<0.0) 
						f0=(-pow(-f1,1.0/(power_b))+1.0)/2.0;
					 else
						f0=(pow(f1,1.0/(power_b))+1.0)/2.0;
					 f1=1-f0;
					 f2=1-(f0*2);
					 f3=1-fabs(pow(f2,power_a));
					 f4=f0;
					 *(v++)=v1[0]+f4*v3[0]+f3*((v2[0]*f1)-(v2[3]*f0));
					 *(v++)=v1[1]+f4*v3[1]+f3*((v2[1]*f1)-(v2[4]*f0));
					 *(v++)=v1[2]+f4*v3[2]+f3*((v2[2]*f1)-(v2[5]*f0));
					 *(vc++)=radius;
					 *(vc++)=*(v-3);
					 *(vc++)=*(v-2);
					 *(vc++)=*(v-1);

					 f0=((float)b+1)/sampling;
					 f1=2*(f0-0.5);
					 if(f1<0.0) 
						f0=(-pow(-f1,1.0/(power_b))+1.0)/2.0;
					 else
						f0=(pow(f1,1.0/(power_b))+1.0)/2.0;
					 f1=1-f0;
					 f2=1-(f0*2);
					 f3=1-fabs(pow(f2,power_a));
					 f4=f0;

					 *(v++)=v1[0]+f4*v3[0]+f3*((v2[0]*f1)-(v2[3]*f0));
					 *(v++)=v1[1]+f4*v3[1]+f3*((v2[1]*f1)-(v2[4]*f0));
					 *(v++)=v1[2]+f4*v3[2]+f3*((v2[2]*f1)-(v2[5]*f0));
					 *(vc++)=*(v-3);
					 *(vc++)=*(v-2);
					 *(vc++)=*(v-1);
					 I->NC++;
					 I->N++;
				  }
			 }
		  v1+=3;
		  v2+=3;
		  v3+=3;
		  atp+=1;
		  s++;
		}
	 FreeP(dv);
	 FreeP(dl);
	 FreeP(tv);
	 FreeP(nv);
  }

  FreeP(at);
  FreeP(seg);
  FreeP(pv);
  
  if(I->N) 
	 I->V=(float*)mrealloc(I->V,sizeof(float)*(v-I->V));
  else
	 I->V=(float*)mrealloc(I->V,1);
  
  return((void*)(struct Rep*)I);
}


