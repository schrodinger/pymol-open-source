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
#include"Vector.h"
#include"ObjectMolecule.h"
#include"RepCylBond.h"
#include"Color.h"
#include"Setting.h"
#include"main.h"

typedef struct RepCylBond {
  Rep R;
  float *V,*VR;
  int N,NR;
  int NEdge;

} RepCylBond;

void RepCylBondRender(RepCylBond *I,CRay *ray,Pickable **pick);
void subdivide( int n, float *x, float *y);
float *RepCylinder(float *v,float *v1,float *v2,int nEdge,int endCap);

void RepCylBondFree(RepCylBond *I);

void RepCylBondFree(RepCylBond *I)
{
  FreeP(I->VR);
  FreeP(I->V);
  OOFreeP(I);
}

void RepCylBondRender(RepCylBond *I,CRay *ray,Pickable **pick)
{
  int a;
  float *v=I->V;
  int c=I->N;

  if(ray) {
	 v=I->VR;
	 c=I->NR;
	 while(c--) {
		ray->fCylinder3fv(ray,v+4,v+7,*(v+3),v,v);
		v+=10;
	 }
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 while(c--)
		{
		  
		  glColor3fv(v);
		  v+=3;
		  
		  glBegin(GL_TRIANGLE_STRIP);
		  a=I->NEdge+1;
		  while(a--) {
			 glNormal3fv(v);
			 v+=3;
			 glVertex3fv(v);
			 v+=3;
			 glVertex3fv(v);
			 v+=3;
		  }
		  glEnd();
		  
		  glBegin(GL_TRIANGLE_FAN);
		  glNormal3fv(v);
		  v+=3;
		  glVertex3fv(v);
		  v+=3;
		  a=I->NEdge+1;
		  while(a--) {
			 glNormal3fv(v);
			 v+=3;
			 glVertex3fv(v);
			 v+=3;
		  }
		  glEnd();
		  
		  if(*(v++)) {
			 
			 glBegin(GL_TRIANGLE_FAN);
			 glNormal3fv(v);
			 v+=3;
			 glVertex3fv(v);
			 v+=3;
			 a=I->NEdge+1;
			 while(a--) {
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

Rep *RepCylBondNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,a2,*b,c1,c2,s1,s2,b1,b2;
  float *v,*vv1,*vv2,*v0,*vr;
  float v1[3],v2[3];
  float radius;
  int nEdge;
  int half_bonds;
  int visFlag;
  OOAlloc(RepCylBond);

  obj = cs->Obj;
  visFlag=false;
  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = *(b++);
      b2 = *(b++);
      b++;
      if(obj->AtomInfo[b1].visRep[cRepCyl]||
         obj->AtomInfo[b2].visRep[cRepCyl]) {
        visFlag=true;
        break;
      }
    }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  nEdge = SettingGet(cSetting_stick_quality);
  radius = SettingGet(cSetting_stick_radius);
  half_bonds = SettingGet(cSetting_half_bonds);  

  RepInit(&I->R);
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepCylBondRender;
  I->R.fFree=(void (*)(struct Rep *))RepCylBondFree;

  I->V = NULL;
  I->VR = NULL;
  I->N = 0;
  I->NR = 0;
  I->R.fRecolor=NULL;

  if(obj->NBond) {
	 I->V = (float*)mmalloc(((obj->NBond)*((nEdge+2)*42)+20)*sizeof(float));
	 ErrChkPtr(I->V);

	 I->VR=(float*)mmalloc(sizeof(float)*obj->NBond*10*3);
	 ErrChkPtr(I->VR);
	 
	 I->NEdge = nEdge;
	 
	 v=I->V;
	 vr=I->VR;
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
				c1=*(cs->Color+a1);
				c2=*(cs->Color+a2);
				
				vv1 = cs->Coord+3*a1;
				vv2 = cs->Coord+3*a2;
				
				s1=obj->AtomInfo[b1].visRep[cRepCyl];
				s2=obj->AtomInfo[b2].visRep[cRepCyl];

				if(!(s1&&s2))
              if(!half_bonds) {
                s1 = 0;
                s2 = 0;
              }
				
				if(s1||s2)
				  {
					 
					 if((c1==c2)&&s1&&s2)
						{
						  
						  v1[0]=vv1[0];
						  v1[1]=vv1[1];
						  v1[2]=vv1[2];
						  
						  v2[0]=vv2[0];
						  v2[1]=vv2[1];
						  v2[2]=vv2[2];
						  
						  v0 = ColorGet(c1);

						  *(vr++)=*(v0);
						  *(vr++)=*(v0+1);
						  *(vr++)=*(v0+2);
						  *(vr++)=radius;						  

						  *(vr++)=*(v1);
						  *(vr++)=*(v1+1);
						  *(vr++)=*(v1+2);

						  *(vr++)=*(v2);
						  *(vr++)=*(v2+1);
						  *(vr++)=*(v2+2);




						  I->NR++;

 						  *(v++)=*(v0++);
						  *(v++)=*(v0++);
						  *(v++)=*(v0++);
						  
						  I->N++;
						  

						  v=RepCylinder(v,v1,v2,nEdge,1);
						  
						}
					 else
						{
						  
						  v1[0]=vv1[0];
						  v1[1]=vv1[1];
						  v1[2]=vv1[2];
						  
						  v2[0]=(vv1[0]+vv2[0])/2;
						  v2[1]=(vv1[1]+vv2[1])/2;
						  v2[2]=(vv1[2]+vv2[2])/2;
						  
						  if(s1) 
							 {
								v0 = ColorGet(c1);

								*(vr++)=*(v0);
								*(vr++)=*(v0+1);
								*(vr++)=*(v0+2);
								*(vr++)=radius;
								
								*(vr++)=*(v1);
								*(vr++)=*(v1+1);
								*(vr++)=*(v1+2);
								
								*(vr++)=*(v2);
								*(vr++)=*(v2+1);
								*(vr++)=*(v2+2);
								
								I->NR++;
								
								*(v++)=*(v0++);
								*(v++)=*(v0++);
								*(v++)=*(v0++);

								I->N++;
								v=RepCylinder(v,v1,v2,nEdge,0);
							 }
						  
						  v1[0]=vv2[0];
						  v1[1]=vv2[1];
						  v1[2]=vv2[2];
						  
						  if(s2) 
							 {
								v0 = ColorGet(c2);


								*(vr++)=*(v0);
								*(vr++)=*(v0+1);
								*(vr++)=*(v0+2);
								*(vr++)=radius;
								
								*(vr++)=*(v1);
								*(vr++)=*(v1+1);
								*(vr++)=*(v1+2);
								
								*(vr++)=*(v2);
								*(vr++)=*(v2+1);
								*(vr++)=*(v2+2);
								
								I->NR++;

								*(v++)=*(v0++);
								*(v++)=*(v0++);
								*(v++)=*(v0++);
								
								I->N++;
								v=RepCylinder(v,v1,v2,nEdge,0);
							 }
						}
				  }
			 }
		}
	 /*	 printf(" RepCylBond: %d triplets\n",(v-I->V)/3);*/
	 I->V = Realloc(I->V,float,(v-I->V));
	 I->VR = Realloc(I->VR,float,(vr-I->VR));
  }
  return((void*)(struct Rep*)I);
}




void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++)
	 {
		x[a]=cos(a*2*PI/n);
		y[a]=sin(a*2*PI/n);
	 }
}

float *RepCylinder(float *v,float *v1,float *v2,int nEdge,int endCap)
{

  float d[3],t[3],p0[3],p1[3],p2[3];
  float x[50],y[50];
  float tube_size;
  float overlap;
  float nub;
  int c;

  tube_size = SettingGet(cSetting_stick_radius);
  overlap = tube_size*SettingGet(cSetting_stick_overlap);
  nub = tube_size*SettingGet(cSetting_stick_nub);

 subdivide(nEdge,x,y);

  /* direction vector */
  
  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  normalize3f(p0);
  
  v1[0]-=p0[0]*overlap;
  v1[1]-=p0[1]*overlap;
  v1[2]-=p0[2]*overlap;
  
  if(endCap) {
	 v2[0]+=p0[0]*overlap;
	 v2[1]+=p0[1]*overlap;
	 v2[2]+=p0[2]*overlap;
  }
  
  d[0] = (v2[0] - v1[0]);
  d[1] = (v2[1] - v1[1]);
  d[2] = (v2[2] - v1[2]);
  
  t[0] = d[1];
  t[1] = d[2];
  t[2] = -d[0];
  
  cross_product3f(d,t,p1);
  
  normalize3f(p1);
  
  cross_product3f(d,p1,p2);
  
  normalize3f(p2);
  
  /* now we have a coordinate system*/
  
  for(c=0;c<=nEdge;c++)
	 {
		v[0] = p1[0]*tube_size*x[c] + p2[0]*tube_size*y[c];
		v[1] = p1[1]*tube_size*x[c] + p2[1]*tube_size*y[c];
		v[2] = p1[2]*tube_size*x[c] + p2[2]*tube_size*y[c];
		
		v[3] = v1[0] + v[0];
		v[4] = v1[1] + v[1];
		v[5] = v1[2] + v[2];
		
		v[6] = v[3] + d[0];
		v[7] = v[4] + d[1];
		v[8] = v[5] + d[2];
		
		normalize3f(v);
		v+=9;			 
	 }
  
  v[0] = -p0[0];
  v[1] = -p0[1];
  v[2] = -p0[2];
  
  v[3] = v1[0] - p0[0]*nub;
  v[4] = v1[1] - p0[1]*nub;
  v[5] = v1[2] - p0[2]*nub;
  
  v+=6;
  
  for(c=0;c<=nEdge;c++)
	 {
		
		v[0] = p1[0]*tube_size*x[c] + p2[0]*tube_size*y[c];
		v[1] = p1[1]*tube_size*x[c] + p2[1]*tube_size*y[c];
		v[2] = p1[2]*tube_size*x[c] + p2[2]*tube_size*y[c];
		
		v[3] = v1[0] + v[0];
		v[4] = v1[1] + v[1];
		v[5] = v1[2] + v[2];
		
		v+=6;
	 }

  if(endCap) 
	 {
		*(v++)=1.0;
		v[0] = p0[0];
		v[1] = p0[1];
		v[2] = p0[2];
		
		v[3] = v2[0] + p0[0]*nub;
		v[4] = v2[1] + p0[1]*nub;
		v[5] = v2[2] + p0[2]*nub;
		
		v+=6;
		
		for(c=0;c<=nEdge;c++)
		  {
			 
			 v[0] = p1[0]*tube_size*x[c] + p2[0]*tube_size*y[c];
			 v[1] = p1[1]*tube_size*x[c] + p2[1]*tube_size*y[c];
			 v[2] = p1[2]*tube_size*x[c] + p2[2]*tube_size*y[c];
			 
			 v[3] = v2[0] + v[0];
			 v[4] = v2[1] + v[1];
			 v[5] = v2[2] + v[2];
			 
			 v+=6;
		  }
	 }
  else 
	 {
		*(v++)=0.0;
	 }
  
  return(v);
}



