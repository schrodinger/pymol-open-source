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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Base.h"
#include"OOMac.h"
#include"Vector.h"
#include"ObjectMolecule.h"
#include"RepCylBond.h"
#include"Color.h"
#include"Setting.h"
#include"main.h"
#include"Feedback.h"

typedef struct RepCylBond {
  Rep R;
  float *V,*VR;
  int N,NR;
  int NEdge;
  float *VP;
  int NP;
} RepCylBond;

void RepCylBondRender(RepCylBond *I,CRay *ray,Pickable **pick);
static void subdivide( int n, float *x, float *y);
float *RepCylinder(float *v,float *v1,float *v2,int nEdge,int endCap,CoordSet *cs,
                   ObjectMolecule *obj);

void RepCylBondFree(RepCylBond *I);

void RepCylBondFree(RepCylBond *I)
{
  FreeP(I->VR);
  FreeP(I->VP);
  FreeP(I->V);
  RepFree(&I->R);
  OOFreeP(I);
}

float *RepCylinderBox(float *v,float *v1,float *v2,CoordSet *cs,ObjectMolecule *obj);

void RepCylBondRender(RepCylBond *I,CRay *ray,Pickable **pick)
{
  int a;
  float *v;
  int c;
  int i,j;
  Pickable *p;


  if(ray) {
    
    PRINTFD(FB_RepCylBond)
      " RepCylBondRender: rendering raytracable...\n"
      ENDFD;

	 v=I->VR;
	 c=I->NR;
	 while(c--) {
		ray->fSausage3fv(ray,v+4,v+7,*(v+3),v,v);
		v+=10;
	 }
  } else if(pick&&PMGUI) {

  PRINTFD(FB_RepCylBond)
    " RepCylBondRender: rendering pickable...\n"
    ENDFD;

	 i=(*pick)->index;

	 v=I->VP;
	 c=I->NP;
	 p=I->R.P;

	 while(c--) {

		i++;

		if(!(*pick)[0].ptr) {
		  /* pass 1 - low order bits */

		  glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
		  VLACheck((*pick),Pickable,i);
		  p++;
		  (*pick)[i] = *p; /* copy object and atom info */
		} else { 
		  /* pass 2 - high order bits */

		  j=i>>12;

		  glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4)); 

		}			 
      
      glBegin(GL_TRIANGLE_STRIP);

		glVertex3fv(v+ 0);
		glVertex3fv(v+ 3);

		glVertex3fv(v+ 6);
		glVertex3fv(v+ 9);

		glVertex3fv(v+12);
		glVertex3fv(v+15);

		glVertex3fv(v+18);
		glVertex3fv(v+21);

		glVertex3fv(v+ 0);
		glVertex3fv(v+ 3);

      glEnd();

      glBegin(GL_TRIANGLE_STRIP);
      
		glVertex3fv(v+ 0);
		glVertex3fv(v+ 6);

		glVertex3fv(v+18);
		glVertex3fv(v+12);

      glEnd();

      glBegin(GL_TRIANGLE_STRIP);
      
		glVertex3fv(v+ 3);
		glVertex3fv(v+ 9);

		glVertex3fv(v+21);
		glVertex3fv(v+15);

      glEnd();

      v+=24;

	 }
	 (*pick)[0].index = i; /* pass the count */

  } else if(PMGUI) {

    {
      
      v=I->V;
      c=I->N;

      PRINTFD(FB_RepCylBond)
        " RepCylBondRender: rendering GL...\n"
        ENDFD;
      
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
      PRINTFD(FB_RepCylBond)
        " RepCylBondRender: done.\n"
        ENDFD;
    }
  }
}

Rep *RepCylBondNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,a2,c1,c2,s1,s2,b1,b2;
  BondType *b;
  float *v,*vv1,*vv2,*v0,*vr;
  float v1[3],v2[3],h[3];
  float radius;
  int nEdge;
  int half_bonds;
  int visFlag;
  unsigned int v_size,vr_size,rp_size,vp_size;
  Pickable *rp;
  AtomInfoType *ai1,*ai2;

  OOAlloc(RepCylBond);

  PRINTFD(FB_RepCylBond)
    " RepCylBondNew-Debug: entered.\n"
    ENDFD;

  obj = cs->Obj;
  visFlag=false;
  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = b->index[0];
      b2 = b->index[1];
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

  nEdge = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_quality);
  radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_radius);
  half_bonds = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_half_bonds);  

  RepInit(&I->R);
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepCylBondRender;
  I->R.fFree=(void (*)(struct Rep *))RepCylBondFree;

  I->V = NULL;
  I->VR = NULL;
  I->N = 0;
  I->NR = 0;
  I->NP = 0;
  I->VP = NULL;

  if(obj->NBond) {

    v_size = ((obj->NBond)*((nEdge+2)*42)+32);
	 I->V = Alloc(float,v_size);
	 ErrChkPtr(I->V);

    vr_size = obj->NBond*10*3;
    I->VR=Alloc(float,vr_size);
	 ErrChkPtr(I->VR);
	 
	 I->NEdge = nEdge;
	 
	 v=I->V;
	 vr=I->VR;
	 b=obj->Bond;
	 for(a=0;a<obj->NBond;a++)
		{
        b1 = b->index[0];
        b2 = b->index[1];
        b++;

        if(obj->DiscreteFlag) {
          if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
            a1=obj->DiscreteAtmToIdx[b1];
            a2=obj->DiscreteAtmToIdx[b2];
          } else {
            a1=-1;
            a2=-1;
          }
        } else {
          a1=cs->AtmToIdx[b1];
          a2=cs->AtmToIdx[b2];
        }
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
						  

						  v=RepCylinder(v,v1,v2,nEdge,1,cs,obj);
						  
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
								v=RepCylinder(v,v1,v2,nEdge,0,cs,obj);
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
								v=RepCylinder(v,v1,v2,nEdge,0,cs,obj);
							 }
						}
				  }
			 }
		}
	 PRINTFD(FB_RepCylBond)
      " RepCylBond-DEBUG: %d triplets\n",(v-I->V)/3
      ENDFD;

    if(v_size<(v-I->V))
      ErrFatal("RepCylBond","V array overrun.");
    if(vr_size<(vr-I->VR))
      ErrFatal("RepCylBond","VR array overrun.");
    
	 I->V = Realloc(I->V,float,(v-I->V));
	 I->VR = Realloc(I->VR,float,(vr-I->VR));

	 if(SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_pickable)) { 

      PRINTFD(FB_RepCylBond)
        " RepCylBondNEW: generating pickable version\n"
        ENDFD;
      /* pickable versions are simply capped boxes, 
         vertices: 8 points * 3 = 32  * 2 = 48 floats per bond
      */

      vp_size = obj->NBond*48;
      I->VP=Alloc(float,vp_size);
		ErrChkPtr(I->VP);
		
      rp_size = 2*obj->NBond+1;
		I->R.P=Alloc(Pickable,rp_size);
		ErrChkPtr(I->R.P);
		rp = I->R.P + 1; /* skip first record! */

		v=I->VP;
		b=obj->Bond;
		for(a=0;a<obj->NBond;a++)
		  {
          b1 = b->index[0];
          b2 = b->index[1];
          b++;
          if(obj->DiscreteFlag) {
            if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
              a1=obj->DiscreteAtmToIdx[b1];
              a2=obj->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=cs->AtmToIdx[b1];
            a2=cs->AtmToIdx[b2];
          }
			 if((a1>=0)&&(a2>=0))
				{
              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;
				  s1=ai1->visRep[cRepCyl];
				  s2=ai2->visRep[cRepCyl];
				  
              if(!(s1&&s2)) {
                if(!half_bonds) {
                  s1 = 0;
                  s2 = 0;
                }
              } 

				  if(s1||s2)
					 {	
						copy3f(cs->Coord+3*a1,v1);
                    copy3f(cs->Coord+3*a2,v2);
						
						h[0]=(v1[0]+v2[0])/2;
						h[1]=(v1[1]+v2[1])/2;
						h[2]=(v1[2]+v2[2])/2;
						
						if(s1&(!ai1->masked))
						  {
							 I->NP++;
                      rp->ptr = (void*)obj;
							 rp->index = b1;
                      rp->bond = a;
                      rp++;

                      v = RepCylinderBox(v,v1,h,cs,obj);
						  }
						if(s2&(!ai2->masked))
						  {
							 I->NP++;
                      rp->ptr = (void*)obj;
							 rp->index = b2;
                      rp->bond = a;
                      rp++;

                      v = RepCylinderBox(v,h,v2,cs,obj);
						  }
					 }
				}
		  }
      
      if(vp_size<(v-I->VP))
        ErrFatal("RepCylBond","VP array overrun.");
      if(rp_size<=(I->NP))
        ErrFatal("RepCylBond","RP array overrun.");

		I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
		I->R.P[0].index = I->NP;
      I->VP = Realloc(I->VP,float,(v-I->VP));

      PRINTFD(FB_RepCylBond)
        " RepCylBondNew: I->NP: %d I->VP: %p\n",I->NP,I->VP
        ENDFD;
	 }
  }

  return((void*)(struct Rep*)I);
}

static void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++)
	 {
		x[a]=cos(a*2*PI/n);
		y[a]=sin(a*2*PI/n);
	 }
}


float *RepCylinderBox(float *v,float *v1,float *v2,CoordSet *cs,ObjectMolecule *obj)
{

  float d[3],t[3],p0[3],p1[3],p2[3],n[3];
  float tube_size;
  float overlap;
  float nub;

  tube_size = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_radius)
    *0.7;

  overlap = tube_size*SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_overlap);
  nub = tube_size*SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_nub);

  overlap+=(nub/2);

  /* direction vector */
  
  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  normalize3f(p0);
  
  v1[0]-=p0[0]*overlap;
  v1[1]-=p0[1]*overlap;
  v1[2]-=p0[2]*overlap;
  
  v2[0]+=p0[0]*overlap;
  v2[1]+=p0[1]*overlap;
  v2[2]+=p0[2]*overlap;
  
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

  n[0] = p1[0]*tube_size*(-1) + p2[0]*tube_size*(-1);
  n[1] = p1[1]*tube_size*(-1) + p2[1]*tube_size*(-1);
  n[2] = p1[2]*tube_size*(-1) + p2[2]*tube_size*(-1);

  v[0] = v1[0] + n[0];
  v[1] = v1[1] + n[1];
  v[2] = v1[2] + n[2];
		
  v[3] = v[0] + d[0];
  v[4] = v[1] + d[1];
  v[5] = v[2] + d[2];

  v+=6;

  n[0] = p1[0]*tube_size*( 1) + p2[0]*tube_size*(-1);
  n[1] = p1[1]*tube_size*( 1) + p2[1]*tube_size*(-1);
  n[2] = p1[2]*tube_size*( 1) + p2[2]*tube_size*(-1);

  v[0] = v1[0] + n[0];
  v[1] = v1[1] + n[1];
  v[2] = v1[2] + n[2];
		
  v[3] = v[0] + d[0];
  v[4] = v[1] + d[1];
  v[5] = v[2] + d[2];

  v+=6;

  n[0] = p1[0]*tube_size*( 1) + p2[0]*tube_size*( 1);
  n[1] = p1[1]*tube_size*( 1) + p2[1]*tube_size*( 1);
  n[2] = p1[2]*tube_size*( 1) + p2[2]*tube_size*( 1);

  v[0] = v1[0] + n[0];
  v[1] = v1[1] + n[1];
  v[2] = v1[2] + n[2];
		
  v[3] = v[0] + d[0];
  v[4] = v[1] + d[1];
  v[5] = v[2] + d[2];

  v+=6;

  n[0] = p1[0]*tube_size*(-1) + p2[0]*tube_size*( 1);
  n[1] = p1[1]*tube_size*(-1) + p2[1]*tube_size*( 1);
  n[2] = p1[2]*tube_size*(-1) + p2[2]*tube_size*( 1);

  v[0] = v1[0] + n[0];
  v[1] = v1[1] + n[1];
  v[2] = v1[2] + n[2];
		
  v[3] = v[0] + d[0];
  v[4] = v[1] + d[1];
  v[5] = v[2] + d[2];

  v+=6;

  return(v);
}


float *RepCylinder(float *v,float *v1,float *v2,int nEdge,int endCap,
                   CoordSet *cs,ObjectMolecule *obj)
{

  float d[3],t[3],p0[3],p1[3],p2[3];
  float x[50],y[50];
  float tube_size;
  float overlap;
  float nub;
  int c;

  tube_size = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_radius);
  overlap = tube_size*SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_overlap);
  nub = tube_size*SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_nub);

  if(nEdge>50)
    nEdge=50;
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
  
  for(c=nEdge;c>=0;c--)
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
  
  for(c=nEdge;c>=0;c--)
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



