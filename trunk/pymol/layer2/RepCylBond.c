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
#include"Sphere.h"

typedef struct RepCylBond {
  Rep R;
  float *V,*VR;
  int N,NR;
  int NEdge;
  float *VP;
  int NP;
  float *VSP,*VSPC;
  SphereRec *SP;
  int NSP,NSPC;
} RepCylBond;

void RepCylBondRender(RepCylBond *I,CRay *ray,Pickable **pick);
static void subdivide( int n, float *x, float *y);
float *RepCylinder(float *v,float *v1,float *v2,int nEdge,
                   int frontCap, int endCap,
                   float tube_size,float overlap,float nub);

void RepCylBondFree(RepCylBond *I);

void RepCylBondFree(RepCylBond *I)
{
  FreeP(I->VR);
  FreeP(I->VP);
  FreeP(I->V);
  FreeP(I->VSP);
  FreeP(I->VSPC);
  RepFree(&I->R);
  OOFreeP(I);
}

float *RepCylinderBox(float *v,float *v1,float *v2,float tube_size,
                      float overlap,float nub);


void RepCylBondRender(RepCylBond *I,CRay *ray,Pickable **pick)
{
  int a;
  float *v;
  int c,cc;
  int i,j;
  Pickable *p;
  float alpha;
  SphereRec *sp;

  alpha = SettingGet_f(I->R.cs->Setting,I->R.obj->Setting,cSetting_stick_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  if(ray) {
    ray->fTransparentf(ray,1.0F-alpha);    
    PRINTFD(FB_RepCylBond)
      " RepCylBondRender: rendering raytracable...\n"
      ENDFD;

	 v=I->VR;
	 c=I->NR;
	 while(c--) {
		ray->fSausage3fv(ray,v+4,v+7,*(v+3),v,v);
		v+=10;
	 }
    if(I->VSPC) {
      v=I->VSPC;
      c=I->NSPC;
      while(c--) {
        ray->fColor3fv(ray,v);
        v+=3;
        ray->fSphere3fv(ray,v,*(v+3));
        v+=4;
      }
    }

    ray->fTransparentf(ray,0.0);
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
    
    int use_dlst;
    use_dlst = (int)SettingGet(cSetting_use_display_lists);
    if(use_dlst&&I->R.displayList) {
      glCallList(I->R.displayList);
    } else { 
      
      if(use_dlst) {
        if(!I->R.displayList) {
          I->R.displayList = glGenLists(1);
          if(I->R.displayList) {
            glNewList(I->R.displayList,GL_COMPILE_AND_EXECUTE);
          }
        }
      }
      
      v=I->V;
      c=I->N;
      
      PRINTFD(FB_RepCylBond)
        " RepCylBondRender: rendering GL...\n"
        ENDFD;
      
      while(c--)
        {
          /* cylinder entry consists of a color, a fan,
             a cylinder, and another fan (if flagged) */
          
          if(alpha==1.0) {
            glColor3fv(v);
          } else {
            glColor4f(v[0],v[1],v[2],alpha);
          }
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

      if(I->VSP) { /* stick spheres, if present */
        
        v = I->VSP;
        c = I->NSP;
        if(alpha==1.0) {
          
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
        } else {
          sp=I->SP;
          while(c--)
            {
              glColor4f(v[0],v[1],v[2],alpha);
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
      

      PRINTFD(FB_RepCylBond)
        " RepCylBondRender: done.\n"
        ENDFD;
      
      if(use_dlst&&I->R.displayList) {
        glEndList();
      }
    }
  }
}

static void RepValence(float **v_ptr,int *n_ptr, /* opengl */
                       float **vr_ptr,int *nr_ptr, /* ray */
                       float *v1,float *v2,int *other,
                       int a1,int a2,float *coord,
                       float *color1,float *color2,int ord,
                       int n_edge,
                       float tube_size,
                       float overlap,
                       float nub,
                       int half_bonds,
                       int fixed_r)
{

  float d[3],t[3],p0[3],p1[3],p2[3],*vv;
  float v1t[3],v2t[3],vh[3];
  float *v = *v_ptr,*vr = *vr_ptr;
  int n = *n_ptr,nr = *nr_ptr;
  int a3,ck;

  /* First, we need to construct a coordinate system */

  /* get direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  copy3f(p0,d);
  normalize3f(p0);
  
  /* need a third atom to get planarity*/

  a3 = -1;
  ck = other[a1*2];
  if((ck>=0)&&(ck!=a2))
    a3=ck;
  else {
    ck = other[a1*2+1];
    if((ck>=0)&&(ck!=a2))
      a3=ck;
    else {
      ck = other[a2*2];
      if((ck>=0)&&(ck!=a1))
        a3=ck;
      else {
        ck = other[a2*2+1];
        if((ck>=0)&&(ck!=a1))
          a3=ck;
        }
      }
  }

  if(a3<0) {    
    t[0] = p0[0];
    t[1] = p0[1];
    t[2] = -p0[2];
  } else {
    vv= coord+3*a3;
    t[0] = *(vv++)-v1[0];
    t[1] = *(vv++)-v1[1];
    t[2] = *(vv++)-v1[2];
    normalize3f(t);
  }
  
  cross_product3f(d,t,p1);
  
  normalize3f(p1);

  if(length3f(p1)==0.0) {
    p1[0]=p0[1];
    p1[1]=p0[2];
    p1[2]=p0[0];
    cross_product3f(p0,p1,p2);
    normalize3f(p2);
  } else {
    cross_product3f(d,p1,p2);
    
    normalize3f(p2);
  }

  /* we have a coordinate system*/

  /* Next, we need to determine how many cylinders */
  
  switch(ord) {
  case 2:
    {
      float radius = tube_size;
      float overlap_r;
      float nub_r;
      if(!fixed_r) {
        radius/=2.5;
      }
      overlap_r = radius*overlap;
      nub_r = radius*nub;

      t[0] = p2[0]*1.5*radius;
      t[1] = p2[1]*1.5*radius;
      t[2] = p2[2]*1.5*radius;
      
      if(!half_bonds) {

        /* opengl */

        copy3f(color1,v);
        v+=3;
        
        add3f(v1,t,v1t);
        add3f(v2,t,v2t);
        
        v=RepCylinder(v,v1t,v2t,n_edge,1,1,radius,
                      overlap_r,nub_r);
        n++;

        /* ray */

        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius;						          
        copy3f(v1t,vr);
        vr+=3;        
        copy3f(v2t,vr);
        vr+=3;
        nr++;

        /* opengl */

        copy3f(color1,v);
        
        v+=3;
        
        subtract3f(v1,t,v1t);
        subtract3f(v2,t,v2t);
        
        v=RepCylinder(v,v1t,v2t,n_edge,1,1,radius,
                      overlap_r,nub_r);

        /* ray */

        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius;						          
        copy3f(v1t,vr);
        vr+=3;        
        copy3f(v2t,vr);
        vr+=3;
        nr++;

        n++;
      } else {
        vh[0] = (v1[0]+v2[0])*0.5F;
        vh[1] = (v1[1]+v2[1])*0.5F;
        vh[2] = (v1[2]+v2[2])*0.5F;

        if(color1) {

        /* opengl */
          copy3f(color1,v);
          v+=3;

          add3f(v1,t,v1t);
          add3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;

        /* opengl */
          copy3f(color1,v);
          v+=3;

          subtract3f(v1,t,v1t);
          subtract3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;
          
        }
        if(color2) {

          copy3f(color2,v);
          v+=3;

          add3f(v2,t,v1t);
          add3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color2,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;

          copy3f(color2,v);
          v+=3;

          subtract3f(v2,t,v1t);
          subtract3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;
          /* ray */
          
          copy3f(color2,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;


        }
      }
    }
    break;
  case 3:
    {
      float radius = tube_size;
      float overlap_r;
      float nub_r;
      if(!fixed_r) {
        radius/=3.5;
      }
      overlap_r = radius*overlap;
      nub_r = radius*nub;

      t[0] = p2[0]*2.5*radius;
      t[1] = p2[1]*2.5*radius;
      t[2] = p2[2]*2.5*radius;

      if(!half_bonds) {
        /* opengl */
        copy3f(color1,v);
        
        v+=3;
        
        v=RepCylinder(v,v1,v2,n_edge,1,1,radius,
                      overlap_r,nub_r);
        n++;

        /* ray */
        
        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius;						          
        copy3f(v1,vr);
        vr+=3;        
        copy3f(v2,vr);
        vr+=3;
        nr++;

        /* opengl */
        copy3f(color1,v);
        
        v+=3;
        
        add3f(v1,t,v1t);
        add3f(v2,t,v2t);
        
        v=RepCylinder(v,v1t,v2t,n_edge,1,1,radius,
                      overlap_r,nub_r);
        n++;
        
        /* ray */
        
        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius;						          
        copy3f(v1t,vr);
        vr+=3;        
        copy3f(v2t,vr);
        vr+=3;
        nr++;

        /* opengl */
        copy3f(color1,v);
        
        v+=3;
        
        subtract3f(v1,t,v1t);
        subtract3f(v2,t,v2t);
        
        v=RepCylinder(v,v1t,v2t,n_edge,1,1,radius,
                      overlap_r,nub_r);
        n++;

        /* ray */
        
        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius;						          
        copy3f(v1t,vr);
        vr+=3;        
        copy3f(v2t,vr);
        vr+=3;
        nr++;

      } else {
        vh[0] = (v1[0]+v2[0])*0.5F;
        vh[1] = (v1[1]+v2[1])*0.5F;
        vh[2] = (v1[2]+v2[2])*0.5F;

        if(color1) {
          /* opengl */
          copy3f(color1,v);
          v+=3;
          
          v=RepCylinder(v,v1,vh,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;
          
          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1,vr);
          vr+=3;        
          copy3f(vh,vr);
          vr+=3;
          nr++;

          /* opengl */
          copy3f(color1,v);
          v+=3;

          add3f(v1,t,v1t);
          add3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;

          /* opengl */
          copy3f(color1,v);
          v+=3;

          subtract3f(v1,t,v1t);
          subtract3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;

          
        }
        if(color2) {
          
          /* opengl */
          copy3f(color2,v);
          v+=3;
          
          v=RepCylinder(v,v2,vh,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color2,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v2,vr);
          vr+=3;        
          copy3f(vh,vr);
          vr+=3;
          nr++;

          /* opengl */
          copy3f(color2,v);
          v+=3;

          add3f(v2,t,v1t);
          add3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color2,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;

          /* opengl */
          copy3f(color2,v);
          v+=3;

          subtract3f(v2,t,v1t);
          subtract3f(vh,t,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
                        overlap_r,nub_r);
          n++;

          /* ray */
          
          copy3f(color2,vr);
          vr+=3;        
          *(vr++)=radius;						          
          copy3f(v1t,vr);
          vr+=3;        
          copy3f(v2t,vr);
          vr+=3;
          nr++;

        }
      }      
    }

    break;
  case 4: 
    
    break;
  }
  *v_ptr = v;
  *n_ptr = n;
  *vr_ptr = vr;
  *nr_ptr = nr;

}


Rep *RepCylBondNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,a2,c1,c2,s1,s2,b1,b2,*o;
  BondType *b;
  float *v,*vv1,*vv2,*v0,*vr,*vsp,*vspc;
  float v1[3],v2[3],h[3];
  float radius;
  int nEdge;
  float valence;
  float overlap,nub,overlap_r,nub_r;
  int half_bonds,*other=NULL;
  int visFlag;
  int maxCyl;
  int ord;
  int stick_ball;
  float stick_ball_ratio=1.0F;
  unsigned int v_size,vr_size,rp_size,vp_size;
  Pickable *rp;
  AtomInfoType *ai1,*ai2;
  SphereRec *sp = NULL;
  float *rgb1,*rgb2,rgb1_buf[3],rgb2_buf[3];
  int fixed_radius = false;
  int caps_req = true;

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

  valence = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_valence);

  maxCyl = 0;

  b=obj->Bond;
  for(a=0;a<obj->NBond;a++)
    {
      b1 = b->index[0];
      b2 = b->index[1];
      ord = b->order;
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
          obj->AtomInfo[a1].temp1 = false; /* use the kludge field for sphere marker */
          obj->AtomInfo[a2].temp1 = false;
          
          if(valence!=0.0F) {
            switch(ord) {
            case 1:
              maxCyl+=2;
              break;
            case 2:
              maxCyl+=4;
              break;
            case 3:
              maxCyl+=6;
              break;
            case 4:
              maxCyl+=8;
              break;
            }
          } else
            maxCyl+=2;
        }
    }


  nEdge = (int)SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_quality);
  radius = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_radius);
  half_bonds = (int)SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_half_bonds);  


  RepInit(&I->R);
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepCylBondRender;
  I->R.fFree=(void (*)(struct Rep *))RepCylBondFree;
  I->R.obj=(CObject*)obj;
  I->R.cs = cs;

  I->V = NULL;
  I->VR = NULL;
  I->N = 0;
  I->NR = 0;
  I->NP = 0;
  I->VP = NULL;
  I->SP = NULL;
  I->VSP = NULL;
  I->NSP = 0;
  I->VSPC = NULL;
  I->NSPC = 0;
  if(obj->NBond) {

    stick_ball = SettingGet_b(cs->Setting,obj->Obj.Setting,cSetting_stick_ball);
    overlap = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_overlap);
    nub = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_nub);

    overlap_r = overlap * radius;
    nub_r = nub * radius;

    if(valence!=0.0) /* build list of up to 2 connected atoms for each atom */
      {
        other=Alloc(int,2*obj->NAtom);
        o=other;
        for(a=0;a<obj->NAtom;a++) {
          *(o++)=-1;
          *(o++)=-1;
        }
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
                o=other+2*a1;
                if(*o!=a2) {
                  if(*o>=0) o++;
                  *o=a2;
                }
                o=other+2*a2;              
                if(*o!=a1) {
                  if(*o>=0) o++;
                  *o=a1;
                }
              }
          }
        
        fixed_radius = SettingGet_b(cs->Setting,obj->Obj.Setting,cSetting_stick_fixed_radius);
      }
    
    /* OpenGL */

    v_size = ((maxCyl)*((nEdge+1)*21)+22);
	 I->V = Alloc(float,v_size);
	 ErrChkPtr(I->V);

    /* RayTrace */

    vr_size = maxCyl*10*3;
    I->VR=Alloc(float,vr_size);
	 ErrChkPtr(I->VR);

    /* spheres for stick & balls */
	 
    if(stick_ball) {
      int ds;
      stick_ball_ratio = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_stick_ball_ratio);

      ds = SettingGet_i(cs->Setting,obj->Obj.Setting,cSetting_sphere_quality);
      if(ds<0) ds=0;
      switch(ds) {
      case 0: sp=Sphere0; break;
      case 1: sp=Sphere1; break;
      case 2: sp=Sphere2; break;
      case 3: sp=Sphere3; break;
      default: sp=Sphere4; break;
      }
      I->SP = sp;
      I->VSP=Alloc(float,maxCyl*2*(3+sp->NVertTot*6));
      I->VSPC=Alloc(float,maxCyl*2*7);
      ErrChkPtr(I->VSP);
    }
	 I->NEdge = nEdge;
	 
	 v=I->V;
	 vr=I->VR;
    vsp = I->VSP;
    vspc = I->VSPC;
	 b=obj->Bond;
	 for(a=0;a<obj->NBond;a++)
		{
        b1 = b->index[0];
        b2 = b->index[1];
        ord = b->order;
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
				
            if(stick_ball) {
              float vdw = stick_ball_ratio * radius;
              int d,e;
              if(stick_ball_ratio>=1.0F) /* don't use caps if spheres are big enough */
                caps_req = false;
              if(s1&&(!obj->AtomInfo[b1].temp1)) { /* just once for each atom... */
                int *q=sp->Sequence;
                int *s=sp->StripLen;
                obj->AtomInfo[b1].temp1=1;
                {
                  if(ColorCheckRamped(c1)) {
                    ColorGetRamped(c1,vv1,rgb2_buf);
                    rgb1 = rgb1_buf;
                  } else {
                    rgb1 = ColorGet(c1);
                  }
                }
                copy3f(rgb1,vsp);
                vsp+=3;
                for(d=0;d<sp->NStrip;d++)
                  {
                    for(e=0;e<(*s);e++)
                      {
                        *(vsp++)=sp->dot[*q][0]; /* normal */
                        *(vsp++)=sp->dot[*q][1];
                        *(vsp++)=sp->dot[*q][2];
                        *(vsp++)=vv1[0]+vdw*sp->dot[*q][0]; /* point */
                        *(vsp++)=vv1[1]+vdw*sp->dot[*q][1];
                        *(vsp++)=vv1[2]+vdw*sp->dot[*q][2];
                        q++;
                      }
                    s++;
                  }
                I->NSP++;
                copy3f(rgb1,vspc);
                vspc+=3;
                copy3f(vv1,vspc);
                vspc+=3;
                *(vspc++)=vdw;
                I->NSPC++;
              }
              if(s2&&!(obj->AtomInfo[b2].temp1)) { /* just once for each atom... */
                int *q=sp->Sequence;
                int *s=sp->StripLen;
                obj->AtomInfo[b2].temp1=1;
                
                if(ColorCheckRamped(c2)) {
                  ColorGetRamped(c2,vv2,rgb2_buf);
                  rgb2 = rgb2_buf;
                } else {
                  rgb2 = ColorGet(c2);
                }
              
                copy3f(rgb2,vsp);
                vsp+=3;

                for(d=0;d<sp->NStrip;d++)
                  {
                    for(e=0;e<(*s);e++)
                      {
                        *(vsp++)=sp->dot[*q][0]; /* normal */
                        *(vsp++)=sp->dot[*q][1];
                        *(vsp++)=sp->dot[*q][2];
                        *(vsp++)=vv2[0]+vdw*sp->dot[*q][0]; /* point */
                        *(vsp++)=vv2[1]+vdw*sp->dot[*q][1];
                        *(vsp++)=vv2[2]+vdw*sp->dot[*q][2];
                        q++;
                      }
                    s++;
                  }
                I->NSP++;

                copy3f(rgb2,vspc);
                vspc+=3;
                copy3f(vv2,vspc);
                vspc+=3;
                *(vspc++)=vdw;
                I->NSPC++;
              }
            }

				if(s1||s2)
				  {
					 
                if((valence!=0.0)&&(ord>1)&&(ord<4)) {
                  
                  if((c1==c2)&&s1&&s2&&(!ColorCheckRamped(c1))) {

                    
                    v0 = ColorGet(c1);

                    RepValence(&v,&I->N,
                               &vr,&I->NR,
                               vv1,vv2,other,
                               a1,a2,cs->Coord,
                               v0,NULL,ord,nEdge,
                               radius,
                               overlap,
                               nub,
                               false,
                               fixed_radius);
                  } else {

                    rgb1 = NULL;
                    if(s1) {
                      if(ColorCheckRamped(c1)) {
                        ColorGetRamped(c1,vv1,rgb1_buf);
                        rgb1 = rgb1_buf;
                      } else {
                        rgb1 = ColorGet(c1);
                      }
                    }

                    rgb2 = NULL;
                    if(s2) 
                      {
                        if(ColorCheckRamped(c2)) {
                          ColorGetRamped(c2,vv2,rgb2_buf);
                          rgb2 = rgb2_buf;
                        } else {
                          rgb2 = ColorGet(c2);
                        }
                      }
                    
                     RepValence(&v,&I->N,
                               &vr,&I->NR,
                               vv1,vv2,other,
                               a1,a2,cs->Coord,
                               rgb1,rgb2,ord,nEdge,
                               radius,
                               overlap,
                               nub,
                               true,
                                fixed_radius);
                  }
                  
                } else {
                  
                  if((c1==c2)&&s1&&s2&&(!ColorCheckRamped(c1)))
                    {
                      
                      copy3f(vv1,v1);
                      copy3f(vv2,v2);
                      
                      v0 = ColorGet(c1);
                      
                      /* ray-tracing */
                    
                      copy3f(v0,vr);
                      vr+=3;
                    
                      *(vr++)=radius;						  
                    
                      copy3f(v1,vr);
                      vr+=3;
                    
                      copy3f(v2,vr);
                      vr+=3;
                    
                      I->NR++;
                    
                      /* store color */
                    
                      copy3f(v0,v);
                      v+=3;
                    
                      I->N++;
                    
                      /* generate a cylinder */
                    
                      v=RepCylinder(v,v1,v2,nEdge,caps_req,caps_req,radius,overlap_r,nub_r);
                    } else {                    
                      v1[0]=vv1[0];
                      v1[1]=vv1[1];
                      v1[2]=vv1[2];
                      
                      v2[0]=(vv1[0]+vv2[0])*0.5F;
                      v2[1]=(vv1[1]+vv2[1])*0.5F;
                      v2[2]=(vv1[2]+vv2[2])*0.5F;
                      
                      if(s1) 
                        {
                          if(ColorCheckRamped(c1)) {
                            ColorGetRamped(c1,v1,vr);
                            v0=vr;
                            vr+=3;
                          } else {
                            v0 = ColorGet(c1);
                            *(vr++)=*(v0);
                            *(vr++)=*(v0+1);
                            *(vr++)=*(v0+2);
                          }
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
                          v=RepCylinder(v,v1,v2,nEdge,caps_req,0,radius,overlap_r,nub_r);
                        }
                      
                      v1[0]=vv2[0];
                      v1[1]=vv2[1];
                      v1[2]=vv2[2];
                      
                      if(s2) 
                        {
                          if(ColorCheckRamped(c2)) {
                            ColorGetRamped(c2,v1,vr);
                            v0=vr;
                            vr+=3;
                          } else {
                            v0 = ColorGet(c2);
                            
                            *(vr++)=*(v0);
                            *(vr++)=*(v0+1);
                            *(vr++)=*(v0+2);
                          }
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
                          v=RepCylinder(v,v1,v2,nEdge,caps_req,0,radius,overlap_r,nub_r);
                        }
                    }
                }
				  }
			 }
		}
	 PRINTFD(FB_RepCylBond)
      " RepCylBond-DEBUG: %d triplets\n",(v-I->V)/3
      ENDFD;

    if((signed)v_size<(v-I->V))
      ErrFatal("RepCylBond","V array overrun.");
    if((signed)vr_size<(vr-I->VR))
      ErrFatal("RepCylBond","VR array overrun.");
    
	 I->V = Realloc(I->V,float,(v-I->V));
	 I->VR = Realloc(I->VR,float,(vr-I->VR));
    if(I->VSP) 
      I->VSP = Realloc(I->VSP,float,(vsp-I->VSP));
    if(I->VSPC) 
      I->VSPC = Realloc(I->VSPC,float,(vspc-I->VSPC));
	 if(SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_pickable)) { 

      PRINTFD(FB_RepCylBond)
        " RepCylBondNEW: generating pickable version\n"
        ENDFD;
      /* pickable versions are simply capped boxes, 
         vertices: 8 points * 3 = 32  * 2 = 48 floats per bond
      */

      vp_size = maxCyl*24;
      I->VP=Alloc(float,vp_size);
		ErrChkPtr(I->VP);
		
      rp_size = maxCyl+1;
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

                      v = RepCylinderBox(v,v1,h,radius,overlap_r,nub_r);
						  }
						if(s2&(!ai2->masked))
						  {
							 I->NP++;
                      rp->ptr = (void*)obj;
							 rp->index = b2;
                      rp->bond = a;
                      rp++;

                      v = RepCylinderBox(v,h,v2,radius,overlap_r,nub_r);
						  }
					 }
				}
		  }
      
      if((signed)vp_size<(v-I->VP))
        ErrFatal("RepCylBond","VP array overrun.");
      if((signed)rp_size<=(I->NP))
        ErrFatal("RepCylBond","RP array overrun.");

		I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
		I->R.P[0].index = I->NP;
      I->VP = Realloc(I->VP,float,(v-I->VP));

      PRINTFD(FB_RepCylBond)
        " RepCylBondNew: I->NP: %d I->VP: %p\n",I->NP,I->VP
        ENDFD;
	 }
  }
  FreeP(other);

  return((void*)(struct Rep*)I);
}

static void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++)
	 {
		x[a]=(float)cos(a*2*PI/n);
		y[a]=(float)sin(a*2*PI/n);
	 }
}


float *RepCylinderBox(float *v,float *vv1,float *vv2,
                      float tube_size,float overlap,float nub)
{
  
  float d[3],t[3],p0[3],p1[3],p2[3],n[3];
  float v1[3],v2[3];

  tube_size *= 0.7F;

  overlap+=(nub/2);

  /* direction vector */

  subtract3f(vv2,vv1,p0);

  normalize3f(p0);

  v1[0]=vv1[0]-p0[0]*overlap;
  v1[1]=vv1[1]-p0[1]*overlap;
  v1[2]=vv1[2]-p0[2]*overlap;
  
  v2[0]=vv2[0]+p0[0]*overlap;
  v2[1]=vv2[1]+p0[1]*overlap;
  v2[2]=vv2[2]+p0[2]*overlap;
    
  d[0] = (v2[0] - v1[0]);
  d[1] = (v2[1] - v1[1]);
  d[2] = (v2[2] - v1[2]);
  
  t[0] = d[0];
  t[1] = d[1];
  t[2] = d[2];

  get_system1f3f(t,p1,p2);

  /*
  cross_product3f(d,t,p1);
  
  normalize3f(p1);
  
  cross_product3f(d,p1,p2);
  
  normalize3f(p2);*/
  
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


float *RepCylinder(float *v,float *v1,float *v2,int nEdge,
                   int frontCap, int endCap,
                   float tube_size,float overlap,float nub)
{

  float d[3],t[3],p0[3],p1[3],p2[3];
  float x[50],y[50];
  int c;


  if(nEdge>50)
    nEdge=50;
 subdivide(nEdge,x,y);
 
  /* direction vector */
  
  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
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
  
  
  t[0] = d[0];
  t[1] = d[1];
  t[2] = d[2];

  get_system1f3f(t,p1,p2);

  /* 
  cross_product3f(d,t,p1);
  
  normalize3f(p1);
  
  cross_product3f(d,p1,p2);
  
  normalize3f(p2);*/

  
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
  
  if(frontCap) {

    *(v++)=1.0;
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
  }
  else 
	 {
		*(v++)=0.0F;
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
		*(v++)=0.0F;
	 }
  
  return(v);
}



