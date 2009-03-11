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
  float *VarAlpha,*VarAlphaRay,*VarAlphaSph;
} RepCylBond;

static void subdivide( int n, float *x, float *y);
float *RepCylinder(float *v,float *v1,float *v2,int nEdge,
                   int frontCap, int endCap,
                   float tube_size,float overlap,float nub);

void RepCylBondFree(RepCylBond *I);

void RepCylBondFree(RepCylBond *I)
{
  FreeP(I->VarAlpha);
  FreeP(I->VarAlphaRay);
  FreeP(I->VarAlphaSph);
  FreeP(I->VR);
  FreeP(I->VP);
  FreeP(I->V);
  FreeP(I->VSP);
  FreeP(I->VSPC);
  RepPurge(&I->R);
  OOFreeP(I);
}

float *RepCylinderBox(float *v,float *v1,float *v2,float tube_size,
                      float overlap,float nub);


static void RepCylBondRender(RepCylBond *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int a;
  float *v,*var_alpha;
  int c,cc;
  int i,j;
  Pickable *p;
  float alpha;
  SphereRec *sp;
  register PyMOLGlobals *G=I->R.G;


  alpha = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_stick_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  if(ray) {
    ray->fTransparentf(ray,1.0F-alpha);    
    PRINTFD(G,FB_RepCylBond)
      " RepCylBondRender: rendering raytracable...\n"
      ENDFD;

    v=I->VR;
    c=I->NR;
    var_alpha = I->VarAlphaRay;
    while(c--) {
      if(var_alpha) {
        ray->fTransparentf(ray,1.0F-*(var_alpha++));
      }
      ray->fSausage3fv(ray,v+4,v+7,*(v+3),v,v);
      v+=10;
    }
    var_alpha = I->VarAlphaSph;
    if(I->VSPC) {
      v=I->VSPC;
      c=I->NSPC;
      while(c--) {
        if(var_alpha) {
          ray->fTransparentf(ray,1.0F-*(var_alpha++));
        }
        ray->fColor3fv(ray,v);
        v+=3;
        ray->fSphere3fv(ray,v,*(v+3));
        v+=4;
      }
    }

    ray->fTransparentf(ray,0.0);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      
      PRINTFD(G,FB_RepCylBond)
        " RepCylBondRender: rendering pickable...\n"
        ENDFD;

      i=(*pick)->src.index;

      v=I->VP;
      c=I->NP;
      p=I->R.P;

      while(c--) {

        i++;

        if(!(*pick)[0].src.bond) {
          /* pass 1 - low order bits */

          glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
          VLACheck((*pick),Picking,i);
          p++;
          (*pick)[i].src = *p; /* copy object and atom info */
          (*pick)[i].context = I->R.context;
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
      (*pick)[0].src.index = i; /* pass the count */

    } else {
      int use_dlst;    

      use_dlst = (int)SettingGet(G,cSetting_use_display_lists);
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
        var_alpha = I->VarAlpha;
        PRINTFD(G,FB_RepCylBond)
          " RepCylBondRender: rendering GL...\n"
          ENDFD;
      
        while(c--) {
          /* cylinder entry consists of a color, a fan,
             a cylinder, and another fan (if flagged) */
          
          if((alpha==1.0)&&(!var_alpha)) {
            glColor3fv(v);
          } else if(var_alpha) { 
            glColor4f(v[0],v[1],v[2],*(var_alpha++));
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
          var_alpha = I->VarAlphaSph;

          if((alpha==1.0)&&!(var_alpha)) {
            
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
            while(c--) {
              if(!var_alpha) {
                glColor4f(v[0],v[1],v[2],alpha);
              } else {
                glColor4f(v[0],v[1],v[2],*(var_alpha++));
              }
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
        PRINTFD(G,FB_RepCylBond)
          " RepCylBondRender: done.\n"
          ENDFD;
      
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
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
                       int fixed_r,
                       float scale_r)
{

  float d[3],t[3],p0[3],p1[3],p2[3],*vv;
  float v1t[3],v2t[3],vh[3];
  float *v = *v_ptr,*vr = *vr_ptr;
  int n = *n_ptr,nr = *nr_ptr;
  int a3;
  int double_sided;
      
  /* First, we need to construct a coordinate system */

  /* get direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  copy3f(p0,d);
  normalize3f(p0);
  
  /* need a third atom to get planarity*/
  a3 = ObjectMoleculeGetPrioritizedOther(other,a1,a2,&double_sided);

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
        radius*=scale_r;
        radius/=2.5;
      }

      overlap_r = radius*overlap;
      nub_r = radius*nub;

      t[0] = p2[0]*1.5F*radius;
      t[1] = p2[1]*1.5F*radius;
      t[2] = p2[2]*1.5F*radius;
      
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
        radius*=scale_r;
        radius/=3.5;
      }

      overlap_r = radius*overlap;
      nub_r = radius*nub;

      t[0] = p2[0]*2.5F*radius;
      t[1] = p2[1]*2.5F*radius;
      t[2] = p2[2]*2.5F*radius;

      if(!half_bonds) {
        /* opengl */
        copy3f(color1,v);
        
        v+=3;
        copy3f(v1,v1t);
        copy3f(v2,v2t);
        v=RepCylinder(v,v1t,v2t,n_edge,1,1,radius,
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
          
          copy3f(v1,v1t);
          v=RepCylinder(v,v1t,vh,n_edge,1,0,radius,
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
          
          copy3f(v2,v2t);
          v=RepCylinder(v,v2t,vh,n_edge,1,0,radius,
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
    {
      float radius = tube_size;
      float radius2 = tube_size;
      float overlap_r,overlap_r2;
      float nub_r,nub_r2;
      float along[3],adj[3],v1tt[3],v2tt[3];
      float inner1a = 0.24F;
      float inner1b = 0.44F;
      float inner2a = 0.5F+(0.5F-inner1b);
      float inner2b = 1.0F-inner1a;

      if(!fixed_r) {
        radius*=scale_r;
        radius2=radius/2.5F;
        t[0] = p2[0]*1.5F*radius;
        t[1] = p2[1]*1.5F*radius;
        t[2] = p2[2]*1.5F*radius;
      } else {
        inner1a -= 0.04F;
        inner2b = 1.0F-inner1a;
        t[0] = p2[0]*3*radius;
        t[1] = p2[1]*3*radius;
        t[2] = p2[2]*3*radius;
      }
      
      overlap_r = radius*overlap;
      nub_r = radius*nub;
      overlap_r2 = radius2*overlap;
      nub_r2 = radius2*nub;

      
      if(!half_bonds) {

        /* opengl */

        copy3f(color1,v);
        v+=3;
        
        copy3f(v1,v1t);
        copy3f(v2,v2t);
        v=RepCylinder(v,v1t,v2t,n_edge,1,1,radius,
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

        subtract3f(v1,t,v1t);
        subtract3f(v2,t,v2t);
        
        subtract3f(v2t,v1t,along);
        scale3f(along,inner1a,adj);
        add3f(adj,v1t,v1tt);
        scale3f(along,inner1b,adj);
        add3f(adj,v1t,v2tt);

        /* ray */

        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius2;						          
        copy3f(v1tt,vr);
        vr+=3;        
        copy3f(v2tt,vr);
        vr+=3;
        nr++;

        /* opengl */

        copy3f(color1,v); 
        v+=3;
        
        v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                      overlap_r2,nub_r2);
        n++;


        scale3f(along,inner2a,adj);
        add3f(adj,v1t,v1tt);
        scale3f(along,inner2b,adj);
        add3f(adj,v1t,v2tt);

        /* ray */

        copy3f(color1,vr);
        vr+=3;        
        *(vr++)=radius2;						          
        copy3f(v1tt,vr);
        vr+=3;        
        copy3f(v2tt,vr);
        vr+=3;
        nr++;

        /* opengl */

        copy3f(color1,v);
        v+=3;

        v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                      overlap_r2,nub_r2);        
        n++;

        if(double_sided) {
          
          
          add3f(v1,t,v1t);
          add3f(v2,t,v2t);
          
          subtract3f(v2t,v1t,along);
          scale3f(along,inner1a,adj);
          add3f(adj,v1t,v1tt);
          scale3f(along,inner1b,adj);
          add3f(adj,v1t,v2tt);
          
          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius2;						          
          copy3f(v1tt,vr);
          vr+=3;        
          copy3f(v2tt,vr);
          vr+=3;
          nr++;
          
          /* opengl */
          
          copy3f(color1,v);
          v+=3;

          v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                        overlap_r2,nub_r2);
          n++;
          
          scale3f(along,inner2a,adj);
          add3f(adj,v1t,v1tt);
          scale3f(along,inner2b,adj);
          add3f(adj,v1t,v2tt);
          
          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius2;						          
          copy3f(v1tt,vr);
          vr+=3;        
          copy3f(v2tt,vr);
          vr+=3;
          nr++;
          
          /* opengl */
          
          copy3f(color1,v);          
          v+=3;
          
          v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                        overlap_r2,nub_r2);        
          n++;
        
        }

      } else {
        vh[0] = (v1[0]+v2[0])*0.5F;
        vh[1] = (v1[1]+v2[1])*0.5F;
        vh[2] = (v1[2]+v2[2])*0.5F;

        if(color1) {

          /* opengl */

          copy3f(color1,v);
          v+=3;

          copy3f(v1,v1t);
          copy3f(vh,v2t);
          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
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
          
          subtract3f(v1,t,v1t);
          subtract3f(v2,t,v2t);
          
          subtract3f(v2t,v1t,along);
          scale3f(along,inner1a,adj);
          add3f(adj,v1t,v1tt);
          scale3f(along,inner1b,adj);
          add3f(adj,v1t,v2tt);

          /* ray */
          
          copy3f(color1,vr);
          vr+=3;        
          *(vr++)=radius2;						          
          copy3f(v1tt,vr);
          vr+=3;        
          copy3f(v2tt,vr);
          vr+=3;
          nr++;
          
          /* open gl */

          copy3f(color1,v);
          v+=3;

          v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                        overlap_r2,nub_r2);
          n++;
         
          if(double_sided) {
            
            add3f(v1,t,v1t);
            add3f(v2,t,v2t);
            
            subtract3f(v2t,v1t,along);
            scale3f(along,inner1a,adj);
            add3f(adj,v1t,v1tt);
            scale3f(along,inner1b,adj);
            add3f(adj,v1t,v2tt);
            
            /* ray */
            
            copy3f(color1,vr);
            vr+=3;        
            *(vr++)=radius2;						          
            copy3f(v1tt,vr);
            vr+=3;        
            copy3f(v2tt,vr);
            vr+=3;
            nr++;
            
            /* open gl */
            
            copy3f(color1,v);
            v+=3;

            v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                          overlap_r2,nub_r2);
            n++;
            
            
          }
        }

        if(color2) {

          copy3f(color2,v);
          v+=3;

          copy3f(v2,v1t);
          copy3f(vh,v2t);

          v=RepCylinder(v,v1t,v2t,n_edge,1,0,radius,
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

          subtract3f(v1,t,v1t);
          subtract3f(v2,t,v2t);
          
          subtract3f(v2t,v1t,along);
          scale3f(along,inner2a,adj);
          add3f(adj,v1t,v1tt);
          scale3f(along,inner2b,adj);
          add3f(adj,v1t,v2tt);

          /* ray */
          
          copy3f(color2,vr);
          vr+=3;        
          *(vr++)=radius2;						          
          copy3f(v1tt,vr);
          vr+=3;        
          copy3f(v2tt,vr);
          vr+=3;
          nr++;

          /* opengl */

          copy3f(color2,v);
          v+=3;

          v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                        overlap_r2,nub_r2);
          n++;
          
          if(double_sided) {
            
            add3f(v1,t,v1t);
            add3f(v2,t,v2t);
            
            subtract3f(v2t,v1t,along);
            scale3f(along,inner2a,adj);
            add3f(adj,v1t,v1tt);
            scale3f(along,inner2b,adj);
            add3f(adj,v1t,v2tt);
            
            /* ray */
            
            copy3f(color2,vr);
            vr+=3;        
            *(vr++)=radius2;						          
            copy3f(v1tt,vr);
            vr+=3;        
            copy3f(v2tt,vr);
            vr+=3;
            nr++;
            
            /* opengl */
            
            copy3f(color2,v);
            v+=3;
            
            v=RepCylinder(v,v1tt,v2tt,n_edge,1,1,radius2,
                          overlap_r2,nub_r2);
            n++;
            
          }
        }
      }
    }
    break;
  }
  *v_ptr = v;
  *n_ptr = n;
  *vr_ptr = vr;
  *nr_ptr = nr;

}


Rep *RepCylBondNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  int a,a1,a2,c1,c2,s1,s2;
  register int b1,b2;
  register BondType *b;
  float *v,*vv1,*vv2,*v0,*vr,*vsp,*vspc;
  float v1[3],v2[3],h[3];
  float radius;
  int nEdge;
  float valence;
  float overlap,nub;
  int half_bonds,*other=NULL;
  int visFlag;
  int maxCyl;
  int ord;
  int stick_ball, stick_ball_color=-1;
  float stick_ball_ratio=1.0F;
  unsigned int v_size,vr_size,rp_size,vp_size;
  Pickable *rp;
  register AtomInfoType *ai1,*ai2;
  SphereRec *sp = NULL;
  float *rgb1,*rgb2,rgb1_buf[3],rgb2_buf[3];
  int fixed_radius = false;
  int caps_req = true;
  int valence_flag = false;
  int hide_long = false;
  int stick_color = 0;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 1;
  int na_mode;
  int *marked = NULL;
  float scale_r = 1.0F;
  int variable_alpha = false;
  int n_var_alpha=0, n_var_alpha_ray=0,n_var_alpha_sph=0;
  float transp,h_scale;
  int valence_found  = false;
  const float _0p9 = 0.9F;
  
  OOAlloc(G,RepCylBond);

  PRINTFD(G,FB_RepCylBond)
    " RepCylBondNew-Debug: entered.\n"
    ENDFD;

  obj = cs->Obj;
  visFlag=false;
  b=obj->Bond;
  ai1=obj->AtomInfo;
  if(obj->RepVisCache[cRepCyl])
    for(a=0;a<obj->NBond;a++)
      {
        b1 = b->index[0];
        b2 = b->index[1];
        
        b++;
        if(ai1[b1].visRep[cRepCyl]||
           ai1[b2].visRep[cRepCyl]) {
          visFlag=true;
          break;
        }
      }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  marked = Calloc(int,obj->NAtom);

  valence = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_valence);
  valence_flag = (valence!=0.0F);
  maxCyl = 0;

  stick_color = SettingGet_color(G,cs->Setting,obj->Obj.Setting,cSetting_stick_color);
  cartoon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                           cSetting_cartoon_side_chain_helper);
  ribbon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                          cSetting_ribbon_side_chain_helper);

  transp = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_transparency);
  hide_long = SettingGet_b(G,cs->Setting,obj->Obj.Setting,cSetting_hide_long_bonds);

  b=obj->Bond;
  for(a=0;a<obj->NBond;a++) {
      b1 = b->index[0];
      b2 = b->index[1];
      ord = b->order;
      
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
      if((a1>=0)&&(a2>=0)) {
        int bd_valence_flag;

        if((!variable_alpha) && AtomInfoCheckBondSetting(G,b,cSetting_stick_transparency))
          variable_alpha = true;
        
        AtomInfoGetBondSetting_b(G,b,cSetting_valence,valence_flag,&bd_valence_flag);
        
        if(bd_valence_flag) {
          valence_found = true;
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
            maxCyl+=6;
            break;
          }
        } else
          maxCyl+=2;
      }
      b++;
  }


  nEdge = (int)SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_quality);
  radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_radius);
  half_bonds = (int)SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_half_bonds);  
  na_mode = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_cartoon_nucleic_acid_mode);
  h_scale = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_h_scale);
  RepInit(G,&I->R);
  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepCylBondRender;
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
  I->VarAlpha = NULL;
  I->VarAlphaRay = NULL;
  I->VarAlphaSph = NULL;

  if(obj->NBond) {
    int draw_mode = SettingGetGlobal_i(G,cSetting_draw_mode);
    int draw_quality = (((draw_mode == 1)||(draw_mode==-2)));
      
    stick_ball = SettingGet_b(G,cs->Setting,obj->Obj.Setting,cSetting_stick_ball);
    overlap = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_overlap);
    nub = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_nub);
    
    if(draw_quality) {
      nEdge = 16;
      overlap = 0.05;
    }

    if(cartoon_side_chain_helper || ribbon_side_chain_helper) {
      /* mark atoms that are bonded to atoms without a
         visible cartoon or ribbon */

      b=obj->Bond;
      for(a=0;a<obj->NBond;a++) {

          b1 = b->index[0];
          b2 = b->index[1];
          ord = b->order;
          
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
          if((a1>=0)&&(a2>=0)) {
            register AtomInfoType *ati1=obj->AtomInfo+b1;
            register AtomInfoType *ati2=obj->AtomInfo+b2;


            if((!ati1->hetatm) && (!ati2->hetatm)) {
              if(((cartoon_side_chain_helper && ati1->visRep[cRepCartoon] && !ati2->visRep[cRepCartoon]) ||
                  (ribbon_side_chain_helper && ati1->visRep[cRepRibbon] && !ati2->visRep[cRepRibbon]))) {
                marked[b1] = 1;
              }
              if(((cartoon_side_chain_helper && ati2->visRep[cRepCartoon] && !ati1->visRep[cRepCartoon]) ||
                  (ribbon_side_chain_helper && ati2->visRep[cRepRibbon] && !ati1->visRep[cRepRibbon]))) {
                marked[b2] = 1;
              }
            }
          }
          b++;
      }

    }


    if(valence_found) {/* build list of up to 2 connected atoms for each atom */
      other=ObjectMoleculeGetPrioritizedOtherIndexList(obj,cs);
      fixed_radius = SettingGet_b(G,cs->Setting,obj->Obj.Setting,cSetting_stick_fixed_radius);
      scale_r = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_valence_scale);
    }
    
    /* OpenGL */

    v_size = maxCyl * ( (9+6+6) * (nEdge+1) + 3 );
    /* each cylinder is 9+6+6+3= 21*(nEdge+1) floats for coordinates
       plus 3 more floats for color */

    /*    printf("debug maxCyl: %d nEdge: %d v_size: %d\n",maxCyl,nEdge,v_size);*/

    if(variable_alpha) 
      I->VarAlpha = Alloc(float,maxCyl);
    I->V = Alloc(float,v_size);
    ErrChkPtr(G,I->V);

    /* RayTrace */

    vr_size = maxCyl*10*3;
    if(variable_alpha) 
      I->VarAlphaRay = Alloc(float,maxCyl);
    I->VR=Alloc(float,vr_size);
    ErrChkPtr(G,I->VR);

    /* spheres for stick & balls */
	if(stick_ball) {
      stick_ball_ratio = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_ball_ratio);
      stick_ball_color = SettingGet_b(G,cs->Setting,obj->Obj.Setting,cSetting_stick_ball_color);
    } else if(draw_quality) {
      stick_ball = true;
    }
    if(stick_ball) {
      int ds;

      ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sphere_quality);
      if(ds<0) ds=0;
      if(ds>4) ds=4;

      if(draw_quality && (ds<3)) ds = 3;

      sp = G->Sphere->Sphere[ds];

      I->SP = sp;
      I->VSP=Alloc(float,maxCyl*2*(3+sp->NVertTot*6));
      I->VSPC=Alloc(float,maxCyl*2*7);
      I->VarAlphaSph = Alloc(float,maxCyl*2);
      ErrChkPtr(G,I->VSP);
    }
    I->NEdge = nEdge;
	 
    v=I->V;
    vr=I->VR;
    vsp = I->VSP;
    vspc = I->VSPC;
    b=obj->Bond;
    for(a=0;a<obj->NBond;a++) {
      
      b1 = b->index[0];
      b2 = b->index[1];
      ord = b->order;
      
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
      if((a1>=0)&&(a2>=0)) {

        register AtomInfoType *ati1=obj->AtomInfo+b1;
        register AtomInfoType *ati2=obj->AtomInfo+b2;
        int bd_stick_color;        
        float bd_radius,bd_radius_full;
        float overlap_r,nub_r;
        float bd_transp;
	
        AtomInfoGetBondSetting_color(G,b,cSetting_stick_color,stick_color,&bd_stick_color);
        AtomInfoGetBondSetting_f(G,b,cSetting_stick_radius,radius,&bd_radius);
        if(variable_alpha) 
          AtomInfoGetBondSetting_f(G,b,cSetting_stick_transparency,transp,&bd_transp);

	bd_radius_full=fabs(bd_radius);
	if(bd_radius<0.0F) {
	  bd_radius = -bd_radius;
	  if((ati1->protons == cAN_H)||(ati2->protons == cAN_H))
	    bd_radius = bd_radius*h_scale; /* scaling for bonds involving hydrogen */

	}
        overlap_r = overlap * bd_radius;
        nub_r = nub * bd_radius;
        
        if(bd_stick_color<0) {
          if(bd_stick_color==cColorObject) {
            c1 = (c2 = obj->Obj.Color);
          } else if(ColorCheckRamped(G,bd_stick_color)) {
            c1 = (c2 = bd_stick_color);            
          } else {
            c1=*(cs->Color+a1);
            c2=*(cs->Color+a2);
          }
        } else {
          c1 = (c2 = bd_stick_color);
        }
        
        vv1 = cs->Coord+3*a1;
        vv2 = cs->Coord+3*a2;
        
        s1=ati1->visRep[cRepCyl];
        s2=ati2->visRep[cRepCyl];
        
        if(!(s1&&s2))
          if(!half_bonds) {
            s1 = 0;
            s2 = 0;
          }

        if(hide_long && (s1||s2)) {
          float cutoff = (ati1->vdw + ati2->vdw) * _0p9;
          ai1 = obj->AtomInfo + b1;
          ai2 = obj->AtomInfo + b2;
          if(!within3f(vv1,vv2,cutoff)) /* atoms separated by more than 90% of the sum of their vdw radii */
            s1 = s2 = 0;
        }
        
        if( (!ati1->hetatm) && (!ati2->hetatm) &&
            ((cartoon_side_chain_helper && ati1->visRep[cRepCartoon] && ati2->visRep[cRepCartoon]) ||
             (ribbon_side_chain_helper && ati1->visRep[cRepRibbon] && ati2->visRep[cRepRibbon]))) {
          
          register char *name1=ati1->name;
          register int prot1=ati1->protons;
          register char *name2=ati2->name;
          register int prot2=ati2->protons;
          
          if(prot1 == cAN_C) { 
            if((name1[1]=='A')&&(name1[0]=='C')&&(!name1[2])) { /* CA */
              if(prot2 == cAN_C) { 
                if((name2[1]=='B')&&(name2[0]=='C')&&(!name2[2]))
                  c1 = c2;  /* CA-CB */
                else if((!name2[1])&&(name2[0]=='C')&&(!marked[b2]))
                  s1 = s2 = 0; /* suppress CA-C */
              } else if(prot2 == cAN_H) 
                s1 = s2 = 0; /* suppress all CA-hydrogens */
            } else if((na_mode==1)&&(prot2 == cAN_C)) {
              if((((name2[3]==0)&&
                   ((name2[2]=='*')||(name2[2]=='\''))&&
                   (name2[1]=='5')&&
                   (name2[0]=='C')))&&
                 (((name1[3]==0)&&
                   ((name1[2]=='*')||(name1[2]=='\''))&&
                   (name1[1]=='4')&&
                   (name1[0]=='C'))))
                s1 = s2 = 0;
            }
          } else if(prot1 == cAN_N) { 
            if((!name1[1])&&(name1[0]=='N')) { /* N */
              if(prot2 == cAN_C) {
                if((name2[1]=='D')&&(name2[0]=='C')&&(!name2[2])) 
                  c1 = c2; /* N->CD in PRO */
                else if((name2[1]=='A')&&(name2[0]=='C')&&(!name2[2])&&(!marked[b1]))
                  {
                    char *resn2 = ati2->resn;
                    if(!((resn2[0]=='P')&&(resn2[1]=='R')&&(resn2[2]=='O')))
                      s1 = s2 = 0; /* suppress N-CA, except in pro */
                    else
                      c1 = c2;
                  }
                else if((!name2[1])&&(name2[0]=='C')&&(!marked[b1]))
                  s1 = s2 = 0; /* suppress N-C */
              } else if(prot2 == cAN_H)
                s1 = s2 = 0; /* suppress all N-hydrogens */
            }
          } else if((prot1 == cAN_O)&&(prot2 == cAN_C)) { 
            if((!name2[1])&&(name2[0]=='C')&&
               (((!name1[1])&&(name1[0]=='O'))||
                ((name1[3]==0)&&(name1[2]=='T')&&(name1[1]=='X')&&(name1[0]=='O')))
               &&(!marked[b2]))
              s1 = s2 = 0; /* suppress C-O,OXT */
            else if(na_mode==1) {
              if((((name2[3]==0)&&
                   ((name2[2]=='*')||(name2[2]=='\''))&&
                   ((name2[1]=='3')||(name2[1]=='5'))&&
                   (name2[0]=='C')))&&
                 (((name1[3]==0)&&
                   ((name1[2]=='*')||(name1[2]=='\''))&&
                   ((name1[1]=='3')||(name1[1]=='5'))&&
                   (name1[0]=='O'))))
                s1 = s2 = 0; 
            } 
          } else if((prot1 == cAN_P)&&(prot2 == cAN_O)) {
            if((!name1[1])&&(name1[0]=='P')&&
               (((name2[3]==0)&&(name2[2]=='P')&&
                 ((name2[1]=='1')||(name2[1]=='2')||(name2[1]=='3'))
                 &&(name2[0]=='O'))))
              s1 = s2 = 0; /* suppress P-O1P,O2P,O3P */
            else if(na_mode==1) {
              if((!name1[1])&&(name1[0]=='P')&&
                 (((name2[3]==0)&&
                   ((name2[2]=='*')||(name2[2]=='\''))&&
                   ((name2[1]=='3')||(name2[1]=='5'))&&
                   (name2[0]=='O'))))
                s1 = s2 = 0;
            }
          } 
          if(prot2 == cAN_C) {
            if((name2[1]=='A')&&(name2[0]=='C')&&(!name2[2])) { /* CA */
              if(prot1 == cAN_C) { 
                if((name1[1]=='B')&&(name1[0]=='C')&&(!name1[2]))
                  c2 = c1; /* CA-CB */
                else if((!name1[1])&&(name1[0]=='C')&&(!marked[b1]))
                  s1 = s2 = 0; /* suppress CA-C */
              } else if(prot1 == cAN_H) 
                s1 = s2 = 0; /* suppress all CA-hydrogens */
            } else if((na_mode==1)&&(prot2 == cAN_C)) {
              if((((name1[3]==0)&&
                   ((name1[2]=='*')||(name1[2]=='\''))&&
                   (name1[1]=='5')&&
                   (name1[0]=='C')))&&
                 (((name2[3]==0)&&
                   ((name2[2]=='*')||(name2[2]=='\''))&&
                   (name2[1]=='4')&&
                   (name2[0]=='C'))))
                s1 = s2 = 0;
            }
          } else if(prot2 == cAN_N) {
            if((!name2[1])&&(name2[0]=='N')) { /* N */
              if(prot1 == cAN_C) { 
                if((name1[1]=='D')&&(name1[0]=='C')&&(!name1[2]))
                  c2 = c1; /* N->CD in PRO */
                else if((name1[1]=='A')&&(name1[0]=='C')&&(marked[b2])) 
                  { 
                    char *resn1 = ati1->resn;
                    if(!((resn1[0]=='P')&&(resn1[1]=='R')&&(resn1[2]=='O')))
                      s1 = s2 = 0; /* suppress N-CA, except in pro */
                    else
                      c2 = c1;
                  }
                else if((!name1[1])&&(name1[0]=='C')&&(!marked[b2]))
                  s1 = s2 = 0; /* suppress N-C */
              } else if(prot1 == cAN_H)
                s1 = s2 = 0; /* suppress all N-hydrogens */
            }
          } else if((prot2 == cAN_O)&&(prot1 == cAN_C)) {
            if((!name1[1])&&(name1[0]=='C')&&
               (((!name2[1])&&(name2[0]=='O'))||
                ((name2[3]==0)&&(name2[2]=='T')&&(name2[1]=='X')&&(name2[0]=='O')))
               &&(!marked[b1]))
              s1 = s2 = 0; /* suppress C-O,OXT */
            else if (na_mode==1) {
              if((((name1[3]==0)&&
                   ((name1[2]=='*')||(name1[2]=='\''))&&
                   ((name1[1]=='3')||(name1[1]=='5'))&&
                   (name1[0]=='C')))&&
                 (((name2[3]==0)&&
                   ((name2[2]=='*')||(name2[2]=='\''))&&
                   ((name2[1]=='3')||(name2[1]=='5'))&&
                   (name2[0]=='O'))))
                s1 = s2 = 0;
            }
          } else if((prot2 == cAN_P)&&(prot1 == cAN_O)) {
            if((!name2[1])&&(name2[0]=='P')&&
               (((name1[3]==0)&&(name1[2]=='P')&&
                 ((name1[1]=='1')||(name1[1]=='2')||(name1[1]=='3'))
                 &&(name1[0]=='O'))))
              s1 = s2 = 0; /* suppress P-O1P,O2P,O3P */
            else if(na_mode==1) {
              if((!name2[1])&&(name2[0]=='P')&&
                 (((name1[3]==0)&&
                   ((name1[2]=='*')||(name1[2]=='\''))&&
                   ((name1[1]=='3')||(name1[1]=='5'))&&
                   (name1[0]=='O'))))
                s1 = s2 = 0;
            }
          }
        }
          
        if(stick_ball) {
          int d,e;
          if(stick_ball_ratio>=1.0F) /* don't use caps if spheres are big enough */
            caps_req = false;
          if(s1&&(!marked[b1])) { /* just once for each atom... */
            int *q=sp->Sequence;
            int *s=sp->StripLen;
            float vdw = stick_ball_ratio * ((ati1->protons == cAN_H) ? bd_radius : bd_radius_full);
            float vdw1 = (vdw>=0) ? vdw : -ati1->vdw * vdw;
            int sbc1 = (stick_ball_color==cColorDefault) ? c1 : stick_ball_color;
            
            if(sbc1==cColorAtomic)
              sbc1 = ati1->color;
            
            marked[b1]=1;
            {
              if(ColorCheckRamped(G,sbc1)) {
                ColorGetRamped(G,sbc1,vv1,rgb2_buf,state);
                rgb1 = rgb1_buf;
              } else {
                rgb1 = ColorGet(G,sbc1);
              }
            }
            copy3f(rgb1,vsp);
            vsp+=3;
            for(d=0;d<sp->NStrip;d++) {
              for(e=0;e<(*s);e++) {
                *(vsp++)=sp->dot[*q][0]; /* normal */
                *(vsp++)=sp->dot[*q][1];
                *(vsp++)=sp->dot[*q][2];
                *(vsp++)=vv1[0]+vdw1*sp->dot[*q][0]; /* point */
                *(vsp++)=vv1[1]+vdw1*sp->dot[*q][1];
                *(vsp++)=vv1[2]+vdw1*sp->dot[*q][2];
                q++;
              }
              s++;
            }
            I->NSP++;
            copy3f(rgb1,vspc);
            vspc+=3;
            copy3f(vv1,vspc);
            vspc+=3;
            *(vspc++)=vdw1;
            I->NSPC++;
          }
          
          if(s2&&!(marked[b2])) { /* just once for each atom... */
            int *q=sp->Sequence;
            int *s=sp->StripLen;
            float vdw = stick_ball_ratio * ((ati2->protons == cAN_H) ? bd_radius : bd_radius_full);
            float vdw2 = (vdw>=0) ? vdw : -ati2->vdw * vdw;
            int sbc2 = (stick_ball_color==cColorDefault) ? c2 : stick_ball_color;
            
            marked[b2]=1;
            
            if(sbc2==cColorAtomic)
              sbc2 = ati2->color;
            
            if(ColorCheckRamped(G,sbc2)) {
              ColorGetRamped(G,sbc2,vv2,rgb2_buf,state);
              rgb2 = rgb2_buf;
            } else {
              rgb2 = ColorGet(G,sbc2);
            }
            
            copy3f(rgb2,vsp);
            vsp+=3;
            
            for(d=0;d<sp->NStrip;d++) {
              for(e=0;e<(*s);e++)  {
                *(vsp++)=sp->dot[*q][0]; /* normal */
                *(vsp++)=sp->dot[*q][1];
                *(vsp++)=sp->dot[*q][2];
                *(vsp++)=vv2[0]+vdw2*sp->dot[*q][0]; /* point */
                *(vsp++)=vv2[1]+vdw2*sp->dot[*q][1];
                *(vsp++)=vv2[2]+vdw2*sp->dot[*q][2];
                q++;
              }
              s++;
            }
            I->NSP++;
            
            copy3f(rgb2,vspc);
            vspc+=3;
            copy3f(vv2,vspc);
            vspc+=3;
            *(vspc++)=vdw2;
            I->NSPC++;
          }
        }

        if(hide_long && (s1||s2)) {
          float cutoff = (ati1->vdw + ati2->vdw) * _0p9;
          ai1 = obj->AtomInfo + b1;
          ai2 = obj->AtomInfo + b2;
          if(!within3f(vv1,vv2,cutoff)) /* atoms separated by more than 90% of the sum of their vdw radii */
            s1 = s2 = 0;
        }

        if(s1||s2)
          {
            int bd_valence_flag;

            AtomInfoGetBondSetting_b(G,b,cSetting_valence,valence_flag,&bd_valence_flag);
            
            if((bd_valence_flag)&&(ord>1)&&(ord<5)) {
                  
              if((c1==c2)&&s1&&s2&&(!ColorCheckRamped(G,c1))) {

                    
                v0 = ColorGet(G,c1);

                RepValence(&v,&I->N,
                           &vr,&I->NR,
                           vv1,vv2,other,
                           a1,a2,cs->Coord,
                           v0,NULL,ord,nEdge,
                           bd_radius,
                           overlap,
                           nub,
                           false,
                           fixed_radius,
                           scale_r);
              } else {

                rgb1 = NULL;
                if(s1) {
                  if(ColorCheckRamped(G,c1)) {
                    ColorGetRamped(G,c1,vv1,rgb1_buf,state);
                    rgb1 = rgb1_buf;
                  } else {
                    rgb1 = ColorGet(G,c1);
                  }
                }

                rgb2 = NULL;
                if(s2) 
                  {
                    if(ColorCheckRamped(G,c2)) {
                      ColorGetRamped(G,c2,vv2,rgb2_buf,state);
                      rgb2 = rgb2_buf;
                    } else {
                      rgb2 = ColorGet(G,c2);
                    }
                  }
                    
                RepValence(&v,&I->N,
                           &vr,&I->NR,
                           vv1,vv2,other,
                           a1,a2,cs->Coord,
                           rgb1,rgb2,ord,nEdge,
                           bd_radius,
                           overlap,
                           nub,
                           true,
                           fixed_radius,
                           scale_r);
              }
                  
            } else {
                  
              if((c1==c2)&&s1&&s2&&(!ColorCheckRamped(G,c1)))
                {
                      
                  copy3f(vv1,v1);
                  copy3f(vv2,v2);
                      
                  v0 = ColorGet(G,c1);
                      
                  /* ray-tracing */
                    
                  copy3f(v0,vr);
                  vr+=3;
                    
                  *(vr++)=bd_radius;						  
                    
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
                  v=RepCylinder(v,v1,v2,nEdge,caps_req,caps_req,bd_radius,overlap_r,nub_r);
                } else {                    
                  v1[0]=vv1[0];
                  v1[1]=vv1[1];
                  v1[2]=vv1[2];
                      
                  v2[0]=(vv1[0]+vv2[0])*0.5F;
                  v2[1]=(vv1[1]+vv2[1])*0.5F;
                  v2[2]=(vv1[2]+vv2[2])*0.5F;
                      
                  if(s1) {
                    if(ColorCheckRamped(G,c1)) {
                      ColorGetRamped(G,c1,v1,vr,state);
                      v0=vr;
                      vr+=3;
                    } else {
                      v0 = ColorGet(G,c1);
                      *(vr++)=*(v0);
                      *(vr++)=*(v0+1);
                      *(vr++)=*(v0+2);
                    }
                    *(vr++)=bd_radius;
                    
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
                    v=RepCylinder(v,v1,v2,nEdge,caps_req,0,bd_radius,overlap_r,nub_r);
                  }
                      
                  v1[0]=vv2[0];
                  v1[1]=vv2[1];
                  v1[2]=vv2[2];
                      
                  if(s2) {
                    if(ColorCheckRamped(G,c2)) {
                      ColorGetRamped(G,c2,v1,vr,state);
                      v0=vr;
                      vr+=3;
                    } else {
                      v0 = ColorGet(G,c2);
                      
                      *(vr++)=*(v0);
                      *(vr++)=*(v0+1);
                      *(vr++)=*(v0+2);
                    }
                    *(vr++)=bd_radius;
                    
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
                    v=RepCylinder(v,v1,v2,nEdge,caps_req,0,bd_radius,overlap_r,nub_r);
                  }
                }
            }

            if(variable_alpha) { /* record alpha values for each */
              float bd_alpha = (1.0F - bd_transp);
              while(n_var_alpha<I->N) {
                I->VarAlpha[n_var_alpha++] = bd_alpha;
              }
              while(n_var_alpha_ray<I->NR) {
                I->VarAlphaRay[n_var_alpha_ray++] = bd_alpha;
              }
              while(n_var_alpha_sph<I->NSP) {
                I->VarAlphaSph[n_var_alpha_sph++] = bd_alpha;
              }
            }
          }
      }
      /*      printf("%d\n",(v-I->V)/( (9+6+6) * (nEdge+1) + 3 ));*/
      b++;
    }

    PRINTFD(G,FB_RepCylBond)
      " RepCylBond-DEBUG: %d triplets\n",(int)(v-I->V)/3
      ENDFD;

    if((signed)v_size<(v-I->V))
      ErrFatal(G,"RepCylBond","V array overrun.");
    if((signed)vr_size<(vr-I->VR))
      ErrFatal(G,"RepCylBond","VR array overrun.");
    
    I->V = ReallocForSure(I->V,float,(v-I->V));
    I->VR = ReallocForSure(I->VR,float,(vr-I->VR));
    if(I->VSP) 
      I->VSP = ReallocForSure(I->VSP,float,(vsp-I->VSP));
    if(I->VSPC) 
      I->VSPC = ReallocForSure(I->VSPC,float,(vspc-I->VSPC));
    if(I->VarAlpha)
      I->VarAlpha = ReallocForSure(I->VarAlpha,float,n_var_alpha);
    if(n_var_alpha_ray) {
      if(I->VarAlphaRay)
        I->VarAlphaRay = ReallocForSure(I->VarAlphaRay,float,n_var_alpha_ray);
    } else {
      FreeP(I->VarAlphaRay);
    }
    if(n_var_alpha_sph) {
      if(I->VarAlphaSph)
        I->VarAlphaSph = ReallocForSure(I->VarAlphaSph,float,n_var_alpha_sph);
    } else {
      FreeP(I->VarAlphaSph);
    }
      
    if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) { 

      PRINTFD(G,FB_RepCylBond)
        " RepCylBondNEW: generating pickable version\n"
        ENDFD;
      /* pickable versions are simply capped boxes, 
         vertices: 8 points * 3 = 32  * 2 = 48 floats per bond
      */
       
      vp_size = maxCyl*24;
      I->VP=Alloc(float,vp_size);
      ErrChkPtr(G,I->VP);
       
      rp_size = maxCyl+1;
      I->R.P=Alloc(Pickable,rp_size);
      ErrChkPtr(G,I->R.P);
      rp = I->R.P + 1; /* skip first record! */
      I->R.context.object = (void*)obj;
      I->R.context.state = state;
       
      v=I->VP;
      b=obj->Bond;
      for(a=0;a<obj->NBond;a++)     {
        b1 = b->index[0];
        b2 = b->index[1];
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
        if((a1>=0)&&(a2>=0))	{
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

          if(s1||s2) {	
            float bd_radius;
            float overlap_r,nub_r;
                    
            AtomInfoGetBondSetting_f(G,b,cSetting_stick_radius,radius,&bd_radius);
                    
            overlap_r = overlap * bd_radius;
            nub_r = nub * bd_radius;
                    
            copy3f(cs->Coord+3*a1,v1);
            copy3f(cs->Coord+3*a2,v2);
						
            h[0]=(v1[0]+v2[0])/2;
            h[1]=(v1[1]+v2[1])/2;
            h[2]=(v1[2]+v2[2])/2;
						
            if(s1&(!ai1->masked)) {
              I->NP++;
              rp->index = b1;
              rp->bond = a;
              rp++;
                             
              v = RepCylinderBox(v,v1,h,bd_radius,overlap_r,nub_r);
            }
            if(s2&(!ai2->masked))  {
              I->NP++;
              rp->index = b2;
              rp->bond = a;
              rp++;
                             
              v = RepCylinderBox(v,h,v2,bd_radius,overlap_r,nub_r);
            }
          }
        }
        b++;
      }
      
      if((signed)vp_size<(v-I->VP))
        ErrFatal(G,"RepCylBond","VP array overrun.");
      if((signed)rp_size<=(I->NP))
        ErrFatal(G,"RepCylBond","RP array overrun.");

      I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
      I->R.P[0].index = I->NP;
      I->VP = ReallocForSure(I->VP,float,(v-I->VP));
      
      PRINTFD(G,FB_RepCylBond)
        " RepCylBondNew: I->NP: %d I->VP: %p\n",I->NP,
        (void*)I->VP
        ENDFD;
    }
  }
  FreeP(other);
  FreeP(marked);
  return((void*)(struct Rep*)I);
}

static void subdivide( int n, float *x, float *y)
{
  int a;
  if(n<3) {n=3;}
  for(a=0;a<=n;a++) {
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
  
  /*  t[0] = d[1];
  t[1] = d[2];
  t[2] = -d[0];
  */
  get_divergent3f(d,t);

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


float *RepCylinder(float *v,float *v1,float *v2,int nEdge,
                   int frontCap, int endCap,
                   float tube_size,float overlap,float nub)
{

  float d[3],t[3],p0[3],p1[3],p2[3];
  float x[51],y[51];
  int c;

  if(nEdge>50)
    nEdge=50;

  subdivide(nEdge,x,y); /* NOTE: lift this out...*/
 
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
  
  /*  t[0] = d[1];
  t[1] = d[2];
  t[2] = -d[0];*/
  
  get_divergent3f(d,t);

  cross_product3f(d,t,p1);
  
  normalize3f(p1);
  
  cross_product3f(d,p1,p2);
  
  normalize3f(p2);
  
  /* now we have a coordinate system*/
  
  for(c=nEdge;c>=0;c--) {
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

static void RepCylinderImmediate(float *v1,float *v2,int nEdge,
                                 int frontCap, int endCap,
                                 float overlap,float nub,
                                 float *x, float *y)
{
  float d[3],t[3],p0[3],p1[3],p2[3];
  float vv1[3],vv2[3];
  float v[3],vv[3],vvv[3];
  int c;

  /* direction vector */
  
  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);
  
  normalize3f(p0);

  vv1[0]=v1[0] - p0[0]*overlap;
  vv1[1]=v1[1] - p0[1]*overlap;
  vv1[2]=v1[2] - p0[2]*overlap;
  
  v1 = vv1;

  if(endCap) {
    vv2[0]=v2[0] + p0[0]*overlap;
    vv2[1]=v2[1] + p0[1]*overlap;
    vv2[2]=v2[2] + p0[2]*overlap;
    v2 = vv2;
  }
  
  d[0] = (v2[0] - v1[0]);
  d[1] = (v2[1] - v1[1]);
  d[2] = (v2[2] - v1[2]);
  
  get_divergent3f(d,t);

  cross_product3f(d,t,p1);
  
  normalize3f(p1);
  
  cross_product3f(d,p1,p2);
  
  normalize3f(p2);
  
  /* now we have a coordinate system*/
  
 glBegin(GL_TRIANGLE_STRIP);
 
 for(c=nEdge;c>=0;c--) {
   v[0] = p1[0]*x[c] + p2[0]*y[c];
   v[1] = p1[1]*x[c] + p2[1]*y[c];
   v[2] = p1[2]*x[c] + p2[2]*y[c];

   vv[0] = v1[0] + v[0];
   vv[1] = v1[1] + v[1];
   vv[2] = v1[2] + v[2];

   glNormal3fv(v);

   vvv[0] = vv[0] + d[0];
   vvv[1] = vv[1] + d[1];
   vvv[2] = vv[2] + d[2];

   glVertex3fv(vv);
   glVertex3fv(vvv);
 }
 glEnd();

 if(frontCap) {
   
    glBegin(GL_TRIANGLE_FAN);

    v[0] = -p0[0];
    v[1] = -p0[1];
    v[2] = -p0[2];
    
    vv[0] = v1[0] - p0[0]*nub;
    vv[1] = v1[1] - p0[1]*nub;
    vv[2] = v1[2] - p0[2]*nub;

    glNormal3fv(v);
    glVertex3fv(vv);

    for(c=nEdge;c>=0;c--) {
      
      v[0] = p1[0]*x[c] + p2[0]*y[c];
      v[1] = p1[1]*x[c] + p2[1]*y[c];
      v[2] = p1[2]*x[c] + p2[2]*y[c];
      
      vv[0] = v1[0] + v[0];
      vv[1] = v1[1] + v[1];
      vv[2] = v1[2] + v[2];
      
      glNormal3fv(v);
      glVertex3fv(vv);
    }

    glEnd();
 }
 
  if(endCap) {

    v[0] = p0[0];
    v[1] = p0[1];
    v[2] = p0[2];
	
    vv[0] = v2[0] + p0[0]*nub;
    vv[1] = v2[1] + p0[1]*nub;
    vv[2] = v2[2] + p0[2]*nub;

    glBegin(GL_TRIANGLE_FAN);
    
    glNormal3fv(v);
    glVertex3fv(vv);
	
    for(c=0;c<=nEdge;c++) {
      
      v[0] = p1[0]*x[c] + p2[0]*y[c];
      v[1] = p1[1]*x[c] + p2[1]*y[c];
      v[2] = p1[2]*x[c] + p2[2]*y[c];
      
      vv[0] = v2[0] + v[0];
      vv[1] = v2[1] + v[1];
      vv[2] = v2[2] + v[2];
      glNormal3fv(v);
      glVertex3fv(vv);
    }
    glEnd();
  }
}

void RepCylBondRenderImmediate(CoordSet *cs, RenderInfo *info)
{
  /* performance optimized, so it does not support the following:

  - anything other than opengl
  - display of bond valences
  - per-bond & per-atom properties
  - half-bonds
  - helper settings such as cartoon_side_chain_helper
  - suppression of long bonds
  - color ramps
  - atom picking
  - display lists
  - transparency 
  
  */

  PyMOLGlobals *G=cs->State.G;
    if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)) )
    return;
  else {
    int active = false;
    ObjectMolecule *obj = cs->Obj;
    int nEdge = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_stick_quality);
    float radius = fabs(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_radius));
    float overlap = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_overlap);
    float nub = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_stick_nub);
    float overlap_r = radius*overlap;
    float nub_r = radius*nub;
    float x[51],y[51];
    
    if(nEdge>50)
      nEdge=50;
    
    subdivide(nEdge,x,y);  
    {
      int c;
      for(c=0;c<=nEdge;c++) {
        x[c] *= radius;
        y[c] *= radius;
      }
    }

    {
      int a;
      int nBond = obj->NBond;
      BondType *bd = obj->Bond;
      AtomInfoType *ai = obj->AtomInfo;
      int *atm2idx = cs->AtmToIdx;
      int discreteFlag = obj->DiscreteFlag;
      int last_color = -9;
      float *coord = cs->Coord;
      const float _pt5 = 0.5F;

      for(a=0;a<nBond;a++) {
        int b1 = bd->index[0];
        int b2 = bd->index[1];
        AtomInfoType *ai1, *ai2;
        bd++;
        
        if( (ai1 = ai+b1)->visRep[cRepCyl] && (ai2 = ai+b2)->visRep[cRepCyl]) {
          int a1, a2;
          active = true;
          if(discreteFlag) {
            /* not optimized */
            if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
              a1=obj->DiscreteAtmToIdx[b1];
              a2=obj->DiscreteAtmToIdx[b2];
            } else {
              a1=-1;
              a2=-1;
            }
          } else {
            a1=atm2idx[b1];
            a2=atm2idx[b2];
          }
          
          if((a1>=0)&&(a2>=0)) {
            int c1 = ai1->color;
            int c2 = ai2->color;
            
            float *v1 = coord+3 * a1;
            float *v2 = coord+3 * a2;
            
            if(c1 == c2) { /* same colors -> one cylinder */
              if(c1!=last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G,c1));
              }
              
              RepCylinderImmediate(v1,v2,nEdge,1,1,
                                   overlap_r,nub_r,x,y);
              
            } else { /* different colors -> two cylinders, no interior */
              float avg[3];
              
              avg[0] = (v1[0]+v2[0]) * _pt5;
              avg[1] = (v1[1]+v2[1]) * _pt5;
              avg[2] = (v1[2]+v2[2]) * _pt5;
              
              if(c1!=last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G,c1));
              }

              RepCylinderImmediate(v1,avg,nEdge,1,0,
                                   overlap_r,nub_r,x,y);
              
              if(c2!=last_color) {
                last_color = c2;
                glColor3fv(ColorGet(G,c2));
              }

              RepCylinderImmediate(v2,avg,nEdge,1,0,
                                   overlap_r,nub_r,x,y);
            }
          }
        }
      }
    }
    if(!active)
      cs->Active[cRepCyl] = false;
  }
}


