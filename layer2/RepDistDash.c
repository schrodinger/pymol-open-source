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

#include"os_std.h"
#include"os_gl.h"


#include"OOMac.h"
#include"RepDistDash.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"Object.h"

typedef struct RepDistDash {
  Rep R;
  float *V;
  int N;
  Object *Obj;
  float linewidth,radius;
} RepDistDash;

#include"ObjectDist.h"

void RepDistDashRender(RepDistDash *I,CRay *ray,Pickable **pick);
void RepDistDashFree(RepDistDash *I);

void RepDistDashFree(RepDistDash *I)
{
  VLAFreeP(I->V);
  RepFree(&I->R);
  OOFreeP(I);
}


void RepDistDashRender(RepDistDash *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  float *vc;

  if(ray) {

    vc = ColorGet(I->Obj->Color);
	 v=I->V;
	 c=I->N;
	 
	 while(c>0) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]);*/
		ray->fCylinder3fv(ray,v,v+3,I->radius,vc,vc);
		v+=6;
      c-=2;
	 }

  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 
	 v=I->V;
	 c=I->N;
	 
    glDisable(GL_LIGHTING);
    glLineWidth(I->linewidth);
	 glBegin(GL_LINES);	 
	 SceneResetNormal(true);
	 while(c>0) {
		glVertex3fv(v);
		v+=3;
		glVertex3fv(v);
		v+=3;
      c-=2;
	 }
	 glEnd();
    glEnable(GL_LIGHTING);
  }
}

Rep *RepDistDashNew(DistSet *ds)
{
  int a;
  int n;
  float *v,*v1,*v2,d[3],d1[3],d2[3];
  float l,ph;
  float dash_len,dash_gap,dash_sum,seg;
  OOAlloc(RepDistDash);

  if(!ds->NIndex) {
    OOFreeP(I);
    return(NULL); 
  }

  RepInit(&I->R);


  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepDistDashRender;
  I->R.fFree=(void (*)(struct Rep *))RepDistDashFree;
  I->R.fRecolor=NULL;

  dash_len = SettingGet(cSetting_dash_length);
  dash_gap = SettingGet(cSetting_dash_gap);
  dash_sum = dash_len+dash_gap;
  if(dash_sum<R_SMALL4) dash_sum=0.5;

  I->linewidth = SettingGet(cSetting_dash_width);
  I->radius = SettingGet(cSetting_dash_radius);

  I->N=0;
  I->V=NULL;
  I->R.P=NULL;
  I->Obj = (Object*)ds->Obj;

  n=0;
  if(ds->NIndex) {
	 I->V=VLAlloc(float,ds->NIndex*10);

	 for(a=0;a<ds->NIndex;a=a+2) {
      v1 = ds->Coord+3*a;
      v2 = ds->Coord+3*(a+1);

      subtract3f(v2,v1,d);

      l = length3f(d);
      
      l -= dash_gap;
      
      ph = dash_sum-fmod((l+dash_gap)/2.0,dash_sum);
      if(l>R_SMALL4) {

        copy3f(v1,d1);
        normalize3f(d);
        scale3f(d,dash_gap,d2);        
        scale3f(d2,0.5,d2);
        add3f(d1,d2,d1);

        while(l>0.0) {
          if(ph<dash_len) {
            seg = dash_len-ph;
            if(l<seg) seg = l;
            scale3f(d,seg,d2);
            if((seg/dash_len)>0.2) {
              VLACheck(I->V,float,(n*3)+5);
              v=I->V+n*3;
              copy3f(d1,v);
              v+=3;
              add3f(d1,d2,d1);
              copy3f(d1,v);
              n+=2;
            } else 
              add3f(d1,d2,d1);
            ph = dash_len;
          } else {
            /* gap */
            seg = dash_gap;
            if(l<seg) seg = l;            
            scale3f(d,seg,d2);
            add3f(d1,d2,d1);
            ph = 0.0;
          }
          l-=seg;
        }
      }
    }
    VLASize(I->V,float,n*3);
    I->N=n;
  }
  return((void*)(struct Rep*)I);
}


