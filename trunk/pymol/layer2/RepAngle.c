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


#include"OOMac.h"
#include"RepAngle.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"

typedef struct RepAngle {
  Rep R;
  float *V;
  int N;
  CObject *Obj;
  DistSet *ds;
  float linewidth,radius;
} RepAngle;

#include"ObjectDist.h"

void RepAngleFree(RepAngle *I);

void RepAngleFree(RepAngle *I)
{
  VLAFreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepAngleRender(RepAngle *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  float *vc;
  int round_ends;
  int color = SettingGet_color(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_angle_color);
  I->linewidth = SettingGet_f(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_width);
  I->radius = SettingGet_f(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_radius);
  round_ends = SettingGet_b(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_round_ends);

  if(ray) {

    float radius;

    if(I->radius==0.0F) {
      radius = ray->PixelRadius*I->linewidth/2.0F;
    } else {
      radius = I->radius;
    }

    if(color<0)
      color = I->Obj->Color;
    vc = ColorGet(G,color);
	 v=I->V;
	 c=I->N;
	 
	 while(c>0) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]);*/
      if(round_ends) {
        ray->fSausage3fv(ray,v,v+3,radius,vc,vc);
      } else {
        ray->fCustomCylinder3fv(ray,v,v+3,radius,vc,vc,cCylCapFlat,cCylCapFlat);
      }
		v+=6;
      c-=2;
	 }

  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      int use_dlst;

      if(info->width_scale_flag) {
        glLineWidth(I->linewidth * info->width_scale);
      } else {
        glLineWidth(I->linewidth);
      }
      if(color>=0)
        glColor3fv(ColorGet(G,color));

      use_dlst = (int)SettingGet(G,cSetting_use_display_lists);
      if(use_dlst&&I->R.displayList) {
        glCallList(I->R.displayList);
      } else { 

        SceneResetNormal(G,true);

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
      
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);	 
        while(c>0) {
          glVertex3fv(v);
          v+=3;
          glVertex3fv(v);
          v+=3;
          c-=2;
        }
        glEnd();
        glEnable(GL_LIGHTING);

        glEnable(GL_LIGHTING);
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
      }
    }
  }
}

Rep *RepAngleNew(DistSet *ds)
{
  PyMOLGlobals *G=ds->State.G;
  int a;
  int n = 0;
  float *v,*v1,*v2,*v3,*v4,d1[3],d2[3],d3[3],n1[3],n3[3],l1,l2,x[3],y[3];
  float length,radius,angle,pos,phase;
  float dash_len,dash_gap,dash_sum;

  OOAlloc(G,RepAngle);

  PRINTFD(G,FB_RepAngle)
    "RepAngleNew: entered.\n"
    ENDFD;
  if(!ds->NAngleIndex) {
    OOFreeP(I);
    return(NULL); 
  }

  RepInit(G,&I->R);
  
  
  I->R.fRender=(void (*)(struct Rep *, RenderInfo *info))RepAngleRender;
  I->R.fFree=(void (*)(struct Rep *))RepAngleFree;
  I->R.fRecolor=NULL;

  dash_len = SettingGet_f(G,ds->Setting,ds->Obj->Obj.Setting,cSetting_dash_length);
  dash_gap = SettingGet_f(G,ds->Setting,ds->Obj->Obj.Setting,cSetting_dash_gap);
  dash_sum = dash_len+dash_gap;
  if(dash_sum<R_SMALL4) dash_sum=0.1F;

  I->N=0;
  I->V=NULL;
  I->R.P=NULL;
  I->Obj = (CObject*)ds->Obj;
  I->ds = ds;

  n=0;
  if(ds->NAngleIndex) {
	 I->V=VLAlloc(float,ds->NAngleIndex*10);

	 for(a=0;a<ds->NAngleIndex;a=a+5) {
      v1 = ds->AngleCoord+3*a;
      v2 = ds->AngleCoord+3*(a+1);
      v3 = ds->AngleCoord+3*(a+2);
      v4 = ds->AngleCoord+3*(a+3);
      subtract3f(v1,v2,d1);
      subtract3f(v3,v2,d2);

      l1 = (float)length3f(d1);
      l2 = (float)length3f(d2);

      if(l1>l2)
        radius = l2;
      else
        radius = l1;
      radius *= SettingGet_f(G,ds->Setting,ds->Obj->Obj.Setting,cSetting_angle_size);
      
      angle = get_angle3f(d1,d2);

      normalize23f(d1,n1);

      remove_component3f(d2,n1,d3);

      if(length3f(d3)<R_SMALL8) {
        d3[0]=1.0F;
        d3[1]=0.0F;
        d3[2]=0.0F;
      } else {
        normalize23f(d3,n3);
      }

      scale3f(n1,radius,x);
      scale3f(n3,radius,y);

      if(v4[0]!=0.0F) { /* line 1 flag */
        
        VLACheck(I->V,float,(n*3)+5);
        v=I->V+n*3;
        copy3f(v1,v);
        v+=3;
        copy3f(v2,v);
        n+=2;
      }

      if(v4[1]!=0.0F) { /* line 2 flag */

        VLACheck(I->V,float,(n*3)+5);
        v=I->V+n*3;
        copy3f(v3,v);
        v+=3;
        copy3f(v2,v);
        n+=2;
      }

      /* now we have a relevant orthogonal axes */

      length = (float)(angle * radius * 2); 

      /* figure out dash/gap phasing that will lead to nicely space dashes and gaps */

      phase = dash_sum - (float)fmod(length/2+(dash_gap/2), dash_sum); 
      pos = -phase;

      if(length>R_SMALL4) {
        
        float mod_pos;
        float vx[3],vy[3];
        float cur_angle;
        float cons_pos1, cons_pos2;

        while(pos<length) {

          mod_pos = (float)fmod(pos + phase, dash_sum);

          VLACheck(I->V,float,(n*3)+5);
          
          cons_pos1 = pos;
          if(cons_pos1<0.0F) cons_pos1 = 0.0F;
          cons_pos2 = pos + dash_len;
          if(cons_pos2>length) cons_pos2 = length;
          
          if(cons_pos1<cons_pos2) {
            cur_angle = angle * cons_pos1/length;
            
            v=I->V+n*3;
            scale3f(x,(float)cos(cur_angle),vx);
            scale3f(y,(float)sin(cur_angle),vy);
            add3f(vx,vy,v);
            add3f(v2,v,v);
            
            cur_angle = angle * cons_pos2/length;
            
            v+=3;
            scale3f(x,(float)cos(cur_angle),vx);
            scale3f(y,(float)sin(cur_angle),vy);
            add3f(vx,vy,v);
            add3f(v2,v,v);
            
            n+=2;
          }
          pos+=dash_sum;
        }
      }
    }
    VLASize(I->V,float,n*3);
    I->N=n;
  }
  return((void*)(struct Rep*)I);
}


