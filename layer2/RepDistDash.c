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
#include"RepDistDash.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"

typedef struct RepDistDash {
  Rep R;
  float *V;
  int N;
  CObject *Obj;
  DistSet *ds;
  float linewidth,radius;
} RepDistDash;

#include"ObjectDist.h"

void RepDistDashFree(RepDistDash *I);

void RepDistDashFree(RepDistDash *I)
{
  VLAFreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}


static void RepDistDashRender(RepDistDash *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  float *vc;
  int round_ends;
  int color = SettingGet_color(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_color);
  float line_width = SettingGet_f(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_width);
                                             
  I->radius = SettingGet_f(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_radius);
  round_ends = SettingGet_b(G,I->ds->Setting,I->ds->Obj->Obj.Setting,cSetting_dash_round_ends);
line_width = SceneGetDynamicLineWidth(info,line_width);
  
  if(ray) {

    float radius;

    if(I->radius<=0.0F) {
      radius = ray->PixelRadius*line_width/2.0F;
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
        glLineWidth(line_width * info->width_scale);
      } else {
        glLineWidth(line_width);
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
      
        if(!info->line_lighting)
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

Rep *RepDistDashNew(DistSet *ds)
{
  PyMOLGlobals *G=ds->State.G;
  int a;
  int n;
  float *v,*v1,*v2,d[3],d1[3];
  float l;
  float dash_len,dash_gap,dash_sum;
#if 0
  float d2[3], ph, seg;
#endif

  OOAlloc(G,RepDistDash);

  if(!ds->NIndex) {
    OOFreeP(I);
    return(NULL); 
  }

  RepInit(G,&I->R);
  
  
  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepDistDashRender;
  I->R.fFree=(void (*)(struct Rep *))RepDistDashFree;
  I->R.fRecolor=NULL;

  dash_len = SettingGet_f(G,ds->Setting,ds->Obj->Obj.Setting,cSetting_dash_length);
  dash_gap = SettingGet_f(G,ds->Setting,ds->Obj->Obj.Setting,cSetting_dash_gap);
  dash_sum = dash_len+dash_gap;
  if(dash_sum<R_SMALL4) dash_sum=0.5;

  I->N=0;
  I->V=NULL;
  I->R.P=NULL;
  I->Obj = (CObject*)ds->Obj;
  I->ds = ds;

  n=0;
  if(ds->NIndex) {
	 I->V=VLAlloc(float,ds->NIndex*10);

	 for(a=0;a<ds->NIndex;a=a+2) {
      v1 = ds->Coord+3*a;
      v2 = ds->Coord+3*(a+1);

      subtract3f(v2,v1,d);

      l = (float)length3f(d);

      if(l>R_SMALL4) {
        
        copy3f(v1,d1);
        normalize3f(d);
        
        if(dash_gap>R_SMALL4) {
          float avg[3],proj1[3],proj2[3];
          float l_left = l/2.0F;
          float l_used = 0.0F;
          float half_dash_gap = dash_gap*0.5;
          
          average3f(v1,v2,avg);
          while(l_left>dash_sum) {
            VLACheck(I->V,float,(n*3)+11);
            v=I->V+n*3;
            scale3f(d,l_used+half_dash_gap,proj1);
            scale3f(d,l_used+dash_len+half_dash_gap,proj2);
            add3f(avg,proj1,v);
            add3f(avg,proj2,v+3);
            subtract3f(avg,proj1,v+6);
            subtract3f(avg,proj2,v+9);
            n+=4;
            l_left-=dash_sum;
            l_used += dash_sum;
          }
          if(l_left>dash_gap) {
            l_left-=dash_gap;
            scale3f(d,l_used+half_dash_gap,proj1);
            scale3f(d,l_used+l_left+half_dash_gap,proj2);
            VLACheck(I->V,float,(n*3)+11);
            v=I->V+n*3;
            add3f(avg,proj1,v);
            add3f(avg,proj2,v+3);
            subtract3f(avg,proj1,v+6);
            subtract3f(avg,proj2,v+9);
            n+=4;
          }
        } else if(dash_len>R_SMALL4) {
          VLACheck(I->V,float,(n*3)+5);
          v=I->V+n*3;
          copy3f(v1,v);
          copy3f(v2,v+3);
          n+=2;
        }
      }
#if 0
      l -= dash_gap;
      
      ph = dash_sum-(float)fmod((l+dash_gap)/2.0,dash_sum);
      if(l>R_SMALL4) {

        copy3f(v1,d1);
        normalize3f(d);
        scale3f(d,dash_gap,d2);        
        scale3f(d2,0.5F,d2);
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
#endif

    }
    VLASize(I->V,float,n*3);
    I->N=n;
  }
  return((void*)(struct Rep*)I);
}


