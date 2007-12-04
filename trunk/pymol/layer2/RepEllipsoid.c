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
#include"RepEllipsoid.h"
#include"Color.h"
#include"Setting.h"
#include"Feedback.h"
#include"Matrix.h"
#include"CGO.h"

typedef struct RepEllipsoid {
  Rep R; /* must be first! */
  CGO *ray,*std;
} RepEllipsoid;

#include"ObjectMolecule.h"

void RepEllipsoidFree(RepEllipsoid *I);

void RepEllipsoidFree(RepEllipsoid *I)
{
  if(I->ray)
    CGOFree(I->ray);
  if(I->std)
    CGOFree(I->std);
  RepPurge(&I->R);
  OOFreeP(I);
}


static void RepEllipsoidRender(RepEllipsoid *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  register PyMOLGlobals *G=I->R.G;
  if(ray) {
    PRINTFD(G,FB_RepEllipsoid)
      " RepEllipsoidRender: rendering ray...\n"
      ENDFD;
    
    if(I->ray)  
      CGORenderRay(I->ray,ray,NULL,I->R.cs->Setting,
                   I->R.obj->Setting);
    else if(I->std)
      CGORenderRay(I->std,ray,NULL,I->R.cs->Setting,
                   I->R.obj->Setting);    
  } else if(G->HaveGUI && G->ValidContext) {

    if(pick) {
      if(I->std) {
        CGORenderGLPicking(I->std,pick,&I->R.context,
                            I->R.cs->Setting,I->R.obj->Setting);
      }
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
        
        PRINTFD(G,FB_RepEllipsoid)
          " RepEllipsoidRender: rendering GL...\n"
          ENDFD;
        if(I->std) 
          CGORenderGL(I->std,NULL,I->R.cs->Setting,
                      I->R.obj->Setting,info);
        
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
      }
    }
  }
}

const double problevel[50] = {0.4299, 0.5479, 0.6334, 0.7035, 0.7644, 
                              0.8192, 0.8694, 0.9162, 0.9605, 1.0026,
                              1.0430, 1.0821, 1.1200, 1.1570, 1.1932,
                              1.2288, 1.2638, 1.2985, 1.3330, 1.3672,
                              1.4013, 1.4354, 1.4695, 1.5037, 1.5382,
                              1.5729, 1.6080, 1.6436, 1.6797, 1.7164,
                              1.7540, 1.7924, 1.8318, 1.8724, 1.9144,
                              1.9580, 2.0034, 2.0510, 2.1012, 2.1544,
                              2.2114, 2.2730, 2.3404, 2.4153, 2.5003,
                              2.5997, 2.7216, 2.8829, 3.1365, 6.0000 };

Rep *RepEllipsoidNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;

  OOCalloc(G,RepEllipsoid); /* allocates & sets I */
  
  obj = cs->Obj;

  {
    int visible_flag = false;
    int a;
    visible_flag=false;
    if(obj->RepVisCache[cRepEllipsoid])
      for(a=0;a<cs->NIndex;a++) {
        if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepEllipsoid]) {
          visible_flag=true;
          break;
        }
      }
    if(!visible_flag) {
      OOFreeP(I);
      return(NULL); /* skip if no dots are visible */
    }
  }

  RepInit(G,&I->R);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo *))RepEllipsoidRender;
  I->R.fFree = (void (*)(struct Rep *))RepEllipsoidFree;
  I->R.cs = cs;
  I->R.obj = (CObject*)obj;
  I->R.context.object = (void*)obj;
  I->R.context.state = state;

  /*  I->R.fSameVis=(int (*)(struct Rep*, struct CoordSet*))RepEllipsoidSameVis;*/

  {
    int ellipsoid_color = SettingGet_color(G,cs->Setting,obj->Obj.Setting,
                                         cSetting_ellipsoid_color);
    
    int cartoon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                                 cSetting_cartoon_side_chain_helper);

    int ribbon_side_chain_helper = SettingGet_b(G,cs->Setting, obj->Obj.Setting,
                                                cSetting_ribbon_side_chain_helper);

    float ellipsoid_scale = SettingGet_f(G,cs->Setting,obj->Obj.Setting,
                                         cSetting_ellipsoid_scale);

    float transp = SettingGet_f(G,cs->Setting,obj->Obj.Setting,
                                cSetting_ellipsoid_transparency);

    int pickable = SettingGet_b(G,cs->Setting,obj->Obj.Setting,
                                cSetting_pickable);

    float prob = SettingGet_f(G,cs->Setting,obj->Obj.Setting,
                                cSetting_ellipsoid_probability);
    double matrix_factor = 0.0F;
    float pradius = 0.0F;
      {

      int iprob = (prob+0.01F)*50.0F - 1;
      if(iprob<0) iprob = 0;
      if(iprob>49) iprob = 49;
      pradius = problevel[iprob];
      matrix_factor = -(1/(pradius*pradius));
    }

    I->ray = CGONew(G); /* describe the ellipsoids analytically */

    if(I->ray) {
      int a, a1;
      int vis_flag;
      AtomInfoType *ai;
      float last_alpha = 1.0F;
      
      for(a=0;a<cs->NIndex;a++) {
        a1 = cs->IdxToAtm[a];
        ai = obj->AtomInfo+a1;
        vis_flag = ai->visRep[cRepEllipsoid];
        
        if(vis_flag &&
           (!ai->hetatm) &&
           ((cartoon_side_chain_helper && ai->visRep[cRepCartoon]) ||
            (ribbon_side_chain_helper && ai->visRep[cRepRibbon]))) {
          
          register char *name1=ai->name;
          register int prot1=ai->protons;


          if(prot1 == cAN_N) { 
            if((!name1[1])&&(name1[0]=='N')) { /* N */
              register char *resn1 = ai->resn;
              if(!((resn1[0]=='P')&&(resn1[1]=='R')&&(resn1[2]=='O')))
                vis_flag=false;
            }
          } else if(prot1 == cAN_O) { 
            if((!name1[1])&&(name1[0]=='O'))
              vis_flag=false;
          } else if(prot1 == cAN_C) {
            if((!name1[1])&&(name1[0]=='C'))
              vis_flag=false;
          }
        }
        
        if(vis_flag) {
          double u11,u22,u33,u12,u13,u23;
          
          u11 = ai->U11;
          u22 = ai->U22;
          u33 = ai->U33;
          u12 = ai->U12;
          u13 = ai->U13;
          u23 = ai->U23;
          if (u11 || u22 || u33 || u12 || u13 || u23) {
            int n_rot;
            double matrix[16];
            double e_val[4];
            double e_vec[16];
              
            matrix[0] = u11;
            matrix[1] = u12;
            matrix[2] = u13;
            matrix[3] = 0.0;
            matrix[4] = u12;
            matrix[5] = u22;
            matrix[6] = u23;
            matrix[7] = 0.0;
            matrix[8] = u13;
            matrix[9] = u23;
            matrix[10] = u33;
            matrix[11] = 0.0;
            matrix[12] = 0.0;
            matrix[13] = 0.0;
            matrix[14] = 0.0;
            matrix[15] = matrix_factor;
            
            if(xx_matrix_jacobi_solve(e_vec, e_val, &n_rot, matrix, 4)) {

              float at_ellipsoid_scale;
              int at_ellipsoid_color;
              float at_transp;
              float *v = cs->Coord+3*a;          

              float mag[3];
              float scale[3];

              float mx;
              float r_el,n0[3],n1[3],n2[3];
              int c1;

              AtomInfoGetSetting_f(G, ai, cSetting_ellipsoid_scale, ellipsoid_scale, &at_ellipsoid_scale);
              AtomInfoGetSetting_f(G, ai, cSetting_ellipsoid_transparency, transp, &at_transp);
              AtomInfoGetSetting_color(G, ai, cSetting_ellipsoid_color, ellipsoid_color, &at_ellipsoid_color);
              
              if(at_ellipsoid_color==-1)
                c1=*(cs->Color+a);
              else
                c1=at_ellipsoid_color;

              n0[0] = e_vec[0];
              n0[1] = e_vec[4];
              n0[2] = e_vec[8];
              n1[0] = e_vec[1];
              n1[1] = e_vec[5];
              n1[2] = e_vec[9];
              n2[0] = e_vec[2];
              n2[1] = e_vec[6];
              n2[2] = e_vec[10];
              
              normalize3f(n0);
              normalize3f(n1);
              normalize3f(n2);
              mag[0] = sqrt1f(e_val[0]);
              mag[1] = sqrt1f(e_val[1]);
              mag[2] = sqrt1f(e_val[2]);
              
              mx = mag[0];
              if( mx < mag[1]) mx = mag[1];
              if( mx < mag[2]) mx = mag[2];
              
              scale[0] = mag[0]/mx;
              scale[1] = mag[1]/mx;
              scale[2] = mag[2]/mx;
              
              scale3f(n0,scale[0],n0);
              scale3f(n1,scale[1],n1);
              scale3f(n2,scale[2],n2);
              
              r_el = mx * pradius * ellipsoid_scale;
              
              {
                float vc[3];
                if(ColorCheckRamped(G,c1)) {
                  ColorGetRamped(G,c1,v,vc,state);
                  CGOColorv(I->ray, vc);
                } else {
                  CGOColorv(I->ray, ColorGet(G,c1));
                }
              }

              { 
                float alpha = 1.0F - at_transp;
                if(alpha != last_alpha) {
                  CGOAlpha(I->ray,alpha);
                  last_alpha = alpha;
                }
              }
              if(pickable && (!ai->masked))
                CGOPickColor(I->ray,a1,cPickableAtom);

              CGOEllipsoid(I->ray,v,r_el,n0,n1,n2);
            }
          }
        }
      }
      CGOStop(I->ray);
      I->std = CGOSimplify(I->ray,0); /* convert analytical to discrete */
    }
  }
  return((void*)(struct Rep*)I);
}



