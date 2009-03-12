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
#include"RepLabel.h"
#include"Color.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"Scene.h"
#include"Text.h"

typedef struct RepLabel {
  Rep R;
  float *V;
  int *L;
  int N;
  int OutlineColor;
} RepLabel;

#include"ObjectMolecule.h"

void RepLabelFree(RepLabel *I);

void RepLabelInit(void)
{
}

void RepLabelFree(RepLabel *I)
{
  FreeP(I->R.P);
  FreeP(I->V);
  FreeP(I->L);
  OOFreeP(I);
}

static void RepLabelRender(RepLabel *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  int *l=I->L;
  int font_id = SettingGet_i(G,I->R.cs->Setting,I->R.obj->Setting,
                             cSetting_label_font_id);
  float font_size = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,
                                 cSetting_label_size);

  if(ray) {
    if(c) {
      char *st;
      TextSetOutlineColor(G,I->OutlineColor);
      while(c--) {
        if(*l) {
          st = OVLexicon_FetchCString(G->Lexicon,*l);
          TextSetPosNColor(G,v+3,v);
          TextRenderRay(G,ray,font_id,st,font_size,v+6);
        }
        v+=9;
        l++;
      }
    }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      Pickable *p = I->R.P;
      int i;
      if(c) {
        char *st;
        int float_text = (int)SettingGet(G,cSetting_float_labels);
        if(float_text) 
          glDisable(GL_DEPTH_TEST);	 

        i=(*pick)->src.index;
        while(c--) {
          if(*l) {
            int first_pass = (!(*pick)[0].src.bond);
            i++;            
            TextSetPosNColor(G,v+3,v);
            TextSetPickColor(G, first_pass, i);
            if(first_pass) {
              VLACheck((*pick),Picking,i);
              p++;
              (*pick)[i].src = *p; /* copy object and atom info */
              (*pick)[i].context = I->R.context;
            }
            st = OVLexicon_FetchCString(G->Lexicon,*l);
            TextRenderOpenGL(G,info,font_id,st,font_size,v+6);
          }
          l++;
          v+=9;
        }
        if(float_text)
          glEnable(GL_DEPTH_TEST);	 
        (*pick)[0].src.index = i; /* pass the count */
      }
    } else {
      if(c) {
        char *st;
        int float_text = (int)SettingGet(G,cSetting_float_labels);
        if(float_text)
          glDisable(GL_DEPTH_TEST);	 
        if(!info->line_lighting)
          glDisable(GL_LIGHTING); 
        TextSetOutlineColor(G,I->OutlineColor);
        while(c--) {
          if(*l) {
            TextSetPosNColor(G,v+3,v);
            st = OVLexicon_FetchCString(G->Lexicon,*l);
            TextRenderOpenGL(G,info,font_id,st,font_size,v+6);
          }
          l++;
          v+=9;
        }
        glEnable(GL_LIGHTING);
        glEnable(GL_BLEND);
        if(float_text)
          glEnable(GL_DEPTH_TEST);	 
      }
    }
  }
}

Rep *RepLabelNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  int a,a1,vFlag,c1;
  float *v,*v0,*vc;
  float *lab_pos;
  int *l;
  int label_color;
  LabPosType *lp = NULL;
  Pickable *rp = NULL;
  AtomInfoType *ai;
  OOAlloc(G,RepLabel);
  
  obj = cs->Obj;
  vFlag=false;
  if(obj->RepVisCache[cRepLabel])
    for(a=0;a<cs->NIndex;a++) {
      if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepLabel])
        {
          vFlag=true;
          break;
        }
    }
  if(!vFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no label are visible */
  }

  label_color = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_label_color);
  RepInit(G,&I->R);
  
  obj = cs->Obj;
  I->R.fRender=(void (*)(struct Rep *,RenderInfo *))RepLabelRender;
  I->R.fFree=(void (*)(struct Rep *))RepLabelFree;
  I->R.fRecolor=NULL;
  I->R.obj=(CObject*)obj;
  I->R.cs = cs;
  I->R.context.object = (void*)obj;
  I->R.context.state = state;

  /* raytracing primitives */

  I->L=Alloc(int,cs->NIndex);
  ErrChkPtr(G,I->L);
  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*9);
  ErrChkPtr(G,I->V);


  I->OutlineColor = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_label_outline_color);

  lab_pos = SettingGet_3fv(G,cs->Setting,obj->Obj.Setting,cSetting_label_position);

  if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
    I->R.P=Alloc(Pickable,cs->NIndex+1);
    ErrChkPtr(G,I->R.P);
    rp = I->R.P + 1; /* skip first record! */
  }

  I->N=0;
  
  v=I->V; 
  l=I->L;
  for(a=0;a<cs->NIndex;a++) {
    a1 = cs->IdxToAtm[a];
    ai = obj->AtomInfo+a1;
    if(cs->LabPos) {
      lp = cs->LabPos + a;
    }
    if(ai->visRep[cRepLabel]&&(ai->label)) {
      int at_label_color;
      AtomInfoGetSetting_color(G, ai, cSetting_label_color, 
                               label_color, &at_label_color);

      /*
      float at_label_pos = lab_pos;

        AtomInfoGetSetting_3fv(G, ai, cSetting_label_position, 
        label_pos, &at_label_pos);
      */

      I->N++;
      if((at_label_color>=0) ||
         (at_label_color == cColorFront) ||
         (at_label_color == cColorBack))
        c1 = at_label_color;
      else
        c1=*(cs->Color+a);
      vc = ColorGet(G,c1); /* save new color */
      *(v++)=*(vc++);
      *(v++)=*(vc++);
      *(v++)=*(vc++);
      v0 = cs->Coord+3*a;			 
      *(v++)=*(v0++);
      *(v++)=*(v0++);
      *(v++)=*(v0++);
      if(lp) {
        switch(lp->mode) {
        case 1:  /* local absolute positioning, global relative */
          add3f(lp->offset, v-3, v-3);
          copy3f(lab_pos,v);
          break;
        default:
          copy3f(lab_pos,v);
          break;
        }
      } else {
        copy3f(lab_pos,v);
      }
              
      v+=3;
      if(rp) {
        rp->index = a1;
        rp->bond = cPickableLabel; /* label indicator */
        rp++;
      }
      *(l++) = ai->label;
    }
  }

  if(I->N) {
    I->V = ReallocForSure(I->V,float,(v-I->V));
    I->L = ReallocForSure(I->L,int,(l-I->L));      
    if(rp) {
      I->R.P = ReallocForSure(I->R.P,Pickable,(rp-I->R.P));
      I->R.P[0].index = I->N;  /* unnec? */
    }
  } else {
    I->V=ReallocForSure(I->V,float,1);
    I->L=ReallocForSure(I->L,int,1);
    if(rp) {
      FreeP(I->R.P);
    }
  }
  return((void*)(struct Rep*)I);
}


