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
  char *L;
  int N;
} RepLabel;

#include"ObjectMolecule.h"

void RepLabelRender(RepLabel *I,CRay *ray,Pickable **pick);
void RepLabelFree(RepLabel *I);

void RepLabelInit(void)
{
}

void RepLabelFree(RepLabel *I)
{
  FreeP(I->V);
  FreeP(I->L);
  OOFreeP(I);
}

void RepLabelRender(RepLabel *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  char *l=I->L;
  int font_id = SettingGet_i(I->R.cs->Setting,I->R.obj->Setting,cSetting_label_font_id);
  
  if(ray) {

    if(c) {
      while(c--) {
        if(*l) {
          TextSetPosNColor(v+3,v);
          l = TextRenderRay(ray,font_id,l);
        }
        v+=6;
      }
    }
    
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {

    if(c) {
      int float_text;
      float_text = (int)SettingGet(cSetting_float_labels);
      if(float_text)
        glDisable(GL_DEPTH_TEST);	 
      glDisable(GL_LIGHTING);
      while(c--) {
        if(*l) {
          TextSetPosNColor(v+3,v);
          l = TextRenderOpenGL(font_id,l);
        }
        v+=6;
      }
      glEnable(GL_LIGHTING);
      if(float_text)
        glEnable(GL_DEPTH_TEST);	 
    }
  }
}

Rep *RepLabelNew(CoordSet *cs)
{
  ObjectMolecule *obj;
  int a,a1,vFlag,c1;
  float *v,*v0,*vc;
  char *p,*l;
  int label_color;
  AtomInfoType *ai;
  OOAlloc(RepLabel);
  
  obj = cs->Obj;
  vFlag=false;
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

  label_color = SettingGet_i(cs->Setting,obj->Obj.Setting,cSetting_label_color);

  
  RepInit(&I->R);
  
  obj = cs->Obj;
  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepLabelRender;
  I->R.fFree=(void (*)(struct Rep *))RepLabelFree;
  I->R.fRecolor=NULL;
  I->R.obj=(CObject*)obj;
  I->R.cs = cs;

  /* raytracing primitives */

  I->L=Alloc(char,sizeof(LabelType)*cs->NIndex);
  ErrChkPtr(I->L);
  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*6);
  ErrChkPtr(I->V);

  I->N=0;
  
  v=I->V; 
  l=I->L;
  for(a=0;a<cs->NIndex;a++)
	 {
		a1 = cs->IdxToAtm[a];
      ai = obj->AtomInfo+a1;
		if(ai->visRep[cRepLabel]&&(ai->label[0]))
		  {
			 I->N++;
          if(label_color>=0) 
            c1 = label_color;
          else
            c1=*(cs->Color+a);
			 vc = ColorGet(c1); /* save new color */
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 v0 = cs->Coord+3*a;			 
			 *(v++)=*(v0++);
			 *(v++)=*(v0++);
			 *(v++)=*(v0++);
          p=ai->label;
          while(*p) 
            *(l++)=*(p++);
          *(l++)=0;
		  }
	 }

  if(I->N) 
	 {
		I->V=ReallocForSure(I->V,float,(v-I->V));
		I->L=ReallocForSure(I->L,char,(l-I->L));      
	 }
  else
	 {
		I->V=ReallocForSure(I->V,float,1);
		I->L=ReallocForSure(I->L,char,1);
	 }
  return((void*)(struct Rep*)I);
}


