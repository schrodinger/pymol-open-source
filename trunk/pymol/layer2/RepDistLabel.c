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
#include"RepDistLabel.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"
#include"Text.h"

typedef char DistLabel[6];

typedef struct RepDistLabel {
  Rep R;
  float *V;
  int N;
  DistLabel *L;
  CObject *Obj;
} RepDistLabel;

#include"ObjectDist.h"

void RepDistLabelRender(RepDistLabel *I,CRay *ray,Pickable **pick);
void RepDistLabelFree(RepDistLabel *I);

void RepDistLabelFree(RepDistLabel *I)
{
  FreeP(I->V);
  FreeP(I->L);
  RepPurge(&I->R);
  OOFreeP(I);
}

void RepDistLabelRender(RepDistLabel *I,CRay *ray,Pickable **pick)
{
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  DistLabel *l = I->L;
  int n = 0;
  int color;
  int font_id = SettingGet_i(G,NULL,I->Obj->Setting,cSetting_label_font_id);

  if(ray) {

    color = SettingGet_color(G,NULL,I->Obj->Setting,cSetting_label_color);
    
    if(color>=0)
      TextSetColor(G,ColorGet(G,color));
    else
      TextSetColor(G,ColorGet(G,I->Obj->Color));

	 while(c--) {
      TextSetPos(G,v);
      TextRenderRay(G,ray,font_id,l[n]);
      v+=3;
      n++;
	 }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      int float_text;
	 
      float_text = (int)SettingGet(G,cSetting_float_labels);
      if(float_text)
        glDisable(GL_DEPTH_TEST);	 
      glDisable(GL_LIGHTING);
    
      color = SettingGet_color(G,NULL,I->Obj->Setting,cSetting_label_color);
    
      if(color>=0)
        TextSetColor(G,ColorGet(G,color));
      else
        TextSetColor(G,ColorGet(G,I->Obj->Color));
      while(c--) {

        TextSetPos(G,v);
        TextRenderOpenGL(G,font_id,l[n]);
        v+=3;
        n++;
      }
      glEnable(GL_LIGHTING);
      if(float_text)
        glEnable(GL_DEPTH_TEST);	 
    }
  }
}

Rep *RepDistLabelNew(DistSet *ds)
{
  PyMOLGlobals *G=ds->G;
  int a;
  int n;
  float *v,*v1,*v2,d[3],di;
  char buffer[255];

  DistLabel *l;
  OOAlloc(G,RepDistLabel);

  if(!ds->NIndex) {
    OOFreeP(I);
    return(NULL); 
  }

  RepInit(G,&I->R);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepDistLabelRender;
  I->R.fFree=(void (*)(struct Rep *))RepDistLabelFree;
  I->R.fRecolor=NULL;

  I->N=0;
  I->V=NULL;
  I->R.P=NULL;
  I->Obj = (CObject*)ds->Obj;

  n=0;
  if(ds->NIndex) {
	 I->V=Alloc(float,3*(ds->NIndex/2+1));
    I->L=Alloc(DistLabel,(ds->NIndex/2)+1);
    v = I->V;
    l = I->L;
	 for(a=0;a<ds->NIndex;a=a+2) {
      v1 = ds->Coord+3*a;
      v2 = ds->Coord+3*(a+1);
      average3f(v2,v1,d);
      di = (float)diff3f(v1,v2);
      sprintf(buffer,"%1.2f",di);
      buffer[5]=0;
      strcpy(l[n],buffer);
      copy3f(d,v);
      v+=3;
      n++;
    }
    I->N=n;
  }
  return((void*)(struct Rep*)I);
}


