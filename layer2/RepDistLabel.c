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


#include<GL/gl.h>
#include<GL/glut.h>
#include<string.h>

#include"OOMac.h"
#include"RepDistLabel.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"Object.h"

typedef char DistLabel[6];

typedef struct RepDistLabel {
  Rep R;
  float *V;
  int N;
  DistLabel *L;
  Object *Obj;
} RepDistLabel;

#include"ObjectDist.h"

void RepDistLabelRender(RepDistLabel *I,CRay *ray,Pickable **pick);
void RepDistLabelFree(RepDistLabel *I);

void RepDistLabelFree(RepDistLabel *I)
{
  FreeP(I->V);
  FreeP(I->L);
  RepFree(&I->R);
  OOFreeP(I);
}

void RepDistLabelRender(RepDistLabel *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  DistLabel *l = I->L;
  char *cc;
  int n = 0;

  if(ray) {
  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 
	 SceneResetNormal(true);
	 while(c--) {
      glRasterPos3fv(v);
		v+=3;
      cc = l[n];
      while(*cc) glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(cc++));
      n++;
	 }
  }
}

Rep *RepDistLabelNew(DistSet *ds)
{
  int a;
  int n;
  float *v,*v1,*v2,d[3],di;
  char buffer[255];

  DistLabel *l;
  OOAlloc(RepDistLabel);

  if(!ds->NIndex) {
    OOFreeP(I);
    return(NULL); 
  }

  RepInit(&I->R);

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepDistLabelRender;
  I->R.fFree=(void (*)(struct Rep *))RepDistLabelFree;
  I->R.fRecolor=NULL;

  I->N=0;
  I->V=NULL;
  I->R.P=NULL;
  I->Obj = (Object*)ds->Obj;

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
      di = diff3f(v1,v2);
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


