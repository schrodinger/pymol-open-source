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
#include"RepNonbonded.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Setting.h"

typedef struct RepNonbonded {
  Rep R;
  float *V,*VP;
  Pickable *P;
  int N,NP;
  float Width;
  float Radius;
} RepNonbonded;

#include"ObjectMolecule.h"

void RepNonbondedFree(RepNonbonded *I);

void RepNonbondedFree(RepNonbonded *I)
{
  FreeP(I->VP);
  FreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepNonbondedRender(RepNonbonded *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  unsigned int i,j;
  Pickable *p;
  float alpha;

  alpha = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_nonbonded_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  if(ray) {
    float radius;
    ray->fTransparentf(ray,1.0F-alpha);
    
    if(I->Radius==0.0F) {
      radius = ray->PixelRadius*I->Width/2.0F;
    } else {
      radius = I->Radius;
    }

	 v=I->V;
	 c=I->N;
	 
	 while(c--) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]);*/
      ray->fSausage3fv(ray,v+3,v+6,radius,v,v);
      ray->fSausage3fv(ray,v+9,v+12,radius,v,v);
      ray->fSausage3fv(ray,v+15,v+18,radius,v,v);
		v+=21;
	 }
    ray->fTransparentf(ray,0.0);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {

      i=(*pick)->src.index;

      v=I->VP;
      c=I->NP;
      p=I->R.P;
	 
      glBegin(GL_LINES);
	 
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

        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;
        glVertex3fv(v);
        v+=3;

      }
      glEnd();

      (*pick)[0].src.index = i;
    } else {

      int use_dlst;
      
      if(info->width_scale_flag) 
        glLineWidth(I->Width*info->width_scale);
      else
        glLineWidth(I->Width);

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
        if(c) {
          glDisable(GL_LIGHTING);
          glBegin(GL_LINES);	 
          SceneResetNormal(G,true);
          while(c--) {

            if(alpha==1.0) {
              glColor3fv(v);
            } else {
              glColor4f(v[0],v[1],v[2],alpha);
            }
            v+=3;
        
            glVertex3fv(v);
            v+=3;
            glVertex3fv(v);
            v+=3;
        
            glVertex3fv(v);
            v+=3;
            glVertex3fv(v);
            v+=3;
        
            glVertex3fv(v);
            v+=3;
            glVertex3fv(v);
            v+=3;
        
          }
          glEnd();
          glEnable(GL_LIGHTING);
        }
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
      }
    }
  }
}

Rep *RepNonbondedNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  int a,a1,c1;
  float nonbonded_size;
  float *v,*v0,*v1;
  int *active;
  AtomInfoType *ai;
  int nAtom = 0;
  float tmpColor[3];
  OOAlloc(G,RepNonbonded);
  obj = cs->Obj;

  active = Alloc(int,cs->NIndex);
  if(obj->RepVisCache[cRepNonbonded])
    for(a=0;a<cs->NIndex;a++) {
      ai = obj->AtomInfo+cs->IdxToAtm[a];
      active[a] =(!ai->bonded) && (ai->visRep[ cRepNonbonded]);/*&& (!ai->masked);*/
      if(active[a]) {
        if(ai->masked)
          active[a]=-1;
        else
          active[a]=1;
      }
      if(active[a]) nAtom++;
    }
  if(!nAtom) {
    OOFreeP(I);
    FreeP(active);
    return(NULL); /* skip if no dots are visible */
  }

  nonbonded_size = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_nonbonded_size);
  RepInit(G,&I->R);

  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepNonbondedRender;
  I->R.fFree=(void (*)(struct Rep *))RepNonbondedFree;

  I->N=0;
  I->NP=0;
  I->V=NULL;
  I->VP=NULL;
  I->R.P=NULL;
  I->R.fRecolor=NULL;
  I->R.obj = (CObject*)(cs->Obj);
  I->R.cs = cs;

  I->Width = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_line_width);
  I->Radius = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_line_radius);
  I->V=(float*)mmalloc(sizeof(float)*nAtom*21);
  ErrChkPtr(G,I->V);
  v=I->V;
  for(a=0;a<cs->NIndex;a++) 
    if(active[a]) {
      c1=*(cs->Color+a);

      v1 = cs->Coord+3*a;

      if(ColorCheckRamped(G,c1)) {
        ColorGetRamped(G,c1,v1,tmpColor,state);
        v0 = tmpColor;
      } else {
        v0 = ColorGet(G,c1);
      }

      *(v++)=*(v0++);
      *(v++)=*(v0++);
      *(v++)=*(v0++);

      *(v++)=v1[0]-nonbonded_size;
      *(v++)=v1[1];
      *(v++)=v1[2];

      *(v++)=v1[0]+nonbonded_size;
      *(v++)=v1[1];
      *(v++)=v1[2];

      *(v++)=v1[0];
      *(v++)=v1[1]-nonbonded_size;
      *(v++)=v1[2];

      *(v++)=v1[0];
      *(v++)=v1[1]+nonbonded_size;
      *(v++)=v1[2];

      *(v++)=v1[0];
      *(v++)=v1[1];
      *(v++)=v1[2]-nonbonded_size;

      *(v++)=v1[0];
      *(v++)=v1[1];
      *(v++)=v1[2]+nonbonded_size;
      I->N++;
    }
  I->V = ReallocForSure(I->V,float,(v-I->V));

  /* now create pickable verson */
  
  if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
    I->VP=(float*)mmalloc(sizeof(float)*nAtom*21);
    ErrChkPtr(G,I->VP);
    
    I->R.P=Alloc(Pickable,cs->NIndex+1);
    ErrChkPtr(G,I->R.P);
    
    v=I->VP;
    
    for(a=0;a<cs->NIndex;a++) 
      if(active[a]>0) {
        
        a1=cs->IdxToAtm[a];

        if(!obj->AtomInfo[a1].masked) {
          I->NP++;
          I->R.P[I->NP].index = a1;
          I->R.P[I->NP].bond = -1;
          v1 = cs->Coord+3*a;
          
          *(v++)=v1[0]-nonbonded_size;
          *(v++)=v1[1];
          *(v++)=v1[2];
          
          *(v++)=v1[0]+nonbonded_size;
          *(v++)=v1[1];
          *(v++)=v1[2];
          
          *(v++)=v1[0];
          *(v++)=v1[1]-nonbonded_size;
          *(v++)=v1[2];
          
          *(v++)=v1[0];
          *(v++)=v1[1]+nonbonded_size;
          *(v++)=v1[2];
          
          *(v++)=v1[0];
          *(v++)=v1[1];
          *(v++)=v1[2]-nonbonded_size;
          
          *(v++)=v1[0];
          *(v++)=v1[1];
          *(v++)=v1[2]+nonbonded_size;
        }
      }
    I->R.P = Realloc(I->R.P,Pickable,I->NP+1);
    I->R.context.object = (void*)obj;
    I->R.context.state = state;
    I->R.P[0].index = I->NP;
    I->VP = ReallocForSure(I->VP,float,(v-I->VP));
  }
  FreeP(active);
  return((void*)(struct Rep*)I);
}

