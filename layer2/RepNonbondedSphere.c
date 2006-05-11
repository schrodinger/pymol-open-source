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
#include"RepNonbondedSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Setting.h"
#include"main.h"

typedef struct RepNonbondedSphere {
  Rep R;
  float *V;
  float *VC;
  SphereRec *SP;
  int N,NC;
  float *VP;
  Pickable *P;
  int NP;

} RepNonbondedSphere;

#include"ObjectMolecule.h"

void RepNonbondedSphereFree(RepNonbondedSphere *I);

void RepNonbondedSphereInit(void)
{
}

void RepNonbondedSphereFree(RepNonbondedSphere *I)
{
  FreeP(I->VP);
  RepPurge(&I->R);
  FreeP(I->VC);
  FreeP(I->V);
  OOFreeP(I);
}

static void RepNonbondedSphereRender(RepNonbondedSphere *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  register PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  int cc=0;
  int a;
  SphereRec *sp;
  int i,j;
  Pickable *p;

  float alpha;

  alpha = SettingGet_f(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_nonbonded_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  
  if(ray) {
    ray->fTransparentf(ray,1.0F-alpha);
	 v=I->VC;
	 c=I->NC;
	 while(c--) {
		ray->fColor3fv(ray,v);
		v+=3;
		ray->fSphere3fv(ray,v,*(v+3));
		v+=4;
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

      sp=I->SP;
      while(c--)
        {
          if(alpha==1.0) {
            glColor3fv(v);
          } else {
            glColor4f(v[0],v[1],v[2],alpha);
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
}

Rep *RepNonbondedSphereNew(CoordSet *cs,int state)
{
  PyMOLGlobals *G=cs->State.G;
  ObjectMolecule *obj;
  int a,c,d,c1;
  float *v,*v0,*vc;
  float nonbonded_size;
  int *q, *s;
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds;
  int *active=NULL;
  AtomInfoType *ai;
  int nSphere = 0;
  int a1;
  float *v1;
  float tmpColor[3];
  OOAlloc(G,RepNonbondedSphere);
  obj = cs->Obj;

  active = Alloc(int,cs->NIndex);
  
  if(obj->RepVisCache[cRepNonbondedSphere])
    for(a=0;a<cs->NIndex;a++) {
      ai = obj->AtomInfo+cs->IdxToAtm[a];
      active[a] =(!ai->bonded) && (ai->visRep[ cRepNonbondedSphere ]);
      if(active[a]) {
        if(ai->masked)
          active[a]=-1;
        else
          active[a]=1;
      }
      if(active[a]) nSphere++;
    }
  if(!nSphere) {
    OOFreeP(I);
    FreeP(active);
    return(NULL); /* skip if no dots are visible */
  }
  
  nonbonded_size = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_nonbonded_size);
  
  /* get current dot sampling */
  ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_dot_density);
  ds=1;
  if(ds<0) ds=0;
  if(ds>3) ds=3;
  sp = G->Sphere->Sphere[ds];

  RepInit(G,&I->R);
  I->R.fRender=(void (*)(struct Rep *, RenderInfo *))RepNonbondedSphereRender;
  I->R.fFree=(void (*)(struct Rep *))RepNonbondedSphereFree;
  I->R.fRecolor=NULL;
  I->R.obj = (CObject*)(cs->Obj);
  I->R.cs = cs;
  I->N=0;
  I->NC=0;
  I->V=NULL;
  I->VC=NULL;
  I->SP=NULL;
  I->NP=0;
  I->VP=NULL;
  I->R.P=NULL;

  /* raytracing primitives */

  I->VC=(float*)mmalloc(sizeof(float)*nSphere*7);
  ErrChkPtr(G,I->VC);
  I->NC=0;

  v=I->VC; 

  for(a=0;a<cs->NIndex;a++)
	 {
      if(active[a])
		  {
			 I->NC++;
			 c1=*(cs->Color+a);
			 v0 = cs->Coord+3*a;			 
          if(ColorCheckRamped(G,c1)) {
            ColorGetRamped(G,c1,v0,tmpColor,state);
            vc = tmpColor;
          } else {
            vc = ColorGet(G,c1);
          }
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(v0++);
			 *(v++)=*(v0++);
			 *(v++)=*(v0++);
			 *(v++)=nonbonded_size;
		  }
	 }

  if(I->NC) 
	 I->VC=ReallocForSure(I->VC,float,(v-I->VC));
  else
	 I->VC=ReallocForSure(I->VC,float,1);

  I->V=(float*)mmalloc(sizeof(float)*nSphere*(3+sp->NVertTot*6));
  ErrChkPtr(G,I->V);

  /* rendering primitives */

  I->N=0;
  I->SP=sp;
  v=I->V;

  for(a=0;a<cs->NIndex;a++)
	 {
		if(active[a])
		  {
			 c1=*(cs->Color+a);
			 v0 = cs->Coord+3*a;
			 vc = ColorGet(G,c1);

          if(ColorCheckRamped(G,c1)) {
            ColorGetRamped(G,c1,v0,tmpColor,state);
            vc = tmpColor;
          } else {
            vc = ColorGet(G,c1);
          }

			 *(v++)=*(vc++);
			 *(v++)=*(vc++);
			 *(v++)=*(vc++);

          q=sp->Sequence;
          s=sp->StripLen;
          
          for(d=0;d<sp->NStrip;d++)
            {
              for(c=0;c<(*s);c++)
                {
                  *(v++)=sp->dot[*q][0]; /* normal */
                  *(v++)=sp->dot[*q][1];
                  *(v++)=sp->dot[*q][2];
                  *(v++)=v0[0]+nonbonded_size*sp->dot[*q][0]; /* point */
                  *(v++)=v0[1]+nonbonded_size*sp->dot[*q][1];
                  *(v++)=v0[2]+nonbonded_size*sp->dot[*q][2];
                  q++;
                }
              s++;
            }
			 I->N++;
		  }
	 }
  
  if(I->N) 
    I->V=ReallocForSure(I->V,float,(v-I->V));
  else
    I->V=ReallocForSure(I->V,float,1);

  /* use pickable representation from nonbonded */
  if(SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_pickable)) {
    I->VP=(float*)mmalloc(sizeof(float)*nSphere*18);
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
    I->VP = Realloc(I->VP,float,I->NP*21);
  }

  FreeP(active);
  return((void*)(struct Rep*)I);
}


