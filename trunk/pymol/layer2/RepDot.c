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
#include"RepDot.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"ObjectMolecule.h"

void RepDotRender(RepDot *I,CRay *ray,Pickable **pick);
void RepDotFree(RepDot *I);

void RepDotInit(void)
{
}

void RepDotFree(RepDot *I)
{
  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->T);
  FreeP(I->F);
  FreeP(I->VN);
  FreeP(I->A);
  FreeP(I->Atom);
  OOFreeP(I);
}

void RepDotRender(RepDot *I,CRay *ray,Pickable **pick)
{
  float *v=I->V;
  int c=I->N;
  int cc=0;

  if(ray) {
	 while(c--)
		{
		  if(!cc) /* load up the current vertex color */
			 {
				cc=(*(v++));
				ray->fColor3fv(ray,v);
				v+=3;
			 }
		  v+=3;
		  ray->fSphere3fv(ray,v,I->dotSize);
		  v+=3;
		  cc--;
		}
	 /*	 v=I->VC;
	 c=I->NC;
	 while(c--) {
		ray->fColor3fv(ray,v);
		v+=3;
		ray->fSphere3fv(ray,v,*(v+3));
		v+=4;
		}*/

  } else if(pick&&PMGUI) {
  } else if(PMGUI) {
	 glBegin(GL_POINTS);
	 while(c--)
		{
		  if(!cc) /* load up the current vertex color */
			 {
				cc=(*(v++));
				glColor3fv(v);
				v+=3;
			 }
		  glNormal3fv(v);
		  v+=3;
		  glVertex3fv(v);
		  v+=3;
		  cc--;
		}
	 glEnd();
  }
}

Rep *RepDotNew(CoordSet *cs)
{
  return(RepDotDoNew(cs,cRepDotNormal));
}

Rep *RepDotDoNew(CoordSet *cs,int mode)
{

  /* this routine does double duty - generating the dot representation,
     but also acting as our surface area computation routine.
     Modes: cRepDotNormal,cRepDotAreaType
  */

  ObjectMolecule *obj;
  int a,b,flag,h,k,l,i,j,c1;
  float *v,*v0,*vc,vdw,*vn;
  float *aa=NULL;
  int *tp=NULL;
  int *tf=NULL;
  float *countPtr = NULL;
  int colorCnt,lastColor;
  Vector3f v1;
  MapType *map;
  SphereRec *sp = Sphere0; 
  int ds;
  float max_vdw = MAX_VDW;
  float solv_rad=0.0;
  int inclH = true;
  int cullByFlag = false;
  int visFlag;
  int atm,*ati=NULL;
  AtomInfoType *ai1,*ai2;

  OOAlloc(RepDot);

  obj = cs->Obj;

  if(mode==cRepDotAreaType) { /* assume all atoms "visible" for area comp. */
    visFlag=true;
  } else {
    visFlag=false;
    for(a=0;a<cs->NIndex;a++) {
      if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepDot])
        {
          visFlag=true;
          break;
        }
    }
  }
  if(!visFlag) {
    OOFreeP(I);
    return(NULL); /* skip if no dots are visible */
  }

  RepInit(&I->R);

  I->dotSize = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_dot_radius);

  I->A=NULL;
  I->T=NULL;
  I->F=NULL;
  I->V=NULL;
  I->VC=NULL;
  I->VN=NULL;
  I->Atom=NULL;
  I->R.fRecolor=NULL;

  cullByFlag = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_trim_dots); /* are we using flags 24 & 25 */
  inclH = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_dot_hydrogens); /* are we ignoring hydrogens? */
  if(SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_dot_mode)>0.0) { /* are we generating a solvent surface? */
    solv_rad = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_solvent_radius); /* if so, get solvent radius */
  }

  /* get current dot sampling */
  ds = (int)SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_dot_density);

  max_vdw+=solv_rad;

  if(ds<0) ds=0;
  switch(ds) { /* Note: significantly affects the accuracy of our area comp. */
  case 0: sp=Sphere0; break;
  case 1: sp=Sphere1; break;
  case 2: sp=Sphere2; break;
  default: sp=Sphere3; break;
  }

  I->R.fRender=(void (*)(struct Rep *, CRay *, Pickable **))RepDotRender;
  I->R.fFree=(void (*)(struct Rep *))RepDotFree;
  
  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*sp->nDot*10);
  ErrChkPtr(I->V);

  if(mode==cRepDotAreaType) { /* in area mode, we need to export save addl. info 
                               * such as the normal vectors, the partial area, 
                               * the originating atom, etc. */
	 I->A=Alloc(float,cs->NIndex*sp->nDot);
	 I->T=Alloc(int,cs->NIndex*sp->nDot);
	 I->F=Alloc(int,cs->NIndex*sp->nDot);
	 I->VN=Alloc(float,cs->NIndex*sp->nDot*3);
    I->Atom=Alloc(int,cs->NIndex*sp->nDot);
	 aa=I->A;
	 tp=I->T;
	 tf=I->F;
    ati=I->Atom;
	 inclH=true;
	 cullByFlag=true;
  }
  vn=I->VN;

  I->N=0;
  lastColor=-1;
  colorCnt=0;
  map=MapNew(max_vdw,cs->Coord,cs->NIndex,NULL);
  v=I->V;
  if(map)
	 {
		MapSetupExpress(map);
		for(a=0;a<cs->NIndex;a++)
		  {
          atm = cs->IdxToAtm[a];
          ai1 = obj->AtomInfo+atm;
			 if(ai1->visRep[cRepDot]||mode==cRepDotAreaType) 
				if((inclH||(!ai1->hydrogen))&&
               ((!cullByFlag)||
                (!(ai1->flags&cAtomFlag_exclude)))) {
              /* If we are culling, flag 24 controls which atoms 
                 will have dot surfaces generated for them.
              */
              if(cs->Color)
                c1=*(cs->Color+a);
              else
                c1 = 0;
				  v0 = cs->Coord+3*a;
				  vdw = ai1->vdw+solv_rad;
				  for(b=0;b<sp->nDot;b++)
					 {
						v1[0]=v0[0]+vdw*sp->dot[b].v[0];
						v1[1]=v0[1]+vdw*sp->dot[b].v[1];
						v1[2]=v0[2]+vdw*sp->dot[b].v[2];
						
						MapLocus(map,v1,&h,&k,&l);

						flag=true;

						i=*(MapEStart(map,h,k,l));
						if(i) {
						  j=map->EList[i++];
						  while(j>=0) {
                      ai2 = obj->AtomInfo+cs->IdxToAtm[j];
							 if((inclH||(!(ai2->hydrogen)))&&
								 ((!cullByFlag)||
                          (!(ai2->flags&cAtomFlag_ignore))))  
                        /* If we are cullilng, flag 25 controls which atoms 
                           are considered "present" in the surface area 
                           calculation (i.e. able to occlude surface) */
								if(j!=a)
								  if(within3f(cs->Coord+3*j,v1,ai2->vdw+solv_rad))
									 {
										flag=false;
										break;
									 }
							 j=map->EList[i++];
						  }
						}
						if(flag)
						  {
							 switch(mode) {
							 case cRepDotNormal:
								if(lastColor!=c1) /* new color */
								  {
									 if(countPtr) /* after first pass */
										*countPtr=(float)colorCnt; /* save count */
									 colorCnt=1;
									 countPtr=v++;
									 vc = ColorGet(c1); /* save new color */
									 lastColor=c1;
									 *(v++)=*(vc++);
									 *(v++)=*(vc++);
									 *(v++)=*(vc++);
								  }
								else 
								  colorCnt++;
								*(v++)=sp->dot[b].v[0];
								*(v++)=sp->dot[b].v[1];
								*(v++)=sp->dot[b].v[2];
								*(v++)=v1[0];
								*(v++)=v1[1];
								*(v++)=v1[2];
								I->N++;
								break;
							 case cRepDotAreaType:
								*(v++)=v1[0];
								*(v++)=v1[1];
								*(v++)=v1[2];
								*(aa++)=vdw*vdw*sp->dot[b].area; /* area */
								*(tp++)=ai1->customType; /* numeric type */
								*(tf++)=ai1->flags; /* flags */
								*(vn++)=sp->dot[b].v[0];
								*(vn++)=sp->dot[b].v[1];
								*(vn++)=sp->dot[b].v[2];
                        *(ati++)=atm;
								I->N++;
								break;
							 }
						  }
					 }
				}
		  }
		if(countPtr) *countPtr=(float)colorCnt; /* save count */
		MapFree(map);
	 }
  
  I->V = Realloc(I->V,float,(v-I->V));
  
  if(mode==cRepDotAreaType) {
	 I->A = Realloc(I->A,float,(aa-I->A));
	 I->T= Realloc(I->T,int,(tp-I->T));
	 I->F= Realloc(I->F,int,(tf-I->F));
	 I->VN= Realloc(I->VN,float,(vn-I->VN));
    I->Atom= Realloc(I->Atom,int,(ati-I->Atom));
  }
  return((void*)(struct Rep*)I);
}




