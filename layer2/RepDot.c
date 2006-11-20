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
#include"Scene.h"

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

static void RepDotRender(RepDot *I,RenderInfo *info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G=I->R.G;
  float *v=I->V;
  int c=I->N;
  int cc=0;
  if(ray) {
    float radius;

    if(I->dotSize<=0.0F) {
      radius = ray->PixelRadius*I->Width/1.4142F;
    } else {
      radius = I->dotSize;
    }

	 while(c--)
		{
		  if(!cc) {/* load up the current vertex color */
            cc=(int)(*(v++));
            ray->fColor3fv(ray,v);
				v+=3;
          }
		  v+=3;
		  ray->fSphere3fv(ray,v,radius);
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

  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      int normals = SettingGet_i(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_dot_normals);
      int lighting = SettingGet_i(G,I->R.cs->Setting,I->R.obj->Setting,cSetting_dot_lighting);
      int use_dlst;

      if(!normals)
        SceneResetNormal(G,true);
      if(!lighting)
        glDisable(GL_LIGHTING);
      
      
      use_dlst = (int)SettingGet(G,cSetting_use_display_lists);

      if(info->width_scale_flag) 
        glPointSize(I->Width*info->width_scale);
      else
        glPointSize(I->Width);

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


        glBegin(GL_POINTS);
        while(c--)
          {
            if(!cc) /* load up the current vertex color */
              {
                cc=(int)(*(v++));
                glColor3fv(v);
                v+=3;
              }
            if(normals) 
              glNormal3fv(v);
            v+=3;
            glVertex3fv(v);
            v+=3;
            cc--;
          }
        glEnd();
      
        if(use_dlst&&I->R.displayList) {
          glEndList();
        }
        if(!lighting)
          glEnable(GL_LIGHTING);
      }
    }
  }
}

Rep *RepDotNew(CoordSet *cs,int state)
{
  return(RepDotDoNew(cs,cRepDotNormal,state));
}

Rep *RepDotDoNew(CoordSet *cs,int mode,int state)
{

  /* this routine does double duty - generating the dot representation,
     but also acting as our surface area computation routine.
     Modes: cRepDotNormal,cRepDotAreaType
  */
  PyMOLGlobals *G=cs->State.G;
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
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds;
  float max_vdw = MAX_VDW;
  float solv_rad=0.0;
  int inclH = true;
  int cullByFlag = false;
  int visFlag;
  int atm,*ati=NULL;
  AtomInfoType *ai1,*ai2;
  int dot_color;
  OOAlloc(G,RepDot);

  obj = cs->Obj;

  if(mode==cRepDotAreaType) { /* assume all atoms "visible" for area comp. */
    visFlag=true;
  } else {
    visFlag=false;
    if(obj->RepVisCache[cRepDot])
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

  RepInit(G,&I->R);

  I->dotSize = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_dot_radius);

  I->A=NULL;
  I->T=NULL;
  I->F=NULL;
  I->V=NULL;
  I->VC=NULL;
  I->VN=NULL;
  I->Atom=NULL;
  I->R.fRecolor=NULL;

  I->Width = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_dot_width);
  cullByFlag = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_trim_dots); /* are we using flags 24 & 25 */

  dot_color = SettingGet_color(G,cs->Setting,obj->Obj.Setting,cSetting_dot_color); /* are we using flags 24 & 25 */
  inclH = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_dot_hydrogens); /* are we ignoring hydrogens? */
  if(SettingGet_b(G,cs->Setting,obj->Obj.Setting,cSetting_dot_solvent)) { /* are we generating a solvent surface? */
    solv_rad = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_solvent_radius); /* if so, get solvent radius */
  }

  /* get current dot sampling */
  ds = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_dot_density);

  max_vdw+=solv_rad;

/* Note: significantly affects the accuracy of our area comp. */
  if(ds<0) ds=0;
  if(ds>4) ds=4;
  sp = G->Sphere->Sphere[ds];

  I->R.fRender=(void (*)(struct Rep *, RenderInfo *info))RepDotRender;
  I->R.fFree=(void (*)(struct Rep *))RepDotFree;
  I->R.obj=(CObject*)obj;
  I->R.cs = cs;

  I->V=(float*)mmalloc(sizeof(float)*cs->NIndex*sp->nDot*10);
  ErrChkPtr(G,I->V);

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
  map=MapNew(G,max_vdw,cs->Coord,cs->NIndex,NULL);
  v=I->V;
  if(map) {
    MapSetupExpress(map);
    for(a=0;a<cs->NIndex;a++) {
      atm = cs->IdxToAtm[a];
      ai1 = obj->AtomInfo+atm;
      if(ai1->visRep[cRepDot]||mode==cRepDotAreaType) 
        if((inclH||(!ai1->hydrogen))&&
           ((!cullByFlag)||
            (!(ai1->flags&cAtomFlag_exfoliate)))) {
          int at_dot_color;
          
          AtomInfoGetSetting_color(G, ai1, cSetting_dot_color, 
                                   dot_color, &at_dot_color);
          

          /* If we are culling, flag 24 controls which atoms 
             will have dot surfaces generated for them.
          */
          if(at_dot_color==-1) {
            if(cs->Color)
              c1=*(cs->Color+a);
            else
              c1 = 0;
          } else {
            c1 = at_dot_color;
          }
          v0 = cs->Coord+3*a;
          vdw = ai1->vdw+solv_rad;
          for(b=0;b<sp->nDot;b++) {
            v1[0]=v0[0]+vdw*sp->dot[b][0];
            v1[1]=v0[1]+vdw*sp->dot[b][1];
            v1[2]=v0[2]+vdw*sp->dot[b][2];
            
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
                    if(within3f(cs->Coord+3*j,v1,ai2->vdw+solv_rad)) {
                      flag=false;
                      break;
                    }
                j=map->EList[i++];
              }
            }
            if(flag) {
              switch(mode) {
              case cRepDotNormal:
                
                if((lastColor!=c1)||ColorCheckRamped(G,c1)) { /* new color */
                  if(countPtr) /* after first pass */
                    *countPtr=(float)colorCnt; /* save count */
                  colorCnt=1;
                  countPtr=v++;
                  vc = ColorGet(G,c1); /* save new color */
                  lastColor=c1;
                  if(ColorCheckRamped(G,c1)) {
                    ColorGetRamped(G,c1,v1,v,state);
                    v+=3;
                  } else {
                    *(v++)=*(vc++);
                    *(v++)=*(vc++);
                    *(v++)=*(vc++);
                  }
                }
                else 
                  colorCnt++;
                *(v++)=sp->dot[b][0];
                *(v++)=sp->dot[b][1];
                *(v++)=sp->dot[b][2];
                *(v++)=v1[0];
                *(v++)=v1[1];
                *(v++)=v1[2];
                I->N++;
                break;
              case cRepDotAreaType:
                *(v++)=v1[0];
                *(v++)=v1[1];
                *(v++)=v1[2];
                *(aa++)=vdw*vdw*sp->area[b]; /* area */
                *(tp++)=ai1->customType; /* numeric type */
                *(tf++)=ai1->flags; /* flags */
                *(vn++)=sp->dot[b][0];
                *(vn++)=sp->dot[b][1];
                *(vn++)=sp->dot[b][2];
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
  
  I->V = ReallocForSure(I->V,float,(v-I->V));
  
  if(mode==cRepDotAreaType) {
    I->A = ReallocForSure(I->A,float,(aa-I->A));
    I->T= ReallocForSure(I->T,int,(tp-I->T));
    I->F= ReallocForSure(I->F,int,(tf-I->F));
    I->VN= ReallocForSure(I->VN,float,(vn-I->VN));
    I->Atom= ReallocForSure(I->Atom,int,(ati-I->Atom));
  }
  return((void*)(struct Rep*)I);
}




