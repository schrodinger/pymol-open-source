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
#include"ObjectSurface.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Debug.h"
#include"Parse.h"
#include"Tetsurf.h"
#include"Vector.h"
#include"Color.h"
#include"main.h"
#include"Scene.h"
#include"Setting.h"

ObjectSurface *ObjectSurfaceNew(void);

static void ObjectSurfaceFree(ObjectSurface *I);
void ObjectSurfaceStateInit(ObjectSurfaceState *ms);
void ObjectSurfaceRecomputeExtent(ObjectSurface *I);

static void ObjectSurfaceStateFree(ObjectSurfaceState *ms)
{
  VLAFreeP(ms->N);
  VLAFreeP(ms->V);
  VLAFreeP(ms->AtomVertex);
  if(ms->UnitCellCGO)
    CGOFree(ms->UnitCellCGO);
}

static void ObjectSurfaceFree(ObjectSurface *I) {
  int a;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active)
      ObjectSurfaceStateFree(I->State+a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  
  OOFreeP(I);
}

void ObjectSurfaceDump(ObjectSurface *I,char *fname,int state)
{
  float *v;
  int *n;
  int c;
  FILE *f;
  f=fopen(fname,"w");
  if(!f) 
    ErrMessage("ObjectSurfaceDump","can't open file for writing");
  else {
    if(state<I->NState) {
      n=I->State[state].N;
      v=I->State[state].V;
      if(n&&v)
        while(*n)
          {
            c=*(n++);
            while(c--) {
              fprintf(f,"%10.4f%10.4f%10.4f\n",v[0],v[1],v[2]);
              v+=3;
            }
          }
    }
    fclose(f);
    PRINTFB(FB_ObjectSurface,FB_Actions)
      " ObjectSurfaceDump: %s written to %s\n",I->Obj.Name,fname
      ENDFB;
  }
}

static void ObjectSurfaceInvalidate(ObjectSurface *I,int rep,int level,int state)
{
  int a;
  for(a=0;a<I->NState;a++) {
    I->State[a].RefreshFlag=true;
  }
}

static void ObjectSurfaceUpdate(ObjectSurface *I) 
{
  int a;
  int c;
  ObjectSurfaceState *ms;
  int *n;
  float *v;
  int *old_n;
  float *old_v;
  int n_cur;
  int n_seg;
  int n_line;
  int flag;
  int last_flag=0;
  int h,k,l;
  int i,j;
  MapType *voxelmap; /* this has nothing to do with isosurfaces... */
  
  for(a=0;a<I->NState;a++) {
    ms = I->State+a;
    if(ms->Active) {
      if(ms->RefreshFlag||ms->ResurfaceFlag) {
        ms->Crystal = *(ms->Map->Crystal);
        if(I->Obj.RepVis[cRepCell]) {
          if(ms->UnitCellCGO)
            CGOFree(ms->UnitCellCGO);
          ms->UnitCellCGO = CrystalGetUnitCellCGO(&ms->Crystal);
        } 
        ms->RefreshFlag=false;
      }

      if(ms->N&&ms->V&&I->Obj.RepVis[cRepSurface]) {
        if(ms->ResurfaceFlag) {
          ms->ResurfaceFlag=false;
          PRINTF " ObjectSurface: updating \"%s\".\n" , I->Obj.Name ENDF;
          if(ms->Map->Field) TetsurfVolume(ms->Map->Field,
                                           ms->Level,
                                           &ms->N,&ms->V,
                                           ms->Range,
                                           ms->Mode); 
          if(ms->CarveFlag&&ms->AtomVertex&&
             VLAGetSize(ms->N)&&VLAGetSize(ms->V)) {
            /* cull my friend, cull */
            voxelmap=MapNew(-ms->CarveBuffer,ms->AtomVertex,VLAGetSize(ms->AtomVertex)/3,NULL);
            if(voxelmap) {
            
              MapSetupExpress(voxelmap);  
            
              old_n = ms->N;
              old_v = ms->V;
              ms->N=VLAlloc(int,VLAGetSize(old_n));
              ms->V=VLAlloc(float,VLAGetSize(old_v));
            
              n = old_n;
              v = old_v;
              n_cur = 0;
              n_seg = 0;
              n_line = 0;
              while(*n)
                {
                  last_flag=false;
                  c=*(n++);
                  while(c--) {
                    flag=false;
                    MapLocus(voxelmap,v,&h,&k,&l);
                    i=*(MapEStart(voxelmap,h,k,l));
                    if(i) {
                      j=voxelmap->EList[i++];
                      while(j>=0) {
                        if(within3f(ms->AtomVertex+3*j,v,ms->CarveBuffer)) {
                          flag=true;
                          break;
                        }
                        j=voxelmap->EList[i++];
                      }
                    }
                    if(flag&&(!last_flag)) {
                      copy3f(v,ms->V+n_line*3);
                      n_cur++;
                      n_line++;
                    }
                    if(flag&&last_flag) { /* continue segment */
                      copy3f(v,ms->V+n_line*3);
                      n_cur++;
                      n_line++;
                    } if((!flag)&&last_flag) { /* terminate segment */
                      ms->N[n_seg]=n_cur;
                      n_seg++;
                      n_cur=0;
                    }
                    last_flag=flag;
                    v+=3;
                  }
                  if(last_flag) { /* terminate segment */
                    ms->N[n_seg]=n_cur;
                    n_seg++;
                    n_cur=0;
                  }
                }
              ms->N[n_seg]=0;
              VLAFreeP(old_n);
              VLAFreeP(old_v);
              MapFree(voxelmap);
            }
          }
        }
      }
    }
  }
  SceneDirty();
}

static void ObjectSurfaceRender(ObjectSurface *I,int state,CRay *ray,Pickable **pick,int pass)
{
  float *v = NULL;
  float *vc;
  int *n = NULL;
  int c;
  int a=0;
  ObjectSurfaceState *ms = NULL;

  if(state<I->NState) {
    if(I->State[state].Active)
      if(I->State[state].V&&I->State[state].N)
        ms=I->State+state;
  }
  while(1) {
    if(state<0) { /* all_states */
      ms = I->State + a;
    } else {
      if(!ms) {
        if(I->NState&&I->State[0].V&&I->State[0].N&&SettingGet(cSetting_static_singletons))
          ms=I->State;
      }
    }
    if(ms) {
      if(ms->Active) {
        v=ms->V;
        n=ms->N;
        if(ray) {
          
          if(ms->UnitCellCGO&&(I->Obj.RepVis[cRepCell]))
            CGORenderRay(ms->UnitCellCGO,ray,ColorGet(I->Obj.Color),
                         I->Obj.Setting,NULL);
          ms->Radius=SettingGet_f(I->Obj.Setting,NULL,cSetting_mesh_radius);
          if(n&&v&&I->Obj.RepVis[cRepSurface]) {
            vc = ColorGet(I->Obj.Color);
            if(ms->Mode) {
              ray->fColor3fv(ray,vc);
              while(*n)
                {
                  c=*(n++);
                  if(c--)
                    {
                      v+=3;
                      while(c--)
                        {
                          ray->fSphere3fv(ray,v,ms->Radius);
                          v+=3;
                        }
                    }
                }
            } else {
              while(*n)
                {
                  c=*(n++);
                  if(c--)
                    {
                      v+=3;
                      while(c--)
                        {
                          ray->fCylinder3fv(ray,v-3,v,ms->Radius,vc,vc);
                          v+=3;
                        }
                    }
                }
            }
          }
        } else if(pick&&PMGUI) {
        } else if(PMGUI) {
          if(!pass) {
            if(ms->UnitCellCGO&&(I->Obj.RepVis[cRepCell]))
              CGORenderGL(ms->UnitCellCGO,ColorGet(I->Obj.Color),
                          I->Obj.Setting,NULL);
            if(n&&v&&I->Obj.RepVis[cRepSurface]) {
              ObjectUseColor(&I->Obj);
              glLineWidth(SettingGet_f(I->Obj.Setting,NULL,cSetting_mesh_width));
              while(*n)
                {
                  c=*(n++);
                  /*                  if(ms->Mode) 
glBegin(GL_LINE_STRIP);
                    else */
                    glBegin(GL_POINTS);



                    
                  SceneResetNormal(false);
                  while(c--) {
                    glVertex3fv(v);
                    v+=3;
                  }
                  glEnd();
                }
            }
          }
        }
      }
    }
    if(state<0) break;
    a = a + 1;
    if(a>=I->NState) break;
  }
}

/*========================================================================*/

static int ObjectSurfaceGetNStates(ObjectSurface *I) 
{
  return(I->NState);
}

/*========================================================================*/
ObjectSurface *ObjectSurfaceNew(void)
{
  OOAlloc(ObjectSurface);
  
  ObjectInit((CObject*)I);
  
  I->NState = 0;
  I->State=VLAMalloc(10,sizeof(ObjectSurfaceState),5,true); /* autozero important */

  I->Obj.type = cObjectSurface;
  
  I->Obj.fFree = (void (*)(struct CObject *))ObjectSurfaceFree;
  I->Obj.fUpdate =  (void (*)(struct CObject *)) ObjectSurfaceUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, int, CRay *, Pickable **,int ))ObjectSurfaceRender;
  I->Obj.fInvalidate =(void (*)(struct CObject *,int,int,int))ObjectSurfaceInvalidate;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectSurfaceGetNStates;
  return(I);
}

/*========================================================================*/
void ObjectSurfaceStateInit(ObjectSurfaceState *ms)
{
  if(!ms->V) {
    ms->V = VLAlloc(float,10000);
  }
  if(!ms->N) {
    ms->N = VLAlloc(int,10000);

  }
  ms->N[0]=0;
  ms->Active=true;
  ms->ResurfaceFlag=true;
  ms->ExtentFlag=false;
  ms->CarveFlag=false;
  ms->AtomVertex=NULL;
  ms->UnitCellCGO=NULL;
}

/*========================================================================*/
ObjectSurface *ObjectSurfaceFromBox(ObjectSurface *obj,ObjectMap *map,
int state,float *mn,float *mx,float level,int mode,
float carve,float *vert_vla)
{
  ObjectSurface *I;
  ObjectSurfaceState *ms;

  if(!obj) {
    I=ObjectSurfaceNew();
  } else {
    I=obj;
  }

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectSurfaceState,state);
    I->NState=state+1;
  }

  ms=I->State+state;
  ObjectSurfaceStateInit(ms);

  ms->Map = map;
  ms->Level = level;
  ms->Mode = mode;
  TetsurfGetRange(ms->Map->Field,ms->Map->Crystal,mn,mx,ms->Range);
  copy3f(mn,ms->ExtentMin); /* this is not exactly correct...should actually take vertex points from range */
  copy3f(mx,ms->ExtentMax);
  ms->ExtentFlag = true;
  if(carve>=0.0) 
  ms->CarveFlag=true;
  ms->CarveBuffer = carve;
  ms->AtomVertex = vert_vla;
  if(I) {
    ObjectSurfaceRecomputeExtent(I);
  }
  I->Obj.ExtentFlag=true;
  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]);*/
  SceneChanged();
  SceneCountFrames();
  return(I);
}

/*========================================================================*/

void ObjectSurfaceRecomputeExtent(ObjectSurface *I)
{
  int extent_flag = false;
  int a;
  ObjectSurfaceState *ms;

  for(a=0;a<I->NState;a++) {
    ms=I->State+a;
    if(ms->Active) {
      if(ms->ExtentFlag) {
        if(!extent_flag) {
          extent_flag=true;
          copy3f(ms->ExtentMax,I->Obj.ExtentMax);
          copy3f(ms->ExtentMin,I->Obj.ExtentMin);
        } else {
          max3f(ms->ExtentMax,I->Obj.ExtentMax,I->Obj.ExtentMax);
          min3f(ms->ExtentMin,I->Obj.ExtentMin,I->Obj.ExtentMin);
        }
      }
    }
  }
  I->Obj.ExtentFlag=extent_flag;
}

