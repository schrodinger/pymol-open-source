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

#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectMesh.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Debug.h"
#include"Parse.h"
#include"Isosurf.h"
#include"Vector.h"
#include"Color.h"
#include"main.h"
#include"Scene.h"
#include"Setting.h"

ObjectMesh *ObjectMeshNew(void);

static void ObjectMeshFree(ObjectMesh *I);
void ObjectMeshStateInit(ObjectMeshState *ms);
void ObjectMeshRecomputeExtent(ObjectMesh *I);

static void ObjectMeshStateFree(ObjectMeshState *ms)
{
  VLAFreeP(ms->N);
  VLAFreeP(ms->V);
  VLAFreeP(ms->AtomVertex);
  if(ms->UnitCellCGO)
    CGOFree(ms->UnitCellCGO);
}

static void ObjectMeshFree(ObjectMesh *I) {
  int a;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active)
      ObjectMeshStateFree(I->State+a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  
  OOFreeP(I);
}

void ObjectMeshDump(ObjectMesh *I,char *fname,int state)
{
  float *v;
  int *n;
  int c;
  FILE *f;
  f=fopen(fname,"w");
  if(!f) 
    ErrMessage("ObjectMeshDump","can't open file for writing");
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
    PRINTFB(FB_ObjectMesh,FB_Actions)
      " ObjectMeshDump: %s written to %s\n",I->Obj.Name,fname
      ENDFB;
  }
}

static void ObjectMeshInvalidate(ObjectMesh *I,int rep,int level,int state)
{
  int a;
  for(a=0;a<I->NState;a++) {
    I->State[a].RefreshFlag=true;
  }
}

static void ObjectMeshUpdate(ObjectMesh *I) 
{
  int a;
  int c;
  ObjectMeshState *ms;
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

      if(ms->N&&ms->V&&I->Obj.RepVis[cRepMesh]) {
        if(ms->ResurfaceFlag) {
          ms->ResurfaceFlag=false;
          PRINTF " ObjectMesh: updating \"%s\".\n" , I->Obj.Name ENDF;
          if(ms->Map->Field) IsosurfVolume(ms->Map->Field,
                                           ms->Level,
                                           &ms->N,&ms->V,
                                           ms->Range,
                                           ms->DotFlag); 
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

static void ObjectMeshRender(ObjectMesh *I,int state,CRay *ray,Pickable **pick,int pass)
{
  float *v = NULL;
  float *vc;
  int *n = NULL;
  int c;
  int a=0;
  ObjectMeshState *ms = NULL;

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
          if(n&&v&&I->Obj.RepVis[cRepMesh]) {
            vc = ColorGet(I->Obj.Color);
            if(ms->DotFlag) {
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
            if(n&&v&&I->Obj.RepVis[cRepMesh]) {
              ObjectUseColor(&I->Obj);
              glLineWidth(SettingGet_f(I->Obj.Setting,NULL,cSetting_mesh_width));
              while(*n)
                {
                  c=*(n++);
                  if(ms->DotFlag) 
                    glBegin(GL_POINTS);
                  else 
                    glBegin(GL_LINE_STRIP);
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

static int ObjectMeshGetNStates(ObjectMesh *I) 
{
  return(I->NState);
}

/*========================================================================*/
ObjectMesh *ObjectMeshNew(void)
{
  OOAlloc(ObjectMesh);
  
  ObjectInit((Object*)I);
  
  I->NState = 0;
  I->State=VLAMalloc(10,sizeof(ObjectMeshState),5,true); /* autozero important */

  I->Obj.type = cObjectMesh;
  
  I->Obj.fFree = (void (*)(struct Object *))ObjectMeshFree;
  I->Obj.fUpdate =  (void (*)(struct Object *)) ObjectMeshUpdate;
  I->Obj.fRender =(void (*)(struct Object *, int, CRay *, Pickable **,int ))ObjectMeshRender;
  I->Obj.fInvalidate =(void (*)(struct Object *,int,int,int))ObjectMeshInvalidate;
  I->Obj.fGetNFrame = (int (*)(struct Object *)) ObjectMeshGetNStates;
  return(I);
}

/*========================================================================*/
void ObjectMeshStateInit(ObjectMeshState *ms)
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
ObjectMesh *ObjectMeshFromBox(ObjectMesh *obj,ObjectMap *map,
int state,float *mn,float *mx,float level,int dotFlag,
float carve,float *vert_vla)
{
  ObjectMesh *I;
  ObjectMeshState *ms;

  if(!obj) {
    I=ObjectMeshNew();
  } else {
    I=obj;
  }

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectMeshState,state);
    I->NState=state+1;
  }

  ms=I->State+state;
  ObjectMeshStateInit(ms);

  ms->Map = map;
  ms->Level = level;
  ms->DotFlag = dotFlag;
  IsosurfGetRange(ms->Map->Field,ms->Map->Crystal,mn,mx,ms->Range);
  copy3f(mn,ms->ExtentMin); /* this is not exactly correct...should actually take vertex points from range */
  copy3f(mx,ms->ExtentMax);
  ms->ExtentFlag = true;
  if(carve>=0.0) 
  ms->CarveFlag=true;
  ms->CarveBuffer = carve;
  ms->AtomVertex = vert_vla;
  if(I) {
    ObjectMeshRecomputeExtent(I);
  }
  I->Obj.ExtentFlag=true;
  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]);*/
  SceneChanged();
  SceneCountFrames();
  return(I);
}

/*========================================================================*/

void ObjectMeshRecomputeExtent(ObjectMesh *I)
{
  int extent_flag = false;
  int a;
  ObjectMeshState *ms;

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

