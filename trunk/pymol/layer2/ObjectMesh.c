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

static void ObjectMeshFree(ObjectMesh *I) {
  int a;
  for(a=0;a<I->NState;a++) {
    VLAFreeP(I->State[a].N);
    VLAFreeP(I->State[a].V);
  }
  VLAFreeP(I->State);

  OOFreeP(I);
}

void ObjectMeshDump(ObjectMesh *I,char *fname,int state)
{
  float *v;
  int *n;
  int c;
  FILE *f;
  OrthoLineType buf;
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
    sprintf(buf,"%s written to %s\n",I->Obj.Name,fname);
    ErrOk("ObjectMeshDump",buf);
  }
}

static void ObjectMeshUpdate(ObjectMesh *I) 
{
  int a;
  ObjectMeshState *ms;

  for(a=0;a<I->NState;a++) {
    ms = I->State+a;
    
    if(ms->N&&ms->V) {
      if(ms->ResurfaceFlag) {
        ms->ResurfaceFlag=false;
        PRINTF " ObjectMesh: updating \"%s\".\n" , I->Obj.Name ENDF
          if(ms->Map->Field) IsosurfVolume(ms->Map->Field,
                                            ms->Level,
                                            &ms->N,&ms->V,
                                            ms->Range,
                                            ms->DotFlag); 
      }
    }
  }
  SceneDirty();
}

static void ObjectMeshRender(ObjectMesh *I,int state,CRay *ray,Pickable **pick)
{
  float *v = NULL;
  float *vc;
  int *n = NULL;
  int c;
  ObjectMeshState *ms = NULL;

  if(state<I->NState) {
    if(I->State[state].V&&I->State[state].N)
      ms=I->State+state;
  }
  if(!ms) {
    if(I->NState&&I->State[0].V&&I->State[0].N&&SettingGet(cSetting_static_singletons))
      ms=I->State;
  }
  if(ms) {
    v=ms->V;
    n=ms->N;
  }
  if(ray) {
    ms->Radius=SettingGet(cSetting_mesh_radius);
	 if(n&&v) {
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
	 if(n&&v) {
      ObjectUseColor(&I->Obj);
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
  I->State=VLAMalloc(10,sizeof(ObjectMeshState),5,true);

  I->Obj.type = cObjectMesh;

  I->Obj.fFree = (void (*)(struct Object *))ObjectMeshFree;
  I->Obj.fUpdate =  (void (*)(struct Object *)) ObjectMeshUpdate;
  I->Obj.fRender =(void (*)(struct Object *, int, CRay *, Pickable **))ObjectMeshRender;

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
  ms->ResurfaceFlag=true;
  ms->ExtentFlag=false;
}

/*========================================================================*/
ObjectMesh *ObjectMeshFromBox(ObjectMesh *obj,ObjectMap *map,int state,float *mn,float *mx,float level,int dotFlag)
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
  I->Obj.ExtentFlag=extent_flag;
}

