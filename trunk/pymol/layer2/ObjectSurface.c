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
#include"Executive.h"
#include"PConv.h"

ObjectSurface *ObjectSurfaceNew(void);

static void ObjectSurfaceFree(ObjectSurface *I);
void ObjectSurfaceStateInit(ObjectSurfaceState *ms);
void ObjectSurfaceRecomputeExtent(ObjectSurface *I);


static PyObject *ObjectSurfaceStateAsPyList(ObjectSurfaceState *I)
{
  PyObject *result = NULL;

  result = PyList_New(15);
  
  PyList_SetItem(result,0,PyInt_FromLong(I->Active));
  PyList_SetItem(result,1,PyString_FromString(I->MapName));
  PyList_SetItem(result,2,PyInt_FromLong(I->MapState));
  PyList_SetItem(result,3,CrystalAsPyList(&I->Crystal));
  PyList_SetItem(result,4,PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result,5,PConvFloatArrayToPyList(I->ExtentMin,3));
  PyList_SetItem(result,6,PConvFloatArrayToPyList(I->ExtentMax,3));
  PyList_SetItem(result,7,PConvIntArrayToPyList(I->Range,6));
  PyList_SetItem(result,8,PyFloat_FromDouble(I->Level));
  PyList_SetItem(result,9,PyFloat_FromDouble(I->Radius));
  PyList_SetItem(result,10,PyInt_FromLong(I->CarveFlag));
  PyList_SetItem(result,11,PyFloat_FromDouble(I->CarveBuffer));
  if(I->CarveFlag&&I->AtomVertex) {
    PyList_SetItem(result,12,PConvFloatVLAToPyList(I->AtomVertex));
  } else {
    PyList_SetItem(result,12,PConvAutoNone(NULL));
  }
  PyList_SetItem(result,13,PyInt_FromLong(I->DotFlag));
  PyList_SetItem(result,14,PyInt_FromLong(I->Mode));
  PyList_SetItem(result,15,PyInt_FromLong(I->Side));

#if 0
  char MapName[ObjNameMax];
  int MapState;
  CCrystal Crystal;
  int Active;
  int *N;
  float *V;
  int Range[6];
  float ExtentMin[3],ExtentMax[3];
  int ExtentFlag;
  float Level,Radius;
  int RefreshFlag;
  int ResurfaceFlag;
  float *AtomVertex;
  int CarveFlag;
  float CarveBuffer;
  int DotFlag;
  CGO *UnitCellCGO;
#endif

  return(PConvAutoNone(result));  
}


static PyObject *ObjectSurfaceAllStatesAsPyList(ObjectSurface *I)
{
  
  PyObject *result=NULL;
  int a;
  result = PyList_New(I->NState);
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result,a,ObjectSurfaceStateAsPyList(I->State+a));
    } else {
      PyList_SetItem(result,a,PConvAutoNone(NULL));
    }
  }
  return(PConvAutoNone(result));  

}

static int ObjectSurfaceStateFromPyList(ObjectSurfaceState *I,PyObject *list)
{
  int ok=true;
  int ll=0;
  PyObject *tmp;
  if(ok) ok=(list!=NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active=false;
    else {
      ObjectSurfaceStateInit(I);
      if(ok) ok=PyList_Check(list);
      if(ok) ll=PyList_Size(list);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->Active);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,1),I->MapName,ObjNameMax);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->MapState);
      if(ok) ok = CrystalFromPyList(&I->Crystal,PyList_GetItem(list,3));
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,4),&I->ExtentFlag);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,5),I->ExtentMin,3);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,6),I->ExtentMax,3);
      if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,7),I->Range,6);
      if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,8),&I->Level);
      if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,9),&I->Radius);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,10),&I->CarveFlag);
      if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,11),&I->CarveBuffer);
      if(ok) {
        tmp = PyList_GetItem(list,12);
        if(tmp == Py_None)
          I->AtomVertex = NULL;
        else 
          ok = PConvPyListToFloatVLA(tmp,&I->AtomVertex);
      }
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,13),&I->DotFlag);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,14),&I->Mode);

      if(ok&&(ll>15)) PConvPyIntToInt(PyList_GetItem(list,15),&I->Side);

      if(ok) {
        I->RefreshFlag=true;
        I->ResurfaceFlag=true;
      }
    }
  }
  return(ok);
}

static int ObjectSurfaceAllStatesFromPyList(ObjectSurface *I,PyObject *list)
{
  int ok=true;
  int a;
  VLACheck(I->State,ObjectSurfaceState,I->NState);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    for(a=0;a<I->NState;a++) {
      ok = ObjectSurfaceStateFromPyList(I->State+a,PyList_GetItem(list,a));
      if(!ok) break;
    }
  }
  return(ok);
}

int ObjectSurfaceNewFromPyList(PyObject *list,ObjectSurface **result)
{
  int ok = true;
  ObjectSurface *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);

  I=ObjectSurfaceNew();
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NState);
  if(ok) ok = ObjectSurfaceAllStatesFromPyList(I,PyList_GetItem(list,2));
  if(ok) {
    (*result) = I;
    ObjectSurfaceRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return(ok);
}

PyObject *ObjectSurfaceAsPyList(ObjectSurface *I)
{
  
  PyObject *result=NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NState));
  PyList_SetItem(result,2,ObjectSurfaceAllStatesAsPyList(I));

  return(PConvAutoNone(result));  
}

static void ObjectSurfaceStateFree(ObjectSurfaceState *ms)
{
  if(PMGUI) {
    if(ms->displayList) {
      glDeleteLists(ms->displayList,1);
      ms->displayList = 0;
    }
  }
  
  
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
  f=fopen(fname,"wb");
  if(!f) 
    ErrMessage("ObjectSurfaceDump","can't open file for writing");
  else {
    if(state<I->NState) {
      n=I->State[state].N;
      v=I->State[state].V;
      if(n&&v)
        while(*n)
          {
            v+=12;
            c=*(n++);
            c-=4;
            while(c>0) {
              fprintf(f,
"%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
                      *(v-9),*(v-8),*(v-7),
                      *(v-12),*(v-11),*(v-10),
                      *(v-3),*(v-2),*(v-1),
                      *(v-6),*(v-5),*(v-4),
                      *(v+3),*(v+4),*(v+5),
                      *(v),  *(v+1),*(v+2));

              v+=6;
              c-=2;
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
  ObjectSurfaceState *ms;
  ObjectMapState *oms = NULL;
  ObjectMap *map = NULL;
  MapType *voxelmap=NULL; /* this has nothing to do with isosurfaces... */
  int ok=true;

  for(a=0;a<I->NState;a++) {
    ms = I->State+a;
    if(ms->Active) {
      map = ExecutiveFindObjectMapByName(ms->MapName);
      if(!map) {
        ok=false;
        PRINTFB(FB_ObjectSurface,FB_Errors)
          "ObjectSurfaceUpdate-Error: map '%s' has been deleted.\n",ms->MapName
          ENDFB;
        ms->ResurfaceFlag=false;
      }
      if(map) {
        oms = ObjectMapGetState(map,ms->MapState);
        if(!oms) ok=false;
      }
      if(oms) {
        if(ms->RefreshFlag||ms->ResurfaceFlag) {
          ms->Crystal = *(oms->Crystal);
          if(I->Obj.RepVis[cRepCell]) {
            if(ms->UnitCellCGO)
              CGOFree(ms->UnitCellCGO);
            ms->UnitCellCGO = CrystalGetUnitCellCGO(&ms->Crystal);
          } 
          ms->RefreshFlag=false;
        }
      }
      if(map&&ms&&oms&&ms->N&&ms->V&&I->Obj.RepVis[cRepSurface]) {
        if(ms->ResurfaceFlag) {


          ms->ResurfaceFlag=false;
          PRINTF " ObjectSurface: updating \"%s\".\n" , I->Obj.Name ENDF;
          if(oms->Field) {
            TetsurfGetRange(oms->Field,oms->Crystal,
                            ms->ExtentMin,ms->ExtentMax,ms->Range);

            if(ms->CarveFlag&&ms->AtomVertex) {
              voxelmap=MapNew(-ms->CarveBuffer,ms->AtomVertex,
                              VLAGetSize(ms->AtomVertex)/3,NULL);
              if(voxelmap)
                MapSetupExpress(voxelmap);  
            }
            TetsurfVolume(oms->Field,
                          ms->Level,
                          &ms->N,&ms->V,
                          ms->Range,
                          ms->Mode,
                          voxelmap,
                          ms->AtomVertex,
                          ms->CarveBuffer,
                          ms->Side); 
            if(voxelmap)
              MapFree(voxelmap);
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
  float *col;
  int *n = NULL;
  int c;
  int a=0;
  ObjectSurfaceState *ms = NULL;
  float alpha;

  ObjectPrepareContext(&I->Obj,ray);

  alpha = SettingGet_f(NULL,I->Obj.Setting,cSetting_transparency);
  alpha=1.0-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0;
  
  if(state>=0) 
    if(state<I->NState) 
      if(I->State[state].Active)
        if(I->State[state].V&&I->State[state].N)
          ms=I->State+state;
  while(1) {
    if(state<0) { /* all_states */
      ms = I->State + a;
    } else {
      if(!ms) {
        if(I->NState&&
           ((SettingGet(cSetting_static_singletons)&&(I->NState==1))))
          ms=I->State;
      }
    }
    if(ms) {
      if(ms->Active&&ms->V&&ms->N) {
        v=ms->V;
        n=ms->N;
        if(ray) {
          ray->fTransparentf(ray,1.0-alpha);       
          if(ms->UnitCellCGO&&(I->Obj.RepVis[cRepCell]))
            CGORenderRay(ms->UnitCellCGO,ray,ColorGet(I->Obj.Color),
                         I->Obj.Setting,NULL);
          ms->Radius=SettingGet_f(I->Obj.Setting,NULL,cSetting_mesh_radius);
          if(n&&v&&I->Obj.RepVis[cRepSurface]) {
            vc = ColorGet(I->Obj.Color);
            while(*n)
              {
                c=*(n++);
                switch(ms->Mode) {
                case 2:
                  v+=12;
                  c-=4;
                  while(c>0) {
                    /*                  dump3f(v-12,"v-12 ");
                                        dump3f(v-9 ,"v-9  ");
                                        dump3f(v-6, "v-6  ");
                                        dump3f(v-3, "v-3  ");
                                        dump3f(v,   "v    ");
                                        dump3f(v+3, "v+3  ");
                    */
                    
                    ray->fTriangle3fv(ray,v-9,v-3,v+3,
                                      v-12,v-6,v,
                                      vc,vc,vc);
                    v+=6;
                    c-=2;
                  }
                  break;
                case 1:
                  c--;
                  v+=3;
                  while(c>0) {
                    ray->fSausage3fv(ray,v-3,v,ms->Radius,vc,vc);
                    v+=3;
                    c--;
                  }
                  break;
                case 0:
                default:
                  while(c>0) {
                    ray->fSphere3fv(ray,v,ms->Radius);
                    v+=3;
                    c--;
                  }
                  break;
                }
              }
            
          }
          ray->fTransparentf(ray,0.0);
        } else if(pick&&PMGUI) {
        } else if(PMGUI) {
          if(!pass) {

            int use_dlst;

            if(ms->UnitCellCGO&&(I->Obj.RepVis[cRepCell]))
              CGORenderGL(ms->UnitCellCGO,ColorGet(I->Obj.Color),
                          I->Obj.Setting,NULL);

            SceneResetNormal(false);
            col = ColorGet(I->Obj.Color);
            glColor4f(col[0],col[1],col[2],alpha);

            use_dlst = (int)SettingGet(cSetting_use_display_lists);
            if(use_dlst&&ms->displayList) {
              glCallList(ms->displayList);
            } else { 
              
              if(use_dlst) {
                if(!ms->displayList) {
                  ms->displayList = glGenLists(1);
                  if(ms->displayList) {
                    glNewList(ms->displayList,GL_COMPILE_AND_EXECUTE);
                  }
                }
              }

            
            if(n&&v&&I->Obj.RepVis[cRepSurface]) {

              glLineWidth(SettingGet_f(I->Obj.Setting,NULL,cSetting_mesh_width));
              while(*n)
                {
                  c=*(n++);
                  switch(ms->Mode) {
                  case 2:
                    glBegin(GL_TRIANGLE_STRIP);
                    while(c>0) {
                      glNormal3fv(v);
                      v+=3;
                      glVertex3fv(v);
                      v+=3;
                      c-=2;
                    }
                    glEnd();
                    break;
                  case 1:
                    glBegin(GL_LINES);
                    while(c>0) {
                      glVertex3fv(v);
                      v+=3;
                      c--;
                    }
                    glEnd();
                    break;
                  case 0:
                  default:
                    glBegin(GL_POINTS);
                    while(c>0) {
                      glVertex3fv(v);
                      v+=3;
                      c--;
                    }
                    glEnd();
                    break;
                  }
                }
            }
            }

            if(use_dlst&&ms->displayList) {
              glEndList();
            }
          }
        }
      }
    }
    if(state>=0) break; /* only rendering one state */
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
  if(ms->AtomVertex) {
    VLAFreeP(ms->AtomVertex);
  }

  ms->N[0]=0;
  ms->Active=true;
  ms->ResurfaceFlag=true;
  ms->ExtentFlag=false;
  ms->CarveFlag=false;
  ms->AtomVertex=NULL;
  ms->UnitCellCGO=NULL;
  ms->Side = 0;
  ms->displayList = 0;
}

/*========================================================================*/
ObjectSurface *ObjectSurfaceFromBox(ObjectSurface *obj,ObjectMap *map,
                                    int map_state,
int state,float *mn,float *mx,float level,int mode,
float carve,float *vert_vla,int side)
{
  ObjectSurface *I;
  ObjectSurfaceState *ms;
  ObjectMapState *oms;

  mode=2;
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

  strcpy(ms->MapName,map->Obj.Name);
  ms->MapState = map_state;
  oms = ObjectMapGetState(map,map_state);

  ms->Level = level;
  ms->Mode = mode;
  ms->Side = side;
  if(oms) {
    TetsurfGetRange(oms->Field,oms->Crystal,mn,mx,ms->Range);
    copy3f(mn,ms->ExtentMin); /* this is not exactly correct...should actually take vertex points from range */
    copy3f(mx,ms->ExtentMax);
    ms->ExtentFlag = true;
  }
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

int ObjectSurfaceSetLevel(ObjectSurface *I,float level,int state)
{
  int a;
  int ok=true;
  int once_flag=true;
  ObjectSurfaceState *ms;
  if(state>=I->NState) {
    ok=false;
  } else {
    for(a=0;a<I->NState;a++) {
      if(state<0) {
        once_flag=false;
      }
      if(!once_flag) {
        state = a;
      }
      ms = I->State + state;
      if(ms->Active) {
        ms->ResurfaceFlag=true;
        ms->RefreshFlag=true;
        ms->Level = level;
      }
      if(once_flag) {
        break;
      }
    }
  }
  return(ok);
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

