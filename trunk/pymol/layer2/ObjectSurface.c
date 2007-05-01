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
#include"P.h"
#include"Util.h"
#include"PyMOLGlobals.h"
#include"Matrix.h"

ObjectSurface *ObjectSurfaceNew(PyMOLGlobals *G);

static void ObjectSurfaceFree(ObjectSurface *I);
void ObjectSurfaceStateInit(PyMOLGlobals *G,ObjectSurfaceState *ms);
void ObjectSurfaceRecomputeExtent(ObjectSurface *I);

#ifndef _PYMOL_NOPY
static PyObject *ObjectSurfaceStateAsPyList(ObjectSurfaceState *I)
{
  PyObject *result = NULL;

  result = PyList_New(16);
  
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
  ObjectNameType MapName;
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

static int ObjectSurfaceStateFromPyList(PyMOLGlobals *G,ObjectSurfaceState *I,PyObject *list)
{
  int ok=true;
  int ll=0;
  PyObject *tmp;
  if(ok) ok=(list!=NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active=false;
    else {
      ObjectSurfaceStateInit(G,I);
      if(ok) ok=PyList_Check(list);
      if(ok) ll=PyList_Size(list);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->Active);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,1),I->MapName,WordLength);
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
      ok = ObjectSurfaceStateFromPyList(I->Obj.G,I->State+a,PyList_GetItem(list,a));
      if(!ok) break;
    }
  }
  return(ok);
}
#endif

int ObjectSurfaceNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectSurface **result)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  ObjectSurface *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);

  I=ObjectSurfaceNew(G);
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectFromPyList(G,PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NState);
  if(ok) ok = ObjectSurfaceAllStatesFromPyList(I,PyList_GetItem(list,2));
  if(ok) {
    (*result) = I;
    ObjectSurfaceRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return(ok);
#endif
}


PyObject *ObjectSurfaceAsPyList(ObjectSurface *I)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  
  PyObject *result=NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NState));
  PyList_SetItem(result,2,ObjectSurfaceAllStatesAsPyList(I));

  return(PConvAutoNone(result));  
#endif
}

static void ObjectSurfaceStateFree(ObjectSurfaceState *ms)
{
  ObjectStatePurge(&ms->State);
  if(ms->State.G->HaveGUI) {
    if(ms->displayList) {
      if(PIsGlutThread()) {
        if(ms->State.G->ValidContext) {
          glDeleteLists(ms->displayList,1);
          ms->displayList = 0;
        }
      } else {
        char buffer[255]; /* pass this off to the main thread */
        sprintf(buffer,"_cmd.gl_delete_lists(%d,%d)\n",ms->displayList,1);
        PParse(ms->State.G,buffer);
      }
    }
  }
  
  
  VLAFreeP(ms->N);
  VLAFreeP(ms->V);
  FreeP(ms->VC);
  FreeP(ms->RC);
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
    ErrMessage(I->Obj.G,"ObjectSurfaceDump","can't open file for writing");
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
    PRINTFB(I->Obj.G,FB_ObjectSurface,FB_Actions)
      " ObjectSurfaceDump: %s written to %s\n",I->Obj.Name,fname
      ENDFB(I->Obj.G);
  }
}


static void ObjectSurfaceInvalidate(ObjectSurface *I,int rep,int level,int state)
{
  int a;
  int once_flag=true;
  if(level>=cRepInvExtents) {
    I->Obj.ExtentFlag=false;
  }
  for(a=0;a<I->NState;a++) {
    if(state<0) once_flag=false;
    if(!once_flag) state=a;
    I->State[state].RefreshFlag=true;
    if(level>=cRepInvAll) {
      I->State[state].ResurfaceFlag=true;     
      SceneChanged(I->Obj.G);
    } else if(level>=cRepInvColor) {
      I->State[state].RecolorFlag=true;
      SceneChanged(I->Obj.G);
    } else {
      SceneInvalidate(I->Obj.G);
    }
    if(once_flag) break;
  }
}

int ObjectSurfaceInvalidateMapName(ObjectSurface *I,char *name)
{
  int a;
  ObjectSurfaceState *ms;
  int result=false;
  for(a=0;a<I->NState;a++) {
    ms=I->State+a;
    if(ms->Active) {
      if(strcmp(ms->MapName,name)==0) {
        ObjectSurfaceInvalidate(I,cRepAll,cRepInvAll,a);
        result=true;
      }
    }
  }
  return result;
}

static void ObjectSurfaceStateUpdateColors(ObjectSurface *I, ObjectSurfaceState *ms)
{ 
  int one_color_flag = true;
  int cur_color = SettingGet_color(I->Obj.G, I->Obj.Setting, NULL, cSetting_surface_color);

  if(cur_color == -1)
    cur_color = I->Obj.Color;

  if(ColorCheckRamped(I->Obj.G,cur_color))
    one_color_flag = false;

  ms->OneColor = cur_color;
  if(ms->V) {
    int ramped_flag = false;
    float *v = ms->V;
    float *vc;
    int *rc;
    int a;
    int state = ms - I->State;
    int base_n_vert = ms->base_n_V;
    switch(ms->Mode) {
    case 3:
    case 2:
      {
        int n_vert = VLAGetSize(ms->V)/6;
        base_n_vert /=6;

        if(ms->VC && (ms->VCsize<n_vert)) {
          FreeP(ms->VC);
          FreeP(ms->RC);
        }
        
        if(!ms->VC) {
          ms->VCsize = n_vert;
          ms->VC = Alloc(float,n_vert*3);
        }
        if(!ms->RC) {
          ms->RC = Alloc(int,n_vert);
        }
        rc = ms->RC;
        vc = ms->VC;
        v+=3;
        if(vc) {
          for(a=0;a<n_vert;a++) {
            if(a==base_n_vert) {
              int new_color = SettingGet_color(I->Obj.G, I->Obj.Setting, 
                                               NULL, cSetting_surface_negative_color);
              if(new_color==-1) new_color=cur_color;
              if(new_color!=cur_color) {
                one_color_flag = false;
                cur_color = new_color;
              }
            }
            if(ColorCheckRamped(I->Obj.G,cur_color)) {
              ColorGetRamped(I->Obj.G,cur_color,v,vc,state);
              *rc = cur_color;
              ramped_flag=true;
            } else {
              float *col = ColorGet(I->Obj.G,cur_color);
              copy3f(col,vc);
            }
            rc++;
            vc+=3;
            v+=6; /* alternates with normals */
          }
        }
      }
      break;
    case 1:
    case 0:
    default:
      {
        int n_vert = VLAGetSize(ms->V)/3;
        base_n_vert /=3;
        if(ms->VC && (ms->VCsize<n_vert)) {
          FreeP(ms->VC);
          FreeP(ms->RC);
        }
        
        if(!ms->VC) {
          ms->VCsize = n_vert;
          ms->VC = Alloc(float,n_vert*3);
        }
        if(!ms->RC) {
          ms->RC = Alloc(int,n_vert);
        }
        rc = ms->RC;
        vc = ms->VC;
        if(vc) {
          for(a=0;a<n_vert;a++) {
            if(a==base_n_vert) {
              int new_color = SettingGet_color(I->Obj.G, I->Obj.Setting, 
                                               NULL, cSetting_surface_negative_color);
              if(new_color==-1) new_color=cur_color;
              if(new_color!=cur_color)
                one_color_flag = false;
                cur_color = new_color;
            }

            if(ColorCheckRamped(I->Obj.G,cur_color)) {
              ColorGetRamped(I->Obj.G,cur_color,v,vc,state);
              *rc = cur_color;
              ramped_flag=true;
            } else {
              float *col = ColorGet(I->Obj.G,cur_color);
              copy3f(col,vc);
            }
            rc++;
            vc+=3;
            v+=3; 
          }
        }
      }
      break;
    }

    if(one_color_flag && (!ramped_flag)) {
      FreeP(ms->VC);
      FreeP(ms->RC);
    } else if( (!ramped_flag) || (!SettingGet_b(I->Obj.G,NULL,I->Obj.Setting,cSetting_ray_color_ramps))) {
      FreeP(ms->RC);
    }
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
  float carve_buffer;
  for(a=0;a<I->NState;a++) {
    ms = I->State+a;
    if(ms->Active) {
      map = ExecutiveFindObjectMapByName(I->Obj.G,ms->MapName);
      if(!map) {
        ok=false;
        PRINTFB(I->Obj.G,FB_ObjectSurface,FB_Errors)
          "ObjectSurfaceUpdate-Error: map '%s' has been deleted.\n",ms->MapName
          ENDFB(I->Obj.G);
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
          
          if(oms->State.Matrix) {
            ObjectStateSetMatrix(&ms->State,oms->State.Matrix);
          } else if(ms->State.Matrix) {
            ObjectStateResetMatrix(&ms->State);
          }
          ms->displayListInvalid = true;          
          ms->RefreshFlag=false;
        }
      }
      if(map&&ms&&oms&&ms->N&&ms->V&&I->Obj.RepVis[cRepSurface]) {
        if(ms->ResurfaceFlag) {
          ms->ResurfaceFlag=false;
          ms->RecolorFlag=true;
          PRINTFB(I->Obj.G,FB_ObjectSurface,FB_Details)
           " ObjectSurface: updating \"%s\".\n" , I->Obj.Name 
            ENDFB(I->Obj.G);
          if(oms->Field) {

            {
              float *min_ext,*max_ext;
              float tmp_min[3],tmp_max[3];
              if(MatrixInvTransformExtentsR44d3f(ms->State.Matrix,
                                                 ms->ExtentMin,ms->ExtentMax,
                                                 tmp_min,tmp_max)) {
                min_ext = tmp_min;
                max_ext = tmp_max;
              } else {
                min_ext = ms->ExtentMin;
                max_ext = ms->ExtentMax;
              }
              
              TetsurfGetRange(I->Obj.G,oms->Field,oms->Crystal,
                              min_ext,max_ext,ms->Range);
            }

            if(ms->CarveFlag&&ms->AtomVertex) {
              carve_buffer = ms->CarveBuffer;
              if(carve_buffer<0.0F) {
                carve_buffer = -carve_buffer;
              }

              voxelmap=MapNew(I->Obj.G,-carve_buffer,ms->AtomVertex,
                              VLAGetSize(ms->AtomVertex)/3,NULL);
              if(voxelmap)
                MapSetupExpress(voxelmap);  
            }

            ms->nT=TetsurfVolume(I->Obj.G,oms->Field,
                                 ms->Level,
                                 &ms->N,&ms->V,
                                 ms->Range,
                                 ms->Mode,
                                 voxelmap,
                                 ms->AtomVertex,
                                 ms->CarveBuffer,
                                 ms->Side);
            
            if(!SettingGet_b(I->Obj.G,I->Obj.Setting,NULL,cSetting_surface_negative_visible)) { 
              ms->base_n_V = VLAGetSize(ms->V);
            } else {
              /* do we want the negative surface too? */

              int nT2;
              int *N2 = VLAlloc(int,10000);
              float *V2 =  VLAlloc(float,10000);

              nT2 = TetsurfVolume(I->Obj.G,oms->Field,
                                  -ms->Level,
                                  &N2,&V2,
                                  ms->Range,
                                  ms->Mode,
                                  voxelmap,
                                  ms->AtomVertex,
                                  ms->CarveBuffer,
                                  ms->Side);
              if(N2&&V2) {
                
                int base_n_N = VLAGetSize(ms->N);
                int base_n_V = VLAGetSize(ms->V);
                int addl_n_N = VLAGetSize(N2);
                int addl_n_V = VLAGetSize(V2);

                ms->base_n_V = base_n_V;

                /* make room */

                VLASize(ms->N, int, base_n_N + addl_n_N);
                VLASize(ms->V, float, base_n_V + addl_n_V);
                
                /* copy vertex data */

                memcpy(((char*)ms->V)+(sizeof(float)*base_n_V),
                       V2, sizeof(float)*addl_n_V);

                /* copy strip counts */

                memcpy(((char*)ms->N)+(sizeof(int)*(base_n_N-1)),
                       N2, sizeof(int)*addl_n_N);
                ms->N[ base_n_N + addl_n_N - 1] = 0;

                ms->nT += nT2;
                VLAFreeP(N2);
                VLAFreeP(V2);
              }

            }
            
            if(voxelmap)
              MapFree(voxelmap);

            if(ms->State.Matrix) { /* in we're in a different reference frame...*/
              double *matrix = ms->State.Matrix;
              float *v;
              int *n,c;
              
              v=ms->V;
              n=ms->N;

              if(n&&v) {

                while(*n) {
                  c=*(n++);
                  switch(ms->Mode) {
                  case 3:
                  case 2:
                    transform44d3fas33d3f(matrix,v,v);
                    transform44d3f(matrix,v+3,v+3);
                    transform44d3fas33d3f(matrix,v+6,v+6);
                    transform44d3f(matrix,v+9,v+9);
                    v+=12;
                    c-=4;
                    while(c>0) {
                      transform44d3fas33d3f(matrix,v,v);
                      transform44d3f(matrix,v+3,v+3);
                      v+=6;
                      c-=2;
                    }
                    break;
                  case 1:
                    transform44d3f(matrix,v,v);
                    c--;
                    v+=3;
                    while(c>0) {
                      transform44d3f(matrix,v,v);
                      v+=3;
                      c--;
                    }
                    break;
                  case 0:
                  default:
                    while(c>0) {
                      transform44d3f(matrix,v,v);
                      v+=3;
                      c--;
                    }
                    break;
                  }
                }
              }
            }
          }
        }
        if(ms->RecolorFlag) {
          ObjectSurfaceStateUpdateColors(I,ms);
          ms->RecolorFlag=false;
        }
      }
    }
  }
  if(!I->Obj.ExtentFlag) {
    ObjectSurfaceRecomputeExtent(I);
  }
  SceneInvalidate(I->Obj.G);
}
#if 0
static int ZOrderFn(float *array,int l,int r)
{
  return (array[l]<=array[r]);
}

static int ZRevOrderFn(float *array,int l,int r)
{
  return (array[l]>=array[r]);
}
#endif
static void ObjectSurfaceRender(ObjectSurface *I,RenderInfo *info)
{
  PyMOLGlobals *G = I->Obj.G;
  int state = info->state;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int pass = info->pass;
  float *v = NULL;
  float *vc = NULL;
  int *rc = NULL;
  float *col;
  int *n = NULL;
  int c;
  int a=0;
  ObjectSurfaceState *ms = NULL;
  float alpha;

  ObjectPrepareContext(&I->Obj,ray);

  alpha = SettingGet_f(G,NULL,I->Obj.Setting,cSetting_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;
  
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
           ((SettingGet(G,cSetting_static_singletons)&&(I->NState==1))))
          ms=I->State;
      }
    }
    if(ms) {
      if(ms->Active&&ms->V&&ms->N) {
        v=ms->V;
        n=ms->N;
        if(ray) {
          if(ms->UnitCellCGO&&(I->Obj.RepVis[cRepCell]))
            CGORenderRay(ms->UnitCellCGO,ray,ColorGet(G,I->Obj.Color),
                         I->Obj.Setting,NULL);

          ray->fTransparentf(ray,1.0F-alpha);       
          ms->Radius=SettingGet_f(G,I->Obj.Setting,NULL,cSetting_mesh_radius);
          if(ms->Radius==0.0F) {
            ms->Radius = ray->PixelRadius*
              SettingGet_f(I->Obj.G,I->Obj.Setting,NULL,cSetting_mesh_width)/2.0F;
          }
          
          if(n&&v&&I->Obj.RepVis[cRepSurface]) {
            float cc[3];
            float colA[3],colB[3],colC[3];
            ColorGetEncoded(G,ms->OneColor,cc);
            vc = ms->VC;

            rc = ms->RC;
            while(*n)
              {
                c=*(n++);
                switch(ms->Mode) {
                case 3:
                case 2:
                  v+=12;
                  if(vc) vc+=6;
                  c-=4;
                  while(c>0) {
                    if(vc) {
                      register float *cA=vc-6, *cB=vc-3, *cC=vc;
                      if(rc) {
                        if(rc[0]<-1) ColorGetEncoded(G,rc[0], (cA=colA));
                        if(rc[1]<-1) ColorGetEncoded(G,rc[1], (cB=colB));
                        if(rc[2]<-1) ColorGetEncoded(G,rc[2], (cC=colC));
                        rc++;
                      }
                      ray->fTriangle3fv(ray,v-9,v-3,v+3,
                                        v-12,v-6,v,
                                        cA,cB,cC);
                      vc+=3;
                    } else {
                      ray->fTriangle3fv(ray,v-9,v-3,v+3,
                                        v-12,v-6,v,
                                        cc,cc,cc);
                    }
                    v+=6;
                    c-=2;
                  }
                  break;
                case 1:
                  c--;
                  v+=3;
                  if(vc) vc+=3;
                  while(c>0) {
                    if(vc) {
                      register float *cA=vc-3, *cB=vc;
                      if(rc) {
                        if(rc[0]<-1) ColorGetEncoded(G,rc[0], (cA=colA));
                        if(rc[1]<-1) ColorGetEncoded(G,rc[1], (cB=colB));
                        rc++;
                      }
                      ray->fSausage3fv(ray,v-3,v,ms->Radius,cA,cB);
                      vc+=3;
                    } else 
                      ray->fSausage3fv(ray,v-3,v,ms->Radius,cc,cc);
                    v+=3;
                    c--;
                  }
                  break;
                case 0:
                default:
                  while(c>0) {
                    if(vc) {
                      ray->fColor3fv(ray,vc);
                      vc+=3;
                    }
                    ray->fSphere3fv(ray,v,ms->Radius);
                    v+=3;
                    c--;
                  }
                  break;
                }
              }
            
          }
          ray->fTransparentf(ray,0.0);
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
          } else {
            
            int render_now = false;

            if(alpha>0.0001) {
              render_now = (pass==-1);
            } else 
              render_now = (!pass);
            if(render_now) {
              int use_dlst;
              
              if(ms->UnitCellCGO&&(I->Obj.RepVis[cRepCell]))
                CGORenderGL(ms->UnitCellCGO,ColorGet(G,I->Obj.Color),
                            I->Obj.Setting,NULL,info);
              
              SceneResetNormal(G,false);
              col = ColorGet(G,ms->OneColor);
              glColor4f(col[0],col[1],col[2],alpha);
              
              use_dlst = (int)SettingGet(G,cSetting_use_display_lists);

              if(use_dlst && ms->displayList && ms->displayListInvalid) {
                glDeleteLists(ms->displayList,1);
                ms->displayList = 0;
                ms->displayListInvalid = false;
              }

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
                  
                  if((ms->Mode>1)&&(alpha!=1.0)) { /* transparent */
                    
                    int t_mode;

                    t_mode  = SettingGet_i(G,NULL,I->Obj.Setting,cSetting_transparency_mode);
                    
                    if(t_mode) { /* high quality (sorted) transparency? */
                      
                      float **t_buf=NULL,**tb;
                      float **c_buf=NULL,**tc;
                      float *z_value=NULL,*zv;
                      int *ix=NULL;
                      int n_tri = 0;
                      float sum[3];
                      float matrix[16];
                      int parity;
                      glGetFloatv(GL_MODELVIEW_MATRIX,matrix);
                      
                      t_buf = Alloc(float*,ms->nT*9);
                      vc = ms->VC;

                      if(vc) {
                        c_buf = Alloc(float*,ms->nT*9);
                      }

                      z_value = Alloc(float,ms->nT);
                      ix = Alloc(int,ms->nT);
                      
                      zv = z_value;
                      tb = t_buf;
                      tc = c_buf;

                      while(*n)
                        {
                          parity=true;
                          c=*(n++);
                          v+=12;
                          if(vc)
                            vc+=6;
                          c-=4;
                          while(c>0) {
                            
                            if(parity) {
                              *(tb++) = v-12;
                              *(tb++) = v-9;
                              *(tb++) = v-6;
                              *(tb++) = v-3;
                              *(tb++) = v;
                              *(tb++) = v+3;
                              if(vc) {
                                *(tc++) = vc-6;
                                *(tc++) = vc-3;
                                *(tc++) = vc;
                              }
                            } else {
                              *(tb++) = v-12;
                              *(tb++) = v-9;
                              *(tb++) = v;
                              *(tb++) = v+3;
                              *(tb++) = v-6;
                              *(tb++) = v-3;
                              if(vc) {
                                *(tc++) = vc-6;
                                *(tc++) = vc;
                                *(tc++) = vc-3;
                              }
                            }
                            
                            parity=!parity;
                            
                            add3f(tb[-1],tb[-3],sum);
                            add3f(sum,tb[-5],sum);
                            
                            *(zv++) = matrix[2]*sum[0]+matrix[6]*sum[1]+matrix[10]*sum[2];
                            n_tri++;
                            
                            if(vc) vc+=3;
                            v+=6;
                            c-=2;
                          }
                        }


                      switch(t_mode) {
                      case 1:
                        UtilSemiSortFloatIndex(n_tri,z_value,ix,true);
                        /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn);*/
                        break;
                      default:
                        UtilSemiSortFloatIndex(n_tri,z_value,ix,false);
                        /*UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
                        break;
                      }
                      
                      c=n_tri;
                      
                      col=ColorGet(G,ms->OneColor);
                      
                      glColor4f(col[0],col[1],col[2],alpha);
                      glBegin(GL_TRIANGLES);
                      for(c=0;c<n_tri;c++) {
                        
                        tb = t_buf+6*ix[c];
                        if(vc) tc = c_buf+3*ix[c];

                        if(vc) { glColor4f(tc[0][0],tc[0][1],tc[0][2],alpha); tc++;}
                        glNormal3fv(*(tb++));
                        glVertex3fv(*(tb++));
                        if(vc) { glColor4f(tc[0][0],tc[0][1],tc[0][2],alpha); tc++;}
                        glNormal3fv(*(tb++));
                        glVertex3fv(*(tb++));
                        if(vc) { glColor4f(tc[0][0],tc[0][1],tc[0][2],alpha); tc++;}
                        glNormal3fv(*(tb++));
                        glVertex3fv(*(tb++));
                      }
                      glEnd();
                      
                      FreeP(ix);
                      FreeP(z_value);
                      FreeP(t_buf);
                      FreeP(c_buf);
                    } else {  /* fast, but unoptimized transparency */
                      vc = ms->VC;
                      while(*n)
                        {
                          c=*(n++);
                          
                          glBegin(GL_TRIANGLE_STRIP);
                          while(c>0) {
                            glNormal3fv(v);
                            v+=3;
                            if(vc) {
                              glColor3fv(vc);
                              vc+=3;
                            }
                            glVertex3fv(v);
                            v+=3;
                            c-=2;
                          }
                          glEnd();
                        }
                    }
                  } else {
                    glLineWidth(SettingGet_f(G,I->Obj.Setting,NULL,cSetting_mesh_width));
                    vc = ms->VC;
                    while(*n)
                      {
                        c=*(n++);
                        switch(ms->Mode) {
                        case 3:
                        case 2:
                          glBegin(GL_TRIANGLE_STRIP);
                          while(c>0) {
                            glNormal3fv(v);
                            v+=3;
                            if(vc) {
                              glColor3fv(vc);
                              vc+=3;
                            }
                            glVertex3fv(v);
                            v+=3;
                            c-=2;
                          }
                          glEnd();
                          break;
                        case 1:
                          glBegin(GL_LINES);
                          while(c>0) {
                            if(vc) {
                              glColor3fv(vc);
                              vc+=3;
                            }
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
                            if(vc) {
                              glColor3fv(vc);
                              vc+=3;
                            }
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
              }
            
              if(use_dlst&&ms->displayList) {
                glEndList();
              }
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
ObjectSurface *ObjectSurfaceNew(PyMOLGlobals *G)
{
  OOAlloc(G,ObjectSurface);
  
  ObjectInit(G,(CObject*)I);
  
  I->NState = 0;
  I->State=VLAMalloc(10,sizeof(ObjectSurfaceState),5,true); /* autozero important */

  I->Obj.type = cObjectSurface;
  
  I->Obj.fFree = (void (*)(struct CObject *))ObjectSurfaceFree;
  I->Obj.fUpdate =  (void (*)(struct CObject *)) ObjectSurfaceUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, RenderInfo *info))ObjectSurfaceRender;
  I->Obj.fInvalidate =(void (*)(struct CObject *,int,int,int))ObjectSurfaceInvalidate;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectSurfaceGetNStates;
  return(I);
}

/*========================================================================*/
void ObjectSurfaceStateInit(PyMOLGlobals *G,ObjectSurfaceState *ms)
{
  if(ms->Active)
    ObjectStatePurge(&ms->State);
  ObjectStateInit(G,&ms->State);
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
  ms->nT=0;
  ms->VC = NULL;
  ms->RC = NULL;
  ms->VCsize = 0;
  ms->Active=true;
  ms->ResurfaceFlag=true;
  ms->RecolorFlag=false;
  ms->ExtentFlag=false;
  ms->CarveFlag=false;
  ms->AtomVertex=NULL;
  ms->UnitCellCGO=NULL;
  ms->Side = 0;
  ms->displayList = 0;
  ms->displayListInvalid = true;
}

/*========================================================================*/
ObjectSurface *ObjectSurfaceFromBox(PyMOLGlobals *G,ObjectSurface *obj,ObjectMap *map,
                                    int map_state,
int state,float *mn,float *mx,float level,int mode,
float carve,float *vert_vla,int side)
{
  ObjectSurface *I;
  ObjectSurfaceState *ms;
  ObjectMapState *oms;

  if(!obj) {
    I=ObjectSurfaceNew(G);
  } else {
    I=obj;
  }

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectSurfaceState,state);
    I->NState=state+1;
  }

  ms=I->State+state;
  ObjectSurfaceStateInit(G,ms);

  strcpy(ms->MapName,map->Obj.Name);
  ms->MapState = map_state;
  oms = ObjectMapGetState(map,map_state);

  ms->Level = level;
  ms->Mode = mode;
  ms->Side = side;
  if(oms) {

    if(oms->State.Matrix) {
      ObjectStateSetMatrix(&ms->State,oms->State.Matrix);
    } else if(ms->State.Matrix) {
      ObjectStateResetMatrix(&ms->State);
    }

    copy3f(mn,ms->ExtentMin); /* this is not exactly correct...should actually take vertex points from range */
    copy3f(mx,ms->ExtentMax);
    
    {
      float *min_ext,*max_ext;
      float tmp_min[3],tmp_max[3];
      if(MatrixInvTransformExtentsR44d3f(ms->State.Matrix,
                                         ms->ExtentMin,ms->ExtentMax,
                                         tmp_min,tmp_max)) {
        min_ext = tmp_min;
        max_ext = tmp_max;
      } else {
        min_ext = ms->ExtentMin;
        max_ext = ms->ExtentMax;
      }
      
      TetsurfGetRange(G,oms->Field,oms->Crystal,min_ext,max_ext,ms->Range);
    }
    ms->ExtentFlag = true;
  }
  if(carve!=0.0) {
    ms->CarveFlag=true;
    ms->CarveBuffer = carve;
    ms->AtomVertex = vert_vla;

    if(ms->State.Matrix) {
      int n=VLAGetSize(ms->AtomVertex)/3;
      float *v = ms->AtomVertex;
      double *matrix = ms->State.Matrix;
      while(n--) {
        /* convert back into original map coordinates 
           for surface carving operation */
        inverse_transform44d3f(matrix, v,v);
        v+=3;
      }
    }

  }
  if(I) {
    ObjectSurfaceRecomputeExtent(I);
  }
  I->Obj.ExtentFlag=true;
  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]);*/
  SceneChanged(G);
  SceneCountFrames(G);
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

  if(I->Obj.TTTFlag && I->Obj.ExtentFlag) {
    float *ttt;
    double tttd[16];
    if(ObjectGetTTT(&I->Obj,&ttt,-1)) {
      convertTTTfR44d(ttt,tttd);
      MatrixTransformExtentsR44d3f(tttd,
                                   I->Obj.ExtentMin,I->Obj.ExtentMax,
                                   I->Obj.ExtentMin,I->Obj.ExtentMax);
    }
  }
}

