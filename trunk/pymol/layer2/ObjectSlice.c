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
-* Filipe Maia 
-*
Z* -------------------------------------------------------------------
*/

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectSlice.h"
#include"Base.h"
#include"Matrix.h"
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
#include"Executive.h"
#include"PConv.h"
#include"P.h"
#include"Text.h"
#include"Util.h"
#include"ButMode.h"
#include"ObjectGadgetRamp.h"

typedef struct ObjRec {
  struct CObject *obj;  
  struct ObjRec *next;
} ObjRec;

ObjectSlice *ObjectSliceNew(void);

static void ObjectSliceFree(ObjectSlice *I);
void ObjectSliceStateInit(ObjectSliceState *ms);
void ObjectSliceRecomputeExtent(ObjectSlice *I);

static PyObject *ObjectSliceStateAsPyList(ObjectSliceState *I)
{
  PyObject *result = NULL;

  result = PyList_New(10);
  
  PyList_SetItem(result,0,PyInt_FromLong(I->Active));
  PyList_SetItem(result,1,PyString_FromString(I->MapName));
  PyList_SetItem(result,2,PyInt_FromLong(I->MapState));
  PyList_SetItem(result,3,PConvFloatArrayToPyList(I->ExtentMin,3));
  PyList_SetItem(result,4,PConvFloatArrayToPyList(I->ExtentMax,3));
  PyList_SetItem(result,5,PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result,6,PConvFloatArrayToPyList(I->origin,3));
  PyList_SetItem(result,7,PConvFloatArrayToPyList(I->system,9));
  PyList_SetItem(result,8,PyFloat_FromDouble(I->MapMean));
  PyList_SetItem(result,9,PyFloat_FromDouble(I->MapStdev));
  
#if 0
  int Active;
  char MapName[ObjNameMax];
  int MapState;
  float ExtentMin[3];
  float ExtentMax[3];
  int ExtentFlag;

  float origin[3]; /* the origin of the plane */
  float system[9]; /* x, y, and z of the system */
  float grid;   /* sampling interval for the map */
#endif

  return(PConvAutoNone(result));  
}

static PyObject *ObjectSliceAllStatesAsPyList(ObjectSlice *I)
{
  
  PyObject *result=NULL;
  int a;
  result = PyList_New(I->NState);
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result,a,ObjectSliceStateAsPyList(I->State+a));
    } else {
      PyList_SetItem(result,a,PConvAutoNone(NULL));
    }
  }
  return(PConvAutoNone(result));  

}

static int ObjectSliceStateFromPyList(ObjectSliceState *I,PyObject *list)
{
  int ok=true;

  int ll;
  if(ok) ok=(list!=NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active=false;
    else {
      ObjectSliceStateInit(I);
      if(ok) ok=(list!=NULL);
      if(ok) ok=PyList_Check(list);
      if(ok) ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */      

      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->Active);
      if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,1),I->MapName,ObjNameMax);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->MapState);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,3),I->ExtentMin,3);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,4),I->ExtentMax,3);
      if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,5),&I->ExtentFlag);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,6),I->origin,3);
      if(ok) ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list,7),I->system,9);
      if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,8),&I->MapMean);
      if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,9),&I->MapStdev);

      I->RefreshFlag=true;
    }
  }

  return(ok);
}

static int ObjectSliceAllStatesFromPyList(ObjectSlice *I,PyObject *list)
{
  int ok=true;
  int a;
  VLACheck(I->State,ObjectSliceState,I->NState);
  if(ok) ok=PyList_Check(list);
  if(ok) {
    for(a=0;a<I->NState;a++) {
      ok = ObjectSliceStateFromPyList(I->State+a,PyList_GetItem(list,a));
      if(!ok) break;
    }
  }
  return(ok);
}

int ObjectSliceNewFromPyList(PyObject *list,ObjectSlice **result)
{
  int ok = true;
  int ll;
  ObjectSlice *I=NULL;
  (*result) = NULL;
  
  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  I=ObjectSliceNew();
  if(ok) ok = (I!=NULL);
  
  if(ok) ok = ObjectFromPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NState);
  if(ok) ok = ObjectSliceAllStatesFromPyList(I,PyList_GetItem(list,2));
  if(ok) {
    (*result) = I;
    ObjectSliceRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return(ok);
}



PyObject *ObjectSliceAsPyList(ObjectSlice *I)
{
  
  PyObject *result=NULL;

  result = PyList_New(3);
  PyList_SetItem(result,0,ObjectAsPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NState));
  PyList_SetItem(result,2,ObjectSliceAllStatesAsPyList(I));

  return(PConvAutoNone(result));  
}

static void ObjectSliceStateFree(ObjectSliceState *oss)
{
  if(PMGUI) {
    if(oss->displayList) {
      if(PIsGlutThread()) {
        glDeleteLists(oss->displayList,1);
        oss->displayList = 0;
      } else {
        char buffer[255]; /* pass this off to the main thread */
        sprintf(buffer,"_cmd.gl_delete_lists(%d,%d)\n",oss->displayList,1);
        PParse(buffer);
      }
    }
  }
  VLAFreeP(oss->normals);
  VLAFreeP(oss->colors);
  VLAFreeP(oss->values);
  VLAFreeP(oss->points);
  VLAFreeP(oss->flags);
  VLAFreeP(oss->strips);
}

static void ObjectSliceFree(ObjectSlice *I) {
  int a;
  for(a=0;a<I->NState;a++) {
    if(I->State[a].Active)
      ObjectSliceStateFree(I->State+a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  
  OOFreeP(I);
}

static void ObjectSliceInvalidate(ObjectSlice *I,int rep,int level,int state)
{
  int a;
  int once_flag=true;
  for(a=0;a<I->NState;a++) {
    if(state<0) once_flag=false;
    if(!once_flag) state=a;
    I->State[state].RefreshFlag=true;
    SceneChanged();
    if(once_flag) break;
  }
}


static void ObjectSliceStateAssignColors(ObjectSliceState *oss, ObjectGadgetRamp *ogr)
{
    /* compute the colors */
  if(oss && oss->values && oss->colors) {
    int x,y;
    int *min = oss->min;
    int *max = oss->max;
    float *value = oss->values;
    int *flag = oss->flags;
    float *color = oss->colors;
    for(y=min[1];y<=max[1];y++) {
      for(x=min[0];x<=max[0];x++) {
        if(*flag) {
          ObjectGadgetRampInterpolate(ogr, *value, color);
        }
        color +=3;
        value++;
        flag++;
      }
    }
  }
}


static void ObjectSliceStateUpdate(ObjectSlice *I,ObjectSliceState *oss, ObjectMapState *oms)
{
  int ok=true;
  int min[2] = {0,0}, max[2] = {0,0}; /* limits of the rectangle */
  int need_normals = false;
  float grid = SettingGet_f(NULL,I->Obj.Setting,cSetting_slice_grid);  
  if(grid<0.01F)
    grid=0.01F;

  /* for the given map, compute a new set of interpolated point with accompanying levels */

  /* first, find the limits of the enclosing rectangle, starting at the slice origin, 
     via a simple brute-force approach... */

  if(ok) {
    int size = 1, minus_size;
    int a;
    int keep_going=true;
    float point[3];

    while(keep_going) {
      keep_going=false;
      minus_size = -size;

      for(a=-size;a<=size;a++) {
        
        if(max[1]!=size) {
          point[0] = grid*a;
          point[1] = grid*size;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms,point)) {
              keep_going = true;
              if(max[1]<size) max[1]=size;
          }
        }

        if(min[1]!=minus_size) {
          point[0] = grid*a;
          point[1] = grid*minus_size;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          
          if(ObjectMapStateContainsPoint(oms,point)) {
            keep_going = true;
            if(min[1]>minus_size) min[1]=minus_size;
          }
        }
        
        if(max[0]!=size) {
          point[0] = grid*size;
          point[1] = grid*a;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms,point)) {
            keep_going = true;
            if(max[0]<size) max[0]=size;
          }
        }
        
        if(min[0]!=minus_size) {
          point[0] = grid*minus_size;
          point[1] = grid*a;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms,point)) {
            keep_going = true;
            if(min[0]>minus_size) min[0]=minus_size;
          }
        }
      }
      size++;
    }
    oss->max[0] = max[0];
    oss->max[1] = max[1];
    oss->min[0] = min[0];
    oss->min[1] = min[1];
  }

  /* now confirm that storage is available */
  if (ok) {
    int n_alloc = (1+oss->max[0]-oss->min[0])*(1+oss->max[1]-oss->min[1]);

    if(!oss->points) {
      oss->points = VLAlloc(float, n_alloc*3);
    } else {
      VLACheck(oss->points, float, n_alloc*3); /* note: this is a macro which reassigns the pointer */
    }

    if(!oss->values) {
      oss->values = VLAlloc(float, n_alloc);
    } else {
      VLACheck(oss->values, float, n_alloc); /* note: this is a macro which reassigns the pointer */
    }

    if(!oss->colors) {
      oss->colors = VLACalloc(float, n_alloc*3);
    } else {
      VLACheck(oss->colors, float, n_alloc*3); /* note: this is a macro which reassigns the pointer */
    }

    if(!oss->flags) {
      oss->flags = VLAlloc(int, n_alloc);
    } else {
      VLACheck(oss->flags, int, n_alloc); /* note: this is a macro which reassigns the pointer */
    }
    if(!(oss->points&&oss->values&&oss->flags)) {
      ok=false;
      PRINTFB(FB_ObjectSlice,FB_Errors)
        "ObjectSlice-Error: allocation failed\n"
        ENDFB;
    }

    if(!oss->strips) /* this is range-checked during use */
      oss->strips = VLAlloc(int,n_alloc);

    oss->n_points = n_alloc;
  }

  /* generate the coordinates */

  if(ok) {
    int x,y;
    float *point = oss->points;
    for(y=min[1];y<=max[1];y++) {
      for(x=min[0];x<=max[0];x++) {
        point[0] = grid*x;
        point[1] = grid*y;
        point[2] = 0.0F;
        transform33f3f(oss->system, point, point);
        add3f(oss->origin, point, point);
        point +=3;
      }
    }
  }

  /* interpolate and flag the points inside the map */

  if(ok) {
    ObjectMapStateInterpolate(oms,oss->points, oss->values, oss->flags, oss->n_points);
  }

  /* apply the height scale (if nonzero) */

  if(ok) {

    if(SettingGet_b(NULL,I->Obj.Setting,cSetting_slice_height_map))
      {
        float height_scale = SettingGet_f(NULL,I->Obj.Setting,cSetting_slice_height_scale);
        float *value = oss->values;
        float up[3],scaled[3],factor;
        int x,y;
        float *point = oss->points;


        need_normals=true;
        up[0] = oss->system[2];
        up[1] = oss->system[5];
        up[2] = oss->system[8];
        
        for(y=min[1];y<=max[1];y++) {
          for(x=min[0];x<=max[0];x++) {
            factor = ((*value-oss->MapMean)/oss->MapStdev) * height_scale;
            scale3f(up,factor,scaled);
            add3f(scaled,point,point);
            point +=3;
            value ++;
          }
        }
        
      }
  }
/* now generate efficient triangle strips based on the points that are present in the map */

  if(ok) {
    int x,y;
    int cols = 1+max[0]-min[0];
    int flag00,flag01,flag10,flag11;
    int offset = 0, offset00, offset01, offset10, offset11;
    int strip_active = false;
    int n = 0;

    for(y=min[1];y<max[1];y++) {
      offset00 = offset;
      for(x=min[0];x<max[0];x++) {
        
        offset01 = offset00+1;
        offset10 = offset00+cols;
        offset11 = offset10+1;

        flag00 = oss->flags[offset00];
        flag01 = oss->flags[offset01];
        flag10 = oss->flags[offset10];
        flag11 = oss->flags[offset11];

        /* first triangle - forward handedness: 10 00 11 */

        if(strip_active) {
          if(flag10 && flag00 && flag11) {
            /* continue current strip */

            VLACheck(oss->strips,int,n);
            oss->strips[n] = offset10;
            n++;
          } else {
            /* terminate current strip */

            VLACheck(oss->strips,int,n);
            oss->strips[n] = -2;
            strip_active=false;
            n++;
          }
        } else if(flag10 & flag00 && flag11) {
          /* start a new strip with correct parity */
          
          VLACheck(oss->strips,int,n+3);
          oss->strips[n] = -1;
          oss->strips[n+1] = offset10;
          oss->strips[n+2] = offset00;
          oss->strips[n+3] = offset11;
          n+=4;
          strip_active=true;
        }

        /* second triangle -- reverse handedness: 00 11 01*/

        if(strip_active) {
          if(flag00 && flag11 && flag01) {
            /* continue current strip */
            VLACheck(oss->strips,int,n);
            oss->strips[n] = offset01;
            n++;
          } else {
            /* terminate current strip */
            VLACheck(oss->strips,int,n);
            oss->strips[n] = -2;
            strip_active=false;
            n++;
          }
        } else if(flag00 & flag11 && flag01) {
          /* singleton triangle -- improper order for strip */
          
          VLACheck(oss->strips,int,n+5);
          oss->strips[n+0] = -1;
          oss->strips[n+1] = offset11;
          oss->strips[n+2] = offset00;
          oss->strips[n+3] = offset01;
          oss->strips[n+4] = -2;
          n+=5;
        }
        offset00++;
      }
      if(strip_active) {
        /* terminate current strip */
        VLACheck(oss->strips,int,n);
        oss->strips[n] = -2;
        strip_active=false;
        n++;
      }
      offset+=cols;
    }
    VLACheck(oss->strips,int,n);    
    n++;
    oss->n_strips = n;

  }

  /* compute triangle normals if we need them */

  if(!need_normals) {
    VLAFreeP(oss->normals);
  } else {
    int *cnt = NULL;

    if(!oss->normals) {
      oss->normals = VLAlloc(float, oss->n_points*3);
    } else {
      VLACheck(oss->normals, float, oss->n_points*3); /* note: this is a macro which reassigns the pointer */
    }
    cnt = Calloc(int,oss->n_points);

    if(cnt&&oss->normals) {
      int *strip = oss->strips;
      float *point = oss->points;
      float *normal = oss->normals;
      int n=oss->n_strips;
      int a;
      int offset0,offset1,offset2,offset;
      int strip_active =false;
      int tri_count = 0;

      float d1[3],d2[3],cp[3];
      UtilZeroMem(oss->normals,sizeof(float)*3*oss->n_points);

      for(a=0;a<n;a++) {
        offset = *(strip++);
        switch(offset) {
        case -1:
          strip_active=true;
          tri_count = 0;
          break;
        case -2:
          strip_active=false;
          break;
        default:
          if(strip_active) {
            tri_count++;
            offset2 = offset1;
            offset1 = offset0;
            offset0 = offset;
            
            if(tri_count>=3) {
              
              if(tri_count&0x1) { /* get the handedness right ... */
                subtract3f(point+3*offset1,point+3*offset0,d1);
                subtract3f(point+3*offset2,point+3*offset1,d2);
              } else {
                subtract3f(point+3*offset0,point+3*offset1,d1);
                subtract3f(point+3*offset2,point+3*offset0,d2);
              }
              cross_product3f(d2,d1,cp);
              normalize3f(cp);
              add3f(cp,normal+3*offset0,normal+3*offset0);
              add3f(cp,normal+3*offset1,normal+3*offset1);
              add3f(cp,normal+3*offset2,normal+3*offset2);
              cnt[offset0]++;
              cnt[offset1]++;
              cnt[offset2]++;
            }
          }
        }
      }

      { /* now normalize the average normals for active vertices */
        int x,y;
        float *normal = oss->normals;
        int *c = cnt;
        for(y=min[1];y<=max[1];y++) {
          for(x=min[0];x<=max[0];x++) {
            if(*c) 
              normalize3f(normal);
            point +=3;
            c++;
          }
        }
      }
    }
    FreeP(cnt);
  }
}

static void ObjectSliceUpdate(ObjectSlice *I) 
{

  ObjectSliceState *oss;
  ObjectMapState *oms = NULL;
  ObjectMap *map = NULL;
  ObjectGadgetRamp *ogr = NULL;

  int ok=true;
  int a;
  for(a=0;a<I->NState;a++) {
    oss = I->State + a;
    if(oss && oss->Active) {      
      map = ExecutiveFindObjectMapByName(oss->MapName);
      if(!map) {
        ok=false;
        PRINTFB(FB_ObjectSlice,FB_Errors)
          "ObjectSliceUpdate-Error: map '%s' has been deleted.\n",oss->MapName
          ENDFB;
      }
      if(map) {
        oms = ObjectMapGetState(map,oss->MapState);
        if(!oms) ok=false;
      }
      if(oms) {
        
        if(oss->RefreshFlag) {
          oss->RefreshFlag=false;
          PRINTFB(FB_ObjectSlice,FB_Blather)
            " ObjectSlice: updating \"%s\".\n" , I->Obj.Name 
            ENDFB;
          if(oms->Field) {
            ObjectSliceStateUpdate(I,oss, oms);
            ogr = ColorGetRamp(I->Obj.Color);
            if(ogr) 
              ObjectSliceStateAssignColors(oss, ogr);
            else { /* solid color */
              float *solid = ColorGet(I->Obj.Color);
              float *color = oss->colors;
              for(a=0;a<oss->n_points;a++) {
                *(color++)=solid[0];
                *(color++)=solid[1];             
                *(color++)=solid[2];
              }
            }
          }
        }
      }
      SceneDirty();
    }
  }
}

void ObjectSliceDrag(ObjectSlice *I, int state, int mode, float *pt, float *mov, float *z_dir)
{
  ObjectSliceState *oss = NULL;

  if(state>=0) 
    if(state<I->NState) 
      if(I->State[state].Active)
        oss=I->State+state;

  if(oss) {
    switch(mode) {
    case cButModeRotFrag: /* rotated about origin */
      {
        float v3[3];
        float n0[3];
        float n1[3];
        float n2[3];
        float cp[3];
        float mat[9];
        float theta;

        copy3f(oss->origin,v3);
        
        subtract3f(pt,v3,n0);
        add3f(pt,mov,n1);
        subtract3f(n1,v3,n1);

        normalize3f(n0);
        normalize3f(n1);
        cross_product3f(n0,n1,cp);
        
        theta = (float)asin(length3f(cp));

        normalize23f(cp,n2);        

        rotation_matrix3f(theta,n2[0], n2[1], n2[2], mat);
        
        multiply33f33f(mat,oss->system,oss->system);
        
        ObjectSliceInvalidate(I,cRepSlice,cRepAll,state);
        SceneDirty();
        
      }
      break;
    case cButModeMovFrag: /* move along "up" direction */
      {
        float up[3],v1[3];
        up[0]=oss->system[2];
        up[1]=oss->system[5];
        up[2]=oss->system[8];
        
        project3f(mov,up,v1);
        add3f(v1,oss->origin,oss->origin);
        ObjectSliceInvalidate(I,cRepSlice,cRepAll,state);        
        SceneDirty();
      }
      break;
    case cButModeTorFrag: 
      break;
    }
  }
}

int ObjectSliceGetVertex(ObjectSlice *I,int index,int base, float *v)
{
  int state=index-1;
  int offset=base-1;
  int result=false;

  ObjectSliceState *oss = NULL;

  if(state>=0) 
    if(state<I->NState) 
      if(I->State[state].Active)
        oss=I->State+state;

  if(oss) {
    if((offset>=0)&&(offset<oss->n_points))
      if(oss->flags[offset]) {
        copy3f(oss->points+3*offset,v);
        result=true;
      }
  }
  return(result);
}

int ObjectSliceGetOrigin(ObjectSlice *I,int state,float *origin)
{
  int ok=false;
  
  int cur_state = 0;
  ObjectSliceState *oss = NULL;

  if(state>=0) 
    if(state<I->NState) 
      if(I->State[state].Active)
        oss=I->State+state;

  while(1) {
    if(state<0) { /* all_states */
      oss = I->State + cur_state;
    } else {
      if(!oss) {
        if(I->NState&&
           ((SettingGet(cSetting_static_singletons)&&(I->NState==1))))
          oss=I->State;
      }
    }
    if(oss) {
      if(oss->Active) {	
        copy3f(oss->origin,origin); 
        ok=true;
      }
    }
    if(state>=0) break;
    cur_state = cur_state + 1;
    if(cur_state>=I->NState) break;
  }
  return ok;
}

static void ObjectSliceRender(ObjectSlice *I,int state,CRay *ray,Pickable **pick,int pass)
{

  int cur_state = 0;
  float alpha;
  
  ObjectSliceState *oss = NULL;

  if(SettingGet_b(NULL,I->Obj.Setting,cSetting_slice_track_camera)) {
    int update_flag = false;

    if(state>=0) 
      if(state<I->NState) 
        if(I->State[state].Active)
          oss=I->State+state;
    
    while(1) {
      if(state<0) { /* all_states */
        oss = I->State + cur_state;
      } else {
        if(oss) {

          SceneViewType view;
          float pos[3];
          
          SceneGetPos(pos);
          SceneGetView(view);
          
          if((diffsq3f(pos,oss->origin)>R_SMALL8) ||
             (diffsq3f(view,oss->system)>R_SMALL8) ||
             (diffsq3f(view+4,oss->system+3)>R_SMALL8) ||
             (diffsq3f(view+8,oss->system+6)>R_SMALL8))
            {
              copy3f(pos,oss->origin);
              
              copy3f(view,oss->system);
              copy3f(view+4,oss->system+3);
              copy3f(view+8,oss->system+6);
              oss->RefreshFlag=true;
              update_flag=true;
            }
        }
        if(state>=0) break;
        cur_state = cur_state + 1;
        if(cur_state>=I->NState) break;
      }
    }
    ObjectSliceUpdate(I);
  }

  ObjectPrepareContext(&I->Obj,ray);
  alpha = SettingGet_f(NULL,I->Obj.Setting,cSetting_transparency);
  alpha=1.0F-alpha;
  if(fabs(alpha-1.0)<R_SMALL4)
    alpha=1.0F;

  if(state>=0) 
    if(state<I->NState) 
      if(I->State[state].Active)
        oss=I->State+state;

  while(1) {
    if(state<0) { /* all_states */
      oss = I->State + cur_state;
    } else {
      if(!oss) {
        if(I->NState&&
           ((SettingGet(cSetting_static_singletons)&&(I->NState==1))))
          oss=I->State;
      }
    }
    if(oss) {
      if(oss->Active) {	
        if(ray){
          
          ray->fTransparentf(ray,1.0F-alpha);       
          if(I->Obj.RepVis[cRepSlice]) {
            float normal[3], *n0,*n1,*n2;
            int *strip = oss->strips;
            float *point = oss->points;
            float *color = oss->colors;
            int n=oss->n_strips;
            int a;
            int offset0,offset1,offset2,offset;
            int strip_active =false;
            int tri_count = 0;
            
            normal[0]=oss->system[2];
            normal[1]=oss->system[5];
            normal[2]=oss->system[8];

            n0=normal;
            n1=normal;
            n2=normal;

            for(a=0;a<n;a++) {
              offset = *(strip++);
              switch(offset) {
              case -1:
                if(!strip_active) {
                  glBegin(GL_TRIANGLES);
                }
                strip_active=true;
                tri_count = 0;
                break;
              case -2:
                if(strip_active)
                  glEnd();
                strip_active=false;
                break;
              default:
                if(strip_active) {
                  tri_count++;
                  offset2 = offset1;
                  offset1 = offset0;
                  offset0 = offset;
                  
                  if(tri_count>=3) {
                    
                    if(oss->normals) {
                      n0 = oss->normals + 3*offset0;
                      n1 = oss->normals + 3*offset1;
                      n2 = oss->normals + 3*offset2;
                    }
                      
                    if(tri_count&0x1) { /* get the handedness right ... */
                      ray->fTriangle3fv(ray,
                                        point+3*offset0,
                                        point+3*offset1,
                                        point+3*offset2,
                                        n0,n1,n2,
                                        color+3*offset0,
                                        color+3*offset1,
                                        color+3*offset2);
                    } else {
                      ray->fTriangle3fv(ray,
                                        point+3*offset1,
                                        point+3*offset0,
                                        point+3*offset2,
                                        n1,n0,n2,
                                        color+3*offset1,
                                        color+3*offset0,
                                        color+3*offset2);
                    }
                  }
                }
                break;
              }
            }
          }
          ray->fTransparentf(ray,0.0);
        } else if(pick&&PMGUI) {
          
          int i=(*pick)->index;
          int j;
          Pickable p;
          
          p.ptr = (void*)I;
          p.index = state+1;
          
          if(I->Obj.RepVis[cRepSlice]) {
            int *strip = oss->strips;
            float *point = oss->points;
            int n=oss->n_strips;
            int a;
            int offset0,offset1,offset2,offset;
            int strip_active =false;
            int tri_count = 0;
            for(a=0;a<n;a++) {
              offset = *(strip++);
              switch(offset) {
              case -1:
                if(!strip_active) {
                  glBegin(GL_TRIANGLES);
                }
                strip_active=true;
                tri_count = 0;
                break;
              case -2:
                if(strip_active)
                  glEnd();
                strip_active=false;
                break;
              default:
                if(strip_active) {
                  tri_count++;
                  offset2 = offset1;
                  offset1 = offset0;
                  offset0 = offset;
                  
                  if(tri_count>=3) {
                    
                    i++;
                    if(!(*pick)[0].ptr) {
                      /* pass 1 - low order bits */
                      glColor3ub((uchar)((i&0xF)<<4),(uchar)((i&0xF0)|0x8),(uchar)((i&0xF00)>>4)); 
                      VLACheck((*pick),Pickable,i);
                      (*pick)[i] = p; /* copy object and atom info */
                    } else { 
                      /* pass 2 - high order bits */
                      
                      j=i>>12;
                      glColor3ub((uchar)((j&0xF)<<4),(uchar)((j&0xF0)|0x8),(uchar)((j&0xF00)>>4)); 
                    }
                    p.bond = offset0 + 1;
                    
                    if(tri_count&0x1) { /* get the handedness right ... */
                      glVertex3fv(point+3*offset0);
                      glVertex3fv(point+3*offset1);
                      glVertex3fv(point+3*offset2);
                    } else {
                      glVertex3fv(point+3*offset1);
                      glVertex3fv(point+3*offset0);
                      glVertex3fv(point+3*offset2);
                    }
                  }
                }
                break;
              }
            }
            if(strip_active) { /* just in case */
              glEnd();
            }
          }
          (*pick)[0].index = i; /* pass the count */
        } else if(PMGUI) {

          int render_now = false;
          if(alpha>0.0001) {
            render_now = (pass==-1);
          } else 
            render_now = (!pass);

          if(render_now) {
            
            int use_dlst;
       
            SceneResetNormal(false);
            ObjectUseColor(&I->Obj);
            use_dlst = (int)SettingGet(cSetting_use_display_lists);
            if(use_dlst&&oss->displayList) {
              glCallList(oss->displayList);
            } else { 
              
              if(use_dlst) {
                if(!oss->displayList) {
                  oss->displayList = glGenLists(1);
                  if(oss->displayList) {
                    glNewList(oss->displayList,GL_COMPILE_AND_EXECUTE);
                  }
                }
              }

              if(I->Obj.RepVis[cRepSlice]) {
                int *strip = oss->strips;
                float *point = oss->points;
                float *color = oss->colors;
                float *col;
                float *vnormal = oss->normals;
                int n=oss->n_strips;
                int a;
                int offset;
                int strip_active =false;

              {
                float normal[3];
                normal[0]=oss->system[2];
                normal[1]=oss->system[5];
                normal[2]=oss->system[8];
                glNormal3fv(normal);
              }

                for(a=0;a<n;a++) {
                  offset = *(strip++);
                  switch(offset) {
                  case -1:
                    if(!strip_active) 
                      glBegin(GL_TRIANGLE_STRIP);
                    strip_active=true;
                    break;
                  case -2:
                    if(strip_active)
                      glEnd();
                    strip_active=false;
                    break;
                  default:
                    if(strip_active) {
                      col = color+3*offset;
                      if(vnormal)
                        glNormal3fv(vnormal+3*offset);
                      glColor4f(col[0],col[1],col[2],alpha);
                      glVertex3fv(point+3*offset);
                    }
                    break;
                  }
                }
                if(strip_active) /* just in case */
                  glEnd();
              }

            }
            
            if(use_dlst&&oss->displayList) {
              glEndList();
            }
            
          }
        }
      }
  }
  if(state>=0) break;
    cur_state = cur_state + 1;
    if(cur_state>=I->NState) break;
  }

}

/*========================================================================*/

static int ObjectSliceGetNStates(ObjectSlice *I) 
{
  return(I->NState);
}


ObjectSliceState *ObjectSliceStateGetActive(ObjectSlice *I,int state)
{
  ObjectSliceState *ms = NULL;
  if(state>=0) {
    if(state<I->NState) {
      ms=&I->State[state];
      if(!ms->Active)
        ms = NULL;
    }
  }
  return(ms);
}

/*========================================================================*/
ObjectSlice *ObjectSliceNew(void)
{
  OOAlloc(ObjectSlice);
  
  ObjectInit((CObject*)I);
  
  I->NState = 0;
  I->State=VLAMalloc(10,sizeof(ObjectSliceState),5,true); /* autozero important */

  I->Obj.type = cObjectSlice;
  
  I->Obj.fFree = (void (*)(struct CObject *))ObjectSliceFree;
  I->Obj.fUpdate =  (void (*)(struct CObject *)) ObjectSliceUpdate;
  I->Obj.fRender =(void (*)(struct CObject *, int, CRay *, Pickable **,int ))ObjectSliceRender;
  I->Obj.fInvalidate =(void (*)(struct CObject *,int,int,int))ObjectSliceInvalidate;
  I->Obj.fGetNFrame = (int (*)(struct CObject *)) ObjectSliceGetNStates;
  return(I);
}

/*========================================================================*/
void ObjectSliceStateInit(ObjectSliceState *oss)
{
  oss->displayList=0;
  oss->Active=true;
  oss->RefreshFlag=true;  
  oss->ExtentFlag=false;

  oss->values=NULL;
  oss->points=NULL;
  oss->flags=NULL;
  oss->colors=NULL;
  oss->strips=NULL;

  oss->n_points = 0;
  oss->n_strips = 0;

  UtilZeroMem(&oss->system,sizeof(float)*9); /* simple orthogonal coordinate system */
  oss->system[0] = 1.0F;
  oss->system[4] = 1.0F;
  oss->system[8] = 1.0F;

  zero3f(oss->origin); 

}

/*========================================================================*/
ObjectSlice *ObjectSliceFromMap(ObjectSlice * obj,ObjectMap *map,int state,int map_state)
{
  ObjectSlice *I;
  ObjectSliceState *oss;
  ObjectMapState *oms;



  if(!obj) {
    I=ObjectSliceNew();
  } else {
    I=obj;
  }

  if(state<0) state=I->NState;
  if(I->NState<=state) {
    VLACheck(I->State,ObjectSliceState,state);
    I->NState=state+1;
  }

  oss=I->State+state;

  ObjectSliceStateInit(oss);
  oss->MapState = map_state;
  oms = ObjectMapGetState(map,map_state);
  if(oms) {
    if(oss->points) {
      VLAFreeP(oss->points);
    }
    if(oss->values) {
      VLAFreeP(oss->points);
    }
    if(oss->flags) {
      VLAFreeP(oss->flags);
    }

    {
      float tmp[3];
      if(ObjectMapStateGetExcludedStats(oms,NULL,0.0F,0.0F,tmp)) {
        oss->MapMean = tmp[1];
        oss->MapStdev = tmp[2]-tmp[1];
      } else {
        oss->MapMean =0.0F;
        oss->MapStdev = 1.0F;
      }
    }
    /* simply copy the extents from the map -- not quite correct, but probably good enough */
    
    memcpy(oss->ExtentMin,oms->ExtentMin,3*sizeof(float));
    memcpy(oss->ExtentMax,oms->ExtentMax,3*sizeof(float));
  }

  strcpy(oss->MapName,map->Obj.Name);
  oss->ExtentFlag = true;

  /* set the origin of the slice to the center of the map */

  average3f(oss->ExtentMin,oss->ExtentMax,oss->origin);

  /* set the slice's system matrix to the current camera rotation matrix */

  {
    
    SceneViewType view;
     
    SceneGetView(view);
    copy3f(view,oss->system);
    copy3f(view+4,oss->system+3);
    copy3f(view+8,oss->system+6);
  }

  oss->RefreshFlag = true;

  if(I) {
    ObjectSliceRecomputeExtent(I);
  }

  I->Obj.ExtentFlag=true;

  SceneChanged();
  SceneCountFrames();
  return(I);
}

/*========================================================================*/

void ObjectSliceRecomputeExtent(ObjectSlice *I)
{
  int extent_flag = false;
  int a;
  ObjectSliceState *ms;

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

