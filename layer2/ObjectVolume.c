
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include<math.h>

#include"OOMac.h"
#include"ObjectVolume.h"
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
#include"Executive.h"
#include"PConv.h"
#include"P.h"
#include"Matrix.h"
#include"ObjectGadgetRamp.h"
#include"ShaderMgr.h"
#include"Field.h"


ObjectVolume *ObjectVolumeNew(PyMOLGlobals * G);

static void ObjectVolumeFree(ObjectVolume * I);
static void ObjectVolumeInvalidate(ObjectVolume * I, int rep, int level, int state);
void ObjectVolumeStateInit(PyMOLGlobals * G, ObjectVolumeState * vs);
void ObjectVolumeRecomputeExtent(ObjectVolume * I);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define HISTOGRAM_SIZE 64
  
#define volume_nColors 360

static ObjectVolumeState * ObjectVolumeGetActiveState(ObjectVolume * I) {
  int a;
  ok_assert(1, I);
  for(a = 0; a < I->NState; a++)
    if(I->State[a].Active)
      return I->State + a;
ok_except1:
  return NULL;
}

ObjectMapState * ObjectVolumeStateGetMapState(ObjectVolumeState * vs) {
  ObjectMapState *oms = NULL;
  ObjectMap *map = NULL;

  PyMOLGlobals * G = vs->State.G;

  map = ExecutiveFindObjectMapByName(G, vs->MapName);
  if(!map) {
    PRINTFB(G, FB_ObjectVolume, FB_Errors)
      "ObjectVolume-Error: map '%s' has been deleted.\n", vs->MapName
      ENDFB(G);
    return NULL;
  }

  return ObjectMapGetState(map, vs->MapState);
}

ObjectMapState * ObjectVolumeGetMapState(ObjectVolume * I) {
  ObjectVolumeState * ovs = ObjectVolumeGetActiveState(I);
  if(ovs)
    return ObjectVolumeStateGetMapState(ovs);
  return NULL;
}

#ifndef _PYMOL_NOPY
static PyObject *ObjectVolumeStateAsPyList(ObjectVolumeState * I)
{
  PyObject *result = NULL;
  result = PyList_New(19);
  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  PyList_SetItem(result, 1, PyString_FromString(I->MapName));
  PyList_SetItem(result, 2, PyInt_FromLong(I->MapState));
  PyList_SetItem(result, 3, PConvAutoNone(NULL) /* CrystalAsPyList(&I->Crystal) */);
  PyList_SetItem(result, 4, PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 7, PConvAutoNone(NULL) /* PConvIntArrayToPyList(I->Range, 6) */);
  PyList_SetItem(result, 8, PyFloat_FromDouble(0.0 /* I->Level */));
  PyList_SetItem(result, 9, PyFloat_FromDouble(0.0 /* I->Radius */));
  PyList_SetItem(result, 10, PyInt_FromLong(I->CarveFlag));
  PyList_SetItem(result, 11, PyFloat_FromDouble(I->CarveBuffer));
  if(I->CarveFlag && I->AtomVertex) {
    PyList_SetItem(result, 12, PConvFloatVLAToPyList(I->AtomVertex));
  } else {
    PyList_SetItem(result, 12, PConvAutoNone(NULL));
  }
  PyList_SetItem(result, 13, PyInt_FromLong(0 /* I->VolumeMode */));
  PyList_SetItem(result, 14, PyFloat_FromDouble(0.0 /* I->AltLevel */));
  PyList_SetItem(result, 15, PyInt_FromLong(1 /* I->quiet */));
  if(I->Field) {
    PyList_SetItem(result, 16, IsosurfAsPyList(I->Field));
  } else {
    PyList_SetItem(result, 16, PConvAutoNone(NULL));
  }
  PyList_SetItem(result,17,PyInt_FromLong(I->RampSize));
  if (I->Ramp) {
    PyList_SetItem(result, 18, PConvFloatArrayToPyList(I->Ramp, 5 * I->RampSize));
  } else {
    PyList_SetItem(result, 18, PConvAutoNone(NULL));
  }
  return (PConvAutoNone(result));
}
#endif

#ifndef _PYMOL_NOPY
static PyObject *ObjectVolumeAllStatesAsPyList(ObjectVolume * I)
{

  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NState);
  for(a = 0; a < I->NState; a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result, a, ObjectVolumeStateAsPyList(I->State + a));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(NULL));
    }
  }
  return (PConvAutoNone(result));
}
#endif

#ifndef _PYMOL_NOPY
static int ObjectVolumeStateFromPyList(PyMOLGlobals * G, ObjectVolumeState * I,
                                     PyObject * list)
{
  int ok = true;
  int ll = 0;
  PyObject *tmp;
  if(ok)
    ok = (list != NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active = false;
    else {
      ObjectVolumeStateInit(G, I);
      if(ok)
        ok = (list != NULL);
      if(ok)
        ok = PyList_Check(list);
      if(ok)
        ll = PyList_Size(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->Active);
      if(ok)
        ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->MapName, WordLength);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->MapState);
#if 0
      if(ok)
        ok = CrystalFromPyList(&I->Crystal, PyList_GetItem(list, 3));
#endif
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 4), &I->ExtentFlag);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 5), I->ExtentMin, 3);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6), I->ExtentMax, 3);
#if 0
      if(ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 7), I->Range, 6);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 8), &I->Level);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 9), &I->Radius);
#endif
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 10), &I->CarveFlag);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 11), &I->CarveBuffer);
      if(ok) {
        tmp = PyList_GetItem(list, 12);
        if(tmp == Py_None)
          I->AtomVertex = NULL;
        else
          ok = PConvPyListToFloatVLA(tmp, &I->AtomVertex);
      }
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 13), &I->VolumeMode);
      if(ok) {
        I->RefreshFlag = true;
        I->ResurfaceFlag = true;
      }
#if 0
      if(ok && (ll > 14)) {
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 14), &I->AltLevel);
      } else {
        I->AltLevel = I->Level;
      }
      if(ok && (ll > 15)) {
        ok = PConvPyIntToInt(PyList_GetItem(list, 15), &I->quiet);
      } else {
        I->quiet = true;
      }
#endif
      if(ok && (ll > 16)) {
        tmp = PyList_GetItem(list, 16);
        if(tmp == Py_None)
          I->Field = NULL;
        else
          ok = ((I->Field = IsosurfNewFromPyList(G, tmp)) != NULL);
      }
      if (ok && (ll > 17)) {
        ok = PConvPyIntToInt(PyList_GetItem(list, 17), &I->RampSize);
      }
      if (ok && (ll > 18)) {
        tmp = PyList_GetItem(list, 18);
        if(tmp == Py_None)
          I->Ramp = NULL;
        else
          ok = PConvPyListToFloatArray(tmp, &I->Ramp);
      }
    }
  }
  return (ok);
}
#endif

#ifndef _PYMOL_NOPY
static int ObjectVolumeAllStatesFromPyList(ObjectVolume * I, PyObject * list)
{

  int ok = true;
  int a;
  VLACheck(I->State, ObjectVolumeState, I->NState);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    for(a = 0; a < I->NState; a++) {
      ok = ObjectVolumeStateFromPyList(I->Obj.G, I->State + a, PyList_GetItem(list, a));
      if(!ok)
        break;
    }
  }
  return (ok);
}
#endif

int ObjectVolumeNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectVolume ** result)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  int ll;
  ObjectVolume *I = NULL;
  (*result) = NULL;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  I = ObjectVolumeNew(G);
  if(ok)
    ok = (I != NULL);
  if(ok)
    ok = ObjectFromPyList(G, PyList_GetItem(list, 0), &I->Obj);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->NState);
  if(ok)
    ok = ObjectVolumeAllStatesFromPyList(I, PyList_GetItem(list, 2));
  if(ok) {
    (*result) = I;
    ObjectVolumeRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
#endif
}

PyObject *ObjectVolumeAsPyList(ObjectVolume * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(&I->Obj));
  PyList_SetItem(result, 1, PyInt_FromLong(I->NState));
  PyList_SetItem(result, 2, ObjectVolumeAllStatesAsPyList(I));
  return (PConvAutoNone(result));
#endif
}

static void ObjectVolumeStateFree(ObjectVolumeState * vs)
{
    int t;
  ObjectStatePurge(&vs->State);
  if(vs->State.G->HaveGUI) {
    for (t=0; t<2; t++) {
      if (vs->textures[t]) {
        glDeleteTextures(1, (const GLuint *) &vs->textures[t]);
	vs->textures[t] = 0;
	/*
        if(PIsGlutThread()) {
          if(vs->State.G->ValidContext) {
            glDeleteTextures(1, (const GLuint *) &vs->textures[t]);
            vs->textures[t] = 0;
          }
        } else {
          char buffer[255];       // pass this off to the main thread
          sprintf(buffer, "_cmd.gl_delete_texture(cmd._COb,%d)\n", vs->textures[t]);
          PParse(vs->State.G, buffer);
          vs->textures[t] = 0;
        }
      */
      }
    }
  }
  if(vs->Field) {
    IsosurfFieldFree(vs->State.G, vs->Field);
    vs->Field = NULL;
  }
  if(vs->volume) {
    FieldFree(vs->volume);
  }
  VLAFreeP(vs->AtomVertex);
  if (vs->Ramp)
    FreeP(vs->Ramp);
  vs->Active = false;
}

static void ObjectVolumeFree(ObjectVolume * I)
{
  int a,i;
  for(a = 0; a < I->NState; a++) {
    for (i=0; i<2; i++){
      if (I->State[a].textures[i]){
       glDeleteTextures(1, (const GLuint *) &I->State[a].textures[i]);
       I->State[a].textures[i] = 0;
      }
    }
    if(I->State[a].Active)
      ObjectVolumeStateFree(I->State + a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);

  OOFreeP(I);
}

int ObjectVolumeInvalidateMapName(ObjectVolume * I, char *name)
{
  int a;
  ObjectVolumeState *vs;
  int result = false;
  for(a = 0; a < I->NState; a++) {
    vs = I->State + a;
    if(vs->Active) {
      if(strcmp(vs->MapName, name) == 0) {
        ObjectVolumeInvalidate(I, cRepAll, cRepInvAll, a);
        result = true;
      }
    }
  }
  return result;
}

static void ObjectVolumeInvalidate(ObjectVolume * I, int rep, int level, int state)
{
  int a;
  int once_flag = true;
  if(level >= cRepInvExtents) {
    I->Obj.ExtentFlag = false;
  }

  PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeInvalidate-Msg: %d states.\n", I->NState
    ENDFB(I->Obj.G);

  if((rep == cRepVolume) || (rep == cRepAll) || (rep == cRepExtent)) {
    for(a = 0; a < I->NState; a++) {
      if(state < 0)
        once_flag = false;
      if(!once_flag)
        state = a;
      if(level == cRepInvColor || level == cRepInvAll) {
        I->State[state].RecolorFlag = true;
      }
      if(level != cRepInvColor) {
        I->State[state].ResurfaceFlag = true;
        I->State[state].RefreshFlag = true;
      }
      SceneChanged(I->Obj.G);
      if(once_flag)
        break;
    }
  }
}

static void ObjectVolumeUpdate(ObjectVolume * I)
{
  int a;
  ObjectVolumeState *vs;
  ObjectMapState *oms = NULL;
  ObjectMap *map = NULL;
  float carve_buffer;
  int avoid_flag = false;
  int flag;
  int h, k, l;
  int i, j;
  float range;
  MapType *voxelmap;            /* this has nothing to do with isosurfaces... */
  PyMOLGlobals * G = I->Obj.G;

  for(a = 0; a < I->NState; a++) {
    vs = I->State + a;
    if(!vs || !vs->Active)
      continue;

    PRINTFD(G, FB_ObjectVolume)
      "ObjectVolumeUpdate: state=%d, refresh=%d, resurface=%d.\n",
      a, vs->RefreshFlag, vs->ResurfaceFlag ENDFD(G);

    oms = ObjectVolumeStateGetMapState(vs);

    if(!oms) {
      vs->ResurfaceFlag = false;
      continue;
    }

    if(vs->RefreshFlag || vs->ResurfaceFlag) {
      if(oms->State.Matrix) {
        ObjectStateSetMatrix(&vs->State, oms->State.Matrix);
      } else if(vs->State.Matrix) {
        ObjectStateResetMatrix(&vs->State);
      }
    }

    // handle legacy or default color ramp
    if(!vs->Ramp || vs->RampSize && vs->Ramp[0] == 0.f
        && vs->Ramp[5 * (vs->RampSize - 1)] == 359.f) {
      float min_max_mean_stdev[4];

      // data min/max
      range =  SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_data_range);
      ObjectMapStateGetHistogram(I->Obj.G, oms, 0, range, min_max_mean_stdev, 0.f, 0.f);

      if(vs->Ramp) {
        // legacy color ramp (0..359)
        range = min_max_mean_stdev[1] - min_max_mean_stdev[0];
        PRINTFB(G, FB_ObjectVolume, FB_Warnings)
          " ObjectVolumeUpdate: detected legacy color ramp\n" ENDFB(G);
        for (i = 0; i < vs->RampSize * 5; i += 5) {
          vs->Ramp[i] = vs->Ramp[i] / 359.f * range + min_max_mean_stdev[0];
        }
      } else {
        // default color ramp (1.0 sigma peak)
        if(!vs->Ramp) {
          float defaultramp[] = {
            min_max_mean_stdev[2] + 0.7 * min_max_mean_stdev[3],
            0., 0., 1., 0.0,
            min_max_mean_stdev[2] + 1.0 * min_max_mean_stdev[3],
            0., 1., 1., 0.2,
            min_max_mean_stdev[2] + 1.3 * min_max_mean_stdev[3],
            0., 0., 1., 0.0
          };
          vs->RecolorFlag = true;
          vs->RampSize = 3;
          vs->Ramp = Alloc(float, 5 * vs->RampSize);
          memcpy(vs->Ramp, defaultramp, 5 * vs->RampSize * sizeof(float));
        }
      }
    }

    if(I->Obj.RepVis[cRepVolume] && vs->ResurfaceFlag) {
      Isofield *field = NULL;
      vs->ResurfaceFlag = false;
      if(vs->Field) {
        field = vs->Field;
      } else if(oms->Field) {
        field = oms->Field;
      } else {
        field = NULL;
      }

      if(field) {
        float *min_ext, *max_ext;
        float tmp_min[3], tmp_max[3];
        if(MatrixInvTransformExtentsR44d3f(vs->State.Matrix,
              vs->ExtentMin, vs->ExtentMax,
              tmp_min, tmp_max)) {
          min_ext = tmp_min;
          max_ext = tmp_max;
        } else {
          min_ext = vs->ExtentMin;
          max_ext = vs->ExtentMax;
        }

        if (vs->volume)
          FieldFree(vs->volume);

        vs->volume = FieldNewCopy(I->Obj.G, field->data);

        IsofieldGetCorners(G, field, vs->Corner);

        // transform corners by state matrix
        if(vs->State.Matrix) {
          for(i = 0; i < 8; i++)
            transform44d3f(vs->State.Matrix,
                vs->Corner + 3 * i,
                vs->Corner + 3 * i);
        }
      }

      if(vs->volume && vs->CarveFlag && vs->AtomVertex) {

        carve_buffer = vs->CarveBuffer;
        if(vs->CarveBuffer < 0.0F) {
          avoid_flag = true;
          carve_buffer = -carve_buffer;
        }

        // cull my friend, cull */ 
        voxelmap = MapNew(I->Obj.G,
            -carve_buffer, vs->AtomVertex,
            VLAGetSize(vs->AtomVertex) / 3, NULL);
        if(voxelmap) {

          int x, y, z;
          int dx, dy, dz;
          float fx, fy, fz;
          float vv[3];
          float *fdata = (float*)vs->volume->data;
          float min_val = fdata[0];
          float idx, idy, idz;

          float vx[3] = {vs->Corner[3]-vs->Corner[0],
            vs->Corner[4]-vs->Corner[1],
            vs->Corner[5]-vs->Corner[2]};
          float vy[3] = {vs->Corner[6]-vs->Corner[0],
            vs->Corner[7]-vs->Corner[1],
            vs->Corner[8]-vs->Corner[2]};
          float vz[3] = {vs->Corner[12]-vs->Corner[0],
            vs->Corner[13]-vs->Corner[1],
            vs->Corner[14]-vs->Corner[2]};

          MapSetupExpress(voxelmap);


          dx = vs->volume->dim[0];
          dy = vs->volume->dim[1];
          dz = vs->volume->dim[2];
          for (x=0;x<dx*dy*dz;x++) {
            if (fdata[x] < min_val) min_val = fdata[x];
          }
          idx = 1.0 / dx;
          idy = 1.0 / dy;
          idz = 1.0 / dz;
          for (z=0; z<dz; z++) {
            fz = z * idz;
            for (y=0; y<dy; y++) {
              fy = y * idy;
              for (x=0; x<dx; x++) {
                fx = x * idx;
                vv[0] = vs->Corner[0] + fx * vx[0] + fy * vy[0] + fz * vz[0];
                vv[1] = vs->Corner[1] + fx * vx[1] + fy * vy[1] + fz * vz[1];
                vv[2] = vs->Corner[2] + fx * vx[2] + fy * vy[2] + fz * vz[2];
                flag = false;
                MapLocus(voxelmap, vv, &h, &k, &l);
                i = *(MapEStart(voxelmap, h, k, l));
                if(i) {
                  j = voxelmap->EList[i++];
                  while(j >= 0) {
                    if(within3f(vs->AtomVertex + 3 * j, vv, carve_buffer)) {
                      flag = true;
                      break;
                    }
                    j = voxelmap->EList[i++];
                  }
                }
                if(avoid_flag)
                  flag = !flag;
                if(!flag) { 
                  *Ffloat3p(vs->volume, x, y, z) = min_val;
                }
              }
            }
          }
          MapFree(voxelmap);
        }
      }
    }
    vs->isUpdated = true;
    SceneInvalidate(I->Obj.G);
  }
  if(!I->Obj.ExtentFlag) {
    ObjectVolumeRecomputeExtent(I);
    if(I->Obj.ExtentFlag)
      SceneInvalidate(I->Obj.G);
  }
}

int ObjectVolumeAddSlicePoint(float *p0, float *p1, float *zaxis, float d, float *slice, float *t0, float *t1, float *tex_coords, float *origin);
void ObjectVolumeDrawSlice(float *points, float *tex_coords, int n_points, float *zaxis);

/*
 * Converting Ramp to `volume_nColors * 4` sized interpolated RGBA color
 * array. Returns allocated memory.
 * Assigns data minimum and range covered by ramp to `ramp_min` and `ramp_range`.
 */
float * ObjectVolumeStateGetColors(PyMOLGlobals * G, ObjectVolumeState * ovs,
    float *ramp_min, float *ramp_range) {
  int i, j, k;
  int lowerId, upperId;
  float mixc, mixcincr, r_min, range, binsize;
  float * colors;

  ok_assert(1, ovs->Ramp && ovs->RampSize > 1);

  r_min = ovs->Ramp[0];
  range = ovs->Ramp[5 * (ovs->RampSize - 1)] - r_min;

  ok_assert(1, range > R_SMALL4);

  binsize = range / volume_nColors;
  r_min -= binsize * 2;
  range += binsize * 4;

  colors = Calloc(float, 4 * volume_nColors);
  ok_assert(1, colors);

  for (i = 0; i < ovs->RampSize; i++) {
    lowerId = upperId;
    upperId = (int) (volume_nColors * (ovs->Ramp[i * 5] - r_min) / range);

    if(i == 0)
      continue;

    mixcincr = 1.f / (upperId - lowerId);

    for (j = lowerId, mixc = 1.f; j < upperId; j++, mixc -= mixcincr){
      if(j < 0 || j >= volume_nColors)
        continue;

      for (k = 0; k < 4; k++)
        colors[j * 4 + k] = ovs->Ramp[i * 5 - 4 + k] * mixc +
                            ovs->Ramp[i * 5 + 1 + k] * (1.f - mixc);
    }
  }

  *ramp_min = r_min;
  *ramp_range = range;
  return colors;
ok_except1:
  PRINTFB(G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeStateGetColors failed\n" ENDFB(G);
  return NULL;
}

/*
 * Render bounding box
 * TODO: duplicate in other reps?
 */
static void ExtentRender(float * corner) {
  int i, j, ci[] = {
     0,  3,  3,  9,  9,  6,  6,  0,
    12, 15, 15, 21, 21, 18, 18, 12,
     0, 12,  3, 15,  9, 21,  6, 18
  };
#ifndef PURE_OPENGL_ES_2
#ifdef _PYMOL_GL_DRAWARRAYS
  GLfloat polyVerts[8 * 3 * 3];
  for(i = 0; i < 8 * 3; i++)
    for(j = 0; j < 3; j++)
      polyVerts[i * 3 + j] = corner[ci[i] + j];
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, polyVerts);
  glDrawArrays(GL_LINES, 0, 24);
  glDisableClientState(GL_VERTEX_ARRAY);
#else
  glBegin(GL_LINES);
  for(i = 0; i < 8 * 3; i++)
    glVertex3fv(corner + ci[i]);
  glEnd();
#endif
#endif
}

static void ObjectVolumeRender(ObjectVolume * I, RenderInfo * info)
{
#ifndef PURE_OPENGL_ES_2
  PyMOLGlobals *G = I->Obj.G;
  int state = info->state;
  CRay *ray = info->ray;
  int pass = info->pass;
  int *n = NULL;
  int a = 0;
  ObjectVolumeState *vs = NULL;
  float volume_layers =  SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_layers);
  /* make this a setting? */
  GLint alpha_func;
  GLfloat alpha_ref;
  float tex_corner[24];
  float *corner, *m;
  float zaxis[3];
  float points[36], tex_coords[36];
  int n_points;
  float d, sliceRange, sliceDelta;
  float origin[3];
  CShaderPrg *shaderPrg;

  if(info->pick || pass != -1)
    return;

  if(!G->HaveGUI || !G->ValidContext)
    return;

  /* bail if no shaders */
  if (G && !(CShaderMgr_ShadersPresent(G->ShaderMgr)))
      return;

  // ViewElem/TTT Matrix
  ObjectPrepareContext(&I->Obj, ray);

  for(a = 0; a < I->NState; ++a) {

    if(state < 0 || state == a) {
      vs = I->State + a;
    } else if(a == 0 && I->NState == 1 && SettingGetGlobal_b(G, cSetting_static_singletons)) {
      vs = I->State;
    } else {
      continue;
    }

    if(!vs || !vs->Active)
      continue;

    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
      "ObjectVolumeRender-Msg: state=%d, pass=%d, refresh=%d, recolor=%d.\n",
      a, pass, vs->RefreshFlag, vs->RecolorFlag ENDFB(I->Obj.G); 

    corner = vs->Corner;

    SceneResetNormal(I->Obj.G, false);

    // render bounding box
    if(I->Obj.RepVis[cRepExtent]) {
      if(!info->line_lighting)
        glDisable(GL_LIGHTING);
      ObjectUseColor(&I->Obj);
      ExtentRender(corner);
    }

    // upload color ramp texture
    if (vs->RecolorFlag) {
      float * colors = ObjectVolumeStateGetColors(G, vs, &vs->ramp_min, &vs->ramp_range);
      if(!colors)
        continue;

      if (vs->textures[1]) {
        glDeleteTextures(1, (const GLuint *) &vs->textures[1]);
        vs->textures[1] = 0;
      }

      glGenTextures(1, (GLuint *) &vs->textures[1]);
      glBindTexture(GL_TEXTURE_1D, vs->textures[1]);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
      glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, volume_nColors, 0, GL_RGBA, GL_FLOAT, colors);

      mfree(colors);
      vs->RecolorFlag = false;
    }

    // upload map data texture
    if (!vs->textures[0] || vs->RefreshFlag) {
      int volume_bit_depth;

      volume_bit_depth = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_volume_bit_depth);
      volume_bit_depth = (volume_bit_depth < 17) ? GL_R16F : GL_R32F;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
      glTexImage3D = getTexImage3D(); 
      glActiveTexture = getActiveTexture();
      if (! glActiveTexture || ! glTexImage3D) {
        PRINTFB(G, FB_ObjectVolume, FB_Errors)
          " ObjectVolumeRender-Error: Could not bind the glActiveTexture or glTexImage3D function.\n"
          ENDFB(G);
        return;
      } 
#endif
/* END PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */

      if (vs->textures[0]) {
        glDeleteTextures(1, (const GLuint *) &vs->textures[0]);
        vs->textures[0] = 0;
      }
      // Create a 3D texture
      glGenTextures(1, (GLuint *) &vs->textures[0]);
      glBindTexture(GL_TEXTURE_3D, vs->textures[0]);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

      glTexImage3D(GL_TEXTURE_3D, 0, volume_bit_depth,
          vs->volume->dim[2], vs->volume->dim[1], vs->volume->dim[0], 0, 
          GL_RED, GL_FLOAT, vs->volume->data);

      vs->RefreshFlag = false;
    }

    // render volume
    if(I->Obj.RepVis[cRepVolume]) {
      int i, j;

      glDisable(GL_LIGHTING);

      // half grid cell inset in texture corners (texture coordinate units)
      for(j = 0; j < 3; j++) {
        float offset = 0.5 / vs->volume->dim[2 - j];
        int bit = 1 << (2 - j);
        for(i = 0; i < 8; i++)
          tex_corner[i * 3 + j] = (i & bit) ? 1.0 - offset : offset;
      }

      // for z-axis
      m = SceneGetMatrix(G);

      for(j = 0; j < 3; j++) {
        // map center
        origin[j] = corner[j] + 0.5 * (corner[21 + j] - corner[j]);
        // z-axis
        zaxis[j] = m[j * 4 + 2];
      }

      // determine number of slices based on max extent
      // and slice option
      sliceRange = 0.5*sqrt(2.0) * 
        fmax(fmax(fabs(corner[21]-corner[0]), fabs(corner[22]-corner[1])), 
            fabs(corner[23]-corner[2]));
      sliceDelta = (sliceRange / volume_layers);

      // load shader
      shaderPrg = CShaderMgr_GetShaderPrg(G->ShaderMgr, "volume");
      CShaderPrg_Enable(shaderPrg);
      CShaderPrg_Set1i(shaderPrg, "volumeTex", 0);
      CShaderPrg_Set1i(shaderPrg, "colorTex", 1);
      CShaderPrg_Set1f(shaderPrg, "volumeScale", 1.0 / vs->ramp_range);
      CShaderPrg_Set1f(shaderPrg, "volumeBias", (-vs->ramp_min) / vs->ramp_range);

      // background and fog stuff
      {
        float fog[4];
        int bg_gradient = SettingGet_b(G, NULL, NULL, cSetting_bg_gradient);
        const char * bg_image_filename = SettingGet_s(G, NULL, NULL, cSetting_bg_image_filename);

        CShaderPrg_Set1f(shaderPrg, "fogIsSolidColor", bg_gradient || bg_image_filename && bg_image_filename[0] ? 0.f : 1.f);
        CShaderPrg_Set3fv(shaderPrg, "fogSolidColor", ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb)));
        CShaderPrg_SetFogUniforms(G, shaderPrg);
        CShaderPrg_Set1f(shaderPrg, "fog_enabled", SettingGetGlobal_b(G, cSetting_depth_cue) ? 1.f : 0.f);

        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, OrthoGetBackgroundTextureID(G));
        if (!(shaderPrg->uniform_set & 4)){
          CShaderPrg_Set1i(shaderPrg, "bgTextureMap", 4);
          shaderPrg->uniform_set |= 4;
        }

        SceneSetFog(G, fog);
      }

      // bind color ramp and map data textures
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_1D, vs->textures[1]);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_3D, vs->textures[0]);

      // alpha: everything passes
      // Not sure if we really need to restore this
      glGetIntegerv(GL_ALPHA_TEST_FUNC, &alpha_func);
      glGetFloatv(GL_ALPHA_TEST_REF, &alpha_ref);
      glAlphaFunc(GL_ALWAYS, 0.0);

      // This is setting used for PyMOL, but just to be on a safe side
      // we set glBlendFunct explicitely here
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      // draw slices
      {
        int i, cornerindices[] = {
          0,   3,    3,  9,    9,  6,    6,  0,
          12, 15,   15, 21,   21, 18,   18, 12,
          0,  12,    3, 15,    9, 21,    6, 18
        };
        for (d=sliceRange; d>=-sliceRange; d -= sliceDelta) {
          // Slice the volume
          n_points = 0;
          for(i = 0; i < 24; i += 2) {
            int j = cornerindices[i], k = cornerindices[i + 1];
            n_points += ObjectVolumeAddSlicePoint(
                corner + j, corner + k, zaxis, d, points + n_points,
                tex_corner + j, tex_corner + k, tex_coords + n_points, origin);
          }
          ObjectVolumeDrawSlice(points, tex_coords, n_points/3, zaxis);
        }
      }

      CShaderPrg_Disable(shaderPrg);

      // restore
      glAlphaFunc(alpha_func, alpha_ref);
    }

    glEnable(GL_LIGHTING);
  }
#endif
}

void ObjectVolumeDrawSlice(float *points, float *tex_coords, int n_points, float *zaxis)
{
  float center[3], v[3], w[3], q[3];
  float angles[12];
  float a, c, s;
  int vertices[12];
  int i, j;

    if (!n_points) return;

    // Calculate the polygon center
    zero3f(center);

    for (i=0; i<3*n_points; i+=3) {
      add3f(center, points + i, center);      // center += (points + i)
    }    

    scale3f(center, 1. / n_points, center);   // center /= n_points
    subtract3f(points, center, v);            // v = points - center
    normalize3f(v);
    
    // Sort vertices by rotation angle around the central axis
    for (i=0; i<n_points; i++) {
      subtract3f(points + 3 * i, center, w);  // w = (points + 3 * i) - center
      normalize3f(w);
      cross_product3f(v, w, q);
      c = dot_product3f(v, w);
      s = dot_product3f(zaxis, q);
      a = atan2(s, c);
      if (a < 0.0f) a += 2.0f * M_PI;
      j = i-1;
      while (j>=0 && angles[j]>a) {
        angles[j+1] = angles[j];
        vertices[j+1] = vertices[j];
        j--;
      }
      angles[j+1] = a;
      vertices[j+1] = i;
    }

    // Now the vertices are sorted so draw the polygon
#ifdef _PYMOL_GL_DRAWARRAYS
    {
      int nverts = n_points, pl = 0;
      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
      ALLOCATE_ARRAY(GLfloat,texVals,nverts*3)
      float *tmp_ptr;
      for (i=0; i<n_points; i++) {
	tmp_ptr = &tex_coords[3*vertices[(i)%n_points]];
	texVals[pl] = tmp_ptr[0]; texVals[pl+1] = tmp_ptr[1]; texVals[pl+2] = tmp_ptr[2];	
	tmp_ptr = &points[3*vertices[(i)%n_points]];
	vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];	
	pl += 3;
      }
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_TEXTURE_COORD_ARRAY);
      glVertexPointer(3, GL_FLOAT, 0, vertVals);
      glTexCoordPointer(3, GL_FLOAT, 0, texVals);
      glDrawArrays(GL_TRIANGLE_FAN, 0, nverts);
      glDisableClientState(GL_VERTEX_ARRAY);
      glDisableClientState(GL_TEXTURE_COORD_ARRAY);
      DEALLOCATE_ARRAY(vertVals)
      DEALLOCATE_ARRAY(texVals)
    }
#else
    glBegin(GL_POLYGON);
    for (i=0; i<n_points; i++) {
      glTexCoord3fv(&tex_coords[3 * vertices[i]]);
      glVertex3fv(&points[3 * vertices[i]]);
    }
    glEnd();
#endif
}

int ObjectVolumeAddSlicePoint(float *pt0, float *pt1, float *zaxis, float d, 
                              float *coords, float *t0, float *t1, float *tex_coords, float *origin)
{
    
  float p0[3];
  float p1[3];
  float u;

    p0[0] = pt0[0] - origin[0];
    p0[1] = pt0[1] - origin[1];
    p0[2] = pt0[2] - origin[2];
    p1[0] = pt1[0] - origin[0];
    p1[1] = pt1[1] - origin[1];
    p1[2] = pt1[2] - origin[2];

    u = (zaxis[0]*p0[0] + zaxis[1]*p0[1] + zaxis[2]*p0[2] + d) /
      (zaxis[0]*(p0[0]-p1[0]) + zaxis[1]*(p0[1]-p1[1]) + zaxis[2]*(p0[2]-p1[2]));

    if (u>=0.0F && u<=1.0F) {
        coords[0] = pt0[0] + (pt1[0]-pt0[0])*u;
        coords[1] = pt0[1] + (pt1[1]-pt0[1])*u;
        coords[2] = pt0[2] + (pt1[2]-pt0[2])*u;
        tex_coords[0] = t0[0] + (t1[0]-t0[0])*u;
        tex_coords[1] = t0[1] + (t1[1]-t0[1])*u;
        tex_coords[2] = t0[2] + (t1[2]-t0[2])*u;
        return 3;
    }
    return 0;
}

/*========================================================================*/

static int ObjectVolumeGetNStates(ObjectVolume * I)
{
  return (I->NState);
}


/*========================================================================*/
ObjectVolume *ObjectVolumeNew(PyMOLGlobals * G)
{
  OOAlloc(G, ObjectVolume);

  ObjectInit(G, (CObject *) I);

  I->NState = 0;
  I->State = VLACalloc(ObjectVolumeState, 10);   /* autozero important */

  I->Obj.type = cObjectVolume;

  I->Obj.fFree = (void (*)(CObject *)) ObjectVolumeFree;
  I->Obj.fUpdate = (void (*)(CObject *)) ObjectVolumeUpdate;
  I->Obj.fRender = (void (*)(CObject *, RenderInfo *)) ObjectVolumeRender;
  I->Obj.fInvalidate = (void (*)(CObject *, int, int, int)) ObjectVolumeInvalidate;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectVolumeGetNStates;
  return (I);
}


/*========================================================================*/
void ObjectVolumeStateInit(PyMOLGlobals * G, ObjectVolumeState * vs)
{
  if(vs->Active)
    ObjectStatePurge(&vs->State);
  if(vs->Field) {
    IsosurfFieldFree(vs->State.G, vs->Field);
    vs->Field = NULL;
  }
  ObjectStateInit(G, &vs->State);
  if(vs->AtomVertex) {
    VLAFreeP(vs->AtomVertex);
  }
  vs->Active = true;
  vs->ResurfaceFlag = true;
  vs->RecolorFlag = true;
  vs->ExtentFlag = false;
  vs->CarveFlag = false;
  vs->CarveBuffer = 0.0;
  vs->AtomVertex = NULL;
  vs->caption[0] = 0;
  vs->Field = NULL;
  vs->volume = NULL;
  vs->textures[0] = 0;
  vs->textures[1] = 0;
  vs->isUpdated = false;
  // Initial ramp
  vs->RampSize = 0;
  vs->Ramp = NULL;
}


/*========================================================================*/
ObjectVolume *ObjectVolumeFromXtalSym(PyMOLGlobals * G, ObjectVolume * obj, ObjectMap * map,
                                  CSymmetry * sym,
                                  int map_state,
                                  int state, float *mn, float *mx,
                                  float level, int meshMode,
                                  float carve, float *vert_vla,
                                  float alt_level, int quiet)
{
  int ok = true;
  ObjectVolume *I;
  ObjectVolumeState *vs;
  ObjectMapState *oms;
  int created = !obj;

  if(created) {
    I = ObjectVolumeNew(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectVolumeState, state);
    I->NState = state + 1;
  }

  vs = I->State + state;
  ObjectVolumeStateInit(G, vs);

  strcpy(vs->MapName, map->Obj.Name);
  vs->MapState = map_state;
  oms = ObjectMapGetState(map, map_state);

#if 0
  vs->VolumeMode = meshMode;
#endif
  if(oms) {

    copy3f(mn, vs->ExtentMin);  /* this is not exactly correct...should actually take vertex points from range */
    copy3f(mx, vs->ExtentMax);

    if(oms->State.Matrix) {
      ObjectStateSetMatrix(&vs->State, oms->State.Matrix);
    } else if(vs->State.Matrix) {
      ObjectStateResetMatrix(&vs->State);
    }

    {
      float *min_ext, *max_ext;
      float tmp_min[3], tmp_max[3];
      if(MatrixInvTransformExtentsR44d3f(vs->State.Matrix,
                                         vs->ExtentMin, vs->ExtentMax,
                                         tmp_min, tmp_max)) {
        min_ext = tmp_min;
        max_ext = tmp_max;
      } else {
        min_ext = vs->ExtentMin;
        max_ext = vs->ExtentMax;
      }

      if(sym) {
        int eff_range[6];

        if(IsosurfGetRange
           (G, oms->Field, oms->Symmetry->Crystal, min_ext, max_ext, eff_range, false)) {
          int fdim[3];
          int expand_result;
          /* need to generate symmetry-expanded temporary map */

          fdim[0] = eff_range[3] - eff_range[0];
          fdim[1] = eff_range[4] - eff_range[1];
          fdim[2] = eff_range[5] - eff_range[2];
          vs->Field = IsosurfFieldAlloc(I->Obj.G, fdim);

          expand_result =
            IsosurfExpand(oms->Field, vs->Field, oms->Symmetry->Crystal, sym, eff_range);

          if(expand_result == 0) {
            ok = false;
            if(!quiet) {
              PRINTFB(G, FB_ObjectVolume, FB_Warnings)
                " ObjectVolume-Warning: no symmetry expanded map points found.\n" ENDFB(G);
            }
          } else {
            if(!quiet) {
              PRINTFB(G, FB_ObjectVolume, FB_Warnings)
                " ObjectVolume-Warning: not all symmetry expanded points covered by map.\n"
                ENDFB(G);
            }
          }
        }
      }
    }
    vs->ExtentFlag = true;
  }
  if(ok) {
    if(carve != 0.0) {
      vs->CarveFlag = true;
      vs->CarveBuffer = carve;
      vs->AtomVertex = vert_vla;
    }
    I->Obj.ExtentFlag = false;
  }
  if(!ok && created) {
    ObjectVolumeStateFree(vs);
    I = NULL;
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


/*========================================================================*/
ObjectVolume *ObjectVolumeFromBox(PyMOLGlobals * G, ObjectVolume * obj, ObjectMap * map,
                              int map_state,
                              int state, float *mn, float *mx,
                              float level, int meshMode,
                              float carve, float *vert_vla, float alt_level, int quiet)
{
  return ObjectVolumeFromXtalSym(G, obj, map, NULL, map_state, state, mn, mx,
      level, meshMode, carve, vert_vla, alt_level, quiet);
}


/*========================================================================*/

void ObjectVolumeRecomputeExtent(ObjectVolume * I)
{
  int extent_flag = false;
  int a;
  ObjectVolumeState *vs;

  for(a = 0; a < I->NState; a++) {
    vs = I->State + a;
    if(vs->Active) {
      if(vs->ExtentFlag) {
        if(!extent_flag) {
          extent_flag = true;
          copy3f(vs->ExtentMax, I->Obj.ExtentMax);
          copy3f(vs->ExtentMin, I->Obj.ExtentMin);
        } else {
          max3f(vs->ExtentMax, I->Obj.ExtentMax, I->Obj.ExtentMax);
          min3f(vs->ExtentMin, I->Obj.ExtentMin, I->Obj.ExtentMin);
        }
      }
    }
  }

  I->Obj.ExtentFlag = extent_flag;

  if(I->Obj.TTTFlag && I->Obj.ExtentFlag) {
    float *ttt;
    double tttd[16];
    if(ObjectGetTTT(&I->Obj, &ttt, -1)) {
      convertTTTfR44d(ttt, tttd);
      MatrixTransformExtentsR44d3f(tttd,
                                   I->Obj.ExtentMin, I->Obj.ExtentMax,
                                   I->Obj.ExtentMin, I->Obj.ExtentMax);
    }
  }
}


/*==============================================================================*/
PyObject * ObjectVolumeGetField(ObjectVolume * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* TODO: Allow for multi-state maps? */
  CField* F;
  int n_elem;
  PyObject * result = NULL;
  ObjectVolumeState *ovs;

  if(I && (ovs = ObjectVolumeGetActiveState(I))) {
    F = ovs->volume;
    n_elem = F->size / F->base_size;
    result = PConvFloatArrayToPyList((float *) F->data, n_elem);
  }

  return (PConvAutoNone(result));
#endif
}

/*==============================================================================*/
PyObject * ObjectVolumeGetRamp(ObjectVolume * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* TODO: Allow for multi-state maps? */
  PyObject * result = NULL;
  ObjectVolumeState *ovs;

  if(I && (ovs = ObjectVolumeGetActiveState(I))) {
    if(!ovs->isUpdated)
      ObjectVolumeUpdate(I);

    result = PConvFloatArrayToPyList(ovs->Ramp, 5 * ovs->RampSize);
  }
  
  return (PConvAutoNone(result));
#endif
}

/*==============================================================================*/
int ObjectVolumeSetRamp(ObjectVolume * I, float *ramp_list, int list_size)
{
  /* TODO: Allow for multi-state maps? */
  float min_val, max_val;
  ObjectVolumeState *ovs = ObjectVolumeGetActiveState(I);

  ok_assert(1, ovs && ramp_list && list_size > 0);

  FreeP(ovs->Ramp);
  ovs->Ramp = ramp_list;
  ovs->RampSize = list_size / 5;
  ovs->RecolorFlag = true;

  SceneChanged(I->Obj.G);

  return true;
ok_except1:
  PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Errors)
    "ObjectVolumeSetRamp failed" ENDFB(I->Obj.G);
  return false;
}

/*==============================================================================*/
int ObjectVolumeGetIsUpdated(ObjectVolume * I)
{
  /* TODO: Allow for multi-state maps? */
  ObjectVolumeState * ovs = ObjectVolumeGetActiveState(I);

  if (!ovs)
    return -1;

  return ovs->isUpdated;
}

