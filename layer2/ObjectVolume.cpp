
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

#include <algorithm>

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

#define clamp(x,l,h) ((x) < (l) ? (l) : (x) > (h) ? (h) : (x))

static ObjectVolume *ObjectVolumeNew(PyMOLGlobals * G);
static void ObjectVolumeFree(ObjectVolume * I);
static void ObjectVolumeInvalidate(ObjectVolume * I, int rep, int level, int state);
static void ObjectVolumeStateInit(PyMOLGlobals * G, ObjectVolumeState * vs);
static void ObjectVolumeRecomputeExtent(ObjectVolume * I);

static ObjectVolumeState * ObjectVolumeGetActiveState(ObjectVolume * I) {
  int a;
  ok_assert(1, I);
  for(a = 0; a < I->NState; a++)
    if(I->State[a].Active)
      return I->State + a;
ok_except1:
  return NULL;
}

static ObjectMapState * ObjectVolumeStateGetMapState(ObjectVolumeState * vs) {
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
  PyList_SetItem(result, 10, PyInt_FromLong(/* I->CarveFlag */ I->AtomVertex != NULL));
  PyList_SetItem(result, 11, PyFloat_FromDouble(I->CarveBuffer));
  PyList_SetItem(result, 12, I->AtomVertex ?
      PConvFloatVLAToPyList(I->AtomVertex) : PConvAutoNone(NULL));
  PyList_SetItem(result, 13, PyInt_FromLong(0 /* I->VolumeMode */));
  PyList_SetItem(result, 14, PyFloat_FromDouble(0.0 /* I->AltLevel */));
  PyList_SetItem(result, 15, PyInt_FromLong(1 /* I->quiet */));
  if(I->Field) {
    PyList_SetItem(result, 16, IsosurfAsPyList(I->State.G, I->Field));
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
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 10), &I->CarveFlag);
#endif
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 11), &I->CarveBuffer);
      if(ok) {
        tmp = PyList_GetItem(list, 12);
        if(tmp == Py_None)
          I->AtomVertex = NULL;
        else
          ok = PConvPyListToFloatVLA(tmp, &I->AtomVertex);
      }
#if 0
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 13), &I->VolumeMode);
#endif
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

int ObjectVolumeNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectVolume ** result)
{
  int ok = true;
  ObjectVolume *I = NULL;
  (*result) = NULL;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
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
}

PyObject *ObjectVolumeAsPyList(ObjectVolume * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(&I->Obj));
  PyList_SetItem(result, 1, PyInt_FromLong(I->NState));
  PyList_SetItem(result, 2, ObjectVolumeAllStatesAsPyList(I));
  return (PConvAutoNone(result));
}

/*
 * Does actually NOT free the instance, only it's fields.
 */
static void ObjectVolumeStateFree(ObjectVolumeState * vs)
{
  // the instance is only "Active" when it has been initialized. Never free
  // uninitialized instances.
  if(!vs->Active)
    return;
  ObjectStatePurge(&vs->State);
  if(vs->State.G->HaveGUI) {
    vs->State.G->ShaderMgr->freeGPUBuffers(vs->textures, 3);
  }
  if(vs->Field) {
    IsosurfFieldFree(vs->State.G, vs->Field);
    vs->Field = NULL;
  }
  FieldFreeP(vs->carvemask);
  VLAFreeP(vs->AtomVertex);
  if (vs->Ramp)
    FreeP(vs->Ramp);
  vs->Active = false;
}

static void ObjectVolumeFree(ObjectVolume * I)
{
  int a;
  for(a = 0; a < I->NState; a++) {
    ObjectVolumeStateFree(I->State + a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);

  OOFreeP(I);
}

int ObjectVolumeInvalidateMapName(ObjectVolume * I, const char *name, const char * new_name)
{
  int a;
  ObjectVolumeState *vs;
  int result = false;
  for(a = 0; a < I->NState; a++) {
    vs = I->State + a;
    if(vs->Active) {
      if(strcmp(vs->MapName, name) == 0) {
        if (new_name)
          strcpy(vs->MapName, new_name);
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

/*
 * Get the field either from the associated map, or from vs->Field in case
 * this is a reduced or symmetry expanded volume.
 */
static CField * ObjectVolumeStateGetField(ObjectVolumeState * vs) {
  if (!vs)
    return NULL;
  if(vs->Field)
    return vs->Field->data;
  return ObjectVolumeStateGetMapState(vs)->Field->data;
}

CField * ObjectVolumeGetField(ObjectVolume * I) {
  return ObjectVolumeStateGetField(ObjectVolumeGetActiveState(I));
}

/*
 * Get a 4x4 (incl. translation) FracToReal from corner array
 */
static void get44FracToRealFromCorner(const float * corner, float * frac2real)
{
  float tmp[16];
  identity44f(tmp);
  subtract3f(corner +  3, corner, tmp);
  subtract3f(corner +  6, corner, tmp + 4);
  subtract3f(corner + 12, corner, tmp + 8);
  copy3f(corner, tmp + 12);
  transpose44f44f(tmp, frac2real);
}

static void ObjectVolumeUpdate(ObjectVolume * I)
{
  int a;
  ObjectVolumeState *vs;
  ObjectMapState *oms = NULL;
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
      a, vs->RefreshFlag, vs->ResurfaceFlag ENDFD;

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

      // data min/max/mean/stdev
      range = SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_data_range);
      ObjectMapStateGetHistogram(I->Obj.G, oms, 0, range, vs->min_max_mean_stdev, 0.f, 0.f);
    }

    // handle legacy or default color ramp
    if(!vs->Ramp || (vs->RampSize
          && vs->Ramp[0] == 0.f
          && vs->Ramp[5 * (vs->RampSize - 1)] == 359.f)) {

      if(vs->Ramp) {
        // legacy color ramp (0..359)
        range = vs->min_max_mean_stdev[1] - vs->min_max_mean_stdev[0];
        PRINTFB(G, FB_ObjectVolume, FB_Warnings)
          " ObjectVolumeUpdate: detected legacy color ramp\n" ENDFB(G);
        for (i = 0; i < vs->RampSize * 5; i += 5) {
          vs->Ramp[i] = vs->Ramp[i] / 359.f * range + vs->min_max_mean_stdev[0];
        }
      } else {
        // default color ramp (1.0 sigma peak)
        if(!vs->Ramp) {
          float defaultramp[] = {
            vs->min_max_mean_stdev[2] + 0.7f * vs->min_max_mean_stdev[3],
            0.f, 0.f, 1.f, 0.0f,
            vs->min_max_mean_stdev[2] + 1.0f * vs->min_max_mean_stdev[3],
            0.f, 1.f, 1.f, 0.2f,
            vs->min_max_mean_stdev[2] + 1.3f * vs->min_max_mean_stdev[3],
            0.f, 0.f, 1.f, 0.0f
          };
          vs->RecolorFlag = true;
          vs->RampSize = 3;
          vs->Ramp = pymol::malloc<float>(5 * vs->RampSize);
          memcpy(vs->Ramp, defaultramp, 5 * vs->RampSize * sizeof(float));
        }
      }
    }

    if((I->Obj.visRep & cRepVolumeBit) && vs->ResurfaceFlag) {
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
        // get bounds and dimension data from field
        copy3(field->data->dim, vs->dim);
        IsofieldGetCorners(G, field, vs->Corner);

        // transform corners by state matrix
        if(vs->State.Matrix) {
          for(i = 0; i < 8; i++)
            transform44d3f(vs->State.Matrix,
                vs->Corner + 3 * i,
                vs->Corner + 3 * i);
        }
      }

      if(/* CarveFlag */ vs->AtomVertex) {

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
          float vv[3];
          float frac2real[16];

          MapSetupExpress(voxelmap);

          dx = vs->dim[0];
          dy = vs->dim[1];
          dz = vs->dim[2];

          get44FracToRealFromCorner(vs->Corner, frac2real);

          // initialize carve mask
          FieldFreeP(vs->carvemask);
          vs->carvemask = FieldNew(G, (int*) vs->dim, 3, sizeof(GLubyte), cFieldOther);

          // loop over voxels
          for (z = 0; z < dz; z++) {
            for (y = 0; y < dy; y++) {
              for (x = 0; x < dx; x++) {
                float frac[3] = {(x + .5f) / dx, (y + .5f) / dy, (z + .5f) / dz};

                transform44f3f(frac2real, frac, vv);
                flag = avoid_flag;

                // loop over close atoms
                MapLocus(voxelmap, vv, &h, &k, &l);
                for(i = *(MapEStart(voxelmap, h, k, l));
                    i && (j = voxelmap->EList[i]) >= 0; i++) {
                  if(within3f(vs->AtomVertex + 3 * j, vv, carve_buffer)) {
                    flag = !flag;
                    break;
                  }
                }

                // 0xFF (masked) or 0 (not masked), will be 1.0 or 0.0 in shader
                *((GLubyte*)F3p(vs->carvemask, x, y, z)) = flag ? 0x0 : 0xFF;
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

static
int ObjectVolumeAddSlicePoint(float *p0, float *p1, float *zaxis, float d, float *slice, float *t0, float *t1, float *tex_coords, float *origin);
static
void ObjectVolumeDrawSlice(float *points, float *tex_coords, int n_points, float *zaxis);

/*
 * Converting Ramp to `count * 4` sized interpolated RGBA color
 * array. Returns allocated memory.
 * Assigns data minimum and range covered by ramp to `ramp_min` and `ramp_range`.
 */
static float * ObjectVolumeStateGetColors(PyMOLGlobals * G, ObjectVolumeState * ovs,
    int count, float *ramp_min, float *ramp_range) {
  int i, j, k;
  int lowerId, upperId = 0;
  float mixc, mixcincr, r_min, range;
  float stdev = ovs->min_max_mean_stdev[3];
  float * colors;

  ok_assert(1, ovs->Ramp && ovs->RampSize > 1);

  r_min = ovs->Ramp[0];
  range = ovs->Ramp[5 * (ovs->RampSize - 1)] - r_min;

  ok_assert(1, range > R_SMALL4);

  r_min -= stdev * 0.5;
  range += stdev;

  colors = pymol::calloc<float>(4 * count);
  ok_assert(1, colors);

  for (i = 0; i < ovs->RampSize; i++) {
    lowerId = upperId;
    upperId = (int) (count * (ovs->Ramp[i * 5] - r_min) / range);

    if(i == 0)
      continue;

    mixcincr = 1.f / (upperId - lowerId);

    for (j = lowerId, mixc = 1.f; j < upperId; j++, mixc -= mixcincr){
      if(j < 0 || j >= count)
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
 * Adjust alpha values in the given RGBA array (in place) by:
 *   alpha_new = 1 - exp(-alpha * factor)
 */
static void ColorsAdjustAlpha(float * colors, int count, float factor) {
  int j;
  for (j = 3; j < count * 4; j += 4) {
    colors[j] = 1. - expf(-colors[j] * factor);
  }
}

/*
 * Render bounding box
 * TODO: duplicate in other reps?
 */
static void ExtentRender(float * corner) {
#ifndef PURE_OPENGL_ES_2
  int i, ci[] = {
     0,  3,  3,  9,  9,  6,  6,  0,
    12, 15, 15, 21, 21, 18, 18, 12,
     0, 12,  3, 15,  9, 21,  6, 18
  };
  glBegin(GL_LINES);
  for(i = 0; i < 8 * 3; i++)
    glVertex3fv(corner + ci[i]);
  glEnd();
#endif
}

static size_t createColorTexture(PyMOLGlobals * G, const float *colors, const int count)
{
  size_t texname = 0;
#ifndef PURE_OPENGL_ES_2
  auto tex = G->ShaderMgr->newGPUBuffer<textureBuffer_t>(
    tex::format::RGBA,
    tex::data_type::FLOAT,
    tex::filter::LINEAR,
    tex::filter::LINEAR,
    tex::wrap::CLAMP
    );
  tex->texture_data_1D(count, colors);
  texname = tex->get_hash_id();
#endif
  return texname;
}

static size_t createPreintegrationTexture(PyMOLGlobals * G, const float *Table, const int count)
{
  float factor, tmp1[4];
  Vector4f *sat = pymol::malloc<Vector4f>(count + 1);
  int i, sb, sf, lookupindex = 0;
  GLfloat *lookupImg = pymol::malloc<GLfloat>(count * count * 4);

  memset(sat[0], 0, sizeof(sat[0]));

  // summed area table
  for (i = 0; i < count; i++) {
    tmp1[3] = Table[i * 4 + 3];
    scale3f(Table + i * 4, tmp1[3], tmp1);
    add4f(tmp1, sat[i], sat[i + 1]);
}

  // make quadratic lookup table
  for (sb = 0; sb < count; sb++) {
    for (sf = 0; sf < count; sf++) {
      GLfloat col[4];
      int smin, smax;

      if (sb < sf) { smin=sb; smax=sf; } else { smin=sf; smax=sb; }

      if (sat[smax + 1][3] != sat[smin][3]) {
        factor = 1.f / (sat[smax + 1][3] - sat[smin][3]);

        for (i = 0; i < 3; i++)
          col[i] = (sat[smax + 1][i] - sat[smin][i]) * factor;

        col[3] = 1. / (factor * (smax + 1 - smin));

      } else {
        for (i = 0; i < 4; i++)
          col[i] = 0.f;
      }
      for (i = 0; i < 4; i++)
        lookupImg[lookupindex++] = clamp(col[i], 0., 1.);
    }
  }

  // upload texture
  auto tex = G->ShaderMgr->newGPUBuffer<textureBuffer_t>(
    tex::format::RGBA,
    tex::data_type::FLOAT,
    tex::filter::NEAREST,
    tex::filter::NEAREST,
    tex::wrap::CLAMP_TO_EDGE,
    tex::wrap::CLAMP_TO_EDGE
    );
  tex->texture_data_2D(count, count, lookupImg);

  mfree(sat);
  mfree(lookupImg);

  return tex->get_hash_id();
}

static void ObjectVolumeRender(ObjectVolume * I, RenderInfo * info)
{
#ifndef PURE_OPENGL_ES_2
  PyMOLGlobals *G = I->Obj.G;
  int state = info->state;
  int pass = info->pass;
  int a = 0;
  ObjectVolumeState *vs = NULL;
  float volume_layers =  SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_layers);
  short volume_mode = SettingGetGlobal_i(G, cSetting_volume_mode);
  short ortho = SettingGetGlobal_i(G, cSetting_ortho);
  /* make this a setting? */
  GLint alpha_func;
  GLfloat alpha_ref;
  float tex_corner[24];
  float *corner;
  const float *ttt;
  float zaxis[3];
  float points[36], tex_coords[36];
  int n_points;
  float d, sliceRange, sliceDelta;
  float origin[3];
  CShaderPrg *shaderPrg;
  bool volume_t = 0;

  if(info->pick || pass != -1)
    return;

  if(!G->HaveGUI || !G->ValidContext)
    return;

  /* bail if no shaders */
  if (G && !(G->ShaderMgr->ShadersPresent()))
      return;

  if (info->pass < 0){
    volume_t = SettingGetGlobal_i(G, cSetting_transparency_mode) == 3;
  }

  // ViewElem/TTT Matrix
  ObjectPrepareContext(&I->Obj, info);

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
    if((I->Obj.visRep & cRepExtentBit)) {
      if(!info->line_lighting)
        glDisable(GL_LIGHTING);
      ObjectUseColor(&I->Obj);
      ExtentRender(corner);
    }

    // upload color ramp texture
    if (vs->RecolorFlag) {
      const int volume_nColors = 512;

      float * colors = ObjectVolumeStateGetColors(G, vs, volume_nColors, &vs->ramp_min, &vs->ramp_range);
      if(!colors)
        continue;

      // volume_layers default is 256, adjust alpha to maintain integrated
      // opacity with different layer numbers
      ColorsAdjustAlpha(colors, volume_nColors, 256. / volume_layers);

      if (vs->textures[1]) {
        G->ShaderMgr->freeGPUBuffer(vs->textures[1]);
      }

      vs->textures[1] =
#ifdef _PYMOL_IP_EXTRAS
#endif
        createColorTexture(G, colors, volume_nColors);

      tex::env(tex::env_name::ENV_MODE, tex::env_param::REPLACE);

      mfree(colors);
      vs->RecolorFlag = false;
    }

    // upload map data texture
    if (!vs->textures[0] || vs->RefreshFlag) {
      tex::data_type volume_bit_depth;
      CField * field = ObjectVolumeStateGetField(vs);

      if(!field) {
        PRINTFB(G, FB_ObjectVolume, FB_Errors)
          " ObjectVolumeRender-Error: Could not get field data.\n" ENDFB(G);
        return;
      }

      int volume_bit_val = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_volume_bit_depth);
      volume_bit_depth = (volume_bit_val < 17) ? tex::data_type::HALF_FLOAT : tex::data_type::FLOAT;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#if 0
      glTexImage3D = getTexImage3D(); 
      if (! glTexImage3D) {
        PRINTFB(G, FB_ObjectVolume, FB_Errors)
          " ObjectVolumeRender-Error: Could not bind the glActiveTexture or glTexImage3D function.\n"
          ENDFB(G);
        return;
      } 
#endif
/* END PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */

      if (vs->textures[0]) {
        G->ShaderMgr->freeGPUBuffer(vs->textures[0]);
        vs->textures[0] = 0;
      }
      if (vs->textures[2]) {
        G->ShaderMgr->freeGPUBuffer(vs->textures[2]);
        vs->textures[2] = 0;
      }

      auto tex3dGenBind = [](PyMOLGlobals * G, tex::data_type dtype) -> size_t {
        auto texture = G->ShaderMgr->newGPUBuffer<textureBuffer_t>(
          tex::format::R,
          dtype,
          tex::filter::LINEAR,
          tex::filter::LINEAR,
          tex::wrap::CLAMP,
          tex::wrap::CLAMP,
          tex::wrap::CLAMP
          );
        tex::env(tex::env_name::ENV_MODE, tex::env_param::REPLACE);
        return texture->get_hash_id();
      };
      // Create a 3D texture
      vs->textures[0] = tex3dGenBind(G, volume_bit_depth);
      auto t0 = G->ShaderMgr->getGPUBuffer<textureBuffer_t>(vs->textures[0]);
      t0->texture_data_3D(field->dim[2], field->dim[1], field->dim[0], field->data);

      // Create 3D carve mask texture
      if(vs->carvemask) {
        vs->textures[2] = tex3dGenBind(G, tex::data_type::UBYTE);
        auto t2 = G->ShaderMgr->getGPUBuffer<textureBuffer_t>(vs->textures[2]);
        t2->texture_data_3D(
          vs->carvemask->dim[2],
          vs->carvemask->dim[1],
          vs->carvemask->dim[0],
          vs->carvemask->data
          );

        // not needed anymore, data now in texture memory
        FieldFreeP(vs->carvemask);
      }

      vs->RefreshFlag = false;
    }

    // render volume
    if((I->Obj.visRep & cRepVolumeBit)) {
      int i, j;

      glDisable(GL_LIGHTING);

      // half grid cell inset in texture corners (texture coordinate units)
      for(j = 0; j < 3; j++) {
        float offset = 0.5 / vs->dim[2 - j];
        int bit = 1 << (2 - j);
        for(i = 0; i < 8; i++)
          tex_corner[i * 3 + j] = (i & bit) ? 1.0 - offset : offset;
      }

      // for z-axis
      SceneGetViewNormal(G, zaxis);

      for(j = 0; j < 3; j++) {
        // map center
        origin[j] = corner[j] + 0.5 * (corner[21 + j] - corner[j]);
      }

      // TTT (movie object motions)
      if(ObjectGetTTT(&I->Obj, &ttt, -1))
        MatrixTransformC44fAs33f3f(ttt, zaxis, zaxis);

      // determine number of slices based on max extent
      // and slice option
      sliceRange = 0.5*sqrt(2.0) * 
        std::max(std::max(fabs(corner[21]-corner[0]), fabs(corner[22]-corner[1])), 
            fabs(corner[23]-corner[2]));
      sliceDelta = (sliceRange / volume_layers);

      // load shader
      shaderPrg = G->ShaderMgr->GetShaderPrg(volume_t ? "volume_t" : "volume");
      if (!shaderPrg)
	return;
      shaderPrg->Enable();
      shaderPrg->Set_Stereo_And_AnaglyphMode();
      shaderPrg->Set1i("volumeTex", 0);
      shaderPrg->Set1i("colorTex1D", 1);
      shaderPrg->Set1i("colorTex2D", 1);
      shaderPrg->Set1i("carvemask", 5);
      shaderPrg->Set1i("carvemaskFlag", vs->textures[2] != 0);
      shaderPrg->Set1f("volumeScale", 1.0 / vs->ramp_range);
      shaderPrg->Set1f("volumeBias", (-vs->ramp_min) / vs->ramp_range);

      // for pre-integrated rendering
#ifdef _PYMOL_IP_EXTRAS
#endif

      // background and fog stuff
      {
        shaderPrg->SetBgUniforms();
      }

      // bind color ramp and map data textures
      glActiveTexture(GL_TEXTURE1);
      G->ShaderMgr->bindGPUBuffer(vs->textures[1]);

      glActiveTexture(GL_TEXTURE0);
      G->ShaderMgr->bindGPUBuffer(vs->textures[0]);

      if (vs->textures[2]) {
        glActiveTexture(GL_TEXTURE5);
        G->ShaderMgr->bindGPUBuffer(vs->textures[2]);
      }

      // alpha: everything passes
      // Not sure if we really need to restore this
      glGetIntegerv(GL_ALPHA_TEST_FUNC, &alpha_func);
      glGetFloatv(GL_ALPHA_TEST_REF, &alpha_ref);
      glAlphaFunc(GL_ALWAYS, 0.0);

      // don't write to the depth buffer
      GLboolean depth_writemask;
      glGetBooleanv(GL_DEPTH_WRITEMASK, &depth_writemask);
      if (depth_writemask)
        glDepthMask(GL_FALSE);

      // Cheap hack, should be replaced with non-immediate calls
      glFlush();
      glFinish();

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

      shaderPrg->Disable();

      // restore
      if (depth_writemask)
        glDepthMask(GL_TRUE);
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
      if (a < 0.0f) a += 2.0f * PI;
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
    glBegin(GL_POLYGON);
    for (i=0; i<n_points; i++) {
      glTexCoord3fv(&tex_coords[3 * vertices[i]]);
      glVertex3fv(&points[3 * vertices[i]]);
    }
    glEnd();
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
  vs->CarveBuffer = 0.0;
  vs->AtomVertex = NULL;
  vs->caption[0] = 0;
  zero3i(vs->dim);
  vs->carvemask = NULL;
  vs->textures[0] = 0; // 3D volume (map)
  vs->textures[1] = 0; // 1D/2D color table
  vs->textures[2] = 0; // 3D carvemask
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
                                  float level, int box_mode,
                                  float carve, float *vert_vla,
                                  float alt_level, int quiet)
{
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

      if(sym && box_mode) {
        int eff_range[6];

        IsosurfGetRange(G, oms->Field, oms->Symmetry->Crystal, min_ext, max_ext, eff_range, false);

        {
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

  vs->CarveBuffer = carve;
  vs->AtomVertex = vert_vla;

  I->Obj.ExtentFlag = false;

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
    const float *ttt;
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
PyObject * ObjectVolumeGetRamp(ObjectVolume * I)
{
  /* TODO: Allow for multi-state maps? */
  PyObject * result = NULL;
  ObjectVolumeState *ovs;

  if(I && (ovs = ObjectVolumeGetActiveState(I))) {
    if(!ovs->isUpdated)
      ObjectVolumeUpdate(I);

    result = PConvFloatArrayToPyList(ovs->Ramp, 5 * ovs->RampSize);
  }
  
  return (PConvAutoNone(result));
}

/*==============================================================================*/
int ObjectVolumeSetRamp(ObjectVolume * I, float *ramp_list, int list_size)
{
  /* TODO: Allow for multi-state maps? */
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
