
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"OOMac.h"
#include"ObjectSlice.h"
#include"Base.h"
#include"Matrix.h"
#include"MemoryDebug.h"
#include"Map.h"
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
#include"CGO.h"
#include"ShaderMgr.h"

#define START_STRIP -1
#define STOP_STRIP -2

static void ObjectSliceRecomputeExtent(ObjectSlice * I);

static PyObject *ObjectSliceStateAsPyList(ObjectSliceState * I)
{

  PyObject *result = NULL;

  result = PyList_New(10);

  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  PyList_SetItem(result, 1, PyString_FromString(I->MapName));
  PyList_SetItem(result, 2, PyInt_FromLong(I->MapState));
  PyList_SetItem(result, 3, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 4, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 5, PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->origin, 3));
  PyList_SetItem(result, 7, PConvFloatArrayToPyList(I->system, 9));
  PyList_SetItem(result, 8, PyFloat_FromDouble(I->MapMean));
  PyList_SetItem(result, 9, PyFloat_FromDouble(I->MapStdev));

  return (PConvAutoNone(result));
}

static PyObject *ObjectSliceAllStatesAsPyList(ObjectSlice * I)
{

  PyObject *result = NULL;
  int a;
  result = PyList_New(I->State.size());
  for(a = 0; a < I->State.size(); a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result, a, ObjectSliceStateAsPyList(I->State.data() + a));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(NULL));
    }
  }
  return (PConvAutoNone(result));

}

static int ObjectSliceStateFromPyList(PyMOLGlobals * G, ObjectSliceState * I,
                                      PyObject * list)
{
  int ok = true;

  if(ok)
    ok = (list != NULL);
  if(ok) {
    if(!PyList_Check(list))
      I->Active = false;
    else {
      if(ok)
        ok = (list != NULL);
      if(ok)
        ok = PyList_Check(list);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */

      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->Active);
      if(ok)
        ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->MapName, WordLength);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->MapState);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 3), I->ExtentMin, 3);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 4), I->ExtentMax, 3);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 5), &I->ExtentFlag);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6), I->origin, 3);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 7), I->system, 9);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 8), &I->MapMean);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 9), &I->MapStdev);

      I->RefreshFlag = true;
    }
  }

  return (ok);
}

static int ObjectSliceAllStatesFromPyList(ObjectSlice * I, PyObject * list, int size)
{
  int ok = true;
  int a;
  if(ok)
    ok = PyList_Check(list);
  assert(size == PyList_Size(list));
  if(ok) {
    for(a = 0; a < size; a++) {
      CPythonVal *val = CPythonVal_PyList_GetItem(I->G, list, a);
      I->State.emplace_back(I->G);
      ok = ObjectSliceStateFromPyList(I->G, &I->State[a], val);
      CPythonVal_Free(val);
      if(!ok)
        break;
    }
  }
  return (ok);
}

int ObjectSliceNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectSlice ** result)
{
  int ok = true;
  ObjectSlice *I = NULL;
  (*result) = NULL;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  I = new ObjectSlice(G);
  if(ok)
    ok = (I != NULL);

  if(ok){
    CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 0);    
    ok = ObjectFromPyList(G, val, I);
    CPythonVal_Free(val);
  }
  int size;
  if(ok)
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 1, &size);
  if(ok){
    CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 2);
    ok = ObjectSliceAllStatesFromPyList(I, val, size);
    CPythonVal_Free(val);
  }
  if(ok) {
    (*result) = I;
    ObjectSliceRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
}

PyObject *ObjectSliceAsPyList(ObjectSlice * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->State.size()));
  PyList_SetItem(result, 2, ObjectSliceAllStatesAsPyList(I));

  return (PConvAutoNone(result));
}

void ObjectSlice::invalidate(cRep_t rep, cRepInv_t level, int state)
{
  int a;
  int once_flag = true;
  for(a = 0; a < State.size(); a++) {
    if(state < 0)
      once_flag = false;
    if(!once_flag)
      state = a;
    State[state].RefreshFlag = true;
    SceneChanged(G);
    if(once_flag)
      break;
  }
}

static void ObjectSliceStateAssignColors(ObjectSliceState * oss, ObjectGadgetRamp * ogr)
{
  /* compute the colors */
  if(oss && oss->values && oss->colors) {
    int x, y;
    int *min = oss->min;
    int *max = oss->max;
    float *value = oss->values.data();
    int *flag = oss->flags.data();
    float *color = oss->colors.data();
    for(y = min[1]; y <= max[1]; y++) {
      for(x = min[0]; x <= max[0]; x++) {
        if(*flag) {
          ObjectGadgetRampInterpolate(ogr, *value, color);
          ColorLookupColor(oss->G, color);
        }
        color += 3;
        value++;
        flag++;
      }
    }
  }
}

static void ObjectSliceStateUpdate(ObjectSlice * I, ObjectSliceState * oss,
                                   ObjectMapState * oms)
{
  int ok = true;
  int min[2] = { 0, 0 }, max[2] = {
  0, 0};                        /* limits of the rectangle */
  int need_normals = false;
  int track_camera = SettingGet_b(I->G, NULL, I->Setting.get(), cSetting_slice_track_camera);
  float grid = SettingGet_f(I->G, NULL, I->Setting.get(), cSetting_slice_grid);
  int min_expand = 1;

  if(SettingGet_b(I->G, NULL, I->Setting.get(), cSetting_slice_dynamic_grid)) {
    float resol = SettingGet_f(I->G, NULL, I->Setting.get(),
                               cSetting_slice_dynamic_grid_resolution);
    float scale = SceneGetScreenVertexScale(I->G, oss->origin);
    oss->last_scale = scale;
    grid = resol * scale;
  }
  oss->shaderCGO.reset();
  if (track_camera){
    oss->outline_n_points = 0;
  }
  if(grid < 0.01F)
    grid = 0.01F;

  /* for the given map, compute a new set of interpolated points with accompanying levels */

  /* first, find the limits of the enclosing rectangle, starting at the slice origin, 
     via a simple brute-force approach... */

  if(oss->ExtentFlag) {         /* how far out do we need to go to be sure we intersect the map? */
    min_expand = (int) (diff3f(oss->ExtentMax, oss->ExtentMin) / grid);
  }
  if(ok) {
    int size = 1, minus_size;
    int a;
    int keep_going = true;
    int n_cycles = 0;
    float point[3];

    while(keep_going || (n_cycles < min_expand)) {
      keep_going = false;
      minus_size = -size;
      n_cycles++;

      for(a = -size; a <= size; a++) {

        if((max[1] != size) || (min[0] > a) || (max[0] < a)) {
          point[0] = grid * a;
          point[1] = grid * size;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms, point)) {
            keep_going = true;
            if(max[1] < size)
              max[1] = size;
            if(min[0] > a)
              min[0] = a;
            if(max[0] < a)
              max[0] = a;
          }

        } else
          keep_going = true;

        if((min[1] != minus_size) || (min[0] > a) || (max[0] < a)) {
          point[0] = grid * a;
          point[1] = grid * minus_size;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms, point)) {
            keep_going = true;
            if(min[1] > minus_size)
              min[1] = minus_size;
            if(min[0] > a)
              min[0] = a;
            if(max[0] < a)
              max[0] = a;

          }
        } else
          keep_going = true;

        if((max[0] != size) || (min[1] > a) || (max[1] < a)) {
          point[0] = grid * size;
          point[1] = grid * a;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms, point)) {
            keep_going = true;
            if(max[0] < size)
              max[0] = size;
            if(min[1] > a)
              min[1] = a;
            if(max[1] < a)
              max[1] = a;
          }
        } else
          keep_going = true;

        if((min[0] != minus_size) || (min[1] > a) || (max[1] < a)) {
          point[0] = grid * minus_size;
          point[1] = grid * a;
          point[2] = 0.0F;
          transform33f3f(oss->system, point, point);
          add3f(oss->origin, point, point);
          if(ObjectMapStateContainsPoint(oms, point)) {
            keep_going = true;
            if(min[0] > minus_size)
              min[0] = minus_size;
            if(min[1] > a)
              min[1] = a;
            if(max[1] < a)
              max[1] = a;
          }
        } else
          keep_going = true;
      }
      if(keep_going)
        min_expand = 0;         /* if we've hit, then don't keep searching blindly */

      size++;
    }
    oss->max[0] = max[0];
    oss->max[1] = max[1];
    oss->min[0] = min[0];
    oss->min[1] = min[1];
  }
  /* now confirm that storage is available */
  if(ok) {
    int n_alloc = (1 + oss->max[0] - oss->min[0]) * (1 + oss->max[1] - oss->min[1]);

    oss->points.reserve(n_alloc * 3);
    oss->values.reserve(n_alloc);
    oss->colors.reserve(n_alloc * 3);
    oss->flags.reserve(n_alloc);

    if(!(oss->points && oss->values && oss->flags)) {
      ok = false;
      PRINTFB(I->G, FB_ObjectSlice, FB_Errors)
        "ObjectSlice-Error: allocation failed\n" ENDFB(I->G);
    }

    if(!oss->strips)            /* this is range-checked during use */
      oss->strips = pymol::vla<int>(n_alloc);

    oss->n_points = n_alloc;
  }

  /* generate the coordinates */

  if(ok) {
    int x, y;
    float *point = oss->points.data();
    for(y = min[1]; y <= max[1]; y++) {
      for(x = min[0]; x <= max[0]; x++) {
        point[0] = grid * x;
        point[1] = grid * y;
        point[2] = 0.0F;
        transform33f3f(oss->system, point, point);
        add3f(oss->origin, point, point);
        point += 3;
      }
    }
  }

  /* interpolate and flag the points inside the map */

  if(ok) {
    ObjectMapStateInterpolate(oms, oss->points.data(), oss->values.data(), oss->flags.data(), oss->n_points);
  }

  /* apply the height scale (if nonzero) */

  if(ok) {

    if(SettingGet_b(I->G, NULL, I->Setting.get(), cSetting_slice_height_map)) {
      float height_scale =
        SettingGet_f(I->G, NULL, I->Setting.get(), cSetting_slice_height_scale);
      float *value = oss->values.data();
      float up[3], scaled[3], factor;
      int x, y;
      float *point = oss->points.data();

      need_normals = true;
      up[0] = oss->system[2];
      up[1] = oss->system[5];
      up[2] = oss->system[8];

      for(y = min[1]; y <= max[1]; y++) {
        for(x = min[0]; x <= max[0]; x++) {
          factor = ((*value - oss->MapMean) / oss->MapStdev) * height_scale;
          scale3f(up, factor, scaled);
          add3f(scaled, point, point);
          point += 3;
          value++;
        }
      }
      /* TODO: For all edge points, move them onto the closest outline line, generated by GenerateOutlineOfSlice */
    }
  }

/* now generate efficient triangle strips based on the points that are present in the map */

  if(ok) {
    int x, y;
    int cols = 1 + max[0] - min[0];
    int flag00, flag01, flag10, flag11;
    int offset = 0, offset00, offset01, offset10, offset11;
    int strip_active = false;
    int n = 0;
    for(y = min[1]; y < max[1]; y++) {
      offset00 = offset;
      for(x = min[0]; x < max[0]; x++) {
        offset01 = offset00 + 1;
        offset10 = offset00 + cols;
        offset11 = offset10 + 1;

        flag00 = oss->flags[offset00];
        flag01 = oss->flags[offset01];
        flag10 = oss->flags[offset10];
        flag11 = oss->flags[offset11];

        /* first triangle - forward handedness: 10 00 11 */
        if(strip_active) {
          if(flag10 && flag00 && flag11) {
            /* continue current strip */

            VLACheck(oss->strips, int, n);
            oss->strips[n] = offset10;
            n++;
          } else {
            /* terminate current strip */

            VLACheck(oss->strips, int, n);
            oss->strips[n] = STOP_STRIP;
            strip_active = false;
            n++;
          }
        } else if(flag10 & flag00 && flag11) {
          /* start a new strip with correct parity */

          VLACheck(oss->strips, int, n + 3);
          oss->strips[n] = START_STRIP;
          oss->strips[n + 1] = offset10;
          oss->strips[n + 2] = offset00;
          oss->strips[n + 3] = offset11;
          n += 4;
          strip_active = true;
        }

        /* second triangle -- reverse handedness: 00 11 01 */

        if(strip_active) {
          if(flag00 && flag11 && flag01) {
            /* continue current strip */
            VLACheck(oss->strips, int, n);
            oss->strips[n] = offset01;
            n++;
          } else {
            /* terminate current strip */
            VLACheck(oss->strips, int, n);
            oss->strips[n] = STOP_STRIP;
            strip_active = false;
            n++;
          }
        } else if(flag00 & flag11 && flag01) {
          /* singleton triangle -- improper order for strip */

          VLACheck(oss->strips, int, n + 5);
          oss->strips[n + 0] = START_STRIP;
          oss->strips[n + 1] = offset11;
          oss->strips[n + 2] = offset00;
          oss->strips[n + 3] = offset01;
          oss->strips[n + 4] = STOP_STRIP;
          n += 5;
        }
        offset00++;
      }
      if(strip_active) {
        /* terminate current strip */
        VLACheck(oss->strips, int, n);
        oss->strips[n] = STOP_STRIP;
        strip_active = false;
        n++;
      }
      offset += cols;
    }
    VLACheck(oss->strips, int, n);
    n++;
    oss->n_strips = n;

  }

  /* compute triangle normals if we need them */

  if(!need_normals) {
    VLAFreeP(oss->normals);
  } else {
    int *cnt = NULL;

    oss->normals.reserve(oss->n_points * 3);
    cnt = pymol::calloc<int>(oss->n_points);

    if(cnt && oss->normals) {
      int *strip = oss->strips.data();
      float *point = oss->points.data();
      float *normal = oss->normals.data();
      int n = oss->n_strips;
      int a;
      int offset0 = 0, offset1 = 0, offset2, offset;
      int strip_active = false;
      int tri_count = 0;

      float d1[3], d2[3], cp[3];
      UtilZeroMem(oss->normals.data(), sizeof(float) * 3 * oss->n_points);

      for(a = 0; a < n; a++) {
        offset = *(strip++);
        switch (offset) {
        case START_STRIP:
          strip_active = true;
          tri_count = 0;
          break;
        case STOP_STRIP:
          strip_active = false;
          break;
        default:
          if(strip_active) {
            tri_count++;
            offset2 = offset1;
            offset1 = offset0;
            offset0 = offset;

            if(tri_count >= 3) {

              if(tri_count & 0x1) {     /* get the handedness right ... */
                subtract3f(point + 3 * offset1, point + 3 * offset0, d1);
                subtract3f(point + 3 * offset2, point + 3 * offset1, d2);
              } else {
                subtract3f(point + 3 * offset0, point + 3 * offset1, d1);
                subtract3f(point + 3 * offset2, point + 3 * offset0, d2);
              }
              cross_product3f(d2, d1, cp);
              normalize3f(cp);
              add3f(cp, normal + 3 * offset0, normal + 3 * offset0);
              add3f(cp, normal + 3 * offset1, normal + 3 * offset1);
              add3f(cp, normal + 3 * offset2, normal + 3 * offset2);
              cnt[offset0]++;
              cnt[offset1]++;
              cnt[offset2]++;
            }
          }
        }
      }

      {                         /* now normalize the average normals for active vertices */
        int x, y;
        float *normal = oss->normals.data();
        int *c = cnt;
        for(y = min[1]; y <= max[1]; y++) {
          for(x = min[0]; x <= max[0]; x++) {
            if(*c)
              normalize3f(normal);
            point += 3;
            c++;
          }
        }
      }
    }
    FreeP(cnt);
  }
}

void ObjectSlice::update()
{

  auto I = this;
  ObjectSliceState *oss;
  ObjectMapState *oms = NULL;
  ObjectMap *map = NULL;
  ObjectGadgetRamp *ogr = NULL;

  int a;
  for(a = 0; a < I->State.size(); a++) {
    oss = &I->State[a];
    if(oss && oss->Active) {
      map = ExecutiveFindObjectMapByName(I->G, oss->MapName);
      if(!map) {
        PRINTFB(I->G, FB_ObjectSlice, FB_Errors)
          "ObjectSliceUpdate-Error: map '%s' has been deleted.\n", oss->MapName
          ENDFB(I->G);
      }
      if(map) {
        oms = ObjectMapGetState(map, oss->MapState);
      }
      if(oms) {

        if(oss->RefreshFlag) {
          oss->RefreshFlag = false;
          PRINTFB(I->G, FB_ObjectSlice, FB_Blather)
            " ObjectSlice: updating \"%s\".\n", I->Name ENDFB(I->G);
          if(oms->Field) {
            ObjectSliceStateUpdate(I, oss, oms);
            ogr = ColorGetRamp(I->G, I->Color);
            if(ogr)
              ObjectSliceStateAssignColors(oss, ogr);
            else {              /* solid color */
              const float *solid = ColorGet(I->G, I->Color);
              float *color = oss->colors.data();
              for(a = 0; a < oss->n_points; a++) {
                *(color++) = solid[0];
                *(color++) = solid[1];
                *(color++) = solid[2];
              }
            }
          }
        }
      }
      SceneInvalidate(I->G);
    }
  }
}

void ObjectSliceDrag(ObjectSlice * I, int state, int mode, float *pt, float *mov,
                     float *z_dir)
{
  ObjectSliceState *oss = NULL;

  if(state >= 0)
    if(state < I->State.size())
      if(I->State[state].Active)
        oss = &I->State[state];

  if(oss) {
    switch (mode) {
    case cButModeRotFrag:      /* rotated about origin */
    case cButModeRotObj:
      {
        float v3[3];
        float n0[3];
        float n1[3];
        float n2[3];
        float cp[3];
        float mat[9];
        float theta;

        copy3f(oss->origin, v3);

        subtract3f(pt, v3, n0);
        add3f(pt, mov, n1);
        subtract3f(n1, v3, n1);

        normalize3f(n0);
        normalize3f(n1);
        cross_product3f(n0, n1, cp);

        theta = (float) asin(length3f(cp));

        normalize23f(cp, n2);

        rotation_matrix3f(theta, n2[0], n2[1], n2[2], mat);

        multiply33f33f(mat, oss->system, oss->system);

        I->invalidate(cRepSlice, cRepInvAll, state);
        SceneInvalidate(I->G);

      }
      break;
    case cButModeMovFrag:      /* move along "up" direction */
    case cButModeMovFragZ:
    case cButModeMovObj:
    case cButModeMovObjZ:
      {
        float up[3], v1[3];
        up[0] = oss->system[2];
        up[1] = oss->system[5];
        up[2] = oss->system[8];

        project3f(mov, up, v1);
        add3f(v1, oss->origin, oss->origin);
        I->invalidate(cRepSlice, cRepInvAll, state);
        SceneInvalidate(I->G);
      }
      break;
    case cButModeTorFrag:
      break;
    }
  }
}

int ObjectSliceGetVertex(ObjectSlice * I, int index, int base, float *v)
{
  int state = index - 1;
  int offset = base - 1;
  int result = false;

  ObjectSliceState *oss = NULL;

  if(state >= 0)
    if(state < I->State.size())
      if(I->State[state].Active)
        oss = &I->State[state];

  if(oss) {
    if((offset >= 0) && (offset < oss->n_points))
      if(oss->flags[offset]) {
        copy3f(oss->points + 3 * offset, v);
        result = true;
      }
  }
  return (result);
}

#if 0
static int ObjectSliceAddSlicePoint(float *pt0, float *pt1, float *zaxis, float d,
			     float *coords, float *origin)
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
    return 3;
  }
  return 0;
}

#ifndef PURE_OPENGL_ES_2
static
void ObjectSliceDrawSlice(CGO *cgo, float *points, int n_points, float *zaxis)
{
  float center[3], v[3], w[3], q[3];
  float angles[12];
  float a, c, s;
  int vertices[12];
  int i, j;

  if (!n_points) return;
  
  // Calculate the polygon center
  center[0] = center[1] = center[2] = 0.0;
  
  for (i=0; i<3*n_points; i+=3) {
    center[0] += points[i];
    center[1] += points[i+1];
    center[2] += points[i+2];
  }    
  
  center[0] /= (float)n_points;
  center[1] /= (float)n_points;
  center[2] /= (float)n_points;
  
  v[0] = points[0]-center[0];
  v[1] = points[1]-center[1];
  v[2] = points[2]-center[2];
  
  normalize3f(v);
  
  // Sort vertices by rotation angle around the central axis
  for (i=0; i<n_points; i++) {
    w[0] = points[3*i]-center[0];
    w[1] = points[3*i+1]-center[1];
    w[2] = points[3*i+2]-center[2];
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
  if (cgo){
    CGOBegin(cgo, GL_LINE_LOOP);
    for (i=0; i<n_points; i++) {
      CGOVertexv(cgo, &points[3*vertices[(i)%n_points]]);
    }
    CGOEnd(cgo);
  } else {
    glBegin(GL_LINE_LOOP);
    for (i=0; i<n_points; i++) {
      glVertex3fv(&points[3*vertices[(i)%n_points]]);
    }
    glEnd();
  }
}

static
void GenerateOutlineOfSlice(PyMOLGlobals *G, ObjectSliceState *oss, CGO *cgo){
  int n_points = oss->outline_n_points;
  float *points = oss->outline_points;
  float *m = SceneGetMatrix(G);
  float *zaxis = oss->outline_zaxis, *origin;//, origin[3];
  float d = 0.f; // not sure what this should be
  if (!oss->outline_n_points){
    zaxis[0] = m[2];
    zaxis[1] = m[6];
    zaxis[2] = m[10];
    /*
    origin[0] = oss->Corner[0] + 0.5*(oss->Corner[21]-oss->Corner[0]);
    origin[1] = oss->Corner[1] + 0.5*(oss->Corner[22]-oss->Corner[1]);
    origin[2] = oss->Corner[2] + 0.5*(oss->Corner[23]-oss->Corner[2]);
    */
    origin = oss->origin;
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[0],&oss->Corner[3],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[3],&oss->Corner[9],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[9],&oss->Corner[6],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[6],&oss->Corner[0],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[12],&oss->Corner[15],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[15],&oss->Corner[21],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[21],&oss->Corner[18],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[18],&oss->Corner[12],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[0],&oss->Corner[12],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[3],&oss->Corner[15],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[9],&oss->Corner[21],zaxis,d,&points[n_points], origin);
    n_points += ObjectSliceAddSlicePoint(&oss->Corner[6],&oss->Corner[18],zaxis,d,&points[n_points], origin);
    oss->outline_n_points = n_points;
  }
  if (cgo){
    CGOColor(cgo, 1.f, 0.f, 0.f );
    CGOSphere(cgo, origin, 1.f);
    CGOColor(cgo, 1.f, 1.f, 1.f );
  } else {
    glColor3f(1.f,1.f,1.f);
  }
  ObjectSliceDrawSlice(cgo, points, n_points/3, zaxis);
}
#endif
#endif


void ObjectSlice::render(RenderInfo * info)
{

  auto I = this;
  int state = info->state;
  CRay *ray = info->ray;
  auto pick = info->pick;
  const RenderPass pass = info->pass;
  int cur_state = 0;
  float alpha;
  int track_camera = SettingGet_b(G, NULL, I->Setting.get(), cSetting_slice_track_camera);
  int dynamic_grid = SettingGet_b(G, NULL, I->Setting.get(), cSetting_slice_dynamic_grid);
  ObjectSliceState *oss = NULL;
  int use_shaders = !track_camera && SettingGet_b(G, NULL, I->Setting.get(), cSetting_use_shaders);
  // just in case, since slice uses immediate mode, but this should never happen
  G->ShaderMgr->Disable_Current_Shader();

  if(track_camera || dynamic_grid) {
    int update_flag = false;

    if(state >= 0)
      if(state < I->State.size())
        if(I->State[state].Active)
          oss = &I->State[state];

    while(1) {
      if(state < 0) {           /* all_states */
        oss = &I->State[cur_state];
      } else {
        if(oss) {

          SceneViewType view;
          float pos[3];

          SceneGetCenter(G, pos);
          SceneGetView(G, view);

          if(track_camera) {
            if((diffsq3f(pos, oss->origin) > R_SMALL8) ||
               (diffsq3f(view, oss->system) > R_SMALL8) ||
               (diffsq3f(view + 4, oss->system + 3) > R_SMALL8) ||
               (diffsq3f(view + 8, oss->system + 6) > R_SMALL8)) {
              copy3f(pos, oss->origin);

              copy3f(view, oss->system);
              copy3f(view + 4, oss->system + 3);
              copy3f(view + 8, oss->system + 6);
              oss->RefreshFlag = true;
              update_flag = true;
            }
          }
          if(dynamic_grid && (!update_flag)) {
            float scale = SceneGetScreenVertexScale(G, oss->origin);

            if(fabs(scale - oss->last_scale) > R_SMALL4) {
              update_flag = true;
              oss->RefreshFlag = true;
            }
          }
        }
        if(state >= 0)
          break;
        cur_state = cur_state + 1;
        if(cur_state >= I->State.size())
          break;
      }
    }
    I->update();
  }

  ObjectPrepareContext(I, info);
  alpha = SettingGet_f(G, NULL, I->Setting.get(), cSetting_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;

  if(state >= 0)
    if(state < I->State.size())
      if(I->State[state].Active)
        oss = &I->State[state];

  while(1) {
    if(state < 0) {             /* all_states */
      oss = &I->State[cur_state];
    } else {
      if(!oss) {
        if(I->State.size() && ((SettingGetGlobal_b(G, cSetting_static_singletons) && (I->State.size() == 1))))
          oss = I->State.data();
      }
    }
    if(oss) {
      if(oss->Active) {
        if(ray) {

          ray->transparentf(1.0F - alpha);
          if((I->visRep & cRepSliceBit)) {
            float normal[3], *n0, *n1, *n2;
            int *strip = oss->strips.data();
            float *point = oss->points.data();
            float *color = oss->colors.data();
            int n = oss->n_strips;
            int a;
            int offset0 = 0, offset1 = 0, offset2, offset;
            int strip_active = false;
            int tri_count = 0;

            normal[0] = oss->system[2];
            normal[1] = oss->system[5];
            normal[2] = oss->system[8];

            n0 = normal;
            n1 = normal;
            n2 = normal;

            for(a = 0; a < n; a++) {
              offset = *(strip++);
              switch (offset) {
              case START_STRIP:
                strip_active = true;
                tri_count = 0;
                break;
              case STOP_STRIP:
                strip_active = false;
                break;
              default:
                if(strip_active) {
                  tri_count++;
                  offset2 = offset1;
                  offset1 = offset0;
                  offset0 = offset;

                  if(tri_count >= 3) {

                    if(oss->normals) {
                      n0 = oss->normals + 3 * offset0;
                      n1 = oss->normals + 3 * offset1;
                      n2 = oss->normals + 3 * offset2;
                    }

                    if(tri_count & 0x1) {       /* get the handedness right ... */
                      ray->triangle3fv(
                                        point + 3 * offset0,
                                        point + 3 * offset1,
                                        point + 3 * offset2,
                                        n0, n1, n2,
                                        color + 3 * offset0,
                                        color + 3 * offset1, color + 3 * offset2);
                    } else {
                      ray->triangle3fv(
                                        point + 3 * offset1,
                                        point + 3 * offset0,
                                        point + 3 * offset2,
                                        n1, n0, n2,
                                        color + 3 * offset1,
                                        color + 3 * offset0, color + 3 * offset2);
                    }
                  }
                }
                break;
              }
            }
          }
          ray->transparentf(0.0);
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
            if (oss->shaderCGO && (I->visRep & cRepSliceBit)){
              CGORenderGLPicking(oss->shaderCGO.get(), info, &I->context, I->Setting.get(), NULL);
            } else {
#ifndef PURE_OPENGL_ES_2
            Picking p;
	    SceneSetupGLPicking(G);
            p.context.object = I;
            p.context.state = 0;
            p.src.index = state + 1;
            p.src.bond = 0;

            if((I->visRep & cRepSliceBit)) {
              int *strip = oss->strips.data();
              float *point = oss->points.data();
              int n = oss->n_strips;
              int a;
              int offset0 = 0, offset1 = 0, offset2, offset;
              int strip_active = false;
              int tri_count = 0;
              for(a = 0; a < n; a++) {
                offset = *(strip++);
                switch (offset) {
                case START_STRIP:
                  if(!strip_active) {
                    glBegin(GL_TRIANGLES);
                  }
                  strip_active = true;
                  tri_count = 0;
                  break;
                case STOP_STRIP:
                  if(strip_active)
                    glEnd();
                  strip_active = false;
                  break;
                default:
                  if(strip_active) {
                    tri_count++;
                    offset2 = offset1;
                    offset1 = offset0;
                    offset0 = offset;

                    if(tri_count >= 3) {
                      unsigned char color[4];
                      AssignNewPickColor(nullptr, pick, color, &I->context,
                          p.src.index, p.src.bond);
                      glColor4ubv(color);

                      if(tri_count & 0x1) {     /* get the handedness right ... */
                        glVertex3fv(point + 3 * offset0);
                        glVertex3fv(point + 3 * offset1);
                        glVertex3fv(point + 3 * offset2);
                      } else {
                        glVertex3fv(point + 3 * offset1);
                        glVertex3fv(point + 3 * offset0);
                        glVertex3fv(point + 3 * offset2);
                      }
                      p.src.bond = offset0 + 1;
                    }
                  }
                  break;
                }
              }
              if(strip_active) {        /* just in case */
                glEnd();
              }
            }
#endif
            }
          } else {  // !pick
            int render_now = false;
            if(alpha > 0.0001) {
              render_now = (pass == RenderPass::Transparent);
            } else
              render_now = pass == RenderPass::Antialias;

            if(render_now) {
	      int already_rendered = false;

              if (oss->shaderCGO){
                CGORenderGL(oss->shaderCGO.get(), NULL, NULL, I->Setting.get(), info, NULL);
                already_rendered = true;
              } else {
                oss->shaderCGO.reset(CGONew(G));
              }

	      if (!already_rendered){
                SceneResetNormalCGO(G, oss->shaderCGO.get(), false);
                ObjectUseColorCGO(oss->shaderCGO.get(), I);
		
                if((I->visRep & cRepSliceBit)) {
                  int *strip = oss->strips.data();
                  float *point = oss->points.data();
                  float *color = oss->colors.data();
                  float *vnormal = oss->normals.data();
                  int n = oss->n_strips;
                  int a;
                  int offset;
                  int strip_active = false;
		  
                  {
                    float normal[3];
                    normal[0] = oss->system[2];
                    normal[1] = oss->system[5];
                    normal[2] = oss->system[8];
		    
                    CGONormalv(oss->shaderCGO.get(), normal);
                  }
		    
                  for(a = 0; a < n; a++) {
                    offset = *(strip++);
                    switch (offset) {
                    case START_STRIP:
                      if(!strip_active){
                        CGOBegin(oss->shaderCGO.get(), GL_TRIANGLE_STRIP);
                      }
                      strip_active = true;
                      break;
                    case STOP_STRIP:
                      if(strip_active)
                        CGOEnd(oss->shaderCGO.get());
                      strip_active = false;
                      break;
                    default:
                      if(strip_active) {
                        float *col;
                        col = color + 3 * offset;
                        if(vnormal)
                          CGONormalv(oss->shaderCGO.get(), vnormal + 3 * offset);
                        CGOAlpha(oss->shaderCGO.get(), alpha);
                        CGOColor(oss->shaderCGO.get(), col[0], col[1], col[2]);
                        CGOPickColor(oss->shaderCGO.get(), state + 1, offset + 1);
                        CGOVertexv(oss->shaderCGO.get(), point + 3 * offset);
                      }
                      break;
                    }
                  }
                  if(strip_active)      /* just in case */
                    CGOEnd(oss->shaderCGO.get());
		}
		
                CGOStop(oss->shaderCGO.get());
                if (use_shaders){
                  oss->shaderCGO.reset(CGOOptimizeToVBONotIndexed(oss->shaderCGO.get()));
                  assert(oss->shaderCGO->use_shader);
                }
                CGORenderGL(oss->shaderCGO.get(), NULL, NULL, I->Setting.get(), info, NULL);
                SceneInvalidatePicking(G);  // any time cgo is re-generated, needs to invalidate so
                // pick colors can be re-assigned
              }
	    }
	  }
        }
      }
    }
    if(state >= 0)
      break;
    cur_state = cur_state + 1;
    if(cur_state >= I->State.size())
      break;
  }

}


/*========================================================================*/

int ObjectSlice::getNFrame() const
{
  return State.size();
}

ObjectSliceState *ObjectSliceStateGetActive(ObjectSlice * I, int state)
{
  ObjectSliceState *ms = NULL;
  if(state >= 0) {
    if(state < I->State.size()) {
      ms = &I->State[state];
      if(!ms->Active)
        ms = NULL;
    }
  }
  return (ms);
}


/*========================================================================*/
ObjectSlice::ObjectSlice(PyMOLGlobals * G) : pymol::CObject(G)
{
  auto I = this;

  I->type = cObjectSlice;

  I->context.object = I;
  I->context.state = 0;
}

pymol::CObject* ObjectSlice::clone() const
{
  return new ObjectSlice(*this);
}
/*========================================================================*/
ObjectSlice *ObjectSliceFromMap(PyMOLGlobals * G, ObjectSlice * obj, ObjectMap * map,
                                int state, int map_state)
{
  ObjectSlice *I;
  ObjectSliceState *oss;
  ObjectMapState *oms;

  if(!obj) {
    I = new ObjectSlice(G);
  } else {
    I = obj;
  }

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, G);
  }

  oss = &I->State[state];

  oss->MapState = map_state;
  oms = ObjectMapGetState(map, map_state);
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
      if(ObjectMapStateGetExcludedStats(G, oms, NULL, 0.0F, 0.0F, tmp)) {
        oss->MapMean = tmp[1];
        oss->MapStdev = tmp[2] - tmp[1];
      } else {
        oss->MapMean = 0.0F;
        oss->MapStdev = 1.0F;
      }
    }
    /* simply copy the extents from the map -- not quite correct, but probably good enough */

    memcpy(oss->ExtentMin, oms->ExtentMin, 3 * sizeof(float));
    memcpy(oss->ExtentMax, oms->ExtentMax, 3 * sizeof(float));
    memcpy(oss->Corner, oms->Corner, 24 * sizeof(float));
  }

  strcpy(oss->MapName, map->Name);
  oss->ExtentFlag = true;

  /* set the origin of the slice to the center of the map */

  average3f(oss->ExtentMin, oss->ExtentMax, oss->origin);

  /* set the slice's system matrix to the current camera rotation matrix */

  {

    SceneViewType view;

    SceneGetView(G, view);
    copy3f(view, oss->system);
    copy3f(view + 4, oss->system + 3);
    copy3f(view + 8, oss->system + 6);
  }

  oss->RefreshFlag = true;

  if(I) {
    ObjectSliceRecomputeExtent(I);
  }

  I->ExtentFlag = true;

  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


/*========================================================================*/

void ObjectSliceRecomputeExtent(ObjectSlice * I)
{
  int extent_flag = false;
  int a;
  ObjectSliceState *ms;

  for(a = 0; a < I->State.size(); a++) {
    ms = &I->State[a];
    if(ms->Active) {
      if(ms->ExtentFlag) {
        if(!extent_flag) {
          extent_flag = true;
          copy3f(ms->ExtentMax, I->ExtentMax);
          copy3f(ms->ExtentMin, I->ExtentMin);
        } else {
          max3f(ms->ExtentMax, I->ExtentMax, I->ExtentMax);
          min3f(ms->ExtentMin, I->ExtentMin, I->ExtentMin);
        }
      }
    }
  }
  I->ExtentFlag = extent_flag;
}
