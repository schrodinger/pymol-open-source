
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

#ifndef _PYMOL_NOPY
static PyObject *ObjectVolumeStateAsPyList(ObjectVolumeState * I)
{
  PyObject *result = NULL;
  result = PyList_New(19);
  PyList_SetItem(result, 0, PyInt_FromLong(I->Active));
  PyList_SetItem(result, 1, PyString_FromString(I->MapName));
  PyList_SetItem(result, 2, PyInt_FromLong(I->MapState));
  PyList_SetItem(result, 3, CrystalAsPyList(&I->Crystal));
  PyList_SetItem(result, 4, PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->ExtentMin, 3));
  PyList_SetItem(result, 6, PConvFloatArrayToPyList(I->ExtentMax, 3));
  PyList_SetItem(result, 7, PConvIntArrayToPyList(I->Range, 6));
  PyList_SetItem(result, 8, PyFloat_FromDouble(I->Level));
  PyList_SetItem(result, 9, PyFloat_FromDouble(I->Radius));
  PyList_SetItem(result, 10, PyInt_FromLong(I->CarveFlag));
  PyList_SetItem(result, 11, PyFloat_FromDouble(I->CarveBuffer));
  if(I->CarveFlag && I->AtomVertex) {
    PyList_SetItem(result, 12, PConvFloatVLAToPyList(I->AtomVertex));
  } else {
    PyList_SetItem(result, 12, PConvAutoNone(NULL));
  }
  PyList_SetItem(result, 13, PyInt_FromLong(I->VolumeMode));
  PyList_SetItem(result, 14, PyFloat_FromDouble(I->AltLevel));
  PyList_SetItem(result, 15, PyInt_FromLong(I->quiet));
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
      if(ok)
        ok = CrystalFromPyList(&I->Crystal, PyList_GetItem(list, 3));
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 4), &I->ExtentFlag);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 5), I->ExtentMin, 3);
      if(ok)
        ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 6), I->ExtentMax, 3);
      if(ok)
        ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 7), I->Range, 6);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 8), &I->Level);
      if(ok)
        ok = PConvPyFloatToFloat(PyList_GetItem(list, 9), &I->Radius);
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
    if(vs->displayList) {
      if(PIsGlutThread()) {
        if(vs->State.G->ValidContext) {
          glDeleteLists(vs->displayList, 1);
          vs->displayList = 0;
        }
      } else {
        char buffer[255];       /* pass this off to the main thread */
        sprintf(buffer, "_cmd.gl_delete_lists(cmd._COb,%d,%d)\n", vs->displayList, 1);
        PParse(vs->State.G, buffer);
        vs->displayList = 0;
      }
    }
    for (t=0; t<2; t++) {
      if (vs->textures[t]) {
        if(PIsGlutThread()) {
          if(vs->State.G->ValidContext) {
            glDeleteTextures(1, (const GLuint *) &vs->textures[t]);
            vs->textures[t] = 0;
          }
        } else {
          char buffer[255];       /* pass this off to the main thread */
          sprintf(buffer, "_cmd.gl_delete_texture(cmd._COb,%d)\n", vs->textures[t]);
          PParse(vs->State.G, buffer);
          vs->textures[t] = 0;
        }
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
  VLAFreeP(vs->N);
  VLAFreeP(vs->V);
  FreeP(vs->VC);
  FreeP(vs->RC);
  VLAFreeP(vs->AtomVertex);
  if(vs->UnitCellCGO) {
    CGOFree(vs->UnitCellCGO);
    vs->UnitCellCGO = NULL;
  }
  if (vs->Histogram)
    free(vs->Histogram);
  if (vs->Ramp)
    free(vs->Ramp);
  vs->Active = false;
}

static void ObjectVolumeFree(ObjectVolume * I)
{
  int a;
  for(a = 0; a < I->NState; a++) {
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

void ObjectVolumeDump(ObjectVolume * I, char *fname, int state)
{
  float *v;
  int *n;
  int c;
  FILE *f;
  f = fopen(fname, "wb");
  if(!f)
    ErrMessage(I->Obj.G, "ObjectVolumeDump", "can't open file for writing");
  else {
    if(state < I->NState) {
      n = I->State[state].N;
      v = I->State[state].V;
      if(n && v)
        while(*n) {
          c = *(n++);
          if(!I->State[state].VolumeMode) {
            fprintf(f, "\n");
          }
          while(c--) {
            fprintf(f, "%10.4f%10.4f%10.4f\n", v[0], v[1], v[2]);
            v += 3;
          }
        }
    }
    fclose(f);
    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Actions)
      " ObjectVolumeDump: %s written to %s\n", I->Obj.Name, fname ENDFB(I->Obj.G);
  }
}

#if 0
static char *ObjectVolumeGetCaption(ObjectVolume * I)
{
  int state = ObjectGetCurrentState((CObject *) I, false);
  if(state >= 0) {
    if(state < I->NState) {
      char *p;
      ObjectVolumeState *vs = I->State + state;
      sprintf(vs->caption, "@ %1.4f\n", vs->Level);
      p = vs->caption + strlen(vs->caption);
      while(p > vs->caption) {
        p--;
        if(*p != '0')
          break;
      }
      return vs->caption;
    }
  }
  return NULL;
}
#endif

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

  if((rep == cRepVolume) || (rep == cRepAll) || (rep == cRepCell)) {
    for(a = 0; a < I->NState; a++) {
      if(state < 0)
        once_flag = false;
      if(!once_flag)
        state = a;
      if (level == cRepInvColor)
        I->State[state].RecolorFlag = true;
      if(level < cRepInvAll) {
        I->State[state].ResurfaceFlag = true;
        I->State[state].RecolorFlag = true;
        SceneChanged(I->Obj.G);
      } else {
        SceneInvalidate(I->Obj.G);
      }
      if(once_flag)
        break;
    }
  }
}


int ObjectVolumeGetLevel(ObjectVolume * I, int state, float *result)
{
  int ok = true;
  ObjectVolumeState *vs;
  if(state >= I->NState) {
    ok = false;
  } else {
    if(state < 0) {
      state = 0;
    }
    vs = I->State + state;
    if(vs->Active && result) {
      *result = vs->Level;
    } else
      ok = false;
  }
  return (ok);
}

int ObjectVolumeSetLevel(ObjectVolume * I, float level, int state, int quiet)
{
  int a;
  int ok = true;
  int once_flag = true;
  ObjectVolumeState *vs;
  if(state >= I->NState) {
    ok = false;
  } else {
    for(a = 0; a < I->NState; a++) {
      if(state < 0) {
        once_flag = false;
      }
      if(!once_flag) {
        state = a;
      }
      vs = I->State + state;
      if(vs->Active) {
        vs->ResurfaceFlag = true;
        vs->RefreshFlag = true;
        vs->Level = level;
        vs->quiet = quiet;
      }
      if(once_flag) {
        break;
      }
    }
  }
  return (ok);
}


static void ObjectVolumeStateUpdateColors(ObjectVolume * I, ObjectVolumeState * vs)
{
  int one_color_flag = true;
  int cur_color = SettingGet_color(I->Obj.G, I->Obj.Setting, NULL, cSetting_mesh_color);

  if(cur_color == -1)
    cur_color = I->Obj.Color;

  if(ColorCheckRamped(I->Obj.G, cur_color))
    one_color_flag = false;

  vs->OneColor = cur_color;
  if(vs->V) {
    int ramped_flag = false;
    float *v = vs->V;
    float *vc;
    int *rc;
    int a;
    int state = vs - I->State;
    int n_vert = VLAGetSize(vs->V) / 3;
    int base_n_vert = vs->base_n_V / 3;

    if(vs->VC && (vs->VCsize < n_vert)) {
      FreeP(vs->VC);
      FreeP(vs->RC);
    }

    if(!vs->VC) {
      vs->VCsize = n_vert;
      vs->VC = Alloc(float, n_vert * 3);
    }
    if(!vs->RC) {
      vs->RC = Alloc(int, n_vert);
    }
    rc = vs->RC;
    vc = vs->VC;
    if(vc) {
      for(a = 0; a < n_vert; a++) {
        if(a == base_n_vert) {
          int new_color = SettingGet_color(I->Obj.G, I->Obj.Setting,
                                           NULL, cSetting_mesh_negative_color);
          if(new_color == -1)
            new_color = cur_color;
          if(new_color != cur_color) {
            one_color_flag = false;
            cur_color = new_color;
          }
        }
        if(ColorCheckRamped(I->Obj.G, cur_color)) {
          ColorGetRamped(I->Obj.G, cur_color, v, vc, state);
          *rc = cur_color;
          ramped_flag = true;
        } else {
          float *col = ColorGet(I->Obj.G, cur_color);
          copy3f(col, vc);
        }
        rc++;
        vc += 3;
        v += 3;
      }
    }

    if(one_color_flag && (!ramped_flag)) {
      FreeP(vs->VC);
      FreeP(vs->RC);
    } else if((!ramped_flag)
              ||
              (!SettingGet_b(I->Obj.G, NULL, I->Obj.Setting, cSetting_ray_color_ramps))) {
      FreeP(vs->RC);
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
  int ok = true;
  // ObjectGadgetRamp *ogr = NULL;
  float range;

  // int mesh_skip = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_mesh_skip);

  MapType *voxelmap;            /* this has nothing to do with isosurfaces... */

  PyMOLGlobals * G = I->Obj.G;

  if (G && !(CShaderMgr_ShadersPresent(G->ShaderMgr)))
      return;

  for(a = 0; a < I->NState; a++) {
    vs = I->State + a;
    if(vs->Active) {

      map = ExecutiveFindObjectMapByName(I->Obj.G, vs->MapName);
      if(!map) {
        ok = false;
        PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Errors)
          "ObjectVolumeUpdate-Error: map '%s' has been deleted.\n", vs->MapName
          ENDFB(I->Obj.G);
        vs->ResurfaceFlag = false;
      }
      if(map) {
        oms = ObjectMapGetState(map, vs->MapState);
        if(!oms)
          ok = false;
      }
      if(oms) {
        if(vs->RefreshFlag || vs->ResurfaceFlag) {
          memcpy(vs->Corner, oms->Corner, sizeof(float) * 24);
          if(!vs->Field) {
            vs->Crystal = *(oms->Crystal);
          }

          if(I->Obj.RepVis[cRepCell]) {
            if(vs->UnitCellCGO)
              CGOFree(vs->UnitCellCGO);
            vs->UnitCellCGO = CrystalGetUnitCellCGO(&vs->Crystal);
          }

          if(oms->State.Matrix) {
            ObjectStateSetMatrix(&vs->State, oms->State.Matrix);
          } else if(vs->State.Matrix) {
            ObjectStateResetMatrix(&vs->State);
          }
          vs->RefreshFlag = false;
          vs->displayListInvalid = true;
        }
        if (!vs->Histogram) {
          vs->Histogram = malloc(sizeof(float) * (HISTOGRAM_SIZE + 4)); 
        }
        range =  SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_data_range);
        ObjectMapStateGetHistogram(I->Obj.G, oms, HISTOGRAM_SIZE, range, vs->Histogram);
      }

      if(map && oms && vs->N && vs->V && I->Obj.RepVis[cRepVolume]) {
        if(vs->ResurfaceFlag) {
          Isofield *field = NULL;
          vs->RecolorFlag = true;
          vs->ResurfaceFlag = false;
          if(!vs->quiet) {
          }
          if(vs->Field) {
            field = vs->Field;
          } else if(oms->Field) {
            field = oms->Field;
          }

          if(field) {
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

              IsosurfGetRange(I->Obj.G, field, oms->Crystal,
                              min_ext, max_ext, vs->Range, true);
            }

            if (vs->volume)
                FieldFree(vs->volume);

            vs->volume = FieldNewCopy(I->Obj.G, field->data);

          }

          if(vs->CarveFlag && vs->AtomVertex) {

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
        if(vs->RecolorFlag) {
          ObjectVolumeStateUpdateColors(I, vs);
          //vs->RecolorFlag = false;
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

int ObjectVolumeColor(ObjectVolume * I, float * colors, int ncolors) {

  int i;
  PyMOLGlobals * G = I->Obj.G;
  ObjectVolumeState * vs = NULL;

  PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeColor-Update: Coloring volume with %d colors.\n", ncolors
    ENDFB(I->Obj.G);  

  if(!G) {
    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeColor-Error: Invalid Globals.\n"
    ENDFB(I->Obj.G);  
    return false;
  }
  
  if(I->NState &&
     ((SettingGet(I->Obj.G, cSetting_static_singletons) && (I->NState == 1))))
    vs = I->State;

  PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeColor-Msg: Found VolumeState.\n"
    ENDFB(I->Obj.G); 

  if (vs && vs->colors) {
    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
      "ObjectVolumeColor-Msg: Found colors, freeing.\n"
      ENDFB(I->Obj.G); 

    free(vs->colors);
    vs->colors = NULL;
  }

    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeColor-Msg: Freed Colors.\n"
    ENDFB(I->Obj.G); 
  /* ncolors actually counts the length of the array, not 
   * real number of RGBA colors */
    vs->colors = (float*) malloc(ncolors * sizeof(float));

    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
      "ObjectVolumeColor-Msg: Alloc'd Colors.\n"
      ENDFB(I->Obj.G); 

    if(vs->colors) {
      for(i=0;i<ncolors;i++) {
	vs->colors[i] = (float) colors[i];
      }

      PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
	"ObjectVolumeColor-Update: Successfully copied the colors.\n"
	ENDFB(I->Obj.G);  
    } else {
      PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
	"ObjectVolumeColor-Update: Copying colors failed--no memory.\n"
	ENDFB(I->Obj.G);  
      /* manual coloring failed; use default */
      I->Obj.Color = 0;
      return false;
    }

    /* setting this to 1 uses the users' colors */
    /* setting this to 0 uses default */
    I->Obj.Color = 1; 
    
    PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
      "ObjectVolumeColor-Msg: Exiting VolumeColor.\n"
      ENDFB(I->Obj.G); 

  return true;
}


static void ObjectVolumeRender(ObjectVolume * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->Obj.G;
  float *v = NULL;
  int state = info->state;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int pass = info->pass;
  int *n = NULL;
  int a = 0;
  ObjectVolumeState *vs = NULL;
  int volume_bit_depth =  (int) SettingGet_i(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_bit_depth);
  float volume_layers =  SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_volume_layers);
  /* make this a setting? */
  int volume_nColors = 360;
  GLint alpha_func;
  GLfloat alpha_ref;
  float min_val, max_val, scale, bias;
  float tex_corner[24];
  float x_offset, y_offset, z_offset;
  float *min_ext, *max_ext, *corner, *m;
  float tmp_min[3], tmp_max[3], zaxis[3];
  float points[36], tex_coords[36];
  int n_points;
  float d, sliceRange, sliceDelta;
  float origin[3];
  CShaderPrg * p;
  float *fog_color;
  float fog_enabled;

  /* bail if no shaders */
  if (G && !(CShaderMgr_ShadersPresent(G->ShaderMgr)))
      return;

  if (volume_bit_depth==4)
    volume_bit_depth = GL_RGBA4;
  else if (volume_bit_depth==8)
    volume_bit_depth = GL_RGBA;
  else if (volume_bit_depth==16)
    volume_bit_depth = GL_RGBA16F_ARB;
  else if (volume_bit_depth==32)
    volume_bit_depth = GL_RGBA32F_ARB;
  /* default back to 16-bit for improper value */
  else
    volume_bit_depth = GL_RGBA16F_ARB;

  ObjectPrepareContext(&I->Obj, ray);

  if(state >= 0)
    if(state < I->NState)
      if(I->State[state].Active)
        if(I->State[state].V && I->State[state].N)
          vs = I->State + state;

  PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
    "ObjectVolumeRender-Msg: Found VolumeState.\n"
    ENDFB(I->Obj.G); 

  while(1) {
    if(state < 0) {             /* all_states */
      vs = I->State + a;
    } else {
      if(!vs) {
        if(I->NState &&
           ((SettingGet(I->Obj.G, cSetting_static_singletons) && (I->NState == 1))))
          vs = I->State;
      }
    }

    if(vs) {
      if(vs->Active && vs->V && vs->N) {

	PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
	  "ObjectVolumeRender-Msg: vs->Active and V and N\n"
	  ENDFB(I->Obj.G);

        v = vs->V;
        n = vs->N;
        if(G->HaveGUI && G->ValidContext) {
          if(pick) {
          } else {
            if(pass == -1) {
              int use_dlst;

PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeRender-Msg: !pass=%d and refresh=%d and recolor=%d.\n", 
  !pass, vs->RefreshFlag, vs->RecolorFlag
	ENDFB(I->Obj.G); 
		
              if(!info->line_lighting)
                glDisable(GL_LIGHTING);

              SceneResetNormal(I->Obj.G, false);
              /*ObjectUseColor(&I->Obj);*/

              // Never use display lists for volume rendering
              use_dlst = 0; 

              if(use_dlst && vs->displayList && vs->displayListInvalid) {
                glDeleteLists(vs->displayList, 1);
                vs->displayList = 0;
                vs->displayListInvalid = false;
              }

              if (!vs->textures[0] || vs->RefreshFlag) {

PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeRender-Msg: texture not inited or we need refresh.\n"
  ENDFB(I->Obj.G); 

                vs->RefreshFlag = false;

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
		/* FIX ME
		 * This is non-ideal, but works (for now). */
		glTexImage3D = getTexImage3D(); 
		glActiveTexture = getActiveTexture();

		if (! glActiveTexture) {
		  printf("Could not bind the glActiveTexture function.  No volume shaders.\n");
		  abort();
		} 
		else 
		  printf("glActiveTexture is: %x", glActiveTexture);
		
		if (! glTexImage3D) {
		  printf("Could not bind the glTexImage3D function.  No volume shaders.\n");
		  abort();
		} 
		else
		  printf("glTexImage3D is: %x", glTexImage3D);
#endif
/* END PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  
PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeRender-Msg: About to create the texture-unit.\n"
  ENDFB(I->Obj.G); 
/*
                // For future experiments with tex coord jittering

                if (vs->textures[2]) {
                  glDeleteTextures(1, &vs->textures[2]);
                  vs->textures[2] = 0;
                }
                unsigned char random_data[3*1024];
                for (i=0;i<3*1024;i++) {
                    random_data[i] = (rand() % 255);
                }
                glGenTextures(1, &vs->textures[2]);
                glBindTexture(GL_TEXTURE_2D, vs->textures[2]);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
                glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 32, 32, 0, GL_RGB, GL_UNSIGNED_BYTE, random_data);
*/
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

                min_val = vs->Histogram[0];
                max_val = vs->Histogram[1];

                scale = 1.0F/(max_val-min_val);
                bias = -min_val*scale;
                
                glPixelTransferf(GL_RED_BIAS,  bias);
                glPixelTransferf(GL_RED_SCALE, scale);
    		        glTexImage3D(GL_TEXTURE_3D, 0, volume_bit_depth,
		      	      vs->volume->dim[2], vs->volume->dim[1], vs->volume->dim[0], 0, 
                  GL_RED, GL_FLOAT, vs->volume->data);
                glPixelTransferf(GL_RED_BIAS, 0.0);
                glPixelTransferf(GL_RED_SCALE, 1.0);
PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeRender-Msg: Created.\n"
  ENDFB(I->Obj.G); 
                // Create a color map - this will be replaced by a proper Volume Ramp object
                if (vs->textures[1]) {
                  glDeleteTextures(1, (const GLuint *) &vs->textures[1]);
                  vs->textures[1] = 0;
                }
              }

          		if (vs->RecolorFlag || 0==I->Obj.Color) {
          		  int i;
PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeColor-Update: Making new color ramp.\n"
  ENDFB(I->Obj.G);
          		  /* default coloring */
          		  vs->RecolorFlag = false;
          		  if (0==I->Obj.Color) {
PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeColor-Update: No colors, making new.\n"
  ENDFB(I->Obj.G);
		    
                  /* don't need to create more space, just overwrite the values */
            		  if (vs->colors) free(vs->colors);
            		  vs->colors = calloc(volume_nColors, sizeof(float)*4);

          		    /* default rep? */
          		    for (i=0; i<volume_nColors; i++) {
            	      float norm;
          		      vs->colors[4*i+1] = (1./((float) volume_nColors-1) )*(1.0-i);
		                vs->colors[4*i] = 0.25;
          		      vs->colors[4*i+2] = 1.0-(1./((float) volume_nColors-1))*(i);
		                norm = (i-(volume_nColors/2.))/(volume_nColors/2.0);
          		      vs->colors[4*i+3] = 0.1*exp(-(norm+0.5)*(norm+0.5)) + 0.5*(exp(-(norm-1.1)*(norm-1.1)/0.4));
		              }
		            } else {
PRINTFB(I->Obj.G, FB_ObjectVolume, FB_Blather)
  "ObjectVolumeColor-Update: Using established color ramp.\n"
  ENDFB(I->Obj.G);  
            		}
			  glGenTextures(1, (GLuint *) &vs->textures[1]);
                glBindTexture(GL_TEXTURE_1D, vs->textures[1]);
                glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
                glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
                glTexImage1D(GL_TEXTURE_1D, 0, volume_bit_depth, volume_nColors, 0, GL_RGBA, GL_FLOAT, vs->colors);
              } // recolor

              if(use_dlst && vs->displayList) {
                glCallList(vs->displayList);
              } else {
                if(use_dlst) {
                  if(!vs->displayList) {
                    vs->displayList = glGenLists(1);
                    if(vs->displayList) {
                      glNewList(vs->displayList, GL_COMPILE_AND_EXECUTE);
                    }
                  }
                }
                if(n && v && I->Obj.RepVis[cRepVolume]) {
                  x_offset = 0.5 / vs->volume->dim[2];
                  y_offset = 0.5 / vs->volume->dim[1];
                  z_offset = 0.5 / vs->volume->dim[0];

                  tex_corner[0] = x_offset;
                  tex_corner[1] = y_offset;
                  tex_corner[2] = z_offset;
                  tex_corner[3] = x_offset;
                  tex_corner[4] = y_offset;
                  tex_corner[5] = 1.0-z_offset;
                  tex_corner[6] = x_offset;
                  tex_corner[7] = 1.0-y_offset;
                  tex_corner[8] = z_offset;
                  tex_corner[9] = x_offset;
                  tex_corner[10] = 1.0-y_offset;
                  tex_corner[11] = 1.0-z_offset;
                  tex_corner[12] = 1.0-x_offset;
                  tex_corner[13] = y_offset;
                  tex_corner[14] = z_offset;
                  tex_corner[15] = 1.0-x_offset;
                  tex_corner[16] = y_offset;
                  tex_corner[17] = 1.0-z_offset;
                  tex_corner[18] = 1.0-x_offset;
                  tex_corner[19] = 1.0-y_offset;
                  tex_corner[20] = z_offset;
                  tex_corner[21] = 1.0-x_offset;
                  tex_corner[22] = 1.0-y_offset;
                  tex_corner[23] = 1.0-z_offset;

                  p = CShaderMgr_GetShaderPrg(G->ShaderMgr, "volume");

                  corner = vs->Corner;

                  origin[0] = corner[0] + 0.5*(corner[21]-corner[0]);
		  origin[1] = corner[1] + 0.5*(corner[22]-corner[1]);
		  origin[2] = corner[2] + 0.5*(corner[23]-corner[2]);

                  if(MatrixInvTransformExtentsR44d3f(vs->State.Matrix,
                                         vs->ExtentMin, vs->ExtentMax,
                                         tmp_min, tmp_max)) {
                    min_ext = tmp_min;
                    max_ext = tmp_max;
                  } else {
                    min_ext = vs->ExtentMin;
                    max_ext = vs->ExtentMax;
                  }
                  glDisable(GL_LIGHTING);
		  /* extents */
		  if(I->Obj.RepVis[cRepExtent] && corner) {
		    glBegin(GL_LINES);
		    glColor3f(1,1,1);
		    glVertex3f(corner[0], corner[1], corner[2]);
		    glVertex3f(corner[3], corner[4], corner[5]);
		    glVertex3f(corner[3], corner[4], corner[5]);
		    glVertex3f(corner[9], corner[10], corner[11]);
		    glVertex3f(corner[9], corner[10], corner[11]);
		    glVertex3f(corner[6], corner[7], corner[8]); 
		    glVertex3f(corner[6], corner[7], corner[8]);
		    glVertex3f(corner[0], corner[1], corner[2]); 

		    glVertex3f(corner[12], corner[13], corner[14]);
		    glVertex3f(corner[15], corner[16], corner[17]);
		    glVertex3f(corner[15], corner[16], corner[17]);
		    glVertex3f(corner[21], corner[22], corner[23]);
		    glVertex3f(corner[21], corner[22], corner[23]);
		    glVertex3f(corner[18], corner[19], corner[20]);
		    glVertex3f(corner[18], corner[19], corner[20]); 
		    glVertex3f(corner[12], corner[13], corner[14]); 

		    glVertex3f(corner[0], corner[1], corner[2]);
		    glVertex3f(corner[12], corner[13], corner[14]);
		    glVertex3f(corner[3], corner[4], corner[5]);
		    glVertex3f(corner[15], corner[16], corner[17]);
		    glVertex3f(corner[9], corner[10], corner[11]);
		    glVertex3f(corner[21], corner[22], corner[23]);
		    glVertex3f(corner[6], corner[7], corner[8]); 
		    glVertex3f(corner[18], corner[19], corner[20]);
		    glEnd();
		  }

                  m = SceneGetMatrix(G);
                  zaxis[0] = m[2];
		  zaxis[1] = m[6];
		  zaxis[2] = m[10];

                  // determine number of slices based on max extent
                  // and slice option
                  sliceRange = 0.5*sqrt(2.0) * 
                    fmax(fmax(fabs(corner[21]-corner[0]), fabs(corner[22]-corner[1])), 
                         fabs(corner[23]-corner[2]));
		  sliceDelta = (sliceRange / volume_layers);


/* // DEBUG BOUNDING CUBE */
/* glColor3f(1,0,1); */
/* glBegin(GL_LINES); */
/* glVertex3f(origin[0]-sliceRange,origin[1]-sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]+sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]+sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]+sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]+sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]-sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]-sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]-sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]-sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]+sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]+sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]+sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]+sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]-sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]-sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]-sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]-sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]-sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]+sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]-sliceRange,origin[1]+sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]+sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]+sliceRange,origin[2]+sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]-sliceRange,origin[2]-sliceRange); */
/* glVertex3f(origin[0]+sliceRange,origin[1]-sliceRange,origin[2]+sliceRange); */
/* glEnd(); */
/* // */

                CShaderPrg_Enable(p);
                CShaderPrg_Set1i(p, "volumeTex", 0);
                CShaderPrg_Set1i(p, "colorTex", 1);
              	CShaderPrg_Set1f(p, "nSlices", volume_layers);

                fog_color = SettingGetfv(G, cSetting_bg_rgb);                
                fog_enabled = SettingGet(G, cSetting_depth_cue) ? 1.0 : 0.0;
           
                CShaderPrg_Set1f(p, "fog_r", fog_color[0]);
                CShaderPrg_Set1f(p, "fog_g", fog_color[1]);
                CShaderPrg_Set1f(p, "fog_b", fog_color[2]);
                CShaderPrg_Set1f(p, "fog_enabled", fog_enabled);

                glActiveTexture(GL_TEXTURE2);
                glBindTexture(GL_TEXTURE_2D, vs->textures[2]);
                glActiveTexture(GL_TEXTURE1);
                glBindTexture(GL_TEXTURE_1D, vs->textures[1]);
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_3D, vs->textures[0]);

                // PyMOL uses different glAlphaFunct and we need to reset it
                // to default value (everything passes)

                // Not sure if we really need to restore this
                glGetIntegerv(GL_ALPHA_TEST_FUNC, &alpha_func);
                glGetFloatv(GL_ALPHA_TEST_REF, &alpha_ref);

                glAlphaFunc(GL_ALWAYS, 0.0);

                // This is setting used for PyMOL, but just to be on a safe side
                // we set glBlendFunct explicitely here
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                for (d=sliceRange; d>=-sliceRange; d -= sliceDelta) {

                    n_points = 0;
                    // Slice the volume
                    n_points += ObjectVolumeAddSlicePoint(&corner[0],&corner[3],zaxis,d,&points[n_points], 
                      &tex_corner[0], &tex_corner[3], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[3],&corner[9],zaxis,d,&points[n_points], 
                      &tex_corner[3], &tex_corner[9], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[9],&corner[6],zaxis,d,&points[n_points], 
                      &tex_corner[9], &tex_corner[6], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[6],&corner[0],zaxis,d,&points[n_points], 
                      &tex_corner[6], &tex_corner[0], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[12],&corner[15],zaxis,d,&points[n_points], 
                      &tex_corner[12], &tex_corner[15], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[15],&corner[21],zaxis,d,&points[n_points], 
                      &tex_corner[15], &tex_corner[21], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[21],&corner[18],zaxis,d,&points[n_points], 
                      &tex_corner[21], &tex_corner[18], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[18],&corner[12],zaxis,d,&points[n_points], 
                      &tex_corner[18], &tex_corner[12], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[0],&corner[12],zaxis,d,&points[n_points], 
                      &tex_corner[0], &tex_corner[12], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[3],&corner[15],zaxis,d,&points[n_points], 
                      &tex_corner[3], &tex_corner[15], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[9],&corner[21],zaxis,d,&points[n_points], 
                      &tex_corner[9], &tex_corner[21], &tex_coords[n_points], origin);
                    n_points += ObjectVolumeAddSlicePoint(&corner[6],&corner[18],zaxis,d,&points[n_points], 
                      &tex_corner[6], &tex_corner[18], &tex_coords[n_points], origin);
                    ObjectVolumeDrawSlice(points, tex_coords, n_points/3, zaxis);
                  }
                  CShaderPrg_Disable(p);
                  glAlphaFunc(alpha_func, alpha_ref);
                  glEnable(GL_LIGHTING);
                }
                if(use_dlst && vs->displayList) {
                  glEndList();
                }

              }
              glEnable(GL_LIGHTING);
            }
          }
        }
      }
    }
    if(state >= 0)
      break;
    a = a + 1;
    if(a >= I->NState)
      break;
  }
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
    glBegin(GL_POLYGON);
    for (i=0; i<n_points; i++) {
      glTexCoord3fv(&tex_coords[3*vertices[(i)%n_points]]);
      glVertex3fv(&points[3*vertices[(i)%n_points]]);
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
  I->State = VLAMalloc(10, sizeof(ObjectVolumeState), 5, true);   /* autozero important */

  I->Obj.type = cObjectVolume;
  I->Obj.Color = 0;  /* 0 = auto; 1 = user */

  I->Obj.fFree = (void (*)(CObject *)) ObjectVolumeFree;
  I->Obj.fUpdate = (void (*)(CObject *)) ObjectVolumeUpdate;
  I->Obj.fRender = (void (*)(CObject *, RenderInfo *)) ObjectVolumeRender;
  I->Obj.fInvalidate = (void (*)(CObject *, int, int, int)) ObjectVolumeInvalidate;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectVolumeGetNStates;
  /*  I->Obj.fGetCaption = (char *(*)(CObject *))ObjectVolumeGetCaption; */
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
  if(!vs->V) {
    vs->V = VLAlloc(float, 10000);
  }
  if(!vs->N) {
    vs->N = VLAlloc(int, 10000);
  }
  if(vs->AtomVertex) {
    VLAFreeP(vs->AtomVertex);
  }
  vs->N[0] = 0;
  vs->Active = true;
  vs->ResurfaceFlag = true;
  vs->RecolorFlag = false;
  vs->ExtentFlag = false;
  vs->CarveFlag = false;
  vs->quiet = true;
  vs->CarveBuffer = 0.0;
  vs->AtomVertex = NULL;
  vs->UnitCellCGO = NULL;
  vs->displayList = 0;
  vs->displayListInvalid = true;
  vs->caption[0] = 0;
  vs->Field = NULL;
  vs->textures[0] = 0;
  vs->textures[1] = 0;
  vs->textures[2] = 0;
//  vs->Histogram = NULL;
  vs->Histogram = (float*)calloc(sizeof(float), HISTOGRAM_SIZE+4);
  vs->isUpdated = false;
  // Initial ramp
  vs->RampSize = 5;
  vs->Ramp = (float*)malloc(5*5*sizeof(float));
  vs->Ramp[0] = 0.0; vs->Ramp[1] = 0.0; vs->Ramp[2] = 0.0; vs->Ramp[3] = 1.0; vs->Ramp[4] = 0.0;
  vs->Ramp[5] = 200.0; vs->Ramp[6] = 0.0; vs->Ramp[7] = 0.0; vs->Ramp[8] = 1.0; vs->Ramp[9] = 0.0;
  vs->Ramp[10] = 210.0; vs->Ramp[11] = 1.0; vs->Ramp[12] = 0.0; vs->Ramp[13] = 0.2; vs->Ramp[14] = 0.2;
  vs->Ramp[15] = 220.0; vs->Ramp[16] = 0.0; vs->Ramp[17] = 0.0; vs->Ramp[18] = 1.0; vs->Ramp[19] = 0.0;
  vs->Ramp[20] = 359.0; vs->Ramp[21] = 0.0; vs->Ramp[22] = 0.0; vs->Ramp[23] = 1.0; vs->Ramp[24] = 0.0;
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

  if(!obj) {
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

  vs->Level = level;
  vs->AltLevel = alt_level;
  vs->VolumeMode = meshMode;
  vs->quiet = quiet;
  if(oms) {

//    if((meshMode == 3) && (vs->AltLevel < vs->Level)) {
      /* gradient object -- need to auto-set range */
      if(!ObjectMapStateGetDataRange(G, oms, &vs->Level, &vs->AltLevel)) {
        vs->Level = -1.0F;
        vs->AltLevel = 1.0F;
      }
//    }

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

      {
        int eff_range[6];

        if(IsosurfGetRange
           (G, oms->Field, oms->Crystal, min_ext, max_ext, eff_range, false)) {
          int fdim[3];
          int expand_result;
          /* need to generate symmetry-expanded temporary map */

          vs->Crystal = *(oms->Crystal);
          fdim[0] = eff_range[3] - eff_range[0];
          fdim[1] = eff_range[4] - eff_range[1];
          fdim[2] = eff_range[5] - eff_range[2];
          vs->Field = IsosurfFieldAlloc(I->Obj.G, fdim);

          expand_result =
            IsosurfExpand(oms->Field, vs->Field, oms->Crystal, sym, eff_range);

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

          vs->Range[0] = 0;
          vs->Range[1] = 0;
          vs->Range[2] = 0;
          vs->Range[3] = fdim[0];
          vs->Range[4] = fdim[1];
          vs->Range[5] = fdim[2];

        } else {
          /* mesh entirely contained within bounds of current map */
          int a;
          for(a = 0; a < 6; a++) {
            vs->Range[a] = eff_range[a];
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
    if(I) {
      ObjectVolumeRecomputeExtent(I);
    }
    I->Obj.ExtentFlag = true;
    /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]); */
  }
  if(!ok) {
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
  ObjectVolume *I;
  ObjectVolumeState *vs;
  ObjectMapState *oms;

  if(!obj) {
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

  vs->Level = level;
  vs->AltLevel = alt_level;
  vs->VolumeMode = meshMode;
  vs->quiet = quiet;
  if(oms) {

    if((meshMode == 3) && (vs->AltLevel < vs->Level)) {
      /* gradient object -- need to auto-set range */
      if(!ObjectMapStateGetDataRange(G, oms, &vs->Level, &vs->AltLevel)) {
        vs->Level = -1.0F;
        vs->AltLevel = 1.0F;
      }
    }

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

      IsosurfGetRange(G, oms->Field, oms->Crystal, min_ext, max_ext, vs->Range, true);
    }
    vs->ExtentFlag = true;
  }
  if(carve != 0.0) {
    vs->CarveFlag = true;
    vs->CarveBuffer = carve;
    vs->AtomVertex = vert_vla;
  }
  if(I) {
    ObjectVolumeRecomputeExtent(I);
  }
  I->Obj.ExtentFlag = true;
  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]); */
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
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
  int a;

  if (!I) return NULL;

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetField Entering" ENDFD;

  for(a = 0; a < I->NState; a++) 
    if(I->State[a].Active) {
      ovs = I->State + a;
      if (I) {
        F = ovs->volume;
        n_elem = F->size / F->base_size;
        result = PConvFloatArrayToPyList((float *) F->data, n_elem);
      }
      
      break;
    }

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetField Exiting" ENDFD;

  return (PConvAutoNone(result));
#endif
}

/*==============================================================================*/
PyObject * ObjectVolumeGetHistogram(ObjectVolume * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* TODO: Allow for multi-state maps? */
  int n_elem;
  PyObject * result = NULL;
  ObjectVolumeState *ovs;
  int a;

  if (!I) return NULL;

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetHistogram Entering" ENDFD;

  for(a = 0; a < I->NState; a++)
    if(I->State[a].Active) {
      ovs = I->State + a;
      if (I) {
        n_elem = HISTOGRAM_SIZE + 4;
        result = PConvFloatArrayToPyList(ovs->Histogram, n_elem);
      }
      break;
    }

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetHistogram Exiting" ENDFD;

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
  int n_elem;
  PyObject * result = NULL;
  ObjectVolumeState *ovs;
  int a;

  if (!I) return NULL;

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetRamp Entering" ENDFD;

  for(a = 0; a < I->NState; a++)
    if(I->State[a].Active) {
      ovs = I->State + a;  
      if (I) { 
        n_elem = 5 * ovs->RampSize;
        result = PConvFloatArrayToPyList(ovs->Ramp, n_elem);
      }
      break;
    }

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetRamp Exiting" ENDFD;
  
  return (PConvAutoNone(result));
#endif
}

/*==============================================================================*/
PyObject * ObjectVolumeSetRamp(ObjectVolume * I, float *ramp_list, int list_size)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* TODO: Allow for multi-state maps? */
  PyObject * result = NULL;
  ObjectVolumeState *ovs;
  int a;

  if (!I) return NULL;

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-SetRamp Entering" ENDFD;

  for(a = 0; a < I->NState; a++) 
    if(I->State[a].Active) {
      ovs = I->State + a;
      if (ramp_list && list_size>0) {
        if (ovs->Ramp) {
          free(ovs->Ramp);
        }
        ovs->Ramp = ramp_list;
        ovs->RampSize = list_size / 5;
      }
      result = PyInt_FromLong(list_size);
      break; 
   }
  
  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-SetRamp Exiting" ENDFD;
  
  return (PConvAutoNone(result));
#endif
}

/*==============================================================================*/
PyObject * ObjectVolumeGetIsUpdated(ObjectVolume * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  /* TODO: Allow for multi-state maps? */
  PyObject * result = NULL;
  ObjectVolumeState *ovs;
  int a;
  
  if (!I) return NULL;
  
  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetIsUpdated Entering" ENDFD;
  
  for(a = 0; a < I->NState; a++)
    if(I->State[a].Active) {
      ovs = I->State + a;  
      if (I) { 
        result = PyInt_FromLong(ovs->isUpdated);
      }
      break;
    }

  PRINTFD(I->Obj.G,FB_ObjectVolume) "ObjectVolume-GetIsUpdated Exiting" ENDFD;
   
  return (PConvAutoNone(result));
#endif
}

