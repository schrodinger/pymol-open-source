
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
#include"Executive.h"
#include"PConv.h"
#include"P.h"
#include"Matrix.h"

ObjectMesh *ObjectMeshNew(PyMOLGlobals * G);

static void ObjectMeshFree(ObjectMesh * I);
static void ObjectMeshInvalidate(ObjectMesh * I, int rep, int level, int state);
void ObjectMeshStateInit(PyMOLGlobals * G, ObjectMeshState * ms);
void ObjectMeshRecomputeExtent(ObjectMesh * I);

#ifndef _PYMOL_NOPY
static PyObject *ObjectMeshStateAsPyList(ObjectMeshState * I)
{
  PyObject *result = NULL;

  result = PyList_New(17);

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
  PyList_SetItem(result, 13, PyInt_FromLong(I->MeshMode));
  PyList_SetItem(result, 14, PyFloat_FromDouble(I->AltLevel));
  PyList_SetItem(result, 15, PyInt_FromLong(I->quiet));
  if(I->Field) {
    PyList_SetItem(result, 16, IsosurfAsPyList(I->Field));
  } else {
    PyList_SetItem(result, 16, PConvAutoNone(NULL));
  }
  return (PConvAutoNone(result));
}
#endif

#ifndef _PYMOL_NOPY
static PyObject *ObjectMeshAllStatesAsPyList(ObjectMesh * I)
{

  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NState);
  for(a = 0; a < I->NState; a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result, a, ObjectMeshStateAsPyList(I->State + a));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(NULL));
    }
  }
  return (PConvAutoNone(result));
}
#endif

#ifndef _PYMOL_NOPY
static int ObjectMeshStateFromPyList(PyMOLGlobals * G, ObjectMeshState * I,
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
      ObjectMeshStateInit(G, I);
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
        ok = PConvPyIntToInt(PyList_GetItem(list, 13), &I->MeshMode);
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
    }
  }
  return (ok);
}
#endif

#ifndef _PYMOL_NOPY
static int ObjectMeshAllStatesFromPyList(ObjectMesh * I, PyObject * list)
{

  int ok = true;
  int a;
  VLACheck(I->State, ObjectMeshState, I->NState);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    for(a = 0; a < I->NState; a++) {
      ok = ObjectMeshStateFromPyList(I->Obj.G, I->State + a, PyList_GetItem(list, a));
      if(!ok)
        break;
    }
  }
  return (ok);
}
#endif

int ObjectMeshNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectMesh ** result)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  int ll;
  ObjectMesh *I = NULL;
  (*result) = NULL;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  I = ObjectMeshNew(G);
  if(ok)
    ok = (I != NULL);

  if(ok)
    ok = ObjectFromPyList(G, PyList_GetItem(list, 0), &I->Obj);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->NState);
  if(ok)
    ok = ObjectMeshAllStatesFromPyList(I, PyList_GetItem(list, 2));
  if(ok) {
    (*result) = I;
    ObjectMeshRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
#endif
}

PyObject *ObjectMeshAsPyList(ObjectMesh * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(&I->Obj));
  PyList_SetItem(result, 1, PyInt_FromLong(I->NState));
  PyList_SetItem(result, 2, ObjectMeshAllStatesAsPyList(I));

  return (PConvAutoNone(result));
#endif
}

static void ObjectMeshStateFree(ObjectMeshState * ms)
{
  ObjectStatePurge(&ms->State);
  if(ms->State.G->HaveGUI) {
    if(ms->displayList) {
      if(PIsGlutThread()) {
        if(ms->State.G->ValidContext) {
          glDeleteLists(ms->displayList, 1);
          ms->displayList = 0;
        }
      } else {
        char buffer[255];       /* pass this off to the main thread */
        sprintf(buffer, "_cmd.gl_delete_lists(cmd._COb,%d,%d)\n", ms->displayList, 1);
        PParse(ms->State.G, buffer);
        ms->displayList = 0;
      }
    }
  }
  if(ms->Field) {
    IsosurfFieldFree(ms->State.G, ms->Field);
    ms->Field = NULL;
  }
  VLAFreeP(ms->N);
  VLAFreeP(ms->V);
  FreeP(ms->VC);
  FreeP(ms->RC);
  VLAFreeP(ms->AtomVertex);
  if(ms->UnitCellCGO) {
    CGOFree(ms->UnitCellCGO);
    ms->UnitCellCGO = NULL;
  }
  ms->Active = false;
}

static void ObjectMeshFree(ObjectMesh * I)
{
  int a;
  for(a = 0; a < I->NState; a++) {
    if(I->State[a].Active)
      ObjectMeshStateFree(I->State + a);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);

  OOFreeP(I);
}

int ObjectMeshInvalidateMapName(ObjectMesh * I, char *name)
{
  int a;
  ObjectMeshState *ms;
  int result = false;
  for(a = 0; a < I->NState; a++) {
    ms = I->State + a;
    if(ms->Active) {
      if(strcmp(ms->MapName, name) == 0) {
        ObjectMeshInvalidate(I, cRepAll, cRepInvAll, a);
        result = true;
      }
    }
  }
  return result;
}

void ObjectMeshDump(ObjectMesh * I, char *fname, int state)
{
  float *v;
  int *n;
  int c;
  FILE *f;
  f = fopen(fname, "wb");
  if(!f)
    ErrMessage(I->Obj.G, "ObjectMeshDump", "can't open file for writing");
  else {
    if(state < I->NState) {
      n = I->State[state].N;
      v = I->State[state].V;
      if(n && v)
        while(*n) {
          c = *(n++);
          if(!I->State[state].MeshMode) {
            fprintf(f, "\n");
          }
          while(c--) {
            fprintf(f, "%10.4f%10.4f%10.4f\n", v[0], v[1], v[2]);
            v += 3;
          }
        }
    }
    fclose(f);
    PRINTFB(I->Obj.G, FB_ObjectMesh, FB_Actions)
      " ObjectMeshDump: %s written to %s\n", I->Obj.Name, fname ENDFB(I->Obj.G);
  }
}

#if 0
static char *ObjectMeshGetCaption(ObjectMesh * I)
{
  int state = ObjectGetCurrentState((CObject *) I, false);
  if(state >= 0) {
    if(state < I->NState) {
      char *p;
      ObjectMeshState *ms = I->State + state;
      sprintf(ms->caption, "@ %1.4f\n", ms->Level);
      p = ms->caption + strlen(ms->caption);
      while(p > ms->caption) {
        p--;
        if(*p != '0')
          break;
      }
      return ms->caption;
    }
  }
  return NULL;
}
#endif

static void ObjectMeshInvalidate(ObjectMesh * I, int rep, int level, int state)
{
  int a;
  int once_flag = true;
  if(level >= cRepInvExtents) {
    I->Obj.ExtentFlag = false;
  }
  if((rep == cRepMesh) || (rep == cRepAll) || (rep == cRepCell)) {
    for(a = 0; a < I->NState; a++) {
      if(state < 0)
        once_flag = false;
      if(!once_flag)
        state = a;
      I->State[state].RefreshFlag = true;
      if(level >= cRepInvAll) {
        I->State[state].ResurfaceFlag = true;
        SceneChanged(I->Obj.G);
      } else if(level >= cRepInvColor) {
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

int ObjectMeshGetLevel(ObjectMesh * I, int state, float *result)
{
  int ok = true;
  ObjectMeshState *ms;
  if(state >= I->NState) {
    ok = false;
  } else {
    if(state < 0) {
      state = 0;
    }
    ms = I->State + state;
    if(ms->Active && result) {
      *result = ms->Level;
    } else
      ok = false;
  }
  return (ok);
}

int ObjectMeshSetLevel(ObjectMesh * I, float level, int state, int quiet)
{
  int a;
  int ok = true;
  int once_flag = true;
  ObjectMeshState *ms;
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
      ms = I->State + state;
      if(ms->Active) {
        ms->ResurfaceFlag = true;
        ms->RefreshFlag = true;
        ms->Level = level;
        ms->quiet = quiet;
      }
      if(once_flag) {
        break;
      }
    }
  }
  return (ok);
}

static void ObjectMeshStateUpdateColors(ObjectMesh * I, ObjectMeshState * ms)
{
  int one_color_flag = true;
  int cur_color = SettingGet_color(I->Obj.G, I->Obj.Setting, NULL, cSetting_mesh_color);

  if(cur_color == -1)
    cur_color = I->Obj.Color;

  if(ColorCheckRamped(I->Obj.G, cur_color))
    one_color_flag = false;

  ms->OneColor = cur_color;
  if(ms->V) {
    int ramped_flag = false;
    float *v = ms->V;
    float *vc;
    int *rc;
    int a;
    int state = ms - I->State;
    int n_vert = VLAGetSize(ms->V) / 3;
    int base_n_vert = ms->base_n_V / 3;

    if(ms->VC && (ms->VCsize < n_vert)) {
      FreeP(ms->VC);
      FreeP(ms->RC);
    }

    if(!ms->VC) {
      ms->VCsize = n_vert;
      ms->VC = Alloc(float, n_vert * 3);
    }
    if(!ms->RC) {
      ms->RC = Alloc(int, n_vert);
    }
    rc = ms->RC;
    vc = ms->VC;
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
      FreeP(ms->VC);
      FreeP(ms->RC);
    } else if((!ramped_flag)
              ||
              (!SettingGet_b(I->Obj.G, NULL, I->Obj.Setting, cSetting_ray_color_ramps))) {
      FreeP(ms->RC);
    }
  }
}

static void ObjectMeshUpdate(ObjectMesh * I)
{
  PyMOLGlobals *G = I->Obj.G;
  int a;
  int c;
  ObjectMeshState *ms;
  ObjectMapState *oms = NULL;
  ObjectMap *map = NULL;

  int *n;
  float *v;
  float carve_buffer;
  int avoid_flag = false;
  int *old_n;
  float *old_v;
  int n_cur;
  int n_seg;
  int n_line;
  int flag;
  int last_flag = 0;
  int h, k, l;
  int i, j;
  int ok = true;
  int mesh_skip = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_mesh_skip);

  MapType *voxelmap;            /* this has nothing to do with isosurfaces... */

  for(a = 0; a < I->NState; a++) {
    ms = I->State + a;
    if(ms->Active) {

      map = ExecutiveFindObjectMapByName(I->Obj.G, ms->MapName);
      if(!map) {
        ok = false;
        PRINTFB(I->Obj.G, FB_ObjectMesh, FB_Errors)
          "ObjectMeshUpdate-Error: map '%s' has been deleted.\n", ms->MapName
          ENDFB(I->Obj.G);
        ms->ResurfaceFlag = false;
      }
      if(map) {
        oms = ObjectMapGetState(map, ms->MapState);
        if(!oms)
          ok = false;
      }
      if(oms) {
        if(ms->RefreshFlag || ms->ResurfaceFlag) {
          if(!ms->Field) {
            ms->Crystal = *(oms->Crystal);
          }

          if(I->Obj.RepVis[cRepCell]) {
            if(ms->UnitCellCGO)
              CGOFree(ms->UnitCellCGO);
            ms->UnitCellCGO = CrystalGetUnitCellCGO(&ms->Crystal);
          }

          if(oms->State.Matrix) {
            ObjectStateSetMatrix(&ms->State, oms->State.Matrix);
          } else if(ms->State.Matrix) {
            ObjectStateResetMatrix(&ms->State);
          }
          ms->RefreshFlag = false;
          ms->displayListInvalid = true;
        }
      }

      if(map && oms && ms->N && ms->V && I->Obj.RepVis[cRepMesh]) {
        if(ms->ResurfaceFlag) {
          Isofield *field = NULL;
          ms->RecolorFlag = true;
          ms->ResurfaceFlag = false;
          if(!ms->quiet) {
            PRINTFB(G, FB_ObjectMesh, FB_Details)
              " ObjectMesh: updating \"%s\".\n", I->Obj.Name ENDFB(G);
          }
          if(ms->Field) {
            field = ms->Field;
          } else if(oms->Field) {
            field = oms->Field;
          }

          if(field) {
            {
              float *min_ext, *max_ext;
              float tmp_min[3], tmp_max[3];
              if(MatrixInvTransformExtentsR44d3f(ms->State.Matrix,
                                                 ms->ExtentMin, ms->ExtentMax,
                                                 tmp_min, tmp_max)) {
                min_ext = tmp_min;
                max_ext = tmp_max;
              } else {
                min_ext = ms->ExtentMin;
                max_ext = ms->ExtentMax;
              }

              IsosurfGetRange(I->Obj.G, field, oms->Crystal,
                              min_ext, max_ext, ms->Range, true);
            }
            /*                      printf("Mesh-DEBUG: %d %d %d %d %d %d\n",
               ms->Range[0],
               ms->Range[1],
               ms->Range[2],
               ms->Range[3],
               ms->Range[4],
               ms->Range[5]); */
            IsosurfVolume(I->Obj.G, I->Obj.Setting, NULL,
                          field,
                          ms->Level,
                          &ms->N, &ms->V,
                          ms->Range, ms->MeshMode, mesh_skip, ms->AltLevel);

            if(!SettingGet_b
               (I->Obj.G, I->Obj.Setting, NULL, cSetting_mesh_negative_visible)) {
              ms->base_n_V = VLAGetSize(ms->V);
            } else if(ms->MeshMode != 3) {
              /* do we want the negative surface too? */

              int *N2 = VLAlloc(int, 10000);
              float *V2 = VLAlloc(float, 10000);

              IsosurfVolume(I->Obj.G, I->Obj.Setting, NULL,
                            field,
                            -ms->Level,
                            &N2, &V2, ms->Range, ms->MeshMode, mesh_skip, ms->AltLevel);

              if(N2 && V2) {

                int base_n_N = VLAGetSize(ms->N);
                int base_n_V = VLAGetSize(ms->V);
                int addl_n_N = VLAGetSize(N2);
                int addl_n_V = VLAGetSize(V2);

                ms->base_n_V = base_n_V;

                /* make room */

                VLASize(ms->N, int, base_n_N + addl_n_N);
                VLASize(ms->V, float, base_n_V + addl_n_V);

                /* copy vertex data */

                memcpy(((char *) ms->V) + (sizeof(float) * base_n_V),
                       V2, sizeof(float) * addl_n_V);

                /* copy strip counts */

                memcpy(((char *) ms->N) + (sizeof(int) * (base_n_N - 1)),
                       N2, sizeof(int) * addl_n_N);
                ms->N[base_n_N + addl_n_N - 1] = 0;

                VLAFreeP(N2);
                VLAFreeP(V2);
              }

            }

            if(ms->State.Matrix && VLAGetSize(ms->N) && VLAGetSize(ms->V)) {
              int count;
              /* take map coordinates back to view coordinates if necessary */
              v = ms->V;
              count = VLAGetSize(ms->V) / 3;
              while(count--) {
                transform44d3f(ms->State.Matrix, v, v);
                v += 3;
              }
            }

          }
          if(ms->CarveFlag && ms->AtomVertex && VLAGetSize(ms->N) && VLAGetSize(ms->V)) {
            carve_buffer = ms->CarveBuffer;
            if(ms->CarveBuffer < 0.0F) {
              avoid_flag = true;
              carve_buffer = -carve_buffer;
            }

            /* cull my friend, cull */
            voxelmap = MapNew(I->Obj.G,
                              -carve_buffer, ms->AtomVertex,
                              VLAGetSize(ms->AtomVertex) / 3, NULL);
            if(voxelmap) {

              MapSetupExpress(voxelmap);

              old_n = ms->N;
              old_v = ms->V;
              ms->N = VLAlloc(int, VLAGetSize(old_n));
              ms->V = VLAlloc(float, VLAGetSize(old_v));

              n = old_n;
              v = old_v;
              n_cur = 0;
              n_seg = 0;
              n_line = 0;
              while(*n) {
                last_flag = false;
                c = *(n++);
                while(c--) {
                  flag = false;
                  MapLocus(voxelmap, v, &h, &k, &l);
                  i = *(MapEStart(voxelmap, h, k, l));
                  if(i) {
                    j = voxelmap->EList[i++];
                    while(j >= 0) {
                      if(within3f(ms->AtomVertex + 3 * j, v, carve_buffer)) {
                        flag = true;
                        break;
                      }
                      j = voxelmap->EList[i++];
                    }
                  }
                  if(avoid_flag)
                    flag = !flag;
                  if(flag && (!last_flag)) {
                    VLACheck(ms->V, float, 3 * (n_line + 1));
                    copy3f(v, ms->V + n_line * 3);
                    n_cur++;
                    n_line++;
                  }
                  if(flag && last_flag) {       /* continue segment */
                    VLACheck(ms->V, float, 3 * (n_line + 1));
                    copy3f(v, ms->V + n_line * 3);
                    n_cur++;
                    n_line++;
                  }
                  if((!flag) && last_flag) {    /* terminate segment */
                    VLACheck(ms->N, int, n_seg);
                    ms->N[n_seg] = n_cur;
                    n_seg++;
                    n_cur = 0;
                  }
                  last_flag = flag;
                  v += 3;
                }
                if(last_flag) { /* terminate segment */
                  VLACheck(ms->N, int, n_seg);
                  ms->N[n_seg] = n_cur;
                  n_seg++;
                  n_cur = 0;
                }
              }
              VLACheck(ms->N, int, n_seg);
              ms->N[n_seg] = 0;
              VLAFreeP(old_n);
              VLAFreeP(old_v);
              MapFree(voxelmap);
            }
          }
        }
        if(ms->RecolorFlag) {
          ObjectMeshStateUpdateColors(I, ms);
          ms->RecolorFlag = false;
        }
      }
    }
    SceneInvalidate(I->Obj.G);
  }
  if(!I->Obj.ExtentFlag) {
    ObjectMeshRecomputeExtent(I);
    if(I->Obj.ExtentFlag)
      SceneInvalidate(I->Obj.G);
  }
}

static void ObjectMeshRender(ObjectMesh * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->Obj.G;
  float *v = NULL;
  float *vc;
  int *rc;
  float radius;
  int state = info->state;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int pass = info->pass;
  int *n = NULL;
  int c;
  int a = 0;
  float line_width = SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_mesh_width);
  ObjectMeshState *ms = NULL;
  line_width = SceneGetDynamicLineWidth(info, line_width);

  ObjectPrepareContext(&I->Obj, ray);

  if(state >= 0)
    if(state < I->NState)
      if(I->State[state].Active)
        if(I->State[state].V && I->State[state].N)
          ms = I->State + state;

  while(1) {
    if(state < 0) {             /* all_states */
      ms = I->State + a;
    } else {
      if(!ms) {
        if(I->NState &&
           ((SettingGet(I->Obj.G, cSetting_static_singletons) && (I->NState == 1))))
          ms = I->State;
      }
    }
    if(ms) {
      if(ms->Active && ms->V && ms->N) {
        v = ms->V;
        n = ms->N;
        if(ray) {

          if(ms->UnitCellCGO && (I->Obj.RepVis[cRepCell]))
            CGORenderRay(ms->UnitCellCGO, ray, ColorGet(I->Obj.G, I->Obj.Color),
                         I->Obj.Setting, NULL);
          if(ms->MeshMode != 1) {
            radius = SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_mesh_radius);

            if(radius == 0.0F) {
              radius = ray->PixelRadius * line_width / 2.0F;
            }
          } else {
            radius = SettingGet_f(I->Obj.G, I->Obj.Setting, NULL, cSetting_dot_radius);
            if(radius == 0.0F) {
              radius =
                ray->PixelRadius * SettingGet_f(I->Obj.G, I->Obj.Setting, NULL,
                                                cSetting_dot_width) / 1.4142F;
            }
          }

          if(n && v && I->Obj.RepVis[cRepMesh]) {
            float cc[3];
            float colA[3], colB[3];
            ColorGetEncoded(G, ms->OneColor, cc);
            vc = ms->VC;
            rc = ms->RC;
            /* vc = ColorGet(I->Obj.G,I->Obj.Color); */
            if(ms->MeshMode == 1) {
              ray->fColor3fv(ray, cc);
              while(*n) {
                c = *(n++);
                while(c--) {
                  if(vc) {
                    register float *cA = vc;
                    if(rc) {
                      if(rc[0] < -1)
                        ColorGetEncoded(G, rc[0], (cA = colA));
                      rc++;
                    }
                    ray->fColor3fv(ray, cA);
                    ray->fSphere3fv(ray, v, radius);
                    vc += 3;
                  } else {
                    ray->fSphere3fv(ray, v, radius);
                  }
                  v += 3;
                }
              }
            } else {
              while(*n) {
                c = *(n++);
                if(c--) {
                  v += 3;
                  if(vc) {
                    vc += 3;
                    if(rc)
                      rc++;
                  }
                  while(c--) {
                    if(vc) {
                      register float *cA = vc - 3, *cB = vc;
                      if(rc) {
                        if(rc[-1] < -1)
                          ColorGetEncoded(G, rc[-1], (cA = colA));
                        if(rc[0] < -1)
                          ColorGetEncoded(G, rc[0], (cB = colB));
                        rc++;
                      }
                      ray->fSausage3fv(ray, v - 3, v, radius, cA, cB);
                      vc += 3;
                    } else {
                      ray->fSausage3fv(ray, v - 3, v, radius, cc, cc);
                    }
                    v += 3;
                  }
                }
              }
            }
          }
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
          } else {
            if(!pass) {
              int use_dlst;
              if(ms->UnitCellCGO && (I->Obj.RepVis[cRepCell]))
                CGORenderGL(ms->UnitCellCGO, ColorGet(I->Obj.G, I->Obj.Color),
                            I->Obj.Setting, NULL, info);

              if(!info->line_lighting)
                glDisable(GL_LIGHTING);
              SceneResetNormal(I->Obj.G, false);
              ObjectUseColor(&I->Obj);
              use_dlst = (int) SettingGet(I->Obj.G, cSetting_use_display_lists);

              if(use_dlst && ms->displayList && ms->displayListInvalid) {
                glDeleteLists(ms->displayList, 1);
                ms->displayList = 0;
                ms->displayListInvalid = false;
              }

              if(use_dlst && ms->displayList) {
                glCallList(ms->displayList);
              } else {

                if(use_dlst) {
                  if(!ms->displayList) {
                    ms->displayList = glGenLists(1);
                    if(ms->displayList) {
                      glNewList(ms->displayList, GL_COMPILE_AND_EXECUTE);
                    }
                  }
                }

                if(n && v && I->Obj.RepVis[cRepMesh]) {
                  vc = ms->VC;
                  if(ms->MeshMode == 1)
                    glPointSize(SettingGet_f
                                (I->Obj.G, I->Obj.Setting, NULL, cSetting_dot_width));
                  else
                    glLineWidth(line_width);
                  while(*n) {
                    c = *(n++);
                    if(ms->MeshMode == 1)
                      glBegin(GL_POINTS);
                    else
                      glBegin(GL_LINE_STRIP);
                    while(c--) {
                      if(vc) {
                        glColor3fv(vc);
                        vc += 3;
                      }
                      glVertex3fv(v);
                      v += 3;
                    }
                    glEnd();
                  }
                }
                if(use_dlst && ms->displayList) {
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


/*========================================================================*/

static int ObjectMeshGetNStates(ObjectMesh * I)
{
  return (I->NState);
}


/*========================================================================*/
ObjectMesh *ObjectMeshNew(PyMOLGlobals * G)
{
  OOAlloc(G, ObjectMesh);

  ObjectInit(G, (CObject *) I);

  I->NState = 0;
  I->State = VLAMalloc(10, sizeof(ObjectMeshState), 5, true);   /* autozero important */

  I->Obj.type = cObjectMesh;

  I->Obj.fFree = (void (*)(CObject *)) ObjectMeshFree;
  I->Obj.fUpdate = (void (*)(CObject *)) ObjectMeshUpdate;
  I->Obj.fRender = (void (*)(CObject *, RenderInfo *)) ObjectMeshRender;
  I->Obj.fInvalidate = (void (*)(CObject *, int, int, int)) ObjectMeshInvalidate;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectMeshGetNStates;
  /*  I->Obj.fGetCaption = (char *(*)(CObject *))ObjectMeshGetCaption; */
  return (I);
}


/*========================================================================*/
void ObjectMeshStateInit(PyMOLGlobals * G, ObjectMeshState * ms)
{
  if(ms->Active)
    ObjectStatePurge(&ms->State);
  if(ms->Field) {
    IsosurfFieldFree(ms->State.G, ms->Field);
    ms->Field = NULL;
  }
  ObjectStateInit(G, &ms->State);
  if(!ms->V) {
    ms->V = VLAlloc(float, 10000);
  }
  if(!ms->N) {
    ms->N = VLAlloc(int, 10000);
  }
  if(ms->AtomVertex) {
    VLAFreeP(ms->AtomVertex);
  }
  ms->N[0] = 0;
  ms->Active = true;
  ms->ResurfaceFlag = true;
  ms->RecolorFlag = false;
  ms->ExtentFlag = false;
  ms->CarveFlag = false;
  ms->quiet = true;
  ms->CarveBuffer = 0.0;
  ms->AtomVertex = NULL;
  ms->UnitCellCGO = NULL;
  ms->displayList = 0;
  ms->displayListInvalid = true;
  ms->caption[0] = 0;
  ms->Field = NULL;
}


/*========================================================================*/
ObjectMesh *ObjectMeshFromXtalSym(PyMOLGlobals * G, ObjectMesh * obj, ObjectMap * map,
                                  CSymmetry * sym,
                                  int map_state,
                                  int state, float *mn, float *mx,
                                  float level, int meshMode,
                                  float carve, float *vert_vla,
                                  float alt_level, int quiet)
{
  int ok = true;
  ObjectMesh *I;
  ObjectMeshState *ms;
  ObjectMapState *oms;

  if(!obj) {
    I = ObjectMeshNew(G);
  } else {
    I = obj;
  }

  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectMeshState, state);
    I->NState = state + 1;
  }

  ms = I->State + state;
  ObjectMeshStateInit(G, ms);

  strcpy(ms->MapName, map->Obj.Name);
  ms->MapState = map_state;
  oms = ObjectMapGetState(map, map_state);

  ms->Level = level;
  ms->AltLevel = alt_level;
  ms->MeshMode = meshMode;
  ms->quiet = quiet;
  if(oms) {

    if((meshMode == 3) && (ms->AltLevel < ms->Level)) {
      /* gradient object -- need to auto-set range */
      if(!ObjectMapStateGetDataRange(G, oms, &ms->Level, &ms->AltLevel)) {
        ms->Level = -1.0F;
        ms->AltLevel = 1.0F;
      }
    }

    copy3f(mn, ms->ExtentMin);  /* this is not exactly correct...should actually take vertex points from range */
    copy3f(mx, ms->ExtentMax);

    if(oms->State.Matrix) {
      ObjectStateSetMatrix(&ms->State, oms->State.Matrix);
    } else if(ms->State.Matrix) {
      ObjectStateResetMatrix(&ms->State);
    }

    {
      float *min_ext, *max_ext;
      float tmp_min[3], tmp_max[3];
      if(MatrixInvTransformExtentsR44d3f(ms->State.Matrix,
                                         ms->ExtentMin, ms->ExtentMax,
                                         tmp_min, tmp_max)) {
        min_ext = tmp_min;
        max_ext = tmp_max;
      } else {
        min_ext = ms->ExtentMin;
        max_ext = ms->ExtentMax;
      }

      {
        int eff_range[6];

        if(IsosurfGetRange
           (G, oms->Field, oms->Crystal, min_ext, max_ext, eff_range, false)) {
          int fdim[3];
          int expand_result;
          /* need to generate symmetry-expanded temporary map */

          ms->Crystal = *(oms->Crystal);
          fdim[0] = eff_range[3] - eff_range[0];
          fdim[1] = eff_range[4] - eff_range[1];
          fdim[2] = eff_range[5] - eff_range[2];
          ms->Field = IsosurfFieldAlloc(I->Obj.G, fdim);

          expand_result =
            IsosurfExpand(oms->Field, ms->Field, oms->Crystal, sym, eff_range);

          if(expand_result == 0) {
            ok = false;
            if(!quiet) {
              PRINTFB(G, FB_ObjectMesh, FB_Warnings)
                " ObjectMesh-Warning: no symmetry expanded map points found.\n" ENDFB(G);
            }
          } else {
            if(!quiet) {
              PRINTFB(G, FB_ObjectMesh, FB_Warnings)
                " ObjectMesh-Warning: not all symmetry expanded points covered by map.\n"
                ENDFB(G);
            }
          }

          ms->Range[0] = 0;
          ms->Range[1] = 0;
          ms->Range[2] = 0;
          ms->Range[3] = fdim[0];
          ms->Range[4] = fdim[1];
          ms->Range[5] = fdim[2];

        } else {
          /* mesh entirely contained within bounds of current map */
          int a;
          for(a = 0; a < 6; a++) {
            ms->Range[a] = eff_range[a];
          }
        }
      }
    }
    ms->ExtentFlag = true;
  }
  if(ok) {
    if(carve != 0.0) {
      ms->CarveFlag = true;
      ms->CarveBuffer = carve;
      ms->AtomVertex = vert_vla;
    }
    if(I) {
      ObjectMeshRecomputeExtent(I);
    }
    I->Obj.ExtentFlag = true;
    /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]); */
  }
  if(!ok) {
    ObjectMeshStateFree(ms);
    I = NULL;
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


/*========================================================================*/
ObjectMesh *ObjectMeshFromBox(PyMOLGlobals * G, ObjectMesh * obj, ObjectMap * map,
                              int map_state,
                              int state, float *mn, float *mx,
                              float level, int meshMode,
                              float carve, float *vert_vla, float alt_level, int quiet)
{
  ObjectMesh *I;
  ObjectMeshState *ms;
  ObjectMapState *oms;

  if(!obj) {
    I = ObjectMeshNew(G);
  } else {
    I = obj;
  }

  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectMeshState, state);
    I->NState = state + 1;
  }

  ms = I->State + state;
  ObjectMeshStateInit(G, ms);

  strcpy(ms->MapName, map->Obj.Name);
  ms->MapState = map_state;
  oms = ObjectMapGetState(map, map_state);

  ms->Level = level;
  ms->AltLevel = alt_level;
  ms->MeshMode = meshMode;
  ms->quiet = quiet;
  if(oms) {

    if((meshMode == 3) && (ms->AltLevel < ms->Level)) {
      /* gradient object -- need to auto-set range */
      if(!ObjectMapStateGetDataRange(G, oms, &ms->Level, &ms->AltLevel)) {
        ms->Level = -1.0F;
        ms->AltLevel = 1.0F;
      }
    }

    copy3f(mn, ms->ExtentMin);  /* this is not exactly correct...should actually take vertex points from range */
    copy3f(mx, ms->ExtentMax);

    if(oms->State.Matrix) {
      ObjectStateSetMatrix(&ms->State, oms->State.Matrix);
    } else if(ms->State.Matrix) {
      ObjectStateResetMatrix(&ms->State);
    }

    {
      float *min_ext, *max_ext;
      float tmp_min[3], tmp_max[3];
      if(MatrixInvTransformExtentsR44d3f(ms->State.Matrix,
                                         ms->ExtentMin, ms->ExtentMax,
                                         tmp_min, tmp_max)) {
        min_ext = tmp_min;
        max_ext = tmp_max;
      } else {
        min_ext = ms->ExtentMin;
        max_ext = ms->ExtentMax;
      }

      IsosurfGetRange(G, oms->Field, oms->Crystal, min_ext, max_ext, ms->Range, true);
    }
    ms->ExtentFlag = true;
  }
  if(carve != 0.0) {
    ms->CarveFlag = true;
    ms->CarveBuffer = carve;
    ms->AtomVertex = vert_vla;
  }
  if(I) {
    ObjectMeshRecomputeExtent(I);
  }
  I->Obj.ExtentFlag = true;
  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]); */
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


/*========================================================================*/

void ObjectMeshRecomputeExtent(ObjectMesh * I)
{
  int extent_flag = false;
  int a;
  ObjectMeshState *ms;

  for(a = 0; a < I->NState; a++) {
    ms = I->State + a;
    if(ms->Active) {
      if(ms->ExtentFlag) {
        if(!extent_flag) {
          extent_flag = true;
          copy3f(ms->ExtentMax, I->Obj.ExtentMax);
          copy3f(ms->ExtentMin, I->Obj.ExtentMin);
        } else {
          max3f(ms->ExtentMax, I->Obj.ExtentMax, I->Obj.ExtentMax);
          min3f(ms->ExtentMin, I->Obj.ExtentMin, I->Obj.ExtentMin);
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
