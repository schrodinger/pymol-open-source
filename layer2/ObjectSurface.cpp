
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
#include"ObjectSurface.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"Map.h"
#include"Parse.h"
#include"Tetsurf.h"
#include"CarveHelper.h"
#include"ContourSurf.h"
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
#include"ShaderMgr.h"
#include"CGO.h"

static void ObjectSurfaceRecomputeExtent(ObjectSurface * I);

static PyObject *ObjectSurfaceStateAsPyList(ObjectSurfaceState * I)
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
  PyList_SetItem(result, 13, PyInt_FromLong(I->DotFlag));
  PyList_SetItem(result, 14, PyInt_FromLong(static_cast<int>(I->Mode)));
  PyList_SetItem(result, 15, PyInt_FromLong(static_cast<int>(I->Side)));
  PyList_SetItem(result, 16, PyInt_FromLong(I->quiet));

  return (PConvAutoNone(result));
}

static PyObject *ObjectSurfaceAllStatesAsPyList(ObjectSurface * I)
{

  auto result = PyList_New(I->State.size());
  for(int a = 0; a < I->State.size(); a++) {
    if(I->State[a].Active) {
      PyList_SetItem(result, a, ObjectSurfaceStateAsPyList(&I->State[a]));
    } else {
      PyList_SetItem(result, a, PConvAutoNone(NULL));
    }
  }
  return (PConvAutoNone(result));

}

static int ObjectSurfaceStateFromPyList(PyMOLGlobals * G, ObjectSurfaceState * I,
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
      *I = ObjectSurfaceState(G);
      if(ok)
        ok = PyList_Check(list);
      if(ok)
        ll = PyList_Size(list);
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
        ok = PConvPyIntToInt(PyList_GetItem(list, 13), &I->DotFlag);

      if(ok)
        PConvFromPyListItem(G, list, 14, I->Mode);

      if(ok && (ll > 15))
        PConvFromPyListItem(G, list, 15, I->Side);

      if(ok && (ll > 16))
        PConvPyIntToInt(PyList_GetItem(list, 16), &I->quiet);

      if(ok) {
        I->RefreshFlag = true;
        I->ResurfaceFlag = true;
      }
    }
  }
  return (ok);
}

static int ObjectSurfaceAllStatesFromPyList(ObjectSurface * I, PyObject * list, int size)
{
  int ok = true;
  I->State.reserve(size);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    for(int a = 0; a < size; a++) {
      CPythonVal *val = CPythonVal_PyList_GetItem(I->G, list, a);
      I->State.emplace_back(I->G);
      ok = ObjectSurfaceStateFromPyList(I->G, &I->State.back(), val);
      CPythonVal_Free(val);
      if(!ok)
        break;
    }
  }
  return (ok);
}

int ObjectSurfaceNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectSurface ** result)
{
  int ok = true;
  ObjectSurface *I = NULL;
  (*result) = NULL;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);

  I = new ObjectSurface(G);
  if(ok)
    ok = (I != NULL);

  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  int size = 0;
  if(ok)
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 1, &size);
  if(ok){
    CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 2);
    ok = ObjectSurfaceAllStatesFromPyList(I, val, size);
    CPythonVal_Free(val);
  }
  if(ok) {
    (*result) = I;
    ObjectSurfaceRecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
}

PyObject *ObjectSurfaceAsPyList(ObjectSurface * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->State.size()));
  PyList_SetItem(result, 2, ObjectSurfaceAllStatesAsPyList(I));

  return (PConvAutoNone(result));
}

void ObjectSurfaceDump(ObjectSurface * I, const char *fname, int state, int quiet)
{
  FILE* f = fopen(fname, "wb");
  if(!f)
    ErrMessage(I->G, "ObjectSurfaceDump", "can't open file for writing");
  else {
    if(state < I->State.size()) {
      auto n = I->State[state].N.data();
      auto v = I->State[state].V.data();
      if(n && v)
        while(*n) {
          v += 12;
          int c = *(n++);
          c -= 4;
          bool backface = true;
          const float *v1, *v2;
          while(c > 0) {
            if ((backface = !backface)) {
              v1 = v - 6;
              v2 = v - 12;
            } else {
              v1 = v - 12;
              v2 = v - 6;
            }
            fprintf(f,
                    "%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
                    *(v1 + 3), *(v1 + 4), *(v1 + 5), *(v1), *(v1 + 1), *(v1 + 2),
                    *(v2 + 3), *(v2 + 4), *(v2 + 5), *(v2), *(v2 + 1), *(v2 + 2),
                    *(v + 3), *(v + 4), *(v + 5), *(v), *(v + 1), *(v + 2));

            v += 6;
            c -= 2;
          }
        }
    }
    fclose(f);
    if (!quiet) {
      PRINTFB(I->G, FB_ObjectSurface, FB_Actions)
        " ObjectSurfaceDump: %s written to %s\n", I->Name, fname ENDFB(I->G);
    }
  }
}

void ObjectSurface::invalidate(cRep_t rep, cRepInv_t level, int state)
{
  auto I = this;
  int once_flag = true;
  if(level >= cRepInvExtents) {
    I->ExtentFlag = false;
  }
  if((rep == cRepSurface) || (rep == cRepMesh) || (rep == cRepAll)) {
    for(int a = 0; a < I->State.size(); a++) {
      if(state < 0)
        once_flag = false;
      if(!once_flag)
        state = a;
      I->State[state].RefreshFlag = true;
      if(level >= cRepInvRep) {
        I->State[state].ResurfaceFlag = true;
	if(I->State[state].shaderCGO){
	  I->State[state].shaderCGO.reset();
	}
        SceneChanged(I->G);
      } else if(level >= cRepInvColor) {
        I->State[state].RecolorFlag = true;
	if(I->State[state].shaderCGO){
	  I->State[state].shaderCGO.reset();
	}
        SceneChanged(I->G);
      } else {
        SceneInvalidate(I->G);
      }
      if(once_flag)
        break;
    }
  }
}

int ObjectSurfaceInvalidateMapName(ObjectSurface * I, const char *name, const char * new_name)
{
  int result = false;
  for(int a = 0; a < I->State.size(); a++) {
    auto ms = &I->State[a];
    if(ms->Active) {
      if(strcmp(ms->MapName, name) == 0) {
        if (new_name)
          strcpy(ms->MapName, new_name);
        I->invalidate(cRepAll, cRepInvAll, a);
        result = true;
      }
    }
  }
  return result;
}

static void ObjectSurfaceStateUpdateColors(ObjectSurface * I, ObjectSurfaceState * ms)
{
  int one_color_flag = true;
  int cur_color =
    SettingGet_color(I->G, I->Setting.get(), NULL, cSetting_surface_color);

  if(cur_color == -1)
    cur_color = I->Color;

  if(ColorCheckRamped(I->G, cur_color))
    one_color_flag = false;

  ms->OneColor = cur_color;
  if(ms->V) {
    int ramped_flag = false;
    float *v = ms->V.data();
    float *vc;
    int *rc;
    int a;
    int state = std::distance(I->State.data(), ms);
    int base_n_vert = ms->base_n_V;
    switch (ms->Mode) {
    case cIsosurfaceMode::triangles_grad_normals:
    case cIsosurfaceMode::triangles_tri_normals:
      {
        int n_vert = VLAGetSize(ms->V) / 6;
        base_n_vert /= 6;

        if(!ms->VC.empty() && (ms->VCsize() < n_vert)) {
          ms->VC.clear();
          ms->RC.clear();
        }

        if(ms->VC.empty()) {
          ms->VC.resize(n_vert * 3);
        }
        if(ms->RC.empty()) {
          ms->RC.resize(n_vert);
        }
        rc = ms->RC.data();
        vc = ms->VC.empty() ? nullptr : ms->VC.data();
        v += 3;
        if(vc) {
          for(a = 0; a < n_vert; a++) {
            if(a == base_n_vert) {
              int new_color = SettingGet_color(I->G, I->Setting.get(),
                                               NULL, cSetting_surface_negative_color);
              if(new_color == -1)
                new_color = cur_color;
              if(new_color != cur_color) {
                one_color_flag = false;
                cur_color = new_color;
              }
            }
            if(ColorCheckRamped(I->G, cur_color)) {
              ColorGetRamped(I->G, cur_color, v, vc, state);
              *rc = cur_color;
              ramped_flag = true;
            } else {
              const float *col = ColorGet(I->G, cur_color);
              copy3f(col, vc);
            }
            rc++;
            vc += 3;
            v += 6;             /* alternates with normals */
          }
        }
      }
      break;
    default:
      {
        int n_vert = VLAGetSize(ms->V) / 3;
        base_n_vert /= 3;
        if(!ms->VC.empty() && (ms->VCsize() < n_vert)) {
          ms->VC.clear();
          ms->RC.clear();
        }

        if(ms->VC.empty()) {
          ms->VC.resize(n_vert * 3);
        }
        if(ms->RC.empty()) {
          ms->RC.resize(n_vert);
        }
        rc = ms->RC.data();
        vc = ms->VC.empty() ? nullptr : ms->VC.data();
        if(vc) {
          for(a = 0; a < n_vert; a++) {
            if(a == base_n_vert) {
              int new_color = SettingGet_color(I->G, I->Setting.get(),
                                               NULL, cSetting_surface_negative_color);
              if(new_color == -1)
                new_color = cur_color;
              if(new_color != cur_color)
                one_color_flag = false;
              cur_color = new_color;
            }

            if(ColorCheckRamped(I->G, cur_color)) {
              ColorGetRamped(I->G, cur_color, v, vc, state);
              *rc = cur_color;
              ramped_flag = true;
            } else {
              const float *col = ColorGet(I->G, cur_color);
              copy3f(col, vc);
            }
            rc++;
            vc += 3;
            v += 3;
          }
        }
      }
      break;
    }

    if(one_color_flag && (!ramped_flag)) {
      ms->VC.clear();
      ms->RC.clear();
    } else if((!ramped_flag)
              ||
              (!SettingGet_b(I->G, NULL, I->Setting.get(), cSetting_ray_color_ramps))) {
      ms->RC.clear();
    }
  }
}

void ObjectSurface::update()
{
  auto I = this;
  for(auto& msref : I->State) {
    ObjectSurfaceState *ms = &msref;
    ObjectMapState *oms = NULL;
    ObjectMap *map = NULL;

    if(ms->Active) {
      map = ExecutiveFindObjectMapByName(I->G, ms->MapName);
      if(!map) {
        PRINTFB(I->G, FB_ObjectSurface, FB_Errors)
          "ObjectSurfaceUpdate-Error: map '%s' has been deleted.\n", ms->MapName
          ENDFB(I->G);
        ms->ResurfaceFlag = false;
      }
      if(map) {
        oms = ObjectMapGetState(map, ms->MapState);
      }
      if(oms) {
        if(!oms->Matrix.empty()) {
          ObjectStateSetMatrix(ms, oms->Matrix.data());
        } else if(!ms->Matrix.empty()) {
          ObjectStateResetMatrix(ms);
        }

        if(I->visRep & cRepCellBit){
          if (!ms->UnitCellCGO || ms->RefreshFlag || ms->ResurfaceFlag) {
            ms->Crystal = oms->Symmetry->Crystal;
            if((I->visRep & cRepCellBit)) {
              ms->UnitCellCGO.reset(CrystalGetUnitCellCGO(&ms->Crystal));
            }
          }
          ms->RefreshFlag = false;
        }
      }
      if(map && ms && oms && ms->N && ms->V && (I->visRep & cRepSurfaceBit)) {
        if(ms->ResurfaceFlag) {
          ms->ResurfaceFlag = false;
          ms->RecolorFlag = true;
          if(!ms->quiet) {
            PRINTFB(I->G, FB_ObjectSurface, FB_Details)
              " ObjectSurface: updating \"%s\".\n", I->Name ENDFB(I->G);
          }

          ms->shaderCGO.reset();

          if(oms->Field) {

            {
              float *min_ext, *max_ext;
              float tmp_min[3], tmp_max[3];
              if(MatrixInvTransformExtentsR44d3f(ms->Matrix.data(),
                                                 ms->ExtentMin, ms->ExtentMax,
                                                 tmp_min, tmp_max)) {
                min_ext = tmp_min;
                max_ext = tmp_max;
              } else {
                min_ext = ms->ExtentMin;
                max_ext = ms->ExtentMax;
              }

              TetsurfGetRange(I->G, oms->Field.get(), &oms->Symmetry->Crystal,
                              min_ext, max_ext, ms->Range);
            }

            std::unique_ptr<CarveHelper> carvehelper;
            if(ms->CarveFlag && ms->AtomVertex) {
              carvehelper.reset(new CarveHelper(G, ms->CarveBuffer,
                  ms->AtomVertex, ms->AtomVertex.size() / 3));
            }

            ms->nT = ContourSurfVolume(I->G, oms->Field.get(),
                                   ms->Level,
                                   ms->N, ms->V,
                                   ms->Range,
                                   ms->Mode,
                                   carvehelper.get(),
                                   ms->Side);

            if(!SettingGet_b
               (I->G, I->Setting.get(), NULL, cSetting_surface_negative_visible)) {
              ms->base_n_V = VLAGetSize(ms->V);
            } else {
              /* do we want the negative surface too? */

              int nT2;
              pymol::vla<int> N2(10000);
              pymol::vla<float> V2(10000);

              nT2 = ContourSurfVolume(I->G, oms->Field.get(),
                                  -ms->Level,
                                  N2, V2,
                                  ms->Range,
                                  ms->Mode,
                                  carvehelper.get(),
                                  ms->Side);
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
                std::copy_n(V2.data(), addl_n_V, ms->V.data() + base_n_V);

                /* copy strip counts */
                std::copy_n(N2.data(), addl_n_N, ms->N.data() + base_n_N - 1);

                ms->N[base_n_N + addl_n_N - 1] = 0;

                ms->nT += nT2;
              }
            }

            if(!ms->Matrix.empty()) {      /* in we're in a different reference frame... */
              double *matrix = ms->Matrix.data();
              auto v = ms->V.data();
              auto n = ms->N.data();

              if(n && v) {

                while(*n) {
                  int c = *(n++);
                  switch (ms->Mode) {
                  case cIsosurfaceMode::triangles_grad_normals:
                  case cIsosurfaceMode::triangles_tri_normals:
                    transform44d3fas33d3f(matrix, v, v);
                    transform44d3f(matrix, v + 3, v + 3);
                    transform44d3fas33d3f(matrix, v + 6, v + 6);
                    transform44d3f(matrix, v + 9, v + 9);
                    v += 12;
                    c -= 4;
                    while(c > 0) {
                      transform44d3fas33d3f(matrix, v, v);
                      transform44d3f(matrix, v + 3, v + 3);
                      v += 6;
                      c -= 2;
                    }
                    break;
                  case cIsosurfaceMode::lines:
                    transform44d3f(matrix, v, v);
                    c--;
                    v += 3;
                    while(c > 0) {
                      transform44d3f(matrix, v, v);
                      v += 3;
                      c--;
                    }
                    break;
                  case cIsosurfaceMode::dots:
                  default:
                    while(c > 0) {
                      transform44d3f(matrix, v, v);
                      v += 3;
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
          ObjectSurfaceStateUpdateColors(I, ms);
          ms->RecolorFlag = false;
        }
      }
    }
  }
  if(!I->ExtentFlag) {
    ObjectSurfaceRecomputeExtent(I);
  }
  SceneInvalidate(I->G);
}

static void ObjectSurfaceRenderGlobalTransparency(PyMOLGlobals * G,
    RenderInfo * info, ObjectSurfaceState *ms, const float *col, float alpha)
{
  auto v = ms->V.data();
  const float* vc = ms->VC.empty() ? nullptr : ms->VC.data();
  auto n = ms->N.data();
  
  while(*n) {
    int parity = 1;
    int c = *(n++);

    v += 6;
    if(vc)
      vc += 3;
    c -= 2;

    v += 6;
    if(vc)
      vc += 3;
    c -= 2;

    while(c > 0) {
      if(vc) {
        CGOAlphaTriangle(info->alpha_cgo,
                         v + (3 - 6), v + (3 - 12), v + 3,
                         v - 6, v - 12, v,
                         vc - 3, vc - 6, vc,
                         alpha, alpha, alpha, parity);
      } else {
        CGOAlphaTriangle(info->alpha_cgo,
                         v + (3 - 6), v + (3 - 12), v + 3,
                         v - 6, v - 12, v,
                         col, col, col,
                         alpha, alpha, alpha, parity);
      }
      v += 6;
      if(vc)
        vc += 3;
      c -= 2;
      parity = !parity;
    }
  }
}

static void ObjectSurfaceRenderUnOptimizedTransparency(ObjectSurfaceState *ms, float alpha){
  auto v = ms->V.data();
  const float* vc = ms->VC.empty() ? nullptr : ms->VC.data();
  auto n = ms->N.data();

  while(*n) {
    int c = *(n++);
    CGOBegin(ms->shaderCGO.get(), GL_TRIANGLE_STRIP);
    while(c > 0) {
      CGONormalv(ms->shaderCGO.get(), v);
      v += 3;
      if(vc) {
        CGOColorv(ms->shaderCGO.get(), vc);
        vc += 3;
      }
      CGOVertexv(ms->shaderCGO.get(), v);
      v += 3;
      c -= 2;
    }
    CGOEnd(ms->shaderCGO.get());
  }
}

static void ObjectSurfaceRenderOpaque(PyMOLGlobals * G, ObjectSurface * I, ObjectSurfaceState *ms){
  auto v = ms->V.data();
  auto n = ms->N.data();
  CGOSpecial(ms->shaderCGO.get(), LINEWIDTH_DYNAMIC_MESH);
  const float* vc = ms->VC.empty() ? nullptr : ms->VC.data();

  while(*n) {
    int c = *(n++);
    switch (ms->Mode) {
    case cIsosurfaceMode::triangles_grad_normals:
    case cIsosurfaceMode::triangles_tri_normals:
      CGOBegin(ms->shaderCGO.get(), GL_TRIANGLE_STRIP);
      while(c > 0) {
        CGONormalv(ms->shaderCGO.get(), v);
        v += 3;
        if(vc) {
          CGOColorv(ms->shaderCGO.get(), vc);
          vc += 3;
        }
        CGOVertexv(ms->shaderCGO.get(), v);
        v += 3;
        c -= 2;
      }
      CGOEnd(ms->shaderCGO.get());
      break;
    case cIsosurfaceMode::lines:
      CGOBegin(ms->shaderCGO.get(), GL_LINES);
      while(c > 0) {
        if(vc) {
          CGOColorv(ms->shaderCGO.get(), vc);
          vc += 3;
        }
        CGOVertexv(ms->shaderCGO.get(), v);
        v += 3;
        c--;
      }
      CGOEnd(ms->shaderCGO.get());
      break;
    case cIsosurfaceMode::dots:
    default:
      CGOBegin(ms->shaderCGO.get(), GL_POINTS);
      while(c > 0) {
        if(vc) {
          CGOColorv(ms->shaderCGO.get(), vc);
          vc += 3;
        }
        CGOVertexv(ms->shaderCGO.get(), v);
        v += 3;
        c--;
      }
      CGOEnd(ms->shaderCGO.get());
    }
  }
}

static void ObjectSurfaceRenderRay(PyMOLGlobals * G, ObjectSurface *I,
    RenderInfo * info, ObjectSurfaceState *ms)
{
  float *v = ms->V.data();
  int c;
  int* n = ms->N.data();
  float alpha = 1.0F - SettingGet_f(G, NULL, I->Setting.get(), cSetting_transparency);
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;


  CRay *ray = info->ray;
  if(ms->UnitCellCGO && (I->visRep & cRepCellBit)){
    int rayok = CGORenderRay(ms->UnitCellCGO.get(), ray, info, ColorGet(G, I->Color),
                             NULL, I->Setting.get(), NULL);
    if (!rayok){
      ms->UnitCellCGO.reset();
    }
  }
  
  ray->transparentf(1.0F - alpha);
  ms->Radius = SettingGet_f(G, I->Setting.get(), NULL, cSetting_mesh_radius);
  if(ms->Radius == 0.0F) {
    ms->Radius = ray->PixelRadius *
      SettingGet_f(I->G, I->Setting.get(), NULL, cSetting_mesh_width) / 2.0F;
  }
  
  if(n && v && (I->visRep & cRepSurfaceBit)) {
    float cc[3];
    float colA[3], colB[3], colC[3];
    ColorGetEncoded(G, ms->OneColor, cc);
    float* vc = ms->VC.empty() ? nullptr : ms->VC.data();
    
    const int* rc = ms->RC.empty() ? nullptr : ms->RC.data();
    while(*n) {
      c = *(n++);
      switch (ms->Mode) {
      case cIsosurfaceMode::triangles_grad_normals:
      case cIsosurfaceMode::triangles_tri_normals:
        v += 12;
        if(vc)
          vc += 6;
        c -= 4;
        while(c > 0) {
          if(vc) {
            float *cA = vc - 6, *cB = vc - 3, *cC = vc;
            if(rc) {
              if(rc[0] < -1)
                ColorGetEncoded(G, rc[0], (cA = colA));
              if(rc[1] < -1)
                ColorGetEncoded(G, rc[1], (cB = colB));
              if(rc[2] < -1)
                ColorGetEncoded(G, rc[2], (cC = colC));
              rc++;
            }
            ray->triangle3fv(v - 9, v - 3, v + 3,
                             v - 12, v - 6, v, cA, cB, cC);
            vc += 3;
          } else {
            ray->triangle3fv(v - 9, v - 3, v + 3,
                             v - 12, v - 6, v, cc, cc, cc);
          }
          v += 6;
          c -= 2;
        }
        break;
      case cIsosurfaceMode::lines:
        c--;
        v += 3;
        if(vc)
          vc += 3;
        while(c > 0) {
          if(vc) {
            float *cA = vc - 3, *cB = vc;
            if(rc) {
              if(rc[0] < -1)
                ColorGetEncoded(G, rc[0], (cA = colA));
              if(rc[1] < -1)
                ColorGetEncoded(G, rc[1], (cB = colB));
              rc++;
            }
            ray->sausage3fv(v - 3, v, ms->Radius, cA, cB);
            vc += 3;
          } else
            ray->sausage3fv(v - 3, v, ms->Radius, cc, cc);
          v += 3;
          c--;
        }
        break;
      case cIsosurfaceMode::dots:
      default:
        while(c > 0) {
          if(vc) {
            ray->color3fv(vc);
            vc += 3;
          }
          ray->sphere3fv(v, ms->Radius);
          v += 3;
          c--;
        }
        break;
      }
    }
  }
  ray->transparentf(0.0);
}

static void ObjectSurfaceRenderCell(PyMOLGlobals *G, ObjectSurface * I,
    RenderInfo * info, ObjectSurfaceState *ms, short use_shader)
{
  const float *color = ColorGet(G, I->Color);
  if (use_shader != ms->UnitCellCGO->has_draw_buffers){
    if (use_shader){
      CGO *convertcgo = CGOOptimizeToVBONotIndexed(ms->UnitCellCGO.get(), 0);
      ms->UnitCellCGO.reset(convertcgo);
      assert(ms->UnitCellCGO->use_shader);
    } else {
      ms->UnitCellCGO.reset(CrystalGetUnitCellCGO(&ms->Crystal));
    }
  }
  CGORenderGL(ms->UnitCellCGO.get(), color,
              I->Setting.get(), NULL, info, NULL);
}

void ObjectSurface::render(RenderInfo * info)
{
  auto I = this;
  int state = info->state;
  CRay *ray = info->ray;
  auto pick = info->pick;
  const RenderPass pass = info->pass;
  const float *col;
  ObjectSurfaceState *ms = NULL;
  float alpha;
  ObjectPrepareContext(I, info);

  alpha = 1.0F - SettingGet_f(G, NULL, I->Setting.get(), cSetting_transparency);
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;

  StateIterator iter(G, I->Setting.get(), state, I->State.size());
  while(iter.next()) {
    ms = &I->State[iter.state];
    if(ms && ms->Active && ms->V && ms->N) {
      if(ray) {
        ObjectSurfaceRenderRay(G, I, info, ms);
      } else if(G->HaveGUI && G->ValidContext) {
        if(!pick) {  // no picking for ObjectSurfaces
          int render_now = false;
          short use_shader;
          use_shader = SettingGetGlobal_b(G, cSetting_surface_use_shader) & 
            SettingGetGlobal_b(G, cSetting_use_shaders);

          if(info && info->alpha_cgo) {
            render_now = (pass == RenderPass::Opaque);
            use_shader = false;
          } else if(alpha < 1.0F) {
            render_now = (pass == RenderPass::Transparent);
          } else {
            render_now = (pass == RenderPass::Opaque);
          }

          if((I->visRep & cRepCellBit) && ms->UnitCellCGO && (pass == RenderPass::Opaque)){
            ObjectSurfaceRenderCell(G, I, info, ms, use_shader);
          }

          if(render_now) {
            if (ms->shaderCGO && use_shader != ms->shaderCGO->has_draw_buffers){
              ms->shaderCGO.reset();
            }

            if (ms->shaderCGO){
              CGORenderGL(ms->shaderCGO.get(), NULL, NULL, NULL, info, NULL);
              continue;
            }
            
            // Generating CGO
            ms->shaderCGO.reset(CGONew(G));
            ms->shaderCGO->use_shader = true;

            CGOResetNormal(ms->shaderCGO.get(), false);

            col = ColorGet(G, ms->OneColor);
            if((alpha != 1.0)) {
              CGOAlpha(ms->shaderCGO.get(), alpha);
            }
            CGOColorv(ms->shaderCGO.get(), col);
            
            if(I->visRep & cRepSurfaceBit) {
              if (static_cast<int>(ms->Mode) > 1 &&
                  alpha != 1.0) {         /* transparent triangles */
                if(info->alpha_cgo) {     /* global transparency */
                  ObjectSurfaceRenderGlobalTransparency(G, info, ms, col, alpha);
                } else {  /* cgo transparency with sorting if needed */
                  ObjectSurfaceRenderUnOptimizedTransparency(ms, alpha);
                }
              } else { /* opaque, triangles */
                ObjectSurfaceRenderOpaque(G, I, ms);
              }
            }
            CGOStop(ms->shaderCGO.get());

            if (use_shader){
              auto convertcgo = CGOOptimizeToVBOIndexed(ms->shaderCGO.get(), 0,
                  nullptr, true, (alpha != 1.0) /* embedTransparency */);
              if (convertcgo){
                ms->shaderCGO.reset(convertcgo);
              }
              ms->shaderCGO->use_shader = true;              
              CGORenderGL(ms->shaderCGO.get(), NULL, NULL, NULL, info, NULL);
            } else {
              if (alpha != 1.0){
                // use_shader = 0
                CGO *convertcgo = CGOConvertTrianglesToAlpha(ms->shaderCGO.get());
                ms->shaderCGO.reset(convertcgo);
                ms->shaderCGO->render_alpha = 1;
              }
              ms->shaderCGO->use_shader = false;
              CGORenderGL(ms->shaderCGO.get(), NULL, NULL, NULL, info, NULL);
            }
          }
        }
      }
    }
  }
}


/*========================================================================*/

int ObjectSurface::getNFrame() const
{
  return State.size();
}


/*========================================================================*/
ObjectSurface::ObjectSurface(PyMOLGlobals* G)
    : pymol::CObject(G)
{
  type = cObjectSurface;
}

/*========================================================================*/
ObjectSurfaceState::ObjectSurfaceState(PyMOLGlobals* G)
    : CObjectState(G)
    , Crystal(G)
    , Active(true)
{
  V = pymol::vla<float>(10000);
  N = pymol::vla<int>(10000);
}

/*========================================================================*/
ObjectSurface *ObjectSurfaceFromBox(PyMOLGlobals * G, ObjectSurface * obj,
                                    ObjectMap * map, int map_state, int state, float *mn,
                                    float *mx, float level, cIsosurfaceMode mode, float carve,
                                    pymol::vla<float>&& vert_vla, cIsosurfaceSide side, int quiet)
{
  ObjectSurface *I;
  ObjectSurfaceState *ms;
  ObjectMapState *oms;

  if(!obj) {
    I = new ObjectSurface(G);
  } else {
    I = obj;
  }

  if(state < 0)
    state = I->State.size();
  if(I->State.size() <= state) {
    VecCheckEmplace(I->State, state, G);
  }

  ms = &I->State[state];
  *ms = ObjectSurfaceState(G);

  strcpy(ms->MapName, map->Name);
  ms->MapState = map_state;
  oms = ObjectMapGetState(map, map_state);

  ms->Level = level;
  ms->Mode = mode;
  ms->Side = side;
  ms->quiet = quiet;
  if(oms) {

    if(!oms->Matrix.empty()) {
      ObjectStateSetMatrix(ms, oms->Matrix.data());
    } else if(!ms->Matrix.empty()) {
      ObjectStateResetMatrix(ms);
    }

    copy3f(mn, ms->ExtentMin);  /* this is not exactly correct...should actually take vertex points from range */
    copy3f(mx, ms->ExtentMax);

    {
      float *min_ext, *max_ext;
      float tmp_min[3], tmp_max[3];
      if(MatrixInvTransformExtentsR44d3f(ms->Matrix.data(),
                                         ms->ExtentMin, ms->ExtentMax,
                                         tmp_min, tmp_max)) {
        min_ext = tmp_min;
        max_ext = tmp_max;
      } else {
        min_ext = ms->ExtentMin;
        max_ext = ms->ExtentMax;
      }

      TetsurfGetRange(G, oms->Field.get(), &oms->Symmetry->Crystal, min_ext, max_ext, ms->Range);
    }
    ms->ExtentFlag = true;
  }
  if(carve != 0.0) {
    ms->CarveFlag = true;
    ms->CarveBuffer = carve;
    ms->AtomVertex = std::move(vert_vla);

    const double *matrix = ObjectStateGetInvMatrix(ms);

    if(matrix) {
      int n = VLAGetSize(ms->AtomVertex) / 3;
      float *v = ms->AtomVertex.data();
      while(n--) {
        /* convert back into original map coordinates 
           for surface carving operation */
        transform44d3f(matrix, v, v);
        v += 3;
      }
    }

  }
  if(I) {
    ObjectSurfaceRecomputeExtent(I);
  }
  I->ExtentFlag = true;
  /*  printf("Brick %d %d %d %d %d %d\n",I->Range[0],I->Range[1],I->Range[2],I->Range[3],I->Range[4],I->Range[5]); */
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}

pymol::Result<float> ObjectSurfaceGetLevel(ObjectSurface * I, int state)
{
  if(state >= int(I->State.size())) {
    return pymol::make_error("Invalid surface state");
  } else {
    if(state < 0) {
      state = 0;
    }
    auto ms = &I->State[state];
    if(ms->Active) {
      return ms->Level;
    } else {
      return pymol::make_error("Invalid Surface state");
    }
  }
}

int ObjectSurfaceSetLevel(ObjectSurface * I, float level, int state, int quiet)
{
  int ok = true;
  int once_flag = true;
  if(state >= int(I->State.size())) {
    ok = false;
  } else {
    for(int a = 0; a < I->State.size(); a++) {
      if(state < 0) {
        once_flag = false;
      }
      if(!once_flag) {
        state = a;
      }
      auto ms = &I->State[state];
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


/*========================================================================*/

void ObjectSurfaceRecomputeExtent(ObjectSurface * I)
{
  int extent_flag = false;

  for(auto& msr : I->State) {
    auto ms = &msr;
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

  if(I->TTTFlag && I->ExtentFlag) {
    const float *ttt;
    double tttd[16];
    if(ObjectGetTTT(I, &ttt, -1)) {
      convertTTTfR44d(ttt, tttd);
      MatrixTransformExtentsR44d3f(tttd,
                                   I->ExtentMin, I->ExtentMax,
                                   I->ExtentMin, I->ExtentMax);
    }
  }
}

/*========================================================================*/

pymol::CObject* ObjectSurface::clone() const
{
  return new ObjectSurface(*this);
}

