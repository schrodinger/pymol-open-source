
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
#include"ObjectCGO.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"VFont.h"
#include"ShaderMgr.h"

static PyObject *ObjectCGOStateAsPyList(ObjectCGOState * I)
{
  PyObject *result = NULL;

  result = PyList_New(1);
  if(I->origCGO)
    PyList_SetItem(result, 0, CGOAsPyList(I->origCGO));
  else
    PyList_SetItem(result, 0, PConvAutoNone(NULL));
  return (PConvAutoNone(result));

}

static PyObject *ObjectCGOAllStatesAsPyList(ObjectCGO * I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NState);
  for(a = 0; a < I->NState; a++) {
    PyList_SetItem(result, a, ObjectCGOStateAsPyList(I->State + a));
  }
  return (PConvAutoNone(result));
}

static int ObjectCGOStateFromPyList(PyMOLGlobals * G, ObjectCGOState * I, PyObject * list,
                                    int version)
{
  int ok = true;
  int ll, pl = 0;
  PyObject *tmp;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  CGOFree(I->origCGO);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok && ll==2) {
    tmp = PyList_GetItem(list, 0);
    if(tmp == Py_None)
      I->origCGO = NULL;
    else {
      ok = ((I->origCGO = CGONewFromPyList(G, tmp, version, 1)) != NULL);
    }
    pl++;
  }
  if(ok && !I->origCGO) {
    tmp = PyList_GetItem(list, pl);
    if(tmp == Py_None)
      I->origCGO = NULL;
    else {
      ok = ((I->origCGO = CGONewFromPyList(G, tmp, version, 0)) != NULL);
    }
  }
  return (ok);
}

static int ObjectCGOAllStatesFromPyList(ObjectCGO * I, PyObject * list, int version)
{
  int ok = true;
  int a;
  VLACheck(I->State, ObjectCGOState, I->NState);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    for(a = 0; a < I->NState; a++) {
      auto *val = PyList_GetItem(list, a);
      ok =
        ObjectCGOStateFromPyList(I->G, I->State + a, val,
                                 version);
      if(!ok)
        break;
    }
  }
  return (ok);

}

int ObjectCGONewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectCGO ** result,
                           int version)
{
  int ok = true;
  ObjectCGO *I = NULL;
  (*result) = NULL;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);

  I = new ObjectCGO(G);
  if(ok)
    ok = (I != NULL);

  if(ok){
    auto *val = PyList_GetItem(list, 0);
    ok = ObjectFromPyList(G, val, I);
  }
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->NState);
  if(ok)
    ok = ObjectCGOAllStatesFromPyList(I, PyList_GetItem(list, 2), version);
  if(ok) {
    (*result) = I;
    ObjectCGORecomputeExtent(I);
  } else {
    /* cleanup? */
  }
  return (ok);
}

PyObject *ObjectCGOAsPyList(ObjectCGO * I)
{
  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(I));
  PyList_SetItem(result, 1, PyInt_FromLong(I->NState));
  PyList_SetItem(result, 2, ObjectCGOAllStatesAsPyList(I));

  return (PConvAutoNone(result));
}


/*========================================================================*/

ObjectCGO::~ObjectCGO()
{
  auto I = this;
  for(int a = 0; a < I->NState; a++) {
    CGOFree(I->State[a].renderCGO);
    CGOFree(I->State[a].origCGO);
  }
  VLAFreeP(I->State);
}


/*========================================================================*/

void ObjectCGORecomputeExtent(ObjectCGO * I)
{
  float mx[3], mn[3];
  int extent_flag = false;
  int a;
  int has_normals = 0;
  CGO *cgo;
  for(a = 0; a < I->NState; a++){
    cgo = I->State[a].origCGO;
    if (!cgo){
      cgo = I->State[a].renderCGO;
    }
    if(cgo) {
      if(CGOGetExtent(cgo, mn, mx)) {
        if(!extent_flag) {
          extent_flag = true;
          copy3f(mx, I->ExtentMax);
          copy3f(mn, I->ExtentMin);
        } else {
          max3f(mx, I->ExtentMax, I->ExtentMax);
          min3f(mn, I->ExtentMin, I->ExtentMin);
        }
      }
      if (!has_normals && cgo && CGOHasNormals(cgo)){
	has_normals = 1;
      }
    }
  }
  I->ExtentFlag = extent_flag;
  SettingCheckHandle(I->G, &I->Setting);
  SettingSet_i(I->Setting, cSetting_cgo_lighting, has_normals);
}


/*========================================================================*/
void ObjectCGO::invalidate(int rep, int level, int state)
{
  auto I = this;
  ObjectCGOState *sobj = NULL;

  if(state < 0) {
    int a;
    for(a = 0; a < I->NState; a++) {
      sobj = I->State + a;
      CGOFree(sobj->renderCGO);
    }
  } else {
    if((state >= 0) && (state < I->NState)) {
      sobj = I->State + state;
      CGOFree(sobj->renderCGO);
    }
  }
}


/*========================================================================*/

void ObjectCGO::update()
{
  for (int a = 0; a < NState; ++a) {
    CGOFree(State[a].renderCGO);
  }
  SceneInvalidate(G);    /* needed ? */
}


/*========================================================================*/

int ObjectCGO::getNFrame() const
{
  return NState;
}

static void ObjectCGORenderState(PyMOLGlobals* G, int pass, CRay* ray,
    ObjectCGO* I, RenderInfo* info, ObjectCGOState* sobj, const float* color,
    ObjectGadgetRamp* ramp, int use_shader, bool cgo_lighting)
{
  if(ray) {
    if(sobj) {
      if(sobj->origCGO){
        CGO *cgo = sobj->origCGO, *cgo_copy = NULL;
        if (cgo_lighting && CGOHasAnyTriangleVerticesWithoutNormals(cgo)) {
          cgo = cgo_copy = CGOGenerateNormalsForTriangles(cgo);
        }
        CGORenderRay(cgo, ray, info, color, ramp, I->Setting, NULL);
        CGOFree(cgo_copy);
      }
    }
  } else if(G->HaveGUI && G->ValidContext && pass) {
    if(info->pick) { // no picking yet
    } else {
      bool pass_is_opaque = (pass > 0);
      if(sobj && ((sobj->hasTransparency ^ pass_is_opaque) || (sobj->hasOpaque == pass_is_opaque))){
	{
	  CShaderPrg *shaderPrg;
	  int two_sided_lighting = SettingGet_i(G, I->Setting, NULL, cSetting_two_sided_lighting);
          bool backface_cull = SettingGet_i(G, I->Setting, NULL, cSetting_backface_cull);
	  if (two_sided_lighting<0){
	    two_sided_lighting = !cgo_lighting;
	  }
	  two_sided_lighting &= cgo_lighting;  // only set two_sided_lighting if cgo_lighting is set
#ifndef PURE_OPENGL_ES_2
	  if (cgo_lighting){
	    glEnable(GL_LIGHTING);
	  } else {
	    glDisable(GL_LIGHTING);
	  }
	  if (two_sided_lighting){
	    if (use_shader)
	      glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
	    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	  } else {
	    if (use_shader)
	      glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);
	    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	  }
#endif
	  if (backface_cull){
	    glCullFace(GL_BACK);
	    glEnable(GL_CULL_FACE);
	  }

	  if (use_shader){
	    shaderPrg = G->ShaderMgr->Enable_DefaultShader(info->pass);
	    if (!shaderPrg) return;
	    shaderPrg->SetLightingEnabled(cgo_lighting);
	    shaderPrg->Set1i("two_sided_lighting_enabled", two_sided_lighting);
	    sobj->renderCGO->use_shader = use_shader;
	    sobj->renderCGO->debug = SettingGetGlobal_i(G, cSetting_cgo_debug);
	    CGORenderGL(sobj->renderCGO, color, I->Setting, NULL, info, NULL);
	    shaderPrg->Disable();
	  } else {
	    sobj->renderCGO->use_shader = use_shader;
	    sobj->renderCGO->debug = SettingGetGlobal_i(G, cSetting_cgo_debug);
	    CGORenderGL(sobj->renderCGO, color, I->Setting, NULL, info, NULL);
	  }

	    if (backface_cull){
	      glDisable(GL_CULL_FACE);
	    }
#ifndef PURE_OPENGL_ES_2
	    if (two_sided_lighting){
	      if (use_shader)
		glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);
	      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	    }
	    if (!cgo_lighting){
	      glEnable(GL_LIGHTING);
	    }
#endif
	}
      }
    }
  }
}

static void ObjectCGOGenerateCGO(PyMOLGlobals * G, ObjectCGO * I, ObjectCGOState *sobj, bool use_shader, bool cgo_lighting, const float *color, ObjectGadgetRamp *ramp, int state)
{
  if (sobj->renderCGO && 
      ((use_shader ^ sobj->renderWithShaders) ||
       (cgo_lighting ^ sobj->cgo_lighting))){
    // if renderWithShaders doesn't match use_shader, clear CGO and re-generate
    CGOFree(sobj->renderCGO);
  }
  if (!sobj->renderCGO){
    float colorWithA[4];
    short someLinesWithoutNormals = 0;
    if (color){
      colorWithA[0] = color[0]; colorWithA[1] = color[1]; colorWithA[2] = color[2];
    } else {
      colorWithA[0] = 1.f; colorWithA[1] = 1.f; colorWithA[2] = 1.f;
    }
    colorWithA[3] = 1.f - SettingGet_f(G, I->Setting, NULL, cSetting_cgo_transparency);

    bool hasTransparency = (colorWithA[3] < 1.f || CGOHasTransparency(sobj->origCGO));
    bool hasOpaque = (colorWithA[3] == 1.f || CGOHasOpaque(sobj->origCGO));

    CGO *allCylinders = NULL;
    CGO *allSpheres = NULL;
    std::unique_ptr<CGO, CGODeleter> preOpt;

    {
      std::unique_ptr<CGO, CGODeleter> inputWithLighting;
      const CGO* inputCGO = sobj->origCGO;

      if (cgo_lighting){
	if (CGOHasAnyTriangleVerticesWithoutNormals(inputCGO)){
	  // we only need normals if cgo_lighting is on
          inputWithLighting.reset(CGOGenerateNormalsForTriangles(inputCGO));
          inputCGO = inputWithLighting.get();
	}
	someLinesWithoutNormals = CGOHasAnyLineVerticesWithoutNormals(inputCGO);
	if (!use_shader && someLinesWithoutNormals){
	  // if some lines without normals, turn lighting off on lines
          inputWithLighting.reset(CGOTurnLightingOnLinesOff(inputCGO, use_shader));
          inputCGO = inputWithLighting.get();
	}
      }

      CGO* convertcgo = CGONew(G);
      CGOColorv(convertcgo, colorWithA);
      CGOAlpha(convertcgo, colorWithA[3]);
      CGOAppend(convertcgo, inputCGO);
      inputWithLighting.reset();

      if (use_shader){
        bool t_mode_3 = SettingGetGlobal_i(G, cSetting_transparency_mode)==3;
        if ((t_mode_3 || !hasTransparency)
            && G->ShaderMgr->Get_DefaultSphereShader(0)
            && G->ShaderMgr->Get_CylinderShader(0))
        {
          if (CGOHasCylinderOperations(convertcgo)){
            allCylinders = CGONew(G);
            CGOEnable(allCylinders, GL_CYLINDER_SHADER);
            CGO* newCGO = CGOConvertShaderCylindersToCylinderShader(convertcgo, allCylinders);
            allCylinders->free_append(newCGO);
            assert(newCGO == nullptr);
            CGODisable(allCylinders, GL_CYLINDER_SHADER);
            CGOStop(allCylinders);
            CGO *allButCylinders = CGONew(G);
            CGOFilterOutCylinderOperationsInto(convertcgo, allButCylinders);
            CGOStop(allButCylinders);
            CGOFree(convertcgo);
            convertcgo = allButCylinders;
          }
          if (CGOHasOperationsOfType(convertcgo, CGO_SPHERE)){
            CGO *allButSpheres = CGONew(G);
            allSpheres = CGOOptimizeSpheresToVBONonIndexed(convertcgo, 0, true, allButSpheres);
            if (allSpheres){
              CGOFree(convertcgo);
              CGOStop(allButSpheres);
              convertcgo = allButSpheres;
            } else {
              CGOFree(allButSpheres);
            }
          }
          preOpt.reset(CGOSimplify(convertcgo, 0));
        } else {
          preOpt.reset(CGOSimplifyNoCompress(convertcgo, 0));
        }
      } else {
        preOpt.reset(CGOSimplifyNoCompress(convertcgo, 0));
      }
      CGOFree(convertcgo);
    }

    if (ramp){
      preOpt.reset(CGOColorByRamp(G, preOpt.get(), ramp, state, I->Setting));
    }

    sobj->hasTransparency = hasTransparency;
    sobj->hasOpaque = hasOpaque;

    if (use_shader){
      if(preOpt && preOpt->has_begin_end){
	preOpt.reset(CGOCombineBeginEnd(preOpt.get(), 0));
      }

      preOpt.reset(CGOOptimizeToVBOIndexedWithColorEmbedTransparentInfo(
          preOpt.get(), 0, colorWithA, false));

      if (someLinesWithoutNormals){
        // if some lines without normals, turn lighting off on lines
        CGO* convertcgo = preOpt.release();
        preOpt.reset(CGOTurnLightingOnLinesOff(convertcgo, use_shader));
        CGOFreeWithoutVBOs(convertcgo);
      }

      if (allCylinders){
        preOpt->free_append(allCylinders);
      }

      if (allSpheres){
        preOpt->free_append(allSpheres);
      }

      sobj->renderCGO = preOpt.release();
    } else {
      assert(sobj->hasTransparency == CGOHasTransparency(preOpt.get()));
      assert(sobj->hasOpaque == CGOHasOpaque(preOpt.get()));

      if (sobj->hasTransparency) {
        sobj->renderCGO = CGOConvertTrianglesToAlpha(preOpt.get());
	sobj->renderCGO->render_alpha = 2;
      } else {
        sobj->renderCGO = CGOSimplify(preOpt.get(), 0);
      }
    }

    assert(allCylinders == nullptr);
    assert(allSpheres == nullptr);

    sobj->renderWithShaders = use_shader;
    sobj->cgo_lighting = cgo_lighting;
  }
}
/*========================================================================*/

void ObjectCGO::render(RenderInfo * info)
{
  auto I = this;
  int state = info->state;
  CRay *ray = info->ray;
  int pass = info->pass;
  ObjectCGOState *sobj = NULL;
  const float *color = NULL;
  bool use_shader = false, cgo_lighting = false;
  ObjectGadgetRamp *ramp = NULL;
  
  use_shader = SettingGetGlobal_b(G, cSetting_cgo_use_shader) &
    SettingGetGlobal_b(G, cSetting_use_shaders);
  cgo_lighting = SettingGet_i(G, I->Setting, NULL, cSetting_cgo_lighting);

  ObjectPrepareContext(I, info);
  ramp = ColorGetRamp(G, I->Color);
  color = ColorGet(G, I->Color);

  if(!I->State)
    return;

  if(pass || info->ray) {
    if((I->visRep & cRepCGOBit)) {
      for(StateIterator iter(G, I->Setting, state, I->NState); iter.next();) {
        sobj = I->State + iter.state;
        if (!sobj->origCGO)
          continue;
        if (!ray)
          ObjectCGOGenerateCGO(G, I, sobj, use_shader, cgo_lighting, color, ramp, iter.state);
          ObjectCGORenderState(G, pass, ray, I, info, sobj, color, ramp, use_shader, cgo_lighting);
      }
    }
  }
}


/*========================================================================*/
ObjectCGO::ObjectCGO(PyMOLGlobals * G) : CObject(G)
{
  State = VLACalloc(ObjectCGOState, 10);
  type = cObjectCGO;
}


/*========================================================================*/
static CGO *ObjectCGOPyListFloatToCGO(PyMOLGlobals * G, PyObject * list)
{
  CGO *cgo = NULL;
  int len;
  int ok = true;
  int result;
  float *raw = NULL;
  if(PyList_Check(list)) {
    len = PConvPyListToFloatArray(list, &raw);
    if(len < 0)
      len = 0;
    if(raw) {
      if(ok) {
        cgo = CGONewSized(G, len);
        if(cgo) {
          result = CGOFromFloatArray(cgo, raw, len);
          if(result) {
            PRINTF " FloatToCGO: error encountered on element %d\n", result ENDF(G);
          }
          CGOStop(cgo);
        }
      }
      FreeP(raw);
    }
  }
  return (cgo);
}


/*========================================================================*/
static CGO *ObjectCGOFloatArrayToCGO(PyMOLGlobals * G, float *raw, int len, int quiet)
{
  CGO *cgo = NULL;
  int ok = true;
  int result;

  if(raw) {
    if(ok) {
      cgo = CGONewSized(G, len);
      if(cgo) {
        result = CGOFromFloatArray(cgo, raw, len);
        if(result && !quiet) {
          PRINTF " FloatToCGO: error encountered on element %d\n", result ENDF(G);
        }
        CGOStop(cgo);
      }
    }
  }
  return (cgo);
}


/*========================================================================*/
ObjectCGO *ObjectCGOFromCGO(PyMOLGlobals * G, ObjectCGO * obj, CGO * cgo, int state)
{
  ObjectCGO *I = NULL;

  if(obj) {
    if(obj->type != cObjectCGO)     /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = new ObjectCGO(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCGOState, state);
    I->NState = state + 1;
  }

  CGOFree(I->State[state].renderCGO);
  CGOFree(I->State[state].origCGO);
  I->State[state].origCGO = cgo;

  if(I) {
    ObjectCGORecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


/*========================================================================*/

ObjectCGO *ObjectCGONewVFontTest(PyMOLGlobals * G, const char *text, float *pos)
{

  ObjectCGO *I = NULL;
  int font_id;
  CGO *cgo = NULL;
  float scale[2] = { 1.0, 1.0 };

  font_id = VFontLoad(G, 1, 1, 1, true);
  cgo = CGONew(G);
  VFontWriteToCGO(G, font_id, cgo, text, pos, scale, NULL, NULL);
  I = ObjectCGOFromCGO(G, NULL, cgo, 0);

  return (I);
}


/*========================================================================*/
ObjectCGO *ObjectCGODefine(PyMOLGlobals * G, ObjectCGO * obj, PyObject * pycgo, int state)
{                               /* assumes blocked interpreter */
  ObjectCGO *I = NULL;

  CGO *cgo, *font_cgo;
  int est;

  if(obj) {
    if(obj->type != cObjectCGO)     /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = new ObjectCGO(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCGOState, state);
    I->NState = state + 1;
  }

  CGOFree(I->State[state].origCGO);

  if(PyList_Check(pycgo)) {
    if(PyList_Size(pycgo)) {
      if(PyFloat_Check(PyList_GetItem(pycgo, 0))) {
        cgo = ObjectCGOPyListFloatToCGO(G, pycgo);
        if(cgo) {
          est = CGOCheckForText(cgo);
          if(est) {
            CGOPreloadFonts(cgo);
            font_cgo = CGODrawText(cgo, est, NULL);
            CGOFree(cgo);
            cgo = font_cgo;
          }
          est = CGOCheckComplex(cgo);
	  I->State[state].origCGO = cgo;
        } else {
          ErrMessage(G, "ObjectCGO", "could not parse CGO List.");
        }
      }
    }
  }
  if(I) {
    ObjectCGORecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}

ObjectCGO *ObjectCGOFromFloatArray(PyMOLGlobals * G, ObjectCGO * obj,
                                   float *array, int size, int state, int quiet)
{
  ObjectCGO *I = NULL;

  CGO *cgo, *font_cgo;
  int est;

  if(obj) {
    if(obj->type != cObjectCGO)     /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = new ObjectCGO(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCGOState, state);
    I->NState = state + 1;
  }

  CGOFree(I->State[state].renderCGO);
  CGOFree(I->State[state].origCGO);

  cgo = ObjectCGOFloatArrayToCGO(G, array, size, quiet);
  if(cgo) {
    est = CGOCheckForText(cgo);
    if(est) {
      CGOPreloadFonts(cgo);
      font_cgo = CGODrawText(cgo, est, NULL);
      CGOFree(cgo);
      cgo = font_cgo;
    }
    est = CGOCheckComplex(cgo);
    I->State[state].origCGO = cgo;
  } else if(!quiet) {
    ErrMessage(G, "ObjectCGO", "could not parse CGO.");
  }
  if(I) {
    ObjectCGORecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


