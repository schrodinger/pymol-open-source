
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

#ifndef _PYMOL_NOPY
static PyObject *ObjectCGOStateAsPyList(ObjectCGOState * I)
{

  PyObject *result = NULL;

  result = PyList_New(1);
  /*  NO LONGER NEED TO SAVE I->std, just ray tracing version 
  if(I->std)
    PyList_SetItem(result, 0, CGOAsPyList(I->std));
  else
  PyList_SetItem(result, 0, PConvAutoNone(NULL));*/
  if(I->ray)
    PyList_SetItem(result, 0, CGOAsPyList(I->ray));
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
  int ll;
  PyObject *tmp;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok) {
    tmp = PyList_GetItem(list, 0);
    if(tmp == Py_None)
      I->std = NULL;
    else
      ok = ((I->std = CGONewFromPyList(G, PyList_GetItem(list, 0), version)) != NULL);
  }
  if(ok) {
    tmp = PyList_GetItem(list, 1);
    if(tmp == Py_None)
      I->ray = NULL;
    else
      ok = ((I->ray = CGONewFromPyList(G, PyList_GetItem(list, 1), version)) != NULL);
    if (!I->std && I->ray){
      I->std = CGOSimplify(I->ray, 0);
#ifdef _PYMOL_CGO_DRAWARRAYS
      {
	CGO *convertcgo = NULL;
	if(I->std && I->std->has_begin_end){
	  convertcgo = CGOCombineBeginEnd(I->std, 0);
	  CGOFree(I->std);
	  I->std = convertcgo;
	}
      }
#endif
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
      ok =
        ObjectCGOStateFromPyList(I->Obj.G, I->State + a, PyList_GetItem(list, a),
                                 version);
      if(!ok)
        break;
    }
  }
  return (ok);

}
#endif

int ObjectCGONewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectCGO ** result,
                           int version)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  ObjectCGO *I = NULL;
  (*result) = NULL;
  if(ok)
    ok = (list != Py_None);
  if(ok)
    ok = PyList_Check(list);

  I = ObjectCGONew(G);
  if(ok)
    ok = (I != NULL);

  if(ok)
    ok = ObjectFromPyList(G, PyList_GetItem(list, 0), &I->Obj);
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
#endif
}

PyObject *ObjectCGOAsPyList(ObjectCGO * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  result = PyList_New(3);
  PyList_SetItem(result, 0, ObjectAsPyList(&I->Obj));
  PyList_SetItem(result, 1, PyInt_FromLong(I->NState));
  PyList_SetItem(result, 2, ObjectCGOAllStatesAsPyList(I));

  return (PConvAutoNone(result));
#endif
}


/*========================================================================*/

void ObjectCGOFree(ObjectCGO * I)
{
  int a;
  for(a = 0; a < I->NState; a++) {
    if(I->State[a].shaderCGO && I->State[a].std!=I->State[a].shaderCGO)
      CGOFree(I->State[a].shaderCGO);
    if(I->State[a].std)
      CGOFree(I->State[a].std);
    if(I->State[a].ray)
      CGOFree(I->State[a].ray);
  }
  VLAFreeP(I->State);
  ObjectPurge(&I->Obj);
  OOFreeP(I);
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
    cgo = I->State[a].std;
    if (!cgo){
      cgo = I->State[a].shaderCGO;
    }
    if(cgo) {
      if(CGOGetExtent(cgo, mn, mx)) {
        if(!extent_flag) {
          extent_flag = true;
          copy3f(mx, I->Obj.ExtentMax);
          copy3f(mn, I->Obj.ExtentMin);
        } else {
          max3f(mx, I->Obj.ExtentMax, I->Obj.ExtentMax);
          min3f(mn, I->Obj.ExtentMin, I->Obj.ExtentMin);
        }
      }
      if (!has_normals && cgo && CGOHasNormals(cgo)){
	has_normals = 1;
      }
    }
  }
  I->Obj.ExtentFlag = extent_flag;
  SettingCheckHandle(I->Obj.G, &I->Obj.Setting);
  SettingSet_i(I->Obj.Setting, cSetting_cgo_lighting, has_normals);
}


/*========================================================================*/
static void ObjectCGOInvalidate(ObjectCGO * I, int rep, int level, int state)
{
  ObjectCGOState *sobj = NULL;
  if(state < 0) {
    int a;
    for(a = 0; a < I->NState; a++) {
      I->State[a].valid = false;
      sobj = I->State + a;
      if (sobj->shaderCGO){
	CGOFree(sobj->shaderCGO);	      
	sobj->shaderCGO = 0;
      }
    }
  } else {
    if((state >= 0) && (state < I->NState)) {
      I->State[state].valid = false;
      sobj = I->State + state;
      if (sobj->shaderCGO){
	CGOFree(sobj->shaderCGO);	      
	sobj->shaderCGO = 0;
      }
    }
  }
}


/*========================================================================*/

static void ObjectCGOUpdate(ObjectCGO * I)
{
  int a;
  for(a = 0; a < I->NState; a++) {
    ObjectCGOState *ocs = I->State + a;
    if (ocs->shaderCGO){
      CGOFree(ocs->shaderCGO);	      
      ocs->shaderCGO = 0;
    }
    if(!ocs->valid) {
      if(ocs->std && ocs->ray) {
        int est = CGOCheckComplex(ocs->ray);
        if(est) {
          if(ocs->std)
            CGOFree(ocs->std);
          ocs->std = CGOSimplify(ocs->ray, est);
#ifdef _PYMOL_CGO_DRAWARRAYS
	  {
	    CGO *convertcgo = NULL;
	    if(ocs->std && ocs->std->has_begin_end){
	      convertcgo = CGOCombineBeginEnd(ocs->std, 0);
	      CGOFree(ocs->std);
	      ocs->std = convertcgo;
	    }
	  }
#endif
        }
      }
      ocs->valid = true;
    }
  }
  SceneInvalidate(I->Obj.G);    /* needed ? */
}


/*========================================================================*/

static int ObjectCGOGetNState(ObjectCGO * I)
{
  return (I->NState);
}


/*========================================================================*/

static void ObjectCGORender(ObjectCGO * I, RenderInfo * info)
{
  register PyMOLGlobals *G = I->Obj.G;
  int state = info->state;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  int pass = info->pass;
  ObjectCGOState *sobj = NULL;
  int a;
  float *color;
  int use_shader = 0;
  
  use_shader = (int) SettingGet(G, cSetting_cgo_use_shader) &
    (int)SettingGet(G, cSetting_use_shaders) & 
    (int)!SettingGet(G, cSetting_transparency_global_sort);

  ObjectPrepareContext(&I->Obj, ray);

  color = ColorGet(G, I->Obj.Color);

  if((pass == 1) || info->ray) {
    if(I->Obj.RepVis[cRepCGO]) {

      if(state < I->NState) {
        sobj = I->State + state;
      }
      if(state < 0) {
        if(I->State) {
          for(a = 0; a < I->NState; a++) {
            sobj = I->State + a;
#ifdef _PYMOL_CGO_DRAWBUFFERS
	    if (use_shader){
	      CGO *convertcgo = NULL;
	      if (!sobj->shaderCGO){
		float colorWithA[4];
		if (color){
		  colorWithA[0] = color[0]; colorWithA[1] = color[1]; colorWithA[2] = color[2];
		} else {
		  colorWithA[0] = 1.f; colorWithA[1] = 1.f; colorWithA[2] = 1.f;
		}
		colorWithA[3] = 1.f - SettingGet_f(G, I->Obj.Setting, NULL, cSetting_cgo_transparency);
		if (CGOHasCylinderOperations(sobj->std)){
		  convertcgo = CGOOptimizeGLSLCylindersToVBOIndexedNoColor(sobj->std, 0);
		  //		  convertcgo->enable_shaders = true;
		} else {
		  convertcgo = CGOOptimizeToVBOIndexedWithColor(sobj->std, 0, colorWithA);
		}
		sobj->shaderCGO = convertcgo;
	      }
	    } else if (sobj->shaderCGO){
	      CGOFree(sobj->shaderCGO);	      
	      sobj->shaderCGO = 0;
	    }
#endif
            if(ray) {
	      int try_std = false;
              if(sobj->ray){
                int rayok = CGORenderRay(sobj->ray, ray, color, I->Obj.Setting, NULL);
		if (!rayok){
		  CGOFree(sobj->ray);
		  sobj->ray = NULL;
		  try_std = true;
		}
	      } else {
		try_std = true;
	      }
	      if (try_std && sobj->std){
		int rayok = CGORenderRay(sobj->std, ray, color, I->Obj.Setting, NULL);
		if (!rayok){
		  CGOFree(sobj->std);
		  sobj->std = NULL;
		}
	      }
            } else if(G->HaveGUI && G->ValidContext) {
              if(pick) {
              } else {
		CShaderPrg *shaderPrg;
		int cgo_lighting, two_sided_lighting;
		cgo_lighting = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_cgo_lighting);
		two_sided_lighting = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_two_sided_lighting);
		if (two_sided_lighting<0){
		  two_sided_lighting = SceneGetTwoSidedLighting(G);
		}
		if (use_shader && sobj->shaderCGO){
		  shaderPrg = CShaderPrg_Enable_DefaultShader(G);
		  CShaderPrg_SetLightingEnabled(shaderPrg, cgo_lighting);
		  CShaderPrg_Set1i(shaderPrg, "two_sided_lighting_enabled", two_sided_lighting);
		  sobj->shaderCGO->use_shader = use_shader;
		  sobj->shaderCGO->debug = SettingGet(G, cSetting_cgo_debug);
		  CGORenderGL(sobj->shaderCGO, color, I->Obj.Setting, NULL, info, NULL);
		  CShaderPrg_Disable(shaderPrg);
		} else {
		  sobj->std->use_shader = use_shader;
		  sobj->std->debug = SettingGet(G, cSetting_cgo_debug);
		  if (cgo_lighting){
		    glEnable(GL_LIGHTING);
		  } else {
		    glDisable(GL_LIGHTING);
		  }
		  if (two_sided_lighting){
		    GLLIGHTMODELI(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
		  } else {
		    GLLIGHTMODELI(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
		  }
		  sobj->std->use_shader = use_shader;
		  sobj->std->debug = SettingGet(G, cSetting_cgo_debug);
		  CGORenderGL(sobj->std, color, I->Obj.Setting, NULL, info, NULL);
		  if (SceneGetTwoSidedLighting(G)){
		    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
		  } else {
		    glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);
		  }
          if (!cgo_lighting){
		    glEnable(GL_LIGHTING);
		  }
		}
              }
            }
          }
        }
      } else {
        if(!sobj) {
          if(I->NState && SettingGet(G, cSetting_static_singletons))
            sobj = I->State;
        }
#ifdef _PYMOL_CGO_DRAWBUFFERS
	if (use_shader){
	  if (!sobj->shaderCGO && sobj->std){
	    float colorWithA[4];
	    if (color){
	      colorWithA[0] = color[0]; colorWithA[1] = color[1]; colorWithA[2] = color[2];
	    } else {
	      colorWithA[0] = 1.f; colorWithA[1] = 1.f; colorWithA[2] = 1.f;
	    }
	    colorWithA[3] = 1.f - SettingGet_f(G, I->Obj.Setting, NULL, cSetting_cgo_transparency);
	    if (CGOHasCylinderOperations(sobj->std)){
	      sobj->shaderCGO = CGOOptimizeGLSLCylindersToVBOIndexedNoColor(sobj->std, 0);
	      //	      sobj->shaderCGO->enable_shaders = true;
	    } else {
	      sobj->shaderCGO = CGOOptimizeToVBOIndexedWithColor(sobj->std, 0, colorWithA);
	    }
	  }
	} else if (sobj->shaderCGO){
	  CGOFree(sobj->shaderCGO);	      
	  sobj->shaderCGO = 0;
	}
#endif
        if(ray) {
          if(sobj) {
	    int try_std = false;
            if(sobj->ray){
              int rayok = CGORenderRay(sobj->ray, ray, color, I->Obj.Setting, NULL);
	      if (!rayok){
		CGOFree(sobj->ray);
		sobj->ray = NULL;
		try_std = true;
	      }
	    } else {
	      try_std = true;
	    }
	    if (try_std && sobj->std){
              int rayok = CGORenderRay(sobj->std, ray, color, I->Obj.Setting, NULL);
	      if (!rayok){
		CGOFree(sobj->std);
		sobj->std = NULL;
	      }
	    }
          }
        } else if(G->HaveGUI && G->ValidContext) {
          if(pick) {
          } else {
            if(sobj)
              if(sobj->std){
		CShaderPrg *shaderPrg;
		int cgo_lighting, two_sided_lighting;
		cgo_lighting = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_cgo_lighting);
		two_sided_lighting = SettingGet_i(G, I->Obj.Setting, NULL, cSetting_two_sided_lighting);
		if (two_sided_lighting<0){
		  two_sided_lighting = SceneGetTwoSidedLighting(G);
		}
		if (use_shader){
		  shaderPrg = CShaderPrg_Enable_DefaultShader(G);
		  if (!shaderPrg) return;
		  CShaderPrg_SetLightingEnabled(shaderPrg, cgo_lighting);
		  CShaderPrg_Set1i(shaderPrg, "two_sided_lighting_enabled", two_sided_lighting);
		  sobj->shaderCGO->use_shader = use_shader;
		  sobj->shaderCGO->debug = SettingGet(G, cSetting_cgo_debug);
		  CGORenderGL(sobj->shaderCGO, color, I->Obj.Setting, NULL, info, NULL);
		  CShaderPrg_Disable(shaderPrg);
		} else {
		  sobj->std->use_shader = use_shader;
		  sobj->std->debug = SettingGet(G, cSetting_cgo_debug);
		  if (cgo_lighting){
		    glEnable(GL_LIGHTING);
		  } else {
		    glDisable(GL_LIGHTING);
		  }
		  if (two_sided_lighting){
		    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
		  } else {
		    glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);
		  }
		  CGORenderGL(sobj->std, color, I->Obj.Setting, NULL, info, NULL);
		  if (SceneGetTwoSidedLighting(G)){
		    glEnable(GL_VERTEX_PROGRAM_TWO_SIDE);
		  } else {
		    glDisable(GL_VERTEX_PROGRAM_TWO_SIDE);
		  }
		  if (!cgo_lighting){
		    glEnable(GL_LIGHTING);
		  }
		}
	      }
          }
        }
      }
    }
  }
}


/*========================================================================*/
ObjectCGO *ObjectCGONew(PyMOLGlobals * G)
{
  OOAlloc(G, ObjectCGO);

  ObjectInit(G, (CObject *) I);

  I->State = VLAMalloc(10, sizeof(ObjectCGOState), 5, true);
  I->NState = 0;
  I->Obj.type = cObjectCGO;
  I->Obj.fFree = (void (*)(CObject *)) ObjectCGOFree;
  I->Obj.fUpdate = (void (*)(CObject *)) ObjectCGOUpdate;
  I->Obj.fInvalidate = (void (*)(CObject *, int rep, int level, int state))
    ObjectCGOInvalidate;
  I->Obj.fRender = (void (*)(CObject *, RenderInfo *)) ObjectCGORender;
  I->Obj.fGetNFrame = (int (*)(CObject *)) ObjectCGOGetNState;

  return (I);
}

#ifndef _PYMOL_NOPY

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
#ifdef _PYMOL_CGO_DRAWARRAYS
	  {
	    CGO *convertcgo = NULL;
	    if(cgo && cgo->has_begin_end){
	      convertcgo = CGOCombineBeginEnd(cgo, 0);
	      CGOFree(cgo);
	      cgo = convertcgo;
	    }
	  }
#endif
        }
      }
      FreeP(raw);
    }
  }
  return (cgo);
}
#endif


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
  int est = 0;

  if(obj) {
    if(obj->Obj.type != cObjectCGO)     /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = ObjectCGONew(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCGOState, state);
    I->NState = state + 1;
  }

  if(I->State[state].shaderCGO && I->State[state].std!=I->State[state].shaderCGO) {
    CGOFree(I->State[state].shaderCGO);
    I->State[state].shaderCGO = 0;
  }
  if(I->State[state].std) {
    CGOFree(I->State[state].std);
  }
  if(I->State[state].ray) {
    CGOFree(I->State[state].ray);
  }
  if (cgo)
    est = CGOCheckComplex(cgo);
#ifdef _PYMOL_CGO_DRAWARRAYS
  {
    CGO *convertcgo = NULL;
    if(cgo && cgo->has_begin_end){
      convertcgo = CGOCombineBeginEnd(cgo, 0);
      CGOFree(cgo);
      cgo = convertcgo;
    }
  }
#endif
  if(est) {
    I->State[state].ray = cgo;
    I->State[state].std = CGOSimplify(cgo, est);
  } else {
    I->State[state].std = cgo;
  }
  I->State[state].valid = true;

  if(I) {
    ObjectCGORecomputeExtent(I);
  }
  SceneChanged(G);
  SceneCountFrames(G);
  return (I);
}


/*========================================================================*/

ObjectCGO *ObjectCGONewVFontTest(PyMOLGlobals * G, char *text, float *pos)
{

  ObjectCGO *I = NULL;
  int font_id;
  CGO *cgo = NULL;
  float scale[2] = { 1.0, 1.0 };

  font_id = VFontLoad(G, 1, 1, 1, true);
  cgo = CGONew(G);
  VFontWriteToCGO(G, font_id, cgo, text, pos, scale, NULL);
  I = ObjectCGOFromCGO(G, NULL, cgo, 0);

  return (I);
}


/*========================================================================*/
ObjectCGO *ObjectCGODefine(PyMOLGlobals * G, ObjectCGO * obj, PyObject * pycgo, int state)
{                               /* assumes blocked interpreter */
#ifdef _PYMOL_NOPY
  return NULL;
#else
  ObjectCGO *I = NULL;

  CGO *cgo, *font_cgo;
  int est;

  if(obj) {
    if(obj->Obj.type != cObjectCGO)     /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = ObjectCGONew(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCGOState, state);
    I->NState = state + 1;
  }

  if(I->State[state].std) {
    CGOFree(I->State[state].std);
    I->State[state].std = NULL;
  }
  if(I->State[state].ray) {
    CGOFree(I->State[state].ray);
    I->State[state].ray = NULL;
  }
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
#ifdef _PYMOL_CGO_DRAWARRAYS
	  {
	    CGO *convertcgo = NULL;
	    if(cgo && cgo->has_begin_end){
	      convertcgo = CGOCombineBeginEnd(cgo, 0);
	      CGOFree(cgo);
	      cgo = convertcgo;
	    }
	  }
#endif
          if(est) {
            I->State[state].ray = cgo;
            I->State[state].std = CGOSimplify(cgo, est);
          } else {
            I->State[state].std = cgo;
          }
          I->State[state].valid = true;
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
#endif
}

ObjectCGO *ObjectCGOFromFloatArray(PyMOLGlobals * G, ObjectCGO * obj,
                                   float *array, int size, int state, int quiet)
{
  ObjectCGO *I = NULL;

  CGO *cgo, *font_cgo;
  int est;

  if(obj) {
    if(obj->Obj.type != cObjectCGO)     /* TODO: handle this */
      obj = NULL;
  }
  if(!obj) {
    I = ObjectCGONew(G);
  } else {
    I = obj;
  }
  if(state < 0)
    state = I->NState;
  if(I->NState <= state) {
    VLACheck(I->State, ObjectCGOState, state);
    I->NState = state + 1;
  }
  if(I->State[state].shaderCGO && I->State[state].std!=I->State[state].shaderCGO) {
    CGOFree(I->State[state].shaderCGO);
    I->State[state].shaderCGO = 0;
  }
  if(I->State[state].std) {
    CGOFree(I->State[state].std);
  }
  if(I->State[state].ray) {
    CGOFree(I->State[state].ray);
  }
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
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      CGO *convertcgo = NULL;
      if(cgo && cgo->has_begin_end){
	convertcgo = CGOCombineBeginEnd(cgo, 0);
	CGOFree(cgo);
	cgo = convertcgo;
      }
    }
#endif
    if(est) {
      I->State[state].ray = cgo;
      I->State[state].std = CGOSimplify(cgo, est);
    } else {
      I->State[state].std = cgo;
    }
    I->State[state].valid = true;
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


