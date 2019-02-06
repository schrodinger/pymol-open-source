
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

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Scene.h"
#include"GadgetSet.h"
#include"Color.h"
#include"PConv.h"
#include"main.h"
#include"CGO.h"
#include"ShaderMgr.h"
#include"Ray.h"

int GadgetSetGetVertex(GadgetSet * I, int index, int base, float *v)
{
  int ok = true;
  float *v0, *v1;
  if(index < I->NCoord) {
    v0 = I->Coord + 3 * index;
    if(base < 0) {
      copy3f(v0, v);
      if(index){
        add3f(I->Coord, v, v);
      }
    } else if(base < I->NCoord) {
      v1 = I->Coord + 3 * base;
      add3f(v1, v0, v);
      if(index)
        add3f(I->Coord, v, v);
    } else {
      ok = false;
    }
  } else
    ok = false;
  return (ok);
}

int GadgetSetSetVertex(GadgetSet * I, int index, int base, float *v)
{
  int ok = true;
  float *v0, *v1;
  if(index < I->NCoord) {
    v0 = I->Coord + 3 * index;
    if(base < 0) {
      copy3f(v, v0);
      if(index){
        subtract3f(v0, I->Coord, v0);
      } else {
	if (I->offsetPtOP)
	  copy3f(v0, &I->StdCGO->op[I->offsetPtOP]);
	if (I->offsetPtOPick)
	  copy3f(v0, &I->PickCGO->op[I->offsetPtOPick]);
      }
    } else if(base < I->NCoord) {
      v1 = I->Coord + 3 * base;
      subtract3f(v, v1, v0);
      if(index)
        subtract3f(v0, I->Coord, v0);
    } else {
      ok = false;
    }
  } else
    ok = false;
  return (ok);
}

int GadgetSetFromPyList(PyMOLGlobals * G, PyObject * list, GadgetSet ** gs, int version)
{
  int ok = true;
  GadgetSet *I = NULL;
  PyObject *tmp = NULL;

  if(*gs) {
    (*gs)->fFree();
    *gs = NULL;
  }

  if(list == Py_None) {         /* allow None for GSet */
    *gs = NULL;
  } else {

    if(ok)
      I = GadgetSetNew(G);
    if(ok)
      ok = (I != NULL);
    if(ok)
      ok = (list != NULL);
    if(ok)
      ok = PyList_Check(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */

    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->NCoord);
    if(ok && I->NCoord)
      ok = PConvPyListToFloatVLA(PyList_GetItem(list, 1), &I->Coord);

    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->NNormal);
    if(ok && I->NNormal)
      ok = PConvPyListToFloatVLA(PyList_GetItem(list, 3), &I->Normal);

    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 4), &I->NColor);
    if(ok && I->NColor)
      ok = PConvPyListToFloatVLA(PyList_GetItem(list, 5), &I->Color);

    if(ok)
      ok = ((tmp = PyList_GetItem(list, 6)) != NULL);
    if(ok && (tmp != Py_None))
      ok = ((I->ShapeCGO = CGONewFromPyList(I->G, tmp, version)) != NULL);

    if(ok)
      ok = ((tmp = PyList_GetItem(list, 7)) != NULL);
    if(ok && (tmp != Py_None))
      ok = ((I->PickShapeCGO = CGONewFromPyList(I->G, tmp, version)) != NULL);

    if(ok && I->ShapeCGO)
      if(CGOCheckForText(I->ShapeCGO)) {
        CGOPreloadFonts(I->ShapeCGO);
      }

    if(!ok) {
      if(I)
        I->fFree();
    } else {
      *gs = I;
    }
  }

  return (ok);
}

PyObject *GadgetSetAsPyList(GadgetSet * I, bool incl_cgos)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(8);

    PyList_SetItem(result, 0, PyInt_FromLong(I->NCoord));

    if(I->NCoord) {
      PyList_SetItem(result, 1, PConvFloatArrayToPyList(I->Coord, I->NCoord * 3));
    } else {
      PyList_SetItem(result, 1, PConvAutoNone(NULL));
    }

    PyList_SetItem(result, 2, PyInt_FromLong(I->NNormal));

    if(I->NNormal) {
      PyList_SetItem(result, 3, PConvFloatArrayToPyList(I->Normal, I->NNormal * 3));
    } else {
      PyList_SetItem(result, 3, PConvAutoNone(NULL));
    }

    PyList_SetItem(result, 4, PyInt_FromLong(I->NColor));

    if(I->NColor) {
      PyList_SetItem(result, 5, PConvFloatArrayToPyList(I->Color, I->NColor));
    } else {
      PyList_SetItem(result, 5, PConvAutoNone(NULL));
    }

    if(incl_cgos && I->ShapeCGO) {
      PyList_SetItem(result, 6, CGOAsPyList(I->ShapeCGO));
    } else {
      PyList_SetItem(result, 6, PConvAutoNone(NULL));
    }

    if(incl_cgos && I->PickShapeCGO) {
      PyList_SetItem(result, 7, CGOAsPyList(I->PickShapeCGO));
    } else {
      PyList_SetItem(result, 7, PConvAutoNone(NULL));
    }

  }
  return (PConvAutoNone(result));
}


/*========================================================================*/
int GadgetSetGetExtent(GadgetSet * I, float *mn, float *mx)
{
  float *v;
  int a;
  v = I->Coord;
  for(a = 0; a < I->NCoord; a++) {
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
  }
  return (I->NCoord);
}


/*========================================================================*/
void GadgetSet::invalidateRep(int type, int level)
{
#if TO_DO
  GadgetSet * I = this;
  int a;
  PRINTFD(I->G, FB_GadgetSet)
    " GadgetSetInvalidateRep: entered.\n" ENDFD;
  if(type >= 0) {
    if(type < I->NRep) {
      SceneChanged(I->G);
      if(I->Rep[type]) {
        I->Rep[type]->fFree(I->Rep[type]);
        I->Rep[type] = NULL;
      }
    }
  } else {
    for(a = 0; a < I->NRep; a++) {
      SceneChanged(I->G);
      if(I->Rep[a]) {
        switch (level) {
        case cRepInvColor:
          if(I->Rep[a]->fRecolor) {
            I->Rep[a]->fInvalidate(I->Rep[a], (struct CoordSet *) I, level);
          } else {
            I->Rep[a]->fFree(I->Rep[a]);
            I->Rep[a] = NULL;
          }
          break;
        default:
          I->Rep[a]->fFree(I->Rep[a]);
          I->Rep[a] = NULL;
          break;
        }
      }
    }
  }
#endif
}


/*========================================================================*/
void GadgetSet::update()
{
  GadgetSet * I = this;
  if(I->StdCGO) {
    CGOFree(I->StdCGO);
    I->offsetPtOP = 0;
    I->StdCGO = NULL;
  }
  if(I->PickCGO) {
    CGOFree(I->PickCGO);
    I->offsetPtOPick = 0;
    I->PickCGO = NULL;
  }
}


/*========================================================================*/
void GadgetSet::render(RenderInfo * info)
{
  GadgetSet * I = this;
  int pass = info->pass;
  CRay *ray = info->ray;
  auto pick = info->pick;
  const float *color;
  PickContext context;

  context.object = I->Obj;
  context.state = I->State;

  color = ColorGet(I->G, I->Obj->Obj.Color);

  if(pass < 0 || ray || pick) {
    PyMOLGlobals *G = I->G;

    if(ray) {
      if(I->ShapeCGO){
	float mat[16] = { 1.f, 0.f, 0.f, I->Coord[0], 
                          0.f, 1.f, 0.f, I->Coord[1], 
                          0.f, 0.f, 1.f, I->Coord[2],
			  0.f, 0.f, 0.f, 1.f };
	RayPushTTT(ray);
	RaySetTTT(ray, true, mat);  /* Used to set the ray-tracing matrix,
				       this works, but is there another way to do this? */
	CGORenderRay(I->ShapeCGO, ray, info, color, NULL, I->Obj->Obj.Setting, NULL);
	RayPopTTT(ray);
      }
    } else if(G->HaveGUI && G->ValidContext) {
      short use_shader = (short) SettingGetGlobal_b(I->G, cSetting_use_shaders);
      if(pick) {
	if (!I->PickCGO && I->PickShapeCGO){
	  CGO *convertcgo;
	  int ok = true;
	  convertcgo = CGOCombineBeginEnd(I->PickShapeCGO, 0);
	  CHECKOK(ok, convertcgo);
	  if (ok){
	    if (use_shader){
	      CGO *tmpCGO;
	      tmpCGO = CGOOptimizeToVBOIndexedNoShader(convertcgo, 0);
	      I->PickCGO = CGONew(G);
	      CGODisable(I->PickCGO, GL_DEPTH_TEST);
	      CGOEnable(I->PickCGO, GL_RAMP_SHADER);
	      I->offsetPtOPick = CGOUniform3f(I->PickCGO, RAMP_OFFSETPT, (const float*)I->Coord);
	      CGOAppendNoStop(I->PickCGO, tmpCGO);
	      CGOFreeWithoutVBOs(tmpCGO);
	      CGODisable(I->PickCGO, GL_RAMP_SHADER);
	      CGOEnable(I->PickCGO, GL_DEPTH_TEST);
	      CGOStop(I->PickCGO);
	      I->PickCGO->use_shader = true;
	      CGOFree(convertcgo);
	    } else {
	      I->PickCGO = convertcgo;
	    }
	  } else {
	    CGOFree(convertcgo);
	  }
	}
        if(I->PickCGO) {
	  if (use_shader){
	    CGORenderGLPicking(I->PickCGO, info, &context, I->Obj->Obj.Setting, NULL);
#ifndef PURE_OPENGL_ES_2
	  } else {
	    glDisable(GL_DEPTH_TEST);
	    glTranslatef(I->Coord[0],I->Coord[1],I->Coord[2]);
	    CGORenderGLPicking(I->PickShapeCGO, info, &context, I->Obj->Obj.Setting, NULL);
	    glTranslatef(-I->Coord[0],-I->Coord[1],-I->Coord[2]);
	    glEnable(GL_DEPTH_TEST);
#endif
	  }
        }
      } else {
	if (!I->StdCGO && I->ShapeCGO){
	  CGO *convertcgo;
	  int ok = true;
	  convertcgo = CGOCombineBeginEnd(I->ShapeCGO, 0);
	  CHECKOK(ok, convertcgo);
	  if (ok){
	    if (use_shader){
	      CGO *tmpCGO;
	      tmpCGO = CGOOptimizeToVBONotIndexedNoShader(convertcgo, 0);
	      I->StdCGO = CGONew(G);
	      CGODisable(I->StdCGO, GL_DEPTH_TEST);
	      CGOEnable(I->StdCGO, GL_RAMP_SHADER);
	      I->offsetPtOP = CGOUniform3f(I->StdCGO, RAMP_OFFSETPT,  (const float*)I->Coord);
	      CGOAppendNoStop(I->StdCGO, tmpCGO);
	      CGOFreeWithoutVBOs(tmpCGO);
	      CGODisable(I->StdCGO, GL_RAMP_SHADER);
	      CGOEnable(I->StdCGO, GL_DEPTH_TEST);
	      CGOStop(I->StdCGO);
	      I->StdCGO->use_shader = true;
	      CGOFree(convertcgo);
	    } else {
	      I->StdCGO = convertcgo;
	    }
	  } else {
	    CGOFree(convertcgo);
	  }
	}
        if(I->StdCGO) {
	  if (use_shader){
	    if (color)
              CGORenderGL(I->StdCGO, NULL, I->Obj->Obj.Setting, NULL, info, NULL);
#ifndef PURE_OPENGL_ES_2
	  } else {
	    glDisable(GL_DEPTH_TEST);
	    glTranslatef(I->Coord[0],I->Coord[1],I->Coord[2]);
	    CGORenderGL(I->ShapeCGO, NULL, I->Obj->Obj.Setting, NULL, info, NULL);
	    glTranslatef(-I->Coord[0],-I->Coord[1],-I->Coord[2]);
	    glEnable(GL_DEPTH_TEST);
#endif
	  }
        }
      }
    }
  }
}


/*========================================================================*/
GadgetSet *GadgetSetNew(PyMOLGlobals * G)
{
  OOAlloc(G, GadgetSet);
  I->G = G;
  I->NCoord = 0;
  I->NColor = 0;
  I->NNormal = 0;
  I->Coord = NULL;
  I->Normal = NULL;
  I->Color = NULL;
  I->Setting = NULL;
  I->PickCGO = NULL;
  I->StdCGO = NULL;
  I->ShapeCGO = NULL;
  I->PickShapeCGO = NULL;
  I->offsetPtOP = 0;
  I->offsetPtOPick = 0;

  return (I);
}


/*========================================================================*/
void GadgetSet::fFree()
{
  GadgetSet * I = this;
  if(I) {
    CGOFree(I->PickCGO);
    CGOFree(I->PickShapeCGO);
    CGOFree(I->StdCGO);
    CGOFree(I->ShapeCGO);
    I->offsetPtOP = 0;
    I->offsetPtOPick = 0;
    VLAFreeP(I->Coord);
    VLAFreeP(I->Normal);
    VLAFreeP(I->Color);
    OOFreeP(I);
  }
}
