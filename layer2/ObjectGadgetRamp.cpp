
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
#include"ObjectGadgetRamp.h"
#include"GadgetSet.h"
#include"Base.h"
#include"MemoryDebug.h"
#include"CGO.h"
#include"Scene.h"
#include"Setting.h"
#include"PConv.h"
#include"main.h"
#include"Color.h"
#include"VFont.h"
#include"ObjectMolecule.h"
#include"Executive.h"
#include"Util.h"
#include"P.h"
#include"PyMOLObject.h"

static void ObjectGadgetRampBuild(ObjectGadgetRamp *);
static int ObjectGadgetRampHandleInputColors(ObjectGadgetRamp *);

ObjectGadgetRamp::~ObjectGadgetRamp()
{
  auto I = this;
  ColorForgetExt(I->G, I->Name);
}

static void ObjectGadgetRampCalculate(ObjectGadgetRamp * I, float v, float *result)
{
  const float _1 = 1.0F;
  const float _0 = 0.0F;
  /* from Filipe Maia */

  /* All of this functions are taken right of the gnuplot manual  */
  if(v > _1)
    v = _1;
  else if(v < _0)
    v = _0;

  switch (I->CalcMode) {
  case cRAMP_TRADITIONAL:
    result[0] = (float) sqrt(v);
    result[1] = v * v * v;
    result[2] = (float) sin(v * 2 * cPI);
    break;
  case cRAMP_SLUDGE:
    result[0] = v;
    result[1] = (float) fabs(v - 0.5F);
    result[2] = v * v * v * v;
    break;
  case cRAMP_OCEAN:
    result[0] = 3 * v - 2;
    result[1] = (float) fabs((3 * v - 1) / 2);
    result[2] = v;
    break;
  case cRAMP_HOT:
    result[0] = 3 * v;
    result[1] = 3 * v - 1;
    result[2] = 3 * v - 2;
    break;
  case cRAMP_GRAYABLE:
    result[0] = v / 0.32F - 0.78125F;
    result[1] = 2 * v - 0.84F;
    result[2] = v / 0.08F - 11.5F;      /* I'm not so sure about this one */
    break;
  case cRAMP_RAINBOW:
    result[0] = (float) fabs(2 * v - 0.5F);
    result[1] = (float) sin(v * cPI);
    result[2] = (float) cos(v * cPI / 2.0F);
    break;
  case cRAMP_AFMHOT:
    result[0] = 2 * v;
    result[1] = 2 * v - 0.5F;
    result[2] = 2 * v - 1.0F;
    break;
  case cRAMP_GRAYSCALE:
    result[0] = v;
    result[1] = v;
    result[2] = v;
    break;
  default:                     /* default is simply white */
    result[0] = 1.0F;
    result[1] = 1.0F;
    result[2] = 1.0F;
    break;
  }
  clamp3f(result);
}

/**
 * Get the "Level" array, eventually spaced out to match number of colors.
 */
static const float * ObjectGadgetRampGetLevel(ObjectGadgetRamp * I) {
  if (!I->Level || !I->Color)
    return I->Level;

  int n_color = VLAGetSize(I->Color) / 3;
  if (n_color == I->NLevel || n_color < 2)
    return I->Level;

  if (!I->LevelTmp) {
    float level0 = I->Level[0], level1 = I->Level[I->NLevel - 1];
    I->LevelTmp = pymol::vla<float>(n_color);
    for (int i = 0; i < n_color; ++i) {
      float a = i / (float) (n_color - 1);
      I->LevelTmp[i] = level0 * (1.f - a) + level1 * a;
    }
  }

  return I->LevelTmp;
}

/**
 * We support storing special color indices as negative R in the RGB color
 */
inline int GetSpecial(const float * rgb) {
  if (rgb[0] < 0.f)
    return (int) rgb[0];
  return 0;
}

static int _ObjectGadgetRampInterpolate(ObjectGadgetRamp * I, float level, float *color,
                                        const float *table)
{
  const float *i_level = ObjectGadgetRampGetLevel(I);
  int n_level = VLAGetSize(i_level);
  const float _0 = 0.0F;
  const float _1 = 1.0F;
  int ok = true;
  if(i_level && table) {
    int level_is_ge = -1;
    int level_is_le = n_level;
    int i = 0;
    i = n_level - 1;
    while(i >= 0) {
      float f = i_level[i];
      if(level >= f) {
        level_is_ge = i;
        break;
      }
      i--;
    }
    i = 0;
    while(i < n_level) {
      float f = i_level[i];
      if(level <= f) {
        level_is_le = i;
        break;
      } else
        i++;
    }
    if(level_is_ge != level_is_le) {
      if(level_is_le == 0) {    /* lower extreme */
        copy3f(table, color);
      } else if(level_is_ge == (n_level - 1)) { /* upper extreme */
        copy3f(table + 3 * (n_level - 1), color);
      } else {
        float d, x0, x1;

        d = i_level[level_is_ge] - i_level[level_is_le];
        if(fabs(d) > R_SMALL8) {
          x0 = (level - i_level[level_is_le]) / d;
          x1 = 1.0F - x0;
          for(i = 0; i < 3; i++) {
            color[i] = x0 * table[3 * level_is_ge + i] + x1 * table[3 * level_is_le + i];
          }
          clamp3f(color);
        } else {
          copy3f(table + 3 * level_is_ge, color);
        }
      }
    } else {                    /* dead on the specified level */
      copy3f(table + 3 * level_is_ge, color);
      clamp3f(color);
    }
  } else {
    float base, range;
    if(n_level && i_level) {
      base = i_level[0];
      range = i_level[n_level - 1] - base;
      if(fabs(range) < R_SMALL8)
        range = _1;
    } else {
      base = _0;
      range = _1;
    }
    level = (level - base) / range;
    ObjectGadgetRampCalculate(I, level, color);
  }
  return (ok);
}

static int _ObjectGadgetRampBlend(ObjectGadgetRamp * I, float *color,
                                  const float *table, int mode)
{
  /* this capability needs to be re-thought */

  const float *i_level = ObjectGadgetRampGetLevel(I);
  int n_level = VLAGetSize(i_level);
  const float _1 = 1.0F;
  float avg[3];
  int ok = true;
  zero3f(avg);

  {
    int cnt = 0;
    switch (mode) {
    case 1:
    case 2:
      break;
    default:
      if(i_level && table) {
        int i;
        for(i = 0; i < n_level; i++) {
          add3f(table + 3 * i, avg, avg);
          cnt++;
        }
        if(cnt) {
          float fact = _1 / cnt;
          scale3f(avg, fact, avg);
        }
        clamp3f(avg);
      }
      copy3f(avg, color);
    }
  }

  switch (mode) {
  case 1:                      /* min components */
  case 3:
    ones3f(color);
    if(i_level && table) {
      int i, j;
      for(i = 0; i < n_level; i++) {
        for(j = 0; j < 3; j++) {
          color[j] = (color[j] < table[3 * i + j]) ? color[j] : table[3 * i + j];
        }
      }
      clamp3f(color);
    }
    if(mode == 3) {             /* average serves as a minimum */
      int j;
      for(j = 0; j < 3; j++) {
        color[j] = (color[j] > avg[j]) ? color[j] : avg[j];
      }
    }
    break;
  case 2:                      /* max components */
    zero3f(color);
    if(i_level && table) {
      int i, j;
      for(i = 0; i < n_level; i++) {
        for(j = 0; j < 3; j++) {
          color[j] = (color[j] > table[3 * i + j]) ? color[j] : table[3 * i + j];
        }
      }
      clamp3f(color);
    }
    break;
  default:                     /* simple average of all colors */
    copy3f(avg, color);
    break;
  }
  return (ok);
}

#define MAX_COLORS 64

static int ObjectGadgetRampInterpolateWithSpecial(ObjectGadgetRamp * I,
                                                  float level,
                                                  float *color,
                                                  const float *atomic,
                                                  const float *object,
                                                  const float *vertex, int state, int blend_all)
{
  /* now thread-safe...via stack copy of colors */

  float stack_color[MAX_COLORS * 3];
  const float *i_level = ObjectGadgetRampGetLevel(I);
  const float *i_color = I->Color;

  if(i_level && i_color) {
    int i = 0;
    int n_level = VLAGetSize(i_level);
    /* mix special coloring into the table */
    const float *src = i_color;
    float *dst = stack_color;
    if((n_level + 2) > MAX_COLORS)
      n_level = MAX_COLORS - 2;
    while(i < n_level) {
      int index = GetSpecial(src);
      switch (index) {
        case 0:
          copy3f(src, dst);
          break;
        case cColorDefault:
        case cColorAtomic:
          copy3f(atomic, dst);
          break;
        case cColorObject:
          copy3f(object, dst);
          break;
        default:               /* allow nested ramps */
          ColorGetRamped(I->G, index, vertex, dst, state);
          break;
      }
      dst += 3;
      src += 3;
      i++;
    }
    i_color = stack_color;
  }

  /* interpolate using static tables */
  if(blend_all)
    return _ObjectGadgetRampBlend(I, color, i_color, I->SrcState);
  return _ObjectGadgetRampInterpolate(I, level, color, i_color);
}

int ObjectGadgetRampInterpolate(ObjectGadgetRamp * I, float level, float *color)
{
  int result = _ObjectGadgetRampInterpolate(I, level, color, I->Color);
  return result;
}

PyObject *ObjectGadgetRampAsPyList(ObjectGadgetRamp * I)
{
  PyObject *result = NULL;

  result = PyList_New(11);

  PyList_SetItem(result, 0, ObjectGadgetPlainAsPyList(I, false));
  PyList_SetItem(result, 1, PyInt_FromLong(I->RampType));
  PyList_SetItem(result, 2, PyInt_FromLong(I->NLevel));
  if(I->Level && I->NLevel) {
    PyList_SetItem(result, 3, PConvFloatVLAToPyList(I->Level));
  } else {
    PyList_SetItem(result, 3, PConvAutoNone(NULL));
  }
  if(I->Color && I->NLevel) {
    PyList_SetItem(result, 4, PConvFloatVLAToPyList(I->Color));
  } else {
    PyList_SetItem(result, 4, PConvAutoNone(NULL));
  }
  PyList_SetItem(result, 5, PyInt_FromLong(I->var_index));
  PyList_SetItem(result, 6, PyString_FromString(I->SrcName));
  PyList_SetItem(result, 7, PyInt_FromLong(I->SrcState));
  PyList_SetItem(result, 8, PyInt_FromLong(I->CalcMode));

  // I->Special, removed in PyMOL 1.8
  bool any = false;
  int* special = NULL;
  int pse_export_version = SettingGetGlobal_f(I->G, cSetting_pse_export_version) * 1000;
  if (I->Color && pse_export_version < 1800) {
    int n_color = VLAGetSize(I->Color) / 3;
    special = VLAlloc(int, n_color);
    for (int a = 0; a < n_color; ++a) {
      any = (special[a] = GetSpecial(I->Color + a * 3)) || any;
  }
  }
  PyList_SetItem(result, 9, any ? PConvIntVLAToPyList(special) : PConvAutoNone(NULL));
  VLAFreeP(special);

  PyList_SetItem(result, 10, PConvAutoNone(NULL) /* I->Extreme, removed in PyMOL 1.8 */);
  return (PConvAutoNone(result));
}

int ObjectGadgetRampNewFromPyList(PyMOLGlobals * G, PyObject * list,
                                  ObjectGadgetRamp ** result, int version)
{
  ObjectGadgetRamp *I = NULL;
  int ok = true;
  int ll = 0;

  if(ok)
    I = new ObjectGadgetRamp(G);
  if(ok)
    ok = (I != NULL);
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  if(ok){
    auto val = PyList_GetItem(list, 0);
    ok = ObjectGadgetInitFromPyList(G, val, I, version);
  }
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->RampType);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 2), &I->NLevel);
  if(ok && I->NLevel)
    ok = PConvPyListToFloatVLA(PyList_GetItem(list, 3), &I->Level);
  if(ok && I->NLevel) {
    PyObject *item = PyList_GetItem(list, 4);
    if(item != Py_None) {
      ok = PConvPyListToFloatVLA(item, &I->Color);
    }
  }
  if(ok)
    ok = CPythonVal_PConvPyStrToStr_From_List(G, list, 6, I->SrcName, WordLength);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 7), &I->SrcState);
  if(ok && (ll > 8))
    ok = CPythonVal_PConvPyIntToInt_From_List(G, list, 8, &I->CalcMode);

  if(ok && I->NLevel && (ll > 10)) {
    CPythonVal *item = CPythonVal_PyList_GetItem(G, list, 10);
    if(!CPythonVal_IsNone(item)) {
      // ObjectGadgetRamp::Extreme removed in PyMOL 1.8
      // Copy extreme colors to the beginning and end of the Color array
      // and repeat the first and last level value
      float *extreme = NULL;
      PConvPyListToFloatVLA(item, &extreme);
      if (extreme) {
        I->NLevel += 2;

        VLASize(I->Level, float, I->NLevel);
        for (int i = I->NLevel - 2; i >= 1; --i)
          I->Level[i] = I->Level[i - 1];
        I->Level[I->NLevel - 1] = I->Level[I->NLevel - 2];

        if (I->Color) {
          VLASize(I->Color, float, I->NLevel * 3);
          for (int i = (I->NLevel - 1) * 3 - 1; i >= 3; --i)
            I->Color[i] = I->Color[i - 3];
          copy3f(extreme, I->Color.data());
          copy3f(extreme + 3, I->Color + (I->NLevel - 1) * 3);
        }
        VLAFreeP(extreme);
      }
    }
    CPythonVal_Free(item);
  }

  ObjectGadgetRampHandleInputColors(I);
  ObjectGadgetRampBuild(I);

  if(ok)
    (*result) = I;
  return (ok);
}

int ObjectGadgetRampInterVertex(ObjectGadgetRamp * I, const float *pos, float *color, int state)
{
  float level;
  int ok = true;
  switch (I->RampType) {
  case cRampMap:
    if(!I->Map)
      I->Map = ExecutiveFindObjectMapByName(I->G, I->SrcName);
    if(!ExecutiveValidateObjectPtr(I->G, I->Map, cObjectMap))
      ok = false;
    else {
      int src_state;
      if(I->SrcState >= 0)
        src_state = I->SrcState;
      else
        src_state = state;
      if(src_state < 0)
        src_state = SceneGetState(I->G);
      if(ok)
        ok = (I->Map != NULL);
      if(ok)
        ok = ObjectMapInterpolate(I->Map, src_state, pos, &level, NULL, 1);
      if(ok)
        ok = ObjectGadgetRampInterpolate(I, level, color);
    }
    break;
  case cRampMol:
    if(!I->Mol)
      I->Mol = ExecutiveFindObjectMoleculeByName(I->G, I->SrcName);
    if(!ExecutiveValidateObjectPtr(I->G, I->Mol, cObjectMolecule))
      ok = false;
    else {
      float cutoff = 1.0F;
      float dist;
      int sub_vdw = false;
      if(state < 0)
        state = SceneGetState(I->G);
      if(I->Level && I->NLevel) {
        cutoff = I->Level[I->NLevel - 1];
        if(I->Level[0] < 0.0F) {
          sub_vdw = true;
          cutoff += MAX_VDW;
        }
      }
      if(ok)
        ok = (I->Mol != NULL);
      if(ok) {
	if (I->Mol->NCSet==1) // if only one state, then set state to 0
	  state = 0;
        if(SettingGet_b
           (I->G, I->Setting.get(), NULL,
            cSetting_ramp_blend_nearby_colors)) {
          float atomic[3];
          int index =
            ObjectMoleculeGetNearestBlendedColor(I->Mol, pos, cutoff, state, &dist,
                                                 atomic, sub_vdw);
          if(index >= 0) {
            const float *object = ColorGetRaw(I->G, I->Mol->Color);

            if(!ObjectGadgetRampInterpolateWithSpecial(I, dist, color, atomic,
                                                       object, pos, state, false)) {
              copy3f(I->Color.data(), color);
            }
          } else {
            float white[3] = { 1.0F, 1.0F, 1.0F };
            if(!ObjectGadgetRampInterpolateWithSpecial(I, cutoff + 1.0F, color, white,
                                                       white, pos, state, false)) {
              copy3f(I->Color.data(), color);
            }
          }
        } else {
          int index =
            ObjectMoleculeGetNearestAtomIndex(I->Mol, pos, cutoff, state, &dist);
          if(index >= 0) {
            const float *atomic = ColorGetRaw(I->G, I->Mol->AtomInfo[index].color);
            const float *object = ColorGetRaw(I->G, I->Mol->Color);

            if(sub_vdw) {
              dist -= I->Mol->AtomInfo[index].vdw;
              if(dist < 0.0F)
                dist = 0.0F;
            }

            if(!ObjectGadgetRampInterpolateWithSpecial(I, dist, color, atomic,
                                                       object, pos, state, false)) {
              copy3f(I->Color.data(), color);
            }
          } else {
            float white[3] = { 1.0F, 1.0F, 1.0F };
            if(!ObjectGadgetRampInterpolateWithSpecial(I, cutoff + 1.0F, color, white,
                                                       white, pos, state, false)) {
              copy3f(I->Color.data(), color);
            }
          }
        }
      }
    }
    break;
  case cRampNone:
    {
      float white[3] = { 1.0F, 1.0F, 1.0F };
      if(!ObjectGadgetRampInterpolateWithSpecial(I, 0.0F, color, white, white, pos, state, true)) {     /* simple blend */
        copy3f(I->Color.data(), color);
      }
    }
    break;
  default:
    ok = false;
    break;
  }
  return (ok);
}

static void ObjectGadgetRampUpdateCGO(ObjectGadgetRamp * I, GadgetSet * gs)
{
  CGO *cgo;
  int a;
  char buffer[255];
  int blocked = false;
  int font_id = 0;
  int n_color = I->Color ? VLAGetSize(I->Color) / 3 : 0;

  blocked = PAutoBlock(I->G);
  font_id = VFontLoad(I->G, 1.0, 1, 1, true);
  if(blocked)
    PUnblock(I->G);

  cgo = CGONewSized(I->G, 100);

  /* behind text */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOColor(cgo, 0.05F, 0.05F, 0.05F);
  CGONormal(cgo, 0.f, 0.f, 1.f); // normal 2
  CGOVertex(cgo, I->border, -(I->border + I->bar_height), I->border); // 9
  CGOVertex(cgo, I->border, -(I->height + I->border), I->border); // 7
  CGOVertex(cgo, I->width + I->border, -(I->border + I->bar_height), I->border); // 10
  CGOVertex(cgo, I->width + I->border, -(I->height + I->border), I->border); // 8
  CGOEnd(cgo);


  CGOColor(cgo, 1.0F, 1.0F, 1.0F);
  CGONormal(cgo, 0.f, 0.f, 1.f); // text/characters need z-normal to include normals in vbo
  if(I->Level && I->NLevel) {
    float exindent = (n_color > 0) ? I->bar_height : 0.f;
    float pos[] = { I->border + I->text_border + exindent,
                    I->text_border - (I->border + I->height),
		    I->border + I->text_raise };
    float scale[] = { I->text_scale_h, I->text_scale_v };
    float axes[] = { 1.0F, 0.0F, 0.0F,
		     0.0F, 1.0F, 0.0F,
		     0.0F, 0.0F, 1.0F };
    float color[] = { 1.f, 1.f, 1.f };
    /* left text for ramp */
    sprintf(buffer, "%0.3f", I->Level[0]);
    VFontWriteToCGO(I->G, font_id, cgo, buffer, pos, scale, axes, color);

    /* right text, right justified for ramp */
    pos[0] = I->width + I->border - exindent;
    pos[1] = I->text_border - (I->border + I->height);
    pos[2] = I->border + I->text_raise ;
    sprintf(buffer, "%0.3f", I->Level[I->NLevel - 1]);
    /* indent for right justification */
    VFontIndent(I->G, font_id, buffer, pos, scale, axes, -1.f);
    VFontWriteToCGO(I->G, font_id, cgo, buffer, pos, scale, axes, color);
  }

  /* center */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 0.f, 0.f, 1.f); // normal 2
  if(n_color > 0) {
    const float *i_level = ObjectGadgetRampGetLevel(I);
    const float *src = I->Color;
    float stack_color[6], stack_level[2] = {0.f, 1.f};

    if(n_color == 1) {
      n_color = 2;
      copy3f(src, stack_color);
      copy3f(src, stack_color + 3);
      src = stack_color;
      i_level = stack_level;
    }

    // can this ever happen?
    if (!i_level) {
      i_level = stack_level;
      n_color = 2;
    }

    float range = i_level[n_color - 1] - i_level[0];

    // avoid devide by zero
    if (fabs(range) < R_SMALL8) {
      range = 1.f;
      i_level = stack_level;
      n_color = 2;
    }

    // repeat first and last color in a `bar_height` wide square
    for(a = -1; a <= n_color; a++) {
      float tmp[3] = {1.f, 1.f, 1.f};
      float v1 = I->border;

      if(!GetSpecial(src)) {
        copy3f(src, tmp);
        ColorLookupColor(I->G, tmp);
      }

      if (a == n_color) {
        v1 += I->width;
      } else if (a != -1) {
        v1 += I->bar_height + (I->width - 2 * I->bar_height) * (i_level[a] - i_level[0]) / range;
        if (a != n_color - 1) {
          src += 3;
        }
      }

      CGOColorv(cgo, tmp);
      CGOVertex(cgo, v1, -I->border,                 I->border);
      CGOVertex(cgo, v1, -I->border - I->bar_height, I->border);
    }
  } else {
    int samples = 20;
    float fxn;
    float color[3];

    for(a = 0; a < samples; a++) {
      fxn = a / (samples - 1.0F);

      ObjectGadgetRampCalculate(I, fxn, color);
      CGOColorv(cgo, color);

      CGOVertex(cgo, I->border + (I->width * fxn), -I->border, I->border);
      CGOVertex(cgo, I->border + (I->width * fxn), -(I->border + I->bar_height), I->border);
    }
  }
  /* center */
  CGOEnd(cgo);

  CGOColor(cgo, 1.F, 1.F, 1.F);

  /* top */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);

  CGONormal(cgo, 0.f, 0.f, 1.f); // normal 2
  CGOVertex(cgo, I->border, -I->border, I->border); // 5
  CGOVertex(cgo, I->width + I->border, -I->border, I->border); // 6
  CGONormal(cgo, 0.f, 1.f, .1f); // normal 1
  CGOVertex(cgo, 0.0, 0.0, 0.0); // 1
  CGOVertex(cgo, I->width + I->border * 2, 0.0, 0.0); // 2
  CGOEnd(cgo);

  /* bottom */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 0.f, -1.f, .1f); // normal 4
  CGOVertex(cgo, 0.f, -(I->height + I->border * 2), 0.f); // 3
  CGOVertex(cgo, I->width + I->border * 2, -(I->height + I->border * 2), 0.f); // 4
  CGONormal(cgo, 0.f, 0.f, 1.f); // normal 2
  CGOVertex(cgo, I->border, -(I->height + I->border), I->border); // 7
  CGOVertex(cgo, I->width + I->border, -(I->height + I->border), I->border); // 8
  CGOEnd(cgo);

  /* left */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, -1.f, 0.f, 0.1f); // normal 3
  CGOVertex(cgo, 0.f, 0.f, 0.f); // 1
  CGOVertex(cgo, 0.f, -(I->height + I->border * 2), 0.f); // 3
  CGONormal(cgo, 0.f, 0.f, 1.f); // normal 2
  CGOVertex(cgo, I->border, -I->border, I->border); // 5
  CGOVertex(cgo, I->border, -(I->height + I->border), I->border); // 7
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGONormal(cgo, 0.f, 0.f, 1.f); // normal 2
  CGOVertex(cgo, I->width + I->border, -I->border, I->border); // 6
  CGOVertex(cgo, I->width + I->border, -(I->height + I->border), I->border); // 8
  CGONormal(cgo, 1.f, 0.f, .1f); // normal 0
  CGOVertex(cgo, I->width + I->border * 2, 0.f, 0.f); // 2
  CGOVertex(cgo, I->width + I->border * 2, -(I->height + I->border * 2), 0.f); // 4
  CGOEnd(cgo);
  CGOStop(cgo);

  CGOFree(gs->ShapeCGO);
  gs->ShapeCGO = cgo;

  //#ifndef _PYMOL_NOPY
  CGOPreloadFonts(gs->ShapeCGO);
  //#endif

  cgo = CGONewSized(I->G, 100);
  CGODotwidth(cgo, 5);
  CGOPickColor(cgo, 0, cPickableGadget);

  /* top */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 0.f, 0.f, 0.f); // 1
  CGOVertex(cgo, I->width + I->border * 2, 0.f, 0.f); // 2
  CGOVertex(cgo, I->border, -I->border, I->border); // 5
  CGOVertex(cgo, I->width + I->border, -I->border, I->border); // 6
  CGOEnd(cgo);

  /* bottom */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 0.f, -(I->height + I->border * 2), 0.f); // 3
  CGOVertex(cgo, I->width + I->border * 2, -(I->height + I->border * 2), 0.f); // 4
  CGOVertex(cgo, I->border, -(I->height + I->border), I->border); // 7
  CGOVertex(cgo, I->width + I->border, -(I->height + I->border), I->border); // 8
  CGOEnd(cgo);

  /* left */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, 0.f, 0.f, 0.f); // 1
  CGOVertex(cgo, 0.f, -(I->height + I->border * 2), 0.f); // 3
  CGOVertex(cgo, I->border, -I->border, I->border); // 5
  CGOVertex(cgo, I->border, -(I->height + I->border), I->border); // 7
  CGOEnd(cgo);

  /* right */
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, I->width + I->border, -I->border, I->border); // 6
  CGOVertex(cgo, I->width + I->border, -(I->height + I->border), I->border); // 8
  CGOVertex(cgo, I->width + I->border * 2, 0.f, 0.f); // 2
  CGOVertex(cgo, I->width + I->border * 2, -(I->height + I->border * 2), 0.f); // 4
  CGOEnd(cgo);

  /* band */
  CGOPickColor(cgo, 1, cPickableGadget);
  CGOBegin(cgo, GL_TRIANGLE_STRIP);
  CGOVertex(cgo, I->border, -I->border, I->border); // 5
  CGOVertex(cgo, I->width + I->border, -I->border, I->border); // 6
  CGOVertex(cgo, I->border, -(I->height + I->border), I->border); // 7
  CGOVertex(cgo, I->width + I->border, -(I->height + I->border), I->border); // 8
  CGOEnd(cgo);
  CGOStop(cgo);

  CGOFree(gs->PickShapeCGO);

  gs->PickShapeCGO = cgo;
}

static void ObjectGadgetRampBuild(ObjectGadgetRamp * I)
{
  GadgetSet *gs = NULL;
  ObjectGadget *og;

  OrthoBusyPrime(I->G);

  og = I;
  gs = GadgetSetNew(I->G);

  gs->NCoord = 2;
  I->var_index = gs->NCoord;
  gs->Coord = VLAlloc(float, gs->NCoord * 3);

  gs->Coord[0] = I->x;
  gs->Coord[1] = I->y;
  gs->Coord[2] = .3f;
  gs->Coord[3] = gs->Coord[4] = gs->Coord[5] = 0.f;
  gs->NNormal = 0;
  gs->Normal = NULL;

  // was a memory leak in < 1.8.3.1
  for (int i = 0; i < og->NGSet; ++i) {
    delete og->GSet[i];
    og->GSet[i] = NULL;
  }

  og->GSet[0] = gs;
  og->NGSet = 1;
  gs->Obj = (ObjectGadget *) I;
  gs->State = 0;

  ObjectGadgetRampUpdateCGO(I, gs);
  gs->update();

}


/*========================================================================*/
void ObjectGadgetRamp::update()
{
  auto I = this;
  float scale;

  if(I->Changed) {
    scale = (1.0F + 5 * I->GSet[0]->Coord[1 * 3]);
    I->GSet[0]->Coord[1 * 3] = 0.0;
    switch (I->RampType) {
    case cRampMol:
      {
        int a;
        for(a = 0; a < I->NLevel; a++) {
          I->Level[a] = I->Level[a] * scale;
        }
        ExecutiveInvalidateRep(I->G, cKeywordAll, cRepAll, cRepInvColor);
      }
      break;
    default:
      if(I->NLevel == 2) {
        {
          float mean = (I->Level[0] + I->Level[1]) / 2.0F;
          I->Level[0] = (I->Level[0] - mean) * scale + mean;
          I->Level[1] = (I->Level[1] - mean) * scale + mean;
          ExecutiveInvalidateRep(I->G, cKeywordAll, cRepAll, cRepInvColor);
        }
        break;
      } else if(I->NLevel == 3) {
        I->Level[0] = (I->Level[0] - I->Level[1]) * scale + I->Level[1];
        I->Level[2] = (I->Level[2] - I->Level[1]) * scale + I->Level[1];
        ExecutiveInvalidateRep(I->G, cKeywordAll, cRepAll, cRepInvColor);
      }
    }
    VLAFreeP(I->LevelTmp);
    if(I->NGSet)
      if(I->GSet[0]) {
        ObjectGadgetRampUpdateCGO(I, I->GSet[0]);
        ObjectGadgetUpdateStates(I);
      }
    ObjectGadgetUpdateExtents(I);
    I->Changed = false;
    SceneChanged(I->G);
  }
}

static int ObjectGadgetRampHandleInputColors(ObjectGadgetRamp * I)
{
  VLAFreeP(I->LevelTmp);

  if(I->NLevel < 1) {
    VLASize(I->Level, float, 1);
    I->NLevel = 1;
    I->Level[0] = 0.0F;
  }
  if(I->Color) {
    int n_color = VLAGetSize(I->Color) / 3;

    /* handle cases where number of colors doesn't make expectation */
    if(!n_color) {
      VLASize(I->Color, float, 3);
      I->Color[0] = I->Color[1] = I->Color[2] = 1.0F;
      n_color++;
    }

    if(n_color != I->NLevel && I->NLevel != 2) {
      PRINTFB(I->G, FB_ObjectGadget, FB_Warnings)
        " GadgetRamp-Warning: number of colors (%d) and number of levels (%d) don't\n"
        " match and n_level != 2. Support for trailing extreme colors dropped in 1.8.",
        n_color, I->NLevel ENDFB(I->G);
    }

    if(n_color < I->NLevel) {
      // repeat last color until n_color == NLevel
      VLASize(I->Color, float, 3 * I->NLevel);
      while(n_color < I->NLevel) {
        copy3f(I->Color + 3 * (n_color - 1), I->Color + 3 * n_color);
        n_color++;
      }
    }
  }

  return true;
}


/*========================================================================*/
ObjectGadgetRamp *ObjectGadgetRampMapNewAsDefined(PyMOLGlobals * G,
                                                  ObjectGadgetRamp *I,
                                                  ObjectMap * map,
                                                  pymol::vla<float>&& level_vla,
                                                  pymol::vla<float>&& color_vla,
                                                  int map_state,
                                                  float *vert_vla, float beyond,
                                                  float within, float sigma, int zero,
                                                  int calc_mode)
{
  if (!I)
    I = new ObjectGadgetRamp(G);

  I->RampType = cRampMap;

  if (color_vla || calc_mode > 0) {
    I->Color = std::move(color_vla);
    I->CalcMode = calc_mode;
  }

  {
    ObjectMapState *ms;
    float tmp_level[3];
    if(map_state < 0)
      map_state = 0;
    if(vert_vla && map && (ms = ObjectMapGetState(map, map_state))) {
      if(ObjectMapStateGetExcludedStats(G, ms, vert_vla, beyond, within, tmp_level)) {
        tmp_level[0] = tmp_level[1] + (tmp_level[0] - tmp_level[1]) * sigma;
        tmp_level[2] = tmp_level[1] + (tmp_level[2] - tmp_level[1]) * sigma;
        if(zero) {
          if(tmp_level[1] < 0.0F) {
            tmp_level[1] = 0.0F;
            tmp_level[2] = -tmp_level[0];
          } else if(tmp_level[1] > 0.0F) {
            tmp_level[1] = 0.0F;
            tmp_level[0] = -tmp_level[2];
          }
        }
      }
      I->Level = pymol::vla<float>(3);
      copy3f(tmp_level, I->Level.data());
      VLAFreeP(level_vla);
    } else if (level_vla) {
      I->Level = std::move(level_vla);
    }
  }

    I->NLevel = VLAGetSize(I->Level);
  ObjectGadgetRampHandleInputColors(I);
  ObjectGadgetRampBuild(I);

  if (map) {
    I->Map = map;
    I->SrcState = map_state;
    UtilNCopy(I->SrcName, map->Name, WordLength);
  }

  return (I);

}


/*========================================================================*/
ObjectGadgetRamp *ObjectGadgetRampMolNewAsDefined(PyMOLGlobals * G,
                                                  ObjectGadgetRamp *I,
                                                  ObjectMolecule * mol,
                                                  pymol::vla<float>&& level_vla,
                                                  pymol::vla<float>&& color_vla,
                                                  int mol_state, int calc_mode)
{
  if (!I) {
    I = new ObjectGadgetRamp(G);
    I->RampType = cRampNone;
    UtilNCopy(I->SrcName, "none", WordLength);
  }

  if(mol) {
    I->RampType = cRampMol;
    I->Mol = mol;
    I->SrcState = mol_state;
    UtilNCopy(I->SrcName, mol->Name, WordLength);
  }

  if (color_vla || calc_mode > 0) {
    I->Color = std::move(color_vla);
    I->CalcMode = calc_mode;
  }

  if(level_vla) {
    I->Level = std::move(level_vla);
    I->NLevel = VLAGetSize(I->Level);
  }

  ObjectGadgetRampHandleInputColors(I);
  ObjectGadgetRampBuild(I);
  return (I);
}

void ObjectGadgetRamp::invalidate(cRep_t rep, cRepInv_t level,
                                       int state)
{
}


/*========================================================================*/
ObjectGadgetRamp::ObjectGadgetRamp(PyMOLGlobals * G) : ObjectGadget(G)
{
  auto I = this;
  I->GadgetType = cGadgetRamp;
  I->SrcName[0] = 0;
}
