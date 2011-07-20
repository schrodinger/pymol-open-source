
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
#include"Ortho.h"
#include"Word.h"

#include"Color.h"
#include"PConv.h"
#include"ObjectGadgetRamp.h"
#include"Util.h"
#include"Executive.h"
#include"MyPNG.h"
#include"Scene.h"

#include"OVContext.h"
#include"OVreturns.h"

static int AutoColor[] = {
  26,                           /* carbon */
  5,                            /* cyan */
  154,                          /* lightmagenta */
  6,                            /* yellow */
  9,                            /* salmon */
  29,                           /* hydrogen */
  11,                           /* slate */
  13,                           /* orange */
  10,                           /* lime */
  5262,                         /* deepteal */
  12,                           /* hotpink */
  36,                           /* yelloworange */
  5271,                         /* violetpurple */
  124,                          /* grey70 */
  17,                           /* marine */
  18,                           /* olive */
  5270,                         /* smudge */
  20,                           /* teal */
  5272,                         /* dirtyviolet */
  52,                           /* wheat */
  5258,                         /* deepsalmon */
  5274,                         /* lightpink */
  5257,                         /* aquamarine */
  5256,                         /* paleyellow */
  15,                           /* limegreen */
  5277,                         /* skyblue */
  5279,                         /* warmpink */
  5276,                         /* limon */
  53,                           /* violet */
  5278,                         /* bluewhite */
  5275,                         /* greencyan */
  5269,                         /* sand */
  22,                           /* forest */
  5266,                         /* lightteal */
  5280,                         /* darksalmon */
  5267,                         /* splitpea */
  5268,                         /* raspberry */
  104,                          /* grey50 */
  23,                           /* deepblue */
  51,                           /* brown */
};

static int nAutoColor = 40;
static void lookup_color(CColor * I, float *in, float *out, int big_endian);

void ColorGetBkrdContColor(PyMOLGlobals * G, float *rgb, int invert_flag)
{
  float *bkrd = SettingGetfv(G, cSetting_bg_rgb);

  if(!invert_flag) {
    if((bkrd[0] + bkrd[1] + bkrd[2]) > 0.5F) {
      rgb[0] = 1.0F;
      rgb[1] = 1.0F;
      rgb[2] = 1.0F;
    } else {
      rgb[0] = 0.0F;
      rgb[1] = 0.0F;
      rgb[2] = 0.0F;
    }
  }

  {
    int a;
    for(a = 0; a < 3; a++)
      if(fabs(bkrd[a] - rgb[a]) < 0.5F) {
        rgb[a] = 1.0F - rgb[a];
        if(fabs(bkrd[a] - rgb[a]) < 0.5F) {
          if(bkrd[a] > 0.5F)
            rgb[a] = 0.0F;
          else
            rgb[a] = 1.0F;
        }
      }
  }
}

unsigned int ColorGet32BitWord(PyMOLGlobals * G, float *rgba)
{
  CColor *I = G->Color;
  unsigned int rc, gc, bc, ac;
  unsigned int result;

  rc = (int) (255 * rgba[0] + 0.49999F);
  gc = (int) (255 * rgba[1] + 0.49999F);
  bc = (int) (255 * rgba[2] + 0.49999F);
  ac = (int) (255 * rgba[3] + 0.49999F);

  if(rc > 255)
    rc = 255;
  if(bc > 255)
    bc = 255;
  if(gc > 255)
    gc = 255;
  if(ac > 255)
    ac = 255;

  if(I->BigEndian) {
    result = (rc << 24) | (gc << 16) | (bc << 8) | ac;
  } else {
    result = (ac << 24) | (bc << 16) | (gc << 8) | rc;
  }
  return result;
}

int ColorGetNext(PyMOLGlobals * G)
{
  int result;
  int next;
  next = (int) SettingGet(G, cSetting_auto_color_next);

  if(next >= nAutoColor)
    next = 0;
  result = AutoColor[next];
  next++;
  if(next >= nAutoColor)
    next = 0;
  SettingSet(G, cSetting_auto_color_next, (float) next);
  return (result);
}

int ColorGetCurrent(PyMOLGlobals * G)
{
  int result;
  int next;
  next = (int) SettingGet(G, cSetting_auto_color_next);
  next--;
  if(next < 0)
    next = (nAutoColor - 1);
  result = AutoColor[next];
  return (result);
}

int ColorCheckRamped(PyMOLGlobals * G, int index)
{
  return (index <= (cColorExtCutoff));
}

ObjectGadgetRamp *ColorGetRamp(PyMOLGlobals * G, int index)
{
  register CColor *I = G->Color;
  ObjectGadgetRamp *result = NULL;
  if(index <= cColorExtCutoff) {
    index = cColorExtCutoff - index;
    if(index < I->NExt) {
      if(!I->Ext[index].Ptr) {
        if(I->Ext[index].Name) {
          char *name = OVLexicon_FetchCString(I->Lex, I->Ext[index].Name);
          I->Ext[index].Ptr = (void *) ExecutiveFindObjectByName(G, name);
        }
      }
      if(I->Ext[index].Ptr)
        result = (ObjectGadgetRamp *) I->Ext[index].Ptr;
    }
  }
  return result;
}

int ColorGetRamped(PyMOLGlobals * G, int index, float *vertex, float *color, int state)
{
  register CColor *I = G->Color;
  int ok = false;
  if(index <= cColorExtCutoff) {
    index = cColorExtCutoff - index;
    if(index < I->NExt) {
      if(!I->Ext[index].Ptr) {
        if(I->Ext[index].Name) {
          char *name = OVLexicon_FetchCString(I->Lex, I->Ext[index].Name);
          I->Ext[index].Ptr = (void *) ExecutiveFindObjectByName(G, name);
        }
      }
      if(I->Ext[index].Ptr)
        ok = ObjectGadgetRampInterVertex((ObjectGadgetRamp *) I->Ext[index].Ptr,
                                         vertex, color, state);
    }

  }
  if(!ok) {
    color[0] = 1.0;
    color[1] = 1.0;
    color[2] = 1.0;
  } else if(I->LUTActive) {
    lookup_color(I, color, color, I->BigEndian);
  }
  return (ok);
}

static int ColorFindExtByName(PyMOLGlobals * G, char *name, int null_okay, int *best)
{
  register CColor *I = G->Color;
  int result = -1;
  int wm;
  int a;
  int mybest;
  if(!best)
    best = &mybest;
  *best = 0;
  for(a = 0; a < I->NExt; a++) {
    int color_lex = I->Ext[a].Name;
    if(color_lex) {
      char *color_name = OVLexicon_FetchCString(I->Lex, color_lex);
      wm = WordMatch(G, name, color_name, true);
      if(wm < 0) {
        if(null_okay || (I->Ext[a].Ptr)) {
          result = a;
          *best = 0;
          break;
        }
      } else if((wm > 0) && ((*best) < wm)) {
        if(null_okay || (I->Ext[a].Ptr)) {
          result = a;
          *best = wm;
        }
      }
    }
  }
  return (result);
}

void ColorRegisterExt(PyMOLGlobals * G, char *name, void *ptr, int type)
{
  register CColor *I = G->Color;
  int a;

  a = ColorFindExtByName(G, name, true, NULL);
  if(a < 0) {
    VLACheck(I->Ext, ExtRec, I->NExt);
    a = I->NExt;
    I->NExt++;
    {
      OVreturn_word result = OVLexicon_GetFromCString(I->Lex, name);
      if(OVreturn_IS_OK(result)) {
        I->Ext[a].Name = result.word;
      } else {
        I->Ext[a].Name = 0;
      }
    }
  }
  if(a >= 0) {
    I->Ext[a].Ptr = ptr;
    I->Ext[a].Type = type;
  }
}

void ColorForgetExt(PyMOLGlobals * G, char *name)
{
  register CColor *I = G->Color;
  int a;
  a = ColorFindExtByName(G, name, true, NULL);

  if(a >= 0) {                  /* currently leaks memory in I->Ext array -- TODO fix */
    if(I->Ext[a].Name) {
      OVLexicon_DecRef(I->Lex, I->Ext[a].Name);
      OVOneToOne_DelForward(I->Idx, I->Ext[a].Name);
    }
    I->Ext[a].Ptr = NULL;
  }
}

PyObject *ColorExtAsPyList(PyMOLGlobals * G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  register CColor *I = G->Color;
  PyObject *result, *list;
  ExtRec *ext;
  int a;

  result = PyList_New(I->NExt);
  ext = I->Ext;
  for(a = 0; a < I->NExt; a++) {
    list = PyList_New(2);
    {
      char *name = OVLexicon_FetchCString(I->Lex, ext->Name);
      PyList_SetItem(list, 0, PyString_FromString(name));
    }
    PyList_SetItem(list, 1, PyInt_FromLong(ext->Type));
    PyList_SetItem(result, a, list);
    ext++;
  }
  return (result);
#endif
}


/*========================================================================*/
PyObject *ColorAsPyList(PyMOLGlobals * G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  register CColor *I = G->Color;
  PyObject *result, *list;
  ColorRec *color;
  int n_custom = 0;
  int a, c;
  color = I->Color;
  for(a = 0; a < I->NColor; a++) {
    if(color->Custom || color->LutColorFlag)
      n_custom++;
    color++;
  }
  result = PyList_New(n_custom);
  c = 0;
  color = I->Color;
  for(a = 0; a < I->NColor; a++) {
    if(color->Custom || color->LutColorFlag) {
      list = PyList_New(7);
      {
        char *name = OVLexicon_FetchCString(I->Lex, color->Name);
        PyList_SetItem(list, 0, PyString_FromString(name));
      }
      PyList_SetItem(list, 1, PyInt_FromLong(a));
      PyList_SetItem(list, 2, PConvFloatArrayToPyList(color->Color, 3));
      PyList_SetItem(list, 3, PyInt_FromLong((int) color->Custom));
      PyList_SetItem(list, 4, PyInt_FromLong((int) color->LutColorFlag));
      PyList_SetItem(list, 5, PConvFloatArrayToPyList(color->LutColor, 3));
      PyList_SetItem(list, 6, PyInt_FromLong((int) color->Fixed));
      PyList_SetItem(result, c, list);
      c++;
    }
    color++;
  }
  return (result);
#endif
}

#if 0

/*========================================================================*/
PyObject *ColorTableAsPyList(PyMOLGlobals * G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  register CColor *I = G->Color;
  PyObject *result;

  result = PyList_New(2);
  PyList_SetItem(result, 0, PyFloat_FromDouble(I->Gamma));
  if(0 && I->ColorTable) {      /* for now, we don't embed the color table due to size... */
    PyList_SetItem(result, 1, PConvIntArrayToPyList((int *) (I->ColorTable), 512 * 512));
  } else {
    PyList_SetItem(result, 1, PConvAutoNone(Py_None));
  }
  return result;
#endif
}
#endif


/*========================================================================*/
int ColorConvertOldSessionIndex(PyMOLGlobals * G, int index)
{
  register CColor *I = G->Color;
  if(index > cColorExtCutoff) {
    if(I->HaveOldSessionColors) {
      ColorRec *col = I->Color + (I->NColor - 1);
      int a;
      for(a = I->NColor - 1; a >= 0; a--) {
        if(index == col->old_session_index) {
          index = a;
          break;
        }
        col--;
      }
    }
  } else if(I->HaveOldSessionExtColors) {
    ExtRec *ext = I->Ext + (I->NExt - 1);
    int a;
    for(a = I->NExt - 1; a >= 0; a--) {
      if(index == ext->old_session_index) {
        index = cColorExtCutoff - a;
        break;
      }
      ext--;
    }
  }
  return index;                 /* failsafe */
}

int ColorExtFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int n_ext = 0;
  int a;
  int ok = true;
  int ll;
  register CColor *I = G->Color;
  PyObject *rec;
  ExtRec *ext;

  if(partial_restore) {
    ext = I->Ext;
    for(a = 0; a < I->NExt; a++) {
      ext->old_session_index = 0;
      ext++;
    }
    I->HaveOldSessionExtColors = true;
  } else {
    I->HaveOldSessionExtColors = false;
  }

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);

  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  if(ok) {
    n_ext = PyList_Size(list);
    if(partial_restore) {
      VLACheck(I->Ext, ExtRec, n_ext + I->NExt);
      ext = I->Ext + I->NExt;
    } else {
      VLACheck(I->Ext, ExtRec, n_ext);
      ext = I->Ext;
    }
    for(a = 0; a < n_ext; a++) {
      rec = PyList_GetItem(list, a);
      if(ok)
        ok = (rec != NULL);
      if(ok)
        ok = PyList_Check(rec);
      if(ok) {
        WordType name;
        OVreturn_word result;
        ok = PConvPyStrToStr(PyList_GetItem(rec, 0), name, sizeof(WordType));
        if(OVreturn_IS_OK(result = OVLexicon_GetFromCString(I->Lex, name))) {
          OVOneToOne_Set(I->Idx, result.word, a);
          ext->Name = result.word;
        } else {
          ext->Name = 0;
        }
      }
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(rec, 1), &ext->Type);
      ext->old_session_index = cColorExtCutoff - a;
      ext++;
    }
    if(ok)
      I->NExt = (ext - I->Ext);

  }
  return (ok);
#endif
}


/*========================================================================*/
int ColorFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int n_custom = 0;
  int a;
  int index = 0, old_session_index = 0;
  int ok = true;
  int ll = 0;
  register CColor *I = G->Color;
  PyObject *rec;
  ColorRec *color = NULL;

  if(partial_restore) {
    color = I->Color;
    for(a = 0; a < I->NColor; a++) {
      color->old_session_index = 0;
      color++;
    }
  }
  I->HaveOldSessionColors = false;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok) {
    n_custom = PyList_Size(list);
    for(a = 0; a < n_custom; a++) {
      rec = PyList_GetItem(list, a);
      if(ok)
        ok = (rec != NULL);
      if(ok)
        ok = PyList_Check(rec);
      if(ok)
        ll = PyList_Size(rec);
      /* TO SUPPORT BACKWARDS COMPATIBILITY...
         Always check ll when adding new PyList_GetItem's */
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(rec, 1), &index);
      if(ok) {
        old_session_index = index;
        if(partial_restore) {
          if(I->NColor > index) {       /* conflicts with an existing color... */
            I->HaveOldSessionColors = true;
            index = I->NColor;
          }
        }
        if(index >= I->NColor) {
          VLACheck(I->Color, ColorRec, index);  /* auto-zeros */
          I->NColor = index + 1;
        }
        color = I->Color + index;
        color->old_session_index = old_session_index;
        if(ok) {
          WordType name;
          OVreturn_word result;
          ok = PConvPyStrToStr(PyList_GetItem(rec, 0), name, sizeof(WordType));
          if(OVreturn_IS_OK(result = OVLexicon_GetFromCString(I->Lex, name))) {
            OVOneToOne_Set(I->Idx, result.word, index);
            color->Name = result.word;
          } else {
            color->Name = 0;
          }
        }
        if(ok)
          ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(rec, 2), color->Color, 3);
        if(PyList_Size(rec) >= 6) {
          if(ok)
            ok = PConvPyIntToChar(PyList_GetItem(rec, 3), &color->Custom);
          if(ok)
            ok = PConvPyIntToChar(PyList_GetItem(rec, 4), &color->LutColorFlag);
          if(ok)
            ok =
              PConvPyListToFloatArrayInPlace(PyList_GetItem(rec, 5), color->LutColor, 3);
        } else {
          if(ok) {
            color->Custom = true;
          }
        }
      }
      if(ok && (ll > 6)) {
        if(ok)
          ok = PConvPyIntToChar(PyList_GetItem(rec, 6), &color->Fixed);
      } else if(ok && color) {
        color->Fixed = false;
      }
      if(!ok)
        break;
    }
  }
  return (ok);
#endif
}


/*========================================================================*/
#if 0
int ColorTableFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  register CColor *I = G->Color;
  int ok = true;
  int ll = 0;
  if(!partial_restore) {
    if(ok)
      ok = (list != NULL);
    if(ok)
      ok = PyList_Check(list);
    if(ok) {
      ll = PyList_Size(list);
      if(ll > 1) {
        float tmp_float;
        if(ok)
          ok = PConvPyFloatToFloat(PyList_GetItem(list, 0), &tmp_float);
        if(ok)
          I->Gamma = tmp_float;
        {
          PyObject *tmp = PyList_GetItem(list, 1);
          if(tmp && (tmp != Py_None)) {
#if 0
            int PConvPyListToIntArrayInPlace(PyObject * obj, int *ff, ov_size ll);

            int PConvPyListToIntArrayInPlace(tmp, int *ff, ov_size ll);
            PyList_SetItem(list, 2, PConvIntArrayToPyList(I->ColorTable, 512 * 512));

#endif

            ColorUpdateFromLut(G, -1);
          }
        }
      }
    }
  }
  return ok;
#endif
}
#endif


/*========================================================================*/
void ColorDef(PyMOLGlobals * G, char *name, float *v, int mode, int quiet)
{
  register CColor *I = G->Color;
  int color = -1;
  int a;
  int best;
  int wm;

  {
    OVreturn_word result;
    if(OVreturn_IS_OK(result = OVLexicon_BorrowFromCString(I->Lex, name)))
      if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->Idx, result.word))) {
        color = result.word;
      }
  }

  if(color < 0) {
    best = 0;
    for(a = 0; a < I->NColor; a++) {
      int color_ext = I->Color[a].Name;
      if(color_ext) {
        char *color_name = OVLexicon_FetchCString(I->Lex, color_ext);
        wm = WordMatch(G, name, color_name, true);
        if(wm < 0) {
          color = a;
          break;
        }
      }
    }
  }

  if(color < 0) {
    OVreturn_word result;
    color = I->NColor;
    VLACheck(I->Color, ColorRec, I->NColor);
    I->NColor++;

    if(OVreturn_IS_OK(result = OVLexicon_GetFromCString(I->Lex, name))) {
      OVOneToOne_Set(I->Idx, result.word, color);
      I->Color[color].Name = result.word;
    } else {
      I->Color[color].Name = 0;
    }
  }

  I->Color[color].Color[0] = v[0];
  I->Color[color].Color[1] = v[1];
  I->Color[color].Color[2] = v[2];

  switch (mode) {
  case 1:
    I->Color[color].Fixed = true;
    break;
  default:
    I->Color[color].Fixed = false;
    break;
  }

  I->Color[color].Custom = true;
  ColorUpdateFromLut(G, color);

  if(!quiet) {
    PRINTFB(G, FB_Executive, FB_Actions)
      " Color: \"%s\" defined as [ %3.3f, %3.3f, %3.3f ].\n", name, v[0], v[1], v[2]
      ENDFB(G);

  }

  PRINTFD(G, FB_Color)
    " Color: and assigned number %d.\n", color ENDFD;
}


/*========================================================================*/
int ColorGetIndex(PyMOLGlobals * G, char *name)
{
  register CColor *I = G->Color;
  int color = -1;               /* default for unknown is white */
  int ext_color;
  int a;
  int i;
  int wm, best = 0;
  int ext_best = 0;
  int is_numeric = true;
  int found = false;

  {
    char *c;
    c = name;
    while(*c) {
      if((((*c) < '0') || ((*c) > '9')) && ((*c) != '-')) {
        is_numeric = false;
        break;
      }
      c++;
    }
  }

  if(is_numeric) {
    if(sscanf(name, "%d", &i)) {
      if((i < I->NColor) && (i >= 0))
        return (i);
      else if(i == cColorNewAuto)
        return (ColorGetNext(G));
      else if(i == cColorCurAuto)
        return (ColorGetCurrent(G));
      else if(i == cColorAtomic)
        return cColorAtomic;
      else if(i == cColorObject)
        return cColorObject;
      else if(i == cColorFront)
        return cColorFront;
      else if(i == cColorBack)
        return cColorBack;
      else if(i == -1)
        return -1;
    }
  }
  if((name[0] == '0') && (name[1] == 'x')) {    /* explicit hex RGB 0x000000 */
    int tmp_color;
    if(sscanf(name + 2, "%x", (unsigned int *) &tmp_color) == 1) {
      tmp_color = (cColor_TRGB_Bits |
                   (tmp_color & 0x00FFFFFF) | ((tmp_color >> 2) & 0x3F000000));
      return tmp_color;
    }
  }
  if(WordMatch(G, name, "default", true))
    return (-1);
  if(WordMatch(G, name, "auto", true))
    return (ColorGetNext(G));
  if(WordMatch(G, name, "current", true))
    return (ColorGetCurrent(G));
  if(WordMatch(G, name, "atomic", true))
    return (cColorAtomic);
  if(WordMatch(G, name, "object", true))
    return (cColorObject);
  if(WordMatch(G, name, "front", true))
    return (cColorFront);
  if(WordMatch(G, name, "back", true))
    return (cColorBack);

  if(I->Lex) {                  /* search for a perfect match (fast!) */
    OVreturn_word result;
    if(OVreturn_IS_OK(result = OVLexicon_BorrowFromCString(I->Lex, name)))
      if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->Idx, result.word))) {
        found = true;
        color = result.word;
      }
  }
  if(!found) {                  /* search for an imperfect match */
    for(a = 0; a < I->NColor; a++) {
      int color_ext = I->Color[a].Name;
      if(color_ext) {
        char *color_name = OVLexicon_FetchCString(I->Lex, color_ext);
        wm = WordMatch(G, name, color_name, true);
        if(wm < 0) {
          color = a;
          best = 0;
          break;
        } else if((wm > 0) && (best < wm)) {
          color = a;
          best = wm;
        }
      }
    }
    if(best || (color < 0)) {
      ext_color = ColorFindExtByName(G, name, false, &ext_best);
      if(ext_color >= 0) {
        ext_color = -10 - ext_color;    /* indicates external */
        if((!ext_best) || (ext_best > best))    /* perfect or better match? */
          color = ext_color;
      }
    }
  }
  return (color);
}


/*========================================================================*/
float *ColorGetNamed(PyMOLGlobals * G, char *name)
{
  return (ColorGet(G, ColorGetIndex(G, name)));
}


/*========================================================================*/
char *ColorGetName(PyMOLGlobals * G, int index)
{
  register CColor *I = G->Color;
  if((index >= 0) && (index < I->NColor)) {
    return OVLexicon_FetchCString(I->Lex, I->Color[index].Name);
  } else if((index & cColor_TRGB_Mask) == cColor_TRGB_Bits) {
    index = (((index & 0xFFFFFF) | ((index << 2) & 0xFC000000) |        /* convert 6 bits of trans into 8 */
              ((index >> 4) & 0x03000000)));
    if(index & 0xFF000000)      /* if transparent */
      sprintf(I->RGBName, "%08x", index);
    else                        /* else */
      sprintf(I->RGBName, "%06x", index);
    return I->RGBName;
  } else if(index <= cColorExtCutoff) {
    int a = cColorExtCutoff - index;
    if(a < I->NExt) {
      return OVLexicon_FetchCString(I->Lex, I->Ext[a].Name);
    } else
      return NULL;
  }
  return (NULL);
}


/*========================================================================*/
int ColorGetStatus(PyMOLGlobals * G, int index)
{
  register CColor *I = G->Color;
  /* return 0 if color is invalid, -1 if hidden; 
     1 otherwise */
  int result = 0;
  if((index >= 0) && (index < I->NColor)) {
    int color_ext = I->Color[index].Name;
    if(color_ext) {
      char *c = OVLexicon_FetchCString(I->Lex, color_ext);
      result = 1;
      while(*c) {
        if(((*c) >= '0') && ((*c) <= '9')) {
          result = -1;
          break;
        }
        c++;
      }
    }
  }
  return (result);
}


/*========================================================================*/
int ColorGetNColor(PyMOLGlobals * G)
{
  register CColor *I = G->Color;
  return (I->NColor);
}


/*========================================================================*/
void ColorFree(PyMOLGlobals * G)
{
  register CColor *I = G->Color;
  if(I->ColorTable) {
    FreeP(I->ColorTable);
  }
  VLAFreeP(I->Color);
  VLAFreeP(I->Ext);
  if(I->Lex)
    OVLexicon_Del(I->Lex);
  if(I->Idx)
    OVOneToOne_Del(I->Idx);
  FreeP(I);
}


/*========================================================================*/

static int reg_name(OVLexicon * lex, OVOneToOne * o2o, int index, char *name)
{
  OVreturn_word result;
  if(OVreturn_IS_OK(result = OVLexicon_GetFromCString(lex, name))) {
    OVOneToOne_Set(o2o, result.word, index);
    return result.word;
  } else {
    return 0;
  }
}

void ColorReset(PyMOLGlobals * G)
{

/* PyMOL core color names

  1   1   1   white
 .5  .5  .5   grey/gray
  0   0   0   black 

  1   0   0   red
  0   1   0   green
  0   0   1   blue

  1   1   0   yellow
  1   0   1   magenta
  0   1   1   cyan

  1   1  .5   paleyellow  .
  1  .5   1   violet      .
 .5   1   1   aquamarine  .

  1  .5  .5   deepsalmon  .
 .5   1  .5   palegreen   .
 .5  .5   1   slate       .

 .75 .75  0   olive       .
 .75  0  .75  purple      .
  0  .75 .75  teal        .

 .6  .6  .1   deepolive   .
 .6  .1  .6   deeppurple  .
 .1  .6  .6   deepteal    .

  1  .5   0   orange      .
  1   0  .5   hotpink     .
 .5   1   0   chartreuse  .
  0   1  .5   limegreen   .
  0  .5   1   marine      .
 .5   0   1   purpleblue  .

*/

  register CColor *I = G->Color;
  register ColorRec *color = I->Color;
  register int n_color = 0;

  char name[10];
  int a;
  int set1;
  float f;
  float spectrumS[13][3] = {
    {1.0, 0.0, 1.0},            /* magenta - 0 */
    {0.5, 0.0, 1.0},
    {0.0, 0.0, 1.0},            /* blue - 166.66  */
    {0.0, 0.5, 1.0},
    {0.0, 1.0, 1.0},            /* cyan - 333.33 */

    {0.0, 1.0, 0.5},
    {0.0, 1.0, 0.0},            /* green - 500 */
    {0.5, 1.0, 0.0},
    {1.0, 1.0, 0.0},            /* yellow - 666.66 */
    {1.0, 0.5, 0.0},

    {1.0, 0.0, 0.0},            /* red - 833.33 */
    {1.0, 0.0, 0.5},
    {1.0, 0.0, 1.0},            /* magenta - 999 */
  };

  float spectrumR[13][3] = {
    {1.0, 1.0, 0.0},            /* yellow - 0 */
    {0.5, 1.0, 0.0},            /* chartreuse */
    {0.0, 1.0, 0.0},            /* green - 166.66 */
    {0.0, 1.0, 0.5},            /* limegreen */
    {0.0, 1.0, 1.0},            /* cyan - 333.33 */

    {0.0, 0.5, 1.0},            /* marine */
    {0.0, 0.0, 1.0},            /* blue - 500 */
    {0.5, 0.0, 1.0},            /* purpleblue */
    {1.0, 0.0, 1.0},            /* magenta - 666.66 */
    {1.0, 0.0, 0.5},            /* hotpink */

    {1.0, 0.0, 0.0},            /* red - 833.33 */
    {1.0, 0.5, 0.0},            /* orange */
    {1.0, 1.0, 0.0},            /* yellow - 999 */
  };

  float spectrumC[][3] = {
    {1.0, 1.0, 0.0},            /* yellow - 0 */
    {0.0, 0.0, 1.0},            /* blue - 83.333 */
    {1.0, 0.0, 0.0},            /* red - 167.67 */
    {0.0, 1.0, 0.0},            /* green - 250.00 */
    {1.0, 0.0, 1.0},            /* magenta - 333.33 */

    {0.0, 1.0, 1.0},            /* cyan - 416.67 */
    {1.0, 1.0, 0.0},            /* yellow - 500.00 */
    {0.0, 1.0, 0.0},            /* green - 583.33 */
    {0.0, 0.0, 1.0},            /* blue - 666.67 */
    {1.0, 0.0, 1.0},            /* magenta - 750.00 */

    {1.0, 1.0, 0.0},            /* yellow - 833.33 */
    {1.0, 0.0, 0.0},            /* red - 916.67 */
    {0.0, 1.0, 1.0},            /* cyan - 999 */
  };

  float spectrumW[][3] = {
    {1.0, 1.0, 0.0},            /* yellow - 0 */
    {1.0, 1.0, 1.0},            /* white */
    {0.0, 0.0, 1.0},            /* blue  - 83.333 */
    {1.0, 1.0, 1.0},            /* white */
    {1.0, 0.0, 0.0},            /* red - 166.67 */

    {1.0, 1.0, 1.0},            /* white */
    {0.0, 1.0, 0.0},            /* green - 250.00 */
    {1.0, 1.0, 1.0},            /* white */
    {1.0, 0.0, 1.0},            /* magenta - 333.33 */
    {1.0, 1.0, 1.0},            /* white */

    {0.0, 1.0, 1.0},            /* cyan - 416.67 */
    {1.0, 1.0, 1.0},            /* white */
    {1.0, 1.0, 0.0},            /* yellow - 500.00 */
    {1.0, 1.0, 1.0},            /* white */
    {0.0, 1.0, 0.0},            /* green - 583.33 */

    {1.0, 1.0, 1.0},            /* white */
    {0.0, 0.0, 1.0},            /* blue - 666.67 */
    {1.0, 1.0, 1.0},            /* white */
    {1.0, 0.0, 1.0},            /* magenta - 750.00 */
    {1.0, 1.0, 1.0},            /* white */

    {1.0, 1.0, 0.0},            /* yellow - 833.33 */
    {1.0, 1.0, 1.0},            /* white */
    {1.0, 0.0, 0.0},            /* red - 916.67 */
    {1.0, 1.0, 1.0},            /* white */
    {0.0, 1.0, 1.0},            /* cyan - 999 */
  };

  float spectrumO[29][3] = {
    /* a rainbow with perceptive color balancing and extra blue/red at the ends */
    {1.0, 0.0, 1.0},            /* violet */
    {0.8F, 0.0, 1.0},

    {0.5F, 0.0, 1.0},           /* blend */

    {0.0, 0.0, 1.0},            /* blue */
    {0.0, 0.0, 1.0},            /* blue */
    {0.0, 0.2F, 1.0},

    {0.0, 0.5F, 1.0},           /* blend */

    {0.0, 0.8F, 1.0},
    {0.0, 1.0, 1.0},            /* cyan */
    {0.0, 1.0, 0.8F},

    {0.0, 1.0, 0.5F},           /* blend */

    {0.0, 1.0, 0.2F},
    {0.0, 1.0, 0.0},            /* green */
    {0.2F, 1.0, 0.0},

    {0.5F, 1.0, 0.0},           /* blend */

    {0.8F, 1.0, 0.0},
    {1.0, 1.0, 0.0},            /* yellow */
    {1.0, 0.9F, 0.0},

    {1.0, 0.75F, 0.0},          /* blend */

    {1.0, 0.6F, 0.0},
    {1.0, 0.5F, 0.0},           /* orange */
    {1.0, 0.4F, 0.0},

    {1.0, 0.3F, 0.0},           /* blend */

    {1.0, 0.2F, 0.0},
    {1.0, 0.0, 0.0},            /* red */
    {1.0, 0.0, 0.0},            /* red */

    {1.0, 0.0, 0.5F},           /* blend */

    {1.0, 0.0, 0.8F},           /* violet */
    {1.0, 0.0, 1.0},            /* violet */
  };

  if(I->Lex)
    OVLexicon_Del(I->Lex);
  I->Lex = OVLexicon_New(G->Context->heap);

  if(I->Idx)
    OVOneToOne_Del(I->Idx);

  I->Idx = OVOneToOne_New(G->Context->heap);

  /* BLUE->VIOLET->RED r546 to r909 */
  /* BLUE->CYAN->GREEN->YELLOW->RED s182 to s909 */
  /* BLUE->WHITE->RED w00 to */

  color->Name = reg_name(I->Lex, I->Idx, n_color, "white");
  color->Color[0] = 1.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "black");
  color->Color[0] = 0.0F;
  color->Color[1] = 0.0F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "blue");
  color->Color[0] = 0.0F;
  color->Color[1] = 0.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "green");
  color->Color[0] = 0.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "red");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.0F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "cyan");
  color->Color[0] = 0.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "yellow");
  color->Color[0] = 1.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "dash");
  color->Color[0] = 1.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "magenta");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "salmon");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.6F;       /* was 0.5 */
  color->Color[2] = 0.6F;       /* wat 0.5 */
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lime");
  color->Color[0] = 0.5F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "slate");
  color->Color[0] = 0.5F;
  color->Color[1] = 0.5F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "hotpink");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.0F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "orange");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.5F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "chartreuse");        /* AKA puke green */
  color->Color[0] = 0.5F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "limegreen");
  color->Color[0] = 0.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "purpleblue");        /* legacy name */
  color->Color[0] = 0.5F;
  color->Color[1] = 0.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "marine");
  color->Color[0] = 0.0F;
  color->Color[1] = 0.5F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "olive");
  color->Color[0] = 0.77F;
  color->Color[1] = 0.70F;
  color->Color[2] = 0.00F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "purple");
  color->Color[0] = 0.75F;
  color->Color[1] = 0.00F;
  color->Color[2] = 0.75F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "teal");
  color->Color[0] = 0.00F;
  color->Color[1] = 0.75F;
  color->Color[2] = 0.75F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "ruby");
  color->Color[0] = 0.6F;
  color->Color[1] = 0.2F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "forest");
  color->Color[0] = 0.2F;
  color->Color[1] = 0.6F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deepblue");  /* was "deep" */
  color->Color[0] = 0.25F;
  color->Color[1] = 0.25F;
  color->Color[2] = 0.65F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "grey");      /* english spelling */
  color->Color[0] = 0.5F;
  color->Color[1] = 0.5F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "gray");      /* american spelling */
  color->Color[0] = 0.5F;
  color->Color[1] = 0.5F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "carbon");
  color->Color[0] = 0.2F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "nitrogen");
  color->Color[0] = 0.2F;
  color->Color[1] = 0.2F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "oxygen");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.3F;
  color->Color[2] = 0.3F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "hydrogen");
  color->Color[0] = 0.9F;
  color->Color[1] = 0.9F;
  color->Color[2] = 0.9F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "brightorange");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.7F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "sulfur");
  color->Color[0] = 0.9F;       /* needs to be far enough from "yellow" */
  color->Color[1] = 0.775F;     /* to be resolved, while still slightly on */
  color->Color[2] = 0.25F;      /* the yellow side of yelloworange */
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tv_red");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.2F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tv_green");
  color->Color[0] = 0.2F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tv_blue");
  color->Color[0] = 0.3F;
  color->Color[1] = 0.3F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tv_yellow");
  color->Color[0] = 1.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "yelloworange");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.87F;
  color->Color[2] = 0.37F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tv_orange");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.55F;
  color->Color[2] = 0.15F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br0");
  color->Color[0] = 0.1F;
  color->Color[1] = 0.1F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br1");
  color->Color[0] = 0.2F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.9F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br2");
  color->Color[0] = 0.3F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.8F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br3");
  color->Color[0] = 0.4F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.7F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br4");
  color->Color[0] = 0.5F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.6F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br5");
  color->Color[0] = 0.6F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br6");
  color->Color[0] = 0.7F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.4F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br7");
  color->Color[0] = 0.8F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.3F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br8");
  color->Color[0] = 0.9F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.2F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "br9");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.1F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "pink");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.65F;
  color->Color[2] = 0.85F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "firebrick");
  color->Color[0] = 0.698F;
  color->Color[1] = 0.13F;
  color->Color[2] = 0.13F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "chocolate");
  color->Color[0] = 0.555F;
  color->Color[1] = 0.222F;
  color->Color[2] = 0.111F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "brown");
  color->Color[0] = 0.65F;
  color->Color[1] = 0.32F;
  color->Color[2] = 0.17F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "wheat");
  color->Color[0] = 0.99F;
  color->Color[1] = 0.82F;
  color->Color[2] = 0.65F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "violet");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.5F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  /* greybow */

  strcpy(name, "grey00");       /* english spelling */
  for(a = 0; a < 100; a = a + 1) {
    name[5] = (a % 10) + '0';
    name[4] = ((a % 100) / 10) + '0';
    /* sprintf(color->Name,"grey%02d",a); */
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = a / 99.0F;
    color->Color[1] = a / 99.0F;
    color->Color[2] = a / 99.0F;
    n_color++;
    color++;
  }

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lightmagenta");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.2F;
  color->Color[2] = 0.8F;
  n_color++;
  color++;

#define A_DIV 83.333333333F

  /* full spectrum (s000-s999) */

  strcpy(name, "s000");
  for(a = 0; a < 1000; a = a + 1) {
    set1 = (int) (a / A_DIV);
    name[3] = (a % 10) + '0';
    name[2] = ((a % 100) / 10) + '0';
    name[1] = ((a % 1000) / 100) + '0';
    /* sprintf(color->Name,"s%03d",a); */
    f = 1.0F - (a - (set1 * A_DIV)) / A_DIV;
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = f * spectrumS[set1][0] + (1.0F - f) * spectrumS[set1 + 1][0];
    color->Color[1] = f * spectrumS[set1][1] + (1.0F - f) * spectrumS[set1 + 1][1];
    color->Color[2] = f * spectrumS[set1][2] + (1.0F - f) * spectrumS[set1 + 1][2];
    n_color++;
    color++;
  }

  /* offset & reversed full spectrum (r000-r999) */

  strcpy(name, "r000");
  for(a = 0; a < 1000; a = a + 1) {
    set1 = (int) (a / A_DIV);
    /* sprintf(color->Name,"r%03d",a); */
    name[3] = (a % 10) + '0';
    name[2] = ((a % 100) / 10) + '0';
    name[1] = ((a % 1000) / 100) + '0';
    f = 1.0F - (a - (set1 * A_DIV)) / A_DIV;
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = f * spectrumR[set1][0] + (1.0F - f) * spectrumR[set1 + 1][0];
    color->Color[1] = f * spectrumR[set1][1] + (1.0F - f) * spectrumR[set1 + 1][1];
    color->Color[2] = f * spectrumR[set1][2] + (1.0F - f) * spectrumR[set1 + 1][2];
    n_color++;
    color++;
  }

  /* complementary spectra (c000-c999) */

  strcpy(name, "c000");
  for(a = 0; a < 1000; a = a + 1) {
    set1 = (int) (a / A_DIV);
    /*     sprintf(color->Name,"c%03d",a); */
    name[3] = (a % 10) + '0';
    name[2] = ((a % 100) / 10) + '0';
    name[1] = ((a % 1000) / 100) + '0';
    f = 1.0F - (a - (set1 * A_DIV)) / A_DIV;
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = f * spectrumC[set1][0] + (1.0F - f) * spectrumC[set1 + 1][0];
    color->Color[1] = f * spectrumC[set1][1] + (1.0F - f) * spectrumC[set1 + 1][1];
    color->Color[2] = f * spectrumC[set1][2] + (1.0F - f) * spectrumC[set1 + 1][2];
    n_color++;
    color++;
  }

#define W_DIV 41.666666667F

  /* complementary spectra separated by white (w000-w999) */

  strcpy(name, "w000");
  for(a = 0; a < 1000; a = a + 1) {
    set1 = (int) (a / W_DIV);
    /* sprintf(color->Name,"w%03d",a); */
    name[3] = (a % 10) + '0';
    name[2] = ((a % 100) / 10) + '0';
    name[1] = ((a % 1000) / 100) + '0';
    f = 1.0F - (a - (set1 * W_DIV)) / W_DIV;
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = f * spectrumW[set1][0] + (1.0F - f) * spectrumW[set1 + 1][0];
    color->Color[1] = f * spectrumW[set1][1] + (1.0F - f) * spectrumW[set1 + 1][1];
    color->Color[2] = f * spectrumW[set1][2] + (1.0F - f) * spectrumW[set1 + 1][2];
    n_color++;
    color++;
  }

  color->Name = reg_name(I->Lex, I->Idx, n_color, "density");
  color->Color[0] = 0.1F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.6F;
  n_color++;
  color++;

  strcpy(name, "gray00");       /* american */
  for(a = 0; a < 100; a = a + 1) {
    name[5] = (a % 10) + '0';
    name[4] = ((a % 100) / 10) + '0';
    /* sprintf(color->Name,"gray%02d",a); */
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = a / 99.0F;
    color->Color[1] = a / 99.0F;
    color->Color[2] = a / 99.0F;
    n_color++;
    color++;
  }

  /* original full spectrum, with extra blue and red at the ends (o000-o999) */

#define B_DIV 35.7143F

  strcpy(name, "o000");
  for(a = 0; a < 1000; a = a + 1) {
    set1 = (int) (a / B_DIV);
    name[3] = (a % 10) + '0';
    name[2] = ((a % 100) / 10) + '0';
    name[1] = ((a % 1000) / 100) + '0';
    /* sprintf(color->Name,"o%03d",a); */
    f = 1.0F - (a - (set1 * B_DIV)) / B_DIV;
    color->Name = reg_name(I->Lex, I->Idx, n_color, name);
    color->Color[0] = f * spectrumO[set1][0] + (1.0F - f) * spectrumO[set1 + 1][0];
    color->Color[1] = f * spectrumO[set1][1] + (1.0F - f) * spectrumO[set1 + 1][1];
    color->Color[2] = f * spectrumO[set1][2] + (1.0F - f) * spectrumO[set1 + 1][2];
    n_color++;
    color++;
  }

  color->Name = reg_name(I->Lex, I->Idx, n_color, "paleyellow");
  color->Color[0] = 1.0F;
  color->Color[1] = 1.0F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "aquamarine");
  color->Color[0] = 0.5F;
  color->Color[1] = 1.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deepsalmon");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.5F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "palegreen");
  color->Color[0] = 0.65F;
  color->Color[1] = 0.9F;
  color->Color[2] = 0.65F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deepolive");
  color->Color[0] = 0.6F;
  color->Color[1] = 0.6F;
  color->Color[2] = 0.1F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deeppurple");
  color->Color[0] = 0.6F;
  color->Color[1] = 0.1F;
  color->Color[2] = 0.6F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deepteal");
  color->Color[0] = 0.1F;
  color->Color[1] = 0.6F;
  color->Color[2] = 0.6F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lightblue");
  color->Color[0] = 0.75F;
  color->Color[1] = 0.75;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lightorange");
  color->Color[0] = 1.0F;
  color->Color[1] = 0.8F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "palecyan");
  color->Color[0] = 0.8F;
  color->Color[1] = 1.0F;
  color->Color[2] = 1.0F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lightteal");
  color->Color[0] = 0.4F;
  color->Color[1] = 0.7F;
  color->Color[2] = 0.7F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "splitpea");
  color->Color[0] = 0.52F;
  color->Color[1] = 0.75F;
  color->Color[2] = 0.00F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "raspberry");
  color->Color[0] = 0.70F;
  color->Color[1] = 0.30F;
  color->Color[2] = 0.40F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "sand");
  color->Color[0] = 0.72F;
  color->Color[1] = 0.55F;
  color->Color[2] = 0.30F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "smudge");
  color->Color[0] = 0.55F;
  color->Color[1] = 0.70F;
  color->Color[2] = 0.40F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "violetpurple");
  color->Color[0] = 0.55F;
  color->Color[1] = 0.25F;
  color->Color[2] = 0.60F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "dirtyviolet");
  color->Color[0] = 0.70F;
  color->Color[1] = 0.50F;
  color->Color[2] = 0.50F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deepsalmon");
  color->Color[0] = 1.00F;
  color->Color[1] = 0.42F;
  color->Color[2] = 0.42F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lightpink");
  color->Color[0] = 1.00F;
  color->Color[1] = 0.75F;
  color->Color[2] = 0.87F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "greencyan");
  color->Color[0] = 0.25F;
  color->Color[1] = 1.00F;
  color->Color[2] = 0.75F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "limon");
  color->Color[0] = 0.75F;
  color->Color[1] = 1.00F;
  color->Color[2] = 0.25F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "skyblue");
  color->Color[0] = 0.20F;
  color->Color[1] = 0.50F;
  color->Color[2] = 0.80F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "bluewhite");
  color->Color[0] = 0.85F;
  color->Color[1] = 0.85F;
  color->Color[2] = 1.00F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "warmpink");
  color->Color[0] = 0.85F;
  color->Color[1] = 0.20F;
  color->Color[2] = 0.50F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "darksalmon");
  color->Color[0] = 0.73F;
  color->Color[1] = 0.55F;
  color->Color[2] = 0.52F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "helium");
  color->Color[0] = 0.850980392F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lithium");
  color->Color[0] = 0.800000000F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "beryllium");
  color->Color[0] = 0.760784314F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "boron");
  color->Color[0] = 1.000000000F;
  color->Color[1] = 0.709803922F;
  color->Color[2] = 0.709803922F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "fluorine");
  color->Color[0] = 0.701960784F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "neon");
  color->Color[0] = 0.701960784F;
  color->Color[1] = 0.890196078F;
  color->Color[2] = 0.960784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "sodium");
  color->Color[0] = 0.670588235F;
  color->Color[1] = 0.360784314F;
  color->Color[2] = 0.949019608F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "magnesium");
  color->Color[0] = 0.541176471F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "aluminum");
  color->Color[0] = 0.749019608F;
  color->Color[1] = 0.650980392F;
  color->Color[2] = 0.650980392F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "silicon");
  color->Color[0] = 0.941176471F;
  color->Color[1] = 0.784313725F;
  color->Color[2] = 0.627450980F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "phosphorus");
  color->Color[0] = 1.000000000F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "chlorine");
  color->Color[0] = 0.121568627F;
  color->Color[1] = 0.941176471F;
  color->Color[2] = 0.121568627F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "argon");
  color->Color[0] = 0.501960784F;
  color->Color[1] = 0.819607843F;
  color->Color[2] = 0.890196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "potassium");
  color->Color[0] = 0.560784314F;
  color->Color[1] = 0.250980392F;
  color->Color[2] = 0.831372549F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "calcium");
  color->Color[0] = 0.239215686F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "scandium");
  color->Color[0] = 0.901960784F;
  color->Color[1] = 0.901960784F;
  color->Color[2] = 0.901960784F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "titanium");
  color->Color[0] = 0.749019608F;
  color->Color[1] = 0.760784314F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "vanadium");
  color->Color[0] = 0.650980392F;
  color->Color[1] = 0.650980392F;
  color->Color[2] = 0.670588235F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "chromium");
  color->Color[0] = 0.541176471F;
  color->Color[1] = 0.600000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "manganese");
  color->Color[0] = 0.611764706F;
  color->Color[1] = 0.478431373F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "iron");
  color->Color[0] = 0.878431373F;
  color->Color[1] = 0.400000000F;
  color->Color[2] = 0.200000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "cobalt");
  color->Color[0] = 0.941176471F;
  color->Color[1] = 0.564705882F;
  color->Color[2] = 0.627450980F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "nickel");
  color->Color[0] = 0.313725490F;
  color->Color[1] = 0.815686275F;
  color->Color[2] = 0.313725490F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "copper");
  color->Color[0] = 0.784313725F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 0.200000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "zinc");
  color->Color[0] = 0.490196078F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 0.690196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "gallium");
  color->Color[0] = 0.760784314F;
  color->Color[1] = 0.560784314F;
  color->Color[2] = 0.560784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "germanium");
  color->Color[0] = 0.400000000F;
  color->Color[1] = 0.560784314F;
  color->Color[2] = 0.560784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "arsenic");
  color->Color[0] = 0.741176471F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 0.890196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "selenium");
  color->Color[0] = 1.000000000F;
  color->Color[1] = 0.631372549F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "bromine");
  color->Color[0] = 0.650980392F;
  color->Color[1] = 0.160784314F;
  color->Color[2] = 0.160784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "krypton");
  color->Color[0] = 0.360784314F;
  color->Color[1] = 0.721568627F;
  color->Color[2] = 0.819607843F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "rubidium");
  color->Color[0] = 0.439215686F;
  color->Color[1] = 0.180392157F;
  color->Color[2] = 0.690196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "strontium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "yttrium");
  color->Color[0] = 0.580392157F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "zirconium");
  color->Color[0] = 0.580392157F;
  color->Color[1] = 0.878431373F;
  color->Color[2] = 0.878431373F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "niobium");
  color->Color[0] = 0.450980392F;
  color->Color[1] = 0.760784314F;
  color->Color[2] = 0.788235294F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "molybdenum");
  color->Color[0] = 0.329411765F;
  color->Color[1] = 0.709803922F;
  color->Color[2] = 0.709803922F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "technetium");
  color->Color[0] = 0.231372549F;
  color->Color[1] = 0.619607843F;
  color->Color[2] = 0.619607843F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "ruthenium");
  color->Color[0] = 0.141176471F;
  color->Color[1] = 0.560784314F;
  color->Color[2] = 0.560784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "rhodium");
  color->Color[0] = 0.039215686F;
  color->Color[1] = 0.490196078F;
  color->Color[2] = 0.549019608F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "palladium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.411764706F;
  color->Color[2] = 0.521568627F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "silver");
  color->Color[0] = 0.752941176F;
  color->Color[1] = 0.752941176F;
  color->Color[2] = 0.752941176F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "cadmium");
  color->Color[0] = 1.000000000F;
  color->Color[1] = 0.850980392F;
  color->Color[2] = 0.560784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "indium");
  color->Color[0] = 0.650980392F;
  color->Color[1] = 0.458823529F;
  color->Color[2] = 0.450980392F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tin");
  color->Color[0] = 0.400000000F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 0.501960784F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "antimony");
  color->Color[0] = 0.619607843F;
  color->Color[1] = 0.388235294F;
  color->Color[2] = 0.709803922F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tellurium");
  color->Color[0] = 0.831372549F;
  color->Color[1] = 0.478431373F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "iodine");
  color->Color[0] = 0.580392157F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.580392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "xenon");
  color->Color[0] = 0.258823529F;
  color->Color[1] = 0.619607843F;
  color->Color[2] = 0.690196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "cesium");
  color->Color[0] = 0.341176471F;
  color->Color[1] = 0.090196078F;
  color->Color[2] = 0.560784314F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "barium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.788235294F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lanthanum");
  color->Color[0] = 0.439215686F;
  color->Color[1] = 0.831372549F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "cerium");
  color->Color[0] = 1.000000000F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "praseodymium");
  color->Color[0] = 0.850980392F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "neodymium");
  color->Color[0] = 0.780392157F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "promethium");
  color->Color[0] = 0.639215686F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "samarium");
  color->Color[0] = 0.560784314F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "europium");
  color->Color[0] = 0.380392157F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "gadolinium");
  color->Color[0] = 0.270588235F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "terbium");
  color->Color[0] = 0.188235294F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "dysprosium");
  color->Color[0] = 0.121568627F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.780392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "holmium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 1.000000000F;
  color->Color[2] = 0.611764706F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "erbium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.901960784F;
  color->Color[2] = 0.458823529F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "thulium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.831372549F;
  color->Color[2] = 0.321568627F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "ytterbium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.749019608F;
  color->Color[2] = 0.219607843F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lutetium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.670588235F;
  color->Color[2] = 0.141176471F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "hafnium");
  color->Color[0] = 0.301960784F;
  color->Color[1] = 0.760784314F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tantalum");
  color->Color[0] = 0.301960784F;
  color->Color[1] = 0.650980392F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "tungsten");
  color->Color[0] = 0.129411765F;
  color->Color[1] = 0.580392157F;
  color->Color[2] = 0.839215686F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "rhenium");
  color->Color[0] = 0.149019608F;
  color->Color[1] = 0.490196078F;
  color->Color[2] = 0.670588235F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "osmium");
  color->Color[0] = 0.149019608F;
  color->Color[1] = 0.400000000F;
  color->Color[2] = 0.588235294F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "iridium");
  color->Color[0] = 0.090196078F;
  color->Color[1] = 0.329411765F;
  color->Color[2] = 0.529411765F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "platinum");
  color->Color[0] = 0.815686275F;
  color->Color[1] = 0.815686275F;
  color->Color[2] = 0.878431373F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "gold");
  color->Color[0] = 1.000000000F;
  color->Color[1] = 0.819607843F;
  color->Color[2] = 0.137254902F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "mercury");
  color->Color[0] = 0.721568627F;
  color->Color[1] = 0.721568627F;
  color->Color[2] = 0.815686275F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "thallium");
  color->Color[0] = 0.650980392F;
  color->Color[1] = 0.329411765F;
  color->Color[2] = 0.301960784F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lead");
  color->Color[0] = 0.341176471F;
  color->Color[1] = 0.349019608F;
  color->Color[2] = 0.380392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "bismuth");
  color->Color[0] = 0.619607843F;
  color->Color[1] = 0.309803922F;
  color->Color[2] = 0.709803922F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "polonium");
  color->Color[0] = 0.670588235F;
  color->Color[1] = 0.360784314F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "astatine");
  color->Color[0] = 0.458823529F;
  color->Color[1] = 0.309803922F;
  color->Color[2] = 0.270588235F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "radon");
  color->Color[0] = 0.258823529F;
  color->Color[1] = 0.509803922F;
  color->Color[2] = 0.588235294F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "francium");
  color->Color[0] = 0.258823529F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.400000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "radium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.490196078F;
  color->Color[2] = 0.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "actinium");
  color->Color[0] = 0.439215686F;
  color->Color[1] = 0.670588235F;
  color->Color[2] = 0.980392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "thorium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.729411765F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "protactinium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.631372549F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "uranium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.560784314F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "neptunium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.501960784F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "plutonium");
  color->Color[0] = 0.000000000F;
  color->Color[1] = 0.419607843F;
  color->Color[2] = 1.000000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "americium");
  color->Color[0] = 0.329411765F;
  color->Color[1] = 0.360784314F;
  color->Color[2] = 0.949019608F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "curium");
  color->Color[0] = 0.470588235F;
  color->Color[1] = 0.360784314F;
  color->Color[2] = 0.890196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "berkelium");
  color->Color[0] = 0.541176471F;
  color->Color[1] = 0.309803922F;
  color->Color[2] = 0.890196078F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "californium");
  color->Color[0] = 0.631372549F;
  color->Color[1] = 0.211764706F;
  color->Color[2] = 0.831372549F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "einsteinium");
  color->Color[0] = 0.701960784F;
  color->Color[1] = 0.121568627F;
  color->Color[2] = 0.831372549F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "fermium");
  color->Color[0] = 0.701960784F;
  color->Color[1] = 0.121568627F;
  color->Color[2] = 0.729411765F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "mendelevium");
  color->Color[0] = 0.701960784F;
  color->Color[1] = 0.050980392F;
  color->Color[2] = 0.650980392F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "nobelium");
  color->Color[0] = 0.741176471F;
  color->Color[1] = 0.050980392F;
  color->Color[2] = 0.529411765F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lawrencium");
  color->Color[0] = 0.780392157F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.400000000F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "rutherfordium");
  color->Color[0] = 0.800000000F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.349019608F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "dubnium");
  color->Color[0] = 0.819607843F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.309803922F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "seaborgium");
  color->Color[0] = 0.850980392F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.270588235F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "bohrium");
  color->Color[0] = 0.878431373F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.219607843F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "hassium");
  color->Color[0] = 0.901960784F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.180392157F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "meitnerium");
  color->Color[0] = 0.921568627F;
  color->Color[1] = 0.000000000F;
  color->Color[2] = 0.149019608F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "deuterium");
  color->Color[0] = 0.9F;
  color->Color[1] = 0.9F;
  color->Color[2] = 0.9F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "lonepair");
  color->Color[0] = 0.5F;
  color->Color[1] = 0.5F;
  color->Color[2] = 0.5F;
  n_color++;
  color++;

  color->Name = reg_name(I->Lex, I->Idx, n_color, "pseudoatom");
  color->Color[0] = 0.9F;
  color->Color[1] = 0.9F;
  color->Color[2] = 0.9F;
  n_color++;
  color++;

  color = I->Color;
  for(a = 0; a < n_color; a++) {
    /* mark all current colors non-custom so that they don't get saved in session files */
    color[a].Custom = false;
  }

  I->NColor = n_color;
  I->NExt = 0;
}

int ColorTableLoad(PyMOLGlobals * G, char *fname, float gamma, int quiet)
{
  register CColor *I = G->Color;
  int ok = true;

  I->Gamma = gamma;
  if(!fname[0]) {
    ColorUpdateFromLut(G, -1);
  } else {
    int width = 512, height = 512;
    unsigned int *table = NULL;

    if(!strcmp(fname, "rgb")) {
      if(I->ColorTable) {
        FreeP(I->ColorTable);
        I->ColorTable = NULL;
        PRINTFB(G, FB_Color, FB_Actions)
          " Color: purged table; restoring RGB colors.\n" ENDFB(G);
      }
      ColorUpdateFromLut(G, -1);
    } else if(!strcmp(fname, "greyscale")) {

      int x, y;
      unsigned int r = 0, g = 0, b = 0;
      unsigned int *pixel, mask, *p;
      unsigned int rc;

      FreeP(I->ColorTable);
      if(I->BigEndian)
        mask = 0x000000FF;
      else
        mask = 0xFF000000;

      table = Alloc(unsigned int, 512 * 512);

      p = (unsigned int *) table;
      for(x = 0; x < width; x++)
        for(y = 0; y < height; y++)
          *(p++) = mask;

      for(y = 0; y < height; y++)
        for(x = 0; x < width; x++) {
          rc = (r + g + b)/3;

          pixel = table + ((width) * y) + x;
          if(I->BigEndian) {
            *(pixel) = mask | (rc << 24) | (rc << 16) | (rc << 8);
          } else {
            *(pixel) = mask | (rc << 16) | (rc << 8) | rc;
          }
          b = b + 4;
          if(!(0xFF & b)) {
            b = 0;
            g = g + 4;
            if(!(0xFF & g)) {
              g = 0;
              r = r + 4;
            }
          }
        }

      I->ColorTable = table;
      if(!quiet) {
        PRINTFB(G, FB_Color, FB_Actions)
          " Color: defined table '%s'.\n", fname ENDFB(G);
      }

      ColorUpdateFromLut(G, -1);
      ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
      SceneChanged(G);


    } else if(!strcmp(fname, "pymol")) {

      int x, y;
      unsigned int r = 0, g = 0, b = 0;
      unsigned int *pixel, mask, *p;
      unsigned int rc, bc, gc;
      unsigned int gf, bf, rf;

      float green_max = 0.75F;
      float red_max = 0.95F;
      float blue_max = 0.97F;

      float min_factor = 0.15F;

      red_max = SettingGet(G, cSetting_pymol_space_max_red);
      green_max = SettingGet(G, cSetting_pymol_space_max_green);
      blue_max = SettingGet(G, cSetting_pymol_space_max_blue);
      min_factor = SettingGet(G, cSetting_pymol_space_min_factor);

      FreeP(I->ColorTable);
      if(I->BigEndian)
        mask = 0x000000FF;
      else
        mask = 0xFF000000;

      table = Alloc(unsigned int, 512 * 512);

      p = (unsigned int *) table;
      for(x = 0; x < width; x++)
        for(y = 0; y < height; y++)
          *(p++) = mask;

      for(y = 0; y < height; y++)
        for(x = 0; x < width; x++) {
          rc = r;
          gc = g;
          bc = b;

          if((r >= g) && (r >= b)) {
            if(rc > 255 * red_max) {
              rc = (unsigned int) (red_max * 255);
              bc = bc * rc / r;
              gc = gc * rc / r;
            }
          } else if((g >= b) && (g >= r)) {
            if(gc > 255 * green_max) {
              gc = (unsigned int) (green_max * 255);
              bc = bc * gc / g;
              rc = rc * gc / g;
            }
          } else if((b >= g) && (b >= r)) {
            if(bc > 255 * blue_max) {
              bc = (unsigned int) (blue_max * 255);
              gc = gc * bc / b;
              rc = rc * bc / b;
            }
          }

          rf = (int) (min_factor * rc + 0.49999F);
          gf = (int) (min_factor * gc + 0.49999F);
          bf = (int) (min_factor * bc + 0.49999F);

          if(rc < gf)
            rc = (int) gf;
          if(bc < gf)
            bc = (int) gf;

          if(rc < bf)
            rc = (int) bf;
          if(gc < bf)
            gc = (int) bf;

          if(gc < rf)
            gc = (int) rf;
          if(bc < rf)
            bc = (int) rf;

          if(rc > 255)
            rc = 255;
          if(bc > 255)
            bc = 255;
          if(gc > 255)
            gc = 255;

          pixel = table + ((width) * y) + x;
          if(I->BigEndian) {
            *(pixel) = mask | (rc << 24) | (gc << 16) | (bc << 8);
          } else {
            *(pixel) = mask | (bc << 16) | (gc << 8) | rc;
          }
          b = b + 4;
          if(!(0xFF & b)) {
            b = 0;
            g = g + 4;
            if(!(0xFF & g)) {
              g = 0;
              r = r + 4;
            }
          }
        }

      I->ColorTable = table;
      if(!quiet) {
        PRINTFB(G, FB_Color, FB_Actions)
          " Color: defined table '%s'.\n", fname ENDFB(G);
      }

      ColorUpdateFromLut(G, -1);
      ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
      SceneChanged(G);

    } else {
      if(strlen(fname)) {

        unsigned int u_width = (unsigned int) width, u_height = (unsigned int) height;
        unsigned char *u_table = (unsigned char *) table;
        if(MyPNGRead(fname, &u_table, &u_width, &u_height)) {
          table = (unsigned int *) u_table;
          width = (signed int) u_width;
          height = (signed int) u_height;
          if((width == 512) && (height == 512)) {
            FreeP(I->ColorTable);
            I->ColorTable = table;
            if(!quiet) {
              PRINTFB(G, FB_Color, FB_Actions)
                " Color: loaded table '%s'.\n", fname ENDFB(G);
            }

            ColorUpdateFromLut(G, -1);

          } else {
            PRINTFB(G, FB_Color, FB_Errors)
              " ColorTableLoad-Error: invalid dimensions w x h  = %d x %d; should be 512 x 512.\n",
              width, height ENDFB(G);

            ok = false;
          }
        } else {
          PRINTFB(G, FB_Color, FB_Errors)
            " ColorTableLoad-Error: unable to load '%s'.\n", fname ENDFB(G);
          ok = false;
        }
      } else {
        PRINTFB(G, FB_Color, FB_Actions)
          " Color: purged table; colors unchanged.\n" ENDFB(G);
        FreeP(I->ColorTable);
      }
    }
    if(!ok) {
      FreeP(table);
    }
  }
  if(ok) {
    ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
    SceneChanged(G);
  }
  return (ok);
}

#if 0
static void unlookup_color(unsigned int *table, float *in, float *out, int big_endian)
{
  /* simple iterative approach to finding closest input for a desired
     output given the current color mapping table */

  float cur_in[3], cur_out[3], diff[3];
  copy3f(in, cur_in);
  for(a = 0; a < 10; a++) {
    lookup_color(table, cur_in, cur_out, big_endian);
    diff3f(cur_out, in, diff);
    scale3f(0.5, diff, diff);
    subtract3f(cur_in, diff, cur_in);
  }
  copy3f(cur_in, out);
}
#endif

static void lookup_color(CColor * I, float *in, float *out, int big_endian)
{
  const float _1 = 1.0F;
  unsigned int *table = I->ColorTable;
  if(table) {
    unsigned int r, g, b, rr, gr, br;
    unsigned int ra, ga, ba;
    unsigned int rc[2][2][2], gc[2][2][2], bc[2][2][2];
    unsigned int *entry;
    int x, y, z;
    float fr, fg, fb, frm1x, fgm1, fbm1, rct, gct, bct;
    const float _2 = 2.0F, _0 = 0.0F, _05 = 0.5F, _04999 = 0.4999F;
    const float inv255 = 1.0F / 255.0F;

    r = ((int) (255 * in[0] + _05)) & 0xFF;
    g = ((int) (255 * in[1] + _05)) & 0xFF;
    b = ((int) (255 * in[2] + _05)) & 0xFF;

    rr = r & 0x3;
    gr = g & 0x3;
    br = b & 0x3;

    r = (r >> 2);
    g = (g >> 2);
    b = (b >> 2);

    /* now for a crude little trilinear */

    for(x = 0; x < 2; x++) {
      ra = r + x;
      if(ra > 63)
        ra = 63;
      for(y = 0; y < 2; y++) {
        ga = g + y;
        if(ga > 63)
          ga = 63;
        for(z = 0; z < 2; z++) {
          ba = b + z;
          if(ba > 63)
            ba = 63;

          entry = table + (ra << 12) + (ga << 6) + ba;

          if(big_endian) {
            rc[x][y][z] = 0xFF & ((*entry) >> 24);
            gc[x][y][z] = 0xFF & ((*entry) >> 16);
            bc[x][y][z] = 0xFF & ((*entry) >> 8);
          } else {
            rc[x][y][z] = 0xFF & ((*entry));
            gc[x][y][z] = 0xFF & ((*entry) >> 8);
            bc[x][y][z] = 0xFF & ((*entry) >> 16);
          }
        }
      }
    }

    frm1x = rr / 4.0F;
    fgm1 = gr / 4.0F;
    fbm1 = br / 4.0F;

    fr = 1.0F - frm1x;
    fg = 1.0F - fgm1;
    fb = 1.0F - fbm1;

    rct = _04999 +
      (fr * fg * fb * rc[0][0][0]) +
      (frm1x * fg * fb * rc[1][0][0]) +
      (fr * fgm1 * fb * rc[0][1][0]) +
      (fr * fg * fbm1 * rc[0][0][1]) +
      (frm1x * fgm1 * fb * rc[1][1][0]) +
      (fr * fgm1 * fbm1 * rc[0][1][1]) +
      (frm1x * fg * fbm1 * rc[1][0][1]) + (frm1x * fgm1 * fbm1 * rc[1][1][1]);

    gct = _04999 +
      (fr * fg * fb * gc[0][0][0]) +
      (frm1x * fg * fb * gc[1][0][0]) +
      (fr * fgm1 * fb * gc[0][1][0]) +
      (fr * fg * fbm1 * gc[0][0][1]) +
      (frm1x * fgm1 * fb * gc[1][1][0]) +
      (fr * fgm1 * fbm1 * gc[0][1][1]) +
      (frm1x * fg * fbm1 * gc[1][0][1]) + (frm1x * fgm1 * fbm1 * gc[1][1][1]);

    bct = _04999 +
      (fr * fg * fb * bc[0][0][0]) +
      (frm1x * fg * fb * bc[1][0][0]) +
      (fr * fgm1 * fb * bc[0][1][0]) +
      (fr * fg * fbm1 * bc[0][0][1]) +
      (frm1x * fgm1 * fb * bc[1][1][0]) +
      (fr * fgm1 * fbm1 * bc[0][1][1]) +
      (frm1x * fg * fbm1 * bc[1][0][1]) + (frm1x * fgm1 * fbm1 * bc[1][1][1]);

    if(r >= 63)
      rct += rr;
    if(g >= 63)
      gct += gr;
    if(b >= 63)
      bct += br;

    if(rct <= _2)
      rct = _0;                 /* make sure black is black */
    if(gct <= _2)
      gct = _0;
    if(bct <= _2)
      bct = _0;

    out[0] = rct * inv255;
    out[1] = gct * inv255;
    out[2] = bct * inv255;
  } else {
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
  }

  if((I->Gamma != 1.0F) && (I->Gamma > R_SMALL4)) {
    float inv_gamma = 1.0F / I->Gamma;
    float inp = (out[0] + out[1] + out[2]) * (1 / 3.0F);
    if(inp >= R_SMALL4) {
      float sig = (float) (pow(inp, inv_gamma)) / inp;
      out[0] *= sig;
      out[1] *= sig;
      out[2] *= sig;
    }
  }

  if(out[0] > _1)
    out[0] = _1;
  if(out[1] > _1)
    out[1] = _1;
  if(out[2] > _1)
    out[2] = _1;

}


/*========================================================================*/
void ColorUpdateFromLut(PyMOLGlobals * G, int index)
{
  int i;
  int once = false;
  register CColor *I = G->Color;
  float *color, *new_color;

  I->LUTActive = (I->ColorTable || (I->Gamma != 1.0F));

  i = index;
  if(index >= 0) {
    once = true;
  }
  for(i = 0; i < I->NColor; i++) {
    if(!once)
      index = i;

    if(index < I->NColor) {
      if(!I->LUTActive) {
        I->Color[index].LutColorFlag = false;
      } else if(!I->Color[index].Fixed) {
        color = I->Color[index].Color;
        new_color = I->Color[index].LutColor;
        lookup_color(I, color, new_color, I->BigEndian);

        PRINTFD(G, FB_Color)
          "%5.3f %5.3f %5.3f -> %5.3f %5.3f %5.3f\n",
          color[0], color[1], color[2], new_color[0], new_color[1], new_color[2]
          ENDFD;

        I->Color[index].LutColorFlag = true;
      }
    }

    if(once)
      break;
  }
}


/*========================================================================*/
int ColorLookupColor(PyMOLGlobals * G, float *color)
{
  register CColor *I = G->Color;
  if(I->LUTActive) {
    lookup_color(I, color, color, I->BigEndian);
    return true;
  } else {
    return false;
  }
}


/*========================================================================*/
int ColorInit(PyMOLGlobals * G)
{
  CColor *I = NULL;

  if((I = (G->Color = Calloc(CColor, 1)))) {
    unsigned int test;
    unsigned char *testPtr;

    test = 0xFF000000;
    testPtr = (unsigned char *) &test;
    I->BigEndian = (*testPtr) && 1;

    I->Color = VLAMalloc(5500, sizeof(ColorRec), 5, true);
    I->Ext = VLAMalloc(2, sizeof(ExtRec), 5, true);
    I->Gamma = 1.0F;

    ColorReset(G);              /* will alloc I->Idx and I->Lex */
    I->Front[0] = 1.0F;
    I->Front[1] = 1.0F;
    I->Front[2] = 1.0F;
    return 1;
  } else {
    return 0;
  }
}

void ColorUpdateFront(PyMOLGlobals * G, float *back)
{
  register CColor *I = G->Color;
  copy3f(back, I->Back);
  I->Front[0] = 1.0F - back[0];
  I->Front[1] = 1.0F - back[1];
  I->Front[2] = 1.0F - back[2];
  if(diff3f(I->Front, back) < 0.5F)
    zero3f(I->Front);
}


/*========================================================================*/
float *ColorGetSpecial(PyMOLGlobals * G, int index)
{
  if(index >= 0)
    return ColorGet(G, index);
  else {
    register CColor *I = G->Color;
    I->RGBColor[0] = (float) index;
    I->RGBColor[1] = -1.0F;
    I->RGBColor[2] = -1.0F;
    return I->RGBColor;
  }
}

float *ColorGet(PyMOLGlobals * G, int index)
{
  register CColor *I = G->Color;
  float *ptr;
  if((index >= 0) && (index < I->NColor)) {
    if(I->Color[index].LutColorFlag && SettingGetGlobal_b(G, cSetting_clamp_colors))
      ptr = I->Color[index].LutColor;
    else
      ptr = I->Color[index].Color;
    return (ptr);
  } else if((index & cColor_TRGB_Mask) == cColor_TRGB_Bits) {   /* a 24-bit RGB color */
    I->RGBColor[0] = ((index & 0x00FF0000) >> 16) / 255.0F;
    I->RGBColor[1] = ((index & 0x0000FF00) >> 8) / 255.0F;
    I->RGBColor[2] = ((index & 0x000000FF)) / 255.0F;
    if(I->LUTActive)
      lookup_color(I, I->RGBColor, I->RGBColor, I->BigEndian);
    return I->RGBColor;
  } else if(index == cColorFront) {
    return I->Front;
  } else if(index == cColorBack) {
    return I->Back;
  } else {
    /* invalid color id, then simply return white */
    return (I->Color[0].Color);
  }
}

float *ColorGetRaw(PyMOLGlobals * G, int index)
{
  register CColor *I = G->Color;
  float *ptr;
  if((index >= 0) && (index < I->NColor)) {
    ptr = I->Color[index].Color;
    return (ptr);
  } else if((index & cColor_TRGB_Mask) == cColor_TRGB_Bits) {   /* a 24-bit RGB color */
    I->RGBColor[0] = ((index & 0x00FF0000) >> 16) / 255.0F;
    I->RGBColor[1] = ((index & 0x0000FF00) >> 8) / 255.0F;
    I->RGBColor[2] = ((index & 0x000000FF)) / 255.0F;
    return I->RGBColor;
  } else {
    /* invalid color id, then simply return white */
    return (I->Color[0].Color);
  }
}

int ColorGetEncoded(PyMOLGlobals * G, int index, float *color)
{
  register CColor *I = G->Color;
  float *ptr;
  if((index >= 0) && (index < I->NColor)) {
    if(I->Color[index].LutColorFlag && SettingGetGlobal_b(G, cSetting_clamp_colors))
      ptr = I->Color[index].LutColor;
    else
      ptr = I->Color[index].Color;
    copy3f(ptr, color);
  } else if((index & cColor_TRGB_Mask) == cColor_TRGB_Bits) {   /* a 24-bit RGB color */
    float rgb_color[3];
    rgb_color[0] = ((index & 0x00FF0000) >> 16) / 255.0F;
    rgb_color[1] = ((index & 0x0000FF00) >> 8) / 255.0F;
    rgb_color[2] = ((index & 0x000000FF)) / 255.0F;
    if(I->LUTActive)
      lookup_color(I, rgb_color, rgb_color, I->BigEndian);
    copy3f(rgb_color, color);
  } else if(index <= cColorExtCutoff) {
    color[0] = (float) index;
    color[1] = 0.0F;
    color[2] = 0.0F;
  } else if(index == cColorFront) {
    copy3f(I->Front, color);
  } else if(index == cColorBack) {
    copy3f(I->Back, color);
  } else {
    color[0] = 1.0F;
    color[1] = 1.0F;
    color[2] = 1.0F;
    /* invalid color id, then simply return white */
    return 0;
  }
  return 1;
}
