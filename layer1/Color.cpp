
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
static void lookup_color(CColor * I, const float *in, float *out, int big_endian);

void ColorGetBkrdContColor(PyMOLGlobals * G, float *rgb, int invert_flag)
{
  const float *bkrd = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb));

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

unsigned int ColorGet32BitWord(PyMOLGlobals * G, const float *rgba)
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
  next = SettingGetGlobal_i(G, cSetting_auto_color_next);

  if(next >= nAutoColor)
    next = 0;
  result = AutoColor[next];
  next++;
  if(next >= nAutoColor)
    next = 0;
  SettingSetGlobal_i(G, cSetting_auto_color_next, next);
  return (result);
}

int ColorGetCurrent(PyMOLGlobals * G)
{
  int result;
  int next;
  next = SettingGetGlobal_i(G, cSetting_auto_color_next);
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
  CColor *I = G->Color;
  if(index <= cColorExtCutoff) {
    index = cColorExtCutoff - index;
    if (index < I->Ext.size()) {
      auto& ext = I->Ext[index];
      if (!ext.Ptr && ext.Name) {
        ext.Ptr = ExecutiveFindObject<ObjectGadgetRamp>(G, ext.Name);
      }
      return ext.Ptr;
    }
  }
  return nullptr;
}

int ColorGetRamped(PyMOLGlobals * G, int index, const float *vertex, float *color, int state)
{
  CColor *I = G->Color;
  int ok = false;
  if (auto* ptr = ColorGetRamp(G, index)) {
    ok = ObjectGadgetRampInterVertex(ptr, vertex, color, state);
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

/**
 * Gets a color as 3 floats from an index and writes it into
 * the color argument.  If the index is a ramp, then it uses the vertex and state arguments to lookup the
 * color value in the ramp.
 * NOTES: does not support index values cColorObject(-5) or cColorAtomic(-4) color since the object
 *        or atom color is not passed in.
 *
 * @param index - color index value
 * @param vertex - x/y/z used for ramp lookup (if color index is a ramp)
 * @param[out] color - output color array of 3 floats
 * @param state - state lookup if ramp
 *
 * @return whether the color index is dependent on a ramp.
 */
bool ColorGetCheckRamped(PyMOLGlobals * G, int index, const float *vertex, float *color, int state)
{
  bool isRamped = false;
  if(ColorCheckRamped(G, index)) {
    ColorGetRamped(G, index, vertex, color, state);
    isRamped = true;
  } else {
    copy3f(ColorGet(G, index), color);
  }
  return isRamped;
}

/**
 * Find a record by case-insensitive name
 *
 * @param seq Indexable container (Color or Ext)
 * @param name Color name
 * @return seq index or -1 if not found
 */
template <typename Sequence>
static int findByCaseInsensitiveName(
    PyMOLGlobals* G, const Sequence& seq, const char* name)
{
  for (int a = 0; a < seq.size(); ++a) {
    auto* color_name = seq[a].Name;
    if (color_name) {
      int wm = WordMatch(G, name, color_name, true);
      if (wm < 0) {
        return a;
      }
    }
  }

  return -1;
}

/**
 * Find a record by case-insensitive and/or partial name
 *
 * @param seq Indexable container (Color or Ext)
 * @param name Color name
 * @param[in,out] best Word match score (0 for perfect match), must not be
 * negative
 * @return seq index or -1 if not found
 */
template <typename Sequence>
static int findByCaseInsensitivePrefix(
    PyMOLGlobals* G, const Sequence& seq, const char* name, int& best)
{
  int best_a = -1;
  assert(best >= 0);

  // search for an imperfect match
  for (int a = 0; a < seq.size(); ++a) {
    auto* color_name = seq[a].Name;
    if (color_name) {
      auto wm = WordMatch(G, name, color_name, true);
      if (wm < 0) {
        // perfect case-insensitive match
        best = 0;
        return a;
      }

      if (best < wm) {
        // prefix match
        best = wm;
        best_a = a;
      }
    }
  }

  return best_a;
}

/**
 * Find a color ramp by case-insensitive name
 *
 * @param name Color name (ramp name)
 * @return Ext index or -1 if not found
 */
static int ColorFindExtByName(PyMOLGlobals* G, const char* name)
{
  return findByCaseInsensitiveName(G, G->Color->Ext, name);
}

/**
 * Map name to index (idx[name] = index)
 *
 * If the name is already in use and the index can't be reused, then clear the
 * name on the existing color or ext record.
 *
 * @param index Color index
 * @param name Color name
 * @param reuse If the name already exists, reuse the existing index if possible
 * (not possible to reuse a ramp index for a color or vice versa)
 *
 * @return pointer to stored name string
 */
static const char* reg_name(CColor* const I, CColor::ColorIdx const index,
    const char* name, bool reuse = false)
{
  auto handle = I->Idx.emplace(name, index);
  auto& handle_name = handle.first->first;
  auto& handle_index = handle.first->second;

  if (handle_index != index &&
      (!reuse || bool(cColorExtCutoff < handle_index) !=
                     bool(cColorExtCutoff < index))) {
    assert(!handle.second);

    // if we're stealing a name to a new index, clear the name on the old record
    if (handle_index <= cColorExtCutoff) {
      auto& ext = I->Ext[cColorExtCutoff - handle_index];
      assert(ext.Name == handle_name.c_str());
      ext.Name = nullptr;
    } else if (handle_index >= 0) {
      auto& col = I->Color[handle_index];
      assert(col.Name == handle_name.c_str());
      col.Name = nullptr;
    }

    handle_index = index;
  }

  return handle_name.c_str();
}

void ColorRegisterExt(PyMOLGlobals* G, const char* name, ObjectGadgetRamp* ptr)
{
  CColor *I = G->Color;
  int a;

  a = ColorFindExtByName(G, name);
  if(a < 0) {
    a = I->Ext.size();

    I->Ext.emplace_back();
    auto& ext = I->Ext.back();

    ext.Name = reg_name(I, cColorExtCutoff - a, name);
    assert(I->Idx[ext.Name] == cColorExtCutoff - a);
  }
  if(a >= 0) {
    I->Ext[a].Ptr = ptr;
  }
}

void ColorForgetExt(PyMOLGlobals * G, const char *name)
{
  CColor *I = G->Color;
  auto a = ColorFindExtByName(G, name);

  if (a < 0)
    return;

  // currently leaks memory in I->Ext array
  auto& ext = I->Ext[a];
  ext.Ptr = nullptr;

  // HaveOldSessionExtColors should only be true while we're loading a partial
  // session, and ColorForgetExt probably means that we're replacing the ramp
  // object with a another one, so we don't want to lose the name+index
  // relationship.
  if (ext.Name && !I->HaveOldSessionExtColors) {
    I->Idx.erase(ext.Name);
    ext.Name = nullptr;
  }
}

PyObject *ColorExtAsPyList(PyMOLGlobals * G)
{
  CColor *I = G->Color;

  auto* result = PyList_New(I->Ext.size());

  size_t a = 0;

  for (const auto& ext : I->Ext) {
    auto* list = PyList_New(2);
    const char* name = ext.Name ? ext.Name : "";
    PyList_SetItem(list, 0, PyString_FromString(name));

    // obsolete since PyMOL 2.5, store for backwards compatibility
    PyList_SetItem(list, 1, PyInt_FromLong(cColorGadgetRamp));

    PyList_SetItem(result, a++, list);
  }

  assert(a == I->Ext.size());

  return result;
}


/*========================================================================*/
PyObject *ColorAsPyList(PyMOLGlobals * G)
{
  CColor *I = G->Color;

  size_t n_custom = 0;
  for (const auto& color : I->Color) {
    if (color.Custom || color.LutColorFlag) {
      n_custom++;
    }
  }

  auto* result = PyList_New(n_custom);

  size_t a = 0;
  size_t c = 0;

  for (const auto& color : I->Color) {
    if (color.Custom || color.LutColorFlag) {
      auto* list = PyList_New(7);
      PyList_SetItem(list, 0, PyString_FromString(color.Name));
      PyList_SetItem(list, 1, PyInt_FromLong(a));
      PyList_SetItem(list, 2, PConvFloatArrayToPyList(color.Color, 3));
      PyList_SetItem(list, 3, PyInt_FromLong(color.Custom));
      PyList_SetItem(list, 4, PyInt_FromLong(color.LutColorFlag));
      PyList_SetItem(list, 5, PConvFloatArrayToPyList(color.LutColor, 3));
      PyList_SetItem(list, 6, PyInt_FromLong(color.Fixed));
      PyList_SetItem(result, c++, list);
    }
    ++a;
  }

  assert(c == n_custom);

  return result;
}

/*========================================================================*/
int ColorConvertOldSessionIndex(PyMOLGlobals * G, int index)
{
  CColor *I = G->Color;
  if(index > cColorExtCutoff) {
    if(I->HaveOldSessionColors) {
      for (int a = int(I->Color.size()) - 1; a >= 0; --a) {
        if (index == I->Color[a].old_session_index) {
          return a;
        }
      }
    }
  } else if(I->HaveOldSessionExtColors) {
    for (int a = int(I->Ext.size()) - 1; a >= 0; --a) {
      if (index == I->Ext[a].old_session_index) {
        return cColorExtCutoff - a;
      }
    }
  }
  return index;                 /* failsafe */
}

#define return_error_if_fail(e) p_return_val_if_fail((e), false);

int ColorExtFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore)
{
  CColor *I = G->Color;
  size_t n_ext = 0;

  assert(!I->HaveOldSessionExtColors);

  if (list && PyList_Check(list)) {
    n_ext = PyList_Size(list);
  }

  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  if (partial_restore) {
    I->HaveOldSessionExtColors = I->Ext.size() > 0;

    for (auto& ext : I->Ext) {
      ext.old_session_index = 0;
    }
  } else {
    I->Ext.clear();
  }

  for (int a = 0; a < n_ext; ++a) {
    auto* rec = PyList_GetItem(list, a);

    return_error_if_fail(rec != nullptr);
    return_error_if_fail(PyList_Check(rec));

    std::string name;
    return_error_if_fail(PConvFromPyListItem(G, rec, 0, name));

    char const* name_ptr =
        reg_name(I, cColorExtCutoff - I->Ext.size(), name.c_str(), true);
    int const a_new = cColorExtCutoff - I->Idx[name];

    assert(a_new >= 0);
    assert(a_new <= I->Ext.size());
    assert(a_new == a || partial_restore);

    if (a_new == I->Ext.size()) {
      I->Ext.emplace_back();
    } else {
      assert(partial_restore);
    }

    auto& ext = I->Ext[a_new];
    ext.Name = name_ptr;
    ext.old_session_index = cColorExtCutoff - a;

    CPythonVal_Free(rec);
  }

  return true;
}


/*========================================================================*/
int ColorFromPyList(PyMOLGlobals * G, PyObject * list, int partial_restore)
{
  CColor* I = G->Color;

  assert(!I->HaveOldSessionColors);

  if (partial_restore) {
    for (auto& color : I->Color) {
      color.old_session_index = 0;
    }
  }

  return_error_if_fail(list != nullptr );
  return_error_if_fail(PyList_Check(list));

  int const n_custom = PyList_Size(list);

  for (int a = 0; a < n_custom; ++a) {
    auto rec = PyList_GetItem(list, a);

    return_error_if_fail(rec && PyList_Check(rec));

    auto const ll = PyList_Size(rec);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */

    int index = 0;
    return_error_if_fail(PConvFromPyListItem(G, rec, 1, index));

    std::string name;
    return_error_if_fail(PConvFromPyListItem(G, rec, 0, name));

    int const old_session_index = index;
    if (partial_restore && I->Color.size() > index) {
      // conflicts with an existing color
      index = I->Color.size();
      I->HaveOldSessionColors = true;
    }

    if (index >= I->Color.size()) {
      assert(I->Color.size() == index);
      I->Color.emplace_back(reg_name(I, index, name.c_str()));
    }

    auto& color = I->Color[index];
    color.old_session_index = old_session_index;

    assert(name == color.Name);
    assert(index == I->Idx[name]);

    return_error_if_fail(CPythonVal_PConvPyListToFloatArrayInPlace_From_List(
        G, rec, 2, color.Color, 3));

    if (PyList_Size(rec) >= 6) {
      return_error_if_fail(PConvFromPyListItem(G, rec, 3, color.Custom));
      return_error_if_fail(PConvFromPyListItem(G, rec, 4, color.LutColorFlag));
      return_error_if_fail(CPythonVal_PConvPyListToFloatArrayInPlace_From_List(
          G, rec, 5, color.LutColor, 3));
    } else {
      color.Custom = true;
    }

    if (ll > 6) {
      PConvFromPyListItem(G, rec, 6, color.Fixed);
    } else {
      color.Fixed = false;
    }

    CPythonVal_Free(rec);
  }

  return true;
}

/*========================================================================*/
void ColorDef(PyMOLGlobals * G, const char *name, const float *v, int mode, int quiet)
{
  CColor *I = G->Color;
  int color = -1;

  // Search for a perfect case-sensitive match
  {
    auto it = I->Idx.find(name);
    if (it != I->Idx.end()) {
      color = it->second;
    }
  }

  if (color < 0) {
    // Not found or ramp index -> do slow search
    color = findByCaseInsensitiveName(G, I->Color, name);

    if (color < 0) {
      // Not found -> new entry
      color = I->Color.size();
      I->Color.emplace_back(reg_name(I, color, name));
      assert(I->Idx[name] == color);
    }
  }

  copy3f(v, I->Color[color].Color);

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
int ColorGetIndex(PyMOLGlobals * G, const char *name)
{
  CColor *I = G->Color;
  int i;
  int is_numeric = true;

  {
    const char *c;
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
      if((i < I->Color.size()) && (i >= 0))
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
      else if(i == cColorDefault)
        return cColorDefault;
      if (i & cColor_TRGB_Bits)
        return i;
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

  // the following block used to allow prefix matches (before PyMOL 2.5)
  if(WordMatch(G, name, "default", true) < 0)
    return cColorDefault;
  if(WordMatch(G, name, "auto", true) < 0)
    return (ColorGetNext(G));
  if(WordMatch(G, name, "current", true) < 0)
    return (ColorGetCurrent(G));
  if(WordMatch(G, name, "atomic", true) < 0)
    return (cColorAtomic);
  if(WordMatch(G, name, "object", true) < 0)
    return (cColorObject);
  if(WordMatch(G, name, "front", true) < 0)
    return (cColorFront);
  if(WordMatch(G, name, "back", true) < 0)
    return (cColorBack);

  // search for a perfect case-sensitive match (fast!)
  {
    auto it = I->Idx.find(name);
    if (it != I->Idx.end()) {
      return it->second;
    }
  }

  // search for case-insensitive or partial match
  // TODO does this even make sense? What's the use case? Should this be
  // restricted to non-ambiguous matches? Note that the Python cmd.color()
  // function does its own prefix lookup and rejects ambiguous matches.
  int best = 0;
  int color = findByCaseInsensitivePrefix(G, I->Color, name, best);
  if (best != 0 || color < 0) {
    int const ext_color = findByCaseInsensitivePrefix(G, I->Ext, name, best);
    if (ext_color >= 0) {
      color = cColorExtCutoff - ext_color;
    }
  }

  return color;
}


/*========================================================================*/
const float *ColorGetNamed(PyMOLGlobals * G, const char *name)
{
  return (ColorGet(G, ColorGetIndex(G, name)));
}


/*========================================================================*/
const char *ColorGetName(PyMOLGlobals * G, int index)
{
  CColor *I = G->Color;
  if((index >= 0) && (index < I->Color.size())) {
    return I->Color[index].Name;
  } else if((index & cColor_TRGB_Mask) == cColor_TRGB_Bits) {
    index = (((index & 0xFFFFFF) | ((index << 2) & 0xFC000000) |        /* convert 6 bits of trans into 8 */
              ((index >> 4) & 0x03000000)));
    if(index & 0xFF000000)      /* if transparent */
      sprintf(I->RGBName, "0x%08x", index);
    else                        /* else */
      sprintf(I->RGBName, "0x%06x", index);
    return I->RGBName;
  } else if(index <= cColorExtCutoff) {
    int a = cColorExtCutoff - index;
    if (a < I->Ext.size()) {
      return I->Ext[a].Name;
    } else
      return NULL;
  }
  return (NULL);
}


/*========================================================================*/
int ColorGetStatus(PyMOLGlobals * G, int index)
{
  CColor *I = G->Color;
  /* return 0 if color is invalid, -1 if hidden; 
     1 otherwise */
  int result = 0;
  if((index >= 0) && (index < I->Color.size())) {
    auto* color_name = I->Color[index].Name;
    if(color_name) {
      const char* c = color_name;
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
  CColor *I = G->Color;
  return (I->Color.size());
}


/*========================================================================*/
void ColorFree(PyMOLGlobals * G)
{
  DeleteP(G->Color);
}


/*========================================================================*/

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

  CColor *I = G->Color;

  I->Idx.clear();
  I->Ext.clear();

  auto& Color = I->Color;
  Color.clear();
  Color.reserve(5500);

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

  /* BLUE->VIOLET->RED r546 to r909 */
  /* BLUE->CYAN->GREEN->YELLOW->RED s182 to s909 */
  /* BLUE->WHITE->RED w00 to */

#define reg_named_color(name, R, G, B)                                         \
  {                                                                            \
    Color.emplace_back(reg_name(I, Color.size(), name));                       \
    set3f(Color.back().Color, R, G, B);                                        \
    assert(I->Idx[name] == Color.size() - 1);                                  \
  }

  reg_named_color("white", 1.F, 1.F, 1.F);
  reg_named_color("black", 0.F, 0.F, 0.F);
  reg_named_color("blue", 0.F, 0.F, 1.F);
  reg_named_color("green", 0.F, 1.F, 0.F);
  reg_named_color("red", 1.F, 0.F, 0.F);
  reg_named_color("cyan", 0.F, 1.F, 1.F);
  reg_named_color("yellow", 1.F, 1.F, 0.F);
  reg_named_color("dash", 1.F, 1.F, 0.F);
  reg_named_color("magenta", 1.F, 0.F, 1.F);
  reg_named_color("salmon", 1.F, 0.6F, 0.6F);
  reg_named_color("lime", 0.5F, 1.F, 0.5F);
  reg_named_color("slate", 0.5F, 0.5F, 1.F);
  reg_named_color("hotpink", 1.F, 0.F, 0.5F);
  reg_named_color("orange", 1.F, 0.5F, 0.F);
  reg_named_color("chartreuse", 0.5F, 1.F, 0.F); /* AKA puke green */
  reg_named_color("limegreen", 0.F, 1.F, 0.5F);
  reg_named_color("purpleblue", 0.5F, 0.F, 1.F); /* legacy name */
  reg_named_color("marine", 0.F, 0.5F, 1.F);
  reg_named_color("olive", 0.77F, 0.7F, 0.F);
  reg_named_color("purple", 0.75F, 0.F, 0.75F);
  reg_named_color("teal", 0.F, 0.75F, 0.75F);
  reg_named_color("ruby", 0.6F, 0.2F, 0.2F);
  reg_named_color("forest", 0.2F, 0.6F, 0.2F);
  reg_named_color("deepblue", 0.25F, 0.25F, 0.65F); /* was "deep" */
  reg_named_color("grey", 0.5F, 0.5F, 0.5F); /* english spelling */
  reg_named_color("gray", 0.5F, 0.5F, 0.5F); /* american spelling */
  reg_named_color("carbon", 0.2F, 1.F, 0.2F);
  reg_named_color("nitrogen", 0.2F, 0.2F, 1.F);
  reg_named_color("oxygen", 1.F, 0.3F, 0.3F);
  reg_named_color("hydrogen", 0.9F, 0.9F, 0.9F);
  reg_named_color("brightorange", 1.F, 0.7F, 0.2F);
  reg_named_color("sulfur", 0.9F, 0.775F, 0.25F);
  reg_named_color("tv_red", 1.F, 0.2F, 0.2F);
  reg_named_color("tv_green", 0.2F, 1.F, 0.2F);
  reg_named_color("tv_blue", 0.3F, 0.3F, 1.F);
  reg_named_color("tv_yellow", 1.F, 1.F, 0.2F);
  reg_named_color("yelloworange", 1.F, 0.87F, 0.37F);
  reg_named_color("tv_orange", 1.F, 0.55F, 0.15F);
  reg_named_color("br0", 0.1F, 0.1F, 1.F);
  reg_named_color("br1", 0.2F, 0.1F, 0.9F);
  reg_named_color("br2", 0.3F, 0.1F, 0.8F);
  reg_named_color("br3", 0.4F, 0.1F, 0.7F);
  reg_named_color("br4", 0.5F, 0.1F, 0.6F);
  reg_named_color("br5", 0.6F, 0.1F, 0.5F);
  reg_named_color("br6", 0.7F, 0.1F, 0.4F);
  reg_named_color("br7", 0.8F, 0.1F, 0.3F);
  reg_named_color("br8", 0.9F, 0.1F, 0.2F);
  reg_named_color("br9", 1.F, 0.1F, 0.1F);
  reg_named_color("pink", 1.F, 0.65F, 0.85F);
  reg_named_color("firebrick", 0.698F, 0.13F, 0.13F);
  reg_named_color("chocolate", 0.555F, 0.222F, 0.111F);
  reg_named_color("brown", 0.65F, 0.32F, 0.17F);
  reg_named_color("wheat", 0.99F, 0.82F, 0.65F);
  reg_named_color("violet", 1.F, 0.5F, 1.F);

  /* greybow */

  strcpy(name, "grey00");       /* english spelling */
  for(a = 0; a < 100; a = a + 1) {
    name[5] = (a % 10) + '0';
    name[4] = ((a % 100) / 10) + '0';
    /* sprintf(color->Name,"grey%02d",a); */
    reg_named_color(name, a / 99.F, a / 99.F, a / 99.F);
  }

  reg_named_color("lightmagenta", 1.F, 0.2F, 0.8F);

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
    reg_named_color(name,
        f * spectrumS[set1][0] + (1.F - f) * spectrumS[set1 + 1][0],
        f * spectrumS[set1][1] + (1.F - f) * spectrumS[set1 + 1][1],
        f * spectrumS[set1][2] + (1.F - f) * spectrumS[set1 + 1][2]);
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
    reg_named_color(name,
        f * spectrumR[set1][0] + (1.F - f) * spectrumR[set1 + 1][0],
        f * spectrumR[set1][1] + (1.F - f) * spectrumR[set1 + 1][1],
        f * spectrumR[set1][2] + (1.F - f) * spectrumR[set1 + 1][2]);
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
    reg_named_color(name,
        f * spectrumC[set1][0] + (1.F - f) * spectrumC[set1 + 1][0],
        f * spectrumC[set1][1] + (1.F - f) * spectrumC[set1 + 1][1],
        f * spectrumC[set1][2] + (1.F - f) * spectrumC[set1 + 1][2]);
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
    reg_named_color(name,
        f * spectrumW[set1][0] + (1.F - f) * spectrumW[set1 + 1][0],
        f * spectrumW[set1][1] + (1.F - f) * spectrumW[set1 + 1][1],
        f * spectrumW[set1][2] + (1.F - f) * spectrumW[set1 + 1][2]);
  }

  reg_named_color("density", 0.1F, 0.1F, 0.6F);

  strcpy(name, "gray00");       /* american */
  for(a = 0; a < 100; a = a + 1) {
    name[5] = (a % 10) + '0';
    name[4] = ((a % 100) / 10) + '0';
    /* sprintf(color->Name,"gray%02d",a); */
    reg_named_color(name, a / 99.F, a / 99.F, a / 99.F);
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
    reg_named_color(name,
        f * spectrumO[set1][0] + (1.F - f) * spectrumO[set1 + 1][0],
        f * spectrumO[set1][1] + (1.F - f) * spectrumO[set1 + 1][1],
        f * spectrumO[set1][2] + (1.F - f) * spectrumO[set1 + 1][2]);
  }

  reg_named_color("paleyellow", 1.F, 1.F, 0.5F);
  reg_named_color("aquamarine", 0.5F, 1.F, 1.F);
  reg_named_color("deepsalmon", 1.F, 0.5F, 0.5F);
  reg_named_color("palegreen", 0.65F, 0.9F, 0.65F);
  reg_named_color("deepolive", 0.6F, 0.6F, 0.1F);
  reg_named_color("deeppurple", 0.6F, 0.1F, 0.6F);
  reg_named_color("deepteal", 0.1F, 0.6F, 0.6F);
  reg_named_color("lightblue", 0.75F, 0.75F, 1.F);
  reg_named_color("lightorange", 1.F, 0.8F, 0.5F);
  reg_named_color("palecyan", 0.8F, 1.F, 1.F);
  reg_named_color("lightteal", 0.4F, 0.7F, 0.7F);
  reg_named_color("splitpea", 0.52F, 0.75F, 0.F);
  reg_named_color("raspberry", 0.7F, 0.3F, 0.4F);
  reg_named_color("sand", 0.72F, 0.55F, 0.3F);
  reg_named_color("smudge", 0.55F, 0.7F, 0.4F);
  reg_named_color("violetpurple", 0.55F, 0.25F, 0.6F);
  reg_named_color("dirtyviolet", 0.7F, 0.5F, 0.5F);

  // was deepsalmon (duplicated name!)
  reg_named_color("_deepsalmon", 1.F, 0.42F, 0.42F);

  reg_named_color("lightpink", 1.F, 0.75F, 0.87F);
  reg_named_color("greencyan", 0.25F, 1.F, 0.75F);
  reg_named_color("limon", 0.75F, 1.F, 0.25F);
  reg_named_color("skyblue", 0.2F, 0.5F, 0.8F);
  reg_named_color("bluewhite", 0.85F, 0.85F, 1.F);
  reg_named_color("warmpink", 0.85F, 0.2F, 0.5F);
  reg_named_color("darksalmon", 0.73F, 0.55F, 0.52F);
  reg_named_color("helium", 0.850980392F, 1.F, 1.F);
  reg_named_color("lithium", 0.8F, 0.501960784F, 1.F);
  reg_named_color("beryllium", 0.760784314F, 1.F, 0.F);
  reg_named_color("boron", 1.F, 0.709803922F, 0.709803922F);
  reg_named_color("fluorine", 0.701960784F, 1.F, 1.F);
  reg_named_color("neon", 0.701960784F, 0.890196078F, 0.960784314F);
  reg_named_color("sodium", 0.670588235F, 0.360784314F, 0.949019608F);
  reg_named_color("magnesium", 0.541176471F, 1.F, 0.F);
  reg_named_color("aluminum", 0.749019608F, 0.650980392F, 0.650980392F);
  reg_named_color("silicon", 0.941176471F, 0.784313725F, 0.62745098F);
  reg_named_color("phosphorus", 1.F, 0.501960784F, 0.F);
  reg_named_color("chlorine", 0.121568627F, 0.941176471F, 0.121568627F);
  reg_named_color("argon", 0.501960784F, 0.819607843F, 0.890196078F);
  reg_named_color("potassium", 0.560784314F, 0.250980392F, 0.831372549F);
  reg_named_color("calcium", 0.239215686F, 1.F, 0.F);
  reg_named_color("scandium", 0.901960784F, 0.901960784F, 0.901960784F);
  reg_named_color("titanium", 0.749019608F, 0.760784314F, 0.780392157F);
  reg_named_color("vanadium", 0.650980392F, 0.650980392F, 0.670588235F);
  reg_named_color("chromium", 0.541176471F, 0.6F, 0.780392157F);
  reg_named_color("manganese", 0.611764706F, 0.478431373F, 0.780392157F);
  reg_named_color("iron", 0.878431373F, 0.4F, 0.2F);
  reg_named_color("cobalt", 0.941176471F, 0.564705882F, 0.62745098F);
  reg_named_color("nickel", 0.31372549F, 0.815686275F, 0.31372549F);
  reg_named_color("copper", 0.784313725F, 0.501960784F, 0.2F);
  reg_named_color("zinc", 0.490196078F, 0.501960784F, 0.690196078F);
  reg_named_color("gallium", 0.760784314F, 0.560784314F, 0.560784314F);
  reg_named_color("germanium", 0.4F, 0.560784314F, 0.560784314F);
  reg_named_color("arsenic", 0.741176471F, 0.501960784F, 0.890196078F);
  reg_named_color("selenium", 1.F, 0.631372549F, 0.F);
  reg_named_color("bromine", 0.650980392F, 0.160784314F, 0.160784314F);
  reg_named_color("krypton", 0.360784314F, 0.721568627F, 0.819607843F);
  reg_named_color("rubidium", 0.439215686F, 0.180392157F, 0.690196078F);
  reg_named_color("strontium", 0.F, 1.F, 0.F);
  reg_named_color("yttrium", 0.580392157F, 1.F, 1.F);
  reg_named_color("zirconium", 0.580392157F, 0.878431373F, 0.878431373F);
  reg_named_color("niobium", 0.450980392F, 0.760784314F, 0.788235294F);
  reg_named_color("molybdenum", 0.329411765F, 0.709803922F, 0.709803922F);
  reg_named_color("technetium", 0.231372549F, 0.619607843F, 0.619607843F);
  reg_named_color("ruthenium", 0.141176471F, 0.560784314F, 0.560784314F);
  reg_named_color("rhodium", 0.039215686F, 0.490196078F, 0.549019608F);
  reg_named_color("palladium", 0.F, 0.411764706F, 0.521568627F);
  reg_named_color("silver", 0.752941176F, 0.752941176F, 0.752941176F);
  reg_named_color("cadmium", 1.F, 0.850980392F, 0.560784314F);
  reg_named_color("indium", 0.650980392F, 0.458823529F, 0.450980392F);
  reg_named_color("tin", 0.4F, 0.501960784F, 0.501960784F);
  reg_named_color("antimony", 0.619607843F, 0.388235294F, 0.709803922F);
  reg_named_color("tellurium", 0.831372549F, 0.478431373F, 0.F);
  reg_named_color("iodine", 0.580392157F, 0.F, 0.580392157F);
  reg_named_color("xenon", 0.258823529F, 0.619607843F, 0.690196078F);
  reg_named_color("cesium", 0.341176471F, 0.090196078F, 0.560784314F);
  reg_named_color("barium", 0.F, 0.788235294F, 0.F);
  reg_named_color("lanthanum", 0.439215686F, 0.831372549F, 1.F);
  reg_named_color("cerium", 1.F, 1.F, 0.780392157F);
  reg_named_color("praseodymium", 0.850980392F, 1.F, 0.780392157F);
  reg_named_color("neodymium", 0.780392157F, 1.F, 0.780392157F);
  reg_named_color("promethium", 0.639215686F, 1.F, 0.780392157F);
  reg_named_color("samarium", 0.560784314F, 1.F, 0.780392157F);
  reg_named_color("europium", 0.380392157F, 1.F, 0.780392157F);
  reg_named_color("gadolinium", 0.270588235F, 1.F, 0.780392157F);
  reg_named_color("terbium", 0.188235294F, 1.F, 0.780392157F);
  reg_named_color("dysprosium", 0.121568627F, 1.F, 0.780392157F);
  reg_named_color("holmium", 0.F, 1.F, 0.611764706F);
  reg_named_color("erbium", 0.F, 0.901960784F, 0.458823529F);
  reg_named_color("thulium", 0.F, 0.831372549F, 0.321568627F);
  reg_named_color("ytterbium", 0.F, 0.749019608F, 0.219607843F);
  reg_named_color("lutetium", 0.F, 0.670588235F, 0.141176471F);
  reg_named_color("hafnium", 0.301960784F, 0.760784314F, 1.F);
  reg_named_color("tantalum", 0.301960784F, 0.650980392F, 1.F);
  reg_named_color("tungsten", 0.129411765F, 0.580392157F, 0.839215686F);
  reg_named_color("rhenium", 0.149019608F, 0.490196078F, 0.670588235F);
  reg_named_color("osmium", 0.149019608F, 0.4F, 0.588235294F);
  reg_named_color("iridium", 0.090196078F, 0.329411765F, 0.529411765F);
  reg_named_color("platinum", 0.815686275F, 0.815686275F, 0.878431373F);
  reg_named_color("gold", 1.F, 0.819607843F, 0.137254902F);
  reg_named_color("mercury", 0.721568627F, 0.721568627F, 0.815686275F);
  reg_named_color("thallium", 0.650980392F, 0.329411765F, 0.301960784F);
  reg_named_color("lead", 0.341176471F, 0.349019608F, 0.380392157F);
  reg_named_color("bismuth", 0.619607843F, 0.309803922F, 0.709803922F);
  reg_named_color("polonium", 0.670588235F, 0.360784314F, 0.F);
  reg_named_color("astatine", 0.458823529F, 0.309803922F, 0.270588235F);
  reg_named_color("radon", 0.258823529F, 0.509803922F, 0.588235294F);
  reg_named_color("francium", 0.258823529F, 0.F, 0.4F);
  reg_named_color("radium", 0.F, 0.490196078F, 0.F);
  reg_named_color("actinium", 0.439215686F, 0.670588235F, 0.980392157F);
  reg_named_color("thorium", 0.F, 0.729411765F, 1.F);
  reg_named_color("protactinium", 0.F, 0.631372549F, 1.F);
  reg_named_color("uranium", 0.F, 0.560784314F, 1.F);
  reg_named_color("neptunium", 0.F, 0.501960784F, 1.F);
  reg_named_color("plutonium", 0.F, 0.419607843F, 1.F);
  reg_named_color("americium", 0.329411765F, 0.360784314F, 0.949019608F);
  reg_named_color("curium", 0.470588235F, 0.360784314F, 0.890196078F);
  reg_named_color("berkelium", 0.541176471F, 0.309803922F, 0.890196078F);
  reg_named_color("californium", 0.631372549F, 0.211764706F, 0.831372549F);
  reg_named_color("einsteinium", 0.701960784F, 0.121568627F, 0.831372549F);
  reg_named_color("fermium", 0.701960784F, 0.121568627F, 0.729411765F);
  reg_named_color("mendelevium", 0.701960784F, 0.050980392F, 0.650980392F);
  reg_named_color("nobelium", 0.741176471F, 0.050980392F, 0.529411765F);
  reg_named_color("lawrencium", 0.780392157F, 0.F, 0.4F);
  reg_named_color("rutherfordium", 0.8F, 0.F, 0.349019608F);
  reg_named_color("dubnium", 0.819607843F, 0.F, 0.309803922F);
  reg_named_color("seaborgium", 0.850980392F, 0.F, 0.270588235F);
  reg_named_color("bohrium", 0.878431373F, 0.F, 0.219607843F);
  reg_named_color("hassium", 0.901960784F, 0.F, 0.180392157F);
  reg_named_color("meitnerium", 0.921568627F, 0.F, 0.149019608F);
  reg_named_color("deuterium", 0.9F, 0.9F, 0.9F);
  reg_named_color("lonepair", 0.5F, 0.5F, 0.5F);
  reg_named_color("pseudoatom", 0.9F, 0.9F, 0.9F);
}

int ColorTableLoad(PyMOLGlobals * G, const char *fname, float gamma, int quiet)
{
  CColor *I = G->Color;
  int ok = true;

  I->Gamma = gamma;
  if(!fname[0]) {
    ColorUpdateFromLut(G, -1);
  } else {
    int width = 512, height = 512;
    if(!strcmp(fname, "rgb")) {
      if(!I->ColorTable.empty()) {
        I->ColorTable.clear();
        PRINTFB(G, FB_Color, FB_Actions)
          " Color: purged table; restoring RGB colors.\n" ENDFB(G);
      }
      ColorUpdateFromLut(G, -1);
    } else if(!strcmp(fname, "greyscale")) {

      int x, y;
      unsigned int r = 0, g = 0, b = 0;
      unsigned int *pixel, mask, *p;
      unsigned int rc;

      if(I->BigEndian)
        mask = 0x000000FF;
      else
        mask = 0xFF000000;

      I->ColorTable.resize(512 * 512);

      p = I->ColorTable.data();
      for(x = 0; x < width; x++)
        for(y = 0; y < height; y++)
          *(p++) = mask;

      for(y = 0; y < height; y++)
        for(x = 0; x < width; x++) {
          rc = (r + g + b)/3;

          pixel = I->ColorTable.data() + ((width) * y) + x;
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

      red_max = SettingGetGlobal_f(G, cSetting_pymol_space_max_red);
      green_max = SettingGetGlobal_f(G, cSetting_pymol_space_max_green);
      blue_max = SettingGetGlobal_f(G, cSetting_pymol_space_max_blue);
      min_factor = SettingGetGlobal_f(G, cSetting_pymol_space_min_factor);

      if(I->BigEndian)
        mask = 0x000000FF;
      else
        mask = 0xFF000000;

      I->ColorTable.resize(512 * 512);

      p = I->ColorTable.data();
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

          pixel = I->ColorTable.data() + ((width) * y) + x;
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

      if(!quiet) {
        PRINTFB(G, FB_Color, FB_Actions)
          " Color: defined table '%s'.\n", fname ENDFB(G);
      }

      ColorUpdateFromLut(G, -1);
      ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
      SceneChanged(G);

    } else {
      if(strlen(fname)) {

        auto image = MyPNGRead(fname);
        if(image) {
          std::tie(width, height) = image->getSize();
          if((width == 512) && (height == 512)) {
            auto imageSize = width * height;
            I->ColorTable.resize(imageSize);
            std::copy(image->pixels(), image->pixels() + imageSize, I->ColorTable.data());

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
        I->ColorTable.clear();
      }
    }
  }
  if(ok) {
    ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
    SceneChanged(G);
  }
  return (ok);
}

static void lookup_color(CColor * I, const float *in, float *out, int big_endian)
{
  const float _1 = 1.0F;
  unsigned int *table = I->ColorTable.data();
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
  CColor *I = G->Color;
  float *color, *new_color;

  I->LUTActive = (!I->ColorTable.empty() || (I->Gamma != 1.0F));

  i = index;
  if(index >= 0) {
    once = true;
  }
  for(i = 0; i < I->Color.size(); i++) {
    if(!once)
      index = i;

    if(index < I->Color.size()) {
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
  CColor *I = G->Color;
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

  if ((G->Color = new CColor())) {
    I = G->Color;
    unsigned int test;
    unsigned char *testPtr;

    test = 0xFF000000;
    testPtr = (unsigned char *) &test;
    I->BigEndian = (*testPtr) & 0x01;

    ColorReset(G);              /* will alloc I->Idx and I->Lex */
    return 1;
  } else {
    return 0;
  }
}

void ColorUpdateFront(PyMOLGlobals * G, const float *back)
{
  CColor *I = G->Color;
  copy3f(back, I->Back);
  I->Front[0] = 1.0F - back[0];
  I->Front[1] = 1.0F - back[1];
  I->Front[2] = 1.0F - back[2];
  if(diff3f(I->Front, back) < 0.5F)
    zero3f(I->Front);
}

void ColorUpdateFrontFromSettings(PyMOLGlobals * G){
  int bg_gradient = SettingGet_b(G, NULL, NULL, cSetting_bg_gradient);
  const char * bg_image_filename = SettingGet_s(G, NULL, NULL, cSetting_bg_image_filename);
  short bg_image = bg_image_filename && bg_image_filename[0];
  
  if (!bg_gradient){
    if (!bg_image && !OrthoBackgroundDataIsSet(*G->Ortho)){
      const float *v = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb));
      ColorUpdateFront(G, v);
    } else {
      float v[] = { 0.f, 0.f, 0.f };
      ColorUpdateFront(G, v);
    }
  } else {
    float vv[3];
    const float *v = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb_bottom));
    const float *vb = ColorGet(G, SettingGet_color(G, NULL, NULL, cSetting_bg_rgb_top));
    average3f(v, vb, vv);
    ColorUpdateFront(G, vv);    
  }
}


/*========================================================================*/
const float *ColorGetSpecial(PyMOLGlobals * G, int index)
{
  if(index >= 0)
    return ColorGet(G, index);
  else {
    CColor *I = G->Color;
    I->RGBColor[0] = (float) index;
    I->RGBColor[1] = -1.0F;
    I->RGBColor[2] = -1.0F;
    return I->RGBColor;
  }
}

const float *ColorGet(PyMOLGlobals * G, int index)
{
  CColor *I = G->Color;
  const float *ptr;
  if((index >= 0) && (index < I->Color.size())) {
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

const float *ColorGetRaw(PyMOLGlobals * G, int index)
{
  CColor *I = G->Color;
  const float *ptr;
  if((index >= 0) && (index < I->Color.size())) {
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
  CColor *I = G->Color;
  float *ptr;
  if((index >= 0) && (index < I->Color.size())) {
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

int Color3fToInt(PyMOLGlobals * G, const float *rgb){
  unsigned int rc, gc, bc;
  rc = pymol_roundf(rgb[0] * 255.);
  gc = pymol_roundf(rgb[1] * 255.);
  bc = pymol_roundf(rgb[2] * 255.);
  return ( ( cColor_TRGB_Bits & 0xFF000000) | 
	   ( ( rc << 16 ) & 0x00FF0000) |
	   ( ( gc << 8 ) & 0x0000FF00) |
	   ( ( bc & 0x000000FF ) ) );
}
