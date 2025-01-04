
/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/
#include "os_python.h"

#include "os_gl.h"
#include "os_predef.h"
#include "os_std.h"

#include "Base.h"
#include "CGO.h"
#include "CGOGL.h"
#include "CGORenderer.h"
#include "CoordSet.h"
#include "Err.h"
#include "Feedback.h"
#include "GadgetSet.h"
#include "Matrix.h"
#include "ObjectGadgetRamp.h"
#include "P.h"
#include "PConv.h"
#include "Picking.h"
#include "PrintUtils.h"
#include "PyMOLGlobals.h"
#include "Ray.h"
#include "Rep.h"
#include "Scene.h"
#include "ScenePicking.h"
#include "Setting.h"
#include "ShaderMgr.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Util.h"
#include "VFont.h"
#include "Vector.h"

#include "pymol/algorithm.h"

#define VAR_FOR_NORMAL pl
#define VERTEX_NORMAL_SIZE 3
#define VAR_FOR_NORMAL_CNT_PLUS

#define VALUES_PER_IMPOSTER_SPACE_COORD 1

constexpr unsigned VERTEX_PICKCOLOR_RGBA_SIZE = 1;  // 4 unsigned bytes
constexpr unsigned VERTEX_PICKCOLOR_INDEX_SIZE = 2; // index + bond
constexpr unsigned VERTEX_PICKCOLOR_SIZE = VERTEX_PICKCOLOR_RGBA_SIZE + //
                                           VERTEX_PICKCOLOR_INDEX_SIZE;
constexpr unsigned VERTEX_ACCESSIBILITY_SIZE = 1;

#include <algorithm>
#include <cassert>
#include <iostream>

template <typename T> inline T CLAMPVALUE(T val, T minimum, T maximum)
{
  return (val < minimum) ? minimum : (val > maximum) ? maximum : val;
}

inline cCylCap cap1_from_cyl_shader_bits(const unsigned char bits)
{
  return (bits & 1) ? ((bits & cCylShaderCap1RoundBit) ? cCylCap::Round
                                                       : cCylCap::Flat)
                    : cCylCap::None;
}
inline cCylCap cap2_from_cyl_shader_bits(const unsigned char bits)
{
  return (bits & 2) ? ((bits & cCylShaderCap2RoundBit) ? cCylCap::Round
                                                       : cCylCap::Flat)
                    : cCylCap::None;
}
inline unsigned char cyl_shader_bits_from_caps(
    const cCylCap cap1, const cCylCap cap2)
{
  return ((cap1 == cCylCap::Flat)
                 ? cCylShaderCap1Flat
                 : ((cap1 == cCylCap::Round) ? cCylShaderCap1Round : 0)) |
         ((cap2 == cCylCap::Flat)
                 ? cCylShaderCap2Flat
                 : ((cap2 == cCylCap::Round) ? cCylShaderCap2Round : 0));
}

#define WARN_UNEXPECTED_OPERATION(G, op)                                       \
  PRINTFB(G, FB_CGO, FB_Warnings)                                              \
  " %s-Warning: unexpected op=0x%x (line %d)\n", __func__, op, __LINE__ ENDFB(G)

// like g_return_val_if_fail from glib
#define RETURN_VAL_IF_FAIL(expr, val)                                          \
  {                                                                            \
    if (!(expr))                                                               \
      return (val);                                                            \
  }

static void set_current_pick_color(CGO* cgo, unsigned int idx, int bnd)
{
  if (cgo) {
    cgo->current_pick_color_index = idx;
    cgo->current_pick_color_bond = bnd;
  }
}

/**
 * Assign a new pick color index for {context, index, bond} (only in the
 * first pass) and convert that to an RGBA color for the current picking pass.
 *
 * @param[out] cgo Set the "current pick color" in this CGO
 * @param[in,out] pickmgr Pick color manager
 * @param[out] color RGBA pick color to use in this pass
 * @param[in] context Object state identifier
 * @param[in] index Primary index
 * @param[in] bond Secondary index
 */
void AssignNewPickColor(CGO* cgo, PickColorManager* pickmgr,
    unsigned char* color, const PickContext* context, unsigned int index,
    int bond)
{
  set_current_pick_color(cgo, index, bond);
  pickmgr->colorNext(color, context, index, bond);
}

std::size_t CGO_sz[] = {
    CGO_NULL_SZ,                           //
    CGO_NULL_SZ,                           //
    CGO_BEGIN_SZ,                          //
    CGO_END_SZ,                            //
    CGO_VERTEX_SZ,                         //
    CGO_NORMAL_SZ,                         //
    fsizeof<cgo::draw::color>(),           //
    fsizeof<cgo::draw::sphere>(),          //
    CGO_TRIANGLE_SZ,                       //
    fsizeof<cgo::draw::cylinder>(),        //
    CGO_LINEWIDTH_SZ,                      //
    CGO_WIDTHSCALE_SZ,                     //
    CGO_ENABLE_SZ,                         //
    CGO_DISABLE_SZ,                        //
    fsizeof<cgo::draw::sausage>(),         //
    fsizeof<cgo::draw::custom_cylinder>(), //
    CGO_DOTWIDTH_SZ,                       //
    CGO_ALPHA_TRIANGLE_SZ,                 //
    CGO_ELLIPSOID_SZ,                      //
    CGO_FONT_SZ,                           //
    CGO_FONT_SCALE_SZ,                     //
    CGO_FONT_VERTEX_SZ,                    //
    CGO_FONT_AXES_SZ,                      //
    CGO_CHAR_SZ,                           //
    CGO_INDENT_SZ,                         //
    CGO_ALPHA_SZ,                          //
    CGO_QUADRIC_SZ,                        //
    CGO_CONE_SZ,                           //
    fsizeof<cgo::draw::arrays>(),          //
    CGO_NULL_SZ, CGO_RESET_NORMAL_SZ,      //
    CGO_PICK_COLOR_SZ,                     //
    CGO_NULL_SZ,                           // CGO_DRAW_BUFFERS_SZ no longer used
    fsizeof<cgo::draw::buffers_indexed>(), //
    CGO_BOUNDING_BOX_SZ,                   //
    fsizeof<cgo::draw::buffers_not_indexed>(),             //
    CGO_SPECIAL_SZ,                                        //
    fsizeof<cgo::draw::cylinder_buffers>(),                //
    fsizeof<cgo::draw::shadercylinder>(),                  //
    fsizeof<cgo::draw::shadercylinder2ndcolor>(),          //
    fsizeof<cgo::draw::sphere_buffers>(),                  //
    CGO_ACCESSIBILITY_SZ,                                  //
    CGO_DRAW_TEXTURE_SZ,                                   //
    fsizeof<cgo::draw::textures>(),                        //
    fsizeof<cgo::draw::screen_textures>(),                 //
    CGO_TEX_COORD_SZ,                                      //
    fsizeof<cgo::draw::label>(),                           //
    fsizeof<cgo::draw::labels>(),                          //
    CGO_DRAW_CONNECTOR_SZ,                                 //
    fsizeof<cgo::draw::connectors>(),                      //
    CGO_DRAW_TRILINES_SZ,                                  //
    CGO_UNIFORM3F_SZ,                                      //
    CGO_SPECIAL_WITH_ARG_SZ,                               //
    fsizeof<cgo::draw::line>(),                            //
    fsizeof<cgo::draw::splitline>(),                       //
    fsizeof<cgo::draw::custom>(),                          //
    fsizeof<cgo::draw::vertex_attribute_3f>(),             //
    fsizeof<cgo::draw::vertex_attribute_4ub>(),            //
    fsizeof<cgo::draw::vertex_attribute_1f>(),             //
    fsizeof<cgo::draw::mask_attribute_if_picking>(),       //
    fsizeof<cgo::draw::bind_vbo_for_picking>(),            //
    CGO_VERTEX_BEGIN_LINE_STRIP_SZ,                        //
    CGO_INTERPOLATED_SZ,                                   //
    CGO_VERTEX_CROSS_SZ,                                   //
    fsizeof<cgo::draw::vertex_attribute_4ub_if_picking>(), //
    fsizeof<cgo::draw::custom_cylinder_alpha>(),           //
    CGO_BEZIER_SZ,                                         //
    fsizeof<cgo::draw::bezier_buffers>(),                  //
    CGO_NULL_SZ                                            //
};

/**
 * Get the number of elements in `CGO_sz`
 */
size_t CGO_sz_size()
{
  return sizeof(CGO_sz) / sizeof(*CGO_sz);
}

static float* CGO_add(CGO* I, unsigned c);
static float* CGO_size(CGO* I, int sz);
static int CGOSimpleCylinder(CGO* I, const float* v1, const float* v2,
    const float tube_size, const float* c1, const float* c2, float a1,
    const float a2, const bool interp, const cCylCap cap1, const cCylCap cap2,
    const Pickable* pickcolor2 = nullptr, const bool stick_round_nub = false);
template <typename CylinderT>
static int CGOSimpleCylinder(CGO* I, const CylinderT& cyl, const float a1,
    const float a2, const bool interp, const cCylCap cap1, const cCylCap cap2,
    const Pickable* pickcolor2 = nullptr, const bool stick_round_nub = false);
static int CGOSimpleEllipsoid(CGO* I, const float* v, float vdw,
    const float* n0, const float* n1, const float* n2);
static int CGOSimpleQuadric(CGO* I, const float* v, float vdw, const float* q);
static int CGOSimpleSphere(
    CGO* I, const float* v, float vdw, short sphere_quality);
static int CGOSimpleCone(CGO* I, const float* v1, const float* v2, float r1,
    float r2, const float* c1, const float* c2, cCylCap cap1, cCylCap cap2);
int CGOSimpleCone(CGO* I, const cgo::draw::cone& cone);

/**
 * Inverse function of CGOArrayFromPyListInPlace
 *
 * I: (input) Primitive CGO (may contain CGO_DRAW_ARRAYS)
 *
 * Return: All-float Python list primitive CGO
 */
static PyObject* CGOArrayAsPyList(const CGO* I)
{
  std::vector<float> flat;
  flat.reserve(I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto op = it.op_code();
    auto pc = it.data();
    auto sz = CGO_sz[op];

    flat.push_back(op);

    switch (op) {
    case CGO_BEGIN:
    case CGO_ENABLE:
    case CGO_DISABLE:
    case CGO_SPECIAL:
      // first member int
      flat.push_back(*reinterpret_cast<const int*>(pc));
      ++pc;
      --sz;
      break;
    case CGO_PICK_COLOR:
      assert(sz == 2);
      flat.push_back(*reinterpret_cast<const int*>(pc++)); // atom index
      flat.push_back(*reinterpret_cast<const int*>(pc++)); // bond index
      sz -= 2;
      break;
    case CGO_DRAW_ARRAYS: {
      auto sp = reinterpret_cast<const cgo::draw::arrays*>(pc);
      flat.push_back(sp->mode);
      flat.push_back(sp->arraybits);
      flat.push_back(sp->narrays); // (redundant)
      flat.push_back(sp->nverts);
      pc = sp->get_data();
      sz = sp->get_data_length();
    }
    }

    // float members
    for (; sz; --sz) {
      flat.push_back(*(pc++));
    }
  }

  return PConvToPyObject(flat);
}

PyObject* CGOAsPyList(CGO* I)
{
  PyObject* result;
  result = PyList_New(2);
  PyObject* list = CGOArrayAsPyList(I);
  PyList_SetItem(result, 0, PyInt_FromLong(PyList_Size(list)));
  PyList_SetItem(result, 1, list);
  return (result);
}

static float CPythonVal_PyFloat_AsDouble_From_List(
    void* G, PyObject* list, size_t i)
{
  float out;
  PConvPyFloatToFloat(PyList_GetItem(list, i), &out);
  return out;
}

/**
 * Inverse function of CGOArrayAsPyList
 *
 * list: (input) All-float Python list primitive CGO (may contain
 * CGO_DRAW_ARRAYS) I: (output) empty CGO
 */
static int CGOArrayFromPyListInPlace(PyObject* list, CGO* I)
{
  // sanity check
  if (!list || !PyList_Check(list))
    return false;

  auto G = I->G;

#define GET_FLOAT(i)                                                           \
  ((float) CPythonVal_PyFloat_AsDouble_From_List(I->G, list, i))
#define GET_INT(i) ((int) CPythonVal_PyFloat_AsDouble_From_List(I->G, list, i))

  for (int i = 0, l = PyList_Size(list); i < l;) {
    unsigned op = GET_INT(i++);
    ok_assert(1, op < CGO_sz_size());
    int sz = CGO_sz[op];
    float* fdata = I->add_to_buffer(sz + 1);
    CGO_write_int(fdata, op);

    switch (op) {
    case CGO_STOP:
      // don't increment size for null terminator
      I->c -= 1;
      return true;
    case CGO_BEGIN:
      I->has_begin_end = true;
    case CGO_ENABLE:
    case CGO_DISABLE:
    case CGO_SPECIAL:
      // first member int
      ok_assert(1, i < l);
      CGO_write_int(fdata, GET_INT(i++));
      sz--;
      break;
    case CGO_PICK_COLOR:
      ok_assert(1, i + 1 < l);
      CGO_write_int(fdata, GET_INT(i++)); // atom index
      CGO_write_int(fdata, GET_INT(i++)); // bond index
      sz -= 2;
      break;
    case CGO_DRAW_ARRAYS: {
      // has abstract superclass, need to be constructed!
      ok_assert(1, i + 3 < l);
      auto sp = new (fdata)
          cgo::draw::arrays(GET_INT(i), GET_INT(i + 1), GET_INT(i + 3));

      // sanity check
      int narrays_check = GET_INT(i + 2);
      if (sp->narrays != narrays_check) {
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGO-Warning: narrays mismatch: %d != %d\n", sp->narrays,
            narrays_check ENDFB(I->G);
      }

      // data
      sz = sp->get_data_length();
      sp->floatdata = fdata = I->allocate_in_data_heap(sz);

      i += 4;
    } break;
    }

    // float members
    for (; sz; --sz) {
      ok_assert(1, i < l);
      *(fdata++) = GET_FLOAT(i++);
    }
  }

#undef GET_FLOAT
#undef GET_INT

  return true;

ok_except1:
  PRINTFB(G, FB_CGO, FB_Errors) " %s-Error: Corrupt data\n", __func__ ENDFB(G);
  return false;
}

CGO* CGONewFromPyList(
    PyMOLGlobals* G, PyObject* list, int version, bool shouldCombine)
{
  int ok = true;
  auto I = CGONew(G);
  if (ok)
    ok = (list != nullptr);
  if (ok)
    ok = PyList_Check(list);
  /* TO ENABLE BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if ((version > 0) && (version <= 86)) {
    if (ok)
      ok = PConvFromPyListItem(G, list, 0, I->c);
    if (ok)
      VLACheck(I->op, float, I->c);
    if (ok)
      ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 1), I->op, I->c);
  } else {
    if (ok)
      ok = CGOArrayFromPyListInPlace(PyList_GetItem(list, 1), I);
  }
  if (!ok) {
    CGOFree(I);
  }
  {
    CGO* cgo = nullptr;
    if (shouldCombine && I && I->has_begin_end) {
      cgo = CGOCombineBeginEnd(I, 0);
      CGOFree(I);
    } else {
      cgo = I;
    }
    return cgo;
  }
}

CGO::CGO(PyMOLGlobals* G, int size)
    : G(G)
{
  op = VLACalloc(float, size + 32);
  cgo_shader_ub_color = SettingGet<bool>(G, cSetting_cgo_shader_ub_color);
  cgo_shader_ub_normal = SettingGet<bool>(G, cSetting_cgo_shader_ub_normal);
}

void CGOSetUseShader(CGO* I, int use_shader)
{
  I->use_shader = use_shader;
  if (use_shader) {
    I->cgo_shader_ub_color =
        SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color);
    I->cgo_shader_ub_normal =
        SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal);
  } else {
    I->cgo_shader_ub_color = 0;
    I->cgo_shader_ub_normal = 0;
  }
}
void CGOReset(CGO* I)
{
  I->c = 0;
  I->z_flag = false;
  I->alpha = 1.f;
  I->has_begin_end = false;
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->normal[0] = 0.f;
  I->normal[1] = 0.f;
  I->normal[2] = 1.f;
  I->color[0] = 0.f;
  I->color[1] = 0.f;
  I->color[2] = 1.f;
  I->pickColor[0] = 0;
  I->pickColor[1] = 0;
  I->pickColor[2] = 0;
  I->pickColor[3] = 255;
  I->current_accessibility = 1.f;
}

void CGOFree(CGO*& I, bool withVBOs)
{
  if (I) {
    if (!withVBOs) {
      I->has_draw_buffers = false;
    }
    DeleteP(I);
  }
}

CGO::~CGO()
{
  if (has_draw_buffers) {
    CGOFreeVBOs(this);
  }
  FreeP(i_start);
  VLAFreeP(op);
}

static float* CGO_add(CGO* I, unsigned c)
{
  float* at;
  VLACheck(I->op, float, I->c + c);
  if (!I->op) {
    return nullptr;
  }
  at = I->op + I->c;
  I->c += c;
  return (at);
}

static float* CGO_size(CGO* I, int sz)
{
  float* at;
  VLASize(I->op, float, sz);
  if (!I->op) {
    return nullptr;
  }
  at = I->op + I->c;
  I->c = sz;
  return (at);
}

/*===== Object Creation Routines =======*/

int CGOFromFloatArray(CGO* I, const float* src, int len)
{
  int iarg;
  int ok;
  int all_ok = true;
  int bad_entry = 0;
  int sz;
  int a;
  int cc = 0;
  float val;
  float *pc, *save_pc, *tf;
  VLACheck(I->op, float, I->c + len + 32);
  save_pc = I->op + I->c;
  while (len-- > 0) {
    cc++;
    const auto op = static_cast<unsigned>(*(src++));

    if (op >= CGO_sz_size()) {
      bad_entry = cc;
      break;
    }

    sz = CGO_sz[op];
    if (len < sz)
      break; /* discard short instruction */
    len -= sz;
    pc = save_pc;
    CGO_write_int(pc, op);
    ok = true;
    for (a = 0; a < sz; a++) {
      cc++;
      val = *(src++);
      if (std::abs(val) <= R_SMALL8) {
        val = 0;
      }
      if ((FLT_MAX - val) > 0.0F) { /* make sure we have a real float */
        *(pc++) = val;
      } else {
        *(pc++) = 0.0;
        ok = false;
      }
    }
    if (ok) {
      switch (op) {
      case CGO_END:
      case CGO_VERTEX:
      case CGO_BEGIN:
        I->has_begin_end = true;
      }
      switch (op) { /* now convert any instructions with int arguments */
      case CGO_BEGIN:
      case CGO_ENABLE:
      case CGO_DISABLE:
      case CGO_SPECIAL:
        tf = save_pc + 1;
        iarg = (int) *(tf);
        CGO_write_int(tf, iarg);
        break;
      case CGO_PICK_COLOR:
        tf = save_pc + 1;
        CGO_write_int(tf, int(save_pc[1])); // atom index
        CGO_write_int(tf, int(save_pc[2])); // bond index
        break;
      }
      save_pc = pc;
      I->c += sz + 1;
    } else { /* discard illegal instructions */
      if (all_ok)
        bad_entry = cc;
      all_ok = false;
    }
  }
  return (bad_entry);
}

int CGOBegin(CGO* I, int mode)
{
  float* pc = CGO_add(I, CGO_BEGIN_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_BEGIN);
  CGO_write_int(pc, mode);
  I->has_begin_end = true;
  I->texture[0] = 0.f;
  I->texture[1] = 0.f;
  return true;
}

int CGOEnd(CGO* I)
{
  float* pc = CGO_add(I, CGO_END_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_END);
  I->has_begin_end = true;
  return true;
}

int CGOEnable(CGO* I, int mode)
{
  float* pc = CGO_add(I, CGO_ENABLE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ENABLE);
  CGO_write_int(pc, mode);
  return true;
}

int CGODisable(CGO* I, int mode)
{
  float* pc = CGO_add(I, CGO_DISABLE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_DISABLE);
  CGO_write_int(pc, mode);
  return true;
}

int CGOLinewidth(CGO* I, float v)
{
  float* pc = CGO_add(I, CGO_LINEWIDTH_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_LINEWIDTH);
  *(pc++) = v;
  return true;
}

/**
 * implements special-case operations inside a CGO
 *
 * v: lookup value defined for each special operation (see CGO.h)
 */
int CGOSpecial(CGO* I, int v)
{
  float* pc = CGO_add(I, CGO_SPECIAL_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_SPECIAL);
  CGO_write_int(pc, v);
  return true;
}

/**
 * implements special-case operations with an argument
 * inside a CGO
 *
 * v: lookup value defined for each special operation (see CGO.h)
 * argval : argument value
 */
int CGOSpecialWithArg(CGO* I, int v, float argval)
{
  float* pc = CGO_add(I, CGO_SPECIAL_WITH_ARG_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_SPECIAL_WITH_ARG);
  CGO_write_int(pc, v);
  *pc = argval;
  return true;
}

int CGODotwidth(CGO* I, float v)
{
  float* pc = CGO_add(I, CGO_DOTWIDTH_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_DOTWIDTH);
  *(pc++) = v;
  return true;
}

/* CGOUniform3f - specifies a 3f uniform variable and
   its value.  This function returns the offset of where
   these values are stored inside the CGO float array
   so that they can be accessed and changed from outside
   the CGO.

 */
int CGOUniform3f(CGO* I, int uniform_id, const float* value)
{
  float* pc = CGO_add(I, CGO_UNIFORM3F_SZ + 1);
  if (!pc)
    return 0;
  CGO_write_int(pc, CGO_UNIFORM3F);
  CGO_write_int(pc, uniform_id);
  copy3f(value, pc);
  return pc - I->op;
}

int CGOBoundingBox(CGO* I, const float* min, const float* max)
{
  float* pc = CGO_add(I, CGO_BOUNDING_BOX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_BOUNDING_BOX);
  *(pc++) = *(min);
  *(pc++) = *(min + 1);
  *(pc++) = *(min + 2);
  *(pc++) = *(max);
  *(pc++) = *(max + 1);
  *(pc++) = *(max + 2);
  return true;
}

int CGOAccessibility(CGO* I, float a)
{
  float* pc = CGO_add(I, CGO_ACCESSIBILITY_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ACCESSIBILITY);
  *(pc++) = a;
  return true;
}

int CGODrawTexture(CGO* I, float* worldPos, float* screenMin, float* screenMax,
    float* textExtent)
{
  float* pc = CGO_add(I, CGO_DRAW_TEXTURE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_DRAW_TEXTURE);
  *(pc++) = worldPos[0];
  *(pc++) = worldPos[1];
  *(pc++) = worldPos[2];
  *(pc++) = screenMin[0];
  *(pc++) = screenMin[1];
  *(pc++) = screenMin[2];
  *(pc++) = screenMax[0];
  *(pc++) = screenMax[1];
  *(pc++) = screenMax[2];
  *(pc++) = textExtent[0];
  *(pc++) = textExtent[1];
  *(pc++) = textExtent[2];
  *(pc++) = textExtent[3];
  return true;
}

int CGODrawConnector(CGO* I, float* targetPt3d, float* labelCenterPt3d,
    float text_width, float text_height, float* indentFactor,
    float* screenWorldOffset, float* connectorColor, short relativeMode,
    int draw_flags, float bkgrd_transp, float* bkgrd_color,
    float rel_ext_length, float connectorWidth)
{
  float* pc = CGO_add(I, CGO_DRAW_CONNECTOR_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_DRAW_CONNECTOR);
  *(pc++) = targetPt3d[0];
  *(pc++) = targetPt3d[1];
  *(pc++) = targetPt3d[2];
  *(pc++) = labelCenterPt3d[0];
  *(pc++) = labelCenterPt3d[1];
  *(pc++) = labelCenterPt3d[2];
  *(pc++) = indentFactor[0];
  *(pc++) = indentFactor[1];
  *(pc++) = rel_ext_length; /* place for ext_length relative to height (i.e.,
                               text_height which is total height */
  *(pc++) = screenWorldOffset[0];
  *(pc++) = screenWorldOffset[1];
  *(pc++) = screenWorldOffset[2];
  *(pc++) = text_width;
  *(pc++) = text_height;
  *(pc++) = connectorColor[0];
  *(pc++) = connectorColor[1];
  *(pc++) = connectorColor[2];
  *(pc++) = (float) relativeMode;
  *(pc++) = (float) draw_flags;
  *(pc++) = bkgrd_color[0];
  *(pc++) = bkgrd_color[1];
  *(pc++) = bkgrd_color[2];
  *(pc++) = bkgrd_transp;
  *(pc++) = connectorWidth; // place for label_connector_width
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGODrawLabel(CGO* I, int texture_id, float* targetPos, float* worldPos,
    float* screenWorldOffset, float* screenMin, float* screenMax,
    float* textExtent, short relativeMode)
{
  float* pc = CGO_add(I, CGO_DRAW_LABEL_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_DRAW_LABEL);
  *(pc++) = worldPos[0];
  *(pc++) = worldPos[1];
  *(pc++) = worldPos[2];
  *(pc++) = screenWorldOffset[0];
  *(pc++) = screenWorldOffset[1];
  *(pc++) = screenWorldOffset[2];
  *(pc++) = screenMin[0];
  *(pc++) = screenMin[1];
  *(pc++) = screenMin[2];
  *(pc++) = screenMax[0];
  *(pc++) = screenMax[1];
  *(pc++) = screenMax[2];
  *(pc++) = textExtent[0];
  *(pc++) = textExtent[1];
  *(pc++) = textExtent[2];
  *(pc++) = textExtent[3];
  *(pc++) = (float) relativeMode;
  *(pc++) = targetPos[0];
  *(pc++) = targetPos[1];
  *(pc++) = targetPos[2];
  return true;
}
#endif

#ifdef WITH_UNUSED_FUNCTIONS
int CGOConev(CGO* I, const float* p1, const float* p2, float r1, float r2,
    const float* c1, const float* c2, float cap1, float cap2)
{
  float* pc = CGO_add(I, CGO_CONE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_CONE);
  *(pc++) = *(p1++);
  *(pc++) = *(p1++);
  *(pc++) = *(p1++);
  *(pc++) = *(p2++);
  *(pc++) = *(p2++);
  *(pc++) = *(p2++);
  *(pc++) = r1;
  *(pc++) = r2;
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
  *(pc++) = cap1;
  *(pc++) = cap2;
  return true;
}
#endif

int CGOPickColor(CGO* I, unsigned int index, int bond)
{
  // check if uchar is -1 since extrude does this for masked atoms
  if (index == (unsigned int) -1) {
    bond = cPickableNoPick;
  }
  if (I->current_pick_color_index == index &&
      I->current_pick_color_bond == bond)
    return true;
  float* pc = CGO_add(I, CGO_PICK_COLOR_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_PICK_COLOR);
  CGO_write_uint(pc, index);
  CGO_write_int(pc, bond);
  I->current_pick_color_index = index;
  I->current_pick_color_bond = bond;
  return true;
}

int CGOAlpha(CGO* I, float alpha)
{
  float* pc = CGO_add(I, CGO_ALPHA_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ALPHA);
  *(pc++) = alpha;
  I->alpha = alpha;
  return true;
}

int CGOSphere(CGO* I, const float* v1, float r)
{
  float* pc = CGO_add(I, CGO_SPHERE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_SPHERE);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = r;
  return true;
}

int CGOEllipsoid(CGO* I, const float* v1, float r, const float* n1,
    const float* n2, const float* n3)
{
  float* pc = CGO_add(I, CGO_ELLIPSOID_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ELLIPSOID);

  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = r;
  *(pc++) = *(n1++);
  *(pc++) = *(n1++);
  *(pc++) = *(n1++);
  *(pc++) = *(n2++);
  *(pc++) = *(n2++);
  *(pc++) = *(n2++);
  *(pc++) = *(n3++);
  *(pc++) = *(n3++);
  *(pc++) = *(n3++);
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOQuadric(CGO* I, const float* v, float r, const float* q)
{
  float* pc = CGO_add(I, CGO_QUADRIC_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_QUADRIC);

  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = r;

  *(pc++) = *(q++);
  *(pc++) = *(q++);
  *(pc++) = *(q++);
  *(pc++) = *(q++);
  *(pc++) = *(q++);

  *(pc++) = *(q++);
  *(pc++) = *(q++);
  *(pc++) = *(q++);
  *(pc++) = *(q++);
  *(pc++) = *(q++);
  return true;
}
#endif

void CGOSetZVector(CGO* I, float z0, float z1, float z2)
{
  I->z_flag = true;
  I->z_vector[0] = z0;
  I->z_vector[1] = z1;
  I->z_vector[2] = z2;
  I->z_min = FLT_MAX;
  I->z_max = -FLT_MAX;
}

const static float one_third = 1.0F / 3.0F;
const static float _0 = 0.0F;

int CGOAlphaTriangle(CGO* I, const float* v1, const float* v2, const float* v3,
    const float* n1, const float* n2, const float* n3, const float* c1,
    const float* c2, const float* c3, float a1, float a2, float a3, int reverse)
{
  if (v1 && v2 && v3) {
    float* pc = CGO_add(I, CGO_ALPHA_TRIANGLE_SZ + 1);
    float z = _0;
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_ALPHA_TRIANGLE);
    CGO_write_int(pc, 0); // this is the place for the next triangle in the bin
    *(pc++) = (v1[0] + v2[0] + v3[0]) * one_third;
    *(pc++) = (v1[1] + v2[1] + v3[1]) * one_third;
    *(pc++) = (v1[2] + v2[2] + v3[2]) * one_third;
    if (I->z_flag) {
      float* zv = I->z_vector;
      z = pc[-3] * zv[0] + pc[-2] * zv[1] + pc[-1] * zv[2];
      if (z > I->z_max)
        I->z_max = z;
      if (z < I->z_min)
        I->z_min = z;
    }
    *(pc++) = z;

    if (reverse) {
      *(pc++) = *(v2++); /* vertices @ +5 */
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
    } else {
      *(pc++) = *(v1++); /* vertices @ +5 */
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
    }

    *(pc++) = *(v3++);
    *(pc++) = *(v3++);
    *(pc++) = *(v3++);

    if (reverse) {
      *(pc++) = *(n2++); /* normals @ +14 */
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
    } else {
      *(pc++) = *(n1++); /* normals @ +14 */
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
    }
    *(pc++) = *(n3++);
    *(pc++) = *(n3++);
    *(pc++) = *(n3++);

    if (reverse) {
      *(pc++) = *(c2++); /* colors @ +23 */
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = a2;
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = a1;
    } else {
      *(pc++) = *(c1++); /* colors @ +23 */
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = a1;
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = a2;
    }
    *(pc++) = *(c3++);
    *(pc++) = *(c3++);
    *(pc++) = *(c3++);
    *(pc++) = a3;
  }
  return true;
}

int CGOVertex(CGO* I, float v1, float v2, float v3)
{
  float* pc = CGO_add(I, CGO_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  return true;
}

int CGOVertexv(CGO* I, const float* v)
{
  float* pc = CGO_add(I, CGO_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

int CGOVertexCrossv(CGO* I, const float* v)
{
  float* pc = CGO_add(I, CGO_VERTEX_CROSS_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_VERTEX_CROSS);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOInterpolated(CGO* I, const bool interp)
{
  float* pc = CGO_add(I, CGO_INTERPOLATED_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_INTERPOLATED);
  *(pc++) = interp ? 1.f : 0.f;
  I->interpolated = interp;
  return true;
}
#endif

int CGOColor(CGO* I, float v1, float v2, float v3)
{
  float* pc = CGO_add(I, CGO_COLOR_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_COLOR);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  I->color[0] = v1;
  I->color[1] = v2;
  I->color[2] = v3;
  return true;
}

int CGOColorv(CGO* I, const float* v)
{
  return CGOColor(I, v[0], v[1], v[2]);
}

int CGOTexCoord2f(CGO* I, float v1, float v2)
{
  float* pc = CGO_add(I, CGO_TEX_COORD_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_TEX_COORD);
  *(pc++) = v1;
  *(pc++) = v2;
  I->texture[0] = v1;
  I->texture[1] = v2;
  return true;
}

int CGONormal(CGO* I, float v1, float v2, float v3)
{
  float* pc = CGO_add(I, CGO_NORMAL_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_NORMAL);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  I->normal[0] = v1;
  I->normal[1] = v2;
  I->normal[2] = v3;
  return true;
}

int CGOResetNormal(CGO* I, int mode)
{
  float* pc = CGO_add(I, CGO_RESET_NORMAL_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_RESET_NORMAL);
  CGO_write_int(pc, mode);
  SceneGetResetNormal(I->G, I->normal, mode);
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOFontVertexv(CGO* I, const float* v)
{
  float* pc = CGO_add(I, CGO_FONT_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_FONT_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

int CGOFontVertex(CGO* I, float x, float y, float z)
{
  float* pc = CGO_add(I, CGO_FONT_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_FONT_VERTEX);
  *(pc++) = x;
  *(pc++) = y;
  *(pc++) = z;
  return true;
}
#endif

int CGOFontScale(CGO* I, float v1, float v2)
{
  float* pc = CGO_add(I, CGO_FONT_SCALE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_FONT_SCALE);
  *(pc++) = v1;
  *(pc++) = v2;
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOChar(CGO* I, char c)
{
  float* pc = CGO_add(I, CGO_CHAR_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_CHAR);
  *(pc++) = (float) c;
  return true;
}

int CGOIndent(CGO* I, char c, float dir)
{
  float* pc = CGO_add(I, CGO_INDENT_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_INDENT);
  *(pc++) = (float) c;
  *(pc++) = dir;
  return true;
}

int CGOWrite(CGO* I, const char* str)
{
  float* pc;

  while (*str) {
    pc = CGO_add(I, CGO_CHAR_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(str++);
  }
  return true;
}

int CGOWriteLeft(CGO* I, const char* str)
{
  float* pc;
  const char* s = str;
  while (*s) {
    pc = CGO_add(I, CGO_INDENT_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_INDENT);
    *(pc++) = (float) *(s++);
    *(pc++) = -1.0F;
  }
  s = str;
  while (*s) {
    pc = CGO_add(I, CGO_CHAR_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(s++);
  }
  return true;
}

int CGOWriteIndent(CGO* I, const char* str, float indent)
{
  float* pc;
  const char* s = str;
  while (*s) {
    pc = CGO_add(I, CGO_INDENT_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_INDENT);
    *(pc++) = (float) *(s++);
    *(pc++) = indent;
  }
  s = str;
  while (*s) {
    pc = CGO_add(I, CGO_CHAR_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(s++);
  }
  return true;
}
#endif

int CGONormalv(CGO* I, const float* v)
{
  float* pc = CGO_add(I, CGO_NORMAL_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_NORMAL);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

/**
 * Add a null terminator to the CGO buffer, but don't increment
 * the size variable (CGO::c).
 */
int CGOStop(CGO* I)
{
#define CGO_STOP_ZEROS 1

  float* pc = CGO_size(I, I->c + CGO_STOP_ZEROS);
  if (!pc)
    return false;
  UtilZeroMem(pc, sizeof(float) * CGO_STOP_ZEROS);
  I->c -= CGO_STOP_ZEROS;
  return true;
}

int CGOCheckComplex(CGO* I)
{
  int fc = 0;
  const SphereRec* sp = I->G->Sphere->Sphere[1];

  /* stick_quality needs to match *every* CGO? */
  auto nEdge = SettingGet<int>(I->G, cSetting_stick_quality);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_CYLINDER:
    case CGO_CONE:
    case CGO_SAUSAGE:
    case CGO_CUSTOM_CYLINDER:
    case CGO_CUSTOM_CYLINDER_ALPHA:
      fc += 3 * (3 + (nEdge + 1) * 9) + 9;
      break;
    case CGO_ELLIPSOID:
    case CGO_QUADRIC:
    case CGO_SPHERE:
      fc += (sp->NVertTot * 6) + (sp->NStrip * 3) + 3;
      break;
    case CGO_DRAW_ARRAYS: {
      cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      fc += sp->nverts;
    } break;
    case CGO_DRAW_BUFFERS_INDEXED: {
      cgo::draw::buffers_indexed* sp = reinterpret_cast<decltype(sp)>(pc);
      switch (sp->mode) {
      case GL_TRIANGLES:
        fc += sp->nindices / 3;
        break;
      case GL_LINES:
        fc += sp->nindices / 2;
        break;
      }
    } break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED: {
      cgo::draw::buffers_not_indexed* sp = reinterpret_cast<decltype(sp)>(pc);
      switch (sp->mode) {
      case GL_TRIANGLES:
        fc += sp->nverts / 3;
        break;
      case GL_LINES:
        fc += sp->nverts / 2;
        break;
      }
    } break;
    case CGO_DRAW_SPHERE_BUFFERS: {
      cgo::draw::sphere_buffers* sp = reinterpret_cast<decltype(sp)>(pc);
      fc += sp->num_spheres * VerticesPerSphere();
    } break;
    case CGO_DRAW_CYLINDER_BUFFERS: {
      cgo::draw::cylinder_buffers* sp = reinterpret_cast<decltype(sp)>(pc);
      fc += sp->num_cyl * NumVerticesPerCylinder();
    } break;
    }
  }
  return (fc);
}

int CGOPreloadFonts(CGO* I)
{
  int ok = true;
  int font_seen = false;
  int font_id;

  auto blocked = PAutoBlock(I->G);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
    case CGO_FONT:
      ok = ok && (VFontLoad(I->G, 1.0, 1, 1, true));
      font_seen = true;
      break;
    case CGO_CHAR:
      if (!font_seen) {
        font_id = VFontLoad(I->G, 1.0, 1, 1, true);
        ok = ok && font_id;
        font_seen = true;
      }
      break;
    }
  }
  if (blocked)
    PUnblock(I->G);
  return (ok);
}

int CGOCheckForText(CGO* I)
{
  int fc = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
    case CGO_FONT:
    case CGO_FONT_AXES:
    case CGO_FONT_SCALE:
      fc++;
      break;
    case CGO_INDENT:
    case CGO_FONT_VERTEX:
      fc++;
      break;
    case CGO_CHAR:
      fc += 3 + 2 * 3 * 10; /* est 10 lines per char */
      break;
    }
  }
  PRINTFD(I->G, FB_CGO)
  " CGOCheckForText-Debug: %d\n", fc ENDFD;

  return (fc);
}

CGO* CGODrawText(const CGO* I, int est, float* camera)
{ /* assumes blocked intepreter */
  CGO* cgo;
  int font_id = 0;
  char text[2] = " ";
  float pos[] = {0.0F, 0.0F, 0.0F};
  float axes[] = {1.0F, 0.0F, 0.0F, 0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 1.0F};
  float scale[2] = {1.0, 1.0};

  cgo = CGONewSized(I->G, I->c + est);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_FONT:
      break;
    case CGO_FONT_AXES:
      break;
    case CGO_FONT_SCALE:
      scale[0] = pc[0];
      scale[1] = pc[1];
      break;
    case CGO_FONT_VERTEX:
      copy3f(pc, pos);
      break;
    case CGO_INDENT:
      text[0] = (unsigned char) *pc;
      VFontIndent(I->G, font_id, text, pos, scale, axes, pc[1]);
      break;
    case CGO_CHAR:
      if (!font_id) {
        font_id = VFontLoad(I->G, 1.0, 1, 1, false);
      }
      text[0] = (unsigned char) *pc;
      VFontWriteToCGO(I->G, font_id, cgo, text, pos, scale, axes, cgo->color);
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  if (cgo && cgo->has_begin_end) {
    /* this is mainly for VFontWriteToCGO() that still creates CGOBegin/CGOEnd
     */
    if (cgo && cgo->has_begin_end) {
      CGO* convertcgo = nullptr;
      convertcgo = CGOCombineBeginEnd(cgo, 0);
      CGOFree(cgo);
      cgo = convertcgo;
    }
  }
  return (cgo);
}

static void CGOAddVertexToDrawArrays(CGO* cgo, int pl, int plc, int pla,
    const float* vertex, short notHaveValue, float* vertexVals,
    float* normalVals, float* colorVals, float* pickColorVals,
    float* accessibilityVals)
{
  float* tmp_ptr;
  if (notHaveValue & CGO_NORMAL_ARRAY) {
    if (pl) {
      tmp_ptr = &normalVals[pl - 3];
      copy3f(tmp_ptr, &normalVals[pl]);
    } else {
      copy3f(cgo->normal, &normalVals[pl]);
    }
  }
  if (notHaveValue & CGO_COLOR_ARRAY) {
    if (plc) {
      tmp_ptr = &colorVals[plc - 4];
      copy4f(tmp_ptr, &colorVals[plc]);
    } else {
      copy3f(&colorVals[plc], cgo->color);
      colorVals[plc + 3] = cgo->alpha;
    }
  }
  if (pickColorVals) {
    CGO_put_uint(pickColorVals + pla * 2, cgo->current_pick_color_index);
    CGO_put_int(pickColorVals + pla * 2 + 1, cgo->current_pick_color_bond);
  }
  if (accessibilityVals) {
    accessibilityVals[pla] = cgo->current_accessibility;
  }
  copy3f(vertex, &vertexVals[pl]);
}

bool CGOCombineBeginEnd(CGO** I, bool do_not_split_lines)
{
  CGO* cgo = CGOCombineBeginEnd(*I, 0, do_not_split_lines);
  CGOFree(*I);
  *I = cgo;
  return (cgo != nullptr);
}

/**
 * Converts Begin/End blocks into CGO_DRAW_ARRAYS operations.
 *
 * @param I Primitive CGO to convert
 * @param est Output CGO size estimate (size of buffer to "reserve")
 * @param do_not_split_lines ???
 *
 * @return New converted CGO with has_begin_end=false
 */
CGO* CGOCombineBeginEnd(const CGO* I, int est, bool do_not_split_lines)
{
  CGO* cgo;

  int ok = true;
  if (!I)
    return nullptr;
  cgo = CGONewSized(I->G, 0);
  ok &= cgo ? true : false;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_END:
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGOCombineBeginEnd: op=0x%02x encountered without CGO_BEGIN\n",
          op ENDFB(I->G);
      break;
    case CGO_BEGIN: {
      float firstColor[3], firstAlpha;
      char hasFirstColor = 0, hasFirstAlpha = 0;
      int nverts = 0, damode = CGO_VERTEX_ARRAY, err = 0;

      // read int argument of the BEGIN operation
      int mode = CGO_get_int(it.data());
      ++it;

      // we want to iterate twice over the BEGIN/END block
      auto it2 = it;

      // first iteration over BEGIN/END block (consumes 'it')
      for (; !err && it != CGO_END; ++it) {
        auto pc = it.data();
        switch (it.op_code()) {
        case CGO_DRAW_ARRAYS:
        case CGO_STOP:
          PRINTFB(I->G, FB_CGO, FB_Errors)
          " CGO-Error: CGOCombineBeginEnd: invalid op=0x%02x inside "
          "BEGIN/END\n",
              it.op_code() ENDFB(I->G);
          err = true;
          continue;
        case CGO_NORMAL:
          damode |= CGO_NORMAL_ARRAY;
          break;
        case CGO_COLOR:
          if (!nverts) {
            hasFirstColor = 1;
            copy3f(pc, firstColor);
          } else {
            hasFirstColor = 0;
            damode |= CGO_COLOR_ARRAY;
          }
          copy3f(pc, cgo->color);
          break;
        case CGO_PICK_COLOR:
          damode |= CGO_PICK_COLOR_ARRAY;
          break;
        case CGO_ACCESSIBILITY:
          damode |= CGO_ACCESSIBILITY_ARRAY;
          break;
        case CGO_VERTEX:
          nverts++;
          break;
        case CGO_ALPHA:
          cgo->alpha = *pc;
          if (!nverts) {
            hasFirstAlpha = 1;
            firstAlpha = cgo->alpha;
          } else {
            hasFirstAlpha = 0;
            damode |= CGO_COLOR_ARRAY;
          }
          break;
        case CGO_LINE:
          nverts += 2;
          break;
        case CGO_SPLITLINE: {
          auto splitline = reinterpret_cast<const cgo::draw::splitline*>(pc);
          if (do_not_split_lines ||
              (splitline->flags & cgo::draw::splitline::equal_colors)) {
            nverts += 2;
          } else {
            nverts += 4;
          }
        } break;
        }
      }
      if (nverts > 0 && !err) {
        int pl = 0, plc = 0, pla = 0;
        float *vertexVals = nullptr, *normalVals = nullptr,
              *colorVals = nullptr, *pickColorVals = nullptr,
              *accessibilityVals = nullptr;

        if (hasFirstAlpha || hasFirstColor) {
          if (hasFirstAlpha) {
            CGOAlpha(cgo, firstAlpha);
          }
          if (hasFirstColor) {
            CGOColorv(cgo, firstColor);
          }
        }

        float* nxtVals = cgo->add<cgo::draw::arrays>(mode, damode, nverts);
        ok &= nxtVals ? true : false;
        if (!ok)
          continue;

        assert(damode & CGO_VERTEX_ARRAY);
        vertexVals = nxtVals;
        nxtVals += VERTEX_POS_SIZE * nverts;

        if (damode & CGO_NORMAL_ARRAY) {
          normalVals = nxtVals;
          nxtVals += nverts * VERTEX_NORMAL_SIZE;
        }

        if (damode & CGO_COLOR_ARRAY) {
          colorVals = nxtVals;
          nxtVals += nverts * VERTEX_COLOR_SIZE;
        }

        if (damode & CGO_PICK_COLOR_ARRAY) {
          pickColorVals = nxtVals + VERTEX_PICKCOLOR_RGBA_SIZE * nverts;
          nxtVals += nverts * VERTEX_PICKCOLOR_SIZE;
        }

        if (damode & CGO_ACCESSIBILITY_ARRAY) {
          accessibilityVals = nxtVals;
          nxtVals += nverts * VERTEX_ACCESSIBILITY_SIZE;
        }

        auto notHaveValue = damode;

        // second iteration (with copy of iterator, doesn't consume 'it')
        for (; ok && it2 != CGO_END; ++it2) {
          auto pc = it2.data();
          switch (it2.op_code()) {
          case CGO_NORMAL:
            copy3f(pc, &normalVals[pl]);
            notHaveValue &= ~CGO_NORMAL_ARRAY;
            break;
          case CGO_COLOR:
            if (colorVals) {
              copy3f(pc, &colorVals[plc]);
              colorVals[plc + 3] = cgo->alpha;
              notHaveValue &= ~CGO_COLOR_ARRAY;
            }
            copy3f(pc, cgo->color);
            break;
          case CGO_PICK_COLOR:
            cgo->current_pick_color_index = CGO_get_uint(pc);
            cgo->current_pick_color_bond = CGO_get_int(pc + 1);
            notHaveValue &= ~CGO_PICK_COLOR_ARRAY;
            break;
          case CGO_ACCESSIBILITY:
            cgo->current_accessibility = pc[0];
            break;
          case CGO_SPLITLINE: {
            auto splitline = reinterpret_cast<const cgo::draw::splitline*>(pc);
            float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
                CONVERT_COLOR_VALUE(splitline->color2[1]),
                CONVERT_COLOR_VALUE(splitline->color2[2])};
            if (do_not_split_lines ||
                (splitline->flags & cgo::draw::splitline::equal_colors)) {
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex1,
                  notHaveValue, vertexVals, normalVals, colorVals,
                  pickColorVals, accessibilityVals);
              pl += 3;
              plc += 4;
              pla++;
              notHaveValue = damode;

              if (!(splitline->flags & cgo::draw::splitline::equal_colors)) {
                if (colorVals) {
                  copy3f(color2, &colorVals[plc]);
                  colorVals[plc + 3] = cgo->alpha;
                  notHaveValue = notHaveValue & ~CGO_COLOR_ARRAY;
                }
                copy3f(color2, cgo->color);
              }
              if (pickColorVals) {
                cgo->current_pick_color_index = splitline->index;
                cgo->current_pick_color_bond = splitline->bond;
                notHaveValue = notHaveValue & ~CGO_PICK_COLOR_ARRAY;
              }
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex2,
                  notHaveValue, vertexVals, normalVals, colorVals,
                  pickColorVals, accessibilityVals);
              pl += 3;
              plc += 4;
              pla++;
              notHaveValue = damode;
            } else {
              float mid[3];
              add3f(splitline->vertex1, splitline->vertex2, mid);
              mult3f(mid, .5f, mid);
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex1,
                  notHaveValue, vertexVals, normalVals, colorVals,
                  pickColorVals, accessibilityVals);
              notHaveValue = damode;
              pl += 3;
              plc += 4;
              pla++;
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, mid, notHaveValue,
                  vertexVals, normalVals, colorVals, pickColorVals,
                  accessibilityVals);
              pl += 3;
              plc += 4;
              pla++;
              if (colorVals) {
                copy3f(color2, &colorVals[plc]);
                colorVals[plc + 3] = cgo->alpha;
                notHaveValue = notHaveValue & ~CGO_COLOR_ARRAY;
              }
              copy3f(color2, cgo->color);
              if (pickColorVals) {
                cgo->current_pick_color_index = splitline->index;
                cgo->current_pick_color_bond = splitline->bond;
                notHaveValue = notHaveValue & ~CGO_PICK_COLOR_ARRAY;
              }
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, mid, notHaveValue,
                  vertexVals, normalVals, colorVals, pickColorVals,
                  accessibilityVals);
              notHaveValue = damode;
              pl += 3;
              plc += 4;
              pla++;
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex2,
                  notHaveValue, vertexVals, normalVals, colorVals,
                  pickColorVals, accessibilityVals);
              pl += 3;
              plc += 4;
              pla++;
              notHaveValue = damode;
            }
          } break;
          case CGO_LINE: {
            auto line = reinterpret_cast<const cgo::draw::line*>(pc);
            CGOAddVertexToDrawArrays(cgo, pl, plc, pla, line->vertex1,
                notHaveValue, vertexVals, normalVals, colorVals, pickColorVals,
                accessibilityVals);
            pl += 3;
            plc += 4;
            pla++;
            notHaveValue = damode;
            CGOAddVertexToDrawArrays(cgo, pl, plc, pla, line->vertex2,
                notHaveValue, vertexVals, normalVals, colorVals, pickColorVals,
                accessibilityVals);
            pl += 3;
            plc += 4;
            pla++;
          } break;
          case CGO_VERTEX:
            CGOAddVertexToDrawArrays(cgo, pl, plc, pla, pc, notHaveValue,
                vertexVals, normalVals, colorVals, pickColorVals,
                accessibilityVals);
            pl += 3;
            plc += 4;
            pla++;
            notHaveValue = damode;
            break;
          case CGO_ALPHA:
            // in case we're before CGO_COLOR
            cgo->alpha = *pc;
            if (colorVals) {
              // in case we're after CGO_COLOR
              colorVals[plc + 3] = *pc;
            }
            break;
          }
        }
      }
    } break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      cgo->add_to_cgo(op, pc);
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  if (ok) {
    ok &= CGOStop(cgo);
    if (ok) {
      cgo->use_shader = I->use_shader;
      if (cgo->use_shader) {
        cgo->cgo_shader_ub_color =
            SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
        cgo->cgo_shader_ub_normal =
            SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
      }
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

/**
 * Release all shader resources from this CGO
 */
void CGOFreeVBOs(CGO* I)
{
  constexpr bool freevbos = true;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
    case CGO_DRAW_TRILINES: {
      unsigned buf = it.cast<cgo::draw::trilines>()->buffer;
      if (freevbos)
        I->G->ShaderMgr->AddVBOToFree(buf);
    } break;
    case CGO_DRAW_CUSTOM: {
      auto sp = it.cast<cgo::draw::custom>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->iboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    } break;
    case CGO_DRAW_SPHERE_BUFFERS: {
      auto sp = it.cast<cgo::draw::sphere_buffers>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    } break;
    case CGO_DRAW_LABELS: {
      auto sp = it.cast<cgo::draw::labels>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    } break;
    case CGO_DRAW_TEXTURES: {
      auto sp = it.cast<cgo::draw::textures>();
      if (freevbos)
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
    } break;
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS: {
      auto sp = it.cast<cgo::draw::screen_textures>();
      if (freevbos)
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
    } break;
    case CGO_DRAW_CYLINDER_BUFFERS: {
      auto sp = it.cast<cgo::draw::cylinder_buffers>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->iboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    } break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED: {
      auto sp = it.cast<cgo::draw::buffers_not_indexed>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    } break;
    case CGO_DRAW_BUFFERS_INDEXED: {
      auto sp = it.cast<cgo::draw::buffers_indexed>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffers({sp->vboid, sp->iboid, sp->pickvboid});
      }
    } break;
    case CGO_DRAW_CONNECTORS: {
      auto sp = it.cast<cgo::draw::connectors>();
      if (freevbos)
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
    } break;
    }
  }
}

#define set_min_max(mn, mx, pt)                                                \
  {                                                                            \
    if (mn[0] > *pt)                                                           \
      mn[0] = *pt;                                                             \
    if (mn[1] > *(pt + 1))                                                     \
      mn[1] = *(pt + 1);                                                       \
    if (mn[2] > *(pt + 2))                                                     \
      mn[2] = *(pt + 2);                                                       \
    if (mx[0] < *pt)                                                           \
      mx[0] = *pt;                                                             \
    if (mx[1] < *(pt + 1))                                                     \
      mx[1] = *(pt + 1);                                                       \
    if (mx[2] < *(pt + 2))                                                     \
      mx[2] = *(pt + 2);                                                       \
  }

struct CGOCount
{
  int num_total_vertices{};
  int num_total_indexes{};
  int num_total_vertices_lines{};
  int num_total_indexes_lines{};
  int num_total_vertices_points{};
};

static CGOCount CGOCountNumVertices(const CGO* I);

#ifdef WITH_UNUSED_FUNCTIONS
void CGOCountNumVerticesDEBUG(const CGO* I)
{
  auto count = CGOCountNumVertices(I);
  CGOCountNumVertices(I);
  printf("CGOCountNumVerticesDEBUG: num_total_vertices=%d num_total_indexes=%d "
         "num_total_vertices_lines=%d num_total_indexes_lines=%d "
         "num_total_vertices_points=%d\n",
      count.num_total_vertices, count.num_total_indexes,
      count.num_total_vertices_lines, count.num_total_indexes_lines,
      count.num_total_vertices_points);
}
#endif

static CGOCount CGOCountNumVertices(const CGO* I)
{
  CGOCount count{};
  short err = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    err = 0;
    switch (op) {
    case CGO_DRAW_ARRAYS: {
      const auto sp = it.cast<cgo::draw::arrays>();
      short shouldCompress = false, shouldCompressLines = false,
            shouldCompressPoints = false;
      switch (sp->mode) {
      case GL_TRIANGLE_FAN:
      case GL_TRIANGLE_STRIP:
      case GL_TRIANGLES:
        shouldCompress = true;
        break;
      case GL_LINES:
      case GL_LINE_STRIP:
      case GL_LINE_LOOP:
        shouldCompressLines = true;
        break;
      case GL_POINTS:
        shouldCompressPoints = true;
        break;
      default:
        break;
      }
      if (!shouldCompress && !shouldCompressLines && !shouldCompressPoints) {
        // skipped
      } else if (shouldCompressLines) {
        count.num_total_vertices_lines += sp->nverts;
        switch (sp->mode) {
        case GL_LINE_LOOP:
          count.num_total_indexes_lines += 2 * sp->nverts;
          break;
        case GL_LINE_STRIP:
          count.num_total_indexes_lines += 2 * (sp->nverts - 1);
          break;
        case GL_LINES:
          count.num_total_indexes_lines += sp->nverts;
          break;
        }
      } else if (shouldCompress) {
        count.num_total_vertices += sp->nverts;
        switch (sp->mode) {
        case GL_TRIANGLE_FAN:
          count.num_total_indexes += 3 * (sp->nverts - 2);
          break;
        case GL_TRIANGLE_STRIP:
          count.num_total_indexes += 3 * (sp->nverts - 2);
          break;
        case GL_TRIANGLES:
          count.num_total_indexes += sp->nverts;
          break;
        }
      } else if (shouldCompressPoints) {
        count.num_total_vertices_points += sp->nverts;
      }
    } break;
    case CGO_END:
      if (!err) {
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOCountNumVertices: CGO_END encountered, should call "
        "CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
        err = true;
      }
    case CGO_VERTEX:
      if (!err) {
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOCountNumVertices: CGO_VERTEX encountered, should call "
        "CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
        err = true;
      }
    case CGO_BEGIN:
      if (!err) {
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOCountNumVertices: CGO_BEGIN encountered, should call "
        "CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
        err = true;
      }
    default:
      break;
    }
  }
  return count;
}

/**
 * @param[in] I Primitive CGO
 * @return number of vertices, or -1 on error
 */
static int CGOCountNumVerticesForScreen(const CGO* I)
{
  auto G = I->G;
  constexpr int MODE_INVALID = -1;
  int begin_mode = MODE_INVALID;
  int num_total_indexes = 0;
  int nverts = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    // sanity check: non-primitive operations
    switch (op) {
    case CGO_DRAW_ARRAYS:
      WARN_UNEXPECTED_OPERATION(G, op);
      return -1;
    }

    if (begin_mode == MODE_INVALID) {
      switch (op) {
      case CGO_BEGIN:
        begin_mode = it.cast<cgo::draw::begin>()->mode;
        break;
      case CGO_VERTEX:
      case CGO_END:
        WARN_UNEXPECTED_OPERATION(G, op);
        return -1;
      }
    } else {
      switch (op) {
      case CGO_BEGIN:
        WARN_UNEXPECTED_OPERATION(G, op);
        return -1;
      case CGO_VERTEX:
        ++nverts;
        break;
      case CGO_END:
        switch (begin_mode) {
        case GL_TRIANGLE_FAN:
        case GL_TRIANGLE_STRIP:
          num_total_indexes += 3 * (nverts - 2);
          break;
        case GL_TRIANGLES:
          num_total_indexes += nverts;
          break;
        default:
          // not implemented
          assert(false);
        }
        begin_mode = MODE_INVALID;
        nverts = 0;
        break;
      }
    }
  }

  return num_total_indexes;
}

static void SetVertexValuesForVBO(PyMOLGlobals* G, CGO* cgo, int pl, int plc,
    int cnt, int incr, const float* vertexValsDA, const float* normalValsDA,
    const float* colorValsDA, const float* pickColorValsDA, float* vertexVals,
    uchar* normalValsC, float* normalVals, uchar* colorValsUC, float* colorVals,
    float* pickColorVals, float* accessibilityVals = nullptr,
    const float* accessibilityValsDA = nullptr)
{
  int pl2 = pl + 1, pl3 = pl + 2;
  int pln1 = VAR_FOR_NORMAL, pln2 = VAR_FOR_NORMAL + 1,
      pln3 = VAR_FOR_NORMAL + 2;
  int plc2 = plc + 1, plc3 = plc + 2, plc4 = plc + 3;
  int c, c2, c3;
  int cc, cc2, cc3, cc4;
  int pcc = incr * 2, pcco = cnt * 2;
  c = cnt * 3;
  c2 = c + 1;
  c3 = c + 2;
  cc = cnt * 4;
  cc2 = cc + 1;
  cc3 = cc + 2;
  cc4 = cc + 3;
  vertexVals[pl] = vertexValsDA[c];
  vertexVals[pl2] = vertexValsDA[c2];
  vertexVals[pl3] = vertexValsDA[c3];

  if (normalValsC) {
    if (normalValsDA) {
      normalValsC[pln1] = CLIP_NORMAL_VALUE(normalValsDA[c]);
      normalValsC[pln2] = CLIP_NORMAL_VALUE(normalValsDA[c2]);
      normalValsC[pln3] = CLIP_NORMAL_VALUE(normalValsDA[c3]);
    } else {
      normalValsC[pln1] = CLIP_NORMAL_VALUE(cgo->normal[0]);
      normalValsC[pln2] = CLIP_NORMAL_VALUE(cgo->normal[1]);
      normalValsC[pln3] = CLIP_NORMAL_VALUE(cgo->normal[2]);
    }

#ifdef ALIGN_VBOS_TO_4_BYTE_ARRAYS
    normalValsC[pln3 + 1] = 127;
#endif
  } else {
    if (normalValsDA) {
      normalVals[pln1] = normalValsDA[c];
      normalVals[pln2] = normalValsDA[c2];
      normalVals[pln3] = normalValsDA[c3];
    } else {
      normalVals[pln1] = cgo->normal[0];
      normalVals[pln2] = cgo->normal[1];
      normalVals[pln3] = cgo->normal[2];
    }
  }

  if (colorValsUC) {
    if (colorValsDA) {
      colorValsUC[plc] = CLIP_COLOR_VALUE(colorValsDA[cc]);
      colorValsUC[plc2] = CLIP_COLOR_VALUE(colorValsDA[cc2]);
      colorValsUC[plc3] = CLIP_COLOR_VALUE(colorValsDA[cc3]);
      colorValsUC[plc4] = CLIP_COLOR_VALUE(colorValsDA[cc4]);
    } else {
      colorValsUC[plc] = CLIP_COLOR_VALUE(cgo->color[0]);
      colorValsUC[plc2] = CLIP_COLOR_VALUE(cgo->color[1]);
      colorValsUC[plc3] = CLIP_COLOR_VALUE(cgo->color[2]);
      colorValsUC[plc4] = CLIP_COLOR_VALUE(cgo->alpha);
    }
  } else {
    if (colorValsDA) {
      colorVals[plc] = colorValsDA[cc];
      colorVals[plc2] = colorValsDA[cc2];
      colorVals[plc3] = colorValsDA[cc3];
      colorVals[plc4] = colorValsDA[cc4];
    } else {
      colorVals[plc] = cgo->color[0];
      colorVals[plc2] = cgo->color[1];
      colorVals[plc3] = cgo->color[2];
      colorVals[plc4] = cgo->alpha;
    }
  }
  if (pickColorValsDA) {
    cgo->current_pick_color_index = CGO_get_uint(pickColorValsDA + pcco);
    cgo->current_pick_color_bond = CGO_get_int(pickColorValsDA + pcco + 1);
  }
  CGO_put_uint(pickColorVals + pcc, cgo->current_pick_color_index);
  CGO_put_int(pickColorVals + pcc + 1, cgo->current_pick_color_bond);
  if (accessibilityValsDA) {
    accessibilityVals[pl / 3] = accessibilityValsDA[cnt];
  }
}

struct NormalColorFormatSize {
  VertexFormat normalFormat{};
  VertexFormat colorFormat{};
  std::size_t normalSize{};
  std::size_t colorSize{};
};

static NormalColorFormatSize GetNormalColorFormatSize(PyMOLGlobals* G)
{
  NormalColorFormatSize fmt{};

  fmt.normalFormat =
      VERTEX_NORMAL_SIZE == 3 ? VertexFormat::Float3 : VertexFormat::Float4;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)) {
    fmt.normalFormat = VERTEX_NORMAL_SIZE == 3 ? VertexFormat::Byte3Norm
                                               : VertexFormat::Byte4Norm;
  }
  fmt.normalSize = GetSizeOfVertexFormat(fmt.normalFormat);

  fmt.colorFormat = VertexFormat::Float4;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
    fmt.colorFormat = VertexFormat::UByte4Norm;
  }
  fmt.colorSize = GetSizeOfVertexFormat(fmt.colorFormat);

  return fmt;
}

static int OptimizePointsToVBO(const CGO* I, CGO* cgo,
    int num_total_vertices_points, float* min, float* max,
    short* has_draw_buffer, bool addshaders)
{
  auto G = I->G;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;
  bool has_normals = false, has_colors = false;
  int ok = true;

  cgo->alpha = 1.f;
  cgo->color[0] = 1.f;
  cgo->color[1] = 1.f;
  cgo->color[2] = 1.f;

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal)
             ? 1
             : VERTEX_NORMAL_SIZE;
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(num_total_vertices_points) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * num_total_vertices_points;
  nxtn = 3;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)) {
    normalValsC = (uchar*) normalVals;
    nxtn = 1;
  }
  auto* colorVals = normalVals + nxtn * num_total_vertices_points;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * num_total_vertices_points);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_BOUNDING_BOX:
    case CGO_DRAW_SPHERE_BUFFERS:
    case CGO_DRAW_LABELS:
    case CGO_DRAW_TEXTURES:
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
    case CGO_DRAW_CYLINDER_BUFFERS:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    case CGO_DRAW_BUFFERS_INDEXED:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "ERROR: OptimizePointsToVBO used with unsupported CGO ops" ENDFB(I->G);
      return 0;
      break;
    case CGO_NORMAL:
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      has_normals = true;
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      has_colors = true;
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_ARRAYS: {
      const auto sp = it.cast<cgo::draw::arrays>();
      short shouldCompress = false;
      switch (sp->mode) {
      case GL_POINTS:
        shouldCompress = true;
      default:
        break;
      }
      /*	TODO : DO WE NEED TO COMPENSATE FOR THIS? if (!has_normals &&
        arrays & CGO_NORMAL_ARRAY){ arrays = arrays ^ CGO_NORMAL_ARRAY; narrays
        -= 1;
        }*/
      if (shouldCompress) {
        int cnt, nxtn = VERTEX_POS_SIZE, incr = 0;
        float *vertexValsDA = nullptr, *nxtVals = nullptr,
              *colorValsDA = nullptr, *normalValsDA = nullptr;
        float *pickColorValsDA = nullptr, *pickColorValsTMP;

        nxtVals = vertexValsDA = sp->floatdata;

        for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
          set_min_max(min, max, &vertexValsDA[cnt]);
        }
        if (sp->arraybits & CGO_NORMAL_ARRAY) {
          has_normals = true;
          nxtVals = normalValsDA = vertexValsDA + (nxtn * sp->nverts);
        }
        if (sp->arraybits & CGO_COLOR_ARRAY) {
          has_colors = true;
          nxtVals = colorValsDA = nxtVals + (nxtn * sp->nverts);
          nxtn = VERTEX_COLOR_SIZE;
        }
        if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
          nxtVals = nxtVals + (nxtn * sp->nverts);
          pickColorValsDA = nxtVals + sp->nverts;
          nxtn = VERTEX_PICKCOLOR_SIZE;
        }
        pickColorValsTMP = pickColorVals + (idxpl * 2);
        switch (sp->mode) {
        case GL_POINTS:
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(I->G, cgo, pl, plc, cnt, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP);
            idxpl++;
            pl += 3;
            plc += 4;
          }
          break;
        }
        vpl += sp->nverts;
      }
    } break;
    }
    ok &= !I->G->Interrupt;
  }
  if (ok) {
    short arrays = CGO_VERTEX_ARRAY | CGO_PICK_COLOR_ARRAY;
    auto fmt = GetNormalColorFormatSize(I->G);

    VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();

    BufferDataDesc bufData = {{"a_Vertex", VertexFormat::Float3,
        sizeof(float) * num_total_vertices_points * 3, vertexVals}};

    if (has_normals) {
      bufData.push_back({"a_Normal", fmt.normalFormat,
          num_total_vertices_points * fmt.normalSize, normalVals});
    }
    if (has_colors) {
      bufData.push_back({"a_Color", fmt.colorFormat,
          num_total_vertices_points * fmt.colorSize, colorVals});
    }
    ok = vbo->bufferData(std::move(bufData));

    if (ok && has_colors) {
      arrays |= CGO_COLOR_ARRAY;
    }

    size_t vboid = vbo->get_hash_id();

    if (ok) {
      float* newPickColorVals;
      if (addshaders)
        CGOEnable(cgo, GL_DEFAULT_SHADER);
      newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(
          GL_POINTS, arrays, num_total_vertices_points, vboid);
      CHECKOK(ok, newPickColorVals);
      if (ok && addshaders)
        ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      if (!newPickColorVals)
        I->G->ShaderMgr->freeGPUBuffer(vboid);
      if (ok)
        memcpy(newPickColorVals + num_total_vertices_points, pickColorVals,
            num_total_vertices_points * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(vboid);
    }
  }
  return ok;
  /* END GL_POINTS */
  //    printf("num_total_vertices_points=%d\n", num_total_vertices_points);
}

static void FixPickColorsForLine(float* pick1, float* pick2)
{
  unsigned int p1 = CGO_get_uint(pick1);
  unsigned int p2 = CGO_get_uint(pick2);
  int b1 = CGO_get_int(pick1 + 1);
  int b2 = CGO_get_int(pick2 + 1);
  if (p1 != p2 || b1 != b2) {
    // if the pick colors are different, then pick the first one
    CGO_put_uint(pick1, p2);
    CGO_put_int(pick1 + 1, b2);
  }
}

static void FixPickColorsForTriangle(float* pick1, float* pick2, float* pick3)
{
  unsigned int p1 = CGO_get_uint(pick1);
  unsigned int p2 = CGO_get_uint(pick2);
  unsigned int p3 = CGO_get_uint(pick3);
  int b1 = CGO_get_int(pick1 + 1);
  int b2 = CGO_get_int(pick2 + 1);
  int b3 = CGO_get_int(pick3 + 1);
  if (p1 != p2 || p1 != p3 || p2 != p3 || b1 != b2 || b1 != b3 || b2 != b3) {
    // right now, if the pick colors are different, then pick majority,
    // otherwise, pick first one
    if (p1 == p2 && b1 == b2) {
      CGO_put_uint(pick3, p1);
      CGO_put_int(pick3 + 1, b1);
    } else if (p1 == p3 && b1 == b3) {
      CGO_put_uint(pick2, p1);
      CGO_put_int(pick2 + 1, b1);
    } else if (p2 == p3 && b2 == b3) {
      CGO_put_uint(pick1, p2);
      CGO_put_int(pick1 + 1, b2);
    } else {
      CGO_put_uint(pick2, p1);
      CGO_put_int(pick2 + 1, b1);
      CGO_put_uint(pick3, p1);
      CGO_put_int(pick3 + 1, b1);
    }
  }
}

static bool CGOProcessCGOtoArrays(const CGO* I, CGO* cgo, CGO* addtocgo,
    float* min, float* max, int* ambient_occlusion, float* vertexVals,
    float* normalVals, uchar* normalValsC, float* colorVals, uchar* colorValsUC,
    float* pickColorVals, float* accessibilityVals, bool& has_normals,
    bool& has_colors, bool& has_accessibility)
{
  auto G = I->G;
  int idxpl = 0;
  int pl = 0, plc = 0, vpl = 0;
  int ok = true;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_BOUNDING_BOX:
    case CGO_DRAW_SPHERE_BUFFERS:
    case CGO_DRAW_LABELS:
    case CGO_DRAW_TEXTURES:
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
    case CGO_DRAW_CYLINDER_BUFFERS:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    case CGO_DRAW_BUFFERS_INDEXED: {
      const float* newpc = pc;
      if (addtocgo)
        addtocgo->add_to_cgo(op, newpc);
    } break;
    case CGO_NORMAL:
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      has_normals = true;
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      has_colors = true;
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_ACCESSIBILITY:
      cgo->current_accessibility = pc[0];
      has_accessibility = true;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_ARRAYS: {
      const auto sp = it.cast<cgo::draw::arrays>();
      short shouldCompress = false;
      switch (sp->mode) {
      case GL_TRIANGLE_FAN:
      case GL_TRIANGLE_STRIP:
      case GL_TRIANGLES:
        shouldCompress = true;
      default:
        break;
      }
      if (shouldCompress) {
        int cnt, incr = 0;
        const float* nxtVals = sp->floatdata;
        const float *vertexValsDA = nullptr, *colorValsDA = nullptr,
                    *normalValsDA = nullptr, *pickColorValsDA = nullptr,
                    *accessibilityValsDA = nullptr;

        assert(sp->arraybits & CGO_VERTEX_ARRAY);
        vertexValsDA = nxtVals;
        nxtVals += sp->nverts * VERTEX_POS_SIZE;

        for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
          set_min_max(min, max, &vertexValsDA[cnt]);
        }
        if (sp->arraybits & CGO_NORMAL_ARRAY) {
          has_normals = true;
          normalValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_NORMAL_SIZE;
        }

        if (sp->arraybits & CGO_COLOR_ARRAY) {
          has_colors = true;
          colorValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_COLOR_SIZE;
        }
        if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
          nxtVals += VERTEX_PICKCOLOR_RGBA_SIZE * sp->nverts;
          pickColorValsDA = nxtVals;
          nxtVals += VERTEX_PICKCOLOR_INDEX_SIZE * sp->nverts;
        }
        float* pickColorValsTMP = pickColorVals + (idxpl * 2);
        if (sp->arraybits & CGO_ACCESSIBILITY_ARRAY) {
          if (!(*ambient_occlusion) && incr) {
            for (cnt = 0; cnt < incr; cnt++) {
              /* if ambient_occlusion, need to fill in the array */
              accessibilityVals[cnt] = 1.f;
            }
          }
          (*ambient_occlusion) = 1;
          accessibilityValsDA = nxtVals;
          nxtVals += VERTEX_ACCESSIBILITY_SIZE * sp->nverts;
          has_accessibility = true;
        } else {
          if (*ambient_occlusion) {
            for (cnt = incr; cnt < incr + sp->nverts; cnt++) {
              /* if ambient_occlusion, need to fill in the array */
              accessibilityVals[cnt] = 1.f;
            }
          }
        }
        switch (sp->mode) {
        case GL_TRIANGLES:
          for (cnt = 0; ok && cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            if (incr && (incr % 3) == 0) {
              FixPickColorsForTriangle(pickColorValsTMP + (incr - 3) * 2,
                  pickColorValsTMP + (incr - 2) * 2,
                  pickColorValsTMP + (incr - 1) * 2);
            }
            idxpl++;
            pl += 3;
            plc += 4;
            ok &= !G->Interrupt;
          }
          break;
        case GL_TRIANGLE_STRIP: {
          short flip = 0;
          for (cnt = 2; ok && cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt - (flip ? 0 : 2), incr++,
                vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
                vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            idxpl++;
            pl += 3;
            plc += 4;
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
                vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
                vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            idxpl++;
            pl += 3;
            plc += 4;
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt - (flip ? 2 : 0), incr++,
                vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
                vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            FixPickColorsForTriangle(pickColorValsTMP + (incr - 3) * 2,
                pickColorValsTMP + (incr - 2) * 2,
                pickColorValsTMP + (incr - 1) * 2);
            idxpl++;
            pl += 3;
            plc += 4;
            ok &= !G->Interrupt;
            flip = !flip;
          }
        } break;
        case GL_TRIANGLE_FAN:
          for (cnt = 2; ok && cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(G, cgo, pl, plc, 0, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            idxpl++;
            pl += 3;
            plc += 4;
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
                vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
                vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            idxpl++;
            pl += 3;
            plc += 4;
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP, accessibilityVals, accessibilityValsDA);
            FixPickColorsForTriangle(pickColorValsTMP + (incr - 3) * 2,
                pickColorValsTMP + (incr - 2) * 2,
                pickColorValsTMP + (incr - 1) * 2);
            idxpl++;
            pl += 3;
            plc += 4;
            ok &= !G->Interrupt;
          }
          break;
        }
        vpl += sp->nverts;
      }
    } break;
    }
    ok &= !G->Interrupt;
  }
  ok &= !G->Interrupt;
  return ok;
}

static bool CGOProcessScreenCGOtoArrays(PyMOLGlobals* G, CGO* cgo,
    float* vertexVals, float* texcoordVals, float* colorVals,
    uchar* colorValsUC)
{
  int pl = 0;
  cgo->alpha = 1.f;

  for (auto it = cgo->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_BOUNDING_BOX:
    case CGO_DRAW_SPHERE_BUFFERS:
    case CGO_DRAW_LABELS:
    case CGO_DRAW_TEXTURES:
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
    case CGO_DRAW_CYLINDER_BUFFERS:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_ARRAYS:
    case CGO_ACCESSIBILITY:
      WARN_UNEXPECTED_OPERATION(G, op);
      return false;
    case CGO_NORMAL:
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      break;
    case CGO_TEX_COORD:
      cgo->texture[0] = *pc;
      cgo->texture[1] = *(pc + 1);
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_BEGIN: {
      int mode = it.cast<cgo::draw::begin>()->mode;
      int nverts = 0, ipl = 0;
      cgo->texture[0] = cgo->texture[1] = 0.f;

      for (++it;; ++it) {
        if (it.is_stop()) {
          WARN_UNEXPECTED_OPERATION(G, CGO_STOP);
          return false;
        }

        const auto op = it.op_code();
        if (op == CGO_END) {
          break;
        }

        const auto pc = it.data();

        switch (op) {
        case CGO_TEX_COORD:
          cgo->texture[0] = *pc;
          cgo->texture[1] = *(pc + 1);
          break;
        case CGO_COLOR:
          cgo->color[0] = *pc;
          cgo->color[1] = *(pc + 1);
          cgo->color[2] = *(pc + 2);
          break;
        case CGO_ALPHA:
          cgo->alpha = *pc;
          break;
        case CGO_VERTEX: {
          switch (mode) {
          case GL_TRIANGLES: {
            int vpl = pl * 3, tpl = pl * 2, cpl = pl * 4;
            vertexVals[vpl] = *pc;
            vertexVals[vpl + 1] = *(pc + 1);
            vertexVals[vpl + 2] = *(pc + 2);
            texcoordVals[tpl] = cgo->texture[0];
            texcoordVals[tpl + 1] = cgo->texture[1];
            if (colorValsUC) {
              colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
              colorValsUC[cpl + 1] = CLIP_COLOR_VALUE(cgo->color[1]);
              colorValsUC[cpl + 2] = CLIP_COLOR_VALUE(cgo->color[2]);
              colorValsUC[cpl + 3] = CLIP_COLOR_VALUE(cgo->alpha);
            } else {
              colorVals[cpl] = cgo->color[0];
              colorVals[cpl + 1] = cgo->color[1];
              colorVals[cpl + 2] = cgo->color[2];
              colorVals[cpl + 3] = cgo->alpha;
            }
            pl++;
          } break;
          case GL_TRIANGLE_STRIP: {
            int vpl = pl * 3, tpl = pl * 2, cpl = pl * 4;
            if (ipl < 3) {
              vertexVals[vpl] = *pc;
              vertexVals[vpl + 1] = *(pc + 1);
              vertexVals[vpl + 2] = *(pc + 2);
              texcoordVals[tpl] = cgo->texture[0];
              texcoordVals[tpl + 1] = cgo->texture[1];
              if (colorValsUC) {
                colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
                colorValsUC[cpl + 1] = CLIP_COLOR_VALUE(cgo->color[1]);
                colorValsUC[cpl + 2] = CLIP_COLOR_VALUE(cgo->color[2]);
                colorValsUC[cpl + 3] = CLIP_COLOR_VALUE(cgo->alpha);
              } else {
                colorVals[cpl] = cgo->color[0];
                colorVals[cpl + 1] = cgo->color[1];
                colorVals[cpl + 2] = cgo->color[2];
                colorVals[cpl + 3] = cgo->alpha;
              }
              pl++;
              ipl++;
            } else {
              vertexVals[vpl] = vertexVals[vpl - 6];
              vertexVals[vpl + 1] = vertexVals[vpl - 5];
              vertexVals[vpl + 2] = vertexVals[vpl - 4];
              texcoordVals[tpl] = texcoordVals[tpl - 4];
              texcoordVals[tpl + 1] = texcoordVals[tpl - 3];
              if (colorValsUC) {
                colorValsUC[cpl] = colorValsUC[cpl - 8];
                colorValsUC[cpl + 1] = colorValsUC[cpl - 7];
                colorValsUC[cpl + 2] = colorValsUC[cpl - 6];
                colorValsUC[cpl + 3] = colorValsUC[cpl - 5];
              } else {
                colorVals[cpl] = colorVals[cpl - 8];
                colorVals[cpl + 1] = colorVals[cpl - 7];
                colorVals[cpl + 2] = colorVals[cpl - 6];
                colorVals[cpl + 3] = colorVals[cpl - 5];
              }
              pl++;
              vpl += 3;
              tpl += 2;
              cpl += 4;
              ipl++;
              vertexVals[vpl] = vertexVals[vpl - 6];
              vertexVals[vpl + 1] = vertexVals[vpl - 5];
              vertexVals[vpl + 2] = vertexVals[vpl - 4];
              texcoordVals[tpl] = texcoordVals[tpl - 4];
              texcoordVals[tpl + 1] = texcoordVals[tpl - 3];
              if (colorValsUC) {
                colorValsUC[cpl] = colorValsUC[cpl - 8];
                colorValsUC[cpl + 1] = colorValsUC[cpl - 7];
                colorValsUC[cpl + 2] = colorValsUC[cpl - 6];
                colorValsUC[cpl + 3] = colorValsUC[cpl - 5];
              } else {
                colorVals[cpl] = colorVals[cpl - 8];
                colorVals[cpl + 1] = colorVals[cpl - 7];
                colorVals[cpl + 2] = colorVals[cpl - 6];
                colorVals[cpl + 3] = colorVals[cpl - 5];
              }
              pl++;
              vpl += 3;
              tpl += 2;
              cpl += 4;
              ipl++;
              vertexVals[vpl] = *pc;
              vertexVals[vpl + 1] = *(pc + 1);
              vertexVals[vpl + 2] = *(pc + 2);
              texcoordVals[tpl] = cgo->texture[0];
              texcoordVals[tpl + 1] = cgo->texture[1];
              if (colorValsUC) {
                colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
                colorValsUC[cpl + 1] = CLIP_COLOR_VALUE(cgo->color[1]);
                colorValsUC[cpl + 2] = CLIP_COLOR_VALUE(cgo->color[2]);
                colorValsUC[cpl + 3] = CLIP_COLOR_VALUE(cgo->alpha);
              } else {
                colorVals[cpl] = cgo->color[0];
                colorVals[cpl + 1] = cgo->color[1];
                colorVals[cpl + 2] = cgo->color[2];
                colorVals[cpl + 3] = cgo->alpha;
              }
              pl++;
              ipl++;
            }
          } break;
          case GL_TRIANGLE_FAN:
          default:
            printf(
                "CGOProcessScreenCGOtoArrays: WARNING: mode=%d not implemented "
                "yet GL_LINES=%d GL_LINE_STRIP=%d GL_LINE_LOOP=%d\n",
                mode, GL_LINES, GL_LINE_STRIP, GL_LINE_LOOP);
            break;
          }
          nverts++;
        }
        }
      }
    } break;
    }
  }
  return true;
}

static bool OptimizeVertsToVBONotIndexed(const CGO* I, CGO* cgo,
    const CGOCount& count, float* min, float* max, short* has_draw_buffer,
    bool addshaders)
{
  auto G = I->G;
  bool ok = true;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;

  cgo->alpha = 1.f;
  cgo->color[0] = 1.f;
  cgo->color[1] = 1.f;
  cgo->color[2] = 1.f;

  unsigned mul =
      VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE + VERTEX_ACCESSIBILITY_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal)
              ? 1
              : VERTEX_NORMAL_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color)
              ? 1
              : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_indexes) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * count.num_total_indexes;
  unsigned nxtn = VERTEX_NORMAL_SIZE;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)) {
    normalValsC = (uchar*) normalVals;
    nxtn = 1;
  }
  auto* colorVals = normalVals + nxtn * count.num_total_indexes;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * count.num_total_indexes);
  nxtn = 3;
  auto* accessibilityVals = pickColorVals + nxtn * count.num_total_indexes;

  bool has_normals = false, has_colors = false, has_accessibility = false;
  int ambient_occlusion{};
  ok = CGOProcessCGOtoArrays(I, cgo, cgo, min, max, &ambient_occlusion,
      vertexVals, normalVals, normalValsC, colorVals, colorValsUC,
      pickColorVals, accessibilityVals, has_normals, has_colors,
      has_accessibility);
  if (!ok) {
    if (!G->Interrupt)
      PRINTFB(G, FB_CGO, FB_Errors)
      "ERROR: CGOProcessCGOtoArrays() could not allocate enough "
      "memory\n" ENDFB(G);
    return false;
  }
  if (ok) {
    auto fmt = GetNormalColorFormatSize(G);

    VertexBuffer* vbo =
        G->ShaderMgr->newGPUBuffer<VertexBuffer>(buffer_layout::SEQUENTIAL);
    BufferDataDesc bufData = {{"a_Vertex", VertexFormat::Float3,
        sizeof(float) * count.num_total_indexes * 3, vertexVals}};
    if (has_normals) {
      bufData.push_back({"a_Normal", fmt.normalFormat,
          count.num_total_indexes * fmt.normalSize, normalVals});
    }
    if (has_colors) {
      bufData.push_back({"a_Color", fmt.colorFormat,
          count.num_total_indexes * fmt.colorSize, colorVals});
    }
    if (has_accessibility) {
      bufData.push_back({"a_Accessibility", VertexFormat::Float,
          sizeof(float) * count.num_total_indexes, accessibilityVals});
    }
    ok = vbo->bufferData(std::move(bufData));

    size_t vboid = vbo->get_hash_id();
    // picking VBO: generate a buffer twice the size needed, for each picking
    // pass
    VertexBuffer* pickvbo = G->ShaderMgr->newGPUBuffer<VertexBuffer>(
        buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
    ok = pickvbo->bufferData({BufferDesc{"a_Color", VertexFormat::UByte4Norm,
                                  sizeof(float) * count.num_total_indexes},
        BufferDesc{"a_Color", VertexFormat::UByte4Norm,
            sizeof(float) * count.num_total_indexes}});
    size_t pickvboid = pickvbo->get_hash_id();

    if (ok) {
      float* newPickColorVals;
      int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
                    CGO_PICK_COLOR_ARRAY;
      if (ambient_occlusion) {
        arrays |= CGO_ACCESSIBILITY_ARRAY;
      }
      if (addshaders)
        CGOEnable(cgo, GL_DEFAULT_SHADER_WITH_SETTINGS);
      newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(
          GL_TRIANGLES, arrays, count.num_total_indexes, vboid, pickvboid);
      if (ok && addshaders)
        ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      CHECKOK(ok, newPickColorVals);
      if (!newPickColorVals) {
        G->ShaderMgr->freeGPUBuffer(pickvboid);
        G->ShaderMgr->freeGPUBuffer(vboid);
      }
      if (!ok) {
        PRINTFB(G, FB_CGO, FB_Errors)
        "CGOOptimizeToVBONotIndexed: ERROR: "
        "CGODrawBuffersNotIndexed() could not allocate enough memory\n" ENDFB(
            G);
        return false;
      }
      memcpy(newPickColorVals + count.num_total_indexes, pickColorVals,
          count.num_total_indexes * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      G->ShaderMgr->freeGPUBuffer(vboid);
      G->ShaderMgr->freeGPUBuffer(pickvboid);
    }
  }
  return ok;
}

static bool OptimizeLinesToVBONotIndexed(const CGO* I, CGO* cgo,
    const CGOCount& count, float* min, float* max, short* has_draw_buffer,
    bool addshaders)
{
  auto G = I->G;
  bool ok = true;
  bool has_color = false, has_normals = false;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;

  cgo->alpha = 1.f;
  cgo->color[0] = 1.f;
  cgo->color[1] = 1.f;
  cgo->color[2] = 1.f;

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal)
             ? 1
             : VERTEX_NORMAL_SIZE;
  mul +=
      SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_indexes_lines) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * count.num_total_indexes_lines;
  nxtn = 3;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)) {
    normalValsC = (uchar*) normalVals;
    nxtn = 1;
  }

  auto* colorVals = normalVals + nxtn * count.num_total_indexes_lines;
  if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * count.num_total_indexes_lines);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
    case CGO_SPECIAL:
    case CGO_RESET_NORMAL: {
      const float* newpc = pc;
      cgo->add_to_cgo(op, newpc);
    } break;
    case CGO_NORMAL:
      has_normals = true;
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      break;
    case CGO_COLOR:
      has_color = true;
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      break;
    case CGO_ACCESSIBILITY:
      cgo->current_accessibility = pc[0];
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_ARRAYS: {
      auto sp = it.cast<cgo::draw::arrays>();
      short shouldCompress = false;
      switch (sp->mode) {
      case GL_LINE_LOOP:
      case GL_LINE_STRIP:
      case GL_LINES:
        shouldCompress = true;
      default:
        break;
      }

      if (shouldCompress) {
        int cnt, incr = 0;
        const float* nxtVals = sp->floatdata;
        const float *vertexValsDA = nullptr, *colorValsDA = nullptr,
                    *normalValsDA = nullptr, *pickColorValsDA = nullptr;

        assert(sp->arraybits & CGO_VERTEX_ARRAY);
        vertexValsDA = nxtVals;
        nxtVals += sp->nverts * VERTEX_POS_SIZE;

        for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
          set_min_max(min, max, &vertexValsDA[cnt]);
        }
        if (sp->arraybits & CGO_NORMAL_ARRAY) {
          has_normals = true;
          normalValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_NORMAL_SIZE;
        }

        if (sp->arraybits & CGO_COLOR_ARRAY) {
          has_color = true;
          colorValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_COLOR_SIZE;
        }
        if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
          nxtVals += VERTEX_PICKCOLOR_RGBA_SIZE * sp->nverts;
          pickColorValsDA = nxtVals;
          nxtVals += VERTEX_PICKCOLOR_INDEX_SIZE * sp->nverts;
        }
        float* pickColorValsTMP = pickColorVals + (idxpl * 2);
        switch (sp->mode) {
        case GL_LINES:
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP);
            if (incr && (incr % 2) == 0) {
              FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                  pickColorValsTMP + (incr - 1) * 2);
            }
            idxpl++;
            pl += 3;
            plc += 4;
          }
          break;
        case GL_LINE_STRIP:
          for (cnt = 1; cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
                vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
                vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP);
            idxpl++;
            pl += 3;
            plc += 4;
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP);
            FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                pickColorValsTMP + (incr - 1) * 2);
            idxpl++;
            pl += 3;
            plc += 4;
          }
          break;
        case GL_LINE_LOOP:
          for (cnt = 1; cnt < sp->nverts; cnt++) {
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
                vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
                vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP);
            idxpl++;
            pl += 3;
            plc += 4;
            SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, vertexValsDA,
                normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
                normalValsC, normalVals, colorValsUC, colorVals,
                pickColorValsTMP);
            FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
                pickColorValsTMP + (incr - 1) * 2);
            idxpl++;
            pl += 3;
            plc += 4;
          }
          SetVertexValuesForVBO(G, cgo, pl, plc, 0, incr++, vertexValsDA,
              normalValsDA, colorValsDA, pickColorValsDA, vertexVals,
              normalValsC, normalVals, colorValsUC, colorVals,
              pickColorValsTMP);
          idxpl++;
          pl += 3;
          plc += 4;
          SetVertexValuesForVBO(G, cgo, pl, plc, sp->nverts - 1, incr++,
              vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
              vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
              pickColorValsTMP);
          FixPickColorsForLine(pickColorValsTMP + (incr - 2) * 2,
              pickColorValsTMP + (incr - 1) * 2);
          idxpl++;
          pl += 3;
          plc += 4;
          break;
        }

        //	  pl += 3 * nverts;
        //	  plc += 4 * nverts;
        vpl += sp->nverts;
      }
    } break;
    }
  }
  auto fmt = GetNormalColorFormatSize(G);

  VertexBuffer* vbo =
      G->ShaderMgr->newGPUBuffer<VertexBuffer>(buffer_layout::SEQUENTIAL);
  BufferDataDesc bufData = {{"a_Vertex", VertexFormat::Float3,
      sizeof(float) * count.num_total_indexes_lines * 3, vertexVals}};

  if (has_normals) {
    bufData.push_back({"a_Normal", fmt.normalFormat,
        count.num_total_indexes_lines * fmt.normalSize, normalVals});
  }
  if (has_color) {
    bufData.push_back({"a_Color", fmt.colorFormat,
        count.num_total_indexes_lines * fmt.colorSize, colorVals});
  }
  ok = vbo->bufferData(std::move(bufData));
  size_t vboid = vbo->get_hash_id();

  // picking VBO: generate a buffer twice the size needed, for each picking
  // pass
  VertexBuffer* pickvbo = G->ShaderMgr->newGPUBuffer<VertexBuffer>(
      buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
  ok &= pickvbo->bufferData({BufferDesc("a_Color", VertexFormat::UByte4Norm,
                                 sizeof(float) * count.num_total_indexes_lines),
      BufferDesc("a_Color", VertexFormat::UByte4Norm,
          sizeof(float) * count.num_total_indexes_lines)});
  size_t pickvboid = pickvbo->get_hash_id();

  if (ok) {
    float* newPickColorVals;
    if (addshaders)
      CGOEnable(cgo, GL_DEFAULT_SHADER_WITH_SETTINGS);
    CGODisable(cgo, GL_SHADER_LIGHTING);
    newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(GL_LINES,
        CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
            CGO_PICK_COLOR_ARRAY,
        count.num_total_indexes_lines, vboid, pickvboid);
    if (ok && addshaders)
      ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
    CHECKOK(ok, newPickColorVals);
    if (!ok) {
      PRINTFB(G, FB_CGO, FB_Errors)
      "CGOOptimizeToVBONotIndexed: ERROR: "
      "CGODrawBuffersNotIndexed() could not allocate enough memory\n" ENDFB(G);
      if (!newPickColorVals) {
        G->ShaderMgr->freeGPUBuffer(pickvboid);
        G->ShaderMgr->freeGPUBuffer(vboid);
      }
      return false;
    }
    memcpy(newPickColorVals + count.num_total_indexes_lines, pickColorVals,
        count.num_total_indexes_lines * 2 * sizeof(float));
    *has_draw_buffer = true;
  } else {
    G->ShaderMgr->freeGPUBuffer(pickvboid);
    G->ShaderMgr->freeGPUBuffer(vboid);
  }
  return ok;
}

bool CGOOptimizeToVBONotIndexed(CGO** I)
{
  CGO* cgo = CGOOptimizeToVBONotIndexed(*I, 0, true);
  CGOFree(*I);
  *I = cgo;
  return (cgo != nullptr);
}

/**
 * Converts primitive operations into VBO operations.
 *
 * If the input CGO has sortable alpha triangles, you should use
 * CGOOptimizeToVBOIndexed() instead.
 *
 * @param I Primitive CGO to convert
 * @param est Output CGO size estimate (size of buffer to "reserve")
 * @param addshaders Add Enable/Disable shader operations
 *
 * @return New converted CGO with use_shader=true
 */
CGO* CGOOptimizeToVBONotIndexed(const CGO* I, int est, bool addshaders)
{
  auto G = I->G;

  std::unique_ptr<CGO> I_begin_end_combined;
  if (I->has_begin_end) {
    I_begin_end_combined.reset(CGOCombineBeginEnd(I));
    I = I_begin_end_combined.get();
    assert(!I->has_begin_end);
  }

  short has_draw_buffer = false;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX},
        max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  int ambient_occlusion = 0;
  int ok = true;
  auto cgo = CGONew(G);

  auto count = CGOCountNumVertices(I);
  if (count.num_total_vertices_points > 0) {
    if (!OptimizePointsToVBO(I, cgo, count.num_total_vertices_points, min, max,
            &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (count.num_total_indexes > 0) {
    if (!OptimizeVertsToVBONotIndexed(
            I, cgo, count, min, max, &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (ok && count.num_total_indexes_lines > 0) {
    if (!OptimizeLinesToVBONotIndexed(
            I, cgo, count, min, max, &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }

  if (ok && (count.num_total_vertices > 0 || count.num_total_vertices_lines > 0 ||
                count.num_total_vertices_points > 0)) {
    ok &= CGOBoundingBox(cgo, min, max);
  }

  if (ok)
    ok &= CGOStop(cgo);

  if (!ok) {
    CGOFree(cgo);
    return nullptr;
  }

  cgo->has_draw_buffers |= has_draw_buffer;
  cgo->use_shader = true;
  cgo->cgo_shader_ub_color = SettingGet<int>(G, cSetting_cgo_shader_ub_color);
  cgo->cgo_shader_ub_normal = SettingGet<int>(G, cSetting_cgo_shader_ub_normal);
  return (cgo);
}

static bool OptimizeVertsToVBOIndexed(const CGO* I, CGO* cgo,
    const CGOCount& count, float* min, float* max, short* has_draw_buffer,
    bool addshaders, bool embedTransparencyInfo)
{
  auto G = I->G;
  bool ok = true;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;
  short ambient_occlusion = 0;
  float* sumarray = nullptr;
  int n_data = 0;

  if (embedTransparencyInfo) {
    int n_tri = count.num_total_indexes / 3;
    int bytes_to_allocate =
        2 * count.num_total_indexes *
            sizeof(VertexIndex_t) + // vertexIndicesOriginal, vertexIndices
        3 * count.num_total_indexes * sizeof(float) + // 3 * for sum
        n_tri * sizeof(float) +
        2 * n_tri * sizeof(int) +
        256 * sizeof(int); // z_value (float * n_tri), ix (n_tri * int),
                            // sort_mem ((n_tri + 256) * int)
    // round to 4 byte words for the length of the CGO
    n_data = bytes_to_allocate / 4 + (((bytes_to_allocate % 4) == 0) ? 0 : 1);
  }
  std::vector<VertexIndex_t> vertexIndices(count.num_total_indexes);

  unsigned mul =
      VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE + VERTEX_ACCESSIBILITY_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal)
              ? 1
              : VERTEX_NORMAL_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color)
              ? 1
              : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_vertices) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* normalVals = vertexVals + 3 * count.num_total_vertices;
  nxtn = 3;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)) {
    normalValsC = (uchar*) normalVals;
    nxtn = 1;
  }

  auto* colorVals = normalVals + nxtn * count.num_total_vertices;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  auto* pickColorVals = (colorVals + nxtn * count.num_total_vertices);
  auto* accessibilityVals = pickColorVals + 3 * count.num_total_vertices;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
    case CGO_NORMAL:
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_ACCESSIBILITY:
      cgo->current_accessibility = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_ARRAYS: {
      auto sp = it.cast<cgo::draw::arrays>();
      short shouldCompress = false;
      switch (sp->mode) {
      case GL_TRIANGLE_FAN:
      case GL_TRIANGLE_STRIP:
      case GL_TRIANGLES:
        shouldCompress = true;
      default:
        break;
      }
      if (shouldCompress) {
        int cnt, nxtn = 3;
        float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0,
              *normalValsDA, *accessibilityValsDA;
        float *pickColorValsDA, *pickColorValsTMP, *accessibilityValsTMP;

        nxtVals = vertexValsDA = sp->floatdata;
        for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
          set_min_max(min, max, &vertexValsDA[cnt]);
        }
        for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
          vertexVals[pl + cnt] = vertexValsDA[cnt];
        }
        if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)) {
          if (sp->arraybits & CGO_NORMAL_ARRAY) {
            nxtVals = normalValsDA = vertexValsDA + (nxtn * sp->nverts);
            for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
              normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] =
                  CLIP_NORMAL_VALUE(normalValsDA[cnt]);
            }
          } else {
            uchar norm[3] = {CLIP_NORMAL_VALUE(cgo->normal[0]),
                CLIP_NORMAL_VALUE(cgo->normal[1]),
                CLIP_NORMAL_VALUE(cgo->normal[2])};
            for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
              normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] =
                  norm[cnt % 3];
            }
          }
        } else {
          if (sp->arraybits & CGO_NORMAL_ARRAY) {
            nxtVals = normalValsDA = vertexValsDA + (nxtn * sp->nverts);
            for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
              normalVals[pl + cnt] = normalValsDA[cnt];
            }
          } else {
            for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
              normalVals[pl + cnt] = cgo->normal[cnt % 3];
            }
          }
        }
        nxtn = 3;
        if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
          if (sp->arraybits & CGO_COLOR_ARRAY) {
            nxtVals = colorValsDA = nxtVals + (nxtn * sp->nverts);
            for (cnt = 0; cnt < sp->nverts * 4; cnt += 4) {
              colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
              colorValsUC[plc + cnt + 1] =
                  CLIP_COLOR_VALUE(colorValsDA[cnt + 1]);
              colorValsUC[plc + cnt + 2] =
                  CLIP_COLOR_VALUE(colorValsDA[cnt + 2]);
              colorValsUC[plc + cnt + 3] =
                  CLIP_COLOR_VALUE(colorValsDA[cnt + 3]);
            }
            nxtn = VERTEX_COLOR_SIZE;
          } else {
            uchar col[4] = {CLIP_COLOR_VALUE(cgo->color[0]),
                CLIP_COLOR_VALUE(cgo->color[1]),
                CLIP_COLOR_VALUE(cgo->color[2]),
                CLIP_COLOR_VALUE(cgo->alpha)};
            for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
              colorValsUC[plc + cnt] = col[cnt % 4];
            }
          }
        } else {
          if (sp->arraybits & CGO_COLOR_ARRAY) {
            nxtVals = colorValsDA = nxtVals + (nxtn * sp->nverts);
            for (cnt = 0; cnt < sp->nverts * 4; cnt += 4) {
              colorVals[plc + cnt] = colorValsDA[cnt];
              colorVals[plc + cnt + 1] = colorValsDA[cnt + 1];
              colorVals[plc + cnt + 2] = colorValsDA[cnt + 2];
              colorVals[plc + cnt + 3] = colorValsDA[cnt + 3];
            }
            nxtn = VERTEX_COLOR_SIZE;
          } else {
            float col[4] = {
                cgo->color[0], cgo->color[1], cgo->color[2], cgo->alpha};
            for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
              colorVals[plc + cnt] = col[cnt % 4];
            }
          }
        }
        if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
          nxtVals = nxtVals + (nxtn * sp->nverts);
          pickColorValsDA = nxtVals + sp->nverts;
          pickColorValsTMP = pickColorVals + (vpl * 2);
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
            CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
          }
          nxtn = VERTEX_PICKCOLOR_SIZE;
        } else {
          pickColorValsTMP = pickColorVals + (vpl * 2);
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            CGO_put_uint(pickColorValsTMP++, cgo->current_pick_color_index);
            CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
          }
        }
        if (sp->arraybits & CGO_ACCESSIBILITY_ARRAY) {
          if (!ambient_occlusion) {
            for (cnt = 0; cnt < vpl; cnt++) {
              accessibilityVals[cnt] = 1.f;
            }
          }
          ambient_occlusion = 1;
          nxtVals = nxtVals + (nxtn * sp->nverts);
          accessibilityValsDA = nxtVals;
          accessibilityValsTMP = accessibilityVals + vpl;
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            accessibilityValsTMP[cnt] = accessibilityValsDA[cnt];
          }
        } else {
          if (ambient_occlusion) {
            accessibilityValsTMP = accessibilityVals + vpl;
            for (cnt = 0; cnt < sp->nverts; cnt++) {
              accessibilityValsTMP[cnt] = 1.f;
            }
          }
        }
        switch (sp->mode) {
        case GL_TRIANGLES:
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            vertexIndices[idxpl++] = vpl + cnt;
          }
          break;
        case GL_TRIANGLE_STRIP: {
          short flip = 0;
          for (cnt = 2; cnt < sp->nverts; cnt++) {
            vertexIndices[idxpl++] = vpl + cnt - (flip ? 0 : 2);
            vertexIndices[idxpl++] = vpl + cnt - 1;
            vertexIndices[idxpl++] = vpl + cnt - (flip ? 2 : 0);
            flip = !flip;
          }
        } break;
        case GL_TRIANGLE_FAN:
          for (cnt = 2; cnt < sp->nverts; cnt++) {
            vertexIndices[idxpl++] = vpl;
            vertexIndices[idxpl++] = vpl + cnt - 1;
            vertexIndices[idxpl++] = vpl + cnt;
          }
          break;
        }
        pl += 3 * sp->nverts;
        plc += 4 * sp->nverts;
        vpl += sp->nverts;
      }
    } break;
    default:
      break;
    }
    ok &= !I->G->Interrupt;
  }
  if (sumarray) {
    for (idxpl = 0; idxpl < count.num_total_indexes; idxpl += 3) {
      add3f(&vertexVals[3 * vertexIndices[idxpl]],
          &vertexVals[3 * vertexIndices[idxpl + 1]], sumarray);
      add3f(&vertexVals[3 * vertexIndices[idxpl + 2]], sumarray, sumarray);
      sumarray += 3;
    }
  }
  if (ok) {
    auto fmt = GetNormalColorFormatSize(I->G);

    VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
    ok &= vbo->bufferData(
        {BufferDesc{"a_Vertex", VertexFormat::Float3,
              sizeof(float) * count.num_total_vertices * 3, vertexVals},
            BufferDesc{"a_Normal", fmt.normalFormat,
                count.num_total_vertices * fmt.normalSize, normalVals},
            BufferDesc{"a_Color", fmt.colorFormat,
                count.num_total_vertices * fmt.colorSize, colorVals},
            BufferDesc{"a_Accessibility", VertexFormat::Float,
                sizeof(float) * count.num_total_vertices, accessibilityVals}});

    IndexBuffer* ibo = I->G->ShaderMgr->newGPUBuffer<IndexBuffer>();
    ok &= ibo->bufferData({BufferDesc{nullptr, VertexFormat::UInt,
        sizeof(VertexIndex_t) * count.num_total_indexes, vertexIndices.data()}});

    size_t vboid = vbo->get_hash_id();
    size_t iboid = ibo->get_hash_id();

    VertexBuffer* pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
        buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
    ok &= pickvbo->bufferData({BufferDesc{"a_Color", VertexFormat::UByte4Norm,
                                    sizeof(float) * count.num_total_indexes},
        BufferDesc{"a_Color", VertexFormat::UByte4Norm,
            sizeof(float) * count.num_total_indexes}});
    size_t pickvboid = pickvbo->get_hash_id();

    if (ok) {
      float* newPickColorVals;
      int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
                    CGO_PICK_COLOR_ARRAY;
      if (ambient_occlusion) {
        arrays |= CGO_ACCESSIBILITY_ARRAY;
      }
      if (addshaders)
        CGOEnable(cgo, GL_DEFAULT_SHADER);
      newPickColorVals = cgo->add<cgo::draw::buffers_indexed>(GL_TRIANGLES,
          arrays, count.num_total_indexes, count.num_total_vertices, vboid, iboid, n_data,
          pickvboid);
      if (embedTransparencyInfo) {
        int n_tri = count.num_total_indexes / 3;
        float* sumarray;
        float* sum = sumarray = newPickColorVals + count.num_total_vertices * 3;
        float* z_value = sum + (count.num_total_indexes * 3);
        int* ix = (int*) z_value + n_tri;
        int* sort_mem = ix + n_tri;
        auto vertexIndicesOriginalTI =
            (VertexIndex_t*) (sort_mem + n_tri + 256);

        for (idxpl = 0; idxpl < count.num_total_indexes; idxpl += 3) {
          add3f(&vertexVals[3 * vertexIndices[idxpl]],
              &vertexVals[3 * vertexIndices[idxpl + 1]], sumarray);
          add3f(
              &vertexVals[3 * vertexIndices[idxpl + 2]], sumarray, sumarray);
          sumarray += 3;
        }
        memcpy(vertexIndicesOriginalTI, vertexIndices.data(),
            sizeof(VertexIndex_t) * count.num_total_indexes);
      }

      if (addshaders && ok)
        ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      CHECKOK(ok, newPickColorVals);
      if (!newPickColorVals) {
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(iboid);
      }
      if (ok)
        memcpy(newPickColorVals + count.num_total_vertices, pickColorVals,
            count.num_total_vertices * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(pickvboid);
      I->G->ShaderMgr->freeGPUBuffer(vboid);
      I->G->ShaderMgr->freeGPUBuffer(iboid);
    }
  }
  return ok;
}

static bool OptimizeLinesToVBOIndexed(const CGO* I, CGO* cgo,
    const CGOCount& count, float* min, float* max, short* has_draw_buffer,
    bool addshaders)
{
  auto G = I->G;
  bool ok = true;
  uchar* colorValsUC = 0;
  uchar* normalValsC = 0;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, sz;
  bool hasNormals = 0;

  hasNormals = !CGOHasAnyLineVerticesWithoutNormals(I);
  std::vector<VertexIndex_t> vertexIndexes(count.num_total_indexes_lines);

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  if (hasNormals) {
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal)
                ? 1
                : VERTEX_NORMAL_SIZE;
  }
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color)
              ? 1
              : VERTEX_COLOR_SIZE;
  auto const tot = size_t(count.num_total_vertices_lines) * mul;

  std::vector<float> vertexValsVec(tot);
  auto* vertexVals = vertexValsVec.data();
  auto* nxtVals = vertexVals + VERTEX_POS_SIZE * count.num_total_vertices_lines;
  float* normalVals = nullptr;
  if (hasNormals) {
    normalVals = nxtVals;
    if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)) {
      normalValsC = (uchar*) normalVals;
      sz = 1;
    } else {
      sz = VERTEX_NORMAL_SIZE;
    }
    nxtVals += sz * count.num_total_vertices_lines;
  }

  auto* colorVals = nxtVals;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
    colorValsUC = (uchar*) colorVals;
    sz = 1;
  } else {
    sz = VERTEX_COLOR_SIZE;
  }
  nxtVals += sz * count.num_total_vertices_lines;

  auto* pickColorVals = nxtVals;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
    case CGO_NORMAL:
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      break;
    case CGO_ACCESSIBILITY:
      cgo->current_accessibility = *pc;
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_ARRAYS: {
      auto sp = it.cast<cgo::draw::arrays>();
      short shouldCompress = false;
      switch (sp->mode) {
      case GL_LINE_LOOP:
      case GL_LINE_STRIP:
      case GL_LINES:
        shouldCompress = true;
      default:
        break;
      }
      if (shouldCompress) {
        int cnt, nxtn = 3;
        float *vertexValsDA = 0, *nxtVals2 = 0, *colorValsDA = 0,
              *normalValsDA = 0;
        float *pickColorValsDA = 0, *pickColorValsTMP;

        nxtVals2 = vertexValsDA = sp->floatdata;
        for (cnt = 0; cnt < sp->nverts * 3; cnt += 3) {
          set_min_max(min, max, &vertexValsDA[cnt]);
        }
        for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
          vertexVals[pl + cnt] = vertexValsDA[cnt];
        }
        if (normalVals) {
          if (sp->arraybits & CGO_NORMAL_ARRAY) {
            if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)) {
              nxtVals2 = normalValsDA = nxtVals2 + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] =
                    CLIP_NORMAL_VALUE(normalValsDA[cnt]);
              }
            } else {
              nxtVals2 = normalValsDA = nxtVals2 + (nxtn * sp->nverts);
              for (cnt = 0; cnt < sp->nverts * 3; cnt++) {
                normalVals[VAR_FOR_NORMAL + cnt] = normalValsDA[cnt];
              }
            }
          }
        }
        if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)) {
          if (sp->arraybits & CGO_COLOR_ARRAY) {
            nxtVals2 = colorValsDA = nxtVals2 + (nxtn * sp->nverts);
            for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
              colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
            }
            nxtn = 4;
          } else {
            uchar col[4] = {CLIP_COLOR_VALUE(cgo->color[0]),
                CLIP_COLOR_VALUE(cgo->color[1]),
                CLIP_COLOR_VALUE(cgo->color[2]),
                CLIP_COLOR_VALUE(cgo->alpha)};
            for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
              colorValsUC[plc + cnt] = col[cnt % 4];
            }
          }
        } else {
          if (sp->arraybits & CGO_COLOR_ARRAY) {
            nxtVals2 = colorValsDA = nxtVals2 + (nxtn * sp->nverts);
            for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
              colorVals[plc + cnt] = colorValsDA[cnt];
            }
            nxtn = 4;
          } else {
            float col[4] = {
                cgo->color[0], cgo->color[1], cgo->color[2], cgo->alpha};
            for (cnt = 0; cnt < sp->nverts * 4; cnt++) {
              colorVals[plc + cnt] = col[cnt % 4];
            }
          }
        }
        if (sp->arraybits & CGO_PICK_COLOR_ARRAY) {
          nxtVals2 = nxtVals2 + (nxtn * sp->nverts);
          pickColorValsDA = nxtVals2 + sp->nverts;
          pickColorValsTMP = pickColorVals + (vpl * 2);
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
            CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
          }
          nxtn = VERTEX_PICKCOLOR_SIZE;
        } else {
          pickColorValsTMP = pickColorVals + (vpl * 2);
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            CGO_put_uint(pickColorValsTMP++, cgo->current_pick_color_index);
            CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
          }
        }
        if (idxpl + sp->nverts > count.num_total_indexes_lines) {
          PRINTFB(I->G, FB_CGO, FB_Errors)
          "ERROR: CGOOptimizeToVBOIndexed() num_total_indexes_lines=%d "
          "mode=%d nverts=%d idxpl=%d\n",
              count.num_total_indexes_lines, sp->mode, sp->nverts,
              idxpl ENDFB(I->G);
        }
        switch (sp->mode) {
        case GL_LINES:
          for (cnt = 0; cnt < sp->nverts; cnt++) {
            vertexIndexes[idxpl++] = vpl + cnt;
          }
          break;
        case GL_LINE_STRIP:
          for (cnt = 1; cnt < sp->nverts; cnt++) {
            vertexIndexes[idxpl++] = vpl + cnt - 1;
            vertexIndexes[idxpl++] = vpl + cnt;
          }
          break;
        case GL_LINE_LOOP:
          for (cnt = 1; cnt < sp->nverts; cnt++) {
            vertexIndexes[idxpl++] = vpl + cnt - 1;
            vertexIndexes[idxpl++] = vpl + cnt;
          }
          vertexIndexes[idxpl++] = vpl;
          vertexIndexes[idxpl++] = vpl + sp->nverts - 1;
          break;
        }

        pl += 3 * sp->nverts;
        plc += 4 * sp->nverts;
        vpl += sp->nverts;
      }
    } break;
    case CGO_SPECIAL:
      CGOSpecial(cgo, CGO_get_int(pc));
    }
    ok &= !I->G->Interrupt;
  }
  if (ok) {
    auto fmt = GetNormalColorFormatSize(I->G);

    VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
    ok &= vbo->bufferData({
        BufferDesc{"a_Vertex", VertexFormat::Float3,
            sizeof(float) * count.num_total_vertices_lines * 3, vertexVals},
        BufferDesc{"a_Normal", fmt.normalFormat,
            count.num_total_vertices_lines * fmt.normalSize, normalVals},
        BufferDesc{"a_Color", fmt.colorFormat,
            count.num_total_vertices_lines * fmt.colorSize, colorVals},
    });

    IndexBuffer* ibo = I->G->ShaderMgr->newGPUBuffer<IndexBuffer>();
    ok &= ibo->bufferData({BufferDesc(nullptr, VertexFormat::UInt,
        sizeof(VertexIndex_t) * count.num_total_indexes_lines,
        vertexIndexes.data())});

    size_t vboid = vbo->get_hash_id();
    size_t iboid = ibo->get_hash_id();

    VertexBuffer* pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
        buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
    ok &= pickvbo->bufferData({BufferDesc("a_Color", VertexFormat::UByte4Norm,
                                    sizeof(float) * count.num_total_indexes),
        BufferDesc("a_Color", VertexFormat::UByte4Norm,
            sizeof(float) * count.num_total_indexes)});
    size_t pickvboid = pickvbo->get_hash_id();

    if (ok) {
      if (addshaders) {
        CGOEnable(cgo, GL_DEFAULT_SHADER);
        // TODO: Check if this next line is supposed to be in this
        // if-statement. VBONotIndexed has it otherwise.
        CGODisable(cgo, GL_SHADER_LIGHTING);
      }
      auto newPickColorVals = cgo->add<cgo::draw::buffers_indexed>(GL_LINES,
          CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY |
              CGO_PICK_COLOR_ARRAY,
          count.num_total_indexes_lines, count.num_total_vertices_lines, vboid,
          iboid, 0, pickvboid);
      CHECKOK(ok, newPickColorVals);
      if (addshaders && ok)
        ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      if (!newPickColorVals) {
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(iboid);
      }
      if (ok)
        memcpy(newPickColorVals + count.num_total_vertices_lines, pickColorVals,
            count.num_total_vertices_lines * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(pickvboid);
      I->G->ShaderMgr->freeGPUBuffer(vboid);
      I->G->ShaderMgr->freeGPUBuffer(iboid);
    }
  }
  return ok;
}

/**
 * Similar to CGOOptimizeToVBONotIndexed(), but works with sortable alpha
 * triangles.
 *
 * @param I Primitive CGO to convert
 * @param est Output CGO size estimate (size of buffer to "reserve")
 * @param color Initial color to use until a CGO_COLOR operation is observed
 * @param addshaders Add Enable/Disable shader operations
 * @param embedTransparencyInfo ???
 *
 * @return New converted CGO with use_shader=true
 */
CGO* CGOOptimizeToVBOIndexed(const CGO* I, int est, const float* color,
    bool addshaders, bool embedTransparencyInfo)
{
  auto G = I->G;
  CGO* cgo;

  std::unique_ptr<CGO> I_begin_end_combined;
  if (I->has_begin_end) {
    I_begin_end_combined.reset(CGOCombineBeginEnd(I));
    I = I_begin_end_combined.get();
    assert(!I->has_begin_end);
  }

  short has_draw_buffer = false;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX},
        max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  int ok = true;

  auto count = CGOCountNumVertices(I);

  cgo = CGONewSized(I->G, I->c + est);
  CHECKOK(ok, cgo);
  if (ok) {
    if (color) {
      cgo->color[0] = color[0];
      cgo->color[1] = color[1];
      cgo->color[2] = color[2];
      cgo->alpha = color[3];
    } else {
      cgo->color[0] = 1.f;
      cgo->color[1] = 1.f;
      cgo->color[2] = 1.f;
      cgo->alpha = 1.f;
    }
  }

  if (count.num_total_vertices_points > 0) {
    /* This does not need to be indexed (for now) */
    if (!OptimizePointsToVBO(I, cgo, count.num_total_vertices_points, min, max,
            &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }

  if (count.num_total_vertices > 0) {
    if (!OptimizeVertsToVBOIndexed(I, cgo, count, min, max, &has_draw_buffer,
            addshaders, embedTransparencyInfo)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (ok && count.num_total_vertices_lines > 0) {
    if (!OptimizeLinesToVBOIndexed(
            I, cgo, count, min, max, &has_draw_buffer, addshaders)) {
      CGOFree(cgo);
      return nullptr;
    }
  }
  if (ok && (count.num_total_vertices > 0 || count.num_total_vertices_lines > 0)) {
    ok &= CGOBoundingBox(cgo, min, max);
  }

  if (ok)
    ok &= CGOStop(cgo);
  if (ok) {
    if (has_draw_buffer) {
      cgo->has_draw_buffers = true;
    }
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader) {
      cgo->cgo_shader_ub_color =
          SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
      cgo->cgo_shader_ub_normal =
          SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

struct OptimizeSphereData {
  std::vector<float> vert;
  std::vector<unsigned char> color;
  std::vector<int> pickColor;
  std::vector<unsigned char> rightUpFlagsUB;
  std::vector<float> rightUpFlags;
  float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
  float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  int total_vert = 0;
  int total_spheres = 0;
};

/**
 * @param[in] I input CGO
 * @param[out] cgo output CGO
 * @param[in] num_total_spheres Number of spheres to optimize
 * @param[out] leftOverCGO CGO data not relevant to spheres.
 * @pre cgo must not be nullptr
 *
 * Note: To use leftOverCGO, it must already be preallocated.
 */

static OptimizeSphereData GetOptimizeSphereData(
    const CGO* I, CGO*& cgo, int num_total_spheres, CGO* leftOverCGO)
{
  OptimizeSphereData sphereData;
  const int tot = VerticesPerSphere() * 4 * num_total_spheres;
#if defined(PURE_OPENGL_ES_2)
  int rightup_flags[6] = {0, 1, 3, 3, 2, 0};
#else
  int rightup_flags[4] = {0, 1, 3, 2};
#endif

  auto cgo_shader_ub_flags =
      SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_flags);

  sphereData.vert.resize(tot);
  auto vertVals = sphereData.vert.data();

  sphereData.color.resize(tot);
  auto colorValsUB = sphereData.color.data();

  unsigned char* rightUpFlagValsUB = nullptr;
  float* rightUpFlagVals = nullptr;
  if (cgo_shader_ub_flags) {
    sphereData.rightUpFlagsUB.resize(VALUES_PER_IMPOSTER_SPACE_COORD *
                                     VerticesPerSphere() * num_total_spheres);
    rightUpFlagValsUB = sphereData.rightUpFlagsUB.data();
  } else {
    sphereData.rightUpFlags.resize(VALUES_PER_IMPOSTER_SPACE_COORD *
                                   VerticesPerSphere() * num_total_spheres);
    rightUpFlagVals = sphereData.rightUpFlags.data();
  }

  bool has_picking = CGOHasOperationsOfType(I, CGO_PICK_COLOR);
  int* pickcolorVals = nullptr;
  if (has_picking) {
    // atom/bond info for picking, 2 ints for each sphere
    sphereData.pickColor.resize(num_total_spheres * 2 * 4);
    pickcolorVals = sphereData.pickColor.data();
  }

  cgo->alpha = 1.f;
  float min_alpha = 1.f;
  bool copyNormalToLeftOver = false;
  bool copyColorToLeftOver = false;
  bool copyPickColorToLeftOver = false;
  bool copyAlphaToLeftOver = false;
  bool ok = true;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_NORMAL:
      cgo->normal[0] = *pc;
      cgo->normal[1] = *(pc + 1);
      cgo->normal[2] = *(pc + 2);
      copyNormalToLeftOver = true;
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc;
      cgo->color[1] = *(pc + 1);
      cgo->color[2] = *(pc + 2);
      copyColorToLeftOver = true;
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      if (cgo->alpha < min_alpha)
        min_alpha = cgo->alpha;
      copyAlphaToLeftOver = true;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      copyPickColorToLeftOver = true;
      break;
    case CGO_SPHERE:
      for (int vv = 0; vv < VerticesPerSphere();
           vv++) { // generate eight vertices of a bounding box for each
                   // cylinder
        vertVals[0] = *(pc);
        vertVals[1] = *(pc + 1);
        vertVals[2] = *(pc + 2);
        vertVals[3] = *(pc + 3);
        set_min_max(sphereData.min, sphereData.max, vertVals);
        if (cgo_shader_ub_flags) {
          rightUpFlagValsUB[0] = rightup_flags[vv];
          rightUpFlagValsUB++;
        } else {
          rightUpFlagVals[0] = rightup_flags[vv];
          rightUpFlagVals++;
        }
        colorValsUB[0] = CLIP_COLOR_VALUE(cgo->color[0]);
        colorValsUB[1] = CLIP_COLOR_VALUE(cgo->color[1]);
        colorValsUB[2] = CLIP_COLOR_VALUE(cgo->color[2]);
        colorValsUB[3] = CLIP_COLOR_VALUE(cgo->alpha);
        colorValsUB += 4;
        vertVals += 4;
        sphereData.total_vert++;
      }
      if (has_picking) {
        *(pickcolorVals++) = cgo->current_pick_color_index;
        *(pickcolorVals++) = cgo->current_pick_color_bond;
      }
      sphereData.total_spheres++;
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      "WARNING: CGOOptimizeSpheresToVBONonIndexed() CGO_DRAW_BUFFERS_INDEXED "
      "or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
          op ENDFB(I->G);
      break;
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      "WARNING: CGOOptimizeCylindersToVBO() "
      "CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS encountered op=0x%X\n",
          op ENDFB(I->G);
      break;
    case CGO_DRAW_LABELS:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_LABELS encountered "
      "op=0x%X\n",
          op ENDFB(I->G);
      break;
    case CGO_DRAW_TEXTURES:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_TEXTURES encountered "
      "op=0x%X\n",
          op ENDFB(I->G);
      break;
    case CGO_DRAW_ARRAYS:
    default:
      if (!leftOverCGO)
        break;
      if (copyAlphaToLeftOver) {
        copyAlphaToLeftOver = false;
        CGOAlpha(leftOverCGO, cgo->alpha);
      }
      if (copyColorToLeftOver) {
        copyColorToLeftOver = false;
        CGOColor(leftOverCGO, cgo->color[0], cgo->color[1], cgo->color[2]);
      }
      if (copyNormalToLeftOver) {
        CGONormalv(leftOverCGO, cgo->normal);
      }
      if (copyPickColorToLeftOver) {
        copyPickColorToLeftOver = false;
        CGOPickColor(leftOverCGO, cgo->current_pick_color_index,
            cgo->current_pick_color_bond);
      }
      leftOverCGO->add_to_cgo(op, pc);
    }
#ifndef _WEBGL
    ok &= !I->G->Interrupt;
#endif
  }
  return sphereData;
}

/**
 * Creates and fills OpenGL buffers with optimized sphere Data
 * @param[in] I input CGO
 * @param[out] cgo output optimized CGO
 * @param[in] addshaders adds shader call to CGO
 * @param[in] num_total_spheres number of total spheres to optimize
 * @param[in] sphereData Optimized Sphere Data
 * @pre cgo must not be nullptr
 */

static bool PopulateGLBufferOptimizedSphereData(const CGO* I, CGO* cgo,
    bool addshaders, int num_total_spheres,
    const OptimizeSphereData& sphereData)
{
  bool ok = true;
  auto rtp = VertexFormat::Float;
  short rsz = sizeof(float);
  const void* radiusptr = sphereData.rightUpFlags.data();
  auto cgo_shader_ub_flags =
      SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_flags);
  if (cgo_shader_ub_flags) {
    rtp = VertexFormat::UByte;
    rsz = sizeof(std::uint8_t);
    radiusptr = sphereData.rightUpFlagsUB.data();
  }

  VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
  ok &= vbo->bufferData({
      BufferDesc("a_vertex_radius", VertexFormat::Float4,
          sizeof(float) * sphereData.total_vert * 4, sphereData.vert.data()),
      BufferDesc("a_Color", VertexFormat::UByte4Norm,
          sizeof(float) * sphereData.total_vert, sphereData.color.data()),
      BufferDesc("a_rightUpFlags", rtp,
          rsz * sphereData.total_vert * VALUES_PER_IMPOSTER_SPACE_COORD,
          radiusptr),
  });
  size_t vboid = vbo->get_hash_id();

  VertexBuffer* pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
      buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
  ok &= pickvbo->bufferData({BufferDesc("a_Color", VertexFormat::UByte4Norm, 0),
                                BufferDesc("a_Color", VertexFormat::UByte4Norm,
                                    sizeof(float) * sphereData.total_vert)},
      0, sizeof(float) * sphereData.total_vert * 2, 0);
  size_t pickvboid = pickvbo->get_hash_id();

  cgo->has_draw_buffers = true;
  cgo->has_draw_sphere_buffers = true;

  auto freebuffers = [vboid, pickvboid, I]() {
    I->G->ShaderMgr->freeGPUBuffer(vboid);
    I->G->ShaderMgr->freeGPUBuffer(pickvboid);
  };
  if (ok) {
    if (addshaders)
      CGOEnable(cgo, GL_SPHERE_SHADER);
    auto pickcolor_data = (int*) cgo->add<cgo::draw::sphere_buffers>(
        sphereData.total_spheres, (cgo_shader_ub_flags ? 3 : 1), vboid,
        pickvboid); // always cgo_shader_ub_color
    CHECKOK(ok, pickcolor_data);
    if (ok && !sphereData.pickColor.empty()) {
      memcpy(pickcolor_data, sphereData.pickColor.data(),
          num_total_spheres * 2 * 4);
    }
    if (ok && addshaders)
      ok &= CGODisable(cgo, GL_SPHERE_SHADER);
    if (!ok) {
      freebuffers();
    }
  } else {
    freebuffers();
  }
  return ok;
}

/**
 * Optimizes Sphere CGO and populates data into GPU buffers
 * @param[in] I input CGO
 * @param[in] est estimated size of output CGO
 * @param[in] addshaders adds shader call to CGO
 * @param[out] leftOverCGO CGO data not relevant to spheres.
 * @return optimized CGO
 *
 * Note: If estimated size is unknown beforehand, provide 0 to est.
 * Note: To use leftOverCGO, it must already be preallocated.
 */

CGO* CGOOptimizeSpheresToVBONonIndexed(
    const CGO* I, int est, bool addshaders, CGO* leftOverCGO)
{
  bool ok = true;
  int num_total_spheres = CGOCountNumberOfOperationsOfType(I, CGO_SPHERE);

  if (num_total_spheres <= 0) {
    return nullptr;
  }
  auto cgo = CGONewSized(I->G, I->c + est);
  auto sphereData =
      GetOptimizeSphereData(I, cgo, num_total_spheres, leftOverCGO);

  if (ok && sphereData.total_spheres > 0) {
    ok = PopulateGLBufferOptimizedSphereData(
        I, cgo, addshaders, num_total_spheres, sphereData);
  }

  if (ok && num_total_spheres > 0) {
    ok &= CGOBoundingBox(cgo, sphereData.min, sphereData.max);
  }

  if (ok)
    ok &= CGOStop(cgo);

  if (ok) {
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader) {
      cgo->cgo_shader_ub_color = true;
      cgo->cgo_shader_ub_normal =
          SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_normal);
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

CGO* CGOOptimizeBezier(const CGO* I)
{
  auto cgo = std::make_unique<CGO>(I->G);
  int num_splines = CGOCountNumberOfOperationsOfType(I, CGO_BEZIER);
  auto vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
  std::vector<float> vertData;
  vertData.reserve(num_splines * CGO_BEZIER_SZ);
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_BEZIER: {
      vertData.resize(vertData.size() + CGO_BEZIER_SZ);
      std::copy_n(pc, CGO_BEZIER_SZ, vertData.end() - CGO_BEZIER_SZ);
    } break;
    }
  }

  std::size_t numDimensions = 3;
  std::size_t numVerts = 4;
  vbo->bufferData({
      BufferDesc("position", VertexFormat::Float3,
          sizeof(float) * numVerts * numDimensions, vertData.data()),
  });
  size_t vboid = vbo->get_hash_id();

  CGOEnable(cgo.get(), GL_BEZIER_SHADER);
  cgo->add<cgo::draw::bezier_buffers>(vboid);
  CGODisable(cgo.get(), GL_BEZIER_SHADER);
  cgo->use_shader = true;
  return cgo.release();
}

/**
 * converts a CGO that has primitives into pure geometry,
 *    and converts CGO_BEGIN/CGO_END blocks into CGO_DRAW_ARRAYS
 *    operations, similar to what CGOCombineBeginEnd() does.
 *
 * I:               input CGO
 * est:             initial size of the newly allocated CGO array
 *                     that is returned by this function
 * sphere_quality:  the quality of the spheres generated by this function
 *                     (if -1, defaults to cgo_sphere_quality)
 * stick_round_nub: if true, a round cap is generated, otherwise, it generates
 *                  the old "pointed" caps
 */
CGO* CGOSimplify(
    const CGO* I, int est, short sphere_quality, bool stick_round_nub)
{
  auto G = I->G;
  int ok = true;
  if (sphere_quality < 0) {
    sphere_quality =
        SettingGet_i(I->G, nullptr, nullptr, cSetting_cgo_sphere_quality);
  }

  std::unique_ptr<CGO> cgo_managed(CGONew(G, I->c + est));
  auto* const cgo = cgo_managed.get();
  RETURN_VAL_IF_FAIL(cgo, nullptr);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_COLOR:
      copy3f(pc, cgo->color);
      CGOColorv(cgo, pc);
      break;
    case CGO_PICK_COLOR:
      CGOPickColor(cgo, CGO_get_uint(pc), CGO_get_int(pc + 1));
      break;
    case CGO_SHADER_CYLINDER: {
      float v2[3];
      int cap = CGO_get_int(pc + 7);
      cCylCap fcap = cap1_from_cyl_shader_bits(cap);
      cCylCap bcap = cap2_from_cyl_shader_bits(cap);
      add3f(pc, pc + 3, v2);
      ok &= CGOSimpleCylinder(cgo, pc, v2, *(pc + 6), 0, 0, cgo->alpha,
          cgo->alpha, (cap & cCylShaderInterpColor), fcap, bcap, nullptr,
          stick_round_nub);
    } break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR: {
      auto cyl = it.cast<cgo::draw::shadercylinder2ndcolor>();
      float v1[3];
      int cap = cyl->cap;
      cCylCap fcap = cap1_from_cyl_shader_bits(cap);
      cCylCap bcap = cap2_from_cyl_shader_bits(cap);
      Pickable pickcolor2 = {cyl->pick_color_index, cyl->pick_color_bond};
      float color1[3] = {cgo->color[0], cgo->color[1], cgo->color[2]};
      add3f(pc, pc + 3, v1);
      float mid[3];
      mult3f(cyl->axis, .5f, mid);
      add3f(cyl->origin, mid, mid);
      float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : cgo->alpha;
      if (cap & cCylShaderInterpColor) {
        ok &= CGOSimpleCylinder(cgo, cyl->origin, v1, cyl->tube_size, color1,
            cyl->color2, cgo->alpha, alpha2, true, bcap, fcap, &pickcolor2,
            stick_round_nub);
      } else {
        ok &= CGOColorv(cgo, color1);
        ok &= CGOSimpleCylinder(cgo, cyl->origin, mid, cyl->tube_size, color1,
            nullptr, cgo->alpha, alpha2, false, fcap, cCylCap::None, nullptr,
            stick_round_nub);
        ok &= CGOColorv(cgo, cyl->color2);
        ok &= CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
        ok &= CGOSimpleCylinder(cgo, mid, v1, cyl->tube_size, cyl->color2,
            nullptr, cgo->alpha, alpha2, false, cCylCap::None, bcap, nullptr,
            stick_round_nub);
      }
    } break;
    case CGO_CYLINDER: {
      auto cyl = it.cast<cgo::draw::cylinder>();
      ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true,
          cCylCap::Flat, cCylCap::Flat, nullptr, stick_round_nub);
    } break;
    case CGO_CONE: {
      auto cone = it.cast<cgo::draw::cone>();
      ok &= CGOSimpleCone(cgo, *cone);
    } break;
    case CGO_SAUSAGE:
      ok &= CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10,
          cgo->alpha, cgo->alpha, true, cCylCap::Round, cCylCap::Round, nullptr,
          stick_round_nub);
      break;
    case CGO_CUSTOM_CYLINDER: {
      auto cyl = it.cast<cgo::draw::custom_cylinder>();
      ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true,
          cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
    } break;
    case CGO_CUSTOM_CYLINDER_ALPHA: {
      auto cyl = it.cast<cgo::draw::custom_cylinder_alpha>();
      ok &= CGOSimpleCylinder(cgo, *cyl, cyl->color1[3], cyl->color2[3], true,
          cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
    } break;
    case CGO_SPHERE:
      ok &= CGOSimpleSphere(cgo, pc, *(pc + 3), sphere_quality);
      break;
    case CGO_ELLIPSOID:
      ok &= CGOSimpleEllipsoid(cgo, pc, *(pc + 3), pc + 4, pc + 7, pc + 10);
      break;
    case CGO_QUADRIC:
      ok &= CGOSimpleQuadric(cgo, pc, *(pc + 3), pc + 4);
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    case CGO_DRAW_SPHERE_BUFFERS:
    case CGO_DRAW_CYLINDER_BUFFERS:
    case CGO_DRAW_LABELS:
    case CGO_DRAW_TEXTURES:
    case CGO_END:
    case CGO_VERTEX:
      WARN_UNEXPECTED_OPERATION(G, op);
      return nullptr;
    case CGO_BEGIN: {
      float firstColor[3], firstAlpha;
      char hasFirstColor = 0, hasFirstAlpha = 0;
      int nverts = 0, damode = CGO_VERTEX_ARRAY, err = 0;
      int mode = it.cast<cgo::draw::begin>()->mode;

      // remember for a second iteration
      auto it2 = it;

      for (++it;; ++it) {
        if (it.is_stop()) {
          WARN_UNEXPECTED_OPERATION(G, CGO_STOP);
          return nullptr;
        }

        const auto op = it.op_code();
        if (op == CGO_END) {
          break;
        }

        const auto pc = it.data();

        switch (op) {
        case CGO_DRAW_ARRAYS:
          WARN_UNEXPECTED_OPERATION(G, op);
          return nullptr;
        case CGO_NORMAL:
          damode |= CGO_NORMAL_ARRAY;
          break;
        case CGO_COLOR:
          if (!nverts) {
            hasFirstColor = 1;
            firstColor[0] = pc[0];
            firstColor[1] = pc[1];
            firstColor[2] = pc[2];
          } else {
            hasFirstColor = 0;
            damode |= CGO_COLOR_ARRAY;
          }
          break;
        case CGO_PICK_COLOR:
          damode |= CGO_PICK_COLOR_ARRAY;
          break;
        case CGO_ACCESSIBILITY:
          damode |= CGO_ACCESSIBILITY_ARRAY;
          break;
        case CGO_VERTEX:
          nverts++;
          break;
        case CGO_ALPHA:
          cgo->alpha = *pc;
          if (!nverts) {
            hasFirstAlpha = 1;
            firstAlpha = cgo->alpha;
          } else {
            hasFirstAlpha = 0;
            damode |= CGO_COLOR_ARRAY;
          }
        }
      }

      if (nverts > 0 && !err) {
        int pl = 0, plc = 0, pla = 0;
        float *vertexVals, *tmp_ptr;
        float *normalVals = 0, *colorVals = 0, *nxtVals = 0, *pickColorVals = 0,
              *accessibilityVals = 0;
        short notHaveValue = 0, nxtn = 3;
        if (hasFirstAlpha || hasFirstColor) {
          if (hasFirstAlpha) {
            CGOAlpha(cgo, firstAlpha);
          }
          if (hasFirstColor) {
            CGOColorv(cgo, firstColor);
          }
        }
        nxtVals = vertexVals =
            cgo->add<cgo::draw::arrays>(mode, damode, nverts);
        RETURN_VAL_IF_FAIL(vertexVals, nullptr);
        if (damode & CGO_NORMAL_ARRAY) {
          nxtVals = normalVals = vertexVals + (nxtn * nverts);
        }
        if (damode & CGO_COLOR_ARRAY) {
          nxtVals = colorVals = nxtVals + (nxtn * nverts);
          nxtn = VERTEX_COLOR_SIZE;
        }
        if (damode & CGO_PICK_COLOR_ARRAY) {
          nxtVals = nxtVals + (nxtn * nverts);
          pickColorVals = nxtVals + nverts;
          nxtn = VERTEX_PICKCOLOR_SIZE;
        }
        if (damode & CGO_ACCESSIBILITY_ARRAY) {
          nxtVals = nxtVals + (nxtn * nverts);
          accessibilityVals = nxtVals;
          nxtn = 1;
        }
        notHaveValue = damode;
        bool skiptoend = false;

        // second iteration
        for (++it2; !skiptoend; ++it2) {
          const auto op = it2.op_code();
          if (op == CGO_END) {
            break;
          }

          const auto pc = it2.data();

          switch (op) {
          case CGO_NORMAL:
            normalVals[pl] = pc[0];
            normalVals[pl + 1] = pc[1];
            normalVals[pl + 2] = pc[2];
            notHaveValue &= ~CGO_NORMAL_ARRAY;
            break;
          case CGO_COLOR:
            if (colorVals) {
              colorVals[plc] = pc[0];
              colorVals[plc + 1] = pc[1];
              colorVals[plc + 2] = pc[2];
              colorVals[plc + 3] = cgo->alpha;
              notHaveValue &= ~CGO_COLOR_ARRAY;
            }
            break;
          case CGO_PICK_COLOR:
            CGOPickColor(cgo, CGO_get_uint(pc), CGO_get_int(pc + 1));
            notHaveValue &= ~CGO_PICK_COLOR_ARRAY;
            break;
          case CGO_ACCESSIBILITY:
            cgo->current_accessibility = pc[0];
            break;
          case CGO_VERTEX:
            if (notHaveValue & CGO_NORMAL_ARRAY) {
              if (pl) {
                tmp_ptr = &normalVals[pl - 3];
                normalVals[pl] = tmp_ptr[0];
                normalVals[pl + 1] = tmp_ptr[1];
                normalVals[pl + 2] = tmp_ptr[2];
              } else {
                copy3f(cgo->normal, &normalVals[pl]);
              }
            }
            if (notHaveValue & CGO_COLOR_ARRAY) {
              if (plc) {
                tmp_ptr = &colorVals[plc - 4];
                colorVals[plc] = tmp_ptr[0];
                colorVals[plc + 1] = tmp_ptr[1];
                colorVals[plc + 2] = tmp_ptr[2];
                colorVals[plc + 3] = tmp_ptr[3];
              } else {
                copy3f(cgo->color, &colorVals[plc]);
                colorVals[plc + 3] = cgo->alpha;
              }
            }
            if (pickColorVals) {
              CGO_put_uint(
                  pickColorVals + pla * 2, cgo->current_pick_color_index);
              CGO_put_int(
                  pickColorVals + pla * 2 + 1, cgo->current_pick_color_bond);
            }
            if (accessibilityVals) {
              accessibilityVals[pla] = cgo->current_accessibility;
            }
            vertexVals[pl++] = pc[0];
            vertexVals[pl++] = pc[1];
            vertexVals[pl++] = pc[2];
            plc += 4;
            pla++;
            if (pla >= nverts) // anything past the last vertex is ignored
              skiptoend = true;
            notHaveValue = damode;
            break;
          case CGO_ALPHA:
            // in case we're before CGO_COLOR
            cgo->alpha = *pc;
            if (colorVals) {
              // in case we're after CGO_COLOR
              colorVals[plc + 3] = *pc;
            }
            break;
          }
        }
      }
    } break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
    default:
      cgo->add_to_cgo(op, pc);
    }

    if (G->Interrupt) {
      return nullptr;
    }

    RETURN_VAL_IF_FAIL(ok, nullptr);
  }

  CGOStop(cgo);
  return cgo_managed.release();
}

/**
 * converts a CGO that has primitives into pure geomtry, just like CGOSimplify
 *    but without converting the CGO_BEGIN/CGO_END blocks.
 *
 */
CGO* CGOSimplifyNoCompress(
    const CGO* I, int est, short sphere_quality, bool stick_round_nub)
{
  CGO* cgo;

  int ok = true;
  if (sphere_quality < 0) {
    sphere_quality =
        SettingGet_i(I->G, nullptr, nullptr, cSetting_cgo_sphere_quality);
  }

  cgo = CGONewSized(I->G, I->c + est);
  CHECKOK(ok, cgo);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_PICK_COLOR:
      CGOPickColor(cgo, CGO_get_uint(pc), CGO_get_int(pc + 1));
      break;
    case CGO_SHADER_CYLINDER: {
      float v2[3];
      int cap = CGO_get_int(pc + 7);
      cCylCap fcap = cap1_from_cyl_shader_bits(cap);
      cCylCap bcap = cap2_from_cyl_shader_bits(cap);
      add3f(pc, pc + 3, v2);
      ok &= CGOSimpleCylinder(cgo, pc, v2, *(pc + 6), 0, 0, cgo->alpha,
          cgo->alpha, (cap & cCylShaderInterpColor), fcap, bcap, nullptr,
          stick_round_nub);
    } break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR: {
      auto cyl = it.cast<cgo::draw::shadercylinder2ndcolor>();
      float v1[3];
      int cap = cyl->cap;
      cCylCap fcap = cap1_from_cyl_shader_bits(cap);
      cCylCap bcap = cap2_from_cyl_shader_bits(cap);
      Pickable pickcolor2 = {cyl->pick_color_index, cyl->pick_color_bond};
      float color1[3] = {cgo->color[0], cgo->color[1], cgo->color[2]};
      add3f(cyl->origin, cyl->axis, v1);
      float mid[3];
      mult3f(cyl->axis, .5f, mid);
      add3f(cyl->origin, mid, mid);
      float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : cgo->alpha;
      if (cap & cCylShaderInterpColor) {
        ok &= CGOSimpleCylinder(cgo, cyl->origin, v1, cyl->tube_size, color1,
            cyl->color2, cgo->alpha, alpha2, true, bcap, fcap, &pickcolor2,
            stick_round_nub);
      } else {
        ok &= CGOColorv(cgo, color1);
        ok &= CGOSimpleCylinder(cgo, cyl->origin, mid, cyl->tube_size, color1,
            nullptr, cgo->alpha, alpha2, false, fcap, cCylCap::None, nullptr,
            stick_round_nub);
        ok &= CGOColorv(cgo, cyl->color2);
        ok &= CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
        ok &= CGOSimpleCylinder(cgo, mid, v1, cyl->tube_size, cyl->color2,
            nullptr, cgo->alpha, alpha2, false, cCylCap::None, bcap, nullptr,
            stick_round_nub);
      }
    } break;
    case CGO_CYLINDER: {
      auto cyl = it.cast<cgo::draw::cylinder>();
      ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true,
          cCylCap::Flat, cCylCap::Flat, nullptr, stick_round_nub);
    } break;
    case CGO_CONE: {
      auto cone = it.cast<cgo::draw::cone>();
      ok &= CGOSimpleCone(cgo, *cone);
    } break;
    case CGO_SAUSAGE:
      ok &= CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10,
          cgo->alpha, cgo->alpha, true, cCylCap::Round, cCylCap::Round, nullptr,
          stick_round_nub);
      break;
    case CGO_CUSTOM_CYLINDER: {
      auto cyl = it.cast<cgo::draw::custom_cylinder>();
      ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true,
          cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
    } break;
    case CGO_CUSTOM_CYLINDER_ALPHA: {
      auto cyl = it.cast<cgo::draw::custom_cylinder_alpha>();
      ok &= CGOSimpleCylinder(cgo, *cyl, cyl->color1[3], cyl->color2[3], true,
          cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
    } break;
    case CGO_SPHERE:
      ok &= CGOSimpleSphere(cgo, pc, *(pc + 3), sphere_quality);
      break;
    case CGO_ELLIPSOID:
      ok &= CGOSimpleEllipsoid(cgo, pc, *(pc + 3), pc + 4, pc + 7, pc + 10);
      break;
    case CGO_QUADRIC:
      ok &= CGOSimpleQuadric(cgo, pc, *(pc + 3), pc + 4);
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "CGOSimplifyNoCompress-Error: CGO_DRAW_BUFFERS_INDEXED "
      "encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "CGOSimplifyNoCompress-Error: CGO_DRAW_BUFFERS_NOT_INDEXED "
      "encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_SPHERE_BUFFERS:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "CGOSimplifyNoCompress-Error: CGO_DRAW_SPHERE_BUFFERS "
      "encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_CYLINDER_BUFFERS:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "CGOSimplifyNoCompress-Error: CGO_DRAW_CYLINDER_BUFFERS "
      "encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_LABELS:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "CGOSimplifyNoCompress-Error: CGO_DRAW_LABELS encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_TEXTURES:
      PRINTFB(I->G, FB_CGO, FB_Errors)
      "CGOSimplifyNoCompress-Error: CGO_DRAW_TEXTURES encountered\n" ENDFB(
          I->G);
      break;
    case CGO_BEGIN:
      cgo->has_begin_end = true;
    default:
      cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  if (ok) {
    ok &= CGOStop(cgo);
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

CGO* CGOOptimizeTextures(const CGO* I, int est)
{
  CGO* cgo = nullptr;
  int ok = true;
  auto num_total_textures = CGOCountNumberOfOperationsOfType(I, CGO_DRAW_TEXTURE);
  //  printf("CGOOptimizeTextures: num_total_textures=%d\n",
  //  num_total_textures);
  if (num_total_textures) {
    int place3 = 0, place2 = 0;
    std::vector<float> worldPos(num_total_textures * 18);
    std::vector<float> screenValues(num_total_textures * 18);
    std::vector<float> textExtents(num_total_textures * 12);
    std::vector<float> pickColorValsVec(num_total_textures * 12); /* pick index and bond */
    auto* pickColorVals = pickColorValsVec.data();

    cgo = CGONewSized(I->G, 0);

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto op = it.op_code();
      const auto pc = it.data();

      switch (op) {
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeTextures() CGO_DRAW_BUFFERS_INDEXED or "
        "CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_TEXTURE: {
        float screenMin[3], screenMax[3], textExtent[4];
        copy3f(pc, &worldPos[place3]);
        copy3f(pc, &worldPos[place3 + 3]);
        copy3f(pc, &worldPos[place3 + 6]);
        copy3f(pc, &worldPos[place3 + 9]);
        copy3f(pc, &worldPos[place3 + 12]);
        copy3f(pc, &worldPos[place3 + 15]);
        copy3f(pc + 3, screenMin);
        copy3f(pc + 6, screenMax);
        copy4f(pc + 9, textExtent);
        copy3f(screenMin, &screenValues[place3]);
        copy3f(screenMin, &screenValues[place3 + 3]);
        copy3f(screenMin, &screenValues[place3 + 6]);
        copy3f(screenMin, &screenValues[place3 + 9]);
        copy3f(screenMin, &screenValues[place3 + 12]);
        copy3f(screenMax, &screenValues[place3 + 15]);
        screenValues[place3 + 4] = screenMax[1];
        screenValues[place3 + 6] = screenMax[0];
        screenValues[place3 + 10] = screenMax[1];
        screenValues[place3 + 12] = screenMax[0];
        screenValues[place3 + 17] = screenMin[2];
        place3 += 18;
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[0];
        textExtents[place2++] = textExtent[1];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[0];
        textExtents[place2++] = textExtent[3];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[2];
        textExtents[place2++] = textExtent[1];
        CGO_put_int(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_uint(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[0];
        textExtents[place2++] = textExtent[3];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[2];
        textExtents[place2++] = textExtent[1];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[2];
        textExtents[place2++] = textExtent[3];
      } break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
          buffer_layout::SEQUENTIAL);
      ok &= vbo->bufferData(
          {BufferDesc("attr_worldpos", VertexFormat::Float3,
               sizeof(float) * num_total_textures * 18, worldPos.data()),
              BufferDesc("attr_screenoffset", VertexFormat::Float3,
                  sizeof(float) * num_total_textures * 18, screenValues.data()),
              BufferDesc("attr_texcoords", VertexFormat::Float3,
                  sizeof(float) * num_total_textures * 18, textExtents.data())});
      size_t vboid = vbo->get_hash_id();

      if (ok) {
        float* pickArray =
            cgo->add<cgo::draw::textures>(num_total_textures, vboid);
        CHECKOK(ok, pickArray);
        if (!pickArray)
          I->G->ShaderMgr->freeGPUBuffer(vboid);
        if (ok)
          memcpy(pickArray + num_total_textures * 6, pickColorVals,
              num_total_textures * 12 * sizeof(float));
        if (ok)
          ok &= CGOStop(cgo);
      } else {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
      }
      if (!ok) {
        CGOFree(cgo);
      }
    }
  }
  return cgo;
}

CGO* CGOConvertToLabelShader(const CGO* I, CGO* addTo)
{
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;

  AttribDataOp world_pos_op = {{CGO_DRAW_LABEL, 1, FLOAT3_TO_FLOAT3,
      offsetof(cgo::draw::label, world_pos), 0}};
  AttribDataOp screen_offset_op = {{CGO_DRAW_LABEL, 2, FLOAT3_TO_FLOAT3,
      offsetof(cgo::draw::label, screen_world_offset), 0}};
  AttribDataOp screen_min_op = {{CGO_DRAW_LABEL, 3, FLOAT3_TO_FLOAT3,
      offsetof(cgo::draw::label, screen_min), 0}};
  AttribDataOp screen_max_op = {{CGO_DRAW_LABEL, 4, FLOAT3_TO_FLOAT3,
      offsetof(cgo::draw::label, screen_max), 0}};
  AttribDataOp text_extent_op = {{CGO_DRAW_LABEL, 5, FLOAT2_TO_FLOAT2,
      offsetof(cgo::draw::label, text_extent), 0}};
  AttribDataOp relative_mode_op = {{CGO_DRAW_LABEL, 6, FLOAT_TO_FLOAT,
      offsetof(cgo::draw::label, relative_mode), 0}};
  AttribDataOp target_pos_op = {{CGO_DRAW_LABEL, 7, FLOAT3_TO_FLOAT3,
      offsetof(cgo::draw::label, target_pos), 6}};

  AttribDataDesc attrDesc = {
      {"attr_worldpos", VertexFormat::Float3, world_pos_op},
      {"attr_targetpos", VertexFormat::Float3, target_pos_op},
      {"attr_screenoffset", VertexFormat::Float3, screen_min_op},
      {"attr_texcoords", VertexFormat::Float2, text_extent_op},
      {"attr_screenworldoffset", VertexFormat::Float3, screen_offset_op},
      {"attr_relative_mode", VertexFormat::Float, relative_mode_op}};

  auto ComputeScreenValues = [](void* varData, const float* pc,
                                 void* screenData, int idx) {
    auto sp = reinterpret_cast<const cgo::draw::label*>(pc);
    const auto& smin = sp->screen_min;
    const auto& smax = sp->screen_max;
    float* v = reinterpret_cast<float*>(varData);
    switch (idx) {
    case 0:
      v[0] = smin[0];
      v[1] = smin[1];
      v[2] = smin[2];
      break;
    case 1:
      v[0] = smin[0];
      v[1] = smax[1];
      v[2] = smin[2];
      break;
    case 2:
      v[0] = smax[0];
      v[1] = smin[1];
      v[2] = smin[2];
      break;
    case 3:
      v[0] = smin[0];
      v[1] = smax[1];
      v[2] = smin[2];
      break;
    case 4:
      v[0] = smax[0];
      v[1] = smin[1];
      v[2] = smin[2];
      break;
    case 5:
      v[0] = smax[0];
      v[1] = smax[1];
      v[2] = smin[2];
      break;
    };
  };

  auto ComputeTexCoords = [](void* varData, const float* pc, void* discard,
                              int idx) {
    auto sp = reinterpret_cast<const cgo::draw::label*>(pc);
    float* v = reinterpret_cast<float*>(varData);
    const auto& te = sp->text_extent;
    static struct {
      int x, y;
    } const idxs[6] = {{0, 1}, {0, 3}, {2, 1}, {0, 3}, {2, 1}, {2, 3}};
    v[0] = te[idxs[idx].x];
    v[1] = te[idxs[idx].y];
  };

  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {ComputeScreenValues, nullptr, "attr_screenoffset"});
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {ComputeTexCoords, nullptr, "attr_texcoords"});

  uchar pickdata[4] = {0, 0, 0, 0};
  addTo->add<cgo::draw::vertex_attribute_4ub>(
      G->ShaderMgr->GetAttributeUID("attr_pickcolor"), pickdata);

  AttribDataOp pickOp = {{CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataDesc pickDesc = {
      {"attr_pickcolor", VertexFormat::UByte4Norm, pickOp}};
  return CGOConvertToShader(
      I, attrDesc, pickDesc, GL_TRIANGLES, buffer_layout::INTERLEAVED, true);
}

CGO* CGOOptimizeLabels(const CGO* I, int est, bool addshaders)
{
  CGO* cgo = nullptr;
  int ok = true;
  auto num_total_labels = CGOCountNumberOfOperationsOfType(I, CGO_DRAW_LABEL);
  if (num_total_labels) {
    int place3 = 0, place2 = 0, place = 0;
    std::vector<float> worldPos(num_total_labels * 6 * 17);
    auto* screenValues = worldPos.data() + (num_total_labels * 18);
    auto* targetPos = screenValues + (num_total_labels * 18);
    auto* screenWorldValues = targetPos + (num_total_labels * 18);
    auto* textExtents = screenWorldValues + (num_total_labels * 18);
    auto* pickColorVals =
        textExtents + (num_total_labels * 12); /* pick index and bond */
    auto relativeMode = (float*) (pickColorVals + (num_total_labels * 12));
    cgo = CGONewSized(I->G, 0);

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto op = it.op_code();
      const auto pc = it.data();

      switch (op) {
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeLabels() CGO_DRAW_BUFFERS_INDEXED or "
        "CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_LABEL: {
        float screenWorldOffset[3], screenMin[3], screenMax[3], textExtent[4];
        copy3f(pc, &worldPos[place3]);
        copy3f(pc, &worldPos[place3 + 3]);
        copy3f(pc, &worldPos[place3 + 6]);
        copy3f(pc, &worldPos[place3 + 9]);
        copy3f(pc, &worldPos[place3 + 12]);
        copy3f(pc, &worldPos[place3 + 15]);
        copy3f(pc + 3, screenWorldOffset);
        copy3f(pc + 6, screenMin);
        copy3f(pc + 9, screenMax);
        copy4f(pc + 12, textExtent);
        copy3f(screenWorldOffset, &screenWorldValues[place3]);
        copy3f(&screenWorldValues[place3], &screenWorldValues[place3 + 3]);
        copy3f(&screenWorldValues[place3], &screenWorldValues[place3 + 6]);
        copy3f(&screenWorldValues[place3], &screenWorldValues[place3 + 9]);
        copy3f(&screenWorldValues[place3], &screenWorldValues[place3 + 12]);
        copy3f(&screenWorldValues[place3], &screenWorldValues[place3 + 15]);
        copy3f(screenMin, &screenValues[place3]);
        copy3f(screenMin, &screenValues[place3 + 3]);
        copy3f(screenMin, &screenValues[place3 + 6]);
        copy3f(screenMin, &screenValues[place3 + 9]);
        copy3f(screenMin, &screenValues[place3 + 12]);
        copy3f(screenMax, &screenValues[place3 + 15]);
        screenValues[place3 + 4] = screenMax[1];
        screenValues[place3 + 6] = screenMax[0];
        screenValues[place3 + 10] = screenMax[1];
        screenValues[place3 + 12] = screenMax[0];
        screenValues[place3 + 17] = screenMin[2];
        copy3f(pc + 17, &targetPos[place3]);
        copy3f(pc + 17, &targetPos[place3 + 3]);
        copy3f(pc + 17, &targetPos[place3 + 6]);
        copy3f(pc + 17, &targetPos[place3 + 9]);
        copy3f(pc + 17, &targetPos[place3 + 12]);
        copy3f(pc + 17, &targetPos[place3 + 15]);
        place3 += 18;
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[0];
        textExtents[place2++] = textExtent[1];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[0];
        textExtents[place2++] = textExtent[3];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[2];
        textExtents[place2++] = textExtent[1];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[0];
        textExtents[place2++] = textExtent[3];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[2];
        textExtents[place2++] = textExtent[1];
        CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
        CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
        textExtents[place2++] = textExtent[2];
        textExtents[place2++] = textExtent[3];
        {
          uchar rM = (uchar) * (pc + 16);
          relativeMode[place++] = rM;
          relativeMode[place++] = rM;
          relativeMode[place++] = rM;
          relativeMode[place++] = rM;
          relativeMode[place++] = rM;
          relativeMode[place++] = rM;
        }
      } break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      // Static Vertex Data
      VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
          buffer_layout::SEQUENTIAL);
      ok &=
          vbo->bufferData({BufferDesc("attr_worldpos", VertexFormat::Float3,
                               sizeof(float) * num_total_labels * 18, worldPos.data()),
              BufferDesc("attr_targetpos", VertexFormat::Float3,
                  sizeof(float) * num_total_labels * 18, targetPos),
              BufferDesc("attr_screenoffset", VertexFormat::Float3,
                  sizeof(float) * num_total_labels * 18, screenValues),
              BufferDesc("attr_texcoords", VertexFormat::Float2,
                  sizeof(float) * num_total_labels * 12, textExtents),
              BufferDesc("attr_screenworldoffset", VertexFormat::Float3,
                  sizeof(float) * num_total_labels * 18, screenWorldValues),
              BufferDesc("attr_relative_mode", VertexFormat::Float,
                  sizeof(float) * num_total_labels * 6, relativeMode)});
      size_t vboid = vbo->get_hash_id();

      VertexBuffer* pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
          buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
      ok &= pickvbo->bufferData(
          {BufferDesc("attr_pickcolor", VertexFormat::UByteNorm, 0),
              BufferDesc("attr_pickcolor", VertexFormat::UByteNorm,
                  sizeof(float) * num_total_labels * 6)},
          0, sizeof(float) * num_total_labels * 12, 0);
      size_t pickvboid = pickvbo->get_hash_id();

      auto freebuffers = [vboid, pickvboid, I]() {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
      };

      if (ok) {
        float* pickArray = nullptr;
        if (addshaders) {
          CGOEnable(cgo, GL_LABEL_SHADER);
        }
        pickArray =
            cgo->add<cgo::draw::labels>(num_total_labels, vboid, pickvboid);
        if (addshaders) {
          CGODisable(cgo, GL_LABEL_SHADER);
        }
        CHECKOK(ok, pickArray);
        if (!pickArray) {
          freebuffers();
        }
        if (ok)
          memcpy(pickArray + num_total_labels * 6, pickColorVals,
              num_total_labels * 12 * sizeof(float));
        if (ok)
          ok &= CGOStop(cgo);
      } else {
        freebuffers();
      }
      if (!ok) {
        CGOFree(cgo);
      }
    }
  }
  return cgo;
}

CGO* CGOOptimizeConnectors(const CGO* I, int est)
{
  CGO* cgo = nullptr;
  int num_total_connectors;
  int ok = true;
  int use_geometry_shaders =
      SettingGetGlobal_b(I->G, cSetting_use_geometry_shaders);
  int factor = (use_geometry_shaders ? 1 : 4);
  num_total_connectors =
      CGOCountNumberOfOperationsOfType(I, CGO_DRAW_CONNECTOR);

  if (num_total_connectors) {
    uchar* isCenterPt = nullptr;
    int place3 = 0, place2 = 0, place = 0;
    std::vector<float> targetPt3d(num_total_connectors * 20 *
                             factor, 0); /* too much, relativeMode only needs 1
                                         byte per vertex, instead of 1 float */
    auto* labelCenterPt3d = targetPt3d.data() + (num_total_connectors * 3 * factor);
    auto* indentFactor = labelCenterPt3d + (num_total_connectors * 3 * factor);
    auto* screenWorldOffset = indentFactor + (num_total_connectors * 2 * factor);
    auto* connectorColor = screenWorldOffset + (num_total_connectors * 3 * factor);
    auto* textSize = connectorColor + (num_total_connectors * factor);
    auto* relativeMode = (uchar*) (textSize + (num_total_connectors * 2 * factor));
    auto* drawBkgrd = (uchar*) (relativeMode + (num_total_connectors * factor));
    auto* bkgrdColor = (float*) (drawBkgrd + (num_total_connectors * factor));
    auto* relExtLength = (float*) (bkgrdColor + (num_total_connectors * factor));
    auto* connectorWidth = (float*) (relExtLength + (num_total_connectors * factor));
    if (!use_geometry_shaders)
      isCenterPt = (uchar*) (connectorWidth + (num_total_connectors * factor));
    else
      isCenterPt = nullptr;
    cgo = CGONewSized(I->G, 0);

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto op = it.op_code();
      const auto pc = it.data();

      switch (op) {
      case CGO_PICK_COLOR:
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        break;
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
        PRINTFB(I->G, FB_CGO, FB_Warnings)
        "WARNING: CGOOptimizeConnectors() CGO_DRAW_BUFFERS_INDEXED or "
        "CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
            op ENDFB(I->G);
        break;
      case CGO_DRAW_CONNECTOR: {
        uchar* uc;
        int f;
        if (!use_geometry_shaders) {
          isCenterPt[place] = 0;
          isCenterPt[place + 1] = 2;
          isCenterPt[place + 2] = 2;
          isCenterPt[place + 3] = 1;
        }
        copy3f(pc, &targetPt3d[place3]);
        copy3f(pc + 3, &labelCenterPt3d[place3]);
        copy2f(pc + 6, &indentFactor[place2]);
        copy3f(pc + 9, &screenWorldOffset[place3]);
        copy2f(pc + 12, &textSize[place2]);
        relativeMode[place] = (uchar) (int) pc[17];
        drawBkgrd[place] = (uchar) (int) pc[18];
        uc = (uchar*) &bkgrdColor[place];
        uc[0] = CLIP_COLOR_VALUE(*(pc + 19));
        uc[1] = CLIP_COLOR_VALUE(*(pc + 20));
        uc[2] = CLIP_COLOR_VALUE(*(pc + 21));
        uc[3] = CLIP_COLOR_VALUE(*(pc + 22));
        uc = (uchar*) &connectorColor[place];
        uc[0] = CLIP_COLOR_VALUE(*(pc + 14));
        uc[1] = CLIP_COLOR_VALUE(*(pc + 15));
        uc[2] = CLIP_COLOR_VALUE(*(pc + 16));
        uc[3] = 255;
        relExtLength[place] = *(pc + 8);
        connectorWidth[place] = *(pc + 23);
        place3 += 3;
        place2 += 2;
        place += 1;
        for (f = 1; f < factor; f++) {
          copy3f(pc, &targetPt3d[place3]);
          copy3f(pc + 3, &labelCenterPt3d[place3]);
          copy2f(pc + 6, &indentFactor[place2]);
          copy3f(pc + 9, &screenWorldOffset[place3]);
          copy2f(pc + 12, &textSize[place2]);
          relativeMode[place] = (uchar) (int) pc[17];
          drawBkgrd[place] = (uchar) (int) pc[18];
          relExtLength[place] = *(pc + 8);
          connectorWidth[place] = *(pc + 23);
          uc = (uchar*) &bkgrdColor[place];
          uc[0] = uc[-4];
          uc[1] = uc[-3];
          uc[2] = uc[-2];
          uc[3] = uc[-1];
          uc = (uchar*) &connectorColor[place];
          uc[0] = uc[-4];
          uc[1] = uc[-3];
          uc[2] = uc[-2];
          uc[3] = uc[-1];
          place3 += 3;
          place2 += 2;
          place += 1;
        }
      } break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      const size_t quant = factor * num_total_connectors;
      VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
      ok = vbo->bufferData({BufferDesc("a_target_pt3d", VertexFormat::Float3,
                                sizeof(float) * 3 * quant, targetPt3d.data()),
          BufferDesc("a_center_pt3d", VertexFormat::Float3,
              sizeof(float) * 3 * quant, labelCenterPt3d),
          BufferDesc("a_indentFactor", VertexFormat::Float2,
              sizeof(float) * 2 * quant, indentFactor),
          BufferDesc("a_screenWorldOffset", VertexFormat::Float3,
              sizeof(float) * 3 * quant, screenWorldOffset),
          BufferDesc("a_textSize", VertexFormat::Float2,
              sizeof(float) * 2 * quant, textSize),
          BufferDesc("a_Color", VertexFormat::UByte4Norm, sizeof(float) * quant,
              connectorColor),
          BufferDesc("a_relative_mode", VertexFormat::UByte,
              sizeof(std::uint8_t) * quant, relativeMode),
          BufferDesc("a_draw_flags", VertexFormat::UByte,
              sizeof(std::uint8_t) * quant, drawBkgrd),
          BufferDesc("a_bkgrd_color", VertexFormat::UByte4Norm,
              sizeof(float) * quant, bkgrdColor),
          BufferDesc("a_rel_ext_length", VertexFormat::Float,
              sizeof(float) * quant, relExtLength),
          BufferDesc("a_con_width", VertexFormat::Float, sizeof(float) * quant,
              connectorWidth),
          BufferDesc("a_isCenterPt", VertexFormat::UByte,
              sizeof(std::uint8_t) * quant, isCenterPt)});
      size_t vboid = vbo->get_hash_id();
      if (ok) {
        cgo->add<cgo::draw::connectors>(num_total_connectors, vboid);
        if (ok)
          ok &= CGOStop(cgo);
      }
      if (!ok) {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        CGOFree(cgo);
      }
    }
  }
  CheckGLErrorOK(I->G, "ERROR: CGOOptimizeConnectors() end returns err=%d\n");
  return cgo;
}

CGO* CGOExpandDrawTextures(const CGO* I, int est)
{
  CGO* cgo = CGONew(I->G);
  int ok = true;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      "WARNING: CGOOptimizeTextures() CGO_DRAW_BUFFERS_INDEXED or "
      "CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n",
          op ENDFB(I->G);
      break;
    case CGO_DRAW_TEXTURE: {
      float screenMin[3], screenMax[3], textExtent[4];
      float alpha = cgo->alpha;
      CGOAlpha(cgo, 0.f);
      CGOColor(cgo, 0.f, 0.f, 0.f);
      copy3f(pc + 3, screenMin);
      copy3f(pc + 6, screenMax);
      copy4f(pc + 9, textExtent);
      CGOBegin(cgo, GL_TRIANGLES);
      CGOTexCoord2f(cgo, textExtent[0], textExtent[1]);
      CGOVertexv(cgo, screenMin);
      CGOTexCoord2f(cgo, textExtent[0], textExtent[3]);
      CGOVertex(cgo, screenMin[0], screenMax[1], screenMin[2]);
      CGOTexCoord2f(cgo, textExtent[2], textExtent[1]);
      CGOVertex(cgo, screenMax[0], screenMin[1], screenMin[2]);
      CGOTexCoord2f(cgo, textExtent[0], textExtent[3]);
      CGOVertex(cgo, screenMin[0], screenMax[1], screenMin[2]);
      CGOTexCoord2f(cgo, textExtent[2], textExtent[1]);
      CGOVertex(cgo, screenMax[0], screenMin[1], screenMin[2]);
      CGOTexCoord2f(cgo, textExtent[2], textExtent[3]);
      CGOVertex(cgo, screenMax[0], screenMax[1], screenMin[2]);
      CGOEnd(cgo);
      CGOAlpha(cgo, alpha);
    } break;
    default:
      cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  CGOStop(cgo);
  return cgo;
}

/* ======== Raytrace Renderer ======== */

int CGOGetExtent(const CGO* I, float* mn, float* mx)
{
  int result = false;

#define check_extent(v, r)                                                     \
  {                                                                            \
    if (!result) {                                                             \
      mn[0] = ((*(v)) - r);                                                    \
      mx[0] = ((*(v)) + r);                                                    \
      mn[1] = ((*(v + 1)) - r);                                                \
      mx[1] = ((*(v + 1)) + r);                                                \
      mn[2] = ((*(v + 2)) - r);                                                \
      mx[2] = ((*(v + 2)) + r);                                                \
      result = true;                                                           \
    } else {                                                                   \
      if (mn[0] > ((*(v)) - r))                                                \
        mn[0] = ((*(v)) - r);                                                  \
      if (mx[0] < ((*(v)) + r))                                                \
        mx[0] = ((*(v)) + r);                                                  \
      if (mn[1] > ((*((v) + 1)) - r))                                          \
        mn[1] = ((*((v) + 1)) - r);                                            \
      if (mx[1] < ((*((v) + 1)) + r))                                          \
        mx[1] = ((*((v) + 1)) + r);                                            \
      if (mn[2] > ((*((v) + 2)) - r))                                          \
        mn[2] = ((*((v) + 2)) - r);                                            \
      if (mx[2] < ((*((v) + 2)) + r))                                          \
        mx[2] = ((*((v) + 2)) + r);                                            \
    }                                                                          \
  }

#define check_extent4(v, r)                                                    \
  {                                                                            \
    if (!result) {                                                             \
      mn[0] = ((*(v)) - r);                                                    \
      mx[0] = ((*(v)) + r);                                                    \
      mn[1] = ((*(v + 1)) - r);                                                \
      mx[1] = ((*(v + 1)) + r);                                                \
      mn[2] = ((*(v + 2)) - r);                                                \
      mx[2] = ((*(v + 2)) + r);                                                \
      mn[3] = ((*(v + 3)) - r);                                                \
      mx[3] = ((*(v + 3)) + r);                                                \
      result = true;                                                           \
    } else {                                                                   \
      if (mn[0] > ((*(v)) - r))                                                \
        mn[0] = ((*(v)) - r);                                                  \
      if (mx[0] < ((*(v)) + r))                                                \
        mx[0] = ((*(v)) + r);                                                  \
      if (mn[1] > ((*((v) + 1)) - r))                                          \
        mn[1] = ((*((v) + 1)) - r);                                            \
      if (mx[1] < ((*((v) + 1)) + r))                                          \
        mx[1] = ((*((v) + 1)) + r);                                            \
      if (mn[2] > ((*((v) + 2)) - r))                                          \
        mn[2] = ((*((v) + 2)) - r);                                            \
      if (mx[2] < ((*((v) + 2)) + r))                                          \
        mx[2] = ((*((v) + 2)) + r);                                            \
      if (mn[3] > ((*((v) + 3)) - r))                                          \
        mn[3] = ((*((v) + 3)) - r);                                            \
      if (mx[3] < ((*((v) + 3)) + r))                                          \
        mx[3] = ((*((v) + 3)) + r);                                            \
    }                                                                          \
  }

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
    case CGO_VERTEX:
      check_extent(pc, 0);
      break;
    case CGO_SPHERE:
    case CGO_ELLIPSOID:
      check_extent(pc, *(pc + 3));
      break;
    case CGO_CYLINDER:
    case CGO_CONE:
    case CGO_SAUSAGE:
    case CGO_CUSTOM_CYLINDER:
    case CGO_CUSTOM_CYLINDER_ALPHA:
      check_extent(pc, *(pc + 6));
      check_extent(pc + 3, *(pc + 6));
      break;
    case CGO_TRIANGLE:
      check_extent(pc, 0);
      check_extent(pc + 3, 0);
      check_extent(pc + 6, 0);
      break;
    case CGO_DRAW_ARRAYS: {
      const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      const float* pct = sp->floatdata;
      int pl;

      if (sp->arraybits & CGO_VERTEX_ARRAY) {
        for (pl = 0; pl < sp->nverts; pl++) {
          check_extent(pct, 0);
          pct += VERTEX_POS_SIZE;
        }
      }
    } break;
    case CGO_BOUNDING_BOX: {
      if (!result) {
        mn[0] = (*pc);
        mn[1] = *(pc + 1);
        mn[2] = *(pc + 2);
        mx[0] = *(pc + 3);
        mx[1] = *(pc + 4);
        mx[2] = *(pc + 5);
        result = true;
      } else {
        if (mn[0] > *pc)
          mn[0] = (*pc);
        if (mn[1] > *(pc + 1))
          mn[1] = *(pc + 1);
        if (mn[2] > *(pc + 2))
          mn[2] = *(pc + 2);
        if (mx[0] < *(pc + 3))
          mx[0] = *(pc + 3);
        if (mx[1] < *(pc + 4))
          mx[1] = *(pc + 4);
        if (mx[2] < *(pc + 5))
          mx[2] = *(pc + 5);
      }
    }
    }
  }
  return (result);
}

int CGOHasNormals(const CGO* I)
{
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
    case CGO_NORMAL:
    case CGO_SPHERE:
    case CGO_ELLIPSOID:
    case CGO_CYLINDER:
    case CGO_CONE:
    case CGO_SAUSAGE:
    case CGO_CUSTOM_CYLINDER:
    case CGO_CUSTOM_CYLINDER_ALPHA:
      return true;
    case CGO_DRAW_ARRAYS:
      if (it.cast<cgo::draw::arrays>()->arraybits & CGO_NORMAL_ARRAY) {
        return true;
      }
      break;
    }
  }
  return false;
}

static int CGOQuadricToEllipsoid(const float* v, float r, const float* q,
    float* r_el, float* n0, float* n1, float* n2)
{
  int ok = false;
  double inp_matrix[16];
  double e_val[4];
  double e_vec[16];
  double inverse[16];

  inp_matrix[0] = q[0];
  inp_matrix[1] = q[3];
  inp_matrix[2] = q[5];
  inp_matrix[3] = q[6];
  inp_matrix[4] = q[3];
  inp_matrix[5] = q[1];
  inp_matrix[6] = q[4];
  inp_matrix[7] = q[7];
  inp_matrix[8] = q[5];
  inp_matrix[9] = q[4];
  inp_matrix[10] = q[2];
  inp_matrix[11] = q[8];
  inp_matrix[12] = q[6];
  inp_matrix[13] = q[7];
  inp_matrix[14] = q[8];
  inp_matrix[15] = q[9];

  if (xx_matrix_invert(inverse, inp_matrix, 4)) {

    /* inverse now contains Uij coefficients */
    float pradius = sqrt1f(-1 / inverse[15]);
    int n_rot;

    if (xx_matrix_jacobi_solve(e_vec, e_val, &n_rot, inverse, 4)) {
      float mag[3];
      float scale[3];
      float mx;
      n0[0] = e_vec[0];
      n0[1] = e_vec[4];
      n0[2] = e_vec[8];
      n1[0] = e_vec[1];
      n1[1] = e_vec[5];
      n1[2] = e_vec[9];
      n2[0] = e_vec[2];
      n2[1] = e_vec[6];
      n2[2] = e_vec[10];

      normalize3f(n0);
      normalize3f(n1);
      normalize3f(n2);
      mag[0] = sqrt1f(e_val[0]);
      mag[1] = sqrt1f(e_val[1]);
      mag[2] = sqrt1f(e_val[2]);

      mx = mag[0];
      if (mx < mag[1])
        mx = mag[1];
      if (mx < mag[2])
        mx = mag[2];

      scale[0] = mag[0] / mx;
      scale[1] = mag[1] / mx;
      scale[2] = mag[2] / mx;

      scale3f(n0, scale[0], n0);
      scale3f(n1, scale[1], n1);
      scale3f(n2, scale[2], n2);

      *r_el = mx * pradius;
      ok = true;
    }
  }
  return ok;
}

static int CGORenderQuadricRay(CRay* ray, float* v, float r, float* q)
{
  float r_el, n0[3], n1[3], n2[3];
  int ok = true;
  if (CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    ok &= ray->ellipsoid3fv(v, r_el, n0, n1, n2);
  return ok;
}

/* ======== Raytrace Renderer ======== */

int CGORenderRay(CGO* I, CRay* ray, RenderInfo* info, const float* color,
    ObjectGadgetRamp* ramp, CSetting* set1, CSetting* set2)
{
#ifdef _PYMOL_NO_RAY
  return 0;
#else
  int vc = 0;
  float linewidth = 1.0F;
  float widthscale = 0.15F;
  float lineradius, dotradius, dotwidth;
  float white[] = {1.0, 1.0, 1.0};
  float zee[] = {0.0, 0.0, 1.0};
  int ok = true;
  const float *n0 = nullptr, *n1 = nullptr, *n2 = nullptr, *v0 = nullptr,
              *v1 = nullptr, *v2 = nullptr, *c0 = nullptr, *c1 = nullptr,
              *c2 = nullptr;
  float rampc0[3], rampc1[3], rampc2[3];
  int mode = -1;
  /* workaround; multi-state ray-trace bug */
  if (!I) {
    assert("TODO investigate" && false);
    return 0; /* not sure if it should return 0 or 1, 0 - fails, but is it a
                 memory issue? might not be since the arg is nullptr */
  }

  I->G->CGORenderer->alpha =
      1.0F - SettingGet_f(I->G, set1, set2, cSetting_cgo_transparency);

  widthscale = SettingGet_f(I->G, set1, set2, cSetting_cgo_ray_width_scale);

  /*  printf("debug %8.9f\n",SceneGetScreenVertexScale(I->G,zee)); */
  linewidth = SettingGet_f(I->G, set1, set2, cSetting_cgo_line_width);
  if (linewidth < 0.0F)
    linewidth = 1.0F;
  lineradius = SettingGet_f(I->G, set1, set2, cSetting_cgo_line_radius);
  dotwidth = SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_width);
  dotradius = SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_radius);
  if (lineradius < 0.0F)
    lineradius = linewidth * ray->PixelRadius / 2.0F;
  if (dotradius < 0.0F)
    dotradius = dotwidth * ray->PixelRadius / 2.0F;
  if (widthscale < 0.0F)
    widthscale = ray->PixelRadius / 2.0F;
  if (color)
    c0 = color;
  else
    c0 = white;
  ray->transparentf(1.0F - I->G->CGORenderer->alpha);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
    case CGO_BEGIN:
      mode = CGO_get_int(pc);
      vc = 0;
      n0 = zee;
      break;
    case CGO_END:
      switch (mode) {
      case GL_LINE_LOOP:
        if (vc > 1)
          ok &= ray->sausage3fv(v0, v2, lineradius, c0, c2);
        break;
      }
      mode = -1;
      break;
    case CGO_WIDTHSCALE:
      widthscale = *pc;
      lineradius = widthscale * linewidth;
      dotradius = widthscale * dotwidth;
      break;
    case CGO_DOTWIDTH:
      dotwidth = *pc;
      dotradius = widthscale * dotwidth;
      break;
    case CGO_LINEWIDTH:
      linewidth = *pc;
      lineradius = widthscale * linewidth;
      break;
    case CGO_NORMAL:
      n0 = pc;
      break;
    case CGO_SPECIAL_WITH_ARG: {
      float argval = *(pc + 1);
      switch (*(int*) pc) {
      case LINEWIDTH_FOR_LINES:
        linewidth = argval;
        lineradius = widthscale * linewidth;
      }
    } break;
    case CGO_SPECIAL: {
      switch ((*(int*) pc)) {
      case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON: {
        float radius = SettingGet_f(I->G, set1, set2, cSetting_ribbon_radius);
        if (radius == 0.0F) {
          float ribbon_width =
              SettingGet_f(I->G, set1, set2, cSetting_ribbon_width);
          float line_width = SceneGetDynamicLineWidth(info, ribbon_width);
          SceneGetDynamicLineWidth(info, line_width);
          radius = ray->PixelRadius * line_width / 2.0F;
        }
        lineradius = radius;
      } break;
      case LINEWIDTH_FOR_LINES: {
        float radius = SettingGet_f(I->G, set1, set2, cSetting_line_radius);
        if (radius <= 0.0F) {
          float line_width =
              SettingGet_f(I->G, set1, set2, cSetting_line_width);
          line_width = SceneGetDynamicLineWidth(info, line_width);
          radius = ray->PixelRadius * line_width / 2.0F;
        }
        lineradius = radius;
      } break;
      case CYLINDER_WIDTH_FOR_NONBONDED:
      case LINEWIDTH_WITH_SCALE: {
        float line_width = SettingGet_f(I->G, set1, set2, cSetting_line_width);
        line_width = SceneGetDynamicLineWidth(info, line_width);
        lineradius = widthscale * line_width / 2.f;
      } break;
      }
    } break;
    case CGO_COLOR:
      c0 = pc;
      ray->color3fv(c0);
      break;
    case CGO_ALPHA:
      I->G->CGORenderer->alpha = *pc;
      ray->transparentf(1.0F - *pc);
      break;
    case CGO_LINE: {
      auto line = reinterpret_cast<cgo::draw::line*>(pc);
      ok &= ray->sausage3fv(line->vertex1, line->vertex2, lineradius, c0, c0);
    } break;
    case CGO_SPLITLINE: {
      auto splitline = reinterpret_cast<cgo::draw::splitline*>(pc);
      float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
          CONVERT_COLOR_VALUE(splitline->color2[1]),
          CONVERT_COLOR_VALUE(splitline->color2[2])};
      if (splitline->flags & cgo::draw::splitline::interpolation) {
        ok &= ray->sausage3fv(
            splitline->vertex1, splitline->vertex2, lineradius, c0, color2);
      } else {
        float mid[3];
        add3f(splitline->vertex1, splitline->vertex2, mid);
        mult3f(mid, .5f, mid);
        ok &= ray->customCylinder3fv(splitline->vertex1, mid, lineradius, c0,
            c0, cCylCap::Round, cCylCap::None);
        ok &= ray->customCylinder3fv(mid, splitline->vertex2, lineradius,
            color2, color2, cCylCap::None, cCylCap::Round);
      }
    } break;
    case CGO_VERTEX_CROSS: {
      float pt1[3], pt2[3];
      float nonbonded_size =
          SettingGet_f(I->G, set1, set2, cSetting_nonbonded_size);
      copy3f(pc, pt1);
      copy3f(pc, pt2);
      pt1[0] -= nonbonded_size;
      pt2[0] += nonbonded_size;
      ok &= ray->sausage3fv(pt1, pt2, lineradius, c0, c0);

      copy3f(pc, pt1);
      copy3f(pc, pt2);
      pt1[1] -= nonbonded_size;
      pt2[1] += nonbonded_size;
      ok &= ray->sausage3fv(pt1, pt2, lineradius, c0, c0);

      copy3f(pc, pt1);
      copy3f(pc, pt2);
      pt1[2] -= nonbonded_size;
      pt2[2] += nonbonded_size;
      ok &= ray->sausage3fv(pt1, pt2, lineradius, c0, c0);
    } break;
    case CGO_VERTEX_BEGIN_LINE_STRIP:
    case CGO_VERTEX:
      v0 = pc;

      if (ramp) {
        if (!ObjectGadgetRampInterVertex(ramp, v0, rampc0, -1)) {
          copy3f(white, rampc0);
        }
        c0 = rampc0;
      }
      switch (mode) {
      case GL_POINTS:
        ok &= ray->sphere3fv(v0, dotradius);
        break;
      case GL_LINES:
        if (vc & 0x1)
          ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
        v1 = v0;
        if (!ramp) {
          c1 = c0;
        }
        break;
      case GL_LINE_STRIP:
        if (vc) {
          ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
        }
        v1 = v0;
        if (!ramp) {
          c1 = c0;
        }
        break;
      case GL_LINE_LOOP:
        if (vc)
          ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
        else {
          v2 = v0;
          c2 = c0;
        }
        v1 = v0;
        if (!ramp)
          c1 = c0;
        break;
      case GL_TRIANGLES:
        if (((vc + 1) % 3) == 0)
          ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
        v2 = v1;
        n2 = n1;
        v1 = v0;
        n1 = n0;
        if (!ramp) {
          c2 = c1;
          c1 = c0;
        }
        break;
      case GL_TRIANGLE_STRIP:
        if (vc > 1)
          ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
        v2 = v1;
        n2 = n1;
        v1 = v0;
        n1 = n0;
        if (!ramp) {
          c2 = c1;
          c1 = c0;
        }
        break;
      case GL_TRIANGLE_FAN:
        if (vc > 1)
          ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
        else if (!vc) {
          n2 = n0;
          v2 = v0;
          if (!ramp)
            c2 = c0;
        }
        v1 = v0;
        n1 = n0;
        if (!ramp)
          c1 = c0;
        break;
      }
      if (ramp) {
        switch (mode) {
        case GL_TRIANGLES:
        case GL_TRIANGLE_STRIP:
        case GL_TRIANGLE_FAN:
          copy3f(rampc1, rampc2);
          c2 = rampc2;
        case GL_LINES:
        case GL_LINE_STRIP:
        case GL_LINE_LOOP:
          copy3f(rampc0, rampc1);
          c1 = rampc1;
          break;
        }
      }
      vc++;
      break;
    case CGO_SPHERE:
      ray->color3fv(c0);
      ok &= ray->sphere3fv(pc, *(pc + 3));
      break;
    case CGO_ELLIPSOID:
      ray->color3fv(c0);
      ok &= ray->ellipsoid3fv(pc, *(pc + 3), pc + 4, pc + 7, pc + 10);
      break;
    case CGO_QUADRIC:
      ray->color3fv(c0);
      ok &= CGORenderQuadricRay(ray, pc, *(pc + 3), pc + 4);
      break;
    case CGO_CONE:
      ok &= ray->cone3fv(pc, pc + 3, *(pc + 6), *(pc + 7), pc + 8, pc + 11,
          static_cast<cCylCap>(int(pc[14])), static_cast<cCylCap>(int(pc[15])));
      break;
    case CGO_CUSTOM_CYLINDER: {
      auto cyl = reinterpret_cast<cgo::draw::custom_cylinder*>(pc);
      ok &= ray->customCylinder3fv(*cyl);
    } break;
    case CGO_CUSTOM_CYLINDER_ALPHA: {
      auto cyl = reinterpret_cast<cgo::draw::custom_cylinder_alpha*>(pc);
      ok &= ray->customCylinderAlpha3fv(*cyl);
    } break;
    case CGO_SHADER_CYLINDER: {
      float p2[3];
      int cap = CGO_get_int(pc + 7);
      const cCylCap cap1 = cap1_from_cyl_shader_bits(cap);
      const cCylCap cap2 = cap2_from_cyl_shader_bits(cap);
      add3f(pc, pc + 3, p2);
      ok &= ray->customCylinder3fv(
          pc, p2, *(pc + 6), ray->CurColor, ray->CurColor, cap1, cap2);
    } break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR: {
      auto cyl = reinterpret_cast<cgo::draw::shadercylinder2ndcolor*>(pc);
      float v1[3];
      int cap = cyl->cap;
      const cCylCap fcap = cap1_from_cyl_shader_bits(cap);
      const cCylCap bcap = cap2_from_cyl_shader_bits(cap);
      int colorinterp = cap & cCylShaderInterpColor;
      const float* color1 = c0;
      const float* color2 = cyl->color2;
      add3f(cyl->origin, cyl->axis, v1);
      float alpha1 = I->G->CGORenderer->alpha;
      float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : alpha1;
      if (colorinterp || equal3f(color1, color2)) {
        ok &= ray->customCylinder3fv(
            pc, v1, cyl->tube_size, color1, color2, fcap, bcap, alpha1, alpha2);
      } else {
        float mid[3];
        mult3f(cyl->axis, .5f, mid);
        add3f(cyl->origin, mid, mid);

        ray->color3fv(c0);
        ok &= ray->customCylinder3fv(cyl->origin, mid, cyl->tube_size, color1,
            color1, fcap, cCylCap::None, alpha1, alpha2);
        ray->color3fv(cyl->color2);
        ok &= ray->customCylinder3fv(mid, v1, cyl->tube_size, color2, color2,
            cCylCap::None, bcap, alpha1, alpha2);
      }
    } break;
    case CGO_CYLINDER: {
      auto* cyl = reinterpret_cast<cgo::draw::cylinder*>(pc);
      ok &= ray->cylinder3fv(*cyl);
    } break;
    case CGO_SAUSAGE:
      ok &= ray->sausage3fv(pc, pc + 3, *(pc + 6), pc + 7, pc + 10);
      break;
    case CGO_TRIANGLE:
      ok &= ray->triangle3fv(pc, pc + 3, pc + 6, pc + 9, pc + 12, pc + 15,
          pc + 18, pc + 21, pc + 24);
      break;
    case CGO_DRAW_ARRAYS: {
      cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      int const mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
      float* vertexVals = sp->floatdata;
      float *normalVals = 0, *colorVals = 0;
      int offset = 0;
      if (arrays & CGO_VERTEX_ARRAY) {
        vertexVals = sp->floatdata;
        offset += nverts * VERTEX_POS_SIZE;
      }
      if (arrays & CGO_NORMAL_ARRAY) {
        normalVals = sp->floatdata + offset;
        offset += nverts * VERTEX_NORMAL_SIZE;
      }
      if (arrays & CGO_COLOR_ARRAY) {
        colorVals = sp->floatdata + offset;
        offset += nverts * VERTEX_COLOR_SIZE;
      }
      vc = 0;
      for (int v = 0, pl = 0, plc = 0; ok && v < nverts;
           v++, pl += 3, plc += 4) {
        if (normalVals) {
          n0 = &normalVals[pl];
        }
        if (colorVals) {
          c0 = &colorVals[plc];
          ray->color3fv(c0);
          ray->transparentf(1.0f - c0[3]);
        }
        if (vertexVals) {
          v0 = &vertexVals[pl];
        }
        switch (mode) {
        case GL_POINTS:
          ok &= ray->sphere3fv(v0, dotradius);
          break;
        case GL_LINES:
          if (vc & 0x1)
            ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
          v1 = v0;
          c1 = c0;
          break;
        case GL_LINE_STRIP:
          if (vc)
            ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
          v1 = v0;
          c1 = c0;
          break;
        case GL_LINE_LOOP:
          if (vc)
            ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
          else {
            v2 = v0;
            c2 = c0;
          }
          v1 = v0;
          c1 = c0;
          break;
        case GL_TRIANGLES:
          if (((vc + 1) % 3) == 0)
            ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
          v2 = v1;
          c2 = c1;
          n2 = n1;
          v1 = v0;
          c1 = c0;
          n1 = n0;
          break;
        case GL_TRIANGLE_STRIP:
          if (vc > 1)
            ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
          v2 = v1;
          c2 = c1;
          n2 = n1;
          v1 = v0;
          c1 = c0;
          n1 = n0;
          break;
        case GL_TRIANGLE_FAN:
          if (vc > 1)
            ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
          else if (!vc) {
            n2 = n0;
            v2 = v0;
            c2 = c0;
          }
          v1 = v0;
          c1 = c0;
          n1 = n0;
          break;
        }
        vc++;
      }
    } break;
    default:
      break;
    }
  }

  if (ok)
    ray->transparentf(0.0F);
  return ok;
#endif
}

/**
 * Add axes-aligned cube to CGO with BEGIN/END GL_TRIANGLE_STRIP.
 *
 * @param center XYZ center position
 * @param radius Radius of a sphere with the same volume as the cube
 */
void CGOSimpleCube(CGO* I, const float* center, float radius)
{
  float const x = center[0];
  float const y = center[1];
  float const z = center[2];

  // match volume of sphere
  float const r = radius * 0.805996 /* (M_PI / 6)^(1/3) */;

  CGOBegin(I, GL_TRIANGLE_STRIP);
  CGONormal(I, 0., 0., 1.);
  CGOVertex(I, x + r, y + r, z + r);
  CGOVertex(I, x - r, y + r, z + r);
  CGOVertex(I, x + r, y - r, z + r);
  CGOVertex(I, x - r, y - r, z + r);
  CGOEnd(I);
  CGOBegin(I, GL_TRIANGLE_STRIP);
  CGONormal(I, 1., 0., 0.);
  CGOVertex(I, x + r, y - r, z - r);
  CGOVertex(I, x + r, y + r, z - r);
  CGOVertex(I, x + r, y - r, z + r);
  CGOVertex(I, x + r, y + r, z + r);
  CGOEnd(I);
  CGOBegin(I, GL_TRIANGLE_STRIP);
  CGONormal(I, 0., 1., 0.);
  CGOVertex(I, x + r, y + r, z - r);
  CGOVertex(I, x - r, y + r, z - r);
  CGOVertex(I, x + r, y + r, z + r);
  CGOVertex(I, x - r, y + r, z + r);
  CGOEnd(I);
  CGOBegin(I, GL_TRIANGLE_STRIP);
  CGONormal(I, 0., 0., -1.);
  CGOVertex(I, x - r, y - r, z - r);
  CGOVertex(I, x - r, y + r, z - r);
  CGOVertex(I, x + r, y - r, z - r);
  CGOVertex(I, x + r, y + r, z - r);
  CGOEnd(I);
  CGOBegin(I, GL_TRIANGLE_STRIP);
  CGONormal(I, -1., 0., 0.);
  CGOVertex(I, x - r, y + r, z + r);
  CGOVertex(I, x - r, y + r, z - r);
  CGOVertex(I, x - r, y - r, z + r);
  CGOVertex(I, x - r, y - r, z - r);
  CGOEnd(I);
  CGOBegin(I, GL_TRIANGLE_STRIP);
  CGONormal(I, 0., -1., 0.);
  CGOVertex(I, x - r, y - r, z + r);
  CGOVertex(I, x - r, y - r, z - r);
  CGOVertex(I, x + r, y - r, z + r);
  CGOVertex(I, x + r, y - r, z - r);
  CGOEnd(I);
}

/**
 * Add axes-aligned tetrahedron to CGO with BEGIN/END GL_TRIANGLES.
 *
 * @param center XYZ center position
 * @param radius Distance between center and (each) edge
 */
void CGOSimpleTetrahedron(CGO* I, const float* center, float radius)
{
  float vertices[][3] = {
      {1.f, 1.f, 1.f},
      {-1.f, -1.f, 1.f},
      {1.f, -1.f, -1.f},
      {-1.f, 1.f, -1.f},
  };

  for (float* const v : vertices) {
    scale3f(v, radius, v);
    add3f(v, center, v);
  }

  constexpr float n = 0.57735027;

  CGOBegin(I, GL_TRIANGLES);
  CGONormal(I, n, -n, n);
  CGOVertexv(I, vertices[0]);
  CGOVertexv(I, vertices[1]);
  CGOVertexv(I, vertices[2]);
  CGONormal(I, n, n, -n);
  CGOVertexv(I, vertices[0]);
  CGOVertexv(I, vertices[2]);
  CGOVertexv(I, vertices[3]);
  CGONormal(I, -n, n, n);
  CGOVertexv(I, vertices[0]);
  CGOVertexv(I, vertices[3]);
  CGOVertexv(I, vertices[1]);
  CGONormal(I, -n, -n, -n);
  CGOVertexv(I, vertices[1]);
  CGOVertexv(I, vertices[3]);
  CGOVertexv(I, vertices[2]);
  CGOEnd(I);
}

/* translation function which turns cylinders and spheres into triangles */

static int CGOSimpleSphere(
    CGO* I, const float* v, float vdw, short sphere_quality)
{
  SphereRec* sp;
  int *q, *s;
  int b, c;
  int ok = true;
  /* cgo_sphere_quality is between 0 and (NUMBER_OF_SPHERE_LEVELS-1) */

  sp = I->G->Sphere->Sphere[CLAMPVALUE<short>(
      sphere_quality, 0, (NUMBER_OF_SPHERE_LEVELS - 1))];

  q = sp->Sequence;
  s = sp->StripLen;

  for (b = 0; b < sp->NStrip; b++) {
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
    for (c = 0; ok && c < (*s); c++) {
      ok &= CGONormalv(I, sp->dot[*q]);
      if (ok)
        ok &= CGOVertex(I, v[0] + vdw * sp->dot[*q][0],
            v[1] + vdw * sp->dot[*q][1], v[2] + vdw * sp->dot[*q][2]);
      q++;
    }
    if (ok)
      ok &= CGOEnd(I);
    s++;
  }
  return ok;
}

static int CGOSimpleQuadric(CGO* I, const float* v, float r, const float* q)
{
  float r_el, n0[3], n1[3], n2[3];
  int ok = true;
  if (CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    ok &= CGOSimpleEllipsoid(I, v, r_el, n0, n1, n2);
  return ok;
}

static int CGOSimpleEllipsoid(CGO* I, const float* v, float vdw,
    const float* n0, const float* n1, const float* n2)
{
  SphereRec* sp;
  int *q, *s;
  int b, c;
  int ds;
  float nn0[3], nn1[3], nn2[3];
  float scale[3], scale_sq[3];
  int ok = true;

  normalize23f(n0, nn0);
  normalize23f(n1, nn1);
  normalize23f(n2, nn2);

  scale[0] = (float) length3f(n0);
  scale[1] = (float) length3f(n1);
  scale[2] = (float) length3f(n2);

  scale_sq[0] = scale[0] * scale[0];
  scale_sq[1] = scale[1] * scale[1];
  scale_sq[2] = scale[2] * scale[2];

  ds = SettingGet_i(I->G, nullptr, nullptr, cSetting_cgo_ellipsoid_quality);
  if (ds < 0)
    ds = SettingGet_i(I->G, nullptr, nullptr, cSetting_ellipsoid_quality);
  if (ds < 0)
    ds = 0;
  if (ds > 3)
    ds = 3;
  sp = I->G->Sphere->Sphere[ds];

  q = sp->Sequence;
  s = sp->StripLen;

  for (b = 0; b < sp->NStrip; b++) {
    ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
    for (c = 0; ok && c < (*s); c++) {
      float* sp_dot_q = sp->dot[*q];
      float s0 = vdw * sp_dot_q[0];
      float s1 = vdw * sp_dot_q[1];
      float s2 = vdw * sp_dot_q[2];
      float d0[3], d1[3], d2[3], vv[3], direction[3];
      float dd0, dd1, dd2, ss0, ss1, ss2;
      float comp0[3], comp1[3], comp2[3];
      float surfnormal[3];
      int i;

      scale3f(n0, s0, d0);
      scale3f(n1, s1, d1);
      scale3f(n2, s2, d2);

      for (i = 0; i < 3; i++) {
        vv[i] = d0[i] + d1[i] + d2[i];
      }
      normalize23f(vv, direction);
      add3f(v, vv, vv);

      dd0 = dot_product3f(direction, nn0);
      dd1 = dot_product3f(direction, nn1);
      dd2 = dot_product3f(direction, nn2);

      if (scale[0] > R_SMALL8) {
        ss0 = dd0 / scale_sq[0];
      } else {
        ss0 = 0.0F;
      }
      if (scale[1] > R_SMALL8) {
        ss1 = dd1 / scale_sq[1];
      } else {
        ss1 = 0.0F;
      }

      if (scale[2] > R_SMALL8) {
        ss2 = dd2 / scale_sq[2];
      } else {
        ss2 = 0.0F;
      }

      scale3f(nn0, ss0, comp0);
      scale3f(nn1, ss1, comp1);
      scale3f(nn2, ss2, comp2);

      for (i = 0; i < 3; i++) {
        surfnormal[i] = comp0[i] + comp1[i] + comp2[i];
      }
      normalize3f(surfnormal);

      ok &= CGONormalv(I, surfnormal);
      if (ok)
        ok &= CGOVertexv(I, vv);
      q++;
    }
    if (ok)
      ok &= CGOEnd(I);
    s++;
  }
  return ok;
}

/**
 * Triangulated round cap (half-globe)
 */
void CGORoundNub(CGO* I,
    const float* v1, // cap center
    const float* p0, // normal along axis
    const float* p1, // x coord in cap space
    const float* p2, // y coord in cap space
    int direction,   // 1 or -1
    int nEdge,       // "quality"
    float size)
{
  const int cmax = (nEdge + 3) / 2;
  const float PI_over_cmax = PI / ((cmax - 1) * 2);
  const float PI_over_nEdge = (PI * 2) / nEdge;
  float z2 = 1.f;

  // z coord in cap space
  float p3[3];
  scale3f(p0, direction, p3);

  CGOBegin(I, GL_TRIANGLE_STRIP);

  // from equator to pole (latitudinal)
  for (int c = 1; c < cmax; c += 1) {
    float z1 = z2;
    z2 = cos(c * PI_over_cmax);

    // around cylinder axis (longitudinal)
    for (int d = (nEdge + 1) * (-direction); d; d += direction) {
      float z3 = z1;

      // 2 vertices
      for (int e = -1; e < 1; ++e) {
        float x = cos(d * PI_over_nEdge) * sin((c + e) * PI_over_cmax);
        float y = sin(d * PI_over_nEdge) * sin((c + e) * PI_over_cmax);
        float normal[3], vertex[3];

        normal[0] = p1[0] * x + p2[0] * y + p3[0] * z3;
        normal[1] = p1[1] * x + p2[1] * y + p3[1] * z3;
        normal[2] = p1[2] * x + p2[2] * y + p3[2] * z3;

        vertex[0] = v1[0] + normal[0] * size;
        vertex[1] = v1[1] + normal[1] * size;
        vertex[2] = v1[2] + normal[2] * size;

        normalize3f(normal);
        CGONormalv(I, normal);
        CGOVertexv(I, vertex);

        z3 = z2;
      }
    }
  }

  CGOEnd(I);
}

static int CGOSimpleCylinder(CGO* I, const float* v1, const float* v2,
    const float tube_size, const float* c1, const float* c2, const float alpha1,
    const float alpha2, const bool interp, const cCylCap cap1,
    const cCylCap cap2, const Pickable* pickcolor2, const bool stick_round_nub)
{
#define MAX_EDGE 50

  float d[3], t[3], p0[3], p1[3], p2[3], vv1[3], vv2[3], v_buf[9], *v;
  float x[MAX_EDGE + 1], y[MAX_EDGE + 1];
  float overlap;
  float nub;
  bool colorFlag, interpColorFlag;
  int nEdge;
  int c;
  int ok = true;
  float midcolor[3];
  float midalpha{alpha1};
  Pickable pickcolor[2];
  pickcolor[0].index = I->current_pick_color_index;
  pickcolor[0].bond = I->current_pick_color_bond;
  if (pickcolor2) {
    pickcolor[1].index = pickcolor2->index;
    pickcolor[1].bond = pickcolor2->bond;
  } else {
    pickcolor[1].index = pickcolor[0].index;
    pickcolor[1].bond = pickcolor[0].bond;
  }
  v = v_buf;
  nEdge = SettingGetGlobal_i(I->G, cSetting_stick_quality);
  overlap = tube_size * SettingGetGlobal_f(I->G, cSetting_stick_overlap);
  nub = tube_size * SettingGetGlobal_f(I->G, cSetting_stick_nub);

  if (nEdge > MAX_EDGE)
    nEdge = MAX_EDGE;
  subdivide(nEdge, x, y);

  colorFlag = (c1 != c2) && c2;
  colorFlag |= alpha1 != alpha2;

  interpColorFlag = c2 && interp && pickcolor2;
  if (interpColorFlag) {
    average3f(c1, c2, midcolor);
    midalpha = (alpha1 + alpha2) / 2.0f;
  }
  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  normalize3f(p0);

  if (cap1 == cCylCapRound && !stick_round_nub) {
    vv1[0] = v1[0] - p0[0] * overlap;
    vv1[1] = v1[1] - p0[1] * overlap;
    vv1[2] = v1[2] - p0[2] * overlap;
  } else {
    vv1[0] = v1[0];
    vv1[1] = v1[1];
    vv1[2] = v1[2];
  }
  if (cap2 == cCylCapRound && !stick_round_nub) {
    vv2[0] = v2[0] + p0[0] * overlap;
    vv2[1] = v2[1] + p0[1] * overlap;
    vv2[2] = v2[2] + p0[2] * overlap;
  } else {
    vv2[0] = v2[0];
    vv2[1] = v2[1];
    vv2[2] = v2[2];
  }

  d[0] = (vv2[0] - vv1[0]);
  d[1] = (vv2[1] - vv1[1]);
  d[2] = (vv2[2] - vv1[2]);
  if (pickcolor2) {
    mult3f(d, .5f, d);
  }
  get_divergent3f(d, t);
  cross_product3f(d, t, p1);
  normalize3f(p1);

  cross_product3f(d, p1, p2);

  normalize3f(p2);

  /* now we have a coordinate system */

  if (ok)
    ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
  for (c = nEdge; ok && c >= 0; c--) {
    v[0] = p1[0] * x[c] + p2[0] * y[c];
    v[1] = p1[1] * x[c] + p2[1] * y[c];
    v[2] = p1[2] * x[c] + p2[2] * y[c];

    v[3] = vv1[0] + v[0] * tube_size;
    v[4] = vv1[1] + v[1] * tube_size;
    v[5] = vv1[2] + v[2] * tube_size;

    v[6] = v[3] + d[0];
    v[7] = v[4] + d[1];
    v[8] = v[5] + d[2];

    ok &= CGONormalv(I, v);
    if (ok && (colorFlag || interpColorFlag)) {
      ok &= CGOColorv(I, c1);
      ok &= CGOAlpha(I, alpha1);
    }
    if (ok)
      ok &= CGOVertexv(I, v + 3);
    if (ok && interpColorFlag) {
      ok &= CGOColorv(I, midcolor);
      ok &= CGOAlpha(I, midalpha);
    } else if (ok && colorFlag && !pickcolor2) {
      ok &= CGOColorv(I, c2);
      ok &= CGOAlpha(I, alpha2);
    }
    if (ok)
      ok &= CGOVertexv(I, v + 6);
  }
  if (ok)
    ok &= CGOEnd(I);
  if (pickcolor2) {
    ok &= CGOColorv(I, c2);
    ok &= CGOAlpha(I, alpha2);
    CGOPickColor(I, pickcolor2->index, pickcolor2->bond);
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
    for (c = nEdge; ok && c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv1[0] + v[0] * tube_size + d[0];
      v[4] = vv1[1] + v[1] * tube_size + d[1];
      v[5] = vv1[2] + v[2] * tube_size + d[2];

      v[6] = v[3] + d[0];
      v[7] = v[4] + d[1];
      v[8] = v[5] + d[2];

      ok &= CGONormalv(I, v);
      if (ok && interpColorFlag) {
        ok &= CGOColorv(I, midcolor);
        ok &= CGOAlpha(I, midalpha);
      }
      if (ok)
        ok &= CGOVertexv(I, v + 3);
      if (ok && interpColorFlag) {
        ok &= CGOColorv(I, c2);
        ok &= CGOAlpha(I, alpha2);
      }
      if (ok)
        ok &= CGOVertexv(I, v + 6);
    }
    if (ok)
      ok &= CGOEnd(I);
  }

  if (ok && cap1 != cCylCap::None) {
    if (ok && colorFlag && c1) {
      ok &= CGOColorv(I, c1);
      ok &= CGOAlpha(I, alpha1);
    }
    if (pickcolor2)
      CGOPickColor(I, pickcolor[0].index, pickcolor[0].bond);

    if (stick_round_nub && cap1 == cCylCapRound) {
      CGORoundNub(I, v1, p0, p1, p2, -1, nEdge, tube_size);
    } else {
      v[0] = -p0[0];
      v[1] = -p0[1];
      v[2] = -p0[2];

      if (cap1 == cCylCapRound) {
        v[3] = vv1[0] - p0[0] * nub;
        v[4] = vv1[1] - p0[1] * nub;
        v[5] = vv1[2] - p0[2] * nub;
      } else {
        v[3] = vv1[0];
        v[4] = vv1[1];
        v[5] = vv1[2];
      }

      if (ok)
        ok &= CGOBegin(I, GL_TRIANGLE_FAN);
      if (ok)
        ok &= CGONormalv(I, v);
      if (ok)
        ok &= CGOVertexv(I, v + 3);

      for (c = nEdge; ok && c >= 0; c--) {
        v[0] = p1[0] * x[c] + p2[0] * y[c];
        v[1] = p1[1] * x[c] + p2[1] * y[c];
        v[2] = p1[2] * x[c] + p2[2] * y[c];

        v[3] = vv1[0] + v[0] * tube_size;
        v[4] = vv1[1] + v[1] * tube_size;
        v[5] = vv1[2] + v[2] * tube_size;

        if (cap1 == cCylCapRound)
          ok &= CGONormalv(I, v);
        if (ok)
          ok &= CGOVertexv(I, v + 3);
      }
      if (ok)
        ok &= CGOEnd(I);
    }
  }

  if (ok && cap2 != cCylCap::None) {
    if (ok && colorFlag && c2) {
      ok &= CGOColorv(I, c2);
      ok &= CGOAlpha(I, alpha2);
    }
    if (pickcolor2)
      CGOPickColor(I, pickcolor2->index, pickcolor2->bond);

    if (stick_round_nub && cap2 == cCylCapRound) {
      CGORoundNub(I, v2, p0, p1, p2, 1, nEdge, tube_size);
    } else {
      v[0] = p0[0];
      v[1] = p0[1];
      v[2] = p0[2];

      if (cap2 == cCylCapRound) {
        v[3] = vv2[0] + p0[0] * nub;
        v[4] = vv2[1] + p0[1] * nub;
        v[5] = vv2[2] + p0[2] * nub;
      } else {
        v[3] = vv2[0];
        v[4] = vv2[1];
        v[5] = vv2[2];
      }

      if (ok)
        ok &= CGOBegin(I, GL_TRIANGLE_FAN);
      if (ok)
        ok &= CGONormalv(I, v);
      if (ok)
        ok &= CGOVertexv(I, v + 3);

      for (c = 0; ok && c <= nEdge; c++) {
        v[0] = p1[0] * x[c] + p2[0] * y[c];
        v[1] = p1[1] * x[c] + p2[1] * y[c];
        v[2] = p1[2] * x[c] + p2[2] * y[c];

        v[3] = vv2[0] + v[0] * tube_size;
        v[4] = vv2[1] + v[1] * tube_size;
        v[5] = vv2[2] + v[2] * tube_size;

        if (cap2 == cCylCapRound)
          ok &= CGONormalv(I, v);
        if (ok)
          ok &= CGOVertexv(I, v + 3);
      }
      if (ok)
        ok &= CGOEnd(I);
    }
  }
  return ok;
}

template <typename CylinderT>
static int CGOSimpleCylinder(CGO* I, const CylinderT& cyl, const float a1,
    const float a2, const bool interp, const cCylCap cap1, const cCylCap cap2,
    const Pickable* pickcolor2, const bool stick_round_nub)
{
  return CGOSimpleCylinder(I, cyl.vertex1, cyl.vertex2, cyl.radius, cyl.color1,
      cyl.color2, a1, a2, interp, cap1, cap2, pickcolor2, stick_round_nub);
}

static int CGOSimpleCone(CGO* I, const float* v1, const float* v2, float r1,
    float r2, const float* c1, const float* c2, cCylCap cap1, cCylCap cap2)
{
#define MAX_EDGE 50

  float d[3], t[3], p0[3], p1[3], p2[3], vv1[3], vv2[3], v_buf[9], *v;
  float x[MAX_EDGE + 1], y[MAX_EDGE + 1], edge_normal[3 * (MAX_EDGE + 1)];
  int colorFlag;
  int nEdge;
  int c;
  int ok = true;

  v = v_buf;
  nEdge = SettingGetGlobal_i(I->G, cSetting_cone_quality);

  if (nEdge > MAX_EDGE)
    nEdge = MAX_EDGE;
  subdivide(nEdge, x, y);

  colorFlag = (c1 != c2) && c2;

  ok &= CGOColorv(I, c1);

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  normalize3f(p0);

  {
    vv1[0] = v1[0];
    vv1[1] = v1[1];
    vv1[2] = v1[2];
  }

  {
    vv2[0] = v2[0];
    vv2[1] = v2[1];
    vv2[2] = v2[2];
  }

  d[0] = (vv2[0] - vv1[0]);
  d[1] = (vv2[1] - vv1[1]);
  d[2] = (vv2[2] - vv1[2]);

  get_divergent3f(d, t);

  cross_product3f(d, t, p1);

  normalize3f(p1);

  cross_product3f(d, p1, p2);

  normalize3f(p2);

  /* now we have a coordinate system */

  {
    float len = diff3f(v1, v2);
    float vt[3], nt[3];
    float slope = 0.0F;

    if (len) {
      slope = (r1 - r2) / len;
    }
    for (c = nEdge; c >= 0; c--) {
      vt[0] = p1[0] * x[c] + p2[0] * y[c];
      vt[1] = p1[1] * x[c] + p2[1] * y[c];
      vt[2] = p1[2] * x[c] + p2[2] * y[c];

      scale3f(p0, slope, nt);
      add3f(nt, vt, vt);
      normalize3f(vt);
      copy3f(vt, edge_normal + 3 * c);
    }
  }

  /* now we have normals */
  if (ok)
    ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
  for (c = nEdge; ok && c >= 0; c--) {
    v[0] = p1[0] * x[c] + p2[0] * y[c];
    v[1] = p1[1] * x[c] + p2[1] * y[c];
    v[2] = p1[2] * x[c] + p2[2] * y[c];

    v[3] = vv1[0] + v[0] * r1;
    v[4] = vv1[1] + v[1] * r1;
    v[5] = vv1[2] + v[2] * r1;

    v[6] = vv1[0] + v[0] * r2 + d[0];
    v[7] = vv1[1] + v[1] * r2 + d[1];
    v[8] = vv1[2] + v[2] * r2 + d[2];

    ok &= CGONormalv(I, edge_normal + 3 * c);
    if (ok && colorFlag)
      CGOColorv(I, c1);
    if (ok)
      CGOVertexv(I, v + 3);
    if (ok && colorFlag)
      CGOColorv(I, c2);
    if (ok)
      CGOVertexv(I, v + 6);
  }
  if (ok)
    ok &= CGOEnd(I);

  if (ok && cap1 != cCylCap::None) {
    v[0] = -p0[0];
    v[1] = -p0[1];
    v[2] = -p0[2];

    {
      v[3] = vv1[0];
      v[4] = vv1[1];
      v[5] = vv1[2];
    }

    if (colorFlag)
      ok &= CGOColorv(I, c1);
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_FAN);
    if (ok)
      ok &= CGONormalv(I, v);
    if (ok)
      ok &= CGOVertexv(I, v + 3);

    for (c = nEdge; ok && c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv1[0] + v[0] * r1;
      v[4] = vv1[1] + v[1] * r1;
      v[5] = vv1[2] + v[2] * r1;

      if (cap1 == cCylCapRound)
        ok &= CGONormalv(I, v);
      if (ok)
        ok &= CGOVertexv(I, v + 3);
    }
    if (ok)
      ok &= CGOEnd(I);
  }

  if (ok && cap2 != cCylCap::None) {

    v[0] = p0[0];
    v[1] = p0[1];
    v[2] = p0[2];

    {
      v[3] = vv2[0];
      v[4] = vv2[1];
      v[5] = vv2[2];
    }

    if (colorFlag)
      ok &= CGOColorv(I, c2);
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_FAN);
    if (ok)
      ok &= CGONormalv(I, v);
    if (ok)
      ok &= CGOVertexv(I, v + 3);

    for (c = 0; ok && c <= nEdge; c++) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv2[0] + v[0] * r2;
      v[4] = vv2[1] + v[1] * r2;
      v[5] = vv2[2] + v[2] * r2;

      if (cap2 == cCylCapRound)
        ok &= CGONormalv(I, v);
      if (ok)
        ok &= CGOVertexv(I, v + 3);
    }
    if (ok)
      ok &= CGOEnd(I);
  }
  return ok;
}

int CGOSimpleCone(CGO* I, const cgo::draw::cone& cone)
{
  return CGOSimpleCone(I, cone.vertex1, cone.vertex2, cone.radius1,
      cone.radius2, cone.color1, cone.color2, static_cast<cCylCap>(cone.cap1),
      static_cast<cCylCap>(cone.cap2));
}

/* CGOGetNextDrawBufferedIndex: This is used by RepSurface to */
/* get the data from the CGO_DRAW_BUFFERS_INDEXED operation so */
/* that it can update the indices for semi-transparent surfaces. */
const cgo::draw::buffers_not_indexed* CGOGetNextDrawBufferedNotIndex(
    const CGO* cgo)
{
  for (auto it = cgo->begin(); !it.is_stop(); ++it) {
    if (it.op_code() == CGO_DRAW_BUFFERS_NOT_INDEXED) {
      return it.cast<cgo::draw::buffers_not_indexed>();
    }
  }
  return nullptr;
}

bool CGO::append(const CGO& source_, bool stopAtEnd)
{
  const CGO* const source = &source_;
  int ok = 1;

  for (auto it = source->begin(); !it.is_stop(); ++it) {
    add_to_cgo(it.op_code(), it.data());
  }

  if (stopAtEnd)
    ok &= CGOStop(this);
  has_draw_buffers |= source->has_draw_buffers;
  has_draw_cylinder_buffers |= source->has_draw_cylinder_buffers;
  return ok;
}

/**
 * Appends `src` to the end of this CGO. Takes ownership of data
 * (incl. VBOs) and leaves `src` as a valid but empty CGO.
 */
void CGO::move_append(CGO&& src_)
{
  CGO* const src = &src_;
  if (!src->c)
    return;

  // copy buffer
  VLACheck(op, float, c + src->c);
  memcpy(op + c, src->op, src->c * sizeof(float));

  // update sizes
  c += src->c;
  src->c = 0;

  // null terminators (CGO_STOP)
  *(op + c) = 0;
  *(src->op) = 0;

  // move heap data
  for (auto& ref : src->_data_heap) {
    _data_heap.emplace_back(std::move(ref));
  }
  src->_data_heap.clear();

  // copy boolean flags
  has_draw_buffers |= src->has_draw_buffers;
  has_draw_cylinder_buffers |= src->has_draw_cylinder_buffers;
  has_draw_sphere_buffers |= src->has_draw_sphere_buffers;
  has_begin_end |= src->has_begin_end;
  use_shader |= src->use_shader;
  render_alpha |= src->render_alpha;
  src->has_draw_buffers = false;
}

/**
 * Appends `src` to the end of this CGO and then free's `src`
 * and sets the pointer to nullptr.
 */
void CGO::free_append(CGO*& src)
{
  free_append(std::move(src));
  assert(src == nullptr);
}
void CGO::free_append(CGO*&& src)
{
  if (!src)
    return;
  move_append(std::move(*src));
  DeleteP(src);
}

int CGOAppend(CGO* dest, const CGO* source, bool stopAtEnd)
{
  return dest->append(*source, stopAtEnd);
}

int CGOCountNumberOfOperationsOfType(const CGO* I, int optype)
{
  std::set<int> ops = {optype};
  return CGOCountNumberOfOperationsOfTypeN(I, ops);
}

int CGOCountNumberOfOperationsOfTypeN(const CGO* I, const std::set<int>& optype)
{
  int numops = 0;
  for (auto cgoit = I->begin(); !cgoit.is_stop(); ++cgoit) {
    if (optype.count(cgoit.op_code()))
      numops++;
  }
  return (numops);
}

int CGOCountNumberOfOperationsOfTypeN(
    const CGO* I, const std::map<int, int>& optype)
{
  int numops = 0;
  for (auto cgoit = I->begin(); !cgoit.is_stop(); ++cgoit) {
    auto it = optype.find(cgoit.op_code());
    if (it != optype.end())
      numops += it->second;
  }
  return (numops);
}

bool CGOHasOperationsOfType(const CGO* I, int optype)
{
  std::set<int> ops = {optype};
  return CGOHasOperationsOfTypeN(I, ops);
}

bool CGOHasOperations(const CGO* I)
{
  return !I->begin().is_stop();
}

bool CGOHasOperationsOfTypeN(const CGO* I, const std::set<int>& optype)
{
  if (!I->op)
    return false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (optype.count(it.op_code()))
      return 1;
  }
  return (0);
}

static bool CGOFilterOutOperationsOfTypeN(
    const CGO* I, CGO* cgo, const std::set<int>& optype)
{
  if (!I->op)
    return false;

  bool ret = false;
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto op = it.op_code();
    if (optype.find(op) == optype.end()) {
      auto pc = it.data();
      cgo->add_to_cgo(op, pc);
    } else {
      ret = true; // returns if filtered anything
    }
  }
  return ret;
}

bool CGOFilterOutCylinderOperationsInto(const CGO* I, CGO* cgo)
{
  static std::set<int> optypes = {CGO_SHADER_CYLINDER,
      CGO_SHADER_CYLINDER_WITH_2ND_COLOR, CGO_SAUSAGE, CGO_CYLINDER,
      CGO_CUSTOM_CYLINDER, CGO_CUSTOM_CYLINDER_ALPHA};
  return CGOFilterOutOperationsOfTypeN(I, cgo, optypes);
}

bool CGOFilterOutBezierOperationsInto(const CGO* I, CGO* cgo)
{
  static std::set<int> optypes = {CGO_BEZIER};
  return CGOFilterOutOperationsOfTypeN(I, cgo, optypes);
}

bool CGOHasCylinderOperations(const CGO* I)
{
  static std::set<int> optypes = {CGO_SHADER_CYLINDER,
      CGO_SHADER_CYLINDER_WITH_2ND_COLOR, CGO_SAUSAGE, CGO_CYLINDER,
      CGO_CUSTOM_CYLINDER, CGO_CUSTOM_CYLINDER_ALPHA};
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOHasSphereOperations(const CGO* I)
{
  static std::set<int> optypes = {CGO_SPHERE};
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOHasBezierOperations(const CGO* I)
{
  static std::set<int> optypes = {CGO_BEZIER};
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOCheckWhetherToFree(PyMOLGlobals* G, CGO* I)
{
  if (I->use_shader) {
    if (I->cgo_shader_ub_color !=
            SettingGetGlobal_i(G, cSetting_cgo_shader_ub_color) ||
        I->cgo_shader_ub_normal !=
            SettingGetGlobal_i(G, cSetting_cgo_shader_ub_normal)) {
      return true;
    }
  }
  return false;
}

CGO* CGOConvertLinesToShaderCylinders(const CGO* I, int est)
{

  int tot_nverts = 0, tot_ncyls = 0;

  CGO* cgo = CGONewSized(I->G, I->c + est);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_DRAW_ARRAYS: {
      const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      float* vals =
          cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
      int nvals = sp->narrays * sp->nverts;
      memcpy(vals, sp->floatdata, nvals);
    } break;
    case CGO_END:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGOConvertLinesToShaderCylinders: CGO_END encountered without "
      "CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);
      break;
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGOConvertLinesToShaderCylinders: CGO_VERTEX encountered without "
      "CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);
      break;
    case CGO_BEGIN: {
      const float *last_vertex = nullptr, *last_color = nullptr,
                  *current_color = nullptr, *color = nullptr;
      unsigned int last_pick_color_idx = 0;
      int last_pick_color_bnd = cPickableNoPick;
      int nverts = 0, err = 0;
      int mode = CGO_get_int(pc);

      for (++it; !err && it != CGO_END; ++it) {
        auto pc = it.data();
        int op = it.op_code();

        switch (op) {
        case CGO_VERTEX:
          if (last_vertex) {
            switch (mode) {
            case GL_LINES:
            case GL_LINE_STRIP: {
              float axis[3];
              bool pick_color_diff = false;
              axis[0] = pc[0] - last_vertex[0];
              axis[1] = pc[1] - last_vertex[1];
              axis[2] = pc[2] - last_vertex[2];
              pick_color_diff =
                  (cgo->current_pick_color_index != last_pick_color_idx ||
                      cgo->current_pick_color_bond != last_pick_color_bnd);
              if (last_color && current_color &&
                  (!equal3f(last_color, current_color) || pick_color_diff)) {
                CGOColorv(cgo, last_color);
                if (pick_color_diff) {
                  Pickable pickcolor2 = {cgo->current_pick_color_index,
                      cgo->current_pick_color_bond};
                  CGOPickColor(cgo, last_pick_color_idx, last_pick_color_bnd);
                  cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, last_vertex,
                      axis, 1.f, cCylShaderBothCapsRound, current_color,
                      &pickcolor2);
                  CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
                } else {
                  cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, last_vertex,
                      axis, 1.f, cCylShaderBothCapsRound, current_color);
                }
                CGOColorv(cgo, current_color);
              } else {
                cgo->add<cgo::draw::shadercylinder>(
                    last_vertex, axis, 1.f, cCylShaderBothCapsRound);
              }
              last_vertex = pc;
              last_pick_color_idx = cgo->current_pick_color_index;
              last_pick_color_bnd = cgo->current_pick_color_bond;
              tot_ncyls++;
            }
              if (mode == GL_LINES) {
                last_vertex = nullptr;
                last_color = nullptr;
              }
            }
          } else {
            last_vertex = pc;
            current_color = color;
            last_pick_color_idx = cgo->current_pick_color_index;
            last_pick_color_bnd = cgo->current_pick_color_bond;
          }
          nverts++;
          break;
        case CGO_LINE: {
          float axis[3];
          auto line = reinterpret_cast<const cgo::draw::line*>(pc);
          subtract3f(line->vertex2, line->vertex1, axis);
          cgo->add<cgo::draw::shadercylinder>(
              line->vertex1, axis, 1.f, cCylShaderBothCapsRound);
          tot_ncyls++;
        } break;
        case CGO_SPLITLINE: {
          float axis[3];
          auto splitline = reinterpret_cast<const cgo::draw::splitline*>(pc);
          Pickable pickcolor2 = {splitline->index, splitline->bond};
          float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
              CONVERT_COLOR_VALUE(splitline->color2[1]),
              CONVERT_COLOR_VALUE(splitline->color2[2])};
          unsigned char flags = splitline->flags;
          subtract3f(splitline->vertex2, splitline->vertex1, axis);
          if ((flags & cgo::draw::splitline::equal_colors) &&
              (flags & cgo::draw::splitline::no_split_for_pick)) {
            cgo->add<cgo::draw::shadercylinder>(
                splitline->vertex1, axis, 1., cCylShaderBothCapsRound);
          } else {
            int cap = cCylShaderBothCapsRound;
            if (flags & splitline->flags &
                cgo::draw::splitline::interpolation) {
              cap |= cCylShaderInterpColor;
            }
            cgo->add<cgo::draw::shadercylinder2ndcolor>(
                cgo, splitline->vertex1, axis, 1., cap, color2, &pickcolor2);
            last_pick_color_idx = splitline->index;
            last_pick_color_bnd = splitline->bond;
          }
          tot_ncyls++;
        } break;
        case CGO_COLOR:
          if (op == CGO_COLOR) {
            last_color = current_color;
            current_color = pc;
            color = pc;
          }
        case CGO_PICK_COLOR:
          if (op == CGO_PICK_COLOR) {
            cgo->current_pick_color_index = CGO_get_uint(pc);
            cgo->current_pick_color_bond = CGO_get_int(pc + 1);
          }
        default:
          cgo->add_to_cgo(op, pc);
        }
      }

      tot_nverts += nverts;
    } break;
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  if (tot_ncyls) {
    return (cgo);
  } else {
    CGOFree(cgo);
    return nullptr;
  }
}

#ifdef WITH_UNUSED_FUNCTIONS
/* CGOSplitUpLinesForPicking: This operation goes through */
/* a CGO and returns a new CGO that has the same lines but */
/* a line that has two different pick colors will get split */
/* at its midpoint into two separate lines so that it can */
/* be used for picking */
CGO* CGOSplitUpLinesForPicking(const CGO* I)
{
  auto G = I->G;

  std::unique_ptr<CGO> cgo_managed(new CGO(G));
  CGO* cgo = cgo_managed.get();
  int tot_nverts = 0;

  CGOBegin(cgo, GL_LINES);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_PICK_COLOR: {
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
    } break;
    case CGO_END:
    case CGO_VERTEX:
      WARN_UNEXPECTED_OPERATION(G, op);
      return nullptr;
    case CGO_BEGIN: {
      const float *last_vertex = nullptr, *last_color = nullptr,
                  *current_color = nullptr, *color = nullptr;
      unsigned int last_pick_color_idx = 0;
      int last_pick_color_bnd = cPickableNoPick;
      int nverts = 0;
      const int mode = it.cast<cgo::draw::begin>()->mode;

      for (++it;; ++it) {
        if (it.is_stop()) {
          WARN_UNEXPECTED_OPERATION(G, CGO_STOP);
          return nullptr;
        }

        const auto op = it.op_code();
        if (op == CGO_END) {
          break;
        }

        const auto pc = it.data();

        switch (op) {
        case CGO_VERTEX:
          if (last_vertex) {
            switch (mode) {
            case GL_LINES:
            case GL_LINE_STRIP: {
              bool pick_color_diff = false;
              pick_color_diff =
                  (cgo->current_pick_color_index != last_pick_color_idx ||
                      cgo->current_pick_color_bond != last_pick_color_bnd);
              if (pick_color_diff ||
                  (last_color && current_color &&
                      (!equal3f(last_color, current_color)))) {
                if (pick_color_diff) {
                  float haxis[3];
                  float mid[3];
                  uint curp_idx = cgo->current_pick_color_index;
                  int curp_bnd = cgo->current_pick_color_bond;
                  haxis[0] = .5f * (pc[0] - last_vertex[0]);
                  haxis[1] = .5f * (pc[1] - last_vertex[1]);
                  haxis[2] = .5f * (pc[2] - last_vertex[2]);
                  add3f(last_vertex, haxis, mid);
                  CGOPickColor(cgo, last_pick_color_idx, last_pick_color_bnd);
                  CGOVertexv(cgo, last_vertex);
                  CGOVertexv(cgo, mid);
                  CGOPickColor(cgo, curp_idx, curp_bnd);
                  CGOVertexv(cgo, mid);
                  CGOVertexv(cgo, pc);
                } else {
                  CGOVertexv(cgo, last_vertex);
                  CGOVertexv(cgo, pc);
                }
              } else {
                CGOVertexv(cgo, last_vertex);
                CGOVertexv(cgo, pc);
              }
              last_vertex = pc;
              last_pick_color_idx = cgo->current_pick_color_index;
              last_pick_color_bnd = cgo->current_pick_color_bond;
            }
              if (mode == GL_LINES) {
                last_vertex = nullptr;
                last_color = nullptr;
              }
            }
          } else {
            last_vertex = pc;
            current_color = color;
            last_pick_color_idx = cgo->current_pick_color_index;
            last_pick_color_bnd = cgo->current_pick_color_bond;
          }
          nverts++;
          break;
        case CGO_COLOR: {
          last_color = current_color;
          current_color = pc;
          color = pc;
        } break;
        case CGO_PICK_COLOR: {
          cgo->current_pick_color_index = CGO_get_uint(pc);
          cgo->current_pick_color_bond = CGO_get_int(pc + 1);
        } break;
        }
      }
      tot_nverts += nverts;
    } break;
    }
  }

  if (!tot_nverts) {
    return nullptr;
  }

  CGOEnd(cgo);
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }

  return cgo_managed.release();
}
#endif

static void trilinesBufferAddVertex(float*& buffer,
    const float* v1,    // vertex
    const float* v2,    // vertex other end of line
    const float* color, // RGB color
    float alpha,        // alpha
    signed char uv)     // uv
{
  // vertex
  (*buffer++) = v1[0];
  (*buffer++) = v1[1];
  (*buffer++) = v1[2];

  // othervertex
  (*buffer++) = v2[0];
  (*buffer++) = v2[1];
  (*buffer++) = v2[2];

  (*buffer++) = (float) uv;

  // RGBA
  unsigned char* byte_view = (unsigned char*) (buffer++);
  (*byte_view++) = CLIP_COLOR_VALUE(color[0]);
  (*byte_view++) = CLIP_COLOR_VALUE(color[1]);
  (*byte_view++) = CLIP_COLOR_VALUE(color[2]);
  (*byte_view++) = CLIP_COLOR_VALUE(alpha);
}

static void trilinesBufferAddVertices(float*& buffer,
    const float* v1,    // vertex
    const float* v2,    // vertex other end of line
    const float* color, // RGB color
    float alpha)        // alpha
{
  // Vertex 1
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 1); //-1, 1);
  // Vertex 3
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 3); // 1, 1);
  // Vertex 2
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 0); //-1, -1);
  // Vertex 4
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 3); // 1, 1);
  // Vertex 3
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 2); // 1, -1);
  // Vertex 1
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 1); //-1, 1);
}

static void CGOTrilines_GetCurrentColor(const float*& current_color,
    const float* colorv, const float* last_color, const float* cc)
{
  if (!current_color) {
    if (colorv) {
      current_color = colorv;
    } else if (last_color) {
      current_color = last_color;
    } else {
      current_color = cc;
    }
  }
}

/**
 * Changes
 *
 *     CGO_ENABLE <frommode>
 *
 * to
 *
 *     CGO_ENABLE <tomode>
 *
 * in place.
 */
void CGOChangeShadersTo(CGO* I, int frommode, int tomode)
{
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (it.op_code() == CGO_ENABLE) {
      auto eo = it.cast<cgo::draw::enable>();
      if (eo->mode == frommode) {
        eo->mode = tomode;
      }
    }
  }
}

CGO* CGOOptimizeScreenTexturesAndPolygons(CGO* I, int est)
{
  auto G = I->G;
  int ok = true;

  int num_total_indices = CGOCountNumVerticesForScreen(I);
  if (num_total_indices <= 0)
    return nullptr;

  std::unique_ptr<CGO> cgo_managed(new CGO(G));
  auto* const cgo = cgo_managed.get();

  {
    float *colorVals = 0, *texcoordVals;
    int tot, nxtn;
    uchar* colorValsUC = 0;
    CGOAlpha(cgo, 1.f);
    cgo->alpha = 1.f;
    cgo->color[0] = 1.f;
    cgo->color[1] = 1.f;
    cgo->color[2] = 1.f;

    {
      int mul =
          6; // 3 - screenoffset/vertex, 2 - texture coordinates, 1 - color
      /*
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
        mul++;
      } else {
        mul += 4;
        }*/
      tot = num_total_indices * mul;
    }

    auto vertexVals_managed = std::vector<float>(tot);
    float* vertexVals = vertexVals_managed.data();
    texcoordVals = vertexVals + 3 * num_total_indices;
    nxtn = 2;
    colorVals = texcoordVals + nxtn * num_total_indices;
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
    ok = CGOProcessScreenCGOtoArrays(
        G, I, vertexVals, texcoordVals, colorVals, colorValsUC);
    RETURN_VAL_IF_FAIL(ok && !G->Interrupt, nullptr);
    if (ok) {
      VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
      ok = vbo->bufferData(
          {BufferDesc("attr_screenoffset", VertexFormat::Float3,
               sizeof(float) * num_total_indices * 3, vertexVals),
              BufferDesc("attr_texcoords", VertexFormat::Float2,
                  sizeof(float) * num_total_indices * 2, texcoordVals),
              BufferDesc("attr_backgroundcolor", VertexFormat::UByte4Norm,
                  sizeof(std::uint8_t) * num_total_indices * 4, colorValsUC)});
      size_t vboid = vbo->get_hash_id();
      if (ok) {
        CGOEnable(cgo, GL_SCREEN_SHADER);
        cgo->add<cgo::draw::screen_textures>(num_total_indices, vboid);
        if (ok)
          ok &= CGODisable(cgo, GL_SCREEN_SHADER);
        RETURN_VAL_IF_FAIL(ok, nullptr);
      } else {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
      }
    }
    cgo->use_shader = true;
  }

  return cgo_managed.release();
}

CGO* CGOColorByRamp(PyMOLGlobals* G, const CGO* I, ObjectGadgetRamp* ramp,
    int state, CSetting* set1)
{
  if (!I) {
    return nullptr;
  }

  auto cgo = CGONewSized(G, 0);
  if (!cgo) {
    return nullptr;
  }

  int ok = true;
  float white[3] = {1.f, 1.f, 1.f};
  float probe_radius = SettingGet_f(G, set1, nullptr, cSetting_solvent_radius);
  float v_above[3], n0[3] = {0.f, 0.f, 0.f};
  int ramp_above =
      SettingGet_i(G, set1, nullptr, cSetting_surface_ramp_above_mode) == 1;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    bool skipCopy = false;
    switch (op) {
    case CGO_DRAW_ARRAYS: {
      const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      float* vals =
          cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
      int nvals = sp->narrays * sp->nverts;
      ok &= vals ? true : false;
      if (ok)
        memcpy(vals, sp->floatdata, nvals);
      skipCopy = true;
    } break;
    case CGO_NORMAL:
      copy3f(pc, n0);
      break;
    case CGO_VERTEX: {
      float color[3];
      copy3f(white, color);
      if (ramp_above) {
        copy3f(n0, v_above);
        scale3f(v_above, probe_radius, v_above);
        add3f(pc, v_above, v_above);
      } else {
        copy3f(pc, v_above);
      }
      if (ObjectGadgetRampInterVertex(ramp, v_above, color, state)) {
        CGOColorv(cgo, color);
      } else {
        CGOColorv(cgo, white);
      }
    } break;
    }
    if (!skipCopy) {
      cgo->add_to_cgo(op, pc);
    }
  }
  if (ok) {
    ok &= CGOStop(cgo);
    if (ok) {
      cgo->use_shader = I->use_shader;
      if (cgo->use_shader) {
        cgo->cgo_shader_ub_color =
            SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
        cgo->cgo_shader_ub_normal =
            SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
      }
    }
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

/**
 * FIXME: This function always returns true for `checkOpaque=ture`
 */
int CGOHasTransparency(const CGO* I, bool checkTransp, bool checkOpaque)
{
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (it.op_code() == CGO_ALPHA) {
      const auto pc = it.data();
      if (checkTransp && *pc < 1.f)
        return 1;
      if (checkOpaque && *pc == 1.f)
        return 1;
    }
  }

  return checkOpaque;
}

CGO* CGOConvertTrianglesToAlpha(const CGO* I)
{
  int tot_nverts = 0;

  CGO* cgo = CGONewSized(I->G, I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_DRAW_ARRAYS: {
      const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      const int mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
      const float* nxtVals = sp->floatdata;

      assert(arrays & CGO_VERTEX_ARRAY);
      const float* vertexValsDA = nxtVals;
      nxtVals += VERTEX_POS_SIZE * nverts;

      const float* normalValsDA = nullptr;
      if (arrays & CGO_NORMAL_ARRAY) {
        normalValsDA = nxtVals;
        nxtVals += VERTEX_NORMAL_SIZE * nverts;
      }

      const float* colorValsDA = nullptr;
      if (arrays & CGO_COLOR_ARRAY) {
        colorValsDA = nxtVals;
        nxtVals += VERTEX_COLOR_SIZE * nverts;
      }

      const float* const vertexVals0 = vertexValsDA;
      const float* const normalVals0 = normalValsDA;
      const float* const colorVals0 = colorValsDA;

      switch (mode) {
      case GL_TRIANGLES: {
        for (int cnt = 0; cnt < nverts; cnt += 3) {
          if (colorVals0) {
            CGOAlphaTriangle(cgo, vertexValsDA, vertexValsDA + 3,
                vertexValsDA + 6, normalValsDA, normalValsDA + 3,
                normalValsDA + 6, colorValsDA, colorValsDA + 4, colorValsDA + 8,
                *(colorValsDA + 3), *(colorValsDA + 7), *(colorValsDA + 11), 0);
          } else {
            CGOAlphaTriangle(cgo, vertexValsDA, vertexValsDA + 3,
                vertexValsDA + 6, normalValsDA, normalValsDA + 3,
                normalValsDA + 6, cgo->color, cgo->color, cgo->color,
                cgo->alpha, cgo->alpha, cgo->alpha, 0);
          }
          vertexValsDA += 9;
          normalValsDA += 9;
          if (colorVals0)
            colorValsDA += 12;
        }
      }
        tot_nverts += nverts;
        break;
      case GL_TRIANGLE_STRIP: {
        short flip = 0;
        vertexValsDA += 6;
        normalValsDA += 6;
        if (colorVals0)
          colorValsDA += 8;

        for (int cnt = 2; cnt < nverts; cnt++) {
          if (colorVals0) {
            CGOAlphaTriangle(cgo, vertexValsDA - 6, vertexValsDA - 3,
                vertexValsDA, normalValsDA - 6, normalValsDA - 3, normalValsDA,
                colorValsDA - 8, colorValsDA - 4, colorValsDA,
                *(colorValsDA - 5), *(colorValsDA - 1), *(colorValsDA + 3),
                flip);
          } else {
            CGOAlphaTriangle(cgo, vertexValsDA - 6, vertexValsDA - 3,
                vertexValsDA, normalValsDA - 6, normalValsDA - 3, normalValsDA,
                cgo->color, cgo->color, cgo->color, cgo->alpha, cgo->alpha,
                cgo->alpha, flip);
          }
          vertexValsDA += 3;
          normalValsDA += 3;
          if (colorVals0)
            colorValsDA += 4;
          flip = !flip;
        }
      }
        tot_nverts += nverts;
        break;
      case GL_TRIANGLE_FAN: {
        vertexValsDA += 6;
        normalValsDA += 6;
        if (colorVals0)
          colorValsDA += 8;
        for (int cnt = 2; cnt < nverts; cnt++) {
          if (colorVals0) {
            CGOAlphaTriangle(cgo, vertexVals0, vertexValsDA - 3, vertexValsDA,
                normalVals0, normalValsDA - 3, normalValsDA, colorVals0,
                colorValsDA - 4, colorValsDA, *(colorVals0 + 3),
                *(colorValsDA - 1), *(colorValsDA + 3), 0);
          } else {
            CGOAlphaTriangle(cgo, vertexVals0, vertexValsDA - 3, vertexValsDA,
                normalVals0, normalValsDA - 3, normalValsDA, cgo->color,
                cgo->color, cgo->color, cgo->alpha, cgo->alpha, cgo->alpha, 0);
          }
          vertexValsDA += 3;
          normalValsDA += 3;
          if (colorVals0)
            colorValsDA += 4;
        }
      }
        tot_nverts += nverts;
        break;
      }
    } break;
    case CGO_END:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGOConvertTrianglesToAlpha: CGO_END encountered without CGO_BEGIN but "
      "skipped for OpenGLES\n" ENDFB(I->G);
      break;
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
      " CGOConvertTrianglesToAlpha: CGO_VERTEX encountered without CGO_BEGIN "
      "but skipped for OpenGLES\n" ENDFB(I->G);
      break;
    case CGO_BEGIN: {
      float vertices[3][3], colors[4][3], normals[4][3], alpha[4];
      short verticespl = 2, colorspl = 2, normalspl = 2, alphapl = 2;
      short hasShifted = 0;
      int nverts = 0, err = 0;
      int mode = CGO_get_int(pc);
      short mode_is_triangles = 0, flip = 0, mode_is_fan = 0;
      copy3f(cgo->color, colors[3]);
      copy3f(cgo->normal, normals[3]);
      alpha[3] = cgo->alpha;
      switch (mode) {
      case GL_TRIANGLE_FAN:
        mode_is_fan = 1;
      case GL_TRIANGLES:
      case GL_TRIANGLE_STRIP:
        mode_is_triangles = 1;
      }
      if (!mode_is_triangles) {
        CGOBegin(cgo, mode);
      }

      for (++it; !err && it != CGO_END; ++it) {
        auto pc = it.data();
        int op = it.op_code();
        short add_to_cgo = 1;
        switch (op) {
        case CGO_VERTEX:
          if (mode_is_triangles) {
            if (!(hasShifted & 1)) { // colors
              if (colorspl >= 0) {
                copy3f(colors[colorspl + 1], colors[colorspl]);
                colorspl--;
              } else {
                if (!mode_is_fan)
                  copy3f(colors[1], colors[2]);
                copy3f(colors[0], colors[1]);
              }
            }
            if (!(hasShifted & 2)) { // normals
              if (normalspl >= 0) {
                copy3f(normals[normalspl + 1], normals[normalspl]);
                normalspl--;
              } else {
                if (!mode_is_fan)
                  copy3f(normals[1], normals[2]);
                copy3f(normals[0], normals[1]);
              }
            }
            if (!(hasShifted & 4)) { // alphas
              if (alphapl >= 0) {
                alpha[alphapl] = alpha[alphapl + 1];
                alphapl--;
              } else {
                if (!mode_is_fan)
                  alpha[2] = alpha[1];
                alpha[1] = alpha[0];
              }
            }
            if (verticespl >= 0) {
              copy3f(pc, vertices[verticespl]);
              verticespl--;
            } else {
              if (!mode_is_fan)
                copy3f(vertices[1], vertices[2]);
              copy3f(vertices[0], vertices[1]);
              copy3f(pc, vertices[0]);
            }

            nverts++;
            switch (mode) {
            case GL_TRIANGLES:
              if (!(nverts % 3)) {
                CGOAlphaTriangle(cgo, vertices[2], vertices[1], vertices[0],
                    normals[2], normals[1], normals[0], colors[2], colors[1],
                    colors[0], alpha[2], alpha[1], alpha[0], 0);
              }
              break;
            case GL_TRIANGLE_STRIP:
              if (verticespl < 0) {
                int off0, off2;
                if (flip) {
                  off0 = 0;
                  off2 = 2;
                } else {
                  off0 = 2;
                  off2 = 0;
                }
                flip = !flip;
                CGOAlphaTriangle(cgo, vertices[off0], vertices[1],
                    vertices[off2], normals[off0], normals[1], normals[off2],
                    colors[off0], colors[1], colors[off2], alpha[off0],
                    alpha[1], alpha[off2], 0);
              }
              break;
            case GL_TRIANGLE_FAN:
              if (verticespl < 0) {
                CGOAlphaTriangle(cgo, vertices[2], vertices[1], vertices[0],
                    normals[2], normals[1], normals[0], colors[2], colors[1],
                    colors[0], alpha[2], alpha[1], alpha[0], 0);
              }
            }
            add_to_cgo = !mode_is_triangles;
            hasShifted = 0;
          } else {
            add_to_cgo = 1;
          }
        case CGO_COLOR:
          if (op == CGO_COLOR) {
            add_to_cgo = !mode_is_triangles;
            if (mode_is_triangles) {
              if (colorspl >= 0) {
                copy3f(pc, colors[colorspl]);
                colorspl--;
              } else {
                if (!mode_is_fan)
                  copy3f(colors[1], colors[2]);
                copy3f(colors[0], colors[1]);
                copy3f(pc, colors[0]);
              }
              hasShifted |= 1;
            }
          }
        case CGO_NORMAL:
          if (op == CGO_NORMAL) {
            add_to_cgo = !mode_is_triangles;
            if (mode_is_triangles) {
              if (normalspl >= 0) {
                copy3f(pc, normals[normalspl]);
                normalspl--;
              } else {
                if (!mode_is_fan)
                  copy3f(normals[1], normals[2]);
                copy3f(normals[0], normals[1]);
                copy3f(pc, normals[0]);
              }
              hasShifted |= 2;
            }
          }
        case CGO_ALPHA:
          if (op == CGO_ALPHA) {
            add_to_cgo = !mode_is_triangles;
            if (mode_is_triangles) {
              if (alphapl >= 0)
                alpha[alphapl--] = *pc;
              else {
                if (!mode_is_fan)
                  alpha[2] = alpha[1];
                alpha[1] = alpha[0];
                alpha[0] = *pc;
              }
              hasShifted |= 4;
            }
          }
        default:
          if (add_to_cgo) {
            cgo->add_to_cgo(op, pc);
          }
        }
      }

      if (!mode_is_triangles) {
        CGOEnd(cgo);
      }

      tot_nverts += nverts;
    } break;
    case CGO_COLOR:
      if (op == CGO_COLOR) {
        copy3f(pc, cgo->color);
      }
    case CGO_NORMAL:
      if (op == CGO_NORMAL) {
        copy3f(pc, cgo->normal);
      }
    case CGO_ALPHA:
      if (op == CGO_ALPHA) {
        cgo->alpha = *pc;
      }
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  if (tot_nverts) {
    return (cgo);
  } else {
    CGOFree(cgo);
    return nullptr;
  }
}

/**
 * Converts TRIANGLE(_STRIP|_FAN) to TRIANGLES and generates
 * normals for all triangles. Discards any existing normals for
 * triangles.
 *
 * I: primitive CGO
 * return: new primitive CGO with normals on triangles
 */
CGO* CGOGenerateNormalsForTriangles(const CGO* I)
{
  auto G = I->G;
  auto cgo = CGONewSized(G, I->c);

  float vertices[3][3];
  float current_color[3] = {0.f, 0.f, 0.f}, colors[3][3];
  float current_normal[3];
  float current_alpha = 0, alphas[3];

  bool has_alpha = false;
  bool has_color = false;

  int mode = 0;
  bool inside_begin_triangles = false;
  int current_i = 0;
  int vertex_count = 0;
  bool flip = false;
  bool emit;

  const int indices_regular[] = {0, 1, 2};
  const int indices_flipped[] = {0, 2, 1};

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    auto op = it.op_code();

    if (op == CGO_BEGIN) {
      mode = *reinterpret_cast<const int*>(pc);

      switch (mode) {
      case GL_TRIANGLE_STRIP:
      case GL_TRIANGLE_FAN:
      case GL_TRIANGLES:
        current_i = 0;
        vertex_count = 0;
        flip = false;
        inside_begin_triangles = true;

        CGOBegin(cgo, GL_TRIANGLES);
        continue; // for-loop, no add_to_cgo
      }

      inside_begin_triangles = false;
    } else if (op == CGO_END) {
      inside_begin_triangles = false;
    }

    if (!inside_begin_triangles) {
      cgo->add_to_cgo(op, pc);
      continue;
    }

    // handle operations inside BEGIN/END TRIANGLE(S|_STRIP|_FAN)
    switch (op) {
    case CGO_VERTEX:
      copy3(reinterpret_cast<const float*>(pc), vertices[current_i]);
      copy3(current_color, colors[current_i]);
      alphas[current_i] = current_alpha;

      ++vertex_count;

      switch (mode) {
      case GL_TRIANGLE_STRIP:
        current_i = vertex_count % 3;
        emit = (vertex_count > 2);
        break;
      case GL_TRIANGLE_FAN:
        current_i = ((vertex_count + 1) % 2) + 1;
        emit = (vertex_count > 2);
        break;
      default:
        current_i = vertex_count % 3;
        emit = (current_i == 0);
      }

      if (emit) {
        auto* indices = flip ? indices_flipped : indices_regular;

        if (mode != GL_TRIANGLES) {
          flip = !flip;
        }

        CalculateTriangleNormal(vertices[0], vertices[indices[1]],
            vertices[indices[2]], current_normal);

        CGONormalv(cgo, current_normal);

        for (int j = 0; j < 3; ++j) {
          int k = indices[j];
          if (has_color)
            CGOColorv(cgo, colors[k]);
          if (has_alpha)
            CGOAlpha(cgo, alphas[k]);
          CGOVertexv(cgo, vertices[k]);
        }
      }

      break;
    case CGO_COLOR:
      copy3(reinterpret_cast<const float*>(pc), current_color);
      has_color = true;
      break;
    case CGO_ALPHA:
      current_alpha = *reinterpret_cast<const float*>(pc);
      has_alpha = true;
      break;
    case CGO_NORMAL:
      // discard, we will generate new normals
      break;
    default:
      PRINTFB(G, FB_CGO, FB_Warnings)
      " CGO-Warning: CGOGenerateNormalsForTriangles: unhandled op=0x%02x "
      "inside BEGIN/END\n",
          op ENDFB(G);
      cgo->add_to_cgo(op, pc);
    }
  }

  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  return (cgo);
}

CGO* CGOTurnLightingOnLinesOff(const CGO* I, bool use_shader)
{
  bool cur_mode_is_lines = false;
  auto cgo = CGONewSized(I->G, I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_DRAW_ARRAYS: {
      const cgo::draw::arrays* sp = reinterpret_cast<decltype(sp)>(pc);
      float* vals;
      int nvals = sp->narrays * sp->nverts;
      switch (sp->mode) {
      case GL_LINES:
      case GL_LINE_STRIP:
        CGODisable(cgo, CGO_GL_LIGHTING);
        cur_mode_is_lines = true;
      }
      vals = cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
      memcpy(vals, sp->floatdata, nvals);
      if (cur_mode_is_lines) {
        CGOEnable(cgo, CGO_GL_LIGHTING);
        cur_mode_is_lines = false;
      }
    } break;
    case CGO_DRAW_BUFFERS_INDEXED: {
      const cgo::draw::buffers_indexed* sp = reinterpret_cast<decltype(sp)>(pc);
      int mode = sp->mode, mode_is_lines = 0;
      switch (mode) {
      case GL_LINES:
      case GL_LINE_STRIP:
        mode_is_lines = true;
      }
      if (mode_is_lines) {
        CGODisable(cgo, CGO_GL_LIGHTING);
      }
      cgo->copy_op_from<cgo::draw::buffers_indexed>(pc);
      if (mode_is_lines) {
        CGOEnable(cgo, CGO_GL_LIGHTING);
      }
    } break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED: {
      const cgo::draw::buffers_not_indexed* sp =
          reinterpret_cast<decltype(sp)>(pc);
      int mode = sp->mode, mode_is_lines = 0;
      switch (mode) {
      case GL_LINES:
      case GL_LINE_STRIP:
        mode_is_lines = true;
      }
      if (mode_is_lines) {
        CGODisable(cgo, CGO_GL_LIGHTING);
      }
      cgo->copy_op_from<cgo::draw::buffers_not_indexed>(pc);
      if (mode_is_lines) {
        CGOEnable(cgo, CGO_GL_LIGHTING);
      }
    } break;
    case CGO_END: {
      CGOEnd(cgo);
      if (cur_mode_is_lines) {
        CGOEnable(cgo, CGO_GL_LIGHTING);
        cur_mode_is_lines = 0;
      }
    } break;
    case CGO_BEGIN: {
      int mode = CGO_get_int(pc);
      switch (mode) {
      case GL_LINES:
      case GL_LINE_STRIP:
        CGODisable(cgo, CGO_GL_LIGHTING);
        cur_mode_is_lines = true;
        break;
      default:
        if (!use_shader) { // no shaders, not lines, turn lighting on
          CGOEnable(cgo, CGO_GL_LIGHTING);
        }
      }
      CGOBegin(cgo, mode);
    } break;
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  cgo->use_shader = use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  return (cgo);
}

bool CGOHasAnyTriangleVerticesWithoutNormals(const CGO* I, bool checkTriangles)
{
  bool inside = false;
  bool hasNormal = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_BEGIN:
      switch (CGO_get_int(pc)) {
      case GL_TRIANGLE_FAN:
      case GL_TRIANGLES:
      case GL_TRIANGLE_STRIP:
        if (checkTriangles)
          inside = 1;
        break;
      case GL_LINE_STRIP:
      case GL_LINES:
        if (!checkTriangles)
          inside = 1;
        break;
      }
      break;
    case CGO_END:
      inside = 0;
      break;
    case CGO_NORMAL:
      hasNormal = 1;
      break;
    case CGO_VERTEX:
      if (inside && !hasNormal)
        return 1;
      break;
    case CGO_DRAW_ARRAYS: {
      const auto sp = it.cast<cgo::draw::arrays>();
      switch (sp->mode) {
      case GL_TRIANGLE_FAN:
      case GL_TRIANGLES:
      case GL_TRIANGLE_STRIP:
        if (checkTriangles) {
          if (!(sp->arraybits & CGO_NORMAL_ARRAY)) {
            return 1;
          }
        }
        break;
      case GL_LINE_STRIP:
      case GL_LINES:
        if (!checkTriangles) {
          if (!(sp->arraybits & CGO_NORMAL_ARRAY)) {
            return 1;
          }
        }
        break;
      }
    } break;
    }
  }
  return 0;
}

void CGO::add_to_cgo(int op, const float* pc)
{
  switch (op) {
  case CGO_STOP:
    // only append to buffer, don't increment size
    CGOStop(this);
    break;
  case CGO_DRAW_ARRAYS:
    copy_op_from<cgo::draw::arrays>(pc);
    break;
  case CGO_DRAW_BUFFERS_INDEXED:
    copy_op_from<cgo::draw::buffers_indexed>(pc);
    break;
  case CGO_DRAW_TEXTURES:
    copy_op_from<cgo::draw::textures>(pc);
    break;
  case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
    copy_op_from<cgo::draw::screen_textures>(pc);
    break;
  case CGO_DRAW_LABELS:
    copy_op_from<cgo::draw::labels>(pc);
    break;
  case CGO_DRAW_CONNECTORS:
    copy_op_from<cgo::draw::connectors>(pc);
    break;
  case CGO_DRAW_BUFFERS_NOT_INDEXED:
    copy_op_from<cgo::draw::buffers_not_indexed>(pc);
    break;
  case CGO_DRAW_SPHERE_BUFFERS:
    copy_op_from<cgo::draw::sphere_buffers>(pc);
    break;
  case CGO_DRAW_CYLINDER_BUFFERS:
    copy_op_from<cgo::draw::cylinder_buffers>(pc);
    break;
  case CGO_DRAW_CUSTOM:
    copy_op_from<cgo::draw::custom>(pc);
    break;
  default:
    int sz = CGO_sz[op];
    std::copy_n(pc - 1, sz + 1, add_to_buffer(sz + 1));
  };
}

void CGO::print_table() const
{
  display_table_t table;
  table.begin_row()
      .insert_cell("#")
      .insert_cell("OP_CODE")
      .insert_cell("OP_SZ")
      .insert_cell("DATA");

  unsigned j = 0;
  for (auto it = begin(); !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op_code = it.op_code();

    table.begin_row().insert_cell(++j).insert_cell(op_code).insert_cell(
        CGO_sz[op_code]);
    int sz = CGO_sz[op_code];
    std::stringstream ss;
    for (int i = 0; i < sz; i++) {
      ss << CGO_get_int(pc + i);
      if (i != (sz - 1))
        ss << ", ";
    }
    table.insert_cell(ss.str());
  }
  table.display();
}

CGO* CGOConvertSpheresToPoints(const CGO* I)
{
  CGO* cgo;

  int ok = true;
  cgo = CGONew(I->G);
  CHECKOK(ok, cgo);
  CGOBegin(cgo, GL_POINTS);

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto pc = it.data();
    const auto op = it.op_code();

    switch (op) {
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      CGOPickColor(
          cgo, cgo->current_pick_color_index, cgo->current_pick_color_bond);
      break;
    case CGO_SHADER_CYLINDER:
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
    case CGO_CYLINDER:
    case CGO_CONE:
    case CGO_SAUSAGE:
    case CGO_CUSTOM_CYLINDER:
    case CGO_CUSTOM_CYLINDER_ALPHA:
    case CGO_END:
    case CGO_VERTEX:
    case CGO_BEGIN:
    case CGO_ELLIPSOID:
    case CGO_QUADRIC:
    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    case CGO_DRAW_SPHERE_BUFFERS:
    case CGO_DRAW_CYLINDER_BUFFERS:
    case CGO_DRAW_LABELS:
      break;
    case CGO_SPHERE:
      CGOVertexv(cgo, pc);
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
    default:
      cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  CGOEnd(cgo);
  if (ok) {
    ok &= CGOStop(cgo);
  }
  if (!ok) {
    CGOFree(cgo);
  }
  return (cgo);
}

// Will create an interleaved VBO with { vertex, otherVertex, uv, and texcoord
// info } Currently, this function does not support/parse lines inside
// CGODrawArrays, i.e., a CGO that CGOCombineBeginEnd was used on
CGO* CGOConvertLinesToTrilines(const CGO* I, bool addshaders)
{
  static std::set<int> lineops = {CGO_VERTEX, CGO_LINE, CGO_SPLITLINE};
  auto G = I->G;
  const int nLines = CGOCountNumberOfOperationsOfTypeN(I, lineops) + 1;

  if (nLines == 0) {
    return nullptr;
  }

  int line_counter = 0;
  GLuint glbuff = 0;
  const float* colorv = nullptr;
  unsigned int buff_size = nLines * 6 * (8 * sizeof(float));

  // VBO memory -- number of lines x 4 vertices per line x (vertex + otherVertex
  // + normal + texCoord + color)
  std::vector<float> buffer_start(buff_size);
  float* buffer = buffer_start.data();

  std::unique_ptr<CGO> cgo(new CGO(G));

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_DRAW_ARRAYS: {
      auto sp = it.cast<cgo::draw::arrays>();
      float* vals =
          cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
      int nvals = sp->narrays * sp->nverts;
      memcpy(vals, sp->floatdata, nvals);
    } break;
    case CGO_END:
      WARN_UNEXPECTED_OPERATION(G, op);
      return nullptr;
    case CGO_BEGIN: {
      const float *last_vertex = nullptr, *last_color = nullptr,
                  *current_color = nullptr, *color = nullptr;
      const int mode = it.cast<cgo::draw::begin>()->mode;

      for (++it;; ++it) {
        if (it.is_stop()) {
          WARN_UNEXPECTED_OPERATION(G, CGO_STOP);
          return nullptr;
        }

        const auto op = it.op_code();
        if (op == CGO_END) {
          break;
        }

        const auto pc = it.data();

        switch (op) {
        case CGO_VERTEX:
          if (last_vertex) {
            switch (mode) {
            case GL_LINES:
            case GL_LINE_STRIP: {
              float cc[3] = {1, 1, 1};
              float alpha = cgo->alpha;
              CGOTrilines_GetCurrentColor(
                  current_color, colorv, last_color, cc);
              trilinesBufferAddVertices(
                  buffer, pc, last_vertex, current_color, alpha);
              line_counter++;
              last_vertex = pc;
            }
              if (mode == GL_LINES) {
                last_vertex = nullptr;
                last_color = nullptr;
              }
            }
          } else {
            last_vertex = pc;
            current_color = color;
          }
          break;
        case CGO_LINE: {
          auto line = it.cast<cgo::draw::line>();
          float cc[3] = {1, 1, 1};
          float alpha = cgo->alpha;
          CGOTrilines_GetCurrentColor(current_color, colorv, last_color, cc);
          trilinesBufferAddVertices(
              buffer, line->vertex1, line->vertex2, current_color, alpha);
          line_counter++;
        } break;
        case CGO_SPLITLINE: {
          auto splitline = it.cast<cgo::draw::splitline>();
          float cc[3] = {1, 1, 1};
          float alpha = cgo->alpha;
          float mid[3];
          float color2[] = {CONVERT_COLOR_VALUE(splitline->color2[0]),
              CONVERT_COLOR_VALUE(splitline->color2[1]),
              CONVERT_COLOR_VALUE(splitline->color2[2])};
          add3f(splitline->vertex1, splitline->vertex2, mid);
          mult3f(mid, .5f, mid);
          CGOTrilines_GetCurrentColor(current_color, colorv, last_color, cc);
          trilinesBufferAddVertices(
              buffer, splitline->vertex1, mid, current_color, alpha);
          trilinesBufferAddVertices(
              buffer, mid, splitline->vertex2, color2, alpha);
          line_counter += 2;
        } break;
        case CGO_COLOR:
          last_color = current_color;
          current_color = pc;
          color = pc;
          break;
        }
      }
    } break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_COLOR:
      colorv = pc;
      break;
    }
  }

  cgo->use_shader = I->use_shader;
  if (cgo->use_shader) {
    cgo->cgo_shader_ub_color =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal =
        SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }

  {
    /* TODO: Rmove GL calls and use Generic Vertex Buffers */
    glGenBuffers(1, &glbuff);
    glBindBuffer(GL_ARRAY_BUFFER, glbuff);
    glBufferData(GL_ARRAY_BUFFER, line_counter * 6 * 8 * sizeof(float),
        buffer_start.data(), GL_STATIC_DRAW);
    CheckGLErrorOK(G, "ERROR: CGOConvertLinesToTriangleStrips() glBindBuffer "
                      "returns err=%d\n");
    if (addshaders)
      cgo->add<cgo::draw::enable>(GL_TRILINES_SHADER);
    cgo->add<cgo::draw::trilines>(line_counter * 6, glbuff);
    cgo->has_draw_buffers = true;
    if (addshaders)
      cgo->add<cgo::draw::disable>(GL_TRILINES_SHADER);
    CGOStop(cgo.get());
  }

  return cgo.release();
}

/**
 * copies data for a particular attribute operation into the array used to load
 * the VBO. this takes into account whether it is interleaved or not.
 *
 * isInterleaved    : whether the VBO is interleaved
 * nvert            : which vertex in the VBO
 * attribOp         : the attribute op
 * vertexDataSize   : total vertex data size in VBO (for interleaved)
 * dataPtrs         : all data pointers for attributes (for interleaved, they
 * are all the pointer to the one array) attrOffset       : offsets of the
 * attributes (for interleaved) pcarg            : pc pointer to CGO operation
 * data pick_data        : pointer to pick data for current vertex (writes to if
 * pick data) has_pick_colorBS : keeps track of which pick attributes have been
 * set
 *
 */
static void copyAttributeForOp(bool isInterleaved, int& nvert,
    AttribOp* attribOp, int vertexDataSize, std::vector<void*>& dataPtrs,
    std::vector<int>& attrOffset, const float* pcarg, float* pick_data,
    int& has_pick_colorBS, int pstride)
{
  auto attrDesc = attribOp->desc;
  int ord = attrDesc->order;
  int copyord = -1;
  void* dataPtr = dataPtrs[ord];
  unsigned char* pc = ((unsigned char*) pcarg) + attribOp->offset;
  if (isInterleaved) {
    dataPtr =
        (unsigned char*) dataPtr + nvert * vertexDataSize + attrOffset[ord];
    if (attribOp->copyAttribDesc) {
      copyord = attribOp->copyAttribDesc->order;
      pc = ((unsigned char*) dataPtrs[ord]) + nvert * vertexDataSize +
           attrOffset[copyord];
    }
  } else {
    auto sz = GetSizeOfVertexFormat(attrDesc->m_format);
    dataPtr = (unsigned char*) dataPtr + nvert * sz;
    if (attribOp->copyAttribDesc) {
      copyord = attribOp->copyAttribDesc->order;
      auto copysz = GetSizeOfVertexFormat(attribOp->copyAttribDesc->m_format);
      pc = (unsigned char*) dataPtr + nvert * copysz;
    }
  }
  switch (attribOp->conv_type) {
  case NO_COPY:
    break;
  case FLOAT_TO_FLOAT:
    *((float*) dataPtr) = *((float*) pc);
    break;
  case FLOAT2_TO_FLOAT2:
    *((float*) dataPtr) = *((float*) pc);
    *((float*) dataPtr + 1) = *((float*) pc + 1);
    break;
  case FLOAT3_TO_FLOAT3:
    copy3f((float*) pc, (float*) dataPtr);
    break;
  case FLOAT4_TO_FLOAT4:
    *((float*) dataPtr) = *((float*) pc);
    *((float*) dataPtr + 1) = *((float*) pc + 1);
    *((float*) dataPtr + 2) = *((float*) pc + 2);
    *((float*) dataPtr + 3) = *((float*) pc + 3);
    break;
  case FLOAT3_TO_UB3: {
    auto dataPtrUB = (unsigned char*) dataPtr;
    float* pcf = (float*) pc;
    dataPtrUB[0] = CLIP_COLOR_VALUE(pcf[0]);
    dataPtrUB[1] = CLIP_COLOR_VALUE(pcf[1]);
    dataPtrUB[2] = CLIP_COLOR_VALUE(pcf[2]);
  } break;
  case FLOAT1_TO_UB_4TH: {
    auto dataPtrUB = (unsigned char*) dataPtr;
    float* pcf = (float*) pc;
    dataPtrUB[3] = CLIP_COLOR_VALUE(pcf[0]);
  } break;
  case UB3_TO_UB3: {
    auto dataPtrUB = (unsigned char*) dataPtr;
    auto pcUB = (unsigned char*) pc;
    dataPtrUB[0] = pcUB[0];
    dataPtrUB[1] = pcUB[1];
    dataPtrUB[2] = pcUB[2];
    break;
  }
  case UINT_INT_TO_PICK_DATA: {
    float* pcf = (float*) pc;
    unsigned int index = CGO_get_uint(pcf);
    int bond = CGO_get_int(pcf + 1);
    CGO_put_uint(ord * 2 + pick_data, index);
    CGO_put_int(ord * 2 + pick_data + 1, bond);
    has_pick_colorBS |= (1 << ord);
  } break;
  case FLOAT4_TO_UB4: {
    auto dataPtrUB = (unsigned char*) dataPtr;
    float* pcf = (float*) pc;
    dataPtrUB[0] = CLIP_COLOR_VALUE(pcf[0]);
    dataPtrUB[1] = CLIP_COLOR_VALUE(pcf[1]);
    dataPtrUB[2] = CLIP_COLOR_VALUE(pcf[2]);
    dataPtrUB[3] = CLIP_COLOR_VALUE(pcf[3]);
  } break;
  case CYL_CAP_TO_CAP: {
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] = *pc;
  } break;
  case CYL_CAPS_ARE_ROUND: {
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] = cCylShaderBothCapsRound | cCylShaderInterpColor;
  } break;
  case CYL_CAPS_ARE_FLAT: {
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] = cCylShaderBothCapsFlat | cCylShaderInterpColor;
  } break;
  case CYL_CAPS_ARE_CUSTOM: {
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    float* pcf = (float*) pc;
    const cCylCap cap1 = static_cast<cCylCap>(int(pcf[0]));
    const cCylCap cap2 = static_cast<cCylCap>(int(pcf[1]));
    dataPtrUB[0] =
        cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
  } break;
  case UB1_TO_INTERP: {
    bool interp = (pc[0] & cgo::draw::splitline::interpolation);
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] = interp ? 1 : 0;
  } break;
  case UB1_INTERP_TO_CAP: {
    bool interp = (pc[0] & cgo::draw::splitline::interpolation);
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] =
        (cCylShaderBothCapsRound | (interp ? cCylShaderInterpColor : 0));
  } break;
  case FLOAT1_TO_INTERP: {
    float interp = *((float*) pc);
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] = (interp > .5f) ? 1 : 0;
  } break;
  case FLOAT1_INTERP_TO_CAP: {
    float interp = *((float*) pc);
    unsigned char* dataPtrUB = (unsigned char*) dataPtr;
    dataPtrUB[0] = (cCylShaderBothCapsRound |
                    ((interp > .5f) ? cCylShaderInterpColor : 0));
  } break;
  case UB4_TO_UB4: {
    auto dataPtrUB = (unsigned char*) dataPtr;
    auto pcUB = (unsigned char*) pc;
    dataPtrUB[0] = pcUB[0];
    dataPtrUB[1] = pcUB[1];
    dataPtrUB[2] = pcUB[2];
    dataPtrUB[3] = pcUB[3];
    break;
  }
  case PICK_DATA_TO_PICK_DATA: {
    float* pcf;
    if (copyord < 0) {
      pcf = (float*) pc;
    } else {
      pcf = (copyord * 2 + pick_data);
      if (nvert) {
        pcf -= pstride;
      }
    }
    unsigned int index = CGO_get_uint(pcf);
    int bond = CGO_get_int(pcf + 1);
    CGO_put_uint(ord * 2 + pick_data, index);
    CGO_put_int(ord * 2 + pick_data + 1, bond);
    has_pick_colorBS |= (1 << ord);
    break;
  }
  }
}

/**
 * copies data for a particular attribute into the array used to load the VBO.
 * this takes into account whether it is interleaved or not, and if an attribute
 * has repeat values
 *
 * isInterleaved  : whether the VBO is interleaved
 * nvert          : which vertex in the VBO
 * attribDesc     : the attribute description
 * vertexDataSize : total vertex data size in VBO (for interleaved)
 * dataPtrs       : all data pointers for attributes (for interleaved, they are
 * all the pointer to the one array) attrOffset     : offsets of the attributes
 * (for interleaved)
 *
 */
static void copyAttributeForVertex(bool isInterleaved, int& nvert,
    AttribDesc& attribDesc, const int vertexDataSize,
    std::vector<void*>& dataPtrs, std::vector<int>& attrOffset)
{
  int ord = attribDesc.order;
  void* dataPtr = dataPtrs[ord];
  unsigned char* pc = nullptr;
  auto attrSize = GetSizeOfVertexFormat(attribDesc.m_format);
  if (isInterleaved) {
    dataPtr =
        (unsigned char*) dataPtr + nvert * vertexDataSize + attrOffset[ord];
    pc = (unsigned char*) dataPtr - vertexDataSize;
  } else {
    dataPtr = (unsigned char*) dataPtr + nvert * attrSize;
    pc = (unsigned char*) dataPtr - attrSize;
  }
  if (attribDesc.repeat_value && attribDesc.repeat_value_length) {
    int pos = (nvert % attribDesc.repeat_value_length);
    pc = attribDesc.repeat_value + pos * attrSize;
    memcpy(dataPtr, pc, attrSize);
  } else {
    memcpy(dataPtr, pc, attrSize);
  }
}
/**
 * check all attributes (pick and non-pick) to see if they are specified in the
 * CGO (I) also checks to see if any picking is specified and sets has_picking
 * argument if any of the attributes are not specified and if a default_value is
 * set (in the AttribDesc) then the associated vertex_attribute CGO operation is
 * inserted into the cgo that is passed in.
 *
 * I:           primitive CGO that is processed
 * attrData:    definition of attributes
 * pickData:    definition of pick attributes
 * cgo:         new cgo that could have vertex_attribute CGO operations added
 * has_picking: if there are any operations in the CGO (I) that specifies
 * different values for picking. if has_picking is set, then a VBO for picking
 * is generated in CGOConvertToShader()
 *
 */
static void CheckAttributesForUsage(const CGO* I, AttribDataDesc& attrData,
    AttribDataDesc& pickData, CGO* cgo, bool& has_picking)
{
  size_t attrIdx = 0;

  // need to check attributes:
  //  - remove any that are not needed
  //  - add glVertexAttrib for those that are removed
  std::map<int, int>
      opToAttrUsed; // bitmask for each op to which attributes are set
  for (auto& attrDesc : attrData) {
    auto attrOps = &attrDesc.attrOps;
    attrDesc.order = attrIdx++;
    for (auto attrOpIt = attrOps->begin(); attrOpIt != attrOps->end();
         ++attrOpIt) {
      auto attrOp = &(*attrOpIt);
      if (opToAttrUsed.find(attrOp->op) == opToAttrUsed.end())
        opToAttrUsed[attrOp->op] = 1 << attrDesc.order;
      else
        opToAttrUsed[attrOp->op] |= 1 << attrDesc.order;
    }
  }
  // add picking ops (1 << attrIdx) i.e., any pick op is the last bit
  int pidx = 0;
  for (auto pickDataIt = pickData.begin(); pickDataIt != pickData.end();
       ++pickDataIt) {
    auto pickDesc = &(*pickDataIt);
    auto pickOps = &pickDesc->attrOps;
    pickDesc->order = pidx++;
    for (auto pickOpIt = pickOps->begin(); pickOpIt != pickOps->end();
         ++pickOpIt) {
      auto pickOp = &(*pickOpIt);
      pickOp->desc = pickDesc;
      if (opToAttrUsed.find(pickOp->op) == opToAttrUsed.end())
        opToAttrUsed[pickOp->op] = 1 << attrIdx;
      else
        opToAttrUsed[pickOp->op] |= 1 << attrIdx;
    }
  }
  size_t totAttrIdx = (1 << (attrIdx + 1)) - 1;
  size_t allAttrIdxUsed = 0;
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    if (opToAttrUsed.find(op) != opToAttrUsed.end()) {
      int attrUsed = opToAttrUsed[op];
      if (attrUsed & (1 << attrIdx)) { // if picking, need to check values, and
                                       // take them out if cPickableNoPick
        switch (op) {
        case cgo::draw::shadercylinder2ndcolor::op_code:
          if (reinterpret_cast<const cgo::draw::shadercylinder2ndcolor*>(pc)
                  ->pick_color_bond == cPickableNoPick)
            attrUsed ^= (1 << attrIdx);
          break;
        case cgo::draw::splitline::op_code:
          if (reinterpret_cast<const cgo::draw::splitline*>(pc)->bond ==
              cPickableNoPick)
            attrUsed ^= (1 << attrIdx);
        }
      }
      allAttrIdxUsed |= attrUsed;
      if (allAttrIdxUsed == totAttrIdx)
        break;
    }
  }
  has_picking =
      allAttrIdxUsed & (1 << attrIdx); // has_picking if the last bit is set

  if (allAttrIdxUsed != totAttrIdx) {
    // go through any attributes that aren't used:
    //   - add associated vertex_attribute type (if default_value is set)
    //   - remove attribute from attrData description so that it isn't included
    //   in VBO
    AttribDataDesc attrDataNew;
    for (auto idx = 0; idx < attrIdx; ++idx) {
      if (!attrData[idx].repeat_value && !(allAttrIdxUsed & (1 << idx))) {
        // attribute not used, need to create glVertexAttrib
        if (attrData[idx].default_value) {
          // need to add glVertexAttrib CGO OP
          int attr_lookup_idx =
              I->G->ShaderMgr->GetAttributeUID(attrData[idx].attr_name);
          switch (GetVertexFormatBaseType(attrData[idx].m_format)) {
          case VertexFormatBaseType::Float:
            switch (VertexFormatToGLSize(attrData[idx].m_format)) {
            case 1:
              cgo->add<cgo::draw::vertex_attribute_1f>(
                  attr_lookup_idx, *(float*) attrData[idx].default_value);
              break;
            case 3:
              cgo->add<cgo::draw::vertex_attribute_3f>(
                  attr_lookup_idx, attrData[idx].default_value);
              break;
            default:
              std::cerr << "\tBAD SIZE: attrData[idx].m_format="
                        << static_cast<int>(attrData[idx].m_format)
                        << std::endl;
            }
            break;
          case VertexFormatBaseType::UByte:
            switch (VertexFormatToGLSize(attrData[idx].m_format)) {
            case 1: {
              float val;
              unsigned char valuc = *attrData[idx].default_value;
              if (VertexFormatIsNormalized(attrData[idx].m_format)) {
                val = CLAMPVALUE(valuc / 255.f, 0.f, 1.f);
              } else {
                val = (float) valuc;
              }
              cgo->add<cgo::draw::vertex_attribute_1f>(attr_lookup_idx, val);
            } break;
            case 4:
              cgo->add<cgo::draw::vertex_attribute_4ub>(
                  attr_lookup_idx, attrData[idx].default_value);
              break;
            default:
              std::cerr << "\tNOT IMPLEMENTED: attrData[idx].m_format="
                        << static_cast<int>(attrData[idx].m_format)
                        << std::endl;
            }
          }
        }
      } else {
        attrDataNew.push_back(attrData[idx]);
      }
    }
    attrData.swap(attrDataNew); // only keep attributes that are used
  }
}

/**
 * Populates two structures and sets vertsperpickinfo
 * opToCntPer           : CGO op to how many vertices are generated for each op.
 * opToOrderedAttribOps : CGO op to an ordered map of AttribOps that define how
 * we operate on the attribute arrays for each CGO op. attrData             :
 * definition of attributes pickData             : definition of pick attributes
 * vertsperpickinfo     : number of vertices per pick (only used when
 *
 */
static void PopulateOpsIntoStructuresForConversion(
    std::map<int, int>& opToCntPer,
    std::map<int, std::map<int, AttribOp*>>& opToOrderedAttribOps,
    AttribDataDesc& attrData, AttribDataDesc& pickData, int& vertsperpickinfo,
    const bool has_picking)
{
  size_t attrIdx = 0;
  for (auto& attrDesc : attrData) {
    auto attrOps = &attrDesc.attrOps;
    attrDesc.order = attrIdx++;
    for (auto attrOpIt = attrOps->begin(); attrOpIt != attrOps->end();
         ++attrOpIt) {
      auto attrOp = &(*attrOpIt);
      attrOp->desc = &attrDesc;
      if (attrOp->copyFromAttr >= 0) {
        attrOp->copyAttribDesc = &attrData[attrOp->copyFromAttr];
      }
      if (attrOp->incr_vertices > 0) {
        if (!vertsperpickinfo) {
          vertsperpickinfo = attrOp->incr_vertices;
        } else {
          if (attrOp->incr_vertices != vertsperpickinfo) {
            std::cerr << "WARNING: attrOp->incr_vertices set to multiple "
                         "values, vertsperpickinfo="
                      << vertsperpickinfo
                      << " attrOp->incr_vertices=" << attrOp->incr_vertices
                      << " : picking might get confused" << std::endl;
          }
        }
        if (opToCntPer.find(attrOp->op) == opToCntPer.end())
          opToCntPer[attrOp->op] = attrOp->incr_vertices;
        else
          opToCntPer[attrOp->op] += attrOp->incr_vertices;
      }
      if (opToOrderedAttribOps.find(attrOp->op) == opToOrderedAttribOps.end())
        opToOrderedAttribOps[attrOp->op] = std::map<int, AttribOp*>({});
      opToOrderedAttribOps[attrOp->op][attrOp->order] = attrOp;
    }
  }
  if (has_picking) {
    for (auto pickDataIt = pickData.begin(); pickDataIt != pickData.end();
         ++pickDataIt) {
      auto pickDesc = &(*pickDataIt);
      auto pickOps = &pickDesc->attrOps;
      for (auto pickOpIt = pickOps->begin(); pickOpIt != pickOps->end();
           ++pickOpIt) {
        auto pickOp = &(*pickOpIt);
        if (pickOp->copyFromAttr >= 0) {
          pickOp->copyAttribDesc = &pickData[pickOp->copyFromAttr];
        }
        if (opToOrderedAttribOps.find(pickOp->op) == opToOrderedAttribOps.end())
          opToOrderedAttribOps[pickOp->op] = std::map<int, AttribOp*>({});
        opToOrderedAttribOps[pickOp->op][pickOp->order] = pickOp;
      }
    }
  }
}

/**
 * converts a "primitive" CGO into a CGO that renders a custom operation
 *
 * I:                   primitive CGO that is processed
 * attrData:            definition of attributes that are accumulated and put
 * into the VBO pickData:            definition of pick attributes that
 * accumulate pick data and put into the interleaved picking VBO mode: which
 * openGL mode to use for rendering (e.g., GL_POINTS, GL_TRIANGLES,
 * GL_LINE_STRIPS, etc.) layout:              SEPARATE, SEQUENTIAL, or
 * INTERLEAVED : how the VBO is layed out in memory (default: INTERLEAVED)
 * check_attr_for_data: if true, this function checks whether all attributes are
 * used, and if any are not specified, and default values are defined in
 * AttribDesc, then glVertexAttrib CGO ops are created on the returned CGO
 * (default: true) idx_array:           if specified, this array is used to
 * specify indices for each fragment in an Indexed Buffer and glDrawElements is
 * used to render (instead of glDrawArrays) nvertsperfrag:       in conjunction
 * with idx_array, how many vertices for each set of indices (default: 0) both
 * idx_array and nvertsperfrag are used to specify vertices and geometry for
 * each fragment, such as the box for cylinders nfragspergroup:      Currently,
 * this represents the number of fragments that are inside of a group.  For
 * example, crosses as cylinders (CGOConvertCrossesToCylinderShader) have 36
 * indices per fragment (nvertsperfrag=36) and 3 fragments per group
 * (nfragspergroup=3) (i.e., one for each line). note: idx_array, nvertsperfrag
 * and nfragspergroup should probably be moved to AttrOp in the future to
 * support different types of fragments.
 *
 * returns a CGO that consists of vertex_attrib_* (e.g., glVertexAttrib) and a
 * custom CGO operation that calls either glDrawArrays or glDrawElements. It
 * also supports picking if specified in the pickData.
 *
 */
CGO* CGOConvertToShader(const CGO* I, AttribDataDesc& attrData,
    AttribDataDesc& pickData, int mode, const buffer_layout layout,
    bool check_attr_for_data, int* idx_array, int nvertsperfrag,
    int nfragspergroup)
{
  CGO* cgo;
  int ok = true;
  bool isInterleaved = (layout == buffer_layout::INTERLEAVED);
  bool has_picking = true;
  std::map<std::string, AttribDesc*> attrToDesc;

  cgo = CGONew(I->G);
  cgo->use_shader = true;

  if (check_attr_for_data) {
    CheckAttributesForUsage(I, attrData, pickData, cgo, has_picking);
  } else {
    // if attributes aren't checked, still need to set pick order and desc
    // pointer
    int pidx = 0;
    for (auto& pickDesc : pickData) {
      pickDesc.order = pidx++;
      for (auto& pickOp : pickDesc.attrOps) {
        pickOp.desc = &pickDesc;
      }
    }
  }
  std::map<int, int> opToCntPer;
  std::map<int, std::map<int, AttribOp*>> opToOrderedAttribOps;

  int vertsperpickinfo = 0;
  PopulateOpsIntoStructuresForConversion(opToCntPer, opToOrderedAttribOps,
      attrData, pickData, vertsperpickinfo, has_picking);

  // Populate these variables used for accumulating and setting VBO data arrays
  int vertexDataSize = 0;
  std::vector<int> attrSizes;
  std::vector<int> attrOffset; // for interleaved
  int curoffset = 0;           // for interleaved
  size_t attrIdx = 0;
  for (auto& attrDesc : attrData) {
    attrDesc.order = attrIdx++;
    auto attrSize = GetSizeOfVertexFormat(attrDesc.m_format);
    attrSizes.push_back(attrSize);
    vertexDataSize += attrSize;

    if (isInterleaved) {
      attrOffset.push_back(curoffset);
      curoffset += attrSize;
    } else {
      attrOffset.push_back(0); // no offset when not interleaved
    }
    attrToDesc[attrDesc.attr_name] = &attrDesc;
  }

  // populate funcData.attrib
  for (auto& attrDesc : attrData)
    for (auto& attrop : attrDesc.attrOps)
      for (auto& funcData : attrop.funcDataConversions) {
        if (attrToDesc.find(funcData.attribName) != attrToDesc.end()) {
          funcData.attrib = attrToDesc[funcData.attribName];
        }
      }

  // Since some cards require word-aligned strides (e.g., ATI)
  // we need to make sure it is word-aligned.  If it isn't, and you
  // want to save memory, maybe SEQUENTIAL is a better option
  if (vertexDataSize % 4) {
    vertexDataSize += 4 - (vertexDataSize % 4);
  }

  int ntotalverts = CGOCountNumberOfOperationsOfTypeN(I, opToCntPer);

  // PYMOL-2668
  if (ntotalverts == 0) {
    CGOStop(cgo);
    return cgo;
  }

  size_t pickvbohash = 0;
  int pickDataSize = pickData.size();
  int pstride = 2 * pickDataSize;
  // Generate VBOs: both for vertex data (vbo) and picking color data (pickvbo)
  // if necessary
  // - pickvbo is interleaved so that for multiple channels, pick data for each
  // vertex
  //   is contiguous
  VertexBuffer* vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(layout);
  VertexBuffer* pickvbo = nullptr;
  if (pickDataSize) {
    pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(
        buffer_layout::SEQUENTIAL, GL_DYNAMIC_DRAW);
    pickvbohash = pickvbo->get_hash_id();
  }

  // adding SPECIAL OPERATIONS (for now) before adding custom OP
  // - for now, this is the only operation that needs to be passed to the new
  // CGO
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    switch (op) {
    case CGO_SPECIAL:
      cgo->add_to_cgo(op, pc);
      break;
    }
  }

  // defines how many passes we have
  // We *could* set npickcolattr=1 if we had access to PickColorConverter here
  // and getTotalBits() == 32
  const int npickcolattr = SHADER_PICKING_PASSES_MAX;
  /* for picking, we need to mask the attributes and bind the buffer in the the
   * picking VBO */
  int pl = 0;
  int npickattr = pickData.size();
  for (auto& pickDesc : pickData) {
    cgo->add<cgo::draw::mask_attribute_if_picking>(
        I->G->ShaderMgr->GetAttributeUID(pickDesc.attr_name),
        vbo->get_hash_id());
    if (has_picking) {
      cgo->add<cgo::draw::bind_vbo_for_picking>(pickvbohash, pl++, npickattr);
    } else {
      /* if no picking, should render black */
      static unsigned char zerocolor[]{0, 0, 0, 0};
      cgo->add<cgo::draw::vertex_attribute_4ub_if_picking>(
          I->G->ShaderMgr->GetAttributeUID(pickDesc.attr_name), zerocolor);
    }
  }
  size_t iboid = 0;
  int num_total_indexes = 0;
  if (nvertsperfrag) {
    int nfrags = nfragspergroup * ntotalverts / vertsperpickinfo;
    int nvertsperindivfrag = vertsperpickinfo / nfragspergroup;
    num_total_indexes = nfrags * nvertsperfrag;
    std::vector<VertexIndex_t> vertexIndices(num_total_indexes);
    int idxpl = 0;
    // using vertsperpickinfo as verts per frag
    for (int cnt = 0, vpl = 0; cnt < nfrags; ++cnt) {
      for (int idx_array_pl = 0; idx_array_pl < nvertsperfrag; ++idx_array_pl) {
        vertexIndices[idxpl] = idx_array[idx_array_pl] + vpl;
        idxpl++;
      }
      vpl += nvertsperindivfrag;
    }
    IndexBuffer* ibo = I->G->ShaderMgr->newGPUBuffer<IndexBuffer>();
    ok &= ibo->bufferData({BufferDesc(nullptr, VertexFormat::UInt,
        sizeof(VertexIndex_t) * num_total_indexes, vertexIndices.data())});
    iboid = ibo->get_hash_id();
  }

  // pick_data is interleaved if more than one attribute
  float* pick_data = cgo->add<cgo::draw::custom>(mode, ntotalverts,
      vbo->get_hash_id(), pickvbohash, vertsperpickinfo, pickDataSize, iboid,
      num_total_indexes);
  std::vector<unsigned char> allData(ntotalverts * vertexDataSize);
  std::vector<void*> dataPtrs;
  std::vector<int> repeat_attr_idx;
  int allAttrBS = 0;

  // Initialize first entry in array(s) with default values and populate
  // dataPtrs
  if (isInterleaved) {
    int pl = 0;
    auto attrDataIt = attrData.begin();
    auto attrOffsetIt = attrOffset.begin();
    for (; attrDataIt != attrData.end() && attrOffsetIt != attrOffset.end();
         ++attrDataIt, ++attrOffsetIt) {
      auto attrDesc = &(*attrDataIt);
      if (attrDesc->repeat_value) {
        repeat_attr_idx.push_back(pl);
      } else {
        allAttrBS |= (1 << attrDesc->order);
      }
      auto attrOffset = *attrOffsetIt;
      unsigned char* first_value = nullptr;
      first_value = (attrDesc->default_value
                         ? attrDesc->default_value
                         : (attrDesc->repeat_value ? attrDesc->repeat_value
                                                   : nullptr));
      if (first_value) {
        auto attrSize = GetSizeOfVertexFormat(attrDesc->m_format);
        memcpy(allData.data() + attrOffset, first_value, attrSize);
      }
      dataPtrs.push_back(allData.data());
      ++pl;
    }
  } else {
    auto* curAllDataPtr = allData.data();
    int pl = 0;
    for (auto& attrDesc : attrData) {
      if (attrDesc.repeat_value) {
        repeat_attr_idx.push_back(pl);
      } else {
        allAttrBS |= (1 << attrDesc.order);
      }
      dataPtrs.push_back(curAllDataPtr);
      unsigned char* first_value = nullptr;
      first_value =
          (attrDesc.default_value
                  ? attrDesc.default_value
                  : (attrDesc.repeat_value ? attrDesc.repeat_value : nullptr));
      if (first_value) {
        memcpy(curAllDataPtr, first_value, attrSizes[pl]);
      }
      curAllDataPtr = curAllDataPtr + ntotalverts * attrSizes[pl];
      ++pl;
    }
  }

  int nvert = 0;
  int attrBS = 0;
  int allPickAttrBS = (1 << pickData.size()) - 1;
  int has_pick_colorBS = allPickAttrBS;
  for (int pi = 0; pi < pickDataSize; ++pi) {
    CGO_put_uint(2 * pi + pick_data, 0);
    CGO_put_int(2 * pi + pick_data + 1, cPickableNoPick);
  }
  bool cont =
      true; // need to break while statement as well if past all vertices

  // This is the loop that goes through the CGO and accumulates all of the data
  // for both the rendering and picking VBOs.
  // - For each OP, go through the list of Attribute OPs
  // - Each Attribute OP:
  //    - copies attribute data from the OP
  //    - can generate vertices (incr_vertices)
  // - Attributes are kept track of (attrBS) and when vertices are generated,
  //   the attributes that have not been newly written for the current vertex
  //   are copied from the previous one.
  for (auto it = I->begin(); cont && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    if (opToOrderedAttribOps.find(op) != opToOrderedAttribOps.end()) {
      std::map<int, AttribOp*>* attribOpsInOrder = &opToOrderedAttribOps[op];
      for (auto attribOpIt : *attribOpsInOrder) {
        AttribOp* attribOp = attribOpIt.second;
        int ord = attribOp->desc->order;
        cont = nvert < ntotalverts;
        if (!cont)
          break;
        copyAttributeForOp(isInterleaved, nvert, attribOp, vertexDataSize,
            dataPtrs, attrOffset, pc, pick_data, has_pick_colorBS, pstride);
        if (ord >= 0) // picking is negative, has_pick_colorBS is used instead
          attrBS |= (1 << ord);
        else
          std::cout << "   ord=%d\n" << ord << std::endl;
        if (attribOp->incr_vertices) {
          if (has_pick_colorBS != allPickAttrBS) {
            // copy pick colors that haven't been set from previous vertex
            for (int pi = 0; pi < pickDataSize; ++pi) {
              if (has_pick_colorBS ^ (1 << pi)) {
                CGO_put_uint(pick_data, CGO_get_uint(pick_data - pstride));
                CGO_put_int(
                    pick_data + 1, CGO_get_int(pick_data - pstride + 1));
              }
              pick_data += 2;
            }
          } else {
            pick_data += pstride;
          }
          has_pick_colorBS = 0;
          if (!nvert && attrBS != allAttrBS) {
            // for the first vertex, all attributes should be set
            for (auto idx = 0; idx < attrData.size(); ++idx) {
              if (!(attrBS & (1 << idx))) {
                if (!attrData[idx].default_value) {
                  std::cerr << "WARNING: attribute #" << idx << " ("
                            << attrData[idx].attr_name
                            << ") not set for first"
                               " vertex and does not have default value"
                            << std::endl;
                }
              }
            }
          }
          if (nvert && attrBS != allAttrBS) {
            // for each vertex that hasn't been written for the current vertex,
            // copy it from the previous vertex
            for (auto idx = 0; idx < attrData.size(); ++idx) {
              if (!(attrBS & (1 << idx))) {
                copyAttributeForVertex(isInterleaved, nvert, attrData[idx],
                    vertexDataSize, dataPtrs, attrOffset);
              }
            }
          }
          attrBS = 0;

          // creating new vertices, all attribute data should be copied into new
          // vertex.
          for (int nxt = 0; nxt < attribOp->incr_vertices; ++nxt) {
            /* for now, always copy values from previous */
            ++nvert;
            if (nvert < ntotalverts) {
              {
                // last should not need to copy into next, since we call
                // copyAttributeForVertex (above) for all attributes that
                // haven't been set
                // - note: it might be faster (especially for interleaved) to
                // copy into the next
                //         vertex, then the above copyAttributeForVertex() would
                //         not be needed
                if (isInterleaved) {
                  auto* dest = allData.data() + vertexDataSize * nvert;
                  memcpy(dest, dest - vertexDataSize, vertexDataSize);
                } else {
                  auto dataPtrIt = dataPtrs.begin();
                  auto attrDataIt = attrData.begin();
                  for (; attrDataIt != attrData.end() &&
                         dataPtrIt != dataPtrs.end();
                       ++attrDataIt, ++dataPtrIt) {
                    auto attrDesc = &(*attrDataIt);
                    auto dataPtr = *dataPtrIt;
                    auto attrSize = GetSizeOfVertexFormat(attrDesc->m_format);
                    void* dest = ((unsigned char*) dataPtr) + attrSize * nvert;
                    memcpy((unsigned char*) dest,
                        ((unsigned char*) dest) - attrSize, attrSize);
                  }
                }
              }
              // always copy repeat attributes
              if (!repeat_attr_idx.empty()) {
                for (auto ridx = repeat_attr_idx.begin();
                     ridx != repeat_attr_idx.end(); ++ridx) {
                  copyAttributeForVertex(isInterleaved, nvert, attrData[*ridx],
                      vertexDataSize, dataPtrs, attrOffset);
                }
              }
            }
            if (!attribOp->funcDataConversions.empty()) {
              // for all attributes, call the funcDataConversion() if defined
              int nvert_m_1 = nvert - 1;
              for (auto funcData : attribOp->funcDataConversions) {
                auto funcAttrib = funcData.attrib;
                auto order = funcAttrib->order;
                if (isInterleaved) {
                  unsigned char* dest = allData.data() +
                                        vertexDataSize * nvert_m_1 +
                                        attrOffset[order];
                  funcData.funcDataConversion(
                      dest, pc, funcData.funcDataGlobalArg, nxt);
                } else {
                  auto dataPtr = dataPtrs[order];
                  void* dest =
                      ((unsigned char*) dataPtr) + attrSizes[order] * nvert_m_1;
                  funcData.funcDataConversion(
                      dest, pc, funcData.funcDataGlobalArg, nxt);
                }
              }
            }
          }
        }
      }
    }
  }

  /* Generate Pick Buffers with all pick attributes (if necessary) */
  if (pickvbo) {
    BufferDataDesc pickBufferData;
    for (int i = 0; i < npickcolattr; i++) {
      for (auto& pickDesc : pickData) {
        auto pickSize = GetSizeOfVertexFormat(pickDesc.m_format);
        pickBufferData.push_back(BufferDesc(
            pickDesc.attr_name, pickDesc.m_format, pickSize * nvert));
      }
    }
    pickvbo->bufferData(std::move(pickBufferData));
  }

  /* Generate VBO Buffers with all pick attributes based on the VertexBuffer
   * type SEPARATE/SEQUENTIAL/INTERLEAVED*/
  BufferDataDesc bufferData;
  switch (layout) {
  case buffer_layout::SEPARATE:
  case buffer_layout::SEQUENTIAL: {
    auto attrDataIt = attrData.begin();
    auto dataPtrIt = dataPtrs.begin();
    auto attrSizeIt = attrSizes.begin();
    for (; attrDataIt != attrData.end() && dataPtrIt != dataPtrs.end() &&
           attrSizeIt != attrSizes.end();
         ++attrDataIt, ++dataPtrIt, ++attrSizeIt) {
      auto attrDesc = &(*attrDataIt);
      auto dataPtr = *dataPtrIt;
      auto attrSize = *attrSizeIt;
      bufferData.push_back(BufferDesc(
          attrDesc->attr_name, attrDesc->m_format, nvert * attrSize, dataPtr));
    }
    vbo->bufferData(std::move(bufferData));
    break;
  } break;
  case buffer_layout::INTERLEAVED: {
    auto attrDataIt = attrData.begin();
    auto attrOffsetIt = attrOffset.begin();
    for (; attrDataIt != attrData.end() && attrOffsetIt != attrOffset.end();
         ++attrDataIt, ++attrOffsetIt) {
      auto attrDesc = &(*attrDataIt);
      auto offset = *attrOffsetIt;
      bufferData.push_back(BufferDesc{attrDesc->attr_name, attrDesc->m_format,
          0, nullptr, (std::uint32_t) offset});
    }
    vbo->bufferData(std::move(bufferData), allData.data(),
        (size_t) (nvert * vertexDataSize), (size_t) vertexDataSize);
    break;
  }
  }

  CGOStop(cgo);
  return cgo;
}

// CGOCheckSplitLineInterpolationIsSame:
//   - returns true if always the same
//   - returns false if not always the same
bool CGOCheckSplitLineInterpolationIsSame(const CGO* I, bool& interp_value)
{
  bool interp_value_first = false;
  bool interp_value_is_set = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
    case cgo::draw::splitline::op_code:
      interp_value = (it.cast<cgo::draw::splitline>()->flags &
                      cgo::draw::splitline::interpolation);
      break;
    case CGO_INTERPOLATED:
      interp_value = it.cast<float>()[0] > 0.5f;
      break;
    default:
      continue;
    }
    if (!interp_value_is_set) {
      interp_value_first = interp_value;
      interp_value_is_set = true;
    } else if (interp_value != interp_value_first) {
      return false;
    }
  }
  return true;
}

// CGOCheckShaderCylinderCapInfoIsSame:
//   - returns true if always the same
//   - returns false if not always the same
static bool CGOCheckShaderCylinderCapInfoIsSame(
    const CGO* I, unsigned char& cap_value)
{
  unsigned char cap_value_first = 0;
  bool cap_value_first_is_set = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
    case cgo::draw::shadercylinder::op_code:
      cap_value = it.cast<cgo::draw::shadercylinder>()->cap;
      break;
    case cgo::draw::shadercylinder2ndcolor::op_code:
      cap_value = it.cast<cgo::draw::shadercylinder2ndcolor>()->cap;
      break;
    case cgo::draw::sausage::op_code:
      cap_value = cCylShaderBothCapsRound | cCylShaderInterpColor;
      break;
    case cgo::draw::cylinder::op_code:
      cap_value = cCylShaderBothCapsFlat | cCylShaderInterpColor;
      break;
    case cgo::draw::custom_cylinder::op_code: {
      auto cc = it.cast<cgo::draw::custom_cylinder>();
      const cCylCap cap1 = cc->get_cap1();
      const cCylCap cap2 = cc->get_cap2();
      cap_value = cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
    } break;
    case cgo::draw::custom_cylinder_alpha::op_code: {
      auto cc = it.cast<cgo::draw::custom_cylinder_alpha>();
      const cCylCap cap1 = cc->get_cap1();
      const cCylCap cap2 = cc->get_cap2();
      cap_value = cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
    } break;
    default:
      continue;
    }

    if (!cap_value_first_is_set) {
      cap_value_first = cap_value;
      cap_value_first_is_set = true;
    } else if (cap_value != cap_value_first) {
      return false;
    }
  }

  return true;
}

CGO* CGOConvertToTrilinesShader(const CGO* I, CGO* addTo, bool add_color)
{
  PyMOLGlobals* G = I->G;

  AttribDataOp vertexOps = {
      {CGO_LINE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex1), 0},
      {CGO_SPLITLINE, 2, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::splitline, vertex1), 0}};
  AttribDataOp vertexOtherOps = {
      {CGO_LINE, 2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex2), 6},
      {CGO_SPLITLINE, 5, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::splitline, vertex2), 6}};
  AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 6, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  AttribDataOp color2Ops = {{CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 3, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 7, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::splitline, index), 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 4, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::splitline, index), 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  AttribDataDesc attrDesc = {{"a_Vertex", VertexFormat::Float3, vertexOps},
      {"a_OtherVertex", VertexFormat::Float3, vertexOtherOps},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"a_UV", VertexFormat::UByte}};

  if (add_color) {
    static unsigned char default_color[] = {
        255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
    attrDesc[2].default_value = default_color;
    attrDesc[3].default_value = default_color;
  }
  AttribDesc* uvdesc = &attrDesc[attrDesc.size() - 1];
  uvdesc->repeat_value_length = 6;
  static unsigned char uv_bits[] = {1, 3, 0, 3, 2, 1};
  uvdesc->repeat_value = uv_bits;

  bool interp_same, interp_value;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_interpolate"),
        interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_TO_INTERP,
        offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_interpolate attribute
    attrDesc.push_back({"a_interpolate", VertexFormat::UByte, interpOps});
  }
  if (!add_color) {
    attrDesc.erase(attrDesc.begin() + 2); // a_Color
    attrDesc.erase(attrDesc.begin() + 2); // a_Color2
  }

  return CGOConvertToShader(
      I, attrDesc, pickDesc, GL_TRIANGLES, buffer_layout::INTERLEAVED);
}

CGO* CGOConvertToLinesShader(const CGO* I, CGO* addTo, bool add_color)
{
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;
  AttribDataOp vertexOps = {
      {CGO_LINE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex1), 1},
      {CGO_SPLITLINE, 2, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::splitline, vertex1), 1},
      {CGO_LINE, 2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex2), 1},
      {CGO_SPLITLINE, 5, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::splitline, vertex2), 1}};
  AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 3, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 4, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::splitline, index), 0}};

  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps}};
  AttribDataDesc attrDesc = {{"a_Vertex", VertexFormat::Float3, vertexOps},
      {"a_Color", VertexFormat::UByte4Norm, colorOps}};
  if (add_color) {
    static unsigned char default_color[] = {
        255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
    attrDesc[1].default_value = default_color;
  }
  bool interp_same, interp_value;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_interpolate"),
        interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_TO_INTERP,
        offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_interpolate attribute
    attrDesc.push_back({"a_interpolate", VertexFormat::UByte, interpOps});
  }
#ifndef PURE_OPENGL_ES_2
  {
    attrDesc.push_back({"a_line_position", VertexFormat::UByte});
    AttribDesc* lpdesc = &attrDesc[attrDesc.size() - 1];
    lpdesc->repeat_value_length = 2;
    static unsigned char flip_bits[] = {0, 1};
    lpdesc->repeat_value = flip_bits;
  }
#endif
  if (!add_color) {
    attrDesc.erase(attrDesc.begin() + 1); // a_Color
  }

  return CGOConvertToShader(
      I, attrDesc, pickDesc, GL_LINES, buffer_layout::INTERLEAVED);
}

CGO* CGOConvertLinesToCylinderShader(const CGO* I, CGO* addTo, bool add_color)
{
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;

  AttribDataOp vertex1Ops = {
      {CGO_LINE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex1), 0},
      {CGO_SPLITLINE, 2, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::splitline, vertex1), 0}};
  AttribDataOp vertex2Ops = {
      {CGO_LINE, 2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::line, vertex2), 8},
      {CGO_SPLITLINE, 5, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::splitline, vertex2), 8}};
  static AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 6, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};
  static AttribDataOp color2Ops = {{CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0},
      {CGO_SPLITLINE, 3, UB3_TO_UB3, offsetof(cgo::draw::splitline, color2)}};

  AttribDataDesc attrDesc = {{"attr_vertex1", VertexFormat::Float3, vertex1Ops},
      {"attr_vertex2", VertexFormat::Float3, vertex2Ops},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"attr_radius", VertexFormat::Float}};

  AttribDesc* fdesc;
  static unsigned char cyl_flags[] = {
      0, 4, 6, 2, 1, 5, 7, 3}; // right(4)/up(2)/out(1)

  attrDesc.push_back({"attr_flags", VertexFormat::UByte});
  fdesc = &attrDesc[attrDesc.size() - 1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  if (add_color) {
    static unsigned char default_color[] = {
        255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
    fdesc = &attrDesc[2];
    fdesc->default_value = default_color;
    fdesc = &attrDesc[3];
    fdesc->default_value = default_color;
  }

  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*) &default_radius;

  int box_indices[36] = {// box indices
      0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 0, 7, 3, 7, 0, 4, 3,
      6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4};
  int* box_indices_ptr = nullptr;
  box_indices_ptr = box_indices;

  bool interp_same, interp_value = false;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_cap"),
        (cCylShaderBothCapsRound | (interp_value ? cCylShaderInterpColor : 0)));
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_INTERP_TO_CAP,
        offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_cap attribute
    attrDesc.push_back({"a_cap", VertexFormat::UByte, interpOps});
  }

  if (!add_color) {
    attrDesc.erase(attrDesc.begin() + 2); // attr_colors
    attrDesc.erase(attrDesc.begin() + 2); // attr_colors2
  }

  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 7, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::splitline, index), 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SPLITLINE, 4, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::splitline, index), 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
      buffer_layout::INTERLEAVED, true, box_indices_ptr, 36);
}

struct CrossSizeData {
  float cross_size;
  bool forward;
  CrossSizeData(float _cross_size, bool _forward)
      : cross_size(_cross_size)
      , forward(_forward)
  {
  }
};

static void CrossVertexConversion(
    void* varData, const float* pc, void* crossData, int idx)
{
  CrossSizeData* csd = (CrossSizeData*) crossData;
  int idxpl = idx / 8; // X Y or Z
  float* varDataF = ((float*) varData);
  varDataF[idxpl] += (csd->forward ? csd->cross_size : -csd->cross_size);
}

CGO* CGOConvertCrossesToCylinderShader(
    const CGO* I, CGO* addTo, float cross_size_arg)
{
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;
  AttribDataOp vertex1Ops = {{CGO_VERTEX_CROSS, 1, FLOAT3_TO_FLOAT3, 0, 0}};
  AttribDataOp vertex2Ops = {
      {CGO_VERTEX_CROSS, 2, FLOAT3_TO_FLOAT3, 0, 3 * 8, 0}};

  static AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0}, {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0}};
  static AttribDataOp color2Ops = {
      {CGO_COLOR, 1, FLOAT3_TO_UB3, 0}, {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0}};

  CrossSizeData crossData[] = {{cross_size_arg, false}, {cross_size_arg, true}};
  AttribDataDesc attrDesc = {{"attr_vertex1", VertexFormat::Float3, vertex1Ops},
      {"attr_vertex2", VertexFormat::Float3, vertex2Ops},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"attr_radius", VertexFormat::Float}};

  attrDesc.reserve(10);
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversion, &crossData[0], "attr_vertex1"});
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversion, &crossData[1], "attr_vertex2"});

  AttribDesc* fdesc;
  static unsigned char cyl_flags[] = {
      0, 4, 6, 2, 1, 5, 7, 3}; // right(4)/up(2)/out(1)

  attrDesc.push_back({"attr_flags", VertexFormat::UByte});
  fdesc = &attrDesc[attrDesc.size() - 1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  unsigned char default_color[] = {
      255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
  fdesc = &attrDesc[2];
  fdesc->default_value = default_color;
  fdesc = &attrDesc[3];
  fdesc->default_value = default_color;
  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*) &default_radius;

  int box_indices[36] = {// box indices
      0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 0, 7, 3, 7, 0, 4, 3,
      6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4};
  int* box_indices_ptr = nullptr;
  box_indices_ptr = box_indices;

  addTo->add<cgo::draw::vertex_attribute_1f>(
      G->ShaderMgr->GetAttributeUID("a_cap"), cCylShaderBothCapsRound);

  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
      buffer_layout::INTERLEAVED, true, box_indices_ptr, 36, 3);
}

struct CrossSizeDataLines {
  float cross_size;
  CrossSizeDataLines(float _cross_size)
      : cross_size(_cross_size)
  {
  }
};

static void CrossVertexConversionLines(
    void* varData, const float* pc, void* crossData, int idx)
{
  CrossSizeDataLines* csd = (CrossSizeDataLines*) crossData;
  int idxpl = idx / 2; // X Y or Z
  bool forward = idx % 2;
  float* varDataF = ((float*) varData);
  varDataF[idxpl] += (forward ? csd->cross_size : -csd->cross_size);
}

CGO* CGOConvertCrossesToLinesShader(
    const CGO* I, CGO* addTo, float cross_size_arg)
{
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;
  AttribDataOp vertexOps = {
      {CGO_VERTEX_CROSS, 1, FLOAT3_TO_FLOAT3, 0, 6}}; // 6 vertices for a cross
  AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0}, {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};

  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps}};
  AttribDataDesc attrDesc = {{"a_Vertex", VertexFormat::Float3, vertexOps},
      {"a_Color", VertexFormat::UByte4Norm, colorOps}};
  unsigned char default_color[] = {
      255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
  attrDesc[1].default_value = default_color;

  CrossSizeDataLines crossData = {cross_size_arg};
  attrDesc[0].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversionLines, &crossData, "a_Vertex"});

  bool interp_same, interp_value = false;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_interpolate"),
        interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps = {{CGO_SPLITLINE, 1, UB1_TO_INTERP,
        offsetof(cgo::draw::splitline, flags), 0}};
    // need to add a_interpolate attribute
    attrDesc.push_back({"a_interpolate", VertexFormat::UByte, interpOps});
  }
#ifndef PURE_OPENGL_ES_2
  {
    attrDesc.push_back({"a_line_position", VertexFormat::UByte});
    AttribDesc* lpdesc = &attrDesc[attrDesc.size() - 1];
    lpdesc->repeat_value_length = 2;
    static unsigned char flip_bits[] = {0, 1};
    lpdesc->repeat_value = flip_bits;
  }
#endif
  return CGOConvertToShader(
      I, attrDesc, pickDesc, GL_LINES, buffer_layout::INTERLEAVED);
}

static void CrossVertexConversionTrilines(
    void* varData, const float* pc, void* crossData, int idx)
{
  CrossSizeData* csd = (CrossSizeData*) crossData;
  int idxpl = idx / 6; // X Y or Z
  float* varDataF = ((float*) varData);
  varDataF[idxpl] += (csd->forward ? csd->cross_size : -csd->cross_size);
}

CGO* CGOConvertCrossesToTrilinesShader(
    const CGO* I, CGO* addTo, float cross_size_arg)
{
  PyMOLGlobals* G = I->G;

  AttribDataOp vertexOps = {{CGO_VERTEX_CROSS, 1, FLOAT3_TO_FLOAT3, 0, 0}};
  AttribDataOp vertexOtherOps = {
      {CGO_VERTEX_CROSS, 2, FLOAT3_TO_FLOAT3, 0, 6 * 3}};
  AttribDataOp colorOps = {
      {CGO_COLOR, 0, FLOAT3_TO_UB3, 0}, {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0}};
  AttribDataOp color2Ops = {
      {CGO_COLOR, 1, FLOAT3_TO_UB3, 0}, {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0}};
  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  AttribDataDesc attrDesc = {{"a_Vertex", VertexFormat::Float3, vertexOps},
      {"a_OtherVertex", VertexFormat::Float3, vertexOtherOps},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"a_UV", VertexFormat::UByte}};

  CrossSizeData crossData[] = {{cross_size_arg, false}, {cross_size_arg, true}};

  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversionTrilines, &crossData[0], "a_Vertex"});
  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {CrossVertexConversionTrilines, &crossData[1], "a_OtherVertex"});

  unsigned char default_color[] = {
      255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
  attrDesc[2].default_value = default_color;
  attrDesc[3].default_value = default_color;

  AttribDesc* uvdesc = &attrDesc[attrDesc.size() - 1];
  uvdesc->repeat_value_length = 6;
  static unsigned char uv_bits[] = {1, 3, 0, 3, 2, 1};
  uvdesc->repeat_value = uv_bits;

  addTo->add<cgo::draw::vertex_attribute_1f>(
      G->ShaderMgr->GetAttributeUID("a_interpolate"), 0.f);

  return CGOConvertToShader(
      I, attrDesc, pickDesc, GL_TRIANGLES, buffer_layout::INTERLEAVED);
}

cgo::draw::shadercylinder2ndcolor::shadercylinder2ndcolor(CGO* I,
    const float* _origin, const float* _axis, const float _tube_size, int _cap,
    const float* _color2, Pickable* pickcolor2, const float _alpha)
    : tube_size(_tube_size)
    , alpha(_alpha)
{
  copy3f(_origin, origin);
  copy3f(_axis, axis);
  cap = _cap;
  copy3f(_color2, color2);
  if (pickcolor2) {
    I->current_pick_color_index = pick_color_index = pickcolor2->index;
    I->current_pick_color_bond = pick_color_bond = pickcolor2->bond;
  } else {
    pick_color_index = I->current_pick_color_index;
    pick_color_bond = I->current_pick_color_bond;
  }
};

static void SetVertexFromOriginAxisForCylinder(
    void* varData, const float* pc, void* np, int idx)
{
  float* varDataF = ((float*) varData);
  add3f(pc, pc + 3, varDataF); // adding origin and axis for both shadercylinder
                               // and shadercylinder2ndcolor
}

/**
 * converts all cylinders in the input CGO to a CGO custom operation, which
 * includes picking information (if it exists)
 *
 * I     - input CGO (includes cylinders)
 * addTo - CGO that vertex_attribute operations are added to (if needed), in
 * this case for caps if values for all cylinders are the same
 *
 */
CGO* CGOConvertShaderCylindersToCylinderShader(const CGO* I, CGO* addTo)
{
  /* Lines that pass in two vertices per line */
  PyMOLGlobals* G = I->G;

  // TODO: NEED TO ADD: CGO_CUSTOM_CYLINDER and CGO_CYLINDER
  AttribDataOp vertex1Ops = {
      {CGO_SHADER_CYLINDER, 1, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::shadercylinder, origin), 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 1, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::shadercylinder2ndcolor, origin), 0},
      {CGO_SAUSAGE, 1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::sausage, vertex1),
          0},
      {CGO_CYLINDER, 1, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::cylinder, vertex1), 0},
      {CGO_CUSTOM_CYLINDER, 1, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::custom_cylinder, vertex1), 0},
      {CGO_CUSTOM_CYLINDER_ALPHA, 1, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::custom_cylinder_alpha, vertex1), 0}};
  AttribDataOp vertex2Ops = {{CGO_SHADER_CYLINDER, 5, FLOAT3_TO_FLOAT3,
                                 offsetof(cgo::draw::shadercylinder, axis), 8},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 6, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::shadercylinder2ndcolor, axis), 8},
      {CGO_SAUSAGE, 6, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::sausage, vertex2),
          8},
      {CGO_CYLINDER, 6, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::cylinder, vertex2), 8},
      {CGO_CUSTOM_CYLINDER, 6, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::custom_cylinder, vertex2), 8},
      {CGO_CUSTOM_CYLINDER_ALPHA, 6, FLOAT3_TO_FLOAT3,
          offsetof(cgo::draw::custom_cylinder_alpha, vertex2), 8}};
  static AttribDataOp colorOps = {{CGO_COLOR, 0, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 0, FLOAT1_TO_UB_4TH, 0},
      {CGO_SAUSAGE, 4, FLOAT3_TO_UB3, offsetof(cgo::draw::sausage, color1)},
      {CGO_CYLINDER, 4, FLOAT3_TO_UB3, offsetof(cgo::draw::cylinder, color1)},
      {CGO_CUSTOM_CYLINDER, 4, FLOAT3_TO_UB3,
          offsetof(cgo::draw::custom_cylinder, color1)},
      {CGO_CUSTOM_CYLINDER_ALPHA, 4, FLOAT4_TO_UB4,
          offsetof(cgo::draw::custom_cylinder_alpha, color1)}};
  static AttribDataOp color2Ops = {{CGO_COLOR, 1, FLOAT3_TO_UB3, 0},
      {CGO_ALPHA, 1, FLOAT1_TO_UB_4TH, 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 2, FLOAT3_TO_UB3,
          offsetof(cgo::draw::shadercylinder2ndcolor, color2)},
      {CGO_SAUSAGE, 5, FLOAT3_TO_UB3, offsetof(cgo::draw::sausage, color2)},
      {CGO_CYLINDER, 5, FLOAT3_TO_UB3, offsetof(cgo::draw::cylinder, color2)},
      {CGO_CUSTOM_CYLINDER, 5, FLOAT3_TO_UB3,
          offsetof(cgo::draw::custom_cylinder, color2)},
      {CGO_CUSTOM_CYLINDER_ALPHA, 5, FLOAT4_TO_UB4,
          offsetof(cgo::draw::custom_cylinder_alpha, color2)}};
  AttribDataOp radiusOps = {
      {CGO_SHADER_CYLINDER, 2, FLOAT_TO_FLOAT,
          offsetof(cgo::draw::shadercylinder, tube_size), 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 3, FLOAT_TO_FLOAT,
          offsetof(cgo::draw::shadercylinder2ndcolor, tube_size), 0},
      {CGO_SAUSAGE, 3, FLOAT_TO_FLOAT, offsetof(cgo::draw::sausage, radius), 0},
      {CGO_CYLINDER, 3, FLOAT_TO_FLOAT, offsetof(cgo::draw::cylinder, radius),
          0},
      {CGO_CUSTOM_CYLINDER, 3, FLOAT_TO_FLOAT,
          offsetof(cgo::draw::custom_cylinder, radius), 0},
      {CGO_CUSTOM_CYLINDER_ALPHA, 3, FLOAT_TO_FLOAT,
          offsetof(cgo::draw::custom_cylinder_alpha, radius), 0}};

  AttribDataDesc attrDesc = {{"attr_vertex1", VertexFormat::Float3, vertex1Ops},
      {"attr_vertex2", VertexFormat::Float3, vertex2Ops},
      {"a_Color", VertexFormat::UByte4Norm, colorOps},
      {"a_Color2", VertexFormat::UByte4Norm, color2Ops},
      {"attr_radius", VertexFormat::Float, radiusOps}};

  AttribDesc* fdesc;
  static unsigned char cyl_flags[] = {
      0, 4, 6, 2, 1, 5, 7, 3}; // right(4)/up(2)/out(1)

  attrDesc[1].attrOps[0].funcDataConversions.push_back(
      {SetVertexFromOriginAxisForCylinder, nullptr, "attr_vertex2"});
  attrDesc[1].attrOps[1].funcDataConversions.push_back(
      {SetVertexFromOriginAxisForCylinder, nullptr, "attr_vertex2"});

  attrDesc.push_back({"attr_flags", VertexFormat::UByte});
  fdesc = &attrDesc[attrDesc.size() - 1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  unsigned char default_color[] = {
      255, 255, 255, 255}; // to write in alpha if CGO doesn't have it
  fdesc = &attrDesc[2];
  fdesc->default_value = default_color;
  fdesc = &attrDesc[3];
  fdesc->default_value = default_color;
  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*) &default_radius;

  int box_indices[36] = {// box indices
      0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 0, 7, 3, 7, 0, 4, 3,
      6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4};
  int* box_indices_ptr = nullptr;
  box_indices_ptr = box_indices;

  bool interp_same;
  unsigned char cap_value = 0;
  if ((interp_same = CGOCheckShaderCylinderCapInfoIsSame(I, cap_value))) {
    addTo->add<cgo::draw::vertex_attribute_1f>(
        G->ShaderMgr->GetAttributeUID("a_cap"), cap_value);
  } else {
    AttribDataOp interpOps = {
        {CGO_SHADER_CYLINDER, 3, CYL_CAP_TO_CAP,
            offsetof(cgo::draw::shadercylinder, cap), 0},
        {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 4, CYL_CAP_TO_CAP,
            offsetof(cgo::draw::shadercylinder2ndcolor, cap), 0},
        {CGO_SAUSAGE, 2, CYL_CAPS_ARE_ROUND, 0, 0},
        {CGO_CYLINDER, 2, CYL_CAPS_ARE_FLAT, 0, 0},
        {CGO_CUSTOM_CYLINDER, 2, CYL_CAPS_ARE_CUSTOM,
            offsetof(cgo::draw::custom_cylinder, cap1), 0},
        {CGO_CUSTOM_CYLINDER_ALPHA, 2, CYL_CAPS_ARE_CUSTOM,
            offsetof(cgo::draw::custom_cylinder_alpha, cap1), 0},
    };
    attrDesc.push_back({"a_cap", VertexFormat::UByte, interpOps});
  }

  AttribDataOp extraPickColorOps = {
      {CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 8, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::shadercylinder2ndcolor, pick_color_index), 0}};
  AttribDataOp extraPickColor2Ops = {
      {CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0},
      {CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 5, UINT_INT_TO_PICK_DATA,
          offsetof(cgo::draw::shadercylinder2ndcolor, pick_color_index), 0}};
  AttribDataDesc pickDesc = {
      {"a_Color", VertexFormat::UByte4Norm, extraPickColorOps},
      {"a_Color2", VertexFormat::UByte4Norm, extraPickColor2Ops}};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES,
      buffer_layout::INTERLEAVED, true, box_indices_ptr, 36);
}

/**
 * CGO iterator increment
 */
CGO::const_iterator& CGO::const_iterator::operator++()
{
  const unsigned op = op_code();

  // Corrupted OP codes should never make it into a CGO, so we don't want to
  // handle this gracefully. Only import of CGOs (from PSEs, from Python, etc.)
  // should handle invalid codes gracefully.
  assert(op < CGO_sz_size());

  m_pc += CGO_sz[op] + 1;
  return *this;
}

void CGORender(CGO* I, const float* color, CSetting* set1, CSetting* set2,
    RenderInfo* info, Rep* rep)
{
  CGORenderGL(I, color, set1, set2, info, rep);
}

void CGORenderPicking(CGO* I, RenderInfo* info, PickContext* context,
    CSetting* set1, CSetting* set2, Rep* rep)
{
  CGORenderGLPicking(I, info, context, set1, set2, rep);
}

void CGORenderAlpha(CGO* I, RenderInfo* info, bool calcDepth)
{
  CGORenderGLAlpha(I, info, calcDepth);
}
