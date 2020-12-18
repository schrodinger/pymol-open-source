
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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Feedback.h"
#include"CGO.h"
#include"Base.h"
#include"OOMac.h"
#include"Setting.h"
#include"Sphere.h"
#include"PConv.h"
#include"GadgetSet.h"
#include"VFont.h"
#include"P.h"
#include"PyMOLGlobals.h"
#include"Ray.h"
#include"Util.h"
#include"Scene.h"
#include"ScenePicking.h"
#include"Matrix.h"
#include"ShaderMgr.h"
#include"CoordSet.h"
#include"Rep.h"
#include"Vector.h"
#include"ObjectGadgetRamp.h"
#include"Triangle.h"
#include "Picking.h"

#include "pymol/algorithm.h"

#if defined(_PYMOL_IOS) && !defined(_WEBGL)
#define ALIGN_VBOS_TO_4_BYTE_ARRAYS
#define VAR_FOR_NORMAL  plc
#define VAR_FOR_NORMAL_CNT_PLUS   + (cnt / 3)
#define VERTEX_NORMAL_SIZE 4
#else
#define VAR_FOR_NORMAL  pl
#define VERTEX_NORMAL_SIZE 3
#define VAR_FOR_NORMAL_CNT_PLUS   
#endif

#define VALUES_PER_IMPOSTER_SPACE_COORD 1

#if defined(PURE_OPENGL_ES_2)
#define VERTICES_PER_SPHERE 6
#else
#define VERTICES_PER_SPHERE 4
#endif

#if defined(_PYMOL_IOS) && !defined(_WEBGL)
#define NUM_VERTICES_PER_CYLINDER 4
#define NUM_TOTAL_VERTICES_PER_CYLINDER 6
#else
#define NUM_VERTICES_PER_CYLINDER 8
#define NUM_TOTAL_VERTICES_PER_CYLINDER 36
#endif

constexpr unsigned VERTEX_PICKCOLOR_RGBA_SIZE = 1;  // 4 unsigned bytes
constexpr unsigned VERTEX_PICKCOLOR_INDEX_SIZE = 2; // index + bond
constexpr unsigned VERTEX_PICKCOLOR_SIZE = VERTEX_PICKCOLOR_RGBA_SIZE + //
                                           VERTEX_PICKCOLOR_INDEX_SIZE;
constexpr unsigned VERTEX_ACCESSIBILITY_SIZE = 1;

#ifdef PURE_OPENGL_ES_2
#define glVertexAttrib4ubv(loc, data) glVertexAttrib4f(loc, \
    (data)[0] / 255.f, (data)[1] / 255.f, (data)[2] / 255.f, (data)[3] / 255.f);
#endif

#include <cassert>
#include <iostream>
#include <algorithm>

#define MAX_INDICES_FOR_IOS 65536

const float g_ones4f[4] = {1.f, 1.f, 1.f, 1.f};

template <typename T>
inline T CLAMPVALUE(T val, T minimum, T maximum) {
  return 
    (val < minimum) ? minimum :
    (val > maximum) ? maximum : val;
}

inline cCylCap cap1_from_cyl_shader_bits(const unsigned char bits)
{
  return (bits & 1)
             ? ((bits & cCylShaderCap1RoundBit) ? cCylCap::Round : cCylCap::Flat)
             : cCylCap::None;
}
inline cCylCap cap2_from_cyl_shader_bits(const unsigned char bits)
{
  return (bits & 2)
             ? ((bits & cCylShaderCap2RoundBit) ? cCylCap::Round : cCylCap::Flat)
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

#if defined(_PYMOL_IOS) && !defined(_WEBGL) 
extern "C" void firePyMOLLimitationWarning();
#define CHECK_GL_ERROR_OK(printstr)			\
  if ((err = glGetError())!=0 || I->G->Interrupt != 0){		\
      if (err)        \
	PRINTFB(I->G, FB_CGO, FB_Errors) printstr, err ENDFB(I->G);	   \
  }
#else
#define CHECK_GL_ERROR_OK(printstr)	\
  if ((err = glGetError()) != 0){						\
     PRINTFB(I->G, FB_CGO, FB_Errors) printstr, err ENDFB(I->G);	   \
  }
#endif

#define WARN_UNEXPECTED_OPERATION(G, op)                                       \
  PRINTFB(G, FB_CGO, FB_Warnings)                                              \
  " %s-Warning: unexpected op=0x%x (line %d)\n", __func__, op, __LINE__ ENDFB(G)

// like g_return_val_if_fail from glib
#define RETURN_VAL_IF_FAIL(expr, val)                                          \
  {                                                                            \
    if (!(expr))                                                               \
      return (val);                                                            \
  }

struct _CCGORenderer {
  PyMOLGlobals *G;
  RenderInfo *info;
  Rep *rep;
  const float *color;
  float alpha;
  short sphere_quality;
  bool isPicking;
  unsigned pick_pass() const { return info->pick->m_pass; }
  bool use_shader; // OpenGL 1.4+, e.g., glEnableVertexAttribArray() (on) vs. glEnableClientState() (off)
  bool debug;
  CSetting *set1, *set2;
};

static
void set_current_pick_color(
    CGO * cgo,
    unsigned int idx,
    int bnd)
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

static
int CGOConvertDebugMode(int debug, int modeArg){
  int mode = modeArg;
  if (debug==1){
    switch (mode){
    case GL_TRIANGLES:
      mode = GL_LINES;
      break;
    case GL_TRIANGLE_STRIP:
      mode = GL_LINE_STRIP;
      break;
    case GL_TRIANGLE_FAN:
      mode = GL_LINES;
      break;
    }
  } else {
    mode = GL_POINTS;
  }
  return mode;
}


int CGORendererInit(PyMOLGlobals * G)
{
  CCGORenderer *I = NULL;

  I = (G->CGORenderer = pymol::calloc<CCGORenderer>(1));
  if(I) {
    I->G = G;
    I->isPicking = false;
    I->alpha = 1.0F;
    return 1;
  } else
    return 0;
}

void CGORendererFree(PyMOLGlobals * G)
{
  FreeP(G->CGORenderer);
}

int CGO_sz[] = {
  CGO_NULL_SZ,
  CGO_NULL_SZ,
  CGO_BEGIN_SZ,
  CGO_END_SZ,

  CGO_VERTEX_SZ,
  CGO_NORMAL_SZ,
  CGO_COLOR_SZ,
  CGO_SPHERE_SZ,

  CGO_TRIANGLE_SZ,
  fsizeof<cgo::draw::cylinder>(),
  CGO_LINEWIDTH_SZ,
  CGO_WIDTHSCALE_SZ,

  CGO_ENABLE_SZ,
  CGO_DISABLE_SZ,
  fsizeof<cgo::draw::sausage>(),
  fsizeof<cgo::draw::custom_cylinder>(),

  CGO_DOTWIDTH_SZ,
  CGO_ALPHA_TRIANGLE_SZ,
  CGO_ELLIPSOID_SZ,
  CGO_FONT_SZ,

  CGO_FONT_SCALE_SZ,
  CGO_FONT_VERTEX_SZ,
  CGO_FONT_AXES_SZ,
  CGO_CHAR_SZ,

  CGO_INDENT_SZ,
  CGO_ALPHA_SZ,
  CGO_QUADRIC_SZ,
  CGO_CONE_SZ,

  fsizeof<cgo::draw::arrays>(),
  CGO_NULL_SZ,
  CGO_RESET_NORMAL_SZ,
  CGO_PICK_COLOR_SZ,

  CGO_NULL_SZ, // CGO_DRAW_BUFFERS_SZ no longer used
  fsizeof<cgo::draw::buffers_indexed>(),
  CGO_BOUNDING_BOX_SZ,
  fsizeof<cgo::draw::buffers_not_indexed>(),
  CGO_SPECIAL_SZ,
  fsizeof<cgo::draw::cylinder_buffers>(),
  fsizeof<cgo::draw::shadercylinder>(),
  fsizeof<cgo::draw::shadercylinder2ndcolor>(),
  fsizeof<cgo::draw::sphere_buffers>(),
  CGO_ACCESSIBILITY_SZ,
  CGO_DRAW_TEXTURE_SZ,
  fsizeof<cgo::draw::textures>(),
  fsizeof<cgo::draw::screen_textures>(),
  CGO_TEX_COORD_SZ,
  fsizeof<cgo::draw::label>(),
  fsizeof<cgo::draw::labels>(),
  CGO_DRAW_CONNECTOR_SZ,
  fsizeof<cgo::draw::connectors>(),
  CGO_DRAW_TRILINES_SZ,  CGO_UNIFORM3F_SZ,
  CGO_SPECIAL_WITH_ARG_SZ,
  fsizeof<cgo::draw::line>(),
  fsizeof<cgo::draw::splitline>(),
  fsizeof<cgo::draw::custom>(),
  fsizeof<cgo::draw::vertex_attribute_3f>(),
  fsizeof<cgo::draw::vertex_attribute_4ub>(),
  fsizeof<cgo::draw::vertex_attribute_1f>(),
  fsizeof<cgo::draw::mask_attribute_if_picking>(),
  fsizeof<cgo::draw::bind_vbo_for_picking>(),
  CGO_VERTEX_BEGIN_LINE_STRIP_SZ,  CGO_INTERPOLATED_SZ,  CGO_VERTEX_CROSS_SZ,
  fsizeof<cgo::draw::vertex_attribute_4ub_if_picking>(),
  fsizeof<cgo::draw::custom_cylinder_alpha>(),
  CGO_NULL_SZ
};

/**
 * Get the number of elements in `CGO_sz`
 */
size_t CGO_sz_size()
{
  return sizeof(CGO_sz) / sizeof(*CGO_sz);
}

// I think CGO rendering functions should not modify CGO's, so the
// data pointer should be const. Current exception: `pickcolorsset`
#define CGO_OP_DATA_CONST const
typedef CGO_OP_DATA_CONST float* const* CGO_op_data;
typedef void CGO_op(CCGORenderer * I, CGO_op_data);
typedef CGO_op *CGO_op_fn;

static float *CGO_add(CGO * I, unsigned c);
static float *CGO_size(CGO * I, int sz);
static int CGOSimpleCylinder(CGO * I, const float *v1, const float *v2, const float tube_size, const float *c1,
                             const float *c2, float a1, const float a2, const bool interp,
                             const cCylCap cap1,
                             const cCylCap cap2,
                             const Pickable *pickcolor2 = nullptr, const bool stick_round_nub = false);
template<typename CylinderT>
static int CGOSimpleCylinder(CGO * I, const CylinderT &cyl, const float a1, const float a2, const bool interp,
                             const cCylCap cap1,
                             const cCylCap cap2,
                             const Pickable *pickcolor2 = nullptr, const bool stick_round_nub = false);
static int CGOSimpleEllipsoid(CGO * I, const float *v, float vdw, const float *n0, const float *n1,
			      const float *n2);
static int CGOSimpleQuadric(CGO * I, const float *v, float vdw, const float *q);
static int CGOSimpleSphere(CGO * I, const float *v, float vdw, short sphere_quality);
static int CGOSimpleCone(CGO * I, const float *v1, const float *v2, float r1, float r2, const float *c1,
			 const float *c2, cCylCap cap1, cCylCap cap2);


/**
 * Inverse function of CGOArrayFromPyListInPlace
 *
 * I: (input) Primitive CGO (may contain CGO_DRAW_ARRAYS)
 *
 * Return: All-float Python list primitive CGO
 */
static PyObject *CGOArrayAsPyList(const CGO * I)
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
    case CGO_DRAW_ARRAYS:
      {
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
    for(; sz; --sz) {
      flat.push_back(*(pc++));
    }
  }

  return PConvToPyObject(flat);
}

PyObject *CGOAsPyList(CGO * I)
{
  PyObject *result;
  result = PyList_New(2);
  PyObject *list = CGOArrayAsPyList(I);
  PyList_SetItem(result, 0, PyInt_FromLong(PyList_Size(list)));
  PyList_SetItem(result, 1, list);
  return (result);
}

static float CPythonVal_PyFloat_AsDouble_From_List(void * G, PyObject * list, size_t i) {
  float out;
  PConvPyFloatToFloat(PyList_GetItem(list, i), &out);
  return out;
}

/**
 * Inverse function of CGOArrayAsPyList
 *
 * list: (input) All-float Python list primitive CGO (may contain CGO_DRAW_ARRAYS)
 * I: (output) empty CGO
 */
static int CGOArrayFromPyListInPlace(PyObject * list, CGO * I)
{
  // sanity check
  if (!list || !PyList_Check(list))
    return false;

  auto G = I->G;

#define GET_FLOAT(i) ((float) CPythonVal_PyFloat_AsDouble_From_List(I->G, list, i))
#define GET_INT(i)   ((int)   CPythonVal_PyFloat_AsDouble_From_List(I->G, list, i))

  for (int i = 0, l = PyList_Size(list); i < l;) {
    unsigned op = GET_INT(i++);
    ok_assert(1, op < CGO_sz_size());
    int sz = CGO_sz[op];
    float * fdata = I->add_to_buffer(sz + 1);
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
    case CGO_DRAW_ARRAYS:
      {
        // has abstract superclass, need to be constructed!
        ok_assert(1, i + 3 < l);
        auto sp = new (fdata) cgo::draw::arrays(
            GET_INT(i),
            GET_INT(i + 1),
            GET_INT(i + 3));

        // sanity check
        int narrays_check = GET_INT(i + 2);
        if (sp->narrays != narrays_check) {
          PRINTFB(I->G, FB_CGO, FB_Warnings)
            " CGO-Warning: narrays mismatch: %d != %d\n",
            sp->narrays, narrays_check ENDFB(I->G);
        }

        // data
        sz = sp->get_data_length();
        sp->floatdata = fdata = I->allocate_in_data_heap(sz);

        i += 4;
      }
      break;
    }

    // float members
    for(; sz; --sz) {
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

CGO *CGONewFromPyList(PyMOLGlobals * G, PyObject * list, int version, bool shouldCombine)
{
  int ok = true;
  auto I = CGONew(G);
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  /* TO ENABLE BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if((version > 0) && (version <= 86)) {
    if(ok)
      ok = PConvFromPyListItem(G, list, 0, I->c);
    if(ok)
      VLACheck(I->op, float, I->c);
    if(ok)
      ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 1), I->op, I->c);
  } else {
    if(ok)
      ok = CGOArrayFromPyListInPlace(PyList_GetItem(list, 1), I);
  }
  if(!ok) {
    CGOFree(I);
  }
  {
    CGO *cgo = NULL;
    if (shouldCombine && I && I->has_begin_end){
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

void CGOSetUseShader(CGO *I, int use_shader){
  I->use_shader = use_shader;
  if (use_shader){
    I->cgo_shader_ub_color = SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color);
    I->cgo_shader_ub_normal = SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal);
  } else {
    I->cgo_shader_ub_color = 0;
    I->cgo_shader_ub_normal = 0;
  }
}
void CGOReset(CGO * I)
{
  I->c = 0;
  I->z_flag = false;
  I->alpha = 1.f;
  I->has_begin_end = false;
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->normal[0] = 0.f; I->normal[1] = 0.f; I->normal[2] = 1.f;
  I->color[0] = 0.f; I->color[1] = 0.f; I->color[2] = 1.f;
  I->pickColor[0] = 0; I->pickColor[1] = 0; I->pickColor[2] = 0; I->pickColor[3] = 255;
  I->current_accessibility = 1.f;
}

void CGOFree(CGO * &I, bool withVBOs)
{
  if(I) {
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

static float *CGO_add(CGO * I, unsigned c)
{
  float *at;
  VLACheck(I->op, float, I->c + c);
  if (!I->op){
    return NULL;
  }
  at = I->op + I->c;
  I->c += c;
  return (at);
}

static float *CGO_size(CGO * I, int sz)
{
  float *at;
  VLASize(I->op, float, sz);
  if (!I->op){
      return NULL;
  }
  at = I->op + I->c;
  I->c = sz;
  return (at);
}


/*===== Object Creation Routines =======*/

int CGOFromFloatArray(CGO * I, const float *src, int len)
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
  while(len-- > 0) {
    cc++;
    const auto op = static_cast<unsigned>(*(src++));

    if (op >= CGO_sz_size()) {
      bad_entry = cc;
      break;
    }

    sz = CGO_sz[op];
    if(len < sz)
      break;                    /* discard short instruction */
    len -= sz;
    pc = save_pc;
    CGO_write_int(pc, op);
    ok = true;
    for(a = 0; a < sz; a++) {
      cc++;
      val = *(src++);
      if((FLT_MAX - val) > 0.0F) {      /* make sure we have a real float */
        *(pc++) = val;
      } else {
        *(pc++) = 0.0;
        ok = false;
      }
    }
    if(ok) {
      switch (op) {
      case CGO_END:
      case CGO_VERTEX:
      case CGO_BEGIN:
	I->has_begin_end = true;
      }
      switch (op) {             /* now convert any instructions with int arguments */
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
    } else {                    /* discard illegal instructions */
      if(all_ok)
        bad_entry = cc;
      all_ok = false;
    }
  }
  return (bad_entry);
}

int CGOBegin(CGO * I, int mode)
{
  float *pc = CGO_add(I, CGO_BEGIN_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_BEGIN);
  CGO_write_int(pc, mode);
  I->has_begin_end = true;
  I->texture[0] = 0.f;
  I->texture[1] = 0.f;
  return true;
}

int CGOEnd(CGO * I)
{
  float *pc = CGO_add(I, CGO_END_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_END);
  I->has_begin_end = true;
  return true;
}

int CGOEnable(CGO * I, int mode)
{
  float *pc = CGO_add(I, CGO_ENABLE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ENABLE);
  CGO_write_int(pc, mode);
  return true;
}

int CGODisable(CGO * I, int mode)
{
  float *pc = CGO_add(I, CGO_DISABLE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_DISABLE);
  CGO_write_int(pc, mode);
  return true;
}

int CGOLinewidth(CGO * I, float v)
{
  float *pc = CGO_add(I, CGO_LINEWIDTH_SZ + 1);
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
int CGOSpecial(CGO * I, int v)
{
  float *pc = CGO_add(I, CGO_SPECIAL_SZ + 1);
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
int CGOSpecialWithArg(CGO * I, int v, float argval)
{
  float *pc = CGO_add(I, CGO_SPECIAL_WITH_ARG_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_SPECIAL_WITH_ARG);
  CGO_write_int(pc, v);
  *pc = argval;
  return true;
}

int CGODotwidth(CGO * I, float v)
{
  float *pc = CGO_add(I, CGO_DOTWIDTH_SZ + 1);
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
int CGOUniform3f(CGO *I, int uniform_id, const float *value){
  float *pc = CGO_add(I, CGO_UNIFORM3F_SZ + 1);
  if (!pc)
    return 0;
  CGO_write_int(pc, CGO_UNIFORM3F);
  CGO_write_int(pc, uniform_id);
  copy3f(value, pc);
  return pc - I->op;
}

int CGOBoundingBox(CGO *I, const float *min, const float *max){
  float *pc = CGO_add(I, CGO_BOUNDING_BOX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_BOUNDING_BOX);
  *(pc++) = *(min);
  *(pc++) = *(min+1);
  *(pc++) = *(min+2);
  *(pc++) = *(max);
  *(pc++) = *(max+1);
  *(pc++) = *(max+2);
  return true;
}

int CGOAccessibility(CGO * I, float a)
{
  float *pc = CGO_add(I, CGO_ACCESSIBILITY_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ACCESSIBILITY);
  *(pc++) = a;
  return true;
}

int CGODrawTexture(CGO *I, int texture_id, float *worldPos, float *screenMin, float *screenMax, float *textExtent)
{
  float *pc = CGO_add(I, CGO_DRAW_TEXTURE_SZ + 1);
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

int CGODrawConnector(CGO *I, float *targetPt3d, float *labelCenterPt3d, float text_width, float text_height, float *indentFactor, float *screenWorldOffset, float *connectorColor, short relativeMode, int draw_flags, float bkgrd_transp, float *bkgrd_color, float rel_ext_length, float connectorWidth)
{
  float *pc = CGO_add(I, CGO_DRAW_CONNECTOR_SZ + 1);
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
  *(pc++) = rel_ext_length; /* place for ext_length relative to height (i.e., text_height which is total height */
  *(pc++) = screenWorldOffset[0];
  *(pc++) = screenWorldOffset[1];
  *(pc++) = screenWorldOffset[2];
  *(pc++) = text_width;
  *(pc++) = text_height;
  *(pc++) = connectorColor[0];
  *(pc++) = connectorColor[1];
  *(pc++) = connectorColor[2];
  *(pc++) = (float)relativeMode;
  *(pc++) = (float)draw_flags;
  *(pc++) = bkgrd_color[0];
  *(pc++) = bkgrd_color[1];
  *(pc++) = bkgrd_color[2];
  *(pc++) = bkgrd_transp;
  *(pc++) = connectorWidth; // place for label_connector_width
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGODrawLabel(CGO *I, int texture_id, float *targetPos, float *worldPos, float *screenWorldOffset, float *screenMin, float *screenMax, float *textExtent, short relativeMode)
{
  float *pc = CGO_add(I, CGO_DRAW_LABEL_SZ + 1);
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
  *(pc++) = (float)relativeMode;
  *(pc++) = targetPos[0];
  *(pc++) = targetPos[1];
  *(pc++) = targetPos[2];
  return true;
}
#endif

#ifdef WITH_UNUSED_FUNCTIONS
int CGOConev(CGO * I,
    const float *p1,
    const float *p2, float r1, float r2,
    const float *c1,
    const float *c2,
              float cap1, float cap2)
{
  float *pc = CGO_add(I, CGO_CONE_SZ + 1);
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

int CGOPickColor(CGO * I, unsigned int index, int bond)
{
  // check if uchar is -1 since extrude does this for masked atoms
  if (index == (unsigned int)-1){
    bond = cPickableNoPick;
  }
  if (I->current_pick_color_index==index &&
      I->current_pick_color_bond==bond)
    return true;
  float *pc = CGO_add(I, CGO_PICK_COLOR_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_PICK_COLOR);
  CGO_write_uint(pc, index);
  CGO_write_int(pc, bond);
  I->current_pick_color_index = index;
  I->current_pick_color_bond = bond;
  return true;
}

int CGOAlpha(CGO * I, float alpha)
{
  float *pc = CGO_add(I, CGO_ALPHA_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_ALPHA);
  *(pc++) = alpha;
  I->alpha = alpha;
  return true;
}

int CGOSphere(CGO * I, const float *v1, float r)
{
  float *pc = CGO_add(I, CGO_SPHERE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_SPHERE);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = r;
  return true;
}

int CGOEllipsoid(CGO * I, const float *v1, float r,
    const float *n1,
    const float *n2,
    const float *n3)
{
  float *pc = CGO_add(I, CGO_ELLIPSOID_SZ + 1);
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
int CGOQuadric(CGO * I, const float *v, float r, const float *q)
{
  float *pc = CGO_add(I, CGO_QUADRIC_SZ + 1);
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

void CGOSetZVector(CGO * I, float z0, float z1, float z2)
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

int CGOAlphaTriangle(CGO * I,
                      const float *v1, const float *v2, const float *v3,
                      const float *n1, const float *n2, const float *n3,
                      const float *c1, const float *c2, const float *c3,
                      float a1, float a2, float a3, int reverse)
{
  if(v1 && v2 && v3) {
    float *pc = CGO_add(I, CGO_ALPHA_TRIANGLE_SZ + 1);
    float z = _0;
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_ALPHA_TRIANGLE);
    CGO_write_int(pc, 0);  // this is the place for the next triangle in the bin
    *(pc++) = (v1[0] + v2[0] + v3[0]) * one_third;
    *(pc++) = (v1[1] + v2[1] + v3[1]) * one_third;
    *(pc++) = (v1[2] + v2[2] + v3[2]) * one_third;
    if(I->z_flag) {
      float *zv = I->z_vector;
      z = pc[-3] * zv[0] + pc[-2] * zv[1] + pc[-1] * zv[2];
      if(z > I->z_max)
        I->z_max = z;
      if(z < I->z_min)
        I->z_min = z;
    }
    *(pc++) = z;

    if(reverse) {
      *(pc++) = *(v2++);        /* vertices @ +5 */
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
    } else {
      *(pc++) = *(v1++);        /* vertices @ +5 */
      *(pc++) = *(v1++);
      *(pc++) = *(v1++);
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
      *(pc++) = *(v2++);
    }

    *(pc++) = *(v3++);
    *(pc++) = *(v3++);
    *(pc++) = *(v3++);

    if(reverse) {
      *(pc++) = *(n2++);        /* normals @ +14 */
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
    } else {
      *(pc++) = *(n1++);        /* normals @ +14 */
      *(pc++) = *(n1++);
      *(pc++) = *(n1++);
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
      *(pc++) = *(n2++);
    }
    *(pc++) = *(n3++);
    *(pc++) = *(n3++);
    *(pc++) = *(n3++);

    if(reverse) {
      *(pc++) = *(c2++);        /* colors @ +23 */
      *(pc++) = *(c2++);
      *(pc++) = *(c2++);
      *(pc++) = a2;
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = *(c1++);
      *(pc++) = a1;
    } else {
      *(pc++) = *(c1++);        /* colors @ +23 */
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

int CGOVertex(CGO * I, float v1, float v2, float v3)
{
  float *pc = CGO_add(I, CGO_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
  return true;
}

int CGOVertexv(CGO * I, const float *v)
{
  float *pc = CGO_add(I, CGO_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

int CGOVertexCrossv(CGO * I, const float *v)
{
  float *pc = CGO_add(I, CGO_VERTEX_CROSS_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_VERTEX_CROSS);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOInterpolated(CGO * I, const bool interp)
{
  float *pc = CGO_add(I, CGO_INTERPOLATED_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_INTERPOLATED);
  *(pc++) = interp ? 1.f : 0.f;
  I->interpolated = interp;
  return true;
}
#endif

int CGOColor(CGO * I, float v1, float v2, float v3)
{
  float *pc = CGO_add(I, CGO_COLOR_SZ + 1);
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

int CGOColorv(CGO * I, const float *v)
{
  return CGOColor(I, v[0], v[1], v[2]);
}

int CGOTexCoord2f(CGO * I, float v1, float v2){
  float *pc = CGO_add(I, CGO_TEX_COORD_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_TEX_COORD);
  *(pc++) = v1;
  *(pc++) = v2;
  I->texture[0] = v1;
  I->texture[1] = v2;
  return true;
}

int CGONormal(CGO * I, float v1, float v2, float v3)
{
  float *pc = CGO_add(I, CGO_NORMAL_SZ + 1);
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

int CGOResetNormal(CGO * I, int mode)
{
  float *pc = CGO_add(I, CGO_RESET_NORMAL_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_RESET_NORMAL);
  CGO_write_int(pc, mode);
  SceneGetResetNormal(I->G, I->normal, mode);
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOFontVertexv(CGO * I, const float *v)
{
  float *pc = CGO_add(I, CGO_FONT_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_FONT_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  return true;
}

int CGOFontVertex(CGO * I, float x, float y, float z)
{
  float *pc = CGO_add(I, CGO_FONT_VERTEX_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_FONT_VERTEX);
  *(pc++) = x;
  *(pc++) = y;
  *(pc++) = z;
  return true;
}
#endif

int CGOFontScale(CGO * I, float v1, float v2)
{
  float *pc = CGO_add(I, CGO_FONT_SCALE_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_FONT_SCALE);
  *(pc++) = v1;
  *(pc++) = v2;
  return true;
}

#ifdef WITH_UNUSED_FUNCTIONS
int CGOChar(CGO * I, char c)
{
  float *pc = CGO_add(I, CGO_CHAR_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_CHAR);
  *(pc++) = (float) c;
  return true;
}

int CGOIndent(CGO * I, char c, float dir)
{
  float *pc = CGO_add(I, CGO_INDENT_SZ + 1);
  if (!pc)
    return false;
  CGO_write_int(pc, CGO_INDENT);
  *(pc++) = (float) c;
  *(pc++) = dir;
  return true;
}

int CGOWrite(CGO * I, const char *str)
{
  float *pc;

  while(*str) {
    pc = CGO_add(I, CGO_CHAR_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(str++);
  }
  return true;
}

int CGOWriteLeft(CGO * I, const char *str)
{
  float *pc;
  const char *s = str;
  while(*s) {
    pc = CGO_add(I, CGO_INDENT_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_INDENT);
    *(pc++) = (float) *(s++);
    *(pc++) = -1.0F;
  }
  s = str;
  while(*s) {
    pc = CGO_add(I, CGO_CHAR_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(s++);
  }
  return true;
}

int CGOWriteIndent(CGO * I, const char *str, float indent)
{
  float *pc;
  const char *s = str;
  while(*s) {
    pc = CGO_add(I, CGO_INDENT_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_INDENT);
    *(pc++) = (float) *(s++);
    *(pc++) = indent;
  }
  s = str;
  while(*s) {
    pc = CGO_add(I, CGO_CHAR_SZ + 1);
    if (!pc)
      return false;
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(s++);
  }
  return true;
}
#endif

int CGONormalv(CGO * I, const float *v)
{
  float *pc = CGO_add(I, CGO_NORMAL_SZ + 1);
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
int CGOStop(CGO * I)
{
#define CGO_STOP_ZEROS 1

  float *pc = CGO_size(I, I->c + CGO_STOP_ZEROS);
  if (!pc)
    return false;
  UtilZeroMem(pc, sizeof(float) * CGO_STOP_ZEROS);
  I->c -= CGO_STOP_ZEROS;
  return true;
}

int CGOCheckComplex(CGO * I)
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
    case CGO_DRAW_ARRAYS:
      {
        cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
	fc += sp->nverts;
      }
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
      {
        cgo::draw::buffers_indexed * sp = reinterpret_cast<decltype(sp)>(pc);
	switch(sp->mode){
	case GL_TRIANGLES:
	  fc += sp->nindices / 3;
	  break;
	case GL_LINES:
	  fc += sp->nindices / 2;
	  break;
	}
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
        cgo::draw::buffers_not_indexed * sp = reinterpret_cast<decltype(sp)>(pc);
	switch(sp->mode){
	case GL_TRIANGLES:
	  fc += sp->nverts / 3;
	  break;
	case GL_LINES:
	  fc += sp->nverts / 2;
	  break;
	}
      }
      break;
    case CGO_DRAW_SPHERE_BUFFERS:
      {
        cgo::draw::sphere_buffers * sp = reinterpret_cast<decltype(sp)>(pc);
        fc += sp->num_spheres * VERTICES_PER_SPHERE;
      }
      break;
    case CGO_DRAW_CYLINDER_BUFFERS:
      {
        cgo::draw::cylinder_buffers * sp = reinterpret_cast<decltype(sp)>(pc);
        fc += sp->num_cyl * NUM_VERTICES_PER_CYLINDER;
      }
      break;
    }
  }
  return (fc);
}

int CGOPreloadFonts(CGO * I)
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
      if(!font_seen) {
        font_id = VFontLoad(I->G, 1.0, 1, 1, true);
        ok = ok && font_id;
        font_seen = true;
      }
      break;
    }
  }
  if(blocked)
    PUnblock(I->G);
  return (ok);
}

int CGOCheckForText(CGO * I)
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
      fc += 3 + 2 * 3 * 10;     /* est 10 lines per char */
      break;
    }
  }
  PRINTFD(I->G, FB_CGO)
    " CGOCheckForText-Debug: %d\n", fc ENDFD;

  return (fc);
}

CGO *CGODrawText(const CGO * I, int est, float *camera)
{                               /* assumes blocked intepreter */
  CGO *cgo;
  int font_id = 0;
  char text[2] = " ";
  float pos[] = { 0.0F, 0.0F, 0.0F };
  float axes[] = { 1.0F, 0.0F, 0.0F,
    0.0F, 1.0F, 0.0F,
    0.0F, 0.0F, 1.0F
  };
  float scale[2] = { 1.0, 1.0 };

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
      if(!font_id) {
        font_id = VFontLoad(I->G, 1.0, 1, 1, false);
      }
      text[0] = (unsigned char) *pc;
      VFontWriteToCGO(I->G, font_id, cgo, text, pos, scale, axes, cgo->color);
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  if (cgo && cgo->has_begin_end){
    /* this is mainly for VFontWriteToCGO() that still creates CGOBegin/CGOEnd */
    if(cgo && cgo->has_begin_end){
      CGO *convertcgo = NULL;
      convertcgo = CGOCombineBeginEnd(cgo, 0);
      CGOFree(cgo);
      cgo = convertcgo;
    }
  }
  return (cgo);
}

static
void CGOAddVertexToDrawArrays(CGO *cgo, int pl, int plc, int pla, const float *vertex,
                              short notHaveValue, float *vertexVals, float *normalVals,
                              float *colorVals, float *pickColorVals, float *accessibilityVals){
  float *tmp_ptr;
  if (notHaveValue & CGO_NORMAL_ARRAY){
    if (pl){
      tmp_ptr = &normalVals[pl-3];
      copy3f(tmp_ptr, &normalVals[pl]);
    } else {
      copy3f(cgo->normal, &normalVals[pl]);
    }
  }
  if (notHaveValue & CGO_COLOR_ARRAY){
    if (plc){
      tmp_ptr = &colorVals[plc-4];
      copy4f(tmp_ptr, &colorVals[plc]);
    } else {
      copy3f(&colorVals[plc], cgo->color);
      colorVals[plc+3] = cgo->alpha;
    }
  }
  if (pickColorVals){
    CGO_put_uint(pickColorVals + pla * 2, cgo->current_pick_color_index);
    CGO_put_int(pickColorVals + pla * 2 + 1, cgo->current_pick_color_bond);
  }
  if (accessibilityVals){
    accessibilityVals[pla] = cgo->current_accessibility;
  }
  copy3f(vertex, &vertexVals[pl]);
}

bool CGOCombineBeginEnd(CGO ** I, bool do_not_split_lines) {
  CGO *cgo = CGOCombineBeginEnd(*I, 0, do_not_split_lines);
  CGOFree(*I);
  *I = cgo;
  return (cgo != NULL);
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
CGO *CGOCombineBeginEnd(const CGO * I, int est, bool do_not_split_lines)
{
  CGO *cgo;

  int ok = true;
  if (!I)
      return NULL;
  cgo = CGONewSized(I->G, 0);
  ok &= cgo ? true : false;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_END:
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings)
        " CGOCombineBeginEnd: op=0x%02x encountered without CGO_BEGIN\n", op
        ENDFB(I->G);
      break;
    case CGO_BEGIN:
      {
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
	      " CGO-Error: CGOCombineBeginEnd: invalid op=0x%02x inside BEGIN/END\n",
	      it.op_code() ENDFB(I->G);
	    err = true;
	    continue;
	  case CGO_NORMAL:
	    damode |= CGO_NORMAL_ARRAY;
	    break;
	  case CGO_COLOR:
	    if (!nverts){
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
	    if (!nverts){
	      hasFirstAlpha = 1;
	      firstAlpha = cgo->alpha;
	    } else {
	      hasFirstAlpha = 0;
	      damode |= CGO_COLOR_ARRAY;
	    }
            break;
          case CGO_LINE:
            nverts+=2;
            break;
          case CGO_SPLITLINE:
            {
              auto splitline = reinterpret_cast<const cgo::draw::splitline *>(pc);
              if (do_not_split_lines || (splitline->flags & cgo::draw::splitline::equal_colors)){
                nverts+=2;
              } else {
                nverts+=4;
              }
            }
            break;
	  }
	}
	if (nverts>0 && !err){
	  int pl = 0, plc = 0, pla = 0;
          float *vertexVals = nullptr, *normalVals = nullptr,
                *colorVals = nullptr, *pickColorVals = nullptr,
                *accessibilityVals = nullptr;

	  if (hasFirstAlpha || hasFirstColor){
	    if (hasFirstAlpha){
	      CGOAlpha(cgo, firstAlpha);
	    }
	    if (hasFirstColor){
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

          if (damode & CGO_NORMAL_ARRAY){
            normalVals = nxtVals;
            nxtVals += nverts * VERTEX_NORMAL_SIZE;
	  }

	  if (damode & CGO_COLOR_ARRAY){
            colorVals = nxtVals;
            nxtVals += nverts * VERTEX_COLOR_SIZE;
	  }

	  if (damode & CGO_PICK_COLOR_ARRAY){
            pickColorVals = nxtVals + VERTEX_PICKCOLOR_RGBA_SIZE * nverts;
            nxtVals += nverts * VERTEX_PICKCOLOR_SIZE;
	  }

	  if (damode & CGO_ACCESSIBILITY_ARRAY){
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
	      if (colorVals){
		copy3f(pc, &colorVals[plc]);
		colorVals[plc+3] = cgo->alpha;
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
	    case CGO_SPLITLINE:
              {
                auto splitline = reinterpret_cast<const cgo::draw::splitline *>(pc);
                float color2[] = { CONVERT_COLOR_VALUE(splitline->color2[0]),
                                   CONVERT_COLOR_VALUE(splitline->color2[1]),
                                   CONVERT_COLOR_VALUE(splitline->color2[2]) };
                if (do_not_split_lines || (splitline->flags & cgo::draw::splitline::equal_colors)){
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex1, notHaveValue, vertexVals,
                                           normalVals, colorVals, pickColorVals, accessibilityVals);
                  pl+=3; plc+=4; pla++;
                  notHaveValue = damode;

                  if (!(splitline->flags & cgo::draw::splitline::equal_colors)){
                    if (colorVals){
                      copy3f(color2, &colorVals[plc]);
                      colorVals[plc+3] = cgo->alpha;
                      notHaveValue = notHaveValue & ~CGO_COLOR_ARRAY;
                    }
                    copy3f(color2, cgo->color);
                  }
                  if (pickColorVals){
                    cgo->current_pick_color_index = splitline->index;
                    cgo->current_pick_color_bond = splitline->bond;
                    notHaveValue = notHaveValue & ~CGO_PICK_COLOR_ARRAY;
                  }
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex2, notHaveValue, vertexVals,
                                           normalVals, colorVals, pickColorVals, accessibilityVals);
                  pl+=3; plc+=4; pla++;
                  notHaveValue = damode;
                } else {
                  float mid[3];
                  add3f(splitline->vertex1, splitline->vertex2, mid);
                  mult3f(mid, .5f, mid);
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex1, notHaveValue, vertexVals,
                                           normalVals, colorVals, pickColorVals, accessibilityVals);
                  notHaveValue = damode;
                  pl+=3; plc+=4; pla++;
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, mid, notHaveValue, vertexVals,
                                           normalVals, colorVals, pickColorVals, accessibilityVals);
                  pl+=3; plc+=4; pla++;
                  if (colorVals){
                    copy3f(color2, &colorVals[plc]);
                    colorVals[plc+3] = cgo->alpha;
                    notHaveValue = notHaveValue & ~CGO_COLOR_ARRAY;
                  }
                  copy3f(color2, cgo->color);
                  if (pickColorVals){
                    cgo->current_pick_color_index = splitline->index;
                    cgo->current_pick_color_bond = splitline->bond;
                    notHaveValue = notHaveValue & ~CGO_PICK_COLOR_ARRAY;
                  }
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, mid, notHaveValue, vertexVals,
                                           normalVals, colorVals, pickColorVals, accessibilityVals);
                  notHaveValue = damode;
                  pl+=3; plc+=4; pla++;
                  CGOAddVertexToDrawArrays(cgo, pl, plc, pla, splitline->vertex2, notHaveValue, vertexVals,
                                           normalVals, colorVals, pickColorVals, accessibilityVals);
                  pl+=3; plc+=4; pla++;
                  notHaveValue = damode;
                }
              }
              break;
	    case CGO_LINE:
              {
                auto line = reinterpret_cast<const cgo::draw::line *>(pc);
                CGOAddVertexToDrawArrays(cgo, pl, plc, pla, line->vertex1, notHaveValue, vertexVals,
                                         normalVals, colorVals, pickColorVals, accessibilityVals);
                pl+=3; plc+=4; pla++;
                notHaveValue = damode;
                CGOAddVertexToDrawArrays(cgo, pl, plc, pla, line->vertex2, notHaveValue, vertexVals,
                                         normalVals, colorVals, pickColorVals, accessibilityVals);
                pl+=3; plc+=4; pla++;
              }
	      break;
	    case CGO_VERTEX:
              CGOAddVertexToDrawArrays(cgo, pl, plc, pla, pc, notHaveValue, vertexVals,
                                       normalVals, colorVals, pickColorVals, accessibilityVals);
	      pl+=3; plc+=4; pla++;
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
      }
      break;
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
  if (ok){
    ok &= CGOStop(cgo);
    if (ok){
      cgo->use_shader = I->use_shader;
      if (cgo->use_shader){
	cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
	cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
      }
    }
  }
  if (!ok){
    CGOFree(cgo);
  }
  return (cgo);
}

/**
 * Release all shader resources from this CGO
 */
void CGOFreeVBOs(CGO * I) {
  constexpr bool freevbos = true;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    switch (op) {
    case CGO_DRAW_TRILINES:
    {
      unsigned buf = it.cast<cgo::draw::trilines>()->buffer;
      if (freevbos)
        I->G->ShaderMgr->AddVBOToFree(buf);
    }
    break;
    case CGO_DRAW_CUSTOM:
    {
      auto sp = it.cast<cgo::draw::custom>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->iboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    }
    break;
    case CGO_DRAW_SPHERE_BUFFERS:
    {
      auto sp = it.cast<cgo::draw::sphere_buffers>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    }
    break;
    case CGO_DRAW_LABELS:
    {
      auto sp = it.cast<cgo::draw::labels>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    }
    break;
    case CGO_DRAW_TEXTURES:
    {
      auto sp = it.cast<cgo::draw::textures>();
      if (freevbos)
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
    }
    break;
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
    {
      auto sp = it.cast<cgo::draw::screen_textures>();
      if (freevbos)
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
    }
    break;
    case CGO_DRAW_CYLINDER_BUFFERS:
    {
      auto sp = it.cast<cgo::draw::cylinder_buffers>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->iboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    }
    break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
    {
      auto sp = it.cast<cgo::draw::buffers_not_indexed>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
        I->G->ShaderMgr->freeGPUBuffer(sp->pickvboid);
      }
    }
    break;
    case CGO_DRAW_BUFFERS_INDEXED:
    {
      auto sp = it.cast<cgo::draw::buffers_indexed>();
      if (freevbos) {
        I->G->ShaderMgr->freeGPUBuffers({ sp->vboid, sp->iboid, sp->pickvboid });
      }
    }
    break;
    case CGO_DRAW_CONNECTORS:
    {
      auto sp = it.cast<cgo::draw::connectors>();
      if (freevbos)
        I->G->ShaderMgr->freeGPUBuffer(sp->vboid);
    }
    break;
    }
  }
}

#define set_min_max(mn, mx, pt) {		\
  if (mn[0]>*pt) mn[0] = *pt; \
  if (mn[1]>*(pt+1)) mn[1] = *(pt+1);		\
  if (mn[2]>*(pt+2)) mn[2] = *(pt+2);		\
  if (mx[0]<*pt) mx[0] = *pt; \
  if (mx[1]<*(pt+1)) mx[1] = *(pt+1);		\
  if (mx[2]<*(pt+2)) mx[2] = *(pt+2);}


static void CGOCountNumVertices(const CGO *I, int *num_total_vertices, int *num_total_indexes,
			 int *num_total_vertices_lines, int *num_total_indexes_lines,
			 int *num_total_vertices_points);

#ifdef WITH_UNUSED_FUNCTIONS
void CGOCountNumVerticesDEBUG(const CGO *I){
  int num_total_vertices=0, num_total_indexes=0, num_total_vertices_lines=0, num_total_indexes_lines=0, num_total_vertices_points=0;
  CGOCountNumVertices(I, &num_total_vertices, &num_total_indexes, &num_total_vertices_lines, &num_total_indexes_lines, &num_total_vertices_points);
  printf("CGOCountNumVerticesDEBUG: num_total_vertices=%d num_total_indexes=%d num_total_vertices_lines=%d num_total_indexes_lines=%d num_total_vertices_points=%d\n", num_total_vertices, num_total_indexes, num_total_vertices_lines, num_total_indexes_lines, num_total_vertices_points);
  
}
#endif

static void CGOCountNumVertices(const CGO *I, int *num_total_vertices, int *num_total_indexes,
			 int *num_total_vertices_lines, int *num_total_indexes_lines,
			 int *num_total_vertices_points){
  int verts_skipped = 0;
  short err = 0;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();

    err = 0;
    switch (op) {
    case CGO_DRAW_ARRAYS:
      {
        const auto sp = it.cast<cgo::draw::arrays>();
	short shouldCompress = false, shouldCompressLines = false, shouldCompressPoints = false;
	switch(sp->mode){
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
	if (!shouldCompress && !shouldCompressLines && !shouldCompressPoints){
	  verts_skipped += sp->nverts;
	} else if (shouldCompressLines) {
	  *num_total_vertices_lines += sp->nverts;
	  switch(sp->mode){
	  case GL_LINE_LOOP:
	    *num_total_indexes_lines += 2 * sp->nverts;
	    break;
	  case GL_LINE_STRIP:
	    *num_total_indexes_lines += 2 * (sp->nverts - 1);
	    break;
	  case GL_LINES:
	    *num_total_indexes_lines += sp->nverts;
	    break;
	  }
	} else if (shouldCompress){
	  *num_total_vertices += sp->nverts;
	  switch(sp->mode){
	  case GL_TRIANGLE_FAN:
	    *num_total_indexes += 3 * (sp->nverts - 2);
	    break;
	  case GL_TRIANGLE_STRIP:
	    *num_total_indexes += 3 * (sp->nverts - 2);
	    break;
	  case GL_TRIANGLES:
	    *num_total_indexes += sp->nverts;
	    break;
	  }
	} else if (shouldCompressPoints){
	  *num_total_vertices_points += sp->nverts;
	}
      }
	break;
    case CGO_END:
      if (!err){
	PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOCountNumVertices: CGO_END encountered, should call CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
	err = true;
      }
    case CGO_VERTEX:
      if (!err){
	PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOCountNumVertices: CGO_VERTEX encountered, should call CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
	err = true;
      }
    case CGO_BEGIN:
      if (!err){
	PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOCountNumVertices: CGO_BEGIN encountered, should call CGOCombineBeginEnd before CGOCountNumVertices\n" ENDFB(I->G);
	err = true;
      }
    default:
      break;
    }
  }
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

static
void SetVertexValuesForVBO(PyMOLGlobals * G, CGO *cgo, int pl, int plc, int cnt, int incr,
                           const float *vertexValsDA,
                           const float *normalValsDA,
                           const float *colorValsDA,
                           const float *pickColorValsDA,
			   float *vertexVals, uchar *normalValsC,
                           float *normalVals, uchar *colorValsUC, float *colorVals,
                           float *pickColorVals, 
                           float *accessibilityVals=NULL,
                           const float *accessibilityValsDA=nullptr)
{
  int pl2 = pl + 1, pl3 = pl + 2;
  int pln1 = VAR_FOR_NORMAL, pln2 = VAR_FOR_NORMAL + 1, pln3 = VAR_FOR_NORMAL + 2;
  int plc2 = plc + 1, plc3 = plc + 2, plc4 = plc + 3;
  int c, c2, c3;
  int cc, cc2, cc3, cc4;
  int pcc = incr * 2, pcco = cnt * 2;
  c = cnt * 3; c2 = c + 1; c3 = c + 2;
  cc = cnt * 4; cc2 = cc + 1; cc3 = cc + 2; cc4 = cc + 3;
  vertexVals[pl] = vertexValsDA[c]; vertexVals[pl2] = vertexValsDA[c2]; vertexVals[pl3] = vertexValsDA[c3];

    if (normalValsC){
      if (normalValsDA){
	normalValsC[pln1] = CLIP_NORMAL_VALUE(normalValsDA[c]); normalValsC[pln2] = CLIP_NORMAL_VALUE(normalValsDA[c2]); normalValsC[pln3] = CLIP_NORMAL_VALUE(normalValsDA[c3]);
      } else {
	normalValsC[pln1] = CLIP_NORMAL_VALUE(cgo->normal[0]); normalValsC[pln2] = CLIP_NORMAL_VALUE(cgo->normal[1]); normalValsC[pln3] = CLIP_NORMAL_VALUE(cgo->normal[2]);
      }
    
#ifdef ALIGN_VBOS_TO_4_BYTE_ARRAYS
      normalValsC[pln3+1] = 127;
#endif
  } else {
    if (normalValsDA){
	normalVals[pln1] = normalValsDA[c]; normalVals[pln2] = normalValsDA[c2]; normalVals[pln3] = normalValsDA[c3];
      } else {
	normalVals[pln1] = cgo->normal[0]; normalVals[pln2] = cgo->normal[1]; normalVals[pln3] = cgo->normal[2];
      }
  }

  if (colorValsUC){
    if (colorValsDA){
      colorValsUC[plc] = CLIP_COLOR_VALUE(colorValsDA[cc]); colorValsUC[plc2] = CLIP_COLOR_VALUE(colorValsDA[cc2]); 
      colorValsUC[plc3] = CLIP_COLOR_VALUE(colorValsDA[cc3]); colorValsUC[plc4] = CLIP_COLOR_VALUE(colorValsDA[cc4]);
    } else {
      colorValsUC[plc] = CLIP_COLOR_VALUE(cgo->color[0]); colorValsUC[plc2] = CLIP_COLOR_VALUE(cgo->color[1]);
      colorValsUC[plc3] = CLIP_COLOR_VALUE(cgo->color[2]); colorValsUC[plc4] = CLIP_COLOR_VALUE(cgo->alpha);
    }
  } else {
    if (colorValsDA){
      colorVals[plc] = colorValsDA[cc]; colorVals[plc2] = colorValsDA[cc2];
      colorVals[plc3] = colorValsDA[cc3]; colorVals[plc4] = colorValsDA[cc4];
    } else {
      colorVals[plc] = cgo->color[0]; colorVals[plc2] = cgo->color[1];
      colorVals[plc3] = cgo->color[2]; colorVals[plc4] = cgo->alpha;
    }
  }
  if (pickColorValsDA){
    cgo->current_pick_color_index = CGO_get_uint(pickColorValsDA + pcco);
    cgo->current_pick_color_bond = CGO_get_int(pickColorValsDA + pcco + 1);
  }
  CGO_put_uint(pickColorVals + pcc, cgo->current_pick_color_index);
  CGO_put_int(pickColorVals + pcc + 1, cgo->current_pick_color_bond);
  if (accessibilityValsDA){
    accessibilityVals[pl/3] = accessibilityValsDA[cnt];
  }
}

static int OptimizePointsToVBO(const CGO *I, CGO *cgo, int num_total_vertices_points, float *min, float *max, short *has_draw_buffer, bool addshaders){
  auto G = I->G;
  float *vertexVals = 0, *colorVals = 0, *normalVals = 0;
  float *pickColorVals;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
  uchar *colorValsUC = 0;
  uchar *normalValsC = 0;
  bool has_normals = false, has_colors = false;
  int ok = true;
  
  cgo->alpha = 1.f;
  cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;

  unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal) ? 1 : VERTEX_NORMAL_SIZE;
  mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
  auto const tot = size_t(num_total_vertices_points) * mul;

  vertexVals = pymol::malloc<float>(tot);
  CHECKOK(ok, vertexVals);
  if (!ok){
    PRINTFB(G, FB_CGO, FB_Errors)
    "%s-Error(%d): vertexVals could not be allocated (tot=%zu)\n", __func__,
        __LINE__, tot ENDFB(G);
    return 0;
  }
  normalVals = vertexVals + 3 * num_total_vertices_points;
  nxtn = 3;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
    normalValsC = (uchar*) normalVals;
    nxtn = 1;
  }
  colorVals = normalVals + nxtn * num_total_vertices_points;
  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  pickColorVals = (colorVals + nxtn * num_total_vertices_points);

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
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: OptimizePointsToVBO used with unsupported CGO ops" ENDFB(I->G);
      return 0;
      break;
    case CGO_NORMAL:
      cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
      has_normals = true;
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
      has_colors = true;
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_DRAW_ARRAYS:
      {
        const auto sp = it.cast<cgo::draw::arrays>();
	short shouldCompress = false;
	switch(sp->mode){
	case GL_POINTS:
	  shouldCompress = true;
	default:
	  break;
	}
	/*	TODO : DO WE NEED TO COMPENSATE FOR THIS? if (!has_normals && arrays & CGO_NORMAL_ARRAY){
	  arrays = arrays ^ CGO_NORMAL_ARRAY;
	  narrays -= 1;
	  }*/
	if (shouldCompress){	
	  int cnt, nxtn = VERTEX_POS_SIZE, incr = 0;
	  float *vertexValsDA = NULL, *nxtVals = NULL, *colorValsDA = NULL, *normalValsDA = NULL;
	  float *pickColorValsDA = NULL, *pickColorValsTMP;

	  nxtVals = vertexValsDA = sp->floatdata;

	  for (cnt=0; cnt<sp->nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  if (sp->arraybits & CGO_NORMAL_ARRAY){
	    has_normals = true;
	    nxtVals = normalValsDA = vertexValsDA + (nxtn*sp->nverts);
	  }
	  if (sp->arraybits & CGO_COLOR_ARRAY){
	    has_colors = true;
	    nxtVals = colorValsDA = nxtVals + (nxtn*sp->nverts);
	    nxtn = VERTEX_COLOR_SIZE;
	  }
	  if (sp->arraybits & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*sp->nverts);
	    pickColorValsDA = nxtVals + sp->nverts;
	    nxtn = VERTEX_PICKCOLOR_SIZE;
	  }
	  pickColorValsTMP = pickColorVals + (idxpl * 2);
	  switch (sp->mode){
	  case GL_POINTS:
	    for (cnt = 0; cnt < sp->nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  }
	  vpl += sp->nverts;
        }
      }
      break;
    }
    ok &= !I->G->Interrupt;
  }
  if (ok){
    short arrays = CGO_VERTEX_ARRAY | CGO_PICK_COLOR_ARRAY;
    short nsz = 12;
    GLenum ntp = GL_FLOAT;
    bool nnorm = GL_FALSE;
    if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
      nsz = 3;
      ntp = GL_BYTE;
      nnorm = GL_TRUE;
    }

    short csz = 4;
    GLenum ctp = GL_FLOAT;
    bool cnorm = GL_FALSE;
    if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
      csz = 1;
      ctp = GL_UNSIGNED_BYTE;
      cnorm = GL_TRUE;
    }

    VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();

    BufferDataDesc bufData =
      { { "a_Vertex", GL_FLOAT, 3, sizeof(float) * num_total_vertices_points * 3, vertexVals, GL_FALSE } };

    if (has_normals){
      bufData.push_back( { "a_Normal", ntp,      3, (size_t)(num_total_vertices_points * nsz), normalVals, nnorm } );
    }
    if (has_colors){
      bufData.push_back( { "a_Color",  ctp,      4, sizeof(float) * num_total_vertices_points * csz, colorVals, cnorm } );
    }
    ok = vbo->bufferData(std::move(bufData));

    if (ok && has_colors){
      arrays |= CGO_COLOR_ARRAY;
    }

    size_t vboid = vbo->get_hash_id();

    if (ok){
      float *newPickColorVals ;
      if (addshaders)
	CGOEnable(cgo, GL_DEFAULT_SHADER);
      newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(GL_POINTS, arrays, num_total_vertices_points, vboid);
      CHECKOK(ok, newPickColorVals);
      if (ok && addshaders)
	ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
      if (!newPickColorVals)
	I->G->ShaderMgr->freeGPUBuffer(vboid);
      if (ok)
	memcpy(newPickColorVals + num_total_vertices_points, pickColorVals, num_total_vertices_points * 2 * sizeof(float));
      *has_draw_buffer = true;
    } else {
      I->G->ShaderMgr->freeGPUBuffer(vboid);
    }
  }
  FreeP(vertexVals);
  return ok;
  /* END GL_POINTS */
  //    printf("num_total_vertices_points=%d\n", num_total_vertices_points);
}

static
void FixPickColorsForLine(float *pick1, float *pick2){
  unsigned int p1 = CGO_get_uint(pick1);
  unsigned int p2 = CGO_get_uint(pick2);
  int b1 = CGO_get_int(pick1 + 1);
  int b2 = CGO_get_int(pick2 + 1);
  if (p1 != p2 || b1 != b2){
    // if the pick colors are different, then pick the first one
    CGO_put_uint(pick1, p2);
    CGO_put_int(pick1 + 1, b2);
  }
}

static
void FixPickColorsForTriangle(float *pick1, float *pick2, float *pick3){
  unsigned int p1 = CGO_get_uint(pick1);
  unsigned int p2 = CGO_get_uint(pick2);
  unsigned int p3 = CGO_get_uint(pick3);
  int b1 = CGO_get_int(pick1 + 1);
  int b2 = CGO_get_int(pick2 + 1);
  int b3 = CGO_get_int(pick3 + 1);
  if (p1 != p2 || p1 != p3 || p2 != p3 ||
      b1 != b2 || b1 != b3 || b2 != b3){
    // right now, if the pick colors are different, then pick majority, otherwise, pick first one
    if (p1 == p2 && b1 == b2){
      CGO_put_uint(pick3, p1);
      CGO_put_int(pick3 + 1, b1);
    } else if (p1 == p3 && b1 == b3){
      CGO_put_uint(pick2, p1);
      CGO_put_int(pick2 + 1, b1);
    } else if (p2 == p3 && b2 == b3){
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
    case CGO_DRAW_BUFFERS_INDEXED:
    {
      const float * newpc = pc;
      if (addtocgo)
        addtocgo->add_to_cgo(op, newpc);
    }
      break;
    case CGO_NORMAL:
      cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
      has_normals = true;
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
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
    case CGO_DRAW_ARRAYS:
      {
        const auto sp = it.cast<cgo::draw::arrays>();
	short shouldCompress = false;
	switch(sp->mode){
	case GL_TRIANGLE_FAN:
	case GL_TRIANGLE_STRIP:
	case GL_TRIANGLES:
	  shouldCompress = true;
	default:
	  break;
	}
	if (shouldCompress){
          int cnt, incr = 0;
          const float* nxtVals = sp->floatdata;
          const float *vertexValsDA = nullptr, *colorValsDA = nullptr,
                      *normalValsDA = nullptr, *pickColorValsDA = nullptr,
                      *accessibilityValsDA = nullptr;

          assert(sp->arraybits & CGO_VERTEX_ARRAY);
          vertexValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_POS_SIZE;

	  for (cnt=0; cnt<sp->nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  if (sp->arraybits & CGO_NORMAL_ARRAY){
            has_normals = true;
	    normalValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_NORMAL_SIZE;
	  }

	  if (sp->arraybits & CGO_COLOR_ARRAY){
            has_colors = true;
	    colorValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_COLOR_SIZE;
	  }
	  if (sp->arraybits & CGO_PICK_COLOR_ARRAY){
            nxtVals += VERTEX_PICKCOLOR_RGBA_SIZE * sp->nverts;
            pickColorValsDA = nxtVals;
            nxtVals += VERTEX_PICKCOLOR_INDEX_SIZE * sp->nverts;
	  }
          float* pickColorValsTMP = pickColorVals + (idxpl * 2);
	  if (sp->arraybits & CGO_ACCESSIBILITY_ARRAY){
	    if (!(*ambient_occlusion) && incr){
	      for (cnt=0; cnt<incr;cnt++){
		/* if ambient_occlusion, need to fill in the array */
		accessibilityVals[cnt] = 1.f;
	      }
	    }
	    (*ambient_occlusion) = 1;
	    accessibilityValsDA = nxtVals;
	    nxtVals += VERTEX_ACCESSIBILITY_SIZE * sp->nverts;
            has_accessibility = true;
	  } else {
	    if (*ambient_occlusion){
	      for (cnt=incr; cnt<incr+sp->nverts;cnt++){
		/* if ambient_occlusion, need to fill in the array */
		accessibilityVals[cnt] = 1.f;
	      }
	    }
	  }
	  switch (sp->mode){
	  case GL_TRIANGLES:
	    for (cnt = 0; ok && cnt < sp->nverts; cnt++){
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++,
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
              if (incr && (incr % 3) == 0){
                FixPickColorsForTriangle(pickColorValsTMP + (incr-3) * 2,
                                         pickColorValsTMP + (incr-2) * 2,
                                         pickColorValsTMP + (incr-1) * 2);
              }
	      idxpl++; pl += 3; plc += 4;
	      ok &= !G->Interrupt;
	    }
	    break;
	  case GL_TRIANGLE_STRIP:
	    {
	      short flip = 0;
	      for (cnt = 2; ok && cnt < sp->nverts; cnt++){
		SetVertexValuesForVBO(G, cgo, pl, plc, cnt - (flip ? 0 : 2), incr++,
				      vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				      vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				      pickColorValsTMP, accessibilityVals, accessibilityValsDA);
		idxpl++; pl += 3; plc += 4;
		SetVertexValuesForVBO(G, cgo, pl, plc, cnt-1, incr++,
				      vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				      vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				      pickColorValsTMP, accessibilityVals, accessibilityValsDA);
		idxpl++; pl += 3; plc += 4;
		SetVertexValuesForVBO(G, cgo, pl, plc, cnt - (flip ? 2 : 0) , incr++,
				      vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				      vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				      pickColorValsTMP, accessibilityVals, accessibilityValsDA);
                FixPickColorsForTriangle(pickColorValsTMP + (incr-3) * 2,
                                         pickColorValsTMP + (incr-2) * 2,
                                         pickColorValsTMP + (incr-1) * 2);
		idxpl++; pl += 3; plc += 4;
		ok &= !G->Interrupt;
		flip = !flip;
	      }
	    }
	    break;
	  case GL_TRIANGLE_FAN:
	    for (cnt = 2; ok && cnt < sp->nverts; cnt++){
	      SetVertexValuesForVBO(G, cgo, pl, plc, 0, incr++,
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt - 1, incr++,
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++,
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA,
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals,
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
              FixPickColorsForTriangle(pickColorValsTMP + (incr-3) * 2,
                                       pickColorValsTMP + (incr-2) * 2,
                                       pickColorValsTMP + (incr-1) * 2);
	      idxpl++; pl += 3; plc += 4;
	      ok &= !G->Interrupt;
	    }
	    break;
	  }
	  vpl += sp->nverts;
	}
      }
      break;
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
      cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
      break;
    case CGO_TEX_COORD:
      cgo->texture[0] = *pc; cgo->texture[1] = *(pc + 1);
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      break;
    case CGO_BEGIN:
      {
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
	    cgo->texture[0] = *pc; cgo->texture[1] = *(pc + 1);
	    break;
	  case CGO_COLOR:
	    cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
	    break;
	  case CGO_ALPHA:
	    cgo->alpha = *pc;
	    break;
	  case CGO_VERTEX:
	    {
	      switch (mode){
	      case GL_TRIANGLES:
		{
		  int vpl = pl * 3, tpl = pl * 2, cpl = pl * 4;
		  vertexVals[vpl] = *pc;
		  vertexVals[vpl + 1] = *(pc + 1);
		  vertexVals[vpl + 2] = *(pc + 2);
		  texcoordVals[tpl] = cgo->texture[0];
		  texcoordVals[tpl+1] = cgo->texture[1];
		  if (colorValsUC){
		    colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]);
		    colorValsUC[cpl+1] = CLIP_COLOR_VALUE(cgo->color[1]);
		    colorValsUC[cpl+2] = CLIP_COLOR_VALUE(cgo->color[2]);
		    colorValsUC[cpl+3] = CLIP_COLOR_VALUE(cgo->alpha);
		  } else {
		    colorVals[cpl] = cgo->color[0];
		    colorVals[cpl+1] = cgo->color[1];
		    colorVals[cpl+2] = cgo->color[2];
		    colorVals[cpl+3] = cgo->alpha;
		  }
		  pl++;
		}
		break;
	      case GL_TRIANGLE_STRIP:
		{
		  int vpl = pl * 3, tpl = pl * 2, cpl = pl * 4;
		  if (ipl < 3){
		    vertexVals[vpl] = *pc; vertexVals[vpl + 1] = *(pc + 1); vertexVals[vpl + 2] = *(pc + 2);
		    texcoordVals[tpl] = cgo->texture[0]; texcoordVals[tpl+1] = cgo->texture[1];
		    if (colorValsUC){
		      colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]); colorValsUC[cpl+1] = CLIP_COLOR_VALUE(cgo->color[1]);
		      colorValsUC[cpl+2] = CLIP_COLOR_VALUE(cgo->color[2]); colorValsUC[cpl+3] = CLIP_COLOR_VALUE(cgo->alpha);
		    } else {
		      colorVals[cpl] = cgo->color[0]; colorVals[cpl+1] = cgo->color[1];
		      colorVals[cpl+2] = cgo->color[2]; colorVals[cpl+3] = cgo->alpha;
		    }
		    pl++; ipl++;
		  } else {
		    vertexVals[vpl] = vertexVals[vpl-6]; vertexVals[vpl + 1] = vertexVals[vpl-5]; vertexVals[vpl + 2] = vertexVals[vpl-4];
		    texcoordVals[tpl] = texcoordVals[tpl-4]; texcoordVals[tpl+1] = texcoordVals[tpl-3];
		    if (colorValsUC){
		      colorValsUC[cpl] = colorValsUC[cpl-8]; colorValsUC[cpl+1] = colorValsUC[cpl-7];
		      colorValsUC[cpl+2] = colorValsUC[cpl-6]; colorValsUC[cpl+3] = colorValsUC[cpl-5];
		    } else {
		      colorVals[cpl] = colorVals[cpl-8]; colorVals[cpl+1] = colorVals[cpl-7];
		      colorVals[cpl+2] = colorVals[cpl-6]; colorVals[cpl+3] = colorVals[cpl-5];
		    }
		    pl++; vpl+=3; tpl+=2; cpl+=4; ipl++;
		    vertexVals[vpl] = vertexVals[vpl-6]; vertexVals[vpl + 1] = vertexVals[vpl-5]; vertexVals[vpl + 2] = vertexVals[vpl-4];
		    texcoordVals[tpl] = texcoordVals[tpl-4]; texcoordVals[tpl+1] = texcoordVals[tpl-3];
		    if (colorValsUC){
		      colorValsUC[cpl] = colorValsUC[cpl-8]; colorValsUC[cpl+1] = colorValsUC[cpl-7];
		      colorValsUC[cpl+2] = colorValsUC[cpl-6]; colorValsUC[cpl+3] = colorValsUC[cpl-5];
		    } else {
		      colorVals[cpl] = colorVals[cpl-8]; colorVals[cpl+1] = colorVals[cpl-7];
		      colorVals[cpl+2] = colorVals[cpl-6]; colorVals[cpl+3] = colorVals[cpl-5];
		    }
		    pl++; vpl+=3; tpl+=2; cpl+=4; ipl++;
		    vertexVals[vpl] = *pc; vertexVals[vpl + 1] = *(pc + 1); vertexVals[vpl + 2] = *(pc + 2);
		    texcoordVals[tpl] = cgo->texture[0]; texcoordVals[tpl+1] = cgo->texture[1];
		    if (colorValsUC){
		      colorValsUC[cpl] = CLIP_COLOR_VALUE(cgo->color[0]); colorValsUC[cpl+1] = CLIP_COLOR_VALUE(cgo->color[1]);
		      colorValsUC[cpl+2] = CLIP_COLOR_VALUE(cgo->color[2]); colorValsUC[cpl+3] = CLIP_COLOR_VALUE(cgo->alpha);
		    } else {
		      colorVals[cpl] = cgo->color[0]; colorVals[cpl+1] = cgo->color[1];
		      colorVals[cpl+2] = cgo->color[2]; colorVals[cpl+3] = cgo->alpha;
		    }
		    pl++; ipl++;
		  }
		}
		break;
	      case GL_TRIANGLE_FAN:
	      default:
		printf("CGOProcessScreenCGOtoArrays: WARNING: mode=%d not implemented yet GL_LINES=%d GL_LINE_STRIP=%d GL_LINE_LOOP=%d\n", mode, GL_LINES, GL_LINE_STRIP, GL_LINE_LOOP);
		break;
	      }
	      nverts++;
	    }
	  }
	}
      }
      break;
    }
  }
  return true;
}

bool CGOOptimizeToVBONotIndexed(CGO ** I) {
  CGO *cgo = CGOOptimizeToVBONotIndexed(*I, 0, true, NULL);
  CGOFree(*I);
  *I = cgo;
  return (cgo != NULL);
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
 * @param returnedData ???
 *
 * @return New converted CGO with use_shader=true
 */
CGO *CGOOptimizeToVBONotIndexed(const CGO * I, int est, bool addshaders, float **returnedData)
{
  auto G = I->G;

  std::unique_ptr<CGO> I_begin_end_combined;
  if (I->has_begin_end) {
    I_begin_end_combined.reset(CGOCombineBeginEnd(I));
    I = I_begin_end_combined.get();
    assert(!I->has_begin_end);
  }

  int num_total_vertices = 0, num_total_indexes = 0, num_total_vertices_lines = 0, num_total_indexes_lines = 0,
    num_total_vertices_points = 0;
  short has_draw_buffer = false;
  float min[3] = { FLT_MAX, FLT_MAX, FLT_MAX }, max[3] = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
  int ambient_occlusion = 0;
  int ok = true;
  auto cgo = CGONew(G);

  CGOCountNumVertices(I, &num_total_vertices, &num_total_indexes,
                         &num_total_vertices_lines, &num_total_indexes_lines,
                         &num_total_vertices_points);
  if (num_total_vertices_points>0){
    if (!OptimizePointsToVBO(I, cgo, num_total_vertices_points, min, max, &has_draw_buffer, addshaders)){
      CGOFree(cgo);
      return NULL;
    }
  }
  if (num_total_indexes>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals;
    float *pickColorVals, *accessibilityVals = 0;
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;

    cgo->alpha = 1.f;
    cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;

    unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE + VERTEX_ACCESSIBILITY_SIZE;
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal) ? 1 : VERTEX_NORMAL_SIZE;
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
    auto const tot = size_t(num_total_indexes) * mul;

    vertexVals = pymol::malloc<float>(tot);
    if (!vertexVals){
      PRINTFB(G, FB_CGO, FB_Errors)
      "%s-Error(%d): vertexVals could not be allocated (tot=%zu)\n", __func__,
          __LINE__, tot ENDFB(G);
      CGOFree(cgo);
      return (NULL);
    }
    normalVals = vertexVals + 3 * num_total_indexes;
    unsigned nxtn = VERTEX_NORMAL_SIZE;
    if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      nxtn = 1;
    }
    colorVals = normalVals + nxtn * num_total_indexes;
    if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
      colorValsUC = (uchar*) colorVals;
      nxtn = 1;
    } else {
      nxtn = 4;
    }
    pickColorVals = (colorVals + nxtn * num_total_indexes);
    nxtn = 3;
    accessibilityVals = pickColorVals + nxtn * num_total_indexes;

    bool has_normals = false, has_colors = false, has_accessibility = false;
    ok = CGOProcessCGOtoArrays(I, cgo, cgo, min, max, &ambient_occlusion,
        vertexVals, normalVals, normalValsC, colorVals, colorValsUC,
        pickColorVals, accessibilityVals, has_normals, has_colors,
        has_accessibility);
    if (!ok){
      if (!G->Interrupt)
        PRINTFB(G, FB_CGO, FB_Errors) "ERROR: CGOProcessCGOtoArrays() could not allocate enough memory\n" ENDFB(G);
      FreeP(vertexVals);      
      CGOFree(cgo);
      return (NULL);
    }
    if (ok){
      short nsz = VERTEX_NORMAL_SIZE * 4;
      GLenum ntp = GL_FLOAT;
      bool nnorm = GL_FALSE;
      if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)) {
        nsz = VERTEX_NORMAL_SIZE;
        ntp = GL_BYTE;
        nnorm = GL_TRUE;
      }

      short csz = 4;
      GLenum ctp = GL_FLOAT;
      bool cnorm = GL_FALSE;
      if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)) {
        csz = 1;
        ctp = GL_UNSIGNED_BYTE;
        cnorm = GL_TRUE;
      }

      VertexBuffer * vbo = G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL);
      BufferDataDesc bufData =
        { { "a_Vertex", GL_FLOAT, 3, sizeof(float) * num_total_indexes * 3, vertexVals, GL_FALSE } };
      if (has_normals){
        bufData.push_back( { "a_Normal", ntp,      VERTEX_NORMAL_SIZE, (size_t)(num_total_indexes * nsz), normalVals, nnorm } );
      }
      if (has_colors){
        bufData.push_back( { "a_Color",  ctp,      4, sizeof(float) * num_total_indexes * csz, colorVals, cnorm } );
      }
      if (has_accessibility){
        bufData.push_back( { "a_Accessibility", GL_FLOAT, 1, sizeof(float) * num_total_indexes, accessibilityVals, GL_FALSE } );
      }
      ok = vbo->bufferData(std::move(bufData));

      size_t vboid = vbo->get_hash_id();
      // picking VBO: generate a buffer twice the size needed, for each picking pass
      VertexBuffer * pickvbo = G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
      ok = pickvbo->bufferData({
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes, 0, GL_TRUE ),
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes, 0, GL_TRUE )
        });
      size_t pickvboid = pickvbo->get_hash_id();

      if (ok){
	float *newPickColorVals ;
	int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY;
	if (ambient_occlusion){
	  arrays |= CGO_ACCESSIBILITY_ARRAY;
	}
	if (addshaders)
	  CGOEnable(cgo, GL_DEFAULT_SHADER_WITH_SETTINGS);
	newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(GL_TRIANGLES, arrays, num_total_indexes, vboid, pickvboid);
	if (ok && addshaders)
	  ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
	CHECKOK(ok, newPickColorVals);
	if (!newPickColorVals){
          G->ShaderMgr->freeGPUBuffer(pickvboid);
          G->ShaderMgr->freeGPUBuffer(vboid);
	}
	if (!ok){
	  PRINTFB(G, FB_CGO, FB_Errors) "CGOOptimizeToVBONotIndexedWithReturnedData: ERROR: CGODrawBuffersNotIndexed() could not allocate enough memory\n" ENDFB(G);
	  FreeP(vertexVals);
	  CGOFree(cgo);
	  return (NULL);
	}
	memcpy(newPickColorVals + num_total_indexes, pickColorVals, num_total_indexes * 2 * sizeof(float));
	has_draw_buffer = true;
      } else {
        G->ShaderMgr->freeGPUBuffer(vboid);
        G->ShaderMgr->freeGPUBuffer(pickvboid);
      }
    }
    if (ok && returnedData){
      returnedData[0] = vertexVals;
    } else {
      FreeP(vertexVals);
    }
  }
  if (ok && num_total_indexes_lines>0){
    bool has_color = false, has_normals = false;
    float *vertexVals = 0, *colorVals = 0, *normalVals;
    float *pickColorVals;
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;

    cgo->alpha = 1.f;
    cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;

    unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal) ? 1 : VERTEX_NORMAL_SIZE;
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
    auto const tot = size_t(num_total_indexes_lines) * mul;

    vertexVals = pymol::malloc<float>(tot);
    if (!vertexVals){
      PRINTFB(G, FB_CGO, FB_Errors)
      "%s-Error(%d): vertexVals could not be allocated (tot=%zu)\n", __func__,
          __LINE__, tot ENDFB(G);
      CGOFree(cgo);
      return (NULL);
    }
    normalVals = vertexVals + 3 * num_total_indexes_lines;
    nxtn = 3;
    if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      nxtn = 1;
    }

    colorVals = normalVals + nxtn * num_total_indexes_lines;
    if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      nxtn = 1;
    } else {
      nxtn = 4;
    }
    pickColorVals = (colorVals + nxtn * num_total_indexes_lines);

    for (auto it = I->begin(); !it.is_stop(); ++it) {
      auto pc = it.data();
      const auto op = it.op_code();

      switch (op) {
      case CGO_SPECIAL:
      case CGO_RESET_NORMAL:
	{
	  const float *newpc = pc;
          cgo->add_to_cgo(op, newpc);
	}
	break;
      case CGO_NORMAL:
        has_normals = true;
	cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
	break;
      case CGO_COLOR:
        has_color = true;
	cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
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
      case CGO_DRAW_ARRAYS:
      {
        auto sp = it.cast<cgo::draw::arrays>();
        short shouldCompress = false;
	switch(sp->mode){
	case GL_LINE_LOOP:
	case GL_LINE_STRIP:
	case GL_LINES:
	  shouldCompress = true;
	default:
	  break;
	}

	if (shouldCompress){
          int cnt, incr = 0;
          const float* nxtVals = sp->floatdata;
          const float *vertexValsDA = nullptr, *colorValsDA = nullptr,
                      *normalValsDA = nullptr, *pickColorValsDA = nullptr;

          assert(sp->arraybits & CGO_VERTEX_ARRAY);
          vertexValsDA = nxtVals;
          nxtVals += sp->nverts * VERTEX_POS_SIZE;

	  for (cnt=0; cnt<sp->nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  if (sp->arraybits & CGO_NORMAL_ARRAY){
            has_normals = true;
            normalValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_NORMAL_SIZE;
	  }

	  if (sp->arraybits & CGO_COLOR_ARRAY){
            has_color = true;
            colorValsDA = nxtVals;
            nxtVals += sp->nverts * VERTEX_COLOR_SIZE;
	  }
	  if (sp->arraybits & CGO_PICK_COLOR_ARRAY){
            nxtVals += VERTEX_PICKCOLOR_RGBA_SIZE * sp->nverts;
            pickColorValsDA = nxtVals;
            nxtVals += VERTEX_PICKCOLOR_INDEX_SIZE * sp->nverts;
	  }
          float* pickColorValsTMP = pickColorVals + (idxpl * 2);
	  switch (sp->mode){
	  case GL_LINES:
	    for (cnt = 0; cnt < sp->nverts; cnt++){
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP);
              if (incr && (incr % 2) == 0){
                FixPickColorsForLine(pickColorValsTMP + (incr-2) * 2,
                                     pickColorValsTMP + (incr-1) * 2);
              }
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  case GL_LINE_STRIP:
	    for (cnt = 1; cnt < sp->nverts; cnt++){
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt-1, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP);
              FixPickColorsForLine(pickColorValsTMP + (incr-2) * 2,
                                   pickColorValsTMP + (incr-1) * 2);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  case GL_LINE_LOOP:
	    for (cnt = 1; cnt < sp->nverts; cnt++){
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt-1, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(G, cgo, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP);
              FixPickColorsForLine(pickColorValsTMP + (incr-2) * 2,
                                   pickColorValsTMP + (incr-1) * 2);
	      idxpl++; pl += 3; plc += 4;
	    }
	    SetVertexValuesForVBO(G, cgo, pl, plc, 0, incr++, 
				  vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				  vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				  pickColorValsTMP);
	    idxpl++; pl += 3; plc += 4;
	    SetVertexValuesForVBO(G, cgo, pl, plc, sp->nverts-1, incr++, 
				  vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				  vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				  pickColorValsTMP);
            FixPickColorsForLine(pickColorValsTMP + (incr-2) * 2,
                                 pickColorValsTMP + (incr-1) * 2);
	    idxpl++; pl += 3; plc += 4;
	    break;
	  }

	  //	  pl += 3 * nverts;
	  //	  plc += 4 * nverts;
	  vpl += sp->nverts;
	}
      }
	break;
      }
    }
    {
      short nsz = VERTEX_NORMAL_SIZE * 4;
      GLenum ntp = GL_FLOAT;
      bool nnorm = GL_FALSE;
      if (SettingGet<int>(G, cSetting_cgo_shader_ub_normal)){
        nsz = VERTEX_NORMAL_SIZE;
        ntp = GL_BYTE;
        nnorm = GL_TRUE;
      }

      short csz = 4;
      GLenum ctp = GL_FLOAT;
      bool cnorm = GL_FALSE;
      if (SettingGet<int>(G, cSetting_cgo_shader_ub_color)){
        csz = 1;
        ctp = GL_UNSIGNED_BYTE;
        cnorm = GL_TRUE;
      }

      VertexBuffer * vbo = G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL);
      BufferDataDesc bufData =
        { { "a_Vertex", GL_FLOAT, 3, sizeof(float) * num_total_indexes_lines * 3, vertexVals, GL_FALSE } };

      if (has_normals){
        bufData.push_back( { "a_Normal", ntp,      VERTEX_NORMAL_SIZE, (size_t)(num_total_indexes_lines * nsz), normalVals, nnorm } );
      }
      if (has_color){
        bufData.push_back( { "a_Color",  ctp,      4, sizeof(float) * num_total_indexes_lines * csz, colorVals, cnorm } );
      }
      ok = vbo->bufferData(std::move(bufData));
      size_t vboid = vbo->get_hash_id();

      // picking VBO: generate a buffer twice the size needed, for each picking pass
      VertexBuffer * pickvbo = G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
      ok &= pickvbo->bufferData({
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes_lines, 0, GL_TRUE ),
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes_lines, 0, GL_TRUE )
        });
      size_t pickvboid = pickvbo->get_hash_id();

      if (ok){
	float *newPickColorVals ;
	if (addshaders)
	  CGOEnable(cgo, GL_DEFAULT_SHADER_WITH_SETTINGS);
	CGODisable(cgo, GL_SHADER_LIGHTING);
	newPickColorVals = cgo->add<cgo::draw::buffers_not_indexed>(GL_LINES, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, num_total_indexes_lines, vboid, pickvboid);
	if (ok && addshaders)
	  ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
	CHECKOK(ok, newPickColorVals);
	if (!ok){
	  PRINTFB(G, FB_CGO, FB_Errors) "CGOOptimizeToVBONotIndexedWithReturnedData: ERROR: CGODrawBuffersNotIndexed() could not allocate enough memory\n" ENDFB(G);
	  FreeP(vertexVals);
	  CGOFree(cgo);
	  if (!newPickColorVals) {
            G->ShaderMgr->freeGPUBuffer(pickvboid);
            G->ShaderMgr->freeGPUBuffer(vboid);
          }
	  return (NULL);
	}
	memcpy(newPickColorVals + num_total_indexes_lines, pickColorVals, num_total_indexes_lines * 2 * sizeof(float));
	has_draw_buffer = true;
      } else {
        G->ShaderMgr->freeGPUBuffer(pickvboid);
        G->ShaderMgr->freeGPUBuffer(vboid);
      }
    }
    if (ok && returnedData){
      returnedData[1] = vertexVals;
    } else {
      FreeP(vertexVals);
    }
  }

  if (ok && (num_total_vertices>0 || num_total_vertices_lines>0 || num_total_vertices_points>0)){
    ok &= CGOBoundingBox(cgo, min, max);
  }

  if (ok)
    ok &= CGOStop(cgo);

  if (!ok){
    CGOFree(cgo);
    return nullptr;
  }

  cgo->has_draw_buffers |= has_draw_buffer;
  cgo->use_shader = true;
  cgo->cgo_shader_ub_color = SettingGet<int>(G, cSetting_cgo_shader_ub_color);
  cgo->cgo_shader_ub_normal = SettingGet<int>(G, cSetting_cgo_shader_ub_normal);
  return (cgo);
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
CGO *CGOOptimizeToVBOIndexed(const CGO * I, int est,
    const float *color, bool addshaders, bool embedTransparencyInfo)
{
  auto G = I->G;
  CGO *cgo;

  std::unique_ptr<CGO> I_begin_end_combined;
  if (I->has_begin_end) {
    I_begin_end_combined.reset(CGOCombineBeginEnd(I));
    I = I_begin_end_combined.get();
    assert(!I->has_begin_end);
  }

  int num_total_vertices = 0, num_total_indexes = 0, num_total_vertices_lines = 0, num_total_indexes_lines = 0,
    num_total_vertices_points = 0;
  short has_draw_buffer = false;
  float min[3] = { FLT_MAX, FLT_MAX, FLT_MAX }, max[3] = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
  int ok = true;

  CGOCountNumVertices(I, &num_total_vertices, &num_total_indexes,
		      &num_total_vertices_lines, &num_total_indexes_lines,
		      &num_total_vertices_points);

  cgo = CGONewSized(I->G, I->c + est);
  CHECKOK(ok, cgo);
  if (ok){
    if (color){
      cgo->color[0] = color[0]; cgo->color[1] = color[1]; cgo->color[2] = color[2];
      cgo->alpha = color[3];
    } else {
      cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;
      cgo->alpha = 1.f;
    }
  }

#if defined(_PYMOL_IOS) && !defined(_WEBGL)
  if (num_total_indexes > MAX_INDICES_FOR_IOS || num_total_indexes_lines > MAX_INDICES_FOR_IOS){
    PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() VBO Memory Limitation: The requested \n       operation requires a larger buffer than PyMOL currently allows. \n       The operation has not entirely completed successfully.\n" ENDFB(I->G);
    firePyMOLLimitationWarning();
    CGOFree(cgo);
    return NULL;
  }
#endif

  if (num_total_vertices_points>0){
    /* This does not need to be indexed (for now) */
    if (!OptimizePointsToVBO(I, cgo, num_total_vertices_points, min, max, &has_draw_buffer, addshaders)){
      CGOFree(cgo);
      return NULL;
    }
  }

  if (num_total_vertices>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals, *accessibilityVals = 0;
    float *pickColorVals;
    GL_C_INT_TYPE *vertexIndices; 
    short vertexIndicesAllocated = 0;
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, nxtn;
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;
    short ambient_occlusion = 0;
    float *sumarray = NULL;
    int n_data = 0;

    if (embedTransparencyInfo){
      int n_tri = num_total_indexes / 3;
      int bytes_to_allocate = 2 * num_total_indexes * sizeof(GL_C_INT_TYPE) + // vertexIndicesOriginal, vertexIndices
	3 * num_total_indexes * sizeof(float) +  // 3 * for sum
	n_tri * sizeof(float) + 2 * n_tri * sizeof(int) + 256 * sizeof(int);    // z_value (float * n_tri), ix (n_tri * int), sort_mem ((n_tri + 256) * int)
      // round to 4 byte words for the length of the CGO
      n_data = bytes_to_allocate / 4 + (((bytes_to_allocate % 4) == 0) ? 0 : 1) ;
    }
    vertexIndices = pymol::calloc<GL_C_INT_TYPE>(num_total_indexes);
    vertexIndicesAllocated = 1;
    if (!vertexIndices){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() vertexIndices could not be allocated\n" ENDFB(I->G);	
      CGOFree(cgo);
      return (NULL);
    }

    unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE + VERTEX_ACCESSIBILITY_SIZE;
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal) ? 1 : VERTEX_NORMAL_SIZE;
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
    auto const tot = size_t(num_total_vertices) * mul;

    vertexVals = pymol::malloc<float>(tot);
    if (!vertexVals){
      PRINTFB(G, FB_CGO, FB_Errors)
      "%s-Error(%d): vertexVals could not be allocated (tot=%zu)\n", __func__,
          __LINE__, tot ENDFB(G);
      CGOFree(cgo);
      return (NULL);
    }
    normalVals = vertexVals + 3 * num_total_vertices;
    nxtn = 3;
    if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      nxtn = 1;
    }

    colorVals = normalVals + nxtn * num_total_vertices;
    if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      nxtn = 1;
    } else {
      nxtn = 4;
    }
    pickColorVals = (colorVals + nxtn * num_total_vertices);
    accessibilityVals = pickColorVals + 3 * num_total_vertices;

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto pc = it.data();
      const auto op = it.op_code();

      switch (op) {
      case CGO_NORMAL:
	cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
	break;
      case CGO_COLOR:
	cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
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
      case CGO_DRAW_ARRAYS:
	{
        auto sp = it.cast<cgo::draw::arrays>();
	short shouldCompress = false;
	switch(sp->mode){
	case GL_TRIANGLE_FAN:
	case GL_TRIANGLE_STRIP:
	case GL_TRIANGLES:
	  shouldCompress = true;
	default:
	  break;
	}
	if (shouldCompress){	
	  int cnt, nxtn = 3;
	  float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0, *normalValsDA, *accessibilityValsDA;
	  float *pickColorValsDA, *pickColorValsTMP, *accessibilityValsTMP;

	  nxtVals = vertexValsDA = sp->floatdata;
	  for (cnt=0; cnt<sp->nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  for (cnt=0; cnt<sp->nverts*3; cnt++){
	    vertexVals[pl + cnt] = vertexValsDA[cnt];
	  }
	  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
	    if (sp->arraybits & CGO_NORMAL_ARRAY){
	      nxtVals = normalValsDA = vertexValsDA + (nxtn*sp->nverts);
	      for (cnt=0; cnt<sp->nverts*3; cnt++){
		normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] = CLIP_NORMAL_VALUE(normalValsDA[cnt]);
	      }
	    } else {
	      uchar norm[3] = { CLIP_NORMAL_VALUE(cgo->normal[0]), CLIP_NORMAL_VALUE(cgo->normal[1]), CLIP_NORMAL_VALUE(cgo->normal[2]) };
	      for (cnt=0; cnt<sp->nverts*3; cnt++){
		normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] = norm[cnt%3];
	      }
	    }
	  } else {
	    if (sp->arraybits & CGO_NORMAL_ARRAY){
	      nxtVals = normalValsDA = vertexValsDA + (nxtn*sp->nverts);
	      for (cnt=0; cnt<sp->nverts*3; cnt++){
		normalVals[pl + cnt] = normalValsDA[cnt];
	      }
	    } else {
	      for (cnt=0; cnt<sp->nverts*3; cnt++){
		normalVals[pl + cnt] = cgo->normal[cnt%3];
	      }
	    }
	  }
	  nxtn = 3;
	  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
	    if (sp->arraybits & CGO_COLOR_ARRAY){
	      nxtVals = colorValsDA = nxtVals + (nxtn*sp->nverts);
	      for (cnt=0; cnt<sp->nverts*4; cnt+=4){
		colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
		colorValsUC[plc + cnt + 1] = CLIP_COLOR_VALUE(colorValsDA[cnt+1]);
		colorValsUC[plc + cnt + 2] = CLIP_COLOR_VALUE(colorValsDA[cnt+2]);
		colorValsUC[plc + cnt + 3] = CLIP_COLOR_VALUE(colorValsDA[cnt+3]);
	      }
	      nxtn = VERTEX_COLOR_SIZE;
	    } else {
	      uchar col[4] = { CLIP_COLOR_VALUE(cgo->color[0]), CLIP_COLOR_VALUE(cgo->color[1]), CLIP_COLOR_VALUE(cgo->color[2]), CLIP_COLOR_VALUE(cgo->alpha) };
	      for (cnt=0; cnt<sp->nverts*4; cnt++){
		colorValsUC[plc + cnt] = col[cnt%4];
	      }
	    }
	  } else {
	    if (sp->arraybits & CGO_COLOR_ARRAY){
	      nxtVals = colorValsDA = nxtVals + (nxtn*sp->nverts);
	      for (cnt=0; cnt<sp->nverts*4; cnt+=4){
		colorVals[plc + cnt] = colorValsDA[cnt];
		colorVals[plc + cnt + 1] = colorValsDA[cnt+1];
		colorVals[plc + cnt + 2] = colorValsDA[cnt+2];
		colorVals[plc + cnt + 3] = colorValsDA[cnt+3];
	      }
	      nxtn = VERTEX_COLOR_SIZE;
	    } else {
	      float col[4] = { cgo->color[0], cgo->color[1], cgo->color[2], cgo->alpha };
	      for (cnt=0; cnt<sp->nverts*4; cnt++){
		colorVals[plc + cnt] = col[cnt%4];
	      }
	    }
	  }
	  if (sp->arraybits & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*sp->nverts);
	    pickColorValsDA = nxtVals + sp->nverts;
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<sp->nverts; cnt++){
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	    }
	    nxtn = VERTEX_PICKCOLOR_SIZE;
	  } else {
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<sp->nverts; cnt++){
	      CGO_put_uint(pickColorValsTMP++, cgo->current_pick_color_index);
	      CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
	    }
	  }
	  if (sp->arraybits & CGO_ACCESSIBILITY_ARRAY){
	    if (!ambient_occlusion){
	      for (cnt=0; cnt<vpl; cnt++){
		accessibilityVals[cnt] = 1.f;
	      }
	    }
	    ambient_occlusion = 1;
	    nxtVals = nxtVals + (nxtn*sp->nverts);
	    accessibilityValsDA = nxtVals;
	    accessibilityValsTMP = accessibilityVals + vpl;
	    for (cnt=0; cnt<sp->nverts; cnt++){
	      accessibilityValsTMP[cnt] = accessibilityValsDA[cnt];
	    }
	  } else {
	    if (ambient_occlusion){
	      accessibilityValsTMP = accessibilityVals + vpl;
	      for (cnt=0; cnt<sp->nverts; cnt++){
		accessibilityValsTMP[cnt] = 1.f;
	      }
	    }
	  }
	  switch (sp->mode){
	  case GL_TRIANGLES:
	    for (cnt = 0; cnt < sp->nverts; cnt++){
	      vertexIndices[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_TRIANGLE_STRIP:
	    {
	      short flip = 0;
	      for (cnt = 2; cnt < sp->nverts; cnt++){
		vertexIndices[idxpl++] = vpl + cnt - (flip ? 0 : 2);
		vertexIndices[idxpl++] = vpl + cnt - 1;
		vertexIndices[idxpl++] = vpl + cnt - (flip ? 2 : 0);
		flip = !flip;
	      }
	    }
	    break;
	  case GL_TRIANGLE_FAN:
	    for (cnt = 2; cnt < sp->nverts; cnt++){
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
	}
	break;
      default:
	break;
      }
      ok &= !I->G->Interrupt;
    }
    if (sumarray){
      for (idxpl = 0; idxpl < num_total_indexes; idxpl+=3){
	add3f(&vertexVals[3 * vertexIndices[idxpl]], &vertexVals[3 * vertexIndices[idxpl+1]], sumarray);
	add3f(&vertexVals[3 * vertexIndices[idxpl+2]], sumarray, sumarray);
	sumarray += 3;
      }
    }
    if (ok) {
      short nsz = VERTEX_NORMAL_SIZE * 4;
      GLenum ntp = GL_FLOAT;
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
        nsz = VERTEX_NORMAL_SIZE;
        ntp = GL_BYTE;
      }

      short csz = 4;
      GLenum ctp = GL_FLOAT;
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
        csz = 1;
        ctp = GL_UNSIGNED_BYTE;
      }

      VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
      ok &= vbo->bufferData({
          BufferDesc( "a_Vertex",        GL_FLOAT, 3, sizeof(float) * num_total_vertices * 3, vertexVals, GL_FALSE ),
          BufferDesc( "a_Normal",        ntp,      VERTEX_NORMAL_SIZE, num_total_vertices * nsz, normalVals, GL_FALSE ),
          BufferDesc( "a_Color",         ctp,      4, sizeof(float) * num_total_vertices * csz, colorVals, GL_TRUE ),
          BufferDesc( "a_Accessibility", GL_FLOAT, 1, sizeof(float) * num_total_vertices, accessibilityVals, GL_FALSE )
        });

      IndexBuffer * ibo = I->G->ShaderMgr->newGPUBuffer<IndexBuffer>();
      ok &= ibo->bufferData({
          BufferDesc( GL_UNSIGNED_INT, sizeof(GL_C_INT_TYPE) * num_total_indexes, vertexIndices )
        });

      size_t vboid = vbo->get_hash_id();
      size_t iboid = ibo->get_hash_id();

      VertexBuffer * pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
      ok &= pickvbo->bufferData({
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes, 0, GL_TRUE ),
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes, 0, GL_TRUE )
        });
      size_t pickvboid = pickvbo->get_hash_id();

      if (ok) {
	float *newPickColorVals ;
	int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY;
	if (ambient_occlusion){
	  arrays |= CGO_ACCESSIBILITY_ARRAY;
	}
	if (addshaders)
	  CGOEnable(cgo, GL_DEFAULT_SHADER);
	newPickColorVals = cgo->add<cgo::draw::buffers_indexed>(GL_TRIANGLES, arrays, num_total_indexes, num_total_vertices, vboid, iboid, n_data, pickvboid);
	if (embedTransparencyInfo){
	  int n_tri = num_total_indexes/3;
	  float *sumarray;
	  float *sum = sumarray = newPickColorVals + num_total_vertices*3;
	  float *z_value = sum + (num_total_indexes*3);
	  int *ix = (int *)z_value + n_tri;
	  int *sort_mem = ix + n_tri;
	  GL_C_INT_TYPE *vertexIndicesOriginalTI = (GL_C_INT_TYPE *)(sort_mem + n_tri + 256);
	  
	  for (idxpl = 0; idxpl < num_total_indexes; idxpl+=3){
	    add3f(&vertexVals[3 * vertexIndices[idxpl]], &vertexVals[3 * vertexIndices[idxpl+1]], sumarray);
	    add3f(&vertexVals[3 * vertexIndices[idxpl+2]], sumarray, sumarray);
	    sumarray += 3;
	  }
	  memcpy(vertexIndicesOriginalTI, vertexIndices, sizeof(GL_C_INT_TYPE) * num_total_indexes);
	}

	if (addshaders && ok)
	  ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
	CHECKOK(ok, newPickColorVals);
	if (!newPickColorVals){
          I->G->ShaderMgr->freeGPUBuffer(pickvboid);
          I->G->ShaderMgr->freeGPUBuffer(vboid);
          I->G->ShaderMgr->freeGPUBuffer(iboid);
	}
	if (ok)
	  memcpy(newPickColorVals + num_total_vertices, pickColorVals, num_total_vertices * 2 * sizeof(float));
	has_draw_buffer = true;
      } else {
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(iboid);
      }
    }
    if (vertexIndicesAllocated)
      FreeP(vertexIndices);
    FreeP(vertexVals);
  }
  if (ok && num_total_vertices_lines>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals = NULL, *nxtVals;
    float *pickColorVals;
    GL_C_INT_TYPE *vertexIndexes; 
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, sz;
    bool hasNormals = 0;

    hasNormals = !CGOHasAnyLineVerticesWithoutNormals(I);
    vertexIndexes = pymol::malloc<GL_C_INT_TYPE>(num_total_indexes_lines);
    if (!vertexIndexes){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() vertexIndexes could not be allocated\n" ENDFB(I->G);	
      CGOFree(cgo);
      return (NULL);
    }

    unsigned mul = VERTEX_POS_SIZE + VERTEX_PICKCOLOR_SIZE;
    if (hasNormals) {
      mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_normal) ? 1 : VERTEX_NORMAL_SIZE;
    }
    mul += SettingGet<bool>(G, cSetting_cgo_shader_ub_color) ? 1 : VERTEX_COLOR_SIZE;
    auto const tot = size_t(num_total_vertices_lines) * mul;

    vertexVals = pymol::malloc<float>(tot);
    if (!vertexVals){
      PRINTFB(G, FB_CGO, FB_Errors)
      "%s-Error(%d): vertexVals could not be allocated (tot=%zu)\n", __func__,
          __LINE__, tot ENDFB(G);
      CGOFree(cgo);
      return (NULL);
    }
    nxtVals = vertexVals + VERTEX_POS_SIZE * num_total_vertices_lines;

    if (hasNormals){
      normalVals = nxtVals;
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
	normalValsC = (uchar*) normalVals;
	sz = 1;
      } else {
        sz = VERTEX_NORMAL_SIZE;
      }
      nxtVals += sz * num_total_vertices_lines;
    }

    colorVals = nxtVals;
    if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      sz = 1;
    } else {
      sz = VERTEX_COLOR_SIZE;
    }
    nxtVals += sz * num_total_vertices_lines;

    pickColorVals = nxtVals;

    for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
      const auto pc = it.data();
      const auto op = it.op_code();

      switch (op) {
      case CGO_NORMAL:
	cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
	break;
      case CGO_COLOR:
	cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
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
      case CGO_DRAW_ARRAYS:
	{
        auto sp = it.cast<cgo::draw::arrays>();
	short shouldCompress = false;
	switch(sp->mode){
	case GL_LINE_LOOP:
	case GL_LINE_STRIP:
	case GL_LINES:
	  shouldCompress = true;
	default:
	  break;
	}
	if (shouldCompress){	
	  int cnt, nxtn = 3;
	  float *vertexValsDA = 0, *nxtVals2 = 0, *colorValsDA = 0, *normalValsDA = 0;
	  float *pickColorValsDA = 0, *pickColorValsTMP;

	  nxtVals2 = vertexValsDA = sp->floatdata;
	  for (cnt=0; cnt<sp->nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  for (cnt=0; cnt<sp->nverts*3; cnt++){
	    vertexVals[pl + cnt] = vertexValsDA[cnt];
	  }
	  if (normalVals){
	    if (sp->arraybits & CGO_NORMAL_ARRAY){
	      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
		nxtVals2 = normalValsDA = nxtVals2 + (nxtn*sp->nverts);
		for (cnt=0; cnt<sp->nverts*3; cnt++){
		  normalValsC[VAR_FOR_NORMAL + cnt VAR_FOR_NORMAL_CNT_PLUS] = CLIP_NORMAL_VALUE(normalValsDA[cnt]);
		}
	      } else {
		nxtVals2 = normalValsDA = nxtVals2 + (nxtn*sp->nverts);
		for (cnt=0; cnt<sp->nverts*3; cnt++){
		  normalVals[VAR_FOR_NORMAL + cnt] = normalValsDA[cnt];
		}
	      }
	    }
	  }
	  if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
	    if (sp->arraybits & CGO_COLOR_ARRAY){
	      nxtVals2 = colorValsDA = nxtVals2 + (nxtn*sp->nverts);
	      for (cnt=0; cnt<sp->nverts*4; cnt++){
		colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
	      }
	      nxtn = 4;
	    } else {
	      uchar col[4] = { CLIP_COLOR_VALUE(cgo->color[0]), CLIP_COLOR_VALUE(cgo->color[1]), CLIP_COLOR_VALUE(cgo->color[2]), CLIP_COLOR_VALUE(cgo->alpha) };
	      for (cnt=0; cnt<sp->nverts*4; cnt++){
		colorValsUC[plc + cnt] = col[cnt%4];
	      }
	    }
	  } else {
	    if (sp->arraybits & CGO_COLOR_ARRAY){
	      nxtVals2 = colorValsDA = nxtVals2 + (nxtn*sp->nverts);
	      for (cnt=0; cnt<sp->nverts*4; cnt++){
		colorVals[plc + cnt] = colorValsDA[cnt];
	      }
	      nxtn = 4;
	    } else {
	      float col[4] = { cgo->color[0], cgo->color[1], cgo->color[2], cgo->alpha };
	      for (cnt=0; cnt<sp->nverts*4; cnt++){
		colorVals[plc + cnt] = col[cnt%4];
	      }
	    }
	  }
	  if (sp->arraybits & CGO_PICK_COLOR_ARRAY){
	    nxtVals2 = nxtVals2 + (nxtn*sp->nverts);
	    pickColorValsDA = nxtVals2 + sp->nverts;
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<sp->nverts; cnt++){
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	    }
	    nxtn = VERTEX_PICKCOLOR_SIZE;
	  } else {
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<sp->nverts; cnt++){
	      CGO_put_uint(pickColorValsTMP++, cgo->current_pick_color_index);
	      CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
	    }
	  }
	  if (idxpl + sp->nverts > num_total_indexes_lines){
	    PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() num_total_indexes_lines=%d mode=%d nverts=%d idxpl=%d\n", num_total_indexes_lines, sp->mode, sp->nverts, idxpl ENDFB(I->G);
	  }
	  switch (sp->mode){
	  case GL_LINES:
	    for (cnt = 0; cnt < sp->nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_LINE_STRIP:
	    for (cnt = 1; cnt < sp->nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt - 1;
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_LINE_LOOP:
	    for (cnt = 1; cnt < sp->nverts; cnt++){
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
	}
	break;
      case CGO_SPECIAL:
	CGOSpecial(cgo, CGO_get_int(pc));
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      short nsz = VERTEX_NORMAL_SIZE * sizeof(float);
      GLenum ntp = GL_FLOAT;
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_normal)){
        nsz = VERTEX_NORMAL_SIZE;
        ntp = GL_BYTE;
      }

      short csz = VERTEX_COLOR_SIZE;
      GLenum ctp = GL_FLOAT;
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
        csz = 1;
        ctp = GL_UNSIGNED_BYTE;
      }

      VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
      ok &= vbo->bufferData({
          BufferDesc( "a_Vertex",        GL_FLOAT, 3, sizeof(float) * num_total_vertices_lines * 3, vertexVals, GL_FALSE ),
          BufferDesc( "a_Normal",        ntp,      VERTEX_NORMAL_SIZE, num_total_vertices_lines * nsz, normalVals, GL_FALSE ),
          BufferDesc( "a_Color",         ctp,      4, sizeof(float) * num_total_vertices_lines * csz, colorVals, GL_TRUE )
        });

      IndexBuffer * ibo = I->G->ShaderMgr->newGPUBuffer<IndexBuffer>();
      ok &= ibo->bufferData({
          BufferDesc( GL_UNSIGNED_INT, sizeof(GL_C_INT_TYPE) * num_total_indexes_lines, vertexIndexes )
        });

      size_t vboid = vbo->get_hash_id();
      size_t iboid = ibo->get_hash_id();

      VertexBuffer * pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
      ok &= pickvbo->bufferData({
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes, 0, GL_TRUE ),
          BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * num_total_indexes, 0, GL_TRUE )
        });
      size_t pickvboid = pickvbo->get_hash_id();

      if (ok){
	float *newPickColorVals ;
	if (addshaders){
	  CGOEnable(cgo, GL_DEFAULT_SHADER);
	  CGODisable(cgo, GL_SHADER_LIGHTING);
	}
	newPickColorVals = cgo->add<cgo::draw::buffers_indexed>(GL_LINES,
                                                                CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY |
                                                                CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY,
                                                                num_total_indexes_lines,
                                                                num_total_vertices_lines, vboid, iboid, 0, pickvboid);
	CHECKOK(ok, newPickColorVals);
	if (addshaders && ok)
	  ok &= CGODisable(cgo, GL_DEFAULT_SHADER);
	if (!newPickColorVals) {
          I->G->ShaderMgr->freeGPUBuffer(pickvboid);
          I->G->ShaderMgr->freeGPUBuffer(vboid);
          I->G->ShaderMgr->freeGPUBuffer(iboid);
        }
	if (ok)
	  memcpy(newPickColorVals + num_total_vertices_lines, 
		 pickColorVals, num_total_vertices_lines * 2 * sizeof(float));
	has_draw_buffer = true;
      } else {
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(iboid);
      }
    }
    FreeP(vertexIndexes);
    FreeP(vertexVals);
  }
  if (ok && (num_total_vertices>0 || num_total_vertices_lines>0)){
    ok &= CGOBoundingBox(cgo, min, max);
  }

  if (ok)
    ok &= CGOStop(cgo);
  if (ok){
    if (has_draw_buffer){
      cgo->has_draw_buffers = true;
    }
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader){
      cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color); 
      cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
    }
  }
  if (!ok){
    CGOFree(cgo);
  }
  return (cgo);
}

struct OptimizeSphereData
{
  std::vector<float> vert;
  std::vector<unsigned char> color;
  std::vector<int> pickColor;
  std::vector<unsigned char> rightUpFlagsUB;
  std::vector<float> rightUpFlags;
  float min[3] = { FLT_MAX, FLT_MAX, FLT_MAX };
  float max[3] = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
  int total_vert = 0;
  int total_spheres = 0;
};

/**
 * @param[in] I input CGO
 * @param[out] cgo output CGO
 * @param[in] num_total_spheres Number of spheres to optimize
 * @param[out] leftOverCGO CGO data not relevant to spheres.
 * @pre cgo must not be NULL
 *
 * Note: To use leftOverCGO, it must already be preallocated.
 */

static OptimizeSphereData GetOptimizeSphereData(const CGO* I, CGO*& cgo, int num_total_spheres, CGO* leftOverCGO)
{
  OptimizeSphereData sphereData;
  const int tot = VERTICES_PER_SPHERE * 4 * num_total_spheres;
#if defined(PURE_OPENGL_ES_2)
  int rightup_flags[6] = { 0, 1, 3, 3, 2, 0 };
#else
  int rightup_flags[4] = { 0, 1, 3, 2 };
#endif

  auto cgo_shader_ub_flags = SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_flags);

  sphereData.vert.resize(tot);
  auto vertVals = sphereData.vert.data();

  sphereData.color.resize(tot);
  auto colorValsUB = sphereData.color.data();

  unsigned char* rightUpFlagValsUB = nullptr;
  float* rightUpFlagVals = nullptr;
  if (cgo_shader_ub_flags){
    sphereData.rightUpFlagsUB.resize(VALUES_PER_IMPOSTER_SPACE_COORD * VERTICES_PER_SPHERE * num_total_spheres);
    rightUpFlagValsUB = sphereData.rightUpFlagsUB.data();
  } else {
    sphereData.rightUpFlags.resize(VALUES_PER_IMPOSTER_SPACE_COORD * VERTICES_PER_SPHERE * num_total_spheres);
    rightUpFlagVals = sphereData.rightUpFlags.data();
  }

  bool has_picking = CGOHasOperationsOfType(I, CGO_PICK_COLOR);
  int* pickcolorVals = nullptr;
  if (has_picking){
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
      cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
      copyNormalToLeftOver = true;
    break;
    case CGO_COLOR:
      cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
      copyColorToLeftOver = true;
    break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      if (cgo->alpha < min_alpha) min_alpha = cgo->alpha;
      copyAlphaToLeftOver = true;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = CGO_get_uint(pc);
      cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      copyPickColorToLeftOver = true;
      break;
    case CGO_SPHERE:
      for (int vv=0; vv<VERTICES_PER_SPHERE; vv++) { // generate eight vertices of a bounding box for each cylinder
        vertVals[0] = *(pc);
        vertVals[1] = *(pc+1);
        vertVals[2] = *(pc+2);
        vertVals[3] = *(pc+3);
        set_min_max(sphereData.min, sphereData.max, vertVals);
        if (cgo_shader_ub_flags){
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
      if (has_picking){
        *(pickcolorVals++) = cgo->current_pick_color_index;
        *(pickcolorVals++) = cgo->current_pick_color_bond;
      }
      sphereData.total_spheres++;
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeSpheresToVBONonIndexed() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n", op ENDFB(I->G);
      break;
    case CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS:
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_SCREEN_TEXTURES_AND_POLYGONS encountered op=0x%X\n", op ENDFB(I->G);
      break;
    case CGO_DRAW_LABELS:
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_LABELS encountered op=0x%X\n", op ENDFB(I->G);
      break;
    case CGO_DRAW_TEXTURES:
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeCylindersToVBO() CGO_DRAW_TEXTURES encountered op=0x%X\n", op ENDFB(I->G);
      break;
    case CGO_DRAW_ARRAYS:
    default:
      if (!leftOverCGO)
        break;
      if (copyAlphaToLeftOver){
        copyAlphaToLeftOver = false;
        CGOAlpha(leftOverCGO, cgo->alpha);
      }
      if (copyColorToLeftOver){
        copyColorToLeftOver = false;
        CGOColor(leftOverCGO, cgo->color[0],  cgo->color[1],  cgo->color[2] );
      }
      if (copyNormalToLeftOver){
        CGONormalv(leftOverCGO, cgo->normal );
      }
      if (copyPickColorToLeftOver){
        copyPickColorToLeftOver = false;
        CGOPickColor(leftOverCGO, cgo->current_pick_color_index, cgo->current_pick_color_bond);
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
 * @pre cgo must not be NULL
 */

static bool PopulateGLBufferOptimizedSphereData(const CGO* I, CGO* cgo, bool addshaders, int num_total_spheres, const OptimizeSphereData& sphereData)
{
  bool ok = true;
  GLenum rtp = GL_FLOAT;
  short rsz  = sizeof(float);
  const void* radiusptr = sphereData.rightUpFlags.data();
  auto cgo_shader_ub_flags = SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_flags);
  if (cgo_shader_ub_flags) {
    rtp = GL_UNSIGNED_BYTE;
    rsz = sizeof(GLubyte);
    radiusptr = sphereData.rightUpFlagsUB.data();
  }

  VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
  ok &= vbo->bufferData({
      BufferDesc( "a_vertex_radius", GL_FLOAT, 4, sizeof(float) * sphereData.total_vert * 4, sphereData.vert.data(), GL_FALSE ),
      BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * sphereData.total_vert, sphereData.color.data(), GL_TRUE ),
      BufferDesc( "a_rightUpFlags", rtp, VALUES_PER_IMPOSTER_SPACE_COORD, rsz * sphereData.total_vert * VALUES_PER_IMPOSTER_SPACE_COORD, radiusptr, GL_FALSE )
    });
  size_t vboid = vbo->get_hash_id();

  VertexBuffer * pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
  ok &= pickvbo->bufferData({
      BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, 0, GL_TRUE ),
      BufferDesc( "a_Color", GL_UNSIGNED_BYTE, 4, sizeof(float) * sphereData.total_vert, GL_TRUE )
    }, 0, sizeof(float) * sphereData.total_vert * 2, 0);
  size_t pickvboid = pickvbo->get_hash_id();

  cgo->has_draw_buffers = true;
  cgo->has_draw_sphere_buffers = true;

  auto freebuffers = [vboid, pickvboid, I]() {
    I->G->ShaderMgr->freeGPUBuffer(vboid);
    I->G->ShaderMgr->freeGPUBuffer(pickvboid);
  };
  if (ok){
    if (addshaders)
      CGOEnable(cgo, GL_SPHERE_SHADER);
    auto pickcolor_data = (int*)cgo->add<cgo::draw::sphere_buffers>(sphereData.total_spheres, (cgo_shader_ub_flags ? 3 : 1), vboid, pickvboid); // always cgo_shader_ub_color
    CHECKOK(ok, pickcolor_data);
    if (ok && !sphereData.pickColor.empty()){
      memcpy(pickcolor_data, sphereData.pickColor.data(), num_total_spheres * 2 * 4);
    }
    if (ok && addshaders)
      ok &= CGODisable(cgo, GL_SPHERE_SHADER);
    if (!ok){
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

CGO *CGOOptimizeSpheresToVBONonIndexed(const CGO * I, int est, bool addshaders, CGO *leftOverCGO)
{
  bool ok = true;
  int num_total_spheres = CGOCountNumberOfOperationsOfType(I, CGO_SPHERE);

  if (num_total_spheres <= 0) {
    return nullptr;
  }
  auto cgo = CGONewSized(I->G, I->c + est);
  auto sphereData = GetOptimizeSphereData(I, cgo, num_total_spheres, leftOverCGO);

  if (ok && sphereData.total_spheres > 0) {
    ok = PopulateGLBufferOptimizedSphereData(I, cgo, addshaders, num_total_spheres, sphereData);
  }

  if (ok && num_total_spheres>0) {
    ok &= CGOBoundingBox(cgo, sphereData.min, sphereData.max);
  }

  if (ok)
    ok &= CGOStop(cgo);

  if (ok){
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader){
      cgo->cgo_shader_ub_color = true;
      cgo->cgo_shader_ub_normal = SettingGet<bool>(cgo->G, cSetting_cgo_shader_ub_normal);
    }
  }
  if (!ok){
    CGOFree(cgo);
  }
  return (cgo);
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
CGO *CGOSimplify(const CGO * I, int est, short sphere_quality, bool stick_round_nub)
{
  auto G = I->G;
  int ok = true;
  if (sphere_quality < 0){
    sphere_quality = SettingGet_i(I->G, NULL, NULL, cSetting_cgo_sphere_quality);
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
    case CGO_SHADER_CYLINDER:
      {
	float v2[3];
        int cap = CGO_get_int(pc + 7);
        cCylCap fcap = cap1_from_cyl_shader_bits(cap);
        cCylCap bcap = cap2_from_cyl_shader_bits(cap);
	add3f(pc, pc + 3, v2);
        ok &= CGOSimpleCylinder(cgo, pc, v2, *(pc + 6), 0, 0, cgo->alpha, cgo->alpha, (cap & cCylShaderInterpColor),
                                fcap, bcap, nullptr, stick_round_nub);
      }
      break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
      {
        auto cyl = it.cast<cgo::draw::shadercylinder2ndcolor>();
	float v1[3];
        int cap = cyl->cap;
        cCylCap fcap = cap1_from_cyl_shader_bits(cap);
        cCylCap bcap = cap2_from_cyl_shader_bits(cap);
        Pickable pickcolor2 = { cyl->pick_color_index, cyl->pick_color_bond };
        float color1[3] = { cgo->color[0], cgo->color[1], cgo->color[2] };
	add3f(pc, pc + 3, v1);
        float mid[3];
        mult3f(cyl->axis, .5f, mid);
        add3f(cyl->origin, mid, mid);
        float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : cgo->alpha;
        if (cap & cCylShaderInterpColor){
          ok &= CGOSimpleCylinder(cgo, cyl->origin, v1, cyl->tube_size, color1, cyl->color2, cgo->alpha, alpha2, true, bcap, fcap, &pickcolor2, stick_round_nub);
        } else {
          ok &= CGOColorv(cgo, color1);
          ok &= CGOSimpleCylinder(cgo, cyl->origin, mid, cyl->tube_size, color1, NULL, cgo->alpha, alpha2, false, fcap, cCylCap::None, nullptr, stick_round_nub);
          ok &= CGOColorv(cgo, cyl->color2);
          ok &= CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
          ok &= CGOSimpleCylinder(cgo, mid, v1, cyl->tube_size, cyl->color2, NULL, cgo->alpha, alpha2, false, cCylCap::None, bcap, nullptr, stick_round_nub);
        }
      }
      break;
    case CGO_CYLINDER:
      {
        auto cyl = it.cast<cgo::draw::cylinder>();
        ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true, cCylCap::Flat, cCylCap::Flat, nullptr, stick_round_nub);
      }
      break;
    case CGO_CONE:
      ok &= CGOSimpleCone(cgo, pc, pc + 3, *(pc + 6), *(pc + 7), pc + 8, pc + 11,
          static_cast<cCylCap>(int(pc[14])),
          static_cast<cCylCap>(int(pc[15])));
      break;
    case CGO_SAUSAGE:
      ok &= CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10, cgo->alpha, cgo->alpha, true, cCylCap::Round, cCylCap::Round, nullptr, stick_round_nub);
      break;
    case CGO_CUSTOM_CYLINDER:
      {
        auto cyl = it.cast<cgo::draw::custom_cylinder>();
        ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true,
            cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
      }
      break;
    case CGO_CUSTOM_CYLINDER_ALPHA:
      {
        auto cyl = it.cast<cgo::draw::custom_cylinder_alpha>();
        ok &= CGOSimpleCylinder(cgo, *cyl, cyl->color1[3], cyl->color2[3], true,
            cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
      }
      break;
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
    case CGO_BEGIN:
      {
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
	    if (!nverts){
	      hasFirstColor = 1;
	      firstColor[0] = pc[0]; firstColor[1] = pc[1]; firstColor[2] = pc[2];
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
	    if (!nverts){
	      hasFirstAlpha = 1;
	      firstAlpha = cgo->alpha;
	    } else {
	      hasFirstAlpha = 0;
	      damode |= CGO_COLOR_ARRAY;
	    }
	  }
	}

	if (nverts>0 && !err){
	  int pl = 0, plc = 0, pla = 0;
	  float *vertexVals, *tmp_ptr;
	  float *normalVals = 0, *colorVals = 0, *nxtVals = 0, *pickColorVals = 0, *accessibilityVals = 0;
	  short notHaveValue = 0, nxtn = 3;
	  if (hasFirstAlpha || hasFirstColor){
	    if (hasFirstAlpha){
	      CGOAlpha(cgo, firstAlpha);
	    }
	    if (hasFirstColor){
	      CGOColorv(cgo, firstColor);
	    }
	  }
	  nxtVals = vertexVals = cgo->add<cgo::draw::arrays>(mode, damode, nverts);
          RETURN_VAL_IF_FAIL(vertexVals, nullptr);
	  if (damode & CGO_NORMAL_ARRAY){
	    nxtVals = normalVals = vertexVals + (nxtn*nverts);
	  }
	  if (damode & CGO_COLOR_ARRAY){
	    nxtVals = colorVals = nxtVals + (nxtn*nverts);
	    nxtn = VERTEX_COLOR_SIZE;
	  }
	  if (damode & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorVals = nxtVals + nverts;
	    nxtn = VERTEX_PICKCOLOR_SIZE;
	  }
	  if (damode & CGO_ACCESSIBILITY_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
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
	      normalVals[pl] = pc[0]; normalVals[pl+1] = pc[1]; normalVals[pl+2] = pc[2];
	      notHaveValue &= ~CGO_NORMAL_ARRAY;
	      break;
	    case CGO_COLOR:
	      if (colorVals){
		colorVals[plc] = pc[0]; colorVals[plc+1] = pc[1]; 
		colorVals[plc+2] = pc[2]; colorVals[plc+3] = cgo->alpha;
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
	      if (notHaveValue & CGO_NORMAL_ARRAY){
		if (pl){
		  tmp_ptr = &normalVals[pl-3];
		  normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];		
		} else {
		  copy3f(cgo->normal, &normalVals[pl]);
		}
	      }
	      if (notHaveValue & CGO_COLOR_ARRAY){
		if (plc){
		  tmp_ptr = &colorVals[plc-4];
		  colorVals[plc] = tmp_ptr[0]; colorVals[plc+1] = tmp_ptr[1]; 
		  colorVals[plc+2] = tmp_ptr[2];	colorVals[plc+3] = tmp_ptr[3];
		} else {
		  copy3f(cgo->color, &colorVals[plc]);
		  colorVals[plc+3] = cgo->alpha;
		}
	      }
	      if (pickColorVals){
		CGO_put_uint(pickColorVals + pla * 2, cgo->current_pick_color_index);
		CGO_put_int(pickColorVals + pla * 2 + 1, cgo->current_pick_color_bond);
	      }
	      if (accessibilityVals){
		accessibilityVals[pla] = cgo->current_accessibility;
	      }
	      vertexVals[pl++] = pc[0]; vertexVals[pl++] = pc[1]; vertexVals[pl++] = pc[2];
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
      }
      break;
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
CGO *CGOSimplifyNoCompress(const CGO * I, int est, short sphere_quality, bool stick_round_nub)
{
  CGO *cgo;

  int ok = true;
  if (sphere_quality < 0){
    sphere_quality = SettingGet_i(I->G, NULL, NULL, cSetting_cgo_sphere_quality);
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
    case CGO_SHADER_CYLINDER:
      {
	float v2[3];
        int cap = CGO_get_int(pc + 7);
        cCylCap fcap = cap1_from_cyl_shader_bits(cap);
        cCylCap bcap = cap2_from_cyl_shader_bits(cap);
	add3f(pc, pc + 3, v2);
        ok &= CGOSimpleCylinder(cgo, pc, v2, *(pc + 6), 0, 0, cgo->alpha, cgo->alpha, (cap & cCylShaderInterpColor),
                                fcap, bcap, nullptr, stick_round_nub);
      }
      break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
      {
        auto cyl = it.cast<cgo::draw::shadercylinder2ndcolor>();
	float v1[3];
        int cap = cyl->cap;
        cCylCap fcap = cap1_from_cyl_shader_bits(cap);
        cCylCap bcap = cap2_from_cyl_shader_bits(cap);
        Pickable pickcolor2 = { cyl->pick_color_index, cyl->pick_color_bond };
        float color1[3] = { cgo->color[0], cgo->color[1], cgo->color[2] };
	add3f(cyl->origin, cyl->axis, v1);
        float mid[3];
        mult3f(cyl->axis, .5f, mid);
        add3f(cyl->origin, mid, mid);
        float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : cgo->alpha;
        if (cap & cCylShaderInterpColor){
          ok &= CGOSimpleCylinder(cgo, cyl->origin, v1, cyl->tube_size, color1, cyl->color2, cgo->alpha, alpha2, true, bcap, fcap, &pickcolor2, stick_round_nub);
        } else {
          ok &= CGOColorv(cgo, color1);
          ok &= CGOSimpleCylinder(cgo, cyl->origin, mid, cyl->tube_size, color1, nullptr, cgo->alpha, alpha2, false, fcap, cCylCap::None, nullptr, stick_round_nub);
          ok &= CGOColorv(cgo, cyl->color2);
          ok &= CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
          ok &= CGOSimpleCylinder(cgo, mid, v1, cyl->tube_size, cyl->color2, nullptr, cgo->alpha, alpha2, false, cCylCap::None, bcap, nullptr, stick_round_nub);
        }
      }
      break;
    case CGO_CYLINDER:
      {
        auto cyl = it.cast<cgo::draw::cylinder>();
        ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true, cCylCap::Flat, cCylCap::Flat, nullptr, stick_round_nub);
      }
      break;
    case CGO_CONE:
      ok &= CGOSimpleCone(cgo, pc, pc + 3, *(pc + 6), *(pc + 7), pc + 8, pc + 11,
          static_cast<cCylCap>(int(pc[14])),
          static_cast<cCylCap>(int(pc[15])));
      break;
    case CGO_SAUSAGE:
      ok &= CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10, cgo->alpha, cgo->alpha, true, cCylCap::Round, cCylCap::Round, nullptr, stick_round_nub);
      break;
    case CGO_CUSTOM_CYLINDER:
      {
        auto cyl = it.cast<cgo::draw::custom_cylinder>();
        ok &= CGOSimpleCylinder(cgo, *cyl, cgo->alpha, cgo->alpha, true,
            cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
      }
      break;
    case CGO_CUSTOM_CYLINDER_ALPHA:
      {
        auto cyl = it.cast<cgo::draw::custom_cylinder_alpha>();
        ok &= CGOSimpleCylinder(cgo, *cyl, cyl->color1[3], cyl->color2[3], true,
            cyl->get_cap1(), cyl->get_cap2(), nullptr, stick_round_nub);
      }
      break;
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
	PRINTFB(I->G, FB_CGO, FB_Errors) "CGOSimplifyNoCompress-Error: CGO_DRAW_BUFFERS_INDEXED encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
	PRINTFB(I->G, FB_CGO, FB_Errors) "CGOSimplifyNoCompress-Error: CGO_DRAW_BUFFERS_NOT_INDEXED encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_SPHERE_BUFFERS:
        PRINTFB(I->G, FB_CGO, FB_Errors) "CGOSimplifyNoCompress-Error: CGO_DRAW_SPHERE_BUFFERS encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_CYLINDER_BUFFERS:
        PRINTFB(I->G, FB_CGO, FB_Errors) "CGOSimplifyNoCompress-Error: CGO_DRAW_CYLINDER_BUFFERS encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_LABELS:
	PRINTFB(I->G, FB_CGO, FB_Errors) "CGOSimplifyNoCompress-Error: CGO_DRAW_LABELS encountered\n" ENDFB(I->G);
      break;
    case CGO_DRAW_TEXTURES:
	PRINTFB(I->G, FB_CGO, FB_Errors) "CGOSimplifyNoCompress-Error: CGO_DRAW_TEXTURES encountered\n" ENDFB(I->G);
      break;
    case CGO_BEGIN:
      cgo->has_begin_end = true;
    default:
      cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  if (ok){
    ok &= CGOStop(cgo);
  } 
  if (!ok){
    CGOFree(cgo);
  }
  return (cgo);
}


CGO *CGOOptimizeTextures(const CGO * I, int est)
{
  CGO *cgo = NULL;
  int num_total_textures;
  int ok = true;
  num_total_textures = CGOCountNumberOfOperationsOfType(I, CGO_DRAW_TEXTURE);
  //  printf("CGOOptimizeTextures: num_total_textures=%d\n", num_total_textures);
  if (num_total_textures){
    float *worldPos, *screenValues, *textExtents, *pickColorVals;
    int place3 = 0, place2 = 0;
    worldPos = pymol::malloc<float>(num_total_textures * 18);
    if (!worldPos){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeTextures() worldPos could not be allocated\n" ENDFB(I->G);
      return NULL;
    }
    screenValues = pymol::malloc<float>(num_total_textures * 18);
    if (!screenValues){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeTextures() screenValues could not be allocated\n" ENDFB(I->G);
      FreeP(worldPos);
      return NULL;
    }
    textExtents = pymol::malloc<float>(num_total_textures * 12);
    if (!textExtents){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeTextures() textExtents could not be allocated\n" ENDFB(I->G);
      FreeP(screenValues);
      FreeP(worldPos);
      return NULL;
    }
    pickColorVals = pymol::malloc<float>(num_total_textures * 12); /* pick index and bond */
    if (!pickColorVals){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeTextures() pickColorVals could not be allocated\n" ENDFB(I->G);
      FreeP(textExtents);
      FreeP(screenValues);
      FreeP(worldPos);
      return NULL;
    }

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
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeTextures() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n", op ENDFB(I->G);	
	break;
      case CGO_DRAW_TEXTURE:
	{
	  float screenMin[3], screenMax[3], textExtent[4];
	  copy3f(pc, &worldPos[place3]);
	  copy3f(pc, &worldPos[place3+3]);
	  copy3f(pc, &worldPos[place3+6]);
	  copy3f(pc, &worldPos[place3+9]);
	  copy3f(pc, &worldPos[place3+12]);
	  copy3f(pc, &worldPos[place3+15]);
	  copy3f(pc + 3, screenMin);
	  copy3f(pc + 6, screenMax);
	  copy4f(pc + 9, textExtent);
	  copy3f(screenMin, &screenValues[place3]);
	  copy3f(screenMin, &screenValues[place3+3]);
	  copy3f(screenMin, &screenValues[place3+6]);
	  copy3f(screenMin, &screenValues[place3+9]);
	  copy3f(screenMin, &screenValues[place3+12]);
	  copy3f(screenMax, &screenValues[place3+15]);
	  screenValues[place3+4] = screenMax[1];
	  screenValues[place3+6] = screenMax[0];
	  screenValues[place3+10] = screenMax[1];
	  screenValues[place3+12] = screenMax[0];
	  screenValues[place3+17] = screenMin[2];
	  place3 += 18;
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[0]; textExtents[place2++] = textExtent[1];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[0]; textExtents[place2++] = textExtent[3];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[2]; textExtents[place2++] = textExtent[1];
	  CGO_put_int(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_uint(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[0]; textExtents[place2++] = textExtent[3];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[2]; textExtents[place2++] = textExtent[1];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[2]; textExtents[place2++] = textExtent[3];
	}
	break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL);
      ok &= vbo->bufferData({
          BufferDesc( "attr_worldpos", GL_FLOAT, 3, sizeof(float) * num_total_textures * 18, worldPos, GL_FALSE ),
          BufferDesc( "attr_screenoffset", GL_FLOAT, 3, sizeof(float) * num_total_textures * 18, screenValues, GL_FALSE ),
          BufferDesc( "attr_texcoords", GL_FLOAT, 3, sizeof(float) * num_total_textures * 18, textExtents, GL_FALSE )
        });
      size_t vboid = vbo->get_hash_id();

      if (ok) {
	float *pickArray = cgo->add<cgo::draw::textures>(num_total_textures, vboid);
	CHECKOK(ok, pickArray);
	if (!pickArray)
	  I->G->ShaderMgr->freeGPUBuffer(vboid);
	if (ok)
	  memcpy(pickArray + num_total_textures * 6, pickColorVals, num_total_textures * 12 * sizeof(float));
	if (ok)
	  ok &= CGOStop(cgo);
      } else {
	I->G->ShaderMgr->freeGPUBuffer(vboid);
      }
      if (!ok){
	CGOFree(cgo);
      }
    }
    FreeP(worldPos);
    FreeP(screenValues);
    FreeP(textExtents);
    FreeP(pickColorVals);
  }
  return cgo;
}

CGO *CGOConvertToLabelShader(const CGO *I, CGO * addTo){
  /* Lines that pass in two vertices per line */
  PyMOLGlobals *G = I->G;

  AttribDataOp world_pos_op =
    { { CGO_DRAW_LABEL,       1, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::label, world_pos),           0 } };
  AttribDataOp screen_offset_op =
    { { CGO_DRAW_LABEL,       2, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::label, screen_world_offset), 0 } };
  AttribDataOp screen_min_op =
    { { CGO_DRAW_LABEL,       3, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::label, screen_min),          0 } };
  AttribDataOp screen_max_op =
    { { CGO_DRAW_LABEL,       4, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::label, screen_max),          0 } };
  AttribDataOp text_extent_op =
    { { CGO_DRAW_LABEL,       5, FLOAT2_TO_FLOAT2, offsetof(cgo::draw::label, text_extent),         0 } };
  AttribDataOp relative_mode_op =
    { { CGO_DRAW_LABEL,       6, FLOAT_TO_FLOAT,   offsetof(cgo::draw::label, relative_mode),       0 } };
  AttribDataOp target_pos_op =
    { { CGO_DRAW_LABEL,       7, FLOAT3_TO_FLOAT3, offsetof(cgo::draw::label, target_pos),          6 } };

  AttribDataDesc attrDesc = { { "attr_worldpos",          GL_FLOAT, 3, GL_FALSE, world_pos_op },
                              { "attr_targetpos",         GL_FLOAT, 3, GL_FALSE, target_pos_op },
                              { "attr_screenoffset",      GL_FLOAT, 3, GL_FALSE, screen_min_op },
                              { "attr_texcoords",         GL_FLOAT, 2, GL_FALSE, text_extent_op },
                              { "attr_screenworldoffset", GL_FLOAT, 3, GL_FALSE, screen_offset_op },
                              { "attr_relative_mode",     GL_FLOAT, 1, GL_FALSE, relative_mode_op } };

  auto ComputeScreenValues = [](void * varData, const float * pc, void * screenData, int idx) {
    auto sp = reinterpret_cast<const cgo::draw::label *>(pc);
    const auto& smin = sp->screen_min;
    const auto& smax = sp->screen_max;
    float * v = reinterpret_cast<float *>(varData);
    switch (idx) {
    case 0:
    v[0] = smin[0]; v[1] = smin[1]; v[2] = smin[2];
    break;
    case 1:
    v[0] = smin[0]; v[1] = smax[1]; v[2] = smin[2];
    break;
    case 2:
    v[0] = smax[0]; v[1] = smin[1]; v[2] = smin[2];
    break;
    case 3:
    v[0] = smin[0]; v[1] = smax[1]; v[2] = smin[2];
    break;
    case 4:
    v[0] = smax[0]; v[1] = smin[1]; v[2] = smin[2];
    break;
    case 5:
    v[0] = smax[0]; v[1] = smax[1]; v[2] = smin[2];
    break;
    };
  };

  auto ComputeTexCoords = [](void * varData, const float * pc, void * discard, int idx) {
    auto sp = reinterpret_cast<const cgo::draw::label *>(pc);
    float * v = reinterpret_cast<float *>(varData);
    const auto& te = sp->text_extent;
    static struct { int x, y; } const idxs[6] = {
      { 0, 1 },
      { 0, 3 },
      { 2, 1 },
      { 0, 3 },
      { 2, 1 },
      { 2, 3 }
    };
    v[0] = te[idxs[idx].x];
    v[1] = te[idxs[idx].y];
  };

  attrDesc[1].attrOps[0].funcDataConversions.push_back({ ComputeScreenValues, nullptr, "attr_screenoffset" });
  attrDesc[1].attrOps[0].funcDataConversions.push_back({ ComputeTexCoords, nullptr, "attr_texcoords" });

  uchar pickdata[4] = { 0, 0, 0, 0 };
  addTo->add<cgo::draw::vertex_attribute_4ub>(G->ShaderMgr->GetAttributeUID("attr_pickcolor"), pickdata);

  AttribDataOp pickOp = { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 } };
  AttribDataDesc pickDesc = { { "attr_pickcolor", GL_UNSIGNED_BYTE, 4, GL_TRUE, pickOp } };
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES, VertexBuffer::INTERLEAVED, true);
}

CGO *CGOOptimizeLabels(const CGO * I, int est, bool addshaders)
{
  CGO *cgo = NULL;
  int num_total_labels;
  int ok = true;
  num_total_labels = CGOCountNumberOfOperationsOfType(I, CGO_DRAW_LABEL);
  if (num_total_labels){
    float *targetPos, *worldPos, *screenValues, *screenWorldValues, *textExtents, *pickColorVals;
    float *relativeMode;
    int place3 = 0, place2 = 0, place = 0;
    worldPos = pymol::malloc<float>(num_total_labels * 6 * 17); 
    if (!worldPos){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeLabels() worldPos could not be allocated\n" ENDFB(I->G);
      return NULL;
    }
    screenValues = worldPos + (num_total_labels * 18);
    targetPos = screenValues + (num_total_labels * 18);
    screenWorldValues = targetPos + (num_total_labels * 18);
    textExtents = screenWorldValues + (num_total_labels * 18);
    pickColorVals = textExtents + (num_total_labels * 12); /* pick index and bond */
    relativeMode = (float *)(pickColorVals + (num_total_labels * 12));
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
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeLabels() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n", op ENDFB(I->G);	
	break;
      case CGO_DRAW_LABEL:
	{
	  float screenWorldOffset[3], screenMin[3], screenMax[3], textExtent[4];
	  copy3f(pc, &worldPos[place3]);
	  copy3f(pc, &worldPos[place3+3]);
	  copy3f(pc, &worldPos[place3+6]);
	  copy3f(pc, &worldPos[place3+9]);
	  copy3f(pc, &worldPos[place3+12]);
	  copy3f(pc, &worldPos[place3+15]);
	  copy3f(pc + 3, screenWorldOffset);
	  copy3f(pc + 6, screenMin);
	  copy3f(pc + 9, screenMax);
	  copy4f(pc + 12, textExtent);
	  copy3f(screenWorldOffset, &screenWorldValues[place3]);
	  copy3f(&screenWorldValues[place3], &screenWorldValues[place3+3]);
	  copy3f(&screenWorldValues[place3], &screenWorldValues[place3+6]);
	  copy3f(&screenWorldValues[place3], &screenWorldValues[place3+9]);
	  copy3f(&screenWorldValues[place3], &screenWorldValues[place3+12]);
	  copy3f(&screenWorldValues[place3], &screenWorldValues[place3+15]);
	  copy3f(screenMin, &screenValues[place3]);
	  copy3f(screenMin, &screenValues[place3+3]);
	  copy3f(screenMin, &screenValues[place3+6]);
	  copy3f(screenMin, &screenValues[place3+9]);
	  copy3f(screenMin, &screenValues[place3+12]);
	  copy3f(screenMax, &screenValues[place3+15]);
	  screenValues[place3+4] = screenMax[1];
	  screenValues[place3+6] = screenMax[0];
	  screenValues[place3+10] = screenMax[1];
	  screenValues[place3+12] = screenMax[0];
	  screenValues[place3+17] = screenMin[2];
	  copy3f(pc + 17, &targetPos[place3]);
	  copy3f(pc + 17, &targetPos[place3+3]);
	  copy3f(pc + 17, &targetPos[place3+6]);
	  copy3f(pc + 17, &targetPos[place3+9]);
	  copy3f(pc + 17, &targetPos[place3+12]);
	  copy3f(pc + 17, &targetPos[place3+15]);
	  place3 += 18;
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[0]; textExtents[place2++] = textExtent[1];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[0]; textExtents[place2++] = textExtent[3];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[2]; textExtents[place2++] = textExtent[1];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[0]; textExtents[place2++] = textExtent[3];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[2]; textExtents[place2++] = textExtent[1];
	  CGO_put_uint(pickColorVals + place2, cgo->current_pick_color_index);
	  CGO_put_int(pickColorVals + place2 + 1, cgo->current_pick_color_bond);
	  textExtents[place2++] = textExtent[2]; textExtents[place2++] = textExtent[3];
	  {
	    uchar rM = (uchar)*(pc + 16);
	    relativeMode[place++] = rM;
	    relativeMode[place++] = rM;
	    relativeMode[place++] = rM;
	    relativeMode[place++] = rM;
	    relativeMode[place++] = rM;
	    relativeMode[place++] = rM;
	  }
	}
	break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      // Static Vertex Data
      VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL);
      ok &= vbo->bufferData({
          BufferDesc( "attr_worldpos", GL_FLOAT, 3, sizeof(float)*18*num_total_labels, worldPos,          GL_FALSE ),
          BufferDesc( "attr_targetpos", GL_FLOAT, 3, sizeof(float)*18*num_total_labels, targetPos,         GL_FALSE ),
          BufferDesc( "attr_screenoffset", GL_FLOAT, 3, sizeof(float)*18*num_total_labels, screenValues,      GL_FALSE ),
          BufferDesc( "attr_texcoords", GL_FLOAT, 2, sizeof(float)*12*num_total_labels, textExtents,       GL_FALSE ),
          BufferDesc( "attr_screenworldoffset", GL_FLOAT, 3, sizeof(float)*18*num_total_labels, screenWorldValues, GL_FALSE ),
          BufferDesc( "attr_relative_mode", GL_FLOAT, 1, sizeof(float)*6*num_total_labels,  relativeMode,      GL_FALSE )
        });
      size_t vboid = vbo->get_hash_id();

      VertexBuffer * pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
      ok &= pickvbo->bufferData({
          BufferDesc( "attr_pickcolor", GL_UNSIGNED_BYTE, VERTEX_COLOR_SIZE, 0, GL_TRUE ),
          BufferDesc( "attr_pickcolor", GL_UNSIGNED_BYTE, VERTEX_COLOR_SIZE, sizeof(float) * num_total_labels * 6, GL_TRUE )
        }, 0, sizeof(float) * num_total_labels * 12, 0);
      size_t pickvboid = pickvbo->get_hash_id();

      auto freebuffers = [vboid, pickvboid, I]() {
        I->G->ShaderMgr->freeGPUBuffer(vboid);
        I->G->ShaderMgr->freeGPUBuffer(pickvboid);
      };

      if (ok) {
	float *pickArray = NULL;
	if (addshaders){
	  CGOEnable(cgo, GL_LABEL_SHADER);
	}
	pickArray = cgo->add<cgo::draw::labels>(num_total_labels, vboid, pickvboid);
	if (addshaders){
	  CGODisable(cgo, GL_LABEL_SHADER);
	}
	CHECKOK(ok, pickArray);
	if (!pickArray) {
          freebuffers();
        }
	if (ok)
	  memcpy(pickArray + num_total_labels * 6, pickColorVals, num_total_labels * 12 * sizeof(float));
	if (ok)
	  ok &= CGOStop(cgo);
      } else {
        freebuffers();
      }
      if (!ok){
	CGOFree(cgo);
      }
    }
    FreeP(worldPos);
  }
  return cgo;
}

CGO *CGOOptimizeConnectors(const CGO * I, int est)
{
  CGO *cgo = NULL;
  int num_total_connectors;
  int ok = true;
  int use_geometry_shaders = SettingGetGlobal_b(I->G, cSetting_use_geometry_shaders);
  int factor = (use_geometry_shaders ? 1 : 4);
  num_total_connectors = CGOCountNumberOfOperationsOfType(I, CGO_DRAW_CONNECTOR);

  if (num_total_connectors){
    float *targetPt3d, *labelCenterPt3d, *indentFactor, *screenWorldOffset, *connectorColor, *textSize;
    float *bkgrdColor, *relExtLength, *connectorWidth;
    uchar *relativeMode, *drawBkgrd;
    uchar *isCenterPt = NULL;
    int place3 = 0, place2 = 0, place = 0;
    targetPt3d = pymol::calloc<float>(num_total_connectors * 20 * factor); /* too much, relativeMode only needs 1 byte per vertex, instead of 1 float */
    if (!targetPt3d){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeConnectors() could not be allocated\n" ENDFB(I->G);
      return NULL;
    }
    labelCenterPt3d = targetPt3d + (num_total_connectors*3 * factor);
    indentFactor = labelCenterPt3d + (num_total_connectors*3 * factor);
    screenWorldOffset = indentFactor + (num_total_connectors*2 * factor);
    connectorColor = screenWorldOffset + (num_total_connectors*3 * factor);
    textSize = connectorColor + (num_total_connectors*factor);
    relativeMode = (uchar *)(textSize + (num_total_connectors*2*factor));
    drawBkgrd = (uchar *)(relativeMode + (num_total_connectors*factor));
    bkgrdColor = (float*)(drawBkgrd + (num_total_connectors*factor));
    relExtLength = (float*)(bkgrdColor + (num_total_connectors*factor));
    connectorWidth = (float*)(relExtLength + (num_total_connectors*factor));
    if (!use_geometry_shaders)
      isCenterPt = (uchar *)(connectorWidth + (num_total_connectors*factor));
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
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeConnectors() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n", op ENDFB(I->G);	
	break;
      case CGO_DRAW_CONNECTOR:
	{
	  uchar *uc;
	  int f;
	  if (!use_geometry_shaders){
	      isCenterPt[place] = 0; isCenterPt[place+1] = 2; isCenterPt[place+2] = 2; isCenterPt[place+3] = 1;
	  }
	  copy3f(pc, &targetPt3d[place3]);
	  copy3f(pc + 3, &labelCenterPt3d[place3]);
	  copy2f(pc + 6, &indentFactor[place2]);
	  copy3f(pc + 9, &screenWorldOffset[place3]);
	  copy2f(pc + 12, &textSize[place2]);
	  relativeMode[place] = (uchar)(int)pc[17];
	  drawBkgrd[place] = (uchar)(int)pc[18];
	  uc = (uchar *)&bkgrdColor[place];
	  uc[0] = CLIP_COLOR_VALUE(*(pc+19));
	  uc[1] = CLIP_COLOR_VALUE(*(pc+20));
	  uc[2] = CLIP_COLOR_VALUE(*(pc+21));
	  uc[3] = CLIP_COLOR_VALUE(*(pc+22));
	  uc = (uchar *)&connectorColor[place];
	  uc[0] = CLIP_COLOR_VALUE(*(pc+14));
	  uc[1] = CLIP_COLOR_VALUE(*(pc+15));
	  uc[2] = CLIP_COLOR_VALUE(*(pc+16));
	  uc[3] = 255;
	  relExtLength[place] = *(pc + 8);
	  connectorWidth[place] = *(pc + 23);
	  place3 += 3; place2 += 2; place += 1;
	  for (f=1;f<factor; f++){
	    copy3f(pc, &targetPt3d[place3]);
	    copy3f(pc + 3, &labelCenterPt3d[place3]);
	    copy2f(pc + 6, &indentFactor[place2]);
	    copy3f(pc + 9, &screenWorldOffset[place3]);
	    copy2f(pc + 12, &textSize[place2]);
	    relativeMode[place] = (uchar)(int)pc[17];
	    drawBkgrd[place] = (uchar)(int)pc[18];
	    relExtLength[place] = *(pc + 8);
	    connectorWidth[place] = *(pc + 23);
	    uc = (uchar *)&bkgrdColor[place];
	    uc[0] = uc[-4];
	    uc[1] = uc[-3];
	    uc[2] = uc[-2];
	    uc[3] = uc[-1];
	    uc = (uchar *)&connectorColor[place];
	    uc[0] = uc[-4];
	    uc[1] = uc[-3];
	    uc[2] = uc[-2];
	    uc[3] = uc[-1];
	    place3 += 3; place2 += 2; place += 1;
	  }
	}
	break;
      }
      ok &= !I->G->Interrupt;
    }
    if (ok) {
      const size_t quant = factor * num_total_connectors;
      VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
      ok = vbo->bufferData({
          BufferDesc( "a_target_pt3d",       GL_FLOAT, 3, sizeof(float) * 3 * quant,     targetPt3d,        GL_FALSE ),
          BufferDesc( "a_center_pt3d",       GL_FLOAT, 3, sizeof(float) * 3 * quant,     labelCenterPt3d,   GL_FALSE ),
          BufferDesc( "a_indentFactor",      GL_FLOAT, 2, sizeof(float) * 2 * quant,     indentFactor,      GL_FALSE ),
          BufferDesc( "a_screenWorldOffset", GL_FLOAT, 3, sizeof(float) * 3 * quant,     screenWorldOffset, GL_FALSE ),
          BufferDesc( "a_textSize",          GL_FLOAT, 2, sizeof(float) * 2 * quant,     textSize,          GL_FALSE ),
          BufferDesc( "a_Color",             GL_UNSIGNED_BYTE, 4, sizeof(float) * quant, connectorColor,    GL_TRUE  ),
          BufferDesc( "a_relative_mode",     GL_UNSIGNED_BYTE, 1, sizeof(uchar) * quant, relativeMode,      GL_FALSE ),
          BufferDesc( "a_draw_flags",        GL_UNSIGNED_BYTE, 1, sizeof(uchar) * quant, drawBkgrd,         GL_FALSE ),
          BufferDesc( "a_bkgrd_color",       GL_UNSIGNED_BYTE, 4, sizeof(float) * quant, bkgrdColor,        GL_TRUE  ),
          BufferDesc( "a_rel_ext_length",    GL_FLOAT, 1, sizeof(float) * quant,         relExtLength,      GL_FALSE ),
          BufferDesc( "a_con_width",         GL_FLOAT, 1, sizeof(float) * quant,         connectorWidth,    GL_FALSE ),
          BufferDesc( "a_isCenterPt",        GL_UNSIGNED_BYTE, 1, sizeof(uchar) * quant, isCenterPt,        GL_FALSE )
        });
      size_t vboid = vbo->get_hash_id();
      if (ok) {
	cgo->add<cgo::draw::connectors>(num_total_connectors, vboid);
	if (ok)
	  ok &= CGOStop(cgo);
      }
      if (!ok){
        I->G->ShaderMgr->freeGPUBuffer(vboid);
	CGOFree(cgo);
      }
    }
    FreeP(targetPt3d);
  }
  {
    GLenum err ;
    CHECK_GL_ERROR_OK("ERROR: CGOOptimizeConnectors() end returns err=%d\n");
  }
  return cgo;
}


CGO *CGOExpandDrawTextures(const CGO * I, int est)
{
  CGO *cgo = CGONew(I->G);
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
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeTextures() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n", op ENDFB(I->G);	
      break;
    case CGO_DRAW_TEXTURE:
      {
	float screenMin[3], screenMax[3], textExtent[4];
	float alpha = cgo->alpha;
	CGOAlpha(cgo, 0.f);
	CGOColor(cgo, 0.f,0.f,0.f);
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
      }
      break;
    default:
      cgo->add_to_cgo(op, pc);
    }
    ok &= !I->G->Interrupt;
  }
  CGOStop(cgo);
  return cgo;
}

/* ======== Raytrace Renderer ======== */

int CGOGetExtent(const CGO * I, float *mn, float *mx)
{
  int result = false;

#define check_extent(v,r) {\
    if(!result) {\
      mn[0]=((*(v  ))-r); \
      mx[0]=((*(v  ))+r);  \
      mn[1]=((*(v+1))-r); \
      mx[1]=((*(v+1))+r); \
      mn[2]=((*(v+2))-r); \
      mx[2]=((*(v+2))+r); \
      result=true; \
  } else {\
       if(mn[0]>((*(v    ))-r)) mn[0]=((*(v    ))-r); \
       if(mx[0]<((*(v    ))+r)) mx[0]=((*(v    ))+r); \
       if(mn[1]>((*((v)+1))-r)) mn[1]=((*((v)+1))-r); \
       if(mx[1]<((*((v)+1))+r)) mx[1]=((*((v)+1))+r); \
       if(mn[2]>((*((v)+2))-r)) mn[2]=((*((v)+2))-r); \
       if(mx[2]<((*((v)+2))+r)) mx[2]=((*((v)+2))+r); }}

#define check_extent4(v,r) {\
    if(!result) {\
      mn[0]=((*(v  ))-r); \
      mx[0]=((*(v  ))+r);  \
      mn[1]=((*(v+1))-r); \
      mx[1]=((*(v+1))+r); \
      mn[2]=((*(v+2))-r); \
      mx[2]=((*(v+2))+r); \
      mn[3]=((*(v+3))-r); \
      mx[3]=((*(v+3))+r); \
      result=true; \
  } else {\
       if(mn[0]>((*(v    ))-r)) mn[0]=((*(v    ))-r); \
       if(mx[0]<((*(v    ))+r)) mx[0]=((*(v    ))+r); \
       if(mn[1]>((*((v)+1))-r)) mn[1]=((*((v)+1))-r); \
       if(mx[1]<((*((v)+1))+r)) mx[1]=((*((v)+1))+r); \
       if(mn[2]>((*((v)+2))-r)) mn[2]=((*((v)+2))-r); \
       if(mx[2]<((*((v)+2))+r)) mx[2]=((*((v)+2))+r); \
       if(mn[3]>((*((v)+3))-r)) mn[3]=((*((v)+3))-r); \
       if(mx[3]<((*((v)+3))+r)) mx[3]=((*((v)+3))+r); }}

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
    case CGO_DRAW_ARRAYS:
      {
        const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
	const float *pct = sp->floatdata;
	int pl;

	if (sp->arraybits & CGO_VERTEX_ARRAY){
	  for (pl = 0; pl < sp->nverts; pl++){
	    check_extent(pct, 0);
	    pct += VERTEX_POS_SIZE;
	  }
	}
      }
      break;
    case CGO_BOUNDING_BOX:
      {
	if (!result){
	  mn[0]=(*pc);
	  mn[1]=*(pc+1);
	  mn[2]=*(pc+2);
	  mx[0]=*(pc+3);
	  mx[1]=*(pc+4);
	  mx[2]=*(pc+5);
	  result = true;
	} else {
	  if(mn[0]>*pc) mn[0]=(*pc);
	  if(mn[1]>*(pc+1)) mn[1]=*(pc+1);
	  if(mn[2]>*(pc+2)) mn[2]=*(pc+2);
	  if(mx[0]<*(pc+3)) mx[0]=*(pc+3);
	  if(mx[1]<*(pc+4)) mx[1]=*(pc+4);
	  if(mx[2]<*(pc+5)) mx[2]=*(pc+5);
	}
      }
    }
  }
  return (result);
}

int CGOHasNormals(const CGO * I)
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

static int CGOQuadricToEllipsoid(const float *v, float r, const float *q,
                                 float *r_el, float *n0, float *n1, float *n2)
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

  if(xx_matrix_invert(inverse, inp_matrix, 4)) {

    /* inverse now contains Uij coefficients */
    float pradius = sqrt1f(-1 / inverse[15]);
    int n_rot;

    if(xx_matrix_jacobi_solve(e_vec, e_val, &n_rot, inverse, 4)) {
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
      if(mx < mag[1])
        mx = mag[1];
      if(mx < mag[2])
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

static int CGORenderQuadricRay(CRay * ray, float *v, float r, float *q)
{
  float r_el, n0[3], n1[3], n2[3];
  int ok = true;
  if(CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    ok &= ray->ellipsoid3fv(v, r_el, n0, n1, n2);
  return ok;
}


/* ======== Raytrace Renderer ======== */

int CGORenderRay(CGO * I, CRay * ray, RenderInfo * info, const float *color, ObjectGadgetRamp *ramp, CSetting * set1, CSetting * set2)
{
#ifdef _PYMOL_NO_RAY
  return 0;
#else
  int vc = 0;
  float linewidth = 1.0F;
  float widthscale = 0.15F;
  float lineradius, dotradius, dotwidth;
  float white[] = { 1.0, 1.0, 1.0 };
  float zee[] = { 0.0, 0.0, 1.0 };
  int ok = true;
  const float *n0 = NULL, *n1 = NULL, *n2 = NULL, *v0 = NULL, *v1 = NULL, *v2 = NULL, *c0 =
    NULL, *c1 = NULL, *c2 = NULL;
  float rampc0[3], rampc1[3], rampc2[3];
  int mode = -1;
  /* workaround; multi-state ray-trace bug */
  if (!I) {
    assert("TODO investigate" && false);
    return 0; /* not sure if it should return 0 or 1, 0 - fails, but is it a memory issue? might not be since the arg is NULL */ 
  }

  I->G->CGORenderer->alpha =
    1.0F - SettingGet_f(I->G, set1, set2, cSetting_cgo_transparency);

  widthscale = SettingGet_f(I->G, set1, set2, cSetting_cgo_ray_width_scale);

  /*  printf("debug %8.9f\n",SceneGetScreenVertexScale(I->G,zee)); */
  linewidth = SettingGet_f(I->G, set1, set2, cSetting_cgo_line_width);
  if(linewidth < 0.0F)
    linewidth = 1.0F;
  lineradius = SettingGet_f(I->G, set1, set2, cSetting_cgo_line_radius);
  dotwidth = SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_width);
  dotradius = SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_radius);
  if(lineradius < 0.0F)
    lineradius = linewidth * ray->PixelRadius / 2.0F;
  if(dotradius < 0.0F)
    dotradius = dotwidth * ray->PixelRadius / 2.0F;
  if(widthscale < 0.0F)
    widthscale = ray->PixelRadius / 2.0F;
  if(color)
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
        if(vc > 1)
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
    case CGO_SPECIAL_WITH_ARG:
      {
        float argval = *(pc + 1);
        switch (*(int*)pc){
        case LINEWIDTH_FOR_LINES:
          linewidth = argval;
          lineradius = widthscale * linewidth;
        }
      }
      break;
    case CGO_SPECIAL:
      {
        switch ((*(int*)pc)){
        case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON:
          {
            float radius = SettingGet_f(I->G, set1, set2, cSetting_ribbon_radius);
            if(radius == 0.0F) {
              float ribbon_width = SettingGet_f(I->G, set1, set2, cSetting_ribbon_width);
              float line_width = SceneGetDynamicLineWidth(info, ribbon_width);
              SceneGetDynamicLineWidth(info, line_width);
              radius = ray->PixelRadius * line_width / 2.0F;
            }
            lineradius = radius;
          }
          break;
        case LINEWIDTH_FOR_LINES:
          {
            float radius = SettingGet_f(I->G, set1, set2, cSetting_line_radius);
            if(radius <= 0.0F) {
              float line_width = SettingGet_f(I->G, set1, set2, cSetting_line_width);
              line_width = SceneGetDynamicLineWidth(info, line_width);
              radius = ray->PixelRadius * line_width / 2.0F;
            }
            lineradius = radius;
          }
          break;
        case CYLINDER_WIDTH_FOR_NONBONDED:
        case LINEWIDTH_WITH_SCALE:
          {
            float line_width = SettingGet_f(I->G, set1, set2, cSetting_line_width);
            line_width = SceneGetDynamicLineWidth(info, line_width);
            lineradius = widthscale * line_width / 2.f;
          }
          break;
        }
      }
      break;
    case CGO_COLOR:
      c0 = pc;
      ray->color3fv(c0);
      break;
    case CGO_ALPHA:
      I->G->CGORenderer->alpha = *pc;
      ray->transparentf(1.0F - *pc);
      break;
    case CGO_LINE:
      {
        auto line = reinterpret_cast<cgo::draw::line *>(pc);
        ok &= ray->sausage3fv(line->vertex1, line->vertex2, lineradius, c0, c0);
      }
      break;
    case CGO_SPLITLINE:
      {
        auto splitline = reinterpret_cast<cgo::draw::splitline *>(pc);
        float color2[] = { CONVERT_COLOR_VALUE(splitline->color2[0]),
                           CONVERT_COLOR_VALUE(splitline->color2[1]),
                           CONVERT_COLOR_VALUE(splitline->color2[2]) };
        if (splitline->flags & cgo::draw::splitline::interpolation){
          ok &= ray->sausage3fv(splitline->vertex1, splitline->vertex2, lineradius, c0, color2);
        } else {
          float mid[3];
          add3f(splitline->vertex1, splitline->vertex2, mid);
          mult3f(mid, .5f, mid);
          ok &= ray->customCylinder3fv(splitline->vertex1, mid, 
                                       lineradius, c0, c0, cCylCap::Round, cCylCap::None);
          ok &= ray->customCylinder3fv(mid, splitline->vertex2, 
                                       lineradius, color2, color2, cCylCap::None, cCylCap::Round);
        }
      }
      break;
    case CGO_VERTEX_CROSS:
      {
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
      }
      break;
    case CGO_VERTEX_BEGIN_LINE_STRIP:
    case CGO_VERTEX:
      v0 = pc;

      if (ramp){
	if (!ObjectGadgetRampInterVertex(ramp, v0, rampc0, -1)){
	  copy3f(white, rampc0);
	}
	c0 = rampc0;
      }
      switch (mode) {
      case GL_POINTS:
        ok &= ray->sphere3fv(v0, dotradius);
        break;
      case GL_LINES:
        if(vc & 0x1)
          ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
        v1 = v0;
	if (!ramp){
	  c1 = c0;
	}
        break;
      case GL_LINE_STRIP:
        if(vc){
          ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
        }
        v1 = v0;
	if (!ramp){
	  c1 = c0;
	}
        break;
      case GL_LINE_LOOP:
        if(vc)
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
	if( ((vc + 1) % 3) == 0)
          ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
        v2 = v1;
        n2 = n1;
        v1 = v0;
        n1 = n0;
	if (!ramp){
	  c2 = c1;
	  c1 = c0;
	}
        break;
      case GL_TRIANGLE_STRIP:
        if(vc > 1)
          ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
        v2 = v1;
        n2 = n1;
        v1 = v0;
        n1 = n0;
	if (!ramp){
	  c2 = c1;
	  c1 = c0;
	}
        break;
      case GL_TRIANGLE_FAN:
        if(vc > 1)
          ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
        else if(!vc) {
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
      if (ramp){
	switch (mode){
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
          static_cast<cCylCap>(int(pc[14])),
          static_cast<cCylCap>(int(pc[15])));
      break;
    case CGO_CUSTOM_CYLINDER:
      {
        auto cyl = reinterpret_cast<cgo::draw::custom_cylinder*>(pc);
        ok &= ray->customCylinder3fv(*cyl);
      }
      break;
    case CGO_CUSTOM_CYLINDER_ALPHA:
      {
        auto cyl = reinterpret_cast<cgo::draw::custom_cylinder_alpha*>(pc);
        ok &= ray->customCylinderAlpha3fv(*cyl);
      }
      break;
    case CGO_SHADER_CYLINDER:
      {
        float p2[3];
        int cap = CGO_get_int(pc + 7);
        const cCylCap cap1 = cap1_from_cyl_shader_bits(cap);
        const cCylCap cap2 = cap2_from_cyl_shader_bits(cap);
        add3f(pc, pc + 3, p2);
        ok &= ray->customCylinder3fv(pc, p2, *(pc + 6), ray->CurColor, ray->CurColor,
                                     cap1, cap2);
      }
      break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
      {
        auto cyl = reinterpret_cast<cgo::draw::shadercylinder2ndcolor*>(pc);
        float v1[3];
        int cap = cyl->cap;
        const cCylCap fcap = cap1_from_cyl_shader_bits(cap);
        const cCylCap bcap = cap2_from_cyl_shader_bits(cap);
        int colorinterp = cap & cCylShaderInterpColor;
        const float *color1 = c0;
        const float *color2 = cyl->color2;
        add3f(cyl->origin, cyl->axis, v1);
        float alpha1 = I->G->CGORenderer->alpha;
        float alpha2 = cyl->alpha >= 0.f ? cyl->alpha : alpha1;
        if (colorinterp || equal3f(color1, color2)) {
          ok &= ray->customCylinder3fv(pc, v1, cyl->tube_size, color1, color2, fcap, bcap, alpha1, alpha2);
        } else {
          float mid[3];
          mult3f(cyl->axis, .5f, mid);
          add3f(cyl->origin, mid, mid);

          ray->color3fv(c0);
          ok &= ray->customCylinder3fv(cyl->origin, mid, cyl->tube_size, color1, color1, fcap, cCylCap::None, alpha1, alpha2);
          ray->color3fv(cyl->color2);
          ok &= ray->customCylinder3fv(mid, v1, cyl->tube_size, color2, color2, cCylCap::None, bcap, alpha1, alpha2);
        }
      }
      break;
    case CGO_CYLINDER:
      {
        auto *cyl = reinterpret_cast<cgo::draw::cylinder*>(pc);
        ok &= ray->cylinder3fv(*cyl);
      }
      break;
    case CGO_SAUSAGE:
      ok &= ray->sausage3fv(pc, pc + 3, *(pc + 6), pc + 7, pc + 10);
      break;
    case CGO_TRIANGLE:
      ok &= ray->triangle3fv(pc, pc + 3, pc + 6, pc + 9, pc + 12, pc + 15, pc + 18,
			      pc + 21, pc + 24);
      break;
    case CGO_DRAW_ARRAYS:
      {
        cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
	int const mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
	float *vertexVals = sp->floatdata;
	float *normalVals = 0, *colorVals = 0;
        int offset = 0;
	if (arrays & CGO_VERTEX_ARRAY){
	  vertexVals = sp->floatdata;
          offset += nverts * VERTEX_POS_SIZE;
	}
	if (arrays & CGO_NORMAL_ARRAY){
	  normalVals = sp->floatdata + offset;
          offset += nverts * VERTEX_NORMAL_SIZE;
	}
	if (arrays & CGO_COLOR_ARRAY){
	  colorVals = sp->floatdata + offset;
          offset += nverts * VERTEX_COLOR_SIZE;
	}
	vc = 0;
	for (int v=0, pl=0, plc=0; ok && v<nverts; v++, pl+=3, plc+=4){
	  if (normalVals){
	    n0 = &normalVals[pl];
	  }
	  if (colorVals){
	    c0 = &colorVals[plc];
	    ray->color3fv(c0);
	    ray->transparentf(1.0f - c0[3]);
	  }
	  if (vertexVals){
	    v0 = &vertexVals[pl];
	  }
	  switch (mode){
	  case GL_POINTS:
	    ok &= ray->sphere3fv(v0, dotradius);
	    break;
	  case GL_LINES:
	    if(vc & 0x1)
	      ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
	    v1 = v0;
	    c1 = c0;
	    break;
	  case GL_LINE_STRIP:
	    if(vc)
	      ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
	    v1 = v0;
	    c1 = c0;
	    break;
	  case GL_LINE_LOOP:
	    if(vc)
	      ok &= ray->sausage3fv(v0, v1, lineradius, c0, c1);
	    else {
	      v2 = v0;
	      c2 = c0;
	    }
	    v1 = v0;
	    c1 = c0;
	    break;
	  case GL_TRIANGLES:
	    if( ((vc + 1) % 3) == 0)
	      ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
	    v2 = v1;
	    c2 = c1;
	    n2 = n1;
	    v1 = v0;
	    c1 = c0;
	    n1 = n0;
	    break;
	  case GL_TRIANGLE_STRIP:
	    if(vc > 1)
	      ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
	    v2 = v1;
	    c2 = c1;
	    n2 = n1;
	    v1 = v0;
	    c1 = c0;
	    n1 = n0;
	    break;
	  case GL_TRIANGLE_FAN:
	    if(vc > 1)
	      ok &= ray->triangle3fv(v0, v1, v2, n0, n1, n2, c0, c1, c2);
	    else if(!vc) {
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
      }
      break;
    default:
      break;
    }
  }

  if (ok)
    ray->transparentf(0.0F);
  return ok;
#endif
}


/* ======== GL Rendering ======== */

static int CGO_gl_begin_WARNING_CALLED = false, CGO_gl_end_WARNING_CALLED = false, CGO_gl_vertex_WARNING_CALLED = false;
static void CGO_gl_begin(CCGORenderer * I, CGO_op_data pc){ 
  if (I->use_shader){
    if (!CGO_gl_begin_WARNING_CALLED) { 
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_begin() is called but not implemented in OpenGLES\n" ENDFB(I->G);
      CGO_gl_begin_WARNING_CALLED = true; 
    }
  } else {
    int mode = CGO_get_int(*pc);
    if (I->debug)
      mode = CGOConvertDebugMode(I->debug, mode);
    glBegin(mode);
  }
}
static void CGO_gl_end(CCGORenderer * I, CGO_op_data){ 
  if (I->use_shader){
    if (!CGO_gl_end_WARNING_CALLED) {
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_end() is called but not implemented in OpenGLES\n" ENDFB(I->G);
      CGO_gl_end_WARNING_CALLED = true; 
    }
  } else {
    glEnd();
  }
}
static void CGO_gl_vertex(CCGORenderer * I, CGO_op_data v){
  if (I->use_shader){
  if (!CGO_gl_vertex_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_vertex_WARNING_CALLED = true;
  }
  } else {
    glVertex3fv(*v);
  }
}

static void CGO_gl_vertex_cross(CCGORenderer * I, CGO_op_data v){
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader){
#endif
  if (!CGO_gl_vertex_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_vertex_WARNING_CALLED = true;
  }
#ifndef PURE_OPENGL_ES_2
  } else {
    CSetting * set1 = NULL, * set2 = NULL;
    if (I->rep&&I->rep->cs) set1 = I->rep->cs->Setting.get();
    if (I->rep&&I->rep->obj) set2 = I->rep->obj->Setting.get();
    float nonbonded_size =
      SettingGet_f(I->G, set1, set2, cSetting_nonbonded_size);
    float pt[3];
    copy3f(*v, pt);
    pt[0] -= nonbonded_size;
    glVertex3fv(pt);
    pt[0] += 2 * nonbonded_size;
    glVertex3fv(pt);
    copy3f(*v, pt);
    pt[1] -= nonbonded_size;
    glVertex3fv(pt);
    pt[1] += 2 * nonbonded_size;
    glVertex3fv(pt);
    copy3f(*v, pt);
    pt[2] -= nonbonded_size;
    glVertex3fv(pt);
    pt[2] += 2 * nonbonded_size;
    glVertex3fv(pt);
  }
#endif
}

static void CGO_gl_line(CCGORenderer * I, CGO_op_data v){
#ifndef PURE_OPENGL_ES_2
  if (!I->use_shader){
    auto line = reinterpret_cast<const cgo::draw::line *>(*v);
    glVertex3fv(line->vertex1);
    glVertex3fv(line->vertex2);
  }
#endif
}

static void CGO_gl_splitline(CCGORenderer * I, CGO_op_data v){
#ifndef PURE_OPENGL_ES_2
  if (!I->use_shader){
    auto splitline = reinterpret_cast<const cgo::draw::splitline *>(*v);
    bool interpolation = splitline->flags & cgo::draw::splitline::interpolation;
    bool equal_colors = splitline->flags & cgo::draw::splitline::equal_colors;
    bool no_split_for_pick = splitline->flags & cgo::draw::splitline::no_split_for_pick;

    if (I->isPicking){
      if (no_split_for_pick){
        glVertex3fv(splitline->vertex1);
        glVertex3fv(splitline->vertex2);
      } else {
        float h[3];
        average3f(splitline->vertex1, splitline->vertex2, h);
        glVertex3fv(splitline->vertex1);
        glVertex3fv(h);
        unsigned char col[4];
        AssignNewPickColor(nullptr, I->info->pick, col, &I->rep->context,
                           splitline->index, splitline->bond);
        glColor4ubv(col);
        glVertex3fv(h);
        glVertex3fv(splitline->vertex2);
      }
    } else if (interpolation || equal_colors){
      glVertex3fv(splitline->vertex1);
      if (!equal_colors)
        glColor4ub(splitline->color2[0], splitline->color2[1], splitline->color2[2], CLIP_COLOR_VALUE(I->alpha));
      glVertex3fv(splitline->vertex2);
    } else {
      float h[3];
      average3f(splitline->vertex1, splitline->vertex2, h);
      glVertex3fv(splitline->vertex1);
      glVertex3fv(h);
      glColor4ub(splitline->color2[0], splitline->color2[1], splitline->color2[2], CLIP_COLOR_VALUE(I->alpha));
      glVertex3fv(h);
      glVertex3fv(splitline->vertex2);
    }
  }
#endif
}


static void CGO_gl_normal(CCGORenderer * I, CGO_op_data varg){
  const float* v = *varg;
  if (I->use_shader){
    glVertexAttrib3fv(VERTEX_NORMAL, v);
  } else {
    glNormal3f(v[0],v[1],v[2]);
  }
}

static void CGO_gl_draw_arrays(CCGORenderer * I, CGO_op_data pc){
  auto sp = reinterpret_cast<const cgo::draw::arrays*>(*pc);
  int mode = sp->mode, arrays = sp->arraybits, narrays = sp->narrays, nverts = sp->nverts;
  const float* data = sp->floatdata;
  (void) narrays;
#ifndef PURE_OPENGL_ES_2
  if (I->use_shader){
#endif

  if (arrays & CGO_VERTEX_ARRAY) glEnableVertexAttribArray(VERTEX_POS);
  if (arrays & CGO_NORMAL_ARRAY) glEnableVertexAttribArray(VERTEX_NORMAL);
  if (I->isPicking){
    if (arrays & CGO_PICK_COLOR_ARRAY){
      glEnableVertexAttribArray(VERTEX_COLOR);
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY)
      glEnableVertexAttribArray(VERTEX_COLOR);
  }

  if (arrays & CGO_VERTEX_ARRAY){
#ifdef _WEBGL
#else
    glVertexAttribPointer(VERTEX_POS, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, data);
#endif
    data += nverts * VERTEX_POS_SIZE;
  }
  if (arrays & CGO_NORMAL_ARRAY){
#ifdef _WEBGL
#else
    glVertexAttribPointer(VERTEX_NORMAL, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, data);
#endif
    data += nverts * VERTEX_NORMAL_SIZE;
  }
  if (I->isPicking){
    if (arrays & CGO_COLOR_ARRAY){
      data += nverts * VERTEX_COLOR_SIZE;
    }
    if (arrays & CGO_PICK_COLOR_ARRAY){
#ifdef _WEBGL
    glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
    glBufferData(GL_ARRAY_BUFFER, nverts * 4, data, GL_STATIC_DRAW);
    glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_FALSE, 0, 0);
#else
      glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_FALSE, 0, data);     
#endif
      data += nverts * VERTEX_PICKCOLOR_SIZE;
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY){
#ifdef _WEBGL
#else
      glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, data);
#endif
      data += nverts * VERTEX_COLOR_SIZE;
    }
    if (arrays & CGO_PICK_COLOR_ARRAY){
      data += nverts * VERTEX_PICKCOLOR_SIZE;
    }
  }
  if (I->debug){
    mode = CGOConvertDebugMode(I->debug, mode);
  }
  glDrawArrays(mode, 0, nverts);

  if (I->isPicking){
    if (arrays & CGO_PICK_COLOR_ARRAY){
      glDisableVertexAttribArray(VERTEX_COLOR);
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY)
      glDisableVertexAttribArray(VERTEX_COLOR);
  }
  if (arrays & CGO_VERTEX_ARRAY) glDisableVertexAttribArray(VERTEX_POS);
  if (arrays & CGO_NORMAL_ARRAY) glDisableVertexAttribArray(VERTEX_NORMAL);

#ifndef PURE_OPENGL_ES_2
  } else {

    int pl, pla, plc;
    const float *vertexVals = nullptr;
    const float *colorVals = 0, *normalVals = 0, *tmp_ptr;
    const uchar *pickColorVals = 0, *tmp_pc_ptr;
    float alpha = I->alpha;
    if (arrays & CGO_VERTEX_ARRAY){
      vertexVals = data;
      data += nverts * VERTEX_POS_SIZE;
    }
    if (arrays & CGO_NORMAL_ARRAY){
      normalVals = data;
      data += nverts * VERTEX_NORMAL_SIZE;
    }
    if (I->isPicking){
      alpha = 1.f;
      if (arrays & CGO_COLOR_ARRAY){
        data += nverts * VERTEX_COLOR_SIZE;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY){
	pickColorVals = (uchar*)data;
        data += nverts * VERTEX_PICKCOLOR_SIZE;
      }
    } else {
      if (arrays & CGO_COLOR_ARRAY){
	colorVals = data;
	data += nverts*4;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY){
        data += nverts * VERTEX_PICKCOLOR_SIZE;
      }
    }
    if (arrays & CGO_ACCESSIBILITY_ARRAY) {
      data += nverts * VERTEX_ACCESSIBILITY_SIZE;
    }

    if (I->debug){
      mode = CGOConvertDebugMode(I->debug, mode);
    }

    glBegin(mode);
    for (pl = 0, pla = 0, plc = 0; pl<nverts; pl++, pla+=3, plc+=4){
      if (pickColorVals){
        tmp_pc_ptr = &pickColorVals[plc]; /* the pick colors are saved with rgba */
        glColor4ub(tmp_pc_ptr[0], tmp_pc_ptr[1], tmp_pc_ptr[2], tmp_pc_ptr[3]);
      } else {
        if (colorVals){
          tmp_ptr = &colorVals[plc];
          glColor4f(tmp_ptr[0], tmp_ptr[1], tmp_ptr[2], alpha);	
        }
        if (normalVals){
          tmp_ptr = &normalVals[pla];
          glNormal3fv(&normalVals[pla]);
        }
      }
      if (vertexVals){
	tmp_ptr = &vertexVals[pla];
	glVertex3fv(&vertexVals[pla]);
      }
    }
    glEnd();
  }
#endif
}

static
void TransparentInfoSortIX(PyMOLGlobals * G, float *sum, float *z_value, 
			   int *ix, int n_tri, int *sort_mem, int t_mode);
static
void CGOReorderIndicesWithTransparentInfo(PyMOLGlobals * G, int nindices, 
					  size_t vbuf, int n_tri, int *ix, 
					  GL_C_INT_TYPE *vertexIndicesOriginal, 
					  GL_C_INT_TYPE *vertexIndices);

static void CGO_gl_draw_buffers_indexed(CCGORenderer * I, CGO_op_data pc){
  auto sp = reinterpret_cast<const cgo::draw::buffers_indexed*>(*pc);
  int mode = sp->mode, nindices = sp->nindices,
    nverts = sp->nverts, n_data = sp->n_data;
  size_t vboid = sp->vboid, iboid = sp->iboid;
  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(vboid);
  IndexBuffer * ibo  = I->G->ShaderMgr->getGPUBuffer<IndexBuffer>(iboid);
  GLenum err ;
  CHECK_GL_ERROR_OK("beginning of CGO_gl_draw_buffers_indexed err=%d\n");

  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();

  if (!shaderPrg){
    return;
  }

  if (I->isPicking){
    int attr_a_Color = shaderPrg->GetAttribLocation("a_Color");
    vbo->maskAttributes({ attr_a_Color });
    shaderPrg->Set1i("fog_enabled", 0);
    shaderPrg->Set1i("lighting_enabled", 0);
    if (I->use_shader){
      if (sp->pickvboid){
        VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
        pickvbo->bind(shaderPrg->id, I->pick_pass());
      } else {
        glEnableVertexAttribArray(attr_a_Color);
        glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, sp->floatdata);
      }
    }
  }
  if (n_data){
    // if transparency data, then sort it
    int n_tri = nindices/3;
    float *sum = sp->floatdata + nverts*3;
    float *z_value = sum + (nindices*3);
    int *ix = (int *)(z_value + n_tri);
    int *sort_mem = ix + n_tri;
    int t_mode;
    CSetting * set1 = NULL, * set2 = NULL;
    if (I->rep&&I->rep->cs) set1 = I->rep->cs->Setting.get();
    if (I->rep&&I->rep->obj) set2 = I->rep->obj->Setting.get();
    t_mode = SettingGet_i(I->G, set1, set2, cSetting_transparency_mode);
    if (t_mode!=3){
      GL_C_INT_TYPE *vertexIndicesOriginalTI = (GL_C_INT_TYPE *)(sort_mem + n_tri + 256);
      GL_C_INT_TYPE *vertexIndicesTI = vertexIndicesOriginalTI + nindices;
      TransparentInfoSortIX(I->G, sum, z_value, ix, n_tri, sort_mem, t_mode);
      CGOReorderIndicesWithTransparentInfo(I->G, nindices, iboid, n_tri, ix,
                                           vertexIndicesOriginalTI, vertexIndicesTI);
    }
  }

  if (I->debug){
    mode = CGOConvertDebugMode(I->debug, mode);
  }
  vbo->bind(shaderPrg->id);
  ibo->bind();

  CHECK_GL_ERROR_OK("CGO_gl_draw_buffers_indexed: before glDrawElements err=%d\n");
  glDrawElements(mode, nindices, GL_C_INT_ENUM, 0);
  CHECK_GL_ERROR_OK("CGO_gl_draw_buffers_indexed: after glDrawElements err=%d\n");

  vbo->unbind();
  ibo->unbind();

  if (I->isPicking) {
    VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
    if (pickvbo)
      pickvbo->unbind();
  }

  CHECK_GL_ERROR_OK("CGO_gl_draw_buffers_indexed: end err=%d\n");
}

static void CGO_gl_draw_buffers_not_indexed(CCGORenderer * I, CGO_op_data pc){
  const cgo::draw::buffers_not_indexed * sp = reinterpret_cast<decltype(sp)>(*pc);
  int mode = sp->mode;

  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg){
    return;
  }
  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vbo)
    return;
  if (I->isPicking){
    int attr_a_Color = shaderPrg->GetAttribLocation("a_Color");
    vbo->maskAttributes({ attr_a_Color });
    shaderPrg->Set1i("fog_enabled", 0);
    shaderPrg->Set1i("lighting_enabled", 0);
    if (I->use_shader){
      if (sp->pickvboid){
        VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
        pickvbo->bind(shaderPrg->id, I->pick_pass());
      } else {
        glEnableVertexAttribArray(attr_a_Color);
        glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, sp->floatdata);
      }
    }
  }

  if (I->debug){
    mode = CGOConvertDebugMode(I->debug, mode);
  }

  vbo->bind(shaderPrg->id);
  glDrawArrays(mode, 0, sp->nverts);
  vbo->unbind();

  if (I->isPicking) {
    VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
    if (pickvbo)
      pickvbo->unbind();
  }
}

static void CGO_gl_mask_attribute_if_picking(CCGORenderer * I, CGO_op_data pc){
  if (I->isPicking){
    const cgo::draw::mask_attribute_if_picking * sp = reinterpret_cast<decltype(sp)>(*pc);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    if (!shaderPrg){
      return;
    }
    VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
    if (!vbo)
      return;
    int loc = shaderPrg->GetAttribLocation(I->G->ShaderMgr->GetAttributeName(sp->attr_lookup_idx));
    vbo->maskAttribute(loc);
  }
}

static void CGO_gl_bind_vbo_for_picking(CCGORenderer * I, CGO_op_data pc){
  if (I->isPicking){
    const cgo::draw::bind_vbo_for_picking * sp = reinterpret_cast<decltype(sp)>(*pc);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    if (!shaderPrg){
      return;
    }
    VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
    if (!vbo)
      return;
    vbo->bind(shaderPrg->id, sp->which_attr_idx + sp->npickattrs * I->pick_pass());
  }
}

static void CGO_gl_draw_custom(CCGORenderer * I, CGO_op_data pc){
  const cgo::draw::custom * sp = reinterpret_cast<decltype(sp)>(*pc);

  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg){
    return;
  }
  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vbo)
    return;
  IndexBuffer * ibo = NULL;
  if (sp->iboid){
    ibo = I->G->ShaderMgr->getGPUBuffer<IndexBuffer>(sp->iboid);
  }
  vbo->bind(shaderPrg->id);
  if (ibo){
    ibo->bind();
    glDrawElements(sp->mode, sp->nindices, GL_C_INT_ENUM, 0);
  } else {
    glDrawArrays(sp->mode, 0, sp->nverts);
  }
  vbo->unbind();
  if (sp->pickvboid) {
    VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
    if (pickvbo)
      pickvbo->unbind();
  }
  if (ibo)
    ibo->unbind();

}

static void CGO_gl_draw_sphere_buffers(CCGORenderer * I, CGO_op_data pc) {
  const cgo::draw::sphere_buffers * sp = reinterpret_cast<decltype(sp)>(*pc);
  int num_spheres = sp->num_spheres;
  int attr_color;
  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);
  CShaderPrg *shaderPrg;
  int pickable = 0;

  shaderPrg = I->G->ShaderMgr->Get_DefaultSphereShader(I->info ? I->info->pass : RenderPass::Antialias);
  if (!shaderPrg){
    return;
  }

  attr_color = shaderPrg->GetAttribLocation("a_Color");

  if (I->isPicking){
    vbo->maskAttributes({ attr_color });
    pickable = SettingGet_i(I->G, I->set1, I->set2, cSetting_pickable);
    shaderPrg->Set1i("lighting_enabled", 0);
    if (pickable){
      pickvbo->bind(shaderPrg->id, I->pick_pass());
    } else {
      assert(I->info->pick);
      unsigned char nopick[4] = {};
      I->info->pick->colorNoPick(nopick);
      glVertexAttrib4ubv(attr_color, nopick);
    }
  }

  vbo->bind(shaderPrg->id);
  glDrawArrays(GL_QUADS, 0, num_spheres * 4);

  vbo->unbind();
}

static void CGO_gl_draw_cylinder_buffers(CCGORenderer * I, CGO_op_data pc) {
  const cgo::draw::cylinder_buffers * sp = reinterpret_cast<decltype(sp)>(*pc);
  int  num_cyl = sp->num_cyl;
  int min_alpha = sp->alpha;
  int attr_colors, attr_colors2;
  CShaderPrg *shaderPrg;
  int pickable = 0;
  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  IndexBuffer * ibo = I->G->ShaderMgr->getGPUBuffer<IndexBuffer>(sp->iboid);
  VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);

  shaderPrg = I->G->ShaderMgr->Get_CylinderShader(I->info ? I->info->pass : RenderPass::Antialias);

  if (!shaderPrg){
    return;
  }
  attr_colors = shaderPrg->GetAttribLocation("a_Color");
  attr_colors2 = shaderPrg->GetAttribLocation("a_Color2");

  if (I->isPicking){
    pickable = SettingGet_i(I->G, I->set1, I->set2, cSetting_pickable);
    shaderPrg->Set1i("lighting_enabled", 0);
  }
  if (I->isPicking){
    vbo->maskAttributes({ attr_colors, attr_colors2 });
    if (pickable){
      // in first pass: 1st half of vbo, in second pass: 2nd half of vbo
      // first color (offset 0)
      pickvbo->bind(shaderPrg->id, I->pick_pass());
      // second color (offset 4)
      pickvbo->bind(shaderPrg->id, I->pick_pass() + SHADER_PICKING_PASSES_MAX);
    } else {
      assert(I->info->pick);
      unsigned char nopick[4] = {};
      I->info->pick->colorNoPick(nopick);
      glVertexAttrib4ubv(attr_colors, nopick);
      glVertexAttrib4ubv(attr_colors2, nopick);
    }
  }

  vbo->bind(shaderPrg->id);
  ibo->bind();

  if (min_alpha < 255) {
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glDrawElements(GL_TRIANGLES, num_cyl * NUM_TOTAL_VERTICES_PER_CYLINDER, GL_C_INT_ENUM, 0);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthFunc(GL_LEQUAL);
  }
  glDrawElements(GL_TRIANGLES, num_cyl * NUM_TOTAL_VERTICES_PER_CYLINDER, GL_C_INT_ENUM, 0);

  if (min_alpha < 255) {
    glDepthFunc(GL_LESS);
  }

  ibo->unbind();
  vbo->unbind();
  if (I->isPicking)
    pickvbo->unbind();
}
#include "Texture.h"

static void CGO_gl_draw_labels(CCGORenderer * I, CGO_op_data pc) {
  const cgo::draw::labels * sp = reinterpret_cast<decltype(sp)>(*pc);

  CShaderPrg * shaderPrg;
  int t_mode = SettingGetGlobal_i(I->G, cSetting_transparency_mode);

  if (t_mode==3 && I->info && I->info->pass != RenderPass::Transparent){
    // in transparency_mode=3, labels are drawn in the transparency pass=-1
    return;
  }
  shaderPrg = I->G->ShaderMgr->Get_LabelShader(I->info ? I->info->pass : RenderPass::Antialias);
  if (I->rep){
    float label_size;
    CSetting * set1 = NULL, * set2 = NULL;
    if (I->rep->cs) set1 = I->rep->cs->Setting.get();
    if (I->rep->obj) set2 = I->rep->obj->Setting.get();
    label_size = SettingGet_f(I->G, set1, set2, cSetting_label_size);
    shaderPrg->Set1f("scaleByVertexScale", label_size < 0.f ? 1.f : 0.f);
    if (label_size<0.f){
      shaderPrg->Set1f("labelTextureSize", (float)-2.f* I->info->texture_font_size/label_size);
    }
  }

  if (!shaderPrg){
    return;
  }

  VertexBuffer * vbo     = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  VertexBuffer * pickvbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->pickvboid);

  if (I->isPicking){
    pickvbo->bind(shaderPrg->id, I->pick_pass());
  }

  if (!vbo)
    return;
  vbo->bind(shaderPrg->id);

  glDrawArrays(GL_TRIANGLES, 0, sp->ntextures*6);

  vbo->unbind();
  pickvbo->unbind();
}

static void CGO_gl_draw_connectors(CCGORenderer * I, CGO_op_data pc) {
  int use_geometry_shaders = SettingGetGlobal_b(I->G, cSetting_use_geometry_shaders);

  const cgo::draw::connectors * sp = reinterpret_cast<decltype(sp)>(*pc);

  GLenum mode = GL_LINES;
  int factor = 2;
  float lineWidth;
  if (I->isPicking){
    return;
  }
  {
    GLenum err ;
    CHECK_GL_ERROR_OK("ERROR: CGO_gl_draw_connectors begin returns err=%d\n");
  }

  if (use_geometry_shaders){
    mode = GL_POINTS;
    factor = 1;
  } else {
    factor = 4;
  }
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg){
    return;
  }
  if (I->rep){
    float label_size;
    CSetting * set1 = NULL, * set2 = NULL;
    float v_scale = SceneGetScreenVertexScale(I->G, NULL);
    if (I->rep->cs) set1 = I->rep->cs->Setting.get();
    if (I->rep->obj) set2 = I->rep->obj->Setting.get();
    label_size = SettingGet_f(I->G, set1, set2, cSetting_label_size);
    shaderPrg->Set1f("scaleByVertexScale", label_size < 0.f ? 1.f : 0.f);
    lineWidth = SettingGet_f(I->G, set1, set2, cSetting_label_connector_width);
    if (label_size<0.f){
      shaderPrg->Set1f("textureToLabelSize", v_scale * (float)I->info->texture_font_size/label_size);
    } else {
      shaderPrg->Set1f("textureToLabelSize", 1.f);
    }
  } else {
    lineWidth = SettingGetGlobal_f(I->G, cSetting_label_connector_width);
  }
#ifndef _WEBGL
  if (!use_geometry_shaders)
    glLineWidth(lineWidth);
#endif

  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vbo)
    return;
  vbo->bind(shaderPrg->id);
  glDrawArrays(mode, 0, sp->nconnectors*factor);
  vbo->unbind();
  {
    GLenum err ;
    CHECK_GL_ERROR_OK("ERROR: CGO_gl_draw_connectors end returns err=%d\n");
  }
}

static void CGO_gl_draw_textures(CCGORenderer * I, CGO_op_data pc) {
  const cgo::draw::textures * sp = reinterpret_cast<decltype(sp)>(*pc);
  int ntextures = sp->ntextures;
  VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  CShaderPrg * shaderPrg;
  int attr_pickcolor = 0;
  shaderPrg = I->G->ShaderMgr->Get_LabelShader(I->info ? I->info->pass : RenderPass::Antialias);
  if (!shaderPrg){
    return;
  }
  if (I->isPicking){
    attr_pickcolor = shaderPrg->GetAttribLocation("attr_pickcolor");    
  }
  if (attr_pickcolor){
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glEnableVertexAttribArray(attr_pickcolor);
    glVertexAttribPointer(attr_pickcolor, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, sp->floatdata);
  }
  vbo->bind(shaderPrg->id);
  glDrawArrays(GL_TRIANGLES, 0, ntextures*6);
  vbo->unbind();
  if (attr_pickcolor){
    glDisableVertexAttribArray(attr_pickcolor);
  }
}

static void CGO_gl_draw_screen_textures_and_polygons(CCGORenderer * I, CGO_op_data pc) {
  const cgo::draw::screen_textures * sp = reinterpret_cast<decltype(sp)>(*pc);
  int nverts = sp->nverts;
  CShaderPrg * shaderPrg;

  shaderPrg = I->G->ShaderMgr->Get_ScreenShader();
  if (!shaderPrg){
    return;
  }

  VertexBuffer * vb = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(sp->vboid);
  if (!vb)
    return;
  vb->bind(shaderPrg->id);

  glDrawArrays(GL_TRIANGLES, 0, nverts);

  vb->unbind();
}

static void CGO_gl_draw_trilines(CCGORenderer * I, CGO_op_data pc) {
  int nverts = CGO_get_int(*pc);
  int buffer = CGO_get_int(*pc+1);
  int a_vertex, a_othervertex, a_uv, a_color, a_color2;
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg){
    return;
  }
  a_vertex = 0; // a_Vertex is bound to 0 (see ShaderMgr) CShaderPrg_GetAttribLocation(shaderPrg, "a_Vertex");
  a_othervertex = shaderPrg->GetAttribLocation("a_OtherVertex");
  a_uv = shaderPrg->GetAttribLocation("a_UV");
  a_color = shaderPrg->GetAttribLocation("a_Color");
  a_color2 = shaderPrg->GetAttribLocation("a_Color2");

  glEnableVertexAttribArray(a_vertex);
  glEnableVertexAttribArray(a_othervertex);
  glEnableVertexAttribArray(a_uv);
  glEnableVertexAttribArray(a_color);
  glEnableVertexAttribArray(a_color2);

  glBindBuffer(GL_ARRAY_BUFFER, buffer);

  glVertexAttribPointer(a_vertex, 3, GL_FLOAT, GL_FALSE, 32, (const void *)0);
  glVertexAttribPointer(a_othervertex, 3, GL_FLOAT, GL_FALSE, 32, (const void *)12);
  glVertexAttribPointer(a_uv, 1, GL_FLOAT, GL_FALSE, 32, (const void *)24);
  glVertexAttribPointer(a_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, 32, (const void *)28);
  glVertexAttribPointer(a_color2, 4, GL_UNSIGNED_BYTE, GL_TRUE, 32, (const void *)28);
  glDrawArrays(GL_TRIANGLES, 0, nverts);

  glDisableVertexAttribArray(a_vertex);
  glDisableVertexAttribArray(a_othervertex);
  glDisableVertexAttribArray(a_uv);
  glDisableVertexAttribArray(a_color);
  glDisableVertexAttribArray(a_color2);
}

/* CGO_gl_uniform3f - this is the implementation for the 
 * CGOUniform3f/CGO_UNIFORM3F operation. From the uniform_id
 * it looks up the uniform location from the current shader,
 * and sets it to the values in this op.
 *
 */
static void CGO_gl_uniform3f(CCGORenderer * I, CGO_op_data pc) {
  int uniform_id = CGO_get_int(*pc);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (!shaderPrg){
    return;
  }
  int loc = shaderPrg->GetUniformLocation(
      shaderPrg->uniformLocations[uniform_id].c_str());
  const float *pcp = *pc + 1;
  glUniform3f(loc, pcp[0], pcp[1], pcp[2]);
}

static void CGO_gl_linewidth(CCGORenderer * I, CGO_op_data pc)
{
#ifndef _WEBGL
  glLineWidth(**pc);
#endif
}

/**
 * call glLineWidth and set the "line_width" uniform
 */
static void glLineWidthAndUniform(float line_width,
    CShaderPrg * shaderPrg=NULL) {
#ifndef _WEBGL
  glLineWidth(line_width);
#endif

  if (shaderPrg && shaderPrg->name == "trilines")
    shaderPrg->Set1f("line_width", line_width);
}

/* CGO_gl_special - this is the implementation function for 
   CGOSpecial/CGO_SPECIAL.  Each op has its own implementation.
 */
static void CGO_gl_special(CCGORenderer * I, CGO_op_data pc)
{
  int mode = CGO_get_int(*pc);
  bool openVR = SceneGetStereo(I->G) == cStereo_openvr;
  char varwidth = 0;
  float vScale = (I->info ? I->info->vertex_scale : SceneGetScreenVertexScale(I->G, NULL));

  CSetting *csSetting = NULL, *objSetting = NULL;
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (I->rep && I->rep->cs){
    csSetting = I->rep->cs->Setting.get();
  }
  if (I->rep && I->rep->obj){
    objSetting = I->rep->obj->Setting.get();
  }
  switch (mode){
  case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON:
    {
      float line_width = SettingGet_f(I->G, NULL, NULL, cSetting_ribbon_width);
      if (!openVR) line_width = SceneGetDynamicLineWidth(I->info, line_width);
      if (I->info && I->info->width_scale_flag){
        line_width *= I->info->width_scale;
      }
      glLineWidthAndUniform(line_width, shaderPrg);
    }
    break;
  case LINEWIDTH_DYNAMIC_WITH_SCALE_DASH:
    {
      float line_width = SettingGet_f(I->G, NULL, NULL, cSetting_dash_width);
      if (!openVR) line_width = SceneGetDynamicLineWidth(I->info, line_width);
      if (I->info && I->info->width_scale_flag){
        line_width *= I->info->width_scale;
      }
      glLineWidthAndUniform(line_width, shaderPrg);
    }
    break;
  case LINEWIDTH_DYNAMIC_WITH_SCALE:
    {
      float line_width = SettingGet_f(I->G, NULL, NULL, cSetting_line_width);
      if (!openVR) line_width = SceneGetDynamicLineWidth(I->info, line_width);
      if (I->info && I->info->width_scale_flag){
        line_width *= I->info->width_scale;
      }
      glLineWidthAndUniform(line_width, shaderPrg);
    }
    break;
  case LINEWIDTH_WITH_SCALE:
    {
      float line_width = SettingGet_f(I->G, NULL, NULL, cSetting_line_width);
      if (I->info && I->info->width_scale_flag){
        line_width *= I->info->width_scale;
      }
      glLineWidthAndUniform(line_width, shaderPrg);
    }
    break;
  case LINEWIDTH_DYNAMIC_MESH:
    {
      float line_width;
      if (I->rep){
        line_width = SettingGet_f(I->G, I->rep->cs->Setting.get(), I->rep->obj->Setting.get(), cSetting_mesh_width);
      } else {
        line_width = SettingGet_f(I->G, NULL, NULL, cSetting_mesh_width);
      }
      if (!openVR) line_width = SceneGetDynamicLineWidth(I->info, line_width);
      glLineWidthAndUniform(line_width, shaderPrg);
    }
    break;
  case POINTSIZE_DYNAMIC_DOT_WIDTH:
    {
      float ps;
      if(I->info && I->info->width_scale_flag){
        ps = SettingGet_f
          (I->G, csSetting, objSetting,
           cSetting_dot_width) * I->info->width_scale;
      }
      else {
        ps = SettingGet_f
          (I->G, csSetting, objSetting, cSetting_dot_width);
      }
      glPointSize(ps);
      break;
    }
  case CYLINDERWIDTH_DYNAMIC_MESH:
    {
      CSetting *setting = NULL;
      float mesh_width;
      if (I && I->rep && I->rep->obj){
        setting = I->rep->obj->Setting.get();
      }
      mesh_width = SettingGet_f(I->G, setting, NULL, cSetting_mesh_width);
      if (shaderPrg) {
        const float * color = I->color ? I->color : g_ones4f;
	shaderPrg->Set1f("uni_radius", SceneGetLineWidthForCylinders(I->G, I->info, mesh_width));
        shaderPrg->SetAttrib4fLocation("a_Color", color[0], color[1], color[2], I->alpha);
        shaderPrg->SetAttrib4fLocation("a_Color2", color[0], color[1], color[2], I->alpha);
      }
    }
    break;
  case DOTSIZE_WITH_SPHERESCALE:
    {
      float radius = SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
      radius *= vScale;
      if (shaderPrg)
	shaderPrg->Set1f("sphere_size_scale", fabs(radius));
    }
    break;
  case MESH_WIDTH_FOR_SURFACES:
    {
      float mesh_width = SettingGet_f(I->G, csSetting, objSetting, cSetting_mesh_width);
      if (shaderPrg)
	shaderPrg->Set1f("uni_radius", SceneGetLineWidthForCylinders(I->G, I->info, mesh_width));
    }
    break;
  case CYLINDER_WIDTH_FOR_DISTANCES:
    {
      float line_width, radius;
      int round_ends;
      round_ends =
        SettingGet_b(I->G, csSetting, objSetting, cSetting_dash_round_ends);
      line_width = 
        SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_width);
      radius =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_radius);
      
      line_width = SceneGetDynamicLineWidth(I->info, line_width);

      if (shaderPrg) {
	if(radius == 0.0F) {
          float dash_size = SettingGet_f(I->G, csSetting, objSetting, cSetting_dash_width);
	  shaderPrg->Set1f("uni_radius", SceneGetLineWidthForCylindersStatic(I->G, I->info, line_width, dash_size));
	} else {
	  shaderPrg->Set1f("uni_radius", radius);
	}
	if (!round_ends){
	  shaderPrg->Set1i("no_flat_caps", 0);
	}
      }
    }
    break;
  case CYLINDER_WIDTH_FOR_RIBBONS:
    {
      float pixel_scale_value = SettingGetGlobal_f(I->G, cSetting_ray_pixel_scale);
      float line_width, radius;
      line_width = 
        SettingGet_f(I->G, csSetting, objSetting, cSetting_ribbon_width);
      radius =
        SettingGet_f(I->G, csSetting, objSetting, cSetting_ribbon_radius);
      
      line_width = SceneGetDynamicLineWidth(I->info, line_width);
      if(pixel_scale_value < 0)
        pixel_scale_value = 1.0F;
      if (shaderPrg) {
	if(radius == 0.0F) {
	  shaderPrg->Set1f("uni_radius", vScale * pixel_scale_value * line_width/ 2.f);
	} else {
	  shaderPrg->Set1f("uni_radius", radius);
	}
      }
    }
    break;
  case DOT_WIDTH_FOR_DOTS:
    {
      float dot_width = SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
      float radius;
      if(I->info && I->info->width_scale_flag)
        radius = (dot_width * I->info->width_scale);
      else
        radius= dot_width;
      if (shaderPrg)
	shaderPrg->Set1f("g_PointSize", radius);
#ifndef _PYMOL_IOS
      glPointSize(radius);
#endif
    }
    break;
  case DOT_WIDTH_FOR_DOT_SPHERES:
    {
      float dotSize = SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_radius);
      float dot_width = SettingGet_f(I->G, csSetting, objSetting, cSetting_dot_width);
      float radius;
      if(I->info && dotSize <= 0.0F) {
        if(I->info->width_scale_flag)
          radius = dot_width * I->info->width_scale * I->info->vertex_scale / 1.4142F;
        else
          radius = dot_width * I->info->vertex_scale;
      } else {
        radius = dotSize;
      }
      if (shaderPrg)
	shaderPrg->Set1f("sphere_size_scale", fabs(radius));
    }
    break;
  case CYLINDER_WIDTH_FOR_NONBONDED:
    {
      if (shaderPrg){
        float line_width = SettingGet_f(I->G, csSetting, objSetting, cSetting_line_width);
        shaderPrg->Set1f("uni_radius", SceneGetLineWidthForCylindersStatic(I->G, I->info, line_width, line_width));
      }
    }
    break;
  case CYLINDER_WIDTH_FOR_REPWIRE_VARWIDTH:
    varwidth = 1;
  case CYLINDER_WIDTH_FOR_REPWIRE:
    {
      float radius = SettingGet_f(I->G, csSetting, objSetting, cSetting_line_radius);
      if (radius < R_SMALL8) {
        float line_width = SettingGet_f(I->G, csSetting, objSetting, cSetting_line_width);
        float pixel_scale_value = SettingGetGlobal_f(I->G, cSetting_ray_pixel_scale);
        float vertex_scale = vScale;
        float scale_bound = SettingGetGlobal_f(I->G, cSetting_field_of_view)  * cPI / 180.0f * 0.018f;
        if (!varwidth){
          line_width = SceneGetDynamicLineWidth(I->info, line_width);
        }
        if (vertex_scale < scale_bound) {
          vertex_scale = scale_bound;
        }
        if(pixel_scale_value < 0)
          pixel_scale_value = 1.0F;
        radius = vertex_scale * pixel_scale_value * line_width / 2.f;
      }
      if (shaderPrg){
	shaderPrg->Set1f("uni_radius", radius);
      }
    }
    break;
  case ENABLE_BACK_FACES_IF_NOT_TWO_SIDED:
    {
      int two_sided_lighting = SettingGet_i(I->G, csSetting, objSetting, cSetting_two_sided_lighting) > 0;
      if (!two_sided_lighting){
        glCullFace(GL_BACK);
        glEnable(GL_CULL_FACE);
      }
    }
    break;
  case DISABLE_BACK_FACES_IF_NOT_TWO_SIDED:
    {
      int two_sided_lighting = SettingGet_i(I->G, csSetting, objSetting, cSetting_two_sided_lighting) > 0;
      if (!two_sided_lighting){
        glDisable(GL_CULL_FACE);
      }
    }
    break;
  case SET_SURFACE_UNIFORMS:
    {
      float ambient_occlusion_scale = 0.f;
      int ambient_occlusion_mode = SettingGet_i(I->G, csSetting, objSetting, cSetting_ambient_occlusion_mode);
      
      if (ambient_occlusion_mode){
        ambient_occlusion_scale = SettingGet_f(I->G, csSetting, objSetting, cSetting_ambient_occlusion_scale);
      }
      if (shaderPrg)
	shaderPrg->Set1f("ambient_occlusion_scale", ambient_occlusion_scale);
    }
    break;
  case SET_ALIGNMENT_UNIFORMS_ATTRIBS:
    {
      float linewidth = SettingGet_f(I->G, csSetting, objSetting, cSetting_cgo_line_width);
      float lineradius = SettingGet_f(I->G, csSetting, objSetting, cSetting_cgo_line_radius);
      float pixel_scale_value = SettingGetGlobal_f(I->G, cSetting_ray_pixel_scale);
      if (linewidth < 0.f){
        linewidth = 1.f;
      }
      if(pixel_scale_value < 0)
        pixel_scale_value = 1.0F;
      if (lineradius < 0.f){
        lineradius = linewidth * vScale * pixel_scale_value / 2.f;
      }
      shaderPrg->Set1f("uni_radius", lineradius);
      if (I->color){
        shaderPrg->SetAttrib4fLocation("a_Color", I->color[0], I->color[1], I->color[2], 1.f);
        shaderPrg->SetAttrib4fLocation("a_Color2", I->color[0], I->color[1], I->color[2], 1.f);
      }
      glLineWidthAndUniform(lineradius*2.f / vScale, shaderPrg);
    }
    break;
  case LINEWIDTH_FOR_LINES:
    {
      float line_width = SceneGetDynamicLineWidth(I->info,
          SettingGet_f(I->G, NULL, NULL, cSetting_line_width));
      if (I->info && I->info->width_scale_flag){
        line_width *= I->info->width_scale;
      }
      glLineWidthAndUniform(line_width, shaderPrg);
    }
    break;
  case SET_LABEL_SCALE_UNIFORMS:
  {
    if (I->rep){
      float label_size;
      CSetting * set1 = NULL, * set2 = NULL;
      if (I->rep->cs) set1 = I->rep->cs->Setting.get();
      if (I->rep->obj) set2 = I->rep->obj->Setting.get();
      label_size = SettingGet_f(I->G, set1, set2, cSetting_label_size);
      shaderPrg->Set1f("scaleByVertexScale", label_size < 0.f ? 1.f : 0.f);
      if (label_size<0.f){
        shaderPrg->Set1f("labelTextureSize", (float)-2.f* I->info->texture_font_size/label_size);
      }
    }


  }
  break;
  default:
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_special(): bad mode=%d\n", mode ENDFB(I->G);
  }
}

/* CGO_gl_special_with_arg - this is the implementation function for 
   CGOSpecialWithArg/CGO_SPECIAL_WITH_ARG.  Each op has its own implementation.
 */
static void CGO_gl_special_with_arg(CCGORenderer * I, CGO_op_data pc)
{
#ifndef PURE_OPENGL_ES_2
  int mode = CGO_get_int(*pc);
  float argval = *((*pc) + 1);
  bool use_shaders = SettingGetGlobal_b(I->G, cSetting_use_shaders);
  bool sphere_use_shaders = use_shaders && SettingGetGlobal_b(I->G, cSetting_use_shaders);
  switch(mode){
  case LINEWIDTH_FOR_LINES:
    {
      if (!use_shaders){
        glEnd();
        glLineWidth(argval);
        glBegin(GL_LINES);
      }
    }
    break;
  case LINE_LIGHTING:
    if (!I->isPicking && !SettingGetGlobal_b(I->G, cSetting_use_shaders)) {
      if (!I->info->line_lighting){
        bool enableLighting = (int)argval;
        if (enableLighting)
          glEnable(GL_LIGHTING);
        else
          glDisable(GL_LIGHTING);
      }
    }
    break;
  case SPHERE_MODE_OPS:
    {
      float pixel_scale = 1.0F / I->info->vertex_scale;
      int sphere_mode = (int)fabs(argval);
      bool enable = argval > 0.f;
      if (enable){
        float pointSize;
        if((sphere_mode == 1) || (sphere_mode == 6)) {
          pointSize = SettingGet_f(I->G, I->set1, I->set2, cSetting_sphere_point_size);
          glDisable(GL_POINT_SMOOTH);
          glDisable(GL_ALPHA_TEST);
          if (!I->isPicking && !sphere_use_shaders){
            glEnable(GL_LIGHTING);
            glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
          }
        } else {
          float sphere_scale = SettingGet_f(I->G, I->set1, I->set2, cSetting_sphere_scale);
          if((sphere_mode == 3) || (sphere_mode == 8)) {
            glEnable(GL_POINT_SMOOTH);
            glAlphaFunc(GL_GREATER, 0.5F);
            glEnable(GL_ALPHA_TEST);
            glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
            pointSize = sphere_scale * pixel_scale * 2.0F;
          } else {
            glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
            glDisable(GL_POINT_SMOOTH);
            glDisable(GL_ALPHA_TEST);
            pointSize = sphere_scale * pixel_scale * 1.4F;
          }
        }
        if(!I->isPicking && ((sphere_mode == 7) || (sphere_mode == 8)))
          glEnable(GL_LIGHTING);
        glPointSize(pointSize);
      } else {
        if(sphere_mode == 3) {
          glDisable(GL_POINT_SMOOTH);
          glAlphaFunc(GL_GREATER, 0.05F);
        } else {
          glEnable(GL_ALPHA_TEST);
        }
      }
    }
  }
#endif
}

static void CGO_gl_dotwidth(CCGORenderer * I, CGO_op_data pc)
{
  glPointSize(**pc);
}

static void CGO_gl_enable(CCGORenderer * I, CGO_op_data pc)
{
  GLenum mode = CGO_get_int(*pc);
  CShaderMgr *shaderMgr = I->G->ShaderMgr;
  CShaderPrg *shaderPrg = shaderMgr->Get_Current_Shader();
  if (I->use_shader){
    if (true){
      switch(mode){
      case CGO_GL_LIGHTING:
        {
          if (shaderPrg){
            shaderPrg->SetLightingEnabled(1);
          }
        }
        break;
      case GL_SHADER_LIGHTING:
        if (!I->isPicking){
          if (shaderPrg){
            shaderPrg->SetLightingEnabled(1);
          }
        }
        break;
      case GL_TWO_SIDED_LIGHTING:
        {
          if (shaderPrg){
            shaderPrg->Set1i("two_sided_lighting_enabled", 1);
          }
        }
        break;
      case GL_MESH_LIGHTING:
        {
          int lighting =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_mesh_lighting);
          if (shaderPrg){
            shaderPrg->SetLightingEnabled(lighting);
          }
        }
        break;
      case GL_DOT_LIGHTING:
        {
          int lighting =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_dot_lighting);
          if (shaderPrg && !I->isPicking){
            shaderPrg->SetLightingEnabled(lighting);
            shaderPrg->Set1i("two_sided_lighting_enabled", 0);
          }
        }
        break;
      case GL_LABEL_FLOAT_TEXT:
        {
          int float_text =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
          if (float_text){
            glDisable(GL_DEPTH_TEST);
          }
        }
        break;
      case GL_DASH_TRANSPARENCY_DEPTH_TEST:
        {
          float dash_transparency =
            SettingGet_f(I->G, I->set1, I->set2, cSetting_dash_transparency);
          short dash_transparency_enabled;
          bool t_mode_3 =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_transparency_mode) == 3;
          dash_transparency = (dash_transparency < 0.f ? 0.f : (dash_transparency > 1.f ? 1.f : dash_transparency));
          dash_transparency_enabled = (dash_transparency > 0.f);
          if (dash_transparency_enabled && !t_mode_3 && !I->isPicking){
            glDisable(GL_DEPTH_TEST);
          }
        }
        break;
      case GL_DEFAULT_SHADER:
        shaderMgr->Enable_DefaultShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_LINE_SHADER:
        shaderMgr->Enable_LineShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_SURFACE_SHADER:
        shaderMgr->Enable_SurfaceShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_CYLINDER_SHADER:
        shaderMgr->Enable_CylinderShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_SPHERE_SHADER:
        shaderMgr->Enable_DefaultSphereShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_RAMP_SHADER:
        shaderMgr->Enable_RampShader();
        break;
      case GL_DEFAULT_SHADER_WITH_SETTINGS:
        shaderMgr->Enable_DefaultShaderWithSettings(I->set1, I->set2, I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_BACKGROUND_SHADER:
        shaderMgr->Enable_BackgroundShader();
        break;
      case GL_LABEL_SHADER:
        shaderMgr->Enable_LabelShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_CONNECTOR_SHADER:
        shaderMgr->Enable_ConnectorShader(I->info ? I->info->pass : RenderPass::Antialias);
        break;
      case GL_SCREEN_SHADER:
        shaderMgr->Enable_ScreenShader();
        break;
      case GL_TRILINES_SHADER:
        shaderMgr->Enable_TriLinesShader();
        break;
#ifndef _PYMOL_NO_AA_SHADERS
#endif
      case GL_OIT_SHADER:
        shaderMgr->Enable_OITShader();
        break;
      case GL_OIT_COPY_SHADER:
        shaderMgr->Enable_OITCopyShader();
        break;
      case GL_BACK_FACE_CULLING:
        glCullFace(GL_BACK);
        glEnable(GL_CULL_FACE);
        break;
      case GL_DEPTH_TEST:
        glEnable(mode);
        break;
      case GL_DEPTH_TEST_IF_FLOATING:
        {
          int float_text = SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
          if(float_text)
            glEnable(GL_DEPTH_TEST);
        }
        break;
      }
    }
  } else {
    if (!I->isPicking){
      if (mode==CGO_GL_LIGHTING){
        glEnable(GL_LIGHTING);
      }
    }
  }
}

static void CGO_gl_disable(CCGORenderer * I, CGO_op_data pc)
{
  GLenum mode = CGO_get_int(*pc);
  auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
  if (I->use_shader){
      switch(mode){
      case GL_SHADER_LIGHTING:
        {
          if (shaderPrg){
            shaderPrg->SetLightingEnabled(0);
          }
        }
        break;
      case GL_CYLINDER_SHADER:
        glDisable(GL_CULL_FACE);
      case GL_RAMP_SHADER:
      case GL_SCREEN_SHADER:
      case GL_LABEL_SHADER:
      case GL_CONNECTOR_SHADER:
      case GL_DEFAULT_SHADER:
      case GL_SURFACE_SHADER:
      case GL_SPHERE_SHADER:
      case GL_TRILINES_SHADER:
      case GL_OIT_COPY_SHADER:
      case GL_LINE_SHADER:
        I->G->ShaderMgr->Disable_Current_Shader();
        break;
      case GL_LABEL_FLOAT_TEXT:
        {
          int float_text =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
          if (float_text){
            glEnable(GL_DEPTH_TEST);
          }
        }
        break;
      case GL_DASH_TRANSPARENCY_DEPTH_TEST:
        {
          float dash_transparency =
            SettingGet_f(I->G, I->set1, I->set2, cSetting_dash_transparency);
          short dash_transparency_enabled;
          bool t_mode_3 =
            SettingGet_i(I->G, I->set1, I->set2, cSetting_transparency_mode) == 3;
          dash_transparency = (dash_transparency < 0.f ? 0.f : (dash_transparency > 1.f ? 1.f : dash_transparency));
          dash_transparency_enabled = (dash_transparency > 0.f);
          if (dash_transparency_enabled && !t_mode_3 && !I->isPicking){
            glEnable(GL_DEPTH_TEST);
          }
        }
        break;
      case CGO_GL_LIGHTING:
        {
          if (shaderPrg){
            shaderPrg->SetLightingEnabled(0);
          }
        }
        break;
      case GL_TWO_SIDED_LIGHTING:
        {
          if (shaderPrg){
            shaderPrg->Set1i("two_sided_lighting_enabled", 0);
          }
        }
        break;
#if !defined(PURE_OPENGL_ES_2) || defined(_WEBGL)
      case GL_OIT_SHADER:
      case GL_SMAA1_SHADER:
      case GL_SMAA2_SHADER:
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, I->G->ShaderMgr->default_framebuffer_id);
        break;
#endif
      case GL_BACK_FACE_CULLING:
        glDisable(GL_CULL_FACE);
        break;
      case GL_DEPTH_TEST:
        glDisable(mode);
        break;
      case GL_DEPTH_TEST_IF_FLOATING:
        {
          int float_text = SettingGet_i(I->G, I->set1, I->set2, cSetting_float_labels);
          if(float_text)
            glDisable(GL_DEPTH_TEST);
        }
        break;
      }
  } else {
    if (mode!=CGO_GL_LIGHTING || !I->isPicking){
        if (mode==CGO_GL_LIGHTING)
            mode = GL_LIGHTING;
      glDisable(mode);
    }
  }
}

static void CGO_gl_alpha(CCGORenderer * I, CGO_op_data pc)
{
  I->alpha = **pc;
}

static void CGO_gl_reset_normal(CCGORenderer * I, CGO_op_data pc)
{
  SceneResetNormalUseShader(I->G, CGO_get_int(*pc), I->use_shader);
}

static void CGO_gl_null(CCGORenderer * I, CGO_op_data pc)
{
}

static void CGO_gl_error(CCGORenderer * I, CGO_op_data pc)
{
  PRINTFB(I->G, FB_CGO, FB_Warnings)
  " CGO_gl_error() is not suppose to be called op=%d\n",
      CGO_get_int((*pc) - 1) ENDFB(I->G);
}

static void CGO_gl_color(CCGORenderer* I, CGO_op_data varg)
{
  auto* v = *varg;
  if (I->use_shader){
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    if (shaderPrg){
      int attr_a_Color = shaderPrg->GetAttribLocation("a_Color");
      glVertexAttrib4f(attr_a_Color, v[0], v[1], v[2], I->alpha);
    }
  } else {
    glColor4f(v[0], v[1], v[2], I->alpha);
  }
}

static void CGO_gl_sphere(CCGORenderer * I, CGO_op_data varg)
{
  auto *v = *varg;
  if (I->isPicking){
    SphereRender(I->G, 0, v, I->color, I->alpha, v[3]);
  } else {
    SphereRender(I->G, I->sphere_quality, v, NULL, I->alpha, v[3]);
  }
}

static void CGO_gl_vertex_attribute_3f(CCGORenderer * I, CGO_op_data varg)
{
    auto vertex_attr = reinterpret_cast<const cgo::draw::vertex_attribute_3f *>(*varg);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    int loc = shaderPrg->GetAttribLocation(I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx));
    if (loc >= 0)
      glVertexAttrib3fv(loc, vertex_attr->values);
}

static void CGO_gl_vertex_attribute_4ub(CCGORenderer * I, CGO_op_data varg)
{
    auto vertex_attr = reinterpret_cast<const cgo::draw::vertex_attribute_4ub *>(*varg);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    int loc = shaderPrg->GetAttribLocation(I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx));
    if (loc >= 0)
      glVertexAttrib4ubv(loc, vertex_attr->ubdata);
}

static void CGO_gl_vertex_attribute_4ub_if_picking(CCGORenderer * I, CGO_op_data varg)
{
  if (I->isPicking){
    auto vertex_attr = reinterpret_cast<const cgo::draw::vertex_attribute_4ub_if_picking *>(*varg);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    int loc = shaderPrg->GetAttribLocation(I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx));
    if (loc >= 0)
      glVertexAttrib4ubv(loc, vertex_attr->ubdata);
  }
}

static void CGO_gl_vertex_attribute_1f(CCGORenderer * I, CGO_op_data varg)
{
    auto vertex_attr = reinterpret_cast<const cgo::draw::vertex_attribute_1f *>(*varg);
    auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
    const char *name = I->G->ShaderMgr->GetAttributeName(vertex_attr->attr_lookup_idx);
    int loc = shaderPrg->GetAttribLocation(name);
    if (loc >= 0)
      glVertexAttrib1f(loc, vertex_attr->value);
}

/* dispatch table for OpenGL */

CGO_op_fn CGO_gl[] = {
  CGO_gl_null,                  /* 0x00 */
  CGO_gl_null,                  /* 0x01 */
  CGO_gl_begin,                 /* 0x02 */
  CGO_gl_end,                   /* 0x03 */
  CGO_gl_vertex,                /* 0x04 */
  CGO_gl_normal,                /* 0x05 */
  CGO_gl_color,                 /* 0x06 */
  CGO_gl_sphere,                /* 0x07 */
  CGO_gl_null,                  /* 0x08 */
  CGO_gl_null,                  /* 0x09 */

  CGO_gl_linewidth,             /* 0x0A */
  CGO_gl_null,                  /* 0x0B */
  CGO_gl_enable,                /* 0x0C */
  CGO_gl_disable,               /* 0x0D */
  CGO_gl_null,                  /* 0x0E */
  CGO_gl_null,                  /* 0x0F */

  CGO_gl_dotwidth,              /* 0X10 */
  CGO_gl_null,                  /* 0x11 */
  CGO_gl_null,                  /* 0x12 */
  CGO_gl_null,                  /* 0X13 */

  CGO_gl_null,                  /* 0X14 */
  CGO_gl_null,                  /* 0x15 */
  CGO_gl_null,                  /* 0x16 */
  CGO_gl_null,                  /* 0X17 */

  CGO_gl_null,                  /* 0X18 */
  CGO_gl_alpha,                 /* 0x19 */
  CGO_gl_null,                  /* 0x1A */
  CGO_gl_null,                  /* 0X1B */
  CGO_gl_draw_arrays,           /* 0x1C DrawArrays() */
  CGO_gl_null,                  /* 0x1D */
  CGO_gl_reset_normal,          /* 0x1E */
  CGO_gl_null,                  /* pick color  0X1F */
  CGO_gl_null,                  /* 0x20 draw buffers REMOVED */
  CGO_gl_draw_buffers_indexed,          /* 0x21 draw buffers indexed */
  CGO_gl_null,                  /* 0x22 bounding box */
  CGO_gl_draw_buffers_not_indexed,          /* 0x23 draw buffers not indexed */
  CGO_gl_special,                /* 0x24 special */
  CGO_gl_draw_cylinder_buffers,  /* 0x25 draw GLSL cylinders */
  CGO_gl_null,                  /* 0x26 shader cylinder */
  CGO_gl_null,                  /* 0x27 shader cylinder with 2nd color */
  CGO_gl_draw_sphere_buffers,   /* 0x28 draw sphere buffers */
  CGO_gl_null,                  /* 0x29 accessibility used for ambient occlusion */
  CGO_gl_error,          /* 0x2A draw texture */
  CGO_gl_draw_textures,          /* 0x2B draw textures */
  CGO_gl_draw_screen_textures_and_polygons,          /* 0x2C draw screen textures and polygons */
  CGO_gl_error,
  CGO_gl_error,  CGO_gl_draw_labels,
  CGO_gl_error,  CGO_gl_draw_connectors,  CGO_gl_draw_trilines,  CGO_gl_uniform3f,  CGO_gl_special_with_arg,
  CGO_gl_line,  CGO_gl_splitline,  CGO_gl_draw_custom,
  CGO_gl_vertex_attribute_3f, CGO_gl_vertex_attribute_4ub,
  CGO_gl_vertex_attribute_1f, 
  CGO_gl_mask_attribute_if_picking, CGO_gl_bind_vbo_for_picking,
  CGO_gl_vertex, 
  CGO_gl_null, // interpolated
  CGO_gl_vertex_cross, // CGO_VERTEX_CROSS
  CGO_gl_vertex_attribute_4ub_if_picking,
  CGO_gl_error
};

#if 0
static
void SetUCColorToPrev(uchar *color){
  color[0] = color[-4];
  color[1] = color[-3];
  color[2] = color[-2];
  color[3] = color[-1];
}

static
void SetUCColorToPrev8(uchar *color){
  color[0] = color[-8];
  color[1] = color[-7];
  color[2] = color[-6];
  color[3] = color[-5];
}
#endif

static
void SetUCColorToPrevN(int n, uchar *color){
  color[0] = color[-n*4];
  color[1] = color[-n*4+1];
  color[2] = color[-n*4+2];
  color[3] = color[-n*4+3];
}

static
int * get_pickcolorsset_ptr(int op, float * pc) {
#define RETURN_PICKCOLORSETPTR_CASE(cls) \
  case cgo::draw::cls::op_code: \
    return &(reinterpret_cast<cgo::draw::cls*>(pc)->pickcolorsset)
  switch (op) {
    RETURN_PICKCOLORSETPTR_CASE(buffers_indexed);
    RETURN_PICKCOLORSETPTR_CASE(buffers_not_indexed);
    RETURN_PICKCOLORSETPTR_CASE(labels);
    RETURN_PICKCOLORSETPTR_CASE(sphere_buffers);
    RETURN_PICKCOLORSETPTR_CASE(cylinder_buffers);
    RETURN_PICKCOLORSETPTR_CASE(custom);
  }
  return NULL;
}

void CGORenderGLPicking(CGO * I, RenderInfo *info, PickContext * context, CSetting * set1,
                        CSetting * set2, Rep *rep)
{
  PyMOLGlobals *G = I->G;

  if (!G->ValidContext)
    return;

  if (!I->c)
    return;

  CCGORenderer *R = G->CGORenderer;
  bool pickable = (!I->no_pick) &&
    SettingGet_b(G, set1, set2, cSetting_pickable);
  auto pick = info->pick;
  bool reset_colors = !pick->pickColorsValid();

  R->use_shader = I->use_shader;
  R->isPicking = true;
  R->set1 = set1;
  R->set2 = set2;
  R->info = info;
  R->rep = rep;

#ifndef _WEBGL
      glLineWidth(SettingGet_f(G, set1, set2, cSetting_cgo_line_width));
#endif

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    CGO_OP_DATA_CONST float* pc = it.data();

    switch (op) {
      case CGO_COLOR:
        continue;

      case CGO_PICK_COLOR:

        if (reset_colors){ // only if picking info is invalid
          unsigned char col[4];
          AssignNewPickColor(I, pick, col, context, CGO_get_uint(pc),
              pickable ? CGO_get_int(pc + 1) : cPickableNoPick);
#ifndef PURE_OPENGL_ES_2
          if (!I->use_shader){
            glColor4ubv(col);
          }
#endif
        } else {
          PRINTFB(G, FB_CGO, FB_Warnings)
          " %s: unexpected CGO_PICK_COLOR with !reset_colors\n",
              __func__ ENDFB(G);
        }
        continue;

      case CGO_DRAW_ARRAYS:
        {
          const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
          int arrays = sp->arraybits;
          if (reset_colors && arrays & CGO_PICK_COLOR_ARRAY){ // only if picking info is invalid
            int nverts = sp->nverts, v, idx = -1, bnd = -1;
            float *pca = sp->floatdata;

            if (arrays & CGO_VERTEX_ARRAY) { pca += nverts * VERTEX_POS_SIZE; }
            if (arrays & CGO_NORMAL_ARRAY) { pca += nverts * VERTEX_NORMAL_SIZE; }
            if (arrays & CGO_COLOR_ARRAY) { pca += nverts * VERTEX_COLOR_SIZE; }

            auto pickColorValsUC = (uchar*)pca;
            auto pickColorVals = (int*)(pca + nverts * VERTEX_PICKCOLOR_RGBA_SIZE);

            for (v=0;v<nverts; v++) {
              bnd = pickable ? pickColorVals[v * 2 + 1] : cPickableNoPick;
              idx = pickColorVals[v * 2];
              AssignNewPickColor(
                  I, pick, pickColorValsUC + (v * 4), context, idx, bnd);
            }
          }
        }
        break;

      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
      case CGO_DRAW_TEXTURES:
      case CGO_DRAW_LABELS:
      case CGO_DRAW_SPHERE_BUFFERS:
      case CGO_DRAW_CYLINDER_BUFFERS:
      case CGO_DRAW_CUSTOM:
        {
          int pickcolors_are_set = true;
          int* pickcolors_are_set_ptr = get_pickcolorsset_ptr(op, const_cast<float*>(pc));
          if (!pickcolors_are_set_ptr)
            pickcolors_are_set_ptr = &pickcolors_are_set;

          // TODO remove `pickcolorsset` fields from CGOs
          // This assert can fail during "Roving Detail" demo. However, I still
          // question the need of the `pickcolorsset` fields.
          // assert(reset_colors || *pickcolors_are_set_ptr);

          if (reset_colors || !*pickcolors_are_set_ptr){ // only if picking info is invalid
            int nverts = 0;
            int nvertsperfrag = 1;
            int v, pl;
            int bnd = cPickableNoPick;
            unsigned int idx = 0;
            int srcp;
            float *pca = nullptr;
            int *pickDataSrc ;
            uchar *pickColorDestUC = NULL;
            bool free_pick_color_dest = false;
            int destOffset = 0, bufsizemult = 1;
            size_t pickvbo = 0;
            switch (op){
              case CGO_DRAW_CUSTOM:
              {
                const cgo::draw::custom * sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->nverts;
                pickvbo = sp->pickvboid;
                if (!pickvbo)
                  continue;
                pca = sp->floatdata;
                nvertsperfrag = sp->vertsperpickinfo;
                bufsizemult = sp->npickbufs;

                pickColorDestUC = new uchar[bufsizemult * nverts * 4];
              }
                break;
              case CGO_DRAW_BUFFERS_INDEXED:
              {
                const cgo::draw::buffers_indexed * sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->nverts;
                pickvbo = sp->pickvboid;
                pca = sp->floatdata;
              }
                break;
              case CGO_DRAW_BUFFERS_NOT_INDEXED:
              {
                const cgo::draw::buffers_not_indexed * sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->nverts;
                pickvbo = sp->pickvboid;
                pca = sp->floatdata;
              }
                break;
              case CGO_DRAW_SPHERE_BUFFERS:
              {
                const cgo::draw::sphere_buffers * sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->num_spheres * VERTICES_PER_SPHERE;
                nvertsperfrag = VERTICES_PER_SPHERE;
                pickvbo = sp->pickvboid;
                pca = sp->floatdata;

                pickColorDestUC = new uchar[nverts * 4];
              }
              break;
              case CGO_DRAW_CYLINDER_BUFFERS:
              {
                const cgo::draw::cylinder_buffers * sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->num_cyl * NUM_VERTICES_PER_CYLINDER;
                nvertsperfrag = NUM_VERTICES_PER_CYLINDER;
                pickvbo = sp->pickvboid;
                pca = sp->floatdata;
                bufsizemult = 2;

                pickColorDestUC = new uchar[bufsizemult * nverts * 4];
              }
              break;
              case CGO_DRAW_TEXTURES:
              {
                const cgo::draw::textures * sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->ntextures * 6;
                pca = sp->floatdata;
              }
              break;
              case CGO_DRAW_LABELS:
              {
                const cgo::draw::labels * sp;
                sp = reinterpret_cast<decltype(sp)>(pc);
                nverts = sp->ntextures * 6;
                pca = sp->floatdata;
                pickvbo = sp->pickvboid;
              }
                break;
              }

            if (pickColorDestUC) {
              free_pick_color_dest = true;
              pickDataSrc = (int*)(pca);
            } else {
              pickColorDestUC = (uchar*)pca;
              pickDataSrc = (int*)(pca + nverts);
            }

            destOffset = R->pick_pass() * sizeof(float) * nverts * bufsizemult;

            if (!pickable){
              for (int i = 0; i < nverts * bufsizemult; ++i) {
                pick->colorNoPick(pickColorDestUC + 4 * i);
              }
            } else {
              int npickbufs = bufsizemult;
              int ploffsetforbuf = 0;
              if (op == CGO_DRAW_CYLINDER_BUFFERS){
                  // disabled 2016-07-19 TH: code looks almost identical to
                  // else branch and CGO_DRAW_CYLINDER_BUFFERS seem to be
                  // not used anymore.
                  PRINTFB(I->G, FB_CGO, FB_Errors)
                    " FIXME: SUPPOSEDLY UNUSED CODE EXECUTED in CGORenderGLPicking!\n"
                    ENDFB(I->G);
              } else {
                if (op == CGO_DRAW_CUSTOM){
                  ploffsetforbuf = sizeof(float) * nverts; // for multiple picking attributes
                }
                for (v=0, pl = 0;v<nverts; v++, pl += 4){
                  if (v % nvertsperfrag){
                    // if same fragment, same color
                    for (int pi = 0; pi < npickbufs; ++pi){
                      int ploffset = ploffsetforbuf * pi;
                      SetUCColorToPrevN(1, &pickColorDestUC[pl+ploffset]);
                    }
                    continue;
                  }

                  int frag = (int)(v / nvertsperfrag);
                  for (int pi = 0; pi < npickbufs; ++pi){
                    int ploffset = ploffsetforbuf * pi;
                    srcp = 2* ((npickbufs * frag) + pi);
                    idx = pickDataSrc[srcp];
                    bnd = pickDataSrc[srcp + 1];

                    AssignNewPickColor(I, pick, &pickColorDestUC[pl + ploffset],
                        context, idx, bnd);
                  }
                }
              }
            }

            if (pickvbo) {
              // reload entire vbo
              VertexBuffer * vbo = I->G->ShaderMgr->getGPUBuffer<VertexBuffer>(pickvbo);
              vbo->bufferReplaceData(destOffset, sizeof(float) * nverts * bufsizemult, pickColorDestUC);
              (*pickcolors_are_set_ptr) = true;
            }

            if (free_pick_color_dest){
              delete[] pickColorDestUC;
              pickColorDestUC = NULL;
              free_pick_color_dest = false;
            }
          }
        }
        break;
    }

    CGO_gl[op] (R, &pc);
  }

  R->isPicking = false;
}

void CGORenderGL(CGO * I, const float *color, CSetting * set1, CSetting * set2,
                 RenderInfo * info, Rep *rep)
/* this should be as fast as you can make it...

 * the ASM loop is about 2X long as raw looped GL calls,

 * but hopefully superscaler processors won't care */
{
  PyMOLGlobals *G = I->G;

  const float zee[] = {0.f, 0.f, 1.f};
  const float color_tmp[] = {1.f, 1.f, 1.f};

  if (I->render_alpha){
    // for now, the render_alpha_only flag calls CGOSetZVector/CGORenderGLAlpha
    float *ModMatrix = SceneGetModMatrix(G);
    CGOSetZVector(I, ModMatrix[2], ModMatrix[6], ModMatrix[10]);
    CGORenderGLAlpha(I, info, 1);
    if (I->render_alpha == 1) // right now, render_alpha 1: renders alpha only, 2: renders both alpha and rest
      return;
  }

  if (!G->ValidContext) {
    return;
  }

  if (!I->c) {
    return;
  }

  {
    CCGORenderer *R = G->CGORenderer;
    R->info = info;
    R->use_shader = I->use_shader;
    R->debug = I->debug;
    R->sphere_quality = I->sphere_quality;
    R->rep = rep;
    R->color = color;
    R->alpha = 1.0F - SettingGet_f(G, set1, set2, cSetting_cgo_transparency);
    R->set1 = set1;
    R->set2 = set2;
    // normals should be initialized to the view vector
    // (changed BB 9/14 from SceneResetNormalUseShader(), to CScene->LinesNormal, which was arbitrary, I believe)
    SceneResetNormalToViewVector(I->G, I->use_shader);  

    if (!color) {
      color = color_tmp;
    }

    {
      auto shaderPrg = I->G->ShaderMgr->Get_Current_Shader();
      if (shaderPrg && I->use_shader) {
          shaderPrg->SetAttrib4fLocation("a_Color", color[0], color[1], color[2], R->alpha);
      }
      else {
          glColor4f(color[0], color[1], color[2], R->alpha);
      }
    }

#ifndef PURE_OPENGL_ES_2
    const float width_scale = (info && info->width_scale_flag) ? info->width_scale : 1.f;
    glLineWidth(SettingGet_f(G, set1, set2, cSetting_cgo_line_width) * width_scale);
    glPointSize(SettingGet_f(G, set1, set2, cSetting_cgo_dot_width) * width_scale);
#endif

    if (!(info && info->alpha_cgo)) {
      // Regular CGO dispatch table rendering
      for (auto it = I->begin(); !it.is_stop(); ++it) {
        const auto op = it.op_code();
        assert(op < CGO_sz_size());
        CGO_OP_DATA_CONST float* const pc = it.data();

        CGO_gl[op](R, &pc);
      }
      return;
    }

    /* we're sorting transparent triangles globally */
    {
      {
        int mode = -1;
        int vc = 0;
        // triangle normals
        const float *n0 = NULL, *n1 = NULL, *n2 = NULL;
        // triangle vertices
        const float *v0 = NULL, *v1 = NULL, *v2 = NULL;
        // triangle colors
        const float *c0 = color, *c1 = NULL, *c2 = NULL;

        for (auto it = I->begin(); !it.is_stop(); ++it) {
          const auto op = it.op_code();
          assert(op < CGO_sz_size());
          CGO_OP_DATA_CONST float* const pc = it.data();

          if((R->alpha != 1.f)) {
            switch (op) {       /* transparency */
            case CGO_BEGIN:
              mode = CGO_get_int(pc);
              CGO_gl_begin(R, &pc);
              vc = 0;
              n0 = zee;
              break;
            case CGO_END:
              CGO_gl_end(R, &pc);
              mode = -1;
              break;
            case CGO_NORMAL:
              switch (mode) {
              case GL_TRIANGLES:
              case GL_TRIANGLE_STRIP:
              case GL_TRIANGLE_FAN:
                n0 = pc;
                break;
              default:
                CGO_gl_normal(R, &pc);
              }
              break;
            case CGO_COLOR:
              c0 = pc;
              CGO_gl_color(R, &pc);
              break;
            case CGO_TRIANGLE:
              CGOAlphaTriangle(info->alpha_cgo,
                               pc, pc + 3, pc + 6, pc + 9, pc + 12, pc + 15, pc + 18,
                               pc + 21, pc + 24, R->alpha, R->alpha, R->alpha, false);
              break;
	    case CGO_DRAW_ARRAYS:
	      {
                const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
		int mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
		float *vertexVals = 0, *nxtVals = 0, *colorVals = 0, *normalVals;
		float *vertexVals_tmp = 0, *colorVals_tmp = 0, *normalVals_tmp = 0;
		int step;
		short nxtn = 3;
		nxtVals = vertexVals = vertexVals_tmp = sp->floatdata;
		if (arrays & CGO_NORMAL_ARRAY){
		  nxtVals = normalVals = normalVals_tmp = vertexVals + (nxtn*nverts);
		}
		if (arrays & CGO_COLOR_ARRAY){
		  nxtVals = colorVals = colorVals_tmp = nxtVals + (nxtn*nverts);
		  nxtn = 4;
		}
		switch (mode){
		case GL_TRIANGLES:
		  {
		    for (step = 0; step < nverts; step += 3){
		      if (colorVals_tmp){
			c0 = colorVals_tmp; c1 = colorVals_tmp+4; c2 = colorVals_tmp+8;
		      } else {
			c1 = c2 = c0;
		      }
		      if (normalVals_tmp){
			n0 = normalVals_tmp; n1 = normalVals_tmp+3; n2 = normalVals_tmp+6; 
		      } else {
			n1 = n2 = n0;
		      }
		      CGOAlphaTriangle(info->alpha_cgo,
				       vertexVals_tmp, vertexVals_tmp+3, vertexVals_tmp+6,
				       n0, n1, n2,
				       c0, c1, c2,
				       R->alpha, R->alpha, R->alpha, false);
		      vertexVals_tmp += 9;
		      if (normalVals_tmp){
			normalVals_tmp += 9;
		      }
		      if (colorVals_tmp){
			colorVals_tmp += 12;
		      }
		    }
		  }
		  break;
		case GL_TRIANGLE_STRIP:
		  {
		    if (colorVals_tmp){
		      c1 = colorVals_tmp; c2 = colorVals_tmp+4;
		      colorVals_tmp += 8;
		    } else {
		      c1 = c2 = c0;
		    }
		    if (normalVals_tmp){
		      n1 = normalVals_tmp; n2 = normalVals_tmp+3;
		      normalVals_tmp+= 6;
		    } else {
		      n1 = n2 = n0;
		    }
		    vertexVals_tmp += 6;
		    for (step = 2; step < nverts; step++){
		      if (colorVals_tmp){
			c0 = c1; c1 = c2; c2 = colorVals_tmp;
		      }
		      if (normalVals_tmp){
			n0 = n1; n1 = n2; n2 = normalVals_tmp;
		      }
		      CGOAlphaTriangle(info->alpha_cgo,
				       vertexVals_tmp-6, vertexVals_tmp-3, vertexVals_tmp,
				       n0, n1, n2,
				       c0, c1, c2,
				       R->alpha, R->alpha, R->alpha, false);
		      vertexVals_tmp += 3;
		      if (normalVals_tmp){
			normalVals_tmp += 3;
		      }
		      if (colorVals_tmp){
			colorVals_tmp += 4;
		      }
		    }
		  }
		  break;
		case GL_TRIANGLE_FAN:
		  {
		    float *firstVertex = vertexVals_tmp;
		    if (colorVals_tmp){
		      c0 = colorVals_tmp;
		      c2 = colorVals_tmp + 4;
		      colorVals_tmp += 8;
		    } else {
		      c1 = c2 = c0;
		    }
		    if (normalVals_tmp){
		      n0 = normalVals_tmp; 
		      n2 = normalVals_tmp + 3;
		      normalVals_tmp += 6;
		    }
		    vertexVals_tmp += 6;
		    for (step = 2; step < nverts; step++){
		      if (colorVals_tmp){
			c1 = c2; c2 = colorVals_tmp;
		      }
		      if (normalVals_tmp){
			n1 = n2; n2 = normalVals_tmp;
		      }
		      CGOAlphaTriangle(info->alpha_cgo,
				       firstVertex, vertexVals_tmp-3, vertexVals_tmp,
				       n0, n1, n2,
				       c0, c1, c2,
				       R->alpha, R->alpha, R->alpha, false);
		      vertexVals_tmp += 3;
		      if (normalVals_tmp){
			normalVals_tmp += 3;
		      }
		      if (colorVals_tmp){
			colorVals_tmp += 4;
		      }
		    }
		  }
		  break;
		}
	      }
	      break;
            case CGO_VERTEX:
              v0 = pc;
              switch (mode) {
              case GL_TRIANGLES:
                if(3 * ((vc + 1) / 3) == vc + 1) {
                  CGOAlphaTriangle(info->alpha_cgo,
                                   v0, v1, v2, n0, n1, n2, c0, c1, c2,
                                   R->alpha, R->alpha, R->alpha, true);
                }
                v2 = v1;
                c2 = c1;
                n2 = n1;
                v1 = v0;
                c1 = c0;
                n1 = n0;
                vc++;
                break;
              case GL_TRIANGLE_STRIP:
                if(vc > 1) {
                  CGOAlphaTriangle(info->alpha_cgo,
                                   v0, v1, v2, n0, n1, n2, c0, c1, c2,
                                   R->alpha, R->alpha, R->alpha, !(vc & 0x1));
                }
                v2 = v1;
                c2 = c1;
                n2 = n1;
                v1 = v0;
                c1 = c0;
                n1 = n0;
                vc++;
                break;
              case GL_TRIANGLE_FAN:
                if(vc > 1) {
                  CGOAlphaTriangle(info->alpha_cgo,
                                   v0, v1, v2, n0, n1, n2, c0, c1, c2,
                                   R->alpha, R->alpha, R->alpha, false);
                } else if(!vc) {
                  n2 = n0;
                  v2 = v0;
                  c2 = c0;
                }
                v1 = v0;
                c1 = c0;
                n1 = n0;
                vc++;
                break;
              default:
                CGO_gl_vertex(R, &pc);
                break;
              }
              break;
            default:
              CGO_gl[op] (R, &pc);
              break;
            }
          } else {              /* opaque */
	    switch(op){
	    case CGO_COLOR:
	      /* Since CGO operations are done in sequence, alpha could happen 
		 after color is set.  In this case, we still need to keep track of the color 
		 in case there is a transparent object */
	      c0 = pc;
	      break;
	    default:
	      break;
	    }
            CGO_gl[op] (R, &pc);
          }
        }
      }
    }
  }
}

void CGORenderGLAlpha(CGO * I, RenderInfo * info, bool calcDepth)
{
  PyMOLGlobals *G = I->G;
  if(G->ValidContext && I->c) {
    int mode = GL_TRIANGLES;
    if (I->debug){
      mode = CGOConvertDebugMode(I->debug, GL_TRIANGLES);
    }
#ifndef PURE_OPENGL_ES_2
    // not sure why shader is set, but disable it for now,
    // since we are doing immediate mode rendering for global transparency
    G->ShaderMgr->Disable_Current_Shader();
#endif
    /* 1. transform and measure range (if not already known) 
       2. bin into linked lists based on Z-centers
       3. render by layer */

    if(I->z_flag) {
      if(!I->i_start) {
        I->i_size = 256;
        I->i_start = pymol::calloc<int>(I->i_size);
      } else {
        UtilZeroMem(I->i_start, sizeof(int) * I->i_size);
      }
      {
        const int i_size = I->i_size;
        const float* base = I->op;
        int *start = I->i_start;
        int delta = 1, ntris = 0;
        /* bin the triangles */
	if (calcDepth){
          for (auto it = I->begin(); !it.is_stop(); ++it) {
            if (it.op_code() == CGO_ALPHA_TRIANGLE) {
              float* const pc = it.data();
              const float z = dot_product3f(pc + 1, I->z_vector);
	      if(z > I->z_max)
		I->z_max = z;
	      if(z < I->z_min)
		I->z_min = z;
	      pc[4] = z;
	      ntris++;
	    }
          }
        }

        const float range_factor = (0.9999F * i_size) / (I->z_max - I->z_min);

        for (auto it = I->begin(); !it.is_stop(); ++it) {
          if (it.op_code() == CGO_ALPHA_TRIANGLE) {
            float* const pc = it.data();
            assert(base < pc && pc < I->op + I->c);
            auto i = pymol::clamp<int>((pc[4] - I->z_min) * range_factor, 0, i_size);
            CGO_put_int(pc, start[i]);
            start[i] = (pc - base);     /* NOTE: will always be > 0 since we have CGO_read_int'd */
          }
        }

        // for single-layer transparency, render front-to-back
        if(SettingGetGlobal_i(G, cSetting_transparency_mode) == 2) {
          delta = -1;
          start += (i_size - 1);
        }

        /* now render by bin */
#ifndef PURE_OPENGL_ES_2
        glBegin(mode);
        for (int i = 0; i < i_size; i++) {
          int ii = *start;
          start += delta;
          while(ii) {
            const float* pc = base + ii;
            glColor4fv(pc + 23);
            glNormal3fv(pc + 14);
            glVertex3fv(pc + 5);
            glColor4fv(pc + 27);
            glNormal3fv(pc + 17);
            glVertex3fv(pc + 8);
            glColor4fv(pc + 31);
            glNormal3fv(pc + 20);
            glVertex3fv(pc + 11);

            ii = CGO_get_int(pc);
          }
        }
        glEnd();
#endif
      }
    } else {
#ifndef PURE_OPENGL_ES_2
      glBegin(mode);
      for (auto it = I->begin(); !it.is_stop(); ++it) {
        if (it.op_code() == CGO_ALPHA_TRIANGLE) {
          float* const pc = it.data();
          glColor4fv(pc + 23);
          glNormal3fv(pc + 14);
          glVertex3fv(pc + 5);
          glColor4fv(pc + 27);
          glNormal3fv(pc + 17);
          glVertex3fv(pc + 8);
          glColor4fv(pc + 31);
          glNormal3fv(pc + 20);
          glVertex3fv(pc + 11);
        }
      }
      glEnd();
#endif
    }
  }
}


/* translation function which turns cylinders and spheres into triangles */

static int CGOSimpleSphere(CGO * I, const float *v, float vdw, short sphere_quality)
{
  SphereRec *sp;
  int *q, *s;
  int b, c;
  int ok = true;
  /* cgo_sphere_quality is between 0 and (NUMBER_OF_SPHERE_LEVELS-1) */

  sp = I->G->Sphere->Sphere[CLAMPVALUE<short>(sphere_quality, 0, (NUMBER_OF_SPHERE_LEVELS-1)) ];

  q = sp->Sequence;
  s = sp->StripLen;

  for(b = 0; b < sp->NStrip; b++) {
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
    for(c = 0; ok && c < (*s); c++) {
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

static int CGOSimpleQuadric(CGO * I, const float *v, float r, const float *q)
{
  float r_el, n0[3], n1[3], n2[3];
  int ok = true;
  if(CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    ok &= CGOSimpleEllipsoid(I, v, r_el, n0, n1, n2);
  return ok;
}

static int CGOSimpleEllipsoid(CGO* I, const float* v, float vdw,
    const float* n0, const float* n1, const float* n2)
{
  SphereRec *sp;
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

  ds = SettingGet_i(I->G, NULL, NULL, cSetting_cgo_ellipsoid_quality);
  if(ds < 0)
    ds = SettingGet_i(I->G, NULL, NULL, cSetting_ellipsoid_quality);
  if(ds < 0)
    ds = 0;
  if(ds > 3)
    ds = 3;
  sp = I->G->Sphere->Sphere[ds];

  q = sp->Sequence;
  s = sp->StripLen;

  for(b = 0; b < sp->NStrip; b++) {
    ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
    for(c = 0; ok && c < (*s); c++) {
      float *sp_dot_q = sp->dot[*q];
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
      
      for(i = 0; i < 3; i++) {
	vv[i] = d0[i] + d1[i] + d2[i];
      }
      normalize23f(vv, direction);
      add3f(v, vv, vv);
      
      dd0 = dot_product3f(direction, nn0);
      dd1 = dot_product3f(direction, nn1);
      dd2 = dot_product3f(direction, nn2);
      
      if(scale[0] > R_SMALL8) {
	ss0 = dd0 / scale_sq[0];
      } else {
	ss0 = 0.0F;
      }
      if(scale[1] > R_SMALL8) {
	ss1 = dd1 / scale_sq[1];
      } else {
	ss1 = 0.0F;
      }
      
      if(scale[2] > R_SMALL8) {
	ss2 = dd2 / scale_sq[2];
      } else {
	ss2 = 0.0F;
      }
      
      scale3f(nn0, ss0, comp0);
      scale3f(nn1, ss1, comp1);
      scale3f(nn2, ss2, comp2);
      
      for(i = 0; i < 3; i++) {
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
void CGORoundNub(CGO * I,
    const float *v1,    // cap center
    const float *p0,    // normal along axis
    const float *p1,    // x coord in cap space
    const float *p2,    // y coord in cap space
    int direction,      // 1 or -1
    int nEdge,          // "quality"
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
  for (int c = 1; c < cmax; c += 1){
    float z1 = z2;
    z2 = cos(c * PI_over_cmax);

    // around cylinder axis (longitudinal)
    for (int d = (nEdge + 1) * (-direction); d; d += direction){
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

static int CGOSimpleCylinder(CGO * I, const float *v1, const float *v2, const float tube_size,
                             const float *c1, const float *c2, const float alpha1,
                             const float alpha2, const bool interp, const cCylCap cap1, const cCylCap cap2,
                             const Pickable *pickcolor2, const bool stick_round_nub)
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
  if (pickcolor2){
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

  if(nEdge > MAX_EDGE)
    nEdge = MAX_EDGE;
  subdivide(nEdge, x, y);

  colorFlag = (c1 != c2) && c2;
  colorFlag |= alpha1 != alpha2;

  interpColorFlag = c2 && interp && pickcolor2;
  if (interpColorFlag){
    average3f(c1, c2, midcolor);
    midalpha = (alpha1 + alpha2) / 2.0f;
  }
  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  normalize3f(p0);

  if(cap1 == cCylCapRound && !stick_round_nub) {
    vv1[0] = v1[0] - p0[0] * overlap;
    vv1[1] = v1[1] - p0[1] * overlap;
    vv1[2] = v1[2] - p0[2] * overlap;
  } else {
    vv1[0] = v1[0];
    vv1[1] = v1[1];
    vv1[2] = v1[2];
  }
  if(cap2 == cCylCapRound && !stick_round_nub) {
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
  if (pickcolor2){
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
  for(c = nEdge; ok && c >= 0; c--) {
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
    if(ok && (colorFlag || interpColorFlag) ){
      ok &= CGOColorv(I, c1);
      ok &= CGOAlpha(I, alpha1);
    }
    if (ok)
      ok &= CGOVertexv(I, v + 3);
    if (ok && interpColorFlag){
      ok &= CGOColorv(I, midcolor);
      ok &= CGOAlpha(I, midalpha);
    } else if(ok && colorFlag && !pickcolor2){
      ok &= CGOColorv(I, c2);
      ok &= CGOAlpha(I, alpha2);
    }
    if (ok)
      ok &= CGOVertexv(I, v + 6);
  }
  if (ok)
    ok &= CGOEnd(I);
  if (pickcolor2){
    ok &= CGOColorv(I, c2);
    ok &= CGOAlpha(I, alpha2);
    CGOPickColor(I, pickcolor2->index, pickcolor2->bond);
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_STRIP);
    for(c = nEdge; ok && c >= 0; c--) {
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
      if (ok && interpColorFlag){
        ok &= CGOColorv(I, midcolor);
        ok &= CGOAlpha(I, midalpha);
      }
      if (ok)
        ok &= CGOVertexv(I, v + 3);
      if (ok && interpColorFlag){
        ok &= CGOColorv(I, c2);
        ok &= CGOAlpha(I, alpha2);
      }
      if (ok)
        ok &= CGOVertexv(I, v + 6);
    }
    if (ok)
      ok &= CGOEnd(I);
  }

  if(ok && cap1 != cCylCap::None) {
    if(ok && colorFlag && c1){
      ok &= CGOColorv(I, c1);
      ok &= CGOAlpha(I, alpha1);
    }
    if (pickcolor2)
      CGOPickColor(I, pickcolor[0].index, pickcolor[0].bond);

    if(stick_round_nub && cap1 == cCylCapRound) {
      CGORoundNub(I, v1, p0, p1, p2, -1, nEdge, tube_size);
    } else {
      v[0] = -p0[0];
      v[1] = -p0[1];
      v[2] = -p0[2];

      if(cap1 == cCylCapRound) {
        v[3] = vv1[0] - p0[0] * nub;
        v[4] = vv1[1] - p0[1] * nub;
        v[5] = vv1[2] - p0[2] * nub;
      } else {
        v[3] = vv1[0];
        v[4] = vv1[1];
        v[5] = vv1[2];
      }

      if (ok)  ok &= CGOBegin(I, GL_TRIANGLE_FAN);
      if (ok)  ok &= CGONormalv(I, v);
      if (ok)  ok &= CGOVertexv(I, v + 3);

      for(c = nEdge; ok && c >= 0; c--) {
        v[0] = p1[0] * x[c] + p2[0] * y[c];
        v[1] = p1[1] * x[c] + p2[1] * y[c];
        v[2] = p1[2] * x[c] + p2[2] * y[c];

        v[3] = vv1[0] + v[0] * tube_size;
        v[4] = vv1[1] + v[1] * tube_size;
        v[5] = vv1[2] + v[2] * tube_size;

        if(cap1 == cCylCapRound)
          ok &= CGONormalv(I, v);
        if (ok)
          ok &= CGOVertexv(I, v + 3);
      }
      if (ok)
        ok &= CGOEnd(I);
    }
  }

  if(ok && cap2 != cCylCap::None) {
    if(ok && colorFlag && c2){
      ok &= CGOColorv(I, c2);
      ok &= CGOAlpha(I, alpha2);
    }
    if (pickcolor2)
      CGOPickColor(I, pickcolor2->index, pickcolor2->bond);

    if(stick_round_nub && cap2 == cCylCapRound) {
      CGORoundNub(I, v2, p0, p1, p2, 1, nEdge, tube_size);
    } else {
      v[0] = p0[0];
      v[1] = p0[1];
      v[2] = p0[2];

      if(cap2 == cCylCapRound) {
        v[3] = vv2[0] + p0[0] * nub;
        v[4] = vv2[1] + p0[1] * nub;
        v[5] = vv2[2] + p0[2] * nub;
      } else {
        v[3] = vv2[0];
        v[4] = vv2[1];
        v[5] = vv2[2];
      }

      if (ok) ok &= CGOBegin(I, GL_TRIANGLE_FAN);
      if (ok) ok &= CGONormalv(I, v);
      if (ok) ok &= CGOVertexv(I, v + 3);

      for(c = 0; ok && c <= nEdge; c++) {
        v[0] = p1[0] * x[c] + p2[0] * y[c];
        v[1] = p1[1] * x[c] + p2[1] * y[c];
        v[2] = p1[2] * x[c] + p2[2] * y[c];

        v[3] = vv2[0] + v[0] * tube_size;
        v[4] = vv2[1] + v[1] * tube_size;
        v[5] = vv2[2] + v[2] * tube_size;

        if(cap2 == cCylCapRound)
          ok &= CGONormalv(I, v);
        if (ok)
          ok &= CGOVertexv(I, v + 3);
      }
      if (ok) ok &= CGOEnd(I);
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

  if(nEdge > MAX_EDGE)
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

    if(len) {
      slope = (r1 - r2) / len;
    }
    for(c = nEdge; c >= 0; c--) {
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
  for(c = nEdge; ok && c >= 0; c--) {
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
    if(ok && colorFlag)
      CGOColorv(I, c1);
    if (ok)
      CGOVertexv(I, v + 3);
    if(ok && colorFlag)
      CGOColorv(I, c2);
    if (ok)
      CGOVertexv(I, v + 6);
  }
  if (ok)
    ok &= CGOEnd(I);

  if(ok && cap1 != cCylCap::None) {
    v[0] = -p0[0];
    v[1] = -p0[1];
    v[2] = -p0[2];

    {
      v[3] = vv1[0];
      v[4] = vv1[1];
      v[5] = vv1[2];
    }

    if(colorFlag)
      ok &= CGOColorv(I, c1);
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_FAN);
    if (ok)
      ok &= CGONormalv(I, v);
    if (ok)
      ok &= CGOVertexv(I, v + 3);

    for(c = nEdge; ok && c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv1[0] + v[0] * r1;
      v[4] = vv1[1] + v[1] * r1;
      v[5] = vv1[2] + v[2] * r1;

      if(cap1 == cCylCapRound)
        ok &= CGONormalv(I, v);
      if (ok)
	ok &= CGOVertexv(I, v + 3);
    }
    if (ok)
      ok &= CGOEnd(I);
  }

  if(ok && cap2 != cCylCap::None) {

    v[0] = p0[0];
    v[1] = p0[1];
    v[2] = p0[2];

    {
      v[3] = vv2[0];
      v[4] = vv2[1];
      v[5] = vv2[2];
    }

    if(colorFlag)
      ok &= CGOColorv(I, c2);
    if (ok)
      ok &= CGOBegin(I, GL_TRIANGLE_FAN);
    if (ok)
      ok &= CGONormalv(I, v);
    if (ok)
      ok &= CGOVertexv(I, v + 3);

    for(c = 0; ok && c <= nEdge; c++) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv2[0] + v[0] * r2;
      v[4] = vv2[1] + v[1] * r2;
      v[5] = vv2[2] + v[2] * r2;

      if(cap2 == cCylCapRound)
        ok &= CGONormalv(I, v);
      if (ok)
	ok &= CGOVertexv(I, v + 3);
    }
    if (ok)
      ok &= CGOEnd(I);
  }
  return ok;
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
  has_draw_buffers            |= src->has_draw_buffers;
  has_draw_cylinder_buffers   |= src->has_draw_cylinder_buffers;
  has_draw_sphere_buffers     |= src->has_draw_sphere_buffers;
  has_begin_end               |= src->has_begin_end;
  use_shader                  |= src->use_shader;
  render_alpha                |= src->render_alpha;
  src->has_draw_buffers = false;
}

/**
 * Appends `src` to the end of this CGO and then free's `src`
 * and sets the pointer to NULL.
 */
void CGO::free_append(CGO * &src) {
  free_append(std::move(src));
  assert(src == nullptr);
}
void CGO::free_append(CGO * &&src) {
  if (!src)
    return;
  move_append(std::move(*src));
  DeleteP(src);
}

int CGOAppend(CGO *dest, const CGO *source, bool stopAtEnd){
  return dest->append(*source, stopAtEnd);
}

int CGOCountNumberOfOperationsOfType(const CGO *I, int optype){
  std::set<int> ops = { optype };
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

int CGOCountNumberOfOperationsOfTypeN(const CGO *I, const std::map<int, int> &optype){
  int numops = 0;
  for (auto cgoit = I->begin(); !cgoit.is_stop(); ++cgoit) {
    auto it = optype.find(cgoit.op_code());
    if (it != optype.end())
      numops += it->second;
  }
  return (numops);
}

bool CGOHasOperationsOfType(const CGO *I, int optype){
  std::set<int> ops = { optype };
  return CGOHasOperationsOfTypeN(I, ops);
}

bool CGOHasOperations(const CGO *I) {
  return !I->begin().is_stop();
}

bool CGOHasOperationsOfTypeN(const CGO *I, const std::set<int> &optype){
  if (!I->op)
    return false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (optype.count(it.op_code()))
      return 1;
  }
  return (0);
}

static
bool CGOFilterOutOperationsOfTypeN(const CGO *I, CGO *cgo, const std::set<int> &optype){
  if (!I->op)
    return false;

  bool ret = false;
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto op = it.op_code();
    if (optype.find(op) == optype.end()){
      auto pc = it.data();
      cgo->add_to_cgo(op, pc);
    } else {
      ret = true; // returns if filtered anything
    }
  }
  return ret;
}

bool CGOFilterOutCylinderOperationsInto(const CGO *I, CGO *cgo){
  static std::set<int> optypes = { CGO_SHADER_CYLINDER, 
                                   CGO_SHADER_CYLINDER_WITH_2ND_COLOR,
                                   CGO_SAUSAGE,
                                   CGO_CYLINDER,
                                   CGO_CUSTOM_CYLINDER,
                                   CGO_CUSTOM_CYLINDER_ALPHA };
  return CGOFilterOutOperationsOfTypeN(I, cgo, optypes);
}

bool CGOHasCylinderOperations(const CGO *I){
  static std::set<int> optypes = { CGO_SHADER_CYLINDER, 
                                   CGO_SHADER_CYLINDER_WITH_2ND_COLOR,
                                   CGO_SAUSAGE,
                                   CGO_CYLINDER,
                                   CGO_CUSTOM_CYLINDER,
                                   CGO_CUSTOM_CYLINDER_ALPHA };
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOHasSphereOperations(const CGO *I){
  static std::set<int> optypes = { CGO_SPHERE };
  return CGOHasOperationsOfTypeN(I, optypes);
}

bool CGOCheckWhetherToFree(PyMOLGlobals * G, CGO *I){
  if (I->use_shader){
    if (I->cgo_shader_ub_color != SettingGetGlobal_i(G, cSetting_cgo_shader_ub_color) || 
	I->cgo_shader_ub_normal != SettingGetGlobal_i(G, cSetting_cgo_shader_ub_normal)){
      return true;
    }
  }
  return false;
}

CGO *CGOConvertLinesToShaderCylinders(const CGO * I, int est){

  int tot_nverts = 0, tot_ncyls = 0;

  CGO *cgo = CGONewSized(I->G, I->c + est);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_DRAW_ARRAYS:
      {
        const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
	float *vals = cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
	int nvals = sp->narrays*sp->nverts;
        memcpy(vals, sp->floatdata, nvals);
      }
      break;
    case CGO_END:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOConvertLinesToShaderCylinders: CGO_END encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOConvertLinesToShaderCylinders: CGO_VERTEX encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_BEGIN:
      {
	const float *last_vertex = NULL, *last_color = NULL, *current_color = NULL, *color = NULL ;
        unsigned int last_pick_color_idx = 0;
	int last_pick_color_bnd = cPickableNoPick ;
	int nverts = 0, err = 0;
	int mode = CGO_get_int(pc);

        for (++it; !err && it != CGO_END; ++it) {
          auto pc = it.data();
          int op = it.op_code();

	  switch (op) {
	  case CGO_VERTEX:
	    if (last_vertex){
	      switch (mode){
	      case GL_LINES:
	      case GL_LINE_STRIP:
		{
		  float axis[3];
                  bool pick_color_diff = false;
		  axis[0] = pc[0] - last_vertex[0];
		  axis[1] = pc[1] - last_vertex[1];
		  axis[2] = pc[2] - last_vertex[2];
                  pick_color_diff = (cgo->current_pick_color_index != last_pick_color_idx ||
                                     cgo->current_pick_color_bond != last_pick_color_bnd);
		  if (last_color && current_color &&
                      (!equal3f(last_color, current_color) || pick_color_diff)){
		    CGOColorv(cgo, last_color);
                    if (pick_color_diff){
                      Pickable pickcolor2 = { cgo->current_pick_color_index, cgo->current_pick_color_bond };
                      CGOPickColor(cgo, last_pick_color_idx, last_pick_color_bnd);
                      cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, last_vertex, axis, 1.f, cCylShaderBothCapsRound, current_color, &pickcolor2);
                      CGOPickColor(cgo, pickcolor2.index, pickcolor2.bond);
                    } else {
                      cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, last_vertex, axis, 1.f, cCylShaderBothCapsRound, current_color);
                    }
		    CGOColorv(cgo, current_color);
		  } else {
		    cgo->add<cgo::draw::shadercylinder>(last_vertex, axis, 1.f, cCylShaderBothCapsRound);
		  }
		  last_vertex = pc;
                  last_pick_color_idx = cgo->current_pick_color_index;
                  last_pick_color_bnd = cgo->current_pick_color_bond;
		  tot_ncyls++;
		}
		if (mode==GL_LINES){
		  last_vertex = NULL;
		  last_color = NULL;
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
          case CGO_LINE:
	    {
              float axis[3];
              auto line = reinterpret_cast<const cgo::draw::line *>(pc);
              subtract3f(line->vertex2, line->vertex1, axis);
              cgo->add<cgo::draw::shadercylinder>(line->vertex1, axis, 1.f, cCylShaderBothCapsRound);
              tot_ncyls++;
            }
            break;
          case CGO_SPLITLINE:
	    {
              float axis[3];
              auto splitline = reinterpret_cast<const cgo::draw::splitline *>(pc);
              Pickable pickcolor2 = { splitline->index, splitline->bond };
              float color2[] = { CONVERT_COLOR_VALUE(splitline->color2[0]),
                                 CONVERT_COLOR_VALUE(splitline->color2[1]),
                                 CONVERT_COLOR_VALUE(splitline->color2[2]) };
              unsigned char flags = splitline->flags;
              subtract3f(splitline->vertex2, splitline->vertex1, axis);
              if ((flags & cgo::draw::splitline::equal_colors) &&
                  (flags & cgo::draw::splitline::no_split_for_pick)){
                cgo->add<cgo::draw::shadercylinder>(splitline->vertex1, axis, 1., cCylShaderBothCapsRound);
              } else {
                int cap = cCylShaderBothCapsRound;
                if (flags & splitline->flags & cgo::draw::splitline::interpolation){
                  cap |= cCylShaderInterpColor;
                }
                cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, splitline->vertex1, axis, 1., cap, color2, &pickcolor2);
                last_pick_color_idx = splitline->index;
                last_pick_color_bnd = splitline->bond;
              }
              tot_ncyls++;
            }
            break;
	  case CGO_COLOR:
	    if (op == CGO_COLOR){
	      last_color = current_color;
	      current_color = pc;
	      color = pc;
	    }
          case CGO_PICK_COLOR:
            if (op == CGO_PICK_COLOR){
              cgo->current_pick_color_index = CGO_get_uint(pc);
              cgo->current_pick_color_bond = CGO_get_int(pc + 1);
            }
	  default:
            cgo->add_to_cgo(op, pc);
	  }
	}

	tot_nverts += nverts;
      }
      break;
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  if (tot_ncyls){
    return (cgo);
  } else {
    CGOFree(cgo);
    return NULL;
  }
}

#ifdef WITH_UNUSED_FUNCTIONS
/* CGOSplitUpLinesForPicking: This operation goes through */
/* a CGO and returns a new CGO that has the same lines but */
/* a line that has two different pick colors will get split */
/* at its midpoint into two separate lines so that it can */
/* be used for picking */
CGO *CGOSplitUpLinesForPicking(const CGO * I){
  auto G = I->G;

  std::unique_ptr<CGO> cgo_managed(new CGO(G));
  CGO* cgo = cgo_managed.get();
  int tot_nverts = 0;

  CGOBegin(cgo, GL_LINES);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_PICK_COLOR:
      {
        cgo->current_pick_color_index = CGO_get_uint(pc);
        cgo->current_pick_color_bond = CGO_get_int(pc + 1);
      }
      break;
    case CGO_END:
    case CGO_VERTEX:
      WARN_UNEXPECTED_OPERATION(G, op);
      return nullptr;
    case CGO_BEGIN:
      {
	const float *last_vertex = nullptr, *last_color = nullptr,
              *current_color = nullptr, *color = nullptr;
        unsigned int last_pick_color_idx = 0;
	int last_pick_color_bnd = cPickableNoPick ;
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
	    if (last_vertex){
	      switch (mode){
	      case GL_LINES:
	      case GL_LINE_STRIP:
		{
                  bool pick_color_diff = false;
                  pick_color_diff = (cgo->current_pick_color_index != last_pick_color_idx ||
                                     cgo->current_pick_color_bond != last_pick_color_bnd);
		  if (pick_color_diff || 
                      (last_color && current_color &&
                       (!equal3f(last_color, current_color)))){
                    if (pick_color_diff){
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
		if (mode==GL_LINES){
		  last_vertex = NULL;
		  last_color = NULL;
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
	  case CGO_COLOR:
	    {
	      last_color = current_color;
	      current_color = pc;
	      color = pc;
	    }
            break;
          case CGO_PICK_COLOR:
            {
              cgo->current_pick_color_index = CGO_get_uint(pc);
              cgo->current_pick_color_bond = CGO_get_int(pc + 1);
            }
	    break;
	  }
	}
	tot_nverts += nverts;
      }
      break;
    }
  }

  if (!tot_nverts) {
    return nullptr;
  }

  CGOEnd(cgo);
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }

  return cgo_managed.release();
}
#endif

static void trilinesBufferAddVertex(float * &buffer,
    const float * v1,           // vertex
    const float * v2,           // vertex other end of line
    const float * color,        // RGB color
    float alpha,                // alpha
    signed char uv)              // uv
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
  unsigned char *byte_view = (unsigned char *)(buffer++);
  (*byte_view++) = CLIP_COLOR_VALUE(color[0]);
  (*byte_view++) = CLIP_COLOR_VALUE(color[1]);
  (*byte_view++) = CLIP_COLOR_VALUE(color[2]);
  (*byte_view++) = CLIP_COLOR_VALUE(alpha);
}

static void trilinesBufferAddVertices(float * &buffer,
    const float * v1,           // vertex
    const float * v2,           // vertex other end of line
    const float * color,        // RGB color
    float alpha)                // alpha
{
  // Vertex 1
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 1); //-1, 1);
  // Vertex 3
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 3); //1, 1);
  // Vertex 2
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 0); //-1, -1);
  // Vertex 4
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 3); //1, 1);
  // Vertex 3
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 2); //1, -1);
  // Vertex 1
  trilinesBufferAddVertex(buffer, v1, v2, color, alpha, 1); //-1, 1);
}

static
void CGOTrilines_GetCurrentColor(const float*& current_color,
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
void CGOChangeShadersTo(CGO *I, int frommode, int tomode){
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    if (it.op_code() == CGO_ENABLE) {
      auto eo = it.cast<cgo::draw::enable>();
      if (eo->mode == frommode) {
        eo->mode = tomode;
      }
    }
  }
}

CGO *CGOOptimizeScreenTexturesAndPolygons(CGO * I, int est)
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
    uchar *colorValsUC = 0;
    CGOAlpha(cgo, 1.f);
    cgo->alpha = 1.f;
    cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;

    {
      int mul = 6; // 3 - screenoffset/vertex, 2 - texture coordinates, 1 - color
      /*
      if (SettingGetGlobal_i(I->G, cSetting_cgo_shader_ub_color)){
	mul++;
      } else {
	mul += 4;
	}*/
      tot = num_total_indices * mul ;
    }

    auto vertexVals_managed = std::vector<float>(tot);
    float* vertexVals = vertexVals_managed.data();
    texcoordVals = vertexVals + 3 * num_total_indices;
    nxtn = 2;
    colorVals = texcoordVals + nxtn * num_total_indices;
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
    ok = CGOProcessScreenCGOtoArrays(G, I, vertexVals, texcoordVals, colorVals, colorValsUC);
    RETURN_VAL_IF_FAIL(ok && !G->Interrupt, nullptr);
    if (ok){
      VertexBuffer * vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>();
      ok = vbo->bufferData({
          BufferDesc( "attr_screenoffset", GL_FLOAT,         3, sizeof(float) * num_total_indices * 3, vertexVals,   GL_FALSE ),
          BufferDesc( "attr_texcoords", GL_FLOAT,         2, sizeof(float) * num_total_indices * 2, texcoordVals, GL_FALSE ),
          BufferDesc( "attr_backgroundcolor", GL_UNSIGNED_BYTE, 4, sizeof(uchar) * num_total_indices * 4, colorValsUC,  GL_TRUE )
        });
      size_t vboid = vbo->get_hash_id();
      if (ok){
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
  float white[3] = { 1.f, 1.f, 1.f};
  float probe_radius = SettingGet_f(G, set1, NULL, cSetting_solvent_radius);
  float v_above[3], n0[3] = { 0.f, 0.f, 0.f };
  int ramp_above = SettingGet_i(G, set1, NULL, cSetting_surface_ramp_above_mode) == 1;

  for (auto it = I->begin(); ok && !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    bool skipCopy = false;
    switch (op) {
    case CGO_DRAW_ARRAYS:
      {
        const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
	float *vals = cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
	int nvals = sp->narrays*sp->nverts;
	ok &= vals ? true : false;
	if (ok)
          memcpy(vals, sp->floatdata, nvals);
	skipCopy = true;
      }
      break;
    case CGO_NORMAL:
      copy3f(pc, n0);
      break;
    case CGO_VERTEX:
      {
	float color[3];
	copy3f(white, color);
	if (ramp_above){
	  copy3f(n0, v_above);
	  scale3f(v_above, probe_radius, v_above);
	  add3f(pc, v_above, v_above);
	} else {
	  copy3f(pc, v_above);
	}
	if (ObjectGadgetRampInterVertex(ramp, v_above, color, state)){
	  CGOColorv(cgo, color);
	} else {
	  CGOColorv(cgo, white);
	}
      }
      break;
    }
    if (!skipCopy){
      cgo->add_to_cgo(op, pc);
    }
  }
  if (ok){
    ok &= CGOStop(cgo);
    if (ok){
      cgo->use_shader = I->use_shader;
      if (cgo->use_shader){
	cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
	cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
      }
    }
  }
  if (!ok){
    CGOFree(cgo);
  }
  return (cgo);
}

/**
 * FIXME: This function always returns true for `checkOpaque=ture`
 */
int CGOHasTransparency(const CGO *I, bool checkTransp, bool checkOpaque){
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

/* TransparentInfoSortIX - This function sorts all n_tri triangle 
 * centroids in the array sum by:
 * 1) computing z-value in array z_value
 * 2) bin sorting z_values and placing indices in ix array (using Util.cpp)
 *
 * - uses sort_mem as pre-allocated memory to sort
 * - t_mode - either forward (1) or backwards (0) sort
 */
void TransparentInfoSortIX(PyMOLGlobals * G, 
			   float *sum, float *z_value, int *ix,
			   int n_tri, int *sort_mem, int t_mode){
  float *zv;
  float *sv;
  float matrix[16];
  int idx;
  
#ifdef PURE_OPENGL_ES_2
  copy44f(SceneGetModelViewMatrix(G), matrix);
#else
  glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
#endif
  zv = z_value;
  sv = sum;
  
  /* for each triangle, computes the z */
  for (idx = 0; idx<n_tri; ++idx){
    *(zv++) = matrix[2] * sv[0] + matrix[6] * sv[1] + matrix[10] * sv[2];
    sv += 3;
  }

  UtilZeroMem(sort_mem, sizeof(int) * (n_tri + 256));

  switch (t_mode) {
  case 1:
    UtilSemiSortFloatIndexWithNBinsImpl(sort_mem, n_tri, 256, z_value, ix, true); // front to back
    /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn); */
    break;
  default:
    UtilSemiSortFloatIndexWithNBinsImpl(sort_mem, n_tri, 256, z_value, ix, false); // back to front
    /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
    break;
  }
}

CGO *CGOConvertTrianglesToAlpha(const CGO * I){
  int tot_nverts = 0;

  CGO *cgo = CGONewSized(I->G, I->c);

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();

    switch (op) {
    case CGO_DRAW_ARRAYS:
      {
        const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
        const int mode = sp->mode, arrays = sp->arraybits, nverts = sp->nverts;
        const float* nxtVals = sp->floatdata;

        assert(arrays & CGO_VERTEX_ARRAY);
        const float* vertexValsDA = nxtVals;
        nxtVals += VERTEX_POS_SIZE * nverts;

        const float* normalValsDA = nullptr;
        if (arrays & CGO_NORMAL_ARRAY){
          normalValsDA = nxtVals;
          nxtVals += VERTEX_NORMAL_SIZE * nverts;
        }

        const float* colorValsDA = nullptr;
        if (arrays & CGO_COLOR_ARRAY){
          colorValsDA = nxtVals;
          nxtVals += VERTEX_COLOR_SIZE * nverts;
        }

        const float* const vertexVals0 = vertexValsDA;
        const float* const normalVals0 = normalValsDA;
        const float* const colorVals0 = colorValsDA;

        switch (mode){
        case GL_TRIANGLES:
          {
            for (int cnt = 0; cnt < nverts; cnt+=3){
              if (colorVals0){
                CGOAlphaTriangle(cgo, vertexValsDA, vertexValsDA+3, vertexValsDA+6,
                                 normalValsDA, normalValsDA+3, normalValsDA+6,
                                 colorValsDA, colorValsDA+4, colorValsDA+8,
                                 *(colorValsDA + 3), *(colorValsDA + 7), *(colorValsDA + 11), 0);
              } else {
                CGOAlphaTriangle(cgo, vertexValsDA, vertexValsDA+3, vertexValsDA+6,
                                 normalValsDA, normalValsDA+3, normalValsDA+6,
                                 cgo->color, cgo->color, cgo->color,
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
        case GL_TRIANGLE_STRIP:
          {
            short flip = 0;
            vertexValsDA += 6;
            normalValsDA += 6;
            if (colorVals0)
              colorValsDA += 8;

            for (int cnt = 2; cnt < nverts; cnt++){
              if (colorVals0){
                CGOAlphaTriangle(cgo, vertexValsDA-6, vertexValsDA-3, vertexValsDA,
                                 normalValsDA-6, normalValsDA-3, normalValsDA,
                                 colorValsDA-8, colorValsDA-4, colorValsDA,
                                 *(colorValsDA - 5), *(colorValsDA - 1), *(colorValsDA + 3), flip);
              } else {
                CGOAlphaTriangle(cgo, vertexValsDA-6, vertexValsDA-3, vertexValsDA,
                                 normalValsDA-6, normalValsDA-3, normalValsDA,
                                 cgo->color, cgo->color, cgo->color,
                                 cgo->alpha, cgo->alpha, cgo->alpha, flip);
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
        case GL_TRIANGLE_FAN:
          {
            vertexValsDA += 6;
            normalValsDA += 6;
            if (colorVals0)
              colorValsDA += 8;
            for (int cnt = 2; cnt < nverts; cnt++){
              if (colorVals0){
                CGOAlphaTriangle(cgo, vertexVals0, vertexValsDA-3, vertexValsDA,
                                 normalVals0, normalValsDA-3, normalValsDA,
                                 colorVals0, colorValsDA-4, colorValsDA,
                                 *(colorVals0 + 3), *(colorValsDA - 1), *(colorValsDA + 3), 0);
              } else {
                CGOAlphaTriangle(cgo, vertexVals0, vertexValsDA-3, vertexValsDA,
                                 normalVals0, normalValsDA-3, normalValsDA,
                                 cgo->color, cgo->color, cgo->color,
                                 cgo->alpha, cgo->alpha, cgo->alpha, 0);
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
      }
      break;
    case CGO_END:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOConvertTrianglesToAlpha: CGO_END encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOConvertTrianglesToAlpha: CGO_VERTEX encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_BEGIN:
      {
        float vertices[3][3], colors[4][3], normals[4][3], alpha[4] ;
        short verticespl = 2, colorspl = 2, normalspl = 2, alphapl = 2;
        short hasShifted = 0;
        int nverts = 0, err = 0;
        int mode = CGO_get_int(pc);
        short mode_is_triangles = 0, flip = 0, mode_is_fan = 0;
        copy3f(cgo->color, colors[3]);
        copy3f(cgo->normal, normals[3]);
        alpha[3] = cgo->alpha;
        switch (mode){
        case GL_TRIANGLE_FAN:
          mode_is_fan = 1;
        case GL_TRIANGLES:
        case GL_TRIANGLE_STRIP:
          mode_is_triangles = 1;
        }
        if (!mode_is_triangles){
          CGOBegin(cgo, mode);
        }

        for (++it; !err && it != CGO_END; ++it) {
          auto pc = it.data();
          int op = it.op_code();
          short add_to_cgo = 1;
          switch (op) {
          case CGO_VERTEX:
            if (mode_is_triangles){
              if (!(hasShifted & 1)){ // colors
                if (colorspl>=0){
                  copy3f(colors[colorspl+1], colors[colorspl]);
                  colorspl--;
                } else {
                  if (!mode_is_fan)
                    copy3f(colors[1], colors[2]);
                  copy3f(colors[0], colors[1]);
                }
              }
              if (!(hasShifted & 2)){ // normals
                if (normalspl>=0){
                  copy3f(normals[normalspl+1], normals[normalspl]);
                  normalspl--;
                } else {
                  if (!mode_is_fan)
                    copy3f(normals[1], normals[2]);
                  copy3f(normals[0], normals[1]);
                }
              }
              if (!(hasShifted & 4)){ // alphas
                if (alphapl>=0){
                  alpha[alphapl] = alpha[alphapl+1];
                  alphapl--;
                } else {
                  if (!mode_is_fan)
                    alpha[2] = alpha[1];
                  alpha[1] = alpha[0];
                }
              }
              if (verticespl>=0){
                copy3f(pc, vertices[verticespl]) ;
                verticespl--;
              } else {
                if (!mode_is_fan)
                  copy3f(vertices[1], vertices[2]);
                copy3f(vertices[0], vertices[1]);
                copy3f(pc, vertices[0]);
              }

              nverts++;
              switch (mode){
              case GL_TRIANGLES:
                if (!(nverts % 3)){
                  CGOAlphaTriangle(cgo, vertices[2], vertices[1], vertices[0],
                                   normals[2], normals[1], normals[0],
                                   colors[2], colors[1], colors[0], 
                                   alpha[2], alpha[1], alpha[0], 0);
                }
                break;
              case GL_TRIANGLE_STRIP:
                if (verticespl<0){
                  int off0, off2;
                  if (flip){ off0 = 0; off2 = 2; } else { off0 = 2; off2 = 0; }
                  flip = !flip;
                  CGOAlphaTriangle(cgo, vertices[off0], vertices[1], vertices[off2],
                                   normals[off0], normals[1], normals[off2],
                                   colors[off0], colors[1], colors[off2], 
                                   alpha[off0], alpha[1], alpha[off2], 0);
                }
                break;
              case GL_TRIANGLE_FAN:
                if (verticespl<0){
                  CGOAlphaTriangle(cgo, vertices[2], vertices[1], vertices[0],
                                   normals[2], normals[1], normals[0],
                                   colors[2], colors[1], colors[0], 
                                   alpha[2], alpha[1], alpha[0], 0);
                }
              }
              add_to_cgo = !mode_is_triangles;
              hasShifted = 0;
            } else {
              add_to_cgo = 1;
            }
          case CGO_COLOR:
            if (op == CGO_COLOR){
              add_to_cgo = !mode_is_triangles;
              if (mode_is_triangles){
                if (colorspl>=0){
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
            if (op == CGO_NORMAL){
              add_to_cgo = !mode_is_triangles;
              if (mode_is_triangles){
                if (normalspl>=0){
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
            if (op == CGO_ALPHA){
              add_to_cgo = !mode_is_triangles;
              if (mode_is_triangles){
                if (alphapl>=0)
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
            if (add_to_cgo){
              cgo->add_to_cgo(op, pc);
            }
          }
        }

        if (!mode_is_triangles) {
          CGOEnd(cgo);
        }

        tot_nverts += nverts;
      }
      break;
    case CGO_COLOR:
      if (op==CGO_COLOR){
        copy3f(pc, cgo->color);
      }
    case CGO_NORMAL:
      if (op==CGO_NORMAL){
        copy3f(pc, cgo->normal);
      }
    case CGO_ALPHA:
      if (op==CGO_ALPHA){
        cgo->alpha = *pc;
      }
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  if (tot_nverts){
    return (cgo);
  } else {
    CGOFree(cgo);
    return NULL;
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
CGO *CGOGenerateNormalsForTriangles(const CGO * I){
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
          auto * indices = flip ? indices_flipped : indices_regular;

          if (mode != GL_TRIANGLES) {
            flip = !flip;
          }

          CalculateTriangleNormal(vertices[0],
              vertices[indices[1]],
              vertices[indices[2]], current_normal);

          CGONormalv(cgo, current_normal);

          for (int j = 0; j < 3; ++j) {
            int k = indices[j];
            if (has_color) CGOColorv(cgo, colors[k]);
            if (has_alpha) CGOAlpha(cgo, alphas[k]);
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
          " CGO-Warning: CGOGenerateNormalsForTriangles: unhandled op=0x%02x inside BEGIN/END\n",
          op ENDFB(G);
        cgo->add_to_cgo(op, pc);
    }
  }

  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
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
    case CGO_DRAW_ARRAYS:
      {
        const cgo::draw::arrays * sp = reinterpret_cast<decltype(sp)>(pc);
        float *vals;
        int nvals = sp->narrays*sp->nverts;
        switch (sp->mode){
        case GL_LINES:
        case GL_LINE_STRIP:
          CGODisable(cgo, CGO_GL_LIGHTING);
          cur_mode_is_lines = true;
        }
        vals = cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
        memcpy(vals, sp->floatdata, nvals);
        if (cur_mode_is_lines){
          CGOEnable(cgo, CGO_GL_LIGHTING);
          cur_mode_is_lines = false;
        }
      }
      break;
    case CGO_DRAW_BUFFERS_INDEXED:
      {
        const cgo::draw::buffers_indexed * sp = reinterpret_cast<decltype(sp)>(pc);
        int mode = sp->mode, mode_is_lines = 0;
        switch (mode){
        case GL_LINES:
        case GL_LINE_STRIP:
          mode_is_lines = true;
        }
        if (mode_is_lines){
          CGODisable(cgo, CGO_GL_LIGHTING);
        }
        cgo->copy_op_from<cgo::draw::buffers_indexed>(pc);
        if (mode_is_lines){
          CGOEnable(cgo, CGO_GL_LIGHTING);
        }
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
        const cgo::draw::buffers_not_indexed * sp = reinterpret_cast<decltype(sp)>(pc);
        int mode = sp->mode, mode_is_lines = 0;
        switch (mode){
        case GL_LINES:
        case GL_LINE_STRIP:
          mode_is_lines = true;
        }
        if (mode_is_lines){
          CGODisable(cgo, CGO_GL_LIGHTING);
        }
        cgo->copy_op_from<cgo::draw::buffers_not_indexed>(pc);
        if (mode_is_lines){
          CGOEnable(cgo, CGO_GL_LIGHTING);
        }
      }
      break;
    case CGO_END:
      {
        CGOEnd(cgo);
        if (cur_mode_is_lines){
          CGOEnable(cgo, CGO_GL_LIGHTING);
          cur_mode_is_lines = 0;
        }
      }
      break;
    case CGO_BEGIN:
      {
        int mode = CGO_get_int(pc);
        switch (mode){
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
      }
      break;
    default:
      cgo->add_to_cgo(op, pc);
    }
  }
  cgo->use_shader = use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  return (cgo);
}

bool CGOHasAnyTriangleVerticesWithoutNormals(const CGO *I, bool checkTriangles){
  bool inside = false;
  bool hasNormal = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_BEGIN:
      switch (CGO_get_int(pc)){
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
    case CGO_DRAW_ARRAYS:
      {
        const auto sp = it.cast<cgo::draw::arrays>();
        switch (sp->mode){
        case GL_TRIANGLE_FAN:
        case GL_TRIANGLES:
        case GL_TRIANGLE_STRIP:
          if (checkTriangles){
            if (!(sp->arraybits & CGO_NORMAL_ARRAY)){
              return 1;
            }
          }
          break;
        case GL_LINE_STRIP:
        case GL_LINES:
          if (!checkTriangles){
            if (!(sp->arraybits & CGO_NORMAL_ARRAY)){
              return 1;
            }
          }
          break;
        }
      }
      break;
    }
  }
  return 0;
}

/**
 * CGOReorderIndicesWithTransparentInfo : This function
 * takes the triangle index array ix (result from TransparentInfoSortIX)
 * and sets the vertices (vertexIndices) for each triangle from the original
 * indices (vertexIndicesOriginal), then uses glBufferData to set the 
 * GL_ELEMENT_ARRAY_BUFFER to these indices.
 *
 */
void CGOReorderIndicesWithTransparentInfo(PyMOLGlobals * G, 
					  int nindices, size_t vbuf, 
					  int n_tri, int *ix, 
					  GL_C_INT_TYPE *vertexIndicesOriginal,
					  GL_C_INT_TYPE *vertexIndices){
  int c, pl, idx;
  IndexBuffer * ibo = G->ShaderMgr->getGPUBuffer<IndexBuffer>( vbuf );
  if (!vertexIndices){
    PRINTFB(G, FB_RepSurface, FB_Errors) "ERROR: RepSurfaceRender() vertexIndices is not set, nindices=%d\n", nindices ENDFB(G);
  }
  /* updates the vertexIndices from the ix array */
  for(c = 0, pl=0; c < n_tri; c++) {
    idx = ix[c] * 3;
    vertexIndices[pl++] = vertexIndicesOriginal[idx];
    vertexIndices[pl++] = vertexIndicesOriginal[idx + 1];
    vertexIndices[pl++] = vertexIndicesOriginal[idx + 2];
  }
  ibo->bufferSubData(0, sizeof(GL_C_INT_TYPE) * nindices, vertexIndices);
}

void CGO::add_to_cgo(int op, const float * pc) {
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

void CGO::print_table() const {
}

CGO* CGOConvertSpheresToPoints(const CGO* I)
{
  CGO *cgo;

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
      CGOPickColor(cgo, cgo->current_pick_color_index, cgo->current_pick_color_bond);
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
  if (ok){
    ok &= CGOStop(cgo);
  } 
  if (!ok){
    CGOFree(cgo);
  }
  return (cgo);
}

// Will create an interleaved VBO with { vertex, otherVertex, uv, and texcoord info }
// Currently, this function does not support/parse lines inside CGODrawArrays, i.e.,
// a CGO that CGOCombineBeginEnd was used on
CGO *CGOConvertLinesToTrilines(const CGO * I, bool addshaders){
  static std::set<int> lineops = { CGO_VERTEX, CGO_LINE, CGO_SPLITLINE };
  auto G = I->G;
  const int nLines = CGOCountNumberOfOperationsOfTypeN(I, lineops ) + 1;

  if (nLines == 0) {
    return nullptr;
  }

  int line_counter = 0;
  GLuint glbuff = 0;
  const float *colorv = NULL;
  unsigned int buff_size = nLines * 6 * (8 * sizeof(float));
  
  // VBO memory -- number of lines x 4 vertices per line x (vertex + otherVertex + normal + texCoord + color)
  std::vector<float> buffer_start(buff_size);
  float *buffer = buffer_start.data();

  std::unique_ptr<CGO> cgo(new CGO(G));

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    const auto op = it.op_code();
    const auto pc = it.data();

    switch (op) {
    case CGO_DRAW_ARRAYS:
    {
        auto sp = it.cast<cgo::draw::arrays>();
	float *vals = cgo->add<cgo::draw::arrays>(sp->mode, sp->arraybits, sp->nverts);
	int nvals = sp->narrays*sp->nverts;
        memcpy(vals, sp->floatdata, nvals);
      }
      break;
    case CGO_END:
      WARN_UNEXPECTED_OPERATION(G, op);
      return nullptr;
    case CGO_BEGIN:
      {
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
	    if (last_vertex){
	      switch (mode){
	      case GL_LINES:
	      case GL_LINE_STRIP:
		{
		  float cc[3] = { 1, 1, 1 };
		  float alpha = cgo->alpha;
                  CGOTrilines_GetCurrentColor(current_color, colorv, last_color, cc);
		  trilinesBufferAddVertices(buffer, pc, last_vertex, current_color, alpha);
		  line_counter++;
		  last_vertex = pc;
		}
		if (mode==GL_LINES){
		  last_vertex = NULL;
		  last_color = NULL;
		}
	      }
	    } else {
	      last_vertex = pc;
	      current_color = color;
	    }
            break;
	  case CGO_LINE:
            {
              auto line = it.cast<cgo::draw::line>();
              float cc[3] = { 1, 1, 1 };
              float alpha = cgo->alpha;
              CGOTrilines_GetCurrentColor(current_color, colorv, last_color, cc);
              trilinesBufferAddVertices(buffer, line->vertex1, line->vertex2, current_color, alpha);
              line_counter++;
            }
            break;
	  case CGO_SPLITLINE:
            {
              auto splitline = it.cast<cgo::draw::splitline>();
              float cc[3] = { 1, 1, 1 };
              float alpha = cgo->alpha;
              float mid[3];
              float color2[] = { CONVERT_COLOR_VALUE(splitline->color2[0]),
                                 CONVERT_COLOR_VALUE(splitline->color2[1]),
                                 CONVERT_COLOR_VALUE(splitline->color2[2]) };
              add3f(splitline->vertex1, splitline->vertex2, mid);
              mult3f(mid, .5f, mid);
              CGOTrilines_GetCurrentColor(current_color, colorv, last_color, cc);
              trilinesBufferAddVertices(buffer, splitline->vertex1, mid, current_color, alpha);
              trilinesBufferAddVertices(buffer, mid, splitline->vertex2, color2, alpha);
              line_counter+=2;
            }
            break;
	  case CGO_COLOR:
	      last_color = current_color;
	      current_color = pc;
	      color = pc;
            break;
	  }
	}
      }
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_COLOR:
      colorv = pc;
      break;
    }
  }

  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGetGlobal_i(cgo->G, cSetting_cgo_shader_ub_normal);
  }

  {
    int err = 0;
    glGenBuffers(1, &glbuff);
    glBindBuffer(GL_ARRAY_BUFFER, glbuff);
    glBufferData(GL_ARRAY_BUFFER, line_counter * 6 * 8 * sizeof(float), buffer_start.data(), GL_STATIC_DRAW);
    CHECK_GL_ERROR_OK("ERROR: CGOConvertLinesToTriangleStrips() glBindBuffer returns err=%d\n");
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
 * copies data for a particular attribute operation into the array used to load the VBO.
 * this takes into account whether it is interleaved or not.
 * 
 * isInterleaved    : whether the VBO is interleaved
 * nvert            : which vertex in the VBO
 * attribOp         : the attribute op
 * vertexDataSize   : total vertex data size in VBO (for interleaved)
 * dataPtrs         : all data pointers for attributes (for interleaved, they are all the pointer to the one array)
 * attrOffset       : offsets of the attributes (for interleaved)
 * pcarg            : pc pointer to CGO operation data
 * pick_data        : pointer to pick data for current vertex (writes to if pick data)
 * has_pick_colorBS : keeps track of which pick attributes have been set
 *
 */
static
void copyAttributeForOp(bool isInterleaved, int &nvert, AttribOp *attribOp, int vertexDataSize, std::vector<void*> &dataPtrs,
                        std::vector<int> &attrOffset, const float *pcarg, float *pick_data, int &has_pick_colorBS, int pstride){
  auto attrDesc = attribOp->desc;
  int ord = attrDesc->order;
  int copyord = -1;
  void *dataPtr = dataPtrs[ord];
  unsigned char *pc = ((unsigned char *)pcarg) + attribOp->offset;
  if (isInterleaved){
    dataPtr = (unsigned char*) dataPtr + nvert * vertexDataSize + attrOffset[ord];
    if (attribOp->copyAttribDesc){
      copyord = attribOp->copyAttribDesc->order;
      pc = ((unsigned char*) dataPtrs[ord]) + nvert * vertexDataSize + attrOffset[copyord];
    }
  } else {
    int sz = gl_sizeof(attrDesc->type_size) * attrDesc->type_dim;
    dataPtr = (unsigned char*) dataPtr + nvert * sz;
    if (attribOp->copyAttribDesc){
      copyord = attribOp->copyAttribDesc->order;
      int copysz = gl_sizeof(attribOp->copyAttribDesc->type_size) * attribOp->copyAttribDesc->type_dim;
      pc = (unsigned char*) dataPtr + nvert * copysz;
    }
  }
  switch (attribOp->conv_type){
  case NO_COPY:
    break;
  case FLOAT_TO_FLOAT:
    *((float *) dataPtr) = *((float *)pc);
    break;
  case FLOAT2_TO_FLOAT2:
    *((float *) dataPtr) = *((float *)pc);
    *((float *) dataPtr + 1) = *((float *)pc + 1);
    break;
  case FLOAT3_TO_FLOAT3:
    copy3f((float*)pc, (float*)dataPtr);
    break;
  case FLOAT4_TO_FLOAT4:
    *((float *) dataPtr) = *((float *)pc);
    *((float *) dataPtr + 1) = *((float *)pc + 1);
    *((float *) dataPtr + 2) = *((float *)pc + 2);
    *((float *) dataPtr + 3) = *((float *)pc + 3);
    break;
  case FLOAT3_TO_UB3:
    {
      auto dataPtrUB = (unsigned char *)dataPtr;
      float *pcf = (float*)pc;
      dataPtrUB[0] = CLIP_COLOR_VALUE(pcf[0]);
      dataPtrUB[1] = CLIP_COLOR_VALUE(pcf[1]);
      dataPtrUB[2] = CLIP_COLOR_VALUE(pcf[2]);
    }
    break;
  case FLOAT1_TO_UB_4TH:
    {
      auto dataPtrUB = (unsigned char *)dataPtr;
      float *pcf = (float*)pc;
      dataPtrUB[3] = CLIP_COLOR_VALUE(pcf[0]);
    }
    break;
  case UB3_TO_UB3:
    {
      auto dataPtrUB = (unsigned char *)dataPtr;
      auto pcUB = (unsigned char *)pc;
      dataPtrUB[0] = pcUB[0];
      dataPtrUB[1] = pcUB[1];
      dataPtrUB[2] = pcUB[2];
      break;
    }
  case UINT_INT_TO_PICK_DATA:
    {
      float *pcf = (float*)pc;
      unsigned int index = CGO_get_uint(pcf);
      int bond = CGO_get_int(pcf+1);
      CGO_put_uint(ord * 2 + pick_data, index);
      CGO_put_int(ord * 2 + pick_data + 1, bond);
      has_pick_colorBS |= (1 << ord) ;
    }
    break;
  case FLOAT4_TO_UB4:
    {
      auto dataPtrUB = (unsigned char *)dataPtr;
      float *pcf = (float*)pc;
      dataPtrUB[0] = CLIP_COLOR_VALUE(pcf[0]);
      dataPtrUB[1] = CLIP_COLOR_VALUE(pcf[1]);
      dataPtrUB[2] = CLIP_COLOR_VALUE(pcf[2]);
      dataPtrUB[3] = CLIP_COLOR_VALUE(pcf[3]);
    }
    break;
  case CYL_CAP_TO_CAP:
    {
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = *pc;
    }
    break;
  case CYL_CAPS_ARE_ROUND:
    {
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = cCylShaderBothCapsRound | cCylShaderInterpColor;
    }
    break;
  case CYL_CAPS_ARE_FLAT:
    {
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = cCylShaderBothCapsFlat | cCylShaderInterpColor;
    }
    break;
  case CYL_CAPS_ARE_CUSTOM:
    {
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      float *pcf = (float*)pc;
      const cCylCap cap1 = static_cast<cCylCap>(int(pcf[0]));
      const cCylCap cap2 = static_cast<cCylCap>(int(pcf[1]));
      dataPtrUB[0] =
          cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
    }
    break;
  case UB1_TO_INTERP:
    {
      bool interp = (pc[0] & cgo::draw::splitline::interpolation);
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = interp ? 1 : 0;
    }
    break;
  case UB1_INTERP_TO_CAP:
    {
      bool interp = (pc[0] & cgo::draw::splitline::interpolation);
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = (cCylShaderBothCapsRound | (interp ? cCylShaderInterpColor : 0));
    }
    break;
  case FLOAT1_TO_INTERP:
    {
      float interp = *((float*)pc);
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = (interp > .5f) ? 1 : 0;
    }
    break;
  case FLOAT1_INTERP_TO_CAP:
    {
      float interp = *((float*)pc);
      unsigned char *dataPtrUB = (unsigned char *)dataPtr;
      dataPtrUB[0] = (cCylShaderBothCapsRound | ((interp > .5f) ? cCylShaderInterpColor : 0));
    }
    break;
  case UB4_TO_UB4:
    {
      auto dataPtrUB = (unsigned char *)dataPtr;
      auto pcUB = (unsigned char *)pc;
      dataPtrUB[0] = pcUB[0];
      dataPtrUB[1] = pcUB[1];
      dataPtrUB[2] = pcUB[2];
      dataPtrUB[3] = pcUB[3];
      break;
    }
  case PICK_DATA_TO_PICK_DATA:
    {
      float *pcf;
      if (copyord < 0){
        pcf = (float*)pc;
      } else {
        pcf = (copyord * 2 + pick_data);
        if (nvert){
          pcf -= pstride;
        }
      }
      unsigned int index = CGO_get_uint(pcf);
      int bond = CGO_get_int(pcf+1);
      CGO_put_uint(ord * 2 + pick_data, index);
      CGO_put_int(ord * 2 + pick_data + 1, bond);
      has_pick_colorBS |= (1 << ord) ;
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
 * dataPtrs       : all data pointers for attributes (for interleaved, they are all the pointer to the one array)
 * attrOffset     : offsets of the attributes (for interleaved)
 *
 */
static
void copyAttributeForVertex(bool isInterleaved, int &nvert, AttribDesc &attribDesc,
                            const int vertexDataSize, std::vector<void*> &dataPtrs, std::vector<int> &attrOffset){
  int ord = attribDesc.order;
  void *dataPtr = dataPtrs[ord];
  unsigned char *pc = NULL;
  int attrSize = gl_sizeof(attribDesc.type_size) * attribDesc.type_dim;
  if (isInterleaved){
    dataPtr = (unsigned char*) dataPtr + nvert * vertexDataSize + attrOffset[ord];
    pc = (unsigned char*) dataPtr - vertexDataSize;
  } else {
    dataPtr = (unsigned char*) dataPtr + nvert * attrSize;
    pc = (unsigned char*) dataPtr - attrSize;
  }
  if (attribDesc.repeat_value && attribDesc.repeat_value_length){
    int pos = (nvert % attribDesc.repeat_value_length);
    pc = attribDesc.repeat_value + pos * attrSize;
    memcpy(dataPtr, pc, attrSize);
  } else {
    memcpy(dataPtr, pc, attrSize);
  }
}
/**
 * check all attributes (pick and non-pick) to see if they are specified in the CGO (I)
 * also checks to see if any picking is specified and sets has_picking argument
 * if any of the attributes are not specified and if a default_value is set (in the 
 * AttribDesc) then the associated vertex_attribute CGO operation is inserted into the
 * cgo that is passed in.
 *
 * I:           primitive CGO that is processed
 * attrData:    definition of attributes
 * pickData:    definition of pick attributes
 * cgo:         new cgo that could have vertex_attribute CGO operations added
 * has_picking: if there are any operations in the CGO (I) that specifies different values for picking.
 *              if has_picking is set, then a VBO for picking is generated in CGOConvertToShader()
 *
 */
static
void CheckAttributesForUsage(const CGO *I, AttribDataDesc &attrData, AttribDataDesc &pickData, CGO *cgo, bool &has_picking)
{
  size_t attrIdx = 0;

  // need to check attributes:
  //  - remove any that are not needed
  //  - add glVertexAttrib for those that are removed
  std::map<int, int> opToAttrUsed; // bitmask for each op to which attributes are set
  for (auto &attrDesc : attrData){
    auto attrOps = &attrDesc.attrOps;
    attrDesc.order = attrIdx++;
    for (auto attrOpIt = attrOps->begin(); attrOpIt!=attrOps->end(); ++attrOpIt){
      auto attrOp = &(*attrOpIt);
      if (opToAttrUsed.find(attrOp->op) == opToAttrUsed.end())
        opToAttrUsed[attrOp->op] = 1 << attrDesc.order;
      else
        opToAttrUsed[attrOp->op] |= 1 << attrDesc.order;
    }
  }
  // add picking ops (1 << attrIdx) i.e., any pick op is the last bit
  int pidx = 0;
  for (auto pickDataIt = pickData.begin(); pickDataIt!=pickData.end(); ++pickDataIt){
    auto pickDesc = &(*pickDataIt);
    auto pickOps = &pickDesc->attrOps;
    pickDesc->order = pidx++;
    for (auto pickOpIt = pickOps->begin(); pickOpIt!=pickOps->end(); ++pickOpIt){
      auto pickOp = &(*pickOpIt);
      pickOp->desc = pickDesc;
      if (opToAttrUsed.find(pickOp->op) == opToAttrUsed.end())
        opToAttrUsed[pickOp->op] = 1 << attrIdx;
      else
        opToAttrUsed[pickOp->op] |= 1 << attrIdx;
    }
  }
  size_t totAttrIdx = (1 << (attrIdx+1)) - 1;
  size_t allAttrIdxUsed = 0;
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    if (opToAttrUsed.find(op) != opToAttrUsed.end()){
      int attrUsed = opToAttrUsed[op];
      if (attrUsed & (1 << attrIdx)){ // if picking, need to check values, and take them out if cPickableNoPick
        switch (op){
        case cgo::draw::shadercylinder2ndcolor::op_code:
          if (reinterpret_cast<const cgo::draw::shadercylinder2ndcolor *>(pc)->pick_color_bond == cPickableNoPick)
            attrUsed ^= (1 << attrIdx);
          break;
        case cgo::draw::splitline::op_code:
          if (reinterpret_cast<const cgo::draw::splitline *>(pc)->bond == cPickableNoPick)
            attrUsed ^= (1 << attrIdx);
        }
      }
      allAttrIdxUsed |= attrUsed;
      if (allAttrIdxUsed == totAttrIdx)
        break;
    }
  }
  has_picking = allAttrIdxUsed & (1 << attrIdx);  // has_picking if the last bit is set

  if (allAttrIdxUsed != totAttrIdx){
    // go through any attributes that aren't used:
    //   - add associated vertex_attribute type (if default_value is set)
    //   - remove attribute from attrData description so that it isn't included in VBO
    AttribDataDesc attrDataNew;
    for (auto idx = 0; idx < attrIdx; ++idx){
      if (!attrData[idx].repeat_value && !(allAttrIdxUsed & (1 << idx))) {
        // attribute not used, need to create glVertexAttrib
        if (attrData[idx].default_value){
          // need to add glVertexAttrib CGO OP
          int attr_lookup_idx = I->G->ShaderMgr->GetAttributeUID(attrData[idx].attr_name);
          switch (attrData[idx].type_size){
          case GL_FLOAT:
            switch (attrData[idx].type_dim){
            case 1:
              cgo->add<cgo::draw::vertex_attribute_1f>(attr_lookup_idx, *(float*)attrData[idx].default_value);
              break;
            case 3:
              cgo->add<cgo::draw::vertex_attribute_3f>(attr_lookup_idx, attrData[idx].default_value);
              break;
            default:
              std::cerr << "\tNOT IMPLEMENTED: attrData[idx].type_size=" << attrData[idx].type_size << " attrData[idx].type_dim=" << attrData[idx].type_dim << std::endl;
            }
            break;
          case GL_UNSIGNED_BYTE:
            switch (attrData[idx].type_dim){
            case 1:
              {
                float val;
                unsigned char valuc = *attrData[idx].default_value;
                if (attrData[idx].data_norm){
                  val = CLAMPVALUE(valuc / 255.f, 0.f, 1.f);
                } else {
                  val = (float)valuc;
                }
                cgo->add<cgo::draw::vertex_attribute_1f>(attr_lookup_idx, val);
              }
              break;
            case 4:
              cgo->add<cgo::draw::vertex_attribute_4ub>(attr_lookup_idx, attrData[idx].default_value);
              break;
            default:
              std::cerr << "\tNOT IMPLEMENTED: attrData[idx].type_size=" << attrData[idx].type_size << " attrData[idx].type_dim=" << attrData[idx].type_dim << std::endl;
            }
          }
        }
      } else {
        attrDataNew.push_back(attrData[idx]);
      }
    }
    attrData.swap(attrDataNew);  // only keep attributes that are used
  }
}

/**
 * Populates two structures and sets vertsperpickinfo
 * opToCntPer           : CGO op to how many vertices are generated for each op.
 * opToOrderedAttribOps : CGO op to an ordered map of AttribOps that define how we operate on
 *                        the attribute arrays for each CGO op.
 * attrData             : definition of attributes
 * pickData             : definition of pick attributes
 * vertsperpickinfo     : number of vertices per pick (only used when 
 * 
 */
static
void PopulateOpsIntoStructuresForConversion(std::map<int,int> &opToCntPer,
                                            std::map< int, std::map<int, AttribOp*> > &opToOrderedAttribOps, 
                                            AttribDataDesc &attrData, AttribDataDesc &pickData,
                                            int &vertsperpickinfo, const bool has_picking){
  size_t attrIdx = 0;
  for (auto &attrDesc : attrData){
    auto attrOps = &attrDesc.attrOps;
    attrDesc.order = attrIdx++;
    for (auto attrOpIt = attrOps->begin(); attrOpIt!=attrOps->end(); ++attrOpIt){
      auto attrOp = &(*attrOpIt);
      attrOp->desc = &attrDesc;
      if (attrOp->copyFromAttr >= 0){
        attrOp->copyAttribDesc = &attrData[attrOp->copyFromAttr];
      }
      if (attrOp->incr_vertices > 0){
        if (!vertsperpickinfo){
          vertsperpickinfo = attrOp->incr_vertices;
        } else {
          if (attrOp->incr_vertices != vertsperpickinfo){
            std::cerr << "WARNING: attrOp->incr_vertices set to multiple values, vertsperpickinfo=" << vertsperpickinfo << " attrOp->incr_vertices=" << attrOp->incr_vertices << " : picking might get confused" << std::endl;
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
  if (has_picking){
    for (auto pickDataIt = pickData.begin(); pickDataIt!=pickData.end(); ++pickDataIt){
      auto pickDesc = &(*pickDataIt);
      auto pickOps = &pickDesc->attrOps;
      for (auto pickOpIt = pickOps->begin(); pickOpIt!=pickOps->end(); ++pickOpIt){
        auto pickOp = &(*pickOpIt);
        if (pickOp->copyFromAttr >= 0){
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
 * attrData:            definition of attributes that are accumulated and put into the VBO
 * pickData:            definition of pick attributes that accumulate pick data and put into the interleaved picking VBO
 * mode:                which openGL mode to use for rendering (e.g., GL_POINTS, GL_TRIANGLES, GL_LINE_STRIPS, etc.)
 * layout:              SEPARATE, SEQUENTIAL, or INTERLEAVED : how the VBO is layed out in memory (default: INTERLEAVED)
 * check_attr_for_data: if true, this function checks whether all attributes are used, and if any are not specified,
 *                      and default values are defined in AttribDesc, then glVertexAttrib CGO ops are created on the
 *                      returned CGO (default: true)
 * idx_array:           if specified, this array is used to specify indices for each fragment in an Indexed Buffer
 *                      and glDrawElements is used to render (instead of glDrawArrays)
 * nvertsperfrag:       in conjunction with idx_array, how many vertices for each set of indices (default: 0)
 *                      both idx_array and nvertsperfrag are used to specify vertices and geometry for each fragment,
 *                      such as the box for cylinders
 * nfragspergroup:      Currently, this represents the number of fragments that are inside of a group.  For example,
 *                      crosses as cylinders (CGOConvertCrossesToCylinderShader) have 36 indices per fragment (nvertsperfrag=36)
 *                      and 3 fragments per group (nfragspergroup=3) (i.e., one for each line).
 *                      note: idx_array, nvertsperfrag and nfragspergroup should probably be moved to AttrOp in
 *                            the future to support different types of fragments.
 *
 * returns a CGO that consists of vertex_attrib_* (e.g., glVertexAttrib) and a custom CGO operation that calls
 * either glDrawArrays or glDrawElements. It also supports picking if specified in the pickData.
 *
 */
CGO *CGOConvertToShader(const CGO *I, AttribDataDesc &attrData, AttribDataDesc &pickData, int mode, 
                        const VertexBuffer::buffer_layout layout, bool check_attr_for_data,
                        int *idx_array, int nvertsperfrag, int nfragspergroup){
  CGO *cgo;
  int ok = true;
  bool isInterleaved = (layout == VertexBuffer::INTERLEAVED);
  bool has_picking = true;
  std::map<std::string, AttribDesc*> attrToDesc;

  cgo = CGONew(I->G);
  cgo->use_shader = true;

  if (check_attr_for_data){
    CheckAttributesForUsage(I, attrData, pickData, cgo, has_picking);
  } else {
    // if attributes aren't checked, still need to set pick order and desc pointer
    int pidx = 0;
    for (auto &pickDesc : pickData){
      pickDesc.order = pidx++;
      for (auto &pickOp : pickDesc.attrOps){
        pickOp.desc = &pickDesc;
      }
    }
  }
  std::map<int,int> opToCntPer;
  std::map< int, std::map<int, AttribOp*> > opToOrderedAttribOps ;

  int vertsperpickinfo = 0;
  PopulateOpsIntoStructuresForConversion(opToCntPer, opToOrderedAttribOps, attrData,
                                         pickData, vertsperpickinfo, has_picking);


  // Populate these variables used for accumulating and setting VBO data arrays
  int vertexDataSize = 0;
  std::vector<int> attrSizes;
  std::vector<int> attrOffset; // for interleaved
  int curoffset = 0;      // for interleaved
  size_t attrIdx = 0;
  for (auto &attrDesc : attrData){
    attrDesc.order = attrIdx++;
    int attrSize = gl_sizeof(attrDesc.type_size) * attrDesc.type_dim;
    attrSizes.push_back(attrSize);
    vertexDataSize += attrSize;

    if (isInterleaved){
      attrOffset.push_back(curoffset);
      curoffset += attrSize;
    } else {
      attrOffset.push_back(0);  // no offset when not interleaved
    }
    attrToDesc[attrDesc.attr_name] = &attrDesc;
  }

  // populate funcData.attrib
  for (auto &attrDesc : attrData)
    for (auto &attrop : attrDesc.attrOps)
      for (auto &funcData : attrop.funcDataConversions){
        if (attrToDesc.find(funcData.attribName) != attrToDesc.end()){
          funcData.attrib = attrToDesc[funcData.attribName];
        }
      }

  // Since some cards require word-aligned strides (e.g., ATI)
  // we need to make sure it is word-aligned.  If it isn't, and you
  // want to save memory, maybe SEQUENTIAL is a better option
  if (vertexDataSize % 4){
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
  // Generate VBOs: both for vertex data (vbo) and picking color data (pickvbo) if necessary
  // - pickvbo is interleaved so that for multiple channels, pick data for each vertex
  //   is contiguous
  VertexBuffer *vbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(layout);
  VertexBuffer *pickvbo = NULL;
  if (pickDataSize){
    pickvbo = I->G->ShaderMgr->newGPUBuffer<VertexBuffer>(VertexBuffer::SEQUENTIAL, GL_DYNAMIC_DRAW);
    pickvbohash = pickvbo->get_hash_id();
  }

  // adding SPECIAL OPERATIONS (for now) before adding custom OP
  // - for now, this is the only operation that needs to be passed to the new CGO
  for (auto it = I->begin(); !it.is_stop(); ++it) {
    auto pc = it.data();
    int op = it.op_code();
    switch (op){
    case CGO_SPECIAL:
      cgo->add_to_cgo(op, pc);
      break;
    }
  }

  // defines how many passes we have
  // We *could* set npickcolattr=1 if we had access to PickColorConverter here and getTotalBits() == 32
  const int npickcolattr = SHADER_PICKING_PASSES_MAX;
  /* for picking, we need to mask the attributes and bind the buffer in the the picking VBO */
  int pl = 0;
  int npickattr = pickData.size();
  for (auto &pickDesc : pickData){
    cgo->add<cgo::draw::mask_attribute_if_picking>(I->G->ShaderMgr->GetAttributeUID(pickDesc.attr_name), vbo->get_hash_id());
    if (has_picking){
      cgo->add<cgo::draw::bind_vbo_for_picking>(pickvbohash, pl++, npickattr);
    } else {
      /* if no picking, should render black */
      static unsigned char zerocolor[] { 0,0,0,0 };
      cgo->add<cgo::draw::vertex_attribute_4ub_if_picking>(I->G->ShaderMgr->GetAttributeUID(pickDesc.attr_name), zerocolor);
    }
  }
  size_t iboid = 0;
  int num_total_indexes = 0;
  if (nvertsperfrag){
    GL_C_INT_TYPE *vertexIndices; 
    int nfrags = nfragspergroup * ntotalverts/vertsperpickinfo;
    int nvertsperindivfrag = vertsperpickinfo/nfragspergroup;
    num_total_indexes = nfrags * nvertsperfrag;
    vertexIndices = pymol::calloc<GL_C_INT_TYPE>(num_total_indexes);
    int idxpl=0;
    // using vertsperpickinfo as verts per frag
    for (int cnt = 0, vpl = 0; cnt < nfrags; ++cnt){
      for (int idx_array_pl = 0; idx_array_pl < nvertsperfrag; ++idx_array_pl){
        vertexIndices[idxpl] = idx_array[idx_array_pl] + vpl;
        idxpl++;
      }
      vpl+=nvertsperindivfrag;
    }
    IndexBuffer * ibo = I->G->ShaderMgr->newGPUBuffer<IndexBuffer>();
    ok &= ibo->bufferData({
        BufferDesc( GL_C_INT_ENUM, sizeof(GL_C_INT_TYPE) * num_total_indexes, vertexIndices )
          });
    FreeP(vertexIndices);
    iboid = ibo->get_hash_id();
  }

  // pick_data is interleaved if more than one attribute
  float * pick_data = cgo->add<cgo::draw::custom>(mode, ntotalverts, vbo->get_hash_id(), pickvbohash, vertsperpickinfo, pickDataSize, iboid, num_total_indexes);
  void *allData = malloc(ntotalverts * vertexDataSize);
  std::vector<void*> dataPtrs;
  std::vector<int> repeat_attr_idx;
  int allAttrBS = 0;

  // Initialize first entry in array(s) with default values and populate dataPtrs
  if (isInterleaved){
    int pl = 0;
    auto attrDataIt = attrData.begin();
    auto attrOffsetIt = attrOffset.begin();
    for (; attrDataIt!=attrData.end() && attrOffsetIt!=attrOffset.end(); ++attrDataIt, ++attrOffsetIt){
      auto attrDesc = &(*attrDataIt);
      if (attrDesc->repeat_value){
        repeat_attr_idx.push_back(pl);
      } else {
        allAttrBS |= (1 << attrDesc->order);
      }
      auto attrOffset = *attrOffsetIt;
      unsigned char *first_value = NULL;
      first_value = (attrDesc->default_value ? attrDesc->default_value : 
                     (attrDesc->repeat_value ? attrDesc->repeat_value : NULL));
      if (first_value){
        int attrSize = gl_sizeof(attrDesc->type_size) * attrDesc->type_dim;
        memcpy(((unsigned char*)allData)+attrOffset, first_value, attrSize);
      }
      dataPtrs.push_back((void*)allData);
      ++pl;
    }
  } else {
    void *curAllDataPtr = (void*)allData;
    int pl = 0;
    for (auto &attrDesc : attrData){
      if (attrDesc.repeat_value){
        repeat_attr_idx.push_back(pl);
      } else {
        allAttrBS |= (1 << attrDesc.order);
      }
      dataPtrs.push_back(curAllDataPtr);
      unsigned char *first_value = NULL;
      first_value = (attrDesc.default_value ? attrDesc.default_value : 
                     (attrDesc.repeat_value ? attrDesc.repeat_value : NULL));
      if (first_value){
        memcpy((unsigned char*)curAllDataPtr, first_value, attrSizes[pl]);
      }
      curAllDataPtr = ((unsigned char*)curAllDataPtr) + ntotalverts * attrSizes[pl];
      ++pl;
    }
  }

  int nvert = 0;
  int attrBS = 0;
  int allPickAttrBS = (1 << pickData.size()) - 1;
  int has_pick_colorBS = allPickAttrBS;
  for (int pi = 0; pi<pickDataSize; ++pi){
    CGO_put_uint(2*pi + pick_data, 0);
    CGO_put_int(2*pi + pick_data+1, cPickableNoPick);
  }
  bool cont = true; // need to break while statement as well if past all vertices

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
    if (opToOrderedAttribOps.find(op) != opToOrderedAttribOps.end()){
      std::map<int, AttribOp*> *attribOpsInOrder = &opToOrderedAttribOps[op];
      for (auto attribOpIt : *attribOpsInOrder){
        AttribOp *attribOp = attribOpIt.second;
        int ord = attribOp->desc->order;
	cont = nvert < ntotalverts;
	if (!cont)
	  break;
        copyAttributeForOp(isInterleaved, nvert, attribOp, vertexDataSize, dataPtrs, attrOffset, pc, pick_data, has_pick_colorBS, pstride);
        if (ord >= 0)  // picking is negative, has_pick_colorBS is used instead
          attrBS |= (1 << ord);
        else
          std::cout << "   ord=%d\n" << ord << std::endl;
        if (attribOp->incr_vertices){
          if (has_pick_colorBS!=allPickAttrBS){
            // copy pick colors that haven't been set from previous vertex
            for (int pi = 0; pi<pickDataSize; ++pi){
              if (has_pick_colorBS ^ (1 << pi)){
                CGO_put_uint(pick_data, CGO_get_uint(pick_data-pstride));
                CGO_put_int(pick_data + 1, CGO_get_int(pick_data-pstride+1));
              }
              pick_data += 2;
            }
          } else {
            pick_data += pstride;
          }
          has_pick_colorBS = 0;
          if (!nvert && attrBS!=allAttrBS){
            // for the first vertex, all attributes should be set
            for (auto idx = 0; idx < attrData.size(); ++idx){
              if (!(attrBS & (1 << idx))) {
                if (!attrData[idx].default_value){
                  std::cerr << "WARNING: attribute #" << idx <<
                    " (" << attrData[idx].attr_name << ") not set for first"
                    " vertex and does not have default value" << std::endl;
                }
              }
            }
          }
          if (nvert && attrBS!=allAttrBS){
            // for each vertex that hasn't been written for the current vertex, copy it from the previous vertex
            for (auto idx = 0; idx < attrData.size(); ++idx){
              if (!(attrBS & (1 << idx))) {
                copyAttributeForVertex(isInterleaved, nvert, attrData[idx], vertexDataSize, dataPtrs, attrOffset);
              }
            }
          }
          attrBS = 0;

          // creating new vertices, all attribute data should be copied into new vertex.
          for (int nxt = 0; nxt < attribOp->incr_vertices; ++nxt){
            /* for now, always copy values from previous */
            ++nvert;
            if (nvert < ntotalverts){
              {
                // last should not need to copy into next, since we call copyAttributeForVertex (above)
                // for all attributes that haven't been set
                // - note: it might be faster (especially for interleaved) to copy into the next
                //         vertex, then the above copyAttributeForVertex() would not be needed
                if (isInterleaved){
                  void *dest = ((unsigned char*)allData)+vertexDataSize*nvert;
                  memcpy(dest, ((unsigned char*)dest) - vertexDataSize, vertexDataSize);
                } else {
                  auto dataPtrIt = dataPtrs.begin();
                  auto attrDataIt = attrData.begin();
                  for (; attrDataIt!=attrData.end() && dataPtrIt!=dataPtrs.end(); ++attrDataIt, ++dataPtrIt){
                    auto attrDesc = &(*attrDataIt);
                    auto dataPtr = *dataPtrIt;
                    int attrSize = gl_sizeof(attrDesc->type_size) * attrDesc->type_dim;
                    void *dest = ((unsigned char*)dataPtr)+attrSize*nvert;
                    memcpy((unsigned char*)dest, ((unsigned char*)dest) - attrSize, attrSize);
                  }
                }
              }
              // always copy repeat attributes
              if (!repeat_attr_idx.empty()){
                for (auto ridx = repeat_attr_idx.begin(); ridx != repeat_attr_idx.end(); ++ridx){
                  copyAttributeForVertex(isInterleaved, nvert, attrData[*ridx], vertexDataSize, dataPtrs, attrOffset);
                }
              }
            }
            if (!attribOp->funcDataConversions.empty()){
              // for all attributes, call the funcDataConversion() if defined 
              int nvert_m_1 = nvert - 1;
              for (auto funcData : attribOp->funcDataConversions){
                auto funcAttrib = funcData.attrib;
                auto order = funcAttrib->order;
                if (isInterleaved){
                  unsigned char *dest = ((unsigned char*)allData)+vertexDataSize*nvert_m_1 + attrOffset[order];
                  funcData.funcDataConversion(dest, pc, funcData.funcDataGlobalArg, nxt);
                } else {
                  auto dataPtr = dataPtrs[order];
                  void *dest = ((unsigned char*)dataPtr)+attrSizes[order]*nvert_m_1;
                  funcData.funcDataConversion(dest, pc, funcData.funcDataGlobalArg, nxt);
                }
              }
            }

          }
        }
      }
    }
  }

  /* Generate Pick Buffers with all pick attributes (if necessary) */
  if (pickvbo){
    BufferDataDesc pickBufferData;
    for (int i=0; i < npickcolattr; i++){
      for (auto &pickDesc : pickData){
	int pickSize = gl_sizeof(pickDesc.type_size) * pickDesc.type_dim;
        pickBufferData.push_back(BufferDesc(pickDesc.attr_name, pickDesc.type_size,
                                            pickDesc.type_dim, pickSize * nvert, NULL, pickDesc.data_norm));
      }
    }
    pickvbo->bufferData(std::move(pickBufferData));
  }

  /* Generate VBO Buffers with all pick attributes based on the VertexBuffer type SEPARATE/SEQUENTIAL/INTERLEAVED*/
  BufferDataDesc bufferData;
  switch (layout){
  case VertexBuffer::SEPARATE:
  case VertexBuffer::SEQUENTIAL:
    {
      auto attrDataIt = attrData.begin();
      auto dataPtrIt = dataPtrs.begin();
      auto attrSizeIt = attrSizes.begin();
        for (; attrDataIt!=attrData.end() && 
               dataPtrIt!=dataPtrs.end() && 
               attrSizeIt!=attrSizes.end(); ++attrDataIt, ++dataPtrIt, ++attrSizeIt){
        auto attrDesc = &(*attrDataIt);
        auto dataPtr = *dataPtrIt;
        auto attrSize = *attrSizeIt;
        bufferData.push_back(BufferDesc(attrDesc->attr_name, attrDesc->type_size, 
                                        attrDesc->type_dim, nvert * attrSize, dataPtr, attrDesc->data_norm));
      }
        vbo->bufferData(std::move(bufferData));
    break;
    }
    break;
  case VertexBuffer::INTERLEAVED:
    {
      auto attrDataIt = attrData.begin();
      auto attrOffsetIt = attrOffset.begin();
      for (; attrDataIt!=attrData.end() && attrOffsetIt!=attrOffset.end(); ++attrDataIt, ++attrOffsetIt){
        auto attrDesc = &(*attrDataIt);
        auto offset = *attrOffsetIt;
        bufferData.push_back(BufferDesc(attrDesc->attr_name, attrDesc->type_size,
                                        attrDesc->type_dim, offset, attrDesc->data_norm));
      }
      vbo->bufferData(std::move(bufferData), (const void *)allData,
                      (size_t)(nvert*vertexDataSize), (size_t)vertexDataSize);
    break;
    }
  }
  free(allData);

  CGOStop(cgo);
  return cgo;
}

// CGOCheckSplitLineInterpolationIsSame: 
//   - returns true if always the same
//   - returns false if not always the same
bool CGOCheckSplitLineInterpolationIsSame(const CGO *I, bool &interp_value){
  bool interp_value_first = false;
  bool interp_value_is_set = false;

  for (auto it = I->begin(); !it.is_stop(); ++it) {
    switch (it.op_code()) {
    case cgo::draw::splitline::op_code:
      interp_value = (it.cast<cgo::draw::splitline>()->flags & cgo::draw::splitline::interpolation);
      break;
    case CGO_INTERPOLATED:
      interp_value = it.cast<float>()[0] > 0.5f;
      break;
    default:
      continue;
    }
    if (!interp_value_is_set){
      interp_value_first = interp_value;
      interp_value_is_set = true;
    } else if (interp_value != interp_value_first){
      return false;
    }
  }
  return true;
}

// CGOCheckShaderCylinderCapInfoIsSame: 
//   - returns true if always the same
//   - returns false if not always the same
static
bool CGOCheckShaderCylinderCapInfoIsSame(const CGO *I, unsigned char &cap_value){
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
    case cgo::draw::custom_cylinder::op_code:
      {
        auto cc = it.cast<cgo::draw::custom_cylinder>();
        const cCylCap cap1 = cc->get_cap1();
        const cCylCap cap2 = cc->get_cap2();
        cap_value =
            cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
      }
      break;
    case cgo::draw::custom_cylinder_alpha::op_code:
      {
        auto cc = it.cast<cgo::draw::custom_cylinder_alpha>();
        const cCylCap cap1 = cc->get_cap1();
        const cCylCap cap2 = cc->get_cap2();
        cap_value =
            cyl_shader_bits_from_caps(cap1, cap2) | cCylShaderInterpColor;
      }
      break;
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

CGO *CGOConvertToTrilinesShader(const CGO *I, CGO *addTo, bool add_color){
  PyMOLGlobals *G = I->G;

  AttribDataOp vertexOps =
    { { CGO_LINE,       1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::line, vertex1), 0 },
      { CGO_SPLITLINE,  2, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::splitline, vertex1), 0 } };
  AttribDataOp vertexOtherOps =
    { { CGO_LINE,       2, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::line, vertex2), 6 },
      { CGO_SPLITLINE,  5, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::splitline, vertex2), 6 } };
  AttribDataOp colorOps =
    { { CGO_COLOR,      0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      0, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SPLITLINE,  6, UB3_TO_UB3,            offsetof(cgo::draw::splitline, color2) } };
  AttribDataOp color2Ops =
    { { CGO_COLOR,      1, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      1, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SPLITLINE,  3, UB3_TO_UB3,            offsetof(cgo::draw::splitline, color2) } };
  AttribDataOp extraPickColorOps =
    { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 },
      { CGO_SPLITLINE,  7, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::splitline, index), 0 } };
  AttribDataOp extraPickColor2Ops =
    { { CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0 },
      { CGO_SPLITLINE,  4, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::splitline, index), 0 } };
  AttribDataDesc pickDesc =
    { { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps },
      { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColor2Ops }};
  AttribDataDesc attrDesc =
    { { "a_Vertex", GL_FLOAT, 3, GL_FALSE, vertexOps },
      { "a_OtherVertex", GL_FLOAT, 3, GL_FALSE, vertexOtherOps },
      { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, colorOps },
      { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, color2Ops },
      { "a_UV", GL_UNSIGNED_BYTE, 1, GL_FALSE } };

  if (add_color){
    static unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
    attrDesc[2].default_value = default_color;
    attrDesc[3].default_value = default_color;
  }
  AttribDesc *uvdesc = &attrDesc[attrDesc.size()-1];
  uvdesc->repeat_value_length = 6;
  static unsigned char uv_bits[] = { 1, 3, 0, 3, 2, 1 };
  uvdesc->repeat_value = uv_bits;
  
  bool interp_same, interp_value;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))){
    addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_interpolate"), interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps =
      { { CGO_SPLITLINE, 1, UB1_TO_INTERP, offsetof(cgo::draw::splitline, flags), 0 } };
    // need to add a_interpolate attribute
    attrDesc.push_back({ "a_interpolate", GL_UNSIGNED_BYTE, 1, GL_FALSE, interpOps } );
  }
  if (!add_color){
    attrDesc.erase(attrDesc.begin()+2); // a_Color
    attrDesc.erase(attrDesc.begin()+2); // a_Color2
  }

  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES, VertexBuffer::INTERLEAVED);
}

CGO *CGOConvertToLinesShader(const CGO *I, CGO *addTo, bool add_color){
  /* Lines that pass in two vertices per line */
  PyMOLGlobals *G = I->G;
  AttribDataOp vertexOps =
    { { CGO_LINE,       1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::line, vertex1), 1 },
      { CGO_SPLITLINE,  2, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::splitline, vertex1), 1 },
      { CGO_LINE,       2, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::line, vertex2), 1 },
      { CGO_SPLITLINE,  5, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::splitline, vertex2), 1 } };
  AttribDataOp colorOps =
    { { CGO_COLOR,      0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      0, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SPLITLINE,  3, UB3_TO_UB3,            offsetof(cgo::draw::splitline, color2) } };
  AttribDataOp extraPickColorOps =
    { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 },
      { CGO_SPLITLINE,  4, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::splitline, index), 0 } };
  
  AttribDataDesc pickDesc =
    { { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps } };
  AttribDataDesc attrDesc =
    { { "a_Vertex", GL_FLOAT, 3, GL_FALSE, vertexOps },
      { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, colorOps } };
  if (add_color){
    static unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
    attrDesc[1].default_value = default_color;
  }
  bool interp_same, interp_value;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))){
    addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_interpolate"), interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps =
      { { CGO_SPLITLINE, 1, UB1_TO_INTERP, offsetof(cgo::draw::splitline, flags), 0 } };
    // need to add a_interpolate attribute
    attrDesc.push_back({ "a_interpolate", GL_UNSIGNED_BYTE, 1, GL_FALSE, interpOps } );
  }
#ifndef PURE_OPENGL_ES_2
  {
    attrDesc.push_back({ "a_line_position", GL_UNSIGNED_BYTE, 1, GL_FALSE } );
    AttribDesc *lpdesc = &attrDesc[attrDesc.size()-1];
    lpdesc->repeat_value_length = 2;
    static unsigned char flip_bits[] = { 0, 1 };
    lpdesc->repeat_value = flip_bits;
  }
#endif
  if (!add_color){
    attrDesc.erase(attrDesc.begin()+1); // a_Color
  }

  return CGOConvertToShader(I, attrDesc, pickDesc, GL_LINES, VertexBuffer::INTERLEAVED);
}

CGO *CGOConvertLinesToCylinderShader(const CGO *I, CGO *addTo, bool add_color){
  /* Lines that pass in two vertices per line */
  PyMOLGlobals *G = I->G;

  AttribDataOp vertex1Ops =
    { { CGO_LINE,       1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::line, vertex1), 0 },
      { CGO_SPLITLINE,  2, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::splitline, vertex1), 0 } };
  AttribDataOp vertex2Ops =
    { { CGO_LINE,       2, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::line, vertex2), 8 },
      { CGO_SPLITLINE,  5, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::splitline, vertex2), 8 } };
  static AttribDataOp colorOps =
    { { CGO_COLOR,      0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      0, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SPLITLINE,  6, UB3_TO_UB3,            offsetof(cgo::draw::splitline, color2) } };
  static AttribDataOp color2Ops =
    { { CGO_COLOR,      1, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      1, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SPLITLINE,  3, UB3_TO_UB3,            offsetof(cgo::draw::splitline, color2) } };

  AttribDataDesc attrDesc = { { "attr_vertex1", GL_FLOAT, 3, GL_FALSE, vertex1Ops },
                              { "attr_vertex2", GL_FLOAT, 3, GL_FALSE, vertex2Ops },
                              { "a_Color",  GL_UNSIGNED_BYTE, 4, GL_TRUE, colorOps },
                              { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, color2Ops },
                              { "attr_radius", GL_FLOAT, 1, GL_FALSE } };
  AttribDesc *fdesc;
  static unsigned char cyl_flags[] = { 0, 4, 6, 2, 1, 5, 7, 3 }; // right(4)/up(2)/out(1)

  attrDesc.push_back( { "attr_flags", GL_UNSIGNED_BYTE, 1, GL_FALSE } ) ;
  fdesc = &attrDesc[attrDesc.size()-1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  if (add_color){
    static unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
    fdesc = &attrDesc[2];
    fdesc->default_value = default_color;
    fdesc = &attrDesc[3];
    fdesc->default_value = default_color;
  }

  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*)&default_radius;

  int box_indices[36] = { // box indices 
    0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 
    0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4 };
  int *box_indices_ptr = NULL;
  box_indices_ptr = box_indices;

  bool interp_same, interp_value = false;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))){
    addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_cap"), (cCylShaderBothCapsRound | (interp_value ? cCylShaderInterpColor : 0)));
  } else {
    AttribDataOp interpOps =
      { { CGO_SPLITLINE, 1, UB1_INTERP_TO_CAP, offsetof(cgo::draw::splitline, flags), 0 } };
    // need to add a_cap attribute
    attrDesc.push_back({ "a_cap", GL_UNSIGNED_BYTE, 1, GL_FALSE, interpOps } );
  }

  if (!add_color){
    attrDesc.erase(attrDesc.begin()+2); // attr_colors
    attrDesc.erase(attrDesc.begin()+2); // attr_colors2
  }

  AttribDataOp extraPickColorOps = { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 },
                                     { CGO_SPLITLINE,  7, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::splitline, index), 0 } };
  AttribDataOp extraPickColor2Ops = { { CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0 },
                                      { CGO_SPLITLINE,  4, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::splitline, index), 0 } };
  AttribDataDesc pickDesc = { { "a_Color",  GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps },
                              { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColor2Ops }};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES, VertexBuffer::INTERLEAVED, true, box_indices_ptr, 36);
}

struct CrossSizeData {
  float cross_size;
  bool forward;
  CrossSizeData(float _cross_size, bool _forward) : cross_size(_cross_size), forward(_forward){ }
};

static void CrossVertexConversion(void *varData, const float * pc, void *crossData, int idx){
  CrossSizeData *csd = (CrossSizeData*)crossData;
  int idxpl = idx / 8; // X Y or Z
  float *varDataF = ((float*)varData);
  varDataF[idxpl] += (csd->forward ? csd->cross_size : -csd->cross_size);
}

CGO *CGOConvertCrossesToCylinderShader(const CGO *I, CGO *addTo, float cross_size_arg){
  /* Lines that pass in two vertices per line */
  PyMOLGlobals *G = I->G;
  AttribDataOp vertex1Ops =
    { { CGO_VERTEX_CROSS,       1, FLOAT3_TO_FLOAT3,      0, 0 } };
  AttribDataOp vertex2Ops =
    { { CGO_VERTEX_CROSS,       2, FLOAT3_TO_FLOAT3,      0, 3 * 8, 0 } };

  static AttribDataOp colorOps =
    { { CGO_COLOR,      0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      0, FLOAT1_TO_UB_4TH,      0 } };
  static AttribDataOp color2Ops =
    { { CGO_COLOR,      1, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      1, FLOAT1_TO_UB_4TH,      0 } };

  CrossSizeData crossData[] = { { cross_size_arg, false }, { cross_size_arg, true } };
  AttribDataDesc attrDesc = { { "attr_vertex1", GL_FLOAT, 3, GL_FALSE, vertex1Ops },
                              { "attr_vertex2", GL_FLOAT, 3, GL_FALSE, vertex2Ops },
                              { "a_Color",  GL_UNSIGNED_BYTE, 4, GL_TRUE, colorOps },
                              { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, color2Ops },
                              { "attr_radius", GL_FLOAT, 1, GL_FALSE } };
  attrDesc.reserve(10);
  attrDesc[1].attrOps[0].funcDataConversions.push_back( { CrossVertexConversion, &crossData[0], "attr_vertex1" } );
  attrDesc[1].attrOps[0].funcDataConversions.push_back( { CrossVertexConversion, &crossData[1], "attr_vertex2" } );

  AttribDesc *fdesc;
  static unsigned char cyl_flags[] = { 0, 4, 6, 2, 1, 5, 7, 3 }; // right(4)/up(2)/out(1)

  attrDesc.push_back( { "attr_flags", GL_UNSIGNED_BYTE, 1, GL_FALSE } ) ;
  fdesc = &attrDesc[attrDesc.size()-1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
  fdesc = &attrDesc[2];
  fdesc->default_value = default_color;
  fdesc = &attrDesc[3];
  fdesc->default_value = default_color;
  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*)&default_radius;

  int box_indices[36] = { // box indices 
    0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 
    0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4 };
  int *box_indices_ptr = NULL;
  box_indices_ptr = box_indices;

  addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_cap"), cCylShaderBothCapsRound);

  AttribDataOp extraPickColorOps = { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 } };
  AttribDataOp extraPickColor2Ops = { { CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0 } };
  AttribDataDesc pickDesc = { { "a_Color",  GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps },
                              { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColor2Ops }};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES, VertexBuffer::INTERLEAVED, true, box_indices_ptr, 36, 3);
}

struct CrossSizeDataLines {
  float cross_size;
  CrossSizeDataLines(float _cross_size) : cross_size(_cross_size){ }
};

static void CrossVertexConversionLines(void *varData, const float * pc, void *crossData, int idx){
  CrossSizeDataLines *csd = (CrossSizeDataLines*)crossData;
  int idxpl = idx / 2; // X Y or Z
  bool forward = idx % 2;
  float *varDataF = ((float*)varData);
  varDataF[idxpl] += (forward ? csd->cross_size : -csd->cross_size);
}

CGO *CGOConvertCrossesToLinesShader(const CGO *I, CGO *addTo, float cross_size_arg){
  /* Lines that pass in two vertices per line */
  PyMOLGlobals *G = I->G;
  AttribDataOp vertexOps =
    { { CGO_VERTEX_CROSS,       1, FLOAT3_TO_FLOAT3,      0, 6 } };  // 6 vertices for a cross
  AttribDataOp colorOps =
    { { CGO_COLOR,      0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      0, FLOAT1_TO_UB_4TH,      0 } };
  AttribDataOp extraPickColorOps =
    { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 } };

  AttribDataDesc pickDesc =
    { { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps } };
  AttribDataDesc attrDesc =
    { { "a_Vertex", GL_FLOAT, 3, GL_FALSE, vertexOps },
      { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, colorOps } };
  unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
  attrDesc[1].default_value = default_color;

  CrossSizeDataLines crossData = { cross_size_arg };
  attrDesc[0].attrOps[0].funcDataConversions.push_back( { CrossVertexConversionLines, &crossData, "a_Vertex" } );

  bool interp_same, interp_value = false;
  if ((interp_same = CGOCheckSplitLineInterpolationIsSame(I, interp_value))){
    addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_interpolate"), interp_value ? 1.f : 0.f);
  } else {
    AttribDataOp interpOps =
      { { CGO_SPLITLINE, 1, UB1_TO_INTERP, offsetof(cgo::draw::splitline, flags), 0 } };
    // need to add a_interpolate attribute
    attrDesc.push_back({ "a_interpolate", GL_UNSIGNED_BYTE, 1, GL_FALSE, interpOps } );
  }
#ifndef PURE_OPENGL_ES_2
  {
    attrDesc.push_back({ "a_line_position", GL_UNSIGNED_BYTE, 1, GL_FALSE } );
    AttribDesc *lpdesc = &attrDesc[attrDesc.size()-1];
    lpdesc->repeat_value_length = 2;
    static unsigned char flip_bits[] = { 0, 1 };
    lpdesc->repeat_value = flip_bits;
  }
#endif
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_LINES, VertexBuffer::INTERLEAVED);
}

static void CrossVertexConversionTrilines(void *varData, const float * pc, void *crossData, int idx){
  CrossSizeData *csd = (CrossSizeData*)crossData;
  int idxpl = idx / 6; // X Y or Z
  float *varDataF = ((float*)varData);
  varDataF[idxpl] += (csd->forward ? csd->cross_size : -csd->cross_size);
}

CGO *CGOConvertCrossesToTrilinesShader(const CGO *I, CGO *addTo, float cross_size_arg){
  PyMOLGlobals *G = I->G;

  AttribDataOp vertexOps =
    { { CGO_VERTEX_CROSS,       1, FLOAT3_TO_FLOAT3,      0, 0 } };
  AttribDataOp vertexOtherOps =
    { { CGO_VERTEX_CROSS,       2, FLOAT3_TO_FLOAT3,      0, 6 * 3 } };
  AttribDataOp colorOps =
    { { CGO_COLOR,      0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      0, FLOAT1_TO_UB_4TH,      0 } };
  AttribDataOp color2Ops =
    { { CGO_COLOR,      1, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,      1, FLOAT1_TO_UB_4TH,      0 } };
  AttribDataOp extraPickColorOps =
    { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 } };
  AttribDataOp extraPickColor2Ops =
    { { CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0 } };
  AttribDataDesc pickDesc =
    { { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps },
      { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColor2Ops } };
  AttribDataDesc attrDesc =
    { { "a_Vertex", GL_FLOAT, 3, GL_FALSE, vertexOps },
      { "a_OtherVertex", GL_FLOAT, 3, GL_FALSE, vertexOtherOps },
      { "a_Color", GL_UNSIGNED_BYTE, 4, GL_TRUE, colorOps },
      { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, color2Ops },
      { "a_UV", GL_UNSIGNED_BYTE, 1, GL_FALSE } };

  CrossSizeData crossData[] = { { cross_size_arg, false }, { cross_size_arg, true } };

  attrDesc[1].attrOps[0].funcDataConversions.push_back( { CrossVertexConversionTrilines, &crossData[0], "a_Vertex" } );
  attrDesc[1].attrOps[0].funcDataConversions.push_back( { CrossVertexConversionTrilines, &crossData[1], "a_OtherVertex" } );

  unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
  attrDesc[2].default_value = default_color;
  attrDesc[3].default_value = default_color;
  
  AttribDesc *uvdesc = &attrDesc[attrDesc.size()-1];
  uvdesc->repeat_value_length = 6;
  static unsigned char uv_bits[] = { 1, 3, 0, 3, 2, 1 };
  uvdesc->repeat_value = uv_bits;
  
  addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_interpolate"), 0.f);

  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES, VertexBuffer::INTERLEAVED);
}

cgo::draw::shadercylinder2ndcolor::shadercylinder2ndcolor(CGO *I, const float *_origin, 
                                                          const float *_axis, const float _tube_size,
                                                          int _cap, const float *_color2, Pickable *pickcolor2,
                                                          const float _alpha) :
  tube_size(_tube_size), alpha(_alpha) {
  copy3f(_origin, origin);
  copy3f(_axis, axis);
  cap = _cap;
  copy3f(_color2, color2);
  if (pickcolor2){
    I->current_pick_color_index = pick_color_index = pickcolor2->index;
    I->current_pick_color_bond = pick_color_bond = pickcolor2->bond;
  } else {
    pick_color_index = I->current_pick_color_index;
    pick_color_bond = I->current_pick_color_bond;
  }
};

static void SetVertexFromOriginAxisForCylinder(void *varData, const float * pc, void *np, int idx){
  float *varDataF = ((float*)varData);
  add3f(pc, pc + 3, varDataF); // adding origin and axis for both shadercylinder and shadercylinder2ndcolor
}

/**
 * converts all cylinders in the input CGO to a CGO custom operation, which includes picking information (if it exists)
 * 
 * I     - input CGO (includes cylinders)
 * addTo - CGO that vertex_attribute operations are added to (if needed), in this case for caps if values for all
 *         cylinders are the same
 *
 */
CGO *CGOConvertShaderCylindersToCylinderShader(const CGO *I, CGO *addTo){
  /* Lines that pass in two vertices per line */
  PyMOLGlobals *G = I->G;

  // TODO: NEED TO ADD: CGO_CUSTOM_CYLINDER and CGO_CYLINDER
  AttribDataOp vertex1Ops =
    { { CGO_SHADER_CYLINDER,                 1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::shadercylinder, origin), 0 },
      { CGO_SHADER_CYLINDER_WITH_2ND_COLOR,  1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::shadercylinder2ndcolor, origin), 0 },
      { CGO_SAUSAGE,                         1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::sausage, vertex1), 0 },
      { CGO_CYLINDER,                        1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::cylinder, vertex1), 0 },
      { CGO_CUSTOM_CYLINDER,                 1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::custom_cylinder, vertex1), 0 },
      { CGO_CUSTOM_CYLINDER_ALPHA,           1, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::custom_cylinder_alpha, vertex1), 0 } };
  AttribDataOp vertex2Ops =
    { { CGO_SHADER_CYLINDER,                 5, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::shadercylinder, axis), 8 },
      { CGO_SHADER_CYLINDER_WITH_2ND_COLOR,  6, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::shadercylinder2ndcolor, axis), 8 },
      { CGO_SAUSAGE,                         6, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::sausage, vertex2), 8 },
      { CGO_CYLINDER,                        6, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::cylinder, vertex2), 8 },
      { CGO_CUSTOM_CYLINDER,                 6, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::custom_cylinder, vertex2), 8 },
      { CGO_CUSTOM_CYLINDER_ALPHA,           6, FLOAT3_TO_FLOAT3,      offsetof(cgo::draw::custom_cylinder_alpha, vertex2), 8 } };
  static AttribDataOp colorOps =
    { { CGO_COLOR,                           0, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,                           0, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SAUSAGE,                         4, FLOAT3_TO_UB3,         offsetof(cgo::draw::sausage, color1) },
      { CGO_CYLINDER,                        4, FLOAT3_TO_UB3,         offsetof(cgo::draw::cylinder, color1) },
      { CGO_CUSTOM_CYLINDER,                 4, FLOAT3_TO_UB3,         offsetof(cgo::draw::custom_cylinder, color1) },
      { CGO_CUSTOM_CYLINDER_ALPHA,           4, FLOAT4_TO_UB4,         offsetof(cgo::draw::custom_cylinder_alpha, color1) } };
  static AttribDataOp color2Ops =
    { { CGO_COLOR,                           1, FLOAT3_TO_UB3,         0 },
      { CGO_ALPHA,                           1, FLOAT1_TO_UB_4TH,      0 },
      { CGO_SHADER_CYLINDER_WITH_2ND_COLOR,  2, FLOAT3_TO_UB3,         offsetof(cgo::draw::shadercylinder2ndcolor, color2) },
      { CGO_SAUSAGE,                         5, FLOAT3_TO_UB3,         offsetof(cgo::draw::sausage, color2) },
      { CGO_CYLINDER,                        5, FLOAT3_TO_UB3,         offsetof(cgo::draw::cylinder, color2) },
      { CGO_CUSTOM_CYLINDER,                 5, FLOAT3_TO_UB3,         offsetof(cgo::draw::custom_cylinder, color2) },
      { CGO_CUSTOM_CYLINDER_ALPHA,           5, FLOAT4_TO_UB4,         offsetof(cgo::draw::custom_cylinder_alpha, color2) } };
  AttribDataOp radiusOps =
    { { CGO_SHADER_CYLINDER,                 2, FLOAT_TO_FLOAT,        offsetof(cgo::draw::shadercylinder, tube_size), 0 },
      { CGO_SHADER_CYLINDER_WITH_2ND_COLOR,  3, FLOAT_TO_FLOAT,        offsetof(cgo::draw::shadercylinder2ndcolor, tube_size), 0 },
      { CGO_SAUSAGE,                         3, FLOAT_TO_FLOAT,        offsetof(cgo::draw::sausage, radius), 0 },
      { CGO_CYLINDER,                        3, FLOAT_TO_FLOAT,        offsetof(cgo::draw::cylinder, radius), 0 },
      { CGO_CUSTOM_CYLINDER,                 3, FLOAT_TO_FLOAT,        offsetof(cgo::draw::custom_cylinder, radius), 0 },
      { CGO_CUSTOM_CYLINDER_ALPHA,           3, FLOAT_TO_FLOAT,        offsetof(cgo::draw::custom_cylinder_alpha, radius), 0 } };

  AttribDataDesc attrDesc = { { "attr_vertex1", GL_FLOAT,         3, GL_FALSE, vertex1Ops },
                              { "attr_vertex2", GL_FLOAT,         3, GL_FALSE, vertex2Ops },
                              { "a_Color",  GL_UNSIGNED_BYTE, 4, GL_TRUE,  colorOps },
                              { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE,  color2Ops },
                              { "attr_radius",  GL_FLOAT,         1, GL_FALSE, radiusOps } };
  AttribDesc *fdesc;
  static unsigned char cyl_flags[] = { 0, 4, 6, 2, 1, 5, 7, 3 }; // right(4)/up(2)/out(1)

  attrDesc[1].attrOps[0].funcDataConversions.push_back( { SetVertexFromOriginAxisForCylinder, NULL, "attr_vertex2" } );
  attrDesc[1].attrOps[1].funcDataConversions.push_back( { SetVertexFromOriginAxisForCylinder, NULL, "attr_vertex2" } );

  attrDesc.push_back( { "attr_flags", GL_UNSIGNED_BYTE, 1, GL_FALSE } ) ;
  fdesc = &attrDesc[attrDesc.size()-1];
  fdesc->repeat_value = cyl_flags;
  fdesc->repeat_value_length = 8;

  unsigned char default_color[] = { 255, 255, 255, 255 }; // to write in alpha if CGO doesn't have it
  fdesc = &attrDesc[2];
  fdesc->default_value = default_color;
  fdesc = &attrDesc[3];
  fdesc->default_value = default_color;
  float default_radius = 1.f;
  attrDesc[4].default_value = (unsigned char*)&default_radius;

  int box_indices[36] = { // box indices 
    0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 
    0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4 };
  int *box_indices_ptr = NULL;
  box_indices_ptr = box_indices;

  bool interp_same;
  unsigned char cap_value = 0;
  if ((interp_same = CGOCheckShaderCylinderCapInfoIsSame(I, cap_value))){
    addTo->add<cgo::draw::vertex_attribute_1f>(G->ShaderMgr->GetAttributeUID("a_cap"), cap_value );
  } else {
    AttribDataOp interpOps =
      { { CGO_SHADER_CYLINDER,                3, CYL_CAP_TO_CAP,      offsetof(cgo::draw::shadercylinder, cap), 0 },
        { CGO_SHADER_CYLINDER_WITH_2ND_COLOR, 4, CYL_CAP_TO_CAP,      offsetof(cgo::draw::shadercylinder2ndcolor, cap), 0 },
        { CGO_SAUSAGE,                        2, CYL_CAPS_ARE_ROUND,  0, 0 },
        { CGO_CYLINDER,                       2, CYL_CAPS_ARE_FLAT,   0, 0 },
        { CGO_CUSTOM_CYLINDER,                2, CYL_CAPS_ARE_CUSTOM, offsetof(cgo::draw::custom_cylinder, cap1), 0 },
        { CGO_CUSTOM_CYLINDER_ALPHA,          2, CYL_CAPS_ARE_CUSTOM, offsetof(cgo::draw::custom_cylinder_alpha, cap1), 0 },
      };
    attrDesc.push_back({ "a_cap", GL_UNSIGNED_BYTE, 1, GL_FALSE, interpOps } );
  }

  AttribDataOp extraPickColorOps = { { CGO_PICK_COLOR, 1, UINT_INT_TO_PICK_DATA, 0, 0 },
                                     { CGO_SHADER_CYLINDER_WITH_2ND_COLOR,  8, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::shadercylinder2ndcolor, pick_color_index), 0 } };
  AttribDataOp extraPickColor2Ops = { { CGO_PICK_COLOR, 2, UINT_INT_TO_PICK_DATA, 0, 0 },
                                      { CGO_SHADER_CYLINDER_WITH_2ND_COLOR,  5, UINT_INT_TO_PICK_DATA, offsetof(cgo::draw::shadercylinder2ndcolor, pick_color_index), 0 } };
  AttribDataDesc pickDesc = { { "a_Color",  GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColorOps },
                              { "a_Color2", GL_UNSIGNED_BYTE, 4, GL_TRUE, extraPickColor2Ops }};
  return CGOConvertToShader(I, attrDesc, pickDesc, GL_TRIANGLES, VertexBuffer::INTERLEAVED, true, box_indices_ptr, 36);
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
