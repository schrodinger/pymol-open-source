
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
#include"Matrix.h"
#include"ShaderMgr.h"
#include"CoordSet.h"
#include"Rep.h"

#define CLIP_COLOR_VALUE(cv)  ((cv>1.f) ? 255 :  (cv < 0.f) ? 0 : pymol_roundf(cv * 255) )
#define CLIP_NORMAL_VALUE(cv)  ((cv>1.f) ? 127 :  (cv < -1.f) ? -128 : pymol_roundf(((cv + 1.f)/2.f) * 255) - 128 )

struct _CCGORenderer {
  PyMOLGlobals *G;
  RenderInfo *info;
  Rep *rep;
  float alpha;
  short isPicking;
  short use_shader;
  short debug;
  short enable_shaders;
};

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
  register CCGORenderer *I = NULL;

  I = (G->CGORenderer = Calloc(CCGORenderer, 1));
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
  CGO_CYLINDER_SZ,
  CGO_LINEWIDTH_SZ,
  CGO_WIDTHSCALE_SZ,

  CGO_ENABLE_SZ,
  CGO_DISABLE_SZ,
  CGO_SAUSAGE_SZ,
  CGO_CUSTOM_CYLINDER_SZ,

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

#ifdef _PYMOL_CGO_DRAWARRAYS
  CGO_DRAW_ARRAYS_SZ,
#else
  CGO_NULL_SZ,
#endif
  CGO_NULL_SZ,
  CGO_RESET_NORMAL_SZ,
  CGO_PICK_COLOR_SZ,

#ifdef _PYMOL_CGO_DRAWBUFFERS
  CGO_DRAW_BUFFERS_SZ,  
  CGO_DRAW_BUFFERS_INDEXED_SZ,  
  CGO_BOUNDING_BOX_SZ,
  CGO_DRAW_BUFFERS_NOT_INDEXED_SZ,  
  CGO_LINEWIDTH_SPECIAL_SZ,
#else
  CGO_NULL_SZ, CGO_NULL_SZ, CGO_NULL_SZ, CGO_NULL_SZ, CGO_NULL_SZ,
#endif
  CGO_DRAW_CYLINDER_BUFFERS_SZ,
  CGO_SHADER_CYLINDER_SZ,
  CGO_SHADER_CYLINDER_WITH_2ND_COLOR_SZ,
  CGO_DRAW_SPHERE_BUFFERS_SZ,  
  CGO_ACCESSIBILITY_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,
  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,
  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ,  CGO_NULL_SZ
};

typedef void CGO_op(CCGORenderer * I, float **);
typedef CGO_op *CGO_op_fn;

static float *CGO_add(CGO * I, int c);
static float *CGO_size(CGO * I, int sz);
static void subdivide(int n, float *x, float *y);
static void CGOSimpleCylinder(CGO * I, float *v1, float *v2, float tube_size, float *c1,
                              float *c2, int cap1, int cap2);
static void CGOSimpleEllipsoid(CGO * I, float *v, float vdw, float *n0, float *n1,
                               float *n2);
static void CGOSimpleQuadric(CGO * I, float *v, float vdw, float *q);
static void CGOSimpleSphere(CGO * I, float *v, float vdw);
static void CGOSimpleCone(CGO * I, float *v1, float *v2, float r1, float r2, float *c1,
                          float *c2, int cap1, int cap2);

CGO *CGOProcessShape(CGO * I, struct GadgetSet *gs, CGO * result)
{
  register float *p, *pc = I->op;
  register float *n, *nc;
  register int op;
  int sz;
  if(!result)
    result = CGONew(I->G);
  CGOReset(result);
  VLACheck(result->op, float, I->c + 32);

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    sz = CGO_sz[op];
    nc = CGO_add(result, sz + 1);
    *(nc++) = *(pc - 1);
    switch (op) {
    case CGO_NORMAL:
      GadgetSetFetchNormal(gs, pc, nc);
      break;
    case CGO_COLOR:
      GadgetSetFetchColor(gs, pc, nc);
      break;
    case CGO_VERTEX:
      GadgetSetFetch(gs, pc, nc);
#ifdef _PYMOL_CGO_DRAWARRAYS
      result->has_begin_end = true;      
#endif
      break;
    case CGO_FONT_VERTEX:
      GadgetSetFetch(gs, pc, nc);
      break;
    case CGO_SPHERE:
      GadgetSetFetch(gs, pc, nc);
      *(nc + 3) = *(pc + 3);
      break;
    case CGO_CUSTOM_CYLINDER:
      GadgetSetFetch(gs, pc, nc);
      GadgetSetFetch(gs, pc + 3, nc + 3);
      GadgetSetFetchColor(gs, pc + 7, nc + 7);
      GadgetSetFetchColor(gs, pc + 10, nc + 10);
      *(nc + 6) = *(pc + 6);
      *(nc + 13) = *(pc + 13);
      *(nc + 14) = *(pc + 14);
      break;
    case CGO_CYLINDER:
      GadgetSetFetch(gs, pc, nc);
      GadgetSetFetch(gs, pc + 3, nc + 3);
      GadgetSetFetchColor(gs, pc + 7, nc + 7);
      GadgetSetFetchColor(gs, pc + 10, nc + 10);
      *(nc + 6) = *(pc + 6);
      break;
    case CGO_SAUSAGE:
      GadgetSetFetch(gs, pc, nc);
      GadgetSetFetch(gs, pc + 3, nc + 3);
      GadgetSetFetchColor(gs, pc + 7, nc + 7);
      GadgetSetFetchColor(gs, pc + 10, nc + 10);
      *(nc + 6) = *(pc + 6);
      break;
    case CGO_TRIANGLE:
      GadgetSetFetch(gs, pc, nc);
      GadgetSetFetch(gs, pc + 3, nc + 3);
      GadgetSetFetch(gs, pc + 6, nc + 6);
      GadgetSetFetchNormal(gs, pc + 9, nc + 9);
      GadgetSetFetchNormal(gs, pc + 12, nc + 12);
      GadgetSetFetchNormal(gs, pc + 15, nc + 15);
      GadgetSetFetchColor(gs, pc + 18, nc + 18);
      GadgetSetFetchColor(gs, pc + 21, nc + 21);
      GadgetSetFetchColor(gs, pc + 24, nc + 24);
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int mode = CGO_read_int(pc), arrays = CGO_read_int(pc), narrays = CGO_read_int(pc), nverts = CGO_read_int(pc), i, tsz;
	sz = narrays*nverts;
	tsz = sz + 4;
	nc = CGO_add(result, tsz);
	CGO_write_int(nc, mode); 
	CGO_write_int(nc, arrays);
	CGO_write_int(nc, narrays); 
	CGO_write_int(nc, nverts);
	if (arrays & CGO_VERTEX_ARRAY){
	  for (i=0; i<nverts; i++){
	    GadgetSetFetch(gs, pc, nc);
	    pc += 3; nc += 3;
	  }
	}
	if (arrays & CGO_NORMAL_ARRAY){
	  for (i=0; i<nverts; i++){
	    GadgetSetFetchNormal(gs, pc, nc);
	    pc += 3; nc += 3;
	  }
	}
	if (arrays & CGO_COLOR_ARRAY){
	  for (i=0; i<nverts; i++){
	    GadgetSetFetchColor(gs, pc, nc);
	    pc += 4; nc += 4;
	  }
	}
	if (arrays & CGO_PICK_COLOR_ARRAY){
	  pc += nverts * 3;
	  nc += nverts * 3;
	}
      }
      break;
    case CGO_BEGIN:
    case CGO_END:
      result->has_begin_end = true;
#endif
    default:
      p = pc;
      n = nc;
      while(sz--)
        *(n++) = *(p++);
      break;
    }
    pc += CGO_sz[op];
    nc += CGO_sz[op];
  }
  CGOStop(result);
  return (result);
}

#ifndef _PYMOL_NOPY

static PyObject *CGOArrayAsPyList(CGO * I)
{

  register float *pc = I->op;
  register int op;
  int i;
  int cc;
  PyObject *result = NULL;

  result = PyList_New(I->c);

  i = 0;
  if(I->c) {
    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      PyList_SetItem(result, i++, PyFloat_FromDouble((float) op));
      cc = CGO_sz[op];
      switch (op) {             /* now convert any instructions with int arguments */
      case CGO_BEGIN:
      case CGO_ENABLE:
      case CGO_DISABLE:
        PyList_SetItem(result, i++, PyFloat_FromDouble((float) CGO_read_int(pc)));
        cc--;
        break;
      case CGO_DRAW_ARRAYS:
	{
	  int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	  cc = floatlength + 4 ;
	  PyList_SetItem(result, i++, PyFloat_FromDouble((float) CGO_read_int(pc)));
	  PyList_SetItem(result, i++, PyFloat_FromDouble((float) CGO_read_int(pc)));
	  PyList_SetItem(result, i++, PyFloat_FromDouble((float) CGO_read_int(pc)));
	  PyList_SetItem(result, i++, PyFloat_FromDouble((float) CGO_read_int(pc)));
	}
      }
     if(cc > 0)
        while(cc--) {
          PyList_SetItem(result, i++, PyFloat_FromDouble(*(pc++)));
        }
    }
  }
  while(i < I->c) {
    PyList_SetItem(result, i++, PyFloat_FromDouble((float) CGO_STOP));
  }
  return (result);
}
#endif

#ifndef _PYMOL_NOPY

static int CGOArrayFromPyListInPlace(PyObject * list, CGO * I)
{
  int a;
  int c = I->c;
  int cc = 0;
  int ok = true;
  float *pc;
  int sz, op;
  int l;
  if(!list) {
    ok = false;
  } else if(!PyList_Check(list)) {
    ok = false;
  } else {
    l = PyList_Size(list);
    if(l != I->c)
      ok = false;
  }
  if(ok) {
    pc = I->op;

    while(c > 0) {
      op = (int) PyFloat_AsDouble(PyList_GetItem(list, cc++));
      op = op & CGO_MASK;
      c--;
      sz = CGO_sz[op];
      CGO_write_int(pc, op);
      ok = true;
#ifdef _PYMOL_CGO_DRAWARRAYS
      switch (op) {
      case CGO_END:
      case CGO_VERTEX:
      case CGO_BEGIN:
	I->has_begin_end = true;
      }
#endif
      switch (op) {             /* now convert any instructions with int arguments */
      case CGO_BEGIN:
      case CGO_ENABLE:
      case CGO_DISABLE:
        CGO_write_int(pc, (int) PyFloat_AsDouble(PyList_GetItem(list, cc++)));
        c--;
        sz--;
        break;
      case CGO_DRAW_ARRAYS:
	{
	  int arrays, narrays, nverts;	  
	  CGO_write_int(pc, (int) PyFloat_AsDouble(PyList_GetItem(list, cc++)));
	  CGO_write_int(pc, arrays = (int) PyFloat_AsDouble(PyList_GetItem(list, cc++)));
	  CGO_write_int(pc, narrays = (int) PyFloat_AsDouble(PyList_GetItem(list, cc++)));
	  CGO_write_int(pc, nverts = (int) PyFloat_AsDouble(PyList_GetItem(list, cc++)));	
	  c -= 4;
	  sz = narrays*nverts;
	}
	break;
      }

      for(a = 0; a < sz; a++) {
        *(pc++) = (float) PyFloat_AsDouble(PyList_GetItem(list, cc++));
        c--;
      }
    }
  }
  return (ok);
}
#endif

PyObject *CGOAsPyList(CGO * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result;
  result = PyList_New(2);
  PyList_SetItem(result, 0, PyInt_FromLong(I->c));
  PyList_SetItem(result, 1, CGOArrayAsPyList(I));
  return (result);
#endif
}

CGO *CGONewFromPyList(PyMOLGlobals * G, PyObject * list, int version)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  int ok = true;
  int ll;
  OOCalloc(G, CGO);
  I->G = G;
  I->op = NULL;
  I->i_start = 0;
  I->debug = 0;
#ifdef _PYMOL_CGO_DRAWARRAYS
  I->has_begin_end = false;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->has_draw_sphere_buffers = false;
#endif
  I->enable_shaders = 0;
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO ENABLE BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->c);
  if(ok){
    ok = ((I->op = VLAlloc(float, I->c + 1)) != NULL);
  }
  if((version > 0) && (version <= 86)) {
    if(ok)
      ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 1), I->op, I->c);
  } else {
    if(ok)
      ok = CGOArrayFromPyListInPlace(PyList_GetItem(list, 1), I);
  }
  if(!ok) {
    CGOFree(I);
    I = NULL;
  }
#ifdef _PYMOL_CGO_DRAWARRAYS
  {
    CGO *cgo = NULL;
    if (I && I->has_begin_end){
      cgo = CGOCombineBeginEnd(I, 0);
      CGOFree(I);
    } else {
      cgo = I;
    }
    return cgo;
  }
#else
  return (I);
#endif
#endif
}

CGO *CGONew(PyMOLGlobals * G)
{
  OOCalloc(G, CGO);
  I->G = G;
  I->op = VLAlloc(float, 33);
  I->i_start = 0;
  I->alpha = 1.f;
#ifdef _PYMOL_CGO_DRAWARRAYS
  I->has_begin_end = false;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->normal[0] = 0.f; I->normal[1] = 0.f; I->normal[2] = 1.f;
  I->color[0] = 0.f; I->color[1] = 0.f; I->color[2] = 1.f;
  I->pickColor[0] = 0; I->pickColor[1] = 0; I->pickColor[2] = 0; I->pickColor[3] = 255;
  I->current_accessibility = 1.f;
#endif
  I->enable_shaders = 0;
  return (I);
}

CGO *CGONewSized(PyMOLGlobals * G, int size)
{
  OOCalloc(G, CGO);
  I->G = G;
  I->op = VLAlloc(float, size + 32);
  I->i_start = 0;
  I->alpha = 1.f;
#ifdef _PYMOL_CGO_DRAWARRAYS
  I->has_begin_end = false;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->normal[0] = 0.f; I->normal[1] = 0.f; I->normal[2] = 1.f;
  I->color[0] = 0.f; I->color[1] = 0.f; I->color[2] = 1.f;
  I->pickColor[0] = 0; I->pickColor[1] = 0; I->pickColor[2] = 0; I->pickColor[3] = 255;
  I->current_accessibility = 1.f;
#endif
  I->enable_shaders = 0;
  return (I);
}

void CGOReset(CGO * I)
{
  I->c = 0;
  I->z_flag = false;
  I->alpha = 1.f;
#ifdef _PYMOL_CGO_DRAWARRAYS
  I->has_begin_end = false;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->has_draw_buffers = false;
  I->has_draw_cylinder_buffers = false;
  I->normal[0] = 0.f; I->normal[1] = 0.f; I->normal[2] = 1.f;
  I->color[0] = 0.f; I->color[1] = 0.f; I->color[2] = 1.f;
  I->pickColor[0] = 0; I->pickColor[1] = 0; I->pickColor[2] = 0; I->pickColor[3] = 255;
  I->current_accessibility = 1.f;
#endif
}

void CGOFreeWithoutVBOs(CGO * I){
  CGOFreeImpl(I, 0);
}

void CGOFree(CGO * I){
  CGOFreeImpl(I, 1);
}

void CGOFreeImpl(CGO * I, short withVBOs)
{
  if(I) {
#ifdef _PYMOL_CGO_DRAWBUFFERS
    if (withVBOs && I->has_draw_buffers){
      CGOFreeVBOs(I);
    }
#endif
    if(I->i_start) {
      FreeP(I->i_start);
    }
    VLAFreeP(I->op);
  }
  OOFreeP(I);
}

static float *CGO_add(CGO * I, int c)
{
  float *at;
  VLACheck(I->op, float, I->c + c);
  at = I->op + I->c;
  I->c += c;
  return (at);
}

float *CGO_add_GLfloat(CGO * I, int c)
{
  GLfloat *at;
  VLACheck(I->op, GLfloat, I->c + c);
  at = I->op + I->c;
  I->c += c;
  return (at);
}

static float *CGO_size(CGO * I, int sz)
{
  float *at;
  VLASize(I->op, float, sz);
  at = I->op + I->c;
  I->c = sz;
  return (at);
}


/*===== Object Creation Routines =======*/

int CGOFromFloatArray(CGO * I, float *src, int len)
{
  int op, iarg;
  int c;
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
    c = 1;
    op = CGO_MASK & ((int) (*(src++)));
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
#ifdef _PYMOL_CGO_DRAWARRAYS
      switch (op) {
      case CGO_END:
      case CGO_VERTEX:
      case CGO_BEGIN:
	I->has_begin_end = true;
      }
#endif
      switch (op) {             /* now convert any instructions with int arguments */
      case CGO_BEGIN:
      case CGO_ENABLE:
      case CGO_DISABLE:
        tf = save_pc + 1;
        iarg = (int) *(tf);
        CGO_write_int(tf, iarg);
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

void CGOBegin(CGO * I, int mode)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_BEGIN);
  CGO_write_int(pc, mode);
#ifdef _PYMOL_CGO_DRAWARRAYS
  I->has_begin_end = true;
#endif
  /*#ifdef _PYMOL_CGO_DRAWARRAYS
  PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOBegin() is called but not implemented in OpenGLES\n" ENDFB(I->G);
  #endif*/
}

void CGOEnd(CGO * I)
{
  float *pc = CGO_add(I, 1);
  CGO_write_int(pc, CGO_END);
#ifdef _PYMOL_CGO_DRAWARRAYS
  I->has_begin_end = true;
#endif
  /*#ifdef _PYMOL_CGO_DRAWARRAYS
    if (!override){
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOEnd() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    }
    #endif*/
}

void CGOEnable(CGO * I, int mode)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_ENABLE);
  CGO_write_int(pc, mode);
}

void CGODisable(CGO * I, int mode)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_DISABLE);
  CGO_write_int(pc, mode);
}

void CGOLinewidth(CGO * I, float v)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_LINEWIDTH);
  *(pc++) = v;
}

void CGOLinewidthSpecial(CGO * I, int v)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_LINEWIDTH_SPECIAL);
  CGO_write_int(pc, v);
}

void CGODotwidth(CGO * I, float v)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_DOTWIDTH);
  *(pc++) = v;
}

GLfloat *CGODrawArrays(CGO *I, GLenum mode, short arrays, int nverts){
#ifdef _PYMOL_CGO_DRAWARRAYS
  int narrays = 0, floatlength;
  short bit;
  float *pc = CGO_add(I, 5);
  CGO_write_int(pc, CGO_DRAW_ARRAYS);
  CGO_write_int(pc, mode);
  CGO_write_int(pc, arrays);
  for (bit = 0; bit < 4; bit++){
    if ((1 << bit) & arrays){
      narrays+=3;
    }
  }
  if (arrays & CGO_ACCESSIBILITY_ARRAY) narrays++;
  if (arrays & CGO_COLOR_ARRAY) narrays++; //floatlength += nverts;
  CGO_write_int(pc, narrays);
  CGO_write_int(pc, nverts);
  floatlength = narrays*nverts;
  return CGO_add_GLfloat(I, floatlength); /* 4 floats for color array, 3 floats for every other array, 1 float for accessibility  */
#else
  return NULL;
#endif
}

#ifdef _PYMOL_CGO_DRAWBUFFERS
void CGODrawBuffers(CGO *I, GLenum mode, short arrays, int nverts, uint *bufs){
  int narrays = 0;
  short bit;

  float *pc = CGO_add(I, 9);
  CGO_write_int(pc, CGO_DRAW_BUFFERS);
  CGO_write_int(pc, mode);
  CGO_write_int(pc, arrays);
  for (bit = 0; bit < 4; bit++){
    if ((1 << bit) & arrays){
      narrays+=3;
    }
  }
  if (arrays & CGO_ACCESSIBILITY_ARRAY) narrays++;
  if (arrays & CGO_COLOR_ARRAY) narrays++; //floatlength += nverts;
  CGO_write_int(pc, narrays);
  CGO_write_int(pc, nverts);
  for (bit = 0; bit<4; bit++){
    CGO_write_int(pc, bufs[bit]);
  }
  return;
}
void CGOBoundingBox(CGO *I, float *min, float *max){
  float *pc = CGO_add(I, 7);
  CGO_write_int(pc, CGO_BOUNDING_BOX);
  *(pc++) = *(min);
  *(pc++) = *(min+1);
  *(pc++) = *(min+2);
  *(pc++) = *(max);
  *(pc++) = *(max+1);
  *(pc++) = *(max+2);
}

void CGOAccessibility(CGO * I, float a)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_ACCESSIBILITY);
  *(pc++) = a;
}
GLfloat *CGODrawBuffersIndexed(CGO *I, GLenum mode, short arrays, int nindices, int nverts, uint *bufs){
  int narrays = 0;
  short bit;

  float *pc = CGO_add(I, 11);
  CGO_write_int(pc, CGO_DRAW_BUFFERS_INDEXED);
  CGO_write_int(pc, mode);
  CGO_write_int(pc, arrays);
  for (bit = 0; bit < 4; bit++){
    if ((1 << bit) & arrays){
      narrays++;
    }
  }
  if (arrays & CGO_ACCESSIBILITY_ARRAY) narrays++;
  if (arrays & CGO_COLOR_ARRAY) narrays++; //floatlength += nverts;
  CGO_write_int(pc, narrays);
  CGO_write_int(pc, nindices);
  CGO_write_int(pc, nverts);
  for (bit = 0; bit<5; bit++){
    CGO_write_int(pc, bufs[bit]);
  }
  return CGO_add_GLfloat(I, nverts*3); /* 3 floats for pick array, 1st 1/3 for color, 2/3 for atom/bond info */
}

GLfloat *CGODrawBuffersNotIndexed(CGO *I, GLenum mode, short arrays, int nverts, uint *bufs){
  int narrays = 0;
  short bit;

  float *pc = CGO_add(I, 9);
  CGO_write_int(pc, CGO_DRAW_BUFFERS_NOT_INDEXED);
  CGO_write_int(pc, mode);
  CGO_write_int(pc, arrays);
  for (bit = 0; bit < 4; bit++){
    if ((1 << bit) & arrays){
      narrays++;
    }
  }
  if (arrays & CGO_ACCESSIBILITY_ARRAY) narrays++;
  if (arrays & CGO_COLOR_ARRAY) narrays++; //floatlength += nverts;
  CGO_write_int(pc, narrays);
  CGO_write_int(pc, nverts);
  for (bit = 0; bit<4; bit++){
    CGO_write_int(pc, bufs[bit]);
  }
  return CGO_add_GLfloat(I, nverts*3); /* 3 floats for pick array, 1st 1/3 for color, 2/3 for atom/bond info */
}

void CGOShaderCylinder(CGO *I, float *origin, float *axis, float tube_size, int cap){
  float *pc = CGO_add(I, 9);
  CGO_write_int(pc, CGO_SHADER_CYLINDER);
  *(pc++) = *(origin);
  *(pc++) = *(origin+1);
  *(pc++) = *(origin+2);
  *(pc++) = *(axis);
  *(pc++) = *(axis+1);
  *(pc++) = *(axis+2);
  *(pc++) = tube_size;
  CGO_write_int(pc, cap);
}
void CGOShaderCylinder2ndColor(CGO *I, float *origin, float *axis, float tube_size, int cap, float *color2){
  float *pc = CGO_add(I, 12);
  CGO_write_int(pc, CGO_SHADER_CYLINDER_WITH_2ND_COLOR);
  *(pc++) = *(origin);
  *(pc++) = *(origin+1);
  *(pc++) = *(origin+2);
  *(pc++) = *(axis);
  *(pc++) = *(axis+1);
  *(pc++) = *(axis+2);
  *(pc++) = tube_size;
  CGO_write_int(pc, cap);
  *(pc++) = *(color2);
  *(pc++) = *(color2+1);
  *(pc++) = *(color2+2);
}

void CGODrawCylinderBuffers(CGO *I, int num_cyl, int alpha, uint *bufs){
  short bit;

  float *pc = CGO_add(I, 8);
  CGO_write_int(pc, CGO_DRAW_CYLINDER_BUFFERS);
  CGO_write_int(pc, num_cyl);
  CGO_write_int(pc, alpha);
  for (bit = 0; bit<5; bit++){
    CGO_write_int(pc, bufs[bit]);
  }
}

void CGODrawSphereBuffers(CGO *I, int num_spheres, int ub_flags, uint *bufs){
  /* ub_flags are for whether the color and flags buffer are UB:
     1 - color buffer is UB
     2 - flags buffer is UB 
  */
  short bit;
  float *pc = CGO_add(I, 6);
  CGO_write_int(pc, CGO_DRAW_SPHERE_BUFFERS);
  CGO_write_int(pc, num_spheres);
  CGO_write_int(pc, ub_flags);
  for (bit = 0; bit<3; bit++){
    CGO_write_int(pc, bufs[bit]);
  }
}

#endif


void CGOCylinderv(CGO * I, float *p1, float *p2, float r, float *c1, float *c2)
{
  float *pc = CGO_add(I, 14);
  CGO_write_int(pc, CGO_CYLINDER);
  *(pc++) = *(p1++);
  *(pc++) = *(p1++);
  *(pc++) = *(p1++);
  *(pc++) = *(p2++);
  *(pc++) = *(p2++);
  *(pc++) = *(p2++);
  *(pc++) = r;
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
}

void CGOCustomCylinderv(CGO * I, float *p1, float *p2, float r, float *c1, float *c2,
                        float cap1, float cap2)
{
  float *pc = CGO_add(I, 16);
  CGO_write_int(pc, CGO_CUSTOM_CYLINDER);
  *(pc++) = *(p1++);
  *(pc++) = *(p1++);
  *(pc++) = *(p1++);
  *(pc++) = *(p2++);
  *(pc++) = *(p2++);
  *(pc++) = *(p2++);
  *(pc++) = r;
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
  *(pc++) = cap1;
  *(pc++) = cap2;
}

void CGOConev(CGO * I, float *p1, float *p2, float r1, float r2, float *c1, float *c2,
              float cap1, float cap2)
{
  float *pc = CGO_add(I, 17);
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
}

void SetCGOPickColor(float *begPickColorVals, int nverts, int pl, int index, int bond){
  float *colorVals = begPickColorVals + nverts;
  int pld = pl / 3;
  /* Picking Color stores atom/bond info in the second 2/3rds of the array */
  CGO_put_int(colorVals + (pld*2), index);
  CGO_put_int(colorVals + (pld*2) + 1, bond);
#ifdef _PYMOL_CGO_DRAWBUFFERS
  /* TODO: Should set I->pickColor but do not have I */
#endif
}
void CGOPickColor(CGO * I, int index, int bond)
{
  float *pc = CGO_add(I, 3);
  CGO_write_int(pc, CGO_PICK_COLOR);
  CGO_write_int(pc, index);
  CGO_write_int(pc, bond);
  I->current_pick_color_index = index;
  I->current_pick_color_bond = bond;
}

void CGOAlpha(CGO * I, float alpha)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_ALPHA);
  *(pc++) = alpha;
  I->alpha = alpha;
}

void CGOSphere(CGO * I, float *v1, float r)
{
  float *pc = CGO_add(I, 5);
  CGO_write_int(pc, CGO_SPHERE);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = r;
}

void CGOEllipsoid(CGO * I, float *v1, float r, float *n1, float *n2, float *n3)
{
  float *pc = CGO_add(I, 14);
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
}

void CGOQuadric(CGO * I, float *v, float r, float *q)
{
  float *pc = CGO_add(I, 15);
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

}

void CGOSausage(CGO * I, float *v1, float *v2, float r, float *c1, float *c2)
{
  float *pc = CGO_add(I, 14);
  CGO_write_int(pc, CGO_SAUSAGE);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v1++);
  *(pc++) = *(v2++);
  *(pc++) = *(v2++);
  *(pc++) = *(v2++);
  *(pc++) = r;
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c1++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
  *(pc++) = *(c2++);
}

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

void CGOAlphaTriangle(CGO * I,
                      float *v1, float *v2, float *v3,
                      float *n1, float *n2, float *n3,
                      float *c1, float *c2, float *c3,
                      float a1, float a2, float a3, int reverse)
{
  if(v1 && v2 && v3) {
    float *pc = CGO_add(I, CGO_ALPHA_TRIANGLE_SZ + 1);
    register float z = _0;

    CGO_write_int(pc, CGO_ALPHA_TRIANGLE);
    CGO_write_int(pc, 0);
    *(pc++) = (v1[0] + v2[0] + v3[0]) * one_third;
    *(pc++) = (v1[1] + v2[1] + v3[1]) * one_third;
    *(pc++) = (v1[2] + v2[2] + v3[2]) * one_third;
    if(I->z_flag) {
      register float *zv = I->z_vector;
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
}

void CGOVertex(CGO * I, float v1, float v2, float v3)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
}

void CGOVertexv(CGO * I, float *v)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
}

void CGOColor(CGO * I, float v1, float v2, float v3)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_COLOR);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->color[0] = v1;
  I->color[1] = v2;
  I->color[2] = v3;
#endif
}

void CGOColorv(CGO * I, float *v)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_COLOR);
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->color[0] = v[0];
  I->color[1] = v[1];
  I->color[2] = v[2];
#endif
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
}

void CGONormal(CGO * I, float v1, float v2, float v3)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_NORMAL);
  *(pc++) = v1;
  *(pc++) = v2;
  *(pc++) = v3;
#ifdef _PYMOL_CGO_DRAWBUFFERS
  I->normal[0] = v1;
  I->normal[1] = v2;
  I->normal[2] = v3;
#endif
}

void CGOResetNormal(CGO * I, int mode)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_RESET_NORMAL);
  CGO_write_int(pc, mode);
#ifdef _PYMOL_CGO_DRAWBUFFERS
  SceneGetResetNormal(I->G, I->normal, mode);
#endif
}

void CGOFontVertexv(CGO * I, float *v)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_FONT_VERTEX);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
}

void CGOFontVertex(CGO * I, float x, float y, float z)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_FONT_VERTEX);
  *(pc++) = x;
  *(pc++) = y;
  *(pc++) = z;
}

void CGOFontScale(CGO * I, float v1, float v2)
{
  float *pc = CGO_add(I, 3);
  CGO_write_int(pc, CGO_FONT_SCALE);
  *(pc++) = v1;
  *(pc++) = v2;
}

void CGOChar(CGO * I, char c)
{
  float *pc = CGO_add(I, 2);
  CGO_write_int(pc, CGO_CHAR);
  *(pc++) = (float) c;
}

void CGOIndent(CGO * I, char c, float dir)
{
  float *pc = CGO_add(I, 3);
  CGO_write_int(pc, CGO_INDENT);
  *(pc++) = (float) c;
  *(pc++) = dir;
}

void CGOWrite(CGO * I, char *str)
{
  float *pc;

  while(*str) {
    pc = CGO_add(I, 2);
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(str++);
  }
}

void CGOWriteLeft(CGO * I, char *str)
{
  float *pc;
  char *s;
  s = str;
  while(*s) {
    pc = CGO_add(I, 3);
    CGO_write_int(pc, CGO_INDENT);
    *(pc++) = (float) *(s++);
    *(pc++) = -1.0F;
  }
  s = str;
  while(*s) {
    pc = CGO_add(I, 2);
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(s++);
  }
}

void CGOWriteIndent(CGO * I, char *str, float indent)
{
  float *pc;
  char *s;
  s = str;
  while(*s) {
    pc = CGO_add(I, 3);
    CGO_write_int(pc, CGO_INDENT);
    *(pc++) = (float) *(s++);
    *(pc++) = indent;
  }
  s = str;
  while(*s) {
    pc = CGO_add(I, 2);
    CGO_write_int(pc, CGO_CHAR);
    *(pc++) = (float) *(s++);
  }
}

void CGONormalv(CGO * I, float *v)
{
  float *pc = CGO_add(I, 4);
  CGO_write_int(pc, CGO_NORMAL);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
  *(pc++) = *(v++);
}

void CGOStop(CGO * I)
{
  /* add enough zeros to prevent overrun in the event of corruption
   * (include more zeros than the longest instruction in the compiler 

   * although this is wasteful, it does prevent crashes...
   */

#define CGO_STOP_ZEROS 16

  float *pc = CGO_size(I, I->c + CGO_STOP_ZEROS);

  UtilZeroMem(pc, sizeof(float) * CGO_STOP_ZEROS);
}

int CGOCheckComplex(CGO * I)
{
  register float *pc = I->op;
  int fc = 0;
  int nEdge;
  int op;
  SphereRec *sp;

  sp = I->G->Sphere->Sphere[1];

  /* stick_quality needs to match *every* CGO? */
  nEdge = (int) SettingGet(I->G, cSetting_stick_quality);

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    switch (op) {
    case CGO_CYLINDER:
    case CGO_CONE:
    case CGO_SAUSAGE:
    case CGO_CUSTOM_CYLINDER:
      fc += 3 * (3 + (nEdge + 1) * 9) + 9;
      break;
    case CGO_ELLIPSOID:
    case CGO_QUADRIC:
    case CGO_SPHERE:
      fc += (sp->NVertTot * 6) + (sp->NStrip * 3) + 3;
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	fc += nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int mode = CGO_get_int(pc), nverts = CGO_get_int(pc + 4), nindices = CGO_get_int(pc + 3);
	switch(mode){
	case GL_TRIANGLES:
	  fc += nindices / 3;	  
	  break;
	case GL_LINES:
	  fc += nindices / 2;	  
	  break;
	}
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int mode = CGO_get_int(pc), nverts = CGO_get_int(pc + 3);
	switch(mode){
	case GL_TRIANGLES:
	  fc += nverts / 3;	  
	  break;
	case GL_LINES:
	  fc += nverts / 2;	  
	  break;
	}
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (fc);
}

int CGOPreloadFonts(CGO * I)
{
  int ok = true;
  register float *pc = I->op;
  int op;
  int font_seen = false;
  int font_id;
  int blocked = false;

  blocked = PAutoBlock(I->G);
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
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
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  if(blocked)
    PUnblock(I->G);
  return (ok);
}

int CGOCheckForText(CGO * I)
{
  register float *pc = I->op;
  int fc = 0;
  int op;

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
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
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  PRINTFD(I->G, FB_CGO)
    " CGOCheckForText-Debug: %d\n", fc ENDFD;

  return (fc);
}

CGO *CGODrawText(CGO * I, int est, float *camera)
{                               /* assumes blocked intepreter */
  CGO *cgo;

  register float *pc = I->op;
  register float *nc;
  register int op;
  float *save_pc;
  int sz;
  int font_id = 0;
  char text[2] = " ";
  float pos[] = { 0.0F, 0.0F, 0.0F };
  float axes[] = { 1.0F, 0.0F, 0.0F,
    0.0F, 1.0F, 0.0F,
    0.0F, 0.0F, 1.0F
  };
  float scale[2] = { 1.0, 1.0 };

  cgo = CGONewSized(I->G, I->c + est);

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
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
      VFontWriteToCGO(I->G, font_id, cgo, text, pos, scale, axes);
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), tsz;
	sz = narrays*nverts + 4 ;
	tsz = sz;
	nc = CGO_add(cgo, sz + 1);
	*(nc++) = *(pc - 1);
	while(sz--)
	  *(nc++) = *(pc++);
	save_pc += tsz ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4), tsz;
	sz = nverts*3 + 10 ;
	tsz = sz;
	nc = CGO_add(cgo, sz + 1);
	*(nc++) = *(pc - 1);
	while(sz--)
	  *(nc++) = *(pc++);
	save_pc += tsz ;
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3), tsz;
	sz = nverts*3 + 8 ;
	tsz = sz;
	nc = CGO_add(cgo, sz + 1);
	*(nc++) = *(pc - 1);
	while(sz--)
	  *(nc++) = *(pc++);
	save_pc += tsz ;
      }
      break;
#endif
    default:
      sz = CGO_sz[op];
      nc = CGO_add(cgo, sz + 1);
      *(nc++) = *(pc - 1);
      while(sz--)
        *(nc++) = *(pc++);
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
  CGOStop(cgo);
#ifdef _PYMOL_CGO_DRAWARRAYS
  if (cgo && cgo->has_begin_end){
    /* this is mainly for VFontWriteToCGO() that still creates CGOBegin/CGOEnd */
    if(cgo && cgo->has_begin_end){
      CGO *convertcgo = NULL;
      convertcgo = CGOCombineBeginEnd(cgo, 0);
      CGOFree(cgo);
      cgo = convertcgo;
    }
  }
#endif
  return (cgo);
}

CGO *CGOCombineBeginEnd(CGO * I, int est)
{
#ifdef _PYMOL_CGO_DRAWARRAYS
  CGO *cgo;

  register float *pc = I->op;
  register float *nc;
  register int op;
  float *save_pc;
  int sz;

  cgo = CGONewSized(I->G, I->c + est);
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    switch (op) {
    case CGO_DRAW_ARRAYS:
      {
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	GLfloat *vals = CGODrawArrays(cgo, mode, arrays, nverts);
	int nvals = narrays*nverts, onvals ;
	onvals = nvals;
	pc += 4;
	while(nvals--)
	  *(vals++) = *(pc++);
	save_pc += onvals + 4 ;
      }
      break;
    case CGO_END:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOCombineBeginEnd: CGO_END encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOCombineBeginEnd: CGO_VERTEX encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_BEGIN:
      {
	float *origpc = pc;
	int nverts = 0, damode = CGO_VERTEX_ARRAY, err = 0, end = 0;
	int mode = CGO_read_int(pc);

	while(!err && !end && (op = (CGO_MASK & CGO_read_int(pc)))) {
	  switch (op) {
	  case CGO_DRAW_ARRAYS:
	    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOCombineBeginEnd: CGO_DRAW_ARRAYS encountered inside CGO_BEGIN/CGO_END\n" ENDFB(I->G);
	    err = true;
	    continue;
	  case CGO_NORMAL:
	    damode |= CGO_NORMAL_ARRAY;
	    break;
	  case CGO_COLOR:
	    damode |= CGO_COLOR_ARRAY;
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
	  case CGO_END:
	    end = 1;
	  default:
	    break;
	  }
	  sz = CGO_sz[op];
	  pc += sz;
	}
	if (nverts>0 && !err){
	  int pl = 0, plc = 0, plp = 0, pla = 0;
	  float *vertexVals, *tmp_ptr;
	  float *normalVals = 0, *colorVals = 0, *nxtVals = 0, *pickColorVals = 0, *accessibilityVals = 0;
	  uchar *pickColorValsUC;
	  short notHaveValue = 0, nxtn = 3;
	  nxtVals = vertexVals = CGODrawArrays(cgo, mode, damode, nverts);	      
	  if (damode & CGO_NORMAL_ARRAY){
	    nxtVals = normalVals = vertexVals + (nxtn*nverts);
	  }
	  if (damode & CGO_COLOR_ARRAY){
	    nxtVals = colorVals = nxtVals + (nxtn*nverts);
	    nxtn = 4;
	  }
	  if (damode & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorVals = nxtVals + nverts;
	    pickColorValsUC = (uchar*)nxtVals;
	    nxtn = 3;
	  }
	  if (damode & CGO_ACCESSIBILITY_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    accessibilityVals = nxtVals;
	    nxtn = 1;
	  }
	  pc = origpc + 1;
	  notHaveValue = damode;
	  end = 0;
	  while(!err && !end && (op = (CGO_MASK & CGO_read_int(pc)))) {
	    switch (op) {
	    case CGO_NORMAL:
	      normalVals[pl] = pc[0]; normalVals[pl+1] = pc[1]; normalVals[pl+2] = pc[2];
	      notHaveValue = notHaveValue ^ CGO_NORMAL_ARRAY;
	      break;
	    case CGO_COLOR:
	      colorVals[plc] = pc[0]; colorVals[plc+1] = pc[1]; colorVals[plc+2] = pc[2]; colorVals[plc+3] = I->alpha;
	      notHaveValue = notHaveValue ^ CGO_COLOR_ARRAY;
	      break;
	    case CGO_PICK_COLOR:
	      /* TODO need to move uchar and index/bond separately into pickColorVals */
	      CGO_put_int(&pickColorVals[pla * 2], CGO_get_int(pc));
	      CGO_put_int(&pickColorVals[pla * 2 + 1], CGO_get_int(pc+1));
	      cgo->current_pick_color_index = (int)pc[0];
	      cgo->current_pick_color_bond = (int)pc[1];
	      notHaveValue = notHaveValue ^ CGO_PICK_COLOR_ARRAY;
	      break;
	    case CGO_ACCESSIBILITY:
	      cgo->current_accessibility = pc[0];
	      break;
	    case CGO_VERTEX:
	      if (notHaveValue & CGO_NORMAL_ARRAY){
		tmp_ptr = &normalVals[pl-3];
		normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];		
	      }
	      if (notHaveValue & CGO_COLOR_ARRAY){
		tmp_ptr = &colorVals[plc-4];
		colorVals[plc] = tmp_ptr[0]; colorVals[plc+1] = tmp_ptr[1]; 
		colorVals[plc+2] = tmp_ptr[2];	colorVals[plc+3] = tmp_ptr[3];
	      }
	      if (notHaveValue & CGO_PICK_COLOR_ARRAY){
		CGO_put_int(pickColorVals + plp, cgo->current_pick_color_index);
		plp++;
		CGO_put_int(pickColorVals + plp, cgo->current_pick_color_bond);
		plp++;
	      }
	      if (accessibilityVals){
		accessibilityVals[pla] = cgo->current_accessibility;
	      }
	      vertexVals[pl++] = pc[0]; vertexVals[pl++] = pc[1]; vertexVals[pl++] = pc[2];
	      plc+=4;
	      pla++;
	      notHaveValue = damode;
	      break;
	    case CGO_END:
	      end = 1;
	    default:
	      break;
	    }
	    sz = CGO_sz[op];
	    pc += sz;
	  }
	}
	save_pc = pc;
	op = CGO_NULL;
      }
      break;
    case CGO_ALPHA:
      I->alpha = *pc;
    default:
      sz = CGO_sz[op];
      nc = CGO_add(cgo, sz + 1);
      *(nc++) = *(pc - 1);
      while(sz--)
        *(nc++) = *(pc++);
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  return (cgo);
#else
  return NULL;
#endif
}

#ifdef _PYMOL_CGO_DRAWBUFFERS
void CGOFreeVBOs(CGO *I){
  register float *pc = I->op;
  register int op;
  float *save_pc;
  int numbufs = 0, bufoffset;

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    numbufs = 0;
    bufoffset = 0;
    switch (op) {
    case CGO_DRAW_SPHERE_BUFFERS:
      numbufs = 3;
      bufoffset = 2;
    case CGO_DRAW_CYLINDER_BUFFERS:
      if (!numbufs){
	numbufs = 5;
	bufoffset = 2;
      }
    case CGO_DRAW_BUFFERS:
      if (!numbufs){
	numbufs = 4;
	bufoffset = 4;
      }
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      if (!numbufs){
	numbufs = 4;
	bufoffset = 4;
      }
    case CGO_DRAW_BUFFERS_INDEXED:
      if (!numbufs){
	numbufs = 5;
	bufoffset = 5;
      }
      {
	int i, buf;
	for (i=0; i<numbufs; i++){
	  buf = CGO_get_int(pc+bufoffset+i);
	  if (buf){
	    CShaderMgr_AddVBOToFree(I->G->ShaderMgr, buf);
	  }
	}
	  switch (op){
	  case CGO_DRAW_BUFFERS_INDEXED:
	    {
	      int nverts = CGO_get_int(pc + 4);
	      pc += nverts*3 + 10;
	      save_pc += nverts*3 + 10 ;
	    }
	    break;
	  case CGO_DRAW_BUFFERS_NOT_INDEXED:
	    {
	      int nverts = CGO_get_int(pc + 3);
	      pc += nverts*3 + 8;
	      save_pc += nverts*3 + 8 ;
	    }
	    break;
	  }
      }
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {      
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	int nvals = narrays*nverts;
	pc += nvals + 4;
	save_pc += nvals + 4 ;
      }
      break;
#endif
    default:
      break;
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
}
CGO *CGOOptimizeToVBO(CGO * I, int est)
{
  CGO *cgo;

  register float *pc = I->op;
  register float *nc;
  register int op;
  float *save_pc;
  int sz;
  int ndrawbuffers = 0;
  short err = 0;
  short has_draw_buffer = false;
  cgo = CGONewSized(I->G, I->c + est);
#ifndef _PYMOL_CGO_DRAWARRAYS
  ndrawbuffers;
#endif
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    err = 0;
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {      
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	if (true){
	  float *vertexVals = 0, *nxtVals = 0, *colorVals = 0, *normalVals;
	  uchar *pickColorVals;
	  short nxtn = 3, bufpl = 0;
	  uint bufs[4], allbufs[4] = { 0, 0, 0, 0 };
	  
	  glGenBuffers(countBitsInt(arrays), bufs);

	  nxtVals = vertexVals = pc + 4;
	  glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
	  allbufs[0] = bufs[bufpl++];
	  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nverts*3, vertexVals, GL_STATIC_DRAW);
	  if (arrays & CGO_NORMAL_ARRAY){
	    nxtVals = normalVals = vertexVals + (nxtn*nverts);
	    glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
	    allbufs[1] = bufs[bufpl++];
	    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nverts*3, normalVals, GL_STATIC_DRAW);
	  }
	  if (arrays & CGO_COLOR_ARRAY){
	    nxtVals = colorVals = nxtVals + (nxtn*nverts);
	    nxtn = 4;
	    glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
	    allbufs[2] = bufs[bufpl++];
	    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*nverts*4, colorVals, GL_STATIC_DRAW);
	  }
	  if (arrays & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorVals = (uchar*)nxtVals;
	    nxtn = 3;
	    glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
	    allbufs[3] = bufs[bufpl++];
	    glBufferData(GL_ARRAY_BUFFER, sizeof(uchar)*nverts*4, pickColorVals, GL_STATIC_DRAW);
	  }
	  //	  printf("CGODrawBuffers: mode=%d allbufs: %d %d %d %d\n", mode, allbufs[0], allbufs[1], allbufs[2], allbufs[3]);
	  ndrawbuffers++;
	  CGODrawBuffers(cgo, mode, arrays, nverts, allbufs);
	  has_draw_buffer = true;
	  {
	    int nvals = narrays*nverts ;
	    pc += nvals + 4;
	    save_pc += nvals + 4 ;
	  }
	} else {
	  GLfloat *vals = CGODrawArrays(cgo, mode, arrays, nverts);
	  int nvals = narrays*nverts, onvals;
	  onvals = nvals;
	  pc += 4;
	  while(nvals--)
	    *(vals++) = *(pc++);
	  save_pc += onvals + 4 ;
	}
      }
	break;
#endif
    case CGO_END:
      if (!err){
	PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOOptimizeToVBO: CGO_END encountered, should call CGOCombineBeginEnd before CGOOptimizeToVBO\n" ENDFB(I->G);
	err = true;
      }
    case CGO_VERTEX:
      if (!err){
	PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOOptimizeToVBO: CGO_VERTEX encountered, should call CGOCombineBeginEnd before CGOOptimizeToVBO\n" ENDFB(I->G);      
	err = true;
      }
    case CGO_BEGIN:
      if (!err){
	PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOOptimizeToVBO: CGO_BEGIN encountered, should call CGOCombineBeginEnd before CGOOptimizeToVBO\n" ENDFB(I->G);      
	err = true;
      }
    case CGO_ALPHA:
      I->alpha = *pc;
    default:
      sz = CGO_sz[op];
      nc = CGO_add(cgo, sz + 1);
      *(nc++) = *(pc - 1);
      while(sz--)
        *(nc++) = *(pc++);
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
  CGOStop(cgo);
  if (has_draw_buffer){
    cgo->has_draw_buffers = true;
  }
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  //  printf("CGOOptimizeToVBO: End: ndrawbuffers=%d\n", ndrawbuffers);
  return (cgo);
}

#define set_min_max(mn, mx, pt) {		\
  if (mn[0]>*pt) mn[0] = *pt; \
  if (mn[1]>*(pt+1)) mn[1] = *(pt+1);		\
  if (mn[2]>*(pt+2)) mn[2] = *(pt+2);		\
  if (mx[0]<*pt) mx[0] = *pt; \
  if (mx[1]<*(pt+1)) mx[1] = *(pt+1);		\
  if (mx[2]<*(pt+2)) mx[2] = *(pt+2);}


void CGOCountNumVertices(CGO *I, int *num_total_vertices, int *num_total_indexes,
			 int *num_total_vertices_lines, int *num_total_indexes_lines,
			 int *num_total_vertices_points){
  register float *pc = I->op;
  register int op;
  float *save_pc;
  int verts_skipped = 0;
  short err = 0;

#ifndef _PYMOL_CGO_DRAWARRAYS
  verts_skipped;
#endif

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    err = 0;
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {      
	int mode = CGO_get_int(pc), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	short shouldCompress = false, shouldCompressLines = false, shouldCompressPoints = false;
	switch(mode){
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
	  verts_skipped += nverts;
	  {
	    int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	    int nvals = narrays*nverts, onvals ;
	    onvals = nvals;
	    pc += 4;
	    save_pc += onvals + 4 ;
	  }
	} else if (shouldCompressLines) {
	  int nvals = narrays*nverts ;
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	  *num_total_vertices_lines += nverts;
	  switch(mode){
	  case GL_LINE_LOOP:
	    *num_total_indexes_lines += 2 * nverts;
	    break;
	  case GL_LINE_STRIP:
	    *num_total_indexes_lines += 2 * (nverts - 1);
	    break;
	  case GL_LINES:
	    *num_total_indexes_lines += nverts;
	    break;
	  }
	} else if (shouldCompress){
	  int nvals = narrays*nverts ;
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	  *num_total_vertices += nverts;
	  switch(mode){
	  case GL_TRIANGLE_FAN:
	    *num_total_indexes += 3 * (nverts - 2);
	    break;
	  case GL_TRIANGLE_STRIP:
	    *num_total_indexes += 3 * (nverts - 2);
	    break;
	  case GL_TRIANGLES:
	    *num_total_indexes += nverts;
	    break;
	  }
	} else if (shouldCompressPoints){
	  int nvals = narrays*nverts ;
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	  *num_total_vertices_points += nverts;
	}
      }
	break;
#endif
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
    case CGO_ALPHA:
      I->alpha = *pc;
    default:
      break;
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
}

void SetVertexValuesForVBO(PyMOLGlobals * G, CGO *cgo, int arrays, int pl, int plc, int cnt, int incr, float *vertexValsDA, float *normalValsDA, float *colorValsDA, float *pickColorValsDA, 
			   float *vertexVals, uchar *normalValsC, float *normalVals, uchar *colorValsUC, float *colorVals, float *pickColorVals, float *accessibilityVals, float *accessibilityValsDA){
  int pl2 = pl + 1, pl3 = pl + 2;
  int plc2 = plc + 1, plc3 = plc + 2, plc4 = plc + 3;
  int c, c2, c3;
  int cc, cc2, cc3, cc4;
  int pcc = incr * 2, pcco = cnt * 2;
  c = cnt * 3; c2 = c + 1; c3 = c + 2;
  cc = cnt * 4; cc2 = cc + 1; cc3 = cc + 2; cc4 = cc + 3;
  vertexVals[pl] = vertexValsDA[c]; vertexVals[pl2] = vertexValsDA[c2]; vertexVals[pl3] = vertexValsDA[c3];
  if (SettingGet(G, cSetting_cgo_shader_ub_normal)){
    if (normalValsC){
      if (arrays & CGO_NORMAL_ARRAY){
	normalValsC[pl] = CLIP_NORMAL_VALUE(normalValsDA[c]); normalValsC[pl2] = CLIP_NORMAL_VALUE(normalValsDA[c2]); normalValsC[pl3] = CLIP_NORMAL_VALUE(normalValsDA[c3]);
      } else {
	normalValsC[pl] = CLIP_NORMAL_VALUE(cgo->normal[0]); normalValsC[pl2] = CLIP_NORMAL_VALUE(cgo->normal[1]); normalValsC[pl3] = CLIP_NORMAL_VALUE(cgo->normal[2]);
      }
    }
  } else {
    if (normalVals){
      if (arrays & CGO_NORMAL_ARRAY){
	normalVals[pl] = normalValsDA[c]; normalVals[pl2] = normalValsDA[c2]; normalVals[pl3] = normalValsDA[c3];
      } else {
	normalVals[pl] = cgo->normal[0]; normalVals[pl2] = cgo->normal[1]; normalVals[pl3] = cgo->normal[2];
      }
    }
  }
  if (SettingGet(G, cSetting_cgo_shader_ub_color)){
    if (arrays & CGO_COLOR_ARRAY){
      colorValsUC[plc] = CLIP_COLOR_VALUE(colorValsDA[cc]); colorValsUC[plc2] = CLIP_COLOR_VALUE(colorValsDA[cc2]); 
      colorValsUC[plc3] = CLIP_COLOR_VALUE(colorValsDA[cc3]); colorValsUC[plc4] = CLIP_COLOR_VALUE(colorValsDA[cc4]);
    } else {
      colorValsUC[plc] = CLIP_COLOR_VALUE(cgo->color[0]); colorValsUC[plc2] = CLIP_COLOR_VALUE(cgo->color[1]);
      colorValsUC[plc3] = CLIP_COLOR_VALUE(cgo->color[2]); colorValsUC[plc4] = CLIP_COLOR_VALUE(cgo->alpha);
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY){
      colorVals[plc] = colorValsDA[cc]; colorVals[plc2] = colorValsDA[cc2];
      colorVals[plc3] = colorValsDA[cc3]; colorVals[plc4] = colorValsDA[cc4];
    } else {
      colorVals[plc] = cgo->color[0]; colorVals[plc2] = cgo->color[1];
      colorVals[plc3] = cgo->color[2]; colorVals[plc4] = cgo->alpha;
    }
  }
  if (arrays & CGO_PICK_COLOR_ARRAY){
    CGO_put_int(pickColorVals + pcc, CGO_get_int(pickColorValsDA + pcco));
    CGO_put_int(pickColorVals + pcc + 1, CGO_get_int(pickColorValsDA + pcco + 1));
  } else {
    CGO_put_int(pickColorVals + pcc, cgo->current_pick_color_index);
    CGO_put_int(pickColorVals + pcc + 1, cgo->current_pick_color_bond);
  }
  if (arrays & CGO_ACCESSIBILITY_ARRAY){
    accessibilityVals[pl/3] = accessibilityValsDA[cnt];
  }
}

void OptimizePointsToVBO(CGO *I, CGO *cgo, int num_total_vertices_points, float *min, float *max, short *has_draw_buffer){
  float *vertexVals = 0, *colorVals = 0, *normalVals = 0;
  float *pickColorVals;
  int pl = 0, plc = 0, idxpl = 0, vpl = 0, tot, nxtn;
  uchar *colorValsUC = 0;
  uchar *normalValsC = 0;
  int numbufs = 0, bufoffset = 0;
  short has_normals = 0;
  float *pc = I->op;
  int op;
  float *save_pc;
  short err = 0;
  
#ifndef _PYMOL_CGO_DRAWARRAYS
  pl; plc; idxpl; vpl;
#endif
  
  cgo->alpha = 1.f;
  cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;
  
  tot = num_total_vertices_points * (3 * 5) ;
  //    tot = num_total_indexes * (3 * 3 + 2) ;
  /* NOTE/TODO: Not sure why 3*5 needs to be used, but 3*3+2, which is the 
     correct length, crashes in glBufferData */
  vertexVals = Alloc(float, tot);
  if (!vertexVals){
    PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() vertexVals could not be allocated\n" ENDFB(I->G);	
  }
  normalVals = vertexVals + 3 * num_total_vertices_points;
  nxtn = 3;
  if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
    normalValsC = (uchar*) normalVals;
    nxtn = 1;
  }
  colorVals = normalVals + nxtn * num_total_vertices_points;
  if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
    colorValsUC = (uchar*) colorVals;
    nxtn = 1;
  } else {
    nxtn = 4;
  }
  pickColorVals = (colorVals + nxtn * num_total_vertices_points);
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    err = 0;
    numbufs = 0;
    switch (op) {
    case CGO_BOUNDING_BOX:
      {
	register float *nc, *newpc = pc;
	int sz;
	sz = CGO_sz[op];
	nc = CGO_add(cgo, sz + 1);
	*(nc++) = *(pc - 1);
	while(sz--)
	  *(nc++) = *(newpc++);
      }
      break;
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_SPHERE_BUFFERS:
      numbufs = 3;
      bufoffset = 2;
    case CGO_DRAW_CYLINDER_BUFFERS:
      if (!numbufs){
	numbufs = 5;
	bufoffset = 2;
      }
    case CGO_DRAW_BUFFERS:
      if (!numbufs){
	numbufs = 4;
	bufoffset = 4;
      }
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      if (!numbufs){
	numbufs = 4;
	bufoffset = 4;
      }
    case CGO_DRAW_BUFFERS_INDEXED:
      if (!numbufs){
	numbufs = 5;
	bufoffset = 5;
      }
      {
	int i, sz;
	register float *nc, *newpc = pc;
	sz = CGO_sz[op];
	nc = CGO_add(cgo, sz + 1);
	*(nc++) = *(pc - 1);
	while(sz--)
	  *(nc++) = *(newpc++);
	
	for (i=0; i<numbufs; i++){
	  *(pc+bufoffset+i) = 0;
	}
      }
      break;
#endif
    case CGO_NORMAL:
      cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
      has_normals = 1;
      break;
    case CGO_COLOR:
      cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
      break;
    case CGO_ALPHA:
      cgo->alpha = *pc;
      break;
    case CGO_PICK_COLOR:
      cgo->current_pick_color_index = (int)pc[0];
      cgo->current_pick_color_bond = (int)pc[1];
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	short shouldCompress = false;
	switch(mode){
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
	  int nvals = narrays*nverts, cnt, nxtn = 3 ,incr=0;
	  float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0, *normalValsDA = 0;
	  float *pickColorValsDA = NULL, *pickColorValsTMP;
	  
	  nxtVals = vertexValsDA = pc + 4;
	  
	  for (cnt=0; cnt<nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  if (arrays & CGO_NORMAL_ARRAY){
	    nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	  }
	  if (arrays & CGO_COLOR_ARRAY){
	    nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	    nxtn = 4;
	  }
	  if (arrays & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorValsDA = nxtVals + nverts;
	    nxtn = 3;
	  }
	  pickColorValsTMP = pickColorVals + (idxpl * 2);
	  switch (mode){
	  case GL_POINTS:
	    for (cnt = 0; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, NULL, NULL);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  }
	  vpl += nverts;
	  
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	} else {
	  {
	    int nvals = narrays*nverts ;
	    pc += nvals + 4;
	    save_pc += nvals + 4 ;
	  }
	}
      }
      break;
#endif
    default:
      break;
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
  {
    uint bufs[4] = {0, 0, 0, 0 }, allbufs[4] = { 0, 0, 0, 0 };
    short bufpl = 0;
    GLenum err ;
    //      printf("CGOOptimizeToVBOIndexed: End: num_total_vertices=%d num_total_indexes=%d verts_skipped=%d\n", num_total_vertices, num_total_indexes, verts_skipped);
    if ((err = glGetError())){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() BEFORE glGenBuffers returns err=%d\n", err ENDFB(I->G);	
    }
    
    glGenBuffers(3, bufs);
    if ((err = glGetError())){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glGenBuffers returns err=%d\n", err ENDFB(I->G);	
    }
    glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
    if ((err = glGetError())){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
    }
    if (!glIsBuffer(bufs[bufpl])){
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
    } else {
      allbufs[0] = bufs[bufpl++];
      glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices_points*3, vertexVals, GL_STATIC_DRAW);      
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
      }
    }
    
    /*      NO NORMALS IN POINTS */
    if (!has_normals){
      if (bufs[bufpl]){
	if (glIsBuffer(bufs[bufpl])){
	  CShaderMgr_AddVBOToFree(I->G->ShaderMgr, bufs[bufpl]);
	}
	bufs[bufpl] = 0;
      }
      bufpl++;
    } else {
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 3;
	allbufs[1] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices_points*sz, normalVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
    }
    glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
    if ((err = glGetError())){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
    }
    if (!glIsBuffer(bufs[bufpl])){
      PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
    } else {
      short sz = 4;
      allbufs[2] = bufs[bufpl++];
      if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	sz = 1;
      }
      glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices_points*sz, colorVals, GL_STATIC_DRAW);      
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
      }
    }
    {
      GLfloat *newPickColorVals ;
      newPickColorVals = CGODrawBuffersNotIndexed(cgo, GL_POINTS, CGO_VERTEX_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, num_total_vertices_points, allbufs);
      memcpy(newPickColorVals + num_total_vertices_points, pickColorVals, num_total_vertices_points * 2 * sizeof(float));
      *has_draw_buffer = true;
    }
  }
  FreeP(vertexVals);
  /* END GL_POINTS */
  //    printf("num_total_vertices_points=%d\n", num_total_vertices_points);
}

CGO *CGOOptimizeToVBONotIndexed(CGO * I, int est)
{
  CGO *cgo;
  register float *pc = I->op;
  register int op;
  float *save_pc;
  int num_total_vertices = 0, num_total_indexes = 0, num_total_vertices_lines = 0, num_total_indexes_lines = 0,
    num_total_vertices_points = 0;
  short err = 0;
  short has_draw_buffer = false;
  float min[3] = { MAXFLOAT, MAXFLOAT, MAXFLOAT }, max[3] = { -MAXFLOAT, -MAXFLOAT, -MAXFLOAT };
  int ambient_occlusion = 0;
  cgo = CGONewSized(I->G, I->c + est);

  CGOCountNumVertices(I, &num_total_vertices, &num_total_indexes,
                         &num_total_vertices_lines, &num_total_indexes_lines,
                         &num_total_vertices_points);

  if (num_total_vertices_points>0){
    OptimizePointsToVBO(I, cgo, num_total_vertices_points, min, max, &has_draw_buffer);
  }
  if (num_total_indexes>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals;
    float *pickColorVals, *accessibilityVals = 0;
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, tot, nxtn;
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;
    int numbufs = 0, bufoffset = 0;
    pc = I->op;

#ifndef _PYMOL_CGO_DRAWARRAYS
    pl; plc; idxpl; vpl;
#endif

    cgo->alpha = 1.f;
    cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;

    tot = num_total_indexes * (3 * 6) ;
    //    tot = num_total_indexes * (3 * 3 + 2) ;
    /* NOTE/TODO: Not sure why 3*5 needs to be used, but 3*3+2, which is the 
       correct length, crashes in glBufferData */
    vertexVals = Alloc(float, tot);
    if (!vertexVals){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() vertexVals could not be allocated\n" ENDFB(I->G);	
    }
    normalVals = vertexVals + 3 * num_total_indexes;
    nxtn = 3;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      nxtn = 1;
    }

    colorVals = normalVals + nxtn * num_total_indexes;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      nxtn = 1;
    } else {
      nxtn = 4;
    }
    pickColorVals = (colorVals + nxtn * num_total_indexes);
    nxtn = 3;
    accessibilityVals = pickColorVals + nxtn * num_total_indexes;
    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      save_pc = pc;
      err = 0;
      numbufs = 0;
      switch (op) {
      case CGO_BOUNDING_BOX:
	{
	  register float *nc, *newpc = pc;
	  int sz;
	  sz = CGO_sz[op];
	  nc = CGO_add(cgo, sz + 1);
	  *(nc++) = *(pc - 1);
	  while(sz--)
	    *(nc++) = *(newpc++);
	}
	break;
#ifdef _PYMOL_CGO_DRAWBUFFERS
      case CGO_DRAW_SPHERE_BUFFERS:
	numbufs = 3;
	bufoffset = 2;
      case CGO_DRAW_CYLINDER_BUFFERS:
	if (!numbufs){
	  numbufs = 5;
	  bufoffset = 2;
	}
      case CGO_DRAW_BUFFERS:
	if (!numbufs){
	  numbufs = 4;
	  bufoffset = 4;
	}
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
	if (!numbufs){
	  numbufs = 4;
	  bufoffset = 4;
	}
      case CGO_DRAW_BUFFERS_INDEXED:
	if (!numbufs){
	  numbufs = 5;
	  bufoffset = 5;
	}

	{
	  int i, sz;
	  register float *nc, *newpc = pc;
	  sz = CGO_sz[op];
	  nc = CGO_add(cgo, sz + 1);
	  *(nc++) = *(pc - 1);
	  while(sz--)
	    *(nc++) = *(newpc++);
	  
	  for (i=0; i<numbufs; i++){
	    *(pc+bufoffset+i) = 0;
	  }
	}
	break;
#endif
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
	cgo->current_accessibility = pc[0];
	break;
      case CGO_PICK_COLOR:
	cgo->current_pick_color_index = (int)pc[0];
	cgo->current_pick_color_bond = (int)pc[1];
	break;
#ifdef _PYMOL_CGO_DRAWARRAYS
      case CGO_DRAW_ARRAYS:
	{
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	short shouldCompress = false;
	switch(mode){
	case GL_TRIANGLE_FAN:
	case GL_TRIANGLE_STRIP:
	case GL_TRIANGLES:
	  shouldCompress = true;
	default:
	  break;
	}
	if (shouldCompress){	
	  int nvals = narrays*nverts, cnt, nxtn = 3 ,incr=0;
	  float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0, *normalValsDA = NULL, *accessibilityValsDA = NULL;
	  float *pickColorValsDA = NULL, *pickColorValsTMP;

	  nxtVals = vertexValsDA = pc + 4;

	  for (cnt=0; cnt<nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  if (arrays & CGO_NORMAL_ARRAY){
	    nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	  }

	  if (arrays & CGO_COLOR_ARRAY){
	    nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	    nxtn = 4;
	  }
	  if (arrays & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorValsDA = nxtVals + nverts;
	    nxtn = 3;
	  }
	  pickColorValsTMP = pickColorVals + (idxpl * 2);
	  if (arrays & CGO_ACCESSIBILITY_ARRAY){
	    if (!ambient_occlusion && incr){
	      for (cnt=0; cnt<incr;cnt++){
		/* if ambient_occlusion, need to fill in the array */
		accessibilityVals[cnt] = 1.f;
	      }
	    }
	    ambient_occlusion = 1;
	    accessibilityValsDA = nxtVals + nxtn*nverts;
	    /*	    printf("nverts=%d\n", nverts);
	    for (cnt=0; cnt<nverts; cnt+=1000){
	      printf("\taccessibilityValsDA[%d]=%f\n", cnt, accessibilityValsDA[cnt]);
	      }*/
	  } else {
	    if (ambient_occlusion){
	      for (cnt=incr; cnt<incr+nverts;cnt++){
		/* if ambient_occlusion, need to fill in the array */
		accessibilityVals[cnt] = 1.f;
	      }
	    }
	  }
	  /*	  if (idxpl + nverts > num_total_indexes){
	    printf("num_total_indexes=%d mode=%d nverts=%d idxpl=%d\n", num_total_indexes, mode, nverts, idxpl);
	    }*/
	  switch (mode){
	  case GL_TRIANGLES:
	    for (cnt = 0; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  case GL_TRIANGLE_STRIP:
	    for (cnt = 2; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt-2, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt-1, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  case GL_TRIANGLE_FAN:
	    for (cnt = 2; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, 0, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt - 1, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, accessibilityVals, accessibilityValsDA);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  }

	  //	  pl += 3 * nverts;
	  //	  plc += 4 * nverts;
	  vpl += nverts;

	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	} else {
	  {
	    int nvals = narrays*nverts ;
	    pc += nvals + 4;
	    save_pc += nvals + 4 ;
	  }
	}
	}
	break;
#endif
      default:
	break;
      }
      pc = save_pc;
      pc += CGO_sz[op];
    }
    
    {
      uint bufs[4] = {0, 0, 0, 0 }, allbufs[4] = { 0, 0, 0, 0 };
      short bufpl = 0, numbufs = 3;
      GLenum err ;
      //      printf("CGOOptimizeToVBOIndexed: End: num_total_vertices=%d num_total_indexes=%d verts_skipped=%d\n", num_total_vertices, num_total_indexes, verts_skipped);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() BEFORE glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }
      if (ambient_occlusion){
	numbufs++;
      }
      glGenBuffers(numbufs, bufs);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	allbufs[0] = bufs[bufpl++];
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes*3, vertexVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 3;
	allbufs[1] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes*sz, normalVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 4;
	allbufs[2] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes*sz, colorVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      if (ambient_occlusion){
	glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
	}
	if (!glIsBuffer(bufs[bufpl])){
	  PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
	} else {
	  allbufs[3] = bufs[bufpl++];
	  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes, accessibilityVals, GL_STATIC_DRAW);      
	  if ((err = glGetError())){
	    PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	  }
	}
      }
      {
	GLfloat *newPickColorVals ;
	int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY;
	if (ambient_occlusion){
	  arrays |= CGO_ACCESSIBILITY_ARRAY;
	}
	newPickColorVals = CGODrawBuffersNotIndexed(cgo, GL_TRIANGLES, arrays, num_total_indexes, allbufs);
	memcpy(newPickColorVals + num_total_indexes, pickColorVals, num_total_indexes * 2 * sizeof(float));
	has_draw_buffer = true;
      }
    }
    FreeP(vertexVals);
  }
  if (num_total_indexes_lines>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals;
    float *pickColorVals;
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, tot, nxtn;
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;
#ifndef _PYMOL_CGO_DRAWARRAYS
    pl; plc; idxpl; vpl;
#endif
    pc = I->op;
    cgo->alpha = 1.f;
    cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;

    tot = num_total_indexes_lines * (3 * 5) ;
    //    tot = num_total_indexes * (3 * 3 + 2) ;
    /* NOTE/TODO: Not sure why 3*5 needs to be used, but 3*3+2, which is the 
       correct length, crashes in glBufferData */
    vertexVals = Alloc(float, tot);
    if (!vertexVals){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() vertexVals could not be allocated\n" ENDFB(I->G);	
    }
    normalVals = vertexVals + 3 * num_total_indexes_lines;
    nxtn = 3;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      nxtn = 1;
    }

    colorVals = normalVals + nxtn * num_total_indexes_lines;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      nxtn = 1;
    } else {
      nxtn = 4;
    }
    pickColorVals = (colorVals + nxtn * num_total_indexes_lines);
    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      save_pc = pc;
      err = 0;
      switch (op) {
      case CGO_NORMAL:
	cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
	break;
      case CGO_COLOR:
	cgo->color[0] = *pc; cgo->color[1] = *(pc + 1); cgo->color[2] = *(pc + 2);
	break;
      case CGO_ACCESSIBILITY:
	cgo->current_accessibility = pc[0];
	break;
      case CGO_ALPHA:
	cgo->alpha = *pc;
	break;
      case CGO_PICK_COLOR:
	cgo->current_pick_color_index = (int)pc[0];
	cgo->current_pick_color_bond = (int)pc[1];
	break;
#ifdef _PYMOL_CGO_DRAWARRAYS
      case CGO_DRAW_ARRAYS:
	{
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	short shouldCompress = false;
	switch(mode){
	case GL_LINE_LOOP:
	case GL_LINE_STRIP:
	case GL_LINES:
	  shouldCompress = true;
	default:
	  break;
	}
	
	if (shouldCompress){	
	  int nvals = narrays*nverts, cnt, nxtn = 3, incr = 0;
	  float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0, *normalValsDA = NULL;
	  float *pickColorValsDA = NULL, *pickColorValsTMP;

	  nxtVals = vertexValsDA = pc + 4;

	  for (cnt=0; cnt<nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  if (arrays & CGO_NORMAL_ARRAY){
	    nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	  }

	  if (arrays & CGO_COLOR_ARRAY){
	    nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	    nxtn = 4;
	  }
	  if (arrays & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorValsDA = nxtVals + nverts;
	    nxtn = 3;
	  }
	  pickColorValsTMP = pickColorVals + (idxpl * 2);
	  /*	  if (idxpl + nverts > num_total_indexes){
	    printf("num_total_indexes=%d mode=%d nverts=%d idxpl=%d\n", num_total_indexes, mode, nverts, idxpl);
	    }*/
	  switch (mode){
	  case GL_LINES:
	    for (cnt = 0; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, NULL, NULL);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  case GL_LINE_STRIP:
	    for (cnt = 1; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt-1, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, NULL, NULL);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, NULL, NULL);
	      idxpl++; pl += 3; plc += 4;
	    }
	    break;
	  case GL_LINE_LOOP:
	    for (cnt = 1; cnt < nverts; cnt++){
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt-1, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, NULL, NULL);
	      idxpl++; pl += 3; plc += 4;
	      SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, cnt, incr++, 
				    vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				    vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				    pickColorValsTMP, NULL, NULL);
	      idxpl++; pl += 3; plc += 4;
	    }
	    SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, 0, incr++, 
				  vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				  vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				  pickColorValsTMP, NULL, NULL);
	    idxpl++; pl += 3; plc += 4;
	    SetVertexValuesForVBO(I->G, cgo, arrays, pl, plc, nverts-1, incr++, 
				  vertexValsDA, normalValsDA, colorValsDA, pickColorValsDA, 
				  vertexVals, normalValsC, normalVals, colorValsUC, colorVals, 
				  pickColorValsTMP, NULL, NULL);
	    idxpl++; pl += 3; plc += 4;
	    break;
	  }

	  //	  pl += 3 * nverts;
	  //	  plc += 4 * nverts;
	  vpl += nverts;
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	} else {
	  {
	    int nvals = narrays*nverts ;
	    pc += nvals + 4;
	    save_pc += nvals + 4 ;
	  }
	}
	}
	break;
#endif
      default:
	break;
      }
      pc = save_pc;
      pc += CGO_sz[op];
    }
    {
      uint bufs[4] = {0, 0, 0, 0 }, allbufs[4] = { 0, 0, 0, 0 };
      short bufpl = 0;
      GLenum err ;
      //      printf("CGOOptimizeToVBOIndexed: End: num_total_vertices=%d num_total_indexes=%d verts_skipped=%d\n", num_total_vertices, num_total_indexes, verts_skipped);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() BEFORE glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }

      glGenBuffers(3, bufs);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	allbufs[0] = bufs[bufpl++];
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes_lines*3, vertexVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 3;
	allbufs[1] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes_lines*sz, normalVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBONotIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 4;
	allbufs[2] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes_lines*sz, colorVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBONotIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }

      {
	GLfloat *newPickColorVals ;
#if defined(OPENGL_ES_2)
	CGODisable(cgo, GL_SHADER_LIGHTING);
#endif
	newPickColorVals = CGODrawBuffersNotIndexed(cgo, GL_LINES, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, num_total_indexes_lines, allbufs);
	memcpy(newPickColorVals + num_total_indexes_lines, pickColorVals, num_total_indexes_lines * 2 * sizeof(float));
	has_draw_buffer = true;
#if defined(OPENGL_ES_2)
	CGOEnable(cgo, GL_SHADER_LIGHTING);
#endif
      }
    }
    FreeP(vertexVals);
  }

  if (num_total_vertices>0 || num_total_vertices_lines>0 || num_total_vertices_points>0){
    CGOBoundingBox(cgo, min, max);
  }

  CGOStop(cgo);
  if (has_draw_buffer){
    cgo->has_draw_buffers = true;
  }
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  return (cgo);
}

CGO *CGOOptimizeToVBOIndexed(CGO * I, int est){
  return (CGOOptimizeToVBOIndexedWithColor(I, est, NULL));
}

CGO *CGOOptimizeToVBOIndexedWithColor(CGO * I, int est, float *color)
{
  CGO *cgo;

  register float *pc = I->op;
  register int op;
  float *save_pc;
  int num_total_vertices = 0, num_total_indexes = 0, num_total_vertices_lines = 0, num_total_indexes_lines = 0,
    num_total_vertices_points = 0;
  short err = 0;
  short has_draw_buffer = false;
  float min[3] = { MAXFLOAT, MAXFLOAT, MAXFLOAT }, max[3] = { -MAXFLOAT, -MAXFLOAT, -MAXFLOAT };

  cgo = CGONewSized(I->G, I->c + est);
  cgo->alpha = 1.f;
  if (color){
    cgo->color[0] = color[0]; cgo->color[1] = color[1]; cgo->color[2] = color[2];
    cgo->alpha = color[3];
  } else {
    cgo->color[0] = 1.f; cgo->color[1] = 1.f; cgo->color[2] = 1.f;
  }

  CGOCountNumVertices(I, &num_total_vertices, &num_total_indexes,
                         &num_total_vertices_lines, &num_total_indexes_lines,
                         &num_total_vertices_points);

  if (num_total_vertices_points>0){
    /* This does not need to be indexed (for now) */
    OptimizePointsToVBO(I, cgo, num_total_vertices_points, min, max, &has_draw_buffer);
  }

  if (num_total_vertices>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals, *accessibilityVals = 0;
    float *pickColorVals;
    uint *vertexIndexes; 
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, tot, nxtn;
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;
    short ambient_occlusion = 0;
#ifndef _PYMOL_CGO_DRAWARRAYS
    pl; plc; idxpl; vpl;
#endif
    pc = I->op;
    vertexIndexes = Alloc(uint, num_total_indexes);
    if (!vertexIndexes){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() vertexIndexes could not be allocated\n" ENDFB(I->G);	
    }
    tot = num_total_vertices * (3 * 6) ;
    //    tot = num_total_vertices * (3 * 3 + 2) ;
    /* NOTE/TODO: Not sure why 3*5 needs to be used, but 3*3+2, which is the 
       correct length, crashes in glBufferData */
    vertexVals = Alloc(float, tot);
    if (!vertexVals){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() vertexVals could not be allocated\n" ENDFB(I->G);	
    }
    normalVals = vertexVals + 3 * num_total_vertices;
    nxtn = 3;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      nxtn = 1;
    }

    colorVals = normalVals + nxtn * num_total_vertices;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      nxtn = 1;
    } else {
      nxtn = 4;
    }
    pickColorVals = (colorVals + nxtn * num_total_vertices);
    accessibilityVals = pickColorVals + 3 * num_total_vertices;

    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      save_pc = pc;
      err = 0;
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
	cgo->current_pick_color_index = (int)pc[0];
	cgo->current_pick_color_bond = (int)pc[1];
	break;
#ifdef _PYMOL_CGO_DRAWARRAYS
      case CGO_DRAW_ARRAYS:
	{
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	short shouldCompress = false;
	switch(mode){
	case GL_TRIANGLE_FAN:
	case GL_TRIANGLE_STRIP:
	case GL_TRIANGLES:
	  shouldCompress = true;
	default:
	  break;
	}
	if (shouldCompress){	
	  int nvals = narrays*nverts, cnt, nxtn = 3;
	  float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0, *normalValsDA, *accessibilityValsDA;
	  float *pickColorValsDA, *pickColorValsTMP, *accessibilityValsTMP;

	  nxtVals = vertexValsDA = pc + 4;
	  for (cnt=0; cnt<nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  for (cnt=0; cnt<nverts*3; cnt++){
	    vertexVals[pl + cnt] = vertexValsDA[cnt];
	  }
	  if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	    if (arrays & CGO_NORMAL_ARRAY){
	      nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalValsC[pl + cnt] = CLIP_NORMAL_VALUE(normalValsDA[cnt]);
	      }
	    } else {
	      uchar norm[3] = { CLIP_NORMAL_VALUE(cgo->normal[0]), CLIP_NORMAL_VALUE(cgo->normal[1]), CLIP_NORMAL_VALUE(cgo->normal[2]) };
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalValsC[pl + cnt] = norm[cnt%3];
	      }
	    }
	  } else {
	    if (arrays & CGO_NORMAL_ARRAY){
	      nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalVals[pl + cnt] = normalValsDA[cnt];
	      }
	    } else {
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalVals[pl + cnt] = cgo->normal[cnt%3];
	      }
	    }
	  }
	  if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	    if (arrays & CGO_COLOR_ARRAY){
	      nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
	      }
	      nxtn = 4;
	    } else {
	      uchar col[4] = { CLIP_COLOR_VALUE(cgo->color[0]), CLIP_COLOR_VALUE(cgo->color[1]), CLIP_COLOR_VALUE(cgo->color[2]), CLIP_COLOR_VALUE(cgo->alpha) };
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorValsUC[plc + cnt] = col[cnt%4];
	      }
	    }
	  } else {
	    if (arrays & CGO_COLOR_ARRAY){
	      nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorVals[plc + cnt] = colorValsDA[cnt];
	      }
	      nxtn = 4;
	    } else {
	      float col[4] = { cgo->color[0], cgo->color[1], cgo->color[2], cgo->alpha };
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorVals[plc + cnt] = col[cnt%4];
	      }
	    }
	  }
	  if (arrays & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorValsDA = nxtVals + nverts;
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<nverts; cnt++){
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	    }
	    nxtn = 3;
	  } else {
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<nverts; cnt++){
	      CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_index);
	      CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
	    }
	  }
	  if (arrays & CGO_ACCESSIBILITY_ARRAY){
	    if (!ambient_occlusion){
	      for (cnt=0; cnt<vpl; cnt++){
		accessibilityVals[cnt] = 1.f;
	      }
	    }
	    ambient_occlusion = 1;
	    nxtVals = nxtVals + (nxtn*nverts);
	    accessibilityValsDA = nxtVals;
	    accessibilityValsTMP = accessibilityVals + vpl;
	    for (cnt=0; cnt<nverts; cnt++){
	      accessibilityValsTMP[cnt] = accessibilityValsDA[cnt];
	    }
	  } else {
	    if (ambient_occlusion){
	      accessibilityValsTMP = accessibilityVals + vpl;
	      for (cnt=0; cnt<nverts; cnt++){
		accessibilityValsTMP[cnt] = 1.f;
	      }
	    }
	  }

	  /*	  if (idxpl + nverts > num_total_indexes){
	    printf("num_total_indexes=%d mode=%d nverts=%d idxpl=%d\n", num_total_indexes, mode, nverts, idxpl);
	    }*/
	  switch (mode){
	  case GL_TRIANGLES:
	    for (cnt = 0; cnt < nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_TRIANGLE_STRIP:
	    for (cnt = 2; cnt < nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt - 2;
	      vertexIndexes[idxpl++] = vpl + cnt - 1;
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_TRIANGLE_FAN:
	    for (cnt = 2; cnt < nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl;
	      vertexIndexes[idxpl++] = vpl + cnt - 1;
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  }

	  pl += 3 * nverts;
	  plc += 4 * nverts;
	  vpl += nverts;
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	} else {
	  {
	    int nvals = narrays*nverts ;
	    pc += nvals + 4;
	    save_pc += nvals + 4 ;
	  }
	}
	}
	break;
#endif
      default:
	break;
      }
      pc = save_pc;
      pc += CGO_sz[op];
    }
    {
      uint bufs[5] = {0, 0, 0, 0, 0 }, allbufs[5] = { 0, 0, 0, 0, 0 };
      short bufpl = 0;
      GLenum err ;
      //      printf("CGOOptimizeToVBOIndexed: End: num_total_vertices=%d num_total_indexes=%d verts_skipped=%d\n", num_total_vertices, num_total_indexes, verts_skipped);

      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() BEFORE glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }

      glGenBuffers(5, bufs);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	allbufs[0] = bufs[bufpl++];
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices*3, vertexVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 3;
	allbufs[1] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices*sz, normalVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 4;
	allbufs[2] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices*sz, colorVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	allbufs[3] = bufs[bufpl++];
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint)*num_total_indexes, vertexIndexes, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      if (ambient_occlusion){
	glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
	}
	if (!glIsBuffer(bufs[bufpl])){
	  PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
	} else {
	  allbufs[4] = bufs[bufpl++];
	  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_indexes, accessibilityVals, GL_STATIC_DRAW);      
	  if ((err = glGetError())){
	    PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	  }
	}
      }
      {
	GLfloat *newPickColorVals ;
	int arrays = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY;
	if (ambient_occlusion){
	  arrays |= CGO_ACCESSIBILITY_ARRAY;
	}
	newPickColorVals = CGODrawBuffersIndexed(cgo, GL_TRIANGLES, arrays, num_total_indexes, num_total_vertices, allbufs);
	memcpy(newPickColorVals + num_total_vertices, pickColorVals, num_total_vertices * 2 * sizeof(float));
	has_draw_buffer = true;
      }
    }
    FreeP(vertexIndexes);
    FreeP(vertexVals);
  }
  if (num_total_vertices_lines>0){
    float *vertexVals = 0, *colorVals = 0, *normalVals;
    float *pickColorVals;
    uint *vertexIndexes; 
    uchar *colorValsUC = 0;
    uchar *normalValsC = 0;
    int pl = 0, plc = 0, idxpl = 0, vpl = 0, tot, sz;
#ifndef _PYMOL_CGO_DRAWARRAYS
    pl; plc; idxpl; vpl;
#endif
    pc = I->op;

    vertexIndexes = Alloc(uint, num_total_indexes_lines);
    if (!vertexIndexes){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() vertexIndexes could not be allocated\n" ENDFB(I->G);	
    }
    tot = num_total_vertices_lines * (3 * 5) ;
    //    tot = num_total_vertices * (3 * 3 + 2) ;
    /* NOTE/TODO: Not sure why 3*5 needs to be used, but 3*3+2, which is the 
       correct length, crashes in glBufferData */
    vertexVals = Alloc(float, tot);
    if (!vertexVals){
      PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() vertexVals could not be allocated\n" ENDFB(I->G);	
    }
    normalVals = vertexVals + 3 * num_total_vertices_lines;
    sz = 3;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      normalValsC = (uchar*) normalVals;
      sz = 1;
    }
    colorVals = normalVals + sz * num_total_vertices_lines;
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      colorValsUC = (uchar*) colorVals;
      sz = 1;
    } else {
      sz = 4;
    }
    pickColorVals = (colorVals + sz * num_total_vertices_lines);
    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      save_pc = pc;
      err = 0;
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
	cgo->current_pick_color_index = (int)pc[0];
	cgo->current_pick_color_bond = (int)pc[1];
	break;
#ifdef _PYMOL_CGO_DRAWARRAYS
      case CGO_DRAW_ARRAYS:
	{
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	short shouldCompress = false;
	switch(mode){
	case GL_LINE_LOOP:
	case GL_LINE_STRIP:
	case GL_LINES:
	  shouldCompress = true;
	default:
	  break;
	}
	if (shouldCompress){	
	  int nvals = narrays*nverts, cnt, nxtn = 3;
	  float *vertexValsDA = 0, *nxtVals = 0, *colorValsDA = 0, *normalValsDA;
	  float *pickColorValsDA, *pickColorValsTMP;

	  nxtVals = vertexValsDA = pc + 4;
	  for (cnt=0; cnt<nverts*3; cnt+=3){
	    set_min_max(min, max, &vertexValsDA[cnt]);
	  }
	  for (cnt=0; cnt<nverts*3; cnt++){
	    vertexVals[pl + cnt] = vertexValsDA[cnt];
	  }
	  if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	    if (arrays & CGO_NORMAL_ARRAY){
	      nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalValsC[pl + cnt] = CLIP_NORMAL_VALUE(normalValsDA[cnt]);
	      }
	    } else {
	      uchar norm[3] = { CLIP_NORMAL_VALUE(cgo->normal[0]), CLIP_NORMAL_VALUE(cgo->normal[1]), CLIP_NORMAL_VALUE(cgo->normal[2]) };
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalValsC[pl + cnt] = norm[cnt%3];
	      }
	    }
	  } else {
	    if (arrays & CGO_NORMAL_ARRAY){
	      nxtVals = normalValsDA = vertexValsDA + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalVals[pl + cnt] = normalValsDA[cnt];
	      }
	    } else {
	      for (cnt=0; cnt<nverts*3; cnt++){
		normalVals[pl + cnt] = I->normal[cnt%3];
	      }
	    }
	  }
	  if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	    if (arrays & CGO_COLOR_ARRAY){
	      nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorValsUC[plc + cnt] = CLIP_COLOR_VALUE(colorValsDA[cnt]);
	      }
	      nxtn = 4;
	    } else {
	      uchar col[4] = { CLIP_COLOR_VALUE(cgo->color[0]), CLIP_COLOR_VALUE(cgo->color[1]), CLIP_COLOR_VALUE(cgo->color[2]), CLIP_COLOR_VALUE(cgo->alpha) };
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorValsUC[plc + cnt] = col[cnt%4];
	      }
	    }
	  } else {
	    if (arrays & CGO_COLOR_ARRAY){
	      nxtVals = colorValsDA = nxtVals + (nxtn*nverts);
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorVals[plc + cnt] = colorValsDA[cnt];
	      }
	      nxtn = 4;
	    } else {
	      float col[4] = { cgo->color[0], cgo->color[1], cgo->color[2], cgo->alpha };
	      for (cnt=0; cnt<nverts*4; cnt++){
		colorVals[plc + cnt] = col[cnt%4];
	      }
	    }
	  }
	  if (arrays & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorValsDA = nxtVals + nverts;
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<nverts; cnt++){
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	      CGO_put_int(pickColorValsTMP++, CGO_get_int(pickColorValsDA++));
	    }
	    nxtn = 3;
	  } else {
	    pickColorValsTMP = pickColorVals + (vpl * 2);
	    for (cnt=0; cnt<nverts; cnt++){
	      CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_index);
	      CGO_put_int(pickColorValsTMP++, cgo->current_pick_color_bond);
	    }
	  }
	  if (idxpl + nverts > num_total_indexes_lines){
	    PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() num_total_indexes_lines=%d mode=%d nverts=%d idxpl=%d\n", num_total_indexes_lines, mode, nverts, idxpl ENDFB(I->G);
	  }
	  switch (mode){
	  case GL_LINES:
	    for (cnt = 0; cnt < nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_LINE_STRIP:
	    for (cnt = 1; cnt < nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt - 1;
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    break;
	  case GL_LINE_LOOP:
	    for (cnt = 1; cnt < nverts; cnt++){
	      vertexIndexes[idxpl++] = vpl + cnt - 1;
	      vertexIndexes[idxpl++] = vpl + cnt;
	    }
	    vertexIndexes[idxpl++] = vpl;
	    vertexIndexes[idxpl++] = vpl + nverts - 1;
	    break;
	  }

	  pl += 3 * nverts;
	  plc += 4 * nverts;
	  vpl += nverts;
	  pc += nvals + 4;
	  save_pc += nvals + 4 ;
	} else {
	  {
	    int nvals = narrays*nverts ;
	    pc += nvals + 4;
	    save_pc += nvals + 4 ;
	  }
	}
	}
	break;
#endif
      default:
	break;
      }
      pc = save_pc;
      pc += CGO_sz[op];
    }
    {
      uint bufs[5] = {0, 0, 0, 0, 0 }, allbufs[5] = { 0, 0, 0, 0, 0 };
      short bufpl = 0;
      GLenum err ;
      //      printf("CGOOptimizeToVBOIndexed: End: num_total_vertices=%d num_total_indexes=%d verts_skipped=%d\n", num_total_vertices, num_total_indexes, verts_skipped);

      glGenBuffers(4, bufs);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glGenBuffers returns err=%d\n", err ENDFB(I->G);	
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	allbufs[0] = bufs[bufpl++];
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices_lines*3, vertexVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 3;
	allbufs[1] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices_lines*sz, normalVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	short sz = 4;
	allbufs[2] = bufs[bufpl++];
	if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	  sz = 1;
	}
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*num_total_vertices_lines*sz, colorVals, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufs[bufpl]);
      if ((err = glGetError())){
	PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);	
      }
      if (!glIsBuffer(bufs[bufpl])){
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeToVBOIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);		  
      } else {
	allbufs[3] = bufs[bufpl++];
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint)*num_total_indexes_lines, vertexIndexes, GL_STATIC_DRAW);      
	if ((err = glGetError())){
	  PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeToVBOIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);	
	}
      }
      {
	GLfloat *newPickColorVals ;
#if defined(OPENGL_ES_2)
	CGODisable(cgo, GL_SHADER_LIGHTING);
#endif
	newPickColorVals = CGODrawBuffersIndexed(cgo, GL_LINES, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY | CGO_COLOR_ARRAY | CGO_PICK_COLOR_ARRAY, num_total_indexes_lines, num_total_vertices_lines, allbufs);
	memcpy(newPickColorVals + num_total_vertices_lines, pickColorVals, num_total_vertices_lines * 2 * sizeof(float));
	has_draw_buffer = true;
#if defined(OPENGL_ES_2)
	CGOEnable(cgo, GL_SHADER_LIGHTING);
#endif
      }
    }
    FreeP(vertexIndexes);
    FreeP(vertexVals);
  }
  if (num_total_vertices>0 || num_total_vertices_lines>0){
    CGOBoundingBox(cgo, min, max);
  }

  CGOStop(cgo);
  if (has_draw_buffer){
    cgo->has_draw_buffers = true;
  }
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  return (cgo);
}

#endif

CGO *CGOOptimizeGLSLCylindersToVBOIndexedImpl(CGO * I, int est, short no_color, CGO *leftOverCGO);

CGO *CGOOptimizeGLSLCylindersToVBOIndexedNoColor(CGO * I, int est){
  return (CGOOptimizeGLSLCylindersToVBOIndexedImpl(I, est, true, NULL));
}

CGO *CGOOptimizeGLSLCylindersToVBOIndexedWithLeftOver(CGO * I, int est, CGO *leftOverCGO){
  return (CGOOptimizeGLSLCylindersToVBOIndexedImpl(I, est, false, leftOverCGO));
}

CGO *CGOOptimizeGLSLCylindersToVBOIndexed(CGO * I, int est){
  return (CGOOptimizeGLSLCylindersToVBOIndexedImpl(I, est, false, NULL));
}

int CGOCountNumberCustomCylinders(CGO *I, int *has_2nd_color);

CGO *CGOOptimizeGLSLCylindersToVBOIndexedImpl(CGO * I, int est, short no_color, CGO *leftOverCGO)
{
  CGO *cgo = NULL;
  register float *pc = I->op;
  register int op;
  int sz;
  int box_indices[6 * 2 * 3] = { // box indices 
    0, 2, 1, 2, 0, 3, 1, 6, 5, 6, 1, 2, 0, 1, 5, 5, 4, 0, 
    0, 7, 3, 7, 0, 4, 3, 6, 2, 6, 3, 7, 4, 5, 6, 6, 7, 4 };
  int right_idx[8] =  { 0, 1, 1, 0, 0, 1, 1, 0 };
  int up_idx[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };
  int out_idx[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };
  short color2nd = 0, customCyl = 0;

  float *save_pc;
  int num_total_cylinders = 0, num_cylinders_with_2nd_color = 0, num_custom_cylinders = 0, num_custom_cylinders_with_2nd_color = 0;
  short err = 0;
  short has_draw_buffer = false;
  float min[3] = { MAXFLOAT, MAXFLOAT, MAXFLOAT }, max[3] = { -MAXFLOAT, -MAXFLOAT, -MAXFLOAT };
  int vv, total_vert = 0, total_cyl = 0;

  num_custom_cylinders = CGOCountNumberCustomCylinders(I, &num_custom_cylinders_with_2nd_color);
  num_cylinders_with_2nd_color = CGOCountNumberOfOperationsOfType(I, CGO_SHADER_CYLINDER_WITH_2ND_COLOR);
  num_total_cylinders = CGOCountNumberOfOperationsOfType(I, CGO_SHADER_CYLINDER) + num_cylinders_with_2nd_color + num_custom_cylinders;
  num_cylinders_with_2nd_color += num_custom_cylinders_with_2nd_color;
/* 
   The structure for cylinders is following:

   first vertex is a cylinder origin
   second vertex is cylinder axis vector
   third vertex is a corner flag (radius, right, up)
*/

  if (num_total_cylinders>0) {
    float *originVals = 0, *axisVals = 0;
    float *colorVals = 0, *color2Vals = 0;
    float axis[3], rad, col2[3] = { 0.f, 0.f, 0.f };
    int capvals;
    int *indexVals = 0;
    int tot = 4 * 4 * 3 * num_total_cylinders;
    short copyToLeftOver, copyColorToLeftOver, copyPickColorToLeftOver, copyAlphaToLeftOver ;
    float *org_originVals;
    float *org_axisVals;
    float *org_colorVals = NULL;
    float *org_color2Vals = NULL;
    int *org_indexVals;
    float min_alpha;

    cgo = CGONewSized(I->G, I->c + est);

    org_originVals = originVals = Alloc(float, tot);
    org_axisVals = axisVals = Alloc(float, tot);
    if (!no_color){
      org_colorVals = colorVals = Alloc(float, tot);
      if (num_cylinders_with_2nd_color){
	org_color2Vals = color2Vals = Alloc(float, tot);
      }
    }
    org_indexVals = indexVals = Alloc(int, tot);

    pc = I->op;
    cgo->alpha = 1.f;
    min_alpha = 1.f;
    copyToLeftOver = copyColorToLeftOver = copyPickColorToLeftOver = copyAlphaToLeftOver = 0;
    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      copyToLeftOver = false;
      save_pc = pc;
      err = 0;
      color2nd = false;
      customCyl = false;
      sz = -1;
      switch (op) {
      case CGO_NORMAL:
        cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
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
        cgo->current_pick_color_index = (int)pc[0];
        cgo->current_pick_color_bond = (int)pc[1];
	copyPickColorToLeftOver = true;
	break;
      case CGO_SAUSAGE:
      case CGO_CYLINDER:
      case CGO_CUSTOM_CYLINDER:
	customCyl = true;
	axis[0] = *(pc+3) - *(pc);
	axis[1] = *(pc+4) - *(pc+1);
	axis[2] = *(pc+5) - *(pc+2);
	if (op==CGO_CUSTOM_CYLINDER){
	  capvals = ((*(pc+13) > 1.5f) ? 5 : (*(pc+13) > 0.5f) ? 1 : 0) | ((*(pc+14) > 1.5f) ? 10 : (*(pc+14) > 0.5f) ? 2 : 0);
	} else if (op==CGO_CYLINDER) {
	  capvals = 3;
	} else {
	  capvals = 15;
	}
	cgo->color[0] = *(pc+7);
	cgo->color[1] = *(pc+8);
	cgo->color[2] = *(pc+9);
	if (*(pc+7) != *(pc+10) || *(pc+8) != *(pc+11) || *(pc+9) != *(pc+12)){
	  color2nd = true;
	  col2[0] = *(pc+10);
	  col2[1] = *(pc+11);
	  col2[2] = *(pc+12);
	}
      case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
	if (!customCyl){
	  color2nd = true;
	  col2[0] = *(pc+8);
	  col2[1] = *(pc+9);
	  col2[2] = *(pc+10);
	}
      case CGO_SHADER_CYLINDER:
	if (!customCyl){
	  axis[0] = *(pc+3);
	  axis[1] = *(pc+4);
	  axis[2] = *(pc+5);
	  capvals = CGO_get_int(pc+7);
	}
	rad = *(pc+6);
	for (vv=0; vv<8; vv++) { // generate eight vertices of a bounding box for each cylinder
	  originVals[0] = *(pc);
	  originVals[1] = *(pc+1);
	  originVals[2] = *(pc+2);
	  set_min_max(min, max, originVals);
	  axisVals[0] = axis[0];
	  axisVals[1] = axis[1];
	  axisVals[2] = axis[2];
	  originVals[3] = rad;
	  // pack the corner + cap flags into a single float
          // start cap = 1, end cap = 2
	  axisVals[3] = ((capvals << 18) |
			 (right_idx[vv] << 12) |
			 (up_idx[vv] << 6) |
			 (out_idx[vv]));
	  if (!no_color){
	    colorVals[0] = cgo->color[0];
	    colorVals[1] = cgo->color[1];
	    colorVals[2] = cgo->color[2];
	    colorVals[3] = cgo->alpha;
	    if (color2Vals) {
	      if (color2nd){
		color2Vals[0] = col2[0];
		color2Vals[1] = col2[1];
		color2Vals[2] = col2[2];
		color2Vals[3] = cgo->alpha;
	      } else {
		color2Vals[0] = cgo->color[0];
		color2Vals[1] = cgo->color[1];
		color2Vals[2] = cgo->color[2];
		color2Vals[3] = cgo->alpha;
	      }
	      color2Vals += 4;
	    }
	    colorVals += 4;
	  }
	  originVals += 4;
	  axisVals += 4;
	  total_vert++;
	}
	for (vv=0; vv<36; vv++) {
	  *(indexVals++) = box_indices[vv] + 8 * total_cyl;
	}
	total_cyl++;
	break;
#ifdef _PYMOL_CGO_DRAWARRAYS
      case CGO_DRAW_ARRAYS:
	{
	  int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	  int nvals = narrays*nverts, onvals;
	  onvals = nvals;
	  pc += 4;
	  save_pc += onvals + 4 ;
	  sz = onvals + 4 ;
	  copyToLeftOver = true;
	}
	break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeGLSLCylindersToVBO() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=0x%X\n", op ENDFB(I->G);	
	break;
#endif
      default:
	copyToLeftOver = true;
      }
      if (leftOverCGO && copyToLeftOver){
	float *npc = save_pc, *nc;
	if (copyAlphaToLeftOver){
	  CGOAlpha(leftOverCGO, cgo->alpha);
	}
	if (copyColorToLeftOver){
	  CGOColor(leftOverCGO, cgo->color[0],  cgo->color[1],  cgo->color[2] );
	}
	if (copyPickColorToLeftOver){
	  CGOPickColor(leftOverCGO, cgo->current_pick_color_index, cgo->current_pick_color_bond);
	}
	if (sz < 0){
	  sz = CGO_sz[op];
	} else {
	  npc -= sz;
	}
	nc = CGO_add(leftOverCGO, sz + 1);
	*(nc++) = *(npc - 1);
	while(sz--)
	  *(nc++) = *(npc++);
	copyToLeftOver = copyColorToLeftOver = copyPickColorToLeftOver = copyAlphaToLeftOver = 0;
      }
      pc = save_pc;
      pc += CGO_sz[op];
    }
    //    printf("total_cyl=%d total_vert=%d\n", total_cyl, total_vert);
    if (total_cyl > 0) {
      uint bufpl, bufs[5] = { 0, 0, 0, 0, 0};
      glGenBuffers(5, bufs);
      if ((err = glGetError())){
        PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeGLSLCylindersToVBO() glGenBuffers returns err=%d\n", err ENDFB(I->G);
      }

      for (bufpl=0; bufpl<5; bufpl++) {
        glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
        if ((err = glGetError())){
          PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeGLSLCylindersToVBO() glBindBuffer returns err=%d\n", err ENDFB(I->G);
        }

        if (!glIsBuffer(bufs[bufpl])){
          PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeGLSLCylindersToVBO() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);
        } else {
          switch(bufpl) {
            case 0: // midpoint
              glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert*4, org_originVals, GL_STATIC_DRAW);
              break;
            case 1: // axis
              glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert*4, org_axisVals, GL_STATIC_DRAW);
              break;
            case 2: // color
	      if (!no_color){
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert*4, org_colorVals, GL_STATIC_DRAW);
	      } else {
		if (bufs[2]){
		  CShaderMgr_AddVBOToFree(I->G->ShaderMgr, bufs[2]);
		  bufs[2] = 0;
		}
	      }
              break;
            case 3: // color2
	      if (!no_color && num_cylinders_with_2nd_color){
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert*4, org_color2Vals, GL_STATIC_DRAW);
	      } else {
		if (bufs[3]){
		  CShaderMgr_AddVBOToFree(I->G->ShaderMgr, bufs[3]);		
		  bufs[3] = 0;
		}
	      }
              break;
          }
          if ((err = glGetError())){
            PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeGLSLCylindersVBO() glBufferData returns err=%d\n", err ENDFB(I->G);
          }
        }
      }

      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufs[4]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int)*total_cyl*36, org_indexVals, GL_STATIC_DRAW);

      has_draw_buffer = true;
      CGODrawCylinderBuffers(cgo, total_cyl, (int)(255 * min_alpha), bufs);
    }

    FreeP(org_axisVals);
    FreeP(org_originVals);
    if (!no_color){
      FreeP(org_colorVals);
      if (org_color2Vals){
	FreeP(org_color2Vals);
      }
    }
    FreeP(org_indexVals);

    CGOBoundingBox(cgo, min, max);
    CGOStop(cgo);
    if (has_draw_buffer){
      cgo->has_draw_buffers = true;
      cgo->has_draw_cylinder_buffers = true;
    }
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader){
      cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
      cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
    }
  }
  return (cgo);
}

CGO *CGOOptimizeSpheresToVBONonIndexed(CGO * I, int est){
  return (CGOOptimizeSpheresToVBONonIndexedImpl(I, est, NULL));
}

CGO *CGOOptimizeSpheresToVBONonIndexedImpl(CGO * I, int est, CGO *leftOverCGO)
{
  CGO *cgo = NULL;
  register float *pc = I->op;
  register int op;
  int sz;

  int rightup_flags[4] = { 0, 1, 3, 2 };
  float *save_pc;
  int num_total_spheres = 0;
  short err = 0;
  short has_draw_buffer = false;
  float min[3] = { MAXFLOAT, MAXFLOAT, MAXFLOAT }, max[3] = { -MAXFLOAT, -MAXFLOAT, -MAXFLOAT };
  int vv, total_vert = 0, total_spheres = 0;

  num_total_spheres = CGOCountNumberOfOperationsOfType(I, CGO_SPHERE);
  if (num_total_spheres>0) {
    float *vertVals = 0;
    GLubyte *rightUpFlagValsUB = 0;
    float *rightUpFlagVals = 0;
    GLubyte *colorValsUB = 0;
    float *colorVals = 0;
    int tot = 4 * 4 * num_total_spheres;
    float *org_vertVals = NULL;
    float *org_colorVals = NULL;
    GLubyte *org_colorValsUB = NULL;
    GLubyte *org_rightUpFlagValsUB = NULL;
    float *org_rightUpFlagVals = NULL;
    float min_alpha;
    short cgo_shader_ub_color, cgo_shader_ub_flags;
    short copyToLeftOver, copyColorToLeftOver, copyPickColorToLeftOver, copyAlphaToLeftOver ;

    cgo = CGONewSized(I->G, I->c + est);

    cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
    cgo_shader_ub_flags = SettingGet(cgo->G, cSetting_cgo_shader_ub_flags);
  
    org_vertVals = vertVals = Alloc(float, tot);
    if(cgo_shader_ub_color){
      org_colorValsUB = colorValsUB = Alloc(GLubyte, tot);
    } else {
      org_colorVals = colorVals = Alloc(float, tot);
    }
    if (cgo_shader_ub_flags){
      org_rightUpFlagValsUB = rightUpFlagValsUB = Alloc(GLubyte, 4 * num_total_spheres);
    } else {
      org_rightUpFlagVals = rightUpFlagVals = Alloc(float, 4 * num_total_spheres);
    }
    pc = I->op;
    cgo->alpha = 1.f;
    min_alpha = 1.f;
    copyToLeftOver = copyColorToLeftOver = copyPickColorToLeftOver = copyAlphaToLeftOver = 0;
    while((op = (CGO_MASK & CGO_read_int(pc)))) {
      copyToLeftOver = false;
      save_pc = pc;
      err = 0;
      sz = -1;
      switch (op) {
      case CGO_NORMAL:
        cgo->normal[0] = *pc; cgo->normal[1] = *(pc + 1); cgo->normal[2] = *(pc + 2);
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
        cgo->current_pick_color_index = (int)pc[0];
        cgo->current_pick_color_bond = (int)pc[1];
	copyPickColorToLeftOver = true;
	break;
      case CGO_SPHERE:
	for (vv=0; vv<4; vv++) { // generate eight vertices of a bounding box for each cylinder
	  vertVals[0] = *(pc);
	  vertVals[1] = *(pc+1);
	  vertVals[2] = *(pc+2);
	  vertVals[3] = *(pc+3);
	  set_min_max(min, max, vertVals);
	  if (cgo_shader_ub_flags){
	    rightUpFlagValsUB[0] = rightup_flags[vv];
	    rightUpFlagValsUB++;
	  } else {
	    rightUpFlagVals[0] = rightup_flags[vv];
	    rightUpFlagVals++;
	  }
	  if (cgo_shader_ub_color){
	    colorValsUB[0] = CLIP_COLOR_VALUE(cgo->color[0]);
	    colorValsUB[1] = CLIP_COLOR_VALUE(cgo->color[1]);
	    colorValsUB[2] = CLIP_COLOR_VALUE(cgo->color[2]);
	    colorValsUB[3] = CLIP_COLOR_VALUE(cgo->alpha);
	    colorValsUB += 4;
	  } else {
	    colorVals[0] = cgo->color[0];
	    colorVals[1] = cgo->color[1];
	    colorVals[2] = cgo->color[2];
	    colorVals[3] = cgo->alpha;
	    colorVals += 4;
	  }
	  vertVals += 4;
	  total_vert++;
	}
	total_spheres++;
	break;
#ifdef _PYMOL_CGO_DRAWARRAYS
      case CGO_DRAW_ARRAYS:
	{
	  int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	  GLfloat *vals = CGODrawArrays(cgo, mode, arrays, nverts);
	  int nvals = narrays*nverts, onvals;
	  onvals = nvals;
	  pc += 4;
	  while(nvals--)
	    *(vals++) = *(pc++);
	  save_pc += onvals + 4 ;
	  sz = onvals + 4 ;
	  copyToLeftOver = true;
	}
	break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
      case CGO_DRAW_BUFFERS_INDEXED:
      case CGO_DRAW_BUFFERS_NOT_INDEXED:
	PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeSpheresToVBONonIndexed() CGO_DRAW_BUFFERS_INDEXED or CGO_DRAW_BUFFERS_INDEXED encountered op=%d\n", op ENDFB(I->G);	
	break;
#endif
      default:
	copyToLeftOver = true;
	sz = CGO_sz[op];
	pc += sz;
	/*	nc = CGO_add(cgo, sz + 1);
	*(nc++) = *(pc - 1);
	while(sz--)
	  *(nc++) = *(pc++);*/
      }
      if (leftOverCGO && copyToLeftOver){
	float *npc = pc, *nc;
	if (copyAlphaToLeftOver){
	  CGOAlpha(leftOverCGO, cgo->alpha);
	}
	if (copyColorToLeftOver){
	  CGOColor(leftOverCGO, cgo->color[0],  cgo->color[1],  cgo->color[2] );
	}
	if (copyPickColorToLeftOver){
	  CGOPickColor(leftOverCGO, cgo->current_pick_color_index, cgo->current_pick_color_bond);
	}
	if (sz < 0){
	  sz = CGO_sz[op];
	} else {
	  npc -= sz;
	}
	nc = CGO_add(leftOverCGO, sz + 1);
	*(nc++) = *(npc - 1);
	while(sz--)
	  *(nc++) = *(npc++);
	copyToLeftOver = copyColorToLeftOver = copyPickColorToLeftOver = copyAlphaToLeftOver = 0;
      }
      pc = save_pc;
      pc += CGO_sz[op];
    }

    if (total_spheres > 0) {
      uint bufpl, bufs[3] = { 0, 0, 0};
      glGenBuffers(5, bufs);
      if ((err = glGetError())){
        PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeSpheresToVBONonIndexed() glGenBuffers returns err=%d\n", err ENDFB(I->G);
      }

      for (bufpl=0; bufpl<3; bufpl++) {
        glBindBuffer(GL_ARRAY_BUFFER, bufs[bufpl]);
        if ((err = glGetError())){
          PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeSpheresToVBONonIndexed() glBindBuffer returns err=%d\n", err ENDFB(I->G);
        }

        if (!glIsBuffer(bufs[bufpl])){
          PRINTFB(I->G, FB_CGO, FB_Warnings) "WARNING: CGOOptimizeSpheresToVBONonIndexed() glGenBuffers created bad buffer bufpl=%d bufs[bufpl]=%d\n", bufpl, bufs[bufpl] ENDFB(I->G);
        } else {
          switch(bufpl) {
            case 0: // vertex xyz + right/up flag in w
              glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert*4, org_vertVals, GL_STATIC_DRAW);
              break;
            case 1: // color in UNSIGNED_BYTE
	      if (cgo_shader_ub_color){
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLubyte)*total_vert*4, org_colorValsUB, GL_STATIC_DRAW);
	      } else {
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert*4, org_colorVals, GL_STATIC_DRAW);
	      }
              break;
            case 2: // radius
	      if (cgo_shader_ub_flags){
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLubyte)*total_vert, org_rightUpFlagValsUB, GL_STATIC_DRAW);
	      } else {
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*total_vert, org_rightUpFlagVals, GL_STATIC_DRAW);
	      }
              break;
          }
          if ((err = glGetError())){
            PRINTFB(I->G, FB_CGO, FB_Errors) "ERROR: CGOOptimizeSpheresToVBONonIndexed() glBufferData returns err=%d\n", err ENDFB(I->G);
          }
        }
      }

      has_draw_buffer = true;

      CGODrawSphereBuffers(cgo, total_spheres, (cgo_shader_ub_color ? 1 : 0) | (cgo_shader_ub_flags ? 2 : 0), bufs);
    }

    FreeP(org_vertVals);
    if (cgo_shader_ub_color){
      FreeP(org_colorValsUB);
    } else {
      FreeP(org_colorVals);
    }
    if (cgo_shader_ub_flags){
      FreeP(org_rightUpFlagValsUB);
    } else {
      FreeP(org_rightUpFlagVals);
    }

    if (num_total_spheres>0){
      CGOBoundingBox(cgo, min, max);
    }

    CGOStop(cgo);
    
    if (has_draw_buffer){
      cgo->has_draw_buffers = true;
      cgo->has_draw_sphere_buffers = true;
    }
    cgo->use_shader = I->use_shader;
    if (cgo->use_shader){
      cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
      cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
    }
  }
  return (cgo);
}

CGO *CGOSimplify(CGO * I, int est)
{
  CGO *cgo;

  register float *pc = I->op;
  register float *nc;
  register int op;
  float *save_pc;
  int sz;

  cgo = CGONewSized(I->G, I->c + est);
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    switch (op) {
    case CGO_SHADER_CYLINDER:
      {
	float v2[3];
	add3f(pc, pc + 3, v2);
	CGOSimpleCylinder(cgo, pc, v2, *(pc + 6), 0, 0, 1, 1);
      }
      break;
    case CGO_SHADER_CYLINDER_WITH_2ND_COLOR:
      {
	float haxis[3], v1[3], v2[3];
	mult3f(pc + 3, .5f, haxis);
	add3f(pc, haxis, v1);
	CGOSimpleCylinder(cgo, pc, v1, *(pc + 6), 0, 0, 1, 0);
	
	add3f(v1, haxis, v2);
	CGOSimpleCylinder(cgo, v1, v2, *(pc + 6), pc+8, pc+8, 0, 1);
      }
      break;
    case CGO_CYLINDER:
      CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10, 1, 1);
      break;
    case CGO_CONE:
      CGOSimpleCone(cgo, pc, pc + 3, *(pc + 6), *(pc + 7), pc + 8, pc + 11,
                    (int) *(pc + 14), (int) *(pc + 15));
      break;
    case CGO_SAUSAGE:
      CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10, 2, 2);
      break;
    case CGO_CUSTOM_CYLINDER:
      CGOSimpleCylinder(cgo, pc, pc + 3, *(pc + 6), pc + 7, pc + 10, (int) *(pc + 13),
                        (int) *(pc + 14));
      break;
    case CGO_SPHERE:
      CGOSimpleSphere(cgo, pc, *(pc + 3));
      break;
    case CGO_ELLIPSOID:
      CGOSimpleEllipsoid(cgo, pc, *(pc + 3), pc + 4, pc + 7, pc + 10);
      break;
    case CGO_QUADRIC:
      CGOSimpleQuadric(cgo, pc, *(pc + 3), pc + 4);
      break;
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ;
	PRINTFB(I->G, FB_CGO, FB_Errors) "WARNING: CGOSimplify: CGO_DRAW_BUFFERS_INDEXED encountered nverts=%d\n", nverts ENDFB(I->G);
      }
      break;      
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ;
	PRINTFB(I->G, FB_CGO, FB_Errors) "WARNING: CGOSimplify: CGO_DRAW_BUFFERS_NOT_INDEXED encountered nverts=%d\n", nverts ENDFB(I->G);
      }
      break;      
#endif
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	GLfloat *vals = CGODrawArrays(cgo, mode, arrays, nverts);
	int nvals = narrays*nverts, onvals;
	onvals = nvals;
	pc += 4;
	while(nvals--)
	  *(vals++) = *(pc++);
	save_pc += onvals + 4 ;
      }
      break;
    case CGO_END:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOSimplify: CGO_END encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_VERTEX:
      PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOSimplify: CGO_VERTEX encountered without CGO_BEGIN but skipped for OpenGLES\n" ENDFB(I->G);      
      break;
    case CGO_BEGIN:
      {
	float *origpc = pc;
	int nverts = 0, damode = CGO_VERTEX_ARRAY, err = 0, end = 0;
	int mode = CGO_read_int(pc);

	while(!err && !end && (op = (CGO_MASK & CGO_read_int(pc)))) {
	  switch (op) {
	  case CGO_DRAW_ARRAYS:
	    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGOSimplify: CGO_DRAW_ARRAYS encountered inside CGO_BEGIN/CGO_END\n" ENDFB(I->G);
	    err = true;
	    continue;
	  case CGO_NORMAL:
	    damode |= CGO_NORMAL_ARRAY;
	    break;
	  case CGO_COLOR:
	    damode |= CGO_COLOR_ARRAY;
	    break;
	  case CGO_PICK_COLOR:
	    damode |= CGO_PICK_COLOR_ARRAY;
	    break;
	  case CGO_VERTEX:
	    nverts++;
	    break;
	  case CGO_END:
	    end = 1;
	  default:
	    break;
	  }
	  sz = CGO_sz[op];
	  pc += sz;
	}
	if (nverts>0 && !err){
	  int pl = 0, plc = 0, plp = 0;
	  float *vertexVals, *tmp_ptr;
	  float *normalVals = 0, *colorVals = 0, *nxtVals = 0, *pickColorVals = 0;
	  uchar *pickColorValsUC;
	  short notHaveValue = 0, nxtn = 3;
	  nxtVals = vertexVals = CGODrawArrays(cgo, mode, damode, nverts);	      
	  if (damode & CGO_NORMAL_ARRAY){
	    nxtVals = normalVals = vertexVals + (nxtn*nverts);
	  }
	  if (damode & CGO_COLOR_ARRAY){
	    nxtVals = colorVals = nxtVals + (nxtn*nverts);
	    nxtn = 4;
	  }
	  if (damode & CGO_PICK_COLOR_ARRAY){
	    nxtVals = nxtVals + (nxtn*nverts);
	    pickColorVals = nxtVals + nverts;
	    pickColorValsUC = (uchar*)nxtVals;
	    nxtn = 3;
	  }
	  pc = origpc + 1;
	  notHaveValue = damode;
	  end = 0;
	  while(!err && !end && (op = (CGO_MASK & CGO_read_int(pc)))) {
	    switch (op) {
	    case CGO_NORMAL:
	      normalVals[pl] = pc[0]; normalVals[pl+1] = pc[1]; normalVals[pl+2] = pc[2];
	      notHaveValue = notHaveValue ^ CGO_NORMAL_ARRAY;
	      break;
	    case CGO_COLOR:
	      colorVals[plc] = pc[0]; colorVals[plc+1] = pc[1]; 
	      colorVals[plc+2] = pc[2]; colorVals[plc+3] = I->alpha;
	      notHaveValue = notHaveValue ^ CGO_COLOR_ARRAY;
	      break;
	    case CGO_PICK_COLOR:
	      cgo->current_pick_color_index = (int)pc[0];
	      cgo->current_pick_color_bond = (int)pc[1];
	      notHaveValue = notHaveValue ^ CGO_PICK_COLOR_ARRAY;
	      break;
	    case CGO_VERTEX:
	      if (notHaveValue & CGO_NORMAL_ARRAY){
		tmp_ptr = &normalVals[pl-3];
		normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];		
	      }
	      if (notHaveValue & CGO_COLOR_ARRAY){
		tmp_ptr = &colorVals[plc-4];
		colorVals[plc] = tmp_ptr[0]; colorVals[plc+1] = tmp_ptr[1]; 
		colorVals[plc+2] = tmp_ptr[2]; colorVals[plc+3] = tmp_ptr[3]; 
	      }
	      if (notHaveValue & CGO_PICK_COLOR_ARRAY){
		CGO_put_int(pickColorVals + plp, cgo->current_pick_color_index);
		plp++;
		CGO_put_int(pickColorVals + plp, cgo->current_pick_color_bond);
		plp++;
	      }
	      vertexVals[pl++] = pc[0]; vertexVals[pl++] = pc[1]; vertexVals[pl++] = pc[2];
	      plc += 4;
	      notHaveValue = damode;
	      break;
	    case CGO_END:
	      end = 1;
	    default:
	      break;
	    }
	    sz = CGO_sz[op];
	    pc += sz;
	  }
	  save_pc = pc;
	} else {
	  save_pc = origpc;
	}
	op = CGO_NULL;
      }
      break;
#endif
    case CGO_ALPHA:
      I->alpha = *pc;
    default:
      sz = CGO_sz[op];
      nc = CGO_add(cgo, sz + 1);
      *(nc++) = *(pc - 1);
      while(sz--)
        *(nc++) = *(pc++);
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
  CGOStop(cgo);
  return (cgo);
}


/* ======== Raytrace Renderer ======== */

int CGOGetExtent(CGO * I, float *mn, float *mx)
{
  register float *pc = I->op;
  register int op;
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

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
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
      check_extent(pc, *(pc + 6));
      check_extent(pc + 3, *(pc + 6));
      break;
    case CGO_TRIANGLE:
      check_extent(pc, 0);
      check_extent(pc + 3, 0);
      check_extent(pc + 6, 0);
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), pl;
	float *pct = pc + 4;
	int nvals = narrays*nverts;

	if (arrays & CGO_VERTEX_ARRAY){
	  for (pl = 0; pl < nverts; pl++){
	    check_extent(pct, 0);
	    pct += 3;
	  }
	}
	if (arrays & CGO_NORMAL_ARRAY){
	  for (pl = 0; pl < nverts; pl++){
	    pct += 3;
	  }
	}
	if (arrays & CGO_COLOR_ARRAY){
	  for (pl = 0; pl < nverts; pl++){
	    pct += 4;
	  }
	}
	if (arrays & CGO_PICK_COLOR_ARRAY){
	  for (pl = 0; pl < nverts; pl++){
	    pct += 3;
	  }
	}
	pc += nvals + 4;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ;
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ;
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
#endif
    }
    pc += CGO_sz[op];
  }
  return (result);
}

int CGOHasNormals(CGO * I)
{
  register float *pc = I->op;
  register int op;
  int result = false;

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    switch (op) {
    case CGO_NORMAL:
    case CGO_SPHERE:
    case CGO_ELLIPSOID:
    case CGO_CYLINDER:
    case CGO_CONE:
    case CGO_SAUSAGE:
    case CGO_CUSTOM_CYLINDER:
      result |= 1;
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	int nvals = narrays*nverts;
	if (arrays & CGO_NORMAL_ARRAY){
	  result |= 1;
	}
	pc += nvals + 4;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ;
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ;
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (result);
}

static int CGOQuadricToEllipsoid(float *v, float r, float *q,
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

static void CGORenderQuadricRay(CRay * ray, float *v, float r, float *q)
{
  float r_el, n0[3], n1[3], n2[3];
  if(CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    ray->fEllipsoid3fv(ray, v, r_el, n0, n1, n2);

}


/* ======== Raytrace Renderer ======== */

void CGORenderRay(CGO * I, CRay * ray, float *color, CSetting * set1, CSetting * set2)
{
  register float *pc;
  register int op;
  int vc = 0;
  float linewidth = 1.0F;
  float widthscale = 0.15F;
  float lineradius, dotradius, dotwidth;
  float white[] = { 1.0, 1.0, 1.0 };
  float zee[] = { 0.0, 0.0, 1.0 };
  
  float *n0 = NULL, *n1 = NULL, *n2 = NULL, *v0 = NULL, *v1 = NULL, *v2 = NULL, *c0 =
    NULL, *c1 = NULL, *c2 = NULL;
  int mode = -1;
  /* workaround; multi-state ray-trace bug */
  if (I)
    pc = I->op;
  else return;

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
  ray->fTransparentf(ray, 1.0F - I->G->CGORenderer->alpha);

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
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
          ray->fSausage3fv(ray, v0, v2, lineradius, c0, c2);
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
    case CGO_COLOR:
      c0 = pc;
      ray->fColor3fv(ray, c0);
      break;
    case CGO_ALPHA:
      I->G->CGORenderer->alpha = *pc;
      ray->fTransparentf(ray, 1.0F - *pc);
      break;
    case CGO_VERTEX:
      v0 = pc;
      switch (mode) {
      case GL_POINTS:
        ray->fSphere3fv(ray, v0, dotradius);
        break;
      case GL_LINES:
        if(vc & 0x1)
          ray->fSausage3fv(ray, v0, v1, lineradius, c0, c1);
        v1 = v0;
        c1 = c0;
        break;
      case GL_LINE_STRIP:
        if(vc)
          ray->fSausage3fv(ray, v0, v1, lineradius, c0, c1);
        v1 = v0;
        c1 = c0;
        break;
      case GL_LINE_LOOP:
        if(vc)
          ray->fSausage3fv(ray, v0, v1, lineradius, c0, c1);
        else {
          v2 = v0;
          c2 = c0;
        }
        v1 = v0;
        c1 = c0;
        break;
      case GL_TRIANGLES:
	if( ((vc + 1) % 3) == 0)
          ray->fTriangle3fv(ray, v0, v1, v2, n0, n1, n2, c0, c1, c2);
        v2 = v1;
        c2 = c1;
        n2 = n1;
        v1 = v0;
        c1 = c0;
        n1 = n0;
        break;
      case GL_TRIANGLE_STRIP:
        if(vc > 1)
          ray->fTriangle3fv(ray, v0, v1, v2, n0, n1, n2, c0, c1, c2);
        v2 = v1;
        c2 = c1;
        n2 = n1;
        v1 = v0;
        c1 = c0;
        n1 = n0;
        break;
      case GL_TRIANGLE_FAN:
        if(vc > 1)
          ray->fTriangle3fv(ray, v0, v1, v2, n0, n1, n2, c0, c1, c2);
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
      break;
    case CGO_SPHERE:
      ray->fColor3fv(ray, c0);
      ray->fSphere3fv(ray, pc, *(pc + 3));
      break;
    case CGO_ELLIPSOID:
      ray->fColor3fv(ray, c0);
      ray->fEllipsoid3fv(ray, pc, *(pc + 3), pc + 4, pc + 7, pc + 10);
      break;
    case CGO_QUADRIC:
      ray->fColor3fv(ray, c0);
      CGORenderQuadricRay(ray, pc, *(pc + 3), pc + 4);
      break;
    case CGO_CONE:
      ray->fCone3fv(ray, pc, pc + 3, *(pc + 6), *(pc + 7), pc + 8, pc + 11,
                    (int) *(pc + 14), (int) *(pc + 15));
      break;
    case CGO_CUSTOM_CYLINDER:
      ray->fCustomCylinder3fv(ray, pc, pc + 3, *(pc + 6), pc + 7, pc + 10,
			      (int) *(pc + 13), (int) *(pc + 14));
      break;
    case CGO_CYLINDER:
      ray->fCylinder3fv(ray, pc, pc + 3, *(pc + 6), pc + 7, pc + 10);
      break;
    case CGO_SAUSAGE:
      ray->fSausage3fv(ray, pc, pc + 3, *(pc + 6), pc + 7, pc + 10);
      break;
    case CGO_TRIANGLE:
      ray->fTriangle3fv(ray, pc, pc + 3, pc + 6, pc + 9, pc + 12, pc + 15, pc + 18,
                        pc + 21, pc + 24);
      break;
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int mode = CGO_read_int(pc), arrays = CGO_read_int(pc), narrays = CGO_read_int(pc), nverts = CGO_read_int(pc), v, pl, plc;
	float *vertexVals = NULL;
	float *normalVals = 0, *colorVals = 0, *pickColorVals;
	(void)narrays;
	if (arrays & CGO_VERTEX_ARRAY){
	  vertexVals = pc; pc += nverts * 3;
	}
	if (arrays & CGO_NORMAL_ARRAY){
	  normalVals = pc; pc += nverts * 3;
	}	
	if (arrays & CGO_COLOR_ARRAY){
	  colorVals = pc; pc += nverts * 4;
	}
	if (arrays & CGO_PICK_COLOR_ARRAY){
	  pickColorVals = pc ; pc += nverts * 3;
	}
	vc = 0;
	for (v=0, pl=0, plc=0; v<nverts; v++, pl+=3, plc+=4){
	  if (normalVals){
	    n0 = &normalVals[pl];
	  }
	  if (colorVals){
	    c0 = &colorVals[plc];
	    ray->fColor3fv(ray, c0);
	  }
	  if (vertexVals){
	    v0 = &vertexVals[pl];
	  }
	  switch (mode){
	  case GL_POINTS:
	    ray->fSphere3fv(ray, v0, dotradius);
	    break;
	  case GL_LINES:
	    if(vc & 0x1)
	      ray->fSausage3fv(ray, v0, v1, lineradius, c0, c1);
	    v1 = v0;
	    c1 = c0;
	    break;
	  case GL_LINE_STRIP:
	    if(vc)
	      ray->fSausage3fv(ray, v0, v1, lineradius, c0, c1);
	    v1 = v0;
	    c1 = c0;
	    break;
	  case GL_LINE_LOOP:
	    if(vc)
	      ray->fSausage3fv(ray, v0, v1, lineradius, c0, c1);
	    else {
	      v2 = v0;
	      c2 = c0;
	    }
	    v1 = v0;
	    c1 = c0;
	    break;
	  case GL_TRIANGLES:
	    if( ((vc + 1) % 3) == 0)
	      ray->fTriangle3fv(ray, v0, v1, v2, n0, n1, n2, c0, c1, c2);
	    v2 = v1;
	    c2 = c1;
	    n2 = n1;
	    v1 = v0;
	    c1 = c0;
	    n1 = n0;
	    break;
	  case GL_TRIANGLE_STRIP:
	    if(vc > 1)
	      ray->fTriangle3fv(ray, v0, v1, v2, n0, n1, n2, c0, c1, c2);
	    v2 = v1;
	    c2 = c1;
	    n2 = n1;
	    v1 = v0;
	    c1 = c0;
	    n1 = n0;
	    break;
	  case GL_TRIANGLE_FAN:
	    if(vc > 1)
	      ray->fTriangle3fv(ray, v0, v1, v2, n0, n1, n2, c0, c1, c2);
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
#endif
    default:
      break;
    }
    pc += CGO_sz[op];
  }

  ray->fTransparentf(ray, 0.0F);
}


/* ======== GL Rendering ======== */

#ifdef _PYMOL_CGO_DRAWARRAYS
static int CGO_gl_begin_WARNING_CALLED = false, CGO_gl_end_WARNING_CALLED = false, CGO_gl_vertex_WARNING_CALLED = false;
static void CGO_gl_begin(CCGORenderer * I, float **pc){ 
  if (!CGO_gl_begin_WARNING_CALLED) { 
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_begin() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_begin_WARNING_CALLED = true; 
  }
}
static void CGO_gl_end(CCGORenderer * I, float **pc){ 
  if (!CGO_gl_end_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_end() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_end_WARNING_CALLED = true; 
  }
}
static void CGO_gl_vertex(CCGORenderer * I, float **v){
  if (!CGO_gl_vertex_WARNING_CALLED) {
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_vertex() is called but not implemented in OpenGLES\n" ENDFB(I->G);
    CGO_gl_vertex_WARNING_CALLED = true;
  }
}
static void CGO_gl_normal(CCGORenderer * I, float **varg){
  float *v = *varg;
#ifdef OPENGL_ES_2
  if (I->use_shader){
    glVertexAttrib3fv(VERTEX_NORMAL, v);
  } else {
#ifndef PURE_OPENGL_ES_2
    glNormal3f(v[0],v[1],v[2]);
#endif
  }
#else
  glNormal3f(v[0],v[1],v[2]);
#endif
}
#else
static void CGO_gl_begin(CCGORenderer * I, float **pc)
{
  glBegin(CGO_get_int(*pc));
}

static void CGO_gl_end(CCGORenderer * I, float **pc)
{
  glEnd();
}

static void CGO_gl_vertex(CCGORenderer * I, float **v)
{
  glVertex3fv(*v);
}

static void CGO_gl_normal(CCGORenderer * I, float **v)
{
#ifdef OPENGL_ES_2
  if (I->use_shader){
    glVertexAttrib3fv(VERTEX_NORMAL, v);
  } else {
    glNormal3fv(*v);
  }
#else
  glNormal3fv(*v);
#endif
}

#endif

static void CGO_gl_draw_arrays(CCGORenderer * I, float **pc){
  int mode = CGO_read_int(*pc), arrays = CGO_read_int(*pc), narrays = CGO_read_int(*pc), nverts = CGO_read_int(*pc);
  //  CShaderPrg * shaderPrg;
  /*  if (I->enable_shaders){
    shaderPrg = CShaderPrg_Enable_DefaultShader(I->G);
    }*/
  (void) narrays;
#if !defined(_PYMOL_PURE_OPENGL_ES) && (defined(OPENGL_ES_1) || defined(OPENGL_ES_2))
  if (I->use_shader){
#endif

#if defined(OPENGL_ES_1) || defined(OPENGL_ES_2)
#ifdef OPENGL_ES_1
  if (arrays & CGO_VERTEX_ARRAY) glEnableClientState(GL_VERTEX_ARRAY);
  if (arrays & CGO_NORMAL_ARRAY) glEnableClientState(GL_NORMAL_ARRAY);
  if (I->isPicking){
    if (arrays & CGO_PICK_COLOR_ARRAY){
      glEnableClientState(GL_COLOR_ARRAY);
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY)
      glEnableClientState(GL_COLOR_ARRAY);
  }
#else
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
#endif

  if (arrays & CGO_VERTEX_ARRAY){
#ifdef OPENGL_ES_1
    glVertexPointer(3, GL_FLOAT, 0, *pc);
#else
    glVertexAttribPointer(VERTEX_POS, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, *pc);
#endif
    *pc += nverts*3;
  }
  if (arrays & CGO_NORMAL_ARRAY){
#ifdef OPENGL_ES_1
    glNormalPointer(GL_FLOAT, 0, *pc);
#else
    glVertexAttribPointer(VERTEX_NORMAL, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, *pc);
#endif
    *pc += nverts*3;
  }
  if (I->isPicking){
    if (arrays & CGO_COLOR_ARRAY){
      *pc += nverts*4;
    }
    if (arrays & CGO_PICK_COLOR_ARRAY){
#ifdef OPENGL_ES_1
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, *pc);
#else
      glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_FALSE, 0, *pc);      
#endif
      *pc += nverts*3;
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY){
#ifdef OPENGL_ES_1
      glColorPointer(4, GL_FLOAT, 0, *pc);
#else
      glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, *pc);
#endif
      *pc += nverts*4;
    }
    if (arrays & CGO_PICK_COLOR_ARRAY){
      *pc += nverts*3;
    }
  }
  if (I->debug){
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
  }
  glDrawArrays(mode, 0, nverts);

#ifdef OPENGL_ES_1
  if (I->isPicking){
    if (arrays & CGO_PICK_COLOR_ARRAY){
      glDisableClientState(GL_COLOR_ARRAY);
    }
  } else {
    if (arrays & CGO_COLOR_ARRAY)
      glDisableClientState(GL_COLOR_ARRAY);
  }
  if (arrays & CGO_VERTEX_ARRAY) glDisableClientState(GL_VERTEX_ARRAY);
  if (arrays & CGO_NORMAL_ARRAY) glDisableClientState(GL_NORMAL_ARRAY);
#else
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
#endif
#ifndef _PYMOL_PURE_OPENGL_ES
  } else {
#endif
#endif

#ifndef _PYMOL_PURE_OPENGL_ES
    int pl, pla, plc;
    float *vertexVals = NULL;
    float *colorVals = 0, *normalVals = 0, *tmp_ptr, alpha ;
    uchar *pickColorVals = 0, *tmp_pc_ptr;
    alpha = I->alpha;
    if (arrays & CGO_VERTEX_ARRAY){
      vertexVals = *pc;
      *pc += nverts*3;
    }
    if (arrays & CGO_NORMAL_ARRAY){
      normalVals = *pc;
      *pc += nverts*3;
    }
    if (I->isPicking){
      alpha = 1.f;
      if (arrays & CGO_COLOR_ARRAY){
	*pc += nverts*4;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY){
	pickColorVals = (uchar*)*pc;
	*pc += nverts*3;
      }
    } else {
      if (arrays & CGO_COLOR_ARRAY){
	colorVals = *pc;
	*pc += nverts*4;
      }
      if (arrays & CGO_PICK_COLOR_ARRAY){
	*pc += nverts*3;
      }
    }
    if (arrays & CGO_ACCESSIBILITY_ARRAY) *pc += nverts;

    glBegin(mode);
    for (pl = 0, pla = 0, plc = 0; pl<nverts; pl++, pla+=3, plc+=4){
      if (colorVals){
	tmp_ptr = &colorVals[plc];
	glColor4f(tmp_ptr[0], tmp_ptr[1], tmp_ptr[2], alpha);	
      }
      if (pickColorVals){
	tmp_pc_ptr = &pickColorVals[plc]; /* the pick colors are saved with rgba */
	glColor3ub(tmp_pc_ptr[0], tmp_pc_ptr[1], tmp_pc_ptr[2]);	
      }
      if (normalVals){
	tmp_ptr = &normalVals[pla];
	glNormal3fv(&normalVals[pla]);
      }
      if (vertexVals){
	tmp_ptr = &vertexVals[pla];
	glVertex3fv(&vertexVals[pla]);
      }
    }
    glEnd();
  }
  /*  if (I->enable_shaders){
    CShaderPrg_Disable(shaderPrg);
    }*/
#endif
}

static void CGO_gl_draw_buffers(CCGORenderer * I, float **pc){
#ifdef _PYMOL_CGO_DRAWBUFFERS
  int mode = CGO_get_int(*pc), arrays = CGO_get_int(*pc+1), narrays = CGO_get_int(*pc+2), nverts = CGO_get_int(*pc+3);
  uint bufs[4] = { CGO_get_int(*pc+4), CGO_get_int(*pc+5), CGO_get_int(*pc+6), CGO_get_int(*pc+7) };
  CShaderPrg * shaderPrg = NULL;
  if (I->enable_shaders){
    shaderPrg = CShaderPrg_Enable_DefaultShader(I->G);
  }
  (void) narrays;
  (void) arrays;
  if (bufs[0]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[0]);
#ifdef OPENGL_ES_1
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);
#else
    glEnableVertexAttribArray(VERTEX_POS);
    glVertexAttribPointer(VERTEX_POS, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#endif
  }
  if (bufs[1]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[1]);
#ifdef OPENGL_ES_1
    glNormalPointer(GL_FLOAT, 0, 0);
    glEnableClientState(GL_NORMAL_ARRAY);
#else
    glEnableVertexAttribArray(VERTEX_NORMAL);
    glVertexAttribPointer(VERTEX_NORMAL, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#endif
  }
  if ((I->isPicking && bufs[3])){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[3]);
#ifdef OPENGL_ES_1
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
    glEnableClientState(GL_COLOR_ARRAY);
#else
    glEnableVertexAttribArray(VERTEX_COLOR);
    glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
#endif
  } else if (bufs[2]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[2]);
#ifdef OPENGL_ES_1
    glColorPointer(4, GL_FLOAT, 0, 0);
    glEnableClientState(GL_COLOR_ARRAY);
#else
    glEnableVertexAttribArray(VERTEX_COLOR);
    glVertexAttribPointer(VERTEX_COLOR, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#endif
  }

  if (I->debug){
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
  }
  glDrawArrays(mode, 0, nverts);

#ifdef OPENGL_ES_1
  if (bufs[0]) glDisableClientState(GL_VERTEX_ARRAY);
  if (bufs[1]) glDisableClientState(GL_NORMAL_ARRAY);
  if ((I->isPicking && bufs[3])){
    glDisableClientState(GL_COLOR_ARRAY);
  } else if (bufs[2]){
    glDisableClientState(GL_COLOR_ARRAY);
  }
#else
  if (bufs[0]) glDisableVertexAttribArray(VERTEX_POS);
  if (bufs[1]) glDisableVertexAttribArray(VERTEX_NORMAL);
  if ((I->isPicking && bufs[3])){
    glDisableVertexAttribArray(VERTEX_COLOR);
  } else if (bufs[2]){
    glDisableVertexAttribArray(VERTEX_COLOR);
  }
#endif
  //  printf("CGO_gl_draw_buffers called bufs: %d %d %d %d\n", bufs[0], bufs[1], bufs[2], bufs[3]);
  if (I->enable_shaders){
    CShaderPrg_Disable(shaderPrg);
  }
#endif
}

static void CGO_gl_draw_buffers_indexed(CCGORenderer * I, float **pc){
#ifdef _PYMOL_CGO_DRAWBUFFERS
  int mode = CGO_get_int(*pc), arrays = CGO_get_int(*pc+1), narrays = CGO_get_int(*pc+2), nindices = CGO_get_int(*pc+3), 
    nverts = CGO_get_int(*pc+4);
  uint bufs[5] = { CGO_get_int(*pc+5), CGO_get_int(*pc+6), CGO_get_int(*pc+7), CGO_get_int(*pc+8), CGO_get_int(*pc+9) };
  CShaderPrg * shaderPrg;
  int attr_a_Vertex, attr_a_Normal, attr_a_Color, attr_a_Accessibility;
  if (I->enable_shaders){
    shaderPrg = CShaderPrg_Enable_DefaultShader(I->G);
  } else {
    shaderPrg = CShaderMgr_GetShaderPrg(I->G->ShaderMgr, "default");
  }
  attr_a_Vertex = CShaderPrg_GetAttribLocation(shaderPrg, "a_Vertex");
  attr_a_Normal = CShaderPrg_GetAttribLocation(shaderPrg, "a_Normal"); 
  attr_a_Color = CShaderPrg_GetAttribLocation(shaderPrg, "a_Color");
  attr_a_Accessibility = CShaderPrg_GetAttribLocation(shaderPrg, "a_Accessibility");
  (void) arrays;
  (void) narrays;
  if (bufs[0]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[0]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Vertex);
    glVertexAttribPointer(attr_a_Vertex, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Vertex);
      glVertexAttribPointer(attr_a_Vertex, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    } else {
      glVertexPointer(3, GL_FLOAT, 0, 0);
      glEnableClientState(GL_VERTEX_ARRAY);
    }
#else
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);
#endif
#endif
  }
  if (bufs[1]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[1]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Normal);
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_BYTE, GL_TRUE, 0, 0);
    } else {
      glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    }
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Normal);
      if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_BYTE, GL_TRUE, 0, 0);
      } else {
	glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
      }
    } else {
      if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	glNormalPointer(GL_BYTE, 0, 0);
      } else {
	glNormalPointer(GL_FLOAT, 0, 0);
      }
      glEnableClientState(GL_NORMAL_ARRAY);
    }
#else
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      glNormalPointer(GL_BYTE, 0, 0);
    } else {
      glNormalPointer(GL_FLOAT, 0, 0);
    }
    glEnableClientState(GL_NORMAL_ARRAY);
#endif
#endif
  }
  if (I->isPicking){
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Color);
    glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, *pc + 9);
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Color);
      glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, *pc + 9);
    } else {
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, *pc + 9);
      glEnableClientState(GL_COLOR_ARRAY);
    }
#else
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, *pc + 9);
    glEnableClientState(GL_COLOR_ARRAY);
#endif
#endif
  } else if (bufs[2]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[2]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Color);
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
    } else {
      glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    }
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Color);
      if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
      } else {
	glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
      }
    } else {
      if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
      } else {
	glColorPointer(4, GL_FLOAT, 0, 0);
      }
      glEnableClientState(GL_COLOR_ARRAY);
    }
#else
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
    } else {
      glColorPointer(4, GL_FLOAT, 0, 0);
    }
    glEnableClientState(GL_COLOR_ARRAY);
#endif
#endif
  }
  {
    if (bufs[4]){
      glBindBuffer(GL_ARRAY_BUFFER, bufs[4]);
#ifdef PURE_OPENGL_ES_2
      glEnableVertexAttribArray(attr_a_Accessibility);
      glVertexAttribPointer(attr_a_Accessibility, 1, GL_FLOAT, GL_FALSE, 0, 0);
#else
#ifdef OPENGL_ES_2
      if (I->use_shader){
	glEnableVertexAttribArray(attr_a_Accessibility);
	glVertexAttribPointer(attr_a_Accessibility, 1, GL_FLOAT, GL_FALSE, 0, 0);
      } else {
	glVertexPointer(1, GL_FLOAT, 0, 0);
	glEnableClientState(GL_VERTEX_ARRAY);
      }
#else
      glVertexPointer(1, GL_FLOAT, 0, 0);
      glEnableClientState(GL_VERTEX_ARRAY);
#endif
#endif
    } else {
      glVertexAttrib1f(attr_a_Accessibility, 1.f);
    }
  }

  if (bufs[3]){
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufs[3]);
  }
  if (I->debug){
    mode = CGOConvertDebugMode(I->debug, mode);
  }

  glDrawElements(mode, nindices, GL_UNSIGNED_INT, 0);

#ifdef PURE_OPENGL_ES_2
  if (bufs[3]) glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  if (bufs[0]) glDisableVertexAttribArray(attr_a_Vertex);
  if (bufs[1]) glDisableVertexAttribArray(attr_a_Normal);
  if (I->isPicking){
    glDisableVertexAttribArray(attr_a_Color);
  } else if (bufs[2]){
    glDisableVertexAttribArray(attr_a_Color);
  }
#else
#ifdef OPENGL_ES_2
  if (I->use_shader){
    if (bufs[3]) glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    if (bufs[4]) glDisableVertexAttribArray(attr_a_Accessibility);
    if (bufs[0]) glDisableVertexAttribArray(attr_a_Vertex);
    if (bufs[1]) glDisableVertexAttribArray(attr_a_Normal);
    if (I->isPicking){
      glDisableVertexAttribArray(attr_a_Color);
    } else if (bufs[2]){
      glDisableVertexAttribArray(attr_a_Color);
    }
  } else {
    if (bufs[3]) glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    if (bufs[4]) glDisableClientState(attr_a_Accessibility);
    if (bufs[0]) glDisableClientState(GL_VERTEX_ARRAY);
    if (bufs[1]) glDisableClientState(GL_NORMAL_ARRAY);
    if (I->isPicking){
      glDisableClientState(GL_COLOR_ARRAY);
    } else if (bufs[2]){
      glDisableClientState(GL_COLOR_ARRAY);
    }
  }
#else
  if (bufs[3]) glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  if (bufs[0]) glDisableVertexAttribArray(attr_a_Vertex);
  if (bufs[1]) glDisableVertexAttribArray(attr_a_Normal);
  if (bufs[4]) glDisableVertexAttribArray(attr_a_Accessibility);
  if (I->isPicking){
    glDisableVertexAttribArray(attr_a_Color);
  } else if (bufs[2]){
    glDisableVertexAttribArray(attr_a_Color);
  }
#endif
#endif
  *pc += nverts*3 + 10;
  if (I->enable_shaders){
    CShaderPrg_Disable(shaderPrg);
  }
#endif
}

static void CGO_gl_draw_buffers_not_indexed(CCGORenderer * I, float **pc){
#ifdef _PYMOL_CGO_DRAWBUFFERS
  int mode = CGO_get_int(*pc), arrays = CGO_get_int(*pc+1), narrays = CGO_get_int(*pc+2),
    nverts = CGO_get_int(*pc+3);
  uint bufs[4] = { CGO_get_int(*pc+4), CGO_get_int(*pc+5), CGO_get_int(*pc+6), CGO_get_int(*pc+7) };
  CShaderPrg * shaderPrg;
  int attr_a_Vertex, attr_a_Normal, attr_a_Color, attr_a_Accessibility;
  if (I->enable_shaders){
    shaderPrg = CShaderPrg_Enable_DefaultShader(I->G);
  } else {
    shaderPrg = CShaderMgr_GetShaderPrg(I->G->ShaderMgr, "default");
  }
  attr_a_Vertex = CShaderPrg_GetAttribLocation(shaderPrg, "a_Vertex");
  attr_a_Normal = CShaderPrg_GetAttribLocation(shaderPrg, "a_Normal"); 
  attr_a_Color = CShaderPrg_GetAttribLocation(shaderPrg, "a_Color");
  attr_a_Accessibility = CShaderPrg_GetAttribLocation(shaderPrg, "a_Accessibility");

  (void) arrays;
  (void) narrays;
  if (bufs[0]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[0]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Vertex);
    glVertexAttribPointer(attr_a_Vertex, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Vertex);
      glVertexAttribPointer(attr_a_Vertex, VERTEX_POS_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    } else {
      glVertexPointer(3, GL_FLOAT, 0, 0);
      glEnableClientState(GL_VERTEX_ARRAY);
    }
#else
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);
#endif
#endif
  }
  if (bufs[1]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[1]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Normal);
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_BYTE, GL_TRUE, 0, 0);
    } else {
      glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    }
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Normal);
      if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_BYTE, GL_TRUE, 0, 0);
      } else {
	glVertexAttribPointer(attr_a_Normal, VERTEX_NORMAL_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
      }
    } else {
      if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
	glNormalPointer(GL_BYTE, 0, 0);
      } else {
	glNormalPointer(GL_FLOAT, 0, 0);
      }
      glEnableClientState(GL_NORMAL_ARRAY);
    }
#else
    if (SettingGet(I->G, cSetting_cgo_shader_ub_normal)){
      glNormalPointer(GL_BYTE, 0, 0);
    } else {
      glNormalPointer(GL_FLOAT, 0, 0);
    }
    glEnableClientState(GL_NORMAL_ARRAY);
#endif
#endif
  }
  if (I->isPicking){
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Color);
    glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, *pc + 7);
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Color);
      glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, *pc + 7);
    } else {
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, *pc + 7);
      glEnableClientState(GL_COLOR_ARRAY);
    }
#else
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, *pc + 7);
    glEnableClientState(GL_COLOR_ARRAY);
#endif
#endif
  } else if (bufs[2]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[2]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Color);
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
    } else {
      glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
    }
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Color);
      if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
      } else {
	glVertexAttribPointer(attr_a_Color, VERTEX_COLOR_SIZE, GL_FLOAT, GL_FALSE, 0, 0);
      }
    } else {
      if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
      } else {
	glColorPointer(4, GL_FLOAT, 0, 0);
      }
      glEnableClientState(GL_COLOR_ARRAY);
    }
#else
    if (SettingGet(I->G, cSetting_cgo_shader_ub_color)){
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);
    } else {
      glColorPointer(4, GL_FLOAT, 0, 0);
    }
    glEnableClientState(GL_COLOR_ARRAY);
#endif
#endif
  }

  if (bufs[3]){
    glBindBuffer(GL_ARRAY_BUFFER, bufs[3]);
#ifdef PURE_OPENGL_ES_2
    glEnableVertexAttribArray(attr_a_Accessibility);
    glVertexAttribPointer(attr_a_Accessibility, 1, GL_FLOAT, GL_FALSE, 0, 0);
#else
#ifdef OPENGL_ES_2
    if (I->use_shader){
      glEnableVertexAttribArray(attr_a_Accessibility);
      glVertexAttribPointer(attr_a_Accessibility, 1, GL_FLOAT, GL_FALSE, 0, 0);
    } else {
      glVertexPointer(1, GL_FLOAT, 0, 0);
      glEnableClientState(GL_VERTEX_ARRAY);
    }
#else
    glVertexPointer(1, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);
#endif
#endif
  } else {
    glVertexAttrib1f(attr_a_Accessibility, 1.f);
  }

  if (I->debug){
    mode = CGOConvertDebugMode(I->debug, mode);
  }
  glDrawArrays(mode, 0, nverts);

#ifdef PURE_OPENGL_ES_2
  if (bufs[0]) glDisableVertexAttribArray(attr_a_Vertex);
  if (bufs[1]) glDisableVertexAttribArray(attr_a_Normal);
  if (I->isPicking){
    glDisableVertexAttribArray(attr_a_Color);
  } else if (bufs[2]){
    glDisableVertexAttribArray(attr_a_Color);
  }
#else
#ifdef OPENGL_ES_2
  if (I->use_shader){
    if (bufs[0]) glDisableVertexAttribArray(attr_a_Vertex);
    if (bufs[1]) glDisableVertexAttribArray(attr_a_Normal);
    if (I->isPicking){
      glDisableVertexAttribArray(attr_a_Color);
    } else if (bufs[2]){
      glDisableVertexAttribArray(attr_a_Color);
    }
  } else {
    if (bufs[0]) glDisableClientState(GL_VERTEX_ARRAY);
    if (bufs[1]) glDisableClientState(GL_NORMAL_ARRAY);
    if (I->isPicking){
      glDisableClientState(GL_COLOR_ARRAY);
    } else if (bufs[2]){
      glDisableClientState(GL_COLOR_ARRAY);
    }
  }
#else
  if (bufs[0]) glDisableVertexAttribArray(attr_a_Vertex);
  if (bufs[1]) glDisableVertexAttribArray(attr_a_Normal);
  if (I->isPicking){
    glDisableVertexAttribArray(attr_a_Color);
  } else if (bufs[2]){
    glDisableVertexAttribArray(attr_a_Color);
  }
#endif
#endif
  if (bufs[3]) glDisableVertexAttribArray(attr_a_Accessibility);
  *pc += nverts*3 + 8;
  if (I->enable_shaders){
    CShaderPrg_Disable(shaderPrg);
  }
#endif
}

static void CGO_gl_draw_sphere_buffers(CCGORenderer * I, float **pc) {
#ifndef _PYMOL_PURE_OPENGL_ES
#ifdef _PYMOL_CGO_DRAWBUFFERS
  int  num_spheres = CGO_get_int(*pc);
  int ub_flags = CGO_get_int(*pc+1);
  int attr_a_vertex_radius;
  int attr_color;
  int attr_rightup;
  uint bufs[3] = { CGO_get_int(*pc+2), CGO_get_int(*pc+3),
                   CGO_get_int(*pc+4) };
  CShaderPrg *shaderPrg;

  if (I->enable_shaders){
    shaderPrg = CShaderPrg_Enable_SphereShader(I->G, "sphere");
  } else {
    shaderPrg = CShaderMgr_GetShaderPrg(I->G->ShaderMgr, "sphere");
  }
  attr_a_vertex_radius = CShaderPrg_GetAttribLocation(shaderPrg, "a_vertex_radius");
  attr_color = CShaderPrg_GetAttribLocation(shaderPrg, "a_Color"); 
  attr_rightup = CShaderPrg_GetAttribLocation(shaderPrg, "a_rightUpFlags");

  glEnableVertexAttribArray(attr_a_vertex_radius);
  glBindBuffer(GL_ARRAY_BUFFER, bufs[0]); // vertex xyz + right/up flag in w
  glVertexAttribPointer(attr_a_vertex_radius, 4, GL_FLOAT, GL_FALSE, 0, 0);

  glEnableVertexAttribArray(attr_color);
  glBindBuffer(GL_ARRAY_BUFFER, bufs[1]); // color in UNSIGNED_BYTE
  if (ub_flags & 1){
    glVertexAttribPointer(attr_color, 4, GL_UNSIGNED_BYTE, GL_TRUE, 0, 0);
  } else { 
    glVertexAttribPointer(attr_color, 4, GL_FLOAT, GL_FALSE, 0, 0);
  }
  glEnableVertexAttribArray(attr_rightup);
  glBindBuffer(GL_ARRAY_BUFFER, bufs[2]); // radius
  if (ub_flags & 2){
    glVertexAttribPointer(attr_rightup, 1, GL_UNSIGNED_BYTE, GL_FALSE, 0, 0);
  } else {
    glVertexAttribPointer(attr_rightup, 1, GL_FLOAT, GL_FALSE, 0, 0);
  }
  glDrawArrays(GL_QUADS, 0, num_spheres * 4);

  glDisableVertexAttribArray(attr_a_vertex_radius);
  glDisableVertexAttribArray(attr_color);
  glDisableVertexAttribArray(attr_rightup);
  if (I->enable_shaders){
    CShaderPrg_Disable(shaderPrg);
  }
#endif
#endif
}

static void CGO_gl_draw_cylinder_buffers(CCGORenderer * I, float **pc) {
#ifdef _PYMOL_CGO_DRAWBUFFERS
  int  num_cyl = CGO_get_int(*pc);
  int min_alpha = CGO_get_int(*pc+1);
  int attr_origin;
  int attr_axis;
  int attr_colors, attr_colors2;
  uint bufs[5] = { CGO_get_int(*pc+2), CGO_get_int(*pc+3),
                   CGO_get_int(*pc+4), CGO_get_int(*pc+5), CGO_get_int(*pc+6) };
  CShaderPrg *shaderPrg;

  if (I->enable_shaders){
    shaderPrg = CShaderPrg_Enable_CylinderShader(I->G);
  } else {
    shaderPrg = CShaderMgr_GetShaderPrg(I->G->ShaderMgr, "cylinder");
  }
  attr_origin = CShaderPrg_GetAttribLocation(shaderPrg, "attr_origin");
  attr_axis = CShaderPrg_GetAttribLocation(shaderPrg, "attr_axis"); 
  attr_colors = CShaderPrg_GetAttribLocation(shaderPrg, "attr_colors");
  attr_colors2 = CShaderPrg_GetAttribLocation(shaderPrg, "attr_colors2");

  glEnableVertexAttribArray(attr_origin);
  glBindBuffer(GL_ARRAY_BUFFER, bufs[0]); // cylinder origin + radius
  glVertexAttribPointer(attr_origin, 4, GL_FLOAT, GL_FALSE, 0, 0);
 
  glEnableVertexAttribArray(attr_axis);
  glBindBuffer(GL_ARRAY_BUFFER, bufs[1]); // cylinder axis + flags
  glVertexAttribPointer(attr_axis, 4, GL_FLOAT, GL_FALSE, 0, 0);

  if (bufs[2]){
    glEnableVertexAttribArray(attr_colors);
    glBindBuffer(GL_ARRAY_BUFFER, bufs[2]); // colors
    glVertexAttribPointer(attr_colors, 4, GL_FLOAT, GL_FALSE, 0, 0);
  }

  if (bufs[2]||bufs[3]){
    glEnableVertexAttribArray(attr_colors2);
    if (bufs[3]){
      glBindBuffer(GL_ARRAY_BUFFER, bufs[3]); // colors2
    } else if (bufs[2]) {
      glBindBuffer(GL_ARRAY_BUFFER, bufs[2]); // colors
    }
    glVertexAttribPointer(attr_colors2, 4, GL_FLOAT, GL_FALSE, 0, 0);
  }
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufs[4]);

  if (min_alpha < 255) {
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glDrawElements(GL_TRIANGLES, num_cyl * 36, GL_UNSIGNED_INT, 0);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthFunc(GL_LEQUAL);
  }

  glDrawElements(GL_TRIANGLES, num_cyl * 36, GL_UNSIGNED_INT, 0);

  if (min_alpha < 255) {
    glDepthFunc(GL_LESS);
  }

  glDisableVertexAttribArray(attr_origin);
  glDisableVertexAttribArray(attr_axis);
  if (bufs[2]||bufs[3]){
    glDisableVertexAttribArray(attr_colors);
    glDisableVertexAttribArray(attr_colors2);
  }
  if (I->enable_shaders){
    CShaderPrg_Disable(shaderPrg);
  }
#endif
}

static void CGO_gl_linewidth(CCGORenderer * I, float **pc)
{
  glLineWidth(**pc);
}

static void CGO_gl_linewidth_special(CCGORenderer * I, float **pc)
{
  int mode = CGO_get_int(*pc);
  switch (mode){
  case LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON:
    {
      float line_width = SceneGetDynamicLineWidth(I->info, SettingGet_f(I->G, NULL, NULL, cSetting_ribbon_width));
      if (I->info->width_scale_flag){
	glLineWidth(line_width * I->info->width_scale);
      } else {
	glLineWidth(line_width);
      }
    }
    break;
  case LINEWIDTH_DYNAMIC_WITH_SCALE_DASH:
    {
      float line_width = SceneGetDynamicLineWidth(I->info, SettingGet_f(I->G, NULL, NULL, cSetting_dash_width));
      if (I->info->width_scale_flag){
	glLineWidth(line_width * I->info->width_scale);
      } else {
	glLineWidth(line_width);
      }
    }
    break;
  case LINEWIDTH_DYNAMIC_WITH_SCALE:
    {
      float line_width = SceneGetDynamicLineWidth(I->info, SettingGet_f(I->G, NULL, NULL, cSetting_line_width));
      if (I->info->width_scale_flag){
	glLineWidth(line_width * I->info->width_scale);
      } else {
	glLineWidth(line_width);
      }
    }
    break;
  case LINEWIDTH_DYNAMIC_MESH:
    {
      float line_width;
      if (I->rep){
	line_width = SettingGet_f(I->G, I->rep->cs->Setting, I->rep->obj->Setting, cSetting_mesh_width);
      } else {
	line_width = SettingGet_f(I->G, NULL, NULL, cSetting_mesh_width);
      }
      line_width = SceneGetDynamicLineWidth(I->info, line_width);
      glLineWidth(line_width);
    }
    break;
  case POINTSIZE_DYNAMIC_DOT_WIDTH:
    {
      CSetting *csSetting = NULL, *objSetting = NULL;
      if (I->rep && I->rep->cs){
	csSetting = I->rep->cs->Setting;
      }
      if (I->rep && I->rep->obj){
	objSetting = I->rep->obj->Setting;
      }
      if(I->info->width_scale_flag){
	glPointSize(SettingGet_f
		    (I->G, csSetting, objSetting,
		     cSetting_dot_width) * I->info->width_scale);
      }
      else {
	glPointSize(SettingGet_f
		    (I->G, csSetting, objSetting, cSetting_dot_width));
      }
      break;
    }
  default:
    PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_linewidth_special(): bad mode=%d\n", mode ENDFB(I->G);
  }
}

static void CGO_gl_dotwidth(CCGORenderer * I, float **pc)
{
#ifdef PURE_OPENGL_ES_2      
      /* TODO */
#else
  glPointSize(**pc);
#endif
}

static void CGO_gl_enable(CCGORenderer * I, float **pc)
{
  GLenum mode = CGO_get_int(*pc);
#ifdef PURE_OPENGL_ES_2
#else
#ifdef OPENGL_ES_2
  if (I->use_shader){
    if (mode==GL_SHADER_LIGHTING && !I->isPicking){
      CShaderPrg * shaderPrg = CShaderMgr_GetShaderPrg(I->G->ShaderMgr, "default");
      CShaderPrg_Set1i(shaderPrg, "lighting_enabled", 1);
    }
  } else {
    if (mode!=GL_LIGHTING || !I->isPicking){
      glEnable(mode);
    }
  }
#else
  if (mode!=GL_LIGHTING || !I->isPicking){
    glEnable(mode);
  }
#endif
#endif
}

static void CGO_gl_disable(CCGORenderer * I, float **pc)
{
  GLenum mode = CGO_get_int(*pc);

#ifdef PURE_OPENGL_ES_2
#else
#ifdef OPENGL_ES_2
  if (I->use_shader){
    if (mode==GL_SHADER_LIGHTING && !I->isPicking){
      CShaderPrg * shaderPrg = CShaderMgr_GetShaderPrg(I->G->ShaderMgr, "default");
      CShaderPrg_Set1i(shaderPrg, "lighting_enabled", 0);
    }
  } else {
    if (mode!=GL_LIGHTING || !I->isPicking){
      glDisable(mode);
    }
  }
#else
  if (mode!=GL_LIGHTING || !I->isPicking){
    glDisable(mode);
  }
#endif
#endif
}

static void CGO_gl_alpha(CCGORenderer * I, float **pc)
{
  I->alpha = **pc;
}

static void CGO_gl_reset_normal(CCGORenderer * I, float **pc)
{
  SceneResetNormalUseShader(I->G, CGO_get_int(*pc), I->use_shader);
}

static void CGO_gl_null(CCGORenderer * I, float **pc)
{
}

static void CGO_gl_error(CCGORenderer * I, float **pc)
{
  PRINTFB(I->G, FB_CGO, FB_Warnings) " CGO_gl_error() is not suppose to be called op=%d with mask=%d\n", CGO_get_int((*pc)-1), CGO_MASK & CGO_get_int((*pc)-1) ENDFB(I->G);
}

static void CGO_gl_color(CCGORenderer * I, float **varg)
{
  float *v = *varg;
#ifdef PURE_OPENGL_ES_2
  {
    if (I->G->ShaderMgr->current_shader){
      int attr_a_Color = CShaderPrg_GetAttribLocation(I->G->ShaderMgr->current_shader, "a_Color");
      glVertexAttrib4f(attr_a_Color, v[0], v[1], v[2], I->alpha);
    }
  }
#else
#ifdef OPENGL_ES_2
  if (I->use_shader){
    if (I->G->ShaderMgr->current_shader){
      int attr_a_Color = CShaderPrg_GetAttribLocation(I->G->ShaderMgr->current_shader, "a_Color");
      glVertexAttrib4f(attr_a_Color, v[0], v[1], v[2], I->alpha);
    }
  } else {
    glColor4f(v[0], v[1], v[2], I->alpha);
  }
#else
  glColor4f(v[0], v[1], v[2], I->alpha);
#endif
#endif
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
  CGO_gl_null,                  /* 0x07 */
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

#ifdef _PYMOL_CGO_DRAWARRAYS
  CGO_gl_draw_arrays,           /* 0x1C DrawArrays() */
#else
  CGO_gl_null,                  /* 0X1C */
#endif
  CGO_gl_null,                  /* 0x1D */
  CGO_gl_reset_normal,          /* 0x1E */
  CGO_gl_null,                  /* pick color  0X1F */
  CGO_gl_draw_buffers,          /* 0x20 draw buffers */
  CGO_gl_draw_buffers_indexed,          /* 0x21 draw buffers indexed */
  CGO_gl_null,                  /* 0x22 bounding box */
  CGO_gl_draw_buffers_not_indexed,          /* 0x23 draw buffers not indexed */
  CGO_gl_linewidth_special,  /* 0x24 linewidth special */
  CGO_gl_draw_cylinder_buffers,  /* 0x25 draw GLSL cylinders */
  CGO_gl_null,                  /* 0x26 shader cylinder */
  CGO_gl_null,                  /* 0x27 shader cylinder with 2nd color */
  CGO_gl_draw_sphere_buffers,   /* 0x28 draw sphere buffers */
  CGO_gl_null,                  /* 0x29 accessibility used for ambient occlusion */
  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,
  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,
  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error,  CGO_gl_error
};

void CGORenderGLPicking(CGO * I, Picking ** pick, PickContext * context, CSetting * set1,
                        CSetting * set2)
{
  register PyMOLGlobals *G = I->G;
  if(G->ValidContext) {
    float *pc = I->op;
    register int op;
    register CCGORenderer *R = G->CGORenderer;
    int i, j;
    Picking *p;
    R->use_shader = I->use_shader;
    R->isPicking = true;
    if(I->c) {
      i = (*pick)[0].src.index;

      glLineWidth(SettingGet_f(G, set1, set2, cSetting_cgo_line_width));
      while((op = (CGO_MASK & CGO_read_int(pc)))) {
        if(op != CGO_PICK_COLOR) {
#ifdef _PYMOL_CGO_DRAWARRAYS
	  if (op == CGO_DRAW_ARRAYS){
	    int arrays = CGO_get_int(pc+1);
	    if (arrays & CGO_PICK_COLOR_ARRAY){
	      int nverts = CGO_get_int(pc+3), v, pl, idx = -1, pidx, bnd = -1, pbnd;
	      float *pca = pc + 4;
	      float *pickColorVals ;
	      uchar *pickColorValsUC;
	      if (arrays & CGO_VERTEX_ARRAY){ pca += nverts * 3; }
	      if (arrays & CGO_NORMAL_ARRAY){ pca += nverts * 3; }
	      if (arrays & CGO_COLOR_ARRAY){ pca += nverts * 4; }
	      pickColorValsUC = (uchar*)pca;
	      pickColorVals = pca + nverts;
	      pl =0;
	      for (v=0;v<nverts; v++){
		pidx = idx;
		pbnd = bnd;
		idx = CGO_get_int(pickColorVals + (v * 2));
		bnd = CGO_get_int(pickColorVals + (v * 2) + 1);
		i++;
		if(!(*pick)[0].src.bond) {
		  pickColorValsUC[pl++] = ((i & 0xF) << 4);
		  pickColorValsUC[pl++] = ((i & 0xF0) | 0x8);
		  pickColorValsUC[pl++] = ((i & 0xF00) >> 4);
		  pickColorValsUC[pl++] = 255;
		  VLACheck((*pick), Picking, i);
		  p = (*pick) + i;
		  p->context = (*context);
		  p->src.index = (int)idx;
		  p->src.bond = (int) bnd;      /* actually holds state information */
		  I->current_pick_color_index = idx;
		  I->current_pick_color_bond = bnd;
		} else {
		  j = i >> 12;
		  pickColorValsUC[pl++] = ((j & 0xF) << 4);
		  pickColorValsUC[pl++] = ((j & 0xF0) | 0x8);
		  pickColorValsUC[pl++] = ((j & 0xF00) >> 4);
		  pickColorValsUC[pl++] = 255;
		}
	      }
	    }
            CGO_gl[op] (R, &pc);
	  }
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
	  else if (op == CGO_DRAW_BUFFERS_INDEXED || op == CGO_DRAW_BUFFERS_NOT_INDEXED){
	    int nverts = 0, v, pl, idx = -1, pidx, bnd = -1, pbnd, chg = 0;
	    float *pca = NULL;
	    float *pickColorVals ;
	    uchar *pickColorValsUC;
	    switch (op){
	    case CGO_DRAW_BUFFERS_INDEXED:
	      nverts = CGO_get_int(pc+4); 
	      pca = pc + 10;
	      break;
	    case CGO_DRAW_BUFFERS_NOT_INDEXED:
	      nverts = CGO_get_int(pc+3);
	      pca = pc + 8;
	      break;
	    }
	    pickColorValsUC = (uchar*)pca;
	    pickColorVals = pca + nverts;
	    pl =0;
	    for (v=0;v<nverts; v++){
	      pidx = idx;
	      pbnd = bnd;
	      idx = CGO_get_int(pickColorVals + v * 2);
	      bnd = CGO_get_int(pickColorVals + v * 2 + 1);
	      chg = idx != pidx || bnd != pbnd;
	      if (chg)
		i++;
	      if(!(*pick)[0].src.bond) {
		pickColorValsUC[pl++] = ((i & 0xF) << 4);
		pickColorValsUC[pl++] = ((i & 0xF0) | 0x8);
		pickColorValsUC[pl++] = ((i & 0xF00) >> 4);
		pickColorValsUC[pl++] = 255;
		if (chg){
		  VLACheck((*pick), Picking, i);
		  p = (*pick) + i;
		  p->context = (*context);
		  p->src.index = (int)idx;
		  p->src.bond = (int) bnd;      /* actually holds state information */
		  I->current_pick_color_index = idx;
		  I->current_pick_color_bond = bnd;
		}
	      } else {
		j = i >> 12;
		pickColorValsUC[pl++] = ((j & 0xF) << 4);
		pickColorValsUC[pl++] = ((j & 0xF0) | 0x8);
		pickColorValsUC[pl++] = ((j & 0xF00) >> 4);
		pickColorValsUC[pl++] = 255;
	      }
	    }
	  } else 
#endif
          if(op != CGO_COLOR) {
            CGO_gl[op] (R, &pc); /* ignore color changes */
	  }
	  pc += CGO_sz[op];
        } else {
          i++;
          if(!(*pick)[0].src.bond) {
            /* pass 1 - low order bits */
#ifdef PURE_OPENGL_ES_2      
      /* TODO */
#else
#ifdef OPENGL_ES_2
	    if (I->use_shader){
	      GLubyte col[] = { (GLubyte) ((i & 0xF) << 4), (GLubyte) ((i & 0xF0) | 0x8), (GLubyte) ((i & 0xF00) >> 4), 255 };
	      if (I->G->ShaderMgr->current_shader){
		int attr_a_Color = CShaderPrg_GetAttribLocation(I->G->ShaderMgr->current_shader, "a_Accessibility");
		glVertexAttrib4ubv(attr_a_Color, col);       /* we're encoding the index into the color */
	      }
	    } else {
	      glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8), (uchar) ((i & 0xF00) >> 4));       /* we're encoding the index into the color */
	    }
#else
            glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8), (uchar) ((i & 0xF00) >> 4));       /* we're encoding the index into the color */
#endif
#endif
            VLACheck((*pick), Picking, i);
            p = (*pick) + i;
            p->context = (*context);
            p->src.index = CGO_get_int(pc);
            p->src.bond = CGO_get_int(pc + 1);      /* actually holds state information */
            I->current_pick_color_index = p->src.index;
            I->current_pick_color_bond = p->src.bond;
          } else {
            /* pass 2 - high order bits */

            j = i >> 12;

 #ifdef PURE_OPENGL_ES_2      
      /* TODO */
#else
#ifdef OPENGL_ES_2
	    {
	      int attr_a_Color = CShaderPrg_GetAttribLocation(I->G->ShaderMgr->current_shader, "a_Accessibility");
	      if (I->use_shader){
		GLubyte col[] = { (GLubyte) ((j & 0xF) << 4), (GLubyte) ((j & 0xF0) | 0x8), (GLubyte) ((j & 0xF00) >> 4), 255 };
		glVertexAttrib4ubv(attr_a_Color, col);       /* we're encoding the index into the color */
	      } else {
		glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8), (uchar) ((j & 0xF00) >> 4));       /* we're encoding the index into the color */
	      }
	    }
#else
           glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                       (uchar) ((j & 0xF00) >> 4));
#endif
#endif
          }
	  pc += CGO_sz[op];
        }
      }
      (*pick)[0].src.index = i; /* pass the count */
    }
    R->isPicking = false;
  }
}

/* This DEBUG_PRINT_OPS preprocessor, if defined, will print all OPS of every CGO rendered
   by CGORenderGL().  This is only to be used in debugging */
//#define DEBUG_PRINT_OPS

void CGORenderGL(CGO * I, float *color, CSetting * set1, CSetting * set2,
                 RenderInfo * info, Rep *rep)

/* this should be as fast as you can make it...

 * the ASM loop is about 2X long as raw looped GL calls,

 * but hopefully superscaler processors won't care */
{
  register PyMOLGlobals *G = I->G;
  if(G->ValidContext) {
    float *pc = I->op;
    register int op;
    register CCGORenderer *R = G->CGORenderer;
    register float _1 = 1.0F;
#ifdef DEBUG_PRINT_OPS
    CGOCountNumberOfOperationsOfType(I, 0);
#endif
    R->info = info;
    R->use_shader = I->use_shader;
    R->enable_shaders = I->enable_shaders;
    R->debug = I->debug;
    R->rep = rep;
    SceneResetNormalUseShader(I->G, true, I->use_shader);
    if(I->c) {
      R->alpha = 1.0F - SettingGet_f(I->G, set1, set2, cSetting_cgo_transparency);
#ifdef PURE_OPENGL_ES_2
      if (I->G->ShaderMgr->current_shader){
	int attr_a_Color = CShaderPrg_GetAttribLocation(I->G->ShaderMgr->current_shader, "a_Color");
	if (attr_a_Color>=0){
	  if(color)
	    glVertexAttrib4f(attr_a_Color, color[0], color[1], color[2], R->alpha);
	  else
	    glVertexAttrib4f(attr_a_Color, 1.0, 1.0, 1.0, R->alpha);
	}
      }
#else
#ifdef OPENGL_ES_2
      if (I->use_shader){
	if (I->G->ShaderMgr->current_shader){
	  int attr_a_Color = CShaderPrg_GetAttribLocation(I->G->ShaderMgr->current_shader, "a_Color");
	  if (attr_a_Color>=0){
	    if(color)
	      glVertexAttrib4f(attr_a_Color, color[0], color[1], color[2], R->alpha);
	    else
	      glVertexAttrib4f(attr_a_Color, 1.0, 1.0, 1.0, R->alpha);
	  }
	}
      } else {
	if(color)
	  glColor4f(color[0], color[1], color[2], R->alpha);
	else
	  glColor4f(1.0, 1.0, 1.0, R->alpha);
      }
#else
      if(color)
        glColor4f(color[0], color[1], color[2], R->alpha);
      else
        glColor4f(1.0, 1.0, 1.0, R->alpha);
#endif
#endif
#ifdef PURE_OPENGL_ES_2
#else
      if(info && info->width_scale_flag) {
        glLineWidth(SettingGet_f(I->G, set1, set2, cSetting_cgo_line_width) *
                    info->width_scale);
        glPointSize(SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_width) *
                    info->width_scale);

      } else {
        glLineWidth(SettingGet_f(I->G, set1, set2, cSetting_cgo_line_width));
        glPointSize(SettingGet_f(I->G, set1, set2, cSetting_cgo_dot_width));
      }
#endif
      if(info->alpha_cgo) {     /* we're sorting transparent triangles globally */
        register int mode = -1;
        float *n0 = NULL, *n1 = NULL, *n2 = NULL, *v0 = NULL, *v1 = NULL, *v2 =
          NULL, *c0 = NULL, *c1 = NULL, *c2 = NULL;
        float zee[] = { 0.0, 0.0, 1.0 }, color_tmp[] = { 1., 1., 1. };
        int vc = 0;
	if (color){
	  I->color[0] = color[0]; I->color[1] = color[1]; I->color[2] = color[2];
	  c0 = I->color;
	} else {
	  c0 = color_tmp;
	}
        while((op = (CGO_MASK & CGO_read_int(pc)))) {
          if((R->alpha != _1)) {
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
#ifdef _PYMOL_CGO_DRAWARRAYS
	    case CGO_DRAW_ARRAYS:
	      {
		int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), nverts = CGO_get_int(pc + 3);
		float *vertexVals = 0, *nxtVals = 0, *colorVals = 0, *normalVals;
		float *vertexVals_tmp = 0, *colorVals_tmp = 0, *normalVals_tmp = 0;
		int step;	      
		short nxtn = 3;
		nxtVals = vertexVals = vertexVals_tmp = pc + 4;
		pc += nverts*3;
		if (arrays & CGO_NORMAL_ARRAY){
		  nxtVals = normalVals = normalVals_tmp = vertexVals + (nxtn*nverts);
		  pc += nverts*3;
		}
		if (arrays & CGO_COLOR_ARRAY){
		  nxtVals = colorVals = colorVals_tmp = nxtVals + (nxtn*nverts);
		  nxtn = 4;
		  pc += nverts*4;
		}
		if (arrays & CGO_PICK_COLOR_ARRAY){
		  pc += nverts*3;
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
#endif
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
          pc += CGO_sz[op];
        }
      } else {
	int nops = 0;
        while((op = (CGO_MASK & CGO_read_int(pc)))) {
          CGO_gl[op] (R, &pc);
          pc += CGO_sz[op];
	  nops++;
        }
	//	printf("nops=%d\n", nops);
      }
    }
  }
}

void CGORenderGLAlpha(CGO * I, RenderInfo * info)
{
  register PyMOLGlobals *G = I->G;
  if(G->ValidContext && I->c) {

    /* 1. transform and measure range (if not already known) 
       2. bin into linked lists based on Z-centers
       3. render by layer */

    if(I->z_flag) {
      if(!I->i_start) {
        I->i_size = 256;
        I->i_start = Calloc(int, I->i_size);
      } else {
        UtilZeroMem(I->i_start, sizeof(int) * I->i_size);
      }
      {
        register float z_min = I->z_min;
        register int i_size = I->i_size;
        register float range_factor = (0.9999F * i_size) / (I->z_max - z_min);
        register float *base = I->op;
        register float *pc = base;
        register int op, i, ii;
        register int *start = I->i_start, *init_start;
        register int delta = 1;
        /* bin the triangles */
        while((op = (CGO_MASK & CGO_read_int(pc)))) {
          switch (op) {
          case CGO_ALPHA_TRIANGLE:
            i = (int) ((pc[4] - z_min) * range_factor);
            if(i < 0)
              i = 0;
            if(i >= i_size)
              i = i_size;
            CGO_put_int(pc, start[i]);
            start[i] = (pc - base);     /* NOTE: will always be > 0 since we have CGO_read_int'd */
            break;
          }
          pc += CGO_sz[op];
        }
        if(SettingGetGlobal_i(G, cSetting_transparency_mode) == 2) {
          delta = -1;
          start += (i_size - 1);
        }

        /* now render by bin */
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
#ifdef _PYMOL_GL_DRAWARRAYS
	{
	  int num_vertices = 0;
	  init_start = start;
	  for(i = 0; i < i_size; i++) {
	    ii = *start;
	    start += delta;
	    while(ii) {
	      pc = base + ii;
	      num_vertices += 3;
	      ii = CGO_get_int(pc);
	    }
	  }
	  {
	    ALLOCATE_ARRAY(GLfloat,colorVals, num_vertices*4)
	    ALLOCATE_ARRAY(GLfloat,normalVals, num_vertices*3)
	    ALLOCATE_ARRAY(GLfloat,vertexVals, num_vertices*3)
	    int array_place_4 = 0, array_place_3 = 0;
	    float *tmp_ptr;
	    start = init_start;
	    for(i = 0; i < i_size; i++) {
	      ii = *start;
	      start += delta;
	      while(ii) {
		pc = base + ii;
		tmp_ptr = pc + 23;
		colorVals[0] = tmp_ptr[0]; colorVals[1] = tmp_ptr[1]; colorVals[2] = tmp_ptr[2]; colorVals[3] = tmp_ptr[3];
		tmp_ptr = pc + 14;
		normalVals[0] = tmp_ptr[0]; normalVals[1] = tmp_ptr[1]; normalVals[2] = tmp_ptr[2];
		tmp_ptr = pc + 5;
		vertexVals[0] = tmp_ptr[0]; vertexVals[1] = tmp_ptr[1]; vertexVals[2] = tmp_ptr[2];
		array_place_4 += 4;
		array_place_3 += 3;
		tmp_ptr = pc + 27;
		colorVals[0] = tmp_ptr[0]; colorVals[1] = tmp_ptr[1]; colorVals[2] = tmp_ptr[2]; colorVals[3] = tmp_ptr[3];
		tmp_ptr = pc + 17;
		normalVals[0] = tmp_ptr[0]; normalVals[1] = tmp_ptr[1]; normalVals[2] = tmp_ptr[2];
		tmp_ptr = pc + 8;
		vertexVals[0] = tmp_ptr[0]; vertexVals[1] = tmp_ptr[1]; vertexVals[2] = tmp_ptr[2];
		array_place_4 += 4;
		array_place_3 += 3;
		tmp_ptr = pc + 31;
		colorVals[0] = tmp_ptr[0]; colorVals[1] = tmp_ptr[1]; colorVals[2] = tmp_ptr[2]; colorVals[3] = tmp_ptr[3];
		tmp_ptr = pc + 20;
		normalVals[0] = tmp_ptr[0]; normalVals[1] = tmp_ptr[1]; normalVals[2] = tmp_ptr[2];
		tmp_ptr = pc + 11;
		vertexVals[0] = tmp_ptr[0]; vertexVals[1] = tmp_ptr[1]; vertexVals[2] = tmp_ptr[2];
		array_place_4 += 4;
		array_place_3 += 3;
		ii = CGO_get_int(pc);
	      }
	    }
	    glEnableClientState(GL_VERTEX_ARRAY);
	    glEnableClientState(GL_COLOR_ARRAY);
	    glEnableClientState(GL_NORMAL_ARRAY);
	    glVertexPointer(3, GL_FLOAT, 0, vertexVals);
	    glColorPointer(4, GL_FLOAT, 0, colorVals);
	    glNormalPointer(GL_FLOAT, 0, normalVals);
	    glDrawArrays(GL_TRIANGLES, 0, num_vertices);
	    glDisableClientState(GL_VERTEX_ARRAY);
	    glDisableClientState(GL_COLOR_ARRAY);
	    glDisableClientState(GL_NORMAL_ARRAY);
	    DEALLOCATE_ARRAY(colorVals)
	    DEALLOCATE_ARRAY(normalVals)
	    DEALLOCATE_ARRAY(vertexVals)
	  }
	}
#else
	(void)init_start;
        glBegin(GL_TRIANGLES);
        for(i = 0; i < i_size; i++) {
          ii = *start;
          start += delta;
          while(ii) {
            pc = base + ii;

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
#endif
      }
    } else {
      register float *pc = I->op;
      register int op;
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
#ifdef _PYMOL_GL_DRAWARRAYS
      {
	int num_vertices = 0;
	while((op = (CGO_MASK & CGO_read_int(pc)))) {
	  switch (op) {
	  case CGO_ALPHA_TRIANGLE:
	    num_vertices += 3;
	    break;
	  }
	  pc += CGO_sz[op];
	}
	{
	  int array_place_4 = 0, array_place_3 = 0;
	  float *tmp_ptr;
	  ALLOCATE_ARRAY(GLfloat,colorVals, num_vertices*4)
	  ALLOCATE_ARRAY(GLfloat,normalVals, num_vertices*3)
	  ALLOCATE_ARRAY(GLfloat,vertexVals, num_vertices*3)
	  pc = I->op;
	  while((op = (CGO_MASK & CGO_read_int(pc)))) {
	    switch (op) {
	    case CGO_ALPHA_TRIANGLE:
	      tmp_ptr = pc + 23;
	      colorVals[array_place_4++] = tmp_ptr[0]; colorVals[array_place_4++] = tmp_ptr[1];
	      colorVals[array_place_4++] = tmp_ptr[2]; colorVals[array_place_4++] = tmp_ptr[3];
	      tmp_ptr = pc + 14;
	      normalVals[array_place_3] = tmp_ptr[0]; normalVals[array_place_3+1] = tmp_ptr[1]; normalVals[array_place_3+2] = tmp_ptr[2];
	      tmp_ptr = pc + 5;
	      vertexVals[array_place_3++] = tmp_ptr[0]; vertexVals[array_place_3++] = tmp_ptr[1]; vertexVals[array_place_3++] = tmp_ptr[2];
	      tmp_ptr = pc + 27;
	      colorVals[array_place_4++] = tmp_ptr[0]; colorVals[array_place_4++] = tmp_ptr[1];
	      colorVals[array_place_4++] = tmp_ptr[2]; colorVals[array_place_4++] = tmp_ptr[3];
	      tmp_ptr = pc + 17;
	      normalVals[array_place_3] = tmp_ptr[0]; normalVals[array_place_3+1] = tmp_ptr[1]; normalVals[array_place_3+2] = tmp_ptr[2];
	      tmp_ptr = pc + 8;
	      vertexVals[array_place_3++] = tmp_ptr[0]; vertexVals[array_place_3++] = tmp_ptr[1]; vertexVals[array_place_3++] = tmp_ptr[2];
	      tmp_ptr = pc + 31;
	      colorVals[array_place_4++] = tmp_ptr[0]; colorVals[array_place_4++] = tmp_ptr[1];
	      colorVals[array_place_4++] = tmp_ptr[2]; colorVals[array_place_4++] = tmp_ptr[3];
	      tmp_ptr = pc + 20;
	      normalVals[array_place_3] = tmp_ptr[0]; normalVals[array_place_3+1] = tmp_ptr[1]; normalVals[array_place_3+2] = tmp_ptr[2];
	      tmp_ptr = pc + 11;
	      vertexVals[array_place_3++] = tmp_ptr[0]; vertexVals[array_place_3++] = tmp_ptr[1]; vertexVals[array_place_3++] = tmp_ptr[2];
	      break;
	    }
	    pc += CGO_sz[op];
	  }
	  glEnableClientState(GL_VERTEX_ARRAY);
	  glEnableClientState(GL_COLOR_ARRAY);
	  glEnableClientState(GL_NORMAL_ARRAY);
	  glVertexPointer(3, GL_FLOAT, 0, vertexVals);
	  glColorPointer(4, GL_FLOAT, 0, colorVals);
	  glNormalPointer(GL_FLOAT, 0, normalVals);
	  glDrawArrays(GL_TRIANGLES, 0, num_vertices);
	  glDisableClientState(GL_VERTEX_ARRAY);
	  glDisableClientState(GL_COLOR_ARRAY);
	  glDisableClientState(GL_NORMAL_ARRAY);
	  DEALLOCATE_ARRAY(colorVals)
	  DEALLOCATE_ARRAY(normalVals)
	  DEALLOCATE_ARRAY(vertexVals)
	}
      }
#else
      glBegin(GL_TRIANGLES);
      while((op = (CGO_MASK & CGO_read_int(pc)))) {
        switch (op) {
        case CGO_ALPHA_TRIANGLE:
          glColor4fv(pc + 23);
          glNormal3fv(pc + 14);
          glVertex3fv(pc + 5);
          glColor4fv(pc + 27);
          glNormal3fv(pc + 17);
          glVertex3fv(pc + 8);
          glColor4fv(pc + 31);
          glNormal3fv(pc + 20);
          glVertex3fv(pc + 11);
          break;
        }
        pc += CGO_sz[op];
      }
      glEnd();
#endif
#endif
    }
  }
}


/* translation function which turns cylinders and spheres into triangles */

static void CGOSimpleSphere(CGO * I, float *v, float vdw)
{
  SphereRec *sp;
  int *q, *s;
  int b, c;
  int ds;

  /* cgo_sphere_quality is between 0 and (NUMBER_OF_SPHERE_LEVELS-1) */
  ds = SettingGet_i(I->G, NULL, NULL, cSetting_cgo_sphere_quality);

  sp = I->G->Sphere->Sphere[ds];

  q = sp->Sequence;
  s = sp->StripLen;

  for(b = 0; b < sp->NStrip; b++) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = (*s), pl = 0;
      float *vertexVals, *normalVals, *tmp_ptr;
      vertexVals = CGODrawArrays(I, GL_TRIANGLE_STRIP, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY, nverts);      
      normalVals = vertexVals + (nverts*3);
      for(c = 0; c < (*s); c++) {
	tmp_ptr = sp->dot[*q];
	normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
	vertexVals[pl++] = v[0] + vdw * sp->dot[*q][0]; 
	vertexVals[pl++] = v[1] + vdw * sp->dot[*q][1];
	vertexVals[pl++] = v[2] + vdw * sp->dot[*q][2];
	q++;
      }
    }
#else
    CGOBegin(I, GL_TRIANGLE_STRIP);
    for(c = 0; c < (*s); c++) {
      CGONormalv(I, sp->dot[*q]);
      CGOVertex(I, v[0] + vdw * sp->dot[*q][0],
                v[1] + vdw * sp->dot[*q][1], v[2] + vdw * sp->dot[*q][2]);
      q++;
    }
    CGOEnd(I);
#endif
    s++;
  }
}

static void CGOSimpleQuadric(CGO * I, float *v, float r, float *q)
{
  float r_el, n0[3], n1[3], n2[3];
  if(CGOQuadricToEllipsoid(v, r, q, &r_el, n0, n1, n2))
    CGOSimpleEllipsoid(I, v, r_el, n0, n1, n2);
}

static void CGOSimpleEllipsoid(CGO * I, float *v, float vdw, float *n0, float *n1,
                               float *n2)
{
  SphereRec *sp;
  int *q, *s;
  int b, c;
  int ds;
  float nn0[3], nn1[3], nn2[3];
  float scale[3], scale_sq[3];

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
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = (*s), pl = 0;
      float *vertexVals, *normalVals;
      vertexVals = CGODrawArrays(I, GL_TRIANGLE_STRIP, CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY, nverts);      
      normalVals = vertexVals + (nverts*3);
      for(c = 0; c < (*s); c++) {
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
	
	normalVals[pl] = surfnormal[0]; normalVals[pl+1] = surfnormal[1]; normalVals[pl+2] = surfnormal[2];
	vertexVals[pl++] = vv[0]; vertexVals[pl++] = vv[1]; vertexVals[pl++] = vv[2];
	q++;
      }
    }
#else
    CGOBegin(I, GL_TRIANGLE_STRIP);
    for(c = 0; c < (*s); c++) {
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

      CGONormalv(I, surfnormal);
      CGOVertexv(I, vv);
      q++;
    }
    CGOEnd(I);
#endif
    s++;
  }
}

static void subdivide(int n, float *x, float *y)
{
  int a;
  if(n < 3) {
    n = 3;
  }
  for(a = 0; a <= n; a++) {
    x[a] = (float) cos(a * 2 * PI / n);
    y[a] = (float) sin(a * 2 * PI / n);
  }
}

static void CGOSimpleCylinder(CGO * I, float *v1, float *v2, float tube_size, float *c1,
                              float *c2, int cap1, int cap2)
{
#define MAX_EDGE 50

  float d[3], t[3], p0[3], p1[3], p2[3], vv1[3], vv2[3], v_buf[9], *v;
  float x[MAX_EDGE + 1], y[MAX_EDGE + 1];
  float overlap;
  float nub;
  int colorFlag;
  int nEdge;
  int c;

  v = v_buf;
  nEdge = (int) SettingGet(I->G, cSetting_stick_quality);
  overlap = tube_size * SettingGet(I->G, cSetting_stick_overlap);
  nub = tube_size * SettingGet(I->G, cSetting_stick_nub);

  if(nEdge > MAX_EDGE)
    nEdge = MAX_EDGE;
  subdivide(nEdge, x, y);

  colorFlag = (c1 != c2) && c2;

  if (c1)
    CGOColorv(I, c1);

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  normalize3f(p0);

  if(cap1 == cCylCapRound) {
    vv1[0] = v1[0] - p0[0] * overlap;
    vv1[1] = v1[1] - p0[1] * overlap;
    vv1[2] = v1[2] - p0[2] * overlap;
  } else {
    vv1[0] = v1[0];
    vv1[1] = v1[1];
    vv1[2] = v1[2];
  }
  if(cap2 == cCylCapRound) {
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

  get_divergent3f(d, t);

  cross_product3f(d, t, p1);

  normalize3f(p1);

  cross_product3f(d, p1, p2);

  normalize3f(p2);

  /* now we have a coordinate system */

#ifdef _PYMOL_CGO_DRAWARRAYS
  {
    int nverts = (2*nEdge) + 2, pl = 0, plc = 0, damode = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY;
    float *vertexVals;
    float *normalVals, *colorVals = 0, *tmp_ptr;
    if (colorFlag){
      damode |= CGO_COLOR_ARRAY;
    }
    vertexVals = CGODrawArrays(I, GL_TRIANGLE_STRIP, damode, nverts);
    normalVals = vertexVals + (nverts*3);
    if (colorFlag){
      colorVals = normalVals + (nverts*3);
    }
    for(c = nEdge; c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];
      
      v[3] = vv1[0] + v[0] * tube_size;
      v[4] = vv1[1] + v[1] * tube_size;
      v[5] = vv1[2] + v[2] * tube_size;
      
      v[6] = v[3] + d[0];
      v[7] = v[4] + d[1];
      v[8] = v[5] + d[2];
      
      tmp_ptr = v;
      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
      if(colorFlag){
	if (c1){
	  colorVals[plc++] = c1[0]; colorVals[plc++] = c1[1]; 
	  colorVals[plc++] = c1[2]; colorVals[plc++] = I->alpha;
	} else {
	  colorVals[plc++] = I->color[0]; colorVals[plc++] = I->color[1]; 
	  colorVals[plc++] = I->color[2]; colorVals[plc++] = I->alpha;
	}
      }
      tmp_ptr = v + 3;
      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
      pl += 3;
      tmp_ptr = &normalVals[pl-3];
      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
      if(colorFlag){
	colorVals[plc++] = c2[0]; colorVals[plc++] = c2[1]; 
	colorVals[plc++] = c2[2]; colorVals[plc++] = I->alpha;
      }
      tmp_ptr = v + 6;
      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
      pl += 3;
    }
  }
#else
  CGOBegin(I, GL_TRIANGLE_STRIP);
  for(c = nEdge; c >= 0; c--) {
    v[0] = p1[0] * x[c] + p2[0] * y[c];
    v[1] = p1[1] * x[c] + p2[1] * y[c];
    v[2] = p1[2] * x[c] + p2[2] * y[c];

    v[3] = vv1[0] + v[0] * tube_size;
    v[4] = vv1[1] + v[1] * tube_size;
    v[5] = vv1[2] + v[2] * tube_size;

    v[6] = v[3] + d[0];
    v[7] = v[4] + d[1];
    v[8] = v[5] + d[2];

    CGONormalv(I, v);
    if(colorFlag)
      CGOColorv(I, c1);
    CGOVertexv(I, v + 3);
    if(colorFlag)
      CGOColorv(I, c2);
    CGOVertexv(I, v + 6);
  }
  CGOEnd(I);
#endif
  if(cap1) {
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

    if(colorFlag && c1)
      CGOColorv(I, c1);
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = nEdge + 2, pl = 0, damode = CGO_VERTEX_ARRAY;
      float *vertexVals;
      float *normalVals = 0;
      if(cap2 == cCylCapRound){
	damode |= CGO_NORMAL_ARRAY;
      } else {
	CGONormalv(I, v);
      }
      vertexVals = CGODrawArrays(I, GL_TRIANGLE_FAN, damode, nverts);
      if(cap1 == cCylCapRound){
	normalVals = vertexVals + (nverts*3);
	normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
      }
      vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
      pl += 3;
      for(c = nEdge; c >= 0; c--) {
	v[0] = p1[0] * x[c] + p2[0] * y[c];
	v[1] = p1[1] * x[c] + p2[1] * y[c];
	v[2] = p1[2] * x[c] + p2[2] * y[c];
	
	v[3] = vv1[0] + v[0] * tube_size;
	v[4] = vv1[1] + v[1] * tube_size;
	v[5] = vv1[2] + v[2] * tube_size;
	
	if(normalVals){
	  normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
	}
	vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
	pl += 3;
      }
    }
#else
    CGOBegin(I, GL_TRIANGLE_FAN);
    CGONormalv(I, v);
    CGOVertexv(I, v + 3);

    for(c = nEdge; c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv1[0] + v[0] * tube_size;
      v[4] = vv1[1] + v[1] * tube_size;
      v[5] = vv1[2] + v[2] * tube_size;

      if(cap1 == cCylCapRound)
        CGONormalv(I, v);
      CGOVertexv(I, v + 3);
    }
    CGOEnd(I);
#endif
  }

  if(cap2) {
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

    if(colorFlag)
      CGOColorv(I, c2);
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = nEdge + 2, pl = 0, damode = CGO_VERTEX_ARRAY;
      float *vertexVals;
      float *normalVals = 0;
      if(cap2 == cCylCapRound){
	damode |= CGO_NORMAL_ARRAY;
      } else {
	CGONormalv(I, v);
      }
      vertexVals = CGODrawArrays(I, GL_TRIANGLE_FAN, damode, nverts);
      if(cap2 == cCylCapRound){
	normalVals = vertexVals + (nverts*3);
	normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
      }
      vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
      pl += 3;
      for(c = 0; c <= nEdge; c++) {
	v[0] = p1[0] * x[c] + p2[0] * y[c];
	v[1] = p1[1] * x[c] + p2[1] * y[c];
	v[2] = p1[2] * x[c] + p2[2] * y[c];
	
	v[3] = vv2[0] + v[0] * tube_size;
	v[4] = vv2[1] + v[1] * tube_size;
	v[5] = vv2[2] + v[2] * tube_size;
	
	if(normalVals){
	  normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
	}
	vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
	pl += 3;
      }
    }
#else
    CGOBegin(I, GL_TRIANGLE_FAN);
    CGONormalv(I, v);
    CGOVertexv(I, v + 3);

    for(c = 0; c <= nEdge; c++) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv2[0] + v[0] * tube_size;
      v[4] = vv2[1] + v[1] * tube_size;
      v[5] = vv2[2] + v[2] * tube_size;

      if(cap2 == cCylCapRound)
        CGONormalv(I, v);
      CGOVertexv(I, v + 3);
    }
    CGOEnd(I);
#endif
  }
}

static void CGOSimpleCone(CGO * I, float *v1, float *v2, float r1, float r2,
                          float *c1, float *c2, int cap1, int cap2)
{
#define MAX_EDGE 50

  float d[3], t[3], p0[3], p1[3], p2[3], vv1[3], vv2[3], v_buf[9], *v;
  float x[MAX_EDGE + 1], y[MAX_EDGE + 1], edge_normal[3 * (MAX_EDGE + 1)];
#if 0
  float overlap1, overlap2;
#endif
  float nub1, nub2;
  int colorFlag;
  int nEdge;
  int c;

  v = v_buf;
  nEdge = (int) SettingGet(I->G, cSetting_cone_quality);
#if 0
  overlap1 = r1 * SettingGet(I->G, cSetting_stick_overlap);
  overlap2 = r2 * SettingGet(I->G, cSetting_stick_overlap);
#endif

  nub1 = r1 * SettingGet(I->G, cSetting_stick_nub);
  nub2 = r2 * SettingGet(I->G, cSetting_stick_nub);

  if(nEdge > MAX_EDGE)
    nEdge = MAX_EDGE;
  subdivide(nEdge, x, y);

  colorFlag = (c1 != c2) && c2;

  CGOColorv(I, c1);

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  normalize3f(p0);

#if 0
  if(cap1 == cCylCapRound)
    ) {
    vv1[0] = v1[0] - p0[0] * overlap1;
    vv1[1] = v1[1] - p0[1] * overlap1;
    vv1[2] = v1[2] - p0[2] * overlap1;
  } else
#endif

  {
    vv1[0] = v1[0];
    vv1[1] = v1[1];
    vv1[2] = v1[2];
  }

#if 0
  if(cap2 == cCylCapRound) {
    vv2[0] = v2[0] + p0[0] * overlap2;
    vv2[1] = v2[1] + p0[1] * overlap2;
    vv2[2] = v2[2] + p0[2] * overlap2;
  } else
#endif

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
#ifdef _PYMOL_CGO_DRAWARRAYS
  {
    int nverts = (2*nEdge) + 2, pl = 0, plc = 0, damode = CGO_VERTEX_ARRAY | CGO_NORMAL_ARRAY;
    float *vertexVals;
    float *normalVals, *colorVals = 0, *tmp_ptr;
    if (colorFlag){
      damode |= CGO_COLOR_ARRAY;
    }
    vertexVals = CGODrawArrays(I, GL_TRIANGLE_STRIP, damode, nverts);
    normalVals = vertexVals + (nverts*3);
    if (colorFlag){
      colorVals = normalVals + (nverts*3);
    }
    for(c = nEdge; c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];
      
      v[3] = vv1[0] + v[0] * r1;
      v[4] = vv1[1] + v[1] * r1;
      v[5] = vv1[2] + v[2] * r1;

      v[6] = vv1[0] + v[0] * r2 + d[0];
      v[7] = vv1[1] + v[1] * r2 + d[1];
      v[8] = vv1[2] + v[2] * r2 + d[2];
      
      tmp_ptr = edge_normal + 3 * c;
      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
      if(colorFlag){
	colorVals[plc++] = c1[0]; colorVals[plc++] = c1[1]; colorVals[plc++] = c1[2]; colorVals[plc++] = I->alpha;
      }
      tmp_ptr = v + 3;
      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
      pl += 3;
      tmp_ptr = &normalVals[pl-3];
      normalVals[pl] = tmp_ptr[0]; normalVals[pl+1] = tmp_ptr[1]; normalVals[pl+2] = tmp_ptr[2];
      if(colorFlag){
	colorVals[plc++] = c2[0]; colorVals[plc++] = c2[1]; colorVals[plc++] = c2[2]; colorVals[plc++] = I->alpha;
      }
      tmp_ptr = v + 6;
      vertexVals[pl] = tmp_ptr[0]; vertexVals[pl+1] = tmp_ptr[1]; vertexVals[pl+2] = tmp_ptr[2];
      pl += 3;
    }
  }
#else
  CGOBegin(I, GL_TRIANGLE_STRIP);
  for(c = nEdge; c >= 0; c--) {
    v[0] = p1[0] * x[c] + p2[0] * y[c];
    v[1] = p1[1] * x[c] + p2[1] * y[c];
    v[2] = p1[2] * x[c] + p2[2] * y[c];

    v[3] = vv1[0] + v[0] * r1;
    v[4] = vv1[1] + v[1] * r1;
    v[5] = vv1[2] + v[2] * r1;

    v[6] = vv1[0] + v[0] * r2 + d[0];
    v[7] = vv1[1] + v[1] * r2 + d[1];
    v[8] = vv1[2] + v[2] * r2 + d[2];

    CGONormalv(I, edge_normal + 3 * c);
    if(colorFlag)
      CGOColorv(I, c1);
    CGOVertexv(I, v + 3);
    if(colorFlag)
      CGOColorv(I, c2);
    CGOVertexv(I, v + 6);
  }
  CGOEnd(I);
#endif

  if(cap1) {
    v[0] = -p0[0];
    v[1] = -p0[1];
    v[2] = -p0[2];

#if 0
    if(cap1 == cCylCapRound) {
      v[3] = vv1[0] - p0[0] * nub1;
      v[4] = vv1[1] - p0[1] * nub1;
      v[5] = vv1[2] - p0[2] * nub1;
    } else
#endif
    {
      v[3] = vv1[0];
      v[4] = vv1[1];
      v[5] = vv1[2];
    }

    if(colorFlag)
      CGOColorv(I, c1);
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = nEdge + 2, pl = 0, damode = CGO_VERTEX_ARRAY;
      float *vertexVals;
      float *normalVals = 0;
      if(cap2 == cCylCapRound){
	damode |= CGO_NORMAL_ARRAY;
      } else {
	CGONormalv(I, v);	
      }

      vertexVals = CGODrawArrays(I, GL_TRIANGLE_FAN, damode, nverts);
      if(cap2 == cCylCapRound){
	normalVals = vertexVals + (nverts*3);
	normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
      }
      vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
      pl += 3;
      for(c = nEdge; c >= 0; c--) {
	v[0] = p1[0] * x[c] + p2[0] * y[c];
	v[1] = p1[1] * x[c] + p2[1] * y[c];
	v[2] = p1[2] * x[c] + p2[2] * y[c];
	
	v[3] = vv1[0] + v[0] * r1;
	v[4] = vv1[1] + v[1] * r1;
	v[5] = vv1[2] + v[2] * r1;
	
	if(normalVals){
	  normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
	}
	vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
	pl += 3;
      }
    }
#else
    CGOBegin(I, GL_TRIANGLE_FAN);
    CGONormalv(I, v);
    CGOVertexv(I, v + 3);

    for(c = nEdge; c >= 0; c--) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv1[0] + v[0] * r1;
      v[4] = vv1[1] + v[1] * r1;
      v[5] = vv1[2] + v[2] * r1;

      if(cap1 == cCylCapRound)
        CGONormalv(I, v);
      CGOVertexv(I, v + 3);
    }
    CGOEnd(I);
#endif
  }

  if(cap2) {

    v[0] = p0[0];
    v[1] = p0[1];
    v[2] = p0[2];

#if 0
    if(cap2 == cCylCapRound) {
      v[3] = vv2[0] + p0[0] * nub2;
      v[4] = vv2[1] + p0[1] * nub2;
      v[5] = vv2[2] + p0[2] * nub2;
    } else
#endif
    {
      v[3] = vv2[0];
      v[4] = vv2[1];
      v[5] = vv2[2];
    }

    if(colorFlag)
      CGOColorv(I, c2);
#ifdef _PYMOL_CGO_DRAWARRAYS
    {
      int nverts = nEdge + 2, pl = 0, damode = CGO_VERTEX_ARRAY;
      float *vertexVals;
      float *normalVals = 0;
      if(cap2 == cCylCapRound){
	damode |= CGO_NORMAL_ARRAY;
      } else {
	CGONormalv(I, v);
      }
      vertexVals = CGODrawArrays(I, GL_TRIANGLE_FAN, damode, nverts);
      if(cap2 == cCylCapRound){
	normalVals = vertexVals + (nverts*3);
	normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
      }
      vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
      pl += 3;
      for(c = 0; c <= nEdge; c++) {
	v[0] = p1[0] * x[c] + p2[0] * y[c];
	v[1] = p1[1] * x[c] + p2[1] * y[c];
	v[2] = p1[2] * x[c] + p2[2] * y[c];
	
	v[3] = vv2[0] + v[0] * r2;
	v[4] = vv2[1] + v[1] * r2;
	v[5] = vv2[2] + v[2] * r2;
	
	if(normalVals){
	  normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
	}
	vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
	pl += 3;
      }
    }
#else
    CGOBegin(I, GL_TRIANGLE_FAN);
    CGONormalv(I, v);
    CGOVertexv(I, v + 3);

    for(c = 0; c <= nEdge; c++) {
      v[0] = p1[0] * x[c] + p2[0] * y[c];
      v[1] = p1[1] * x[c] + p2[1] * y[c];
      v[2] = p1[2] * x[c] + p2[2] * y[c];

      v[3] = vv2[0] + v[0] * r2;
      v[4] = vv2[1] + v[1] * r2;
      v[5] = vv2[2] + v[2] * r2;

      if(cap2 == cCylCapRound)
        CGONormalv(I, v);
      CGOVertexv(I, v + 3);
    }
    CGOEnd(I);
#endif
  }
}

/* CGOGetNextDrawBufferedIndex: This is used by RepSurface to */
/* get the data from the CGO_DRAW_BUFFERS_INDEXED operation so */
/* that it can update the indices for semi-transparent surfaces. */
float *CGOGetNextDrawBufferedIndex(float *cgo_op)
{
  register float *pc = cgo_op;
  int op;

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	//	int nverts = CGO_get_int(pc + 4);
	//	pc += nverts*3 + 9 ;
	return pc;
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (0);
}

float *CGOGetNextOp(float *cgo_op, int optype)
{
  register float *pc = cgo_op;
  int op;

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    if (op==optype)
      return pc;
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (0);
}

int CGOGetSizeWithoutStops(CGO *I){
  register float *pc = I->op;
  int op;
  while((op = (CGO_MASK & CGO_get_int(pc)))) {
    (void)CGO_read_int(pc);
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (pc - I->op);
}

void CGOAppendImpl(CGO *dest, CGO *source, int stopAtEnd);

void CGOAppendNoStop(CGO *dest, CGO *source){
  CGOAppendImpl(dest, source, 0);
}

void CGOAppend(CGO *dest, CGO *source){
  CGOAppendImpl(dest, source, 1);
}

void CGOAppendImpl(CGO *dest, CGO *source, int stopAtEnd){
  register float *pc = source->op;
  int sz, szd;
  register float *nc;

  sz = CGOGetSizeWithoutStops(source);
  szd = dest->c; 
  if (szd && !(CGO_MASK & CGO_get_int(dest->op + szd-1 ))){
    szd = CGOGetSizeWithoutStops(dest);
  }
  VLASizeForSure(dest->op, float, szd + sz);
  //  VLACheck(dest->op, float, szd + sz);
  dest->c = szd + sz;
  nc = dest->op + szd;
  while(sz--)
    *(nc++) = *(pc++);
  if (stopAtEnd)
    CGOStop(dest);
  return;
}

int CGOCountNumberOfOperationsOfTypeDEBUG(CGO *I, int optype){
  register float *pc = I->op;
  int op, numops = 0, totops = 0;
  if (!optype){
    printf("CGOCountNumberOfOperationsOfType: ");
  }
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    if (!optype){
      printf(" %02X ", op);
    }
    totops++;
    if (op == optype)
      numops++;
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  if (!optype){
    printf("\n");
  }
  //  printf("\n\ttotops=%d\n", totops);
  if(optype){
    return (numops);
  } else {
    return (totops);
  }
}

int CGOCountNumberOfOperationsOfType(CGO *I, int optype){
  register float *pc = I->op;
  int op, numops = 0, totops = 0;
#ifdef DEBUG_PRINT_OPS
  if (!optype){
    printf("CGOCountNumberOfOperationsOfType: ");
  }
#endif
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
#ifdef DEBUG_PRINT_OPS
    if (!optype){
      printf(" 0x%02X ", op);
    }
#endif
    totops++;
    if (op == optype)
      numops++;
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
#ifdef DEBUG_PRINT_OPS
  if (!optype){
    printf("\n");
  }
#endif
  //  printf("\n\ttotops=%d\n", totops);
  if(optype){
    return (numops);
  } else {
    return (totops);
  }
}

short CGOHasOperationsOfType(CGO *I, int optype){
  register float *pc = I->op;
  int op;

  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    //    printf("%X ", op);
    if (op == optype){
      return (1);
    }
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (0);
}
short CGOCheckWhetherToFree(PyMOLGlobals * G, CGO *I){
  if (I->use_shader){
    if (I->cgo_shader_ub_color != SettingGet(G, cSetting_cgo_shader_ub_color) || 
	I->cgo_shader_ub_normal != SettingGet(G, cSetting_cgo_shader_ub_normal)){
      return true;
    }
  }
  return false;
}
CGO *CGOConvertLinesToShaderCylinders(CGO * I, int est){
#ifdef _PYMOL_CGO_DRAWARRAYS
  CGO *cgo;

  register float *pc = I->op;
  register float *nc;
  register int op;
  float *save_pc;
  int sz, tot_nverts = 0, tot_ncyls = 0;

  cgo = CGONewSized(I->G, I->c + est);
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    save_pc = pc;
    switch (op) {
    case CGO_DRAW_ARRAYS:
      {
	int mode = CGO_get_int(pc), arrays = CGO_get_int(pc + 1), narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3);
	GLfloat *vals = CGODrawArrays(cgo, mode, arrays, nverts);
	int nvals = narrays*nverts, onvals ;
	onvals = nvals;
	pc += 4;
	while(nvals--)
	  *(vals++) = *(pc++);
	save_pc += onvals + 4 ;
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
	float *last_vertex = NULL, *last_color = NULL, *current_color = NULL, *color = NULL ;
	int nverts = 0, err = 0, end = 0;
	int mode = CGO_read_int(pc);
	while(!err && !end && (op = (CGO_MASK & CGO_read_int(pc)))) {
	  short add_to_cgo = 1;
	  end = (op == CGO_END);
	  switch (op) {
	  case CGO_VERTEX:
	    if (last_vertex){
	      switch (mode){
	      case GL_LINES:
	      case GL_LINE_STRIP:
		{
		  float axis[3];
		  axis[0] = pc[0] - last_vertex[0];
		  axis[1] = pc[1] - last_vertex[1];
		  axis[2] = pc[2] - last_vertex[2];
		  if (last_color && current_color && 
		      (last_color[0]!=current_color[0] || 
		       last_color[1]!=current_color[1] ||
		       last_color[2]!=current_color[2])){
		    CGOColorv(cgo, last_color);
		    CGOShaderCylinder2ndColor(cgo, last_vertex, axis, 1.f, 15, current_color);		    
		    CGOColorv(cgo, current_color);
		  } else {
		    CGOShaderCylinder(cgo, last_vertex, axis, 1.f, 15);
		  }
		  last_vertex = pc;
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
	      //	      current_color = yellow;
	      //	      CGOColorv(cgo, yellow);
	    }
	    nverts++;
	    add_to_cgo = 0;
	  case CGO_END:
	    if (op == CGO_END){
	      switch (mode){
	      case GL_LINES:
	      case GL_LINE_STRIP:
		add_to_cgo = 0;		
		break;
	      default:
		add_to_cgo = 1;
	      }
	    }
	  case CGO_COLOR:
	    if (op == CGO_COLOR){
	      last_color = current_color;
	      current_color = pc;
	      color = pc;
	    }
	  default:
	    sz = CGO_sz[op];
	    if (add_to_cgo){
	      nc = CGO_add(cgo, sz + 1);
	      *(nc++) = *(pc - 1);
	      while(sz--)
		*(nc++) = *(pc++);
	    } else {
	      pc += sz;
	    }
	  }
	  if (end){
	    break;
	  }
	}
	tot_nverts += nverts;
	save_pc = pc;
      }
      break;
    case CGO_ALPHA:
      I->alpha = *pc;
    default:
      sz = CGO_sz[op];
      nc = CGO_add(cgo, sz + 1);
      *(nc++) = *(pc - 1);
      while(sz--)
        *(nc++) = *(pc++);
    }
    pc = save_pc;
    pc += CGO_sz[op];
  }
  CGOStop(cgo);
  cgo->use_shader = I->use_shader;
  if (cgo->use_shader){
    cgo->cgo_shader_ub_color = SettingGet(cgo->G, cSetting_cgo_shader_ub_color);
    cgo->cgo_shader_ub_normal = SettingGet(cgo->G, cSetting_cgo_shader_ub_normal);
  }
  if (tot_nverts){
    return (cgo);
  } else {
    return NULL;
  }
#else
  return NULL;
#endif

}

int CGOCountNumberCustomCylinders(CGO *I, int *has_2nd_color){
  register float *pc = I->op;
  int op, numops = 0, totops = 0;
  *has_2nd_color = 0;
  while((op = (CGO_MASK & CGO_read_int(pc)))) {
    totops++;
    if (op == CGO_CUSTOM_CYLINDER){
      numops++;
      if (has_2nd_color){
	if (*(pc+7) != *(pc+10) || *(pc+8) != *(pc+11) || *(pc+9) != *(pc+12)){
	  (*has_2nd_color)++;
	}
      }
    }
    switch (op) {
#ifdef _PYMOL_CGO_DRAWARRAYS
    case CGO_DRAW_ARRAYS:
      {
	int narrays = CGO_get_int(pc + 2), nverts = CGO_get_int(pc + 3), floatlength = narrays*nverts;
	pc += floatlength + 4 ;
      }
      break;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
    case CGO_DRAW_BUFFERS_INDEXED:
      {
	int nverts = CGO_get_int(pc + 4);
	pc += nverts*3 + 10 ; 
      }
      break;
    case CGO_DRAW_BUFFERS_NOT_INDEXED:
      {
	int nverts = CGO_get_int(pc + 3);
	pc += nverts*3 + 8 ; 
      }
      break;
#endif
    }
    pc += CGO_sz[op];
  }
  return (totops);
}
