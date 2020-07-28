
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

#include"Base.h"
#include"MemoryDebug.h"
#include"OOMac.h"
#include"RepSurface.h"
#include"Map.h"
#include"Scene.h"
#include"Sphere.h"
#include"Setting.h"
#include"Color.h"
#include"ObjectMolecule.h"
#include"Triangle.h"
#include"Vector.h"
#include"Feedback.h"
#include"main.h"
#include"Util.h"
#include"CGO.h"
#include"P.h"
#include"PConv.h"
#include"Selector.h"
#include"ShaderMgr.h"

#ifdef NT
#undef NT
#endif

struct RepSurface : Rep {
  using Rep::Rep;

  ~RepSurface() override;

  cRep_t type() const override { return cRepSurface; }
  void render(RenderInfo* info) override;
  void invalidate(cRepInv_t level) override;
  Rep* recolor() override;
  bool sameVis() const override;
  bool sameColor() const override;

  int N = 0;
  int NT = 0; //!< number of triangles (?)
  bool proximity = false; //!< cSetting_surface_proximity
  float* V = nullptr;
  float* VN = nullptr;  //!< Normals
  float* VC = nullptr;  //!< Colors
  float* VA = nullptr;  //!< Alpha
  float* VAO = nullptr; //!< Ambient Occlusion per vertex
  int* RC = nullptr;
  int* Vis = nullptr;
  int* T = nullptr;  //!< vertices (of triangles?)
  int* S = nullptr;  //!< strips
  int* AT = nullptr; //!< closest atom for vertices
  bool oneColorFlag = false;
  int oneColor;
  bool allVisibleFlag = false;
  char* LastVisib = nullptr;
  int* LastColor = nullptr;
  bool ColorInvalidated = false;
  int Type;
  float max_vdw;

  /* These variables are for using the shader.  All of them */
  /* are allocated/set when generate_shader_cgo to minimize */
  /* allocation during the rendering loop. */
  CGO *shaderCGO = nullptr;
  CGO *pickingCGO = nullptr;
  bool dot_as_spheres = false;

  int surface_mode = cRepSurface_by_flags;
};

static
void RepSurfaceSmoothEdges(RepSurface * I);

static void setShaderCGO(RepSurface * I, CGO * cgo) {
  if (I->shaderCGO != I->pickingCGO) {
    CGOFree(I->shaderCGO);
  }
  I->shaderCGO = cgo;
}

static void setPickingCGO(RepSurface * I, CGO * cgo) {
  if (I->shaderCGO != I->pickingCGO) {
    CGOFree(I->pickingCGO);
  }
  I->pickingCGO = cgo;
}

RepSurface::~RepSurface()
{
  auto I = this;
  VLAFreeP(I->V);
  VLAFreeP(I->VN);
  setPickingCGO(I, NULL);
  setShaderCGO(I, NULL);
  FreeP(I->VC);
  FreeP(I->VA);
  if (I->VAO){
    VLAFreeP(I->VAO);
    I->VAO = 0;
  }
  FreeP(I->RC);
  FreeP(I->Vis);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  VLAFreeP(I->T);
  VLAFreeP(I->S);
  VLAFreeP(I->AT);
}

typedef struct {
  int nDot;
  float *dot;
  float *dotNormal;
  int *dotCode;
} SolventDot;

typedef struct {
  float vdw;
  int flags;
} SurfaceJobAtomInfo;

static SolventDot *SolventDotNew(PyMOLGlobals * G,
                                 float *coord,
                                 SurfaceJobAtomInfo * atom_info,
                                 float probe_radius, SphereRec * sp,
                                 int *present,
                                 int circumscribe, int surface_mode,
                                 int surface_solvent, int cavity_cull,
                                 int all_visible_flag, float max_vdw,
                                 int cavity_mode, float cavity_radius, 
                                 float cavity_cutoff);

static void SolventDotFree(SolventDot * I)
{
  if(I) {
    VLAFreeP(I->dot);
    VLAFreeP(I->dotNormal);
    VLAFreeP(I->dotCode);
  }
  OOFreeP(I);
}

#ifndef PURE_OPENGL_ES_2
static
void immediate_draw_masked_vertices(
    const float * vc,   // colors
    const float * vn,   // normals
    const float * v,    // vertices
    const int * mask,   // mask
    int count)
{
  for (int i = 0; i < count; ++i) {
    if (!mask[i])
      continue;
    int i3 = i * 3;
    if (vc) glColor3fv(vc + i3);
    if (vn) glNormal3fv(vn + i3);
    glVertex3fv(v + i3);
  }
}

static
void immediate_draw_indexed_vertices(
    const float * vc,     // colors
    const float * vn,     // normals
    const float * v,      // vertices
    const int * indices,  // indices
    int count)
{
  for (int i = 0; i < count; ++i) {
    int i3 = indices[i] * 3;
    if (vc) glColor3fv(vc + i3);
    if (vn) glNormal3fv(vn + i3);
    glVertex3fv(v + i3);
  }
}

static
void immediate_draw_indexed_vertices_alpha(
    const float * vc,     // colors
    const float * va,     // alpha array
    float alpha,          // alpha value if va is NULL
    const float * vn,     // normals
    const float * v,      // vertices
    const int * indices,  // indices
    int count)
{
  for (int i = 0; i < count; ++i) {
    int i3 = indices[i] * 3;
    if (vc)
      glColor4f(vc[i3], vc[i3 + 1], vc[i3 + 2],
          va ? va[indices[i]] : alpha);
    if (vn) glNormal3fv(vn + i3);
    glVertex3fv(v + i3);
  }
}
#endif

static
bool visibility_test(
    bool proximityFlag,
    const int * vi,     // visibility array
    const int * t)      // indices
{
  if (proximityFlag)
    return (vi[t[0]] || vi[t[1]] || vi[t[2]]);
  return   (vi[t[0]] && vi[t[1]] && vi[t[2]]);
}

static int check_and_add(int *cache, int spacing, int t0, int t1)
{
  int *rec;
  int cnt;
  t0++;
  t1++;

  rec = cache + spacing * t0;
  cnt = spacing;
  while(cnt > 0) {
    if(*rec == t1)
      return 1;
    if(!*rec) {
      *rec = t1;
      break;
    }
    rec++;
    cnt--;
  }
  rec = cache + spacing * t1;
  cnt = spacing;
  while(cnt > 0) {
    if(*rec == t0)
      return 1;
    if(!*rec) {
      *rec = t0;
      break;
    }
    rec++;
    cnt--;
  }
  return 0;
}

#define CLAMP_VALUE(v)  ((v>1.f) ? 1.f :  (v < 0.f) ? 0.f : v)

static int AtomInfoIsMasked(ObjectMolecule *obj, int atm){
  AtomInfoType *ait;
  if (atm < 0)
    return cPickableNoPick;
  ait = &obj->AtomInfo[atm];
  return (ait->masked ? cPickableNoPick : cPickableAtom);
}

static int RepSurfaceCGOGenerate(RepSurface * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;
  float *v = I->V;
  float *vn = I->VN;
  float *vc = I->VC;
  float *va = I->VA;
  int *t = I->T;
  int *s = I->S;
  int c = I->N;
  int *vi = I->Vis;
  int *at = I->AT;
  int ok = true;
  float alpha;
  int t_mode;
  CGO *convertcgo = NULL;

  auto* const cs = I->cs;
  auto* const obj = I->cs->Obj;

  bool pick_surface = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_pick_surface);
  short dot_as_spheres = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_as_spheres);
  alpha = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;

  setPickingCGO(I, nullptr);
  setShaderCGO(I, CGONew(G));

  if (!I->shaderCGO)
    return false;

  I->shaderCGO->use_shader = true;
  I->dot_as_spheres = dot_as_spheres;

  if (alpha < 1) {
    I->setHasTransparency();
  }

  if (I->Type == 1) {
    /* no triangle information, so we're rendering dots only */
    int normals =
        SettingGet<int>(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_normals);
    if(!normals){
      CGOResetNormal(I->shaderCGO, true);
    }
      if((alpha != 1.0)) {
	CGOAlpha(I->shaderCGO, alpha);
      }
      if (dot_as_spheres){
	if(c) {
	  ok &= CGOColor(I->shaderCGO, 1.0, 0.0, 0.0);
	  if(ok && I->oneColorFlag) {
	    ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
	  }
	  while(ok && c--) {
	    if(*vi) {
	      if(!I->oneColorFlag) {
		ok &= CGOColorv(I->shaderCGO, vc);
	      }
	      if(ok && normals)
		ok &= CGONormalv(I->shaderCGO, vn);
	      if (ok && pick_surface)
		ok &= CGOPickColor(I->shaderCGO, *at, AtomInfoIsMasked(obj, *at));

	      if (ok)
		ok &= CGOSphere(I->shaderCGO, v, 1.f);
	    }
	    vi++;
	    vc += 3;
	    vn += 3;
	    v += 3;
	    at++;
	  }
	}
      } else {
	ok &= CGODotwidth(I->shaderCGO, SettingGet_f
			  (G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_width));
	if(ok && c) {
	  ok &= CGOColor(I->shaderCGO, 1.0, 0.0, 0.0);
	  ok &= CGOBegin(I->shaderCGO, GL_POINTS);
	  if(ok && I->oneColorFlag) {
	    ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
	  }
	      
	  while(ok && c--) {
	    if(*vi) {
	      if(!I->oneColorFlag) {
		ok &= CGOColorv(I->shaderCGO, vc);
	      }
	      if(normals){
		CGONormalv(I->shaderCGO, vn);
	      }
	      if (ok && pick_surface)
		ok &= CGOPickColor(I->shaderCGO, *at, AtomInfoIsMasked(obj, *at));
	      if (ok)
		ok &= CGOVertexv(I->shaderCGO, v);
	    }
	    vi++;
	    vc += 3;
	    vn += 3;
	    v += 3;
	    at++;
	  }
	  if (ok)
	    ok &= CGOEnd(I->shaderCGO);
	}
      }
  } else if (I->Type == 2) { /* rendering triangle mesh */
    int normals =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_normals);
    if(ok && !normals){
	ok &= CGOResetNormal(I->shaderCGO, true);
    }
    if (ok) {
      float line_width =
	SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_width);
          line_width = SceneGetDynamicLineWidth(info, line_width);

	ok &= CGOSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_MESH);

          c = I->NT;
          if(ok && c) {
            if(I->oneColorFlag) {
		  ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		  while(ok && c--) {
              if (visibility_test(I->proximity, vi, t)) {
		      if(normals) {
			int idx;
			ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
			idx = (*(t + 2));
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);

			idx = (*t);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			if (ok)
			  ok &= CGOEnd(I->shaderCGO);
		      } else {
			int idx;
			ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);

			idx = (*(t + 2));
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			idx = *t;
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = *t;
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = *t;
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			if (ok)
			  ok &= CGOEnd(I->shaderCGO);
		      }
                } else
                  t += 3;
              }
            } else { /* not oneColorFlag */
		  while(ok && c--) {
              if (visibility_test(I->proximity, vi, t)) {
		      if(normals) {
			int idx;
			ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
			
			idx = (*(t + 2));
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			if (ok)
			  ok &= CGOEnd(I->shaderCGO);
		      } else {
			int idx;
			ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);

			idx = (*(t + 2));
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);

			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			if (ok)
			  ok &= CGOEnd(I->shaderCGO);
		      }
		    } else
		      t += 3;
	    }
          }
      }
    }
  } else {
    /* we're rendering triangles */
    if((alpha != 1.0) || va) {
      t_mode =
	SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_transparency_mode);

          if(info && info->alpha_cgo) {
            t_mode = 0;
          }

          if(t_mode) {

            float **t_buf = NULL, **tb;
            float *z_value = NULL, *zv;
            int *ix = NULL;
	int n_tri = 0;
            float sum[3];
            float matrix[16];

            glGetFloatv(GL_MODELVIEW_MATRIX, matrix);

            if(I->oneColorFlag) {
	  t_buf = pymol::malloc<float *>(I->NT * 6);
            } else {
	  t_buf = pymol::malloc<float *>(I->NT * 12);
            }
	    CHECKOK(ok, t_buf);
	    if (ok){
	  z_value = pymol::malloc<float>(I->NT);
	      CHECKOK(ok, z_value);
	    }
	    if (ok){
	  ix = pymol::malloc<int>(I->NT);
	      CHECKOK(ok, ix);
	    }
            zv = z_value;
            tb = t_buf;
            c = I->NT;
	    if (ok){
	      if(I->oneColorFlag) {
		while(c--) {
              if (visibility_test(I->proximity, vi, t)) {
		    *(tb++) = vn + (*t) * 3;
		    *(tb++) = v + (*t) * 3;
		    *(tb++) = vn + (*(t + 1)) * 3;
		    *(tb++) = v + (*(t + 1)) * 3;
		    *(tb++) = vn + (*(t + 2)) * 3;
		    *(tb++) = v + (*(t + 2)) * 3;
		    
		    add3f(tb[-1], tb[-3], sum);
		    add3f(sum, tb[-5], sum);
		    
		    *(zv++) = matrix[2] * sum[0] + matrix[6] * sum[1] + matrix[10] * sum[2];
		    n_tri++;
		  }
		  t += 3;
		}
	      } else {
		while(c--) {
              if (visibility_test(I->proximity, vi, t)) {
		      
		      if(va)
			*(tb++) = va + (*t);
		      else
			*(tb++) = &alpha;
		      *(tb++) = vc + (*t) * 3;
		      *(tb++) = vn + (*t) * 3;
		      *(tb++) = v + (*t) * 3;
		      
		      if(va)
			*(tb++) = va + (*(t + 1));
		      else
			*(tb++) = &alpha;
		      *(tb++) = vc + (*(t + 1)) * 3;
		      *(tb++) = vn + (*(t + 1)) * 3;
		      *(tb++) = v + (*(t + 1)) * 3;
		      
		      if(va)
			*(tb++) = va + (*(t + 2));
		      else
			*(tb++) = &alpha;
		      *(tb++) = vc + (*(t + 2)) * 3;
		      *(tb++) = vn + (*(t + 2)) * 3;
		      *(tb++) = v + (*(t + 2)) * 3;
		      
		      add3f(tb[-1], tb[-5], sum);
		      add3f(sum, tb[-9], sum);
		      *(zv++) =
			matrix[2] * sum[0] + matrix[6] * sum[1] + matrix[10] * sum[2];
		      n_tri++;
		    }
		  t += 3;
		}
	      }
	    }
            switch (t_mode) {
            case 1:
              ok &= UtilSemiSortFloatIndex(n_tri, z_value, ix, true);
              /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn); */
              break;
            default:
              ok &= UtilSemiSortFloatIndex(n_tri, z_value, ix, false);
              /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
              break;
            }
            c = n_tri;
            if(I->oneColorFlag) {
              float col[3];
              ColorGetEncoded(G, I->oneColor, col);
	      if (ok){
		  CGOAlpha(I->shaderCGO, alpha);
		  if (ok)
		    ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
		  if (ok)
		    ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
		  for(c = 0; ok && c < n_tri; c++) {
		    int idx;
		    tb = t_buf + 6 * c;
		    //		    tb = t_buf + 6 * ix[c];
		  if (ok)
		    ok &= CGONormalv(I->shaderCGO, *(tb++));
		  idx = ((*tb - v)/3);
		if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok)
		      ok &= CGOVertexv(I->shaderCGO, *(tb++));
		    if (ok)
		      ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
		if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok)
		      ok &= CGOVertexv(I->shaderCGO, *(tb++));
		    if (ok)
		      ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
		if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
		    if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok)
		      ok &= CGOVertexv(I->shaderCGO, *(tb++));
		  }
		  if (ok)
		    ok &= CGOEnd(I->shaderCGO);
	      }
            } else { /* else I->oneColorFlag */
		  if (ok)
		    ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
		  for(c = 0; ok && c < n_tri; c++) {
		    float *vv, *v_alpha;
		    int idx;
		    tb = t_buf + 12 * c; /* need to update index every frame */
		    //		    tb = t_buf + 12 * ix[c];
		    v_alpha = *(tb++);
		    vv = *(tb++);
		    ok &= CGOAlpha(I->shaderCGO, *v_alpha);
		    if (ok) ok &= CGOColor(I->shaderCGO, vv[0], vv[1], vv[2]);
		    if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
	      if (ok && pick_surface)
		ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
		    
		    v_alpha = *(tb++);
		    vv = *(tb++);
		    if (ok) ok &= CGOAlpha(I->shaderCGO, *v_alpha);
		    if (ok) ok &= CGOColor(I->shaderCGO, vv[0], vv[1], vv[2]);
		    if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
	      if (ok && pick_surface)
		ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
		    
		    v_alpha = *(tb++);
		    vv = *(tb++);
		    if (ok) ok &= CGOAlpha(I->shaderCGO, *v_alpha);
		    if (ok) ok &= CGOColor(I->shaderCGO, vv[0], vv[1], vv[2]);
		    if (ok) ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
	      if (ok && pick_surface)
		ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked(obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
		  }
		  if (ok) ok &= CGOEnd(I->shaderCGO);
		}
	      FreeP(ix);
	      FreeP(z_value);
	      FreeP(t_buf);
          } else if (ok) {
            if(info->alpha_cgo) {       /* global transparency sort */
              if(I->allVisibleFlag) {
                if(I->oneColorFlag) {
                  float col[3];
                  ColorGetEncoded(G, I->oneColor, col);

                  c = *(s++);
                  while(c) {
                    int parity = 0;
                    s += 2;
                    while(ok && c--) {
                      ok &= CGOAlphaTriangle(info->alpha_cgo,
					     v + s[-2] * 3, v + s[-1] * 3, v + (*s) * 3,
					     vn + s[-2] * 3, vn + s[-1] * 3, vn + (*s) * 3,
					     col, col, col, alpha, alpha, alpha, parity);
                      s++;
                      parity = !parity;
                    }
                    c = *(s++);
                  }
                } else {
                  c = *(s++);
                  while(ok && c) {
                    float *col0, *col1, *col2;
                    int parity = 0;
                    float alpha0, alpha1, alpha2;
                    col0 = vc + (*s) * 3;
                    if(va) {
                      alpha0 = va[(*s)];
                    } else {
                      alpha0 = alpha;
                    }
                    s++;
                    col1 = vc + (*s) * 3;
                    if(va) {
                      alpha1 = va[(*s)];
                    } else {
                      alpha1 = alpha;
                    }
                    s++;
                    while(ok && c--) {
                      col2 = vc + (*s) * 3;
                      if(va) {
                        alpha2 = va[(*s)];
                      } else {
                        alpha2 = alpha;
                      }
                      ok &= CGOAlphaTriangle(info->alpha_cgo,
					     v + s[-2] * 3, v + s[-1] * 3, v + (*s) * 3,
					     vn + s[-2] * 3, vn + s[-1] * 3, vn + (*s) * 3,
					     col0, col1, col2, alpha0, alpha1, alpha2, parity);
                      alpha0 = alpha1;
                      alpha1 = alpha2;
                      col0 = col1;
                      col1 = col2;
                      s++;
                      parity = !parity;
                    }
                    c = *(s++);
                  }
                }
              } else if (ok){          /* subset s */
                c = I->NT;
                if(c) {
                  if(I->oneColorFlag) {
                    float color[3];
                    ColorGetEncoded(G, I->oneColor, color);
                    while(ok && c--) {
                  if (visibility_test(I->proximity, vi, t)) {

                        ok &= CGOAlphaTriangle(info->alpha_cgo,
					       v + t[0] * 3, v + t[1] * 3, v + t[2] * 3,
					       vn + t[0] * 3, vn + t[1] * 3, vn + t[2] * 3,
					       color, color, color, alpha, alpha, alpha, 0);
                      }
                      t += 3;
                    }
                  } else {
                    while(ok && c--) {
                  if (visibility_test(I->proximity, vi, t)) {

                        if(va) {
                          ok &= CGOAlphaTriangle(info->alpha_cgo,
						 v + t[0] * 3, v + t[1] * 3, v + t[2] * 3,
						 vn + t[0] * 3, vn + t[1] * 3, vn + t[2] * 3,
						 vc + t[0] * 3, vc + t[1] * 3, vc + t[2] * 3,
						 va[t[0]], va[t[1]], va[t[2]], 0);
                        } else {
                          ok &= CGOAlphaTriangle(info->alpha_cgo,
						 v + t[0] * 3, v + t[1] * 3, v + t[2] * 3,
						 vn + t[0] * 3, vn + t[1] * 3, vn + t[2] * 3,
						 vc + t[0] * 3, vc + t[1] * 3, vc + t[2] * 3,
						 alpha, alpha, alpha, 0);
                        }
                      }
                      t += 3;
                    }
                  }
		  if (ok)
		    CGOEnd(info->alpha_cgo);
                }
              }
            } else if (ok){

              /* fast and ugly */
              /*          glCullFace(GL_BACK);
                 glEnable(GL_CULL_FACE);
                 glDepthMask(GL_FALSE); */
              if(I->allVisibleFlag) {
                if(I->oneColorFlag) {
                  float col[3];
                  ColorGetEncoded(G, I->oneColor, col);

		      if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
		      if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
		      c = *(s++);
		      while(ok && c) {
			if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
                  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
                  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			while(ok && c--) {
			  ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
                    if (ok && pick_surface)
		      ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		} else { /* I->oneColorFlag */
		      c = *(s++);
		      while(ok && c) {
			float *col;
			if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			col = vc + (*s) * 3;
			if(va) {
			  if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*s)]);
			  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
			} else {
			  if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
			  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
			}
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
                  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			col = vc + (*s) * 3;
			if(va) {
			  if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*s)]);
			  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
			} else {
			  if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
			  if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
			}
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
                  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			while(ok && c--) {
			  col = vc + (*s) * 3;
			  if(va) {
			    ok &= CGOAlpha(I->shaderCGO, va[(*s)]);
			    if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
			  } else {
			    ok &= CGOAlpha(I->shaderCGO, alpha);
			    if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
			  }
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
                    if (ok && pick_surface)
		      ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		  }
              } else {          /* subset s */
		    c = I->NT;
		    if(ok && c) {
		      if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
		      if(I->oneColorFlag) {
			float color[3];
			ColorGetEncoded(G, I->oneColor, color);
			if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
			if (ok) ok &= CGOColor(I->shaderCGO, color[0], color[1], color[2]);
			while(ok && c--) {
                    if (visibility_test(I->proximity, vi, t))
			    {
			      CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
                      if (ok && pick_surface)
			ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
                      if (ok && pick_surface)
			ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
                      if (ok && pick_surface)
			ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			    } else
			    t += 3;
			}
		      } else {
			float *col;
			while(ok && c--) {
                    if (visibility_test(I->proximity, vi, t))
			    {
			      
			      col = vc + (*t) * 3;
			      if(va) {
				ok &= CGOAlpha(I->shaderCGO, va[(*t)]);
			      } else {
				ok &= CGOAlpha(I->shaderCGO, alpha);
			      }
			      if (ok) ok &= CGOColorv(I->shaderCGO, col);
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
                      if (ok && pick_surface)
			ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			      col = vc + (*t) * 3;
			      if(va) {
				if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*t)]);
			      } else {
				if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
			      }
			      if (ok) ok &= CGOColorv(I->shaderCGO, col);
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
                      if (ok && pick_surface)
			ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			      col = vc + (*t) * 3;
			      if(va) {
				if (ok) ok &= CGOAlpha(I->shaderCGO, va[(*t)]);
			      } else {
				if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
			      }
			      if (ok) ok &= CGOColorv(I->shaderCGO, col);
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
                      if (ok && pick_surface)
			ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			    } else
			    t += 3;
			}
		      }
		      if (ok) ok &= CGOEnd(I->shaderCGO);
		    }
                }
            }
            /*          glDisable(GL_CULL_FACE);
               glDepthMask(GL_TRUE); */
          }
        } else if (ok) {                /* opaque */
            if(I->allVisibleFlag) {
              if(I->oneColorFlag) {
                if (ok) {
		      CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		      c = *(s++);
		      while(ok && c) {
			if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
                if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			while(ok && c--) {
			  ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
                  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		  }
	      } else {          /* not one color */
                    {
		      c = *(s++);
		      while(ok && c) {
			ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
                if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
                if (ok && pick_surface)
		  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			while(ok && c--) {
			  ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked(obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO );
			c = *(s++);
		      }
		    }
              }                 /* one color */
            } else if (ok) {            /* subsets */
		  c = I->NT;
		  if(ok && c) {
		    ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
		    if(I->oneColorFlag) {
		      if (ok) ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		      while(ok && c--) {
                if (visibility_test(I->proximity, vi, t)) {
			  ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			} else
			  t += 3;
		      }
		    } else {
		      while(ok && c--) {
                if (visibility_test(I->proximity, vi, t)) {
			  ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
		  if (ok && pick_surface)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			} else
			  t += 3;
		      }
		    }
		    if (ok) ok &= CGOEnd(I->shaderCGO);
		  }
            }
	  /*          if (use_shader) {
	    CShaderPrg_Disable(shaderPrg);
	    }*/
        }
      }
      if(ok && SettingGetGlobal_i(G, cSetting_surface_debug)) {
        t = I->T;
        c = I->NT;
	    if(c) {
	      ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
	      while(ok && c--) {
		if(I->allVisibleFlag
              || visibility_test(I->proximity, vi, t)) {
		  ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
		  }
	    if (ok && pick_surface)
	      ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
		  }
	    if (ok && pick_surface)
	      ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
		  }
	    if (ok && pick_surface)
	      ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		} else {
		  t += 3;
		}
	      }
	      if (ok) ok &= CGOEnd(I->shaderCGO);
	    }
        t = I->T;
        c = I->NT;

	    if(ok && c) {
	      ok &= CGOColor(I->shaderCGO, 0.0, 1.0, 0.0);
	      if (ok) ok &= CGODotwidth(I->shaderCGO, 1.0F);
	      while(ok && c--) {
		ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
		if(I->allVisibleFlag
              || visibility_test(I->proximity, vi, t)) {
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
	    if (ok && pick_surface)
	      ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
	    if (ok && pick_surface)
	      ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
	    if (ok && pick_surface)
	      ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked(obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		} else {
		  t += 3;
		}
		if (ok) ok &= CGOEnd(I->shaderCGO);
	      }
	    }

    c = I->N;
      if(ok && c) {
	ok &= CGOColor(I->shaderCGO, 1.0, 0.0, 0.0);
	if (ok) ok &= CGOResetNormal(I->shaderCGO, true);
	if (ok) ok &= CGOBegin(I->shaderCGO, GL_LINES);
	while(ok && c--) {
	  ok &= CGOVertexv(I->shaderCGO, v);
	  if (ok) ok &= CGOVertex(I->shaderCGO, v[0] + vn[0] / 2, v[1] + vn[1] / 2, v[2] + vn[2] / 2);
	  v += 3;
	  vn += 3;
	}
	if (ok) ok &= CGOEnd(I->shaderCGO);
      }
  }


  
  if (ok) ok &= CGOStop(I->shaderCGO);
  if (I->Type != 2){
    CGOCombineBeginEnd(&I->shaderCGO);
  }
  if(I->Type == 1){
    if (dot_as_spheres) {
      CGO *tmpCGO = CGONew(G);
      CHECKOK(ok, tmpCGO);
      ok &= CGOEnable(tmpCGO, GL_SPHERE_SHADER);
      ok &= CGOEnable(tmpCGO, GL_DOT_LIGHTING);  // TODO: this needs normals in sphere shader PYMOL-1870
      ok &= CGOSpecial(tmpCGO, DOTSIZE_WITH_SPHERESCALE);
      convertcgo = CGOOptimizeSpheresToVBONonIndexedNoShader(I->shaderCGO,
          CGO_BOUNDING_BOX_SZ + fsizeof<cgo::draw::sphere_buffers>() + 2);
      ok &= CGOAppendNoStop(tmpCGO, convertcgo); 
      CGOFreeWithoutVBOs(convertcgo);
      ok &= CGODisable(tmpCGO, GL_SPHERE_SHADER);
      CGOStop(tmpCGO);
      convertcgo = tmpCGO;
    } else {
      convertcgo = CGOOptimizeToVBONotIndexedNoShader(I->shaderCGO);
    }
    CHECKOK(ok, convertcgo);
    if (ok)
      convertcgo->use_shader = true;
  } else if (I->Type == 2) {
    CGO *tmpCGO = CGONew(G);
    CHECKOK(ok, tmpCGO);
    auto convertcgo2 = CGOConvertLinesToShaderCylinders(I->shaderCGO, 0);
    CHECKOK(ok, convertcgo2);
    CGOEnable(tmpCGO, GL_CYLINDER_SHADER);
    CGOSpecial(tmpCGO, MESH_WIDTH_FOR_SURFACES);
    convertcgo = CGOConvertShaderCylindersToCylinderShader(convertcgo2, tmpCGO);
    CGOAppendNoStop(tmpCGO, convertcgo);
    CGOFreeWithoutVBOs(convertcgo);
    CGODisable(tmpCGO, GL_CYLINDER_SHADER);
    CGOStop(tmpCGO);
    convertcgo = tmpCGO;
    convertcgo->use_shader = true;
    CGOFree(convertcgo2);
  } else {
    if((alpha != 1.0) || va) { // semi-transparent
      if (ok)
	convertcgo = CGOOptimizeToVBOIndexedWithColorEmbedTransparentInfo(I->shaderCGO, 0, 0, 0);
      CHECKOK(ok, convertcgo);
#ifdef _PYMOL_IOS
#endif
    } else {
      if (ok)
	convertcgo = CGOOptimizeToVBONotIndexedWithReturnedData(I->shaderCGO, 0, false, NULL);
      CHECKOK(ok, convertcgo);
    }
    if (ok)
      convertcgo->use_shader = true;
    {
      CGO *tmpCGO = NULL;
      tmpCGO = CGONew(G);
      CGOEnable(tmpCGO, GL_SURFACE_SHADER);
      //CGOEnable(tmpCGO, CGO_GL_LIGHTING); // do we need this?
      CGOSpecial(tmpCGO, SET_SURFACE_UNIFORMS);
      CGOAppendNoStop(tmpCGO, convertcgo);
      CGODisable(tmpCGO, GL_SURFACE_SHADER);
      CGOStop(tmpCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = tmpCGO;
      convertcgo->use_shader = true;
    }

  }

  if(ok && I->Type == 1 && !dot_as_spheres) {
    setShaderCGO(I, CGONew(G));
    CHECKOK(ok, I->shaderCGO);
    if (ok){
      I->shaderCGO->use_shader = true;
      if (ok) ok &= CGOResetNormal(I->shaderCGO, true);
      if (ok) ok &= CGOEnable(I->shaderCGO, GL_SURFACE_SHADER);
      if (ok) ok &= CGOSpecial(I->shaderCGO, POINTSIZE_DYNAMIC_DOT_WIDTH);
      if (ok) CGOAppendNoStop(I->shaderCGO, convertcgo);
      if (ok) ok &= CGODisable(I->shaderCGO, GL_SURFACE_SHADER);
      if (ok) ok &= CGOStop(I->shaderCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = NULL;
    }
  } else {
    setShaderCGO(I, convertcgo);
    convertcgo = NULL;
  }

  if (I->shaderCGO) {
    I->shaderCGO->no_pick = !pick_surface;
    setPickingCGO(I, I->shaderCGO);
  }

  return ok;
}

void RepSurface::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  float *v = I->V;
  float *vn = I->VN;
  float *vc = I->VC;
  float *va = I->VA;
  int *rc = I->RC;
  int *t = I->T;
  int *s = I->S;
  int c = I->N;
  int *vi = I->Vis;
  int ok = true;
  float alpha;
  int t_mode;
  float ambient_occlusion_scale = 0.f;
  int ambient_occlusion_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ambient_occlusion_mode);
  int ambient_occlusion_mode_div_4 = 0;
  if (ambient_occlusion_mode){
    ambient_occlusion_scale = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ambient_occlusion_scale);
    ambient_occlusion_mode_div_4 = ambient_occlusion_mode / 4;
  }

  if((I->Type != 1) && (!s)) {
    return;
  }

  alpha = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  if(ray) {
#ifndef _PYMOL_NO_RAY
    ray->transparentf(1.0F - alpha);
    if(I->Type == 1) {
      /* dot surface */
      float radius;
      radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_radius);
      if(radius == 0.0F) {
        radius = ray->PixelRadius * SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                                 cSetting_dot_width) / 1.4142F;
      }

      if(I->oneColorFlag) {
        float col[3];
        ColorGetEncoded(G, I->oneColor, col);
        ray->color3fv(col);
      }

      if(c)
        while(ok && c--) {
          if(*vi) {
            if(!I->oneColorFlag) {
              ray->color3fv(vc);
            }
            ok &= ray->sphere3fv(v, radius);
          }
          vi++;
          vc += 3;
          v += 3;
        }
    } else if((I->Type == 0) || (I->Type == 3) || (I->Type == 4) || (I->Type == 5)) {   /* solid surface */
      c = I->NT;

      if(I->oneColorFlag) {
        float col[3], col1[3], col2[3], col3[3];
        ColorGetEncoded(G, I->oneColor, col);
        while(ok && c--) {
          if (visibility_test(I->proximity, vi, t)) {
	    copy3f(col, col1);
	    copy3f(col, col2);
	    copy3f(col, col3);
	    if (I->VAO){
	      float ao1, ao2, ao3;
	      switch (ambient_occlusion_mode_div_4){
	      case 1:
		ao1 = 1.f-ambient_occlusion_scale*(*(I->VAO + *t));
		ao2 = 1.f-ambient_occlusion_scale*(*(I->VAO + *(t+1)));
		ao3 = 1.f-ambient_occlusion_scale*(*(I->VAO + *(t+2)));
		break;
	      case 2:
		ao1 = cos(.5f * PI * CLAMP_VALUE(ambient_occlusion_scale*(*(I->VAO + *t))));
		ao2 = cos(.5f * PI * CLAMP_VALUE(ambient_occlusion_scale*(*(I->VAO + *(t+1)))));
		ao3 = cos(.5f * PI * CLAMP_VALUE(ambient_occlusion_scale*(*(I->VAO + *(t+2)))));
		break;
	      default:
		ao1 = CLAMP_VALUE(1.f / (1.f + exp(.5f*((ambient_occlusion_scale*(*(I->VAO + *t))) - 10.f))));
		ao2 = CLAMP_VALUE(1.f / (1.f + exp(.5f*((ambient_occlusion_scale*(*(I->VAO + *(t+1)))) - 10.f))));
		ao3 = CLAMP_VALUE(1.f / (1.f + exp(.5f*((ambient_occlusion_scale*(*(I->VAO + *(t+2)))) - 10.f))));
	      }
	      mult3f(col1, ao1, col1);
	      mult3f(col2, ao2, col2);
	      mult3f(col3, ao3, col3);
	    }
	    ok &= ray->triangle3fv(v + (*t) * 3, v + (*(t + 1)) * 3, v + (*(t + 2)) * 3,
				    vn + (*t) * 3, vn + (*(t + 1)) * 3, vn + (*(t + 2)) * 3,
				    col1, col2, col3);
	  }
          t += 3;
        }
      } else {
        while(ok && c--) {
          int ttA = *t, ttB = *(t + 1), ttC = *(t + 2);
          if (visibility_test(I->proximity, vi, t)) {
            int ttA3 = ttA * 3, ttB3 = ttB * 3, ttC3 = ttC * 3;
            float cA[3], cB[3], cC[3];
	    copy3f(vc + ttA3, cA);
	    copy3f(vc + ttB3, cB);
	    copy3f(vc + ttC3, cC);
	    //            register float *cA = vc + ttA3, *cB = vc + ttB3, *cC = vc + ttC3;
            if(rc) {
              if(rc[ttA] < -1)
                ColorGetEncoded(G, rc[ttA], cA);
	      //                ColorGetEncoded(G, rc[ttA], (cA = colA));
              if(rc[ttB] < -1)
                ColorGetEncoded(G, rc[ttB], cB);
	      //                ColorGetEncoded(G, rc[ttB], (cB = colB));
              if(rc[ttC] < -1)
                ColorGetEncoded(G, rc[ttC], cC);
	      //                ColorGetEncoded(G, rc[ttC], (cC = colC));
            }
            if((*(vi + ttA)) || (*(vi + ttB)) || (*(vi + ttC))) {
	      if (I->VAO){
		float ao1, ao2, ao3;
		switch (ambient_occlusion_mode_div_4){
		case 1:
		  ao1 = 1.f-ambient_occlusion_scale*(*(I->VAO + *t));
		  ao2 = 1.f-ambient_occlusion_scale*(*(I->VAO + *(t+1)));
		  ao3 = 1.f-ambient_occlusion_scale*(*(I->VAO + *(t+2)));
		  break;
		case 2:
		  ao1 = cos(.5f * PI * CLAMP_VALUE(ambient_occlusion_scale*(*(I->VAO + *t))));
		  ao2 = cos(.5f * PI * CLAMP_VALUE(ambient_occlusion_scale*(*(I->VAO + *(t+1)))));
		  ao3 = cos(.5f * PI * CLAMP_VALUE(ambient_occlusion_scale*(*(I->VAO + *(t+2)))));
		  break;
		default:
		  ao1 = CLAMP_VALUE(1.f / (1.f + exp(.5f*((ambient_occlusion_scale*(*(I->VAO + *t))) - 10.f))));
		  ao2 = CLAMP_VALUE(1.f / (1.f + exp(.5f*((ambient_occlusion_scale*(*(I->VAO + *(t+1)))) - 10.f))));
		  ao3 = CLAMP_VALUE(1.f / (1.f + exp(.5f*((ambient_occlusion_scale*(*(I->VAO + *(t+2)))) - 10.f))));
		}
		mult3f(cA, ao1, cA);
		mult3f(cB, ao2, cB);
		mult3f(cC, ao3, cC);
	      }
              if(va) {
                ok &= ray->triangleTrans3fv(v + ttA3, v + ttB3, v + ttC3,
					     vn + ttA3, vn + ttB3, vn + ttC3,
					     cA, cB, cC,
					     1.0F - va[ttA], 1.0F - va[ttB], 1.0F - va[ttC]);
              } else {
                ok &= ray->triangle3fv(v + ttA3, v + ttB3, v + ttC3,
					vn + ttA3, vn + ttB3, vn + ttC3, cA, cB, cC);
              }
            }
          }
          t += 3;
        }
      }
    } else if(I->Type == 2) {   /* triangle mesh surface */

      float radius;
      int t0, t1, t2;
      int spacing = 10;
      int *cache = pymol::calloc<int>(spacing * (I->N + 1));
      CHECKOK(ok, cache);

      radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_radius);

      if(ok && radius == 0.0F) {
        float line_width =
          SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_width);
        line_width = SceneGetDynamicLineWidth(info, line_width);

        radius = ray->PixelRadius * line_width / 2.0F;
      }
      if (ok){
        float col[3], *col0, *col1, *col2;
	c = I->NT;
	if(I->oneColorFlag) {
	  ColorGetEncoded(G, I->oneColor, col);
          col0 = col1 = col2 = col;
        }
        while(ok && c--) {
	    t0 = (*t);
	    t1 = (*(t + 1));
	    t2 = (*(t + 2));
            if (visibility_test(I->proximity, vi, t)) {
              if (!I->oneColorFlag) {
                col0 = vc + t0 * 3;
                col1 = vc + t1 * 3;
                col2 = vc + t2 * 3;
              }
	      if(!check_and_add(cache, spacing, t0, t1))
		ok &= ray->sausage3fv(v + t0 * 3, v + t1 * 3, radius, col0, col1);
	      if(!check_and_add(cache, spacing, t1, t2))
		ok &= ray->sausage3fv(v + t1 * 3, v + t2 * 3, radius, col1, col2);
	      if(!check_and_add(cache, spacing, t2, t0))
		ok &= ray->sausage3fv(v + t2 * 3, v + t0 * 3, radius, col2, col0);
	    }
	    t += 3;
	}
	FreeP(cache);
      }
    }
    if (ok){
      ray->transparentf(0.0);
    } else {
      /* If not ok, then Clear Entire RepSurface, not just the ray object */
      
    }
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    /* Not ray tracing, but rendering */
    if(pick) {
      // Don't render transparent unpickable surfaces. Do render solid but
      // unpickable surfaces to write the depth buffer and prevent
      // through-picking.
      if (I->pickingCGO && !(I->pickingCGO->no_pick && alpha < 1.F)) {
	CGORenderGLPicking(I->pickingCGO, info, &I->context, cs->Setting.get(), obj->Setting.get());
      }
    } else {

#ifndef PURE_OPENGL_ES_2
      bool use_shader = SettingGetGlobal_b(G, cSetting_surface_use_shader) &&
        SettingGetGlobal_b(G, cSetting_use_shaders);

      if (use_shader && !info->alpha_cgo)
#endif
      {
        bool dot_as_spheres = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_as_spheres);

        if (!I->shaderCGO || CGOCheckWhetherToFree(G, I->shaderCGO) ||
            (I->Type == 1 && I->dot_as_spheres != dot_as_spheres)) {
          ok &= RepSurfaceCGOGenerate(I, info);
        }

        if (ok && I->shaderCGO) {
          const float *color = ColorGet(G, obj->Color);
          CGORenderGL(I->shaderCGO, color, NULL, NULL, info, I);
          return;
        }

        PRINTFB(G, FB_RepSurface, FB_Errors)
          " RepSurfaceCGOGenerate failed\n" ENDFB(G);
      }

      setShaderCGO(I, NULL);

#ifndef PURE_OPENGL_ES_2
      bool two_sided_lighting =
        SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
            cSetting_two_sided_lighting) > 0;
      if (two_sided_lighting){
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
      }

      if(I->Type == 1) {
        /* no triangle information, so we're rendering dots only */
        {
          int normals =
            SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_normals);
          int lighting = info->line_lighting ||
            SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_lighting);

          if(!normals){
            SceneResetNormal(G, true);
            vn = NULL;
          }

          if(!lighting)
            glDisable(GL_LIGHTING);

          glPointSize(SettingGet_f
              (G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_width));
          if(c) {
            glBegin(GL_POINTS);
            if(I->oneColorFlag) {
              glColor3fv(ColorGet(G, I->oneColor));
              vc = NULL;
            } else {
              glColor3f(1.0, 0.0, 0.0);
            }

            immediate_draw_masked_vertices(vc, vn, v, vi, c);
            glEnd();
          }

          if(!lighting)
            glEnable(GL_LIGHTING);
        } /* else use shader */
      } else if(I->Type == 2) { /* rendering triangle mesh */
        if (ok) {

          c = I->NT;
          if(ok && c) {
              int normals =
                SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_normals);
              int lighting = info->line_lighting ||
                SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_lighting);

              if(!normals){
                SceneResetNormal(G, true);
                vn = NULL;
              }

              if(!lighting)
                glDisable(GL_LIGHTING);

              float line_width =
                SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_mesh_width);
              glLineWidth(SceneGetDynamicLineWidth(info, line_width));

              if(I->oneColorFlag) {
                glColor3fv(ColorGet(G, I->oneColor));
                vc = NULL;
              }

              glBegin(GL_LINES);
              for (; c--; t += 3) {
                if (visibility_test(I->proximity, vi, t)) {
                  int indices[] = {t[0], t[1], t[1], t[2], t[2], t[0]};
                  immediate_draw_indexed_vertices(vc, vn, v, indices, 6);
                }
              }
              glEnd();

              if(!lighting)
                glEnable(GL_LIGHTING);
          }
	}
      } else {
        /* we're rendering triangles */

        if((alpha != 1.0) || va) {

          t_mode =
            SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
                         cSetting_transparency_mode);

          if(info && info->alpha_cgo) {
            t_mode = 0;
          }

          if(t_mode) {

            float **t_buf = NULL, **tb;
            float *z_value = NULL, *zv;
            int *ix = NULL;
            int n_tri = 0;
            float sum[3];
            float matrix[16];

            glGetFloatv(GL_MODELVIEW_MATRIX, matrix);

            if(I->oneColorFlag) {
              t_buf = pymol::malloc<float *>(I->NT * 6);
            } else {
              t_buf = pymol::malloc<float *>(I->NT * 12);
            }
	    CHECKOK(ok, t_buf);
	    if (ok){
	      z_value = pymol::malloc<float>(I->NT);
	      CHECKOK(ok, z_value);
	    }
	    if (ok){
	      ix = pymol::malloc<int>(I->NT);
	      CHECKOK(ok, ix);
	    }
            zv = z_value;
            tb = t_buf;
            c = I->NT;
	    if (ok){
	      if(I->oneColorFlag) {
		while(c--) {
                  if (visibility_test(I->proximity, vi, t)) {
		    *(tb++) = vn + (*t) * 3;
		    *(tb++) = v + (*t) * 3;
		    *(tb++) = vn + (*(t + 1)) * 3;
		    *(tb++) = v + (*(t + 1)) * 3;
		    *(tb++) = vn + (*(t + 2)) * 3;
		    *(tb++) = v + (*(t + 2)) * 3;
		    
		    add3f(tb[-1], tb[-3], sum);
		    add3f(sum, tb[-5], sum);
		    
		    *(zv++) = matrix[2] * sum[0] + matrix[6] * sum[1] + matrix[10] * sum[2];
		    n_tri++;
		  }
		  t += 3;
		}
	      } else {
		while(c--) {
                  if (visibility_test(I->proximity, vi, t)) {
		      if(va)
			*(tb++) = va + (*t);
		      else
			*(tb++) = &alpha;
		      *(tb++) = vc + (*t) * 3;
		      *(tb++) = vn + (*t) * 3;
		      *(tb++) = v + (*t) * 3;
		      
		      if(va)
			*(tb++) = va + (*(t + 1));
		      else
			*(tb++) = &alpha;
		      *(tb++) = vc + (*(t + 1)) * 3;
		      *(tb++) = vn + (*(t + 1)) * 3;
		      *(tb++) = v + (*(t + 1)) * 3;
		      
		      if(va)
			*(tb++) = va + (*(t + 2));
		      else
			*(tb++) = &alpha;
		      *(tb++) = vc + (*(t + 2)) * 3;
		      *(tb++) = vn + (*(t + 2)) * 3;
		      *(tb++) = v + (*(t + 2)) * 3;
		      
		      add3f(tb[-1], tb[-5], sum);
		      add3f(sum, tb[-9], sum);

		      *(zv++) =
			matrix[2] * sum[0] + matrix[6] * sum[1] + matrix[10] * sum[2];
		      n_tri++;
		    }
		  t += 3;
		}
	      }
	    }
            switch (t_mode) {
            case 1:
              ok &= UtilSemiSortFloatIndex(n_tri, z_value, ix, true);
              /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn); */
              break;
            default:
              ok &= UtilSemiSortFloatIndex(n_tri, z_value, ix, false);
              /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
              break;
            }
            c = n_tri;
            if(I->oneColorFlag) {
              float col[3];
              ColorGetEncoded(G, I->oneColor, col);
	      if (ok){
		glColor4f(col[0], col[1], col[2], alpha);
		glBegin(GL_TRIANGLES);
		for(c = 0; c < n_tri; c++) {
		  tb = t_buf + 6 * ix[c];
		  glNormal3fv(*(tb++));
		  glVertex3fv(*(tb++));
		  glNormal3fv(*(tb++));
		  glVertex3fv(*(tb++));
		  glNormal3fv(*(tb++));
		  glVertex3fv(*(tb++));
		}
		glEnd();
	      }
            } else { /* else I->oneColorFlag */
		glBegin(GL_TRIANGLES);
		for(c = 0; c < n_tri; c++) {
		  float *vv, *v_alpha;
		  
		  tb = t_buf + 12 * ix[c];
		  
		  v_alpha = *(tb++);
		  vv = *(tb++);
		  glColor4f(vv[0], vv[1], vv[2], *v_alpha);
		  glNormal3fv(*(tb++));
		  glVertex3fv(*(tb++));
		  
		  v_alpha = *(tb++);
		  vv = *(tb++);
		  glColor4f(vv[0], vv[1], vv[2], *v_alpha);
		  glNormal3fv(*(tb++));
		  glVertex3fv(*(tb++));
		  
		  v_alpha = *(tb++);
		  vv = *(tb++);
		  glColor4f(vv[0], vv[1], vv[2], *v_alpha);
		  glNormal3fv(*(tb++));
		  glVertex3fv(*(tb++));
		}
		glEnd();
	      }
	    {
	      FreeP(ix);
	      FreeP(z_value);
	      FreeP(t_buf);
	    }
          } else if (ok) {
            if(info->alpha_cgo) {       /* global transparency sort */
              if(I->allVisibleFlag) {
                if(I->oneColorFlag) {
                  float col[3];
                  ColorGetEncoded(G, I->oneColor, col);

                  glColor4f(col[0], col[1], col[2], alpha);
                  c = *(s++);
                  while(c) {
                    int parity = 0;
                    s += 2;
                    while(ok && c--) {
                      ok &= CGOAlphaTriangle(info->alpha_cgo,
					     v + s[-2] * 3, v + s[-1] * 3, v + (*s) * 3,
					     vn + s[-2] * 3, vn + s[-1] * 3, vn + (*s) * 3,
					     col, col, col, alpha, alpha, alpha, parity);
                      s++;
                      parity = !parity;
                    }
                    c = *(s++);
                  }
                } else {
                  c = *(s++);
                  while(ok && c) {
                    float *col0, *col1, *col2;
                    int parity = 0;
                    float alpha0, alpha1, alpha2;
                    col0 = vc + (*s) * 3;
                    if(va) {
                      alpha0 = va[(*s)];
                    } else {
                      alpha0 = alpha;
                    }
                    s++;
                    col1 = vc + (*s) * 3;
                    if(va) {
                      alpha1 = va[(*s)];
                    } else {
                      alpha1 = alpha;
                    }
                    s++;
                    while(ok && c--) {
                      col2 = vc + (*s) * 3;
                      if(va) {
                        alpha2 = va[(*s)];
                      } else {
                        alpha2 = alpha;
                      }
                      ok &= CGOAlphaTriangle(info->alpha_cgo,
					     v + s[-2] * 3, v + s[-1] * 3, v + (*s) * 3,
					     vn + s[-2] * 3, vn + s[-1] * 3, vn + (*s) * 3,
					     col0, col1, col2, alpha0, alpha1, alpha2, parity);
                      alpha0 = alpha1;
                      alpha1 = alpha2;
                      col0 = col1;
                      col1 = col2;
                      s++;
                      parity = !parity;
                    }
                    c = *(s++);
                  }
                }
              } else if (ok){          /* subset s */
                c = I->NT;
                if(c) {
                  if(I->oneColorFlag) {
                    float color[3];
                    ColorGetEncoded(G, I->oneColor, color);
                    while(ok && c--) {
                      if (visibility_test(I->proximity, vi, t))
                      {

                        ok &= CGOAlphaTriangle(info->alpha_cgo,
					       v + t[0] * 3, v + t[1] * 3, v + t[2] * 3,
					       vn + t[0] * 3, vn + t[1] * 3, vn + t[2] * 3,
					       color, color, color, alpha, alpha, alpha, 0);
                      }
                      t += 3;
                    }
                  } else {
                    while(ok && c--) {
                      if (visibility_test(I->proximity, vi, t))
                      {

                        if(va) {
                          ok &= CGOAlphaTriangle(info->alpha_cgo,
						 v + t[0] * 3, v + t[1] * 3, v + t[2] * 3,
						 vn + t[0] * 3, vn + t[1] * 3, vn + t[2] * 3,
						 vc + t[0] * 3, vc + t[1] * 3, vc + t[2] * 3,
						 va[t[0]], va[t[1]], va[t[2]], 0);
                        } else {
                          ok &= CGOAlphaTriangle(info->alpha_cgo,
						 v + t[0] * 3, v + t[1] * 3, v + t[2] * 3,
						 vn + t[0] * 3, vn + t[1] * 3, vn + t[2] * 3,
						 vc + t[0] * 3, vc + t[1] * 3, vc + t[2] * 3,
						 alpha, alpha, alpha, 0);
                        }
                      }
                      t += 3;
                    }
                  }
		  if (ok)
		    CGOEnd(info->alpha_cgo);
                }
              }
            } else if (ok){

              /* fast and ugly */
                if(I->oneColorFlag) {
                  float col[3];
                  ColorGetEncoded(G, I->oneColor, col);
                  glColor4f(col[0], col[1], col[2], alpha);
                  vc = NULL;
                }

                if(I->allVisibleFlag) {
                  for (; (c = *(s++)); s += c + 2) {
                    glBegin(GL_TRIANGLE_STRIP);
                    immediate_draw_indexed_vertices_alpha(vc, va, alpha, vn, v, s, c + 2);
                    glEnd();
                  }
                } else {
                  c = I->NT;
                  if(c) {
                    glBegin(GL_TRIANGLES);
                    for (; c--; t += 3) {
                      if (visibility_test(I->proximity, vi, t)) {
                        immediate_draw_indexed_vertices_alpha(vc, va, alpha, vn, v, t, 3);
                      }
                    }
                    glEnd();
                  }
                }
            }
          }
        } else if (ok) {                /* opaque */
            if(I->oneColorFlag) {
              glColor3fv(ColorGet(G, I->oneColor));
              vc = NULL;
            }

            if(I->allVisibleFlag) {
              for (; (c = *(s++)); s += c + 2) {
                glBegin(GL_TRIANGLE_STRIP);
                immediate_draw_indexed_vertices(vc, vn, v, s, c + 2);
                glEnd();
              }
            } else {
              c = I->NT;
              if(c) {
                glBegin(GL_TRIANGLES);
                for (; c--; t += 3) {
                  if (visibility_test(I->proximity, vi, t)) {
                    immediate_draw_indexed_vertices(vc, vn, v, t, 3);
                  }
                }
                glEnd();
              }
            }
        }
      }
      if(ok && SettingGetGlobal_i(G, cSetting_surface_debug)) {
        t = I->T;
        c = I->NT;

        glLineWidth(1.0F);

        // draw green triangles, either filled or as line edges,
        // depending on surface_type
        if(c) {
          glColor3f(0.0, 1.0, 0.0);
          if (I->Type == 2) {
            // filled green triangles for surface as mesh
            glBegin(GL_TRIANGLES);
            for (; c--; t += 3) {
              if(I->allVisibleFlag
                  || visibility_test(I->proximity, vi, t)) {
                immediate_draw_indexed_vertices(NULL, vn, v, t, 3);
              }
            }
            glEnd();
          } else {
            // line edges as green lines
            glBegin(GL_LINES);
            for (; c--; t += 3) {
              if(I->allVisibleFlag
                  || visibility_test(I->proximity, vi, t)) {
                int indices[] = {t[0], t[1], t[1], t[2], t[2], t[0]};
                immediate_draw_indexed_vertices(NULL, vn, v, indices, 6);
              }
            }
            glEnd();
          }
        }

        // draw normals as red lines
        c = I->N;
        if(c) {
          SceneResetNormal(G, true);
          glColor3f(1.0, 0.0, 0.0);
          glBegin(GL_LINES);
          while(c--) {
            glVertex3fv(v);
            glVertex3f(v[0] + vn[0] / 2, v[1] + vn[1] / 2, v[2] + vn[2] / 2);
            v += 3;
            vn += 3;
          }
          glEnd();
        }
      }

      if (two_sided_lighting){
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
      }
#endif
    }
  }
  if (!ok){
    I->invalidate(cRepInvPurge);
    I->cs->Active[cRepSurface] = false;
  }
}

bool RepSurface::sameVis() const
{
  for (int idx = 0; idx < cs->NIndex; ++idx) {
    auto* ai = cs->getAtomInfo(idx);
    if (LastVisib[idx] != GET_BIT(ai->visRep, cRepSurface)) {
      return false;
    }
  }

  return true;
}

void RepSurface::invalidate(cRepInv_t level)
{
  Rep::invalidate(level);

  if (level >= cRepInvColor)
    ColorInvalidated = true;
}

bool RepSurface::sameColor() const
{
  if (ColorInvalidated)
    return false;

  const auto* lc = LastColor;
  for (int idx = 0; idx < cs->NIndex; idx++) {
    auto* ai = cs->getAtomInfo(idx);
      if(ai->visRep & cRepSurfaceBit) {
        if(*(lc++) != ai->color) {
          return false;
        }
      }
    }

  return true;
}

Rep* RepSurface::recolor()
{
  auto const I = this;
  assert(cs == this->cs);

  MapType *map = NULL, *ambient_occlusion_map = NULL;
  int a, i0, i, j, c1;
  float *v0, *vc, *va;
  const float *c0;
  float *n0;
  int *lc;
  char *lv;
  int first_color;
  float *v_pos, v_above[3];
  int ramp_above;
  ObjectMolecule *obj;
  float probe_radius;
  float dist;
  float cutoff;
  int inclH;
  int inclInvis;
  int cullByFlag = false;
  int surface_mode, ambient_occlusion_mode;
  int surface_color;
  int *present = NULL;
  int *rc;
  int ramped_flag = false;

  int carve_state = 0;
  int carve_flag = false;
  float carve_cutoff;
  float carve_normal_cutoff;
  int carve_normal_flag;
  const char *carve_selection = NULL;
  float *carve_vla = NULL;
  MapType *carve_map = NULL;

  int clear_state = 0;
  int clear_flag = false;
  float clear_cutoff;
  const char *clear_selection = NULL;
  float *clear_vla = NULL;

  int state = getState();

  float transp;
  int variable_alpha = false;

  MapType *clear_map = NULL;

  AtomInfoType *ai2 = NULL, *ai1;

  obj = cs->Obj;
  ambient_occlusion_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ambient_occlusion_mode);
  surface_mode = I->surface_mode;
  ramp_above =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_ramp_above_mode);
  surface_color =
    SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_color);
  cullByFlag = (surface_mode == cRepSurface_by_flags);
  inclH = !((surface_mode == cRepSurface_heavy_atoms)
            || (surface_mode == cRepSurface_vis_heavy_only));
  inclInvis = !((surface_mode == cRepSurface_vis_only)
                || (surface_mode == cRepSurface_vis_heavy_only));
  probe_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_solvent_radius);
  I->proximity =
    SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_proximity);
  carve_cutoff =
    SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_carve_cutoff);
  clear_cutoff =
    SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_clear_cutoff);
  transp = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_transparency);
  carve_normal_cutoff =
    SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_carve_normal_cutoff);
  carve_normal_flag = carve_normal_cutoff > (-1.0F);

  cutoff = I->max_vdw + 2 * probe_radius;

  if(!I->LastVisib)
    I->LastVisib = pymol::malloc<char>(cs->NIndex);
  if(!I->LastColor)
    I->LastColor = pymol::malloc<int>(cs->NIndex);
  lv = I->LastVisib;
  lc = I->LastColor;
  for(a = 0; a < cs->NIndex; a++) {
    ai2 = cs->getAtomInfo(a);
    *(lv++) = GET_BIT(ai2->visRep, cRepSurface);
    *(lc++) = ai2->color;
  }

  if(I->N) {
    if(carve_cutoff > 0.0F) {
      carve_state =
        SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_carve_state) - 1;
      carve_selection =
        SettingGet_s(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_carve_selection);
      if(carve_selection)
        carve_map = SelectorGetSpacialMapFromSeleCoord(G,
                                                       SelectorIndexByName(G,
                                                                           carve_selection),
                                                       carve_state, carve_cutoff,
                                                       &carve_vla);
      if(carve_map)
        MapSetupExpress(carve_map);
      carve_flag = true;
    }

    if(clear_cutoff > 0.0F) {
      clear_state =
        SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_clear_state) - 1;
      clear_selection =
        SettingGet_s(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_clear_selection);
      if(clear_selection)
        clear_map = SelectorGetSpacialMapFromSeleCoord(G,
                                                       SelectorIndexByName(G,
                                                                           clear_selection),
                                                       clear_state, clear_cutoff,
                                                       &clear_vla);
      if(clear_map)
        MapSetupExpress(clear_map);
      clear_flag = true;
    }

    if(!I->VC)
      I->VC = pymol::malloc<float>(3 * I->N);
    vc = I->VC;
    if(!I->VA)
      I->VA = pymol::malloc<float>(I->N);
    va = I->VA;
    if(!I->RC)
      I->RC = pymol::malloc<int>(I->N);
    rc = I->RC;
    if(!I->Vis)
      I->Vis = pymol::malloc<int>(I->N);
    if(ColorCheckRamped(G, surface_color)) {
      I->oneColorFlag = false;
    } else {
      I->oneColorFlag = true;
    }
    first_color = -1;

    present = pymol::malloc<int>(cs->NIndex);
    {
      int *ap = present;
      for(a = 0; a < cs->NIndex; a++) {
        ai1 = obj->AtomInfo + cs->IdxToAtm[a];
        if((ai1->visRep & cRepSurfaceBit) &&
           (inclH || (!ai1->isHydrogen())) &&
           ((!cullByFlag) || (!(ai1->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
          *ap = 2;
        else
          *ap = 0;
        ap++;
      }
    }

    if(inclInvis) {
      float probe_radiusX2 = probe_radius * 2;
      map =
	MapNewFlagged(G, 2 * I->max_vdw + probe_radius, cs->Coord, cs->NIndex, NULL,
		      present);
      MapSetupExpress(map);
      /* add in nearby invisibles */
      for(a = 0; a < cs->NIndex; a++){
        if(!present[a]) {
          ai1 = obj->AtomInfo + cs->IdxToAtm[a];
          if((!cullByFlag) || !(ai1->flags & cAtomFlag_ignore)) {
            v0 = cs->coordPtr(a);
            i = *(MapLocusEStart(map, v0));
            if(i && map->EList) {
              j = map->EList[i++];
              while(j >= 0) {
                if(present[j] > 1) {
                  ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                  if(within3f
                     (cs->coordPtr(j), v0, ai1->vdw + ai2->vdw + probe_radiusX2)) {
                    present[a] = 1;
                    break;
                  }
                }
                j = map->EList[i++];
              }
            }
          }
        }
      }
      MapFree(map);
      map = NULL;
    }

    /**
     * Update the ambient occlusion accessibility array (VAO)
     */
    auto update_VAO = [&]() {
      if (!ambient_occlusion_mode) {
        VLAFreeP(I->VAO);
        return;
      }

      float maxDist = 0.f, maxDistA =0.f;
      int level_min = 64, level_max = 0;
      double start_time, cur_time;
      short vertex_map = 0; /* vertex or atom map */
      float map_cutoff = cutoff;

      switch (ambient_occlusion_mode % 4){
      case 1:
      case 3:
	vertex_map = 0; /* ambient_occlusion_mode - 1 atoms in map (default), 2 - vertices in map */
	break;
      case 2:
	vertex_map = 1;
      }
      
      if (!I->VAO){
	I->VAO = VLAlloc(float, I->N);
      } else {
	VLASize(I->VAO, float, I->N);
      }
      start_time = UtilGetSeconds(G);

      if (vertex_map){
	ambient_occlusion_map = MapNewFlagged(G, map_cutoff, I->V, I->N, NULL, NULL);
      } else {
	ambient_occlusion_map = MapNewFlagged(G, map_cutoff, cs->Coord, cs->NIndex, NULL, present);
      }      
      MapSetupExpress(ambient_occlusion_map);

      if (ambient_occlusion_mode==3){
	/* per atom */
	float *VAO = pymol::malloc<float>(cs->NIndex);
	short *nVAO = pymol::malloc<short>(cs->NIndex);
	memset(VAO, 0, sizeof(float)*cs->NIndex);
	memset(nVAO, 0, sizeof(short)*cs->NIndex);

	for(a = 0; a < I->N; a++) {
	  int nbits = 0;
	  short level1, level2, has;
	  unsigned long bits = 0L, bit;
	  float d[3], *vn0, v0mod[3];
	  int closeA = -1;
	  float closeDist = FLT_MAX;
	  has = 0;
	  
	  v0 = I->V + 3 * a;
	  vn0 = I->VN + 3 * a;
	  mult3f(vn0, .01f, v0mod);
	  add3f(v0, v0mod, v0mod);

	  i = *(MapLocusEStart(ambient_occlusion_map, v0));
	  if(i && map->EList) {
	    j = ambient_occlusion_map->EList[i++];
	    while(j >= 0) {
	      subtract3f(cs->coordPtr(j), v0, d);
	      dist = (float) length3f(d);
	      if (dist < closeDist){
		closeA = j;
		closeDist = dist;
	      }
	      j = ambient_occlusion_map->EList[i++];
	    }
	  }
	  if (closeA >= 0){
	    if (nVAO[closeA]){
	      I->VAO[a] = VAO[closeA];
	    } else {
	      v0 = cs->coordPtr(closeA);
	      i = *(MapLocusEStart(ambient_occlusion_map, v0));
	      if (i){
		j = ambient_occlusion_map->EList[i++];
		while(j >= 0) {
		  if (closeA==j){
		    j = ambient_occlusion_map->EList[i++];
		    continue;
		  }
		  subtract3f(cs->coordPtr(j), v0, d);
		  dist = (float) length3f(d);
		  if (dist > 12.f){
		    j = ambient_occlusion_map->EList[i++];
		    continue;
		  }
		  has = 1;
		  level1 = (d[2] < 0.f) ? 4 : 0;
		  level1 |= (d[1] < 0.f) ? 2 : 0;
		  level1 |= (d[0] < 0.f) ? 1 : 0;
		  d[0] = fabs(d[0]); d[1] = fabs(d[1]); d[2] = fabs(d[2]);
		  level2 = (d[0] <= d[1]) ? 4 : 0;
		  level2 |= (d[1] <= d[2]) ? 2 : 0;
		  level2 |= (d[0] <= d[2]) ? 1 : 0;
		  
		  bit = level1* 8 + level2;
		  bits |= (1L << bit);
		  j = ambient_occlusion_map->EList[i++];
		}
	      }
	      if (has){
		nbits = countBits(bits);
		if (nbits > level_max) level_max = nbits;
		if (nbits < level_min) level_min = nbits;
		I->VAO[a] = nbits;
	      } else {
		level_min = 0;
		I->VAO[a] = 0.f;
	      }
	      VAO[closeA] = I->VAO[a];
	      nVAO[closeA] = 1;
	    }
	  }
	}
	FreeP(VAO);
	FreeP(nVAO);
      } else {
	float natomsL = 0;
	for(a = 0; a < I->N; a++) {
	  int natoms = 0, nbits = 0;
	  short level1, level2, has;
	  unsigned long bits = 0L, bit;
	  float d[3], *vn0, pt[3], v0mod[3];
	  
	  if (a%1000==0){
	    PRINTFB(G, FB_RepSurface, FB_Debugging) "RepSurfaceColor():  Ambient Occlusion computing mode=%d #vertices=%d done=%d\n", ambient_occlusion_mode, I->N, a ENDFB(G);
	  }
	  v0 = I->V + 3 * a;
	  vn0 = I->VN + 3 * a;
	  mult3f(vn0, .01f, v0mod);
	  add3f(v0, v0mod, v0mod);
	  i = *(MapLocusEStart(ambient_occlusion_map, v0));
	  if(i && ambient_occlusion_map->EList) {
	    j = ambient_occlusion_map->EList[i++];
	    maxDistA = 0.f;
	    has = 0;
	    while(j >= 0) {
	      natomsL++;
	      if (vertex_map && a==j){
		j = ambient_occlusion_map->EList[i++];
		continue;
	      }
	      if (vertex_map){
		copy3f(I->V + j * 3, pt);
		subtract3f(I->V + j * 3, v0mod, d);
	      } else {
		copy3f(cs->coordPtr(j), pt);
		subtract3f(cs->coordPtr(j), v0mod, d);
	      }
	      dist = (float) length3f(d);
	      normalize3f(d);
	      if (get_angle3f(vn0, d) >= (PI/2.f)){
		j = ambient_occlusion_map->EList[i++];
		continue;
	      }
	      if (dist <= .0001f || dist > 12.f){
		j = ambient_occlusion_map->EList[i++];
		continue;
	      }
	      has = 1;
	      (dist > maxDistA) ? maxDistA = dist : 0;
	      level1 = (d[2] < 0.f) ? 4 : 0;
	      level1 |= (d[1] < 0.f) ? 2 : 0;
	      level1 |= (d[0] < 0.f) ? 1 : 0;
	      d[0] = fabs(d[0]); d[1] = fabs(d[1]); d[2] = fabs(d[2]);
	      level2 = (d[0] <= d[1]) ? 4 : 0;
	      level2 |= (d[1] <= d[2]) ? 2 : 0;
	      level2 |= (d[0] <= d[2]) ? 1 : 0;
	      
	      bit = level1* 8 + level2;
	      bits |= (1L << bit);
	      j = ambient_occlusion_map->EList[i++];
	      natoms++;
	    }
	    if(has){
	      nbits = countBits(bits);
	      if (nbits > level_max) level_max = nbits;
	      if (nbits < level_min) level_min = nbits;
	      (maxDistA > maxDist) ? maxDist = maxDistA : 0;
	      I->VAO[a] = nbits;
	    } else {
	      level_min = 0;
	      I->VAO[a] = 0.f;
	    }
	    maxDistA=0.f;
	  }
	}
	PRINTFB(G, FB_RepSurface, FB_Debugging) "RepSurfaceColor():  #vertices=%d Ambient Occlusion average #atoms looked at per vertex = %f\n", I->N, (natomsL/I->N) ENDFB(G);
      }
      {
	int ambient_occlusion_smooth = 
	  SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ambient_occlusion_smooth);
	int surface_quality = 
	  SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_quality);
	float min_max_diff = level_max - level_min;
	if (surface_quality>0)
	  ambient_occlusion_smooth *= surface_quality;
	/* Now we should set VAO from min/max */
	if (min_max_diff){
	  for(a = 0; a < I->N; a++) {
	    I->VAO[a] = ((I->VAO[a] - level_min)/min_max_diff);
	  }
	} else {
	  for(a = 0; a < I->N; a++) {
	    I->VAO[a] = 0.f;
	  }
	}

	/* SMOOTH Accessibility VALUES */
	if (ambient_occlusion_smooth && I->T){
	  int i, j, pt1, pt2, pt3;
	  float ave;
	  float *tmpVAO = pymol::malloc<float>(I->N);
	  int *nVAO = pymol::malloc<int>(I->N), c, *t;
	  
	  for (j=0; j<ambient_occlusion_smooth; j++){
	    memset(nVAO, 0, sizeof(int)*I->N);
	    memset(tmpVAO, 0, sizeof(float)*I->N);

	    t = I->T;
	    c = I->NT;
	    while (c--){
              if (I->allVisibleFlag ||
                  visibility_test(I->proximity, I->Vis, t)) {
		pt1 = *t; pt2 = *(t+1); pt3 = *(t+2);
		nVAO[pt1] += 1; nVAO[pt2] += 1; nVAO[pt3] += 1;

		ave = ave3(I->VAO[pt1], I->VAO[pt2], I->VAO[pt3]);

		tmpVAO[pt1] += ave;
		tmpVAO[pt2] += ave;
		tmpVAO[pt3] += ave;
	      }
	      t +=3;
	    }

	    for (i=0; i<I->N;i++){
	      if (nVAO[i]){
		/* only update if added and greater than original */
		//		if ((tmpVAO[i]/(float)nVAO[i]) > I->VAO[i])
		I->VAO[i] = tmpVAO[i]/(float)nVAO[i];
	      }
	    }
	  }
	  FreeP(tmpVAO);
	  FreeP(nVAO);
	}

      }
      MapFree(ambient_occlusion_map);
      cur_time = UtilGetSeconds(G);

      PRINTFB(G, FB_RepSurface, FB_Debugging) "RepSurfaceColor():  Ambient Occlusion computed #atoms=%d #vertices=%d time=%lf seconds\n", cs->NIndex, I->N, (cur_time-start_time) ENDFB(G);
      ambient_occlusion_map = NULL;

    }; // update_VAO

    /* now, assign colors to each point */
    map = MapNewFlagged(G, cutoff, cs->Coord, cs->NIndex, NULL, present);
    if(map) {
      short color_smoothing = SettingGetGlobal_i(G, cSetting_surface_color_smoothing);
      float color_smoothing_threshold = SettingGetGlobal_f(G, cSetting_surface_color_smoothing_threshold);
      int atm, ok = true;
      MapSetupExpress(map);
      ok &= !G->Interrupt;
      if (ok && !I->AT)
	I->AT = VLACalloc(int, I->N);
      for(a = 0; ok && a < I->N; a++) {
        float at_transp = transp;

        AtomInfoType *ai0 = NULL;
        float minDist = FLT_MAX, minDist2 = FLT_MAX, distDiff = FLT_MAX;
	int pi = -1, catm = -1; /* variables for color smoothing */
        AtomInfoType *pai = NULL, *pai2 = NULL; /* variables for color smoothing */
        c1 = 1;
        i0 = -1;
        v0 = I->V + 3 * a;
        n0 = I->VN + 3 * a;
        auto vi = I->Vis + a;
        /* colors */
        i = *(MapLocusEStart(map, v0));
        if(i && map->EList) {
          j = map->EList[i++];
          while(j >= 0) {
	    atm = cs->IdxToAtm[j];
            ai2 = obj->AtomInfo + atm;
            if((inclH || (!ai2->isHydrogen())) &&
               ((!cullByFlag) || (!(ai2->flags & cAtomFlag_ignore)))) {
              dist = (float) diff3f(v0, cs->coordPtr(j)) - ai2->vdw;
              if(color_smoothing){
		if (dist < minDist){
		  /* switching closest to 2nd closest */
		  pai2 = pai;
		  minDist2 = minDist;
		  pi = j;
		  catm = atm;
		  pai = ai2;
		  minDist = dist;
		} else if (dist < minDist2){
		  /* just setting second closest */
		  pai2 = ai2;
		  minDist2 = dist;
		}
	      } else if (dist < minDist) {
                i0 = j;
		catm = atm;
                ai0 = ai2;
                minDist = dist;
              }
            }
            j = map->EList[i++];
          }
        }
	I->AT[a] = catm;
        if (color_smoothing){
          i0 = pi;
          ai0 = pai;
          /* TODO: should find point closest to v0 between
                   atoms points (cs->coordPtr(pi)) and (cs->coordPtr(pi2)) (including vdw), then set this
                   distance to the distance between the vertex v0
                   and that point. We might want to use the normal
                   to compute this distance.
           */
          distDiff = fabs(minDist2-minDist);
        }
        if(i0 >= 0) {
          int at_surface_color = AtomSettingGetWD(G, ai0, cSetting_surface_color, surface_color);
          at_transp = AtomSettingGetWD(G, ai0, cSetting_transparency, transp);

          if(at_surface_color != -1) {
            c1 = at_surface_color;
            distDiff = FLT_MAX;
          } else {
            c1 = ai0->color;
          }

          if(I->oneColorFlag) {
            if(first_color >= 0) {
              if(first_color != c1)
                I->oneColorFlag = false;
            } else
              first_color = c1;
          }
          if(I->allVisibleFlag)
            *vi = 1;
          else {
            ai2 = obj->AtomInfo + cs->IdxToAtm[i0];
            if((ai2->visRep & cRepSurfaceBit) &&
               (inclH || (!ai2->isHydrogen())) &&
               ((!cullByFlag) || (!(ai2->flags &
                                    (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
              *vi = 1;
            else
              *vi = 0;
          }
        } else {
          *vi = 0;
        }
        if(carve_flag && (*vi)) {       /* is point visible, and are we carving? */
          *vi = 0;

          if(carve_map) {

            i = *(MapLocusEStart(carve_map, v0));
            if(i && carve_map->EList) {
              j = carve_map->EList[i++];
              while(j >= 0) {
                float *v_targ = carve_vla + 3 * j;
                if(within3f(v_targ, v0, carve_cutoff)) {
                  if(!carve_normal_flag) {
                    *vi = 1;
                    break;
                  } else {
                    float v_to[3];
                    subtract3f(v_targ, v0, v_to);
                    if(dot_product3f(v_to, n0) >= carve_normal_cutoff) {
                      *vi = 1;
                      break;
                    }
                  }
                }
                j = carve_map->EList[i++];
              }
            }
          }
        }
        if(clear_flag && (*vi)) {       /* is point visible, and are we clearing? */
          if(clear_map) {
            i = *(MapLocusEStart(clear_map, v0));
            if(i && clear_map->EList) {
              j = clear_map->EList[i++];
              while(j >= 0) {
                if(within3f(clear_vla + 3 * j, v0, clear_cutoff)) {
                  *vi = 0;
                  break;
                }
                j = clear_map->EList[i++];
              }
            }
          }
        }
        /*
           if(ColorCheckRamped(G,surface_color)) {
           c1 = surface_color;
           }
         */

        if(ColorCheckRamped(G, c1)) {
          I->oneColorFlag = false;
          switch (ramp_above) {
          case 1:
            copy3f(n0, v_above);
            scale3f(v_above, probe_radius, v_above);
            add3f(v0, v_above, v_above);
            v_pos = v_above;
            rc[0] = -1;
            break;
          default:
            v_pos = v0;
            rc[0] = c1;
            ramped_flag = true;
            break;
          }
          ColorGetRamped(G, c1, v_pos, vc, state);
          vc += 3;
          rc++;
        } else {
          if (color_smoothing && distDiff < color_smoothing_threshold && pai2) {
            const float *c2;
            float weight, weight2;
            if (color_smoothing==1){
              weight = 1.f + sin(.5f * PI * (distDiff / color_smoothing_threshold));
            } else {
              weight = 1.f + (distDiff / color_smoothing_threshold);
            }
            weight2 = 2.f - weight;
            c0 = ColorGet(G, c1);
            c2 = ColorGet(G, pai2->color);
            *(rc++) = c1;
            *(vc++) = ((weight*(*(c0++))) + (weight2*(*(c2++)))) / 2.f;
            *(vc++) = ((weight*(*(c0++))) + (weight2*(*(c2++)))) / 2.f;
            *(vc++) = ((weight*(*(c0++))) + (weight2*(*(c2++)))) / 2.f;
          } else {
            c0 = ColorGet(G, c1);
            *(rc++) = c1;
            *(vc++) = *(c0++);
            *(vc++) = *(c0++);
            *(vc++) = *(c0++);
          }
        }
        if(at_transp != transp)
          variable_alpha = true;
        *(va++) = 1.0F - at_transp;

        if (at_transp > 0) {
          I->setHasTransparency();
        }

        if(!*vi)
          I->allVisibleFlag = false;
      }
      MapFree(map);
    }
    if(variable_alpha)
      I->oneColorFlag = false;

    if(I->oneColorFlag) {
      I->oneColor = first_color;
    }

    // ambient occlusion
    update_VAO();
  }
  /*
     if(surface_color>=0) {
     I->oneColorFlag=true;
     I->oneColor=surface_color;
     }
   */

  if(G->HaveGUI) {
    setShaderCGO(I, NULL);
  }

  if(I->VA && (!variable_alpha)) {
    FreeP(I->VA);
    I->VA = NULL;
  }

  if(carve_map)
    MapFree(carve_map);
  VLAFreeP(carve_vla);
  if(clear_map)
    MapFree(clear_map);
  VLAFreeP(clear_vla);
  if((!ramped_flag)
     || (!SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_ray_color_ramps)))
    FreeP(I->RC);
  I->ColorInvalidated = false;
  FreeP(present);

  return this;
}

typedef struct {
  /* input */
  float *coord;
  SurfaceJobAtomInfo *atomInfo;

  float maxVdw;
  int allVisibleFlag;

  int nPresent;
  int *presentVla;

  int solventSphereIndex, sphereIndex;

  int surfaceType;
  int circumscribe;
  float probeRadius;
  float carveCutoff;
  float *carveVla;

  int surfaceMode;
  int surfaceSolvent;
  int cavityCull;
  float pointSep;
  float trimCutoff;
  float trimFactor;

  int cavityMode;
  float cavityRadius;
  float cavityCutoff;

  /* results */
  float *V, *VN;
  int N, *T, *S, NT;

} SurfaceJob;

static void SurfaceJobPurgeResult(PyMOLGlobals * G, SurfaceJob * I)
{
  I->N = 0;
  I->NT = 0;
  VLAFreeP(I->V);
  I->V = NULL;
  VLAFreeP(I->VN);
  I->VN = NULL;
  VLAFreeP(I->T);
  I->T = NULL;
  VLAFreeP(I->S);
  I->S = NULL;
}

#ifndef _PYMOL_NOPY
OV_INLINE PyObject *SurfaceJobAtomInfoAsPyTuple(SurfaceJobAtomInfo * atom_info)
{
  PyObject *result = NULL;
  if(atom_info) {
    ov_size size = 2 * VLAGetSize(atom_info) + 1;
    result = PyTuple_New(size);
    if(result) {
      ov_size i;
      PyTuple_SetItem(result, 0, PyInt_FromLong(2));    /* width of array */
      for(i = 1; i < size; i += 2) {
        PyTuple_SetItem(result, i, PyFloat_FromDouble(atom_info->vdw));
        PyTuple_SetItem(result, i + 1, PyInt_FromLong(atom_info->flags));
        atom_info++;
      }
    }
  }
  return (PConvAutoNone(result));
}

#if 0
OV_INLINE SurfaceJobAtomInfo *SurfaceJobAtomInfoVLAFromPyTuple(PyObject * tuple)
{
  SurfaceJobAtomInfo *result = NULL;
  if(tuple && PyTuple_Check(tuple)) {
    ov_size size = PyTuple_Size(tuple);
    if(size) {
      ov_size width = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
      if(width == 2) {
        ov_size vla_size = (size - 1) / 2;
        result = VLAlloc(SurfaceJobAtomInfo, vla_size);
        if(result) {
          SurfaceJobAtomInfo *atom_info = result;
          ov_size i;
          for(i = 1; i < size; i += 2) {
            atom_info->vdw = (float) PyFloat_AsDouble(PyTuple_GetItem(tuple, i));
            atom_info->flags = PyInt_AsLong(PyTuple_GetItem(tuple, i + 1));
            atom_info++;
          }
        }
      }
    }
  }
  return (result);
}
#endif

OV_INLINE PyObject *SurfaceJobInputAsTuple(PyMOLGlobals * G, SurfaceJob * I)
{
  PyObject *result = PyTuple_New(24);
  if(result) {
    PyTuple_SetItem(result, 0, PyString_FromString("SurfaceJob"));
    PyTuple_SetItem(result, 1, PyInt_FromLong(1));      /* version */
    PyTuple_SetItem(result, 2, PConvFloatVLAToPyTuple(I->coord));
    PyTuple_SetItem(result, 3, SurfaceJobAtomInfoAsPyTuple(I->atomInfo));
    PyTuple_SetItem(result, 4, PyFloat_FromDouble(I->maxVdw));
    PyTuple_SetItem(result, 5, PyInt_FromLong(I->allVisibleFlag));
    PyTuple_SetItem(result, 6, PyInt_FromLong(I->nPresent));
    PyTuple_SetItem(result, 7, PConvIntVLAToPyTuple(I->presentVla));
    PyTuple_SetItem(result, 8, PyInt_FromLong(I->solventSphereIndex));
    PyTuple_SetItem(result, 9, PyInt_FromLong(I->sphereIndex));
    PyTuple_SetItem(result, 10, PyInt_FromLong(I->surfaceType));
    PyTuple_SetItem(result, 11, PyInt_FromLong(I->circumscribe));
    PyTuple_SetItem(result, 12, PyFloat_FromDouble(I->probeRadius));
    PyTuple_SetItem(result, 13, PyFloat_FromDouble(I->carveCutoff));
    PyTuple_SetItem(result, 14, PConvFloatVLAToPyTuple(I->carveVla));
    PyTuple_SetItem(result, 15, PyInt_FromLong(I->surfaceMode));
    PyTuple_SetItem(result, 16, PyInt_FromLong(I->surfaceSolvent));
    PyTuple_SetItem(result, 17, PyInt_FromLong(I->cavityCull));
    PyTuple_SetItem(result, 18, PyFloat_FromDouble(I->pointSep));
    PyTuple_SetItem(result, 19, PyFloat_FromDouble(I->trimCutoff));
    PyTuple_SetItem(result, 20, PyFloat_FromDouble(I->trimFactor));

    PyTuple_SetItem(result, 21, PyInt_FromLong(I->cavityMode));
    PyTuple_SetItem(result, 22, PyFloat_FromDouble(I->cavityRadius));
    PyTuple_SetItem(result, 23, PyFloat_FromDouble(I->cavityCutoff));
    
  }
  return result;
}

OV_INLINE PyObject *SurfaceJobResultAsTuple(PyMOLGlobals * G, SurfaceJob * I)
{
  PyObject *result = PyTuple_New(6);
  if(result) {
    PyTuple_SetItem(result, 0, PyInt_FromLong(I->N));
    PyTuple_SetItem(result, 1, PConvFloatVLAToPyTuple(I->V));
    PyTuple_SetItem(result, 2, PConvFloatVLAToPyTuple(I->VN));
    PyTuple_SetItem(result, 3, PyInt_FromLong(I->NT));
    PyTuple_SetItem(result, 4, PConvIntVLAToPyTuple(I->T));
    PyTuple_SetItem(result, 5, PConvIntVLAToPyTuple(I->S));
  }
  return result;
}

OV_INLINE ov_status SurfaceJobResultFromTuple(PyMOLGlobals * G,
                                                     SurfaceJob * I, PyObject * tuple)
{
  ov_status status = OV_STATUS_FAILURE;
  SurfaceJobPurgeResult(G, I);
  if(tuple && PyTuple_Check(tuple)) {
    ov_size size = PyTuple_Size(tuple);
    if(size >= 6) {
      status = OV_STATUS_SUCCESS;

      I->N = PyInt_AsLong(PyTuple_GetItem(tuple, 0));
      if(OV_OK(status))
        status = PConvPyTupleToFloatVLA(&I->V, PyTuple_GetItem(tuple, 1));
      if(OV_OK(status))
        status = PConvPyTupleToFloatVLA(&I->VN, PyTuple_GetItem(tuple, 2));
      I->NT = PyInt_AsLong(PyTuple_GetItem(tuple, 3));
      if(OV_OK(status))
        status = PConvPyTupleToIntVLA(&I->T, PyTuple_GetItem(tuple, 4));
      if(OV_OK(status))
        status = PConvPyTupleToIntVLA(&I->S, PyTuple_GetItem(tuple, 5));
    }
    if(OV_ERR(status))
      SurfaceJobPurgeResult(G, I);
  }
  return status;
}
#endif

static SurfaceJob *SurfaceJobNew(PyMOLGlobals * G)
{
  OOCalloc(G, SurfaceJob);
  return I;
}

static void SurfaceJobFree(PyMOLGlobals * G, SurfaceJob * I)
{
  SurfaceJobPurgeResult(G, I);
  VLAFreeP(I->coord);
  VLAFreeP(I->presentVla);
  VLAFreeP(I->atomInfo);
  VLAFreeP(I->carveVla);
  OOFreeP(I);
}

static int SurfaceJobEliminateCloseDotsType3orMore(PyMOLGlobals * G,
    SurfaceJob * I, int *repeat_flag, int *dot_flag)
{
  int ok = true;
  int jj;
  float dist;
  float nearest;
  float point_sep = I->pointSep;
  float min_sep2 = point_sep * point_sep;
  float diff[3];
  {
    int a;
    for(a = 0; a < I->N; a++)
      dot_flag[a] = 1;
  }
  {
    MapType *map = MapNew(G, point_sep + 0.05F, I->V, I->N, NULL);
    int a;
    float *v = I->V;
    float *vn = I->VN;
    float min_dot = 0.1F;
    CHECKOK(ok, map);
    if (ok)
      ok &= MapSetupExpress(map);
    for(a = 0; ok && a < I->N; a++) {
      if(dot_flag[a]) {
	int i = *(MapLocusEStart(map, v));
	if(i && map->EList) {
	  int j = map->EList[i++];
	  jj = I->N;
	  nearest = point_sep + 1.0F;
	  while(j >= 0) {
	    if(j > a) {
	      if(dot_flag[j]) {
		if(dot_product3f(I->VN + (3 * j), vn) > min_dot) {
		  if(within3fret
		     (I->V + (3 * j), v, point_sep, min_sep2, diff, &dist)) {
		    *repeat_flag = true;
		    if(dist < nearest) {
		      /* try to be as determinstic as possible
			 in terms of how we collapse points */
		      jj = j;
		      nearest = dist;
		    } else if((j < jj) && (fabs(dist - nearest) < R_SMALL4)) {
		      jj = j;
		      nearest = dist;
		    }
		  }
		}
	      }
	    }
	    j = map->EList[i++];
	  }
	  if(jj < I->N) {
	    dot_flag[jj] = 0;
	    add3f(vn, I->VN + (3 * jj), vn);
	    average3f(I->V + (3 * jj), v, v);
	    *repeat_flag = true;
	  }
	}
      }
      v += 3;
      vn += 3;
      ok &= !G->Interrupt;
    }
    MapFree(map);
  }
  return ok;
}

static int SurfaceJobEliminateCloseDotsTypeLessThan3(PyMOLGlobals * G,
    SurfaceJob * I, int *repeat_flag, int *dot_flag)
{
  int ok = true;
  int a;
  float point_sep = I->pointSep;
  MapType *map = MapNew(G, -point_sep, I->V, I->N, NULL);
  float *v = I->V;
  float *vn = I->VN;

  CHECKOK(ok, map);
  if (ok){
    for(a = 0; a < I->N; a++)
      dot_flag[a] = 1;
    ok &= MapSetupExpress(map);
  }
  for(a = 0; ok && a < I->N; a++) {
    if(dot_flag[a]) {
      int i = *(MapLocusEStart(map, v));
      if(i && map->EList) {
	int j = map->EList[i++];
	while(j >= 0) {
	  if(j != a) {
	    if(dot_flag[j]) {
	      if(within3f(I->V + (3 * j), v, point_sep)) {
		dot_flag[j] = 0;
		add3f(vn, I->VN + (3 * j), vn);
		average3f(I->V + (3 * j), v, v);
		*repeat_flag = true;
	      }
	    }
	  }
	  j = map->EList[i++];
	}
      }
    }
    v += 3;
    vn += 3;
    ok &= !G->Interrupt;
  }
  MapFree(map);
  return ok;
}

static void SurfaceJobEliminateFromVArrays(PyMOLGlobals * G, SurfaceJob * I,
    int *dot_flag, short normalize)
{
  float *v0 = I->V;
  float *vn0 = I->VN;
  int *p = dot_flag;
  int c = I->N;
  int a;
  float *v = I->V;
  float *vn = I->VN;
  I->N = 0;
  for(a = 0; a < c; a++) {
    if(*(p++)) {
      *(v0++) = *(v++);
      *(v0++) = *(v++);
      *(v0++) = *(v++);
      if (normalize)
	normalize3f(vn);
      *(vn0++) = *(vn++);
      *(vn0++) = *(vn++);
      *(vn0++) = *(vn++);
      I->N++;
    } else {
      v += 3;
      vn += 3;
    }
  }
}

static int SurfaceJobEliminateCloseDots(PyMOLGlobals * G, SurfaceJob * I){
  int ok = true;
  if(I->N) {
    int repeat_flag = true;
    int *dot_flag = pymol::malloc<int>(I->N);
    CHECKOK(ok, dot_flag);
    while(ok && repeat_flag) {
      repeat_flag = false;
      if(I->surfaceType >= 3) {
	ok = SurfaceJobEliminateCloseDotsType3orMore(G, I, &repeat_flag, dot_flag);
      } else {            /* surface types < 3 */
	ok = SurfaceJobEliminateCloseDotsTypeLessThan3(G, I, &repeat_flag, dot_flag);
      }
      if(ok) {
	SurfaceJobEliminateFromVArrays(G, I, dot_flag, true);
      }
      ok &= !G->Interrupt;
    }
    FreeP(dot_flag);
  }
  return ok;
}

/* For each vertex, lookup all vertices within the neighborhood, and sum the dot_product of the normals.  
   If the average of the dot_products of the normals is less than the trim_cutoff,
   then the middle vertex is eliminated. */
static int SurfaceJobEliminateTroublesomeVerticesMark(PyMOLGlobals * G,
    SurfaceJob * I, int *repeat_flag, MapType *map, int *dot_flag,
    float neighborhood, float trim_cutoff)
{
  int ok = true;
  int a;
  float *v = I->V;
  float *vn = I->VN;
  for(a = 0; ok && a < I->N; a++) {
    if(dot_flag[a]) {
      int i = *(MapLocusEStart(map, v));
      if(i && map->EList) {
	int j = map->EList[i++];
	int n_nbr = 0;
	float dot_sum = 0.0F;
	while(j >= 0) {
	  if(j != a) {
	    if(dot_flag[j]) {
	      float *v0 = I->V + 3 * j;
	      if(within3f(v0, v, neighborhood)) {
		float *n0 = I->VN + 3 * j;
		dot_sum += dot_product3f(n0, vn);
		n_nbr++;
	      }
	    }
	  }
	  j = map->EList[i++];
	}
	if(n_nbr) {
	  dot_sum /= n_nbr;
	  if(dot_sum < trim_cutoff) {
	    dot_flag[a] = false;
	    *repeat_flag = true;
	  }
	}
      }
    }
    v += 3;
    vn += 3;
    ok &= !G->Interrupt;
  }
  return ok;
}

static int SurfaceJobEliminateTroublesomeVertices(PyMOLGlobals * G, SurfaceJob * I){
  int ok = true;
  if((I->surfaceType != 3) &&
     I->N && (I->trimCutoff > 0.0F) && (I->trimFactor > 0.0F)) {
    float trim_cutoff = I->trimCutoff;
    float trim_factor = I->trimFactor;
    int repeat_flag = true;
    float point_sep = I->pointSep;
    float neighborhood = trim_factor * point_sep;
    int *dot_flag = pymol::malloc<int>(I->N);
    CHECKOK(ok, dot_flag);
    if(ok && I->surfaceType == 6) {       /* emprical tweaks */
      trim_factor *= 2.5;
      trim_cutoff *= 1.5;
    }
    while(ok && repeat_flag) {
      MapType *map = MapNew(G, neighborhood, I->V, I->N, NULL);
      CHECKOK(ok, map);
      if (ok){
	int a;
	for(a = 0; a < I->N; a++)
	  dot_flag[a] = 1;
	ok &= MapSetupExpress(map);
      }
      repeat_flag = false;
      ok &= SurfaceJobEliminateTroublesomeVerticesMark(G, I, &repeat_flag, map, dot_flag, neighborhood, trim_cutoff);
      if(ok) {
	SurfaceJobEliminateFromVArrays(G, I, dot_flag, true);
      }
      MapFree(map);
      ok &= !G->Interrupt;
    }
    FreeP(dot_flag);
  }
  return ok;
}

static int SurfaceJobAtomProximityCleanupPass(PyMOLGlobals * G, SurfaceJob * I,
    int *dot_flag, int *present_vla, float probe_radius)
{
  int ok = true;
  float cutoff = 0.5 * probe_radius;
  float *I_coord = I->coord;
  SurfaceJobAtomInfo *I_atom_info = I->atomInfo;
  int n_index = VLAGetSize(I->atomInfo);
  MapType *map =
    MapNewFlagged(G, I->maxVdw + probe_radius, I_coord, n_index, NULL,
		  present_vla);
  int a;
  float *v;
  CHECKOK(ok, map);
  if (ok)
    ok &= MapSetupExpress(map);
  v = I->V;
  for(a = 0; ok && a < I->N; a++) {
    int i = *(MapLocusEStart(map, v));
    if(i && map->EList) {
      int j = map->EList[i++];
      while(j >= 0) {
	SurfaceJobAtomInfo *atom_info = I_atom_info + j;
	if((!present_vla) || present_vla[j]) {
	  if(within3f(I_coord + 3 * j, v, atom_info->vdw + cutoff)) {
	    dot_flag[a] = true;
	  }
	}
	j = map->EList[i++];
      }
    }
    v += 3;
    ok &= !G->Interrupt;
  }
  MapFree(map);
  return ok;
}

static int SurfaceJobRefineCopyNewPoints(SurfaceJob * I, float *new_dot, int n_new){
  int ok = true;
  float *n1 = new_dot + 3;
  float *v1 = new_dot;
  float *v, *vn;
  VLASize(I->V, float, 3 * (I->N + n_new));
  CHECKOK(ok, I->V);
  if (ok)
    VLASize(I->VN, float, 3 * (I->N + n_new));
  CHECKOK(ok, I->VN);
  if (ok){
    v = I->V + 3 * I->N;
    vn = I->VN + 3 * I->N;
    I->N += n_new;
  }
  while(ok && n_new--) {
    copy3f(v1, v);
    copy3f(n1, vn);
    v += 3;
    vn += 3;
    v1 += 6;
    n1 += 6;
  }
  return ok;
}
static int SurfaceJobRefineAddNewVerticesCheckPoint(SurfaceJob * I, MapType *map,
    int *n_new, float **new_dot, int j, float *v, float *vn, float map_cutoff,
    float neighborhood, float insert_cutoff)
{
  int ok = true;
  float *v0 = I->V + 3 * j;
  if(within3f(v0, v, map_cutoff)) {
    int add_new = false;
    float *n0 = I->VN + 3 * j;
    VLACheck(*new_dot, float, (*n_new) * 6 + 5);
    CHECKOK(ok, *new_dot);
    if (ok){
      float *v1 = (*new_dot) + (*n_new) * 6;
      average3f(v, v0, v1);
      if((dot_product3f(n0, vn) < 0.666 /* dot_cutoff, was hardcoded as a variable */ )
	 && (within3f(v0, v, neighborhood))){
	// if the normals are further than dot_cutoff apart
	// and the related points are close to each other, than add new
	add_new = true;
      } else {
	/* if points are too far apart, insert a new one, i.e., 
	   search for any point within insert_cutoff, if not, add */
	int ii = *(MapLocusEStart(map, v1));
	if(ii) {
	  int found = false;
	  int jj = map->EList[ii++];
	  while(jj >= 0) {
	    if(jj != j) {
	      float *vv0 = I->V + 3 * jj;
	      if(within3f(vv0, v1, insert_cutoff)) {
		found = true;
		break;
	      }
	    }
	    jj = map->EList[ii++];
	  }
	  if(!found)
	    add_new = true;
	}
      }
      if(add_new) {
	/* highly divergent, add dot in-between v and v0
	   (averaged, v1 set above, compute the normal (averaged below))
	   and n_new incremented */
	float *n1 = v1 + 3;
	(*n_new)++;
	average3f(vn, n0, n1);
	normalize3f(n1);
      }
    }
  }
  return ok;
}

static int SurfaceJobRefineAddNewVertices(PyMOLGlobals * G, SurfaceJob * I){
  int ok = true;
  float point_sep = I->pointSep;
  int n_new = 0;
  float neighborhood = 2.6 * point_sep; /* these constants need more tuning... */
  float insert_cutoff = 1.1 * point_sep;
  float map_cutoff = neighborhood;
  float *new_dot = VLAlloc(float, 1000);
  float *v, *vn;
  if(map_cutoff < (2.9 * point_sep)) {  /* these constants need more tuning... */
    map_cutoff = 2.9 * point_sep;
  }
  {
    MapType *map = NULL;
    int a;
    map = MapNew(G, map_cutoff, I->V, I->N, NULL);
    CHECKOK(ok, map);
    if (ok)
      ok &= MapSetupExpress(map);
    v = I->V;
    vn = I->VN;
    for(a = 0; ok && a < I->N; a++) {
      int i = *(MapLocusEStart(map, v));
      if(i && map->EList) {
	int j = map->EList[i++];
	while(ok && j >= 0) {
	  if(j > a) {
	    SurfaceJobRefineAddNewVerticesCheckPoint(I, map, &n_new, &new_dot, j, v, vn, map_cutoff, neighborhood, insert_cutoff);
	  }
	  j = map->EList[i++];
	  ok &= !G->Interrupt;
	}
      }
      v += 3;
      vn += 3;
      ok &= !G->Interrupt;
    }
    MapFree(map);
  }
  if(ok && n_new) {
    ok = SurfaceJobRefineCopyNewPoints(I, new_dot, n_new);
  }
  VLAFreeP(new_dot);
  return ok;
}

static void SurfaceJobCheckInteriorSolventSurface(MapType *solv_map, float *v,
    SolventDot *sol_dot, float probe_rad_less, float dist2, int a, int *flag){
  int ii;
  ii = *(MapLocusEStart(solv_map, v));
  if(ii && solv_map->EList) {
    float *i_dot = sol_dot->dot;
    float dist = probe_rad_less;
    int *elist_ii = solv_map->EList + ii;
    float v_0 = v[0];
    int jj_next, jj = *(elist_ii++);
    float v_1 = v[1];
    float *v1 = i_dot + 3 * jj;
    float v_2 = v[2];
    while(jj >= 0) {
      /* huge bottleneck -- optimized for superscaler processors */
      float dx = v1[0], dy, dz;
      jj_next = *(elist_ii++);
      dx -= v_0;
      if(jj != a) {
	dx = (dx < 0.0F) ? -dx : dx;
	dy = v1[1] - v_1;
	if(!(dx > dist)) {
	  dy = (dy < 0.0F) ? -dy : dy;
	  dz = v1[2] - v_2;
	  if(!(dy > dist)) {
	    dx = dx * dx;
	    dz = (dz < 0.0F) ? -dz : dz;
	    dy = dy * dy;
	    if(!(dz > dist)) {
	      dx = dx + dy;
	      dz = dz * dz;
	      if(!(dx > dist2))
		if((dx + dz) <= dist2) {
		  *flag = false;
		  break;
		}
	    }
	  }
	}
      }
      v1 = i_dot + 3 * jj_next;
      jj = jj_next;
    }
  }
}

static void SurfaceJobCheckPresentAndWithin(MapType *map, SurfaceJob * I,
    int *present_vla, float *v, float probe_rad_more, int *flag){
  SurfaceJobAtomInfo *I_atom_info = I->atomInfo;
  float *I_coord = I->coord;

  int i = *(MapLocusEStart(map, v));
  if(i && map->EList) {
    int j = map->EList[i++];
    while(j >= 0) {
      SurfaceJobAtomInfo *atom_info = I_atom_info + j;
      if((!present_vla) || present_vla[j]) {
	if(within3f
	   (I_coord + 3 * j, v,
	    atom_info->vdw + probe_rad_more)) {
	  *flag = false;
	  break;
	}
      }
      j = map->EList[i++];
    }
  }
}

static void SurfaceJobSetProbeRadius(int surface_type, float point_sep,
    float *probe_radius, float *probe_rad_more, float *probe_rad_less,
    float *probe_rad_less2)
{
  float solv_tole = point_sep * 0.04F;

  if(*probe_radius < (2.5F * point_sep)) { /* minimum probe radius allowed */
    *probe_radius = 2.5F * point_sep;
  }
  *probe_rad_more = *probe_radius * (1.0F + solv_tole);
  switch (surface_type) {
  case 0:                /* solid */
  case 3:
  case 4:
  case 5:
  case 6:
    *probe_rad_less = *probe_radius;
    break;
  default:
    *probe_rad_less = *probe_radius * (1.0F - solv_tole);
    break;
  }
  *probe_rad_less2 = (*probe_rad_less) * (*probe_rad_less);
}

static int SurfaceJobRun(PyMOLGlobals * G, SurfaceJob * I)
{
  int ok = true;
  int MaxN;
  int n_present = I->nPresent;
  SphereRec *sp = G->Sphere->Sphere[I->sphereIndex];
  SphereRec *ssp = G->Sphere->Sphere[I->solventSphereIndex];

  SurfaceJobPurgeResult(G, I);

  {
    /* compute limiting storage requirements */
    int tmp = n_present;
    if(tmp < 1)
      tmp = 1;
    if(sp->nDot < ssp->nDot)
      MaxN = tmp * ssp->nDot;
    else
      MaxN = tmp * sp->nDot;
  }
  MaxN = n_present;
  I->V = VLAlloc(float, (MaxN + 1) * 3);
  CHECKOK(ok, I->V);
  if (ok)
    I->VN = VLAlloc(float, (MaxN + 1) * 3);
  CHECKOK(ok, I->VN);

  ok &= !G->Interrupt;

  if(!(I->V && I->VN) || (!ok)) {       /* bail out point -- try to reduce crashes */
    PRINTFB(G, FB_RepSurface, FB_Errors)
      "Error-RepSurface: insufficient memory to calculate surface at this quality.\n"
      ENDFB(G);
    VLAFreeP(I->V);
    I->V = NULL;
    VLAFreeP(I->VN);
    I->VN = NULL;
    ok = false;
  } else {
    SolventDot *sol_dot = NULL;
    float *v = I->V;
    float *vn = I->VN;
    MapType *carve_map = NULL;
    float probe_radius = I->probeRadius;
    int circumscribe = I->circumscribe;
    int surface_type = I->surfaceType;
    float point_sep = I->pointSep;
    int *present_vla = I->presentVla;

    I->N = 0;

    sol_dot = SolventDotNew(G, I->coord, I->atomInfo, probe_radius,
                            ssp, present_vla,
                            circumscribe, I->surfaceMode, I->surfaceSolvent,
                            I->cavityCull, I->allVisibleFlag, I->maxVdw,
                            I->cavityMode, I->cavityRadius, I->cavityCutoff);
    CHECKOK(ok, sol_dot);
    ok &= !G->Interrupt;
    if(ok) {
      if(!I->surfaceSolvent) {
        float probe_rad_more, probe_rad_less, probe_rad_less2;

	SurfaceJobSetProbeRadius(surface_type, point_sep, &probe_radius, &probe_rad_more, &probe_rad_less, &probe_rad_less2);

        if(surface_type >= 5) { /* effectively double-weights atom points */
          if(sol_dot->nDot) {
            int a;
            float *v0 = sol_dot->dot;
            float *n0 = sol_dot->dotNormal;
            for(a = 0; a < sol_dot->nDot; a++) {
              scale3f(n0, -probe_radius, v);
              add3f(v0, v, v);
              copy3f(n0, vn);
              v += 3;
              vn += 3;
              n0 += 3;
              v0 += 3;
              I->N++;
            }
          }
        }
	ok &= !G->Interrupt;
        if(ok) {
          MapType *map, *solv_map = NULL;
          map = MapNewFlagged(G, I->maxVdw + probe_rad_more,
                              I->coord, VLAGetSize(I->coord) / 3, NULL, NULL);
	  CHECKOK(ok, map);
	  if (ok)
	    solv_map = MapNew(G, probe_rad_less, sol_dot->dot, sol_dot->nDot, NULL);
	  CHECKOK(ok, solv_map);
          if(ok) {
            ok &= MapSetupExpress(solv_map);
	    if (ok)
	      ok &= MapSetupExpress(map);
	    ok &= !G->Interrupt;
	    ok &= map->EList && solv_map->EList;
            if(sol_dot->nDot && ok) {
              Vector3f *dot = pymol::malloc<Vector3f>(sp->nDot);
              float *v0, *n0;
	      CHECKOK(ok, dot);
              if (ok){
                int b;
                for(b = 0; b < sp->nDot; b++) {
                  scale3f(sp->dot[b], probe_radius, dot[b]);
                }
              }
              v0 = sol_dot->dot;
              n0 = sol_dot->dotNormal;
              if (ok) {
                int a, b;
                int sp_nDot = sp->nDot;
                for(a = 0; ok && a < sol_dot->nDot; a++) {
                  if(sol_dot->dotCode[a] || (surface_type < 6)) {     /* surface type 6 is completely scribed */
                    OrthoBusyFast(G, a + sol_dot->nDot * 2, sol_dot->nDot * 5); /* 2/5 to 3/5 */
                    for(b = 0; ok && b < sp_nDot; b++) {
                      float *dot_b = dot[b];
                      v[0] = v0[0] + dot_b[0];
                      v[1] = v0[1] + dot_b[1];
                      v[2] = v0[2] + dot_b[2];
                      {
                        int flag = true;
			SurfaceJobCheckInteriorSolventSurface(solv_map, v, sol_dot, probe_rad_less, probe_rad_less2, a, &flag);
                        /* at this point, we have points on the interior of the solvent surface,
                           so now we need to further trim that surface to cover atoms that are present */
                        if(flag) {
			  SurfaceJobCheckPresentAndWithin(map, I, present_vla, v, probe_rad_more, &flag);
                          if(!flag) {   /* compute the normals */
                            vn[0] = -sp->dot[b][0];
                            vn[1] = -sp->dot[b][1];
                            vn[2] = -sp->dot[b][2];
                              I->N++;
			    VLACheck(I->V, float, 3 * (I->N + 1));
			    VLACheck(I->VN, float, 3 * (I->N + 1));
			      CHECKOK(ok, I->V);
			      CHECKOK(ok, I->VN);
			    v = I->V + I->N * 3;
			    vn = I->VN + I->N * 3;
                          }
                        }
                      }
		      ok &= !G->Interrupt;
                    }
                  }
                  v0 += 3;
                  n0 += 3;
		  ok &= !G->Interrupt;
                }
              }
              FreeP(dot);
            }
          }
          MapFree(solv_map);
          MapFree(map);
        }
      } else {
        float *v0 = sol_dot->dot;
        float *n0 = sol_dot->dotNormal;
        int a;
        circumscribe = 0;
        if(sol_dot->nDot) {
	  VLACheck(I->V, float, 3 * (I->N + sol_dot->nDot));
	  VLACheck(I->VN, float, 3 * (I->N + sol_dot->nDot));
	  v = I->V;
	  vn = I->VN;
          for(a = 0; a < sol_dot->nDot; a++) {
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            *(v++) = *(v0++);
            *(vn++) = *(n0++);
            I->N++;
          }
        }
      }
    }
    SolventDotFree(sol_dot);
    sol_dot = NULL;
    ok &= !G->Interrupt;
    if(ok) {
      int refine, ref_count = 1;

      if((surface_type == 0) && (circumscribe)) {
        ref_count = 2;          /* these constants need more tuning... */
      }
      for(refine = 0; ok && refine < ref_count; refine++) {
        /* add new vertices in regions where curvature is very high
           or where there are gaps with no points */
        if(I->N && (surface_type == 0) && (circumscribe)) {
	  ok = SurfaceJobRefineAddNewVertices(G, I);
        }

        if(ok && I->N && (surface_type == 0) && (circumscribe)) {
          /* combine scribing with an atom proximity cleanup pass */
          int *dot_flag = pymol::calloc<int>(I->N);
	  CHECKOK(ok, dot_flag);
	  ok &= SurfaceJobAtomProximityCleanupPass(G, I, dot_flag, present_vla, probe_radius);
	  /* purge unused dots */
	  if (ok)
	    SurfaceJobEliminateFromVArrays(G, I, dot_flag, false); // not normalize?
          FreeP(dot_flag);
        }

        /* now, eliminate dots that are too close to each other */
	ok &= SurfaceJobEliminateCloseDots(G, I);

        /* now eliminate troublesome vertices in regions of extremely high curvature */
	ok &= SurfaceJobEliminateTroublesomeVertices(G, I);
	ok &= !G->Interrupt;
      }
    }

    if(ok && I->N && I->V && I->VN) {
      VLASizeForSure(I->V, float, 3 * I->N);
      CHECKOK(ok, I->V);
      if (ok)
	VLASizeForSure(I->VN, float, 3 * I->N);
      CHECKOK(ok, I->VN);
    }

    PRINTFB(G, FB_RepSurface, FB_Blather)
      " RepSurface: %i surface points.\n", I->N ENDFB(G);

    ok &= !G->Interrupt;

    OrthoBusyFast(G, 3, 5);
    if(I->N) {
      if(ok && surface_type != 1) {   /* not a dot surface... */
        float cutoff = point_sep * 5.0F;
        if((cutoff > probe_radius) && (!I->surfaceSolvent))
          cutoff = probe_radius;
        I->T = TrianglePointsToSurface(G, I->V, I->VN, I->N, cutoff, &I->NT, &I->S, NULL, 
                                       I->cavityMode);
	CHECKOK(ok, I->T);
        PRINTFB(G, FB_RepSurface, FB_Blather)
          " RepSurface: %i triangles.\n", I->NT ENDFB(G);
      }
    } else {
      if (ok)
	VLASizeForSure(I->V, float, 1);
      CHECKOK(ok, I->V);
      if (ok)
	VLASizeForSure(I->VN, float, 1);
      CHECKOK(ok, I->VN);
    }
    if(carve_map)
      MapFree(carve_map);
  }
  return ok;
}

static void RepSurfaceSetSettings(PyMOLGlobals * G, CoordSet * cs,
    ObjectMolecule *obj, int surface_quality, int surface_type, float *point_sep,
    int *sphere_idx, int *solv_sph_idx, int *circumscribe)
{
  if(surface_quality >= 4) {        /* totally impractical */
    *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_best) / 4.f;
    *sphere_idx = 4;
    *solv_sph_idx = 4;
    if(*circumscribe < 0)
      *circumscribe = 91;
  } else {
    switch (surface_quality) {
    case 3:                /* nearly impractical */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_best) / 3.f;
      *sphere_idx = 4;
      *solv_sph_idx = 3;
      if(*circumscribe < 0)
	*circumscribe = 71;
      break;
    case 2:
      /* nearly perfect */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_best) / 2.f;
      *sphere_idx = 3;
      *solv_sph_idx = 3;
      if(*circumscribe < 0)
	*circumscribe = 41;
      break;
    case 1:
      /* good */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_best);
      *sphere_idx = 2;
      *solv_sph_idx = 3;
      if((*circumscribe < 0) && (surface_type == 6))
	*circumscribe = 40;
      break;
    case 0:
      /* 0 - normal */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_normal);
      *sphere_idx = 1;
      *solv_sph_idx = 2;
      if((*circumscribe < 0) && (surface_type == 6))
	*circumscribe = 30;
      break;
    case -1:
      /* -1 */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_poor);
      *sphere_idx = 1;
      *solv_sph_idx = 2;
      if((*circumscribe < 0) && (surface_type == 6))
	*circumscribe = 10;
      break;
    case -2:
      /* -2 god awful */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_poor) * 1.5F;
      *sphere_idx = 1;
      *solv_sph_idx = 1;
      break;
    case -3:
      /* -3 miserable */
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_miserable);
      *sphere_idx = 1;
      *solv_sph_idx = 1;
      break;
    default:
      *point_sep = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_miserable) * 1.18F;
      *sphere_idx = 0;
      *solv_sph_idx = 1;
    }
  }
  /* Fixed problem with surface holes when surface_quality>2, it seems like circumscribe can only be
     used with surface_solvent */
  if((*circumscribe < 0) || (!SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_solvent)))
    *circumscribe = 0;
}

static int RepSurfacePrepareSurfaceJob(PyMOLGlobals * G, SurfaceJob *surf_job,
    RepSurface *I, CoordSet *cs, ObjectMolecule *obj, SurfaceJobAtomInfo **atom_info,
    float *carve_vla, int n_present, int *present_vla, int optimize, int sphere_idx,
    int solv_sph_idx, int surface_type, int circumscribe, float probe_radius,
    float point_sep, float carve_cutoff)
{
  int ok = true;
  surf_job->maxVdw = I->max_vdw;
  surf_job->allVisibleFlag = I->allVisibleFlag;
  
  surf_job->atomInfo = *atom_info;
  (*atom_info) = NULL;
  
  surf_job->nPresent = n_present;
  if(present_vla && optimize) {
    /* implies that n_present < cs->NIndex, so eliminate
       irrelevant atoms & coordinates if we are optimizing subsets */
    surf_job->coord = VLAlloc(float, n_present * 3);
    CHECKOK(ok, surf_job->coord);
    if (ok) {
      int *p = present_vla;
      SurfaceJobAtomInfo *ai_src = surf_job->atomInfo;
      SurfaceJobAtomInfo *ai_dst = surf_job->atomInfo;
      const float *v_src = cs->Coord.data();
      float *v_dst = surf_job->coord;
      int a;
      
      for(a = 0; a < cs->NIndex; a++) {
	if(*(p++)) {
	  copy3f(v_src, v_dst);
	  v_dst += 3;
	  if(ai_dst != ai_src)
	    *ai_dst = *ai_src;
	  ai_dst++;
	}
	v_src += 3;
	ai_src++;
      }
    }
    VLASize(surf_job->atomInfo, SurfaceJobAtomInfo, n_present);
    CHECKOK(ok, surf_job->atomInfo);
  } else {
    surf_job->presentVla = present_vla;
    present_vla = NULL;
    surf_job->coord = VLAlloc(float, cs->NIndex * 3);
    CHECKOK(ok, surf_job->coord);
    if(ok)
      UtilCopyMem(surf_job->coord, cs->Coord, sizeof(float) * 3 * cs->NIndex);
  }
  if (ok){
    surf_job->sphereIndex = sphere_idx;
    surf_job->solventSphereIndex = solv_sph_idx;
    
    surf_job->surfaceType = surface_type;
    surf_job->circumscribe = circumscribe;
    surf_job->probeRadius = probe_radius;
    surf_job->pointSep = point_sep;
    
    surf_job->trimCutoff = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_trim_cutoff);
    surf_job->trimFactor = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_trim_factor);
    
    surf_job->cavityMode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_cavity_mode);
    surf_job->cavityRadius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_cavity_radius);
    surf_job->cavityCutoff = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_cavity_cutoff);
    if(carve_vla)
      surf_job->carveVla = VLACopy(carve_vla, float);
    surf_job->carveCutoff = carve_cutoff;
    surf_job->surfaceMode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_mode);
    surf_job->surfaceSolvent = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_solvent);
    surf_job->cavityCull = SettingGet_i(G, cs->Setting.get(),
					obj->Setting.get(), cSetting_cavity_cull);
  }
  return ok;
}

#ifndef _PYMOL_NOPY
static
void RepSurfaceConvertSurfaceJobToPyObject(PyMOLGlobals *G, SurfaceJob *surf_job, CoordSet *cs, ObjectMolecule *obj, PyObject **entry, PyObject **input, PyObject **output, int *found){
  int cache_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_cache_mode);
  
  if(cache_mode > 0) {
    int blocked = PAutoBlock(G);
    *input = SurfaceJobInputAsTuple(G, surf_job);
    
    if(PCacheGet(G, output, entry, *input) == OV_STATUS_YES) {
      if(OV_OK(SurfaceJobResultFromTuple(G, surf_job, *output))) {
	*found = true;
	PXDecRef(*input);
	*input = NULL;
	PXDecRef(*entry);
	*entry = NULL;
      }
      PXDecRef(*output);
      *output = NULL;
    }
    if(PyErr_Occurred())
      PyErr_Print();
    PAutoUnblock(G, blocked);
  }
}
#endif

static void RepSurfaceFindAllPresentAtoms(ObjectMolecule *obj, CoordSet *cs, int *present_vla, int inclH, int cullByFlag){
  int *ap = present_vla;
  const int *idx_to_atm = cs->IdxToAtm.data();
  const AtomInfoType *obj_AtomInfo = obj->AtomInfo.data();
  int a, cs_NIndex = cs->NIndex;
  for(a = 0; a < cs_NIndex; a++) {
    const AtomInfoType *ai1 = obj_AtomInfo + *(idx_to_atm++);
    if((ai1->visRep & cRepSurfaceBit) &&
       (inclH || (!ai1->isHydrogen())) &&
       ((!cullByFlag) || (!(ai1->flags &
			    (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
      *ap = 2;
    else
      *ap = 0;
    ap++;
  }
}

static int RepSurfaceAddNearByAtomsIfNotSurfaced(PyMOLGlobals *G, MapType *map,
    ObjectMolecule *obj, CoordSet *cs, int *present_vla, int inclH,
    int cullByFlag, float probe_radius, int optimize)
{
  int ok = true;
  float probe_radiusX2 = probe_radius * 2;
  int a;
  for(a = 0; ok && a < cs->NIndex; a++){
    if(!present_vla[a]) {
      AtomInfoType *ai1 = obj->AtomInfo + cs->IdxToAtm[a];
      if((inclH || (!ai1->isHydrogen())) &&
	 ((!cullByFlag) || 
	  !(ai1->flags & cAtomFlag_ignore))) {
        const float* v0 = cs->coordPtr(a);
        int i = *(MapLocusEStart(map, v0));
	if(optimize) {
	  if(i && map->EList) {
	    int j = map->EList[i++];
	    while(j >= 0) {
	      if(present_vla[j] > 1) {
		AtomInfoType *ai2 = obj->AtomInfo + cs->IdxToAtm[j];
		if(within3f
		   (cs->coordPtr(j), v0,
		    ai1->vdw + ai2->vdw + probe_radiusX2)) {
		  present_vla[a] = 1;
		  break;
		}
	      }
	      j = map->EList[i++];
	    }
	  }
	} else
	  present_vla[a] = 1;
      }
    }
    ok &= !G->Interrupt;
  }
  return ok;
}

static int RepSurfaceRemoveAtomsNotWithinCutoff(PyMOLGlobals *G,
    MapType *carve_map, float *carve_vla, CoordSet *cs, int *present_vla,
    float carve_cutoff)
{
  int ok = true;
  int a;
  for(a = 0; ok && a < cs->NIndex; a++) {
    int include_flag = false;
    if(carve_map) {
      const float* v0 = cs->coordPtr(a);
      int i = *(MapLocusEStart(carve_map, v0));
      if(i && carve_map->EList) {
	int j = carve_map->EList[i++];
	while(j >= 0) {
	  if(within3f(carve_vla + 3 * j, v0, carve_cutoff)) {
	    include_flag = true;
	    break;
	  }
	  j = carve_map->EList[i++];
	}
      }
    }
    if(!include_flag)
      present_vla[a] = 0;
    ok &= !G->Interrupt;
  }
  return ok;
}

Rep *RepSurfaceNew(CoordSet * cs, int state)
{
  int ok = true;
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj = cs->Obj;

    int surface_mode =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_mode);
    int cullByFlag = (surface_mode == cRepSurface_by_flags);
    int inclH = !((surface_mode == cRepSurface_heavy_atoms)
                  || (surface_mode == cRepSurface_vis_heavy_only));
    int inclInvis = !((surface_mode == cRepSurface_vis_only)
                      || (surface_mode == cRepSurface_vis_heavy_only));
    int visFlag = false;

    if(GET_BIT(obj->RepVisCache,cRepSurface)) {
      const int *idx_to_atm = cs->IdxToAtm.data();
      const AtomInfoType *obj_AtomInfo = obj->AtomInfo.data();
      int a, cs_NIndex = cs->NIndex;

      // If all atoms have "flag ignore", then don't cull by flags
      if (cullByFlag) {
        bool all_ignore = true;

        for (int idx = 0; idx < cs_NIndex; ++idx) {
          if (!(obj_AtomInfo[idx_to_atm[idx]].flags & cAtomFlag_ignore)) {
            all_ignore = false;
            break;
          }
        }

        if (all_ignore) {
          cullByFlag = false;
          surface_mode = cRepSurface_all;
        }
      }

      for(a = 0; a < cs_NIndex; a++) {
        const AtomInfoType *ai1 = obj_AtomInfo + *(idx_to_atm++);
        if((ai1->visRep & cRepSurfaceBit) &&
           (inclH || (!ai1->isHydrogen())) &&
           ((!cullByFlag) || (!(ai1->flags & (cAtomFlag_exfoliate | cAtomFlag_ignore))))) {
          visFlag = true;
          break;
        }
      }
    }
    if(!visFlag) {
      return (NULL);            /* skip if no thing visible */
    }

    auto I = new RepSurface(cs, state);
    I->surface_mode = surface_mode;

    {
      int surface_flag = false;
      int surface_type = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_type);        
      int surface_quality = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_quality);
      float probe_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_solvent_radius);
      int optimize = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_optimize_subsets);
      int circumscribe = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_circumscribe);
      int sphere_idx = 0, solv_sph_idx = 0;
      MapType *map = NULL;
      float point_sep;
      int *present_vla = NULL;
      int n_present = 0;
      int carve_state = 0;
      int carve_flag = false;
      float carve_cutoff;
      const char *carve_selection = NULL;
      float *carve_vla = NULL;
      MapType *carve_map = NULL;
      bool smooth_edges = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_smooth_edges);

      I->Type = surface_type;

      I->max_vdw = ObjectMoleculeGetMaxVDW(obj);

      RepSurfaceSetSettings(G, cs, obj, surface_quality, surface_type, &point_sep, &sphere_idx, &solv_sph_idx, &circumscribe);

      I->allVisibleFlag = true;
      obj = cs->Obj;

      /* don't waist time computing a Surface unless we need it!! */
      {
        const int *idx_to_atm = cs->IdxToAtm.data();
        const AtomInfoType *obj_AtomInfo = obj->AtomInfo.data();
        int a, cs_NIndex = cs->NIndex;
        for(a = 0; a < cs_NIndex; a++) {
          const AtomInfoType *ai1 = obj_AtomInfo + *(idx_to_atm++);
          if((ai1->visRep & cRepSurfaceBit) &&
             (inclH || (!ai1->isHydrogen())) &&
             ((!cullByFlag) || (!(ai1->flags &
                                  (cAtomFlag_exfoliate | cAtomFlag_ignore)))))
            surface_flag = true;
          else {
            I->allVisibleFlag = false;
            if(surface_flag)
              break;
          }
        }
      }
      if(surface_flag) {
        SurfaceJobAtomInfo *atom_info = VLACalloc(SurfaceJobAtomInfo, cs->NIndex);
	CHECKOK(ok, atom_info);
        if(ok && atom_info) {
          const AtomInfoType *i_ai, *obj_atom_info = obj->AtomInfo.data();
          const int *idx_to_atm = cs->IdxToAtm.data();
          int n_index = cs->NIndex;
          SurfaceJobAtomInfo *i_atom_info = atom_info;
          int i;
	  /* fill in surfacing flags into SurfaceJobAtomInfo array */
          for(i = 0; i < n_index; i++) {
            i_ai = obj_atom_info + idx_to_atm[i];
            i_atom_info->flags = i_ai->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate);
            i_atom_info->vdw = i_ai->vdw;
            i_atom_info++;
          }
        }

        OrthoBusyFast(G, 0, 1);

	if (ok){
	  n_present = cs->NIndex;
	  carve_selection =
	    SettingGet_s(G, cs->Setting.get(), obj->Setting.get(),
			 cSetting_surface_carve_selection);
	  carve_cutoff =
	    SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_surface_carve_cutoff);
	  if((!carve_selection) || (!carve_selection[0]))
	    carve_cutoff = 0.0F;
	  if(carve_cutoff > 0.0F) {
	    carve_state =
	      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
			   cSetting_surface_carve_state) - 1;
	    carve_cutoff += 2 * I->max_vdw + probe_radius;
	    
	    if(carve_selection)
	      carve_map = SelectorGetSpacialMapFromSeleCoord(G,
							     SelectorIndexByName(G,
										 carve_selection),
							     carve_state, carve_cutoff,
							     &carve_vla);
	    if(carve_map)
	      ok &= MapSetupExpress(carve_map);
	    carve_flag = true;
	    I->allVisibleFlag = false;
	  }
	}
        if(ok && !I->allVisibleFlag) {
          /* optimize the space over which we calculate a surface */
          /* first find out which atoms are actually to be surfaced */
          present_vla = VLAlloc(int, cs->NIndex);
	  CHECKOK(ok, present_vla);
          if (ok){
	    RepSurfaceFindAllPresentAtoms(obj, cs, present_vla, inclH, cullByFlag);
          }

	  if (ok)
	    map =
	      MapNewFlagged(G, 2 * I->max_vdw + probe_radius, cs->Coord, cs->NIndex, NULL,
			    present_vla);
	  CHECKOK(ok, map);
	  if (ok)
	    ok &= MapSetupExpress(map);

          if(ok && inclInvis) {
            /* then add in the nearby atoms which are not surfaced and not ignored */
	    ok &= RepSurfaceAddNearByAtomsIfNotSurfaced(G, map, obj, cs, present_vla, inclH, cullByFlag, probe_radius, optimize);
          }

          if(ok && carve_flag && (!optimize)) {
            /* and optimize for carved region */
	    ok &= RepSurfaceRemoveAtomsNotWithinCutoff(G, carve_map, carve_vla, cs, present_vla, carve_cutoff);
          }
          MapFree(map);
          map = NULL;

          /* now count how many atoms we actually need to think about */
          n_present = 0;
          if (ok) {
            int a;
            for(a = 0; a < cs->NIndex; a++) {
              if(present_vla[a]) {
                n_present++;
              }
            }
          }
        }

        if (ok) {
          SurfaceJob *surf_job = SurfaceJobNew(G);
	  CHECKOK(ok, surf_job);
	  ok &= RepSurfacePrepareSurfaceJob(G, surf_job, I, cs, obj, &atom_info, carve_vla, n_present, present_vla, optimize, sphere_idx, solv_sph_idx, surface_type, circumscribe, probe_radius, point_sep, carve_cutoff);

          ok &= !G->Interrupt;

          if(ok) {
            int found = false;
#ifndef _PYMOL_NOPY
            PyObject *entry = NULL;
            PyObject *output = NULL;
            PyObject *input = NULL;
	    int cache_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_cache_mode);
	    RepSurfaceConvertSurfaceJobToPyObject(G, surf_job, cs, obj, &entry, &input, &output, &found);
#endif
            if(ok && !found) {

              ok &= SurfaceJobRun(G, surf_job);

#ifndef _PYMOL_NOPY
              if(cache_mode > 1) {
                int blocked = PAutoBlock(G);
                output = SurfaceJobResultAsTuple(G, surf_job);
                PCacheSet(G, entry, output);
                PXDecRef(entry); entry = NULL;
                PXDecRef(output); output = NULL;
                PXDecRef(input); input = NULL;
                PAutoUnblock(G, blocked);
              }
#endif
            }
#ifndef _PYMOL_NOPY
            if(entry || input || output) {
              int blocked = PAutoBlock(G);
              PXDecRef(entry);
              PXDecRef(input);
              PXDecRef(output);
              PAutoUnblock(G, blocked);
            }
#endif
          }
          /* surf_job must be valid at this point */
	  if (ok){
	    I->N = surf_job->N;
	    surf_job->N = 0;
	    I->V = surf_job->V;
	    surf_job->V = NULL;
	    I->VN = surf_job->VN;
	    surf_job->VN = NULL;
	    I->NT = surf_job->NT;
	    surf_job->NT = 0;
	    I->T = surf_job->T;
	    surf_job->T = NULL;
	    I->S = surf_job->S;
	    surf_job->S = NULL;
	  }
          SurfaceJobPurgeResult(G, surf_job);
          SurfaceJobFree(G, surf_job);
        }
        VLAFreeP(atom_info);

	ok &= !G->Interrupt;
        if(ok)
          I->recolor();

        if (ok && smooth_edges)
          RepSurfaceSmoothEdges(I);
      }
      if(carve_map)
        MapFree(carve_map);
      VLAFreeP(carve_vla);
      VLAFreeP(present_vla);
      OrthoBusyFast(G, 4, 4);
    }

  if(!ok) {
    delete I;
    I = NULL;
  }
  return (Rep *) I;
}

void RepSurfaceSmoothEdges(RepSurface * I)
{
  if (I->allVisibleFlag)
    return;

  std::vector<std::vector<int>> edges(I->N);
  // Find all edges.
  // Uses the following edge structure:
  // edge[v0] -> { v1, t0, t1 }
  // where (v0, v1) defines an edge and (t0, t1)
  // are triangles the edge belongs to.
  // The edge is unique when t0 == t1.
  int *t = I->T;
  int *vi = I->Vis;
  for (auto i = 0; i < I->NT; i++) {
    if (visibility_test(I->proximity, vi, t)) {
      int v0 = t[0];
      int v1 = t[1];
      int v2 = t[2];
      if (v0 < v1) {
        edges[v0].push_back(v1);
        edges[v0].push_back(i);
        edges[v0].push_back(i);
      } else {
        edges[v1].push_back(v0);
        edges[v1].push_back(i);
        edges[v1].push_back(i);
      }
      if (v1 < v2) {
        edges[v1].push_back(v2);
        edges[v1].push_back(i);
        edges[v1].push_back(i);
      } else {
        edges[v2].push_back(v1);
        edges[v2].push_back(i);
        edges[v2].push_back(i);
      }
      if (v2 < v0) {
        edges[v2].push_back(v0);
        edges[v2].push_back(i);
        edges[v2].push_back(i);
      } else {
        edges[v0].push_back(v2);
        edges[v0].push_back(i);
        edges[v0].push_back(i);
      }
    }
    t += 3;
  }
  // Assign triangles to every edge
  for (auto i = 0; i < I->N; i++) {
    for (auto j = 0; j < edges[i].size(); j += 3) {
      for (auto k = 0; k < j; k += 3) {
        if ((j != k) && (edges[i][j] == edges[i][k])) {
          edges[i][j + 2] = edges[i][k + 1];
          edges[i][k + 2] = edges[i][j + 1];
        }
      }
    }
  }
  std::vector<int> unique_edges;
  // Build a list of unique edges
  for (auto i = 0; i < I->N; i++) {
    for (auto j = 0; j < edges[i].size(); j += 3) {
      if (edges[i][j + 1] == edges[i][j + 2]) {
        unique_edges.push_back(i);
        unique_edges.push_back(edges[i][j]);
      }
    }
  }
  if (unique_edges.size() > 0) {
    std::vector<float> new_vertices(3 * I->N);
    UtilCopyMem(&new_vertices[0], I->V, sizeof(float) * 3 * I->N);
    // Average positions of vertices shared by consecutive unique edges
    for (auto i = 0; i < unique_edges.size(); i += 2) {
      for (auto j = 0; j < i; j += 2) {
        auto v0 = -1, v1 = -1, v2 = -1;
        if (unique_edges[i] == unique_edges[j]) {
          v0 = unique_edges[i];
          v1 = unique_edges[i + 1];
          v2 = unique_edges[j + 1];
        } else if (unique_edges[i + 1] == unique_edges[j]) {
          v0 = unique_edges[i + 1];
          v1 = unique_edges[i];
          v2 = unique_edges[j + 1];
        } else if (unique_edges[i] == unique_edges[j + 1]) {
          v0 = unique_edges[i];
          v1 = unique_edges[i + 1];
          v2 = unique_edges[j];
        } else if (unique_edges[i + 1] == unique_edges[j + 1]) {
          v0 = unique_edges[i + 1];
          v1 = unique_edges[i];
          v2 = unique_edges[j];
        }
        if (v0 >= 0) {
          for (auto k = 0; k < 3; k++) {
            new_vertices[3 * v0 + k] =
              (1. / 3.) * (I->V[3 * v0 + k] +
                           I->V[3 * v1 + k] +
                           I->V[3 * v2 + k]);
          }
        }
      }
    }
    UtilCopyMem(I->V, &new_vertices[0], sizeof(float) * 3 * I->N);
  }
}

static int SolventDotFilterOutSameXYZ(PyMOLGlobals * G, MapType *map,
    SurfaceJobAtomInfo * atom_info, SurfaceJobAtomInfo *a_atom_info,
    float *coord, int a, int *present, int *skip_flag)
{
  int ok = true;
  float *v0 = coord + 3 * a;
  int i = *(MapLocusEStart(map, v0));
  if(i && map->EList) {
    int j = map->EList[i++];
    while(ok && j >= 0) {
      SurfaceJobAtomInfo *j_atom_info = atom_info + j;
      if(j > a)       /* only check if this is atom trails */
	if((!present) || present[j]) {
	  if(j_atom_info->vdw == a_atom_info->vdw) {        /* handle singularities */
	    float *v1 = coord + 3 * j;
	    if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2])){
	      // is this necessary? if some different atoms have exact same x/y/z
	      *skip_flag = true;
	    }
	  }
	}
      j = map->EList[i++];
      ok &= !G->Interrupt;
    }
  }
  return ok;
}

static int SolventDotGetDotsAroundVertexInSphere(PyMOLGlobals * G, SolventDot *I,
    MapType *map, SurfaceJobAtomInfo * atom_info, SurfaceJobAtomInfo *a_atom_info,
    float *coord, int a, int *present, SphereRec * sp, float radius, int *dotCnt,
    int stopDot, float *dotPtr, float *dotNormal, int *nDot)
{
  float vdw = a_atom_info->vdw + radius;
  float *v0 = coord + 3 * a;
  int b, ok = true;
  float *v = dotPtr + (*nDot) * 3;
  float *n = NULL;
  Vector3f *sp_dot = sp->dot;
  if (dotNormal){
    n = dotNormal + (*nDot) * 3;
  }
  for(b = 0; ok && b < sp->nDot; b++) {
    float *sp_dot_b = (float*)(sp_dot + b);
    int i;
    int flag = true;
    if (n){
      n[0] = sp_dot_b[0]; n[1] = sp_dot_b[1]; n[2] = sp_dot_b[2];
    }
    v[0] = v0[0] + vdw * sp_dot_b[0];
    v[1] = v0[1] + vdw * sp_dot_b[1];
    v[2] = v0[2] + vdw * sp_dot_b[2];
    i = *(MapLocusEStart(map, v));
    if(i) {
      int j = map->EList[i++];
      while(ok && j >= 0) {
	SurfaceJobAtomInfo *j_atom_info = atom_info + j;
	if((!present) || present[j]) {
	  if(j != a) {
	    int skip_flag = false;
	    if(j_atom_info->vdw == a_atom_info->vdw) {      /* handle singularities */
	      float *v1 = coord + 3 * j;
	      if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2])){
		skip_flag = true;
	      }
	    }
	    if(!skip_flag)
	      if(within3f(coord + 3 * j, v, j_atom_info->vdw + radius)) {
		flag = false;
		break;
	      }
	  }
	}
	j = map->EList[i++];
	ok &= !G->Interrupt;
      }
    }
    if(ok && flag && (*dotCnt < stopDot)) {
      (*dotCnt)++;
      v += 3;
      if (n)
	n += 3;
      (*nDot)++;
    }
  }
  return ok;
}

static int SolventDotCircumscribeAroundVertex(PyMOLGlobals * G, SolventDot *I,
    MapType *map, float *vdw, float dist, float *v0, float *v2, int circumscribe,
    SurfaceJobAtomInfo * atom_info, SurfaceJobAtomInfo *a_atom_info,
    SurfaceJobAtomInfo *jj_atom_info, int *present, int a, int jj, float *coord,
    float probe_radius, int *dotCnt, int stopDot)
{
  int ok = true;
  float vz[3], vx[3], vy[3], vp[3];
  float tri_a = vdw[0], tri_b = vdw[2], tri_c = dist;
  float tri_s = (tri_a + tri_b + tri_c) * 0.5F;
  float area = (float) sqrt1f(tri_s * (tri_s - tri_a) *
			      (tri_s - tri_b) * (tri_s - tri_c));
  float radius = (2 * area) / dist;
  float adj = (float) sqrt1f(vdw[1] - radius * radius);
  int b;
  float *v = I->dot + I->nDot * 3;
  float *n = I->dotNormal + I->nDot * 3;

  subtract3f(v2, v0, vz);
  get_system1f3f(vz, vx, vy);
  
  copy3f(vz, vp);
  scale3f(vp, adj, vp);
  add3f(v0, vp, vp);
  
  for(b = 0; ok && b <= circumscribe; b++) {
    float xcos = (float) cos((b * 2 * cPI) / circumscribe);
    float ysin = (float) sin((b * 2 * cPI) / circumscribe);
    float xcosr = xcos * radius;
    float ysinr = ysin * radius;
    int flag = true, i;
    v[0] = vp[0] + vx[0] * xcosr + vy[0] * ysinr;
    v[1] = vp[1] + vx[1] * xcosr + vy[1] * ysinr;
    v[2] = vp[2] + vx[2] * xcosr + vy[2] * ysinr;
    
    i = *(MapLocusEStart(map, v));
    if(i && map->EList) {
      int j = map->EList[i++];
      while(ok && j >= 0) {
	SurfaceJobAtomInfo *j_atom_info = atom_info + j;
	if((!present) || present[j])
	  if((j != a) && (j != jj)) {
	    int skip_flag = false;
	    if(a_atom_info->vdw == j_atom_info->vdw) {  /* handle singularities */
	      float *v1 = coord + 3 * j;
	      if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
		skip_flag = true;
	    }
	    if(jj_atom_info->vdw == j_atom_info->vdw) { /* handle singularities */
	      float *v1 = coord + 3 * j;
	      if((v2[0] == v1[0]) && (v2[1] == v1[1]) && (v2[2] == v1[2]))
		skip_flag = true;
	    }
	    if(!skip_flag)
	      if(within3f(coord + 3 * j, v,j_atom_info->vdw + probe_radius)) {
		flag = false;
		break;
	      }
	  }
	j = map->EList[i++];
	ok &= !G->Interrupt;
      }
    }
    if(ok && flag && (*dotCnt < stopDot)) {
      float vt0[3], vt2[3];
      subtract3f(v0, v, vt0);
      subtract3f(v2, v, vt2);
      normalize3f(vt0);
      normalize3f(vt2);
      add3f(vt0, vt2, n);
      invert3f(n);
      normalize3f(n);
      /*
	n[0] = vx[0] * xcos + vy[0] * ysin;
	n[1] = vx[1] * xcos + vy[1] * ysin;
	n[2] = vx[2] * xcos + vy[2] * ysin;
      */
      I->dotCode[I->nDot] = 1;  /* mark as exempt */
      (*dotCnt)++;
      v += 3;
      n += 3;
      I->nDot++;
    }
  }
  return ok;
}

static int SolventDotMarkDotsWithinCutoff(PyMOLGlobals * G, SolventDot *I,
    MapType *map, float *I_dot, int nDot, float *cavityDot, int *dot_flag, float cutoff){
  int ok = true, *p = dot_flag;
  float *v = I->dot;
  int a;
  for(a = 0; ok && a < I->nDot; a++) {
    int i = *(MapLocusEStart(map, v));
    if(i && map->EList) {
      int j = map->EList[i++];
      while(j >= 0) {
	if(within3f(cavityDot + (3 * j), v, cutoff)) {
	  *p = true;
	  break;
	}
	j = map->EList[i++];
      }
    }
    v += 3;
    p++;
    if(G->Interrupt) {
      ok = false;
    }
  }
  return ok;
}


static void SolventDotSlideDotsAndInfo(PyMOLGlobals * G, SolventDot *I, int *dot_flag, int flag_value){
  float *v0 = I->dot;
  float *n0 = I->dotNormal;
  int *dc0 = I->dotCode;
  int *p = dot_flag;
  int c = I->nDot;
  float *n = n0;
  float *v = v0;
  int *dc = dc0;
  int a;
  I->nDot = 0;
  for(a = 0; a < c; a++) {
    if(flag_value ? *(p++) : !*(p++)) {
      *(v0++) = *(v++);
      *(n0++) = *(n++);
      *(v0++) = *(v++);
      *(n0++) = *(n++);
      *(v0++) = *(v++);
      *(n0++) = *(n++);
      *(dc0++) = *(dc++);
      I->nDot++;
    } else {
      v += 3;
      n += 3;
    }
  }
  PRINTFD(G, FB_RepSurface)
    " SolventDotNew-DEBUG: %d->%d\n", c, I->nDot ENDFD;
}

static int SolventDotMarkDotsWithinProbeRadius(PyMOLGlobals * G, SolventDot *I,
    MapType *map, int cavity_cull, float probe_radius_plus, int *dot_flag, int *flag){
  int ok = true;
  int *p = dot_flag;
  float *v = I->dot;
  int a;
  for(a = 0; ok && a < I->nDot; a++) {
    if(!dot_flag[a]) {
      int i = *(MapLocusEStart(map, v));
      int cnt = 0;
      if(i && map->EList) {
	int j = map->EList[i++];
	while(j >= 0) {
	  if(j != a) {
	    if(within3f(I->dot + (3 * j), v, probe_radius_plus)) {
	      if(dot_flag[j]) {
		*p = true;
		(*flag) = true;
		break;
	      }
	      cnt++;
	      if(cnt > cavity_cull) {
		*p = true;
		(*flag) = true;
		break;
	      }
	    }
	  }
	  j = map->EList[i++];
	}
      }
    }
    v += 3;
    p++;
    ok &= !G->Interrupt;
  }
  ok &= !G->Interrupt;
  return ok;
}

static SolventDot *SolventDotNew(PyMOLGlobals * G,
                                 float *coord,
                                 SurfaceJobAtomInfo * atom_info,
                                 float probe_radius, SphereRec * sp,
                                 int *present,
                                 int circumscribe, int surface_mode,
                                 int surface_solvent, int cavity_cull,
                                 int all_visible_flag, float max_vdw,
                                 int cavity_mode, float cavity_radius, 
                                 float cavity_cutoff)
{
  int ok = true;
  int stopDot;
  int n_coord = VLAGetSize(atom_info);
  OOCalloc(G, SolventDot);
  CHECKOK(ok, I);
  /*  printf("%p %p %p %f %p %p %p %d %d %d %d %d %f\n",
     G,
     coord,
     atom_info,
     probe_radius,sp,
     extent,present,
     circumscribe,  surface_mode, 
     surface_solvent,  cavity_cull,
     all_visible_flag, max_vdw);
   */

  stopDot = n_coord * sp->nDot + 2 * circumscribe;
  if (ok)
    I->dot = VLAlloc(float, (stopDot + 1) * 3);
  CHECKOK(ok, I->dot);
  if (ok){
    I->dotNormal = VLAlloc(float, (stopDot + 1) * 3);
    CHECKOK(ok, I->dotNormal);
  }
  if (ok){
    I->dotCode = VLACalloc(int, stopDot + 1);
    CHECKOK(ok, I->dotCode);
  }
  I->nDot = 0;
  if (ok) {
    int dotCnt = 0;
    MapType *map = MapNewFlagged(G, max_vdw + probe_radius, coord, n_coord, NULL, present);
    CHECKOK(ok, map);
    ok &= !G->Interrupt;
    if(map && ok) {
      ok &= MapSetupExpress(map);
      if (ok) {
        int a;
        int skip_flag;
        SurfaceJobAtomInfo *a_atom_info = atom_info;
        for(a = 0; ok && a < n_coord; a++) {
          OrthoBusyFast(G, a, n_coord * 5);
          if((!present) || (present[a])) {
            skip_flag = false;
	    ok = SolventDotFilterOutSameXYZ(G, map, atom_info, a_atom_info, coord, a, present, &skip_flag);
            if(ok && !skip_flag) {
	      ok = SolventDotGetDotsAroundVertexInSphere(G, I, map, atom_info, a_atom_info, coord, a, present, sp, probe_radius, &dotCnt, stopDot, I->dot, I->dotNormal, &I->nDot);
            }
          }
          a_atom_info++;
        }
      }

      /* for each pair of proximal atoms, circumscribe a circle for their intersection */
      if (ok) {
        MapType *map2 = NULL;
        if(circumscribe && (!surface_solvent)){
          map2 = MapNewFlagged(G, 2 * (max_vdw + probe_radius), coord, n_coord, NULL, present);
	  CHECKOK(ok, map2);
	}
	ok &= !G->Interrupt;
        if(ok && map2) {
          int a;
          int skip_flag;

          SurfaceJobAtomInfo *a_atom_info = atom_info;
          ok &= MapSetupExpress(map2);
          for(a = 0; ok && a < n_coord; a++) {
            if((!present) || present[a]) {
              float *v0 = coord + 3 * a;
              skip_flag = false;

	      ok = SolventDotFilterOutSameXYZ(G, map2, atom_info, a_atom_info, coord, a, present, &skip_flag);
              if(ok && !skip_flag) {
                int ii = *(MapLocusEStart(map2, v0));
                if(ii) {
                  int jj = map2->EList[ii++];
		  float vdw[3];
		  vdw[0] = a_atom_info->vdw + probe_radius;
		  vdw[1] = vdw[0] * vdw[0];
                  while(ok && jj >= 0) {
                    SurfaceJobAtomInfo *jj_atom_info = atom_info + jj;
                    float dist;
                    if(jj > a)  /* only check if this is atom trails */
                      if((!present) || present[jj]) {
                        float *v2 = coord + 3 * jj;
			vdw[2] = jj_atom_info->vdw + probe_radius;
                        dist = (float) diff3f(v0, v2);
			if((dist > R_SMALL4) && (dist < (vdw[0] + vdw[2]))) {
			  ok = SolventDotCircumscribeAroundVertex(G, I, map, vdw, dist, v0, v2, circumscribe, 
								  atom_info, a_atom_info, jj_atom_info, 
								  present, a, jj, coord, probe_radius, &dotCnt, stopDot);
                        }
                      }
                    jj = map2->EList[ii++];
                  }
                }
              }
            }
            a_atom_info++;
	    ok &= !G->Interrupt;
          }
        }
        MapFree(map2);
      }
    }
    MapFree(map);
  }

  if(ok && cavity_mode) {
    int nCavityDot = 0;
    int dotCnt = 0;
    float *cavityDot = VLAlloc(float, (stopDot + 1) * 3);
    CHECKOK(ok, cavityDot);
    if(cavity_radius<0.0F) {
      cavity_radius = - probe_radius * cavity_radius;
    }
    if(cavity_cutoff<0.0F) {
      cavity_cutoff = cavity_radius - cavity_cutoff * probe_radius;
    }
    {
      MapType *map = MapNewFlagged(G, max_vdw + cavity_radius, coord, n_coord, NULL, present);
      CHECKOK(ok, map);
      if(G->Interrupt)
        ok = false;
      if(ok && map) {
        ok &= MapSetupExpress(map);
        if (ok) {
          int a;
          int skip_flag;
          SurfaceJobAtomInfo *a_atom_info = atom_info;
          for(a = 0; a < n_coord; a++) {
            if((!present) || (present[a])) {
              skip_flag = false;
	      ok = SolventDotFilterOutSameXYZ(G, map, atom_info, a_atom_info, coord, a, present, &skip_flag);
              if(ok && !skip_flag) {
		ok = SolventDotGetDotsAroundVertexInSphere(G, I, map, atom_info, a_atom_info, coord, a, present, sp, cavity_radius, &dotCnt, stopDot, cavityDot, NULL, &nCavityDot);
              }
            }
            a_atom_info++;
          }
        }
      }
      MapFree(map);
    }
    {
      int *dot_flag = pymol::calloc<int>(I->nDot);
      ErrChkPtr(G, dot_flag);
      {
        MapType *map = MapNew(G, cavity_cutoff, cavityDot, nCavityDot, NULL);
        if(map) {
          MapSetupExpress(map);
	  ok = SolventDotMarkDotsWithinCutoff(G, I, map, I->dot, I->nDot, cavityDot, dot_flag, cavity_cutoff);
        }
        MapFree(map);
      }
      SolventDotSlideDotsAndInfo(G, I, dot_flag, false);
      FreeP(dot_flag);
    }
    VLAFreeP(cavityDot);
  }

  if(ok && (cavity_mode != 1) && (cavity_cull > 0) && 
     (probe_radius > 0.75F) && (!surface_solvent)) {
    int *dot_flag = pymol::calloc<int>(I->nDot);
    float probe_radius_plus;
    probe_radius_plus = probe_radius * 1.5F;

    ErrChkPtr(G, dot_flag);

    {
      MapType *map = MapNew(G, probe_radius_plus, I->dot, I->nDot, NULL);
      if(map) {
        int flag = true;
        MapSetupExpress(map);
        while(ok && flag) {
          flag = false;
	  ok = SolventDotMarkDotsWithinProbeRadius(G, I, map, cavity_cull, probe_radius_plus, dot_flag, &flag);
        }
      }
      MapFree(map);
    }

    if (ok) {
      SolventDotSlideDotsAndInfo(G, I, dot_flag, true);
    }

    FreeP(dot_flag);
  }
  if(!ok) {
    SolventDotFree(I);
    I = NULL;
  }
  return I;
}
