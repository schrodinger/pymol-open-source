
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

typedef struct RepSurface {
  Rep R;
  int N;
  int NT;
  int proximity;
  float *V, *VN, *VC, *VA, *VAO; /* VAO - Ambient Occlusion per vertex */
  int *RC;
  int *Vis;
  int *T, *S, *AT;                   /* T=vertices, S=strips, AT=closest atom for vertices */
  int solidFlag;
  int oneColorFlag, oneColor;
  int allVisibleFlag;
  char *LastVisib;
  int *LastColor;
  int ColorInvalidated;
  int Type;
  float max_vdw;
  CGO *debug;

  /* These variables are for using the shader.  All of them */
  /* are allocated/set when generate_shader_cgo to minimize */
  /* allocation during the rendering loop. */
  CGO *shaderCGO, *pickingCGO;
  short dot_as_spheres;
  uint *vertexIndices; /* the vertex index order for transparency, computed each frame and passed to VBO */
  float *sum;  /* the sum of the x,y,z for each triangle. used to compute zvalue. */
  float *z_value; /* the z value for each triangle, computed each frame. */
  int n_tri; /* the number of triangles in the transparent surface */ 
  int *ix; /* the triangle index order for transparency, computed each frame. */
} RepSurface;

void RepSurfaceFree(RepSurface * I);
int RepSurfaceSameVis(RepSurface * I, CoordSet * cs);
void RepSurfaceColor(RepSurface * I, CoordSet * cs);

void RepSurfaceFree(RepSurface * I)
{
  VLAFreeP(I->V);
  VLAFreeP(I->VN);
  if (I->pickingCGO && I->pickingCGO!=I->shaderCGO){
    CGOFree(I->pickingCGO);
    I->pickingCGO = NULL;
  }
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
  }
  if (I->vertexIndices){
    FreeP(I->vertexIndices);
  }
  if (I->sum){
    FreeP(I->sum);
  }
  if (I->z_value){
    FreeP(I->z_value);
  }
  if (I->ix){
    FreeP(I->ix);
  }
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
  CGOFree(I->debug);
  VLAFreeP(I->T);
  VLAFreeP(I->S);
  VLAFreeP(I->AT);
  RepPurge(&I->R);              /* unnecessary, but a good idea */
  OOFreeP(I);
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

void RepSurfaceSortIX(PyMOLGlobals * G, RepSurface *I, int t_mode){
  float *z_value = NULL, *zv;
  float *sum, *sv;
  float matrix[16];
  int idx;
  int *ix = NULL;
  
  glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
  sum = I->sum;
  z_value = I->z_value;
  ix = I->ix;
  zv = z_value;
  sv = sum;
  
  /* for each triangle, computes the z */
  for (idx = 0; idx<I->n_tri; idx++){
    *(zv++) = matrix[2] * sv[0] + matrix[6] * sv[1] + matrix[10] * sv[2];
    sv += 3;
  }
  
  switch (t_mode) {
  case 1:
    UtilSemiSortFloatIndex(I->n_tri, z_value, ix, true);
    /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZOrderFn); */
    break;
  default:
    UtilSemiSortFloatIndex(I->n_tri, z_value, ix, false);
    /* UtilSortIndex(n_tri,z_value,ix,(UtilOrderFn*)ZRevOrderFn); */
    break;
  }
}

int AtomInfoIsMasked(ObjectMolecule *obj, int atm){
  AtomInfoType *ait = &obj->AtomInfo[atm];
  return (ait->masked ? cPickableNoPick : cPickableAtom);
}

static void RepSurfaceRender(RepSurface * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  float *vn = I->VN;
  float *vc = I->VC;
  float *va = I->VA;
  int *rc = I->RC;
  int *t = I->T;
  int *s = I->S;
  int c = I->N;
  int *vi = I->Vis;
  int *at = I->AT;
  int ok = true;
  float alpha;
  int t_mode;
  short dot_as_spheres = SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_as_spheres);
  float ambient_occlusion_scale = 0.f;
  int ambient_occlusion_mode = ambient_occlusion_mode = SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_ambient_occlusion_mode);
  int ambient_occlusion_mode_div_4 = 0;
  if (ambient_occlusion_mode){
    ambient_occlusion_scale = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_ambient_occlusion_scale);
    ambient_occlusion_mode_div_4 = ambient_occlusion_mode / 4;
  }

  if((I->Type != 1) && (!s)) {
    return;
  }

  alpha = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  if(ray) {
    ray->fTransparentf(ray, 1.0F - alpha);
    if(I->Type == 1) {
      /* dot surface */
      float radius;
      radius = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_radius);
      if(radius == 0.0F) {
        radius = ray->PixelRadius * SettingGet_f(G, I->R.cs->Setting,
                                                 I->R.obj->Setting,
                                                 cSetting_dot_width) / 1.4142F;
      }

      if(I->oneColorFlag) {
        float col[3];
        ColorGetEncoded(G, I->oneColor, col);
        ray->fColor3fv(ray, col);
      }

      if(c)
        while(ok && c--) {
          if(*vi) {
            if(!I->oneColorFlag) {
              ray->fColor3fv(ray, vc);
            }
            ok &= ray->fSphere3fv(ray, v, radius);
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
          if((I->proximity
              && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
             || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))){
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
	    ok &= ray->fTriangle3fv(ray, v + (*t) * 3, v + (*(t + 1)) * 3, v + (*(t + 2)) * 3,
				    vn + (*t) * 3, vn + (*(t + 1)) * 3, vn + (*(t + 2)) * 3,
				    col1, col2, col3);
	  }
          t += 3;
        }
      } else {
        while(ok && c--) {
          register int ttA = *t, ttB = *(t + 1), ttC = *(t + 2);
          if((I->proximity && ((*(vi + ttA)) || (*(vi + ttB)) || (*(vi + ttC)))) ||
             ((*(vi + ttA)) && (*(vi + ttB)) && (*(vi + ttC)))) {
            register int ttA3 = ttA * 3, ttB3 = ttB * 3, ttC3 = ttC * 3;
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
                ok &= ray->fTriangleTrans3fv(ray, v + ttA3, v + ttB3, v + ttC3,
					     vn + ttA3, vn + ttB3, vn + ttC3,
					     cA, cB, cC,
					     1.0F - va[ttA], 1.0F - va[ttB], 1.0F - va[ttC]);
              } else {
                ok &= ray->fTriangle3fv(ray, v + ttA3, v + ttB3, v + ttC3,
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
      int *cache = Calloc(int, spacing * (I->N + 1));
      CHECKOK(ok, cache);

      radius = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_mesh_radius);

      if(ok && radius == 0.0F) {
        float line_width =
          SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_mesh_width);
        line_width = SceneGetDynamicLineWidth(info, line_width);

        radius = ray->PixelRadius * line_width / 2.0F;
      }
      if (ok){
	c = I->NT;
	if(I->oneColorFlag) {
	  float col[3];
	  ColorGetEncoded(G, I->oneColor, col);
	  while(ok && c--) {
	    t0 = (*t);
	    t1 = (*(t + 1));
	    t2 = (*(t + 2));
	    if((I->proximity && ((*(vi + t0)) || (*(vi + t1)) || (*(vi + t2)))) ||
	       ((*(vi + t0)) && (*(vi + t1)) && (*(vi + t2)))) {
	      if(!check_and_add(cache, spacing, t0, t1))
		ok &= ray->fSausage3fv(ray, v + t0 * 3, v + t1 * 3, radius, col, col);
	      if(!check_and_add(cache, spacing, t1, t2))
		ok &= ray->fSausage3fv(ray, v + t1 * 3, v + t2 * 3, radius, col, col);
	      if(!check_and_add(cache, spacing, t2, t0))
		ok &= ray->fSausage3fv(ray, v + t2 * 3, v + t0 * 3, radius, col, col);
	    }
	    t += 3;
	  }
	} else {
	  while(ok && c--) {
	    t0 = (*t);
	    t1 = (*(t + 1));
	    t2 = (*(t + 2));
	    
	    if((I->proximity && ((*(vi + t0)) || (*(vi + t1)) || (*(vi + t2)))) ||
	       ((*(vi + t0)) && (*(vi + t1)) && (*(vi + t2))))
	      if((*(vi + t0)) || (*(vi + t1)) || (*(vi + t2))) {
		if(!check_and_add(cache, spacing, t0, t1))
		  ok &= ray->fSausage3fv(ray, v + t0 * 3, v + t1 * 3, radius, vc + t0 * 3,
					 vc + t1 * 3);
		if(!check_and_add(cache, spacing, t1, t2))
		  ok &= ray->fSausage3fv(ray, v + t1 * 3, v + t2 * 3, radius, vc + t1 * 3,
					 vc + t2 * 3);
		if(!check_and_add(cache, spacing, t2, t0))
		  ok &= ray->fSausage3fv(ray, v + t2 * 3, v + t0 * 3, radius, vc + t2 * 3,
					 vc + t0 * 3);
	      }
	    t += 3;
	  }
	}
	FreeP(cache);
      }
    }
    if (ok){
      ray->fTransparentf(ray, 0.0);
    } else {
      /* If not ok, then Clear Entire RepSurface, not just the ray object */
      
    }
  } else if(G->HaveGUI && G->ValidContext) {
    /* Not ray tracing, but rendering */
    if(pick) {
      int pick_surface = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_pick_surface);
      int no_pick_but_write_to_depth_buffer = (!pick_surface && (alpha == 1.0));
      if (I->pickingCGO && (pick_surface || no_pick_but_write_to_depth_buffer)){
	I->pickingCGO->use_shader = false;
	I->pickingCGO->no_pick = no_pick_but_write_to_depth_buffer;
	CGORenderGLPicking(I->pickingCGO, pick, &I->R.context, I->R.cs->Setting, I->R.obj->Setting);
      }
    } else {
      short use_shader, generate_shader_cgo = 0, use_display_lists = 0;
      use_shader = (int) SettingGet(G, cSetting_surface_use_shader) & 
                           (int) SettingGet(G, cSetting_use_shaders);
      use_display_lists = (int) SettingGet(G, cSetting_use_display_lists);

      if (I->shaderCGO && (!use_shader || CGOCheckWhetherToFree(G, I->shaderCGO) ||
			   (I->Type == 1 && I->dot_as_spheres != dot_as_spheres))){
	CGOFree(I->shaderCGO);
	I->shaderCGO = NULL;
      }
#ifdef _PYMOL_GL_CALLLISTS
        if(use_display_lists && I->R.displayList) {
          glCallList(I->R.displayList);
	  return;
	}
#endif
      
      if (use_shader){
	if (!I->shaderCGO){
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok){
	    I->shaderCGO->use_shader = true;
	    I->dot_as_spheres = dot_as_spheres;
	    generate_shader_cgo = 1;
	  }
	} else {
	  CShaderPrg * shaderPrg = 0;
	  if (I->Type == 1){
	    if (dot_as_spheres){
	      float radius = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_width) 
		             * info->vertex_scale;
	      shaderPrg = CShaderPrg_Enable_DefaultSphereShader(G);
	      CShaderPrg_Set1f(shaderPrg, "sphere_size_scale", fabs(radius));
	    } else {
	      shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	      SceneResetNormalUseShaderAttribute(G, 0, true, CShaderPrg_GetAttribLocation(shaderPrg, "a_Normal"));
	    }
	  } else if (I->Type == 2) {
	    float mesh_width = SettingGet_f(G, I->R.obj->Setting, NULL, cSetting_mesh_width);
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    CShaderPrg_Set1f(shaderPrg, "uni_radius", SceneGetLineWidthForCylinders(G, info, mesh_width));
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_SetLightingEnabled(shaderPrg, 1);
	    CShaderPrg_Set1f(shaderPrg, "ambient_occlusion_scale", ambient_occlusion_scale);
	    CShaderPrg_Set1i(shaderPrg, "use_interior_color_threshold", 1);
	    if((alpha != 1.0) || va) {
	      /* Updating indices if alpha */
	      t_mode =
		SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting,
			     cSetting_transparency_mode);
	      if(info && info->alpha_cgo) {
		t_mode = 0;
	      }
	      if(t_mode) {
		RepSurfaceSortIX(G, I, t_mode);
		/* Now once we have ix, we can update the vertexIndices that are stored
		   as VBOs */
		if (I->ix){
		  float *pc = CGOGetNextDrawBufferedIndex(I->shaderCGO->op);
		  int *ix = I->ix;
		  if (pc){
		    int nindices = CGO_get_int(pc+3), c, pl, idx;
		    uint vbuf = CGO_get_int(pc+8);
		    uint *vertexIndices;
		    vertexIndices = I->vertexIndices;
		    //		  vertexIndices = Alloc(uint, nindices);	      
		    if (!vertexIndices){
		      PRINTFB(I->R.G, FB_RepSurface, FB_Errors) "ERROR: RepSurfaceRender() vertexIndices is not set, nindices=%d\n", nindices ENDFB(I->R.G);
		    }
		    /* updates the vertexIndices from the ix array */
		    for(c = 0, pl=0; c < I->n_tri; c++) {
		      idx = ix[c] * 3;
		      vertexIndices[pl++] = idx;
		      vertexIndices[pl++] = idx + 1;
		      vertexIndices[pl++] = idx + 2;
		    }
		    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbuf);
		    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint)*nindices, vertexIndices, GL_STATIC_DRAW);
		  }
		}
	      }
	    }
	  }
	  {
	    float *color;
	    color = ColorGet(G, I->R.obj->Color);
	    I->shaderCGO->enable_shaders = shaderPrg ? 0 : 1;
	    CGORenderGL(I->shaderCGO, color, NULL, NULL, info, &I->R);
	  }
	  if (shaderPrg)
	    CShaderPrg_Disable(shaderPrg);
	  return;
	}
      }
#ifdef _PYMOL_GL_CALLLISTS
      if(use_display_lists) {
	if(!I->R.displayList) {
	  I->R.displayList = glGenLists(1);
	  if(I->R.displayList) {
	    glNewList(I->R.displayList, GL_COMPILE_AND_EXECUTE);
	  }
	}
      }
#endif


      if(I->debug)
        CGORenderGL(I->debug, NULL, NULL, NULL, info, &I->R);
      if(I->Type == 1) {
        /* no triangle information, so we're rendering dots only */
        int normals =
          SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_normals);
        int lighting =
          SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_lighting);
        if(!normals){
	  if (generate_shader_cgo){
	    CGOResetNormal(I->shaderCGO, true);
	  } else {
	    SceneResetNormal(G, true);
	  }
	}
        if(!lighting)
	  if(!info->line_lighting)
	    glDisable(GL_LIGHTING);
	if (generate_shader_cgo){
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
		  if (ok)
		    ok &= CGOPickColor(I->shaderCGO, *at, AtomInfoIsMasked((ObjectMolecule*)I->R.obj, *at));

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
			      (G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_width));
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
		  /*		  if(normals){
		    CGONormalv(I->shaderCGO, vn);
		    }*/
		  if (ok)
		    ok &= CGOPickColor(I->shaderCGO, *at, AtomInfoIsMasked((ObjectMolecule*)I->R.obj, *at));
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
	} else {
          glPointSize(SettingGet_f
                      (G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_width));
          if(c) {
            glColor3f(1.0, 0.0, 0.0);
#ifdef _PYMOL_GL_DRAWARRAYS
            if(I->oneColorFlag) {
              glColor3fv(ColorGet(G, I->oneColor));
            }
	    {
	      int nverts = 0, cinit = c;
	      {
		while(c--) {
		  if(*vi) {
		    nverts++;
		  }
		}
	      }
	      c = cinit; 
	      {
		int pl = 0, plc = 0;
		ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
		ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
		ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
		while(c--) {
		  if(*vi) {
		    if(!I->oneColorFlag) {
		      colorVals[plc++] = vc[0]; colorVals[plc++] = vc[1]; colorVals[plc++] = vc[2]; colorVals[plc++] = 1.f;
		    }
		    if(normals){
		      normalVals[pl] = vn[0]; normalVals[pl+1] = vn[1]; normalVals[pl+2] = vn[2];
		    }
		    ptsVals[pl] = v[0]; ptsVals[pl+1] = v[1]; ptsVals[pl+2] = v[2];
		  }
		  vi++;
		  vc += 3;
		  vn += 3;
		  v += 3;
		  pl += 3;
		}
		glEnableClientState(GL_VERTEX_ARRAY);
		if(!I->oneColorFlag) {
		  glEnableClientState(GL_COLOR_ARRAY);
		}
		if (normals){
		  glEnableClientState(GL_NORMAL_ARRAY);
		  glNormalPointer(GL_FLOAT, 0, normalVals);
		}
		glVertexPointer(3, GL_FLOAT, 0, ptsVals);
		if(!I->oneColorFlag) {
		  glColorPointer(4, GL_FLOAT, 0, colorVals);
		}
		glDrawArrays(GL_POINTS, 0, nverts);
		glDisableClientState(GL_VERTEX_ARRAY);
		if(!I->oneColorFlag) {
		  glDisableClientState(GL_COLOR_ARRAY);
		}
		if (normals){
		  glDisableClientState(GL_NORMAL_ARRAY);
		}
		DEALLOCATE_ARRAY(ptsVals)
		DEALLOCATE_ARRAY(colorVals)
		DEALLOCATE_ARRAY(normalVals)
	      }
	    }
#else
            glBegin(GL_POINTS);
            if(I->oneColorFlag) {
              glColor3fv(ColorGet(G, I->oneColor));
            }

            while(c--) {
              if(*vi) {
                if(!I->oneColorFlag) {
                  glColor3fv(vc);
                }
                if(normals)
                  glNormal3fv(vn);
                glVertex3fv(v);
              }
              vi++;
              vc += 3;
              vn += 3;
              v += 3;
            }
            glEnd();
#endif
          }
	  } /* else use shader */
	  if(!lighting)
	    glEnable(GL_LIGHTING);
      } else if(I->Type == 2) { /* rendering triangle mesh */
        int normals =
          SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_mesh_normals);
        int lighting =
          SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_mesh_lighting);
        if(ok && !normals){
	  if (generate_shader_cgo){
	    ok &= CGOResetNormal(I->shaderCGO, true);
	  } else {
	    SceneResetNormal(G, true);
	  }
	}
        if(!lighting)
          if(!info->line_lighting)
            glDisable(GL_LIGHTING);
	if (ok) {
          float line_width =
            SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_mesh_width);
          line_width = SceneGetDynamicLineWidth(info, line_width);

	  if (generate_shader_cgo){
	    ok &= CGOLinewidthSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_MESH);
	  } else {
	    glLineWidth(line_width);
	  }

          c = I->NT;
          if(ok && c) {
            if(I->oneColorFlag) {
		if (generate_shader_cgo){
		  ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		  while(ok && c--) {
		    if((I->proximity
			&& ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		       || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))) {
		      if(normals) {
			int idx;
			ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
			idx = (*(t + 2));
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			
			idx = (*t);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
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
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			idx = *t;
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = *t;
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = *t;
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			if (ok)
			  ok &= CGOEnd(I->shaderCGO);
		      }
		    } else
		      t += 3;
		  }
	      } else {
              glColor3fv(ColorGet(G, I->oneColor));
              while(c--) {
                if((I->proximity
                    && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
                   || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))) {
                  if(normals) {
#ifdef _PYMOL_GL_DRAWARRAYS
		    {
		      float *tmp_ptr;
		      ALLOCATE_ARRAY(GLfloat,vertVals,4*3)
		      ALLOCATE_ARRAY(GLfloat,normVals,4*3)
		      tmp_ptr = v + (*(t + 2)) * 3;
		      vertVals[0] = tmp_ptr[0]; vertVals[1] = tmp_ptr[1]; vertVals[2] = tmp_ptr[2];
		      tmp_ptr = vn + (*(t + 2)) * 3;
		      normVals[0] = tmp_ptr[0]; normVals[1] = tmp_ptr[1]; normVals[2] = tmp_ptr[2];
		      tmp_ptr = v + (*t) * 3;
		      vertVals[3] = tmp_ptr[0]; vertVals[4] = tmp_ptr[1]; vertVals[5] = tmp_ptr[2];
		      tmp_ptr = vn + (*t) * 3;
		      normVals[3] = tmp_ptr[0]; normVals[4] = tmp_ptr[1]; normVals[5] = tmp_ptr[2];
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[6] = tmp_ptr[0]; vertVals[7] = tmp_ptr[1]; vertVals[8] = tmp_ptr[2];
		      tmp_ptr = vn + (*t) * 3;
		      normVals[6] = tmp_ptr[0]; normVals[7] = tmp_ptr[1]; normVals[8] = tmp_ptr[2];
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[9] = tmp_ptr[0]; vertVals[10] = tmp_ptr[1]; vertVals[11] = tmp_ptr[2];
		      tmp_ptr = vn + (*t) * 3;
		      normVals[9] = tmp_ptr[0]; normVals[10] = tmp_ptr[1]; normVals[11] = tmp_ptr[2];
		      t++;
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_NORMAL_ARRAY);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glNormalPointer(GL_FLOAT, 0, normVals);
		      glDrawArrays(GL_LINE_STRIP, 0, 4);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_NORMAL_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		      DEALLOCATE_ARRAY(normVals)
		    }
#else
                    glBegin(GL_LINE_STRIP);

                    glNormal3fv(vn + (*(t + 2)) * 3);
                    glVertex3fv(v + (*(t + 2)) * 3);

                    glNormal3fv(vn + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glNormal3fv(vn + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glNormal3fv(vn + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glEnd();
#endif
                  } else {
#ifdef _PYMOL_GL_DRAWARRAYS
		    {
		      float *tmp_ptr;
		      ALLOCATE_ARRAY(GLfloat,vertVals,4*3)
		      tmp_ptr = v + (*(t + 2)) * 3;
		      vertVals[0] = tmp_ptr[0]; vertVals[1] = tmp_ptr[1]; vertVals[2] = tmp_ptr[2];
		      tmp_ptr = v + (*t) * 3;
		      vertVals[3] = tmp_ptr[0]; vertVals[4] = tmp_ptr[1]; vertVals[5] = tmp_ptr[2];
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[6] = tmp_ptr[0]; vertVals[7] = tmp_ptr[1]; vertVals[8] = tmp_ptr[2];
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[9] = tmp_ptr[0]; vertVals[10] = tmp_ptr[1]; vertVals[11] = tmp_ptr[2];
		      t++;
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glDrawArrays(GL_LINE_STRIP, 0, 4);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		    }
#else
                    glBegin(GL_LINE_STRIP);

                    glVertex3fv(v + (*(t + 2)) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glEnd();
#endif
                  }
                } else
                  t += 3;
              }
	      } /* end else use_shader */ 
            } else { /* not oneColorFlag */
		if (generate_shader_cgo){
		  while(ok && c--) {
		    if((I->proximity
			&& ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		       || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))) {
		      if(normals) {
			int idx;
			ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
			
			idx = (*(t + 2));
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGONormalv(I->shaderCGO, vn + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
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
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);

			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			idx = (*t);
			if (ok)
			  ok &= CGOColorv(I->shaderCGO, vc + idx * 3);
			if (ok)
			  ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
			if (ok)
			  ok &= CGOVertexv(I->shaderCGO, v + idx * 3);
			t++;
			if (ok)
			  ok &= CGOEnd(I->shaderCGO);
		      }
		    } else
		      t += 3;
		  }
	      } else {
              while(c--) {
                if((I->proximity
                    && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
                   || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))) {
                  if(normals) {
#ifdef _PYMOL_GL_DRAWARRAYS
		    {
		      float *tmp_ptr;
		      ALLOCATE_ARRAY(GLfloat,colorVals,4*4)
		      ALLOCATE_ARRAY(GLfloat,normVals,4*3)
		      ALLOCATE_ARRAY(GLfloat,vertVals,4*3)
		      tmp_ptr = v + (*(t + 2)) * 3;
		      vertVals[0] = tmp_ptr[0]; vertVals[1] = tmp_ptr[1]; vertVals[2] = tmp_ptr[2];
		      tmp_ptr = vn + (*(t + 2)) * 3;
		      normVals[0] = tmp_ptr[0]; normVals[1] = tmp_ptr[1]; normVals[2] = tmp_ptr[2];
		      tmp_ptr = vc + (*(t + 2)) * 3;
		      colorVals[0] = tmp_ptr[0]; colorVals[1] = tmp_ptr[1]; colorVals[2] = tmp_ptr[2]; colorVals[3] = 1.f;

		      tmp_ptr = v + (*t) * 3;
		      vertVals[3] = tmp_ptr[0]; vertVals[4] = tmp_ptr[1]; vertVals[5] = tmp_ptr[2];
		      tmp_ptr = vn + (*t) * 3;
		      normVals[3] = tmp_ptr[0]; normVals[4] = tmp_ptr[1]; normVals[5] = tmp_ptr[2];
		      tmp_ptr = vc + (*t) * 3;
		      colorVals[4] = tmp_ptr[0]; colorVals[5] = tmp_ptr[1]; colorVals[6] = tmp_ptr[2]; colorVals[7] = 1.f;
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[6] = tmp_ptr[0]; vertVals[7] = tmp_ptr[1]; vertVals[8] = tmp_ptr[2];
		      tmp_ptr = vn + (*t) * 3;
		      normVals[6] = tmp_ptr[0]; normVals[7] = tmp_ptr[1]; normVals[8] = tmp_ptr[2];
		      tmp_ptr = vc + (*t) * 3;
		      colorVals[8] = tmp_ptr[0]; colorVals[9] = tmp_ptr[1]; colorVals[10] = tmp_ptr[2]; colorVals[11] = 1.f;
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[9] = tmp_ptr[0]; vertVals[10] = tmp_ptr[1]; vertVals[11] = tmp_ptr[2];
		      tmp_ptr = vn + (*t) * 3;
		      normVals[9] = tmp_ptr[0]; normVals[10] = tmp_ptr[1]; normVals[11] = tmp_ptr[2];
		      tmp_ptr = vc + (*t) * 3;
		      colorVals[12] = tmp_ptr[0]; colorVals[13] = tmp_ptr[1]; colorVals[14] = tmp_ptr[2]; colorVals[15] = 1.f;
		      t++;
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_NORMAL_ARRAY);
		      glEnableClientState(GL_COLOR_ARRAY);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glColorPointer(4, GL_FLOAT, 0, colorVals);
		      glNormalPointer(GL_FLOAT, 0, normVals);
		      glDrawArrays(GL_LINE_STRIP, 0, 4);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_NORMAL_ARRAY);
		      glDisableClientState(GL_COLOR_ARRAY);
		      DEALLOCATE_ARRAY(colorVals)
		      DEALLOCATE_ARRAY(normVals)
		      DEALLOCATE_ARRAY(vertVals)
		    }
#else
                    glBegin(GL_LINE_STRIP);

                    glColor3fv(vc + (*(t + 2)) * 3);
                    glNormal3fv(vn + (*(t + 2)) * 3);
                    glVertex3fv(v + (*(t + 2)) * 3);

                    glColor3fv(vc + (*t) * 3);
                    glNormal3fv(vn + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glColor3fv(vc + (*t) * 3);
                    glNormal3fv(vn + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glColor3fv(vc + (*t) * 3);
                    glNormal3fv(vn + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glEnd();
#endif
                  } else {
#ifdef _PYMOL_GL_DRAWARRAYS
		    {
		      float *tmp_ptr;
		      ALLOCATE_ARRAY(GLfloat,vertVals,4*3)
		      ALLOCATE_ARRAY(GLfloat,colorVals,4*4)
		      tmp_ptr = v + (*(t + 2)) * 3;
		      vertVals[0] = tmp_ptr[0]; vertVals[1] = tmp_ptr[1]; vertVals[2] = tmp_ptr[2];
		      tmp_ptr = vc + (*(t + 2)) * 3;
		      colorVals[0] = tmp_ptr[0]; colorVals[1] = tmp_ptr[1]; colorVals[2] = tmp_ptr[2]; colorVals[3] = 1.f;

		      tmp_ptr = v + (*t) * 3;
		      vertVals[3] = tmp_ptr[0]; vertVals[4] = tmp_ptr[1]; vertVals[5] = tmp_ptr[2];
		      tmp_ptr = vc + (*t) * 3;
		      colorVals[4] = tmp_ptr[0]; colorVals[5] = tmp_ptr[1]; colorVals[6] = tmp_ptr[2]; colorVals[7] = 1.f;
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[6] = tmp_ptr[0]; vertVals[7] = tmp_ptr[1]; vertVals[8] = tmp_ptr[2];
		      tmp_ptr = vc + (*t) * 3;
		      colorVals[8] = tmp_ptr[0]; colorVals[9] = tmp_ptr[1]; colorVals[10] = tmp_ptr[2]; colorVals[11] = 1.f;
		      t++;
		      tmp_ptr = v + (*t) * 3;
		      vertVals[9] = tmp_ptr[0]; vertVals[10] = tmp_ptr[1]; vertVals[11] = tmp_ptr[2];
		      tmp_ptr = vc + (*t) * 3;
		      colorVals[12] = tmp_ptr[0]; colorVals[13] = tmp_ptr[1]; colorVals[14] = tmp_ptr[2]; colorVals[15] = 1.f;
		      t++;
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_COLOR_ARRAY);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glColorPointer(4, GL_FLOAT, 0, colorVals);
		      glDrawArrays(GL_LINE_STRIP, 0, 4);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_COLOR_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		      DEALLOCATE_ARRAY(colorVals)
		    }
#else
                    glBegin(GL_LINE_STRIP);

                    glColor3fv(vc + (*(t + 2)) * 3);
                    glVertex3fv(v + (*(t + 2)) * 3);

                    glColor3fv(vc + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glColor3fv(vc + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glColor3fv(vc + (*t) * 3);
                    glVertex3fv(v + (*t) * 3);
                    t++;
                    glEnd();
#endif
                  }
                } else
                  t += 3;
              }
            }
	    } /* end else use_shader */ 
          }
	}
        if(!lighting)
          glEnable(GL_LIGHTING);
      } else {
        /* we're rendering triangles */

        if((alpha != 1.0) || va) {

          t_mode =
            SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting,
                         cSetting_transparency_mode);

          if(info && info->alpha_cgo) {
            t_mode = 0;
          }

          if(t_mode) {

            float **t_buf = NULL, **tb;
            float *z_value = NULL, *zv;
            int *ix = NULL;
	    float *sumarray = NULL;
            int n_tri = 0, sumarraypl = 0;
            float sum[3];
            float matrix[16];

            glGetFloatv(GL_MODELVIEW_MATRIX, matrix);

            if(I->oneColorFlag) {
              t_buf = Alloc(float *, I->NT * 6);
            } else {
              t_buf = Alloc(float *, I->NT * 12);
            }
	    CHECKOK(ok, t_buf);
	    if (ok){
	      z_value = Alloc(float, I->NT);
	      CHECKOK(ok, z_value);
	    }
	    if (ok){
	      ix = Alloc(int, I->NT);
	      CHECKOK(ok, ix);
	    }
	    if (ok && use_shader && generate_shader_cgo){
	      sumarray = Alloc(float, I->NT * 3);
	      CHECKOK(ok, sumarray);
	    }
            zv = z_value;
            tb = t_buf;
            c = I->NT;
	    if (ok){
	      if(I->oneColorFlag) {
		while(c--) {
		  if((I->proximity
		      && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		     || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))) {
		    
		    *(tb++) = vn + (*t) * 3;
		    *(tb++) = v + (*t) * 3;
		    *(tb++) = vn + (*(t + 1)) * 3;
		    *(tb++) = v + (*(t + 1)) * 3;
		    *(tb++) = vn + (*(t + 2)) * 3;
		    *(tb++) = v + (*(t + 2)) * 3;
		    
		    add3f(tb[-1], tb[-3], sum);
		    add3f(sum, tb[-5], sum);
		    
		    if (sumarray){
		      sumarray[sumarraypl++] = sum[0];
		      sumarray[sumarraypl++] = sum[1];
		      sumarray[sumarraypl++] = sum[2];
		    }
		    *(zv++) = matrix[2] * sum[0] + matrix[6] * sum[1] + matrix[10] * sum[2];
		    n_tri++;
		  }
		  t += 3;
		}
	      } else {
		while(c--) {
		  if((I->proximity
		      && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		     || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))
		    if((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))) {
		      
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
		      if (sumarray){
			sumarray[sumarraypl++] = sum[0];
			sumarray[sumarraypl++] = sum[1];
			sumarray[sumarraypl++] = sum[2];
		      }
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
		if (generate_shader_cgo){
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
		  if (ok)
		    ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok)
		      ok &= CGOVertexv(I->shaderCGO, *(tb++));
		    if (ok)
		      ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
		    if (ok)
		      ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok)
		      ok &= CGOVertexv(I->shaderCGO, *(tb++));
		    if (ok)
		      ok &= CGONormalv(I->shaderCGO, *(tb++));
		    idx = ((*tb - v)/3);
		    if (ok)
		      ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
		    if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok)
		      ok &= CGOVertexv(I->shaderCGO, *(tb++));
		  }
		  if (ok)
		    ok &= CGOEnd(I->shaderCGO);
	      } else {
		glColor4f(col[0], col[1], col[2], alpha);
#ifdef _PYMOL_GL_DRAWARRAYS
		{
		  int nverts = n_tri*3;
		  int pl;
		  ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
		  ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
		  float *tmp_ptr;
		  pl = 0;
		  for(c = 0; c < n_tri; c++) {
		    tb = t_buf + 6 * ix[c];
		    tmp_ptr = *(tb++);
		    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		    tmp_ptr = *(tb++);
		    vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		    pl += 3;
		    tmp_ptr = *(tb++);
		    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		    tmp_ptr = *(tb++);
		    vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		    pl += 3;
		    tmp_ptr = *(tb++);
		    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		    tmp_ptr = *(tb++);
		    vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		    pl += 3;
		  }
		  
		  glEnableClientState(GL_VERTEX_ARRAY);
		  glEnableClientState(GL_NORMAL_ARRAY);
		  glVertexPointer(3, GL_FLOAT, 0, vertVals);
		  glNormalPointer(GL_FLOAT, 0, normVals);
		  glDrawArrays(GL_TRIANGLES, 0, nverts);
		  glDisableClientState(GL_VERTEX_ARRAY);
		  glDisableClientState(GL_NORMAL_ARRAY);
		  DEALLOCATE_ARRAY(vertVals)
		  DEALLOCATE_ARRAY(normVals)
		}
#else
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
#endif
	      } /* end if else use_shader */
	      }
            } else { /* else I->oneColorFlag */
		if (generate_shader_cgo){
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
		    if (ok) 
		      ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
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
		    if (ok) 
		      ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
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
		    if (ok) 
		      ok &= CGOPickColor(I->shaderCGO, I->AT[idx], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[idx]));
		    if (ok && I->VAO){
		      ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + idx));
		    }
		    if (ok) ok &= CGOVertexv(I->shaderCGO, *(tb++));
		  }
		  if (ok) ok &= CGOEnd(I->shaderCGO);
	      } else {
#ifdef _PYMOL_GL_DRAWARRAYS
		{
		  int nverts = n_tri*3;
		  int pl, plc;
		  float *tmp_ptr;
		  ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
		  ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
		  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
		  pl = 0; plc = 0;
		  for(c = 0; c < n_tri; c++) {
		    tb = t_buf + 12 * ix[c];
		    colorVals[plc+3] = **(tb++);
		    tmp_ptr = *(tb++);
		    colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; 
		    plc++;
		    tmp_ptr = *(tb++);
		    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		    tmp_ptr = *(tb++);
		    vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		    pl += 3;
		    
		    colorVals[plc+3] = **(tb++);
		    tmp_ptr = *(tb++);
		    colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; 
		    plc++;
		    tmp_ptr = *(tb++);
		    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		    tmp_ptr = *(tb++);
		    vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		    pl += 3;
		    
		    colorVals[plc+3] = **(tb++);
		    tmp_ptr = *(tb++);
		    colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; 
		    plc++;
		    tmp_ptr = *(tb++);
		    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		    tmp_ptr = *(tb++);
		    vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		    pl += 3;
		  }
		  
		  glEnableClientState(GL_VERTEX_ARRAY);
		  glEnableClientState(GL_NORMAL_ARRAY);
		  glEnableClientState(GL_COLOR_ARRAY);
		  glVertexPointer(3, GL_FLOAT, 0, vertVals);
		  glColorPointer(4, GL_FLOAT, 0, colorVals);
		  glNormalPointer(GL_FLOAT, 0, normVals);
		  glDrawArrays(GL_TRIANGLES, 0, nverts);
		  glDisableClientState(GL_VERTEX_ARRAY);
		  glDisableClientState(GL_COLOR_ARRAY);
		  glDisableClientState(GL_NORMAL_ARRAY);
		  DEALLOCATE_ARRAY(vertVals)
		  DEALLOCATE_ARRAY(normVals)
		  DEALLOCATE_ARRAY(colorVals)
		}
#else
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
#endif
	      }
	    } /* end if else use_shader */
	    if (ok && use_shader && generate_shader_cgo){
	      I->ix = ix;
	      I->z_value = z_value;
	      I->n_tri = n_tri;
	      I->sum = sumarray;
	    } else {
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
                      if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
                                           || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 1))))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 2))))))
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
                      if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
                                           || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 1))))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 2))))))
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
              /*          glCullFace(GL_BACK);
                 glEnable(GL_CULL_FACE);
                 glDepthMask(GL_FALSE); */
              if(I->allVisibleFlag) {
                if(I->oneColorFlag) {
                  float col[3];
                  ColorGetEncoded(G, I->oneColor, col);

		    if (generate_shader_cgo){
		      if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
		      if (ok) ok &= CGOColor(I->shaderCGO, col[0], col[1], col[2]);
		      c = *(s++);
		      while(ok && c) {
			if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			while(ok && c--) {
			  ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		  } else {
		    glColor4f(col[0], col[1], col[2], alpha);
		    c = *(s++);
		    while(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		      {
			int nverts = c + 2;
			int pl = 0;
			float *tmp_ptr;
			ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
			ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			while(c--) {
			  tmp_ptr = vn + (*s) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			  tmp_ptr = v + (*s) * 3;
			  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			  s++;
			}
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, vertVals);
			glNormalPointer(GL_FLOAT, 0, normVals);
			glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			DEALLOCATE_ARRAY(vertVals)
			DEALLOCATE_ARRAY(normVals)
		      }
#else
		      glBegin(GL_TRIANGLE_STRIP);
		      glNormal3fv(vn + (*s) * 3);
		      glVertex3fv(v + (*s) * 3);
		      s++;
		      glNormal3fv(vn + (*s) * 3);
		      glVertex3fv(v + (*s) * 3);
		      s++;
		      while(c--) {
			glNormal3fv(vn + (*s) * 3);
			glVertex3fv(v + (*s) * 3);
			s++;
		      }
		      glEnd();
#endif
		      c = *(s++);
		    }
		  } /* end if else use_shader */
		} else { /* I->oneColorFlag */
		    if (generate_shader_cgo){
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
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
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
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
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
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		  } else {
		    c = *(s++);
		    while(c) {
		      float *col;
#ifdef _PYMOL_GL_DRAWARRAYS
		      {
			int nverts = 2 + c;
			int pl = 0, plc = 0;
			float *tmp_ptr;
			ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
			ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
			ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
			col = vc + (*s) * 3;
			colorVals[plc++] = col[0]; colorVals[plc++] = col[1]; colorVals[plc++] = col[2]; 
			if(va) {
			  colorVals[plc++] = va[(*s)];
			} else {
			  colorVals[plc++] = alpha;
			}
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			col = vc + (*s) * 3;
			colorVals[plc++] = col[0]; colorVals[plc++] = col[1]; colorVals[plc++] = col[2]; 
			if(va) {
			  colorVals[plc++] = va[(*s)];
			} else {
			  colorVals[plc++] = alpha;
			}
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			while(c--) {
			  col = vc + (*s) * 3;
			  colorVals[plc++] = col[0]; colorVals[plc++] = col[1]; colorVals[plc++] = col[2]; 
			  if(va) {
			    colorVals[plc++] = va[(*s)];
			  } else {
			    colorVals[plc++] = alpha;
			  }
			  tmp_ptr = vn + (*s) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			  tmp_ptr = v + (*s) * 3;
			  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			  s++;
			}
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, vertVals);
			glColorPointer(4, GL_FLOAT, 0, colorVals);
			glNormalPointer(GL_FLOAT, 0, normVals);
			glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			DEALLOCATE_ARRAY(vertVals)
			DEALLOCATE_ARRAY(normVals)
			DEALLOCATE_ARRAY(colorVals)
		      }
#else
		      glBegin(GL_TRIANGLE_STRIP);
		      col = vc + (*s) * 3;
		      if(va) {
			glColor4f(col[0], col[1], col[2], va[(*s)]);
		      } else {
			glColor4f(col[0], col[1], col[2], alpha);
		      }
		      glNormal3fv(vn + (*s) * 3);
		      glVertex3fv(v + (*s) * 3);
		      s++;
		      col = vc + (*s) * 3;
		      if(va) {
			glColor4f(col[0], col[1], col[2], va[(*s)]);
		      } else {
			glColor4f(col[0], col[1], col[2], alpha);
		      }
		      glNormal3fv(vn + (*s) * 3);
		      glVertex3fv(v + (*s) * 3);
		      s++;
		      while(c--) {
			col = vc + (*s) * 3;
			if(va) {
			  glColor4f(col[0], col[1], col[2], va[(*s)]);
			} else {
			  glColor4f(col[0], col[1], col[2], alpha);
			}
			glNormal3fv(vn + (*s) * 3);
			glVertex3fv(v + (*s) * 3);
			s++;
		      }
		      glEnd();
#endif
		      c = *(s++);
		    }
		  }
		} /* end if else use_shader */
              } else {          /* subset s */
		  if (generate_shader_cgo){
		    c = I->NT;
		    if(ok && c) {
		      if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
		      if(I->oneColorFlag) {
			float color[3];
			float *col;
			ColorGetEncoded(G, I->oneColor, color);
			if (ok) ok &= CGOAlpha(I->shaderCGO, alpha);
			if (ok) ok &= CGOColor(I->shaderCGO, color[0], color[1], color[2]);
			while(ok && c--) {
			  if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					       || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									     &&
									     (*
									      (vi +
									       (*(t + 1))))
									     &&
									     (*
									      (vi +
									       (*(t + 2))))))
			    {
			      col = vc + (*t) * 3;
			      CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
			      if (ok) 
				ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			      col = vc + (*t) * 3;
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
			      if (ok) 
				ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			      col = vc + (*t) * 3;
			      if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			      if (ok && I->VAO){
				ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			      }
			      if (ok) 
				ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			    } else
			    t += 3;
			}
		      } else {
			float *col;
			while(ok && c--) {
			  if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					       || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									     &&
									     (*
									      (vi +
									       (*(t + 1))))
									     &&
									     (*
									      (vi +
									       (*(t + 2))))))
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
			      if (ok) 
				ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
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
			      if (ok) 
				ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
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
			      if (ok) 
				ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			      if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			      t++;
			    } else
			    t += 3;
			}
		      }
		      if (ok) ok &= CGOEnd(I->shaderCGO);
		    }
		} else {  /* end generate_shader_cgo */
                c = I->NT;
                if(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		  {
		    int nverts = 0, cinit = c;
		    {
		      while(c--) {
			if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					     || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									   &&
									   (*
									    (vi +
									     (*(t + 1))))
									   &&
									   (*
									    (vi +
									     (*(t + 2))))))
			  {
			    nverts += 3;
			  }
			t += 3;
		      }
		    }
		    {
		      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
		      int pl = 0, plc = 0;
		      c = cinit;
		      if(I->oneColorFlag) {
			float color[3];
			float *col, *tmp_ptr;
			ColorGetEncoded(G, I->oneColor, color);
			
			while(c--) {
			  if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					       || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									     &&
									     (*
									      (vi +
									       (*(t + 1))))
									     &&
									     (*
									      (vi +
									       (*(t + 2)))))) {
			      col = vc + (*t) * 3;
			      tmp_ptr = vn + (*t) * 3;
			      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			      tmp_ptr = v + (*t) * 3;
			      vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			      colorVals[plc++] = color[0]; colorVals[plc++] = color[1]; colorVals[plc++] = color[2]; 
			      colorVals[plc++] = alpha;
			      t++;
			      col = vc + (*t) * 3;
			      tmp_ptr = vn + (*t) * 3;
			      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			      tmp_ptr = v + (*t) * 3;
			      vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			      colorVals[plc++] = color[0]; colorVals[plc++] = color[1]; colorVals[plc++] = color[2]; 
			      colorVals[plc++] = alpha;
			      t++;
			      col = vc + (*t) * 3;
			      tmp_ptr = vn + (*t) * 3;
			      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			      tmp_ptr = v + (*t) * 3;
			      vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			      colorVals[plc++] = color[0]; colorVals[plc++] = color[1]; colorVals[plc++] = color[2]; 
			      colorVals[plc++] = alpha;
			      t++;
			  } else {
			    t += 3;
			  }
			}
		      } else {
			float *col, *tmp_ptr;
			while(c--) {
			  if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					       || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									     &&
									     (*
									      (vi +
									       (*(t + 1))))
									     &&
									     (*
									      (vi +
									       (*(t + 2)))))){
			    col = vc + (*t) * 3;
			    tmp_ptr = vn + (*t) * 3;
			    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			    tmp_ptr = v + (*t) * 3;
			    vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			    colorVals[plc++] = col[0]; colorVals[plc++] = col[1]; colorVals[plc++] = col[2];
			    if(va) {
			      colorVals[plc++] = va[(*t)];
			    } else {
			      colorVals[plc++] = alpha;
			    }
			    t++;
			    col = vc + (*t) * 3;
			    tmp_ptr = vn + (*t) * 3;
			    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			    tmp_ptr = v + (*t) * 3;
			    vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			    colorVals[plc++] = col[0]; colorVals[plc++] = col[1]; colorVals[plc++] = col[2];
			    if(va) {
			      colorVals[plc++] = va[(*t)];
			    } else {
			      colorVals[plc++] = alpha;
			    }
			    t++;

			    col = vc + (*t) * 3;
			    tmp_ptr = vn + (*t) * 3;
			    normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			    tmp_ptr = v + (*t) * 3;
			    vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			    colorVals[plc++] = col[0]; colorVals[plc++] = col[1]; colorVals[plc++] = col[2];
			    if(va) {
			      colorVals[plc++] = va[(*t)];
			    } else {
			      colorVals[plc++] = alpha;
			    }
			    t++;
			  } else {
			    t += 3;
			  }
			}
		      }
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_NORMAL_ARRAY);
		      glEnableClientState(GL_COLOR_ARRAY);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glColorPointer(4, GL_FLOAT, 0, colorVals);
		      glNormalPointer(GL_FLOAT, 0, normVals);
		      glDrawArrays(GL_TRIANGLES, 0, nverts);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_NORMAL_ARRAY);
		      glDisableClientState(GL_COLOR_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		      DEALLOCATE_ARRAY(normVals)
		      DEALLOCATE_ARRAY(colorVals)
		    }
		  }
#else
                  glBegin(GL_TRIANGLES);

                  if(I->oneColorFlag) {
                    float color[3];
                    float *col;
                    ColorGetEncoded(G, I->oneColor, color);

                    glColor4f(color[0], color[1], color[2], alpha);
                    while(c--) {
                      if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
                                           || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 1))))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 2))))))
                      {

                        col = vc + (*t) * 3;
                        glNormal3fv(vn + (*t) * 3);
                        glVertex3fv(v + (*t) * 3);
                        t++;
                        col = vc + (*t) * 3;
                        glNormal3fv(vn + (*t) * 3);
                        glVertex3fv(v + (*t) * 3);
                        t++;
                        col = vc + (*t) * 3;
                        glNormal3fv(vn + (*t) * 3);
                        glVertex3fv(v + (*t) * 3);
                        t++;
                      } else
                        t += 3;
                    }
                  } else {
                    float *col;
                    while(c--) {
                      if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
                                           || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 1))))
                                                                         &&
                                                                         (*
                                                                          (vi +
                                                                           (*(t + 2))))))
                      {

                        col = vc + (*t) * 3;
                        if(va) {
                          glColor4f(col[0], col[1], col[2], va[(*t)]);
                        } else {
                          glColor4f(col[0], col[1], col[2], alpha);
                        }
                        glNormal3fv(vn + (*t) * 3);
                        glVertex3fv(v + (*t) * 3);
                        t++;
                        col = vc + (*t) * 3;
                        if(va) {
                          glColor4f(col[0], col[1], col[2], va[(*t)]);
                        } else {
                          glColor4f(col[0], col[1], col[2], alpha);
                        }
                        glNormal3fv(vn + (*t) * 3);
                        glVertex3fv(v + (*t) * 3);
                        t++;
                        col = vc + (*t) * 3;
                        if(va) {
                          glColor4f(col[0], col[1], col[2], va[(*t)]);
                        } else {
                          glColor4f(col[0], col[1], col[2], alpha);
                        }
                        glNormal3fv(vn + (*t) * 3);
                        glVertex3fv(v + (*t) * 3);
                        t++;
                      } else
                        t += 3;
                    }
                  }
                  glEnd();
#endif
                }
              }
	      } /* end else use_shader */ 
            }
            /*          glDisable(GL_CULL_FACE);
               glDepthMask(GL_TRUE); */
          }
        } else if (ok) {                /* opaque */
          int simplify = 0; //, use_shader;
          simplify = (int) SettingGet(G, cSetting_simplify_display_lists);
            if(I->allVisibleFlag) {
              if(I->oneColorFlag) {
                if(use_display_lists && simplify) {      /* simplify: try to help display list optimizer */
		    if (generate_shader_cgo){
		      ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		      c = *(s++);
		      while(ok && c) {
			ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
			s += 2;
			while(ok && c--) {
			  s -= 2;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok)
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
		      }
		  } else {
		    glColor3fv(ColorGet(G, I->oneColor));
		    c = *(s++);
		    while(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		      {
			int nverts = c*3;
			int pl = 0;
			float *tmp_ptr;
			ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
			ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
			s += 2;
			while(c--) {
			  s -= 2;
			  tmp_ptr = vn + (*s) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			  tmp_ptr = v + (*s) * 3;
			  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			  s++;
			  tmp_ptr = vn + (*s) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			  tmp_ptr = v + (*s) * 3;
			  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			  pl += 3;
			  s++;
			  tmp_ptr = vn + (*s) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			  tmp_ptr = v + (*s) * 3;
			  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			  pl += 3;
			  s++;
			}		      
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, vertVals);
			glNormalPointer(GL_FLOAT, 0, normVals);
			glDrawArrays(GL_TRIANGLES, 0, nverts);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			DEALLOCATE_ARRAY(vertVals)
			DEALLOCATE_ARRAY(normVals)
		      }
#else
		      glBegin(GL_TRIANGLES);
		      s += 2;
		      while(c--) {
			s -= 2;
			glNormal3fv(vn + (*s) * 3);
			glVertex3fv(v + (*s) * 3);
			s++;
			glNormal3fv(vn + (*s) * 3);
			glVertex3fv(v + (*s) * 3);
			s++;
			glNormal3fv(vn + (*s) * 3);
			glVertex3fv(v + (*s) * 3);
			s++;
		      }
		      glEnd();
#endif
		      c = *(s++);
		    }
		  } /* end if else use_shader */
                } else if (ok) {
		    if (generate_shader_cgo){
		      CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		      c = *(s++);
		      while(ok && c) {
			if (ok) ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
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
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		  } else {
		    glColor3fv(ColorGet(G, I->oneColor));
		    c = *(s++);
		    while(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		      {
			int nverts = c + 2;
			int pl = 0;
			float *tmp_ptr;
			ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
			ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			while(c--) {
			  tmp_ptr = vn + (*s) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];		      
			  tmp_ptr = v + (*s) * 3;
			  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			  s++;
			}		      
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, vertVals);
			glNormalPointer(GL_FLOAT, 0, normVals);
			glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			DEALLOCATE_ARRAY(vertVals)
			DEALLOCATE_ARRAY(normVals)
		      }
#else
		      glBegin(GL_TRIANGLE_STRIP);
		      glNormal3fv(vn + (*s) * 3);
		      glVertex3fv(v + (*s) * 3);
		      s++;
		      glNormal3fv(vn + (*s) * 3);
		      glVertex3fv(v + (*s) * 3);
		      s++;
		      while(c--) {
			glNormal3fv(vn + (*s) * 3);
			glVertex3fv(v + (*s) * 3);
			s++;
		      }
		      glEnd();
#endif
		      c = *(s++);
		    }
		  }               /* use_display_lists&&simplify */
		} /* end if else generate_shader_cgo */
	      } else {          /* not one color */
		  if (generate_shader_cgo){
		    if(use_display_lists && simplify) {      /* simplify: try to help display list optimizer */
		      c = *(s++);
		      while(ok && c) {
			ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
			s += 2;
			while(ok && c--) {
			  s -= 2;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO);
			c = *(s++);
		      }
		    } else {
		      c = *(s++);
		      while(ok && c) {
			ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
			if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			if (ok && I->VAO){
			  ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			}
			if (ok) 
			  ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			s++;
			while(ok && c--) {
			  ok &= CGOColorv(I->shaderCGO, vc + (*s) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*s) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *s));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*s], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*s]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*s) * 3);
			  s++;
			}
			if (ok) ok &= CGOEnd(I->shaderCGO );
			c = *(s++);
		      }
		    }
		  } else {
                if(use_display_lists && simplify) {      /* simplify: try to help display list optimizer */
                  c = *(s++);
                  while(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		    {
		      int nverts = c*3;
		      int pl = 0, plc = 0;
		      float *tmp_ptr;
		      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
		      s += 2;
		      while(c--) {
			s -= 2;
			tmp_ptr = vc + (*s) * 3;
			colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			tmp_ptr = vc + (*s) * 3;
			colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
			tmp_ptr = vc + (*s) * 3;
			colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			tmp_ptr = v + (*s) * 3;
			vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
			s++;
		      }
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_NORMAL_ARRAY);
		      glEnableClientState(GL_COLOR_ARRAY);
		      glColorPointer(4, GL_FLOAT, 0, colorVals);
		      glNormalPointer(GL_FLOAT, 0, normVals);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glDrawArrays(GL_TRIANGLES, 0, nverts);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_NORMAL_ARRAY);
		      glDisableClientState(GL_COLOR_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		      DEALLOCATE_ARRAY(normVals)
		      DEALLOCATE_ARRAY(colorVals)
		    }
#else
                    glBegin(GL_TRIANGLES);
                    s += 2;
                    while(c--) {
                      s -= 2;
                      glColor3fv(vc + (*s) * 3);
                      glNormal3fv(vn + (*s) * 3);
                      glVertex3fv(v + (*s) * 3);
                      s++;
                      glColor3fv(vc + (*s) * 3);
                      glNormal3fv(vn + (*s) * 3);
                      glVertex3fv(v + (*s) * 3);
                      s++;
                      glColor3fv(vc + (*s) * 3);
                      glNormal3fv(vn + (*s) * 3);
                      glVertex3fv(v + (*s) * 3);
                      s++;
                    }
                    glEnd();
#endif
                    c = *(s++);
                  }
                } else {
                  c = *(s++);
                  while(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		    {
		      int nverts = c + 2;
		      int pl = 0, plc = 0;
		      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
		      float *tmp_ptr;
		      tmp_ptr = vc + (*s) * 3;
		      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1];
		      colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
		      tmp_ptr = vn + (*s) * 3;
		      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		      tmp_ptr = v + (*s) * 3;
		      vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		      s++; pl += 3;
		      tmp_ptr = vc + (*s) * 3;
		      colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1];
		      colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
		      tmp_ptr = vn + (*s) * 3;
		      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		      tmp_ptr = v + (*s) * 3;
		      vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
		      s++; pl += 3;
		      while(c--) {
			tmp_ptr = vc + (*s) * 3;
			colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1];
			colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = 1.f;
			tmp_ptr = vn + (*s) * 3;
			normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			tmp_ptr = v + (*s) * 3;
			vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
			pl += 3;
			s++;
		      }		      
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_NORMAL_ARRAY);
		      glEnableClientState(GL_COLOR_ARRAY);
		      glColorPointer(4, GL_FLOAT, 0, colorVals);
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glNormalPointer(GL_FLOAT, 0, normVals);
		      glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
		      glDisableClientState(GL_COLOR_ARRAY);
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_NORMAL_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		      DEALLOCATE_ARRAY(normVals)
		      DEALLOCATE_ARRAY(colorVals)
		    }

#else
                    glBegin(GL_TRIANGLE_STRIP);
                    glColor3fv(vc + (*s) * 3);
                    glNormal3fv(vn + (*s) * 3);
                    glVertex3fv(v + (*s) * 3);
                    s++;
                    glColor3fv(vc + (*s) * 3);
                    glNormal3fv(vn + (*s) * 3);
                    glVertex3fv(v + (*s) * 3);
                    s++;
                    while(c--) {
                      glColor3fv(vc + (*s) * 3);
                      glNormal3fv(vn + (*s) * 3);
                      glVertex3fv(v + (*s) * 3);
                      s++;
                    }
                    glEnd();
#endif
                    c = *(s++);
                  }
                }
	      } /* end if else generate_shader_cgo */
              }                 /* one color */
            } else if (ok) {            /* subsets */
		if (generate_shader_cgo){
		  c = I->NT;
		  if(ok && c) {
		    ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
		    if(I->oneColorFlag) {
		      if (ok) ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->oneColor));
		      while(ok && c--) {
			if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					     || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									   &&
									   (*
									    (vi + (*(t + 1))))
									   &&
									   (*
									    (vi +
									     (*(t + 2)))))) {
			  
			  ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			} else
			  t += 3;
		      }
		    } else {
		      while(ok && c--) {
			if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					     || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									   &&
									   (*
									    (vi + (*(t + 1))))
									   &&
									   (*
									    (vi +
									     (*(t + 2)))))) {
			  ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			  if (ok) ok &= CGOColorv(I->shaderCGO, vc + (*t) * 3);
			  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
			  if (ok && I->VAO){
			    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
			  }
			  if (ok) 
			    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
			  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
			  t++;
			} else
			  t += 3;
		      }
		    }
		    if (ok) ok &= CGOEnd(I->shaderCGO);
		  }
	      } else {  /* if else generate_shader_cgo */
              c = I->NT;
              if(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
		{
		  int nverts = 0, cinit = c, *tinit = t;
		  {
		    while(c--) {
		      if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					   || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									 &&
									 (*
									  (vi + (*(t + 1))))
									 &&
									 (*
									  (vi +
									   (*(t + 2)))))) {
			nverts += 3;
		      }
		      t += 3;
		    }
		    {
		      int pl = 0, plc = 0;
		      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
		      ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
		      float *tmp_ptr;
		      c = cinit;
		      t = tinit;
		      pl = plc = 0;
		      if(I->oneColorFlag) {
			glColor3fv(ColorGet(G, I->oneColor));
		      }
		      while(c--) {
			if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
					     || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
									   &&
									   (*
									    (vi + (*(t + 1))))
									   &&
									   (*
									    (vi +
									     (*(t + 2)))))) {
			  if(!I->oneColorFlag) {
			    tmp_ptr = vc + (*t) * 3;
			    colorVals[plc] = tmp_ptr[0]; colorVals[plc+1] = tmp_ptr[1]; colorVals[plc+2] = tmp_ptr[2];
			    colorVals[plc+3] = 1.;
			    plc+=4;
			  }
			  tmp_ptr = vn + (*t) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			  tmp_ptr = v + (*t) * 3;
			  vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
			  t++; pl += 3;
			  if(!I->oneColorFlag) {
			    tmp_ptr = vc + (*t) * 3;
			    colorVals[plc] = tmp_ptr[0]; colorVals[plc+1] = tmp_ptr[1]; colorVals[plc+2] = tmp_ptr[2];
			    colorVals[plc+3] = 1.;
			    plc+=4;
			  }
			  tmp_ptr = vn + (*t) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			  tmp_ptr = v + (*t) * 3;
			  vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
			  t++; pl += 3;
			  if(!I->oneColorFlag) {
			    tmp_ptr = vc + (*t) * 3;
			    colorVals[plc] = tmp_ptr[0]; colorVals[plc+1] = tmp_ptr[1]; colorVals[plc+2] = tmp_ptr[2];
			    colorVals[plc+3] = 1.;
			    plc+=4;
			  }
			  tmp_ptr = vn + (*t) * 3;
			  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
			  tmp_ptr = v + (*t) * 3;
			  vertVals[pl] = tmp_ptr[0]; vertVals[pl+1] = tmp_ptr[1]; vertVals[pl+2] = tmp_ptr[2];
			  t++; pl += 3;
			} else {
			  t += 3;
			}
		      }
		      glEnableClientState(GL_VERTEX_ARRAY);
		      glEnableClientState(GL_NORMAL_ARRAY);
		      if(!I->oneColorFlag) {
			glEnableClientState(GL_COLOR_ARRAY);
			glColorPointer(4, GL_FLOAT, 0, colorVals);
		      }
		      glVertexPointer(3, GL_FLOAT, 0, vertVals);
		      glNormalPointer(GL_FLOAT, 0, normVals);
		      glDrawArrays(GL_TRIANGLES, 0, nverts);
		      if(!I->oneColorFlag) {
			glDisableClientState(GL_COLOR_ARRAY);
		      }
		      glDisableClientState(GL_VERTEX_ARRAY);
		      glDisableClientState(GL_NORMAL_ARRAY);
		      DEALLOCATE_ARRAY(vertVals)
		      DEALLOCATE_ARRAY(normVals)
		      DEALLOCATE_ARRAY(colorVals)
		    }
		  }
		}
#else
                glBegin(GL_TRIANGLES);
                if(I->oneColorFlag) {
                  glColor3fv(ColorGet(G, I->oneColor));
                  while(c--) {
                    if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
                                         || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
                                                                       &&
                                                                       (*
                                                                        (vi + (*(t + 1))))
                                                                       &&
                                                                       (*
                                                                        (vi +
                                                                         (*(t + 2)))))) {

                      glNormal3fv(vn + (*t) * 3);
                      glVertex3fv(v + (*t) * 3);
                      t++;
                      glNormal3fv(vn + (*t) * 3);
                      glVertex3fv(v + (*t) * 3);
                      t++;
                      glNormal3fv(vn + (*t) * 3);
                      glVertex3fv(v + (*t) * 3);
                      t++;
                    } else
                      t += 3;
                  }
                } else {
                  while(c--) {
                    if((I->proximity && ((*(vi + (*t))) || (*(vi + (*(t + 1))))
                                         || (*(vi + (*(t + 2)))))) || ((*(vi + (*t)))
                                                                       &&
                                                                       (*
                                                                        (vi + (*(t + 1))))
                                                                       &&
                                                                       (*
                                                                        (vi +
                                                                         (*(t + 2)))))) {

                      glColor3fv(vc + (*t) * 3);
                      glNormal3fv(vn + (*t) * 3);
                      glVertex3fv(v + (*t) * 3);
                      t++;
                      glColor3fv(vc + (*t) * 3);
                      glNormal3fv(vn + (*t) * 3);
                      glVertex3fv(v + (*t) * 3);
                      t++;
                      glColor3fv(vc + (*t) * 3);
                      glNormal3fv(vn + (*t) * 3);
                      glVertex3fv(v + (*t) * 3);
                      t++;
                    } else
                      t += 3;
                  }
                }
                glEnd();
#endif
              }
	      } /* end if else use_shader */
            }
	  /*          if (use_shader) {
	    CShaderPrg_Disable(shaderPrg);
	    }*/
        }
      }
      if(ok && SettingGet(G, cSetting_surface_debug)) {
        t = I->T;
        c = I->NT;
	  if (generate_shader_cgo){
	    if(c) {
	      ok &= CGOBegin(I->shaderCGO, GL_TRIANGLES);
	      while(ok && c--) {
		if(I->allVisibleFlag
		   ||
		   ((I->proximity
		     && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		    || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {
		  ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
		  }
		  if (ok) 
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
		  }
		  if (ok) 
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok && I->VAO){
		    ok &= CGOAccessibility(I->shaderCGO, *(I->VAO + *t));
		  }
		  if (ok) 
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		} else {
		  t += 3;
		}
	      }
	      if (ok) ok &= CGOEnd(I->shaderCGO);
	    }
	} else {

        if(c) {
#ifdef _PYMOL_GL_DRAWARRAYS
	  {
	    int nverts = 0, cinit = c, *tinit = t;
	    while(c--) {
	      if(I->allVisibleFlag
		 ||
		 ((I->proximity
		   && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		  || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {
		nverts += 3;
	      }
	      t += 3;
	    }
	    {
	      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
	      ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
	      int pl = 0;
	      float *tmp_ptr;
	      c = cinit;
	      t = tinit;
	      while(c--) {
		if(I->allVisibleFlag
		   ||
		   ((I->proximity
		     && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		    || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {
		  tmp_ptr = vn + (*t) * 3;
		  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		  tmp_ptr = v + (*t) * 3;
		  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
		  t++;
		  tmp_ptr = vn + (*t) * 3;
		  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		  tmp_ptr = v + (*t) * 3;
		  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
		  t++;
		  tmp_ptr = vn + (*t) * 3;
		  normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
		  tmp_ptr = v + (*t) * 3;
		  vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
		  t++;
		} else {
		  t += 3;
		}
	      }
	      glEnableClientState(GL_VERTEX_ARRAY);
	      glEnableClientState(GL_NORMAL_ARRAY);
	      glVertexPointer(3, GL_FLOAT, 0, vertVals);
	      glNormalPointer(GL_FLOAT, 0, normVals);
	      glDrawArrays(GL_TRIANGLES, 0, nverts);
	      glDisableClientState(GL_VERTEX_ARRAY);
	      glDisableClientState(GL_NORMAL_ARRAY);
	      DEALLOCATE_ARRAY(vertVals)
	      DEALLOCATE_ARRAY(normVals)
	    }
	  }
#else
          glBegin(GL_TRIANGLES);
          while(c--) {

            if(I->allVisibleFlag
               ||
               ((I->proximity
                 && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
                || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {

              glNormal3fv(vn + (*t) * 3);
              glVertex3fv(v + (*t) * 3);
              t++;
              glNormal3fv(vn + (*t) * 3);
              glVertex3fv(v + (*t) * 3);
              t++;
              glNormal3fv(vn + (*t) * 3);
              glVertex3fv(v + (*t) * 3);
              t++;
            } else {
              t += 3;
            }
          }
          glEnd();
#endif
        }
	} /* end else use_shader */ 
        t = I->T;
        c = I->NT;

	  if (generate_shader_cgo){
	    if(ok && c) {
	      ok &= CGOColor(I->shaderCGO, 0.0, 1.0, 0.0);
	      if (ok) ok &= CGODotwidth(I->shaderCGO, 1.0F);
	      while(ok && c--) {
		ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
		if(I->allVisibleFlag
		   ||
		   ((I->proximity
		     && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		    || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok) 
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok) 
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		  if (ok) ok &= CGONormalv(I->shaderCGO, vn + (*t) * 3);
		  if (ok) 
		    ok &= CGOPickColor(I->shaderCGO, I->AT[*t], AtomInfoIsMasked((ObjectMolecule*)I->R.obj, I->AT[*t]));
		  if (ok) ok &= CGOVertexv(I->shaderCGO, v + (*t) * 3);
		  t++;
		} else {
		  t += 3;
		}
		if (ok) ok &= CGOEnd(I->shaderCGO);
	      }
	    }
	} else {/* end generate_shader_cgo */
        if(c) {
          glColor3f(0.0, 1.0, 0.0);
          glLineWidth(1.0F);
          while(c--) {
#ifdef _PYMOL_GL_DRAWARRAYS
            if(I->allVisibleFlag
               ||
               ((I->proximity
                 && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
                || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {
	      float *tmp_ptr;
	      int pl = 0;
	      ALLOCATE_ARRAY(GLfloat,vertVals,3*3)
	      ALLOCATE_ARRAY(GLfloat,normVals,3*3)
	      tmp_ptr = vn + (*t) * 3;
	      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = v + (*t) * 3;
	      vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
	      t++;
	      tmp_ptr = vn + (*t) * 3;
	      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = v + (*t) * 3;
	      vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
	      t++;
	      tmp_ptr = vn + (*t) * 3;
	      normVals[pl] = tmp_ptr[0]; normVals[pl+1] = tmp_ptr[1]; normVals[pl+2] = tmp_ptr[2];
	      tmp_ptr = v + (*t) * 3;
	      vertVals[pl++] = tmp_ptr[0]; vertVals[pl++] = tmp_ptr[1]; vertVals[pl++] = tmp_ptr[2];
	      t++;
	      glEnableClientState(GL_VERTEX_ARRAY);
	      glEnableClientState(GL_NORMAL_ARRAY);
	      glVertexPointer(3, GL_FLOAT, 0, vertVals);
	      glNormalPointer(GL_FLOAT, 0, normVals);
	      glDrawArrays(GL_LINE_STRIP, 0, 3);
	      glDisableClientState(GL_VERTEX_ARRAY);
	      glDisableClientState(GL_NORMAL_ARRAY);
	      DEALLOCATE_ARRAY(vertVals)
	      DEALLOCATE_ARRAY(normVals)
            } else {
              t += 3;
            }
#else
	    /* glBegin/glEnd should be inside the if statement, right? - BB */
            glBegin(GL_LINE_STRIP);

            if(I->allVisibleFlag
               ||
               ((I->proximity
                 && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
                || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2))))))) {

              glNormal3fv(vn + (*t) * 3);
              glVertex3fv(v + (*t) * 3);
              t++;
              glNormal3fv(vn + (*t) * 3);
              glVertex3fv(v + (*t) * 3);
              t++;
              glNormal3fv(vn + (*t) * 3);
              glVertex3fv(v + (*t) * 3);
              t++;
            } else {
              t += 3;
            }
            glEnd();
#endif
          }
        }
	} /* end else use_shader */ 

        c = I->N;
	if (generate_shader_cgo){
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
	} else {
        if(c) {
          SceneResetNormal(G, true);
          glColor3f(1.0, 0.0, 0.0);
#ifdef _PYMOL_GL_DRAWARRAYS
	  {
	    int nverts = c * 2, pl = 0;
	    ALLOCATE_ARRAY(GLfloat,lineVerts,nverts)
	    while(c--) {
	      lineVerts[pl++] = v[0]; lineVerts[pl++] = v[1]; lineVerts[pl++] = v[2];
	      lineVerts[pl++] = v[0] + vn[0] / 2;
	      lineVerts[pl++] = v[1] + vn[1] / 2;
	      lineVerts[pl++] = v[2] + vn[2] / 2;
	      v += 3;
	      vn += 3;
	    }
	    glEnableClientState(GL_VERTEX_ARRAY);
	    glVertexPointer(3, GL_FLOAT, 0, lineVerts);
	    glDrawArrays(GL_LINES, 0, nverts);
	    glDisableClientState(GL_VERTEX_ARRAY);
	    DEALLOCATE_ARRAY(lineVerts)
	  }
#else
          glBegin(GL_LINES);
          while(c--) {
            glVertex3fv(v);
            glVertex3f(v[0] + vn[0] / 2, v[1] + vn[1] / 2, v[2] + vn[2] / 2);
            v += 3;
            vn += 3;
          }
          glEnd();
#endif
        }
      } /* end else use_shader */ 
      }
      /* end of rendering, if using shaders, then render CGO */
      if (ok && use_shader) {
	CShaderPrg * shaderPrg = 0;
	if (generate_shader_cgo){
	  CGO *convertcgo = NULL;
	  if (ok) ok &= CGOStop(I->shaderCGO);
#ifdef _PYMOL_CGO_DRAWARRAYS


	  if (I->Type != 2){
	    if (ok)
	      convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);
	    CHECKOK(ok, convertcgo);
	    CGOFree(I->shaderCGO);
	    I->shaderCGO = convertcgo;
	    convertcgo = NULL;
	    I->pickingCGO = I->shaderCGO;
	    if (I->Type == 1){
	      /* Only needed to simplify spheres in surface_type=1 */
	      CGO *simple = CGOSimplify(I->shaderCGO, 0);
	      CHECKOK(ok, simple);
	      if (ok)
		I->pickingCGO = CGOCombineBeginEnd(simple ,0);
	      CHECKOK(ok, I->pickingCGO);
	      CGOFree(simple);
	    }
	  }
#else
	  (void)convertcgo;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
	  if(I->Type == 1){
	    if (dot_as_spheres) {
	      convertcgo = CGOOptimizeSpheresToVBONonIndexed(I->shaderCGO, CGO_BOUNDING_BOX_SZ + CGO_DRAW_SPHERE_BUFFERS_SZ);
	    } else {
	      convertcgo = CGOOptimizeToVBONotIndexedWithReturnedData(I->shaderCGO, 0, false, NULL);
	    }
	    CHECKOK(ok, convertcgo);
	    if (ok)
	      convertcgo->use_shader = true;
	  } else if (I->Type == 2) {
	    CGO *convertcgo2, *simple;
	    convertcgo2 = CGOConvertLinesToShaderCylinders(I->shaderCGO, 0);
	    CHECKOK(ok, convertcgo2);
	    if (ok)
	      convertcgo = CGOOptimizeGLSLCylindersToVBOIndexed(convertcgo2, 0);
	    CHECKOK(ok, convertcgo);
	    if (ok)
	      convertcgo->use_shader = true;
	    if (ok)
	      simple = CGOCombineBeginEnd(convertcgo2, 0);
	    CHECKOK(ok, simple);
	    if (ok)
	      I->pickingCGO = CGOSimplify(simple, 0);
	    CHECKOK(ok, I->pickingCGO);
	    CGOFree(simple);
	    CGOFree(convertcgo2);
	  } else {
	    if(I->Type != 1 && ((alpha != 1.0) || va)) {
	      if (ok)
		convertcgo = CGOOptimizeToVBOIndexedNoShader(I->shaderCGO, 0);
	      CHECKOK(ok, convertcgo);
	    } else {
	      if (ok)
		convertcgo = CGOOptimizeToVBONotIndexedWithReturnedData(I->shaderCGO, 0, false, NULL);
	      CHECKOK(ok, convertcgo);
	    }
	    if (ok)
	      convertcgo->use_shader = true;
	  }	  
	  if(ok && I->Type == 1 && !dot_as_spheres) {
	    I->shaderCGO = CGONew(G);
	    CHECKOK(ok, I->shaderCGO);
	    if (ok){
	      I->shaderCGO->use_shader = true;
	      if (ok) ok &= CGOResetNormal(I->shaderCGO, true);
	      if (ok) ok &= CGOLinewidthSpecial(I->shaderCGO, POINTSIZE_DYNAMIC_DOT_WIDTH);
	      if (ok) ok &= CGOStop(I->shaderCGO);
	      if (ok) CGOAppend(I->shaderCGO, convertcgo);
	      CGOFreeWithoutVBOs(convertcgo);
	      convertcgo = NULL;
	    }
	  } else {
	    I->shaderCGO = convertcgo;
	    convertcgo = NULL;
	  }
#else
	  (void)convertcgo;
#endif
	}

	if (I->Type == 1){
	  if (dot_as_spheres){
	    float radius = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_width)
                           * info->vertex_scale;
	    shaderPrg = CShaderPrg_Enable_DefaultSphereShader(G);
	    CShaderPrg_Set1f(shaderPrg, "sphere_size_scale", fabs(radius));
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    SceneResetNormalUseShaderAttribute(G, 0, true, CShaderPrg_GetAttribLocation(shaderPrg, "a_Normal"));
	  }
	} else if (I->Type == 2) {
	  float mesh_width = SettingGet_f(G, I->R.obj->Setting, NULL, cSetting_mesh_width);
	  shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	  CShaderPrg_Set1f(shaderPrg, "uni_radius", SceneGetLineWidthForCylinders(G, info, mesh_width));
	} else {
	  shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	  CShaderPrg_Set1f(shaderPrg, "ambient_occlusion_scale", ambient_occlusion_scale);
	  CShaderPrg_Set1i(shaderPrg, "use_interior_color_threshold", 1);
	  {
	    if((alpha != 1.0) || va) {
	      int t_mode =
		SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting,
			     cSetting_transparency_mode);
	      if(info && info->alpha_cgo) {
		t_mode = 0;
	      }
	      if (t_mode){
		RepSurfaceSortIX(G, I, t_mode);
		if (I->ix){
		  float *pc = CGOGetNextDrawBufferedIndex(I->shaderCGO->op);
		  if (pc){
		    int nindices = CGO_get_int(pc+3), c, pl, idx;
		    uint vbuf = CGO_get_int(pc+8);
		    uint *vertexIndices;
		    vertexIndices = Alloc(uint, nindices);	      
		    if (!vertexIndices){
		      PRINTFB(I->R.G, FB_RepSurface, FB_Errors) "ERROR: RepSurfaceRender() vertexIndices could not be allocated: nindices=%d\n", nindices ENDFB(I->R.G);
		    }
		    for(c = 0, pl=0; c < I->n_tri; c++) {
		      idx = I->ix[c] * 3;
		      vertexIndices[pl++] = idx;
		      vertexIndices[pl++] = idx + 1;
		      vertexIndices[pl++] = idx + 2;
		    }
		    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbuf);
		    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint)*nindices, vertexIndices, GL_STATIC_DRAW);
		    I->vertexIndices = vertexIndices;
		  }
		}
	      }
	    }
	  }
	}
    if (I->shaderCGO){
        I->shaderCGO->enable_shaders = shaderPrg ? 0 : 1;
        CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
        if (shaderPrg){
            CShaderPrg_Disable(shaderPrg);
        }
    }
      }
#ifdef _PYMOL_GL_CALLLISTS
    if (use_display_lists && I->R.displayList){
      glEndList();
      glCallList(I->R.displayList);      
    }
#endif
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
    I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
    I->R.cs->Active[cRepSurface] = false;
  }
}

int RepSurfaceSameVis(RepSurface * I, CoordSet * cs)
{
  int same = true;
  char *lv;
  int a;
  AtomInfoType *ai;

  ai = cs->Obj->AtomInfo;
  lv = I->LastVisib;

  for(a = 0; a < cs->NIndex; a++) {
    if(*(lv++) != (ai + cs->IdxToAtm[a])->visRep[cRepSurface]) {
      same = false;
      break;
    }
  }
  return (same);
}

void RepSurfaceInvalidate(struct Rep *I, struct CoordSet *cs, int level){
  RepInvalidate(I, cs, level);
  if (level >= cRepInvColor)
    ((RepSurface*)I)->ColorInvalidated = true;
}

int RepSurfaceSameColor(RepSurface * I, CoordSet * cs)
{
  int same = true;
  int *lc, *cc;
  int a;
  AtomInfoType *ai;

  same = !I->ColorInvalidated;
  if (same){
    ai = cs->Obj->AtomInfo;
    lc = I->LastColor;
    cc = cs->Color;
    for(a = 0; a < cs->NIndex; a++) {
      if((ai + cs->IdxToAtm[a])->visRep[cRepSurface]) {
	if(*(lc++) != *(cc++)) {
	  same = false;
	  break;
	}
      }
    }
  }
  return (same);
}

void RepSurfaceColor(RepSurface * I, CoordSet * cs)
{
  PyMOLGlobals *G = cs->State.G;
  MapType *map = NULL, *ambient_occlusion_map = NULL;
  int a, i0, i, j, c1;
  float *v0, *vc, *c0, *va;
  float *n0;
  int *vi, *lc, *cc;
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
  char *carve_selection = NULL;
  float *carve_vla = NULL;
  MapType *carve_map = NULL;

  int clear_state = 0;
  int clear_flag = false;
  float clear_cutoff;
  char *clear_selection = NULL;
  float *clear_vla = NULL;
  int state = I->R.context.state;
  float transp;
  int variable_alpha = false;

  MapType *clear_map = NULL;

  AtomInfoType *ai2 = NULL, *ai1;

  obj = cs->Obj;
  ambient_occlusion_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ambient_occlusion_mode);
  surface_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_mode);
  ramp_above =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_ramp_above_mode);
  surface_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_surface_color);
  cullByFlag = (surface_mode == cRepSurface_by_flags);
  inclH = !((surface_mode == cRepSurface_heavy_atoms)
            || (surface_mode == cRepSurface_vis_heavy_only));
  inclInvis = !((surface_mode == cRepSurface_vis_only)
                || (surface_mode == cRepSurface_vis_heavy_only));
  probe_radius = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_solvent_radius);
  I->proximity =
    SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_surface_proximity);
  carve_cutoff =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_carve_cutoff);
  clear_cutoff =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_clear_cutoff);
  transp = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_transparency);
  carve_normal_cutoff =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_carve_normal_cutoff);
  carve_normal_flag = carve_normal_cutoff > (-1.0F);

  cutoff = I->max_vdw + 2 * probe_radius;

  if(!I->LastVisib)
    I->LastVisib = Alloc(char, cs->NIndex);
  if(!I->LastColor)
    I->LastColor = Alloc(int, cs->NIndex);
  lv = I->LastVisib;
  lc = I->LastColor;
  cc = cs->Color;
  ai2 = obj->AtomInfo;
  for(a = 0; a < cs->NIndex; a++) {
    *(lv++) = (ai2 + cs->IdxToAtm[a])->visRep[cRepSurface] ? 1 : 0;
    *(lc++) = *(cc++);
  }

  if(I->N) {
    if(carve_cutoff > 0.0F) {
      carve_state =
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_carve_state) - 1;
      carve_selection =
        SettingGet_s(G, cs->Setting, obj->Obj.Setting, cSetting_surface_carve_selection);
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
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_clear_state) - 1;
      clear_selection =
        SettingGet_s(G, cs->Setting, obj->Obj.Setting, cSetting_surface_clear_selection);
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
      I->VC = Alloc(float, 3 * I->N);
    vc = I->VC;
    if(!I->VA)
      I->VA = Alloc(float, I->N);
    va = I->VA;
    if(!I->RC)
      I->RC = Alloc(int, I->N);
    rc = I->RC;
    if(!I->Vis)
      I->Vis = Alloc(int, I->N);
    vi = I->Vis;
    if(ColorCheckRamped(G, surface_color)) {
      I->oneColorFlag = false;
    } else {
      I->oneColorFlag = true;
    }
    first_color = -1;

    present = Alloc(int, cs->NIndex);
    {
      int *ap = present;
      for(a = 0; a < cs->NIndex; a++) {
        ai1 = obj->AtomInfo + cs->IdxToAtm[a];
        if(ai1->visRep[cRepSurface] &&
           (inclH || (!ai1->hydrogen)) &&
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
            v0 = cs->Coord + 3 * a;
            i = *(MapLocusEStart(map, v0));
            if(i && map->EList) {
              j = map->EList[i++];
              while(j >= 0) {
                if(present[j] > 1) {
                  ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                  if(within3f
                     (cs->Coord + 3 * j, v0, ai1->vdw + ai2->vdw + probe_radiusX2)) {
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

    if (ambient_occlusion_mode){
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
	float *VAO = Alloc(float, cs->NIndex);
	short *nVAO = Alloc(short, cs->NIndex);
	memset(VAO, 0, sizeof(float)*cs->NIndex);
	memset(nVAO, 0, sizeof(short)*cs->NIndex);

	for(a = 0; a < I->N; a++) {
	  int nbits = 0;
	  short level1, level2, has;
	  unsigned long bits = 0L, bit;
	  float d[3], *vn0, planew, v0mod[3];
	  int closeA = -1;
	  float closeDist = MAXFLOAT;
	  has = 0;
	  
	  v0 = I->V + 3 * a;
	  vn0 = I->VN + 3 * a;
	  mult3f(vn0, .01f, v0mod);
	  add3f(v0, v0mod, v0mod);
	  planew = -(vn0[0] * v0mod[0] + vn0[1] * v0mod[1] + vn0[2] * v0mod[2]);

	  i = *(MapLocusEStart(ambient_occlusion_map, v0));
	  if(i && map->EList) {
	    j = ambient_occlusion_map->EList[i++];
	    while(j >= 0) {
	      subtract3f(cs->Coord + j * 3, v0, d);
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
	      v0 = cs->Coord + 3 * closeA;
	      i = *(MapLocusEStart(ambient_occlusion_map, v0));
	      if (i){
		j = ambient_occlusion_map->EList[i++];
		while(j >= 0) {
		  if (closeA==j){
		    j = ambient_occlusion_map->EList[i++];
		    continue;
		  }
		  subtract3f(cs->Coord + j * 3, v0, d);
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
	  float d[3], *vn0, pt[3], planew, v0mod[3];
	  
	  if (a%1000==0){
	    PRINTFB(I->R.G, FB_RepSurface, FB_Debugging) "RepSurfaceColor():  Ambient Occlusion computing mode=%d #vertices=%d done=%d\n", ambient_occlusion_mode, I->N, a ENDFB(I->R.G);      
	  }
	  v0 = I->V + 3 * a;
	  vn0 = I->VN + 3 * a;
	  mult3f(vn0, .01f, v0mod);
	  add3f(v0, v0mod, v0mod);
	  planew = -(vn0[0] * v0mod[0] + vn0[1] * v0mod[1] + vn0[2] * v0mod[2]);
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
		copy3f(cs->Coord + j * 3, pt);
		subtract3f(cs->Coord + j * 3, v0mod, d);
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
	PRINTFB(I->R.G, FB_RepSurface, FB_Debugging) "RepSurfaceColor():  #vertices=%d Ambient Occlusion average #atoms looked at per vertex = %f\n", I->N, (natomsL/I->N) ENDFB(I->R.G);
      }
      {
	int ambient_occlusion_smooth = 
	  SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ambient_occlusion_smooth);
	int surface_quality = 
	  SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_quality);
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
	  float *tmpVAO = Alloc(float, I->N);
	  int *nVAO = Alloc(int, I->N), c, *t;
	  
	  for (j=0; j<ambient_occlusion_smooth; j++){
	    memset(nVAO, 0, sizeof(int)*I->N);
	    memset(tmpVAO, 0, sizeof(float)*I->N);

	    t = I->T;
	    c = I->NT;
	    while (c--){
	      if((I->proximity
		  && ((*(vi + (*t))) || (*(vi + (*(t + 1)))) || (*(vi + (*(t + 2))))))
		 || ((*(vi + (*t))) && (*(vi + (*(t + 1)))) && (*(vi + (*(t + 2)))))) {
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

      PRINTFB(I->R.G, FB_RepSurface, FB_Debugging) "RepSurfaceColor():  Ambient Occlusion computed #atoms=%d #vertices=%d time=%lf seconds\n", cs->NIndex, I->N, (cur_time-start_time) ENDFB(I->R.G);      
      ambient_occlusion_map = NULL;
    } else {
      if (I->VAO){
	VLAFreeP(I->VAO);
	I->VAO = 0;
      }
    }
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
        float minDist = MAXFLOAT, minDist2 = MAXFLOAT, distDiff = MAXFLOAT;
	int pi = -1, pi2 = -1, catm = -1; /* variables for color smoothing */
        AtomInfoType *pai = NULL, *pai2 = NULL; /* variables for color smoothing */
        c1 = 1;
        i0 = -1;
        v0 = I->V + 3 * a;
        n0 = I->VN + 3 * a;
        vi = I->Vis + a;
        /* colors */
        i = *(MapLocusEStart(map, v0));
        if(i && map->EList) {
          j = map->EList[i++];
          while(j >= 0) {
	    atm = cs->IdxToAtm[j];
            ai2 = obj->AtomInfo + atm;
            if((inclH || (!ai2->hydrogen)) &&
               ((!cullByFlag) || (!(ai2->flags & cAtomFlag_ignore)))) {
              dist = (float) diff3f(v0, cs->Coord + j * 3) - ai2->vdw;
              if(color_smoothing){
		if (dist < minDist){
		  /* switching closest to 2nd closest */
		  pi2 = pi;
		  pai2 = pai;
		  minDist2 = minDist;
		  pi = j;
		  catm = atm;
		  pai = ai2;
		  minDist = dist;
		} else if (dist < minDist2){
		  /* just setting second closest */
		  pi2 = j;
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
                   atoms points (cs->Coord + pi * 3) and (cs->Coord + pi2 * 3) (including vdw), then set this
                   distance to the distance between the vertex v0
                   and that point. We might want to use the normal
                   to compute this distance.
           */
          distDiff = fabs(minDist2-minDist);
        }
        if(i0 >= 0) {
          int at_surface_color;
          transp = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_transparency);

          AtomInfoGetSetting_f(G, ai0, cSetting_transparency, transp, &at_transp);

          AtomInfoGetSetting_color(G, ai0, cSetting_surface_color,
                                   surface_color, &at_surface_color);

          if(at_surface_color != -1) {
            c1 = at_surface_color;
            distDiff = MAXFLOAT;
          } else {
            c1 = *(cs->Color + i0);
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
            if(ai2->visRep[cRepSurface] &&
               (inclH || (!ai2->hydrogen)) &&
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
                register float *v_targ = carve_vla + 3 * j;
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
          if (color_smoothing && distDiff < color_smoothing_threshold){
            float *c2;
            float weight, weight2;
            if (color_smoothing==1){
              weight = 1.f + sin(.5f * PI * (distDiff / color_smoothing_threshold));
            } else {
              weight = 1.f + (distDiff / color_smoothing_threshold);
            }
            weight2 = 2.f - weight;
            c0 = ColorGet(G, c1);
            c2 = ColorGet(G, *(cs->Color + pi2));
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

        if(!*vi)
          I->allVisibleFlag = false;
        vi++;
      }
      MapFree(map);
    }
    if(variable_alpha)
      I->oneColorFlag = false;

    if(I->oneColorFlag) {
      I->oneColor = first_color;
    }
  }
  /*
     if(surface_color>=0) {
     I->oneColorFlag=true;
     I->oneColor=surface_color;
     }
   */

  if(G->HaveGUI) {
    if (I->shaderCGO){
      CGOFree(I->shaderCGO);
      I->shaderCGO = NULL;
    }
    if (I->vertexIndices){
      FreeP(I->vertexIndices);
      I->vertexIndices = 0;
    }
    if (I->sum){
      FreeP(I->sum);
      I->sum = 0;
    }
    if (I->z_value){
      FreeP(I->z_value);
      I->z_value = 0;
    }
    if (I->ix){
      FreeP(I->ix);
      I->ix = 0;
    }
    I->n_tri = 0;
#ifdef _PYMOL_GL_CALLLISTS
    if(I->R.displayList) {
      if(PIsGlutThread()) {
        if(G->ValidContext) {
          glDeleteLists(I->R.displayList, 1);
          I->R.displayList = 0;
        }
      } else {
        char buffer[255];       /* pass this off to the main thread */
        sprintf(buffer, "_cmd.gl_delete_lists(cmd._COb,%d,%d)\n", I->R.displayList, 1);
        PParse(G, buffer);
        I->R.displayList = 0;
      }
    }
#endif
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
     || (!SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_ray_color_ramps)))
    FreeP(I->RC);
  I->ColorInvalidated = false;
  FreeP(present);
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

static int SurfaceJobRun(PyMOLGlobals * G, SurfaceJob * I)
{
  int ok = true;
  int MaxN;
  int n_index = VLAGetSize(I->atomInfo);
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
    float *I_coord = I->coord;
    int *present_vla = I->presentVla;
    SurfaceJobAtomInfo *I_atom_info = I->atomInfo;

    I->N = 0;

    sol_dot = SolventDotNew(G, I->coord, I->atomInfo, probe_radius,
                            ssp, present_vla,
                            circumscribe, I->surfaceMode, I->surfaceSolvent,
                            I->cavityCull, I->allVisibleFlag, I->maxVdw,
                            I->cavityMode, I->cavityRadius, I->cavityCutoff);
    CHECKOK(ok, sol_dot);
    ok &= !G->Interrupt;
    if(ok && sol_dot) {
      if(!I->surfaceSolvent) {

        float solv_tole = point_sep * 0.04F;
        float probe_rad_more;
        float probe_rad_less;
        float probe_rad_less2;

        if(probe_radius < (2.5F * point_sep)) { /* minimum probe radius allowed */
          probe_radius = 2.5F * point_sep;
        }

        probe_rad_more = probe_radius * (1.0F + solv_tole);

        switch (surface_type) {
        case 0:                /* solid */
        case 3:
        case 4:
        case 5:
        case 6:
          probe_rad_less = probe_radius;
          break;
        default:
          probe_rad_less = probe_radius * (1.0F - solv_tole);
          break;
        }
        probe_rad_less2 = probe_rad_less * probe_rad_less;

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
          MapType *map, *solv_map;
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
              int *dc = sol_dot->dotCode;
              Vector3f *dot = Alloc(Vector3f, sp->nDot);
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
                float dist2 = probe_rad_less2;
                int sp_nDot = sp->nDot;
                for(a = 0; ok && a < sol_dot->nDot; a++) {
                  if(dc[a] || (surface_type < 6)) {     /* surface type 6 is completely scribed */
                    OrthoBusyFast(G, a + sol_dot->nDot * 2, sol_dot->nDot * 5); /* 2/5 to 3/5 */
                    for(b = 0; ok && b < sp_nDot; b++) {
                      float *dot_b = dot[b];
                      v[0] = v0[0] + dot_b[0];
                      v[1] = v0[1] + dot_b[1];
                      v[2] = v0[2] + dot_b[2];
                      {
                        int flag = true;
                        int ii;
			ii = *(MapLocusEStart(solv_map, v));
                        if(ii && solv_map->EList) {
                          register float *i_dot = sol_dot->dot;
                          register float dist = probe_rad_less;
                          register int *elist_ii = solv_map->EList + ii;
                          register float v_0 = v[0];
                          register int jj_next, jj = *(elist_ii++);
                          register float v_1 = v[1];
                          register float *v1 = i_dot + 3 * jj;
                          register float v_2 = v[2];
                          while(jj >= 0) {
                            /* huge bottleneck -- optimized for superscaler processors */
                            register float dx = v1[0], dy, dz;
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
                                        flag = false;
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

                        /* at this point, we have points on the interior of the solvent surface,
                           so now we need to further trim that surface to cover atoms that are present */

                        if(flag) {
                          register int i = *(MapLocusEStart(map, v));
                          if(i && map->EList) {
                            register int j = map->EList[i++];
                            while(j >= 0) {
                              SurfaceJobAtomInfo *atom_info = I_atom_info + j;
                              if((!present_vla) || present_vla[j]) {
                                if(within3f
                                   (I_coord + 3 * j, v,
                                    atom_info->vdw + probe_rad_more)) {
                                  flag = false;
                                  break;
                                }
                              }
                              j = map->EList[i++];
                            }
                          }
                          if(!flag) {   /* compute the normals */
                            vn[0] = -sp->dot[b][0];
                            vn[1] = -sp->dot[b][1];
                            vn[2] = -sp->dot[b][2];
                            if(I->N < MaxN) {
                              I->N++;
                              v += 3;
                              vn += 3;
                            } else {
                              int v_offset = v - I->V;
                              int vn_offset = vn - I->VN;
                              MaxN = MaxN * 2;
                              VLASize(I->V, float, (MaxN + 1) * 3);
			      CHECKOK(ok, I->V);
			      if (ok)
				VLASize(I->VN, float, (MaxN + 1) * 3);
			      CHECKOK(ok, I->VN);
			      if (ok){
				v = I->V + v_offset;
				vn = I->VN + vn_offset;
			      }
                            }
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
          int n_new = 0;

          float neighborhood = 2.6 * point_sep; /* these constants need more tuning... */
          float dot_cutoff = 0.666;
          float insert_cutoff = 1.1 * point_sep;

          float map_cutoff = neighborhood;
          float *new_dot = VLAlloc(float, 1000);

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
              register int i = *(MapLocusEStart(map, v));
              if(i && map->EList) {
                register int j = map->EList[i++];
                while(ok && j >= 0) {
                  if(j > a) {
                    float *v0 = I->V + 3 * j;
                    if(within3f(v0, v, map_cutoff)) {
                      int add_new = false;
                      float *n0 = I->VN + 3 * j;
                      VLACheck(new_dot, float, n_new * 6 + 5);
		      CHECKOK(ok, new_dot);
                      if (ok){
                        float *v1 = new_dot + n_new * 6;
                        average3f(v, v0, v1);
                        if((dot_product3f(n0, vn) < dot_cutoff)
                           && (within3f(v0, v, neighborhood)))
                          add_new = true;
                        else {
                          /* if points are too far apart, insert a new one */
                          register int ii = *(MapLocusEStart(map, v1));
                          if(ii) {
                            int found = false;
                            register int jj = map->EList[ii++];
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
                          /* highly divergent */
                          float *n1 = v1 + 3;
                          n_new++;
                          average3f(vn, n0, n1);
                          normalize3f(n1);
                        }
                      }
                    }
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
            float *n1 = new_dot + 3;
            float *v1 = new_dot;
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
          }
          VLAFreeP(new_dot);
        }

        if(ok && I->N && (surface_type == 0) && (circumscribe)) {

          float cutoff = 0.5 * probe_radius;

          /* combine scribing with an atom proximity cleanup pass */

          int *dot_flag = Calloc(int, I->N);
          MapType *map =
            MapNewFlagged(G, I->maxVdw + probe_radius, I_coord, n_index, NULL,
                          present_vla);
          int a;
	  CHECKOK(ok, map);
          if (ok)
	    ok &= MapSetupExpress(map);
          v = I->V;
          for(a = 0; ok && a < I->N; a++) {
            register int i = *(MapLocusEStart(map, v));
            if(i && map->EList) {
              register int j = map->EList[i++];
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
          map = NULL;

          if(ok) {
            /* purge unused dots */

            float *v0 = I->V;
            float *vn0 = I->VN;
            int *p = dot_flag;
            int c = I->N;
            int a;
            v = I->V;
            vn = I->VN;
            I->N = 0;
            for(a = 0; a < c; a++) {
              if(*(p++)) {
                *(v0++) = *(v++);
                *(v0++) = *(v++);
                *(v0++) = *(v++);
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
          FreeP(dot_flag);
        }

        /* now, eliminate dots that are too close to each other */

        /*    CGOColor(I->debug,0.0,1.0,0.0);
           CGOBegin(I->debug,GL_POINTS);
           for(a=0;a<I->N;a++)
           CGOVertexv(I->debug,I->V+3*a);
           CGOEnd(I->debug);
         */

        if(ok && I->N) {
          int repeat_flag = true;
          float min_dot = 0.1F;
          int *dot_flag = Alloc(int, I->N);
	  CHECKOK(ok, dot_flag);
          while(ok && repeat_flag) {
            repeat_flag = false;

            if(surface_type >= 3) {
              register int jj;
              float dist;
              register float nearest;
              register float min_sep2 = point_sep * point_sep;
              float diff[3];
              {
                int a;
                for(a = 0; a < I->N; a++)
                  dot_flag[a] = 1;
              }
              {
                MapType *map = MapNew(G, point_sep + 0.05F, I->V, I->N, NULL);
                int a;
		CHECKOK(ok, map);
		if (ok)
		  ok &= MapSetupExpress(map);
                v = I->V;
                vn = I->VN;
                for(a = 0; ok && a < I->N; a++) {
                  if(dot_flag[a]) {
                    register int i = *(MapLocusEStart(map, v));
                    if(i && map->EList) {
                      register int j = map->EList[i++];
                      jj = I->N;
                      nearest = point_sep + 1.0F;
                      while(j >= 0) {
                        if(j > a) {
                          if(dot_flag[j]) {
                            if(dot_product3f(I->VN + (3 * j), vn) > min_dot) {
                              if(within3fret
                                 (I->V + (3 * j), v, point_sep, min_sep2, diff, &dist)) {
                                repeat_flag = true;
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
                        repeat_flag = true;
                      }
                    }
                  }
                  v += 3;
                  vn += 3;
		  ok &= !G->Interrupt;
                }
                MapFree(map);
              }
            } else {            /* surface types < 3 */
              int a;
              MapType *map = MapNew(G, -point_sep, I->V, I->N, NULL);
	      CHECKOK(ok, map);
	      if (ok){
		for(a = 0; a < I->N; a++)
		  dot_flag[a] = 1;
		ok &= MapSetupExpress(map);
	      }
              v = I->V;
              vn = I->VN;
              for(a = 0; ok && a < I->N; a++) {
                if(dot_flag[a]) {
                  register int i = *(MapLocusEStart(map, v));
                  if(i && map->EList) {
                    register int j = map->EList[i++];
                    while(j >= 0) {
                      if(j != a) {
                        if(dot_flag[j]) {
                          if(within3f(I->V + (3 * j), v, point_sep)) {
                            dot_flag[j] = 0;
                            add3f(vn, I->VN + (3 * j), vn);
                            average3f(I->V + (3 * j), v, v);
                            repeat_flag = true;
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
            }

            if(ok) {
              float *v0 = I->V;
              float *vn0 = I->VN;
              int *p = dot_flag;
              int c = I->N;
              int a;
              v = I->V;
              vn = I->VN;
              I->N = 0;
              for(a = 0; a < c; a++) {
                if(*(p++)) {
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
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
	    ok &= !G->Interrupt;
          }
          FreeP(dot_flag);
        }

        /* now eliminate troublesome vertices in regions of extremely high curvature */

        if(ok && (surface_type != 3) &&
           I->N && (I->trimCutoff > 0.0F) && (I->trimFactor > 0.0F)) {
          float trim_cutoff = I->trimCutoff;
          float trim_factor = I->trimFactor;
          int repeat_flag = true;
          float neighborhood = trim_factor * point_sep;
          float dot_sum;
          int n_nbr;
          int *dot_flag = Alloc(int, I->N);
	  CHECKOK(ok, dot_flag);
          if(ok && surface_type == 6) {       /* emprical tweaks */
            trim_factor *= 2.5;
            trim_cutoff *= 1.5;
          }
          while(ok && repeat_flag) {
            int a;
            MapType *map = MapNew(G, neighborhood, I->V, I->N, NULL);
	    CHECKOK(ok, map);
	    if (ok){
	      repeat_flag = false;
	      for(a = 0; a < I->N; a++)
		dot_flag[a] = 1;
	      ok &= MapSetupExpress(map);
	    }
            v = I->V;
            vn = I->VN;
            for(a = 0; ok && a < I->N; a++) {
              if(dot_flag[a]) {
                register int i = *(MapLocusEStart(map, v));
                if(i && map->EList) {
                  register int j = map->EList[i++];
                  n_nbr = 0;
                  dot_sum = 0.0F;
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
                      repeat_flag = true;
                    }
                  }
                }
              }
              v += 3;
              vn += 3;
	      ok &= !G->Interrupt;
            }

            if(ok) {
              float *v0 = I->V;
              float *vn0 = I->VN;
              int *p = dot_flag;
              int c = I->N;
              v = I->V;
              vn = I->VN;
              I->N = 0;
              for(a = 0; a < c; a++) {
                if(*(p++)) {
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
                  *(v0++) = *(v++);
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
            MapFree(map);
	    ok &= !G->Interrupt;
          }
          FreeP(dot_flag);
        }
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

Rep *RepSurfaceNew(CoordSet * cs, int state)
{
  int ok = true;
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj = cs->Obj;
  OOCalloc(G, RepSurface);
  CHECKOK(ok, I);
  if (!ok)
    return NULL;
  I->pickingCGO = I->shaderCGO = 0;
  I->vertexIndices = 0;
  I->sum = 0;
  I->z_value = 0;
  I->ix = 0;
  I->n_tri = 0;
  I->AT = 0;
  I->ColorInvalidated = false;
  {
    int surface_mode =
      SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_mode);
    int cullByFlag = (surface_mode == cRepSurface_by_flags);
    int inclH = !((surface_mode == cRepSurface_heavy_atoms)
                  || (surface_mode == cRepSurface_vis_heavy_only));
    int inclInvis = !((surface_mode == cRepSurface_vis_only)
                      || (surface_mode == cRepSurface_vis_heavy_only));
    int visFlag = false;

    if(obj->RepVisCache[cRepSurface]) {
      register int *idx_to_atm = cs->IdxToAtm;
      register AtomInfoType *obj_AtomInfo = obj->AtomInfo;
      register int a, cs_NIndex = cs->NIndex;
      for(a = 0; a < cs_NIndex; a++) {
        register AtomInfoType *ai1 = obj_AtomInfo + *(idx_to_atm++);
        if(ai1->visRep[cRepSurface] &&
           (inclH || (!ai1->hydrogen)) &&
           ((!cullByFlag) || (!(ai1->flags &
                                (cAtomFlag_exfoliate | cAtomFlag_ignore))))) {
          visFlag = true;
          break;
        }
      }
    }
    if(!visFlag) {
      OOFreeP(I);
      return (NULL);            /* skip if no thing visible */
    }

    {
      int surface_flag = false;

      int surface_type =
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_type);
      int surface_solvent =
        SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_surface_solvent);
      int surface_quality =
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_quality);
      float probe_radius =
        SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_solvent_radius);
      int optimize =
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_optimize_subsets);

      int circumscribe =
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_circumscribe);
      float trim_cutoff =
        SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_trim_cutoff);

      float trim_factor =
        SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_trim_factor);
      int sphere_idx = 0, solv_sph_idx = 0;
      MapType *map = NULL;
      float point_sep;
      int *present_vla = NULL;
      int n_present = 0;
      int carve_state = 0;
      int carve_flag = false;
      float carve_cutoff;
      char *carve_selection = NULL;
      float *carve_vla = NULL;
      MapType *carve_map = NULL;

      int cavity_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_surface_cavity_mode);
      float cavity_radius = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_cavity_radius);
      float cavity_cutoff = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_cavity_cutoff);

      I->Type = surface_type;

      I->max_vdw = ObjectMoleculeGetMaxVDW(obj);

      if(surface_quality >= 4) {        /* totally impractical */
        point_sep =
          SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_best) / 4;
        sphere_idx = 4;
        solv_sph_idx = 4;
        if(circumscribe < 0)
          circumscribe = 91;
      } else {
        switch (surface_quality) {
        case 3:                /* nearly impractical */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_best) / 3;
          sphere_idx = 4;
          solv_sph_idx = 3;
          if(circumscribe < 0)
            circumscribe = 71;
          break;
        case 2:
          /* nearly perfect */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_best) / 2;
          sphere_idx = 3;
          solv_sph_idx = 3;
          if(circumscribe < 0)
            circumscribe = 41;
          break;
        case 1:
          /* good */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_best);
          sphere_idx = 2;
          solv_sph_idx = 3;
          if((circumscribe < 0) && (surface_type == 6))
            circumscribe = 40;
          break;
        case 0:
          /* 0 - normal */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_normal);
          sphere_idx = 1;
          solv_sph_idx = 2;
          if((circumscribe < 0) && (surface_type == 6))
            circumscribe = 30;
          break;
        case -1:
          /* -1 */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_poor);
          sphere_idx = 1;
          solv_sph_idx = 2;
          if((circumscribe < 0) && (surface_type == 6))
            circumscribe = 10;
          break;
        case -2:
          /* -2 god awful */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_poor) * 1.5F;
          sphere_idx = 1;
          solv_sph_idx = 1;
          break;
        case -3:
          /* -3 miserable */
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_miserable);
          sphere_idx = 1;
          solv_sph_idx = 1;
          break;
        default:
          point_sep =
            SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                         cSetting_surface_miserable) * 1.18F;
          sphere_idx = 0;
          solv_sph_idx = 1;
        }
      }

      if((circumscribe < 0) || (!surface_solvent))
        circumscribe = 0;

      RepInit(G, &I->R);
      I->R.context.object = (void *) obj;
      I->R.context.state = state;
      I->R.fRender = (void (*)(struct Rep *, RenderInfo * info)) RepSurfaceRender;
      I->R.fFree = (void (*)(struct Rep *)) RepSurfaceFree;
      I->R.fRecolor = (void (*)(struct Rep *, struct CoordSet *)) RepSurfaceColor;
      I->R.fSameVis = (int (*)(struct Rep *, struct CoordSet *)) RepSurfaceSameVis;
      I->R.fSameColor = (int (*)(struct Rep *, struct CoordSet *)) RepSurfaceSameColor;
      I->R.fInvalidate = (void (*)(struct Rep *, struct CoordSet *, int)) RepSurfaceInvalidate;
      I->R.obj = (CObject *) (cs->Obj);
      I->R.cs = cs;
      I->allVisibleFlag = true;
      obj = cs->Obj;

      /* don't waist time computing a Surface unless we need it!! */
      {
        register int *idx_to_atm = cs->IdxToAtm;
        register AtomInfoType *obj_AtomInfo = obj->AtomInfo;
        register int a, cs_NIndex = cs->NIndex;
        for(a = 0; a < cs_NIndex; a++) {
          register AtomInfoType *ai1 = obj_AtomInfo + *(idx_to_atm++);
          if(ai1->visRep[cRepSurface] &&
             (inclH || (!ai1->hydrogen)) &&
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
          AtomInfoType *i_ai, *obj_atom_info = obj->AtomInfo;
          int *idx_to_atm = cs->IdxToAtm;
          int n_index = cs->NIndex;
          SurfaceJobAtomInfo *i_atom_info = atom_info;
          int i;
          for(i = 0; i < n_index; i++) {
            i_ai = obj_atom_info + idx_to_atm[i];
            /* just surfacing flags */
            i_atom_info->flags = i_ai->flags & (cAtomFlag_ignore | cAtomFlag_exfoliate);
            i_atom_info->vdw = i_ai->vdw;
            i_atom_info++;
          }
        }

        OrthoBusyFast(G, 0, 1);

	if (ok){
	  n_present = cs->NIndex;
	  
	  carve_selection =
	    SettingGet_s(G, cs->Setting, obj->Obj.Setting,
			 cSetting_surface_carve_selection);
	  carve_cutoff =
	    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_surface_carve_cutoff);
	  if((!carve_selection) || (!carve_selection[0]))
	    carve_cutoff = 0.0F;
	  if(carve_cutoff > 0.0F) {
	    carve_state =
	      SettingGet_i(G, cs->Setting, obj->Obj.Setting,
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
            register int *ap = present_vla;
            register int *idx_to_atm = cs->IdxToAtm;
            register AtomInfoType *obj_AtomInfo = obj->AtomInfo;
            register int a, cs_NIndex = cs->NIndex;
            for(a = 0; a < cs_NIndex; a++) {
              register AtomInfoType *ai1 = obj_AtomInfo + *(idx_to_atm++);
              if(ai1->visRep[cRepSurface] &&
                 (inclH || (!ai1->hydrogen)) &&
                 ((!cullByFlag) || (!(ai1->flags &
                                      (cAtomFlag_ignore | cAtomFlag_exfoliate)))))
                *ap = 2;
              else
                *ap = 0;
              ap++;
            }
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
            float probe_radiusX2 = probe_radius * 2;
            int a;
            for(a = 0; ok && a < cs->NIndex; a++){
              if(!present_vla[a]) {
                AtomInfoType *ai1 = obj->AtomInfo + cs->IdxToAtm[a];
                if((inclH || (!ai1->hydrogen)) &&
                   ((!cullByFlag) || 
                    !(ai1->flags & cAtomFlag_ignore))) {
                  float *v0 = cs->Coord + 3 * a;
                  int i = *(MapLocusEStart(map, v0));
                  if(optimize) {
                    if(i && map->EList) {
                      int j = map->EList[i++];
                      while(j >= 0) {
                        if(present_vla[j] > 1) {
                          AtomInfoType *ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                          if(within3f
                             (cs->Coord + 3 * j, v0,
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
          }

          if(ok && carve_flag && (!optimize)) {
            /* and optimize for carved region */
            int a;
            for(a = 0; ok && a < cs->NIndex; a++) {
              int include_flag = false;
              if(carve_map) {
                float *v0 = cs->Coord + 3 * a;
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
          if(ok && surf_job) {

            surf_job->maxVdw = I->max_vdw;
            surf_job->allVisibleFlag = I->allVisibleFlag;

            surf_job->atomInfo = atom_info;
            atom_info = NULL;

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
                float *v_src = cs->Coord;
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
	      
	      surf_job->trimCutoff = trim_cutoff;
	      surf_job->trimFactor = trim_factor;
	      
	      surf_job->cavityMode = cavity_mode;
	      surf_job->cavityRadius = cavity_radius;
	      surf_job->cavityCutoff = cavity_cutoff;
	      
	      if(carve_vla)
		surf_job->carveVla = VLACopy(carve_vla, float);
	      surf_job->carveCutoff = carve_cutoff;
	      
	      surf_job->surfaceMode = surface_mode;
	      surf_job->surfaceSolvent = surface_solvent;
	      
	      surf_job->cavityCull = SettingGet_i(G, cs->Setting,
						  obj->Obj.Setting, cSetting_cavity_cull);
	    }
	  }

          ok &= !G->Interrupt;

          if(ok) {
            int found = false;
#ifndef _PYMOL_NOPY
            PyObject *entry = NULL;
            PyObject *output = NULL;
            PyObject *input = NULL;
            int cache_mode =
              SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cache_mode);

            if(cache_mode > 0) {
              int blocked = PAutoBlock(G);
              input = SurfaceJobInputAsTuple(G, surf_job);

              if(PCacheGet(G, &output, &entry, input) == OV_STATUS_YES) {
                if(OV_OK(SurfaceJobResultFromTuple(G, surf_job, output))) {
                  found = true;
                  PXDecRef(input);
                  input = NULL;
                  PXDecRef(entry);
                  entry = NULL;
                }
                PXDecRef(output);
                output = NULL;
              }
              if(PyErr_Occurred())
                PyErr_Print();
              PAutoUnblock(G, blocked);
            }
#endif
            if(ok && !found) {

              ok &= SurfaceJobRun(G, surf_job);

#ifndef _PYMOL_NOPY
              if(cache_mode > 1) {
                int blocked = PAutoBlock(G);
                output = SurfaceJobResultAsTuple(G, surf_job);
                PCacheSet(G, entry, output);
                PXDecRef(entry);
                entry = NULL;
                PXDecRef(output);
                output = NULL;
                PXDecRef(input);
                input = NULL;
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
          RepSurfaceColor(I, cs);
      }
      if(carve_map)
        MapFree(carve_map);
      VLAFreeP(carve_vla);
      VLAFreeP(present_vla);
      if(ok && I->debug)
        ok &= CGOStop(I->debug);
      OrthoBusyFast(G, 4, 4);
    }
  }
  if(!ok) {
    RepSurfaceFree(I);
    I = NULL;
  }
  return ((void *) (struct Rep *) I);
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
  int b;
  float vdw;
  float probe_radius_plus;
  int maxDot = 0;
  int stopDot;
  int n_coord = VLAGetSize(atom_info);
  Vector3f *sp_dot = sp->dot;
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

  probe_radius_plus = probe_radius * 1.5F;

  I->nDot = 0;
  if (ok) {
    int dotCnt = 0;
    MapType *map = MapNewFlagged(G, max_vdw + probe_radius, coord, n_coord, NULL, present);
    CHECKOK(ok, map);
    ok &= !G->Interrupt;
    if(map && ok) {
      float *v = I->dot;
      float *n = I->dotNormal;
      int *dc = I->dotCode;
      int maxCnt = 0;

      ok &= MapSetupExpress(map);
      if (ok) {
        int a;
        int skip_flag;

        SurfaceJobAtomInfo *a_atom_info = atom_info;
        for(a = 0; ok && a < n_coord; a++) {
          OrthoBusyFast(G, a, n_coord * 5);
          if((!present) || (present[a])) {
            register int i;
            float *v0 = coord + 3 * a;
            vdw = a_atom_info->vdw + probe_radius;

            skip_flag = false;

            i = *(MapLocusEStart(map, v0));
            if(i && map->EList) {
              register int j = map->EList[i++];
              while(ok && j >= 0) {
                SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                if(j > a)       /* only check if this is atom trails */
                  if((!present) || present[j]) {
                    if((j_atom_info->vdw == a_atom_info->vdw)) {        /* handle singularities */
                      float *v1 = coord + 3 * j;
                      if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                        skip_flag = true;
                    }
                  }
                j = map->EList[i++];
		ok &= !G->Interrupt;
              }
            }
            if(ok && !skip_flag) {
              for(b = 0; ok && b < sp->nDot; b++) {
                float *sp_dot_b = (float*)(sp_dot + b);
                register int i;
                int flag = true;
                v[0] = v0[0] + vdw * (n[0] = sp_dot_b[0]);
                v[1] = v0[1] + vdw * (n[1] = sp_dot_b[1]);
                v[2] = v0[2] + vdw * (n[2] = sp_dot_b[2]);
                i = *(MapLocusEStart(map, v));
                if(i) {
                  register int j = map->EList[i++];
                  while(ok && j >= 0) {
                    SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                    if((!present) || present[j]) {
                      if(j != a) {
                        skip_flag = false;
                        if(j_atom_info->vdw == a_atom_info->vdw) {      /* handle singularities */
                          float *v1 = coord + 3 * j;
                          if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                            skip_flag = true;
                        }
                        if(!skip_flag)
                          if(within3f(coord + 3 * j, v, j_atom_info->vdw + probe_radius)) {
                            flag = false;
                            break;
                          }
                      }
                    }
                    j = map->EList[i++];
		    ok &= !G->Interrupt;
                  }
                }
                if(ok && flag && (dotCnt < stopDot)) {
                  dotCnt++;
                  v += 3;
                  n += 3;
                  dc++;
                  I->nDot++;
                }
              }
            }
            if(dotCnt > maxCnt) {
              maxCnt = dotCnt;
              maxDot = I->nDot - 1;
            }
          }
          a_atom_info++;
        }
      }

      /* for each pair of proximal atoms, circumscribe a circle for their intersection */

      /*    CGOReset(G->DebugCGO); */

      if (ok) {
        MapType *map2 = NULL;
        if(circumscribe && (!surface_solvent)){
          map2 = MapNewFlagged(G, 2 * (max_vdw + probe_radius), coord, n_coord, NULL, present);
	  CHECKOK(ok, map2);
	}
	ok &= !G->Interrupt;
        if(ok && map2) {
          /*        CGOBegin(G->DebugCGO,GL_LINES); */
          int a;
          int skip_flag;

          SurfaceJobAtomInfo *a_atom_info = atom_info;
          ok &= MapSetupExpress(map2);
          for(a = 0; ok && a < n_coord; a++) {
            if((!present) || present[a]) {
              register int i;
              float vdw2;
              float *v0 = coord + 3 * a;
              vdw = a_atom_info->vdw + probe_radius;
              vdw2 = vdw * vdw;

              skip_flag = false;

              i = *(MapLocusEStart(map2, v0));
              if(i) {
                register int j = map2->EList[i++];
                while(ok && j >= 0) {
                  SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                  if(j > a)     /* only check if this is atom trails */
                    if((!present) || present[j]) {
                      if((j_atom_info->vdw == a_atom_info->vdw)) {      /* handle singularities */
                        float *v2 = coord + 3 * j;
                        if((v0[0] == v2[0]) && (v0[1] == v2[1]) && (v0[2] == v2[2]))
                          skip_flag = true;
                      }
                    }
                  j = map2->EList[i++];
		  ok &= !G->Interrupt;
                }
              }

              if(ok && !skip_flag) {
                register int ii = *(MapLocusEStart(map2, v0));
                if(ii) {
                  register int jj = map2->EList[ii++];
                  while(jj >= 0) {
                    SurfaceJobAtomInfo *jj_atom_info = atom_info + jj;
                    float dist;
                    if(jj > a)  /* only check if this is atom trails */
                      if((!present) || present[jj]) {
                        float vdw3 = jj_atom_info->vdw + probe_radius;

                        float *v2 = coord + 3 * jj;
                        dist = (float) diff3f(v0, v2);
                        if((dist > R_SMALL4) && (dist < (vdw + vdw3))) {
                          float vz[3], vx[3], vy[3], vp[3];
                          float tri_a = vdw, tri_b = vdw3, tri_c = dist;
                          float tri_s = (tri_a + tri_b + tri_c) * 0.5F;
                          float area = (float) sqrt1f(tri_s * (tri_s - tri_a) *
                                                      (tri_s - tri_b) * (tri_s - tri_c));
                          float radius = (2 * area) / dist;
                          float adj = (float) sqrt1f(vdw2 - radius * radius);

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
                            int flag = true;
                            v[0] = vp[0] + vx[0] * xcosr + vy[0] * ysinr;
                            v[1] = vp[1] + vx[1] * xcosr + vy[1] * ysinr;
                            v[2] = vp[2] + vx[2] * xcosr + vy[2] * ysinr;

                            i = *(MapLocusEStart(map, v));
                            if(i && map->EList) {
                              register int j = map->EList[i++];
                              while(ok && j >= 0) {
                                SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                                if((!present) || present[j])
                                  if((j != a) && (j != jj)) {
                                    skip_flag = false;
                                    if(a_atom_info->vdw == j_atom_info->vdw) {  /* handle singularities */
                                      float *v1 = coord + 3 * j;
                                      if((v0[0] == v1[0]) &&
                                         (v0[1] == v1[1]) && (v0[2] == v1[2]))
                                        skip_flag = true;
                                    }
                                    if(jj_atom_info->vdw == j_atom_info->vdw) { /* handle singularities */
                                      float *v1 = coord + 3 * j;
                                      if((v2[0] == v1[0]) &&
                                         (v2[1] == v1[1]) && (v2[2] == v1[2]))
                                        skip_flag = true;
                                    }
                                    if(!skip_flag)
                                      if(within3f
                                         (coord + 3 * j, v,
                                          j_atom_info->vdw + probe_radius)) {
                                        flag = false;
                                        break;
                                      }
                                  }
                                j = map->EList[i++];
				ok &= !G->Interrupt;
                              }
                            }
                            if(ok && flag && (dotCnt < stopDot)) {
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

                              *dc = 1;  /* mark as exempt */

                              dotCnt++;
                              v += 3;
                              n += 3;
                              dc++;
                              I->nDot++;
                            }
                          }
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
    /*    CGOEnd(G->DebugCGO); */
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
        float *v = cavityDot;
        ok &= MapSetupExpress(map);
        if (ok) {
          int a;
          int skip_flag;
          SurfaceJobAtomInfo *a_atom_info = atom_info;
          for(a = 0; a < n_coord; a++) {
            if((!present) || (present[a])) {
              register int i;
              float *v0 = coord + 3 * a;
              vdw = a_atom_info->vdw + cavity_radius;
              skip_flag = false;
              i = *(MapLocusEStart(map, v0));
              if(i && map->EList) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                  if(j > a)       /* only check if this is atom trails */
                    if((!present) || present[j]) {
                      if((j_atom_info->vdw == a_atom_info->vdw)) {        /* handle singularities */
                        float *v1 = coord + 3 * j;
                        if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                          skip_flag = true;
                      }
                    }
                  j = map->EList[i++];
                }
              }
              if(!skip_flag) {
                for(b = 0; b < sp->nDot; b++) {
                  float *sp_dot_b = (float*)(sp_dot + b);
                  register int i;
                  int flag = true;
                  v[0] = v0[0] + vdw * (sp_dot_b[0]);
                  v[1] = v0[1] + vdw * (sp_dot_b[1]);
                  v[2] = v0[2] + vdw * (sp_dot_b[2]);
                  i = *(MapLocusEStart(map, v));
                  if(i) {
                    register int j = map->EList[i++];
                    while(j >= 0) {
                      SurfaceJobAtomInfo *j_atom_info = atom_info + j;
                      if((!present) || present[j]) {
                        if(j != a) {
                          skip_flag = false;
                          if(j_atom_info->vdw == a_atom_info->vdw) {
                            /* handle singularities */
                            float *v1 = coord + 3 * j;
                            if((v0[0] == v1[0]) && (v0[1] == v1[1]) && (v0[2] == v1[2]))
                              skip_flag = true;
                          }
                          if(!skip_flag) {
                            if(within3f(coord + 3 * j, v, j_atom_info->vdw + cavity_radius)) {
                              flag = false;
                              break;
                            }
                          }
                        }
                      }
                      j = map->EList[i++];
                    }
                  }
                  if(flag && (dotCnt < stopDot)) {
                    v += 3;
                    nCavityDot++;
                    dotCnt++;
                  }
                }
              }
            }
            a_atom_info++;
          }
        }
      }
      MapFree(map);
    }
    {
      int *dot_flag = Calloc(int, I->nDot);
      ErrChkPtr(G, dot_flag);
      {
        MapType *map = MapNew(G, cavity_cutoff, cavityDot, nCavityDot, NULL);
        if(map) {
          MapSetupExpress(map);
          {
            int *p = dot_flag;
            float *v = I->dot;
            int a;
            for(a = 0; a < I->nDot; a++) {
              register int i = *(MapLocusEStart(map, v));
              if(i && map->EList) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  if(within3f(cavityDot + (3 * j), v, cavity_cutoff)) {
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
                break;
              }
            }
          }
        }
        MapFree(map);
      }
      
      {
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
          if(!*(p++)) {
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
      FreeP(dot_flag);
    }
    VLAFreeP(cavityDot);
  }
  if(ok && (cavity_mode != 1) && (cavity_cull > 0) && 
     (probe_radius > 0.75F) && (!surface_solvent)) {
    int *dot_flag = Calloc(int, I->nDot);
    ErrChkPtr(G, dot_flag);

    {
      MapType *map = MapNew(G, probe_radius_plus, I->dot, I->nDot, NULL);
      if(map) {
        int flag = true;
        MapSetupExpress(map);
        while(ok && flag) {
          int *p = dot_flag;
          float *v = I->dot;
          int a;
          flag = false;
          for(a = 0; ok && a < I->nDot; a++) {
            if(!dot_flag[a]) {
              register int i = *(MapLocusEStart(map, v));
              int cnt = 0;

              if(i && map->EList) {
                register int j = map->EList[i++];
                while(j >= 0) {
                  if(j != a) {
                    if(within3f(I->dot + (3 * j), v, probe_radius_plus)) {
                      if(dot_flag[j]) {
                        *p = true;
                        flag = true;
                        break;
                      }
                      cnt++;
                      if(cnt > cavity_cull) {
                        *p = true;
                        flag = true;
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
        }
      }
      MapFree(map);
    }

    if (ok) {
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
        if(*(p++)) {
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

    FreeP(dot_flag);
  }
  if(!ok) {
    SolventDotFree(I);
    I = NULL;
  }
  return I;
}
