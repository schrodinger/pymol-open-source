
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
#include"OOMac.h"
#include"RepSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"Util.h"
#include"Feedback.h"
#include "ShaderMgr.h"
#include "Scene.h"
#include"CGO.h"
#include"ObjectMolecule.h"
#include "Lex.h"

#include "ShaderText.h"

#ifdef NT
#undef NT
#endif

/* defer_builds_mode = 5 : Immediate mode for any sphere_mode

   sphere_mode :

   0) Geometry shaders (quality based on sphere_quality, default 1)
   1) rectangular points with the same size, that can be changed with sphere_point_size
   2) rectangles with constant size relative to vdw and scene scale (i.e., changes when zoomed)
      maxed by a multiple of sphere_point_max_size (set it below 1 to see it influence, 3*pixel_scale max)
   3) same as 2 but with circles
   4) Draw multiple points for each sphere to mimic specular reflection of one light on the sphere
   5) Uses the fast ARB Shader that approximates spheres as circles
   9) GLSL Shader Spheres (only when shaders are available)

 */
typedef struct RepSphere {
  Rep R;
  float *V;    /* triangle vertices (if any) */
  float *VC;   /* 8 floats per sphere: 3 color, 1 transparency, 3 coordinates, 1 radius */
  float *VN;   /* normals (if any) computed in RepSphereComputeSphereNormals() */
  SphereRec *SP;
  SphereRec *SSP;  /* sphere record used for picking (if set, otherwise 0) and direct mode */
  int *NT;
  int N, NC, NP;   /* N - number of vertices in V triangles, 
		      NC - number of spheres stored in VC, 
		      NP - number of pickable atoms */
  int cullFlag;
  int spheroidFlag;
  int *LastVisib;
  int *LastColor;
  float LastVertexScale;
  int VariableAlphaFlag;
  CGO *shaderCGO;
} RepSphere;

void RepSphereFree(RepSphere * I);
int RepSphereSameVis(RepSphere * I, CoordSet * cs);

void RepSphereFree(RepSphere * I)
{
  if (I->shaderCGO ){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }

  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->VN);
  FreeP(I->NT);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  RepPurge(&I->R);
  OOFreeP(I);
}

/* MULTI-INSTSANCE TODO:  isn't this a conflict? */
static CShaderPrg *sphereARBShaderPrg = NULL;

#ifdef _PYMOL_ARB_SHADERS
static void RepSphereRenderOneSphere_ARB(PyMOLGlobals *G, RenderInfo *info,
    float *color, float *last_radius, float *cur_radius, float *fog_info,
    float *v)
{
  static const float _00[2] = { 0.0F, 0.0F };
  static const float _01[2] = { 0.0F, 1.0F };
  static const float _11[2] = { 1.0F, 1.0F };
  static const float _10[2] = { 1.0F, 0.0F };

  float v3 = v[3];
  if((*last_radius) != ( (*cur_radius) = v3)) {
    glEnd();
    glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
	  0, 0.0F, 0.0F, v3, 0.0F);
    glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
			       0, fog_info[0], fog_info[1], 0.0F, 0.0F);
    glBegin(GL_QUADS);
    (*last_radius) = (*cur_radius);
  }
  glColor3fv(color);
  glTexCoord2fv(_00);
  glVertex3fv(v);
  glTexCoord2fv(_10);
  glVertex3fv(v);
  glTexCoord2fv(_11);
  glVertex3fv(v);
  glTexCoord2fv(_01);
  glVertex3fv(v);
}
#endif

static void RenderSpherePopulateVariables(PyMOLGlobals *G, RenderInfo *info,
    float *nv, float *fog_info, float *z_front, float *z_back)
{
  /* compute -Ze = (Wc) of fog start */
  nv[3] =
    (info->front +
     (info->back - info->front) * SettingGetGlobal_f(G, cSetting_fog_start));
  /* compute Zc of fog start using std. perspective transformation */
  nv[2] =
    (nv[3] * (info->back + info->front) -
     2 * (info->back * info->front)) / (info->back - info->front);
  /* compute Zc/Wc to get normalized depth coordinate of fog start */
  nv[0] = (nv[2] / nv[3]);
  fog_info[0] = (nv[0] * 0.5) + 0.5;
  
  fog_info[1] = 1.0F / (1.0 - fog_info[0]);     /* effective range of fog */
  
  (*z_front) = info->stereo_front;
  (*z_back) = info->back + ((info->back + info->front) * 0.25);
}

#ifdef _PYMOL_ARB_SHADERS
void RenderSphereMode_Immediate_5(PyMOLGlobals *G, RenderInfo *info, CoordSet *cs, ObjectMolecule *obj, int *repActive, float sphere_scale){
  if (!sphereARBShaderPrg){
    sphereARBShaderPrg = CShaderPrg_NewARB(G, "sphere_arb", sphere_arb_vs, sphere_arb_fs);
  }
  if(sphereARBShaderPrg){
    float fog_info[3];
    float nv[4];
    float z_front, z_back;

    RenderSpherePopulateVariables(G, info, nv, fog_info, &z_front, &z_back);
    
    CShaderPrg_Enable_SphereShaderARB(G);
    
    glNormal3fv(info->view_normal);
    glBegin(GL_QUADS);
    {
      float last_radius = -1.0F, cur_radius;
      int a;
      int nIndex = cs->NIndex;
      AtomInfoType *atomInfo = obj->AtomInfo;
      int *i2a = cs->IdxToAtm;
      float *v = cs->Coord;
      for(a = 0; a < nIndex; a++) {
	AtomInfoType *ai = atomInfo + *(i2a++);
	if(GET_BIT(ai->visRep,cRepSphere)) {
	  float vr[4];
	  copy3f(v, vr);
	  vr[3] = ai->vdw * sphere_scale;
	  (*repActive) = true;
	  RepSphereRenderOneSphere_ARB(G, info, ColorGet(G, ai->color), &last_radius, &cur_radius, fog_info, vr);
	}
	v += 3;
      }
      glEnd();
    }
    CShaderPrg_DisableARB(sphereARBShaderPrg);
  }
}
#endif

static void RenderSphereMode_Immediate_4(PyMOLGlobals *G, RenderInfo *info,
    CoordSet *cs, ObjectMolecule *obj, int *repActive, float pixel_scale)
{
#ifndef PURE_OPENGL_ES_2
        int repeat = true;
        const float _1 = 1.0F;
        const float _2 = 2.0F;
        float x_add = 0.0F, y_add = 0.0F, z_add = 0.0F;
        float z_factor = 0.0F, r_factor = 1.0F;
        float s_factor = 0.0F;
        int pass = 0;
        float max_size = SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                                               cSetting_sphere_point_max_size);
        int clamp_size_flag = (max_size >= 0.0F);
        while(repeat) {

          int a;
          int nIndex = cs->NIndex;
          AtomInfoType *atomInfo = obj->AtomInfo;
          int *i2a = cs->IdxToAtm;
          float *v = cs->Coord;
          float last_radius = -1.0F;
          float last_size = -1.0;
          float largest = 0.0F;

          float zz_factor = _1 - (float) pow(_1 - z_factor, 2);
          if(zz_factor < 0.45F)
            zz_factor = 0.45F;

          repeat = false;

          glBegin(GL_POINTS);
	  
          for(a = 0; a < nIndex; a++) {
            AtomInfoType *ai = atomInfo + *(i2a++);
            if(GET_BIT(ai->visRep,cRepSphere)) {
              float cur_radius = ai->vdw;
              (*repActive) = true;

              if(last_radius != cur_radius) {
                float clamp_radius = cur_radius;
                float size = cur_radius * pixel_scale;

                if(clamp_size_flag)
                  if(size > max_size) {
                    size = max_size;
                    clamp_radius = size / pixel_scale;
                  }
                size *= r_factor;
                if(size != last_size) {
                  glEnd();
                  if(size > largest)
                    largest = size;
                  if(size < _2) {
                    if(!pass) {
                      zz_factor = 1.0F;
                      s_factor = 0.0F;
                    }
                  }
                  if(size < _1) {
                    size = _1;
                    glDisable(GL_POINT_SMOOTH);
                    glDisable(GL_ALPHA_TEST);
                  } else {
                    glEnable(GL_POINT_SMOOTH);
                    glEnable(GL_ALPHA_TEST);
                  }
                  glPointSize(size);
                  glBegin(GL_POINTS);
                }

                x_add = z_factor * clamp_radius * info->view_normal[0];
                y_add = z_factor * clamp_radius * info->view_normal[1];
                z_add = z_factor * clamp_radius * info->view_normal[2];
                last_radius = cur_radius;
                last_size = size;
              }

              {
                float *vc = ColorGet(G, ai->color);
                float r = zz_factor * vc[0] + s_factor;
                float g = zz_factor * vc[1] + s_factor;
                float b = zz_factor * vc[2] + s_factor;

                glColor3f(r > _1 ? _1 : r, g > _1 ? _1 : g, b > _1 ? _1 : b);
                glVertex3f(v[0] + x_add, v[1] + y_add, v[2] + z_add);
              }
            }
            v += 3;
          }

          glEnd();

          if(largest > 2.0F) {
            float reduce = (largest - 2.0F) / largest;
            r_factor *= reduce;
            z_factor = (float) sqrt1f(1.0F - (r_factor * r_factor));
            s_factor = (float) pow(z_factor, 20.0F) * 0.5F;
            repeat = true;
            pass++;
          }
        }
        glDisable(GL_POINT_SMOOTH);
	printf("pass=%d\n", pass);
#endif
}

static void RenderSphereMode_Immediate_Triangles(PyMOLGlobals *G, CoordSet *cs,
    ObjectMolecule *obj, int *repActive, float sphere_scale)
{
  /* triangle-based spheres */
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_quality);
  if(ds < 0) {
    sp = NULL;
  } else {
    if(ds > 4)
      ds = 4;
    sp = G->Sphere->Sphere[ds];
  }
  {
    int a;
    int nIndex = cs->NIndex;
    AtomInfoType *atomInfo = obj->AtomInfo;
    int *i2a = cs->IdxToAtm;
    int last_color = -1;
    float *v = cs->Coord;
    int *sp_Sequence = sp->Sequence;
    int *sp_StripLen = sp->StripLen;
    int sp_NStrip = sp->NStrip;
    Vector3f *sp_dot = sp->dot;
    
    for(a = 0; a < nIndex; a++) {
      AtomInfoType *ai = atomInfo + *(i2a++);
      if(GET_BIT(ai->visRep,cRepSphere)) {
	float vdw = ai->vdw * sphere_scale;
	int c = ai->color;
	float v0 = v[0];
	float v1 = v[1];
	float v2 = v[2];
	(*repActive) = true;
	
	if(c != last_color) {
	  last_color = c;
	  glColor3fv(ColorGet(G, c));
	}
	
	{
	  int *s = sp_StripLen;
	  int *q = sp_Sequence;
	  int b;
	  for(b = 0; b < sp_NStrip; b++) {
	    int nc = *(s++);
#ifdef PURE_OPENGL_ES_2
	    /* TODO */
#else
	    glBegin(GL_TRIANGLE_STRIP);
	    for(c = 0; c < nc; c++) {
	      float *sp_dot_q = &sp_dot[*(q++)][0];
	      glNormal3fv(sp_dot_q);        /* normal */
	      glVertex3f(v0 + vdw * sp_dot_q[0],
			 v1 + vdw * sp_dot_q[1], v2 + vdw * sp_dot_q[2]);
	    }
	    glEnd();
#endif
	  }
	}
      }
      v += 3;
    }
  }
}

static void RenderSphereMode_Immediate_1_2_3(PyMOLGlobals *G, RenderInfo *info,
    CoordSet *cs, ObjectMolecule *obj, int *repActive, float pixel_scale,
    int sphere_mode)
{
  /* sphere_mode is 1, 2, or 3 */
  float max_radius = SettingGet_f(G, cs->Setting, obj->Obj.Setting,
				  cSetting_sphere_point_max_size) * 3 * pixel_scale;
  int clamp_size_flag = (max_radius >= 0.0F);
  
  int a;
  int nIndex = cs->NIndex;
  AtomInfoType *atomInfo = obj->AtomInfo;
  int *i2a = cs->IdxToAtm;
  int last_color = -1;
  float *v = cs->Coord;
  float last_radius = -1.0F;
  
  if(!info->line_lighting)
    glDisable(GL_LIGHTING);
  
  glBegin(GL_POINTS);
  for(a = 0; a < nIndex; a++) {
    AtomInfoType *ai = atomInfo + *(i2a++);
    if(GET_BIT(ai->visRep,cRepSphere)) {
      int c = ai->color;
      (*repActive) = true;
      if(c != last_color) {
	last_color = c;
	glColor3fv(ColorGet(G, c));
      }
      switch (sphere_mode) {
      case 1:
      case 6:
	glVertex3fv(v);
	break;
      case 2:
      case 3:
      case 7:
      case 8:
	{
	  float cur_radius = ai->vdw * pixel_scale;
	  if(last_radius != cur_radius) {
	    glEnd();
	    if(clamp_size_flag)
	      if(cur_radius > max_radius)
		cur_radius = max_radius;
	    glPointSize(cur_radius);
	    glBegin(GL_POINTS);
	    last_radius = cur_radius;
	  }
	  glVertex3fv(v);
	}
	break;
      }
    }
    v += 3;
  }
  glEnd();

  glEnable(GL_LIGHTING);
  
  if(sphere_mode == 3) {
    glDisable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.05F);
  } else {
    glEnable(GL_ALPHA_TEST);
  }
}

void RenderImmediate_DoPreGL(PyMOLGlobals *G, int sphere_mode, float *pixel_scale, CoordSet *cs, ObjectMolecule *obj, float sphere_scale){
  switch (sphere_mode) {
  case 2:
  case 7:
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_ALPHA_TEST);
    (*pixel_scale) *= 1.4F;
    glPointSize(1.0F);
    break;
  case 3:
  case 8:
    glEnable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.5F);
    glEnable(GL_ALPHA_TEST);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glPointSize(1.0F);
    (*pixel_scale) *= 2.0F;
    break;
  case 4:
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_ALPHA_TEST);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glPointSize(1.0F);
    (*pixel_scale) *= 2.0F;
    break;
  default:
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_ALPHA_TEST);
    glPointSize(SettingGet_f
		(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_point_size));
    break;
  }
}

void RepSphereRenderImmediate(CoordSet * cs, RenderInfo * info)
{
  PyMOLGlobals *G = cs->State.G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    int repActive = false;
    ObjectMolecule *obj = cs->Obj;
    int sphere_mode =
      SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_mode);
    float sphere_scale =
      SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_scale);

    if(sphere_mode > 0) {       /* point-based modees */
      float pixel_scale = 1.0F / info->vertex_scale;
      RenderImmediate_DoPreGL(G, sphere_mode, &pixel_scale, cs, obj, sphere_scale);
      switch (sphere_mode){
#ifdef _PYMOL_ARB_SHADERS
      case 5:
	RenderSphereMode_Immediate_5(G, info, cs, obj, &repActive, sphere_scale);
	break;
#endif
      case 4:
	RenderSphereMode_Immediate_4(G, info, cs, obj, &repActive, pixel_scale);
	break;
      default:
	RenderSphereMode_Immediate_1_2_3(G, info, cs, obj, &repActive, pixel_scale, sphere_mode);
      }
    } else {
      RenderSphereMode_Immediate_Triangles(G, cs, obj, &repActive, sphere_scale);
    }

    if(!repActive)              /* didn't draw a single sphere, so we can skip this representation next time around */
      cs->Active[cRepSphere] = false;
  }
}
#ifdef PURE_OPENGL_ES_2
void RepSphereRenderPointsES(int nvertsarg, float *varg, float *vnarg);
void RepSphereRenderMode5PointsES(int nvertsarg, float *varg, float zz_factor, float s_factor, float *dim_add, float _1);
void RepSphereRenderPointsDefaultES(RepSphere * I, Picking **pick, int nvertsarg, int iarg, Pickable *parg, float *varg);
void RepSphereRenderPointsES(int nvertsarg, float *varg, float *vnarg){
    /* TODO */
}
void RepSphereRenderMode5PointsES(int nvertsarg, float *varg, float zz_factor, float s_factor, float *dim_add, float _1){
}
void RepSphereRenderPointsDefaultES(RepSphere * I, Picking **pick, int nvertsarg, int iarg, Pickable *parg, float *varg){
}
#endif


static int RenderSphereMode_Direct(PyMOLGlobals *G, RepSphere *I,
    RenderInfo * info, int carg, float **vptr, float alpha,
    SphereRec *sphereRecPtr)
{
  short use_shader, generate_shader_cgo = 0;
  float *v = *vptr;
  int c = carg;
  int ok = true;
  use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) & 
               SettingGetGlobal_b(G, cSetting_use_shaders);
  if (I->shaderCGO && !use_shader){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }

  if (use_shader){
    if (!I->shaderCGO){
      I->shaderCGO = CGONew(G);
      CHECKOK(ok, I->shaderCGO);
      if (ok)
	I->shaderCGO->use_shader = true;
      generate_shader_cgo = 1;
    } else {
      I->shaderCGO->enable_shaders = true;
      CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
      return true;
    }
  }

  if (generate_shader_cgo){
    if(ok && sphereRecPtr) {
      float last_vdw = -1.0F;
      int variable_alpha = I->VariableAlphaFlag;
      SphereRec *sp = sphereRecPtr;
      while(ok && c--) {
	Vector3f *sp_dot = sp->dot;
	int b, *q, *s;
	if(variable_alpha){
	  ok &= CGOAlpha(I->shaderCGO, v[3]);
	} else {
	  ok &= CGOAlpha(I->shaderCGO, alpha);
	}
	if (ok)
	  ok &= CGOColorv(I->shaderCGO, v);
	(*vptr)+=4; v = *vptr;

	if (ok){
	  float vdw = v[3];
	  last_vdw = vdw;
	  q = sp->Sequence;
	  s = sp->StripLen;
	  for(b = 0; ok && b < sp->NStrip; b++) {
	    int d;
	    ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
	    for(d = 0; ok && d < (*s); d++) {
	      float *norm = sp_dot[*(q++)];
	      ok &= CGONormalv(I->shaderCGO, norm);
	      if (ok)
		ok &= CGOVertex(I->shaderCGO, v[0] + vdw * norm[0], v[1] + vdw * norm[1], v[2] + vdw * norm[2]);
	    }
	    if (ok)
	      ok &= CGOEnd(I->shaderCGO);
	    s++;
	  }
	}
	(*vptr)+=4; v = *vptr;
      }
      if (ok)
	ok &= CGOStop(I->shaderCGO);
    }
    if (ok) {
      CGO *convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0), *convertcgo2;
      CHECKOK(ok, convertcgo);
      if (ok){
	CGOFree(I->shaderCGO);
	I->shaderCGO = convertcgo;
	convertcgo2 = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0);
	CHECKOK(ok, convertcgo2);
	if (ok){
	  CGOFree(I->shaderCGO);
	  I->shaderCGO = convertcgo2;
	}
      }
    }
    if (ok){
      I->shaderCGO->enable_shaders = true;
      CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
    }
  } else {
    if(sphereRecPtr) {
      float last_vdw = -1.0F;
      int dlist = 0;
      int variable_alpha = I->VariableAlphaFlag;
      SphereRec *sp = sphereRecPtr;
      while(c--) {
	Vector3f *sp_dot = sp->dot;
	int b, *q, *s;
	if(variable_alpha)
	  glColor4f(v[0], v[1], v[2], v[3]);
	else
	  glColor4f(v[0], v[1], v[2], alpha);
	(*vptr)+=4; v = *vptr;
	{
	  float vdw = v[3];
	  glTranslatef(v[0], v[1], v[2]);
	  if((vdw != last_vdw) || (!dlist)) {
	    last_vdw = vdw;
	    q = sp->Sequence;
	    s = sp->StripLen;

#ifndef PURE_OPENGL_ES_2
	    for(b = 0; b < sp->NStrip; b++) {
	      int d;
	      glBegin(GL_TRIANGLE_STRIP);
	      for(d = 0; d < (*s); d++) {
		float *norm = sp_dot[*(q++)];
		glNormal3fv(norm);
		glVertex3f(vdw * norm[0], vdw * norm[1], vdw * norm[2]);
	      }
	      glEnd();
	      s++;
	    }
#endif
	  }
	  glTranslatef(-v[0], -v[1], -v[2]);
	}
	(*vptr)+=4; v = *vptr;
      }
    }
  }
  
  if (!ok){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
    I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
    I->R.cs->Active[cRepSphere] = false;
  }
  return ok;
}
   //, float radius, int carg, float pixel_scale, int clamp_size_flag, float max_size){

static void RenderSphereMode_Sprites(PyMOLGlobals *G, RepSphere *I,
    RenderInfo *info, int sphere_mode, int carg, float **vptr, float **vnptr)
{
   int c = carg;
   float pixel_scale = 1.0F / info->vertex_scale;
   float last_radius = -1.0F, cur_radius;
   float size;
   float *v = *vptr, *vn = *vnptr;
   float max_size =
     SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting,
		  cSetting_sphere_point_max_size);
   int clamp_size_flag = (max_size >= 0.0F);

  if((sphere_mode == 3) || (sphere_mode == 8)) {
    pixel_scale *= 2.0F;
    glEnable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.5F);
    glEnable(GL_ALPHA_TEST);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glPointSize(1.0F);
  } else {
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_ALPHA_TEST);
    pixel_scale *= 1.4F;
  }
  if((sphere_mode == 7) || (sphere_mode == 8))
    glEnable(GL_LIGHTING);

  glBegin(GL_POINTS);
  while(c--) {
    if(last_radius != (cur_radius = v[7])) {
      size = cur_radius * pixel_scale;
      glEnd();
      if(clamp_size_flag)
	if(size > max_size)
	  size = max_size;
      glPointSize(size);
      glBegin(GL_POINTS);
      last_radius = cur_radius;
    }
    glColor3fv(v);
    (*vptr)+=4; v = *vptr;
    if(vn) {
      glNormal3fv(vn);
      (*vnptr)+=3; vn = *vnptr;
    }
    glVertex3fv(v);
    (*vptr)+=4; v = *vptr;
  }
  glEnd();

  if(sphere_mode == 3) {
    glDisable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.05F);
  } else {
    glEnable(GL_ALPHA_TEST);
  }
}
 
static void RenderSphereMode_Points(PyMOLGlobals *G, RepSphere *I, RenderInfo *info, int carg){
#ifndef PURE_OPENGL_ES_2
  float _1 = 1.0F;
  float _2 = 2.0F;
  float pixel_scale = 1.0F / info->vertex_scale;
  int repeat = true;
  float x_add = 0.0F, y_add = 0.0F, z_add = 0.0F;
  float z_factor = 0.0F, r_factor = 1.0F;
  float largest;
  float r, g, b;
  float s_factor = 0.0F;
  float zz_factor;
  float clamp_radius;
  float last_size;
  int pass = 0;
  float *v, last_radius, cur_radius, size;
  int c;
  float max_size =
    SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting,
		 cSetting_sphere_point_max_size);
  int clamp_size_flag = (max_size >= 0.0F);

  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_ALPHA_TEST);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glPointSize(1.0F);
  
  pixel_scale *= 2.0F;
  while(repeat) {
    v = I->VC;
    c = I->NC;
    largest = 0.0F;
    zz_factor = _1 - (float) pow(_1 - z_factor, 2);
    if(zz_factor < 0.45F)
      zz_factor = 0.45F;
    
    last_radius = -1.0F;
    last_size = -1.0F;
    repeat = false;

  glBegin(GL_POINTS);
  while(c--) {
    if(last_radius != (cur_radius = v[7])) {
      size = cur_radius * pixel_scale;
      clamp_radius = cur_radius;
      if(clamp_size_flag)
	if(size > max_size) {
	  size = max_size;
	  clamp_radius = size / pixel_scale;
	}
      size *= r_factor;
      if(size != last_size) {
	glEnd();
	if(size > largest)
	  largest = size;
	if(size < _2) {
	  if(!pass) {
	    zz_factor = 1.0F;
	    s_factor = 0.0F;
	  }
	}
	if(size < _1) {
	  size = _1;
	  glDisable(GL_POINT_SMOOTH);
	  glDisable(GL_ALPHA_TEST);
	} else {
	  glEnable(GL_POINT_SMOOTH);
	  glEnable(GL_ALPHA_TEST);
	}
	glPointSize(size);
	glBegin(GL_POINTS);
      }
      x_add = z_factor * clamp_radius * info->view_normal[0];
      y_add = z_factor * clamp_radius * info->view_normal[1];
      z_add = z_factor * clamp_radius * info->view_normal[2];
      last_radius = cur_radius;
      last_size = size;
    }
    r = zz_factor * v[0] + s_factor;
    g = zz_factor * v[1] + s_factor;
    b = zz_factor * v[2] + s_factor;
    
    glColor3f(r > _1 ? _1 : r, g > _1 ? _1 : g, b > _1 ? _1 : b);
    
    v += 4;
    glVertex3f(v[0] + x_add, v[1] + y_add, v[2] + z_add);
    v += 4;
  }
  glEnd();

  if(largest > 2.0F) {
    float reduce = (largest - 2.0F) / largest;
    r_factor *= reduce;
    z_factor = (float) sqrt1f(1.0F - (r_factor * r_factor));
    s_factor = (float) pow(z_factor, 20.0F) * 0.5F;
    repeat = true;
    pass++;
  }
}
 glDisable(GL_POINT_SMOOTH);
#endif
}

static void RenderSphereMode_9(PyMOLGlobals *G, RepSphere *I, RenderInfo *info, float **vptr, int carg){
  int c = carg;
  short use_shader;
  float *v = *vptr;
  use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) & 
               SettingGetGlobal_b(G, cSetting_use_shaders);

  if (I->shaderCGO && !use_shader){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  if (use_shader){
    if (!I->shaderCGO){
      I->shaderCGO = CGONew(G);
      I->shaderCGO->use_shader = true;

      // generating shader
#ifndef PURE_OPENGL_ES_2
      CGOEnable(I->shaderCGO, GL_LIGHTING);
#endif
      while (c--) {
	CGOAlpha(I->shaderCGO, v[3]);
	CGOColorv(I->shaderCGO, v);
	CGOSphere(I->shaderCGO, v+4, v[7]);
	(*vptr)+=8; v = *vptr;
      }
      CGOStop(I->shaderCGO);
      {
	CGO *convertcgo = NULL;
	convertcgo = CGOOptimizeSpheresToVBONonIndexed(I->shaderCGO, 0, true);
	if (convertcgo){
	  CGOFree(I->shaderCGO);    
	  I->shaderCGO = convertcgo;
	}
      }
    }
    if (I->shaderCGO){
      I->shaderCGO->enable_shaders = true;
      CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
    }
  }
}

#ifdef _PYMOL_ARB_SHADERS
void RenderSphereMode_ARB(PyMOLGlobals *G, RenderInfo *info, float **vptr, int carg){
  int c = carg;
  float fog_info[3];
  float nv[4];
  float z_front, z_back;
  float *v = *vptr;
  float last_radius, cur_radius;
  
  RenderSpherePopulateVariables(G, info, nv, fog_info, &z_front, &z_back);
  
  if(Feedback(G, FB_OpenGL, FB_Debugging))
    PyMOLCheckOpenGLErr("before shader");
  CShaderPrg_Enable_SphereShaderARB(G);
  {
    glNormal3fv(info->view_normal);
    (*vptr)+=4; v = *vptr;
    last_radius = -1.f;
    
    glBegin(GL_QUADS);
    while(c--) {
      RepSphereRenderOneSphere_ARB(G, info, v - 4, &last_radius, &cur_radius, fog_info, v);
      (*vptr)+=8; v = *vptr;
    }
    glEnd();
    CShaderPrg_DisableARB(sphereARBShaderPrg);
    if(Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("after shader");
  }
}
#endif

/* simple, default point width points -- modes 1 or 6 */
static void RenderSphereMode_1_or_6(PyMOLGlobals *G, RepSphere *I,
    RenderInfo *info, float **vptr, float **vnptr, int carg, float alpha)
{
#ifndef PURE_OPENGL_ES_2
  int c = carg;
  float *v = *vptr, *vn = *vnptr;

  glPointSize(SettingGet_f
	      (G, I->R.cs->Setting, I->R.obj->Setting,
	       cSetting_sphere_point_size));
  glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_ALPHA_TEST);
  
  glBegin(GL_POINTS);
  if(alpha == 1.0) {
    if(vn) {
      glEnd();
      glEnable(GL_LIGHTING);
      glBegin(GL_POINTS);
      while(c--) {
	glColor3fv(v);
	(*vptr)+=4; v = *vptr;
	glNormal3fv(vn);
	(*vnptr)+=3; vn = *vnptr;
	glVertex3fv(v);
	(*vptr)+=4; v = *vptr;
      }
    } else {
      while(c--) {
	glColor3fv(v);
	(*vptr)+=4; v = *vptr;
	glVertex3fv(v);
	(*vptr)+=4; v = *vptr;
      }
    }
  } else {
    if(vn) {
      glEnd();
      glEnable(GL_LIGHTING);
      glBegin(GL_POINTS);
      while(c--) {
	glColor4f(v[0], v[1], v[2], alpha);
	(*vptr)+=4; v = *vptr;
	glNormal3fv(vn);
	(*vnptr)+=3; vn = *vnptr;
	glVertex3fv(v);
	(*vptr)+=4; v = *vptr;
      }
    } else {
      while(c--) {
	glColor4f(v[0], v[1], v[2], alpha);
	(*vptr)+=4; v = *vptr;
	glVertex3fv(v);
	(*vptr)+=4; v = *vptr;
      }
    }
  }
  glEnd();

  glEnable(GL_ALPHA_TEST);
#endif
}

static void RepSpheresPrepPickingIfNoSphereGeometry(RepSphere * I,
    int sphere_mode, float *pixel_scale)
{
  PyMOLGlobals *G = I->R.G;
  switch (sphere_mode) {
  case 2:
  case 7:
    (*pixel_scale) *= 1.4F;
    glPointSize(1.0F);
    break;
  case 3:
  case 8:
    (*pixel_scale) *= 2.0F;
    glPointSize(1.0F);
    break;
  default:
    {
      float ptsize = SettingGet_f
	(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_sphere_point_size);
      glPointSize(ptsize);
    }
    break;
  }
}

static void RepSpheresSetColorForPicking(RepSphere * I, Picking **pick, int *i, int *j, Pickable **p){
  (*i)++;
  if(!(*pick)[0].src.bond) {
    /* pass 1 - low order bits *            */
    glColor3ub((uchar) ((*i & 0xF) << 4), (uchar) ((*i & 0xF0) | 0x8),
	       (uchar) ((*i & 0xF00) >> 4));
    VLACheck((*pick), Picking, *i);
    (*p)++;
    (*pick)[*i].src = **p;    /* copy object and atom info */
    (*pick)[*i].context = I->R.context;
  } else {
    /* pass 2 - high order bits */
    (*j) = *i >> 12;
    glColor3ub((uchar) ((*j & 0xF) << 4), (uchar) ((*j & 0xF0) | 0x8),
	       (uchar) ((*j & 0xF00) >> 4));
  }
}

static void RepSpheresRenderSphereGeometryForPicking(SphereRec *sp, float *v0, float vdw){
  int *s, *q, b, cc;
  q = sp->Sequence;
  s = sp->StripLen;
  for(b = 0; b < sp->NStrip; b++) {
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else

    glBegin(GL_TRIANGLE_STRIP);
    for(cc = 0; cc < (*s); cc++) {
      glNormal3f(sp->dot[*q][0], sp->dot[*q][1], sp->dot[*q][2]);
      glVertex3f(v0[0] + vdw * sp->dot[*q][0],
		 v0[1] + vdw * sp->dot[*q][1],
		 v0[2] + vdw * sp->dot[*q][2]);
      q++;
    }
    glEnd();
#endif
    s++;
  }
}

static void RepSpheresRenderSphereRecAtVertex(SphereRec *sp, float *v0, float vdw){
#ifndef PURE_OPENGL_ES_2
  Vector3f *sp_dot = sp->dot;
  int b, *q, *s;
  {
    glTranslatef(v0[0], v0[1], v0[2]);
    q = sp->Sequence;
    s = sp->StripLen;
    for(b = 0; b < sp->NStrip; b++) {
      int d;
      glBegin(GL_TRIANGLE_STRIP);
      for(d = 0; d < (*s); d++) {
	float *norm = sp_dot[*(q++)];
	glNormal3fv(norm);
	glVertex3f(vdw * norm[0], vdw * norm[1], vdw * norm[2]);
      }
      glEnd();
      s++;
    }
    glTranslatef(-v0[0], -v0[1], -v0[2]);
  }
#endif
}

static void RepSpheresRenderPointForPicking(RepSphere * I, float vdw, float *v,
    int sphere_mode, float *last_radius, float *cur_radius, float pixel_scale,
    int clamp_size_flag, float max_size, short *hasBegun)
{
#ifndef PURE_OPENGL_ES_2
  float size;
  float *v0 = &v[4];
  SphereRec *sp = I->R.G->Sphere->Sphere[0];
  switch (sphere_mode) {
  case -1:
  case 0:          /* memory-efficient sphere rendering */
    if(I->SSP) {
      sp = I->SSP;
    }
    RepSpheresRenderSphereRecAtVertex(sp, v0, vdw);
    break;
  case 2:
  case 3:
  case 4:
  case 5:               
  case 7:
  case 8:
    {
      {
	(*cur_radius) = v[7];
	size = (*cur_radius) * pixel_scale;
	if (*hasBegun){
	  glEnd();
	  (*hasBegun) = 0;
	}
	if(clamp_size_flag)
	  if(size > max_size)
	    size = max_size;
	glPointSize(size);
	glBegin(GL_POINTS);
	(*hasBegun) = 1;
	*last_radius = *cur_radius;
      }
    }
    glVertex3fv(v0);
    break;
  default:         /* simple, default point width points */
    glVertex3fv(v0);
    break;
  }
#endif
}

static void RepSpheresRenderEndOfPicking(int sphere_mode){
  switch (sphere_mode) {
  case -1:
  case 0:
    break;
  case 3:
  case 4:
  case 8:
    glDisable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.05F);
    break;
  default:
    glEnable(GL_ALPHA_TEST);
    break;
  }
}

static int RepSphereRenderRay(RepSphere * I, RenderInfo * info, float alpha)
{
  CRay *ray = info->ray;
  SphereRec *sp = I->SP;
  int ok = true;
  int c = I->N;
  int cc = 0;
  float *v = I->V, *vc;
  int a;
  ray->transparentf(1.0 - alpha);
  if(I->spheroidFlag) {
    if(sp) {
      while(c--) {
	vc = v;
	v += 3;
	for(a = 0; ok && a < sp->NStrip; a++) {
	  cc = sp->StripLen[a];
	  while(ok && (cc--) > 2) {
	    ok &= ray->triangle3fv(v + 3, v + 9, v + 15, v, v + 6, v + 12, vc, vc, vc);
	    v += 6;
            }
	  v += 12;
	}
      }
    }
  } else {
    int variable_alpha = I->VariableAlphaFlag;
    v = I->VC;
    c = I->NC;
    while(ok && c--) {
      if(variable_alpha) {
	ray->transparentf(1.0F - v[3]);
      }
      ray->color3fv(v);
      v += 4;
      ok &= ray->sphere3fv(v, *(v + 3));
      v += 4;
    }
  }
  ray->transparentf(0.0);
  return ok;
}

static void RepSphereRenderPick(RepSphere * I, RenderInfo * info, float alpha, int sphere_mode)
{
  PyMOLGlobals *G = I->R.G;
  SphereRec *sp = NULL;
  Picking **pick = info->pick;
  int c = I->N;
  float *v = I->V;
  int a;
  int cc = 0;
  int trans_pick_mode = SettingGet_i(G, I->R.cs->Setting,
				     I->R.obj->Setting,
				     cSetting_transparency_picking_mode);
  
  switch (sphere_mode){
  case 0: case 4: case 5: case 9:
    // for picking these modes, use geometry
    sp = I->R.G->Sphere->Sphere[0];
  }
  SceneSetupGLPicking(G);
  
  if(I->R.P && ((trans_pick_mode == 1) || ((trans_pick_mode == 2) && (alpha > 0.9F)))) {
    int i, j;
    Pickable *p;
    i = (*pick)->src.index;
    
    p = I->R.P;
    
    if(I->spheroidFlag && sp) {
      while(c--) {
	int skip = (p[1].index < 0);
	if(!skip) {
	  RepSpheresSetColorForPicking(I, pick, &i, &j, &p);
	} else {
	  p++;
	}
	
	v += 4;
	for(a = 0; a < sp->NStrip; a++) {
	  cc = sp->StripLen[a];
	  if(!skip) {
	    glBegin(GL_TRIANGLE_STRIP);
	    while((cc--) > 0) {
	      glNormal3fv(v);
	      glVertex3fv(v + 3);
	      v += 6;
	    }
	    glEnd();
	  } else {
	    while((cc--) > 0) {
	      v += 6;
	    }
	  }
	}
      }
    } else {
      float last_radius = -1.0F;
      float cur_radius;
      float pixel_scale = 1.0F / info->vertex_scale;
      float max_size = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting,
					     cSetting_sphere_point_max_size) * 3;
      int clamp_size_flag = (max_size >= 0.0F);

      short hasBegun = 0;
      if(!sp) {
	RepSpheresPrepPickingIfNoSphereGeometry(I, sphere_mode, &pixel_scale);
	glBegin(GL_POINTS);
	hasBegun = 1;
      }
      v = I->VC;
      c = I->NC;
      while(c--) {
	int skip = (p[1].index < 0);
	if(!skip) {
	  RepSpheresSetColorForPicking(I, pick, &i, &j, &p);
	  if(sp) {
	    RepSpheresRenderSphereGeometryForPicking(sp,  v + 4, v[7] );
	  } else {
	    RepSpheresRenderPointForPicking(I, v[7], v, sphere_mode, &last_radius, &cur_radius, pixel_scale, clamp_size_flag, max_size, &hasBegun);
	  }
	} else {
	  p++;
	}
	v += 8;
      }
      if(!sp) {
	glEnd();

	RepSpheresRenderEndOfPicking(sphere_mode);
      }
    }
    (*pick)[0].src.index = i;
  }
}

static void RepSphereRender(RepSphere * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V, *vn = I->VN;
  int c = I->N;
  SphereRec *sp = I->SP;
  float alpha;
  int ok = true;
  short use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) & 
                     SettingGetGlobal_b(G, cSetting_use_shaders);
  int sphere_mode = SettingGet_i(G, I->R.cs->Setting,
				 I->R.obj->Setting,
				 cSetting_sphere_mode);

  if (!ray) {
    switch (sphere_mode) {
    case 5:
#ifdef _PYMOL_ARB_SHADERS
      if (!sphereARBShaderPrg && G->HaveGUI && G->ValidContext) {
	sphereARBShaderPrg = CShaderPrg_NewARB(G, "sphere_arb", sphere_arb_vs, sphere_arb_fs);
      }
      if (!sphereARBShaderPrg)
#endif
      {
        PRINTFB(G, FB_ShaderMgr, FB_Warnings)
          " Warning: ARB shaders (sphere_mode=5) not supported.\n" ENDFB(G);
        sphere_mode = 9;
      }
      break;
    case -1:
      sphere_mode = 9;
    case 9:
      if (!use_shader || !CShaderMgr_ShaderPrgExists(G->ShaderMgr, "sphere")) {
        sphere_mode = 0;
      }
    }
  }
  alpha =
    SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_sphere_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  if(ray) {
    ok &= RepSphereRenderRay(I, info, alpha);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      RepSphereRenderPick(I, info, alpha, sphere_mode);
    } else {                    /* not pick */      
      if(!sp) {
        /* no sp -- we're rendering as points */
        v = I->VC;
        c = I->NC;

        I->LastVertexScale = info->vertex_scale;

          if((sphere_mode > 0) && (!info->line_lighting))
            glDisable(GL_LIGHTING);

          switch (sphere_mode) {
          case -1:
          case 0:              /* memory-efficient sphere rendering */
	    if (ok)
	      ok &= RenderSphereMode_Direct(G, I, info, c, &v, alpha, I->SSP);
            break;
          case 2:
          case 3:
          case 7:
          case 8:
	    RenderSphereMode_Sprites(G, I, info, sphere_mode, c, &v, &vn);
	    break;
          case 4:
	    RenderSphereMode_Points(G, I, info, c);
	    break;
	  case 5:          /* use vertex/fragment program */
#ifdef _PYMOL_ARB_SHADERS
	    RenderSphereMode_ARB(G, info, &v, c);
	    break;
#endif
	  case 9: // use GLSL shader
	    RenderSphereMode_9(G, I, info, &v, c);
	    break;
	  default:
	    RenderSphereMode_1_or_6(G, I, info, &v, &vn, c, alpha);
            break;
          }
          glEnable(GL_LIGHTING);
      } else {                  /* real spheres, drawn with triangles -- not points or impostors */
	ok &= RenderSphereMode_Direct(G, I, info, c, &v, alpha, I->SSP);
      }
    }
  }
}

int RepSphereSameVis(RepSphere * I, CoordSet * cs)
{
  int *lv, *lc;
  int a;
  AtomInfoType *ai;
  if(I->LastVisib && I->LastColor) {
    lv = I->LastVisib;
    lc = I->LastColor;

    for(a = 0; a < cs->NIndex; a++) {
      ai = cs->getAtomInfo(a);
      if(*(lv++) != GET_BIT(ai->visRep, cRepSphere)) {
        return false;
      }
      if(*(lc++) != ai->color) {
        return false;
      }
    }
  } else {
    return false;
  }
  return true;
}

static int RadiusOrder(float *list, int a, int b)
{
  return (list[a * 8 + 7] <= list[b * 8 + 7]);
}

static bool RepSphereDetermineAtomVisibility(PyMOLGlobals *G,
    AtomInfoType *ati1, int cartoon_side_chain_helper, int ribbon_side_chain_helper)
{
  if (!(ati1->flags & cAtomFlag_polymer))
    return true;

  bool sc_helper =
    (GET_BIT(ati1->visRep, cRepCartoon) &&
     AtomSettingGetWD(G, ati1, cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper)) ||
    (GET_BIT(ati1->visRep, cRepRibbon) &&
     AtomSettingGetWD(G, ati1, cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper));

  if (sc_helper) {
    int prot1 = ati1->protons;

    if(prot1 == cAN_N) {
      if(ati1->name == G->lex_const.N) {
        if(ati1->resn != G->lex_const.PRO)
	  return false;
      }
    } else if(prot1 == cAN_O) {
      if(ati1->name == G->lex_const.O)
	return false;
    } else if(prot1 == cAN_C) {
      if(ati1->name == G->lex_const.C)
	return false;
    }
  }
  return true;
}

static void RepSphereAddAtomVisInfoToStoredVC(RepSphere *I, ObjectMolecule *obj,
    CoordSet * cs, int state, float *varg, int a1, AtomInfoType *ati1, int a,
    int *mf, float sphere_scale, int sphere_color, float transp,
    int *variable_alpha, float sphere_add)
{
  PyMOLGlobals *G = cs->State.G;
  float at_transp = transp;
  int c1;
  float *v0, *vc;
  float *v = varg;

  float at_sphere_scale = AtomSettingGetWD(G, ati1, cSetting_sphere_scale, sphere_scale);
  int at_sphere_color = AtomSettingGetWD(G, ati1, cSetting_sphere_color, sphere_color);

  if(AtomSettingGetIfDefined(G, ati1, cSetting_sphere_transparency, &at_transp))
    *variable_alpha = true;
  
  if(I->R.P) {
    I->NP++;
    if(!ati1->masked) {
      I->R.P[I->NP].index = a1;
    } else {
      I->R.P[I->NP].index = -1;
    }
    I->R.P[I->NP].bond = -1;
  }
  
  *mf = true;
  I->NC++;
  if(at_sphere_color == -1)
    c1 = ati1->color;
  else
    c1 = at_sphere_color;
  v0 = cs->Coord + 3 * a;
  if(ColorCheckRamped(G, c1)) {
    ColorGetRamped(G, c1, v0, v, state);
    v += 3;
  } else {
    vc = ColorGet(G, c1);   /* save new color */
    *(v++) = *(vc++);
    *(v++) = *(vc++);
    *(v++) = *(vc++);
  }
  *(v++) = 1.0F - at_transp;
  
  *(v++) = *(v0++);         /* coordinate */
  *(v++) = *(v0++);
  *(v++) = *(v0++);
  *(v++) = obj->AtomInfo[a1].vdw * at_sphere_scale + sphere_add;    /* radius */
}

static int RepSphereComputeSphereNormals(RepSphere *I){
  PyMOLGlobals *G = I->R.G;
  float range = 6.0F;
  float *vc = I->VC;
  float *dot = G->Sphere->Sphere[1]->dot[0];
  int n_dot = G->Sphere->Sphere[1]->nDot;
  int nc = I->NC;
  int *active = NULL;
  float *v_tmp = NULL;
  int ok = true;
  int a;
  MapType *map = NULL;
  float *v0, *v;
  active = Alloc(int, 2 * n_dot);
  CHECKOK(ok, active);
  if (ok)
    v_tmp = Alloc(float, 3 * nc);
  CHECKOK(ok, v_tmp);
  if (ok) {
    float *src = vc + 4;
    float *dst = v_tmp;
    for(a = 0; a < nc; a++) {       /* create packed array of sphere centers */
      *(dst++) = *(src++);
      *(dst++) = *(src++);
      *(dst++) = *(src++);
      src += 5;
    }
    {
      map = MapNew(G, range, v_tmp, nc, NULL);
      CHECKOK(ok, map);
      if (ok)
	I->VN = Alloc(float, I->NC * 3);
      CHECKOK(ok, I->VN);
      if(ok && map && I->VN) {
	float dst;
	float *vv;
	int nbr_flag;
	int n_dot_active, *da;
	float cut_mult = -1.0F;
	float range2 = range * range;

	ok &= MapSetupExpress(map);
	if (ok){
	  v = vc + 4;
	  v0 = I->VN;
	  for(a = 1; a < n_dot; a++) {
	    float t_dot = dot_product3f(dot, dot + a * 3);
	    if(cut_mult < t_dot)
	      cut_mult = t_dot;
	  }
	}
	for(a = 0; ok && a < nc; a++) {
	  int h, k, l, b, i, j;
	  float v1[3];
	  nbr_flag = false;
	  MapLocus(map, v, &h, &k, &l);
	  da = active;
	  for(b = 0; b < n_dot; b++) {
	    *(da++) = b * 3;
	  }
	  n_dot_active = n_dot;
	  i = *(MapEStart(map, h, k, l));
	  if(i) {
	    j = map->EList[i++];
	    while(ok && j >= 0) {
	      if(j != a) {
		vv = v_tmp + 3 * j;
		if(within3fret(vv, v, range, range2, v1, &dst)) {
		  float cutoff = dst * cut_mult;
		  b = 0;
		  while(b < n_dot_active) {
		    vv = dot + active[b];
		    if(dot_product3f(v1, vv) > cutoff) {
		      n_dot_active--;
		      active[b] = active[n_dot_active];
		    }
		    b++;
		  }
		}
	      }
	      j = map->EList[i++];
	      ok &= !G->Interrupt;
	    }
	  }
	  if (ok){
	    if(!n_dot_active) {
	      v0[0] = 0.0F;
	      v0[1] = 0.0F;
	      v0[2] = 1.0F;
	    } else {
	      zero3f(v0);
	      b = 0;
	      while(b < n_dot_active) {
		vv = dot + active[b];
		add3f(vv, v0, v0);
		b++;
	      }
	      normalize3f(v0);
	    }
	    v += 8;
	    v0 += 3;
	  }
	  ok &= !G->Interrupt;
	}
      }
    }
  }
  MapFree(map);
  map = NULL;
  FreeP(v_tmp);
  map = NULL;
  FreeP(active);
  return ok;
}

static int RepSphereWriteSphereRecIntoArray(SphereRec *sp, int spheroidFlag,
    CoordSet * cs, float **varg, int a1, float *v0, float vdw,
    float spheroid_scale)
{
  PyMOLGlobals *G = cs->State.G;
  float *v = *varg;
  int b, *q, *s, c, ok = true;
  q = sp->Sequence;
  s = sp->StripLen;
  if(spheroidFlag) {
    for(b = 0; ok && b < sp->NStrip; b++) {
      float *sphLen = cs->Spheroid + (sp->nDot * a1);
      float *sphNorm = cs->SpheroidNormal + (3 * sp->nDot * a1);
      for(c = 0; c < (*s); c++) {
	float sphTmp, *sphTmpN = sphNorm + 3 * (*q);
	*(v++) = *(sphTmpN++);
	*(v++) = *(sphTmpN++);
	*(v++) = *(sphTmpN++);
	sphTmp = (*(sphLen + (*q))) * spheroid_scale;
	*(v++) = v0[0] + sphTmp * sp->dot[*q][0];       /* point */
	*(v++) = v0[1] + sphTmp * sp->dot[*q][1];
	*(v++) = v0[2] + sphTmp * sp->dot[*q][2];
	q++;
      }
      s++;
      ok &= !G->Interrupt;
    }
  } else {
    for(b = 0; ok && b < sp->NStrip; b++) {
      for(c = 0; ok && c < (*s); c++) {
	*(v++) = sp->dot[*q][0];        /* normal */
	*(v++) = sp->dot[*q][1];
	*(v++) = sp->dot[*q][2];
	*(v++) = v0[0] + vdw * sp->dot[*q][0];  /* point */
	*(v++) = v0[1] + vdw * sp->dot[*q][1];
	*(v++) = v0[2] + vdw * sp->dot[*q][2];
	q++;
	ok &= !G->Interrupt;
      }
      s++;
      ok &= !G->Interrupt;
    }
  }
  return ok;
}
/* 
 When Culling, triangles that are within other close spheres (i.e., within vdw plus solvent_radius),
 then the triangles are not rendered.  This saves geometry, but adds to calculating these overlaps.
 right now, it looks like these new strips are generated on the fly, this whole RepSphere.c file could be
 enhanced by CGOs (other than sphere_mode 9 which already uses them)

 sphere_add - the setting solvent_radius 
*/
static int RepSphereGenerateGeometryCullForSphere(SphereRec *sp,
    ObjectMolecule *obj, CoordSet * cs, float **varg, MapType *map,
    float vdw, float *v0, int *visFlag, int *marked, float sphere_scale,
    int a, float sphere_add, int *nt)
{
  PyMOLGlobals *G = cs->State.G;
  float *v = *varg;
  int ok = true;
  int b;
  int h, k, l, i, j, a2;
  int *q, *s, c;
  float v1[3];
  short restart;

  /* go through each vertex in the SphereRec, flag it to 0 if there are any coordinates
     that are within the vdw * sphere_scale + solvent radius */
  for(b = 0; ok && b < sp->nDot; b++) {       /* Sphere culling mode - more strips, but many fewer atoms */
    v1[0] = v0[0] + vdw * sp->dot[b][0];
    v1[1] = v0[1] + vdw * sp->dot[b][1];
    v1[2] = v0[2] + vdw * sp->dot[b][2];
    MapLocus(map, v1, &h, &k, &l);
    visFlag[b] = 1;
    i = *(MapEStart(map, h, k, l));
    if(i) {
      j = map->EList[i++];
      while(j >= 0) {
	a2 = cs->IdxToAtm[j];
	if(marked[a2]) {
	  AtomInfoType *ati2 = obj->AtomInfo + a2;

          float at2_sphere_scale = AtomSettingGetWD(G, ati2, cSetting_sphere_scale, sphere_scale);

	  if(j != a)
	    if(within3f(cs->Coord + 3 * j, v1,
			ati2->vdw * at2_sphere_scale + sphere_add)) {
	      /* do not render atoms if there are atoms within this distance */
	      visFlag[b] = 0;
	      break;
	    }
	}
	j = map->EList[i++];
      }
    }
    ok &= !G->Interrupt;
  }
  q = sp->Sequence;
  s = sp->StripLen;
  for(b = 0; ok && b < sp->NStrip; b++) {
    /* this is an attempt to fill in *some* of the cracks
     * by checking to see if the center of the triangle is visible 
     * IMHO - the increase in framerates is worth missing a triangle
     * here or there, and the user can always turn off sphere culling */
    q += 2;
    for(c = 2; c < (*s); c++) {
      float v1[3];
      int q0, q1, q2, flag;
      q0 = *q;
      q1 = *(q - 1);
      q2 = *(q - 2);
      
      if((!visFlag[q0]) && (!visFlag[q1]) && (!visFlag[q2]))
	v1[0] = v0[0] + vdw * sp->dot[q0][0];
      v1[1] = v0[1] + vdw * sp->dot[q0][1];
      v1[2] = v0[2] + vdw * sp->dot[q0][2];
      
      v1[0] += v0[0] + vdw * sp->dot[q1][0];
      v1[1] += v0[1] + vdw * sp->dot[q1][1];
      v1[2] += v0[2] + vdw * sp->dot[q1][2];
      
      v1[0] += v0[0] + vdw * sp->dot[q2][0];
      v1[1] += v0[1] + vdw * sp->dot[q2][1];
      v1[2] += v0[2] + vdw * sp->dot[q2][2];
      
      v1[0] /= 3;
      v1[1] /= 3;
      v1[2] /= 3;
      
      flag = true;
      i = *(MapEStart(map, h, k, l));
      if(i) {
	j = map->EList[i++];
	while(j >= 0) {
	  a2 = cs->IdxToAtm[j];
	  if(marked[a2]) {
	    if(j != a)
	      if(within3f
		 (cs->Coord + 3 * j, v1,
		  cs->Obj->AtomInfo[a2].vdw * sphere_scale + sphere_add)) {
		flag = false;
		break;
	      }
	  }
	  j = map->EList[i++];
	}
      }
      if(flag) {
	visFlag[q0] = 1;
	visFlag[q1] = 1;
	visFlag[q2] = 1;
      }
      q++;
    }
    s++;
    ok &= !G->Interrupt;
  }
  
  *(nt) = 0;            /* how many passes through the triangle renderer? */
  q = sp->Sequence;
  s = sp->StripLen;
  
  for(b = 0; ok && b < sp->NStrip; b++) {
    restart = true;      /* startin a new strip */
    for(c = 0; c < (*s); c++) {
      if(c > 1) {       /* on third vertex or better */
	int q0, q1, q2;
	q0 = *q;        /* get the indices of the triangle in this strip */
	q1 = *(q - 1);
	q2 = *(q - 2);
	if(visFlag[q0] || (visFlag[q1]) || (visFlag[q2])) {     /* visible? */
	  *(v++) = restart;     /* store continuing string flag */
	  
	  if(restart) { /* not continuing...this is a new strip */
	    if(c & 0x1) /* make sure strip starts off "right" */
	      *(v - 1) = 2.0;
	    *(v++) = sp->dot[q2][0];    /* normal */
	    *(v++) = sp->dot[q2][1];
	    *(v++) = sp->dot[q2][2];
	    *(v++) = v0[0] + vdw * sp->dot[q2][0];      /* point */
	    *(v++) = v0[1] + vdw * sp->dot[q2][1];
	    *(v++) = v0[2] + vdw * sp->dot[q2][2];
	    *(v++) = sp->dot[q1][0];    /* normal */
	    *(v++) = sp->dot[q1][1];
	    *(v++) = sp->dot[q1][2];
	    *(v++) = v0[0] + vdw * sp->dot[q1][0];      /* point */
	    *(v++) = v0[1] + vdw * sp->dot[q1][1];
	    *(v++) = v0[2] + vdw * sp->dot[q1][2];
	    *(v++) = sp->dot[q0][0];    /* normal */
	    *(v++) = sp->dot[q0][1];
	    *(v++) = sp->dot[q0][2];
	    *(v++) = v0[0] + vdw * sp->dot[q0][0];      /* point */
	    *(v++) = v0[1] + vdw * sp->dot[q0][1];
	    *(v++) = v0[2] + vdw * sp->dot[q0][2];
	  } else {      /* continue strip */
	    *(v++) = sp->dot[q0][0];    /* normal */
	    *(v++) = sp->dot[q0][1];
	    *(v++) = sp->dot[q0][2];
	    *(v++) = v0[0] + vdw * sp->dot[q0][0];      /* point */
	    *(v++) = v0[1] + vdw * sp->dot[q0][1];
	    *(v++) = v0[2] + vdw * sp->dot[q0][2];
	  }
	  restart = false;
	  (*nt)++;
	} else {
	  restart = true;        /* next triangle is a new strip */
	}
      }
      q++;
    }
    s++;
    ok &= !G->Interrupt;
  }
  return ok;
}

static int RepSphereGenerateGeometryForSphere(RepSphere *I, ObjectMolecule *obj,
    CoordSet * cs, int state, int a1, AtomInfoType *ati1, int a,
    float sphere_scale, int sphere_color, float spheroid_scale, float transp,
    int *variable_alpha, float sphere_add, int spheroidFlag, SphereRec *sp,
    int *visFlag, int *marked, MapType *map, int *nt, float **varg)
{
  PyMOLGlobals *G = cs->State.G;
  float *v = *varg;
  float at_transp = transp;
  int ok = true;
  int c1;
  float *v0;
  float vdw;

  float at_sphere_scale = AtomSettingGetWD(G, ati1, cSetting_sphere_scale, sphere_scale);
  int at_sphere_color = AtomSettingGetWD(G, ati1, cSetting_sphere_color, sphere_color);

  if(AtomSettingGetIfDefined(G, ati1, cSetting_sphere_transparency, &at_transp))
    (*variable_alpha) = true;
  
  if(at_sphere_color == -1)
    c1 = ati1->color;
  else
    c1 = at_sphere_color;
  v0 = cs->Coord + 3 * a;
  vdw = ati1->vdw * at_sphere_scale + sphere_add;
  if(ColorCheckRamped(G, c1)) {
    ColorGetRamped(G, c1, v0, v, state);
    v += 3;
  } else {
    float *vc = ColorGet(G, c1);
    *(v++) = *(vc++);
    *(v++) = *(vc++);
    *(v++) = *(vc++);
  }
  *(v++) = 1.0F - at_transp;      /* alpha */
  if(I->cullFlag && (!spheroidFlag) && (sp)) {
    ok &= RepSphereGenerateGeometryCullForSphere(sp, obj, cs, &v, map, vdw, v0, visFlag, marked, sphere_scale, a, sphere_add, nt);
  } else if(sp) {
    ok &= RepSphereWriteSphereRecIntoArray(sp, spheroidFlag, cs, &v, a1, v0, vdw, spheroid_scale);
  } else if (ok) {                /* if sp is null, then we're simply drawing points */
    *(v++) = v0[0];
    *(v++) = v0[1];
    *(v++) = v0[2];
  }
  I->N++;
  *varg = v;
  return ok;
}


Rep *RepSphereNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int ok = true;
  int a, a1;
  float *v;
  int *lv, *lc;
  SphereRec *sp = G->Sphere->Sphere[0];
  int sphere_quality, *nt;
  int *visFlag = NULL;
  MapType *map = NULL;
  AtomInfoType *ai2;
  int spheroidFlag = false;
  float spheroid_scale;
  float sphere_scale, sphere_add = 0.f;
  int sphere_color;
  int *map_flag = NULL, *mf;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 0;
  AtomInfoType *ati1;
  int vis_flag;
  int sphere_mode;
  int *marked = NULL;
  float transp;
  int variable_alpha = false;
#ifdef _this_code_is_not_used
  float vv0[3], vv1[3], vv2[3];
  float tn[3], vt1[3], vt2[3], xtn[3], *tn0, *tn1, *tn2;
#endif
  int draw_mode = SettingGetGlobal_i(G, cSetting_draw_mode);
  int draw_quality = (((draw_mode == 1) || (draw_mode == -2) || (draw_mode == 2)));
  short use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) & 
                     SettingGetGlobal_b(G, cSetting_use_shaders);

  // skip if not visible
  if(!cs->hasRep(cRepSphereBit))
    return NULL;

  OOCalloc(G, RepSphere);
  CHECKOK(ok, I);
  if (!ok)
    return NULL;
  obj = cs->Obj;

  marked = Calloc(int, obj->NAtom);
  CHECKOK(ok, marked);
  if (ok)
    RepInit(G, &I->R);
  I->shaderCGO = NULL;

  if (ok){
    sphere_quality = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_quality);
    sphere_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_mode);
    if (!use_shader && (sphere_mode == 5 || sphere_mode == 9)){
      sphere_mode = 0;
    }
    if(sphere_mode > 0)
      sphere_quality = -1;
    
    if(sphere_quality < 0) {
      sp = NULL;
    } else {
      if(draw_quality && (sphere_quality < 3))
	sphere_quality = 3;
      if(sphere_quality > 4)
	sphere_quality = 4;
      sp = G->Sphere->Sphere[sphere_quality];
    }
  }
  if (ok){
    sphere_color =
      SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_color);
    cartoon_side_chain_helper =
      SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_side_chain_helper);
    ribbon_side_chain_helper =
      SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_side_chain_helper);
    transp = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_transparency);
    spheroid_scale =
      SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_spheroid_scale);
    if(sp && spheroid_scale && cs->Spheroid)
      spheroidFlag = 1;
    else
      spheroidFlag = 0;
    sphere_scale = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_scale);
  }

  if (ok){
    I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepSphereRender;
    I->R.fFree = (void (*)(struct Rep *)) RepSphereFree;
    I->R.fSameVis = (int (*)(struct Rep *, struct CoordSet *)) RepSphereSameVis;
    I->LastVertexScale = -1.0F;
    I->R.obj = (CObject *) obj;
    I->R.cs = cs;
    I->R.context.object = (void *) obj;
    I->R.context.state = state;
  }
  /* raytracing primitives */

  if (ok)
    I->VC = Alloc(float, cs->NIndex * 8);
  CHECKOK(ok, I->VC);
  if (ok){
    I->NC = 0;
    map_flag = Calloc(int, cs->NIndex);
  }
  CHECKOK(ok, map_flag);

  if (ok){
    I->NT = NULL;
    nt = NULL;
    v = I->VC;
    mf = map_flag;
    
    if(SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_solvent)) { /* are we generating a solvent surface? */
      sphere_add = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_solvent_radius);       /* if so, get solvent radius */
    }
    
    if(SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_pickable)) {
      I->R.P = Alloc(Pickable, cs->NIndex + 1);
      CHECKOK(ok, I->R.P);
    }
  }

  I->spheroidFlag = spheroidFlag;
  for(a = 0; ok && a < cs->NIndex; a++) {
    a1 = cs->IdxToAtm[a];
    ati1 = obj->AtomInfo + a1;
    /* store temporary visibility information */
    marked[a1] = GET_BIT(ati1->visRep,cRepSphere) &&
        RepSphereDetermineAtomVisibility(G, ati1, 
                                         cartoon_side_chain_helper, ribbon_side_chain_helper);
    if(marked[a1]) {
      RepSphereAddAtomVisInfoToStoredVC(I, obj, cs, state, v, a1, ati1, a, mf, sphere_scale, sphere_color, transp, &variable_alpha, sphere_add);
      v += 8;
    }
    mf++;
    ok &= !G->Interrupt;
  }
  if (ok){
    I->VariableAlphaFlag = variable_alpha;
    if(I->NC)
      I->VC = ReallocForSure(I->VC, float, (v - I->VC));
    else
      I->VC = ReallocForSure(I->VC, float, 1);
    CHECKOK(ok, I->VC);
    if(ok && I->R.P) {
      I->R.P = Realloc(I->R.P, Pickable, I->NP + 1);
      CHECKOK(ok, I->R.P);
      if (ok)
	I->R.P[0].index = I->NP;
    }
  }
  if (ok){
    if(variable_alpha)
      I->cullFlag = false;
    else {
      I->cullFlag = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cull_spheres);
      if(I->cullFlag < 0) { /* if negative, then set if more than one state */
	I->cullFlag = !(obj->NCSet > 1);
      }
    }
    if(draw_quality)
      I->cullFlag = false;
    if(spheroidFlag || (!sp))
      I->cullFlag = false;
    if((I->cullFlag < 2) &&
       (SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_sculpting)))
      /* optimize build-time performance when sculpting */
      I->cullFlag = false;
    if((I->cullFlag < 2) && (SettingGetGlobal_f(G, cSetting_roving_spheres) != 0.0F))
      I->cullFlag = false;
    if(sp && ((I->cullFlag < 2) ) ) {
      /* don't cull unless forced */
      I->SSP = sp;
      sp = NULL;
    }
  }
  if(!sp) {                     /* if sp==null, then we're drawing a point-based sphere rep */
    /* sort the vertices by radius */
    if(ok && I->NC && I->VC && (!spheroidFlag) && (sphere_mode != 1)) {
      int *ix = NULL;
      float *vc_tmp = NULL;
      Pickable *pk_tmp = NULL;
      int a;
      ix = Alloc(int, I->NC);
      CHECKOK(ok, ix);
      if (ok)
	vc_tmp = Alloc(float, I->NC * 8);
      CHECKOK(ok, vc_tmp);
      if (ok)
	pk_tmp = Alloc(Pickable, I->NP + 1);
      CHECKOK(ok, pk_tmp);

      if(ok && vc_tmp && pk_tmp && ix) {
        UtilCopyMem(vc_tmp, I->VC, sizeof(float) * 8 * I->NC);
        if (I->R.P)
          UtilCopyMem(pk_tmp, I->R.P, sizeof(Pickable) * (I->NP + 1));

        UtilSortIndex(I->NC, I->VC, ix, (UtilOrderFn *) RadiusOrder);

        if (I->R.P)
          UtilCopyMem(I->R.P, pk_tmp, sizeof(Pickable));
        for(a = 0; a < I->NC; a++) {
          UtilCopyMem(I->VC + (a * 8), vc_tmp + (8 * ix[a]), sizeof(float) * 8);
          if (I->R.P)
            UtilCopyMem(I->R.P + (a + 1), pk_tmp + ix[a] + 1, sizeof(Pickable));
        }
      }
      FreeP(vc_tmp);
      FreeP(ix);
      FreeP(pk_tmp);
    }
    if(ok && (sphere_mode >= 6) && (sphere_mode < 9) && I->NC) {
      /* compute sphere normals in VN to approximate a surface */
      ok = RepSphereComputeSphereNormals(I);
    }
    I->cullFlag = false;
    I->V = NULL;
    I->NT = NULL;
    I->N = 0;
    I->SP = NULL;
  } else {
    /* if sp, drawing geometry-based sphere rep */
    if(I->cullFlag && sp) {
      if (ok)
	I->V = Alloc(float, I->NC * (sp->NVertTot * 31));    /* double check 31 */
      CHECKOK(ok, I->V);
      if (ok)
	I->NT = Alloc(int, cs->NIndex);
      CHECKOK(ok, I->NT);
      if (ok)
	visFlag = Alloc(int, sp->nDot);
      CHECKOK(ok, visFlag);
      /* hmm...need to compute max(sphere_scale) for all atoms... */
      if (ok)
	map =
	  MapNewFlagged(G, MAX_VDW * sphere_scale + sphere_add, cs->Coord, cs->NIndex, NULL,
			map_flag);
      CHECKOK(ok, map);
      if (ok)
	ok &= MapSetupExpress(map);
    } else if (ok){
      if(sp)
        I->V = Alloc(float, I->NC * (4 + sp->NVertTot * 6));
      else
        I->V = Alloc(float, I->NC * 7);    /* one color, one alpha, one vertex per spheres */
      CHECKOK(ok, I->V);
    }
    /* rendering primitives */
    I->N = 0;
    I->SP = sp;
    v = I->V;
    nt = I->NT;

    for(a = 0; ok && a < cs->NIndex; a++) {
      a1 = cs->IdxToAtm[a];
      ati1 = obj->AtomInfo + a1;
      vis_flag = marked[a1];

      /* don't show backbone atoms if side_chain_helper is on */

      if(vis_flag) {
	ok &= RepSphereGenerateGeometryForSphere(I, obj, cs, state, a1, ati1, a, sphere_scale, sphere_color, spheroid_scale, transp, &variable_alpha, sphere_add, spheroidFlag, sp, visFlag, marked, map, nt, &v);
	if(nt)
	  nt++;
	ok &= !G->Interrupt;
	if (!ok)
	  break;
      }
    }
  }
  if(ok) {
    if(!I->LastVisib)
      I->LastVisib = Alloc(int, cs->NIndex);
    CHECKOK(ok, I->LastVisib);
    if(ok && !I->LastColor)
      I->LastColor = Alloc(int, cs->NIndex);
    CHECKOK(ok, I->LastColor);
    if (ok){
      lv = I->LastVisib;
      lc = I->LastColor;
      obj = cs->Obj;
      ai2 = obj->AtomInfo;
      if(sphere_color == -1){
	for(a = 0; a < cs->NIndex; a++) {
          int at = cs->IdxToAtm[a];
	  *(lv++) = marked[at];
	  *(lc++) = (ai2 + at)->color;
	}
      } else {
	for(a = 0; a < cs->NIndex; a++) {
	  *(lv++) = marked[cs->IdxToAtm[a]];
	  *(lc++) = sphere_color;
	}
      }
    }
  }

  if(ok && I->V) {
    if(I->N) {
      I->V = ReallocForSure(I->V, float, (v - I->V));
      CHECKOK(ok, I->V);
      if(ok && I->NT){
        I->NT = ReallocForSure(I->NT, int, (nt - I->NT));
	CHECKOK(ok, I->NT);
      }
    } else {
      I->V = ReallocForSure(I->V, float, 1);
      CHECKOK(ok, I->V);
      if(ok && I->NT){
        I->NT = ReallocForSure(I->NT, int, 1);
	CHECKOK(ok, I->NT);
      }
    }
  }
  FreeP(marked);
  FreeP(visFlag);
  FreeP(map_flag);
  if(map)
    MapFree(map);
  if(!ok) {
    RepSphereFree(I);
    I = NULL;
  }
  return (Rep *) I;
}
