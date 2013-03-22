
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

#ifdef NT
#undef NT
#endif

typedef struct RepSphere {
  Rep R;
  float *V;                     /* triangle vertices (if any) */
  float *VC;                    /* sphere centers, colors, alpha, and radii */
  float *VN;                    /* normals (if any) */
  SphereRec *SP, *SSP;
  int *NT;
  int N, NC, NP;
  int cullFlag, spheroidFlag;
  int *LastVisib;
  int *LastColor;
  float LastVertexScale;
  int VariableAlphaFlag;
  CGO *shaderCGO;
} RepSphere;

#include"ObjectMolecule.h"

#ifdef _PYMOL_OPENGL_SHADERS

#ifndef GL_FRAGMENT_PROGRAM_ARB
#define GL_FRAGMENT_PROGRAM_ARB                         0x8804
#endif


#include "ShaderMgr.h"

/* END PROPRIETARY CODE SEGMENT */


/* NOTE -- right now this shader program only runs in perspective mode */

#include "ShaderText.h"

/*
  normal depth routine...does not work!  why?
  "#MAD_SAT fogFactor.x, fogInfo.x, fragment.texcoord.w, fogInfo.y;\n",
  "#LRP color.xyz, fogFactor.x, color, fogColor;\n",
*/

#endif

void RepSphereFree(RepSphere * I);
int RepSphereSameVis(RepSphere * I, CoordSet * cs);

void RepSphereFree(RepSphere * I)
{
#ifdef _PYMOL_OPENGL_SHADERS
  if(I->R.G->HaveGUI && I->R.G->ValidContext) {
  }
#endif
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

#ifdef _PYMOL_OPENGL_SHADERS

/* MULTI-INSTSANCE TODO:  isn't this a conflict? */
static CShaderPrg *sphereARBShaderPrg = NULL;
#endif
#ifdef _PYMOL_GL_DRAWARRAYS
void RepSphereRenderImmediatePointsES(PyMOLGlobals *G, int sphere_mode, AtomInfoType *atomInfo, int starta, int enda, int nverts, int *i2aarg, float *varg){
  int a, *i2a = i2aarg, pl = 0, plc = 0;
  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
  float *cur_color, *v = varg;
  int last_color = -1;

  for(a = starta; a < enda; a++) {
    AtomInfoType *ai = atomInfo + *(i2a++);
    if(ai->visRep[cRepSphere]) {
      int c = ai->color;
      if(c != last_color) {
	last_color = c;
	cur_color = ColorGet(G, c);
      }
      switch (sphere_mode) {
      case 1:
      case 2:
      case 3:
	colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];
      }
    }
    v += 3;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
  glColorPointer(4, GL_FLOAT, 0, colorVals);
  glDrawArrays(GL_LINES, 0, nverts);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  DEALLOCATE_ARRAY(ptsVals)
  DEALLOCATE_ARRAY(colorVals)
}

void RepSphereRenderImmediateMode4PointsES(PyMOLGlobals *G, AtomInfoType *atomInfo, int starta, int enda, int nverts, int *i2aarg, float *varg, 
					   float zz_factor, float s_factor, float *dim_add, float _1){
  int a, *i2a = i2aarg, pl = 0, plc = 0;
  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
  float *v = varg;
  for(a = starta; a < enda; a++) {
    AtomInfoType *ai = atomInfo + *(i2a++);
    if(ai->visRep[cRepSphere]) {
      float *vc = ColorGet(G, ai->color);
      float r = zz_factor * vc[0] + s_factor;
      float g = zz_factor * vc[1] + s_factor;
      float b = zz_factor * vc[2] + s_factor;
      colorVals[plc++] = r > _1 ? _1 : r; colorVals[plc++] = g > _1 ? _1 : g; colorVals[plc++] = b > _1 ? _1 : b; colorVals[plc++] = 1.f;
      ptsVals[pl++] = v[0] + dim_add[0]; ptsVals[pl++] = v[1] + dim_add[1]; ptsVals[pl++] = v[2] + dim_add[2];
    }
    v += 3;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
  glColorPointer(4, GL_FLOAT, 0, colorVals);
  glDrawArrays(GL_POINTS, 0, nverts);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  DEALLOCATE_ARRAY(ptsVals)
  DEALLOCATE_ARRAY(colorVals)
}

void RepSphereRenderImmediateQuadsES(PyMOLGlobals *G, AtomInfoType *atomInfo, int starta, int enda, int nverts, int *i2aarg, 
				     float *varg, float *_00, float *_10, float *_11, float *_01, float sphere_scale, 
				     float *fog_info, float *m, float cur_radius, float z_front, float z_back, float _1, float cutoff, float m_cutoff){
  int a, *i2a = i2aarg;
  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
  ALLOCATE_ARRAY(GLfloat,texVals,nverts*2)
  int pl = 0, plt = 0, plc = 0;
  float *cur_color, *v = varg;
  AtomInfoType *fai = atomInfo + *(i2a);
  register float v0, v1, v2, nv0, nv1, nv3, v3;
  v3 = fai->vdw * sphere_scale;
#ifndef _PYMOL_GL_DRAWARRAYS
  glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
			     0, 0.0F, 0.0F, v3, 0.0F);
  glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
			     0, fog_info[0], fog_info[1], 0.0F, 0.0F);
#endif
  for(a = starta; a < enda; a++) {
    AtomInfoType *ai = atomInfo + *(i2a++);
    if(ai->visRep[cRepSphere]) {
      v0 = v[0];
      v1 = v[1];
      v2 = v[2];
      nv3 = m[3] * v0 + m[7] * v1 + m[11] * v2 + m[15];       /* compute Wc */
      if(((nv3 - cur_radius) > z_front) && (nv3 < z_back)) {  /* is it within the viewing volume? */
	nv0 = m[0] * v0 + m[4] * v1 + m[8] * v2 + m[12];
	nv3 = _1 / nv3;
	nv1 = m[1] * v0 + m[5] * v1 + m[9] * v2 + m[13];
	nv0 *= nv3;
	nv1 *= nv3;
	if((nv0 < cutoff) && (nv0 > m_cutoff) &&
	   (nv1 < cutoff) && (nv1 > m_cutoff)) {
	  cur_color = ColorGet(G, ai->color);
	  colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	  texVals[plt++] = _00[0]; texVals[plt++] = _00[1];
	  ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];

	  colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	  texVals[plt++] = _10[0]; texVals[plt++] = _10[1];
	  ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];

	  colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	  texVals[plt++] = _11[0]; texVals[plt++] = _11[1];
	  ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];


	  colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	  texVals[plt++] = _11[0]; texVals[plt++] = _11[1];
	  ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];

	  colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	  texVals[plt++] = _01[0]; texVals[plt++] = _01[1];
	  ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];

	  colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	  texVals[plt++] = _00[0]; texVals[plt++] = _00[1];
	  ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];
	}
      }
    }
    v += 3;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
  glColorPointer(4, GL_FLOAT, 0, colorVals);
  glTexCoordPointer(2, GL_FLOAT, 0, texVals);
  glDrawArrays(GL_TRIANGLES, 0, nverts);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  DEALLOCATE_ARRAY(ptsVals)
  DEALLOCATE_ARRAY(colorVals)
  DEALLOCATE_ARRAY(texVals)
}
#endif

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
      register float pixel_scale = 1.0F / info->vertex_scale;

      switch (sphere_mode) {
      case 2:
        glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_ALPHA_TEST);
        pixel_scale *= 1.4F;
        glPointSize(1.0F);
        break;
      case 3:
        glEnable(GL_POINT_SMOOTH);
        glAlphaFunc(GL_GREATER, 0.5F);
        glEnable(GL_ALPHA_TEST);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        glPointSize(1.0F);
        pixel_scale *= 2.0F;
        break;
      case 4:
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_ALPHA_TEST);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        glPointSize(1.0F);
        pixel_scale *= 2.0F;
        break;
      default:
        glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
        glDisable(GL_POINT_SMOOTH);
        glDisable(GL_ALPHA_TEST);
        glPointSize(SettingGet_f
                    (G, cs->Setting, obj->Obj.Setting, cSetting_sphere_point_size));
        break;
      }

#ifdef _PYMOL_OPENGL_SHADERS
      if(sphere_mode == 5) {
	if (!sphereARBShaderPrg){
	  sphereARBShaderPrg = CShaderPrg_NewARB(G, "sphere_arb", sphere_arb_vs, sphere_arb_fs);
	}
        if(sphereARBShaderPrg){
          float fog_info[3];
          float _00[2] = { 0.0F, 0.0F };
          float _01[2] = { 0.0F, 1.0F };
          float _11[2] = { 1.0F, 1.0F };
          float _10[2] = { 1.0F, 0.0F };
          const float _1 = 1.0F;
          register float v0, v1, v2, nv0, nv1, nv2, nv3, v3;
          register float *m = info->pmv_matrix;
          register float cutoff = 1.2F;
          register float m_cutoff = -cutoff;
          register float z_front, z_back;
          /* compute -Ze = (Wc) of fog start */
          nv3 =
            (info->front +
             (info->back - info->front) * SettingGetGlobal_f(G, cSetting_fog_start));
          /* compute Zc of fog start using std. perspective transformation */
          nv2 =
            (nv3 * (info->back + info->front) -
             2 * (info->back * info->front)) / (info->back - info->front);
          /* compute Zc/Wc to get normalized depth coordinate of fog start */
          nv0 = (nv2 / nv3);
          fog_info[0] = (nv0 * 0.5) + 0.5;

          fog_info[1] = 1.0F / (1.0 - fog_info[0]);     /* effective range of fog */

          z_front = info->stereo_front;
          z_back = info->back + ((info->back + info->front) * 0.25);

	  CShaderPrg_Enable_SphereShaderARB(G);

          glNormal3fv(info->view_normal);
#ifdef _PYMOL_GL_DRAWARRAYS
	  {
            float last_radius = -1.0F, cur_radius;
            int a;
            int nIndex = cs->NIndex;
            AtomInfoType *atomInfo = obj->AtomInfo;
            int *i2a = cs->IdxToAtm;
            float *v = cs->Coord;
	    int nverts = 0;
	    int starta = 0, *starti2a = i2a;
	    float *startv = v;
            for(a = 0; a < nIndex; a++) {
              AtomInfoType *ai = atomInfo + *(i2a++);
              if(ai->visRep[cRepSphere]) {
                repActive = true;
                v3 = ai->vdw * sphere_scale;
                v0 = v[0];
                v1 = v[1];
                v2 = v[2];
                if(last_radius != (cur_radius = v3)) {
		  if (nverts>0){
		    RepSphereRenderImmediateQuadsES(G, atomInfo, starta, a, nverts, starti2a, startv, _00, _10, _11, _01, 
						    sphere_scale, fog_info, m, last_radius, z_front, z_back, _1, cutoff, m_cutoff);
		    starta = a;
		    nverts = 0;
		    starti2a = i2a;
		    startv = v;
		  }
                  last_radius = cur_radius;
                }
                /*  MatrixTransformC44f4f(info->pmv_matrix, v, nv); */
                nv3 = m[3] * v0 + m[7] * v1 + m[11] * v2 + m[15];       /* compute Wc */
                if(((nv3 - cur_radius) > z_front) && (nv3 < z_back)) {  /* is it within the viewing volume? */
                  nv0 = m[0] * v0 + m[4] * v1 + m[8] * v2 + m[12];
                  nv3 = _1 / nv3;
                  nv1 = m[1] * v0 + m[5] * v1 + m[9] * v2 + m[13];
                  nv0 *= nv3;
                  nv1 *= nv3;
                  if((nv0 < cutoff) && (nv0 > m_cutoff) &&
                     (nv1 < cutoff) && (nv1 > m_cutoff)) {
		    nverts += 6;
                  }
                }
              }
              v += 3;
            }
	    if (nverts>0){
	      RepSphereRenderImmediateQuadsES(G, atomInfo, starta, a, nverts, starti2a, startv, _00, _10, _11, _01, 
					      sphere_scale, fog_info, m, last_radius, z_front, z_back, _1, cutoff, m_cutoff);
	    }
#else
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
              if(ai->visRep[cRepSphere]) {
                repActive = true;
                v3 = ai->vdw * sphere_scale;
                v0 = v[0];
                v1 = v[1];
                v2 = v[2];
                if(last_radius != (cur_radius = v3)) {
                  glEnd();
                  glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
                                             0, 0.0F, 0.0F, v3, 0.0F);
                  glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
                                             0, fog_info[0], fog_info[1], 0.0F, 0.0F);
                  glBegin(GL_QUADS);
                  last_radius = cur_radius;
                }
                /*  MatrixTransformC44f4f(info->pmv_matrix, v, nv); */
                nv3 = m[3] * v0 + m[7] * v1 + m[11] * v2 + m[15];       /* compute Wc */
                if(((nv3 - cur_radius) > z_front) && (nv3 < z_back)) {  /* is it within the viewing volume? */
                  nv0 = m[0] * v0 + m[4] * v1 + m[8] * v2 + m[12];
                  nv3 = _1 / nv3;
                  nv1 = m[1] * v0 + m[5] * v1 + m[9] * v2 + m[13];
                  nv0 *= nv3;
                  nv1 *= nv3;
                  if((nv0 < cutoff) && (nv0 > m_cutoff) &&
                     (nv1 < cutoff) && (nv1 > m_cutoff)) {
                    glColor3fv(ColorGet(G, ai->color));
                    glTexCoord2fv(_00);
                    glVertex3fv(v);
                    glTexCoord2fv(_10);
                    glVertex3fv(v);
                    glTexCoord2fv(_11);
                    glVertex3fv(v);
                    glTexCoord2fv(_01);
                    glVertex3fv(v);
                  }
                }
              }
              v += 3;
            }
            glEnd();
#endif
          }
	  CShaderPrg_DisableARB(sphereARBShaderPrg);
        }
      } else
#endif
      if(sphere_mode == 4) {
        int repeat = true;
        const float _1 = 1.0F;
        const float _2 = 2.0F;
        float x_add = 0.0F, y_add = 0.0F, z_add = 0.0F;
#ifdef _PYMOL_GL_DRAWARRAYS
	float dim_add[3] = { 0.0F, 0.0F, 0.0F };
#endif
        float z_factor = 0.0F, r_factor = 1.0F;
        float s_factor = 0.0F;
        int pass = 0;
        register float max_size = SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                                               cSetting_sphere_point_max_size);
        register int clamp_size_flag = (max_size >= 0.0F);

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
#ifdef _PYMOL_GL_DRAWARRAYS
	  {
	    int *starti2a = i2a, starta = 0, nverts = 0;
	    float *startv = v;  
	    for(a = 0; a < nIndex; a++) {
	      AtomInfoType *ai = atomInfo + *(i2a++);
	      if(ai->visRep[cRepSphere]) {
		float cur_radius = ai->vdw;
		repActive = true;
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
		    if (nverts>0){
		      RepSphereRenderImmediateMode4PointsES(G, atomInfo, starta, a, nverts, starti2a, startv, zz_factor, s_factor, dim_add, _1);
		      nverts = 0;
		      starta = a;
		      startv = v;
		      starti2a = i2a;
		    }
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
		  }
		  x_add = z_factor * clamp_radius * info->view_normal[0];
		  y_add = z_factor * clamp_radius * info->view_normal[1];
		  z_add = z_factor * clamp_radius * info->view_normal[2];
		  dim_add[0] = x_add; dim_add[1] = y_add; dim_add[2] = z_add; 
		  last_radius = cur_radius;
		  last_size = size;
		}
		nverts++;
	      }
	      v += 3;
	    }
	    if (nverts>0){
	      RepSphereRenderImmediateMode4PointsES(G, atomInfo, starta, a, nverts, starti2a, startv, zz_factor, s_factor, dim_add, _1);
	    }
	  }
#else
          glBegin(GL_POINTS);
	  
          for(a = 0; a < nIndex; a++) {
            AtomInfoType *ai = atomInfo + *(i2a++);
            if(ai->visRep[cRepSphere]) {
              float cur_radius = ai->vdw;
              repActive = true;

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
#endif
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
      } else {
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

#ifdef _PYMOL_GL_DRAWARRAYS
	{
	  int *starti2a = i2a, starta = 0, nverts = 0;
	  float *startv = v;
	  (void)last_color;
	  for(a = 0; a < nIndex; a++) {
	    AtomInfoType *ai = atomInfo + *(i2a++);
	    if(ai->visRep[cRepSphere]) {
	      repActive = true;
	      switch (sphere_mode) {
	      case 1:
		nverts++;
		break;
	      case 2:
	      case 3:
		{
		  float cur_radius = ai->vdw * pixel_scale;
		  if(last_radius != cur_radius) {
		    glPointSize(last_radius);
		    if (nverts>0){
		      RepSphereRenderImmediatePointsES(G, sphere_mode, atomInfo, starta, a, nverts, starti2a, startv);
		    }
		    if(clamp_size_flag)
		      if(cur_radius > max_radius)
			cur_radius = max_radius;
		    glPointSize(cur_radius);
		    last_radius = cur_radius;
		    starta = a;
		    startv = v;
		    nverts = 0;
		  }
		  nverts++;
		}
		break;
	      }
	    }
	    v += 3;
	  }
	  if (nverts>0){
	    RepSphereRenderImmediatePointsES(G, sphere_mode, atomInfo, starta, a, nverts, starti2a, startv);
	  }
	}
#else
        glBegin(GL_POINTS);
        for(a = 0; a < nIndex; a++) {
          AtomInfoType *ai = atomInfo + *(i2a++);
          if(ai->visRep[cRepSphere]) {
            int c = ai->color;
            repActive = true;
            if(c != last_color) {
              last_color = c;
              glColor3fv(ColorGet(G, c));
            }
            switch (sphere_mode) {
            case 1:
              glVertex3fv(v);
              break;
            case 2:
            case 3:
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
#endif
        glEnable(GL_LIGHTING);

        if(sphere_mode == 3) {
          glDisable(GL_POINT_SMOOTH);
          glAlphaFunc(GL_GREATER, 0.05F);
        } else {
          glEnable(GL_ALPHA_TEST);
        }
      }
    } else {                    /* triangle-based spheres */

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
          if(ai->visRep[cRepSphere]) {
            float vdw = ai->vdw * sphere_scale;
            int c = ai->color;
            float v0 = v[0];
            float v1 = v[1];
            float v2 = v[2];
            repActive = true;

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
#ifdef _PYMOL_GL_DRAWARRAYS
		{
		  int nverts = nc;
		  ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
		  ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
		  int pl = 0;
		  for(c = 0; c < nc; c++) {
		    float *sp_dot_q = &sp_dot[*(q++)][0];
		    normalVals[pl] = sp_dot_q[0]; normalVals[pl+1] = sp_dot_q[1]; normalVals[pl+2] = sp_dot_q[2];
		    vertexVals[pl++] = v0 + vdw * sp_dot_q[0];
		    vertexVals[pl++] = v1 + vdw * sp_dot_q[1];
		    vertexVals[pl++] = v2 + vdw * sp_dot_q[2];
		  }
		  glEnableClientState(GL_VERTEX_ARRAY);
		  glEnableClientState(GL_NORMAL_ARRAY);
		  glVertexPointer(3, GL_FLOAT, 0, vertexVals);
		  glNormalPointer(GL_FLOAT, 0, normalVals);
		  glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
		  glDisableClientState(GL_VERTEX_ARRAY);
		  glDisableClientState(GL_NORMAL_ARRAY);
		  DEALLOCATE_ARRAY(vertexVals)
		  DEALLOCATE_ARRAY(normalVals)
		}
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

    if(!repActive)              /* didn't draw a single sphere, so we can skip this representation next time around */
      cs->Active[cRepSphere] = false;
  }
}
#ifdef _PYMOL_GL_DRAWARRAYS
void RepSphereRenderTriangleStripsES(int nvertsarg, float *varg){
  int nverts = nvertsarg;
  float *v = varg;
  ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
  ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
  int pl = 0;
  float restart;
  while (nverts>0){
    restart = *(v++);
    if (restart){
      if (restart == 2.0){
	normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
	vertexVals[pl] = v[3]; vertexVals[pl+1] = v[4]; vertexVals[pl+2] = v[5];
	nverts--;
	pl += 3;
      }
      normalVals[pl] = v[0]; normalVals[pl++] = v[1]; normalVals[pl++] = v[2];
      v += 3;
      vertexVals[pl] = v[0]; vertexVals[pl+1] = v[1]; vertexVals[pl+2] = v[2];
      v += 3;
      pl += 3;
      normalVals[pl] = v[0]; normalVals[pl++] = v[1]; normalVals[pl++] = v[2];
      v += 3;
      vertexVals[pl] = v[0]; vertexVals[pl+1] = v[1]; vertexVals[pl+2] = v[2];
      v += 3;
      pl += 3;
      nverts -= 2;
    }
    normalVals[pl] = v[0]; normalVals[pl++] = v[1]; normalVals[pl++] = v[2];
    v += 3;
    vertexVals[pl] = v[0]; vertexVals[pl+1] = v[1]; vertexVals[pl+2] = v[2];
    nverts--;
    v += 3;
    pl += 3;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, normalVals);
  glVertexPointer(3, GL_FLOAT, 0, vertexVals);
  glDrawArrays(GL_TRIANGLE_STRIP, 0, nvertsarg);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  DEALLOCATE_ARRAY(vertexVals)
  DEALLOCATE_ARRAY(normalVals)
}
void RepSphereRenderPointsES(int nvertsarg, float *varg, float *vnarg){
  float *v = varg, *vn = vnarg;
  int nverts = nvertsarg, plc = 0;
  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
  ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
  int pl = 0;

  while (nverts>0){
    colorVals[plc++] = v[0]; colorVals[plc++] = v[1]; colorVals[plc++] = v[2]; colorVals[plc++] = 1.f;
    v += 4;
    if (vn){
      normalVals[pl] = vn[0]; normalVals[pl+1] = vn[1]; normalVals[pl+2] = vn[2];
      vn += 3;
    }
    ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];    
    v += 4;
    nverts--;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
  glColorPointer(4, GL_FLOAT, 0, colorVals);
  glNormalPointer(GL_FLOAT, 0, normalVals);
  glDrawArrays(GL_POINTS, 0, nvertsarg);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  DEALLOCATE_ARRAY(ptsVals)
  DEALLOCATE_ARRAY(colorVals)
  DEALLOCATE_ARRAY(normalVals)
}

void RepSphereRenderMode5PointsES(int nvertsarg, float *varg, float zz_factor, float s_factor, float *dim_add, float _1){
  int nverts = nvertsarg;
  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
    int pl = 0, plc = 0;
  float r,g,b, *v = varg;
  while (nverts>0){
    r = zz_factor * v[0] + s_factor;
    g = zz_factor * v[1] + s_factor;
    b = zz_factor * v[2] + s_factor;
    colorVals[plc++] = r > _1 ? _1 : r; colorVals[plc++] = g > _1 ? _1 : g; colorVals[plc++] = b > _1 ? _1 : b; colorVals[plc++] = 1.f;
    v += 4;
    ptsVals[pl++] = v[0] + dim_add[0]; ptsVals[pl++] = v[1] + dim_add[1]; ptsVals[pl++] = v[2] + dim_add[2];
    v += 4;
    nverts--;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
  glColorPointer(4, GL_FLOAT, 0, colorVals);
  glDrawArrays(GL_POINTS, 0, nvertsarg);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  DEALLOCATE_ARRAY(ptsVals)
  DEALLOCATE_ARRAY(colorVals)
}
void RepSphereRenderPointsDefaultES(RepSphere * I, Picking **pick, int nvertsarg, int iarg, Pickable *parg, float *varg){
  int nverts = nvertsarg;
  float cur_color[3] = { 0.F, 0.F, 0.F };
  int i = iarg;
  Pickable *p = parg;
  float *v = varg;
  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
  ALLOCATE_ARRAY(GLubyte,colorVals,nverts*3)
  int pl = 0, j, plc = 0;
  while(nverts>0){
    int skip = (p[1].index < 0);
    if(!skip) {
      i++;
      if(!(*pick)[0].src.bond) {
	/* pass 1 - low order bits *            */
	cur_color[0] = (uchar) ((i & 0xF) << 4);
	cur_color[1] = (uchar) ((i & 0xF0) | 0x8);
	cur_color[2] = (uchar) ((i & 0xF00) >> 4);
	VLACheck((*pick), Picking, i);
	p++;
	(*pick)[i].src = *p;    /* copy object and atom info */
	(*pick)[i].context = I->R.context;
      } else {
	/* pass 2 - high order bits */
	j = i >> 12;
	cur_color[0] = (uchar) ((j & 0xF) << 4);
	cur_color[1] = (uchar) ((j & 0xF0) | 0x8);
	cur_color[2] = (uchar) ((j & 0xF00) >> 4);
      }
    } else {
      p++;
    }
    v += 4;
    if (!skip){
      colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
      ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];
      nverts--;
    }
    v += 4;
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
  glColorPointer(4, GL_UNSIGNED_BYTE, 0, colorVals);
  glDrawArrays(GL_POINTS, 0, nvertsarg);
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  DEALLOCATE_ARRAY(ptsVals)
  DEALLOCATE_ARRAY(colorVals)
}

#endif


int RenderSphereMode_Direct(PyMOLGlobals *G, RepSphere *I, RenderInfo * info, int carg, float **vptr, float alpha, SphereRec *sphereRecPtr){
  short use_shader, generate_shader_cgo = 0;
  float *v = *vptr;
  int c = carg;
  int ok = true;
  use_shader = (int) SettingGet(G, cSetting_sphere_use_shader) & 
    (int) SettingGet(G, cSetting_use_shaders);
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
	  register float vdw = v[3];
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
	  register float vdw = v[3];
	  glTranslatef(v[0], v[1], v[2]);
	  if((vdw != last_vdw) || (!dlist)) {
	    last_vdw = vdw;
	    q = sp->Sequence;
	    s = sp->StripLen;
#ifdef _PYMOL_GL_CALLLISTS
	    if(!dlist)
	      dlist = glGenLists(1);
	    if(dlist) {
	      glNewList(dlist, GL_COMPILE_AND_EXECUTE);
	    }
#endif
	    
	    for(b = 0; b < sp->NStrip; b++) {
	      int d;
#ifdef _PYMOL_GL_DRAWARRAYS
	      {
		int nverts = *s;
		ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
		ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
		int pl = 0;
		
		for(d = 0; d < (*s); d++) {
		  float *norm = sp_dot[*(q++)];
		  normalVals[pl] = norm[0]; normalVals[pl+1] = norm[1]; normalVals[pl+2] = norm[2];
		  vertexVals[pl++] = vdw * norm[0]; vertexVals[pl++] = vdw * norm[1];
		  vertexVals[pl++] = vdw * norm[2];
		}
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, vertexVals);
		glNormalPointer(GL_FLOAT, 0, normalVals);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		DEALLOCATE_ARRAY(vertexVals)
		DEALLOCATE_ARRAY(normalVals)
	      }
#else
	      glBegin(GL_TRIANGLE_STRIP);
	      for(d = 0; d < (*s); d++) {
		float *norm = sp_dot[*(q++)];
		glNormal3fv(norm);
		glVertex3f(vdw * norm[0], vdw * norm[1], vdw * norm[2]);
	      }
	      glEnd();
#endif
	      s++;
	    }
#ifdef _PYMOL_GL_CALLLISTS
	    if(dlist) {
	      glEndList();
	    }
	  } else {
	    glCallList(dlist);
#endif
	  }
	  glTranslatef(-v[0], -v[1], -v[2]);
	}
	(*vptr)+=4; v = *vptr;
      }
#ifdef _PYMOL_GL_CALLLISTS
      if(dlist)
	glDeleteLists(dlist, 1);
#endif
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

void RenderSphereMode_Sprites(PyMOLGlobals *G, RepSphere *I, RenderInfo *info, int sphere_mode, int carg, float **vptr, float **vnptr){
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
#ifdef _PYMOL_GL_DRAWARRAYS
  {
    float *startv = v, *startvn = vn;
    int nverts = 0;
    while(c--) {
      if(last_radius != (cur_radius = v[7])) {
	size = cur_radius * pixel_scale;
	if (nverts>0){
	  RepSphereRenderPointsES(nverts, startv, startvn);
	  nverts = 0;
	  startv = v;
	  startvn = vn;
	}
	if(clamp_size_flag)
	  if(size > max_size)
	    size = max_size;
	glPointSize(size);
	last_radius = cur_radius;
      }
      if(vn) {
	(*vnptr)+=3; vn = *vnptr;
      }
      nverts++;
      (*vptr)+=8; v = *vptr;
    }
    if (nverts>0){
      RepSphereRenderPointsES(nverts, startv, startvn);
    }
  }
#else
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
#endif
  if(sphere_mode == 3) {
    glDisable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.05F);
  } else {
    glEnable(GL_ALPHA_TEST);
  }
}
 
void RenderSphereMode_Points(PyMOLGlobals *G, RepSphere *I, RenderInfo *info, float radius, int carg){
  register float _1 = 1.0F;
  register float _2 = 2.0F;
  float pixel_scale = 1.0F / info->vertex_scale;
  int repeat = true;
  register float x_add = 0.0F, y_add = 0.0F, z_add = 0.0F;
#ifdef _PYMOL_GL_DRAWARRAYS
  float dim_add[3] = { 0.0F, 0.0F, 0.0F };
#endif
  register float z_factor = 0.0F, r_factor = 1.0F;
  register float largest;
  register float r, g, b;
  register float s_factor = 0.0F;
  register float zz_factor;
  register float clamp_radius;
  register float last_size;
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
#ifdef _PYMOL_GL_DRAWARRAYS
    {
      int nverts;
      float *startv = v;  
      (void)r; (void)g; (void)b;
      nverts = 0;
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
	    if (nverts>0){
	      RepSphereRenderMode5PointsES(nverts, startv, zz_factor, s_factor, dim_add, _1);
	      startv = v;
	      nverts = 0;
	    }
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
	}
	x_add = z_factor * clamp_radius * info->view_normal[0];
	y_add = z_factor * clamp_radius * info->view_normal[1];
	z_add = z_factor * clamp_radius * info->view_normal[2];
	dim_add[0] = x_add; dim_add[1] = y_add; dim_add[2] = z_add;
	last_radius = cur_radius;
	last_size = size;
      }
      nverts++;
      v += 8;
    }
    if (nverts>0){
      RepSphereRenderMode5PointsES(nverts, startv, zz_factor, s_factor, dim_add, _1);
    }
  }
#else
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
#endif
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
}

void RenderSphereMode_9(PyMOLGlobals *G, RepSphere *I, RenderInfo *info, int n_quad_verts, float **vptr, float radius, int carg){
  int c = carg;
  int vc = 0;
  int cc = 0;
  int ac = 0;
  int attr;
  CShaderPrg *shaderPrg;
  short use_shader, generate_shader_cgo = 0;
  float *v = *vptr;
  use_shader = (int) SettingGet(G, cSetting_sphere_use_shader) & 
    (int) SettingGet(G, cSetting_use_shaders);

  if (I->shaderCGO && !use_shader){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  if (use_shader){
    if (!I->shaderCGO){
      I->shaderCGO = CGONew(G);
      I->shaderCGO->use_shader = true;
      generate_shader_cgo = 1;
    } else {
      I->shaderCGO->enable_shaders = true;
      CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
      return;
    }
  }
  
  if (generate_shader_cgo){
    CGOEnable(I->shaderCGO, GL_LIGHTING);
    while (c--) {
      CGOAlpha(I->shaderCGO, v[3]);
      CGOColorv(I->shaderCGO, v);
      CGOSphere(I->shaderCGO, v+4, v[7]);
      (*vptr)+=8; v = *vptr;
    }
    CGOStop(I->shaderCGO);
    {
      CGO *convertcgo = NULL;
      convertcgo = CGOOptimizeSpheresToVBONonIndexed(I->shaderCGO, 0);
      if (convertcgo){
	CGOFree(I->shaderCGO);    
	I->shaderCGO = convertcgo;
      }
    }
    
    {
      I->shaderCGO->enable_shaders = true;
      CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
      return;
    }
  } else {
    ALLOCATE_ARRAY(GLfloat,colorVals,c*4*4)
    ALLOCATE_ARRAY(GLfloat,vertexVals,c*4*3)
    ALLOCATE_ARRAY(GLfloat,attribVals,c*4*3)
    n_quad_verts = c * 4;
    if(Feedback(G, FB_OpenGL, FB_Debugging)) {
      PRINTF "GLSL Sphere Shader: n_quad_verts: %d\n", 
	n_quad_verts ENDF(G);
    }
    shaderPrg = CShaderPrg_Enable_SphereShader(G, "spheredirect");
    
    attr = CShaderPrg_GetAttribLocation(shaderPrg, "sphere_attributes");
    
    while (c--) {
      radius = v[7];
      attribVals[ac++] = -1.0;
      attribVals[ac++] = -1.0;
      attribVals[ac++] = radius;
      colorVals[cc++] = v[0];
      colorVals[cc++] = v[1];
      colorVals[cc++] = v[2];
      colorVals[cc++] = v[3];
      vertexVals[vc++] = v[4];
      vertexVals[vc++] = v[5];
      vertexVals[vc++] = v[6];
      
      attribVals[ac++] =  1.0;
      attribVals[ac++] = -1.0;
      attribVals[ac++] = radius;
      colorVals[cc++] = v[0];
      colorVals[cc++] = v[1];
      colorVals[cc++] = v[2];
      colorVals[cc++] = v[3];
      vertexVals[vc++] = v[4];
      vertexVals[vc++] = v[5];
      vertexVals[vc++] = v[6];
      
      attribVals[ac++] = 1.0;
      attribVals[ac++] = 1.0;
      attribVals[ac++] = radius;
      colorVals[cc++] = v[0];
      colorVals[cc++] = v[1];
      colorVals[cc++] = v[2];
      colorVals[cc++] = v[3];
      vertexVals[vc++] = v[4];
      vertexVals[vc++] = v[5];
      vertexVals[vc++] = v[6];
      
      attribVals[ac++] = -1.0;
      attribVals[ac++] =  1.0;
      attribVals[ac++] = radius;
      colorVals[cc++] = v[0];
      colorVals[cc++] = v[1];
      colorVals[cc++] = v[2];
      colorVals[cc++] = v[3];
      vertexVals[vc++] = v[4];
      vertexVals[vc++] = v[5];
      vertexVals[vc++] = v[6];
      
      glBegin(GL_QUADS);
      glColor4f(v[0], v[1], v[2], v[3]);
      glVertexAttrib3f(attr, -1.0, -1.0, radius);
      glVertex3f(v[4], v[5], v[6]);
      glVertexAttrib3f(attr,  1.0, -1.0, radius);
      glVertex3f(v[4], v[5], v[6]);
      glVertexAttrib3f(attr,  1.0,  1.0, radius);
      glVertex3f(v[4], v[5], v[6]);
      glVertexAttrib3f(attr, -1.0,  1.0, radius);
      glVertex3f(v[4], v[5], v[6]);
      glEnd();
      
      (*vptr)+=8; v = *vptr;
    }
    CShaderPrg_Disable(shaderPrg);
    
    DEALLOCATE_ARRAY(colorVals)
    DEALLOCATE_ARRAY(vertexVals)
    DEALLOCATE_ARRAY(attribVals)
  }
}

void RenderSphereMode_ARB(PyMOLGlobals *G, RenderInfo *info, float **vptr, int carg){
  if(sphereARBShaderPrg) {
    int c = carg;
    float fog_info[3];
    float _1 = 1.0F;
    float _00[2] = { 0.0F, 0.0F };
    float _01[2] = { 0.0F, 1.0F };
    float _11[2] = { 1.0F, 1.0F };
    float _10[2] = { 1.0F, 0.0F };
    register float v0, v1, v2, nv0, nv1, nv2, nv3, v3;
    register float *m = info->pmv_matrix;
    register float cutoff = 1.2F;
    register float m_cutoff = -cutoff;
    register float z_front, z_back;
    float *v = *vptr;
    float last_radius, cur_radius;
    /* compute -Ze = (Wc) of fog start */
    nv3 =
      (info->front +
       (info->back - info->front) * SettingGetGlobal_f(G,
						       cSetting_fog_start));
    /* compute Zc of fog start using std. perspective transformation */
    nv2 =
      (nv3 * (info->back + info->front) -
       2 * (info->back * info->front)) / (info->back - info->front);
    /* compute Zc/Wc to get normalized depth coordinate of fog start */
    nv0 = (nv2 / nv3);
    fog_info[0] = (nv0 * 0.5) + 0.5;
    /* printf("%8.3f %8.3f %8.3f %8.3f\n", nv3, nv2, nv0, fog_info[0]); */
    fog_info[1] = 1.0F / (1.0 - fog_info[0]);     /* effective range of fog */
    
    z_front = info->stereo_front;
    z_back = info->back + ((info->back + info->front) * 0.25);
    
    if(Feedback(G, FB_OpenGL, FB_Debugging))
      PyMOLCheckOpenGLErr("before shader");
    
    CShaderPrg_Enable_SphereShaderARB(G);
    
    {
      glNormal3fv(info->view_normal);
      
      (*vptr)+=4; v = *vptr;
      last_radius = -1.f;
      glBegin(GL_QUADS);
      while(c--) {
	v3 = v[3];
	v0 = v[0];
	v1 = v[1];
	v2 = v[2];
	
	if(last_radius != (cur_radius = v3)) {
	  glEnd();
	  glProgramEnvParameter4fARB(GL_VERTEX_PROGRAM_ARB,
				     0, 0.0F, 0.0F, v3, 0.0F);
	  glProgramEnvParameter4fARB(GL_FRAGMENT_PROGRAM_ARB,
				     0, fog_info[0], fog_info[1], 0.0F,
				     0.0F);
	  glBegin(GL_QUADS);
	  last_radius = cur_radius;
	}
	
	/*  MatrixTransformC44f4f(info->pmv_matrix, v, nv); */
	
	nv3 = m[3] * v0 + m[7] * v1 + m[11] * v2 + m[15]; /* compute Wc */
	
	if(((nv3 - cur_radius) > z_front) && (nv3 < z_back)) {    /* is it within the viewing volume? */
	  nv0 = m[0] * v0 + m[4] * v1 + m[8] * v2 + m[12];
	  
	  nv3 = _1 / nv3;
	  nv1 = m[1] * v0 + m[5] * v1 + m[9] * v2 + m[13];
	  nv0 *= nv3;
	  nv1 *= nv3;
	  
	  if((nv0 < cutoff) && (nv0 > m_cutoff) &&
	     (nv1 < cutoff) && (nv1 > m_cutoff)) {
	    
	    glColor3fv(v - 4);
	    
	    glTexCoord2fv(_00);
	    glVertex3fv(v);
	    
	    glTexCoord2fv(_10);
	    glVertex3fv(v);
	    
	    glTexCoord2fv(_11);
	    glVertex3fv(v);
	    
	    glTexCoord2fv(_01);
	    glVertex3fv(v);
	  }
	}
	(*vptr)+=8; v = *vptr;
      }
      glEnd();
      
      CShaderPrg_DisableARB(sphereARBShaderPrg);
      if(Feedback(G, FB_OpenGL, FB_Debugging))
	PyMOLCheckOpenGLErr("after shader");
    }
  }
}

/* simple, default point width points -- modes 1 or 6 */
void RenderSphereMode_1_or_6(PyMOLGlobals *G, RepSphere *I, RenderInfo *info, float **vptr, float **vnptr, int carg, float alpha){
  int c = carg;
  float *v = *vptr, *vn = *vnptr;
  glPointSize(SettingGet_f
	      (G, I->R.cs->Setting, I->R.obj->Setting,
	       cSetting_sphere_point_size));
  glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_ALPHA_TEST);
  
#ifdef _PYMOL_GL_DRAWARRAYS
  {
    int nverts = c;
    ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
    ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
    ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
    int pl =0, plc = 0;
    while (c--){
      colorVals[plc++] = v[0]; colorVals[plc++] = v[1]; 
      colorVals[plc++] = v[2]; colorVals[plc++] = alpha;
      (*vptr)+=4; v = *vptr;
      if (vn){
	normalVals[pl] = vn[0]; normalVals[pl+1] = vn[1]; normalVals[pl+2] = vn[2];
	(*vnptr)+=3; vn = *vnptr;
      }
      vertexVals[pl++] = v[0]; vertexVals[pl++] = v[1]; 
      vertexVals[pl++] = v[2]; 
      (*vptr)+=4; v = *vptr;
    }
    glEnable(GL_LIGHTING);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    if (vn){
      glEnableClientState(GL_NORMAL_ARRAY);
      glNormalPointer(GL_FLOAT, 0, normalVals);
    }
    glVertexPointer(3, GL_FLOAT, 0, vertexVals);
    glColorPointer(4, GL_FLOAT, 0, colorVals);
    glDrawArrays(GL_POINTS, 0, nverts);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    if (vn){
      glDisableClientState(GL_NORMAL_ARRAY);
    }
    DEALLOCATE_ARRAY(colorVals)
    DEALLOCATE_ARRAY(normalVals)
    DEALLOCATE_ARRAY(vertexVals)
  }
#else
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
#endif
  glEnable(GL_ALPHA_TEST);
}

void RenderSphereMode_Default(PyMOLGlobals *G, RepSphere *I, int carg, float **vptr, float alpha, SphereRec *sp){
  int cc, c = carg, flag, a;
  int variable_alpha = I->VariableAlphaFlag;
  int use_dlst;
  int *nt;
  float *v = *vptr;
  float restart;
  use_dlst = (int) SettingGet(G, cSetting_use_display_lists);
  
#ifdef _PYMOL_GL_CALLLISTS
  if(use_dlst && I->R.displayList) {
    glCallList(I->R.displayList);
  } else {                /* display list */
    if(use_dlst) {
      if(!I->R.displayList) {
	I->R.displayList = glGenLists(1);
	if(I->R.displayList) {
	  glNewList(I->R.displayList, GL_COMPILE_AND_EXECUTE);
	}
      }
    }
#endif
    if(I->cullFlag) {
      if((alpha == 1.0) && (!variable_alpha)) {
	nt = I->NT;       /* number of passes for each sphere */
	while(c--) {      /* iterate through all atoms */
	  glColor3fv(v);
	  (*vptr)+=4; v = *vptr;
	  cc = *(nt++);
	  flag = 0;
#ifdef _PYMOL_GL_DRAWARRAYS
	  {
	    int nverts=0;
	    float *vstart = v;
	    while(cc--) {   /* execute loop this many times */
	      restart = *(v++);
	      if(restart) {
		if(flag) {
		  RepSphereRenderTriangleStripsES(nverts, vstart);
		  nverts = 0;
		  vstart = v - 1;
		}
		if(restart == 2.0) {        /* swap triangle polarity */
		  nverts++;
		}
		nverts += 2;
		(*vptr)+=12; v = *vptr;
	      }
	      nverts++;
	      (*vptr)+=6; v = *vptr;
	      flag = 1;
	    }
	    if (nverts>0){
	      RepSphereRenderTriangleStripsES(nverts, vstart);
	    }
	  }
#else
	  glBegin(GL_TRIANGLE_STRIP);
	  while(cc--) {   /* execute loop this many times */
	    restart = *v;
	    (*vptr)++; v = *vptr;
	    if(restart) {
	      if(flag) {
		glEnd();
		glBegin(GL_TRIANGLE_STRIP);
	      }
	      if(restart == 2.0) {        /* swap triangle polarity */
		glNormal3fv(v);
		glVertex3fv(v + 3);
	      }
	      glNormal3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glVertex3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glNormal3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glVertex3fv(v);
	      (*vptr)+=3; v = *vptr;
	    }
	    glNormal3fv(v);
	    (*vptr)+=3; v = *vptr;
	    glVertex3fv(v);
	    (*vptr)+=3; v = *vptr;
	    flag = 1;
	  }
	  glEnd();
#endif
	}
      } else {
	nt = I->NT;       /* number of passes for each sphere */
	while(c--) {      /* iterate through all atoms */
	  glColor4f(v[0], v[1], v[2], v[3]);
	  (*vptr)+=4; v = *vptr;
	  cc = *(nt++);
	  flag = 0;
#ifdef _PYMOL_GL_DRAWARRAYS
	  {
	    int nverts=0;
	    float *vstart = v;
	    while(cc--) {   /* execute loop this many times */
	      restart = *(v++);
	      if(restart) {
		if(flag) {
		  RepSphereRenderTriangleStripsES(nverts, vstart);
		  nverts = 0;
		  vstart = v - 1;
		}
		if(restart == 2.0) {        /* swap triangle polarity */
		  nverts++;
		}
		nverts += 2;
		(*vptr)+=12; v = *vptr;
	      }
	      nverts++;
	      (*vptr)+=6; v = *vptr;
	      flag = 1;
	    }
	    if (nverts>0){
	      RepSphereRenderTriangleStripsES(nverts, vstart);
	    }
	  }
#else
	  glBegin(GL_TRIANGLE_STRIP);
	  while(cc--) {   /* execute loop this many times */
	    restart = *v;
	    (*vptr)++; v = *vptr;
	    if(restart) {
	      if(flag) {
		glEnd();
		glBegin(GL_TRIANGLE_STRIP);
	      }
	      if(restart == 2.0) {        /* swap triangle polarity */
		glNormal3fv(v);
		glVertex3fv(v + 3);
	      }
	      glNormal3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glVertex3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glNormal3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glVertex3fv(v);
	      (*vptr)+=3; v = *vptr;
	    }
	    glNormal3fv(v);
	    (*vptr)+=3; v = *vptr;
	    glVertex3fv(v);
	    (*vptr)+=3; v = *vptr;
	    flag = 1;
	  }
	  glEnd();
#endif
	}
      }
    } else if(sp) {
      if((alpha == 1.0) && !variable_alpha) {
	while(c--) {
	  glColor3fv(v);
	  (*vptr)+=4; v = *vptr;
	  for(a = 0; a < sp->NStrip; a++) {
#ifdef _PYMOL_GL_DRAWARRAYS
	    cc = sp->StripLen[a];
	    {
	      int nverts = cc;
	      ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
	      ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
	      int pl = 0;
	      
	      while(cc--) {
		normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
		(*vptr)+=3; v = *vptr;
		vertexVals[pl] = v[0]; vertexVals[pl+1] = v[1]; vertexVals[pl+2] = v[2];
		(*vptr)+=3; v = *vptr;
		pl += 3;
	      }
	      glEnableClientState(GL_VERTEX_ARRAY);
	      glEnableClientState(GL_NORMAL_ARRAY);
	      glVertexPointer(3, GL_FLOAT, 0, vertexVals);
	      glNormalPointer(GL_FLOAT, 0, normalVals);
	      glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
	      glDisableClientState(GL_VERTEX_ARRAY);
	      glDisableClientState(GL_NORMAL_ARRAY);
	      DEALLOCATE_ARRAY(normalVals)
	      DEALLOCATE_ARRAY(vertexVals)
	    }
#else
	    glBegin(GL_TRIANGLE_STRIP);
	    cc = sp->StripLen[a];
	    while(cc--) {
	      glNormal3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glVertex3fv(v);
	      (*vptr)+=3; v = *vptr;
	    }
	    glEnd();
#endif
	  }
	}
      } else {
	while(c--) {
	  glColor4f(v[0], v[1], v[2], v[3]);
	  (*vptr)+=4; v = *vptr;
	  for(a = 0; a < sp->NStrip; a++) {
#ifdef _PYMOL_GL_DRAWARRAYS
	    {
	      int nverts = cc;
	      ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
	      ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
	      int pl = 0;
	      
	      while(cc--) {
		normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
		(*vptr)+=3; v = *vptr;
		vertexVals[pl] = v[0]; vertexVals[pl+1] = v[1]; vertexVals[pl+2] = v[2];
		(*vptr)+=3; v = *vptr;
		pl += 3;
	      }
	      glEnableClientState(GL_VERTEX_ARRAY);
	      glEnableClientState(GL_NORMAL_ARRAY);
	      glVertexPointer(3, GL_FLOAT, 0, vertexVals);
	      glNormalPointer(GL_FLOAT, 0, normalVals);
	      glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
	      glDisableClientState(GL_VERTEX_ARRAY);
	      glDisableClientState(GL_NORMAL_ARRAY);
	      DEALLOCATE_ARRAY(normalVals)
	      DEALLOCATE_ARRAY(vertexVals)
	    }
#else
	    glBegin(GL_TRIANGLE_STRIP);
	    cc = sp->StripLen[a];
	    while(cc--) {
	      glNormal3fv(v);
	      (*vptr)+=3; v = *vptr;
	      glVertex3fv(v);
	      (*vptr)+=3; v = *vptr;
	    }
	    glEnd();
#endif
	  }
	}
      }
    }
#ifdef _PYMOL_GL_CALLLISTS
    if(use_dlst && I->R.displayList) {
      glEndList();
    }
  }
#endif
}

static void RepSphereRender(RepSphere * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V, *vc, *vn = I->VN;
  int c = I->N;
  int cc = 0;
  int a;
  SphereRec *sp = I->SP;
  float alpha;
  int n_quad_verts;
  float radius;
  int ok = true;

#ifdef _PYMOL_OPENGL_SHADERS
  /* TO DO -- garbage collect -- IMPORTANT! */
  {
    int sphere_mode = SettingGet_i(G, I->R.cs->Setting,
                                   I->R.obj->Setting,
                                   cSetting_sphere_mode);
    if (!ray && sphere_mode == 5 && G->HaveGUI && G->ValidContext && !sphereARBShaderPrg){
      sphereARBShaderPrg = CShaderPrg_NewARB(G, "sphere_arb", sphere_arb_vs, sphere_arb_fs);      
    }
  }
#endif
  alpha =
    SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_sphere_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  if(ray) {
    ray->fTransparentf(ray, 1.0 - alpha);
    if(I->spheroidFlag) {
      if(sp) {
        while(c--) {
          vc = v;
          v += 3;
          for(a = 0; ok && a < sp->NStrip; a++) {
            cc = sp->StripLen[a];
            while(ok && (cc--) > 2) {
              ok &= ray->fTriangle3fv(ray, v + 3, v + 9, v + 15, v, v + 6, v + 12, vc, vc, vc);
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
          ray->fTransparentf(ray, 1.0F - v[3]);
        }
        ray->fColor3fv(ray, v);
        v += 4;

	ok &= ray->fSphere3fv(ray, v, *(v + 3));
        v += 4;
      }
    }
    ray->fTransparentf(ray, 0.0);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      int trans_pick_mode = SettingGet_i(G, I->R.cs->Setting,
                                         I->R.obj->Setting,
                                         cSetting_transparency_picking_mode);

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
              i++;
              if(!(*pick)[0].src.bond) {
                /* pass 1 - low order bits *            */
                glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8),
                           (uchar) ((i & 0xF00) >> 4));
                VLACheck((*pick), Picking, i);
                p++;
                (*pick)[i].src = *p;    /* copy object and atom info */
                (*pick)[i].context = I->R.context;
              } else {
                /* pass 2 - high order bits */
                j = i >> 12;
                glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                           (uchar) ((j & 0xF00) >> 4));
              }
            } else {
              p++;
            }

            v += 4;
            for(a = 0; a < sp->NStrip; a++) {
              cc = sp->StripLen[a];
              if(!skip) {
#ifdef _PYMOL_GL_DRAWARRAYS
		{
		  int nverts = cc-1;
		  ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
		  ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
		  int pl = 0;
		  while((cc--) > 0) {
		    normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
		    vertexVals[pl++] = v[3]; vertexVals[pl++] = v[4]; vertexVals[pl++] = v[5];
		    v += 6;
		  }
		  glEnableClientState(GL_VERTEX_ARRAY);
		  glEnableClientState(GL_NORMAL_ARRAY);
		  glVertexPointer(3, GL_FLOAT, 0, vertexVals);
		  glNormalPointer(GL_FLOAT, 0, normalVals);
		  glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
		  glDisableClientState(GL_VERTEX_ARRAY);
		  glDisableClientState(GL_NORMAL_ARRAY);
		  DEALLOCATE_ARRAY(vertexVals)
		  DEALLOCATE_ARRAY(normalVals)
		}
#else
                glBegin(GL_TRIANGLE_STRIP);
                while((cc--) > 0) {
                  glNormal3fv(v);
                  glVertex3fv(v + 3);
                  v += 6;
                }
                glEnd();
#endif
              } else {
                while((cc--) > 0) {
                  v += 6;
                }
              }
            }
          }
        } else {
          register float last_radius = -1.0F;
          register float cur_radius;
          register float pixel_scale = 1.0F / info->vertex_scale;
          register float max_size = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting,
                                                 cSetting_sphere_point_max_size) * 3;
          register int clamp_size_flag = (max_size >= 0.0F);
          register float size;
          int sphere_mode = SettingGet_i(G, I->R.cs->Setting,
                                         I->R.obj->Setting,
                                         cSetting_sphere_mode);
#ifdef _PYMOL_GL_DRAWARRAYS
	  int starti, nverts = 0;
	  float *startv;
	  Pickable *startp;
	  
#endif

          if(!sp) {
            switch (sphere_mode) {
            case -1:
            case 0:
              break;
            case 1:
            case 6:
#ifndef _PYMOL_GL_DRAWARRAYS
              glBegin(GL_POINTS);
#endif
              break;
            case 5:
            case 4:
            case 3:
            case 8:
              glEnable(GL_POINT_SMOOTH);
              glAlphaFunc(GL_GREATER, 0.5F);
              glEnable(GL_ALPHA_TEST);
              glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
              glPointSize(1.0F);
              pixel_scale *= 2.0F;
#ifndef _PYMOL_GL_DRAWARRAYS
              glBegin(GL_POINTS);
#endif
              break;
            case 2:
            case 7:
              glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
              glDisable(GL_POINT_SMOOTH);
              glDisable(GL_ALPHA_TEST);
              pixel_scale *= 1.4F;
#ifndef _PYMOL_GL_DRAWARRAYS
              glBegin(GL_POINTS);
#endif
              break;
            default:
              glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
              glDisable(GL_POINT_SMOOTH);
              glDisable(GL_ALPHA_TEST);
              glPointSize(SettingGet_f
                          (G, I->R.cs->Setting, I->R.obj->Setting,
                           cSetting_sphere_point_size));
#ifndef _PYMOL_GL_DRAWARRAYS
              glBegin(GL_POINTS);
#endif
              break;
            }
          }

          v = I->VC;
          c = I->NC;
#ifdef _PYMOL_GL_DRAWARRAYS
	  starti = i;
	  startp = p;
	  startv = v;
#endif
          while(c--) {
            int skip = (p[1].index < 0);
            if(!skip) {
              i++;
              if(!(*pick)[0].src.bond) {
                /* pass 1 - low order bits *            */
                glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8),
                           (uchar) ((i & 0xF00) >> 4));
                VLACheck((*pick), Picking, i);
                p++;
#ifndef _PYMOL_GL_DRAWARRAYS
                (*pick)[i].src = *p;    /* copy object and atom info */
                (*pick)[i].context = I->R.context;
#endif
              } else {
                /* pass 2 - high order bits */
                j = i >> 12;
                glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                           (uchar) ((j & 0xF00) >> 4));
              }
            } else {
              p++;
            }

            if(sp) {
              if(!skip) {
                int *s, *q, b;
                float *v0, vdw;

                v0 = v + 4;
                vdw = v[7];
                q = sp->Sequence;
                s = sp->StripLen;
                for(b = 0; b < sp->NStrip; b++) {
#ifdef _PYMOL_GL_DRAWARRAYS
		  {
		    int nverts = *s;
		    ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
		    ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
		    int pl = 0;
		    
		    for(cc = 0; cc < (*s); cc++) {
		      normalVals[pl] = sp->dot[*q][0]; normalVals[pl+1] = sp->dot[*q][1]; normalVals[pl+2] = sp->dot[*q][2];
		      vertexVals[pl++] = v0[0] + vdw * sp->dot[*q][0]; vertexVals[pl++] = v0[1] + vdw * sp->dot[*q][1]; 
		      vertexVals[pl++] = v0[2] + vdw * sp->dot[*q][2];
		      q++;
		    }
		    glEnableClientState(GL_VERTEX_ARRAY);
		    glEnableClientState(GL_NORMAL_ARRAY);
		    glVertexPointer(3, GL_FLOAT, 0, vertexVals);
		    glNormalPointer(GL_FLOAT, 0, normalVals);
		    glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
		    glDisableClientState(GL_VERTEX_ARRAY);
		    glDisableClientState(GL_NORMAL_ARRAY);
		    DEALLOCATE_ARRAY(vertexVals)
		    DEALLOCATE_ARRAY(normalVals)
		  }
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
              v += 8;
            } else {
              switch (sphere_mode) {
              case -1:
              case 0:          /* memory-efficient sphere rendering */
                if(I->SSP) {
                  SphereRec *sp = I->SSP;
                  Vector3f *sp_dot = sp->dot;
                  int b, *q, *s;
                  v += 4;
                  if(!skip) {
                    register float vdw = v[3];
                    glTranslatef(v[0], v[1], v[2]);
                    q = sp->Sequence;
                    s = sp->StripLen;
                    for(b = 0; b < sp->NStrip; b++) {
                      int d;
#ifdef _PYMOL_GL_DRAWARRAYS
		      {
			int nverts = *s;
			ALLOCATE_ARRAY(GLfloat,vertexVals,nverts*3)
			ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
			int pl = 0;
			
			for(d = 0; d < (*s); d++) {
			  float *norm = sp_dot[*(q++)];
			  normalVals[pl] = norm[0]; normalVals[pl+1] = norm[1]; normalVals[pl+2] = norm[2];
			  vertexVals[pl++] = vdw * norm[0]; vertexVals[pl++] = vdw * norm[1];
			  vertexVals[pl++] = vdw * norm[2];
			}
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glVertexPointer(3, GL_FLOAT, 0, vertexVals);
			glNormalPointer(GL_FLOAT, 0, normalVals);
			glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			DEALLOCATE_ARRAY(vertexVals)
			DEALLOCATE_ARRAY(normalVals)
		      }
#else
                      glBegin(GL_TRIANGLE_STRIP);
                      for(d = 0; d < (*s); d++) {
                        float *norm = sp_dot[*(q++)];
                        glNormal3fv(norm);
                        glVertex3f(vdw * norm[0], vdw * norm[1], vdw * norm[2]);
                      }
                      glEnd();
#endif
                      s++;
                    }
                    glTranslatef(-v[0], -v[1], -v[2]);
                  }
                  v += 4;
                }
                break;
              case 2:
              case 3:
              case 4:
              case 5:               
              case 7:
              case 8:
              case 9:
                if(!skip) {
                  if(last_radius != (cur_radius = v[7])) {
                    size = cur_radius * pixel_scale;
#ifdef _PYMOL_GL_DRAWARRAYS
		    /* TODO Need to draw points */
		    if (nverts>0){
		      RepSphereRenderPointsDefaultES(I, pick, nverts, starti, startp, startv);
		      starti = i;
		      startp = p;
		      startv = v;
		      nverts = 0;
		    }
#else
                    glEnd();
#endif
                    if(clamp_size_flag)
                      if(size > max_size)
                        size = max_size;
                    glPointSize(size);
#ifndef _PYMOL_GL_DRAWARRAYS
		    glBegin(GL_POINTS);
#endif
                    last_radius = cur_radius;
                  }
                }
                v += 4;
#ifdef _PYMOL_GL_DRAWARRAYS
                if(!skip)
		  nverts++;
#else
                if(!skip)
                  glVertex3fv(v);
#endif
                v += 4;
                break;
              default:         /* simple, default point width points */
                v += 4;
#ifdef _PYMOL_GL_DRAWARRAYS
                if(!skip)
                  nverts++;
#else
                if(!skip)
                  glVertex3fv(v);
#endif
                v += 4;
                break;
              }
            }
          }
          if(!sp) {
            switch (sphere_mode) {
            case -1:
            case 0:
              break;
            case 3:
            case 4:
            case 8:
              glDisable(GL_POINT_SMOOTH);
              glAlphaFunc(GL_GREATER, 0.05F);

#ifdef _PYMOL_GL_DRAWARRAYS
	      /* TODO Need to draw points */
	      if (nverts>0){
		RepSphereRenderPointsDefaultES(I, pick, nverts, starti, startp, startv);
	      }
#else
              glEnd();
#endif
              break;
            default:
#ifdef _PYMOL_GL_DRAWARRAYS
	      /* TODO Need to draw points */
	      if (nverts>0){
		RepSphereRenderPointsDefaultES(I, pick, nverts, starti, startp, startv);
	      }
#else
              glEnd();
#endif
              glEnable(GL_ALPHA_TEST);
              break;
            }
          }
        }
        (*pick)[0].src.index = i;
      }
    } else {                    /* not pick */

      if(!sp) {
        /* no sp -- we're rendering as points */
        int use_dlst;
        int sphere_mode = SettingGet_i(G, I->R.cs->Setting,
                                       I->R.obj->Setting,
                                       cSetting_sphere_mode);
        v = I->VC;
        c = I->NC;

#ifdef _PYMOL_GL_CALLLISTS
	switch (sphere_mode){
	case 2: case 3: case 7: case 8:
	  /* scaleable reps... */
          if(I->R.displayList) {
            if(I->LastVertexScale != info->vertex_scale) {
              glDeleteLists(I->R.displayList, 1);
              I->R.displayList = 0;
            }
          }
        }
#endif
        I->LastVertexScale = info->vertex_scale;

        use_dlst = (int) SettingGet(G, cSetting_use_display_lists);
        switch (sphere_mode) {
        case -1: case 0: case 4: case 5: case 9:
          use_dlst = 0;
          break;
        }
#ifdef _PYMOL_GL_CALLLISTS
        if(use_dlst && I->R.displayList) {
          glCallList(I->R.displayList);
        } else {                /* display list */

          if(use_dlst) {
            if(!I->R.displayList) {
              I->R.displayList = glGenLists(1);
              if(I->R.displayList) {
                glNewList(I->R.displayList, GL_COMPILE_AND_EXECUTE);
              }
            }
          }
#endif
          if((sphere_mode > 0) && (!info->line_lighting))
            glDisable(GL_LIGHTING);
	  if((sphere_mode == 5)
#ifdef _PYMOL_OPENGL_SHADERS
	     && (!sphereARBShaderPrg)
#endif
	     )
	    sphere_mode = 4;

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
	    RenderSphereMode_Points(G, I, info, radius, c);
	    break;
	  case 5:          /* use vertex/fragment program */
	    RenderSphereMode_ARB(G, info, &v, c);
	    break;
	  case 9: // use GLSL shader
	    RenderSphereMode_9(G, I, info, n_quad_verts, &v, radius, c);
	    break;
	  default:
	    RenderSphereMode_1_or_6(G, I, info, &v, &vn, c, alpha);
            break;
          }
          glEnable(GL_LIGHTING);

#ifdef _PYMOL_GL_CALLLISTS
          if(use_dlst && I->R.displayList) {
            glEndList();
          }
        }
#endif
      } else {                  /* real spheres, drawn with triangles -- not points or impostors */
	RenderSphereMode_Default(G, I, c, &v, alpha, sp);
      }
    }
  }
}

int RepSphereSameVis(RepSphere * I, CoordSet * cs)
{
  int same = true;
  int *lv, *lc, *cc;
  int a;
  AtomInfoType *ai;
  if(I->LastVisib && I->LastColor) {
    ai = cs->Obj->AtomInfo;
    lv = I->LastVisib;
    lc = I->LastColor;
    cc = cs->Color;

    for(a = 0; a < cs->NIndex; a++) {
      if(*(lv++) != (ai + cs->IdxToAtm[a])->visRep[cRepSphere]) {
        same = false;
        break;
      }
      if(*(lc++) != *(cc++)) {
        same = false;
        break;
      }
    }
  } else {
    same = false;
  }
  return (same);
}

static int RadiusOrder(float *list, int a, int b)
{
  return (list[a * 8 + 7] <= list[b * 8 + 7]);
}

Rep *RepSphereNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int ok = true;
  int a, b, c, a1, c1, a2, i, j, k, h, l;
  float *v, *v0, *vc, vdw, v1[3];
  float restart;
  int *q, *s, q0, q1, q2;
  int *lv, *lc, *cc;
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds, *nt, flag;
  int *visFlag = NULL;
  MapType *map = NULL;
  int vFlag;
  AtomInfoType *ai2;
  int spheroidFlag = false;
  float spheroid_scale;
  float *sphLen, sphTmp, *sphNorm, *sphTmpN;
  float sphere_scale, sphere_add = 0.0;
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

  OOCalloc(G, RepSphere);
  CHECKOK(ok, I);
  if (!ok)
    return NULL;
  obj = cs->Obj;
  vFlag = false;
  if(obj->RepVisCache[cRepSphere])
    for(a = 0; a < cs->NIndex; a++) {
      if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepSphere]) {
        vFlag = true;
        break;
      }
    }
  if(!vFlag) {
    OOFreeP(I);
    return (NULL);              /* skip if no dots are visible */
  }
  marked = Calloc(int, obj->NAtom);
  CHECKOK(ok, marked);
  if (ok)
    RepInit(G, &I->R);
  I->shaderCGO = NULL;

  if (ok){
    ds = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_quality);
    sphere_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_mode);
    if(sphere_mode > 0)
      ds = -1;
    
    if(ds < 0) {
      sp = NULL;
    } else {
      if(draw_quality && (ds < 3))
	ds = 3;
      if(ds > 4)
	ds = 4;
      sp = G->Sphere->Sphere[ds];
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
    I->VC = (float *) mmalloc(sizeof(float) * cs->NIndex * 8);
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
    
    if(SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_pickable)) {
      I->R.P = Alloc(Pickable, cs->NIndex + 1);
      CHECKOK(ok, I->R.P);
    }
  }

  I->spheroidFlag = spheroidFlag;
  for(a = 0; ok && a < cs->NIndex; a++) {
    a1 = cs->IdxToAtm[a];
    ati1 = obj->AtomInfo + a1;
    vis_flag = ati1->visRep[cRepSphere];

    if(vis_flag &&
       (!ati1->hetatm) &&
       (!(ati1->flags & cAtomFlag_solvent)) &&
       ((cartoon_side_chain_helper && ati1->visRep[cRepCartoon]) ||
        (ribbon_side_chain_helper && ati1->visRep[cRepRibbon]))) {

      register char *name1 = ati1->name;
      register int prot1 = ati1->protons;

      if(prot1 == cAN_N) {
        if((!name1[1]) && (name1[0] == 'N')) {  /* N */
          register char *resn1 = ati1->resn;
          if(!((resn1[0] == 'P') && (resn1[1] == 'R') && (resn1[2] == 'O')))
            vis_flag = false;
        }
      } else if(prot1 == cAN_O) {
        if((!name1[1]) && (name1[0] == 'O'))
          vis_flag = false;
      } else if(prot1 == cAN_C) {
        if((!name1[1]) && (name1[0] == 'C'))
          vis_flag = false;
      }
    }

    marked[a1] = vis_flag;      /* store temporary visibility information */

    if(vis_flag) {
      float at_sphere_scale;
      int at_sphere_color;
      float at_transp;

      AtomInfoGetSetting_f(G, ati1, cSetting_sphere_scale, sphere_scale,
                           &at_sphere_scale);
      if(AtomInfoGetSetting_f(G, ati1, cSetting_sphere_transparency, transp, &at_transp))
        variable_alpha = true;
      AtomInfoGetSetting_color(G, ati1, cSetting_sphere_color, sphere_color,
                               &at_sphere_color);

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
        c1 = *(cs->Color + a);
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
      
      if(I->cullFlag < 0) {
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
    
    if((I->cullFlag < 2) && (SettingGet(G, cSetting_roving_spheres) != 0.0F))
      I->cullFlag = false;
    
    if(sp && (I->cullFlag < 2) && (!spheroidFlag)) {
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
        UtilCopyMem(pk_tmp, I->R.P, sizeof(Pickable) * (I->NP + 1));

        UtilSortIndex(I->NC, I->VC, ix, (UtilOrderFn *) RadiusOrder);

        UtilCopyMem(I->R.P, pk_tmp, sizeof(Pickable));
        for(a = 0; a < I->NC; a++) {
          UtilCopyMem(I->VC + (a * 8), vc_tmp + (8 * ix[a]), sizeof(float) * 8);
          UtilCopyMem(I->R.P + (a + 1), pk_tmp + ix[a] + 1, sizeof(Pickable));
        }
      }
      FreeP(vc_tmp);
      FreeP(ix);
      FreeP(pk_tmp);
    }

    if(ok && (sphere_mode >= 6) && (sphere_mode < 9) && I->NC) {
      /* compute sphere normals to approximate a surface */
      register float range = 6.0F;
      float *vc = I->VC;
      float *dot = G->Sphere->Sphere[1]->dot[0];
      int n_dot = G->Sphere->Sphere[1]->nDot;
      int nc = I->NC;
      int *active = NULL;
      float *v_tmp = NULL;
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
            register float *vv;
            register int nbr_flag;
            register int n_dot_active, *da;
            register float cut_mult = -1.0F;
            register float range2 = range * range;

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
                      register float cutoff = dst * cut_mult;
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
    }

    I->cullFlag = false;
    I->V = NULL;
    I->NT = NULL;
    I->N = 0;
    I->SP = NULL;

  } else {

    if(I->cullFlag && sp) {
      if (ok)
	I->V = (float *) mmalloc(sizeof(float) * I->NC * (sp->NVertTot * 31));    /* double check 31 */
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
        I->V = (float *) mmalloc(sizeof(float) * I->NC * (4 + sp->NVertTot * 6));
      else
        I->V = (float *) mmalloc(sizeof(float) * I->NC * 7);    /* one color, one alpha, one vertex per spheres */
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
        float at_sphere_scale;
        int at_sphere_color;
        float at_transp;

        AtomInfoGetSetting_f(G, ati1, cSetting_sphere_scale, sphere_scale,
                             &at_sphere_scale);
        AtomInfoGetSetting_color(G, ati1, cSetting_sphere_color, sphere_color,
                                 &at_sphere_color);
        if(AtomInfoGetSetting_f
           (G, ati1, cSetting_sphere_transparency, transp, &at_transp))
          variable_alpha = true;

        if(at_sphere_color == -1)
          c1 = *(cs->Color + a);
        else
          c1 = at_sphere_color;
        v0 = cs->Coord + 3 * a;
        vdw = ati1->vdw * at_sphere_scale + sphere_add;
        if(ColorCheckRamped(G, c1)) {
          ColorGetRamped(G, c1, v0, v, state);
          v += 3;
        } else {
          vc = ColorGet(G, c1);
          *(v++) = *(vc++);
          *(v++) = *(vc++);
          *(v++) = *(vc++);
        }

        *(v++) = 1.0F - at_transp;      /* alpha */

        if(I->cullFlag && (!spheroidFlag) && (sp)) {
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
                  float at2_sphere_scale;
                  AtomInfoType *ati2 = obj->AtomInfo + a2;
                  AtomInfoGetSetting_f(G, ati2,
                                       cSetting_sphere_scale, sphere_scale,
                                       &at2_sphere_scale);

                  if(j != a)
                    if(within3f(cs->Coord + 3 * j, v1,
                                ati2->vdw * at2_sphere_scale + sphere_add)) {
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
            restart = 1.0;      /* startin a new strip */
            for(c = 0; c < (*s); c++) {
              if(c > 1) {       /* on third vertex or better */
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
                  restart = 0.0;
                  (*nt)++;
                } else {
                  restart = 1.0;        /* next triangle is a new strip */
                }
              }
              q++;
            }
            s++;
	    ok &= !G->Interrupt;
          }
        } else if(sp) {
          q = sp->Sequence;
          s = sp->StripLen;
          if(spheroidFlag) {
            for(b = 0; ok && b < sp->NStrip; b++) {
              sphLen = cs->Spheroid + (sp->nDot * a1);
              sphNorm = cs->SpheroidNormal + (3 * sp->nDot * a1);
              for(c = 0; c < (*s); c++) {
                sphTmpN = sphNorm + 3 * (*q);
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
        } else if (ok) {                /* if sp is null, then we're simply drawing points */
          *(v++) = v0[0];
          *(v++) = v0[1];
          *(v++) = v0[2];
        }
        I->N++;
        if(nt)
          nt++;
      }
      ok &= !G->Interrupt;
      if (!ok)
	break;
    }
  }

  if(ok && sp) {                      /* don't do this if we're trying to conserve RAM */

    if(!I->LastVisib)
      I->LastVisib = Alloc(int, cs->NIndex);
    if(!I->LastColor)
      I->LastColor = Alloc(int, cs->NIndex);
    lv = I->LastVisib;
    lc = I->LastColor;
    cc = cs->Color;
    obj = cs->Obj;
    ai2 = obj->AtomInfo;
    if(sphere_color == -1)
      for(a = 0; a < cs->NIndex; a++) {
        *(lv++) = marked[cs->IdxToAtm[a]];
        *(lc++) = *(cc++);
    } else
      for(a = 0; a < cs->NIndex; a++) {
        *(lv++) = marked[cs->IdxToAtm[a]];
        *(lc++) = sphere_color;
      }
  }

  if(ok && I->V) {
    if(I->N) {
      I->V = ReallocForSure(I->V, float, (v - I->V));
      CHECKOK(ok, I->V);
      if(ok && I->NT)
        I->NT = ReallocForSure(I->NT, int, (nt - I->NT));
      CHECKOK(ok, I->NT);
    } else {
      I->V = ReallocForSure(I->V, float, 1);
      CHECKOK(ok, I->V);
      if(ok && I->NT)
        I->NT = ReallocForSure(I->NT, int, 1);
      CHECKOK(ok, I->NT);
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
  return ((void *) (struct Rep *) I);
}
