
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
#include"os_gl.h"

#include"OOMac.h"
#include"RepWireBond.h"
#include"SideChainHelper.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Setting.h"
#include"ShaderMgr.h"
#include"CGO.h"

typedef struct RepWireBond {
  Rep R;
  float *V, *VP;
  /*  Pickable *P; */
  int N, NP;
  float Width, *VarWidth;
  float Radius;
  CGO *shaderCGO;
} RepWireBond;

#include"ObjectMolecule.h"

void RepWireBondFree(RepWireBond * I);
static void RepValence(float *v, float *v1, float *v2, int *other, int a1,
                       int a2, float *coord, float *color, int ord,
                       float tube_size, int half_state, int fancy);

static void RepAromatic(float *v1, float *v2, int *other,
                        int a1, int a2, float *coord, float *color,
                        float tube_size, int half_state, float **v_ptr, int *n_ptr)
{
  float d[3], t[3], p0[3], p1[3], p2[3], *vv;
  int a3;
  float *v = *v_ptr;
  float f, f_1;
  int n = *n_ptr;
  int double_sided;

  v[0] = color[0];
  v[1] = color[1];
  v[2] = color[2];

  v[9] = color[0];
  v[10] = color[1];
  v[11] = color[2];

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  copy3f(p0, d);
  normalize3f(p0);

  /* need a prioritized third atom to get planarity */

  a3 = ObjectMoleculeGetPrioritizedOther(other, a1, a2, &double_sided);

  if(a3 < 0) {
    t[0] = p0[0];
    t[1] = p0[1];
    t[2] = -p0[2];
  } else {
    vv = coord + 3 * a3;
    t[0] = *(vv++) - v1[0];
    t[1] = *(vv++) - v1[1];
    t[2] = *(vv++) - v1[2];
    normalize3f(t);
  }

  cross_product3f(d, t, p1);

  normalize3f(p1);

  if(length3f(p1) == 0.0) {
    p1[0] = p0[1];
    p1[1] = p0[2];
    p1[2] = p0[0];
    cross_product3f(p0, p1, p2);
    normalize3f(p2);
  } else {
    cross_product3f(d, p1, p2);

    normalize3f(p2);
  }

  switch (half_state) {
  case 0:                      /* full bond */

    t[0] = p2[0] * tube_size * 2;
    t[1] = p2[1] * tube_size * 2;
    t[2] = p2[2] * tube_size * 2;

    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    v[3] = v1[0];
    v[4] = v1[1];
    v[5] = v1[2];

    v[6] = v2[0];
    v[7] = v2[1];
    v[8] = v2[2];

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];

    f = 0.14F;
    f_1 = 1.0F - f;

    v[12] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[13] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[14] = (f_1 * v1[2] + f * v2[2]) - t[2];

    f = 0.4F;
    f_1 = 1.0F - f;

    v[15] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[16] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[17] = (f_1 * v1[2] + f * v2[2]) - t[2];

    v[18] = color[0];
    v[19] = color[1];
    v[20] = color[2];

    f = 0.6F;
    f_1 = 1.0F - f;

    v[21] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[22] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[23] = (f_1 * v1[2] + f * v2[2]) - t[2];

    f = 0.86F;
    f_1 = 1.0F - f;

    v[24] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[25] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[26] = (f_1 * v1[2] + f * v2[2]) - t[2];

    v += 27;
    n += 3;

    if(double_sided) {

      v[0] = color[0];
      v[1] = color[1];
      v[2] = color[2];

      f = 0.14F;
      f_1 = 1.0F - f;

      v[3] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[4] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[5] = (f_1 * v1[2] + f * v2[2]) + t[2];

      f = 0.4F;
      f_1 = 1.0F - f;

      v[6] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[7] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[8] = (f_1 * v1[2] + f * v2[2]) + t[2];

      v[9] = color[0];
      v[10] = color[1];
      v[11] = color[2];

      f = 0.6F;
      f_1 = 1.0F - f;

      v[12] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[13] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[14] = (f_1 * v1[2] + f * v2[2]) + t[2];

      f = 0.86F;
      f_1 = 1.0F - f;

      v[15] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[16] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[17] = (f_1 * v1[2] + f * v2[2]) + t[2];

      v += 18;
      n += 2;

    }

    break;
  case 1:

    t[0] = p2[0] * tube_size * 2;
    t[1] = p2[1] * tube_size * 2;
    t[2] = p2[2] * tube_size * 2;

    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    v[3] = v1[0];
    v[4] = v1[1];
    v[5] = v1[2];

    v[6] = (v2[0] + v1[0]) / 2.0F;
    v[7] = (v2[1] + v1[1]) / 2.0F;
    v[8] = (v2[2] + v1[2]) / 2.0F;

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];

    f = 0.14F;
    f_1 = 1.0F - f;

    v[12] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[13] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[14] = (f_1 * v1[2] + f * v2[2]) - t[2];

    f = 0.4F;
    f_1 = 1.0F - f;

    v[15] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[16] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[17] = (f_1 * v1[2] + f * v2[2]) - t[2];

    v += 18;
    n += 2;

    if(double_sided) {

      v[0] = color[0];
      v[1] = color[1];
      v[2] = color[2];

      f = 0.14F;
      f_1 = 1.0F - f;

      v[3] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[4] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[5] = (f_1 * v1[2] + f * v2[2]) + t[2];

      f = 0.4F;
      f_1 = 1.0F - f;

      v[6] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[7] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[8] = (f_1 * v1[2] + f * v2[2]) + t[2];

      v += 9;
      n++;

    }
    break;
  case 2:

    t[0] = p2[0] * tube_size * 2;
    t[1] = p2[1] * tube_size * 2;
    t[2] = p2[2] * tube_size * 2;

    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    v[3] = (v2[0] + v1[0]) / 2.0F;
    v[4] = (v2[1] + v1[1]) / 2.0F;
    v[5] = (v2[2] + v1[2]) / 2.0F;

    v[6] = v2[0];
    v[7] = v2[1];
    v[8] = v2[2];

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];

    f = 0.60F;
    f_1 = 1.0F - f;

    v[12] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[13] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[14] = (f_1 * v1[2] + f * v2[2]) - t[2];

    f = 0.86F;
    f_1 = 1.0F - f;

    v[15] = (f_1 * v1[0] + f * v2[0]) - t[0];
    v[16] = (f_1 * v1[1] + f * v2[1]) - t[1];
    v[17] = (f_1 * v1[2] + f * v2[2]) - t[2];

    v += 18;
    n += 2;

    if(double_sided) {

      v[0] = color[0];
      v[1] = color[1];
      v[2] = color[2];

      f = 0.60F;
      f_1 = 1.0F - f;

      v[3] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[4] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[5] = (f_1 * v1[2] + f * v2[2]) + t[2];

      f = 0.86F;
      f_1 = 1.0F - f;

      v[6] = (f_1 * v1[0] + f * v2[0]) + t[0];
      v[7] = (f_1 * v1[1] + f * v2[1]) + t[1];
      v[8] = (f_1 * v1[2] + f * v2[2]) + t[2];

      v += 9;
      n++;

    }

    break;
  }
  *v_ptr = v;
  *n_ptr = n;

}

void RepWireBondFree(RepWireBond * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  FreeP(I->VarWidth);
  FreeP(I->VP);
  FreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}


/* lower memory use and higher performance for
   display of large trajectories, etc. */

void RepWireBondRenderImmediate(CoordSet * cs, RenderInfo * info)
{
  /* performance optimized, so it does not support the following:

     - anything other than opengl
     - display of bond valences
     - per-bond & per-atom properties, including color and line-width
     - half-bonds
     - helper settings such as cartoon_side_chain_helper
     - suppression of long bonds
     - color ramps
     - atom picking
     - display lists
     - transparency 

   */
  PyMOLGlobals *G = cs->State.G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    int active = false;
    ObjectMolecule *obj = cs->Obj;
    float line_width, line_width_setting =
      SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_line_width);
    line_width = SceneGetDynamicLineWidth(info, line_width_setting);

    if(info->width_scale_flag)
      glLineWidth(line_width * info->width_scale);
    else
      glLineWidth(line_width);

    SceneResetNormal(G, true);
    if(!info->line_lighting)
      glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    {
      int a;
      int nBond = obj->NBond;
      BondType *bd = obj->Bond;
      AtomInfoType *ai = obj->AtomInfo;
      int *atm2idx = cs->AtmToIdx;
      int discreteFlag = obj->DiscreteFlag;
      int last_color = -9;
      float *coord = cs->Coord;
      const float _pt5 = 0.5F;

      for(a = 0; a < nBond; a++) {
        int b1 = bd->index[0];
        int b2 = bd->index[1];
        AtomInfoType *ai1, *ai2;
        bd++;
        if(GET_BIT((ai1 = ai + b1)->visRep,cRepLine) && GET_BIT((ai2 = ai + b2)->visRep,cRepLine)) {
          int a1, a2;
          active = true;
          if(discreteFlag) {
            /* not optimized */
            if((cs == obj->DiscreteCSet[b1]) && (cs == obj->DiscreteCSet[b2])) {
              a1 = obj->DiscreteAtmToIdx[b1];
              a2 = obj->DiscreteAtmToIdx[b2];
            } else {
              a1 = -1;
              a2 = -1;
            }
          } else {
            a1 = atm2idx[b1];
            a2 = atm2idx[b2];
          }
          if((a1 >= 0) && (a2 >= 0)) {
            int c1 = ai1->color;
            int c2 = ai2->color;

            float *v1 = coord + 3 * a1;
            float *v2 = coord + 3 * a2;

            if(c1 == c2) {      /* same colors -> one line */
              if(c1 != last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G, c1));
              }
              glVertex3fv(v1);
              glVertex3fv(v2);  /* we done */
            } else {            /* different colors -> two lines */
              float avg[3];

              avg[0] = (v1[0] + v2[0]) * _pt5;
              avg[1] = (v1[1] + v2[1]) * _pt5;
              avg[2] = (v1[2] + v2[2]) * _pt5;

              if(c1 != last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G, c1));
              }
              glVertex3fv(v1);
              glVertex3fv(avg);

              if(c2 != last_color) {
                last_color = c2;
                glColor3fv(ColorGet(G, c2));
              }
              glVertex3fv(avg);
              glVertex3fv(v2);
            }
          }
        }
      }
    }
    glEnd();
    glEnable(GL_LIGHTING);
    if(!active)
      cs->Active[cRepLine] = false;
  }
}

static void RepWireBondRender(RepWireBond * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->R.G;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  float *v = I->V, *vw = I->VarWidth;
  float last_width = -1.0F;
  int c = I->N;
  unsigned int i, j;
  Pickable *p;
  int ok = true;
  float line_width = SceneGetDynamicLineWidth(info, I->Width);
  float line_width_setting =
    SettingGetGlobal_f(G, cSetting_line_width);
  // 0.018f is found by trial and error
  // TODO: this is not sufficient to solve the problem of disappearing cylinders
  float scale_bound = SettingGetGlobal_f(G, cSetting_field_of_view)  * cPI / 180.0f * 0.018f;

  if(ray) {

    float radius;
    float pixel_radius = ray->PixelRadius;
    if (pixel_radius < scale_bound) {
      pixel_radius = scale_bound;
    }
    if(I->Radius <= 0.0F) {
      radius = ray->PixelRadius * line_width / 2.0F;
    } else {
      vw = NULL;
      radius = I->Radius;
    }

    v = I->V;
    c = I->N;

    while(ok && c--) {
      if(vw) {
        if(last_width != *vw) {
          last_width = *vw;
          radius = ray->PixelRadius * last_width / 2.0F;
        }
        vw++;
      }
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]); */
      ok &= ray->sausage3fv(v + 3, v + 6, radius, v, v);
      v += 9;
    }

  } else if(G->HaveGUI && G->ValidContext) {
    int nvidia_bugs = SettingGetGlobal_i(G, cSetting_nvidia_bugs);

    if(pick) {

      i = (*pick)->src.index;

      v = I->VP;
      c = I->NP;
      p = I->R.P;

#ifdef PURE_OPENGL_ES_2
      (void) nvidia_bugs;
#else
      glBegin(GL_LINES);

      while(c--) {

        i++;

        if(!(*pick)[0].src.bond) {
          /* pass 1 - low order bits */

          glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8), (uchar) ((i & 0xF00) >> 4)); /* we're encoding the index into the color */
          VLACheck((*pick), Picking, i);
          p++;
          (*pick)[i].src = *p;  /* copy object and atom info */
          (*pick)[i].context = I->R.context;

        } else {
          /* pass 2 - high order bits */

          j = i >> 12;

          glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                     (uchar) ((j & 0xF00) >> 4));

        }
        if(nvidia_bugs) {
          glFlush();
        }
        glVertex3fv(v);
        v += 3;
        glVertex3fv(v);
        v += 3;

      }
      glEnd();
#endif
      (*pick)[0].src.index = i; /* pass the count */
    } else { /* else not pick i.e., when rendering */
      short use_shader, generate_shader_cgo = 0;
      short line_as_cylinders ;
      int nvidia_bugs = SettingGetGlobal_i(G, cSetting_nvidia_bugs);
      use_shader = SettingGetGlobal_b(G, cSetting_line_use_shader) & 
                   SettingGetGlobal_b(G, cSetting_use_shaders);

      line_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_line_as_cylinders);
      if (!use_shader && I->shaderCGO){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (I->shaderCGO && (line_as_cylinders ^ I->shaderCGO->has_draw_cylinder_buffers)){
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
	  CShaderPrg *shaderPrg;
	  if (line_as_cylinders){
	    // vertex scale is bound so that cylinders cannot disappear when it gets too low
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if (!shaderPrg) return;
	    if (vw){
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width_setting/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    if (!shaderPrg) return;
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  
	  CShaderPrg_Disable(shaderPrg);
	  return;
	}
      }

      v = I->V;
      c = I->N;

      if (ok && generate_shader_cgo){
	ok &= CGOLinewidthSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE);
	
	if(ok && !info->line_lighting)
	  ok &= CGODisable(I->shaderCGO, GL_LIGHTING);
	ok &= CGOResetNormal(I->shaderCGO, true);
      } else {

	if(info->width_scale_flag)
	  glLineWidth(line_width * info->width_scale);
	else
	  glLineWidth(line_width);
	
	if(!info->line_lighting)
	  glDisable(GL_LIGHTING);
	SceneResetNormal(G, true);
      }

      
      if (generate_shader_cgo){
	float curColor[3];
	while(ok && c--) {
	  //	  float cylinder_width = line_width;
	  float cylinder_width = line_width_setting;
	  if(vw) {
	    if(last_width != *vw) {
	      last_width = *vw;
	      ok &= CGOLinewidth(I->shaderCGO, last_width);
	    }
	    cylinder_width = *vw;
	    vw++;
	  }
	  if (ok){
	    ok &= CGOColorv(I->shaderCGO, v);
	    copy3f(v, curColor);
	    v += 3;
	  }
	  if (ok){
	    if (line_as_cylinders){
	      float *origin, axis[3];
	      origin = v;
	      v += 3;
	      axis[0] = v[0] - origin[0];
	      axis[1] = v[1] - origin[1];
	      axis[2] = v[2] - origin[2];
	      v += 3;

	      {
		if (c && equal3f(&v[-3], &v[3]) && !equal3f(curColor, v)){
		  /* if successive bonds share midpoint, then draw one cylinder with two colors */
		  origin = &v[-6];
		  axis[0] = v[6] - origin[0];
		  axis[1] = v[7] - origin[1];
		  axis[2] = v[8] - origin[2];
		  // if next bond has same half-point, then draw one cylinder with two colors
		  ok &= CGOShaderCylinder2ndColor(I->shaderCGO, origin, axis, cylinder_width/line_width_setting, 15, v);
		  v += 9;
		  c--;
		} else {
		  /* Storing the cylinder_width divided by the current line_width setting */
		  ok &= CGOShaderCylinder(I->shaderCGO, origin, axis, cylinder_width/line_width_setting, 15);
		}
	      }
	    } else {
	      ok &= CGOBegin(I->shaderCGO, GL_LINES);
	      if (ok){
		ok &= CGOVertexv(I->shaderCGO, v);
		v += 3;
	      }
	      if (ok){
		ok &= CGOVertexv(I->shaderCGO, v);
		v += 3;
	      }
	      if (ok)
		ok &= CGOEnd(I->shaderCGO);
	    }
	  }
	}
      } else {
	while(c--) {
	  if(vw) {
	    if(last_width != *vw) {
	      last_width = *vw;
	      glLineWidth(last_width);
	    }
	    vw++;
	  }
#ifdef PURE_OPENGL_ES_2
    /* TODO */
#else
	  glBegin(GL_LINES);
	  glColor3fv(v);
	  v += 3;
	  if(nvidia_bugs) {
	    glFlush();
	  }
	  glVertex3fv(v);
	  v += 3;
	  glVertex3fv(v);
	  v += 3;
	  glEnd();
#endif
	}
      }
      if (generate_shader_cgo){
	if (ok)
	  ok &= CGOEnable(I->shaderCGO, GL_LIGHTING);
      } else {
	glEnable(GL_LIGHTING);
      }
      if (use_shader) {
	if (ok && generate_shader_cgo){
	  CGO *convertcgo = NULL;
	  if (ok)
	    ok &= CGOStop(I->shaderCGO);
	  if (ok){
	    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);    
	    CGOFree(I->shaderCGO);    
	    I->shaderCGO = convertcgo;
	    CHECKOK(ok, I->shaderCGO);
	    convertcgo = NULL;
	  }
	  if (ok && I->shaderCGO){
	    if (line_as_cylinders){
              convertcgo = CGOOptimizeGLSLCylindersToVBOIndexed(I->shaderCGO, 0);
	    } else {
              convertcgo = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0);
	    }
	  }
      CGOFree(I->shaderCGO);
      I->shaderCGO = convertcgo;
      CHECKOK(ok, I->shaderCGO);
	}
	
	if (ok){
	  CShaderPrg *shaderPrg;
	  if (line_as_cylinders){
	    // vertex scale is bound so that cylinders cannot disappear when it gets too low
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if (!shaderPrg) return;
	    if (vw){
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width_setting/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    if (!shaderPrg) return;
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }	 
	  
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  
	  CShaderPrg_Disable(shaderPrg);
	}
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
    I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
    I->R.cs->Active[cRepLine] = false;
  }
}

static
bool IsBondTerminal(ObjectMolecule *obj, int b1, int b2){
  int *neighbor = obj->Neighbor;
  if(neighbor) {
    int mem, nbr;
    int heavy1 = 0, heavy2 = 0;
    AtomInfoType *atomInfo = obj->AtomInfo;
    nbr = neighbor[b1] + 1;
    while(((mem = neighbor[nbr]) >= 0)) {
      if(atomInfo[mem].protons > 1) {
        heavy1++;
      }
      nbr += 2;
    }
    nbr = neighbor[b2] + 1;
    while(((mem = neighbor[nbr]) >= 0)) {
      if(atomInfo[mem].protons > 1) {
        heavy2++;
      }
      nbr += 2;
    }
    if((heavy1 < 2) || (heavy2 < 2))
      return true;
  }
  return false;
}

Rep *RepWireBondNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj = cs->Obj;
  int a1, a2, b1, b2;
  int a, c1, c2, s1, s2, ord;
  BondType *b;
  int half_bonds, *other = NULL;
  float valence;
  float *v, *v0, *v1, *v2, h[3];
  int visFlag;
  int maxSegment = 0;
  int maxBond = 0;
  float tmpColor[3];
  float line_width;
  int valence_flag = false;
  Pickable *rp;
  AtomInfoType *ai1, *ai2;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 0;
  int line_stick_helper = 0;
  int na_mode;
  bool *marked = NULL;
  int valence_found = false;
  int variable_width = false;
  int n_line_width = 0;
  int line_color;
  int hide_long = false;
  int fancy;
  const float _0p9 = 0.9F;
  int ok = true;

  OOAlloc(G, RepWireBond);
  CHECKOK(ok, I);
  PRINTFD(G, FB_RepWireBond)
    " RepWireBondNew-Debug: entered.\n" ENDFD;

  visFlag = false;
  b = obj->Bond;
  if(ok && GET_BIT(obj->RepVisCache,cRepLine)){
    for(a = 0; a < obj->NBond; a++) {
      b1 = b->index[0];
      b2 = b->index[1];
      b++;
      if(GET_BIT(obj->AtomInfo[b1].visRep,cRepLine) || GET_BIT(obj->AtomInfo[b2].visRep,cRepLine)) {
        visFlag = true;
        break;
      }
    }
  }
  if(!visFlag) {
    OOFreeP(I);
    return (NULL);              /* skip if no dots are visible */
  }
  marked = Calloc(bool, obj->NAtom);
  CHECKOK(ok, marked);
  if (!ok){
    OOFreeP(I);
    return (NULL);
  }
  
  valence_flag = SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_valence);
  valence = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_valence_size);
  cartoon_side_chain_helper = SettingGet_b(G, cs->Setting, obj->Obj.Setting,
					   cSetting_cartoon_side_chain_helper);
  ribbon_side_chain_helper = SettingGet_b(G, cs->Setting, obj->Obj.Setting,
					  cSetting_ribbon_side_chain_helper);
  line_stick_helper = SettingGet_b(G, cs->Setting, obj->Obj.Setting,
				   cSetting_line_stick_helper);
  line_color = SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_line_color);
  line_width = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_line_width);
  
  if(line_stick_helper && (SettingGet_f(G, cs->Setting, obj->Obj.Setting,
					cSetting_stick_transparency) > R_SMALL4))
    line_stick_helper = false;
  half_bonds = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_half_bonds);
  hide_long = SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_hide_long_bonds);
  na_mode =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_nucleic_acid_mode);
  int na_mode_ribbon =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_nucleic_acid_mode);
  fancy = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_valence_mode) == 1;
  
  b = obj->Bond;
  
  for(a = 0; a < obj->NBond; a++) {
    b1 = b->index[0];
    b2 = b->index[1];

    if(obj->DiscreteFlag) {
      if((cs == obj->DiscreteCSet[b1]) && (cs == obj->DiscreteCSet[b2])) {
        a1 = obj->DiscreteAtmToIdx[b1];
        a2 = obj->DiscreteAtmToIdx[b2];
      } else {
        a1 = -1;
        a2 = -1;
      }
    } else {
      a1 = cs->AtmToIdx[b1];
      a2 = cs->AtmToIdx[b2];
    }
    if((a1 >= 0) && (a2 >= 0)) {
      if((!variable_width) && AtomInfoCheckBondSetting(G, b, cSetting_line_width))
        variable_width = true;
      auto bd_valence_flag = BondSettingGetWD(G, b, cSetting_valence, valence_flag);
      if(bd_valence_flag) {

        valence_found = true;
        if((b->order > 0) && (b->order < 4)) {
          maxSegment += 2 * b->order;
        } else if(b->order == 4) {      /* aromatic */
          maxSegment += 10;
        } else {
          maxSegment += 2;
        }
      } else
        maxSegment += 2;
      maxBond++;
    }
    b++;
  }

  RepInit(G, &I->R);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo * info)) RepWireBondRender;
  I->R.fFree = (void (*)(struct Rep *)) RepWireBondFree;
  I->Width = line_width;
  I->Radius = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_line_radius);

  I->shaderCGO = 0;
  I->N = 0;
  I->NP = 0;
  I->V = NULL;
  I->VP = NULL;
  I->VarWidth = NULL;
  I->R.P = NULL;
  I->R.fRecolor = NULL;
  I->R.context.object = (void *) obj;
  I->R.context.state = state;
  I->R.cs = cs;

  if(obj->NBond) {

    if(valence_found)           /* build list of up to 2 connected atoms for each atom */
      other = ObjectMoleculeGetPrioritizedOtherIndexList(obj, cs);

    if(variable_width) {
      I->VarWidth = Alloc(float, maxSegment);
      CHECKOK(ok, I->VarWidth);
    }

    if (ok)
      I->V = Alloc(float, maxSegment * 9);
    CHECKOK(ok, I->V);

    if(ok && (cartoon_side_chain_helper || ribbon_side_chain_helper)) {
      SideChainHelperMarkNonCartoonBonded(marked, obj, cs,
          cartoon_side_chain_helper,
          ribbon_side_chain_helper);
    }

    v = I->V;
    b = obj->Bond;

    for(a = 0; ok && a < obj->NBond; a++) {

      b1 = b->index[0];
      b2 = b->index[1];
      ord = b->order;

      /*
         b1 = *(b++);
         b2 = *(b++);
         ord = (*(b++));
       */
      if(obj->DiscreteFlag) {
        if((cs == obj->DiscreteCSet[b1]) && (cs == obj->DiscreteCSet[b2])) {
          a1 = obj->DiscreteAtmToIdx[b1];
          a2 = obj->DiscreteAtmToIdx[b2];
        } else {
          a1 = -1;
          a2 = -1;
        }
      } else {
        a1 = cs->AtmToIdx[b1];
        a2 = cs->AtmToIdx[b2];
      }
      if((a1 >= 0) && (a2 >= 0)) {

        AtomInfoType *ati1 = obj->AtomInfo + b1;
        AtomInfoType *ati2 = obj->AtomInfo + b2;

        s1 = GET_BIT(ati1->visRep,cRepLine);
        s2 = GET_BIT(ati2->visRep,cRepLine);

        if((s1 || s2) && !(s1 && s2))
          if(!half_bonds) {
            if(line_stick_helper &&
               (((!s1) && GET_BIT(ati1->visRep,cRepCyl) && (!GET_BIT(ati2->visRep,cRepCyl))) ||
                ((!s2) && GET_BIT(ati2->visRep,cRepCyl) && (!GET_BIT(ati1->visRep,cRepCyl)))))
              s1 = s2 = 1;      /* turn on line when both stick and line are alternately shown */
            else {
              s1 = 0;
              s2 = 0;
            }
          }

        if(hide_long && (s1 || s2)) {
          float cutoff = (ati1->vdw + ati2->vdw) * _0p9;
          v1 = cs->Coord + 3 * a1;
          v2 = cs->Coord + 3 * a2;
          ai1 = obj->AtomInfo + b1;
          ai2 = obj->AtomInfo + b2;
          if(!within3f(v1, v2, cutoff)) /* atoms separated by more than 90% of the sum of their vdw radii */
            s1 = s2 = 0;
        }

        if(s1 || s2) {
          float bd_line_width = line_width;
          int terminal = false;

          auto bd_valence_flag = BondSettingGetWD(G, b, cSetting_valence, valence_flag);
          auto bd_line_color = BondSettingGetWD(G, b, cSetting_line_color, line_color);

          if(fancy && bd_valence_flag && (b->order > 1)) {
            terminal = IsBondTerminal(obj, b1, b2);
          }

          if(variable_width) {
            bd_line_width = BondSettingGetWD(G, b, cSetting_line_width, line_width);
          }

          if(bd_line_color < 0) {
            if(bd_line_color == cColorObject) {
              c1 = (c2 = obj->Obj.Color);
            } else if(ColorCheckRamped(G, bd_line_color)) {
              c1 = (c2 = bd_line_color);
            } else {
              c1 = ati1->color;
              c2 = ati2->color;
            }
          } else {
            c1 = (c2 = bd_line_color);
          }

          v1 = cs->Coord + 3 * a1;
          v2 = cs->Coord + 3 * a2;

          if (line_stick_helper && (ati1->visRep & ati2->visRep & cRepCylBit)) {
            s1 = s2 = 0;
          } else if ((ati1->flags & ati2->flags & cAtomFlag_polymer)) {
            // side chain helpers
            if ((cRepCartoonBit & ati1->visRep & ati2->visRep)) {
              bool sc_helper = AtomSettingGetWD(G, ati1,
                  cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper);

              if (!sc_helper)
                AtomSettingGetIfDefined(G, ati2, cSetting_cartoon_side_chain_helper, &sc_helper);

              if (sc_helper &&
                  SideChainHelperFilterBond(G, marked, ati1, ati2, b1, b2, na_mode, &c1, &c2))
                s1 = s2 = 0;
            }

            if ((s1 || s2) && (cRepRibbonBit & ati1->visRep & ati2->visRep)) {
              bool sc_helper = AtomSettingGetWD(G, ati1,
                  cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper);

              if (!sc_helper)
                AtomSettingGetIfDefined(G, ati2, cSetting_ribbon_side_chain_helper, &sc_helper);

              if (sc_helper &&
                  SideChainHelperFilterBond(G, marked, ati1, ati2, b1, b2, na_mode_ribbon, &c1, &c2))
                s1 = s2 = 0;
            }
          }

          if((c1 == c2) && s1 && s2 && (!ColorCheckRamped(G, c1))) {

            v0 = ColorGet(G, c1);

            if((bd_valence_flag) && (ord > 1) && (ord < 4)) {
              RepValence(v, v1, v2, other, a1, a2, cs->Coord, v0, ord, valence, 0, fancy
                         && !terminal);
              v += ord * 9;
              I->N += ord;
            } else if(bd_valence_flag && (ord == 4)) {  /* aromatic */
              RepAromatic(v1, v2, other, a1, a2, cs->Coord, v0, valence, 0, &v, &I->N);
            } else {
              I->N++;
              *(v++) = *(v0++);
              *(v++) = *(v0++);
              *(v++) = *(v0++);

              *(v++) = *(v1++);
              *(v++) = *(v1++);
              *(v++) = *(v1++);

              *(v++) = *(v2++);
              *(v++) = *(v2++);
              *(v++) = *(v2++);
            }
          } else {

            h[0] = (v1[0] + v2[0]) / 2;
            h[1] = (v1[1] + v2[1]) / 2;
            h[2] = (v1[2] + v2[2]) / 2;

            if(s1) {

              if(ColorCheckRamped(G, c1)) {
                ColorGetRamped(G, c1, v1, tmpColor, state);
                v0 = tmpColor;
              } else {
                v0 = ColorGet(G, c1);
              }

              if((bd_valence_flag) && (ord > 1) && (ord < 4)) {
                RepValence(v, v1, h, other, a1, a2, cs->Coord, v0, ord, valence, 1, fancy
                           && !terminal);
                v += ord * 9;
                I->N += ord;
              } else if(bd_valence_flag && (ord == 4)) {
                RepAromatic(v1, v2, other, a1, a2, cs->Coord, v0, valence, 1, &v, &I->N);
              } else {

                I->N++;
                *(v++) = *(v0++);
                *(v++) = *(v0++);
                *(v++) = *(v0++);

                *(v++) = *(v1++);
                *(v++) = *(v1++);
                *(v++) = *(v1++);

                *(v++) = h[0];
                *(v++) = h[1];
                *(v++) = h[2];
              }
            }
            if(s2) {
              if(ColorCheckRamped(G, c2)) {
                ColorGetRamped(G, c2, v2, tmpColor, state);
                v0 = tmpColor;
              } else {
                v0 = ColorGet(G, c2);
              }
              if((bd_valence_flag) && (ord > 1) && (ord < 4)) {
                RepValence(v, h, v2, other, a1, a2, cs->Coord, v0, ord, valence, 2, fancy
                           && !terminal);
                v += ord * 9;
                I->N += ord;
              } else if(bd_valence_flag && (ord == 4)) {
                RepAromatic(v1, v2, other, a1, a2, cs->Coord, v0, valence, 2, &v, &I->N);
              } else {
                I->N++;
                *(v++) = *(v0++);
                *(v++) = *(v0++);
                *(v++) = *(v0++);

                *(v++) = h[0];
                *(v++) = h[1];
                *(v++) = h[2];

                *(v++) = *(v2++);
                *(v++) = *(v2++);
                *(v++) = *(v2++);
              }

            }
          }

          /* record effective line_widths for these segments */

          if(variable_width) {
            while(n_line_width < I->N) {
              I->VarWidth[n_line_width] = bd_line_width;
              n_line_width++;
            }
          }
        }
      }
      b++;
      ok &= !G->Interrupt;
    }
    if (ok)
      I->V = ReallocForSure(I->V, float, (v - I->V));
    CHECKOK(ok, I->V);
    if(ok && I->VarWidth) {
      I->VarWidth = ReallocForSure(I->VarWidth, float, n_line_width);
      CHECKOK(ok, I->VarWidth);
    }

    /* now create pickable verson */

    if(ok && SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_pickable)) {
      I->VP = Alloc(float, maxBond * 6 * 2);
      CHECKOK(ok, I->VP);

      if (ok)
	I->R.P = Alloc(Pickable, 2 * maxBond + 1);
      CHECKOK(ok, I->R.P);
      if (ok){
	rp = I->R.P + 1;          /* skip first record! */

	v = I->VP;
	b = obj->Bond;
      }
      for(a = 0; ok && a < obj->NBond; a++) {

        b1 = b->index[0];
        b2 = b->index[1];
        b++;
        if(obj->DiscreteFlag) {
          if((cs == obj->DiscreteCSet[b1]) && (cs == obj->DiscreteCSet[b2])) {
            a1 = obj->DiscreteAtmToIdx[b1];
            a2 = obj->DiscreteAtmToIdx[b2];
          } else {
            a1 = -1;
            a2 = -1;
          }
        } else {
          a1 = cs->AtmToIdx[b1];
          a2 = cs->AtmToIdx[b2];
        }
        if((a1 >= 0) && (a2 >= 0)) {

          ai1 = obj->AtomInfo + b1;
          ai2 = obj->AtomInfo + b2;
          s1 = GET_BIT(ai1->visRep,cRepLine);
          s2 = GET_BIT(ai2->visRep,cRepLine);

          if(!(s1 && s2)) {
            if(!half_bonds) {
              s1 = 0;
              s2 = 0;
            }
          }

          if(hide_long && (s1 || s2)) {
            float cutoff = (ai1->vdw + ai2->vdw) * _0p9;
            v1 = cs->Coord + 3 * a1;
            v2 = cs->Coord + 3 * a2;
            ai1 = obj->AtomInfo + b1;
            ai2 = obj->AtomInfo + b2;
            if(!within3f(v1, v2, cutoff))       /* atoms separated by more than 90% of the sum of their vdw radii */
              s1 = s2 = 0;
          }

          if(s1 || s2) {

            v1 = cs->Coord + 3 * a1;
            v2 = cs->Coord + 3 * a2;

            h[0] = (v1[0] + v2[0]) / 2;
            h[1] = (v1[1] + v2[1]) / 2;
            h[2] = (v1[2] + v2[2]) / 2;

            if(s1 & (!ai1->masked)) {

              I->NP++;
              rp->index = b1;
              rp->bond = a;
              rp++;

              *(v++) = *(v1++);
              *(v++) = *(v1++);
              *(v++) = *(v1++);

              *(v++) = h[0];
              *(v++) = h[1];
              *(v++) = h[2];
            }
            if(s2 & (!ai2->masked)) {

              I->NP++;
              rp->index = b2;
              rp->bond = a;
              rp++;

              *(v++) = h[0];
              *(v++) = h[1];
              *(v++) = h[2];

              *(v++) = *(v2++);
              *(v++) = *(v2++);
              *(v++) = *(v2++);
            }
          }
        }
	ok &= !G->Interrupt;
      }
      if (ok){
	I->R.P = Realloc(I->R.P, Pickable, I->NP + 1);
	CHECKOK(ok, I->R.P);
	if (ok)
	  I->R.P[0].index = I->NP;
      }
      if (ok)
	I->VP = ReallocForSure(I->VP, float, (v - I->VP));
      CHECKOK(ok, I->VP);
    }
  }
  FreeP(marked);
  FreeP(other);
  if (!ok){
    RepWireBondFree(I);
    I = NULL;
  }
  return (Rep *) I;
}

static void RepValence(float *v, float *v1, float *v2, int *other,
                       int a1, int a2, float *coord, float *color, int ord,
                       float tube_size, int half_state, int fancy)
{

  float d[3], t[3], p0[3], p1[3], p2[3], *vv;
  int a3;
  const float indent = tube_size;

  v[0] = color[0];
  v[1] = color[1];
  v[2] = color[2];

  v[9] = color[0];
  v[10] = color[1];
  v[11] = color[2];

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  copy3f(p0, d);
  normalize3f(p0);

  /* need a prioritized third atom to get planarity */

  a3 = ObjectMoleculeGetPrioritizedOther(other, a1, a2, NULL);

  if(a3 < 0) {
    t[0] = p0[0];
    t[1] = p0[1];
    t[2] = -p0[2];
  } else {
    vv = coord + 3 * a3;
    t[0] = *(vv++) - v1[0];
    t[1] = *(vv++) - v1[1];
    t[2] = *(vv++) - v1[2];
    normalize3f(t);
  }

  cross_product3f(d, t, p1);

  normalize3f(p1);

  if(length3f(p1) == 0.0) {
    p1[0] = p0[1];
    p1[1] = p0[2];
    p1[2] = p0[0];
    cross_product3f(p0, p1, p2);
    normalize3f(p2);
  } else {
    cross_product3f(d, p1, p2);

    normalize3f(p2);
  }

  /* now we have a coordinate system */

  t[0] = p2[0] * tube_size;
  t[1] = p2[1] * tube_size;
  t[2] = p2[2] * tube_size;

  switch (ord) {
  case 2:
    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];

    if(fancy) {

      float f, f_1;
      v[3] = v1[0];
      v[4] = v1[1];
      v[5] = v1[2];

      v[6] = v2[0];
      v[7] = v2[1];
      v[8] = v2[2];

      if(half_state == 2) {
        v[12] = v1[0] - 2 * t[0];
        v[13] = v1[1] - 2 * t[1];
        v[14] = v1[2] - 2 * t[2];
      } else {
        if(half_state == 1)
          f = indent * 2;
        else
          f = indent;

        f_1 = 1.0F - f;

        v[12] = (f_1 * v1[0] + f * v2[0]) - 2 * t[0];
        v[13] = (f_1 * v1[1] + f * v2[1]) - 2 * t[1];
        v[14] = (f_1 * v1[2] + f * v2[2]) - 2 * t[2];
      }

      if(half_state == 1) {
        v[15] = v2[0] - 2 * t[0];
        v[16] = v2[1] - 2 * t[1];
        v[17] = v2[2] - 2 * t[2];

      } else {
        if(half_state == 2)
          f = 1.0 - 2 * indent;
        else
          f = 1.0 - indent;
        f_1 = 1.0F - f;
        v[15] = (f_1 * v1[0] + f * v2[0]) - 2 * t[0];
        v[16] = (f_1 * v1[1] + f * v2[1]) - 2 * t[1];
        v[17] = (f_1 * v1[2] + f * v2[2]) - 2 * t[2];
      }

    } else {
      v[3] = v1[0] - t[0];
      v[4] = v1[1] - t[1];
      v[5] = v1[2] - t[2];

      v[6] = v2[0] - t[0];
      v[7] = v2[1] - t[1];
      v[8] = v2[2] - t[2];

      v[12] = v1[0] + t[0];
      v[13] = v1[1] + t[1];
      v[14] = v1[2] + t[2];

      v[15] = v2[0] + t[0];
      v[16] = v2[1] + t[1];
      v[17] = v2[2] + t[2];
    }
    break;
  case 3:
    t[0] = t[0] * 2;
    t[1] = t[1] * 2;
    t[2] = t[2] * 2;

    v[0] = color[0];
    v[1] = color[1];
    v[2] = color[2];

    if(fancy) {
      float f, f_1;

      if(half_state == 2) {
        v[3] = v1[0] - t[0];
        v[4] = v1[1] - t[1];
        v[5] = v1[2] - t[2];
      } else {
        if(half_state == 1)
          f = indent * 2;
        else
          f = indent;

        f_1 = 1.0F - f;

        v[3] = (f_1 * v1[0] + f * v2[0]) - t[0];
        v[4] = (f_1 * v1[1] + f * v2[1]) - t[1];
        v[5] = (f_1 * v1[2] + f * v2[2]) - t[2];
      }

      if(half_state == 1) {
        v[6] = v2[0] - t[0];
        v[7] = v2[1] - t[1];
        v[8] = v2[2] - t[2];

      } else {
        if(half_state == 2)
          f = 1.0 - 2 * indent;
        else
          f = 1.0 - indent;
        f_1 = 1.0F - f;
        v[6] = (f_1 * v1[0] + f * v2[0]) - t[0];
        v[7] = (f_1 * v1[1] + f * v2[1]) - t[1];
        v[8] = (f_1 * v1[2] + f * v2[2]) - t[2];
      }

      if(half_state == 2) {
        v[12] = v1[0] + t[0];
        v[13] = v1[1] + t[1];
        v[14] = v1[2] + t[2];
      } else {
        if(half_state == 1)
          f = indent * 2;
        else
          f = indent;

        f_1 = 1.0F - f;

        v[12] = (f_1 * v1[0] + f * v2[0]) + t[0];
        v[13] = (f_1 * v1[1] + f * v2[1]) + t[1];
        v[14] = (f_1 * v1[2] + f * v2[2]) + t[2];
      }

      if(half_state == 1) {
        v[15] = v2[0] + t[0];
        v[16] = v2[1] + t[1];
        v[17] = v2[2] + t[2];
      } else {
        if(half_state == 2)
          f = 1.0 - 2 * indent;
        else
          f = 1.0 - indent;
        f_1 = 1.0F - f;
        v[15] = (f_1 * v1[0] + f * v2[0]) + t[0];
        v[16] = (f_1 * v1[1] + f * v2[1]) + t[1];
        v[17] = (f_1 * v1[2] + f * v2[2]) + t[2];
      }
    } else {

      v[3] = v1[0] - t[0];
      v[4] = v1[1] - t[1];
      v[5] = v1[2] - t[2];

      v[6] = v2[0] - t[0];
      v[7] = v2[1] - t[1];
      v[8] = v2[2] - t[2];

      v[12] = v1[0] + t[0];
      v[13] = v1[1] + t[1];
      v[14] = v1[2] + t[2];

      v[15] = v2[0] + t[0];
      v[16] = v2[1] + t[1];
      v[17] = v2[2] + t[2];

    }

    v[9] = color[0];
    v[10] = color[1];
    v[11] = color[2];

    v[18] = color[0];
    v[19] = color[1];
    v[20] = color[2];

    v[21] = v1[0];
    v[22] = v1[1];
    v[23] = v1[2];

    v[24] = v2[0];
    v[25] = v2[1];
    v[26] = v2[2];
    break;
  }
}
