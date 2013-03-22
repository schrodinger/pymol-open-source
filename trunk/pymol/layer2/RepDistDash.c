
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

#include"OOMac.h"
#include"RepDistDash.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"
#include"CGO.h"
#include"ShaderMgr.h"
#include"CoordSet.h"

typedef struct RepDistDash {
  Rep R;
  float *V;
  int N;
  CObject *Obj;
  DistSet *ds;
  float linewidth, radius;
  CGO *shaderCGO;
} RepDistDash;

#include"ObjectDist.h"

void RepDistDashFree(RepDistDash * I);

void RepDistDashFree(RepDistDash * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  VLAFreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepDistDashRender(RepDistDash * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  float *vc;
  int round_ends;
  int ok = true;
  int color =
    SettingGet_color(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_color);
  float line_width =
    SettingGet_f(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_width);

  I->radius =
    SettingGet_f(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_radius);
  round_ends =
    SettingGet_b(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_round_ends);
  line_width = SceneGetDynamicLineWidth(info, line_width);

  if(ray) {

    float radius;

    if(I->radius <= 0.0F) {
      radius = ray->PixelRadius * line_width / 2.0F;
    } else {
      radius = I->radius;
    }

    if(color < 0)
      color = I->Obj->Color;
    vc = ColorGet(G, color);
    v = I->V;
    c = I->N;

    while(ok && c > 0) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]); */
      if(round_ends) {
        ok &= ray->fSausage3fv(ray, v, v + 3, radius, vc, vc);
      } else {
        ok &= ray->fCustomCylinder3fv(ray, v, v + 3, radius, vc, vc, cCylCapFlat, cCylCapFlat);
      }
      v += 6;
      c -= 2;
    }

  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      short use_shader, generate_shader_cgo = 0, use_display_lists = 0, dash_as_cylinders = 0;

      use_shader = (int) SettingGet(G, cSetting_dash_use_shader) & 
                           (int) SettingGet(G, cSetting_use_shaders);
      use_display_lists = (int) SettingGet(G, cSetting_use_display_lists);
      dash_as_cylinders = (int) SettingGet(G, cSetting_render_as_cylinders) && SettingGet(G, cSetting_dash_as_cylinders);

      if (!use_shader && I->shaderCGO){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }
      if (I->shaderCGO && (dash_as_cylinders ^ I->shaderCGO->has_draw_cylinder_buffers)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
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
	  if (ok)
	    I->shaderCGO->use_shader = true;
	  generate_shader_cgo = 1;
	} else if (ok) {
	  CShaderPrg *shaderPrg;
	  if (dash_as_cylinders){
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if(I->radius == 0.0F) {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", I->radius);
	    }
	    if (!round_ends){
	      CShaderPrg_Set1f(shaderPrg, "no_flat_caps", 0.f);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }
	  if (!shaderPrg) return;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);

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
#else
      (void) use_display_lists;
#endif

      if (generate_shader_cgo){
	if (ok)
	  ok &= CGOLinewidthSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
	if (ok)
	  ok &= CGOResetNormal(I->shaderCGO, true);
      } else {
	if(info->width_scale_flag) {
	  glLineWidth(line_width * info->width_scale);
	} else {
	  glLineWidth(line_width);
	}
        SceneResetNormal(G, true);
      }

      if (generate_shader_cgo){
	if (ok){
	  if(color >= 0){
	    ok &= CGOColorv(I->shaderCGO, ColorGet(G, color));
	  } else if (I->Obj && I->Obj->Color >= 0){
	    ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->Obj->Color));
	  }
	}
	v = I->V;
	c = I->N;
	if (dash_as_cylinders){
	  float *origin = NULL, axis[3];
	  while(ok && c > 0) {
	    origin = v;
	    v += 3;
	    axis[0] = v[0] - origin[0];
	    axis[1] = v[1] - origin[1];
	    axis[2] = v[2] - origin[2];
	    v += 3;
	    ok &= CGOShaderCylinder(I->shaderCGO, origin, axis, 1.f, 15);
	    c -= 2;
	  }
	} else {
	  ok &= CGOBegin(I->shaderCGO, GL_LINES);
	  while(ok && c > 0) {
	    ok &= CGOVertexv(I->shaderCGO, v);
	    v += 3;
	    if (ok)
	      ok &= CGOVertexv(I->shaderCGO, v);
	    v += 3;
	    c -= 2;
	  }
	  if (ok)
	    ok &= CGOEnd(I->shaderCGO);
	}
      } else {
	if(color >= 0){
	  glColor3fv(ColorGet(G, color));
	}
	v = I->V;
	c = I->N;
	if(!info->line_lighting)
	  glDisable(GL_LIGHTING);
#ifdef _PYMOL_GL_DRAWARRAYS
	{
	  int nverts = c, pl;	  
	  ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
	  pl = 0;
	  while(c > 0) {
	    vertVals[pl++] = v[0]; vertVals[pl++] = v[1]; vertVals[pl++] = v[2];
	    v += 3;
	    vertVals[pl++] = v[0]; vertVals[pl++] = v[1]; vertVals[pl++] = v[2];
	    v += 3;
	    c -= 2;
	  }
	  glEnableClientState(GL_VERTEX_ARRAY);
	  glVertexPointer(3, GL_FLOAT, 0, vertVals);
	  glDrawArrays(GL_LINES, 0, nverts);
	  glDisableClientState(GL_VERTEX_ARRAY);
	  DEALLOCATE_ARRAY(vertVals)
	}
#else
	glBegin(GL_LINES);
	while(c > 0) {
	  glVertex3fv(v);
	  v += 3;
	  glVertex3fv(v);
	  v += 3;
	  c -= 2;
	}
	glEnd();
#endif
	glEnable(GL_LIGHTING);
      }
      if (use_shader) {
	if (generate_shader_cgo){
	  CGO *convertcgo = NULL;
	  if (ok)
	    ok &= CGOStop(I->shaderCGO);
#ifdef _PYMOL_CGO_DRAWARRAYS
	  if (ok)
	    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);    
	  CHECKOK(ok, convertcgo);
	  CGOFree(I->shaderCGO);    
	  I->shaderCGO = convertcgo;
	  convertcgo = NULL;
#else
	  (void)convertcgo;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
	  if (ok){
	    if (dash_as_cylinders){
	      convertcgo = CGOOptimizeGLSLCylindersToVBOIndexed(I->shaderCGO, 0);
	    } else {
	      convertcgo = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0);
	    }
	    CHECKOK(ok, convertcgo);
	  }
	  if (convertcgo){
	    CGOFree(I->shaderCGO);
	    I->shaderCGO = convertcgo;
	    convertcgo = NULL;
	  }
#else
	  (void)convertcgo;
#endif
	}
	
	if (ok) {
	  CShaderPrg *shaderPrg;
	  if (dash_as_cylinders){
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if(I->radius == 0.0F) {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", I->radius);
	    }
	    if (!round_ends){
	      CShaderPrg_Set1f(shaderPrg, "no_flat_caps", 0.f);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }	 

	  if (!shaderPrg)
	    return;
        
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  
	  CShaderPrg_Disable(shaderPrg);
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
    I->ds->Rep[cRepDash] = NULL;
    RepDistDashFree(I);
  }
}

Rep *RepDistDashNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->State.G;
  int a;
  int n;
  float *v, *v1, *v2, d[3], d1[3];
  float l;
  float dash_len, dash_gap, dash_sum;
  int ok = true;

  OOAlloc(G, RepDistDash);
  CHECKOK(ok, I);

  if(!ok || !ds->NIndex) {
    OOFreeP(I);
    return (NULL);
  }

  RepInit(G, &I->R);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepDistDashRender;
  I->R.fFree = (void (*)(struct Rep *)) RepDistDashFree;
  I->R.fRecolor = NULL;
  I->R.context.state = state;
  dash_len = SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dash_length);
  dash_gap = SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if(dash_sum < R_SMALL4)
    dash_sum = 0.5;

  I->shaderCGO = 0;
  I->N = 0;
  I->V = NULL;
  I->R.P = NULL;
  I->Obj = (CObject *) ds->Obj;
  I->ds = ds;

  n = 0;
  if(ds->NIndex) {
    I->V = VLAlloc(float, ds->NIndex * 10);
    CHECKOK(ok, I->V);
    for(a = 0; ok && a < ds->NIndex; a = a + 2) {
      v1 = ds->Coord + 3 * a;
      v2 = ds->Coord + 3 * (a + 1);

      /* vector from v2->v1 */
      subtract3f(v2, v1, d);

      l = (float) length3f(d);

      if(l > R_SMALL4) {

        copy3f(v1, d1);
	/* this makes d the direction vector of the distance measure from v2->v1 */
        normalize3f(d);

        if(dash_gap > R_SMALL4) {
          float avg[3], proj1[3], proj2[3];
          float l_left = l / 2.0F;
          float l_used = 0.0F;
          float half_dash_gap = dash_gap * 0.5;

          average3f(v1, v2, avg);
          while(ok && l_left > dash_sum) {
            VLACheck(I->V, float, (n * 3) + 11);
	    CHECKOK(ok, I->V);
            v = I->V + n * 3;
            scale3f(d, l_used + half_dash_gap, proj1);
            scale3f(d, l_used + dash_len + half_dash_gap, proj2);
            add3f(avg, proj1, v);
            add3f(avg, proj2, v + 3);
            subtract3f(avg, proj1, v + 6);
            subtract3f(avg, proj2, v + 9);
            n += 4;
            l_left -= dash_sum;
            l_used += dash_sum;
          }
          if(ok && l_left > dash_gap) {
            l_left -= dash_gap;
            scale3f(d, l_used + half_dash_gap, proj1);
            scale3f(d, l_used + l_left + half_dash_gap, proj2);
            VLACheck(I->V, float, (n * 3) + 11);
            v = I->V + n * 3;
            add3f(avg, proj1, v);
            add3f(avg, proj2, v + 3);
            subtract3f(avg, proj1, v + 6);
            subtract3f(avg, proj2, v + 9);
            n += 4;
          }
        } else if(dash_len > R_SMALL4) {
          VLACheck(I->V, float, (n * 3) + 5);
	  CHECKOK(ok, I->V);
	  if (ok){
	    v = I->V + n * 3;
	    copy3f(v1, v);
	    copy3f(v2, v + 3);
	    n += 2;
	  }
        }
      }
    }
    if (ok)
      VLASize(I->V, float, n * 3);
    CHECKOK(ok, I->V);
    if (ok)
      I->N = n;
  }
  if (!ok){
    RepDistDashFree(I);
    I = NULL;
  }
  return ((void *) (struct Rep *) I);
}
