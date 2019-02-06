
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
#ifdef _WEBGL
#include "WebPyMOLLibrary.h"
#endif

typedef struct RepDistDash {
  Rep R;
  float *V;
  int N;
  CObject *Obj;
  DistSet *ds;
  float linewidth, radius;
  CGO *shaderCGO;
  bool shaderCGO_has_cylinders, shaderCGO_has_trilines;
} RepDistDash;

#include"ObjectDist.h"

static
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

/* Has no prototype */
static void RepDistDashCGOGenerate(RepDistDash * I)
{
  int ok = true;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  int color =
    SettingGet_color(G, NULL, I->ds->Obj->Obj.Setting, cSetting_dash_color);
  short dash_as_cylinders = 0;

  dash_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_dash_as_cylinders);

  if (ok)
    ok &= CGOSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
  if (ok)
    ok &= CGOResetNormal(I->shaderCGO, true);
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
      ok &= (bool)I->shaderCGO->add<cgo::draw::shadercylinder>(origin, axis, 1.f, 15);
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
}

static void RepDistDashRender(RepDistDash * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  auto pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  const float *vc;
  int round_ends;
  int ok = true;
  int color =
    SettingGet_color(G, NULL, I->ds->Obj->Obj.Setting, cSetting_dash_color);
  float line_width =
    SettingGet_f(G, NULL, I->ds->Obj->Obj.Setting, cSetting_dash_width);
  float dash_transparency =
    SettingGet_f(G, NULL, I->ds->Obj->Obj.Setting, cSetting_dash_transparency);
  bool t_mode_3 =
    SettingGet_i(G, NULL, I->ds->Obj->Obj.Setting, cSetting_transparency_mode) == 3;
  short dash_transparency_enabled;
  dash_transparency = (dash_transparency < 0.f ? 0.f : (dash_transparency > 1.f ? 1.f : dash_transparency));
  dash_transparency_enabled = (dash_transparency > 0.f);

  if (!(ray || pick) && (!info->pass || (info->pass > 0) == dash_transparency_enabled))
    return;

  if(color < 0)
    color = I->Obj->Color;

  I->radius =
    SettingGet_f(G, NULL, I->ds->Obj->Obj.Setting, cSetting_dash_radius);
  round_ends =
    SettingGet_b(G, NULL, I->ds->Obj->Obj.Setting, cSetting_dash_round_ends);
  line_width = SceneGetDynamicLineWidth(info, line_width);

  if(ray) {
    float radius;
    if (dash_transparency_enabled){
      ray->transparentf(dash_transparency);      
    }
    if(I->radius <= 0.0F) {
      radius = ray->PixelRadius * line_width / 2.0F;
    } else {
      radius = I->radius;
    }

    vc = ColorGet(G, color);
    v = I->V;
    c = I->N;

    while(ok && c > 0) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]); */
      if(round_ends) {
        ok &= ray->sausage3fv(v, v + 3, radius, vc, vc);
      } else {
        ok &= ray->customCylinder3fv(v, v + 3, radius, vc, vc, cCylCapFlat, cCylCapFlat);
      }
      v += 6;
      c -= 2;
    }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      short use_shader, generate_shader_cgo = 0, dash_as_cylinders = 0;

      use_shader = SettingGetGlobal_b(G, cSetting_dash_use_shader) & 
                   SettingGetGlobal_b(G, cSetting_use_shaders);
      dash_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_dash_as_cylinders);
      if (!GET_FRAGDEPTH_SUPPORT() && dash_as_cylinders)
        dash_as_cylinders = false;
      if (!use_shader && I->shaderCGO){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }
      if (I->shaderCGO && (dash_as_cylinders ^ I->shaderCGO_has_cylinders)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }
      if (I->shaderCGO && !dash_as_cylinders && I->shaderCGO_has_trilines != SettingGetGlobal_b(G, cSetting_trilines)){
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
	  if (dash_transparency_enabled){
	    CGOAlpha(I->shaderCGO, 1.f-dash_transparency);
	  }
	  RepDistDashCGOGenerate(I);
	} else if (ok) {
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  return;
	}
      }
#ifndef PURE_OPENGL_ES_2
      if (!generate_shader_cgo) {
	if(info->width_scale_flag) {
	  glLineWidth(line_width * info->width_scale);
	} else {
	  glLineWidth(line_width);
	}
        SceneResetNormal(G, true);
      }

      {
	if(color >= 0){
	  if (dash_transparency_enabled){
	    const float *col = ColorGet(G, color);
	    glColor4f(col[0], col[1], col[2], 1.f-dash_transparency);
	  } else {
	    glColor3fv(ColorGet(G, color));
	  }
	} else if (dash_transparency_enabled){
	  float col[4];
	  copy3f(ColorGet(I->Obj->G, I->Obj->Color), col);
	  col[3] = 1.f-dash_transparency;
	  glColor4fv(col);
	}
	v = I->V;
	c = I->N;
	if (dash_transparency_enabled && !t_mode_3)
	  glDisable(GL_DEPTH_TEST);	
	if(!info->line_lighting)
	  glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	while(c > 0) {
	  glVertex3fv(v);
	  v += 3;
	  glVertex3fv(v);
	  v += 3;
	  c -= 2;
	}
	glEnd();
	glEnable(GL_LIGHTING);
        if (dash_transparency_enabled && !t_mode_3)
          glEnable(GL_DEPTH_TEST);	
      }
#endif
      if (use_shader) {
	if (generate_shader_cgo){
	  CGO *convertcgo = NULL;
	  if (ok)
	    ok &= CGOStop(I->shaderCGO);
	  {
            bool trilines = SettingGetGlobal_b(G, cSetting_trilines);
            if (dash_as_cylinders || !trilines) {
	  if (ok)
	    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);    
	  CHECKOK(ok, convertcgo);
	  CGOFree(I->shaderCGO);    
	  I->shaderCGO = convertcgo;
	  convertcgo = NULL;
            }
	  if (ok){
	    if (dash_as_cylinders){
		CGO *tmpCGO = CGONew(G);
		if (ok) ok &= CGOEnable(tmpCGO, GL_CYLINDER_SHADER);
		if (ok) ok &= CGOSpecial(tmpCGO, CYLINDER_WIDTH_FOR_DISTANCES);
		convertcgo = CGOConvertShaderCylindersToCylinderShader(I->shaderCGO, tmpCGO);
		if (ok) ok &= CGOEnable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
		if (ok) ok &= CGOAppendNoStop(tmpCGO, convertcgo);
		if (ok) ok &= CGODisable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
		if (ok) ok &= CGODisable(tmpCGO, GL_CYLINDER_SHADER);
		if (ok) ok &= CGOStop(tmpCGO);
		CGOFreeWithoutVBOs(convertcgo);
		convertcgo = tmpCGO;
                I->shaderCGO_has_cylinders = true;
                I->shaderCGO_has_trilines = false;
	      } else {
		CGO *tmpCGO = CGONew(G);
                int shader = trilines ? GL_TRILINES_SHADER : GL_DEFAULT_SHADER;
                if (ok) ok &= CGOEnable(tmpCGO, shader);
		if (ok) ok &= CGODisable(tmpCGO, CGO_GL_LIGHTING);
		if (trilines) {
                  if (ok) ok &= CGOSpecial(tmpCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
                  convertcgo = CGOConvertLinesToTrilines(I->shaderCGO, false);
		} else {
		  convertcgo = CGOOptimizeToVBONotIndexedNoShader(I->shaderCGO, 0);
		}
                I->shaderCGO_has_trilines = trilines;
		if (ok) ok &= CGOEnable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
		if (ok) ok &= CGOAppendNoStop(tmpCGO, convertcgo);
		if (ok) ok &= CGODisable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
		if (ok) ok &= CGODisable(tmpCGO, shader);
		if (ok) ok &= CGOStop(tmpCGO);
		CGOFreeWithoutVBOs(convertcgo);
		convertcgo = tmpCGO;
                I->shaderCGO_has_cylinders = false;
	      }
	      convertcgo->use_shader = true;
	    }
	  }
	  if (convertcgo){
	    CGOFree(I->shaderCGO);
	    I->shaderCGO = convertcgo;
	    convertcgo = NULL;
	  }
	}
	
	if (ok) {
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	}
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->ds->Rep[cRepDash] = NULL;
    RepDistDashFree(I);
  }
}

Rep *RepDistDashNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->State.G;
  int a;
  int n;
  float *v, *v1, *v2, d[3];
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
  I->R.obj = &ds->Obj->Obj;
  I->R.context.state = state;
  dash_len = SettingGet_f(G, NULL, ds->Obj->Obj.Setting, cSetting_dash_length);
  dash_gap = SettingGet_f(G, NULL, ds->Obj->Obj.Setting, cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if(dash_sum < R_SMALL4)
    dash_sum = 0.5;

  I->shaderCGO = 0;
  I->shaderCGO_has_cylinders = false;
  I->shaderCGO_has_trilines = false;
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
  return (Rep *) I;
}
