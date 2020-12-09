
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

#include"Err.h"
#include"RepAngle.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"
#include"ShaderMgr.h"
#include"CGO.h"
#include"CoordSet.h"

struct RepAngle : Rep {
  using Rep::Rep;

  ~RepAngle() override;

  cRep_t type() const override { return cRepAngle; }
  void render(RenderInfo* info) override;

  pymol::vla<float> V;
  int N = 0;

  DistSet *ds;
  float linewidth, radius;
  CGO* shaderCGO = nullptr;
};

#include"ObjectDist.h"

RepAngle::~RepAngle()
{
  CGOFree(shaderCGO);
}

static int RepAngleCGOGenerate(RepAngle * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;
  float line_width;
  int ok = true;
  CGO *convertcgo = NULL;
  int dash_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_dash_as_cylinders);
  int color =
    SettingGet_color(G, NULL, I->ds->Obj->Setting.get(), cSetting_angle_color);
  I->linewidth = line_width = 
    SettingGet_f(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_width);
  I->radius =
    SettingGet_f(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_radius);
  line_width = SceneGetDynamicLineWidth(info, line_width);
  ok &= CGOSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
  if (ok)
    ok &= CGOResetNormal(I->shaderCGO, true);
  if (ok){
    if (color < 0) {
      color = I->getObj()->Color;
    }
    if(color >= 0){
      ok &= CGOColorv(I->shaderCGO, ColorGet(G, color));
    }
  }
  const float* v = I->V.data();
  int c = I->N;
  if (dash_as_cylinders){
    float axis[3];
    while(ok && c > 0) {
      const float* origin = v;
      v += 3;
      axis[0] = v[0] - origin[0];
      axis[1] = v[1] - origin[1];
      axis[2] = v[2] - origin[2];
      v += 3;
      ok &= (bool)I->shaderCGO->add<cgo::draw::shadercylinder>(origin, axis, 1.f, 15);
      c -= 2;
    }
  } else if (ok) {
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
  if (ok)
    ok &= CGOStop(I->shaderCGO);
  if (ok)
    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);    
  CHECKOK(ok, convertcgo);
  CGOFree(I->shaderCGO);    
  I->shaderCGO = convertcgo;
  convertcgo = NULL;
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
    } else {
      CGO *tmpCGO = CGONew(G);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DEFAULT_SHADER);
      if (ok) ok &= CGODisable(tmpCGO, CGO_GL_LIGHTING);
      convertcgo = CGOOptimizeToVBONotIndexedNoShader(I->shaderCGO);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
      if (ok) ok &= CGOAppendNoStop(tmpCGO, convertcgo);
      if (ok) ok &= CGODisable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
      if (ok) ok &= CGODisable(tmpCGO, GL_DEFAULT_SHADER);
      if (ok) ok &= CGOStop(tmpCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = tmpCGO;
    }
    convertcgo->use_shader = true;
  }
  CHECKOK(ok, convertcgo);
  if (convertcgo){
    CGOFree(I->shaderCGO);
    I->shaderCGO = convertcgo;
    convertcgo = NULL;
  }
  return ok;
}

static void RepAngleRenderImmediate(RepAngle * I, RenderInfo * info, int color, 
                                    short dash_transparency_enabled, float dash_transparency)
{
#ifndef PURE_OPENGL_ES_2
  PyMOLGlobals *G = I->G;
  const float* v = I->V.data();
  int c = I->N;
  float line_width;
  bool t_mode_3 =
    SettingGet_i(G, NULL, I->ds->Obj->Setting.get(), cSetting_transparency_mode) == 3;
  line_width = 
    SettingGet_f(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_width);
  line_width = SceneGetDynamicLineWidth(info, line_width);

  if(info->width_scale_flag) {
    glLineWidth(line_width * info->width_scale);
  } else {
    glLineWidth(line_width);
  }
  SceneResetNormal(G, true);
  
  if (dash_transparency_enabled){
    const float *col = ColorGet(G, color);
    glColor4f(col[0], col[1], col[2], 1.f-dash_transparency);
  } else {
    glColor3fv(ColorGet(G, color));
  }
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
#endif
}

void RepAngle::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  const float* v = I->V.data();
  int c = I->N;
  const float *vc;
  int round_ends;
  float line_width;
  int ok = true;
  float dash_transparency;
  short dash_transparency_enabled;
  int color =
    SettingGet_color(G, NULL, I->ds->Obj->Setting.get(), cSetting_angle_color);
  if(color < 0)
    color = getObj()->Color;
  I->linewidth = line_width = 
    SettingGet_f(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_width);
  I->radius =
    SettingGet_f(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_radius);
  round_ends =
    SettingGet_b(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_round_ends);
  line_width = SceneGetDynamicLineWidth(info, line_width);
  dash_transparency =
    SettingGet_f(G, NULL, I->ds->Obj->Setting.get(), cSetting_dash_transparency);
  dash_transparency = (dash_transparency < 0.f ? 0.f : (dash_transparency > 1.f ? 1.f : dash_transparency));
  dash_transparency_enabled = (dash_transparency > 0.f);

  if (!(ray || pick) && (info->pass == RenderPass::Antialias || (info->pass == RenderPass::Opaque) == dash_transparency_enabled))
    return;

  if(ray) {

    float radius;
    if (dash_transparency_enabled){
      ray->transparentf(dash_transparency);      
    }

    if(I->radius == 0.0F) {
      radius = ray->PixelRadius * line_width / 2.0F;
    } else {
      radius = I->radius;
    }

    vc = ColorGet(G, color);
    v = I->V.data();
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
      short use_shader, generate_shader_cgo = 0;

      use_shader = SettingGetGlobal_b(G, cSetting_dash_use_shader) & 
	           SettingGetGlobal_b(G, cSetting_use_shaders);

      if (!use_shader && I->shaderCGO){
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
	  ok &= RepAngleCGOGenerate(I, info);
   	} else {
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, I);
	  return;
	}
      }

      if (!generate_shader_cgo) {
	RepAngleRenderImmediate(I, info, color, dash_transparency_enabled, dash_transparency);
      } else {
	CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, I);
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->ds->Rep[cRepAngle] = NULL;
    delete I;
  }
}

Rep *RepAngleNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->G;
  int a;
  int n = 0;
  float *v, *v1, *v2, *v3, *v4, d1[3], d2[3], d3[3], n1[3], n3[3], l1, l2, x[3], y[3];
  float length, radius, angle, pos, phase;
  float dash_len, dash_gap, dash_sum;
  int ok = true;
  float dash_transparency;
  dash_transparency =
    SettingGet_f(G, NULL, ds->Obj->Setting.get(), cSetting_dash_transparency);
  dash_transparency = (dash_transparency < 0.f ? 0.f : (dash_transparency > 1.f ? 1.f : dash_transparency));

  PRINTFD(G, FB_RepAngle)
    "RepAngleNew: entered.\n" ENDFD;
  if(!ok || !ds->NAngleIndex) {
    return (NULL);
  }

  auto I = new RepAngle(ds->Obj, state);

  dash_len = SettingGet_f(G, NULL, ds->Obj->Setting.get(), cSetting_dash_length);
  dash_gap = SettingGet_f(G, NULL, ds->Obj->Setting.get(), cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if(dash_sum < R_SMALL4)
    dash_sum = 0.1F;

  I->ds = ds;

  n = 0;
  if(ds->NAngleIndex) {
    I->V.resize(ds->NAngleIndex * 10);
    CHECKOK(ok, I->V);
    for(a = 0; ok && a < ds->NAngleIndex; a = a + 5) {
      v1 = ds->AngleCoord + 3 * a;
      v2 = ds->AngleCoord + 3 * (a + 1);
      v3 = ds->AngleCoord + 3 * (a + 2);
      v4 = ds->AngleCoord + 3 * (a + 3);
      subtract3f(v1, v2, d1);
      subtract3f(v3, v2, d2);

      l1 = (float) length3f(d1);
      l2 = (float) length3f(d2);

      if(l1 > l2)
        radius = l2;
      else
        radius = l1;
      radius *= SettingGet_f(G, NULL, ds->Obj->Setting.get(), cSetting_angle_size);

      angle = get_angle3f(d1, d2);

      normalize23f(d1, n1);

      remove_component3f(d2, n1, d3);

      if(length3f(d3) < R_SMALL8) {
        d3[0] = 1.0F;
        d3[1] = 0.0F;
        d3[2] = 0.0F;
      } else {
        normalize23f(d3, n3);
      }

      scale3f(n1, radius, x);
      scale3f(n3, radius, y);

      if(v4[0] != 0.0F) {       /* line 1 flag */
        I->V.check((n * 3) + 5);
	CHECKOK(ok, I->V);
	if (ok){
	  v = I->V + n * 3;
	  copy3f(v1, v);
	  v += 3;
	  copy3f(v2, v);
	  n += 2;
	}
      }
      
      if(ok && v4[1] != 0.0F) {       /* line 2 flag */
	I->V.check((n * 3) + 5);
	CHECKOK(ok, I->V);
	if (ok){
	  v = I->V + n * 3;
	  copy3f(v3, v);
	  v += 3;
	  copy3f(v2, v);
	  n += 2;
	}
      }
      if (!ok)
	break;
      /* now we have a relevant orthogonal axes */

      length = (float) (angle * radius * 2);

      /* figure out dash/gap phasing that will lead to nicely spaced dashes and gaps */

      phase = dash_sum - (float) fmod(length / 2 + (dash_gap / 2), dash_sum);
      pos = -phase;

      if(length > R_SMALL4) {

        float vx[3], vy[3];
        float cur_angle;
        float cons_pos1, cons_pos2;

        while(ok && pos < length) {

          I->V.check((n * 3) + 5);
	  CHECKOK(ok, I->V);

	  if (!ok)
	    break;
          cons_pos1 = pos;
          if(cons_pos1 < 0.0F)
            cons_pos1 = 0.0F;
          cons_pos2 = pos + dash_len;
          if(cons_pos2 > length)
            cons_pos2 = length;

          if(cons_pos1 < cons_pos2) {
            cur_angle = angle * cons_pos1 / length;

            v = I->V + n * 3;
            scale3f(x, (float) cos(cur_angle), vx);
            scale3f(y, (float) sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(v2, v, v);

            cur_angle = angle * cons_pos2 / length;

            v += 3;
            scale3f(x, (float) cos(cur_angle), vx);
            scale3f(y, (float) sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(v2, v, v);

            n += 2;
          }
          pos += dash_sum;
        }
      }
    }
    if (ok)
      I->V.resize(n * 3);
    CHECKOK(ok, I->V);
    if (ok)
      I->N = n;
  }
  if (!ok){
    delete I;
    I = NULL;
  }
  return (Rep *) I;
}
