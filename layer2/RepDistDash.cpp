
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

#include <optional>

struct RepDistDash : Rep {
  using Rep::Rep;

  ~RepDistDash() override;

  cRep_t type() const override { return cRepDash; }
  void render(RenderInfo* info) override;

  float* V = nullptr;
  int N = 0;
  DistSet *ds;
  float linewidth, radius;
  CGO* shaderCGO = nullptr;
  bool shaderCGO_has_cylinders = false;
  bool shaderCGO_has_trilines = false;
};

#include"ObjectDist.h"

RepDistDash::~RepDistDash()
{
  CGOFree(shaderCGO);
  VLAFreeP(V);
}

/* Has no prototype */
static void RepDistDashCGOGenerate(RepDistDash* I, std::optional<float> dash_transparency)
{
  int ok = true;
  PyMOLGlobals* G = I->G;
  float* v = I->V;
  int c = I->N;
  int color = SettingGet_color(
      G, nullptr, I->ds->Obj->Setting.get(), cSetting_dash_color);

  bool dash_as_cylinders = SettingGet<bool>(G, cSetting_render_as_cylinders) &&
                           SettingGet<bool>(G, cSetting_dash_as_cylinders);

  I->shaderCGO = CGONew(G);
  CHECKOK(ok, I->shaderCGO);
  if (ok)
    I->shaderCGO->use_shader = true;

  if (dash_transparency) {
    CGOAlpha(I->shaderCGO, 1.f - *dash_transparency);
  }

  if (ok)
    ok &= CGOSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
  if (ok)
    ok &= CGOResetNormal(I->shaderCGO, true);
  if (ok) {
    if (color < 0) {
      color = I->getObj()->Color;
    }
    if (color >= 0) {
      ok &= CGOColorv(I->shaderCGO, ColorGet(G, color));
    }
  }
  v = I->V;
  c = I->N;
  if (dash_as_cylinders) {
    float *origin = nullptr, axis[3];
    while (ok && c > 0) {
      origin = v;
      v += 3;
      axis[0] = v[0] - origin[0];
      axis[1] = v[1] - origin[1];
      axis[2] = v[2] - origin[2];
      v += 3;
      ok &= (bool) I->shaderCGO->add<cgo::draw::shadercylinder>(
          origin, axis, 1.f, 15);
      c -= 2;
    }
  } else {
    ok &= CGOBegin(I->shaderCGO, GL_LINES);
    while (ok && c > 0) {
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

  CGO* convertcgo = nullptr;
  if (ok)
    ok &= CGOStop(I->shaderCGO);
  bool trilines = SettingGetGlobal_b(G, cSetting_trilines);
  if (dash_as_cylinders || !trilines) {
    if (ok)
      convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);
    CHECKOK(ok, convertcgo);
    CGOFree(I->shaderCGO);
    I->shaderCGO = convertcgo;
    convertcgo = nullptr;
  }
  if (ok) {
    if (dash_as_cylinders) {
      CGO* tmpCGO = CGONew(G);
      if (ok)
        ok &= CGOEnable(tmpCGO, GL_CYLINDER_SHADER);
      if (ok)
        ok &= CGOSpecial(tmpCGO, CYLINDER_WIDTH_FOR_DISTANCES);
      convertcgo =
          CGOConvertShaderCylindersToCylinderShader(I->shaderCGO, tmpCGO);
      if (ok)
        ok &= CGOEnable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
      if (ok)
        ok &= CGOAppendNoStop(tmpCGO, convertcgo);
      if (ok)
        ok &= CGODisable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
      if (ok)
        ok &= CGODisable(tmpCGO, GL_CYLINDER_SHADER);
      if (ok)
        ok &= CGOStop(tmpCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = tmpCGO;
      I->shaderCGO_has_cylinders = true;
      I->shaderCGO_has_trilines = false;
    } else {
      CGO* tmpCGO = CGONew(G);
      int shader = trilines ? GL_TRILINES_SHADER : GL_DEFAULT_SHADER;
      if (ok)
        ok &= CGOEnable(tmpCGO, shader);
      if (ok)
        ok &= CGODisable(tmpCGO, CGO_GL_LIGHTING);
      if (trilines) {
        if (ok)
          ok &= CGOSpecial(tmpCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
        convertcgo = CGOConvertLinesToTrilines(I->shaderCGO, false);
      } else {
        convertcgo = CGOOptimizeToVBONotIndexedNoShader(I->shaderCGO);
      }
      I->shaderCGO_has_trilines = trilines;
      if (ok)
        ok &= CGOEnable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
      if (ok)
        ok &= CGOAppendNoStop(tmpCGO, convertcgo);
      if (ok)
        ok &= CGODisable(tmpCGO, GL_DASH_TRANSPARENCY_DEPTH_TEST);
      if (ok)
        ok &= CGODisable(tmpCGO, shader);
      if (ok)
        ok &= CGOStop(tmpCGO);
      CGOFreeWithoutVBOs(convertcgo);
      convertcgo = tmpCGO;
      I->shaderCGO_has_cylinders = false;
    }
    convertcgo->use_shader = true;
  }
  if (convertcgo) {
    CGOFree(I->shaderCGO);
    I->shaderCGO = convertcgo;
    convertcgo = nullptr;
  }
}

static void RepDistDashRay(RepDistDash* I, RenderInfo* info,
    std::optional<float> dash_transparency, float line_width, int color)
{
  auto G = I->G;
  bool round_ends = SettingGet<bool>(
      G, nullptr, I->ds->Obj->Setting.get(), cSetting_dash_round_ends);
  auto ray = info->ray;
  float radius;
  if (dash_transparency){
    ray->transparentf(*dash_transparency);      
  }
  if(I->radius <= 0.0F) {
    radius = ray->PixelRadius * line_width / 2.0F;
  } else {
    radius = I->radius;
  }

  auto vc = ColorGet(G, color);
  auto v = I->V;
  auto c = I->N;
  bool ok = true;

  while (ok && c > 0) {
    /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]); */
    if (round_ends) {
      ok &= ray->sausage3fv(v, v + 3, radius, vc, vc);
    } else {
      ok &= ray->customCylinder3fv(v, v + 3, radius, vc, vc, cCylCapFlat, cCylCapFlat);
    }
    v += 6;
    c -= 2;
  }
}

static void RepDistDashRenderRaster(
    RepDistDash* I, RenderInfo* info, std::optional<float> dash_transparency)
{
  auto G = I->G;
  bool dash_as_cylinders = SettingGet<bool>(G, cSetting_render_as_cylinders) &&
                           SettingGet<bool>(G, cSetting_dash_as_cylinders);

  if (!GET_FRAGDEPTH_SUPPORT() && dash_as_cylinders)
    dash_as_cylinders = false;
  if (I->shaderCGO && (dash_as_cylinders ^ I->shaderCGO_has_cylinders)) {
    CGOFree(I->shaderCGO);
    I->shaderCGO = nullptr;
  }
  if (I->shaderCGO && !dash_as_cylinders &&
      I->shaderCGO_has_trilines != SettingGetGlobal_b(G, cSetting_trilines)) {
    CGOFree(I->shaderCGO);
    I->shaderCGO = nullptr;
  }

  if (!I->shaderCGO) {
    RepDistDashCGOGenerate(I, dash_transparency);
  }

  CGORender(I->shaderCGO, nullptr, nullptr, nullptr, info, I);
}

static void RepDistDashRenderPick(RepDistDash* I, RenderInfo* info)
{
  /* NO-OP */
}

static void RepDistDashRenderImmediate(RepDistDash* I, RenderInfo* info,
    std::optional<float> dash_transparency, float line_width, int color)
{
  auto G = I->G;
  bool t_mode_3 = SettingGet<int>(G, nullptr, I->ds->Obj->Setting.get(),
                      cSetting_transparency_mode) == 3;

  if (I->shaderCGO) {
    CGOFree(I->shaderCGO);
    I->shaderCGO = nullptr;
  }

#ifndef PURE_OPENGL_ES_2
  if (info->width_scale_flag) {
    glLineWidth(line_width * info->width_scale);
  } else {
    glLineWidth(line_width);
  }
  SceneResetNormal(G, true);

  if (color >= 0) {
    if (dash_transparency) {
      const float* col = ColorGet(G, color);
      glColor4f(col[0], col[1], col[2], 1.f - *dash_transparency);
    } else {
      glColor3fv(ColorGet(G, color));
    }
  }
  auto v = I->V;
  auto c = I->N;
  if (dash_transparency && !t_mode_3)
    glDisable(GL_DEPTH_TEST);
  if (!info->line_lighting)
    glDisable(GL_LIGHTING);
  glBegin(GL_LINES);
  while (c > 0) {
    glVertex3fv(v);
    v += 3;
    glVertex3fv(v);
    v += 3;
    c -= 2;
  }
  glEnd();
  glEnable(GL_LIGHTING);
  if (dash_transparency && !t_mode_3)
    glEnable(GL_DEPTH_TEST);
#endif
}

void RepDistDash::render(RenderInfo* info)
{
  auto I = this;
  CRay* ray = info->ray;
  auto pick = info->pick;
  int color =
    SettingGet_color(G, nullptr, I->ds->Obj->Setting.get(), cSetting_dash_color);


  float dash_transparency =
    SettingGet<float>(G, nullptr, I->ds->Obj->Setting.get(), cSetting_dash_transparency);
  dash_transparency = std::clamp(dash_transparency, 0.0f, 1.0f);
  std::optional<float> dash_transparency_opt;
  if (dash_transparency > 0.f) {
    dash_transparency_opt = dash_transparency;
  }

  if (!(ray || pick) &&
      (info->pass == RenderPass::Antialias ||
          (info->pass == RenderPass::Opaque) == dash_transparency_opt.has_value()))
    return;

  if (color < 0)
    color = getObj()->Color;

  I->radius = SettingGet<float>(
      G, nullptr, I->ds->Obj->Setting.get(), cSetting_dash_radius);

  float line_width =
    SettingGet<float>(G, nullptr, I->ds->Obj->Setting.get(), cSetting_dash_width);
  line_width = SceneGetDynamicLineWidth(info, line_width);

  if (ray) {
    RepDistDashRay(I, info, dash_transparency_opt, line_width, color);
    return;
  }
  if (!(G->HaveGUI && G->ValidContext)) {
    return;
  }
  if (pick) {
    RepDistDashRenderPick(I, info);
    return;
  }

  bool use_shader = SettingGet<bool>(G, cSetting_dash_use_shader) &
                    SettingGet<bool>(G, cSetting_use_shaders);

  if (use_shader) {
    RepDistDashRenderRaster(I, info, dash_transparency_opt);
    return;
  }
  RepDistDashRenderImmediate(I, info, dash_transparency_opt, line_width, color);
}

Rep *RepDistDashNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->G;
  int a;
  int n;
  float *v, *v1, *v2, d[3];
  float l;
  float dash_len, dash_gap, dash_sum;
  int ok = true;

  if(!ok || !ds->NIndex) {
    return (nullptr);
  }

  auto I = new RepDistDash(ds->Obj, state);

  dash_len = SettingGet_f(G, nullptr, ds->Obj->Setting.get(), cSetting_dash_length);
  dash_gap = SettingGet_f(G, nullptr, ds->Obj->Setting.get(), cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if(dash_sum < R_SMALL4)
    dash_sum = 0.5;

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
    delete I;
    I = nullptr;
  }
  return (Rep *) I;
}
