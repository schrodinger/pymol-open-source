
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
#include"RepNonbonded.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Setting.h"
#include"ShaderMgr.h"
#include"CGO.h"

struct RepNonbonded : Rep {
  using Rep::Rep;

  ~RepNonbonded() override;

  cRep_t type() const override { return cRepNonbonded; }
  void render(RenderInfo* info) override;

  CGO *primitiveCGO;
  CGO *shaderCGO;
  bool shaderCGO_has_cylinders;
};

#include"ObjectMolecule.h"

RepNonbonded::~RepNonbonded()
{
  CGOFree(primitiveCGO);
  CGOFree(shaderCGO);
}

void RepNonbondedRenderImmediate(CoordSet * cs, RenderInfo * info)
{
#ifndef PURE_OPENGL_ES_2
  PyMOLGlobals *G = cs->G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    int active = false;
    ObjectMolecule *obj = cs->Obj;
    float line_width =
      SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_line_width);
    float nonbonded_size =
      SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_nonbonded_size);

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
      int nIndex = cs->NIndex;
      const AtomInfoType* atomInfo = obj->AtomInfo.data();
      const int* i2a = cs->IdxToAtm.data();
      int last_color = -1;
      const float *v = cs->Coord;

      for(a = 0; a < nIndex; a++) {
        const AtomInfoType* ai = atomInfo + *(i2a++);
        if((!ai->bonded) && (ai->visRep & cRepNonbondedBit)) {
          int c = ai->color;
          float v0 = v[0];
          float v1 = v[1];
          float v2 = v[2];
          active = true;
          if(c != last_color) {
            last_color = c;
            glColor3fv(ColorGet(G, c));
          }

          glVertex3f(v0 - nonbonded_size, v1, v2);
          glVertex3f(v0 + nonbonded_size, v1, v2);

          glVertex3f(v0, v1 - nonbonded_size, v2);
          glVertex3f(v0, v1 + nonbonded_size, v2);

          glVertex3f(v0, v1, v2 - nonbonded_size);
          glVertex3f(v0, v1, v2 + nonbonded_size);
        }
        v += 3;
      }
    }
    glEnd();
    glEnable(GL_LIGHTING);
    if(!active)
      cs->Active[cRepNonbonded] = false;
  }
#endif
}

static int RepNonbondedCGOGenerate(RepNonbonded * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;
  float alpha;
  int ok = true;
  CGO *convertcgo = NULL;
  short nonbonded_as_cylinders ;
  short use_shader;
  float nonbonded_size =
    SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_nonbonded_size);

  nonbonded_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_nonbonded_as_cylinders);
  use_shader = SettingGetGlobal_b(G, cSetting_nonbonded_use_shader) & 
    SettingGetGlobal_b(G, cSetting_use_shaders);

  alpha =
    SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_nonbonded_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  
  if (use_shader){
    if (ok && I->shaderCGO){
      CGOFree(I->shaderCGO);
    }
    if (ok){
      if (nonbonded_as_cylinders){
        CGO *tmpCGO = CGONew(G);
        if (ok) ok &= CGOEnable(tmpCGO, GL_CYLINDER_SHADER);
        if (ok) ok &= CGOSpecial(tmpCGO, CYLINDER_WIDTH_FOR_REPWIRE);
        convertcgo = CGOConvertCrossesToCylinderShader(I->primitiveCGO, tmpCGO, nonbonded_size);
        if (ok) ok &= CGOAppendNoStop(tmpCGO, convertcgo);
        if (ok) ok &= CGODisable(tmpCGO, GL_CYLINDER_SHADER);
        if (ok) ok &= CGOStop(tmpCGO);
        CGOFreeWithoutVBOs(convertcgo);
        I->shaderCGO_has_cylinders = true;
        convertcgo = tmpCGO;
      } else {
        bool trilines = SettingGetGlobal_b(G, cSetting_trilines);
        CGO *tmpCGO = CGONew(G), *tmp2CGO;
        int shader = trilines ? GL_TRILINES_SHADER : GL_LINE_SHADER;

        if (ok) ok &= CGOEnable(tmpCGO, shader);
        if (ok) ok &= CGODisable(tmpCGO, CGO_GL_LIGHTING);
        if (trilines) {
          if (ok) ok &= CGOSpecial(tmpCGO, LINEWIDTH_DYNAMIC_WITH_SCALE);
          tmp2CGO = CGOConvertCrossesToTrilinesShader(I->primitiveCGO, tmpCGO, nonbonded_size);
        } else {
          tmp2CGO = CGOConvertCrossesToLinesShader(I->primitiveCGO, tmpCGO, nonbonded_size);
        }
        if (ok) ok &= CGOAppendNoStop(tmpCGO, tmp2CGO);
        if (ok) ok &= CGODisable(tmpCGO, shader);
        if (ok) ok &= CGOStop(tmpCGO);
        CGOFreeWithoutVBOs(tmp2CGO);
        convertcgo = tmpCGO;
        I->shaderCGO_has_cylinders = false;
      }
      convertcgo->use_shader = true;
    }
    CHECKOK(ok, convertcgo);
    if (!ok || convertcgo){
      CGOFree(I->shaderCGO);
      I->shaderCGO = convertcgo;
      I->shaderCGO->use_shader = use_shader;
      convertcgo = NULL;
    }
  } else {
    // no shaders
    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);
    CGOFree(I->shaderCGO);
    I->shaderCGO = convertcgo;
    I->shaderCGO->use_shader = use_shader;
    convertcgo = NULL;
  }
  return ok;
}

void RepNonbonded::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  int ok = true;
  float alpha =
    SettingGet_f(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_nonbonded_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  if(ray) {
#ifndef _PYMOL_NO_RAY
    CGORenderRay(I->primitiveCGO, ray, info, NULL, NULL, I->cs->Setting.get(), I->cs->Obj->Setting.get());
    ray->transparentf(0.0);
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      CGORenderGLPicking(I->shaderCGO ? I->shaderCGO : I->primitiveCGO, info, &I->context, I->cs->Setting.get(), I->obj->Setting.get());
    } else {
      /* not pick, but render */
      bool use_shader = SettingGetGlobal_b(G, cSetting_nonbonded_use_shader) && SettingGetGlobal_b(G, cSetting_use_shaders);
      if (!use_shader){
        CGORenderGL(I->primitiveCGO, NULL, NULL, NULL, info, I);
        return;
      }
      bool nonbonded_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_nonbonded_as_cylinders);
      if (I->shaderCGO && use_shader != I->shaderCGO->use_shader){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (I->shaderCGO && (nonbonded_as_cylinders ^ I->shaderCGO_has_cylinders)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (!I->shaderCGO){
        I->shaderCGO = CGONew(G);
        CHECKOK(ok, I->shaderCGO);
        if (ok){
          I->shaderCGO->use_shader = use_shader;
        }
        ok &= RepNonbondedCGOGenerate(I, info);
      }
      CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, I);
    }
  }
}

Rep *RepNonbondedNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  bool hasNonbondedAtoms = false;

  ObjectMolecule *obj = cs->Obj;

  if((obj->RepVisCache & cRepNonbondedBit)){
    for(int a = 0; a < cs->NIndex; a++) {
      AtomInfoType *ai = obj->AtomInfo + cs->IdxToAtm[a];
      if (!ai->bonded && (ai->visRep & cRepNonbondedBit)){
        hasNonbondedAtoms = true;
        break;
      }
    }
  }
  if(!hasNonbondedAtoms) {
    return (NULL);              /* skip if no dots are visible */
  }

  auto I = new RepNonbonded(cs, state);

  I->shaderCGO = NULL;

  I->primitiveCGO = CGONew(G);

  CGOSpecialWithArg(I->primitiveCGO, LINE_LIGHTING, 0.f);
  CGOSpecial(I->primitiveCGO, LINEWIDTH_FOR_LINES);
  CGOBegin(I->primitiveCGO, GL_LINES); // for immediate mode
  bool first = true;
  int a1, c1;
  float tmpColor[3];
  for(int a = 0; a < cs->NIndex; a++){
    a1 = cs->IdxToAtm[a];
    AtomInfoType *ai = obj->AtomInfo + a1;
    if(!ai->bonded && (ai->visRep & cRepNonbondedBit)) {
      c1 = ai->color;
      const float* v1 = cs->coordPtr(a);
      ColorGetCheckRamped(G, c1, v1, tmpColor, state);
      if (first || !equal3f(I->primitiveCGO->color, tmpColor)){
        CGOColorv(I->primitiveCGO, tmpColor);
      }
      CGOPickColor(I->primitiveCGO, a1, (ai->masked) ? cPickableNoPick : cPickableAtom);
      CGOVertexCrossv(I->primitiveCGO, v1);
      first = false;
    }
  }
  CGOEnd(I->primitiveCGO); // for immediate mode
  CGOSpecialWithArg(I->primitiveCGO, LINE_LIGHTING, 1.f);
  return (Rep *) I;
}
