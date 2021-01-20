
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

struct RepWireBond : Rep {
  using Rep::Rep;

  ~RepWireBond() override;

  cRep_t type() const override { return cRepLine; }
  void render(RenderInfo* info) override;

  CGO *shaderCGO = nullptr;
  CGO *primitiveCGO = nullptr;
  bool shaderCGO_has_cylinders = false;
};

#include"ObjectMolecule.h"

static int RepLine(CGO *cgo, bool s1, bool s2, bool isRamped,
    const float *v1,
    const float *v2,
    const float *v1color,
    unsigned int b1, unsigned int b2, int a,
    const float *v2color,
    bool b1masked, bool b2masked)
{
  int ok = true;
  if (s1 && s2){
    CGOColorv(cgo, v1color);
    CGOPickColor(cgo, b1, b1masked ? cPickableNoPick : a);
    {
      // if not ramped, then insert vertices so colors are not interpolated since lines interpolate by default
      bool eq = equal3f(v1color, v2color);
      bool split = !eq || b1 != b2;
      if (split){
        cgo->add<cgo::draw::splitline>(v1, v2, v2color, b2, b2masked ? cPickableNoPick : a, isRamped, b1==b2, eq);
        cgo->current_pick_color_index = b2;
        cgo->current_pick_color_bond = b2masked ? cPickableNoPick : a;
      } else {
        cgo->add<cgo::draw::line>(v1, v2);
      }
    }
  } else {
    // if half bond, then split for either s1 or s2
    float h[3];
    average3f(v1, v2, h);
    if (s1){
      CGOColorv(cgo, v1color);
      CGOPickColor(cgo, b1, b1masked ? cPickableNoPick : a);
      cgo->add<cgo::draw::line>(v1, h);
    } else {
      if (v2color)
        CGOColorv(cgo, v2color);
      if (b2)
        CGOPickColor(cgo, b2, b2masked ? cPickableNoPick : a);
      cgo->add<cgo::draw::line>(h, v2);
    }
  }
  return ok;
}

static void RepValence(CGO *cgo, bool s1, bool s2, bool isRamped,
    const float *v1,
    const float *v2,
    const int *other,
    int a1, int a2,
    const float *coord,
    const float *color1,
    const float *color2, int ord,
    float tube_size, int fancy, unsigned int b1, unsigned int b2, int a, bool b1masked, bool b2masked)
{

  float d[3], t[3], p0[3], p1[3], p2[3];
  int a3;
  const float indent = tube_size;
  /* direction vector */
  subtract3f(v2, v1, p0);
  copy3f(p0, d);
  normalize3f(p0);
  /* need a prioritized third atom to get planarity */
  a3 = ObjectMoleculeGetPrioritizedOther(other, a1, a2, NULL);
  if(a3 < 0) {
    t[0] = p0[0];
    t[1] = p0[1];
    t[2] = -p0[2];
  } else {
    subtract3f(coord + 3 * a3, v1, t);
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
  mult3f(p2, tube_size, t);

  bool ord3 = (ord == 3);
  if (ord3)
    mult3f(t, 2.f, t);
  if (fancy || ord3){
    RepLine(cgo, s1, s2, isRamped, v1, v2, color1, b1, b2, a, color2, b1masked, b2masked);
  }
  if(fancy) {
    float f[] = { indent, 1.f - indent };
    float f_1[] = { 1.f - f[0], 1.f - f[1] };
    float vv1[] = { (f_1[0] * v1[0] + f[0] * v2[0]) - 2 * t[0],
                    (f_1[0] * v1[1] + f[0] * v2[1]) - 2 * t[1],
                    (f_1[0] * v1[2] + f[0] * v2[2]) - 2 * t[2] };
    float vv2[] = { (f_1[1] * v1[0] + f[1] * v2[0]) - 2 * t[0],
                    (f_1[1] * v1[1] + f[1] * v2[1]) - 2 * t[1],
                    (f_1[1] * v1[2] + f[1] * v2[2]) - 2 * t[2] };
    RepLine(cgo, s1, s2, isRamped, vv1, vv2, color1, b1, b2, a, color2, b1masked, b2masked);
  } else {
    float vv1[][3] = { { v1[0] - t[0], v1[1] - t[1], v1[2] - t[2] },
                       { v1[0] + t[0], v1[1] + t[1], v1[2] + t[2] } };
    float vv2[][3] = { { v2[0] - t[0], v2[1] - t[1], v2[2] - t[2] },
                       { v2[0] + t[0], v2[1] + t[1], v2[2] + t[2] } };
    RepLine(cgo, s1, s2, isRamped, vv1[0], vv2[0], color1, b1, b2, a, color2, b1masked, b2masked);
    RepLine(cgo, s1, s2, isRamped, vv1[1], vv2[1], color1, b1, b2, a, color2, b1masked, b2masked);
  }
}

static void RepAromatic(CGO *cgo, bool s1, bool s2, bool isRamped,
    const float *v1,
    const float *v2,
    const int *other,
    int a1, int a2,
    const float *coord,
    const float *color1,
    const float *color2,
                        float tube_size, int half_state, unsigned int b1, unsigned int b2, int a, bool b1masked, bool b2masked)
{
  float d[3], t[3], p0[3], p1[3], p2[3];
  int a3;
  int double_sided;

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
    const float* vv = coord + 3 * a3;
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

    t[0] = p2[0] * tube_size * 2;
    t[1] = p2[1] * tube_size * 2;
    t[2] = p2[2] * tube_size * 2;
  RepLine(cgo, s1, s2, isRamped, v1, v2, color1, b1, b2, a, color2, b1masked, b2masked);
  if (s1){
    CGOColorv(cgo, color1);
    CGOPickColor(cgo, b1, b1masked ? cPickableNoPick : a);
    float f[] = { 0.14F, 0.4F } ;
    float f_1[] = { 1.0F - f[0], 1.0F - f[1] };
    float pt1[] = { (f_1[0] * v1[0] + f[0] * v2[0]),
                    (f_1[0] * v1[1] + f[0] * v2[1]),
                    (f_1[0] * v1[2] + f[0] * v2[2]) };
    float pt2[] = { (f_1[1] * v1[0] + f[1] * v2[0]),
                    (f_1[1] * v1[1] + f[1] * v2[1]),
                    (f_1[1] * v1[2] + f[1] * v2[2]) };
    float p1[3], p2[3];
    subtract3f(pt1, t, p1);
    subtract3f(pt2, t, p2);
    cgo->add<cgo::draw::line>(p1, p2);
    if(double_sided) {
      add3f(pt1, t, p1);
      add3f(pt2, t, p2);
      cgo->add<cgo::draw::line>(p1, p2);
    }
  }
  if (s2){
    CGOColorv(cgo, color2);
    CGOPickColor(cgo, b2, b2masked ? cPickableNoPick : a);
    float f[] = { 0.6F, 0.86F } ;
    float f_1[] = { 1.0F - f[0], 1.0F - f[1] };
    float pt1[] = { (f_1[0] * v1[0] + f[0] * v2[0]),
                    (f_1[0] * v1[1] + f[0] * v2[1]),
                    (f_1[0] * v1[2] + f[0] * v2[2]) };
    float pt2[] = { (f_1[1] * v1[0] + f[1] * v2[0]),
                    (f_1[1] * v1[1] + f[1] * v2[1]),
                    (f_1[1] * v1[2] + f[1] * v2[2]) };
    float p1[3], p2[3];
    subtract3f(pt1, t, p1);
    subtract3f(pt2, t, p2);
    cgo->add<cgo::draw::line>(p1, p2);
    if(double_sided) {
      add3f(pt1, t, p1);
      add3f(pt2, t, p2);
      cgo->add<cgo::draw::line>(p1, p2);
    }
  }
}

RepWireBond::~RepWireBond()
{
  CGOFree(shaderCGO);
  CGOFree(primitiveCGO);
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
  PyMOLGlobals *G = cs->G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    int active = false;
    ObjectMolecule *obj = cs->Obj;
    float line_width, line_width_setting =
      SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_line_width);
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
      const BondType *bd = obj->Bond.data();
      const AtomInfoType *ai = obj->AtomInfo.data();
      int last_color = -9;
      const float *coord = cs->Coord.data();

      for(a = 0; a < nBond; a++) {
        int b1 = bd->index[0];
        int b2 = bd->index[1];
        const AtomInfoType *ai1, *ai2;
        bd++;
        if(GET_BIT((ai1 = ai + b1)->visRep,cRepLine) && GET_BIT((ai2 = ai + b2)->visRep,cRepLine)) {
          active = true;
          int a1 = cs->atmToIdx(b1);
          int a2 = cs->atmToIdx(b2);
          if((a1 >= 0) && (a2 >= 0)) {
            int c1 = ai1->color;
            int c2 = ai2->color;

            const float *v1 = coord + 3 * a1;
            const float *v2 = coord + 3 * a2;

            if(c1 == c2) {      /* same colors -> one line */
              if(c1 != last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G, c1));
              }
              glVertex3fv(v1);
              glVertex3fv(v2);  /* we done */
            } else {            /* different colors -> two lines */
              float avg[3];
              average3f(v1, v2, avg);

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

static int RepWireBondCGOGenerate(RepWireBond * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;
  CGO *convertcgo = NULL;
  int ok = true;
  short line_as_cylinders = 0;
  line_as_cylinders = SettingGetGlobal_b(G, cSetting_use_shaders) && SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_line_as_cylinders);

  {
    if (ok && I->primitiveCGO){
      if (line_as_cylinders){
        CGO *tmpCGO = CGONew(G);

        if (ok) ok &= CGOEnable(tmpCGO, GL_CYLINDER_SHADER);
        if (ok) ok &= CGOSpecial(tmpCGO, CYLINDER_WIDTH_FOR_REPWIRE);
        convertcgo = CGOConvertLinesToCylinderShader(I->primitiveCGO, tmpCGO);
        I->shaderCGO_has_cylinders = true;

        if (ok) ok &= CGOAppendNoStop(tmpCGO, convertcgo);

        if (ok) ok &= CGODisable(tmpCGO, GL_CYLINDER_SHADER);
        if (ok) ok &= CGOStop(tmpCGO);
        CGOFreeWithoutVBOs(convertcgo);
        convertcgo = tmpCGO;
      } else {
        bool trilines = SettingGetGlobal_b(G, cSetting_trilines);
        CGO *tmpCGO = CGONew(G), *tmp2CGO;
        int shader = trilines ? GL_TRILINES_SHADER : GL_LINE_SHADER;

        if (ok) ok &= CGOEnable(tmpCGO, shader);
        if (ok) ok &= CGODisable(tmpCGO, CGO_GL_LIGHTING);
        if (trilines) {
          if (ok) ok &= CGOSpecial(tmpCGO, LINEWIDTH_DYNAMIC_WITH_SCALE);
          tmp2CGO = CGOConvertToTrilinesShader(I->primitiveCGO, tmpCGO);
        } else {
          tmp2CGO = CGOConvertToLinesShader(I->primitiveCGO, tmpCGO);
        }
        if (ok) ok &= CGOAppendNoStop(tmpCGO, tmp2CGO);
        if (ok) ok &= CGODisable(tmpCGO, shader);
        if (ok) ok &= CGOStop(tmpCGO);
        CGOFreeWithoutVBOs(tmp2CGO);
        convertcgo = tmpCGO;
      }
      convertcgo->use_shader = true;
    }
    CGOFree(I->shaderCGO);
    I->shaderCGO = convertcgo;
    CHECKOK(ok, I->shaderCGO);
  }
  return ok;
}
	  
void RepWireBond::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  int ok = true;

  if(ray) {
#ifndef _PYMOL_NO_RAY
    CGORenderRay(primitiveCGO, ray, info, NULL, NULL, cs->Setting.get(), cs->Obj->Setting.get());
    ray->transparentf(0.0);
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    bool use_shader = SettingGetGlobal_b(G, cSetting_line_use_shader) &&
                      SettingGetGlobal_b(G, cSetting_use_shaders);
    if(pick) {
      CGORenderGLPicking(use_shader ? I->shaderCGO : I->primitiveCGO, info, &I->context, NULL, NULL, I);
    } else { /* else not pick i.e., when rendering */
      short line_as_cylinders ;
      line_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_line_as_cylinders);
      if (I->shaderCGO && (!use_shader || (line_as_cylinders ^ I->shaderCGO_has_cylinders))){
      CGOFree(I->shaderCGO);
        I->shaderCGO_has_cylinders = 0;
	}
	if (ok){
        if (use_shader) {
          if (!I->shaderCGO)
            ok &= RepWireBondCGOGenerate(I, info);
          CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, I);
	    } else {
          CGORenderGL(I->primitiveCGO, NULL, NULL, NULL, info, I);
	}
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->invalidate(cRepInvPurge);
    I->cs->Active[cRepLine] = false;
  }
}

static
bool IsBondTerminal(ObjectMolecule *obj, int b1, int b2){
  auto const has_heavy_neighbors = [obj](int atm, int atm_other) {
    for (auto const& neighbor : AtomNeighbors(obj, atm)) {
      if (neighbor.atm != atm_other &&
          obj->AtomInfo[neighbor.atm].protons > 1) {
        return true;
      }
    }
    return false;
  };

  return !has_heavy_neighbors(b1, b2) || !has_heavy_neighbors(b2, b1);
}

static int RepWireZeroOrderBond(CGO *cgo, bool s1, bool s2,
    const float *v1,
    const float *v2,
    const float *rgb1,
    const float *rgb2,
                                unsigned int b1, unsigned int b2, int a, float dash_gap, float dash_length, bool b1masked, bool b2masked)
{
  int ok = true;
  float axis[3], naxis[3];
  subtract3f(v2, v1, axis);
  copy3f(axis, naxis);
  normalize3f(naxis);
  float blen = length3f(axis);
  float dash_tot = dash_gap + dash_length;
  int ndashes = blen / dash_tot;

  // only do even number of dashes
  if (ndashes < 2) {
    ndashes = 2;
  } else if (ndashes % 2) {
    --ndashes;
  }

  float remspace = blen - (ndashes * dash_length); // remaining space for first gaps
  float dgap = remspace / (ndashes - 1.f); // endpoints at each vertex, therefore only account for ndashes-1 spaces
  float placep[3], placep2[3], adddlen[3], adddtot[3];
  float dplace;
  int ndashes_drawn = 0;
  bool color2_set = false;
  mult3f(naxis, dash_length, adddlen); // adddlen - length of dash as x/y/z vector
  mult3f(naxis, dash_length + dgap, adddtot); // adddtot - length of dash plus gap as x/y/z vector

  copy3f(v1, placep);
  if (s1){
    ok &= CGOColorv(cgo, rgb1);
    ok &= CGOPickColor(cgo, b1, b1masked ? cPickableNoPick : a);
    for (dplace = 0.f; (dplace+dash_length) < blen / 2.f; ){
      add3f(placep, adddlen, placep2);
      cgo->add<cgo::draw::line>(placep, placep2);
      add3f(placep, adddtot, placep);
      dplace += dash_length + dgap;
      ++ndashes_drawn;
    }
    if (!s2){
      if (dplace < blen / 2.f){
        // if we are behind the mid-point, only s1, so draw a half-bond
        add3f(placep, adddlen, placep2);
        cgo->add<cgo::draw::line>(placep, placep2);
        add3f(placep, adddtot, placep);
        dplace += dash_length + dgap;
        ++ndashes_drawn;
      }
    }
  } else {
    float tmpp[3];
    dplace = (ndashes/2) * (dash_length + dgap);
    mult3f(naxis, dplace, tmpp);
    add3f(v1, tmpp, placep);
    ndashes_drawn = ndashes/2;
    // if !s1, then definitely s2, so draw half-bond
    if (dplace <= blen / 2.f){
      // if no s1, and we are behind the mid-point, draw half-bond with only s2
      add3f(placep, adddlen, placep2);
      ok &= CGOColorv(cgo, rgb2);
      ok &= CGOPickColor(cgo, b2, b2masked ? cPickableNoPick : a);
      color2_set = true;
      cgo->add<cgo::draw::line>(placep, placep2);
      add3f(placep, adddtot, placep);
      dplace += dash_length + dgap;
      ++ndashes_drawn;
    }
  }
  if (s2){
    if (dplace < blen / 2.f){
      // if we are behind the mid-point, draw a split cylinder with both colors
      float tmpp[3];
      mult3f(axis, .5f, tmpp);
      add3f(v1, tmpp, tmpp);
      cgo->add<cgo::draw::line>(placep, tmpp);
      add3f(placep, adddlen, placep2);
      if (!color2_set){
        ok &= CGOColorv(cgo, rgb2);
        ok &= CGOPickColor(cgo, b2, b2masked ? cPickableNoPick : a);
      }
      cgo->add<cgo::draw::line>(tmpp, placep2);
      add3f(placep, adddtot, placep);
      dplace += dash_length + dgap;
      ++ndashes_drawn;
    } else if (!color2_set){
      ok &= CGOColorv(cgo, rgb2);
      ok &= CGOPickColor(cgo, b2, b2masked ? cPickableNoPick : a);
    }
    while (ndashes_drawn < ndashes){
      add3f(placep, adddlen, placep2);
      cgo->add<cgo::draw::line>(placep, placep2);
      add3f(placep, adddtot, placep);
      dplace += dash_length + dgap;
      ++ndashes_drawn;
    }
  }

  return ok;
}

Rep *RepWireBondNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj = cs->Obj;
  int a1, a2;
  unsigned int b1, b2;
  int a, c1, c2, ord;
  bool s1, s2;
  const BondType *b;
  int half_bonds, *other = NULL;
  float valence;
  const float *v1, *v2;
  int visFlag;
  float line_width;
  int valence_flag = false;
  const AtomInfoType *ai1;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 0;
  int line_stick_helper = 0;
  int na_mode;
  bool *marked = NULL;
  int valence_found = false;
  int variable_width = false;
  float last_line_width = -1.f;
  int line_color;
  int hide_long = false;
  int fancy;
  const float _0p9 = 0.9F;
  int ok = true;
  unsigned int line_counter = 0;
  PRINTFD(G, FB_RepWireBond)
    " RepWireBondNew-Debug: entered.\n" ENDFD;

  visFlag = false;
  b = obj->Bond;
  ai1 = obj->AtomInfo;
  if(ok && GET_BIT(obj->RepVisCache,cRepLine)){
    for(a = 0; a < obj->NBond; a++) {
      b1 = b->index[0];
      b2 = b->index[1];
      if((cRepLineBit & ai1[b1].visRep & ai1[b2].visRep)) {
        visFlag = true;
        break;
      }
      b++;
    }
  }
  if(!visFlag) {
    return (NULL);              /* skip if no dots are visible */
  }
  marked = pymol::calloc<bool>(obj->NAtom);
  CHECKOK(ok, marked);
  if (!ok){
    return (NULL);
  }
  
  valence_flag = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence);
  valence = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence_size);
  cartoon_side_chain_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
					   cSetting_cartoon_side_chain_helper);
  ribbon_side_chain_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
					  cSetting_ribbon_side_chain_helper);
  line_stick_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
				   cSetting_line_stick_helper);
  line_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_line_color);
  line_width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_line_width);
  
  if(line_stick_helper && (SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
					cSetting_stick_transparency) > R_SMALL4))
    line_stick_helper = false;
  half_bonds = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_half_bonds);
  hide_long = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_hide_long_bonds);
  na_mode =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_cartoon_nucleic_acid_mode);
  int na_mode_ribbon =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_nucleic_acid_mode);
  fancy = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence_mode) == 1;
  auto valence_zero_mode =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence_zero_mode);
  
  b = obj->Bond;
  
  for(a = 0; a < obj->NBond; a++) {
    b1 = b->index[0];
    b2 = b->index[1];

    a1 = cs->atmToIdx(b1);
    a2 = cs->atmToIdx(b2);
    if((a1 >= 0) && (a2 >= 0)) {
      if(!variable_width)
        if (AtomInfoCheckBondSetting(G, b, cSetting_line_width)){
          variable_width = true;
          if (valence_found) break;
        }
      auto bd_valence_flag = BondSettingGetWD(G, b, cSetting_valence, valence_flag);
      if(bd_valence_flag) {
        valence_found = true;
        if (variable_width) break;
      }
    }
    b++;
  }

  auto I = new RepWireBond(cs, state);

  I->primitiveCGO = CGONew(G);

  CGOSpecialWithArg(I->primitiveCGO, LINE_LIGHTING, 0.f);

  CGOSpecial(I->primitiveCGO, LINEWIDTH_FOR_LINES);
  CGOBegin(I->primitiveCGO, GL_LINES);

  if(obj->NBond) {

    if(valence_found)           /* build list of up to 2 connected atoms for each atom */
      other = ObjectMoleculeGetPrioritizedOtherIndexList(obj, cs);

    if(ok && (cartoon_side_chain_helper || ribbon_side_chain_helper)) {
      SideChainHelperMarkNonCartoonBonded(marked, obj, cs,
          cartoon_side_chain_helper,
          ribbon_side_chain_helper);
    }


    b = obj->Bond;

    for(a = 0; ok && a < obj->NBond; ++a, ++b) {
      b1 = b->index[0];
      b2 = b->index[1];
      ord = b->order;

      if (ord == 0 && valence_zero_mode == 0)
        continue;

      a1 = cs->atmToIdx(b1);
      a2 = cs->atmToIdx(b2);
      if((a1 >= 0) && (a2 >= 0)) {

        AtomInfoType *ati1 = obj->AtomInfo + b1;
        AtomInfoType *ati2 = obj->AtomInfo + b2;

        s1 = (ati1->visRep & cRepLineBit);
        s2 = (ati2->visRep & cRepLineBit);

        if(s1 ^ s2){
          if(!half_bonds) {
            if(line_stick_helper &&
               (((!s1) && (cRepCylBit & ati1->visRep) && !(cRepCylBit & ati2->visRep)) ||
                ((!s2) && (cRepCylBit & ati2->visRep) && !(cRepCylBit & ati1->visRep))))
              s1 = s2 = 1;      /* turn on line when both stick and line are alternately shown */
            else {
              s1 = 0;
              s2 = 0;
            }
          }
        }

        auto const s1_before_symop = s1;
        auto const s2_before_symop = s2;
        int symop_pass = 0;

        pymol::SymOp symop[2] = {pymol::SymOp(), b->symop_2};
        assert(!symop[0]);
        float vv_buf[2][3];

      inv_sym_bond:

        v1 = cs->coordPtrSym(a1, symop[0], vv_buf[0], symop_pass);
        v2 = cs->coordPtrSym(a2, symop[1], vv_buf[1], symop_pass);

        if (!v1 || !v2) {
          PRINTFB(G, FB_RepCylBond, FB_Warnings)
          " %s-Warning: Failed to get symmetry coordiantes\n",
              __func__ ENDFB(G);
          continue;
        }

        // show half-bond for atom which connects to a symmetry mate
        s1 = s1_before_symop && !symop[0];
        s2 = s2_before_symop && !symop[1];

        if(hide_long && (s1 || s2)) {
          float cutoff = (ati1->vdw + ati2->vdw) * _0p9;
          if(!within3f(v1, v2, cutoff)) /* atoms separated by more than 90% of the sum of their vdw radii */
            s1 = s2 = 0;
        }

        if(s1 || s2) {
          float rgb1[3], rgb2[3];
          int terminal = false;

          auto bd_valence_flag = BondSettingGetWD(G, b, cSetting_valence, valence_flag);
          auto bd_line_color = BondSettingGetWD(G, b, cSetting_line_color, line_color);

          if(fancy && bd_valence_flag && (b->order > 1)) {
            terminal = IsBondTerminal(obj, b1, b2);
          }

          if(variable_width) {
            auto bd_line_width = BondSettingGetWD(G, b, cSetting_line_width, line_width);
            if (last_line_width!=bd_line_width){
              CGOSpecialWithArg(I->primitiveCGO, LINEWIDTH_FOR_LINES, bd_line_width);
              last_line_width = bd_line_width;
            }
          }

          if(bd_line_color < 0) {
            if(bd_line_color == cColorObject) {
              c1 = (c2 = obj->Color);
            } else if(ColorCheckRamped(G, bd_line_color)) {
              c1 = (c2 = bd_line_color);
            } else {
              c1 = ati1->color;
              c2 = ati2->color;
            }
          } else {
            c1 = (c2 = bd_line_color);
          }

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

          bool isRamped = false;
          isRamped = ColorGetCheckRamped(G, c1, v1, rgb1, state);
          isRamped = ColorGetCheckRamped(G, c2, v2, rgb2, state) | isRamped;
          if (s1 || s2){
            if (ord == 0 && valence_zero_mode == 2) {
              ord = 1;
            }

            if (!ord){
              RepWireZeroOrderBond(I->primitiveCGO, s1, s2, v1, v2, rgb1, rgb2, b1, b2, a, .15f, .15f, ati1->masked, ati2->masked);
            } else if (!bd_valence_flag || ord <= 1){
              RepLine(I->primitiveCGO, s1, s2, isRamped, v1, v2, rgb1, b1, b2, a, rgb2, ati1->masked, ati2->masked);
            } else {
              if (ord == 4){
                RepAromatic(I->primitiveCGO, s1, s2, isRamped, v1, v2, other, a1, a2, cs->Coord, rgb1, rgb2, valence, 0, b1, b2, a, ati1->masked, ati2->masked);
              } else {
                RepValence(I->primitiveCGO, s1, s2, isRamped, v1, v2, other, a1, a2, cs->Coord, rgb1, rgb2, ord, valence, fancy && !terminal, b1, b2, a, ati1->masked, ati2->masked);
              }
            }
            line_counter++;
          }
        }

        // If this was a half-bond to a symmetry mate, do another pass and
        // render the other half.
        if (symop_pass == 0 && ati1 != ati2 && symop[1]) {
          symop_pass = 1;
          std::swap(symop[0], symop[1]);
          goto inv_sym_bond;
        }
      }
      ok &= !G->Interrupt;
    }
  }
  CGOEnd(I->primitiveCGO);
  CGOSpecialWithArg(I->primitiveCGO, LINE_LIGHTING, 1.f);
  CGOStop(I->primitiveCGO);
  FreeP(marked);
  FreeP(other);
  if (!ok || !line_counter){
    delete I;
    I = NULL;
  }
  return (Rep *) I;
}

