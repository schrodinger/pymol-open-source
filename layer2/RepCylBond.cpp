
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
#include"Err.h"
#include"Vector.h"
#include"ObjectMolecule.h"
#include"RepCylBond.h"
#include"SideChainHelper.h"
#include"Color.h"
#include"Setting.h"
#include"Feedback.h"
#include"ShaderMgr.h"
#include"Scene.h"
#include"CGO.h"
#include "Lex.h"

#include <iostream>

#ifdef _PYMOL_IOS
extern "C" void fireMemoryWarning();
#endif

struct RepCylBond : Rep {
  using Rep::Rep;

  ~RepCylBond() override;

  cRep_t type() const override { return cRepCyl; }
  void render(RenderInfo* info) override;

  CGO* primitiveCGO = nullptr;
  CGO* renderCGO = nullptr;
};

/* RepCylinder -- This function is a helper function that generates a cylinder for RepCylBond.
 *   Depending on the s1 and s2 arguments, it will either draw a bond or a half bond.  Half bonds
 *   always have flat caps on the inside.  The isRamped argument specifies whether the colors are
 *   interpolated only when s1 and s2 are set.
 *
 * PARAMS
 *
 * cgo - CGO that cylinder operation is added
 * s1, s2 - whether first(s1)/second(s2) half of cylinder is rendered
 * isRamped - only when s1 and s2, whether color interpolation is used
 * v1, v2 - x/y/z of both vertices of cylinder
 * frontCap/endCap - whether front/end of cylinder should be round or not
 * tube_size - radius of cylinder
 * v2color - second color specified as 3 floats (optional)
 * v2pickcolor - second pick color info (ptr to structure Pickable) specified (optional)
 *
 * RETURN VALUE: returns ok (whether adding operation(s) were successful
 *
 */
static int RepCylinder(CGO* cgo, bool s1, bool s2, bool isRamped,
    float const* v1, float const* v2, bool frontCap, bool endCap,
    float tube_size, float const* v2color = nullptr,
    Pickable* v2pickcolor = nullptr)
{
  float axis[3];
  int ok = true;
  subtract3f(v2, v1, axis);
  if (s1 && s2){
    short cap = (frontCap ? cCylShaderCap1Round : 0) |
                (endCap ? cCylShaderCap2Round : 0) |
                (isRamped ? cCylShaderInterpColor : 0);
    if (v2color){
      ok &= (bool)cgo->add<cgo::draw::shadercylinder2ndcolor>(cgo, v1, axis, tube_size, cap, v2color, v2pickcolor);
    } else {
      ok &= (bool)cgo->add<cgo::draw::shadercylinder>(v1, axis, tube_size, cap);
    }
  } else {
    // if either s1 or s2 is 0, then draw half bond
    mult3f(axis, .5f, axis);
    if (s1) {
      short cap = (frontCap ? cCylShaderCap1Round : 0) | cCylShaderCap2Flat;
      ok &= (bool)cgo->add<cgo::draw::shadercylinder>(v1, axis, tube_size, cap);
    } else if (s2){
      short cap = (endCap ? cCylShaderCap2Round : 0) | cCylShaderCap1Flat;
      float v1new[3];
      add3f(v1, axis, v1new);
      if (v2color){
        ok &= CGOColorv(cgo, v2color);
      }
      if (v2pickcolor){
        ok &= CGOPickColor(cgo, v2pickcolor->index, v2pickcolor->bond);
      }
      ok &= (bool)cgo->add<cgo::draw::shadercylinder>(v1new, axis, tube_size, cap);
    }
  }
  return ok;
}

RepCylBond::~RepCylBond()
{
  auto I = this;
  CGOFree(I->primitiveCGO);
  CGOFree(I->renderCGO);
}

static int RepCylBondCGOGenerate(RepCylBond * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->G;

  const CGO* input = I->primitiveCGO;
  assert(input);

  bool const use_shader = info->use_shaders && //
                          SettingGet<bool>(*I->cs, cSetting_stick_use_shader);
  bool const as_cylinders =
      use_shader && //
      SettingGet<bool>(*I->cs, cSetting_stick_as_cylinders) &&
      SettingGet<bool>(*I->cs, cSetting_render_as_cylinders) &&
      G->ShaderMgr->ShaderPrgExists("cylinder");

  std::unique_ptr<CGO> convertcgo;

  if (as_cylinders) {
    convertcgo.reset(CGONew(G));

    CGOEnable(convertcgo.get(), GL_CYLINDER_SHADER);
    std::unique_ptr<CGO> cylindercgo(
        CGOConvertShaderCylindersToCylinderShader(input, convertcgo.get()));
    convertcgo->move_append(std::move(*cylindercgo));
    CGODisable(convertcgo.get(), GL_CYLINDER_SHADER);

    std::unique_ptr<CGO> spherescgo(
        CGOOptimizeSpheresToVBONonIndexed(input, 0, true));
    if (spherescgo) {
      convertcgo->move_append(std::move(*spherescgo));
    }
  } else {
    std::unique_ptr<CGO> simplified(CGOSimplify(input, 0, //
        SettingGet<int>(G, cSetting_cgo_sphere_quality),
        SettingGet<int>(G, cSetting_stick_round_nub)));
    p_return_val_if_fail(simplified, false);
    if (use_shader) {
      convertcgo.reset(CGOOptimizeToVBONotIndexed(simplified.get()));
    } else {
      convertcgo.reset(CGOCombineBeginEnd(simplified.get()));
    }
  }

  p_return_val_if_fail(convertcgo, false);

  assert(!I->renderCGO);
  I->renderCGO = convertcgo.release();
  CGOSetUseShader(I->renderCGO, use_shader);

  return true;
}

void RepCylBond::render(RenderInfo * info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  PyMOLGlobals *G = I->G;
  int ok = true;

  if(ray) {
#ifndef _PYMOL_NO_RAY
    CGORenderRay(I->primitiveCGO, ray, info, NULL, NULL, I->cs->Setting.get(), I->obj->Setting.get());
    ray->transparentf(0.0);
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    bool use_shader = SettingGetGlobal_b(G, cSetting_stick_use_shader)
      && SettingGetGlobal_b(G, cSetting_use_shaders);

    if (I->renderCGO && (CGOCheckWhetherToFree(G, I->renderCGO) || ((bool)I->renderCGO->use_shader) != use_shader)){
      CGOFree(I->renderCGO);
      I->renderCGO = 0;
    }

    if(pick) {
      PRINTFD(G, FB_RepCylBond)
        " RepCylBondRender: rendering pickable...\n" ENDFD;

      if (I->renderCGO){
        CGORenderGLPicking(I->renderCGO, info, &I->context, I->cs->Setting.get(), I->obj->Setting.get());
      }
    } else { /* else not pick, i.e., when rendering */
      if (!I->renderCGO){
        ok &= RepCylBondCGOGenerate(I, info);
        assert(I->renderCGO);
      }
      const float *color = ColorGet(G, I->obj->Color);
      I->renderCGO->debug = SettingGetGlobal_i(G, cSetting_stick_debug);
      CGORenderGL(I->renderCGO, color, NULL, NULL, info, I);
    }
  }
}

static int RepZeroOrderBond(RepCylBond* I, CGO* cgo, bool s1, bool s2,
    const float* vv1, const float* vv2, float zradius, const float* rgb1,
    const float* rgb2, unsigned int b1, unsigned int b2, int a, bool b1masked,
    bool b2masked)
{
  float axis[3], naxis[3];
  subtract3f(vv2, vv1, axis);
  copy3f(axis, naxis);
  normalize3f(naxis);
  float blen = length3f(axis);
  float dgap = zradius*6.f, dlen = zradius*3.f; // dlen - dash length, dgap - gap length
  float placep[3], placep2[3], adddlen[3], adddtot[3];
  float dplace;
  int ndashes = blen / (dlen + dgap);

  // only do even number of dashes
  if (ndashes < 2) {
    ndashes = 2;
  } else if (ndashes % 2) {
    --ndashes;
  }

  float remspace = blen - (ndashes * dlen); // remaining space for first gaps
  dgap = remspace / (ndashes + 1.f);

  int ok = true;
  mult3f(naxis, dlen, adddlen); // adddlen - length of dash as x/y/z vector
  mult3f(naxis, dlen + dgap, adddtot); // adddtot - length of dash plus gap as x/y/z vector
  mult3f(naxis, dgap, placep);
  add3f(vv1, placep, placep);
  if (s1){
    ok &= CGOColorv(I->primitiveCGO, rgb1);
    ok &= CGOPickColor(I->primitiveCGO, b1, b1masked ? cPickableNoPick : a);
    for (dplace = dgap; (dplace+dlen) < blen / 2.f; ){
      add3f(placep, adddlen, placep2);
      ok &= RepCylinder(I->primitiveCGO, true, true, false, placep, placep2, true, true, zradius);
      add3f(placep, adddtot, placep);
      dplace += dlen + dgap;
    }
    if (!s2){
      if (dplace < blen / 2.f){
        // if we are behind the mid-point, only s1, so draw a half-bond
        add3f(placep, adddlen, placep2);
        ok &= RepCylinder(I->primitiveCGO, s1, s2, false, placep, placep2, true, false, zradius);
        add3f(placep, adddtot, placep);
        dplace += dlen + dgap;
      }
    }
  } else {
    float tmpp[3];
    dplace = dgap + (ndashes/2) * (dlen + dgap);
    mult3f(naxis, dplace, tmpp);
    add3f(vv1, tmpp, placep);
    // if !s1, then definitely s2, so draw half-bond
    if (dplace < blen / 2.f){
      // if no s1, and we are behind the mid-point, draw half-bond with only s2
      add3f(placep, adddlen, placep2);
      ok &= CGOColorv(I->primitiveCGO, rgb2);
      ok &= CGOPickColor(I->primitiveCGO, b2, b2masked ? cPickableNoPick : a);
      ok &= RepCylinder(I->primitiveCGO, s1, s2, false, placep, placep2, false, true, zradius);
      add3f(placep, adddtot, placep);
      dplace += dlen + dgap;
    }
  }
  if (s2){
    if (dplace < blen / 2.f){
      // if we are behind the mid-point, draw a split cylinder with both colors
      add3f(placep, adddlen, placep2);
      Pickable pickcolor2 = { b2, b2masked ? cPickableNoPick : a };
      ok &= RepCylinder(I->primitiveCGO, true, true, false, placep, placep2, true, true, zradius, rgb2, &pickcolor2);
      add3f(placep, adddtot, placep);
      dplace += dlen + dgap;
    }
    ok &= CGOColorv(I->primitiveCGO, rgb2);
    ok &= CGOPickColor(I->primitiveCGO, b2, b2masked ? cPickableNoPick : a);
    for (; (dplace + dlen) < blen; ){
      add3f(placep, adddlen, placep2);
      ok &= RepCylinder(I->primitiveCGO, true, true, false, placep, placep2, true, true, zradius);
      add3f(placep, adddtot, placep);
      dplace += dlen + dgap;
    }
  }
  return ok;
}

static int RepValence(RepCylBond *I, CGO *cgo, bool s1, bool s2, bool isRamped,
		      const float *v1, const float *v2, const int *other,
		      int a1, int a2, const float *coord,
		      const float *color1, const float *color2, int ord,
		      float tube_size,
		      bool fixed_r, float scale_r,
		      Pickable pickdata[] = NULL)
{
  float d[3], t[3], p0[3], p1[3], p2[3];
  const float* vv;
  float v1t[3], v2t[3];
  int a3;
  int double_sided;
  int ok = true;

  /* First, we need to construct a coordinate system */

  /* get direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  copy3f(p0, d);
  normalize3f(p0);

  /* need a third atom to get planarity */
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

  /* we have a coordinate system */

  /* Next, we need to determine how many cylinders */
  switch (ord) {
  case 2:
    {
      float radius = tube_size;
      if(!fixed_r) {
        radius *= scale_r;
        radius /= 2.5;
      }

      t[0] = p2[0] * 1.5F * radius;
      t[1] = p2[1] * 1.5F * radius;
      t[2] = p2[2] * 1.5F * radius;

      add3f(v1, t, v1t);
      add3f(v2, t, v2t);
      if (ok)
        ok &= CGOColorv(cgo, color1);
      if (ok && pickdata)
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);
      Pickable *pickdataptr = NULL;
      if (pickdata)
        pickdataptr = &pickdata[1];
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1t, v2t, true, true, radius, color2, pickdataptr);
      subtract3f(v1, t, v1t);
      subtract3f(v2, t, v2t);
      if (ok)
        ok &= CGOColorv(cgo, color1);
      if (ok && pickdata)
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1t, v2t, true, true, radius, color2, pickdataptr);
    }
    break;
  case 3:
    {
      float radius = tube_size;
      if(!fixed_r) {
        radius *= scale_r;
        radius /= 3.5;
      }

      t[0] = p2[0] * 2.5F * radius;
      t[1] = p2[1] * 2.5F * radius;
      t[2] = p2[2] * 2.5F * radius;

      if (ok)
        ok &= CGOColorv(cgo, color1);
      if (ok && pickdata)
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);

      copy3f(v1, v1t);
      copy3f(v2, v2t);
      Pickable *pickdataptr = NULL;
      if (pickdata)
        pickdataptr = &pickdata[1];
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1t, v2t, true, true, radius, color2, pickdataptr);
      add3f(v1, t, v1t);
      add3f(v2, t, v2t);
      if (ok && pickdata){
        ok &= CGOColorv(cgo, color1);
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);
      }
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1t, v2t, true, true, radius, color2, pickdataptr);
      subtract3f(v1, t, v1t);
      subtract3f(v2, t, v2t);
      if (ok && pickdata){
        ok &= CGOColorv(cgo, color1);
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);
      }
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1t, v2t, true, true, radius, color2, pickdataptr);
    }
    break;
  case 4:
    {
      float radius = tube_size;
      float radius2 = tube_size;
      float along[3], adj[3], v1tt[3], v2tt[3];
      float inner1a = 0.24F;
      float inner1b = 0.44F;
      float inner2a = 0.5F + (0.5F - inner1b);
      float inner2b = 1.0F - inner1a;

      if(!fixed_r) {
        radius *= scale_r;
        radius2 = radius / 2.5F;
        t[0] = p2[0] * 1.5F * radius;
        t[1] = p2[1] * 1.5F * radius;
        t[2] = p2[2] * 1.5F * radius;
      } else {
        inner1a -= 0.04F;
        inner2b = 1.0F - inner1a;
        t[0] = p2[0] * 3 * radius;
        t[1] = p2[1] * 3 * radius;
        t[2] = p2[2] * 3 * radius;
      }

      if (ok)
        ok &= CGOColorv(cgo, color1);
      if (ok && pickdata)
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);

      subtract3f(v1, t, v1t);
      subtract3f(v2, t, v2t);
      subtract3f(v2t, v1t, along);
      scale3f(along, inner1a, adj);
      add3f(adj, v1t, v1tt);
      scale3f(along, inner1b, adj);
      add3f(adj, v1t, v2tt);
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1tt, v2tt, true, true, radius2);
      if(double_sided) {
        add3f(v1, t, v1t);
        add3f(v2, t, v2t);
        subtract3f(v2t, v1t, along);
        scale3f(along, inner1a, adj);
        add3f(adj, v1t, v1tt);
        scale3f(along, inner1b, adj);
        add3f(adj, v1t, v2tt);
        if (ok)
          ok &= RepCylinder(cgo, s1, s2, isRamped, v1tt, v2tt, true, true, radius2);
      }

      if (ok)
        ok &= CGOColorv(cgo, color2);
      if (ok && pickdata)
        ok &= CGOPickColor(cgo, pickdata[1].index, pickdata[1].bond);
      subtract3f(v1, t, v1t);
      subtract3f(v2, t, v2t);
      subtract3f(v2t, v1t, along);
      scale3f(along, inner2a, adj);
      add3f(adj, v1t, v1tt);
      scale3f(along, inner2b, adj);
      add3f(adj, v1t, v2tt);
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1tt, v2tt, true, true, radius2);
      if(double_sided) {
        add3f(v1, t, v1t);
        add3f(v2, t, v2t);
        subtract3f(v2t, v1t, along);
        scale3f(along, inner2a, adj);
        add3f(adj, v1t, v1tt);
        scale3f(along, inner2b, adj);
        add3f(adj, v1t, v2tt);
        if (ok)
          ok &= RepCylinder(cgo, s1, s2, isRamped, v1tt, v2tt, true, true, radius2);
      }

      Pickable *pickdataptr = NULL;
      if (ok)
        ok &= CGOColorv(cgo, color1);
      if (ok && pickdata){
        ok &= CGOPickColor(cgo, pickdata[0].index, pickdata[0].bond);
        pickdataptr = &pickdata[1];
      }
      if (ok)
        ok &= RepCylinder(cgo, s1, s2, isRamped, v1, v2, true, true, radius, color2, pickdataptr);
    }
    break;
  }
  return ok;
}

Rep *RepCylBondNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj;
  int a1, a2; // can be -1 for missing atoms
  int a;
  int c1, c2, s1, s2;
  unsigned int b1, b2;
  const BondType *b;
  float radius;
  float valence;
  int half_bonds, *other = NULL;
  int visFlag;
  int ord;
  int stick_ball, stick_ball_color = -1;
  float stick_ball_ratio = 1.0F;
  const AtomInfoType *ai1;
  bool fixed_radius = false;
  int valence_flag = false;
  int hide_long = false;
  int stick_color = 0;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 1;
  int na_mode;
  bool *marked = NULL;
  float *capdrawn = NULL;
  float scale_r = 1.0F;
  float transp, h_scale;
  float prev_transp = -1;
  int valence_found = false;
  const float _0p9 = 0.9F;
  short shader_mode = 0;
  int ok = true;

  PRINTFD(G, FB_RepCylBond)
    " RepCylBondNew-Debug: entered.\n" ENDFD;
  obj = cs->Obj;
  visFlag = false;
  b = obj->Bond;
  ai1 = obj->AtomInfo;
  if(obj->RepVisCache & cRepCylBit)
    for(a = 0; a < obj->NBond; a++) {
      b1 = b->index[0];
      b2 = b->index[1];
      if((cRepCylBit & ai1[b1].visRep & ai1[b2].visRep)) {
	visFlag = true;
	break;
      }
      b++;
    }
  if(!visFlag) {
    return (NULL);              /* skip if no sticks are visible */
  }

  capdrawn = pymol::calloc<float>(obj->NAtom); // max radius of caps
  marked = pymol::calloc<bool>(obj->NAtom);
  CHECKOK(ok, marked);
  if (!ok){
    return NULL;
  }

  valence = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence);
  valence_flag = (valence != 0.0F);

  stick_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_color);
  cartoon_side_chain_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                           cSetting_cartoon_side_chain_helper);
  ribbon_side_chain_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                          cSetting_ribbon_side_chain_helper);

  transp = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_transparency);
  hide_long = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_hide_long_bonds);

  std::set<int> all_zero_order_bond_atoms;
  b = obj->Bond;
  for(a = 0; ok && a < obj->NBond; a++) {
    b1 = b->index[0];
    b2 = b->index[1];
    ord = b->order;
    a1 = cs->atmToIdx(b1);
    a2 = cs->atmToIdx(b2);

    if((a1 >= 0) && (a2 >= 0)) {
      AtomInfoType *ati1 = obj->AtomInfo + b1;
      AtomInfoType *ati2 = obj->AtomInfo + b2;
      s1 = GET_BIT(ati1->visRep, cRepCyl);
      s2 = GET_BIT(ati2->visRep, cRepCyl);

      if (s1 && s2){
        if (!valence_found)
          valence_found = BondSettingGetWD(G, b, cSetting_valence, valence_flag);
        
        if (!ord){
          all_zero_order_bond_atoms.insert(b1);
          all_zero_order_bond_atoms.insert(b2);
        }
      }
    }
    b++;
    ok &= !G->Interrupt;
  }
  radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_radius);
  half_bonds = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_half_bonds);
  na_mode =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_cartoon_nucleic_acid_mode);
  int na_mode_ribbon =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_nucleic_acid_mode);
  h_scale = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_h_scale);
  auto valence_zero_scale =
    SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence_zero_scale);
  auto valence_zero_mode =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_valence_zero_mode);

  auto I = new RepCylBond(cs, state);

  I->primitiveCGO = CGONew(G);
  if(ok && obj->NBond) {
    stick_ball = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_ball);

    shader_mode = SettingGetGlobal_b(G, cSetting_use_shaders) && SettingGetGlobal_b(G, cSetting_stick_as_cylinders) 
      && SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_stick_use_shader);

    if(cartoon_side_chain_helper || ribbon_side_chain_helper) {
      SideChainHelperMarkNonCartoonBonded(marked, obj, cs,
          cartoon_side_chain_helper,
          ribbon_side_chain_helper);
    }

    if(valence_found) {         /* build list of up to 2 connected atoms for each atom */
      other = ObjectMoleculeGetPrioritizedOtherIndexList(obj, cs);
      CHECKOK(ok, other);
      if (ok){
        fixed_radius = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_fixed_radius);
        scale_r = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_valence_scale);
      }
    }

    /* spheres for stick & balls */
    stick_ball_ratio = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_ball_ratio);
    stick_ball_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_ball_color);

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
        float bd_radius_full;

        auto bd_stick_color = BondSettingGetWD(G, b, cSetting_stick_color, stick_color);
        auto bd_radius = BondSettingGetWD(G, b, cSetting_stick_radius, radius);

        // version <=1.8.2 used negative stick_radius to turn on
        // stick_h_scale, which had a default of 0.4 (now: 1.0)
        if(bd_radius < 0.0F) {
          bd_radius = -bd_radius;

          // legacy behavior
          if (h_scale == 1.0F)
            h_scale = 0.4F;
        }

        // before H scaling
        bd_radius_full = bd_radius;

        // scaling for bonds involving hydrogen
        if (ati1->isHydrogen() || ati2->isHydrogen())
          bd_radius *= h_scale;

        if(bd_stick_color < 0) {
          if(bd_stick_color == cColorObject) {
            c1 = (c2 = obj->Color);
          } else if(ColorCheckRamped(G, bd_stick_color)) {
            c1 = (c2 = bd_stick_color);
          } else {
            c1 = ati1->color;
            c2 = ati2->color;
          }
        } else {
          c1 = (c2 = bd_stick_color);
        }

        s1 = GET_BIT(ati1->visRep, cRepCyl);
        s2 = GET_BIT(ati2->visRep, cRepCyl);

        if (!(s1 || s2)) {
          continue;
        }

        if (!(s1 && s2) && !half_bonds) {
          continue;
        }

        auto const s1_before_symop = s1;
        auto const s2_before_symop = s2;
        int symop_pass = 0;

        pymol::SymOp symop[2] = {pymol::SymOp(), b->symop_2};
        assert(!symop[0]);
        float vv_buf[2][3];
        float const *vv1, *vv2;

      inv_sym_bond:

        vv1 = cs->coordPtrSym(a1, symop[0], vv_buf[0], symop_pass);
        vv2 = cs->coordPtrSym(a2, symop[1], vv_buf[1], symop_pass);

        if (!vv1 || !vv2) {
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
          if(!within3f(vv1, vv2, cutoff))       /* atoms separated by more than 90% of the sum of their vdw radii */
            s1 = s2 = 0;
        }

        // side chain helpers
        if ((s1 || s2) && (ati1->flags & ati2->flags & cAtomFlag_polymer)) {
          if ((cRepCartoonBit & ati1->visRep & ati2->visRep)) {
            bool sc_helper =
              AtomSettingGetWD(G, ati1,
                  cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper) ||
              AtomSettingGetWD(G, ati2,
                  cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper);
            if (sc_helper &&
                SideChainHelperFilterBond(G, marked, ati1, ati2, b1, b2, na_mode, &c1, &c2))
              s1 = s2 = 0;
          }

          if ((s1 || s2) && (cRepRibbonBit & ati1->visRep & ati2->visRep)) {
            bool sc_helper =
              AtomSettingGetWD(G, ati1,
                  cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper) ||
              AtomSettingGetWD(G, ati2,
                  cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper);
            if (sc_helper &&
                SideChainHelperFilterBond(G, marked, ati1, ati2, b1, b2, na_mode_ribbon, &c1, &c2))
              s1 = s2 = 0;
          }
        }

          /* This means that if stick_ball gets changed, the RepCylBond needs to be completely invalidated */

        auto stick_ball_impl = [&](AtomInfoType * ati1, int b1, int c1, float const* vv1) {
          int stick_ball_1 = AtomSettingGetWD(G, ati1, cSetting_stick_ball, stick_ball);
          if(stick_ball_1) {
            float vdw = stick_ball_ratio * ((ati1->protons == cAN_H) ? bd_radius : bd_radius_full);
            float vdw1 = (vdw >= 0) ? vdw : -ati1->vdw * vdw;

            // only draw cap if larger than any previously drawn (PYMOL-2527)
            if (vdw1 < capdrawn[b1])
              return;

            int sbc1 = (stick_ball_color == cColorDefault) ? c1 : stick_ball_color;
            float rgb1[3];
            if(sbc1 == cColorAtomic)
              sbc1 = ati1->color;
            capdrawn[b1] = vdw1;
            ColorGetCheckRamped(G, sbc1, vv1, rgb1, state);
            CGOColorv(I->primitiveCGO, rgb1);
            CGOPickColor(I->primitiveCGO, b1, ati1->masked ? cPickableNoPick : a);
            CGOSphere(I->primitiveCGO, vv1, vdw1);
          }
        };

        if(s1 || s2) {
          auto const bd_transp =
              BondSettingGetWD(G, b, cSetting_stick_transparency, transp);

          if (prev_transp != bd_transp) {
            prev_transp = bd_transp;
            CGOAlpha(I->primitiveCGO, 1.0F - bd_transp);

            if (bd_transp > 0) {
              I->setHasTransparency();
            }
          }

          if (s1) stick_ball_impl(ati1, b1, c1, vv1);
          if (s2) stick_ball_impl(ati2, b2, c2, vv2);

          float rgb1[3], rgb2[3];
          bool isRamped = false;
          isRamped = ColorGetCheckRamped(G, c1, vv1, rgb1, state);
          isRamped = ColorGetCheckRamped(G, c2, vv2, rgb2, state) | isRamped;

          if (ord == 0) {
            bd_radius *= valence_zero_scale;
            if (valence_zero_mode == 2) {
              ord = 1;
            }
          }

          if (!ord){
            // zero order bonds
            ok &= RepZeroOrderBond(I, I->primitiveCGO, s1, s2, vv1, vv2, bd_radius, rgb1, rgb2, b1, b2, a, ati1->masked, ati2->masked);
          } else {
            all_zero_order_bond_atoms.erase(b1);
            all_zero_order_bond_atoms.erase(b2);

            bool bd_valence_flag = (ord > 1) && (ord < 5) &&
              BondSettingGetWD(G, b, cSetting_valence, valence_flag);

            if(bd_valence_flag) {
              Pickable pickdata[] = { { b1, ati1->masked ? cPickableNoPick : a },
                                      { b2, ati2->masked ? cPickableNoPick : a } };
              ok &= RepValence(I, I->primitiveCGO, s1, s2, isRamped, vv1, vv2, other,
                               a1, a2, cs->Coord, rgb1, rgb2, ord,
                               bd_radius, fixed_radius, scale_r, pickdata);
            } else {
              ok &= CGOColorv(I->primitiveCGO, rgb1);
              ok &= CGOPickColor(I->primitiveCGO, b1, ati1->masked ? cPickableNoPick : a);
              /* generate a cylinder */
              if (ok){
                Pickable pickdata = { b2, ati2->masked ? cPickableNoPick : a };

                bool drawcap1 = bd_radius > capdrawn[b1] && s1;
                bool drawcap2 = bd_radius > capdrawn[b2] && s2;

                ok &= RepCylinder(I->primitiveCGO, s1, s2, isRamped, vv1, vv2,
                    drawcap1, drawcap2, bd_radius, rgb2, &pickdata);

                if (shader_mode) {
                  // don't render caps twice with the cylinder shader
                  if (drawcap1) capdrawn[b1] = bd_radius;
                  if (drawcap2) capdrawn[b2] = bd_radius;
                }
              }
            }
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
    }
    /* for all zero order bond atoms that do not have any other ordered bonds,
       a sphere should be rendered for it (a cylinder with vertices that are almost
       exactly the same is used to render a sphere so that we won't need to use 
       the sphere shader excessively. */
    for (auto at : all_zero_order_bond_atoms){
      ai1 = obj->AtomInfo + at;
      c1 = ai1->color;
      float *v1 = cs->coordPtr(cs->atmToIdx(at));
      float v2[3];
      float rgb1[3];
      float v[3] = {R_SMALL4,0.,0.};
      // instead of using the sphere shader, we use a cylinder with 0 length (or very close to zero,
      // since the shader uses the axis for many computations/checks
      add3f(v1, v, v2);

      ColorGetCheckRamped(G, c1, v1, rgb1, state);
      ok &= CGOColorv(I->primitiveCGO, rgb1);
      ok &= CGOPickColor(I->primitiveCGO, at, ai1->masked ? cPickableNoPick : cPickableAtom);
      ok &= RepCylinder(I->primitiveCGO, true, true, false, v1, v2, true, true, radius);
    }
  }
  FreeP(other);
  FreeP(marked);
  FreeP(capdrawn);

  CGOStop(I->primitiveCGO);
  if (!ok){
    delete I;
    I = NULL;
  }

  return (Rep *) I;
}

#ifndef PURE_OPENGL_ES_2
static void RepCylinderImmediate(const float *v1arg, const float *v2arg, int nEdge,
                                 int frontCapArg, int endCapArg,
                                 float overlap, float nub, float radius, float **dir)
{
  float d[3], t[3], p0[3], p1[3], p2[3], v1ptr[3], v2ptr[3], *v1, *v2;
  float v[3], vv[3], vvv[3];
  float x, y;
  int c, frontCap = frontCapArg, endCap = endCapArg, tmpCap;

  p0[0] = (v2arg[0] - v1arg[0]);
  p0[1] = (v2arg[1] - v1arg[1]);
  p0[2] = (v2arg[2] - v1arg[2]);

  normalize3f(p0);

  v1ptr[0] = v1arg[0]; v1ptr[1] = v1arg[1]; v1ptr[2] = v1arg[2];
  v2ptr[0] = v2arg[0]; v2ptr[1] = v2arg[1]; v2ptr[2] = v2arg[2];

  v1ptr[0] -= p0[0] * overlap;
  v1ptr[1] -= p0[1] * overlap;
  v1ptr[2] -= p0[2] * overlap;

  if(endCap) {
    v2ptr[0] += p0[0] * overlap;
    v2ptr[1] += p0[1] * overlap;
    v2ptr[2] += p0[2] * overlap;
  }

  v1 = v1ptr;
  v2 = v2ptr;

  d[0] = (v2[0] - v1[0]);
  d[1] = (v2[1] - v1[1]);
  d[2] = (v2[2] - v1[2]);

  if (dir){
    if (!*dir){
      *dir = pymol::malloc<float>(3);
      (*dir)[0] = d[0]; (*dir)[1] = d[1]; (*dir)[2] = d[2];
    } else {
      if (get_angle3f(d, *dir)>=(cPI/2.)){
	v1 = v2ptr;
	v2 = v1ptr;
	d[0] = -d[0];     d[1] = -d[1];     d[2] = -d[2];
	tmpCap = frontCap;
	frontCap = endCap;
	endCap = tmpCap;
      }
    }
  }

  /* direction vector */

  p0[0] = (v2[0] - v1[0]);
  p0[1] = (v2[1] - v1[1]);
  p0[2] = (v2[2] - v1[2]);

  normalize3f(p0);
  get_divergent3f(d, t);
  cross_product3f(d, t, p1);
  normalize3f(p1);
  cross_product3f(d, p1, p2);
  normalize3f(p2);

  /* now we have a coordinate system */
  glBegin(GL_TRIANGLE_STRIP);

  for(c = nEdge; c >= 0; c--) {
    x = (float) radius * cos(c * 2 * PI / nEdge);
    y = (float) radius * sin(c * 2 * PI / nEdge);
    v[0] = p1[0] * x + p2[0] * y;
    v[1] = p1[1] * x + p2[1] * y;
    v[2] = p1[2] * x + p2[2] * y;

    vv[0] = v1[0] + v[0];
    vv[1] = v1[1] + v[1];
    vv[2] = v1[2] + v[2];

    glNormal3fv(v);

    vvv[0] = vv[0] + d[0];
    vvv[1] = vv[1] + d[1];
    vvv[2] = vv[2] + d[2];

    glVertex3fv(vv);
    glVertex3fv(vvv);
  }
  glEnd();

  if(frontCap) {
    v[0] = -p0[0];
    v[1] = -p0[1];
    v[2] = -p0[2];

    vv[0] = v1[0] - p0[0] * nub;
    vv[1] = v1[1] - p0[1] * nub;
    vv[2] = v1[2] - p0[2] * nub;

    glBegin(GL_TRIANGLE_FAN);

    glNormal3fv(v);
    glVertex3fv(vv);

    for(c = nEdge; c >= 0; c--) {
      x = (float) radius * cos(c * 2 * PI / nEdge);
      y = (float) radius * sin(c * 2 * PI / nEdge);
      v[0] = p1[0] * x + p2[0] * y;
      v[1] = p1[1] * x + p2[1] * y;
      v[2] = p1[2] * x + p2[2] * y;

      vv[0] = v1[0] + v[0];
      vv[1] = v1[1] + v[1];
      vv[2] = v1[2] + v[2];

      glNormal3fv(v);
      glVertex3fv(vv);
    }

    glEnd();
  }

  if(endCap) {

    v[0] = p0[0];
    v[1] = p0[1];
    v[2] = p0[2];

    vv[0] = v2[0] + p0[0] * nub;
    vv[1] = v2[1] + p0[1] * nub;
    vv[2] = v2[2] + p0[2] * nub;

    glBegin(GL_TRIANGLE_FAN);

    glNormal3fv(v);
    glVertex3fv(vv);

    for(c = 0; c <= nEdge; c++) {
      x = (float) radius * cos(c * 2 * PI / nEdge);
      y = (float) radius * sin(c * 2 * PI / nEdge);
      v[0] = p1[0] * x + p2[0] * y;
      v[1] = p1[1] * x + p2[1] * y;
      v[2] = p1[2] * x + p2[2] * y;

      vv[0] = v2[0] + v[0];
      vv[1] = v2[1] + v[1];
      vv[2] = v2[2] + v[2];
      glNormal3fv(v);
      glVertex3fv(vv);
    }
    glEnd();
  }
}
#endif

void RepCylBondRenderImmediate(CoordSet * cs, RenderInfo * info)
{
#ifndef PURE_OPENGL_ES_2
  /* performance optimized, so it does not support the following:

     - anything other than opengl
     - display of bond valences
     - per-bond & per-atom properties
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
    int nEdge = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_quality);
    float radius =
      fabs(SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_radius));
    float overlap =
      SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_overlap);
    float nub = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_stick_nub);
    float overlap_r = radius * overlap;
    float nub_r = radius * nub;

    {
      int a;
      int nBond = obj->NBond;
      const BondType *bd = obj->Bond.data();
      const AtomInfoType *ai = obj->AtomInfo.data();
      int last_color = -9;
      const float *coord = cs->Coord.data();
      const float _pt5 = 0.5F;

      for(a = 0; a < nBond; a++) {
        int b1 = bd->index[0];
        int b2 = bd->index[1];
        const AtomInfoType *ai1, *ai2;
        bd++;

        if( ((ai1 = ai + b1)->visRep & cRepCylBit) &&
            ((ai2 = ai + b2)->visRep & cRepCylBit)) {
          int a1, a2;
          active = true;
          a1 = cs->atmToIdx(b1);
          a2 = cs->atmToIdx(b2);

          if((a1 >= 0) && (a2 >= 0)) {
            int c1 = ai1->color;
            int c2 = ai2->color;

            const float *v1 = coord + 3 * a1;
            const float *v2 = coord + 3 * a2;

            if(c1 == c2) {      /* same colors -> one cylinder */
              if(c1 != last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G, c1));
              }

	      /* overlap is half since it is one cylinder representing both halfs of a bond */
              RepCylinderImmediate(v1, v2, nEdge, 1, 1, overlap_r, nub_r, radius, NULL);

            } else {            /* different colors -> two cylinders, no interior */
              float avg[3], *dir = NULL;

              avg[0] = (v1[0] + v2[0]) * _pt5;
              avg[1] = (v1[1] + v2[1]) * _pt5;
              avg[2] = (v1[2] + v2[2]) * _pt5;

              if(c1 != last_color) {
                last_color = c1;
                glColor3fv(ColorGet(G, c1));
              }

              RepCylinderImmediate(v1, avg, nEdge, 1, 0, overlap_r, nub_r, radius, &dir);

              if(c2 != last_color) {
                last_color = c2;
                glColor3fv(ColorGet(G, c2));
              }

              RepCylinderImmediate(v2, avg, nEdge, 1, 0, overlap_r, nub_r, radius, &dir);
	      if (dir){
		FreeP(dir);
		dir = 0;
	      }
            }
          }
        }
      }
    }
    if(!active)
      cs->Active[cRepCyl] = false;
  }
#endif
}
