
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
#include"RepNonbondedSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Setting.h"
#include"main.h"
#include"ShaderMgr.h"
#include"Scene.h"
#include"CGO.h"

struct RepNonbondedSphere : Rep {
  using Rep::Rep;

  ~RepNonbondedSphere() override;

  cRep_t type() const override { return cRepNonbondedSphere; }
  void render(RenderInfo* info) override;

  CGO *shaderCGO, *primitiveCGO;
};

#include"ObjectMolecule.h"

RepNonbondedSphere::~RepNonbondedSphere()
{
  CGOFree(shaderCGO);
  CGOFree(primitiveCGO);
}

void RepNonbondedSphere::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;

  if(ray) {
#ifndef _PYMOL_NO_RAY
    CGORenderRay(I->primitiveCGO, ray, info, NULL, NULL, I->cs->Setting.get(), I->obj->Setting.get());
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      if (I->shaderCGO){
	CGORenderGLPicking(I->shaderCGO, info, &I->context, I->cs->Setting.get(), I->obj->Setting.get());
      } else if (I->primitiveCGO){
	CGORenderGLPicking(I->primitiveCGO, info, &I->context, I->cs->Setting.get(), I->obj->Setting.get());
      }
    } else { /* rendering */
      short use_shader, use_sphere_shader;
      use_shader = SettingGetGlobal_i(G, cSetting_nb_spheres_use_shader) &&
                   SettingGetGlobal_b(G, cSetting_use_shaders);
      use_sphere_shader = (SettingGetGlobal_i(G, cSetting_nb_spheres_use_shader)==1) &&
	                   SettingGetGlobal_b(G, cSetting_use_shaders);

      if (I->shaderCGO){
	if (!use_shader || use_sphere_shader ^ I->shaderCGO->has_draw_sphere_buffers){
	  CGOFree(I->shaderCGO);
	  I->shaderCGO = 0;
	}
      }
      if (use_shader){
	if (!I->shaderCGO){
          if (use_sphere_shader){
            I->shaderCGO = CGOOptimizeSpheresToVBONonIndexed(I->primitiveCGO, 0, true);
          } else {
            ok_assert(1, I->shaderCGO = CGOSimplify(I->primitiveCGO, 0, SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_nb_spheres_quality)));
            ok_assert(1, CGOOptimizeToVBONotIndexed(&I->shaderCGO));
          }
          I->shaderCGO->use_shader = true;
        }
        CGORenderGL(I->shaderCGO, NULL, I->cs->Setting.get(), I->obj->Setting.get(), info, I);
      } else {
        CGORenderGL(I->primitiveCGO, NULL, I->cs->Setting.get(), I->obj->Setting.get(), info, I);
      }
    }
  }
  return;
ok_except1:
  CGOFree(I->shaderCGO);
  I->invalidate(cRepInvPurge);
  I->cs->Active[cRepNonbondedSphere] = false;
}

Rep *RepNonbondedSphereNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj = cs->Obj;

  unsigned char *active = NULL;
  int nSphere = 0;

  float prev_transp = -1;
  float const transp = SettingGet<float>(G, cs->Setting.get(),
      obj->Setting.get(), cSetting_nonbonded_transparency);

  int ok = true;

  if (ok)
    active = pymol::malloc<unsigned char>(cs->NIndex);
  CHECKOK(ok, active);

  if((obj->RepVisCache & cRepNonbondedSphereBit)){
    for(int a = 0; a < cs->NIndex; a++) {
      AtomInfoType *ai = obj->AtomInfo + cs->IdxToAtm[a];
      active[a] = (!ai->bonded && (ai->visRep & cRepNonbondedSphereBit));
      if (active[a])
        nSphere++;
    }
  }
  if(!nSphere) {
    FreeP(active);
    return (NULL);
  }
  float nb_spheres_size =
    SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_nb_spheres_size);

  auto I = new RepNonbondedSphere(cs, state);
  I->shaderCGO = NULL;
  I->primitiveCGO = NULL;

  /* Generate primitiveCGO */
  int NP = 0;
  I->primitiveCGO = CGONew(G);
  for(int a = 0; ok && a < cs->NIndex; a++){
    if(active[a]) {
      int a1 = cs->IdxToAtm[a];
      AtomInfoType *ai = obj->AtomInfo + a1;
      NP++;
      const float* v1 = cs->coordPtr(a);
      int c1 = ai->color;
      const float *vc;
      if(ColorCheckRamped(G, c1)) {
        float tmpColor[3];
        ColorGetRamped(G, c1, v1, tmpColor, state);
        vc = tmpColor;
      } else {
        vc = ColorGet(G, c1);
      }
      CGOPickColor(I->primitiveCGO, a1, (ai->masked ? cPickableNoPick : cPickableAtom));

      auto const at_transp =
          AtomSettingGetWD(G, ai, cSetting_nonbonded_transparency, transp);

      if (prev_transp != at_transp) {
        CGOAlpha(I->primitiveCGO, 1.f - at_transp);

        if (at_transp > 0) {
          I->setHasTransparency();
        }
      }
      CGOColorv(I->primitiveCGO, vc);
      CGOSphere(I->primitiveCGO, v1, nb_spheres_size);
    }
    ok &= !G->Interrupt;
  }
  CGOStop(I->primitiveCGO);
  I->primitiveCGO->sphere_quality = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_nb_spheres_quality);
  FreeP(active);
  if (!ok){
    delete I;
    I = NULL;
  }
  return (Rep *) I;
}
