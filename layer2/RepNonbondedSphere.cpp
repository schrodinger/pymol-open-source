
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

typedef struct RepNonbondedSphere {
  Rep R;
  CGO *shaderCGO, *primitiveCGO;
} RepNonbondedSphere;

#include"ObjectMolecule.h"

static
void RepNonbondedSphereFree(RepNonbondedSphere * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  if (I->primitiveCGO){
    CGOFree(I->primitiveCGO);
    I->primitiveCGO = 0;
  }
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepNonbondedSphereRender(RepNonbondedSphere * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  auto pick = info->pick;
  PyMOLGlobals *G = I->R.G;

  if(ray) {
#ifndef _PYMOL_NO_RAY
    CGORenderRay(I->primitiveCGO, ray, info, NULL, NULL, I->R.cs->Setting, I->R.obj->Setting);
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      if (I->shaderCGO){
	CGORenderGLPicking(I->shaderCGO, info, &I->R.context, I->R.cs->Setting, I->R.obj->Setting);
      } else if (I->primitiveCGO){
	CGORenderGLPicking(I->primitiveCGO, info, &I->R.context, I->R.cs->Setting, I->R.obj->Setting);
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
            ok_assert(1, I->shaderCGO = CGOSimplify(I->primitiveCGO, 0, SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_nb_spheres_quality)));
            ok_assert(1, CGOCombineBeginEnd(&I->shaderCGO));
            ok_assert(1, CGOOptimizeToVBONotIndexed(&I->shaderCGO));
          }
          I->shaderCGO->use_shader = true;
        }
        CGORenderGL(I->shaderCGO, NULL, I->R.cs->Setting, I->R.obj->Setting, info, &I->R);
      } else {
        CGORenderGL(I->primitiveCGO, NULL, I->R.cs->Setting, I->R.obj->Setting, info, &I->R);
      }
    }
  }
  return;
ok_except1:
  CGOFree(I->shaderCGO);
  I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
  I->R.cs->Active[cRepNonbondedSphere] = false;
}

Rep *RepNonbondedSphereNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj = cs->Obj;

  unsigned char *active = NULL;
  int nSphere = 0;
  float transp =
    1.f - SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_nonbonded_transparency);
  int ok = true;

  OOAlloc(G, RepNonbondedSphere);
  CHECKOK(ok, I);

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
    OOFreeP(I);
    FreeP(active);
    return (NULL);
  }
  float nb_spheres_size =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_nb_spheres_size);

  RepInit(G, &I->R);
  I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepNonbondedSphereRender;
  I->R.fFree = (void (*)(struct Rep *)) RepNonbondedSphereFree;
  I->R.fRecolor = NULL;
  I->R.obj = (CObject *) (cs->Obj);
  I->R.cs = cs;
  I->shaderCGO = NULL;
  I->primitiveCGO = NULL;

  /* Generate primitiveCGO */
  float at_transp;
  int NP = 0;
  bool alpha_set = false;
  I->primitiveCGO = CGONew(G);
  CGOAlpha(I->primitiveCGO, transp);
  for(int a = 0; ok && a < cs->NIndex; a++){
    if(active[a]) {
      int a1 = cs->IdxToAtm[a];
      AtomInfoType *ai = obj->AtomInfo + a1;
      NP++;
      const float *v1 = cs->Coord + 3 * a;
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
      if(AtomSettingGetIfDefined(G, ai, cSetting_nonbonded_transparency, &at_transp)){
        CGOAlpha(I->primitiveCGO, 1.f - at_transp);
        alpha_set = true;
      } else if (alpha_set){
        /* if atom level transparency is set, and this atom doesn't have transparency,
           then set back to object-level transparency */
        CGOAlpha(I->primitiveCGO, transp);
        alpha_set = false;
      }
      CGOColorv(I->primitiveCGO, vc);
      CGOSphere(I->primitiveCGO, v1, nb_spheres_size);
    }
    ok &= !G->Interrupt;
  }
  CGOStop(I->primitiveCGO);
  I->primitiveCGO->sphere_quality = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_nb_spheres_quality);
  if (ok){
    I->R.context.object = (void *) obj;
    I->R.context.state = state;
  }
  FreeP(active);
  if (!ok){
    RepNonbondedSphereFree(I);
    I = NULL;
  }
  return (Rep *) I;
}
