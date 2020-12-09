/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#include "RepSphere.h"
#include "RepSphereGenerate.h"
#include "CGO.h"
#include "Feedback.h"
#include "ShaderMgr.h"
#include "Err.h"

void RepSphere_Generate_Triangles(PyMOLGlobals *G, RepSphere *I,
                                  RenderInfo *info) {
  short use_shader;
  int ok = true;
  int sphere_quality = SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(),
                                    cSetting_sphere_quality);

  use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) &&
               SettingGetGlobal_b(G, cSetting_use_shaders);

  // generate the CGO
  if (use_shader) {
    CGO *convertcgo = CGOSimplify(I->primitiveCGO, 0, sphere_quality);
    CHECKOK(ok, convertcgo);
    if (ok){
      I->renderCGO = CGOOptimizeToVBONotIndexed(convertcgo, 0);
      assert(I->renderCGO->use_shader);
    }
    CGOFree(convertcgo);
  } else {
    I->renderCGO = I->primitiveCGO;
  }
  CHECKOK(ok, I->renderCGO);

  if (!ok) {
    CGOFree(I->renderCGO);
    I->invalidate(cRepInvPurge);
    I->cs->Active[cRepSphere] = false;
  } else {
    I->renderCGO->sphere_quality = sphere_quality;
  }
}

void RepSphere_Generate_Impostor_Spheres(PyMOLGlobals *G, RepSphere *I,
                                         RenderInfo *info) {
  if (!I->renderCGO) {
    CGO *convertcgo = NULL;
    convertcgo = CGOOptimizeSpheresToVBONonIndexed(I->primitiveCGO, 0, true);
    if (convertcgo) {
      I->renderCGO = convertcgo;
      I->renderCGO->use_shader = true;
    }
  }
}

/* simple, default point width points -- modes 1 or 6 */
void RepSphere_Generate_Point_Sprites(PyMOLGlobals *G, RepSphere *I,
                                      RenderInfo *info, int sphere_mode) {
  short use_shader;
  use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) &
               SettingGetGlobal_b(G, cSetting_use_shaders);

  CGO *pointCGO = CGOConvertSpheresToPoints(I->primitiveCGO);
  // generate the CGO
  if (use_shader) {
    I->renderCGO = CGOOptimizeToVBONotIndexed(pointCGO, 0);

    CGO *newcgo = CGONew(G);

    CGOSpecialWithArg(newcgo, SPHERE_MODE_OPS, sphere_mode);
    CGOAppendNoStop(newcgo, I->renderCGO);
    CGOSpecialWithArg(newcgo, SPHERE_MODE_OPS, -sphere_mode);

    CGOStop(newcgo);
    CGOFreeWithoutVBOs(I->renderCGO);
    I->renderCGO = newcgo;
    I->renderCGO->use_shader = true;
    CGOFree(pointCGO);
  } else {
    CGO *newcgo = CGONew(G);

    CGOSpecialWithArg(newcgo, SPHERE_MODE_OPS, sphere_mode);
    CGOAppendNoStop(newcgo, pointCGO);
    CGOSpecialWithArg(newcgo, SPHERE_MODE_OPS, -sphere_mode);

    CGOStop(newcgo);
    I->renderCGO = newcgo;
    CGOFree(pointCGO);
  }
  return;
}
