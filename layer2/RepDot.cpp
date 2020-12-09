
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
#include"RepDot.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"ObjectMolecule.h"
#include"Scene.h"
#include"ShaderMgr.h"
#include"CGO.h"

#include "pymol/algorithm.h"

RepDot::~RepDot()
{
  auto I = this;
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->T);
  FreeP(I->F);
  FreeP(I->VN);
  FreeP(I->A);
  FreeP(I->Atom);
}

static int RepDotCGOGenerate(RepDot * I)
{
  PyMOLGlobals *G = I->G;
  float *v = I->V;
  int c = I->N;
  int cc = 0;
  int ok = true;
  CGO *cgo = NULL;

  int normals =
    SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_dot_normals);
  bool const dot_as_spheres = SettingGet<bool>(
      G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_dot_as_spheres);

  cgo = CGONew(G);
  CHECKOK(ok, cgo);
  if (dot_as_spheres){
    while(ok && c--) {
      if(!cc) {             /* load up the current vertex color */
	cc = (int) (*(v++));
	ok &= CGOColorv(cgo, v);
	v += 3;
      }
      if(ok && normals)
	ok &= CGONormalv(cgo, v);
      v += 3;
      if (ok)
	ok &= CGOSphere(cgo, v, 1.f);
      v += 3;
      cc--;
    }
  } else {
    if (ok)
      ok &= CGOBegin(cgo, GL_POINTS);
    while(ok && c--) {
      if(!cc) {             /* load up the current vertex color */
	cc = (int) (*(v++));
	ok &= CGOColorv(cgo, v);
	v += 3;
      }
      if(normals)
	CGONormalv(cgo, v);
      v += 3;
      if (ok)
	ok &= CGOVertexv(cgo, v);
      v += 3;
      cc--;
    }
    if (ok)
      ok &= CGOEnd(cgo);
  }
  if (ok)
    ok &= CGOStop(cgo);
  if (ok) {
    if (dot_as_spheres){
      CGO *tmpCGO = CGONew(G);
      if (ok) ok &= CGOEnable(tmpCGO, GL_SPHERE_SHADER);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DOT_LIGHTING);
      if (ok) ok &= CGOSpecial(tmpCGO, DOT_WIDTH_FOR_DOT_SPHERES);
      if (ok) {
        tmpCGO->free_append(CGOOptimizeSpheresToVBONonIndexedNoShader(cgo,
            CGO_BOUNDING_BOX_SZ + fsizeof<cgo::draw::sphere_buffers>() + 2));
      }
      if (ok) ok &= CGODisable(tmpCGO, GL_SPHERE_SHADER);
      if (ok) ok &= CGOStop(tmpCGO);
      I->shaderCGO = tmpCGO;
    } else {
      CGO *tmpCGO = CGONew(G);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DEFAULT_SHADER);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DOT_LIGHTING);
      if (ok) ok &= CGOSpecial(tmpCGO, DOT_WIDTH_FOR_DOTS);
      if (ok) {
        tmpCGO->free_append(CGOOptimizeToVBONotIndexedNoShader(cgo));
      }
      if (ok) ok &= CGODisable(tmpCGO, GL_DEFAULT_SHADER);
      if (ok) ok &= CGOStop(tmpCGO);
      I->shaderCGO = tmpCGO;
    }
  }
  if (ok){
    I->shaderCGO->use_shader = true;
    I->shaderCGO_as_spheres = dot_as_spheres;
  }
  CGOFree(cgo);

  return ok;

}

void RepDot::render(RenderInfo * info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  float *v = I->V;
  int c = I->N;
  int cc = 0;
  int ok = true;

  if(ray) {
#ifndef _PYMOL_NO_RAY
    float radius;

    if(I->dotSize <= 0.0F) {
      radius = ray->PixelRadius * I->Width / 1.4142F;
    } else {
      radius = I->dotSize;
    }

    while(ok && c--) {
      if(!cc) {                 /* load up the current vertex color */
        cc = (int) (*(v++));
        ray->color3fv(v);
        v += 3;
      }
      v += 3;
      ok &= ray->sphere3fv(v, radius);
      v += 3;
      cc--;
    }
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else { /* else not pick, i.e., when rendering */
      int normals =
        SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_dot_normals);
      bool const dot_as_spheres = SettingGet<bool>(
          G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_as_spheres);

      bool const use_shader = SettingGet<bool>(G, cSetting_dot_use_shader) &&
                              SettingGet<bool>(G, cSetting_use_shaders);

      if (I->shaderCGO && ((!use_shader || CGOCheckWhetherToFree(G, I->shaderCGO)) ||
			   I->shaderCGO_as_spheres!= dot_as_spheres)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (use_shader){
	if (!I->shaderCGO){
	  ok &= RepDotCGOGenerate(I);
	}

	if (ok) {
	  const float *color;
	  color = ColorGet(G, I->obj->Color);
	  CGORenderGL(I->shaderCGO, color, NULL, NULL, info, I);
	  return; /* should not do any other rendering after shaderCGO has
		    been rendered */
	}
      } else {
	if(!normals)
	  SceneResetNormal(G, true);
        int lighting =
          SettingGet_i(G, I->cs->Setting.get(), I->obj->Setting.get(), cSetting_dot_lighting);
	if(!lighting) {
	  if(!info->line_lighting)
	    glDisable(GL_LIGHTING);
	}

	if(info->width_scale_flag)
	  glPointSize(I->Width * info->width_scale);
	else
	  glPointSize(I->Width);

        glBegin(GL_POINTS);
        while(c--) {
          if(!cc) {             /* load up the current vertex color */
            cc = (int) (*(v++));
            glColor3fv(v);
            v += 3;
          }
          if(normals)
            glNormal3fv(v);
          v += 3;
          glVertex3fv(v);
          v += 3;
          cc--;
        }
        glEnd();

        if(!lighting)
          glEnable(GL_LIGHTING);
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->invalidate(cRepInvPurge);
    I->cs->Active[cRepDot] = false;
  }
}

Rep *RepDotNew(CoordSet * cs, int state)
{
  return (RepDotDoNew(cs, cRepDotNormal, state));
}

/**
 * This routine does double duty - generating the dot representation, but also
 * acting as our surface area computation routine.
 *
 * @param mode cRepDotNormal for dot representation, or cRepDotAreaType for
 * `cmd.get_area()` computation.
 * @param state Object state for ramped coloring
 */
Rep *RepDotDoNew(CoordSet * cs, cRepDot_t mode, int state)
{
  PyMOLGlobals *G = cs->G;
  float *v, *vn;
  float *aa = NULL;
  int *tp = NULL;
  int *tf = NULL;
  float *countPtr = NULL;
  int* ati = nullptr;
  auto obj = cs->Obj;

  // skip if no dots are visible (assume all atoms "visible" for area comp.)
  if (mode != cRepDotAreaType && !cs->hasRep(cRepDotBit)) {
    return nullptr;
  }

  // are we using flags 24 & 25
  auto cullByFlag =
      SettingGet<bool>(G, cs->Setting.get(), obj->Setting.get(), cSetting_trim_dots);

  auto dot_color =
      SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_color);

  // are we ignoring hydrogens?
  auto inclH =
      SettingGet<bool>(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_hydrogens);

  float solv_rad = 0.f;
  if(SettingGet_b(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_solvent)) {    /* are we generating a solvent surface? */
    solv_rad = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_solvent_radius); /* if so, get solvent radius */
  }

  std::unique_ptr<MapType> map(
      MapNew(G, MAX_VDW + solv_rad, cs->Coord, cs->NIndex, nullptr));
  if (!map) {
    return nullptr;
  }

  // get current dot sampling
  // Note: significantly affects the accuracy of our area comp.
  auto ds = SettingGet<int>(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_density);
  SphereRec const* sp = G->Sphere->Sphere[pymol::clamp(ds, 0, 4)];

  int lastColor = cColorDefault;
  int colorCnt = 0;

  auto I = new RepDot(cs, state);

  I->dotSize = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_radius);
  I->Width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_dot_width);

  I->V = pymol::malloc<float>(cs->NIndex * sp->nDot * 10);
  ok_assert(1, I->V);
  v = I->V;

  // in area mode, we need to save additional info such as the normal vectors,
  // the partial area, the originating atom, etc.
  if (mode == cRepDotAreaType) {
    I->A = pymol::malloc<float>(cs->NIndex * sp->nDot);
    ok_assert(1, I->A);
    I->T = pymol::malloc<int>(cs->NIndex * sp->nDot);
    ok_assert(1, I->T);
    I->F = pymol::malloc<int>(cs->NIndex * sp->nDot);
    ok_assert(1, I->F);
    I->VN = pymol::malloc<float>(cs->NIndex * sp->nDot * 3);
    ok_assert(1, I->VN);
    I->Atom = pymol::malloc<int>(cs->NIndex * sp->nDot);
    ok_assert(1, I->Atom);

    aa = I->A;
    tp = I->T;
    tf = I->F;
    vn = I->VN;
    ati = I->Atom;
    inclH = true;
    cullByFlag = true;
  }

  for (int a = 0; a < cs->NIndex; ++a) {
    auto const atm = cs->IdxToAtm[a];
    auto const& ai1 = obj->AtomInfo[atm];

    if (mode != cRepDotAreaType && !(ai1.visRep & cRepDotBit)) {
      continue;
    }

    if (!inclH && ai1.isHydrogen()) {
      continue;
    }

    // If we are culling, flags control which atoms will have dot surfaces
    // generated for them.
    if (cullByFlag && (ai1.flags & (cAtomFlag_exfoliate | cAtomFlag_ignore))) {
      continue;
    }

    auto c1 = AtomSettingGetWD(G, &ai1, cSetting_dot_color, dot_color);

    if (c1 == cColorDefault) {
      c1 = ai1.color;
    }

    const float* v0 = cs->coordPtr(a);
    const float vdw = ai1.vdw + solv_rad;
    for (int b = 0; b < sp->nDot; b++) {
      const float v1[] = {
          v0[0] + vdw * sp->dot[b][0],
          v0[1] + vdw * sp->dot[b][1],
          v0[2] + vdw * sp->dot[b][2],
      };

      bool flag = true;

      for (int j : MapEIter(*map, v1, false)) {
        if (j == a) {
          continue;
        }

        auto const& ai2 = obj->AtomInfo[cs->IdxToAtm[j]];
        if (!inclH && ai2.isHydrogen()) {
          continue;
        }

        // If we are cullilng, flag 25 controls which atoms are considered
        // "present" in the surface area calculation (i.e. able to occlude
        // surface)
        if (cullByFlag && (ai2.flags & cAtomFlag_ignore)) {
          continue;
        }

        if (within3f(cs->coordPtr(j), v1, ai2.vdw + solv_rad)) {
          flag = false;
          break;
        }
      }

      if (!flag) {
        continue;
      }

      switch (mode) {
      case cRepDotNormal:
        if (c1 == lastColor && !ColorCheckRamped(G, c1)) {
          colorCnt++;
        } else {
          /* new color */
          if (countPtr)                   /* after first pass */
            *countPtr = (float) colorCnt; /* save count */
          colorCnt = 1;
          countPtr = v++;
          lastColor = c1;
          // save new color
          if (ColorCheckRamped(G, c1)) {
            ColorGetRamped(G, c1, v1, v, state);
          } else {
            copy3(ColorGet(G, c1), v);
          }
          v += 3;
        }
        *(v++) = sp->dot[b][0];
        *(v++) = sp->dot[b][1];
        *(v++) = sp->dot[b][2];
        *(v++) = v1[0];
        *(v++) = v1[1];
        *(v++) = v1[2];
        I->N++;
        break;
      case cRepDotAreaType:
        *(v++) = v1[0];
        *(v++) = v1[1];
        *(v++) = v1[2];
        *(aa++) = vdw * vdw * sp->area[b]; /* area */
        *(tp++) = ai1.customType;         /* numeric type */
        *(tf++) = ai1.flags;              /* flags */
        *(vn++) = sp->dot[b][0];
        *(vn++) = sp->dot[b][1];
        *(vn++) = sp->dot[b][2];
        *(ati++) = atm;
        I->N++;
        break;

      default:
        assert(false);
      }
    }

    ok_assert(1, !G->Interrupt);
  }

  // save count
  if (countPtr)
    *countPtr = (float) colorCnt;

  I->V = ReallocForSure(I->V, float, (v - I->V));
  ok_assert(1, I->V);

  if (mode == cRepDotAreaType) {
    I->A = ReallocForSure(I->A, float, (aa - I->A));
    ok_assert(1, I->A);
    I->T = ReallocForSure(I->T, int, (tp - I->T));
    ok_assert(1, I->T);
    I->F = ReallocForSure(I->F, int, (tf - I->F));
    ok_assert(1, I->F);
    I->VN = ReallocForSure(I->VN, float, (vn - I->VN));
    ok_assert(1, I->VN);
    I->Atom = ReallocForSure(I->Atom, int, (ati - I->Atom));
    ok_assert(1, I->Atom);
  }

  return (Rep*) I;

ok_except1:
  delete I;
  return nullptr;
}
