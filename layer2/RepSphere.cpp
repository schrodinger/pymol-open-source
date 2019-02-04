
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
#include"RepSphere.h"
#include"RepSphereImmediate.h"
#include"RepSphereGenerate.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"Util.h"
#include"Feedback.h"
#include "ShaderMgr.h"
#include "Scene.h"
#include"CGO.h"
#include"ObjectMolecule.h"
#include "Lex.h"

#define SPHERE_NORMAL_RANGE 6.f
#define SPHERE_NORMAL_RANGE2 (SPHERE_NORMAL_RANGE*SPHERE_NORMAL_RANGE)

/* defer_builds_mode = 5 : Immediate mode for any sphere_mode

   sphere_mode :

   0) Geometry shaders (quality based on sphere_quality, default 1)
   1) rectangular points with the same size, that can be changed with sphere_point_size
   2) rectangles with constant size relative to vdw and scene scale (i.e., changes when zoomed)
      maxed by a multiple of sphere_point_max_size (set it below 1 to see it influence, 3*pixel_scale max)
   3) same as 2 but with circles
   4) no longer available
   5) Uses the fast ARB Shader that approximates spheres as circles
   6-8) same as 1-3 but with normals computed from close atoms to mimic nice lighting
   9) GLSL Shader Spheres (only when shaders are available)

 */

static
void RepSphereFree(RepSphere * I)
{
  if (I->primitiveCGO == I->renderCGO) {
    I->primitiveCGO = 0;
  }
  CGOFree(I->primitiveCGO);
  CGOFree(I->renderCGO);
  CGOFree(I->spheroidCGO);
  FreeP(I->LastColor);
  FreeP(I->LastVisib);
  RepPurge(&I->R);
  OOFreeP(I);
}

/* MULTI-INSTSANCE TODO:  isn't this a conflict? */
CShaderPrg *sphereARBShaderPrg = NULL;

void RenderSphereComputeFog(PyMOLGlobals *G, RenderInfo *info, float *fog_info)
{
  /* compute -Ze = (Wc) of fog start */
  float nv[4];
  nv[3] =
    (info->front +
     (info->back - info->front) * SettingGetGlobal_f(G, cSetting_fog_start));
  /* compute Zc of fog start using std. perspective transformation */
  nv[2] =
    (nv[3] * (info->back + info->front) -
     2 * (info->back * info->front)) / (info->back - info->front);
  /* compute Zc/Wc to get normalized depth coordinate of fog start */
  nv[0] = (nv[2] / nv[3]);
  fog_info[0] = (nv[0] * 0.5) + 0.5;
  
  fog_info[1] = 1.0F / (1.0 - fog_info[0]);     /* effective range of fog */
  
}

#ifndef _PYMOL_NO_RAY
static int RepSphereRenderRay(PyMOLGlobals *G, RepSphere * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  float alpha = 1.0F - 
    SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_sphere_transparency);
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;
  ray->transparentf(1.0 - alpha);
  if (I->spheroidCGO){
    CGORenderRay(I->spheroidCGO, ray, info, NULL, NULL, I->R.cs->Setting, I->R.obj->Setting);
  } else {
    CGORenderRay(I->primitiveCGO, ray, info, NULL, NULL, I->R.cs->Setting, I->R.obj->Setting);
  }
  ray->transparentf(0.0);
  return true;
}
#endif

static void RepSphereRenderPick(RepSphere * I, RenderInfo * info, int sphere_mode)
{
  PyMOLGlobals *G = I->R.G;

  if (!I->renderCGO){
    // only for sphere_mode 5, where we don't use a renderCGO (yet) ARB: immediate mode GL_QUADS
    short use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) &&
      SettingGetGlobal_b(G, cSetting_use_shaders);
    CGO *convertcgo = CGOSimplify(I->primitiveCGO, 0, 0);
    CGO *convertcgo2 = CGOCombineBeginEnd(convertcgo, 0);
    if (use_shader){
      I->renderCGO = CGOOptimizeToVBONotIndexed(convertcgo2, 0);
      CGOFree(convertcgo2);
    } else {
      I->renderCGO = convertcgo2;
    }
    I->renderCGO->use_shader = use_shader;
    CGOFree(convertcgo);
  }
  CGORenderGLPicking(I->renderCGO, info, &I->R.context, I->R.cs->Setting, I->R.obj->Setting);
}

static int RepGetSphereMode(PyMOLGlobals *G, RepSphere * I, bool use_shader){
  int sphere_mode = SettingGet_i(G, I->R.cs->Setting,
				 I->R.obj->Setting,
				 cSetting_sphere_mode);
  if (sphere_mode == 4) // sphere_mode 4 no longer exists, use default
    sphere_mode = -1;
    switch (sphere_mode) {
    case 5:
#ifdef _PYMOL_ARB_SHADERS
      if (!sphereARBShaderPrg && G->HaveGUI && G->ValidContext) {
      sphereARBShaderPrg = CShaderPrg::NewARB(G, "sphere_arb",
                                              G->ShaderMgr->GetShaderSource("sphere_arb_vs.vs"),
                                              G->ShaderMgr->GetShaderSource("sphere_arb_fs.fs"));
      }
      if (!sphereARBShaderPrg)
#endif
      {
        PRINTFB(G, FB_ShaderMgr, FB_Warnings)
          " Warning: ARB shaders (sphere_mode=5) not supported.\n" ENDFB(G);
        if (!use_shader || !G->ShaderMgr->ShaderPrgExists("sphere")) {
        sphere_mode = 9;
        } else {
          sphere_mode = 0;
        }
      }
      break;
    case -1:
      sphere_mode = 9;
    case 9:
    if (!use_shader || !G->ShaderMgr->ShaderPrgExists("sphere")) {
        sphere_mode = 0;
      }
    }
  return sphere_mode;
}

static void RepSphereRender(RepSphere * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  auto pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  int ok = true;
  bool use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) &&
                    SettingGetGlobal_b(G, cSetting_use_shaders);
  if(ray) {
#ifndef _PYMOL_NO_RAY
    ok &= RepSphereRenderRay(G, I, info);
#endif
    return;
  }
  int sphere_mode = RepGetSphereMode(G, I, use_shader);
  if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      RepSphereRenderPick(I, info, sphere_mode);
    } else {                    /* not pick, render! */
#ifdef _PYMOL_ARB_SHADERS
      if (sphere_mode == 5){
        // we need to check sphere_mode 5 here until
        // we implement the CGO ARB operation, since 
        // picking uses the I->renderCGO
        RepSphere_Generate_ARB_Spheres(G, I, info);
        return; // sphere_mode 5 does not use I->renderCGO yet
      }
#endif
      if (I->renderCGO){
        if (I->renderCGO->use_shader != use_shader){
          CGOFree(I->renderCGO);
          I->renderCGO = 0;
        } else {
          CGORenderGL(I->renderCGO, NULL, NULL, NULL, info, &I->R);
          return;
        }
      }
      // Generate renderCGO for sphere_mode
      switch (sphere_mode) {
      case 0:              /* memory-efficient sphere rendering */
        RepSphere_Generate_Triangles(G, I, info);
        break;
      case 9: // use GLSL impostor shader
        RepSphere_Generate_Impostor_Spheres(G, I, info);
        break;
      default:
        // sphere_modes 1,2,3,6,7,8
        RepSphere_Generate_Point_Sprites(G, I, info, sphere_mode);
        break;
      }

      CHECKOK(ok, I->renderCGO);
      if (!ok){
        CGOFree(I->renderCGO);
        I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
        I->R.cs->Active[cRepSphere] = false;
      }

      if (I->renderCGO)
        CGORenderGL(I->renderCGO, NULL, NULL, NULL, info, &I->R);
    }
  }
}

static
int RepSphereSameVis(RepSphere * I, CoordSet * cs)
{
  bool *lv;
  int *lc;
  int a;
  AtomInfoType *ai;
  if(I->LastVisib && I->LastColor) {
    lv = I->LastVisib;
    lc = I->LastColor;

    for(a = 0; a < cs->NIndex; a++) {
      ai = cs->getAtomInfo(a);
      if(*(lv++) != GET_BIT(ai->visRep, cRepSphere)) {
        return false;
      }
      if(*(lc++) != ai->color) {
        return false;
      }
    }
  } else {
    return false;
  }
  return true;
}

static bool RepSphereDetermineAtomVisibility(PyMOLGlobals *G,
    AtomInfoType *ati1, int cartoon_side_chain_helper, int ribbon_side_chain_helper)
{
  if (!(ati1->flags & cAtomFlag_polymer))
    return true;

  bool sc_helper =
    (GET_BIT(ati1->visRep, cRepCartoon) &&
     AtomSettingGetWD(G, ati1, cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper)) ||
    (GET_BIT(ati1->visRep, cRepRibbon) &&
     AtomSettingGetWD(G, ati1, cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper));

  if (sc_helper) {
    int prot1 = ati1->protons;

    if(prot1 == cAN_N) {
      if(ati1->name == G->lex_const.N) {
        if(ati1->resn != G->lex_const.PRO)
	  return false;
      }
    } else if(prot1 == cAN_O) {
      if(ati1->name == G->lex_const.O)
	return false;
    } else if(prot1 == cAN_C) {
      if(ati1->name == G->lex_const.C)
	return false;
    }
  }
  return true;
}

static void RepSphereAddAtomVisInfoToStoredVC(RepSphere *I, ObjectMolecule *obj,
    CoordSet * cs, int state, int a1, AtomInfoType *ati1, int a,
    float sphere_scale, int sphere_color, float transp,
    int *variable_alpha, float sphere_add)
{
  PyMOLGlobals *G = cs->State.G;
  float at_transp = transp;
  int c1;
  float *v0, vc[3];
  const float *vcptr;

  float at_sphere_scale = AtomSettingGetWD(G, ati1, cSetting_sphere_scale, sphere_scale);
  int at_sphere_color = AtomSettingGetWD(G, ati1, cSetting_sphere_color, sphere_color);

  if(AtomSettingGetIfDefined(G, ati1, cSetting_sphere_transparency, &at_transp))
    *variable_alpha = true;
  
  CGOPickColor(I->primitiveCGO, a1, ati1->masked ? cPickableNoPick : cPickableAtom);
  
  if(at_sphere_color == -1)
    c1 = ati1->color;
  else
    c1 = at_sphere_color;
  v0 = cs->Coord + 3 * a;

  if(ColorCheckRamped(G, c1)) {
    ColorGetRamped(G, c1, v0, vc, state);
    vcptr = vc;
  } else {
    vcptr = ColorGet(G, c1);   /* save new color */
  }
  float alpha = 1.0F - at_transp;
  CGOAlpha(I->primitiveCGO, alpha);
  CGOColorv(I->primitiveCGO, vcptr);
  float radius = obj->AtomInfo[a1].vdw * at_sphere_scale + sphere_add;
  CGOSphere(I->primitiveCGO, v0, radius);
}

/* This function is extraneous to do every time 
   we need to compute normals for RepSphere */
static float SphereComputeCutMultiplier(SphereRec *sr){
  int a;
  float *dot = sr->dot[0];
  int n_dot = sr->nDot;
	float cut_mult = -1.0F;
	  for(a = 1; a < n_dot; a++) {
	    float t_dot = dot_product3f(dot, dot + a * 3);
	    if(cut_mult < t_dot)
	      cut_mult = t_dot;
	  }
  return cut_mult;
}

/*
 * for the spheroid implementation, this function sets color
 * and pickcolor in the CGO given the atom idx
 *
 */
static void RepSphereCGOSetSphereColorAndPick(ObjectMolecule *obj, CoordSet * cs, CGO *cgo, int idx, int state, float transp, int sphere_color){
  PyMOLGlobals *G = obj->Obj.G;
  int a1 = cs->IdxToAtm[idx];
  float at_transp;
  AtomInfoType *ati1 = obj->AtomInfo + a1;

  int at_sphere_color = AtomSettingGetWD(G, ati1, cSetting_sphere_color, sphere_color);

  if(AtomSettingGetIfDefined(G, ati1, cSetting_sphere_transparency, &at_transp)) {
    float alpha = 1.0F - at_transp;
    CGOAlpha(cgo, alpha);
  }
  int c1;
  if(at_sphere_color == -1)
    c1 = ati1->color;
  else
    c1 = at_sphere_color;
  if(ColorCheckRamped(G, c1)) {
    float *v0 = cs->Coord + 3 * idx;
    float color[3];
    ColorGetRamped(G, c1, v0, color, state);
    CGOColorv(cgo, color);
  } else {
    const float *color = ColorGet(G, c1);   /* save new color */
    CGOColorv(cgo, color);
  }

}

static
CGO *RepSphereGeneratespheroidCGO(ObjectMolecule * I, CoordSet *cs, SphereRec *sp, int state){
  int idx, a, b, c;
  int *q, *s;
  bool ok = true;
  float spheroid_scale =
    SettingGet_f(I->Obj.G, cs->Setting, I->Obj.Setting, cSetting_spheroid_scale);
  int sphere_color =
    SettingGet_color(I->Obj.G, cs->Setting, I->Obj.Setting, cSetting_sphere_color);
  float transp = SettingGet_f(I->Obj.G, cs->Setting, I->Obj.Setting, cSetting_sphere_transparency);

  CGO *cgo = CGONew(I->Obj.G);
  for(idx = 0; idx < cs->NIndex; idx++) {
    float *v0 = &cs->Coord[3 * idx];
    a = cs->IdxToAtm[idx];
    q = sp->Sequence;
    s = sp->StripLen;
    RepSphereCGOSetSphereColorAndPick(I, cs, cgo, idx, state, transp, sphere_color);
    for(b = 0; ok && b < sp->NStrip; b++) {
      float *sphLen = cs->Spheroid + (sp->nDot * a);
      float *sphNorm = cs->SpheroidNormal + (3 * sp->nDot * a);
      CGOBegin(cgo, GL_TRIANGLE_STRIP);
      for(c = 0; c < (*s); c++) {
        float sphTmp, *sphTmpN = sphNorm + 3 * (*q);
        CGONormalv(cgo, sphTmpN);
        sphTmp = (*(sphLen + (*q))) * spheroid_scale;
        // point
        CGOVertex(cgo, v0[0] + sphTmp * sp->dot[*q][0],
                       v0[1] + sphTmp * sp->dot[*q][1],
                       v0[2] + sphTmp * sp->dot[*q][2]);
        q++;
      }
      CGOEnd(cgo);
      s++;
      ok &= !I->Obj.G->Interrupt;
    }
  }
  CGOStop(cgo);
  if (!ok){
    CGOFree(cgo);
  }
  return cgo;
}

/*
 * when normals are needed (sphere_mode 6-8), 
 * this function computes the normal for an atom
 * and pickcolor in the CGO given the atom idx
 *
 */
static void RepSphereSetNormalForSphere(RepSphere *I, MapType *map, float *v_tmp, 
                                        float *v, float cut_mult, int a,
                                        int *active, float *dot, int n_dot){
	  int h, k, l, b, i, j;
	  float v1[3];
  int n_dot_active, *da;
  float *vv;
  float dst;

	  MapLocus(map, v, &h, &k, &l);
	  da = active;
	  for(b = 0; b < n_dot; b++) {
	    *(da++) = b * 3;
	  }
	  n_dot_active = n_dot;
	  i = *(MapEStart(map, h, k, l));
	  if(i) {
	    j = map->EList[i++];
    while(j >= 0) {
	      if(j != a) {
		vv = v_tmp + 3 * j;
        if(within3fret(vv, v, SPHERE_NORMAL_RANGE, SPHERE_NORMAL_RANGE2, v1, &dst)) {
		  float cutoff = dst * cut_mult;
		  b = 0;
		  while(b < n_dot_active) {
		    vv = dot + active[b];
		    if(dot_product3f(v1, vv) > cutoff) {
		      n_dot_active--;
		      active[b] = active[n_dot_active];
		    }
		    b++;
		  }
		}
	      }
	      j = map->EList[i++];
	    }
	  }
  float v0[3];
	    if(!n_dot_active) {
	      v0[0] = 0.0F;
	      v0[1] = 0.0F;
	      v0[2] = 1.0F;
	    } else {
	      zero3f(v0);
	      b = 0;
	      while(b < n_dot_active) {
		vv = dot + active[b];
		add3f(vv, v0, v0);
		b++;
	      }
	      normalize3f(v0);
    CGONormalv(I->primitiveCGO, v0);
	    }
}

Rep *RepSphereNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int ok = true;
  int a, a1;
  bool *lv;
  int *lc;
  float sphere_scale, sphere_add = 0.f;
  int sphere_color;
  int cartoon_side_chain_helper = 0;
  int ribbon_side_chain_helper = 0;
  AtomInfoType *ati1;
  int sphere_mode = 0;
  bool *marked = NULL;
  float transp;
  int variable_alpha = false;
  short use_shader = SettingGetGlobal_b(G, cSetting_sphere_use_shader) &&
                     SettingGetGlobal_b(G, cSetting_use_shaders);
  // skip if not visible
  if(!cs->hasRep(cRepSphereBit))
    return NULL;

  OOCalloc(G, RepSphere);
  CHECKOK(ok, I);
  if (!ok)
    return NULL;
  obj = cs->Obj;

  marked = pymol::calloc<bool>(obj->NAtom);
  CHECKOK(ok, marked);
  if (ok)
    RepInit(G, &I->R);
  I->renderCGO = NULL;
  I->primitiveCGO = NULL;
  if (cs->Spheroid)
    I->spheroidCGO = RepSphereGeneratespheroidCGO(obj, cs, G->Sphere->Sphere[1], state);

  if (ok){
    sphere_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_mode);
    if (!use_shader && (sphere_mode == 5 || sphere_mode == 9)){
      sphere_mode = 0;
    }
  }
  if (ok){
    sphere_color =
      SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_color);
    cartoon_side_chain_helper =
      SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_side_chain_helper);
    ribbon_side_chain_helper =
      SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_side_chain_helper);
    transp = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_transparency);
    sphere_scale = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_scale);
  }

  if (ok){
    I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepSphereRender;
    I->R.fFree = (void (*)(struct Rep *)) RepSphereFree;
    I->R.fSameVis = (int (*)(struct Rep *, struct CoordSet *)) RepSphereSameVis;
    I->R.obj = (CObject *) obj;
    I->R.cs = cs;
    I->R.context.object = (void *) obj;
    I->R.context.state = state;
  }
  /* raytracing primitives */

  if (ok){
    if(SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_sphere_solvent)) { /* are we generating a solvent surface? */
      sphere_add = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_solvent_radius);       /* if so, get solvent radius */
    }
    
    if(SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_pickable)) {
      I->R.P = pymol::malloc<Pickable>(cs->NIndex + 1);
      CHECKOK(ok, I->R.P);
    }
  }
  I->primitiveCGO = CGONew(G);

  bool needNormals = (sphere_mode >= 6) && (sphere_mode < 9);
  int nspheres = 0;
  if (needNormals){
    float *v_tmp = VLAlloc(float, 1024);
  for(a = 0; ok && a < cs->NIndex; a++) {
    a1 = cs->IdxToAtm[a];
    ati1 = obj->AtomInfo + a1;
    /* store temporary visibility information */
    marked[a1] = GET_BIT(ati1->visRep,cRepSphere) &&
        RepSphereDetermineAtomVisibility(G, ati1, 
                                         cartoon_side_chain_helper, ribbon_side_chain_helper);
    if(marked[a1]) {
        int cnc = nspheres * 3;
        nspheres++;
        VLACheck(v_tmp, float, cnc + 3);
        copy3f(cs->Coord + 3 * a, &v_tmp[cnc]);
    }
    ok &= !G->Interrupt;
  }
    MapType *map = MapNew(G, SPHERE_NORMAL_RANGE, v_tmp, nspheres, NULL);
    float cut_mult = SphereComputeCutMultiplier(G->Sphere->Sphere[1]);
    float *dot = G->Sphere->Sphere[1]->dot[0];
    int n_dot = G->Sphere->Sphere[1]->nDot;
    int *active = pymol::malloc<int>(2 * n_dot);

	ok &= MapSetupExpress(map);
    for(a = 0; ok && a < cs->NIndex; a++) {
      a1 = cs->IdxToAtm[a];
      if(marked[a1]) {
      ati1 = obj->AtomInfo + a1;
        RepSphereSetNormalForSphere(I, map, v_tmp, &v_tmp[a * 3], cut_mult, a, active, dot, n_dot);
        RepSphereAddAtomVisInfoToStoredVC(I, obj, cs, state, a1, ati1, a, sphere_scale, sphere_color, transp, &variable_alpha, sphere_add);
      }
	ok &= !G->Interrupt;
      }
    FreeP(active);
    VLAFreeP(v_tmp);
    MapFree(map);
  } else { // no normals necessary
    for(a = 0; ok && a < cs->NIndex; a++) {
      a1 = cs->IdxToAtm[a];
      ati1 = obj->AtomInfo + a1;
      /* store temporary visibility information */
      marked[a1] = GET_BIT(ati1->visRep,cRepSphere) && 
        RepSphereDetermineAtomVisibility(G, ati1, 
                                         cartoon_side_chain_helper, ribbon_side_chain_helper);
      if(marked[a1]) {
        nspheres++;
        RepSphereAddAtomVisInfoToStoredVC(I, obj, cs, state, a1, ati1, a, sphere_scale, sphere_color, transp, &variable_alpha, sphere_add);
      }
      ok &= !G->Interrupt;
    }
  }
  CGOStop(I->primitiveCGO);

  if(ok) {
    if(!I->LastVisib)
      I->LastVisib = pymol::malloc<bool>(cs->NIndex);
    CHECKOK(ok, I->LastVisib);
    if(ok && !I->LastColor)
      I->LastColor = pymol::malloc<int>(cs->NIndex);
    CHECKOK(ok, I->LastColor);
    if (ok){
      lv = I->LastVisib;
      lc = I->LastColor;
      obj = cs->Obj;
      AtomInfoType *ai2 = obj->AtomInfo;
      if(sphere_color == -1){
	for(a = 0; a < cs->NIndex; a++) {
          int at = cs->IdxToAtm[a];
	  *(lv++) = marked[at];
	  *(lc++) = (ai2 + at)->color;
	}
      } else {
	for(a = 0; a < cs->NIndex; a++) {
	  *(lv++) = marked[cs->IdxToAtm[a]];
	  *(lc++) = sphere_color;
	}
      }
    }
  }

  FreeP(marked);
  if(nspheres == 0 || !ok) {
    RepSphereFree(I);
    I = NULL;
  }
  return (Rep *) I;
}
