
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

static void RepDotRender(RepDot * I, RenderInfo * info);

static
void RepDotFree(RepDot * I)
{
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
  OOFreeP(I);
}

static int RepDotCGOGenerate(RepDot * I, RenderInfo * info)
{
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  int cc = 0;
  int ok = true;
  CGO *cgo = NULL;

  int normals =
    SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_normals);
  short dot_as_spheres = SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_as_spheres);

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
      CGO *tmpCGO = CGONew(G), *tmp2CGO = NULL;
      if (ok) ok &= CGOEnable(tmpCGO, GL_SPHERE_SHADER);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DOT_LIGHTING);
      if (ok) ok &= CGOSpecial(tmpCGO, DOT_WIDTH_FOR_DOT_SPHERES);
      tmp2CGO = CGOOptimizeSpheresToVBONonIndexedNoShader(cgo,
          CGO_BOUNDING_BOX_SZ + fsizeof<cgo::draw::sphere_buffers>() + 2);
      if (ok)
	ok &= CGOAppendNoStop(tmpCGO, tmp2CGO);
      CGOFreeWithoutVBOs(tmp2CGO);
      if (ok) ok &= CGODisable(tmpCGO, GL_SPHERE_SHADER);
      if (ok) ok &= CGOStop(tmpCGO);
      I->shaderCGO = tmpCGO;
    } else {
      CGO *convertcgo = CGOCombineBeginEnd(cgo, 0), *tmp2CGO = NULL;
      CGO *tmpCGO = CGONew(G);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DEFAULT_SHADER);
      if (ok) ok &= CGOEnable(tmpCGO, GL_DOT_LIGHTING);
      if (ok) ok &= CGOSpecial(tmpCGO, DOT_WIDTH_FOR_DOTS);
      CHECKOK(ok, convertcgo);
      if (ok)
	tmp2CGO = CGOOptimizeToVBONotIndexedNoShader(convertcgo, CGO_BOUNDING_BOX_SZ + I->N * 3 + 7);
      CHECKOK(ok, tmp2CGO);
      if (ok)
	ok &= CGOAppendNoStop(tmpCGO, tmp2CGO);
      CGOFreeWithoutVBOs(tmp2CGO);
      CGOFree(convertcgo);
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

  /* now that the shaderCGO is created, we can just call RepDotRender to render it */
  if (ok)
    RepDotRender(I, info);
      
  return ok;

}

static void RepDotRender(RepDot * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  auto pick = info->pick;
  PyMOLGlobals *G = I->R.G;
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
      short use_shader, generate_shader_cgo = 0;
      int normals =
        SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_normals);
      short dot_as_spheres = SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_as_spheres);

      use_shader = SettingGetGlobal_b(G, cSetting_dot_use_shader) & 
                   SettingGetGlobal_b(G, cSetting_use_shaders);

      if (I->shaderCGO && ((!use_shader || CGOCheckWhetherToFree(G, I->shaderCGO)) ||
			   I->shaderCGO_as_spheres!= dot_as_spheres)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (use_shader){
	if (!I->shaderCGO){
	  generate_shader_cgo = 1;
	  ok &= RepDotCGOGenerate(I, info);
	} else {
	  const float *color;
	  color = ColorGet(G, I->R.obj->Color);
	  CGORenderGL(I->shaderCGO, color, NULL, NULL, info, &I->R);
	  return; /* should not do any other rendering after shaderCGO has
		    been rendered */
	}
      }

      if (!generate_shader_cgo) {
	if(!normals)
	  SceneResetNormal(G, true);
        int lighting =
          SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_lighting);
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
    I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
    I->R.cs->Active[cRepDot] = false;
  }
}

Rep *RepDotNew(CoordSet * cs, int state)
{
  return (RepDotDoNew(cs, cRepDotNormal, state));
}

Rep *RepDotDoNew(CoordSet * cs, int mode, int state)
{

  /* this routine does double duty - generating the dot representation,
     but also acting as our surface area computation routine.
     Modes: cRepDotNormal,cRepDotAreaType
   */
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int a, b, flag, h, k, l, i, j, c1;
  float *v, *v0, vdw, *vn;
  const float *vc;
  float *aa = NULL;
  int *tp = NULL;
  int *tf = NULL;
  float *countPtr = NULL;
  int colorCnt, lastColor;
  Vector3f v1;
  MapType *map = NULL;
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds;
  float max_vdw = MAX_VDW;
  float solv_rad = 0.0;
  int inclH = true;
  int cullByFlag = false;
  int visFlag;
  int atm, *ati = NULL;
  AtomInfoType *ai1, *ai2;
  int dot_color;
  int ok = true;
  OOAlloc(G, RepDot);
  CHECKOK(ok, I);
  if (ok)
    obj = cs->Obj;

  if (ok){
    if(mode == cRepDotAreaType) { /* assume all atoms "visible" for area comp. */
      visFlag = true;
    } else {
      visFlag = cs->hasRep(cRepDotBit);
    }
  }
  if(!ok || !visFlag) {
    OOFreeP(I);
    return (NULL);              /* skip if no dots are visible */
  }

  RepInit(G, &I->R);

  I->dotSize = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_dot_radius);

  I->A = NULL;
  I->T = NULL;
  I->F = NULL;
  I->V = NULL;
  I->VC = NULL;
  I->VN = NULL;
  I->Atom = NULL;
  I->R.fRecolor = NULL;
  I->shaderCGO = 0;

  I->Width = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_dot_width);
  cullByFlag = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_trim_dots);      /* are we using flags 24 & 25 */

  dot_color = SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_dot_color);   /* are we using flags 24 & 25 */
  inclH = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_dot_hydrogens);       /* are we ignoring hydrogens? */
  if(SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_dot_solvent)) {    /* are we generating a solvent surface? */
    solv_rad = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_solvent_radius); /* if so, get solvent radius */
  }

  /* get current dot sampling */
  ds = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_dot_density);

  max_vdw += solv_rad;


/* Note: significantly affects the accuracy of our area comp. */
  if(ds < 0)
    ds = 0;
  if(ds > 4)
    ds = 4;
  sp = G->Sphere->Sphere[ds];

  I->R.fRender = (void (*)(struct Rep *, RenderInfo * info)) RepDotRender;
  I->R.fFree = (void (*)(struct Rep *)) RepDotFree;
  I->R.obj = (CObject *) obj;
  I->R.cs = cs;

  I->V = pymol::malloc<float>(cs->NIndex * sp->nDot * 10);
  CHECKOK(ok, I->V);

  if(ok && mode == cRepDotAreaType) { /* in area mode, we need to export save addl. info 
                                 * such as the normal vectors, the partial area, 
                                 * the originating atom, etc. */
    I->A = pymol::malloc<float>(cs->NIndex * sp->nDot);
    CHECKOK(ok, I->A);
    if (ok)
      I->T = pymol::malloc<int>(cs->NIndex * sp->nDot);
    CHECKOK(ok, I->T);
    if (ok)
      I->F = pymol::malloc<int>(cs->NIndex * sp->nDot);
    CHECKOK(ok, I->F);
    if (ok)
      I->VN = pymol::malloc<float>(cs->NIndex * sp->nDot * 3);
    CHECKOK(ok, I->VN);
    if (ok)
      I->Atom = pymol::malloc<int>(cs->NIndex * sp->nDot);
    CHECKOK(ok, I->Atom);
    if (ok){
      aa = I->A;
      tp = I->T;
      tf = I->F;
      ati = I->Atom;
      inclH = true;
      cullByFlag = true;
    }
  }
  vn = I->VN;
  I->N = 0;
  lastColor = -1;
  colorCnt = 0;
  if (ok)
    map = MapNew(G, max_vdw, cs->Coord, cs->NIndex, NULL);
  CHECKOK(ok, map);
  v = I->V;
  if(ok && map) {
    ok &= MapSetupExpress(map);
    for(a = 0; ok && a < cs->NIndex; a++) {
      atm = cs->IdxToAtm[a];
      ai1 = obj->AtomInfo + atm;
      if((ai1->visRep & cRepDotBit) || mode == cRepDotAreaType)
        if((inclH || (!ai1->isHydrogen())) &&
           ((!cullByFlag) || (!(ai1->flags & cAtomFlag_exfoliate)))) {
          c1 = AtomSettingGetWD(G, ai1, cSetting_dot_color, dot_color);

          /* If we are culling, flag 24 controls which atoms 
             will have dot surfaces generated for them.
           */
          if(c1 == -1) {
            c1 = ai1->color;
          }
          v0 = cs->Coord + 3 * a;
          vdw = ai1->vdw + solv_rad;
          for(b = 0; b < sp->nDot; b++) {
            v1[0] = v0[0] + vdw * sp->dot[b][0];
            v1[1] = v0[1] + vdw * sp->dot[b][1];
            v1[2] = v0[2] + vdw * sp->dot[b][2];

            MapLocus(map, v1, &h, &k, &l);

            flag = true;

            i = *(MapEStart(map, h, k, l));
            if(i) {
              j = map->EList[i++];
              while(j >= 0) {
                ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                if((inclH || (!(ai2->isHydrogen()))) &&
                   ((!cullByFlag) || (!(ai2->flags & cAtomFlag_ignore))))
                  /* If we are cullilng, flag 25 controls which atoms 
                     are considered "present" in the surface area 
                     calculation (i.e. able to occlude surface) */
                  if(j != a)
                    if(within3f(cs->Coord + 3 * j, v1, ai2->vdw + solv_rad)) {
                      flag = false;
                      break;
                    }
                j = map->EList[i++];
              }
            }
            if(flag) {
              switch (mode) {
              case cRepDotNormal:

                if((lastColor != c1) || ColorCheckRamped(G, c1)) {      /* new color */
                  if(countPtr)  /* after first pass */
                    *countPtr = (float) colorCnt;       /* save count */
                  colorCnt = 1;
                  countPtr = v++;
                  vc = ColorGet(G, c1); /* save new color */
                  lastColor = c1;
                  if(ColorCheckRamped(G, c1)) {
                    ColorGetRamped(G, c1, v1, v, state);
                    v += 3;
                  } else {
                    *(v++) = *(vc++);
                    *(v++) = *(vc++);
                    *(v++) = *(vc++);
                  }
                } else
                  colorCnt++;
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
                *(aa++) = vdw * vdw * sp->area[b];      /* area */
                *(tp++) = ai1->customType;      /* numeric type */
                *(tf++) = ai1->flags;   /* flags */
                *(vn++) = sp->dot[b][0];
                *(vn++) = sp->dot[b][1];
                *(vn++) = sp->dot[b][2];
                *(ati++) = atm;
                I->N++;
                break;
              }
            }
          }
        }
      ok &= !G->Interrupt;
    }
    if(countPtr)
      *countPtr = (float) colorCnt;     /* save count */
    MapFree(map);
  }
  if (ok)
    I->V = ReallocForSure(I->V, float, (v - I->V));
  CHECKOK(ok, I->V);
  if(ok && mode == cRepDotAreaType) {
    I->A = ReallocForSure(I->A, float, (aa - I->A));
    CHECKOK(ok, I->A);
    if (ok)
      I->T = ReallocForSure(I->T, int, (tp - I->T));
    CHECKOK(ok, I->T);
    if (ok)
      I->F = ReallocForSure(I->F, int, (tf - I->F));
    CHECKOK(ok, I->F);
    if (ok)
      I->VN = ReallocForSure(I->VN, float, (vn - I->VN));
    CHECKOK(ok, I->VN);
    if (ok)
      I->Atom = ReallocForSure(I->Atom, int, (ati - I->Atom));
    CHECKOK(ok, I->Atom);
  }
  if (!ok){
    RepDotFree(I);
    I = NULL;
  }
  return (Rep *) I;
}
