
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
#include"RepRibbon.h"
#include"Color.h"
#include"Setting.h"
#include"Word.h"
#include"Scene.h"
#include"main.h"
#include"Feedback.h"
#include"ShaderMgr.h"
#include"CGO.h"
#include "Lex.h"

struct RepRibbon : Rep {
  using Rep::Rep;

  ~RepRibbon() override;

  cRep_t type() const override { return cRepRibbon; }
  void render(RenderInfo* info) override;

  float ribbon_width;
  float radius;
  CGO *shaderCGO;
  CGO *primitiveCGO;
  bool shaderCGO_has_cylinders;
};

#include"ObjectMolecule.h"

RepRibbon::~RepRibbon()
{
  CGOFree(primitiveCGO);
  CGOFree(shaderCGO);
}

void RepRibbon::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  int ok = true;
  short use_shader = SettingGetGlobal_b(G, cSetting_ribbon_use_shader) &&
                     SettingGetGlobal_b(G, cSetting_use_shaders);
  bool ribbon_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) &&
                             SettingGet<bool>(G, I->cs->Setting.get(),
                                                 I->obj->Setting.get(),
                                                 cSetting_ribbon_as_cylinders);

  if(ray) {
#ifndef _PYMOL_NO_RAY
    CGORenderRay(I->primitiveCGO, ray, info, NULL, NULL, I->cs->Setting.get(), I->obj->Setting.get());
#endif
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      CGORenderGLPicking(I->shaderCGO ? I->shaderCGO : I->primitiveCGO, info, &I->context, I->cs->Setting.get(), I->obj->Setting.get(), I);
    } else {
      if (!use_shader && I->shaderCGO){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }
      if (I->shaderCGO && (ribbon_as_cylinders ^ I->shaderCGO_has_cylinders)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (use_shader){
	if (!I->shaderCGO){
          CGO *convertcgo = NULL;
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok)
	    I->shaderCGO->use_shader = true;

          if (ok)
            ok &= CGOResetNormal(I->shaderCGO, true);
          if (ribbon_as_cylinders){
            if (ok) ok &= CGOEnable(I->shaderCGO, GL_CYLINDER_SHADER);
            if (ok) ok &= CGOSpecial(I->shaderCGO, CYLINDER_WIDTH_FOR_RIBBONS);
            convertcgo = CGOConvertLinesToCylinderShader(I->primitiveCGO, I->shaderCGO);
            if (ok) ok &= CGOAppendNoStop(I->shaderCGO, convertcgo);
            if (ok) ok &= CGODisable(I->shaderCGO, GL_CYLINDER_SHADER);
            if (ok) ok &= CGOStop(I->shaderCGO);
          } else {
            int trilines = SettingGetGlobal_b(G, cSetting_trilines);
            int shader = trilines ? GL_TRILINES_SHADER : GL_LINE_SHADER;
            if (ok) ok &= CGOEnable(I->shaderCGO, shader);
            if (ok) ok &= CGODisable(I->shaderCGO, CGO_GL_LIGHTING);
            if (trilines) {
              if (ok) ok &= CGOSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON);
              convertcgo = CGOConvertToTrilinesShader(I->primitiveCGO, I->shaderCGO);
            } else {
              convertcgo = CGOConvertToLinesShader(I->primitiveCGO, I->shaderCGO);
            }
            if (ok) ok &= CGOAppendNoStop(I->shaderCGO, convertcgo);
            if (ok) ok &= CGODisable(I->shaderCGO, shader);
            if (ok) ok &= CGOStop(I->shaderCGO);
          }
          I->shaderCGO_has_cylinders = ribbon_as_cylinders;
          CGOFreeWithoutVBOs(convertcgo);
          I->shaderCGO->use_shader = true;
        }
        CGORenderGL(I->shaderCGO, NULL, I->cs->Setting.get(), I->obj->Setting.get(), info, I);
        return;
      } else {
        CGORenderGL(I->primitiveCGO, NULL, I->cs->Setting.get(), I->obj->Setting.get(), info, I);
        return;
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->invalidate(cRepInvPurge);
    I->cs->Active[cRepRibbon] = false;
  }
}

Rep *RepRibbonNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj;
  int a, b, a1, a2, *i, *s, *at, *seg, nAt, *atp;
  float *v, *v1, *v2, *v3;
  float *pv = NULL;
  float *dv = NULL;
  float *nv = NULL;
  float *tv = NULL;

  float f0, f1, f2, f3, f4;
  float *d;
  float *dl = NULL;
  int nSeg;
  int sampling;
  float power_a = 5;
  float power_b = 5;
  float throw_;
  float dev;
  int trace, trace_mode;
  int ribbon_color;
  int na_mode;
  AtomInfoType *ai, *last_ai = NULL;
  AtomInfoType *trailing_O3p_ai = NULL, *leading_O5p_ai = NULL;
  int trailing_O3p_a = 0, leading_O5p_a = 0, leading_O5p_a1 = 0;

  // skip if not visible
  if(!cs->hasRep(cRepRibbonBit))
    return NULL;

  auto I = new RepRibbon(cs, state);

  obj = cs->Obj;

  power_a = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_power);
  power_b = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_power_b);
  throw_ = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_throw);
  int trace_ostate = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_trace_atoms);
  trace_mode = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_trace_atoms_mode);
  na_mode =
    SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_nucleic_acid_mode);

  ribbon_color =
    SettingGet_color(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_color);

  sampling = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_sampling);
  if(sampling < 1)
    sampling = 1;
  I->radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_radius);
  I->ribbon_width = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_width);

  /* find all of the CA points */

  auto const nAtIndex = cs->getNIndex(); // was NAtIndex
  at = pymol::malloc<int>(nAtIndex * 2);
  pv = pymol::malloc<float>(nAtIndex * 6);
  seg = pymol::malloc<int>(nAtIndex * 2);

  i = at;
  v = pv;
  s = seg;

  nAt = 0;
  nSeg = 0;
  a2 = -1;
  for(a1 = 0; a1 < obj->NAtom; ++a1) {
    a = cs->atmToIdx(a1);
    if(a >= 0) {
      ai = obj->AtomInfo + a1;
      if(GET_BIT(obj->AtomInfo[a1].visRep,cRepRibbon)) {
        trace = AtomSettingGetWD(G, ai, cSetting_ribbon_trace_atoms, trace_ostate);

        if(trace || ((obj->AtomInfo[a1].protons == cAN_C) &&
                     (WordMatchExact(G, G->lex_const.CA, obj->AtomInfo[a1].name, true)) &&
                     !AtomInfoSameResidueP(G, last_ai, ai))) {
          PRINTFD(G, FB_RepRibbon)
            " RepRibbon: found atom in %d; a1 %d a2 %d\n", obj->AtomInfo[a1].resv, a1, a2
            ENDFD;

          // auto-detect CA-only models
          if (!ai->bonded)
            trace = true;

          if(a2 >= 0) {
            if(trace) {
              if(!AtomInfoSequential
                 (G, obj->AtomInfo + a2, obj->AtomInfo + a1, trace_mode))
                a2 = -1;
            } else {
              if(!ObjectMoleculeCheckBondSep(obj, a1, a2, 3))   /* CA->N->C->CA = 3 bonds */
                a2 = -1;
            }
          }
          PRINTFD(G, FB_RepRibbon)
            " RepRibbon: found atom in %d; a1 %d a2 %d\n", obj->AtomInfo[a1].resv, a1, a2
            ENDFD;
          last_ai = ai;
          if(a2 < 0)
            nSeg++;
          *(s++) = nSeg;
          nAt++;
          *(i++) = a;
          v1 = cs->coordPtr(a);
          *(v++) = *(v1++);
          *(v++) = *(v1++);
          *(v++) = *(v1++);

          a2 = a1;
        } else if((((na_mode != 1) && (ai->protons == cAN_P) &&
                    (WordMatchExact(G, G->lex_const.P, ai->name, true))) ||
                   ((na_mode == 1) && (ai->protons == cAN_C) &&
                    (WordMatchExact(G, "C4*", LexStr(G, ai->name), 1) ||
                     WordMatchExact(G, "C4'", LexStr(G, ai->name), 1)))) &&
                  !AtomInfoSameResidueP(G, last_ai, ai)) {
          if(a2 >= 0) {
            if(!ObjectMoleculeCheckBondSep(obj, a1, a2, 6)) {   /* six bonds between phosphates */
              if(trailing_O3p_ai && ((na_mode == 2) || (na_mode == 4))) {
                /* 3' end of nucleic acid */
                *(s++) = nSeg;
                nAt++;
                *(i++) = trailing_O3p_a;
                v1 = cs->coordPtr(trailing_O3p_a);
                *(v++) = *(v1++);
                *(v++) = *(v1++);
                *(v++) = *(v1++);
              }
              a2 = -1;
            }
          }

          trailing_O3p_ai = NULL;

          if(leading_O5p_ai && (a2 < 0) && ((na_mode == 3) || (na_mode == 4))) {
            if((!AtomInfoSameResidueP(G, ai, leading_O5p_ai)) &&
               ObjectMoleculeCheckBondSep(obj, a1, leading_O5p_a1, 5)) {
              nSeg++;
              *(s++) = nSeg;
              nAt++;
              *(i++) = leading_O5p_a;
              v1 = cs->coordPtr(leading_O5p_a);
              *(v++) = *(v1++);
              *(v++) = *(v1++);
              *(v++) = *(v1++);
              a2 = leading_O5p_a1;
            }
          }
          leading_O5p_ai = NULL;
          last_ai = ai;
          if(a2 < 0)
            nSeg++;
          *(s++) = nSeg;
          nAt++;
          *(i++) = a;
          v1 = cs->coordPtr(a);
          *(v++) = *(v1++);
          *(v++) = *(v1++);
          *(v++) = *(v1++);

          a2 = a1;
        } else if((a2 >= 0) &&
                  last_ai &&
                  (ai->protons == cAN_O) &&
                  (last_ai->protons == cAN_P) &&
                  ((na_mode == 2) || (na_mode == 4)) &&
                  (WordMatchExact(G, "O3'", LexStr(G, ai->name), 1) ||
                   WordMatchExact(G, "O3*", LexStr(G, ai->name), 1)) &&
                  AtomInfoSameResidueP(G, last_ai, ai) &&
                  ObjectMoleculeCheckBondSep(obj, a1, a2, 5)) {
          trailing_O3p_ai = ai;
          trailing_O3p_a = a;
        } else if((ai->protons == cAN_O) &&
                  ((na_mode == 3) || (na_mode == 4)) &&
                  (WordMatchExact(G, "O5'", LexStr(G, ai->name), 1) ||
                   WordMatchExact(G, "O5*", LexStr(G, ai->name), 1))) {
          leading_O5p_ai = ai;
          leading_O5p_a = a;
          leading_O5p_a1 = a1;
        }
      }
    }

  }
  if(trailing_O3p_ai && ((na_mode == 2) || (na_mode == 4))) {
    /* 3' end of nucleic acid */
    *(s++) = nSeg;
    nAt++;
    *(i++) = trailing_O3p_a;
    v1 = cs->coordPtr(trailing_O3p_a);
    *(v++) = *(v1++);
    *(v++) = *(v1++);
    *(v++) = *(v1++);
  }
  PRINTFD(G, FB_RepRibbon)
    " RepRibbon: nAt %d\n", nAt ENDFD;

  if(nAt) {
    /* compute differences and normals */

    s = seg;
    v = pv;

    dv = pymol::malloc<float>(nAt * 6);
    nv = pymol::malloc<float>(nAt * 6);
    dl = pymol::malloc<float>(nAt * 2);
    v1 = dv;
    v2 = nv;
    d = dl;

    for(a = 0; a < (nAt - 1); a++) {
      if(*s == *(s + 1)) {
        float d_1;
        subtract3f(v + 3, v, v1);
        *d = (float) length3f(v1);
        if(*d > R_SMALL4) {
          d_1 = 1.0F / (*d);
          scale3f(v1, d_1, v2);
        } else if(a) {
          copy3f(v2 - 3, v2);
        } else {
          zero3f(v2);
        }
      }
      d++;
      v += 3;
      v1 += 3;
      v2 += 3;
      s++;
    }

    /* compute tangents */

    s = seg;
    v = nv;

    tv = pymol::malloc<float>(nAt * 6 + 6);
    v1 = tv;

    *(v1++) = *(v++);           /* first segment */
    *(v1++) = *(v++);
    *(v1++) = *(v++);
    s++;

    for(a = 1; a < (nAt - 1); a++) {
      if((*s == *(s - 1)) && (*s == *(s + 1))) {
        add3f(v, (v - 3), v1);
        normalize3f(v1);
      } else if(*s == *(s - 1)) {
        *(v1) = *(v - 3);       /* end a segment */
        *(v1 + 1) = *(v - 2);
        *(v1 + 2) = *(v - 1);
      } else if(*s == *(s + 1)) {
        *(v1) = *(v);           /* new segment */
        *(v1 + 1) = *(v + 1);
        *(v1 + 2) = *(v + 2);
      }
      v += 3;
      v1 += 3;
      s++;
    }

    *(v1++) = *(v - 3);         /* last segment */
    *(v1++) = *(v - 2);
    *(v1++) = *(v - 1);

  }

  /* okay, we now have enough info to generate smooth interpolations */

  I->shaderCGO = 0;
  I->primitiveCGO = 0;

  if(nAt) {
    v1 = pv;                    /* points */
    v2 = tv;                    /* tangents */
    v3 = dv;                    /* direction vector */
    d = dl;
    s = seg;
    atp = at;

    I->primitiveCGO = CGONew(G);
    CGOSpecialWithArg(I->primitiveCGO, LINE_LIGHTING, 0.f);

    float alpha = 1.f - SettingGet_f(G, NULL, I->obj->Setting.get(), cSetting_ribbon_transparency);
    if(fabs(alpha-1.0) < R_SMALL4)
      alpha = 1.0F;
    CGOAlpha(I->primitiveCGO, alpha);  // would be good to set these at render time instead
    CGOSpecial(I->primitiveCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON);

    if (alpha < 1) {
      I->setHasTransparency();
    }

    // This is required for immediate mode rendering
    CGOBegin(I->primitiveCGO, GL_LINES);

    auto const get_color = [&](AtomInfoType const* ai) {
      auto c = AtomSettingGetWD(G, ai, cSetting_ribbon_color, ribbon_color);
      return (c != cColorDefault) ? c : ai->color;
    };

    float origV1[14], origV2[14], *origV = origV2;
    bool origVis1 = false;
    for(a = 0; a < (nAt - 1); a++) {
      if(*s == *(s + 1)) {
        int const atm[2] = {cs->IdxToAtm[atp[0]], cs->IdxToAtm[atp[1]]};

        AtomInfoType const *ai1 = obj->AtomInfo + atm[0];
        AtomInfoType const *ai2 = obj->AtomInfo + atm[1];

        int const color[2] = {get_color(ai1), get_color(ai2)};
        int const atmpk[2] = {
            ai1->masked ? cPickableNoPick : cPickableAtom,
            ai2->masked ? cPickableNoPick : cPickableAtom,
        };

        dev = throw_ * (*d);
        for(b = 0; b < sampling; b++) { /* needs optimization */
          origVis1 = !origVis1;
          if (origVis1){
            origV = origV1;
          } else {
            origV = origV2;
          }

          size_t const i1 = (b + 0) * 2 < sampling ? 0 : 1;
          size_t const i2 = (b + 1) * 2 > sampling ? 1 : 0;
          assert(i1 <= i2);

          /* 
             14 floats per v:
             v[0] = index1
             v[1-3] = color1
             v[4-6] = vertex1
             v[7] = index2
             v[8-10] = color2
             v[11-13] = vertex2
          */


          f0 = ((float) b) / sampling;  /* fraction of completion */
          f0 = smooth(f0, power_a);     /* bias sampling towards the center of the curve */

          // start of line/cylinder 

          f1 = 1.0F - f0;
          f2 = smooth(f0, power_b);
          f3 = smooth(f1, power_b);
          f4 = dev * f2 * f3;   /* displacement magnitude */

          /* store vertex1 v[4-6] */
          origV[4] = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);
          origV[5] = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);
          origV[6] = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);

          bool isRamped =
              ColorGetCheckRamped(G, color[i1], origV + 4, origV + 1, state);

          f0 = ((float) b + 1) / sampling;
          f0 = smooth(f0, power_a);

          /* end of line/cylinder */

          f1 = 1.0F - f0;
          f2 = smooth(f0, power_b);
          f3 = smooth(f1, power_b);
          f4 = dev * f2 * f3;   /* displacement magnitude */

          /* store vertex2 v[11-13] */
          origV[11] = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);
          origV[12] = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);
          origV[13] = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);

          if (ColorGetCheckRamped(G, color[i2], origV + 11, origV + 8, state)) {
            isRamped = true;
          }

          CGOPickColor(I->primitiveCGO, atm[i1], atmpk[i1]);
          CGOColorv(I->primitiveCGO, origV + 1);
          I->primitiveCGO->add<cgo::draw::splitline>(origV + 4, origV + 11,
              origV + 8, atm[i2], atmpk[i2], isRamped, false, false);
        }
      }
      v1 += 3;
      v2 += 3;
      v3 += 3;
      d++;
      atp += 1;
      s++;
    }
    CGOEnd(I->primitiveCGO);
    CGOSpecialWithArg(I->primitiveCGO, LINE_LIGHTING, 1.f);
    CGOStop(I->primitiveCGO);

    FreeP(dv);
    FreeP(dl);
    FreeP(tv);
    FreeP(nv);
  } else {
    delete I;
    I = NULL;
  }

  FreeP(at);
  FreeP(seg);
  FreeP(pv);

  return (Rep *) I;
}

void RepRibbonRenderImmediate(CoordSet * cs, RenderInfo * info)
{
#ifndef PURE_OPENGL_ES_2
  /* performance optimized to provide a simple C-alpha trace -- no smoothing */

  PyMOLGlobals *G = cs->G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    ObjectMolecule *obj = cs->Obj;
    int active = false;
    int nAtIndex = obj->NAtom;
    const AtomInfoType *obj_AtomInfo = obj->AtomInfo.data();
    const AtomInfoType *ai, *last_ai = NULL;
    int trace, trace_ostate =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_trace_atoms);
    int trace_mode =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_trace_atoms_mode);
    int na_mode =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_nucleic_acid_mode);
    float ribbon_width =
      SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_ribbon_width);
    int a1, a2 = -1;
    int color, last_color = -9;

    glLineWidth(ribbon_width);
    SceneResetNormal(G, true);
    if(!info->line_lighting)
      glDisable(GL_LIGHTING);
    glBegin(GL_LINE_STRIP);
    for(a1 = 0; a1 < nAtIndex; a1++) {
      auto a = cs->atmToIdx(a1);
      if(a >= 0) {
        ai = obj_AtomInfo + a1;
        if(GET_BIT(ai->visRep,cRepRibbon)) {
          trace = AtomSettingGetWD(G, ai, cSetting_ribbon_trace_atoms, trace_ostate);

          if(trace || ((ai->protons == cAN_C) &&
                       (WordMatchExact(G, G->lex_const.CA, ai->name, true)) &&
                       !AtomInfoSameResidueP(G, last_ai, ai))) {
            if(a2 >= 0) {
              if(trace) {
                if(!AtomInfoSequential
                   (G, obj_AtomInfo + a2, obj_AtomInfo + a1, trace_mode))
                  a2 = -1;
              } else {
                if(!ObjectMoleculeCheckBondSep(obj, a1, a2, 3)) /* CA->N->C->CA = 3 bonds */
                  a2 = -1;
              }
            }
            if(a2 == -1) {
              glEnd();
              glBegin(GL_LINE_STRIP);
            }
            color = ai->color;
            if(color != last_color) {
              last_color = color;
              glColor3fv(ColorGet(G, color));
            }
            glVertex3fv(cs->coordPtr(a));
            active = true;
            last_ai = ai;
            a2 = a1;
          } else if((((na_mode != 1) && (ai->protons == cAN_P) &&
                      (WordMatchExact(G, G->lex_const.P, ai->name, true))) ||
                     ((na_mode == 1) && (ai->protons == cAN_C) &&
                      (WordMatchExact(G, "C4*", LexStr(G, ai->name), 1) ||
                       WordMatchExact(G, "C4'", LexStr(G, ai->name), 1)))) &&
                    !AtomInfoSameResidueP(G, last_ai, ai)) {
            if(a2 >= 0) {
              if(!ObjectMoleculeCheckBondSep(obj, a1, a2, 6)) { /* six bonds between phosphates */
                a2 = -1;
              }
            }

            if(a2 == -1) {
              glEnd();
              glBegin(GL_LINE_STRIP);
            }
            color = ai->color;
            if(color != last_color) {
              last_color = color;
              glColor3fv(ColorGet(G, color));
            }
            glVertex3fv(cs->coordPtr(a));
            active = true;
            last_ai = ai;
            a2 = a1;
          }
        }
      }
    }
    glEnd();
    glEnable(GL_LIGHTING);
    if(!active)
      cs->Active[cRepRibbon] = false;
  }
#endif
}
