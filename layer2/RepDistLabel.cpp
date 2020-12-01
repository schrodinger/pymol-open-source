
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
#include"RepDistLabel.h"
#include"RepLabel.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"
#include"Text.h"

#include"Word.h"
#include"CGO.h"
#include"CoordSet.h"
#include"Util.h"

typedef char DistLabel[12];

struct RepDistLabel : Rep {
  using Rep::Rep;

  ~RepDistLabel() override;

  cRep_t type() const override { return cRepLabel; }
  void render(RenderInfo* info) override;

  float* V = nullptr;
  int N = 0;
  DistLabel *L;
  DistSet *ds;
  int OutlineColor;
  CGO* shaderCGO = nullptr;
  int texture_font_size = 0;
};

#define SHADERCGO I->shaderCGO

#include"ObjectDist.h"

RepDistLabel::~RepDistLabel()
{
  auto I = this;
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
  }
  VLAFreeP(I->V);
  VLAFreeP(I->L);
}

void RepDistLabel::render(RenderInfo* info)
{
  auto I = this;
  auto* obj = getObj();
  CRay *ray = info->ray;
  auto pick = info->pick;
  float *v = I->V;
  int c = I->N;
  DistLabel *l = I->L;
  int n = 0;
  int font_id = SettingGet_i(G, NULL, obj->Setting.get(), cSetting_label_font_id);
  float font_size = SettingGet_f(G, NULL, obj->Setting.get(), cSetting_label_size);
  int float_text = SettingGet_i(G, NULL, obj->Setting.get(), cSetting_float_labels);
  int ok = true;
  short use_shader = SettingGetGlobal_b(G, cSetting_use_shaders);
  if (I->MaxInvalid >= cRepInvRep)
    return;
  font_id = SettingCheckFontID(G, NULL, obj->Setting.get(), font_id);

  if (I->shaderCGO && font_size < 0.f){
    int size;
    if (InvalidateShaderCGOIfTextureNeedsUpdate(G, font_size, I->texture_font_size, &size)){
      CGOFree(I->shaderCGO);
      I->texture_font_size = size;
    }
  }

  auto color = SettingGet_color(G, nullptr, obj->Setting.get(), cSetting_label_color);
  if (color < 0               //
      && color != cColorFront //
      && color != cColorBack) {
    color = obj->Color;
  }

  if(ray) {
    TextSetOutlineColor(G, I->OutlineColor);
      TextSetColor(G, ColorGet(G, color));

    while(c--) {
      TextSetPos(G, v);
      TextRenderRay(G, ray, font_id, l[n], font_size, v + 3, false, 0);
      v += 6;
      n++;
    }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      if (I->shaderCGO){
        if(float_text)
          glDisable(GL_DEPTH_TEST);
	CGORenderGLPicking(I->shaderCGO, info, &I->context, NULL, NULL);
        if(float_text)
          glEnable(GL_DEPTH_TEST);
	return;
    } else {
	Pickable *p = I->P;
        TextSetIsPicking(G, true);
	SceneSetupGLPicking(G);
	if(c) {
	  if(float_text)
	    glDisable(GL_DEPTH_TEST);
	  
	  while(c--) {
	    if(*l) {
	      TextSetPos(G, v);
              p++;
              AssignNewPickColor(nullptr, pick, TextGetColorUChar4uv(G),
                  &I->context, p->index, p->bond);
              TextSetColorFromUColor(G);
	      TextSetLabelBkgrdInfo(G, 1.f, 1.2f, NULL);
	      TextSetLabelPosIsSet(G, 0);
	      if (!TextRenderOpenGL(G, info, font_id, l[n], font_size, v + 3, false, 0, 1, 0)){
                TextSetIsPicking(G, false);
		return;
              }
	      n++;
	    }
	    v += 6;
	  }
	  if(float_text)
	    glEnable(GL_DEPTH_TEST);
	}
        TextSetIsPicking(G, false);
      }
    } else {
	Pickable *p = I->P;
	
        if (use_shader){
	if (!I->shaderCGO){
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok){
	    I->shaderCGO->use_shader = true;
	  }
	} else {
	    info->texture_font_size = I->texture_font_size;
	    CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, I);
	  return;
	}
	} else if (I->shaderCGO){
	  CGOFree(I->shaderCGO);
	  if(float_text)
	    glDisable(GL_DEPTH_TEST);
	}

      TextSetOutlineColor(G, I->OutlineColor);
        TextSetColor(G, ColorGet(G, color));
      while(c--) {
	p++;
	if (ok && I->shaderCGO)
	  ok &= CGOPickColor(I->shaderCGO, p->index, p->bond);
        TextSetPos(G, v);
	TextSetLabelBkgrdInfo(G, 1.f, 1.2f, NULL);
	TextSetLabelPosIsSet(G, 0);
	if (!TextRenderOpenGL(G, info, font_id, l[n], font_size, v + 3, false, 0, 1, SHADERCGO))
	  return;
        v += 6;
        n++;
      }
      if (ok && I->shaderCGO){
	ok &= CGOStop(I->shaderCGO);
	if (ok){
	  CGO *tmpCGO = CGONew(G);
          CGOEnable(tmpCGO, GL_LABEL_SHADER);
	  CGODisable(tmpCGO, GL_DEPTH_TEST_IF_FLOATING);
          CGOSpecial(tmpCGO, SET_LABEL_SCALE_UNIFORMS);
	  CGO *convertcgo = CGOConvertToLabelShader(I->shaderCGO, tmpCGO);
	  if (!convertcgo) {
	    CGOFree(tmpCGO);
	    CGOFree(I->shaderCGO);
	    return;
	  }
	  CGOAppendNoStop(tmpCGO, convertcgo);
	  CGOFreeWithoutVBOs(convertcgo);
	  CGOEnable(tmpCGO, GL_DEPTH_TEST_IF_FLOATING);
	  CGODisable(tmpCGO, GL_LABEL_SHADER);
          CGOStop(tmpCGO);
	  convertcgo = tmpCGO;
	  CHECKOK(ok, convertcgo);
	  CGOFree(I->shaderCGO);
	  I->shaderCGO = convertcgo;
	  convertcgo = NULL;
	}
	if (ok && I->shaderCGO){
	  I->shaderCGO->use_shader = true;
	  I->render(info); // recursion !?
	  return;
	}
      }
      if(float_text)
        glEnable(GL_DEPTH_TEST);
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->ds->Rep[cRepLabel] = NULL;
    delete I;
  }
}

Rep *RepDistLabelNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->G;
  int a;
  int n = 0;
  float *v, *v1, *v2, *v3, d[3], di;
  char buffer[255];
  const float *lab_pos =
    SettingGet_3fv(G, NULL, ds->Obj->Setting.get(), cSetting_label_position);
  int default_digits =
    SettingGet_i(G, NULL, ds->Obj->Setting.get(), cSetting_label_digits);
  Pickable *rp = NULL;

  if(!(ds->NIndex || ds->NAngleIndex || ds->NDihedralIndex)) {
    ds->LabCoord.clear();
    ds->LabPos.clear();
    return nullptr;
  }

  default_digits = pymol::clamp(default_digits, 0, 10);

  auto I = new RepDistLabel(ds->Obj, state);

  I->ds = ds;

  I->OutlineColor =
    SettingGet_i(G, NULL, ds->Obj->Setting.get(), cSetting_label_outline_color);

  int nLabel = 0;
  int ok = true;
  if(ds->NIndex || ds->NAngleIndex || ds->NDihedralIndex) {
    nLabel = (ds->NIndex / 2 + ds->NAngleIndex / 5 + ds->NDihedralIndex / 6);

    ds->LabCoord.resize(nLabel);
    ds->LabPos.resize(nLabel);

    if(SettingGet_b(G, NULL, ds->Obj->Setting.get(), cSetting_pickable)) {
      I->P = pymol::malloc<Pickable>(nLabel + 1);
      CHECKOK(ok, I->P);
      if (ok)
	rp = I->P + 1;          /* skip first record! */
    }

    if (ok)
      I->V = VLAlloc(float, 3 * (ds->NIndex / 2 + ds->NAngleIndex / 5) + 1);
    CHECKOK(ok, I->V);
    if (ok)
      I->L = VLAlloc(DistLabel, (ds->NIndex / 2 + ds->NAngleIndex / 5) + 1);
    CHECKOK(ok, I->L);

    n = 0;
    auto* lc = ds->LabCoord.data();

    if(ds->NIndex) {
      int digits = SettingGet_i(G, NULL, ds->Obj->Setting.get(),
                                cSetting_label_distance_digits);
      WordType format;
      if(digits < 0)
        digits = default_digits;
      if(digits > 10)
        digits = 10;
      sprintf(format, "%c0.%df", '%', digits);
      for(a = 0; ok && a < ds->NIndex; a = a + 2) {
        v1 = ds->Coord + 3 * a;
        v2 = ds->Coord + 3 * (a + 1);
        average3f(v2, v1, d);
        di = (float) diff3f(v1, v2);
        sprintf(buffer, format, di);

        VLACheck(I->V, float, 6 * n + 5);
	CHECKOK(ok, I->V);
	if (ok)
	  VLACheck(I->L, DistLabel, n);
	CHECKOK(ok, I->L);

	if (!ok)
	  break;
        v = I->V + 6 * n;
        UtilNCopy(I->L[n], buffer, sizeof(DistLabel));
        copy3f(d, v);
        std::copy_n(v, 3, lc->data());
        lc++;
        if(!ds->LabPos.empty()) {
          const auto& lp = ds->LabPos[n];
          switch (lp.mode) {
          case 1:
            add3f(lp.offset, v, v);
            copy3f(lab_pos, v + 3);
            break;
          default:
            copy3f(lab_pos, v + 3);
            break;
          }
        } else {
          copy3f(lab_pos, v + 3);
        }

        if(rp) {
          rp->index = n;        /* label index */
          rp->bond = cPickableLabel;    /* label indicator */
          rp++;
        }

        n++;
      }
    }

    if(ok && ds->NAngleIndex) {

      float d1[3], d2[3], n1[3], n2[3];
      float avg[3];

      float l1, l2;
      float radius;
      int digits = SettingGet_i(G, NULL, ds->Obj->Setting.get(),
                                cSetting_label_angle_digits);
      WordType format;
      if(digits < 0)
        digits = default_digits;
      if(digits > 10)
        digits = 10;
      sprintf(format, "%c0.%df", '%', digits);

      for(a = 0; ok && a < ds->NAngleIndex; a = a + 5) {
        v1 = ds->AngleCoord + 3 * a;
        v2 = ds->AngleCoord + 3 * (a + 1);
        v3 = ds->AngleCoord + 3 * (a + 2);
        subtract3f(v1, v2, d1);
        subtract3f(v3, v2, d2);

        normalize23f(d1, n1);
        normalize23f(d2, n2);

        average3f(n1, n2, avg);

        l1 = (float) length3f(d1);
        l2 = (float) length3f(d2);

        if(l1 > l2)
          radius = l2;
        else
          radius = l1;
        radius *=
          SettingGet_f(G, NULL, ds->Obj->Setting.get(),
                       cSetting_angle_size) * SettingGet_f(G, NULL,
                                                           ds->Obj->Setting.get(),
                                                           cSetting_angle_label_position);

        normalize3f(avg);
        if((avg[0] == 0.0F) && (avg[1] == 0.0F) && (avg[2] == 0.0F))
          avg[0] = 1.0F;

        scale3f(avg, radius, avg);
        add3f(v2, avg, avg);

        di = (float) (180.0F * get_angle3f(d1, d2) / PI);
        sprintf(buffer, format, di);

        VLACheck(I->V, float, 6 * n + 5);
	CHECKOK(ok, I->V);
	if (ok)
	  VLACheck(I->L, DistLabel, n);
	CHECKOK(ok, I->L);
	if (!ok)
	  break;
        v = I->V + 6 * n;
        UtilNCopy(I->L[n], buffer, sizeof(DistLabel));
        copy3f(avg, v);
        std::copy_n(v, 3, lc->data());
        lc++;
        if(!ds->LabPos.empty()) {
          const auto& lp = ds->LabPos[n];
          switch (lp.mode) {
          case 1:
            add3f(lp.offset, v, v);
            copy3f(lab_pos, v + 3);
            break;
          default:
            copy3f(lab_pos, v + 3);
            break;
          }
        } else {
          copy3f(lab_pos, v + 3);
        }
        if(rp) {
          rp->index = n;        /* label index */
          rp->bond = cPickableLabel;    /* label indicator */
          rp++;
        }
        n++;
      }
    }

    if(ok && ds->NDihedralIndex) {

      float d12[3], d32[3], d43[3], n32[3];
      float p12[3], p43[3], np12[3], np43[3];

      float a32[3];
      float l1, l2;
      float radius;
      float dihedral_size =
        SettingGet_f(G, NULL, ds->Obj->Setting.get(), cSetting_dihedral_size);
      float dihedral_label_position = SettingGet_f(G, NULL, ds->Obj->Setting.get(),
                                                   cSetting_dihedral_label_position);

      float *v4;
      float avg[3];
      int digits = SettingGet_i(G, NULL, ds->Obj->Setting.get(),
                                cSetting_label_dihedral_digits);
      WordType format;
      if(digits < 0)
        digits = default_digits;
      if(digits > 10)
        digits = 10;
      sprintf(format, "%c0.%df", '%', digits);

      for(a = 0; ok && a < ds->NDihedralIndex; a = a + 6) {
        v1 = ds->DihedralCoord + 3 * a;
        v2 = v1 + 3;
        v3 = v1 + 6;
        v4 = v1 + 9;

        subtract3f(v1, v2, d12);
        subtract3f(v3, v2, d32);
        subtract3f(v4, v3, d43);

        normalize23f(d32, n32);

        remove_component3f(d12, n32, p12);
        remove_component3f(d43, n32, p43);

        average3f(v2, v3, a32);

        l1 = (float) length3f(p12);
        l2 = (float) length3f(p43);

        if(l1 > l2)
          radius = l2;
        else
          radius = l1;
        radius *= dihedral_size * dihedral_label_position;

        normalize23f(p12, np12);
        normalize23f(p43, np43);

        average3f(np12, np43, avg);

        normalize3f(avg);
        if((avg[0] == 0.0F) && (avg[1] == 0.0F) && (avg[2] == 0.0F))
          copy3f(np12, avg);

        scale3f(avg, radius, avg);
        add3f(a32, avg, avg);

        di = (float) (180.0F * get_dihedral3f(v1, v2, v3, v4) / PI);
        sprintf(buffer, format, di);

        VLACheck(I->V, float, 6 * n + 5);
	CHECKOK(ok, I->V);
	if (ok)
	  VLACheck(I->L, DistLabel, n);
	CHECKOK(ok, I->L);
        v = I->V + 6 * n;
        UtilNCopy(I->L[n], buffer, sizeof(DistLabel));
        copy3f(avg, v);
        std::copy_n(v, 3, lc->data());
        lc++;
        if(!ds->LabPos.empty()) {
          const auto& lp = ds->LabPos[n];
          switch (lp.mode) {
          case 1:
            add3f(lp.offset, v, v);
            copy3f(lab_pos, v + 3);
            break;
          default:
            copy3f(lab_pos, v + 3);
            break;
          }
        } else {
          copy3f(lab_pos, v + 3);
        }
        copy3f(lab_pos, v + 3);

        if(rp) {
          rp->index = n;        /* label index */
          rp->bond = cPickableLabel;    /* label indicator */
          rp++;
        }
        n++;
      }
    }
  }

  I->N = n;

  if(ok && rp) {
    I->P = ReallocForSure(I->P, Pickable, (rp - I->P));
    CHECKOK(ok, I->P);
    I->P[0].index = I->N;     /* unnec? */
  }
  if (!ok){
    delete I;
    I = NULL;
  }
  return (Rep *) I;
}
