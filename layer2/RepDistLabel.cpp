
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

typedef struct RepDistLabel {
  Rep R;
  float *V;
  int N;
  DistLabel *L;
  CObject *Obj;
  DistSet *ds;
  int OutlineColor;
  CGO *shaderCGO;
  int texture_font_size;
} RepDistLabel;

#define SHADERCGO I->shaderCGO

#include"ObjectDist.h"

static
void RepDistLabelFree(RepDistLabel * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
  }
  VLAFreeP(I->V);
  VLAFreeP(I->L);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepDistLabelRender(RepDistLabel * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  auto pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  DistLabel *l = I->L;
  int n = 0;
  int color;
  int font_id = SettingGet_i(G, NULL, I->Obj->Setting, cSetting_label_font_id);
  float font_size = SettingGet_f(G, NULL, I->Obj->Setting, cSetting_label_size);
  int float_text = SettingGet_i(G, NULL, I->Obj->Setting, cSetting_float_labels);
  int ok = true;
  short use_shader = SettingGetGlobal_b(G, cSetting_use_shaders);
  if (I->R.MaxInvalid >= cRepInvRep)
    return;
  font_id = SettingCheckFontID(G, NULL, I->Obj->Setting, font_id);

  if (I->shaderCGO && font_size < 0.f){
    int size;
    if (InvalidateShaderCGOIfTextureNeedsUpdate(G, font_size, I->texture_font_size, &size)){
      CGOFree(I->shaderCGO);
      I->texture_font_size = size;
    }
  }

  if(ray) {
    TextSetOutlineColor(G, I->OutlineColor);
    color = SettingGet_color(G, NULL, I->Obj->Setting, cSetting_label_color);

    if((color >= 0) || (color == cColorFront) || (color == cColorBack))
      TextSetColor(G, ColorGet(G, color));
    else
      TextSetColor(G, ColorGet(G, I->Obj->Color));

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
	CGORenderGLPicking(I->shaderCGO, info, &I->R.context, NULL, NULL);
        if(float_text)
          glEnable(GL_DEPTH_TEST);
	return;
    } else {
	Pickable *p = I->R.P;
	unsigned int i;
        TextSetIsPicking(G, true);
	SceneSetupGLPicking(G);
	if(c) {
	  if(float_text)
	    glDisable(GL_DEPTH_TEST);
	  
	  i = pick->begin()->src.index;
	  while(c--) {
	    if(*l) {
	      TextSetPos(G, v);
              p++;
              AssignNewPickColor(NULL, i, pick, &I->R.context, TextGetColorUChar4uv(G), p->index, p->bond);
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
	  (*pick)[0].src.index = i;       /* pass the count */
	}
        TextSetIsPicking(G, false);
      }
    } else {
	Pickable *p = I->R.P;
	
        if (use_shader){
	if (!I->shaderCGO){
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok){
	    I->shaderCGO->use_shader = true;
	  }
	} else {
	    info->texture_font_size = I->texture_font_size;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  return;
	}
	} else if (I->shaderCGO){
	  CGOFree(I->shaderCGO);
	  if(float_text)
	    glDisable(GL_DEPTH_TEST);
	}

      TextSetOutlineColor(G, I->OutlineColor);
      color = SettingGet_color(G, NULL, I->Obj->Setting, cSetting_label_color);

      if((color >= 0) || (color == cColorFront) || (color == cColorBack))
        TextSetColor(G, ColorGet(G, color));
      else
        TextSetColor(G, ColorGet(G, I->Obj->Color));
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
	  RepDistLabelRender(I, info);
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
    RepDistLabelFree(I);
  }
}

Rep *RepDistLabelNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->State.G;
  int a;
  int n = 0;
  float *v, *v1, *v2, *v3, d[3], di;
  char buffer[255];
  const float *lab_pos =
    SettingGet_3fv(G, NULL, ds->Obj->Obj.Setting, cSetting_label_position);
  int default_digits =
    SettingGet_i(G, NULL, ds->Obj->Obj.Setting, cSetting_label_digits);
  Pickable *rp = NULL;
  int ok = true;

  OOAlloc(G, RepDistLabel);
  CHECKOK(ok, I);

  if(!ok || !(ds->NIndex || ds->NAngleIndex || ds->NDihedralIndex)) {
    ds->NLabel = 0;
    VLAFreeP(ds->LabCoord);
    VLAFreeP(ds->LabPos);
    OOFreeP(I);
    return (NULL);
  }

  if(default_digits < 0)
    default_digits = 0;
  if(default_digits > 10)
    default_digits = 10;

  RepInit(G, &I->R);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepDistLabelRender;
  I->R.fFree = (void (*)(struct Rep *)) RepDistLabelFree;
  I->R.fRecolor = NULL;

  I->N = 0;
  I->V = NULL;
  I->R.P = NULL;
  I->Obj = (CObject *) ds->Obj;
  I->R.obj = I->Obj;
  I->ds = ds;
  I->R.context.object = (void *) ds->Obj;
  I->R.context.state = state;

  I->shaderCGO = NULL;
  I->texture_font_size = 0;

  I->OutlineColor =
    SettingGet_i(G, NULL, I->Obj->Setting, cSetting_label_outline_color);

  if(ds->NIndex || ds->NAngleIndex || ds->NDihedralIndex) {
    float *lc;
    ds->NLabel = (ds->NIndex / 2 + ds->NAngleIndex / 5 + ds->NDihedralIndex / 6);

    if(!ds->LabCoord) {         /* store label coordinates */
      ds->LabCoord = VLAlloc(float, 3 * ds->NLabel);
    } else {
      VLACheck(ds->LabCoord, float, 3 * ds->NLabel);
    }
    CHECKOK(ok, ds->LabCoord);
      
    if(ok && ds->LabPos) {            /* make sure this VLA covers all labels */
      VLACheck(ds->LabPos, LabPosType, ds->NLabel);
      CHECKOK(ok, ds->LabPos);
    }

    if(ok && SettingGet_b(G, NULL, ds->Obj->Obj.Setting, cSetting_pickable)) {
      I->R.P = pymol::malloc<Pickable>(ds->NLabel + 1);
      CHECKOK(ok, I->R.P);
      if (ok)
	rp = I->R.P + 1;          /* skip first record! */
    }

    if (ok)
      I->V = VLAlloc(float, 3 * (ds->NIndex / 2 + ds->NAngleIndex / 5) + 1);
    CHECKOK(ok, I->V);
    if (ok)
      I->L = VLAlloc(DistLabel, (ds->NIndex / 2 + ds->NAngleIndex / 5) + 1);
    CHECKOK(ok, I->L);

    n = 0;
    lc = ds->LabCoord;

    if(ds->NIndex) {
      int digits = SettingGet_i(G, NULL, ds->Obj->Obj.Setting,
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
        *(lc++) = v[0];
        *(lc++) = v[1];
        *(lc++) = v[2];
        if(ds->LabPos) {
          LabPosType *lp = ds->LabPos + n;
          switch (lp->mode) {
          case 1:
            add3f(lp->offset, v, v);
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
      int digits = SettingGet_i(G, NULL, ds->Obj->Obj.Setting,
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
          SettingGet_f(G, NULL, ds->Obj->Obj.Setting,
                       cSetting_angle_size) * SettingGet_f(G, NULL,
                                                           ds->Obj->Obj.Setting,
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
        *(lc++) = v[0];
        *(lc++) = v[1];
        *(lc++) = v[2];
        if(ds->LabPos) {
          LabPosType *lp = ds->LabPos + n;
          switch (lp->mode) {
          case 1:
            add3f(lp->offset, v, v);
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
        SettingGet_f(G, NULL, ds->Obj->Obj.Setting, cSetting_dihedral_size);
      float dihedral_label_position = SettingGet_f(G, NULL, ds->Obj->Obj.Setting,
                                                   cSetting_dihedral_label_position);

      float *v4;
      float avg[3];
      int digits = SettingGet_i(G, NULL, ds->Obj->Obj.Setting,
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
        *(lc++) = v[0];
        *(lc++) = v[1];
        *(lc++) = v[2];
        if(ds->LabPos) {
          LabPosType *lp = ds->LabPos + n;
          switch (lp->mode) {
          case 1:
            add3f(lp->offset, v, v);
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
    I->R.P = ReallocForSure(I->R.P, Pickable, (rp - I->R.P));
    CHECKOK(ok, I->R.P);
    I->R.P[0].index = I->N;     /* unnec? */
  }
  if (!ok){
    RepDistLabelFree(I);
    I = NULL;
  }
  return (Rep *) I;
}
