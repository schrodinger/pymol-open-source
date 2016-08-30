
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

typedef struct RepRibbon {
  Rep R;
  float *V;
  float ribbon_width;
  float radius;
  int N;
  int NS;
  int NP;
  CGO *shaderCGO;
} RepRibbon;

#include"ObjectMolecule.h"

void RepRibbonFree(RepRibbon * I);

void RepRibbonInit(void)
{
}

void RepRibbonFree(RepRibbon * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  FreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepRibbonRender(RepRibbon * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  Pickable *p;
  int i, j, ip;
  int last;
  float line_width = SceneGetDynamicLineWidth(info, I->ribbon_width);
  int ok = true;

  /* 

     v[0] = index1
     v[1-3] = color1
     v[4-6] = vertex1
     v[7] = index2
     v[8-10] = color2
     v[11-13] = vertex2
     v[14] = radius

   */

  if(ray) {

    float radius;
    /*    float alpha = SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_ribbon_transparency); */
    float alpha = SettingGet_f(G, NULL, I->R.obj->Setting, cSetting_ribbon_transparency);
    alpha = 1.0F - alpha;
    if(fabs(alpha-1.0) < R_SMALL4)
      alpha = 1.0F;

    if(I->radius == 0.0F) {
      radius = ray->PixelRadius * line_width / 2.0F;
    } else {
      radius = I->radius;
    }

    PRINTFD(G, FB_RepRibbon)
      " RepRibbonRender: rendering raytracable...\n" ENDFD;

    if(c > 0) {
      while(ok && c--) {
	ray->transparentf(1.0F - alpha);
        ok &= ray->sausage3fv(v + 4, v + 11, radius, v + 1, v + 8);
        v += 18;
      }
    }
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {

      PRINTFD(G, FB_RepRibbon)
        " RepRibbonRender: rendering pickable...\n" ENDFD;

      if(c) {
        i = (*pick)->src.index;
        p = I->R.P;
        last = -1;
#ifdef PURE_OPENGL_ES_2
#else
	SceneSetupGLPicking(G);
        glBegin(GL_LINES);
        while(c--) {
          ip = (int) *(v);
          if(ip != last) {
            i++;
            last = ip;
            if(!(*pick)[0].src.bond) {
              /* pass 1 - low order bits */

              glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8), (uchar) ((i & 0xF00) >> 4));     /* we're encoding the index into the color */
              VLACheck((*pick), Picking, i);
              (*pick)[i].src = p[ip];   /* copy object and atom info */
              (*pick)[i].context = I->R.context;
            } else {
              /* pass 2 - high order bits */
              j = i >> 12;
              glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                         (uchar) ((j & 0xF00) >> 4));
            }
          }
          glVertex3fv(v + 4);
          ip = (int) *(v + 7);
          if(ip != last) {
            glVertex3fv(v + 15);        /* switch colors at midpoint */
            glVertex3fv(v + 15);
            i++;
            last = ip;
            if(!(*pick)[0].src.bond) {
              /* pass 1 - low order bits */

              glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8), (uchar) ((i & 0xF00) >> 4));     /* we're encoding the index into the color */
              VLACheck((*pick), Picking, i);
              (*pick)[i].src = p[ip];   /* copy object and atom info */
              (*pick)[i].context = I->R.context;
            } else {
              /* pass 2 - high order bits */
              j = i >> 12;
              glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                         (uchar) ((j & 0xF00) >> 4));
            }
          }
          glVertex3fv(v + 11);
          v += 18;
        }
        glEnd();
#endif
        (*pick)[0].src.index = i;       /* pass the count */
      }
    } else {
      short use_shader, generate_shader_cgo = 0;
      int ribbon_smooth;
      short ribbon_as_cylinders ;
      float alpha = SettingGet_f(G, NULL, I->R.obj->Setting, cSetting_ribbon_transparency);

      use_shader = SettingGetGlobal_b(G, cSetting_ribbon_use_shader) & 
                   SettingGetGlobal_b(G, cSetting_use_shaders);
      ribbon_as_cylinders = SettingGetGlobal_b(G, cSetting_render_as_cylinders) && SettingGetGlobal_b(G, cSetting_ribbon_as_cylinders);
      if (!use_shader && I->shaderCGO){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }
      if (I->shaderCGO && (ribbon_as_cylinders ^ I->shaderCGO->has_draw_cylinder_buffers)){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }

      if (use_shader){
	if (!I->shaderCGO){
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok)
	    I->shaderCGO->use_shader = true;
	  generate_shader_cgo = 1;
	} else if (ok) {
	  CShaderPrg *shaderPrg;
	  if (ribbon_as_cylinders){
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if(I->radius == 0.0F) {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", I->radius);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }
	  if (!shaderPrg) return;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  
	  CShaderPrg_Disable(shaderPrg);
	  return;
	}
      }

      alpha = 1.0F - alpha;
      if(fabs(alpha-1.0) < R_SMALL4)
	alpha = 1.0F;

      ribbon_smooth = SettingGet_i(G, NULL, I->R.obj->Setting, cSetting_ribbon_smooth);
      if(!ribbon_smooth)
        glDisable(GL_LINE_SMOOTH);

      if (generate_shader_cgo){
	ok &= CGOLinewidthSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_RIBBON);
	if (ok)
	  ok &= CGOResetNormal(I->shaderCGO, true);
      } else {
	if(info->width_scale_flag)
	  glLineWidth(line_width * info->width_scale);
	else
	  glLineWidth(line_width);
	SceneResetNormal(G, true);
      }
      PRINTFD(G, FB_RepRibbon)
	" RepRibbonRender: rendering GL...\n" ENDFD;
      
      if(ok && c) {
	int first = true;
	if (generate_shader_cgo){
	  if (ribbon_as_cylinders){
	    float *origin = NULL, *vertex = NULL, *color = NULL, *color2 = NULL, axis[3];
	    while(ok && c--) {
	      if(first) {
		/* glColor3fv(v + 1); */
		ok &= CGOAlpha(I->shaderCGO, alpha);
		if (ok)
		  ok &= CGOColorv(I->shaderCGO, v+1);
		color2 = (v + 1);
		vertex = (v + 4);
		first = false;
	      } else if((v[4] != v[-7]) || (v[5] != v[-6]) || (v[6] != v[-5])) {
		ok &= CGOAlpha(I->shaderCGO, alpha);
		if (ok)
		  ok &= CGOColorv(I->shaderCGO, v+1);
		origin = color = NULL;
		color2 = (v+1);
		vertex = (v+4);
	      }
	      /* glColor3fv(v + 8); */
	      origin = vertex;
	      color = color2;
	      color2 = (v+8);
	      vertex = (v+11);
	      axis[0] = vertex[0] - origin[0];
	      axis[1] = vertex[1] - origin[1];
	      axis[2] = vertex[2] - origin[2];
	      if (ok){
		if (*(color) != *(color2) || *(color+1) != *(color2+1) || *(color+2) != *(color2+2)){
		  ok &= CGOShaderCylinder2ndColor(I->shaderCGO, origin, axis, 1.f, 15, color2);
		} else {
		  ok &= CGOShaderCylinder(I->shaderCGO, origin, axis, 1.f, 15);
		}
	      }
	      if (ok)
		ok &= CGOAlpha(I->shaderCGO, alpha);
	      if (ok)
		ok &= CGOColorv(I->shaderCGO, v+8);
	      v += 18;
	    }
	  } else if (ok) {
	    ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
	    while(ok && c--) {
	      if(first) {
		/* glColor3fv(v + 1); */
		ok &= CGOAlpha(I->shaderCGO, alpha);
		if (ok)
		  ok &= CGOColorv(I->shaderCGO, v+1);
		if (ok)
		  ok &= CGOVertexv(I->shaderCGO, v + 4);
		first = false;
	      } else if((v[4] != v[-7]) || (v[5] != v[-6]) || (v[6] != v[-5])) {
		ok &= CGOEnd(I->shaderCGO);
		if (ok)
		  ok &= CGOBegin(I->shaderCGO, GL_LINE_STRIP);
		/* glColor3fv(v + 1); */
		if (ok)
		  ok &= CGOAlpha(I->shaderCGO, alpha);
		if (ok)
		  ok &= CGOColorv(I->shaderCGO, v+1);
		if (ok)
		  ok &= CGOVertexv(I->shaderCGO, v + 4);
	      }
	      /* glColor3fv(v + 8); */
	      if (ok)
		ok &= CGOAlpha(I->shaderCGO, alpha);
	      if (ok)
		ok &= CGOColorv(I->shaderCGO, v+8);
	      if (ok)
		ok &= CGOVertexv(I->shaderCGO, v + 11);
	      v += 18;
	    }
	    if (ok)
	      ok &= CGOEnd(I->shaderCGO);
	  }
	} else {
	  if(!info->line_lighting)
	    glDisable(GL_LIGHTING);
	  glBegin(GL_LINE_STRIP);
	  while(c--) {
	    if(first) {
	      /* glColor3fv(v + 1); */
	      glColor4f( (v+1)[0], (v+1)[1], (v+1)[2], alpha);
	      glVertex3fv(v + 4);
	      first = false;
	    } else if((v[4] != v[-7]) || (v[5] != v[-6]) || (v[6] != v[-5])) {
	      glEnd();
	      glBegin(GL_LINE_STRIP);
	      /* glColor3fv(v + 1); */
	      glColor4f( (v+1)[0], (v+1)[1], (v+1)[2], alpha);
	      glVertex3fv(v + 4);
	    }
	    /* glColor3fv(v + 8); */
	    glColor4f( (v+8)[0], (v+8)[1], (v+8)[2], alpha);
	    glVertex3fv(v + 11);
	    v += 18;
	  }
	  glEnd();
	  glEnable(GL_LIGHTING);
	}
      }
      if(SettingGetGlobal_b(G, cSetting_line_smooth))
        glEnable(GL_LINE_SMOOTH);

      if (use_shader) {
	if (ok && generate_shader_cgo){
	  CGO *convertcgo = NULL;
	  ok &= CGOStop(I->shaderCGO);
	  if (ok)
	    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);    
	  CHECKOK(ok, convertcgo);
	  CGOFree(I->shaderCGO);    
	  I->shaderCGO = convertcgo;
	  convertcgo = NULL;
	  if (ok){
	    if (ribbon_as_cylinders){
	      convertcgo = CGOOptimizeGLSLCylindersToVBOIndexed(I->shaderCGO, 0);
	    } else {
	      convertcgo = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0);
	    }
	    CHECKOK(ok, convertcgo);
	  }
	  if (convertcgo){
	    CGOFree(I->shaderCGO);
	    I->shaderCGO = convertcgo;
	    convertcgo = NULL;
	  }
	}
	
	if (ok){
	  CShaderPrg *shaderPrg;
	  if (ribbon_as_cylinders){
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if(I->radius == 0.0F) {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", I->radius);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }
	  if (!shaderPrg) return;

	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  
	  CShaderPrg_Disable(shaderPrg);
	}
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
    I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
    I->R.cs->Active[cRepRibbon] = false;
  }
}

Rep *RepRibbonNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int a, b, a1, a2, c1, c2, *i, *s, *at, *seg, nAt, *atp;
  float *v, *v0, *v1, *v2, *v3;
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

  Pickable *rp = NULL;
  OOAlloc(G, RepRibbon);

  obj = cs->Obj;

  RepInit(G, &I->R);
  power_a = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_power);
  power_b = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_power_b);
  throw_ = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_throw);
  int trace_ostate = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_trace_atoms);
  trace_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_trace_atoms_mode);
  na_mode =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_nucleic_acid_mode);

  ribbon_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_color);

  sampling = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_sampling);
  if(sampling < 1)
    sampling = 1;
  I->radius = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_radius);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepRibbonRender;
  I->R.fFree = (void (*)(struct Rep *)) RepRibbonFree;
  I->R.fRecolor = NULL;
  I->R.obj = (CObject *) obj;
  I->R.cs = cs;
  I->ribbon_width = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_width);
  I->R.context.object = (void *) obj;
  I->R.context.state = state;

  /* find all of the CA points */

  at = Alloc(int, cs->NAtIndex * 2);
  pv = Alloc(float, cs->NAtIndex * 6);
  seg = Alloc(int, cs->NAtIndex * 2);

  i = at;
  v = pv;
  s = seg;

  nAt = 0;
  nSeg = 0;
  a2 = -1;
  for(a1 = 0; a1 < cs->NAtIndex; a1++) {
    if(obj->DiscreteFlag) {
      if(cs == obj->DiscreteCSet[a1])
        a = obj->DiscreteAtmToIdx[a1];
      else
        a = -1;
    } else
      a = cs->AtmToIdx[a1];
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
            /*                                                if((abs(obj->AtomInfo[a1].resv-obj->AtomInfo[a2].resv)>1)||
               (obj->AtomInfo[a1].chain[0]!=obj->AtomInfo[a2].chain[0])||
               (!WordMatch(G,obj->AtomInfo[a1].segi,obj->AtomInfo[a2].segi,1))) */
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
          v1 = cs->Coord + 3 * a;
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
                v1 = cs->Coord + 3 * trailing_O3p_a;
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
              v1 = cs->Coord + 3 * leading_O5p_a;
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
          v1 = cs->Coord + 3 * a;
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
    v1 = cs->Coord + 3 * trailing_O3p_a;
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

    dv = Alloc(float, nAt * 6);
    nv = Alloc(float, nAt * 6);
    dl = Alloc(float, nAt * 2);
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

    tv = Alloc(float, nAt * 6 + 6);
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

  if(nAt) {
    I->R.P = Alloc(Pickable, 2 * nAt + 2);
    ErrChkPtr(G, I->R.P);
    I->R.P[0].index = nAt;
    rp = I->R.P + 1;            /* skip first record! */
  }

  I->V = (float *) mmalloc(sizeof(float) * 2 * cs->NIndex * 18 * sampling);
  ErrChkPtr(G, I->V);

  I->shaderCGO = 0;
  I->N = 0;
  I->NP = 0;
  v = I->V;
  if(nAt) {
    v1 = pv;                    /* points */
    v2 = tv;                    /* tangents */
    v3 = dv;                    /* direction vector */
    d = dl;
    s = seg;
    atp = at;
    rp->index = cs->IdxToAtm[*atp];
    if(obj->AtomInfo[cs->IdxToAtm[*atp]].masked)
      rp->index = -1;
    rp->bond = -1;
    I->NP++;
    rp++;
    for(a = 0; a < (nAt - 1); a++) {

      rp->index = cs->IdxToAtm[*(atp + 1)];     /* store pickable for n+2 */
      if(obj->AtomInfo[cs->IdxToAtm[*atp]].masked)
        rp->index = -1;
      rp->bond = -1;
      rp++;
      I->NP++;

      PRINTFD(G, FB_RepRibbon)
        " RepRibbon: seg %d *s %d , *(s+1) %d\n", a, *s, *(s + 1)
        ENDFD;

      if(*s == *(s + 1)) {
        AtomInfoType *ai1 = obj->AtomInfo + cs->IdxToAtm[*atp];
        AtomInfoType *ai2 = obj->AtomInfo + cs->IdxToAtm[*(atp + 1)];

        c1 = AtomSettingGetWD(G, ai1, cSetting_ribbon_color, ribbon_color);
        c2 = AtomSettingGetWD(G, ai2, cSetting_ribbon_color, ribbon_color);

        if (c1 < 0) c1 = ai1->color;
        if (c2 < 0) c2 = ai2->color;

        dev = throw_ * (*d);

        for(b = 0; b < sampling; b++) { /* needs optimization */

          f0 = ((float) b) / sampling;  /* fraction of completion */
          f0 = smooth(f0, power_a);     /* bias sampling towards the center of the curve */

          if(f0 < 0.5) {
            v0 = ColorGet(G, c1);
          } else {
            v0 = ColorGet(G, c2);
          }

          /* store index */
          if(f0 < 0.5F)
            *(v++) = (float) I->NP - 1;
          else
            *(v++) = (float) I->NP;

          /* store colors */

          *(v++) = *(v0++);
          *(v++) = *(v0++);
          *(v++) = *(v0);

          /* start of line/cylinder */

          f1 = 1.0F - f0;
          f2 = smooth(f0, power_b);
          f3 = smooth(f1, power_b);
          f4 = dev * f2 * f3;   /* displacement magnitude */

          *(v++) = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);

          *(v++) = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);

          *(v++) = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);

          f0 = ((float) b + 1) / sampling;
          f0 = smooth(f0, power_a);

          if(f0 < 0.5) {
            v0 = ColorGet(G, c1);
          } else {
            v0 = ColorGet(G, c2);
          }

          /* store index */
          if(f0 < 0.5)
            *(v++) = (float) I->NP - 1;
          else
            *(v++) = (float) I->NP;

          /* store colors */

          *(v++) = *(v0++);
          *(v++) = *(v0++);
          *(v++) = *(v0);

          /* end of line/cylinder */

          f1 = 1.0F - f0;
          f2 = smooth(f0, power_b);
          f3 = smooth(f1, power_b);
          f4 = dev * f2 * f3;   /* displacement magnitude */

          *(v++) = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);

          *(v++) = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);

          *(v++) = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);

          v++;                  /* radius no longer stored here... */

          average3f(v - 4, v - 11, v);

          v += 3;

          I->N++;
        }

      }
      v1 += 3;
      v2 += 3;
      v3 += 3;
      d++;
      atp += 1;
      s++;
    }

    FreeP(dv);
    FreeP(dl);
    FreeP(tv);
    FreeP(nv);
  }

  FreeP(at);
  FreeP(seg);
  FreeP(pv);

  if(I->N)
    I->V = ReallocForSure(I->V, float, (v - I->V));
  else
    I->V = ReallocForSure(I->V, float, 1);

  return (Rep *) I;
}

void RepRibbonRenderImmediate(CoordSet * cs, RenderInfo * info)
{
#ifndef PURE_OPENGL_ES_2
  /* performance optimized to provide a simple C-alpha trace -- no smoothing */

  PyMOLGlobals *G = cs->State.G;
  if(info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    ObjectMolecule *obj = cs->Obj;
    int active = false;
    int nAtIndex = cs->NAtIndex;
    int a;
    AtomInfoType *obj_AtomInfo = obj->AtomInfo;
    AtomInfoType *ai, *last_ai = NULL;
    int trace, trace_ostate =
      SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_trace_atoms);
    int trace_mode =
      SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_trace_atoms_mode);
    int na_mode =
      SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_nucleic_acid_mode);
    float ribbon_width =
      SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_ribbon_width);
    int a1, a2 = -1;
    int color, last_color = -9;

    glLineWidth(ribbon_width);
    SceneResetNormal(G, true);
    if(!info->line_lighting)
      glDisable(GL_LIGHTING);
    glBegin(GL_LINE_STRIP);
    for(a1 = 0; a1 < nAtIndex; a1++) {
      if(obj->DiscreteFlag) {
        if(cs == obj->DiscreteCSet[a1])
          a = obj->DiscreteAtmToIdx[a1];
        else
          a = -1;
      } else
        a = cs->AtmToIdx[a1];
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
            glVertex3fv(cs->Coord + 3 * a);
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
            glVertex3fv(cs->Coord + 3 * a);
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
