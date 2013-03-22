
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
#include"RepAngle.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"
#include"ShaderMgr.h"
#include"CGO.h"
#include"CoordSet.h"

typedef struct RepAngle {
  Rep R;
  float *V;
  int N;
  CObject *Obj;
  DistSet *ds;
  float linewidth, radius;
  CGO *shaderCGO;
} RepAngle;

#include"ObjectDist.h"

void RepAngleFree(RepAngle * I);

void RepAngleFree(RepAngle * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  VLAFreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepAngleRender(RepAngle * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  float *vc;
  int round_ends;
  float line_width;
  int ok = true;
  int color =
    SettingGet_color(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_angle_color);
  I->linewidth = line_width = 
    SettingGet_f(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_width);
  I->radius =
    SettingGet_f(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_radius);
  round_ends =
    SettingGet_b(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_round_ends);
  line_width = SceneGetDynamicLineWidth(info, line_width);

  if(ray) {

    float radius;

    if(I->radius == 0.0F) {
      radius = ray->PixelRadius * line_width / 2.0F;
    } else {
      radius = I->radius;
    }

    if(color < 0)
      color = I->Obj->Color;
    vc = ColorGet(G, color);
    v = I->V;
    c = I->N;

    while(ok && c > 0) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]); */
      if(round_ends) {
        ok &= ray->fSausage3fv(ray, v, v + 3, radius, vc, vc);
      } else {
        ok &= ray->fCustomCylinder3fv(ray, v, v + 3, radius, vc, vc, cCylCapFlat, cCylCapFlat);
      }
      v += 6;
      c -= 2;
    }

  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      short use_shader, generate_shader_cgo = 0, use_display_lists = 0, dash_as_cylinders = 0;

      use_shader = (int) SettingGet(G, cSetting_dash_use_shader) & 
                           (int) SettingGet(G, cSetting_use_shaders);
      use_display_lists = (int) SettingGet(G, cSetting_use_display_lists);
      dash_as_cylinders = (int) SettingGet(G, cSetting_render_as_cylinders) && SettingGet(G, cSetting_dash_as_cylinders);

      if (!use_shader && I->shaderCGO){
	CGOFree(I->shaderCGO);
	I->shaderCGO = 0;
      }
#ifdef _PYMOL_GL_CALLLISTS
      if(use_display_lists && I->R.displayList) {
	glCallList(I->R.displayList);
	return;
      }
#endif

      if (use_shader){
	if (!I->shaderCGO){
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok)
	    I->shaderCGO->use_shader = true;
	  generate_shader_cgo = 1;
	} else {
	  CShaderPrg *shaderPrg;
	  if (dash_as_cylinders){
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if(I->radius == 0.0F) {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", I->radius);
	    }
	    if (!round_ends){
	      CShaderPrg_Set1f(shaderPrg, "no_flat_caps", 0.f);
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
#ifdef _PYMOL_GL_CALLLISTS
      if(use_display_lists) {
	if(!I->R.displayList) {
	  I->R.displayList = glGenLists(1);
	  if(I->R.displayList) {
	    glNewList(I->R.displayList, GL_COMPILE_AND_EXECUTE);
	  }
	}
      }
#else
      (void) use_display_lists;
#endif

      if (generate_shader_cgo){
	if (ok)
	  ok &= CGOLinewidthSpecial(I->shaderCGO, LINEWIDTH_DYNAMIC_WITH_SCALE_DASH);
	if (ok)
	  ok &= CGOResetNormal(I->shaderCGO, true);
      } else {
	if(info->width_scale_flag) {
	  glLineWidth(line_width * info->width_scale);
	} else {
	  glLineWidth(line_width);
	}
        SceneResetNormal(G, true);
      }

      if (generate_shader_cgo){
	if (ok){
	  if(color >= 0){
	    ok &= CGOColorv(I->shaderCGO, ColorGet(G, color));
	  } else if (I->Obj && I->Obj->Color >= 0){
	    ok &= CGOColorv(I->shaderCGO, ColorGet(G, I->Obj->Color));
	  }
	}
	v = I->V;
	c = I->N;
	if (dash_as_cylinders){
	  float *origin = NULL, axis[3];
	  while(ok && c > 0) {
	    origin = v;
	    v += 3;
	    axis[0] = v[0] - origin[0];
	    axis[1] = v[1] - origin[1];
	    axis[2] = v[2] - origin[2];
	    v += 3;
	    ok &= CGOShaderCylinder(I->shaderCGO, origin, axis, 1.f, 15);
	    c -= 2;
	  }
	} else if (ok) {
	  ok &= CGOBegin(I->shaderCGO, GL_LINES);
	  while(ok && c > 0) {
	    ok &= CGOVertexv(I->shaderCGO, v);
	    v += 3;
	    if (ok)
	      ok &= CGOVertexv(I->shaderCGO, v);
	    v += 3;
	    c -= 2;
	  }
	  if (ok)
	    ok &= CGOEnd(I->shaderCGO);
	}
      } else {
	if(color >= 0){
	  glColor3fv(ColorGet(G, color));
	}
        v = I->V;
        c = I->N;

        if(!info->line_lighting)
          glDisable(GL_LIGHTING);
#ifdef _PYMOL_GL_DRAWARRAYS
	{
	  int pl, cinit;
	  ALLOCATE_ARRAY(GLfloat,lineVerts,c)
	  pl = 0;
	  cinit = c;
	  while(c > 0) {
	    lineVerts[pl++] = v[0]; lineVerts[pl++] = v[1]; lineVerts[pl++] = v[2];
	    v += 3;
	    lineVerts[pl++] = v[0]; lineVerts[pl++] = v[1]; lineVerts[pl++] = v[2];
	    v += 3;
	    c -= 2;
	  }
	  glEnableClientState(GL_VERTEX_ARRAY);
	  glVertexPointer(3, GL_FLOAT, 0, lineVerts);
	  glDrawArrays(GL_LINE_LOOP, 0, cinit);
	  glDisableClientState(GL_VERTEX_ARRAY);
	  DEALLOCATE_ARRAY(lineVerts)
	}
#else
        glBegin(GL_LINES);
        while(c > 0) {
          glVertex3fv(v);
          v += 3;
          glVertex3fv(v);
          v += 3;
          c -= 2;
        }
        glEnd();
#endif
        glEnable(GL_LIGHTING);
      }
      if (use_shader) {
	if (generate_shader_cgo){
	  CGO *convertcgo = NULL;
	  if (ok)
	    ok &= CGOStop(I->shaderCGO);
#ifdef _PYMOL_CGO_DRAWARRAYS
	  if (ok)
	    convertcgo = CGOCombineBeginEnd(I->shaderCGO, 0);    
	  CHECKOK(ok, convertcgo);
	  CGOFree(I->shaderCGO);    
	  I->shaderCGO = convertcgo;
	  convertcgo = NULL;
#else
	  (void)convertcgo;
#endif
#ifdef _PYMOL_CGO_DRAWBUFFERS
	  if (ok){
	    if (dash_as_cylinders){
	      convertcgo = CGOOptimizeGLSLCylindersToVBOIndexed(I->shaderCGO, 0);
	    } else {
	      convertcgo = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0);
	    }
	  }
	  CHECKOK(ok, convertcgo);
	  if (convertcgo){
	    CGOFree(I->shaderCGO);
	    I->shaderCGO = convertcgo;
	    convertcgo = NULL;
	  }
#else
	  (void)convertcgo;
#endif
	}
	
	if (ok) {
	  CShaderPrg *shaderPrg;
	  if (dash_as_cylinders){
	    float pixel_scale_value = SettingGetGlobal_f(G, cSetting_ray_pixel_scale);
	    if(pixel_scale_value < 0)
	      pixel_scale_value = 1.0F;
	    shaderPrg = CShaderPrg_Enable_CylinderShader(G);
	    if(I->radius == 0.0F) {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", info->vertex_scale * pixel_scale_value * line_width/ 2.f);
	    } else {
	      CShaderPrg_Set1f(shaderPrg, "uni_radius", I->radius);
	    }
	    if (!round_ends){
	      CShaderPrg_Set1f(shaderPrg, "no_flat_caps", 0.f);
	    }
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_SetLightingEnabled(shaderPrg, 0);
	  }	 

	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  
	  CShaderPrg_Disable(shaderPrg);
	}
      }
#ifdef _PYMOL_GL_CALLLISTS
      if (use_display_lists && I->R.displayList){
	glEndList();
	glCallList(I->R.displayList);      
      }
#endif
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
    I->ds->Rep[cRepAngle] = NULL;
    RepAngleFree(I);
  }
}

Rep *RepAngleNew(DistSet * ds, int state)
{
  PyMOLGlobals *G = ds->State.G;
  int a;
  int n = 0;
  float *v, *v1, *v2, *v3, *v4, d1[3], d2[3], d3[3], n1[3], n3[3], l1, l2, x[3], y[3];
  float length, radius, angle, pos, phase;
  float dash_len, dash_gap, dash_sum;
  int ok = true;

  OOAlloc(G, RepAngle);
  CHECKOK(ok, I);
  
  PRINTFD(G, FB_RepAngle)
    "RepAngleNew: entered.\n" ENDFD;
  if(!ok || !ds->NAngleIndex) {
    OOFreeP(I);
    return (NULL);
  }

  RepInit(G, &I->R);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo * info)) RepAngleRender;
  I->R.fFree = (void (*)(struct Rep *)) RepAngleFree;
  I->R.fRecolor = NULL;

  dash_len = SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dash_length);
  dash_gap = SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if(dash_sum < R_SMALL4)
    dash_sum = 0.1F;

  I->shaderCGO = 0;
  I->N = 0;
  I->V = NULL;
  I->R.P = NULL;
  I->Obj = (CObject *) ds->Obj;
  I->ds = ds;

  n = 0;
  if(ds->NAngleIndex) {
    I->V = VLAlloc(float, ds->NAngleIndex * 10);
    CHECKOK(ok, I->V);
    for(a = 0; ok && a < ds->NAngleIndex; a = a + 5) {
      v1 = ds->AngleCoord + 3 * a;
      v2 = ds->AngleCoord + 3 * (a + 1);
      v3 = ds->AngleCoord + 3 * (a + 2);
      v4 = ds->AngleCoord + 3 * (a + 3);
      subtract3f(v1, v2, d1);
      subtract3f(v3, v2, d2);

      l1 = (float) length3f(d1);
      l2 = (float) length3f(d2);

      if(l1 > l2)
        radius = l2;
      else
        radius = l1;
      radius *= SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_angle_size);

      angle = get_angle3f(d1, d2);

      normalize23f(d1, n1);

      remove_component3f(d2, n1, d3);

      if(length3f(d3) < R_SMALL8) {
        d3[0] = 1.0F;
        d3[1] = 0.0F;
        d3[2] = 0.0F;
      } else {
        normalize23f(d3, n3);
      }

      scale3f(n1, radius, x);
      scale3f(n3, radius, y);

      if(v4[0] != 0.0F) {       /* line 1 flag */
        VLACheck(I->V, float, (n * 3) + 5);
	CHECKOK(ok, I->V);
	if (ok){
	  v = I->V + n * 3;
	  copy3f(v1, v);
	  v += 3;
	  copy3f(v2, v);
	  n += 2;
	}
      }
      
      if(ok && v4[1] != 0.0F) {       /* line 2 flag */
	VLACheck(I->V, float, (n * 3) + 5);
	CHECKOK(ok, I->V);
	if (ok){
	  v = I->V + n * 3;
	  copy3f(v3, v);
	  v += 3;
	  copy3f(v2, v);
	  n += 2;
	}
      }
      if (!ok)
	break;
      /* now we have a relevant orthogonal axes */

      length = (float) (angle * radius * 2);

      /* figure out dash/gap phasing that will lead to nicely spaced dashes and gaps */

      phase = dash_sum - (float) fmod(length / 2 + (dash_gap / 2), dash_sum);
      pos = -phase;

      if(length > R_SMALL4) {

        float mod_pos;
        float vx[3], vy[3];
        float cur_angle;
        float cons_pos1, cons_pos2;

        while(ok && pos < length) {

          mod_pos = (float) fmod(pos + phase, dash_sum);

          VLACheck(I->V, float, (n * 3) + 5);
	  CHECKOK(ok, I->V);

	  if (!ok)
	    break;
          cons_pos1 = pos;
          if(cons_pos1 < 0.0F)
            cons_pos1 = 0.0F;
          cons_pos2 = pos + dash_len;
          if(cons_pos2 > length)
            cons_pos2 = length;

          if(cons_pos1 < cons_pos2) {
            cur_angle = angle * cons_pos1 / length;

            v = I->V + n * 3;
            scale3f(x, (float) cos(cur_angle), vx);
            scale3f(y, (float) sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(v2, v, v);

            cur_angle = angle * cons_pos2 / length;

            v += 3;
            scale3f(x, (float) cos(cur_angle), vx);
            scale3f(y, (float) sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(v2, v, v);

            n += 2;
          }
          pos += dash_sum;
        }
      }
    }
    if (ok)
      VLASize(I->V, float, n * 3);
    CHECKOK(ok, I->V);
    if (ok)
      I->N = n;
  }
  if (!ok){
    RepAngleFree(I);
    I = NULL;
  }
  return ((void *) (struct Rep *) I);
}
