
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
#include"RepDihedral.h"
#include"Color.h"
#include"Scene.h"
#include"main.h"
#include"Vector.h"
#include"Setting.h"
#include"PyMOLObject.h"

typedef struct RepDihedral {
  Rep R;
  float *V;
  int N;
  CObject *Obj;
  DistSet *ds;
  float linewidth, radius;
} RepDihedral;

#include"ObjectDist.h"

void RepDihedralFree(RepDihedral * I);

void RepDihedralFree(RepDihedral * I)
{
  VLAFreeP(I->V);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepDihedralRender(RepDihedral * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  float *vc;
  int round_ends;
  int color =
    SettingGet_color(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dihedral_color);
  I->linewidth =
    SettingGet_f(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_width);
  I->radius =
    SettingGet_f(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_radius);
  round_ends =
    SettingGet_b(G, I->ds->Setting, I->ds->Obj->Obj.Setting, cSetting_dash_round_ends);

  if(ray) {

    float radius;

    if(I->radius == 0.0F) {
      radius = ray->PixelRadius * I->linewidth / 2.0F;
    } else {
      radius = I->radius;
    }

    if(color < 0)
      color = I->Obj->Color;
    vc = ColorGet(G, color);
    v = I->V;
    c = I->N;

    while(c > 0) {
      /*      printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f \n",v[3],v[4],v[5],v[6],v[7],v[8]); */
      if(round_ends) {
        ray->fSausage3fv(ray, v, v + 3, radius, vc, vc);
      } else {
        ray->fCustomCylinder3fv(ray, v, v + 3, radius, vc, vc, cCylCapFlat, cCylCapFlat);
      }
      v += 6;
      c -= 2;
    }

  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else {
      int use_dlst;

      if(info->width_scale_flag) {
        glLineWidth(I->linewidth * info->width_scale);
      } else {
        glLineWidth(I->linewidth);
      }
      if(color >= 0)
        glColor3fv(ColorGet(G, color));

      use_dlst = (int) SettingGet(G, cSetting_use_display_lists);
      if(use_dlst && I->R.displayList) {
        glCallList(I->R.displayList);
      } else {

        SceneResetNormal(G, true);

        if(use_dlst) {
          if(!I->R.displayList) {
            I->R.displayList = glGenLists(1);
            if(I->R.displayList) {
              glNewList(I->R.displayList, GL_COMPILE_AND_EXECUTE);
            }
          }
        }

        v = I->V;
        c = I->N;

        if(!info->line_lighting)
          glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        while(c > 0) {
          glVertex3fv(v);
          v += 3;
          glVertex3fv(v);
          v += 3;
          c -= 2;
        }
        glEnd();

        glEnable(GL_LIGHTING);
        if(use_dlst && I->R.displayList) {
          glEndList();
        }
      }
    }
  }
}

Rep *RepDihedralNew(DistSet * ds)
{
  PyMOLGlobals *G = ds->State.G;
  int a;
  int n;
  float *v;
  float dash_len, dash_gap, dash_sum;

  OOAlloc(G, RepDihedral);

  if(!ds->NDihedralIndex) {
    OOFreeP(I);
    return (NULL);
  }

  RepInit(G, &I->R);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo * info)) RepDihedralRender;
  I->R.fFree = (void (*)(struct Rep *)) RepDihedralFree;
  I->R.fRecolor = NULL;

  dash_len = SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dash_length);
  dash_gap = SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dash_gap);
  dash_sum = dash_len + dash_gap;
  if(dash_sum < R_SMALL4)
    dash_sum = 0.5;

  I->N = 0;
  I->V = NULL;
  I->R.P = NULL;
  I->Obj = (CObject *) ds->Obj;
  I->ds = ds;

  n = 0;
  if(ds->NDihedralIndex) {

    float *v1, *v2, *v3, *v4, *v5, *v6;

    float d12[3], d32[3], d43[3], n12[3], n32[3], n43[3];
    float p12[3], p43[3], np12[3], np43[3], v12[3], v43[3];
    float s12[3], s43[3];

    float a32[3];
    float l1, l2;
    float d3[3], n1[3], n3[3], x[3], y[3];
    float radius, length, angle, phase, pos;
    float dihedral_size =
      SettingGet_f(G, ds->Setting, ds->Obj->Obj.Setting, cSetting_dihedral_size);

    I->V = VLAlloc(float, ds->NDihedralIndex * 10);

    for(a = 0; a < ds->NDihedralIndex; a = a + 6) {
      v1 = ds->DihedralCoord + 3 * a;
      v2 = v1 + 3;
      v3 = v1 + 6;
      v4 = v1 + 9;
      v5 = v1 + 12;
      v6 = v1 + 15;

      subtract3f(v1, v2, d12);
      subtract3f(v3, v2, d32);
      subtract3f(v4, v3, d43);

      normalize23f(d12, n12);
      normalize23f(d32, n32);
      normalize23f(d43, n43);

      remove_component3f(d12, n32, p12);
      remove_component3f(d43, n32, p43);

      average3f(v2, v3, a32);

      l1 = (float) length3f(p12);
      l2 = (float) length3f(p43);

      if(l1 > l2)
        radius = l2;
      else
        radius = l1;
      radius *= dihedral_size;

      normalize23f(p12, np12);
      normalize23f(p43, np43);

      scale3f(np12, radius, v12);
      scale3f(np43, radius, v43);

      extrapolate3f(v12, n12, s12);
      add3f(s12, v2, s12);
      extrapolate3f(v43, n43, s43);
      add3f(s43, v3, s43);

      add3f(a32, v12, v12);
      add3f(a32, v43, v43);

      angle = get_angle3f(p12, p43);

      normalize23f(p12, n1);

      remove_component3f(p43, n1, d3);

      if(length3f(d3) < R_SMALL8) {
        d3[0] = 1.0F;
        d3[1] = 0.0F;
        d3[2] = 0.0F;
      } else {
        normalize23f(d3, n3);
      }

      scale3f(n1, radius, x);
      scale3f(n3, radius, y);

      VLACheck(I->V, float, (n * 3) + 5);
      v = I->V + n * 3;
      copy3f(v12, v);
      v += 3;
      copy3f(a32, v);
      n += 2;

      VLACheck(I->V, float, (n * 3) + 5);
      v = I->V + n * 3;
      copy3f(v43, v);
      v += 3;
      copy3f(a32, v);
      n += 2;

#if 0
      VLACheck(I->V, float, (n * 3) + 5);
      v = I->V + n * 3;
      copy3f(v12, v);
      v += 3;
      copy3f(s12, v);
      n += 2;

      VLACheck(I->V, float, (n * 3) + 5);
      v = I->V + n * 3;
      copy3f(v43, v);
      v += 3;
      copy3f(s43, v);
      n += 2;
#endif

      if(v5[0] != 0.0F) {       /* line 1 flag */

        VLACheck(I->V, float, (n * 3) + 5);
        v = I->V + n * 3;
        copy3f(v1, v);
        v += 3;
        copy3f(v2, v);
        n += 2;
      }

      if(v5[1] != 0.0F) {       /* line 2 flag */

        VLACheck(I->V, float, (n * 3) + 5);
        v = I->V + n * 3;
        copy3f(v3, v);
        v += 3;
        copy3f(v2, v);
        n += 2;
      }

      if(v5[2] != 0.0F) {       /* line 3 flag */

        VLACheck(I->V, float, (n * 3) + 5);
        v = I->V + n * 3;
        copy3f(v3, v);
        v += 3;
        copy3f(v4, v);
        n += 2;
      }

      /* now we have a relevant orthogonal axes */

      length = (float) (angle * radius * 2);

      /* figure out dash/gap phasing that will lead to nicely space dashes and gaps */

      phase = dash_sum - (float) fmod(length / 2 + (dash_gap / 2), dash_sum);
      pos = -phase;

      if(length > R_SMALL4) {

        float mod_pos;
        float vx[3], vy[3];
        float cur_angle;
        float cons_pos1, cons_pos2;

        while(pos < length) {

          mod_pos = (float) fmod(pos + phase, dash_sum);

          VLACheck(I->V, float, (n * 3) + 5);

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
            add3f(a32, v, v);

            cur_angle = angle * cons_pos2 / length;

            v += 3;
            scale3f(x, (float) cos(cur_angle), vx);
            scale3f(y, (float) sin(cur_angle), vy);
            add3f(vx, vy, v);
            add3f(a32, v, v);

            n += 2;
          }
          pos += dash_sum;
        }
      }
    }

    VLASize(I->V, float, n * 3);
    I->N = n;
  }
  return ((void *) (struct Rep *) I);
}
