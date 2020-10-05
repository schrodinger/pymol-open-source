
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
#include "RepSphereImmediate.h"
#include "ShaderMgr.h"
#include "Sphere.h"

#ifndef PURE_OPENGL_ES_2

static void RenderSphereMode_Immediate_Triangles(PyMOLGlobals *G, CoordSet *cs,
                                                 ObjectMolecule *obj,
                                                 int *repActive,
                                                 float sphere_scale) {
  /* triangle-based spheres */
  int ds =
      SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_sphere_quality);
  if (ds < 0) ds = 0;
  if (ds > 4) ds = 4;
  SphereRec *sp = G->Sphere->Sphere[ds];
  {
    int a;
    int nIndex = cs->NIndex;
    const AtomInfoType *atomInfo = obj->AtomInfo.data();
    const int *i2a = cs->IdxToAtm.data();
    int last_color = -1;
    const float *v = cs->Coord.data();
    int *sp_Sequence = sp->Sequence;
    int *sp_StripLen = sp->StripLen;
    int sp_NStrip = sp->NStrip;
    Vector3f *sp_dot = sp->dot;

    for (a = 0; a < nIndex; a++) {
      const AtomInfoType *ai = atomInfo + *(i2a++);
      if (GET_BIT(ai->visRep, cRepSphere)) {
        float vdw = ai->vdw * sphere_scale;
        int c = ai->color;
        float v0 = v[0];
        float v1 = v[1];
        float v2 = v[2];
        (*repActive) = true;

        if (c != last_color) {
          last_color = c;
          glColor3fv(ColorGet(G, c));
        }

        {
          int *s = sp_StripLen;
          int *q = sp_Sequence;
          int b;
          for (b = 0; b < sp_NStrip; b++) {
            int nc = *(s++);
            glBegin(GL_TRIANGLE_STRIP);
            for (c = 0; c < nc; c++) {
              float *sp_dot_q = &sp_dot[*(q++)][0];
              glNormal3fv(sp_dot_q); /* normal */
              glVertex3f(v0 + vdw * sp_dot_q[0], v1 + vdw * sp_dot_q[1],
                         v2 + vdw * sp_dot_q[2]);
            }
            glEnd();
          }
        }
      }
      v += 3;
    }
  }
}

static void RenderSphereMode_Immediate_1_2_3(PyMOLGlobals *G, RenderInfo *info,
                                             CoordSet *cs, ObjectMolecule *obj,
                                             int *repActive, float pixel_scale,
                                             int sphere_mode) {
  /* sphere_mode is 1, 2, or 3 */
  float max_radius = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                  cSetting_sphere_point_max_size) *
                     3 * pixel_scale;
  int clamp_size_flag = (max_radius >= 0.0F);

  int a;
  int nIndex = cs->NIndex;
  const AtomInfoType *atomInfo = obj->AtomInfo.data();
  const int* i2a = cs->IdxToAtm.data();
  int last_color = -1;
  const float *v = cs->Coord.data();
  float last_radius = -1.0F;

  if (!info->line_lighting) glDisable(GL_LIGHTING);

  glBegin(GL_POINTS);
  for (a = 0; a < nIndex; a++) {
    const AtomInfoType *ai = atomInfo + *(i2a++);
    if (GET_BIT(ai->visRep, cRepSphere)) {
      int c = ai->color;
      (*repActive) = true;
      if (c != last_color) {
        last_color = c;
        glColor3fv(ColorGet(G, c));
      }
      switch (sphere_mode) {
        case 1:
        case 6:
          glVertex3fv(v);
          break;
        case 2:
        case 3:
        case 7:
        case 8: {
          float cur_radius = ai->vdw * pixel_scale;
          if (last_radius != cur_radius) {
            glEnd();
            if (clamp_size_flag)
              if (cur_radius > max_radius) cur_radius = max_radius;
            glPointSize(cur_radius);
            glBegin(GL_POINTS);
            last_radius = cur_radius;
          }
          glVertex3fv(v);
        } break;
      }
    }
    v += 3;
  }
  glEnd();

  glEnable(GL_LIGHTING);

  if (sphere_mode == 3) {
    glDisable(GL_POINT_SMOOTH);
    glAlphaFunc(GL_GREATER, 0.05F);
  } else {
    glEnable(GL_ALPHA_TEST);
  }
}

static
void RenderImmediate_DoPreGL(PyMOLGlobals *G, int sphere_mode,
                             float *pixel_scale, CoordSet *cs,
                             ObjectMolecule *obj, float sphere_scale) {
  switch (sphere_mode) {
    case 2:
    case 7:
      glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_ALPHA_TEST);
      (*pixel_scale) *= 1.4F;
      glPointSize(1.0F);
      break;
    case 3:
    case 8:
      glEnable(GL_POINT_SMOOTH);
      glAlphaFunc(GL_GREATER, 0.5F);
      glEnable(GL_ALPHA_TEST);
      glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
      glPointSize(1.0F);
      (*pixel_scale) *= 2.0F;
      break;
    case 4:
      glEnable(GL_POINT_SMOOTH);
      glEnable(GL_ALPHA_TEST);
      glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
      glPointSize(1.0F);
      (*pixel_scale) *= 2.0F;
      break;
    default:
      glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_ALPHA_TEST);
      glPointSize(SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                               cSetting_sphere_point_size));
      break;
  }
}
#endif

void RepSphereRenderImmediate(CoordSet *cs, RenderInfo *info) {
#ifndef PURE_OPENGL_ES_2
  PyMOLGlobals *G = cs->G;
  if (info->ray || info->pick || (!(G->HaveGUI && G->ValidContext)))
    return;
  else {
    int repActive = false;
    ObjectMolecule *obj = cs->Obj;
    int sphere_mode =
        SettingGet_i(G, cs->Setting.get(), obj->Setting.get(), cSetting_sphere_mode);
    float sphere_scale =
        SettingGet_f(G, cs->Setting.get(), obj->Setting.get(), cSetting_sphere_scale);

    if (sphere_mode > 0) { /* point-based modees */
      float pixel_scale = 1.0F / info->vertex_scale;
      RenderImmediate_DoPreGL(G, sphere_mode, &pixel_scale, cs, obj,
                              sphere_scale);
      switch (sphere_mode) {
        case 5:
          // removed (was ARB shader)
          break;
        case 4:
          // sphere_mode 4 taken out: many points per sphere to make them look
          // good, one specular light
          break;
        default:
          RenderSphereMode_Immediate_1_2_3(G, info, cs, obj, &repActive,
                                           pixel_scale, sphere_mode);
      }
    } else {
      RenderSphereMode_Immediate_Triangles(G, cs, obj, &repActive,
                                           sphere_scale);
    }

    if (!repActive) /* didn't draw a single sphere, so we can skip this
                       representation next time around */
      cs->Active[cRepSphere] = false;
  }
#endif
}
