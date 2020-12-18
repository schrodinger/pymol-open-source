
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
#include"RepEllipsoid.h"
#include"Color.h"
#include"Setting.h"
#include"Feedback.h"
#include"Matrix.h"
#include"CGO.h"

struct RepEllipsoid : Rep {
  using Rep::Rep;

  ~RepEllipsoid() override;

  cRep_t type() const override { return cRepEllipsoid; }
  void render(RenderInfo* info) override;

  CGO* ray = nullptr;
  CGO* std = nullptr;
  CGO* shaderCGO = nullptr;
};

#include"ObjectMolecule.h"

RepEllipsoid::~RepEllipsoid()
{
  auto I = this;
  CGOFree(I->ray);
  CGOFree(I->std);
  CGOFree(I->shaderCGO);
}

void RepEllipsoid::render(RenderInfo* info)
{
  auto I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  int ok = true;

  PyMOLGlobals *G = I->G;
  if(ray) {
    int try_std = false;
    PRINTFD(G, FB_RepEllipsoid)
      " RepEllipsoidRender: rendering ray...\n" ENDFD;

    if(I->ray){
      int rayok = CGORenderRay(I->ray, ray, info, NULL, NULL, I->cs->Setting.get(), I->obj->Setting.get());
      if (!rayok){
	CGOFree(I->ray);
	try_std = true;
      }
    } else {
      try_std = true;
    }
    if(try_std && I->std){
      ok &= CGORenderRay(I->std, ray, info, NULL, NULL, I->cs->Setting.get(), I->obj->Setting.get());
      if (!ok){
	CGOFree(I->std);
      }
    }
    CHECKOK(ok, I->std);
  } else if(G->HaveGUI && G->ValidContext) {

    if(pick) {
      if(I->shaderCGO) {
        CGORenderGLPicking(I->shaderCGO, info, &I->context,
                           I->cs->Setting.get(), I->obj->Setting.get());
      } else if(I->std) {
        CGORenderGLPicking(I->std, info, &I->context,
                           I->cs->Setting.get(), I->obj->Setting.get());
      }
    } else {
      int use_shaders;
      use_shaders = SettingGetGlobal_b(G, cSetting_use_shaders);
      
        PRINTFD(G, FB_RepEllipsoid)
          " RepEllipsoidRender: rendering GL...\n" ENDFD;

	if (use_shaders){
	  if (!I->shaderCGO){
            I->shaderCGO = CGOOptimizeToVBONotIndexed(I->std, 0);
            assert(I->shaderCGO->use_shader);
	  }
	} else {
	  CGOFree(I->shaderCGO);	  
	}
	if (I->shaderCGO){
          CGORenderGL(I->shaderCGO, NULL, I->cs->Setting.get(), I->obj->Setting.get(), info, I);
	} else if(I->std){
          CGORenderGL(I->std, NULL, I->cs->Setting.get(), I->obj->Setting.get(), info, I);
	}
    }
  }
}

const double problevel[50] = { 0.4299, 0.5479, 0.6334, 0.7035, 0.7644,
  0.8192, 0.8694, 0.9162, 0.9605, 1.0026,
  1.0430, 1.0821, 1.1200, 1.1570, 1.1932,
  1.2288, 1.2638, 1.2985, 1.3330, 1.3672,
  1.4013, 1.4354, 1.4695, 1.5037, 1.5382,
  1.5729, 1.6080, 1.6436, 1.6797, 1.7164,
  1.7540, 1.7924, 1.8318, 1.8724, 1.9144,
  1.9580, 2.0034, 2.0510, 2.1012, 2.1544,
  2.2114, 2.2730, 2.3404, 2.4153, 2.5003,
  2.5997, 2.7216, 2.8829, 3.1365, 6.0000
};

/**
 * Return true if backbone atom that should be hidden with side_chain_helper
 */
static bool is_sidechainhelper_hidden(PyMOLGlobals * G, const AtomInfoType * ai) {
  if (!(ai->flags & cAtomFlag_polymer))
    return false;

  switch (ai->protons) {
    case cAN_C:
      return ai->name == G->lex_const.C;
    case cAN_N:
      return ai->name == G->lex_const.N && ai->resn != G->lex_const.PRO;
    case cAN_O:
      return ai->name == G->lex_const.O;
  }

  return false;
}

Rep *RepEllipsoidNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->G;
  ObjectMolecule *obj;
  int ok = true;

  // skip if no dots are visible
  if(!cs->hasRep(cRepEllipsoidBit))
    return NULL;

  auto I = new RepEllipsoid(cs, state);
  CHECKOK(ok, I);
  if (!ok)
    return NULL;

  obj = cs->Obj;

  {
    int ellipsoid_color = SettingGet_color(G, cs->Setting.get(), obj->Setting.get(),
                                           cSetting_ellipsoid_color);

    int cartoon_side_chain_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                                 cSetting_cartoon_side_chain_helper);

    int ribbon_side_chain_helper = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                                cSetting_ribbon_side_chain_helper);

    float ellipsoid_scale = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                         cSetting_ellipsoid_scale);

    float transp = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                                cSetting_ellipsoid_transparency);

    int pickable = SettingGet_b(G, cs->Setting.get(), obj->Setting.get(),
                                cSetting_pickable);

    float prob = SettingGet_f(G, cs->Setting.get(), obj->Setting.get(),
                              cSetting_ellipsoid_probability);
    double matrix_factor = 0.0F;
    float pradius = 0.0F;
    {

      int iprob = (prob + 0.01F) * 50.0F - 1;
      if(iprob < 0)
        iprob = 0;
      if(iprob > 49)
        iprob = 49;
      pradius = problevel[iprob];
      matrix_factor = -(1 / (pradius * pradius));
    }

    I->ray = CGONew(G);         /* describe the ellipsoids analytically */
    CHECKOK(ok, I->ray);
    if(I->ray) {
      int a, a1;
      AtomInfoType *ai;
      float last_alpha = 1.0F;

      double *csmatrix = SettingGet_i(G, cs->Setting.get(), obj->Setting.get(),
            cSetting_matrix_mode) > 0 ? NULL : cs->Matrix.data();

      for(a = 0; a < cs->NIndex; a++) {
        a1 = cs->IdxToAtm[a];
        ai = obj->AtomInfo + a1;
        if (!ai->anisou || !(ai->visRep & cRepEllipsoidBit))
          continue;

        if (is_sidechainhelper_hidden(G, ai)) {
          if ((ai->visRep & cRepCartoonBit) && AtomSettingGetWD(G, ai,
                cSetting_cartoon_side_chain_helper, /* d= */ cartoon_side_chain_helper))
            continue;

          if ((ai->visRep & cRepRibbonBit) && AtomSettingGetWD(G, ai,
                cSetting_ribbon_side_chain_helper, /* d= */ ribbon_side_chain_helper))
            continue;
        }

        {
          {
            int n_rot;
            double matrix[16];
            double e_val[4];
            double e_vec[16];

            matrix[0] = ai->anisou[0];  // U11
            matrix[1] = ai->anisou[3];  // U12
            matrix[2] = ai->anisou[4];  // U13
            matrix[3] = 0.0;
            matrix[4] = ai->anisou[3];  // U12
            matrix[5] = ai->anisou[1];  // U22
            matrix[6] = ai->anisou[5];  // U23
            matrix[7] = 0.0;
            matrix[8] = ai->anisou[4];  // U13
            matrix[9] = ai->anisou[5];  // U23
            matrix[10] = ai->anisou[2]; // U33
            matrix[11] = 0.0;
            matrix[12] = 0.0;
            matrix[13] = 0.0;
            matrix[14] = 0.0;
            matrix[15] = matrix_factor;

            if(xx_matrix_jacobi_solve(e_vec, e_val, &n_rot, matrix, 4)) {

              const float* v = cs->coordPtr(a);

              float mag[3];
              float scale[3];

              float mx;
              float r_el, n0[3], n1[3], n2[3];

              float at_ellipsoid_scale  = AtomSettingGetWD(G, ai, cSetting_ellipsoid_scale, ellipsoid_scale);
              float at_transp           = AtomSettingGetWD(G, ai, cSetting_ellipsoid_transparency, transp);
              int c1                    = AtomSettingGetWD(G, ai, cSetting_ellipsoid_color, ellipsoid_color);

              if(c1 == -1)
                c1 = ai->color;

              if(csmatrix)
                left_multiply44d44d(csmatrix, e_vec);

              n0[0] = e_vec[0];
              n0[1] = e_vec[4];
              n0[2] = e_vec[8];
              n1[0] = e_vec[1];
              n1[1] = e_vec[5];
              n1[2] = e_vec[9];
              n2[0] = e_vec[2];
              n2[1] = e_vec[6];
              n2[2] = e_vec[10];

              normalize3f(n0);
              normalize3f(n1);
              normalize3f(n2);
              mag[0] = sqrt1f(e_val[0]);
              mag[1] = sqrt1f(e_val[1]);
              mag[2] = sqrt1f(e_val[2]);

              mx = mag[0];
              if(mx < mag[1])
                mx = mag[1];
              if(mx < mag[2])
                mx = mag[2];

              scale[0] = mag[0] / mx;
              scale[1] = mag[1] / mx;
              scale[2] = mag[2] / mx;

              scale3f(n0, scale[0], n0);
              scale3f(n1, scale[1], n1);
              scale3f(n2, scale[2], n2);

              r_el = mx * pradius * at_ellipsoid_scale;

              {
                float vc[3];
                if(ColorCheckRamped(G, c1)) {
                  ColorGetRamped(G, c1, v, vc, state);
                  ok &= CGOColorv(I->ray, vc);
                } else {
                  ok &= CGOColorv(I->ray, ColorGet(G, c1));
                }
              }

              if (ok) {
                float alpha = 1.0F - at_transp;
                if(alpha != last_alpha) {
                  ok &= CGOAlpha(I->ray, alpha);
                  last_alpha = alpha;

                  if (at_transp > 0) {
                    I->setHasTransparency();
                  }
                }
              }
              if(ok && pickable && (!ai->masked))
                ok &= CGOPickColor(I->ray, a1, cPickableAtom);

              if (ok)
		ok &= CGOEllipsoid(I->ray, v, r_el, n0, n1, n2);
            }
          }
        }
      }
      if (ok)
	ok &= CGOStop(I->ray);
      I->std = CGOSimplify(I->ray, 0);  /* convert analytical to discrete */
      CHECKOK(ok, I->std);
    }
  }
  if (!ok){
    delete I;
    I = NULL;
  }
  return (Rep *) I;
}
