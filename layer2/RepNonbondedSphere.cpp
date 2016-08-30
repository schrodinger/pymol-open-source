
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
#include"RepNonbondedSphere.h"
#include"Color.h"
#include"Sphere.h"
#include"Setting.h"
#include"main.h"
#include"ShaderMgr.h"
#include"Scene.h"
#include"CGO.h"

typedef struct RepNonbondedSphere {
  Rep R;
  float *V;
  float *VC;
  SphereRec *SP;
  int N, NC;
  float *VP;
  Pickable *P;
  int NP;
  int VariableAlphaFlag;
  CGO *shaderCGO;
} RepNonbondedSphere;

#include"ObjectMolecule.h"

void RepNonbondedSphereFree(RepNonbondedSphere * I);

void RepNonbondedSphereInit(void)
{
}

void RepNonbondedSphereFree(RepNonbondedSphere * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  FreeP(I->VP);
  RepPurge(&I->R);
  FreeP(I->VC);
  FreeP(I->V);
  OOFreeP(I);
}

static void RepNonbondedSphereRender(RepNonbondedSphere * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  int cc = 0;
  int a;
  SphereRec *sp;
  int i, j;
  Pickable *p;
  int ok = true;
  float alpha;

  alpha =
    SettingGet_f(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_nonbonded_transparency);
  alpha = 1.0F - alpha;
  if(fabs(alpha - 1.0) < R_SMALL4)
    alpha = 1.0F;

  if(ray) {
    int variable_alpha = I->VariableAlphaFlag;
    ray->transparentf(1.0F - alpha);
    v = I->VC;
    c = I->NC;
    while(ok && c--) {
      if(variable_alpha) {
        ray->transparentf(1.0F - v[3]);
      }
      ray->color3fv(v);
      v += 4;
      ok &= ray->sphere3fv(v, *(v + 3));
      v += 4;
    }
    ray->transparentf(0.0);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {

      i = (*pick)->src.index;

      v = I->VP;
      c = I->NP;
      p = I->R.P;

      SceneSetupGLPicking(G);
      glBegin(GL_LINES);
      while(c--) {
        i++;
        if(!(*pick)[0].src.bond) {
          /* pass 1 - low order bits */
          glColor3ub((uchar) ((i & 0xF) << 4), (uchar) ((i & 0xF0) | 0x8),
                     (uchar) ((i & 0xF00) >> 4));
          VLACheck((*pick), Picking, i);
          p++;
          (*pick)[i].src = *p;  /* copy object and atom info */
          (*pick)[i].context = I->R.context;
        } else {
          /* pass 2 - high order bits */
          j = i >> 12;
          glColor3ub((uchar) ((j & 0xF) << 4), (uchar) ((j & 0xF0) | 0x8),
                     (uchar) ((j & 0xF00) >> 4));
        }
        glVertex3fv(v);
        v += 3;
        glVertex3fv(v);
        v += 3;
        glVertex3fv(v);
        v += 3;
        glVertex3fv(v);
        v += 3;
        glVertex3fv(v);
        v += 3;
        glVertex3fv(v);
        v += 3;
      }
      glEnd();
      (*pick)[0].src.index = i;

    } else { /* rendering */
      int variable_alpha = I->VariableAlphaFlag;
      short use_shader, use_default_shader, use_sphere_shader, generate_shader_cgo = 0;
      use_shader = SettingGetGlobal_i(G, cSetting_nb_spheres_use_shader) &&
                   SettingGetGlobal_b(G, cSetting_use_shaders);
      use_sphere_shader = (SettingGetGlobal_i(G, cSetting_nb_spheres_use_shader)==1) &&
	                   SettingGetGlobal_b(G, cSetting_use_shaders);
      use_default_shader = (SettingGetGlobal_i(G, cSetting_nb_spheres_use_shader)==2) && 
	                    SettingGetGlobal_b(G, cSetting_use_shaders);

      if (I->shaderCGO){
	if (!use_shader || use_sphere_shader ^ I->shaderCGO->has_draw_sphere_buffers){
	  CGOFree(I->shaderCGO);
	  I->shaderCGO = 0;
	}
      }

      if (use_shader){
	if (!I->shaderCGO){
	  I->shaderCGO = CGONew(G);
	  CHECKOK(ok, I->shaderCGO);
	  if (ok)
	    I->shaderCGO->use_shader = true;
	  generate_shader_cgo = 1;
	} else if (ok) {
	  I->shaderCGO->enable_shaders = true;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  return;
	}
      }
      sp = I->SP;

      if (generate_shader_cgo){
	if (ok){
	  if (use_sphere_shader){
	    /* Go through the ray tracing data, its easier (and available!) */
	    int variable_alpha = I->VariableAlphaFlag;
	    ok &= CGOAlpha(I->shaderCGO, alpha);
	    v = I->VC;
	    c = I->NC;
	    while(ok && c--) {
	      if(variable_alpha) {
		ok &= CGOAlpha(I->shaderCGO, v[3]);
	      }
	      if (ok){
		ok &= CGOColorv(I->shaderCGO, v);
		v += 4;
	      }
	      if (ok){
		ok &= CGOSphere(I->shaderCGO, v, *(v + 3));
		v += 4;
	      }
	    }
	    if (ok)
	      ok &= CGOAlpha(I->shaderCGO, 1.);
	  } else {
	    while(ok && c--) {
	      if((alpha == 1.0) && (!variable_alpha)) {
		ok &= CGOAlpha(I->shaderCGO, 1.f);
		if (ok)
		  ok &= CGOColorv(I->shaderCGO, v);
	      } else {
		if(variable_alpha)
		  ok &= CGOAlpha(I->shaderCGO, v[3]);
		else
		  ok &= CGOAlpha(I->shaderCGO, alpha);
		if (ok)
		  ok &= CGOColor(I->shaderCGO, v[0], v[1], v[2]);
	      }
	      v += 4;
	      for(a = 0; ok && a < sp->NStrip; a++) {
		cc = sp->StripLen[a];
		ok &= CGOBegin(I->shaderCGO, GL_TRIANGLE_STRIP);
		while(ok && cc--) {
		  ok &= CGONormalv(I->shaderCGO, v);
		  v += 3;
		  if (ok){
		    ok &= CGOVertexv(I->shaderCGO, v);
		    v += 3;
		  }
		}
		if (ok)
		  ok &= CGOEnd(I->shaderCGO);
	      }
	    }
	  }
	}
      } else {
	while(c--) {
	  if((alpha == 1.0) && (!variable_alpha)) {
	    glColor3fv(v);
	  } else {
	    if(variable_alpha)
	      glColor4f(v[0], v[1], v[2], v[3]);
	    else
	      glColor4f(v[0], v[1], v[2], alpha);
	  }
	  v += 4;
	  for(a = 0; a < sp->NStrip; a++) {
	    cc = sp->StripLen[a];
#ifdef PURE_OPENGL_ES_2
            /* TODO */
#else
	    glBegin(GL_TRIANGLE_STRIP);
	    while(cc--) {
	      glNormal3fv(v);
	      v += 3;
	      glVertex3fv(v);
	      v += 3;
	    }
	    glEnd();
#endif
	  }
        }
      }

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
	    if (use_sphere_shader){
	      convertcgo = CGOOptimizeSpheresToVBONonIndexed(I->shaderCGO, 0);
	    } else {
	      convertcgo = CGOOptimizeToVBONotIndexed(I->shaderCGO, 0);
	    }
	    CHECKOK(ok, convertcgo);
	  }
	  if(convertcgo){
	    CGOFree(I->shaderCGO);
	    I->shaderCGO = convertcgo;
	    I->shaderCGO->use_shader = true;
	    convertcgo = NULL;
	  }
	}
	
	if (ok) {
	  I->shaderCGO->enable_shaders = true;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  return;
	}
      }
    }
  }
  if (!ok){
    CGOFree(I->shaderCGO);
    I->shaderCGO = NULL;
    I->R.fInvalidate(&I->R, I->R.cs, cRepInvPurge);
    I->R.cs->Active[cRepNonbondedSphere] = false;
  }
}

Rep *RepNonbondedSphereNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj = cs->Obj;
  int a, c, d, c1;
  float *v, *v0, *vc;
  float nb_spheres_size;
  int *q, *s;
  SphereRec *sp = G->Sphere->Sphere[0];
  int nb_spheres_quality;
  int *active = NULL;
  AtomInfoType *ai;
  int nSphere = 0;
  int a1;
  float *v1;
  float tmpColor[3];
  int variable_alpha = false;
  float transp =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_nonbonded_transparency);
  int ok = true;

  OOAlloc(G, RepNonbondedSphere);
  CHECKOK(ok, I);

  if (ok)
    active = Alloc(int, cs->NIndex);
  CHECKOK(ok, active);

  if((obj->RepVisCache & cRepNonbondedSphereBit))
    for(a = 0; a < cs->NIndex; a++) {
      ai = obj->AtomInfo + cs->IdxToAtm[a];
      active[a] = (!ai->bonded && (ai->visRep & cRepNonbondedSphereBit));
      if(active[a]) {
        active[a] = (ai->masked) ? -1 : 1;
        nSphere++;
      }
    }
  if(!nSphere) {
    OOFreeP(I);
    FreeP(active);
    return (NULL);              /* skip if no dots are visible */
  }

  nb_spheres_size =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_nb_spheres_size);

  /* get current dot sampling */
  nb_spheres_quality = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_nb_spheres_quality);
  if(nb_spheres_quality < 0)
    nb_spheres_quality = 0;
  if(nb_spheres_quality > (NUMBER_OF_SPHERE_LEVELS-1))
    nb_spheres_quality = NUMBER_OF_SPHERE_LEVELS-1;
  sp = G->Sphere->Sphere[nb_spheres_quality];

  RepInit(G, &I->R);
  I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepNonbondedSphereRender;
  I->R.fFree = (void (*)(struct Rep *)) RepNonbondedSphereFree;
  I->R.fRecolor = NULL;
  I->R.obj = (CObject *) (cs->Obj);
  I->R.cs = cs;
  I->shaderCGO = NULL;
  I->N = 0;
  I->NC = 0;
  I->V = NULL;
  I->VC = NULL;
  I->SP = NULL;
  I->NP = 0;
  I->VP = NULL;
  I->R.P = NULL;

  /* raytracing primitives */

  I->VC = (float *) mmalloc(sizeof(float) * nSphere * 8);
  CHECKOK(ok, I->VC);
  I->NC = 0;

  v = I->VC;

  for(a = 0; ok && a < cs->NIndex; a++) {
    if(active[a]) {
      float at_transp = transp;
      ai = obj->AtomInfo + cs->IdxToAtm[a];

      if(AtomSettingGetIfDefined(G, ai, cSetting_nonbonded_transparency, &at_transp))
	variable_alpha = true;
      
      I->NC++;
      c1 = ai->color;
      v0 = cs->Coord + 3 * a;
      if(ColorCheckRamped(G, c1)) {
	ColorGetRamped(G, c1, v0, tmpColor, state);
	vc = tmpColor;
      } else {
	vc = ColorGet(G, c1);
      }
      *(v++) = *(vc++);
      *(v++) = *(vc++);
      *(v++) = *(vc++);
      *(v++) = 1.0F - at_transp;
      *(v++) = *(v0++);
      *(v++) = *(v0++);
      *(v++) = *(v0++);
      *(v++) = nb_spheres_size;
    }
    ok &= !G->Interrupt;
  }

  if (ok){
    I->VariableAlphaFlag = variable_alpha;
    if(I->NC)
      I->VC = ReallocForSure(I->VC, float, (v - I->VC));
    else
      I->VC = ReallocForSure(I->VC, float, 1);
    CHECKOK(ok, I->VC);
  }
  if (ok)
    I->V = (float *) mmalloc(sizeof(float) * nSphere * (4 + sp->NVertTot * 6));
  CHECKOK(ok, I->V);

  /* rendering primitives */

  I->N = 0;
  I->SP = sp;
  v = I->V;

  for(a = 0; ok && a < cs->NIndex; a++) {
    if(active[a]) {
      float at_transp = transp;
      ai = obj->AtomInfo + cs->IdxToAtm[a];
      c1 = ai->color;
      v0 = cs->Coord + 3 * a;
      vc = ColorGet(G, c1);

      if(AtomSettingGetIfDefined(G, ai, cSetting_nonbonded_transparency, &at_transp))
	variable_alpha = true;
      
      if(ColorCheckRamped(G, c1)) {
	ColorGetRamped(G, c1, v0, tmpColor, state);
	vc = tmpColor;
      } else {
	vc = ColorGet(G, c1);
      }
      
      *(v++) = *(vc++);
      *(v++) = *(vc++);
      *(v++) = *(vc++);
      *(v++) = 1.0F - at_transp;
      
      q = sp->Sequence;
      s = sp->StripLen;
      
      for(d = 0; ok && d < sp->NStrip; d++) {
	for(c = 0; c < (*s); c++) {
	  *(v++) = sp->dot[*q][0];      /* normal */
	  *(v++) = sp->dot[*q][1];
	  *(v++) = sp->dot[*q][2];
	  *(v++) = v0[0] + nb_spheres_size * sp->dot[*q][0];     /* point */
	  *(v++) = v0[1] + nb_spheres_size * sp->dot[*q][1];
	  *(v++) = v0[2] + nb_spheres_size * sp->dot[*q][2];
	  q++;
	}
	s++;
	ok &= !G->Interrupt;
      }
      I->N++;
    }
  }

  if (ok){
    if(I->N)
      I->V = ReallocForSure(I->V, float, (v - I->V));
    else
      I->V = ReallocForSure(I->V, float, 1);
    CHECKOK(ok, I->V);
  }

  /* use pickable representation from nonbonded */
  if(ok && SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_pickable)) {
    I->VP = (float *) mmalloc(sizeof(float) * nSphere * 18);
    CHECKOK(ok, I->VP);

    if (ok)
      I->R.P = Alloc(Pickable, cs->NIndex + 1);
    CHECKOK(ok, I->R.P);

    v = I->VP;

    for(a = 0; ok && a < cs->NIndex; a++){
      if(active[a] > 0) {

        a1 = cs->IdxToAtm[a];

        if(!obj->AtomInfo[a1].masked) {
          I->NP++;

          I->R.P[I->NP].index = a1;
          I->R.P[I->NP].bond = -1;
          v1 = cs->Coord + 3 * a;

          *(v++) = v1[0] - nb_spheres_size;
          *(v++) = v1[1];
          *(v++) = v1[2];
          *(v++) = v1[0] + nb_spheres_size;
          *(v++) = v1[1];
          *(v++) = v1[2];
          *(v++) = v1[0];
          *(v++) = v1[1] - nb_spheres_size;
          *(v++) = v1[2];
          *(v++) = v1[0];
          *(v++) = v1[1] + nb_spheres_size;
          *(v++) = v1[2];
          *(v++) = v1[0];
          *(v++) = v1[1];
          *(v++) = v1[2] - nb_spheres_size;
          *(v++) = v1[0];
          *(v++) = v1[1];
          *(v++) = v1[2] + nb_spheres_size;
        }
      }
      ok &= !G->Interrupt;
    }
    if (ok)
      I->R.P = Realloc(I->R.P, Pickable, I->NP + 1);
    CHECKOK(ok, I->R.P);
    if (ok){
      I->R.context.object = (void *) obj;
      I->R.context.state = state;
      
      I->R.P[0].index = I->NP;
    }
    if (ok)
      I->VP = Realloc(I->VP, float, I->NP * 21);
    CHECKOK(ok, I->VP);
  }
  
  FreeP(active);
  if (!ok){
    RepNonbondedSphereFree(I);
    I = NULL;
  }
  return (Rep *) I;
}
