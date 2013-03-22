
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
  register PyMOLGlobals *G = I->R.G;
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
    ray->fTransparentf(ray, 1.0F - alpha);
    v = I->VC;
    c = I->NC;
    while(ok && c--) {
      if(variable_alpha) {
        ray->fTransparentf(ray, 1.0F - v[3]);
      }
      ray->fColor3fv(ray, v);
      v += 4;
      ok &= ray->fSphere3fv(ray, v, *(v + 3));
      v += 4;
    }
    ray->fTransparentf(ray, 0.0);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {

      i = (*pick)->src.index;

      v = I->VP;
      c = I->NP;
      p = I->R.P;

      SceneSetupGLPicking(G);
#ifdef _PYMOL_GL_DRAWARRAYS
      {
	int nverts = c * 6, pl, plc = 0;
	GLubyte *tmp_ptr;
	ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
	ALLOCATE_ARRAY(GLubyte,colorVals,nverts*4)
	pl = 0;
	while(c--) {
	  i++;
	  if(!(*pick)[0].src.bond) {
	    /* pass 1 - low order bits */
	    colorVals[plc++] = (uchar) ((i & 0xF) << 4);
	    colorVals[plc++] = (uchar) ((i & 0xF0) | 0x8);
	    colorVals[plc++] = (uchar) ((i & 0xF00) >> 4);
	    colorVals[plc++] = (uchar) 255;
	    VLACheck((*pick), Picking, i);
	    p++;
	    (*pick)[i].src = *p;  /* copy object and atom info */
	    (*pick)[i].context = I->R.context;
	  } else {
	    /* pass 2 - high order bits */
	    j = i >> 12;
	    colorVals[plc++] = (uchar) ((j & 0xF) << 4);
	    colorVals[plc++] = (uchar) ((j & 0xF0) | 0x8);
	    colorVals[plc++] = (uchar) ((j & 0xF00) >> 4);
	    colorVals[plc++] = (uchar) 255;
	  }
	  memcpy(&vertVals[pl], v, 18*sizeof(GLfloat));
	  v += 18;
	  pl += 18;
	  tmp_ptr = &colorVals[plc-4];
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = tmp_ptr[3];
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = tmp_ptr[3];
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = tmp_ptr[3];
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = tmp_ptr[3];
	  colorVals[plc++] = tmp_ptr[0]; colorVals[plc++] = tmp_ptr[1]; colorVals[plc++] = tmp_ptr[2]; colorVals[plc++] = tmp_ptr[3];
	}
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertVals);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, colorVals);
	glDrawArrays(GL_LINES, 0, nverts);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	DEALLOCATE_ARRAY(vertVals)
	DEALLOCATE_ARRAY(colorVals)
      }
#else
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
#endif
      (*pick)[0].src.index = i;

    } else { /* rendering */
      int variable_alpha = I->VariableAlphaFlag;
      short use_shader, use_default_shader, use_sphere_shader, generate_shader_cgo = 0, use_display_lists = 0;
      use_shader = (int) SettingGet(G, cSetting_nb_spheres_use_shader) &&
	           (int) SettingGet(G, cSetting_use_shaders);
      use_sphere_shader = (int) (SettingGet(G, cSetting_nb_spheres_use_shader)==1) &&
	                  (int) SettingGet(G, cSetting_use_shaders);
      use_default_shader = (int) (SettingGet(G, cSetting_nb_spheres_use_shader)==2) && 
	                   (int) SettingGet(G, cSetting_use_shaders);
      use_display_lists = (int) SettingGet(G, cSetting_use_display_lists);

      if (I->shaderCGO){
	if (!use_shader || use_sphere_shader ^ I->shaderCGO->has_draw_sphere_buffers){
	  CGOFree(I->shaderCGO);
	  I->shaderCGO = 0;
	}
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
	} else if (ok) {
	  I->shaderCGO->enable_shaders = true;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
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
#ifdef _PYMOL_GL_DRAWARRAYS
	    {
	      int nverts = cc, pl;
	      ALLOCATE_ARRAY(GLfloat,vertVals,nverts*3)
	      ALLOCATE_ARRAY(GLfloat,normVals,nverts*3)
	      pl = 0;
	      while(cc--) {
		normVals[pl] = v[0]; normVals[pl+1] = v[1]; normVals[pl+2] = v[2];
		v += 3;
		vertVals[pl++] = v[0]; vertVals[pl++] = v[1]; vertVals[pl++] = v[2];
		v += 3;
	      }
	      glEnableClientState(GL_VERTEX_ARRAY);
	      glEnableClientState(GL_NORMAL_ARRAY);
	      glVertexPointer(3, GL_FLOAT, 0, vertVals);
	      glNormalPointer(GL_FLOAT, 0, normVals);
	      glDrawArrays(GL_TRIANGLE_STRIP, 0, nverts);
	      glDisableClientState(GL_NORMAL_ARRAY);
	      glDisableClientState(GL_VERTEX_ARRAY);
	      DEALLOCATE_ARRAY(vertVals)
	      DEALLOCATE_ARRAY(normVals)
	    }
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
#else
	  (void)convertcgo;
#endif
	}
	
	if (ok) {
	  I->shaderCGO->enable_shaders = true;
	  CGORenderGL(I->shaderCGO, NULL, NULL, NULL, info, &I->R);
	  return;
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

  if(obj->RepVisCache[cRepNonbondedSphere])
    for(a = 0; a < cs->NIndex; a++) {
      ai = obj->AtomInfo + cs->IdxToAtm[a];
      active[a] = (!ai->bonded) && (ai->visRep[cRepNonbondedSphere]);
      if(active[a]) {
        if(ai->masked)
          active[a] = -1;
        else
          active[a] = 1;
      }
      if(active[a]) {
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
      float at_transp;
      ai = obj->AtomInfo + cs->IdxToAtm[a];
      if(AtomInfoGetSetting_f(G, ai, cSetting_nonbonded_transparency, transp, &at_transp))
	variable_alpha = true;
      
      I->NC++;
      c1 = *(cs->Color + a);
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
      float at_transp;
      ai = obj->AtomInfo + cs->IdxToAtm[a];
      c1 = *(cs->Color + a);
      v0 = cs->Coord + 3 * a;
      vc = ColorGet(G, c1);
      if(AtomInfoGetSetting_f(G, ai, cSetting_nonbonded_transparency, transp, &at_transp))
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
  return ((void *) (struct Rep *) I);
}
