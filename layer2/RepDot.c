
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
#include"RepDot.h"
#include"Color.h"
#include"Sphere.h"
#include"Map.h"
#include"Setting.h"
#include"main.h"
#include"ObjectMolecule.h"
#include"Scene.h"
#include"ShaderMgr.h"

void RepDotFree(RepDot * I);

void RepDotInit(void)
{
}

void RepDotFree(RepDot * I)
{
  if (I->shaderCGO){
    CGOFree(I->shaderCGO);
    I->shaderCGO = 0;
  }
  FreeP(I->VC);
  FreeP(I->V);
  FreeP(I->T);
  FreeP(I->F);
  FreeP(I->VN);
  FreeP(I->A);
  FreeP(I->Atom);
  OOFreeP(I);
}

static void RepDotRender(RepDot * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  PyMOLGlobals *G = I->R.G;
  float *v = I->V;
  int c = I->N;
  int cc = 0;
  if(ray) {
    float radius;

    if(I->dotSize <= 0.0F) {
      radius = ray->PixelRadius * I->Width / 1.4142F;
    } else {
      radius = I->dotSize;
    }

    while(c--) {
      if(!cc) {                 /* load up the current vertex color */
        cc = (int) (*(v++));
        ray->fColor3fv(ray, v);
        v += 3;
      }
      v += 3;
      ray->fSphere3fv(ray, v, radius);
      v += 3;
      cc--;
    }
    /*         v=I->VC;
       c=I->NC;
       while(c--) {
       ray->fColor3fv(ray,v);
       v+=3;
       ray->fSphere3fv(ray,v,*(v+3));
       v+=4;
       } */

  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
    } else { /* else not pick, i.e., when rendering */
      short use_shader, generate_shader_cgo = 0, use_display_lists = 0;
      int normals =
        SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_normals);
      int lighting =
        SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_lighting);
      short dot_as_spheres = SettingGet_i(G, I->R.cs->Setting, I->R.obj->Setting, cSetting_dot_as_spheres);

      use_shader = (int) SettingGet(G, cSetting_dot_use_shader) & 
                   (int) SettingGet(G, cSetting_use_shaders);
      use_display_lists = (int) SettingGet(G, cSetting_use_display_lists);

      if (I->shaderCGO && ((!use_shader || CGOCheckWhetherToFree(G, I->shaderCGO)) ||
			   I->shaderCGO_as_spheres!= dot_as_spheres)){
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
	  generate_shader_cgo = 1;
	} else {
	  CShaderPrg *shaderPrg;
	  float *color;
	  color = ColorGet(G, I->R.obj->Color);

	  I->shaderCGO->enable_shaders = false;
	  if (dot_as_spheres){
	    float radius;
	    if(I->dotSize <= 0.0F) {
	      if(info->width_scale_flag)
		radius = I->Width * info->width_scale * info->vertex_scale / 1.4142F;
	      else
		radius = I->Width * info->vertex_scale;
	    } else {
	      radius = I->dotSize;
	    }
	    shaderPrg = CShaderPrg_Enable_SphereShader(G, "sphere");
	    CShaderPrg_Set1f(shaderPrg, "sphere_size_scale", fabs(radius));
	    CGORenderGL(I->shaderCGO, color, NULL, NULL, info, &I->R);
	    CShaderPrg_Disable(shaderPrg);
	  } else {
	    shaderPrg = CShaderPrg_Enable_DefaultShader(G);
	    CShaderPrg_Set1i(shaderPrg, "lighting_enabled", 0);
	    SceneResetNormalUseShaderAttribute(G, 0, true, CShaderPrg_GetAttribLocation(shaderPrg, "a_Normal"));
	    CGORenderGL(I->shaderCGO, color, NULL, NULL, info, &I->R);
	    CShaderPrg_Disable(shaderPrg);
	  }
	  return; /* should not do any other rendering after shaderCGO has
		    been rendered */
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
#endif

      if (generate_shader_cgo){
	CGO *cgo = CGONew(G);
	I->shaderCGO = CGONew(G);
	if(!normals)
	  CGOResetNormal(I->shaderCGO, true);
	if (dot_as_spheres){
	  while(c--) {
	    if(!cc) {             /* load up the current vertex color */
	      cc = (int) (*(v++));
	      CGOColorv(cgo, v);
	      v += 3;
	    }
	    if(normals)
	      CGONormalv(cgo, v);
	    v += 3;
	    CGOSphere(cgo, v, 1.f);
	    v += 3;
	    cc--;
	  }
	} else {
	  CGOLinewidthSpecial(I->shaderCGO, POINTSIZE_DYNAMIC_DOT_WIDTH);
	  CGOBegin(cgo, GL_POINTS);
	  while(c--) {
	    if(!cc) {             /* load up the current vertex color */
	      cc = (int) (*(v++));
	      CGOColorv(cgo, v);
	      v += 3;
	    }
	    v += 3;
	    CGOVertexv(cgo, v);
	    v += 3;
	    cc--;
	  }
	  CGOEnd(cgo);
	}
	CGOStop(cgo);
	{
	  if (dot_as_spheres){
	    I->shaderCGO = CGOOptimizeSpheresToVBONonIndexed(cgo, CGO_BOUNDING_BOX_SZ + CGO_DRAW_SPHERE_BUFFERS_SZ);
	  } else {
	    CGO *convertcgo = CGOCombineBeginEnd(cgo, 0), *tmpCGO;
	    tmpCGO = CGOOptimizeToVBONotIndexed(convertcgo, CGO_BOUNDING_BOX_SZ + I->N * 3 + 7);
	    CGOAppend(I->shaderCGO, tmpCGO);
	    CGOFreeWithoutVBOs(tmpCGO);
	    CGOFree(convertcgo);
	  }
	  CGOStop(I->shaderCGO);
	}
	I->shaderCGO->use_shader = true;
	I->shaderCGO_as_spheres = dot_as_spheres;
	CGOFree(cgo);

	/* now that the shaderCGO is created, we can just call RepDotRender to render it */
	RepDotRender(I, info);
      } else {
	if(!normals)
	  SceneResetNormal(G, true);
#ifdef PURE_OPENGL_ES_2
	/* TODO */
#else
	if(!lighting) {
	  if(!info->line_lighting)
	    glDisable(GL_LIGHTING);
	}
#endif

#ifdef PURE_OPENGL_ES_2
	/* TODO */
#else
	if(info->width_scale_flag)
	  glPointSize(I->Width * info->width_scale);
	else
	  glPointSize(I->Width);
#endif

#ifdef PURE_OPENGL_ES_2
	/* TODO */
#else
#ifdef _PYMOL_GL_DRAWARRAYS
	{
	  int pl = 0, nverts = c, plc = 0;
	  ALLOCATE_ARRAY(GLfloat,ptsVals,nverts*3)
	  ALLOCATE_ARRAY(GLfloat,colorVals,nverts*4)
	  ALLOCATE_ARRAY(GLfloat,normalVals,nverts*3)
	  float *cur_color;
	  while(c--) {
	    if(!cc) {             /* load up the current vertex color */
	      cc = (int) (*(v++));
	      cur_color = v;
	      v += 3;
	    }
	    colorVals[plc++] = cur_color[0]; colorVals[plc++] = cur_color[1]; colorVals[plc++] = cur_color[2]; colorVals[plc++] = 1.f;
	    if(normals){
	      normalVals[pl] = v[0]; normalVals[pl+1] = v[1]; normalVals[pl+2] = v[2];
	    }
	    v += 3;
	    ptsVals[pl++] = v[0]; ptsVals[pl++] = v[1]; ptsVals[pl++] = v[2];
	    v += 3;
	    cc--;
	  }
	  glEnableClientState(GL_VERTEX_ARRAY);
	  glEnableClientState(GL_COLOR_ARRAY);
	  if (normals){
	    glEnableClientState(GL_NORMAL_ARRAY);
	    glNormalPointer(GL_FLOAT, 0, normalVals);
	  }
	  glColorPointer(4, GL_FLOAT, 0, colorVals);
	  glVertexPointer(3, GL_FLOAT, 0, ptsVals);
	  glDrawArrays(GL_POINTS, 0, nverts);
	  glDisableClientState(GL_COLOR_ARRAY);
	  glDisableClientState(GL_VERTEX_ARRAY);
	  if (normals){
	    glDisableClientState(GL_NORMAL_ARRAY);
	  }	  
	  DEALLOCATE_ARRAY(ptsVals)
	  DEALLOCATE_ARRAY(colorVals)
	  DEALLOCATE_ARRAY(normalVals)
	}
#else
        glBegin(GL_POINTS);
        while(c--) {
          if(!cc) {             /* load up the current vertex color */
            cc = (int) (*(v++));
            glColor3fv(v);
            v += 3;
          }
          if(normals)
            glNormal3fv(v);
          v += 3;
          glVertex3fv(v);
          v += 3;
          cc--;
        }
        glEnd();
#endif
#endif

#ifdef PURE_OPENGL_ES_2
	/* TODO */
#else
        if(!lighting)
          glEnable(GL_LIGHTING);
#endif
#ifdef _PYMOL_GL_CALLLISTS
	if (use_display_lists && I->R.displayList){
	  glEndList();
	  glCallList(I->R.displayList);      
	}
#endif
      }
    }
  }
}

Rep *RepDotNew(CoordSet * cs, int state)
{
  return (RepDotDoNew(cs, cRepDotNormal, state));
}

Rep *RepDotDoNew(CoordSet * cs, int mode, int state)
{

  /* this routine does double duty - generating the dot representation,
     but also acting as our surface area computation routine.
     Modes: cRepDotNormal,cRepDotAreaType
   */
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int a, b, flag, h, k, l, i, j, c1;
  float *v, *v0, *vc, vdw, *vn;
  float *aa = NULL;
  int *tp = NULL;
  int *tf = NULL;
  float *countPtr = NULL;
  int colorCnt, lastColor;
  Vector3f v1;
  MapType *map;
  SphereRec *sp = G->Sphere->Sphere[0];
  int ds;
  float max_vdw = MAX_VDW;
  float solv_rad = 0.0;
  int inclH = true;
  int cullByFlag = false;
  int visFlag;
  int atm, *ati = NULL;
  AtomInfoType *ai1, *ai2;
  int dot_color;
  OOAlloc(G, RepDot);

  obj = cs->Obj;

  if(mode == cRepDotAreaType) { /* assume all atoms "visible" for area comp. */
    visFlag = true;
  } else {
    visFlag = false;
    if(obj->RepVisCache[cRepDot])
      for(a = 0; a < cs->NIndex; a++) {
        if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepDot]) {
          visFlag = true;
          break;
        }
      }
  }
  if(!visFlag) {
    OOFreeP(I);
    return (NULL);              /* skip if no dots are visible */
  }

  RepInit(G, &I->R);

  I->dotSize = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_dot_radius);

  I->A = NULL;
  I->T = NULL;
  I->F = NULL;
  I->V = NULL;
  I->VC = NULL;
  I->VN = NULL;
  I->Atom = NULL;
  I->R.fRecolor = NULL;
  I->shaderCGO = 0;

  I->Width = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_dot_width);
  cullByFlag = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_trim_dots);      /* are we using flags 24 & 25 */

  dot_color = SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_dot_color);   /* are we using flags 24 & 25 */
  inclH = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_dot_hydrogens);       /* are we ignoring hydrogens? */
  if(SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_dot_solvent)) {    /* are we generating a solvent surface? */
    solv_rad = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_solvent_radius); /* if so, get solvent radius */
  }

  /* get current dot sampling */
  ds = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_dot_density);

  max_vdw += solv_rad;


/* Note: significantly affects the accuracy of our area comp. */
  if(ds < 0)
    ds = 0;
  if(ds > 4)
    ds = 4;
  sp = G->Sphere->Sphere[ds];

  I->R.fRender = (void (*)(struct Rep *, RenderInfo * info)) RepDotRender;
  I->R.fFree = (void (*)(struct Rep *)) RepDotFree;
  I->R.obj = (CObject *) obj;
  I->R.cs = cs;

  I->V = (float *) mmalloc(sizeof(float) * cs->NIndex * sp->nDot * 10);
  ErrChkPtr(G, I->V);

  if(mode == cRepDotAreaType) { /* in area mode, we need to export save addl. info 
                                 * such as the normal vectors, the partial area, 
                                 * the originating atom, etc. */
    I->A = Alloc(float, cs->NIndex * sp->nDot);
    I->T = Alloc(int, cs->NIndex * sp->nDot);
    I->F = Alloc(int, cs->NIndex * sp->nDot);
    I->VN = Alloc(float, cs->NIndex * sp->nDot * 3);
    I->Atom = Alloc(int, cs->NIndex * sp->nDot);
    aa = I->A;
    tp = I->T;
    tf = I->F;
    ati = I->Atom;
    inclH = true;
    cullByFlag = true;
  }
  vn = I->VN;

  I->N = 0;
  lastColor = -1;
  colorCnt = 0;
  map = MapNew(G, max_vdw, cs->Coord, cs->NIndex, NULL);
  v = I->V;
  if(map) {
    MapSetupExpress(map);
    for(a = 0; a < cs->NIndex; a++) {
      atm = cs->IdxToAtm[a];
      ai1 = obj->AtomInfo + atm;
      if(ai1->visRep[cRepDot] || mode == cRepDotAreaType)
        if((inclH || (!ai1->hydrogen)) &&
           ((!cullByFlag) || (!(ai1->flags & cAtomFlag_exfoliate)))) {
          int at_dot_color;

          AtomInfoGetSetting_color(G, ai1, cSetting_dot_color, dot_color, &at_dot_color);

          /* If we are culling, flag 24 controls which atoms 
             will have dot surfaces generated for them.
           */
          if(at_dot_color == -1) {
            if(cs->Color)
              c1 = *(cs->Color + a);
            else
              c1 = 0;
          } else {
            c1 = at_dot_color;
          }
          v0 = cs->Coord + 3 * a;
          vdw = ai1->vdw + solv_rad;
          for(b = 0; b < sp->nDot; b++) {
            v1[0] = v0[0] + vdw * sp->dot[b][0];
            v1[1] = v0[1] + vdw * sp->dot[b][1];
            v1[2] = v0[2] + vdw * sp->dot[b][2];

            MapLocus(map, v1, &h, &k, &l);

            flag = true;

            i = *(MapEStart(map, h, k, l));
            if(i) {
              j = map->EList[i++];
              while(j >= 0) {
                ai2 = obj->AtomInfo + cs->IdxToAtm[j];
                if((inclH || (!(ai2->hydrogen))) &&
                   ((!cullByFlag) || (!(ai2->flags & cAtomFlag_ignore))))
                  /* If we are cullilng, flag 25 controls which atoms 
                     are considered "present" in the surface area 
                     calculation (i.e. able to occlude surface) */
                  if(j != a)
                    if(within3f(cs->Coord + 3 * j, v1, ai2->vdw + solv_rad)) {
                      flag = false;
                      break;
                    }
                j = map->EList[i++];
              }
            }
            if(flag) {
              switch (mode) {
              case cRepDotNormal:

                if((lastColor != c1) || ColorCheckRamped(G, c1)) {      /* new color */
                  if(countPtr)  /* after first pass */
                    *countPtr = (float) colorCnt;       /* save count */
                  colorCnt = 1;
                  countPtr = v++;
                  vc = ColorGet(G, c1); /* save new color */
                  lastColor = c1;
                  if(ColorCheckRamped(G, c1)) {
                    ColorGetRamped(G, c1, v1, v, state);
                    v += 3;
                  } else {
                    *(v++) = *(vc++);
                    *(v++) = *(vc++);
                    *(v++) = *(vc++);
                  }
                } else
                  colorCnt++;
                *(v++) = sp->dot[b][0];
                *(v++) = sp->dot[b][1];
                *(v++) = sp->dot[b][2];
                *(v++) = v1[0];
                *(v++) = v1[1];
                *(v++) = v1[2];
                I->N++;
                break;
              case cRepDotAreaType:
                *(v++) = v1[0];
                *(v++) = v1[1];
                *(v++) = v1[2];
                *(aa++) = vdw * vdw * sp->area[b];      /* area */
                *(tp++) = ai1->customType;      /* numeric type */
                *(tf++) = ai1->flags;   /* flags */
                *(vn++) = sp->dot[b][0];
                *(vn++) = sp->dot[b][1];
                *(vn++) = sp->dot[b][2];
                *(ati++) = atm;
                I->N++;
                break;
              }
            }
          }
        }
    }
    if(countPtr)
      *countPtr = (float) colorCnt;     /* save count */
    MapFree(map);
  }

  I->V = ReallocForSure(I->V, float, (v - I->V));

  if(mode == cRepDotAreaType) {
    I->A = ReallocForSure(I->A, float, (aa - I->A));
    I->T = ReallocForSure(I->T, int, (tp - I->T));
    I->F = ReallocForSure(I->F, int, (tf - I->F));
    I->VN = ReallocForSure(I->VN, float, (vn - I->VN));
    I->Atom = ReallocForSure(I->Atom, int, (ati - I->Atom));
  }
  return ((void *) (struct Rep *) I);
}
