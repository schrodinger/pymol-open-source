

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
-* Cameron Mura
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
#include"RepCartoon.h"
#include"Color.h"
#include"Setting.h"
#include"Word.h"
#include"Scene.h"
#include"main.h"
#include"Feedback.h"
#include"CGO.h"
#include"Extrude.h"
#include"ShaderMgr.h"

typedef struct RepCartoon {
  Rep R;                        /* must be first! */
  CGO *ray, *std;

} RepCartoon;

#include"ObjectMolecule.h"

#define ESCAPE_MAX 500

void RepCartoonFree(RepCartoon * I);

void RepCartoonFree(RepCartoon * I)
{
  if(I->ray)
    CGOFree(I->ray);
  if(I->std)
    CGOFree(I->std);
  RepPurge(&I->R);
  OOFreeP(I);
}

static void RepCartoonRender(RepCartoon * I, RenderInfo * info)
{
  float *fog_color, fog_enabled;
  CRay *ray = info->ray;
  Picking **pick = info->pick;
  register PyMOLGlobals *G = I->R.G;
  CShaderPrg * p = CShaderMgr_GetShaderPrg(G->ShaderMgr, "default");

  if(ray) {
    PRINTFD(G, FB_RepCartoon)
      " RepCartoonRender: rendering raytracable...\n" ENDFD;

    if(I->ray)
      CGORenderRay(I->ray, ray, NULL, I->R.cs->Setting, I->R.obj->Setting);
    else if(I->std)
      CGORenderRay(I->std, ray, NULL, I->R.cs->Setting, I->R.obj->Setting);
  } else if(G->HaveGUI && G->ValidContext) {
    if(pick) {
      if(I->std) {
        CGORenderGLPicking(I->std, pick, &I->R.context,
                           I->R.cs->Setting, I->R.obj->Setting);
      }
    } else {
      int use_dlst, use_shader;
      use_dlst = (int) SettingGet(G, cSetting_use_display_lists);

      use_shader = (int) SettingGet(G, cSetting_cartoon_use_shader) & 
	               (int) SettingGet(G, cSetting_use_shaders) ;
      if (use_shader) {
	CShaderPrg_Enable(p);
	/*	ShaderEnable(G); */
        fog_color = SettingGetfv(G, cSetting_bg_rgb);
        fog_enabled = SettingGet(G, cSetting_depth_cue) ? 1.0 : 0.0;
        CShaderPrg_Set1f(p, "fog_r", fog_color[0]);
        CShaderPrg_Set1f(p, "fog_g", fog_color[1]);
        CShaderPrg_Set1f(p, "fog_b", fog_color[2]);
        CShaderPrg_Set1f(p, "fog_enabled", fog_enabled);
      }

      if(use_dlst && I->R.displayList) {
        glCallList(I->R.displayList);
      } else {

        if(use_dlst) {
          if(!I->R.displayList) {
            I->R.displayList = glGenLists(1);
            if(I->R.displayList) {
              glNewList(I->R.displayList, GL_COMPILE_AND_EXECUTE);
            }
          }
        }

        PRINTFD(G, FB_RepCartoon)
          " RepCartoonRender: rendering GL...\n" ENDFD;

        if(I->std)
          CGORenderGL(I->std, NULL, I->R.cs->Setting, I->R.obj->Setting, info);

        if(use_dlst && I->R.displayList) {
          glEndList();
        }
      }
      if (use_shader) {
        CShaderPrg_Disable(p);
        /*ShaderDisable(G);*/
      }
    }
  }
}

static float smooth(float x, float power)
{

  if(x <= 0.5) {
    if(x <= 0.0)
      x = 0.0;
    return ((float) (0.5 * pow(2.0 * x, power)));
  } else {
    if(x >= 1.0)
      x = 1.0;
    return ((float) (1.0 - (0.5 * pow(2 * (1.0 - x), power))));
  }
}

#define NUCLEIC_NORMAL0 "C2"
#define NUCLEIC_NORMAL1 "C3*"
#define NUCLEIC_NORMAL2 "C3'"

#define MAX_RING_ATOM 10


/* atix must contain n_atom + 1 elements, with the first atom repeated at the end */

static void do_ring(PyMOLGlobals * G, int n_atom, int *atix, ObjectMolecule * obj,
                    CoordSet * cs, float ring_width, CGO * cgo, int ring_color,
                    int ring_mode, float ladder_radius, int ladder_color, int ladder_mode,
                    int finder, int sc_helper, int *nuc_flag, int na_mode,
                    float ring_alpha, float alpha, int *marked, float *moved,
                    float ring_radius)
{
  float *v_i[MAX_RING_ATOM];
  float *col[MAX_RING_ATOM];
  float n_up[MAX_RING_ATOM][3];
  float n_dn[MAX_RING_ATOM][3];
  AtomInfoType *ai_i[MAX_RING_ATOM];
  int have_all = true;
  int all_marked = true;
  AtomInfoType *ai;
  int have_C4 = -1;
  int have_C4_prime = -1;
  int have_C_number = -1;
  int nf = false;

  /* first, make sure all atoms have known coordinates */
  {
    int a, i;
    for(i = 0; i <= n_atom; i++) {
      int a1 = atix[i];
      int have_atom = false;
      if(nuc_flag[a1])
        nf = true;
      if(obj->DiscreteFlag) {
        if(cs == obj->DiscreteCSet[a1])
          a = obj->DiscreteAtmToIdx[a1];
        else
          a = -1;
      } else
        a = cs->AtmToIdx[a1];
      if(a >= 0) {
        ai = obj->AtomInfo + a1;
        if(ai->visRep[cRepCartoon]) {
          ai_i[i] = ai;

          {
            int atom_ring_mode = ring_mode;
            int atom_ring_color = ladder_color;
            float atom_ring_radius = ring_radius;
            float atom_ring_width = ring_width;
            int atom_ladder_mode = ladder_mode;
            int atom_ladder_color = ladder_color;
            float atom_ladder_radius = ladder_radius;

            if(AtomInfoGetSetting_i
               (G, ai, cSetting_cartoon_ring_mode, ring_mode, &atom_ring_mode)) {
              ring_mode = atom_ring_mode;
            }
            if(AtomInfoGetSetting_color
               (G, ai, cSetting_cartoon_ring_color, ring_color, &atom_ring_color)) {
              ring_color = atom_ring_color;
            }
            if(AtomInfoGetSetting_f
               (G, ai, cSetting_cartoon_ring_radius, ring_radius, &atom_ring_radius)) {
              ring_radius = atom_ring_radius;
            }
            if(AtomInfoGetSetting_f
               (G, ai, cSetting_cartoon_ring_width, ring_width, &atom_ring_width)) {
              ring_width = atom_ring_width;
            }
            if(AtomInfoGetSetting_i
               (G, ai, cSetting_cartoon_ladder_mode, ladder_mode, &atom_ladder_mode)) {
              ladder_mode = atom_ladder_mode;
            }
            if(AtomInfoGetSetting_color
               (G, ai, cSetting_cartoon_ladder_color, ladder_color, &atom_ladder_color)) {
              ladder_color = atom_ladder_color;
            }
            if(AtomInfoGetSetting_f
               (G, ai, cSetting_cartoon_ladder_radius, ladder_radius,
                &atom_ladder_radius)) {
              ladder_radius = atom_ladder_radius;
            }
          }

          col[i] = ColorGet(G, ai->color);
          v_i[i] = cs->Coord + 3 * a;
          have_atom = true;
          if(WordMatchExact(G, "C4", ai->name, 1))
            have_C4 = a1;
          if(WordMatchExact(G, "C4'", ai->name, 1) ||
             WordMatchExact(G, "C4*", ai->name, 1))
            have_C4_prime = a1;
          if(((ai->name[0] == 'C') || (ai->name[0] == 'c')) &&
             ((ai->name[1] >= '0') && (ai->name[1] <= '9')))
            have_C_number = a1;
        }
        if(!marked[a1])
          all_marked = false;
      }
      if(!have_atom) {
        have_all = false;
        break;
      }
    }
  }

  if(!ring_mode)
    finder = 0;

  ring_width *= 0.5F;

  if(n_atom && have_all && (!all_marked)) {
    if(ladder_mode) {
      int i;
      int a1;
      AtomInfoType *ai2;
      register AtomInfoType *atomInfo = obj->AtomInfo;
      register int mem0, mem1, mem2, mem3, mem4, mem5, mem6, mem7;
      register int *neighbor = obj->Neighbor;
      int nbr[7];
      int sugar_at = -1, base_at = -1;
      int phos3_at = -1, phos5_at = -1;
      int o3_at = -1, o5_at = -1, c1_at = -1;
      int purine_flag = false;
      int c1 = -1, c1_linked = -1;
      int c2 = -1, c2_linked = -1;
      int c3 = -1, c3_linked = -1;
      int c4 = -1, c4_linked = -1;
      int c5 = -1, c5_linked = -1;

      if(n_atom == 5) {         /* going from the sugar to the base */
        for(i = 0; i < n_atom; i++) {
          a1 = atix[i];
          if(!nf)
            nf = nuc_flag[a1];
          ai = atomInfo + a1;
          if((ai->protons == cAN_C) && (!marked[a1]) &&
             (WordMatchExact(G, "C3*", ai->name, 1) ||
              WordMatchExact(G, "C3'", ai->name, 1))) {
            sugar_at = a1;
            mem0 = a1;
            nbr[0] = neighbor[mem0] + 1;
            while((mem1 = neighbor[nbr[0]]) >= 0) {
              if((atomInfo[mem1].protons == cAN_O) && (!marked[mem1])) {
                ai = atomInfo + mem1;
                if(WordMatchExact(G, "O3*", ai->name, 1) ||
                   WordMatchExact(G, "O3'", ai->name, 1))
                  o3_at = mem1;
                nbr[1] = neighbor[mem1] + 1;
                while((mem2 = neighbor[nbr[1]]) >= 0) {
                  if((mem2 != mem0) && (!marked[mem2]) &&
                     (atomInfo[mem2].protons == cAN_P)) {
                    phos3_at = mem2;
                  }
                  nbr[1] += 2;
                }
              }

              if((atomInfo[mem1].protons == cAN_C) && (!marked[mem1])) {
                ai2 = atomInfo + mem1;
                if(WordMatchExact(G, NUCLEIC_NORMAL1, ai2->name, 1) ||
                   WordMatchExact(G, NUCLEIC_NORMAL2, ai2->name, 1))
                  sugar_at = mem1;

                nbr[1] = neighbor[mem1] + 1;
                while((mem2 = neighbor[nbr[1]]) >= 0) {
                  if((mem2 != mem0) && (!marked[mem2]) &&
                     (atomInfo[mem2].protons == cAN_C)) {
                    ai = atomInfo + mem2;
                    if(WordMatchExact(G, "C1*", ai->name, 1) ||
                       WordMatchExact(G, "C1'", ai->name, 1))
                      c1_at = mem2;
                    nbr[2] = neighbor[mem2] + 1;
                    while((mem3 = neighbor[nbr[2]]) >= 0) {
                      if((mem3 != mem1) && (mem3 != mem0)) {
                        if((atomInfo[mem3].protons == cAN_O) && (!marked[mem3])) {
                          ai = atomInfo + mem3;
                          if(WordMatchExact(G, "O5*", ai->name, 1) ||
                             WordMatchExact(G, "O5'", ai->name, 1))
                            o5_at = mem3;

                          nbr[3] = neighbor[mem3] + 1;
                          while((mem4 = neighbor[nbr[3]]) >= 0) {
                            if((mem4 != mem2) && (!marked[mem4]) &&
                               (atomInfo[mem4].protons == cAN_P)) {
                              phos5_at = mem4;
                            }
                            nbr[3] += 2;
                          }
                        }
                        if(atomInfo[mem3].protons == cAN_N) {

                          if(ring_mode) {
                            ai2 = atomInfo + mem3;
                            if((!marked[mem3])
                               && (WordMatchExact(G, "N1", ai2->name, 1)
                                   || WordMatchExact(G, "N9", ai2->name, 1))) {
                              base_at = mem3;
                              if(ring_mode != 3) {
                                ladder_radius = ring_width * 1.5;
                              } else {
                                ladder_radius = ring_width * 3;
                              }
                            }
                          } else {
                            nbr[3] = neighbor[mem3] + 1;
                            while((mem4 = neighbor[nbr[3]]) >= 0) {
                              if((mem4 != mem2) && (mem4 != mem1) && (mem4 != mem0) &&
                                 (atomInfo[mem4].protons == cAN_C)) {

                                nbr[4] = neighbor[mem4] + 1;
                                while((mem5 = neighbor[nbr[4]]) >= 0) {
                                  if((mem5 != mem3) && (mem5 != mem2) && (mem5 != mem1)
                                     && (mem5 != mem0)
                                     && (atomInfo[mem5].protons == cAN_N)
                                     && (marked[mem5])) {
                                    /* must be in a mapped ring */
                                    /* clear flag here */
                                    purine_flag = false;

                                    nbr[5] = neighbor[mem5] + 1;
                                    while((mem6 = neighbor[nbr[5]]) >= 0) {
                                      if((mem6 != mem4) && (mem6 != mem3)
                                         && (mem6 != mem2) && (mem6 != mem1)
                                         && (mem6 != mem0)
                                         && (atomInfo[mem6].protons == cAN_C)
                                         && (marked[mem6])) {
                                        nbr[6] = neighbor[mem6] + 1;
                                        while((mem7 = neighbor[nbr[6]]) >= 0) {
                                          ai2 = atomInfo + mem7;
                                          if((mem7 != mem5) && (mem7 != mem4)
                                             && (mem7 != mem3) && (mem7 != mem2)
                                             && (mem7 != mem2) && (mem7 != mem1)
                                             && (mem7 != mem0) && (ai2->protons == cAN_N)
                                             && (marked[mem7])) {
                                            if(WordMatchExact(G, "N1", ai2->name, 1)) {
                                              /* and set flag */
                                              base_at = mem7;
                                              purine_flag = true;
                                            }
                                          }
                                          nbr[6] += 2;
                                        }
                                      }
                                      nbr[5] += 2;
                                    }

                                    if(!purine_flag) {
                                      ai2 = atomInfo + mem5;
                                      if(marked[mem5]
                                         && WordMatchExact(G, "N3", ai2->name, 1)) {
                                        base_at = mem5;
                                      }
                                    }
                                  }
                                  nbr[4] += 2;
                                }
                              }
                              nbr[3] += 2;
                            }
                          }
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              }
              nbr[0] += 2;
            }
          }
        }
      } else if((!ring_mode) && (n_atom == 6)) {        /* going from the base to the sugar */

        for(i = 0; i < n_atom; i++) {
          a1 = atix[i];
          if(!nf)
            nf = nuc_flag[a1];
          ai = atomInfo + a1;

          /* base-hunting */

          if((ai->protons == cAN_N) && (!marked[a1]) &&
             (WordMatchExact(G, "N1", ai->name, 1) ||
              WordMatchExact(G, "N3", ai->name, 1))) {
            mem0 = a1;
            nbr[0] = neighbor[mem0] + 1;
            while((mem1 = neighbor[nbr[0]]) >= 0) {
              if((atomInfo[mem1].protons == cAN_C) && (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while((mem2 = neighbor[nbr[1]]) >= 0) {
                  if((mem2 != mem0) && (!marked[mem2]) &&
                     (atomInfo[mem2].protons == cAN_N)) {
                    nbr[2] = neighbor[mem2] + 1;
                    while((mem3 = neighbor[nbr[2]]) >= 0) {
                      if((mem3 != mem1) && (mem3 != mem0) &&
                         (atomInfo[mem3].protons == cAN_C)) {
                        nbr[3] = neighbor[mem3] + 1;
                        while((mem4 = neighbor[nbr[3]]) >= 0) {
                          if((mem4 != mem2) && (mem4 != mem1) && (mem4 != mem0)) {
			    
			    if((atomInfo[mem4].protons == cAN_N)||
			       (WordMatchExact(G, "C5", atomInfo[mem4].name,1))) {     /* purine case */
                              nbr[4] = neighbor[mem4] + 1;
                              while((mem5 = neighbor[nbr[4]]) >= 0) {
                                if((mem5 != mem3) && (mem5 != mem2) && (mem5 != mem1)
                                   && (mem5 != mem0) && (marked[mem5])
                                   && (atomInfo[mem5].protons == cAN_C)) {
                                  nbr[5] = neighbor[mem5] + 1;
                                  while((mem6 = neighbor[nbr[5]]) >= 0) {
                                    if((mem6 != mem4) && (mem6 != mem3) && (mem6 != mem2)
                                       && (mem6 != mem1) && (mem6 != mem0)
                                       && ((atomInfo[mem6].protons == cAN_C)
                                           || (atomInfo[mem6].protons == cAN_O))
                                       && (marked[mem6])) {
                                      nbr[6] = neighbor[mem6] + 1;
                                      while((mem7 = neighbor[nbr[6]]) >= 0) {
                                        ai2 = atomInfo + mem7;
                                        if((mem7 != mem5) && (mem7 != mem4)
                                           && (mem7 != mem3) && (mem7 != mem2)
                                           && (mem7 != mem2) && (mem7 != mem1)
                                           && (mem7 != mem0) && (ai2->protons == cAN_C)
                                           && (marked[mem7])) {
                                          if(WordMatchExact
                                             (G, NUCLEIC_NORMAL1, ai2->name, 1)
                                             || WordMatchExact(G, NUCLEIC_NORMAL2,
                                                               ai2->name, 1)) {
                                            base_at = a1;
                                            sugar_at = mem7;
                                          }
                                        }
                                        nbr[6] += 2;
                                      }
                                    }
                                    nbr[5] += 2;
                                  }
                                }
                                nbr[4] += 2;
                              }
                            } else if(((atomInfo[mem4].protons == cAN_C)
                                       || (atomInfo[mem4].protons == cAN_O))
                                      && (marked[mem4])) {
                              /* pyrimidine case */

                              nbr[4] = neighbor[mem4] + 1;
                              while((mem5 = neighbor[nbr[4]]) >= 0) {
                                ai2 = atomInfo + mem5;
                                if((mem5 != mem3) && (mem5 != mem2) && (mem5 != mem1)
                                   && (mem5 != mem0) && (ai2->protons == cAN_C)
                                   && (marked[mem5])) {
                                  if(WordMatchExact(G, NUCLEIC_NORMAL1, ai2->name, 1)
                                     || WordMatchExact(G, NUCLEIC_NORMAL2, ai2->name, 1)) {
                                    base_at = a1;
                                    sugar_at = mem5;
                                  }
                                }
                                nbr[4] += 2;
                              }
                            }
                          }
                          nbr[3] += 2;
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              }
              nbr[0] += 2;
            }
          }
        }
      }
      if(n_atom == 6) {

        for(i = 0; i < n_atom; i++) {
          a1 = atix[i];
          ai = atomInfo + a1;

          /* glycosylation hunting */

          if((ai->protons == cAN_C) && (!marked[a1])) { /* unmarked C in ring */
            mem0 = a1;
            nbr[0] = neighbor[mem0] + 1;
            while((mem1 = neighbor[nbr[0]]) >= 0) {
              if((atomInfo[mem1].protons == cAN_C) &&   /* exocyclic C */
                 (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while((mem2 = neighbor[nbr[1]]) >= 0) {
                  if((mem2 != mem0) && (!marked[mem2])
                     && (atomInfo[mem2].protons == cAN_O)) {
                    /* exocyclic O */
                    nbr[2] = neighbor[mem2] + 1;
                    while((mem3 = neighbor[nbr[2]]) >= 0) {
                      if((mem3 != mem1) && (mem3 != mem0) && (marked[mem3])
                         && (atomInfo[mem3].protons == cAN_C)) {
                        /* cyclic C */
                        if(WordMatchExact(G, "C5", ai->name, 1) &&
                           WordMatchExact(G, "C6", atomInfo[mem1].name, 1)) {
                          c5_linked = mem3;
                          c5 = mem0;
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              } else if((atomInfo[mem1].protons == cAN_O) &&    /* exocyclic O */
                        (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while((mem2 = neighbor[nbr[1]]) >= 0) {
                  if((mem2 != mem0) && (marked[mem2])
                     && (atomInfo[mem2].protons == cAN_C)) {
                    /* cyclic C */
                    if(WordMatchExact(G, "C1", ai->name, 1)) {
                      c1_linked = mem2;
                      c1 = mem0;
                    } else if(WordMatchExact(G, "C2", ai->name, 1)) {
                      c2_linked = mem2;
                      c2 = mem0;
                    } else if(WordMatchExact(G, "C3", ai->name, 1)) {
                      c3_linked = mem2;
                      c3 = mem0;
                    } else if(WordMatchExact(G, "C4", ai->name, 1)) {
                      c4_linked = mem2;
                      c4 = mem0;
                    }
                  } else if((mem2 != mem0) && (!marked[mem2])
                            && (atomInfo[mem2].protons == cAN_C)) {
                    /* exocyclic C */
                    nbr[2] = neighbor[mem2] + 1;
                    while((mem3 = neighbor[nbr[2]]) >= 0) {
                      if((mem3 != mem1) && (mem3 != mem0) && (marked[mem3])
                         && (atomInfo[mem3].protons == cAN_C)) {
                        /* cyclic C */
                        if(WordMatchExact(G, "C5", atomInfo[mem3].name, 1) &&
                           WordMatchExact(G, "C6", atomInfo[mem2].name, 1)) {
                          c5 = mem0;
                          c5_linked = mem3;
                        }
                      } else if((mem3 != mem1) && (mem3 != mem0) && (!marked[mem3])
                                && (atomInfo[mem3].protons == cAN_C)) {
                        /* exocyclic */
                        if(WordMatchExact(G, "C1", ai->name, 1) &&
                           WordMatchExact(G, "CA", atomInfo[mem3].name, 1)) {
                          c1 = mem0;
                          c1_linked = mem3;
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              } else if((atomInfo[mem1].protons == cAN_N) &&    /* exocyclic N */
                        (!marked[mem1])) {
                nbr[1] = neighbor[mem1] + 1;
                while((mem2 = neighbor[nbr[1]]) >= 0) {
                  if((mem2 != mem0) && (!marked[mem2])
                     && (atomInfo[mem2].protons == cAN_C)) {
                    /* exocyclic C */
                    nbr[2] = neighbor[mem2] + 1;
                    while((mem3 = neighbor[nbr[2]]) >= 0) {
                      if((mem3 != mem1) && (mem3 != mem0) && (!marked[mem3])
                         && (atomInfo[mem3].protons == cAN_C)) {
                        /* exocyclic */
                        nbr[3] = neighbor[mem3] + 1;
                        while((mem4 = neighbor[nbr[3]]) >= 0) {
                          if((mem4 != mem2) && (mem3 != mem1) && (!marked[mem4])
                             && (atomInfo[mem4].protons == cAN_C)) {
                            /* exocyclic */
                            if(WordMatchExact(G, "C1", ai->name, 1) &&
                               WordMatchExact(G, "CA", atomInfo[mem4].name, 1)) {
                              c1 = mem0;
                              c1_linked = mem4;
                            }
                          }
                          nbr[3] += 2;
                        }
                      }
                      nbr[2] += 2;
                    }
                  }
                  nbr[1] += 2;
                }
              }
              nbr[0] += 2;
            }
          }
        }
      }

      if(sugar_at >= 0) {       /* still need to find the phosphates... */
        int c3_index = -1;
        ai = atomInfo + sugar_at;
        if(WordMatchExact(G, "C3*", ai->name, 1) || WordMatchExact(G, "C3'", ai->name, 1)) {
          c3_index = sugar_at;
        } else {
          mem0 = sugar_at;
          nbr[0] = neighbor[mem0] + 1;
          while((mem1 = neighbor[nbr[0]]) >= 0) {
            if((atomInfo[mem1].protons == cAN_C) && (!marked[mem1])) {
              ai = atomInfo + mem1;
              if(!(WordMatchExact(G, "C3*", ai->name, 1) ||
                   WordMatchExact(G, "C3'", ai->name, 1))) {
                c3_index = mem1;
              }
            }
            nbr[0] += 2;
          }
        }

        if(c3_index >= 0) {     /* now we know where we are... */

          mem0 = c3_index;
          nbr[0] = neighbor[mem0] + 1;
          while((mem1 = neighbor[nbr[0]]) >= 0) {
            if((atomInfo[mem1].protons == cAN_O) && (!marked[mem1])) {
              ai = atomInfo + mem1;
              if(WordMatchExact(G, "O3*", ai->name, 1) ||
                 WordMatchExact(G, "O3'", ai->name, 1))
                o3_at = mem1;
              nbr[1] = neighbor[mem1] + 1;
              while((mem2 = neighbor[nbr[1]]) >= 0) {
                if((mem2 != mem0) && (!marked[mem2]) && (atomInfo[mem2].protons == cAN_P)) {
                  phos3_at = mem2;
                }
                nbr[1] += 2;
              }
            }
            if((atomInfo[mem1].protons == cAN_C) && (marked[mem1])) {

              nbr[1] = neighbor[mem1] + 1;
              while((mem2 = neighbor[nbr[1]]) >= 0) {

                if((mem2 != mem0) && (!marked[mem2]) && (atomInfo[mem2].protons == cAN_C)) {
                  ai = atomInfo + mem2;
                  if(WordMatchExact(G, "C1*", ai->name, 1) ||
                     WordMatchExact(G, "C1'", ai->name, 1))
                    c1_at = mem2;

                  nbr[2] = neighbor[mem2] + 1;
                  while((mem3 = neighbor[nbr[2]]) >= 0) {
                    if((mem3 != mem1) && (mem3 != mem0) &&
                       (atomInfo[mem3].protons == cAN_O) && (!marked[mem3])) {
                      ai = atomInfo + mem3;
                      if(WordMatchExact(G, "O5*", ai->name, 1) ||
                         WordMatchExact(G, "O5'", ai->name, 1))
                        o5_at = mem3;
                      nbr[3] = neighbor[mem3] + 1;
                      while((mem4 = neighbor[nbr[3]]) >= 0) {
                        if((mem4 != mem2) && (!marked[mem4]) &&
                           (atomInfo[mem4].protons == cAN_P)) {
                          phos5_at = mem4;
                        }
                        nbr[3] += 2;
                      }
                    }
                    nbr[2] += 2;
                  }
                }
                nbr[1] += 2;
              }
            }
            nbr[0] += 2;
          }
        }
      }
      /* glycosylation connectors */
      if(1) {
        if((ring_mode) && (finder >= 3)) {
          if(finder >= 3) {
            int b;
            float glyco_radius;
            switch (ring_mode) {
            case 3:
            case 4:
              glyco_radius = ring_width * 3;
              break;
            case 5:
              glyco_radius = ladder_radius;
              break;
            default:
              glyco_radius = ring_width * 1.5;
              break;
            }

            if((alpha != 1.0F) || (ring_alpha != alpha))
              CGOAlpha(cgo, alpha);

            for(b = 0; b < 5; b++) {
              int g1 = -1, g2 = -1;
              switch (b) {
              case 0:
                g1 = c1;
                g2 = c1_linked;
                break;
              case 1:
                g1 = c2;
                g2 = c2_linked;
                break;
              case 2:
                g1 = c3;
                g2 = c3_linked;
                break;
              case 3:
                g1 = c4;
                g2 = c4_linked;
                break;
              case 4:
                g1 = c5;
                g2 = c5_linked;
                break;
              }

              if((g1 >= 0) && (g2 >= 0)) {
                AtomInfoType *g1_ai = atomInfo + g1;
                AtomInfoType *g2_ai = atomInfo + g2;
                if((g1_ai->visRep[cRepCartoon]) &&
                   (g2_ai->visRep[cRepCartoon]) &&
                   ((!sc_helper) || !(g1_ai->visRep[cRepLine] ||
                                      g1_ai->visRep[cRepCyl] ||
                                      g1_ai->visRep[cRepSphere]))) {

                  float *g1p, *g2p;
                  float *color;
                  float avg[3];

                  {
                    int g1_x, g2_x;

                    if(obj->DiscreteFlag) {
                      if(cs == obj->DiscreteCSet[g1] && cs == obj->DiscreteCSet[g2]) {
                        g1_x = obj->DiscreteAtmToIdx[g1];
                        g2_x = obj->DiscreteAtmToIdx[g2];
                      } else {
                        g1_x = -1;
                        g2_x = -1;
                      }
                    } else {
                      g1_x = cs->AtmToIdx[g1];
                      g2_x = cs->AtmToIdx[g2];
                    }
                    g1p = cs->Coord + 3 * g1_x;
                    g2p = cs->Coord + 3 * g2_x;
                  }

                  if(!((!((ring_mode == 0) || (ring_mode == 4) || (ring_mode == 5))) ||
                       (!marked[g2]))) {
                    g2p = moved + 3 * g2;
                  }
                  if((ring_mode == 0) || (ring_mode == 4) || (ring_mode == 5)) {
                    /* ring center */
                    int i;
                    /* compute average coordinate and mark atoms so that ring is only drawn once */
                    zero3f(avg);
                    for(i = 0; i < n_atom; i++) {
                      add3f(avg, v_i[i], avg);
                    }
                    scale3f(avg, 1.0F / n_atom, avg);
                    g1p = avg;
                  }

                  if(!g1_ai->masked)
                    CGOPickColor(cgo, g1, cPickableAtom);

                  if(ladder_color >= 0) {
                    color = ColorGet(G, ladder_color);
                    CGOCustomCylinderv(cgo, g1p,
                                       g2p, glyco_radius, color, color, 2.0F, 2.0F);
                  } else {
                    CGOCustomCylinderv(cgo, g1p,
                                       g2p,
                                       glyco_radius,
                                       ColorGet(G, g1_ai->color),
                                       ColorGet(G, g2_ai->color), 2.0F, 2.0F);
                  }
                }
              }
            }
          }
        }
      }

      /* see if any of the neighbors are confirmed nucleic acids... */
      if(sugar_at >= 0) {
        if(!nf) {
          mem0 = sugar_at;
          nbr[0] = neighbor[mem0] + 1;
          while((!nf) && (mem1 = neighbor[nbr[0]]) >= 0) {
            if(!nf)
              nf = nuc_flag[mem1];
            nbr[1] = neighbor[mem1] + 1;
            while((!nf) && (mem2 = neighbor[nbr[1]]) >= 0) {
              if(mem2 != mem0) {
                if(!nf)
                  nf = nuc_flag[mem2];
                nbr[2] = neighbor[mem2] + 1;
                while((!nf) && (mem3 = neighbor[nbr[2]]) >= 0) {
                  if((mem3 != mem1) && (mem3 != mem0)) {
                    if(!nf)
                      nf = nuc_flag[mem3];
                    nbr[3] = neighbor[mem3] + 1;
                    while((mem4 = neighbor[nbr[3]]) >= 0) {
                      if(mem4 != mem2) {
                        if(!nf)
                          nf = nuc_flag[mem4];
                        nbr[4] = neighbor[mem4] + 1;
                        while((mem5 = neighbor[nbr[4]]) >= 0) {
                          if(mem5 != mem3) {
                            if(!nf)
                              nf = nuc_flag[mem5];
                          }
                          nbr[4] += 2;
                        }
                      }
                      nbr[3] += 2;
                    }
                  }
                  nbr[2] += 2;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }
        }

        if(nf) {
          if((ring_mode) && ((finder == 1) || (finder >= 3))) {
            if((c1_at >= 0) && (base_at >= 0)) {
              int save_at = sugar_at;
              sugar_at = c1_at;
              {
                AtomInfoType *sug_ai = atomInfo + sugar_at;
                AtomInfoType *bas_ai = atomInfo + base_at;
                if((sug_ai->visRep[cRepCartoon]) &&
                   (bas_ai->visRep[cRepCartoon]) &&
                   ((!sc_helper) || !(bas_ai->visRep[cRepLine] ||
                                      bas_ai->visRep[cRepCyl] ||
                                      bas_ai->visRep[cRepSphere]))) {

                  int sug, bas;
                  float *color;
                  if(obj->DiscreteFlag) {
                    if(cs == obj->DiscreteCSet[sugar_at] &&
                       cs == obj->DiscreteCSet[base_at]) {
                      sug = obj->DiscreteAtmToIdx[sugar_at];
                      bas = obj->DiscreteAtmToIdx[base_at];
                    } else {
                      sug = -1;
                      bas = -1;
                    }
                  } else {
                    sug = cs->AtmToIdx[sugar_at];
                    bas = cs->AtmToIdx[base_at];
                  }

                  if((sug >= 0) && (bas >= 0)) {
                    if(!bas_ai->masked)
                      CGOPickColor(cgo, base_at, cPickableAtom);
                    if(ladder_color >= 0) {
                      color = ColorGet(G, ladder_color);
                      CGOCustomCylinderv(cgo, cs->Coord + 3 * sug,
                                         cs->Coord + 3 * bas,
                                         ladder_radius, color, color, 2.0F, 2.0F);
                    } else {
                      CGOCustomCylinderv(cgo, cs->Coord + 3 * sug,
                                         cs->Coord + 3 * bas,
                                         ladder_radius,
                                         ColorGet(G, sug_ai->color),
                                         ColorGet(G, bas_ai->color), 2.0F, 2.0F);
                    }
                  }
                }
              }
              base_at = save_at;
              sugar_at = save_at;
            }
          }
          if((base_at >= 0) && (sugar_at >= 0)) {
            AtomInfoType *sug_ai = atomInfo + sugar_at;
            AtomInfoType *bas_ai = atomInfo + base_at;
            if((sug_ai->visRep[cRepCartoon]) &&
               (bas_ai->visRep[cRepCartoon]) &&
               ((!sc_helper) || !(bas_ai->visRep[cRepLine] ||
                                  bas_ai->visRep[cRepCyl] ||
                                  bas_ai->visRep[cRepSphere]))) {

              int sug, bas;
              float *color;
              float *v_outer, tmp[3], outer[3];
              if(obj->DiscreteFlag) {
                if(cs == obj->DiscreteCSet[sugar_at] && cs == obj->DiscreteCSet[base_at]) {
                  sug = obj->DiscreteAtmToIdx[sugar_at];
                  bas = obj->DiscreteAtmToIdx[base_at];
                } else {
                  sug = -1;
                  bas = -1;
                }
              } else {
                sug = cs->AtmToIdx[sugar_at];
                bas = cs->AtmToIdx[base_at];
              }

              if((sug >= 0) && (bas >= 0)) {
                int p3, p5;
                v_outer = cs->Coord + 3 * sug;

                if((o3_at >= 0) && (phos3_at < 0))
                  phos3_at = o3_at;
                if((o5_at >= 0) && (phos5_at < 0))
                  phos5_at = o5_at;
                if((na_mode != 1) && (phos3_at >= 0) && (phos5_at >= 0)) {

                  if(obj->DiscreteFlag) {
                    if(cs == obj->DiscreteCSet[phos3_at] &&
                       cs == obj->DiscreteCSet[phos5_at]) {
                      p3 = obj->DiscreteAtmToIdx[phos3_at];
                      p5 = obj->DiscreteAtmToIdx[phos5_at];
                    } else {
                      p3 = -1;
                      p5 = -1;
                    }
                  } else {
                    p3 = cs->AtmToIdx[phos3_at];
                    p5 = cs->AtmToIdx[phos5_at];
                  }
                  if((p3 >= 0) && (p5 >= 0)) {
                    if(ring_mode) {
                      scale3f(cs->Coord + 3 * p5, 0.333333F, outer);
                      scale3f(cs->Coord + 3 * p3, 0.666667F, tmp);
                    } else {
                      scale3f(cs->Coord + 3 * p3, 0.5F, outer);
                      scale3f(cs->Coord + 3 * p5, 0.5F, tmp);
                    }
                    add3f(tmp, outer, outer);
                    v_outer = outer;
                  }
                }
                if(!bas_ai->masked)
                  CGOPickColor(cgo, base_at, cPickableAtom);
                if(ladder_color >= 0) {
                  color = ColorGet(G, ladder_color);
                  CGOCustomCylinderv(cgo, v_outer,
                                     cs->Coord + 3 * bas,
                                     ladder_radius, color, color, 2.0F, 2.0F);
                } else {
                  CGOCustomCylinderv(cgo, v_outer,
                                     cs->Coord + 3 * bas,
                                     ladder_radius,
                                     ColorGet(G, sug_ai->color),
                                     ColorGet(G, bas_ai->color), 2.0F, 2.0F);
                }
              }
            }
          }
        }
      }
    }
    if((!ring_mode) || (finder == 2)) {
      if(ladder_mode) {         /* mark sure all rings traversed are mark */
        int i;
        for(i = 0; i <= n_atom; i++) {
          marked[atix[i]] = true;
        }
      }
    }
    if((!nf) && ((have_C4_prime >= 0) || (have_C4 >= 0))) {
      int nbr[9];
      register int *neighbor = obj->Neighbor;
      register int mem0, mem1, mem2, mem3, mem4, mem5, mem6, mem7, mem8, mem9;
      /* see if any of the neighbors are confirmed nucleic acids... */
      if(have_C4_prime >= 0)
        mem0 = have_C4_prime;
      else if(have_C4 >= 0)
        mem0 = have_C4;
      else
        mem0 = -1;
      if(mem0 >= 1) {
        nbr[0] = neighbor[mem0] + 1;
        while((!nf) && (mem1 = neighbor[nbr[0]]) >= 0) {
          if(!nf)
            nf = nuc_flag[mem1];
          nbr[1] = neighbor[mem1] + 1;
          while((!nf) && (mem2 = neighbor[nbr[1]]) >= 0) {
            if(mem2 != mem0) {
              if(!nf)
                nf = nuc_flag[mem2];
              nbr[2] = neighbor[mem2] + 1;
              while((!nf) && (mem3 = neighbor[nbr[2]]) >= 0) {
                if((mem3 != mem1) && (mem3 != mem0)) {
                  if(!nf)
                    nf = nuc_flag[mem3];
                  nbr[3] = neighbor[mem3] + 1;
                  while((mem4 = neighbor[nbr[3]]) >= 0) {
                    if(mem4 != mem2) {
                      if(!nf)
                        nf = nuc_flag[mem4];
                      nbr[4] = neighbor[mem4] + 1;
                      while((mem5 = neighbor[nbr[4]]) >= 0) {
                        if(mem5 != mem3) {
                          if(!nf)
                            nf = nuc_flag[mem5];
                          nbr[5] = neighbor[mem5] + 1;
                          while((mem6 = neighbor[nbr[5]]) >= 0) {
                            if(mem6 != mem4) {
                              if(!nf)
                                nf = nuc_flag[mem6];
                              nbr[6] = neighbor[mem6] + 1;
                              while((mem7 = neighbor[nbr[6]]) >= 0) {
                                if(mem7 != mem5) {
                                  if(!nf)
                                    nf = nuc_flag[mem7];
                                  nbr[7] = neighbor[mem7] + 1;
                                  while((mem8 = neighbor[nbr[7]]) >= 0) {
                                    if(mem8 != mem6) {
                                      if(!nf)
                                        nf = nuc_flag[mem8];
                                      nbr[8] = neighbor[mem8] + 1;
                                      while((mem9 = neighbor[nbr[8]]) >= 0) {
                                        if(mem9 != mem7) {
                                          if(!nf)
                                            nf = nuc_flag[mem9];
                                        }
                                        nbr[8] += 2;
                                      }
                                    }
                                    nbr[7] += 2;
                                  }
                                }
                                nbr[6] += 2;
                              }
                            }
                            nbr[5] += 2;
                          }
                        }
                        nbr[4] += 2;
                      }
                    }
                    nbr[3] += 2;
                  }
                }
                nbr[2] += 2;
              }
            }
            nbr[1] += 2;
          }
          nbr[0] += 2;
        }
      }
    }
    if(n_atom) {                /* store center of ring */

      float avg[3];
      float avg_col[3];
      int i;
      float up[3], upi[3];
      float vc0[3], vc1[3];
      float *color = NULL;
      /* compute average coordinate and mark atoms so that ring is only drawn once */
      zero3f(avg);
      zero3f(avg_col);

      for(i = 0; i < n_atom; i++) {
        add3f(avg, v_i[i], avg);
        add3f(avg_col, col[i], avg_col);
        marked[atix[i]] = true;
      }

      scale3f(avg, 1.0F / n_atom, avg);
      scale3f(avg_col, 1.0F / n_atom, avg_col);

      for(i = 0; i < n_atom; i++) {
        float *v_moved = moved + atix[i] * 3;
        copy3f(avg, v_moved);   /* store ring center for later use */
      }

      if((nf || (!ladder_mode) || (finder >= 3)) &&
         ring_mode &&
         (((finder == 1) && ((have_C4 >= 0) || (have_C4_prime >= 0))) ||
          ((finder == 2) && ((have_C4 >= 0))) ||
          ((finder == 3) && ((have_C_number >= 0))) || ((finder == 4)))) {

        if((alpha != 1.0F) || (ring_alpha != alpha))
          CGOAlpha(cgo, ring_alpha);

        if(ring_color >= 0) {
          color = ColorGet(G, ring_color);
        } else {
          color = avg_col;
        }

        CGOColorv(cgo, color);

        if((ring_mode == 4) || (ring_mode == 5)) {      /* spherical ring  */
          float radius = ring_radius;
          if(radius < 0.0F) {
            radius = 0.0F;
            if(ring_mode == 4) {
              for(i = 0; i < n_atom; i++) {
                float dist = diff3f(avg, v_i[i]);
                if(radius < dist)
                  radius = dist;
              }
            } else {
              radius = ladder_radius;
            }
          }

          if(n_atom) {
            add3f(avg_col, col[i], avg_col);
            if(!ai_i[0]->masked)
              CGOPickColor(cgo, atix[0], cPickableAtom);
            CGOSphere(cgo, avg, radius);
          }

        } else {

          /* clear the normals */

          for(i = 0; i <= n_atom; i++) {
            zero3f(n_up[i]);
            zero3f(n_dn[i]);
          }

          /* compute average normals */

          {
            float acc[3];
            int ii;

            zero3f(acc);
            for(i = 0; i < n_atom; i++) {
              ii = i + 1;
              subtract3f(v_i[i], avg, vc0);
              subtract3f(v_i[ii], avg, vc1);
              cross_product3f(vc0, vc1, up);
              add3f(up, n_up[i], n_up[i]);
              add3f(up, n_up[ii], n_dn[ii]);
              if(!i) {
                add3f(up, n_up[n_atom], n_up[n_atom]);
              } else if(ii == n_atom) {
                add3f(up, n_up[0], n_up[0]);
              }
              add3f(up, acc, acc);
            }
            normalize3f(up);
            scale3f(up, -1.0F, upi);
          }

          for(i = 0; i <= n_atom; i++) {
            normalize3f(n_up[i]);
            scale3f(n_up[i], -1.0F, n_dn[i]);
          }

          {
            int ii;
            float mid[3];
            float up_add[3];
            float ct[3], cb[3];
            float v0t[3], v0b[3];
            float v1t[3], v1b[3];
            float out[3];

            CGOBegin(cgo, GL_TRIANGLES);

            for(i = 0; i < n_atom; i++) {
              ii = i + 1;
              average3f(v_i[ii], v_i[i], mid);

              subtract3f(mid, avg, out);        /* compute outward-facing normal */
              normalize3f(out);

              scale3f(up, ring_width, up_add);

              add3f(avg, up_add, ct);
              subtract3f(avg, up_add, cb);

              add3f(v_i[i], up_add, v0t);
              subtract3f(v_i[i], up_add, v0b);

              add3f(v_i[ii], up_add, v1t);
              subtract3f(v_i[ii], up_add, v1b);

              CGONormalv(cgo, up);
              if(ring_color < 0)
                CGOColorv(cgo, color);
              CGOPickColor(cgo, atix[i], cPickableAtom);        /* TODO: masking for cartoons! */
              CGOVertexv(cgo, ct);
              CGONormalv(cgo, n_up[i]);
              if(ring_color < 0)
                CGOColorv(cgo, col[i]);
              CGOPickColor(cgo, atix[i], cPickableAtom);
              CGOVertexv(cgo, v0t);
              CGONormalv(cgo, n_up[ii]);
              if(ring_color < 0)
                CGOColorv(cgo, col[ii]);
              CGOPickColor(cgo, atix[ii], cPickableAtom);
              CGOVertexv(cgo, v1t);

              if(ring_mode > 1) {
                CGONormalv(cgo, out);

                if(ring_color < 0)
                  CGOColorv(cgo, col[i]);
                CGOPickColor(cgo, atix[i], cPickableAtom);
                CGOVertexv(cgo, v0t);
                CGOVertexv(cgo, v0b);
                if(ring_color < 0)
                  CGOColorv(cgo, col[ii]);
                CGOPickColor(cgo, atix[ii], cPickableAtom);
                CGOVertexv(cgo, v1t);
                CGOVertexv(cgo, v1t);
                if(ring_color < 0)
                  CGOColorv(cgo, col[i]);
                CGOPickColor(cgo, atix[i], cPickableAtom);
                CGOVertexv(cgo, v0b);
                if(ring_color < 0)
                  CGOColorv(cgo, col[ii]);
                CGOPickColor(cgo, atix[ii], cPickableAtom);
                CGOVertexv(cgo, v1b);
              }

              CGONormalv(cgo, upi);
              if(ring_color < 0)
                CGOColorv(cgo, color);
              CGOPickColor(cgo, atix[i], cPickableAtom);
              CGOVertexv(cgo, cb);
              CGONormalv(cgo, n_dn[ii]);
              if(ring_color < 0)
                CGOColorv(cgo, col[ii]);
              CGOPickColor(cgo, atix[ii], cPickableAtom);
              CGOVertexv(cgo, v1b);
              CGONormalv(cgo, n_dn[i]);
              if(ring_color < 0)
                CGOColorv(cgo, col[i]);
              CGOPickColor(cgo, atix[i], cPickableAtom);
              CGOVertexv(cgo, v0b);

            }
            CGOEnd(cgo);

            if((alpha != 1.0F) || (ring_alpha != alpha))
              CGOAlpha(cgo, alpha);

            if(ring_mode == 1) {
              for(i = 0; i < n_atom; i++) {
                ii = i + 1;
                CGOPickColor(cgo, atix[i], cPickableAtom);
                if(ring_color < 0) {
                  CGOSausage(cgo, v_i[i], v_i[ii], ring_width, col[i], col[ii]);
                } else {
                  CGOSausage(cgo, v_i[i], v_i[ii], ring_width, color, color);
                }
              }
            } else if(ring_mode == 3) {
              for(i = 0; i < n_atom; i++) {
                ii = i + 1;
                CGOPickColor(cgo, atix[i], cPickableAtom);
                if(ring_color < 0) {
                  CGOSausage(cgo, v_i[i], v_i[ii], 3 * ring_width, col[i], col[ii]);
                } else {
                  CGOSausage(cgo, v_i[i], v_i[ii], 3 * ring_width, color, color);
                }
              }
            }
          }
        }
      }
    }
  }
}

static void nuc_acid(PyMOLGlobals * G, int a, int a1, AtomInfoType * ai, CoordSet * cs,
                     ObjectMolecule * obj, int na_mode, int *nuc_flag, int set_flags,
                     int *p_a2, int *p_nSeg, float **p_v_o_last,
                     int **p_s, int **p_i, int **p_cc,
                     int *p_nAt, int *p_cur_car, int **p_ss, int *p_putty_flag,
                     float **p_v, float **p_vo)
{
  int a2 = *p_a2;
  int nSeg = *p_nSeg;
  float *v_o_last = *p_v_o_last;
  int *s = *p_s;
  int *i = *p_i;
  int *cc = *p_cc;
  int nAt = *p_nAt;
  int cur_car = *p_cur_car;
  int *ss = *p_ss;
  int putty_flag = *p_putty_flag;
  float *vo = *p_vo;
  float *v = *p_v;

  int a3, a4, st, nd;
  float *v_o, *v_c, *v_n, t0[3];
  float *v1;

  if(a2 < 0) {
    nSeg++;
    v_o_last = NULL;
  }
  *(s++) = nSeg;
  nAt++;
  *(i++) = a;
  cur_car = ai->cartoon;
  if(cur_car == cCartoon_auto)
    cur_car = cCartoon_tube;
  *ss = 3;                      /* DNA/RNA */

  if(cur_car == cCartoon_putty)
    putty_flag = true;

  *(cc++) = cur_car;
  v1 = cs->Coord + 3 * a;
  *(v++) = *(v1++);
  *(v++) = *(v1++);
  *(v++) = *(v1++);

  if(a2 >= 0) {
    if(set_flags) {
      if((obj->AtomInfo[a2].protons == cAN_P) && (!nuc_flag[a2])) {
        int *nf = NULL;
        AtomInfoBracketResidueFast(G, obj->AtomInfo, obj->NAtom, a2, &st, &nd);

        nf = nuc_flag + st;
        for(a3 = st; a3 <= nd; a3++) {
          *(nf++) = true;
        }
      }
    } else if((na_mode >= 2) && (!nuc_flag[a2])) {      /* just a single nucleotide -- skip */
      cur_car = cCartoon_skip;
      *(cc - 2) = cCartoon_skip;
      *(cc - 1) = cCartoon_skip;
    }
  }

  a2 = a1;

  ss++;

  v_c = NULL;
  v_n = NULL;
  v_o = NULL;

  AtomInfoBracketResidueFast(G, obj->AtomInfo, obj->NAtom, a1, &st, &nd);

  {
    int *nf = NULL;
    if(set_flags && v_o_last)
      nf = nuc_flag + st;
    for(a3 = st; a3 <= nd; a3++) {
      if(nf)
        *(nf++) = true;         /* mark this residue as being part of a nucleic acid chain */
      if(obj->DiscreteFlag) {
        if(cs == obj->DiscreteCSet[a3])
          a4 = obj->DiscreteAtmToIdx[a3];
        else
          a4 = -1;
      } else
        a4 = cs->AtmToIdx[a3];
      if(a4 >= 0) {
        if(na_mode == 1) {
          if(WordMatchExact(G, NUCLEIC_NORMAL1, obj->AtomInfo[a3].name, 1) ||
             WordMatchExact(G, NUCLEIC_NORMAL2, obj->AtomInfo[a3].name, 1)) {
            v_c = cs->Coord + 3 * a4;
          }
        } else if(a3 == a1) {
          v_c = cs->Coord + 3 * a4;
        }
        if(WordMatchExact(G, NUCLEIC_NORMAL0, obj->AtomInfo[a3].name, 1)) {
          v_o = cs->Coord + 3 * a4;
        }
      }
    }
  }
  if(!(v_c && v_o)) {
    zero3f(vo);
    v_o_last = NULL;
  } else {
    if(v_o_last) {
      add3f(v_o, v_o_last, t0);
      add3f(v_o_last, t0, t0);
      scale3f(t0, 0.333333F, t0);
      subtract3f(v_c, t0, vo);
    } else {
      subtract3f(v_c, v_o, vo);
    }
    v_o_last = v_o;
    normalize3f(vo);
  }
  vo += 3;

  *p_a2 = a2;
  *p_nSeg = nSeg;
  *p_v_o_last = v_o_last;
  *p_s = s;
  *p_i = i;
  *p_cc = cc;
  *p_nAt = nAt;
  *p_cur_car = cur_car;
  *p_ss = ss;
  *p_putty_flag = putty_flag;
  *p_vo = vo;
  *p_v = v;
}

Rep *RepCartoonNew(CoordSet * cs, int state)
{
  PyMOLGlobals *G = cs->State.G;
  ObjectMolecule *obj;
  int a, b, c, f, e, a1, a2, c1, c2, i0, *i, *s, *at, *seg, nAt, *atp, a3, a4 =
    0, *car, *cc, *sstype;
  float *v, *v0, *v1, *v2, *v3, *v4, *v5, *vo, *vn, *va;
  float *p0, *p1, *p2, *p3;
  float *pv = NULL;
  float *pvo = NULL, *pva = NULL;
  float *dv = NULL;
  float *nv = NULL;
  float *tv = NULL;
  float *vc = NULL;
  float *tmp = NULL;
  int last, first, end_flag;
  int *vi, atom_index1, atom_index2;
  float f0, f1, f2, f3, f4, dev;
  float *d, dp;
  float *dl = NULL;
  int nSeg;
  int sampling;
  int *ss, *fp;
  float power_a = 5;
  float power_b = 5;
  float loop_radius;
  float tube_radius;
  float putty_radius;

  int visFlag;
  CExtrude *ex = NULL, *ex1;
  int n_p, n_pm1, n_pm2;
  int loop_quality, oval_quality, tube_quality, putty_quality;
  float oval_width, oval_length;
  float dumbbell_radius, dumbbell_width, dumbbell_length;
  float throw;
  int st, nd;
  float *v_c, *v_n, *v_o, *v_o_last = NULL;
  float t0[3], t1[3], t2[3], t3[3], t4[3], o0[12], o1[12];
  float max_dot;
  float length, width;
  int cur_car;
  int contFlag, extrudeFlag;
  int cartoon_debug;
  int fancy_helices;
  int fancy_sheets;
  int refine;
  int contigFlag;
  int discrete_colors;
  int cylindrical_helices;
  int last_color, uniform_color;
  int cartoon_color, highlight_color;
  int cartoon_side_chain_helper;
  int ladder_mode, ladder_color;
  float ladder_radius, ring_radius;
  int round_helices;
  int smooth_loops;
  int na_mode;
  int parity;
  float refine_tips;
  float helix_radius;
  float *h_start = NULL, *h_end = NULL;
  float *sampling_tmp;
  int *flag_tmp;
  int smooth_first, smooth_last, smooth_cycles, flat_cycles;
  int trace, trace_mode;
  int skip_to;
  AtomInfoType *ai, *last_ai = NULL;
  float alpha;
  int putty_flag = false;
  float putty_mean = 10.0F, putty_stdev = 0.0F;
  float putty_max = -FLT_MAX, putty_min = FLT_MAX;
  AtomInfoType *trailing_O3p_ai = NULL;
  int trailing_O3p_a = 0, trailing_O3p_a1 = 0;
  AtomInfoType *leading_O5p_ai = NULL;
  int leading_O5p_a = 0, leading_O5p_a1 = 0;
  int *ring_anchor = NULL;
  int ring_mode, ring_finder, ring_finder_eff;
  int n_ring = 0;
  float ring_width;
  int ring_color;
  int loop_cap, tube_cap;
  int *nuc_flag = NULL;
  int nucleic_color = 0;
  float ring_alpha;

  /* THIS IS BY FAR THE WORST ROUTINE IN PYMOL!
   * DEVELOP ON IT ONLY AT EXTREME RISK TO YOUR MENTAL HEALTH */

  OOAlloc(G, RepCartoon);

  PRINTFD(G, FB_RepCartoon)
    " RepCartoonNew-Debug: entered.\n" ENDFD;

  obj = cs->Obj;
  visFlag = false;
  for(a = 0; a < cs->NIndex; a++) {
    if(obj->AtomInfo[cs->IdxToAtm[a]].visRep[cRepCartoon]) {
      visFlag = true;
      break;
    }
  }
  if(!visFlag) {
    OOFreeP(I);
    return (NULL);              /* skip if not visible */
  }

  RepInit(G, &I->R);
  power_a = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_power);
  power_b = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_power_b);

  cartoon_debug = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_debug);
  length = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_rect_length);
  width = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_rect_width);
  trace = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_trace_atoms);
  trace_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_trace_atoms_mode);

  alpha =
    1.0F - SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_transparency);
  throw = SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_throw);

  sampling = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_sampling);
  if(sampling < 1)
    sampling = 1;
  loop_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_loop_radius);
  if(loop_radius < 0.01F)
    loop_radius = 0.01F;
  loop_quality =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_loop_quality);
  if(loop_quality < 3)
    loop_quality = 3;
  if(SettingGetGlobal_i(G, cSetting_ray_trace_mode) > 0)
    if(loop_quality < 12)
      loop_quality *= 2;

  tube_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_tube_radius);
  if(tube_radius < 0.01F)
    tube_radius = 0.01F;
  tube_quality =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_tube_quality);
  if(tube_quality < 3)
    tube_quality = 3;

  putty_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_putty_radius);
  /* WLD removed: if(putty_radius<0.01F) putty_radius=0.01F; --
     should not constrain what is effectively a scale factor */

  putty_quality =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_putty_quality);
  if(putty_quality < 3)
    putty_quality = 3;

  cartoon_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_color);
  cartoon_side_chain_helper =
    SettingGet_b(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_side_chain_helper);

  highlight_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_highlight_color);

  oval_length =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_oval_length);
  oval_width =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_oval_width);
  oval_quality =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_oval_quality);
  if(oval_quality < 3)
    tube_quality = 3;

  dumbbell_length =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_dumbbell_length);
  dumbbell_width =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_dumbbell_width);
  dumbbell_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_dumbbell_radius);
  if(dumbbell_radius < 0.01F)
    dumbbell_radius = 0.01F;

  fancy_helices =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_fancy_helices);
  fancy_sheets =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_fancy_sheets);
  cylindrical_helices =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_cylindrical_helices);
  refine = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_refine);
  refine_tips =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_refine_tips);

  discrete_colors =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_discrete_colors);
  round_helices =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_round_helices);
  smooth_loops =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_smooth_loops);
  helix_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_helix_radius);

  ring_width =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ring_width);
  if(ring_width < 0.0F) {
    ring_width =
      fabs(SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_stick_radius)) * 0.5F;
  }

  ring_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ring_color);

  if(ring_color == -1)
    ring_color = cartoon_color;

  smooth_first =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_smooth_first);
  smooth_last =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_smooth_last);
  smooth_cycles =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_smooth_cycles);
  flat_cycles =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_flat_cycles);

  na_mode =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_nucleic_acid_mode);
  nucleic_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting,
                     cSetting_cartoon_nucleic_acid_color);

  if(nucleic_color == -1)
    nucleic_color = cartoon_color;

  ring_mode = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ring_mode);
  ring_finder =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ring_finder);
  ring_alpha =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ring_transparency);

  if(ring_alpha < 0.0F)
    ring_alpha = alpha;
  else
    ring_alpha = 1.0F - ring_alpha;

  alpha =
    1.0F - SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_transparency);
  ladder_mode =
    SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ladder_mode);
  ladder_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ladder_radius);
  ladder_color =
    SettingGet_color(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ladder_color);
  ring_radius =
    SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_ring_radius);

  if(ladder_color == -1)
    ladder_color = cartoon_color;

  ring_finder_eff = ring_finder;
  if((!ring_mode) || (ring_finder == 2))
    ring_finder_eff = 1;

  tube_cap = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_tube_cap);
  loop_cap = SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_loop_cap);

  I->R.fRender = (void (*)(struct Rep *, RenderInfo *)) RepCartoonRender;
  I->R.fFree = (void (*)(struct Rep *)) RepCartoonFree;
  I->R.fRecolor = NULL;
  I->R.obj = &obj->Obj;
  I->R.cs = cs;
  I->ray = NULL;
  I->std = NULL;
  I->R.context.object = (void *) obj;
  I->R.context.state = state;

  /* find all of the CA points */

  at = Alloc(int, cs->NAtIndex);        /* cs index pointers */
  pv = Alloc(float, cs->NAtIndex * 3);
  tmp = Alloc(float, cs->NAtIndex * 3);
  pvo = Alloc(float, cs->NAtIndex * 3); /* orientation vector */
  pva = Alloc(float, cs->NAtIndex * 6); /* alternative orientation vectors, two per atom */
  seg = Alloc(int, cs->NAtIndex);
  car = Alloc(int, cs->NAtIndex);
  sstype = Alloc(int, cs->NAtIndex);
  sampling_tmp = Alloc(float, sampling * 3);
  flag_tmp = Calloc(int, cs->NAtIndex);
  nuc_flag = Calloc(int, cs->NAtIndex);
  if(ring_mode || ladder_mode) {
    ring_anchor = VLAlloc(int, cs->NAtIndex / 10 + 1);
  }
  i = at;
  v = pv;
  vo = pvo;
  s = seg;
  cc = car;
  ss = sstype;
  fp = flag_tmp;
  nAt = 0;
  nSeg = 0;
  a2 = -1;
  parity = 1;

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

      if(ai->visRep[cRepCartoon]) {
        if(ring_anchor && (ai->protons != cAN_H) && ((ring_finder_eff >= 3) ||  /* all 5-7 atom rings */
                                                     ((ring_finder_eff <= 2) && /*  C4-containing rings */
                                                      (WordMatchExact
                                                       (G, "C4", ai->name, 1)))
                                                     || ((ring_finder_eff == 1)
                                                         &&
                                                         ((WordMatchExact
                                                           (G, "C4*", ai->name, 1)
                                                           || WordMatchExact(G, "C4'",
                                                                             ai->name,
                                                                             1)))))) {
          VLACheck(ring_anchor, int, n_ring);
          ring_anchor[n_ring] = a1;
          n_ring++;
        }
        /*                        if(!obj->AtomInfo[a1].hetatm) */
        if((!ai->alt[0]) || (ai->alt[0] == 'A')) {
          if(trace || (((ai->protons == cAN_C) &&
                        (WordMatch(G, "CA", ai->name, 1) < 0)) &&
                       !AtomInfoSameResidueP(G, last_ai, ai))) {
            PRINTFD(G, FB_RepCartoon)
              " RepCartoon: found CA in %s; a2 %d\n", ai->resi, a2 ENDFD;

            if(trailing_O3p_ai && ((na_mode == 2) || (na_mode == 4))) {

              /*  3' nucleic acid cap */

              nuc_acid(G, trailing_O3p_a, trailing_O3p_a1, trailing_O3p_ai,
                       cs, obj, na_mode, nuc_flag, false, &a2, &nSeg, &v_o_last, &s, &i,
                       &cc, &nAt, &cur_car, &ss, &putty_flag, &v, &vo);
              a2 = -1;
              trailing_O3p_ai = NULL;
            }

            if(!trace) {
              if(a2 >= 0) {
                /*
                   if((abs(obj->AtomInfo[a1].resv-obj->AtomInfo[a2].resv)>1)||
                   (obj->AtomInfo[a1].chain[0]!=obj->AtomInfo[a2].chain[0])||
                   (!WordMatch(G,obj->AtomInfo[a1].segi,obj->AtomInfo[a2].segi,1))) */
                if(!ObjectMoleculeCheckBondSep(obj, a1, a2, 3)) /* CA->N->C->CA = 3 bonds */
                  a2 = -1;

              }
            } else {
              if(a2 >= 0) {
                if(!AtomInfoSequential
                   (G, obj->AtomInfo + a2, obj->AtomInfo + a1, trace_mode))
                  a2 = -1;
              }
            }
            last_ai = ai;

            PRINTFD(G, FB_RepCartoon)
              " RepCartoon: found CA in %s; a2 %d\n", ai->resi, a2 ENDFD;

            if(a2 < 0)
              nSeg++;
            *(s++) = nSeg;
            nAt++;
            *(i++) = a;
            cur_car = ai->cartoon;

            *fp = ai->flags;    /* store atom flags */

            if(cartoon_side_chain_helper) {
              if(ai->visRep[cRepLine] || ai->visRep[cRepCyl] || ai->visRep[cRepSphere])
                *fp |= cAtomFlag_no_smooth;
            }

            switch (ai->ssType[0]) {
            case 'H':
            case 'h':
              if(cur_car == cCartoon_auto) {
                if(cylindrical_helices)
                  cur_car = cCartoon_skip_helix;
                else if(fancy_helices)
                  cur_car = cCartoon_dumbbell;
                else
                  cur_car = cCartoon_oval;
              }
              *ss = 1;          /* helix */
              parity = 0;
              break;
            case 'S':
            case 's':
              if(cur_car == cCartoon_auto) {
                if(fancy_sheets)
                  cur_car = cCartoon_arrow;
                else
                  cur_car = cCartoon_rect;
              }
              *ss = 2;
              parity = !parity;
              break;
            default:           /* 'L', 'T', 0, etc. */
              if(cur_car == cCartoon_auto) {
                cur_car = cCartoon_loop;
              }
              parity = 0;
              *ss = 0;
              break;
            }

            *(cc++) = cur_car;

            if(cur_car == cCartoon_putty)
              putty_flag = true;

            v1 = cs->Coord + 3 * a;
            *(v++) = *(v1++);
            *(v++) = *(v1++);
            *(v++) = *(v1++);
            a2 = a1;

            ss++;

            fp++;

            v_c = NULL;
            v_n = NULL;
            v_o = NULL;

            AtomInfoBracketResidueFast(G, obj->AtomInfo, obj->NAtom, a1, &st, &nd);

            if(obj->DiscreteFlag) {
              if(cs == obj->DiscreteCSet[nd])
                skip_to = obj->DiscreteAtmToIdx[nd];
            } else
              skip_to = cs->AtmToIdx[nd];

            for(a3 = st; a3 <= nd; a3++) {

              if(obj->DiscreteFlag) {
                if(cs == obj->DiscreteCSet[a3])
                  a4 = obj->DiscreteAtmToIdx[a3];
                else
                  a4 = -1;
              } else
                a4 = cs->AtmToIdx[a3];
              if(a4 >= 0) {
                if(WordMatch(G, "C", obj->AtomInfo[a3].name, 1) < 0) {
                  v_c = cs->Coord + 3 * a4;
                } else if(WordMatch(G, "N", obj->AtomInfo[a3].name, 1) < 0) {
                  v_n = cs->Coord + 3 * a4;
                } else if(WordMatch(G, "O", obj->AtomInfo[a3].name, 1) < 0) {
                  v_o = cs->Coord + 3 * a4;
                }
              }
            }
            if(!(v_c && v_n && v_o)) {
              vo[0] = 0.0;
              vo[1] = 0.0;
              vo[2] = 0.0;
              vo += 3;
            } else {
              /* generate orientation vectors... */

              subtract3f(v_n, v_c, t0); /* t0 = N<---C */
              normalize3f(t0);
              subtract3f(v_n, v_o, t1); /* t1 = N<---O */
              normalize3f(t1);
              cross_product3f(t0, t1, vo);
              normalize3f(vo);
              if(parity) {
                invert3f(vo);
              }
              vo += 3;
            }
          } else if((((na_mode != 1) && (ai->protons == cAN_P) &&
                      (WordMatch(G, "P", ai->name, 1) < 0)) ||
                     ((na_mode == 1) && (ai->protons == cAN_C) &&
                      (WordMatchExact(G, NUCLEIC_NORMAL1, ai->name, 1) ||
                       WordMatchExact(G, NUCLEIC_NORMAL2, ai->name, 1))))
                    && !AtomInfoSameResidueP(G, last_ai, ai)) {
            if(a2 >= 0) {
              if(!ObjectMoleculeCheckBondSep(obj, a1, a2, 6)) { /* six bonds between phosphates */

                /*  3' cap */

                if(trailing_O3p_ai && ((na_mode == 2) || (na_mode == 4))) {

                  nuc_acid(G, trailing_O3p_a, trailing_O3p_a1, trailing_O3p_ai,
                           cs, obj, na_mode, nuc_flag, false, &a2, &nSeg, &v_o_last, &s,
                           &i, &cc, &nAt, &cur_car, &ss, &putty_flag, &v, &vo);
                }
                a2 = -1;
              }
            }
            last_ai = ai;
            trailing_O3p_ai = NULL;

            /*  5' cap */

            if(leading_O5p_ai && (a2 < 0) && ((na_mode == 3) || (na_mode == 4))) {
              if((!AtomInfoSameResidueP(G, ai, leading_O5p_ai)) &&
                 ObjectMoleculeCheckBondSep(obj, a1, leading_O5p_a1, 5)) {

                nuc_acid(G, leading_O5p_a, leading_O5p_a1, leading_O5p_ai,
                         cs, obj, na_mode, nuc_flag, false, &a2, &nSeg, &v_o_last, &s, &i,
                         &cc, &nAt, &cur_car, &ss, &putty_flag, &v, &vo);
              }
            }

            leading_O5p_ai = NULL;

            /* this is the main nucleic acid cartoon section... */

            nuc_acid(G, a, a1, ai, cs, obj, na_mode, nuc_flag, true, &a2, &nSeg,
                     &v_o_last, &s, &i, &cc, &nAt, &cur_car, &ss, &putty_flag, &v, &vo);

          } else if((a2 >= 0) &&
                    last_ai &&
                    (ai->protons == cAN_O) &&
                    (last_ai->protons == cAN_P) &&
                    ((na_mode == 2) || (na_mode == 4)) &&
                    (WordMatchExact(G, "O3'", ai->name, 1) ||
                     WordMatchExact(G, "O3*", ai->name, 1)) &&
                    AtomInfoSameResidueP(G, last_ai, ai) &&
                    ObjectMoleculeCheckBondSep(obj, a1, a2, 5)) {
            trailing_O3p_ai = ai;
            trailing_O3p_a = a;
            trailing_O3p_a1 = a1;
          } else if((ai->protons == cAN_O) &&
                    ((na_mode == 3) || (na_mode == 4)) &&
                    (WordMatchExact(G, "O5'", ai->name, 1) ||
                     WordMatchExact(G, "O5*", ai->name, 1))) {
            leading_O5p_ai = ai;
            leading_O5p_a = a;
            leading_O5p_a1 = a1;
          }

        }
      }
    }
  }

  /* BEGIN 3' cap */
  if(trailing_O3p_ai && ((na_mode == 2) || (na_mode == 4))) {

    nuc_acid(G, trailing_O3p_a, trailing_O3p_a1, trailing_O3p_ai,
             cs, obj, na_mode, nuc_flag, false, &a2, &nSeg, &v_o_last, &s, &i, &cc,
             &nAt, &cur_car, &ss, &putty_flag, &v, &vo);
    a2 = -1;
    trailing_O3p_ai = NULL;

  }

  if(nAt && putty_flag) {
    double sum = 0.0, sumsq = 0.0;
    float value;
    int cnt = 0;

    for(a = 0; a < obj->NAtom; a++) {
      ai = obj->AtomInfo + a;

      if(ai->visRep[cRepCartoon]) {
        value = ai->b;
        sum += value;
        sumsq += (value * value);
        if(value < putty_min)
          putty_min = value;
        if(value > putty_max)
          putty_max = value;
        cnt++;
      }
    }

    if(cnt) {
      putty_mean = (float) (sum / cnt);
      putty_stdev = (float) sqrt1d((sumsq - (sum * sum / cnt)) / (cnt));
    } else {
      /* aren't these assignments unnecessary? */
      putty_mean = 10.0F;
      putty_stdev = 10.0F;
      putty_min = 0.0F;
      putty_max = 10.0F;
    }
  }

  PRINTFD(G, FB_RepCartoon)
    " RepCartoon-Debug: path outlined, interpolating...\n" ENDFD;

  if(nAt) {

    /* compute differences and normals */

    s = seg;
    v = pv;

    dv = Alloc(float, nAt * 3);
    nv = Alloc(float, nAt * 3);
    dl = Alloc(float, nAt);
    v1 = dv;
    v2 = nv;
    d = dl;
    for(a = 0; a < (nAt - 1); a++) {
      PRINTFD(G, FB_RepCartoon)
        " RepCartoon: seg %d *s %d , *(s+1) %d\n", a, *s, *(s + 1)
        ENDFD;

      if(*s == *(s + 1)) {
        subtract3f(v + 3, v, v1);
        *d = (float) length3f(v1);
        if(*d > R_SMALL4) {
          float d_1;
          d_1 = 1.0F / (*d);
          scale3f(v1, d_1, v2);
        } else if(a) {
          copy3f(v2 - 3, v2);
        } else {
          zero3f(v2);
        }
      } else {
        zero3f(v2);
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

    tv = Alloc(float, nAt * 3 + 6);
    v1 = tv;

    *(v1++) = *(v++);           /* first segment */
    *(v1++) = *(v++);
    *(v1++) = *(v++);
    s++;

    for(a = 1; a < (nAt - 1); a++) {
      if((*s == *(s - 1)) && (*s == *(s + 1))) {
        add3f(v, (v - 3), v1);  /* tangent vectors are head-to-tail sums within a segment */
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

    PRINTFD(G, FB_RepCartoon)
      " RepCartoon-Debug: generating coordinate systems...\n" ENDFD;

    if(round_helices) {
      v1 = NULL;
      v2 = NULL;
      v3 = NULL;
      v4 = NULL;
      v5 = NULL;
      s = seg;
      v = pv;
      ss = sstype;
      vo = pvo;
      v0 = tv;
      last = 0;
      if(nAt > 1) {
        for(a = 0; a < nAt; a++) {
          if(a) {
            if(*s != *(s - 1)) {        /* contiguous helices in disconnected segments */
              v1 = NULL;
              v2 = NULL;
              v3 = NULL;
              v4 = NULL;
              v5 = NULL;
              last = 0;
            }
          }
          v5 = v4;
          v4 = v3;
          v3 = v2;
          v2 = v1;
          if(*ss == 1)          /* helix */
            v1 = v;
          else {                /* early termination ? */
            if(last < 2) {
              zero3f(t0);
              if(v2 && v3) {
                subtract3f(v2, v, t0);
                normalize3f(t0);
                subtract3f(v3, v2, t1);
                normalize3f(t1);
                add3f(t1, t0, t0);
                if(v4) {
                  subtract3f(v4, v3, t1);
                  normalize3f(t1);
                  add3f(t1, t0, t0);
                }
                if(v5) {
                  subtract3f(v5, v4, t1);
                  normalize3f(t1);
                  add3f(t1, t0, t0);
                }
                normalize3f(t0);
                cross_product3f(t0, v0 - 3, vo - 3);
                normalize3f(vo - 3);
                cross_product3f(t0, v0 - 6, vo - 6);
                normalize3f(vo - 6);
                if(v4) {
                  cross_product3f(t0, v0 - 9, vo - 9);
                  normalize3f(vo - 9);
                }
                if(v5) {
                  cross_product3f(t0, v0 - 12, vo - 12);
                  normalize3f(vo - 12);
                }

                if(v4 && v5) {
                  /* now make sure there's no goofy flip on the end...
                     of a short, tight helix */
                  if(dot_product3f(vo - 9, vo - 12) < -0.8F)
                    invert3f(vo - 12);
                }
              }
            }
            v1 = NULL;
            v2 = NULL;
            v3 = NULL;
            v4 = NULL;
            v5 = NULL;
            last = 0;
          }
          if(v1 && v2 && v3 && v4) {
            add3f(v1, v4, t0);
            add3f(v2, v3, t1);
            scale3f(t0, 0.2130F, t0);
            scale3f(t1, 0.2870F, t1);

            add3f(t0, t1, t0);
            if(last) {          /* 5th CA or later... */
              subtract3f(t2, t0, t1);
              normalize3f(t1);
              cross_product3f(t1, v0, vo);
              normalize3f(vo);
              cross_product3f(t1, v0 - 3, vo - 3);
              normalize3f(vo - 3);
              cross_product3f(t1, v0 - 6, vo - 6);
              normalize3f(vo - 6);
              if(last == 1) {   /* 5th */
                cross_product3f(t1, v0 - 9, vo - 9);
                normalize3f(vo - 9);
                cross_product3f(t1, v0 - 12, vo - 12);
                normalize3f(vo - 12);
              }
            }
            last++;
            copy3f(t0, t2);
          }
          v += 3;
          ss++;
          vo += 3;
          v0 += 3;
          s++;
        }
      }
    }

    {
      int refine_normals =
        SettingGet_i(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_refine_normals);
      if(refine_normals < 0) {
        if(obj->NCSet > 1) {
          int i, n_set = 0;
          for(i = 0; i < obj->NCSet; i++)
            if(obj->CSet[i]) {
              n_set++;
              if(n_set > 1)
                refine_normals = 0;
              /* default behavior is to not refine normals for multi-state objects */
            }
        }
      }

      if(refine_normals) {

        /* first, make sure orientiation vectors are orthogonal to the tangent */

        v1 = tv + 3;
        vo = pvo + 3;
        s = seg + 1;
        for(a = 1; a < (nAt - 1); a++) {
          if((*s == *(s - 1)) && (*s == *(s + 1))) {
            /* only operate on vectors within the cartoon itself --
               not the end vectors */

            remove_component3f(vo, v1, t0);
            normalize23f(t0, vo);

            /* go on to next vertex */
          }
          v1 += 3;
          vo += 3;
          s++;
        }

        /* now generate alternative inverted orientation vectors */

        va = pva;
        vo = pvo;
        ss = sstype;
        for(a = 0; a < nAt; a++) {

          /* original */
          copy3f(vo, va);
          va += 3;

          /* inverse */
          copy3f(vo, va);
          if(*ss != 1) {
            invert3f(va);
            /* for helix, don't allow inversion of normals, since that
               would confuse the inside & outside of the helix  */
          }
          va += 3;

          /* go on to next vertex */

          vo += 3;
          ss++;
        }

        /* now iterate forward through pairs */

        vo = pvo + 3;
        va = pva + 6;
        v = nv + 3;             /* normals in direction of chain */
        s = seg + 1;

        for(a = 1; a < (nAt - 1); a++) {

          if((*s == *(s + 1)) && (*s == *(s - 1))) {    /* only operate within a segment */
            remove_component3f(vo - 3, v - 3, o0);      /* previous orientation vector */
            normalize3f(o0);    /* is now perp to chain direction */

            v1 = va;            /* candidate orientation vectors for current CA */

            remove_component3f(v1, v - 3, o1);  /* removes chain direction from the two candidates */
            remove_component3f(v1 + 3, v - 3, o1 + 3);
            normalize3f(o1);
            normalize3f(o1 + 3);

            max_dot = dot_product3f(o0, o1);
            v0 = v1;

            dp = dot_product3f(o0, o1 + 3);
            if(dp > max_dot) {
              v0 = v1 + 3;
              max_dot = dp;
            }
            copy3f(v0, vo);     /* updates atom with optimal orientation vector */
          }
          vo += 3;
          va += 6;              /* candidate orientation vectors */
          v += 3;               /* normal */
          s++;
        }

        /* now soften up the kinks */
        v1 = tv + 3;
        va = pva + 6;
        vo = pvo + 3;
        s = seg + 1;
        for(a = 1; a < (nAt - 1); a++) {
          if((*s == *(s - 1)) && (*s == *(s + 1))) {
            dp = (dot_product3f(vo, vo + 3) * dot_product3f(vo, vo - 3));
            if(dp < -0.10F) {   /* threshold value -- could be a setting */
              add3f(vo + 3, vo - 3, t0);
              scale3f(vo, 0.001, t1);
              add3f(t1, t0, t0);
              remove_component3f(t0, v1, t0);
              normalize3f(t0);
              if(dot_product3f(vo, t0) < 0.0F) {
                subtract3f(vo, t0, t2);
              } else {
                add3f(vo, t0, t2);
              }
              normalize3f(t2);
              dp = 2 * (-0.10F - dp);
              if(dp > 1.0F)
                dp = 1.0F;
              mix3f(vo, t2, dp, t3);
              copy3f(t3, va);   /* store modified vector */
              invert3f3f(va, va + 3);
            } else {
              copy3f(vo, va);   /* keep as is */
            }
          }
          v1 += 3;
          vo += 3;
          va += 6;
          s++;
        }

        /* now update */
        va = pva + 6;
        vo = pvo + 3;
        s = seg + 1;
        for(a = 1; a < (nAt - 1); a++) {
          if((*s == *(s - 1)) && (*s == *(s + 1))) {
            copy3f(va, vo);
          }
          vo += 3;
          va += 6;
          s++;
        }
      }
    }

    if(SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_flat_sheets)) {
      s = seg;
      cc = car;
      v = pv;
      ss = sstype;
      vo = pvo;
      v0 = tv;
      last = 0;
      first = -1;
      cur_car = *cc;
      end_flag = false;
      if(nAt > 1) {
        for(a = 0; a < nAt; a++) {
          if(a) {
            if(*s != *(s - 1)) {
              end_flag = true;
            } else if(*ss != 2) {
              end_flag = true;
            }
            if(a == (nAt - 1)) {
              end_flag = 1;
            }
          }
          if(end_flag) {
            if((cur_car != cCartoon_loop) && (cur_car != cCartoon_tube)) {

              f = 1;
              for(c = 0; c < flat_cycles; c++) {
                for(b = first + f; b <= last - f; b++) {        /* iterative averaging */
                  zero3f(t0);
                  for(e = -f; e <= f; e++) {
                    add3f(pv + 3 * (b + e), t0, t0);
                  }
                  scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
                }
                for(b = first + f; b <= last - f; b++) {
                  if(!((*(flag_tmp + b) & cAtomFlag_no_smooth))) {
                    copy3f(tmp + b * 3, pv + b * 3);
                  }
                }
                for(b = first + f; b <= last - f; b++) {
                  zero3f(t0);
                  for(e = -f; e <= f; e++) {
                    add3f(pvo + 3 * (b + e), t0, t0);
                  }
                  scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
                }
                for(b = first + f; b <= last - f; b++) {
                  copy3f(tmp + b * 3, pvo + b * 3);
                  /*                  normalize3f(pvo+b*3); */
                }
                for(b = first + f; b <= last - f; b++) {
                  subtract3f(pv + (b + 1) * 3, pv + (b - 1) * 3, tmp + b * 3);
                  normalize3f(tmp + b * 3);
                  remove_component3f(pvo + b * 3, tmp + b * 3, pvo + b * 3);
                  normalize3f(pvo + b * 3);
                }
              }
            }
            first = -1;
            last = -1;
            end_flag = false;
          }
          if(*ss == 2) {
            if(first < 0)
              first = a;
            cur_car = *cc;
            last = a;
          }
          v += 3;
          ss++;
          vo += 3;
          v0 += 3;
          s++;
          cc++;
        }
      }
    }

    if(smooth_loops) {

      s = seg;
      v = pv;
      ss = sstype;
      vo = pvo;
      v0 = tv;
      last = 0;
      first = -1;
      end_flag = false;
      if(nAt > 1) {
        for(a = 0; a < nAt; a++) {

          if(a) {
            if(*s != *(s - 1)) {
              end_flag = true;
            } else if(*ss != 0) {
              end_flag = true;
            }
            if(a == (nAt - 1))
              end_flag = 1;
          }
          if(end_flag) {

            if(a)
              if(first > 0)     /* 011130 WLD */
                if(*(seg + first) == *(seg + first - 1))
                  first--;

            if(last > 0)
              if(*s == *(s - 1))
                if(last < (nAt - 1))
                  last++;

            for(f = smooth_first; f <= smooth_last; f++) {
              for(c = 0; c < smooth_cycles; c++) {
                for(b = first + f; b <= last - f; b++) {        /* iterative averaging */
                  zero3f(t0);
                  for(e = -f; e <= f; e++) {
                    add3f(pv + 3 * (b + e), t0, t0);
                  }
                  scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
                }
                for(b = first + f; b <= last - f; b++) {
                  if(!(*(flag_tmp + b) & cAtomFlag_no_smooth)) {
                    copy3f(tmp + b * 3, pv + b * 3);
                  }
                }
                for(b = first + f; b <= last - f; b++) {
                  zero3f(t0);
                  for(e = -f; e <= f; e++) {
                    add3f(pvo + 3 * (b + e), t0, t0);
                  }
                  scale3f(t0, 1.0F / (f * 2 + 1), tmp + b * 3);
                }
                for(b = first + f; b <= last - f; b++) {
                  copy3f(tmp + b * 3, pvo + b * 3);
                  normalize3f(pvo + b * 3);
                }
              }

            }
            first = -1;
            last = -1;
            end_flag = false;
          }
          if(*ss == 0) {
            if(first < 0)
              first = a;
            last = a;
          }
          v += 3;
          ss++;
          vo += 3;
          v0 += 3;
          s++;
        }
      }
    }

    if(smooth_loops ||
       SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_flat_sheets)) {

      /* recompute differences and normals */

      s = seg;
      v = pv;
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
        } else {
          zero3f(v2);
        }
        d++;
        v += 3;
        v1 += 3;
        v2 += 3;
        s++;
      }

      /* recompute tangents */

      s = seg;
      v = nv;

      v1 = tv;

      *(v1++) = *(v++);         /* first segment */
      *(v1++) = *(v++);
      *(v1++) = *(v++);
      s++;

      for(a = 1; a < (nAt - 1); a++) {
        if((*s == *(s - 1)) && (*s == *(s + 1))) {
          add3f(v, (v - 3), v1);

          normalize3f(v1);
        } else if(*s == *(s - 1)) {
          *(v1) = *(v - 3);     /* end a segment */
          *(v1 + 1) = *(v - 2);
          *(v1 + 2) = *(v - 1);
        } else if(*s == *(s + 1)) {
          *(v1) = *(v);         /* new segment */
          *(v1 + 1) = *(v + 1);
          *(v1 + 2) = *(v + 2);
        }
        v += 3;
        v1 += 3;
        s++;
      }

      *(v1++) = *(v - 3);       /* last segment */
      *(v1++) = *(v - 2);
      *(v1++) = *(v - 1);

      if(SettingGet_f(G, cs->Setting, obj->Obj.Setting, cSetting_cartoon_flat_sheets)) {
        s = seg + 1;
        ss = sstype + 1;
        v2 = tv + 3;            /* normal */
        for(a = 1; a < (nAt - 1); a++) {
          if((*ss == 2) && (*s == *(s + 1)) && (*s == *(s - 1))) {      /* sheet in same segment */
            if((*ss == *(ss + 1)) && (*ss != *(ss - 1))) {      /* start, bias forwards */
              scale3f(v2 + 3, refine_tips, t0);
              add3f(t0, v2, v2);
              normalize3f(v2);
            } else if((*ss != *(ss + 1)) && (*ss == *(ss - 1))) {       /* end, bias backwards */
              scale3f(v2 - 3, refine_tips, t0);
              add3f(t0, v2, v2);
              normalize3f(v2);
            }
          }
          v2 += 3;
          s++;
          ss++;
        }

      }

    }
  }

  I->ray = CGONew(G);

  if(alpha != 1.0F)
    CGOAlpha(I->ray, alpha);

  /* debugging output */
  if(round_helices) {
    if((cartoon_debug > 0.5) && (cartoon_debug < 2.5)) {
      CGOColor(I->ray, 1.0, 1.0, 1.0);
      CGODisable(I->ray, GL_LIGHTING);
      CGOBegin(I->ray, GL_LINE_STRIP);

      v1 = NULL;
      v2 = NULL;
      v3 = NULL;
      v4 = NULL;
      v = pv;
      if(nAt > 1) {
        CGOBegin(I->ray, GL_LINE_STRIP);
        for(a = 0; a < nAt; a++) {
          v4 = v3;
          v3 = v2;
          v2 = v1;
          v1 = v;
          if(v1 && v2 && v3 && v4) {
            add3f(v1, v4, t0);
            add3f(v2, v3, t1);
            /*            scale3f(t0,0.2024,t0);
               scale3f(t1,0.2727,t1); */

            scale3f(t0, 0.2130F, t0);
            scale3f(t1, 0.2870F, t1);

            add3f(t0, t1, t0);
            CGOVertexv(I->ray, t0);
          }
          v += 3;
        }
        CGOEnd(I->ray);
      }
    }
  }

  PRINTFD(G, FB_RepCartoon)
    " RepCartoon-Debug: creating 3D scaffold...\n" ENDFD;

  /* okay, we now have enough info to generate smooth interpolations */

  if(nAt > 1) {
    ex = ExtrudeNew(G);
    ExtrudeAllocPointsNormalsColors(ex, cs->NIndex * (3 * sampling + 3));
  }

  /* process cylindrical helices first */

  if((nAt > 1) && cylindrical_helices) {

    /* this is confusing because we're borrowing Extrude's arrays 
     * for convenient storage, but not actually calling Extrude */

    n_p = 0;
    v = ex->p;
    vc = ex->c;
    vn = ex->n;
    vi = ex->i;

    last_color = -1;
    uniform_color = true;

    v1 = pv;                    /* points */
    v2 = tv;                    /* tangents */
    vo = pvo;
    d = dl;
    s = seg;
    cc = car;
    atp = at;                   /* cs index pointer */
    a = 0;
    contFlag = true;
    cur_car = cCartoon_skip;
    extrudeFlag = false;
    contigFlag = false;

    while(contFlag) {
      if((*cc) != cur_car) {    /* new cartoon type */
        if(n_p) {               /* any cartoon points? */
          extrudeFlag = true;
        } else {
          cur_car = *(cc);      /* now: go ahead and switch cartoons */
          n_p = 0;
          v = ex->p;
          vc = ex->c;
          vi = ex->i;
          vn = ex->n;
          last_color = -1;
          uniform_color = true;
        }
      }
      if(a && !extrudeFlag) {
        if((*s) != *(s - 1)) {  /* new segment */
          contigFlag = false;
          if(n_p) {             /* any cartoon points? */
            extrudeFlag = true;
          } else {
            n_p = 0;
            v = ex->p;
            vc = ex->c;
            vi = ex->i;
            vn = ex->n;
            last_color = -1;
            uniform_color = true;
          }
        }
      }
      if(!extrudeFlag) {
        if((a < (nAt - 1)) && (*s == *(s + 1))) {       /* working in the same segment... */
          atom_index1 = cs->IdxToAtm[*atp];
          atom_index2 = cs->IdxToAtm[*(atp + 1)];
          c1 = *(cs->Color + *atp);
          c2 = *(cs->Color + *(atp + 1));

          if(cartoon_color >= 0) {
            c1 = (c2 = cartoon_color);
          }

          AtomInfoGetSetting_color(G, obj->AtomInfo + atom_index1, cSetting_cartoon_color,
                                   c1, &c1);
          AtomInfoGetSetting_color(G, obj->AtomInfo + atom_index2, cSetting_cartoon_color,
                                   c2, &c2);

          if(discrete_colors) {
            if(n_p == 0) {
              if(contigFlag) {
                if(cur_car != cCartoon_loop)
                  c2 = c1;
                else {
                  if((*cc + 1) == cur_car)
                    c2 = c1;
                  else
                    c1 = c2;
                }
              } else if((cur_car == cCartoon_loop) && (*(cc + 1) != cCartoon_loop)) {
                c2 = c1;
              }
            } else {
              if((cur_car == cCartoon_loop) && (*(cc + 1) != cCartoon_loop)) {
                c2 = c1;
              }
            }                   /* not contig */

          }

          if((*(cc) == *(cc + 1)) && (c1 != c2))
            uniform_color = false;
          if(last_color >= 0) {
            if(c1 != last_color)
              uniform_color = false;
          }
          last_color = c1;

          v0 = ColorGet(G, c1);
          *(vc++) = *(v0++);
          *(vc++) = *(v0++);
          *(vc++) = *(v0++);
          *(vi++) = atom_index1;

          v0 = ColorGet(G, c2); /* kludge */
          *(vc) = *(v0++);
          *(vc + 1) = *(v0++);
          *(vc + 2) = *(v0++);
          *(vi) = atom_index2;
        } else {
          vc += 3;              /* part of kludge */
          vi++;
        }
        if(cur_car == cCartoon_skip_helix) {
          if(!n_p) {
            h_start = v1;
            h_end = v1;
          } else {
            h_end = v1;
          }
          copy3f(v1, v);        /* just store coordinates until we have a complete cylinder */
          v += 3;
          n_p++;
        }
        v1 += 3;
        v2 += 3;
        v3 += 3;
        vo += 3;
        d++;
        atp++;
        s++;
        cc++;
      }

      a++;
      if(a == nAt) {
        contFlag = false;
        if(n_p)
          extrudeFlag = true;
      }

      if(extrudeFlag) {         /* generate cylinder */
        contigFlag = true;
        if((a < nAt) && extrudeFlag) {
          if(*(s - 1) != *(s))
            contigFlag = false;
        }

        if(n_p > 1) {
          atom_index1 = cs->IdxToAtm[*(atp - 1)];
          c1 = *(cs->Color + *(atp - 1));

          if(cartoon_color >= 0) {
            c1 = cartoon_color;
          }

          AtomInfoGetSetting_color(G, obj->AtomInfo + atom_index1, cSetting_cartoon_color,
                                   c1, &c1);

          if(n_p < 5) {
            copy3f(ex->p, t3);
            copy3f(v - 3, t4);
          } else {
            add3f(ex->p, ex->p + 9, t0);
            add3f(ex->p + 3, ex->p + 6, t1);
            scale3f(t0, 0.2130F, t0);
            scale3f(t1, 0.2870F, t1);
            add3f(t0, t1, t3);

            add3f(v - 3, v - 12, t0);
            add3f(v - 6, v - 9, t1);
            scale3f(t0, 0.2130F, t0);
            scale3f(t1, 0.2870F, t1);
            add3f(t0, t1, t4);

            /* extend helix to line up with CA */
            subtract3f(t4, t3, t0);
            normalize3f(t0);
            subtract3f(v - 3, t3, t1);
            project3f(t1, t0, t4);
            add3f(t3, t4, t4);
            invert3f(t0);
            subtract3f(ex->p, t4, t1);
            project3f(t1, t0, t3);
            add3f(t3, t4, t3);

            /* relocate terminal CA to touch helix, if necessary */

            if(h_start && h_end) {
              subtract3f(h_start, t3, t0);
              f0 = helix_radius - loop_radius * 2;
              if(length3f(t0) > f0) {
                normalize3f(t0);
                scale3f(t0, f0, t1);
                add3f(t1, t3, h_start);
              }

              subtract3f(h_end, t4, t0);
              if(length3f(t0) > f0) {
                normalize3f(t0);
                scale3f(t0, f0, t1);
                add3f(t1, t4, h_end);
              }
            }
          }

          /* push helix out a tad to consume loop */

          subtract3f(t4, t3, t0);
          normalize3f(t0);
          scale3f(t0, loop_radius * 2, t0);
          add3f(t0, t4, t4);
          invert3f(t0);
          add3f(t0, t3, t3);

          if(uniform_color) {
            CGOCylinderv(I->ray, t3, t4, helix_radius, ex->c, ex->c);

          } else {
            subtract3f(t4, t3, t0);
            n_pm1 = n_p - 1;
            n_pm2 = n_p - 2;
            for(b = 0; b < n_pm1; b++) {
              if(!b) {
                scale3f(t0, ((float) b - 0.005F) / n_pm1, t1);  /* add small overlap */
              } else {
                scale3f(t0, ((float) b) / n_pm1, t1);
              }
              if(b < n_pm2) {
                scale3f(t0, ((float) b + 1.005F) / n_pm1, t2);
              } else {
                scale3f(t0, ((float) b + 1) / n_pm1, t2);
              }
              add3f(t3, t1, t1);
              add3f(t3, t2, t2);
              CGOCustomCylinderv(I->ray, t1, t2, helix_radius, ex->c + (b * 3),
                                 ex->c + (b + 1) * 3, (float) (b ? 0 : cCylCapFlat),
                                 (float) (b == n_pm2 ? cCylCapFlat : 0));
            }
          }
        }
        a--;                    /* undo above... */
        extrudeFlag = false;
        n_p = 0;
        v = ex->p;
        vc = ex->c;
        vi = ex->i;
        vn = ex->n;
        uniform_color = true;
        last_color = -1;
      }
    }
  }

  if(nAt > 1) {
    n_p = 0;
    v = ex->p;
    vc = ex->c;
    vn = ex->n;
    vi = ex->i;

    v1 = pv;                    /* points */
    v2 = tv;                    /* tangents */
    vo = pvo;
    d = dl;
    s = seg;
    cc = car;
    atp = at;                   /* cs index pointer */
    a = 0;
    contFlag = true;
    cur_car = cCartoon_skip;
    extrudeFlag = false;
    contigFlag = false;

    while(contFlag) {

      if((*cc) != cur_car) {    /* new cartoon type */
        if(n_p) {               /* any cartoon points? */
          extrudeFlag = true;
        } else {
          cur_car = *(cc);      /* no: go ahead and switch cartoons */
          ExtrudeTruncate(ex, 0);
          n_p = 0;
          v = ex->p;
          vc = ex->c;
          vn = ex->n;
          vi = ex->i;
        }
      }

      /* CONFUSION ALERT -- I don't understand the following code (anymore) */

      if(a < (nAt - 1)) {
        /* put a setting controlled conditional here.. */
        if(((*(cc + 1)) != cur_car) && (cur_car != cCartoon_loop)) {    /* end of segment */
          if(n_p) {             /* any cartoon points? */
            extrudeFlag = true;
          } else {
            cur_car = cCartoon_loop;    /* no: go ahead and switch cartoons */
            ExtrudeTruncate(ex, 0);
            n_p = 0;
            v = ex->p;
            vc = ex->c;
            vn = ex->n;
            vi = ex->i;
          }
        }
      }
      if((a < (nAt - 1)) && !extrudeFlag) {
        if((*s) != *(s + 1)) {  /* new segment */
          contigFlag = false;
          if(n_p) {             /* any cartoon points? */
            extrudeFlag = true;
          } else {
            ExtrudeTruncate(ex, 0);
            n_p = 0;
            v = ex->p;
            vc = ex->c;
            vn = ex->n;
            vi = ex->i;
          }
        }
      }
      if(!extrudeFlag) {
        if((a < (nAt - 1)) && (*s == *(s + 1))) {       /* working in the same segment... */
          c1 = *(cs->Color + *atp);
          c2 = *(cs->Color + *(atp + 1));
          atom_index1 = cs->IdxToAtm[*atp];
          atom_index2 = cs->IdxToAtm[*(atp + 1)];

          if(cartoon_color >= 0) {
            c1 = (c2 = cartoon_color);
          }

          AtomInfoGetSetting_color(G, obj->AtomInfo + atom_index1, cSetting_cartoon_color,
                                   c1, &c1);
          AtomInfoGetSetting_color(G, obj->AtomInfo + atom_index2, cSetting_cartoon_color,
                                   c2, &c2);

          if(nuc_flag[*atp] || nuc_flag[*(atp + 1)]) {  /* this is a nucleic acid ribbon */
            if(nucleic_color >= 0) {
              c1 = (c2 = nucleic_color);
            }
          }

          if(discrete_colors) {
            if(n_p == 0) {
              if(contigFlag) {
                if(cur_car != cCartoon_loop)
                  c2 = c1;
                else {
                  if((*cc + 1) == cur_car)
                    c2 = c1;
                  else
                    c1 = c2;
                }
              } else if((cur_car == cCartoon_loop) && (*(cc + 1) != cCartoon_loop)) {
                c2 = c1;
              }
            } else {
              if((cur_car == cCartoon_loop) && (*(cc + 1) != cCartoon_loop)) {
                c2 = c1;
              }
            }                   /* not contig */
          }

          dev = throw * (*d);
          for(b = 0; b < sampling; b++) {       /* needs optimization */

            if(n_p == 0) {

              /* provide starting point on first point in segment only... */

              f0 = ((float) b) / sampling;      /* fraction of completion */
              if(f0 <= 0.5) {
                v0 = ColorGet(G, c1);
                i0 = atom_index1;
              } else {
                v0 = ColorGet(G, c2);
                i0 = atom_index2;
              }
              f0 = smooth(f0, power_a); /* bias sampling towards the center of the curve */

              /* store colors */

              *(vc++) = *(v0++);
              *(vc++) = *(v0++);
              *(vc++) = *(v0++);
              *(vi++) = i0;
              /* start of line/cylinder */

              f1 = 1.0F - f0;
              f2 = smooth(f0, power_b);
              f3 = smooth(f1, power_b);
              f4 = dev * f2 * f3;       /* displacement magnitude */

              *(v++) = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);

              *(v++) = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);

              *(v++) = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);

              vn += 9;

              copy3f(vo, vn - 6);       /* starter... */

              n_p++;

            }

            f0 = ((float) b + 1) / sampling;
            if(f0 <= 0.5) {
              v0 = ColorGet(G, c1);
              i0 = atom_index1;
            } else {
              v0 = ColorGet(G, c2);
              i0 = atom_index2;
            }
            f0 = smooth(f0, power_a);   /* bias sampling towards the center of the curve */

            /* store colors */

            *(vc++) = *(v0++);
            *(vc++) = *(v0++);
            *(vc++) = *(v0++);
            *(vi++) = i0;

            /* end of line/cylinder */

            f1 = 1.0F - f0;
            f2 = smooth(f0, power_b);
            f3 = smooth(f1, power_b);
            f4 = dev * f2 * f3; /* displacement magnitude */

            *(v++) = f1 * v1[0] + f0 * v1[3] + f4 * (f3 * v2[0] - f2 * v2[3]);

            *(v++) = f1 * v1[1] + f0 * v1[4] + f4 * (f3 * v2[1] - f2 * v2[4]);

            *(v++) = f1 * v1[2] + f0 * v1[5] + f4 * (f3 * v2[2] - f2 * v2[5]);

            /*                remove_component3f(vo,v2,o0);
               remove_component3f(vo+3,v2,o0+3); */

            vn += 3;
            *(vn++) = f1 * (vo[0] * f2) + f0 * (vo[3] * f3);
            *(vn++) = f1 * (vo[1] * f2) + f0 * (vo[4] * f3);
            *(vn++) = f1 * (vo[2] * f2) + f0 * (vo[5] * f3);
            vn += 3;

            if(b == sampling - 1)
              copy3f(vo + 3, vn - 6);   /* starter... */
            n_p++;

          }

          /* now do a smoothing pass along orientation 
             vector to smooth helices, etc... */

          c = refine;
          cross_product3f(vn + 3 - (sampling * 9), vn + 3 - 9, t0);

          cross_product3f(vo, vo + 3, t0);
          if((sampling > 1) && length3f(t0) > R_SMALL4) {

            normalize3f(t0);
            while(c--) {
              p0 = v - (sampling * 3) - 3;
              p1 = v - (sampling * 3);
              p2 = v - (sampling * 3) + 3;
              for(b = 0; b < (sampling - 1); b++) {
                f0 = dot_product3f(t0, p0);
                f1 = dot_product3f(t0, p1);
                f2 = dot_product3f(t0, p2);

                f3 = (f2 + f0) / 2.0F;
                scale3f(t0, f3 - f1, t1);
                p3 = sampling_tmp + b * 3;
                add3f(t1, p1, p3);

                p0 = p1;
                p1 = p2;
                p2 += 3;
              }
              p1 = v - (sampling * 3);
              for(b = 0; b < (sampling - 1); b++) {
                p3 = sampling_tmp + b * 3;
                copy3f(p3, p1);
                p1 += 3;
              }
            }
          }
        }
        v1 += 3;
        v2 += 3;
        v3 += 3;
        vo += 3;
        d++;
        atp += 1;
        s++;
        cc++;

      }

      a++;
      if(a == nAt) {
        contFlag = false;
        if(n_p)
          extrudeFlag = true;
      }
      if(extrudeFlag) {
        contigFlag = true;
        if((a < nAt) && extrudeFlag) {
          if(*(s - 1) != *(s))
            contigFlag = false;
        }

        if((cur_car != cCartoon_skip) && (cur_car != cCartoon_skip_helix)) {

          if((cartoon_debug > 0.5) && (cartoon_debug < 2.5)) {
            CGOColor(I->ray, 0.0, 1.0, 0.0);

            v = ex->p;
            vn = ex->n + 3;
            CGODisable(I->ray, GL_LIGHTING);
            CGOBegin(I->ray, GL_LINES);
            for(b = 0; b < n_p; b++) {
              CGOVertexv(I->ray, v);
              add3f(v, vn, t0);
              CGOVertexv(I->ray, t0);
              v += 3;
              vn += 9;
            }
            CGOEnd(I->ray);
            CGOEnable(I->ray, GL_LIGHTING);
          }

          ExtrudeTruncate(ex, n_p);
          ExtrudeComputeTangents(ex);

          /* set up shape */
          switch (cur_car) {
          case cCartoon_tube:
            ExtrudeCircle(ex, tube_quality, tube_radius);
            ExtrudeBuildNormals1f(ex);
            ExtrudeCGOSurfaceTube(ex, I->ray, tube_cap, NULL);
            break;
          case cCartoon_putty:
            ExtrudeCircle(ex, putty_quality, putty_radius);
            ExtrudeBuildNormals1f(ex);
            ExtrudeComputePuttyScaleFactors(ex, obj,
                                            SettingGet_i(G, cs->Setting, obj->Obj.Setting,
                                                         cSetting_cartoon_putty_transform),
                                            putty_mean, putty_stdev, putty_min, putty_max,
                                            SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                                                         cSetting_cartoon_putty_scale_power),
                                            SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                                                         cSetting_cartoon_putty_range),
                                            SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                                                         cSetting_cartoon_putty_scale_min),
                                            SettingGet_f(G, cs->Setting, obj->Obj.Setting,
                                                         cSetting_cartoon_putty_scale_max),
                                            sampling / 2);

            ExtrudeCGOSurfaceVariableTube(ex, I->ray, 1);
            break;
          case cCartoon_loop:
            ExtrudeCircle(ex, loop_quality, loop_radius);
            ExtrudeBuildNormals1f(ex);
            ExtrudeCGOSurfaceTube(ex, I->ray, loop_cap, NULL);
            break;
          case cCartoon_rect:
            if(highlight_color < 0) {
              ExtrudeRectangle(ex, width, length, 0);
              ExtrudeBuildNormals2f(ex);
              ExtrudeCGOSurfacePolygon(ex, I->ray, 1, NULL);
            } else {
              ExtrudeRectangle(ex, width, length, 1);
              ExtrudeBuildNormals2f(ex);
              ExtrudeCGOSurfacePolygon(ex, I->ray, 0, NULL);
              ExtrudeRectangle(ex, width, length, 2);
              ExtrudeBuildNormals2f(ex);
              ExtrudeCGOSurfacePolygon(ex, I->ray, 1, ColorGet(G, highlight_color));
            }
            break;
          case cCartoon_oval:
            ExtrudeOval(ex, oval_quality, oval_width, oval_length);
            ExtrudeBuildNormals2f(ex);
            if(highlight_color < 0)
              ExtrudeCGOSurfaceTube(ex, I->ray, 1, NULL);
            else
              ExtrudeCGOSurfaceTube(ex, I->ray, 1, ColorGet(G, highlight_color));
            break;
          case cCartoon_arrow:
            ExtrudeRectangle(ex, width, length, 0);
            ExtrudeBuildNormals2f(ex);
            if(highlight_color < 0)
              ExtrudeCGOSurfaceStrand(ex, I->ray, sampling, NULL);
            else
              ExtrudeCGOSurfaceStrand(ex, I->ray, sampling, ColorGet(G, highlight_color));

            /* for PLY files      
               ExtrudeCircle(ex,loop_quality,loop_radius);
               ExtrudeBuildNormals1f(ex);
               ExtrudeCGOSurfaceTube(ex,I->ray,loop_cap,NULL);
             */
            break;
          case cCartoon_dumbbell:
            if(highlight_color < 0) {
              ExtrudeDumbbell1(ex, dumbbell_width, dumbbell_length, 0);
              ExtrudeBuildNormals2f(ex);
              ExtrudeCGOSurfacePolygonTaper(ex, I->ray, sampling, NULL);
            } else {

              ExtrudeDumbbell1(ex, dumbbell_width, dumbbell_length, 1);
              ExtrudeBuildNormals2f(ex);
              ExtrudeCGOSurfacePolygonTaper(ex, I->ray, sampling, NULL);

              ExtrudeDumbbell1(ex, dumbbell_width, dumbbell_length, 2);
              ExtrudeBuildNormals2f(ex);
              ExtrudeCGOSurfacePolygonTaper(ex, I->ray, sampling,
                                            ColorGet(G, highlight_color));
            }
            /*
               ExtrudeCGOSurfacePolygonX(ex,I->ray,1); */

            ex1 = ExtrudeCopyPointsNormalsColors(ex);
            ExtrudeDumbbellEdge(ex1, sampling, -1, dumbbell_length);
            ExtrudeComputeTangents(ex1);
            ExtrudeCircle(ex1, loop_quality, dumbbell_radius);
            ExtrudeBuildNormals1f(ex1);

            ExtrudeCGOSurfaceTube(ex1, I->ray, 1, NULL);
            ExtrudeFree(ex1);

            ex1 = ExtrudeCopyPointsNormalsColors(ex);
            ExtrudeDumbbellEdge(ex1, sampling, 1, dumbbell_length);
            ExtrudeComputeTangents(ex1);
            ExtrudeCircle(ex1, loop_quality, dumbbell_radius);
            ExtrudeBuildNormals1f(ex1);
            ExtrudeCGOSurfaceTube(ex1, I->ray, 1, NULL);
            ExtrudeFree(ex1);

            break;
          }
        }
        a--;                    /* undo above... */
        extrudeFlag = false;
        ExtrudeTruncate(ex, 0);
        n_p = 0;
        v = ex->p;
        vc = ex->c;
        vn = ex->n;
      }
    }
  }

  if(nAt > 1) {
    if((cartoon_debug > 0.5) && (cartoon_debug < 2.5)) {
      CGOColor(I->ray, 1.0, 1.0, 1.0);
      CGODisable(I->ray, GL_LIGHTING);
      CGOBegin(I->ray, GL_LINES);
      v1 = pv;
      v2 = pvo;
      v3 = tv;
      for(a = 0; a < nAt; a++) {
        CGOVertexv(I->ray, v1);
        add3f(v1, v2, t0);
        add3f(v2, t0, t0);
        CGOVertexv(I->ray, t0);
        subtract3f(v1, v3, t0);
        CGOVertexv(I->ray, t0);
        add3f(v1, v3, t0);
        CGOVertexv(I->ray, t0);
        v1 += 3;
        v2 += 3;
        v3 += 3;
      }
      CGOEnd(I->ray);

      CGOEnable(I->ray, GL_LIGHTING);
    }
  }
  if(ex) {
    ExtrudeFree(ex);
  }
  /* draw the rings */

  if(ring_anchor && n_ring) {
    int ring_i;
    int mem[8];
    int nbr[7];
    int *neighbor;
    int *marked = Calloc(int, obj->NAtom);
    float *moved = Calloc(float, obj->NAtom * 3);

    register int escape_count;
    register int *atmToIdx = NULL;

    if(!obj->DiscreteFlag)
      atmToIdx = cs->AtmToIdx;

    ObjectMoleculeUpdateNeighbors(obj);
    neighbor = obj->Neighbor;

    escape_count = ESCAPE_MAX;  /* don't get bogged down with structures 
                                   that have unreasonable connectivity */
    for(ring_i = 0; ring_i < n_ring; ring_i++) {
      mem[0] = ring_anchor[ring_i];
      nbr[0] = neighbor[mem[0]] + 1;
      while(((mem[1] = neighbor[nbr[0]]) >= 0) &&
            ((!atmToIdx) || (atmToIdx[mem[0]] >= 0))) {
        nbr[1] = neighbor[mem[1]] + 1;
        while(((mem[2] = neighbor[nbr[1]]) >= 0) &&
              ((!atmToIdx) || (atmToIdx[mem[1]] >= 0))) {
          if(mem[2] != mem[0]) {
            nbr[2] = neighbor[mem[2]] + 1;
            while(((mem[3] = neighbor[nbr[2]]) >= 0) &&
                  ((!atmToIdx) || (atmToIdx[mem[2]] >= 0))) {
              if(mem[3] != mem[1]) {
                nbr[3] = neighbor[mem[3]] + 1;
                while(((mem[4] = neighbor[nbr[3]]) >= 0) &&
                      ((!atmToIdx) || (atmToIdx[mem[3]] >= 0))) {
                  if((mem[4] != mem[2]) && (mem[4] != mem[1]) && (mem[4] != mem[0])) {
                    nbr[4] = neighbor[mem[4]] + 1;
                    while(((mem[5] = neighbor[nbr[4]]) >= 0) &&
                          ((!atmToIdx) || (atmToIdx[mem[4]] >= 0))) {
                      if(!(escape_count--))
                        goto escape;
                      if((mem[5] != mem[3]) && (mem[5] != mem[2]) && (mem[5] != mem[1])) {
                        if(mem[5] == mem[0]) {  /* five-cycle */
                          /*    printf(" 5: %s(%d) %s(%d) %s(%d) %s(%d) %s(%d)\n",
                             obj->AtomInfo[mem[0]].name,mem[0],
                             obj->AtomInfo[mem[1]].name,mem[1],
                             obj->AtomInfo[mem[2]].name,mem[2],
                             obj->AtomInfo[mem[3]].name,mem[3],
                             obj->AtomInfo[mem[4]].name,mem[4]); */
                          do_ring(G, 5, mem, obj, cs, ring_width, I->ray, ring_color,
                                  ring_mode, ladder_radius, ladder_color, ladder_mode,
                                  ring_finder, cartoon_side_chain_helper, nuc_flag,
                                  na_mode, ring_alpha, alpha, marked, moved, ring_radius);

                        }

                        nbr[5] = neighbor[mem[5]] + 1;
                        while(((mem[6] = neighbor[nbr[5]]) >= 0) &&
                              ((!atmToIdx) || (atmToIdx[mem[5]] >= 0))) {
                          if((mem[6] != mem[4]) && (mem[6] != mem[3])
                             && (mem[6] != mem[2]) && (mem[6] != mem[1])) {
                            if(mem[6] == mem[0]) {      /* six-cycle */
                              /* printf(" 6: %s %s %s %s %s %s\n",
                                 obj->AtomInfo[mem[0]].name,
                                 obj->AtomInfo[mem[1]].name,
                                 obj->AtomInfo[mem[2]].name,
                                 obj->AtomInfo[mem[3]].name,
                                 obj->AtomInfo[mem[4]].name,
                                 obj->AtomInfo[mem[5]].name
                                 ); */
                              do_ring(G, 6, mem, obj, cs, ring_width, I->ray, ring_color,
                                      ring_mode, ladder_radius, ladder_color, ladder_mode,
                                      ring_finder, cartoon_side_chain_helper, nuc_flag,
                                      na_mode, ring_alpha, alpha, marked, moved,
                                      ring_radius);
                            }
                            nbr[6] = neighbor[mem[6]] + 1;
                            while(((mem[7] = neighbor[nbr[6]]) >= 0) &&
                                  ((!atmToIdx) || (atmToIdx[mem[6]] >= 0))) {
                              if((mem[7] != mem[5]) && (mem[7] != mem[4])
                                 && (mem[7] != mem[3]) && (mem[7] != mem[2])
                                 && (mem[7] != mem[1])) {
                                if(mem[7] == mem[0]) {
                                  do_ring(G, 7, mem, obj, cs, ring_width, I->ray,
                                          ring_color, ring_mode, ladder_radius,
                                          ladder_color, ladder_mode, ring_finder,
                                          cartoon_side_chain_helper, nuc_flag, na_mode,
                                          ring_alpha, alpha, marked, moved, ring_radius);
                                }
                              }
                              nbr[6] += 2;
                            }
                          }
                          nbr[5] += 2;
                        }
                      }
                      nbr[4] += 2;
                    }
                  }
                  nbr[3] += 2;
                }
              }
              nbr[2] += 2;
            }
          }
          nbr[1] += 2;
        }
        nbr[0] += 2;
      }
    escape:
      escape_count = ESCAPE_MAX;        /* don't get bogged down with structures 
                                           that have unreasonable connectivity */
    }
    FreeP(marked);
    FreeP(moved);
  }
  CGOStop(I->ray);
  I->std = CGOSimplify(I->ray, 0);
  FreeP(dv);
  FreeP(dl);
  FreeP(tv);
  FreeP(nv);
  FreeP(at);
  FreeP(seg);
  FreeP(pv);
  FreeP(pvo);
  FreeP(pva);
  FreeP(car);
  FreeP(tmp);
  FreeP(sstype);
  FreeP(sampling_tmp);
  FreeP(flag_tmp);
  FreeP(nuc_flag);
  VLAFreeP(ring_anchor);
  return ((void *) (struct Rep *) I);
}
