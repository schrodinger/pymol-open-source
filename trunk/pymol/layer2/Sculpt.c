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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"
#include"OOMac.h"
#include"Feedback.h"
#include"Util.h"
#include"Sculpt.h"
#include"SculptCache.h"

#ifndef R_SMALL8
#define R_SMALL8 0.00000001
#endif

#define NB_HASH_SIZE 262144
#define EX_HASH_SIZE 65536

#define nb_hash(v) \
(((((int)*(v  ))>> 2)&0x0003F)|\
 ((((int)*(v+1))<< 4)&0x00FC0)|\
 ((((int)*(v+2))<<10)&0x3F000))

#define nb_hash_off(v,d,e,f) \
(((((d)+(int)*(v  ))>> 2)&0x0003F)|\
 ((((e)+(int)*(v+1))<< 4)&0x00FC0)|\
 ((((f)+(int)*(v+2))<<10)&0x3F000))

#define ex_hash(a,b) \
((((a)    )&0x00FF)|\
 (((b)<< 8)&0xFF00))

int SculptCheckBump(float *v1,float *v2,float *diff,float *dist,float cutoff);
int SculptDoBump(float target,float actual,float *d,float *d0to1,float *d1to0,float wt);

CSculpt *SculptNew(void)
{
  OOAlloc(CSculpt);
  I->Shaker = ShakerNew();
  I->NBList = VLAlloc(int,150000);
  I->NBHash = Alloc(int,NB_HASH_SIZE);
  I->EXList = VLAlloc(int,100000);
  I->EXHash = Alloc(int,EX_HASH_SIZE);
  I->Don = VLAlloc(int,1000);
  I->Acc = VLAlloc(int,1000);
  return(I);
}

void SculptMeasureObject(CSculpt *I,ObjectMolecule *obj,int state)
{
  int a,a0,a1,a2,a3,b0,b1,b2,b3;
  BondType *b;
  float *v0,*v1,*v2,*v3,d,dummy;
  CoordSet *cs;
  int n0,n1,n2;
  int *planer = NULL;
  int *linear = NULL;
  int nex = 1;
  int *j,*k,xhash;
  int ex_type;
  AtomInfoType *ai,*ai1,*ai2,*oai;
  int xoffset;
  int use_cache = 1;

  PRINTFD(FB_Sculpt)
    " SculptMeasureObject-Debug: entered.\n"
    ENDFD;

  ShakerReset(I->Shaker);

  UtilZeroMem(I->NBHash,NB_HASH_SIZE*sizeof(int));
  UtilZeroMem(I->EXHash,EX_HASH_SIZE*sizeof(int));

  if(state<obj->NCSet)
    if(obj->CSet[state])
      {
        oai = obj->AtomInfo;

        VLACheck(I->Don,int,obj->NAtom);
        VLACheck(I->Acc,int,obj->NAtom);
        ai = obj->AtomInfo;
        for(a=0;a<obj->NAtom;a++) {
          I->Don[a]=false;
          I->Acc[a]=false;
          if(!ai->sculpt_id) /* insure all atoms have unique sculpt IDs */
            ai->sculpt_id=SculptCacheNewID();
          ai++;
        }
        
        ObjectMoleculeVerifyChemistry(obj);
        ObjectMoleculeUpdateNeighbors(obj);

        cs = obj->CSet[state];

        use_cache = SettingGet_i(cs->Setting,obj->Obj.Setting,cSetting_sculpt_memory);
        if(obj->NBond) {

          planer=Alloc(int,obj->NAtom);
          linear=Alloc(int,obj->NAtom);
          ai = obj->AtomInfo;
          for(a=0;a<obj->NAtom;a++) {
            planer[a]=(ai->geom==cAtomInfoPlaner);
            linear[a]=(ai->geom==cAtomInfoLinear);
            ai++;
          }

          b=obj->Bond;
          for(a=0;a<obj->NBond;a++)
            {
              b1 = b->index[0];
              b2 = b->index[1];
              
              /* brain-dead donor/acceptor assignment
               * REPLACE later on with pattern-based system */

              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;
                            
              if(ai1->protons==cAN_H) {

                /* donors: any H attached to O, N */
                switch(ai2->protons) {
                case cAN_O: 
                  I->Don[b1]=true; 
                  I->Don[b2]=true; /* mark heavy atom too... */
                  break;
                case cAN_N: 
                  I->Don[b1]=true; 
                  I->Don[b2]=true;
                  break;
                }
              } else if(ai2->protons==cAN_H) {
                switch(ai1->protons) {
                case cAN_O: 
                  I->Don[b1]=true; 
                  I->Don[b2]=true; /* mark heavy atom too... */
                  break;
                case cAN_N: 
                  I->Don[b1]=true; 
                  I->Don[b2]=true; /* mark heavy atom too... */
                  break;
                }
              } else {
                
                /* assume O is always an acceptor...*/

                if(ai1->protons==cAN_O) I->Acc[b1]=true;
                if(ai2->protons==cAN_O) I->Acc[b2]=true;

                /* nitrogens with lone pairs are acceptors */

                if(b->order==2) {
                  if(ai1->protons==cAN_N) {
                    if(b->order==2) {
                      n1 = obj->Neighbor[b1];
                      if(obj->Neighbor[n1]<3) { /* N with L.P. */
                        I->Acc[b1]=true;
                      }
                    }
                  }
                }
                if(ai2->protons==cAN_N) {
                  if(b->order==2) {
                    n2 = obj->Neighbor[b2];
                    if(obj->Neighbor[n2]<3) { /* N with L.P. */
                      I->Acc[b2]=true;
                    }
                  }
                }
              }

              /*              if(Feedback(FB_Sculpt,FB_Debugging)) {
                              for(a=0;a<obj->NAtom;a++) {
                              I->Don[a]=false;
                              I->Acc[a]=false;
                              }
                              
                              }*/


              /* exclusions */

              xhash = ( (b2>b1) ? ex_hash(b1,b2) : ex_hash(b2,b1));
              VLACheck(I->EXList,int,nex+3);
              j = I->EXList+nex;
              *(j++)=*(I->EXHash+xhash);
              if(b2>b1) {
                *(j++)=b1;
                *(j++)=b2;
              } else {
                *(j++)=b2;
                *(j++)=b1;
              }
              *(j++)=2; /* 1-2 exclusion */
              *(I->EXHash+xhash)=nex;
              nex+=4;


              b++;
              if(obj->DiscreteFlag) {
                if((cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
                  a1=obj->DiscreteAtmToIdx[b1];
                  a2=obj->DiscreteAtmToIdx[b2];
                } else {
                  a1=-1;
                  a2=-1;
                }
              } else {
                a1=cs->AtmToIdx[b1];
                a2=cs->AtmToIdx[b2];
              }
              if((a1>=0)&&(a2>=0))
                {
                  v1 = cs->Coord+3*a1;
                  v2 = cs->Coord+3*a2;
                  d = (float)diff3f(v1,v2);
                  if(use_cache) {
                    if(!SculptCacheQuery(cSculptBond,
                                         oai[b1].sculpt_id,
                                         oai[b2].sculpt_id,0,0,&d))
                      SculptCacheStore(cSculptBond,
                                       oai[b1].sculpt_id,
                                       oai[b2].sculpt_id,0,0,d);
                  }
                  ShakerAddDistCon(I->Shaker,b1,b2,d,cShakerDistBond); 
                  /* NOTE: storing atom indices, not coord. ind.! */
                }
            }
          
          /* now pick up those 1-3 interations */
          
          /* b1-b0-b2 */

          for(b0=0;b0<obj->NAtom;b0++) {
            n0 = obj->Neighbor[b0]+1;
            while(obj->Neighbor[n0]>=0) {
              b1 = obj->Neighbor[n0];
              n1 = n0+2;
              while(obj->Neighbor[n1]>=0) {
                b2 = obj->Neighbor[n1];


                xhash = ( (b2>b1) ? ex_hash(b1,b2) : ex_hash(b2,b1));
                VLACheck(I->EXList,int,nex+3);
                j = I->EXList+nex;
                *(j++)=*(I->EXHash+xhash);
                if(b2>b1) {
                  *(j++)=b1;
                  *(j++)=b2;
                } else {
                  *(j++)=b2;
                  *(j++)=b1;
                }
                *(j++)=3; /* 1-3 exclusion */
                *(I->EXHash+xhash)=nex;
                nex+=4;

                if(obj->DiscreteFlag) {
                  if((cs==obj->DiscreteCSet[b0])&&(cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
                    a0=obj->DiscreteAtmToIdx[b0];
                    a1=obj->DiscreteAtmToIdx[b1];
                    a2=obj->DiscreteAtmToIdx[b2];
                  } else {
                    a0=-1;
                    a1=-1;
                    a2=-1;
                  }
                } else {
                  a0=cs->AtmToIdx[b0];
                  a1=cs->AtmToIdx[b1];
                  a2=cs->AtmToIdx[b2];
                }

                if((a0>=0)&&(a1>=0)&&(a2>=0)) {
                  v1 = cs->Coord+3*a1;
                  v2 = cs->Coord+3*a2;
                  d = (float)diff3f(v1,v2);
                  if(use_cache) {
                    if(!SculptCacheQuery(cSculptAngl,
                                         oai[b0].sculpt_id,
                                         oai[b1].sculpt_id,
                                         oai[b2].sculpt_id,0,&d))
                      SculptCacheStore(cSculptAngl,
                                       oai[b0].sculpt_id,
                                       oai[b1].sculpt_id,
                                       oai[b2].sculpt_id,0,d);
                  }
                  ShakerAddDistCon(I->Shaker,b1,b2,d,cShakerDistAngle); 


                  if(linear[b0]&&(linear[b1]||linear[b2])) {
                    
                    if(use_cache) {
                      if(!SculptCacheQuery(cSculptLine,
                                           oai[b1].sculpt_id,
                                           oai[b0].sculpt_id,
                                           oai[b2].sculpt_id,0,&dummy))
                        SculptCacheStore(cSculptLine,
                                         oai[b1].sculpt_id,
                                         oai[b0].sculpt_id,
                                         oai[b2].sculpt_id,0,0.0);
                    }
                    ShakerAddLineCon(I->Shaker,b1,b0,b2); 
                  }
                }
                n1+=2;
              }
              n0+=2;
            }
          }
          
#if 0
          /* add additional limiting distance restraints into object */

          /* b1-b0-b2-b3-b4 */

          for(b0=0;b0<obj->NAtom;b0++) {
            n0 = obj->Neighbor[b0]+1;
            while(obj->Neighbor[n0]>=0) {
              b1 = obj->Neighbor[n0];
              n1 = obj->Neighbor[b0]+1;
              while(obj->Neighbor[n1]>=0) {
                b2 = obj->Neighbor[n1];
                if(b1!=b2) {
                  n2 =  obj->Neighbor[b2]+1;
                  while(obj->Neighbor[n2]>=0) {
                    b3 = obj->Neighbor[n2];
                    if(b3!=b0) {
                      n3 =  obj->Neighbor[b3]+1;
                      while(obj->Neighbor[n3]>=0) {
                        b4 = obj->Neighbor[n3];
                        if(b2!=b4) {
                          
                          if(obj->DiscreteFlag) {
                            if((cs==obj->DiscreteCSet[b0])&&
                               (cs==obj->DiscreteCSet[b1])&&
                               (cs==obj->DiscreteCSet[b2])&&
                               (cs==obj->DiscreteCSet[b3])&&
                               (cs==obj->DiscreteCSet[b4])
                               ) {
                              a0=obj->DiscreteAtmToIdx[b0];
                              a1=obj->DiscreteAtmToIdx[b1];
                              a2=obj->DiscreteAtmToIdx[b2];
                              a3=obj->DiscreteAtmToIdx[b3];
                              a4=obj->DiscreteAtmToIdx[b4];
                            } else {
                              a0=-1;
                              a1=-1;
                              a2=-1;
                              a3=-1;
                              a4=-1;
                            }
                          } else {
                            a0=cs->AtmToIdx[b0];
                            a1=cs->AtmToIdx[b1];
                            a2=cs->AtmToIdx[b2];
                            a3=cs->AtmToIdx[b3];
                            a4=cs->AtmToIdx[b4];
                          }
                          if((a0>=0)&&(a1>=0)&&(a2>=0)&&(a3>=0)&&(a4>=0)) {
                            v0 = cs->Coord+3*a0;
                            v1 = cs->Coord+3*a1;
                            v2 = cs->Coord+3*a2;
                            v3 = cs->Coord+3*a3;
                            v4 = cs->Coord+3*a4;
                            d = diff3f(v0,v1)+diff3f(v0,v2)+diff3f(v2,v3)+diff3f(v3,v4);
                            ShakerAddDistCon(I->Shaker,b0,b1,d,cShakerDistLimit); 
                          }   
                        }
                        n3+=2;
                      }                
                    }
                    n2+=2;
                  }
                }
                n1+=2;
              }
              n0+=2;
            }
          }
#endif

          /* and record the pyramidal and planer geometries */
          
          for(b0=0;b0<obj->NAtom;b0++) {
            n0 = obj->Neighbor[b0]+1;
            while(obj->Neighbor[n0]>=0) {
              b1 = obj->Neighbor[n0];
              n1 = n0+2;
              while(obj->Neighbor[n1]>=0) {
                b2 = obj->Neighbor[n1];
                n2 = n1+2;
                while(obj->Neighbor[n2]>=0) {
                  b3 = obj->Neighbor[n2];
                  
                                
                  if(obj->DiscreteFlag) {
                    if((cs==obj->DiscreteCSet[b0])&&
                       (cs==obj->DiscreteCSet[b1])&&
                       (cs==obj->DiscreteCSet[b2])) {
                      a0=obj->DiscreteAtmToIdx[b0];
                      a1=obj->DiscreteAtmToIdx[b1];
                      a2=obj->DiscreteAtmToIdx[b2];
                      a3=obj->DiscreteAtmToIdx[b3];
                    } else {
                      a0=-1;
                      a1=-1;
                      a2=-1;
                      a3=-1;
                    }
                  } else {
                    a0=cs->AtmToIdx[b0];
                    a1=cs->AtmToIdx[b1];
                    a2=cs->AtmToIdx[b2];
                    a3=cs->AtmToIdx[b3];
                  }
                  if((a0>=0)&&(a1>=0)&&(a2>=0)&&(a3>=0)) {
                    v0 = cs->Coord+3*a0;
                    v1 = cs->Coord+3*a1;
                    v2 = cs->Coord+3*a2;
                    v3 = cs->Coord+3*a3;
                    d = ShakerGetPyra(v0,v1,v2,v3);
                    if(fabs(d)<0.05) {
                      planer[b0]=true;
                    }
                    if(planer[b0])
                      d=0.0;
                    if(use_cache) {
                      if(!SculptCacheQuery(cSculptPyra,
                                           oai[b1].sculpt_id,
                                           oai[b0].sculpt_id,
                                           oai[b2].sculpt_id,
                                           oai[b3].sculpt_id,
                                           &d))
                        SculptCacheStore(cSculptPyra,
                                         oai[b1].sculpt_id,
                                         oai[b0].sculpt_id,
                                         oai[b2].sculpt_id,
                                         oai[b3].sculpt_id,
                                         d);
                    }
                    ShakerAddPyraCon(I->Shaker,b0,b1,b2,b3,d); 
                  }                
                  n2+=2;
                }
                n1+=2;
              }
              n0+=2;
            }
          }

          /* b1\b0_b2/b3 */
        
          for(b0=0;b0<obj->NAtom;b0++) {
            n0 = obj->Neighbor[b0]+1;
            while(obj->Neighbor[n0]>=0) {
              b1 = obj->Neighbor[n0];

              n1 = obj->Neighbor[b0]+1;
              while(obj->Neighbor[n1]>=0) {
                b2 = obj->Neighbor[n1];

                if(b1!=b2) {
                  n2 =  obj->Neighbor[b2]+1;
                  while(obj->Neighbor[n2]>=0) {
                    b3 = obj->Neighbor[n2];
                    if((b3!=b0)&&(b3>b1)) {
                      
                      /* check 1-4 exclusion */
                      xhash = ex_hash(b1,b3);
                      
                      ex_type = 4; 
                      
                      xoffset = *(I->EXHash+xhash);
                      while(xoffset) {
                        k = I->EXList + xoffset;
                        if((abs(*(k+3))==4)&&
                           (*(k+1)==b1)&&
                           (*(k+2)==b3)) {
                          if((b0!=*(k+4))&&
                             (b2!=*(k+5))) {
                            if(planer[b0]&&planer[b2]&&
                               planer[*(k+4)]&&planer[*(k+5)]) {
                              /* aromatic */ 
                              *(k+3)=-4;
                            }
                          }
                          ex_type = 0; /* duplicate, skip */
                          break;
                        }
                        xoffset = *k;
                      }
                      if(ex_type) {
                        VLACheck(I->EXList,int,nex+5);
                        j = I->EXList+nex;
                        *(j++)=*(I->EXHash+xhash);
                        *(j++)=b1;
                        *(j++)=b3;
                        *(j++)=ex_type;
                        *(j++)=b0;
                        *(j++)=b2;
                        *(I->EXHash+xhash)=nex;
                        nex+=6;
                      }
                      
                      /* planarity */
                      
                      if(obj->DiscreteFlag) {
                        if((cs==obj->DiscreteCSet[b0])&&
                           (cs==obj->DiscreteCSet[b1])&&
                           (cs==obj->DiscreteCSet[b2])&&
                           (cs==obj->DiscreteCSet[b3])) {
                          a0=obj->DiscreteAtmToIdx[b0];
                          a1=obj->DiscreteAtmToIdx[b1];
                          a2=obj->DiscreteAtmToIdx[b2];
                          a3=obj->DiscreteAtmToIdx[b3];
                        } else {
                          a0=-1;
                          a1=-1;
                          a2=-1;
                          a3=-1;
                        }
                      } else {
                        a0=cs->AtmToIdx[b0];
                        a1=cs->AtmToIdx[b1];
                        a2=cs->AtmToIdx[b2];
                        a3=cs->AtmToIdx[b3];
                      }
                      if((a0>=0)&&(a1>=0)&&(a2>=0)&&(a3>=0)) {
                        v0 = cs->Coord+3*a0;
                        v1 = cs->Coord+3*a1;
                        v2 = cs->Coord+3*a2;
                        v3 = cs->Coord+3*a3;
                        
                        d = 0.0;
                        if(fabs(get_dihedral3f(v1,v0,v2,v3))<deg_to_rad(10.0))
                          if(planer[b0]&&planer[b2]) {
                            d = 1.0; 
                          }
                        if(use_cache) {
                          if(!SculptCacheQuery(cSculptPlan,
                                               oai[b1].sculpt_id,
                                               oai[b0].sculpt_id,
                                               oai[b2].sculpt_id,
                                               oai[b3].sculpt_id,
                                               &d))
                            SculptCacheStore(cSculptPlan,
                                             oai[b1].sculpt_id,
                                             oai[b0].sculpt_id,
                                             oai[b2].sculpt_id,
                                             oai[b3].sculpt_id,
                                             d);
                        }
                        if(d>0.0) {
                          ShakerAddPlanCon(I->Shaker,b1,b0,b2,b3); 
                        }

                      }
                    }
                    n2+=2;
                  }
                }
                n1+=2;
              }
              n0+=2;
            }
          }
          FreeP(planer);
          FreeP(linear);
        }
      }

  PRINTFB(FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NDistCon %d\n",I->Shaker->NDistCon
    ENDFB;
  PRINTFB(FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NPyraCon %d\n",I->Shaker->NPyraCon
    ENDFB;
  PRINTFB(FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NPlanCon %d\n",I->Shaker->NPlanCon
    ENDFB;


 PRINTFD(FB_Sculpt)
    " SculptMeasureObject-Debug: leaving...\n"
    ENDFD;


}

int SculptCheckBump(float *v1,float *v2,float *diff,float *dist,float cutoff)
{
  register float d2;
  diff[0] = (v1[0]-v2[0]);
  if(fabs(diff[0])>cutoff) return(false);
  diff[1] = (v1[1]-v2[1]);
  if(fabs(diff[1])>cutoff) return(false);
  diff[2] = (v1[2]-v2[2]);
  if(fabs(diff[2])>cutoff) return(false);
  d2 = (diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
  if(d2<(cutoff*cutoff)) {
    *dist = (float)sqrt(d2);
    return(true);
  }
  return(false);
}

int SculptDoBump(float target,float actual,float *d,float *d0to1,float *d1to0,float wt)
{
  float push[3];
  float dev,dev_2,sc;

  dev = target-actual;
  if(fabs(dev)>R_SMALL8) {
    dev_2 = wt*dev/2.0F;
    if(actual>R_SMALL8) { /* nonoverlapping */
      sc = dev_2/actual;
      scale3f(d,sc,push);
      add3f(push,d0to1,d0to1);
      subtract3f(d1to0,push,d1to0);
    } else { /* overlapping, so just push along X */
      d0to1[0]-=dev_2;
      d1to0[0]+=dev_2;
    }
    return 1;
  }
  return 0;
}
       
void SculptIterateObject(CSculpt *I,ObjectMolecule *obj,int state,int n_cycle)
{
  CShaker *shk;
  int a,aa,a0,a1,a2,a3,b0,b1,b2,b3;
  CoordSet *cs;
  ShakerDistCon *sdc;
  ShakerPyraCon *spc;
  ShakerPlanCon *snc;
  ShakerLineCon *slc;

  float *disp = NULL;
  float *v,*v0,*v1,*v2,*v3;
  float diff[3],len;
  int *atm2idx = NULL;
  int *cnt = NULL;
  int *i,*j;
  float sc;
  int hash;
  int nb_next;
  int h,k,l;
  int offset,xoffset;
  float cutoff;
  int ex,ex1;
  int eval_flag;
  int mask;
  float wt;
  float vdw;
  float vdw14;
  float vdw_wt;
  float vdw_wt14;
  float bond_wt;
  float angl_wt;
  float pyra_wt;
  float plan_wt;
  float line_wt;
  int active_flag=false;
  float hb_overlap,hb_overlap_base;
  int *active,n_active;
  AtomInfoType *ai0,*ai1;
  double task_time;

  PRINTFD(FB_Sculpt)
    " SculptIterateObject-Debug: entered state=%d n_cycle=%d\n",state,n_cycle
    ENDFD;

  if(state<obj->NCSet)
    if(obj->CSet[state]&&n_cycle)
      {
        disp = Alloc(float,3*obj->NAtom);
        atm2idx = Alloc(int,obj->NAtom);
        cnt = Alloc(int,obj->NAtom);
        active = Alloc(int,obj->NAtom);
        shk=I->Shaker;

        PRINTFD(FB_Sculpt)
          " SIO-Debug: NDistCon %d\n",shk->NDistCon
          ENDFD;

        cs = obj->CSet[state];

        vdw = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_scale);
        vdw14 = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_scale14);
        vdw_wt = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_weight);
        vdw_wt14 = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_weight14);
        bond_wt =  SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_bond_weight);
        angl_wt =  SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_angl_weight);
        pyra_wt =  SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_pyra_weight);
        plan_wt =  SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_plan_weight);
        line_wt =  SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_line_weight);
        mask = SettingGet_i(cs->Setting,obj->Obj.Setting,cSetting_sculpt_field_mask);
        hb_overlap = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_hb_overlap);
        hb_overlap_base = SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_hb_overlap_base);

        n_active = 0;
        ai0=obj->AtomInfo;
        for(a=0;a<obj->NAtom;a++) {
          if(ai0->flags&cAtomFlag_exclude) {
            a1=-1;
          } else {
            if(obj->DiscreteFlag) {
              if(cs==obj->DiscreteCSet[a]) {
                a1=obj->DiscreteAtmToIdx[a];
              } else {
                a1=-1;
              }
            } else {
              a1=cs->AtmToIdx[a];
            }
          }
          if(a1>=0) {
            active_flag=true;
            active[n_active] = a;
            n_active++;
          }
          atm2idx[a] = a1;
          ai0++;
        }
        
        if(active_flag) {

          /* first, create coordinate -> vertex mapping */
          /* and count number of constraints */

          task_time = UtilGetSeconds();
          while(n_cycle--) {

            /* initialize displacements to zero */
        
            v = disp;
            i = cnt;
            for(aa=0;aa<n_active;aa++) {
              a = active[aa];
              v=disp+a*3;
              *(v++)=0.0;
              *(v++)=0.0;          
              *(v++)=0.0;
              cnt[a]=0;
            }
          
            /* apply distance constraints */

            sdc=shk->DistCon;
            for(a=0;a<shk->NDistCon;a++) {
              b1 = sdc->at0;
              b2 = sdc->at1;
              
              switch(sdc->type) {
              case cShakerDistBond:
                wt = bond_wt;
                eval_flag = cSculptBond & mask;
                break;
              case cShakerDistAngle:
                wt = angl_wt;
                eval_flag = cSculptAngl & mask;
                break;
              case cShakerDistLimit:
                wt = 2.0F;
                eval_flag = true;
                break;
              default:
                eval_flag = false;
                wt=0.0F;
                break;
              }
              
              if(eval_flag) {
                a1 = atm2idx[b1]; /* coordinate set indices */
                a2 = atm2idx[b2];
                cnt[b1]++;
                cnt[b2]++;
              
                if((a1>=0)&&(a2>=0))
                  {
                    v1 = cs->Coord+3*a1;
                    v2 = cs->Coord+3*a2;
                    if(sdc->type!=cShakerDistLimit) {
                      ShakerDoDist(sdc->targ,v1,v2,disp+b1*3,disp+b2*3,wt);
                    } else {
                      ShakerDoDistLimit(sdc->targ,v1,v2,disp+b1*3,disp+b2*3,wt);
                    }
                  }
              }
              sdc++;
            }

            /* apply line constraints */
            
            if(cSculptLine & mask) {
              slc=shk->LineCon;
              
              for(a=0;a<shk->NLineCon;a++) {
                b0 = slc->at0;
                b1 = slc->at1;
                b2 = slc->at2;
                a0 = atm2idx[b0]; /* coordinate set indices */
                a1 = atm2idx[b1];
                a2 = atm2idx[b2];
                
                cnt[b0]++;
                cnt[b1]++;
                cnt[b2]++;
                
                if((a0>=0)&&(a1>=0)&&(a2>=0))
                  {
                    v0 = cs->Coord+3*a0;
                    v1 = cs->Coord+3*a1;
                    v2 = cs->Coord+3*a2;
                    ShakerDoLine(v0,v1,v2,disp+b0*3,disp+b1*3,disp+b2*3,line_wt);
                  }
                slc++;
              }
            }

            /* apply pyramid constraints */

            if(cSculptPyra & mask) {
              spc=shk->PyraCon;
              for(a=0;a<shk->NPyraCon;a++) {
              
                b0 = spc->at0;
                b1 = spc->at1;
                b2 = spc->at2;
                b3 = spc->at3;
                a0 = atm2idx[b0];
                a1 = atm2idx[b1];
                a2 = atm2idx[b2];
                a3 = atm2idx[b3];
              
                if((a0>=0)&&(a1>=0)&&(a2>=0)&&(a3>=0)) {
                  v0 = cs->Coord+3*a0;
                  v1 = cs->Coord+3*a1;
                  v2 = cs->Coord+3*a2;
                  v3 = cs->Coord+3*a3;
                  ShakerDoPyra(spc->targ,
                               v0,v1,v2,v3,
                               disp+b0*3,
                               disp+b1*3,
                               disp+b2*3,
                               disp+b3*3,
                               pyra_wt);
                
                
                  cnt[b0]++;
                  cnt[b1]++;
                  cnt[b2]++;
                  cnt[b3]++;
                }
                spc++;
              }
            }

            if(cSculptPlan & mask) {

              /* apply planarity constraints */
            
              snc=shk->PlanCon;
              for(a=0;a<shk->NPlanCon;a++) {
              
                b0 = snc->at0;
                b1 = snc->at1;
                b2 = snc->at2;
                b3 = snc->at3;
                a0 = atm2idx[b0];
                a1 = atm2idx[b1];
                a2 = atm2idx[b2];
                a3 = atm2idx[b3];
              
                if((a0>=0)&&(a1>=0)&&(a2>=0)&&(a3>=0)) {
                  v0 = cs->Coord+3*a0;
                  v1 = cs->Coord+3*a1;
                  v2 = cs->Coord+3*a2;
                  v3 = cs->Coord+3*a3;
                  ShakerDoPlan(v0,v1,v2,v3,
                               disp+b0*3,
                               disp+b1*3,
                               disp+b2*3,
                               disp+b3*3,
                               plan_wt);
                  cnt[b0]++;
                  cnt[b1]++;
                  cnt[b2]++;
                  cnt[b3]++;
                }
              
                snc++;
              }
            }

            if((cSculptVDW|cSculptVDW14)&mask) {
              /* compute non-bonded interations */
            
              /* construct nonbonded hash */
            
              nb_next = 1;
              for(aa=0;aa<n_active;aa++) {
                b0 = active[aa];
                a0 = atm2idx[b0];
                VLACheck(I->NBList,int,nb_next+2);
                v0 = cs->Coord+3*a0;
                hash = nb_hash(v0);
                i = I->NBList+nb_next;
                *(i++)=*(I->NBHash+hash);
                *(i++)=hash;
                *(i++)=b0;
                *(I->NBHash+hash)=nb_next;
                nb_next+=3;
              }
            
              /* find neighbors for each atom */
            
              for(aa=0;aa<n_active;aa++) {
                b0 = active[aa];
                a0 = atm2idx[b0];
                ai0=obj->AtomInfo+b0;
                v0 = cs->Coord+3*a0;
                for(h=-4;h<5;h+=4)
                  for(k=-4;k<5;k+=4)
                    for(l=-4;l<5;l+=4) 
                      {
                        offset = *(I->NBHash+nb_hash_off(v0,h,k,l));
                        while(offset) {
                          i = I->NBList + offset;
                          b1 = *(i+2);
                          if(b1>b0) { 
                            /* determine exclusion (if any) */
                            xoffset = *(I->EXHash+ex_hash(b0,b1));
                            ex = 10;
                            while(xoffset) {
                              j = I->EXList + xoffset;
                              if((*(j+1)==b0)&&
                                 (*(j+2)==b1)) {
                                ex1 = *(j+3);
                                if(ex1<ex) {
                                  ex=ex1;
                                }
                              }
                              xoffset = (*j);
                            }
                            if(ex>3) {
                              ai1=obj->AtomInfo+b1;
                              cutoff = ai0->vdw+ai1->vdw;
                              if(ex==4) { /* 1-4 interation */
                                cutoff*=vdw14;
                                wt = vdw_wt14;
                                eval_flag = cSculptVDW14 & mask;
                              } else { /* standard interaction */
                                if(I->Don[b0]&&I->Acc[b1]) { /* h-bond */
                                  if(ai0->protons==cAN_H) {
                                    cutoff-=hb_overlap;
                                  } else {
                                    cutoff-=hb_overlap_base;
                                  }
                                } else if(I->Acc[b0]&&I->Don[b1]) { /* h-bond */
                                  if(ai1->protons==cAN_H) {
                                    cutoff-=hb_overlap;
                                  } else {
                                    cutoff-=hb_overlap_base;
                                  } 
                                }
                                cutoff=cutoff*vdw;
                                wt = vdw_wt;
                                eval_flag = cSculptVDW & mask;
                              }
                              if(eval_flag) {
                                a1 = atm2idx[b1];
                                v1 = cs->Coord+3*a1;
                                if(SculptCheckBump(v0,v1,diff,&len,cutoff))
                                  if(SculptDoBump(cutoff,len,diff,
                                                  disp+b0*3,disp+b1*3,wt)) {
                                    cnt[b0]++;
                                    cnt[b1]++;
                                  }
                              }
                            }
                          }
                          offset=(*i);
                        }
                      }
              }
              
              /* clean up nonbonded hash */
              
              i = I->NBList+2;
              while(nb_next>1) {
                *(I->NBHash+*i)=0; 
                i+=3;
                nb_next-=3;
              }
              
            }
            /* average the displacements */
            
            for(aa=0;aa<n_active;aa++) {
              a = active[aa];
              if(cnt[a]) {
                if(!obj->AtomInfo[a].protekted) {
                  v1 = disp+3*a;
                  sc = 0.99F/cnt[a];
                  scale3f(v1,sc,v1);
                  v2 = cs->Coord+3*atm2idx[a];
                  add3f(v1,v2,v2);
                }
              }
            }
            if(cs->fInvalidateRep) {
              cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
            } else {
              ObjectMoleculeInvalidate(obj,cRepAll,cRepInvCoord);
            }
          
          }
          
          task_time = UtilGetSeconds() - task_time;
          PRINTFB(FB_Sculpt,FB_Blather)
            " Sculpt: %2.5f seconds\n",task_time
            ENDFB;
        }
        FreeP(active);
        FreeP(cnt);
        FreeP(disp);
        FreeP(atm2idx);
      }

  PRINTFD(FB_Sculpt)
    " SculptIterateObject-Debug: leaving...\n"
    ENDFD;

}

void SculptFree(CSculpt *I)
{
  VLAFreeP(I->Don);
  VLAFreeP(I->Acc);
  VLAFreeP(I->NBList);
  VLAFreeP(I->EXList);

  FreeP(I->NBHash);
  FreeP(I->EXHash);
  ShakerFree(I->Shaker);
  OOFreeP(I);
}

