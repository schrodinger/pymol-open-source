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

#include"os_std.h"
#include"os_gl.h"
#include"OOMac.h"
#include"Feedback.h"
#include"Util.h"
#include"Sculpt.h"

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
((((d+(int)*(v  ))>> 2)&0x0003F)|\
 (((e+(int)*(v+1))<< 4)&0x00FC0)|\
 (((f+(int)*(v+2))<<10)&0x3F000))

#define ex_hash(a,b) \
((((a)    )&0x00FF)|\
 (((b)<< 8)&0xFF00))

CSculpt *SculptNew(void)
{
  OOAlloc(CSculpt);
  I->Shaker = ShakerNew();
  I->NBList = VLAlloc(int,150000);
  I->NBHash = Alloc(int,NB_HASH_SIZE);
  I->EXList = VLAlloc(int,100000);
  I->EXHash = Alloc(int,EX_HASH_SIZE);

  UtilZeroMem(I->NBHash,NB_HASH_SIZE*sizeof(int));
  UtilZeroMem(I->EXHash,EX_HASH_SIZE*sizeof(int));
  return(I);
}

void SculptMeasureObject(CSculpt *I,ObjectMolecule *obj,int state)
{
  int a,a0,a1,a2,a3,b0,b1,b2,b3;
  BondType *b;
  float *v0,*v1,*v2,*v3,d;
  CoordSet *cs;
  int n0,n1,n2;
  int *planer = NULL;
  int nex = 1;
  int *j,xhash;

  AtomInfoType *ai;

  ShakerReset(I->Shaker);

  if(state<obj->NCSet)
    if(obj->CSet[state])
      {

        ObjectMoleculeVerifyChemistry(obj);
        ObjectMoleculeUpdateNeighbors(obj);

        cs = obj->CSet[state];

        if(obj->NBond) {

          planer=Alloc(int,obj->NAtom);
          ai = obj->AtomInfo;
          for(a=0;a<obj->NAtom;a++) {
            planer[a]=(ai->geom==3);
            ai++;
          }
          
          b=obj->Bond;
          for(a=0;a<obj->NBond;a++)
            {
              b1 = b->index[0];
              b2 = b->index[1];

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
                  d = diff3f(v1,v2);
                  /*
                  d = AtomInfoGetBondLength(obj->AtomInfo+b1,obj->AtomInfo+b2);
                  */
                  ShakerAddDistCon(I->Shaker,b1,b2,d,1); 
                  /* NOTE: storing atom indices, not coord. ind.! */
                }
            }
          
          /* now pick up those 1-3 interations */
          
          
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
                  d = diff3f(v1,v2);
                  ShakerAddDistCon(I->Shaker,b1,b2,d,0); 
                }
                n1+=2;
              }
              n0+=2;
            }
          }
          
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
              n1 = n0+2;
              while(obj->Neighbor[n1]>=0) {
                b2 = obj->Neighbor[n1];

                n2 =  obj->Neighbor[b2]+1;
                while(obj->Neighbor[n2]>=0) {
                  b3 = obj->Neighbor[n2];
                  if(b3!=b0) {
                    
                    xhash = ( (b3>b1) ? ex_hash(b1,b3) : ex_hash(b3,b1));
                    VLACheck(I->EXList,int,nex+3);
                    j = I->EXList+nex;
                    *(j++)=*(I->EXHash+xhash);
                    if(b2>b1) {
                      *(j++)=b1;
                      *(j++)=b3;
                    } else {
                      *(j++)=b3;
                      *(j++)=b1;
                    }
                    *(j++)=4; /* 1-4 exclusion */
                    *(I->EXHash+xhash)=nex;
                    nex+=4;


                  
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
                    
                      if(fabs(get_dihedral3f(v1,v0,v2,v3))<deg_to_rad(10.0))
                        if(planer[b0]&&planer[b2]) {
                          ShakerAddPlanCon(I->Shaker,b1,b0,b2,b3); 

                        }
                    }
                  }
                  n2+=2;
                }
                n1+=2;
              }
              n0+=2;
            }
          }
        FreeP(planer);
        }
      }

  
  PRINTFB(FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NDistCon %d\n",I->Shaker->NDistCon
    ENDFD;
  PRINTFB(FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NPyraCon %d\n",I->Shaker->NPyraCon
    ENDFD;
  PRINTFB(FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NPlanCon %d\n",I->Shaker->NPlanCon
    ENDFD;
}
int SculptCheckBump(float *v1,float *v2,float *diff,float *dist,float cutoff);
int SculptDoBump(float target,float actual,float *d,float *d0to1,float *d1to0,float wt);

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
    *dist = sqrt(d2);
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
    dev_2 = wt*dev/2.0;
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
  int a,a0,a1,a2,a3,b0,b1,b2,b3;
  CoordSet *cs;
  ShakerDistCon *sdc;
  ShakerPyraCon *spc;
  ShakerPlanCon *snc;

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

  AtomInfoType *ai0;

  PRINTFD(FB_Sculpt)
    " SculptIterateObject-Debug: entered state=%d n_cycle=%d\n",state,n_cycle
    ENDFD;

  if(state<obj->NCSet)
    if(obj->CSet[state]&&n_cycle)
      {
        disp = Alloc(float,3*obj->NAtom);
        atm2idx = Alloc(int,obj->NAtom);
        cnt = Alloc(int,obj->NAtom);
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
        mask =  SettingGet_f(cs->Setting,obj->Obj.Setting,cSetting_sculpt_field_mask);

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
          atm2idx[a] = a1;
          ai0++;
        }
        
        /* first, create coordinate -> vertex mapping */
        /* and count number of constraints */

        while(n_cycle--) {

          /* initialize displacements to zero */
        
          v = disp;
          i = cnt;
          for(a=0;a<obj->NAtom;a++) {
            if(atm2idx[a]) {
              *(v++)=0.0;
              *(i++)=0;
              *(v++)=0.0;          
              *(v++)=0.0;
            } else {
              i++;
              v+=3;
            }
          }
          
          /* apply distance constraints */

          sdc=shk->DistCon;
          for(a=0;a<shk->NDistCon;a++) {
            b1 = sdc->at0;
            b2 = sdc->at1;
            if(sdc->type) {
              wt = bond_wt;
              eval_flag = cSculptBond & mask;
            } else {
              wt = angl_wt;
              eval_flag = cSculptAngl & mask;
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
                  ShakerDoDist(sdc->targ,v1,v2,disp+b1*3,disp+b2*3,wt);
                }
            }
            sdc++;
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

          if((cSculptVDW||cSculptVDW14)&mask) {
            /* compute non-bonded interations */
            
            /* construct nonbonded hash */
            
            nb_next = 1;
            for(b0=0;b0<obj->NAtom;b0++) {
              a0 = atm2idx[b0];
              if(a0>=0) {
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
            }
            
            /* find neighbors for each atom */
            
            for(b0=0;b0<obj->NAtom;b0++) {
              a0 = atm2idx[b0];
              if(a0>=0) {
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
                              cutoff = ai0->vdw+obj->AtomInfo[b1].vdw;
                              if(ex==4) {
                                cutoff*=vdw14;
                                wt = vdw_wt14;
                                eval_flag = cSculptVDW14 & mask;
                              } else {
                                cutoff=cutoff*vdw;
                                wt = vdw_wt;
                                eval_flag = cSculptVDW & mask;
                              }
                              if(eval_flag) {
                                a1 = atm2idx[b1];
                                v1 = cs->Coord+3*a1;
                                if(SculptCheckBump(v0,v1,diff,&len,cutoff))
                                  if(SculptDoBump(cutoff,len,diff,disp+b0*3,disp+b1*3,wt)) {
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
            
          for(a=0;a<obj->NAtom;a++) {
            if(atm2idx[a]) 
              if(cnt[a]) {
                if(!obj->AtomInfo[a].protected) {
                  v1 = disp+3*a;
                  sc = 0.99/cnt[a];
                  scale3f(v1,sc,v1);
                  v2 = cs->Coord+3*atm2idx[a];
                  add3f(v1,v2,v2);
                }
              }
          }
          
          ObjectMoleculeInvalidate(obj,cRepAll,cRepInvCoord);
          
          
          
        }
        FreeP(cnt);
        FreeP(disp);
        FreeP(atm2idx);
      }
}

void SculptFree(CSculpt *I)
{
  VLAFreeP(I->NBList);
  VLAFreeP(I->EXList);
  FreeP(I->NBHash);
  FreeP(I->EXHash);
  ShakerFree(I->Shaker);
  OOFreeP(I);
}

