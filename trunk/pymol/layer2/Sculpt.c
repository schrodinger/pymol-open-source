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

#include"Map.h"

#include"Sculpt.h"

#ifndef R_SMALL8
#define R_SMALL8 0.00000001
#endif

CSculpt *SculptNew(void)
{
  OOAlloc(CSculpt);
  I->Shaker = ShakerNew();
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
  AtomInfoType *ai;

  ShakerReset(I->Shaker);

  if(state<obj->NCSet)
    if(obj->CSet[state])
      if(obj->NBond) {

        planer=Alloc(int,obj->NAtom);
        ai = obj->AtomInfo;
        for(a=0;a<obj->NAtom;a++) {
          planer[a]=(ai->geom==3);
          ai++;
        }

        cs = obj->CSet[state];
        b=obj->Bond;
        for(a=0;a<obj->NBond;a++)
          {
            b1 = b->index[0];
            b2 = b->index[1];
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
                ShakerAddDistCon(I->Shaker,b1,b2,d); 
                /* NOTE: storing atom indices, not coord. ind.! */
              }
          }

        /* now pick up those 1-3 interations */

        ObjectMoleculeVerifyChemistry(obj);
        ObjectMoleculeUpdateNeighbors(obj);

        for(b0=0;b0<obj->NAtom;b0++) {
          n0 = obj->Neighbor[b0]+1;
          while(obj->Neighbor[n0]>=0) {
            b1 = obj->Neighbor[n0];
            n1 = n0+2;
            while(obj->Neighbor[n1]>=0) {
              b2 = obj->Neighbor[n1];
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
                ShakerAddDistCon(I->Shaker,b1,b2,d); 
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
                  if((cs==obj->DiscreteCSet[b0])&&(cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
                    a0=obj->DiscreteAtmToIdx[b0];
                    a1=obj->DiscreteAtmToIdx[b1];
                    a2=obj->DiscreteAtmToIdx[b2];
                    a2=obj->DiscreteAtmToIdx[b3];
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
                  
                  if(obj->DiscreteFlag) {
                    if((cs==obj->DiscreteCSet[b0])&&(cs==obj->DiscreteCSet[b1])&&(cs==obj->DiscreteCSet[b2])) {
                      a0=obj->DiscreteAtmToIdx[b0];
                      a1=obj->DiscreteAtmToIdx[b1];
                      a2=obj->DiscreteAtmToIdx[b2];
                      a2=obj->DiscreteAtmToIdx[b3];
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
                      if(planer[b0]&&planer[b2])
                        ShakerAddPlanCon(I->Shaker,b1,b0,b2,b3); 
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
  
  PRINTFD(FB_Sculpt)
    " Sculpt-Debug: I->Shaker->NDistCon %d\n",I->Shaker->NDistCon
    ENDFD;
  PRINTFD(FB_Sculpt)
    " Sculpt-Debug: I->Shaker->NPyraCon %d\n",I->Shaker->NPyraCon
    ENDFD;
  PRINTFD(FB_Sculpt)
    " Sculpt-Debug: I->Shaker->NPlanCon %d\n",I->Shaker->NPlanCon
    ENDFD;
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
  int *atm2idx = NULL;
  int *cnt = NULL;
  float sc;

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

        for(a=0;a<obj->NAtom;a++) {
          atm2idx[a]=-1;
          cnt[a]=0;
        }
        
        /* first, create coordinate -> vertex mapping */
        /* and count number of constraints */

        sdc=shk->DistCon;
        for(a=0;a<shk->NDistCon;a++) {
          b1 = sdc->at0;
          b2 = sdc->at1;
          
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
          atm2idx[b1] = a1;
          atm2idx[b2] = a2;
          cnt[b1]++;
          cnt[b2]++;
          sdc++;
        }

        spc=shk->PyraCon;
        for(a=0;a<shk->NPyraCon;a++) {
          cnt[spc->at0]++; 
          cnt[spc->at1]++;
          cnt[spc->at2]++;
          cnt[spc->at3]++;
          spc++;
        }
        
        while(n_cycle--) {

          /* initialize displacements to zero */
        
          v = disp;
          for(a=0;a<obj->NAtom;a++) {
            *(v++)=0.0;
            *(v++)=0.0;          
            *(v++)=0.0;
          }
          
          /* apply distance constraints */

          sdc=shk->DistCon;
          for(a=0;a<shk->NDistCon;a++) {
            b1 = sdc->at0;
            b2 = sdc->at1;
            a1 = atm2idx[b1]; /* coordinate set indices */
            a2 = atm2idx[b2];
            
            if((a1>=0)&&(a2>=0))
              {
                v1 = cs->Coord+3*a1;
                v2 = cs->Coord+3*a2;
                ShakerDoDist(sdc->targ,v1,v2,disp+b1*3,disp+b2*3);
              }
            sdc++;
          }

          /* apply pyramid constraints */

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
                           disp+b3*3);
            }

            spc++;
          }

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
                           disp+b3*3);
            }

            snc++;
          }

          /* average the displacements */

          for(a=0;a<obj->NAtom;a++) {
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
        }
        ObjectMoleculeInvalidate(obj,cRepAll,cRepInvCoord);
        
        FreeP(cnt);
        FreeP(disp);
        FreeP(atm2idx);
      }
}

void SculptFree(CSculpt *I)
{
  ShakerFree(I->Shaker);
  OOFreeP(I);
}

