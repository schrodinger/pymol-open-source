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
  int a,a1,a2,b1,b2;
  BondType *b;
  float *v1,*v2,d;
  CoordSet *cs;

  ShakerReset(I->Shaker);

  if(state<obj->NCSet)
    if(obj->CSet[state])
      if(obj->NBond) {
        cs = obj->CSet[state];
        b=obj->Bond;
        printf("%d\n",obj->NBond); 
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
                ShakerAddCons(I->Shaker,b1,b2,d); 
                /* NOTE: storing atom indices, not coord. ind.! */
              }
          }
      }
  PRINTFD(FB_Sculpt)
    " Sculpt-Debug: I->Shaker->NDistCon %d\n",I->Shaker->NDistCon
    ENDFD;
}
       
void SculptIterateObject(CSculpt *I,ObjectMolecule *obj,int state,int n_cycle)
{
  int n_dst;
  CShaker *shk;
  int a,a1,a2,b1,b2;
  CoordSet *cs;
  ShakerDistCon *sdc;
  float *disp = NULL;
  float *v,*v1,*v2;
  int *atm2idx = NULL;
  int *cnt = NULL;
  float dev;
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
        
        /* next, initialize displacements to zero */
        
        while(n_cycle--) {
          v = disp;
          for(a=0;a<obj->NAtom;a++) {
            *(v++)=0.0;
            *(v++)=0.0;          
          *(v++)=0.0;
          }
          
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
          
          for(a=0;a<obj->NAtom;a++) {
            if(cnt[a]) {
              if(!obj->AtomInfo[a].protected) {
                v1 = disp+3*a;
                sc = 0.1/cnt[a];
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

