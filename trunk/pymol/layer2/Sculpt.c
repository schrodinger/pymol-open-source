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

#include"CGO.h"

#ifndef R_SMALL8
#define R_SMALL8 0.00000001
#endif

#define NB_HASH_SIZE 262144
#define EX_HASH_SIZE 65536

#define nb_hash(v) \
(((((int)*(v  ))>> 2)&0x0003F)|\
 ((((int)*(v+1))<< 4)&0x00FC0)|\
 ((((int)*(v+2))<<10)&0x3F000))

#define nb_hash_off_i0(v0i,d) \
  ((((d)+v0i)>> 2)&0x0003F)

#define nb_hash_off_i1(v1i,e) \
 ((((e)+v1i)<< 4)&0x00FC0)

#define nb_hash_off_i2(v2i,f) \
 ((((f)+v2i)<<10)&0x3F000)

#define nb_hash_off(v,d,e,f) \
(((((d)+(int)*(v  ))>> 2)&0x0003F)|\
 ((((e)+(int)*(v+1))<< 4)&0x00FC0)|\
 ((((f)+(int)*(v+2))<<10)&0x3F000))

/* below are empirically optimized */

#define ex_hash_i0(a) \
 (((a)^((a)>>5))&0x00FF)

#define ex_hash_i1(b) \
 ((  ((b)<<5))&0xFF00)

#define ex_hash(a,b) \
(((((a)^((a)>>5)))&0x00FF)|\
 (((    ((b)<<5)))&0xFF00))

#ifdef _PYMOL_INLINE
__inline__
#endif
static float ShakerDoDist(float target,float *v0,float *v1,float *d0to1,float *d1to0,float wt)
{
  float d[3],push[3];
  float len,dev,dev_2,sc,result;

  subtract3f(v0,v1,d);
  len = (float)length3f(d);
  dev = target-len;
  if((result=(float)fabs(dev))>R_SMALL8) {
    dev_2 = wt*dev/2.0F;
    if(len>R_SMALL8) { /* nonoverlapping */
      sc = dev_2/len;
      scale3f(d,sc,push);
      add3f(push,d0to1,d0to1);
      subtract3f(d1to0,push,d1to0);
    } else { /* overlapping, so just push along X */
      d0to1[0]-=dev_2;
      d1to0[0]+=dev_2;
    }
  } else
    result = 0.0;
  return result;
}

#ifdef _PYMOL_INLINE
__inline__
#endif
static float ShakerDoTors(int type,float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,float tole,float wt)
{
  
  float push0[3],push3[3];
  float axis[3],seg0[3],seg1[3],perp0[3],perp1[3];
  float dir[3];
  float sc;
  float sign,dp;
  float result = 0.0F;
  
  /* v0       v3
     \      /
     v1__v2 */
  
  subtract3f(v2,v1,axis);    
  subtract3f(v0,v1,seg0);
  subtract3f(v3,v2,seg1);
  cross_product3f(seg0,axis,perp0);
  cross_product3f(axis,seg1,perp1);
  
  normalize3f(perp0);
  normalize3f(perp1);

#if 0
  {
    float vert[3];
    /* debuggin */
    
    CGOColor(DebugCGO,1.0,0.0,0.0);
    CGOBegin(DebugCGO,GL_LINES);
    CGOVertexv(DebugCGO,v1);
    add3f(perp0,v1,vert);
    CGOVertexv(DebugCGO,vert);
    
    CGOColor(DebugCGO,0.0,1.0,1.0);
    CGOVertexv(DebugCGO,v2);
    add3f(perp1,v2,vert);
    CGOVertexv(DebugCGO,vert);
    CGOEnd(DebugCGO);
  }
#endif

  dp = dot_product3f(perp0,perp1);

  switch(type) {
  case cShakerTorsSP3SP3:
    
    if(dp>-0.5F)
      return 0.0F;
    result = ((float)fabs(dp))-0.5F;
    
    if(result<tole) /* discontinuous low bottom well */
      result = result/25.F;       
    break;
  case cShakerTorsDisulfide:
    if(fabs(dp)<tole)
      return 0.0F;
    result = -dp;
    if(result<tole)
      result = result/25.F;             
    break;
  }


    cross_product3f(perp0,perp1,dir);
    sign = dot_product3f(axis,dir);

  if(sign<0.0F)
    sc=wt*result;
  else
    sc=-wt*result;
  
  scale3f(perp0,sc,push0);
  scale3f(perp1,sc,push3);
  
  add3f(p0,push0,p0);
  add3f(p3,push3,p3);
  subtract3f(p1,push0,p1);
  subtract3f(p2,push3,p2);

#if 0
    {
      float vert[3];
      /* debuggin */
      
      CGOColor(DebugCGO,1.0,0.5,0.0);
      CGOBegin(DebugCGO,GL_LINES);
      CGOLinewidth(DebugCGO,3.0);
      
      /* draw from */
      
      CGOVertexv(DebugCGO,v0);
      add3f(v0,push0,vert);
      CGOVertexv(DebugCGO,vert);
      
      /* draw to */
      
      CGOVertexv(DebugCGO,v3);
      add3f(v3,push3,vert);
      CGOVertexv(DebugCGO,vert);
      CGOEnd(DebugCGO);              
    }    
#endif

  return result;
  
}

CSculpt *SculptNew(PyMOLGlobals *G)
{
  OOAlloc(G,CSculpt);
  I->G=G;
  I->Shaker = ShakerNew(G);
  I->NBList = VLAlloc(int,150000);
  I->NBHash = Alloc(int,NB_HASH_SIZE);
  I->EXList = VLAlloc(int,100000);
  I->EXHash = Alloc(int,EX_HASH_SIZE);
  I->Don = VLAlloc(int,1000);
  I->Acc = VLAlloc(int,1000);
  { 
    int a;
    for(a=1;a<256;a++)
      I->inverse[a] = 1.0F/a;
  }
  return(I);
}

void SculptMeasureObject(CSculpt *I,ObjectMolecule *obj,int state)
{
  PyMOLGlobals *G=I->G;
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

  PRINTFD(G,FB_Sculpt)
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
            ai->sculpt_id=SculptCacheNewID(G);
          ai++;
        }
        
        ObjectMoleculeVerifyChemistry(obj);
        ObjectMoleculeUpdateNeighbors(obj);

        cs = obj->CSet[state];

        use_cache = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_memory);
        if(obj->NBond) {

          planer=Alloc(int,obj->NAtom);
          linear=Alloc(int,obj->NAtom);
          ai = obj->AtomInfo;
          for(a=0;a<obj->NAtom;a++) {
            planer[a]=(ai->geom==cAtomInfoPlaner);
            linear[a]=(ai->geom==cAtomInfoLinear);
            ai++;
          }

          /* brain-dead donor/acceptor assignment
           * REPLACE later on with pattern-based system */


          /* pass 1 */

          b=obj->Bond;
          for(a=0;a<obj->NBond;a++)
            {
              b1 = b->index[0];
              b2 = b->index[1];
              
              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;

              /* make blanket assumption that all nitrogens with 
                 <3 bonds are donors -- we qualify this below...*/
              
              if(ai1->protons==cAN_N) {
                n1 = obj->Neighbor[b1];
                if(obj->Neighbor[n1]<3) { /* N with L.P. */
                  I->Don[b1]=true;
                }
              }

              if(ai2->protons==cAN_N) {
                n2 = obj->Neighbor[b2];
                if(obj->Neighbor[n2]<3) { /* N with L.P. */
                  I->Don[b2]=true;
                }
              }

              /* assume O is always an acceptor...*/
              
              if(ai1->protons==cAN_O) I->Acc[b1]=true;
              if(ai2->protons==cAN_O) I->Acc[b2]=true;
              b++;
            }

          /* pass 2 */
          b=obj->Bond;             
          for(a=0;a<obj->NBond;a++)
            {
              b1 = b->index[0];
              b2 = b->index[1];

              /* nitrogens with lone pairs are acceptors 
                 (not donors as assumed above) */
              
              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;
              
              if(ai1->protons==cAN_N) {
                if(b->order==2) {
                  n1 = obj->Neighbor[b1];
                  if(obj->Neighbor[n1]<3) { /* N with L.P. */
                    I->Acc[b1]=true;
                    I->Don[b1]=false;
                  }
                }
              }
              if(ai2->protons==cAN_N) {
                if(b->order==2) {
                  n2 = obj->Neighbor[b2];
                  if(obj->Neighbor[n2]<3) { /* N with L.P. */
                    I->Acc[b2]=true;
                    I->Don[b2]=false;
                  }
                }
              }
              b++;
            }

          /* pass 3 */
          b=obj->Bond;
          for(a=0;a<obj->NBond;a++)
            {
              b1 = b->index[0];
              b2 = b->index[1];
              
              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;
                     
              /* however, every NH is a donor, 
                 even if it's SP2 */
              
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
              }

              b++;
            }

          /* atom pass */
          ai1 = obj->AtomInfo;
          for(a=0;a<obj->NAtom;a++) {
            /* make sure all nonbonded atoms get categorized */

            n0 = obj->Neighbor[a];
            if(obj->Neighbor[n0]==0) { /* nonbonded */
              if(ai1->protons==cAN_O) {
                I->Don[a] = true;
                I->Acc[a] = true;
              } else if(ai1->protons==cAN_N) {
                I->Don[a] = true;
              } 
            }
            /*            
            if(I->Acc[a]) {
              printf("ACC %s %s %s\n",ai1->chain,ai1->resi,ai1->name);
            }
            if(I->Don[a]) {
              printf("DON %s %s %s\n",ai1->chain,ai1->resi,ai1->name);
              }*/

            ai1++;
          }
          
          /*  exclusions */
          b=obj->Bond;          
          for(a=0;a<obj->NBond;a++)
            {
              b1 = b->index[0];
              b2 = b->index[1];
              
              ai1=obj->AtomInfo+b1;
              ai2=obj->AtomInfo+b2;
          
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
                    if(!SculptCacheQuery(G,cSculptBond,
                                         oai[b1].sculpt_id,
                                         oai[b2].sculpt_id,0,0,&d))
                      SculptCacheStore(G,cSculptBond,
                                       oai[b1].sculpt_id,
                                       oai[b2].sculpt_id,0,0,d);
                  }
                  ShakerAddDistCon(I->Shaker,b1,b2,d,cShakerDistBond); 
                  /* NOTE: storing atom indices, not coord. ind.! */
                }
              b++;
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
                    if(!SculptCacheQuery(G,cSculptAngl,
                                         oai[b0].sculpt_id,
                                         oai[b1].sculpt_id,
                                         oai[b2].sculpt_id,0,&d))
                      SculptCacheStore(G,cSculptAngl,
                                       oai[b0].sculpt_id,
                                       oai[b1].sculpt_id,
                                       oai[b2].sculpt_id,0,d);
                  }
                  ShakerAddDistCon(I->Shaker,b1,b2,d,cShakerDistAngle); 


                  if(linear[b0]&&(linear[b1]||linear[b2])) {
                    
                    if(use_cache) {
                      if(!SculptCacheQuery(G,cSculptLine,
                                           oai[b1].sculpt_id,
                                           oai[b0].sculpt_id,
                                           oai[b2].sculpt_id,0,&dummy))
                        SculptCacheStore(G,cSculptLine,
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
                      if(!SculptCacheQuery(G,cSculptPyra,
                                           oai[b1].sculpt_id,
                                           oai[b0].sculpt_id,
                                           oai[b2].sculpt_id,
                                           oai[b3].sculpt_id,
                                           &d))
                        SculptCacheStore(G,cSculptPyra,
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
                      
                      if(!(planer[b0]||planer[b2]||linear[b0]||linear[b2])) {
                        int type;
                        if((oai[b0].protons == cAN_S)&&
                           (oai[b2].protons == cAN_S))
                          type = cShakerTorsDisulfide;
                        else
                          type = cShakerTorsSP3SP3;
                        ShakerAddTorsCon(I->Shaker,b1,b0,b2,b3,type);
                      }
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
                        if(planer[b0]&&planer[b2])
                          *(j++)=-4;
                        else
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
                          if(!SculptCacheQuery(G,cSculptPlan,
                                               oai[b1].sculpt_id,
                                               oai[b0].sculpt_id,
                                               oai[b2].sculpt_id,
                                               oai[b3].sculpt_id,
                                               &d))
                            SculptCacheStore(G,cSculptPlan,
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

  PRINTFB(G,FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NDistCon %d\n",I->Shaker->NDistCon
    ENDFB(G);
  PRINTFB(G,FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NPyraCon %d\n",I->Shaker->NPyraCon
    ENDFB(G);
  PRINTFB(G,FB_Sculpt,FB_Blather)
    " Sculpt: I->Shaker->NPlanCon %d\n",I->Shaker->NPlanCon
    ENDFB(G);


 PRINTFD(G,FB_Sculpt)
    " SculptMeasureObject-Debug: leaving...\n"
    ENDFD;


}
#ifdef _PYMOL_INLINE
__inline__
#endif
static int SculptCheckBump(float *v1,float *v2,float *diff,float *dist,float cutoff)
{
  register float d2;
  diff[0] = (v1[0]-v2[0]);
  diff[1] = (v1[1]-v2[1]);
  if(fabs(diff[0])>cutoff) return(false);
  diff[2] = (v1[2]-v2[2]);
  if(fabs(diff[1])>cutoff) return(false);
  if(fabs(diff[2])>cutoff) return(false);
  d2 = (diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
  if(d2<(cutoff*cutoff)) {
    *dist = (float)sqrt(d2);
    return(true);
  }
  return(false);
}

#ifdef _PYMOL_INLINE
__inline__
#endif
static int SculptDoBump(float target,float actual,float *d,
                 float *d0to1,float *d1to0,float wt,float *strain)
{
  float push[3];
  float dev,dev_2,sc,abs_dev;

  dev = target-actual;
  if((abs_dev=(float)fabs(dev))>R_SMALL8) {
    dev_2 = wt*dev/2.0F;
    (*strain) += abs_dev;
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


       
float SculptIterateObject(CSculpt *I,ObjectMolecule *obj,int state,int n_cycle)
{
  PyMOLGlobals *G=I->G;
  CShaker *shk;
  int a,aa,a0,a1,a2,a3,b0,b1,b2,b3;
  CoordSet *cs;
  ShakerDistCon *sdc;
  ShakerPyraCon *spc;
  ShakerPlanCon *snc;
  ShakerLineCon *slc;
  ShakerTorsCon *stc;
  float *disp = NULL;
  float *v,*v0,*v1,*v2,*v3;
  float diff[3],len;
  int *atm2idx = NULL;
  int *cnt = NULL;
  int *i,*j;
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
  float tors_wt;
  float tors_tole;
  int active_flag=false;
  float hb_overlap,hb_overlap_base;
  int *active,n_active;
  AtomInfoType *ai0,*ai1;
  double task_time;
  float vdw_magnify;
  int nb_skip,nb_skip_count;
  float total_strain=0.0F;
  int total_count=1;

  PRINTFD(G,FB_Sculpt)
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

        PRINTFD(G,FB_Sculpt)
          " SIO-Debug: NDistCon %d\n",shk->NDistCon
          ENDFD;

        cs = obj->CSet[state];

        nb_skip = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_nb_interval);
        if(nb_skip<1) nb_skip=1;
        vdw = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_scale);
        vdw14 = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_scale14);
        vdw_wt = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_weight);
        vdw_wt14 = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_vdw_weight14);
        bond_wt =  SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_bond_weight);
        angl_wt =  SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_angl_weight);
        pyra_wt =  SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_pyra_weight);
        plan_wt =  SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_plan_weight);
        line_wt =  SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_line_weight);
        mask = SettingGet_i(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_field_mask);
        hb_overlap = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_hb_overlap);
        hb_overlap_base = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_hb_overlap_base);
        tors_tole = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_tors_tolerance);
        tors_wt = SettingGet_f(G,cs->Setting,obj->Obj.Setting,cSetting_sculpt_tors_weight);
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

          task_time = UtilGetSeconds(G);
          vdw_magnify = 1.0F;
          nb_skip_count = n_cycle - nb_skip * (n_cycle/nb_skip);
          if(!nb_skip_count) nb_skip_count = nb_skip;
          while(n_cycle--) {
            
            CGOReset(G->DebugCGO);

            total_strain = 0.0F;
            total_count = 0;
            /* initialize displacements to zero */
        
            v = disp;
            i = cnt;
            for(aa=0;aa<n_active;aa++) {
              a = active[aa];
              v=disp+a*3;
              cnt[a]=0;
              *(v  )=0.0F;
              *(v+1)=0.0F;          
              *(v+2)=0.0F;
            }
          
            /* apply distance constraints */

            sdc=shk->DistCon;
            for(a=0;a<shk->NDistCon;a++) {
              b1 = sdc->at0;
              b2 = sdc->at1;
              
              switch(sdc->type) {
              case cShakerDistBond:
                eval_flag = cSculptBond & mask;
                wt = bond_wt;
                break;
              case cShakerDistAngle:
                eval_flag = cSculptAngl & mask;
                wt = angl_wt;
                break;
              case cShakerDistLimit:
                eval_flag = true;
                wt = 2.0F;
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
                      total_strain+=ShakerDoDist(sdc->targ,v1,v2,disp+b1*3,disp+b2*3,wt);
                      total_count++;
                    } else {
                      total_strain+=ShakerDoDistLimit(sdc->targ,v1,v2,disp+b1*3,disp+b2*3,wt);
                      total_count++;
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
                    total_strain+=ShakerDoLine(v0,v1,v2,disp+b0*3,disp+b1*3,disp+b2*3,line_wt);
                    total_count++;
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
                  total_strain+=ShakerDoPyra(spc->targ,
                               v0,v1,v2,v3,
                               disp+b0*3,
                               disp+b1*3,
                               disp+b2*3,
                               disp+b3*3,
                               pyra_wt);
                  total_count++;
                
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
                  total_strain+=ShakerDoPlan(v0,v1,v2,v3,
                                             disp+b0*3,
                                             disp+b1*3,
                                             disp+b2*3,
                                             disp+b3*3,
                                             plan_wt);
                  total_count++;
                  cnt[b0]++;
                  cnt[b1]++;
                  cnt[b2]++;
                  cnt[b3]++;
                }
              
                snc++;
              }
            }
            
            /* apply torsion constraints */

            if(cSculptTors & mask) {

              /* apply planarity constraints */
            
              stc=shk->TorsCon;
              for(a=0;a<shk->NTorsCon;a++) {
              
                b0 = stc->at0;
                b1 = stc->at1;
                b2 = stc->at2;
                b3 = stc->at3;
                a0 = atm2idx[b0];
                a1 = atm2idx[b1];
                a2 = atm2idx[b2];
                a3 = atm2idx[b3];
              
                if((a0>=0)&&(a1>=0)&&(a2>=0)&&(a3>=0)) {
                  v0 = cs->Coord+3*a0;
                  v1 = cs->Coord+3*a1;
                  v2 = cs->Coord+3*a2;
                  v3 = cs->Coord+3*a3;
                  total_strain+=ShakerDoTors(stc->type,
                                             v0,v1,v2,v3,
                                             disp+b0*3,
                                             disp+b1*3,
                                             disp+b2*3,
                                             disp+b3*3,
                                             tors_tole,
                                             tors_wt);
                  total_count++;
                  cnt[b0]++;
                  cnt[b1]++;
                  cnt[b2]++;
                  cnt[b3]++;
                }
                stc++;
              }
            }

            

            /* apply nonbonded interactions */

            if(nb_skip_count>1) { /* don't do nonbonded each round -- skip and then weight extra */
              nb_skip_count--;
              vdw_magnify+=1.0F;
            } else {
              int nb_off0,nb_off1;
              int v0i,v1i,v2i;
              int x0i;
              int don_b0;
              int acc_b0;
              nb_skip_count = nb_skip;
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
                  don_b0 = I->Don[b0];
                  acc_b0 = I->Acc[b0];
                  v0i = (int)(*v0);
                  v1i = (int)(*(v0+1));
                  v2i = (int)(*(v0+2));
                  x0i = ex_hash_i0(b0);
                  for(h=-4;h<5;h+=4) {
                    nb_off0 = nb_hash_off_i0(v0i,h);
                    for(k=-4;k<5;k+=4) {
                      nb_off1 = nb_off0 | nb_hash_off_i1(v1i,k);
                      for(l=-4;l<5;l+=4) { 
                        {
                          /*  offset = *(I->NBHash+nb_hash_off(v0,h,k,l));*/
                          offset = *(I->NBHash + (nb_off1 | nb_hash_off_i2(v2i,l)));
                          while(offset) {
                            i = I->NBList + offset;
                            b1 = *(i+2);
                            if(b1>b0) { 
                              /* determine exclusion (if any) */
                              xoffset = *(I->EXHash+ (x0i | ex_hash_i1(b1)));
                              ex = 10;
                              while(xoffset) {
                                xoffset = (*(j = I->EXList + xoffset));
                                if((*(j+1)==b0)&&(*(j+2)==b1)) {
                                  ex1 = *(j+3);
                                  if(ex1<ex) {
                                    ex=ex1;
                                  }
                                }
                              }
                              if(ex>3) {
                                ai1=obj->AtomInfo+b1;
                                cutoff = ai0->vdw+ai1->vdw;

                                if(ex==4) { /* 1-4 interation */
                                  cutoff*=vdw14;
                                  wt = vdw_wt14 * vdw_magnify;

                                  if(cSculptVDW14 & mask) {
                                    a1 = atm2idx[b1];
                                    v1 = cs->Coord+3*a1;
                                    if(SculptCheckBump(v0,v1,diff,&len,cutoff))
                                      if(SculptDoBump(cutoff,len,diff,
                                                      disp+b0*3,disp+b1*3,wt,&total_strain)) {
                                        cnt[b0]++;
                                        cnt[b1]++;
                                        total_count++;
                                      }
                                  }
                                } else { /* standard interaction */
                                  if(don_b0&&I->Acc[b1]) { /* h-bond */
                                    if(ai0->protons==cAN_H) {
                                      cutoff-=hb_overlap;
                                    } else {
                                      cutoff-=hb_overlap_base;
                                    }
                                  } else if(acc_b0&&I->Don[b1]) { /* h-bond */
                                    if(ai1->protons==cAN_H) {
                                      cutoff-=hb_overlap;
                                    } else {
                                      cutoff-=hb_overlap_base;
                                    } 
                                  }
                                  cutoff=cutoff*vdw;
                                  wt = vdw_wt * vdw_magnify;
                                  if(cSculptVDW & mask) {
                                    a1 = atm2idx[b1];
                                    v1 = cs->Coord+3*a1;
                                    if(SculptCheckBump(v0,v1,diff,&len,cutoff))
                                      if(SculptDoBump(cutoff,len,diff,
                                                      disp+b0*3,disp+b1*3,wt,&total_strain)) {
                                        cnt[b0]++;
                                        cnt[b1]++;
                                        total_count++;
                                      }
                                  }
                                }
                              }
                            }
                            offset=(*i);
                          }
                        }
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
              vdw_magnify = 1.0F;
            }
            /* average the displacements */
               
            {
              int cnt_a;
              float _1 = 1.0F;
              register float inv_cnt;
              int *a_ptr = active;
              register float *lookup_inverse = I->inverse;
              for(aa=0;aa<n_active;aa++) {
                if( (cnt_a = cnt[(a = *(a_ptr++))]) ) {
                  if(!obj->AtomInfo[a].protekted) {
                    v1 = disp+3*a;
                    v2 = cs->Coord+3*atm2idx[a];
                    if(!(cnt_a&0xFFFFFF00)) /* don't divide -- too slow */
                      inv_cnt = lookup_inverse[cnt_a];
                    else
                      inv_cnt = _1/cnt_a;
                    *(v2  )+=(*(v1  ))*inv_cnt;
                    *(v2+1)+=(*(v1+1))*inv_cnt;
                    *(v2+2)+=(*(v1+2))*inv_cnt;
                  }
                }
              }
            }
            if(cs->fInvalidateRep) {
              cs->fInvalidateRep(cs,cRepAll,cRepInvCoord);
            } else {
              ObjectMoleculeInvalidate(obj,cRepAll,cRepInvCoord);
            }
          
          }
          
          task_time = UtilGetSeconds(G) - task_time;
          PRINTFB(G,FB_Sculpt,FB_Blather)
            " Sculpt: %2.5f seconds %8.3f %d %8.3f\n",task_time,total_strain,total_count,
            100*total_strain/total_count
            ENDFB(G);

          if(total_count) 
            total_strain = (1000*total_strain)/total_count;
        }
        FreeP(active);
        FreeP(cnt);
        FreeP(disp);
        FreeP(atm2idx);
      }

  PRINTFD(G,FB_Sculpt)
    " SculptIterateObject-Debug: leaving...\n"
    ENDFD;

  CGOStop(G->DebugCGO);
  

  return total_strain;
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

#if 0
 

                                    if(ex1==4) { /* special handling for 1-4's */
                                      ex14_b0=*(j+4); /* b1\b0_b2/b3 atoms */
                                      ex14_b2=*(j+5);
                                    }

                                      if(SculptDoBump14(cutoff,len,diff,
                                                        disp+b0*3,disp+b1*3,
                                                        wt,&total_strain,
                                                        v0,
                                                        v1,
                                                        cs->Coord+3*atm2idx[ex14_b0],
                                                        cs->Coord+3*atm2idx[ex14_b2],
                                                        disp+ex14_b0*3,
                                                        disp+ex14_b2*3)) {
                                        cnt[b0]++;
                                        cnt[b1]++;
                                        total_count++;
                                      }
                                  }

#ifdef _PYMOL_INLINE
__inline__
#endif
static int SculptDoBump14(float target,float actual,float *d,
                          float *d0to1,float *d1to0,float wt,float *strain,
                          float *v0,float *v1,float *v2,float *v3,
                          float *d2to3,float *d3to2)
{
  /* v0\v2_v3/v1 */

  float push0[3],push1[3],balance[3];
  float axis[3],seg0[3],seg1[3],perp0[3],perp1[3];
  float dir[3];
  float dev,dev_2,sc,abs_dev;
  float sign,dp;

  dev = target-actual;
  if((abs_dev=fabs(dev))>R_SMALL8) {

    /* v0        v1
         \      /
         v2__v3 */

    subtract3f(v3,v2,axis);    
    subtract3f(v0,v2,seg0);
    normalize3f(axis);
    subtract3f(v1,v3,seg1);
    cross_product3f(seg0,axis,perp0);
    cross_product3f(axis,seg1,perp1);

    normalize3f(perp0);
    normalize3f(perp1);

#if 0
    {
      float vert[3];
    /* debuggin */

      CGOColor(DebugCGO,1.0,0.0,0.0);
      CGOBegin(DebugCGO,GL_LINES);
      CGOVertexv(DebugCGO,v2);
      add3f(perp0,v2,vert);
      CGOVertexv(DebugCGO,vert);
      
      CGOColor(DebugCGO,0.0,1.0,1.0);
      CGOVertexv(DebugCGO,v3);
      add3f(perp1,v3,vert);
      CGOVertexv(DebugCGO,vert);
      CGOEnd(DebugCGO);
    }
#endif
    
    if((dp = dot_product3f(perp0,perp1))>-0.5F)
      return 0;
    
    
    dev_2 = wt*dev/2.0F;
    (*strain) += abs_dev;
    if(actual>R_SMALL8) { /* nonoverlapping */

      cross_product3f(perp0,perp1,dir);
      sign = dot_product3f(axis,dir);
      
      if(sign<0.0F)
        sc=wt*(fabs(dp)-0.5);
      else
        sc=wt*(0.5-fabs(dp));
      
      scale3f(perp0,sc,push0);
      scale3f(perp1,sc,push1);
      add3f(push0,push1,balance);
      scale3f(balance,-0.5,balance);

#if 0
      {
        float vert[3];
      /* debuggin */

        CGOColor(DebugCGO,1.0,0.5,0.0);
        CGOBegin(DebugCGO,GL_LINES);
        CGOLinewidth(DebugCGO,3.0);
        
        /* draw from */
        
        CGOVertexv(DebugCGO,v0);
        add3f(v0,push0,vert);
        CGOVertexv(DebugCGO,vert);
        
        /* draw to */
        
        CGOVertexv(DebugCGO,v1);
        add3f(v1,push1,vert);
        CGOVertexv(DebugCGO,vert);

        CGOVertexv(DebugCGO,v2);
        add3f(v2,balance,vert);
        CGOVertexv(DebugCGO,vert);

        CGOVertexv(DebugCGO,v3);
        add3f(v3,balance,vert);
        CGOVertexv(DebugCGO,vert);
        
        CGOEnd(DebugCGO);        
        
      }
#endif

      add3f(d0to1,push0,d0to1);
      add3f(d1to0,push1,d1to0);
      subtract3f(d2to3,push0,d2to3);
      subtract3f(d3to2,push1,d3to2);

    } else { /* overlapping, so just push along X */
      d0to1[0]-=dev_2;
      d1to0[0]+=dev_2;
    }
    return 1;
  }
  return 0;
}

#endif


#if 0
 { int a4,b4,n3;
  float *v4;

          /* add additional long-range distance restraints into object (INEFFICIENT) */

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
 }
#endif
