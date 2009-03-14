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
#include"Base.h"
#include"OOMac.h"
#include"CGO.h"

#include"Map.h"

#include"Shaker.h"

CShaker *ShakerNew(PyMOLGlobals *G)
{
  OOAlloc(G,CShaker);
  I->G=G;
  I->DistCon = VLAlloc(ShakerDistCon,1000);
  I->PyraCon = VLAlloc(ShakerPyraCon,1000);
  I->PlanCon = VLAlloc(ShakerPlanCon,1000);
  I->TorsCon = VLAlloc(ShakerTorsCon,1000);
  I->LineCon = VLAlloc(ShakerLineCon,100);
  I->NDistCon = 0;
  I->NPyraCon = 0;
  I->NPlanCon = 0;
  I->NLineCon = 0;
  I->NTorsCon = 0;
  return(I);
}

void ShakerReset(CShaker *I)
{
  I->NDistCon = 0;
  I->NPyraCon = 0;
  I->NPlanCon = 0;
  I->NLineCon = 0;
  I->NTorsCon = 0;
}

#if 0
/* this code moved into Sculpt.c */

float ShakerDoDist(float target,float *v0,float *v1,float *d0to1,float *d1to0,float wt)
{
  float d[3],push[3];
  float len,dev,dev_2,sc,result;

  subtract3f(v0,v1,d);
  len = (float)length3f(d);
  dev = target-len;
  if((result=fabs(dev))>R_SMALL8) {
    dev_2 = wt*dev*0.5F;
    if(len>R_SMALL8) { /* nonoverlapping */
      sc = dev_2/len;
      scale3f(d,sc,push);
      add3f(push,d0to1,d0to1);
      subtract3f(d1to0,push,d1to0);
    } else { /* overlapping, so just push in a random direction */
      float rd[3];
      get_random3f(rd);
      d0to1[0]-=rd[0]*dev_2;
      d1to0[0]+=rd[0]*dev_2;
      d0to1[1]-=rd[1]*dev_2;
      d1to0[1]+=rd[1]*dev_2;
      d0to1[2]-=rd[2]*dev_2;
      d1to0[2]+=rd[2]*dev_2;
    }
  } else
    result = 0.0;
  return result;
}

float ShakerDoDistLimit(float target,float *v0,float *v1,float *d0to1,float *d1to0,float wt)
{
  float d[3],push[3];
  float len,dev,dev_2,sc;

  if(wt==0.0F) return 0.0F;
  subtract3f(v0,v1,d);
  len = (float)length3f(d);
  dev = target-len;
  if(dev<0.0F) { /* assuming len is non-zero since it is above target */
    dev_2 = wt*dev*0.5F;
    sc = dev_2/len;
    scale3f(d,sc,push);
    add3f(push,d0to1,d0to1);
    subtract3f(d1to0,push,d1to0);
    return -dev;
  } else
    return 0.0F;
}

#endif

void ShakerAddDistCon(CShaker *I,int atom0,int atom1,float target,int type,float wt)
{
  ShakerDistCon *sdc;

  VLACheck(I->DistCon,ShakerDistCon,I->NDistCon);
  sdc = I->DistCon+I->NDistCon;
  sdc->at0=atom0;
  sdc->at1=atom1;
  sdc->targ = target;
  sdc->type = type;
  sdc->weight = wt;
  I->NDistCon++;
}

float ShakerGetPyra(float *targ2, float *v0,float *v1,float *v2,float *v3)
{
  float d0[3],cp[3],d2[3],d3[3];
  float av[3],t0[3];

  add3f(v1,v2,av);
  subtract3f(v2,v1,d2);
  add3f(v3,av,av);
  subtract3f(v3,v1,d3);
  subtract3f(av,v0,t0);
  cross_product3f(d2,d3,cp);
  scale3f(av,0.33333333F,av);
  normalize3f(cp);
  subtract3f(av,v0,d0);

  (*targ2) = length3f(d0);
  return(dot_product3f(d0,cp));
}

float ShakerDoPyra(float targ1,float targ2, 
                   float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,
                   float wt, float inv_wt)
{
  float d0[3],cp[3],d2[3],d3[3];
  float av[3],t0[3],push[3];

  float cur,dev,sc,result1,result2 = 0.0F;

  add3f(v1,v2,av);
  subtract3f(v2,v1,d2);
  add3f(v3,av,av);
  subtract3f(v3,v1,d3);
  subtract3f(av,v0,t0);
  cross_product3f(d2,d3,cp);
  scale3f(av,0.33333333F,av);
  normalize3f(cp);
  subtract3f(av,v0,d0);

  cur = dot_product3f(d0,cp);
  dev = cur-targ1;
  result1 = (float)fabs(dev);
  if(result1>R_SMALL8) {
    sc = wt*dev;
    if((cur*targ1)<0.0) /* inverted */
      sc = sc * inv_wt; /* inversion fixing weight */
    scale3f(cp,sc,push);
    add3f(push,p0,p0);
    scale3f(push,0.333333F,push);
    subtract3f(p1,push,p1);
    subtract3f(p2,push,p2);
    subtract3f(p3,push,p3);
  }
  
  if((targ2>=0.0F) && ((cur*targ1>0.0)||(fabs(targ1)<0.1))) {
    /* so long as we're not inverted...
       also make sure v0 is the right distance from the average point */
    cur = length3f(d0);
    normalize3f(d0);
    dev = cur-targ2;
    result2 = (float)fabs(dev);
    if(result2>R_SMALL4) {
      sc = wt*dev*2.0F;
      scale3f(d0,sc,push);
      add3f(push,p0,p0);
      scale3f(push,0.333333F,push);
      subtract3f(p1,push,p1);
      subtract3f(p2,push,p2);
      subtract3f(p3,push,p3);
    }
  }

  return result1+result2;
}


float ShakerDoLine(float *v0,float *v1,float *v2,
                   float *p0,float *p1,float *p2,float wt)
{
  /* v0-v1-v2 */

  float d0[3],d1[3],cp[3],d2[3],d3[3],d4[3],push[3];
  float dev,sc,lcp,result;

  subtract3f(v2,v1,d2);
  subtract3f(v0,v1,d1);
  normalize3f(d2);
  normalize23f(d1,d0);

  cross_product3f(d2,d0,cp); 
  lcp = (float)length3f(cp);
  if(lcp>R_SMALL4) {
    lcp = 1.0F/lcp;
    scale3f(cp,lcp,cp); /* axis 0 */

    subtract3f(v2,v0,d3);
    normalize3f(d3);  /* axis 1 */
    
    cross_product3f(cp,d3,d4); 
    normalize3f(d4); /* displacement direction */

    dev = dot_product3f(d1,d4); /* current deviation */

    if((result = (float)fabs(dev))>R_SMALL8) {
      sc = wt*dev;
      scale3f(d4,sc,push);
      add3f(push,p1,p1);
      scale3f(push,0.5F,push);
      subtract3f(p0,push,p0);
      subtract3f(p2,push,p2);
    } else {
      result = 0.0;
    }
  } else
    result = 0.0;
  return result;

}



float ShakerDoPlan(float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,
                   float target, int fixed, float wt)
{
  
  float result;
  
  float d01[3],d12[3],d23[3],d03[3],cp0[3],cp1[3],dp,sc,dev,d0[3],push[3];
  double s01,s12,s23,s03;
  
  subtract3f(v0,v1,d01);
  subtract3f(v1,v2,d12);
  subtract3f(v2,v3,d23);
  subtract3f(v0,v3,d03);    
  
  s03 = lengthsq3f(d03);
  s01 = lengthsq3f(d01);
  s12 = lengthsq3f(d12);
  s23 = lengthsq3f(d23);
  
  if( (s03<s01) || (s03<s12) || (s03<s23))
    return 0.0F;
  
  cross_product3f(d01,d12,cp0);
  cross_product3f(d12,d23,cp1);
  
  normalize3f(cp0);
  normalize3f(cp1);
  
  dp = dot_product3f(cp0,cp1);
  
  result = (dev = 1.0F - (float)fabs(dp));
  
  if(dev>R_SMALL4) {
    
    /*
      add3f(cp0,cp1,d0);
      normalize3f(d0);
      
      cross_product3f(cp0,d12,pos);
      dp2 = dot_product3f(cp1,pos);
    */
    
    if(fixed && (dp*target<0.0F)) {

      /* fixed & backwards... */

      if(dp<0.0F) {
	sc = -wt*dev*0.5F;
      } else {
	sc = wt*dev*0.5F;
      }
      sc *= 0.02F; /* weaken considerably to allow resolution of
		      inconsistencies (folded rings, etc.) */

    } else if(dp>0) {
      sc = -wt*dev*0.5F;
    } else {
      sc = wt*dev*0.5F;
    }
    
    if(fixed && (fixed<7)) {
      /* in small rings, ramp up the planarity factor */
      sc *= 8;
    } else {
      sc *= 0.2F;
    }
    
    /* pair-wise nudges */
    
    subtract3f(v0,v3,d0);
    normalize3f(d0);
    scale3f(d0,sc,push);
    add3f(push,p0,p0); 
    subtract3f(p3,push,p3);
    
    subtract3f(v1,v2,d0);
    normalize3f(d0);
    scale3f(d0,sc,push);
    add3f(push,p1,p1); 
    subtract3f(p2,push,p2);
    
    sc = -sc;
    subtract3f(v0,v2,d0);
    normalize3f(d0);
    scale3f(d0,sc,push);
    add3f(push,p0,p0); 
    subtract3f(p2,push,p2);
    
    subtract3f(v1,v3,d0);
    normalize3f(d0);
    scale3f(d0,sc,push);
    add3f(push,p1,p1); 
    subtract3f(p3,push,p3);
    
  } else {
    result = 0.0;
  }
  return result;

}

void ShakerAddPyraCon(CShaker *I,int atom0,int atom1,int atom2,int atom3,float targ1,float targ2)
{
  ShakerPyraCon *spc;
  
  VLACheck(I->PyraCon,ShakerPyraCon,I->NPyraCon);
  spc = I->PyraCon+I->NPyraCon;
  spc->at0=atom0;
  spc->at1=atom1;
  spc->at2=atom2;
  spc->at3=atom3;
  spc->targ1 = targ1;
  spc->targ2 = targ2;
  I->NPyraCon++;
}


void ShakerAddTorsCon(CShaker *I,int atom0,int atom1,int atom2,int atom3,int type)
{
  ShakerTorsCon *stc;
  
  VLACheck(I->TorsCon,ShakerTorsCon,I->NTorsCon);
  stc = I->TorsCon+I->NTorsCon;
  stc->at0=atom0;
  stc->at1=atom1;
  stc->at2=atom2;
  stc->at3=atom3;
  stc->type = type;
  I->NTorsCon++;

}

void ShakerAddPlanCon(CShaker *I,int atom0,int atom1,int atom2,int atom3,float target, int fixed)
{
  ShakerPlanCon *spc;
  
  VLACheck(I->PlanCon,ShakerPlanCon,I->NPlanCon);
  spc = I->PlanCon+I->NPlanCon;
  spc->at0=atom0;
  spc->at1=atom1;
  spc->at2=atom2;
  spc->at3=atom3;
  spc->fixed=fixed;
  spc->target=target;
  I->NPlanCon++;
}

void ShakerAddLineCon(CShaker *I,int atom0,int atom1,int atom2)
{
  ShakerLineCon *slc;
  
  VLACheck(I->LineCon,ShakerLineCon,I->NLineCon);
  slc = I->LineCon+I->NLineCon;
  slc->at0=atom0;
  slc->at1=atom1;
  slc->at2=atom2;
  I->NLineCon++;
}

void ShakerFree(CShaker *I)
{
  VLAFreeP(I->PlanCon);
  VLAFreeP(I->PyraCon);
  VLAFreeP(I->DistCon);
  VLAFreeP(I->LineCon);
  VLAFreeP(I->TorsCon);
  OOFreeP(I);
}





