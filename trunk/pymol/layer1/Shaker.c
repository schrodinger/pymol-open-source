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

#include"Map.h"

#include"Shaker.h"

CShaker *ShakerNew(void)
{
  OOAlloc(CShaker);
  I->DistCon = VLAlloc(ShakerDistCon,1000);
  I->PyraCon = VLAlloc(ShakerPyraCon,1000);
  I->PlanCon = VLAlloc(ShakerPlanCon,1000);
  I->LineCon = VLAlloc(ShakerLineCon,100);
  I->NDistCon = 0;
  I->NPyraCon = 0;
  I->NPlanCon = 0;
  I->NLineCon = 0;
  return(I);
}

void ShakerReset(CShaker *I)
{
  I->NDistCon = 0;
}

float ShakerDoDist(float target,float *v0,float *v1,float *d0to1,float *d1to0,float wt)
{
  float d[3],push[3];
  float len,dev,dev_2,sc;

  subtract3f(v0,v1,d);
  len = (float)length3f(d);
  dev = target-len;
  if(fabs(dev)>R_SMALL8) {
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
    dev = 0.0;
  return dev;
}

void ShakerAddDistCon(CShaker *I,int atom0,int atom1,float target,int type)
{
  ShakerDistCon *sdc;

  VLACheck(I->DistCon,ShakerDistCon,I->NDistCon);
  sdc = I->DistCon+I->NDistCon;
  sdc->at0=atom0;
  sdc->at1=atom1;
  sdc->targ = target;
  sdc->type = type;
  I->NDistCon++;
}

float ShakerGetPyra(float *v0,float *v1,float *v2,float *v3)
{
  float d0[3],cp[3],d2[3],d3[3];
  subtract3f(v2,v1,d2);
  normalize3f(d2);
  subtract3f(v3,v1,d3);
  normalize3f(d3);
  cross_product3f(d2,d3,cp);
  normalize3f(cp);
  subtract3f(v1,v0,d0);
  return(dot_product3f(d0,cp));
}

float ShakerDoPyra(float target,float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,float wt)
{
  float d0[3],cp[3],d2[3],d3[3],push[3];
  float cur,dev,sc;
  subtract3f(v2,v1,d2);
  normalize3f(d2);
  subtract3f(v3,v1,d3);
  normalize3f(d3);
  cross_product3f(d2,d3,cp); 
  normalize3f(cp); /* this is our axis */
  subtract3f(v1,v0,d0);
  cur = dot_product3f(d0,cp);

  dev = cur-target;
  if(fabs(dev)>R_SMALL8) {
    sc = wt*dev;
    scale3f(cp,sc,push);
    add3f(push,p0,p0);
    scale3f(push,1.0F/3,push);
    subtract3f(p1,push,p1);
    subtract3f(p2,push,p2);
    subtract3f(p3,push,p3);
  } else
    dev = 0.0;
  return dev;

}


float ShakerDoLine(float *v0,float *v1,float *v2,
                   float *p0,float *p1,float *p2,float wt)
{
  /* v0-v1-v2 */

  float d0[3],d1[3],cp[3],d2[3],d3[3],d4[3],push[3];
  float dev,sc,lcp;

  subtract3f(v2,v1,d2);
  normalize3f(d2);
  subtract3f(v0,v1,d1);
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

    if(fabs(dev)>R_SMALL8) {
      sc = wt*dev;
      scale3f(d4,sc,push);
      add3f(push,p1,p1);
      scale3f(push,0.5F,push);
      subtract3f(p0,push,p0);
      subtract3f(p2,push,p2);
    } else {
      dev = 0.0;
    }
  } else
    dev = 0.0;
  return dev;

}



float ShakerDoPlan(float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,float wt)
{
  float vc[3],d0[3],d1[3],d2[3],cp[3];
  float push[3];
  float cur,dev,sc;

  average3f(v0,v3,vc);

  subtract3f(v1,vc,d1);
  normalize3f(d1);

  subtract3f(v2,vc,d2);
  normalize3f(d2);

  cross_product3f(d1,d2,cp); 
  normalize3f(cp); /* this is our axis */

  subtract3f(v0,vc,d0);
  cur = dot_product3f(d0,cp);

  dev = (float)fabs(cur);
  if(fabs(dev)>R_SMALL8) {

    sc = -wt*dev/2.0F;

    subtract3f(v0,v3,d0);
    normalize3f(d0);
    scale3f(d0,sc,push);
    add3f(push,p0,p0); 
    subtract3f(p3,push,p3);

    sc = -2.0F*sc ;

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

    
  } else
    dev = 0.0;
  return dev;

}

void ShakerAddPyraCon(CShaker *I,int atom0,int atom1,int atom2,int atom3,float target)
{
  ShakerPyraCon *spc;
  
  VLACheck(I->PyraCon,ShakerPyraCon,I->NPyraCon);
  spc = I->PyraCon+I->NPyraCon;
  spc->at0=atom0;
  spc->at1=atom1;
  spc->at2=atom2;
  spc->at3=atom3;
  spc->targ = target;
  I->NPyraCon++;

}

void ShakerAddPlanCon(CShaker *I,int atom0,int atom1,int atom2,int atom3)
{
  ShakerPlanCon *spc;
  
  VLACheck(I->PlanCon,ShakerPlanCon,I->NPlanCon);
  spc = I->PlanCon+I->NPlanCon;
  spc->at0=atom0;
  spc->at1=atom1;
  spc->at2=atom2;
  spc->at3=atom3;
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
  OOFreeP(I);
}





