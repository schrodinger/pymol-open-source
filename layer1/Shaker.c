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

#include"Map.h"

#include"Shaker.h"

#ifndef R_SMALL8
#define R_SMALL8 0.00000001
#endif

CShaker *ShakerNew(void)
{
  OOAlloc(CShaker);
  I->DistCon = VLAlloc(ShakerDistCon,1000);
  I->PyraCon = VLAlloc(ShakerPyraCon,1000);
  I->PlanCon = VLAlloc(ShakerPlanCon,1000);
  I->NDistCon = 0;
  I->NPyraCon = 0;
  I->NPlanCon = 0;
  return(I);
}

void ShakerReset(CShaker *I)
{
  I->NDistCon = 0;
}

float ShakerDoDist(float target,float *v0,float *v1,float *d0to1,float *d1to0)
{
  float d[3],push[3];
  float len,dev,dev_2,sc;

  subtract3f(v0,v1,d);
  len = length3f(d);
  dev = target-len;
  if(fabs(dev)>R_SMALL8) {
    dev_2 = dev/2.0;
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

void ShakerAddDistCon(CShaker *I,int atom0,int atom1,float target)
{
  ShakerDistCon *sdc;

  VLACheck(I->DistCon,ShakerDistCon,I->NDistCon);
  sdc = I->DistCon+I->NDistCon;
  sdc->at0=atom0;
  sdc->at1=atom1;
  sdc->targ = target;
  I->NDistCon++;
}

float ShakerGetPyra(float *v0,float *v1,float *v2,float *v3)
{
  float d0[0],cp[3],d2[3],d3[3];
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
                   float *p0,float *p1,float *p2,float *p3)
{
  float d0[0],cp[3],d2[3],d3[3],push[3];
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
    sc =  dev;
    scale3f(cp,sc,push);
    add3f(push,p0,p0);
    scale3f(push,1.0/3,push);
    subtract3f(p1,push,p1);
    subtract3f(p2,push,p2);
    subtract3f(p3,push,p3);
  } else
    dev = 0.0;
  return dev;

}

float ShakerDoPlan(float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3)
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

  dev = fabs(cur);
  if(fabs(dev)>R_SMALL8) {

    sc = -dev/2.0;

    subtract3f(v0,v3,d0);
    normalize3f(d0);
    scale3f(d0,sc,push);
    add3f(push,p0,p0); 
    subtract3f(p3,push,p3);

    sc = -2.0*sc ;

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

void ShakerFree(CShaker *I)
{
  VLAFreeP(I->PlanCon);
  VLAFreeP(I->PyraCon);
  VLAFreeP(I->DistCon);
  OOFreeP(I);
}





