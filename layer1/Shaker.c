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

void ShakerAddCons(CShaker *I,int atom0,int atom1,float target)
{
  ShakerDistCon *sdc;

  VLACheck(I->DistCon,ShakerDistCon,I->NDistCon);
  sdc = I->DistCon+I->NDistCon;
  sdc->at0=atom0;
  sdc->at1=atom1;
  sdc->targ = target;
  I->NDistCon++;
}

void ShakerFree(CShaker *I)
{
  VLAFreeP(I->DistCon);
  OOFreeP(I);
}

