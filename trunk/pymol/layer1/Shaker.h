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
#ifndef _H_Shaker
#define _H_Shaker

typedef struct {
  int at0,at1;
  float targ;
} ShakerDistCon;

typedef struct {
  ShakerDistCon *DistCon;
  int NDistCon;
} CShaker;

CShaker *ShakerNew(void);
void ShakerReset(CShaker *I);
void ShakerAddCons(CShaker *I,int atom0,int atom1,float dist);

float ShakerDoDist(float target,float *v0,float *v1,float *d0to1,float *d1to0);

void ShakerFree(CShaker *I);

#endif


