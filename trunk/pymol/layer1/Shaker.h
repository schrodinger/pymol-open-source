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
  int type;
  float targ;
} ShakerDistCon;

typedef struct {
  int at0,at1,at2,at3;
  float targ;
} ShakerPyraCon;

typedef struct {
  int at0,at1,at2,at3;
} ShakerPlanCon;

typedef struct {
  ShakerDistCon *DistCon;
  int NDistCon;
  ShakerPyraCon *PyraCon;
  int NPyraCon;
  ShakerPlanCon *PlanCon;
  int NPlanCon;
} CShaker;

CShaker *ShakerNew(void);
void ShakerReset(CShaker *I);
void ShakerAddDistCon(CShaker *I,int atom0,int atom1,float dist,int type);

void ShakerAddPyraCon(CShaker *I,int atom0,int atom1,int atom2,int atom3,float target);
void ShakerAddPlanCon(CShaker *I,int atom0,int atom1,int atom2,int atom3);

float ShakerGetPyra(float *v0,float *v1,float *v2,float *v3);

float ShakerDoDist(float target,float *v0,float *v1,float *d0to1,float *d1to0,float wt);
float ShakerDoPyra(float target,float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,float wt);

float ShakerDoPlan(float *v0,float *v1,float *v2,float *v3,
                   float *p0,float *p1,float *p2,float *p3,float wt);

void ShakerFree(CShaker *I);

#endif


