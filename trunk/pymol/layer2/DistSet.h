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
#ifndef _H_DistSet
#define _H_DistSet

#include"Rep.h"
#include"Setting.h"

typedef struct DistSet {
  void (*fUpdate)(struct DistSet *I);
  void (*fRender)(struct DistSet *I,CRay *ray,Pickable **pick,int pass);
  void (*fFree)(struct DistSet *I);
  void (*fInvalidateRep)(struct DistSet *I,int type,int level);
  struct ObjectDist *Obj;
  float *Coord;
  int NIndex;
  Rep **Rep; /* an array of pointers to representations */
  int NRep;
  CSetting *Setting;
} DistSet;

#include"ObjectDist.h"

DistSet *DistSetNew(void);

int DistSetGetExtent(DistSet *I,float *mn,float *mx);

#endif

