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
#ifndef _H_Sphere
#define _H_Sphere

#include"Vector.h"

typedef struct DotRec {
  Vector3f v;
  float area;
} DotRec;

typedef struct SphereRec {
  DotRec *dot;
  int *StripLen;
  int *Sequence;
  int NStrip,NVertTot;
  int nDot;
} SphereRec;

extern SphereRec *Sphere0;
extern SphereRec *Sphere1;
extern SphereRec *Sphere2;
extern SphereRec *Sphere3;


void SphereInit(void);
void SphereDone(void);

SphereRec *MakeDotSphere(int level);
void SphereFree(SphereRec *sp);

#endif
