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
#ifndef _H_Rep
#define _H_Rep

#include"Base.h"
#include"Ray.h"

#define cRepAll    -1
#define cRepLine    0
#define cRepCyl     1
#define cRepDot     2
#define cRepMesh    3
#define cRepSphere  4
#define cRepRibbon  5

#define cRepCnt  6

#define cRepInvVisib 10
#define cRepInvColor 20
#define cRepInvCoord 30

typedef struct Rep {
  void (*fRender)(struct Rep *I,CRay *ray,Pickable **pick);  
  void (*fFree)(struct Rep* I);
  Pickable *P;
} Rep;

void RepInit(Rep *I);
void RepFree(Rep *I);

#endif
