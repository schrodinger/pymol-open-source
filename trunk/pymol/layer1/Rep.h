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

#define cRepAll       -1
#define cRepLine       0
#define cRepCyl        1
#define cRepDot        2
#define cRepMesh       3
#define cRepSphere     4
#define cRepRibbon     5
#define cRepSurface    6
#define cRepDash       7
#define cRepLabel      8
#define cRepNonbonded  9
#define cRepNonbondedSphere 10

#define cRepCnt       11

#define cRepInvStatus 10
#define cRepInvColor  15
#define cRepInvText   16
#define cRepInvVisib  20
#define cRepInvPick   21
#define cRepInvCoord  30
#define cRepInvBonds  40
#define cRepInvAtoms  50
#define cRepInvAll    100

struct CoordSet;

typedef struct Rep {
  void            (*fRender)(struct Rep *I,CRay *ray,Pickable **pick);  
  struct Rep     *(*fUpdate)(struct Rep *I,struct CoordSet *cs,int rep);
  void        (*fInvalidate)(struct Rep *I,struct CoordSet *cs,int level);
  void              (*fFree)(struct Rep* I);
  int MaxInvalid,Active;
  Pickable *P;
  /* private */
  void        (*fRecolor)(struct Rep *I,struct CoordSet *cs);
  int         (*fSameVis)(struct Rep *I,struct CoordSet *cs);
  struct Rep *(*fRebuild)(struct Rep *I,struct CoordSet *cs,int rep);
  struct Rep *(*fNew)(struct CoordSet *cs);
} Rep;

void RepInit(Rep *I);
void RepFree(Rep *I);

#endif
