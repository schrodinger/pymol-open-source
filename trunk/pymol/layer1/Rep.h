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


/* Hierarchical invalidation scheme - 
 * each higher level event implies all of the lower levels 
 * These used to be used just for graphics, but are now
 * used by the molecular editor as well */

/* invalidate atom colors */
#define cRepInvColor  15
/* invalidate label text */
#define cRepInvText   16
/* invalidate visible atoms */
#define cRepInvVisib  20
/* invalidate picked atoms */
#define cRepInvPick   21
/* invalidate atomic properties */
#define cRepInvProp   22
/* invalidate coordinates */
#define cRepInvCoord  30
/* invalidate graphic representation */
#define cRepInvRep    35
/* invalidate bond structure */
#define cRepInvBonds  40
/* invalidate atomic structure */
#define cRepInvAtoms  50
/* invalidate everything about a structure */
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
