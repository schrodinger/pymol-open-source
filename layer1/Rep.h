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

#define cCartoon_skip -1
#define cCartoon_auto 0
#define cCartoon_loop 1
#define cCartoon_rect 2
#define cCartoon_oval 3
#define cCartoon_tube 4
#define cCartoon_arrow 5
#define cCartoon_dumbbell 6
#define cCartoon_putty 7

#define cCartoon_skip_helix -2
#define cRepAll       -1
#define cRepNone      -2

/* WARNING: don't change these -- you'll break sessions!
   (you can add to them however, I think) */

#define cRepCyl        0
#define cRepSphere     1
#define cRepSurface    2
#define cRepLabel      3
#define cRepNonbondedSphere 4
#define cRepCartoon    5
#define cRepRibbon     6
#define cRepLine       7
#define cRepMesh       8
#define cRepDot        9
#define cRepDash       10
#define cRepNonbonded  11
#define cRepCell       12
#define cRepCGO        13
#define cRepCallback   14
#define cRepExtent     15
#define cRepSlice      16
#define cRepAngle      17
#define cRepDihedral   18
#define cRepEllipsoid  19

#define cRepCnt        20

/* Hierarchical invalidation scheme - 
 * each higher level event implies all of the lower levels 
 * These used to be used just for graphics, but are now
 * used by the molecular editor as well */

/* invalite display (list) */

#define cRepInvDisplay 1
/* precomputed extents (can change if matrix changes) */
#define cRepInvExtents 5
/* invalidate pickable atoms */
#define cRepInvPick  9
/* invalidate external atom colors */
#define cRepInvExtColor  10
/* invalidate atom colors */
#define cRepInvColor  15
/* invalidate label text */
#define cRepInvText   16
/* invalidate visible atoms */
#define cRepInvVisib  20
#define cRepInvVisib2 21
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
/* invalidate and furthermore, purge existing representations */
#define cRepInvPurge  110

struct CoordSet;
struct Object;



typedef struct Rep {
  PyMOLGlobals *G;
  void            (*fRender)(struct Rep *I,RenderInfo *info);
  struct Rep     *(*fUpdate)(struct Rep *I,struct CoordSet *cs,int state,int rep);
  void        (*fInvalidate)(struct Rep *I,struct CoordSet *cs,int level);
  void              (*fFree)(struct Rep* I);
  int MaxInvalid,Active;
  CObject *obj;
  struct CoordSet *cs;
  Pickable *P;
  PickContext context;
  /* private */
  void        (*fRecolor)(struct Rep *I,struct CoordSet *cs);
  int         (*fSameVis)(struct Rep *I,struct CoordSet *cs);
  struct Rep *(*fRebuild)(struct Rep *I,struct CoordSet *cs,int state,int rep);
  struct Rep *(*fNew)(struct CoordSet *cs,int state);
  int displayList;
  int displayListInvalid;
} Rep;

void RepInit(PyMOLGlobals *G,Rep *I);
void RepPurge(Rep *I);

#endif
