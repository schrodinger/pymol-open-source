
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
#define cCartoon_dash 8

#define cCartoon_skip_helix -2
#define cRepAll       -1
#define cRepNone      -2


/* WARNING: don't change these -- you'll break sessions!
   (you can add to them however, I think) */

enum {
  cRepCyl,             // 0
  cRepSphere,          // 1
  cRepSurface,         // 2
  cRepLabel,           // 3
  cRepNonbondedSphere, // 4
  cRepCartoon,         // 5
  cRepRibbon,          // 6
  cRepLine,            // 7
  cRepMesh,            // 8
  cRepDot,             // 9
  cRepDash,            // 10
  cRepNonbonded,       // 11
  cRepCell,            // 12
  cRepCGO,             // 13
  cRepCallback,        // 14
  cRepExtent,          // 15
  cRepSlice,           // 16
  cRepAngle,           // 17
  cRepDihedral,        // 18
  cRepEllipsoid,       // 19
  cRepVolume,          // 20
  // rep count
  cRepCnt
};

#define cRepCylBit             (1 << 0)
#define cRepSphereBit          (1 << 1)
#define cRepSurfaceBit         (1 << 2)
#define cRepLabelBit           (1 << 3)
#define cRepNonbondedSphereBit (1 << 4)
#define cRepCartoonBit         (1 << 5)
#define cRepRibbonBit          (1 << 6)
#define cRepLineBit            (1 << 7)
#define cRepMeshBit            (1 << 8)
#define cRepDotBit             (1 << 9)
#define cRepDashBit            (1 << 10)
#define cRepNonbondedBit       (1 << 11)
#define cRepCellBit            (1 << 12)
#define cRepCGOBit             (1 << 13)
#define cRepCallbackBit        (1 << 14)
#define cRepExtentBit          (1 << 15)
#define cRepSliceBit           (1 << 16)
#define cRepAngleBit           (1 << 17)
#define cRepDihedralBit        (1 << 18)
#define cRepEllipsoidBit       (1 << 19)
#define cRepVolumeBit          (1 << 20)

/* Add other reps here.  Don't forget to
 * update modules/constants.py::repres{}
 * update modules/constants.py::fb_module, if needed
 * update modules/viewing.py::rep_list
 * create your RepXYZ.h and RepXYZ.c
 */

#define cRepBitmask        ((1 << cRepCnt) - 1)

// all reps which can be shown for atoms
const int cRepsAtomMask = (cRepCylBit | cRepSphereBit | cRepSurfaceBit | \
    cRepLabelBit | cRepNonbondedSphereBit | cRepCartoonBit | cRepRibbonBit | \
    cRepLineBit | cRepMeshBit | cRepDotBit | cRepNonbondedBit | cRepEllipsoidBit);

// all reps which can be shown for objects
const int cRepsObjectMask = (cRepSurfaceBit | cRepMeshBit | cRepDotBit | \
    cRepCellBit | cRepCGOBit | cRepCallbackBit | cRepExtentBit | cRepSliceBit | \
    cRepAngleBit | cRepDihedralBit | cRepVolumeBit | cRepDashBit);

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
  void (*fRender) (struct Rep * I, RenderInfo * info);
  struct Rep *(*fUpdate) (struct Rep * I, struct CoordSet * cs, int state, int rep);
  void (*fInvalidate) (struct Rep * I, struct CoordSet * cs, int level);
  void (*fFree) (struct Rep * I);
  int MaxInvalid, Active;
  CObject *obj;
  struct CoordSet *cs;
  Pickable *P;
  PickContext context;
  /* private */
  void (*fRecolor) (struct Rep * I, struct CoordSet * cs);
  int (*fSameVis) (struct Rep * I, struct CoordSet * cs);
  int (*fSameColor) (struct Rep * I, struct CoordSet * cs);
  struct Rep *(*fRebuild) (struct Rep * I, struct CoordSet * cs, int state, int rep);
  struct Rep *(*fNew) (struct CoordSet * cs, int state);
} Rep;

void RepInit(PyMOLGlobals * G, Rep * I);
void RepPurge(Rep * I);
void RepInvalidate(struct Rep *I, struct CoordSet *cs, int level);

int RepGetAutoShowMask(PyMOLGlobals * G);

class RepIterator {
  int end;

public:
  int rep;

  RepIterator(PyMOLGlobals * G, int rep_);

  bool next() {
    return (++rep < end);
  };
};

#endif
