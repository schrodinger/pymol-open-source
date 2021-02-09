
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

#define cCartoon_skip_helix -2
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
#define cCartoon_cylinder 9

// show/hide/... codes
enum {
  cVis_HIDE,            // 0
  cVis_SHOW,            // 1
  cVis_AS,              // 2
  cVis_TOGGLE,          // 3
};

/* WARNING: don't change these -- you'll break sessions!
   (you can add to them however, I think) */

enum cRep_t {
  cRepNone = -2,
  cRepAll = -1,
  cRepCyl = 0,         // 0
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

inline cRep_t& operator++(cRep_t& rep)
{
  assert(0 <= rep && rep < cRepCnt);
  return (rep = cRep_t(rep + 1));
}

using cRepBitmask_t = int;

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
constexpr cRepBitmask_t cRepsAtomMask = (cRepCylBit | cRepSphereBit | cRepSurfaceBit |
    cRepLabelBit | cRepNonbondedSphereBit | cRepCartoonBit | cRepRibbonBit | \
    cRepLineBit | cRepMeshBit | cRepDotBit | cRepNonbondedBit | cRepEllipsoidBit);

// all reps which can be shown for objects
constexpr cRepBitmask_t cRepsObjectMask = (cRepSurfaceBit | cRepMeshBit | cRepDotBit |
    cRepCellBit | cRepCGOBit | cRepCallbackBit | cRepExtentBit | cRepSliceBit | \
    cRepAngleBit | cRepDihedralBit | cRepVolumeBit | cRepDashBit);

/* Hierarchical invalidation scheme - 
 * each higher level event implies all of the lower levels 
 * These used to be used just for graphics, but are now
 * used by the molecular editor as well */


/* invalite display (list) */

enum cRepInv_t {
  cRepInvNone = 0,
  cRepInvDisplay = 1,

/* precomputed extents (can change if matrix changes) */
  cRepInvExtents = 5,

/* invalidate pickable atoms */
  cRepInvPick = 9,

/* invalidate external atom colors */
  cRepInvExtColor = 10,

/* invalidate atom colors */
  cRepInvColor = 15,

/* invalidate label text */
  cRepInvText = 16,

/* invalidate visible atoms */
  cRepInvVisib = 20,
  cRepInvVisib2 = 21,

/* invalidate atomic properties */
  cRepInvProp = 22,

/* invalidate coordinates */
  cRepInvCoord = 30,

/* invalidate graphic representation */
  cRepInvRep = 35,

/* don't call ObjectMoleculeUpdateNonbonded */
  cRepInvBondsNoNonbonded = 38,

/* invalidate bond structure */
  cRepInvBonds = 40,

/* invalidate atomic structure */
  cRepInvAtoms = 50,

/* invalidate everything about a structure */
  cRepInvAll = 100,

/* invalidate and furthermore, purge existing representations */
  cRepInvPurgeMask = 0x80,
  cRepInvPurgeRep = cRepInvRep | cRepInvPurgeMask,
  cRepInvPurgeAll = cRepInvAll | cRepInvPurgeMask,

/* (alias) */
  cRepInvPurge = cRepInvPurgeRep,
};

struct CoordSet;
struct Object;

struct Rep {
  PyMOLGlobals *G;

  virtual cRep_t type() const = 0;
  virtual void render(RenderInfo* info);
  virtual void invalidate(cRepInv_t level);

  virtual ~Rep();

  pymol::CObject* obj = nullptr; // TODO redundant, use getObj()
  CoordSet* cs = nullptr;

  Pickable* P = nullptr; //!< only used by labels
  PickContext context{};

  //! Object state (>=0, 0-indexed) for picking and ramp colors
  int getState() const { return context.state; }
  pymol::CObject* getObj() const { return context.object; }

  //! True if this rep should be rendered in RenderPass::Transparent. Default is
  //! false, change the flag with setHasTransparency().
  bool hasTransparency() const { return m_has_transparency; }
  //! Sets the hasTransparency() flag.
  void setHasTransparency(bool has = true) { m_has_transparency = has; }

protected:
  cRepInv_t MaxInvalid = cRepInvNone;

private:
  Rep* rebuild();
  virtual Rep* recolor() { return rebuild(); }
  virtual bool sameVis() const { return false; }
  virtual bool sameColor() const { return false; }

  bool m_has_transparency = false;

public:
  Rep* update();

  /** Pointer to static factory function (Only used with molecular
   * representations, DistSet e.g. doesn't use it)
   * @param state Object state for picking and ramp colors
   */
  Rep* (*fNew)(CoordSet* cs, int state) = nullptr;

  Rep(pymol::CObject*, int state);
  Rep(CoordSet*, int state);
};

cRepBitmask_t RepGetAutoShowMask(PyMOLGlobals * G);

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
