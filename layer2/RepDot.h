
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
#ifndef _H_RepDot
#define _H_RepDot

#include"Rep.h"
#include"CoordSet.h"

enum cRepDot_t {
  cRepDotNormal = 0,
  cRepDotAreaType = 1,
};

struct RepDot : Rep {
  using Rep::Rep;

  ~RepDot() override;

  cRep_t type() const override { return cRepDot; }
  void render(RenderInfo* info) override;

  float dotSize;
  float *V = nullptr;
  float *VC = nullptr;
  float *A = nullptr;                     //!< area
  float *VN = nullptr;                    //!< vector normal
  int *T = nullptr;                       //!< custom type
  int *F = nullptr;                       //!< flags
  int N = 0;
  int *Atom = nullptr;                    //!< atom
  float Width;
  CGO* shaderCGO = nullptr;
  bool shaderCGO_as_spheres = false;
};

Rep *RepDotNew(CoordSet * cset, int state);
Rep *RepDotDoNew(CoordSet * cs, cRepDot_t mode, int state);

#endif
