
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
#ifndef _H_RepSphere
#define _H_RepSphere

#include"Rep.h"
#include"CoordSet.h"

struct RepSphere : Rep {
  using Rep::Rep;

  ~RepSphere() override;

  cRep_t type() const override { return cRepSphere; }
  void render(RenderInfo* info) override;
  bool sameVis() const override;

  bool* LastVisib = nullptr;
  int* LastColor = nullptr;
  CGO* renderCGO = nullptr;
  CGO* primitiveCGO = nullptr;
  CGO* spheroidCGO = nullptr;
};

Rep *RepSphereNew(CoordSet * cset, int state);
void RenderSphereComputeFog(PyMOLGlobals *G, RenderInfo *info, float *fog_info);

#endif
