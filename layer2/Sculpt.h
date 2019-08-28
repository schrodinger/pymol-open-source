
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
#ifndef _H_Sculpt
#define _H_Sculpt

#include"Shaker.h"
#include"vla.h"

#include <memory>

#define cSculptBond  0x001
#define cSculptAngl  0x002
#define cSculptPyra  0x004
#define cSculptPlan  0x008
#define cSculptLine  0x010
#define cSculptVDW   0x020
#define cSculptVDW14 0x040
#define cSculptTors  0x080
#define cSculptTri   0x100
#define cSculptMin   0x200
#define cSculptMax   0x400
#define cSculptAvoid 0x800

struct CSculpt {
  PyMOLGlobals *G;
  std::unique_ptr<CShaker> Shaker;
  ObjectMolecule *Obj;
  std::vector<int> NBHash;
  pymol::vla<int> NBList;
  std::vector<int> EXHash;
  pymol::vla<int> EXList;
  pymol::vla<int> Don;
  pymol::vla<int> Acc;
  float inverse[256];
  CSculpt(PyMOLGlobals * G);
};

struct ObjectMolecule;

void SculptMeasureObject(CSculpt * I, ObjectMolecule * obj, int state, int match_state,
                         int match_by_segment);
float SculptIterateObject(CSculpt * I, ObjectMolecule * obj, int state, int n_cycle,
                          float *center);


#endif
