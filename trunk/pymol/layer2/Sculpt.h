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
#include"ObjectMolecule.h"

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

typedef struct CSculpt {
  PyMOLGlobals *G;
  CShaker *Shaker;
  ObjectMolecule *Obj;
  int *NBHash;
  int *NBList;
  int *EXHash;
  int *EXList;
  int *Don;
  int *Acc;
  float inverse[256];
} CSculpt;

CSculpt *SculptNew(PyMOLGlobals *G);
void SculptMeasureObject(CSculpt *I,ObjectMolecule *obj,int state,int match_state,int match_by_segment);
float SculptIterateObject(CSculpt *I,ObjectMolecule *obj,
                          int state,int n_cycle,float *center);

void SculptFree(CSculpt *I);

#endif


