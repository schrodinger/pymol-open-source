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

#define cSculptBond  0x01
#define cSculptAngl  0x02
#define cSculptPyra  0x04
#define cSculptPlan  0x08
#define cSculptLine  0x10
#define cSculptVDW   0x20
#define cSculptVDW14 0x40
#define cSculptTors  0x80

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
void SculptMeasureObject(CSculpt *I,ObjectMolecule *obj,int state);
float SculptIterateObject(CSculpt *I,ObjectMolecule *obj,int state,int n_cycle);

void SculptFree(CSculpt *I);

#endif


