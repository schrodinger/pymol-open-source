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

typedef struct CSculpt {
  CShaker *Shaker;
  ObjectMolecule *Obj;
  int *NBHash;
  int *NBList;
  int *EXHash;
  int *EXList;
} CSculpt;

CSculpt *SculptNew(void);
void SculptMeasureObject(CSculpt *I,ObjectMolecule *obj,int state);
void SculptIterateObject(CSculpt *I,ObjectMolecule *obj,int state,int n_cycle);

void SculptFree(CSculpt *I);

#endif


