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
#ifndef _H_Crystal
#define _H_Crystal

#include"Vector.h"

typedef struct { 
  float Dim[3];
  float Angle[3]; /* stored in degrees for convenience */
  float RealToFrac[9];
  float FracToReal[9];
  float UnitCellVolume;
  float Norm[3];
  float RecipDim[3];
} CCrystal;

void CrystalFree(CCrystal *I);
void CrystalInit(CCrystal *I);
CCrystal *CrystalNew(void);
void CrystalUpdate(CCrystal *I);
void CrystalDump(CCrystal *I);

#endif









