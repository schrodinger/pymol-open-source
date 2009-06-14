
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
#ifndef _H_RepNonbonded
#define _H_RepNonbonded

#include"Rep.h"
#include"CoordSet.h"

Rep *RepNonbondedNew(CoordSet * cset, int state);
void RepNonbondedRenderImmediate(CoordSet * cs, RenderInfo * info);

#endif
