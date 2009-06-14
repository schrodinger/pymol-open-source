
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
#ifndef _H_RepRibbon
#define _H_RepRibbon

#include"Rep.h"
#include"CoordSet.h"

Rep *RepRibbonNew(CoordSet * cset, int state);
void RepRibbonInit(void);
void RepRibbonRenderImmediate(CoordSet * cs, RenderInfo * info);

#endif
