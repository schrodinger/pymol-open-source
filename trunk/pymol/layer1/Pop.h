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
#ifndef _H_Pop
#define _H_Pop

#include"Block.h"

int PopInit(PyMOLGlobals *G);
void PopFree(PyMOLGlobals *G);
void PopFitBlock(Block *block);
Block *PopGetBlock(PyMOLGlobals *G);
int PopPlaceChild(Block *block,int left_x,int right_x,int row_y,int affinity);

#endif
