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
#ifndef _H_Control
#define _H_Control

#include"Ortho.h"

typedef struct {
  Block *Block;
  int Rocking;
}  CControl;

extern CControl Control;

void ControlInit(void);
void ControlFree(void);
Block *ControlGetBlock(void);
int ControlIdling(void);
void ControlInterrupt(void);

#endif
