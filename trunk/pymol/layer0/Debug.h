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
#ifndef _H_DEBUG
#define _H_DEBUG

/* OBSOLETE -- USE THE FEEDBACK FACILITY */

#define DebugSelector (1)
#define DebugParser   (1<<1)
#define DebugMolecule (1<<2)
#define DebugPython   (1<<3)
#define DebugMap      (1<<4)

extern unsigned int DebugState;

#endif
