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
#ifndef _H_main
#define _H_main

void MainFree(void);
void MainDoReshape(int width, int height);
void MainDirty(void);
void MainSwapBuffers(void);
void MainRefreshNow(void);

extern int PyMOLReady;
extern int PMGUI;

#ifdef _PYMOL_MODULE
void was_main(void);
#endif

#endif
