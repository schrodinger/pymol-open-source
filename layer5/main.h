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

#include"os_python.h"

void MainFree(void);
void MainDoReshape(int width, int height);
void MainDirty(void);
void MainResetIdle(void);
void MainSwapBuffers(void);
void MainRefreshNow(void);
int MainSavingUnderWhileIdle(void);
int MainFromPyList(PyObject *list);
PyObject *MainAsPyList(void);

extern int PyMOLReady;
extern int PyMOLTerminating; /* flag to help prevent crashes on shutdown with Windows */
extern int PMGUI;
extern int StereoCapable;
extern int Security;

#ifdef _PYMOL_MODULE
int was_main(void);
#endif


#endif
