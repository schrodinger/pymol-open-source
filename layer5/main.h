
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

int MainSavingUnderWhileIdle(void);

#include"os_python.h"
#include"PyMOLGlobals.h"

PyObject *MainAsPyList(PyMOLGlobals *G);
int MainFromPyList(PyMOLGlobals *G, PyObject * list);

#ifndef _PYMOL_NO_MAIN

void MainFree(void);
void MainDoReshape(int width, int height);
void MainRefreshNow(void);
void MainFlush(void);
void MainFlushAsync(void);
int MainSavingUnderWhileIdle(void);

void MainSetWindowVisibility(int mode);
void MainMaximizeWindow(PyMOLGlobals * G);
void MainSetWindowSize(PyMOLGlobals * G, int w, int h);
void MainSetWindowPosition(PyMOLGlobals * G, int x, int y);
void MainCheckWindowFit(PyMOLGlobals * G);

#endif

#define PYMOL_MAX_OPT_STR  1025

int main_shared(int);

#endif /* _H_main header */
