/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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
void MainFlush(void);
int MainSavingUnderWhileIdle(void);
int MainFromPyList(PyObject *list);
PyObject *MainAsPyList(void);
int MainCheckRedundantOpen(char *file);
  
extern int PyMOLReady;
extern int PyMOLTerminating; /* flag to help prevent crashes on shutdown with Windows */

/* note that the following four options are also contained in
   PyMOLOption global -- they exist here as independent globals only
   because changes haven't yet been made throughout code */

extern int PMGUI;
extern int StereoCapable;
extern int Security;
extern int ExternalGUI;

#define PYMOL_MAX_OPT_STR  1025

typedef struct PyMOLOptionRec {
  /* WARNING: don't delete items or change order unless you also update
     main.c where this global structure is initialized */
  int pmgui, internal_gui, show_splash,
    internal_feedback, security, game_mode,
    force_stereo, winX, winY, blue_line,
    winPX, winPY, external_gui, siginthand,
    reuse_helper, auto_reinitialize;
  char after_load_script[PYMOL_MAX_OPT_STR];
} PyMOLOptionRec;

extern PyMOLOptionRec *PyMOLOption;

#ifdef _PYMOL_MODULE
int was_main(void);
#endif


#endif
