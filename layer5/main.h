
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

#ifdef _PYMOL_NO_MAIN


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _MACPYMOL_XCODE
#include"os_python.h"
#include"PyMOLGlobals.h"
PyObject *MainAsPyList(void);
#endif

/* END PROPRIETARY CODE SEGMENT */

#else

#include"os_python.h"
#include"PyMOLGlobals.h"

void MainFree(void);
void MainDoReshape(int width, int height);
void MainRefreshNow(void);
void MainFlush(void);
void MainFlushAsync(void);
int MainSavingUnderWhileIdle(void);

int MainFromPyList(PyObject * list);
PyObject *MainAsPyList(void);

void MainSetWindowVisibility(int mode);
void MainMaximizeWindow(PyMOLGlobals * G);
void MainSetWindowSize(PyMOLGlobals * G, int w, int h);
void MainSetWindowPosition(PyMOLGlobals * G, int x, int y);
void MainCheckWindowFit(PyMOLGlobals * G);

#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _MACPYMOL_XCODE
void MainUnblock(void);
int MainCheckRedundantOpen(char *file);
void MainRunString(const char *str);
PyObject *MainGetStringResult(const char *str);
void MainDoCommand(const char *str1);
int MainFeedbackOut(char *st);
void MainRunCommand(const char *str1);
void MainMoviePrepareCopy(int *width, int *height, int *length);
int MainMovieCopy(int frame, int width, int height, int rowbytes, void *ptr);
void MainMovieCopyPrepare(int *width, int *height, int *length);
int MainMovieCopyFrame(int frame, int width, int height, int rowbytes, void *ptr);
int MainMoviePurgeFrame(int frame);
void MainMovieCopyFinish(void);
void MainSceneGetSize(int *width, int *height);
int MainSceneCopy(int width, int height, int rowbytes, void *ptr);
#endif

/* END PROPRIETARY CODE SEGMENT */

#define PYMOL_MAX_OPT_STR  1025

int main_exec(int, char **);
int main_shared(int);

#endif /* _H_main header */
