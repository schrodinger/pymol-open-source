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
#ifndef _H_PUtils
#define _H_PUtils

#include<Python.h>

#include"AtomInfo.h"

void PInit(void);
void PFree(void);
void PExit(int code);
void PParse(char *str);

#define cLockAPI 1
#define cLockInbox 2
#define cLockOutbox 3

int PAlterAtom(AtomInfoType *at,char *expr);
void PSleep(int usec);

void PLockAPIAsGlut(void);
void PUnlockAPIAsGlut(void);

void PBlock(void);
void PUnblock(void);

void PBlockAndUnlockAPI(void);
void PLockAPIAndUnblock(void);

void PFlush(void);

void PStereoOff(void);
void PDefineFloat(char *name,float value);

PyObject *PFloatVLAToPyList(float *f);

extern PyObject *P_pm;
extern PyObject *P_pmm;
extern PyObject *P_pmx;

extern PyThreadState *P_glut_thread_state; /* this is the state for the main GUI thread */
extern PyThreadState *P_api_thread_state; /* this is the thread state for a non-glut API thread */
extern int P_glut_thread_active; 
extern int P_glut_thread_keep_out;

#endif







