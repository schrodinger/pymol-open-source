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
#ifndef _H_P
#define _H_P

#include<Python.h>

#include"AtomInfo.h"

void PInit(void);
void PInitEmbedded(int argc,char **argv);
void PGetOptions(int *pmgui,int *internal_gui,int *stereo_capable);

void PFree(void);
void PExit(int code);
void PParse(char *str);

#define cLockAPI 1
#define cLockInbox 2
#define cLockOutbox 3

int PAlterAtom(AtomInfoType *at,char *expr);
int PLabelAtom(AtomInfoType *at,char *expr);
int PAlterAtomState(float *v,char *expr);

void PSleep(int usec);

void PLockAPIAsGlut(void);
void PUnlockAPIAsGlut(void);

void PBlock(void);
void PUnblock(void);

void PBlockAndUnlockAPI(void);
void PLockAPIAndUnblock(void);

void PFlush(void);
void PFlushFast(void);
void PXDecRef(PyObject *obj);

void PStereoOff(void);
void PDefineFloat(char *name,float value);

void PRunString(char *str);
void PDumpTraceback(PyObject *err);

int PTruthCallStr(PyObject *object,char *method,char *argument);

extern PyObject *P_globals;

extern PyObject *P_cmd;
extern PyObject *P_menu;
extern PyObject *P_xray;
extern PyObject *P_chempy;
extern PyObject *P_models;

extern PyThreadState *P_glut_thread_state; /* this is the state for the main GUI thread */
extern PyThreadState *P_api_thread_state; /* this is the thread state for a non-glut API thread */
extern int P_glut_thread_active; 
extern int P_glut_thread_keep_out;

#endif







