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
void PLock(int lock,PyThreadState **save);
void PUnlock(int lock,PyThreadState **save);
void PFlush(PyThreadState **save);
void PStereoOff(void);
void PDefineFloat(char *name,float value);

#endif
