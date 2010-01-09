

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
#ifndef _H_Deferred
#define _H_Deferred

#include "PyMOLGlobals.h"

typedef struct _CDeferred CDeferred;

typedef int DeferredFn(CDeferred * D);

void DeferredInit(PyMOLGlobals * G, CDeferred * I);
void DeferredFree(CDeferred * I);
CDeferred *DeferredExec(CDeferred * I);

struct _CDeferred {
  PyMOLGlobals *G;
  DeferredFn *fn;
  CDeferred *next;
};

#endif
