

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

struct CDeferred {
  PyMOLGlobals *m_G { nullptr };
  int (*fn)(CDeferred *) = nullptr;
  CDeferred *next { nullptr };
  void exec();
  CDeferred(PyMOLGlobals * G) : m_G(G){};
  virtual ~CDeferred() = default;
};

typedef int DeferredFn(CDeferred * D);

#endif
