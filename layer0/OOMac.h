

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
#ifndef _H_OOMac
#define _H_OOMac

#include"Err.h"
#include"MemoryDebug.h"

#define OOAlloc(G,type) \
auto* I = new type; \
ErrChkPtr(G,I);

#define OOCalloc(G,type) \
auto* I = new type(); \
ErrChkPtr(G,I);

#define OOFreeP(ptr) \
{if(ptr) {delete(ptr);ptr=NULL;}}

#endif
