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
#ifndef _H_OOMac
#define _H_OOMac

#include"Err.h"
#include"MemoryDebug.h"

#define OOAlloc(G,type) \
type *I;\
I = (type*)mmalloc(sizeof(type));\
ErrChkPtr(G,I);

#define OOCalloc(G,type) \
type *I;\
  I = (type*)mcalloc(sizeof(type),1);            \
ErrChkPtr(G,I);

#define OOFreeP(ptr) \
{if(ptr) {mfree(ptr);ptr=NULL;}}

#endif
