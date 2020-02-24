

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
#ifndef _H_ListMacros
#define _H_ListMacros


/* simplest possible single-linked list */

#include"MemoryDebug.h"
#include"Err.h"

#define ListInit(List) List = NULL

#define ListAppend(List,Elem,Link,ElemType) \
{ \
  ElemType *current = (List); \
  ElemType *previous = NULL; \
  while(current) \
	 { \
		previous = current; \
		current = current->Link; \
	 } \
  if(previous) \
	 previous->Link = Elem; \
  else \
	 (List) = Elem; \
  (Elem)->Link = NULL; \
}

#define ListFree(List,Link,ElemType) \
{ \
  ElemType *current = List; \
  ElemType *previous = NULL; \
  while(current) \
	 { \
      if(previous) \
		  mfree(previous); \
		previous = current; \
		current = current->Link; \
	 } \
  if(previous) \
	 mfree(previous); \
  (List) = NULL; \
}

#define ListDetach(List,Elem,Link,ElemType) \
{ \
  ElemType *current = List; \
  ElemType *previous = NULL; \
  while(current) \
	 { \
		if(current == (Elem)) \
		  break; \
		previous = current; \
		current = current->Link; \
	 } \
  if(current) \
	 { \
		if(previous) \
		  previous->Link = current->Link; \
        else \
          (List) = current->Link; \
	  (Elem)->Link = NULL; \
	 } \
}

#define ListDelete(List,Elem,Link,ElemType) \
{ \
   ElemType *copy = (Elem); \
   ListDetach(List,copy,Link,ElemType); \
   mfree(copy); \
}

#define ListIterate(List,Counter,Link) \
   ( (Counter) = ((List) ? (((Counter) ? (Counter)->Link : (List))) : NULL))


/* Elem handling routines */

#define ListElemAlloc(G,Elem,ElemType)		\
{						\
    if(!(Elem))					\
      {							\
	(Elem) = pymol::malloc<ElemType>(1);		\
	ErrChkPtr(G,Elem);				\
      }							\
}

#define ListElemCalloc(G,Elem,ElemType) \
{						\
  if(!(Elem))						  \
    {							  \
      (Elem) = pymol::calloc<ElemType>(1);		  \
      ErrChkPtr(G,Elem);				  \
    }							  \
}

#define ListElemFree(Elem) { mfree(Elem); Elem = NULL; }

#endif /* _H_ListMacros */
