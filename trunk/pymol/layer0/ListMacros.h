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
#ifndef _H_ListMacros
#define _H_ListMacros

#include"MemoryDebug.h"
#include"Err.h"

#define ListVarDeclare(ListType,ElemType) \
   static struct { ElemType *List,*current,*previous,*copy,*storage; int count,flag; } ListType

#define ListInit(List) List = NULL

#define ListAppend(List,Elem,Link,ListType) \
{ \
  ListType.current = (List); \
  ListType.previous = NULL; \
  while(ListType.current) \
	 { \
		ListType.previous = ListType.current; \
		ListType.current = ListType.current->Link; \
	 } \
  if(ListType.previous) \
	 ListType.previous->Link = Elem; \
  else \
	 (List) = Elem; \
  (Elem)->Link = NULL; \
}

#define ListIndex(List,Elem,c,Link,ListType) \
{ \
  ListType.current = (List); \
  ListType.count = (c);\
  while(ListType.count--) \
	 { \
      if(!ListType.current) {break;} \
		ListType.current = ListType.current->Link; \
	 } \
  Elem=ListType.current;\
}

#define ListSort(List,fOrdered,Link,ListType) \
{ \
  ListType.flag = 1; \
  while(ListType.flag) \
	 { \
		ListType.flag = 0; \
		ListType.current = List; \
		ListType.previous = NULL; \
		while(ListType.current) \
		  { \
			 if(ListType.current->Link) \
				{ \
				  if(!fOrdered(ListType.current,ListType.current->Link)) \
					 { \
						ListType.flag = 1; \
						if(ListType.previous) \
						  { \
							 ListType.previous->Link = ListType.current->Link; \
							 ListType.current->Link = ListType.previous->Link->Link; \
							 ListType.previous->Link->Link = ListType.current; \
						  } \
						else \
						  { \
							 (List) = ListType.current->Link; \
							 ListType.current->Link = (List)->Link; \
							 (List)->Link->Link = ListType.current; \
						  } \
					 } \
				} \
			 ListType.previous = ListType.current; \
			 ListType.current = ListType.current->Link; \
		  } \
	 } \
}

#define ListInsert(List,Elem,After,Link,ListType) \
{ \
  ListType.copy = After; \
  ListType.current = List; \
  ListType.previous = NULL; \
  while(ListType.current) \
		{ \
		  if(ListType.previous == (After)) \
			 break; \
   	  ListType.previous = ListType.current; \
		  ListType.current = ListType.current->Link; \
		} \
  if(ListType.previous) \
	 { \
		(Elem)->Link = ListType.current; \
		ListType.previous->Link = Elem; \
	 } \
  else \
	 { \
		(Elem)->Link = List; \
		(List) = Elem; \
	 } \
}

#define ListInsertSorted(List,Elem,fOrdered,Link,ListType) \
{ \
  ListType.current = (List); \
  ListType.previous = NULL; \
  while(ListType.current) \
		{ \
		  if(fOrdered(Elem,ListType.current)) \
			 break; \
   	  ListType.previous = ListType.current; \
		  ListType.current = ListType.current->Link; \
		} \
  if(ListType.previous) \
	 { \
		(Elem)->Link = ListType.current; \
		ListType.previous->Link = Elem; \
	 } \
  else \
	 { \
		(Elem)->Link = List; \
		(List) = Elem; \
	 } \
}


#define ListFree(List,Link,ListType) \
{ \
  ListType.current = List; \
  ListType.previous = NULL; \
  while(ListType.current) \
	 { \
      if(ListType.previous) \
		  mfree(ListType.previous); \
		ListType.previous = ListType.current; \
		ListType.current = ListType.current->Link; \
	 } \
  if(ListType.previous) \
	 mfree(ListType.previous); \
  (List) = NULL; \
}

#define ListDetach(List,Elem,Link,ListType) \
{ \
  ListType.copy = Elem; \
  ListType.current = List; \
  ListType.previous = NULL; \
  while(ListType.current) \
	 { \
		if(ListType.current == (Elem)) \
		  break; \
		ListType.previous = ListType.current; \
		ListType.current = ListType.current->Link; \
	 } \
  if(ListType.current) \
	 { \
		if(ListType.previous) \
        ListType.previous->Link = ListType.current->Link; \
      else \
        (List) = ListType.current->Link; \
		(Elem)->Link = ListType.current; \
	 } \
}


#define ListDelete(List,Elem,Link,ListType) \
{ \
   ListType.copy = (Elem); \
   ListDetach(List,ListType.copy,Link,ListType); \
   mfree(ListType.copy); \
}

#define ListIterate(List,Counter,Link,ListType) \
   ( (Counter) = ((Counter) ? (Counter)->Link : (List)))

/* Elem handling routines */

#define ListElemAlloc(Elem,ElemType) \
{ \
if(!(Elem)) \
  { \
	 (Elem) = (ElemType*)mmalloc(sizeof(ElemType)); \
	 ErrChkPtr(Elem); \
  } \
}


#define ListElemInit(List,Link) (List)->Link = NULL

#define ListElemFree(Elem) { mfree(Elem); Elem = NULL; }


#endif
