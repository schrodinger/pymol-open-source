

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
  register ElemType *current = (List); \
  register ElemType *previous = NULL; \
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

#define ListInsert(List,Elem,After,Link,ElemType) \
{ \
  register ElemType *current = List; \
  register ElemType *previous = NULL; \
  while(current) \
		{ \
		  if(previous == (After)) \
			 break; \
   	  previous = current; \
		  current = current->Link; \
		} \
  if(previous) \
	 { \
		(Elem)->Link = current; \
		previous->Link = Elem; \
	 } \
  else \
	 { \
		(Elem)->Link = List; \
		(List) = Elem; \
	 } \
}

#define ListFree(List,Link,ElemType) \
{ \
  register ElemType *current = List; \
  register ElemType *previous = NULL; \
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
  register ElemType *current = List; \
  register ElemType *previous = NULL; \
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
   register ElemType *copy = (Elem); \
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
	(Elem) = (ElemType*)mmalloc(sizeof(ElemType));	\
	ErrChkPtr(G,Elem);				\
      }							\
}

#define ListElemCalloc(G,Elem,ElemType) \
{						\
  if(!(Elem))						  \
    {							  \
      (Elem) = (ElemType*)mcalloc(sizeof(ElemType),1);	  \
      ErrChkPtr(G,Elem);				  \
    }							  \
}

#define ListElemInit(List,Link) (List)->Link = NULL

#define ListElemFree(Elem) { mfree(Elem); Elem = NULL; }



/* -- JV
 * Doubly linked list macros
 * -- JV
 */

/* Create the circular, doubly linked list w/a sentinel node */
#define DListInit(List, Pre, Post, ElemType)	    \
do { \
  List = (ElemType*)malloc(sizeof(ElemType));	    \
  (List)->Pre = (List)->Post = List;		    \
} while(0)

/* DListInsert -- Insert Elem at head of list
 * List -- any structure with previous and next pointers (a doubly linked list)
 * Elem -- the Element to add
 * Pre  -- pointer to previous element in list
 * Post -- pointer to next element in list
 */
#define DListInsert(List,Elem,Pre,Post) \
do {\
  (Elem)->Post = List;			\
  (Elem)->Pre = (List)->Pre;		\
  (List)->Pre = (Elem);		\
  (Elem)->Pre->Post = (Elem);	\
} while(0)

/* DListRemove -- remove Element from the list, do not delete it
 * Elem -- the element to remove
 * Pre  -- the link to the previous element 
 * Post -- the link to the next element
 */
#define DListRemove(Elem,Pre,Post) \
do { \
  if ((Elem)->Pre && (Elem)->Post) {			  \
      (Elem)->Pre->Post = (Elem)->Post;			  \
      (Elem)->Post->Pre = (Elem)->Pre;			  \
  }							  \
  (Elem)->Pre = (Elem)->Post = NULL;			  \
} while (0)

/* DListIterate -- Iterate Elem across all items in List using Post as the link to next */
#define DListIterate(List,Elem,Post) \
    for((Elem) = (List)->Post; (Elem) != (List); (Elem) = (Elem)->Post) 

/* Can similarly do reverse iteration w/Post=Pre */

/* For all these ElemAlloc macros, it calls if(!Elem)
 * indicating that all blank incoming Elem's must be initialized
 * to NULL.  Just calling :ElemType* foo;" won't do.
 */
#define DListElemAlloc(G,Elem,ElemType) \
do {						\
 if(!(Elem))						\
   {							\
     (Elem) = (ElemType*)mmalloc(sizeof(ElemType));	\
     ErrChkPtr(G,Elem);					\
   }							\
 } while (0)

#define DListElemCalloc(G,Elem,ElemType) \
do {							  \
  if(!(Elem))						  \
    {							  \
      (Elem) = (ElemType*)mcalloc(sizeof(ElemType),1);	  \
      ErrChkPtr(G,Elem);				  \
    }							  \
} while (0)

#define DListElemInit(Elem,Pre,Post) (Elem)->Pre = (Elem)->Post = NULL

#define DListElemFree(Elem) { mfree(Elem); Elem = NULL; }

#endif /* _H_ListMacros */















