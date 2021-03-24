

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
#include"Result.h"

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

/**
 * @brief appends element to list
 * @param list list to append element to
 * @param ele element to append
 * @tparam ElemType type of the element
 */

template<typename ElemType>
void ListAppendT(ElemType*& list, ElemType* ele){
  ElemType* current = list;
  ElemType* previous = nullptr;
  while(current){
    previous = current;
    current = current->next;
  }
  if(previous){
    previous->next = ele;
  } else{
    list = ele;
  }
  ele->next = nullptr;
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

/**
 * @brief removes element from list
 * @param list list to remove element from
 * @param ele element to remove
 * @tparam ElemType type of the element
 * @return removed element
 * @note Removed element should be same as ele input.
 */

template<typename ElemType>
ElemType* ListDetachT(ElemType*& list, ElemType* ele){
  ElemType* current = list;
  ElemType* previous = nullptr;
  while(current){
    if(current == ele){
      break;
    }
    previous = current;
    current = current->next;
  }
  if(current){
    if(previous){
      previous->next = current->next;
    }
    else{
      list = current->next;
    }
    ele->next = nullptr;
  }
  return ele;
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

/**
 * Retrives position of element in list
 * @param ele target element
 * @return position in list
 */

template<typename ElemType>
pymol::Result<std::size_t> ListGetPosition(ElemType*& list, ElemType* ele)
{
  ElemType* current = list;
  std::size_t i = 0;
  for (; current; i++, current = current->next) {
    if (current == ele) {
      return i;
    }
  }
  return pymol::make_error("Element not found");
}

/**
 * Inserts an element at given position
 * @param ele target element
 * @param pos position in list
 */

template<typename ElemType>
pymol::Result<> ListInsertAt(ElemType*& list, ElemType* ele, std::size_t pos)
{
  ElemType* current = list;
  ElemType* previous = nullptr;
  std::size_t i = 0;
  while (current) {
    if (i++ == pos) {
      ele->next = current;
      if (previous) {
        previous->next = ele;
      }
      return {};
    }
    previous = current;
    current = current->next;
  }
  if (i == pos) {
    previous->next = ele;
    return {};
  }
  return pymol::make_error("Invalid pos");
}

#endif /* _H_ListMacros */
