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
#ifndef _H_List
#define _H_List

/* All data structures which use List must have this header */

typedef struct {
  int link;
} ListElem;

typedef struct { 
  /* NOTE: happy lists must have record sizes equal or larger 
     than the length of this struct */
  int rec_size;
  int next_avail;
} List;

typedef struct {
  int link; 
  int value;
} ListInt;

typedef struct {
  int link; 
  int value[2];
} ListInt2;

typedef struct {
  int link; 
  int value[3];
} ListInt3;

typedef struct {
  int link; 
  int value[4];
} ListInt4;

typedef struct {
  int link; 
  int value[5];
} ListInt5;

void *ListNew(int init_size,int rec_size); 
int   ListElemNewZero(void *list_ptr);
int   ListElemNew(void *list_ptr_ptr);
void  ListFree(void *list);
void  ListElemFree(void *list,int elem);
void  ListElemFreeChain(void *list,int start);
int   ListLen(void *list,int start);
int   ListGetNAlloc(void *list);

/* shortcut routines for managing tops of lists */
int   ListElemPop(void *list,int elem);
int   ListElemPush(void *list_ptr_ptr,int elem);

int   ListElemPushInt(ListInt **list,int elem,int value);
int   ListElemPopInt(ListInt *list,int elem,int *value);
int   ListElemGetInt(ListInt *list,int elem,int *value);

int   ListElemPurgeInt(ListInt *list,int start, int value); /* slow! */
#endif









