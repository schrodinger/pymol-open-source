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
#ifndef _H_mac
#define _H_mac

#include "os_memory.h"

#ifndef NULL
#define NULL ((void*)0)
#endif

#define mac_calloc(type) ((type*)os_calloc(1,sizeof(type)))
/* allocate a new, blank object of given type */

#define mac_calloc_array(type,cnt) ((type*)os_calloc(cnt,sizeof(type)))
/* allocate a new, blank array of given type */

#define mac_malloc(type) ((type*)os_malloc(sizeof(type)))
/* allocate a block of given type (not blanked) */

#define mac_malloc_array(type,cnt) ((type*)os_malloc(sizeof(type)*(cnt)))
/* allocate an array of blocks of given type (not blanked) */

#define mac_free(ptr) {if(ptr) {os_free(ptr);(ptr)=NULL;}}

/* free a pointer if it is non-null, then nullify */

#endif

