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
#ifndef _H_os_memory
#define _H_os_memory

/* This file can be included by C programs for debugging of malloc, realloc, free in C 

  In order to use the debugging system, you must do the following...

  *use os_malloc for malloc
  *use os_realloc for realloc
  *use os_calloc for calloc
  *use os_free for free

*/
#include<stdlib.h>

/* ==================== Master Switch ============================= 
 * Define _os_memory_debug_on to enable the debugging system...
 * Note: the debugging system is not thread-safe, so turn this
 *       flag off when you are using threads
*/

#define _os_memory_debug_on

/* ================================================================ 
 * Don't touch below unless you know what you are doing */


#ifndef _os_memory_debug_on

#define os_malloc malloc
#define os_realloc realloc
#define os_calloc calloc
#define os_free free

#define os_memory_dump()

#else

#define _OSMemoryPointer 1
#define _OSMemoryVLA     2

void *OSMemoryMalloc(unsigned int size,const char *file,int line,int type);
void *OSMemoryCalloc(unsigned int count,unsigned int size,const char *file,int line,int type);
void *OSMemoryRealloc(void *ptr,unsigned int size,const char *file,int line,int type);
void OSMemoryFree(void *ptr,const char *file,int line,int type);
void OSMemoryDump(void);

#define os_malloc(x) OSMemoryMalloc(x,__FILE__,__LINE__,_OSMemoryPointer)
#define os_calloc(n,x) OSMemoryCalloc(n,x,__FILE__,__LINE__,_OSMemoryPointer)
#define os_realloc(x,y) OSMemoryRealloc(x,y,__FILE__,__LINE__,_OSMemoryPointer)
#define os_free(x) OSMemoryFree(x,__FILE__,__LINE__,_OSMemoryPointer)
#define os_memory_dump() OSMemoryDump()


#endif

void OSMemoryZero(char *p,char *q);

#define os_zero(p,q) OSMemoryZero(p,q)

#endif







