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
#ifndef _H_MemoryDebug
#define _H_MemoryDebug

#include<stdlib.h>

/* This file can be included by C and C++ programs for
   debugging of malloc, realloc, free in C and in addition,
   of new and delete in C++ 

  it also includes a two functions for variable length arrays

  In order to use the debugging system, you must do the following...

  *use mmalloc for malloc
  *use mrealloc for mrealloc
  *use mfree for free
  *use mnew for new

  *use mregister to remember some address or pointer
  *use mforget to forget some address or pointer

  *note calloc isn't supported because I don't use it

  Notice that because of these requirements, it isn't easy
  to use this memory debuggin system with existing code...

*/

/* ==================== Master Switch ============================= 
 * Define _MemoryDebug_ON to enable the debugging system...*/

#define _MemoryDebug_ON

/* ================================================================ 
 * Don't touch below unless you know what you are doing */

typedef struct VLARec {
  unsigned int nAlloc;
  unsigned int recSize;
  unsigned int growFactor;
  int autoZero;
} VLARec;

#define VLACheck(ptr,type,rec) (ptr=(type*)((((rec)>=((VLARec*)(ptr))[-1].nAlloc) ? VLAExpand(ptr,(rec)) : (ptr))))

 /* NOTE: rec is index (total-1) */

#define VLAlloc(type,initSize) (type*)VLAMalloc(initSize,sizeof(type),2,0)
#define VLAFreeP(ptr) {if(ptr) {VLAFree(ptr);ptr=NULL;}}
#define Alloc(type,size) (type*)mmalloc(sizeof(type)*(size))
#define Realloc(ptr,type,size) (type*)mrealloc(ptr,sizeof(type)*(size))

#define FreeP(ptr) {if(ptr) {mfree(ptr);ptr=NULL;}}

void *VLAExpand(void *ptr,unsigned int rec); /* NOTE: rec is index (total-1) */

#ifndef _MemoryDebug_ON
void *VLAMalloc(unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero); /*growfactor 1-10*/
#else
#define VLAMalloc(a,b,c,d) _VLAMalloc(__FILE__,__LINE__,a,b,c,d)
void *_VLAMalloc(const char *file,int line,unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero); /*growfactor 1-10*/
#endif

void VLAFree(void *ptr);
void *VLASetSize(void *ptr,unsigned int newSize);

#ifndef _MemoryDebug_ON
/* _MemoryDebug_ON not defined */

#define mmalloc malloc
#define mrealloc realloc
#define mfree free
#define mregister(x,y) 
#define mforget(x)

#define MemoryDebugDump()

#ifdef __cplusplus
#define mnew new
#endif

#else
/* _MemoryDebug_ON is defined */

#include<stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define _MDPointer 1
#define _MDObject 2
#define _MDMarker 3

#define mmalloc(x) MemoryDebugMalloc(x,__FILE__,__LINE__,_MDPointer)
#define mrealloc(x,y) MemoryDebugRealloc(x,y,__FILE__,__LINE__,_MDPointer)
#define mfree(x) MemoryDebugFree(x,__FILE__,__LINE__,_MDPointer)
#define mregister(x,y) MemoryDebugRegister((void*)x,y,__FILE__,__LINE__)
#define mforget(x) MemoryDebugForget((void*)x,__FILE__,__LINE__)

void MemoryDebugRegister(void *addr,const char *note,
			 const char *file,int line);
void MemoryDebugForget(void *addr,const char *file,int line);

void *MemoryDebugMalloc(size_t size,const char *file,int line,int type);
void *MemoryDebugRealloc(void *ptr,size_t size,
			 const char *file,int line,int type);
void MemoryDebugFree(void *ptr,const char *file,int line,int type);
void MemoryDebugQuietFree(void *ptr,int type);

void MemoryDebugDump(void);

#ifdef __cplusplus
}

void *operator new(size_t size, const char *file,int line);

#define mnew new(__FILE__,__LINE__)

#endif

#endif
#endif







