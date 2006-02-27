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

#include "os_std.h"
#include "PyMOLGlobals.h"

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

  Notice that because of these requirements, it isn't easy
  to use this memory debuggin system with existing code...

*/

/* ==================== Master Switch ============================= 
 * Define _MemoryDebug_ON to enable the debugging system...*/

/* WARNING!!! MemoryDebug is not thread safe...it must be disabled
   for stable multi-threaded operation within the PyMOL core */

#define _MemoryDebug_OFF

/* ================================================================ 
 * Don't touch below unless you know what you are doing */


#define CopyArray(dst,src,type,count) memcpy(dst,src,sizeof(type)*(count))

void UtilMemCpy(void *dst,void *src,unsigned int *size);

typedef struct VLARec {
  unsigned int nAlloc;
  unsigned int recSize;
  float growFactor;
  int autoZero;
} VLARec;

/* NOTE: in VLACheck, rec is a zero based array index, not a record count */
#define VLACheck(ptr,type,rec) (ptr=(type*)(((((unsigned)rec)>=((VLARec*)(ptr))[-1].nAlloc) ? VLAExpand(ptr,(rec)) : (ptr))))

#define VLAlloc(type,initSize) (type*)VLAMalloc(initSize,sizeof(type),5,0)
#define VLACalloc(type,initSize) (type*)VLAMalloc(initSize,sizeof(type),5,1)
#define VLAFreeP(ptr) {if(ptr) {VLAFree(ptr);ptr=NULL;}}
#define VLASize(ptr,type,size) {ptr=(type*)VLASetSize(ptr,size);}
#define VLASizeForSure(ptr,type,size) {ptr=(type*)VLASetSizeForSure(ptr,size);}
#define VLACopy(ptr,type) (type*)VLANewCopy(ptr);

#define Alloc(type,size) (type*)mmalloc(sizeof(type)*(size))
#define Calloc(type,size) (type*)mcalloc(sizeof(type),size)
#define Realloc(ptr,type,size) (type*)mrealloc(ptr,sizeof(type)*(size))

#define FreeP(ptr) {if(ptr) {mfree(ptr);ptr=NULL;}}

void *VLAExpand(void *ptr,unsigned int rec); /* NOTE: rec is index (total-1) */
void *MemoryReallocForSure(void *ptr, unsigned int newSize);
void *MemoryReallocForSureSafe(void *ptr, unsigned int newSize, unsigned int oldSize);

#ifndef _MemoryDebug_ON
void *VLAMalloc(unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero); /*growfactor 1-10*/

#else
#define VLAMalloc(a,b,c,d) _VLAMalloc(__FILE__,__LINE__,a,b,c,d)

void *_VLAMalloc(const char *file,int line,unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero); /*growfactor 1-10*/
#endif


void VLAFree(void *ptr);
void *VLASetSize(void *ptr,unsigned int newSize);
void *VLASetSizeForSure(void *ptr,unsigned int newSize);

unsigned int VLAGetSize(void *ptr);
void *VLANewCopy(void *ptr);
void MemoryZero(char *p,char *q);

#ifndef _MemoryDebug_ON
/* _MemoryDebug_ON not defined */

#define mcalloc calloc
#define mmalloc malloc
#define mrealloc realloc
#define mfree free
#define mregister(x,y) 
#define mforget(x)
#define ReallocForSure(ptr,type,size) (type*)MemoryReallocForSure(ptr,sizeof(type)*(size))
#define ReallocForSureSafe(ptr,type,size,old_size) (type*)MemoryReallocForSure(ptr,sizeof(type)*(size),sizeof(type)*(old_size))

#define MemoryDebugDump()

#ifdef __cplusplus
#define mnew new
#endif

#define MD_FILE_LINE_Call
#define MD_FILE_LINE_Decl
#define MD_FILE_LINE_Nest
#define MD_FILE_LINE_PTR_Call


#else
/* _MemoryDebug_ON is defined */

#ifdef __cplusplus
extern "C" {
#endif

#define _MDPointer 1
#define _MDObject 2
#define _MDMarker 3

#define MD_FILE_LINE_Call ,__FILE__,__LINE__
#define MD_FILE_LINE_Decl ,const char *file,int line
#define MD_FILE_LINE_Nest ,file,line
#define MD_FILE_LINE_PTR_Call ,__FILE__,__LINE__,_MDPointer

#define mmalloc(x) MemoryDebugMalloc(x,__FILE__,__LINE__,_MDPointer)
#define mcalloc(x,y) MemoryDebugCalloc(x,y,__FILE__,__LINE__,_MDPointer)
#define mrealloc(x,y) MemoryDebugRealloc(x,y,__FILE__,__LINE__,_MDPointer)
#define mfree(x) MemoryDebugFree(x,__FILE__,__LINE__,_MDPointer)
#define mregister(x,y) MemoryDebugRegister((void*)x,y,__FILE__,__LINE__)
#define mforget(x) MemoryDebugForget((void*)x,__FILE__,__LINE__)
#define ReallocForSure(ptr,type,size) (type*)MemoryDebugReallocForSure(ptr,sizeof(type)*(size),__FILE__,__LINE__,_MDPointer)
#define ReallocForSureSafe(ptr,type,size,old_size) (type*)MemoryDebugReallocForSureSafe(ptr,sizeof(type)*(size),\
                    sizeof(type)*(old_size),__FILE__,__LINE__,_MDPointer)

void MemoryDebugRegister(void *addr,const char *note,
                         const char *file,int line);
void MemoryDebugForget(void *addr,const char *file,int line);

void *MemoryDebugMalloc(size_t size,const char *file,int line,int type);
void *MemoryDebugCalloc(size_t nmemb,size_t size,const char *file,int line,int type);
void *MemoryDebugRealloc(void *ptr,size_t size,
			 const char *file,int line,int type);
void *MemoryDebugReallocForSure(void *ptr,size_t size,const char *file,
                                int line,int type);

void MemoryDebugFree(void *ptr,const char *file,int line,int type);
void MemoryDebugQuietFree(void *ptr,int type);

void MemoryDebugDump(void);
int MemoryDebugUsage(void);

#ifdef __cplusplus
}

void *operator new(size_t size, const char *file,int line);

#define mnew new(__FILE__,__LINE__)

#endif

#endif

#endif







