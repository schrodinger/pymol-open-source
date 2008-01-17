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
#ifndef _H_vla
#define _H_vla

#include "os_memory.h"

/* WARNING, almost all vla macros have the side effect of changing the pointer...
   so be careful not to ever use temporary copies of vla pointers with these operations */

typedef struct VLARec {
  unsigned int nAlloc;
  unsigned int recSize;
  unsigned int growFactor;
  int autoZero;
} VLARec;

/* NOTE: in vla_check, rec is a zero based array index, not a record count */

#define vla_check(ptr,type,rec) (ptr=(type*)(((((unsigned int)(rec))>=((VLARec*)(ptr))[-1].nAlloc) ? VLAExpand(ptr,(rec)) : (ptr))))
#define vla_malloc(ptr,type,initSize) (ptr=(type*)VLAMalloc(initSize,sizeof(type),5,0))
#define vla_calloc(ptr,type,initSize) (ptr=(type*)VLAMalloc(initSize,sizeof(type),5,1))
#define vla_free(ptr) {if(ptr) {VLAFree(ptr);ptr=NULL;}}
#define vla_set_size(ptr,type,size) {ptr=(type*)VLASetSize(ptr,size);}
#define vla_get_size(ptr) VLAGetSize2(ptr)

unsigned int VLAGetSize2(void *ptr);

#ifndef _os_memory_debug_on

void *VLAExpand(void *ptr,unsigned int rec); 
void *VLAMalloc(unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero); /*growfactor 1-10*/
void VLAFree(void *ptr);
void *VLASetSize(void *ptr,unsigned int newSize);

#else

#define VLAMalloc(a,b,c,d) _champVLAMalloc(__FILE__,__LINE__,a,b,c,d)
#define VLAExpand(a,b) _champVLAExpand(__FILE__,__LINE__,a,b)
#define VLAFree(a) _champVLAFree(__FILE__,__LINE__,a)
#define VLASetSize(a,b) _champVLASetSize(__FILE__,__LINE__,a,b)

void *_champVLAExpand(const char *file,int line,
                 void *ptr,unsigned int rec); 
void _champVLAFree(const char *file,int line,
              void *ptr);
void *_champVLASetSize(const char *file,int line,
                  void *ptr,unsigned int newSize);
void *_champVLAMalloc(const char *file,int line,
                 unsigned int initSize,unsigned int recSize,
                 unsigned int growFactor,int autoZero); 
#endif

#endif







