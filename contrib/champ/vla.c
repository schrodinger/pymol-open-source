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
/* This file can be compiled under C as a .c file, or under C++ as a .cc file*/

#include<stdio.h>

#include"vla.h"


#ifndef _os_memory_debug_on
void *VLAExpand(void *ptr,unsigned int rec)
#else
void *_VLAExpand(const char *file,int line,
                void *ptr,unsigned int rec)
#endif
{
  VLARec *vla;
  char *start,*stop;
  unsigned int soffset=0;
  vla = &(((VLARec*)ptr)[-1]);
  if(rec>=vla->nAlloc)
	 {
		if(vla->autoZero)
		  soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		vla->nAlloc = (rec*(vla->growFactor+10)/10)+1;
#ifndef _os_memory_debug_on
		vla=(void*)os_realloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec));
#else
		vla=(void*)OSMemoryRealloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec),file,line,_OSMemoryVLA);
#endif

		if(!vla)
		  {
			 printf("VLAExpand-ERR: realloc failed\n");
			 exit(EXIT_FAILURE);
		  }
		if(vla->autoZero)
		  {
			 start = ((char*)vla) + soffset;
			 stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
			 os_zero(start,stop);
		  }
	 }
  return((void*)&(vla[1]));
}

#ifndef _os_memory_debug_on
void *VLAMalloc(unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero)
#else
void *_VLAMalloc(const char *file,int line,unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero)
#endif
{
  VLARec *vla;
#ifndef _os_memory_debug_on
  if(autoZero)
    vla=(void*)os_calloc((initSize*recSize)+sizeof(VLARec));
  else
    vla=(void*)os_malloc((initSize*recSize)+sizeof(VLARec));
#else
  if(autoZero)
    vla=OSMemoryCalloc(1,(initSize*recSize)+sizeof(VLARec),file,line,_OSMemoryVLA);
  else
    vla=OSMemoryMalloc((initSize*recSize)+sizeof(VLARec),file,line,_OSMemoryVLA);
#endif

  if(!vla)
	 {
		printf("VLAMalloc-ERR: realloc failed\n");
		exit(EXIT_FAILURE);
	 }
  vla->nAlloc=initSize;
  vla->recSize=recSize;
  vla->growFactor=growFactor;
  vla->autoZero=autoZero;
  return((void*)&(vla[1]));
}

#ifndef _os_memory_debug_on
void VLAFree(void *ptr)
#else
void _VLAFree(const char *file,int line,void *ptr)
#endif
{
  VLARec *vla;
  vla = &(((VLARec*)ptr)[-1]);
#ifndef _os_memory_debug_on
  os_free(vla);
#else
  OSMemoryFree(vla,file,line,_OSMemoryVLA);
#endif
}

unsigned int VLAGetSize2(void *ptr)
{
  VLARec *vla;
  vla = &((VLARec*)ptr)[-1];
  return(vla->nAlloc);
}

#ifndef _os_memory_debug_on
void *VLASetSize(void *ptr,unsigned int newSize)
#else
void *_VLASetSize(const char *file,int line,void *ptr,unsigned int newSize)
#endif
{
  VLARec *vla;
  unsigned int orig_size = 0;
  char *start;
  char *stop;
  vla = &((VLARec*)ptr)[-1];
  if(vla->autoZero)
	 orig_size = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
  vla->nAlloc = newSize;
#ifndef _os_memory_debug_on
  vla=(void*)os_realloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec));
#else
  vla=OSMemoryRealloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec),file,line,_OSMemoryVLA);
#endif

  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed\n");
		exit(EXIT_FAILURE);
	 }
  if(vla->autoZero)
	 {
      start = ((char*)vla)+orig_size;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		if(start<stop)
		  os_zero(start,stop);
	 }
  return((void*)&(vla[1]));
}











