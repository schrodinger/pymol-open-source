#ifndef _H_OVHeapArray
#define _H_OVHeapArray

#include "OVHeap.h"

/* private implementation */


typedef struct {
  ov_port_size_t size, rec_size;
  OVHeap *heap;
  ov_boolean auto_zero;
} _OVHeapArray;


void *_OVHeapArray_Alloc(OVHeap *heap,ov_port_size_t rec_size,
                         ov_port_size_t size,int zero
#ifdef OVHeap_TRACKING
                         ,const char *file,int line
#endif
                         );


void *_OVHeapArray_Check(void *ptr,ov_port_size_t index
#ifdef OVHeap_TRACKING
                         ,const char *file,int line
#endif
                         );


void  _OVHeapArray_Free(void *ptr
#ifdef OVHeap_TRACKING
                        ,const char *file,int line
#endif
                        );


void *_OVHeapArray_SetSize(void *ptr, ov_port_size_t new_size
#ifdef OVHeap_TRACKING
                        ,const char *file,int line
#endif
                        );

ov_port_size_t OVHeapArray_GetSize(void *ptr);


/* begin interface */

#ifndef OVHeap_TRACKING

#define OVHeapArray_MALLOC(heap,type,size)  (type*)_OVHeapArray_Alloc(heap,sizeof(type),\
                                                                   size,0)
#define OVHeapArray_CALLOC(heap,type,size)  (type*)_OVHeapArray_Alloc(heap,sizeof(type),\
                                                                   size,1)

#define OVHeapArray_CHECK(ptr,type,index) (((index)>=(((_OVHeapArray*)(ptr))[-1].size)) ? \
    (((index)<((_OVHeapArray*)((ptr)=\
        (type*)_OVHeapArray_Check((void*)(ptr),index)))[-1].size))\
           : OV_TRUE)

#define OVHeapArray_SET_SIZE(ptr,type,new_size) \
 (((new_size)== \
      ((_OVHeapArray*)((ptr)= \
                       (type*)_OVHeapArray_SetSize((void*)(ptr),(new_size))))[-1].size))


#define OVHeapArray_FREE(ptr) _OVHeapArray_Free((void*)(ptr))


#else

#define OVHeapArray_MALLOC(heap,type,size)  (type*)_OVHeapArray_Alloc(heap,sizeof(type),\
                                                                   size,0,\
                                                                   __FILE__,__LINE__)
#define OVHeapArray_CALLOC(heap,type,size)  (type*)_OVHeapArray_Alloc(heap,sizeof(type),\
                                                                   size,1,\
                                                                   __FILE__,__LINE__)

#define OVHeapArray_CHECK(ptr,type,index) (((index)>=(((_OVHeapArray*)(ptr))[-1].size)) ? \
    (((index)<((_OVHeapArray*)((ptr)=\
        (type*)_OVHeapArray_Check((void*)(ptr),index,__FILE__,__LINE__)))[-1].size))\
           : OV_TRUE)

#define OVHeapArray_SET_SIZE(ptr,type,new_size) \
 (((new_size)== \
      ((_OVHeapArray*)((ptr)= \
                       (type*)_OVHeapArray_SetSize((void*)(ptr),(new_size),__FILE__,__LINE__)))[-1].size))


#define OVHeapArray_FREE(ptr) _OVHeapArray_Free((void*)(ptr),__FILE__,__LINE__)

#endif

#define OVHeapArray_FREE_AUTO_NULL(ptr) { if(ptr) { OVHeapArray_FREE(ptr); (ptr)=OV_NULL;}}

#endif







