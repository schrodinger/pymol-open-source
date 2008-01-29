#ifndef  _H_OVHeap
#define _H_OVHeap

#include "ov_port.h"
#include "ov_types.h"

#ifdef OV_JENARIX

/* NEW Jenarix-based OVHeap wrapper */

#ifdef XX_HEAP_TRACKER
#ifndef OVHeap_TRACKING
#define OVHeap_TRACKING
#endif
#endif

#define OVHeap_DUMP_FILES_TOO      OV_HEAP_DUMP_FILES_TOO
#define OVHeap_DUMP_NO_ADDRESSES   OV_HEAP_DUMP_NO_ADDRESSES

typedef int *OVHeap;
#define OVHeap_New()   ((OVHeap*)1)
#define OVHeap_Del(x) 
#define OVHeap_Malloc(I,size) OV_HEAP_MALLOC_RAW_VOID(size)
#define OVHeap_Calloc(I,num,size) OV_HEAP_CALLOC_RAW_VOID((ov_size)(num)*(ov_size)(size))
#define OVHeap_Realloc(I,ptr,size) OV_HEAP_REALLOC_RAW_VOID(ptr,size)
#define OVHeap_Free(I,ptr) OV_HEAP_FREE_RAW(ptr)
#define OVHeap_Dump(I,flags) ov_heap_dump(flags)
#define OVHeap_Usage(I) ov_heap_usage()

/* convenience macros */

#define OVHeap_ALLOC(I,type) OVHeap_Calloc(I,1,sizeof(type))

#define OVHeap_MALLOC(I,type,num) OVHeap_Malloc(I,sizeof(type)*(num))
#define OVHeap_CALLOC(I,type,num) OVHeap_Calloc(I,num,sizeof(type))
#define OVHeap_REALLOC(I,ptr,type,num) OVHeap_Realloc(I,ptr,sizeof(type)*(num))
#define OVHeap_FREE_AUTO_NULL(I,ptr) { if(ptr) {OVHeap_Free(I,ptr); ptr = NULL;}}

#else
/* OLD proven OVHeap implementation */

#ifdef OVHeap_TRACKING

struct _OVHeap;  typedef struct _OVHeap OVHeap;

/* what kinds of blocks do we keep track of? */

#define OVHeap_BLOCK_TYPE 0
#define OVHeap_ARRAY_TYPE 1

/* how do we want the dump printed? */

#define OVHeap_DUMP_FILES_TOO      0x1
#define OVHeap_DUMP_NO_ADDRESSES   0x2

/* constructor/destructor */

#define OVHeap_New()   _OVHeap_New()
#define OVHeap_Del(x)  _OVHeap_Del(x)

/* implementation -- available for use by derivative memory managers */

OVHeap *_OVHeap_New(void);
void    _OVHeap_Del(OVHeap *I);

void   *_OVHeap_Malloc(OVHeap *I,ov_size size,const char *file,
                     int line,int type);
void   *_OVHeap_Calloc(OVHeap *I,ov_size num,ov_size size,
                     const char *file,int line,int type);
void   *_OVHeap_Realloc(OVHeap *I,void *ptr,ov_size size,
                      const char *file,int line,int type);
void    _OVHeap_Free(OVHeap *I,void *ptr,const char *file,int line,int type);

/* macros for tracking */

#define OVHeap_Malloc(I,size) \
       _OVHeap_Malloc(I,size,__FILE__,__LINE__,OVHeap_BLOCK_TYPE)
#define OVHeap_Calloc(I,num,size) \
       _OVHeap_Calloc(I,num,size,__FILE__,__LINE__,OVHeap_BLOCK_TYPE)
#define OVHeap_Realloc(I,ptr,size) \
       _OVHeap_Realloc(I,ptr,size,__FILE__,__LINE__,OVHeap_BLOCK_TYPE)
#define OVHeap_Free(I,ptr) \
       _OVHeap_Free(I,ptr,__FILE__,__LINE__,OVHeap_BLOCK_TYPE)


void OVHeap_Dump(OVHeap *I,ov_uint32 flags);
int  OVHeap_Usage(OVHeap *I);

#else

typedef int *OVHeap;

/* macros without tracking */

#define OVHeap_New()            ((void*)1)
#define OVHeap_Del(I)      

#define OVHeap_Malloc(I,size)    malloc(size)
#define OVHeap_Calloc(I,num,size)  calloc(num,size)
#define OVHeap_Realloc(I,ptr,size) realloc(ptr,size)
#define OVHeap_Free(I,ptr)      free(ptr)

#define OVHeap_Dump(I,log_to_files)

#endif

/* convenience macros */

#define OVHeap_ALLOC(I,type) OVHeap_Calloc(I,1,sizeof(type))

#define OVHeap_MALLOC(I,type,num) OVHeap_Malloc(I,sizeof(type)*(num))
#define OVHeap_CALLOC(I,type,num) OVHeap_Calloc(I,num,sizeof(type))
#define OVHeap_REALLOC(I,ptr,type,num) OVHeap_Realloc(I,ptr,sizeof(type)*(num))
#define OVHeap_FREE_AUTO_NULL(I,ptr) { if(ptr) {OVHeap_Free(I,ptr); ptr = NULL;}}

#endif
#endif







