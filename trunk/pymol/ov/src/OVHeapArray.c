#include"OVHeapArray.h"
#ifndef OV_JENARIX
#include"ov_utility.h"

#define OVHeapArray_COMPLAIN
#define OVHeapArray_NULL_CHECK

#ifndef OVHeap_TRACKING

#define ov_array_malloc(I,size) OVHeap_Malloc(I,size)
#define ov_array_calloc(I,num,size) OVHeap_Calloc(I,num,size)
#define ov_array_realloc(I,ptr,size) OVHeap_Realloc(I,ptr,size)
#define ov_array_free(I,ptr) OVHeap_Free(I,ptr)

#else

#define ov_array_malloc(I,size) \
       _OVHeap_Malloc(I,size,file,line,OVHeap_ARRAY_TYPE)
#define ov_array_calloc(I,num,size) \
       _OVHeap_Calloc(I,num,size,file,line,OVHeap_ARRAY_TYPE)
#define ov_array_realloc(I,ptr,size) \
       _OVHeap_Realloc(I,ptr,size,file,line,OVHeap_ARRAY_TYPE)
#define ov_array_free(I,ptr) \
       _OVHeap_Free(I,ptr,file,line,OVHeap_ARRAY_TYPE)

/* TO DO */

#endif


void *_OVHeapArray_Check(void *ptr,ov_size index
#ifdef OVHeap_TRACKING
                         ,const char *file,int line
#endif
)
{ /* array-growth routine */
  _OVHeapArray *vla = &(((_OVHeapArray*)ptr)[-1]);
  if(index>=vla->size) {
    ov_size new_size = (index+(index>>1)+1); 
    _OVHeapArray *new_vla;
    new_vla = (_OVHeapArray*)ov_array_realloc(vla->heap,vla,
                                              sizeof(_OVHeapArray)+(vla->unit_size*new_size));
    if(!new_vla) {
#ifdef OVHeapArray_COMPLAIN
      fprintf(stderr,"_OVHeapArray_Check-Error: realloc failed\n");
#endif
    } else {
      vla = new_vla;
      if(vla->auto_zero) {
        char *start = ((char*)vla) + sizeof(_OVHeapArray)+(vla->unit_size*vla->size);
        char *stop = ((char*)vla) + sizeof(_OVHeapArray)+(vla->unit_size*new_size);
        ov_utility_zero_range(start,stop);
      }
      vla->size = new_size;
    }
  }
  return((void*)&(vla[1]));
}

void *_OVHeapArray_Alloc(OVHeap *heap,ov_size unit_size,
                         ov_size size,int zero
#ifdef OVHeap_TRACKING
                         ,const char *file,int line
#endif
)
{
  _OVHeapArray *vla;
  if(zero) {
    vla = ov_array_calloc(heap,1,sizeof(_OVHeapArray)+(unit_size*size));
  } else {
    vla = ov_array_malloc(heap,sizeof(_OVHeapArray)+(unit_size*size));
  }
  if(!vla) {
#ifdef OVHeapArray_COMPLAIN
    fprintf(stderr,"_OVHeapArray: realloc failed\n");
#endif
    return OV_NULL;
  } else {
    vla->heap=heap;
    vla->size=size;
    vla->unit_size=unit_size;
    vla->auto_zero = zero;
    return((void*)&(vla[1]));
  }
}

void  _OVHeapArray_Free(void *ptr
#ifdef OVHeap_TRACKING
                        ,const char *file,int line
#endif
                        )
{
  _OVHeapArray *vla;
#ifdef OVHeapArray_NULL_CHECK
  if(!ptr) {
#ifdef OVHeapArray_COMPLAIN
		fprintf(stderr,"_OVHeapArray_Free-Error: tried to free NULL pointer!\n");
#endif
  } else
#endif
 {
    vla = &(((_OVHeapArray*)ptr)[-1]);
    ov_array_free(vla->heap,vla);
  }
}

ov_size OVHeapArray_GetSize(void *ptr
#ifdef OVHeap_TRACKING
                            ,const char *file,int line
#endif
                            )
{
  _OVHeapArray *vla;
  vla = &((_OVHeapArray*)ptr)[-1];
  return(vla->size);
}

/*
void *VLANewCopy(void *ptr)
{
  _OVHeapArray *vla,*new_vla;
  unsigned int size;
  vla = &((_OVHeapArray*)ptr)[-1];
  size = (vla->unit_size*vla->size)+sizeof(_OVHeapArray);
  new_vla=(void*)mmalloc(size);
  if(!new_vla)
	 {
		printf("VLACopy-ERR: mmalloc failed\n");
		exit(EXIT_FAILURE);
	 }
  else
    {
      memcpy(new_vla,vla,size);
    }
  return((void*)&(new_vla[1]));
}
*/

void *_OVHeapArray_SetSize(void *ptr, ov_size new_size
#ifdef OVHeap_TRACKING
                        ,const char *file,int line
#endif
                        )
{
  _OVHeapArray *vla,*new_vla;
  vla = &((_OVHeapArray*)ptr)[-1];
  new_vla=(void*)ov_array_realloc(vla->heap,vla,(vla->unit_size*new_size)+sizeof(_OVHeapArray));
  if(!new_vla)	 {
#ifdef OVHeapArray_COMPLAIN      
		fprintf(stderr,"VLASetSize-ERR: realloc failed.\n");
#endif
  } else {
    vla = new_vla;
    if(new_size>vla->size && vla->auto_zero) {
      char *start = ((char*)vla) + sizeof(_OVHeapArray)+(vla->unit_size*vla->size);
      char *stop = ((char*)vla) + sizeof(_OVHeapArray)+(vla->unit_size*new_size);
      ov_utility_zero_range(start,stop);
    }
    vla->size = new_size;
  }
  return((void*)&(vla[1]));
}

#endif
