

#include "OVOneToAny.h"
#include "OVHeapArray.h"
#include "ov_utility.h"

#define HASH(value,mask) (((value^(value>>24))^((value>>8)^(value>>16)))&mask)

/* FYI: "up" stands for UniquePair -- a precursor to OneToAny */

typedef struct {
  int active;
  ov_word forward_value, reverse_value;
  ov_size forward_next;
} up_element;

struct _OVOneToAny {
  OVHeap *heap;
  ov_uword mask;
  ov_size size,n_inactive;
  ov_word next_inactive;
  up_element *elem;
  ov_word *forward; 
};


OVstatus OVOneToAny_Init(OVOneToAny *up,OVHeap *heap)
{
  ov_utility_zero_range(up,up+1);
  up->heap = heap;
  return_OVstatus_SUCCESS;
}

OVOneToAny *OVOneToAny_New(OVHeap *heap)
{
  OVOneToAny *up;
  up = OVHeap_ALLOC(heap,OVOneToAny);
  up->heap = heap;
  return up;
}


void OVOneToAny_Purge(OVOneToAny *up)
{
  if(up) {    
    OVHeapArray_FREE_AUTO_NULL(up->elem);
    OVHeap_FREE_AUTO_NULL(up->heap,up->forward);
  }
}

void OVOneToAny_Del(OVOneToAny *up)
{
  if(up) {
    OVOneToAny_Purge(up);
    OVHeap_FREE_AUTO_NULL(up->heap,up);
  }  
}

void OVOneToAny_Reset(OVOneToAny *up) 
{
  OVOneToAny_Purge(up);
  OVOneToAny_Init(up,up->heap);
}

void OVOneToAny_Dump(OVOneToAny *up)
{
  ov_uword a;
  ov_boolean empty = OV_TRUE;
  if(up&&up->mask) {
    for(a=0;a<=up->mask;a++) {
      if(up->forward[a]) {
        fprintf(stderr,
" OVOneToAny_Dump: Hashes forward[0x%02x]->%d\n",
                (unsigned int)a,(int)up->forward[a]);
        empty = OV_FALSE;
      }
    }

    for(a=0;a<up->size;a++)
      if(up->elem[a].active) {
        fprintf(stderr,
" OVOneToAny_Dump: Elements %d:    %d (->%d)    %d \n",
                a+1,
                (int)up->elem[a].forward_value,
                (int)up->elem[a].forward_next,
                (int)up->elem[a].reverse_value                );
        empty=OV_FALSE;
      }
  }
  if(empty) {
    fprintf(stderr," OVOneToAny_Dump: Empty.\n");
  }
}


static void Reload(OVOneToAny *up)
{ /* assumes hash tables are clean and initialized to zero */
#ifdef DEBUG_UP
  fprintf(stderr,"Reload-Debug: entered\n");
#endif
  register ov_uword mask = up->mask;

  if(up->elem&&mask) {
    { 
      register up_element *elem = up->elem;
      register ov_uword a;
      for(a=0;a<up->size;a++) {
        if(elem->active) {
          elem->forward_next = 0; /* 0 is the sentinel for end of list */
        }
        elem++;
      }
    }
    
    { 
      register ov_uword a;
      register ov_word *forward = up->forward;
      register up_element *elem = up->elem;
      {
        register ov_word fwd;
        register ov_word fwd_val;
        for(a=0;a<up->size;a++) {
          if(elem->active) {
            fwd_val = elem->forward_value;
            fwd = HASH(fwd_val, mask);
            elem->forward_next = forward[fwd];
            forward[fwd] = a+1; /* NOTE: 1 based indices */
          }
          elem++;
        }
      }
    }
  }

#ifdef DEBUG_UP
  {
    ov_uword a;
    for(a=0;a<=up->mask;a++) {
      fprintf(stderr,"Reload-Debug: forward[%d]=%d\n",
              a,up->forward[a],);
    }
  }
#endif

}

OVreturn_word OVOneToAny_GetKey(OVOneToAny *up,ov_word forward_value)
{
#ifdef DEBUG_OVOneToAny
  fprintf(stderr,"OVOneToAnyGetKey-Debug: %d\n",forward_value);
#endif
  if(!up) {
    OVreturn_word result = { OVstatus_NULL_PTR };
    return result;
  } else {
    register ov_uword mask = up->mask;
    if(mask) {
      ov_word hash = HASH(forward_value,mask);
      register up_element *elem = up->elem;
      register ov_word index = up->forward[hash];
      register up_element *cur_elem = elem+(index-1);
#ifdef DEBUG_OVOneToAny
      fprintf(stderr,"OVOneToAnyGetKey-Debug: hash %d index %d\n",hash, index);
#endif

      while(index) {
#ifdef DEBUG_OVOneToAny
        fprintf(stderr,"OVOneToAnyGetKey-Debug: index %d forward_value %d\n", index, cur_elem->forward_value);
#endif

        if(cur_elem->forward_value==forward_value) {
          OVreturn_word result = { OVstatus_SUCCESS };
          result.word = cur_elem->reverse_value;
          return result;
        }
        index = cur_elem->forward_next;
        cur_elem = elem+(index-1);
      }
    }
    {
      OVreturn_word result = { OVstatus_NOT_FOUND };
      return result;
    }
  }
}


static OVstatus Recondition(OVOneToAny *up,ov_uword size,int force)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    ov_uword mask = up->mask;
#ifdef DEBUG_UP
    fprintf(stderr,"Recondition-Debug: entered for size %d.\n",size);
#endif
    if((size>mask)||((size<<2)<mask)||force) {
      
      while((size<<2)<mask) {
        mask = mask>>1;
        if(mask<2) break;
      }
      
      while(size>mask) {
        mask=(mask<<1)+1;
      }
      

#ifdef DEBUG_UP
      fprintf(stderr,"Recondition-Debug: mask %d\n",mask);
#endif
      {
        if(!up->elem) {
          up->elem    = OVHeapArray_CALLOC(up->heap,up_element,size);
          if(!up->elem) {
            return_OVstatus_OUT_OF_MEMORY;
          }
        } 
        if(mask!=up->mask) {
          ov_word *tmp_forward = OVHeap_CALLOC(up->heap,ov_word,mask+1);
          if(!tmp_forward) { /* validate */
            OVHeap_FREE_AUTO_NULL(up->heap,tmp_forward);
            /* being unable to condition is not an error */          
          } else {
            /* impossible to fail after here... */
            OVHeap_FREE_AUTO_NULL(up->heap,up->forward);
            up->forward = tmp_forward;
            up->mask = mask;
          } 
        } else {
          ov_utility_zero_range(up->forward,up->forward+(up->mask+1));
        }
        Reload(up);
      }
    }
  }
  return_OVstatus_SUCCESS;
}


OVstatus OVOneToAny_Pack(OVOneToAny *up) 
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    if(up->n_inactive&&up->elem) {
      ov_uword new_size = 0;
      up_element *src=up->elem,*dst=up->elem;
      ov_uword a;

      for(a=0;a<up->size;a++) {
        if(src->active) {
          if(src>dst) {
            *dst = *src;
          } 
          dst++;
          new_size++;
        }
        src++;
      }
      up->n_inactive = 0;
      up->next_inactive = 0;
      if(new_size<up->size) {
        if(!OVHeapArray_SET_SIZE(up->elem, up_element, new_size))
          ov_utility_zero_range(up->elem+new_size,up->elem+up->size);
      }
      up->size = new_size;
      return Recondition(up,new_size,OV_TRUE);
    }
    return_OVstatus_SUCCESS;
  }
}

 
OVreturn_size OVOneToAny_GetSize(OVOneToAny *up)
{
  if(!up) {
    OVreturn_size result = { OVstatus_NULL_PTR };
    return result;
  } else {
    OVreturn_size result = { OVstatus_SUCCESS };
    result.size = up->size - up->n_inactive;
    return result;
  }
}


OVstatus OVOneToAny_DelKey(OVOneToAny *up,ov_word forward_value)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    register ov_word mask = up->mask;
    if(mask) {
      register ov_word fwd_hash = HASH(forward_value,mask);
      register ov_word fwd = up->forward[fwd_hash];
      if(!fwd) {
        return_OVstatus_NOT_FOUND;
      } else {
        register up_element *fwd_elem = NULL;
        register up_element *elem = up->elem;
        register ov_word fwd_last = 0;
        
        while(fwd) {
          fwd_elem = elem+(fwd-1);
          if(fwd_elem->forward_value == forward_value) 
            break;
          fwd_last = fwd;
          fwd = fwd_elem->forward_next;
        }
        
        if(fwd_elem) {
          if(fwd) {
            
            /* excise elements */
            
            if(fwd_last)
              up->elem[fwd_last-1].forward_next = fwd_elem->forward_next;
            else
              up->forward[fwd_hash] = fwd_elem->forward_next;
            
            /* store as inactive */
            
            fwd_elem->active = OV_FALSE;
            fwd_elem->forward_next = up->next_inactive;
            up->next_inactive = fwd;
            up->n_inactive++;
            if(up->n_inactive>(up->size>>1)) /* over half of bits are inactive */
              OVOneToAny_Pack(up);          
            return_OVstatus_SUCCESS;
          }
        } 
      }
    }
    return_OVstatus_NOT_FOUND;
  }
}

void OVOneToAny_Stats(OVOneToAny *up)
{
  if(up&&up->mask) {
    int max_len=0;
    ov_uword a;
    for(a=0;a<up->mask;a++) {
      {
        ov_word index = up->forward[a];
        up_element *elem = up->elem;
        int cnt = 0;
        if(index) {
          up_element *cur_elem;
          while(index) {
            cur_elem = elem+(index-1);
            index = cur_elem->forward_next;
            cnt++;
          }
          if(cnt>max_len)
            max_len=cnt;
        }
      }
     
    }
    fprintf(stderr," OVOneToAny_Stats: MaxLen=%d ",(int)max_len);
    fprintf(stderr,"active=%d n_inactive=%d ",(int)(up->size-up->n_inactive),
            (int)up->n_inactive);
    fprintf(stderr,"mask=0x%x n_alloc=%lu\n",(unsigned int)up->mask,
            (unsigned long)OVHeapArray_GET_SIZE(up->elem));
  }
}


OVstatus OVOneToAny_SetKey(OVOneToAny *up, ov_word forward_value, ov_word reverse_value)
{
#ifdef DEBUG_OVOneToAny
      fprintf(stderr,"OVOneToAnySetKey-Debug: %d,%d\n",forward_value,reverse_value);
#endif
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    register ov_word mask = up->mask;
    register ov_word fwd_hash = HASH(forward_value,mask);
    register up_element *fwd_elem = NULL;
    register up_element *rev_elem = NULL;
    register ov_word fwd;
    if(!mask) {
      fwd = 0;
    } else {
    
      fwd = up->forward[fwd_hash];
      
#ifdef DEBUG_OVOneToAny
      fprintf(stderr,"OVOneToAnySet-Debug: fwd_hash %d mask %d size %d\n",
              fwd_hash,up->mask,up->size);
      fprintf(stderr,"OVOneToAnySet-Debug: before search fwd %d \n",fwd);
#endif
      
      { /* find elements if they exist, and detect erroneous conditions */
        
        register up_element *elem = up->elem;
        
#ifdef DEBUG_OVOneToAny
        {
          int a;
          for(a=0;a<up->size;a++) {
            fprintf(stderr,"OVOneToAnySet-Debug: on entry %d forward_next: %d\n",
                 a+1,elem[a].forward_next);
        }
        
        for(a=0;a<=up->mask;a++) {
          fprintf(stderr,
"OVOneToAnySet-Debug: on entry %d forward[%d]=%d:\n",
                 a,a,up->forward[a]);
        }
      }
#endif
      
        while(fwd) {
          fwd_elem = elem+(fwd-1);
          if(fwd_elem->forward_value==forward_value)
            break;
          fwd = fwd_elem->forward_next;
        }
      }
      
      if(fwd) {
        return_OVstatus_DUPLICATE;
      }
    }
    if(!(fwd)) {  /* new entry */
      ov_size new_index;
      /* new pair */
#ifdef DEBUG_OVOneToAny
      fprintf(stderr,"OVOneToAnySet-Debug: New pair.\n");
#endif
      if(up->n_inactive) {
        new_index = up->next_inactive;
        up->next_inactive = up->elem[new_index-1].forward_next;
        up->n_inactive--;
      } else {
        if(up->elem&&(!OVHeapArray_CHECK(up->elem, up_element, up->size))) {
          return_OVstatus_OUT_OF_MEMORY;
        } else {
          OVstatus result;
          if(OVreturn_IS_ERROR(result = Recondition(up,up->size+1,OV_FALSE))) {
            return result;
          } else {
            /* guaranteed to succeed past this point, so we can increase size */
            new_index = ++up->size;
          }
        }
      }
      {
        up_element *elem = up->elem + (new_index-1);
        elem->forward_value = forward_value;
        elem->reverse_value = reverse_value;
        elem->active = OV_TRUE;
        
        /* regenerate new hashes */
        mask = up->mask;
        fwd_hash = HASH(forward_value, mask);
        
        {
          ov_word *forward_start_index = up->forward + fwd_hash;
          
          elem->forward_next = *forward_start_index;
          *forward_start_index = new_index; /* note the +1 offset */
        }
        
#ifdef DEBUG_OVOneToAny
        {
          int a;
          for(a=0;a<=up->mask;a++) {
            fprintf(stderr,"OVOneToAnySet-Debug: forward[%d]=%d\n",
                    a,up->forward[a]);
          }
        }
#endif
        
      }
      
    } else if(fwd_elem!=rev_elem) { 
      return_OVstatus_MISMATCH;
    } else { 
      return_OVstatus_NO_EFFECT;
      /* exists and matched, so do nothing */
    }
  }
  return_OVstatus_SUCCESS; 
}


