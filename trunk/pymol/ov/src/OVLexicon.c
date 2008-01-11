
#include "OVLexicon.h"
#include "OVOneToOne.h"
#include "OVHeapArray.h"
#include "OVreturns.h"

/* should only be accessed by special methods */

typedef struct {
  ov_word offset;
  ov_word next; /* NOTE: 1-based, 0 is the sentinel */
  ov_word ref_cnt;
  ov_uword hash;
  ov_size size;
} lex_entry;

struct _OVLexicon {
  OVHeap *heap;
  OVOneToOne *up; /* maps hash_key to index of first zero-based entry */
  lex_entry *entry; /* maintained as a 1-based index */
  ov_uword n_entry, n_active;
  ov_char8 *data;
  ov_uword data_size;
  ov_uword data_unused;
  ov_word free_index; /* NOTE: 1-based, 0 is the sentinel */
};

ov_uword OVLexicon_GetNActive(OVLexicon *uk)
{
  return uk->n_active;
}

OVLexicon *OVLexicon_New(OVHeap *heap)
{
  OVLexicon *I = NULL;
  if(heap) {
    I = OVHeap_ALLOC(heap,OVLexicon);
    if(I) {
      I->heap = heap;
      I->up = OVOneToOne_New(heap);
      if(!(I->up)) {
        OVLexicon_DEL_AUTO_NULL(I);
      }
    }
  }
  return I;
}

void OVLexicon_Del(OVLexicon *I)
{
  if(I) {
    OVOneToOne_DEL_AUTO_NULL(I->up);
    if(I->entry) {
      I->entry++; /* allow for 1-based indexing */
      OVHeapArray_FREE_AUTO_NULL(I->entry);
    }
    OVHeapArray_FREE_AUTO_NULL(I->data);
    OVHeap_Free(I->heap,I);
  }
}

static OVstatus OVLexicon_CheckStorage(OVLexicon *uk, ov_word entry_size, ov_size data_size)
{
  if(!uk->entry) { /* auto zero to initialize ref_cnt fields */
    uk->entry = OVHeapArray_CALLOC(uk->heap, lex_entry, entry_size);
    if(!uk->entry) {
      return_OVstatus_OUT_OF_MEMORY;
    }
    uk->entry--; /* allow for 1-based indexing */

  } else {
    uk->entry++; /* allow for 1-based indexing */
    if(!OVHeapArray_CHECK(uk->entry, lex_entry, (ov_size)(entry_size-1))) {
      return_OVstatus_OUT_OF_MEMORY;
    }
    uk->entry--; /* allow for 1-based indexing */
  }
    
  if(!uk->data) {
    uk->data = OVHeapArray_MALLOC(uk->heap, ov_char8, data_size);
    if(!uk->data) {
      return_OVstatus_OUT_OF_MEMORY;
    }
  } else if(!OVHeapArray_CHECK(uk->data, ov_char8, data_size-1)) {
    return_OVstatus_OUT_OF_MEMORY;
  }
  return_OVstatus_SUCCESS;
}

static ov_word _GetCStringHash(ov_uchar8 *str)
{
  register ov_uchar8 *p = str;
  register ov_word x;
  register ov_size len = 0;
  register ov_uchar8 c;

  x = *p << 7;
  while ( (c=*(p++)) ) {
#if 0
    x = (1000003*x) + c;    /* PYTHON (time: G5 = 3.2, P3 = 19.3, P4=11.1)*/
#endif
#if 0
    x = (x << 6) + (x << 16) - x + c; /* aho (G5 = 3.6 sec, P3 = 15.2?!, P4=3.6) */
#endif
#if 1
    x = (x << 5) + x + c;   /*  djb2 (G5 = 2.8, P3 = 18.7, P4=2.9) FASTEST OVERALL */
#endif
    len++;
  }
  x ^= len;
  return x;
}

#if 0
ov_word OVLexicon_GetCStringHash(ov_char8 *str) {
  return _GetCStringHash((ov_uchar8*)str);
}
#endif

OVstatus OVLexicon_Pack(OVLexicon *uk)
{
  
  if(uk->entry && uk->data && uk->n_entry && uk->data_unused) {
    ov_size new_count = 0;
    ov_size new_size = 0;

    { /* compute storage requirements */
      ov_size a;
      ov_size n_entry = uk->n_entry;
      lex_entry *cur_entry = uk->entry + 1; /* NOTE: 1-based array */
      for(a=0;a<n_entry;a++) {
        if(cur_entry->ref_cnt>0) {
          new_size+=cur_entry->size;
          new_count++;
        }
       cur_entry++;
      }
    }
    if((!new_count)&&(!new_size)) { /* if lexicon is completely empty, then purge both */
      uk->entry++; /* alloc for 1-based indexing */
      OVHeapArray_FREE_AUTO_NULL(uk->entry);
      OVHeapArray_FREE_AUTO_NULL(uk->data);
      OVOneToOne_Reset(uk->up);
      uk->n_entry = 0;
      uk->n_active = 0;
      uk->data_unused = 0;
      uk->data_size = 0;
      uk->free_index = 0;
    } else { /* otherwise, pack the string fields, and track the free entries */
      OVstatus status;
      ov_char8 *old_data = uk->data;
      uk->data = NULL;
      
      if(OVreturn_IS_ERROR(status = OVLexicon_CheckStorage(uk,uk->n_entry,new_size))) {
        uk->data = old_data;
        return status;
      }
      
      { /* copy data */
        register ov_word a;
        register ov_size n_entry = uk->n_entry, new_size = 0;
        register lex_entry *cur_entry = uk->entry + 1; /* NOTE: 1-based array */
        register ov_char8 *data = old_data;
        register ov_char8 *new_data = uk->data;
        register ov_word free_index = 0;
        register ov_size entry_size;
        
        for(a=1;a<=(ov_word)n_entry;a++) {
          if(cur_entry->ref_cnt>0) {
            entry_size = cur_entry->size;
            memcpy(new_data, data+cur_entry->offset, entry_size);
            cur_entry->offset = new_size;
            new_size += entry_size;
            new_data += entry_size;
          } else { 
            /*  remember for later */
            cur_entry->next = free_index;
            cur_entry->ref_cnt = 0;
            free_index = a;
          }
          cur_entry++;
        }
        
        OVHeapArray_FREE(old_data);
        uk->data_unused = 0;
        uk->data_size = new_size;
        uk->free_index = free_index;

      }
    }
  }
  return_OVstatus_SUCCESS;
}

OVstatus OVLexicon_DecRef(OVLexicon *uk, ov_word id)
{
  if((!uk->entry)||(id<1)||(id>(ov_word)uk->n_entry)) { /* range checking */
    return_OVstatus_NOT_FOUND;
  } else {
    register lex_entry *cur_entry = uk->entry + id;
    ov_word ref_cnt = (--cur_entry->ref_cnt);
    if(ref_cnt<0) {
      return_OVstatus_INVALID_REF_CNT;
    } else if(!ref_cnt) {
      OVreturn_word result = OVOneToOne_GetForward(uk->up,cur_entry->hash);
      if(OVreturn_IS_OK(result)) {
        register ov_word index = result.word;
        if(index!=id) {
          register lex_entry *entry = uk->entry;
          register ov_word next;
          while(index&&(next=entry[index].next)!=id)
            index = next;
          if(index) 
            entry[index].next = entry[id].next; /* excise from list */
        } else {
          /* remove entry from OneToOne */
          OVOneToOne_DelReverse(uk->up,id);
          /* if non-terminal, then add next entry to OneToOne */
          if(cur_entry->next) { 
            OVOneToOne_Set(uk->up,cur_entry->hash,cur_entry->next); /* NOTE: no error checking performed! */
          }
        }
      }
      uk->data_unused+=cur_entry->size;
      uk->n_active--;
      if(uk->data_unused>=(uk->data_size>>1)) /* if 50% unutilized, then pack */
        OVLexicon_Pack(uk);
    }
    return_OVstatus_SUCCESS;
  }
}

OVstatus OVLexicon_IncRef(OVLexicon *uk, ov_word id)
{
  if((!uk->entry)||(id<1)||(id>(ov_word)uk->n_entry)) { /* range checking */
    return_OVstatus_NOT_FOUND;
  } else {
    register lex_entry *entry = uk->entry + id;
    ov_word ref_cnt = (++entry->ref_cnt);
    if(ref_cnt<2) { /* was reference count zero or less? */
      /* safety precauations */
      entry->ref_cnt = 0; 
      entry->size = 0;
      entry->offset = 0; 
      return_OVstatus_INVALID_REF_CNT;
    } else {
      return_OVstatus_SUCCESS;
    }
  }
}


OVreturn_word OVLexicon_GetFromCString(OVLexicon *uk,ov_char8 *str)
{
  ov_word hash = _GetCStringHash((ov_uchar8*)str);
  OVreturn_word search = OVOneToOne_GetForward(uk->up, hash);
  register ov_word index = 0;
  ov_word cur_index = 0;
  ov_char8 *c;

  if(OVreturn_IS_OK(search)) {
    /* the hash key has been registered, so follow the chain
     * looking for a match... */
    register ov_char8 *data = uk->data;
    register lex_entry *entry = uk->entry,*entry_ptr;
    cur_index = (index = search.word);
    while(index) { /* found */
      c = data + (entry_ptr = (entry + index))->offset;
      if(strcmp(c,str)!=0) {/* verify match */
        index = entry_ptr->next; /* not a match? keep looking */
      } else {
        break;
      }
    }
  }
  if(!index) { 
    /* no match was found (or hash is new) so add this string to the lexicon */
    ov_size st_size = strlen(str) + 1;
    register lex_entry *entry, *entry_ptr, *cur_entry_ptr;
    OVstatus status;
    {
      ov_size new_size = uk->data_size + st_size;
      ov_size new_n_entry = uk->n_entry;
      if(!uk->free_index)
        new_n_entry++;
      /* allocate storage as necessary */
      if(OVreturn_IS_ERROR(status = OVLexicon_CheckStorage(uk,new_n_entry,new_size))) {
        OVreturn_word result;
        result.status = status.status;
        result.word = 0;
        return result;
      }
    }
    /* WARNING: uk->entry, uk->data_size, and uk->data may have changed! */

    if(!uk->free_index)
      index = (++uk->n_entry); /* NOTE: 1-based indices */
    else {
      index = uk->free_index;
      uk->free_index = uk->entry[index].next;
    }
    uk->n_active++;

    if(!cur_index) {
      /* if the hash key was new, add it to the OneToOne, and setup entry */      
      if(OVreturn_IS_ERROR(status = OVOneToOne_Set(uk->up, hash, index))) {
        OVreturn_word result;
        uk->entry[index].next = uk->free_index; /* record this as a free entry */
        uk->free_index = index; 
        uk->n_active--;
        result.status = status.status;
        result.word = 0;
        return result;
      }
      entry_ptr = uk->entry + index;
      entry_ptr->next = 0;

    } else { 
      /* otherwise, simply add this entry on to the list in position 2 */

      cur_entry_ptr = (entry = uk->entry) + cur_index;
      entry_ptr = entry + index;
      entry_ptr->next = cur_entry_ptr->next;
      cur_entry_ptr->next = index;
    }
    /* increase the counts and sizes to accomodate this new entry & copy info */
    entry_ptr->hash = hash;
    entry_ptr->size = st_size;
    entry_ptr->offset = uk->data_size;
    entry_ptr->ref_cnt++; /* increase ref_cnt */
    strcpy(uk->data + uk->data_size, str);
    uk->data_size += st_size;
  } else {
    lex_entry *entry_ptr;
    entry_ptr = uk->entry + index; 
    entry_ptr->ref_cnt++; /* increase ref_cnt */
  }
  {
    OVreturn_word result = { OVstatus_SUCCESS };
    result.word = index;
    return result;
  }
}



OVreturn_word OVLexicon_BorrowFromCString(OVLexicon *uk,ov_char8 *str)
{
  ov_word hash = _GetCStringHash((ov_uchar8*)str);
  OVreturn_word search = OVOneToOne_GetForward(uk->up, hash);
  register ov_word index = 0;
  ov_word cur_index = 0;
  ov_char8 *c;

  if(OVreturn_IS_ERROR(search)) {
    return search;
  } else {
    /* the hash key has been registered, so follow the chain
     * looking for a match... */
    register ov_char8 *data = uk->data;
    register lex_entry *entry = uk->entry,*entry_ptr;
    cur_index = (index = search.word);
    while(index) { /* found */
      c = data + (entry_ptr = (entry + index))->offset;
      if(strcmp(c,str)!=0) {/* verify match */
        index = entry_ptr->next; /* not a match? keep looking */
      } else {
        break;
      }
    }
  }

  if(!index) { 
    OVreturn_word result = { OVstatus_NOT_FOUND };
    result.word = index;
    return result;
  } else {
    OVreturn_word result = { OVstatus_SUCCESS };
    result.word = index;
    return result;
  }
}

ov_char8 *OVLexicon_FetchCString(OVLexicon *uk,ov_word id)
{
  if(id<=(ov_word)uk->n_entry)
    return uk->data + uk->entry[id].offset;
  return NULL;
}

