
#include "OVOneToOne.h"
#include "OVHeapArray.h"
#include "ov_utility.h"

#define HASH(value,mask) (((value^(value>>24))^((value>>8)^(value>>16)))&mask)


/* FYI: "up" stands for UniquePair -- a precursor to OneToOne */

typedef struct {
  int active;
  ov_word forward_value, reverse_value;
  ov_size forward_next, reverse_next;
} up_element;

struct _OVOneToOne {
  OVHeap *heap;
  ov_uword mask;
  ov_size size, n_inactive;
  ov_word next_inactive;
  up_element *elem;
  ov_word *forward;
  ov_word *reverse;
};

OVstatus OVOneToOne_Init(OVOneToOne * up, OVHeap * heap)
{
  ov_utility_zero_range(up, up + 1);
  up->heap = heap;
  return_OVstatus_SUCCESS;
}

OVOneToOne *OVOneToOne_New(OVHeap * heap)
{
  OVOneToOne *up;
  up = OVHeap_ALLOC(heap, OVOneToOne);
  up->heap = heap;
  return up;
}

void OVOneToOne_Purge(OVOneToOne * up)
{
  if(up) {
    OVHeapArray_FREE_AUTO_NULL(up->elem);
    OVHeap_FREE_AUTO_NULL(up->heap, up->forward);
    OVHeap_FREE_AUTO_NULL(up->heap, up->reverse);
  }
}

void OVOneToOne_Del(OVOneToOne * up)
{
  if(up) {
    OVOneToOne_Purge(up);
    OVHeap_FREE_AUTO_NULL(up->heap, up);
  }
}

void OVOneToOne_Reset(OVOneToOne * up)
{
  OVOneToOne_Purge(up);
  OVOneToOne_Init(up, up->heap);
}

void OVOneToOne_Dump(OVOneToOne * up)
{
  ov_uword a;
  ov_boolean empty = OV_TRUE;
  if(up && up->mask) {
    for(a = 0; a <= up->mask; a++) {
      if(up->forward[a] || up->reverse[a]) {
        fprintf(stderr,
                " OVOneToOne_Dump: Hashes forward[0x%02x]->%d    reverse[0x%02x]->%d\n",
                (unsigned int) a, (int) up->forward[a],
                (unsigned int) a, (int) up->reverse[a]);
        empty = OV_FALSE;
      }
    }

    for(a = 0; a < up->size; a++)
      if(up->elem[a].active) {
        fprintf(stderr,
                " OVOneToOne_Dump: Elements %d:    %d (->%d)    %d (->%d)\n",
                (int) a + 1,
                (int) up->elem[a].forward_value,
                (int) up->elem[a].forward_next,
                (int) up->elem[a].reverse_value, (int) up->elem[a].reverse_next);
        empty = OV_FALSE;
      }
  }
  if(empty) {
    fprintf(stderr, " OVOneToOne_Dump: Empty. \n");
  }
}

static void Reload(OVOneToOne * up)
{                               /* assumes hash tables are clean and initialized to zero */
#ifdef DEBUG_UP
  fprintf(stderr, "Reload-Debug: entered\n");
#endif
  register ov_uword mask = up->mask;

  if(up->elem && mask) {
    {
      register up_element *elem = up->elem;
      register ov_uword a;
      for(a = 0; a < up->size; a++) {
        if(elem->active) {
          elem->forward_next = 0;       /* 0 is the sentinel for end of list */
          elem->reverse_next = 0;
        }
        elem++;
      }
    }

    {
      register ov_uword a;
      register ov_word *forward = up->forward;
      register ov_word *reverse = up->reverse;
      register up_element *elem = up->elem;
      {
        register ov_word fwd, rev;
        register ov_word fwd_val;
        register ov_word rev_val;
        for(a = 0; a < up->size; a++) {
          if(elem->active) {
            fwd_val = elem->forward_value;
            rev_val = elem->reverse_value;
            fwd = HASH(fwd_val, mask);
            rev = HASH(rev_val, mask);
            elem->forward_next = forward[fwd];
            forward[fwd] = a + 1;       /* NOTE: 1 based indices */
            elem->reverse_next = reverse[rev];
            reverse[rev] = a + 1;
          }
          elem++;
        }
      }
    }
  }
#ifdef DEBUG_UP
  {
    ov_uword a;
    for(a = 0; a <= up->mask; a++) {
      fprintf(stderr, "Reload-Debug: forward[%d]=%d, reverse[%d]=%d\n",
              a, up->forward[a], a, up->reverse[a]);
    }
  }
#endif

}

OVreturn_word OVOneToOne_GetReverse(OVOneToOne * up, ov_word reverse_value)
{
  if(!up) {
    OVreturn_word result = { OVstatus_NULL_PTR };
    return result;
  } else {
    if(up->mask) {
      ov_word hash = HASH(reverse_value, up->mask);
      register up_element *elem = up->elem;
      register ov_word index = up->reverse[hash];
      register up_element *cur_elem = elem + (index - 1);

      while(index) {
        if(cur_elem->reverse_value == reverse_value) {
          OVreturn_word result = { OVstatus_SUCCESS };
          result.word = cur_elem->forward_value;
          return result;
        }
        index = cur_elem->reverse_next;
        cur_elem = elem + (index - 1);
      }
    }
    {
      OVreturn_word result = { OVstatus_NOT_FOUND };
      return result;
    }
  }
}

OVreturn_word OVOneToOne_IterateForward(OVOneToOne * up, ov_word * hidden)
{
  if(!up) {
    OVreturn_word result = { OVstatus_NULL_PTR };
    return result;
  } else {
    OVreturn_word result = { OVstatus_YES };
    register unsigned int a;
    register up_element *cur_elem = up->elem + (*hidden);
    for(a = *hidden; a < up->size; a++) {
      if(cur_elem->active) {
        result.word = cur_elem->forward_value;
        *hidden = a + 1;
        return result;
      }
      cur_elem++;
    }
    *hidden = 0;
    result.status = OVstatus_NO;
    return result;
  }
}

OVreturn_word OVOneToOne_GetForward(OVOneToOne * up, ov_word forward_value)
{
  if(!up) {
    OVreturn_word result = { OVstatus_NULL_PTR };
    return result;
  } else {
    register ov_uword mask = up->mask;
    if(mask) {
      ov_word hash = HASH(forward_value, mask);
      register up_element *elem = up->elem;
      register ov_word index = up->forward[hash];
      register up_element *cur_elem = elem + (index - 1);
      while(index) {
        if(cur_elem->forward_value == forward_value) {
          OVreturn_word result = { OVstatus_SUCCESS };
          result.word = cur_elem->reverse_value;
          return result;
        }
        index = cur_elem->forward_next;
        cur_elem = elem + (index - 1);
      }
    }
    {
      OVreturn_word result = { OVstatus_NOT_FOUND };
      return result;
    }
  }
}

static OVstatus Recondition(OVOneToOne * up, ov_uword size, int force)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    ov_uword mask = up->mask;
#ifdef DEBUG_UP
    fprintf(stderr, "Recondition-Debug: entered for size %d.\n", size);
#endif
    if((size > mask) || ((size << 2) < mask) || force) {

      while((size << 2) < mask) {
        mask = mask >> 1;
        if(mask < 2)
          break;
      }

      while(size > mask) {
        mask = (mask << 1) + 1;
      }

#ifdef DEBUG_UP
      fprintf(stderr, "Recondition-Debug: mask %d\n", mask);
#endif
      {
        if(!up->elem) {
          up->elem = OVHeapArray_CALLOC(up->heap, up_element, size);
          if(!up->elem) {
            return_OVstatus_OUT_OF_MEMORY;
          }
        }
        if(mask != up->mask) {
          ov_word *tmp_forward = OVHeap_CALLOC(up->heap, ov_word, mask + 1);
          ov_word *tmp_reverse = OVHeap_CALLOC(up->heap, ov_word, mask + 1);
          if(!(tmp_forward && tmp_reverse)) {   /* validate */
            OVHeap_FREE_AUTO_NULL(up->heap, tmp_forward);
            OVHeap_FREE_AUTO_NULL(up->heap, tmp_reverse);
            /* being unable to condition is not an error */
          } else {
            /* impossible to fail after here... */
            OVHeap_FREE_AUTO_NULL(up->heap, up->forward);
            OVHeap_FREE_AUTO_NULL(up->heap, up->reverse);
            up->forward = tmp_forward;
            up->reverse = tmp_reverse;
            up->mask = mask;
          }
        } else {
          ov_utility_zero_range(up->forward, up->forward + (up->mask + 1));
          ov_utility_zero_range(up->reverse, up->reverse + (up->mask + 1));
        }
        Reload(up);
      }
    }
  }
  return_OVstatus_SUCCESS;
}

OVstatus OVOneToOne_Pack(OVOneToOne * up)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    if(up->n_inactive && up->elem) {
      ov_uword new_size = 0;
      up_element *src = up->elem, *dst = up->elem;
      ov_uword a;

      for(a = 0; a < up->size; a++) {
        if(src->active) {
          if(src > dst) {
            *dst = *src;
          }
          dst++;
          new_size++;
        }
        src++;
      }
      up->n_inactive = 0;
      up->next_inactive = 0;
      if(new_size < up->size) {
        if(!OVHeapArray_SET_SIZE(up->elem, up_element, new_size))
          ov_utility_zero_range(up->elem + new_size, up->elem + up->size);
      }
      up->size = new_size;
      return Recondition(up, new_size, OV_TRUE);
    }
    return_OVstatus_SUCCESS;
  }
}

OVreturn_size OVOneToOne_GetSize(OVOneToOne * up)
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

OVstatus OVOneToOne_DelReverse(OVOneToOne * up, ov_word reverse_value)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    register ov_word mask = up->mask;
    if(mask) {
      register ov_word rev_hash = HASH(reverse_value, mask);
      register ov_word rev = up->reverse[rev_hash];

      if(!rev) {
        return_OVstatus_NOT_FOUND;
      } else {
        register up_element *rev_elem = NULL;
        register up_element *elem = up->elem;
        register ov_word rev_last = 0;

        while(rev) {
          rev_elem = elem + (rev - 1);
          if(rev_elem->reverse_value == reverse_value)
            break;
          rev_last = rev;
          rev = rev_elem->reverse_next;
        }

        if(rev_elem) {
          register ov_word rev_fwd_val = rev_elem->forward_value;
          register ov_word fwd_hash = HASH(rev_fwd_val, mask);
          register ov_word fwd = up->forward[fwd_hash];
          register ov_word fwd_last = 0;
          register up_element *fwd_elem = NULL;

          while(fwd) {
            fwd_elem = elem + (fwd - 1);
            if(fwd_elem == rev_elem)
              break;
            fwd_last = fwd;
            fwd = fwd_elem->forward_next;
          }

          if(rev && (rev == fwd)) {

            /* excise element */

            if(rev_last)
              up->elem[rev_last - 1].reverse_next = rev_elem->reverse_next;
            else
              up->reverse[rev_hash] = rev_elem->reverse_next;

            if(fwd_last)
              up->elem[fwd_last - 1].forward_next = fwd_elem->forward_next;
            else
              up->forward[fwd_hash] = fwd_elem->forward_next;

            /* store as inactive */

            rev_elem->active = OV_FALSE;
            rev_elem->forward_next = up->next_inactive;
            up->next_inactive = rev;
            up->n_inactive++;
            if(up->n_inactive > (up->size >> 1))        /* over half of bits are inactive */
              OVOneToOne_Pack(up);
            return_OVstatus_SUCCESS;
          }
        }
      }
    }
    return_OVstatus_NOT_FOUND;
  }
}

OVstatus OVOneToOne_DelForward(OVOneToOne * up, ov_word forward_value)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    register ov_word mask = up->mask;
    if(mask) {
      register ov_word fwd_hash = HASH(forward_value, mask);
      register ov_word fwd = up->forward[fwd_hash];
      if(!fwd) {
        return_OVstatus_NOT_FOUND;
      } else {
        register up_element *fwd_elem = NULL;
        register up_element *elem = up->elem;
        register ov_word fwd_last = 0;

        while(fwd) {
          fwd_elem = elem + (fwd - 1);
          if(fwd_elem->forward_value == forward_value)
            break;
          fwd_last = fwd;
          fwd = fwd_elem->forward_next;
        }

        if(fwd_elem) {
          register ov_word fwd_rev_val = fwd_elem->reverse_value;
          register ov_word rev_hash = HASH(fwd_rev_val, mask);
          register ov_word rev = up->reverse[rev_hash];
          register ov_word rev_last = 0;
          up_element *rev_elem = NULL;

          while(rev) {
            rev_elem = elem + (rev - 1);
            if(rev_elem == fwd_elem)
              break;
            rev_last = rev;
            rev = rev_elem->reverse_next;
          }

          if(fwd && (fwd == rev)) {

            /* excise elements */

            if(fwd_last)
              up->elem[fwd_last - 1].forward_next = fwd_elem->forward_next;
            else
              up->forward[fwd_hash] = fwd_elem->forward_next;

            if(rev_last)
              up->elem[rev_last - 1].reverse_next = rev_elem->reverse_next;
            else
              up->reverse[rev_hash] = rev_elem->reverse_next;

            /* store as inactive */

            fwd_elem->active = OV_FALSE;
            fwd_elem->forward_next = up->next_inactive;
            up->next_inactive = fwd;
            up->n_inactive++;
            if(up->n_inactive > (up->size >> 1))        /* over half of bits are inactive */
              OVOneToOne_Pack(up);
            return_OVstatus_SUCCESS;
          }
        }
      }
    }
    return_OVstatus_NOT_FOUND;
  }
}

void OVOneToOne_Stats(OVOneToOne * up)
{
  if(up && up->mask) {
    int max_len = 0;
    ov_uword a;
    for(a = 0; a < up->mask; a++) {
      {
        ov_word index = up->forward[a];
        up_element *elem = up->elem;
        int cnt = 0;
        if(index) {
          up_element *cur_elem;
          while(index) {
            cur_elem = elem + (index - 1);
            index = cur_elem->forward_next;
            cnt++;
          }
          if(cnt > max_len)
            max_len = cnt;
        }
      }

      {
        ov_word index = up->reverse[a];
        up_element *elem = up->elem;
        int cnt = 0;
        if(index) {
          up_element *cur_elem;
          while(index) {
            cur_elem = elem + (index - 1);
            index = cur_elem->reverse_next;
            cnt++;
          }
          if(cnt > max_len)
            max_len = cnt;
        }
      }

    }
    fprintf(stderr, " OVOneToOne_Stats: MaxLen=%d ", max_len);
    fprintf(stderr, "active=%d n_inactive=%d ", (int)( up->size - up->n_inactive),
            (int) up->n_inactive);
    fprintf(stderr, "mask=0x%x n_alloc=%lu\n", (unsigned int) up->mask,
            (unsigned long) OVHeapArray_GET_SIZE(up->elem));
  }
}

OVstatus OVOneToOne_Set(OVOneToOne * up, ov_word forward_value, ov_word reverse_value)
{
  if(!up) {
    return_OVstatus_NULL_PTR;
  } else {
    register ov_word mask = up->mask;
    register ov_word fwd_hash = HASH(forward_value, mask);
    register ov_word rev_hash = HASH(reverse_value, mask);
    register up_element *fwd_elem = NULL;
    register up_element *rev_elem = NULL;
    register ov_word fwd;
    register ov_word rev;
    if(!mask) {
      fwd = 0;
      rev = 0;
    } else {

      fwd = up->forward[fwd_hash];
      rev = up->reverse[rev_hash];

#ifdef DEBUG_OVOneToOne
      fprintf(stderr, "OVOneToOneSet-Debug: set %d,%d\n", forward_value, reverse_value);
      fprintf(stderr, "OVOneToOneSet-Debug: fwd_hash %d rev_hash %d mask %d size %d\n",
              fwd_hash, rev_hash, up->mask, up->size);
      fprintf(stderr, "OVOneToOneSet-Debug: before search fwd rev %d %d\n", fwd, rev);
#endif

      {                         /* find elements if they exist, and detect erroneous conditions */

        register up_element *elem = up->elem;

#ifdef DEBUG_OVOneToOne
        {
          int a;
          for(a = 0; a < up->size; a++) {
            fprintf(stderr, "OVOneToOneSet-Debug: on entry %d forward_next: %d\n",
                    a + 1, elem[a].forward_next);
          }

          for(a = 0; a <= up->mask; a++) {
            fprintf(stderr,
                    "OVOneToOneSet-Debug: on entry %d forward hash %d: reverse_hash: %d\n",
                    a, up->forward[a], up->reverse[a]);
          }
        }
#endif

        while(fwd) {
          fwd_elem = elem + (fwd - 1);
          if(fwd_elem->forward_value == forward_value)
            break;
          fwd = fwd_elem->forward_next;
        }
        while(rev) {
          rev_elem = elem + (rev - 1);
          if(rev_elem->reverse_value == reverse_value)
            break;
          rev = rev_elem->reverse_next;
        }
      }

      if((fwd && (!rev)) || (rev && (!fwd))) {
        return_OVstatus_DUPLICATE;
      }
    }
    if(!(fwd || rev)) {
      ov_size new_index;
      /* new pair */
#ifdef DEBUG_OVOneToOne
      fprintf(stderr, "OVOneToOneSet-Debug: New pair.\n");
#endif
      if(up->n_inactive) {
        new_index = up->next_inactive;
        up->next_inactive = up->elem[new_index - 1].forward_next;
        up->n_inactive--;
      } else {
        if(up->elem && (!OVHeapArray_CHECK(up->elem, up_element, up->size))) {
          return_OVstatus_OUT_OF_MEMORY;
        } else {
          OVstatus result;
          if(OVreturn_IS_ERROR(result = Recondition(up, up->size + 1, OV_FALSE))) {
            return result;
          } else {
            /* guaranteed to succeed past this point, so we can increase size */
            new_index = ++up->size;
          }
        }
      }
      {
        up_element *elem = up->elem + (new_index - 1);
        elem->forward_value = forward_value;
        elem->reverse_value = reverse_value;
        elem->active = OV_TRUE;

        /* regenerate new hashes */
        mask = up->mask;
        fwd_hash = HASH(forward_value, mask);
        rev_hash = HASH(reverse_value, mask);

        {
          ov_word *forward_start_index = up->forward + fwd_hash;
          ov_word *reverse_start_index = up->reverse + rev_hash;

          elem->forward_next = *forward_start_index;
          *forward_start_index = new_index;     /* note the +1 offset */
          elem->reverse_next = *reverse_start_index;
          *reverse_start_index = new_index;     /* note the +1 offset */
        }

#ifdef DEBUG_OVOneToOne
        {
          int a;
          for(a = 0; a <= up->mask; a++) {
            fprintf(stderr, "OVOneToOneSet-Debug: forward[%d]=%d, reverse[%d]=%d\n",
                    a, up->forward[a], a, up->reverse[a]);
          }
        }
#endif

      }

    } else if(fwd_elem != rev_elem) {
      return_OVstatus_MISMATCH;
    } else {
      return_OVstatus_NO_EFFECT;
      /* exists and matched, so do nothing */
    }
  }
  return_OVstatus_SUCCESS;
}
