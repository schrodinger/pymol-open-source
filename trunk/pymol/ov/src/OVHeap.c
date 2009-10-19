#include"OVHeap.h"

#ifndef OV_JX
#ifdef OVHeap_TRACKING


/* all of the code in this module is only used when OVHeap_TRACKING is
 * active otherwise, OVHeap is simply a Macro-wrapper around stdlib's
 * memory management routines */

#define OVHeap_SORT

#define OVHeap_ERROR_LOG stderr

static void HandleOutOfMemory(void)
{
  fprintf(OVHeap_ERROR_LOG,
          "OVHeap-Error: *************************************************\n");
  fprintf(OVHeap_ERROR_LOG,
          "OVHeap-Error: *** EEK! This program just ran out of memory! ***\n");
  fprintf(OVHeap_ERROR_LOG,
          "OVHeap-Error: *************************************************\n");
  fflush(OVHeap_ERROR_LOG);
#ifdef OVHeap_ABORT_ON_ERROR
  OVHeap_Dump(false);
  abort();
#endif

}

static void HandleError(void)
{
  fprintf(OVHeap_ERROR_LOG,
          "OVHeap-Error: *************************************************\n");
  fprintf(OVHeap_ERROR_LOG,
          "OVHeap-Error: *** EEK! Memory corruption may have occurred! ***\n");
  fprintf(OVHeap_ERROR_LOG,
          "OVHeap-Error: *************************************************\n");
  fflush(OVHeap_ERROR_LOG);

#ifdef OVHeap_ABORT_ON_ERROR
  OVHeap_Dump(false);
  abort();
#endif
}

typedef struct DebugRec {
  struct DebugRec *next;
  char file[32], note[64];
  int line;
  ov_size size;
  int type;
} DebugRec;

#define HASH(x) ((x>>11)&0x3FF)

struct _OVHeap {
  DebugRec *HashTable[1024];
  int Count;
  int MaxCount;
};

OVHeap *_OVHeap_New(void)
{
  OVHeap *I = ov_os_calloc(1, sizeof(OVHeap));
  return I;
}

void _OVHeap_Del(OVHeap * I)
{
#ifndef OV_JX
#ifdef OVHeap_TRACKING
  OVHeap_Dump(I, 0);
#endif
#endif
  ov_os_free((void *) I);
}

int OVHeap_Usage(OVHeap * I)
{
  int a;
  unsigned int tot = 0;
  DebugRec *rec;
  for(a = 0; a < 1024; a++) {
    rec = I->HashTable[a];
    while(rec) {
      tot += rec->size;
      rec = rec->next;
    }
  }
  return (tot);
}

typedef char OutputLine[1280];

typedef struct {
  void *p;
  char text[1280];
} DebugOutput;

static int OutputInOrder(DebugOutput * output, int line1, int line2)
{
  return (output[line1].p <= output[line2].p);
}

static void OVHeap_SortOutput(int n, DebugOutput * array, int *x)
{
  int l, a, r, t, i;

  if(n < 1)
    return;
  else if(n == 1) {
    x[0] = 0;
    return;
  }
  x--;
  for(a = 1; a <= n; a++)
    x[a] = a;
  l = (n >> 1) + 1;
  r = n;
  while(1) {
    if(l > 1)
      t = x[--l];
    else {
      t = x[r];
      x[r] = x[1];
      if(--r == 1) {
        x[1] = t;
        break;
      }
    }
    i = l;
    a = l << 1;
    while(a <= r) {
      if(a < r && (!OutputInOrder(array, x[a + 1] - 1, x[a] - 1)))
        a++;
      if(!OutputInOrder(array, x[a] - 1, t - 1)) {
        x[i] = x[a];
        a += (i = a);
      } else
        a = r + 1;
    }
    x[i] = t;
  }
  x++;
  for(a = 0; a < n; a++)
    x[a]--;
}

void OVHeap_Dump(OVHeap * I, ov_uint32 flags)
{
  int a;
  int cnt = 0;
  unsigned int tot = 0;
  DebugRec *rec;
  char type[] = "FV";
#ifdef OVHeap_SORT
  DebugOutput *output;
  int *index;
#endif

#ifdef OVHeap_SORT
  for(a = 0; a < 1024; a++) {
    rec = I->HashTable[a];
    while(rec) {
      rec = rec->next;
      cnt++;
    }
  }
  output = (DebugOutput *) malloc(cnt * sizeof(DebugOutput));
  index = (int *) malloc(cnt * sizeof(int));
  cnt = 0;
#endif
  fprintf(OVHeap_ERROR_LOG,
          "================================= HEAP =================================\n");
  fflush(OVHeap_ERROR_LOG);
  for(a = 0; a < 1024; a++) {
    rec = I->HashTable[a];
    while(rec) {
      {
        tot += rec->size;

        if(flags & OVHeap_DUMP_NO_ADDRESSES) {
#ifdef OVHeap_SORT
          output[cnt].p = rec;
          sprintf(output[cnt].text,
#else
          fprintf(OVHeap_ERROR_LOG,
#endif
                  "OVHeap: (%7x) %c %s", (unsigned int)
                  rec->size, type[rec->type], rec->file);

        } else {

#ifdef OVHeap_SORT
          output[cnt].p = rec;
          sprintf(output[cnt].text,
#else
          fprintf(OVHeap_ERROR_LOG,
#endif
                  "OVHeap:%8x -%8x (%7x) %c %s:%-4d",
                  (unsigned int) rec + 1,
                  (unsigned int) ((char *) (rec + 1) + rec->size), (unsigned int)
                  rec->size, type[rec->type], rec->file, rec->line);
        }

        if(flags & OVHeap_DUMP_FILES_TOO) {
          FILE *f;
          int line = rec->line;
          char buffer[1024], *c;
          f = fopen(rec->file, "r");
          while(line--)
            fgets(buffer, 2048, f);
          c = buffer;
          while(*c && *c < 33)
            c++;
#ifdef OVHeap_SORT
          strcat(output[cnt].text, c);
#else
          fprintf(OVHeap_ERROR_LOG, "%s", c);
#endif
          fclose(f);
        } else {
#ifdef OVHeap_SORT
          strcat(output[cnt].text, "\n");
#else
          fprintf(OVHeap_ERROR_LOG, "\n");
#endif
        }
      }
      rec = rec->next;
      cnt++;
    }
  }
#ifdef OVHeap_SORT
  OVHeap_SortOutput(cnt, output, index);
  for(a = 0; a < cnt; a++) {
    fprintf(OVHeap_ERROR_LOG, "%s", output[index[a]].text);
  }
#endif

  fprintf(OVHeap_ERROR_LOG,
          "OVHeap: Summary: Blocks expected %d, found %d, peaked at %d.\n", I->Count, cnt,
          I->MaxCount);
  fprintf(OVHeap_ERROR_LOG, "OVHeap: Summary: Total bytes allocated 0x%x (%0.3f MB).\n",
          tot, tot / (1024.0 * 1024));
  fprintf(OVHeap_ERROR_LOG,
          "========================================================================\n");
  fflush(OVHeap_ERROR_LOG);
}

static void OVHeap_HashAdd(OVHeap * I, DebugRec * rec)
{
  int hash;

  hash = (int) rec;
  hash = HASH(hash);
  rec->next = I->HashTable[hash];
  I->HashTable[hash] = rec;
}

static DebugRec *OVHeap_HashRemove(OVHeap * I, void *ptr)
{
  DebugRec *rec, *cur, *last;
  int hash;

  rec = (DebugRec *) ptr;
  rec--;
  hash = (int) rec;
  hash = HASH(hash);
  last = NULL;
  cur = I->HashTable[hash];
  while(cur) {
    if(cur == rec) {
      if(last)
        last->next = cur->next;
      else
        I->HashTable[hash] = cur->next;
      break;
    }
    last = cur;
    cur = cur->next;
  }
  return (cur);
}


/* 
   Ustatus OVHeap_Wrap(Ustatus input,void **ptrptr,char *file, int line, int type)
   {
   DebugRec *rec,*cur;
   int hash;
  
   rec=(DebugRec*)(*ptrptr);
   rec--;
   hash=(int)rec;
   hash=HASH(hash);
   cur=I->HashTable[hash];
   while(cur)
   {
   if(cur==rec)
   {
   strcpy(rec->file,file);
   rec->line=line;
   rec->type=type;
   break;
   }
   cur=cur->next;
   }  
   return input;
   }
*/

void *_OVHeap_Malloc(OVHeap * I, ov_size size, const char *file, int line, int type)
{
  DebugRec *rec;

  rec = (DebugRec *) malloc(sizeof(DebugRec) + size);
  if(!rec) {
    if(size) {
      fprintf(OVHeap_ERROR_LOG, "OVHeap_Malloc-Error: malloc failed %lu.\n",
              (unsigned long) size);
      HandleOutOfMemory();
    }
  } else {
    strcpy(rec->file, file);
    rec->line = line;
    rec->size = size;
    rec->type = type;
    OVHeap_HashAdd(I, rec);
    rec++;
    I->Count++;
    if(I->MaxCount < I->Count)
      I->MaxCount = I->Count;
  }
  return ((void *) rec);
}

void *_OVHeap_Calloc(OVHeap * I, ov_size num, ov_size size,
                     const char *file, int line, int type)
{
  DebugRec *rec;

  rec = (DebugRec *) calloc(1, sizeof(DebugRec) + size * num);
  if(!rec) {
    if(size) {
      fprintf(OVHeap_ERROR_LOG, "OVHeap_Malloc-Error: calloc failed %lu x %lu.\n",
              (unsigned long) size, (unsigned long) num);
      HandleOutOfMemory();
    }
  } else {
    strcpy(rec->file, file);
    rec->line = line;
    rec->size = size * num;
    rec->type = type;
    OVHeap_HashAdd(I, rec);
    rec++;
    I->Count++;
    if(I->MaxCount < I->Count)
      I->MaxCount = I->Count;
  }
  return ((void *) rec);
}

void *_OVHeap_Realloc(OVHeap * I, void *ptr, ov_size size,
                      const char *file, int line, int type)
{
  DebugRec *rec, *new_rec;

  if(!ptr)
    return (_OVHeap_Malloc(I, size, file, line, type));
  else {
    rec = OVHeap_HashRemove(I, ptr);
    if(!rec) {
      fprintf(OVHeap_ERROR_LOG,
              "OVHeap_Realloc-Error: realloc() corrupted heap or bad ptr! (%s:%i @%p)\n",
              file, line, ptr);
#ifdef OVHeap_ABORT_ON_ERROR
      OVHeap_Dump(false);
      abort();
#endif
      HandleError();
    } else if(rec->type != type) {
      fprintf(OVHeap_ERROR_LOG,
              "OVHeap_Realloc-Error: ptr %p of wrong type: %i!=%i (%s:%i)\n", ptr,
              rec->type, type, file, line);
#ifdef OVHeap_ABORT_ON_ERROR
      OVHeap_Dump(false);
      abort();
#endif
      HandleError();
    } else {
      new_rec = (DebugRec *) realloc(rec, size + sizeof(DebugRec));
      if(!new_rec) {
        fprintf(OVHeap_ERROR_LOG, "OVHeap_Realloc-Error: realloc() failed! (%s:%i)\n",
                file, line);
        HandleOutOfMemory();
        OVHeap_HashAdd(I, rec); /* put it back */
        rec++;
      } else {
        rec = new_rec;
        OVHeap_HashAdd(I, rec);
        rec->size = size;
        rec++;
        return ((void *) rec);
      }
    }
  }
  return (NULL);
}

void _OVHeap_Free(OVHeap * I, void *ptr, const char *file, int line, int type)
{
  DebugRec *rec;

  if(!ptr) {
    fprintf(OVHeap_ERROR_LOG,
            "OVHeap_Free-Error: free() called with NULL pointer (%s:%i)\n", file, line);
#ifdef OVHeap_ABORT_ON_ERROR
    OVHeap_Dump(false);
    abort();
#endif
    HandleError();
  }
  rec = OVHeap_HashRemove(I, ptr);
  if(rec) {
    if(rec->type != type) {
      fprintf(OVHeap_ERROR_LOG,
              "OVHeap_Free-Error: ptr %p type mismatch: %i!=%i (%s:%i)\n", ptr, rec->type,
              type, file, line);
#ifdef OVHeap_ABORT_ON_ERROR
      OVHeap_Dump(false);
      abort();
#endif
      HandleError();
    }
    free(rec);
  } else {

    fprintf(OVHeap_ERROR_LOG,
            "OVHeap_Free-Error: free(): corrupted tree or bad ptr! (%s:%i @%p)\n",
            file, line, ptr);
    HandleError();
  }
  I->Count--;
}

#endif
#endif
