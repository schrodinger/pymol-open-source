/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
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

#include"MemoryCache.h"
#include"MemoryDebug.h"
#include"Setting.h"

typedef struct {
  int *ptr;
  unsigned int size;
} MemoryCacheRec;

typedef MemoryCacheRec MemoryCacheThread[cMemoryCache_max_block];
static MemoryCacheThread MemoryCache[cMemoryCache_max_group];

void MemoryCacheInit(void)
{
  MemoryZero((void*)MemoryCache, 
             (void*)(((&MemoryCache[0][0] + 
                       (cMemoryCache_max_block*cMemoryCache_max_group)))));

}

void *MemoryCacheMalloc(unsigned int size,int group_id,int block_id)
{
  MemoryCacheRec *rec = &MemoryCache[group_id][block_id];
  
  if((group_id<0)||(!SettingGet(cSetting_cache_memory)))
    return(mmalloc(size));
  if(!rec->ptr) {
    rec->size = size;
    rec->ptr = mmalloc(size);
  } else if(rec->size<size) {
    rec->size = size;
    mfree(rec->ptr);
    rec->ptr = mmalloc(size);
  }
  return(rec->ptr);
}

void *MemoryCacheCalloc(unsigned int number, unsigned int size,int group_id,int block_id)
{
  MemoryCacheRec *rec = &MemoryCache[group_id][block_id];
  unsigned int true_size = number * size;

  if((group_id<0)||(!(int)SettingGet(cSetting_cache_memory)))
    return(mcalloc(number,size));
  if(!rec->ptr) {
    rec->size = true_size;
    rec->ptr = mcalloc(number,size);
  } else if(rec->size<true_size) {
    mfree(rec->ptr);
    rec->size = true_size;
    rec->ptr = mcalloc(number,size);
  } else {
    mfree(rec->ptr);
    rec->size = true_size;
    rec->ptr = mcalloc(number,size);
  }
  return(rec->ptr);
}

void *MemoryCacheRealloc(void *ptr, unsigned int size,int group_id, int block_id)
{
  /* no checking done */
  MemoryCacheRec *rec = &MemoryCache[group_id][block_id];

  if((group_id<0)||(!(int)SettingGet(cSetting_cache_memory)))
    return(mrealloc(ptr,size));
  if(ptr!=rec->ptr)
    printf("Error: Memory Cache Mismatch 2 %d %d\n",group_id,block_id);
  if(!rec->ptr) {
    rec->size= size;
    rec->ptr = mrealloc(ptr,size);
  } else if(rec->size<size) {
    rec->size = size;
    rec->ptr = mrealloc(ptr,size);
  }
  return(rec->ptr);
}

void MemoryCacheFree(void *ptr,int group_id, int block_id,int force)
{
  MemoryCacheRec *rec = &MemoryCache[group_id][block_id];
  if((group_id<0)||(!(int)SettingGet(cSetting_cache_memory)))
    return(mfree(ptr));
  if(ptr!=rec->ptr)
    printf("Error: Memory Cache Mismatch 2 %d %d\n",group_id,block_id);
  if(force) {
    mfree(rec->ptr);
    rec->ptr = NULL;
  }
}

void MemoryCacheDone(void)
{
  int a,b;

  for(a=0;a<cMemoryCache_max_group;a++) {
    for(b=0;b<cMemoryCache_max_block;b++) {
      MemoryCacheRec *rec = &MemoryCache[a][b];
      if(rec->ptr)
        mfree(rec->ptr);
    }
  }
}
