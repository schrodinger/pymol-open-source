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

#ifdef _MemoryCache_ON

#include"MemoryDebug.h"
#include"Setting.h"


typedef struct {
  int *ptr;
  unsigned int size;
} MemoryCacheRec;

typedef MemoryCacheRec MemoryCacheThread[cMemoryCache_max_block];

struct _CMemoryCache {
  MemoryCacheThread Cache[cMemoryCache_max_group];
};

void MemoryCacheInit(PyMOLGlobals *G)
{
  G->MemoryCache=Calloc(CMemoryCache,1);
}

void *MemoryCacheMalloc(PyMOLGlobals *G,unsigned int size,int group_id,int block_id)
{
  register CMemoryCache *I = G->MemoryCache;
  register MemoryCacheRec *rec = &I->Cache[group_id][block_id];
  
  if((group_id<0)||(!SettingGet(G,cSetting_cache_memory)))
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

void *MemoryCacheCalloc(PyMOLGlobals *G,unsigned int number, unsigned int size,int group_id,int block_id)
{
  register CMemoryCache *I = G->MemoryCache;
  register MemoryCacheRec *rec = &I->Cache[group_id][block_id];
  unsigned int true_size = number * size;

  /* interesting result: calloc is faster than cacheing */

  if((group_id<0)||(!(int)SettingGet(G,cSetting_cache_memory)))
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

void *MemoryCacheRealloc(PyMOLGlobals *G,void *ptr, unsigned int size,int group_id, int block_id)
{
  /* no checking done */
  register CMemoryCache *I = G->MemoryCache;
  register MemoryCacheRec *rec = &I->Cache[group_id][block_id];

  if((group_id<0)||(!(int)SettingGet(G,cSetting_cache_memory)))
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

void MemoryCacheFree(PyMOLGlobals *G,void *ptr,int group_id, int block_id,int force)
{
  register CMemoryCache *I = G->MemoryCache;
  MemoryCacheRec *rec = &I->Cache[group_id][block_id];
  if((group_id<0)||(!(int)SettingGet(G,cSetting_cache_memory))) {
    mfree(ptr);
    return;
  }
  if(rec->ptr&&(ptr!=rec->ptr))
    printf("Error: Memory Cache Mismatch 2 %d %d\n",group_id,block_id);
  if(force) {
    if(rec->ptr) 
      mfree(rec->ptr);
    rec->ptr = NULL;
  }
}

void MemoryCacheDone(PyMOLGlobals *G)
{
  int a,b;
  register CMemoryCache *I = G->MemoryCache;
  for(a=0;a<cMemoryCache_max_group;a++) {
    for(b=0;b<cMemoryCache_max_block;b++) {
      MemoryCacheRec *rec = &I->Cache[a][b];
      if(rec->ptr)
        mfree(rec->ptr);
    }
  }
  FreeP(G->MemoryCache);
}
#endif
