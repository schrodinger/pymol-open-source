
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
#include"os_python.h"
#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"
#include"OOMac.h"
#include"Feedback.h"
#include"Util.h"
#include"SculptCache.h"

#define CACHE_HASH_SIZE 65536

#define cache_hash(d,e,f,g) \
((((d)      )&0x003F)|\
 (((e+g)<< 6)&0x0FC0)|\
 (((f-g)<<12)&0xF000))

static void SculptCacheCheck(PyMOLGlobals * G)
{
  register CSculptCache *I = G->SculptCache;
  if(!I->Hash) {
    I->Hash = Calloc(int, CACHE_HASH_SIZE);
  }
}

int SculptCacheInit(PyMOLGlobals * G)
{
  register CSculptCache *I = NULL;
  if((I = (G->SculptCache = Calloc(CSculptCache, 1)))) {
    I->Hash = NULL;             /* don't allocate until we need it */
    I->List = VLAlloc(SculptCacheEntry, 16);
    I->NCached = 1;
    return 1;
  } else {
    return 0;
  }
}

void SculptCachePurge(PyMOLGlobals * G)
{
  register CSculptCache *I = G->SculptCache;
  if(I->Hash) {
    FreeP(I->Hash);
    I->Hash = NULL;
  }
  I->NCached = 1;
}

void SculptCacheFree(PyMOLGlobals * G)
{
  register CSculptCache *I = G->SculptCache;
  FreeP(I->Hash);
  VLAFreeP(I->List);
  FreeP(G->SculptCache);
}

int SculptCacheQuery(PyMOLGlobals * G, int rest_type, int id0, int id1, int id2, int id3,
                     float *value)
{
  register CSculptCache *I = G->SculptCache;
  int *v, i;
  SculptCacheEntry *e;
  int found = false;
  if(!I->Hash)
    SculptCacheCheck(G);
  if(I->Hash) {
    v = I->Hash + cache_hash(id0, id1, id2, id3);
    i = (*v);
    while(i) {
      e = I->List + i;
      if((e->rest_type == rest_type) && (e->id0 == id0) && (e->id1 == id1)
         && (e->id2 == id2) && (e->id3 == id3)) {
        found = true;
        *value = e->value;
        break;
      }
      i = e->next;
    }
  }
  return (found);
}

void SculptCacheStore(PyMOLGlobals * G, int rest_type, int id0, int id1, int id2, int id3,
                      float value)
{
  register CSculptCache *I = G->SculptCache;
  int *v, i;
  SculptCacheEntry *e;
  int found = false;
  if(!I->Hash)
    SculptCacheCheck(G);
  if(I->Hash) {
    v = I->Hash + cache_hash(id0, id1, id2, id3);
    i = (*v);
    while(i) {
      e = I->List + i;;
      if((e->rest_type == rest_type) && (e->id0 == id0) && (e->id1 == id1)
         && (e->id2 == id2) && (e->id3 == id3)) {
        e->value = value;
        found = true;
        break;
      }
      i = e->next;
    }
    if(!found) {
      VLACheck(I->List, SculptCacheEntry, I->NCached);
      v = I->Hash + cache_hash(id0, id1, id2, id3);
      e = I->List + I->NCached;
      e->next = *v;
      *v = I->NCached;          /* become new head of list */

      e->rest_type = rest_type;
      e->id0 = id0;
      e->id1 = id1;
      e->id2 = id2;
      e->id3 = id3;
      e->value = value;

      I->NCached++;
    }
  }
}
