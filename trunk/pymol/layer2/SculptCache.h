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
#ifndef _H_SculptCache
#define _H_SculptCache

#include"Sculpt.h"

typedef struct SculptCacheEntry {
  int rest_type;
  int id0,id1,id2,id3;
  float value;
  int next;
} SculptCacheEntry;


struct _CSculptCache {
  int NCached;
  int *Hash;
  SculptCacheEntry *List;
};

int SculptCacheInit(PyMOLGlobals *G);
void SculptCacheFree(PyMOLGlobals *G);
void SculptCachePurge(PyMOLGlobals *G);

int SculptCacheQuery(PyMOLGlobals *G,int rest_type,int id0,int id1,int id2,int id3,float *value);
void SculptCacheStore(PyMOLGlobals *G,int rest_type,int id0,int id1,int id2,int id3,float value);

#endif


