

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#ifndef _H_MemoryCache
#define _H_MemoryCache

#define _MemoryCache_OFF


/* NOTE: at the present time, MemoryCaching is not compatible with a
   running PyMOL free of global state -- it's basically just a
   performace hack for systems with slow memory management.
*/


/* cacheable memory blocks (really just for the ray-tracer)  */

#define cCache_no_cache                                   0
#define cCache_ray_antialias_buffer                       1
#define cCache_ray_vert2prim                              2
#define cCache_ray_primitive                              3

#define cCache_basis_vertex                               4

#define cCache_basis_radius                               5

#define cCache_basis_radius2                              6
#define cCache_basis_vert2normal                          7
#define cCache_basis_normal                               8

#define cCache_basis_precomp                              9

#define cCache_ray_basis                                 10

#define cCache_ray_map                                   11

#define cCache_map_head_offset                           1
#define cCache_map_link_offset                           2
#define cCache_map_ehead_offset                          3
#define cCache_map_emask_offset                          4
#define cCache_map_elist_offset                          5
#define cCache_map_ehead_new_offset                      6
#define cCache_map_elist_new_offset                      7

#define cCache_basis_tempVertex                          31
#define cCache_basis_tempRef                             32

#define cCache_basis_site                                33
#define cCache_basis_value                               34

#define cCache_map_scene_cache                           40
#define cCache_map_shadow_cache                          44
#define cCache_map_cache_offset                          1
#define cCache_map_cache_link_offset                     2

#define cCache_ray_edging_buffer                         50


/* other constants */

#define cMemoryCache_max_block 100
#define cMemoryCache_max_group 16

#ifdef _MemoryCache_ON

#include "PyMOLGlobals.h"
#include "MemoryDebug.h"

void MemoryCacheInit(PyMOLGlobals * G);
void MemoryCacheDone(PyMOLGlobals * G);
void MemoryCacheReplaceBlock(PyMOLGlobals * G, int group_id, int old_block_id,
                             int new_block_id);

void *_MemoryCacheMalloc(PyMOLGlobals * G, unsigned int size, int group_id,
                         int block_id MD_FILE_LINE_Decl);
void *_MemoryCacheCalloc(PyMOLGlobals * G, unsigned int number, unsigned int size,
                         int group_id, int block_id MD_FILE_LINE_Decl);
void *_MemoryCacheRealloc(PyMOLGlobals * G, void *ptr, unsigned int size, int group_id,
                          int block_id MD_FILE_LINE_Decl);
void *_MemoryCacheShrinkForSure(PyMOLGlobals * G, void *ptr, unsigned int size,
                                int group_id, int block_id MD_FILE_LINE_Decl);
void _MemoryCacheFree(PyMOLGlobals * G, void *ptr, int group_id, int block_id,
                      int force MD_FILE_LINE_Decl);

#define CacheAlloc(G,type,size,thread,id) (type*)_MemoryCacheMalloc(G,sizeof(type)*(size),thread,id MD_FILE_LINE_Call)
#define CacheCalloc(G,type,size,thread,id) (type*)_MemoryCacheCalloc(G,sizeof(type),size,thread,id MD_FILE_LINE_Call)
#define CacheRealloc(G,ptr,type,size,thread,id) (type*)_MemoryCacheRealloc(G,ptr,sizeof(type)*(size),thread,id MD_FILE_LINE_Call)
#define CacheFreeP(G,ptr,thread,id,force) {if(ptr) {_MemoryCacheFree(G,ptr,thread,id,force MD_FILE_LINE_Call);ptr=NULL;}}

#define VLACacheCheck(G,ptr,type,rec,t,i) (ptr=(type*)(((((unsigned)rec)>=((VLARec*)(ptr))[-1].nAlloc) ? VLACacheExpand(G,ptr,(rec),t,i) : (ptr))))
#define VLACacheAlloc(G,type,initSize,t,i) (type*)VLACacheMalloc(G,initSize,sizeof(type),3,0,t,i)
#define VLACacheFreeP(G,ptr,t,i,f) {if(ptr) {VLACacheFree(G,ptr,t,i,f);ptr=NULL;}}
#define VLACacheSize(G,ptr,type,size,t,i) {ptr=(type*)VLACacheSetSize(G,ptr,size,t,i);}
#define VLACacheSizeForSure(G,ptr,type,size,t,i) {ptr=(type*)VLACacheSetSizeForSure(G,ptr,size,t,i);}

#ifndef _MemoryDebug_ON
void *VLACacheMalloc(PyMOLGlobals * G, unsigned int initSize, unsigned int recSize, unsigned int growFactor, int autoZero, int thread, int index);      /*growfactor 1-10 */
#else
#define VLACacheMalloc(G,a,b,c,d,t,i) _VLACacheMalloc(G,__FILE__,__LINE__,a,b,c,d,t,i)
void *_VLACacheMalloc(PyMOLGlobals * G, const char *file, int line, unsigned int initSize, unsigned int recSize, unsigned int growFactor, int autoZero, int thread, int index); /*growfactor 1-10 */
#endif

void VLACacheFree(PyMOLGlobals * G, void *ptr, int thread, int id, int force);
void *VLACacheSetSize(PyMOLGlobals * G, void *ptr, unsigned int newSize, int group_id,
                      int block_id);
void *VLACacheSetSizeForSure(PyMOLGlobals * G, void *ptr, unsigned int newSize,
                             int group_id, int block_id);
void *VLACacheExpand(PyMOLGlobals * G, void *ptr, unsigned int rec, int thread_index,
                     int block_id);

#else


/* memory cache off */

#define VLACacheCheck(G,ptr,type,rec,t,i) VLACheck(ptr,type,rec)
#define VLACacheAlloc(G,type,initSize,t,i) VLAlloc(type,initSize)
#define VLACacheFreeP(G,ptr,t,i,f) VLAFreeP(ptr)
#define VLACacheSize(G,ptr,type,size,t,i) VLASize(ptr,type,size)
#define VLACacheSizeForSure(G,ptr,type,size,t,i) VLASizeForSure(ptr,type,size)
#define VLACacheExpand(G,ptr,rec,thread_index,i) VLAExpand(ptr,rec)

#define MemoryCacheInit(x)
#define MemoryCacheDone(x)
#define MemoryCacheReplaceBlock(G,g,o,n)

#define VLACacheMalloc(G,a,b,c,d,t,i) VLAMalloc(a,b,c,d)
#define VLACacheFree(G,p,t,i,f) VLAFree(p)

#define CacheAlloc(G,type,size,thread,id) (type*)mmalloc(sizeof(type)*(size))
#define CacheCalloc(G,type,size,thread,id) (type*)mcalloc(sizeof(type),size)
#define CacheRealloc(G,ptr,type,size,thread,id) (type*)mrealloc(sizeof(type)*(size))
#define CacheFreeP(G,ptr,thread,id,force) {if(ptr) {mfree(ptr);ptr=NULL;}}

#endif

#endif
