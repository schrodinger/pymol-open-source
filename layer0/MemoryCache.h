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

void MemoryCacheInit(void);
void *MemoryCacheMalloc(unsigned int size,int group_id,int block_id);
void *MemoryCacheCalloc(unsigned int number, unsigned int size,int group_id,int block_id);
void *MemoryCacheRealloc(void *ptr, unsigned int size,int group_id, int block_id);
void MemoryCacheFree(void *ptr,int group_id, int block_id,int force);
void MemoryCacheDone(void);

#else
/* memory cache off */

#define MemoryCacheInit(x)
#define MemoryCacheMalloc(size,group_id,block_id) mmalloc(size)
#define MemoryCacheCalloc(number, size,group_id,block_id) mcalloc(number,size)
#define MemoryCacheRealloc(ptr,size,group_id,block_id) mrealloc(ptr,size)
#define MemoryCacheFree(ptr,group_id, block_id,force) mfree(ptr)
#define MemoryCacheDone()

#endif

#define CacheAlloc(type,size,thread,id) (type*)MemoryCacheMalloc(sizeof(type)*(size),thread,id)
#define CacheCalloc(type,size,thread,id) (type*)MemoryCacheCalloc(sizeof(type),size,thread,id)
#define CacheRealloc(ptr,type,size,thread,id) (type*)MemoryCacheRealloc(ptr,sizeof(type)*(size),thread,id)
#define CacheFreeP(ptr,thread,id,force) {if(ptr) {MemoryCacheFree(ptr,thread,id,force);ptr=NULL;}}

#endif
