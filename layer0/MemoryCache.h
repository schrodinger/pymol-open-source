

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

/* memory cache off */

#define VLACacheCheck(G,ptr,type,rec,t,i) VLACheck(ptr,type,rec)
#define VLACacheAlloc(G,type,initSize,t,i) VLAlloc(type,initSize)
#define VLACacheFreeP(G,ptr,t,i,f) VLAFreeP(ptr)
#define VLACacheSize(G,ptr,type,size,t,i) VLASize(ptr,type,size)
#define VLACacheSizeForSure(G,ptr,type,size,t,i) VLASizeForSure(ptr,type,size)
#define VLACacheExpand(G,ptr,rec,thread_index,i) VLAExpand(ptr,rec)

#define MemoryCacheReplaceBlock(G,g,o,n)

#define VLACacheMalloc(G,a,b,c,d,t,i) VLAMalloc(a,b,c,d)
#define VLACacheFree(G,p,t,i,f) VLAFree(p)

#define CacheAlloc(G,type,size,thread,id) pymol::malloc<type>(size)
#define CacheCalloc(G,type,size,thread,id) pymol::calloc<type>(size)
#define CacheRealloc(G,ptr,type,size,thread,id) pymol::realloc<type>(size)
#define CacheFreeP(G,ptr,thread,id,force) {if(ptr) {mfree(ptr);ptr=NULL;}}

#endif
