

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
#ifndef _H_MemoryDebug
#define _H_MemoryDebug

#include <vector>

#include "os_std.h"
#include "PyMOLGlobals.h"


/* ================================================================ 
 * Don't touch below unless you know what you are doing */

typedef struct VLARec {
  ov_size size, unit_size;
  float grow_factor;
  bool auto_zero;
} VLARec;


/* NOTE: in VLACheck, rec is a zero based array index, not a record count */
#define VLACheck(ptr,type,rec) VLACheck2<type>(ptr, rec)

#define VLAlloc(type,init_size) (type*)VLAMalloc(init_size,sizeof(type),5,0)
#define VLACalloc(type,init_size) (type*)VLAMalloc(init_size,sizeof(type),5,1)
#define VLASize(ptr,type,size) VLASize2<type>(ptr,size)
#define VLASizeForSure(ptr,type,size) VLASizeForSure2<type>(ptr,size)

#define VLACopy(ptr,type) (type*)VLANewCopy(ptr);
#define VLAInsert(ptr,type,index,count) {ptr=(type*)VLAInsertRaw(ptr,index,count);}
#define VLADelete(ptr,type,index,count) {ptr=(type*)VLADeleteRaw(ptr,index,count);}

#if defined(_PYMOL_IOS) && !defined(_WEBGL)
#include "MemoryUsage.h"
extern "C" void fireMemoryWarning();
#define PYMOL_IOS_MEMORY_CHECK                                                 \
  if (pymol::memory_available() < sizeof(T) * num) {                           \
    fireMemoryWarning();                                                       \
    return nullptr;                                                            \
  }
#else
#define PYMOL_IOS_MEMORY_CHECK
#endif

namespace pymol
{
using ::free;

struct default_free {
  void operator()(void* ptr) const { free(ptr); }
};

template <typename T> T* malloc(size_t num)
{
  return (T*) ::malloc(num * sizeof(T));
}

template <typename T> T* calloc(size_t num)
{
  return (T*) ::calloc(num, sizeof(T));
}

template <typename T> T* realloc(T* ptr, size_t num)
{
  return (T*) ::realloc(ptr, num * sizeof(T));
}
} // namespace pymol

#define FreeP(ptr) {if(ptr) {mfree(ptr);ptr=NULL;}}
#define DeleteP(ptr) {if(ptr) {delete ptr;ptr=NULL;}}
#define DeleteAP(ptr) {if(ptr) {delete[] ptr;ptr=NULL;}}

void *VLAExpand(void *ptr, ov_size rec);        /* NOTE: rec is index (total-1) */
void *MemoryReallocForSure(void *ptr, size_t newSize);
void *MemoryReallocForSureSafe(void *ptr, size_t newSize, size_t oldSize);

void *VLADeleteRaw(void *ptr, int index, unsigned int count);
void *VLAInsertRaw(void *ptr, int index, unsigned int count);

void *VLAMalloc(ov_size init_size, ov_size unit_size, unsigned int grow_factor, int auto_zero); /*growfactor 1-10 */

void VLAFree(void *ptr);
void *VLASetSize(void *ptr, size_t newSize);
void *VLASetSizeForSure(void *ptr, size_t newSize);

size_t VLAGetSize(const void *ptr);
void *VLANewCopy(const void *ptr);
void MemoryZero(char *p, char *q);


#define mfree pymol::free
#define mstrdup strdup
#define ReallocForSure(ptr,type,size) (type*)MemoryReallocForSure(ptr,sizeof(type)*(size))


inline size_t VLAGetByteSize(const void *ptr) {
  const VLARec *vla = ((const VLARec *) ptr) - 1;
  return vla->size * vla->unit_size;
}

/**
 * Templated version of the `VLACopy` macro
 */
template <typename T>
T * VLACopy2(const T * vla) {
  return VLACopy((void*)vla, T);
}

/**
 * @brief std::vector version of VLACheck. Checks to see if index i is valid for insertion.
 *        If not, a resize will be attempted.
 * @param vec: vector whose size will be check for valid insertion at index i
 * @param i: index for position where an element may be inserted into vec
 * Note: Use of this function should be limited. Used for safe replacement of VLACheck
 * Note: This function can throw.
 */
template <typename T>
void VecCheck(std::vector<T> &vec, std::size_t i){
  if(i >= vec.size()){
    vec.resize(i + 1);
  }
}

/**
 * @brief Similar to VecCheck but constructs objects with arguments. Useful when
 * no default constructor available.
 * @tparam T type of elements in vector and also created during resize
 * @tparam Ts types of T's constructor arguments
 * @param vec: vector whose size will be check for valid insertion at index i
 * @param i: index for position where an element may be inserted into vec
 * @param args: constructor arguments of T
 */
template <typename T, typename... Ts>
void VecCheckEmplace(std::vector<T>& vec, std::size_t i, Ts... args)
{
  vec.reserve(i + 1);
  for (auto s = vec.size(); s <= i; s++) {
    vec.emplace_back(args...);
  }
}

template <typename T>
T* VLACheck2(T*& ptr, size_t pos) {
  if (pos >= ((VLARec*) ptr)[-1].size) {
    ptr = static_cast<T*>(VLAExpand(ptr, pos));
  }
  return ptr;
}

template <typename T>
void VLASize2(T*& ptr, size_t size) {
  ptr = static_cast<T*>(VLASetSize(ptr, size));
}

template <typename T>
void VLASizeForSure2(T*& ptr, size_t size) {
  ptr = static_cast<T*>(VLASetSizeForSure(ptr, size));
}

template <typename T>
void VLAFreeP(T*& ptr) {
  if (ptr) {
    VLAFree(ptr);
    ptr = nullptr;
  }
}

#endif

// vi:sw=2:expandtab
