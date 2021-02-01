/*
 * (c) Schrodinger, Inc.
 */

#pragma once

#include "MemoryDebug.h"
#include "pymol/memory.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <type_traits>

namespace pymol
{

/*
 * NOTE: Use this only for existing VLAs
 * otherwise, Use std::vector<T> instead!
 *
 * NOTE: safe RAII version of our PyMOL VLAs
 *
 * NOTE: allocation exceptions handled by VLA API
 */

/**
 * Variable sized array.
 *
 * Compatible with legacy VLA API (MemoryDebug.h) through implicit casts.
 *
 * Conceptually (mostly) equivalent to std::vector<T>. Differences are:
 * - no capacity(), but size() will grow in a similar fashion when calling
 *   check()
 */
template <typename T> class vla
{
#if defined(__clang__) || !defined(__GNUC__) || __GNUC__ >= 5
  // Not available in GCC 4.8
  static_assert(std::is_trivially_copyable<T>::value, "");
#endif

  T* m_vla = nullptr;

  void swap(vla<T>& other) noexcept { std::swap(m_vla, other.m_vla); }

public:
  // implicit NULL constructor
  vla() = default;

  // NULL assignment
  vla<T>& operator=(std::nullptr_t)
  {
    freeP();
    return *this;
  }

  /**
   * Takes ownership of a legacy VLA pointer
   * @param ptr Pointer to a legacy VLA which was constructed with VLAlloc (or a related function)
   */
  template <typename U> friend vla<U> vla_take_ownership(U* ptr);

  /**
   * Construct a new zero-initialized container
   * @param size the size of the container
   */
  explicit vla(std::size_t size) { m_vla = VLACalloc(T, size); }

  // constructor with initializer list
  // Empty list constructs a NULL VLA to be consistent with default constructor
  vla(std::initializer_list<T> init)
  {
    if (init.size() == 0)
      return;
    resize(init.size());
    std::copy(init.begin(), init.end(), this->begin());
  }

  // copy constructor
  vla(const vla<T>& other) { m_vla = VLACopy2<T>(other.m_vla); }

  // copy assignment
  vla<T>& operator=(const vla<T>& other)
  {
    vla<T> tmp(other);
    swap(tmp);
    return *this;
  }

  // move constructor
  vla(vla<T>&& other) noexcept { swap(other); }

  // move assignment
  vla<T>& operator=(vla<T>&& other) noexcept
  {
    swap(other);
    return *this;
  }
  // destructor
  ~vla() { freeP(); }

  /**
   * Returns pointer to the underlying array serving as element storage.
   */
  const T* data() const { return m_vla; }
  T* data() { return m_vla; }

  /// legacy VLA cast
  operator const T*() const { return m_vla; }

  /// legacy address-of VLA cast
  T** operator&() { return &m_vla; }

  // note: VS2015 fails in various situations if this is not a template
  /// legacy pointer arithmetic
  template <typename S> const T* operator+(S i) const { return m_vla + i; }
  /// legacy pointer arithmetic
  template <typename S> T* operator+(S i) { return m_vla + i; }

  // note: VS2015 32bit fails with "overloads have similar conversions" if this is not a template
  /**
   * Returns a reference to the element at specified location @a i. No bounds
   * checking is performed.
   * @param i position of the element to return
   */
  template <typename S> T& operator[](S i) { return m_vla[i]; }
  template <typename S> const T& operator[](S i) const { return m_vla[i]; }

  /// legacy VLA member access
  T* operator->() { return m_vla; }
  /// legacy VLA member access
  const T* operator->() const { return m_vla; }

  /**
   * checks whether the owned pointer is not nullptr
   */
  explicit operator bool() const { return m_vla; }

  /**
   * Resizes the container to contain @a newSize elements.
   * If this didn't manage any data yet, a new zero-initialized buffer will be
   * allocated. Otherwise new elements will only be zero-initialized if the
   * first allocation was done with VLACalloc.
   * @param newSize new size of the container
   */
  void resize(std::size_t newSize)
  {
    if (m_vla == nullptr) {
      m_vla = VLACalloc(T, newSize);
    } else {
      VLASize(m_vla, T, newSize);
    }
  }

  /**
   * Returns the number of elements in the container
   */
  std::size_t size() const
  {
    if (m_vla == nullptr) {
      return 0;
    }
    return VLAGetSize(m_vla);
  }

  /**
   * Grow the container size if index @a i is out-of-bounds.
   * @post size() > i
   * @return pointer to element @a i
   * @pre operator bool() != false
   */
  T* check(std::size_t i)
  {
    assert(m_vla != nullptr);
    VLACheck(m_vla, T, i);
    return m_vla + i;
  }

  /**
   * Erases all elements from the container and frees the underlying buffer.
   * @post size() == 0
   * @post operator bool() == false
   */
  void freeP()
  {
    if (m_vla != nullptr) {
      VLAFreeP(m_vla);
    }
  }

  T* begin() { return m_vla; }
  T* end() { return m_vla + size(); }
  const T* begin() const { return m_vla; }
  const T* end() const { return m_vla + size(); }

  /**
   * Unlike resize(), won't shrink the size, and unlike check(), can be
   * called on empty container.
   * @post size() >= count
   * @post operator bool() != false
   */
  void reserve(std::size_t count)
  {
    if (m_vla == nullptr) {
      m_vla = VLACalloc(T, count);
    } else if (count != 0) {
      VLACheck(m_vla, T, count - 1);
    }
  }

  /**
   * Inserts 'count' number of elements to the vla at a given index
   * @param index position of the vla where elements are added
   * @param count number of elements added
   */
  void insert(int index, int count)
  {
    m_vla = static_cast<T*>(VLAInsertRaw(m_vla, index, count));
  }

  /**
   * Removes 'count' number of elements to the vla at a given index
   * @param index position of the vla where elements are removed
   * @param count number of elements removed
   */
  void erase(int index, int count)
  {
    m_vla = static_cast<T*>(VLADeleteRaw(m_vla, index, count));
  }
};

template <typename U> vla<U> vla_take_ownership(U* ptr)
{
  vla<U> instance;
  instance.m_vla = ptr;
  return instance;
}
} // namespace pymol

template <typename T> pymol::vla<T> VLACopy2(const pymol::vla<T>& v)
{
  return v; // calls copy constructor
}

template <typename T> void VLACheck2(pymol::vla<T>& v, size_t pos)
{
  v.check(pos);
}

template <typename T> void VLASize2(pymol::vla<T>& v, size_t size)
{
  v.resize(size);
}

template <typename T> void VLAFreeP(pymol::vla<T>& v)
{
  v.freeP();
}

// vi:sw=2:expandtab
