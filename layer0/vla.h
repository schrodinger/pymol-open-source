/*
 * (c) Schrodinger, Inc.
 */

#pragma once

#include "MemoryDebug.h"
#include "LangUtil.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

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
  T* m_vla = nullptr;

  void swap(vla<T>& other) noexcept { std::swap(m_vla, other.m_vla); }

public:
  // implicit NULL constructor
  vla(std::nullptr_t) {}

  /**
   * Takes ownership of a legacy VLA pointer
   * @param vla Pointer to a legacy VLA which was constructed with VLAlloc (or a related function)
   */
  explicit vla(T* vla = nullptr) : m_vla(vla) {}

  /**
   * Construct a new zero-initialized container
   * @param size the size of the container
   */
  explicit vla(std::size_t size) { m_vla = VLACalloc(T, size); }

  /**
   * Constructs the container with @a size copies of elements with value @a
   * value.
   * @param size the size of the container
   * @param value the value to initialize elements of the container with
   */
  vla(std::size_t size, T value)
  {
    m_vla = VLAlloc(T, size);
    for (size_t i = 0; i < size; ++i) {
      new (m_vla + i) T(value);
    }
  }

  // constructor with initializer list
  vla(std::initializer_list<T> init) : vla(init.size())
  {
    std::copy(init.begin(), init.end(), this->begin());
  }

  // constructor from std::vector
  explicit vla(std::vector<T>& vec) : vla(vec.size())
  {
    std::copy(vec.begin(), vec.end(), this->begin());
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

  // data access - Allows for some usability with VLA APIs
  const T* data() const { return m_vla; }
  T*& data() { return m_vla; }

  /// legacy VLA cast
  operator const T*() const { return m_vla; }
  /// legacy VLA cast
  operator T*&() { return m_vla; }

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
   */
  T* check(std::size_t i)
  {
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
      pymol::destroy(begin(), end());
      VLAFreeP(m_vla);
    }
  }

  // Util functions
  std::vector<T> toStdVector() const
  {
    return std::vector<T>(m_vla, m_vla + size());
  }

  T* begin() { return m_vla; }
  T* end() { return m_vla + size(); }
  const T* begin() const { return m_vla; }
  const T* end() const { return m_vla + size(); }
};
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
  VLAFreeP(v.data());
}

// vi:sw=2:expandtab
