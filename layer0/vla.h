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
 *
 * Conceptually equivalent to std::vector<T>
 */

template <typename T> class vla
{
  T* m_vla = nullptr;

  void swap(vla<T>& other) noexcept { std::swap(m_vla, other.m_vla); }

public:
  // implicit NULL constructor
  vla(std::nullptr_t) {}

  // constructor -- takes ownership of pointer
  explicit vla(T* vla = nullptr) : m_vla(vla) {}

  // constructor with size
  explicit vla(std::size_t size) { m_vla = VLACalloc(T, size); }

  // constructor with size and default value
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
  operator const T*() const { return m_vla; }
  operator T*&() { return m_vla; }

  T** operator&() { return &m_vla; }

  const T* operator+(size_t i) const { return m_vla + i; }
  T* operator+(size_t i) { return m_vla + i; }

  T& operator[](std::size_t i) { return m_vla[i]; }
  const T& operator[](std::size_t i) const { return m_vla[i]; }
  T* operator->() { return m_vla; }
  const T* operator->() const { return m_vla; }

  // checks whether the owned pointer is not nullptr
  explicit operator bool() const { return m_vla; }

  // Memory Management
  void resize(std::size_t newSize)
  {
    if (m_vla == nullptr) {
      m_vla = VLACalloc(T, newSize);
    } else {
      VLASize(m_vla, T, newSize);
    }
  }
  std::size_t size() const
  {
    if (m_vla == nullptr) {
      return 0;
    }
    return VLAGetSize(m_vla);
  }

  T* check(std::size_t i)
  {
    VLACheck(m_vla, T, i);
    return m_vla + i;
  }

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
