/*
 * (c) Schrodinger, Inc.
 */

#pragma once

#include <algorithm>

#include "MemoryDebug.h"

/*
 * Class to handle ownership transfer of VLAs
 *
 * Conceptually equivalent to std::unique_ptr<T, VLAFree>
 */
template <typename T>
class unique_vla_ptr {
  T * m_vla;

  void swap(unique_vla_ptr<T> &other) {
    std::swap(m_vla, other.m_vla);
  }

public:
  // destructor
  ~unique_vla_ptr() {
    VLAFreeP(m_vla);
  }

  // constructor
  unique_vla_ptr() : m_vla(NULL) {}
  unique_vla_ptr(T * vla) : m_vla(vla) {}

  // move constructor
  unique_vla_ptr(unique_vla_ptr<T> && other) : m_vla(NULL) {
    swap(other);
  }

  // move assignment
  unique_vla_ptr<T> & operator= (unique_vla_ptr<T> && other) {
    swap(other);
    return *this;
  }

  // data access
  const T * get() const { return m_vla; }
  operator const T * () const { return m_vla; }

  // checks whether the owned pointer is not NULL
  operator bool() const {
    return m_vla;
  }
};
