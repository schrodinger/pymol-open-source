#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <cstring>
#include <cmath>
#include "os_python.h"
#include "PConv.h"
#include <catch2/catch.hpp>

namespace pymol {
namespace test {

// Checks whether obj is zero'd out (Struct of all PoD Types without non-default
// values are 0)
template <typename T> static bool isStructZero(const T &obj) {
  const auto size = sizeof(T);
  std::vector<char> buffer{size, 0};
  if (std::memcmp(buffer.data(), &obj, size)) {
    return false;
  }
  return true;
}

// Checks whether array arr is zeroed
template <typename T>
static bool isArrayZero(const T *arr, const std::size_t len) {
  auto bytelen = len * sizeof(T);
  std::vector<char> buffer(bytelen, 0);
  return std::memcmp(buffer.data(), arr, bytelen) == 0;
}

// Checks whether arrays are equal
template <typename T>
static bool isArrayEqual(const T *arr1, const T *arr2, const std::size_t len) {
  return std::equal(arr1, arr1 + len, arr2);
}

// Checks whether type has all special member functions
template <typename T> static bool hasAllSpecialMemberFunctions()
{
  return std::is_default_constructible<T>::value &&
         std::is_copy_constructible<T>::value &&
         std::is_copy_assignable<T>::value &&
         std::is_move_constructible<T>::value &&
         std::is_move_assignable<T>::value;
}

// Checks whether ptr is equal to nullptr
template <typename T> static bool isNullptr(const T *ptr) {
  return ptr == nullptr;
}
struct PYMOL_TEST_API {
  static PyObject *PYMOL_TEST_SUCCESS;
  static PyObject *PYMOL_TEST_FAILURE;
};
}; // namespace test
}; // namespace pymol

