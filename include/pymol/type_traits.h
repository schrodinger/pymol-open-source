#pragma once

#include <cstddef>
#include <type_traits>

namespace pymol
{

// Non-STL
//! Casts a pointer of type T to a pointer of type T[N]
template <size_t N, typename T> std::remove_reference_t<T (*)[N]> reshape(T* flat)
{
  return reinterpret_cast<T(*)[N]>(flat);
}

//! Casts a pointer of type T[N] to a pointer of type T
template <typename T> std::remove_extent_t<T>* flatten(T* shaped)
{
  return reinterpret_cast<std::remove_extent_t<T>*>(shaped);
}
}
