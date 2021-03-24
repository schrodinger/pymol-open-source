#pragma once

#include <cstddef>
#include <type_traits>

namespace pymol
{

template<typename T>
using decay_t = typename std::decay<T>::type;

template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

template <typename T>
using remove_reference_t = typename std::remove_reference<T>::type;

template <typename T>
using remove_cv_t = typename std::remove_cv<T>::type;

template <typename T>
using remove_cvref_t = remove_cv_t<remove_reference_t<T>>;

template <typename T>
using remove_extent_t = typename std::remove_extent<T>::type;

template <typename T, typename U>
using common_type_t = typename std::common_type<T, U>::type;

template <typename T>
using result_of_t = typename std::result_of<T>::type;

// Non-STL
template <typename T, typename U>
using forward_check_t = pymol::enable_if_t<std::is_same<pymol::remove_reference_t<T>, U>::value>;

//! Casts a pointer of type T to a pointer of type T[N]
template <size_t N, typename T> remove_reference_t<T (*)[N]> reshape(T* flat)
{
  return reinterpret_cast<T(*)[N]>(flat);
}

//! Casts a pointer of type T[N] to a pointer of type T
template <typename T> remove_extent_t<T>* flatten(T* shaped)
{
  return reinterpret_cast<remove_extent_t<T>*>(shaped);
}
}
