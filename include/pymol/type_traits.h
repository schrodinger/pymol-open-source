#pragma once

#include <type_traits>

namespace pymol
{

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

// Non-STL
template <typename T, typename U>
using forward_check_t = pymol::enable_if_t<std::is_same<pymol::remove_reference_t<T>, U>::value>;

}
