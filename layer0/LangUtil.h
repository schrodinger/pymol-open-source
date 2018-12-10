#pragma once

#include <memory>
#include <type_traits>

namespace pymol {

/*
 * @brief C++14's std::make_unique<T>
 * @param args arguments for constructor of Type T
 * @return A std::unique_ptr for type T
 */

// C++14 Type Traits Helper
template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

template <typename T>
using remove_reference_t = typename std::remove_reference<T>::type;

template <typename T>
using remove_extent_t = typename std::remove_extent<T>::type;

template <typename T, typename U>
using forward_check_t = pymol::enable_if_t<std::is_same<pymol::remove_reference_t<T>, U>::value>;

template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <typename T, typename = pymol::enable_if_t<std::is_array<T>::value>>
std::unique_ptr<T> make_unique(std::size_t size){
  return std::unique_ptr<T>(new pymol::remove_extent_t<T>[size]());
}

/*
 * @brief C++17's std::destroy_at
 */
template <typename T> void destroy_at(T* p) { p->~T(); }

/*
 * @brief C++17's std::destroy
 */

template <typename T> void destroy(T* iter, T* end)
{
  for(; iter != end; ++iter) {
    destroy_at(std::addressof(*iter));
  }
}

/*
 * @brief C++17's version of std::clamp
 */
template<typename T>
const T& clamp(const T& value, const T& low, const T& high){
  return std::max(low, std::min(value, high));
}

} // namespace pymol
