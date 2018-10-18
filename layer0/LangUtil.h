#pragma once

#include <memory>
#include <type_traits>

namespace pymol {

/*
 * @brief C++14's std::make_unique<T>
 * @param args arguments for constructor of Type T
 * @return A std::unique_ptr for type T
 */

template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <typename T, typename = typename std::enable_if<std::is_array<T>::value>::type>
std::unique_ptr<T> make_unique(std::size_t size){
  return std::unique_ptr<T>(new typename std::remove_extent<T>::type[size]());
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

}// namespace schrodinger


/*
 * Allows for operator |= on std::bitset::reference
 */

template<typename T, typename U, typename =
    typename std::enable_if<std::is_convertible<T, bool>::value &&
                            std::is_convertible<U, bool>::value>::type>
void operator|= (T l, const U r)
  { l = l | r; }
