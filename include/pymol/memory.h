#pragma once

#include <memory>

#include "pymol/type_traits.h"

namespace pymol
{

#if __cplusplus >= 201402L
using std::make_unique;
#else
/**
 * @brief C++14's std::make_unique<T>
 * @param args arguments for constructor of Type T
 * @return A std::unique_ptr for type T
 */

template <typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <typename T, typename = pymol::enable_if_t<std::is_array<T>::value>>
std::unique_ptr<T> make_unique(std::size_t size){
  return std::unique_ptr<T>(new pymol::remove_extent_t<T>[size]());
}
#endif

#if __cplusplus >= 201703L
using std::destroy_at;
using std::destroy;
#else
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
#endif

} // namespace pymol
