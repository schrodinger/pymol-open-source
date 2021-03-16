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

/**
 * A copyable unique pointer. Will make a copy of the managed object with "new".
 */
template <typename T, class D = std::default_delete<T>>
class copyable_ptr : public std::unique_ptr<T, D>
{
  T* copy_ptr() const { return *this ? new T(**this) : nullptr; }

public:
  using std::unique_ptr<T, D>::unique_ptr;

  copyable_ptr() = default;
  // move
  copyable_ptr(copyable_ptr&& other) : copyable_ptr(other.release()) {}
  copyable_ptr& operator=(copyable_ptr&& other)
  {
    this->reset(other.release());
    return *this;
  }

  // copy
  copyable_ptr(const copyable_ptr &other) : copyable_ptr(other.copy_ptr()) {}
  copyable_ptr& operator=(const copyable_ptr& other)
  {
    this->reset(other.copy_ptr());
    return *this;
  }
};

template <typename T, typename... Args>
copyable_ptr<T> make_copyable(Args &&... args) {
  return copyable_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * A unique pointer which copies to nullptr.
 */
template <typename T, class D = std::default_delete<T>>
class cache_ptr : public std::unique_ptr<T, D>
{
public:
  using std::unique_ptr<T, D>::unique_ptr;

  cache_ptr() = default;
  // move
  cache_ptr(cache_ptr&& other) noexcept : cache_ptr(other.release()) {}
  cache_ptr& operator=(cache_ptr&& other) noexcept
  {
    this->reset(other.release());
    return *this;
  }

  // copy
  cache_ptr(const cache_ptr& other) {}
  cache_ptr& operator=(const cache_ptr& other)
  {
    this->reset();
    return *this;
  }
};

template <typename T, typename... Args>
cache_ptr<T> make_cache(Args &&... args) {
  return cache_ptr<T>(new T(std::forward<Args>(args)...));
}

/**
 * Take ownership of a raw pointer with a custom delete function.
 *
 * Example:
 *
 *     auto s = unique_ptr_take_ownership(strdup("Hello"), free);
 */
template <typename T, typename Deleter>
std::unique_ptr<T, Deleter> unique_ptr_take_ownership(T* ptr, Deleter func)
{
  return std::unique_ptr<T, Deleter>(ptr, func);
}

} // namespace pymol
