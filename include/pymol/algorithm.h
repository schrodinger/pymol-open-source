#pragma once

#include <algorithm>
#include <cassert>

namespace pymol
{

/**
 * C++20's std::erase_if
 */
template <class C, class Pred> void erase_if(C& c, Pred pred)
{
#if __cplusplus >= 202002L
  std::erase_if(c, pred);
#else
  c.erase(std::remove_if(c.begin(), c.end(), pred), c.end());
#endif
}

/**
 * @brief C++17's version of std::clamp
 */
template<typename T>
const T& clamp(const T& value, const T& low, const T& high){
  assert(low <= high);
  return std::max(low, std::min(value, high));
}

/**
 * @brief C++14's std::equal
 */

template <typename InIter1, typename InIter2>
bool equal(InIter1 first1, InIter1 last1, InIter2 first2)
{
  for (; first1 != last1; ++first1, ++first2) {
    if (*first1 != *first2) {
      return false;
    }
  }
  return true;
}

} // namespace pymol
