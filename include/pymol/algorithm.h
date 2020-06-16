#pragma once

#include <algorithm>
#include <cassert>

namespace pymol
{

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
