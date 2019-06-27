#pragma once

#include <algorithm>

namespace pymol
{

/**
 * @brief C++17's version of std::clamp
 */
template<typename T>
const T& clamp(const T& value, const T& low, const T& high){
  return std::max(low, std::min(value, high));
}

} // namespace pymol
