#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include "pymol/type_traits.h"
#include <cstdlib>
#include <functional>
#if __cplusplus >= 201703L
#include <numeric>
#endif

#include "pymol/type_traits.h"

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
 * @brief Checks whether two floating point values are nearly equal
 * @param first floating point number
 * @param second floating point number
 * @param epsilon rounding cutoff
 * @tparam T first floating point type
 * @tparam U second floating point type
 */

template <typename T, typename U, typename CommonT = std::common_type_t<T, U>>
bool almost_equal(T a, U b, CommonT epsilon = 1e-6)
{
  return std::abs(a - b) <= epsilon;
}

/**
 * @brief Checks whether two floating point ranges are nearly equal
 * @param first iterator to beginning of first f.p. range
 * @param second iterator to beginning of first f.p. range
 * @param epsilon rounding cutoff
 * @tparam InIter1 first floating point iterator type
 * @tparam InIter2 second floating point iterator type
 * @return true if all elements in the range are almost equal
 */

template <typename InIter1, typename InIter2>
bool almost_equal_n(
    InIter1 first1, std::size_t n, InIter2 first2, float epsilon = 1e-6)
{
  for (std::size_t n_count{0u}; n_count != n; ++first1, ++first2, ++n_count) {
    if (!almost_equal(*first1, *first2, epsilon)) {
      return false;
    }
  }
  return true;
}

namespace ranges
{
/**
 * @brief Checks whether an element in a container satisfies a given predicate
 * @param cont container to check.
 * @param pred predicate to satisfy
 * @tparam ContainerT type of container
 * @tparam Pred type of predicate
 */

template <typename RangeT, typename ValueT>
bool contains(const RangeT& cont, const ValueT& value)
{
  auto it = std::find(std::begin(cont), std::end(cont), value);
  return it != std::end(cont);
}

/**
 * @brief Checks whether an element in a container satisfies a given predicate
 * @param cont container to check.
 * @param pred predicate to satisfy
 * @tparam ContainerT type of container
 * @tparam Pred type of predicate
 */

template <typename RangeT, typename Pred>
bool contains_if(const RangeT& cont, Pred pred)
{
  auto it = std::find_if(std::begin(cont), std::end(cont), pred);
  return it != std::end(cont);
}

/**
 * @brief checks to see if two ranges are equal
 * @param first reference range
 * @param second other range
 * @tparam Range type
 * @return true if both ranges are equal
 */

template <typename RangeT1, typename RangeT2>
bool equal(const RangeT1& first, const RangeT2& second)
{
  auto range1Size = std::distance(std::begin(first), std::end(first));
  auto range2Size = std::distance(std::begin(second), std::end(second));
  if (range1Size != range2Size) {
    return false;
  }
  return std::equal(std::begin(first), std::end(first), std::begin(second));
}

/**
 * @brief checks to see if two ranges are equal
 * @param first reference range
 * @param second other range
 * @tparam Range type
 * @return true if both ranges are equal
 */

template <typename RangeT1, typename RangeT2, typename Pred>
bool equal(const RangeT1& first, const RangeT2& second, Pred p)
{
  auto range1Size = std::distance(std::begin(first), std::end(first));
  auto range2Size = std::distance(std::begin(second), std::end(second));
  if (range1Size != range2Size) {
    return false;
  }
  auto firstStart = std::begin(first);
  auto firstEnd = std::end(first);
  auto secondStart = std::begin(second);
  for (; firstStart != firstEnd; ++firstStart, ++secondStart) {
    if (!p(*firstStart, *secondStart)) {
      return false;
    }
  }
  return true;
  //return std::equal(std::begin(first), std::end(first), std::begin(second));
}

/**
 * @brief Performs a left fold over a ranges
 * @param range range to fold
 * @param init initial accumulation value
 * @param op binary operation to call per element in range
 * @return result accmulated value over range
 * @note modeled after C++23's proposed std::ranges::fold
 */

template<typename RangeT, typename T, typename BinaryOp>
T left_fold(const RangeT& range, T init = T{}, BinaryOp op = std::plus<T>())
{
#if __cplusplus >= 201703L
  return std::accumulate(std::begin(range), std::end(range), init, op);
#else
  for(auto it = std::begin(range); it != std::end(range); ++it)
  {
    init = op(init, *it);
  }
  return init;
#endif
}

template<typename RangeT, typename T = typename RangeT::value_type>
T left_fold(const RangeT& range, T init = T{})
{
  return left_fold(range, init, std::plus<T>());
}

} // namespace ranges
} // namespace pymol
