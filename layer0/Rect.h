#pragma once

namespace pymol
{

template<typename T>
struct Rect
{
  T x1{};
  T x2{};
  T y1{};
  T y2{};

  Rect() = default;
  Rect(T x1_, T x2_, T y1_, T y2_)
    : x1(x1_)
    , x2(x2_)
    , y1(y1_)
    , y2(y2_)
    {}

  /**
   * Determines if point is inside rect
   * @param x 2D point x coord
   * @param y 2D point y coord
   * @param proper disallows point to be on edge of rect
   */
  bool contains (T x, T y, bool proper = true) const noexcept {
    if (proper) {
      return x > x1 && x < x2 && y > y1 && y < y2;
    }
    return x >= x1 && x <= x2 && y >= y1 && y <= y2;
  }
};

} // namespace pymol

