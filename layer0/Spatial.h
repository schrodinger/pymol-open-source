#pragma once

#include <cstdint>

struct Offset2D
{
  std::int32_t x;
  std::int32_t y;
};

struct Extent2D
{
  std::uint32_t width;
  std::uint32_t height;
};

struct Rect2D
{
  Offset2D offset;
  Extent2D extent;
};
