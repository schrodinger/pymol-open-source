#pragma once

#include <array>

/**
 * Temporary location to house these typedefs
 */

namespace pymol
{
  using Mat3 = std::array<float, 9>;
  using Mat4 = std::array<float, 16>;
  using Vec3 = std::array<float, 3>;
}
