#pragma once

#include <cstdint>
#include <type_traits>

namespace pymol
{

/**
 * C++ 23's std::to_underlying
 * @param e The enum to convert
 * @return The underlying type of the enum
 */
template <class EnumT>
constexpr std::underlying_type_t<EnumT> to_underlying(EnumT e) noexcept
{
  return static_cast<std::underlying_type_t<EnumT>>(e);
}

} // namespace pymol
