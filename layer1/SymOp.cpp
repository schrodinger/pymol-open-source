/**
 * @file Symmetry operation stuff
 *
 * (c) Schrodinger, Inc.
 */

#include "SymOp.h"

#include <cassert>
#include <cstdio>

namespace pymol
{

bool SymOp::reset(const char* code)
{
  assert(code);

  auto n = std::sscanf(code, "%hhu_%c%c%c", &index, &x, &y, &z);

  index = (n > 0) ? (index - 1) : 0;

  if (n < 4) {
    x = y = z = 0;
  } else {
    x -= static_cast<signed char>('5');
    y -= static_cast<signed char>('5');
    z -= static_cast<signed char>('5');
  }

  return true;
}

std::string SymOp::to_string() const
{
  char code[8];
  std::snprintf(code, 8, "%u_%d%d%d", index + 1, x + 5, y + 5, z + 5);
  return code;
}

} // namespace pymol
