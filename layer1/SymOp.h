/**
 * @file Symmetry operation stuff
 *
 * (c) Schrodinger, Inc.
 */

#pragma once

#include <string>

namespace pymol
{

/**
 * Symmetry operation corresponding to records like
 * "_struct_conn.ptnr1_symmetry".
 *
 * This is a trivial type and should be calloc'd, for compatibility with
 * pymol::vla<BondType>.
 */
struct SymOp {
  unsigned char index; ///< 0-based symmetry operation index (0-191)
  signed char x, y, z; ///< Offset in fractional space

  SymOp() = default;
  explicit SymOp(const char* code) { reset(code); }

  /// Update from a string. "1_555", "1" or "" will set all members to zero.
  bool reset(const char* code);

  /// Get string representation. The default is "1_555".
  std::string to_string() const;

  /// False for default, true for everything else
  operator bool() const { return index || x || y || z; }
};

} // namespace pymol
