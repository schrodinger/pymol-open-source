/**
 * @file Tools for manipulating SymOp bonds
 *
 * (c) Schrodinger, Inc.
 */

#pragma once

#include "PyMOLGlobals.h"
#include "SymOp.h"

#include <vector>

struct CoordSet;

namespace pymol
{
std::vector<std::pair<SymOp, float>> find_bond_symops(
    CoordSet const&, unsigned, unsigned, float tolerance = 1.05);
} // namespace pymol
