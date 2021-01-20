/**
 * @file Tools for manipulating SymOp bonds
 *
 * (c) Schrodinger, Inc.
 */

#include "SymOpTools.h"

#include "CoordSet.h"
#include "SymOp.h"

#include "pymol/algorithm.h"

#include <algorithm>

namespace pymol
{
/**
 * Find all symmetry operations of a PBC (periodic boundary condition) bond.
 *
 * Guesses the correct symmetry operations based on the distance to
 * all possible symmetry mates. This will include all images that are within 5%
 * of the shortest distance (if tolerance is 1.05).
 *
 * @param cs Coordinate set which provides atom coordinates and symmetry
 * information
 * @param atm1 Atom index "from"
 * @param atm2 Atom index "to"
 * @param tolerance Distance tolerance (>= 1.0)
 * @return List of symmetry operations and distances for atm2, given that atm1
 * has symmetry operation 1_555.
 */
std::vector<std::pair<SymOp, float>> find_bond_symops(
    CoordSet const& cs, unsigned atm1, unsigned atm2, float tolerance)
{
  assert(tolerance >= 1);

  auto const* sym = cs.getSymmetry();

  if (!sym) {
    return {};
  }

  int const idx[] = {
      cs.atmToIdx(atm1),
      cs.atmToIdx(atm2),
  };

  if (idx[0] == -1 || idx[1] == -1) {
    return {};
  }

  auto const n_symmat = sym->getNSymMat();
  float const* const v_fr = cs.coordPtr(idx[0]);
  float v_to_buf[3];
  float dist_min = FLT_MAX;
  float dist_min_tol = FLT_MAX;

  std::vector<std::pair<SymOp, float>> symops;

  for (auto symop = SymOp("1_444"); symop.x < 2; ++symop.x) {
    for (symop.y = -1; symop.y < 2; ++symop.y) {
      for (symop.z = -1; symop.z < 2; ++symop.z) {
        for (symop.index = 0; symop.index < n_symmat; ++symop.index) {
          // atom which binds to its own symmetry mate must not be 1_555
          if (idx[0] == idx[1] && !symop) {
            continue;
          }

          auto const v_to = cs.coordPtrSym(idx[1], symop, v_to_buf);
          assert(v_to);

          auto const dist = diffsq3f(v_fr, v_to);
          if (dist_min_tol < dist) {
            continue;
          }

          if (dist_min > dist) {
            dist_min = dist;
            dist_min_tol = dist_min * tolerance * tolerance;

            using cr_t = decltype(symops)::const_reference;
            pymol::erase_if(symops,
                [dist_min_tol](cr_t sd) { return dist_min_tol < sd.second; });
          }

          symops.emplace_back(symop, dist);
        }
      }
    }
  }

  return symops;
}
} // namespace pymol
