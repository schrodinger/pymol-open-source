/*
 * (c) Schrodinger, Inc.
 */

#include "os_predef.h"

#include "CoordSet.h"
#include "Feedback.h"
#include "ObjectMolecule.h"
#include "ObjectMolecule3.h"
#include "Vector.h"

#include <cassert>
#include <cmath>
#include <unordered_map>
#include <vector>

/**
 * Find and enumerate molecules and return the mapping of atom indices to
 * molecule identifiers.
 *
 * A molecule is a set of atoms which is connected by bonds of non-zero order.
 */
static std::vector<int> ObjectMoleculeGetMolMappingVec(
    ObjectMolecule const& objmol)
{
  std::vector<int> mapping(size_t(objmol.NAtom));

  auto const inv = [](int x) { return -(x + 1); };

  auto const get_mol = [&mapping](int atm) {
    while (atm >= 0)
      atm = mapping[atm];
    return atm;
  };

  for (size_t i = 0; i != mapping.size(); ++i) {
    mapping[i] = inv(i);
    assert(mapping[i] < 0);
  }

  for (size_t b = 0; b < objmol.NBond; ++b) {
    auto const& bnd = objmol.Bond[b];

    // Exclude zero-order bonds, but include symop bonds
    if (bnd.order < 1) {
      continue;
    }

    auto mol0 = get_mol(bnd.index[0]);
    auto mol1 = get_mol(bnd.index[1]);

    assert(mol0 < 0);
    assert(mol1 < 0);
    assert(mapping[inv(mol0)] == mol0);

    if (mol0 != mol1) {
      mapping[inv(mol1)] = inv(mol0);
    }
  }

  for (auto& atm : mapping) {
    atm = get_mol(atm);
  }

  return mapping;
}

static std::unordered_map<int, std::vector<unsigned>>
ObjectMoleculeGetMolMappingMap(ObjectMolecule const& objmol)
{
  std::unordered_map<int, std::vector<unsigned>> molecules;
  auto const molmapvec = ObjectMoleculeGetMolMappingVec(objmol);
  for (unsigned atm = 0; atm != molmapvec.size(); ++atm) {
    molecules[molmapvec[atm]].push_back(atm);
  }
  return molecules;
}

/**
 * Unwrap periodic boundary conditions so that molecules don't jump.
 */
void ObjectMoleculePBCUnwrap(ObjectMolecule& objmol, bool const bymol)
{
  auto G = objmol.G;
  CoordSet const* cs_prev = nullptr;
  CoordSet* cs_curr = nullptr;

  bool sg_warning_shown = false;

  auto const molecules = ObjectMoleculeGetMolMappingMap(objmol);

  for (StateIndex_t state = 0; state != objmol.NCSet;
       ++state, cs_prev = cs_curr) {
    cs_curr = objmol.CSet[state];
    if (!cs_curr)
      continue;

    auto const* sym = cs_curr->getSymmetry();
    if (!sym || sym->Crystal.isSuspicious())
      continue;

    if (!sg_warning_shown) {
      auto sg = pymol::zstring_view(sym->spaceGroup());
      if (sg != "" && sg != "P 1" && sg != "P1") {
        PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
        " %s-Warning: Space group is not 'P 1'.\n", __func__ ENDFB(G);
        sg_warning_shown = true;
      }
    }

    CoordSetRealToFrac(cs_curr, &(sym->Crystal));

    if (!cs_prev)
      continue;

    if (bymol) /* by molecule */ {
      for (auto const& mol : molecules) {
        double center_prev[4] = {};
        double center_curr[4] = {};

        for (auto atm : mol.second) {
          auto const idx_prev = cs_prev->atmToIdx(atm);
          auto const idx_curr = cs_curr->atmToIdx(atm);
          if (idx_prev != -1) {
            pymol::add3(center_prev, cs_prev->coordPtr(idx_prev), center_prev);
            center_prev[3] += 1;
          }
          if (idx_curr != -1) {
            pymol::add3(center_curr, cs_curr->coordPtr(idx_curr), center_curr);
            center_curr[3] += 1;
          }
        }

        float offset[3] = {};
        for (int i = 0; i != 3; ++i) {
          center_prev[i] /= center_prev[3];
          center_curr[i] /= center_curr[3];
          offset[i] = std::round(center_curr[i] - center_prev[i]);
        }

        for (auto atm : mol.second) {
          auto const idx = cs_curr->atmToIdx(atm);
          if (idx != -1) {
            auto* coord = cs_curr->coordPtr(idx);
            pymol::subtract3(coord, offset, coord);
          }
        }
      }
    } else /* by atom */ {
      for (unsigned atm = 0; atm != objmol.NAtom; ++atm) {
        auto const idx_prev = cs_prev->atmToIdx(atm);
        auto const idx_curr = cs_curr->atmToIdx(atm);

        if (idx_prev == -1 || idx_curr == -1)
          continue;

        auto* coord_prev = cs_prev->coordPtr(idx_prev);
        auto* coord_curr = cs_curr->coordPtr(idx_curr);

        for (int i = 0; i != 3; ++i) {
          coord_curr[i] -= std::round(coord_curr[i] - coord_prev[i]);
        }
      }
    }
  }

  for (StateIndex_t state = 0; state != objmol.NCSet; ++state) {
    cs_curr = objmol.CSet[state];
    if (!cs_curr)
      continue;

    auto const* sym = cs_curr->getSymmetry();
    if (!sym || sym->Crystal.isSuspicious())
      continue;

    CoordSetFracToReal(cs_curr, &(sym->Crystal));
  }

  objmol.invalidate(cRepAll, cRepInvCoord, cStateAll);
}

/**
 * Wrap molecules into periodic boundary box.
 *
 * @param center Target box center in real space
 */
void ObjectMoleculePBCWrap(ObjectMolecule& objmol, float const* center)
{
  float center_auto[3];

  auto const molecules = ObjectMoleculeGetMolMappingMap(objmol);

  for (StateIndex_t state = 0; state != objmol.NCSet; ++state) {
    auto* cs = objmol.CSet[state];
    if (!cs)
      continue;

    auto const* sym = cs->getSymmetry();
    if (!sym || sym->Crystal.isSuspicious())
      continue;

    // Default center is the coordinate average of the first state
    if (!center) {
      pymol::meanNx3(cs->Coord.data(), cs->NIndex, center_auto);
      center = center_auto;
    }

    CoordSetRealToFrac(cs, &(sym->Crystal));

    // Apply inverse state matrix to center
    float centerX[3];
    if (cs->getPremultipliedMatrix()) {
      transform44d3f(ObjectStateGetInvMatrix(cs), center, centerX);
    } else {
      copy3(center, centerX);
    }

    // Center in fractional coordinates
    transform33f3f(sym->Crystal.realToFrac(), centerX, centerX);

    for (auto const& mol : molecules) {
      double molcenter[4] = {};

      for (auto atm : mol.second) {
        auto const idx = cs->atmToIdx(atm);
        if (idx != -1) {
          pymol::add3(molcenter, cs->coordPtr(idx), molcenter);
          molcenter[3] += 1;
        }
      }

      for (int i = 0; i != 3; ++i) {
        molcenter[i] = std::round(molcenter[i] / molcenter[3] - centerX[i]);
      }

      for (auto atm : mol.second) {
        auto const idx = cs->atmToIdx(atm);
        if (idx != -1) {
          auto* coord = cs->coordPtr(idx);
          pymol::subtract3(coord, molcenter, coord);
        }
      }
    }

    CoordSetFracToReal(cs, &(sym->Crystal));
  }

  objmol.invalidate(cRepAll, cRepInvCoord, cStateAll);
}
