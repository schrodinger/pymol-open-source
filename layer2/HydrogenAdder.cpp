/*
 * (c) Schrodinger, Inc.
 */

#include "os_std.h"

#include "CoordSet.h"
#include "ObjectMolecule.h"
#include "Selector.h"
#include "HydrogenAdder.h"
#include "Err.h"

#include <cassert>

/**
 * Add coordinates for atom `atm`.
 *
 * @param atm Atom index (in ObjectMolecule::AtomInfo)
 * @param v Atom coordinates (`float[3]`)
 *
 * @pre `atm` doesn't have coordinates yet in given coord set
 */
static
void AppendAtomVertex(CoordSet* cs, unsigned atm, const float* v)
{
  assert(cs->atmToIdx(atm) == -1);

  int const idx = cs->NIndex;
  cs->setNIndex(idx + 1);

  cs->IdxToAtm[idx] = atm;

  if (cs->Obj->DiscreteFlag) {
    cs->Obj->DiscreteAtmToIdx[atm] = idx;
    cs->Obj->DiscreteCSet[atm] = cs;
  } else {
    cs->AtmToIdx[atm] = idx;
  }

  copy3f(v, cs->coordPtr(idx));
}

/**
 * If `atm` has a planar (sp2) configuration, then write the plane's normal
 * vector to the `normal` out pointer and return true.
 *
 * @param atm Atom index (in ObjectMolecule::AtomInfo)
 * @param[out] normal Normal vector (`float[3]` with length 1)
 * @param h_fix If true, then ignore hydrogen neighbors
 *
 * @pre Neighbors up-to-date
 */
static
bool get_planer_normal_cs(
    const ObjectMolecule* I,
    const CoordSet* cs, unsigned atm, float *normal,
    bool h_fix)
{
  int nOcc = 0;
  float occ[3 * 3];

  if (I->AtomInfo[atm].geom != cAtomInfoPlanar)
    return false;

  int idx = cs->atmToIdx(atm);
  if (idx == -1)
    return false;

  const float* center_coord = cs->coordPtr(idx);

  for (auto const& neighbor : AtomNeighbors(I, atm)) {
    if (h_fix && I->AtomInfo[neighbor.atm].isHydrogen())
      continue;

    // get neighbor coordinate
    int neighbor_idx = cs->atmToIdx(neighbor.atm);
    if (neighbor_idx == -1)
      continue;

    const float* neighbor_coord = cs->coordPtr(neighbor_idx);

    // points away from center
    float* vvec = occ + 3 * nOcc;
    subtract3f(neighbor_coord, center_coord, vvec);
    normalize3f(vvec);

    if (++nOcc == 3)
      // more doesn't make sence for a planar system
      break;
  }

  if (nOcc < 2)
    return false;

  cross_product3f(occ, occ + 3, normal);

  // avg of all three cross products
  if (nOcc > 2) {
    float v2[3];
    for (int offset = 0; offset < 6; offset += 3) {
      cross_product3f(occ + offset, occ + 6, v2);
      if (dot_product3f(normal, v2) < 0) {
        scale3f(v2, -1.f, v2);
      }
      add3f(normal, v2, normal);
    }
  }

  normalize3f(normal);

  return true;
}

/**
 * Calculate plausible coordinates for those neighbors of `atm` which don't
 * have coordinates yet.
 *
 * @param h_fix also reposition hydrogens with existing coordinates.
 *
 * @return Number of added/updated coordinates.
 *
 * @pre Neighbors up-to-date
 *
 * @note Similar to ::ObjectMoleculeFindOpenValenceVector ("Evolutionary
 * descendant", code duplication event)
 */
int ObjectMoleculeSetMissingNeighborCoords(
    ObjectMolecule* I, CoordSet* cs, unsigned atm, bool h_fix)
{
  auto G = I->G;
  int n_present = 0;
  float cbuf[4 * 3];
  int present_atm = -1;
  int missing_atm[4];
  int n_missing = 0;

  const AtomInfoType* ai = I->AtomInfo + atm;

  int idx = cs->atmToIdx(atm);
  if (idx == -1)
    return 0;

  const float* center_coord = cs->coordPtr(idx);

  for (auto const& neighbor : AtomNeighbors(I, atm)) {
    if (n_present == 4)
      break;

    // get neighbor coordinate
    int neighbor_idx = cs->atmToIdx(neighbor.atm);
    if (neighbor_idx == -1 ||
        (h_fix && I->AtomInfo[neighbor.atm].isHydrogen())) {
      missing_atm[n_missing++] = neighbor.atm;
      continue;
    }

    const float* neighbor_coord = cs->coordPtr(neighbor_idx);

    // points away from center
    float* vvec = cbuf + 3 * n_present;
    subtract3f(neighbor_coord, center_coord, vvec);
    normalize3f(vvec);

    present_atm = neighbor.atm;
    ++n_present;
  }

  if (n_missing == 0)
    // nothing to do
    return 0;

  int n_system = n_present;
  if (n_system == 0) {
    get_random3f(cbuf);
    ++n_system;
  }

  switch (ai->geom) {
    float t[3], z[3];

    // Tetrahedral system: 109.5 degree angles
    case cAtomInfoTetrahedral:
      switch (n_system) {
        case 1:
          get_system1f3f(cbuf, t, z);
          scale3f(cbuf, -0.334F, t);    // cos(-109.5)
          scale3f(z, 0.943F, z);        // sin( 109.5)
          add3f(z, t, cbuf + 3);
          normalize3f(cbuf + 3);
        case 2:
          add3f(cbuf, cbuf + 3, t);
          normalize3f(t);
          scale3f(t, -1.0F, t);
          cross_product3f(cbuf, cbuf + 3, z);
          normalize3f(z);
          scale3f(z, 1.41F, z);         // tan(109.5 / 2.0)
          add3f(t, z, cbuf + 6);
          normalize3f(cbuf + 6);
        case 3:
          add3f(cbuf, cbuf + 3, t);
          add3f(cbuf + 6, t, t);
          scale3f(t, -1.0F, cbuf + 9);
          normalize3f(cbuf + 9);
      }
      n_system = 4;
      break;

    // Planar system: 120.0 degree angles
    case cAtomInfoPlanar:
      switch (n_system) {
        case 1:
          if (present_atm >= 0 &&
              get_planer_normal_cs(I, cs, present_atm, t, h_fix)) {
            get_system2f3f(cbuf, t, z);
          } else {
            get_system1f3f(cbuf, t, z);
          }
          scale3f(cbuf, -0.500F, t);
          scale3f(z, 0.866F, z);        // sin( 120.0)
          add3f(z, t, cbuf + 3);
          normalize3f(cbuf + 3);
        case 2:
          add3f(cbuf, cbuf + 3, t);
          scale3f(t, -1.0F, cbuf + 6);
          normalize3f(cbuf + 6);
      }
      n_system = 3;
      break;

    // Linear system: 180.0 degree angles
    case cAtomInfoLinear:
      switch (n_system) {
        case 1:
          scale3f(cbuf, -1.0F, cbuf + 3);
          normalize3f(cbuf + 3);
      }
      n_system = 2;
      break;
  }

  if (n_missing > n_system - n_present) {
    n_missing = n_system - n_present;
  }

  // adding coordinates will invalidate pointer
  float center_coord_copy[3];
  copy3f(center_coord, center_coord_copy);
  center_coord = nullptr;

  for (int i = 0; i < n_missing; ++i) {
    float bondlength = AtomInfoGetBondLength(G,
        I->AtomInfo + atm,
        I->AtomInfo + missing_atm[i]);
    float* coord = cbuf + (n_present + i) * 3;
    scale3f(coord, bondlength, coord);
    add3f(coord, center_coord_copy, coord);

    if (h_fix && (idx = cs->atmToIdx(missing_atm[i])) != -1) {
      copy3f(coord, cs->coordPtr(idx));
    } else {
      AppendAtomVertex(cs, missing_atm[i], coord);
    }
  }

  return n_missing;
}

/**
 * Add hydrogens to selection
 *
 * @param sele Valid atom selection
 * @param state Object state (can be all (-1) or current (-2))
 *
 * @return False if `I` has no atoms in the selection or if the chemistry (atom
 * geometry and valence) can't be determined.
 */
int ObjectMoleculeAddSeleHydrogensRefactored(ObjectMolecule* I, int sele, int state)
{
  auto G = I->G;
  auto const n_atom_old = I->NAtom;

  bool seleFlag = false;
  for (unsigned atm = 0; atm < n_atom_old; atm++) {
    const auto ai = I->AtomInfo + atm;
    if (SelectorIsMember(G, ai->selEntry, sele)) {
      seleFlag = true;
      break;
    }
  }

  if (!seleFlag) {
    return true;
  }

  if (!ObjectMoleculeVerifyChemistry(I, state)) {
    ErrMessage(G, " AddHydrogens", "missing chemical geometry information.");
    return false;
  }

  // add hydrogens (without coordinates)
  for (unsigned atm = 0; atm < n_atom_old; ++atm) {
    const auto ai = I->AtomInfo + atm;

    if (ai->isMetal())
      continue;

    if (!SelectorIsMember(G, ai->selEntry, sele))
      continue;

    int nneighbors = AtomNeighbors(I, atm).size();
    int nimplicit = ai->valence - nneighbors;

    if (nimplicit <= 0)
      continue;

    I->AtomInfo.reserve(I->NAtom + nimplicit);
    I->Bond.reserve(I->NBond + nimplicit);

    for (int i = 0; i < nimplicit; ++i) {
      // bond
      auto bond = I->Bond + I->NBond++;
      BondTypeInit2(bond, atm, I->NAtom, 1);

      // atom
      auto atom = I->AtomInfo + I->NAtom++;
      atom->protons = cAN_H;
      atom->geom = cAtomInfoSingle;
      atom->valence = 1;
      ObjectMoleculePrepareAtom(I, atm, atom, /* uniquefy */ false);
    }
  }

  // grow index arrays
  ObjectMoleculeExtendIndices(I, cSelectorUpdateTableAllStates);

  I->invalidate(cRepAll, cRepInvBonds, state);

  AtomInfoUniquefyNames(G,
      I->AtomInfo, n_atom_old,
      I->AtomInfo + n_atom_old, nullptr,
      I->NAtom - n_atom_old);

  // fill coordinates
  for (StateIterator iter(I, state); iter.next();) {
    CoordSet* cs = I->CSet[iter.state];
    if (!cs)
      continue;

    for (unsigned idx = 0; idx < cs->NIndex; ++idx) {
      auto atm = cs->IdxToAtm[idx];
      if (atm >= n_atom_old)
        continue;

      const auto ai = I->AtomInfo + atm;
      if (!SelectorIsMember(G, ai->selEntry, sele))
        continue;

      ObjectMoleculeSetMissingNeighborCoords(I, cs, atm);
    }
  }

  I->invalidate(cRepAll, cRepInvAtoms, state);
  ObjectMoleculeSort(I);
  ObjectMoleculeUpdateIDNumbers(I);

  return true;
}
