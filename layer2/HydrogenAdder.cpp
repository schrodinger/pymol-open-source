/*
 * (c) Schrodinger, Inc.
 */

#include "os_std.h"

#include "CoordSet.h"
#include "ObjectMolecule.h"
#include "Selector.h"
#include "HydrogenAdder.h"
#include "Err.h"

/*
 * Add coordinates for atom `atm`.
 *
 * Pre-conditions:
 * - `atm` doesn't have coordinates yet in given coord set
 */
static
void AppendAtomVertex(CoordSet* cs, unsigned atm, const float* v)
{
  int idx = cs->NIndex++;
  VLACheck(cs->Coord, float, idx * 3 + 2);
  VLACheck(cs->IdxToAtm, int, idx);

  cs->IdxToAtm[idx] = atm;

  if (cs->Obj->DiscreteFlag) {
    cs->Obj->DiscreteAtmToIdx[atm] = idx;
    cs->Obj->DiscreteCSet[atm] = cs;
  } else {
    cs->AtmToIdx[atm] = idx;
  }

  copy3f(v, cs->coordPtr(idx));
}

/*
 * If `atm` has a planar (sp2) configuration, then write the plane's normal
 * vector to the `normal` out pointer and return true.
 *
 * Pre-conditions:
 * - Neighbors up-to-date
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

  int neighbor_atm, tmp;
  ITERNEIGHBORATOMS(I->Neighbor, atm, neighbor_atm, tmp) {
    if (h_fix && I->AtomInfo[neighbor_atm].isHydrogen())
      continue;

    // get neighbor coordinate
    int neighbor_idx = cs->atmToIdx(neighbor_atm);
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

/*
 * Calculate plausible coordinates for those neighbors of `atm` which don't
 * have coordinates yet.
 *
 * h_fix: also reposition hydrogens with existing coordinates.
 *
 * Returns the number of added/updated coordinates.
 *
 * Pre-conditions:
 * - Neighbors up-to-date
 *
 * Note: Similar to ObjectMoleculeFindOpenValenceVector ("Evolutionary
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

  int neighbor_atm, tmp;
  ITERNEIGHBORATOMS(I->Neighbor, atm, neighbor_atm, tmp) {
    if (n_present == 4)
      break;

    // get neighbor coordinate
    int neighbor_idx = cs->atmToIdx(neighbor_atm);
    if (neighbor_idx == -1 ||
        (h_fix && I->AtomInfo[neighbor_atm].isHydrogen())) {
      missing_atm[n_missing++] = neighbor_atm;
      continue;
    }

    const float* neighbor_coord = cs->coordPtr(neighbor_idx);

    // points away from center
    float* vvec = cbuf + 3 * n_present;
    subtract3f(neighbor_coord, center_coord, vvec);
    normalize3f(vvec);

    present_atm = neighbor_atm;
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

/*
 * Add hydrogens to selection
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

  ObjectMoleculeUpdateNeighbors(I);

  // add hydrogens (without coordinates)
  for (unsigned atm = 0; atm < n_atom_old; ++atm) {
    const auto ai = I->AtomInfo + atm;

    if (ai->isMetal())
      continue;

    if (!SelectorIsMember(G, ai->selEntry, sele))
      continue;

    int nneighbors = I->Neighbor[I->Neighbor[atm]];
    int nimplicit = ai->valence - nneighbors;

    if (nimplicit <= 0)
      continue;

    VLACheck(I->AtomInfo, AtomInfoType, I->NAtom + nimplicit - 1);
    VLACheck(I->Bond,     BondType,     I->NBond + nimplicit - 1);

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
  for (StateIterator iter(G, nullptr, cSelectorUpdateTableAllStates, I->NCSet);
      iter.next();) {
    CoordSet* cs = I->CSet[iter.state];
    if (cs)
      cs->extendIndices(I->NAtom);
  }

  I->invalidate(cRepAll, cRepInvBonds, state);
  ObjectMoleculeUpdateNeighbors(I);

  AtomInfoUniquefyNames(G,
      I->AtomInfo, n_atom_old,
      I->AtomInfo + n_atom_old, nullptr,
      I->NAtom - n_atom_old);

  // fill coordinates
  for (StateIterator iter(G, I->Setting, state, I->NCSet); iter.next();) {
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
