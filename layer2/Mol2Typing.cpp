/*
 * Mol2 atom typing
 *
 * (c) Schrodinger, Inc.
 */

#include "os_std.h"

#include "Mol2Typing.h"
#include "AtomInfo.h"
#include "ObjectMolecule.h"

/**
 * atm: Atom index of a carbon atom with geom=3
 *
 * Return: True if atm has 3 neighbors which are all nitrogens with geom=3
 */
static bool isGuanidiniumCarbon(ObjectMolecule * obj, int atm) {
  int neighbor_count = 0;
  int charge = 0;

  for (auto const& item : AtomNeighbors(obj, atm)) {
    AtomInfoType const* neighbor = obj->AtomInfo.data() + item.atm;
    if (neighbor->protons != cAN_N || neighbor->geom != 3)
      return false;
    ++neighbor_count;
    charge += neighbor->formalCharge;
  }

  return neighbor_count == 3 && charge > 0;
}

/**
 * TODO: bond order 4 seems to be no guarantee for a ring, which is
 * required for aromaticity
 *
 * Return: always false, since PyMOL doesn't know enough about aromaticity
 */
static bool isAromaticAtom(ObjectMolecule * obj, int atm) {
#if 0
  for (auto const& neighbor : AtomNeighbors(obj, atm)) {
    if (obj->Bond[neighbor.bond].order == 4)
      return true;
  }
#endif

  return false;
}

/**
 * atm: Atom index of an oxygen atom
 *
 * Return: True if atom is part of a carboxylate or phosphate group
 */
static bool isCarboxylateOrPhosphateOxygen(ObjectMolecule * obj, int atm) {
  int o_count = 0, other_count = 0;

  auto const neighbors = AtomNeighbors(obj, atm);

  // must have only one neighbor
  if (neighbors.size() != 1)
    return false;

  // get that one neighbor as center of the acidic group
  atm = neighbors[0].atm;

  // check center atom
  AtomInfoType * ai = obj->AtomInfo + atm;
  if (!(ai->protons == cAN_C && ai->geom == 3) &&
      !(ai->protons == cAN_P && ai->geom == 4))
    return false;

  // iterate over neighbors of center atom
  for (auto const& item : AtomNeighbors(obj, atm)) {
    AtomInfoType const* neighbor = obj->AtomInfo.data() + item.atm;
    if (neighbor->protons == cAN_O)
      ++o_count;
    else
      ++other_count;
  }

  // carboxylate
  if (ai->protons == cAN_C)
    return (o_count == 2 && other_count == 1);

  // phosphate
  return (o_count == 4 && other_count == 0);
}

/**
 * atm: Atom index of a sulfur atom
 *
 * Return: Number of bound Oxygens if bound to two non-Oxygen atoms. Otherwise 0.
 */
static int sulfurCountOxygenNeighbors(ObjectMolecule * obj, int atm) {
  int o_count = 0, other_count = 0;

  for (auto const& item : AtomNeighbors(obj, atm)) {
    AtomInfoType const* neighbor = obj->AtomInfo.data() + item.atm;
    if (neighbor->protons == cAN_O)
      ++o_count;
    else
      ++other_count;
  }

  return (other_count == 2) ? o_count : 0;
}

/**
 * Get the Tripos Mol2 atom type
 *
 * Pre-condition: ObjectMoleculeVerifyChemistry
 */
const char * getMOL2Type(ObjectMolecule * obj, int atm) {
  auto G = obj->G;
  auto ai = obj->AtomInfo + atm;

  switch (ai->protons) {
    case cAN_C:
      switch (ai->geom) {
        case cAtomInfoLinear:
          return "C.1";
        case cAtomInfoPlanar:
          if (isAromaticAtom(obj, atm))
            return "C.ar";
          if (isGuanidiniumCarbon(obj, atm))
            return "C.cat";
          return "C.2";
        case cAtomInfoTetrahedral:
          return "C.3";
      }
      break;

    case cAN_N:
      switch (ai->geom) {
        case cAtomInfoLinear:
          return "N.1";
        case cAtomInfoPlanar:
          if ((ai->flags & cAtomFlag_polymer)
              && ai->name == G->lex_const.N)
            return "N.am";
          if (isAromaticAtom(obj, atm))
            return "N.ar";
          if (ai->valence == 2 && ai->formalCharge == 0)
            return "N.2";
          return "N.pl3";
        case cAtomInfoTetrahedral:
          return (ai->formalCharge == 1) ? "N.4" : "N.3";
      }
      break;

    case cAN_O:
      if (isCarboxylateOrPhosphateOxygen(obj, atm))
        return "O.co2";
      switch (ai->geom) {
        case cAtomInfoPlanar:       return "O.2";
        case cAtomInfoTetrahedral:  return "O.3";
      }
      break;

    case cAN_S:
      switch (sulfurCountOxygenNeighbors(obj, atm)) {
        case 1: return "S.o";
        case 2: return "S.o2";
      }
      switch (ai->geom) {
        case cAtomInfoPlanar:       return "S.2";
        case cAtomInfoTetrahedral:  return "S.3";
      }
      break;

    case cAN_P:
      if (ai->geom == 4)
        return "P.3";
      break;

    case cAN_Cr:
      if (ai->geom == 4)
        return "Cr.th";
      return "Cr.oh";

    case cAN_Co:
      return "Co.oh";
  }

  if (ai->protons >= 0 && ai->protons < ElementTableSize)
    return ElementTable[ai->protons].symbol;

  return "Du";
}
