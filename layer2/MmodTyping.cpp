/*
 * MacroModel atom typing
 *
 * (c) Schrodinger, Inc.
 */

#include "MmodTyping.h"
#include "AtomInfo.h"

/**
 * Get the MacroModel atom type.
 *
 * Pre-condition: ObjectMoleculeVerifyChemistry
 */
int getMacroModelAtomType(const AtomInfoType * ai) {

  // generic ions
  switch (ai->formalCharge) {
    case -2:
      switch (ai->protons) {
        case cAN_S:  return 114; // Sm Sulfide anion, S2-
        case cAN_O:  return 115; // Om Oxide anion, O2-
      }
      break;
    case -1:
      switch (ai->protons) {
        case cAN_C:  return  10; // CM Carbanion (C-)
        case cAN_O:  return  18; // OM O- (alkoxide, carboxylate)
        case cAN_H:  return  45; // H5 H- (Anion)
        case cAN_S:  return  51; // SM S- (thiolate anion)
        case cAN_Cl: return 102; // Cm Cl-
        case cAN_F:  return 104; // Fm F-
        case cAN_Br: return 105; // Bm Br-
        case cAN_I:  return 106; // Im I-
      }
      break;
    case 0:
      switch (ai->protons) {
        case cAN_Li: return  93; // L0 Lithium neutral
        case cAN_Mg: return  94; // M0 Magnesium neutral
      }
      break;
    case 1:
      switch (ai->protons) {
        case cAN_C:  return  11; // CP Carbocation (C+)
        case cAN_H:  return  44; // H4 H+ (Cation)
        case cAN_Li: return  65; // Li Li+
        case cAN_Na: return  66; // Na Na+
        case cAN_K:  return  67; // K0 K+
        case cAN_Rb: return  68; // Rb Rb+
        case cAN_Cs: return  69; // Cs Cs+
        case cAN_Cu: return  85; // c1 Cu+
        case cAN_S:  return 100; // SP S+
      }
      break;
    case 2:
      switch (ai->protons) {
        case cAN_Ca: return  70; // Ca Ca+2
        case cAN_Ba: return  71; // Ba Ba+2
        case cAN_Mg: return  72; // Mg Mg+2
        case cAN_Fe: return  79; // f2 Fe+2
        case cAN_Co: return  81; // o2 Co+2
        case cAN_Ni: return  83; // n2 Ni+2
        case cAN_Cu: return  86; // c2 Cu+2
        case cAN_Zn: return  87; // Zn Zn+2
      }
      break;
    case 3:
      switch (ai->protons) {
        case cAN_Fe: return  80; // f3 Fe+3
        case cAN_Co: return  82; // o3 Co+3
        case cAN_Ni: return  84; // n3 Ni+3
      }
      break;
  }

  // others
  switch (ai->protons) {
    case cAN_C:
      switch (ai->geom) {
        case cAtomInfoLinear:       return 1; // C1 Carbon - sp
        case cAtomInfoPlanar:       return 2; // C2 Carbon - sp2
        case cAtomInfoTetrahedral:  return 3; // C3 Carbon - sp3
      }
      return 14; // C0 Any Carbon
    case cAN_O:
      if ((ai->flags & cAtomFlag_solvent) && !ai->bonded) {
        return 19; // OW United atom H2O - Water
      }
      switch (ai->geom) {
        case cAtomInfoPlanar:       return 15; // O2 Oxygen - double bond
        case cAtomInfoTetrahedral:  return 16; // O3 Oxygen - single bond
      }
      return 23; // O0 Any oxygen;
    case cAN_N:
      switch (ai->geom) {
        case cAtomInfoLinear:
          return 24; // N1 Nitrogen - sp
        case cAtomInfoPlanar:
          switch (ai->formalCharge) {
            case -1: return 38; // NM N- - sp2
            case  1: return 31; // N4 N+ - sp2
          }
          return 25; // N2 Nitrogen - sp2
        case cAtomInfoTetrahedral:
          switch (ai->formalCharge) {
            case -1: return 39; // NP N- - sp2
            case  1: return 32; // N5 N+ - sp2
          }
          return 26; // N3 Nitrogen - sp3
      }
      return 40; // N0 Any Nitrogen
    case cAN_H:
      return 48; // H0 Any Hydrogen
    case cAN_S:
      switch (ai->geom) {
        case cAtomInfoPlanar: return 101; // S2 Sulfur - sp2
      }
      return 52; // S0 Any Sulfur
    case cAN_P:
      if (ai->geom == cAtomInfoTetrahedral) {
        switch (ai->valence) {
          case 3: return  53; // P3 Phosphorus, trivalent
          case 4: return 107; // P5 Phosphorus, pentavalent tetrahedral
        }
      }
      return 108; // P0 Any phosphorus
    case cAN_B:
      switch (ai->geom) {
        case cAtomInfoPlanar:       return 54; // B2 Boron - sp2
        case cAtomInfoTetrahedral:  return 55; // B3 Boron - sp3
      }
      return 103; // B0 Any boron
    case cAN_F:
      return 56; // F0 Fluorine
    case cAN_Cl:
      return 57; // Cl Chlorine
    case cAN_Br:
      return 58; // Br Bromine
    case cAN_I:
      return 59; // I0 Iodine
    case cAN_Si:
      return 60; // Si Silicon
    case cAN_Mn:
      switch (ai->formalCharge) {
        case 2: return 73; // M2 Mn+2
        case 3: return 74; // M3 Mn+3
        case 4: return 75; // M4 Mn+4
        case 5: return 76; // M5 Mn+5
        case 6: return 77; // M6 Mn+6
        case 7: return 78; // M7 Mn+7
      }
      break;
    case cAN_Se:
      return 112; // Se Selenium
    case 0:
      if (strcmp(ai->elem, "LP") == 0)
        return 63; // Lp Lone electron pair
      return 61; // Du Dummy atom
  }

  return 64; // 00 Any Atom
}
