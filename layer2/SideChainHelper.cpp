/*
 * This file contains source code for the PyMOL computer program
 * Copyright (c) Schrodinger, LLC.
 */

#include <algorithm>

#include "os_std.h"

#include "SideChainHelper.h"
#include "Lex.h"

/**
 * Mark side_chain_helper atoms that are bonded to atoms without a visible
 * cartoon or ribbon
 */
void SideChainHelperMarkNonCartoonBonded(bool * marked,
    const ObjectMolecule * obj,
    const CoordSet * cs,
    bool cartoon_side_chain_helper,
    bool ribbon_side_chain_helper)
{
  auto G = obj->G;

  auto b = obj->Bond.data();
  auto b_end = b + obj->NBond;

  for(; b < b_end; ++b) {
    auto b1 = b->index[0];
    auto b2 = b->index[1];

    auto a1 = cs->atmToIdx(b1);
    auto a2 = cs->atmToIdx(b2);

    if (a1 < 0 || a2 < 0)
      continue;

    auto ati1 = obj->AtomInfo + b1;
    auto ati2 = obj->AtomInfo + b2;

    if (!(ati1->flags & ati2->flags & cAtomFlag_polymer))
      continue;

    if (!marked[b1]) {
      marked[b1] =
        ((ati1->visRep & cRepCartoonBit) && !(ati2->visRep & cRepCartoonBit) &&
         AtomSettingGetWD(G, ati1, cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper)) ||
        ((ati1->visRep & cRepRibbonBit) && !(ati2->visRep & cRepRibbonBit) &&
         AtomSettingGetWD(G, ati1, cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper));
    }

    if (!marked[b2]) {
      marked[b2] =
        ((ati2->visRep & cRepCartoonBit) && !(ati1->visRep & cRepCartoonBit) &&
         AtomSettingGetWD(G, ati2, cSetting_cartoon_side_chain_helper, cartoon_side_chain_helper)) ||
        ((ati2->visRep & cRepRibbonBit) && !(ati1->visRep & cRepRibbonBit) &&
         AtomSettingGetWD(G, ati2, cSetting_ribbon_side_chain_helper, ribbon_side_chain_helper));
    }
  }
}

// matches c0[35][*']
inline bool is_35prime(const char * p, char c0) {
  return p[0] == c0 &&
    (p[1] == '3' || p[1] == '5') &&
    (p[2] == '*' || p[2] == '\'') && !p[3];
}

// matches C[45][*']
inline bool is_C45prime(const char * p) {
  return p[0] == 'C' &&
    (p[1] == '4' || p[1] == '5') &&
    (p[2] == '*' || p[2] == '\'') && !p[3];
}

/**
 * Return true if bond is hidden with side_chain_helper.
 * c1/c2 are in-out variables for color transfer.
 */
bool SideChainHelperFilterBond(PyMOLGlobals * G,
    const bool *marked,
    const AtomInfoType *ati1,
    const AtomInfoType *ati2,
    int b1, int b2, int na_mode, int *c1, int *c2)
{
  if (ati1->protons == cAN_H ||
      ati2->protons == cAN_N ||
      ati2->protons == cAN_O ||
      (ati1->protons == cAN_C && ati2->protons == cAN_C && ati2->name == G->lex_const.CA)
     ) {
    std::swap(ati1, ati2);
    std::swap(b1, b2);
    std::swap(c1, c2);
  }

  const char *name1 = LexStr(G, ati1->name);
  int prot1 = ati1->protons;
  const char *name2 = LexStr(G, ati2->name);
  int prot2 = ati2->protons;

  switch (prot1) {
    case cAN_C:
      if (ati1->name == G->lex_const.CA) {  /* CA */
        if(prot2 == cAN_C) {
          if(ati2->name == G->lex_const.CB)
            *c1 = *c2;      /* CA-CB */
          else if(ati2->name == G->lex_const.C && (!marked[b2]))
            return true;  /* suppress CA-C */
        } else if(prot2 == cAN_H)
          return true;    /* suppress all CA-hydrogens */
      } else if((na_mode == 1) && (prot2 == cAN_C)) {
        if (is_C45prime(name2) && is_C45prime(name1))
          /* supress C[45][*']-C[45][*'] */
          return true;
      }
      break;
    case cAN_N:
      if (ati1->name == G->lex_const.N) {   /* N */
        if(prot2 == cAN_C) {
          if (ati2->name == G->lex_const.CD)
            *c1 = *c2;      /* N->CD in PRO */
          else if (ati2->name == G->lex_const.CA
              && (!marked[b1])) {
            if(ati2->resn != G->lex_const.PRO)
              return true;        /* suppress N-CA, except in pro */
            *c1 = *c2;
          } else if (ati2->name == G->lex_const.C && (!marked[b1]))
            return true;  /* suppress N-C */
        } else if(prot2 == cAN_H)
          return true;    /* suppress all N-hydrogens */
      }
      break;
    case cAN_O:
      if(prot2 == cAN_C) {
        if ( ati2->name == G->lex_const.C &&
            (ati1->name == G->lex_const.O ||
             ati1->name == G->lex_const.OXT)
            && (!marked[b2]))
          return true;      /* suppress C-O,OXT */
        else if(na_mode == 1) {
          if (is_35prime(name2, 'C') && is_35prime(name1, 'O'))
            /* supress O[35][*']-C[35][*'] */
            return true;
        }
      } else if(prot2 == cAN_P) {
        if (ati2->name == G->lex_const.P &&
            strlen(name1) == 3 && name1[0] == 'O' && (
              (name1[2] == 'P' && name1[1] >= '1' && name1[1] <= '3') ||
              (name1[1] == 'P' && name1[2] >= '1' && name1[2] <= '3')))
          /* suppress P-O([123]P|P[123]) */
          return true;
        else if(na_mode == 1) {
          if (ati2->name == G->lex_const.P && is_35prime(name1, 'O'))
            /* supress P-O[35][*'] */
            return true;
        }
      }
      break;
  }

  return false;
}
