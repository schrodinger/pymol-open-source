/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"Base.h"
#include"Parse.h"
#include"OOMac.h"
#include"Vector.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Map.h"
#include"Selector.h"
#include"ObjectMolecule.h"
#include"Ortho.h"
#include"Util.h"
#include"Matrix.h"
#include"Scene.h"
#include"P.h"
#include"PConv.h"
#include"Executive.h"
#include"Setting.h"
#include"Sphere.h"
#include"main.h"
#include"CGO.h"
#include"Editor.h"
#include"Sculpt.h"
#include"OVContext.h"
#include"OVOneToOne.h"
#include"OVLexicon.h"
#include"ListMacros.h"
#include"File.h"
#include "Lex.h"
#include "MolV3000.h"
#include "HydrogenAdder.h"

#ifdef _WEBGL
#endif

#define cMaxNegResi 100

#define ntrim ParseNTrim
#define nextline ParseNextLine
#define ncopy ParseNCopy
#define nskip ParseNSkip

#include <algorithm>
#include <array>
#include <cassert>
#include <set>
#include <unordered_map>
#include <vector>

#ifndef NO_MMLIBS
#include "mmpymolx.h"
#endif

static
CoordSet *ObjectMoleculeMMDStr2CoordSet(PyMOLGlobals * G, const char *buffer,
                                        AtomInfoType ** atInfoPtr, const char **restart);

static
void ObjectMoleculeTransformTTTf(ObjectMolecule * I, float *ttt, int state);

static
int ObjectMoleculeGetAtomGeometry(const ObjectMolecule * I, int state, int at);

static
void ObjectMoleculeInferHBondFromChem(ObjectMolecule * I);

/**
 * Order function for sorting bonds by atom indices (ignores other properties
 * like bond order).
 *
 * Pre-condition: index pairs are in order (index[0] < index[1])
 */
static
int BondTypeInOrder(PyMOLGlobals * G, const BondType * bonds, int i1, int i2) {
  auto lhs = bonds[i1].index;
  auto rhs = bonds[i2].index;
  return lhs[0] < rhs[0] || (lhs[0] == rhs[0] && lhs[1] < rhs[1]);
}

/**
 * Remove duplicated bonds. Disregards if two atoms would be bonded by two
 * bonds of different bond order (only one will be kept).
 */
static
void ObjectMoleculeRemoveDuplicateBonds(PyMOLGlobals * G, ObjectMolecule * I) {
  // make sure index pairs are in order
  for (int i = 0; i < I->NBond; ++i) {
    auto& bond = I->Bond[i];
    if (bond.index[0] > bond.index[1] && !bond.symop_2) {
      std::swap(bond.index[0], bond.index[1]);
    }
  }

  // get sorted indices
  int * sorted = pymol::malloc<int>(I->NBond);
  UtilSortIndexGlobals(G, I->NBond, I->Bond.data(), sorted,
      (UtilOrderFnGlobals *) BondTypeInOrder);

  // purge duplicated bonds and mark for removal
  for (int i = 0, j = -1; i < I->NBond; ++i) {
    int b = sorted[i];
    auto ptr = I->Bond + sorted[i];

    bool equal = false;
    if (j != -1) {
      const auto& other = I->Bond[j];
      if (ptr->index[0] == other.index[0] &&
          ptr->index[1] == other.index[1])
        equal = true;
    }

    if (!equal) {
      j = b;
    } else {
      AtomInfoPurgeBond(G, ptr);

      // mark for removal
      ptr->index[0] = ptr->index[1] = 0;
    }
  }

  FreeP(sorted);

  // remove purged bonds from array
  int j = 0;
  for (int i = 0; i < I->NBond; ++i) {
    auto& bond = I->Bond[i];
    if (bond.index[0] != bond.index[1]) {
      if (j != i)
        std::swap(I->Bond[j], bond);
      ++j;
    }
  }
  I->NBond = j;
  VLASize(I->Bond, BondType, I->NBond);
}

/**
 * Convert a discrete object to non-discrete. Will merge atoms from states
 * by identifiers. Atom level properties like color, b, q, etc. are lost
 * on the merged atoms.
 */
static
bool ObjectMoleculeSetNotDiscrete(PyMOLGlobals * G, ObjectMolecule * I) {
  if (!I->DiscreteFlag)
    return true;

  VLAFreeP(I->DiscreteAtmToIdx);
  VLAFreeP(I->DiscreteCSet);
  I->DiscreteFlag = false;

  for (int atm = 0; atm < I->NAtom; ++atm) {
    I->AtomInfo[atm].discrete_state = 0;
  }

  int * outdex;
  int * aindex = AtomInfoGetSortedIndex(G, I, I->AtomInfo, I->NAtom, &outdex);

  // merge duplicate atoms
  for (int i = 0, j = -1; i < I->NAtom; ++i) {
    int atm = aindex[i];
    auto ai = I->AtomInfo + atm;

    if (j == -1 || !AtomInfoMatch(G, ai, I->AtomInfo + j, false, false)) {
      j = atm;
    } else {
      ai->deleteFlag = true;
    }

    // remapping of merged indices
    outdex[atm] = j;
  }

  // point coordsets to merged atoms
  for (int state = 0; state < I->NCSet; ++state) {
    auto cs = I->CSet[state];
    if (!cs)
      continue;

    for (int idx = 0; idx < cs->NIndex; ++idx) {
      int atm = cs->IdxToAtm[idx];
      cs->IdxToAtm[idx] = outdex[atm];
    }
  }

  // point bonds to merged atoms
  for (int i = 0; i < I->NBond; ++i) {
    auto& bond = I->Bond[i];
    bond.index[0] = outdex[bond.index[0]];
    bond.index[1] = outdex[bond.index[1]];
  }

  AtomInfoFreeSortedIndexes(G, &aindex, &outdex);

  ObjectMoleculeRemoveDuplicateBonds(G, I);

  I->updateAtmToIdx();
  ObjectMoleculePurge(I);

  return true;
}

/**
 * Make a non-discrete object discrete, or vice versa.
 */
int ObjectMoleculeSetDiscrete(PyMOLGlobals * G, ObjectMolecule * I, int discrete)
{
  int state, idx, ao, an, ao1, ao2, an1, an2;
  int maxnatom, natom = I->NAtom, nbond = I->NBond;
  int *aostate2an = NULL;
  char *bondseen = NULL;
  BondType *bond;
  CoordSet *cs;

  if (!discrete) {
    if (!I->DiscreteFlag)
      return true;

    return ObjectMoleculeSetNotDiscrete(G, I);
  }

  if (I->DiscreteFlag)
    return true;

  // upper bound for number of discrete atoms
  maxnatom = I->NAtom * I->NCSet;

  // mapping (for bonds): atom_old -> atom_new
  ok_assert(1, aostate2an = pymol::malloc<int>(I->NAtom));
  ok_assert(1, bondseen = pymol::calloc<char >(I->NBond));

  // discrete setup
  I->DiscreteFlag = discrete;
  ok_assert(1, I->DiscreteAtmToIdx = pymol::vla<int>(maxnatom));
  ok_assert(1, I->DiscreteCSet = pymol::vla<CoordSet*>(maxnatom));

  // for all coordinate sets
  for (state = 0; state < I->NCSet; state++) {
    cs = I->CSet[state];

    if (!cs)
      continue;

    // init the atom_old -> atom_new array
    for (ao = 0; ao < I->NAtom; ao++)
      aostate2an[ao] = -1;

    // for all atoms in coordinate set
    for (idx = 0; idx < cs->NIndex; idx++) {
      ao = an = cs->IdxToAtm[idx];

      if (I->DiscreteCSet[ao]) {
        // seen before, have to copy
        an = natom++;

        VLACheck(I->AtomInfo, AtomInfoType, an);
        ok_assert(1, I->AtomInfo);

        AtomInfoCopy(G,
            I->AtomInfo + ao,
            I->AtomInfo + an);
        cs->IdxToAtm[idx] = an;
      }

      I->AtomInfo[an].discrete_state = state + 1; // 1-based :-(
      I->DiscreteCSet[an] = cs;
      I->DiscreteAtmToIdx[an] = idx;

      // mapping (for bonds): (atom_old, state) -> atom_new
      aostate2an[ao] = an;
    }

    // unused in discrete coordsets
    cs->AtmToIdx.clear();

    // for all old bonds
    for (idx = 0; idx < I->NBond; idx++) {
      bond = I->Bond + idx;

      // old atom indices
      ao1 = bond->index[0];
      ao2 = bond->index[1];

      // new atom indices
      an1 = aostate2an[ao1];
      an2 = aostate2an[ao2];

      // do both atoms exist in this coord set?
      if (an1 == -1 || an2 == -1)
        continue;

      if (bondseen[idx]) {
        // seen before, have to copy
        VLACheck(I->Bond, BondType, nbond);
        ok_assert(1, I->Bond);

        bond = I->Bond + (nbond++);
        AtomInfoBondCopy(G, I->Bond + idx, bond);
      } else {
        bondseen[idx] = true;
      }

      bond->index[0] = an1;
      bond->index[1] = an2;
    }
  }

  // not needed anymore
  mfree(aostate2an);
  mfree(bondseen);

  // update N* count fields
  I->NAtom = natom;
  I->NBond = nbond;

  // trim VLAs memory to actual size
  if (I->NBond)
    VLASize(I->Bond, BondType, I->NBond);
  if (I->NAtom)
    VLASize(I->AtomInfo, AtomInfoType, I->NAtom);

  I->setNDiscrete(I->NAtom);

  I->invalidate(cRepAll, cRepInvAll, -1);

  return true;

ok_except1:
  PRINTFB(G, FB_ObjectMolecule, FB_Errors)
    " %s: memory allocation failed\n", __func__ ENDFB(G);
  return false;
}

int ObjectMoleculeCheckFullStateSelection(ObjectMolecule * I, int sele, int state)
{
  PyMOLGlobals *G = I->G;
  int result = false;
  if((state >= 0) && (state < I->NCSet)) {
    const AtomInfoType *ai = I->AtomInfo.data();
    CoordSet *cs = I->CSet[state];
    if(cs) {
      int a;
      int at;
      result = true;
      for(a = 0; a < cs->NIndex; a++) {
        at = cs->IdxToAtm[a];
        if(!SelectorIsMember(G, ai[at].selEntry, sele)) {
          result = false;
          break;
        }
      }
    }
  }
  return result;
}

/* PARAMS
 *   ch -- ptr to empty str of length 'len'
 *  len -- str len 
 * RETURNS
 *   ch, with the caption in place
 * NOTES
 *   User owns the buffer so must clean up after it
 */
char *ObjectMolecule::getCaption(char * ch, int len) const
{
  auto I = this;
  int objState;
  int n = 0;
  int show_state = 0;
  int show_as_fraction = 0;
  const char *frozen_str = "";

  int state = ObjectGetCurrentState(I, false);
  int counter_mode = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_state_counter_mode);
  int frozen = SettingGetIfDefined_i(I->G, I->Setting.get(), cSetting_state, &objState);

  /* if frozen print (blue) STATE / NSTATES
   * if not frozen, print STATE/NSTATES
   * if beyond NSTATES, print * /NSTATES.
   */

  if(frozen) { /* frozen color */
    frozen_str = "\\789";
  } else {
    if (I->DiscreteFlag) { /* discrete states */
      frozen_str = "\\993";
    } else { /* normal case */
      frozen_str = "";
    }
  }

  switch(counter_mode) {
  case 0: /* off */
    show_state = show_as_fraction = 0;
    break;
  case 2: /* just state */
    show_state = 1;
    show_as_fraction = 0;
    break;
  case -1: /* fraction, full-on */
  case 1:
  default:
    show_state = show_as_fraction = 1;
    break;
  }

  /* bail on null string or no room */
  if (!ch || len==0)
    return NULL;

  ch[0] = 0;

  /* if the state is valid, setup the label */
  if(state >= 0) {
    if (state < I->NCSet) {
      const CoordSet *cs = I->CSet[state];
      if(cs) {
	if(show_state) {
	  if (show_as_fraction) {
	    if (strlen(cs->Name)) { 	  /* NAME */
	      n = snprintf(ch, len, "%s %s%d/%d", cs->Name, frozen_str, state+1, I->NCSet);
	    } 
	    else { /* no name */
	      n = snprintf(ch, len, "%s%d/%d", frozen_str, state+1, I->NCSet);
	    }
	  } else { /* not fraction */
	    if (strlen(cs->Name)) {
	      n = snprintf(ch, len, "%s %s%d", cs->Name, frozen_str, state+1);
	    } else { /* no name */
	      n = snprintf(ch, len, "%s%d", frozen_str, state+1);
	    }
	  }
	} else { /* no state */
	  n = snprintf(ch, len, "%s", cs->Name);
	}
      } else { /* no coord set, can be an N-state object missing a CSet */
      } 
    } else { /* state > NCSet, out of range due to other object or global setting */
      if(show_state) {
	if(show_as_fraction) {
	  n = snprintf(ch, len, "%s--/%d", frozen_str, I->NCSet);
	} else { /* no fraction */
	  n = snprintf(ch, len, "%s--", frozen_str);
	}
      }
    }
  } else if (state == -1) {
    // all states
    n = snprintf(ch, len, "%s*/%d", frozen_str, I->NCSet);
  } else {
    /* blank out the title if outside the valid # of states */
  }

  if (n > len)
    return NULL;

  return ch;
}


#define MAX_BOND_DIST 50

/* find sets of atoms with identical skeletons */

struct match_info {
  AtomInfoType *ai_a;
  AtomInfoType *ai_b;
  BondType *bi_a;
  BondType *bi_b;
  const int *nbr_a, *nbr_b;
  int *matched;

  // values: -1 .. 2
  using mark_type = signed char;

  std::vector<mark_type> atom_mark_a;
  std::vector<mark_type> atom_mark_b;
  std::vector<mark_type> bond_mark_a;
  std::vector<mark_type> bond_mark_b;
};

#define recmat3(x,y,z) \
  (recursive_match(u_a[0],u_b[x],x_a[0],x_b[x],mi) &&   \
   recursive_match(u_a[1],u_b[y],x_a[1],x_b[y],mi) &&   \
   recursive_match(u_a[2],u_b[z],x_a[2],x_b[z],mi))

#define recmat4(w,x,y,z)  \
  (recursive_match(u_a[0],u_b[w],x_a[0],x_b[w],mi) && \
   recursive_match(u_a[1],u_b[x],x_a[1],x_b[x],mi) && \
   recursive_match(u_a[2],u_b[y],x_a[2],x_b[y],mi) && \
   recursive_match(u_a[3],u_b[z],x_a[3],x_b[z],mi))

static void undo_match(int *start, match_info * mi)
{
  /* remove the matched atom set */
  int *match = mi->matched;
  while(match > start) {
    int a = match[-4];
    int b = match[-3];
    int aa = match[-2];
    int bb = match[-1];
    mi->atom_mark_a[a] = false;
    mi->atom_mark_b[b] = false;
    mi->bond_mark_a[aa] = false;
    mi->bond_mark_b[bb] = false;
    match -= 4;
  }
  mi->matched = start;
}

static int recursive_match(int a, int b, int aa, int bb, match_info * mi)
{
  if (mi->atom_mark_a[a] && mi->atom_mark_b[b]) {
    if (aa >= 0 && bb >= 0 && !(mi->bond_mark_a[aa] || mi->bond_mark_b[bb])) {
      /* note bond */
      mi->atom_mark_a[a] = true;
      mi->atom_mark_b[b] = true;
      mi->bond_mark_a[aa] = true;
      mi->bond_mark_b[bb] = true;
      mi->matched[0] = a;       /* atoms */
      mi->matched[1] = b;
      mi->matched[2] = aa;      /* bonds */
      mi->matched[3] = bb;
      mi->matched += 4;
    }
    return true;
  } else if (mi->atom_mark_a[a] != mi->atom_mark_b[b])
    return false;
  else if(mi->ai_a[a].protons != mi->ai_b[b].protons)
    return false;
  else {
    int ni_a = mi->nbr_a[a];
    int ni_b = mi->nbr_b[b];
    int num_a = mi->nbr_a[ni_a++];
    int num_b = mi->nbr_b[ni_b++];

    if((num_a != num_b) || (num_a > 4)) /* must have the same number of neighbors (four or less) */
      return false;
    else {
      int match_flag = false;
      int u_a[4], u_b[4], x_a[4], x_b[4];
      int n_u_a = 0;
      int n_u_b = 0;
      int *before = mi->matched;
      mi->atom_mark_a[a] = true;
      mi->atom_mark_b[b] = true;
      mi->bond_mark_a[aa] = true;
      mi->bond_mark_b[bb] = true;
      mi->matched[0] = a;       /* atoms */
      mi->matched[1] = b;
      mi->matched[2] = aa;      /* bonds */
      mi->matched[3] = bb;
      mi->matched += 4;

      /* collect unvisited bonds */

      while(num_a--) {
        int i_a = mi->nbr_a[ni_a + 1];
        if(!mi->bond_mark_a[i_a]) {      /* bond not yet visited */
          u_a[n_u_a] = mi->nbr_a[ni_a]; /* atom index */
          x_a[n_u_a++] = i_a;
        }
        ni_a += 2;
      }

      while(num_b--) {
        int i_b = mi->nbr_b[ni_b + 1];
        if(!mi->bond_mark_b[i_b]) {
          u_b[n_u_b] = mi->nbr_b[ni_b];
          x_b[n_u_b++] = i_b;
        }
        ni_b += 2;
      }
      /* NOTE: implicitly relying upon C short-circuit evaluation */
      if(n_u_a == n_u_b) {
        if(!n_u_a)
          match_flag = true;
        else if(n_u_a == 1) {
          match_flag = recursive_match(u_a[0], u_b[0], x_a[0], x_b[0], mi);
        } else if(n_u_a == 2) {
          match_flag =
            (recursive_match(u_a[0], u_b[0], x_a[0], x_b[0], mi) &&
             recursive_match(u_a[1], u_b[1], x_a[1], x_b[1], mi)) ||
            (recursive_match(u_a[0], u_b[1], x_a[0], x_b[1], mi) &&
             recursive_match(u_a[1], u_b[0], x_a[1], x_b[0], mi));
        } else if(n_u_a == 3) {
          match_flag =
            recmat3(0, 1, 2) ||
            recmat3(0, 2, 1) ||
            recmat3(1, 0, 2) || recmat3(1, 2, 0) || recmat3(2, 0, 1) || recmat3(2, 1, 0);
        } else if(n_u_a == 4) {
          match_flag =
            recmat4(0, 1, 2, 3) ||
            recmat4(0, 2, 1, 3) ||
            recmat4(1, 0, 2, 3) ||
            recmat4(1, 2, 0, 3) || recmat4(2, 0, 1, 3) || recmat4(2, 1, 0, 3);
          if(!match_flag) {
            match_flag =
              recmat4(0, 1, 3, 2) ||
              recmat4(0, 2, 3, 1) ||
              recmat4(1, 0, 3, 2) ||
              recmat4(1, 2, 3, 0) || recmat4(2, 0, 3, 1) || recmat4(2, 1, 3, 0);
            if(!match_flag) {
              match_flag =
                recmat4(0, 3, 1, 2) ||
                recmat4(0, 3, 2, 1) ||
                recmat4(1, 3, 0, 2) ||
                recmat4(1, 3, 2, 0) || recmat4(2, 3, 0, 1) || recmat4(2, 3, 1, 0);
              if(!match_flag) {
                match_flag =
                  recmat4(3, 0, 1, 2) ||
                  recmat4(3, 0, 2, 1) ||
                  recmat4(3, 1, 0, 2) ||
                  recmat4(3, 1, 2, 0) || recmat4(3, 2, 0, 1) || recmat4(3, 2, 1, 0);
              }
            }
          }
        }
      }
      if(!match_flag) {         /* clean the match tree */
        undo_match(before, mi);
      }
      return match_flag;
    }
  }
}

int ObjectMoleculeXferValences(ObjectMolecule * Ia, int sele1, int sele2,
                               int target_state, ObjectMolecule * Ib, int sele3,
                               int source_state, int quiet)
{
  int *matched = NULL;
  int match_found = false;
  PyMOLGlobals *G = Ia->G;
  if(Ia == Ib)
    return false;

  {
    int max_match = Ia->NAtom + Ia->NBond;
    if(max_match < (Ib->NAtom + Ib->NBond))
      max_match = (Ib->NAtom + Ib->NBond);
    matched = pymol::calloc<int>(max_match * 4);
  }

  {
    int a, b;
    AtomInfoType *ai_a = Ia->AtomInfo.data();
    AtomInfoType *ai_b = Ib->AtomInfo.data();
    BondType *bi_a = Ia->Bond.data();
    BondType *bi_b = Ib->Bond.data();

    match_info mi;

    mi.atom_mark_a.resize(Ia->NAtom);
    mi.atom_mark_b.resize(Ib->NAtom);
    mi.bond_mark_a.resize(Ia->NBond);
    mi.bond_mark_b.resize(Ib->NBond);

    // confirm zero-initialization
    assert(std::none_of(mi.atom_mark_a.begin(), mi.atom_mark_a.end(),
        [](bool m) { return m; }));

    mi.ai_a = ai_a;
    mi.ai_b = ai_b;
    mi.bi_a = bi_a;
    mi.bi_b = bi_b;
    mi.nbr_a = Ia->getNeighborArray();
    mi.nbr_b = Ib->getNeighborArray();
    mi.matched = matched;
    for(a = 0; a < Ia->NAtom; a++) {
      if(!mi.atom_mark_a[a]) {
        int a_entry = ai_a[a].selEntry;
        if(SelectorIsMember(G, a_entry, sele1) || SelectorIsMember(G, a_entry, sele2)) {
          for(b = 0; b < Ib->NAtom; b++) {
            if(SelectorIsMember(G, ai_b[b].selEntry, sele3)) {
              if(recursive_match(a, b, -1, -1, &mi)) {
                int *match = mi.matched;
                match_found = true;
                /* now graft bond orders for bonds in matching atoms */

                while(match > matched) {
                  int at_a = match[-4];
                  int at_b = match[-3];
                  int bd_a = match[-2];
                  int bd_b = match[-1];

                  if((bd_a >= 0) && (bd_a >= 0)) {
                    int a1 = bi_a[bd_a].index[0];
                    int a2 = bi_a[bd_a].index[1];

                    int a1_entry = ai_a[a1].selEntry;
                    int a2_entry = ai_a[a2].selEntry;

                    if((SelectorIsMember(G, a1_entry, sele1) &&
                        SelectorIsMember(G, a2_entry, sele2)) ||
                       (SelectorIsMember(G, a2_entry, sele1) &&
                        SelectorIsMember(G, a1_entry, sele2))) {
                      /* only update bonds which actually sit between the two selections */
                      bi_a[bd_a].order = bi_b[bd_b].order;
                      ai_a[at_a].chemFlag = false;
                    }
                  }
                  mi.atom_mark_b[at_b] = false;     /* release source for future matching */
                  if(bd_b >= 0)
                    mi.bond_mark_b[bd_b] = false;
                  match -= 4;
                }
                break;
              }
            }
          }
        }
      }
    }
  }
  FreeP(matched);
  return match_found;
}

void ObjectMoleculeTransformState44f(ObjectMolecule * I, int state, const float *matrix,
                                     int log_trans, int homogenous, int transformed)
{
  int a;
  float tmp_matrix[16];
  CoordSet *cs;
  int use_matrices = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);
  if(use_matrices<0) use_matrices = 0;
  if(!use_matrices) {
    ObjectMoleculeTransformSelection(I, state, -1, matrix, log_trans, I->Name,
                                     homogenous, true);
  } else {
    double dbl_matrix[16];
    if(state == -2)
      state = ObjectGetCurrentState(I, false);
    /* ensure homogenous matrix to preserve programmer sanity */
    if(!homogenous) {
      convertTTTfR44d(matrix, dbl_matrix);
      copy44d44f(dbl_matrix, tmp_matrix);
      matrix = tmp_matrix;
    } else {
      copy44f44d(matrix, dbl_matrix);
    }

    if(state < 0) {             /* all states */
      for(a = 0; a < I->NCSet; a++) {
        cs = I->CSet[a];
        if(cs)
          ObjectStateLeftCombineMatrixR44d(cs, dbl_matrix);
      }
    } else if(state < I->NCSet) {       /* single state */
      cs = I->CSet[state];
      if(cs)
        ObjectStateLeftCombineMatrixR44d(cs, dbl_matrix);
    } else if(I->NCSet == 1) {  /* static singleton state */
      cs = I->CSet[0];
      if(cs && SettingGet_b(I->G, I->Setting.get(), NULL, cSetting_static_singletons)) {
        ObjectStateLeftCombineMatrixR44d(cs, dbl_matrix);
      }
    }
  }
}


/*========================================================================*/
static int ObjectMoleculeFixSeleHydrogens(ObjectMolecule * I, int sele, int state)
{
  int a;
  int seleFlag = false;
  const AtomInfoType *ai0;
  int ok = true;

  ai0 = I->AtomInfo;
  for(a = 0; a < I->NAtom; a++) {
    if(SelectorIsMember(I->G, ai0->selEntry, sele)) {
      seleFlag = true;
      break;
    }
    ai0++;
  }
  if(seleFlag) {
    seleFlag = false;
    if(!ObjectMoleculeVerifyChemistry(I, state)) {
      ErrMessage(I->G, " AddHydrogens", "missing chemical geometry information.");
    } else {
      ai0 = I->AtomInfo;
      for(a = 0; a < I->NAtom; a++) {
        if(!ai0->isHydrogen()) {    /* only do heavies */
          if(SelectorIsMember(I->G, ai0->selEntry, sele)) {
            for(StateIterator iter(I->G, I->Setting.get(), state, I->NCSet);
                iter.next();) {
              auto cs = I->CSet[iter.state];
              if (!cs)
                continue;

              seleFlag |= ObjectMoleculeSetMissingNeighborCoords(I, cs, a, true);
            }
          }
        }
        ai0++;
      }
    }
    if(seleFlag)
      I->invalidate(cRepAll, cRepInvAll, -1);
  }
  return ok;
}

static const char *skip_fortran(int num, int per_line, const char *p)
{
  int a, b;
  b = 0;
  for(a = 0; a < num; a++) {
    if((++b) == per_line) {
      b = 0;
      p = nextline(p);
    }
  }
  if(b || (!num))
    p = nextline(p);
  return (p);
}

void ObjectMoleculeOpRecInit(ObjectMoleculeOpRec * op)
{
  UtilZeroMem((char *) op, sizeof(ObjectMoleculeOpRec));
}

/**
 * Returns the most proton-rich element with the lowest priority value (OG1
 * before OG2, CG, HB1)
 */
int ObjectMoleculeGetTopNeighbor(PyMOLGlobals * G,
                                 ObjectMolecule * I, int start, int excluded)
{
  int highest_at = -1, highest_prot = 0, lowest_pri = 9999;

  for (auto const& neighbor : AtomNeighbors(I, start)) {
    auto const at = neighbor.atm;
    auto const *ai = I->AtomInfo.data() + at;
    if((highest_at < 0) && (at != excluded)) {
      highest_prot = ai->protons;
      lowest_pri = ai->priority;
      highest_at = at;
    } else if(((ai->protons > highest_prot) ||
               ((ai->protons == highest_prot) && (ai->priority < lowest_pri)))
              && (at != excluded)) {
      highest_prot = ai->protons;
      highest_at = at;
      lowest_pri = ai->priority;
    }
  }
  return highest_at;
}

/**
 * Selection mapping for `load_traj` to only load coordinates for a subset of
 * atoms. Returns a mapping of old to new coordinate indices and modifies the
 * index-related members of `cs`.
 *
 * @param[in] obj Object Molecule for evaluation of the selection
 * @param[in,out] cs Candidate coordinate set (with valid index arrays)
 * @param[in] selection Named selection
 * @return index mapping, or NULL if `selection` is not valid
 */
std::unique_ptr<int[]> LoadTrajSeleHelper(
    const ObjectMolecule* obj, CoordSet* cs, const char* selection)
{
  auto G = obj->G;
  int sele0 = SelectorIndexByName(G, selection);

  // We can skip "all" (0) here since it's a trivial case
  if (sele0 <= 0) {
    return nullptr;
  }

  auto xref = std::unique_ptr<int[]>(new int[cs->NIndex]);

  int idx_new = 0;

  for (int idx = 0; idx < cs->NIndex; ++idx) {
    auto atm = cs->IdxToAtm[idx];
    if (SelectorIsMember(G, obj->AtomInfo[atm].selEntry, sele0)) {
      cs->IdxToAtm[idx_new] = atm;
      cs->AtmToIdx[atm] = idx_new;
      xref[idx] = idx_new;
      ++idx_new;
    } else {
      cs->AtmToIdx[atm] = -1;
      xref[idx] = -1;
    }
  }

  cs->NIndex = idx_new;
  cs->IdxToAtm.resize(cs->NIndex);
  cs->Coord.resize(cs->NIndex * 3);

  return xref;
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadTRJFile(PyMOLGlobals * G, ObjectMolecule * I,
                                          const char *fname, int frame, int interval,
                                          int average, int start, int stop, int max,
                                          const char *sele, int image, const float *shift, int quiet)
{
  FILE *f;
  char *buffer;
  const char *p;
  char cc[MAXLINELEN];
  int n_read;
  int to_go;
  int skip_first_line = true;
  int periodic = false;
  int angles = true;
  float f0, f1, f2, *fp;
  float box[3], angle[3];
  float r_cent[3], r_trans[3];
  int r_act, r_val, r_cnt;
  float *r_fp_start = NULL, *r_fp_stop = NULL;
  int a, b, c, i;
  int zoom_flag = false;
  int cnt = 0;
  int n_avg = 0;
  int icnt;
  int ncnt = 0;
  float zerovector[3] = { 0.0, 0.0, 0.0 };
  CoordSet *cs = NULL;

  if (I->DiscreteFlag) {
    PRINTFB(G, FB_ObjectMolecule, FB_Errors)
      " %s: Discrete objects not supported\n", __func__ ENDFB(G);
    return I;
  }

  if(!shift)
    shift = zerovector;
  if(interval < 1)
    interval = 1;

  icnt = interval;
#define BUFSIZE 4194304
#define GETTING_LOW 10000

  f = pymol_fopen(fname, "rb");
  if(!f) {
    ErrMessage(G, __func__, "Unable to open file!");
  } else {
    if(I->CSTmpl) {
      cs = CoordSetCopy(I->CSTmpl);
    } else if (I->NCSet > 0) {
      cs = CoordSetCopy(I->CSet[0]);
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " ObjMolLoadTRJFile: Missing topology" ENDFB(G);
      return (I);
    }

    auto xref = LoadTrajSeleHelper(I, cs, sele);

    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " ObjMolLoadTRJFile: Loading from \"%s\".\n", fname ENDFB(G);
    buffer = pymol::malloc<char>(BUFSIZE + 1);     /* 1 MB read buffer */
    p = buffer;
    buffer[0] = 0;
    n_read = 0;
    to_go = 0;
    a = 0;
    b = 0;
    c = 0;
    f1 = 0.0;
    f2 = 0.0;
    while(1) {
      to_go = n_read - (p - buffer);
      if(to_go < GETTING_LOW)
        if(!feof(f)) {
          if(to_go)
            memcpy(buffer, p, to_go);
          n_read = fread(buffer + to_go, 1, BUFSIZE - to_go, f);
          n_read = to_go + n_read;
          buffer[n_read] = 0;
          p = buffer;
          if(skip_first_line) {
            p = nextline(p);
            skip_first_line = false;
          }
          to_go = n_read - (p - buffer);
        }
      if(!to_go)
        break;
      p = ncopy(cc, p, 8);
      if((++b) == 10) {
        b = 0;
        p = nextline(p);
      }
      f0 = f1;
      f1 = f2;
      if(sscanf(cc, "%f", &f2) == 1) {
        if((++c) == 3) {
          c = 0;
          if((cnt + 1) >= start) {
            if(icnt <= 1) {
              if(xref) {
                if(xref[a] >= 0)
                  fp = cs->coordPtr(xref[a]);
                else
                  fp = NULL;
              } else {
                fp = cs->coordPtr(a);
              }
              if(fp) {
                if(n_avg) {
                  *(fp++) += f0;
                  *(fp++) += f1;
                  *(fp++) += f2;
                } else {
                  *(fp++) = f0;
                  *(fp++) = f1;
                  *(fp++) = f2;
                }
              }
            }
          }
          if((++a) == I->NAtom) {

            cnt++;
            a = 0;
            if(b)
              p = nextline(p);
            b = 0;

            if(cs->PeriodicBoxType != CoordSet::NoPeriodicity) {
              /* read periodic box */

              c = 0;
              periodic = true;
              angles = true;

              p = ncopy(cc, p, 8);
              if(sscanf(cc, "%f", &box[0]) != 1)
                periodic = false;
              p = ncopy(cc, p, 8);
              if(sscanf(cc, "%f", &box[1]) != 1)
                periodic = false;
              p = ncopy(cc, p, 8);
              if(sscanf(cc, "%f", &box[2]) != 1)
                periodic = false;

              p = ncopy(cc, p, 8);
              if(sscanf(cc, "%f", &angle[0]) != 1)
                angles = false;

              p = ncopy(cc, p, 8);
              if(sscanf(cc, "%f", &angle[1]) != 1)
                angles = false;

              p = ncopy(cc, p, 8);
              if(sscanf(cc, "%f", &angle[2]) != 1)
                angles = false;
              if(periodic) {
                cs->Symmetry = pymol::make_unique<CSymmetry>(G);
                cs->Symmetry->Crystal.setDims(box);
                if(angles) {
                  cs->Symmetry->Crystal.setAngles(angle);
                }
                p = nextline(p);
                b = 0;
              }

              if(cs->PeriodicBoxType == CoordSet::Octahedral)
                periodic = false;       /* can't handle this yet... */
            }

            if((stop > 0) && (cnt >= stop))
              break;
            if(cnt >= start) {
              icnt--;
              if(icnt > 0) {
                PRINTFB(G, FB_ObjectMolecule, FB_Details)
                  " ObjectMolecule: skipping set %d...\n", cnt ENDFB(G);
              } else {
                icnt = interval;
                n_avg++;
              }

              if(icnt == interval) {
                if(n_avg < average) {
                  PRINTFB(G, FB_ObjectMolecule, FB_Details)
                    " ObjectMolecule: averaging set %d...\n", cnt ENDFB(G);
                } else {

                  /* compute average */
                  if(n_avg > 1) {
                    fp = cs->Coord.data();
                    for(i = 0; i < cs->NIndex; i++) {
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                      *(fp++) /= n_avg;
                    }
                  }
                  if(periodic && image) {       
                    /* Perform residue-based period image transformation */
                    i = 0;
                    r_cnt = 0;
                    r_act = 0;  /* 0 unspec, 1=load, 2=image, 3=leave */
                    r_val = -1;
                    while(r_act != 3) {
                      if(i >= cs->NIndex) {
                        if(r_cnt)
                          r_act = 2;
                        else
                          r_act = 3;
                      }
                      if(r_act == 0) {
                        /* start new residue */
                        r_cnt = 0;
                        r_act = 1;      /* now load */
                      }
                      if(r_act == 1) {
                        if(i < cs->NIndex) {

                          /* is there a coordinate for atom? */
                          if(xref) {
                            if(xref[i] >= 0)
                              fp = cs->coordPtr(xref[i]);
                            else
                              fp = NULL;
                          } else {
                            fp = cs->coordPtr(i);
                          }
                          if(fp) {      /* yes there is... */
                            if(r_cnt) {
                              if(r_val != I->AtomInfo[cs->IdxToAtm[i]].resv) {
                                r_act = 2;      /* end of residue-> time to image */
                              } else {
                                r_cnt++;
                                r_cent[0] += *(fp++);
                                r_cent[1] += *(fp++);
                                r_cent[2] += *(fp++);
                                r_fp_stop = fp; /* stop here */
                                i++;
                              }
                            } else {
                              r_val = I->AtomInfo[cs->IdxToAtm[i]].resv;
                              r_cnt++;
                              r_fp_start = fp;  /* start here */
                              r_cent[0] = *(fp++);
                              r_cent[1] = *(fp++);
                              r_cent[2] = *(fp++);
                              r_fp_stop = fp;   /* stop here */
                              i++;
                            }
                          } else {
                            i++;
                          }
                        } else {
                          r_act = 2;    /* image */
                        }
                      }

                      if(r_act == 2) {  /* time to image */
                        if(r_cnt) {
                          r_cent[0] /= r_cnt;
                          r_cent[1] /= r_cnt;
                          r_cent[2] /= r_cnt;
                          const auto& periodicbox = cs->getSymmetry()->Crystal;
                          transform33f3f(periodicbox.realToFrac(), r_cent, r_cent);
                          r_trans[0] = fmodf(1000.0F + shift[0] + r_cent[0], 1.0F);
                          r_trans[1] = fmodf(1000.0F + shift[1] + r_cent[1], 1.0F);
                          r_trans[2] = fmodf(1000.0F + shift[2] + r_cent[2], 1.0F);
                          r_trans[0] -= r_cent[0];
                          r_trans[1] -= r_cent[1];
                          r_trans[2] -= r_cent[2];
                          transform33f3f(periodicbox.fracToReal(), r_trans, r_trans);
                          fp = r_fp_start;
                          while(fp < r_fp_stop) {
                            *(fp++) += r_trans[0];
                            *(fp++) += r_trans[1];
                            *(fp++) += r_trans[2];
                          }
                        }
                        r_act = 0;      /* reset */
                        r_cnt = 0;
                      }
                    }
                  }

                  /* add new coord set */
                  cs->invalidateRep(cRepAll, cRepInvRep);
                  if(frame < 0)
                    frame = I->NCSet;
                  if(!I->NCSet) {
                    zoom_flag = true;
                  }

                  VLACheck(I->CSet, CoordSet *, frame);
                  if(I->NCSet <= frame)
                    I->NCSet = frame + 1;
                  delete I->CSet[frame];
                  I->CSet[frame] = cs;
                  ncnt++;

                  if(average < 2) {
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                      " ObjectMolecule: read set %d into state %d...\n", cnt, frame + 1
                      ENDFB(G);
                  } else {
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                      " ObjectMolecule: averaging set %d...\n", cnt ENDFB(G);
                    PRINTFB(G, FB_ObjectMolecule, FB_Details)
                      " ObjectMolecule: average loaded into state %d...\n", frame + 1
                      ENDFB(G);
                  }
                  frame++;
                  cs = CoordSetCopy(cs);
                  n_avg = 0;
                  if((stop > 0) && (cnt >= stop))
                    break;
                  if((max > 0) && (ncnt >= max))
                    break;
                }
              }
            } else {
              PRINTFB(G, FB_ObjectMolecule, FB_Details)
                " ObjectMolecule: skipping set %d...\n", cnt ENDFB(G);
            }
          }
        }
      } else {
        PRINTFB(G, FB_ObjectMolecule, FB_Errors)
          " ObjMolLoadTRJFile-Error: Failed to read expected coordinate value.\n  This traj. does not match the loaded parameter/topology file.\n  Likely cause: either the atom count or the periodic box settings\n  are inconsistent between the two files.\n"
          ENDFB(G);
        break;
      }
    }
    mfree(buffer);
  }
  delete cs;
  SceneChanged(G);
  SceneCountFrames(G);
  if(zoom_flag)
    if(SettingGetGlobal_i(G, cSetting_auto_zoom)) {
      ExecutiveWindowZoom(G, I->Name, 0.0, -1, 0, 0, quiet);        /* auto zoom (all states) */
    }

  return (I);
}

ObjectMolecule *ObjectMoleculeLoadRSTFile(PyMOLGlobals * G, ObjectMolecule * I,
                                          const char *fname, int frame, int quiet, char mode)
{
  /*
   * mode = 0: AMBER coordinate/restart file (one frame only)
   * mode = 1: AMBER trajectory
   * mode = 2: AMBER trajectory with box
   *           http://ambermd.org/formats.html
   */
  int ok = true;
  char *buffer, *p;
  char cc[MAXLINELEN];
  float f0, f1, f2, *fp;
  int a, b, c;
  int zoom_flag = false;
  CoordSet *cs = NULL;
  char ncolumn = 6; // number of coordinates per line
  char nbyte = 12;  // width of one coordinate

  if(mode > 0) {
    ncolumn = 10;
    nbyte = 8;
  }

#define BUFSIZE 4194304
#define GETTING_LOW 10000

  else {
    if(I->CSTmpl) {
      cs = CoordSetCopy(I->CSTmpl);
    } else if (I->NCSet > 0) {
      cs = CoordSetCopy(I->CSet[0]);
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " ObjMolLoadRSTFile: Missing topology" ENDFB(G);
      return (I);
    }
    CHECKOK(ok, cs);
    if (ok){
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
	" ObjMolLoadRSTFile: Loading from \"%s\".\n", fname ENDFB(G);
      p = buffer = FileGetContents(fname, NULL);
      if(!buffer)
        ok = ErrMessage(G, __func__, "Unable to open file!");
    }
    if (ok){
      p = nextline(p);
      if (mode == 0) // skip NATOM,TIME
        p = nextline(p);
    }
    a = 0;
    b = 0;
    c = 0;
    f1 = 0.0;
    f2 = 0.0;
    while(ok && *p) {
      p = ncopy(cc, p, nbyte);
      if((++b) == ncolumn) {
        b = 0;
        p = nextline(p);
      }
      f0 = f1;
      f1 = f2;
      if(sscanf(cc, "%f", &f2) == 1) {
        if((++c) == 3) {
          c = 0;
          fp = cs->coordPtr(a);
          *(fp++) = f0;
          *(fp++) = f1;
          *(fp++) = f2;

          if((++a) == I->NAtom) {
            a = 0;
            if(b)
              p = nextline(p);
            if(mode == 2) // skip box
              p = nextline(p);
            b = 0;
            /* add new coord set */
            cs->invalidateRep(cRepAll, cRepInvRep);
            if(frame < 0)
              frame = I->NCSet;
            if(!I->NCSet) {
              zoom_flag = true;
            }

            VLACheck(I->CSet, CoordSet *, frame);
	    CHECKOK(ok, I->CSet);
	    if (ok){
	      if(I->NCSet <= frame)
		I->NCSet = frame + 1;
	      delete I->CSet[frame];
	      I->CSet[frame] = cs;
	    }
            PRINTFB(G, FB_ObjectMolecule, FB_Details)
              " ObjectMolecule: read coordinates into state %d...\n", frame + 1 ENDFB(G);

            if (ok)
              cs = CoordSetCopy(cs);
            CHECKOK(ok, cs);

            if (mode == 0) // restart file has only one frame
              break;

            frame += 1;
          }
        }
      } else {
        PRINTFB(G, FB_ObjectMolecule, FB_Errors)
          " ObjMolLoadRSTFile: atom/coordinate mismatch.\n" ENDFB(G);
        break;
      }
    }
    mfree(buffer);
  }
  delete cs;

  SceneChanged(G);
  SceneCountFrames(G);
  if(zoom_flag){
    if(SettingGetGlobal_i(G, cSetting_auto_zoom)) {
      ExecutiveWindowZoom(G, I->Name, 0.0, -1, 0, 0, quiet);        /* auto zoom (all states) */
    }
  }
  return (I);
}

static const char *findflag(PyMOLGlobals * G, const char *p, const char *flag, const char *format)
{

  char cc[MAXLINELEN];
  char pat[MAXLINELEN] = "%FLAG ";
  int l;

  PRINTFD(G, FB_ObjectMolecule)
    " findflag: flag %s format %s\n", flag, format ENDFD;

  strcat(pat, flag);
  l = strlen(pat);
  while(*p) {
    p = ncopy(cc, p, l);
    if(WordMatch(G, cc, pat, true) < 0) {
      p = nextline(p);
      break;
    }
    p = nextline(p);
    if(!*p) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " ObjectMolecule-Error: Unrecognized file format (can't find \"%s\").\n", pat
        ENDFB(G);
    }
  }

  strcpy(pat, "%FORMAT(");
  strcat(pat, format);
  strcat(pat, ")");
  l = strlen(pat);
  while(*p) {
    p = ncopy(cc, p, l);
    if(WordMatch(G, cc, pat, true) < 0) {
      p = nextline(p);
      break;
    }
    p = nextline(p);
    if(!*p) {
      PRINTFB(G, FB_ObjectMolecule, FB_Errors)
        " ObjectMolecule-Error: Unrecognized file format (can't find \"%s\").\n", pat
        ENDFB(G);
    }

  }
  return (p);
}


/*========================================================================*/
static CoordSet *ObjectMoleculeTOPStr2CoordSet(PyMOLGlobals * G, const char *buffer,
                                        AtomInfoType ** atInfoPtr)
{
  const char *p;
  int nAtom;
  int a, b, c, bi, last_i, at_i, aa, rc;
  float *coord = NULL;
  float *f;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL, *ai;
  BondType *bond = NULL, *bd;
  int nBond = 0;
  int auto_show = RepGetAutoShowMask(G);
  int amber7 = false;

  WordType title;
  ResName *resn;

  char cc[MAXLINELEN];
  int ok = true;
  int i0, i1, i2;

  /* trajectory parameters */

  int NTYPES, NBONH, MBONA, NTHETH, MTHETA;
  int NPHIH, MPHIA, NHPARM, NPARM, NNB, NRES;
  int NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA;
  int NATYP, NPHB, IFPERT, NBPER, NGPER, NDPER;
  int MBPER, MGPER, MDPER, IFBOX = 0, NMXRS, IFCAP;
  int NEXTRA, IPOL = 0;
  int wid, col;
  float BETA;
  float BOX1, BOX2, BOX3;

  cset = CoordSetNew(G);

  p = buffer;
  nAtom = 0;
  if(atInfoPtr)
    atInfo = *atInfoPtr;
  if(!atInfo)
    ErrFatal(G, "TOPStr2CoordSet", "need atom information record!");
  /* failsafe for old version.. */

  ncopy(cc, p, 8);
  if(strcmp(cc, "%VERSION") == 0) {
    amber7 = true;
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
      " ObjectMolecule: Attempting to read Amber7 topology file.\n" ENDFB(G);
  } else {
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
      " ObjectMolecule: Assuming this is an Amber6 topology file.\n" ENDFB(G);
  }

  /* read title */
  if(amber7) {
    p = findflag(G, buffer, "TITLE", "20a4");
  }

  p = ncopy(cc, p, 20);
  title[0] = 0;
  sscanf(cc, "%s", title);
  p = nextline(p);

  if(amber7) {

    p = findflag(G, buffer, "POINTERS", "10I8");

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &nAtom) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NTYPES) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NBONH) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MBONA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NTHETH) == 1);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MTHETA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPHIH) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MPHIA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NHPARM) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPARM) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NNB) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NRES) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NBONA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NTHETA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPHIA) == 1);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NUMBND) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NUMANG) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPTRA) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NATYP) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NPHB) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &IFPERT) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NBPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NGPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NDPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MBPER) == 1);

    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MGPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &MDPER) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &IFBOX) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NMXRS) == 1);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &IFCAP) == 1);

    p = nextline(p);
    p = ncopy(cc, p, 8);
    ok = ok && (sscanf(cc, "%d", &NEXTRA) == 1);

  } else {

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &nAtom) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NTYPES) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NBONH) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MBONA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NTHETH) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MTHETA) == 1);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPHIH) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MPHIA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NHPARM) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPARM) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NNB) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NRES) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NBONA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NTHETA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPHIA) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NUMBND) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NUMANG) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPTRA) == 1);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NATYP) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NPHB) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &IFPERT) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NBPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NGPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NDPER) == 1);

    p = nextline(p);

    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MBPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MGPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &MDPER) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &IFBOX) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &NMXRS) == 1);
    p = ncopy(cc, p, 6);
    ok = ok && (sscanf(cc, "%d", &IFCAP) == 1);

    p = ncopy(cc, p, 6);
    if(sscanf(cc, "%d", &NEXTRA) != 1)
      NEXTRA = 0;

  }

  if(!ok) {
    ErrMessage(G, "TOPStrToCoordSet", "Error reading counts lines");
    ok_raise(1);
  } else {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " TOPStr2CoordSet: read counts line nAtom %d NBONA %d NBONH %d\n",
      nAtom, NBONA, NBONH ENDFB(G);
  }

  switch (IFBOX) {
  case 2:
    cset->PeriodicBoxType = CoordSet::Octahedral;
    PRINTFB(G, FB_ObjectMolecule, FB_Details)
      " TOPStrToCoordSet: Warning: can't currently image a truncated octahedron...\n"
      ENDFB(G);
    break;
  case 1:
    cset->PeriodicBoxType = CoordSet::Orthogonal;
    break;
  case 0:
  default:
    cset->PeriodicBoxType = CoordSet::NoPeriodicity;
    break;
  }

  p = nextline(p);

  if(ok) {
    VLACheck(atInfo, AtomInfoType, nAtom);

    if(amber7) {
      p = findflag(G, buffer, "ATOM_NAME", "20a4");
    }
    /* read atoms */

    b = 0;
    for(a = 0; a < nAtom; a++) {
      p = ntrim(cc, p, 4);
      atInfo[a].name = LexIdx(G, cc);
      if((++b) == 20) {
        b = 0;
        p = nextline(p);
      }
    }

    if(b)
      p = nextline(p);

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading atom names");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read atom names.\n" ENDFB(G);
    }

    /* read charges */

    if(amber7) {
      p = findflag(G, buffer, "CHARGE", "5E16.8");
    }

    b = 0;
    for(a = 0; a < nAtom; a++) {
      p = ncopy(cc, p, 16);
      ai = atInfo + a;
      if(!sscanf(cc, "%f", &ai->partialCharge))
        ok = false;
      else {
        ai->partialCharge /= 18.2223F;  /* convert to electron charge */
      }
      if((++b) == 5) {
        b = 0;
        p = nextline(p);
      }
    }

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading charges");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read charges.\n" ENDFB(G);
    }
    if(b)
      p = nextline(p);

    if(!amber7) {
      /* skip masses */

      p = skip_fortran(nAtom, 5, p);
    }

    /* read LJ atom types */

    if(amber7) {
      p = findflag(G, buffer, "ATOM_TYPE_INDEX", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    b = 0;
    for(a = 0; a < nAtom; a++) {
      p = ncopy(cc, p, wid);
      ai = atInfo + a;
      if(!sscanf(cc, "%d", &ai->customType))
        ok = false;
      if((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if(b)
      p = nextline(p);

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error LJ atom types");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read LJ atom types.\n" ENDFB(G);
    }

    if(!amber7) {
      /* skip excluded atom counts */

      p = skip_fortran(nAtom, 12, p);

      /* skip NB param arrays */

      p = skip_fortran(NTYPES * NTYPES, 12, p);
    }

    /* read residue labels */

    if(amber7) {
      p = findflag(G, buffer, "RESIDUE_LABEL", "20a4");
    }

    resn = pymol::malloc<ResName>(NRES);

    b = 0;
    for(a = 0; a < NRES; a++) {
      p = ncopy(cc, p, 4);
      if(!sscanf(cc, "%s", resn[a]))
        resn[a][0] = 0;
      if((++b) == 20) {
        b = 0;
        p = nextline(p);
      }
    }
    if(b)
      p = nextline(p);

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading residue labels");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read residue labels.\n" ENDFB(G);
    }

    /* read residue assignments */

    if(amber7) {
      p = findflag(G, buffer, "RESIDUE_POINTER", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    b = 0;
    last_i = 0;
    rc = 0;
    for(a = 0; a < NRES; a++) {
      p = ncopy(cc, p, wid);
      if(sscanf(cc, "%d", &at_i)) {
        if(last_i)
          for(aa = (last_i - 1); aa < (at_i - 1); aa++) {
            ai = atInfo + aa;
            ai->resn = LexIdx(G, resn[a - 1]);
            ai->resv = rc;
          }
        rc++;
        last_i = at_i;
      }
      if((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if(b)
      p = nextline(p);
    if(last_i)
      for(aa = (last_i - 1); aa < nAtom; aa++) {
        ai = atInfo + aa;
        ai->resn = LexIdx(G, resn[NRES - 1]);
        ai->resv = rc;
      }
    rc++;

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading residues");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read residues.\n" ENDFB(G);
    }

    FreeP(resn);

    if(!amber7) {
      /* skip bond force constants */

      p = skip_fortran(NUMBND, 5, p);

      /* skip bond lengths */

      p = skip_fortran(NUMBND, 5, p);

      /* skip angle force constant */

      p = skip_fortran(NUMANG, 5, p);

      /* skip angle eq */

      p = skip_fortran(NUMANG, 5, p);

      /* skip dihedral force constant */

      p = skip_fortran(NPTRA, 5, p);

      /* skip dihedral periodicity */

      p = skip_fortran(NPTRA, 5, p);

      /* skip dihedral phases */

      p = skip_fortran(NPTRA, 5, p);

      /* skip SOLTYs */

      p = skip_fortran(NATYP, 5, p);

      /* skip LJ terms r12 */

      p = skip_fortran((NTYPES * (NTYPES + 1)) / 2, 5, p);

      /* skip LJ terms r6 */

      p = skip_fortran((NTYPES * (NTYPES + 1)) / 2, 5, p);

    }

    /* read bonds */

    if(amber7) {
      p = findflag(G, buffer, "BONDS_INC_HYDROGEN", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    nBond = NBONH + NBONA;

    bond = VLACalloc(BondType, nBond);

    bi = 0;

    b = 0;
    c = 0;
    i0 = 0;
    i1 = 0;
    for(a = 0; a < 3 * NBONH; a++) {
      p = ncopy(cc, p, wid);
      i2 = i1;
      i1 = i0;
      if(!sscanf(cc, "%d", &i0))
        ok = false;
      if((++c) == 3) {
        c = 0;
        bd = bond + bi;
        bd->index[0] = (abs(i2) / 3);
        bd->index[1] = (abs(i1) / 3);
        bd->order = 1;
        bi++;
      }
      if((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if(b)
      p = nextline(p);

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error hydrogen containing bonds");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read %d hydrogen containing bonds.\n", NBONH ENDFB(G);
    }

    if(amber7) {
      p = findflag(G, buffer, "BONDS_WITHOUT_HYDROGEN", "10I8");
      col = 10;
      wid = 8;
    } else {
      col = 12;
      wid = 6;
    }

    b = 0;
    c = 0;
    for(a = 0; a < 3 * NBONA; a++) {
      p = ncopy(cc, p, wid);
      i2 = i1;
      i1 = i0;
      if(!sscanf(cc, "%d", &i0))
        ok = false;
      if((++c) == 3) {
        c = 0;
        bd = bond + bi;
        bd->index[0] = (abs(i2) / 3);
        bd->index[1] = (abs(i1) / 3);
        bd->order = 1; // PYMOL-2707
        bi++;
      }
      if((++b) == col) {
        b = 0;
        p = nextline(p);
      }
    }
    if(b)
      p = nextline(p);

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error hydrogen free bonds");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read %d hydrogen free bonds.\n", NBONA ENDFB(G);
    }

    if(!amber7) {

      /* skip hydrogen angles */

      p = skip_fortran(4 * NTHETH, 12, p);

      /* skip non-hydrogen angles */

      p = skip_fortran(4 * NTHETA, 12, p);

      /* skip hydrogen dihedrals */

      p = skip_fortran(5 * NPHIH, 12, p);

      /* skip non hydrogen dihedrals */

      p = skip_fortran(5 * NPHIA, 12, p);

      /* skip nonbonded exclusions */

      p = skip_fortran(NNB, 12, p);

      /* skip hydrogen bonds ASOL */

      p = skip_fortran(NPHB, 5, p);

      /* skip hydrogen bonds BSOL */

      p = skip_fortran(NPHB, 5, p);

      /* skip HBCUT */

      p = skip_fortran(NPHB, 5, p);

    }
    /* read AMBER atom types */

    if(amber7) {
      p = findflag(G, buffer, "AMBER_ATOM_TYPE", "20a4");
    }

    b = 0;
    for(a = 0; a < nAtom; a++) {
      OrthoLineType temp;
      p = ncopy(cc, p, 4);
      ai = atInfo + a;
      if(sscanf(cc, "%s", temp) != 1)
        ok = false;
      else {
        ai->textType = LexIdx(G, temp);
      }
      if((++b) == 20) {
        b = 0;
        p = nextline(p);
      }
    }
    if(b)
      p = nextline(p);

    if(!ok) {
      ErrMessage(G, "TOPStrToCoordSet", "Error reading atom types");
    } else {
      PRINTFB(G, FB_ObjectMolecule, FB_Blather)
        " TOPStr2CoordSet: read atom types.\n" ENDFB(G);
    }

    if(!amber7) {
      /* skip TREE classification */

      p = skip_fortran(nAtom, 20, p);

      /* skip tree joining information */

      p = skip_fortran(nAtom, 12, p);

      /* skip last atom rotated blah blah blah */

      p = skip_fortran(nAtom, 12, p);

    }

    if(IFBOX > 0) {

      int IPTRES, NSPM = 0, NSPSOL;

      if(amber7) {
        p = findflag(G, buffer, "SOLVENT_POINTERS", "3I8");
        wid = 8;
      } else {
        wid = 6;
      }
      p = ncopy(cc, p, wid);
      ok = ok && (sscanf(cc, "%d", &IPTRES) == 1);
      p = ncopy(cc, p, wid);
      ok = ok && (sscanf(cc, "%d", &NSPM) == 1);
      p = ncopy(cc, p, wid);
      ok = ok && (sscanf(cc, "%d", &NSPSOL) == 1);

      ok_assert(1, ok);

      p = nextline(p);

      if(amber7) {
        p = findflag(G, buffer, "ATOMS_PER_MOLECULE", "10I8");
        col = 10;
      } else {
        col = 12;
      }

      /* skip num atoms per box */

      p = skip_fortran(NSPM, col, p);

      if(amber7) {
        p = findflag(G, buffer, "BOX_DIMENSIONS", "5E16.8");
      }
      wid = 16;

      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BETA) == 1);
      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BOX1) == 1);
      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BOX2) == 1);
      p = ncopy(cc, p, 16);
      ok = ok && (sscanf(cc, "%f", &BOX3) == 1);

      if(ok) {
        float angle[3] = {90.f, BETA, 90.f};

        if((BETA > 109.47) && (BETA < 109.48)) {
          cset->PeriodicBoxType = CoordSet::Octahedral;
          angle[0] = 109.47122;
          angle[1] = 109.47122;
          angle[2] = 109.47122;
        } else if(BETA == 60.0) {
          angle[0] = 60.0;   /* rhombic dodecahedron (from ptraj.c) */
          angle[1] = 90.0;
          angle[2] = 60.0;
        }

        cset->Symmetry = pymol::make_unique<CSymmetry>(G);
        cset->Symmetry->Crystal.setDims(BOX1, BOX2, BOX3);
        cset->Symmetry->Crystal.setAngles(angle);
      }
      /* skip periodic box */

      p = nextline(p);

    }

    if(!amber7) {

      if(IFCAP > 0) {
        p = nextline(p);
        p = nextline(p);
        p = nextline(p);
      }

      if(IFPERT > 0) {

        /* skip perturbed bond atoms */

        p = skip_fortran(2 * NBPER, 12, p);

        /* skip perturbed bond atom pointers */

        p = skip_fortran(2 * NBPER, 12, p);

        /* skip perturbed angles */

        p = skip_fortran(3 * NGPER, 12, p);

        /* skip perturbed angle pointers */

        p = skip_fortran(2 * NGPER, 12, p);

        /* skip perturbed dihedrals */

        p = skip_fortran(4 * NDPER, 12, p);

        /* skip perturbed dihedral pointers */

        p = skip_fortran(2 * NDPER, 12, p);

        /* skip residue names */

        p = skip_fortran(NRES, 20, p);

        /* skip atom names */

        p = skip_fortran(nAtom, 20, p);

        /* skip atom symbols */

        p = skip_fortran(nAtom, 20, p);

        /* skip unused field */

        p = skip_fortran(nAtom, 5, p);

        /* skip perturbed flags */

        p = skip_fortran(nAtom, 12, p);

        /* skip LJ atom flags */

        p = skip_fortran(nAtom, 12, p);

        /* skip perturbed charges */

        p = skip_fortran(nAtom, 5, p);

      }

      if(IPOL > 0) {

        /* skip atomic polarizabilities */

        p = skip_fortran(nAtom, 5, p);

      }

      if((IPOL > 0) && (IFPERT > 0)) {

        /* skip atomic polarizabilities */

        p = skip_fortran(nAtom, 5, p);

      }
    }
    /* for future reference 

       %FLAG LES_NTYP
       %FORMAT(10I8)
       %FLAG LES_TYPE
       %FORMAT(10I8)
       %FLAG LES_FAC
       %FORMAT(5E16.8)
       %FLAG LES_CNUM
       %FORMAT(10I8)
       %FLAG LES_ID
       %FORMAT(10I8)

       Here is the additional information for LES topology formats:
       First, if NPARM ==1, LES entries are in topology (NPARM is the 10th
       entry in the initial list of control parameters); otherwise the standard
       format applies.
       So, with NPARM=1, you just need to read a few more things at the very
       end of topology file:
       LES_NTYP (format: I6) ... one number, number of LES types
       and four arrays:
       LES_TYPE (12I6) ... NATOM integer entries
       LES_FAC (E16.8) ... LES_NTYPxLES_NTYP float entries
       LES_CNUM (12I6) ... NATOM integer entries
       LES_ID (12I6)   ... NATOM integer entries

       and that's it. Your parser must have skipped this information because it
       was at the end of the file. Maybe that's good enough.

     */

    coord = VLAlloc(float, 3 * nAtom);

    f = coord;
    for(a = 0; a < nAtom; a++) {
      *(f++) = 0.0;
      *(f++) = 0.0;
      *(f++) = 0.0;
      ai = atInfo + a;
      ai->id = a + 1;           /* assign 1-based identifiers */
      ai->rank = a;
      AtomInfoAssignParameters(G, ai);
      AtomInfoAssignColors(G, ai);
      ai->visRep = auto_show;
    }
  }
  if(ok) {
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->TmpBond = pymol::vla_take_ownership(bond);
    cset->NTmpBond = nBond;
  } else {
ok_except1:
    delete cset;
    cset = NULL;
    ErrMessage(G, __func__, "failed");
  }
  if(atInfoPtr)
    *atInfoPtr = atInfo;

  return (cset);
}


/*========================================================================*/
static ObjectMolecule *ObjectMoleculeReadTOPStr(PyMOLGlobals * G, ObjectMolecule * I,
                                         char *TOPStr, int frame, int discrete)
{
  CoordSet *cset = NULL;
  pymol::vla<AtomInfoType> atInfo(1);
  int ok = true;
  int isNew = true;
  unsigned int nAtom = 0;

  if(!I)
    isNew = true;
  else
    isNew = false;

  if(ok) {
    if(isNew) {
      I = (ObjectMolecule *) new ObjectMolecule(G, discrete);
      CHECKOK(ok, I);
      if (ok)
        std::swap(atInfo, I->AtomInfo);
    }
    if(ok && isNew) {
      I->Color = AtomInfoUpdateAutoColor(G);
    }

    if (ok)
      cset = ObjectMoleculeTOPStr2CoordSet(G, TOPStr, &atInfo);
    CHECKOK(ok, cset);
  }

  /* include coordinate set */
  if(ok) {
    nAtom = cset->NIndex;

    if(I->DiscreteFlag && atInfo) {
      unsigned int a;
      int fp1 = frame + 1;
      AtomInfoType *ai = atInfo.data();
      for(a = 0; a < nAtom; a++) {
        (ai++)->discrete_state = fp1;
      }
    }

    cset->Obj = I;
    cset->enumIndices();
    cset->invalidateRep(cRepAll, cRepInvRep);
    if(isNew) {
      std::swap(I->AtomInfo, atInfo);
    } else if (ok){
      ok &= ObjectMoleculeMerge(I, std::move(atInfo), cset, false, cAIC_AllMask, true);
    }
    if(isNew)
      I->NAtom = nAtom;
    /* 
       if(frame<0) frame=I->NCSet;
       VLACheck(I->CSet,CoordSet*,frame);
       if(I->NCSet<=frame) I->NCSet=frame+1;
       if(I->CSet[frame]) I->CSet[frame]->fFree(I->CSet[frame]);
       I->CSet[frame] = cset;
     */

    if(ok && isNew)
      ok &= ObjectMoleculeConnect(I, cset, false);
    if(cset->Symmetry && (!I->Symmetry)) {
      I->Symmetry.reset(new CSymmetry(*cset->Symmetry));
      CHECKOK(ok, I->Symmetry);
    }

    delete I->CSTmpl;
    I->CSTmpl = cset;           /* save template coordinate set */

    SceneCountFrames(G);
    if (ok)
      ok &= ObjectMoleculeExtendIndices(I, -1);
    if (ok)
      ok &= ObjectMoleculeSort(I);
    if (ok){
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
    }
  }
  if (!ok){
    DeleteP(I)
  }
  return (I);
}

ObjectMolecule *ObjectMoleculeLoadTOPFile(PyMOLGlobals * G, ObjectMolecule * obj,
                                          const char *fname, int frame, int discrete)
{
  ObjectMolecule *I = NULL;
  char *buffer;

  buffer = FileGetContents(fname, NULL);

  if(!buffer)
    ErrMessage(G, __func__, "Unable to open file!");
  else {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " %s: Loading from %s.\n", __func__, fname ENDFB(G);

    I = ObjectMoleculeReadTOPStr(G, obj, buffer, frame, discrete);
    mfree(buffer);
  }

  return (I);
}

void ObjectMoleculeSculptClear(ObjectMolecule * I)
{
  PRINTFD(I->G, FB_ObjectMolecule)
    " %s: entered.\n", __func__ ENDFD;
  if(I->Sculpt)
    DeleteP(I->Sculpt);
}

void ObjectMoleculeSculptImprint(ObjectMolecule * I, int state, int match_state,
                                 int match_by_segment)
{
  PRINTFD(I->G, FB_ObjectMolecule)
    " %s: entered.\n", __func__ ENDFD;

  if(!I->Sculpt)
    I->Sculpt = new CSculpt(I->G);
  SculptMeasureObject(I->Sculpt, I, state, match_state, match_by_segment);
}

float ObjectMoleculeSculptIterate(ObjectMolecule * I, int state, int n_cycle,
                                  float *center)
{
  PRINTFD(I->G, FB_ObjectMolecule)
    " %s: entered.\n", __func__ ENDFD;
  if(I->Sculpt) {
    return SculptIterateObject(I->Sculpt, I, state, n_cycle, center);
  } else
    return 0.0F;
}

/**
 * - Assigns new id to all atoms with AtomInfoType.id == -1
 * - Assigns ObjectMolecule.AtomCounter if -1
 *
 * Cost: O(NAtom + NBond)
 */
void ObjectMoleculeUpdateIDNumbers(ObjectMolecule * I)
{
  int a;
  int max;
  AtomInfoType *ai;

  if(I->AtomCounter < 0) {
    max = -1;
    ai = I->AtomInfo.data();
    for(a = 0; a < I->NAtom; a++) {
      if(ai->id > max)
        max = ai->id;
      ai++;
    }
    I->AtomCounter = max + 1;
  }
  ai = I->AtomInfo.data();
  for(a = 0; a < I->NAtom; a++) {
    if(ai->id < 0)
      ai->id = I->AtomCounter++;
    ai++;
  }
}

/*========================================================================*/
int ObjectMoleculeGetPhiPsi(ObjectMolecule * I, int ca, float *phi, float *psi, int state)
{
  int np = -1;
  int cm = -1;
  int c = -1;
  int n = -1;
  int result = false;
  const AtomInfoType *ai;
  float v_ca[3];
  float v_n[3];
  float v_c[3];
  float v_cm[3];
  float v_np[3];
  auto G = I->G;

  ai = I->AtomInfo;

  if(ai[ca].name == G->lex_const.CA) {
    /* find C */
    for (auto const& neighbor : AtomNeighbors(I, ca)) {
      if (ai[neighbor.atm].name == G->lex_const.C) {
        c = neighbor.atm;
        break;
      }
    }

    /* find N */
    for (auto const& neighbor : AtomNeighbors(I, ca)) {
      if (ai[neighbor.atm].name == G->lex_const.N) {
        n = neighbor.atm;
        break;
      }
    }

    /* find NP */
    if(c >= 0) {
      for (auto const& neighbor : AtomNeighbors(I, c)) {
        if (ai[neighbor.atm].name == G->lex_const.N) {
          np = neighbor.atm;
          break;
        }
      }
    }

    /* find CM */
    if(n >= 0) {
      for (auto const& neighbor : AtomNeighbors(I, n)) {
        if (ai[neighbor.atm].name == G->lex_const.C) {
          cm = neighbor.atm;
          break;
        }
      }
    }
    if((ca >= 0) && (np >= 0) && (c >= 0) && (n >= 0) && (cm >= 0)) {
      if(ObjectMoleculeGetAtomVertex(I, state, ca, v_ca) &&
         ObjectMoleculeGetAtomVertex(I, state, n, v_n) &&
         ObjectMoleculeGetAtomVertex(I, state, c, v_c) &&
         ObjectMoleculeGetAtomVertex(I, state, cm, v_cm) &&
         ObjectMoleculeGetAtomVertex(I, state, np, v_np)) {

        (*phi) = rad_to_deg(get_dihedral3f(v_c, v_ca, v_n, v_cm));
        (*psi) = rad_to_deg(get_dihedral3f(v_np, v_c, v_ca, v_n));
        result = true;
      }
    }
  }
  return (result);
}


/*========================================================================*/
int ObjectMoleculeCheckBondSep(ObjectMolecule * I, int a0, int a1, int dist)
{
  int result = false;
  int n0;
  int stack[MAX_BOND_DIST + 1];
  int history[MAX_BOND_DIST + 1];
  int depth = 0;
  int distinct;
  int a;
  if(dist > MAX_BOND_DIST)
    return false;

  auto* const Neighbor = I->getNeighborArray();

  depth = 1;
  history[depth] = a0;
  stack[depth] = Neighbor[a0] + 1;   /* go to first neighbor */
  while(depth) {                /* keep going until we've traversed tree */
    while(Neighbor[stack[depth]] >= 0) {     /* end of branches? go back up one bond */
      n0 = Neighbor[stack[depth]];   /* get current neighbor index */
      stack[depth] += 2;        /* set up next neighbor */
      distinct = true;          /* check to see if current candidate is distinct from ancestors */
      for(a = 1; a < depth; a++) {
        if(history[a] == n0)
          distinct = false;
      }
      if(distinct) {
        if(depth < dist) {      /* are not yet at the proper distance? */
          if(distinct) {
            depth++;
            stack[depth] = Neighbor[n0] + 1; /* then keep moving outward */
            history[depth] = n0;
          }
        } else if(n0 == a1)     /* otherwise, see if we have a match */
          result = true;
      }
    }
    depth--;
  }
  return result;
}


/*========================================================================*/
void ObjectGotoState(pymol::CObject* I, int state)
{
  auto nstates = I->getNFrame();
  if (nstates > 1 || !SettingGet<bool>(I->G, cSetting_static_singletons)) {
    if(state > nstates)
      state = nstates - 1;
    if(state < 0)
      state = nstates - 1;
    SceneSetFrame(I->G, 0, state);
  }
}


/*========================================================================*/
CObjectState* ObjectMolecule::_getObjectState(int state)
{
  return CSet[state];
}


/*========================================================================*/
pymol::copyable_ptr<CSetting>* ObjectMolecule::getSettingHandle(int state)
{
  auto I = this;

  if (state < -1) {
    state = I->getCurrentState();
  }

  if(state < 0) {
    return (&I->Setting);
  } else if(state < I->NCSet) {
    if(I->CSet[state]) {
      return (&I->CSet[state]->Setting);
    } else {
      return (NULL);
    }
  } else {
    return (NULL);
  }
}

/*========================================================================*/
CSymmetry const* ObjectMolecule::getSymmetry(int state) const
{
  auto cs = getCoordSet(state);
  if (cs && cs->Symmetry) {
    return cs->Symmetry.get();
  }

  return Symmetry.get();
}

/**
 * @return False if state does not exist
 */
bool ObjectMolecule::setSymmetry(CSymmetry const& symmetry, int state)
{
  bool const all_states = (state == cSelectorUpdateTableAllStates);
  bool success = false;

  if (all_states) {
    Symmetry.reset(new CSymmetry(symmetry));
    success = true;
  }

  for (StateIterator iter(G, Setting.get(), state, NCSet); iter.next();) {
    auto cs = CSet[iter.state];
    if (cs) {
      cs->Symmetry.reset(all_states ? nullptr : new CSymmetry(symmetry));
      cs->UnitCellCGO.reset();
      cs->invalidateRep(cRepCell, cRepInvRep);
      success = true;
    }
  }

  return success;
}

/*========================================================================*/
pymol::Result<> ObjectMoleculeSetStateTitle(ObjectMolecule * I, int state, const char *text)
{
  auto cs = I->getCoordSet(state);
  if (!cs) {
    return pymol::make_error("Invalid state ", state + 1);
  }
  cs->setTitle(text);
  return {};
}


/*========================================================================*/
const char *ObjectMoleculeGetStateTitle(const ObjectMolecule * I, int state)
{
  auto cs = I->getCoordSet(state);
  if (!cs) {
    PRINTFB(I->G, FB_ObjectMolecule, FB_Errors)
      "Error: invalid state %d\n", state + 1 ENDFB(I->G);
    return nullptr;
  }
  return cs->Name;
}


/*========================================================================*/
void ObjectMoleculeRenderSele(ObjectMolecule * I, int curState, int sele, int vis_only SELINDICATORARG)
{

  PyMOLGlobals *G = I->G;
  CoordSet *cs;
  int a, nIndex;
  const float *coord, *v;
  int flag = true;
  int all_vis = !vis_only;
  int visRep;
  float tmp_matrix[16], v_tmp[3], *matrix = NULL;
  int use_matrices =
    SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);

  if(use_matrices<0) use_matrices = 0;

  if (SettingGetIfDefined_i(G, I->Setting.get(), cSetting_all_states, &a)) {
    curState = a ? -1 : SettingGet_i(G, I->Setting.get(), NULL, cSetting_state);
  } else if (SettingGetIfDefined_i(G, I->Setting.get(), cSetting_state, &a)) {
    curState = a - 1;
  }

  if(G->HaveGUI && G->ValidContext) {
    const AtomInfoType *atInfo = I->AtomInfo.data();

    for(StateIterator iter(G, I->Setting.get(), curState, I->NCSet);
        iter.next();) {
      if((cs = I->CSet[iter.state])) {
	    const auto* idx2atm = cs->IdxToAtm.data();
	    nIndex = cs->NIndex;
	    coord = cs->Coord;
	    if(use_matrices && !cs->Matrix.empty()) {
	      copy44d44f(cs->Matrix.data(), tmp_matrix);
	      matrix = tmp_matrix;
	    } else
	      matrix = NULL;
	    
	    if(I->TTTFlag) {
	      if(!matrix) {
		convertTTTfR44f(I->TTT, tmp_matrix);
	      } else {
		float ttt[16];
		convertTTTfR44f(I->TTT, ttt);
		left_multiply44f44f(ttt, tmp_matrix);
	      }
	      matrix = tmp_matrix;
	    }
	    
	    for(a = 0; a < nIndex; a++) {
	      if(SelectorIsMember(G, atInfo[*(idx2atm++)].selEntry, sele)) {
		if(all_vis)
		  flag = true;
		else {
		  visRep = atInfo[idx2atm[-1]].visRep;
		  flag = false;
		  if(visRep & (cRepCylBit | cRepSphereBit | cRepSurfaceBit |
			       cRepLabelBit | cRepNonbondedSphereBit | cRepCartoonBit |
			       cRepRibbonBit | cRepLineBit | cRepMeshBit |
			       cRepDotBit | cRepNonbondedBit)){
		    flag = true;
		  }
		}
		if(flag) {
		  v = coord + a + a + a;
		  if(matrix) {
		    transform44f3f(matrix, v, v_tmp);
		    if (SELINDICATORVAR)
		      CGOVertexv(SELINDICATORVAR, v_tmp);
		    else
		      glVertex3fv(v_tmp);
		  } else {
		    if (SELINDICATORVAR)
		      CGOVertexv(SELINDICATORVAR, v);
		    else
		      glVertex3fv(v);
		  }
		}
	      }
	    }
	  }
    }
  }
}

/*========================================================================*/
static CoordSet *ObjectMoleculeXYZStr2CoordSet(PyMOLGlobals * G, const char *buffer,
                                               AtomInfoType ** atInfoPtr, const char **restart)
{
  const char *p, *p_store;
  int nAtom;
  int a, c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL, *ai;
  char cc[MAXLINELEN];
  int atomCount;
  BondType *bond = NULL;
  int nBond = 0;
  int b1, b2;
  WordType tmp_name;
  int auto_show = RepGetAutoShowMask(G);
  int tinker_xyz = true;
  int valid_atom;
  int have_n_atom = false;
  BondType *ii;

  p = buffer;
  nAtom = 0;
  atInfo = *atInfoPtr;

  p_store = p;
  p = ncopy(cc, p, MAXLINELEN - 1);
  if(sscanf(cc, "%d", &nAtom) != 1) {
    nAtom = 0;
    tinker_xyz = false;
    p = p_store;
  } else {
    have_n_atom = true;
    p = nskip(p, 2);
    p = ncopy(tmp_name, p, sizeof(WordType) - 1);
    p = nextline(p);
  }

  if(tinker_xyz && nAtom) {     /* test Tinker XYZ formatting assumption */
    const char *pp = p;
    int dummy_int;
    float dummy_float;
    AtomName dummy_name;

    pp = ncopy(cc, pp, 6);
    if(!sscanf(cc, "%d", &dummy_int))
      tinker_xyz = false;       /* id */
    pp = nskip(pp, 2);
    pp = ncopy(cc, pp, 3);
    if(sscanf(cc, "%s", dummy_name) != 1)
      tinker_xyz = false;       /* name */
    pp = ncopy(cc, pp, 12);
    if(sscanf(cc, "%f", &dummy_float) != 1)
      tinker_xyz = false;       /* x */
    pp = ncopy(cc, pp, 12);
    if(sscanf(cc, "%f", &dummy_float) != 1)
      tinker_xyz = false;       /* y */
    pp = ncopy(cc, pp, 12);
    if(sscanf(cc, "%f", &dummy_float) != 1)
      tinker_xyz = false;       /* z */
    pp = ncopy(cc, pp, 6);
    if(sscanf(cc, "%d", &dummy_int) != 1)
      tinker_xyz = false;       /* numeric type */
  }

  if(!tinker_xyz) {
    const char *pp = p;
    int have_atom_line = true;
    float dummy_float;
    AtomName dummy_name;

    pp = ParseWordCopy(cc, pp, sizeof(AtomName) - 1);
    if(sscanf(cc, "%s", dummy_name) != 1)
      have_atom_line = false;   /* name */
    pp = ParseWordCopy(cc, pp, MAXLINELEN - 1);
    if(sscanf(cc, "%f", &dummy_float) != 1)
      have_atom_line = false;   /* x */
    pp = ParseWordCopy(cc, pp, MAXLINELEN - 1);
    if(sscanf(cc, "%f", &dummy_float) != 1)
      have_atom_line = false;   /* y */
    pp = ParseWordCopy(cc, pp, MAXLINELEN - 1);
    if(sscanf(cc, "%f", &dummy_float) != 1)
      have_atom_line = false;   /* z */
    if(!have_atom_line) {       /* copy the comment line into the title field */
      p = ncopy(tmp_name, p, sizeof(WordType) - 1);
      p = nextline(p);
    }
  }

  if(nAtom) {
    coord = VLAlloc(float, 3 * nAtom);
    if(atInfo)
      VLACheck(atInfo, AtomInfoType, nAtom);
  } else {
    coord = VLAlloc(float, 3);
  }

  if(tinker_xyz) {
    nBond = 0;
    bond = VLACalloc(BondType, 6 * nAtom);      /* is this a safe assumption? */
    ii = bond;
  }

  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " ObjectMoleculeReadXYZ: Found %i atoms...\n", nAtom ENDFB(G);

  a = 0;
  atomCount = 0;
  while(*p) {
    VLACheck(atInfo, AtomInfoType, atomCount);
    ai = atInfo + atomCount;

    if(!tinker_xyz) {
      valid_atom = true;

      p = ParseWordCopy(cc, p, MAXLINELEN - 1);
      UtilCleanStr(cc);
      if(!cc[0])
        valid_atom = false;
      if(valid_atom) {
        strncpy(ai->elem, cc, cElemNameLen);
        ai->name = LexIdx(G, cc);

        ai->rank = atomCount;
        ai->id = atomCount + 1;

        VLACheck(coord, float, a * 3 + 2);
        p = ParseWordCopy(cc, p, MAXLINELEN - 1);
        if(sscanf(cc, "%f", coord + a) != 1)
          valid_atom = false;
        p = ParseWordCopy(cc, p, MAXLINELEN - 1);
        if(sscanf(cc, "%f", coord + a + 1) != 1)
          valid_atom = false;
        p = ParseWordCopy(cc, p, MAXLINELEN - 1);
        if(sscanf(cc, "%f", coord + a + 2) != 1)
          valid_atom = false;

        ai->resn = LexIdx(G, "UNK");

        ai->alt[0] = 0;
        ai->chain = 0;
        ai->resv = atomCount + 1;

        ai->q = 1.0;
        ai->b = 0.0;

        ai->segi = 0;

        ai->visRep = auto_show;

        /* in the absense of external tinker information, assume hetatm */

        ai->hetatm = 1;

        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);
      }
    } else {                    /* tinker XYZ */

      valid_atom = true;

      p = ncopy(cc, p, 6);
      if(!sscanf(cc, "%d", &ai->id))
        break;
      ai->rank = atomCount;

      p = nskip(p, 2);          /* to 12 */
      p = ntrim(cc, p, 3);
      ai->name = LexIdx(G, cc);

      ai->resn = LexIdx(G, "UNK");

      ai->alt[0] = 0;
      ai->chain = 0;

      ai->resv = atomCount + 1;

      valid_atom = true;

      p = ncopy(cc, p, 12);
      sscanf(cc, "%f", coord + a);
      p = ncopy(cc, p, 12);
      sscanf(cc, "%f", coord + (a + 1));
      p = ncopy(cc, p, 12);
      sscanf(cc, "%f", coord + (a + 2));

      ai->q = 1.0;
      ai->b = 0.0;

      ai->segi = 0;
      ai->elem[0] = 0;          /* let atom info guess/infer atom type */

      ai->visRep = auto_show;

      p = ncopy(cc, p, 6);
      sscanf(cc, "%d", &ai->customType);

      /* in the absense of external tinker information, assume hetatm */

      ai->hetatm = 1;

      AtomInfoAssignParameters(G, ai);
      AtomInfoAssignColors(G, ai);

      b1 = atomCount;
      for(c = 0; c < 6; c++) {
        p = ncopy(cc, p, 6);
        if(!cc[0])
          break;
        if(!sscanf(cc, "%d", &b2))
          break;
        if(b1 < (b2 - 1)) {
          VLACheck(bond, BondType, nBond);
          ii = bond + nBond;
          nBond++;
          ii->index[0] = b1;
          ii->index[1] = b2 - 1;
          ii->order = 1;        /* missing bond order information */
          ii++;
        }
      }
    }

    if(valid_atom) {
      PRINTFD(G, FB_ObjectMolecule)
        " ObjectMolecule-DEBUG: %s %s %d %s %8.3f %8.3f %8.3f %6.2f %6.2f %s\n",
         LexStr(G, ai->name), LexStr(G, ai->resn), ai->resv, LexStr(G, ai->chain),
         *(coord + a), *(coord + a + 1), *(coord + a + 2), ai->b, ai->q, LexStr(G, ai->segi) ENDFD;

      a += 3;
      atomCount++;

    }
    p = nextline(p);
    if(have_n_atom && (atomCount >= nAtom)) {
      int dummy;
      ncopy(cc, p, MAXLINELEN - 1);
      if(sscanf(cc, "%d", &dummy) == 1)
        *restart = p;
      break;
    }
  }

  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " XYZStr2CoordSet: Read %d bonds.\n", nBond ENDFB(G);

  if(!tinker_xyz)
    nAtom = atomCount;          /* use number of atoms actually read */

  cset = CoordSetNew(G);
  cset->NIndex = nAtom;
  cset->Coord = pymol::vla_take_ownership(coord);
  cset->TmpBond = pymol::vla_take_ownership(bond);
  cset->NTmpBond = nBond;
  strcpy(cset->Name, tmp_name);
  if(atInfoPtr)
    *atInfoPtr = atInfo;
  return (cset);
}


/*========================================================================*/
int ObjectMoleculeAreAtomsBonded(ObjectMolecule * I, int i0, int i1)
{
  int result = false;
  int a;
  const BondType *b;
  b = I->Bond;
  for(a = 0; a < I->NBond; a++) {
    if(i0 == b->index[0]) {
      if(i1 == b->index[1]) {
        result = true;
        break;
      }
    }
    if(i1 == b->index[0]) {
      if(i0 == b->index[1]) {
        result = true;
        break;
      }
    }
    b++;
  }
  return (result);
}


/*========================================================================*/
int ObjectMoleculeRenameAtoms(ObjectMolecule * I, int *flag, int force)
{
  PyMOLGlobals * G = I->G;
  AtomInfoType *ai;
  int a;
  int result;
  if(force) {
    ai = I->AtomInfo.data();
    if(!flag) {
      for(a = 0; a < I->NAtom; a++) {
        LexAssign(G, ai->name, 0);
        ai++;
      }
    } else {
      for(a = 0; a < I->NAtom; a++) {
        if(flag[a])
          LexAssign(G, ai->name, 0);
        ai++;
      }
    }
  }
  result = AtomInfoUniquefyNames(I->G, NULL, 0, I->AtomInfo.data(), flag, I->NAtom);
  return result;
}


/*========================================================================*/
/**
 * Transform coordinates of `cs` in place and then append it to `tcs`.
 *
 * @param coord_orig Original cs coordinates
 * @param at0 Anchor atom index
 * @param mat1_inv Source coordinate system
 * @param valence_vec Valence vector. If NULL, don't move. If null-vector, don't merge.
 */
static bool AddCoordinateIntoCoordSet(CoordSet* tcs, CoordSet* cs,
    const float* coord_orig, int at0, const float* mat1_inv,
    const float* valence_vec)
{
  if (!tcs) {
    return true;
  }

  if (valence_vec) {
    if (lengthsq3f(valence_vec) < 0.1) {
      // null-vector
      return true;
    }

    // anchor
    auto idx0 = tcs->atmToIdx(at0);
    assert(idx0 >= 0);
    const float* va0 = tcs->coordPtr(idx0);

    // Target coordinate system
    float matT[16], mat[16];
    identity44f(matT);
    copy3f(valence_vec, matT);
    get_system1f3f(matT, matT + 4, matT + 8);
    copy3(va0, matT + 12);
    transpose44f44f(matT, mat);

    // combine with source coordinate system
    right_multiply44f44f(mat, mat1_inv);

    // transform coordset
    for (int idx = 0; idx < cs->NIndex; idx++) {
      transform44f3f(mat, &coord_orig[3 * idx], cs->coordPtr(idx));
    }
  }

  return CoordSetMerge(tcs->Obj, tcs, cs);
}

/**
 * @param at0 Atom index
 * @param index0 If different from at0: Index of atom to replace (e.g. Hydrogen)
 * @param[out] out (3x1) Normalized open valence vector
 * @return True if valid open valence vector found
 */
static bool GetTargetValenceVector(
    const CoordSet* tcs, int at0, int index0, float* out)
{
  if (!tcs) {
    return false;
  }

  if (at0 != index0) {
    // anchor
    auto ca0 = tcs->atmToIdx(at0);
    if (ca0 < 0) {
      return false;
    }

    // hydrogen (or other replaced atom)
    auto ch0 = tcs->atmToIdx(index0);
    if (ch0 < 0) {
      return false;
    }

    const float* vh0 = tcs->coordPtr(ch0);
    const float* va0 = tcs->coordPtr(ca0);
    subtract3f(vh0, va0, out);
    return true;
  }

  return CoordSetFindOpenValenceVector(tcs, at0, out);
}

/*========================================================================*/
/**
 * Merge src into I.
 *
 * Optionally creates a bond between index0 and index1, or the atoms bound to
 * index0 and index1 if they are "single geometry" atoms like hydrogens.
 *
 * @param index0 Atom index in I
 * @param index1 Atom index in src
 * @param create_bond If false, then don't create a bond, just combine the
 * objects. Implies move_flag=false.
 * @param move_flag If true, then transform the source coordinates to form a
 * reasonable bond geometry.
 */
pymol::Result<> ObjectMoleculeFuse(ObjectMolecule* I, int const index0,
    const ObjectMolecule* src, int const index1, bool create_bond,
    bool move_flag)
{
  constexpr int state1 = 0;
  constexpr bool edit = true;

  auto* G = I->G;

  if (!create_bond) {
    move_flag = false;
  }

  auto const* scs = src->getCoordSet(state1);
  if (!scs) {
    return pymol::Error("no source coordset");
  }

  auto* ai0 = I->AtomInfo.data();
  auto const* ai1 = src->AtomInfo.data();

  /**
   * If we want a bond and atm has cAtomInfoSingle geometry, then return the
   * atom index of its (only) neighbor. Otherwise, return atm.
   */
  auto const get_anchor_atm = [create_bond](
                                  const ObjectMolecule* obj, int atm) {
    if (create_bond) {
      if (obj->AtomInfo[atm].geom == cAtomInfoSingle) {
        auto neighbors = AtomNeighbors(obj, atm);
        if (neighbors.size() == 1) {
          atm = neighbors[0].atm;
        }
      }
    }
    return atm;
  };

  int const at0 = get_anchor_atm(I, index0);
  int const at1 = get_anchor_atm(src, index1);

  assert(!(at0 < 0 || at1 < 0));

  auto const anch1 = scs->atmToIdx(at1);
  if (anch1 < 0) {
    return pymol::Error("no coordinate for source anchor atom");
  }

  /* copy atoms and atom info into a 1:1 direct mapping */
  auto cs = pymol::make_unique<CoordSet>(G);
  cs->setNIndex(scs->NIndex);
  cs->enumIndices();

  auto nai = pymol::vla<AtomInfoType>(scs->NIndex);

  for (int idx = 0; idx < scs->NIndex; ++idx) {
    copy3f(scs->coordPtr(idx), cs->coordPtr(idx));
    auto atm1 = scs->IdxToAtm[idx];

    AtomInfoCopy(G, ai1 + atm1, nai + idx);

    if (edit) {
      if (atm1 == at1) {
        nai[idx].temp1 = 2; /* clear marks */
      } else {
        nai[idx].temp1 = 0; /* clear marks */
      }
    }
  }

  /* copy internal bond information */
  cs->TmpBond.reserve(src->NBond);
  cs->NTmpBond = 0;
  for (int b = 0; b != src->NBond; ++b) {
    auto const& b1 = src->Bond[b];
    auto a0 = scs->atmToIdx(b1.index[0]);
    auto a1 = scs->atmToIdx(b1.index[1]);
    if (a0 >= 0 && a1 >= 0) {
      auto& b0 = cs->TmpBond[cs->NTmpBond++];
      b0 = b1;
      b0.index[0] = a0;
      b0.index[1] = a1;
    }
  }

  // source coordinate system
  float mat1_inv[16];

  if (create_bond) {
    const float * va1 = scs->coordPtr(anch1);

    if (at0 != index0) {
      ai0[index0].deleteFlag = true;
    }

    identity44f(mat1_inv);

    if (at1 != index1) {
      auto const hydr1 = scs->atmToIdx(index1);
      if (hydr1 == -1) {
        return pymol::Error("no source attachment vector found");
      }
      nai[hydr1].deleteFlag = true;
      const float* vh1 = scs->coordPtr(hydr1);
      subtract3f(va1, vh1, mat1_inv);
    } else {
      auto found = CoordSetFindOpenValenceVector(scs, at1, mat1_inv);
      if (!found) {
        return pymol::Error("no source attachment vector found");
      }
      scale3f(mat1_inv, -1.0F, mat1_inv);
    }

    // bond length
    float const d = AtomInfoGetBondLength(G, ai0 + at0, ai1 + at1);

    // rotation matrix
    get_system1f3f(mat1_inv, mat1_inv + 4, mat1_inv + 8);

    // translation vector
    float mat1_trans[16];
    identity44f(mat1_trans);
    mat1_trans[3] = (mat1_inv[0] * d - va1[0]);
    mat1_trans[7] = (mat1_inv[1] * d - va1[1]);
    mat1_trans[11] = (mat1_inv[2] * d - va1[2]);

    right_multiply44f44f(mat1_inv, mat1_trans);

    /* set up the linking bond */
    cs->TmpLinkBond.resize(1);
    cs->NTmpLinkBond = 1;
    auto* bond = cs->TmpLinkBond.data();
    BondTypeInit2(bond, at0, anch1);
  }

  AtomInfoUniquefyNames(I, nai.data(), cs->NIndex);

  /* set up tags which will enable use to continue editing bond */

  if (edit) {
    for (int atm = 0; atm < I->NAtom; ++atm) {
      ai0[atm].temp1 = 0;
    }
    ai0[at0].temp1 = 1;
  }

  int const state0 = I->DiscreteFlag ? (ai0[at0].discrete_state - 1)
                                     : cSelectorUpdateTableAllStates;
  assert(!I->DiscreteFlag || I->DiscreteCSet[at0] == I->CSet[state0]);

  // Get target valence vectors for all relevant coordinate sets. Do this before
  // merging atoms to not leave the object molecule in an half-done state in
  // case we can't find valence vectors.
  std::vector<std::array<float, 3>> target_vectors;
  if (move_flag) {
    bool found_any_target_vector = false;
    target_vectors.resize(I->NCSet);

    for (StateIterator iter(G, nullptr, state0, I->NCSet); iter.next();) {
      float* vec = target_vectors[iter.state].data();
      if (GetTargetValenceVector(I->CSet[iter.state], at0, index0, vec)) {
        found_any_target_vector = true;
      } else {
        zero3(vec);
      }
    }

    if (!found_any_target_vector) {
      return pymol::Error("no target attachment vector found");
    }
  }

  // will free nai, cs->TmpBond and cs->TmpLinkBond
  // invalidates the ai0 pointer!
  bool ok = ObjectMoleculeMerge(
                I, std::move(nai), cs.get(), false, cAIC_AllMask, true) &&
            ObjectMoleculeExtendIndices(I, state0);
  p_return_val_if_fail(ok, pymol::Error::MEMORY);

  // Get untransformed copy of source coordinates
  pymol::vla<float> coord_orig;
  if (move_flag) {
    coord_orig = cs->Coord;
    p_return_val_if_fail(coord_orig, pymol::Error::MEMORY);
  }

  // Add coordinates into the coordinate set(s)
  for (StateIterator iter(G, nullptr, state0, I->NCSet); iter.next();) {
    const float* const vec =
        move_flag ? target_vectors[iter.state].data() : nullptr;
    ok = AddCoordinateIntoCoordSet(I->CSet[iter.state], cs.get(),
        coord_orig.data(), at0, mat1_inv, vec);
    p_return_val_if_fail(ok, pymol::Error::MEMORY);
  }

  ok = ObjectMoleculeSort(I);
  p_return_val_if_fail(ok, pymol::Error::MEMORY);

  ObjectMoleculeUpdateIDNumbers(I);

  if (at0 != index0 || at1 != index1) {
    ObjectMoleculePurge(I);
  }

  // edit the resulting bond
  if (edit) {
    int atm_pk1 = -1;
    int atm_pk2 = -1;
    for (int atm = 0; atm < I->NAtom; ++atm) {
      if (I->AtomInfo[atm].temp1 == 1) {
        atm_pk2 = atm;
      } else if (I->AtomInfo[atm].temp1 == 2) {
        atm_pk1 = atm;
      }
    }
    if (atm_pk2 >= 0 && atm_pk1 >= 0) {
      auto sele1 = pymol::string_format("%s`%d", I->Name, atm_pk1 + 1);
      auto sele2 = pymol::string_format("%s`%d", I->Name, atm_pk2 + 1);
      EditorSelect(G, sele1.data(), sele2.data(), "", "", false, true, true);
    }
  }

  return {};
}

/*========================================================================*/
int ObjectMoleculeVerifyChemistry(ObjectMolecule * I, int state)
{
  int result = false;
  const AtomInfoType *ai;
  int a;
  int flag;

  if(state < 0) {
    /* use the first defined state */
    for(a = 0; a < I->NCSet; a++) {
      if(I->CSet[a]) {
        state = a;
        break;
      }
    }
  }
  ai = I->AtomInfo;
  flag = true;
  for(a = 0; a < I->NAtom; a++) {
    if(!ai->chemFlag) {
      flag = false;
    }
    ai++;
  }
  if((!flag) && (state >= 0) && (state < I->NCSet)) {
    if(I->CSet[state]) {
      ObjectMoleculeInferChemFromBonds(I, state);
      ObjectMoleculeInferChemFromNeighGeom(I, state);
      ObjectMoleculeInferHBondFromChem(I);
      /*      ObjectMoleculeInferChemForProtein(I,0); */
    }
    flag = true;
    ai = I->AtomInfo;
    for(a = 0; a < I->NAtom; a++) {
      if(!ai->chemFlag) {
        flag = false;
        break;
      }
      ai++;
    }
  }
  if(flag)
    result = true;
  return (result);
}


/*========================================================================*/
int ObjectMoleculeAttach(ObjectMolecule * I, int index,
    pymol::vla<AtomInfoType>&& nai)
{
  int a;
  AtomInfoType *ai;
  BondType* bond;
  float v[3], v0[3], d;
  CoordSet *cs = NULL;
  int ok = false;

  ai = I->AtomInfo + index;

  ok_assert(1, cs = CoordSetNew(I->G));
  ok_assert(1, cs->Coord = pymol::vla<float>(3));

  cs->NIndex = 1;
  ok_assert(1, cs->TmpLinkBond = pymol::vla<BondType>(1));

  cs->NTmpLinkBond = 1;
  bond = cs->TmpLinkBond.data();
  BondTypeInit2(bond, index, 0);
  cs->enumIndices();

  ok_assert(1, ObjectMoleculePrepareAtom(I, index, nai.data()));
  d = AtomInfoGetBondLength(I->G, ai, nai);

  ok_assert(1, ObjectMoleculeMerge(I, std::move(nai),
        cs, false, cAIC_AllMask, true)); // will free nai and cs->TmpLinkBond
  ok_assert(1, ObjectMoleculeExtendIndices(I, -1));

  for(a = 0; a < I->NCSet; a++) {       /* add atom to each coordinate set */
    if (const auto* cs_a = I->CSet[a]) {
      CoordSetGetAtomVertex(cs_a, index, v0);
      CoordSetFindOpenValenceVector(cs_a, index, v);
      scale3f(v, d, v);
      add3f(v0, v, cs->Coord.data());
      ok_assert(1, CoordSetMerge(I, I->CSet[a], cs));
    }
  }

  ok_assert(1, ObjectMoleculeSort(I));
  ObjectMoleculeUpdateIDNumbers(I);

  ok = true;
ok_except1:
  delete cs;
  return ok;
}


/*========================================================================*/
int ObjectMoleculeFillOpenValences(ObjectMolecule * I, int index)
{
  int a;
  AtomInfoType *ai;
  int result = 0;
  int flag = true;
  float v[3], v0[3], d;
  CoordSet *cs = NULL;
  int ok = true;

  if((index >= 0) && (index <= I->NAtom)) {
    while(ok) {
      ai = I->AtomInfo + index;
      auto const nn = AtomNeighbors(I, index).size();

      assert(flag);

      if((nn >= ai->valence) || (!flag))
        break;
      flag = false;

      if (ok)
	cs = CoordSetNew(I->G);
      CHECKOK(ok, cs);
      if (ok){
	cs->Coord = pymol::vla<float>(3);
	CHECKOK(ok, cs->Coord);

	cs->NIndex = 1;
	if (ok)
	  cs->TmpLinkBond = pymol::vla<BondType>(1);
	CHECKOK(ok, cs->TmpLinkBond);
	if (ok){
	  cs->NTmpLinkBond = 1;
          auto* bond = cs->TmpLinkBond.data();
          BondTypeInit2(bond, index, 0);
	}
      }
      
      if(ok)
        cs->enumIndices();
      auto atInfo = pymol::vla<AtomInfoType>(1);
      AtomInfoType* nai = atInfo.data();
      if (ok){
	UtilNCopy(nai->elem, "H", 2);
	nai->geom = cAtomInfoSingle;
	nai->valence = 1;
	ok &= ObjectMoleculePrepareAtom(I, index, nai);
	d = AtomInfoGetBondLength(I->G, ai, nai);
	if (ok)
          ok &= ObjectMoleculeMerge(I, std::move(atInfo),
              cs, false, cAIC_AllMask, true);       /* will free nai and cs->TmpLinkBond  */
      }
      if (ok)
	ok &= ObjectMoleculeExtendIndices(I, -1);
      for(a = 0; ok &&  a < I->NCSet; a++) {   /* add atom to each coordinate set */
        if (auto* cs_a = I->CSet[a]) {
          CoordSetGetAtomVertex(cs_a, index, v0);
          CoordSetFindOpenValenceVector(cs_a, index, v);
          scale3f(v, d, v);
          add3f(v0, v, cs->Coord.data());
          ok &= CoordSetMerge(I, cs_a, cs);
        }
      }
      delete cs;
      result++;
      flag = true;
    }
  }
  ObjectMoleculeUpdateIDNumbers(I);
  return (result);
}

#define MaxOcc 100


/*========================================================================*/
/**
 * @param[out] normal (3x1) Normal vector of planar geometry at atom `atm`
 *
 * @return True if the atom has planar geometry and the normal could be computed
 *
 * @pre neighbors are defined
 */
static bool get_planer_normal(const CoordSet* cs, int const atm, float* normal)
{
  int found = false;
  int nOcc = 0;
  float occ[MaxOcc * 3];
  float v0[3], v1[3], v2[3], n0[3];

  if (CoordSetGetAtomVertex(cs, atm, v0)) {
    const auto* I = cs->Obj;

    // look for an attached non-hydrogen as a base
    for (auto const& neighbor : AtomNeighbors(I, atm)) {
      if (CoordSetGetAtomVertex(cs, neighbor.atm, v1)) {
        subtract3f(v1, v0, n0);
        normalize3f(n0);        /* n0's point away from center atom */
        copy3f(n0, occ + 3 * nOcc);
        nOcc++;
        if(nOcc == MaxOcc)      /* safety valve */
          break;
      }
    }
    switch (I->AtomInfo[atm].geom) {
    case cAtomInfoPlanar:
      if(nOcc > 1) {
        cross_product3f(occ, occ + 3, normal);
        if(nOcc > 2) {
          cross_product3f(occ, occ + 6, v2);
          if(dot_product3f(normal, v2) < 0) {
            subtract3f(normal, v2, normal);
          } else {
            add3f(normal, v2, normal);
          }
          cross_product3f(occ + 3, occ + 6, v2);
          if(dot_product3f(normal, v2) < 0) {
            subtract3f(normal, v2, normal);
          } else {
            add3f(normal, v2, normal);
          }
        }
        normalize3f(normal);
        found = true;
      }
      break;
    }
  }
  return found;
}

/*========================================================================*/
/**
 * @param atm Atom index
 * @param[out] v (3x1) Normalized open valence vector
 * @param seek (3x1) Primer for the direction, or NULL for random priming.
 * @param ignore_index Atom index of neighbor to ignore
 * @return True if valid open valence vector found
 */
bool CoordSetFindOpenValenceVector(
    const CoordSet* cs, int atm, float* v, const float* seek, int ignore_index)
{
  int nOcc = 0;
  float occ[MaxOcc * 3];
  int last_occ = -1;
  float v0[3], v1[3], n0[3] = {0.0F,0.0F,0.0F}, t[3];
  int result = false;
  float y[3], z[3];

  /* default is +X */
  v[0] = 1.0;
  v[1] = 0.0;
  v[2] = 0.0;

  if(cs) {
    const auto* I = cs->Obj;

    if (atm >= 0 && atm <= I->NAtom) {
      if (CoordSetGetAtomVertex(cs, atm, v0)) {
        const auto* ai = I->AtomInfo.data() + atm;

        // look for an attached non-hydrogen as a base
        for (auto const& neighbor : AtomNeighbors(I, atm)) {
          auto const a1 = neighbor.atm;
          if(a1 != ignore_index) {
            if (CoordSetGetAtomVertex(cs, a1, v1)) {
              last_occ = a1;
              subtract3f(v1, v0, n0);
              normalize3f(n0);  /* n0's point away from center atom */
              copy3f(n0, occ + 3 * nOcc);
              nOcc++;
              if(nOcc == MaxOcc)        /* safety valve */
                break;
            }
          }
        }
        if((!nOcc) || (nOcc > 4) || (ai->geom == cAtomInfoNone)) {
          if(!seek)
            get_random3f(v);
          else
            copy3f(seek, v);
          result = true;
        } else {
          switch (nOcc) {
          case 1:              /* only one current occupied position */
            switch (ai->geom) {
            case cAtomInfoTetrahedral:
              if(!seek) {
                get_system1f3f(occ, y, z);
                scale3f(occ, -0.334F, v);
                scale3f(z, 0.943F, t);
                add3f(t, v, v);
              } else {          /* point hydrogen towards sought vector */
                copy3f(seek, z);
                get_system2f3f(occ, z, y);
                scale3f(occ, -0.334F, v);
                scale3f(z, 0.943F, t);
                add3f(t, v, v);
              }
              result = true;
              break;
            case cAtomInfoPlanar:
              {
                if(!seek) {
                  if (last_occ >= 0 && get_planer_normal(cs, last_occ, n0)) {
                    copy3f(n0, y);
                    get_system2f3f(occ, y, z);
                  } else {
                    get_system1f3f(occ, y, z);
                  }
                  scale3f(occ, -0.500F, v);
                  scale3f(z, 0.866F, t);
                  add3f(t, v, v);
                } else {
                  copy3f(seek, z);
                  get_system2f3f(occ, z, y);
                  scale3f(occ, -0.500F, v);
                  scale3f(z, 0.866F, t);
                  add3f(t, v, v);
                }
                result = true;
              }
              break;
            case cAtomInfoLinear:
              scale3f(occ, -1.0F, v);
              result = true;
              break;
            default:
              if(!seek)
                get_random3f(v);
              else
                copy3f(seek, v);
              result = false;
              break;
            }
            break;
          case 2:              /* only two current occupied positions */
            switch (ai->geom) {
            case cAtomInfoTetrahedral:
              add3f(occ, occ + 3, t);
              get_system2f3f(t, occ, z);
              scale3f(t, -1.0F, v);
              if(seek) {
                if(dot_product3f(z, seek) < 0.0F) {
                  invert3f(z);
                }
              }
              scale3f(z, 1.41F, t);
              add3f(t, v, v);
              result = true;
              break;
            case cAtomInfoPlanar:
              add3f(occ, occ + 3, t);
              scale3f(t, -1.0F, v);
              result = true;
              break;
            default:
              if(!seek) {
                add3f(occ, occ + 3, t);
                scale3f(t, -1.0F, v);
                if(length3f(t) < 0.1)
                  get_random3f(v);
              } else
                copy3f(seek, v);
              /* hypervalent */
              result = false;
              break;
            }
            break;
          case 3:              /* only three current occupied positions */
            switch (ai->geom) {
            case cAtomInfoTetrahedral:
              add3f(occ, occ + 3, t);
              add3f(occ + 6, t, t);
              scale3f(t, -1.0F, v);
              result = true;
              break;
            default:
              if(!seek) {
                add3f(occ, occ + 3, t);
                add3f(occ + 6, t, t);
                scale3f(t, -1.0F, v);
                if(length3f(t) < 0.1)
                  get_random3f(v);
              } else
                copy3f(seek, v);
              /* hypervalent */
              result = false;
              break;
            }
            break;
          case 4:
            if(!seek)
              get_random3f(v);
            else
              copy3f(seek, v);
            /* hypervalent */
            result = false;
            break;
          }
        }
      }
    }
  }
  normalize3f(v);
  return (result);
#undef MaxOcc

}


/*========================================================================*/
void ObjectMoleculeCreateSpheroid(ObjectMolecule * I, int average)
{
  CoordSet *cs;
  int a, b, c, a0;
  SphereRec *sp;
  float *v, *v0, *s, *f, ang, min_dist, *max_sq;
  int *i;
  float *center = NULL;
  float d0[3], n0[3], d1[3], d2[3];
  float p0[3], p1[3], p2[3];
  int t0, t1, t2, bt0, bt1, bt2;
  float dp, l, *fsum = NULL;
  float spheroid_smooth;
  float spheroid_fill;
  float spheroid_ratio = 0.1F;  /* minimum ratio of width over length */
  float spheroid_minimum = 0.02F;       /* minimum size - to insure valid normals */
  int row, *count = NULL, base;
  int nRow;
  int first = 0;
  int last = 0;
  int current;
  int cscount;
  int n_state = 0;
  sp = GetSpheroidSphereRec(I->G);

  nRow = I->NAtom * sp->nDot;

  center = pymol::malloc<float>(I->NAtom * 3);
  count = pymol::malloc<int>(I->NAtom);
  fsum = pymol::malloc<float>(nRow);
  max_sq = pymol::malloc<float>(I->NAtom);

  spheroid_smooth = SettingGetGlobal_f(I->G, cSetting_spheroid_smooth);
  spheroid_fill = SettingGetGlobal_f(I->G, cSetting_spheroid_fill);
  /* first compute average coordinate */

  if(average < 1)
    average = I->NCSet;
  current = 0;
  cscount = 0;
  while(current < I->NCSet) {
    if(I->CSet[current]) {
      if(!cscount)
        first = current;
      cscount++;
      last = current + 1;
    }

    if(cscount == average || current == I->NCSet - 1) {
      PRINTFB(I->G, FB_ObjectMolecule, FB_Details)
        " ObjectMolecule: computing spheroid from states %d to %d.\n",
        first + 1, last ENDFB(I->G);

      auto spheroid = std::vector<float>(nRow);

      v = center;
      i = count;
      for(a = 0; a < I->NAtom; a++) {
        *(v++) = 0.0;
        *(v++) = 0.0;
        *(v++) = 0.0;
        *(i++) = 0;
      }

      for(b = first; b < last; b++) {
        cs = I->CSet[b];
        if(cs) {
          v = cs->Coord.data();
          for(a = 0; a < cs->NIndex; a++) {
            a0 = cs->IdxToAtm[a];
            v0 = center + 3 * a0;
            add3f(v, v0, v0);
            (*(count + a0))++;
            v += 3;
          }
        }
      }

      i = count;
      v = center;
      for(a = 0; a < I->NAtom; a++)
        if(*i) {
          (*(v++)) /= (*i);
          (*(v++)) /= (*i);
          (*(v++)) /= (*i++);
        } else {
          v += 3;
          i++;
        }

      /* now go through and compute radial distances */

      f = fsum;
      s = spheroid.data();
      for(a = 0; a < nRow; a++) {
        *(f++) = 0.0;
        *(s++) = 0.0;
      }

      v = max_sq;
      for(a = 0; a < I->NAtom; a++)
        *(v++) = 0.0;

      for(b = first; b < last; b++) {
        cs = I->CSet[b];
        if(cs) {
          v = cs->Coord.data();
          for(a = 0; a < cs->NIndex; a++) {
            a0 = cs->IdxToAtm[a];
            base = (a0 * sp->nDot);
            v0 = center + (3 * a0);
            subtract3f(v, v0, d0);      /* subtract from average */
            l = lengthsq3f(d0);
            if(l > max_sq[a0])
              max_sq[a0] = l;
            if(l > 0.0) {
              float isq = (float) (1.0 / sqrt1d(l));
              scale3f(d0, isq, n0);
              for(c = 0; c < sp->nDot; c++) {   /* average over spokes */
                dp = dot_product3f(sp->dot[c], n0);
                row = base + c;
                if(dp >= 0.0) {
                  ang = (float) ((acos(dp) / spheroid_smooth) * (cPI / 2.0));
                  if(ang > spheroid_fill)
                    ang = spheroid_fill;
                  /* take envelop to zero over that angle */
                  if(ang <= (cPI / 2.0)) {
                    dp = (float) cos(ang);
                    fsum[row] += dp * dp;
                    spheroid[row] += l * dp * dp * dp;
                  }
                }
              }
            }
            v += 3;
          }
        }
      }

      f = fsum;
      s = spheroid.data();
      for(a = 0; a < I->NAtom; a++) {
        min_dist = (float) (spheroid_ratio * sqrt(max_sq[a]));
        if(min_dist < spheroid_minimum)
          min_dist = spheroid_minimum;
        for(b = 0; b < sp->nDot; b++) {
          if(*f > R_SMALL4) {
            (*s) = (float) (sqrt1d((*s) / (*(f++))));   /* we put the "rm" in "rms" */
          } else {
            f++;
          }
          if(*s < min_dist)
            *s = min_dist;
          s++;
        }
      }

      /* set frame 0 coordinates to the average */

      cs = I->CSet[first];
      if(cs) {
        v = cs->Coord.data();
        for(a = 0; a < cs->NIndex; a++) {
          a0 = cs->IdxToAtm[a];
          v0 = center + 3 * a0;
          copy3f(v0, v);
          v += 3;
        }
      }

      /* now compute surface normals */

      auto norm = std::vector<float>(nRow * 3);
      for(a = 0; a < nRow; a++) {
        zero3f(norm.data() + a * 3);
      }
      for(a = 0; a < I->NAtom; a++) {
        base = a * sp->nDot;
        for(b = 0; b < sp->NTri; b++) {
          t0 = sp->Tri[b * 3];
          t1 = sp->Tri[b * 3 + 1];
          t2 = sp->Tri[b * 3 + 2];
          bt0 = base + t0;
          bt1 = base + t1;
          bt2 = base + t2;
          copy3f(sp->dot[t0], p0);
          copy3f(sp->dot[t1], p1);
          copy3f(sp->dot[t2], p2);
          /*      scale3f(sp->dot[t0].v,spheroid[bt0],p0);
             scale3f(sp->dot[t1].v,spheroid[bt1],p1);
             scale3f(sp->dot[t2].v,spheroid[bt2],p2); */
          subtract3f(p1, p0, d1);
          subtract3f(p2, p0, d2);
          cross_product3f(d1, d2, n0);
          normalize3f(n0);
          v = norm.data() + bt0 * 3;
          add3f(n0, v, v);
          v = norm.data() + bt1 * 3;
          add3f(n0, v, v);
          v = norm.data() + bt2 * 3;
          add3f(n0, v, v);
        }
      }

      f = norm.data();
      for(a = 0; a < I->NAtom; a++) {
        base = a * sp->nDot;
        for(b = 0; b < sp->nDot; b++) {
          normalize3f(f);
          f += 3;
        }
      }

      if(I->CSet[first]) {
        I->CSet[first]->Spheroid = std::move(spheroid);
        I->CSet[first]->SpheroidNormal = std::move(norm);
      }

      for(b = first + 1; b < last; b++) {
        delete I->CSet[b];
        I->CSet[b] = NULL;
      }

      if(n_state != first) {
        I->CSet[n_state] = I->CSet[first];
        I->CSet[first] = NULL;
      }
      n_state++;

      cscount = 0;
    }
    current++;
  }
  I->NCSet = n_state;
  FreeP(center);
  FreeP(count);
  FreeP(fsum);
  FreeP(max_sq);

  I->invalidate(cRepSphere, cRepInvProp, -1);
}


/*========================================================================*/
void ObjectMoleculeReplaceAtom(ObjectMolecule * I, int index, AtomInfoType&& ai)
{
  if((index >= 0) && (index <= I->NAtom)) {
    I->AtomInfo[index] = std::move(ai);
    I->invalidate(cRepAll, cRepInvAtoms, -1);
    /* could we put in a refinement step here? */
  }
}


/*========================================================================*/
int ObjectMoleculePrepareAtom(ObjectMolecule * I, int index, AtomInfoType * ai,
    bool uniquefy)
{
  /* match existing properties of the old atom */
  AtomInfoType *ai0;
  int ok = true;

  if((index >= 0) && (index <= I->NAtom)) {
    ai0 = I->AtomInfo + index;
    ai->resv = ai0->resv;
    ai->hetatm = ai0->hetatm;
    ai->flags = ai0->flags;

    if (!ai->geom)
      ai->geom = ai0->geom;

    ai->discrete_state = ai0->discrete_state;
    ai->q = ai0->q;
    ai->b = ai0->b;
    strcpy(ai->alt, ai0->alt);
    ai->inscode = ai0->inscode;
    LexAssign(I->G, ai->segi, ai0->segi);
    LexAssign(I->G, ai->chain, ai0->chain);
    LexAssign(I->G, ai->resn, ai0->resn);
    ai->visRep = ai0->visRep;
    ai->id = -1;
    ai->rank = -1;

    AtomInfoAssignParameters(I->G, ai);

    if (uniquefy) {
      AtomInfoUniquefyNames(I->G, I->AtomInfo, I->NAtom, ai, NULL, 1);
    }

    if((ai->elem[0] == ai0->elem[0]) && (ai->elem[1] == ai0->elem[1]))
      ai->color = ai0->color;
    else if((ai->elem[0] == 'C') && (ai->elem[1] == 0)) {
      int found = false;
      for (auto const& neighbor : AtomNeighbors(I, index)) {
        AtomInfoType const* ai1 = I->AtomInfo.data() + neighbor.atm;
        if(ai1->protons == cAN_C) {
          ai->color = ai1->color;
          found = true;
          break;
        }
      }
      if(ok && !found) {
        /* if no carbon nearby, then color according to the object color */
        ai->color = I->Color;
      }
    } else {
      AtomInfoAssignColors(I->G, ai);
    }
  }
  return ok;
}


/*========================================================================*/
int ObjectMoleculePreposReplAtom(ObjectMolecule * I, int index, AtomInfoType * ai)
{
  float v0[3], v1[3], v[3];
  float d0[3], d, n0[3];
  int cnt;
  float t[3], sum[3];
  int a;
  int ncycle;
  int ok = true;
  if (ok){
    for(a = 0; a < I->NCSet; a++) {
      if(I->CSet[a]) {
	if(ObjectMoleculeGetAtomVertex(I, a, index, v0)) {
	  copy3f(v0, v);          /* default is direct superposition */
	  ncycle = -1;
	  while(ncycle) {
	    cnt = 0;
	    zero3f(sum);
            /* look for an attached non-hydrogen as a base */
            for (auto const& neighbor : AtomNeighbors(I, index)) {
	      auto const a1 = neighbor.atm;
              auto const* ai1 = I->AtomInfo.data() + a1;
              if(ai1->protons != 1)
		if(ObjectMoleculeGetAtomVertex(I, a, a1, v1)) {
		  d = AtomInfoGetBondLength(I->G, ai, ai1);
		  subtract3f(v0, v1, n0);
		  normalize3f(n0);
		  scale3f(n0, d, d0);
		  add3f(d0, v1, t);
		  add3f(t, sum, sum);
		  cnt++;
		}
	    }
	    if(cnt) {
	      scale3f(sum, 1.0F / cnt, sum);
	      copy3f(sum, v0);
	      if((cnt > 1) && (ncycle < 0))
		ncycle = 5;
	    }
	    ncycle = abs(ncycle) - 1;
	  }
	  if(cnt)
	    copy3f(sum, v);
	  ObjectMoleculeSetAtomVertex(I, a, index, v);
	}
      }
    }
  }
  return ok;
}


/*========================================================================*/

#if 1
void ObjectMoleculeSaveUndo(ObjectMolecule * I, int state, int log)
{
  CoordSet *cs;
  PyMOLGlobals *G = I->G;
  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter] = -1;
  if(state < 0)
    state = 0;
  if(I->NCSet == 1)
    state = 0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    I->UndoCoord[I->UndoIter] = pymol::malloc<float>(cs->NIndex * 3);
    memcpy(I->UndoCoord[I->UndoIter], cs->Coord, sizeof(float) * cs->NIndex * 3);
    I->UndoState[I->UndoIter] = state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }
  I->UndoIter = cUndoMask & (I->UndoIter + 1);
  ExecutiveSetLastObjectEdited(G, I);
  if(log) {
    OrthoLineType line;
    if(SettingGetGlobal_i(I->G, cSetting_logging)) {
      sprintf(line, "cmd.push_undo(\"%s\",%d)\n", I->Name, state + 1);
      PLog(G, line, cPLog_no_flush);
    }
  }

}


/*========================================================================*/
void ObjectMoleculeUndo(ObjectMolecule * I, int dir)
{
  CoordSet *cs;
  int state;

  FreeP(I->UndoCoord[I->UndoIter]);
  I->UndoState[I->UndoIter] = -1;
  state = SceneGetState(I->G);
  if(state < 0)
    state = 0;
  if(I->NCSet == 1)
    state = 0;
  state = state % I->NCSet;
  cs = I->CSet[state];
  if(cs) {
    I->UndoCoord[I->UndoIter] = pymol::malloc<float>(cs->NIndex * 3);
    memcpy(I->UndoCoord[I->UndoIter], cs->Coord, sizeof(float) * cs->NIndex * 3);
    I->UndoState[I->UndoIter] = state;
    I->UndoNIndex[I->UndoIter] = cs->NIndex;
  }

  I->UndoIter = cUndoMask & (I->UndoIter + dir);
  if(!I->UndoCoord[I->UndoIter])
    I->UndoIter = cUndoMask & (I->UndoIter - dir);

  if(I->UndoState[I->UndoIter] >= 0) {
    state = I->UndoState[I->UndoIter];
    if(state < 0)
      state = 0;

    if(I->NCSet == 1)
      state = 0;
    state = state % I->NCSet;
    cs = I->CSet[state];
    if(cs) {
      if(cs->NIndex == I->UndoNIndex[I->UndoIter]) {
        memcpy(cs->Coord.data(), I->UndoCoord[I->UndoIter], sizeof(float) * cs->NIndex * 3);
        I->UndoState[I->UndoIter] = -1;
        FreeP(I->UndoCoord[I->UndoIter]);
        cs->invalidateRep(cRepAll, cRepInvCoord);
        SceneChanged(I->G);
      }
    }
  }
}
#endif

int ObjectMoleculeAddBond(ObjectMolecule * I, int sele0, int sele1, int order, pymol::zstring_view symop)
{
  int a1, a2;
  AtomInfoType *ai1, *ai2;
  int s1, s2;
  int c = 0;

  /* TO DO: optimize for performance -- we shouldn't be doing full
     table scans */

  ai1 = I->AtomInfo.data();
  for(a1 = 0; a1 < I->NAtom; a1++) {
    s1 = ai1->selEntry;
    if(SelectorIsMember(I->G, s1, sele0)) {
      ai2 = I->AtomInfo.data();
      for(a2 = 0; a2 < I->NAtom; a2++) {
        s2 = ai2->selEntry;
        if(SelectorIsMember(I->G, s2, sele1)) {
          if(!I->Bond){
            I->Bond = pymol::vla<BondType>(1);
	  }
          if(I->Bond) {
            BondType* bnd = I->Bond.check(I->NBond);
            BondTypeInit2(bnd, a1, a2, order);

            assert(!bnd->symop_2);
            if (!symop.empty()) {
              bnd->symop_2.reset(symop.c_str());
            }

            I->NBond++;
            c++;
            I->AtomInfo[a1].chemFlag = false;
            I->AtomInfo[a2].chemFlag = false;
            I->AtomInfo[a1].bonded = true;
            I->AtomInfo[a2].bonded = true;
          }
        }
        ai2++;
      }
    }
    ai1++;
  }
  if(c) {
    I->invalidate(cRepAll, cRepInvBondsNoNonbonded, -1);
  }
  return (c);
}

/*========================================================================*/
pymol::Result<> ObjectMoleculeAddBondByIndices(
    ObjectMolecule* I, unsigned atm1, unsigned atm2, int order)
{
  if (atm1 >= I->NAtom || atm2 >= I->NAtom) {
    return pymol::make_error("atom index out of bounds");
  }

  if (!I->Bond) {
    I->Bond = pymol::vla<BondType>(1);
  } else {
    I->Bond.check(I->NBond);
  }

  if (!I->Bond) {
    return pymol::Error::MEMORY;
  }

  auto& bnd = I->Bond[I->NBond++];
  BondTypeInit2(&bnd, atm1, atm2, order);

  I->AtomInfo[atm1].chemFlag = false;
  I->AtomInfo[atm2].chemFlag = false;
  I->AtomInfo[atm1].bonded = true;
  I->AtomInfo[atm2].bonded = true;

  I->invalidate(cRepAll, cRepInvBondsNoNonbonded, -1);

  return {};
}

/*========================================================================*/
int ObjectMoleculeAdjustBonds(ObjectMolecule * I, int sele0, int sele1, int mode,
                              int order, pymol::zstring_view symop)
{
  auto const G = I->G;
  int a0, a1;
  int cnt = 0;
  BondType *b0;
  int both;
  int a;

  if(I->Bond) {
    b0 = I->Bond.data();
    for(a = 0; a < I->NBond; a++) {
      a0 = b0->index[0];
      a1 = b0->index[1];

      both = 0;

      if (SelectorIsMember(G, I->AtomInfo[a0].selEntry, sele0) &&
          SelectorIsMember(G, I->AtomInfo[a1].selEntry, sele1)) {
        both = 2;
      } else if ( //
          SelectorIsMember(G, I->AtomInfo[a1].selEntry, sele0) &&
          SelectorIsMember(G, I->AtomInfo[a0].selEntry, sele1)) {
        std::swap(a0, a1);
        both = 2;
      }

      if(both == 2) {
        cnt++;
        switch (mode) {
        case 0:                /* cycle */
          switch(SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_editor_bond_cycle_mode)) {
          case 1: /* 1 arom 2 3 */
            switch(b0->order) {
            case 1:
              b0->order = 4;
              break;
            case 4:
              b0->order = 2;
              break;
            case 2:
              b0->order = 3;
              break;
            default:
              b0->order = 1;
              break;
            }
            break;
          case 2: /* 1 2 3 arom */
            b0->order++;
            if(b0->order > 4)
              b0->order = 1;
            break;
          default: /* old way -> 1 2 3 */
            b0->order++;
            if(b0->order > 3)
              b0->order = 1;
            break;
          }
          I->AtomInfo[a0].chemFlag = false;
          I->AtomInfo[a1].chemFlag = false;
          break;
        case 1:                /* set */
          b0->order = order;
          I->AtomInfo[a0].chemFlag = false;
          I->AtomInfo[a1].chemFlag = false;
          break;
        }

        if (!symop.empty()) {
          b0->symop_2.reset(symop.c_str());
        }
      }
      b0++;
    }
    if(cnt) {
      I->invalidate(cRepLine, cRepInvBonds, -1);
      I->invalidate(cRepCyl, cRepInvBonds, -1);
      I->invalidate(cRepNonbonded, cRepInvBonds, -1);
      I->invalidate(cRepNonbondedSphere, cRepInvBonds, -1);
      I->invalidate(cRepRibbon, cRepInvBonds, -1);
      I->invalidate(cRepCartoon, cRepInvBonds, -1);
    }
  }

  return (cnt);
}


/*========================================================================*/
int ObjectMoleculeRemoveBonds(ObjectMolecule * I, int sele0, int sele1)
{
  int a0, a1;
  int offset = 0;
  BondType *b0, *b1;
  int both;
  int s;
  int a;

  if(I->Bond) {
    offset = 0;
    b0 = I->Bond.data();
    b1 = I->Bond.data();
    for(a = 0; a < I->NBond; a++) {
      a0 = b0->index[0];
      a1 = b0->index[1];

      both = 0;
      s = I->AtomInfo[a0].selEntry;
      if(SelectorIsMember(I->G, s, sele0))
        both++;
      s = I->AtomInfo[a1].selEntry;
      if(SelectorIsMember(I->G, s, sele1))
        both++;
      if(both < 2) {            /* reverse combo */
        both = 0;
        s = I->AtomInfo[a1].selEntry;
        if(SelectorIsMember(I->G, s, sele0))
          both++;
        s = I->AtomInfo[a0].selEntry;
        if(SelectorIsMember(I->G, s, sele1))
          both++;
      }

      if(both == 2) {
        AtomInfoPurgeBond(I->G, b0);
        offset--;
        b0++;
        I->AtomInfo[a0].chemFlag = false;
        I->AtomInfo[a1].chemFlag = false;
      } else if(offset) {
        *(b1++) = *(b0++);      /* copy bond info */
      } else {
        *(b1++) = *(b0++);      /* copy bond info */
      }
    }
    if(offset) {
      I->NBond += offset;
      VLASize(I->Bond, BondType, I->NBond);
      I->invalidate(cRepLine, cRepInvBonds, -1);
      I->invalidate(cRepCyl, cRepInvBonds, -1);
      I->invalidate(cRepNonbonded, cRepInvBonds, -1);
      I->invalidate(cRepNonbondedSphere, cRepInvBonds, -1);
      I->invalidate(cRepRibbon, cRepInvBonds, -1);
      I->invalidate(cRepCartoon, cRepInvBonds, -1);
    }
  }

  return (-offset);
}


/*========================================================================*/
void ObjectMoleculePurge(ObjectMolecule * I)
{
  PyMOLGlobals *G = I->G;
  int offset = 0;

  // remove the object selection and free up any selection entries.
  // note that we don't delete atom selection members -- those may be needed in the new object.
  SelectorDelete(G, I->Name);

  // old-to-new mapping
  auto oldToNew = std::vector<int>(size_t(I->NAtom), -1);
  for (int atm = 0; atm < I->NAtom; ++atm) {
    auto &ai0 = I->AtomInfo[atm];
    if(ai0.deleteFlag) {
      AtomInfoPurge(G, &ai0);
      offset--;
      assert(oldToNew[atm] == -1);
    } else {
      if(offset) {
        I->AtomInfo[atm + offset] = std::move(ai0);
      }
      oldToNew[atm] = atm + offset;
    }
  }
  if(offset) {
    I->NAtom += offset;
    I->AtomInfo.resize(I->NAtom);
    for (int a = 0; a < I->NCSet; ++a) {
      if (auto cs = I->CSet[a])
        CoordSetAdjustAtmIdx(cs, oldToNew.data());
    }
    if (I->CSTmpl) {
      CoordSetAdjustAtmIdx(I->CSTmpl, oldToNew.data());
    }
  }

  I->updateAtmToIdx();

  // step 4, bonds

  offset = 0;
  auto b0 = I->Bond.data();
  auto b1 = I->Bond.data();
  for (int b = 0; b < I->NBond; ++b) {
    auto a0 = b0->index[0];
    auto a1 = b0->index[1];
    if(a0 < 0 || a1 < 0 || (oldToNew[a0] < 0) || (oldToNew[a1] < 0)) {
      /* deleting bond */
      AtomInfoPurgeBond(I->G, b0);
      offset--;
      b0++;
    } else {
      if(offset) {
	*b1 = *b0;
      }
      b1->index[0] = oldToNew[a0];      /* copy bond info */
      b1->index[1] = oldToNew[a1];
      b0++;
      b1++;
    }
  }
  if(offset) {
    I->NBond += offset;
    VLASize(I->Bond, BondType, I->NBond);
  }

  I->invalidate(cRepAll, cRepInvAtoms, -1);
}


/*========================================================================*/
/**
 * Determines hybridization from coordinates in those few cases where it is
 * unambiguous.
 */
int ObjectMoleculeGetAtomGeometry(const ObjectMolecule* I, int state, int at)
{

  int result = -1;
  float v0[3], v1[3], v2[3], v3[3];
  float d1[3], d2[3], d3[3];
  float cp1[3], cp2[3], cp3[3];
  float avg;
  float dp;
  auto const neighbors = AtomNeighbors(I, at);
  auto const nn = neighbors.size();
  if(nn == 4)
    result = cAtomInfoTetrahedral;
  else if(nn == 3) {
    /* check cross products */
    ObjectMoleculeGetAtomVertex(I, state, at, v0);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[0].atm, v1);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[1].atm, v2);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[2].atm, v3);
    subtract3f(v1, v0, d1);
    subtract3f(v2, v0, d2);
    subtract3f(v3, v0, d3);
    cross_product3f(d1, d2, cp1);
    cross_product3f(d2, d3, cp2);
    cross_product3f(d3, d1, cp3);
    normalize3f(cp1);
    normalize3f(cp2);
    normalize3f(cp3);
    avg = (dot_product3f(cp1, cp2) +
           dot_product3f(cp2, cp3) + dot_product3f(cp3, cp1)) / 3.0F;
    if(avg > 0.75)
      result = cAtomInfoPlanar;
    else
      result = cAtomInfoTetrahedral;
  } else if(nn == 2) {
    ObjectMoleculeGetAtomVertex(I, state, at, v0);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[0].atm, v1);
    ObjectMoleculeGetAtomVertex(I, state, neighbors[1].atm, v2);
    subtract3f(v1, v0, d1);
    subtract3f(v2, v0, d2);
    normalize3f(d1);
    normalize3f(d2);
    dp = dot_product3f(d1, d2);
    if(dp < -0.75)
      result = cAtomInfoLinear;
  }
  return (result);
}


/*========================================================================*/
static float compute_avg_center_dot_cross_fn(ObjectMolecule * I, CoordSet * cs,
                                             int n_atom, int *atix)
{
  float result = 0.0F;
  float *v[9];
  int missing_flag = false;
  {
    int a, i;
    for(i = 0; i < n_atom; i++) {
      int a1 = atix[i];
      a = cs->atmToIdx(a1);
      if(a < 0) {
        missing_flag = true;
        break;
      } else {
        v[i] = cs->coordPtr(a);
      }
    }
  }
  if(!missing_flag) {
    int i;
    float d10[3], d20[3], cp[8][3];
    float avg = 0.0;

    v[n_atom] = v[1];
    for(i = 1; i < n_atom; i++) {
      subtract3f(v[i], v[0], d10);
      subtract3f(v[i + 1], v[0], d20);
      normalize3f(d10);
      normalize3f(d20);
      cross_product3f(d10, d20, cp[i]);
      normalize3f(cp[i]);
      if(i > 1) {
        if(dot_product3f(cp[i - 1], cp[i]) < 0.0)
          invert3f(cp[i]);
      }
    }
    copy3f(cp[1], cp[n_atom]);
    for(i = 1; i < n_atom; i++) {
      avg += dot_product3f(cp[i], cp[i + 1]);
    }
    result = avg / (n_atom - 1);
  }
  return result;
}

static int verify_planer_bonds(const ObjectMolecule* I, const CoordSet* cs,
    int n_atom, const int* atix, const int* neighbor, const float* dir,
    float cutoff)
{
  int a, i;
  for(i = 0; i < n_atom; i++) {
    int a1 = atix[i];
    a = cs->atmToIdx(a1);
    if(a >= 0) {
      const float* v0 = cs->coordPtr(a);
      int n = neighbor[a1] + 1;
      int a2;
      while((a2 = neighbor[n]) >= 0) {
        n += 2;
        a = cs->atmToIdx(a2);
        if(a >= 0) {
          const float* v1 = cs->coordPtr(a);
          float d10[] = { 0.f, 0.f, 0.f };
          subtract3f(v1, v0, d10);
          normalize3f(d10);
          {
	    float dot = fabs(dot_product3f(d10, dir));
	    if(dot > cutoff) {
             switch (I->AtomInfo[a1].protons) {
              case cAN_C:
              case cAN_N:
              case cAN_O:
              case cAN_S:
                switch (I->AtomInfo[a2].protons) {
                case cAN_C:
                case cAN_N:
                case cAN_O:
                case cAN_S:
                  return 0;
                  break;
                }
                break;
              }
            }
          }
        }
      }
    }
  }
  return 1;
}

static float compute_avg_ring_dot_cross_fn(ObjectMolecule * I, CoordSet * cs,
                                           int n_atom, int *atix, float *dir)
{
  float result = 0.0F;
  float *v[9];
  int missing_flag = false;
  {
    int a, i;
    for(i = 0; i < n_atom; i++) {
      int a1 = atix[i];
      a = cs->atmToIdx(a1);
      if(a < 0) {
        missing_flag = true;
        break;
      } else {
        v[i] = cs->coordPtr(a);
      }
    }
  }
  if(!missing_flag) {
    int i;
    float d10[3], d21[3], cp[8][3];
    float avg = 0.0;

    v[n_atom] = v[0];
    v[n_atom + 1] = v[1];
    for(i = 0; i < n_atom; i++) {
      subtract3f(v[i + 0], v[i + 1], d10);
      subtract3f(v[i + 2], v[i + 1], d21);
      normalize3f(d10);
      normalize3f(d21);
      cross_product3f(d10, d21, cp[i]);
      normalize3f(cp[i]);
      if(i) {
        if(dot_product3f(cp[i - 1], cp[i]) < 0.0)
          invert3f(cp[i]);
      }
      add3f(cp[i], dir, dir);
    }
    copy3f(cp[0], cp[n_atom]);
    for(i = 0; i < n_atom; i++) {
      avg += dot_product3f(cp[i], cp[i + 1]);
    }
    result = avg / n_atom;
  }
  normalize3f(dir);
  return result;
}

typedef struct {
  int cyclic, planer, aromatic;
} ObservedInfo;

void ObjectMoleculeGuessValences(ObjectMolecule * I, int state, int *flag1, int *flag2,
                                 int reset)
{
  /* this a hacked 80% solution ...it will get things wrong, but it is
     better than nothing! */

  const float planer_cutoff = 0.96F;
  CoordSet *cs = NULL;
  ObservedInfo *obs_atom = NULL;
  ObservedInfo *obs_bond = NULL;
  int *flag = NULL;
  int warning1 = 0, warning2 = 0;

/* WORKAROUND of a possible -funroll-loops inlining optimizer bug in gcc 3.3.3 */

  float (*compute_avg_center_dot_cross) (ObjectMolecule *, CoordSet *, int, int *);
  float (*compute_avg_ring_dot_cross) (ObjectMolecule *, CoordSet *, int, int *, float *);

  compute_avg_center_dot_cross = compute_avg_center_dot_cross_fn;
  compute_avg_ring_dot_cross = compute_avg_ring_dot_cross_fn;


/* end WORKAROUND */

  if((state >= 0) && (state < I->NCSet)) {
    cs = I->CSet[state];
  }
  if(cs) {
    obs_atom = pymol::calloc<ObservedInfo>(I->NAtom);
    obs_bond = pymol::calloc<ObservedInfo>(I->NBond);
  }
  flag = pymol::calloc<int>(I->NAtom);
  if(flag) {
    if(!flag1) {
      int a, *flag_a = flag;
      const AtomInfoType *ai = I->AtomInfo.data();
      /* default behavior: only reset hetatm valences */
      for(a = 0; a < I->NAtom; a++) {
        *(flag_a++) = (ai++)->hetatm;
      }
    } else if(flag1 && flag2) {
      int a, *flag_a = flag, *flag1_a = flag1, *flag2_a = flag2;
      for(a = 0; a < I->NAtom; a++) {
        *(flag_a++) = (*(flag1_a++) || *(flag2_a++));
      }
    } else if(flag1) {
      int a, *flag_a = flag, *flag1_a = flag1;
      for(a = 0; a < I->NAtom; a++) {
        *(flag_a++) = *(flag1_a++);
      }
    }
  }
  if(!flag1)
    flag1 = flag;
  if(!flag2)
    flag2 = flag;
  if(reset) {
    /* reset chemistry information and bond orders for selected atoms */
    {
      int a;
      AtomInfoType *ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        if(flag[a])
          ai->chemFlag = 0;
        ai++;
      }
    }
    {
      int b;
      BondType *bi = I->Bond.data();
      for(b = 0; b < I->NBond; b++) {
        int at0 = bi->index[0];
        int at1 = bi->index[1];
        if((flag1[at0] && flag2[at1]) || (flag1[at1] && flag2[at0]))
          bi->order = 1;
        bi++;
      }
    }
  }
  if(cs && obs_bond && obs_atom && flag && flag1 && flag2) {
    int a;
    int const* const neighbor = I->getNeighborArray();
    AtomInfoType *atomInfo = I->AtomInfo.data();
    BondType *bondInfo = I->Bond.data();

    for(a = 0; a < I->NAtom; a++) {
      AtomInfoType *ai = atomInfo + a;
      if((!ai->chemFlag) && (flag[a])) {
        /*  determine whether or not atom participates in a planer system with 5 or 6 atoms */

        {
          int mem[9];
          int nbr[7];
          const int ESCAPE_MAX = 500;

          const int* atmToIdx = I->DiscreteFlag ? nullptr : cs->AtmToIdx.data();

          int escape_count = ESCAPE_MAX;    /* don't get bogged down with structures 
                                           that have unreasonable connectivity */
          mem[0] = a;
          nbr[0] = neighbor[mem[0]] + 1;
          while(((mem[1] = neighbor[nbr[0]]) >= 0) &&
                ((!atmToIdx) || (atmToIdx[mem[0]] >= 0))) {
            nbr[1] = neighbor[mem[1]] + 1;
            while(((mem[2] = neighbor[nbr[1]]) >= 0) &&
                  ((!atmToIdx) || (atmToIdx[mem[1]] >= 0))) {
              if(mem[2] != mem[0]) {
                nbr[2] = neighbor[mem[2]] + 1;
                while(((mem[3] = neighbor[nbr[2]]) >= 0) &&
                      ((!atmToIdx) || (atmToIdx[mem[2]] >= 0))) {
                  if(mem[3] != mem[1]) {
                    nbr[3] = neighbor[mem[3]] + 1;
                    while(((mem[4] = neighbor[nbr[3]]) >= 0) &&
                          ((!atmToIdx) || (atmToIdx[mem[3]] >= 0))) {
                      if((mem[4] != mem[2]) && (mem[4] != mem[1]) && (mem[4] != mem[0])) {
                        nbr[4] = neighbor[mem[4]] + 1;
                        while(((mem[5] = neighbor[nbr[4]]) >= 0) &&
                              ((!atmToIdx) || (atmToIdx[mem[4]] >= 0))) {
                          if(!(escape_count--))
                            goto escape;
                          if((mem[5] != mem[3]) && (mem[5] != mem[2])
                             && (mem[5] != mem[1])) {
                            if(mem[5] == mem[0]) {      /* five-cycle */
                              int i;
                              float dir[] = { 0.f, 0.f, 0.f } ;
                              float avg_dot_cross =
                                compute_avg_ring_dot_cross(I, cs, 5, mem, dir);
                              for(i = 0; i < 5; i++) {
                                obs_atom[mem[i]].cyclic = true;
                                obs_bond[neighbor[nbr[0] + 1]].cyclic = true;
                              }
                              if(avg_dot_cross > planer_cutoff) {
                                if(verify_planer_bonds
                                   (I, cs, 5, mem, neighbor, dir, 0.35F)) {
                                  for(i = 0; i < 5; i++) {
                                    obs_atom[mem[i]].planer = true;
                                    obs_bond[neighbor[nbr[i] + 1]].planer = true;
                                  }
                                }
                              }
                            }

                            nbr[5] = neighbor[mem[5]] + 1;
                            while(((mem[6] = neighbor[nbr[5]]) >= 0) &&
                                  ((!atmToIdx) || (atmToIdx[mem[5]] >= 0))) {
                              if((mem[6] != mem[4]) && (mem[6] != mem[3])
                                 && (mem[6] != mem[2]) && (mem[6] != mem[1])) {
                                if(mem[6] == mem[0]) {  /* six-cycle */
                                  int i;
                                  float dir[3] = { 0.f, 0.f, 0.f } ;
                                  float avg_dot_cross =
                                    compute_avg_ring_dot_cross(I, cs, 6, mem, dir);
                                  for(i = 0; i < 6; i++) {
                                    obs_atom[mem[i]].cyclic = true;
                                    obs_bond[neighbor[nbr[i] + 1]].cyclic = true;
                                  }
                                  if(avg_dot_cross > planer_cutoff) {
                                    if(verify_planer_bonds
                                       (I, cs, 6, mem, neighbor, dir, 0.35F)) {
                                      for(i = 0; i < 6; i++) {
                                        obs_atom[mem[i]].planer = true;
                                        obs_bond[neighbor[nbr[i] + 1]].planer = true;
                                      }
                                    }
                                  }
                                }
                              }
                              nbr[5] += 2;
                            }
                          }
                          nbr[4] += 2;
                        }
                      }
                      nbr[3] += 2;
                    }
                  }
                  nbr[2] += 2;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }

        escape:
          escape_count = ESCAPE_MAX;    /* don't get bogged down with structures 
                                           that have unreasonable connectivity */
	  warning1 = 1;
        }
      }
    }

    for(a = 0; a < I->NAtom; a++) {
      AtomInfoType *ai = atomInfo + a;
      if((!ai->chemFlag) && flag[a]) {
        ObservedInfo *ob_at = obs_atom + a;

        {
          switch (ai->protons) {
          case cAN_P:
          case cAN_S:
          case cAN_N:
          case cAN_C:
            {
              int atm[5];
              int bnd[5];
              int n = neighbor[a];
              int nn = neighbor[n++];
              if(nn > 4)
                nn = 4;
              atm[0] = a;
              {
                int i;
                for(i = 1; i <= nn; i++) {
                  atm[i] = neighbor[n];
                  bnd[i] = neighbor[n + 1];
                  n += 2;
                }
              }
              {
                int o1_at = -1, o2_at = -1, o3_at = -1, o4_at = -1;
                int o1_bd = 0, o2_bd = 0, o3_bd = 0, o4_bd = 0;
                float o1_len = 0.0F, o2_len = 0.0F;
                float n1_v[3] = { 0.0F, 0.0F, 0.0F };
                int n1_at = -1, n2_at = -1, n3_at = -1;
                int n1_bd = 0, n2_bd = 0, n3_bd = 0;
                float n1_len = 0.0F, n2_len = 0.0F, n3_len = 0.0F;
                int c1_at = -1, c2_at = -1;
                float c1_v[3] = { 0.0F, 0.0F, 0.0F };
                float *v0 = NULL;

                {
                  int idx0 = cs->atmToIdx(a);
                  if(idx0 >= 0)
                    v0 = cs->coordPtr(idx0);
                }
                {
                  int i;
                  for(i = 1; i <= nn; i++) {
                    float *v1 = NULL;

                    {
                      int idx1 = cs->atmToIdx(atm[i]);
                      if(idx1 >= 0)
                        v1 = cs->coordPtr(idx1);
                    }
                    if(v0 && v1) {
                      float diff[3];
                      subtract3f(v1, v0, diff);
                      {

                        float bond_len = length3f(diff);

                        switch (atomInfo[atm[i]].protons) {
                        case cAN_C:
                          if(c1_at < 0) {
                            c1_at = atm[i];
                            copy3f(diff, c1_v);
                          } else if(c2_at < 0) {
                            c2_at = atm[i];
                          }
                          break;
                        case cAN_O:
                          if(o1_at < 0) {
                            o1_at = atm[i];
                            o1_len = bond_len;
                            o1_bd = bnd[i];
                          } else if(o2_at < 0) {
                            o2_at = atm[i];
                            o2_len = bond_len;
                            o2_bd = bnd[i];
                          } else if(o3_at < 0) {
                            o3_at = atm[i];
                            o3_bd = bnd[i];
                          } else if(o4_at < 0) {
                            o4_at = atm[i];
                            o4_bd = bnd[i];
                          }
                          break;
                        case cAN_N:
                          if(n1_at < 0) {
                            copy3f(diff, n1_v);
                            n1_at = atm[i];
                            n1_len = bond_len;
                            n1_bd = bnd[i];
                          } else if(n2_at < 0) {
                            n2_at = atm[i];
                            n2_len = bond_len;
                            n2_bd = bnd[i];
                          } else if(n3_at < 0) {
                            n3_at = atm[i];
                            n3_len = bond_len;
                            n3_bd = bnd[i];
                          }
                          break;
                        }
                      }
                    }
                  }
                }
                {
                  float avg_dot_cross = 0.0F;

                  switch (ai->protons) {
                  case cAN_C:  /* planer carbons */
                    if(nn == 3) {
                      avg_dot_cross = compute_avg_center_dot_cross(I, cs, 4, atm);

                      if(avg_dot_cross > planer_cutoff) {

                        if((n1_at >= 0) && (o1_at >= 0) && (o2_at < 0)) {
                          /* simple amide? */
                          if(o1_len < 1.38F) {
                            if(neighbor[neighbor[o1_at]] == 1)
                              if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                                bondInfo[o1_bd].order = 2;
                          }
                        } else if((n1_at >= 0) && (o1_at >= 0) && (o2_at >= 0)) {
                          /* carbamyl */
                          if((o1_len < 1.38F) && (neighbor[neighbor[o1_at]] == 1) &&
                             (o2_len < 1.38F) && (neighbor[neighbor[o2_at]] == 1)) {
                            if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                              bondInfo[o1_bd].order = 4;
                            if((flag1[a] && flag2[o2_at]) || (flag2[a] && flag1[o2_at]))
                              bondInfo[o2_bd].order = 4;
                          } else if((o1_len < 1.38F) && (neighbor[neighbor[o1_at]] == 1)) {
                            if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                              bondInfo[o1_bd].order = 2;
                          } else if((o2_len < 1.38F) && (neighbor[neighbor[o2_at]] == 1))
                            if((flag1[a] && flag2[o2_at]) || (flag2[a] && flag1[o2_at]))
                              bondInfo[o2_bd].order = 2;
                        } else if((n1_at < 0) && (o1_at >= 0) && (o2_at < 0)) {
                          /* ketone */
                          if((o1_len < 1.31F) && (neighbor[neighbor[o1_at]] == 1)) {
                            if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                              bondInfo[o1_bd].order = 2;
                          }
                        } else if((o1_at >= 0) && (o2_at >= 0) && (n1_at < 0)) {
                          /* simple carboxylate? */
                          if((o1_len < 1.38F) && (o2_len < 1.38F) &&
                             (neighbor[neighbor[o1_at]] == 1) &&
                             (neighbor[neighbor[o2_at]] == 1)) {
                            if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                              bondInfo[o1_bd].order = 4;
                            if((flag1[a] && flag2[o2_at]) || (flag2[a] && flag1[o2_at]))
                              bondInfo[o2_bd].order = 4;
                          } else if((o1_len < 1.38F) && (neighbor[neighbor[o1_at]] == 1)) {     /* esters */
                            if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                              bondInfo[o1_bd].order = 2;
                          } else if((o2_len < 1.38F) && (neighbor[neighbor[o2_at]] == 1)) {
                            if((flag1[a] && flag2[o2_at]) || (flag2[a] && flag1[o2_at]))
                              bondInfo[o2_bd].order = 2;
                          }
                        } else if((n1_at >= 0) && (n2_at >= 0) && (n3_at < 0)
                                  && (c1_at >= 0) && (n1_len < 1.43F) && (n2_len < 1.43F)
                                  && obs_atom[c1_at].planer && obs_atom[c1_at].cyclic
                                  && (!ob_at->cyclic) && (!obs_atom[n1_at].cyclic)
                                  && (!obs_atom[n2_at].cyclic)) {
                          if((flag1[a] && flag2[n1_at]) || (flag2[a] && flag1[n1_at]))
                            bondInfo[n1_bd].order = 4;
                          atomInfo[n1_at].valence = 3;
                          atomInfo[n1_at].geom = cAtomInfoPlanar;
                          atomInfo[n1_at].chemFlag = 2;
                          if((flag1[a] && flag2[n2_at]) || (flag2[a] && flag1[n2_at]))
                            bondInfo[n2_bd].order = 4;
                          atomInfo[n2_at].valence = 3;
                          atomInfo[n2_at].geom = cAtomInfoPlanar;
                          atomInfo[n2_at].chemFlag = 2;
                        } else if((n1_at >= 0) && (n2_at >= 0) && (n3_at >= 0)) {
                          /* guanido with no hydrogens */
                          if((n1_len < 1.44F) && (n2_len < 1.44F) && (n3_len < 1.44F)) {
                            if((neighbor[neighbor[n1_at]] == 1) &&
                               (neighbor[neighbor[n2_at]] == 1) &&
                               (neighbor[neighbor[n3_at]] >= 2)) {
                              if((flag1[a] && flag2[n1_at]) || (flag2[a] && flag1[n1_at]))
                                bondInfo[n1_bd].order = 4;
                              atomInfo[n1_at].valence = 3;
                              atomInfo[n1_at].geom = cAtomInfoPlanar;
                              atomInfo[n1_at].chemFlag = 2;
                              if((flag1[a] && flag2[n2_at]) || (flag2[a] && flag1[n2_at]))
                                bondInfo[n2_bd].order = 4;
                              atomInfo[n2_at].valence = 3;
                              atomInfo[n2_at].geom = cAtomInfoPlanar;
                              atomInfo[n2_at].chemFlag = 2;
                            } else if((neighbor[neighbor[n1_at]] == 1) &&
                                      (neighbor[neighbor[n2_at]] >= 2) &&
                                      (neighbor[neighbor[n3_at]] == 1)) {
                              if((flag1[a] && flag2[n1_at]) || (flag2[a] && flag1[n1_at]))
                                bondInfo[n1_bd].order = 4;
                              atomInfo[n1_at].valence = 3;
                              atomInfo[n1_at].geom = cAtomInfoPlanar;
                              atomInfo[n1_at].chemFlag = 2;
                              if((flag1[a] && flag2[n3_at]) || (flag2[a] && flag1[n3_at]))
                                bondInfo[n3_bd].order = 4;
                              atomInfo[n3_at].valence = 3;
                              atomInfo[n3_at].geom = cAtomInfoPlanar;
                              atomInfo[n3_at].chemFlag = 2;
                            } else if((neighbor[neighbor[n1_at]] >= 2) &&
                                      (neighbor[neighbor[n2_at]] == 1) &&
                                      (neighbor[neighbor[n3_at]] == 1)) {
                              if((flag1[a] && flag2[n2_at]) || (flag2[a] && flag1[n2_at]))
                                bondInfo[n2_bd].order = 4;
                              atomInfo[n2_at].valence = 3;
                              atomInfo[n2_at].geom = cAtomInfoPlanar;
                              atomInfo[n2_at].chemFlag = 2;
                              if((flag1[a] && flag2[n3_at]) || (flag2[a] && flag1[n3_at]))
                                bondInfo[n3_bd].order = 4;
                              atomInfo[n3_at].valence = 3;
                              atomInfo[n3_at].geom = cAtomInfoPlanar;
                              atomInfo[n3_at].chemFlag = 2;
                            }
                          }
                        }
                      }
                    }
                    /* any carbon */

                    /* handle imines and nitriles */

                    if((nn >= 2) && (nn <= 3) && (n1_at >= 0) && (o1_at < 0) &&
                       (n2_at < 0) && (n1_len < 1.36F) &&
                       (!ob_at->cyclic) && (!obs_bond[n1_bd].cyclic)
                       && (!obs_atom[n1_at].planer) && ((nn == 2)
                                                        || ((nn == 3)
                                                            && (avg_dot_cross >
                                                                planer_cutoff)))) {

                      float n1_dot_cross = 1.0F;
                      int n2 = neighbor[n1_at];
                      int nn2 = neighbor[n2++];

                      {         /* check nitrogen planarity */
                        int atm2[5];
                        if(nn2 > 2) {
                          nn2 = 3;
                          atm2[0] = n1_at;
                          {
                            int i2;
                            for(i2 = 1; i2 <= nn2; i2++) {
                              atm2[i2] = neighbor[n2];
                              n2 += 2;
                            }
                          }
                          n1_dot_cross = compute_avg_center_dot_cross(I, cs, 4, atm2);
                        }
                      }
                      if(n1_dot_cross > planer_cutoff) {
                        if((flag1[a] && flag2[n1_at]) || (flag2[a] && flag1[n1_at]))
                          bondInfo[n1_bd].order = 2;
                        if((n1_len < 1.24F) && (c1_at >= 0) && (nn2 == 1)) {
                          normalize3f(n1_v);
                          normalize3f(c1_v);
                          if(dot_product3f(n1_v, c1_v) < -0.9) {
                            if((flag1[a] && flag2[n1_at]) || (flag2[a] && flag1[n1_at]))
                              bondInfo[n1_bd].order = 3;
                          }
                        }
                      }
                    }
                    break;
                  case cAN_N:
                    if(nn == 3) {
                      avg_dot_cross = compute_avg_center_dot_cross(I, cs, 4, atm);

                      if((avg_dot_cross > planer_cutoff)) {
                        if((o1_at >= 0) && (o2_at >= 0) && (o3_at < 0)) {
                          /* nitro */
                          if(neighbor[neighbor[o1_at]] == 1) {
                            if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                              bondInfo[o1_bd].order = 4;
                          }
                          if(neighbor[neighbor[o2_at]] == 1) {
                            if((flag1[a] && flag2[o2_at]) || (flag2[a] && flag1[o2_at]))
                              bondInfo[o2_bd].order = 4;
                          }
                        }
                      }
                    }
                    break;
                  case cAN_S:
                  case cAN_P:
                    if((o1_at >= 0) && (o2_at >= 0) && (o3_at >= 0) && (o4_at >= 0)) {
                      /* sulfate, phosphate */
                      int o1 = -1, o2 = -1, o3 = -1;
                      int a1 = 0;
                      if(neighbor[neighbor[o1_at]] == 1) {
                        o1 = o1_bd;
                        a1 = o1_at;
                      }
                      if(neighbor[neighbor[o2_at]] == 1) {
                        if(o1 < 0) {
                          o1 = o2_bd;
                          a1 = o2_at;
                        } else if(o2 < 0) {
                          o2 = o2_bd;
                        }
                      }
                      if(neighbor[neighbor[o3_at]] == 1) {
                        if(o1 < 0) {
                          o1 = o3_bd;
                          a1 = o3_at;
                        } else if(o2 < 0)
                          o2 = o3_bd;
                        else if(o3 < 0)
                          o3 = o3_bd;
                      }
                      if(neighbor[neighbor[o4_at]] == 1) {
                        if(o1 < 0) {
                          o1 = o4_bd;
                          a1 = o4_at;
                        } else if(o2 < 0)
                          o2 = o4_bd;
                        else if(o3 < 0)
                          o3 = o4_bd;
                      }
                      if(o2 >= 0) {
                        if((flag1[a] && flag2[a1]) || (flag2[a] && flag1[a1]))
                          bondInfo[o1].order = 2;
                        if(o2 == o2_bd) {
                          atomInfo[o2_at].formalCharge = -1;
                        } else if(o2 == o3_bd) {
                          atomInfo[o3_at].formalCharge = -1;
                        } else if(o2 == o4_bd) {
                          atomInfo[o4_at].formalCharge = -1;
                        }
                      }
                    } else if((o1_at >= 0) && (o2_at >= 0) && (o3_at >= 0) && (o4_at < 0)) {
                      /* sulfonamide */
                      int o1 = -1, o2 = -1;
                      int a1 = 0, a2 = 0;
                      if(neighbor[neighbor[o1_at]] == 1) {
                        o1 = o1_bd;
                        a1 = o1_at;
                      }
                      if(neighbor[neighbor[o2_at]] == 1) {
                        if(o1 < 0) {
                          o1 = o2_bd;
                          a1 = o2_at;
                        } else if(o2 < 0) {
                          o2 = o2_bd;
                          a2 = o2_at;
                        }
                      }
                      if(neighbor[neighbor[o3_at]] == 1) {
                        if(o1 < 0) {
                          o1 = o3_bd;
                          a1 = o3_at;
                        } else if(o2 < 0) {
                          o2 = o3_bd;
                          a2 = o3_at;
                        }
                      }
                      if(o1 >= 0) {
                        if((flag1[a] && flag2[a1]) || (flag2[a] && flag1[a1]))
                          bondInfo[o1].order = 2;
                      }
                      if(o2 >= 0) {
                        if((flag1[a] && flag2[a2]) || (flag2[a] && flag1[a2]))
                          bondInfo[o2].order = 2;
                      }
                    } else if((o1_at >= 0) && (o2_at >= 0) && (o3_at < 0)) {
                      /* sulphone */
                      if(neighbor[neighbor[o1_at]] == 1) {
                        if((flag1[a] && flag2[o1_at]) || (flag2[a] && flag1[o1_at]))
                          bondInfo[o1_bd].order = 2;
                      }
                      if(neighbor[neighbor[o2_at]] == 1) {
                        if((flag1[a] && flag2[o2_at]) || (flag2[a] && flag1[o2_at]))
                          bondInfo[o2_bd].order = 2;
                      }
                    }
                    break;
                  }
                }
              }
            }
          }
          /* this sets aromatic bonds for cyclic planer systems */

          if(ob_at->cyclic && ob_at->planer) {

            switch (ai->protons) {
            case cAN_C:
            case cAN_N:
            case cAN_O:
            case cAN_S:
              {
                int n, a0, b0;
                n = neighbor[a] + 1;
                while(1) {
                  a0 = neighbor[n];
                  if(a0 < 0)
                    break;
                  b0 = neighbor[n + 1];
                  n += 2;
                  if(obs_atom[a0].cyclic && obs_atom[a0].planer && obs_bond[b0].cyclic) {
                    obs_atom[a0].aromatic = true;
                    switch (I->AtomInfo[a0].protons) {
                    case cAN_C:
                    case cAN_N:
                    case cAN_O:
                    case cAN_S:
                      if((flag1[a] && flag2[a0]) || (flag2[a] && flag1[a0]))
                        I->Bond[b0].order = 4;
                      break;
                    }
                  }
                }
              }
              break;
            }
          }
        }
      }
    }

    /* now try to address some simple cases with aromatic nitrogens */
    for(a = 0; a < I->NAtom; a++) {
      AtomInfoType *ai = atomInfo + a;
      if((!ai->chemFlag) && (ai->protons == cAN_N) &&
         (ai->formalCharge == 0) && flag[a] &&
         obs_atom[a].cyclic && obs_atom[a].aromatic) {

        int n = neighbor[a];
        int nn = neighbor[n++];

        if(nn == 2) {           /* only two explicit neighbors */

          int mem[9];
          int nbr[7];
          const int ESCAPE_MAX = 500;

          const int* atmToIdx = I->DiscreteFlag ? nullptr : cs->AtmToIdx.data();

          int escape_count = ESCAPE_MAX;    /* don't get bogged down with structures 
                                           that have unreasonable connectivity */
          mem[0] = a;
          nbr[0] = neighbor[mem[0]] + 1;
          while(((mem[1] = neighbor[nbr[0]]) >= 0) &&
                ((!atmToIdx) || (atmToIdx[mem[0]] >= 0))) {
            nbr[1] = neighbor[mem[1]] + 1;
            while(((mem[2] = neighbor[nbr[1]]) >= 0) &&
                  ((!atmToIdx) || (atmToIdx[mem[1]] >= 0))) {
              if(mem[2] != mem[0]) {
                nbr[2] = neighbor[mem[2]] + 1;
                while(((mem[3] = neighbor[nbr[2]]) >= 0) &&
                      ((!atmToIdx) || (atmToIdx[mem[2]] >= 0))) {
                  if(mem[3] != mem[1]) {
                    nbr[3] = neighbor[mem[3]] + 1;
                    while(((mem[4] = neighbor[nbr[3]]) >= 0) &&
                          ((!atmToIdx) || (atmToIdx[mem[3]] >= 0))) {
                      if((mem[4] != mem[2]) && (mem[4] != mem[1]) && (mem[4] != mem[0])) {
                        nbr[4] = neighbor[mem[4]] + 1;
                        while(((mem[5] = neighbor[nbr[4]]) >= 0) &&
                              ((!atmToIdx) || (atmToIdx[mem[4]] >= 0))) {
                          if(!(escape_count--))
                            goto escape2; /* BUG FIX: need a new escape2, instead of 
					     mistakenly going back to the first escape,
					     which is in the loop above */
                          if((mem[5] != mem[3]) && (mem[5] != mem[2])
                             && (mem[5] != mem[1])) {
                            if(mem[5] == mem[0] && (!ai->chemFlag)) {

                              /* unassigned aromatic nitrogen-containing five-cycle */

                              /* c1ccnc1 becomes c1ccn[H]c1 */

                              if((atomInfo[mem[1]].protons == cAN_C) &&
                                 (atomInfo[mem[2]].protons == cAN_C) &&
                                 (atomInfo[mem[3]].protons == cAN_C) &&
                                 (atomInfo[mem[4]].protons == cAN_C) &&
                                 obs_atom[mem[1]].aromatic &&
                                 obs_atom[mem[2]].aromatic &&
                                 obs_atom[mem[3]].aromatic && obs_atom[mem[4]].aromatic) {
                                ai->valence = 3;
                                ai->chemFlag = 2;
                                ai->geom = cAtomInfoPlanar;
                              }

                              /* c1ncnc1 becomes c1n[H]cnc1 */

                              if((atomInfo[mem[1]].protons == cAN_C) &&
                                 (atomInfo[mem[2]].protons == cAN_N) &&
                                 (atomInfo[mem[2]].formalCharge == 0) &&
                                 (!atomInfo[mem[2]].chemFlag) &&
                                 (atomInfo[mem[3]].protons == cAN_C) &&
                                 (atomInfo[mem[4]].protons == cAN_C) &&
                                 obs_atom[mem[1]].aromatic &&
                                 obs_atom[mem[2]].aromatic &&
                                 obs_atom[mem[3]].aromatic && obs_atom[mem[4]].aromatic) {

                                int n2 = neighbor[mem[2]];
                                int nn2 = neighbor[n2++];
                                if(nn2 == 2) {  /* second nitrogen also ambiguous */
                                  ai->valence = 3;
                                  ai->chemFlag = 2;
                                  ai->geom = cAtomInfoPlanar;
                                }
                              }

                              /* c1cnnc1 becomes c1cn[H]nc1 */

                              if((atomInfo[mem[1]].protons == cAN_N) &&
                                 (atomInfo[mem[1]].formalCharge == 0) &&
                                 (!atomInfo[mem[1]].chemFlag) &&
                                 (atomInfo[mem[2]].protons == cAN_C) &&
                                 (atomInfo[mem[3]].protons == cAN_C) &&
                                 (atomInfo[mem[4]].protons == cAN_C) &&
                                 obs_atom[mem[1]].aromatic &&
                                 obs_atom[mem[2]].aromatic &&
                                 obs_atom[mem[3]].aromatic && obs_atom[mem[4]].aromatic) {

                                int n2 = neighbor[mem[1]];
                                int nn2 = neighbor[n2++];
                                if(nn2 == 2) {  /* second nitrogen also ambiguous */
                                  ai->valence = 3;
                                  ai->chemFlag = 2;
                                  ai->geom = cAtomInfoPlanar;
                                }
                              }

                            }
                          }
                          nbr[4] += 2;
                        }
                      }
                      nbr[3] += 2;
                    }
                  }
                  nbr[2] += 2;
                }
              }
              nbr[1] += 2;
            }
            nbr[0] += 2;
          }
        escape2:           /* BUG FIX: Need separate escape for this loop */
          escape_count = ESCAPE_MAX;    /* don't get bogged down with structures 
                                           that have unreasonable connectivity */
	  warning2 = 1;
        }
      }
    }
  }
  if (warning1 || warning2){
	  PRINTFB(I->G, FB_ObjectMolecule, FB_Blather)
	    " %s(%d,%d): Unreasonable connectivity in heteroatom,\n  unsuccessful in guessing valences.\n", __func__, warning1, warning2
	     ENDFB(I->G);
  }
  FreeP(obs_bond);
  FreeP(obs_atom);
  FreeP(flag);
}


/*========================================================================*/

void ObjectMoleculeInferChemForProtein(ObjectMolecule * I, int state)
{
  /* Infers chemical relations for a molecules under protein assumptions.
   * 
   * NOTE: this routine needs an all-atom model (with hydrogens!)
   * and it will make mistakes on non-protein atoms (if they haven't
   * already been assigned)
   */

  int changedFlag = true;

  /* first, try to find all amids and acids */
  while(changedFlag) {
    changedFlag = false;
    for (int a = 0; a < I->NAtom; a++) {
      auto const* const ai = I->AtomInfo.data() + a;
      if(ai->chemFlag) {
        if(ai->geom == cAtomInfoPlanar)
          if(ai->protons == cAN_C) {
            auto const neighbors = AtomNeighbors(I, a);
            if (neighbors.size() > 1) {
              AtomInfoType* ai1 = nullptr;
              for (auto const& neighbor : neighbors) {
                auto* const ai0 = I->AtomInfo.data() + neighbor.atm;
                if((ai0->protons == cAN_O) && (!ai0->chemFlag)) {
                  ai1 = ai0;    /* found candidate carbonyl */
                  break;
                }
              }
              if (ai1) {
                for (auto const& neighbor : neighbors) {
                  auto* const ai0 = I->AtomInfo.data() + neighbor.atm;
                  if (ai0 != ai1) {
                    if(ai0->protons == cAN_O) {
                      if(!ai0->chemFlag) {
                        ai0->chemFlag = true;   /* acid */
                        ai0->geom = cAtomInfoPlanar;
                        ai0->valence = 1;
                        ai1->chemFlag = true;
                        ai1->geom = cAtomInfoPlanar;
                        ai1->valence = 1;
                        changedFlag = true;
                        break;
                      }
                    } else if(ai0->protons == cAN_N) {
                      if(!ai0->chemFlag) {
                        ai0->chemFlag = true;   /* amide N */
                        ai0->geom = cAtomInfoPlanar;
                        ai0->valence = 3;
                        ai1->chemFlag = true;   /* amide =O */
                        ai1->geom = cAtomInfoPlanar;
                        ai1->valence = 1;
                        changedFlag = true;
                        break;
                      } else if(ai0->geom == cAtomInfoPlanar) {
                        ai1->chemFlag = true;   /* amide =O */
                        ai1->geom = cAtomInfoPlanar;
                        ai1->valence = 1;
                        changedFlag = true;
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
      }
    }
  }
  /* then handle aldehydes and amines (partial amides - both missing a valence) */

  changedFlag = true;
  while(changedFlag) {
    changedFlag = false;
    for (int a = 0; a < I->NAtom; a++) {
      auto* const ai = I->AtomInfo.data() + a;
      if(!ai->chemFlag) {
        if(ai->protons == cAN_C) {
          auto const neighbors = AtomNeighbors(I, a);
          if (neighbors.size() > 1) {
            for (auto const& neighbor : neighbors) {
              auto* const ai0 = I->AtomInfo.data() + neighbor.atm;
              if((ai0->protons == cAN_O) && (!ai0->chemFlag)) { /* =O */
                ai->chemFlag = true;
                ai->geom = cAtomInfoPlanar;
                ai->valence = 1;
                ai0->chemFlag = true;
                ai0->geom = cAtomInfoPlanar;
                ai0->valence = 3;
                changedFlag = true;
                break;
              }
            }
          }
        } else if(ai->protons == cAN_N) {
          if((!ai->chemFlag) || ai->geom != cAtomInfoLinear) {
            if(ai->formalCharge == 0) {
              ai->chemFlag = true;
              ai->geom = cAtomInfoPlanar;
              ai->valence = 3;
            }
          }
        }
      }
    }
  }

}


/*========================================================================*/
/**
 * Assigns:
 * - geom
 * - valence
 * - chemFlag
 */
void ObjectMoleculeInferChemFromNeighGeom(ObjectMolecule * I, int state)
{
  /* infers chemical relations from neighbors and geometry 
   * NOTE: very limited in scope */

  int a;
  int changedFlag = true;
  int geom;
  int carbonVal[10];

  AtomInfoType *ai, *ai2;

  carbonVal[cAtomInfoTetrahedral] = 4;
  carbonVal[cAtomInfoPlanar] = 3;
  carbonVal[cAtomInfoLinear] = 2;

  while(changedFlag) {
    changedFlag = false;
    for(a = 0; a < I->NAtom; a++) {
      ai = I->AtomInfo + a;
      if(!ai->chemFlag) {
        geom = ObjectMoleculeGetAtomGeometry(I, state, a);
        switch (ai->protons) {
        case cAN_K:
          ai->chemFlag = 1;
          ai->geom = cAtomInfoNone;
          ai->valence = 0;
          break;
        case cAN_H:
        case cAN_F:
        case cAN_I:
        case cAN_Br:
          ai->chemFlag = 1;
          ai->geom = cAtomInfoSingle;
          ai->valence = 1;
          break;
        case cAN_O: {
          auto const neighbors = AtomNeighbors(I, a);
          auto const nn = neighbors.size();
          if(nn != 1) {         /* water, hydroxy, ether */
            ai->chemFlag = 1;
            ai->geom = cAtomInfoTetrahedral;
            ai->valence = 2;
          } else {              /* hydroxy or carbonyl? check carbon geometry */
            ai2 = I->AtomInfo.data() + neighbors[0].atm;
            if(ai2->chemFlag) {
              if((ai2->geom == cAtomInfoTetrahedral) || (ai2->geom == cAtomInfoLinear)) {
                ai->chemFlag = 1;       /* hydroxy */
                ai->geom = cAtomInfoTetrahedral;
                ai->valence = 2;
              }
            }
          }
          break;
        }
        case cAN_C:
          if(geom >= 0) {
            ai->geom = geom;
            ai->valence = carbonVal[geom];
            ai->chemFlag = true;
          } else {
            auto const neighbors = AtomNeighbors(I, a);
            auto const nn = neighbors.size();
            if(nn == 1) {       /* only one neighbor */
              ai2 = I->AtomInfo.data() + neighbors[0].atm;
              if(ai2->chemFlag && (ai2->geom == cAtomInfoTetrahedral)) {
                ai->chemFlag = true;    /* singleton carbon bonded to tetC must be tetC */
                ai->geom = cAtomInfoTetrahedral;
                ai->valence = 4;
              }
            }
          }
          break;
        case cAN_N:
          if(geom == cAtomInfoPlanar) {
            ai->chemFlag = true;
            ai->geom = cAtomInfoPlanar;
            ai->valence = 3;
          } else if(geom == cAtomInfoTetrahedral) {
            ai->chemFlag = true;
            ai->geom = cAtomInfoTetrahedral;
            ai->valence = 4;
          }
          break;
        case cAN_S: {
          auto const nn = AtomNeighbors(I, a).size();
          if(nn == 4) {         /* sulfone */
            ai->chemFlag = true;
            ai->geom = cAtomInfoTetrahedral;
            ai->valence = 4;
          } else if(nn == 3) {  /* suloxide */
            ai->chemFlag = true;
            ai->geom = cAtomInfoTetrahedral;
            ai->valence = 3;
          } else if(nn == 2) {  /* thioether */
            ai->chemFlag = true;
            ai->geom = cAtomInfoTetrahedral;
            ai->valence = 2;
          }
          break;
        }
        case cAN_Cl:
          ai->chemFlag = 1;
          if(ai->formalCharge == 0) {
            ai->geom = cAtomInfoSingle;
            ai->valence = 1;
          } else {
            ai->geom = cAtomInfoNone;
            ai->valence = 0;
          }
          break;
        }
        if(ai->chemFlag)
          changedFlag = true;
      }
    }
  }
}


/*========================================================================*/
/**
 * Assigns, based on valence, geom, formalCharge, and bonded atoms:
 * - hb_donor
 * - hb_acceptor
 */
void ObjectMoleculeInferHBondFromChem(ObjectMolecule * I)
{
  int a;
  AtomInfoType *ai;
  int a1;
  int has_hydro;
  /* initialize accumulators on uncategorized atoms */

  const lexborrow_t lex_pseudo = LexBorrow(I->G, "pseudo");

  ai = I->AtomInfo.data();
  for(a = 0; a < I->NAtom; a++) {
    auto const neighbors = AtomNeighbors(I, a);
    int nn = neighbors.size();
    ai->hb_donor = false;
    ai->hb_acceptor = false;

    has_hydro = (nn < ai->valence);     /* implicit hydrogens? */

    if(!has_hydro) {
      /* explicit hydrogens? */
      switch (ai->protons) {
      case cAN_N:
      case cAN_O:
        for (auto const& neighbor : neighbors) {
          a1 = neighbor.atm;
          if(I->AtomInfo[a1].protons == 1) {
            has_hydro = true;
            break;
          }
          if (I->AtomInfo[a1].name == lex_pseudo && --nn < ai->valence) {
            has_hydro = true;
            break;
          }
        }
        break;
      }
    }

    switch (ai->protons) {
      /* cat-ions, lewis acids etc. */
    case cAN_Fe:
    case cAN_Ca:
    case cAN_Cu:
    case cAN_K:
    case cAN_Na:
    case cAN_Mg:
    case cAN_Zn:
    case cAN_Hg:
    case cAN_Sr:
    case cAN_Ba:
      ai->hb_donor = true;
      break;
    case cAN_N:
      if(has_hydro)
        ai->hb_donor = true;
      else {
        int delocalized = false;
        int has_double_bond = false;
        int neighbor_has_double_bond = false;
        int mem[3];
        int nbr[3];
        int const* neighbor = I->getNeighborArray();

        mem[0] = a;
        nbr[0] = neighbor[mem[0]] + 1;
        while(((mem[1] = neighbor[nbr[0]]) >= 0)) {
          int b_order = I->Bond[neighbor[nbr[0] + 1]].order;
          if(b_order > 1) {     /* any double/triple/aromatic bonds? */
            delocalized = true;
          }
          if(b_order == 2) {
            has_double_bond = true;
          }
          nbr[1] = neighbor[mem[1]] + 1;
          while(((mem[2] = neighbor[nbr[1]]) >= 0)) {
            if(mem[2] != mem[0]) {
              int b_order2 = I->Bond[neighbor[nbr[1] + 1]].order;
              if(b_order2 == 2) {
                neighbor_has_double_bond = true;
              }
            }
            nbr[1] += 2;
          }
          nbr[0] += 2;
        }
        if((ai->formalCharge <= 0) && delocalized && (nn < 3)) {
          /* delocalized nitrogen can likely serve as an acceptor */
          ai->hb_acceptor = true;
        }
        if(delocalized && (neighbor_has_double_bond) && (!has_double_bond) &&
           (ai->geom == cAtomInfoPlanar) && (nn == 2) && (ai->formalCharge >= 0)) {
          /* there's a fair chance of a resonance structure with this nitrogen as a donor */
          ai->hb_donor = true;
        }
        if((ai->geom != cAtomInfoPlanar) && (nn == 3) && (ai->formalCharge >= 0)
           && (!delocalized)) {
          /* tertiary amine case -- assume potential donor */
          ai->hb_donor = true;
        }
      }
      break;
    case cAN_O:
      if(ai->formalCharge <= 0)
        ai->hb_acceptor = true;
      if(has_hydro)
        ai->hb_donor = true;
      else {
        int has_double_bond = false;
        int neighbor_has_aromatic_bond = false;
        int mem[3];
        int nbr[3];
        auto* const neighbor = I->getNeighborArray();

        mem[0] = a;
        nbr[0] = neighbor[mem[0]] + 1;
        while(((mem[1] = neighbor[nbr[0]]) >= 0)) {
          int b_order = I->Bond[neighbor[nbr[0] + 1]].order;
          if(b_order == 2) {
            has_double_bond = true;
          }
          nbr[1] = neighbor[mem[1]] + 1;
          while(((mem[2] = neighbor[nbr[1]]) >= 0)) {
            if(mem[2] != mem[0]) {
              int b_order2 = I->Bond[neighbor[nbr[1] + 1]].order;
              if(b_order2 == 4) {
                neighbor_has_aromatic_bond = true;
              }
            }
            nbr[1] += 2;
          }
          nbr[0] += 2;
        }

        if(has_double_bond && neighbor_has_aromatic_bond && (ai->formalCharge >= 0)) {
          /* allow for phenolic resonance structures (and the like) */
          ai->hb_donor = true;
        }
      }
      break;
    }
    ai++;
  }

}


/*========================================================================*/
/**
 * Assigns:
 * - geom
 * - valence
 * - chemFlag
 */
void ObjectMoleculeInferChemFromBonds(ObjectMolecule * I, int state)
{

  int a, b;
  const BondType *b0;
  AtomInfoType *ai, *ai0, *ai1 = NULL;
  int a0, a1;
  int expect, order;
  int changedFlag;
  /* initialize accumulators on uncategorized atoms */

  ai = I->AtomInfo.data();
  for(a = 0; a < I->NAtom; a++) {
    if(!ai->chemFlag) {
      ai->geom = 0;
      ai->valence = 0;
    }
    ai++;
  }

  // Ignore "pseudo" atoms (Desmond virtual sites)
  const lexborrow_t lex_pseudo = LexBorrow(I->G, "pseudo");
  std::vector<unsigned char> pseudo_neighbor_count(
      lex_pseudo == LEX_BORROW_NOTFOUND ? 0 : I->NAtom);

  /* find maximum bond order for each atom */

  b0 = I->Bond;
  for(b = 0; b < I->NBond; b++) {
    a0 = b0->index[0];
    a1 = b0->index[1];
    ai0 = I->AtomInfo + a0;
    ai1 = I->AtomInfo + a1;
    order = b0->order;
    b0++;

    // count "pseudo" neighbors
    if (!pseudo_neighbor_count.empty()) {
      if (ai0->name == lex_pseudo) {
        pseudo_neighbor_count[a1] += 1;
      }
      if (ai1->name == lex_pseudo) {
        pseudo_neighbor_count[a0] += 1;
      }
    }

    if(!ai0->chemFlag) {
      if(order > ai0->geom)
        ai0->geom = order;
      ai0->valence += order;
    }
    if(!ai1->chemFlag) {
      if(order > ai1->geom)
        ai1->geom = order;
      ai1->valence += order;
    }
    if(order == 3) {
      /* override existing chemistry * this is a temp fix to a pressing problem...
         we need to rethink the chemisty assignment ordering (should bond
         information come first? */
      ai0->geom = cAtomInfoLinear;
      ai1->geom = cAtomInfoLinear;
      if(ai0->chemFlag != 2) {
        switch (ai0->protons) {
        case cAN_C:
          ai0->valence = 2;
          break;
        default:
          ai0->valence = 1;
        }
        ai0->chemFlag = true;
      }
      if(ai1->chemFlag != 2) {
        switch (ai1->protons) {
        case cAN_C:
          ai1->valence = 2;
          break;
        default:
          ai1->valence = 1;
        }
        ai1->chemFlag = true;
      }
    } else if(order == 4) {
      ai0->geom = cAtomInfoPlanar;
      ai1->geom = cAtomInfoPlanar;
      if(ai0->chemFlag != 2) {
        switch (ai0->protons) {
        case cAN_O:
          ai0->valence = 1;
          break;
        case cAN_N:
          ai0->valence = (ai0->formalCharge == 1) ? 3 : 2; // TODO wrong for neutral -NH-
          break;
        case cAN_C:
          ai0->valence = 3;
          break;
        case cAN_S:
          ai0->valence = 2;
          break;
        default:
          ai0->valence = 4;
        }
        ai0->chemFlag = true;
      }
      if(ai1->chemFlag != 2) {
        switch (ai1->protons) {
        case cAN_O:
          ai1->valence = 1;
          break;
        case cAN_N:
          ai1->valence = (ai1->formalCharge == 1) ? 3 : 2; // TODO wrong for neutral -NH-
          break;
        case cAN_C:
          ai1->valence = 3;
          break;
        default:
          ai1->valence = 1;
        }
        ai1->chemFlag = true;
      }
    }
  }

  /* now set up valences and geometries */

  ai = I->AtomInfo.data();
  for(a = 0; a < I->NAtom; a++) {
    if(!ai->chemFlag) {
      expect = AtomInfoGetExpectedValence(I->G, ai);
      int nn = AtomNeighbors(I, a).size();

      // don't count "pseudo" atom neighbors
      if (!pseudo_neighbor_count.empty()) {
        nn -= pseudo_neighbor_count[a];
      }

      if(ai->geom == 3) {
        ai->geom = cAtomInfoLinear;
        switch (ai->protons) {
        case cAN_C:
          ai->valence = 2;
          break;
        default:
          ai->valence = 1;
        }
        ai->chemFlag = true;
      } else {
        if(expect < 0)
          expect = -expect;     /* for now, just ignore this issue */
        /*      printf("%d %d %d %d\n",ai->geom,ai->valence,nn,expect); */
        if(ai->valence == expect) {     /* sum of bond orders equals valence */
          ai->chemFlag = true;
          ai->valence = nn;
          switch (ai->geom) {   /* max bond order observed */
          case 0:
            ai->geom = cAtomInfoNone;
            break;
          case 2:
            ai->geom = cAtomInfoPlanar;
            break;
          case 3:
            ai->geom = cAtomInfoLinear;
            break;
          default:
            if(expect == 1)
              ai->geom = cAtomInfoSingle;
            else
              ai->geom = cAtomInfoTetrahedral;
            break;
          }
        } else if(ai->valence < expect) {       /* missing a bond */
          ai->chemFlag = true;
          ai->valence = nn + (expect - ai->valence);
          switch (ai->geom) {
          case 2:
            ai->geom = cAtomInfoPlanar;
            break;
          case 3:
            ai->geom = cAtomInfoLinear;
            break;
          default:
            if(expect == 1)
              ai->geom = cAtomInfoSingle;
            else
              ai->geom = cAtomInfoTetrahedral;
            break;
          }
        } else if(ai->valence > expect) {
          ai->chemFlag = true;
          ai->valence = nn;
          switch (ai->geom) {
          case 2:
            ai->geom = cAtomInfoPlanar;
            break;
          case 3:
            ai->geom = cAtomInfoLinear;
            break;
          default:
            if(expect == 1)
              ai->geom = cAtomInfoSingle;
            else
              ai->geom = cAtomInfoTetrahedral;
            break;
          }
          if(nn > 3)
            ai->geom = cAtomInfoTetrahedral;
        }
      }
    }
    ai++;
  }

  /* now go through and make sure conjugated amines are planer */
  changedFlag = true;
  while(changedFlag) {
    changedFlag = false;
    ai = I->AtomInfo.data();
    for(a = 0; a < I->NAtom; a++) {
      if(ai->chemFlag) {
        if(ai->protons == cAN_N)
          if(ai->formalCharge < 1)
            if(ai->geom == cAtomInfoTetrahedral) {
              /* search for uncharged tetrahedral nitrogen */
              for (auto const& neighbor : AtomNeighbors(I, a)) {
                auto const* ai0 = I->AtomInfo.data() + neighbor.atm;
                if((ai0->chemFlag) && (ai0->geom == cAtomInfoPlanar) &&
                   ((ai0->protons == cAN_C) || (ai0->protons == cAN_N))) {
                  ai->geom = cAtomInfoPlanar;   /* found probable delocalization */
                  if(ai->formalCharge == 0)
                    ai->valence = 3;      /* just in case... */
                  changedFlag = true;
                  break;
                }
              }
            }
      }
      ai++;
    }
  }

  /* now go through and make sure conjugated anions are planer */
  changedFlag = true;
  while(changedFlag) {
    changedFlag = false;
    ai = I->AtomInfo.data();
    for(a = 0; a < I->NAtom; a++) {
      if(ai->chemFlag) {
        if(ai->protons == cAN_O)
          if(ai->formalCharge == -1)
            if((ai->geom == cAtomInfoTetrahedral) || (ai->geom == cAtomInfoSingle)) {
              /* search for anionic tetrahedral oxygen */
              for (auto const& neighbor : AtomNeighbors(I, a)) {
                auto const* ai0 = I->AtomInfo + neighbor.atm;
                if((ai0->chemFlag) && (ai0->geom == cAtomInfoPlanar) &&
                   ((ai0->protons == cAN_C) || (ai0->protons == cAN_N))) {
                  ai->geom = cAtomInfoPlanar;   /* found probable delocalization */
                  changedFlag = true;
                  break;
                }
              }
            }
      }
      ai++;
    }
  }

}


/*========================================================================*/
int ObjectMoleculeTransformSelection(ObjectMolecule * I, int state,
                                     int sele, const float *matrix, int log,
                                     const char *sname, int homogenous, int global)
{
  /* called from "translate [5,5,5], objSele" */
  /* if sele == -1, then the whole object state is transformed */
  PyMOLGlobals *G = I->G;
  int a, s;
  int flag = false;
  CoordSet *cs;
  const AtomInfoType *ai;
  int logging;
  int all_states = false, inp_state;
  int ok = true;
  float homo_matrix[16], tmp_matrix[16];
  const float* input_matrix = matrix;

  inp_state = state;
  if(state == -2)
    state = ObjectGetCurrentState(I, false);
  if(state < 0) {
    all_states = true;
    state = -1;
  }
  PRINTFD(G, FB_ObjectMolecule)
    "ObjMolTransSele-Debug: state %d\n", state ENDFD;
  while(1) {
    if(all_states) {
      state++;
      if(state >= I->NCSet)
        break;
    }
    if(state < I->NCSet) {
      cs = I->CSet[state];
      if(cs) {
        int use_matrices = SettingGet_i(G, I->Setting.get(),
                                        NULL, cSetting_matrix_mode);
        if(use_matrices<0) use_matrices = 0;

        if(global &&!homogenous) {      /* convert matrix to homogenous */
          convertTTTfR44f(matrix, homo_matrix);
          matrix = homo_matrix;
          input_matrix = homo_matrix;
          homogenous = true;
        }

        if(global &&((use_matrices && !cs->Matrix.empty()) || I->TTTFlag)) {
          /* if input coordinates are in the global system,
             they may need to be converted to local coordinates */

          matrix = input_matrix;

          /* global to object */

          if(I->TTTFlag) {
            float ttt[16];
            if(matrix != tmp_matrix) {
              copy44f(matrix, tmp_matrix);
            }
            convertTTTfR44f(I->TTT, ttt);
            {
              float ttt_inv[16];
              invert_special44f44f(ttt, ttt_inv);
              left_multiply44f44f(ttt_inv, tmp_matrix);
              right_multiply44f44f(tmp_matrix, ttt);
            }
            matrix = tmp_matrix;
          }

          /* object to state */

          if(use_matrices) {
            if(!cs->Matrix.empty()) {
              double tmp_double[16], tmp_inv[16];
              copy44f44d(matrix, tmp_double);
              invert_special44d44d(cs->Matrix.data(), tmp_inv);
              left_multiply44d44d(tmp_inv, tmp_double);
              right_multiply44d44d(tmp_double, cs->Matrix.data());
              copy44d44f(tmp_double, tmp_matrix);
              matrix = tmp_matrix;
            }
          }
        }

        if(sele >= 0) {         /* transforming select atoms */
          ai = I->AtomInfo;
          for(a = 0; a < I->NAtom; a++) {
            s = ai->selEntry;
            if(!(ai->protekted == 1))
              if(SelectorIsMember(G, s, sele)) {
                if(homogenous)
                  CoordSetTransformAtomR44f(cs, a, matrix);
                else
                  CoordSetTransformAtomTTTf(cs, a, matrix);
                flag = true;
              }
            ai++;
          }
        } else {                /* transforming whole coordinate set */
          if(!use_matrices) {
            ai = I->AtomInfo;
            for(a = 0; a < I->NAtom; a++) {
              if(!(ai->protekted == 1)) {
                if(homogenous)
                  CoordSetTransformAtomR44f(cs, a, matrix);
                else
                  CoordSetTransformAtomTTTf(cs, a, matrix);
              }
              ai++;
            }
            flag = true;
            CoordSetRecordTxfApplied(cs, matrix, homogenous);
          } else {
            ObjectMoleculeTransformState44f(I, state, matrix, false, homogenous, false);
          }
        }
        if(flag) {
          cs->invalidateRep(cRepAll, cRepInvCoord);
          ExecutiveUpdateCoordDepends(G, I);
        }
      }
    }
    if(!all_states)
      break;
  }

  if(log) {
    OrthoLineType line;
    WordType sele_str = ",'";
    logging = SettingGetGlobal_i(G, cSetting_logging);
    if(sele >= 0) {
      strcat(sele_str, sname);
    }
    strcat(sele_str, "'");
    switch (logging) {
    case cPLog_pml:
      sprintf(line,
              "_ cmd.transform_object('%s',[\\\n_ %15.9f,%15.9f,%15.9f,%15.9f,\\\n_ %15.9f,%15.9f,%15.9f,%15.9f,\\\n_ %15.9f,%15.9f,%15.9f,%15.9f,\\\n_ %15.9f,%15.9f,%15.9f,%15.9f\\\n_     ],%d,%d%s,%d)\n",
              I->Name,
              matrix[0], matrix[1], matrix[2], matrix[3],
              matrix[4], matrix[5], matrix[6], matrix[7],
              matrix[8], matrix[9], matrix[10], matrix[11],
              matrix[12], matrix[13], matrix[14], matrix[15],
              inp_state + 1, 0, sele_str, homogenous);
      PLog(G, line, cPLog_no_flush);
      break;
    case cPLog_pym:

      sprintf(line,
              "cmd.transform_object('%s',[\n%15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f,\n%15.9f,%15.9f,%15.9f,%15.9f\n],%d,%d%s,%d)\n",
              I->Name,
              matrix[0], matrix[1], matrix[2], matrix[3],
              matrix[4], matrix[5], matrix[6], matrix[7],
              matrix[8], matrix[9], matrix[10], matrix[11],
              matrix[12], matrix[13], matrix[14], matrix[15],
              inp_state + 1, 0, sele_str, homogenous);
      PLog(G, line, cPLog_no_flush);
      break;
    default:
      break;
    }
  }
  return (ok);
}


/*========================================================================*/
int ObjectMoleculeGetAtomIndex(const ObjectMolecule* I, SelectorID_t sele)
{
  if(sele < 0)
    return (-1);

  for (int a = 0; a < I->NAtom; a++) {
    auto const s = I->AtomInfo[a].selEntry;
    if(SelectorIsMember(I->G, s, sele))
      return (a);
  }
  return (-1);
}


/*========================================================================*/
/**
 * Update all AtomInfoType.bonded
 *
 * Cost: O(NAtom + NBond)
 */
void ObjectMoleculeUpdateNonbonded(ObjectMolecule * I)
{
  int a;
  const BondType *b;
  AtomInfoType *ai;
  int nAtom = I->NAtom;
  int nBond = I->NBond;

  ai = I->AtomInfo.data();

  for(a = 0; a < nAtom; a++)
    (ai++)->bonded = false;

  b = I->Bond;
  ai = I->AtomInfo.data();
  for(a = 0; a < nBond; a++) {
    ai[b->index[0]].bonded = true;
    ai[b->index[1]].bonded = true;
    b++;
  }
}


/**
 * Get the raw ObjectMolecule::Neighbor storage structure, which is generated
 * from the ObjectMolecule::Bond array. This structure is cached, to clear the
 * cache, call ObjectMolecule::invalidate(level=cRepInvBonds).
 *
 * Direct usage of this array should be avoided, use the AtomNeighbors() class
 * instead which is a higher-level and more readable interface to the same
 * data.
 *
 * Changed in PyMOL 2.1.1: Ignore zero-order bonds (PYMOL-3025)
 *
 * @return Pointer to cached array, or NULL if memory allocation failed
 *
 * @verbatim
     0       list offset for atom 0 = n
     1       list offset for atom 1 = n + m + 1
     ...
     n-1     list offset for atom n-1

     n       count for atom 0 
     n+1     neighbor of atom 0
     n+2     bond index
     n+3     neighbor of atom 0
     n+4     bond index
     ...
     n+m     -1 terminator for atom 0

     n+m+1   count for atom 1
     n+m+2   neighbor of atom 1
     n+m+3   bond index
     n+m+4   neighbor of atom 1
     n+m+5   bond index
     etc.

     NOTE: all atoms have an offset and a terminator whether or not they have any bonds:

     Here's an example of an ALA
     [ offet for atom1, offset for atom2, ..., offset for atom9, \
     {10, 16, 26, 32, 36, 46, 50, 54, 58, 62, \
     (atom1)  nbonds=2, neighbor index1=5, bond index=1, neighbor index2=1, bond index=0, null terminator=-1 \
     2, 5, 1, 1, 0, -1, 
     (atom2) nbonds=4, nbr index1=6, bond id=4; nbr idx2=4, bond idx=3; nbr index3=2, bond idx=2; nbr index=0, bond index=0, null=-1 \
     4, 6, 4, 4, 3, 2, 2, 0, 0, -1,
     (atom3) ...
     2, 3, 5, 1, 2, -1, 1, 2, 5, -1,
     4, 9, 8, 8, 7, 7, 6, 1, 3, -1,
     1, 0, 1, -1,
     1, 1, 4, -1,
     1, 4, 6, -1,
     1, 4, 7, -1,
     1, 4, 8, -1 }
   @endverbatim

   */
int const* ObjectMolecule::getNeighborArray() const
{
  auto const I = this;

  /* If no neighbors have been calculated, fill in the table */
  if(!I->Neighbor) {

    /* Create/check the VLA */
    auto const size = (I->NAtom * 3) + (I->NBond * 4);

    const_cast<ObjectMolecule*>(this)->Neighbor.reset(new int[size]);
    auto* const neighbor = this->Neighbor.get();

    p_return_val_if_fail(neighbor, nullptr);

    /* initialize; zero out neighbors */
    std::fill_n(neighbor, I->NAtom, 0);

    /* count neighbors for each atom */
    const auto* bnd = I->Bond.data();
    for(int b = 0; b < I->NBond; b++) {
      if (bnd->order && !bnd->hasSymOp()) {
        neighbor[bnd->index[0]]++;
        neighbor[bnd->index[1]]++;
      }
      bnd++;
    }

    /* set up offsets and list terminators */
    auto c = I->NAtom;
    for (int a = 0; a < I->NAtom; a++) {
      auto const d = neighbor[a];       /* get number of neighbors */
      neighbor[c] = d;       /* store neighbor count */
      neighbor[a] = c + d + d + 1;   /* set initial position to end of list, we'll fill backwards */
      neighbor[neighbor[a]] = -1; /* store terminator */
      c += d + d + 2;
    }

    /* now load neighbors in a sequential list for each atom (reverse order) */
    bnd = I->Bond.data();
    for(int b = 0; b < I->NBond; b++, ++bnd) {
      auto const l0 = bnd->index[0];
      auto const l1 = bnd->index[1];

      if (!bnd->order || bnd->hasSymOp())
        continue;

      neighbor[l0]--;
      neighbor[neighbor[l0]] = b; /* store bond indices (for I->Bond) */
      neighbor[l0]--;
      neighbor[neighbor[l0]] = l1;        /* store neighbor references (I->AtomInfo, etc.) */

      neighbor[l1]--;
      neighbor[neighbor[l1]] = b; /* store bond indices (for I->Bond) */
      neighbor[l1]--;
      neighbor[neighbor[l1]] = l0;        /* store neighbor references (I->AtomInfo, etc.) */
    }

    // adjust down to point to the count, not the first entry
    for (int a = 0; a < I->NAtom; a++) {
      if (neighbor[a] >= 0)
        --neighbor[a];
    }
  }
  return this->Neighbor.get();
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
static CoordSet *ObjectMoleculeChemPyModel2CoordSet(PyMOLGlobals * G,
                                                    PyObject * model,
                                                    AtomInfoType ** atInfoPtr)
{
  int nAtom, nBond;
  int a, c;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL, *ai;
  float *f;
  BondType *ii, *bond = NULL;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);
  int hetatm;
  int ignore_ids;
  PyObject *atomList = NULL;
  PyObject *bondList = NULL;
  PyObject *atom = NULL;
  PyObject *bnd = NULL;
  PyObject *index = NULL;
  PyObject *crd = NULL;
  PyObject *tmp = NULL;

  ignore_ids = !SettingGetGlobal_b(G, cSetting_preserve_chempy_ids);

  nAtom = 0;
  nBond = 0;
  if(atInfoPtr)
    atInfo = *atInfoPtr;

  atomList = PyObject_GetAttrString(model, "atom");
  if(atomList && PyList_Check(atomList))
    nAtom = PyList_Size(atomList);
  else
    ok = ErrMessage(G, __func__, "can't get atom list");

  if(ok) {
    coord = VLAlloc(float, 3 * nAtom);
    if(atInfo)
      VLACheck(atInfo, AtomInfoType, nAtom);
  }

  if(ok) {

    f = coord;
    for(a = 0; a < nAtom; a++) {

      atom = PyList_GetItem(atomList, a);
      if(!atom)
        ok = ErrMessage(G, __func__, "can't get atom");
      crd = PyObject_GetAttrString(atom, "coord");
      if(!crd)
        ok = ErrMessage(G, __func__, "can't get coordinates");
      else {
        for(c = 0; c < 3; c++) {
          tmp = PySequence_GetItem(crd, c);
          if(tmp)
            ok = PConvPyObjectToFloat(tmp, f++);
          Py_XDECREF(tmp);
          if(!ok) {
            ErrMessage(G, __func__, "can't read coordinates");
            break;
          }
        }
      }
      Py_XDECREF(crd);

      ai = atInfo + a;
      ai->id = a;               /* chempy models are zero-based */
      if(!ignore_ids) {
        if(ok) {                /* get chempy atom id if extant */
          if(PTruthCallStr(atom, "has", "id")) {
            tmp = PyObject_GetAttrString(atom, "id");
            if(tmp)
              ok = PConvPyObjectToInt(tmp, &ai->id);
            if(!ok)
              ErrMessage(G, __func__,
                         "can't read atom identifier");
            Py_XDECREF(tmp);
          } else {
            ai->id = -1;
          }
        }
      }
      ai->rank = a;             /* ranks are always zero-based */

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "name");
        if(tmp) {
          WordType tmp_word;

          ok = PConvPyObjectToStrMaxClean(tmp, tmp_word, sizeof(WordType) - 1);
          AtomInfoCleanAtomName(tmp_word);
          ai->name = LexIdx(G, tmp_word);
        }
        if(!ok)
          ErrMessage(G, __func__, "can't read name");
        Py_XDECREF(tmp);
      }

      if(ok) {
        ai->textType = 0;
        if(PTruthCallStr(atom, "has", "text_type")) {
          tmp = PyObject_GetAttrString(atom, "text_type");
          if(tmp) {
            OrthoLineType temp;
            ok = PConvPyObjectToStrMaxClean(tmp, temp, sizeof(OrthoLineType) - 1);
            ai->textType = LexIdx(G, temp);
          }
          if(!ok)
            ErrMessage(G, __func__, "can't read text_type");
          Py_XDECREF(tmp);
        }
      }

      if(ok) {
        ai->custom = 0;
        if(PTruthCallStr(atom, "has", "custom")) {
          tmp = PyObject_GetAttrString(atom, "custom");
          if(tmp) {
            OrthoLineType temp;
            ok = PConvPyObjectToStrMaxClean(tmp, temp, sizeof(OrthoLineType) - 1);
            ai->custom = LexIdx(G, temp);
          }
          if(!ok)
            ErrMessage(G, __func__, "can't read custom");
          Py_XDECREF(tmp);
        }
      }

      if(ok) {
        if(PTruthCallStr(atom, "has", "vdw")) {
          tmp = PyObject_GetAttrString(atom, "vdw");
          if(tmp)
            ok = PConvPyObjectToFloat(tmp, &ai->vdw);
          if(!ok)
            ErrMessage(G, __func__, "can't read vdw radius");
          Py_XDECREF(tmp);
        } else {
          ai->vdw = 0.0f;
        }
      }
      if(ok) {
        if(PTruthCallStr(atom, "has", "bohr")) {
          tmp = PyObject_GetAttrString(atom, "bohr");
          if(tmp)
            ok = PConvPyObjectToFloat(tmp, &ai->elec_radius);
          if(!ok)
            ErrMessage(G, __func__,
                       "can't read elec. radius");
          Py_XDECREF(tmp);
        } else {
          ai->elec_radius = 0.0F;
        }
      }

      if(ok) {
        if(PTruthCallStr(atom, "has", "stereo")) {
          tmp = PyObject_GetAttrString(atom, "stereo");
          if(tmp) {
            char tmp_char;
            ok = PConvPyObjectToChar(tmp, &tmp_char);
            ai->stereo = tmp_char;
          }
          if(!ok)
            ErrMessage(G, __func__, "can't read stereo");
          Py_XDECREF(tmp);
        } else {
          ai->stereo = 0;
        }
      }

      if(ok) {
        if(PTruthCallStr(atom, "has", "numeric_type")) {
          tmp = PyObject_GetAttrString(atom, "numeric_type");
          if(tmp)
            ok = PConvPyObjectToInt(tmp, &ai->customType);
          if(!ok)
            ErrMessage(G, __func__,
                       "can't read numeric_type");
          Py_XDECREF(tmp);
        } else {
          ai->customType = cAtomInfoNoType;
        }
      }

      if(ok) {
        if(PTruthCallStr(atom, "has", "formal_charge")) {
          tmp = PyObject_GetAttrString(atom, "formal_charge");
          if(tmp)
            ok = PConvPyObjectToChar(tmp, (char *) &ai->formalCharge);
          if(!ok)
            ErrMessage(G, __func__,
                       "can't read formal_charge");
          Py_XDECREF(tmp);
        } else {
          ai->formalCharge = 0;
        }
      }

      if(ok) {
        if(PTruthCallStr(atom, "has", "partial_charge")) {
          tmp = PyObject_GetAttrString(atom, "partial_charge");
          if(tmp)
            ok = PConvPyObjectToFloat(tmp, &ai->partialCharge);
          if(!ok)
            ErrMessage(G, __func__,
                       "can't read partial_charge");
          Py_XDECREF(tmp);
        } else {
          ai->partialCharge = 0.0;
        }
      }

      if(ok) {
        if(PTruthCallStr(atom, "has", "flags")) {
          tmp = PyObject_GetAttrString(atom, "flags");
          if(tmp)
            ok = PConvPyObjectToInt(tmp, (int *) &ai->flags);
          if(!ok)
            ErrMessage(G, __func__, "can't read flags");
          Py_XDECREF(tmp);
        } else {
          ai->flags = 0;
        }
      }

      if(ok) {
        char buf[8] = "";
        tmp = PyObject_GetAttrString(atom, "resn");
        if(tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, buf, sizeof(buf) - 1);
        ai->resn = LexIdx(G, buf);
        if(!ok)
          ErrMessage(G, __func__, "can't read resn");
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "ins_code");
        if(tmp) {
          ResIdent tmp_ins_code;
          ok = PConvPyObjectToStrMaxClean(tmp, tmp_ins_code, sizeof(ResIdent) - 1);
          if(!ok)
            ErrMessage(G, __func__, "can't read ins_code");
          else if(tmp_ins_code[0] != '?') {
            ai->setInscode(tmp_ins_code[0]);
          }
        }
        Py_XDECREF(tmp);
      }
      if(ok) {
        if(PTruthCallStr(atom, "has", "resi_number")) {
          tmp = PyObject_GetAttrString(atom, "resi_number");
          if(tmp)
            ok = PConvPyObjectToInt(tmp, &ai->resv);
          if(!ok)
            ErrMessage(G, __func__, "can't read resi_number");
          Py_XDECREF(tmp);
        }
      }

      if(ok) {
        OrthoLineType temp = "";
        tmp = PyObject_GetAttrString(atom, "segi");
        if(tmp)
          PConvPyObjectToStrMaxClean(tmp, temp, sizeof(OrthoLineType) - 1);
        ai->segi = LexIdx(G, temp);
        if(!ok)
          ErrMessage(G, __func__, "can't read segi");
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "b");
        if(tmp)
          ok = PConvPyObjectToFloat(tmp, &ai->b);
        if(!ok)
          ErrMessage(G, __func__, "can't read b value");
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "u");
        if(tmp) {
          ok = PConvPyObjectToFloat(tmp, &ai->b);
          if(!ok)
            ErrMessage(G, __func__, "can't read u value");
          else
            ai->b *= (8 * cPI * cPI);   /* B-factor = 8 pi^2 U-factor */
        } else {
          assert(PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_AttributeError));
          PyErr_Clear();
        }
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "u_aniso");
        if(tmp) {
          float u[6];
          if(PConvPyListToFloatArrayInPlace(tmp, u, 6)) {
            // only allocate if not all zero
            if(u[0] || u[1] || u[2] || u[3] || u[4] || u[5])
              std::copy_n(u, 6, ai->get_anisou());
          }
          Py_DECREF(tmp);
        } else {
          assert(PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_AttributeError));
          PyErr_Clear();
        }
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "q");
        if(tmp)
          ok = PConvPyObjectToFloat(tmp, &ai->q);
        if(!ok)
          ErrMessage(G, __func__, "can't read occupancy");
        Py_XDECREF(tmp);
      }

      if(ok) {
        OrthoLineType temp = "";
        tmp = PyObject_GetAttrString(atom, "chain");
        if(tmp)
          PConvPyObjectToStrMaxClean(tmp, temp, sizeof(OrthoLineType) - 1);
        ai->chain = LexIdx(G, temp);
        if(!ok)
          ErrMessage(G, __func__, "can't read chain");
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "hetatm");
        if(tmp)
          ok = PConvPyObjectToInt(tmp, &hetatm);
        if(!ok)
          ErrMessage(G, __func__, "can't read hetatm");
        else {
          ai->hetatm = hetatm;
          if(!PTruthCallStr(atom, "has", "flags")) {
            if(ai->hetatm)
              ai->flags = cAtomFlag_ignore;
          }
        }
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "alt");
        if(tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, ai->alt, sizeof(Chain) - 1);
        if(!ok)
          ErrMessage(G, __func__,
                     "can't read alternate conformation");
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "symbol");
        if(tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, ai->elem, sizeof(ElemName) - 1);
        if(!ok)
          ErrMessage(G, __func__, "can't read symbol");
        Py_XDECREF(tmp);
      }

      if(ok) {
        tmp = PyObject_GetAttrString(atom, "ss");
        if(tmp)
          ok = PConvPyObjectToStrMaxClean(tmp, ai->ssType, sizeof(SSType) - 1);
        if(!ok)
          ErrMessage(G, __func__,
                     "can't read secondary structure");
        Py_XDECREF(tmp);
      }

      if(ok && PyObject_HasAttrString(atom, "label")) {
        tmp = PyObject_GetAttrString(atom, "label");
        if(tmp) {
          if (PyString_Check(tmp))
            ai->label = LexIdx(G, PyString_AsString(tmp));
        }
        Py_XDECREF(tmp);
      }

      if(ok && PyObject_HasAttrString(atom, "sphere_scale")) {
        tmp = PyObject_GetAttrString(atom, "sphere_scale");
        if(tmp) {
          float value;
          if(PConvPyFloatToFloat(tmp, &value)) {
            SettingSet(G, cSetting_sphere_scale, value, ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if(ok && PyObject_HasAttrString(atom, "cartoon_color")) {
        tmp = PyObject_GetAttrString(atom, "cartoon_color");
        if(tmp) {
          int color_index;
          ok = PConvPyObjectToInt(tmp, &color_index);
          if(!ok)
            ErrMessage(G, __func__, "bad cartoon color info");
          else {
            SettingSet(G, cSetting_cartoon_color, color_index, ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if(ok && PyObject_HasAttrString(atom, "cartoon_transparency")) {
        tmp = PyObject_GetAttrString(atom, "cartoon_transparency");
        if(tmp) {
          float alpha_val;
          ok = PConvPyObjectToFloat(tmp, &alpha_val);
          if(!ok)
            ErrMessage(G, __FUNCTION__, "bad alpha value");
          else {
            SettingSet(G, cSetting_cartoon_transparency, alpha_val, ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if(ok && PyObject_HasAttrString(atom, "cartoon_trgb")) {
        tmp = PyObject_GetAttrString(atom, "cartoon_trgb");
        if(tmp) {
          unsigned int trgb;
          ok = PConvPyObjectToInt(tmp, (signed int *) &trgb);
          if(!ok)
            ErrMessage(G, __func__, "bad cartoon color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            SettingSet(G, cSetting_cartoon_color, ColorGetIndex(G, color_name), ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if(ok && PyObject_HasAttrString(atom, "label_trgb")) {
        tmp = PyObject_GetAttrString(atom, "label_trgb");
        if(tmp) {
          unsigned int trgb;
          ok = PConvPyObjectToInt(tmp, (signed int *) &trgb);
          if(!ok)
            ErrMessage(G, __func__, "bad label color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            SettingSet(G, cSetting_label_color, ColorGetIndex(G, color_name), ai);
          }
        }
        Py_XDECREF(tmp);
      }
      if(ok && PyObject_HasAttrString(atom, "ribbon_color")) {
        tmp = PyObject_GetAttrString(atom, "ribbon_color");
        if(tmp) {
          int color_index;
          ok = PConvPyObjectToInt(tmp, &color_index);
          if(!ok)
            ErrMessage(G, __func__, "bad ribbon color info");
          else {
            SettingSet(G, cSetting_ribbon_color, color_index, ai);
          }
        }
        Py_XDECREF(tmp);
      }

      if(ok && PyObject_HasAttrString(atom, "ribbon_trgb")) {
        tmp = PyObject_GetAttrString(atom, "ribbon_trgb");
        if(tmp) {
          unsigned int trgb;
          ok = PConvPyObjectToInt(tmp, (signed int *) &trgb);
          if(!ok)
            ErrMessage(G, __func__, "bad cartoon color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
            SettingSet(G, cSetting_ribbon_color, ColorGetIndex(G, color_name), ai);
          }
        }
        Py_XDECREF(tmp);
      }

      atInfo[a].visRep = auto_show;

      if(ok && PyObject_HasAttrString(atom, "visible")) {
        unsigned int vis;
        tmp = PyObject_GetAttrString(atom, "visible");
        if(tmp) {
          ok = PConvPyObjectToInt(tmp, (signed int *) &vis);
          if(!ok)
            ErrMessage(G, __func__, "bad visibility info");
          else {
            atInfo[a].visRep = vis;

            if (!(ai->flags & cAtomFlag_class)) {
              ai->flags |= cAtomFlag_inorganic; // suppress auto_show_classified
            }
          }
        }
        Py_XDECREF(tmp);
      }

      if(ok && atInfo) {
        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);
      }

      if(ok && PyObject_HasAttrString(atom, "trgb")) {
        /* were we given a Transparency-Red-Blue-Green value? */
        unsigned int trgb;
        tmp = PyObject_GetAttrString(atom, "trgb");
        if(tmp) {
          ok = PConvPyObjectToInt(tmp, (signed int *) &trgb);
          if(!ok)
            ErrMessage(G, __func__, "bad color info");
          else {
            char color_name[24];
            sprintf(color_name, "0x%08x", trgb);
	    atInfo[a].color = ColorGetIndex(G, color_name);
          }
          Py_DECREF(tmp);
        }
      }

      if(!ok)
        break;

      assert(!PyErr_Occurred());
    }
  }

  if(nAtom) {
    bondList = PyObject_GetAttrString(model, "bond");
    if(bondList && PyList_Check(bondList))
      nBond = PyList_Size(bondList);
    else
      ok = ErrMessage(G, __func__, "can't get bond list");

    if(ok) {
      bond = VLACalloc(BondType, nBond);
      ii = bond;
      for(a = 0; a < nBond; a++) {
        bnd = PyList_GetItem(bondList, a);
        if(!bnd)
          ok = ErrMessage(G, __func__, "can't get bond");
        index = PyObject_GetAttrString(bnd, "index");
        if(!index)
          ok =
            ErrMessage(G, __func__, "can't get bond indices");
        else {
          for(c = 0; c < 2; c++) {
            tmp = PyList_GetItem(index, c);
            if(tmp)
              ok = PConvPyObjectToInt(tmp, &ii->index[c]);
            if(!ok) {
              ErrMessage(G, __func__,
                         "can't read coordinates");
              break;
            }
          }
        }
        if(ok) {
          int order = 0;
          tmp = PyObject_GetAttrString(bnd, "order");
          if(tmp)
            ok = PConvPyObjectToInt(tmp, &order);
          if(!ok)
            ErrMessage(G, __func__, "can't read bond order");
          Py_XDECREF(tmp);
          ii->order = order;
        }

        auto const set_symop = [&](const char* key) {
          if (ok && PyObject_HasAttrString(bnd, key)) {
            ok = false;
            if (auto tmp = PyObject_GetAttrString(bnd, key)) {
              if (auto const* code = PyUnicode_AsUTF8(tmp)) {
                ii->symop_2.reset(code);
                ok = true;
              }
              Py_DECREF(tmp);
            }
          }
        };

        set_symop("symmetry_2");

        if(ok && PyObject_HasAttrString(bnd, "stick_radius")) {
          tmp = PyObject_GetAttrString(bnd, "stick_radius");
          if(tmp) {
            float value;
            if(PConvPyFloatToFloat(tmp, &value)) {
              SettingSet(G, cSetting_stick_radius, value, ii);
            }
          }
          Py_XDECREF(tmp);
        }
        if(ok && PyObject_HasAttrString(bnd, "valence")) {
          tmp = PyObject_GetAttrString(bnd, "valence");
          if(tmp) {
            int value;
            if(PConvPyIntToInt(tmp, &value)) {
              SettingSet(G, cSetting_valence, value, ii);
            }
          }
          Py_XDECREF(tmp);
        }

        Py_XDECREF(index);
        ii++;
      }
    }
  }

  Py_XDECREF(atomList);
  Py_XDECREF(bondList);

  if(ok) {
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);

#ifdef _PYMOL_IP_PROPERTIES
#endif
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
  }
  if(atInfoPtr)
    *atInfoPtr = atInfo;

  if(PyErr_Occurred())
    PyErr_Print();
  return (cset);
}
#endif


/*========================================================================*/
ObjectMolecule *ObjectMoleculeLoadChemPyModel(PyMOLGlobals * G,
                                              ObjectMolecule * I,
                                              PyObject * model, int frame, int discrete)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  CoordSet *cset = NULL;
  auto atInfo = pymol::vla<AtomInfoType>(10);
  int ok = true;
  int isNew = true;
  unsigned int nAtom = 0;
  int fractional = false;
  int connect_mode = -1;
  int auto_bond = false;
  PyObject *tmp, *mol;

  if(!I)
    isNew = true;
  else
    isNew = false;

  if(ok) {

    if(isNew) {
      I = (ObjectMolecule *) new ObjectMolecule(G, discrete);
      std::swap(atInfo, I->AtomInfo);
      isNew = true;
    } else {
      if (discrete)
	ObjectMoleculeSetDiscrete(G, I, true);
    }

    if(isNew) {
      I->Color = AtomInfoUpdateAutoColor(G);
    }
    cset = ObjectMoleculeChemPyModel2CoordSet(G, model, &atInfo);

    if(!cset)
      ok = false;
    else {
      mol = PyObject_GetAttrString(model, "molecule");
      if(mol) {
        if(PyObject_HasAttrString(mol, "title")) {
          tmp = PyObject_GetAttrString(mol, "title");
          if(tmp) {
            UtilNCopy(cset->Name, PyString_AsString(tmp), sizeof(WordType));
            Py_DECREF(tmp);
            if(!strcmp(cset->Name, "untitled")) /* ignore untitled */
              cset->Name[0] = 0;
          }
        }
        Py_DECREF(mol);
      }
      if(PyObject_HasAttrString(model, "spheroid") &&
         PyObject_HasAttrString(model, "spheroid_normals")) {
        tmp = PyObject_GetAttrString(model, "spheroid");
        if(tmp) {
          PConvFromPyObject(G, tmp, cset->Spheroid);
          Py_DECREF(tmp);
        }
        tmp = PyObject_GetAttrString(model, "spheroid_normals");
        if(tmp) {
          PConvFromPyObject(G, tmp, cset->SpheroidNormal);
          Py_DECREF(tmp);
        }
      }
      if(PyObject_HasAttrString(model, "spacegroup") &&
         PyObject_HasAttrString(model, "cell")) {
        CSymmetry *symmetry = new CSymmetry(G);
        if(symmetry) {
          tmp = PyObject_GetAttrString(model, "spacegroup");
          if(tmp) {
            const char *tmp_str = NULL;
            if(PConvPyStrToStrPtr(tmp, &tmp_str)) {
              symmetry->setSpaceGroup(tmp_str);
            }
            Py_DECREF(tmp);
          }
          tmp = PyObject_GetAttrString(model, "cell");
          if(tmp) {
            float cell[6];
            if(PConvPyListToFloatArrayInPlace(tmp, cell, 6)) {
              symmetry->Crystal.setDims(cell);
              symmetry->Crystal.setAngles(cell + 3);
            }
            Py_DECREF(tmp);
          }
          cset->Symmetry = std::unique_ptr<CSymmetry>(symmetry);
        }
      }
      if(PyObject_HasAttrString(model, "fractional")) {
        tmp = PyObject_GetAttrString(model, "fractional");
        if(tmp) {
          int tmp_int = 0;
          if(PConvPyIntToInt(tmp, &tmp_int)) {
            fractional = tmp_int;
          }
          Py_DECREF(tmp);
        }
      }
      if(PyObject_HasAttrString(model, "connect_mode")) {
        tmp = PyObject_GetAttrString(model, "connect_mode");
        if(tmp) {
          int tmp_int = 0;
          if(PConvPyIntToInt(tmp, &tmp_int)) {
            auto_bond = true;
            connect_mode = tmp_int;
          }
          Py_DECREF(tmp);
        }
      }
      nAtom = cset->NIndex;
    }
  }
  /* include coordinate set */
  if(ok) {
    if(frame < 0)
      frame = I->NCSet;

    if(I->DiscreteFlag && atInfo) {
      unsigned int a;
      int fp1 = frame + 1;
      AtomInfoType *ai = atInfo.data();
      for(a = 0; a < nAtom; a++) {
        (ai++)->discrete_state = fp1;
      }
    }

    cset->Obj = I;
    cset->enumIndices();
    cset->invalidateRep(cRepAll, cRepInvRep);
    if(isNew) {
      std::swap(I->AtomInfo, atInfo);
    } else {
      ObjectMoleculeMerge(I, std::move(atInfo), cset, false, cAIC_AllMask, true);
    }
    if(isNew)
      I->NAtom = nAtom;
    VLACheck(I->CSet, CoordSet *, frame);
    if(I->NCSet <= frame)
      I->NCSet = frame + 1;
    delete I->CSet[frame];
    I->CSet[frame] = cset;
    if (fractional && cset->Symmetry) {
      CoordSetFracToReal(cset, &cset->Symmetry->Crystal);
    }
    if(ok && isNew)
      ok &= ObjectMoleculeConnect(I, cset, auto_bond, connect_mode);
    if(cset->Symmetry && (!I->Symmetry)) {
      I->Symmetry.reset(new CSymmetry(*cset->Symmetry));
    }
    SceneCountFrames(G);
    if (ok)
      ok &= ObjectMoleculeExtendIndices(I, frame);
    if (ok)
      ok &= ObjectMoleculeSort(I);
    if (ok){
      ObjectMoleculeUpdateIDNumbers(I);
      ObjectMoleculeUpdateNonbonded(I);
    }
  }
  return (I);
#endif
}


/*========================================================================*/
/**
 * Update coordinates of an exisiting coordset or insert/append a new
 * coordset.
 *
 * coords:      flat coordinate array of length NIndex * 3
 * coords_len:  must be NIndex * 3
 * frame:       coordset index
 */
ObjectMolecule *ObjectMoleculeLoadCoords(PyMOLGlobals * G, ObjectMolecule * I,
                                         const float * coords, int coords_len, int frame)
{
  CoordSet *cset = NULL;
  int a;
  bool is_new = false;

  if(frame < 0) {
    frame = I->NCSet;
  } else if (frame < I->NCSet) {
    cset = I->CSet[frame];
  }

  if (!cset) {
    // template coordinate set, if available
    cset = I->CSTmpl;

    // find any coordinate set
    for(a = 0; !cset && a < I->NCSet; ++a)
      cset = I->CSet[a];
    ok_assert(1, cset);
    cset = CoordSetCopy(cset);
    is_new = true;
  }

  // check atom count
  if(coords_len != cset->NIndex * 3) {
    ErrMessage(G, "LoadCoords", "atom count mismatch");
    ok_raise(1);
  }

  // copy coordinates
  for(a = 0; a < coords_len; ++a) {
    cset->Coord[a] = coords[a];
  }

  cset->invalidateRep(cRepAll, cRepInvRep);

  // include coordinate set
  if (is_new) {
    VLACheck(I->CSet, CoordSet *, frame);
    if(I->NCSet <= frame)
      I->NCSet = frame + 1;
    I->CSet[frame] = cset;
    SceneCountFrames(G);
  }

  // success
  return (I);

  // error handling
ok_except1:
  if(is_new && cset)
    delete cset;
  ErrMessage(G, "LoadCoords", "failed");
  return NULL;
}

/*========================================================================*/
/* see above... but look up object by name
 */
ObjectMolecule *ObjectMoleculeLoadCoords(PyMOLGlobals * G, const char * name,
                                         const float * coords, int coords_len, int frame)
{
  pymol::CObject * cobj = ExecutiveFindObjectByName(G, name);

  if(!cobj || cobj->type != cObjectMolecule) {
    ErrMessage(G, "LoadCoords", "named object molecule not found.");
    return NULL;
  }

  return ObjectMoleculeLoadCoords(G, (ObjectMolecule *) cobj,
      coords, coords_len, frame);
}


/*========================================================================*/
/* see above... but take coordinates from Python list
 *
 * coords: 2d Python float sequence with shape (NIndex, 3)
 */
ObjectMolecule *ObjectMoleculeLoadCoords(PyMOLGlobals * G, ObjectMolecule * I,
                                         PyObject * coords, int frame)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  CoordSet *cset = NULL;
  int a, b, l;
  PyObject *v, *w;
  float *f;
  bool is_new = false;

  if(!PySequence_Check(coords)) {
    ErrMessage(G, "LoadCoords", "passed argument is not a sequence");
    ok_raise(1);
  }

  if(frame < 0) {
    frame = I->NCSet;
  } else if (frame < I->NCSet) {
    cset = I->CSet[frame];
  }

  if (!cset) {
    // template coordinate set, if available
    cset = I->CSTmpl;

    // find any coordinate set
    for(a = 0; !cset && a < I->NCSet; ++a)
      cset = I->CSet[a];
    ok_assert(1, cset);
    cset = CoordSetCopy(cset);
    is_new = true;
  }

  // check atom count
  l = PySequence_Size(coords);
  if(l != cset->NIndex) {
    ErrMessage(G, "LoadCoords", "atom count mismatch");
    ok_raise(1);
  }

  // copy coordinates
  f = cset->Coord.data();
  for(a = 0; a < l; a++) {
    v = PySequence_ITEM(coords, a);

    for(b = 0; b < 3; b++) {
      if(!(w = PySequence_GetItem(v, b)))
        break;

      f[a * 3 + b] = (float) PyFloat_AsDouble(w);
      Py_DECREF(w);
    }

    Py_DECREF(v);
    ok_assert(2, !PyErr_Occurred());
  }

  cset->invalidateRep(cRepAll, cRepInvRep);

  // include coordinate set
  if (is_new) {
    VLACheck(I->CSet, CoordSet *, frame);
    if(I->NCSet <= frame)
      I->NCSet = frame + 1;
    I->CSet[frame] = cset;
    SceneCountFrames(G);
  }

  // success
  return (I);

  // error handling
ok_except2:
  PyErr_Print();
ok_except1:
  if(is_new && cset)
    delete cset;
  ErrMessage(G, "LoadCoords", "failed");
  return NULL;
#endif
}


/*========================================================================*/
int ObjectMoleculeExtendIndices(ObjectMolecule * I, int state)
{
  int a;
  CoordSet *cs;

  if(I->DiscreteFlag && (state >= 0)) {
    /* if object is discrete, then we don't need to extend each state,
       just the current one (which updates object DiscreteAtmToIdx) */
    cs = I->CSTmpl;
    ok_assert(1, (!cs) || cs->extendIndices(I->NAtom));
    if(state < I->NCSet) {
      cs = I->CSet[state];
      ok_assert(1, (!cs) || cs->extendIndices(I->NAtom));
    }
  } else {                      /* do all states */
    for(a = -1; a < I->NCSet; a++) {
      cs = (a < 0) ? I->CSTmpl : I->CSet[a];
      ok_assert(1, (!cs) || cs->extendIndices(I->NAtom));
    }
  }
  return true;
ok_except1:
  return false;
}


/*========================================================================*/

static CoordSet *ObjectMoleculeMOLStr2CoordSet(PyMOLGlobals * G, const char *buffer,
                                               AtomInfoType ** atInfoPtr, const char **restart)
{
  const char *p;
  int nAtom, nBond;
  int a, cnt, atm, chg;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  char cc[MAXLINELEN], cc1[MAXLINELEN], resn[MAXLINELEN] = "UNK";
  float *f;
  BondType *ii;
  BondType *bond = NULL;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);
  WordType nameTmp;

  p = buffer;
  nAtom = 0;
  if(atInfoPtr)
    atInfo = *atInfoPtr;

  /*  p=ParseWordCopy(nameTmp,p,sizeof(WordType)-1); */
  p = ParseNCopy(nameTmp, p, sizeof(WordType) - 1);
  PRINTFB(G, FB_ObjectMolecule, FB_Blather)
    " ObjMolMOLStr2CoordSet: title '%s'\n", nameTmp ENDFB(G)
    p = nextline(p);
  p = nextline(p);
  p = nextline(p);

  if(ok) {
    p = ncopy(cc, p, 3);
    if(sscanf(cc, "%d", &nAtom) != 1)
      ok = ErrMessage(G, "ReadMOLFile", "bad atom count");
  }

  if(ok) {
    p = ncopy(cc, p, 3);
    if(sscanf(cc, "%d", &nBond) != 1)
      ok = ErrMessage(G, "ReadMOLFile", "bad bond count");
  }

  if(ok) {
    coord = VLAlloc(float, 3 * nAtom);
    if(atInfo)
      VLACheck(atInfo, AtomInfoType, nAtom);
  }

  p = nextline(p);

  /* read coordinates and atom names */

  if(ok) {
    f = coord;
    for(a = 0; a < nAtom; a++) {
      if(ok) {
        p = ncopy(cc, p, 10);
        if(sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad coordinate");
      }
      if(ok) {
        p = ncopy(cc, p, 10);
        if(sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad coordinate");
      }
      if(ok) {
        p = ncopy(cc, p, 10);
        if(sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad coordinate");
      }
      if(ok) {
        p = nskip(p, 1);
        p = ntrim(cc, p, 3);
        strncpy(atInfo[a].elem, cc, cElemNameLen);
        atInfo[a].name = LexIdx(G, cc);
        atInfo[a].visRep = auto_show;
      }
      if(ok) {
        int tmp_int;
        p = nskip(p, 2);
        p = ncopy(cc, p, 3);
        if(sscanf(cc, "%hhi", &atInfo[a].formalCharge) == 1) {
          if(atInfo[a].formalCharge) {
            atInfo[a].formalCharge = 4 - atInfo[a].formalCharge;
          }
        }
        p = ncopy(cc, p, 3);
        if(sscanf(cc, "%d", &tmp_int) != 1)
          atInfo[a].stereo = 0;
        else
          atInfo[a].stereo = tmp_int;
      }
      if(ok && atInfo) {
        atInfo[a].id = a + 1;
        atInfo[a].rank = a;
        atInfo[a].resn = LexIdx(G, resn);
        atInfo[a].hetatm = true;
        AtomInfoAssignParameters(G, atInfo + a);
        AtomInfoAssignColors(G, atInfo + a);
        atInfo[a].alt[0] = 0;
        atInfo[a].segi = 0;
        atInfo[a].inscode = 0;
      }
      p = nextline(p);
      if(!ok)
        break;
    }
  }
  if(ok) {
    bond = VLACalloc(BondType, nBond);
    ii = bond;
    for(a = 0; a < nBond; a++) {
      if(ok) {
        p = ncopy(cc, p, 3);
        if(sscanf(cc, "%d", &ii->index[0]) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad bond atom");
      }

      if(ok) {
        p = ncopy(cc, p, 3);
        if(sscanf(cc, "%d", &ii->index[1]) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad bond atom");
      }

      if(ok) {
        int order = 0;
        p = ncopy(cc, p, 3);
        if(sscanf(cc, "%d", &order) != 1)
          ok = ErrMessage(G, "ReadMOLFile", "bad bond order");
        ii->order = order;
      }
      if(ok) {
        p = ncopy(cc, p, 3);
        // stereo
      }
      ii++;
      if(!ok)
        break;
      p = nextline(p);
    }
    ii = bond;
    for(a = 0; a < nBond; a++) {
      ii->index[0]--;           /* adjust bond indexs down one */
      ii->index[1]--;
      ii++;
    }
  }
  while(*p) {                   /* read M  CHG records */
    auto p_line = p;
    p = ncopy(cc, p, 6);
    if(cc[5] == ' ')
      cc[5] = 0;
    if((!strcmp(cc, "M  END")) || (!strcmp(cc, "M END"))) {
      /* denotes end of MOL block */
      break;
    }
    if((!strcmp(cc, "M  CHG")) || (!strcmp(cc, "M CHG"))) {
      p = ncopy(cc, p, 3);
      if(sscanf(cc, "%d", &cnt) == 1) {
        while(cnt--) {
          p = ncopy(cc, p, 4);
          p = ncopy(cc1, p, 4);
          if(!((*cc) || (*cc1)))
            break;
          if((sscanf(cc, "%d", &atm) == 1) && (sscanf(cc1, "%d", &chg) == 1)) {
            atm--;
            if((atm >= 0) && (atm < nAtom))
              atInfo[atm].formalCharge = chg;
          }
        }
      }
    } else if (!strcmp(cc, "M  V30")) {
      p = MOLV3000Parse(G, p_line, atInfo, bond, coord, nAtom, nBond);
      if (!p) {
        ok = false;
        break;
      }
      continue;
    }
    p = nextline(p);
  }
  if(ok) {
    (*restart) = p;
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);
    strcpy(cset->Name, nameTmp);
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
    (*restart) = NULL;
  }
  if(atInfoPtr)
    *atInfoPtr = atInfo;
  return (cset);
}


/*========================================================================*/

static CoordSet *ObjectMoleculeSDF2Str2CoordSet(PyMOLGlobals * G, const char *buffer,
                                                AtomInfoType ** atInfoPtr,
                                                const char **next_mol)
{
  char cc[MAXLINELEN];
  const char *p;
  CoordSet *result = NULL;
  result = ObjectMoleculeMOLStr2CoordSet(G, buffer, atInfoPtr, next_mol);
  p = *next_mol;
  if(p) {
    while(*p) {                 /* we simply need to skip until we've read past the end of the SDF record */
      p = ncopy(cc, p, 4);
      p = nextline(p);
      if(!strcmp(cc, "$$$$")) { /* SDF record separator... */
        break;
      }
    }
    if(!*p)
      p = NULL;
  }
  *next_mol = p;
  return result;
}

/**
 * Get the sum of bond orders for every atom. For this purpose, aromatic
 * bonds will be considered either single or double bond, based on
 * the entire system (heuristic method, room for improvement).
 */
static
std::vector<signed char> get_bond_order_sums(ObjectMolecule * obj) {
  std::vector<signed char> valences(obj->NAtom);
  std::vector<signed char> freevalences(obj->NAtom);
  std::vector<signed char> orders(obj->NBond);

  // bond order sums as if all aromatic bonds were single bonds
  for (int b = 0; b < obj->NBond; ++b) {
    auto bond = obj->Bond + b;
    auto order = bond->order;
    orders[b] = order;
    if (order > 3)
      order = 1;
    valences[bond->index[0]] += order;
    valences[bond->index[1]] += order;
  }

  // determine free valences as if all aromatic bonds were single bonds
  for (int atm = 0; atm < obj->NAtom; ++atm) {
    int tmp = 0;
    switch (obj->AtomInfo[atm].protons) {
      case cAN_C:             tmp = 4; break;
      case cAN_N: case cAN_P: tmp = 5; break;
      case cAN_O: case cAN_S: tmp = 6; break;
    }
    if (tmp) {
      tmp -= valences[atm];
      if (tmp) {
        freevalences[atm] = (tmp - 1) % 2 + 1;
      }
    }
  }

  // (This is a heuristic, gets most of the ZINC dataset correct)
  // Do two passes over all aromatic bonds.
  // 1st pass: assign double bonds between atoms with `freevalence == 1`
  // 2nd pass: assign double bonds between atoms with `freevalence >= 1`
  for (int secondpass = 0; secondpass < 2; ++secondpass) {
    for (int b = 0; b < obj->NBond; ++b) {
      if (orders[b] == 4) {
        auto atm1 = obj->Bond[b].index[0];
        auto atm2 = obj->Bond[b].index[1];
        if (secondpass ?
            (freevalences[atm1] >= 1 && freevalences[atm2] >= 1 ) :
            (freevalences[atm1] == 1 && freevalences[atm2] == 1)) {
          freevalences[atm1] = 0;
          freevalences[atm2] = 0;
          valences[atm1] += 1;
          valences[atm2] += 1;
          orders[b] = 2;
        }
      }
    }
  }

  return valences;
}

/**
 * For each atom, set `formalCharge` based on mol2 `textType`
 *
 * See also: getMOL2Type() in layer2/Mol2Typing.cpp
 */
static void ObjectMoleculeMOL2SetFormalCharges(PyMOLGlobals *G, ObjectMolecule *obj){

  // check if structure has explicit hydrogens
  bool has_hydrogens = false;
  for (int at = 0; at < obj->NAtom; ++at) {
    auto ai = obj->AtomInfo + at;
    if (ai->isHydrogen()) {
      has_hydrogens = true;
      break;
    }
  }

  if (!has_hydrogens) {
    // could eventually do PDB nomenclature charge assignment here...
    return;
  }

  // Unfortunately, aromatic bonds (type 4) can't be used to determine
  // formal charges. The following uses a heuristic to guess bond orders
  // for aromatic bonds.
  auto valences = get_bond_order_sums(obj);

  // (period currently incompatible with G->lex_const)
  auto lex_N_4 = LexBorrow(G, "N.4");

  for (int at = 0; at < obj->NAtom; ++at) {
    int fcharge = 0;
    auto ai = obj->AtomInfo + at;

    if (ai->protons == cAN_N) {
      if (ai->textType == lex_N_4) {
        fcharge = 1;
      } else if (valences[at] == 2) {
        fcharge = -1;
      } else if (valences[at] == 4) {
        fcharge = 1;
      }
    } else if (ai->protons == cAN_O) {
      if (valences[at] == 1) {
        fcharge = -1;
      }
    }

    ai->formalCharge = fcharge;
  }
}

/*========================================================================*/

static CoordSet *ObjectMoleculeMOL2Str2CoordSet(PyMOLGlobals * G,
                                                const char *buffer,
                                                AtomInfoType ** atInfoPtr,
                                                const char **next_mol)
{
  const char *p;
  int nAtom, nBond, nSubst, nFeat, nSet;
  int a;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL;
  char cc[MAXLINELEN];
  float *f;
  const char *last_p;
  BondType *ii;
  BondType *bond = NULL;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);
  int have_molecule = false;
  WordType nameTmp;
  bool has_non_hetatm = false;

  // Maestro writes <0> as subst_name for waters
  const lexidx_t empty_subst_name = LexIdx(G, "<0>");

  p = buffer;
  nAtom = 0;
  if(atInfoPtr)
    atInfo = *atInfoPtr;

  while((*p) && ok) {

    last_p = p;
    /* top level -- looking for an RTI, assuming p points to the beginning of a line */
    p = ParseWordCopy(cc, p, MAXLINELEN);
    if(cc[0] == '@') {
      if(WordMatchExact(G, cc, "@<TRIPOS>MOLECULE", true)
         || WordMatchExact(G, cc, "@MOLECULE", true)) {
        if(have_molecule) {
          *next_mol = last_p;
          break;                /* next record of multi-mol2 */
        }
        p = ParseNextLine(p);
        p = ParseNTrim(nameTmp, p, sizeof(WordType) - 1);       /* get mol name */
        p = ParseNextLine(p);
        p = ParseWordCopy(cc, p, MAXLINELEN);
        if(sscanf(cc, "%d", &nAtom) != 1)
          ok = ErrMessage(G, "ReadMOL2File", "bad atom count");
        else {
          coord = VLAlloc(float, 3 * nAtom);
          if(atInfo)
            VLACheck(atInfo, AtomInfoType, nAtom);

          p = ParseWordCopy(cc, p, MAXLINELEN);
          if(sscanf(cc, "%d", &nBond) != 1) {
            nBond = 0;
          }
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if(sscanf(cc, "%d", &nSubst) != 1) {
            nSubst = 0;
          }
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if(sscanf(cc, "%d", &nFeat) != 1) {
            nFeat = 0;
          }
          p = ParseWordCopy(cc, p, MAXLINELEN);
          if(sscanf(cc, "%d", &nSet) != 1) {
            nSet = 0;
          }

        }
        p = ParseNextLine(p);
        p = ParseNextLine(p);
        have_molecule = true;
      } else if(WordMatchExact(G, cc, "@<TRIPOS>ATOM", true)
                || WordMatchExact(G, cc, "@ATOM", true)) {
        if(!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@ATOM before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);
        f = coord;
        for(a = 0; a < nAtom; a++) {
          AtomInfoType *ai = atInfo + a;

          if(ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if(sscanf(cc, "%d", &ai->id) != 1) {
              ok = ErrMessage(G, "ReadMOL2File", "bad atom id");
	    }
          }
          if(ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            cc[cAtomNameLen] = 0;
            UtilCleanStr(cc);
            if(cc[0])
              LexAssign(G, ai->name, cc);
            else
              ok = ErrMessage(G, "ReadMOL2File", "bad atom name");
          }
          if(ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if(sscanf(cc, "%f", f++) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad x coordinate");
          }
          if(ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if(sscanf(cc, "%f", f++) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad y coordinate");
          }
          if(ok) {
            p = ParseWordCopy(cc, p, MAXLINELEN);
            if(sscanf(cc, "%f", f++) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad z coordinate");
          }
          if(ok) {
            OrthoLineType temp;
            p = ParseWordCopy(cc, p, OrthoLineLength - 1);
            if(sscanf(cc, "%s", temp) != 1)
              ok = ErrMessage(G, "ReadMOL2File", "bad atom type");
            else {              /* convert atom type to elem symbol */
              char *tt = temp;
              char *el = ai->elem;
              int elc = 0;
              while(*tt && ((*tt) != '.')) {
                *(el++) = *(tt++);
                elc++;
                if(elc > cElemNameLen)
                  break;
              }
              *el = 0;
              if(el[2])
                el[0] = 0;

              ai->textType = LexIdx(G, temp);
            }
          }
          if(ok) {
            char cc_resi[8];
            p = ParseWordCopy(cc_resi, p, sizeof(cc_resi) - 1);
            if(cc_resi[0]) {         /* subst_id is residue identifier */
              ai->setResi(cc_resi);
              p = ParseWordCopy(cc, p, MAXLINELEN);
              if(cc[0]) {
                if (nSubst == 0) {
                  // without substructure information (e.g. OpenBabel exported MOL2):
                  // if subst_name includes the subst_id (e.g. 5 ALA5) then strip the number
                  int len_resi = strlen(cc_resi);
                  int len_resn = strlen(cc);
                  if (len_resn > len_resi) {
                    if (strcmp(cc_resi, cc + len_resn - len_resi) == 0) {
                      cc[len_resn - len_resi] = 0;
                    }
                  }
                }

                LexAssign(G, ai->resn, cc);

                p = ParseWordCopy(cc, p, MAXLINELEN);
                if(cc[0]) {
                  if(sscanf(cc, "%f", &ai->partialCharge) != 1)
                    ok = ErrMessage(G, "ReadMOL2File", "bad atom charge");
                }

                // status_bit
                p = ParseWordCopy(cc, p, MAXLINELEN);
                if (cc[0]) {
                  // for water, substitute resn "<0>" with "HOH" for
                  // correct classification as "solvent"
                  if (ai->resn == empty_subst_name && strstr(cc, "WATER")) {
                    LexAssign(G, ai->resn, G->lex_const.HOH);
                  }
                }
              }
            }
          }
          p = ParseNextLine(p);

          ai->visRep = auto_show;

          ai->id = a + 1;
          ai->rank = a;
          ai->hetatm = true;
          AtomInfoAssignParameters(G, atInfo + a);
          AtomInfoAssignColors(G, atInfo + a);
          ai->alt[0] = 0;
          ai->segi = 0;

        }
      } else if(WordMatchExact(G, cc, "@<TRIPOS>SET", true)
                || WordMatchExact(G, cc, "@SET", true)) {
        if(!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@SET before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);
        if(ok) {
          for(a = 0; a < nSet; a++) {
            if(ok) {
              char ss = 0;
              p = ParseWordCopy(cc, p, 50);
              cc[5] = 0;
              if(!strcmp("HELIX", cc))
                ss = 'H';
              if(!strcmp("SHEET", cc))
                ss = 'S';
              p = ParseWordCopy(cc, p, 50);
              if(!strcmp("STATIC", cc)) {
                int nMember;

                p = ParseNextLine(p);
                p = ParseWordCopy(cc, p, 20);
                if(sscanf(cc, "%d", &nMember) != 1)
                  ok = ErrMessage(G, "ReadMOL2File", "bad member count");
                else {
                  while(ok && (nMember > 0)) {
                    p = ParseWordCopy(cc, p, 20);
                    if((!cc[0]) && (!*p)) {
                      ok = false;
                      break;
                    }
                    if(cc[0] != '\\') {
                      int atom_id;
                      if(sscanf(cc, "%d", &atom_id) != 1) {
                        ok = ErrMessage(G, "ReadMOL2File", "bad member");
                      } else {
                        atom_id--;
                        if((atom_id >= 0) && (atom_id < nAtom)) {
                          atInfo[atom_id].ssType[0] = ss;
                          atInfo[atom_id].ssType[1] = 0;
                        }
                        nMember--;
                      }
                    } else {
                      p = ParseNextLine(p);
                    }
                  }
                }
                p = ParseNextLine(p);
              } else {
                p = ParseNextLine(p);
                p = ParseNextLine(p);
              }
            }
          }
        }
      } else if(WordMatchExact(G, cc, "@<TRIPOS>BOND", true)
                || WordMatchExact(G, cc, "@BOND", true)) {
        if(!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@BOND before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);

        if(ok) {
          bond = VLACalloc(BondType, nBond);
          ii = bond;
          for(a = 0; a < nBond; a++) {
            if(ok) {
              p = ParseWordCopy(cc, p, 20);
            }

            if(ok) {
              p = ParseWordCopy(cc, p, 20);
              if(sscanf(cc, "%d", &ii->index[0]) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad source atom id");
            }

            if(ok) {
              p = ParseWordCopy(cc, p, 20);
              if(sscanf(cc, "%d", &ii->index[1]) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad target atom id");
            }

            if(ok) {
              p = ParseWordCopy(cc, p, 20);
              if(!cc[1]) {
                switch (cc[0]) {
                case '1':
                  ii->order = 1;
                  break;
                case '2':
                  ii->order = 2;
                  break;
                case '3':
                  ii->order = 3;
                  break;
                case '4':
                  ii->order = 4;
                  break;
                }
              } else if(WordMatchExact(G, "ar", cc, true)) {
                ii->order = 4;
              } else if(WordMatchExact(G, "am", cc, true)) {
                ii->order = 1;
              } else if(WordMatchExact(G, "un", cc, true)) {
                ii->order = 1;
              } else if(WordMatchExact(G, "nc", cc, true)) {
                ii->order = 0;  /* is this legal in PyMOL? */
              } else if(WordMatchExact(G, "du", cc, true)) {
                ii->order = 1;
              } else
                ok = ErrMessage(G, "ReadMOL2File", "bad bond type");
            }
            ii++;
            if(!ok)
              break;
            p = ParseNextLine(p);
          }
          ii = bond;
          for(a = 0; a < nBond; a++) {
            ii->index[0]--;     /* adjust bond indexs down one */
            ii->index[1]--;
            ii++;
          }
        }
      } else if(WordMatchExact(G, cc, "@<TRIPOS>SUBSTRUCTURE", true)
                || WordMatchExact(G, cc, "@SUBSTRUCTURE", true)) {
        if(!have_molecule) {
          ok = ErrMessage(G, "ReadMOL2File", "@SUBSTSRUCTURE before @MOLECULE!");
          break;
        }
        p = ParseNextLine(p);
        if(ok) {
          WordType subst_name;
          SegIdent segment;     /* what MOL2 calls chain */
          lexidx_t chain = 0;
          WordType subst_type;
          ResIdent resi;
          ResName resn;
          int id, dict_type, root, resv;
          int end_line, seg_flag, subst_flag, resi_flag;
          int chain_flag, resn_flag;
          std::unordered_map<int, std::vector<int>> resv_map;

          {
            int b;
            AtomInfoType *ai = atInfo;
            for(b = 0; b < nAtom; b++) {
              resv_map[ai->resv].push_back(b);
              ai++;
            }
          }

          for(a = 0; a < nSubst; a++) {
            bool hetatm = true;
            segment[0] = 0;
            subst_name[0] = 0;
            LexAssign(G, chain, 0);
            resn[0] = 0;
            resi[0] = 0;
            end_line = false;
            seg_flag = false;
            subst_flag = false;

            resi_flag = false;
            chain_flag = false;
            resn_flag = false;

            if(ok) {
              const char *save_p = p;
              p = ParseWordCopy(cc, p, 20);
              if(sscanf(cc, "%d", &id) != 1) {
                if(cc[0] != '@')
                  ok = ErrMessage(G, "ReadMOL2File", "bad substructure id");
                else {
                  p = save_p;
                  break;        /* missing substructure... */
                }
              }
            }
            if(ok) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              if(sscanf(cc, "%s", subst_name) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad substructure name");
            }
            if(ok) {
              p = ParseWordCopy(cc, p, 20);
              if(sscanf(cc, "%d", &root) != 1)
                ok = ErrMessage(G, "ReadMOL2File", "bad target root atom id");
              else
                root--;         /* convert 1-based to 0-based */
            }

            /* optional data */
            if(ok && (!end_line)) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              if(sscanf(cc, "%s", subst_type) != 1) {
                end_line = true;
                subst_type[0] = 0;
              } else {
                subst_flag = 1;
              }
            }

            if(ok && (!end_line)) {
              p = ParseWordCopy(cc, p, 20);
              if(sscanf(cc, "%d", &dict_type) != 1)
                end_line = true;
            }

            if(ok && (!end_line)) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              cc[cSegiLen] = 0;
              if(sscanf(cc, "%s", segment) != 1) {
                end_line = true;
                segment[0] = 0;
              } else if(strcmp(segment, "****")) {
                seg_flag = true;
                if(!segment[1]) {       /* if segment is single letter, then also assign to chain field */
                  LexAssign(G, chain, segment);
                  chain_flag = true;
                }
              }
            }

            if(ok && (!end_line)) {
              p = ParseWordCopy(cc, p, sizeof(WordType) - 1);
              cc[cResnLen] = 0;
              if(!strcmp(cc, "****")) { /* whoops -- no residue name... */
                char *pp;
                UtilNCopy(resn, subst_name, sizeof(ResName));
                pp = resn;
                /* any numbers? */
                while(*pp) {
                  if(((*pp) >= '0') && ((*pp) <= '9'))
                    break;
                  pp++;
                }
                /* trim them off */
                while(pp >= resn) {
                  if(((*pp) >= '0') && ((*pp) <= '9'))
                    *pp = 0;
                  else
                    break;
                  pp--;
                }

                if(resn[0])
                  resn_flag = true;
              } else {
                if(sscanf(cc, "%s", resn) != 1) {
                  end_line = true;
                  resn[0] = 0;
                } else {
                  resn_flag = true;
                }
              }
            }

            if(ok && (root >= 0) && (root < nAtom)) {

              if((strlen(subst_name) == 4) &&
                 ((subst_name[0] >= '1') && (subst_name[0] <= '9')) &&
                 ((((subst_name[1] >= 'A') && (subst_name[1] <= 'Z')) ||
                   ((subst_name[1] >= 'a') && (subst_name[1] <= 'z'))) ||
                  (((subst_name[2] >= 'A') && (subst_name[2] <= 'Z')) ||
                   ((subst_name[2] >= 'a') && (subst_name[2] <= 'z'))) ||
                  (((subst_name[3] >= 'A') && (subst_name[3] <= 'Z')) ||
                   ((subst_name[3] >= 'a') && (subst_name[3] <= 'z'))))) {

                /* if subst_name looks like a PDB code, then it isn't a residue identifier
                   so do nothing... */
              } else {

                if(!subst_flag) {

                  /* Merck stuffs the chain ID and the residue ID in the
                     subst_name field, so handle that case first */

                  /* get the residue value */

                  ParseIntCopy(cc, subst_name, 20);
                  if(sscanf(cc, "%d", &resv) == 1) {
                    /* we have the resv, so now get the chain if there is one */
                    char *pp = subst_name;
                    if(!((pp[0] >= '0') && (pp[0] <= '9'))) {
                      char tmp[2] = {pp[0], 0};
                      LexAssign(G, chain, tmp);
                      chain_flag = true;
                      pp++;
                    }
                    UtilNCopy(resi, pp, cResiLen);      /* now get the textual residue identifier */

                    if(resi[0]) {
                      resi_flag = true;
                    }
                  }
                } else {
                  /* normal mode: assume that the substructure ID is a vanilla residue identifier */

                  char *pp = subst_name;

                  if(resn_flag) {       /* if we have a resn, then remove it from the subst_name */
                    char *qq = resn;
                    while((*pp) && (*qq)) {
                      if(*pp == *qq) {
                        pp++;
                        qq++;
                      } else
                        break;
                    }
                  }

                  ParseIntCopy(cc, pp, 20);
                  if(sscanf(cc, "%d", &resv) == 1) {
                    UtilNCopy(resi, pp, cResiLen);      /* now get the textual residue identifier */
                    if(resi[0]) {
                      resi_flag = true;
                    }
                  }

                  if (strcmp(subst_type, "RESIDUE") == 0) {
                    hetatm = false;
                    has_non_hetatm = true;
                  }
                }
              }

              /* for now, assume that atom ids are 1-based and sequential */

              if(ok) {
                if(resi_flag || chain_flag || resn_flag || seg_flag) {
                  auto it = resv_map.find(id);
                  if (it != resv_map.end()) {
                    for (auto b : it->second) {
                      assert(b >= 0 && b < nAtom);
                  {
                        AtomInfoType *ai = atInfo + b;
                        if(resi_flag)
                          ai->setResi(resi);
                        if(chain_flag) {
                          LexAssign(G, ai->chain, chain);
                        }
                        if(resn_flag)
                          LexAssign(G, ai->resn, resn);
                        if(seg_flag)
                          LexAssign(G, ai->segi, segment);
                        ai->hetatm = hetatm;
                      }
                    }
                  }
                }
              }
            }
            if(!ok)
              break;
            p = ParseNextLine(p);
          }
          LexDec(G, chain);
        }
      } else
        p = ParseNextLine(p);
    } else
      p = ParseNextLine(p);
  }

  LexDec(G, empty_subst_name);

  if (has_non_hetatm) {
    // contains "residues" (assume polymer) so ignore hetatm for surfacing
    for (auto ai = atInfo, ai_end = atInfo + nAtom; ai != ai_end; ++ai) {
      if (ai->hetatm) {
        ai->flags |= cAtomFlag_ignore;
      }
    }
  }

  if(ok) {
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);
    strcpy(cset->Name, nameTmp);
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
  }
  if(atInfoPtr)
    *atInfoPtr = atInfo;
  return (cset);
}

/**
 * Read one molecule in MOL2, MOL, SDF, MMD or XYZ format from the string pointed
 * to by (*next_entry). All these formats (except MOL) support multiple
 * concatenated entries in one file. If multiplex=1, then read only one
 * molecule (one state) and set the next_entry pointer to the beginning of the
 * next entry. Otherwise, read a multi-state molecule.
 */
ObjectMolecule *ObjectMoleculeReadStr(PyMOLGlobals * G, ObjectMolecule * I,
                                      const char **next_entry,
                                      cLoadType_t content_format, int frame,
                                      int discrete, int quiet, int multiplex,
                                      char *new_name,
				      short loadpropertiesall, OVLexicon *loadproplex)
{
  int ok = true;
  CoordSet *cset = NULL;
  pymol::vla<AtomInfoType> atInfo;
  int isNew;
  int nAtom;
  const char *restart = NULL, *start = *next_entry;
  int repeatFlag = true;
  int successCnt = 0;
  char tmpName[WordLength];
  int deferred_tasks = false;
  int skip_out;
  int connect = false;
  int set_formal_charges = false;
  *next_entry = NULL;
  int aic_mask = cAIC_MOLMask;

  while(repeatFlag) {
    repeatFlag = false;
    skip_out = false;

    if(!I)
      isNew = true;
    else
      isNew = false;

    if(isNew) {
      I = (ObjectMolecule *) new ObjectMolecule(G, (discrete > 0));
      std::swap(atInfo, I->AtomInfo);
    } else {
      atInfo = pymol::vla<AtomInfoType>(10);
    }

    if(isNew) {
      I->Color = AtomInfoUpdateAutoColor(G);
    }

    restart = NULL;
    switch (content_format) {
    case cLoadTypeMOL2:
    case cLoadTypeMOL2Str:
      cset = ObjectMoleculeMOL2Str2CoordSet(G, start, &atInfo, &restart);
      if (cset){
	set_formal_charges = true;
      }
      break;
    case cLoadTypeMOL:
    case cLoadTypeMOLStr:
      cset = ObjectMoleculeMOLStr2CoordSet(G, start, &atInfo, &restart);
      restart = NULL;
      break;
    case cLoadTypeSDF2:
    case cLoadTypeSDF2Str:
      cset = ObjectMoleculeSDF2Str2CoordSet(G, start, &atInfo, &restart);
      break;
    case cLoadTypeXYZ:
    case cLoadTypeXYZStr:
      cset = ObjectMoleculeXYZStr2CoordSet(G, start, &atInfo, &restart);
      if(!cset->TmpBond)
        connect = true;
      break;
    case cLoadTypeMMD:
    case cLoadTypeMMDStr:
      cset = ObjectMoleculeMMDStr2CoordSet(G, start, &atInfo, &restart);
      aic_mask = cAIC_MMDMask;
      break;
    }


    if(!cset) {
      if(!isNew)
        VLAFreeP(atInfo);
      if(!successCnt) {
	if (isNew)
          std::swap(I->AtomInfo, atInfo);
        DeleteP(I);
        ok = false;
      } else {
        skip_out = true;
      }
    }

    if(ok && !skip_out) {

      if((discrete>0 && !I->DiscreteFlag) ||  // should set to discrete if explicitly defined
	 ((discrete < 0) && (restart) && isNew && (multiplex <= 0))) {
        /* if default discrete behavior and new object, with
           multi-coordinate set file, and not multiplex, then set
           discrete */

        ObjectMoleculeSetDiscrete(G, I, true);
      }

      if(frame < 0)
        frame = I->NCSet;
      if(I->NCSet <= frame)
        I->NCSet = frame + 1;
      VLACheck(I->CSet, CoordSet *, frame);

      nAtom = cset->NIndex;

      if(I->DiscreteFlag && atInfo) {
        int a;
        int fp1 = frame + 1;
        AtomInfoType *ai = atInfo.data();
        for(a = 0; a < nAtom; a++) {
          (ai++)->discrete_state = fp1;
        }
      }

      if(multiplex > 0)
        UtilNCopy(tmpName, cset->Name, WordLength);

      cset->Obj = I;
      cset->enumIndices();
      cset->invalidateRep(cRepAll, cRepInvRep);
      if(isNew) {
        std::swap(I->AtomInfo, atInfo);
      } else {
        ObjectMoleculeMerge(I, std::move(atInfo), cset, false, aic_mask, false);
        /* NOTE: will release atInfo */
      }

      if(isNew)
        I->NAtom = nAtom;
      if(frame < 0)
        frame = I->NCSet;
      VLACheck(I->CSet, CoordSet *, frame);
      if(I->NCSet <= frame)
        I->NCSet = frame + 1;
      delete I->CSet[frame];
      I->CSet[frame] = cset;

      if(ok && isNew)
        ok &= ObjectMoleculeConnect(I, cset, connect);

      if (ok)
	ok &= ObjectMoleculeExtendIndices(I, frame);
      if (ok)
	ok &= ObjectMoleculeSort(I);

      deferred_tasks = true;
      successCnt++;
      if(!quiet) {
        if(successCnt > 1) {
          if(successCnt == 2) {
            PRINTFB(G, FB_ObjectMolecule, FB_Actions)
              " %s: read through molecule %d.\n", __func__, 1 ENDFB(G);
          }
          if(UtilShouldWePrintQuantity(successCnt)) {
            PRINTFB(G, FB_ObjectMolecule, FB_Actions)
              " %s: read through molecule %d.\n", __func__, successCnt ENDFB(G);
          }
        }
      }
      if(multiplex > 0) {
        UtilNCopy(new_name, tmpName, WordLength);
        if(restart) {
          *next_entry = restart;
        }
      } else if(restart) {
        repeatFlag = true;
        start = restart;
        frame = frame + 1;
      }
    }
  }
  if(deferred_tasks && I) {     /* defer time-consuming tasks until all states have been loaded */
    if (set_formal_charges){
      ObjectMoleculeMOL2SetFormalCharges(G, I);
    }
    SceneCountFrames(G);
    I->invalidate(cRepAll, cRepInvAll, -1);
    ObjectMoleculeUpdateIDNumbers(I);
    ObjectMoleculeUpdateNonbonded(I);
  }
  return (I);
}


/*========================================================================*/
typedef int CompareFn(PyMOLGlobals *, const AtomInfoType *, const AtomInfoType *);
int ObjectMoleculeMerge(ObjectMolecule * I, pymol::vla<AtomInfoType>&& ai,
			CoordSet * cs, int bondSearchFlag, int aic_mask, int invalidate)
{
  PyMOLGlobals *G = I->G;
  int *index, *outdex;
  int a, b, lb = 0, ac;
  int c, nb, a1, a2;
  int found;
  int nAt, nBd, nBond;
  int expansionFlag = false;
  int ok = true;

  const auto n_index = cs->getNIndex();
  const auto oldNAtom = I->NAtom;
  const auto oldNBond = I->NBond;

  /* first, sort the coodinate set */

  index = AtomInfoGetSortedIndex(G, I, ai, n_index, &outdex);
  CHECKOK(ok, index);
  if (!ok)
    return false;

  auto ai2 = pymol::vla<AtomInfoType>(n_index);
  if (ok){
    for(a = 0; a < n_index; a++)
      ai2[a] = std::move(ai[index[a]]);      /* creates a sorted list of atom info records */
  }
  ai = std::move(ai2);

  /* now, match it up with the current object's atomic information */

  if (ok){
    for(a = 0; a < n_index; a++) {
      index[a] = -1;
    }
  }

  if (ok) {
    const auto n_atom = I->NAtom;
    const AtomInfoType* const atInfo = I->AtomInfo.data();
    CompareFn *fCompare;

    if(SettingGetGlobal_b(G, cSetting_pdb_hetatm_sort)) {
      fCompare = AtomInfoCompareIgnoreRank;
    } else {
      fCompare = AtomInfoCompareIgnoreRankHet;
    }

    c = 0;
    b = 0;
    lb = 0;

    for(a = 0; a < n_index; a++) {
      int reverse = false;
      const auto* const ai_a = ai + a;
      found = false;
      if(!I->DiscreteFlag) {    /* don't even try matching for discrete objects */
        while(b < n_atom) {
          ac = (fCompare(G, ai_a, atInfo + b));
          if(!ac) {
            found = true;
            break;
          } else if(ac < 0) {   /* atom is before current, so try going the other way */
            reverse = true;
            break;
          } else if(!b) {
            int idx;
            ac = (fCompare(G, ai_a, atInfo + n_atom - 1));
            if(ac > 0) {        /* atom is after all atoms in list, so don't even bother searching */
              break;
            }
            /* next, try to jump to an appropriate position in the list */
            idx = ((a * n_atom) / n_index) - 1;
            if((idx > 0) && (b != idx) && (idx < n_atom)) {
              ac = (fCompare(G, ai_a, atInfo + idx));
              if(ac > 0)
                b = idx;
            }
          }
          lb = b;               /* last b before atom */
          b++;
        }
        if(reverse && !found) { /* searching going backwards... */
          while(b >= 0) {
            ac = (fCompare(G, ai_a, atInfo + b));
            if(!ac) {
              found = true;
              break;
            } else if(ac > 0) { /* atom after current -- no good */
              break;
            }
            lb = b;
            b--;
          }
          if(b < 0)
            b = 0;
        }
      }
      if(found) {
        index[a] = b;           /* store real atom index b for a in index[a] */
        b++;
      } else {
        index[a] = I->NAtom + c;        /* otherwise, this is a new atom */
        c++;
        b = lb;
      }
    }
  }

  /* first, reassign atom info for matched atoms */

  /* allocate additional space */
  if (ok){
    if(c) {
      expansionFlag = true;
      nAt = I->NAtom + c;
    } else {
      nAt = I->NAtom;
    }
    if(expansionFlag) {
      VLACheck(I->AtomInfo, AtomInfoType, nAt);
      CHECKOK(ok, I->AtomInfo);
    }
  }

  if (ok){
    for(a = 0; a < n_index; a++) {     /* a is in original file space */
      a1 = outdex[cs->IdxToAtm[a]];    /* a1 is in sorted atom info space */
      a2 = index[a1];
      cs->IdxToAtm[a] = a2;       /* a2 is in object space */
      if(a2 < oldNAtom)
        AtomInfoCombine(G, I->AtomInfo + a2, std::move(ai[a1]), aic_mask);
      else
	I->AtomInfo[a2] = std::move(ai[a1]);
    }
  }

  if(ok && I->DiscreteFlag) {
    ok = I->setNDiscrete(nAt);
  }

  I->NAtom = nAt;

  if (ok){
    cs->updateNonDiscreteAtmToIdx(nAt);
  }

  VLAFreeP(ai);                 /* note that we're trusting AtomInfoCombine to have 
                                   already purged all atom-dependent storage (such as strings) */

  AtomInfoFreeSortedIndexes(G, &index, &outdex);

  /* now find and integrate and any new bonds */
  if(ok && expansionFlag) {           /* expansion flag means we have introduced at least 1 new atom */
    pymol::vla<BondType> bond;
    ok &= ObjectMoleculeConnect(I, nBond, bond, cs, bondSearchFlag, -1);
    if(nBond) {
      index = pymol::malloc<int>(nBond);
      CHECKOK(ok, index);
      c = 0;
      b = 0;
      nb = 0;
      if (ok){
	for(a = 0; a < nBond; a++) {      /* iterate over new bonds */
	  found = false;
	  if(!I->DiscreteFlag) {  /* don't even try matching for discrete objects */
	    b = nb;               /* pick up where we left off */
	    while(b < I->NBond) {
	      ac = BondCompare(bond + a, I->Bond + b);
	      if(!ac) {           /* zero is a match */
		found = true;
		break;
	      } else if(ac < 0) { /* gone past position of this bond */
		break;
	      } else if(!b) {
		int idx;
		ac = BondCompare(bond + a, I->Bond + I->NBond - 1);
		
		if(ac > 0)        /* bond is after all bonds in list, so don't even bother searching */
		  break;
		/* next, try to jump to an appropriate position in the list */
		idx = ((a * I->NBond) / nBond) - 1;
		if((idx > 0) && (b != idx) && (idx < I->NBond)) {
		  ac = BondCompare(bond + a, I->Bond + idx);
		  if(ac > 0)
		    b = idx;
		}
	      }
	      b++;                /* no match yet, keep looking */
	    }
	  }
	  if(found) {
	    index[a] = b;         /* existing bond... */
	    nb = b + 1;
	  } else {                /* this is a new bond, save index and increment */
	    index[a] = I->NBond + c;
	    c++;
	  }
	}
      }
      /* first, reassign atom info for matched atoms */
      if(c) {
        /* allocate additional space */
        nBd = I->NBond + c;

        if(!I->Bond)
          I->Bond = pymol::vla<BondType>(1);
	CHECKOK(ok, I->Bond);
	if (ok)
	  VLACheck(I->Bond, BondType, nBd);
	CHECKOK(ok, I->Bond);
        for(a = 0; a < nBond; a++) {    /* copy the new bonds */
          a2 = index[a];
          if(a2 >= I->NBond) {
            I->Bond[a2] = bond[a];
          }
        }
        I->NBond = nBd;
      }
      FreeP(index);
    }
  }
  if(invalidate) {
    if(oldNAtom) {
      if(oldNAtom == I->NAtom) {
        if(oldNBond != I->NBond) {
          I->invalidate(cRepAll, cRepInvBonds, -1);
        }
      } else {
        I->invalidate(cRepAll, cRepInvAtoms, -1);
      }
    }
  }
  return ok;
}


/*========================================================================*/
/**
 * @param state Object state or -2 for current state
 * @return NULL if there is no CoordSet for the given state
 */
CoordSet* ObjectMolecule::getCoordSet(int state)
{
  return static_cast<CoordSet*>(getObjectState(state));
}
const CoordSet* ObjectMolecule::getCoordSet(int state) const
{
  return static_cast<const CoordSet*>(getObjectState(state));
}


/*========================================================================*/
void ObjectMoleculeTransformTTTf(ObjectMolecule * I, float *ttt, int frame)
{
  int b;
  CoordSet *cs;
  for(b = 0; b < I->NCSet; b++) {
    if((frame < 0) || (frame == b)) {
      cs = I->CSet[b];
      if(cs) {
        cs->invalidateRep(cRepAll, cRepInvCoord);
        MatrixTransformTTTfN3f(cs->NIndex, cs->Coord.data(), ttt, cs->Coord.data());
        CoordSetRecordTxfApplied(cs, ttt, false);
      }
    }
  }
}


/*========================================================================*/
bool ObjectMoleculeSeleOp(ObjectMolecule * I, int sele, ObjectMoleculeOpRec * op)
{
  float *coord;
  int a, b, s;
  int c, t_i;
  int a1 = 0, ind;
  float rms;
  float v1[3], v2, *vv1, *vv2, *vt, *vt1, *vt2;
  int hit_flag = false;
  int ok = true;
  int cnt;
  int skip_flag;
  int match_flag = false;
  int offset;
  int priority;
  int use_matrices = false;
  CoordSet *cs;
  AtomInfoType *ai, *ai0;
  PyMOLGlobals *G = I->G;
#ifndef _PYMOL_NOPY
  PyObject* expr_co = nullptr;
  int compileType = Py_single_input;
#endif
#ifdef _WEBGL
#endif
  PRINTFD(G, FB_ObjectMolecule)
    " %s-DEBUG: sele %d op->code %d\n", __func__, sele, op->code ENDFD;
  if(sele >= 0) {
    const char *errstr = "Alter";
    /* always run on entry */
    switch (op->code) {
    case OMOP_LABL:
      errstr = "Label";
      if (op->i2 != cExecutiveLabelEvalOn){
	break;
      }
#ifndef _PYMOL_NOPY
      compileType = Py_eval_input;
#endif
    case OMOP_ALTR:
    case OMOP_AlterState:
      // assume blocked interpreter

      if (op->s1 && op->s1[0]){
#ifndef _PYMOL_NOPY
	expr_co = Py_CompileString(op->s1, "", compileType);
	if(expr_co == NULL) {
	  ok = ErrMessage(G, errstr, "failed to compile expression");
	}
#else
#ifndef _WEBGL
        ok = ErrMessage(G, errstr, "failed to compile expression");
#endif
#endif
      }
      /* PBlockAndUnlockAPI() is not safe.
       * what if "v" is invalidated by another thread? */
      break;
    }
    switch (op->code) {
    case OMOP_ReferenceStore:
    case OMOP_ReferenceRecall:
    case OMOP_ReferenceValidate:
    case OMOP_ReferenceSwap:
      for(b = 0; b < I->NCSet; b++) {
        if(I->CSet[b]) {
          if((b == op->i1) || (op->i1 < 0)) {
            int inv_flag = false;
            cs = I->CSet[b];
            if(cs && CoordSetValidateRefPos(cs)) {
              for(a = 0; a < I->NAtom; a++) {
                s = I->AtomInfo[a].selEntry;
                if(SelectorIsMember(G, s, sele)) {
                  ind = cs->atmToIdx(a);
                  if(ind >= 0) {
                    float *v = cs->coordPtr(ind);
                    RefPosType *rp = cs->RefPos + ind;
                    switch (op->code) {
                    case OMOP_ReferenceStore:
                      copy3f(v, rp->coord);
                      rp->specified = true;
                      break;
                    case OMOP_ReferenceRecall:
                      if(rp->specified) {
                        copy3f(rp->coord, v);
                        inv_flag = true;
                      }
                      break;
                    case OMOP_ReferenceValidate:
                      if(!rp->specified) {
                        copy3f(v, rp->coord);
                        rp->specified = true;
                      }
                      break;
                    case OMOP_ReferenceSwap:
                      if(rp->specified) {
                        copy3f(rp->coord, v1);
                        copy3f(v, rp->coord);
                        copy3f(v1, v);
                        inv_flag = true;
                      }
                      break;
                    }
                    op->i2++;
                  }
                }
              }
            }
            if(inv_flag && cs) {
              cs->invalidateRep(cRepAll, cRepInvRep);
            }
          }
        }
      }
      break;
    case OMOP_AddHydrogens:
      if (ok) {
          ok &= ObjectMoleculeAddSeleHydrogensRefactored(I, sele, op->i1);
      }
      break;
    case OMOP_FixHydrogens:
      if (ok)
	ok &= ObjectMoleculeFixSeleHydrogens(I, sele, -1);      /* state? */
      break;
    case OMOP_RevalenceFromSource:
    case OMOP_RevalenceByGuessing:
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        if(SelectorIsMember(G, ai->selEntry, sele)) {
          hit_flag = true;
          break;
        }
        ai++;
      }
      break;
    case OMOP_PrepareFromTemplate:
      ai0 = op->ai;             /* template atom */
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          ai = I->AtomInfo + a;
          if(op->i1 != 3) {
            ai->hetatm = ai0->hetatm;
            ai->flags = ai0->flags;
            LexAssign(G, ai->chain, ai0->chain);
            strcpy(ai->alt, ai0->alt);
            LexAssign(G, ai->segi, ai0->segi);
          }
          if(op->i1 == 1) {     /* mode 1, merge residue information */
            ai->inscode = ai0->inscode;
            ai->resv = ai0->resv;
            LexAssign(G, ai->resn, ai0->resn);
          }
          AtomInfoAssignColors(G, ai);
          if(op->i3) {
            if((ai->elem[0] == ai0->elem[0]) && (ai->elem[1] == ai0->elem[1]))
              ai->color = ai0->color;
            else if(ai->protons == cAN_C) {
              ai->color = op->i4;
            }
          }
          ai->visRep = ai0->visRep;
          ai->id = -1;
          ai->rank = -1;
          op->i2++;
        }
      }
      break;
    case OMOP_Sort:
      for(a = 0; ok && a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          if (ok)
	    ok &= ObjectMoleculeSort(I);
          break;
        }
      }
      break;
    case OMOP_SetAtomicSetting:
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        if(SelectorIsMember(G, ai->selEntry, sele)) {
          int uid = AtomInfoCheckUniqueID(G, ai);
          ai->has_setting = true;
          if (SettingUniqueSetTypedValue(G, uid, op->i1, op->i2, op->ii1))
            op->i4++;
        }
        ai++;
      }
      break;
    case OMOP_Pop:
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          if(SelectorMoveMember(G, s, sele, op->i1)) {
            op->i2--;
            op->i3++;
          }
          if(!op->i2)
            break;
        }
      }
      break;
    case OMOP_AVRT:            /* average vertex coordinate */
    case OMOP_StateVRT:        /* state vertex coordinate */
      {
        int op_i2 = op->i2;
        int obj_TTTFlag = I->TTTFlag;
        int b_end = I->NCSet;
        if (op->code == OMOP_StateVRT && op->i1 < b_end) {
          b_end = op->i1 + 1;
        }
        if(op_i2) {
          use_matrices =
            SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);
          if(use_matrices<0) use_matrices = 0;
        }
        for(a = 0; a < I->NAtom; a++) {
          s = I->AtomInfo[a].selEntry;
          if(!(priority = SelectorIsMember(G, s, sele)))
            continue;

          cnt = 0;

          // all states for AVRT, one state for StateVRT (don't use
          // StateIterator which depends on settings)
          for(b = op->i1; b < b_end; ++b) {
            if(!(cs = I->CSet[b]))
              continue;

            if((a1 = cs->atmToIdx(a)) == -1)
              continue;

            if(!cnt) {
              VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
            }
            cnt++;
            vv2 = cs->coordPtr(a1);

            if(op_i2) {   /* do we want transformed coordinates? */
              if(use_matrices) {
                if(!cs->Matrix.empty()) {    /* state transformation */
                  transform44d3f(cs->Matrix.data(), vv2, v1);
                  vv2 = v1;
                }
              }
              if(obj_TTTFlag) {
                transformTTT44f3f(I->TTT, vv2, v1);
                vv2 = v1;
              }
            }

            // sum up coordinates over states
            vv1 = op->vv1 + (op->nvv1 * 3);
            add3f(vv1, vv2, vv1);
          }

          // number of summed up coordinates (states) for this atom
          VLACheck(op->vc1, int, op->nvv1);
          op->vc1[op->nvv1] = cnt;

          // ordered_selections
          if(op->vp1) {
            VLACheck(op->vp1, int, op->nvv1);
            op->vp1[op->nvv1] = priority;
          }

          // atom pointer VLA
          if(op->ai1VLA) {
            VLACheck(op->ai1VLA, AtomInfoType *, op->nvv1);
            op->ai1VLA[op->nvv1] = I->AtomInfo + a;
            I->AtomInfo[a].temp1 = a;
            /* KLUDGE ALERT!!! storing atom index in the temp1 field... */
          }

          // number of atoms in selection (incl. the ones with no coordinates)
          op->nvv1++;
        }
      }
      break;
    case OMOP_SFIT:            /* state fitting within a single object */
      vt = pymol::malloc<float>(3 * op->nvv2);  /* temporary (matching) target vertex pointers */
      cnt = 0;
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          cnt++;
          break;
        }
      }
      if(cnt) {                 /* only perform action for selected object */

        for(b = 0; b < I->NCSet; b++) {
          rms = -1.0;
          vt1 = vt;             /* reset target vertex pointers */
          vt2 = op->vv2;
          t_i = 0;              /* original target vertex index */
          if(I->CSet[b] && (b != op->i2)) {
            op->nvv1 = 0;
            for(a = 0; a < I->NAtom; a++) {
              s = I->AtomInfo[a].selEntry;
              if(SelectorIsMember(G, s, sele)) {
                a1 = I->CSet[b]->atmToIdx(a);
                if(a1 >= 0) {

                  match_flag = false;
                  while(t_i < op->nvv2) {
                    if(op->i1VLA[t_i] == a) {   /* same atom? */
                      match_flag = true;
                      break;
                    }
                    if(op->i1VLA[t_i] < a) {    /* catch up? */
                      t_i++;
                      vt2 += 3;
                    } else
                      break;
                  }
                  if(match_flag) {
                    VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                    vv2 = I->CSet[b]->coordPtr(a1);
                    vv1 = op->vv1 + (op->nvv1 * 3);
                    *(vv1++) = *(vv2++);
                    *(vv1++) = *(vv2++);
                    *(vv1++) = *(vv2++);
                    *(vt1++) = *(vt2);
                    *(vt1++) = *(vt2 + 1);
                    *(vt1++) = *(vt2 + 2);
                    op->nvv1++;
                  }
                }
              }
            }
            if(op->nvv1 != op->nvv2) {
              PRINTFB(G, FB_Executive, FB_Warnings)
                "Executive-Warning: Missing atoms in state %d (%d instead of %d).\n",
                b + 1, op->nvv1, op->nvv2 ENDFB(G);
            }
            if(op->nvv1) {
              if(op->i1 != 0)   /* fitting flag */
                rms = MatrixFitRMSTTTf(G, op->nvv1, op->vv1, vt, NULL, op->ttt);
              else
                rms = MatrixGetRMS(G, op->nvv1, op->vv1, vt, NULL);
              if(op->i1 == 2) {
                ObjectMoleculeTransformTTTf(I, op->ttt, b);

                if(op->i3) {
                  const float divisor = (float) op->i3;
                  const float premult = (float) op->i3 - 1.0F;

                  /* mix flag is set, so average the prior target
                     coordinates with these coordinates */

                  vt2 = op->vv2;
                  t_i = 0;      /* original target vertex index */
                  for(a = 0; a < I->NAtom; a++) {
                    s = I->AtomInfo[a].selEntry;
                    if(SelectorIsMember(G, s, sele)) {
                      a1 = I->CSet[b]->atmToIdx(a);
                      if(a1 >= 0) {

                        match_flag = false;
                        while(t_i < op->nvv2) {
                          if(op->i1VLA[t_i] == a) {     /* same atom? */
                            match_flag = true;
                            break;
                          }
                          if(op->i1VLA[t_i] < a) {      /* catch up? */
                            t_i++;
                            vt2 += 3;
                          } else
                            break;
                        }
                        if(match_flag) {
                          vv2 = I->CSet[b]->coordPtr(a1);
                          *(vt2) = ((premult * (*vt2)) + *(vv2++)) / divisor;
                          *(vt2 + 1) = ((premult * (*(vt2 + 1))) + *(vv2++)) / divisor;
                          *(vt2 + 2) = ((premult * (*(vt2 + 2))) + *(vv2++)) / divisor;
                        }
                      }
                    }
                  }
                }
              }
            } else {
              PRINTFB(G, FB_Executive, FB_Warnings)
                "Executive-Warning: No matches found for state %d.\n", b + 1 ENDFB(G);
            }
          }
          VLACheck(op->f1VLA, float, b);
          op->f1VLA[b] = rms;
        }
        VLASize(op->f1VLA, float, I->NCSet);    /* NOTE this action is object-specific! */
      }
      FreeP(vt);
      break;
    case OMOP_OnOff:
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          hit_flag = true;
          break;
        }
      }
      break;
    case OMOP_SaveUndo:        /* save undo */
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          hit_flag = true;
          break;
        }
      }
      break;
    case OMOP_IdentifyObjects: /* identify atoms */
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          VLACheck(op->i1VLA, int, op->i1);
          op->i1VLA[op->i1] = I->AtomInfo[a].id;
          if (op->obj1VLA != nullptr) {
            VLACheck(op->obj1VLA, ObjectMolecule*, op->i1);
            op->obj1VLA[op->i1] = I;
          }
          op->i1++;
        }
      }
      break;
    case OMOP_Index:           /* identify atoms */
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          VLACheck(op->i1VLA, int, op->i1);
          op->i1VLA[op->i1] = a;        /* NOTE: need to incr by 1 before python */
          VLACheck(op->obj1VLA, ObjectMolecule *, op->i1);
          op->obj1VLA[op->i1] = I;
          op->i1++;
        }
      }
      break;
    case OMOP_GetObjects:      /* identify atoms */
      for(a = 0; a < I->NAtom; a++) {
        s = I->AtomInfo[a].selEntry;
        if(SelectorIsMember(G, s, sele)) {
          VLACheck(op->obj1VLA, ObjectMolecule *, op->i1);
          op->obj1VLA[op->i1] = I;
          op->i1++;
          break;
        }
      }
      break;
    case OMOP_CountAtoms:      /* count atoms in object, in selection */
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele))
          op->i1++;
        ai++;
      }
      break;
    case OMOP_PhiPsi:
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          VLACheck(op->i1VLA, int, op->i1);
          op->i1VLA[op->i1] = a;
          VLACheck(op->obj1VLA, ObjectMolecule *, op->i1);
          op->obj1VLA[op->i1] = I;
          VLACheck(op->f1VLA, float, op->i1);
          VLACheck(op->f2VLA, float, op->i1);
          if(ObjectMoleculeGetPhiPsi
             (I, a, op->f1VLA + op->i1, op->f2VLA + op->i1, op->i2))
            op->i1++;
        }
        ai++;
      }
      break;
    case OMOP_Cartoon:         /* adjust cartoon type */
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          if (ai->cartoon!=op->i1){
            op->i3++;
          }
          ai->cartoon = op->i1;
          op->i2++;
        }
        ai++;
      }
      break;
    case OMOP_Protect:         /* protect atoms from movement */
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          ai->protekted = op->i1;
          op->i2++;
        }
        ai++;
      }
      break;
    case OMOP_Mask:            /* protect atoms from selection */
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          ai->masked = op->i1;
          op->i2++;
        }
        ai++;
      }
      break;
    case OMOP_SetB:            /* set B-value */
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          ai->b = op->f1;
          op->i2++;
        }
        ai++;
      }
      break;
    case OMOP_Remove:          /* flag atoms for deletion */
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
	ai->deleteFlag = false;
	s = ai->selEntry;
	if(SelectorIsMember(G, s, sele)) {
	  ai->deleteFlag = true;
	  op->i1++;
	}
	ai++;
      }
      break;
    case OMOP_GetChains:
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          // pointer hack
          ((std::set<lexidx_t> *) (void*) op->ii1)->insert(ai->chain);
          op->i1++;
        }
        ai++;
      }
      break;

    case OMOP_Spectrum:
      ai = I->AtomInfo.data();
      ai0 = NULL;
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          skip_flag = false;
          if(op->i4 && ai0)     /* byres and we've done a residue */
            if(AtomInfoSameResidue(G, ai, ai0))
              skip_flag = true;
          if(!skip_flag) {
            c = (int) (0.49999 + op->i1 * (op->ff1[op->i3] - op->f1) / op->f2);
            if(c < 0)
              c = 0;
            if(c >= op->i1)
              c = op->i1 - 1;
            ai->color = op->ii1[c];

            /*               printf("%8.3 %8.3\n",ai->partial_charge, */
            if(op->i4) {        /* byres */
              offset = -1;
              while((a + offset) >= 0) {
                ai0 = I->AtomInfo + a + offset;
                if(AtomInfoSameResidue(G, ai, ai0)) {
                  ai0->color = op->ii1[c];
                  hit_flag = true;
                } else
                  break;
                offset--;
              }
              offset = 1;
              while((a + offset) < I->NAtom) {
                ai0 = I->AtomInfo + a + offset;
                if(AtomInfoSameResidue(G, ai, ai0)) {
                  ai0->color = op->ii1[c];
                  hit_flag = true;
                } else
                  break;
                offset++;
              }
            }
            ai0 = ai;
          }
          op->i3++;

        }
        ai++;
      }
      break;

    case OMOP_SingleStateVertices:     /* same as OMOP_VERT for a single state */
      ai = I->AtomInfo.data();
      if(op->cs1 < I->NCSet) {
        if(I->CSet[op->cs1]) {
          b = op->cs1;
          for(a = 0; a < I->NAtom; a++) {
            s = ai->selEntry;
            if(SelectorIsMember(G, s, sele)) {
              op->i1++;
              a1 = I->CSet[b]->atmToIdx(a);
              if(a1 >= 0) {
                VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                vv2 = I->CSet[b]->coordPtr(a1);
                vv1 = op->vv1 + (op->nvv1 * 3);
                *(vv1++) = *(vv2++);
                *(vv1++) = *(vv2++);
                *(vv1++) = *(vv2++);
                op->nvv1++;
              }
            }
            ai++;
          }
        }
      }
      break;
    case OMOP_CSetIdxGetAndFlag:
      ai = I->AtomInfo.data();
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          for(b = op->cs1; b <= op->cs2; b++) {
            offset = b - op->cs1;
            if(b < I->NCSet) {
              if(I->CSet[b]) {
                a1 = I->CSet[b]->atmToIdx(a);
                if(a1 >= 0) {
                  op->ii1[op->i1 * offset + op->i2] = 1;        /* presence flag */
                  vv1 = op->vv1 + 3 * (op->i1 * offset + op->i2);       /* atom-based offset */
                  vv2 = I->CSet[b]->coordPtr(a1);
                  *(vv1++) = *(vv2++);
                  *(vv1++) = *(vv2++);
                  *(vv1++) = *(vv2++);
                  op->nvv1++;
                }
              }
            }
          }
          op->i2++;             /* atom index field for atoms within selection... */
        }
        ai++;
      }
      break;
    case OMOP_CSetIdxSetFlagged:
      ai = I->AtomInfo.data();
      hit_flag = false;
      for(a = 0; a < I->NAtom; a++) {
        s = ai->selEntry;
        if(SelectorIsMember(G, s, sele)) {
          for(b = op->cs1; b <= op->cs2; b++) {
            offset = b - op->cs1;
            if(b < I->NCSet) {
              if(I->CSet[b]) {
                a1 = I->CSet[b]->atmToIdx(a);
                if(a1 >= 0) {
                  if(op->ii1[op->i1 * offset + op->i2]) {       /* copy flag */
                    vv1 = op->vv1 + 3 * (op->i1 * offset + op->i2);     /* atom-based offset */
                    vv2 = I->CSet[b]->coordPtr(a1);
                    *(vv2++) = *(vv1++);
                    *(vv2++) = *(vv1++);
                    *(vv2++) = *(vv1++);
                    op->nvv1++;
                    hit_flag = true;
                  }
                }
              }
            }
          }
          op->i2++;             /* atom index field for atoms within selection... */
        }
        ai++;
      }
      break;
    case OMOP_SUMC:            /* performance optimized to speed center & zoom actions */
      {
				/* given a selection, sum up all the coordinates (for centering) */
        float *op_v1 = op->v1;
        int op_i1 = op->i1;
        int op_i2 = op->i2;
        int obj_TTTFlag = I->TTTFlag;
        int i_NCSet = I->NCSet;
        int i_NAtom = I->NAtom;
        int i_DiscreteFlag = I->DiscreteFlag;
        CoordSet* const* i_CSet = I->CSet.data();
        if(op_i2) {
          use_matrices =
            SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);
          if(use_matrices<0) use_matrices = 0;
        }
        ai = I->AtomInfo.data();
        for(a = 0; a < i_NAtom; a++) {
          s = ai->selEntry;
					/* for each atom, if this current atom is in the selection */
          if(SelectorIsMember(G, s, sele)) {
						/* loop over all atoms; regardless of state */
            for(b = 0; b < i_NCSet; b++) {
              if(i_DiscreteFlag) {
                if((cs = I->DiscreteCSet[a]))
                  a1 = I->DiscreteAtmToIdx[a];
              } else {
                if((cs = i_CSet[b]))
                  a1 = cs->AtmToIdx[a];
              }
							/* if valid coordinate set and atom info for this atom */
              if(cs && (a1 >= 0)) {
                coord = cs->coordPtr(a1);
                if(op_i2) {     /* do we want transformed coordinates? */
                  if(use_matrices) {
                    if(!cs->Matrix.empty()) {      /* state transformation */
                      transform44d3f(cs->Matrix.data(), coord, v1);
                      coord = v1;
                    }
                  }
                  if(obj_TTTFlag) {
                    transformTTT44f3f(I->TTT, coord, v1);
                    coord = v1;
                  }
                }
								/* op_v1 += coord */
                add3f(op_v1, coord, op_v1);
								/* count += 1 */
                op_i1++;
              }
              if(i_DiscreteFlag)
                break;
            }
          }
					/* next atom */
          ai++;
        }
        op->i1 = op_i1;
      }
      break;
    case OMOP_MNMX:            /* performance optimized to speed center & zoom actions */
      {
        float *op_v1 = op->v1;
        float *op_v2 = op->v2;
        int op_i1 = op->i1;
        int op_i2 = op->i2;
        int obj_TTTFlag = I->TTTFlag;
        int i_NCSet = I->NCSet;
        int i_NAtom = I->NAtom;
        int i_DiscreteFlag = I->DiscreteFlag;
        CoordSet* const* i_CSet = I->CSet.data();
        if(op_i2) {
          use_matrices =
            SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);
          if(use_matrices<0) use_matrices = 0;
        }
        ai = I->AtomInfo.data();
        for(a = 0; a < i_NAtom; a++) {
          s = ai->selEntry;
          if(SelectorIsMember(G, s, sele)) {
            for(b = 0; b < i_NCSet; b++) {
              if(i_DiscreteFlag) {
                if((cs = I->DiscreteCSet[a]))
                  a1 = I->DiscreteAtmToIdx[a];
              } else {
                if((cs = i_CSet[b]))
                  a1 = cs->AtmToIdx[a];
              }
              if(cs && (a1 >= 0)) {
                coord = cs->coordPtr(a1);
                if(op_i2) {     /* do we want transformed coordinates? */
                  if(use_matrices) {
                    if(!cs->Matrix.empty()) {      /* state transformation */
                      transform44d3f(cs->Matrix.data(), coord, v1);
                      coord = v1;
                    }
                  }
                  if(obj_TTTFlag) {
                    transformTTT44f3f(I->TTT, coord, v1);
                    coord = v1;
                  }
                }
                if(op_i1) {
                  if(op_v1[0] > coord[0])
                    op_v1[0] = coord[0];
                  if(op_v1[1] > coord[1])
                    op_v1[1] = coord[1];
                  if(op_v1[2] > coord[2])
                    op_v1[2] = coord[2];
                  if(op_v2[0] < coord[0])
                    op_v2[0] = coord[0];
                  if(op_v2[1] < coord[1])
                    op_v2[1] = coord[1];
                  if(op_v2[2] < coord[2])
                    op_v2[2] = coord[2];
                } else {
                  op_v1[0] = coord[0];
                  op_v1[1] = coord[1];
                  op_v1[2] = coord[2];
                  op_v2[0] = coord[0];
                  op_v2[1] = coord[1];
                  op_v2[2] = coord[2];
                }
                op_i1++;
              }
              if(i_DiscreteFlag)
                break;
            }
          }
          ai++;
        }
        op->i1 = op_i1;
      }
      break;
    default:
      {
        int inv_flag;
#ifdef _PYMOL_IP_EXTRAS
	int use_stereo = 0, use_text_type = 0;
#endif

        switch (op->code) {
        case OMOP_INVA:
          /* set up an important optimization... */
          for(b = 0; b < I->NCSet; b++) {
            cs = I->CSet[b];
            if(cs)
              cs->objMolOpInvalidated = false;
          }
          break;
        }
        use_matrices = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);
        if(use_matrices<0) use_matrices = 0;
        ai = I->AtomInfo.data();

#ifdef _PYMOL_IP_EXTRAS
        // use stereo or text_type ?
        // only do this for "label2" command (better logic in WrapperObjectSubScript)
	if (op->code == OMOP_LABL && op->i2 == cExecutiveLabelEvalAlt) {
	  use_stereo = PLabelExprUsesVariable(G, op->s1, "stereo");
	  use_text_type = PLabelExprUsesVariable(G, op->s1, "text_type");
	}

#ifdef NO_MMLIBS
        if (use_stereo) {
          PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
            " NO_MMLIBS-Warning: stereochemistry not supported in this PyMOL build.\n" ENDFB(G);
        }
        if (use_text_type) {
          PRINTFB(G, FB_ObjectMolecule, FB_Warnings)
            " NO_MMLIBS-Warning: automatic 'text_type' assignment not supported in this PyMOL build.\n" ENDFB(G);
        }
#endif
#endif

        for(a = 0; a < I->NAtom; a++) {
          switch (op->code) {
          case OMOP_Flag:
            ai->flags &= op->i2;        /* clear flag using mask */
            op->i4++;
            /* no break here is intentional!  */
          case OMOP_FlagSet:
          case OMOP_FlagClear:
          case OMOP_COLR:      /* normal atom based loops */
          case OMOP_VISI:
          case OMOP_CheckVis:
          case OMOP_TTTF:
          case OMOP_ALTR:
          case OMOP_LABL:
          case OMOP_AlterState:
            s = ai->selEntry;
            if(SelectorIsMember(G, s, sele)) {
              switch (op->code) {
              case OMOP_Flag:
              case OMOP_FlagSet:
                ai->flags |= op->i1;    /* set flag */
                op->i3++;
                break;
              case OMOP_FlagClear:
                ai->flags &= op->i2;    /* clear flag */
                op->i3++;
                break;
              case OMOP_VISI:
                switch (op->i2) {
                case cVis_HIDE:
                  ai->visRep &= ~(op->i1);
                  I->visRep &= ~(op->i1); // cell
                  break;
                case cVis_SHOW:
                  ai->visRep |= op->i1 & cRepsAtomMask;
                  I->visRep |= op->i1 & cRepsObjectMask; // cell
                  break;
                case cVis_AS:
                  ai->visRep = op->i1 & cRepsAtomMask;
                  I->visRep = op->i1 & cRepsObjectMask; // cell
                  break;
                }
                break;
              case OMOP_CheckVis:
                if((ai->visRep & op->i1)) {
                  op->i2 = true;
                }
                break;
              case OMOP_COLR:
                if(op->i1 == cColorAtomic)
                  ai->color = AtomInfoGetColor(G, ai);
                else
                  ai->color = op->i1;
                hit_flag = true;
                op->i2++;
                break;
              case OMOP_TTTF:
                hit_flag = true;
                break;
              case OMOP_LABL:
                if(ok) {
                  if(!op->s1[0]) {
                    if(ai->label) {
			op->i1--; /* negative if unlabelling */
                      LexAssign(G, ai->label, 0);
                    }
                    ai->visRep &= ~cRepLabelBit;
                    hit_flag = true;
                  } else {
                    switch (op->i2) {
                    case cExecutiveLabelEvalOn:
		      {
			/* python label expression evaluation */
			CoordSet *cs = NULL;
			if(I->DiscreteFlag && I->DiscreteCSet) {
			  cs = I->DiscreteCSet[a];
			} else if (I->NCSet == 1){
			  cs = I->CSet[0];
			}
#ifndef _PYMOL_NOPY
			if(PLabelAtom(I->G, I, cs, expr_co, a)) {
			  if (ai->label){
			    op->i1++; /* only if the string has been set, report labelled */
			  }
			  ai->visRep |= cRepLabelBit;
			  hit_flag = true;
			} else {
			  ok = false;
			}
#else
			ok = false;
#endif
		      }
                      break;
                    case cExecutiveLabelEvalAlt:
		      {
			if(PLabelAtomAlt(I->G, &I->AtomInfo[a], I->Name, op->s1, a)) {
			  if (ai->label){
			    op->i1++; /* only if the string has been set, report labelled */
			  }
			  ai->visRep |= cRepLabelBit;
			  hit_flag = true;
			} else {
			  ok = false;
			}
		      }
                      break;
                    case cExecutiveLabelEvalOff:
                      {
                        /* simple string label text */
                        AtomInfoType *ai = I->AtomInfo + a;
                        LexDec(G, ai->label);
                        ai->label = LexIdx(G, op->s1);
                      }
                      break;
                    }
                  }
                }
                break;
              case OMOP_ALTR:
                if(ok) {
		  CoordSet *cs = NULL;
		  if(I->DiscreteFlag && I->DiscreteCSet) {
		    cs = I->DiscreteCSet[a];
		  } else if (I->NCSet == 1){
		    cs = I->CSet[0];
		  }
#ifndef _PYMOL_NOPY
                  if(PAlterAtom
                     (I->G, I, cs, expr_co, op->i2, a,
                      op->py_ob1))
                    op->i1++;
                  else
#endif
#ifdef _WEBGL
#endif
                    ok = false;
                }
                break;
              case OMOP_AlterState:
                if(ok) {
                  if(op->i2 < I->NCSet) {
                    cs = I->CSet[op->i2];
                    if(cs) {
                      a1 = cs->atmToIdx(a);
                      if(a1 >= 0) {
#ifndef _PYMOL_NOPY
                        if(PAlterAtomState(I->G, expr_co, op->i3,
                                           I, cs, a, a1, op->i2, op->py_ob1)) {
                          op->i1++;
                          hit_flag = true;
                        } else
#endif
#ifdef _WEBGL
#endif
                          ok = false;
                      }
                    }
                  }
                }
                break;
              }
              break;
            }
            break;

            /* coord-set based properties, iterating only a single coordinate set */
          case OMOP_CSetMinMax:
          case OMOP_CSetCameraMinMax:
          case OMOP_CSetMaxDistToPt:
          case OMOP_CSetSumSqDistToPt:
          case OMOP_CSetSumVertices:
          case OMOP_CSetMoment:
            cs = NULL;
            if((op->cs1 >= 0) && (op->cs1 < I->NCSet)) {
              cs = I->CSet[op->cs1];
            } else if(op->include_static_singletons) {
              if((I->NCSet == 1)
                 && (SettingGet_b(G, NULL, I->Setting.get(), cSetting_static_singletons))) {
                cs = I->CSet[0];        /*treat static singletons as present in each state */
              }
            }

            if(cs) {
              s = ai->selEntry;
              if(SelectorIsMember(G, s, sele)) {
                switch (op->code) {
                case OMOP_CSetSumVertices:
                  a1 = cs->atmToIdx(a);
                  if(a1 >= 0) {
                    coord = cs->coordPtr(a1);
                    if(op->i2) {        /* do we want transformed coordinates? */
                      if(use_matrices) {
                        if(!cs->Matrix.empty()) {  /* state transformation */
                          transform44d3f(cs->Matrix.data(), coord, v1);
                          coord = v1;
                        }
                      }
                      if(I->TTTFlag) {
                        transformTTT44f3f(I->TTT, coord, v1);
                        coord = v1;
                      }
                    }
                    add3f(op->v1, coord, op->v1);
                    op->i1++;
                  }
                  break;
                case OMOP_CSetMinMax:
                  a1 = cs->atmToIdx(a);
                  if(a1 >= 0) {
                    coord = cs->coordPtr(a1);
                    if(op->i2) {        /* do we want transformed coordinates? */
                      if(use_matrices) {
                        if(!cs->Matrix.empty()) {  /* state transformation */
                          transform44d3f(cs->Matrix.data(), coord, v1);
                          coord = v1;
                        }
                      }
                      if(I->TTTFlag) {
                        transformTTT44f3f(I->TTT, coord, v1);
                        coord = v1;
                      }
                    }
                    if(op->i1) {
                      for(c = 0; c < 3; c++) {
                        if(*(op->v1 + c) > *(coord + c))
                          *(op->v1 + c) = *(coord + c);
                        if(*(op->v2 + c) < *(coord + c))
                          *(op->v2 + c) = *(coord + c);
                      }
                    } else {
                      for(c = 0; c < 3; c++) {
                        *(op->v1 + c) = *(coord + c);
                        *(op->v2 + c) = *(coord + c);
                      }
                    }
                    op->i1++;
                  }
                  break;
                case OMOP_CSetCameraMinMax:
                  a1 = cs->atmToIdx(a);
                  if(a1 >= 0) {
                    coord = cs->coordPtr(a1);
                    if(op->i2) {        /* do we want transformed coordinates? */
                      if(use_matrices) {
                        if(!cs->Matrix.empty()) {  /* state transformation */
                          transform44d3f(cs->Matrix.data(), coord, v1);
                          coord = v1;
                        }
                      }
                      if(I->TTTFlag) {
                        transformTTT44f3f(I->TTT, coord, v1);
                        coord = v1;
                      }
                    }
                    MatrixTransformC44fAs33f3f(op->mat1, coord, v1);
                    /* convert to view-space */
                    coord = v1;
                    if(op->i1) {
                      for(c = 0; c < 3; c++) {
                        if(*(op->v1 + c) > *(coord + c))
                          *(op->v1 + c) = *(coord + c);
                        if(*(op->v2 + c) < *(coord + c))
                          *(op->v2 + c) = *(coord + c);
                      }
                    } else {
                      for(c = 0; c < 3; c++) {
                        *(op->v1 + c) = *(coord + c);
                        *(op->v2 + c) = *(coord + c);
                      }
                    }
                    op->i1++;
                  }
                  break;
                case OMOP_CSetSumSqDistToPt:
                  a1 = cs->atmToIdx(a);
                  if(a1 >= 0) {
                    float dist;
                    coord = cs->coordPtr(a1);
                    if(op->i2) {        /* do we want transformed coordinates? */
                      if(use_matrices) {
                        if(!cs->Matrix.empty()) {  /* state transformation */
                          transform44d3f(cs->Matrix.data(), coord, v1);
                          coord = v1;
                        }
                      }
                      if(I->TTTFlag) {
                        transformTTT44f3f(I->TTT, coord, v1);
                        coord = v1;
                      }
                    }
                    dist = (float) diff3f(op->v1, coord);
                    op->d1 += dist * dist;
                    op->i1++;
                  }
                  break;
                case OMOP_CSetMaxDistToPt:
                  a1 = cs->atmToIdx(a);
                  if(a1 >= 0) {
                    float dist;
                    coord = cs->coordPtr(a1);
                    if(op->i2) {        /* do we want transformed coordinates? */
                      if(use_matrices) {
                        if(!cs->Matrix.empty()) {  /* state transformation */
                          transform44d3f(cs->Matrix.data(), coord, v1);
                          coord = v1;
                        }
                      }
                      if(I->TTTFlag) {
                        transformTTT44f3f(I->TTT, coord, v1);
                        coord = v1;
                      }
                    }
                    dist = (float) diff3f(op->v1, coord);
                    if(dist > op->f1)
                      op->f1 = dist;
                    op->i1++;
                  }
                  break;
                case OMOP_CSetMoment:
                  a1 = cs->atmToIdx(a);
                  if(a1 >= 0) {
                    subtract3f(cs->coordPtr(a1), op->v1, v1);
                    v2 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
                    op->d[0][0] += v2 - v1[0] * v1[0];
                    op->d[0][1] += -v1[0] * v1[1];
                    op->d[0][2] += -v1[0] * v1[2];
                    op->d[1][0] += -v1[1] * v1[0];
                    op->d[1][1] += v2 - v1[1] * v1[1];
                    op->d[1][2] += -v1[1] * v1[2];
                    op->d[2][0] += -v1[2] * v1[0];
                    op->d[2][1] += -v1[2] * v1[1];
                    op->d[2][2] += v2 - v1[2] * v1[2];
                  }
                  break;

                }
              }
            }
            break;
          default:
            /* coord-set based properties, iterating as all coordsets within atoms */
            for(b = 0; b < I->NCSet; b++) {
              if(I->DiscreteFlag) {
                cs = I->DiscreteCSet[a];
              } else {
                cs = I->CSet[b];
              }
              if(cs) {
                s = ai->selEntry;
                inv_flag = false;
                if(SelectorIsMember(G, s, sele)) {
                  switch (op->code) {
                  case OMOP_CameraMinMax:
                    a1 = cs->atmToIdx(a);
                    if(a1 >= 0) {
                      coord = cs->coordPtr(a1);
                      if(op->i2) {      /* do we want transformed coordinates? */
                        if(use_matrices) {
                          if(!cs->Matrix.empty()) {        /* state transformation */
                            transform44d3f(cs->Matrix.data(), coord, v1);
                            coord = v1;
                          }
                        }
                        if(I->TTTFlag) {
                          transformTTT44f3f(I->TTT, coord, v1);
                          coord = v1;
                        }
                      }
                      MatrixTransformC44fAs33f3f(op->mat1, coord, v1);
                      /* convert to view-space */
                      coord = v1;
                      if(op->i1) {
                        for(c = 0; c < 3; c++) {
                          if(*(op->v1 + c) > *(coord + c))
                            *(op->v1 + c) = *(coord + c);
                          if(*(op->v2 + c) < *(coord + c))
                            *(op->v2 + c) = *(coord + c);
                        }
                      } else {
                        for(c = 0; c < 3; c++) {
                          *(op->v1 + c) = *(coord + c);
                          *(op->v2 + c) = *(coord + c);
                        }
                      }
                      op->i1++;
                    }
                    break;
                  case OMOP_MaxDistToPt:
                    a1 = cs->atmToIdx(a);
                    if(a1 >= 0) {
                      float dist;
                      coord = cs->coordPtr(a1);
                      if(op->i2) {      /* do we want transformed coordinates? */
                        if(use_matrices) {
                          if(!cs->Matrix.empty()) {        /* state transformation */
                            transform44d3f(cs->Matrix.data(), coord, v1);
                            coord = v1;
                          }
                        }
                        if(I->TTTFlag) {
                          transformTTT44f3f(I->TTT, coord, v1);
                          coord = v1;
                        }
                      }
                      dist = (float) diff3f(coord, op->v1);
                      if(dist > op->f1)
                        op->f1 = dist;
                      op->i1++;
                    }
                    break;
                  case OMOP_INVA:
                    if(!cs->objMolOpInvalidated) {
                      /* performance optimization: avoid repeatedly
                         calling invalidate on the same coord sets */
                      a1 = cs->atmToIdx(a);
                      if(a1 >= 0)       /* selection touches this coordinate set */
                        inv_flag = true;        /* so set the invalidation flag */
                    }
                    break;
                  case OMOP_VERT:
										/* get the atom index whether it's discrete or not */
                    a1 = cs->atmToIdx(a);
                    if(a1 >= 0) {
											/* if a1 is a valid atom index, then copy it's xyz coordinates
											 * into vv1; increment the counter, nvv1 */
                      VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                      vv2 = cs->coordPtr(a1);
                      vv1 = op->vv1 + (op->nvv1 * 3);
                      *(vv1++) = *(vv2++);
                      *(vv1++) = *(vv2++);
                      *(vv1++) = *(vv2++);
                      op->nvv1++;
                    }
                    break;
                  case OMOP_SVRT:
                    /* gives us only vertices for a specific coordinate set */
                    if(b == op->i1) {
                      a1 = cs->atmToIdx(a);
                      if(a1 >= 0) {
                        VLACheck(op->vv1, float, (op->nvv1 * 3) + 2);
                        VLACheck(op->i1VLA, int, op->nvv1);
                        op->i1VLA[op->nvv1] = a;        /* save atom index for later comparisons */
                        vv2 = cs->coordPtr(a1);
                        vv1 = op->vv1 + (op->nvv1 * 3);
                        *(vv1++) = *(vv2++);
                        *(vv1++) = *(vv2++);
                        *(vv1++) = *(vv2++);
                        op->nvv1++;
                      }
                    }
                    break;
                  case OMOP_MOME:
                    /* Moment of inertia tensor - unweighted - assumes v1 is center of molecule */
                    a1 = cs->atmToIdx(a);
                    if(a1 >= 0) {
                      subtract3f(cs->coordPtr(a1), op->v1, v1);
                      v2 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
                      op->d[0][0] += v2 - v1[0] * v1[0];
                      op->d[0][1] += -v1[0] * v1[1];
                      op->d[0][2] += -v1[0] * v1[2];
                      op->d[1][0] += -v1[1] * v1[0];
                      op->d[1][1] += v2 - v1[1] * v1[1];
                      op->d[1][2] += -v1[1] * v1[2];
                      op->d[2][0] += -v1[2] * v1[0];
                      op->d[2][1] += -v1[2] * v1[1];
                      op->d[2][2] += v2 - v1[2] * v1[2];
                    }
                    break;
                  }
                }
                switch (op->code) {
                  /* full coord-set based */
                case OMOP_INVA:
                  /* shouldn't this be calling the object invalidation routine instead? */
                  if(inv_flag) {
                    cs->objMolOpInvalidated = true;
                    for (auto d = cRep_t(0); d < cRepCnt; ++d) {
                      if ((1 << d) & op->i1) {
                        cs->invalidateRep(d, cRepInv_t(op->i2));
                      }
                    }
                  }
                  break;
                }
                if(I->DiscreteFlag) {
                  /* don't iterate every coordinate set for discrete objects! */
                  break;
                }
              }
            }                   /* end coordset section */
            break;
          }
          ai++;
        }
      }                         /* case: default */
      break;
    }
    if(hit_flag) {
      switch (op->code) {
      case OMOP_COLR:
        ExecutiveUpdateColorDepends(I->G, I);
        break;
      case OMOP_TTTF:
        ObjectMoleculeTransformTTTf(I, op->ttt, -1);
        break;
      case OMOP_LABL:
        I->invalidate(cRepLabel, cRepInvText, -1);
        break;
      case OMOP_AlterState:    /* overly coarse - doing all states, could do just 1 */
        if(!op->i3) {           /* not read_only? */
          I->invalidate(cRepAll, cRepInvRep, -1);
          SceneChanged(G);
        }
        break;
      case OMOP_CSetIdxSetFlagged:
        I->invalidate(cRepAll, cRepInvRep, -1);
        SceneChanged(G);
        break;
      case OMOP_SaveUndo:
        op->i2 = true;
        ObjectMoleculeSaveUndo(I, op->i1, false);
        break;
      case OMOP_OnOff:
        ExecutiveSetObjVisib(G, I->Name, op->i1, false);
        break;
      case OMOP_RevalenceFromSource:
        if(ObjectMoleculeXferValences(I, op->i1, op->i2,
                                      op->i3, op->obj3, op->i4, op->i5, op->i6)) {
          ObjectMoleculeVerifyChemistry(I, op->i3);
          I->invalidate(cRepAll, cRepInvBonds, op->i3);
        }
        break;
      case OMOP_RevalenceByGuessing:
        {
          int *flag1 = pymol::calloc<int>(I->NAtom);
          int *flag2 = pymol::calloc<int>(I->NAtom);
          if(flag1 && flag2) {
            int a;
            int *f1 = flag1;
            int *f2 = flag2;
            const AtomInfoType *ai = I->AtomInfo.data();
            for(a = 0; a < I->NAtom; a++) {
              *(f1++) = SelectorIsMember(G, ai->selEntry, op->i1);
              *(f2++) = SelectorIsMember(G, ai->selEntry, op->i2);
              ai++;
            }
            {
              int target_state = op->i3;
              if(target_state < 0)
                target_state = 0;       /* TO DO */
              ObjectMoleculeGuessValences(I, target_state, flag1, flag2, op->i4);
              ObjectMoleculeVerifyChemistry(I, target_state);
              I->invalidate(cRepAll, cRepInvBonds, target_state);
            }
            FreeP(flag1);
            FreeP(flag2);
          }
        }
        break;
      }
    }

    /* always run on exit... */
    switch (op->code) {
    case OMOP_LABL:
      if (op->i2 != cExecutiveLabelEvalOn){
	break;
      }
    case OMOP_ALTR:
    case OMOP_AlterState:
#ifndef _PYMOL_NOPY
      Py_XDECREF(expr_co);
#endif
      break;
    }
    /* */
  }

  return ok;
}


/*========================================================================*/
void ObjectMoleculeGetAtomSele(const ObjectMolecule * I, int index, char *buffer)
{
  PyMOLGlobals * G = I->G;
  assert(index < I->NAtom);
  auto* ai = I->AtomInfo + index;
  char inscode_str[2] = { ai->inscode, '\0' };

  snprintf(buffer, OrthoLineLength, "/%s/%s/%s/%s`%d%s/%s`%s", I->Name,
      LexStr(G, ai->segi),
      LexStr(G, ai->chain),
      LexStr(G, ai->resn), ai->resv, inscode_str,
      LexStr(G, ai->name), ai->alt);
}

/**
 * Get Atom selection string
 * @param I target Object Molecule
 * @param index atom index of I
 */

std::string ObjectMoleculeGetAtomSele(const ObjectMolecule * I, int index)
{
  PyMOLGlobals * G = I->G;
  assert(index < I->NAtom);
  auto* ai = I->AtomInfo + index;
  char inscode_str[2] = { ai->inscode, '\0' };
  std::string buffer = pymol::string_format("/%s/%s/%s/%s`%d%s/%s`%s", I->Name,
      LexStr(G, ai->segi),
      LexStr(G, ai->chain),
      LexStr(G, ai->resn), ai->resv, inscode_str,
      LexStr(G, ai->name), ai->alt);
  return buffer;
}

/*========================================================================*/

/**
 * Get Atom selection string
 * @param I target Object Molecule
 * @param index atom index of I
 */

std::string ObjectMolecule::describeElement(int index) const
{
  auto I = this;
  auto buffer = ObjectMoleculeGetAtomSele(I, index);
  if(!I->AtomInfo[index].alt[0]) {
    // don't include the trailing backtick
    buffer.pop_back();
  }
  return buffer;
}


/*========================================================================*/

std::string ObjectMoleculeGetAtomSeleLog(const ObjectMolecule* I, int index, int quote)
{
  OrthoLineType buffer;
  ObjectMoleculeGetAtomSeleLog(I, index, buffer, quote);
  return buffer;
}

void ObjectMoleculeGetAtomSeleLog(const ObjectMolecule * I, int index, char *buffer, int quote)
{
  char *p = quote ? buffer + 1 : buffer;

  if(SettingGetGlobal_b(I->G, cSetting_robust_logs)) {
    ObjectMoleculeGetAtomSele(I, index, p);
  } else {
    sprintf(p, "(%s`%d)", I->Name, index + 1);
  }

  if (quote) {
    int len = strlen(p);
    buffer[0] = buffer[len + 1] = '"';
    buffer[len + 2] = 0;
  }
}

std::string ObjectMoleculeGetAtomSeleFast(const ObjectMolecule* I, int index)
{
  auto G = I->G;
  auto* ai = I->AtomInfo + index;
  char inscodestr[2] = {ai->inscode, 0};
  return pymol::string_format("(/'%s'/'%s'/'%s'/'%s'`%d%s/'%s'`'%s')", I->Name,
      LexStr(G, ai->segi), LexStr(G, ai->chain), LexStr(G, ai->resn), ai->resv,
      inscodestr, ai->name, ai->alt);
}

/*========================================================================*/
int ObjectMolecule::getNFrame() const
{
  return NCSet;
}

struct _CCoordSetUpdateThreadInfo {
  CoordSet *cs;
  int a;
};

void CoordSetUpdateThread(CCoordSetUpdateThreadInfo * T)
{
  if(T->cs) {
    T->cs->update(T->a);
  }
}

#ifndef _PYMOL_NOPY
static void ObjMolCoordSetUpdateSpawn(PyMOLGlobals * G,
                                      CCoordSetUpdateThreadInfo * Thread, int n_thread,
                                      int n_total)
{
  if(n_total == 1) {
    CoordSetUpdateThread(Thread);
  } else if(n_total) {
    int blocked;
    PyObject *info_list;
    int a, n = 0;
    blocked = PAutoBlock(G);

    PRINTFB(G, FB_Scene, FB_Blather)
      " Scene: updating coordinate sets with %d threads...\n", n_thread ENDFB(G);
    info_list = PyList_New(n_total);
    for(a = 0; a < n_total; a++) {
      PyList_SetItem(info_list, a, PyCapsule_New(Thread + a, nullptr, nullptr));
      n++;
    }
    PXDecRef(PYOBJECT_CALLMETHOD
             (G->P_inst->cmd, "_coordset_update_spawn", "Oi", info_list, n_thread));
    Py_DECREF(info_list);
    PAutoUnblock(G, blocked);
  }
}
#endif


/*========================================================================*/
void ObjectMolecule::update()
{
  auto I = this;
  int a; /*, ok; */

  OrthoBusyPrime(G);
  /* if the cached representation is invalid, reset state */
  if(!I->RepVisCacheValid) {
    /* note which representations are active */
    /* for each atom in each coordset, blank out the representation cache */
    if(I->NCSet > 1) {
      const AtomInfoType *ai = I->AtomInfo.data();
      I->RepVisCache = 0;
      for(a = 0; a < I->NAtom; a++) {
        I->RepVisCache |= ai->visRep;
        ai++;
      }
    } else {
      I->RepVisCache = cRepBitmask;     /* if only one coordinate set, then
                                         * there's no benefit to pre-filtering
                                         * the representations... */
    }
    I->RepVisCacheValid = true;
  }
  {
    /* determine the start/stop states */
    int start = 0;
    int stop = I->NCSet;
    /* set start and stop given an object */
    ObjectAdjustStateRebuildRange(I, &start, &stop);
    if((I->NCSet == 1)
       && (SettingGet_b(G, I->Setting.get(), NULL, cSetting_static_singletons))) {
      start = 0;
      stop = 1;
    }
    if(stop > I->NCSet)
      stop = I->NCSet;

    /* single and multithreaded coord set updates */
    {
#ifndef _PYMOL_NOPY
      int n_thread = SettingGetGlobal_i(G, cSetting_max_threads);
      int multithread = SettingGetGlobal_i(G, cSetting_async_builds);

      if(multithread && (n_thread) && (stop - start) > 1) {
        int cnt = 0;

        /* must precalculate to avoid race-condition since this isn't
           mutexed yet and neighbors are needed by cartoons */
        this->getNeighborArray();

        for(a = start; a < stop; a++)
          if((a<I->NCSet) && I->CSet[a])
            cnt++;
        {
          CCoordSetUpdateThreadInfo *thread_info = pymol::malloc<CCoordSetUpdateThreadInfo>(cnt);
          if(thread_info) {
            cnt = 0;
            for(a = start; a < stop; a++) {
              if((a<I->NCSet) && I->CSet[a]) {
                thread_info[cnt].cs = I->CSet[a];
                thread_info[cnt].a = a;
                cnt++;
              }
            }
            ObjMolCoordSetUpdateSpawn(G, thread_info, n_thread, cnt);
            FreeP(thread_info);
          }

        }

      } else
#endif
      {                         /* single thread */
        for(a = start; a < stop; a++) {
          if((a<I->NCSet) && I->CSet[a] && (!G->Interrupt)) {
	    /* status bar */
            OrthoBusySlow(G, a, I->NCSet);
            PRINTFB(G, FB_ObjectMolecule, FB_Blather)
              " ObjectMolecule-DEBUG: updating representations for state %d of \"%s\".\n",
              a + 1, I->Name ENDFB(G);
            I->CSet[a]->update(a);
          }
        }
      }
    }
  } /* end block */

  PRINTFD(G, FB_ObjectMolecule)
    " ObjectMolecule: updates complete for object %s.\n", I->Name ENDFD;
}

/*========================================================================*/
void ObjectMolecule::invalidate(cRep_t rep, cRepInv_t level, int state)
{
  auto I = this;
  int a;
  PRINTFD(I->G, FB_ObjectMolecule)
    " %s: entered. rep: %d level: %d\n", __func__, rep, level ENDFD;

  auto const level_actual = level;

  // Remove the "purge" bit
  level = static_cast<decltype(level)>(level & ~cRepInvPurgeMask);

  if(level >= cRepInvVisib) {
    I->RepVisCacheValid = false;
  }

  if (level >= cRepInvBondsNoNonbonded) {
    if (level < cRepInvBonds) {
      level = cRepInvBonds;
    } else {
      ObjectMoleculeUpdateNonbonded(I);
    }
  }

  if(level >= cRepInvBonds) {
    this->Neighbor.reset();
    if(I->Sculpt) {
      DeleteP(I->Sculpt);
    }
    if(level >= cRepInvAtoms) {
      SelectorUpdateObjectSele(I->G, I);
    }
  }
  PRINTFD(I->G, FB_ObjectMolecule)
    " %s: invalidating representations...\n", __func__ ENDFD;

  if ( level >= cRepInvColor ) { 
    /* after label, this gets called, so we shouldn't invalidate types b/c PYMOL-317
       not sure if that is the exact level we should cut this off at
       1/28/13 BB: changed from >cRepInvText to >=cRepInvColor for invalidating surface
       colors during gadget invalidation */
    
    int start = 0;
    int stop = I->NCSet;

    if(state >= 0) {
      start = state;
      stop = state + 1;
    }
    if(stop > I->NCSet)
      stop = I->NCSet;
    for(a = start; a < stop; a++) {
      CoordSet *cset = 0;
      cset = I->CSet[a];
      if(cset) {
        cset->invalidateRep(rep, level_actual);
      }
    }
  }
  
  PRINTFD(I->G, FB_ObjectMolecule)
    " %s: leaving...\n", __func__ ENDFD;

}

void ObjectMoleculeInvalidateAtomType(ObjectMolecule *I, int state){
  CoordSet *cset = 0;
  int ai, atm;
  AtomInfoType *at;
  cset = I->CSet[state];
  if (state < 0){
    for (ai=0; ai < I->NAtom; ai++){
      at = &I->AtomInfo[ai];
      at->textType = 0;
    }
  } else {
    for (ai=0; ai < cset->NIndex; ai++){
      atm = cset->IdxToAtm[ai];
      if (atm>=0){
	at = &I->AtomInfo[ai];
	at->textType = 0;
      }
    }
  }
}

/*========================================================================*/
int ObjectMoleculeMoveAtom(ObjectMolecule * I, int state, int index, const float *v, int mode,
                           int log)
{
  int result = 0;
  PyMOLGlobals *G = I->G;
  CoordSet *cs;
  if(!(I->AtomInfo[index].protekted == 1)) {
    if(state < 0)
      state = 0;
    if(I->NCSet == 1)
      state = 0;
    state = state % I->NCSet;
    if((!I->CSet[state]) && (SettingGet_b(G, I->Setting.get(), NULL, cSetting_all_states)))
      state = 0;
    cs = I->CSet[state];
    if(cs) {
      result = CoordSetMoveAtom(I->CSet[state], index, v, mode);
      cs->invalidateRep(cRepAll, cRepInvCoord);
      ExecutiveUpdateCoordDepends(G, I);
    }
  }
  if(log) {
    OrthoLineType line, buffer;
    if(SettingGetGlobal_i(G, cSetting_logging)) {
      ObjectMoleculeGetAtomSele(I, index, buffer);
      sprintf(line, "cmd.translate_atom(\"%s\",%15.9f,%15.9f,%15.9f,%d,%d,%d)\n",
              buffer, v[0], v[1], v[2], state + 1, mode, 0);
      PLog(G, line, cPLog_no_flush);
    }
  }
  return (result);
}


/*========================================================================*/
int ObjectMoleculeMoveAtomLabel(ObjectMolecule * I, int state, int index, float *v, int log, float *diff)
{
  int result = 0;
  CoordSet *cs;
  if(!(I->AtomInfo[index].protekted == 1)) {
    if(state < 0)
      state = 0;
    if(I->NCSet == 1)
      state = 0;
    state = state % I->NCSet;
    if((!I->CSet[state])
       && (SettingGet_b(I->G, I->Setting.get(), NULL, cSetting_all_states)))
      state = 0;
    cs = I->CSet[state];
    if(cs) {
      result = CoordSetMoveAtomLabel(I->CSet[state], index, v, diff);
      cs->invalidateRep(cRepLabel, cRepInvCoord);
    }
  }
  return (result);
}

				 
/*========================================================================*/
int ObjectMoleculeInitBondPath(ObjectMolecule * I, ObjectMoleculeBPRec * bp)
{
  int a;
  bp->dist = pymol::malloc<int>(I->NAtom);
  bp->list = pymol::malloc<int>(I->NAtom);
  for(a = 0; a < I->NAtom; a++)
    bp->dist[a] = -1;
  bp->n_atom = 0;
  return 1;
}


/*========================================================================*/
int ObjectMoleculePurgeBondPath(ObjectMolecule * I, ObjectMoleculeBPRec * bp)
{
  FreeP(bp->dist);
  FreeP(bp->list);
  return 1;
}


/*========================================================================*/
int ObjectMoleculeGetBondPaths(ObjectMolecule * I, int atom,
                               int max, ObjectMoleculeBPRec * bp)
{
  /* returns list of bond counts from atom to all others 
     dist and list must be vla array pointers or NULL */

  int b_cnt = 0;

  /* reinitialize dist array (if we've done at least one pass) */

  for (int a = 0; a < bp->n_atom; a++)
    bp->dist[bp->list[a]] = -1;

  bp->n_atom = 0;
  bp->dist[atom] = 0;
  bp->list[bp->n_atom] = atom;
  bp->n_atom++;

  int cur = 0;
  while(1) {
    b_cnt++;
    if(b_cnt > max)
      break;

    unsigned n_cur = bp->n_atom - cur;

    /* iterate through all current atoms */

    if(!n_cur)
      break;
    while(n_cur--) {
      int const a1 = bp->list[cur++];
      for (auto const& neighbor : AtomNeighbors(I, a1)) {
        auto const a2 = neighbor.atm;
        if(bp->dist[a2] < 0) {  /* for each atom not yet sampled... */
          bp->dist[a2] = b_cnt;
          bp->list[bp->n_atom] = a2;
          bp->n_atom++;
        }
      }
    }
  }
  return (bp->n_atom);
}


/*========================================================================*/
int ***ObjectMoleculeGetBondPrint(ObjectMolecule * I, int max_bond, int max_type,
                                  int *dim)
{
  int a, b, i, c;
  int at1, at2;
  int ***result = NULL;
  ObjectMoleculeBPRec bp;

  dim[0] = max_type + 1;
  dim[1] = max_type + 1;
  dim[2] = max_bond + 1;

  result = (int ***) UtilArrayCalloc((unsigned int *) dim, 3, sizeof(int));

  ObjectMoleculeInitBondPath(I, &bp);
  for(a = 0; a < I->NAtom; a++) {
    at1 = I->AtomInfo[a].customType;
    if((at1 >= 0) && (at1 <= max_type)) {
      ObjectMoleculeGetBondPaths(I, a, max_bond, &bp);
      for(b = 0; b < bp.n_atom; b++) {
        i = bp.list[b];
        at2 = I->AtomInfo[i].customType;
        if((at2 >= 0) && (at2 <= max_type)) {
          c = bp.dist[i];
          result[at1][at2][c]++;
        }
      }
    }
  }
  ObjectMoleculePurgeBondPath(I, &bp);
  return (result);
}


/*========================================================================*/
float ObjectMoleculeGetAvgHBondVector(ObjectMolecule * I, int atom,
                                      int state, float *v, float *incoming)


/* computes average / optima hydrogen bonding vector for an atom */
{
  float result = 0.0;
  int vec_cnt = 0;
  float v_atom[3], v_neigh[3], v_diff[3], v_acc[3] = { 0.0, 0.0, 0.0 };
  int sp2_flag = false;
  int order;

  auto const* cs = I->getCoordSet(state);
  if(cs) {
    if (CoordSetGetAtomVertex(cs, atom, v_atom)) { /* atom exists in this C-set */
      for (auto const& neighbor : AtomNeighbors(I, atom)) {
        auto const a2 = neighbor.atm;
        order = I->Bond[neighbor.bond].order;
        if((order == 2) || (order == 4)) {
          sp2_flag = true;
        }

        if(I->AtomInfo[a2].protons != 1) {      /* ignore hydrogens */
          if(CoordSetGetAtomVertex(cs, a2, v_neigh)) {
            subtract3f(v_atom, v_neigh, v_diff);
            normalize3f(v_diff);
            add3f(v_diff, v_acc, v_acc);
            vec_cnt++;
          }
        }
      }
      if(vec_cnt) {
        result = (float) length3f(v_acc);
        result = result / vec_cnt;
        normalize23f(v_acc, v);
      } else {
        copy3f(v_acc, v);
      }

      if(incoming && (vec_cnt == 1) && (fabs(dot_product3f(v, incoming)) < 0.99F)) {
        /* if we know where the donor is, and the acceptor can
           potentially rotate the lone pair, then we should optimally
           orient the acceptor, if possible */
        AtomInfoType *ai = I->AtomInfo + atom;
        float v_perp[3];
        float v_tmp1[3], v_tmp2[3];
        if(((ai->protons == cAN_O) && (!sp2_flag)) ||   /* C-O-H */
           ((ai->protons == cAN_N) && (sp2_flag))) {    /* C=N-H */

          remove_component3f(incoming, v, v_perp);
          normalize3f(v_perp);
          scale3f(v, 0.333644F, v_tmp1);
          scale3f(v_perp, 0.942699F, v_tmp2);
          add3f(v_tmp1, v_tmp2, v_tmp2);
          subtract3f(v, v_tmp2, v);
          normalize3f(v);
        }
      }
    }
  }
  return (result);
}


/*========================================================================*/
int ObjectMoleculeGetAtomVertex(const ObjectMolecule * I, int state, int index, float *v)
{
  int result = 0;
  if(state < 0)
    state = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_state) - 1;
  if(state < 0)
    state = SceneGetState(I->G);
  if(I->NCSet == 1)
    state = 0;                  /* static singletons always active here it seems */
  state = state % I->NCSet;
  if((!I->CSet[state])
     && (SettingGet_b(I->G, I->Setting.get(), NULL, cSetting_all_states)))
    state = 0;
  if(I->CSet[state])
    result = CoordSetGetAtomVertex(I->CSet[state], index, v);

  return (result);
}


/*========================================================================*/
int ObjectMoleculeGetAtomTxfVertex(const ObjectMolecule * I, int state, int index, float *v)
{
  int result = 0;
  const CoordSet* cs = nullptr;
  if (I->DiscreteFlag){
    cs = I->DiscreteCSet[index];
  }
  if(state < 0){
    state = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_state) - 1;
  }
  if(state < 0)
    state = SceneGetState(I->G);
  if(I->NCSet == 1)
    state = 0;                  /* static singletons always active here it seems */
  state = state % I->NCSet;
  {
    if (!cs)
      cs = I->CSet[state];
    if((!cs) && (SettingGet_b(I->G, I->Setting.get(), NULL, cSetting_all_states))) {
      state = 0;
      cs = I->CSet[state];
    }
    if(cs) {
      result = CoordSetGetAtomTxfVertex(cs, index, v);
    }
  }
  return (result);
}

/*========================================================================*/
int ObjectMoleculeSetAtomVertex(ObjectMolecule * I, int state, int index, float *v)
{
  int result = 0;
  if(state < 0)
    state = SettingGet_i(I->G, NULL, I->Setting.get(), cSetting_state) - 1;
  if(state < 0)
    state = SceneGetState(I->G);
  if(I->NCSet == 1)
    state = 0;
  state = state % I->NCSet;
  if((!I->CSet[state])
     && (SettingGet_b(I->G, I->Setting.get(), NULL, cSetting_all_states)))
    state = 0;
  if(I->CSet[state])
    result = CoordSetSetAtomVertex(I->CSet[state], index, v);
  return (result);
}


/*========================================================================*/
void ObjectMolecule::render(RenderInfo * info)
{
  auto I = this;
  int state = info->state;
  const RenderPass pass = info->pass;
  CoordSet *cs;
  int pop_matrix = false;
  int use_matrices = SettingGet_i(I->G, I->Setting.get(), NULL, cSetting_matrix_mode);
  if(use_matrices<0) use_matrices = 0;
  PRINTFD(I->G, FB_ObjectMolecule)
    " ObjectMolecule: rendering %s pass %d...\n", I->Name, static_cast<int>(pass) ENDFD;

  ObjectPrepareContext(I, info);

  for(StateIterator iter(G, I->Setting.get(), state, I->NCSet); iter.next();) {
    cs = I->CSet[iter.state];
    if(cs) {
      if(use_matrices)
        pop_matrix = ObjectStatePushAndApplyMatrix(cs, info);
      cs->render(info);
      if(pop_matrix)
        ObjectStatePopMatrix(cs, info);
    }
  }
  PRINTFD(I->G, FB_ObjectMolecule)
    " ObjectMolecule: rendering complete for object %s.\n", I->Name ENDFD;
}


/*========================================================================*/
void ObjectMoleculeDummyUpdate(ObjectMolecule * I, int mode)
{
  switch (mode) {
  case cObjectMoleculeDummyOrigin:
    SceneOriginGet(I->G, I->CSet[0]->Coord.data());
    break;
  case cObjectMoleculeDummyCenter:
    SceneGetCenter(I->G, I->CSet[0]->Coord.data());
    break;
  }
}


/*========================================================================*/
ObjectMolecule *ObjectMoleculeDummyNew(PyMOLGlobals * G, int type)
{
  auto I = new ObjectMolecule(G, false);
  I->AtomInfo.resize(I->NAtom = 1);

  auto* cset = new CoordSet(G);
  cset->Obj = I;
  cset->setNIndex(1);
  cset->enumIndices();

  I->CSet.resize(I->NCSet = 1);
  I->CSet[0] = cset;

  return I;
}


/*========================================================================*/

ObjectMolecule::ObjectMolecule(PyMOLGlobals * G, int discreteFlag) : pymol::CObject(G)
{
  auto I = this;
  int a;
  I->type = cObjectMolecule;
  I->CSet = pymol::vla<CoordSet*>(10); /* auto-zero */
  I->DiscreteFlag = discreteFlag;
  if(I->DiscreteFlag) {         /* discrete objects don't share atoms between states */
    I->DiscreteAtmToIdx = pymol::vla<int>(0);
    I->DiscreteCSet = pymol::vla<CoordSet*>(0);
  } else {
    I->DiscreteAtmToIdx = NULL;
    I->DiscreteCSet = NULL;
  }
  I->AtomInfo = pymol::vla<AtomInfoType>(10);
  for(a = 0; a <= cUndoMask; a++) {
    I->UndoCoord[a] = NULL;
    I->UndoState[a] = -1;
  }
  I->UndoIter = 0;
}


/*========================================================================*/
void ObjectMoleculeCopyNoAlloc(const ObjectMolecule* obj, ObjectMolecule* I)
{
  PyMOLGlobals * G = const_cast<PyMOLGlobals*>(obj->G);

  int a;
  BondType *i0;
  const BondType *i1;
  (*I) = (*obj);
  I->Sculpt = NULL;
  I->Setting.reset(SettingCopyAll(G, obj->Setting.get(), nullptr));

  I->ViewElem = NULL;
  I->gridSlotSelIndicatorsCGO = NULL;

  for(a = 0; a <= cUndoMask; a++)
    I->UndoCoord[a] = NULL;
  I->CSet = pymol::vla<CoordSet*>(I->NCSet);   /* auto-zero */
  for(a = 0; a < I->NCSet; a++) {
    I->CSet[a] = CoordSetCopy(obj->CSet[a]);
    if (I->CSet[a])
      I->CSet[a]->Obj = I;
  }

  if(obj->CSTmpl)
    I->CSTmpl = CoordSetCopy(obj->CSTmpl);

  if (obj->DiscreteFlag){
    int sz = VLAGetSize(obj->DiscreteAtmToIdx);
    I->DiscreteAtmToIdx = VLACopy2(obj->DiscreteAtmToIdx);
    I->DiscreteCSet = pymol::vla<CoordSet*>(sz);
    I->updateAtmToIdx();
  }
  I->Bond = pymol::vla<BondType>(I->NBond);
  i0 = I->Bond.data();
  i1 = obj->Bond.data();
  for(a = 0; a < I->NBond; a++) {
    AtomInfoBondCopy(G, i1++, i0++);
  }

  // unfortunately, the AtomInfoType is not copy-constructable yet
  // assert copy done correct
  if (I->AtomInfo.size() != obj->AtomInfo.size()) {
    throw "AtomInfo copy failed";
  }
  AtomInfoType *a0 = I->AtomInfo.data();
  const AtomInfoType* a1 = obj->AtomInfo.data();
  memset(a0, 0, sizeof(AtomInfoType) * I->NAtom);
  for(a = 0; a < I->NAtom; a++)
    AtomInfoCopy(G, a1++, a0++);
}

ObjectMolecule *ObjectMoleculeCopy(const ObjectMolecule * obj)
{
  PyMOLGlobals * G = const_cast<PyMOLGlobals*>(obj->G);
  auto I = new ObjectMolecule(G, obj->DiscreteFlag);
  ObjectMoleculeCopyNoAlloc(obj, I);
  return (I);

}

/*========================================================================*/
/**
 * Set the order of coordinate sets with an index array
 */
int ObjectMoleculeSetStateOrder(ObjectMolecule * I, int * order, int len) {
  int a;
  CoordSet ** csets = VLAlloc(CoordSet *, I->NCSet);

  ok_assert(1, len == I->NCSet);

  // invalidate
  I->invalidate(cRepAll, cRepInvAll, -1);

  // new coord set array
  for(a = 0; a < I->NCSet; a++) {
    int i = order[a];
    ok_assert(1, 0 <= i && i < I->NCSet);
    csets[a] = I->CSet[i];
  }

  VLAFreeP(I->CSet);
  I->CSet = pymol::vla_take_ownership(csets);

  return true;
ok_except1:
  ErrMessage(I->G, "ObjectMoleculeSetStateOrder", "failed");
  VLAFreeP(csets);
  return false;
}

/*========================================================================*/
ObjectMolecule::~ObjectMolecule()
{
  auto I = this;
  int a;
  SelectorPurgeObjectMembers(I->G, I);
  for(a = 0; a < I->NCSet; a++){
    if(I->CSet[a]) {
      delete I->CSet[a];
      I->CSet[a] = NULL;
    }
  }
  VLAFreeP(I->DiscreteAtmToIdx);
  VLAFreeP(I->DiscreteCSet);
  VLAFreeP(I->CSet);
  I->m_ciffile.reset(); // free data

  {
    int nAtom = I->NAtom;
    AtomInfoType *ai = I->AtomInfo.data();

    for(a = 0; a < nAtom; a++) {
      AtomInfoPurge(I->G, ai);
      ai++;
    }
    VLAFreeP(I->AtomInfo);
  }
  {
    int nBond = I->NBond;
    BondType *bi = I->Bond.data();

    for(a = 0; a < nBond; a++) {
      AtomInfoPurgeBond(I->G, bi);
      bi++;
    }
    VLAFreeP(I->Bond);
  }
  for(a = 0; a <= cUndoMask; a++)
    FreeP(I->UndoCoord[a]);
  if(I->Sculpt)
    DeleteP(I->Sculpt);
  delete I->CSTmpl;
}

/*========================================================================*/
ObjectMolecule *ObjectMoleculeReadPDBStr(PyMOLGlobals * G, ObjectMolecule * I,
                                         const char *PDBStr, int state, int discrete,
                                         char *pdb_name,
                                         const char **next_pdb, PDBInfoRec * pdb_info,
                                         int quiet, int *model_number)
{
  CoordSet *cset = NULL;
  pymol::vla<AtomInfoType> atInfo;
  int ok = true;
  int isNew = true;
  unsigned int nAtom = 0;
  const char *start, *restart = NULL;
  int repeatFlag = true;
  int successCnt = 0;
  unsigned int aic_mask = cAIC_PDBMask;

  SegIdent segi_override = "";  /* saved segi for corrupted NMR pdb files */

  start = PDBStr;
  while(repeatFlag) {
    repeatFlag = false;

    if(!I)
      isNew = true;
    else
      isNew = false;

    if(ok) {
      if(isNew) {
        I = (ObjectMolecule *) new ObjectMolecule(G, discrete);
	CHECKOK(ok, I);
	if (ok)
          std::swap(atInfo, I->AtomInfo);
        isNew = true;
      } else {
        atInfo = pymol::vla<AtomInfoType>(10);
	CHECKOK(ok, atInfo);
        isNew = false;
      }
      if(ok && isNew) {
        I->Color = AtomInfoUpdateAutoColor(G);

        if (pdb_info->variant == PDB_VARIANT_VDB ||
            pdb_info->variant == PDB_VARIANT_PQR) {
          // pqr files have no chain identifier by default
          // vdb files have same chain identifiers in all symmetry copies
          SettingSet(cSetting_retain_order, 1, I);
        }
      }
      if (ok)
	cset = ObjectMoleculePDBStr2CoordSet(G, start, &atInfo, &restart,
					     segi_override, pdb_name,
					     next_pdb, pdb_info, quiet, model_number);
      CHECKOK(ok, cset);
      if (ok){
	nAtom = cset->NIndex;
      }
    }
    if(ok && pdb_name && (*next_pdb)) {
      /* problematic scenario? */
    }

    /* include coordinate set */
    if(ok) {
      if(I->DiscreteFlag && atInfo) {
        unsigned int a;
        int fp1 = state + 1;
        AtomInfoType *ai = atInfo.data();
        for(a = 0; a < nAtom; a++) {
          (ai++)->discrete_state = fp1;
        }
      }

      cset->Obj = I;
      cset->enumIndices();
      cset->invalidateRep(cRepAll, cRepInvRep);
      if(isNew) {
        std::swap(I->AtomInfo, atInfo);
      } else {
        ok &= ObjectMoleculeMerge(I, std::move(atInfo), cset, true, aic_mask, true);
        /* NOTE: will release atInfo */
      }
      if(isNew)
        I->NAtom = nAtom;
      if(state < 0)
        state = I->NCSet;
      if(*model_number > 0) {
        if(SettingGetGlobal_b(G, cSetting_pdb_honor_model_number))
          state = *model_number - 1;
      }
      VLACheck(I->CSet, CoordSet *, state);
      CHECKOK(ok, I->CSet);
      if(ok){
	if(I->NCSet <= state)
	  I->NCSet = state + 1;
	delete I->CSet[state];
	I->CSet[state] = cset;
      }
      if(ok && isNew)
        ok &= ObjectMoleculeConnect(I, cset);
      if(ok && cset->Symmetry) {
        I->Symmetry.reset(new CSymmetry(*cset->Symmetry));
      }
      if (I->Symmetry) {
          /* check scale records */
          if(pdb_info &&
             pdb_info->scale.flag[0] &&
             pdb_info->scale.flag[1] && pdb_info->scale.flag[2]) {

            float *sca = pdb_info->scale.matrix;
            sca[15] = 1.0F;

            CoordSetInsureOrthogonal(G, cset, sca, &I->Symmetry->Crystal, quiet);
          }
      }
      SceneCountFrames(G);
      if (ok)
	ok &= ObjectMoleculeExtendIndices(I, state);
      if (ok)
	ok &= ObjectMoleculeSort(I);
      if (ok){
	ObjectMoleculeUpdateIDNumbers(I);
	ObjectMoleculeUpdateNonbonded(I);
	ObjectMoleculeAutoDisableAtomNameWildcard(I);
      }
      if(SettingGetGlobal_b(G, cSetting_pdb_hetatm_guess_valences)) {
        ObjectMoleculeGuessValences(I, state, NULL, NULL, false);
      }

      successCnt++;
      if(!quiet) {
        if(successCnt > 1) {
          if(successCnt == 2) {
            PRINTFB(G, FB_ObjectMolecule, FB_Actions)
              " %s: read MODEL %d\n", __func__, 1 ENDFB(G);
          }
          PRINTFB(G, FB_ObjectMolecule, FB_Actions)
            " %s: read MODEL %d\n", __func__, successCnt ENDFB(G);
        }
      }
    }
    if(restart) {
      repeatFlag = true;
      start = restart;
      state = state + 1;
    }
  }
  if (!ok && isNew){
    DeleteP(I);
  }
  return (I);
}


/*========================================================================*/
static
CoordSet *ObjectMoleculeMMDStr2CoordSet(PyMOLGlobals * G, const char *buffer,
                                        AtomInfoType ** atInfoPtr, const char **restart)
{
  const char *p;
  int nAtom, nBond;
  int a, c, bPart, bOrder;
  float *coord = NULL;
  CoordSet *cset = NULL;
  AtomInfoType *atInfo = NULL, *ai;
  char cc[MAXLINELEN];
  WordType title;

  float *f;
  BondType *ii, *bond = NULL;
  int ok = true;
  int auto_show = RepGetAutoShowMask(G);

  p = buffer;
  nAtom = 0;
  if(atInfoPtr)
    atInfo = *atInfoPtr;

  if(ok) {
    p = ncopy(cc, p, 6);  // could this be too small for large molecules with #atoms > 999999 ?
    if(sscanf(cc, "%d", &nAtom) != 1)
      ok = ErrMessage(G, "ReadMMDFile", "bad atom count");
  }

  if(ok) {
    coord = VLAlloc(float, 3 * nAtom);
    if(atInfo)
      VLACheck(atInfo, AtomInfoType, nAtom);
  }

  if(!atInfo) {
    ErrFatal(G, "MMDStr2CoordSet", "need atom information record!");
    /* failsafe for old version.. */
  }

  nBond = 0;
  if(ok) {
    bond = VLACalloc(BondType, 6 * nAtom);
  }

  p = ParseWordCopy(title, p, sizeof(WordType) - 1);
  UtilCleanStr(title);

  p = nextline(p);

  /* read coordinates and atom names */

  if(ok) {
    f = coord;
    ii = bond;
    for(a = 0; a < nAtom; a++) {
      ai = atInfo + a;

      ai->id = a + 1;
      ai->rank = a;
      if(ok) {
        p = ncopy(cc, p, 4);
        if(sscanf(cc, "%d", &ai->customType) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad atom type");
      }
      if(ok) {
        if(ai->customType <= 14)
          strcpy(ai->elem, "C");
        else if(ai->customType <= 23)
          strcpy(ai->elem, "O");
        else if(ai->customType <= 40)
          strcpy(ai->elem, "N");
        else if(ai->customType <= 48)
          strcpy(ai->elem, "H");
        else if(ai->customType <= 52)
          strcpy(ai->elem, "S");
        else if(ai->customType <= 53)
          strcpy(ai->elem, "P");
        else if(ai->customType <= 55)
          strcpy(ai->elem, "B");
        else if(ai->customType <= 56)
          strcpy(ai->elem, "F");
        else if(ai->customType <= 57)
          strcpy(ai->elem, "Cl");
        else if(ai->customType <= 58)
          strcpy(ai->elem, "Br");
        else if(ai->customType <= 59)
          strcpy(ai->elem, "I");
        else if(ai->customType <= 60)
          strcpy(ai->elem, "Si");
        else if(ai->customType <= 61)
          strcpy(ai->elem, "Du");
        else if(ai->customType <= 62)
          strcpy(ai->elem, "Z0");
        else if(ai->customType <= 63)
          strcpy(ai->elem, "Lp");
        else
          ai->elem[0] = 0;
        /* else strcpy(ai->elem,"?"); WLD 090305 -- guess instead. */
      }
      for(c = 0; c < 6; c++) {
        if(ok) {
          p = ncopy(cc, p, 8);
          if(sscanf(cc, "%d%d", &bPart, &bOrder) != 2)
            ok = ErrMessage(G, "ReadMMDFile", "bad bond record");
          else {
            if(bPart && bOrder && (a < (bPart - 1))) {
              nBond++;
              ii->index[0] = a;
              ii->index[1] = bPart - 1;
              ii->order = bOrder;
              ii++;
            }
          }
        }
      }
      if(ok) {
        p = ncopy(cc, p, 12);
        if(sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad coordinate");
      }
      if(ok) {
        p = ncopy(cc, p, 12);
        if(sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad coordinate");
      }
      if(ok) {
        p = ncopy(cc, p, 12);
        if(sscanf(cc, "%f", f++) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad coordinate");
      }
      if(ok) {
        p = nskip(p, 1);
        p = ncopy(cc, p, 5);
        ai->setResi(cc);

        // single-letter code
        p = nskip(p, 1);

        // chain
        p = ncopy(cc, p, 1);
        LexAssign(G, ai->chain, cc);
      }
      if(ok) {
        p = nskip(p, 4);
        p = ncopy(cc, p, 9);
        if(sscanf(cc, "%f", &ai->partialCharge) != 1)
          ok = ErrMessage(G, "ReadMMDFile", "bad charge");
      }
      if(ok) {
        p = nskip(p, 10);
        p = ncopy(cc, p, 3);
        UtilCleanStr(cc);
        LexAssign(G, ai->resn, cc);
        ai->hetatm = true;
      }

      ai->segi = 0;
      ai->alt[0] = 0;

      if(ok) {
        p = nskip(p, 2);
        p = ntrim(cc, p, 4);
        if(!cc[0]) {
          sprintf(cc, "%s%02d", ai->elem, a + 1);
        }
        ai->name = LexIdx(G, cc);

        ai->visRep = auto_show;
      }
      if(ok) {
        AtomInfoAssignParameters(G, ai);
        AtomInfoAssignColors(G, ai);
      }
      if(!ok)
        break;
      p = nextline(p);
    }
  }
  if(ok)
    VLASize(bond, BondType, nBond);
  if(ok) {
    cset = CoordSetNew(G);
    cset->NIndex = nAtom;
    cset->Coord = pymol::vla_take_ownership(coord);
    cset->NTmpBond = nBond;
    cset->TmpBond = pymol::vla_take_ownership(bond);
    strcpy(cset->Name, title);
  } else {
    VLAFreeP(bond);
    VLAFreeP(coord);
  }
  if(atInfoPtr)
    *atInfoPtr = atInfo;
  *restart = *p ? p : NULL;
  return (cset);
}

#ifdef _PYMOL_IP_EXTRAS
#endif

#ifdef _PYMOL_IP_PROPERTIES
#endif

void AtomInfoSettingGenerateSideEffects(PyMOLGlobals * G, ObjectMolecule *obj, int index, int id){
  switch(index){
  case cSetting_label_position:
  case cSetting_label_placement_offset:
  case cSetting_label_screen_point:
  case cSetting_label_relative_mode:
    obj->invalidate(cRepLabel, cRepInvCoord, -1);
  }
}

static int AtomInfoInOrder(PyMOLGlobals * G, const AtomInfoType * atom, int atom1, int atom2)
{
  return (AtomInfoCompare(G, atom + atom1, atom + atom2) <= 0);
}

static int AtomInfoInOrderIgnoreHet(PyMOLGlobals * G, const AtomInfoType * atom,
    int atom1, int atom2)
{
  return (AtomInfoCompareIgnoreHet(G, atom + atom1, atom + atom2) <= 0);
}

static int AtomInfoInOrigOrder(PyMOLGlobals * G, const AtomInfoType * atom,
    int atom1, int atom2)
{
  if(atom[atom1].rank == atom[atom2].rank)
    return (AtomInfoCompare(G, atom + atom1, atom + atom2) <= 0);
  return (atom[atom1].rank < atom[atom2].rank);
}

int *AtomInfoGetSortedIndex(PyMOLGlobals * G,
    const ObjectMolecule* obj,
    const AtomInfoType* rec, int n, int **outdex)
{
  int *index;
  int a;
  const CSetting* setting = nullptr;

  ok_assert(1, index = pymol::malloc<int>(n + 1));
  ok_assert(1, (*outdex) = pymol::malloc<int>(n + 1));

  if(obj && obj->DiscreteFlag) {
    for(a = 0; a < n; a++)
      index[a] = a;
  } else {
    if(obj)
      setting = obj->Setting.get();

    UtilSortIndexGlobals(G, n, rec, index, (UtilOrderFnGlobals *) (
        SettingGet_b(G, setting, NULL, cSetting_retain_order) ?
          AtomInfoInOrigOrder :
        SettingGet_b(G, setting, NULL, cSetting_pdb_hetatm_sort) ?
          AtomInfoInOrder :
          AtomInfoInOrderIgnoreHet));
  }

  for(a = 0; a < n; a++)
    (*outdex)[index[a]] = a;

  return index;

ok_except1:
  FreeP(index);
  return NULL;
}

/**
 * Set the size of the DiscreteAtmToIdx and DiscreteCSet VLAs and pad
 * them with -1/NULL if necessary.
 */
bool ObjectMolecule::setNDiscrete(int natom) {
  int n = VLAGetSize(DiscreteAtmToIdx);

  if (n == natom)
    return true;

  VLASize(DiscreteAtmToIdx, int, natom);
  VLASize(DiscreteCSet, CoordSet*, natom);

  if (!DiscreteAtmToIdx || !DiscreteCSet)
    return false;

  for (int i = n; i < natom; ++i) {
    DiscreteAtmToIdx[i] = -1;
    DiscreteCSet[i] = NULL;
  }

  return true;
}

/**
 * Update the AtmToIdx or DiscreteAtmToIdx/DiscreteCSet VLAs from the
 * IdxToAtm arrays
 */
bool ObjectMolecule::updateAtmToIdx() {
  if (DiscreteFlag) {
    ok_assert(1, setNDiscrete(NAtom));
  }

  for (int i = -1; i < NCSet; ++i) {
    CoordSet * cset = (i < 0) ? CSTmpl : CSet[i];

    if (!cset)
      continue;

    if (!DiscreteFlag) {
      cset->updateNonDiscreteAtmToIdx(NAtom);
    } else {
      for (int idx = 0; idx < cset->NIndex; ++idx) {
        int atm = cset->IdxToAtm[idx];
        assert(atm < NAtom);
        DiscreteAtmToIdx[atm] = idx;
        DiscreteCSet[atm] = cset;
        AtomInfo[atm].discrete_state = i + 1;
      }
    }
  }

  return true;
ok_except1:
  return false;
}

/**
 * Check if this atom has coordinates in any state
 */
bool ObjectMolecule::atomHasAnyCoordinates(size_t atm) const
{
  for (size_t i = 0; i < NCSet; ++i) {
    auto cset = CSet[i];
    if (cset && cset->atmToIdx(atm) != -1) {
      return true;
    }
  }

  return false;
}

pymol::CObject* ObjectMolecule::clone() const
{
  return ObjectMoleculeCopy(this);
}

AtomNeighbors::AtomNeighbors(const ObjectMolecule* I, int atm)
{
  auto* Neighbor = I->getNeighborArray();

  m_neighbor = Neighbor + Neighbor[atm];
}
