/*
 * Atom iterators
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include <algorithm>

#include "AtomIterators.h"
#include "Selector.h"
#include "SelectorDef.h"

/*========================================================================*/
bool CoordSetAtomIterator::next() {
  for (++atm; atm < obj->NAtom; ++atm) {
    idx = cs->atmToIdx(atm);

    if(idx < 0)
      continue;

    return true;
  }

  return false;
}

/*========================================================================*/
/**
 * Make a new (owned) atom selection and update the global selector table for it.
 */
SeleAtomIterator::SeleAtomIterator(PyMOLGlobals * G_, const char * sele_) {
  G = G_;

  stmp = new char[1024];
  SelectorGetTmp(G, (char*) sele_, stmp);
  sele = SelectorIndexByName(G, stmp);

  SelectorUpdateTable(G, cSelectorUpdateTableAllStates, -1);

  reset();
}

SeleAtomIterator::~SeleAtomIterator() {
  if (stmp) {
    SelectorFreeTmp(G, stmp);
    delete[] stmp;
  }
}

void SeleAtomIterator::reset() {
  a = cNDummyAtoms - 1;
}

/**
 * advance the internal state to the next atom, return false if there is no
 * next atom
 */
bool SeleAtomIterator::next() {
  CSelector *I = G->Selector;

  while ((++a) < I->Table.size()) {
    atm = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];

    if(!SelectorIsMember(G, obj->AtomInfo[atm].selEntry, sele))
      continue;

    return true;
  }

  return false;
}

/*========================================================================*/
/**
 * @param sele_ Atom selection to iterate over
 * @param state_ Object state to iterate over (can be current (-2) or all (-1))
 * @param update_table If true, then update the table (once) and do not call
 * `SelectorIsMember` during iteration, assuming that the table stays valid and
 * contains exactly the selected atoms. If false, then assume the table is
 * up-to-date with a selection different to `sele_` (e.g. with all atoms) and
 * ::SelectorIsMember needs to be called during iteration.
 */
SeleCoordIterator::SeleCoordIterator(
    PyMOLGlobals* G_, int sele_, int state_, bool update_table)
{
  G = G_;
  statearg = state_;

  // current state (use -3 for "effective" state)
  if (statearg == cStateCurrent) {
    statearg = SettingGetGlobal_i(G, cSetting_state) - 1;
  }

  // safety check
  if (statearg < cStateAll) {
    statearg = cSelectorUpdateTableEffectiveStates;
  }

  if (update_table) {
    SelectorUpdateTable(G, statearg, sele_);
  } else {
    sele = sele_;
  }

  setPerObject(false);
  reset();
}

void SeleCoordIterator::reset() {
  a = cNDummyAtoms - 1;
  state = statearg;
  prev_obj = nullptr;
  cs = nullptr;

  if (isMultistate()) {
    state = 0;
    statemax = 0;
  }
}

/**
 * advance the internal state to the next atom, return false if there is no
 * next atom
 */
bool SeleCoordIterator::next() {
  CSelector *I = G->Selector;

  for (a++; a < I->Table.size(); a++) {
    obj = I->Obj[I->Table[a].model];

    if (isMultistate()) {
      if (isPerObject()) {
        if (obj != prev_obj) {
          if (nextStateInPrevObject())
            continue;

          // first cs of next object
          prev_obj = obj;
          state = 0;
        }
      } else if(statemax < obj->NCSet) {
        statemax = obj->NCSet;
      }
    } else if (statearg == cSelectorUpdateTableEffectiveStates &&
               obj != prev_obj) {
      // "effective" state (no support here for settings all_states=1 or state=0)
      state = std::max(0, obj->getCurrentState());
      prev_obj = obj;
    }

    if(state >= obj->NCSet || !(cs = obj->CSet[state]))
      continue;

    atm = I->Table[a].atom;
    idx = cs->atmToIdx(atm);

    if(idx < 0)
      continue;

    if (sele > 0 && !SelectorIsMember(G, getAtomInfo()->selEntry, sele))
      continue;

    return true;
  }

  if (isMultistate()) {
    if (isPerObject()) {
      if (nextStateInPrevObject())
        return next();
    } else if ((++state) < statemax) {
      a = cNDummyAtoms - 1;
      return next();
    }
  }

  return false;
}

// vi:sw=2:expandtab
