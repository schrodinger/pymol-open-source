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
  for (atm++; atm < cs->NAtIndex; atm++) {
    idx = cs->atmToIdx(atm);

    if(idx < 0)
      continue;

    return true;
  }

  return false;
}

/*========================================================================*/
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

/*
 * advance the internal state to the next atom, return false if there is no
 * next atom
 */
bool SeleAtomIterator::next() {
  CSelector *I = G->Selector;

  while ((++a) < I->NAtom) {
    atm = I->Table[a].atom;
    obj = I->Obj[I->Table[a].model];

    if(!SelectorIsMember(G, obj->AtomInfo[atm].selEntry, sele))
      continue;

    return true;
  }

  return false;
}

/*========================================================================*/
/*
 * Quasi constructor, call this if `SeleCoordIterator` has been constructed
 * with the default constructor.
 */
void SeleCoordIterator::init(PyMOLGlobals * G_, int sele_, int state_) {
  G = G_;
  statearg = state_;

  // current state (use -3 for "effective" state)
  if (statearg == -2) {
    statearg = SettingGetGlobal_i(G, cSetting_state) - 1;
  }

  // safety check
  if (statearg < -1) {
    statearg = -3;
  }

  SelectorUpdateTable(G, statearg, sele_);
  setPerObject(false);
  reset();
}

void SeleCoordIterator::reset() {
  a = cNDummyAtoms - 1;
  state = statearg;
  prev_obj = NULL;

  if (isMultistate()) {
    state = 0;
    statemax = 0;
  }
}

/*
 * advance the internal state to the next atom, return false if there is no
 * next atom
 */
bool SeleCoordIterator::next() {
  CSelector *I = G->Selector;

  for (a++; a < I->NAtom; a++) {
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
    } else if (statearg == -3 && obj != prev_obj) {
      // "effective" state (no support here for settings all_states=1 or state=0)
      state = std::max(0, obj->getState());
      prev_obj = obj;
    }

    if(state >= obj->NCSet || !(cs = obj->CSet[state]))
      continue;

    atm = I->Table[a].atom;
    idx = cs->atmToIdx(atm);

    if(idx < 0)
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
