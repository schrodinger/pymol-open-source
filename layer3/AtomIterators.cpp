/*
 * Atom iterators
 *
 * (c) 2014 Schrodinger, Inc.
 */

#include "AtomIterators.h"
#include "Selector.h"

bool CoordSetAtomIterator::next() {
  for (atm++; atm < cs->NAtIndex; atm++) {
    idx = cs->atmToIdx(atm);

    if(idx < 0)
      continue;

    return true;
  }

  return false;
}

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

// vi:sw=2:expandtab
