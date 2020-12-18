/*
 * Biological Assembly Helpers
 *
 * (c) 2017 Schrodinger, Inc.
 */

#include <vector>

#include "os_std.h"

#include "AssemblyHelpers.h"
#include "MemoryDebug.h"

/**
 * Create a coordset for a segi (chain) selection
 */
CoordSet * CoordSetCopyFilterChains(
    const CoordSet * other,
    const AtomInfoType * atInfo,
    const std::set<lexborrow_t> & chains_set)
{
  std::vector<int> idxmap;
  idxmap.reserve(other->NIndex);

  for (int idx = 0; idx < other->NIndex; ++idx)
    if (chains_set.count(atInfo[other->IdxToAtm[idx]].segi) > 0)
      idxmap.push_back(idx);

  CoordSet* cset = new CoordSet(other->G);

  cset->setNIndex(idxmap.size());
  cset->Obj = other->Obj;

  for (int idx = 0; idx < cset->NIndex; ++idx) {
    cset->IdxToAtm[idx] = other->IdxToAtm[idxmap[idx]];
    copy3f(other->coordPtr(idxmap[idx]), cset->coordPtr(idx));
  }

  return cset;
}

/**
 * Replace coordinate sets and set all_states
 */
void ObjectMoleculeSetAssemblyCSets(
    ObjectMolecule * I,
    CoordSet ** assembly_csets)
{
  if (!assembly_csets)
    return;

  if (I->DiscreteFlag) {
    printf("error/TODO: can't make discrete assembly\n");
    return;
  }

  // remove asymetric unit coordinate sets
  for (int i = 0; i < I->NCSet; ++i)
    delete I->CSet[i];

  VLAFreeP(I->CSet);

  // get assembly coordinate sets into ObjectMolecule
  I->CSet = pymol::vla_take_ownership(assembly_csets);
  I->NCSet = VLAGetSize(assembly_csets);
  I->updateAtmToIdx();

  // all_states for multi-model assembly
  if (I->NCSet > 1) {
    SettingSet(cSetting_all_states, true, I);
  }
}
