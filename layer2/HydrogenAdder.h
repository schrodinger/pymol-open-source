#pragma once

#include "ObjectMolecule.h"
#include "CoordSet.h"

int ObjectMoleculeAddSeleHydrogensRefactored(ObjectMolecule * I, int sele, int state);

int ObjectMoleculeSetMissingNeighborCoords(
    ObjectMolecule* I, CoordSet* cs, unsigned atm,
    bool h_fix=false);
