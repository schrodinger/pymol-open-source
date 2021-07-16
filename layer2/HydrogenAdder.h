#pragma once

struct ObjectMolecule;
struct CoordSet;

int ObjectMoleculeAddSeleHydrogensRefactored(ObjectMolecule * I, int sele, int state);

int ObjectMoleculeSetMissingNeighborCoords(
    ObjectMolecule* I, CoordSet* cs, unsigned atm,
    bool h_fix=false);
