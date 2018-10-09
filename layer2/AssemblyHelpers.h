/*
 * Biological Assembly Helpers
 *
 * (c) 2017 Schrodinger, Inc.
 */

#pragma once

#include <set>

#include "AtomInfo.h"
#include "CoordSet.h"
#include "Lex.h"
#include "ObjectMolecule.h"

CoordSet * CoordSetCopyFilterChains(
    const CoordSet * other,
    const AtomInfoType * atInfo,
    const std::set<lexborrow_t> & chains_set);

void ObjectMoleculeSetAssemblyCSets(
    ObjectMolecule * I,
    CoordSet ** assembly_csets);
