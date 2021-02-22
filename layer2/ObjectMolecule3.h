#pragma once

#include "ObjectMolecule.h"

void ObjectMoleculePBCUnwrap(ObjectMolecule&, bool bymol = true);
void ObjectMoleculePBCWrap(ObjectMolecule&, float const* center = nullptr);
