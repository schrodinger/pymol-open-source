/**
 * @file
 * Isosurface with VTKm
 *
 * (c) 2020 Schrodinger, Inc.
 */

#pragma once

#include "PyMOLEnums.h"
#include "vla.h"

struct PyMOLGlobals;
struct Isofield;
class CarveHelper;

int ContourSurfVolume(PyMOLGlobals* G, Isofield* field, float level,
    pymol::vla<int>& num,    //
    pymol::vla<float>& vert, //
    const int* range,        //
    cIsosurfaceMode mode,    //
    const CarveHelper*,      //
    cIsosurfaceSide side);
