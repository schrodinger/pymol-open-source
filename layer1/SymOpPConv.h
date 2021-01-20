/**
 * @file Python serialization of SymOp
 *
 * (c) Schrodinger, Inc.
 */

#pragma once

#include "os_python.h"

struct PyMOLGlobals;

namespace pymol
{
struct SymOp;
}

PyObject* PConvToPyObject(pymol::SymOp const&);
bool PConvFromPyObject(PyMOLGlobals*, PyObject*, pymol::SymOp&);
