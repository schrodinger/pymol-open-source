#pragma once

#include "Executive.h"
#include "Result.h"

#ifndef _PYMOL_NO_PY

pymol::Result<> ExecutiveLoadObject(PyMOLGlobals* G,
    const char* oname, PyObject* model, int frame, int type, int finish,
    int discrete, int quiet, int zoom);

pymol::Result<> ExecutiveSetRawAlignment(PyMOLGlobals* G,
    pymol::zstring_view alnname, PyObject* raw, pymol::zstring_view guidename,
    int state, int quiet);

pymol::Result<float> ExecutiveFitPairs(
    PyMOLGlobals* G, PyObject* list, int quiet);

#endif //_PYMOL_NO_PY
