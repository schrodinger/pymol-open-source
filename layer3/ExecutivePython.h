#pragma once

#ifndef _PYMOL_NO_PY

#include "os_python.h"
#include "Result.h"
#include "pymol/zstring_view.h"

struct ObjectAlignment;

pymol::Result<> ExecutiveLoadObject(PyMOLGlobals* G,
    const char* oname, PyObject* model, int frame, int type, int finish,
    int discrete, int quiet, int zoom);

pymol::Result<> ExecutiveSetRawAlignment(PyMOLGlobals* G,
    pymol::zstring_view alnname, PyObject* raw, pymol::zstring_view guidename,
    int state, int quiet);

pymol::Result<float> ExecutiveFitPairs(
    PyMOLGlobals* G, PyObject* list, int quiet);

/**
 * @param name name of alignment object
 * @param active_only only consider active alignments
 * @param state state of alignment object
 * @return a list of lists of (object, index) tuples containing the
 * raw per-atom alignment relationships
 */
pymol::Result<PyObject*> ExecutiveGetRawAlignment(PyMOLGlobals* G,
    pymol::null_safe_zstring_view name, bool active_only, int state);

#endif //_PYMOL_NO_PY
