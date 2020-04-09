#pragma once

#include "Executive.h"
#include "Result.h"

#ifndef _PYMOL_NO_PY

pymol::Result<> ExecutiveLoadObject(PyMOLGlobals* G,
    const char* oname, PyObject* model, int frame, int type, int finish,
    int discrete, int quiet, int zoom);

#endif //_PYMOL_NO_PY
