/*
 * Molecular file formats export
 *
 * (c) 2016 Schrodinger, Inc.
 */

#include "os_python.h"
#include "os_std.h"
#include "vla.h"

#include "PyMOLGlobals.h"

pymol::vla<char> MoleculeExporterGetStr(PyMOLGlobals * G,
    const char *format,
    const char *sele="all",
    int state = cStateCurrent, // current (-1 in Python API)
    const char *ref_object="",
    int ref_state = cStateAll,
    int multi=-1,
    bool quiet=true);

PyObject *MoleculeExporterGetPyBonds(PyMOLGlobals * G,
    const char *selection, int state);
