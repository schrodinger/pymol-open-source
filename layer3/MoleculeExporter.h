/*
 * Molecular file formats export
 *
 * (c) 2016 Schrodinger, Inc.
 */

#include "os_std.h"
#include "unique_vla_ptr.h"

#include "PyMOLGlobals.h"

unique_vla_ptr<char> MoleculeExporterGetStr(PyMOLGlobals * G,
    const char *format,
    const char *sele="all",
    int state=-2, // current (-1 in Python API)
    const char *ref_object="",
    int ref_state=-1,
    int multi=-1,
    bool quiet=true);
