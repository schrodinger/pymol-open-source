#pragma once

#include "pymol/zstring_view.h"

struct PyMOLGlobals;

void CheckGLErrorOK(PyMOLGlobals* G, pymol::zstring_view errString);
