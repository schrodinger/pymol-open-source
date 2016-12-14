#ifndef _H_os_numpy
#define _H_os_numpy

#include "os_python.h"

#ifdef _PYMOL_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif

#endif
