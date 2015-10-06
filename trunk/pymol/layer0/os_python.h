

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_os_python
#define _H_os_python

#ifdef PYMOL_ICC_LINUX
#include"/usr/include/bits/types.h"
#endif

#ifdef _PYMOL_NOPY
typedef int PyObject;
#undef _PYMOL_NUMPY
#else

// Python.h will redefine those, undef to avoid compiler warning
#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE

#include"Python.h"
#include<pythread.h>

#ifdef _PYMOL_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#endif

#include <string.h>

/*
 * For compatibility with the pickletools, this type represents
 * an optionally owned C string and has to be returned by value.
 */
class SomeString {
  const char * m_str;
public:
  SomeString(const char * s) : m_str(s) {}
  inline const char * data()    const { return m_str; }
  inline const char * c_str()   const { return m_str; }
  operator const char * ()      const { return m_str; } // allows assignment to std::string
  inline size_t length()        const { return m_str ? strlen(m_str) : 0; }
};

inline SomeString PyString_AsSomeString(PyObject * o) {
  return PyString_AsString(o);
}

#endif

#include "os_predef.h"

#define PYOBJECT_CALLMETHOD(o, m, ...) PyObject_CallMethod(o, (char*)m, (char*)__VA_ARGS__)
#define PYOBJECT_CALLFUNCTION(o, ...) PyObject_CallFunction(o, (char*)__VA_ARGS__)

#endif
