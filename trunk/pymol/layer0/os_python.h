

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

#include <string.h>

#if PY_MAJOR_VERSION >= 3
# define PyInt_Check            PyLong_Check
# define PyInt_FromLong         PyLong_FromLong
# define PyInt_AsLong           PyLong_AsLong
# define PyInt_AS_LONG          PyLong_AS_LONG

# define PyNumber_Int           PyNumber_Long

# define PyString_Check                 PyUnicode_Check
# define PyString_Size                  PyUnicode_GetLength
# define PyString_GET_SIZE              PyUnicode_GetLength
# define PyString_FromString            PyUnicode_FromString
# define PyString_FromStringAndSize     PyUnicode_FromStringAndSize
# define PyString_InternFromString      PyUnicode_InternFromString
# define PyString_AsString              PyUnicode_AsUTF8
# define PyString_AS_STRING             PyUnicode_AsUTF8

# define PyCObject_AsVoidPtr(capsule)   PyCapsule_GetPointer(capsule, "name")
# define PyCObject_FromVoidPtr(p, d)    PyCapsule_New(p, "name", (PyCapsule_Destructor) d)
# define PyCObject_Check                PyCapsule_CheckExact

# define PyEval_EvalCode(o, ...)        PyEval_EvalCode((PyObject*)o, __VA_ARGS__)

# define Py_TPFLAGS_HAVE_ITER   0
#endif

/*
 * For compatibility with the pickletools, this type represents
 * an optionally owned C string and has to be returned by value.
 */
class SomeString {
  const char * m_str;
  mutable int m_length;
public:
  SomeString(const char * s, int len=-1) : m_str(s), m_length(len) {}
  inline const char * data()    const { return m_str; }
  inline const char * c_str()   const { return m_str; }
  operator const char * ()      const { return m_str; } // allows assignment to std::string
  inline size_t length()        const {
    if (m_length == -1) {
      m_length = m_str ? strlen(m_str) : 0;
    }
    return m_length;
  }
};

inline SomeString PyString_AsSomeString(PyObject * o) {
  return PyString_AsString(o);
}

inline SomeString PyBytes_AsSomeString(PyObject * o) {
  return SomeString(PyBytes_AsString(o), PyBytes_Size(o));
}

#endif

#include "os_predef.h"

#define PYOBJECT_CALLMETHOD(o, m, ...) PyObject_CallMethod(o, (char*)m, (char*)__VA_ARGS__)
#define PYOBJECT_CALLFUNCTION(o, ...) PyObject_CallFunction(o, (char*)__VA_ARGS__)

#endif
