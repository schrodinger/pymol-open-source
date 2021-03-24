

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

#include "os_predef.h"

#include <memory>

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

/**
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

namespace pymol {
/**
 * Destruction policy for unique_ptr<PyObject, pymol::pyobject_delete>
 *
 * Must only be used if the GIL is guaranteed when operator() is called.
 */
struct pyobject_delete {
  void operator()(PyObject* o) const { Py_DECREF(o); }
};

/**
 * Destruction policy for unique_ptr<PyObject, pymol::pyobject_delete_auto_gil>
 * Does not require GIL to be held.
 */
struct pyobject_delete_auto_gil {
  void operator()(PyObject* o) const
  {
    if (o) {
      auto gstate = PyGILState_Ensure();
      Py_DECREF(o);
      PyGILState_Release(gstate);
    }
  }
};

/**
 * RAII helper to ensure the Python GIL
 */
class GIL_Ensure
{
  PyGILState_STATE state;

public:
  GIL_Ensure();
  ~GIL_Ensure();
};
} // namespace pymol

namespace std
{
/**
 * Destruction policy which ensures the GIL before operator() is called.
 */
template <> struct default_delete<PyObject> {
  void operator()(PyObject* o) const
  {
    pymol::GIL_Ensure gil;
    Py_DECREF(o);
  }
};
} // namespace std

/**
 * Unique pointer which must only be used if the GIL is guaranteed when it goes
 * out of scope or is reset.
 */
using unique_PyObject_ptr = std::unique_ptr<PyObject, pymol::pyobject_delete>;

#endif

#define PYOBJECT_CALLMETHOD PyObject_CallMethod
#define PYOBJECT_CALLFUNCTION PyObject_CallFunction

#endif
