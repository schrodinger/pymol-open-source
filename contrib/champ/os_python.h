/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

#include"Python.h"

# define PyInt_Check            PyLong_Check
# define PyInt_FromLong         PyLong_FromLong
# define PyInt_AsLong           PyLong_AsLong
# define PyInt_AS_LONG          PyLong_AS_LONG

# define PyNumber_Int           PyNumber_Long

# define PyString_Check                 PyUnicode_Check
# define PyString_Size                  PyUnicode_GetLength
# define PyString_FromString            PyUnicode_FromString
# define PyString_FromStringAndSize     PyUnicode_FromStringAndSize
# define PyString_InternFromString      PyUnicode_InternFromString
# define PyString_AsString              PyUnicode_AsUTF8
# define PyString_AS_STRING             PyUnicode_AsUTF8

# define PyCObject_AsVoidPtr(capsule)   PyCapsule_GetPointer(capsule, NULL)
# define PyCObject_FromVoidPtr(p, d)    Error_Function_not_available_PyCapsule_New_destructor_differs()
# define PyCObject_Check                PyCapsule_CheckExact

# define PyEval_EvalCode(o, ...)        PyEval_EvalCode((PyObject*)o, __VA_ARGS__)

#endif



