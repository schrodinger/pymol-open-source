/*
 * PConvArgsToPyList
 * PConvArgsFromPyList
 *
 * Templated convenience function to create a Python list from
 * arbitrary number and vice versa
 *
 */

#include "os_python.h"
#include "PyMOLGlobals.h"
#include "PConv.h"

/*
 * PConvArgsToPyList
 */

inline void _PConvArgsToPyList_SetItem(PyObject * list, int i) {
}

template <typename T, typename... Args>
inline void _PConvArgsToPyList_SetItem(PyObject * list, int i, const T& value, const Args&... rest) {
  PyList_SET_ITEM(list, i, PConvToPyObject(value));
  _PConvArgsToPyList_SetItem(list, i + 1, rest...);
}

template <typename... Args>
inline PyObject * PConvArgsToPyList(const Args&... args) {
  PyObject * list = PyList_New(sizeof...(Args));
  _PConvArgsToPyList_SetItem(list, 0, args...);
  return list;
}

/*
 * PConvArgsFromPyList
 */

inline bool _PConvArgsFromPyList_GetItem(PyMOLGlobals * G, PyObject * list, int len, int i) {
  return (len == i);
}

template <typename T, typename... Args>
inline bool _PConvArgsFromPyList_GetItem(PyMOLGlobals * G, PyObject * list, int len, int i, T& out, Args&... rest) {
  if (len <= i)
    return false;

  PConvFromPyObject(G, PyList_GetItem(list, i), out);
  return _PConvArgsFromPyList_GetItem(G, list, len, i + 1, rest...);
}

template <typename... Args>
inline bool PConvArgsFromPyList(PyMOLGlobals * G, PyObject * list, Args&... args) {
  int len = PyList_Size(list);
  return _PConvArgsFromPyList_GetItem(G, list, len, 0, args...);
}

