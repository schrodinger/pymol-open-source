/*
 * PConvArgsToPyList
 * PConvArgsFromPyList
 *
 * Templated convenience function to create a Python list from
 * 2 to 6 (arbitrary number with C++11) arguments and vice versa
 *
 * TODO: remove non-variadic part when full C++11 supports becomes available
 */

#include "os_python.h"
#include "PyMOLGlobals.h"
#include "PConv.h"

/*
 * PConvArgsToPyList
 */

#if __cplusplus >= 201103L
// variadic templates

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

#else
// no variadic templates

template <class A, class B>
inline PyObject * PConvArgsToPyList(const A& a, const B& b) {
  PyObject * obj = PyList_New(2);
  PyList_SET_ITEM(obj, 0, PConvToPyObject(a));
  PyList_SET_ITEM(obj, 1, PConvToPyObject(b));
  return obj;
}

template <class A, class B, class C>
inline PyObject * PConvArgsToPyList(const A& a, const B& b, const C& c) {
  PyObject * obj = PyList_New(3);
  PyList_SET_ITEM(obj, 0, PConvToPyObject(a));
  PyList_SET_ITEM(obj, 1, PConvToPyObject(b));
  PyList_SET_ITEM(obj, 2, PConvToPyObject(c));
  return obj;
}

template <class A, class B, class C, class D>
inline PyObject * PConvArgsToPyList(const A& a, const B& b, const C& c, const D& d) {
  PyObject * obj = PyList_New(4);
  PyList_SET_ITEM(obj, 0, PConvToPyObject(a));
  PyList_SET_ITEM(obj, 1, PConvToPyObject(b));
  PyList_SET_ITEM(obj, 2, PConvToPyObject(c));
  PyList_SET_ITEM(obj, 3, PConvToPyObject(d));
  return obj;
}

template <class A, class B, class C, class D, class E>
inline PyObject * PConvArgsToPyList(const A& a, const B& b, const C& c, const D& d, const E& e) {
  PyObject * obj = PyList_New(5);
  PyList_SET_ITEM(obj, 0, PConvToPyObject(a));
  PyList_SET_ITEM(obj, 1, PConvToPyObject(b));
  PyList_SET_ITEM(obj, 2, PConvToPyObject(c));
  PyList_SET_ITEM(obj, 3, PConvToPyObject(d));
  PyList_SET_ITEM(obj, 4, PConvToPyObject(e));
  return obj;
}

template <class A, class B, class C, class D, class E, class F>
inline PyObject * PConvArgsToPyList(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f) {
  PyObject * obj = PyList_New(6);
  PyList_SET_ITEM(obj, 0, PConvToPyObject(a));
  PyList_SET_ITEM(obj, 1, PConvToPyObject(b));
  PyList_SET_ITEM(obj, 2, PConvToPyObject(c));
  PyList_SET_ITEM(obj, 3, PConvToPyObject(d));
  PyList_SET_ITEM(obj, 4, PConvToPyObject(e));
  PyList_SET_ITEM(obj, 5, PConvToPyObject(f));
  return obj;
}

#endif

/*
 * PConvArgsFromPyList
 */

#if __cplusplus >= 201103L
// variadic templates

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

#else
// no variadic templates

template <class A, class B>
inline bool PConvArgsFromPyList(PyMOLGlobals * G, PyObject * obj, A& a, B& b) {
  int len = PyList_Size(obj);
  if (len > 0) PConvFromPyObject(G, PyList_GetItem(obj, 0), a);
  if (len > 1) PConvFromPyObject(G, PyList_GetItem(obj, 1), b);
  return true;
}

template <class A, class B, class C>
inline bool PConvArgsFromPyList(PyMOLGlobals * G, PyObject * obj, A& a, B& b, C& c) {
  int len = PyList_Size(obj);
  if (len > 0) PConvFromPyObject(G, PyList_GetItem(obj, 0), a);
  if (len > 1) PConvFromPyObject(G, PyList_GetItem(obj, 1), b);
  if (len > 2) PConvFromPyObject(G, PyList_GetItem(obj, 2), c);
  return true;
}

template <class A, class B, class C, class D>
inline bool PConvArgsFromPyList(PyMOLGlobals * G, PyObject * obj, A& a, B& b, C& c, D& d) {
  int len = PyList_Size(obj);
  if (len > 0) PConvFromPyObject(G, PyList_GetItem(obj, 0), a);
  if (len > 1) PConvFromPyObject(G, PyList_GetItem(obj, 1), b);
  if (len > 2) PConvFromPyObject(G, PyList_GetItem(obj, 2), c);
  if (len > 3) PConvFromPyObject(G, PyList_GetItem(obj, 3), d);
  return true;
}

template <class A, class B, class C, class D, class E>
inline bool PConvArgsFromPyList(PyMOLGlobals * G, PyObject * obj, A& a, B& b, C& c, D& d, E& e) {
  int len = PyList_Size(obj);
  if (len > 0) PConvFromPyObject(G, PyList_GetItem(obj, 0), a);
  if (len > 1) PConvFromPyObject(G, PyList_GetItem(obj, 1), b);
  if (len > 2) PConvFromPyObject(G, PyList_GetItem(obj, 2), c);
  if (len > 3) PConvFromPyObject(G, PyList_GetItem(obj, 3), d);
  if (len > 4) PConvFromPyObject(G, PyList_GetItem(obj, 4), e);
  return true;
}

template <class A, class B, class C, class D, class E, class F>
inline bool PConvArgsFromPyList(PyMOLGlobals * G, PyObject * obj, A& a, B& b, C& c, D& d, E& e, F& f) {
  int len = PyList_Size(obj);
  if (len > 0) PConvFromPyObject(G, PyList_GetItem(obj, 0), a);
  if (len > 1) PConvFromPyObject(G, PyList_GetItem(obj, 1), b);
  if (len > 2) PConvFromPyObject(G, PyList_GetItem(obj, 2), c);
  if (len > 3) PConvFromPyObject(G, PyList_GetItem(obj, 3), d);
  if (len > 4) PConvFromPyObject(G, PyList_GetItem(obj, 4), e);
  if (len > 5) PConvFromPyObject(G, PyList_GetItem(obj, 5), f);
  return true;
}

#endif
