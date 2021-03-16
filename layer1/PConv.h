
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
#ifndef _H_PConv
#define _H_PConv

#include"os_python.h"

#include"PyMOLGlobals.h"
#include"Base.h"
#include"OVLexicon.h"
#include "Result.h"

#include <array>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <vector>
#include <algorithm>

/* Convenient conversion routines for C<->Python data interchange
   
   Note that all of these routines assume that we have the global
   interpreter lock - blocking all other threads.
   
   There are three ways to get it:
   
   - call PBlock() [followe by PUnblock() when done]
   
   - call PBlockAndUnlockAPI - [followed by PLockAPIAndUnblock() when
   done]
   
   - or in response to a call to the PM API, you will have the main
   Python thread by default.  [Note that within an
   APIEntry(),APIExit() block the lock is released, so these
   functions should be called outside of that block].

*/

// CPythonVal macros
#define CPythonVal PyObject
#define CPythonVal_PyString_Check                       PyString_Check
#define CPythonVal_PyList_Check                         PyList_Check
#define CPythonVal_PyList_Size                          PyList_Size
#define CPythonVal_PyList_GetItem(G, list, i)           PyList_GetItem(list, i)
#define CPythonVal_PyDict_GetItemString(G, p, key)      PyDict_GetItemString(p, key)
#define CPythonVal_PConvPyIntToInt                                              PConvPyIntToInt
#define CPythonVal_PConvPyIntToInt_From_List(G, list, i, ptr)                   PConvPyIntToInt(PyList_GetItem(list, i), ptr)
#define CPythonVal_PConvPyFloatToFloat_From_List(G, list, i, ptr)               PConvPyFloatToFloat(PyList_GetItem(list, i), ptr)
#define CPythonVal_PConvPyListToIntArrayInPlace(G, obj, ff, ll)                 PConvPyListToIntArrayInPlace(obj, ff, ll)
#define CPythonVal_PConvPyListToIntArrayInPlace_From_List(G, list, i, ff, ll)   PConvPyListToIntArrayInPlace(PyList_GetItem(list, i), ff, ll)
#define CPythonVal_PConvPyListToFloatArrayInPlaceAutoZero_From_List(G, list, i, ...) \
                   PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list, i), __VA_ARGS__)
#define CPythonVal_PConvPyListToFloatArrayInPlace_From_List(G, list, i, ...) \
                   PConvPyListToFloatArrayInPlace(PyList_GetItem(list, i), __VA_ARGS__)
#define CPythonVal_PConvPyListToFloatVLANoneOkay_From_List(G, list, i, f)       PConvPyListToFloatVLANoneOkay(PyList_GetItem(list, i), f)
#define CPythonVal_PConvPyListToLabPosVec(G, obj)                               PConvPyListToLabPosVec(obj)
#define CPythonVal_PConvPyStrToStr_From_List(G, list, i, ptr, l)                PConvPyStrToStr(PyList_GetItem(list, i), ptr, l)

#define CPythonVal_Free(obj)
#define CPythonVal_FreeAll(PYOBJECT)

#define CPythonVal_New(G, PYOBJECT)                     PYOBJECT
#define CPythonVal_Append_List(LIST, ITEM)              PyList_Append(LIST, ITEM)
#define CPythonVal_New_List()                           PyList_New(0)
#define CPythonVal_New_Tuple(SIZE)                      PyTuple_New(SIZE)
#define CPythonVal_Tuple_SetItem(TUPLE, ITEM, VAL)      PyTuple_SetItem(TUPLE, ITEM, VAL)
#define CPythonVal_New_String(BUF, LEN)                 PyString_FromStringAndSize(BUF, LEN)
#define CPythonVal_New_Boolean(VAL)                     (VAL ? PyBool_FromLong(1) : PyBool_FromLong(0))
#define CPythonVal_New_Integer(VAL)                     PyInt_FromLong(VAL)
#define CPythonVal_New_Float(VAL)                       PyFloat_FromDouble(VAL)

#define CPythonVal_IsNone(PYOBJECT)                     (PYOBJECT == Py_None)

/* == error-checking routines: true = success, false = failure. */


/* NOTE: the string routines will write strings up to the specified
 * length, PLUS a NULL...so watch out for array overruns */

int PConvAttrToStrMaxLen(PyObject * obj, const char *attr, char *str, ov_size ll);

int PConvPyListToBitmask(PyObject * obj, int *bitmask, ov_size ll);
int PConvPyListToExtent(PyObject * obj, float *mn, float *mx);

int PConvAttrToFloatArrayInPlace(PyObject * obj, const char *attr, float *ff, ov_size ll);
int PConvAttrToIntArrayInPlace(PyObject * obj, const char *attr, int *ff, ov_size ll);
int PConvAttrToPtr(PyObject * obj, const char *name, void **cobj);

int PConvCObjectToPtr(PyObject * obj, void **ptr);
int PConvPyListToStrVLAList(PyObject * obj, char **vla, int *n_str);

int PConvPyListToStringVLA(PyObject * obj, char **vla_ptr);
#define PConvPyListToIntVLA(obj, f)     PConvPyListToIntArrayImpl(obj, f, true)

int PConvPyStrToStr(PyObject * obj, char *ptr, int l);
#ifndef _PYMOL_NOPY
int PConvPyStrToStrPtr(PyObject * obj, const char **ptr);
#endif
int PConvPyStrToLexRef(PyObject * obj, OVLexicon * lex, int *lex_ref);
int PConvPyFloatToFloat(PyObject * obj, float *ptr);
int PConvPyIntToChar(PyObject * obj, char *ptr);
int PConvPyIntToInt(PyObject * obj, int *ptr);
int PConvPyBoolToInt(PyObject * obj, int *ptr);
pymol::Result<std::vector<LabPosType>> PConvPyListToLabPosVec(PyObject* obj);


/* Jenarix conventions -- returns before args */

ov_status PConvPyTupleToIntVLA(int **result, PyObject * tuple);
ov_status PConvPyTupleToFloatVLA(float **result, PyObject * tuple);


/* === end === */


/* categories below... */

PyObject *PConvFloatVLAToPyList(const float *vla);
PyObject *PConvFloatVLAToPyTuple(float *vla);
PyObject *PConvIntVLAToPyList(const int *vla);
PyObject *PConvIntVLAToPyTuple(int *vla);
PyObject *PConvIntArrayToPyList(const int *f, int l, bool dump_binary=false);
PyObject *PConvSIntArrayToPyList(const short int *f, int l);
PyObject *PConvSCharArrayToPyList(const signed char *f, int l);
PyObject *PConvLabPosVecToPyList(const std::vector<LabPosType>& vec);

void PConvFloat3ToPyObjAttr(PyObject * obj, const char *attr, const float *v);
void PConvFloatToPyObjAttr(PyObject * obj, const char *attr, float f);
void PConvIntToPyObjAttr(PyObject * obj, const char *attr, int i);
void PConvInt2ToPyObjAttr(PyObject * obj, const char *attr, const int *v);
void PConvStringToPyObjAttr(PyObject * obj, const char *attr, const char *f);

int PConvPyObjectToFloat(PyObject * object, float *value);
int PConvPyObjectToInt(PyObject * object, int *value);
int PConvPyObjectToChar(PyObject * object, char *value);


/* NOTE: the string routines will write strings up to the specified
 * length, PLUS a NULL...so watch out for array overruns */

int PConvPyObjectToStrMaxLen(PyObject * object, char *value, int ln);
int PConvPyObjectToStrMaxClean(PyObject * object, char *value, int ln);

PyObject *PConvStringListToPyList(int l, const char * const *str);
PyObject *PConvStringVLAToPyList(const char *str);

void PConv44PyListTo44f(PyObject * src, float *dest);   /* note loss of precision */

#define PConvPyListToFloatVLA(obj, f)   PConvPyListToFloatArrayImpl(obj, f, true)
int PConvPyListToFloatVLANoneOkay(PyObject * obj, float **f);
int PConvPyList3ToFloatVLA(PyObject * obj, float **f);
#define PConvPyListToFloatArray(obj, f) PConvPyListToFloatArrayImpl(obj, f, false)
int PConvPyListToFloatArrayImpl(PyObject * obj, float **f, bool as_vla);
int PConvPyListToDoubleArray(PyObject * obj, double **f);
int PConvPyListToFloatArrayInPlace(PyObject * obj, float *ff, ov_size ll);
int PConvPyListOrTupleToFloatArrayInPlace(PyObject * obj, float *ff, ov_size ll);
int PConvPyListToFloatArrayInPlaceAutoZero(PyObject * obj, float *ii, ov_size ll);

int PConvPyListToDoubleArrayInPlace(PyObject * obj, double *ff, ov_size ll);

PyObject *PConvFloatArrayToPyList(const float *f, int l, bool dump_binary=false);
PyObject *PConvFloatArrayToPyListNullOkay(const float *f, int l);
PyObject *PConvDoubleArrayToPyList(const double *f, int l);

#define PConvPyListToIntArray(obj, f)   PConvPyListToIntArrayImpl(obj, f, false)
int PConvPyListToIntArrayImpl(PyObject * obj, int **f, bool as_vla);
int PConvPyListToIntArrayInPlace(PyObject * obj, int *ff, ov_size ll);
int PConvPyListToIntArrayInPlaceAutoZero(PyObject * obj, int *ii, ov_size ll);

int PConvPyListToSIntArrayInPlaceAutoZero(PyObject * obj, short int *ii, ov_size ll);
int PConvPyListToSCharArrayInPlaceAutoZero(PyObject * obj, signed char *ii, ov_size ll);

PyObject *PConv3DIntArrayTo3DPyList(int ***array, int *dim);

PyObject *PConvPickleLoads(PyObject * str);
PyObject *PConvPickleDumps(PyObject * obj);
PyObject *PConvAutoNone(PyObject * result);     /* automatically own Py_None */
PyObject *PConvIntToPyDictItem(PyObject * dict, const char *key, int i);

/* ============================================================ */
/*
 * PConvToPyObject: Convert any standart type (primitives
 * and c++ std library) to a python object.
 *
 * Return value: New reference.
 */
inline PyObject * PConvToPyObject(PyObject * v) {
  return v;
}

inline PyObject * PConvToPyObject(int v) {
  return PyInt_FromLong(v);
}

inline PyObject * PConvToPyObject(std::size_t v) {
#ifndef _PYMOL_NOPY
  return PyLong_FromSize_t(v);
#else
  return PyInt_FromLong(v);
#endif
}

inline PyObject * PConvToPyObject(float v) {
  return PyFloat_FromDouble(v);
}

inline PyObject * PConvToPyObject(double v) {
  return PyFloat_FromDouble(v);
}

inline PyObject * PConvToPyObject(const std::string &v) {
  return PyString_FromString(v.c_str());
}

inline PyObject * PConvToPyObject(const char * v) {
#ifndef _PYMOL_NOPY
  if (!v) {
    Py_RETURN_NONE;
  }
#endif
  return PyString_FromString(v);
}

inline PyObject * PConvToPyObject(const float * v, int n) {
  return PConvFloatArrayToPyList((float*)v, n);
}

template <class T>
PyObject * PConvToPyObject(const std::vector<T> &v) {
  int n = v.size();
  PyObject * o = PyList_New(n);

  for (int i = 0; i < n; ++i) {
    PyList_SetItem(o, i, PConvToPyObject(v[i]));
  }

  return o;
}

template <class T, std::size_t N>
PyObject * PConvToPyObject(const std::array<T, N> &arr) {
  PyObject * o = PyList_New(N);

  for (int i = 0; i < N; ++i) {
    PyList_SetItem(o, i, PConvToPyObject(arr[i]));
  }

  return o;
}

/**
 * Convert a set to a Python list
 */
template <class T>
PyObject * PConvToPyObject(const std::set<T> &v) {
  size_t i = 0, n = v.size();
  PyObject * o = PyList_New(n);

  for (auto it = v.begin(); it != v.end(); ++it) {
    PyList_SET_ITEM(o, i++, PConvToPyObject(*it));
  }

  return o;
}

/**
 * Convert a map to a flat Python list
 *
 *     {k1: v1, k2: v2, ...} -> [k1, v1, k2, v2, ...]
 */
template <class K, class V>
PyObject * PConvToPyObject(const std::map<K, V> &v) {
  size_t i = 0, n = v.size();
  PyObject * o = PyList_New(n * 2);

  for (auto it = v.begin(); it != v.end(); ++it) {
    PyList_SET_ITEM(o, i++, PConvToPyObject(it->first));
    PyList_SET_ITEM(o, i++, PConvToPyObject(it->second));
  }

  return o;
}

/**
 * Convert a pair to a Python tuple
 */
template <class T1, class T2>
PyObject* PConvToPyObject(const std::pair<T1, T2> &v) {
  PyObject* o = PyTuple_New(2);
  PyTuple_SET_ITEM(o, 0, PConvToPyObject(v.first));
  PyTuple_SET_ITEM(o, 1, PConvToPyObject(v.second));
  return o;
}

namespace pymol
{
struct Void;
}

inline PyObject* PConvToPyObject(const pymol::Void&)
{
#ifndef _PYMOL_NOPY
  Py_RETURN_NONE;
#else
  return nullptr;
#endif
}

/* ============================================================ */
/**
 * PConvFromPyObject: Templated conversion of a python object to a
 * standart type (primitives and c++ std library).
 */
template <typename Int,
    typename std::enable_if<std::is_integral<Int>::value ||
                            std::is_enum<Int>::value>::type* = nullptr>
inline bool PConvFromPyObject(PyMOLGlobals*, PyObject* obj, Int& out)
{
  auto const value = PyInt_AsLong(obj);
  out = static_cast<Int>(value);
  return value != -1 || !PyErr_Occurred();
}

template <typename Float,
    typename std::enable_if<std::is_floating_point<Float>::value>::type* =
        nullptr>
inline bool PConvFromPyObject(PyMOLGlobals*, PyObject* obj, Float& out)
{
  out = PyFloat_AsDouble(obj);
  return out != -1.0 || !PyErr_Occurred();
}

inline bool PConvFromPyObject(PyMOLGlobals *, PyObject * obj, std::string &out) {
#ifdef _PYMOL_NOPY
  // this would also work with real Python but without error handling (could
  // segfault if `obj` is not of type `str`)
  out = PyString_AsSomeString(obj);
#else
  const char* buffer = PyUnicode_AsUTF8(obj);
  if (!buffer) {
    return false;
  }
  out = buffer;
#endif
  return true;
}

inline bool PConvFromPyObject(PyMOLGlobals *, PyObject * obj, float * out) {
  return PConvPyListToFloatArrayInPlace(obj, out, 0);
}

template <class T>
bool PConvFromPyObject(PyMOLGlobals * G, PyObject * obj, std::vector<T> &out) {
  if (PyBytes_Check(obj)) {
    // binary_dump
    size_t slen = PyBytes_Size(obj);

    if (slen % sizeof(T)) {
      return false;
    }

    out.resize(slen / sizeof(T));

    auto strval = PyBytes_AsSomeString(obj);
    std::copy_n(strval.data(), slen, reinterpret_cast<char*>(out.data()));
    return true;
  }

  if (!PyList_Check(obj))
    return false;

  int n = PyList_Size(obj);

  out.clear();
  out.reserve(n);

  for (int i = 0; i < n; ++i) {
    PyObject *item = PyList_GET_ITEM(obj, i);

    T t;
    if (!PConvFromPyObject(G, item, t))
      return false;

    out.push_back(t);
  }

  return true;
}

/**
 * Convert a Python list to a set
 */
template <class T>
bool PConvFromPyObject(PyMOLGlobals * G, PyObject * obj, std::set<T> &out) {
  if (!PyList_Check(obj))
    return false;

  int n = PyList_Size(obj);

  out.clear();

  for (int i = 0; i < n; ++i) {
    PyObject *item = PyList_GET_ITEM(obj, i);

    T t;
    if (!PConvFromPyObject(G, item, t))
      return false;

    out.insert(t);
  }

  return true;
}

/**
 * Convert a flat Python list to a map, even indices are keys and odd
 * indices are values.
 *
 *     [a, b, c, d, ...] -> {a: b, c: d, ...}
 */
template <class K, class V>
bool PConvFromPyObject(PyMOLGlobals * G, PyObject * obj, std::map<K, V> &out) {
  if (!PyList_Check(obj))
    return false;

  int n = PyList_Size(obj);

  out.clear();

  for (int i = 0; i < n - 1;) {
    PyObject *key   = PyList_GET_ITEM(obj, i++);
    PyObject *value = PyList_GET_ITEM(obj, i++);

    K k;
    if (!PConvFromPyObject(G, key, k))
      return false;

    if (!PConvFromPyObject(G, value, out[k]))
      return false;
  }

  return true;
}

/* ============================================================ */

template <class T>
bool PConvFromPyListItem(PyMOLGlobals* G, PyObject* list, size_t i, T& out)
{
  auto item = PyList_GetItem(list, i);
  auto ok = PConvFromPyObject(G, item, out);
  CPythonVal_Free(item);
  return ok;
}

#endif
