
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
#include"os_python.h"

#include"os_predef.h"
#include"os_python.h"

#include"os_std.h"

#include"Err.h"
#include"MemoryDebug.h"
#include"Base.h"
#include"PConv.h"
#include"P.h"
#include"Util.h"

#include <vector>

#define pickle_mod_name "pickle"

/* Return value: New reference.
 * Load a pickle from the given string
 */
PyObject *PConvPickleLoads(PyObject * str)
{
  PyObject *picklemod = NULL, *obj = NULL;
  ok_assert(1, picklemod = PyImport_ImportModule(pickle_mod_name));
  obj = PYOBJECT_CALLMETHOD(picklemod, "loads", "O", str);
ok_except1:
  Py_XDECREF(picklemod);
  return obj;
}

/* Return value: New reference.
 * Return a string containing an object in pickle format.
 */
PyObject *PConvPickleDumps(PyObject * obj)
{
  PyObject *picklemod = NULL, *str = NULL;
  ok_assert(1, picklemod = PyImport_ImportModule(pickle_mod_name));
  str = PYOBJECT_CALLMETHOD(picklemod, "dumps", "Oi", obj, 1);
ok_except1:
  Py_XDECREF(picklemod);
  return str;
}

PyObject *PConvAutoNone(PyObject * result)
{                               /* automatically own Py_None */
  if(result == Py_None)
    Py_INCREF(result);
  else if(result == NULL) {
    result = Py_None;
    Py_INCREF(result);
  }
  return (result);
}


/* Error-checked utility routines */

int PConvPyListToBitmask(PyObject * obj, int *bitmask, ov_size ll)
{
  std::vector<signed char> visRepArr(ll, 0);

  if (ll > 0)
    ok_assert(1, PConvPyListToSCharArrayInPlaceAutoZero(obj, &visRepArr[0], ll));

  *bitmask = 0;
  for(int i = 0; i < ll; i++)
    if(visRepArr[i])
      SET_BIT(*bitmask, i);

  return true;
ok_except1:
  return false;
}

int PConvPyListToExtent(PyObject * obj, float *mn, float *mx)
{                               /* [[min_x,min_y,min_z],
                                   [max_x,max_y,max_z]] */
  int ok = false;
  PyObject *t1, *t2;
  if(!obj) {
    ok = false;
  } else if(PyList_Check(obj))
    if(PyList_Size(obj) == 2) {
      t1 = PyList_GetItem(obj, 0);
      t2 = PyList_GetItem(obj, 1);
      if(PConvPyListToFloatArrayInPlace(t1, mn, 3) &&
         PConvPyListToFloatArrayInPlace(t2, mx, 3))
        ok = true;
    }
  return (ok);
}

int PConvPyListToStrVLAList(PyObject * obj, char **vla, int *n_str)
{
  int ok = false;
  PyObject *t;
  int n_st = 0, n_ch = 0, nn_ch, l, i;
  if(!*vla)
    *vla = VLAlloc(char, 10);
  if((!obj) || (!*vla)) {
    ok = false;
  } else if(PyList_Check(obj)) {
    n_st = PyList_Size(obj);
    ok = true;
    for(i = 0; i < n_st; i++) {
      t = PyList_GetItem(obj, i);
      if(PyString_Check(t)) {
        l = PyString_Size(t);
        nn_ch = n_ch + l + 1;
        VLACheck(*vla, char, nn_ch);
        auto strval = PyString_AsSomeString(t);
        UtilNCopy((*vla) + n_ch, strval.c_str(), l + 1);
        n_ch = nn_ch;
      } else {
        VLACheck(*vla, char, n_ch + 1);
        (*vla)[n_ch] = 0;
        n_ch++;
      }
    }
  }
  *n_str = n_st;
  return (ok);
}

int PConvAttrToIntArrayInPlace(PyObject * obj, const char *attr, int *f, ov_size ll)
{
  int ok = true;
  PyObject *tmp;
  if(!obj) {
    ok = false;
  } else if(PyObject_HasAttrString(obj, attr)) {
    tmp = PyObject_GetAttrString(obj, attr);
    ok = PConvPyListToIntArrayInPlace(tmp, f, ll);
    Py_DECREF(tmp);
  } else {
    ok = false;
  }
  return (ok);
}

int PConvAttrToFloatArrayInPlace(PyObject * obj, const char *attr, float *f, ov_size ll)
{
  int ok = true;
  PyObject *tmp;
  if(!obj) {
    ok = false;
  } else if(PyObject_HasAttrString(obj, attr)) {
    tmp = PyObject_GetAttrString(obj, attr);
    ok = PConvPyListToFloatArrayInPlace(tmp, f, ll);
    Py_DECREF(tmp);
  } else {
    ok = false;
  }
  return (ok);
}

int PConvAttrToStrMaxLen(PyObject * obj, const char *attr, char *str, ov_size ll)
{
  int ok = true;
  PyObject *tmp;
  if(!obj) {
    ok = false;
  } else if(PyObject_HasAttrString(obj, attr)) {
    tmp = PyObject_GetAttrString(obj, attr);
    ok = PConvPyObjectToStrMaxLen(tmp, str, ll);
    Py_DECREF(tmp);
  } else {
    ok = false;
  }
  return (ok);
}

int PConvAttrToPtr(PyObject * obj, const char *attr, void **cobj)
{
  PyObject *tmp;
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(PyObject_HasAttrString(obj, attr)) {
    tmp = PyObject_GetAttrString(obj, attr);
    ok = PConvCObjectToPtr(tmp, cobj);
    Py_DECREF(tmp);
  } else {
    ok = false;
  }
  return (ok);
}

/**
 * Used via `chempy.map` with
 * https://github.com/cctbx/cctbx_project/blob/master/cctbx/maptbx/boost_python/pymol_interface.cpp
 */
int PConvCObjectToPtr(PyObject * obj, void **ptr)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(!PyCapsule_CheckExact(obj))
    ok = false;
  else
    (*ptr) = PyCapsule_GetPointer(obj, nullptr);
  return (ok);
}

int PConvPyStrToLexRef(PyObject * obj, OVLexicon * lex, int *lex_ref)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(!PyString_Check(obj)) {
    ok = false;
  } else {
    auto strval = PyString_AsSomeString(obj);
    if(!strval.c_str()) {
      ok = false;
    } else {
      OVreturn_word result = OVLexicon_GetFromCString(lex, strval.c_str());
      if(OVreturn_IS_OK(result)) {
        *lex_ref = result.word;
      } else {
        ok = false;
      }
    }
  }
  return ok;
}

#ifndef _PYMOL_NOPY
int PConvPyStrToStrPtr(PyObject * obj, const char **ptr)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(!PyString_Check(obj)) {
    ok = false;
  }
  if(ok)
    *ptr = PyString_AsString(obj);
  return (ok);
}
#endif

int PConvPyStrToStr(PyObject * obj, char *ptr, int size)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(PyBytes_Check(obj)) {
    auto strval = PyBytes_AsSomeString(obj);
    UtilNCopy(ptr, strval.c_str(), size);
  } else if(!PyString_Check(obj)) {
    ok = false;
    if(size)
      *ptr = 0;
  } else {
    auto strval = PyString_AsSomeString(obj);
    UtilNCopy(ptr, strval.c_str(), size);
  }
  return (ok);
}

int PConvPyIntToInt(PyObject * obj, int *ptr)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(PyLong_Check(obj)) {
    *ptr = (int) PyLong_AsLongLong(obj);
  } else if(PyInt_Check(obj)) {
    *ptr = PyInt_AsLong(obj);
  } else {
    ok = false;
  }
  return (ok);
}

int PConvPyFloatToFloat(PyObject * obj, float *ptr)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(!PyFloat_Check(obj)) {
    ok = false;
  } else {
    *ptr = (float) PyFloat_AsDouble(obj);
  }
  return (ok);
}

int PConvPyIntToChar(PyObject * obj, char *ptr)
{
  int ok = true;
  if(!obj) {
    ok = false;
  } else if(!PyInt_Check(obj)) {
    if(!PyLong_Check(obj)) {
      ok = false;
    } else {
      *ptr = (char) PyLong_AsLongLong(obj);
    }
  } else {
    *ptr = (char) PyInt_AsLong(obj);
  }
  return (ok);
}


/* == end == */

int PConvPyObjectToFloat(PyObject * object, float *value)
{
  int result = true;
  PyObject *tmp;
  if(!object)
    result = false;
  else if(PyFloat_Check(object)) {
    (*value) = (float) PyFloat_AsDouble(object);
  } else if(PyLong_Check(object)) {
    (*value) = (float) PyLong_AsLongLong(object);
  } else if(PyInt_Check(object)) {
    (*value) = (float) PyInt_AsLong(object);
  } else {
    tmp = PyNumber_Float(object);
    if(tmp) {
      (*value) = (float) PyFloat_AsDouble(tmp);
      Py_DECREF(tmp);
    } else
      result = false;
  }
  return (result);
}

int PConvPyObjectToInt(PyObject * object, int *value)
{
  int result = true;
  PyObject *tmp;
  if(!object) {
    result = false;
  } else if(PyLong_Check(object)) {
    (*value) = (int) PyLong_AsLongLong(object);
  } else if(PyInt_Check(object)) {
    (*value) = (int) PyInt_AsLong(object);
  } else {
    tmp = PyNumber_Int(object);
    if(tmp) {
      (*value) = (int) PyInt_AsLong(tmp);
      Py_DECREF(tmp);
    } else
      result = false;
  }
  return (result);
}

int PConvPyObjectToChar(PyObject * object, char *value)
{
  int result = true;
  PyObject *tmp;
  if(!object)
    result = false;
  else if(PyInt_Check(object)) {
    (*value) = (char) PyInt_AsLong(object);
  } else if(PyLong_Check(object)) {
    (*value) = (char) PyLong_AsLongLong(object);
  } else {
    tmp = PyNumber_Int(object);
    if(tmp) {
      (*value) = (char) PyInt_AsLong(tmp);
      Py_DECREF(tmp);
    } else
      result = false;
  }
  return (result);
}

int PConvPyObjectToStrMaxLen(PyObject * object, char *value, int ln)
{
  PyObject *tmp;
  int result = true;
  if(!object) {
    result = false;
  } else if(PyBytes_Check(object)) {
    auto strval = PyBytes_AsSomeString(object);
    strncpy(value, strval.c_str(), ln);
  } else if(PyString_Check(object)) {
    auto strval = PyString_AsSomeString(object);
    strncpy(value, strval.c_str(), ln);
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      auto strval = PyString_AsSomeString(tmp);
      strncpy(value, strval.c_str(), ln);
      Py_DECREF(tmp);
    } else
      result = 0;
  }
  if(ln > 0)
    value[ln] = 0;
  else
    value[0] = 0;
  return (result);
}

int PConvPyObjectToStrMaxClean(PyObject * object, char *value, int ln)
{
  PyObject *tmp;
  int result = true;
  if(!object)
    result = false;
  else if(PyString_Check(object)) {
    auto strval = PyString_AsSomeString(object);
    strncpy(value, strval.c_str(), ln);
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      auto strval = PyString_AsSomeString(tmp);
      strncpy(value, strval.c_str(), ln);
      Py_DECREF(tmp);
    } else
      result = 0;
  }
  if(ln > 0)
    value[ln] = 0;
  else
    value[0] = 0;
  UtilCleanStr(value);
  return (result);
}

PyObject *PConvFloatToPyDictItem(PyObject * dict, const char *key, float f)
{
  PyObject *tmp;
  tmp = PyFloat_FromDouble((double) f);
  PyDict_SetItemString(dict, key, tmp);
  Py_XDECREF(tmp);
  return (tmp);
}

PyObject *PConvIntToPyDictItem(PyObject * dict, const char *key, int i)
{
  PyObject *tmp;
  tmp = PyInt_FromLong(i);
  PyDict_SetItemString(dict, key, tmp);
  Py_XDECREF(tmp);
  return (tmp);
}

void PConvFloat3ToPyObjAttr(PyObject * obj, const char *attr, const float *v)
{
  PyObject *t1, *t2, *t3, *tmp;

  t1 = PyFloat_FromDouble((double) v[0]);
  t2 = PyFloat_FromDouble((double) v[1]);
  t3 = PyFloat_FromDouble((double) v[2]);
  tmp = PyList_New(3);
  if(t1 && t2 && t3 && tmp) {
    PyList_SetItem(tmp, 0, t1); /* steals reference */
    PyList_SetItem(tmp, 1, t2); /* steals reference */
    PyList_SetItem(tmp, 2, t3); /* steals reference */
    PyObject_SetAttrString(obj, attr, tmp);
  }
  Py_XDECREF(tmp);
}

void PConvInt2ToPyObjAttr(PyObject * obj, const char *attr, const int *v)
{
  PyObject *t1, *t2, *tmp;

  t1 = PyInt_FromLong((long) v[0]);
  t2 = PyInt_FromLong((long) v[1]);
  tmp = PyList_New(2);
  if(t1 && t2 && tmp) {
    PyList_SetItem(tmp, 0, t1); /* steals reference */
    PyList_SetItem(tmp, 1, t2); /* steals reference */
    PyObject_SetAttrString(obj, attr, tmp);
  }
  Py_XDECREF(tmp);
}

void PConvFloatToPyObjAttr(PyObject * obj, const char *attr, float f)
{
  PyObject *tmp;
  tmp = PyFloat_FromDouble((double) f);
  PyObject_SetAttrString(obj, attr, tmp);
  Py_DECREF(tmp);
}

void PConvIntToPyObjAttr(PyObject * obj, const char *attr, int i)
{
  PyObject *tmp;
  tmp = PyInt_FromLong(i);
  PyObject_SetAttrString(obj, attr, tmp);
  Py_DECREF(tmp);
}

void PConvStringToPyObjAttr(PyObject * obj, const char *attr, const char *f)
{
  PyObject *tmp;
  tmp = PyString_FromString(f);
  PyObject_SetAttrString(obj, attr, tmp);
  Py_DECREF(tmp);
}

int PConvPyListToFloatArrayImpl(PyObject * obj, float **f, bool as_vla)
{
  int a, l;
  int ok = true;
  float *ff;
  if(!obj) {
    *f = NULL;
    ok = false;
  } else if (PyBytes_Check(obj)){
    // binary_dump
    int slen = PyBytes_Size(obj);
    l = slen / sizeof(float);

    if (as_vla) {
      (*f) = VLAlloc(float, l);
    } else {
      (*f) = pymol::malloc<float>(l);
    }

    auto strval = PyBytes_AsSomeString(obj);
    memcpy(*f, strval.data(), slen);
  } else if(!PyList_Check(obj)) {
    *f = NULL;
    ok = false;
  } else {
    l = (int) PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;

    if (as_vla) {
      (*f) = VLAlloc(float, l);
    } else {
      (*f) = pymol::malloc<float>(l);
    }

    ff = (*f);
    for(a = 0; a < l; a++)
      *(ff++) = (float) PyFloat_AsDouble(PyList_GetItem(obj, a));

  }
  return (ok);
}

int PConvPyListToFloatVLANoneOkay(PyObject * obj, float **f)
{
  int a, l;
  float *ff;
  int ok = true;
  if(!obj) {
    *f = NULL;
    ok = false;
  } else if(obj == Py_None) {
    *f = NULL;
    ok = true;
  } else if(!PyList_Check(obj)) {
    *f = NULL;
    ok = false;
  } else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    (*f) = VLAlloc(float, l);
    ff = (*f);
    for(a = 0; a < l; a++)
      *(ff++) = (float) PyFloat_AsDouble(PyList_GetItem(obj, a));
    VLASize((*f), float, l);
  }
  return (ok);

}

int PConvPyList3ToFloatVLA(PyObject * obj, float **f)
{
  int a, b, l;
  float *ff;
  PyObject *triple;
  int ok = true;
  if(!obj) {
    *f = NULL;
    ok = false;
  } else if(!PyList_Check(obj)) {
    *f = NULL;
    ok = false;
  } else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    (*f) = VLAlloc(float, l * 3);
    ff = (*f);
    for(a = 0; a < l; a++) {
      triple = PyList_GetItem(obj, a);
      ok = PyList_Check(triple);
      if(ok)
        ok = (PyList_Size(triple) == 3);
      if(ok) {
        for(b = 0; b < 3; b++)
          *(ff++) = (float) PyFloat_AsDouble(PyList_GetItem(triple, b));
      } else {
        ok = false;
        break;
      }
    }
    VLASize((*f), float, l * 3);
  }
  return (ok);
}

int PConvPyListToDoubleArray(PyObject * obj, double **f)
{
  int a, l;
  double *ff;
  int ok = true;
  if(!obj) {
    *f = NULL;
    l = 0;
  } else if(!PyList_Check(obj)) {
    *f = NULL;
    ok = false;
  } else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    (*f) = pymol::malloc<double>(l);
    ff = (*f);
    for(a = 0; a < l; a++)
      *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj, a));
  }
  return (ok);
}

int PConvPyListToIntArrayImpl(PyObject * obj, int **f, bool as_vla)
{
  int a, l;
  int *ff;
  int ok = true;
  if(!obj) {
    *f = NULL;
    ok = false;
  } else if (PyBytes_Check(obj)){
    // binary_dump
    int slen = PyBytes_Size(obj);
    l = slen / sizeof(int);

    if (as_vla) {
      (*f) = VLAlloc(int, l);
    } else {
      (*f) = pymol::malloc<int>(l);
    }

    auto strval = PyBytes_AsSomeString(obj);
    memcpy(*f, strval.data(), slen);
  } else if(!PyList_Check(obj)) {
    *f = NULL;
    ok = false;
  } else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;

    if (as_vla) {
      (*f) = VLAlloc(int, l);
    } else {
      (*f) = pymol::malloc<int>(l);
    }

    ff = (*f);
    for(a = 0; a < l; a++)
      *(ff++) = PyInt_AsLong(PyList_GetItem(obj, a));
  }
  return (ok);
}

ov_status PConvPyTupleToIntVLA(int **result, PyObject * tuple)
{
  ov_status status = OV_STATUS_FAILURE;
  if(!(tuple && PyTuple_Check(tuple))) {
    *result = NULL;
  } else {
    int *vla = NULL;
    ov_size size = PyTuple_Size(tuple);
    vla = VLAlloc(int, size);
    if(vla) {
      ov_size i;
      int *ptr = vla;
      status = OV_STATUS_SUCCESS;
      for(i = 0; i < size; i++)
        *(ptr++) = PyInt_AsLong(PyTuple_GetItem(tuple, i));
    }
    *result = vla;
  }
  return status;
}

ov_status PConvPyTupleToFloatVLA(float **result, PyObject * tuple)
{
  ov_status status = OV_STATUS_FAILURE;
  if(!(tuple && PyTuple_Check(tuple))) {
    *result = NULL;
  } else {
    float *vla = NULL;
    ov_size size = PyTuple_Size(tuple);
    vla = VLAlloc(float, size);
    if(vla) {
      ov_size i;
      float *ptr = vla;
      status = OV_STATUS_SUCCESS;
      for(i = 0; i < size; i++)
        *(ptr++) = (float) PyFloat_AsDouble(PyTuple_GetItem(tuple, i));
    }
    *result = vla;
  }
  return status;
}

int PConvPyListToDoubleArrayInPlace(PyObject * obj, double *ff, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj) {
    ok = false;
  } else if(!PyList_Check(obj)) {
    ok = false;
  } else {
    l = PyList_Size(obj);
    if(l != ll)
      ok = false;
    else {
      if(!l)
        ok = -1;
      else
        ok = l;
      for(a = 0; a < l; a++)
        *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj, a));
    }
    /* NOTE ASSUMPTION! */
  }
  return (ok);
}

int PConvPyListToFloatArrayInPlace(PyObject * obj, float *ff, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj) {
    ok = false;
  } else if(!PyList_Check(obj)) {
    ok = false;
  } else {
    l = PyList_Size(obj);
    if(ll > 0 && l != ll)
      ok = false;
    else {
      if(!l)
        ok = -1;
      else
        ok = l;
      for(a = 0; a < l; a++)
        *(ff++) = (float) PyFloat_AsDouble(PyList_GetItem(obj, a));
    }
    /* NOTE ASSUMPTION! */
  }
  return (ok);
}

int PConvPyListOrTupleToFloatArrayInPlace(PyObject * obj, float *ff, ov_size ll)
{
  int ok = true, isTuple = false;
  ov_size a, l;
  if(!obj) {
    ok = false;
  } else if(!PyList_Check(obj) && !(isTuple=PyTuple_Check(obj))) {
    ok = false;
  } else {
    if (isTuple)
      l = PyTuple_Size(obj);
    else
      l = PyList_Size(obj);
    if(l != ll)
      ok = false;
    else {
      if(!l)
        ok = -1;
      else
        ok = l;
      if (isTuple)
	for(a = 0; a < l; a++){
	  *(ff++) = (float) PyFloat_AsDouble(PyTuple_GetItem(obj, a));
	}
      else
	for(a = 0; a < l; a++){
	  *(ff++) = (float) PyFloat_AsDouble(PyList_GetItem(obj, a));
	}
    }
    /* NOTE ASSUMPTION! */
  }
  return (ok);
}

int PConvPyListToIntArrayInPlace(PyObject * obj, int *ii, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj)
    ok = false;
  else if(!PyList_Check(obj))
    ok = false;
  else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    if(l != ll)
      ok = false;
    else
      for(a = 0; a < l; a++)
        *(ii++) = PyInt_AsLong(PyList_GetItem(obj, a));
    /* NOTE ASSUMPTION! */
  }
  return (ok);
}

int PConvPyListToIntArrayInPlaceAutoZero(PyObject * obj, int *ii, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj)
    ok = false;
  else if(!PyList_Check(obj))
    ok = false;
  else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    for(a = 0; (a < l) && (a < ll); a++)
      *(ii++) = PyInt_AsLong(PyList_GetItem(obj, a));
    while(l < ll) {
      *(ii++) = 0;
      l++;
    }
  }
  return (ok);
}

int PConvPyListToSIntArrayInPlaceAutoZero(PyObject * obj, short int *ii, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj)
    ok = false;
  else if(!PyList_Check(obj))
    ok = false;
  else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    for(a = 0; (a < l) && (a < ll); a++)
      *(ii++) = (short int) PyInt_AsLong(PyList_GetItem(obj, a));
    while(l < ll) {
      *(ii++) = 0;
      l++;
    }
  }
  return (ok);
}

int PConvPyListToSCharArrayInPlaceAutoZero(PyObject * obj, signed char *ii, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj)
    ok = false;
  else if(!PyList_Check(obj))
    ok = false;
  else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    for(a = 0; (a < l) && (a < ll); a++)
      *(ii++) = (signed char) PyInt_AsLong(PyList_GetItem(obj, a));
    while(l < ll) {
      *(ii++) = 0;
      l++;
    }
  }
  return (ok);
}

int PConvPyListToFloatArrayInPlaceAutoZero(PyObject * obj, float *ii, ov_size ll)
{
  int ok = true;
  ov_size a, l;
  if(!obj)
    ok = false;
  else if(!PyList_Check(obj))
    ok = false;
  else {
    l = PyList_Size(obj);
    if(!l)
      ok = -1;
    else
      ok = l;
    for(a = 0; (a < l) && (a < ll); a++)
      *(ii++) = (float) PyFloat_AsDouble(PyList_GetItem(obj, a));
    while(l < ll) {
      *(ii++) = 0.0f;
      l++;
    }
  }
  return (ok);
}

PyObject *PConvFloatArrayToPyList(const float *f, int l, bool dump_binary)
{
#ifndef PICKLETOOLS
  if (dump_binary){
    return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(f), l * sizeof(float));
  } 
#endif
  int a;
  PyObject *result = PyList_New(l);
  for(a = 0; a < l; a++)
    PyList_SetItem(result, a, PyFloat_FromDouble((double) *(f++)));
  return (PConvAutoNone(result));
}

PyObject *PConvFloatArrayToPyListNullOkay(const float *f, int l)
{
  int a;
  PyObject *result = NULL;
  if(f) {
    result = PyList_New(l);
    for(a = 0; a < l; a++)
      PyList_SetItem(result, a, PyFloat_FromDouble((double) *(f++)));
  }
  return (PConvAutoNone(result));
}

PyObject *PConvDoubleArrayToPyList(const double *f, int l)
{
  int a;
  PyObject *result = PyList_New(l);
  for(a = 0; a < l; a++)
    PyList_SetItem(result, a, PyFloat_FromDouble(*(f++)));
  return (PConvAutoNone(result));
}

PyObject *PConvFloatVLAToPyList(const float *f)
{
  int a, l;
  PyObject *result = NULL;
  l = VLAGetSize(f);
  result = PyList_New(l);
  for(a = 0; a < l; a++) {
    PyList_SetItem(result, a, PyFloat_FromDouble((double) *(f++)));     /* set item steals ref */
  }
  return (PConvAutoNone(result));
}

PyObject *PConvFloatVLAToPyTuple(float *vla)
{
  PyObject *result = NULL;
  if(vla) {
    ov_size size = VLAGetSize(vla);
    result = PyTuple_New(size);
    if(result) {
      ov_size i;
      for(i = 0; i < size; i++) {
        PyTuple_SetItem(result, i, PyFloat_FromDouble((double) *(vla++)));      /* set item steals ref */
      }
    }
  }
  return (PConvAutoNone(result));
}

PyObject *PConvIntVLAToPyList(const int *f)
{
  int a, l;
  PyObject *result = NULL;
  l = VLAGetSize(f);
  result = PyList_New(l);
  for(a = 0; a < l; a++)
    PyList_SetItem(result, a, PyInt_FromLong(*(f++)));
  return (PConvAutoNone(result));
}

PyObject *PConvIntVLAToPyTuple(int *vla)
{
  PyObject *result = NULL;
  if(vla) {
    ov_size size = VLAGetSize(vla);
    result = PyTuple_New(size);
    if(result) {
      ov_size i;
      for(i = 0; i < size; i++) {
        PyTuple_SetItem(result, i, PyInt_FromLong(*(vla++)));   /* set item steals ref */
      }
    }
  }
  return (PConvAutoNone(result));
}

PyObject *PConvIntArrayToPyList(const int *f, int l, bool dump_binary)
{
#ifndef PICKLETOOLS
  if (dump_binary){
    return PyBytes_FromStringAndSize(reinterpret_cast<const char*>(f), l * sizeof(int));
  }
#endif
  int a;
  PyObject *result = PyList_New(l);
  for(a = 0; a < l; a++)
    PyList_SetItem(result, a, PyInt_FromLong(*(f++)));
  return (PConvAutoNone(result));
}

PyObject *PConvSIntArrayToPyList(const short int *f, int l)
{
  int a;
  PyObject *result = PyList_New(l);
  for(a = 0; a < l; a++)
    PyList_SetItem(result, a, PyInt_FromLong(*(f++)));
  return (PConvAutoNone(result));
}

PyObject *PConvSCharArrayToPyList(const signed char *f, int l)
{
  int a;
  PyObject *result = PyList_New(l);
  for(a = 0; a < l; a++)
    PyList_SetItem(result, a, PyInt_FromLong(*(f++)));
  return (PConvAutoNone(result));
}

PyObject *PConv3DIntArrayTo3DPyList(int ***array, int *dim)
{
  int a, b, c;
  PyObject *result, *pyB, *pyC;
  result = PyList_New(dim[0]);
  for(a = 0; a < dim[0]; a++) {
    pyB = PyList_New(dim[1]);
    PyList_SetItem(result, a, pyB);
    for(b = 0; b < dim[1]; b++) {
      pyC = PyList_New(dim[2]);
      PyList_SetItem(pyB, b, pyC);
      for(c = 0; c < dim[2]; c++) {
        PyList_SetItem(pyC, c, PyInt_FromLong(array[a][b][c]));
      }
    }
  }
  return (PConvAutoNone(result));
}

PyObject *PConvStringListToPyList(int l, const char * const *str)
{
  int a;
  PyObject *result = PyList_New(l);
  for(a = 0; a < l; a++) {
    PyList_SetItem(result, a, PyString_FromString(str[a]));
  }
  return (PConvAutoNone(result));
}

/**
 * Converts concatenated null-terminated strings to a list of strings.
 */
PyObject *PConvStringVLAToPyList(const char *vla)
{
  int a, c, n = 0;
  const char *p;
  PyObject *result = NULL;
  p = vla;
  c = VLAGetSize(vla);
  while(c--) {                  /* count strings */
    if(!*(p++))
      n++;
  }

  result = PyList_New(n);
  p = vla;
  for(a = 0; a < n; a++) {
    PyList_SetItem(result, a, PyString_FromString(p));
    while(*(p++));
  }
  return (PConvAutoNone(result));
}

/**
 * Concatenates a list of null-terminated strings (the null-terminator serves as
 * the delimiter).
 */
int PConvPyListToStringVLA(PyObject * obj, char **vla_ptr)
{
  int a, l, ll;
  char *vla = NULL, *q;
  PyObject *i;
  if(obj)
    if(PyList_Check(obj)) {
      l = PyList_Size(obj);
      ll = 0;
      for(a = 0; a < l; a++) {
        i = PyList_GetItem(obj, a);
        if(PyString_Check(i)) {
          ll += PyString_Size(i) + 1;
        }
      }
      vla = VLAlloc(char, ll);
      VLASize(vla, char, ll);
      q = vla;
      for(a = 0; a < l; a++) {
        i = PyList_GetItem(obj, a);
        if(PyString_Check(i)) {
          auto strval = PyString_AsSomeString(i);
          auto p = strval.c_str();
          while(*p)
            *(q++) = *(p++);
          *(q++) = 0;
        }
      }
    }
  (*vla_ptr) = vla;
  return (vla && 1);
}

pymol::Result<std::vector<LabPosType>> PConvPyListToLabPosVec(PyObject* obj)
{
  std::vector<LabPosType> result;
  if(obj)
    if(PyList_Check(obj)) {
      auto l = PyList_Size(obj);
      result.resize(l);
      for(int a = 0; a < l; a++) {
        auto i = PyList_GetItem(obj, a);
        auto q = &result[a];
        if(PyList_Check(i) && (PyList_Size(i) == 7)) {
          auto ok = PConvPyIntToInt(PyList_GetItem(i, 0), &q->mode) &&
            PConvPyFloatToFloat(PyList_GetItem(i, 1), q->pos) &&
            PConvPyFloatToFloat(PyList_GetItem(i, 2), q->pos + 1) &&
            PConvPyFloatToFloat(PyList_GetItem(i, 3), q->pos + 2) &&
            PConvPyFloatToFloat(PyList_GetItem(i, 4), q->offset) &&
            PConvPyFloatToFloat(PyList_GetItem(i, 5), q->offset + 1) &&
            PConvPyFloatToFloat(PyList_GetItem(i, 6), q->offset + 2);
          if (!ok) {
            return pymol::make_error("Invalid subitem.");
          }
        } else {
          return pymol::make_error("Invalid sublist.");
        }
      }
    }
  return result;
}

PyObject* PConvLabPosVecToPyList(const std::vector<LabPosType>& vec)
{                               /* TO DO error handling */
  PyObject* result = nullptr;
  if(!vec.empty()) {
    result = PyList_New(vec.size());
    for (int a = 0; a < vec.size(); a++) {
      const auto& p = vec[a];
      auto item = PyList_New(7);
      PyList_SetItem(item, 0, PyInt_FromLong(p.mode));
      PyList_SetItem(item, 1, PyFloat_FromDouble((double) p.pos[0]));
      PyList_SetItem(item, 2, PyFloat_FromDouble((double) p.pos[1]));
      PyList_SetItem(item, 3, PyFloat_FromDouble((double) p.pos[2]));
      PyList_SetItem(item, 4, PyFloat_FromDouble((double) p.offset[0]));
      PyList_SetItem(item, 5, PyFloat_FromDouble((double) p.offset[1]));
      PyList_SetItem(item, 6, PyFloat_FromDouble((double) p.offset[2]));
      PyList_SetItem(result, a, item);
    }
  }
  return (PConvAutoNone(result));
}

void PConv44PyListTo44f(PyObject * src, float *dest)
{                               /* note lost of precision */
  PyObject *row;
  if(src && dest && PyList_Check(src)) {
    row = PyList_GetItem(src, 0);
    if(row && PyList_Check(row)) {
      dest[0] = (float) PyFloat_AsDouble(PyList_GetItem(row, 0));
      dest[1] = (float) PyFloat_AsDouble(PyList_GetItem(row, 1));
      dest[2] = (float) PyFloat_AsDouble(PyList_GetItem(row, 2));
      dest[3] = (float) PyFloat_AsDouble(PyList_GetItem(row, 3));
    }
    row = PyList_GetItem(src, 1);
    if(row && PyList_Check(row)) {
      dest[4] = (float) PyFloat_AsDouble(PyList_GetItem(row, 0));
      dest[5] = (float) PyFloat_AsDouble(PyList_GetItem(row, 1));
      dest[6] = (float) PyFloat_AsDouble(PyList_GetItem(row, 2));
      dest[7] = (float) PyFloat_AsDouble(PyList_GetItem(row, 3));
    }
    row = PyList_GetItem(src, 2);
    if(row && PyList_Check(row)) {
      dest[8] = (float) PyFloat_AsDouble(PyList_GetItem(row, 0));
      dest[9] = (float) PyFloat_AsDouble(PyList_GetItem(row, 1));
      dest[10] = (float) PyFloat_AsDouble(PyList_GetItem(row, 2));
      dest[11] = (float) PyFloat_AsDouble(PyList_GetItem(row, 3));
    }
    row = PyList_GetItem(src, 3);
    if(row && PyList_Check(row)) {
      dest[12] = (float) PyFloat_AsDouble(PyList_GetItem(row, 0));
      dest[13] = (float) PyFloat_AsDouble(PyList_GetItem(row, 1));
      dest[14] = (float) PyFloat_AsDouble(PyList_GetItem(row, 2));
      dest[15] = (float) PyFloat_AsDouble(PyList_GetItem(row, 3));
    }
  }
}
