
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



#include<stdlib.h>
#include<Python.h>
#include<signal.h>
#include<string.h>
#include<sys/types.h>
#include<sys/time.h>
#include<unistd.h>

#include"MemoryDebug.h"
#include"Base.h"
#include"PConv.h"
#include"P.h"

int PConvPyObjectToFloat(PyObject *object,float *value)
{
  int result = true;
  PyObject *tmp;
  if(PyFloat_Check(object)) {
    (*value) = PyFloat_AsDouble(object);
  } else if(PyInt_Check(object)) {
    (*value) = (float)PyInt_AsLong(object);
  } else {
    tmp = PyNumber_Float(object);
    if(tmp) {
      (*value) = PyFloat_AsDouble(tmp);
      Py_DECREF(tmp);
    } else 
      result=false;
  }
  return(result);
}

int PConvPyObjectToInt(PyObject *object,int *value)
{
  int result = true;
  PyObject *tmp;
  if(PyInt_Check(object)) {
    (*value) = (int)PyInt_AsLong(object);
  } else {
    tmp = PyNumber_Int(object);
    if(tmp) {
      (*value) = (int)PyInt_AsLong(tmp);
      Py_DECREF(tmp);
    } else 
      result=false;
  }
  return(result);
}

int PConvPyObjectToStrMaxLen(PyObject *object,char *value,int ln)
{
  char *st;
  PyObject *tmp;
  int result=true;
  if(PyString_Check(object)) {
    st = PyString_AsString(object);
    strncpy(value,st,ln);
    value[ln]=0;
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      st = PyString_AsString(tmp);
      strncpy(value,st,ln);
      value[ln]=0;
      Py_DECREF(tmp);
    } else
      result=0;
  }
  return(result);
}

int PConvPyObjectToStrMaxClean(PyObject *object,char *value,int ln)
{
  char *st;
  PyObject *tmp;
  int result=true;
  if(PyString_Check(object)) {
    st = PyString_AsString(object);
    strncpy(value,st,ln);
    value[ln]=0;
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      st = PyString_AsString(tmp);
      strncpy(value,st,ln);
      value[ln]=0;
      Py_DECREF(tmp);
    } else
      result=0;
  }
  UtilCleanStr(value);
  return(result);
}

void PConvFloatToPyDictItem(PyObject *dict,char *key,float f)
{
  PyObject *fo;
  fo = PyFloat_FromDouble((double)f);
  PyDict_SetItemString(dict,key,fo);
  Py_DECREF(fo); 
}

void PConvIntToPyDictItem(PyObject *dict,char *key,int i)
{
  PyObject *fo;
  fo = PyInt_FromLong(i);
  PyDict_SetItemString(dict,key,fo);
  Py_DECREF(fo); 
}

void PConvStringToPyDictItem(PyObject *dict,char *key,char *f)
{
  PyObject *tmp;
  tmp = PyString_FromString(f);
  PyDict_SetItemString(dict,key,tmp);
  Py_DECREF(tmp); 
}

PyObject *PConvFloatVLAToPyList(float *f)
{
  int a,l;
  PyObject *result = Py_None;
  l=VLAGetSize(f);
  result=PyList_New(l);
  for(a=0;a<l;a++) {
    PyList_SetItem(result,a,PyFloat_FromDouble((double)f[a]));
  }
  return(result);
}

PyObject *PConvStringListToPyList(int l,char **str)
{
  int a;
  PyObject *result = Py_None;
  result=PyList_New(l);
  for(a=0;a<l;a++) {
    PyList_SetItem(result,a,PyString_FromString(str[a]));
  }
  return(result);
}

void PConv44PyListTo44f(PyObject *src,float *dest) /* note lost of precision */
{
  PyObject *row;

  row = PyList_GetItem(src,0);
  dest[ 0]=PyFloat_AsDouble(PyList_GetItem(row,0));
  dest[ 1]=PyFloat_AsDouble(PyList_GetItem(row,1));
  dest[ 2]=PyFloat_AsDouble(PyList_GetItem(row,2));
  dest[ 3]=PyFloat_AsDouble(PyList_GetItem(row,3));
  row = PyList_GetItem(src,1);
  dest[ 4]=PyFloat_AsDouble(PyList_GetItem(row,0));
  dest[ 5]=PyFloat_AsDouble(PyList_GetItem(row,1));
  dest[ 6]=PyFloat_AsDouble(PyList_GetItem(row,2));
  dest[ 7]=PyFloat_AsDouble(PyList_GetItem(row,3));
  row = PyList_GetItem(src,2);
  dest[ 8]=PyFloat_AsDouble(PyList_GetItem(row,0));
  dest[ 9]=PyFloat_AsDouble(PyList_GetItem(row,1));
  dest[10]=PyFloat_AsDouble(PyList_GetItem(row,2));
  dest[11]=PyFloat_AsDouble(PyList_GetItem(row,3));
  row = PyList_GetItem(src,3);    
  dest[12]=PyFloat_AsDouble(PyList_GetItem(row,0));
  dest[13]=PyFloat_AsDouble(PyList_GetItem(row,1));
  dest[14]=PyFloat_AsDouble(PyList_GetItem(row,2));
  dest[15]=PyFloat_AsDouble(PyList_GetItem(row,3));
}

