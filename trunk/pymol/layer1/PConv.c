
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

#include<Python.h>

#include"os_std.h"

#include"MemoryDebug.h"
#include"Base.h"
#include"PConv.h"
#include"P.h"
#include"Util.h"

int PConvPyObjectToFloat(PyObject *object,float *value)
{
  int result = true;
  PyObject *tmp;
  if(!object)
    result=false;
  else if(PyFloat_Check(object)) {
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
  if(!object)
    result=false;
  else   if(PyInt_Check(object)) {
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
  if(!object)
    result=false;
  else   if(PyString_Check(object)) {
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
  if(!object)
    result=false;
  else if(PyString_Check(object)) {
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
  PyObject *tmp;
  tmp = PyFloat_FromDouble((double)f);
  PyDict_SetItemString(dict,key,tmp);
  Py_XDECREF(tmp); 
}

void PConvIntToPyDictItem(PyObject *dict,char *key,int i)
{
  PyObject *tmp;
  tmp = PyInt_FromLong(i);
  PyDict_SetItemString(dict,key,tmp);
  Py_XDECREF(tmp); 
}

void PConvStringToPyDictItem(PyObject *dict,char *key,char *f)
{
  PyObject *tmp;
  tmp = PyString_FromString(f);
  PyDict_SetItemString(dict,key,tmp);
  Py_XDECREF(tmp); 
}

void PConvFloat3ToPyObjAttr(PyObject *obj,char *attr,float *v)
{
  PyObject *t1,*t2,*t3,*tmp;

  t1 = PyFloat_FromDouble((double)v[0]);
  t2 = PyFloat_FromDouble((double)v[1]);
  t3 = PyFloat_FromDouble((double)v[2]);
  tmp = PyList_New(3);
  if(t1&&t2&&t3&&tmp) {
    PyList_SetItem(tmp,0,t1); /* steals reference */
    PyList_SetItem(tmp,1,t2); /* steals reference */
    PyList_SetItem(tmp,2,t3); /* steals reference */
    PyObject_SetAttrString(obj,attr,tmp);
  }
  Py_XDECREF(tmp); 
}

void PConvInt2ToPyObjAttr(PyObject *obj,char *attr,int *v)
{
  PyObject *t1,*t2,*tmp;

  t1 = PyInt_FromLong((long)v[0]);
  t2 = PyInt_FromLong((long)v[1]);
  tmp = PyList_New(2);
  if(t1&&t2&&tmp) {
    PyList_SetItem(tmp,0,t1); /* steals reference */
    PyList_SetItem(tmp,1,t2); /* steals reference */
    PyObject_SetAttrString(obj,attr,tmp);
  }
  Py_XDECREF(tmp); 
}


void PConvFloatToPyObjAttr(PyObject *obj,char *attr,float f)
{
  PyObject *tmp;
  tmp = PyFloat_FromDouble((double)f);
  PyObject_SetAttrString(obj,attr,tmp);
  Py_DECREF(tmp); 
}

void PConvIntToPyObjAttr(PyObject *obj,char *attr,int i)
{
  PyObject *tmp;
  tmp = PyInt_FromLong(i);
  PyObject_SetAttrString(obj,attr,tmp);
  Py_DECREF(tmp); 
}

void PConvStringToPyObjAttr(PyObject *obj,char *attr,char *f)
{
  PyObject *tmp;
  tmp = PyString_FromString(f);
  PyObject_SetAttrString(obj,attr,tmp);
  Py_DECREF(tmp); 
}

int PConvPyListToFloatArray(PyObject *obj,float **f)
{
  int a,l;
  float *ff;
  l=PyList_Size(obj);
  (*f) = Alloc(float,l);
  ff = (*f);
  for(a=0;a<l;a++)
    *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj,a));
  return(l);
}

int PConvPyListToIntArray(PyObject *obj,int **f)
{
  int a,l;
  int *ff;
  l=PyList_Size(obj);
  (*f) = Alloc(int,l);
  ff = (*f);
  for(a=0;a<l;a++)
    *(ff++) = PyInt_AsLong(PyList_GetItem(obj,a));
  return(l);
}

PyObject *PConvFloatArrayToPyList(float *f,int l)
{
  int a;
  PyObject *result = Py_None;
  result=PyList_New(l);
  for(a=0;a<l;a++) 
    PyList_SetItem(result,a,PyFloat_FromDouble((double)*(f++)));
  return(result);
}

PyObject *PConvFloatVLAToPyList(float *f)
{
  int a,l;
  PyObject *result = Py_None;
  l=VLAGetSize(f);
  result=PyList_New(l);
  for(a=0;a<l;a++) {
    PyList_SetItem(result,a,PyFloat_FromDouble((double)*(f++))); /* set item steals ref */
  }
  return(result);
}

PyObject *PConvIntVLAToPyList(int *f)
{
  int a,l;
  PyObject *result = Py_None;
  l=VLAGetSize(f);
  result=PyList_New(l);
  for(a=0;a<l;a++) 
    PyList_SetItem(result,a,PyInt_FromLong(*(f++)));
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

PyObject *PConvStringVLAToPyList(char *vla)
{
  int a,c,n=0;
  char *p;
  PyObject *result = Py_None;
  p=vla;
  c = VLAGetSize(vla);
  while(c--) { /* count strings */
    if(!*(p++))
      n++;
  }

  result=PyList_New(n); 
  p=vla;
  for(a=0;a<n;a++) {
    PyList_SetItem(result,a,PyString_FromString(p));
    while(*(p++));
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

