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

#include"os_predef.h"
#include"os_python.h"

#include"os_std.h"

#include"MemoryDebug.h"
#include"Base.h"
#include"PConv.h"
#include"P.h"
#include"Util.h"

PyObject *PConvAutoNone(PyObject *result) /* automatically own Py_None */
{
  if(result==Py_None)
    Py_INCREF(result);
  else if(result==NULL) {
    result=Py_None;
    Py_INCREF(result);
  } 
  return(result);
}

/* Error-checked utility routines */

int PConvPyListToExtent(PyObject *obj,float *mn,float *mx) /* [[min_x,min_y,min_z],
                                                              [max_x,max_y,max_z]] */
{
  int ok=false;
  PyObject *t1,*t2;
  if(!obj) {
    ok=false;
  } else if(PyList_Check(obj))
    if(PyList_Size(obj)==2) {
      t1 = PyList_GetItem(obj,0);
      t2 = PyList_GetItem(obj,1);
      if(PConvPyListToFloatArrayInPlace(t1,mn,3)&&
         PConvPyListToFloatArrayInPlace(t2,mx,3))
        ok=true;
    }
  return(ok);
}

int PConvAttrToIntArrayInPlace(PyObject *obj,char *attr,int *f,int ll)
{
  int ok=true;
  PyObject *tmp;
  if(!obj) {
    ok=false;
  } else if(PyObject_HasAttrString(obj,attr)) {
    tmp = PyObject_GetAttrString(obj,attr);
    ok = PConvPyListToIntArrayInPlace(tmp,f,ll);
    Py_DECREF(tmp);
  } else {
    ok=false;
  }
  return(ok);
}

int PConvAttrToFloatArrayInPlace(PyObject *obj,char *attr,float *f,int ll)
{
  int ok=true;
  PyObject *tmp;
  if(!obj) {
    ok=false;
  } else if(PyObject_HasAttrString(obj,attr)) {
    tmp = PyObject_GetAttrString(obj,attr);
    ok = PConvPyListToFloatArrayInPlace(tmp,f,ll);
    Py_DECREF(tmp);
  } else {
    ok=false;
  }
  return(ok);
}

int PConvAttrToStrMaxLen(PyObject *obj,char *attr,char *str,int ll)
{
  int ok=true;
  PyObject *tmp;
  if(!obj) {
    ok=false;
  } else if(PyObject_HasAttrString(obj,attr)) {
    tmp = PyObject_GetAttrString(obj,attr);
    ok = PConvPyObjectToStrMaxLen(tmp,str,ll);
    Py_DECREF(tmp);
  } else {
    ok=false;
  }
  return(ok);
}

int PConvAttrToPtr(PyObject *obj,char *attr,void **cobj)
{
  PyObject *tmp;
  int ok=true;
  if(!obj) {
    ok=false;
  } else if(PyObject_HasAttrString(obj,attr)) {
    tmp = PyObject_GetAttrString(obj,attr);
    ok = PConvCObjectToPtr(tmp,cobj);
    Py_DECREF(tmp);
  } else {
    ok = false;
  }
  return(ok);
}

int PConvCObjectToPtr(PyObject *obj,void **ptr) {
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyCObject_Check(obj))
    ok=false;
  else
    (*ptr) = PyCObject_AsVoidPtr(obj);
  return(ok);
}

int PConvPyStrToStrPtr(PyObject *obj,char **ptr)
{
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyString_Check(obj)) {
    ok=false;
  }
  if(ok) 
    *ptr=PyString_AsString(obj);
  return(ok);
}

int PConvPyStrToStr(PyObject *obj,char *ptr,int size)
{
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyString_Check(obj)) {
    ok=false;
    if(size) *ptr=0;
  }
  else {
    UtilNCopy(ptr,PyString_AsString(obj),size);
  }
  return(ok);
}

int PConvPyIntToInt(PyObject *obj,int *ptr)
{
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyInt_Check(obj)) {
    ok=false;
  }
  else {
    *ptr = PyInt_AsLong(obj);
  }
  return(ok);
}

int PConvPyFloatToFloat(PyObject *obj,float *ptr)
{
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyFloat_Check(obj)) {
    ok=false;
  } else {
    *ptr = (float)PyFloat_AsDouble(obj);
  }
  return(ok);
}

int PConvPyIntToChar(PyObject *obj,char *ptr)
{
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyInt_Check(obj)) {
    ok=false;
  } else {
    *ptr = (char)PyInt_AsLong(obj);
  }
  return(ok);
}

/* == end == */


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
  else if(PyString_Check(object)) {
    st = PyString_AsString(object);
    strncpy(value,st,ln);
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      st = PyString_AsString(tmp);
      strncpy(value,st,ln);
      Py_DECREF(tmp);
    } else
      result=0;
  }
  if(ln>0)
    value[ln]=0;
  else
    value[0]=0;
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
  } else {
    tmp = PyObject_Str(object);
    if(tmp) {
      st = PyString_AsString(tmp);
      strncpy(value,st,ln);
      Py_DECREF(tmp);
    } else
      result=0;
  }
  if(ln>0)
    value[ln]=0;
  else
    value[0]=0;
  UtilCleanStr(value);
  return(result);
}

PyObject *PConvFloatToPyDictItem(PyObject *dict,char *key,float f)
{
  PyObject *tmp;
  tmp = PyFloat_FromDouble((double)f);
  PyDict_SetItemString(dict,key,tmp);
  Py_XDECREF(tmp); 
  return(tmp);
}

PyObject *PConvIntToPyDictItem(PyObject *dict,char *key,int i)
{
  PyObject *tmp;
  tmp = PyInt_FromLong(i);
  PyDict_SetItemString(dict,key,tmp);
  Py_XDECREF(tmp); 
  return(tmp);
}

PyObject *PConvStringToPyDictItem(PyObject *dict,char *key,char *f)
{
  PyObject *tmp;
  tmp = PyString_FromString(f);
  PyDict_SetItemString(dict,key,tmp);
  Py_XDECREF(tmp); 
  return(tmp);
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
  int ok=true;
  float *ff;
  if(!obj) {
    *f=NULL;
    ok=false;
  } else if(!PyList_Check(obj)) {
    *f=NULL;
    ok=false;
  } else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else 
      ok=l;
    (*f) = Alloc(float,l);
    ff = (*f);
    for(a=0;a<l;a++)
      *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj,a));
    
  }
  return(ok);
}

int PConvPyListToFloatVLA(PyObject *obj,float **f)
{
  int a,l;
  float *ff;
  int ok=true;
  if(!obj) {
    *f=NULL;
    ok=false;
  } else if(!PyList_Check(obj)) {
    *f=NULL;
    ok=false;
  } else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else 
      ok=l;
    (*f) = VLAlloc(float,l);
    ff = (*f);
    for(a=0;a<l;a++)
      *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj,a));
    VLASize((*f),float,l);
  }
  return(ok);
}

int PConvPyList3ToFloatVLA(PyObject *obj,float **f)
{
  int a,b,l;
  float *ff;
  PyObject *triple;
  int ok=true;
  if(!obj) {
    *f=NULL;
    ok=false;
  } else if(!PyList_Check(obj)) {
    *f=NULL;
    ok=false;
  } else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else 
      ok=l;
    (*f) = VLAlloc(float,l*3);
    ff = (*f);
    for(a=0;a<l;a++) {
      triple = PyList_GetItem(obj,a);
      ok = PyList_Check(triple);
      if(ok) ok = (PyList_Size(triple)==3);
      if(ok) {
        for(b=0;b<3;b++)
          *(ff++) = PyFloat_AsDouble(PyList_GetItem(triple,b));
      } else {
        ok=false;
        break;
      }
    }
    VLASize((*f),float,l*3);
  }
  return(ok);
}

int PConvPyListToIntArray(PyObject *obj,int **f)
{
  int a,l;
  int *ff;
  int ok=true;
  if(!obj) {
    *f=NULL;
    l=0;
  } else if(!PyList_Check(obj)) {
    *f=NULL;
    ok=false;
  } else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else
      ok=l;
    (*f) = Alloc(int,l);
    ff = (*f);
    for(a=0;a<l;a++)
      *(ff++) = PyInt_AsLong(PyList_GetItem(obj,a));
  }
  return(ok);
}

int PConvPyListToIntVLA(PyObject *obj,int **f)
{
  int a,l;
  int *ff;
  int ok=true;
  if(!obj) {
    *f=NULL;
    l=0;
  } else if(!PyList_Check(obj)) {
    *f=NULL;
    ok=false;
  } else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else
      ok=l;
    (*f) = VLAlloc(int,l);
    ff = (*f);
    for(a=0;a<l;a++)
      *(ff++) = PyInt_AsLong(PyList_GetItem(obj,a));
  }
  return(ok);
}

int PConvPyListToFloatArrayInPlace(PyObject *obj,float *ff,int ll)
{
  int ok = true;
  int a,l;
  if(!obj) { 
    ok=false;
  } else if(!PyList_Check(obj)) {
    ok=false;
  } else {
    l=PyList_Size(obj);
    if (l!=ll) 
      ok=false;
    else {
      if(!l)
        ok=-1;
      else
        ok=l;
      for(a=0;a<l;a++)
        *(ff++) = (float)PyFloat_AsDouble(PyList_GetItem(obj,a));
    }
    /* NOTE ASSUMPTION! */
  }
  return(ok);
}

int PConvPyListToIntArrayInPlace(PyObject *obj,int *ii,int ll)
{
  int ok = true;
  int a,l;
  if(!obj) 
    ok=false;
  else if(!PyList_Check(obj)) 
    ok=false;
  else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else
      ok=l;
    if (l!=ll) 
      ok=false;
    else 
      for(a=0;a<l;a++) 
        *(ii++) = PyInt_AsLong(PyList_GetItem(obj,a)); 
    /* NOTE ASSUMPTION! */
  }
  return(ok);
}

int PConvPyListToIntArrayInPlaceAutoZero(PyObject *obj,int *ii,int ll)
{
  int ok = true;
  int a,l;
  if(!obj) 
    ok=false;
  else if(!PyList_Check(obj)) 
    ok=false;
  else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else
      ok=l;
    for(a=0;a<l;a++) 
      *(ii++) = PyInt_AsLong(PyList_GetItem(obj,a)); 
    while(l<ll) {
      *(ii++)=0;
      l++;
    }
  }
  return(ok);
}

int PConvPyListToSIntArrayInPlaceAutoZero(PyObject *obj,short int *ii,int ll)
{
  int ok = true;
  int a,l;
  if(!obj) 
    ok=false;
  else if(!PyList_Check(obj)) 
    ok=false;
  else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else
      ok=l;
    for(a=0;a<l;a++) 
      *(ii++) = (short int)PyInt_AsLong(PyList_GetItem(obj,a)); 
    while(l<ll) {
      *(ii++)=0;
      l++;
    }
  }
  return(ok);
}

int PConvPyListToFloatArrayInPlaceAutoZero(PyObject *obj,float *ii,int ll)
{
  int ok = true;
  int a,l;
  if(!obj) 
    ok=false;
  else if(!PyList_Check(obj)) 
    ok=false;
  else {
    l=PyList_Size(obj);
    if(!l)
      ok=-1;
    else
      ok=l;
    for(a=0;a<l;a++) 
      *(ii++) = (float)PyFloat_AsDouble(PyList_GetItem(obj,a)); 
    while(l<ll) {
      *(ii++)=0.0f;
      l++;
    }
  }
  return(ok);
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

PyObject *PConvIntArrayToPyList(int *f,int l)
{
  int a;
  PyObject *result = Py_None;
  result=PyList_New(l);
  for(a=0;a<l;a++) 
    PyList_SetItem(result,a,PyInt_FromLong(*(f++)));
  return(result);
}

PyObject *PConvSIntArrayToPyList(short int *f,int l)
{
  int a;
  PyObject *result = Py_None;
  result=PyList_New(l);
  for(a=0;a<l;a++) 
    PyList_SetItem(result,a,PyInt_FromLong(*(f++)));
  return(result);
}

PyObject *PConv3DIntArrayTo3DPyList(int ***array,int *dim)
{
  int a,b,c;
  PyObject *result,*pyB,*pyC;
  result = PyList_New(dim[0]);
  for(a=0;a<dim[0];a++) {
    pyB = PyList_New(dim[1]);
    PyList_SetItem(result,a,pyB);
    for(b=0;b<dim[1];b++) {
      pyC = PyList_New(dim[2]);
      PyList_SetItem(pyB,b,pyC);
      for(c=0;c<dim[2];c++) {
        PyList_SetItem(pyC,c,PyInt_FromLong(array[a][b][c]));
      }
    }
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

int PConvPyListToStringVLA(PyObject *obj,char **vla_ptr)
{
  int a,l,ll;
  char *vla = NULL,*p,*q;
  PyObject *i;
  if(obj) 
    if(PyList_Check(obj)) {
      l=PyList_Size(obj);
      ll=0;
      for(a=0;a<l;a++) {
        i = PyList_GetItem(obj,a);
        if (PyString_Check(i)) {
          ll+=strlen(PyString_AsString(i))+1;
        }
      }
      vla=VLAlloc(char,ll);
      VLASize(vla,char,ll);
      q=vla;
      for(a=0;a<l;a++) {
        i = PyList_GetItem(obj,a);
        if (PyString_Check(i)) {
          p=PyString_AsString(i);
          while(*p) 
            *(q++)=*(p++);
          *(q++)=0;
        }
      }
    }
  (*vla_ptr)=vla;
  return(vla&&1);
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



