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

#ifndef _PYMOL_NOPY

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

int PConvPyListToStrVLAList(PyObject *obj,char **vla, int *n_str)
{
  int ok=false;
  PyObject *t;
  int n_st = 0, n_ch = 0, nn_ch, l, i;
  if(!*vla)
    *vla = VLAlloc(char,10);
  if((!obj)||(!*vla)) {
    ok=false;
  } else if(PyList_Check(obj)) {
    n_st = PyList_Size(obj);
    ok = true;
    for(i=0;i<n_st;i++) {
      t = PyList_GetItem(obj,i);
      if(PyString_Check(t)) {
        l = PyString_Size(t);
        nn_ch = n_ch+l+1;
        VLACheck(*vla, char, nn_ch);
        UtilNCopy((*vla)+n_ch,PyString_AsString(t),l+1);
        n_ch = nn_ch;
      } else {
        VLACheck(*vla, char, n_ch+1);
        (*vla)[n_ch] = 0;
        n_ch++;
      }
    }
  }
  *n_str = n_st;
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

int PConvPyStrToLexRef(PyObject *obj,OVLexicon *lex,int *lex_ref)
{
  int ok=true;
  if(!obj) {
    ok=false;
  } else if (!PyString_Check(obj)) {
    ok=false;
  } else {
    char *ptr = PyString_AsString(obj);
    if(!ptr) {
      ok=false;
    } else {
      OVreturn_word result = OVLexicon_GetFromCString(lex,ptr);
      if(OVreturn_IS_OK(result)) {
        *lex_ref = result.word;
      } else {
        ok=false;
      }
    }
  }
  return ok;
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
	  if(!PyLong_Check(obj)) {
          ok=false;	
	  } else {
		  *ptr = (int)PyLong_AsLongLong(obj);
	  }
  } else {
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
    if(!PyLong_Check(obj)) {
          ok=false;	
	  } else {
		  *ptr = (char)PyLong_AsLongLong(obj);
	  }
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
    (*value) = (float)PyFloat_AsDouble(object);
  } else if(PyInt_Check(object)) {
    (*value) = (float)PyInt_AsLong(object);
  } else if(PyLong_Check(object)) {
    (*value) = (float)PyLong_AsLongLong(object);
  }else {
    tmp = PyNumber_Float(object);
    if(tmp) {
      (*value) = (float)PyFloat_AsDouble(tmp);
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
  } else if(PyLong_Check(object)) {
    (*value) = (int)PyLong_AsLongLong(object);
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

int PConvPyObjectToChar(PyObject *object,char *value)
{
  int result = true;
  PyObject *tmp;
  if(!object)
    result=false;
  else   if(PyInt_Check(object)) {
    (*value) = (char)PyInt_AsLong(object);
  } else if(PyLong_Check(object)) {
    (*value) = (char)PyLong_AsLongLong(object);
  } else {
    tmp = PyNumber_Int(object);
    if(tmp) {
      (*value) = (char)PyInt_AsLong(tmp);
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
      *(ff++) = (float)PyFloat_AsDouble(PyList_GetItem(obj,a));
    
  }
  return(ok);
}

int PConvPyListToFloatVLANoneOkay(PyObject *obj,float **f)
{
  int a,l;
  float *ff;
  int ok=true;
  if(!obj) {
    *f=NULL;
    ok=false;
  } else if(obj==Py_None) {
    *f=NULL;
    ok=true;
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
      *(ff++) = (float)PyFloat_AsDouble(PyList_GetItem(obj,a));
    VLASize((*f),float,l);
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
      *(ff++) = (float)PyFloat_AsDouble(PyList_GetItem(obj,a));
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
          *(ff++) = (float)PyFloat_AsDouble(PyList_GetItem(triple,b));
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

int PConvPyListToDoubleArray(PyObject *obj,double **f)
{
  int a,l;
  double *ff;
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
    (*f) = Alloc(double,l);
    ff = (*f);
    for(a=0;a<l;a++)
      *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj,a));
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

int PConvPyListToDoubleArrayInPlace(PyObject *obj,double *ff,int ll)
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
        *(ff++) = PyFloat_AsDouble(PyList_GetItem(obj,a));
    }
    /* NOTE ASSUMPTION! */
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
    for(a=0;(a<l)&&(a<ll);a++) 
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
    for(a=0;(a<l)&&(a<ll);a++) 
      *(ii++) = (short int)PyInt_AsLong(PyList_GetItem(obj,a)); 
    while(l<ll) {
      *(ii++)=0;
      l++;
    }
  }
  return(ok);
}

int PConvPyListToSCharArrayInPlaceAutoZero(PyObject *obj,signed char *ii,int ll)
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
    for(a=0;(a<l)&&(a<ll);a++) 
      *(ii++) = (signed char)PyInt_AsLong(PyList_GetItem(obj,a)); 
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
    for(a=0;(a<l)&&(a<ll);a++) 
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

PyObject *PConvFloatArrayToPyListNullOkay(float *f,int l)
{
  int a;
  PyObject *result = Py_None;
  if(!f) {
    result = PConvAutoNone(NULL);
  } else {
    result=PyList_New(l);
    for(a=0;a<l;a++) 
      PyList_SetItem(result,a,PyFloat_FromDouble((double)*(f++)));
  }
  return(result);
}

PyObject *PConvDoubleArrayToPyList(double *f,int l)
{
  int a;
  PyObject *result = Py_None;
  result=PyList_New(l);
  for(a=0;a<l;a++) 
    PyList_SetItem(result,a,PyFloat_FromDouble(*(f++)));
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

PyObject *PConvSCharArrayToPyList(signed char *f,int l)
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
  return(PConvAutoNone(result));
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

int PConvPyListToLabPosVLA(PyObject *obj,LabPosType **vla_ptr)
{
  int a,l;
  int ok=true;
  LabPosType *vla = NULL,*q;
  PyObject *i;
  if(obj) 
    if(PyList_Check(obj)) {
      l=PyList_Size(obj);
      vla=VLACalloc(LabPosType,l);
      q=vla;
      for(a=0;a<l;a++) {
        i = PyList_GetItem(obj,a);
        if (PyList_Check(i) && (PyList_Size(i)==7)) {
          ok = ok && PConvPyIntToInt(PyList_GetItem(i,0),&q->mode) && 
            PConvPyFloatToFloat(PyList_GetItem(i,1),q->pos) && 
            PConvPyFloatToFloat(PyList_GetItem(i,2),q->pos+1) && 
            PConvPyFloatToFloat(PyList_GetItem(i,3),q->pos+2) && 
            PConvPyFloatToFloat(PyList_GetItem(i,4),q->offset) && 
            PConvPyFloatToFloat(PyList_GetItem(i,5),q->offset+1) &&
            PConvPyFloatToFloat(PyList_GetItem(i,6),q->offset+2);
        } else {
          VLAFreeP(vla); /* just in case... */
          vla = NULL;
          break;
        }
        q++;
      }
    }
  if(!ok&&(!vla)) {
    VLAFreeP(vla);
  }
  (*vla_ptr)=vla;
  return(ok);
}

PyObject *PConvLabPosVLAToPyList(LabPosType *vla,int l)
{ /* TO DO error handling */
  int a;
  LabPosType *p = vla;
  PyObject *result = Py_None;
  if(p) {
    PyObject *item;
    result=PyList_New(l); 
    for(a=0;a<l;a++) {
      item = PyList_New(7);
      if(item) {
        PyList_SetItem(item,0,PyInt_FromLong(p->mode));
        PyList_SetItem(item,1,PyFloat_FromDouble((double)p->pos[0]));
        PyList_SetItem(item,2,PyFloat_FromDouble((double)p->pos[1]));
        PyList_SetItem(item,3,PyFloat_FromDouble((double)p->pos[2]));
        PyList_SetItem(item,4,PyFloat_FromDouble((double)p->offset[0]));
        PyList_SetItem(item,5,PyFloat_FromDouble((double)p->offset[1]));
        PyList_SetItem(item,6,PyFloat_FromDouble((double)p->offset[2]));
        PyList_SetItem(result,a,item);
      }
      
      p++;
    }
  }
  return(PConvAutoNone(result));
}

void PConv44PyListTo44f(PyObject *src,float *dest) /* note lost of precision */
{
  PyObject *row;
  if(src&&dest&&PyList_Check(src)) {
      row = PyList_GetItem(src,0);
      if(row&&PyList_Check(row)) {
        dest[ 0]=(float)PyFloat_AsDouble(PyList_GetItem(row,0));
        dest[ 1]=(float)PyFloat_AsDouble(PyList_GetItem(row,1));
        dest[ 2]=(float)PyFloat_AsDouble(PyList_GetItem(row,2));
        dest[ 3]=(float)PyFloat_AsDouble(PyList_GetItem(row,3));
      }
      row = PyList_GetItem(src,1); 
      if(row&&PyList_Check(row)) {
        dest[ 4]=(float)PyFloat_AsDouble(PyList_GetItem(row,0));
        dest[ 5]=(float)PyFloat_AsDouble(PyList_GetItem(row,1));
        dest[ 6]=(float)PyFloat_AsDouble(PyList_GetItem(row,2));
        dest[ 7]=(float)PyFloat_AsDouble(PyList_GetItem(row,3));
      }
      row = PyList_GetItem(src,2);
      if(row&&PyList_Check(row)) {
        dest[ 8]=(float)PyFloat_AsDouble(PyList_GetItem(row,0));
        dest[ 9]=(float)PyFloat_AsDouble(PyList_GetItem(row,1));
        dest[10]=(float)PyFloat_AsDouble(PyList_GetItem(row,2));
        dest[11]=(float)PyFloat_AsDouble(PyList_GetItem(row,3));
      }
      row = PyList_GetItem(src,3);    
      if(row&&PyList_Check(row)) {
        dest[12]=(float)PyFloat_AsDouble(PyList_GetItem(row,0));
        dest[13]=(float)PyFloat_AsDouble(PyList_GetItem(row,1));
        dest[14]=(float)PyFloat_AsDouble(PyList_GetItem(row,2));
        dest[15]=(float)PyFloat_AsDouble(PyList_GetItem(row,3));
      }
  }
}

#else
typedef int this_file_is_no_longer_empty;
#endif


