
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
#include"PUtils.h"

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

