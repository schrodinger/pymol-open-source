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
#include"os_std.h"
#include"Base.h"

#include"OOMac.h"
#include"PConv.h"

#include"Field.h"

PyObject *FieldAsPyList(CField *I)
{
 PyObject *result = NULL;
 int n_elem;

  /* first, dump the atoms */

  result = PyList_New(7);
  PyList_SetItem(result,0,PyInt_FromLong(I->type));
  PyList_SetItem(result,1,PyInt_FromLong(I->n_dim));
  PyList_SetItem(result,2,PyInt_FromLong(I->base_size));
  PyList_SetItem(result,3,PyInt_FromLong(I->size));
  PyList_SetItem(result,4,PConvIntArrayToPyList((int*)I->dim,I->n_dim));
  PyList_SetItem(result,5,PConvIntArrayToPyList((int*)I->stride,I->n_dim));
  n_elem = I->size/I->base_size;
  switch(I->type) {
  case cFieldInt:
    PyList_SetItem(result,6,PConvIntArrayToPyList((int*)I->data,n_elem));
    break;
  case cFieldFloat:
    PyList_SetItem(result,6,PConvFloatArrayToPyList((float*)I->data,n_elem));
    break;
  default:
    PyList_SetItem(result,6,PConvAutoNone(Py_None));
    break;
  }

#if 0
  int type;
  char *data;
  unsigned int *dim;
  unsigned int *stride;
  int n_dim;
  unsigned int size;
  unsigned int base_size;
#endif

  return(PConvAutoNone(result));  
}

CField *FieldNewFromPyList(PyObject *list)
{
  int ok=true;
  unsigned int n_elem;
  OOAlloc(CField);

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,0),&I->type);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,1),&I->n_dim);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,2),(int*)&I->base_size);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,3),(int*)&I->size);
  if(ok) ok=PConvPyListToIntArray(PyList_GetItem(list,4),(int**)&I->dim);
  if(ok) ok=PConvPyListToIntArray(PyList_GetItem(list,5),(int**)&I->stride);

  if(ok) {
    n_elem = I->size/I->base_size;
    switch(I->type) {
    case cFieldInt:
      ok=PConvPyListToIntArray(PyList_GetItem(list,6),(int**)&I->data);
      break;
    case cFieldFloat:
      ok=PConvPyListToFloatArray(PyList_GetItem(list,6),(float**)&I->data);
      break;
    default:
      I->data=(char*)mmalloc(I->size);
      break;
    }
  }
  if(!ok) {
    OOFreeP(I);
    I=NULL;
  }
  return(I);
}

float FieldInterpolatef(CField *I,int a,int b,int c,float x,float y,float z)
{
  /* basic trilinear interpolation */

  float x1,y1,z1;
  x1=1.0-x;
  y1=1.0-y;
  z1=1.0-z;
  /*
  printf("%8.5f %8.5f %8.5f %8.3f\n %8.5f %8.5f %8.5f %8.5f \n",
         (Ffloat3(I,a  ,b  ,c  )),
         (Ffloat3(I,a+1,b  ,c  )),
         (Ffloat3(I,a  ,b+1,c  )),
         (Ffloat3(I,a  ,b  ,c+1)),
         (Ffloat3(I,a+1,b+1,c  )),
         (Ffloat3(I,a  ,b+1,c+1)),
         (Ffloat3(I,a+1,b  ,c+1)),
         (Ffloat3(I,a+1,b+1,c+1)));*/

  return(
         (Ffloat3(I,a  ,b  ,c  ) *x1*y1*z1) +
         (Ffloat3(I,a+1,b  ,c  ) *x *y1*z1) +
         (Ffloat3(I,a  ,b+1,c  ) *x1*y *z1) +
         (Ffloat3(I,a  ,b  ,c+1) *x1*y1*z ) +
         (Ffloat3(I,a+1,b+1,c  ) *x *y *z1) +
         (Ffloat3(I,a  ,b+1,c+1) *x1*y *z ) +
         (Ffloat3(I,a+1,b  ,c+1) *x *y1*z ) +
         (Ffloat3(I,a+1,b+1,c+1) *x *y *z ));

}

void FieldZero(CField *I)
{
  char *p,*q;
  p=(char*)I->data;
  q=p+I->size;
  MemoryZero(p,q);
}

CField *FieldNew(int *dim,int n_dim,unsigned int base_size,int type)
{
  unsigned int stride;
  int a;

  OOAlloc(CField);
  I->type=type;
  I->base_size=base_size;
  I->stride=(unsigned int*)Alloc(int,n_dim);
  I->dim=(unsigned int*)Alloc(int,n_dim);
  
  stride = base_size;
  for(a=n_dim-1;a>=0;a--) {
    I->stride[a] = stride;
    I->dim[a] = dim[a];
    stride *= dim[a];
  }
  I->data=(char*)mmalloc(stride);
  I->n_dim=n_dim;
  I->size=stride;
  return(I);
}

void FieldFree(CField *I)
{
  if(I) {
    FreeP(I->dim);
    FreeP(I->stride);
    FreeP(I->data);
  }
  OOFreeP(I);
}
