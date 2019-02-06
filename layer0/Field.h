

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#ifndef _H_Field
#define _H_Field

#include"os_python.h"
#include"PyMOLGlobals.h"

#define cFieldFloat 0
#define cFieldInt 1
#define cFieldOther 2

typedef struct {
  int type;
  char *data;
  unsigned int *dim;
  unsigned int *stride;
  int n_dim;
  unsigned int size;
  unsigned int base_size;
} CField;

/* accessors for getting data from a field */

#define F3p(f,a,b,c) ((f)->data + \
        (a)*(f)->stride[0] + \
        (b)*(f)->stride[1] + \
        (c)*(f)->stride[2])

#define F4p(f,a,b,c,d) ((f)->data + \
        (a)*(f)->stride[0] + \
        (b)*(f)->stride[1] + \
        (c)*(f)->stride[2] + \
        (d)*(f)->stride[3])

#define Ffloat3p(f,a,b,c) ((float*)F3p(f,a,b,c))
#define Ffloat3(f,a,b,c) (*(Ffloat3p(f,a,b,c)))

#define Ffloat4p(f,a,b,c,d) ((float*)F4p(f,a,b,c,d))
#define Ffloat4(f,a,b,c,d) (*(Ffloat4p(f,a,b,c,d)))

#define Fint3p(f,a,b,c) ((int*)F3p(f,a,b,c))
#define Fint3(f,a,b,c) (*(Fint3p(f,a,b,c)))

#define Fint4p(f,a,b,c,d) ((int*)F4p(f,a,b,c,d))
#define Fint4(f,a,b,c,d) (*(Fint4p(f,a,b,c,d)))

#define Fvoid4p(f,a,b,c,d) ((void*)F4p(f,a,b,c,d))

CField *FieldNew(PyMOLGlobals * G, int *dim, int n_dim, unsigned int base_size, int type);
void FieldZero(CField * I);
void FieldFree(CField * I);
float FieldInterpolatef(CField * I, int a, int b, int c, float x, float y, float z);
void FieldInterpolate3f(CField * I, int *locus, float *fract, float *result);

#define FieldFreeP(ptr) {if(ptr){FieldFree(ptr);ptr=NULL;}}

PyObject *FieldAsNumPyArray(CField * I, short copy);
PyObject *FieldAsPyList(PyMOLGlobals * G, CField * I);
CField *FieldNewFromPyList(PyMOLGlobals * G, PyObject * list);
CField *FieldNewFromPyList_From_List(PyMOLGlobals * G, PyObject * list, int);
CField *FieldNewCopy(PyMOLGlobals * G, const CField * src);
int FieldSmooth3f(CField * I);

#endif
