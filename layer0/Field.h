

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
#define Ffloat3(f,a,b,c) (*((float*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2])))

#define Ffloat3p(f,a,b,c) ((float*)((f)->data + \
                                   (a)*(f)->stride[0] + \
                                   (b)*(f)->stride[1] + \
                                   (c)*(f)->stride[2]))

#define Ffloat4(f,a,b,c,d) (*((float*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2] + \
                                     (d)*(f)->stride[3])))

#define Ffloat4p(f,a,b,c,d) ((float*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2] + \
                                     (d)*(f)->stride[3]))

#define Fint3(f,a,b,c) (*((int*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2])))

#define Fint3p(f,a,b,c) ((int*)((f)->data + \
                                   (a)*(f)->stride[0] + \
                                   (b)*(f)->stride[1] + \
                                   (c)*(f)->stride[2]))

#define Fint4(f,a,b,c,d) (*((int*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2] + \
                                     (d)*(f)->stride[3])))

#define Fint4p(f,a,b,c,d) ((int*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2] + \
                                     (d)*(f)->stride[3]))

#define Fvoid4p(f,a,b,c,d) ((void*)((f)->data + \
                                     (a)*(f)->stride[0] + \
                                     (b)*(f)->stride[1] + \
                                     (c)*(f)->stride[2] + \
                                     (d)*(f)->stride[3]))

CField *FieldNew(PyMOLGlobals * G, int *dim, int n_dim, unsigned int base_size, int type);
void FieldZero(CField * I);
void FieldFree(CField * I);
float FieldInterpolatef(CField * I, int a, int b, int c, float x, float y, float z);
void FieldInterpolate3f(CField * I, int *locus, float *fract, float *result);

PyObject *FieldAsPyList(CField * I);
CField *FieldNewFromPyList(PyMOLGlobals * G, PyObject * list);
CField *FieldNewCopy(PyMOLGlobals * G, CField * src);
int FieldSmooth3f(CField * I);

float* FieldSample(CField * I, int skip);

#endif
