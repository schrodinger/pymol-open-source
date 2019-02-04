

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
#include"os_python.h"
#include"os_numpy.h"
#include"os_std.h"
#include"Base.h"

#include"OOMac.h"
#include"PConv.h"

#include"Field.h"
#include"Vector.h"
#include "Setting.h"

/*
 * Get a field as NumPy array. If copy is false, then return an array wrapper
 * around the internal data of field. USE WITH CAUTION, the data pointer will
 * be invalid if the field is freed (e.g. its map object is deleted). If copy
 * is true, then the returned array will have it's own memory (safe).
 */
PyObject *FieldAsNumPyArray(CField * field, short copy)
{
#ifndef _PYMOL_NUMPY
  printf("No numpy support\n");
  return NULL;
#else

  PyObject *result;
  int typenum = -1;
  npy_intp * dims = NULL;

  import_array1(NULL);

  if(field->type == cFieldFloat) {
    switch(field->base_size) {
#ifdef NPY_FLOAT16
      case 2: typenum = NPY_FLOAT16; break;
#endif
      case 4: typenum = NPY_FLOAT32; break;
      case 8: typenum = NPY_FLOAT64; break;
    }
  } else {
    switch(field->base_size) {
      case 1: typenum = NPY_INT8; break;
      case 2: typenum = NPY_INT16; break;
      case 4: typenum = NPY_INT32; break;
      case 8: typenum = NPY_INT64; break;
    }
  }

  if(typenum == -1) {
    printf("error: no typenum for type %d and base_size %d\n",
        field->type, field->base_size);
    return NULL;
  }

  ok_assert(1, dims = pymol::malloc<npy_intp>(field->n_dim));
  copyN(field->dim, dims, field->n_dim);

  if(copy) {
    if((result = PyArray_SimpleNew(field->n_dim, dims, typenum)))
      memcpy(PyArray_DATA((PyArrayObject *)result), field->data, field->size);
  } else {
    result = PyArray_SimpleNewFromData(field->n_dim, dims, typenum, field->data);
  }

  mfree(dims);
  return result;
ok_except1:
  printf("FieldAsNumPyArray failed\n");
  return NULL;
#endif
}

PyObject *FieldAsPyList(PyMOLGlobals * G, CField * I)
{
  PyObject *result = NULL;
  int n_elem;

  int pse_export_version = SettingGetGlobal_f(G, cSetting_pse_export_version) * 1000;
  bool dump_binary = (!pse_export_version || pse_export_version > 1776) && SettingGetGlobal_b(G, cSetting_pse_binary_dump);

  /* first, dump the atoms */

  result = PyList_New(7);
  PyList_SetItem(result, 0, PyInt_FromLong(I->type));
  PyList_SetItem(result, 1, PyInt_FromLong(I->n_dim));
  PyList_SetItem(result, 2, PyInt_FromLong(I->base_size));
  PyList_SetItem(result, 3, PyInt_FromLong(I->size));
  PyList_SetItem(result, 4, PConvIntArrayToPyList((int *) I->dim, I->n_dim));
  PyList_SetItem(result, 5, PConvIntArrayToPyList((int *) I->stride, I->n_dim));
  n_elem = I->size / I->base_size;
  switch (I->type) {
  case cFieldInt:
    PyList_SetItem(result, 6, PConvIntArrayToPyList((int *) I->data, n_elem, dump_binary));
    break;
  case cFieldFloat:
    PyList_SetItem(result, 6, PConvFloatArrayToPyList((float *) I->data, n_elem, dump_binary));
    break;
  default:
    PyList_SetItem(result, 6, PConvAutoNone(Py_None));
    break;
  }

  return (PConvAutoNone(result));

}

CField *FieldNewCopy(PyMOLGlobals * G, const CField * src)
{
  int ok = true;
  OOAlloc(G, CField);

  I->type = src->type;
  I->n_dim = src->n_dim;
  I->base_size = src->base_size;
  I->size = src->size;

  {
    int a;
    I->dim = pymol::malloc<unsigned int>(src->n_dim);
    I->stride = pymol::malloc<unsigned int>(src->n_dim);
    ok = I->dim && I->stride;
    if(ok)
      for(a = 0; a < src->n_dim; a++) {
        I->dim[a] = src->dim[a];
        I->stride[a] = src->stride[a];
      }
  }

  if(ok) {
      ok = ((I->data = pymol::malloc<char>(I->size)) != nullptr);
      if(ok)
        memcpy(I->data, src->data, I->size);
  }
  if(!ok) {
    if(I) {
      FreeP(I->data);
      FreeP(I->dim);
      FreeP(I->stride);
      OOFreeP(I);
    }
    I = NULL;
  }
  return I;
}

CField *FieldNewFromPyList(PyMOLGlobals * G, PyObject * list)
{
  int ok = true;
  int *I_dim = NULL;
  int *I_stride = NULL;

  OOAlloc(G, CField);

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->type);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->n_dim);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 2), (int *) &I->base_size);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 3), (int *) &I->size);
  if(ok)
    ok = PConvPyListToIntArray(PyList_GetItem(list, 4), &I_dim);
  if(ok)
    I->dim = (unsigned int *) I_dim;
  if(ok)
    ok = PConvPyListToIntArray(PyList_GetItem(list, 5), &I_stride);
  if(ok)
    I->stride = (unsigned int *) I_stride;

  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  if(ok) {
    switch (I->type) {
    case cFieldInt:
      {
        int *I_data;
        ok = PConvPyListToIntArray(PyList_GetItem(list, 6), &I_data);
        I->data = (char *) I_data;
      }
      break;
    case cFieldFloat:
      {
        float *I_data;
        ok = PConvPyListToFloatArray(PyList_GetItem(list, 6), &I_data);
        I->data = (char *) I_data;
      }
      break;
    default:
      I->data = pymol::malloc<char>(I->size);
      break;
    }
  }
  if(!ok) {
    OOFreeP(I);
    I = NULL;
  }
  return (I);
}

float FieldInterpolatef(CField * I, int a, int b, int c, float x, float y, float z)
{
  /* basic trilinear interpolation */

  float x1, y1, z1;
  float result1 = 0.0F, result2 = 0.0F;
  float product1, product2;
  x1 = 1.0F - x;
  y1 = 1.0F - y;
  z1 = 1.0F - z;

  {
    char *data = I->data;
    int a_st = I->stride[0];
    int b_st = I->stride[1];
    int c_st = I->stride[2];
    int a_bs = a * a_st;
    int b_bs = b * b_st;
    int c_bs = c * c_st;

    if((product1 = x1 * y1 * z1) != 0.0F)
      result1 += product1 * (*((float *) (data + a_bs + b_bs + c_bs)));
    if((product2 = x * y1 * z1) != 0.0F)
      result2 += product2 * (*((float *) (data + a_bs + a_st + b_bs + c_bs)));
    if((product1 = x1 * y * z1) != 0.0F)
      result1 += product1 * (*((float *) (data + a_bs + b_bs + b_st + c_bs)));
    if((product2 = x1 * y1 * z) != 0.0F)
      result2 += product2 * (*((float *) (data + a_bs + b_bs + c_bs + c_st)));
    if((product1 = x * y * z1) != 0.0F)
      result1 += product1 * (*((float *) (data + a_bs + a_st + b_bs + b_st + c_bs)));
    if((product2 = x1 * y * z) != 0.0F)
      result2 += product2 * (*((float *) (data + a_bs + b_bs + b_st + c_bs + c_st)));
    if((product1 = x * y1 * z) != 0.0F)
      result1 += product1 * (*((float *) (data + a_bs + a_st + b_bs + c_bs + c_st)));
    if((product2 = x * y * z) != 0.0F)
      result2 +=
        product2 * (*((float *) (data + a_bs + a_st + b_bs + b_st + c_bs + c_st)));
  }

  /*
     printf("%8.5f %8.5f %8.5f %8.3f\n %8.5f %8.5f %8.5f %8.5f \n",
     (Ffloat3(I,a  ,b  ,c  )),
     (Ffloat3(I,a+1,b  ,c  )),
     (Ffloat3(I,a  ,b+1,c  )),
     (Ffloat3(I,a  ,b  ,c+1)),
     (Ffloat3(I,a+1,b+1,c  )),
     (Ffloat3(I,a  ,b+1,c+1)),
     (Ffloat3(I,a+1,b  ,c+1)),
     (Ffloat3(I,a+1,b+1,c+1))); */

  return (result1 + result2);

}


/*
float FieldExtrapolatef(CField *I,int a,int b,int c,float x,float y,float z)
{
  float ix = (x > 1.0F) ? 1.0F : 0.0F;
  float iy = (y > 1.0F) ? 1.0F : 0.0F;
  float iz = (z > 1.0F) ? 1.0F : 0.0F;

  x -= ix;
  y -= iy;
  z -= iz;

  {
    float base = FieldInterpolatef(I,a,b,c,ix,iy,iz) - FieldInterpolatef(I,a,b,c,0.0F,0.0F,0.0F);
    float offset = FieldInterpolatef(I,a,b,c,x,y,z);
    return base + offset;
  }
}
*/

void FieldInterpolate3f(CField * I, int *locus, float *fract, float *result)
{
  /* basic trilinear interpolation */

  float x = fract[0];
  float y = fract[1];
  float z = fract[2];
  float x1, y1, z1;
  float result1, result2;
  float product1, product2;
  int a = locus[0];
  int b = locus[1];
  int c = locus[2];
  int d;

  x1 = 1.0F - x;
  y1 = 1.0F - y;
  z1 = 1.0F - z;

  {
    char *data = I->data;
    int a_st = I->stride[0];
    int b_st = I->stride[1];
    int c_st = I->stride[2];
    int d_st = I->stride[3];
    int a_bs = a * a_st;
    int b_bs = b * b_st;
    int c_bs = c * c_st;

    for(d = 0; d < 3; d++) {
      int d_bs = d * d_st;
      result1 = 0.0F;
      result2 = 0.0F;

      if((product1 = x1 * y1 * z1) != 0.0F)
        result1 += product1 * (*((float *) (data + a_bs + b_bs + c_bs + d_bs)));
      if((product2 = x * y1 * z1) != 0.0F)
        result2 += product2 * (*((float *) (data + a_bs + a_st + b_bs + c_bs + d_bs)));
      if((product1 = x1 * y * z1) != 0.0F)
        result1 += product1 * (*((float *) (data + a_bs + b_bs + b_st + c_bs + d_bs)));
      if((product2 = x1 * y1 * z) != 0.0F)
        result2 += product2 * (*((float *) (data + a_bs + b_bs + c_bs + c_st + d_bs)));
      if((product1 = x * y * z1) != 0.0F)
        result1 +=
          product1 * (*((float *) (data + a_bs + a_st + b_bs + b_st + c_bs + d_bs)));
      if((product2 = x1 * y * z) != 0.0F)
        result2 +=
          product2 * (*((float *) (data + a_bs + b_bs + b_st + c_bs + c_st + d_bs)));
      if((product1 = x * y1 * z) != 0.0F)
        result1 +=
          product1 * (*((float *) (data + a_bs + a_st + b_bs + c_bs + c_st + d_bs)));
      if((product2 = x * y * z) != 0.0F)
        result2 +=
          product2 *
          (*((float *) (data + a_bs + a_st + b_bs + b_st + c_bs + c_st + d_bs)));

      result[d] = result1 + result2;
    }
  }
}

int FieldSmooth3f(CField * I)
{
  int a, b, c;
  int na = I->dim[0], nb = I->dim[1], nc = I->dim[2];
  int n_pts = na * nb * nc;
  auto data = (char*) pymol::malloc<float>(n_pts);
  int x, y, z;
  int da, db, dc;
  double tot;
  float cur;
  double inp_sum = 0.0;
  double inp_sumsq = 0.0;
  double out_sum = 0.0;
  double out_sumsq = 0.0;
  int mult, cnt;
  float inp_mean, out_mean, inp_stdev, out_stdev;

  if(data) {
    for(a = 0; a < na; a++)
      for(b = 0; b < nb; b++)
        for(c = 0; c < nc; c++) {
          cur = Ffloat3(I, a, b, c);
          inp_sum += cur;
          inp_sumsq += cur * cur;
          tot = 0.0;
          cnt = 0;
          for(x = -1; x < 2; x++)
            for(y = -1; y < 2; y++)
              for(z = -1; z < 2; z++) {
                da = a + x;
                db = b + y;
                dc = c + z;
                if((da >= 0) && (da < na) &&
                   (db >= 0) && (db < nb) && (dc >= 0) && (dc < nc)) {

                  mult = 1;
                  if(!x)
                    mult = (mult << 1);
                  if(!y)
                    mult = (mult << 1);
                  if(!z)
                    mult = (mult << 1);

                  cur = Ffloat3(I, da, db, dc);
                  tot += mult * cur;
                  cnt += mult;
                }
              }
          tot = tot / cnt;
          *((float *) (data +
                       (a) * I->stride[0] +
                       (b) * I->stride[1] + (c) * I->stride[2])) = (float) tot;
          out_sum += tot;
          out_sumsq += (tot * tot);
        }
    mfree(I->data);
    I->data = data;

    inp_mean = (float) (inp_sum / n_pts);
    inp_stdev = (float) sqrt1d((inp_sumsq - (inp_sum * inp_sum / n_pts)) / (n_pts - 1));

    out_mean = (float) (out_sum / n_pts);
    out_stdev = (float) sqrt1d((out_sumsq - (out_sum * out_sum / n_pts)) / (n_pts - 1));

    if(out_stdev != 0.0F) {
      float scale = inp_stdev / out_stdev;

      for(a = 0; a < na; a++)
        for(b = 0; b < nb; b++)
          for(c = 0; c < nc; c++) {
            cur = Ffloat3(I, a, b, c);
            cur = (cur - out_mean) * scale + inp_mean;
            Ffloat3(I, a, b, c) = cur;
          }
    }
    return 1;
  }
  return 0;
}

void FieldZero(CField * I)
{
  char *p, *q;
  p = (char *) I->data;
  q = p + I->size;
  MemoryZero(p, q);
}

CField *FieldNew(PyMOLGlobals * G, int *dim, int n_dim, unsigned int base_size, int type)
{
  unsigned int stride;
  int a;

  OOAlloc(G, CField);
  I->type = type;
  I->base_size = base_size;
  I->stride = pymol::malloc<unsigned int>(n_dim);
  I->dim = pymol::malloc<unsigned int>(n_dim);

  stride = base_size;
  for(a = n_dim - 1; a >= 0; a--) {
    I->stride[a] = stride;
    I->dim[a] = dim[a];
    stride *= dim[a];
  }
  I->data = pymol::malloc<char>(stride);
  I->n_dim = n_dim;
  I->size = stride;
  return (I);
}

void FieldFree(CField * I)
{
  if(I) {
    FreeP(I->dim);
    FreeP(I->stride);
    FreeP(I->data);
  }
  OOFreeP(I);
}
