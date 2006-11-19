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
#include"Vector.h"

PyObject *FieldAsPyList(CField *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

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
#endif

}

CField *FieldNewCopy(PyMOLGlobals *G,CField *src)
{  
  int ok=true;
  OOAlloc(G,CField);
  
  I->type = src->type;
  I->n_dim = src->n_dim;
  I->base_size = src->base_size;
  I->size = src->size;

  {
    int a;
    I->dim = Alloc(unsigned int,src->n_dim);
    I->stride = Alloc(unsigned int,src->n_dim);
    ok = I->dim && I->stride;
    if(ok) 
      for(a=0;a<src->n_dim;a++) {
        I->dim[a] = src->dim[a];
        I->stride[a] = src->stride[a];
      }
  }

  if(ok) {
    unsigned int n_elem = I->size/I->base_size;
    switch(I->type) {
    case cFieldInt:
      ok = ((I->data = (char*)Alloc(int,n_elem)) != NULL);
      if(ok) memcpy(I->data, src->data, sizeof(int)*n_elem);
      break;
    case cFieldFloat:
      ok = ((I->data = (char*)Alloc(float,n_elem)) != NULL);
      if(ok) memcpy(I->data, src->data, sizeof(float)*n_elem);
      break;
    default:
      ok = ((I->data = (char*)Alloc(char,I->size)) != NULL);
      if(ok) memcpy(I->data, src->data, I->size);
      break;
    }
  }
  if(!ok) {
    if(I) {
      FreeP(I->data);
      FreeP(I->dim);
      FreeP(I->stride);
      OOFreeP(I);
    }
    I=NULL;
  }
  return I;
}

CField *FieldNewFromPyList(PyMOLGlobals *G,PyObject *list)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  int ok=true;
  unsigned int n_elem;
  int ll;
  OOAlloc(G,CField);

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,0),&I->type);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,1),&I->n_dim);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,2),(int*)&I->base_size);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,3),(int*)&I->size);
  if(ok) ok=PConvPyListToIntArray(PyList_GetItem(list,4),(int**)&I->dim);
  if(ok) ok=PConvPyListToIntArray(PyList_GetItem(list,5),(int**)&I->stride);
  
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */

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
#endif
}


float FieldInterpolatef(CField *I,int a,int b,int c,float x,float y,float z)
{
  /* basic trilinear interpolation */

  register float x1,y1,z1;
  register float result1=0.0F,result2=0.0F;
  register float product1,product2;
  x1=1.0F-x;
  y1=1.0F-y;
  z1=1.0F-z;
  
#if 0
  if((product1 = x1*y1*z1)!=0.0F) result1 += product1 * Ffloat3(I,a  ,b  ,c  );
  if((product2 = x *y1*z1)!=0.0F) result2 += product2 * Ffloat3(I,a+1,b  ,c  );
  if((product1 = x1*y *z1)!=0.0F) result1 += product1 * Ffloat3(I,a  ,b+1,c  );
  if((product2 = x1*y1*z )!=0.0F) result2 += product2 * Ffloat3(I,a  ,b  ,c+1);
  if((product1 = x *y *z1)!=0.0F) result1 += product1 * Ffloat3(I,a+1,b+1,c  );
  if((product2 = x1*y *z )!=0.0F) result2 += product2 * Ffloat3(I,a  ,b+1,c+1);
  if((product1 = x *y1*z )!=0.0F) result1 += product1 * Ffloat3(I,a+1,b  ,c+1);
  if((product2 = x *y *z )!=0.0F) result2 += product2 * Ffloat3(I,a+1,b+1,c+1);
#else
 {
   register char *data = I->data;
   register int a_st = I->stride[0];
   register int b_st = I->stride[1];
   register int c_st = I->stride[2];
   register int a_bs = a * a_st;
   register int b_bs = b * b_st;
   register int c_bs = c * c_st;

   if((product1 = x1*y1*z1)!=0.0F) result1 += product1 * (*((float*)(data + a_bs +        b_bs +        c_bs        )));
   if((product2 = x *y1*z1)!=0.0F) result2 += product2 * (*((float*)(data + a_bs + a_st + b_bs +        c_bs        ))); 
   if((product1 = x1*y *z1)!=0.0F) result1 += product1 * (*((float*)(data + a_bs +        b_bs + b_st + c_bs        ))); 
   if((product2 = x1*y1*z )!=0.0F) result2 += product2 * (*((float*)(data + a_bs +        b_bs +        c_bs + c_st ))); 
   if((product1 = x *y *z1)!=0.0F) result1 += product1 * (*((float*)(data + a_bs + a_st + b_bs + b_st + c_bs        ))); 
   if((product2 = x1*y *z )!=0.0F) result2 += product2 * (*((float*)(data + a_bs +        b_bs + b_st + c_bs + c_st ))); 
   if((product1 = x *y1*z )!=0.0F) result1 += product1 * (*((float*)(data + a_bs + a_st + b_bs +        c_bs + c_st ))); 
   if((product2 = x *y *z )!=0.0F) result2 += product2 * (*((float*)(data + a_bs + a_st + b_bs + b_st + c_bs + c_st ))); 
 }
#endif

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
  
  return(result1+result2);

}

void FieldInterpolate3f(CField *I,int *locus,float *fract, float *result)
{
  /* basic trilinear interpolation */

  register float x = fract[0];
  register float y = fract[1];
  register float z = fract[2];
  register float x1,y1,z1;
  register float result1,result2;
  register float product1,product2;
  register int a = locus[0];
  register int b = locus[1];
  register int c = locus[2];
  register int d;

  x1=1.0F-x;
  y1=1.0F-y;
  z1=1.0F-z;
  
#if 0
  for(d=0;d<3;d++) {
    result1 = 0.0F;
    result2 = 0.0F;

    if((product1 = x1*y1*z1)!=0.0F) result1 += product1 * Ffloat4(I,a  ,b  ,c  , d);
    if((product2 = x *y1*z1)!=0.0F) result2 += product2 * Ffloat4(I,a+1,b  ,c  , d);
    if((product1 = x1*y *z1)!=0.0F) result1 += product1 * Ffloat4(I,a  ,b+1,c  , d);
    if((product2 = x1*y1*z )!=0.0F) result2 += product2 * Ffloat4(I,a  ,b  ,c+1, d);
    if((product1 = x *y *z1)!=0.0F) result1 += product1 * Ffloat4(I,a+1,b+1,c  , d);
    if((product2 = x1*y *z )!=0.0F) result2 += product2 * Ffloat4(I,a  ,b+1,c+1, d);
    if((product1 = x *y1*z )!=0.0F) result1 += product1 * Ffloat4(I,a+1,b  ,c+1, d);
    if((product2 = x *y *z )!=0.0F) result2 += product2 * Ffloat4(I,a+1,b+1,c+1, d);
    
    result[d] = result1+result2;
  }
#else
  {
    register char *data = I->data;
    register int a_st = I->stride[0];
    register int b_st = I->stride[1];
    register int c_st = I->stride[2];
    register int d_st = I->stride[3];
    register int a_bs = a * a_st;
    register int b_bs = b * b_st;
    register int c_bs = c * c_st;


    for(d=0;d<3;d++) {
      register int d_bs = d * d_st;
      result1 = 0.0F;
      result2 = 0.0F;
      
      if((product1 = x1*y1*z1)!=0.0F) result1 += product1 * (*((float*)(data + a_bs +        b_bs +        c_bs        + d_bs )));
      if((product2 = x *y1*z1)!=0.0F) result2 += product2 * (*((float*)(data + a_bs + a_st + b_bs +        c_bs        + d_bs ))); 
      if((product1 = x1*y *z1)!=0.0F) result1 += product1 * (*((float*)(data + a_bs +        b_bs + b_st + c_bs        + d_bs ))); 
      if((product2 = x1*y1*z )!=0.0F) result2 += product2 * (*((float*)(data + a_bs +        b_bs +        c_bs + c_st + d_bs ))); 
      if((product1 = x *y *z1)!=0.0F) result1 += product1 * (*((float*)(data + a_bs + a_st + b_bs + b_st + c_bs        + d_bs ))); 
      if((product2 = x1*y *z )!=0.0F) result2 += product2 * (*((float*)(data + a_bs +        b_bs + b_st + c_bs + c_st + d_bs ))); 
      if((product1 = x *y1*z )!=0.0F) result1 += product1 * (*((float*)(data + a_bs + a_st + b_bs +        c_bs + c_st + d_bs ))); 
      if((product2 = x *y *z )!=0.0F) result2 += product2 * (*((float*)(data + a_bs + a_st + b_bs + b_st + c_bs + c_st + d_bs ))); 

     result[d] = result1+result2;
    }
  }
#endif
}

int FieldSmooth3f(CField *I)
{
  register int a,b,c;
  register int na = I->dim[0],nb = I->dim[1],nc = I->dim[2];
  int n_pts = na*nb*nc;
  char *data = (char*)mmalloc(sizeof(float)*n_pts);
  register int x,y,z;
  register int da,db,dc;
  register double tot;
  register float cur;
  register double inp_sum = 0.0;
  register double inp_sumsq = 0.0;
  register double out_sum = 0.0;
  register double out_sumsq = 0.0;
  register int mult,cnt;
  float inp_mean,out_mean,inp_stdev,out_stdev;

  if(data) {
    for(a=0;a<na;a++)
      for(b=0;b<nb;b++)
        for(c=0;c<nc;c++) {
          cur = Ffloat3(I,a,b,c);
          inp_sum += cur;
          inp_sumsq += cur*cur;
          tot = 0.0;
          cnt = 0;
          for(x=-1;x<2;x++) 
            for(y=-1;y<2;y++)
              for(z=-1;z<2;z++) {
                da = a+x;
                db = b+y;
                dc = c+z;
                if( (da>=0) && (da<na) &&
                    (db>=0) && (db<nb) &&
                    (dc>=0) && (dc<nc)) {
                  
                  mult=1;
                  if(!x) mult = (mult<<1);
                  if(!y) mult = (mult<<1);
                  if(!z) mult = (mult<<1);
                  
                  cur = Ffloat3(I,da,db,dc);
                  tot += mult*cur;
                  cnt += mult;
                }
              }
          tot = tot/cnt;
          *((float*)(data + 
                     (a)*I->stride[0] + \
                     (b)*I->stride[1] + \
                     (c)*I->stride[2])) = (float)tot;
          out_sum += tot;
          out_sumsq += (tot*tot);
        }
    mfree(I->data);
    I->data = data;
  
    inp_mean = (float)(inp_sum/n_pts);
    inp_stdev = (float)sqrt1d((inp_sumsq - (inp_sum*inp_sum/n_pts))/(n_pts-1));

    out_mean = (float)(out_sum/n_pts);
    out_stdev = (float)sqrt1d((out_sumsq - (out_sum*out_sum/n_pts))/(n_pts-1));

    if(out_stdev!=0.0F) {
      float scale = inp_stdev / out_stdev;
      
      for(a=0;a<na;a++)
        for(b=0;b<nb;b++)
          for(c=0;c<nc;c++) {
            cur = Ffloat3(I,a,b,c);
            cur = (cur-out_mean)*scale + inp_mean;
            Ffloat3(I,a,b,c) = cur;
          }
    }
    return 1;
  }
  return 0;
}

void FieldZero(CField *I)
{
  char *p,*q;
  p=(char*)I->data;
  q=p+I->size;
  MemoryZero(p,q);
}

CField *FieldNew(PyMOLGlobals *G,int *dim,int n_dim,unsigned int base_size,int type)
{
  unsigned int stride;
  int a;

  OOAlloc(G,CField);
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
