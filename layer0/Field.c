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
#include"OOMac.h"

#include"Field.h"

float FieldInterpolatef(CField *I,int a,int b,int c,float x,float y,float z)
{
  /* basic trilinear interpolation */

  float x1,y1,z1;
  x1=1.0-x;
  y1=1.0-y;
  z1=1.0-z;
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

CField *FieldNew(int *dim,int n_dim,unsigned int base_size)
{
  unsigned int stride;
  int a;

  OOAlloc(CField);
  
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
