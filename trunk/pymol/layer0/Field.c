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
