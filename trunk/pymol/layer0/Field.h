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
#ifndef _H_Field
#define _H_Field

#define cFieldFloat 0

typedef struct {
  char *data;
  unsigned int *dim;
  unsigned int *stride;
} CField;

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

CField *FieldNew(unsigned int *dim,int n_dim,unsigned int base_size);
void FieldFree(CField *I);

#endif
