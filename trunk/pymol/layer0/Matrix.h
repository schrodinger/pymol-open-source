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
#ifndef _H_Matrix
#define _H_Matrix

#include"PyMOLGlobals.h"


/* WARNING: PYMOL's internal matrix handling is an partially-mitigated
   disaster.  There are many different kinds of matrices in use in
   PyMOL, and until recently, it wasn't clear which conventions were
   being used where.  We're trying to change that...

   Note the coding used on routines:

   (1) C44f: A column-major homogenous float 4x4 (OpenGL compatible)

   (2) C33f: A column-major float 3x3 (OpenGL-like)

   (3) TTTf: A row-major float 4x4 TTT matrix (translate, transform,
   translate) which create nothing but trouble and confusion...

   (4) 33f, R33f: A row-major 3x3 compatible with programmer brains as well
   as with the code in vector.c

   (5) R44d: row-major homogenous 4x4 double precision, compatible
       with programmer brains and with PyMOL's long-term future.   

   (6) R44f: as above, but floating point -- for performance or just
       coding convenience

   Other informative codes in function calls:

   N3f = array of floats, grouped in sets of 3
   
   33Tf = uses the transpose of a 3x3 float matrix provided
*/

/* (old notice)

WARNING - MAJOR GOTCHA!  

   A state of confusion has arisen because of mixed row-major and
   column-major ordering in this module.  THIS NEEDS TO BE
   STRAIGHTENED OUT before any further use is made of the following
   routines.

   The problem is that OpenGL uses one matrix storage format, while C
   and Python follow the other convention.  The situation
   arose because some of the routines below are used to emulate OpenGL
   transformations...(and are in fact based on Mesa code).
   
   Proposed Solution: 

   (1) Assume the C-PYTHON convention of fast-column (row-major?)
   ordering - and specifically rename routines which don't conform as
   fast row (column-major?) "CM44f" matrices.  Clean up and resolve
   problems with existing code that call these routines.
   
   (2) Move all of the 4x4 transformation code into Vector.h and rename this
   module Algebra.h/.c (for linear algebra and fitting).

*/

void MatrixGetRotationC44f( float *m44, const float angle, 
                             const float x,const float y,const float z);

void MatrixTransformC44f3f(float *m, float *q,float *p);
void MatrixTransformC44f4f(float *m, float *q,float *p);
void MatrixInvTransformC44fAs33f3f(float *m, float *q,float *p);

void MatrixTransformC44fAs33f3f(float *p, float *m, float *q);
void MatrixTransformC44fn( unsigned int n, float *q, float *m, float *p);
int MatrixEigensolveC33d(PyMOLGlobals *G,double *a, double *wr, double *wi, double *v);

void MatrixTranslateC44f( float *m, const float x,const float y,const float z);
void MatrixRotateC44f( float *m, const float angle, const float x,const float y,const float z);
void MatrixMultiplyC44f( const float *mat, float *product );

void MatrixTransformTTTfN3f(unsigned int n, float *q, float *m, float *p );
float MatrixFitRMSTTTf(PyMOLGlobals *G,int n,float *v1,float *v2,float *wt,float *ttt);

float MatrixGetRMS(PyMOLGlobals *G,int n,float *v1,float *v2,float *wt);
int *MatrixFilter(float cutoff,int window,int n_pass,int nv,float *v1,float *v2);

void MatrixTransformR44fN3f( unsigned int n, float *q, float *m, float *p );

int MatrixInvTransformExtentsR44d3f(double *matrix, 
                                  float *old_min, float *old_max,
                                  float *new_min, float *new_max);
int MatrixTransformExtentsR44d3f(double *matrix, 
                                 float *old_min, float *old_max,
                                 float *new_min, float *new_max);

typedef long int integer;
typedef double doublereal;

int pymol_rg_(integer *nm, integer *n, doublereal *a, doublereal *wr, 
              doublereal *wi, integer *matz,doublereal *z__,integer *iv1,
              doublereal  *fv1,integer  *ierr);

#endif

