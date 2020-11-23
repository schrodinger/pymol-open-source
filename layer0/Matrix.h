

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
#ifndef _H_Matrix
#define _H_Matrix

#include"PyMOLGlobals.h"


/* WARNING: PYMOL's internal matrix handling is a partially-mitigated
   disaster.  There are several different kinds of matrices in use in
   PyMOL, and until recently, it wasn't clear which conventions were
   being used where.  Over time, we WILL fix that.

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

/*========================================================================*/
#define MatrixSetScaleC44f(m, scale) { m[0] = m[5] = m[10] = scale; }

void MatrixGetRotationC44f(float *m44, float angle,
                           float x, float y, float z);

void MatrixTransformC44f3f(const float *m, const float *q, float *p);
void MatrixTransformC44f4f(const float *m, const float *q, float *p);
void MatrixInvTransformC44fAs33f3f(const float *m, const float *q, float *p);

void MatrixTransformC44fAs33f3f(const float *p, const float *m, float *q);
int MatrixEigensolveC33d(PyMOLGlobals * G, const double *a, double *wr, double *wi, double *v);

void MatrixTranslateC44f(float *m, float x, float y, float z);
void MatrixRotateC44f(float *m, float angle, float x, float y, float z);
void MatrixMultiplyC44f(const float *mat, float *product);

void MatrixTransformTTTfN3f(unsigned int n, float *q, const float *m, const float *p);
float MatrixFitRMSTTTf(PyMOLGlobals * G, int n, const float *v1, const float *v2, const float *wt,
                       float *ttt);

float MatrixGetRMS(PyMOLGlobals * G, int n, const float *v1, const float *v2, float *wt);

void MatrixTransformR44fN3f(unsigned int n, float *q, const float *m, const float *p);

int MatrixInvTransformExtentsR44d3f(const double *matrix,
                                    const float *old_min, const float *old_max,
                                    float *new_min, float *new_max);
int MatrixTransformExtentsR44d3f(const double *matrix,
                                 const float *old_min, const float *old_max,
                                 float *new_min, float *new_max);

int MatrixInvertC44f(const float *m, float *out);

int xx_matrix_invert(double *result, const double *input, int size);
int xx_matrix_jacobi_solve(double *e_vec, double *e_val, int *n_rot,
                           const double *input, int size);

#endif
