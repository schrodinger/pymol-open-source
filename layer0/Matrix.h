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

/* WARNING - MAJOR GOTCHA!  

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
   module to Algebra.h/.c (for linear algebra and fitting).

*/

void MatrixRotation44f( float *m44, const float angle, const float x,const float y,const float z);

void MatrixDump44f(float *m,char *prefix);

void MatrixTransform44f3f(float *m, float *q,float *p);
void MatrixTransform44f4f(float *m, float *q,float *p);
void MatrixInvRotate44f3f(float *m, float *q,float *p);

void MatrixTransform44fAs33f3f(float *p, float *m, float *q);
void MatrixTransform44fn( unsigned int n, float *q, const float m[16], float *p);
int MatrixEigensolve33d(double *a, double *wr, double *wi, double *v);

void MatrixLoadIdentity44f(float *m);
void MatrixTranslate44f3f( float *m, const float x,const float y,const float z);
void MatrixRotate44f3f( float *m, const float angle, const float x,const float y,const float z);
void MatrixMultiply44f( const float *mat, float *product );
int  MatrixInvert44f( const float *m, float *out );

#define MatrixInvTransform3f MatrixInvRotate44f3f
#define MatrixTransform3f MatrixTransform44fAs33f3f

/* WARNING - routines below use a "non-conventional" 4x4 
  pre-translation - tranformation - and post-translation (TTT) matrix */

void MatrixApplyTTTfn3f(unsigned int n, float *q, const float m[16], float *p );
float MatrixFitRMS(int n,float *v1,float *v2,float *wt,float *ttt);
float MatrixGetRMS(int n,float *v1,float *v2,float *wt);
int *MatrixFilter(float cutoff,int window,int n_pass,int nv,float *v1,float *v2);

#endif

