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

void MatrixDump44f(float *m,char *prefix);
void MatrixInvTransform3f(float *vector,float *matrix,float *result);
void MatrixTransform3f(float *vector,float *matrix,float *result);
void MatrixTransform44fn( unsigned int n, float *q, const float m[16], float *p);
int MatrixEigensolve33d(double *a, double *wr, double *wi, double *v);

void MatrixLoadIdentity44f(float *m);
void MatrixTranslate44f3f( float *m, const float x,const float y,const float z);
void MatrixRotate44f3f( float *m, const float angle, const float x,const float y,const float z);
void MatrixMultiply44f( const float *mat, float *product );

/* WARNING - routines below use a "non-conventional" 4x4 
  pre-translation - tranformation - and post-translation (TTT) matrix */

void MatrixApplyTTTfn3f(unsigned int n, float *q, const float m[16], float *p );
float MatrixFitRMS(int n,float *v1,float *v2,float *wt,float *ttt);
float MatrixGetRMS(int n,float *v1,float *v2,float *wt);

#endif

