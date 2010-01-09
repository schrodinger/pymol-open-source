

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
#ifndef _H_Triangle
#define _H_Triangle

#include"Vector.h"
#include"PyMOLGlobals.h"

int *TrianglePointsToSurface(PyMOLGlobals * G, float *v, float *vn, int n,
                             float cutoff, int *nTriPtr, int **stripPtr, float *extent,
                             int cavity_mode);

int TriangleDegenerate(float *v1, float *n1, float *v2, float *n2, float *v3, float *n3);

#endif
