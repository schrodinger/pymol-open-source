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
#ifndef _H_Basis
#define _H_Basis

#include"Map.h"
#include"Vector.h"

#define cPrimSphere 1
#define cPrimCylinder 2

typedef struct {
  int (*Intersect)(float *p,float *n,float **d);
  int type;
  float v1[3],v2[3];
  float c1[3];
  float r1,l1;
} CPrimitive;

typedef struct {
  MapType *Map;
  float *Vertex,*Normal;
  float *Radius,*Radius2,MaxRadius,MinVoxel;
  int *Vert2Normal;
  int NVertex;
  int NNormal;
  float LightNormal[3]; /* for lights - this is the direction of the light rays */
  float Color[3]; /* for lights */
  Matrix33f Matrix;
} CBasis;

void BasisInit(CBasis *I);
void BasisFinish(CBasis *I);
void BasisMakeMap(CBasis *I,int *vert2prim,CPrimitive *prim,float *clipBox);
void BasisSetupMatrix(CBasis *I);

int BasisHit(CBasis *I,float *v,float *minDist,int except,
				 int *vert2prim,CPrimitive *prim,float *sphere,
				 int shadow,float front,float back);
#endif


