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
#define cPrimTriangle 3
#define cPrimSausage 4

#define cCylCapNone 0
#define cCylCapFlat 1
#define cCylCapRound 2

typedef struct {
  int type,vert,cull;
  float v1[3],v2[3],v3[3];
  float n0[3],n1[3],n2[3],n3[3];
  float c1[3],c2[3],c3[3];
  float r1,l1; 
  float trans;
  int cap1,cap2;
  int texture;
  float texture_param[3];
} CPrimitive;

typedef struct {
  MapType *Map;
  float *Vertex,*Normal,*Precomp;
  float *Radius,*Radius2,MaxRadius,MinVoxel;
  int *Vert2Normal;
  int NVertex;
  int NNormal;
  float LightNormal[3]; /* for lights - this is the direction of the light rays */
  float Color[3]; /* for lights */
  Matrix33f Matrix;
} CBasis;

typedef struct {
  float base[3];
  CPrimitive *prim;
  float impact[3];
  float tri1,tri2;
  float sphere[3]; /* sphere location if reflecting off of one */
  float surfnormal[3]; /* normal of reflecting surface */
  float dist;
  float dotgle;
  float reflect[3];
} RayInfo;

void BasisInit(CBasis *I);
void BasisFinish(CBasis *I);
void BasisMakeMap(CBasis *I,int *vert2prim,CPrimitive *prim,float *clipBox);
void BasisSetupMatrix(CBasis *I);
void BasisGetTriangleNormal(CBasis *I,RayInfo *r,int i,float *fc);
void BasisTrianglePrecompute(float *v1,float *v2,float *v3,float *pre);

int BasisHit(CBasis *I,RayInfo *r,int except,
				 int *vert2prim,CPrimitive *prim,
				 int shadow,float front,float back,
             float excl_trans,int trans_shadows,
             float fudge);

#endif


