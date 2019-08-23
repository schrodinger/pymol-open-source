

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
#ifndef _H_Sphere
#define _H_Sphere

#include"Vector.h"
#include"PyMOLGlobals.h"

typedef struct SphereRec {
  float *area;
  Vector3f *dot;
  int *StripLen;
  int *Sequence;
  int NStrip, NVertTot;
  int nDot;
  int *Tri;
  int NTri;
  int *Mesh;
  int NMesh;
} SphereRec, *SphereRecPtr;

#define NUMBER_OF_SPHERE_LEVELS 5

struct _CSphere {
  SphereRecPtr Sphere[NUMBER_OF_SPHERE_LEVELS];
  SphereRec *Array;
};

void SphereInit(PyMOLGlobals * G);
void SphereFree(PyMOLGlobals * G);
void SphereRender(PyMOLGlobals * G, int level, const float *centroid, const float *color, float alpha, float radius);

SphereRec* GetSpheroidSphereRec(PyMOLGlobals*);

#endif
