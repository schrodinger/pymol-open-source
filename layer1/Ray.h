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
#ifndef _H_Ray
#define _H_Ray

#include"Basis.h"

#define cRayMaxBasis 10


typedef struct CRay {
  void (*fSphere3fv)(struct CRay *ray,float *v,float r);
  void (*fCylinder3fv)(struct CRay *ray,float *v1,float *v2,float r,float *c1,float *c2);
  void (*fColor3fv)(struct CRay *ray,float *c);
  void (*fTriangle3fv)(struct CRay *ray,
							  float *v1,float *v2,float *v3,
							  float *n1,float *n2,float *n3,
							  float *c1,float *c2,float *c3);
  void (*fTexture)(struct CRay *ray,int mode,float *par);
  void (*fTransparentf)(struct CRay *ray,float t);
  CPrimitive *Primitive;
  int NPrimitive;
  CBasis *Basis;
  int NBasis;
  int *Vert2Prim;
  float CurColor[3];
  float ModelView[16];
  float Volume[6];
  float Range[3];
  int BigEndian;
  int Texture;
  float TextureParam[3];
  float Trans;
  float Random[256];
} CRay;

CRay *RayNew(void);
void RayFree(CRay *I);
void RayPrepare(CRay *I,float v0,float v1,float v2,float v3,float v4,float v5,float *mat);
void RayRender(CRay *I,int width,int height,unsigned int *image,float front,float back,double timing);
void RayRenderPOV(CRay *I,int width,int height,char **headerVLA,char **charVLA,float front,float back,float fov);

#endif



