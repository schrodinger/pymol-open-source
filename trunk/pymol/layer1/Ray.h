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
  void (*fCustomCylinder3fv)(struct CRay *ray,float *v1,float *v2,float r,float *c1,
                             float *c2,int cap1,int cap2);
  void (*fSausage3fv)(struct CRay *ray,float *v1,float *v2,float r,float *c1,float *c2);
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
  int TTTFlag;
  float TTT[16];
  int Context;
  float AspRatio;
  float PixelRadius;
} CRay;

typedef struct {
  CRay *ray;
  int width,height;
  unsigned int *image;
  float front,back;
  unsigned int fore_mask;
  float *bkrd;
  unsigned int background;
  int phase, n_thread;
  float spec_vector[3];
} CRayThreadInfo;

typedef struct {
  CBasis *basis;
  int *vert2prim;
  CPrimitive *prim;
  float *clipBox;
} CRayHashThreadInfo;


CRay *RayNew(void);
void RayFree(CRay *I);
void RayPrepare(CRay *I,float v0,float v1,float v2,
                float v3,float v4,float v5,float *mat,
                float aspRat,int ray_width);
void RayRender(CRay *I,int width,int height,unsigned int *image,
               float front,float back,double timing,float angle);
void RayRenderPOV(CRay *I,int width,int height,char **headerVLA,
                  char **charVLA,float front,float back,float fov,float angle);
void RayRenderTest(CRay *I,int width,int height,float front,float back,float fov);
void RaySetTTT(CRay *I,int flag,float *ttt);
void RaySetContext(CRay *I,int context);
void RayApplyContexToNormal(CRay *I,float *v);
void RayApplyContextToVertex(CRay *I,float *v);
void RayRenderColorTable(CRay *I,int width,int height,int *image);
int RayTraceThread(CRayThreadInfo *T);


int RayHashThread(CRayHashThreadInfo *T);

#endif









