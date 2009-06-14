
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

#include"Base.h"
#include"Basis.h"
#include"PyMOLGlobals.h"

#define cRayMaxBasis 10

typedef struct _CRayAntiThreadInfo CRayAntiThreadInfo;
typedef struct _CRayHashThreadInfo CRayHashThreadInfo;
typedef struct _CRayThreadInfo CRayThreadInfo;

CRay *RayNew(PyMOLGlobals * G, int antialias);
void RayFree(CRay * I);
void RayPrepare(CRay * I, float v0, float v1, float v2,
                float v3, float v4, float v5,
                float fov, float *pos,
                float *mat, float *rotMat,
                float aspRat, int width, int height,
                float pixel_scale, int ortho, float pixel_ratio,
                float back_ratio, float magnified);
void RayRender(CRay * I, unsigned int *image,
               double timing, float angle, int antialias, unsigned int *return_bg);
void RayRenderPOV(CRay * I, int width, int height, char **headerVLA,
                  char **charVLA, float front, float back, float fov, float angle,
                  int antialias);

void RayRenderIDTF(CRay * I, char **node_vla, char **rsrc_vla);

void RayRenderVRML1(CRay * I, int width, int height,
                    char **vla_ptr, float front, float back,
                    float fov, float angle, float z_corr);
void RayRenderVRML2(CRay * I, int width, int height,
                    char **vla_ptr, float front, float back,
                    float fov, float angle, float z_corr);
void RayRenderObjMtl(CRay * I, int width, int height, char **objVLA_ptr,
                     char **mtlVLA_ptr, float front, float back, float fov,
                     float angle, float z_corr);
void RayRenderTest(CRay * I, int width, int height, float front, float back, float fov);
void RaySetTTT(CRay * I, int flag, float *ttt);
void RayGetTTT(CRay * I, float *ttt);
void RayPushTTT(CRay * I);
void RayPopTTT(CRay * I);
void RaySetContext(CRay * I, int context);
void RayApplyContexToNormal(CRay * I, float *v);
void RayApplyContextToVertex(CRay * I, float *v);
void RayRenderColorTable(CRay * I, int width, int height, int *image);
int RayTraceThread(CRayThreadInfo * T);
int RayGetNPrimitives(CRay * I);
void RayGetScaledAxes(CRay * I, float *xn, float *yn);

int RayHashThread(CRayHashThreadInfo * T);
int RayAntiThread(CRayAntiThreadInfo * T);

typedef struct {
  int op;
  int x1, y1, z1;
  int x2, y2, z2;
  int x3, y3, z3;
  int c;
  int r;
} G3dPrimitive;

G3dPrimitive *RayRenderG3d(CRay * I, int width, int height, float front,
                           float back, float fov, int quiet);

struct _CRay {
  void (*fSphere3fv) (CRay * ray, float *v, float r);
  void (*fCylinder3fv) (CRay * ray, float *v1, float *v2, float r, float *c1, float *c2);
  void (*fCustomCylinder3fv) (CRay * ray, float *v1, float *v2, float r, float *c1,
                              float *c2, int cap1, int cap2);
  void (*fCone3fv) (CRay * ray, float *v1, float *v2, float r1, float r2, float *c1,
                    float *c2, int cap1, int cap2);
  void (*fSausage3fv) (CRay * ray, float *v1, float *v2, float r, float *c1, float *c2);
  void (*fColor3fv) (CRay * ray, float *c);
  void (*fTriangle3fv) (CRay * ray,
                        float *v1, float *v2, float *v3,
                        float *n1, float *n2, float *n3, float *c1, float *c2, float *c3);
  void (*fTriangleTrans3fv) (CRay * ray,
                             float *v1, float *v2, float *v3,
                             float *n1, float *n2, float *n3,
                             float *c1, float *c2, float *c3,
                             float t1, float t2, float t3);
  void (*fWobble) (CRay * ray, int mode, float *par);
  void (*fTransparentf) (CRay * ray, float t);
  void (*fCharacter) (CRay * ray, int char_id);
  void (*fInteriorColor3fv) (CRay * ray, float *v, int passive);
  void (*fEllipsoid3fv) (CRay * ray, float *v, float r, float *n1, float *n2, float *n3);
  /* everything below should be private */
  PyMOLGlobals *G;
  CPrimitive *Primitive;
  int NPrimitive;
  CBasis *Basis;
  int NBasis;
  int *Vert2Prim;
  float CurColor[3], IntColor[3];
  float ModelView[16];
  float Rotation[16];
  float Volume[6];
  float Range[3];
  int BigEndian;
  int Wobble;
  float WobbleParam[3];
  float Trans;
  float Random[256];
  int TTTFlag;
  float TTT[16];
  float *TTTStackVLA;
  int TTTStackDepth;
  int Context;
  int CheckInterior;
  float AspRatio;
  int Width, Height;
  float PixelRadius;
  int Ortho;
  float min_box[3];
  float max_box[3];
  int Sampling;
  float PixelRatio;
  float Magnified;              /* ray pixels to screen pixels */
  float FrontBackRatio;
  double PrimSize;
  int PrimSizeCnt;
  float Fov, Pos[3];
};

#endif
