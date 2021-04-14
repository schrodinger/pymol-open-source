
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

#include <memory>
#include <vector>

#include"Base.h"
#include"Basis.h"
#include"PyMOLGlobals.h"
#include"Image.h"
#include"RenderContext.h"

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
void RayRenderCOLLADA(CRay * I, int width, int height,
                    char **vla_ptr, float front, float back, float fov);
void RayRenderObjMtl(CRay * I, int width, int height, char **objVLA_ptr,
                     char **mtlVLA_ptr, float front, float back, float fov,
                     float angle, float z_corr);
void RayRenderTest(CRay * I, int width, int height, float front, float back, float fov);
void RaySetTTT(CRay * I, int flag, float *ttt);
void RayGetTTT(CRay * I, float *ttt);
void RayPushTTT(CRay * I);
void RayPopTTT(CRay * I);
void RaySetContext(CRay * I, pymol::RenderContext context);
void RayRenderColorTable(CRay * I, int width, int height, int *image);
int RayTraceThread(CRayThreadInfo * T);
int RayGetNPrimitives(CRay * I);
void RayGetScaledAxes(CRay * I, float *xn, float *yn);

int RayHashThread(CRayHashThreadInfo * T);
int RayAntiThread(CRayAntiThreadInfo * T);

// JMS: added so they can be used in COLLADA.c
int RayExpandPrimitives(CRay * I);
int RayTransformFirst(CRay * I, int perspective, int identity);
void RayComputeBox(CRay * I);
int TriangleReverse(CPrimitive * p);


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

namespace cgo{
namespace draw{
  struct cylinder;
  struct custom_cylinder;
  struct custom_cylinder_alpha;
}
};

struct _CRay {

  // methods
  int sphere3fv(const float *v, float r);
  int cylinder3fv(const cgo::draw::cylinder &cyl);
  int cylinder3fv(const cgo::draw::cylinder &cyl, const float alpha1, const float alpha2);
  int customCylinder3fv(const cgo::draw::custom_cylinder &cyl, const float alpha1, const float alpha2);
  int customCylinder3fv(const cgo::draw::custom_cylinder &cyl);
  int customCylinder3fv(const float *v1, const float *v2, float r, const float *c1,
                        const float *c2, const cCylCap cap1, const cCylCap cap2,
                        const float alpha1, const float alpha2);
  int customCylinder3fv(const float *v1, const float *v2, float r, const float *c1,
                        const float *c2, const cCylCap cap1, const cCylCap cap2);
  int customCylinderAlpha3fv(const cgo::draw::custom_cylinder_alpha &cyl);
  int cone3fv(const float *v1, const float *v2, float r1, float r2, const float *c1,
		   const float *c2, cCylCap cap1, cCylCap cap2);
  int sausage3fv(const float *v1, const float *v2, float r, const float *c1, const float *c2);
  void color3fv(const float *c);
  int triangle3fv(
		       const float *v1, const float *v2, const float *v3,
		       const float *n1, const float *n2, const float *n3, const float *c1, const float *c2, const float *c3);
  int triangleTrans3fv(
			    const float *v1, const float *v2, const float *v3,
			    const float *n1, const float *n2, const float *n3,
			    const float *c1, const float *c2, const float *c3,
			    float t1, float t2, float t3);
  void wobble(int mode, const float *par);
  void transparentf(float t);
  int character(int char_id);
  void interiorColor3fv(const float *v, int passive);
  int ellipsoid3fv(const float *v, float r, const float *n1, const float *n2, const float *n3);
  int setLastToNoLighting(char no_lighting);

  /* everything below should be private */
  PyMOLGlobals *G;
  CPrimitive *Primitive;
  int NPrimitive;
  CBasis *Basis;
  int NBasis;
  int *Vert2Prim;
  float CurColor[3], IntColor[3];
  float ModelView[16];
  float ProMatrix[16];
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
  pymol::RenderContext context;
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
  std::shared_ptr<pymol::Image> bkgrd_data;

private:
  int cylinder3fv(const float *v1, const float *v2, float r, const float *c1, const float *c2,
                  const float alpha1, const float alpha2);
};

void RayGetScreenVertex(CRay * I, float *v, float *dest);
float RayGetScreenVertexScale(CRay * I, float *v1);
void RayAdjustZtoScreenZ(CRay * ray, float *pos, float z);
void RayAdjustZtoScreenZofPoint(CRay * ray, float *pos, float *zpoint);
void RaySetPointToWorldScreenRelative(CRay * ray, float *pos, float *screenPt);
float RayGetScaledAllAxesAtPoint(CRay * I, float *pt, float *xn, float *yn, float *zn);
float* RayGetProMatrix(CRay * I);

#endif
