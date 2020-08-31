
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
#define cPrimCharacter 5
#define cPrimEllipsoid 6
#define cPrimCone 7


/* proposed 

#define cPrimNodePush 8
#define cPrimNodePop 9
#define cPrimNodeAdd 10
#define cPrimBreak 11

*/

enum class cCylCap {
  None = 0,
  Flat = 1,
  Round = 2,
};

// legacy names
constexpr auto cCylCapNone = cCylCap::None;
constexpr auto cCylCapFlat = cCylCap::Flat;
constexpr auto cCylCapRound = cCylCap::Round;

#define cCylShaderCap1Flat 0x01
#define cCylShaderCap2Flat 0x02
#define cCylShaderCap1RoundBit 0x04
#define cCylShaderCap2RoundBit 0x08
#define cCylShaderCap1Round (cCylShaderCap1Flat | cCylShaderCap1RoundBit)
#define cCylShaderCap2Round (cCylShaderCap2Flat | cCylShaderCap2RoundBit)

#define cCylShaderInterpColor 0x10
#define cCylShaderBothCapsRound (cCylShaderCap1Round | cCylShaderCap2Round)
#define cCylShaderBothCapsFlat (cCylShaderCap1Flat | cCylShaderCap2Flat)

#define cCylShaderMask 0x1F

typedef struct {
  int vert;
  float v1[3], v2[3], v3[3];
  float n0[3], n1[3], n2[3], n3[3];
  float c1[3], c2[3], c3[3], ic[3], tr[3];      /* ic = interior color, tr = transparency */
  float r1, r2, l1;
  float trans;
  int char_id;
  char type;
  cCylCap cap1, cap2;
  int cull;
  char wobble, ramped, no_lighting;
  /* float wobble_param[3] eliminated to save space */
} CPrimitive;                   /* currently 172 bytes -> appoximately 6.5 million primitives per gigabyte */

typedef struct {
  PyMOLGlobals *G;
  MapType *Map;
  float *Vertex, *Normal, *Precomp;
  float *Radius, *Radius2, MaxRadius, MinVoxel;
  int *Vert2Normal;
  int NVertex;
  int NNormal;
  float LightNormal[3];         /* for lights - this is the direction of the light rays */
  float SpecNormal[3];          /* for computing specular reflections */
  float Color[3];               /* for lights */
  Matrix33f Matrix;
} CBasis;

typedef struct {
  float base[3];                /* where is this light ray starting from */
  CPrimitive *prim;
  float impact[3];
  float tri1, tri2;
  float sphere[3];              /* sphere location if reflecting off of one */
  float surfnormal[3];          /* normal of reflecting surface */
  float dist;
  float dotgle, flat_dotgle;
  float reflect[3];
  float trans;
  float dir[3];                 /* what is the direction of this light ray? */
  float skip[3];
} RayInfo;

typedef struct {
  CBasis *Basis;
  RayInfo *rr;
  int except1, except2;         /* primitives to avoid */
  int *vert2prim;
  int shadow;
  float front;
  float back;
  float excl_trans;
  int trans_shadows;
  int nearest_shadow;
  int check_interior;
  int label_shadow_mode;
  CPrimitive *prim;
  MapCache cache;
  float fudge0, fudge1;
  /* returns */
  int interior_flag;
  int pass;
  float back_dist;
} BasisCallRec;

int BasisInit(PyMOLGlobals * G, CBasis * I, int group_id);
void BasisFinish(CBasis * I, int group_id);
int BasisMakeMap(CBasis * I, int *vert2prim, CPrimitive * prim, int n_prim,
		 float *volume,
		 int group_id, int block_base,
		 int perspective, float front, float size_hint);

void BasisSetupMatrix(CBasis * I);
void BasisGetTriangleNormal(CBasis * I, RayInfo * r, int i, float *fc, int perspective);
void BasisGetEllipsoidNormal(CBasis * I, RayInfo * r, int i, int perspective);
void BasisTrianglePrecompute(float *v1, float *v2, float *v3, float *pre);
void BasisTrianglePrecomputePerspective(float *v1, float *v2, float *v3, float *pre);

int BasisHitPerspective(BasisCallRec * BC);
int BasisHitOrthoscopic(BasisCallRec * BC);
int BasisHitShadow(BasisCallRec * BC);

void BasisGetTriangleFlatDotgle(CBasis * I, RayInfo * r, int i);
void BasisGetTriangleFlatDotglePerspective(CBasis * I, RayInfo * r, int i);

void BasisCylinderSausagePrecompute(float *dir, float *pre);

#define PROFILE_BASIS_OFF

#endif
