

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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"

#include"CarveHelper.h"
#include"Isosurf.h"
#include"Tetsurf.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Crystal.h"
#include"Vector.h"
#include"Feedback.h"
#include"P.h"

#define Trace_OFF

#define O3(field,P1,P2,P3,offs) ((field)->get<float>(P1+offs[0],P2+offs[1],P3+offs[2]))

#define O4Ptr(field,P1,P2,P3,P4,offs) ((field)->ptr<float>(P1+offs[0],P2+offs[1],P3+offs[2],P4))

#define I3(field,P1,P2,P3) ((field)->get<int>(P1,P2,P3))

typedef struct {
  float Point[3];
  float Normal[3];
  int NormalFlag;
  int Link;
} PointType;

typedef struct {
  PointType *p[3];
  float n[3];
  int done;
} TriangleType;

typedef struct {
  int link;
  int tri;
} PointLinkType;

#define EdgePtPtr(field,P2,P3,P4,P5) ((field)->ptr<PointType>(P2,P3,P4,P5))

struct _CTetsurf {
  PyMOLGlobals *G;
  TriangleType *Tri;
  PointLinkType *PtLink;

  CField *VertexCodes;
  CField *ActiveEdges;
  CField *Point;

  int AbsDim[3], CurDim[3], CurOff[3];
  int Max[3];
  CField *Coord, *Data, *Grad;
  float Level;
  int Edge[6020];               /* 6017 */
  int EdgeStart[256];
  int TotPrim;
};

static int TetsurfAlloc(CTetsurf * II);
static void TetsurfPurge(CTetsurf * II);
static int TetsurfCodeVertices(CTetsurf * II);

static int TetsurfFindActiveBoxes(CTetsurf * II, cIsosurfaceMode, int& n_strip, int n_vert,
                                  pymol::vla<int>& strip_l, pymol::vla<float>& vert,
                                  const CarveHelper*, cIsosurfaceSide);

#define TetsurfSubSize		50

static void copy3fn(float *v1, float *v2)
{
  v2[0] = v1[0];
  v2[1] = v1[1];
  v2[2] = v1[2];
}


/*===========================================================================*/
static int ProcessTetrahedron(int *edge, int nv, int v0, int v1, int v2, int v3,
                              int e01, int e02, int e03, int e12, int e13, int e23,
                              int reflect)
{
  int bits = (v3 << 3) + (v2 << 2) + (v1 << 1) + v0;
  if(reflect)
    bits = 0xF - bits;
  switch (bits) {
  case 0x0:
    break;
  case 0x1:
    edge[nv++] = e01;
    edge[nv++] = e02;
    edge[nv++] = e03;
    break;
  case 0x2:
    edge[nv++] = e01;
    edge[nv++] = e13;
    edge[nv++] = e12;
    break;
  case 0x3:
    edge[nv++] = e13;
    edge[nv++] = e12;
    edge[nv++] = e02;
    edge[nv++] = e03;
    edge[nv++] = e13;
    edge[nv++] = e02;
    break;
  case 0x4:
    edge[nv++] = e12;
    edge[nv++] = e23;
    edge[nv++] = e02;
    break;
  case 0x5:
    edge[nv++] = e01;
    edge[nv++] = e12;
    edge[nv++] = e03;
    edge[nv++] = e12;
    edge[nv++] = e23;
    edge[nv++] = e03;
    break;
  case 0x6:
    edge[nv++] = e01;
    edge[nv++] = e13;
    edge[nv++] = e02;
    edge[nv++] = e13;
    edge[nv++] = e23;
    edge[nv++] = e02;
    break;
  case 0x7:
    edge[nv++] = e03;
    edge[nv++] = e13;
    edge[nv++] = e23;
    break;
  case 0x8:
    edge[nv++] = e03;
    edge[nv++] = e23;
    edge[nv++] = e13;
    break;
  case 0x9:
    edge[nv++] = e13;
    edge[nv++] = e01;
    edge[nv++] = e02;
    edge[nv++] = e02;
    edge[nv++] = e23;
    edge[nv++] = e13;
    break;
  case 0xA:
    edge[nv++] = e01;
    edge[nv++] = e03;
    edge[nv++] = e12;
    edge[nv++] = e03;
    edge[nv++] = e23;
    edge[nv++] = e12;
    break;
  case 0xB:
    edge[nv++] = e23;
    edge[nv++] = e12;
    edge[nv++] = e02;
    break;
  case 0xC:
    edge[nv++] = e13;
    edge[nv++] = e02;
    edge[nv++] = e12;
    edge[nv++] = e03;
    edge[nv++] = e02;
    edge[nv++] = e13;
    break;
  case 0xD:
    edge[nv++] = e01;
    edge[nv++] = e12;
    edge[nv++] = e13;
    break;
  case 0xE:
    edge[nv++] = e01;
    edge[nv++] = e03;
    edge[nv++] = e02;
    break;
  case 0xF:
    break;
  }
  return nv;
}


/*===========================================================================*/
static CTetsurf *TetsurfNew(PyMOLGlobals * G)
{

  /* there are six tetrahedrons in each cube, and there are 
     sixteen different types of tetrahedrons */

  /* one will encounter 2^8 different kinds of cubes. */

  /* bits -> tetrahedral vertices 
     zyx
     0: 000
     1: 001
     2: 010
     3: 011
     4: 100
     5: 101
     6: 110
     7: 111
   */

  /* edges of the tetrahedron 

     vertex
     3210
     0:0011
     1:0101
     2:0110
     3:1001
     4:1010
     5:1100
   */

  /* edges of the cube 

     76543210
     0:00000011 edge 000-001
     1:00000101 edge 000-010
     2:00001001 edge 000-011
     3:00010001 edge 000-100
     4:00100001 edge 000-101
     5:01000001 edge 000-110
     6:10000001 edge 000-111
     7:00001010 edge 001-011
     8:00100010 edge 001-101
     9:10000010 edge 001-111
     0xA:00001100 edge 010-011
     0xB:01000100 edge 010-110
     0xC:10000100 edge 010-111
     0xD:00110000 edge 100-101
     0xE:01010000 edge 100-110
     0xF:10010000 edge 100-111
     0x10:10001000 edge 011-111
     0x11:10100000 edge 101-111
     0x12:11000000 edge 110-111
   */

#define cE_000_001 0x00
#define cE_000_010 0x01
#define cE_000_011 0x02
#define cE_000_100 0x03
#define cE_000_101 0x04
#define cE_000_110 0x05
#define cE_000_111 0x06
#define cE_001_011 0x07
#define cE_001_101 0x08
#define cE_001_111 0x09
#define cE_010_011 0x0A
#define cE_010_110 0x0B
#define cE_010_111 0x0C
#define cE_100_101 0x0D
#define cE_100_110 0x0E
#define cE_100_111 0x0F
#define cE_011_111 0x10
#define cE_101_111 0x11
#define cE_110_111 0x12

#define cM_000_001 0x00001
#define cM_000_010 0x00002
#define cM_000_011 0x00004
#define cM_000_100 0x00008
#define cM_000_101 0x00010
#define cM_000_110 0x00020
#define cM_000_111 0x00040
#define cM_001_011 0x00080
#define cM_001_101 0x00100
#define cM_001_111 0x00200
#define cM_010_011 0x00400
#define cM_010_110 0x00800
#define cM_010_111 0x01000
#define cM_100_101 0x02000
#define cM_100_110 0x04000
#define cM_100_111 0x08000
#define cM_011_111 0x10000
#define cM_101_111 0x20000
#define cM_110_111 0x40000

  CTetsurf *I = pymol::calloc<CTetsurf>(1);
  int c;
  int nv = 1;
  int last_nv;
  int v000, v100, v010, v110, v001, v101, v011, v111;

  I->G = G;
  I->Tri = NULL;
  I->PtLink = NULL;

  I->VertexCodes = NULL;
  I->ActiveEdges = NULL;
  I->Point = NULL;

  last_nv = nv;

  for(c = 0; c < 256; c++) {

    v000 = (c & 0x01) ? 1 : 0;
    v001 = (c & 0x02) ? 1 : 0;
    v010 = (c & 0x04) ? 1 : 0;
    v011 = (c & 0x08) ? 1 : 0;
    v100 = (c & 0x10) ? 1 : 0;
    v101 = (c & 0x20) ? 1 : 0;
    v110 = (c & 0x40) ? 1 : 0;
    v111 = (c & 0x80) ? 1 : 0;

    /* tetrahedron 0: 000, 001, 011, 111 */
    nv = ProcessTetrahedron(I->Edge, nv, v000, v001, v011, v111,
                            cE_000_001,
                            cE_000_011,
                            cE_000_111, cE_001_011, cE_001_111, cE_011_111, 0);
    /* tetrahedron 1: 000, 001, 101, 111 */
    nv = ProcessTetrahedron(I->Edge, nv, v000, v001, v101, v111,
                            cE_000_001,
                            cE_000_101,
                            cE_000_111, cE_001_101, cE_001_111, cE_101_111, 1);

    /* tetrahedron 2: 000, 010, 011, 111 */
    nv = ProcessTetrahedron(I->Edge, nv, v000, v010, v011, v111,
                            cE_000_010,
                            cE_000_011,
                            cE_000_111, cE_010_011, cE_010_111, cE_011_111, 1);
    /* tetrahedron 3: 000, 010, 110, 111 */
    nv = ProcessTetrahedron(I->Edge, nv, v000, v010, v110, v111,
                            cE_000_010,
                            cE_000_110,
                            cE_000_111, cE_010_110, cE_010_111, cE_110_111, 0);

    /* tetrahedron 4: 000, 100, 101, 111 */
    nv = ProcessTetrahedron(I->Edge, nv, v000, v100, v101, v111,
                            cE_000_100,
                            cE_000_101,
                            cE_000_111, cE_100_101, cE_100_111, cE_101_111, 0);

    /* tetrahedron 5: 000, 100, 110, 111 */
    nv = ProcessTetrahedron(I->Edge, nv, v000, v100, v110, v111,
                            cE_000_100,
                            cE_000_110,
                            cE_000_111, cE_100_110, cE_100_111, cE_110_111, 1);

    I->Edge[nv++] = (-1);       /* sentinel */
    I->EdgeStart[c] = last_nv;
    last_nv = nv;
  }
  /*   printf("%d\n",nv);
     for(c=1;c<nv;c++) {
     if(Edge[c]<0) {
     printf("\n");
     } else {
     printf("%02X ",Edge[c]);
     }
     }
   */
  return I;
}

int TetsurfInit(PyMOLGlobals * G)
{
  G->Tetsurf = TetsurfNew(G);
  return 1;
}


/*===========================================================================*/
static void _TetsurfFree(CTetsurf * I)
{
  FreeP(I);
}

void TetsurfFree(PyMOLGlobals * G)
{
  _TetsurfFree(G->Tetsurf);
  G->Tetsurf = NULL;
}


/*===========================================================================*/
/**
 * @param mn Minimum point in real space (3f)
 * @param mx Maximum point in real space (3f)
 * @param[out] range Minimum and maximum field indices (6i)
 */
void TetsurfGetRange(PyMOLGlobals * G,
    const Isofield* field,
    const CCrystal* cryst, const float* mn, const float* mx, int* range)
{
  float rmn[3], rmx[3];
  float imn[3], imx[3];
  float mix[24], imix[24];
  int a, b;
  PRINTFD(G, FB_Isosurface)
    " IsosurfGetRange: entered mn: %4.2f %4.2f %4.2f mx: %4.2f %4.2f %4.2f\n",
    mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]
    ENDFD;

  for(a = 0; a < 3; a++) {
    rmn[a] = F4(field->points, 0, 0, 0, a);
    rmx[a] = F4(field->points,
                field->dimensions[0] - 1,
                field->dimensions[1] - 1, field->dimensions[2] - 1, a);
  }

  /* get min/max extents of map in fractional space */

  transform33f3f(cryst->realToFrac(), rmn, imn);
  transform33f3f(cryst->realToFrac(), rmx, imx);

  mix[0] = mn[0];
  mix[1] = mn[1];
  mix[2] = mn[2];

  mix[3] = mx[0];
  mix[4] = mn[1];
  mix[5] = mn[2];

  mix[6] = mn[0];
  mix[7] = mx[1];
  mix[8] = mn[2];

  mix[9] = mn[0];
  mix[10] = mn[1];
  mix[11] = mx[2];

  mix[12] = mx[0];
  mix[13] = mx[1];
  mix[14] = mn[2];

  mix[15] = mx[0];
  mix[16] = mn[1];
  mix[17] = mx[2];

  mix[18] = mn[0];
  mix[19] = mx[1];
  mix[20] = mx[2];

  mix[21] = mx[0];
  mix[22] = mx[1];
  mix[23] = mx[2];

  /* compute min/max of query in fractional space */

  for(b = 0; b < 8; b++) {
    transform33f3f(cryst->realToFrac(), mix + 3 * b, imix + 3 * b);
  }

  for(a = 0; a < 3; a++) {
    if(imx[a] != imn[a]) {      /* protect against div by zero */
      int b;
      int mini = 0, maxi = 0, tst_min, tst_max;
      float cur;
      for(b = 0; b < 8; b++) {
        cur =
          ((field->dimensions[a] - 1) * (imix[a + 3 * b] - imn[a]) / (imx[a] - imn[a]));
        tst_min = (int) floor(cur);
        tst_max = ((int) ceil(cur)) + 1;

        if(!b) {
          mini = tst_min;
          maxi = tst_max;
        } else {
          if(mini > tst_min)
            mini = tst_min;
          if(maxi <= tst_max)
            maxi = tst_max;
        }
      }

      range[a] = mini;

      range[a + 3] = maxi;
    } else {
      range[a] = 0;
      range[a + 3] = 1;
    }
    if(range[a] < 0)
      range[a] = 0;
    if(range[a] > field->dimensions[a])
      range[a] = field->dimensions[a];
    if(range[a + 3] < 0)
      range[a + 3] = 0;
    if(range[a + 3] > field->dimensions[a])
      range[a + 3] = field->dimensions[a];
  }
  PRINTFD(G, FB_Isosurface)
    " IsosurfGetRange: returning range: %d %d %d %d %d %d\n",
    range[0], range[1], range[2], range[3], range[4], range[5]
    ENDFD;
}


/*===========================================================================*/
/**
 * Compute an isosurface using the "marching tetrahedra" algorithm.
 *
 * @param[in,out] field     Map data (field->gradients will be generated if
 *                          missing for cIsosurfaceMode::triangles_grad_normals)
 * @param[in] level         Contour level
 * @param[out] num          Number of vertices + normals per strip (e.g. 6 for
 *                          triangle strip with a single triangle)
 * @param[out] vert         Normals (triangle modes only) and vertices for
 *                          strips according to `num`
 * @param[in] range         Min and max indices of box (6i - can be NULL,
 *                          then use entire field)
 * @param[in] mode          Geometry mode (points, lines, triangles)
 * @param[in] carvehelper   For carving (optional)
 * @param[in] side          Front or back face (triangle winding order and
 *                          normal)
 * @return Number of primitives
 */
int TetsurfVolume(PyMOLGlobals* G, Isofield* field, float level,
    pymol::vla<int>& num,    //
    pymol::vla<float>& vert, //
    const int* range,        //
    cIsosurfaceMode mode,    //
    const CarveHelper* carvehelper, //
    cIsosurfaceSide side)
{

  CTetsurf *I;
  if(PIsGlutThread()) {
    I = G->Tetsurf;
  } else {
    I = TetsurfNew(G);
  }

  {
    int ok = true;
    int Steps[3];
    int c, i, j, k;
    int range_store[6];
    int n_strip = 0;
    int n_vert = 0;
    int tot_prim = 0;

    if(mode == cIsosurfaceMode::triangles_grad_normals)
      IsofieldComputeGradients(G, field);

    I->TotPrim = 0;
    if(range) {
      for(c = 0; c < 3; c++) {
        I->AbsDim[c] = field->dimensions[c];
        I->CurDim[c] = TetsurfSubSize + 1;
        Steps[c] = 1 + ((range[3 + c] - range[c]) - 1) / TetsurfSubSize;
      }
    } else {
      range = range_store;
      for(c = 0; c < 3; c++) {
        range_store[c] = 0;
        range_store[3 + c] = field->dimensions[c];
        I->AbsDim[c] = field->dimensions[c];
        I->CurDim[c] = TetsurfSubSize + 1;
        Steps[c] = 1 + (I->AbsDim[c] - 1) / TetsurfSubSize;
      }
    }

    /*   for(c=0;c<3;c++) {
       printf("range %d %d %d\n",c,range[c],range[c+3]);
       printf("steps %d\n",Steps[c]);
       }
     */

    I->Coord = field->points.get();
    I->Grad = field->gradients.get();
    I->Data = field->data.get();
    I->Level = level;
    if(ok)
      ok = TetsurfAlloc(I);

    if(ok) {

      for(i = 0; i < Steps[0]; i++)
        for(j = 0; j < Steps[1]; j++)
          for(k = 0; k < Steps[2]; k++) {
            I->CurOff[0] = TetsurfSubSize * i;
            I->CurOff[1] = TetsurfSubSize * j;
            I->CurOff[2] = TetsurfSubSize * k;
            for(c = 0; c < 3; c++)
              I->CurOff[c] += range[c];
            for(c = 0; c < 3; c++) {
              I->Max[c] = (range[3 + c] - I->CurOff[c]);
              if(I->Max[c] > (TetsurfSubSize + 1))
                I->Max[c] = (TetsurfSubSize + 1);
            }
            /*         
               for(c=0;c<3;c++)
               printf(" TetsurfVolume: c: %i I->CurOff[c]: %i I->Max[c] %i\n",c,I->CurOff[c],I->Max[c]); 
             */

            if(ok) {
              if(TetsurfCodeVertices(I))
                n_vert = TetsurfFindActiveBoxes(I, mode, n_strip, n_vert, num, vert,
                                                carvehelper, side);
            }
          }
      TetsurfPurge(I);
    }

    if(Feedback(G, FB_Isosurface, FB_Blather)) {
      if(static_cast<int>(mode) < 2) {
        printf(" TetsurfVolume: Surface generated using %d vertices.\n", n_vert);
      } else {
        printf(" TetsurfVolume: Surface generated using %d triangles.\n", I->TotPrim);
      }
    }

    /* sentinel strip (0 length) */

    num.resize(++n_strip);
    num[n_strip - 1] = 0;

    /* shrinks sizes for more efficient RAM usage */

    vert.resize(n_vert * 3);

    tot_prim = I->TotPrim;

    if(!PIsGlutThread()) {
      _TetsurfFree(I);
    }
    return (tot_prim);
  }
}


/*===========================================================================*/
static int TetsurfAlloc(CTetsurf * II)
{
  CTetsurf *I = II;
  PyMOLGlobals *G = I->G;

  int ok = true;
  int dim4[4];
  int a;
  for(a = 0; a < 3; a++)
    dim4[a] = I->CurDim[a];
  dim4[3] = 3;

  I->VertexCodes = CField::make<int>(G, I->CurDim, 3);
  ErrChkPtr(G, I->VertexCodes);
  I->ActiveEdges = CField::make<int>(G, I->CurDim, 3);
  ErrChkPtr(G, I->ActiveEdges);

  dim4[3] = 7;                  /* seven different ways now... */
  I->Point = CField::make<PointType>(G, dim4, 4);
  ErrChkPtr(G, I->Point);

  I->Tri = VLAlloc(TriangleType, 50000);
  I->PtLink = VLAlloc(PointLinkType, 50000);

  if(!(I->VertexCodes && I->ActiveEdges && I->Point)) {
    TetsurfPurge(I);
    ok = false;
  }
#ifdef Trace
  printf(" TetsurfAlloc: ok: %i\n", ok);
#endif
  return (ok);
}


/*===========================================================================*/
static void TetsurfPurge(CTetsurf * II)
{
  CTetsurf *I = II;
  if(I->Tri) {
    VLAFreeP(I->Tri);
  }
  if(I->PtLink) {
    VLAFreeP(I->PtLink);
  }
  DeleteP(I->VertexCodes);
  DeleteP(I->ActiveEdges);
  DeleteP(I->Point);
}


/*===========================================================================*/
inline
static void TetsurfInterpolate2(float *pt, float *v0, float l0, float *v1, float l1,
                                float level)
{
  float ratio;
  ratio = (level - l0) / (l1 - l0);
  pt[0] = v0[0] + (v1[0] - v0[0]) * ratio;
  pt[1] = v0[1] + (v1[1] - v0[1]) * ratio;
  pt[2] = v0[2] + (v1[2] - v0[2]) * ratio;
}


/*===========================================================================*/
static void TetsurfInterpolate4(float *pt, float *v0, float l0, float *v1, float l1,
                                float l2, float l3, float level)
{
  float ratio;
  float v[3], l;
  average3f(v0, v1, v);
  l = (l0 + l1 + l2 + l3) * 0.25F;
  if(((l > level) && (l1 > level)) || ((l <= level) && (l0 > level))) { /* l0 vs l */
    ratio = (level - l0) / (l - l0);
    pt[0] = v0[0] + (v[0] - v0[0]) * ratio;
    pt[1] = v0[1] + (v[1] - v0[1]) * ratio;
    pt[2] = v0[2] + (v[2] - v0[2]) * ratio;
  } else {
    ratio = (level - l1) / (l - l1);
    pt[0] = v1[0] + (v[0] - v1[0]) * ratio;
    pt[1] = v1[1] + (v[1] - v1[1]) * ratio;
    pt[2] = v1[2] + (v[2] - v1[2]) * ratio;
  }
}


/*===========================================================================*/
static void TetsurfInterpolate8(float *pt, float *v0, float l0, float *v1, float l1,
                                float l2, float l3, float l4,
                                float l5, float l6, float l7, float level)
{
  float ratio;
  float v[3], l;
  average3f(v0, v1, v);
  l = (l0 + l1 + l2 + l3 + l4 + l5 + l6 + l7) * 0.125F;
  if(((l > level) && (l1 > level)) || ((l <= level) && (l0 > level))) { /* l0 vs l */
    ratio = (level - l0) / (l - l0);
    pt[0] = v0[0] + (v[0] - v0[0]) * ratio;
    pt[1] = v0[1] + (v[1] - v0[1]) * ratio;
    pt[2] = v0[2] + (v[2] - v0[2]) * ratio;
  } else {
    ratio = (level - l1) / (l - l1);
    pt[0] = v1[0] + (v[0] - v1[0]) * ratio;
    pt[1] = v1[1] + (v[1] - v1[1]) * ratio;
    pt[2] = v1[2] + (v[2] - v1[2]) * ratio;
  }
}


/*===========================================================================*/
static int TetsurfFindActiveBoxes(CTetsurf * II, cIsosurfaceMode mode, int &n_strip, int n_vert,
                                  pymol::vla<int>& strip_l, pymol::vla<float>& vert_,
                                  const CarveHelper* carvehelper, cIsosurfaceSide side)
{
  float** const vert = &vert_;
  CTetsurf *I = II;
  int a, b, i, j, k;
#ifdef Trace
  int ECount = 0;
#endif
  int i000, i001, i010, i011, i100, i101, i110, i111;
  float *c000, *c001, *c010, *c011, *c100, *c101, *c110, *c111;
  float d000, d001, d010, d011, d100, d101, d110, d111;
  float *g000 = NULL, *g001 = NULL, *g010 = NULL, *g011 = NULL, *g100 = NULL, *g101 =
    NULL, *g110 = NULL, *g111 = NULL;

  int active;
  int n_active = 0;
  int n_start = 0;
  PointType *e[19], *p0, *p1, *p2;
  int code;
  int eidx;
  int idx;
  TriangleType *tt;
  int n_tri = 0;
  int n_link = 1;

  FieldZero(I->Point);          /* sets initial links to zero */
  FieldZero(I->ActiveEdges);
  n_start = n_vert;
  for(i = 0; i < (I->Max[0] - 1); i++)
    for(j = 0; j < (I->Max[1] - 1); j++)
      for(k = 0; k < (I->Max[2] - 1); k++) {
        active = 0;

        i000 = I3(I->VertexCodes, i, j, k);
        i001 = I3(I->VertexCodes, i, j, k + 1);
        i010 = I3(I->VertexCodes, i, j + 1, k);
        i011 = I3(I->VertexCodes, i, j + 1, k + 1);
        i100 = I3(I->VertexCodes, i + 1, j, k);
        i101 = I3(I->VertexCodes, i + 1, j, k + 1);
        i110 = I3(I->VertexCodes, i + 1, j + 1, k);
        i111 = I3(I->VertexCodes, i + 1, j + 1, k + 1);

        if((i000 != i001) || (i001 != i010) || (i010 != i011) || (i011 != i100) || (i100 != i101) || (i101 != i110) || (i110 != i111)) {        /* this is an active box */

          c000 = O4Ptr(I->Coord, i, j, k, 0, I->CurOff);
          c001 = O4Ptr(I->Coord, i, j, k + 1, 0, I->CurOff);
          c010 = O4Ptr(I->Coord, i, j + 1, k, 0, I->CurOff);
          c011 = O4Ptr(I->Coord, i, j + 1, k + 1, 0, I->CurOff);
          c100 = O4Ptr(I->Coord, i + 1, j, k, 0, I->CurOff);
          c101 = O4Ptr(I->Coord, i + 1, j, k + 1, 0, I->CurOff);
          c110 = O4Ptr(I->Coord, i + 1, j + 1, k, 0, I->CurOff);
          c111 = O4Ptr(I->Coord, i + 1, j + 1, k + 1, 0, I->CurOff);

          if (mode == cIsosurfaceMode::triangles_grad_normals) {
            g000 = O4Ptr(I->Grad, i, j, k, 0, I->CurOff);
            g001 = O4Ptr(I->Grad, i, j, k + 1, 0, I->CurOff);
            g010 = O4Ptr(I->Grad, i, j + 1, k, 0, I->CurOff);
            g011 = O4Ptr(I->Grad, i, j + 1, k + 1, 0, I->CurOff);
            g100 = O4Ptr(I->Grad, i + 1, j, k, 0, I->CurOff);
            g101 = O4Ptr(I->Grad, i + 1, j, k + 1, 0, I->CurOff);
            g110 = O4Ptr(I->Grad, i + 1, j + 1, k, 0, I->CurOff);
            g111 = O4Ptr(I->Grad, i + 1, j + 1, k + 1, 0, I->CurOff);
          }

          d000 = O3(I->Data, i, j, k, I->CurOff);
          d001 = O3(I->Data, i, j, k + 1, I->CurOff);
          d010 = O3(I->Data, i, j + 1, k, I->CurOff);
          d011 = O3(I->Data, i, j + 1, k + 1, I->CurOff);
          d100 = O3(I->Data, i + 1, j, k, I->CurOff);
          d101 = O3(I->Data, i + 1, j, k + 1, I->CurOff);
          d110 = O3(I->Data, i + 1, j + 1, k, I->CurOff);
          d111 = O3(I->Data, i + 1, j + 1, k + 1, I->CurOff);

          e[cE_000_001] = EdgePtPtr(I->Point, i, j, k, cE_000_001);
          e[cE_000_010] = EdgePtPtr(I->Point, i, j, k, cE_000_010);
          e[cE_000_011] = EdgePtPtr(I->Point, i, j, k, cE_000_011);
          e[cE_000_100] = EdgePtPtr(I->Point, i, j, k, cE_000_100);
          e[cE_000_101] = EdgePtPtr(I->Point, i, j, k, cE_000_101);

          e[cE_000_110] = EdgePtPtr(I->Point, i, j, k, cE_000_110);
          e[cE_000_111] = EdgePtPtr(I->Point, i, j, k, cE_000_111);
          e[cE_001_011] = EdgePtPtr(I->Point, i, j, k + 1, cE_000_010);
          e[cE_001_101] = EdgePtPtr(I->Point, i, j, k + 1, cE_000_100);
          e[cE_001_111] = EdgePtPtr(I->Point, i, j, k + 1, cE_000_110);

          e[cE_010_011] = EdgePtPtr(I->Point, i, j + 1, k, cE_000_001);
          e[cE_010_110] = EdgePtPtr(I->Point, i, j + 1, k, cE_000_100);
          e[cE_010_111] = EdgePtPtr(I->Point, i, j + 1, k, cE_000_101);
          e[cE_100_101] = EdgePtPtr(I->Point, i + 1, j, k, cE_000_001);
          e[cE_100_110] = EdgePtPtr(I->Point, i + 1, j, k, cE_000_010);

          e[cE_100_111] = EdgePtPtr(I->Point, i + 1, j, k, cE_000_011);
          e[cE_011_111] = EdgePtPtr(I->Point, i, j + 1, k + 1, cE_000_100);
          e[cE_101_111] = EdgePtPtr(I->Point, i + 1, j, k + 1, cE_000_010);
          e[cE_110_111] = EdgePtPtr(I->Point, i + 1, j + 1, k, cE_000_001);

          /* Generate interpolated coordinates for all active edges */

          if(i000 != i001) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_001)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_001;
              TetsurfInterpolate2(e[cE_000_001]->Point, c000, d000, c001, d001, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_000_001]->Normal, g000, d000, g001, d001,
                                    I->Level);
            }
            active |= cM_000_001;
          }
          if(i000 != i010) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_010)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_010;
              TetsurfInterpolate2(e[cE_000_010]->Point, c000, d000, c010, d010, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_000_010]->Normal, g000, d000, g010, d010,
                                    I->Level);

            }
            active |= cM_000_010;
          }
          if(i000 != i011) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_011)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_011;
              TetsurfInterpolate4(e[cE_000_011]->Point, c000, d000, c011, d011, d001,
                                  d010, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate4(e[cE_000_011]->Normal, g000, d000, g011, d011, d001,
                                    d010, I->Level);
            }
            active |= cM_000_011;
          }
          if(i000 != i100) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_100)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_100;
              TetsurfInterpolate2(e[cE_000_100]->Point, c000, d000, c100, d100, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_000_100]->Normal, g000, d000, g100, d100,
                                    I->Level);

            }
            active |= cM_000_100;
          }
          if(i000 != i101) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_101)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_101;
              TetsurfInterpolate4(e[cE_000_101]->Point, c000, d000, c101, d101, d100,
                                  d001, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate4(e[cE_000_101]->Normal, g000, d000, g101, d101, d100,
                                    d001, I->Level);
            }
            active |= cM_000_101;
          }
          if(i000 != i110) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_110)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_110;
              TetsurfInterpolate4(e[cE_000_110]->Point, c000, d000, c110, d110, d100,
                                  d010, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate4(e[cE_000_110]->Normal, g000, d000, g110, d110, d100,
                                    d010, I->Level);
            }
            active |= cM_000_110;
          }
          if(i000 != i111) {
            if(!(I3(I->ActiveEdges, i, j, k) & cM_000_111)) {
              I3(I->ActiveEdges, i, j, k) |= cM_000_111;
              TetsurfInterpolate8(e[cE_000_111]->Point,
                                  c000, d000, c111, d111,
                                  d001, d010, d011, d100, d101, d110, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate8(e[cE_000_111]->Normal,
                                    g000, d000, g111, d111,
                                    d001, d010, d011, d100, d101, d110, I->Level);
            }
            active |= cM_000_111;
          }
          if(i001 != i011) {
            if(!(I3(I->ActiveEdges, i, j, k + 1) & cM_000_010)) {
              I3(I->ActiveEdges, i, j, k + 1) |= cM_000_010;
              TetsurfInterpolate2(e[cE_001_011]->Point, c001, d001, c011, d011, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_001_011]->Normal, g001, d001, g011, d011,
                                    I->Level);
            }
            active |= cM_001_011;
          }
          if(i001 != i101) {
            if(!(I3(I->ActiveEdges, i, j, k + 1) & cM_000_100)) {
              I3(I->ActiveEdges, i, j, k + 1) |= cM_000_100;
              TetsurfInterpolate2(e[cE_001_101]->Point, c001, d001, c101, d101, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_001_101]->Normal, g001, d001, g101, d101,
                                    I->Level);
            }
            active |= cM_001_101;
          }
          if(i001 != i111) {
            if(!(I3(I->ActiveEdges, i, j, k + 1) & cM_000_110)) {
              I3(I->ActiveEdges, i, j, k + 1) |= cM_000_110;
              TetsurfInterpolate4(e[cE_001_111]->Point, c001, d001, c111, d111, d101,
                                  d011, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate4(e[cE_001_111]->Normal, g001, d001, g111, d111, d101,
                                    d011, I->Level);
            }
            active |= cM_001_111;
          }
          if(i010 != i011) {
            if(!(I3(I->ActiveEdges, i, j + 1, k) & cM_000_001)) {
              I3(I->ActiveEdges, i, j + 1, k) |= cM_000_001;
              TetsurfInterpolate2(e[cE_010_011]->Point, c010, d010, c011, d011, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_010_011]->Normal, g010, d010, g011, d011,
                                    I->Level);
            }
            active |= cM_010_011;
          }
          if(i010 != i110) {
            if(!(I3(I->ActiveEdges, i, j + 1, k) & cM_000_100)) {
              I3(I->ActiveEdges, i, j + 1, k) |= cM_000_100;
              TetsurfInterpolate2(e[cE_010_110]->Point, c010, d010, c110, d110, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_010_110]->Normal, g010, d010, g110, d110,
                                    I->Level);
            }
            active |= cM_010_110;
          }
          if(i010 != i111) {
            if(!(I3(I->ActiveEdges, i, j + 1, k) & cM_000_101)) {
              I3(I->ActiveEdges, i, j + 1, k) |= cM_000_101;
              TetsurfInterpolate4(e[cE_010_111]->Point, c010, d010, c111, d111, d110,
                                  d011, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate4(e[cE_010_111]->Normal, g010, d010, g111, d111, d110,
                                    d011, I->Level);
            }
            active |= cM_010_111;
          }
          if(i100 != i101) {
            if(!(I3(I->ActiveEdges, i + 1, j, k) & cM_000_001)) {
              I3(I->ActiveEdges, i + 1, j, k) |= cM_000_001;
              TetsurfInterpolate2(e[cE_100_101]->Point, c100, d100, c101, d101, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_100_101]->Normal, g100, d100, g101, d101,
                                    I->Level);
            }
            active |= cM_100_101;
          }
          if(i100 != i110) {
            if(!(I3(I->ActiveEdges, i + 1, j, k) & cM_000_010)) {
              I3(I->ActiveEdges, i + 1, j, k) |= cM_000_010;
              TetsurfInterpolate2(e[cE_100_110]->Point, c100, d100, c110, d110, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_100_110]->Normal, g100, d100, g110, d110,
                                    I->Level);
            }
            active |= cM_100_110;
          }
          if(i100 != i111) {
            if(!(I3(I->ActiveEdges, i + 1, j, k) & cM_000_011)) {
              I3(I->ActiveEdges, i + 1, j, k) |= cM_000_011;
              TetsurfInterpolate4(e[cE_100_111]->Point, c100, d100, c111, d111, d101,
                                  d110, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate4(e[cE_100_111]->Normal, g100, d100, g111, d111, d101,
                                    d110, I->Level);
            }
            active |= cM_100_111;
          }
          if(i011 != i111) {
            if(!(I3(I->ActiveEdges, i, j + 1, k + 1) & cM_000_100)) {
              I3(I->ActiveEdges, i, j + 1, k + 1) |= cM_000_100;
              TetsurfInterpolate2(e[cE_011_111]->Point, c011, d011, c111, d111, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_011_111]->Normal, g011, d011, g111, d111,
                                    I->Level);
            }
            active |= cM_011_111;
          }
          if(i101 != i111) {
            if(!(I3(I->ActiveEdges, i + 1, j, k + 1) & cM_000_010)) {
              I3(I->ActiveEdges, i + 1, j, k + 1) |= cM_000_010;
              TetsurfInterpolate2(e[cE_101_111]->Point, c101, d101, c111, d111, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_101_111]->Normal, g101, d101, g111, d111,
                                    I->Level);
            }
            active |= cM_101_111;
          }
          if(i110 != i111) {
            if(!(I3(I->ActiveEdges, i + 1, j + 1, k) & cM_000_001)) {
              I3(I->ActiveEdges, i + 1, j + 1, k) |= cM_000_001;
              TetsurfInterpolate2(e[cE_110_111]->Point, c110, d110, c111, d111, I->Level);
              if (mode == cIsosurfaceMode::triangles_grad_normals)
                TetsurfInterpolate2(e[cE_110_111]->Normal, g110, d110, g111, d111,
                                    I->Level);
            }
            active |= cM_110_111;
          }

          if(active) {
            switch (mode) {
            case cIsosurfaceMode::triangles_tri_normals:
            case cIsosurfaceMode::triangles_grad_normals:
              code =
                (i000 ? 0x01 : 0) |
                (i001 ? 0x02 : 0) |
                (i010 ? 0x04 : 0) |
                (i011 ? 0x08 : 0) |
                (i100 ? 0x10 : 0) |
                (i101 ? 0x20 : 0) | (i110 ? 0x40 : 0) | (i111 ? 0x80 : 0);
              eidx = I->EdgeStart[code];
              while(1) {
                idx = I->Edge[eidx];
                if(idx < 0)
                  break;

                /* assemble a triangle from these three points */

                VLACheck(I->Tri, TriangleType, n_tri);
                tt = I->Tri + n_tri;
                tt->p[0] = e[idx];
                tt->p[1] = e[I->Edge[eidx + 1]];
                tt->p[2] = e[I->Edge[eidx + 2]];

                VLACheck(I->PtLink, PointLinkType, n_link + 3);

                /* link this triangle into the points */

                I->PtLink[n_link].tri = n_tri;
                I->PtLink[n_link].link = tt->p[0]->Link;
                tt->p[0]->Link = n_link;
                n_link++;

                I->PtLink[n_link].tri = n_tri;
                I->PtLink[n_link].link = tt->p[1]->Link;
                tt->p[1]->Link = n_link;
                n_link++;

                I->PtLink[n_link].tri = n_tri;
                I->PtLink[n_link].link = tt->p[2]->Link;
                tt->p[2]->Link = n_link;
                n_link++;
                n_tri++;
                eidx += 3;
              }
              break;
            case cIsosurfaceMode::lines:
              VLACheck(*vert, float, (n_vert * 3) + 200);

              code =
                (i000 ? 0x01 : 0) |
                (i001 ? 0x02 : 0) |
                (i010 ? 0x04 : 0) |
                (i011 ? 0x08 : 0) |
                (i100 ? 0x10 : 0) |
                (i101 ? 0x20 : 0) | (i110 ? 0x40 : 0) | (i111 ? 0x80 : 0);
              eidx = I->EdgeStart[code];
              while(1) {
                idx = I->Edge[eidx];
                if(idx < 0)
                  break;
                copy3fn(e[idx]->Point, (*vert) + (n_vert * 3));
                n_vert++;
                copy3fn(e[I->Edge[eidx + 1]]->Point, (*vert) + (n_vert * 3));
                n_vert++;
                copy3fn(e[I->Edge[eidx + 1]]->Point, (*vert) + (n_vert * 3));
                n_vert++;
                copy3fn(e[I->Edge[eidx + 2]]->Point, (*vert) + (n_vert * 3));
                n_vert++;
                copy3fn(e[I->Edge[eidx + 2]]->Point, (*vert) + (n_vert * 3));
                n_vert++;
                copy3fn(e[idx]->Point, (*vert) + (n_vert * 3));
                n_vert++;
                eidx += 3;
              }
              break;
            case cIsosurfaceMode::dots:
            default:           /* dots */
              VLACheck(*vert, float, (n_vert * 3) + 200);

              if(active & cM_000_001) {
                copy3fn(e[cE_000_001]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_000_010) {
                copy3fn(e[cE_000_010]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_000_011) {
                copy3fn(e[cE_000_011]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_000_100) {
                copy3fn(e[cE_000_100]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_000_101) {
                copy3fn(e[cE_000_101]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_000_110) {
                copy3fn(e[cE_000_110]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_000_111) {
                copy3fn(e[cE_000_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }

              if(active & cM_001_011) {
                copy3fn(e[cE_001_011]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_001_101) {
                copy3fn(e[cE_001_101]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_001_111) {
                copy3fn(e[cE_001_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }

              if(active & cM_010_011) {
                copy3fn(e[cE_010_011]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_010_110) {
                copy3fn(e[cE_010_011]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_010_111) {
                copy3fn(e[cE_010_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }

              if(active & cM_100_101) {
                copy3fn(e[cE_100_101]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_100_110) {
                copy3fn(e[cE_100_110]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_100_111) {
                copy3fn(e[cE_100_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }

              if(active & cM_011_111) {
                copy3fn(e[cE_011_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_101_111) {
                copy3fn(e[cE_101_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }
              if(active & cM_110_111) {
                copy3fn(e[cE_101_111]->Point, (*vert) + (n_vert * 3));
                n_vert += 1;
              }

              break;
            }
            n_active++;
          }
        }
      }

  switch (mode) {
  case cIsosurfaceMode::triangles_tri_normals:
  case cIsosurfaceMode::triangles_grad_normals:
    /* compute area-weighted normal */
    for(a = 0; a < n_tri; a++) {
      float *v0, *v1, *v2;
      float vt1[3], vt2[3];

      tt = I->Tri + a;
      v0 = tt->p[0]->Point;
      v1 = tt->p[1]->Point;
      v2 = tt->p[2]->Point;
      tt->done = false;         /* init */

      subtract3f(v0, v2, vt1);
      subtract3f(v1, v2, vt2);
      if (side != cIsosurfaceSide::back) {
        cross_product3f(vt2, vt1, tt->n);
      } else {
        cross_product3f(vt1, vt2, tt->n);
      }
    }
    /* compute (or lookup) normals at active points */
    for(a = 0; a < n_tri; a++) {
      float v[3];
      float dp;
      tt = I->Tri + a;
      for(b = 0; b < 3; b++) {
        if(!tt->p[b]->NormalFlag) {
          zero3f(v);
          idx = tt->p[b]->Link;
          while(idx > 0) {
            add3f(I->Tri[I->PtLink[idx].tri].n, v, v);
            idx = I->PtLink[idx].link;
          }
          if (mode == cIsosurfaceMode::triangles_grad_normals) {
            /* gradient-based normals */
            dp = dot_product3f(v, tt->p[b]->Normal);
            if(dp < 0.0F) {
              invert3f(tt->p[b]->Normal);
              normalize3f(tt->p[b]->Normal);
            } else if(dp == 0.0F) {     /* fall back on triangle normal */
              normalize23f(v, tt->p[b]->Normal);
            } else {
              normalize3f(tt->p[b]->Normal);
            }
          } else { /* triangle-based normals */
            normalize23f(v, tt->p[b]->Normal);
          }
          tt->p[b]->NormalFlag = true;
        }
      }
    }
    if (mode == cIsosurfaceMode::triangles_tri_normals) {
      /* if we're using triangle normals, then 
         do an additional averaging cycle with no weighting */
      for(a = 0; a < n_tri; a++) {
        tt = I->Tri + a;
        add3f(tt->p[0]->Normal, tt->p[1]->Normal, tt->n);
        add3f(tt->p[2]->Normal, tt->n, tt->n);
        normalize3f(tt->n);
        tt->p[0]->NormalFlag = false;
        tt->p[1]->NormalFlag = false;
        tt->p[2]->NormalFlag = false;
      }
      /* compute normals at active points */
      for(a = 0; a < n_tri; a++) {
        float v[3];
        tt = I->Tri + a;
        for(b = 0; b < 3; b++) {
          if(!tt->p[b]->NormalFlag) {
            zero3f(v);
            idx = tt->p[b]->Link;
            while(idx > 0) {
              add3f(I->Tri[I->PtLink[idx].tri].n, v, v);
              idx = I->PtLink[idx].link;
            }
            normalize23f(v, tt->p[b]->Normal);
            tt->p[b]->NormalFlag = true;
          }
        }
      }
    }

    /* Need to move the points now, right? */

    /* if we are carving, then exclude triangles outside region */
    if (carvehelper) {
      for(a = 0; a < n_tri; a++) {
        tt = I->Tri + a;
        if (carvehelper->is_excluded( //
                tt->p[0]->Point,      //
                tt->p[1]->Point,      //
                tt->p[2]->Point)) {
          tt->done = true;    /* exclude this triangle from the surface */
        }
      }
    }

    /* now create triangle strips (not yet optimal) */
    for(a = 0; a < n_tri; a++) {
      tt = I->Tri + a;
      n_start = n_vert;
      if(!tt->done) {

        VLACheck(*vert, float, (n_vert * 3) + 200);

        if (side != cIsosurfaceSide::back) {
          /* switch order around to get "correct" triangles */

          copy3fn(tt->p[1]->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(tt->p[1]->Point, (*vert) + (n_vert * 3));
          n_vert++;

          copy3fn(tt->p[0]->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(tt->p[0]->Point, (*vert) + (n_vert * 3));
          n_vert++;

          copy3fn(tt->p[2]->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(tt->p[2]->Point, (*vert) + (n_vert * 3));
          n_vert++;

          p0 = tt->p[0];
          p1 = tt->p[2];

        } else {

          copy3fn(tt->p[0]->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(tt->p[0]->Point, (*vert) + (n_vert * 3));
          n_vert++;

          copy3fn(tt->p[1]->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(tt->p[1]->Point, (*vert) + (n_vert * 3));
          n_vert++;

          copy3fn(tt->p[2]->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(tt->p[2]->Point, (*vert) + (n_vert * 3));
          n_vert++;

          p0 = tt->p[1];
          p1 = tt->p[2];
        }

        tt->done = true;

        while(1) {
          p2 = NULL;
          idx = p1->Link;
          while(idx > 0) {
            tt = I->Tri + I->PtLink[idx].tri;

            if(!tt->done) {
              if((tt->p[0] == p0) && (tt->p[1] == p1)) {
                p2 = tt->p[2];
                break;
              }

              if((tt->p[1] == p0) && (tt->p[2] == p1)) {
                p2 = tt->p[0];
                break;
              }

              if((tt->p[2] == p0) && (tt->p[0] == p1)) {
                p2 = tt->p[1];
                break;
              }

              if((tt->p[1] == p0) && (tt->p[0] == p1)) {
                p2 = tt->p[2];
                break;
              }

              if((tt->p[2] == p0) && (tt->p[1] == p1)) {
                p2 = tt->p[0];
                break;
              }

              if((tt->p[0] == p0) && (tt->p[2] == p1)) {
                p2 = tt->p[1];
                break;
              }
            }
            idx = I->PtLink[idx].link;
          }

          if(!p2)
            break;
          tt->done = true;
          VLACheck(*vert, float, (n_vert * 3) + 200);
          copy3fn(p2->Normal, (*vert) + (n_vert * 3));
          n_vert++;
          copy3fn(p2->Point, (*vert) + (n_vert * 3));
          n_vert++;
          p0 = p1;
          p1 = p2;
        }
      }
      if(n_vert > n_start) {
        *strip_l.check(n_strip++) = n_vert - n_start;
      }
    }
    I->TotPrim += n_tri;
    break;

  case cIsosurfaceMode::lines:
    /* Need to move the points now, right? */

    if(n_vert > n_start) {
      *strip_l.check(n_strip++) = n_vert - n_start;
    }
    break;
  case cIsosurfaceMode::dots:
  default:

    /* Need to move the points now, right? */

    if(n_vert > n_start) {
      *strip_l.check(n_strip++) = n_vert - n_start;
    }
    break;
  }
  return (n_vert);
}


/*===========================================================================*/
static int TetsurfCodeVertices(CTetsurf * II)
{
  CTetsurf *I = II;
  int i, j, k;
  int b0, b1;
  int flag1 = false;
  int flag2 = false;
  b0 = 1;
  if(I->Level < 0.0F)
    b0 = 0;
  b1 = 1 - b0;

  for(i = 0; i < I->Max[0]; i++)
    for(j = 0; j < I->Max[1]; j++)
      for(k = 0; k < I->Max[2]; k++) {
        if((O3(I->Data, i, j, k, I->CurOff) > I->Level)) {
          I3(I->VertexCodes, i, j, k) = b0;
          flag1 = true;
        } else {
          I3(I->VertexCodes, i, j, k) = b1;
          flag2 = true;
        }
      }
  return (flag1 && flag2);
}
