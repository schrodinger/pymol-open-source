

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

#include <random>

#include"os_python.h"
#include"os_predef.h"
#include"os_std.h"

#include"Isosurf.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Symmetry.h"
#include"Vector.h"
#include"Feedback.h"
#include"PConv.h"
#include"P.h"
#include"Util.h"

#define Trace_OFF

#define O3(field,P1,P2,P3,offs) ((field)->get<float>((P1)+offs[0],(P2)+offs[1],(P3)+offs[2]))

#define O3Ptr(field,P1,P2,P3,offs) ((field)->ptr<float>((P1)+offs[0],(P2)+offs[1],(P3)+offs[2]))

#define O4Ptr(field,P1,P2,P3,P4,offs) ((field)->ptr<float>((P1)+offs[0],(P2)+offs[1],(P3)+offs[2],P4))

#define I3(field,P1,P2,P3) ((field)->get<int>(P1,P2,P3))

#define I4(field,P1,P2,P3,P4) ((field)->get<int>(P1,P2,P3,P4))

typedef struct PointType {
  float Point[3];
  int NLink;
  struct PointType* Link[4];
} PointType;

#define EdgePtPtr(field,P2,P3,P4,P5) ((field)->ptr(P2,P3,P4,P5))

#define EdgePt(field,P2,P3,P4,P5) ((field)->get(P2,P3,P4,P5))

struct CIsosurf {
  PyMOLGlobals *G;
  CFieldTyped<int> *VertexCodes;
  CFieldTyped<int> *ActiveEdges;
  CFieldTyped<PointType> *Point;
  int NLine;
  int Skip;
  int AbsDim[3], CurDim[3], CurOff[3];
  int Max[3];
  CField *Coord, *Data;
  float Level;
  int Code[256];

  pymol::vla<int>* Num = nullptr;
  int NSeg;
  pymol::vla<float>* Line = nullptr;

};

static int IsosurfAlloc(PyMOLGlobals * G, CIsosurf * II);
static void IsosurfPurge(CIsosurf * II);
static int IsosurfCurrent(CIsosurf * II);
static int IsosurfCodeVertices(CIsosurf * II);
static int IsosurfFindActiveEdges(CIsosurf * II);
static int IsosurfFindLines(CIsosurf * II);
static int IsosurfDrawLines(CIsosurf * II);
static void IsosurfCode(CIsosurf * II, const char *bits1, const char *bits2);
static int IsosurfDrawPoints(CIsosurf * II);
static int IsosurfPoints(CIsosurf * II);
static int IsosurfGradients(PyMOLGlobals * G, CSetting * set1, CSetting * set2,
                            CIsosurf * II, Isofield * field,
                            int *range, float min_level, float max_level);

#define IsosurfSubSize		64

static void _IsosurfFree(CIsosurf * I)
{
  FreeP(I);
}

void IsosurfFree(PyMOLGlobals * G)
{
  _IsosurfFree(G->Isosurf);
  G->Isosurf = nullptr;
}


/*===========================================================================*/
PyObject *IsosurfAsPyList(PyMOLGlobals * G, Isofield * field)
{
  PyObject *result = nullptr;

  result = PyList_New(4);

  PyList_SetItem(result, 0, PConvIntArrayToPyList(field->dimensions, 3));
  PyList_SetItem(result, 1, PyInt_FromLong(field->save_points));
  PyList_SetItem(result, 2, FieldAsPyList(G, field->data.get()));
  if(field->save_points)
    PyList_SetItem(result, 3, FieldAsPyList(G, field->points.get()));
  else
    PyList_SetItem(result, 3, PConvAutoNone(NULL));
  return (PConvAutoNone(result));
}


/*===========================================================================*/
inline
static void IsosurfInterpolate(CIsosurf * I, float *v1, float *l1, float *v2, float *l2,
                               float *pt)
{
  float ratio;
  ratio = (I->Level - *l1) / (*l2 - *l1);
  pt[0] = v1[0] + (v2[0] - v1[0]) * ratio;
  pt[1] = v1[1] + (v2[1] - v1[1]) * ratio;
  pt[2] = v1[2] + (v2[2] - v1[2]) * ratio;
}


/*===========================================================================*/
Isofield *IsosurfNewFromPyList(PyMOLGlobals * G, PyObject * list)
{
  int ok = true;
  int dim4[4];
  int a;

  Isofield *result = nullptr;
  if(ok)
    ok = (list != nullptr);
  if(ok)
    ok = PyList_Check(list);
  /* TO ENABLE BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok)
    ok = (result = new Isofield()) != nullptr;
  if(ok)
    ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 0), result->dimensions, 3);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &result->save_points);
  if(ok){
    result->data.reset(FieldNewFromPyList_From_List(G, list, 2));
    ok = result->data != nullptr;
  }
  if(ok) {
    if(result->save_points) {
      result->points.reset(FieldNewFromPyList_From_List(G, list, 3));
      ok = result->points != nullptr;
    }
    else {
      for(a = 0; a < 3; a++)
        dim4[a] = result->dimensions[a];
      dim4[3] = 3;
      result->points.reset(new CFieldTyped<float>(dim4, 4));
      ok = result->points != nullptr;
    }
  }
  if(!ok) {
    DeleteP(result);
  }
  return (result);
}

/*===========================================================================*/
void IsofieldComputeGradients(PyMOLGlobals * G, Isofield * field)
{
  int dim[4];
  int a, b, c;
  CField *data = field->data.get();

  if(!field->gradients) {

    /* compute gradients relative to grid axis spacing */

    for(a = 0; a < 3; a++)
      dim[a] = field->dimensions[a];
    dim[3] = 3;
    field->gradients.reset(new CFieldTyped<float>(dim, 4));
    auto gradients = field->gradients.get();
    dim[3] = 3;

    /* bulk internal */

    for(a = 1; a < (dim[0] - 1); a++) {
      for(b = 1; b < (dim[1] - 1); b++) {
        for(c = 1; c < (dim[2] - 1); c++) {
          F4(gradients, a, b, c, 0) =
            (F3(data, a + 1, b, c) - F3(data, a - 1, b, c)) / 2.0F;
          F4(gradients, a, b, c, 1) =
            (F3(data, a, b + 1, c) - F3(data, a, b - 1, c)) / 2.0F;
          F4(gradients, a, b, c, 2) =
            (F3(data, a, b, c + 1) - F3(data, a, b, c - 1)) / 2.0F;
        }
      }
    }

    for(a = 0; a < dim[0]; a += (dim[0] - 1)) {

      /* 'a' faces */
      for(b = 1; b < (dim[1] - 1); b++) {
        for(c = 1; c < (dim[2] - 1); c++) {
          if(!a) {
            F4(gradients, a, b, c, 0) = (F3(data, a + 1, b, c) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 0) = (F3(data, a, b, c) - F3(data, a - 1, b, c));
          }
          F4(gradients, a, b, c, 1) =
            (F3(data, a, b + 1, c) - F3(data, a, b - 1, c)) / 2.0F;
          F4(gradients, a, b, c, 2) =
            (F3(data, a, b, c + 1) - F3(data, a, b, c - 1)) / 2.0F;
        }
      }

      /* 'c' edges and all eight corners */
      for(b = 0; b < dim[1]; b += (dim[1] - 1)) {
        for(c = 0; c < dim[2]; c++) {

          if(!a) {
            F4(gradients, a, b, c, 0) = (F3(data, a + 1, b, c) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 0) = (F3(data, a, b, c) - F3(data, a - 1, b, c));
          }

          if(!b) {
            F4(gradients, a, b, c, 1) = (F3(data, a, b + 1, c) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 1) = (F3(data, a, b, c) - F3(data, a, b - 1, c));
          }

          if(!c) {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c + 1) - F3(data, a, b, c));
          } else if(c < (dim[2] - 1)) {
            F4(gradients, a, b, c, 2) =
              (F3(data, a, b, c + 1) - F3(data, a, b, c - 1)) / 2.0F;
          } else {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c) - F3(data, a, b, c - 1));
          }
        }
      }

      /* 'b' edges  */
      for(c = 0; c < dim[2]; c += (dim[2] - 1)) {
        for(b = 1; b < (dim[1] - 1); b++) {
          if(!a) {
            F4(gradients, a, b, c, 0) = (F3(data, a + 1, b, c) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 0) = (F3(data, a, b, c) - F3(data, a - 1, b, c));
          }

          F4(gradients, a, b, c, 1) =
            (F3(data, a, b + 1, c) - F3(data, a, b - 1, c)) / 2.0F;

          if(!c) {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c + 1) - F3(data, a, b, c));
          } else if(c < (dim[2] - 1)) {
            F4(gradients, a, b, c, 2) =
              (F3(data, a, b, c + 1) - F3(data, a, b, c - 1)) / 2.0F;
          } else {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c) - F3(data, a, b, c - 1));
          }
        }
      }
    }

    for(b = 0; b < dim[1]; b += (dim[1] - 1)) {

      for(a = 1; a < (dim[0] - 1); a++) {

        /* 'b' faces */

        for(c = 1; c < (dim[2] - 1); c++) {
          F4(gradients, a, b, c, 0) =
            (F3(data, a + 1, b, c) - F3(data, a - 1, b, c)) / 2.0F;
          if(!b) {
            F4(gradients, a, b, c, 1) = (F3(data, a, b + 1, c) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 1) = (F3(data, a, b, c) - F3(data, a, b - 1, c));
          }
          F4(gradients, a, b, c, 2) =
            (F3(data, a, b, c + 1) - F3(data, a, b, c - 1)) / 2.0F;
        }

        /* 'a' edges */

        for(c = 0; c < dim[2]; c += (dim[2] - 1)) {
          F4(gradients, a, b, c, 0) =
            (F3(data, a + 1, b, c) - F3(data, a - 1, b, c)) / 2.0F;
          if(!b) {
            F4(gradients, a, b, c, 1) = (F3(data, a, b + 1, c) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 1) = (F3(data, a, b, c) - F3(data, a, b - 1, c));
          }
          if(!c) {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c + 1) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c) - F3(data, a, b, c - 1));
          }
        }
      }
    }

    /* 'c' faces */

    for(c = 0; c < dim[2]; c += (dim[2] - 1)) {
      for(a = 1; a < (dim[0] - 1); a++) {
        for(b = 1; b < (dim[1] - 1); b++) {
          F4(gradients, a, b, c, 0) =
            (F3(data, a + 1, b, c) - F3(data, a - 1, b, c)) / 2.0F;
          F4(gradients, a, b, c, 1) =
            (F3(data, a, b + 1, c) - F3(data, a, b - 1, c)) / 2.0F;
          if(!c) {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c + 1) - F3(data, a, b, c));
          } else {
            F4(gradients, a, b, c, 2) = (F3(data, a, b, c) - F3(data, a, b, c - 1));
          }
        }
      }
    }
  }
}


/*===========================================================================*/
Isofield::Isofield(PyMOLGlobals * G, const int * const dims)
{
  int dim4[4];
  std::copy_n(dims, 3, dim4);
  dim4[3] = 3;

  /* Warning: ...FromPyList also allocs and inits from the heap */

  data.reset(new CFieldTyped<float>(dims, 3));
  points.reset(new CFieldTyped<float>(dim4, 4));
  std::copy_n(dims, 3, dimensions);
}

/*===========================================================================*/
static void IsosurfCode(CIsosurf * II, const char *bits1, const char *bits2)
{
  CIsosurf *I = II;
  int c;
  int b;
  int sum1, sum2;

  c = 0;
  while(bits1[c])
    c++;
  c--;
  sum1 = 0;
  b = 1;
  while(c >= 0) {
    if(bits1[c] == '1')
      sum1 = sum1 + b;
    b = b + b;
    c--;
  }

  c = 0;
  while(bits2[c])
    c++;
  c--;
  sum2 = 0;
  b = 1;
  while(c >= 0) {
    if(bits2[c] == '1')
      sum2 = sum2 + b;
    b = b + b;
    c--;
  }

  I->Code[sum1] = sum2;
#ifdef Trace
  printf("IsosurfCode: %s (%i) -> %s (%i)\n", bits1, sum1, bits2, sum2);
#endif
}


/*===========================================================================*/
static CIsosurf *IsosurfNew(PyMOLGlobals * G)
{
  int c;
  CIsosurf *I = pymol::calloc<CIsosurf>(1);
  I->G = G;
  I->VertexCodes = nullptr;
  I->ActiveEdges = nullptr;
  I->Point = nullptr;
  I->Line = nullptr;
  I->Skip = 0;
  for(c = 0; c < 256; c++)
    I->Code[c] = -1;


/*___  
 | / |
 |/  |
 |___|
 32
*/
  IsosurfCode(I, "10000010", "100000");
  IsosurfCode(I, "01000001", "100000");


/*___  
 | \ |
 |  \|
 |___|
 16
*/
  IsosurfCode(I, "10010000", "010000");
  IsosurfCode(I, "01100000", "010000");


/*___  
 |   |
 |  /|
 |_/_|
 8
*/
  IsosurfCode(I, "00101000", "001000");
  IsosurfCode(I, "00010100", "001000");


/*___  
 |   |
 |\  |
 |_\_|
 4
*/
  IsosurfCode(I, "00001001", "000100");
  IsosurfCode(I, "00000110", "000100");


/*___  
 | \ |
 |\ \|
 |_\_|
 16+4=20
*/

  IsosurfCode(I, "01101001", "010100");


/*___  
 | / |
 |/ /|
 |_/_|
 32+8=40
*/
  IsosurfCode(I, "10010110", "101000");


/*___  
 | | |
 | | |
 |_|_|
 2
*/
  IsosurfCode(I, "10001000", "000010");
  IsosurfCode(I, "01000100", "000010");


/*___  
 |   |
 |---|
 |___|
 1
*/
  IsosurfCode(I, "00100010", "000001");
  IsosurfCode(I, "00010001", "000001");

  return (I);
}

int IsosurfInit(PyMOLGlobals * G)
{
  G->Isosurf = IsosurfNew(G);
  return 1;
}


/*===========================================================================*/
int IsosurfExpand(Isofield * field1, Isofield * field2, CCrystal * cryst,
                  CSymmetry * sym, int *range)
{
  float rmn[3], rmx[3];
  float imn[3], imx[3];
  float fstep[3], frange[3];
  int field1max[3];
  int expanded = false;
  int missing = false;

  field1max[0] = field1->dimensions[0] - 1;
  field1max[1] = field1->dimensions[1] - 1;
  field1max[2] = field1->dimensions[2] - 1;
  {
    int a;
    for(a = 0; a < 3; a++) {
      rmn[a] = F4(field1->points, 0, 0, 0, a);
      rmx[a] = F4(field1->points, field1max[0], field1max[1], field1max[2], a);
    }
  }

  /* get min/max extents of map1 in fractional space */

  transform33f3f(cryst->realToFrac(), rmn, imn);
  transform33f3f(cryst->realToFrac(), rmx, imx);

  /* compute step size */

  subtract3f(imx, imn, frange);

  fstep[0] = frange[0] / field1max[0];
  fstep[1] = frange[1] / field1max[1];
  fstep[2] = frange[2] / field1max[2];

  /* compute coordinate points for second field */

  if (int nMat = sym->getNSymMat()) {
    int i, j, k;
    int i_stop, j_stop, k_stop;
    float frac[3];

    i_stop = field2->dimensions[0];
    j_stop = field2->dimensions[1];
    k_stop = field2->dimensions[2];
    for(i = 0; i < i_stop; i++) {
      frac[0] = imn[0] + fstep[0] * (i + range[0]);
      for(j = 0; j < j_stop; j++) {
        frac[1] = imn[1] + fstep[1] * (j + range[1]);
        for(k = 0; k < k_stop; k++) {
          float average = 0.0F;
          float extrapolate_average = 0.0F;
          int cnt = 0;
          int extrapolate_cnt = 0;

          /* first compute the coordinate */

          float *ptr = F4Ptr(field2->points, i, j, k, 0);
          frac[2] = imn[2] + fstep[2] * (k + range[2]);
          transform33f3f(cryst->fracToReal(), frac, ptr);

          /* then compute the value at the coordinate */

          for(int n = nMat - 1; n >= 0; n--) {
            const float *matrix = sym->getSymMat(n);
            float test_frac[3];

            transform44f3f(matrix, frac, test_frac);

            /* we're assuming that the identity matrix appears in the list */

            test_frac[0] -= imn[0];
            test_frac[1] -= imn[1];
            test_frac[2] -= imn[2];

            test_frac[0] -= (int) floor(test_frac[0] + R_SMALL4);
            test_frac[1] -= (int) floor(test_frac[1] + R_SMALL4);
            test_frac[2] -= (int) floor(test_frac[2] + R_SMALL4);

            {
              int a, b, c;
              float x, y, z;
              a = (int) (test_frac[0] / fstep[0]);
              b = (int) (test_frac[1] / fstep[1]);
              c = (int) (test_frac[2] / fstep[2]);
              x = (test_frac[0] / fstep[0]) - a;
              y = (test_frac[1] / fstep[1]) - b;
              z = (test_frac[2] / fstep[2]) - c;

              if((a >= 0) && (b >= 0) && (c >= 0) &&
                 (a <= (field1max[0] + 1)) && (b <= (field1max[1] + 1))
                 && (c <= (field1max[2] + 1))) {
                while(a >= field1max[0]) {
                  a--;
                  x += 1.0F;
                }
                while(b >= field1max[1]) {
                  b--;
                  y += 1.0F;
                }
                while(c >= field1max[2]) {
                  c--;
                  z += 1.0F;
                }

                {
                  const float sloppy_1 = 1.0F + R_SMALL4;
                  if((x <= sloppy_1) && (y <= sloppy_1) && (z <= sloppy_1)) {
                    if(!expanded) {
                      if((matrix[0] != 1.0F) || /* not identity matrix */
                         (matrix[5] != 1.0F) ||
                         (matrix[10] != 1.0F) || (matrix[15] != 1.0F) ||
                         /* and not inside source map */
                         ((imn[0] - frac[0]) > R_SMALL4)
                         || ((frac[0] - imx[0]) > R_SMALL4)
                         || ((imn[1] - frac[1]) > R_SMALL4)
                         || ((frac[1] - imx[1]) > R_SMALL4)
                         || ((imn[2] - frac[2]) > R_SMALL4)
                         || ((frac[2] - imx[2]) > R_SMALL4)) {
                        expanded = true;
                      }
                    }
                    /* found at least one point obtained through symmetry or priodicity */
                    if(x > 1.0F)
                      x = 1.0F;
                    if(y > 1.0F)
                      y = 1.0F;
                    if(z > 1.0F)
                      z = 1.0F;
                    average += FieldInterpolatef(field1->data.get(), a, b, c, x, y, z);
                    cnt++;
                  } else {
                    /* allow 1 cell of extrapolation -- this saves us
                       when someone issues map_double and is then technically
                       missing a plane of data at the edge of the cell */
                    if(((x-1.0F) < sloppy_1) && ((y-1.0F) < sloppy_1) && ((z-1.0F) < sloppy_1)) {
                      if(x > 1.0F)
                        x = 1.0F;
                      if(y > 1.0F)
                        y = 1.0F;
                      if(z > 1.0F)
                        z = 1.0F;
                      extrapolate_average += FieldInterpolatef(field1->data.get(), a, b, c, x, y, z);
                      extrapolate_cnt++;
                    }
                  }
                }
              }
            }
          }
          if(cnt) {
            F3(field2->data, i, j, k) = average / cnt;
          } else if(extrapolate_cnt) {
            F3(field2->data, i, j, k) = extrapolate_average / extrapolate_cnt;
          } else {
            missing = true;
            F3(field2->data, i, j, k) = 0.0F;   /* complain? */
          }
        }
      }
    }
  }
  if(expanded) {
    if(missing) {
      return -1;
    } else {
      return 1;
    }
  } else {
    return 0;
  }
}

int IsosurfGetRange(PyMOLGlobals * G, Isofield * field,
                    CCrystal * cryst, float *mn, float *mx, int *range, int clamp)
{
  float rmn[3], rmx[3];
  float imn[3], imx[3];
  float mix[24], imix[24];
  int a, b;
  int clamped = false;          /* clamped? */
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
    if(range[a] < 0) {
      if(clamp)
        range[a] = 0;
      clamped = true;
    }
    if(range[a] > field->dimensions[a]) {
      if(clamp)
        range[a] = field->dimensions[a];
      clamped = true;
    }
    if(range[a + 3] < 0) {
      if(clamp)
        range[a + 3] = 0;
      clamped = true;
    }
    if(range[a + 3] > field->dimensions[a]) {
      if(clamp)
        range[a + 3] = field->dimensions[a];
      clamped = true;
    }
  }
  PRINTFD(G, FB_Isosurface)
    " IsosurfGetRange: returning range: %d %d %d %d %d %d\n",
    range[0], range[1], range[2], range[3], range[4], range[5]
    ENDFD;
  return clamped;
}


/*===========================================================================*/
int IsosurfVolume(PyMOLGlobals* G, CSetting* set1, CSetting* set2,
    Isofield* field, float level, pymol::vla<int>& num, pymol::vla<float>& vert,
    int* range, cIsomeshMode mode, int skip, float alt_level)
{
  int ok = true;
  CIsosurf *I;
  if(PIsGlutThread()) {
    I = G->Isosurf;
  } else {
    I = IsosurfNew(G);
  }
  CHECKOK(ok, I);
  {
    int Steps[3];
    int c, i, j, k;
    int x, y, z;
    int range_store[6];
    I->Num = std::addressof(num);
    I->Line = std::addressof(vert);
    I->Skip = skip;
    if(range) {
      for(c = 0; c < 3; c++) {
        I->AbsDim[c] = field->dimensions[c];
        I->CurDim[c] = IsosurfSubSize + 1;
        Steps[c] = ((range[3 + c] - range[c]) - 2) / IsosurfSubSize + 1;
      }
    } else {
      range = range_store;
      for(c = 0; c < 3; c++) {
        range[c] = 0;
        range[3 + c] = field->dimensions[c];
        I->AbsDim[c] = field->dimensions[c];
        I->CurDim[c] = IsosurfSubSize + 1;
        Steps[c] = (I->AbsDim[c] - 2) / IsosurfSubSize + 1;
      }
    }

    I->Coord = field->points.get();
    I->Data = field->data.get();
    I->Level = level;
    if(ok)
      ok = IsosurfAlloc(G, I);

    I->NLine = 0;
    I->NSeg = 0;
    I->Num->check(I->NSeg);
    (*I->Num)[I->NSeg] = I->NLine;

    if(ok) {
      switch (mode) {
      case cIsomeshMode::gradient:
        ok = IsosurfGradients(G, set1, set2, I, field, range, level, alt_level);
        IsosurfPurge(I);
        break;
      default:
        for(i = 0; i < Steps[0]; i++) {
          for(j = 0; j < Steps[1]; j++) {
            for(k = 0; k < Steps[2]; k++) {
              if(ok) {
                I->CurOff[0] = IsosurfSubSize * i;
                I->CurOff[1] = IsosurfSubSize * j;
                I->CurOff[2] = IsosurfSubSize * k;
                for(c = 0; c < 3; c++)
                  I->CurOff[c] += range[c];
                for(c = 0; c < 3; c++) {
                  I->Max[c] = range[3 + c] - I->CurOff[c];
                  if(I->Max[c] > (IsosurfSubSize + 1))
                    I->Max[c] = (IsosurfSubSize + 1);
                }
                if(!(i || j || k)) {
                  for(x = 0; x < I->Max[0]; x++)
                    for(y = 0; y < I->Max[1]; y++)
                      for(z = 0; z < I->Max[2]; z++)
                        for(c = 0; c < 3; c++)
                          EdgePt(I->Point, x, y, z, c).NLink = 0;
                }
#ifdef Trace
                for(c = 0; c < 3; c++)
                  printf(" IsosurfVolume: c: %i CurOff[c]: %i Max[c] %i\n", c,
                         I->CurOff[c], I->Max[c]);
#endif

                if(ok)
                  switch (mode) {
                  case cIsomeshMode::isomesh:      /* standard mode - want lines */
                    ok = IsosurfCurrent(I);
                    break;
                  case cIsomeshMode::isodot:      /* point mode - just want points on the isosurface */
                    ok = IsosurfPoints(I);
                    break;
                  }
                if(G->Interrupt) {
                  ok = false;
                }
              }
            }
          }
        }
        IsosurfPurge(I);
        break;
      }
    }

    if(mode != cIsomeshMode::isomesh) {
      PRINTFB(G, FB_Isomesh, FB_Blather)
        " IsosurfVolume: Surface generated using %d dots.\n", I->NLine ENDFB(G);
    } else {
      PRINTFB(G, FB_Isomesh, FB_Blather)
        " IsosurfVolume: Surface generated using %d lines.\n", I->NLine ENDFB(G);
    }

    if(!ok) {
      I->NLine = 0;
      I->NSeg = 0;
    }
    /* shrinks sizes for more efficient RAM usage */

    I->Line->resize(I->NLine * 3);
    I->Num->resize(I->NSeg + 1);
    (*I->Num)[I->NSeg] = 0;        /* important - must terminate the segment list */

    if(!PIsGlutThread()) {
      _IsosurfFree(I);
    }
  }
  return (ok);
}


/*===========================================================================*/
static int IsosurfAlloc(PyMOLGlobals * G, CIsosurf * II)
{
  CIsosurf *I = II;

  int ok = true;
  int dim4[4];
  int a;
  for(a = 0; a < 3; a++)
    dim4[a] = I->CurDim[a];
  dim4[3] = 3;

  I->VertexCodes = new CFieldTyped<int>(I->CurDim, 3);
  I->ActiveEdges = new CFieldTyped<int>(dim4, 4);
  I->Point = new CFieldTyped<PointType>(dim4, 4);
  if(!(I->VertexCodes && I->ActiveEdges && I->Point)) {
    IsosurfPurge(I);
    ok = false;
  }
#ifdef Trace
  printf(" IsosurfAlloc: ok: %i\n", ok);
#endif
  return (ok);
}


/*===========================================================================*/
static void IsosurfPurge(CIsosurf * II)
{
  CIsosurf *I = II;
  DeleteP(I->VertexCodes);
  DeleteP(I->ActiveEdges);
  DeleteP(I->Point);
}


/*===========================================================================*/
static int IsosurfCurrent(CIsosurf * II)
{
  CIsosurf *I = II;
  int ok = true;
  if(IsosurfCodeVertices(I)) {
    if(ok)
      ok = IsosurfFindActiveEdges(I);
    if(ok)
      ok = IsosurfFindLines(I);
    if(ok)
      ok = IsosurfDrawLines(I);
  }
  return (ok);
}


/*===========================================================================*/
static int IsosurfPoints(CIsosurf * II)
{
  CIsosurf *I = II;
  int ok = true;
  if(IsosurfCodeVertices(I)) {
    if(ok)
      ok = IsosurfFindActiveEdges(I);
    if(ok)
      ok = IsosurfDrawPoints(I);
  }
  return (ok);
}


/*===========================================================================*/

static int IsosurfGradients(PyMOLGlobals * G, CSetting * set1, CSetting * set2,
                            CIsosurf * II, Isofield * field,
                            int *range, float min_level, float max_level)
{
  CIsosurf *I = II;
  int ok = true;

  /* use local copies for performance reasons */

  int n_seg = I->NSeg;
  int n_line = I->NLine;
  auto& i_line = *I->Line;
  auto i_data = I->Data;
  auto& i_num = *I->Num;

  /* get cascaded state, object, or global settings */

  int spacing = SettingGet_i(G, set1, set2, cSetting_gradient_spacing);
  float step_size = SettingGet_f(G, set1, set2, cSetting_gradient_step_size);
  float max_walk = SettingGet_f(G, set1, set2, cSetting_gradient_max_length);
  float min_walk = SettingGet_f(G, set1, set2, cSetting_gradient_min_length);
  float min_slope = SettingGet_f(G, set1, set2, cSetting_gradient_min_slope);
  float min_dot = SettingGet_f(G, set1, set2, cSetting_gradient_normal_min_dot);
  float symmetry = SettingGet_f(G, set1, set2, cSetting_gradient_symmetry);

  int symmetry_flag = false;    /* are we searching for symmetric segments? */
  if(symmetry != 0.0F)
    symmetry_flag = true;       /* (very slow process) */
  if(symmetry > 1.0F)
    symmetry = 1.0F / symmetry;

  /* clamp dangerous parameters */

  if(step_size < 0.01F)
    step_size = 0.01F;

  if(min_slope < 0.00001F)
    min_slope = 0.00001F;

  /* make sure we have gradients available for map */

  if(!field->gradients)
    IsofieldComputeGradients(G, field);

  /* and that map has a minimum size */

  if(field->gradients) {

    /* locals for performance */

    CField *gradients = field->gradients.get();
    CField *points = field->points.get();

    /* flags marking excluded regions to avoid (currently wasteful) */
    int *flag = nullptr;

    /* variable length array for recording segment paths */
    int *active_cell = VLAlloc(int, 1000);

    /* ordered list of coordinates for processing */
    int *order = nullptr;

    int range_size;             /* total points in region being drawn */
    int range_dim[3];           /* dimension of drawn region */
    int flag_stride[3];         /* stride values for flag array */

    range_dim[0] = (range[3] - range[0]);
    range_dim[1] = (range[4] - range[1]);
    range_dim[2] = (range[5] - range[2]);

    range_size = range_dim[0] * range_dim[1] * range_dim[2];
    if (ok)
      flag = pymol::calloc<int>(range_size);
    CHECKOK(ok, flag);
    flag_stride[0] = 1;
    flag_stride[1] = range_dim[0];
    flag_stride[2] = range_dim[0] * range_dim[1];

    if (ok)
      order = pymol::calloc<int>(3 * range_size);
    CHECKOK(ok, order);
    if(order && flag && (range_dim[0] > 1) && (range_dim[1] > 1) && (range_dim[2] > 1)) {

      {
        /* compute approximate cell spacing */

        float average_cell_axis_dist;
        float *pos[4];
        pos[0] = Ffloat4p(I->Coord, 0, 0, 0, 0);
        pos[1] = Ffloat4p(I->Coord, 1, 0, 0, 0);
        pos[2] = Ffloat4p(I->Coord, 0, 1, 0, 0);
        pos[3] = Ffloat4p(I->Coord, 0, 0, 1, 0);

        average_cell_axis_dist = (float) ((diff3f(pos[0], pos[1]) +
                                           diff3f(pos[0], pos[2]) +
                                           diff3f(pos[0], pos[3])) / 3.0);

        /* scale parameters into cell units */

        max_walk /= average_cell_axis_dist;
        min_walk /= average_cell_axis_dist;
        step_size /= average_cell_axis_dist;
        min_slope *= average_cell_axis_dist;
      }

      {
        /* generate randomized list of cell coordinates */

        /* always use same seed for same volume */
        std::mt19937 mt(range_size);
        std::uniform_real_distribution<float> dist{};
        {
          /* fill */
          int i, j, k, *p = order;
          for(k = range[2]; k < range[5]; k++) {
            for(j = range[1]; j < range[4]; j++) {
              for(i = range[0]; i < range[3]; i++) {
                p[0] = i;
                p[1] = j;
                p[2] = k;
                p += 3;
              }
            }
          }
        }
        {
          /* shuffle */
          int a;
          for(a = 0; a < range_size; a++) {
            int *p = order + 3 * (int) (range_size * dist(mt));
            int *q = order + 3 * (int) (range_size * dist(mt));
            int t0 = p[0], t1 = p[1], t2 = p[2];
            p[0] = q[0];
            p[1] = q[1];
            p[2] = q[2];
            q[0] = t0;
            q[1] = t1;
            q[2] = t2;
          }
        }
      }

      {
        /* now draw our lines */

        int a;
        int *start_locus = order;
        float prev_grad_normal[3] = { 0.0F, 0.0F, 0.0F };
        for(a = 0; a < range_size; a++) {
          int n_active_cell = 0;        /* how many cells have we traversed */
          float walk = max_walk;        /* distance remaining to travel */

          int abort_n_line = n_line;    /* for backtracking */
          int abort_n_seg = n_seg;

          int pass;             /* the pass are we on */

          float symmetry_max = FLT_MIN; /* if we're trying to get symmetric segments */
          float symmetry_min = FLT_MAX;

          for(pass = 0; pass < 2; pass++) {     /* one pass down the gradient, one up */

            int have_prev = false;      /* flag & storage for previous gradient & locus */
            int *prev_locus = nullptr;

            int locus[3];       /* what cell are we in? */
            float fract[3] = { 0.0F, 0.0F, 0.0F };      /* where in the cell are we? */
            int n_vert = 0;

            locus[0] = start_locus[0];
            locus[1] = start_locus[1];
            locus[2] = start_locus[2];

            for(;;) {

              /* master segment extension loop, using "break" to exit */

              {
                /* normalize locus and fract before each new step */
                int done = false;
                int b;
                for(b = 0; b < 3; b++) {

                  while(fract[b] < 0.0F) {      /* force fract >= 0.0 */
                    fract[b] += 1.0F;
                    locus[b]--;
                  }
                  while(fract[b] >= 1.0F) {     /* force fract into [0.0-1.0) */
                    fract[b] -= 1.0F;
                    locus[b]++;
                  }
                  while(locus[b] > (range[b + 3] - 2)) {        /* above range? done */
                    if(fract[b] <= 0.0F) {
                      locus[b]--;
                      fract[b] += 1.0F;
                      if(locus[b] < range[b]) { /* below range? done */
                        done = true;
                        break;
                      }
                    } else {
                      done = true;
                      break;
                    }
                  }
                  while(locus[b] < range[b]) {  /* below range? done */
                    if(fract[b] > 1.0F) {
                      locus[b]++;
                      fract[b] -= 1.0F;
                      if(locus[b] > (range[b + 3] - 2)) {       /* above range? done */
                        done = true;
                        break;
                      }
                    } else {
                      done = true;
                      break;
                    }
                  }

                }
                if(done)
                  break;
              }

              /* have we moved cells? */

              if((!have_prev) || (have_prev && ((locus[0] != prev_locus[0]) ||
                                                (locus[1] != prev_locus[1]) ||
                                                (locus[2] != prev_locus[2])))) {
                /* above: prev_locus may be nullptr, so relying upon shortcut logic eval */

                /* stop if we hit a flagged cell (flag always in lower corner) */

                if(*(flag + (((locus[0] - range[0]) * flag_stride[0]) +
                             ((locus[1] - range[1]) * flag_stride[1]) +
                             ((locus[2] - range[2]) * flag_stride[2])))) {
                  break;
                }

              }

              {
                /* stop if level exceeds desired ranges */
                float level = FieldInterpolatef(i_data, locus[0], locus[1], locus[2],
                                                fract[0], fract[1], fract[2]);
                if((level < min_level) || (level > max_level))
                  break;
                if(symmetry_flag) {
                  if(symmetry_min > level)
                    symmetry_min = level;
                  if(symmetry_max < level)
                    symmetry_max = level;
                }
              }

              {
                /* interpolate gradient relative to grid */

                float interp_gradient[3];

                FieldInterpolate3f(gradients, locus, fract, interp_gradient);

                if(length3f(interp_gradient) < min_slope) {
                  /* if region is too flat, then bail */
                  break;
                }

                {
                  /* add a line vertex at this point */

                  {
                    float *f;
                    VLACheck(i_line, float, n_line * 3 + 2);
                    f = i_line + (n_line * 3);
                    FieldInterpolate3f(points, locus, fract, f);
                    n_line++;
                    n_vert++;
                  }

                  /* record locus for subsequent oblation */

                  if((!have_prev) || (have_prev && ((locus[0] != prev_locus[0]) ||
                                                    (locus[1] != prev_locus[1]) ||
                                                    (locus[2] != prev_locus[2])))) {

                    VLACheck(active_cell, int, n_active_cell * 3 + 2);
                    {
                      int *xrd = active_cell + (n_active_cell * 3);
                      xrd[0] = locus[0];
                      xrd[1] = locus[1];
                      xrd[2] = locus[2];
                      n_active_cell++;
                      prev_locus = xrd; /* warning: volatile pointer */
                    }
                  }

                  /* adjust length of gradient vector */

                  normalize3f(interp_gradient);

                  /* make sure gradient isn't too divergent to take another step */

                  if(have_prev) {
                    if(dot_product3f(interp_gradient, prev_grad_normal) < min_dot)
                      break;
                  }

                  /* take another step */

                  copy3f(interp_gradient, prev_grad_normal);

                  /* scale and flip sign */

                  if(pass) {
                    scale3f(interp_gradient, -step_size, interp_gradient);
                  } else {
                    scale3f(interp_gradient, step_size, interp_gradient);
                  }

                  /* record progress */

                  walk -= step_size;

                  /* leave if max_walk is reached */

                  if(walk < 0.0F) {
                    break;
                  } else {
                    /* otherwise move */
                    add3f(interp_gradient, fract, fract);
                  }
                  have_prev = true;
                }
              }
            }                   /* for */

            if(n_vert < 2) {    /* quash isolated vertices */
              if(n_vert) {
                n_line = i_num[n_seg];
              }
            } else if(n_line != i_num[n_seg]) { /* or retain count of new line segment */
              VLACheck(i_num, int, n_seg + 1);
              i_num[n_seg] = n_line - i_num[n_seg];
              n_seg++;
              i_num[n_seg] = n_line;
            }

          }
          {
            int abort_segment = false;
            if(symmetry_flag) {
              if((symmetry_max * symmetry_min) >= 0.0F) /* abort if not both +/- pot. sampled */
                abort_segment = true;
              else {
                float symmetry_ratio = (float) (fabs(symmetry_max) / fabs(symmetry_min));
                if(symmetry_ratio > 1.0F)
                  symmetry_ratio = 1.0F / symmetry_ratio;
                if(symmetry_ratio < symmetry)   /* abort if +/- weren't close enough in magnitude */
                  abort_segment = true;
              }
            }
            if((max_walk - walk) < min_walk) {  /* ignore too-short segments */
              abort_segment = true;
            }

            if(abort_segment) {
              n_seg = abort_n_seg;
              n_line = abort_n_line;
              i_num[n_seg] = n_line;
            } else {
              /* otherwise, keep line and oblate neighborhood */

              int *ac = active_cell;
              int b;
              int cutoff_sq = spacing * spacing;
              for(b = 0; b < n_active_cell; b++) {
                int ii = ac[0], jj = ac[1], kk = ac[2];
                int i0 = ii - spacing;
                int j0 = jj - spacing;
                int k0 = kk - spacing;

                int i1 = ii + spacing + 1;
                int j1 = jj + spacing + 1;
                int k1 = kk + spacing + 1;

                if(i0 < range[0])
                  i0 = range[0];
                if(i1 >= range[3])
                  i1 = range[3] - 1;
                if(j0 < range[1])
                  j0 = range[1];
                if(j1 >= range[4])
                  j1 = range[4] - 1;
                if(k0 < range[2])
                  k0 = range[2];
                if(k1 >= range[5])
                  k1 = range[5] - 1;

                {
                  int i, j, k;
                  int *flag1 = flag + (((i0 - range[0]) * flag_stride[0]) +
                                       ((j0 - range[1]) * flag_stride[1]) +
                                       ((k0 - range[2]) * flag_stride[2]));

                  /* highly optimized spherical flag-fill routine */

                  for(k = k0; k < k1; k++) {
                    int *flag2 = flag1;
                    int kk_sq = (kk - k);
                    kk_sq = kk_sq * kk_sq;

                    for(j = j0; j < j1; j++) {
                      int *flag3 = flag2;
                      int jj_sq = (jj - j);
                      jj_sq = (jj_sq * jj_sq) + kk_sq;

                      if(!(jj_sq > cutoff_sq)) {
                        for(i = i0; i < i1; i++) {
                          if(!*flag3) {
                            int tot_sq = (ii - i);
                            tot_sq = (tot_sq * tot_sq) + jj_sq;
                            if(!(tot_sq > cutoff_sq)) {
                              *flag3 = true;
                            }
                          }
                          flag3++;
                        }       /* for i */
                      }
                      flag2 += flag_stride[1];
                    }           /* for j */
                    flag1 += flag_stride[2];
                  }             /* for k */
                }

                ac += 3;        /* advance to next active cell */
              }                 /* for b in active_cell */
            }
          }
          start_locus += 3;
        }                       /* for a in range_size */
      }

    }
    /* purge memory */
    VLAFreeP(active_cell);
    FreeP(order);
    FreeP(flag);
  }

  /* restore modified local copies */
  I->NLine = n_line;
  I->NSeg = n_seg;
  return (ok);
}


/*===========================================================================*/
static int IsosurfDrawPoints(CIsosurf * II)
{
  CIsosurf *I = II;
  float *a, *b;
  int i, j, k;
  int ok = true;

  if(ok) {
    for(i = 0; i < (I->Max[0] - 1); i++) {
      for(j = 0; j < I->Max[1]; j++) {
        for(k = 0; k < I->Max[2]; k++) {
          if((I3(I->VertexCodes, i, j, k)) && (!I3(I->VertexCodes, i + 1, j, k))) {
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i + 1, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i + 1, j, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 0).Point[0]));

            I->Line->check(I->NLine * 3 + 2);
            a = I->Line->data() + (I->NLine * 3);
            b = &(EdgePt(I->Point, i, j, k, 0).Point[0]);
            *(a++) = *(b++);
            *(a++) = *(b++);
            *a = *b;
            I->NLine++;
          } else if(!(I3(I->VertexCodes, i, j, k)) && (I3(I->VertexCodes, i + 1, j, k))) {
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i + 1, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i + 1, j, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 0).Point[0]));

            I->Line->check(I->NLine * 3 + 2);
            a = I->Line->data() + (I->NLine * 3);
            b = &(EdgePt(I->Point, i, j, k, 0).Point[0]);
            *(a++) = *(b++);
            *(a++) = *(b++);
            *a = *b;
            I->NLine++;
          } else
            I4(I->ActiveEdges, i, j, k, 0) = 0;
        }
      }
      if(I->G->Interrupt) {
        ok = false;
        break;
      }
    }
  }

  if(ok) {
    for(i = 0; i < I->Max[0]; i++) {
      for(j = 0; j < (I->Max[1] - 1); j++) {
        for(k = 0; k < I->Max[2]; k++) {
          if((I3(I->VertexCodes, i, j, k)) && (!I3(I->VertexCodes, i, j + 1, k))) {
            I4(I->ActiveEdges, i, j, k, 1) = 2;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i, j + 1, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j + 1, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 1).Point[0]));

            I->Line->check(I->NLine * 3 + 2);
            a = I->Line->data() + (I->NLine * 3);
            b = &(EdgePt(I->Point, i, j, k, 1).Point[0]);
            *(a++) = *(b++);
            *(a++) = *(b++);
            *a = *b;
            I->NLine++;

          } else if(!(I3(I->VertexCodes, i, j, k)) && (I3(I->VertexCodes, i, j + 1, k))) {
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i, j + 1, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j + 1, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 1).Point[0]));

            I->Line->check(I->NLine * 3 + 2);
            a = I->Line->data() + (I->NLine * 3);
            b = &(EdgePt(I->Point, i, j, k, 1).Point[0]);
            *(a++) = *(b++);
            *(a++) = *(b++);
            *a = *b;
            I->NLine++;

          }
        }
      }
      if(I->G->Interrupt) {
        ok = false;
        break;
      }
    }
  }

  for(i = 0; i < I->Max[0]; i++) {
    for(j = 0; j < I->Max[1]; j++) {
      for(k = 0; k < (I->Max[2] - 1); k++) {
        if((I3(I->VertexCodes, i, j, k)) && (!I3(I->VertexCodes, i, j, k + 1))) {
          IsosurfInterpolate(I,
                             O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                             O3Ptr(I->Data, i, j, k, I->CurOff),
                             O4Ptr(I->Coord, i, j, k + 1, 0, I->CurOff),
                             O3Ptr(I->Data, i, j, k + 1, I->CurOff),
                             &(EdgePt(I->Point, i, j, k, 2).Point[0]));

          I->Line->check(I->NLine * 3 + 2);
          a = I->Line->data() + (I->NLine * 3);
          b = &(EdgePt(I->Point, i, j, k, 2).Point[0]);
          *(a++) = *(b++);
          *(a++) = *(b++);
          *a = *b;
          I->NLine++;

        } else if(!(I3(I->VertexCodes, i, j, k)) && (I3(I->VertexCodes, i, j, k + 1))) {
          IsosurfInterpolate(I,
                             O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                             O3Ptr(I->Data, i, j, k, I->CurOff),
                             O4Ptr(I->Coord, i, j, k + 1, 0, I->CurOff),
                             O3Ptr(I->Data, i, j, k + 1, I->CurOff),
                             &(EdgePt(I->Point, i, j, k, 2).Point[0]));

          I->Line->check(I->NLine * 3 + 2);
          a = I->Line->data() + (I->NLine * 3);
          b = &(EdgePt(I->Point, i, j, k, 2).Point[0]);
          *(a++) = *(b++);
          *(a++) = *(b++);
          *a = *b;
          I->NLine++;

        }
      }
      if(I->G->Interrupt) {
        ok = false;
        break;
      }
    }
  }
  if(ok) {
    if(I->NLine != (*I->Num)[I->NSeg]) {   /* any new points? */
      I->Num->check(I->NSeg + 1);
      (*I->Num)[I->NSeg] = I->NLine - (*I->Num)[I->NSeg];
      I->NSeg++;
      (*I->Num)[I->NSeg] = I->NLine;
    }
  }
  return (ok);
}


/*===========================================================================*/
static int IsosurfDrawLines(CIsosurf * II)
{
  CIsosurf *I = II;
  int c, i, j, k;
  float *a, *b;
  int ok = true;
  PointType *Cur, *Start, *Next;
  int MaxLinks, MaxL, Cnt;
  int NLink;
#ifdef Trace
  int LCount = 0;
#endif

  for(i = 0; i < I->Max[0]; i++) {
    for(j = 0; j < I->Max[1]; j++) {
      for(k = 0; k < I->Max[2]; k++) {
        for(c = 0; c < 3; c++) {
          Start = EdgePtPtr(I->Point, i, j, k, c);
          while(Start->NLink) {
            Cur = Start;
            I->Line->check(I->NLine * 3 + 2);
            a = I->Line->data() + (I->NLine * 3);
            b = Cur->Point;
            *(a++) = *(b++);
            *(a++) = *(b++);
            *a = *b;
            I->NLine++;

            while(Cur) {
              if(Cur->NLink) {
                Cur->NLink--;
                NLink = Cur->NLink;
                /* Choose point which has most links */
                MaxL = NLink;
                MaxLinks = Cur->Link[MaxL]->NLink;
                Cnt = MaxL - 1;
                while(Cnt >= 0) {
                  if((Cur->Link[Cnt]->NLink) > MaxLinks) {
                    MaxL = Cnt;
                    MaxLinks = Cur->Link[Cnt]->NLink;
                  }
                  Cnt--;
                }
                Next = Cur->Link[MaxL];
                if(MaxL != NLink)
                  Cur->Link[MaxL] = Cur->Link[NLink];
                /* Remove double link */
                Next->NLink--;
                NLink = Next->NLink;
                Cnt = NLink;
                while(Cnt >= 0) {
                  if(Next->Link[Cnt] == Cur)
                    break;
                  else
                    Cnt--;
                }
                if(Cnt >= 0) {
                  if(Cnt != NLink)
                    Next->Link[Cnt] = Next->Link[NLink];
                }
#ifdef Trace
                else
                  printf(" error: IsosurfDrawLines:  can't find double link\n");
#endif

                Cur = Next;
                I->Line->check(I->NLine * 3 + 2);
                a = I->Line->data() + (I->NLine * 3);
                b = Cur->Point;
                *(a++) = *(b++);
                *(a++) = *(b++);
                *a = *b;
                I->NLine++;
              } else {
#ifdef Trace
                LCount++;
#endif
                Cur = nullptr;
                if(I->NLine != (*I->Num)[I->NSeg]) {       /* any new lines? */
                  I->Num->check(I->NSeg + 1);
                  (*I->Num)[I->NSeg] = I->NLine - (*I->Num)[I->NSeg];
                  I->NSeg++;
                  I->Num->check(I->NSeg);
                  (*I->Num)[I->NSeg] = I->NLine;
                }
              }
            }
          }
        }
      }
    }
    if(I->G->Interrupt) {
      ok = false;
      break;
    }
  }
#ifdef Trace
  if(ok)
    printf(" DrawLineCount: %i\n", LCount);
#endif
  return (ok);
}


/*===========================================================================*/
static int IsosurfFindLines(CIsosurf * II)
{
  CIsosurf *I = II;
  int i, j, k, ip1, jp1, kp1;
  int ok = true;
  int index, cod;
  int Max0m1, Max1m1, Max2m1;
  int skip = I->Skip;
#ifdef Trace
  int LCount = 0;
#endif
  PointType *p1, *p2;

  Max0m1 = I->Max[0] - 1;
  Max1m1 = I->Max[1] - 1;
  Max2m1 = I->Max[2] - 1;
  for(i = 0; i < I->Max[0]; i++) {
    for(j = 0; j < I->Max[1]; j++) {
      for(k = 0; k < I->Max[2]; k++) {
        ip1 = i + 1;
        jp1 = j + 1;
        kp1 = k + 1;
        if((j < Max1m1) && (k < Max2m1) && ((!skip) || !(i % skip))) {  /* i-plane */
          index = I4(I->ActiveEdges, i, j, k, 1) << 2;
          index = (index + I4(I->ActiveEdges, i, jp1, k, 2)) << 2;
          index = (index + I4(I->ActiveEdges, i, j, kp1, 1)) << 2;
          index = index + I4(I->ActiveEdges, i, j, k, 2);
          if(index) {
            cod = I->Code[index];
#ifdef Trace
            if(index && (cod < 0))
              printf("IsosurfFindLines: bad index: %i \n", index);
#endif
            while(cod > 0) {
              p1 = nullptr;
              p2 = nullptr;
              switch (cod) {
              case 40:
              case 32:
                cod = cod - 32;
                p1 = EdgePtPtr(I->Point, i, j, k, 1);
                p2 = EdgePtPtr(I->Point, i, j, k, 2);
                break;
              case 20:
              case 16:
                cod = cod - 16;
                p1 = EdgePtPtr(I->Point, i, j, k, 1);
                p2 = EdgePtPtr(I->Point, i, jp1, k, 2);
                break;
              case 8:
                cod = cod - 8;
                p1 = EdgePtPtr(I->Point, i, j, kp1, 1);
                p2 = EdgePtPtr(I->Point, i, jp1, k, 2);
                break;
              case 4:
                cod = cod - 4;
                p1 = EdgePtPtr(I->Point, i, j, kp1, 1);
                p2 = EdgePtPtr(I->Point, i, j, k, 2);
                break;
              case 2:
                cod = cod - 2;
                p1 = EdgePtPtr(I->Point, i, j, k, 1);
                p2 = EdgePtPtr(I->Point, i, j, kp1, 1);
                break;
              case 1:
                cod = cod - 1;
                p1 = EdgePtPtr(I->Point, i, j, k, 2);
                p2 = EdgePtPtr(I->Point, i, jp1, k, 2);
                break;
              default:
                cod = 0;
                p1 = nullptr;
                p2 = nullptr;
                break;
              }
              if(p1 && p2) {
                p1->Link[p1->NLink] = p2;
                p1->NLink++;
                p2->Link[p2->NLink] = p1;
                p2->NLink++;
#ifdef Trace
                LCount++;
#endif
              }
            }
          }
        }
        if((i < Max0m1) && (j < Max1m1) && ((!skip) || !(k % skip))) {  /* k-plane */
          index = I4(I->ActiveEdges, i, j, k, 0) << 2;
          index = (index + I4(I->ActiveEdges, ip1, j, k, 1)) << 2;
          index = (index + I4(I->ActiveEdges, i, jp1, k, 0)) << 2;
          index = index + I4(I->ActiveEdges, i, j, k, 1);
          if(index) {
            cod = I->Code[index];
#ifdef Trace
            if(index && (cod < 0))
              printf("IsosurfFindLines: bad index: %i \n", index);
#endif
            while(cod > 0) {
              switch (cod) {
              case 40:
              case 32:
                cod = cod - 32;
                p1 = EdgePtPtr(I->Point, i, j, k, 0);
                p2 = EdgePtPtr(I->Point, i, j, k, 1);
                break;
              case 20:
              case 16:
                cod = cod - 16;
                p1 = EdgePtPtr(I->Point, i, j, k, 0);
                p2 = EdgePtPtr(I->Point, ip1, j, k, 1);
                break;
              case 8:
                cod = cod - 8;
                p1 = EdgePtPtr(I->Point, i, jp1, k, 0);
                p2 = EdgePtPtr(I->Point, ip1, j, k, 1);
                break;
              case 4:
                cod = cod - 4;
                p1 = EdgePtPtr(I->Point, i, jp1, k, 0);
                p2 = EdgePtPtr(I->Point, i, j, k, 1);
                break;
              case 2:
                cod = cod - 2;
                p1 = EdgePtPtr(I->Point, i, j, k, 0);
                p2 = EdgePtPtr(I->Point, i, jp1, k, 0);
                break;
              case 1:
                cod = cod - 1;
                p1 = EdgePtPtr(I->Point, i, j, k, 1);
                p2 = EdgePtPtr(I->Point, ip1, j, k, 1);
                break;
              default:
                cod = 0;
                p1 = nullptr;
                p2 = nullptr;
                break;
              }
              if(p1 && p2) {
                p1->Link[p1->NLink] = p2;
                p1->NLink++;
                p2->Link[p2->NLink] = p1;
                p2->NLink++;
#ifdef Trace
                LCount++;
#endif
              }
            }
          }
        }
        if((i < Max0m1) && (k < Max2m1) && ((!skip) || !(j % skip))) {  /* j-plane */
          index = I4(I->ActiveEdges, i, j, k, 0) << 2;
          index = (index + I4(I->ActiveEdges, ip1, j, k, 2)) << 2;
          index = (index + I4(I->ActiveEdges, i, j, kp1, 0)) << 2;
          index = index + I4(I->ActiveEdges, i, j, k, 2);
          if(index) {
            cod = I->Code[index];
#ifdef Trace
            if(index && (cod < 0))
              printf("IsosurfFindLines: bad index: %i \n", index);
#endif
            while(cod > 0) {
              switch (cod) {
              case 40:
              case 32:
                cod = cod - 32;
                p1 = EdgePtPtr(I->Point, i, j, k, 0);
                p2 = EdgePtPtr(I->Point, i, j, k, 2);
                break;
              case 20:
              case 16:
                cod = cod - 16;
                p1 = EdgePtPtr(I->Point, i, j, k, 0);
                p2 = EdgePtPtr(I->Point, ip1, j, k, 2);
                break;
              case 8:
                cod = cod - 8;
                p1 = EdgePtPtr(I->Point, i, j, k + 1, 0);
                p2 = EdgePtPtr(I->Point, ip1, j, k, 2);
                break;
              case 4:
                cod = cod - 4;
                p1 = EdgePtPtr(I->Point, i, j, kp1, 0);
                p2 = EdgePtPtr(I->Point, i, j, k, 2);
                break;
              case 2:
                cod = cod - 2;
                p1 = EdgePtPtr(I->Point, i, j, k, 0);
                p2 = EdgePtPtr(I->Point, i, j, kp1, 0);
                break;
              case 1:
                cod = cod - 1;
                p1 = EdgePtPtr(I->Point, i, j, k, 2);
                p2 = EdgePtPtr(I->Point, ip1, j, k, 2);
                break;
              default:
                cod = 0;
                p1 = nullptr;
                p2 = nullptr;
                break;
              }
              if(p1 && p2) {
                p1->Link[p1->NLink] = p2;
                p1->NLink++;
                p2->Link[p2->NLink] = p1;
                p2->NLink++;
#ifdef Trace
                LCount++;
#endif
              }
            }
          }
        }
      }
    }
    if(I->G->Interrupt) {
      ok = false;
      break;
    }

  }
#ifdef Trace
  printf(" IsosurfFindLines: %i lines found\n", LCount);
#endif
  return (ok);
}


/*===========================================================================*/
static int IsosurfFindActiveEdges(CIsosurf * II)
{
  CIsosurf *I = II;
  int i, j, k;
  int ok = true;
#ifdef Trace
  int ECount = 0;
#endif

  if(ok) {
    for(i = 0; i < (I->Max[0] - 1); i++) {
      for(j = 0; j < I->Max[1]; j++) {
        for(k = 0; k < I->Max[2]; k++) {
          if((I3(I->VertexCodes, i, j, k)) && (!I3(I->VertexCodes, i + 1, j, k))) {
#ifdef Trace
            ECount++;
#endif
            I4(I->ActiveEdges, i, j, k, 0) = 2;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i + 1, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i + 1, j, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 0).Point[0]));
          } else if(!(I3(I->VertexCodes, i, j, k)) && (I3(I->VertexCodes, i + 1, j, k))) {
#ifdef Trace
            ECount++;
#endif
            I4(I->ActiveEdges, i, j, k, 0) = 1;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i + 1, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i + 1, j, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 0).Point[0]));
          } else
            I4(I->ActiveEdges, i, j, k, 0) = 0;
        }
      }
      if(I->G->Interrupt) {
        ok = false;
        break;
      }
    }
  }
  if(ok) {
    for(i = 0; i < I->Max[0]; i++) {
      for(j = 0; j < (I->Max[1] - 1); j++) {
        for(k = 0; k < I->Max[2]; k++) {
          if((I3(I->VertexCodes, i, j, k)) && (!I3(I->VertexCodes, i, j + 1, k))) {
#ifdef Trace
            ECount++;
#endif
            I4(I->ActiveEdges, i, j, k, 1) = 2;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i, j + 1, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j + 1, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 1).Point[0]));
          } else if(!(I3(I->VertexCodes, i, j, k)) && (I3(I->VertexCodes, i, j + 1, k))) {
#ifdef Trace
            ECount++;
#endif
            I4(I->ActiveEdges, i, j, k, 1) = 1;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i, j + 1, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j + 1, k, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 1).Point[0]));
          } else {
            I4(I->ActiveEdges, i, j, k, 1) = 0;
          }
        }
      }
      if(I->G->Interrupt) {
        ok = false;
        break;
      }
    }
  }
  if(ok) {
    for(i = 0; i < I->Max[0]; i++) {
      for(j = 0; j < I->Max[1]; j++) {
        for(k = 0; k < (I->Max[2] - 1); k++) {
          if((I3(I->VertexCodes, i, j, k)) && (!I3(I->VertexCodes, i, j, k + 1))) {
#ifdef Trace
            ECount++;
#endif
            I4(I->ActiveEdges, i, j, k, 2) = 2;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i, j, k + 1, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k + 1, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 2).Point[0]));
          } else if(!(I3(I->VertexCodes, i, j, k)) && (I3(I->VertexCodes, i, j, k + 1))) {
#ifdef Trace
            ECount++;
#endif
            I4(I->ActiveEdges, i, j, k, 2) = 1;
            IsosurfInterpolate(I,
                               O4Ptr(I->Coord, i, j, k, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k, I->CurOff),
                               O4Ptr(I->Coord, i, j, k + 1, 0, I->CurOff),
                               O3Ptr(I->Data, i, j, k + 1, I->CurOff),
                               &(EdgePt(I->Point, i, j, k, 2).Point[0]));
          } else {
            I4(I->ActiveEdges, i, j, k, 2) = 0;
          }
        }
      }
      if(I->G->Interrupt) {
        ok = false;
        break;
      }
    }
  }
#ifdef Trace
  printf(" IsosurfFindActiveEdges: %i active edges found\n", ECount);
#endif
  return (ok);
}


/*===========================================================================*/
static int IsosurfCodeVertices(CIsosurf * II)
{
  CIsosurf *I = II;
  int i, j, k;
  int VCount = 0;
  int ok = true;
  for(i = 0; i < I->Max[0]; i++) {
    for(j = 0; j < I->Max[1]; j++) {
      for(k = 0; k < I->Max[2]; k++) {
        if((O3(I->Data, i, j, k, I->CurOff) > I->Level)) {
          I3(I->VertexCodes, i, j, k) = 1;
          VCount++;
        } else {
          I3(I->VertexCodes, i, j, k) = 0;
        }
      }
    }
    if(I->G->Interrupt) {
      ok = false;
      break;
    }
  }
#ifdef Trace
  printf(" IsosurfCodeVertices: %i of %i vertices above level\n", VCount,
         I->CurDim[0] * I->CurDim[1] * I->CurDim[2]);
#endif
  if(!ok)
    VCount = 0;
  return (VCount);
}

/**
 * corner: output buffer of size 8 * 3
 */
void IsofieldGetCorners(PyMOLGlobals * G, Isofield * field, float * corner) {
  CField* points = field->points.get();
  for(int a = 0; a < 8; a++) {
    int i = (a & 1) ? (points->dim[0] - 1) : 0;
    int j = (a & 2) ? (points->dim[1] - 1) : 0;
    int k = (a & 4) ? (points->dim[2] - 1) : 0;
    memcpy(corner + a * 3, F3Ptr(points, i, j, k), 3 * sizeof(float));
  }
}
