

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
#include"Matrix.h"
#include"CGO.h"
#include"Err.h"
#include"Crystal.h"
#include"Feedback.h"
#include"PConv.h"
#include"Vector.h"

#include <algorithm>

PyObject *CrystalAsPyList(const CCrystal * I)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(2);
    PyList_SetItem(result, 0, PConvFloatArrayToPyList(I->Dim, 3));
    PyList_SetItem(result, 1, PConvFloatArrayToPyList(I->Angle, 3));
  }
  return (PConvAutoNone(result));
}

int CrystalFromPyList(CCrystal * I, PyObject * list)
{
  int ok = true, rok = true;
  int ll = 0;
  if(ok)
    ok = (I != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  rok = ok;
  if(ok && (ll > 0))
    ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 0), I->Dim, 3);
  if(ok && (ll > 1))
    ok = PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 1), I->Angle, 3);

  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */

  return (rok);
}

CCrystal::CCrystal(PyMOLGlobals* GParam) : G(GParam)
{
}

void CrystalDump(const CCrystal * I)
{
  PyMOLGlobals *G = I->G;
  int i;

  PRINTF
    " Crystal: Unit Cell         %8.3f %8.3f %8.3f\n", I->Dim[0], I->Dim[1], I->Dim[2]
    ENDF(G);
  PRINTF
    " Crystal: Alpha Beta Gamma  %8.3f %8.3f %8.3f\n",
    I->Angle[0], I->Angle[1], I->Angle[2]
    ENDF(G);
  PRINTF " Crystal: RealToFrac Matrix\n" ENDF(G);
  for(i = 0; i < 3; i++) {
    PRINTF " Crystal: %9.4f %9.4f %9.4f\n",
      I->realToFrac()[i * 3], I->realToFrac()[i * 3 + 1], I->realToFrac()[i * 3 + 2]
      ENDF(G);
  }
  PRINTF " Crystal: FracToReal Matrix\n" ENDF(G);
  for(i = 0; i < 3; i++) {
    PRINTF
      " Crystal: %9.4f %9.4f %9.4f\n",
      I->fracToReal()[i * 3], I->fracToReal()[i * 3 + 1], I->fracToReal()[i * 3 + 2]
      ENDF(G);
  }
  PRINTF " Crystal: Unit Cell Volume %8.0f.\n", I->unitCellVolume() ENDF(G);

}

static float unitCellVertices[][3] = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, // bottom 4 vertices
				       {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} }; // top 4 vertices
static int unitCellLineIndices[] = { 0, 1, 1, 2, 2, 3, 3, 0,   // bottom 4 lines
				     4, 5, 5, 6, 6, 7, 7, 4,   // top 4 lines
				     0, 4, 1, 5, 2, 6, 3, 7 }; // 4 connector lines

CGO *CrystalGetUnitCellCGO(const CCrystal * I)
{
  PyMOLGlobals *G = I->G;
  float v[3];
  CGO *cgo = NULL;
  if(I) {
    cgo = CGONew(G);
    CGODisable(cgo, GL_LIGHTING);

    float *vertexVals = CGODrawArrays(cgo, GL_LINES, CGO_VERTEX_ARRAY, 24);
    for (int pl = 0; pl < 24 ; pl++){
      transform33f3f(I->fracToReal(), unitCellVertices[unitCellLineIndices[pl]], v);
      copy3f(v, &vertexVals[pl*3]);
    }

    CGOEnable(cgo, GL_LIGHTING);
    CGOStop(cgo);
  }
  return (cgo);
}

float CCrystal::unitCellVolume() const
{
  return determinant33f(fracToReal());
}

void CCrystal::setDims(const float* dims)
{
  setDims(dims[0], dims[1], dims[2]);
}

void CCrystal::setDims(float a, float b, float c)
{
  invalidateMatrices();

  Dim[0] = a;
  Dim[1] = b;
  Dim[2] = c;
}

void CCrystal::setAngles(const float* angles)
{
  setAngles(angles[0], angles[1], angles[2]);
}

void CCrystal::setAngles(float alpha, float beta, float gamma)
{
  invalidateMatrices();

  Angle[0] = alpha ? alpha : 90;
  Angle[1] = beta ? beta : 90;
  Angle[2] = gamma ? gamma : 90;
}

const float* CCrystal::realToFrac() const
{
  if (!m_RealToFracValid) {
    double f2r[9], r2f[9];
    std::copy_n(fracToReal(), 9, f2r);
    xx_matrix_invert(r2f, f2r, 3);

    auto this_mutable = const_cast<CCrystal*>(this);
    this_mutable->m_RealToFracValid = true;
    std::copy_n(r2f, 9, this_mutable->RealToFrac);
  }

  return RealToFrac;
}

const float* CCrystal::fracToReal() const
{
  if (!m_FracToRealValid) {
    auto this_mutable = const_cast<CCrystal*>(this);
    this_mutable->m_FracToRealValid = true;
    auto* f2r = this_mutable->FracToReal;

    identity33f(f2r);

    if (Dim[0] && Dim[1] && Dim[2] && Angle[0] && Angle[1] && Angle[2]) {
      float cabg[3]; // Cosine of axis angle
      float sabg[3]; // Singe of axis angle

      for (int i = 0; i < 3; ++i) {
        cabg[i] = cos(Angle[i] * PI / 180.0);
        sabg[i] = sin(Angle[i] * PI / 180.0);
      }

      float const cabgs[3] = {
          (cabg[1] * cabg[2] - cabg[0]) / (sabg[1] * sabg[2]),
          (cabg[2] * cabg[0] - cabg[1]) / (sabg[2] * sabg[0]),
          (cabg[0] * cabg[1] - cabg[2]) / (sabg[0] * sabg[1]),
      };

      auto const sabgs1 = sqrt1d(1.0 - cabgs[0] * cabgs[0]);

      f2r[0] = Dim[0];
      f2r[1] = cabg[2] * Dim[1];
      f2r[2] = cabg[1] * Dim[2];
      f2r[4] = sabg[2] * Dim[1];
      f2r[5] = -sabg[1] * cabgs[0] * Dim[2];
      f2r[8] = sabg[1] * sabgs1 * Dim[2];
    }
  }

  return FracToReal;
}

void CCrystal::setFracToReal(const float* f2c)
{
  m_FracToRealValid = true;
  m_RealToFracValid = false;

  std::copy_n(f2c, 9, FracToReal);

  float f2c_t[9];
  transpose33f33f(f2c, f2c_t);

  Dim[0] = length3f(f2c_t + 0);
  Dim[1] = length3f(f2c_t + 3);
  Dim[2] = length3f(f2c_t + 6);

  Angle[0] = rad_to_deg(get_angle3f(f2c_t + 3, f2c_t + 6));
  Angle[1] = rad_to_deg(get_angle3f(f2c_t + 0, f2c_t + 6));
  Angle[2] = rad_to_deg(get_angle3f(f2c_t + 0, f2c_t + 3));
}

bool CCrystal::isSuspicious() const
{
  return is_identityf(3, fracToReal(), R_SMALL4) || //
         unitCellVolume() < R_SMALL4;
}
