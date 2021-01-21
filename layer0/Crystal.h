
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
#ifndef _H_Crystal
#define _H_Crystal

#include "os_python.h"

class CGO;

/**
 * Representation of a unit cell.
 *
 * Calling setDims() or setAngles() will define an x-axis and xy-plane aligned
 * coordinate system. Calling setFracToReal() can define an arbitrarily oriented
 * coordinate system.
 */
struct CCrystal {
  PyMOLGlobals* G;

private:
  float Dim[3] = {1.0, 1.0, 1.0};
  float Angle[3] = {90.0, 90.0, 90.0}; /* stored in degrees for convenience */
  float RealToFrac[9];
  float FracToReal[9];

  bool m_RealToFracValid = false;
  bool m_FracToRealValid = false;

  void invalidateMatrices()
  {
    m_RealToFracValid = false;
    m_FracToRealValid = false;
  }

public:
  /// Cell edge lengths in Angstroms (a, b, c)
  const float* dims() const { return Dim; }
  void setDims(const float*);
  void setDims(float, float, float);

  /// Cell angles in degrees (alpha, beta, gamma)
  const float* angles() const { return Angle; }
  void setAngles(const float*);
  void setAngles(float, float, float);

  /// 3x3 transformation matrix from cartesian space to fractional space
  const float* realToFrac() const;

  /// 3x3 transformation matrix from fractional space to cartesian space
  const float* fracToReal() const;

  /// Set the fracToReal() matrix. Also updates all other parameters
  /// (realToFrac(), dims(), angles()).
  void setFracToReal(const float*);

  /// Unit cell volume in Angstrom^3
  float unitCellVolume() const;

  /// True if this cell is an orthogonal 1x1x1 or is singular
  bool isSuspicious() const;

  CCrystal(PyMOLGlobals* GParam);

  friend void CrystalDump(const CCrystal*);
  friend int CrystalFromPyList(CCrystal*, PyObject*);
  friend PyObject* CrystalAsPyList(const CCrystal*);
};

CGO* CrystalGetUnitCellCGO(const CCrystal* I);

#endif
