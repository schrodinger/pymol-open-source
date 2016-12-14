#ifndef TNT_LAPACK_HPP
#define TNT_LAPACK_HPP

extern "C"
{
#include "f2c.h"
#include "clapack.h"
}


/**
  This file is an exmaple of how one integrates TNT arrays with external libraries.
  
  C-Lapack is a translation of the LAPACK Fortran code into C.  This package contains
  methods for solving linear systems of equations and eigenvalue problems.

  Lapack_LinearSolve() calls one of the C-Lapack drivers to solve the equation
  AX = B, where A and B are given.  The function returns X as the solution.

  To compile and link this file, you NEED the C-Lapack library installed.
  See http://www.netlib.org/clapack/ for details.  In particular, your
  link needs to include LAPACK library, as well as the BLAS library,
  TMGLIB and libF77.a library of F2C.
  

*/
namespace TNT
{

  Fortran_Array2D<double> Lapack_LinearSolve( const Fortran_Array2D<double> &A,  
      const Fortran_Array2D<double> &B)
  {
      
      integer m = A.dim1();
      integer n = A.dim2();

      if (m != n)
        return Fortran_Array2D<double>();

      integer nrhs = B.dim2();

      if (n != B.dim1())
        return Fortran_Array2D<double>();


      Fortran_Array2D<double> X(B.copy());
      Fortran_Array2D<double> lu(A.copy());
      Fortran_Array1D<integer>  ipiv(m, integer(0));

      integer  info=0;

      dgesv_( &n, &nrhs, &lu(1,1), &m, &ipiv(1), &X(1,1), &n, &info);

      return X;
  }



}


#endif
// TNT_LAPACK_HPP

