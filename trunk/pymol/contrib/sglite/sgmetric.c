/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef SG_GLOBAL
#include "sglite.h"


/* BEGIN: same code as in unitcell.c
 */

static void MxTranspose(const double *M, int nr, int nc, double *Mt)
{
  int  ir, ic;
  range1(ir, nr)
    range1(ic, nc)
      Mt[ic * nr + ir] = M[ir * nc + ic];
}

static void MxMultiply(const double *A, const double *B,
                       int ma, int na, int nb,
                       double *AB)
{
  /* AB[ma, nb] = A[ma, na] * B[na, nb] */
  int  i, j, k;
  range1(i, ma) {
    range1(k, nb) {
      *AB = 0.;
      range1(j, na) *AB += A[i * na + j] * B[j * nb + k];
      AB++;
    }
  }
}

static void getRtGR(const double G[9], const double R[9], double RtGR[9])
{
  double  Rt[9], GR[9];

  MxTranspose(R, 3, 3, Rt);
  MxMultiply(G, R, 3, 3, 3, GR);
  MxMultiply(Rt, GR, 3, 3, 3, RtGR);
}

/* END: same code as in unitcell.c
 */


int CheckMetricalMatrix(const T_SgOps *SgOps, const double *G,
                        double tolerance)
{
  /* for all R in the representative set of SgOps,
     assert Transpose[R].G.R == G
   */

  int     iSMx, i;
  double  R[9], RtGR[9], delta;

  if (tolerance < 0.) tolerance = 1.e-4;

  range2(iSMx, 1, SgOps->nSMx) {
    rangei(9) R[i] = (double) SgOps->SMx[iSMx].s.R[i];
    getRtGR(G, R, RtGR);
    rangei(9) {
      delta = RtGR[i] - G[i];
      if (delta < 0.) delta *= -1.;
      if (delta > tolerance) {
        SetSgError(
          "Error: metrical matrix is incompatible with symmetry operations");
        return -1;
      }
    }
  }

  return 0;
}
