/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif


#undef SG_GLOBAL
#include "sglite.h"


int SenseOfRotation(const int *R, int Rtype, const int *EV)
{
  /* M.B. Boisen, Jr. & G.V. Gibbs
     Mathematical Crystallography, Revised Edition 1990
     pp. 348-349, 354-356
   */

  int  f, trace;


  f = 1; if (Rtype < 0) f = -1;
  trace = f * (R[0] + R[4] + R[8]);

  if (trace == 3 || trace == -1) return 0; /* 1-fold or 2-fold */

  if (EV[1] == 0 && EV[2] == 0) {
    if (EV[0] * f * R[7] > 0)
      return 1;
  }
  else
    if (f * (R[3] * EV[2] - R[6] * EV[1]) > 0)
      return 1;

  return -1;
}


int GetRtype(const int *RotMx)
{
  int deter = deterRotMx(RotMx);

  if (deter == -1 || deter == 1)
  {
    switch (traceRotMx(RotMx))
    {
      case -3:                  return -1;
      case -2:                  return -6;
      case -1: if (deter == -1) return -4;
               else             return  2;
      case  0: if (deter == -1) return -3;
               else             return  3;
      case  1: if (deter == -1) return -2;
               else             return  4;
      case  2:                  return  6;
      case  3:                  return  1;
    }
  }

  return 0;
}


int SetRotMxInfo(const int *R, T_RMxI *RI)
{
  int        Rtype, ProperOrder, i;
  int        M_ProperR[9], RmI[9];
  const int   *ProperR;


  if (RI) {
    RI->Rtype = 0;
    rangei(3) RI->EV[i] = 0;
    RI->SenseOfRotation = 0;
  }

      Rtype = GetRtype(R);
  if (Rtype == 0)
    return 0;

  if (RI)
  {
    ProperR = R;

        ProperOrder = Rtype;
    if (ProperOrder < 0) {
        ProperOrder *= -1;
      rangei(9) M_ProperR[i] = -R[i];
      ProperR = M_ProperR;
    }

    if (ProperOrder > 1) {
      SetRminusI(ProperR, RmI, 0);
      if (iRowEchelonFormT(RmI, 3, 3, NULL, 0) != 2) return 0;
      if (SolveHomRE2(RmI, RI->EV) != 0) return 0;
      RI->SenseOfRotation = SenseOfRotation(R, Rtype, RI->EV);
    }

    RI->Rtype = Rtype;
  }

  return Rtype;
}


int OrderOfRtype(int Rtype)
{
  if (Rtype > 0) return  Rtype;
  if (Rtype % 2) return -Rtype * 2;
                 return -Rtype;
}


int MakeCumRMx(const int *R, int Rtype, int *CumRMx)
{
  int  Order;
  int  MxA[9], MxB[9];
  int  *RR, *RRR, *Swp, iO, i;


  InitRotMx(CumRMx, 1);

      Order = OrderOfRtype(Rtype);
  if (Order > 1)
  {
    RR = (int *) R;
    RRR = MxA;

    for (iO = 1;;)
    {
      rangei(9) CumRMx[i] += RR[i];

      if (++iO == Order)
        break;

      RotMxMultiply(RRR, R, RR);
      if (RR == R) RR = MxB;
      Swp = RR; RR = RRR; RRR = Swp;
    }
  }

  return Order;
}


int Set_wI_Tr(const int *R, const int *T, const T_RMxI *RI,
              int wI[3] /* STBF */,
              int Tr[3] /* CTBF */)
{
  int     Mul, Mx[9], wl[3], i;
  int     P[9], Pwl[3];
  T_RMxI  BufRI[1];


  if (T == NULL) T = &R[9];

           rangei(3) wI[i] = 0;
  if (Tr)  rangei(3) Tr[i] = 0;

  if (RI == NULL) {
    if (SetRotMxInfo(R, BufRI) == 0) return -1;
    RI = BufRI;
  }

  Mul = MakeCumRMx(R, RI->Rtype, Mx);
  RotMx_t_Vector(wI, Mx, T, 0);
  if (ChangeBaseFactor(wI, Mul, wI, 1, 3) != 0) return 1;

  if (Tr == NULL) return 0;

  rangei(3) wl[i] = -(T[i] - wI[i]) * (CTBF / STBF);

  SetRminusI(R, Mx, 0);
  (void) IdentityMat(P, 3);
  (void) iRowEchelonFormT(Mx, 3, 3, P, 3);
  iMxMultiply(Pwl, P, wl, 3, 3, 1);
  i = iREBacksubst(Mx, Pwl, 3, 3, Tr, NULL);
  if (i < 1) return -1;
  if (i > 1) return  1;

  return 0;
}
