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


int CB_IT(const int SignI, const int T[3],
          const T_RTMx *CBMx, const T_RTMx *InvCBMx,
          int BC_T[3])
{
  /* (C|V)( I|T)(C^-1|W)=( C|CT+V)(C^-1|W)=( I|CT+V+CW)=( I|C(T+W)+V)
     (C|V)(-I|T)(C^-1|W)=(-C|CT+V)(C^-1|W)=(-I|CT+V-CW)=(-I|C(T-W)+V)
   */

  int  TpW[3], i;

  rangei(3) TpW[i] = T[i] * (CTBF / STBF) + SignI * InvCBMx->s.T[i];

  RotMx_t_Vector(BC_T, CBMx->s.R, TpW, 0);

  rangei(3) {
        BC_T[i] += CBMx->s.T[i] * CRBF;
    if (BC_T[i] %  (CRBF * (CTBF / STBF))) return IE(-1);
        BC_T[i] /= (CRBF * (CTBF / STBF));
  }

  return 0;
}


int CB_RMx(int CRiC[9],
           const int CBMxR[9], const int RMx[9], const int InvCBMxR[9])
{
  int  BufMx[9];

  RotMxMultiply(BufMx, RMx,   InvCBMxR);
  RotMxMultiply(CRiC,  CBMxR, BufMx);

  if (ChangeBaseFactor(CRiC, CRBF * CRBF, CRiC, 1, 9) != 0) {
    SetSgError(
        "Error: Change-of-basis -> out of rotation-base-factor range");
    return -1;
  }

  return 0;
}


int CB_SMx(T_RTMx *CSiC,
           const T_RTMx *CBMx, const T_RTMx *SMx, const T_RTMx *InvCBMx)
{
  T_RTMx  BufMx;


  RTMxMultiply(&BufMx, SMx,  InvCBMx, CTBF / STBF, 0);
  RTMxMultiply(CSiC,   CBMx, &BufMx,  CRBF,        CRBF * CTBF);

  if (ChangeBaseFactor(CSiC->s.R, CRBF * CRBF,
                       CSiC->s.R, 1, 9) != 0) {
    SetSgError(
        "Error: Change-of-basis -> out of rotation-base-factor range");
    return -1;
  }

  if (ChangeBaseFactor(CSiC->s.T, CRBF * (CTBF / STBF),
                       CSiC->s.T, 1, 3) != 0) {
    SetSgError(
        "Error: Change-of-basis -> out of translation-base-factor range");
    return -1;
  }

  return 0;
}


int CB_SgLTr(const T_SgOps *SgOps,
             const T_RTMx *CBMx, const T_RTMx *InvCBMx,
             T_SgOps *BC_SgOps)
{
  int     iBV, iLLTr, i;
  int     BVT[3], BC_T[3];


  range1(iBV, 3) {
    rangei(3) BVT[i] = (i == iBV ? STBF : 0);
    if (CB_IT(1, BVT, CBMx, InvCBMx, BC_T) != 0) return -1;
    if (ExpSgLTr(BC_SgOps, BC_T) < 0) return -1;
  }

  range1(iLLTr, SgOps->nLTr) {
    if (CB_IT(1, SgOps->LTr[iLLTr].v, CBMx, InvCBMx, BC_T) != 0) return -1;
    if (ExpSgLTr(BC_SgOps, BC_T) < 0) return -1;
  }

  return 0;
}


static int CB_SgSMx(const T_SgOps *SgOps,
                    const T_RTMx *CBMx, const T_RTMx *InvCBMx,
                    T_SgOps *BC_SgOps)
{
  int     iLSMx;
  int     BC_T[3];
  T_RTMx  BC_SMx[1];


  if (SgOps->fInv == 2) {
    if (CB_IT(-1, SgOps->InvT, CBMx, InvCBMx, BC_T) != 0) return -1;
    if (ExpSgInv(BC_SgOps, BC_T) < 0) return -1;
  }

  range2(iLSMx, 1, SgOps->nSMx) {
    if (CB_SMx(BC_SMx, CBMx, &SgOps->SMx[iLSMx], InvCBMx) != 0) return -1;
    if (ExpSgSMx(BC_SgOps, BC_SMx) < 0) return -1;
  }

  return 0;
}


int CB_SgOps(const T_SgOps *SgOps,
             const T_RTMx *CBMx, const T_RTMx *InvCBMx,
             T_SgOps *BC_SgOps)
{
  if (CB_SgLTr(SgOps, CBMx, InvCBMx, BC_SgOps) != 0) return -1;
  if (CB_SgSMx(SgOps, CBMx, InvCBMx, BC_SgOps) != 0) return -1;

  return 0;
}
