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


static void InvT_SMx(const int InvT[3], T_RTMx *SMx)
{
  int  i;

  rangei(12) SMx->a[i] *= -1;
  rangei( 3) SMx->s.T[i] += InvT[i];
}


int CmpiVect(const int *a, const int *b, int n)
{
  int  n0a, n0b, i;


  n0a = 0; rangei(n) if (a[i] == 0) n0a++;
  n0b = 0; rangei(n) if (b[i] == 0) n0b++;

  if (n0a > n0b) return -1;
  if (n0a < n0b) return  1;

  rangei(n)
  {
    if (a[i] != 0 && b[i] == 0) return -1;
    if (a[i] == 0 && b[i] != 0) return  1;
  }

  rangei(n)
  {
    if (abs(a[i]) < abs(b[i])) return -1;
    if (abs(a[i]) > abs(b[i])) return  1;
  }

  rangei(n)
  {
    if (a[i] > b[i]) return -1;
    if (a[i] < b[i]) return  1;
  }

  return 0;
}


static int CmpLTr(const T_LTr *a, const T_LTr *b)
{
  return MemCmp(a, b, 1);
}


static int CmpSMx(const T_RTMx *a, const T_RTMx *b)
{
  int     i;
  T_RMxI  RI_a[1], RI_b[1];


  SetRotMxInfo(a->s.R, RI_a);
  SetRotMxInfo(b->s.R, RI_b);

  if (abs(RI_a->Rtype) > abs(RI_b->Rtype)) return -1;
  if (abs(RI_a->Rtype) < abs(RI_b->Rtype)) return  1;

  if (RI_a->Rtype > RI_b->Rtype) return -1;
  if (RI_a->Rtype < RI_b->Rtype) return  1;

      i = CmpiVect(RI_a->EV, RI_b->EV, 3);
  if (i) return i;

  if (RI_a->SenseOfRotation > RI_b->SenseOfRotation) return -1;
  if (RI_a->SenseOfRotation < RI_b->SenseOfRotation) return  1;

      i = CmpiVect(a->s.T, b->s.T, 3);
  if (i) return i;

  return MemCmp(a, b, 1);
}


static int TidyT(int nLTr, const T_LTr *LTr, int LTBF, int T[3], int TBF)
{
  int  iLTr, i;
  int  BestT[3], TrialT[3];
  int  LCMTBF, fTBF, fLTBF;


  LCMTBF = iLCM(LTBF, TBF);
   fTBF = LCMTBF /  TBF;
  fLTBF = LCMTBF / LTBF;

  rangei(3) BestT[i] = T[i] * fTBF;
  ViModShort(BestT, 3, LCMTBF);

  range2(iLTr, 1, nLTr) {
    rangei(3)
      TrialT[i] = iModShort(T[i] * fTBF + LTr[iLTr].v[i] * fLTBF, LCMTBF);
    if (CmpiVect(BestT, TrialT, 3) > 0) MemCpy(BestT, TrialT, 3);
  }

  if (ChangeBaseFactor(BestT, LCMTBF, T, TBF, 3) != 0) return -1;
  ViModPositive(T, 3, TBF);

  return 0;
}


int TidySgOps(T_SgOps *SgOps)
{
  int  iLTr, iSMx;
  int  Rtype;


  if (SgOps->fInv == 2)
  {
    if (TidyT(SgOps->nLTr, SgOps->LTr, STBF, SgOps->InvT, STBF) != 0)
      return IE(-1);

    range2(iSMx, 1, SgOps->nSMx) {
          Rtype = GetRtype(SgOps->SMx[iSMx].s.R);
      if (Rtype == 0) return IE(-1);
      if (Rtype < 0) InvT_SMx(SgOps->InvT, &SgOps->SMx[iSMx]);
    }
  }

  range2(iSMx, 1, SgOps->nSMx)
    if (TidyT(SgOps->nLTr, SgOps->LTr, STBF, SgOps->SMx[iSMx].s.T, STBF) != 0)
      return IE(-1);

  if (SgOps->nLTr > 2)
    qsort((void *) &SgOps->LTr[1], SgOps->nLTr - 1, sizeof (*SgOps->LTr),
          (int (*)(const void *, const void *)) CmpLTr);

  range2(iLTr, SgOps->nLTr, SgOps_mLTr)
    IntSetZero(SgOps->LTr[iLTr].v, 3);

  if (SgOps->nSMx > 2)
    qsort((void *) &SgOps->SMx[1], SgOps->nSMx - 1, sizeof (*SgOps->SMx),
          (int (*)(const void *, const void *)) CmpSMx);

  range2(iSMx, SgOps->nSMx, SgOps_mSMx)
    InitRTMx(&SgOps->SMx[iSMx], -1);

  return 0;
}
