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


void ResetLLTr(T_LTr *LLTr, int *nLLTr)
{
  int  i;

  rangei(3) LLTr[0].v[i] = 0;
  *nLLTr = 1;
}


static int AddLLTr(int LTBF, int mLLTr,
                   T_LTr *LLTr, int *nLLTr, const int *NewLTr)
{
  int    NLTr[3], iLLTr, i;


  rangei(3) NLTr[i] = iModPositive(NewLTr[i], LTBF);

  for (iLLTr = 0; iLLTr < *nLLTr; iLLTr++, LLTr++)
    if (MemCmp(LLTr->v, NLTr, 3) == 0)
      return 0;

  if (*nLLTr >= mLLTr) return -1;

  MemCpy(LLTr->v, NLTr, 3);
  (*nLLTr)++;

  return 1;
}


int ExpLLTr(int LTBF, int mLLTr,
            T_LTr *LLTr, int *nLLTr, const int *NewLTr)
{
  int           TrialLTr[3], i;
  int           iLLTr,  jLLTr;
  const T_LTr   *LLTri, *LLTrj;


  iLLTr  = *nLLTr;
   LLTri =  &LLTr[iLLTr];

  jLLTr  = 1;
   LLTrj =  &LLTr[jLLTr];

  for (;;)
  {
    if (NewLTr) {
      if (AddLLTr(LTBF, mLLTr, LLTr, nLLTr, NewLTr) < 0)
        return -1;
    }

    if (jLLTr > iLLTr)
    {
      iLLTr++;
       LLTri++;

      jLLTr  = 1;
       LLTrj = &LLTr[jLLTr];
    }

    if (iLLTr == *nLLTr)
      break;

    rangei(3) TrialLTr[i] = LLTrj->v[i] + LLTri->v[i];
    NewLTr = TrialLTr;

    jLLTr++;
     LLTrj++;
  }

  return 0;
}


void ResetSgOps(T_SgOps *SgOps)
{
  SgOps->NoExpand = 0;
  SgOps->nLSL = 1;
  SgOps->nSSL = 1;
  MemSet(SgOps->LTr, SgOps_mLTr, 0);
  ResetLLTr(SgOps->LTr, &SgOps->nLTr);
  SgOps->fInv = 1;
  IntSetZero(SgOps->InvT, 3);
  SgOps->nSMx = 1;
  MemSet(SgOps->SMx, SgOps_mSMx, 0);
  InitRTMx(&SgOps->SMx[0], 1);
}


static int AddSgLTr(T_SgOps *SgOps, const int *NewLTr)
{
  int  Stat;

      Stat = AddLLTr(STBF, SgOps_mLTr, SgOps->LTr, &SgOps->nLTr, NewLTr);
  if (Stat < 0) {
    SetSgError("Internal Error: SgOps_mLTr too small");
    return -1;
  }

  return Stat;
}


static int AddLtrDueToInvT(T_SgOps *SgOps, const T_RTMx *LSMx)
{
  int  NewLTr[3], i;


  RotMx_t_Vector(NewLTr, LSMx->s.R, SgOps->InvT, 0);
  for (i = 0; i < 3; i++) NewLTr[i] += 2 * LSMx->s.T[i] - SgOps->InvT[i];
  return AddSgLTr(SgOps, NewLTr);
}


static int AddSgInv(T_SgOps *SgOps, const int *InvT)
{
  int     NewLTr[3], i;
  int     iLSMx;
  T_RTMx  *LSMx;

  const int  NNN[] = { 0, 0, 0};


  if (! InvT) InvT = NNN;

  if (SgOps->fInv == 2) /* there is a centre of inversion already */
  {
    for (i = 0; i < 3; i++) NewLTr[i] = SgOps->InvT[i] - InvT[i];
    return AddSgLTr(SgOps, NewLTr);
  }

  for (i = 0; i < 3; i++) SgOps->InvT[i] = iModPositive(InvT[i], STBF);

  SgOps->fInv = 2;

  if (! SgOps->NoExpand)
  {
    LSMx = &SgOps->SMx[1];

    for (iLSMx = 1; iLSMx < SgOps->nSMx; iLSMx++, LSMx++)
      if (AddLtrDueToInvT(SgOps, LSMx) < 0)
        return -1;
  }

  return 1;
}


static int AddSgSMx(T_SgOps *SgOps, const T_RTMx *NewSMx)
{
  int     mR[9], NewLTr[3], InvT[3], i;
  int     iLSMx;
  T_RTMx  *LSMx;


  for (i = 0; i < 9; i++) mR[i] = -NewSMx->s.R[i];

  LSMx = SgOps->SMx;

  for (iLSMx = 0; iLSMx < SgOps->nSMx; iLSMx++, LSMx++)
  {
    if (MemCmp(LSMx->s.R, NewSMx->s.R, 9) == 0) {
      for (i = 0; i < 3; i++)   NewLTr[i] = LSMx->s.T[i] - NewSMx->s.T[i];
      return AddSgLTr(SgOps, NewLTr);
    }

    if (MemCmp(LSMx->s.R, mR, 9) == 0) {
      for (i = 0; i < 3; i++)     InvT[i] = LSMx->s.T[i] + NewSMx->s.T[i];
      return AddSgInv(SgOps, InvT);
    }
  }

  if (SgOps->nSMx >= SgOps_mSMx) {
    SetSgError("Error: Non-crystallographic rotation matrix encountered");
    return -1;
  }

  MemCpy(LSMx->s.R, NewSMx->s.R, 9);
  for (i = 0; i < 3; i++) LSMx->s.T[i] = iModPositive(NewSMx->s.T[i], STBF);

  SgOps->nSMx++;

  if (   ! SgOps->NoExpand
      &&   SgOps->fInv == 2
      && AddLtrDueToInvT(SgOps, LSMx) < 0)
    return -1;

  return 1;
}


static int DoMulSMxLTr(T_SgOps *SgOps, int iLSMx, int iLLTr, int OldOnly)
{
  const T_RTMx  *LSMxi;
  const T_LTr   *LLTri;
  int           NewLTr[3];
  const int     iLLTr0 = iLLTr;


  LSMxi = &SgOps->SMx[iLSMx];

  for (; iLSMx < SgOps->nSMx; iLSMx++, LSMxi++)
  {
    LLTri = &SgOps->LTr[iLLTr0];

    for (iLLTr = iLLTr0; iLLTr < (OldOnly ? SgOps->nLSL : SgOps->nLTr);
         iLLTr++, LLTri++)
    {
      RotMx_t_Vector(NewLTr, LSMxi->s.R, LLTri->v, 0);
      if (AddSgLTr(SgOps, NewLTr) < 0)
        return -1;
    }
  }

  return 0;
}


int ExpSgLTr(T_SgOps *SgOps, const int *NewLTr)
{
  int           TrialLTr[3], i;
  int           iLLTr,  jLLTr;
  const T_LTr   *LLTri, *LLTrj;


  if (SgOps->NoExpand) {
    if (NewLTr) return AddSgLTr(SgOps, NewLTr);
    return 0;
  }

  if (DoMulSMxLTr(SgOps, SgOps->nSSL, 1, 1) < 0) return -1;
  SgOps->nSSL = SgOps->nSMx;

  iLLTr  =  SgOps->nLSL;
   LLTri = &SgOps->LTr[iLLTr];

  jLLTr  = 1;
   LLTrj = &SgOps->LTr[jLLTr];

  for (;;)
  {
    if (NewLTr) {
      if (AddSgLTr(SgOps, NewLTr) < 0)
        return -1;
    }

    if (DoMulSMxLTr(SgOps, 1, SgOps->nLSL, 0) < 0) return -1;
    SgOps->nLSL = SgOps->nLTr;

    if (jLLTr > iLLTr)
    {
      iLLTr++;
       LLTri++;

      jLLTr  = 1;
       LLTrj = &SgOps->LTr[jLLTr];
    }

    if (iLLTr == SgOps->nLTr)
      break;

    for (i = 0; i < 3; i++)
      TrialLTr[i] = LLTrj->v[i] + LLTri->v[i];

    NewLTr = TrialLTr;

    jLLTr++;
     LLTrj++;
  }

  return 0;
}


int ExpSgInv(T_SgOps *SgOps, const int *InvT)
{
  if (AddSgInv(SgOps, InvT) < 0)
    return -1;

  return ExpSgLTr(SgOps, NULL);
}


int ExpSgSMx(T_SgOps *SgOps, const T_RTMx *NewSMx)
{
  int           iLSMx,  jLSMx;
  const T_RTMx  *LSMxi, *LSMxj;
  T_RTMx        TrialSMx[1];


  if (SgOps->NoExpand) {
    if (NewSMx) return AddSgSMx(SgOps, NewSMx);
    return 0;
  }

  iLSMx  =  SgOps->nSMx;
   LSMxi = &SgOps->SMx[iLSMx];

  jLSMx  = 1;
   LSMxj = &SgOps->SMx[jLSMx];

  for (;;)
  {
    if (NewSMx && AddSgSMx(SgOps, NewSMx) < 0)
      return -1;

    if (jLSMx > iLSMx)
    {
      iLSMx++;
       LSMxi++;

      jLSMx  = 1;
       LSMxj = &SgOps->SMx[jLSMx];
    }

    if (iLSMx == SgOps->nSMx)
      break;

    SeitzMxMultiply(TrialSMx, LSMxj, LSMxi);

    NewSMx = TrialSMx;

    jLSMx++;
     LSMxj++;
  }

  return ExpSgLTr(SgOps, NULL);
}


int ExpSgRMx(T_SgOps *SgOps, const int NewRMx[9])
{
  T_RTMx  SMx[1];

  MemCpy(SMx->s.R, NewRMx, 9);
  IntSetZero(SMx->s.T, 3);
  return ExpSgSMx(SgOps, SMx);
}
