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
#include <ctype.h>


#undef SG_GLOBAL
#include "sglite.h"


#define T(i) ((i) * (STBF / 12))
#define V(i, j, k) {{ T(i), T(j), T(k) }}

static const T_LTr  LTr_P[] = { V(0, 0, 0)
                              };
static const T_LTr  LTr_A[] = { V(0, 0, 0),
                                V(0, 6, 6)
                              };
static const T_LTr  LTr_B[] = { V(0, 0, 0),
                                V(6, 0, 6)
                              };
static const T_LTr  LTr_C[] = { V(0, 0, 0),
                                V(6, 6, 0)
                              };
static const T_LTr  LTr_I[] = { V(0, 0, 0),
                                V(6, 6, 6)
                              };
static const T_LTr  LTr_R[] = { V(0, 0, 0),
                                V(8, 4, 4),
                                V(4, 8, 8)
                              };
static const T_LTr  LTr_Q[] = { V(0, 0, 0), /* reverse setting,      */
                                V(4, 8, 4), /* for internal use only */
                                V(8, 4, 8)
                              };
static const T_LTr  LTr_H[] = { V(0, 0, 0),
                                V(8, 4, 0),
                                V(4, 8, 0)
                              };
static const T_LTr  LTr_F[] = { V(0, 0, 0),
                                V(0, 6, 6),
                                V(6, 0, 6),
                                V(6, 6, 0)
                              };
#undef  T
#undef  V


static const T_ConvCType  TabConvCType[] =
  {
    { 'P', 1, LTr_P },
    { 'A', 2, LTr_A },
    { 'B', 2, LTr_B },
    { 'C', 2, LTr_C },
    { 'I', 2, LTr_I },
    { 'R', 3, LTr_R },
    { 'Q', 3, LTr_Q }, /* reverse setting, for internal use only */
    { 'H', 3, LTr_H },
    { 'F', 4, LTr_F }
  };

#define nTabConvCType (sizeof TabConvCType / sizeof (*TabConvCType))


static const int CCMx_PP[] = {  1,  0,  0,  /* Change of Basis Matrices     */
                                0,  1,  0,  /* (coordinate transformations) */
                                0,  0,  1
                             };
static const int CCMx_AP[] = { -1,  0,  0,
                                0, -1,  1,
                                0,  1,  1
                             };
static const int CCMx_BP[] = { -1,  0,  1,
                                0, -1,  0,
                                1,  0,  1
                             };
static const int CCMx_CP[] = {  1,  1,  0,
                                1, -1,  0,
                                0,  0, -1
                             };
static const int CCMx_IP[] = {  0,  1,  1,
                                1,  0,  1,
                                1,  1,  0
                             };
static const int CCMx_RP[] = {  1,  0,  1,
                               -1,  1,  1,
                                0, -1,  1
                             };
static const int CCMx_HP[] = {  1,  1,  0,
                               -1,  2,  0,
                                0,  0,  1
                             };
static const int CCMx_FP[] = { -1,  1,  1,
                                1, -1,  1,
                                1,  1, -1
                             };


static const T_ConvCType *GetConvCType(int Sym)
{
  int                iCType;
  const T_ConvCType  *CType;

  Sym = toupper(Sym);

  if (Sym == 'Q') return NULL;

  CType = TabConvCType;

  for (iCType = 0; iCType < nTabConvCType; iCType++, CType++)
    if (CType->Sym == Sym) return CType;

  return NULL;
}


static int ExpSgConvCType(T_SgOps *SgOps, const T_ConvCType *CType)
{
  int  RetVal, iLTr, RV;

  RetVal = 0;

  for (iLTr = 0; iLTr < CType->nLTr; iLTr++)
  {
        RV = ExpSgLTr(SgOps, CType->LTr[iLTr].v);
    if (RV < 0)
      return -1;
    if (RV != 0)
      RetVal++;
  }

  return RetVal;
}


int ExpSgSymCType(T_SgOps *SgOps, int Sym)
{
  const T_ConvCType  *CType;

      CType = GetConvCType(Sym);
  if (CType == NULL) {
    SetSgError("Error: Illegal symbol for centring type of cell");
    return -1;
  }

  return ExpSgConvCType(SgOps, CType);
}


int GetSymCType(int nLTr, const T_LTr *LTr)
{
  int                iLTr, jLTr;
  int                nMatch, Match[4];
  int                iCType;
  const T_ConvCType  *CType;

  CType = TabConvCType;

  for (iCType = 0; iCType < nTabConvCType; iCType++, CType++)
  {
    if (CType->nLTr == nLTr)
    {
      nMatch = 0;
      range1(iLTr, nLTr) Match[iLTr] = 0;

      range1(iLTr, nLTr) {
        range1(jLTr, nLTr) {
          if (   Match[jLTr] == 0
              && MemCmp(CType->LTr[iLTr].v, LTr[jLTr].v, 3) == 0) {
            Match[jLTr] = 1;
            nMatch++;
            break;
          }
        }
      }

      if (nMatch == nLTr)
        return CType->Sym;
    }
  }

  return '\0';
}


static const int *GetCCMxSymCTypeToPrimitive(int SymCType)
{
  switch (SymCType) {
    case 'P': return CCMx_PP;
    case 'A': return CCMx_AP;
    case 'B': return CCMx_BP;
    case 'C': return CCMx_CP;
    case 'I': return CCMx_IP;
    case 'R': return CCMx_RP;
    case 'H': return CCMx_HP;
    case 'F': return CCMx_FP;
    default: break;
  }

  return NULL;
}


static int GetStdZ2PCBMx(int nLTr, const T_LTr *LTr, T_RTMx Z2PCBMx[2])
{
  int        SymCType, DetF, i;
  const int  *CCMx;

  SymCType = GetSymCType(nLTr, LTr);

      CCMx = GetCCMxSymCTypeToPrimitive(SymCType);
  if (CCMx == NULL) return 0;

  rangei(9) Z2PCBMx[0].s.R[i] = CRBF * CCMx[i];

      DetF = InverseRotMx(Z2PCBMx[0].s.R, Z2PCBMx[1].s.R, CRBF);
  if (DetF != nLTr * (CRBF * CRBF * CRBF)) return IE(-1);

  rangei(3) Z2PCBMx[0].s.T[i] = 0;
  rangei(3) Z2PCBMx[1].s.T[i] = 0;

  return 1;
}


static int CmpTLT(const T_LTr *a, const T_LTr *b)
{
  return CmpiVect(a->v, b->v, 3);
}


static int BuildListTotLTr(const int nLTr, const T_LTr *LTr,
                           const int mTLT,       T_LTr *TLT)
{
  int  nTLT, iTLT, iLTr, iUTr, LinDep, i;
  int  UnitTr[3], nUTr[3], V[3];

  nTLT = 0;

  range2(iLTr, 1, nLTr)
  {
    rangei(3)                     nUTr[i] = 1;
    rangei(3) if (LTr[iLTr].v[i]) nUTr[i] = 2;

    range1(UnitTr[0], nUTr[0])
    range1(UnitTr[1], nUTr[1])
    range1(UnitTr[2], nUTr[2])
    {
      rangei(3) {
            V[i] = LTr[iLTr].v[i] - UnitTr[i] * STBF;
            V[i] *= CRBF;
        if (V[i] %  STBF) return IE(-1);
            V[i] /= STBF;
      }

      range1(iTLT, nTLT) {
              LinDep = AreLinDepV(TLT[iTLT].v, V);
        if   (LinDep) {
          if (LinDep > 0) (void) MemCpy(TLT[iTLT].v, V, 3);
          break;
        }
      }

      if (iTLT == nTLT) {
        if (iTLT == mTLT) return IE(-1);
        (void) MemCpy(TLT[iTLT].v, V, 3);
        nTLT++;
      }
    }
  }

  qsort((void *) TLT, nTLT, sizeof (*TLT),
        (int (*)(const void *, const void *)) CmpTLT);

  if (nTLT + 3 > mTLT) return IE(-1);

  range1(iUTr, 3) {
    rangei(3) TLT[nTLT].v[i] = (i == iUTr ? CRBF : 0);
    nTLT++;
  }

  return nTLT;
}


static int IsLTrBasis(const int nLTr, const T_LTr *LTr,
                      int Basis[2][9])
{
  int  DetF, iLTr, V[3], i;

      DetF = deterRotMx(Basis[0]);
  if (DetF == 0) return 0;
  if (DetF < 0) {
      DetF *= -1;
    rangei(3) Basis[0][i * 3] *= -1;
  }

  if (DetF * nLTr != CRBF * CRBF * CRBF) return 0;

  iCoFactorMxTp(Basis[0], Basis[1]);

  rangei(9) Basis[1][i] *= (CRBF * CRBF);

  rangei(9) {
    if (Basis[1][i] %  DetF) return 0;
        Basis[1][i] /= DetF;
  }

  range2(iLTr, 1, nLTr)
  {
    RotMx_t_Vector(V, Basis[1], LTr[iLTr].v, 0);

    rangei(3)
      if (V[i] % (CRBF * STBF))
        return 0;
  }

  return 1;
}


static int CheckLTrBasis(const T_SgOps *SgOps, const int Basis[2][9],
                         T_RTMx CBMx[2])
{
  int      iMx, i;
  T_SgOps  TstSgOps[1];

  range1(iMx, 2) {
    rangei(9) CBMx[1 - iMx].s.R[i] = Basis[iMx][i];
    rangei(3) CBMx[1 - iMx].s.T[i] = 0;
  }

  ResetSgOps(TstSgOps);
  if (CB_SgOps(SgOps, &CBMx[0], &CBMx[1], TstSgOps) != 0) {
    ClrSgError();
    return 0;
  }

  return 1;
}


static int ConstructZ2PCBMx(const T_SgOps *SgOps, T_RTMx Z2PCBMx[2])
{
  int  stat, i;
  int  Basis[2][9];

#define mTLT 320
  T_LTr  TLT[mTLT];
  int   nTLT, iTLT[3];

      nTLT = BuildListTotLTr(SgOps->nLTr, SgOps->LTr, mTLT, TLT);
  if (nTLT < 0)
    return IE(-1);
#undef mTLT

  range2    (iTLT[0],           0, nTLT - 2)
  {
        rangei(3) Basis[0][i * 3 + 0] = TLT[iTLT[0]].v[i];

    range2  (iTLT[1], iTLT[0] + 1, nTLT - 1)
    {
        rangei(3) Basis[0][i * 3 + 1] = TLT[iTLT[1]].v[i];

      range2(iTLT[2], iTLT[1] + 1, nTLT    )
      {
        rangei(3) Basis[0][i * 3 + 2] = TLT[iTLT[2]].v[i];

        if (IsLTrBasis(SgOps->nLTr, SgOps->LTr, Basis)) {
              stat = CheckLTrBasis(SgOps, (const int (*)[9]) Basis, Z2PCBMx);
          if (stat < 0) return IE(-1);
          if (stat) return 0;
        }
      }
    }
  }

  return IE(-1);
}


int GetZ2PCBMx(const T_SgOps *SgOps, T_RTMx Z2PCBMx[2])
{
  int  stat;

      stat = GetStdZ2PCBMx(SgOps->nLTr, SgOps->LTr, Z2PCBMx);
  if (stat < 0) return IE(-1);
  if (stat == 0) {
    if (ConstructZ2PCBMx(SgOps, Z2PCBMx) != 0) return IE(-1);
  }

  return 0;
}
