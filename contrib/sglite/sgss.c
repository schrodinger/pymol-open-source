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
#include "sgconst.h"


#define  mDiscrGr  8


typedef struct
  {
    int  P[3];
    int  Z[3];
  }
  T_DiscrGr;


int Is_ss(const T_ssVM *ssVM, int n_ssVM, int h, int k, int l)
{
  int  i_ssVM, u;

  range1(i_ssVM, n_ssVM)
  {
    u =  ssVM[i_ssVM].V[0] * h;
    u += ssVM[i_ssVM].V[1] * k;
    u += ssVM[i_ssVM].V[2] * l;

    if (ssVM[i_ssVM].M) {
      if (u % ssVM[i_ssVM].M) return 0; }
    else {
      if (u)                  return 0; }
  }

  return 1;
}


void Set_uvw(const T_ssVM *ssVM, int n_ssVM, int h, int k, int l, int uvw[3])
{
  int  i_ssVM, u;

  range1(i_ssVM, n_ssVM)
  {
    u =  ssVM[i_ssVM].V[0] * h;
    u += ssVM[i_ssVM].V[1] * k;
    u += ssVM[i_ssVM].V[2] * l;

    if (ssVM[i_ssVM].M) u %= ssVM[i_ssVM].M;

    uvw[i_ssVM] = u;
  }
}


/* sgss.c:SetAnyIxGen() and sginfo2/metric.c:SetIxGen() should be consolidated.
 */
static int SetAnyIxGen(const T_SgOps *SgOps, int PG_MGC, int IxGen[2])
{
  int     nGen, iLSMx, i;
  int     PrincipalProperOrder;
  T_RMxI  PrincipalRI[1], RI[1];

  rangei(2) IxGen[i] = -1;

  if (PG_MGC == MGC_Undefined) PG_MGC = GetPG(SgOps);
  if (PG_MGC == MGC_Unknown) return IE(-1);

  nGen = 0;
  PrincipalProperOrder = 0;

  switch (ixXS(PG_MGC))
  {
    case XS_Triclinic:
      if (SgOps->fInv == 1) {
        IxGen[0] = 0;
        nGen = 1;
      }
      break;

    case XS_Monoclinic:
      IxGen[0] = 1;
      nGen = 1;
      break;

    case XS_Orthorhombic:
      IxGen[0] = 1;
      IxGen[1] = 2;
      nGen = 2;
      break;

    case XS_Tetragonal:
                                  PrincipalProperOrder = 4;
    case XS_Trigonal:
      if (! PrincipalProperOrder) PrincipalProperOrder = 3;
    case XS_Hexagonal:
      if (! PrincipalProperOrder) PrincipalProperOrder = 6;

      range2(iLSMx, 1, SgOps->nSMx)
      {
        if (SetRotMxInfo(SgOps->SMx[iLSMx].s.R, PrincipalRI) == 0)
          return IE(-1);

        if (abs(PrincipalRI->Rtype) == PrincipalProperOrder) {
          if (PrincipalRI->SenseOfRotation > 0) {
            IxGen[0] = iLSMx;
            nGen++;
            break;
          }
        }
      }

      if (nGen < 1) return IE(-1);

      range2(iLSMx, 1, SgOps->nSMx)
      {
        if (iLSMx == IxGen[0]) continue;

        if (SetRotMxInfo(SgOps->SMx[iLSMx].s.R, RI) == 0)
          return IE(-1);

        if (   abs(RI->Rtype) == 2
            && MemCmp(RI->EV, PrincipalRI->EV, 3) != 0) {
          IxGen[1] = iLSMx;
          nGen++;
          break;
        }
      }

      break;

    case XS_Cubic:
      range2(iLSMx, 1, SgOps->nSMx)
      {
        if (SetRotMxInfo(SgOps->SMx[iLSMx].s.R, RI) == 0)
          return IE(-1);

        if      (   abs(RI->Rtype) == 3
                 && RI->SenseOfRotation > 0
                 && IxGen[0] < 0) {
                    IxGen[0] = iLSMx;
          nGen++;
        }
        else if (   abs(RI->Rtype) == SgOps->nSMx / 6
                 && RI->SenseOfRotation >= 0
                 && IxGen[1] < 0) {
                    IxGen[1] = iLSMx;
          nGen++;
        }
      }

      if (nGen != 2) return IE(-1);

      break;

    default:
      return IE(-1);
  }

  return nGen;
}


static int CmpOLen2(const int a[3], const int b[3])
{
  int OLen2a, OLen2b, i;

  OLen2a = 0; rangei(3) OLen2a += a[i] * a[i];
  OLen2b = 0; rangei(3) OLen2b += b[i] * b[i];

  if (OLen2a < OLen2b) return -1;
  if (OLen2a > OLen2b) return  1;

  return CmpiVect(a, b, 3);
}


static int ConstructGenRmI(const T_SgOps *SgOps, const T_RTMx Z2PCBMx[2],
                           const int IxGen[2], int nGen,
                           int GenRmI[3 * 3 * 3])
{
  int  nrGenRmI;
  int  iRmI, iGen, i;

      nrGenRmI = (SgOps->fInv - 1 + nGen) * 3;
  if (nrGenRmI > 9) return IE(-1);

  iRmI = 0;

  if (SgOps->fInv == 2) {
    SetRminusI(SgOps->SMx[0].s.R, GenRmI, 1);
    iRmI++;
  }

  if (Z2PCBMx == NULL) {
    range1(iGen, nGen) {
      SetRminusI(SgOps->SMx[IxGen[iGen]].s.R, &GenRmI[iRmI * 9], 0);
      iRmI++;
    }
  }
  else {
    range1(iGen, nGen) {
      if (CB_RMx(&GenRmI[iRmI * 9],
                 Z2PCBMx[0].s.R,
                 SgOps->SMx[IxGen[iGen]].s.R,
                 Z2PCBMx[1].s.R) != 0) return -1;
      range3(i, 0, 9, 4) GenRmI[iRmI * 9 + i] -= 1;
      iRmI++;
    }
  }

  if (iRmI * 3 != nrGenRmI) return IE(-1);

  return nrGenRmI;
}


static int GetContNullSpace(const T_SgOps *SgOps, int IxGen[2], int nGen,
                            T_ssVM ssVM[3])
{
  int  GenRmI[3 * 3 * 3];
  int  nrGenRmI, RankGenRmI;
  int  nIndep, iIndep, IxIndep[3], Sol[4][3];
  int  n_ssVM;

      nrGenRmI = ConstructGenRmI(SgOps, NULL, IxGen, nGen, GenRmI);
  if (nrGenRmI < 0) return IE(-1);

      RankGenRmI = iRowEchelonFormT(GenRmI, nrGenRmI, 3, NULL, 0);
  if (RankGenRmI < 0 || RankGenRmI > 3)
    return IE(-1);

  n_ssVM = 3 - RankGenRmI;

      nIndep = iRESetIxIndep(GenRmI, RankGenRmI, 3, IxIndep, 3);
  if (nIndep < 0) return IE(-1);

  if (nIndep != 2)
  {
    range1(iIndep, nIndep)
    {
      ssVM[iIndep].V[IxIndep[iIndep]] = 1;

      if (iREBacksubst(GenRmI, NULL, RankGenRmI, 3,
                       ssVM[iIndep].V, NULL) < 1)
        return IE(-1);

      ssVM[iIndep].M = 0;
    }
  }
  else
  {
    if (SolveHomRE1(GenRmI, IxIndep, Sol) != 0) return -1;

    qsort((void *) Sol, 4, sizeof (*Sol),
          (int (*)(const void *, const void *)) CmpOLen2);

    range1(iIndep, 2) {
      MemCpy(ssVM[iIndep].V, Sol[iIndep], 3);
      ssVM[iIndep].M = 0;
    }
  }

  return n_ssVM;
}


static int nDLoopStep(int *i, int n, int Low, int High)
{
  int  p, l;

  p = l = n - 1;

  for (; p >= 0;)
  {
        i[p]++;
    if (i[p] > High)
      p--;
    else if (p < l)
      i[++p] = Low - 1;
    else
      return 1;
  }

  return 0;
}


static void UpdateBestZ(int OrigZf[mDiscrGr][3], int nDiscrGr,
                        int BestZf[mDiscrGr][3], int BestM[mDiscrGr],
                        int BestZc[mDiscrGr][3],
                        int Shift[3], int LTBF)
{
  int  iDG, c, i;
  int  Zf[3], M;
  int  Zc[3];

  range2(iDG, 1, nDiscrGr)
  {
    rangei(3) Zf[i] = iModPositive(Shift[i] + OrigZf[iDG][i], LTBF);
    MemCpy(Zc, Zf, 3);
    M = CancelBFGCD(Zc, 3, LTBF);

    rangei(3) {
      if (Zf[i]) {
        c = CmpOLen2(BestZc[iDG], Zc);
        if (c > 0 || (c == 0 && BestM[iDG] > M)) {
          MemCpy(BestZf[iDG], Zf, 3);
          MemCpy(BestZc[iDG], Zc, 3);
          BestM[iDG] = M;
        }
        break;
      }
    }
  }
}


static int BestVect(const T_SgOps *SgOps,
                    const T_ssVM *ssVM, int n_ssVM,
                    int DTBF, T_DiscrGr DiscrGr[mDiscrGr], int nDiscrGr)
{
  int  iLTr, i_ssVM, iDG, gcd, i, j;
  int  fGrd, LTBF, f[2];
  int  OrigZf[mDiscrGr][3];
  int  BestZf[mDiscrGr][3], BestM[mDiscrGr];
  int  BestZc[mDiscrGr][3];
  int  LTr[3][3];

  fGrd = 1;
  LTBF = 1;

  range2(iDG, 1, nDiscrGr) {
    rangei(3) {
      gcd = iGCD(DiscrGr[iDG].Z[i], DTBF * CRBF);
      LTBF = iLCM(LTBF, (DTBF * CRBF) / gcd);
    }
  }

  range2(iLTr, 1, SgOps->nLTr) {
    rangei(3) {
      gcd = iGCD(SgOps->LTr[iLTr].v[i], STBF);
      LTBF = iLCM(LTBF, STBF / gcd);
    }
  }

  range1(i_ssVM, n_ssVM)
    rangei(3) fGrd = iLCM(fGrd, ssVM[i_ssVM].V[i]);

  LTBF *= fGrd;

  if (LTBF > 6) LTBF = iLCM(LTBF,  6);
  else          LTBF = iLCM(LTBF, 12);

  if (SgOps->nLTr == 1 && n_ssVM == 0) return 0;

  range2(iDG, 1, nDiscrGr) {
    if (ChangeBaseFactor(DiscrGr[iDG].Z, DTBF * CRBF,
                            OrigZf[iDG],        LTBF, 3) != 0)
      return IE(-1);
    rangei(3) OrigZf[iDG][i] = iModPositive(OrigZf[iDG][i], LTBF);
    MemCpy(BestZf[iDG], OrigZf[iDG], 3);
    MemCpy(BestZc[iDG], OrigZf[iDG], 3);
    BestM[iDG] = CancelBFGCD(BestZc[iDG], 3, LTBF);
  }

  if (n_ssVM > 2) return IE(-1);

  range1(iLTr, SgOps->nLTr)
  {
    if (ChangeBaseFactor(SgOps->LTr[iLTr].v, STBF, LTr[0], LTBF, 3) != 0)
      return IE(-1);

    rangei(n_ssVM) f[i] = 0;

    do
    {
      rangei(n_ssVM)
        range1(j, 3) LTr[i + 1][j] = LTr[i][j] + f[i] * ssVM[i].V[j];

      UpdateBestZ(OrigZf, nDiscrGr, BestZf, BestM, BestZc,
                  LTr[n_ssVM], LTBF);
    }
    while (nDLoopStep(f, n_ssVM, 0, LTBF - 1));
  }

  range2(iDG, 1, nDiscrGr)
    if (ChangeBaseFactor(BestZf[iDG],           LTBF,
                         DiscrGr[iDG].Z, DTBF * CRBF, 3) != 0)
      return IE(-1);

  return 0;
}


static int SelectDiscrete(int LTBF, int nDiscrGr, T_DiscrGr DiscrGr[mDiscrGr],
                          int mIx, int Ix[3])
{
  int    nIx, iIx;
  T_LTr  LLTr[mDiscrGr];
  int    nLLTr;

  if (nDiscrGr == 1) return 0;

  for (nIx = 1; nIx <= nDiscrGr - 1 && nIx <= mIx; nIx++)
  {
    range1(iIx, nIx) Ix[iIx] = iIx;

    do
    {
      ResetLLTr(LLTr, &nLLTr);
      range1(iIx, nIx)
        if (ExpLLTr(LTBF, mDiscrGr,
                    LLTr, &nLLTr, DiscrGr[Ix[iIx] + 1].P) < 0)
          return IE(-1);

      if (nLLTr >  nDiscrGr) return IE(-1);
      if (nLLTr == nDiscrGr) return nIx;
    }
    while (NextOf_n_from_m(nDiscrGr - 1, nIx, Ix) != 0);
  }

  return IE(-1);
}


static int CmpDiscr(const T_DiscrGr *a, const T_DiscrGr *b)
{
  return CmpiVect(a->Z, b->Z, 3);
}


static int Cmp_ssVM(const T_ssVM *a, const T_ssVM *b)
{
  return CmpiVect(a->V, b->V, 3);
}


int Set_ss(const T_SgOps *SgOps, T_ssVM ssVM[3])
{
  int        ir, i;
  int        nGen, IxGen[2];
  int        n_ssVM;
  T_RTMx     Z2PCBMx[2];
  int        nrSNF, DTBF, nd, id, d, f;
  int        SNF[3 * 3 * 3], Q[3 * 3], xp[3], x[3];
  int        nDiscrGr, iDG;
  T_LTr      LLTr[mDiscrGr];
  T_DiscrGr  DiscrGr[mDiscrGr];
  int        nIx, Ix[3];

  range1(ir, 3) rangei(3) ssVM[ir].V[i] = 0;
  range1(ir, 3) ssVM[ir].M = -1;

      nGen = SetAnyIxGen(SgOps, MGC_Undefined, IxGen);
  if (nGen < 0 || nGen > 2) return IE(-1);

      n_ssVM = GetContNullSpace(SgOps, IxGen, nGen, ssVM);
  if (n_ssVM < 0) return -1;
  if (n_ssVM == 3) return n_ssVM;

  if (GetZ2PCBMx(SgOps, Z2PCBMx) != 0) return -1;

      nrSNF = ConstructGenRmI(SgOps, Z2PCBMx, IxGen, nGen, SNF);
  if (nrSNF < 0) return IE(-1);

      nd = SmithNormalForm(SNF, nrSNF, 3, NULL, Q);
  if (nd < 0 || nd > 3) return IE(-1);

  DTBF = 1; range1(id, 3) DTBF = iLCM(DTBF, SNF[(nd + 1) * id]);

  ResetLLTr(LLTr, &nDiscrGr);

  range1(id, nd) {
    d = SNF[(nd + 1) * id];
    range2(f, 1, d) {
      rangei(3) xp[i] = 0;
      xp[id] = f * DTBF / d;
      iMxMultiply(x, xp, Q, 1, 3, 3);
      if (ExpLLTr(DTBF, mDiscrGr, LLTr, &nDiscrGr, x) < 0)
        return IE(-1);
    }
  }

  range1(iDG, nDiscrGr) {
    MemCpy(DiscrGr[iDG].P, LLTr[iDG].v, 3);
    RotMx_t_Vector(DiscrGr[iDG].Z, Z2PCBMx[1].s.R, DiscrGr[iDG].P, 0);
    rangei(3)
      DiscrGr[iDG].Z[i] = iModPositive(DiscrGr[iDG].Z[i], DTBF * CRBF);
  }

  if (BestVect(SgOps, ssVM, n_ssVM, DTBF, DiscrGr, nDiscrGr) != 0)
    return IE(-1);

  qsort((void *) DiscrGr, nDiscrGr, sizeof (*DiscrGr),
        (int (*)(const void *, const void *)) CmpDiscr);

      nIx = SelectDiscrete(DTBF, nDiscrGr, DiscrGr, 3 - n_ssVM, Ix);
  if (nIx < 0) return IE(-1);

  rangei(nIx) {
    if (n_ssVM >= 3) return IE(-1);
    MemCpy(ssVM[n_ssVM].V, DiscrGr[Ix[i] + 1].Z, 3);
    ssVM[n_ssVM].M = CancelBFGCD(ssVM[n_ssVM].V, 3, DTBF * CRBF);
    n_ssVM++;
  }

  qsort((void *) ssVM, n_ssVM, sizeof (*ssVM),
        (int (*)(const void *, const void *)) Cmp_ssVM);

  return n_ssVM;
}
