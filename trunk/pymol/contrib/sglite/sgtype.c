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
#include "sgrefset.h"


static int SetCountRtype(const T_SgOps *SgOps, int *CountRtype)
{
  int  iLSMx, Rtype, AbsRtype;


  for (Rtype = -6; Rtype <=6; Rtype++) CountRtype[Rtype + 6] = 0;

  for (iLSMx = 1; iLSMx < SgOps->nSMx; iLSMx++)
  {
        Rtype = GetRtype(SgOps->SMx[iLSMx].s.R);
    if (Rtype == 0) return IE(-1);

        AbsRtype = abs(Rtype);
    if (AbsRtype < 2 || AbsRtype == 5 || AbsRtype > 6)
      return IE(-1);

    CountRtype[Rtype + 6]++;
  }

  return 0;
}


int GetPG(const T_SgOps *SgOps)
{
  int  CountRtype[13];


  if (SetCountRtype(SgOps, CountRtype) != 0)
    return MGC_Unknown;

#define CountRtype(i)  CountRtype[(i) + 6]
#define CountProperOrder(i)  (CountRtype(-i) + CountRtype(i))

  if      (CountProperOrder(3) == 8)
  {
    if      (SgOps->nSMx == 12) {
      if      (SgOps->fInv == 1) return MGC_23;
      else if (SgOps->fInv == 2) return MGC_m3b;
    }
    else if (SgOps->nSMx == 24) {
      if      (SgOps->fInv == 1) {
        if (CountRtype( 4) == 6) return MGC_432;
        if (CountRtype(-4) == 6) return MGC_4b3m;
      }
      else if (SgOps->fInv == 2) return MGC_m3bm;
    }
  }
  else if (CountProperOrder(6) == 2)
  {
    if      (SgOps->nSMx ==  6) {
      if      (SgOps->fInv == 1) {
        if (CountRtype( 6) == 2) return MGC_6;
        if (CountRtype(-6) == 2) return MGC_6b;
      }
      else if (SgOps->fInv == 2) return MGC_6_m;
    }
    else if (SgOps->nSMx == 12) {
      if      (SgOps->fInv == 1) {
        if      (CountRtype( 6) == 2) {
          if (CountRtype( 2) == 7) return MGC_622;
          if (CountRtype(-2) == 6) return MGC_6mm;
        }
        else if (CountRtype(-6) == 2) return MGC_6bm2;
      }
      else if (SgOps->fInv == 2) return MGC_6_mmm;
    }
  }
  else if (CountProperOrder(3) == 2)
  {
    if      (SgOps->nSMx ==  3) {
      if      (SgOps->fInv == 1) return MGC_3;
      else if (SgOps->fInv == 2) return MGC_3b;
    }
    else if (SgOps->nSMx ==  6) {
      if      (SgOps->fInv == 1) {
        if (CountRtype( 2) == 3) return MGC_32;
        if (CountRtype(-2) == 3) return MGC_3m;
      }
      else if (SgOps->fInv == 2) return MGC_3bm;
    }
  }
  else if (CountProperOrder(4) == 2)
  {
    if      (SgOps->nSMx ==  4) {
      if      (SgOps->fInv == 1) {
        if (CountRtype( 4) == 2) return MGC_4;
        if (CountRtype(-4) == 2) return MGC_4b;
      }
      else if (SgOps->fInv == 2) return MGC_4_m;
    }
    else if (SgOps->nSMx ==  8) {
      if      (SgOps->fInv == 1) {
        if      (CountRtype( 4) == 2) {
          if (CountRtype( 2) == 5) return MGC_422;
          if (CountRtype(-2) == 4) return MGC_4mm;
        }
        else if (CountRtype(-4) == 2) return MGC_4bm2;
      }
      else if (SgOps->fInv == 2) return MGC_4_mmm;
    }
  }
  else if (CountProperOrder(2) == 3)
  {
    if      (SgOps->fInv == 1) {
      if (CountRtype( 2) == 3) return MGC_222;
      if (CountRtype(-2) == 2) return MGC_mm2;
    }
    else if (SgOps->fInv == 2) return MGC_mmm;
  }
  else if (CountProperOrder(2) == 1)
  {
    if      (SgOps->fInv == 1) {
      if (CountRtype( 2) == 1) return MGC_2;
      if (CountRtype(-2) == 1) return MGC_m;
    }
    else if (SgOps->fInv == 2) return MGC_2_m;
  }
  else if (SgOps->nSMx == 1)
  {
    if      (SgOps->fInv == 1) return MGC_1;
    else if (SgOps->fInv == 2) return MGC_1b;
  }

#undef CountRtype
#undef CountProperOrder

  return IE(MGC_Unknown);
}


static int GetMG(const T_SgOps *StdSgOps, int PG_MGC)
{
  int     TwoFold, Mirror;
  int     iLSMx, Rtype;
  T_RMxI  RI[1];


  if (PG_MGC == MGC_Undefined) PG_MGC = GetPG(StdSgOps);
  if (PG_MGC == MGC_Unknown) return MGC_Unknown;

  TwoFold = 0;
  Mirror  = 0;

  if      (   PG_MGC == MGC_4bm2
           || PG_MGC == MGC_6bm2)
    TwoFold = 1;
  else if (StdSgOps->nLTr == 1) {
    if      (PG_MGC == MGC_32)  TwoFold = 1;
    else if (PG_MGC == MGC_3m)  Mirror  = 1;
    else if (PG_MGC == MGC_3bm) TwoFold = Mirror = 1;
  }

  if (! (TwoFold || Mirror))
    return PG_MGC;

  range2(iLSMx, 1, StdSgOps->nSMx)
  {
        Rtype = GetRtype(StdSgOps->SMx[iLSMx].s.R);
    if (Rtype == 0) return IE(MGC_Unknown);

    if (   (Rtype ==  2 && TwoFold)
        || (Rtype == -2 && Mirror))
    {
      if (SetRotMxInfo(StdSgOps->SMx[iLSMx].s.R, RI) == 0)
        return IE(MGC_Unknown);

      if (MemCmp(RI->EV, EV_100, 3) == 0)
      {
        if (PG_MGC == MGC_4bm2) return MGC_4b2m;
        if (PG_MGC == MGC_32)   return MGC_321;
        if (PG_MGC == MGC_3m)   return MGC_3m1;
        if (PG_MGC == MGC_3bm)  return MGC_3bm1;
        if (PG_MGC == MGC_6bm2) return MGC_6b2m;
      }
    }
  }

  if (PG_MGC == MGC_4bm2) return MGC_4bm2;
  if (PG_MGC == MGC_32)   return MGC_312;
  if (PG_MGC == MGC_3m)   return MGC_31m;
  if (PG_MGC == MGC_3bm)  return MGC_3b1m;
  if (PG_MGC == MGC_6bm2) return MGC_6bm2;

  return IE(MGC_Unknown);
}


static int StartBasis(const T_SgOps *SgOps,
                      int nWanted, int *IxS, int *Ord, int (*EVs)[3])
{
  int     Restart, iLSMx, iWtd, jWtd, RMxO, nFound, UseThisSMx;
  T_RMxI  RI[1];


  for (iWtd = 0; iWtd < nWanted; iWtd++) IxS[iWtd] = -1;

  nFound = 0;

  for (;;)
  {
    Restart = 0;

    for (iLSMx = 0; iLSMx < SgOps->nSMx; iLSMx++)
    {
          RMxO = GetRtype(SgOps->SMx[iLSMx].s.R);
      if (RMxO == 0) return IE(-1);

      for (iWtd = 0; iWtd < nWanted; iWtd++)
      {
        if (IxS[iWtd] < 0 && (RMxO == Ord[iWtd] || -RMxO == Ord[iWtd]))
        {
          if (SetRotMxInfo(SgOps->SMx[iLSMx].s.R, RI) == 0)
            return IE(-1);

          (void) MemCpy(EVs[iWtd], RI->EV, 3);

          if (RI->SenseOfRotation >= 0)
          {
            UseThisSMx = 1;

            for (jWtd = 0; jWtd < nWanted; jWtd++)
            {
              if (   IxS[jWtd] >= 0
                  && MemCmp(EVs[iWtd], EVs[jWtd], 3) == 0)
              {
                if (abs(Ord[jWtd]) >= abs(RMxO))
                  UseThisSMx = 0;
                else {
                  IxS[jWtd] = -1;
                  nFound--;
                  Restart = 1;
                }

                break;
              }
            }

            if (UseThisSMx)
            {
              Ord[iWtd] = RMxO;
              IxS[iWtd] = iLSMx;

                  nFound++;
              if (nFound == nWanted)
                return 0;
            }
          }
          break;
        }
      }
    }
    if (Restart == 0) break;
  }

  return IE(-1);
}


static int SetBasis(const int *R, int Rtype, int Basis[9])
{
  int        i;
  int        ProperOrder, CumMx[9];
  int        M_ProperR[9];
  const int   *ProperR;
  T_RMxI     RI[1];
  int        nIndep, IxIndep[2], Sol[4][3];
  int        nIx, Ix[2];
  int        Det, MinDet, TrialBasis[9];


      ProperOrder = abs(Rtype);
  if (ProperOrder == 1) {
    rangei(9)    Basis[i] = 0;
    range3(i, 0, 9, 4) Basis[i] = 1;
    return 0;
  }

  ProperR = R;

  if (Rtype < 0) {
    rangei(9) M_ProperR[i] = -R[i];
    ProperR = M_ProperR;
  }

  if (SetRotMxInfo(ProperR, RI) < 0) return IE(-1);

  (void) MakeCumRMx(ProperR, ProperOrder, CumMx);
  if (iRowEchelonFormT(CumMx, 3, 3, NULL, 0) != 1)
    return IE(-1);

      nIndep = iRESetIxIndep(CumMx, 1, 3, IxIndep, 2);
  if (nIndep != 2) return IE(-1);

  if (SolveHomRE1(CumMx, IxIndep, Sol) != 0) return -1;

  MinDet = 0;

  (void) MemCpy(&TrialBasis[6], RI->EV, 3);

  nIx = 1; if (ProperOrder == 2) nIx++;

  rangei(nIx) Ix[i] = i;

  do
  {
    rangei(nIx) (void) MemCpy(&TrialBasis[i * 3], Sol[Ix[i]], 3);

    if (nIx == 1) RotMx_t_Vector(&TrialBasis[3], ProperR, &TrialBasis[0], 0);

        Det = deterRotMx(TrialBasis);
    if (Det != 0 && (MinDet == 0 || abs(MinDet) > abs(Det))) {
      MinDet = Det;
      (void) MemCpy(Basis, TrialBasis, 9);
    }
  }
  while (NextOf_n_from_m(4, nIx, Ix) != 0);

  if (MinDet < 0) IntSwap(&Basis[0], &Basis[3], 3);

  return 0;
}


static int StdBasis(const T_SgOps *SgOps, int MGC, int Basis[9])
{
  int  nWtd, Ord[3], IxS[3], EVs[3][3], i;


  if (MGC == MGC_Undefined) MGC = GetPG(SgOps);
  if (MGC == MGC_Unknown) return -1;

  switch (ixLG(MGC)) {
    case ixPG(MGC_1b):    nWtd = 1; Ord[0] = 1;                         break;
    case ixPG(MGC_2_m):   nWtd = 1; Ord[0] = 2;                         break;
    case ixPG(MGC_mmm):   nWtd = 3; Ord[0] = 2; Ord[1] = 2; Ord[2] = 2; break;
    case ixPG(MGC_4_m):   nWtd = 1; Ord[0] = 4;                         break;
    case ixPG(MGC_4_mmm): nWtd = 2; Ord[0] = 4; Ord[1] = 2;             break;
    case ixPG(MGC_3b):    nWtd = 1; Ord[0] = 3;                         break;
    case ixPG(MGC_3bm):   nWtd = 2; Ord[0] = 3; Ord[1] = 2;             break;
    case ixPG(MGC_6_m):   nWtd = 1; Ord[0] = 3;                         break;
    case ixPG(MGC_6_mmm): nWtd = 2; Ord[0] = 3; Ord[1] = 2;             break;
    case ixPG(MGC_m3b):   nWtd = 2; Ord[0] = 3; Ord[1] = 2;             break;
    case ixPG(MGC_m3bm):  nWtd = 2; Ord[0] = 3; Ord[1] = 4;             break;
    default:
      return IE(-1);
  }

  if (StartBasis(SgOps, nWtd, IxS, Ord, EVs) != 0)
    return -1;

  if (nWtd == 1)
  {
    if (SetBasis(SgOps->SMx[IxS[0]].s.R, Ord[0], Basis) != 0)
      return -1;
  }
  else
  {
    (void) MemCpy(&Basis[0], EVs[1], 3);

    if (ixLG(MGC) == ixPG(MGC_m3b) || ixLG(MGC) == ixPG(MGC_m3bm))
    {
      range3(i, 0, 6, 3)
        RotMx_t_Vector(&Basis[i + 3], SgOps->SMx[IxS[0]].s.R, &Basis[i], 0);
    }
    else
    {
      (void) MemCpy(&Basis[6], EVs[0], 3);

      if (nWtd == 3) {
        (void) MemCpy(&Basis[3], EVs[2], 3);
      }
      else {
        RotMx_t_Vector(&Basis[3], SgOps->SMx[IxS[0]].s.R, &Basis[0], 0);
        if (Ord[0] < 0) rangei(3) Basis[3 + i] *= -1;
      }
    }

    if (deterRotMx(Basis) < 0) IntSwap(&Basis[0], &Basis[3], 3);
  }

  return 0;
}


static int Basis2CBMx(const int Basis[9], int BF, T_RTMx *CBMx, T_RTMx *InvCBMx)
{
  int     i;
  T_RTMx  CBMxBuf[2];


  if (   CBMx == NULL)    CBMx = &CBMxBuf[0];
  if (InvCBMx == NULL) InvCBMx = &CBMxBuf[1];

  MemCpy(InvCBMx->s.R, Basis, 9);
  if (TransposedMat(InvCBMx->s.R, 3, 3) == NULL) return -1;
  if (ChangeBaseFactor(InvCBMx->s.R, BF, InvCBMx->s.R, CRBF, 9) != 0) {
    SetSgError("Error: Out of change-of-basis rotation-base-factor range");
    return 0;
  }

  rangei(3) InvCBMx->s.T[i] = 0;

  if ((i = InverseRTMx(InvCBMx, CBMx, CRBF)) == 0) {
    SetSgError("Error: Change-of-basis operator is not invertible");
    return 0;
  }

  return i;
}


static int SetStdIxGen(const T_SgOps *StdSgOps, int PG_MGC, int IxGen[2])
{
  int     nGen, iLSMx, i;
  int     PrincipalProperOrder;
  T_RMxI  RI[1];


  rangei(2) IxGen[i] = -1;

  if (PG_MGC == MGC_Undefined) PG_MGC = GetPG(StdSgOps);
  if (PG_MGC == MGC_Unknown) return -1;

  nGen = 0;
  PrincipalProperOrder = 0;

  switch (ixXS(PG_MGC))
  {
    case XS_Triclinic:
      if (StdSgOps->fInv == 1) {
        IxGen[0] = 0;
        nGen = 1;
      }
      break;

    case XS_Monoclinic:
      IxGen[0] = 1;
      nGen = 1;
      break;

    case XS_Orthorhombic:
      range2(iLSMx, 1, StdSgOps->nSMx)
      {
        if (SetRotMxInfo(StdSgOps->SMx[iLSMx].s.R, RI) == 0)
          return IE(-1);

        if      (MemCmp(RI->EV, EV_001, 3) == 0) IxGen[0] = iLSMx;
        else if (MemCmp(RI->EV, EV_100, 3) == 0) IxGen[1] = iLSMx;
      }

      rangei(2) if (IxGen[i] < 0) return IE(-1);

      nGen = 2;
      break;

    case XS_Tetragonal:
                                  PrincipalProperOrder = 4;
    case XS_Trigonal:
      if (! PrincipalProperOrder) PrincipalProperOrder = 3;
    case XS_Hexagonal:
      if (! PrincipalProperOrder) PrincipalProperOrder = 6;

      range2(iLSMx, 1, StdSgOps->nSMx)
      {
        if (SetRotMxInfo(StdSgOps->SMx[iLSMx].s.R, RI) == 0)
          return IE(-1);

        if (abs(RI->Rtype) == PrincipalProperOrder) {
            if (RI->SenseOfRotation > 0)      IxGen[0] = iLSMx;
        }
        else if (PrincipalProperOrder == 4) {
          if (MemCmp(RI->EV, EV_100, 3) == 0) IxGen[1] = iLSMx;
        }
        else if (PrincipalProperOrder == 3) {
          if (MemCmp(RI->EV, EV_m10, 3) == 0) IxGen[1] = iLSMx;
          if (MemCmp(RI->EV, EV_110, 3) == 0) IxGen[1] = iLSMx;
        }
        else {
          if (MemCmp(RI->EV, EV_m10, 3) == 0) IxGen[1] = iLSMx;
        }
      }

      if (IxGen[0] < 0) return IE(-1);

      nGen++; if (IxGen[1] > 0) nGen++;
      break;

    case XS_Cubic:
      range2(iLSMx, 1, StdSgOps->nSMx)
      {
        if (SetRotMxInfo(StdSgOps->SMx[iLSMx].s.R, RI) == 0)
          return IE(-1);

        if      (abs(RI->Rtype) == 4 && RI->SenseOfRotation > 0) {
          if (MemCmp(RI->EV, EV_001, 3) == 0) IxGen[0] = iLSMx;
        }
        else if (abs(RI->Rtype) == 2 && IxGen[0] < 0) {
          if (MemCmp(RI->EV, EV_001, 3) == 0) IxGen[0] = iLSMx;
        }
        else if (abs(RI->Rtype) == 3 && RI->SenseOfRotation > 0) {
          if (MemCmp(RI->EV, EV_111, 3) == 0) IxGen[1] = iLSMx;
        }
      }

      rangei(2) if (IxGen[i] < 0) return IE(-1);

      nGen = 2;
      break;

    default:
      return IE(-1);
  }

  return nGen;
}


static void MvGenFirst(T_SgOps *StdSgOps, int IxGen[2])
{
  int     iIx, jIx, iLSMx;
  T_RTMx  SMx[1];


  range1(iIx, 2)
  {
    if (IxGen[iIx] < 1) break;

        iLSMx = iIx + 1;
    if (iLSMx != IxGen[iIx])
    {
      (void) MemCpy(SMx, &StdSgOps->SMx[iLSMx], 1);
      (void) MemCpy(&StdSgOps->SMx[iLSMx], &StdSgOps->SMx[IxGen[iIx]], 1);
      (void) MemCpy(&StdSgOps->SMx[IxGen[iIx]], SMx, 1);

      range2(jIx, iIx + 1, 2) {
        if (iLSMx == IxGen[jIx]) {
          IxGen[jIx] = IxGen[iIx];
          break;
        }
      }

      IxGen[iIx] = iLSMx;
    }
  }
}


static int MkGenRStd(T_SgOps *StdSgOps, int nGen)
{
  int     i, iLSMx, Rtype;
  T_RTMx  *SMx;


  if (StdSgOps->nSMx > 1 && StdSgOps->fInv == 2)
  {
    for (iLSMx = 1; iLSMx < nGen + 1; iLSMx++)
    {
                           SMx = &StdSgOps->SMx[iLSMx];
          Rtype = GetRtype(SMx->s.R);
      if (Rtype == 0) return IE(-1);
      if (Rtype < 0)
        SMx_t_InvT(SMx, StdSgOps->InvT, SMx);
        rangei(3) SMx->s.T[i] = iModPositive(SMx->s.T[i], STBF);
    }
  }

  return 0;
}


static int TidyGen(T_SgOps *StdSgOps, int PG_MGC)
{
  int  IxGen[2], nGen;


      nGen = SetStdIxGen(StdSgOps, PG_MGC, IxGen);
  if (nGen < 0) return -1;

  MvGenFirst(StdSgOps, IxGen);

  if (MkGenRStd(StdSgOps, nGen) != 0) return -1;

  return nGen;
}


static int UpdateCBMxT(T_RTMx CBMx[2], const int CBT[3])
{
  int  i;


  rangei(3) CBMx[0].s.T[i] = iModPositive(CBT[i], CTBF);

  if (InverseRTMx(&CBMx[0], &CBMx[1], CRBF) == 0)
    return IE(-1);

  rangei(3) CBMx[1].s.T[i] = iModPositive(CBMx[1].s.T[i], CTBF);

  return 1;
}


static int PrimitiveGenerators(const T_SgOps *SgOps, const int nGen,
                               const T_RTMx Z2PCBMx[2],
                               T_RTMx *PSMx)
{
  int  iGen, n, i;

  iGen = 0;

  if (SgOps->nSMx > 1) {
    range1(iGen, nGen) {
      if (CB_SMx(&PSMx[iGen],
                 &Z2PCBMx[0], &SgOps->SMx[iGen + 1], &Z2PCBMx[1]) != 0)
        return -1;
    }
  }

  if (SgOps->fInv == 2)
  {
    InitRotMx(PSMx[iGen].s.R, -1);
    if (CB_IT(-1, SgOps->InvT, &Z2PCBMx[0], &Z2PCBMx[1], PSMx[iGen].s.T) != 0)
      return -1;

    iGen++;
  }

  n = iGen;

  range1(iGen, n)
    rangei(3) PSMx[iGen].s.T[i] = iModPositive(PSMx[iGen].s.T[i], STBF);

  return n;
}


static int PrimitiveSMxT(const T_SgOps *SgOps, const int nGen,
                         const int Z2PCBMxR[9],
                         int PSMxT[3][3])
{
  int  iGen, n, i;


  iGen = 0;

  if (SgOps->nSMx > 1)
    range1(iGen, nGen)
      RotMx_t_Vector(PSMxT[iGen], Z2PCBMxR, SgOps->SMx[iGen + 1].s.T, 0);

  if (SgOps->fInv == 2) {
    RotMx_t_Vector(PSMxT[iGen], Z2PCBMxR, SgOps->InvT, 0);
    iGen++;
  }

  n = iGen;

  range1(iGen, n) {
    rangei(3) {
      if (PSMxT[iGen][i] %  CRBF) return IE(-1);
          PSMxT[iGen][i] /= CRBF;
          PSMxT[iGen][i] = iModPositive(PSMxT[iGen][i], STBF);
    }
  }

  return n;
}


static int SolveInhomModZ(int *M, int nr, int nc, int *b, int BF, int *x)
{
#define maxr 9
#define maxc 6

  int  nd, d, i;
  int  P[maxr * maxr], Q[maxc * maxc];
  int  Pb[maxr], xp[maxc];


  if (nr > maxr) return IE(-1);
  if (nc > maxc) return IE(-1);

      nd = SmithNormalForm(M, nr, nc, P, Q);
  if (nd < 0 || nd > nc) return IE(-1);

  iMxMultiply(Pb, P, b, nr, nr, 1);

  range2(i, nd, nr) if (Pb[i] % BF != 0) return 0;

  if (x)
  {
    rangei(nc) {
      xp[i] = 0;
          d = M[i * nd + i];
      if (d) {
            xp[i] = Pb[i];
        if (xp[i] %  d) return -1;
            xp[i] /= d;
      }
    }

    iMxMultiply(x, xp, Q, 1, nc, nc);
  }

  return nd + 1;

#undef maxr
#undef maxc
}


static int FindOShift(const T_SgOps *TstSgOps, const int nGen,
                      const T_RTMx Z2PCBMx[2], const T_RTMx TabPSMx[3],
                      int CBT[3])
{
  /*    (I|K)(R|T)(I|-K)=(R|S)
     => S=-RK+T+K=-(R-I)K+T
     => S=-(R-I)K+T
     => (R-I)K=T-S
     => (R-I)^-1(T-S)=K
   */

  int  nPSMx, iPSMx, i;
  int  TmSV[9], x[3];
  int  SNF[9 * 3], nrSNF;


      nPSMx = PrimitiveSMxT(TstSgOps, nGen, Z2PCBMx[0].s.R, (int(*)[3]) TmSV);
  if (nPSMx < 1) return IE(-1);

  range1(iPSMx, nPSMx)
    rangei(3) TmSV[iPSMx * 3 + i] -= TabPSMx[iPSMx].s.T[i];

  nrSNF = nPSMx * 3;

  rangei(nrSNF) TmSV[i] *= (CTBF / STBF);

  range1(iPSMx, nPSMx)
    SetRminusI(TabPSMx[iPSMx].s.R, &SNF[iPSMx * 9], 0);

      i = SolveInhomModZ(SNF, nrSNF, 3, TmSV, CTBF, x);
  if (i < 0) return IE(-1);
  if (i == 0) return 0;

  RotMx_t_Vector(CBT, Z2PCBMx[1].s.R, x, 0);

  if (ChangeBaseFactor(CBT, CRBF, CBT, 1, 3) != 0)
    return IE(-1);

  return 1;
}


static int MatchGenerators(T_SgOps *StdSgOps,
                           T_SgOps *TabSgOps, const int MGC,
                           T_RTMx CBMx[2])
{
  int           nGen, Stat;
  int                i2fold,      i3fold;
  const T_RTMx  *CBMx_2fold, *CBMx_3fold;
  int           TabCType, TstCType;
  T_SgOps       TstSgOps[1];
  T_RTMx        TabZ2PCBMx[2], TabPSMx[3];
  int           CBT[3];


  InitRTMx(&CBMx[0], CRBF);
  InitRTMx(&CBMx[1], CRBF);

  if (TabSgOps->nSMx == 1 && TabSgOps->fInv == 1) return 1;

      nGen = TidyGen(TabSgOps, MGC);
  if (nGen < 0 || nGen > 2) return IE(-1);

  if (GetZ2PCBMx(TabSgOps, TabZ2PCBMx) != 0) return -1;
  if (PrimitiveGenerators(TabSgOps, nGen, TabZ2PCBMx, TabPSMx) < 1)
    return IE(-1);

  CBMx_2fold = NULL;
  CBMx_3fold = NULL;

  if      (ixXS(MGC) == XS_Monoclinic) {
    CBMx_2fold = CBMx_2_101;
    CBMx_3fold = CBMx_3_010;
  }
  else if (ixXS(MGC) == XS_Orthorhombic) {
    CBMx_2fold = CBMx_2_110;
    CBMx_3fold = CBMx_3_111;
  }

  if (CBMx_2fold)
  {
        TabCType = GetSymCType(TabSgOps->nLTr, TabSgOps->LTr);
    if (TabCType == '\0' || TabCType == 'Q')
      return IE(-1);

    range1(i2fold, 2)
    {
      if (i2fold)
        (void) MemCpy(&CBMx[0], CBMx_2fold, 1);

      range1(i3fold, 3)
      {
        if (i3fold) {
          if (CBMxMultiply(&CBMx[0], CBMx_3fold, &CBMx[0]) != 0) return -1;
        }

        if (InverseRTMx(&CBMx[0], &CBMx[1], CRBF) == 0)
          return IE(-1);

        ResetSgOps(TstSgOps);
        if (CB_SgOps(StdSgOps, &CBMx[0], &CBMx[1], TstSgOps) != 0) return -1;

            TstCType = GetSymCType(TstSgOps->nLTr, TstSgOps->LTr);
        if (TstCType == '\0' || TstCType == 'Q')
          return IE(-1);

        if (TstCType != TabCType)
          continue;

        if (nGen != TidyGen(TstSgOps, MGC)) return IE(-1);

        if (    nGen != 2
            ||  (   TabSgOps->SMx[1].s.R[8] == TstSgOps->SMx[1].s.R[8]
                 && TabSgOps->SMx[2].s.R[0] == TstSgOps->SMx[2].s.R[0]))
        {
              Stat = FindOShift(TstSgOps, nGen, TabZ2PCBMx, TabPSMx, CBT);
          if (Stat < 0) return -1;
          if (Stat > 0) return UpdateCBMxT(CBMx, CBT);
        }
      }
    }
  }
  else
  {
    if (nGen != TidyGen(StdSgOps, MGC)) return IE(-1);

        Stat = FindOShift(StdSgOps, nGen, TabZ2PCBMx, TabPSMx, CBT);
    if (Stat < 0) return -1;
    if (Stat > 0) return UpdateCBMxT(CBMx, CBT);
  }

  return 0;
}


static int OShSMxT(const T_RTMx *SMx, const int CBMxT[3], int OshT[3])
{
  int  i;


  RotMx_t_Vector(OshT, SMx->s.R, CBMxT, 0);

  rangei(3)
  {
        OshT[i] -= CBMxT[i];
    if (OshT[i] % (CTBF / STBF)) return IE(-1);
        OshT[i] = iModPositive(SMx->s.T[i] - OshT[i] / (CTBF / STBF), STBF);
  }

  return 0;
}


static int m3bWrongGlide(const T_SgOps *StdSgOps)
{
  int     iSMx, iInv;
  int     Rtype, wI[3];
  T_RMxI  RI[1];
  T_RTMx  LISMx[1];


  if (StdSgOps->fInv != 2) return IE(-1);

  range2(iSMx, 1, StdSgOps->nSMx)
  {
            Rtype = GetRtype(StdSgOps->SMx[iSMx].s.R);
    if (    Rtype == 0) return IE(-1);
    if (abs(Rtype) == 2) {
      if (SetRotMxInfo(StdSgOps->SMx[iSMx].s.R, RI) == 0)
        return IE(-1);
      if (MemCmp(RI->EV, EV_100, 3) == 0)
      {
        iInv = 0; if (Rtype == 2) iInv = 1;
        SetLISMx(StdSgOps, 0, iInv, iSMx, LISMx);
        if (Set_wI_Tr(LISMx->a, NULL, NULL, wI, NULL) != 0)
          return IE(-1);
        if (wI[2] % STBF != 0) return 1;
        return 0;
      }
    }
  }

  return IE(-1);
}


static int FindPreShiftSgOps(const T_SgOps *SgOps, T_RTMx CBMx[2])
{
  /* Shifting the symmetry operators of trigonal and hexagonal
     space groups to the origin before calling StdBasis() allows
     us to work with a STBF < 36.
   */

  int  iSMx, i;
  int  CBT[3], SMxT[3], wI[3], Tr[3];


  rangei(3) CBT[i] = 0;

  if (SgOps->fInv == 2)
  {
    rangei(3) CBT[i] = -SgOps->InvT[i] * (CTBF / STBF / 2);
  }
  else
  {
    range2(iSMx, 1, SgOps->nSMx)
    {
      if (OShSMxT(&SgOps->SMx[iSMx], CBT, SMxT) != 0) return -1;

      if (Set_wI_Tr(SgOps->SMx[iSMx].s.R, SMxT, NULL, wI, Tr) != 0)
        return IE(-1);

      rangei(3) CBT[i] -= Tr[i];
    }
  }

  InitRotMx(CBMx[0].s.R, CRBF);
  (void) UpdateCBMxT(CBMx, CBT);

  return 0;
}


int GetSpaceGroupType(const T_SgOps *SgOps, T_RTMx *CBMx, T_RTMx *InvCBMx)
{
  int           PG_MGC, MGC, SgNumber, Stat, RunAwayCounter;
  int           SymCType, MatchSymCType;
  const char    *HallSymbol;
  T_SgOps       StdSgOps[1], TabSgOps[1];
  int           Basis[9];
  T_RTMx        TotCBMx[2], AddCBMx[2];
  const T_RTMx  *AdjCBMx;


  if (   CBMx) InitRTMx(   CBMx, 0);
  if (InvCBMx) InitRTMx(InvCBMx, 0);

  InitRTMx(&TotCBMx[0], CRBF);
  InitRTMx(&TotCBMx[1], CRBF);
  SgOpsCpy(StdSgOps, SgOps);

  RunAwayCounter = 0;

  do
  {
    if (RunAwayCounter++ > 10) return IE(-1);

    if (GetZ2PCBMx(StdSgOps, AddCBMx) != 0) return -1;
    if (CBMx2Update(TotCBMx, AddCBMx) != 0) return -1;
    ResetSgOps(StdSgOps);
    if (CB_SgOps(SgOps, &TotCBMx[0], &TotCBMx[1], StdSgOps) != 0) return -1;

    if (StdSgOps->nLTr != 1) return IE(-1);

        PG_MGC = GetPG(SgOps);
    if (PG_MGC == MGC_Unknown) return -1;

    if (ixXS(PG_MGC) == XS_Trigonal || ixXS(PG_MGC) == XS_Hexagonal)
    {
      if (FindPreShiftSgOps(StdSgOps, AddCBMx) != 0) return -1;
      if (CBMx2Update(TotCBMx, AddCBMx) != 0) return -1;
      ResetSgOps(StdSgOps);
      if (CB_SgOps(SgOps, &TotCBMx[0], &TotCBMx[1], StdSgOps) != 0) return -1;
    }

    if (StdBasis(StdSgOps, PG_MGC, Basis) != 0) return -1;
    if (Basis2CBMx(Basis, 1, &AddCBMx[0], &AddCBMx[1]) == 0) return -1;
    if (CBMx2Update(TotCBMx, AddCBMx) != 0) return -1;
    ResetSgOps(StdSgOps);
    if (CB_SgOps(SgOps, &TotCBMx[0], &TotCBMx[1], StdSgOps) != 0) return -1;

    SymCType = GetSymCType(StdSgOps->nLTr, StdSgOps->LTr);

    AdjCBMx = NULL;

    if      (    ixLG(PG_MGC) == ixPG(MGC_2_m))
      AdjCBMx = CBMxMon_c_b;
    else if (   (ixLG(PG_MGC) == ixPG(MGC_4_m)
              || ixLG(PG_MGC) == ixPG(MGC_4_mmm))
             && SymCType == 'C')
      AdjCBMx = CBMxCP;
    else if (   (ixLG(PG_MGC) == ixPG(MGC_4_m)
              || ixLG(PG_MGC) == ixPG(MGC_4_mmm))
             && SymCType == 'F')
      AdjCBMx = CBMxFI;
    else if (   (ixLG(PG_MGC) == ixPG(MGC_3b)
              || ixLG(PG_MGC) == ixPG(MGC_3bm))
             && SymCType == 'Q')
      AdjCBMx = CBMxRevObv;
    else if (   (ixLG(PG_MGC) == ixPG(MGC_3bm)
              || ixLG(PG_MGC) == ixPG(MGC_6_mmm))
             && SymCType == 'H')
      AdjCBMx = CBMxHP;
    else if (    ixPG(PG_MGC) == ixPG(MGC_m3b)
             && SymCType == 'P') {
          Stat = m3bWrongGlide(StdSgOps);
      if (Stat < 0) return -1;
      if (Stat) AdjCBMx = CBMx_4_001;
    }

    if (AdjCBMx)
    {
      if (CBMx2Update(TotCBMx, AdjCBMx) != 0) return -1;
      ResetSgOps(StdSgOps);
      if (CB_SgOps(SgOps, &TotCBMx[0], &TotCBMx[1], StdSgOps) != 0) return -1;
      SymCType = GetSymCType(StdSgOps->nLTr, StdSgOps->LTr);
    }
  }
  while (SymCType == '\0');

  if (SymCType == 'Q') return IE(-1);

      MGC = GetMG(StdSgOps, PG_MGC);
  if (MGC == MGC_Unknown) return -1;

  if (   ixXS(PG_MGC) != ixXS(MGC)
      || ixLG(PG_MGC) != ixLG(MGC)
      || ixPG(PG_MGC) != ixPG(MGC)) return IE(-1);

  MatchSymCType = (       ixXS(MGC) != XS_Monoclinic
                   && (   ixXS(MGC) != XS_Orthorhombic
                       || (SymCType == 'I' || SymCType == 'F')));

  for (SgNumber = 1; SgNumber <= 230; SgNumber++)
  {
    HallSymbol = RefSetHallSymbols[SgNumber];

    if (MatchSymCType && SymCType != HallSymbol[1])
      continue;

    if ((SymCType == 'P') != (HallSymbol[1] == 'P'))
      continue;

    if (RefSetMGC[SgNumber] != MGC)
      continue;

    ResetSgOps(TabSgOps);
    if (ParseHallSymbol(HallSymbol, TabSgOps, PHSymOptPedantic) < 0)
      return -1;

    if (TabSgOps->nLTr != StdSgOps->nLTr)
      continue;

        Stat = MatchGenerators(StdSgOps, TabSgOps, MGC, AddCBMx);
    if (Stat < 0) return -1;
    if (Stat == 1)
    {
      if (CBMx2Update(TotCBMx, AddCBMx) != 0) return -1;

      if (deterRotMx(TotCBMx[0].s.R) <= 0) return IE(-1);
      if (deterRotMx(TotCBMx[1].s.R) <= 0) return IE(-1);

      if (   CBMx) MemCpy(   CBMx, &TotCBMx[0], 1);
      if (InvCBMx) MemCpy(InvCBMx, &TotCBMx[1], 1);

      return SgNumber;
    }
  }

  return IE(-1);
}


static int TidyCBMxT(const T_SgOps *RefSgOps,
                     const T_SgOps *SgOps, T_RTMx CBMx[2])
{
  int      MGC, nGen, Stat;
  T_RTMx   TabZ2PCBMx[2], TabPSMx[3];
  int      CBT[3];
  T_SgOps  BufSgOps[1], BC_SgOps[1];


  SgOpsCpy(BufSgOps, RefSgOps);

  IntSetZero(CBMx[0].s.T, 3);
  IntSetZero(CBMx[1].s.T, 3);

  /* done if space group is P 1 */
  if (BufSgOps->nSMx == 1 && BufSgOps->fInv == 1) return 0;

      MGC = GetMG(BufSgOps, MGC_Undefined);
  if (MGC == MGC_Unknown) return IE(-1);

      nGen = TidyGen(BufSgOps, MGC);
  if (nGen < 0 || nGen > 2) return IE(-1);

  if (GetZ2PCBMx(BufSgOps, TabZ2PCBMx) != 0) return -1;
  if (PrimitiveGenerators(BufSgOps, nGen, TabZ2PCBMx, TabPSMx) < 1)
    return IE(-1);

  ResetSgOps(BC_SgOps);
  if (CB_SgOps(SgOps, &CBMx[0], &CBMx[1], BC_SgOps) != 0)
    return IE(-1);

  if (nGen != TidyGen(BC_SgOps, MGC)) return IE(-1);

      Stat = FindOShift(BC_SgOps, nGen, TabZ2PCBMx, TabPSMx, CBT);
  if (Stat <= 0) return IE(-1);

  if (UpdateCBMxT(CBMx, CBT) != 1) return -1;

  return 0;
}


static int CmpCBMx(const T_RTMx *a, const T_RTMx *b)
{
  int  na, nb, i;

  na = (MemCmp(a->s.R, CBMx_1_000->s.R, 9) == 0);
  nb = (MemCmp(b->s.R, CBMx_1_000->s.R, 9) == 0);
  if (  na && ! nb) return -1;
  if (! na &&   nb) return  1;

  na = IntIsZero(a->s.T, 3);
  nb = IntIsZero(b->s.T, 3);
  if (  na && ! nb) return -1;
  if (! na &&   nb) return  1;

  na = 0; rangei(9) if (a->s.R[i] == 0) na++;
  nb = 0; rangei(9) if (b->s.R[i] == 0) nb++;
  if (na > nb) return -1;
  if (na < nb) return  1;

  na = 0; rangei(9) if (abs(a->s.R[i]) == CRBF) na++;
  nb = 0; rangei(9) if (abs(b->s.R[i]) == CRBF) nb++;
  if (na > nb) return -1;
  if (na < nb) return  1;

  na = 0; rangei(9) if (a->s.R[i] > 0) na++;
  nb = 0; rangei(9) if (b->s.R[i] > 0) nb++;
  if (na > nb) return -1;
  if (na < nb) return  1;

      i = CmpiVect(a->s.T, b->s.T, 3);
  if (i != 0) return i;

  return CmpiVect(a->s.R, b->s.R, 9);
}


static void GetMonoRefSetAffNormTrialRanges(const int CBMxR[9],
                                            int r00[1], int r22[1])
{
  /* International Tables Volume A, chapter 15, tables 15.3.3 & 15.3.4.

   M.C = n00 * c00 + n02 * c20,  n00 * c01 + n02 * c21,  n00 * c02 + n02 * c22,
         c10,                    c11,                    c12,
         n20 * c00 + n22 * c20,  n20 * c01 + n22 * c21,  n20 * c02 + n22 * c22

   Determine trial range for n00 and n20:
     max(lcm(c00, c20) / c00,
         lcm(c01, c21) / c01,
         lcm(c02, c22) / c02)
   Determine trial range for n02 and n22:
     max(lcm(c00, c20) / c20,
         lcm(c01, c21) / c21,
         lcm(c02, c22) / c22)
   */

  int  i, l, n;

  r00[0] = 1;
  r22[0] = 1;
  rangei(3) {
    l = iLCM(CBMxR[i], CBMxR[6 + i]);
    if (CBMxR[i]) {
      n = abs(l / CBMxR[i]);
      if (r00[0] < n) r00[0] = n;
    }
    if (CBMxR[i + 6]) {
      n = abs(l / CBMxR[6 + i]);
      if (r22[0] < n) r22[0] = n;
    }
  }
  r00[0]++;
  r22[0]++;
}


static int getBestCBMx(const T_SgOps *SgOps, int SgNumber,
                       const T_SgOps *TdRefSgOps,
                       T_RTMx CBMx[2])
{
  int      nAddlG, iAddlG;
  T_RTMx    AddlG[3];
  T_SgOps  NormSgOps[1];
  int      iInv, iSMx, icmp, det, r00, r22, f;
  T_RTMx   SMx[1], LISMx[2], TrialCBMx[2], M[2], M_TrialCBMx[2], BestCBMx[2];

      nAddlG = GetRefSetNormAddlG(SgNumber, 1, 1, 1, AddlG);
  if (nAddlG < 0)
    return IE(-1);

  SgOpsCpy(NormSgOps, SgOps);
  range1(iAddlG, nAddlG) {
    if (CB_SMx(SMx, &CBMx[1], &AddlG[iAddlG], &CBMx[0]) != 0) return IE(-1);
    if (ExpSgSMx(NormSgOps, SMx) < 0) return IE(-1);
  }

  MemCpy(BestCBMx, CBMx, 2);
  icmp = 0;
  if (deterRotMx(CBMx[0].s.R) < deterRotMx(CBMx[1].s.R))
    icmp = 1;

  range1(iInv, NormSgOps->fInv)
  range1(iSMx, NormSgOps->nSMx)
  {
    SetLISMx(NormSgOps, 0, iInv, iSMx, &LISMx[0]);

        det = deterRotMx(LISMx[0].s.R);
    if (det == 0) return IE(-1);
    if (det < 0) continue;

    if (ChangeBaseFactor(LISMx[0].s.R,    1, LISMx[0].s.R, CRBF, 9) != 0)
      return IE(-1);
    if (ChangeBaseFactor(LISMx[0].s.T, STBF, LISMx[0].s.T, CTBF, 3) != 0)
      return IE(-1);
    if (InverseRTMx(&LISMx[0], &LISMx[1], CRBF) == 0)
      return IE(-1);

    if (CBMx2Multiply(TrialCBMx, CBMx, LISMx) != 0) return IE(-1);

    if (SgNumber < 3 || SgNumber > 15) {
      if (TidyCBMxT(TdRefSgOps, SgOps, TrialCBMx) != 0) return IE(-1);
      if (CmpCBMx(&BestCBMx[icmp], &TrialCBMx[icmp]) > 0)
        MemCpy(BestCBMx, TrialCBMx, 2);
    }
    else {
      IntSetZero(M[0].a, 12);
      IntSetZero(M[1].a, 12);
      GetMonoRefSetAffNormTrialRanges(TrialCBMx[0].s.R, &r00, &r22);
#define loop(i, r) \
      for (M[0].s.R[i] = -r*CRBF; M[0].s.R[i] <= r*CRBF; M[0].s.R[i] += CRBF)
      loop(0, r00)
      loop(2, r22)
      loop(6, r00)
      loop(8, r22) {
        if (CheckMonoRefSetAffNormRestrictions(SgNumber, M[0].s.R, CRBF) != 0)
          continue;
        M[0].s.R[4] = M[0].s.R[0] * M[0].s.R[8] - M[0].s.R[2] * M[0].s.R[6];
        if (   M[0].s.R[4] != -CRBF * CRBF
            && M[0].s.R[4] !=  CRBF * CRBF) continue;
        M[0].s.R[4] /= CRBF;
        /* set M[1] = inverse M[0] */
        f = M[0].s.R[4] / CRBF;
        M[1].s.R[0] =  f * M[0].s.R[8];
        M[1].s.R[2] = -f * M[0].s.R[2];
        M[1].s.R[4] =      M[0].s.R[4];
        M[1].s.R[6] = -f * M[0].s.R[6];
        M[1].s.R[8] =  f * M[0].s.R[0];
        if (CBMx2Multiply(M_TrialCBMx, M, TrialCBMx) != 0) return IE(-1);
        if (TidyCBMxT(TdRefSgOps, SgOps, M_TrialCBMx) != 0) return IE(-1);
        if (CmpCBMx(&BestCBMx[icmp], &M_TrialCBMx[icmp]) > 0)
          MemCpy(BestCBMx, M_TrialCBMx, 2);
      }
#undef loop
    }
  }

  MemCpy(CBMx, BestCBMx, 2);
  if (deterRotMx(CBMx[0].s.R) <= 0) return IE(-1);
  if (deterRotMx(CBMx[1].s.R) <= 0) return IE(-1);

  return 0;
}


int TidyCBMx(const T_SgOps *SgOps, int SgNumber, T_RTMx CBMx[2])
{
  T_SgOps  TdRefSgOps[1];

  if (SgNumber < 1 || SgNumber > 230) return IE(-1);

  ResetSgOps(TdRefSgOps);
  if (ParseHallSymbol(RefSetHallSymbols[SgNumber],
                      TdRefSgOps, PHSymOptPedantic) < 0) return IE(-1);
  if (TidySgOps(TdRefSgOps) != 0) return IE(-1);

  return getBestCBMx(SgOps, SgNumber, TdRefSgOps, CBMx);
}


int BuildHallSymbol(const T_SgOps *SgOps, int SgNumber, const T_RTMx CBMx[2],
                    char *HallSymbol, int sizeHallSymbol)
{
  int         HaveCBMx, iHS;
  char        xyz[128];
  const char  *RefHS;
  T_RTMx      RefCBMx[2], HallCBMx[2];
  T_SgOps     TdRefSgOps[1];

  if (SgNumber < 1 || SgNumber > 230) return IE(-1);
  RefHS = RefSetHallSymbols[SgNumber];

  ResetSgOps(TdRefSgOps);
  if (ParseHallSymbolCBMx(RefHS, TdRefSgOps, PHSymOptPedantic,
                          RefCBMx, &HaveCBMx) < 0) return IE(-1);
  if (TidySgOps(TdRefSgOps) != 0) return IE(-1);

  if (HaveCBMx == 0)
    MemCpy(HallCBMx, CBMx, 2);
  else {
    IntSwap(RefCBMx[0].a, RefCBMx[1].a, 12);
    if (CBMx2Multiply(HallCBMx, RefCBMx, CBMx) != 0) return IE(-1);
  }

  if (getBestCBMx(SgOps, SgNumber, TdRefSgOps, HallCBMx) != 0) return IE(-1);

  for (iHS = 0; RefHS[iHS]; iHS++) {
    if (RefHS[iHS] == ' ' && RefHS[iHS + 1] == '(') break;
    if (iHS >= sizeHallSymbol) return IE(-1);
    HallSymbol[iHS] = RefHS[iHS];
  }
  HallSymbol[iHS] = '\0';

  if (MemCmp(&HallCBMx[1], CBMx_1_000, 1) != 0) {
    if (RTMx2XYZ(&HallCBMx[1], CRBF, CTBF, 0, 0, 1, NULL,
                 xyz, sizeof xyz) == NULL) return IE(-1);
    if (sizeHallSymbol < iHS + (int) strlen(xyz) + 4) return IE(-1);
    strcat(HallSymbol, " (");
    strcat(HallSymbol, xyz);
    strcat(HallSymbol, ")");
  }

  return 0;
}
