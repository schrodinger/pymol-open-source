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


static void ShowCBMx(const T_RTMx *CBMx, const T_RTMx *InvCBMx)
{
  fprintf(stdout, "   CBMx = %s\n",
    RTMx2XYZ(   CBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));
  fprintf(stdout, "InvCBMx = %s\n",
    RTMx2XYZ(InvCBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));
}


static int Test_GetRefSetNormAddlG(const T_SgOps *SgOps)
{
  int      SgNumber;
  int      nAddlG, iAddlG;
  T_RTMx    AddlG[3];
  T_SgOps  TidyOps[1], RefSgOps[1], NormSgOps[1], BC_SgOps[1];
  int      iLTr, iInv, iSMx;
  T_RTMx   CBMx[2], SMx[1], LISMx[2];


  SgNumber = GetSpaceGroupType(SgOps, &CBMx[0], &CBMx[1]);
  printf("  SgNumber = %d\n", SgNumber);
  if (SgNumber < 1) return IE(-1);

      nAddlG = GetRefSetNormAddlG(SgNumber, 1, 1, 1, AddlG);
  if (nAddlG < 0)
    return IE(-1);

  ResetSgOps(RefSgOps);
  if (ParseHallSymbol(RefSetHallSymbols[SgNumber],
                      RefSgOps, PHSymOptPedantic) < 0) return IE(-1);
  if (TidySgOps(RefSgOps) != 0) return IE(-1);

  /* check if additional generators leave RefSgOps invariant */
  range1(iAddlG, nAddlG) {
    MemCpy(&LISMx[0], &AddlG[iAddlG], 1);

    if (ChangeBaseFactor(LISMx[0].s.R,    1, LISMx[0].s.R, CRBF, 9) != 0)
      return IE(-1);
    if (ChangeBaseFactor(LISMx[0].s.T, STBF, LISMx[0].s.T, CTBF, 3) != 0)
      return IE(-1);
    if (InverseRTMx(&LISMx[0], &LISMx[1], CRBF) == 0)
      return IE(-1);

    ResetSgOps(BC_SgOps);
    if (CB_SgOps(RefSgOps, &LISMx[0], &LISMx[1], BC_SgOps) != 0)
      return IE(-1);
    if (TidySgOps(BC_SgOps) != 0) return IE(-1);
    if (SgOpsCmp(RefSgOps, BC_SgOps) != 0) return IE(-1);
  }

  /* expand SgOps with transformed additional generators */
  SgOpsCpy(NormSgOps, SgOps);
  range1(iAddlG, nAddlG) {
    if (CB_SMx(SMx, &CBMx[1], &AddlG[iAddlG], &CBMx[0]) != 0) return IE(-1);
    if (ExpSgSMx(NormSgOps, SMx) < 0) return IE(-1);
  }

  /* check that all matrices in NormSgOps leave SgOps invariant */
  SgOpsCpy(TidyOps, SgOps);
  if (TidySgOps(TidyOps) != 0) return IE(-1);
  range1(iLTr, NormSgOps->nLTr)
  range1(iInv, NormSgOps->fInv)
  range1(iSMx, NormSgOps->nSMx)
  {
    SetLISMx(NormSgOps, iLTr, iInv, iSMx, &LISMx[0]);

    if (ChangeBaseFactor(LISMx[0].s.R,    1, LISMx[0].s.R, CRBF, 9) != 0)
      return IE(-1);
    if (ChangeBaseFactor(LISMx[0].s.T, STBF, LISMx[0].s.T, CTBF, 3) != 0)
      return IE(-1);
    if (InverseRTMx(&LISMx[0], &LISMx[1], CRBF) == 0)
      return IE(-1);

    ResetSgOps(BC_SgOps);
    if (CB_SgOps(SgOps, &LISMx[0], &LISMx[1], BC_SgOps) != 0) return IE(-1);
    if (TidySgOps(BC_SgOps) != 0) return IE(-1);
    if (SgOpsCmp(TidyOps, BC_SgOps) != 0) return IE(-1);
  }

  return 0;
}


static int Test_TidyCBMx(const T_SgOps *SgOps)
{
  int      SgNumber;
  T_RTMx   CBMx[2];
  T_SgOps  RefSgOps[1], BC_SgOps[1];

  SgNumber = GetSpaceGroupType(SgOps, &CBMx[0], &CBMx[1]);
  printf("  SgNumber = %d\n", SgNumber);
  if (SgNumber < 1) return IE(-1);
  if (TidyCBMx(SgOps, SgNumber, CBMx) != 0) return IE(-1);
  ShowCBMx(&CBMx[0], &CBMx[1]);

  ResetSgOps(RefSgOps);
  if (ParseHallSymbol(RefSetHallSymbols[SgNumber],
                      RefSgOps, PHSymOptPedantic) < 0) return IE(-1);
  if (TidySgOps(RefSgOps) != 0) return IE(-1);

  ResetSgOps(BC_SgOps);
  if (CB_SgOps(SgOps, &CBMx[0], &CBMx[1], BC_SgOps) != 0)
    return IE(-1);
  if (TidySgOps(BC_SgOps) != 0) return IE(-1);

  if (SgOpsCmp(RefSgOps, BC_SgOps) != 0) return IE(-1);

  return 0;
}


static int Test_BuildHallSymbol(const T_SgOps *SgOps)
{
  int      SgNumber;
  T_RTMx   CBMx[2];
  T_SgOps  TdSgOps[1], HSSgOps[1];
  char     HallSymbol[128];

  SgNumber = GetSpaceGroupType(SgOps, &CBMx[0], &CBMx[1]);
  printf("  SgNumber = %d\n", SgNumber);
  if (SgNumber < 1) return IE(-1);
  if (BuildHallSymbol(SgOps, SgNumber, CBMx,
                      HallSymbol, sizeof HallSymbol) != 0) return IE(-1);
  printf("  %s\n", HallSymbol);

  ResetSgOps(HSSgOps);
  if (ParseHallSymbol(HallSymbol,
                      HSSgOps, PHSymOptPedantic) < 0) return IE(-1);
  if (TidySgOps(HSSgOps) != 0) return IE(-1);

  SgOpsCpy(TdSgOps, SgOps);
  if (TidySgOps(TdSgOps) != 0) return IE(-1);
  if (SgOpsCmp(TdSgOps, HSSgOps) != 0) {
    SgNumber = GetSpaceGroupType(HSSgOps, &CBMx[0], &CBMx[1]);
    printf("  TdSgOps:\n");
    DumpSgOps(TdSgOps, stdout);
    printf("  HSSgNumber = %d\n", SgNumber);
    printf("  HSSgOps:\n");
    DumpSgOps(HSSgOps, stdout);
    return IE(-1);
  }

  return 0;
}


static int Test_Set_ss(const T_SgOps *SgOps)
{
  int     n_ssVM;
  T_ssVM  ssVM[3];

  n_ssVM = Set_ss(SgOps, ssVM);
  if (n_ssVM < 0) return IE(-1);
  printf("  n_ssVM = %d\n", n_ssVM);

  return 0;
}


static int CountRotMxOrder(const int *R)
{
  int  MxA[9], MxB[9];
  int  *RR, *RRR, *Swp, iO;

  const int IdentityMx[] = { 1, 0, 0,
                             0, 1, 0,
                             0, 0, 1 };

  int  nIdentity = 0;

  RR = (int *) R;
  RRR = MxA;

  for (iO = 1; iO < 99; iO++)
  {
    RotMxMultiply(RRR, R, RR);
    if (RR == R) RR = MxB;
    Swp = RR; RR = RRR; RRR = Swp;

    if (MemCmp(IdentityMx, RR, 9) == 0) nIdentity++;

    if (MemCmp(R, RR, 9) == 0)
      break;
  }

  if (nIdentity != 1) return -iO;

  return iO;
}


static int IsFiniteOrderRotMx(const int R[9], int Rtype)
{
  int  iO, i;
  int  ProperR[9], MxA[9], MxB[9];
  int  *RR, *RRR, *Swp;

  const int IdentityMx[] = { 1, 0, 0,
                             0, 1, 0,
                             0, 0, 1 };

  MemCpy(ProperR, R, 9);
  if (Rtype < 0) rangei(9) ProperR[i] *= -1;

  MemCpy(MxA, ProperR, 9);
  RR = MxA;
  RRR = MxB;

  range2(iO, 1, abs(Rtype)) {
    if (MemCmp(IdentityMx, RR, 9) == 0) return 0;
    RotMxMultiply(RRR, ProperR, RR);
    Swp = RR; RR = RRR; RRR = Swp;
  }

  if (MemCmp(IdentityMx, RR, 9) != 0) return 0;

  return 1;
}


static int BuildListRotMx(int Range, int *LVRMx, int mLVRMx,
                          int ProperOnly, int PositiveSenseOnly)
{
  int     Rtype, Order;
  int     nRtypes, nFinite, nList;
  int     R[9];
  T_RMxI  RI[1];


  nRtypes = 0;
  nFinite = 0;
  nList   = 0;

#define loop(i) for (R[i] = -Range; R[i] <= Range; R[i]++)
  loop(0) loop(1) loop(2)
  loop(3) loop(4) loop(5)
  loop(6) loop(7) loop(8)
#undef loop
  {
        Rtype = GetRtype(R);
    if (Rtype == 0) continue;
    nRtypes++;

    if (IsFiniteOrderRotMx(R, Rtype) == 0) continue;
    nFinite++;

    Order = CountRotMxOrder(R);
    if (OrderOfRtype(Rtype) != Order) {
      printf("OrderMismatch %d %d\n",
        OrderOfRtype(Rtype),
        CountRotMxOrder(R));
      return IE(-1);
    }

    if (SetRotMxInfo(R, RI) != Rtype)
      return IE(-1);

    if (   (ProperOnly == 0 || RI->Rtype > 0)
        && (PositiveSenseOnly == 0 || RI->SenseOfRotation >= 0))
    {
      if (nList == mLVRMx) return IE(-1);

      (void) MemCpy(&LVRMx[nList * 9], R, 9);
                           nList++;
    }
  }

  printf("nRtypes=%d\n", nRtypes);
  printf("nFinite=%d\n", nFinite);
  fflush(stdout);

  return nList;
}


static int TestAll(int Range)
{
  int        F_All, F_i, F_j;
  int        mList, nList, iList, jList, i;
  int        *ListRotMx;
  const int  *Ri, *Rj;
  T_RTMx     SMxi[1], SMxj[1];
  int        Rtype[2];
  T_SgOps    SgOps[1];
  int        nGoodComb, nBadComb;
  int        MGC;


  F_All = 1;
  F_i = -1;
  F_j = -1;

  mList = 20000;

      ListRotMx = malloc(mList * 9 * sizeof (*ListRotMx));
  if (ListRotMx == NULL)
    return IE(-1);

  if (F_All)
      nList = BuildListRotMx(Range, ListRotMx, mList, 0, 0);
  else
      nList = BuildListRotMx(Range, ListRotMx, mList, 1, 1);
  if (nList < 0)
    return IE(-1);

  printf("nList = %d\n", nList);

  mList = nList;

      ListRotMx = realloc(ListRotMx, mList * 9 * sizeof (*ListRotMx));
  if (ListRotMx == NULL)
    return IE(-1);

  for (i = 0; i < 3; i++) SMxi->s.T[i] = 0;
  for (i = 0; i < 3; i++) SMxj->s.T[i] = 0;

  nGoodComb = 0;
  nBadComb = 0;

  if (F_i < 0) iList = 0;
  else         iList = F_i;

  for (; iList < nList; iList++)
  {
    Ri = &ListRotMx[iList * 9];
    (void) MemCpy(SMxi->s.R, Ri, 9);

    if (F_j < 0) jList = iList;
    else         jList = F_j;

    for (; jList < nList; jList++)
    {
      Rj = &ListRotMx[jList * 9];
      (void) MemCpy(SMxj->s.R, Rj, 9);

      ResetSgOps(SgOps);

      if (ExpSgSMx(SgOps, SMxi) < 0)
        return IE(-1);

      if (ExpSgSMx(SgOps, SMxj) < 0) {
        nBadComb++;
        ClrSgError();
      }
      else
      {
        nGoodComb++;

        Rtype[0] = GetRtype(Ri);
        Rtype[1] = GetRtype(Rj);
        printf("%d %d %d", Rtype[0], Rtype[1], SgOps->nSMx);
        if (Rtype[0] == 0 || Rtype[1] == 0)
          return IE(-1);

            MGC = GetPG(SgOps);
        if (MGC == MGC_Unknown)
          return IE(-1);

        printf(" %s %s\n", XS_Name[ixXS(MGC)], MG_Names[ixPG(MGC)]);

        if (SgOps->nLTr != 1) return IE(-1);
        if (F_All == 0 && SgOps->fInv != 1) return IE(-1);

        printf("i,jList = %d %d  %s %s\n",
          iList, jList,
          XS_Name[ixXS(MGC)], MG_Names[ixPG(MGC)]);
        fflush(stdout);

#ifdef JUNK
        i = GetSpaceGroupType(SgOps, NULL, NULL);
        printf("SgNumber = %d\n", i);
        if (i < 0) {
#endif
#ifdef JUNK
        {
          int  CutP[3];
          i = GetCutParamMIx(SgOps, 0, CutP);
          if (i == 0) printf("CutP %d %d %d\n", CutP[0], CutP[1], CutP[2]);
        }
        if (i < 0) {
#endif
#ifdef JUNK
        if (Test_GetRefSetNormAddlG(SgOps) != 0) {
#endif
#ifdef JUNK
        if (Test_TidyCBMx(SgOps) != 0) {
#endif
#ifndef JUNK
        if (Test_BuildHallSymbol(SgOps) != 0) {
#endif
#ifdef JUNK
        if (Test_Set_ss(SgOps) != 0) {
#endif
          printf("i,jList = %d %d  %s %s\n",
            iList, jList,
            XS_Name[ixXS(MGC)], MG_Names[ixPG(MGC)]);
          fflush(stdout);
          fprintf(stderr, "%s\n", SgError);
          fflush(stderr);
          ClrSgError();
        }
      }
      if (F_j >= 0) break;
    }
    if (F_i >= 0) break;
  }

  printf("nGoodComb = %d\n", nGoodComb);
  printf("nBadComb  = %d\n", nBadComb);

  free(ListRotMx);

  return 0;
}


int RunSgLiteTests(const char *HallSymbol, const char *Mode, int Range)
{
  T_SgOps  SgOps[1], BufSgOps[1];
  T_RTMx   Z2PCBMx[2];


  if (strcmp(Mode, "TestAll") == 0) {
    if (TestAll(Range) != 0) return IE(-1);
    return 0;
  }

  ResetSgOps(SgOps);
  if (ParseHallSymbol(HallSymbol, SgOps, PHSymOptPedantic) < 0) return IE(-1);

  if (strcmp(Mode, "Primitive") == 0) {
    if (GetZ2PCBMx(SgOps, Z2PCBMx) != 0) return IE(-1);
    printf("  Primitive setting: CBMx = %s\n",
      RTMx2XYZ(Z2PCBMx, CRBF, CTBF, 0, 0, 1, ", ", NULL, 0));

    ResetSgOps(BufSgOps);
    if (CB_SgOps(SgOps, &Z2PCBMx[0], &Z2PCBMx[1], BufSgOps) != 0)
      return IE(-1);
    SgOpsCpy(SgOps, BufSgOps);
    if (SgOps->nLTr != 1) return IE(-1);
  }

#ifdef JUNK
  {
    int  CutP[3];
    if (GetCutParamMIx(SgOps, 0, CutP) != 0) return IE(-1);
    printf("CutP %d %d %d\n", CutP[0], CutP[1], CutP[2]);
  }
#endif
#ifdef JUNK
  if (Test_GetRefSetNormAddlG(SgOps) != 0) return IE(-1);
#endif
#ifdef JUNK
  if (Test_TidyCBMx(SgOps) != 0) return IE(-1);
#endif
#ifndef JUNK
  if (Test_BuildHallSymbol(SgOps) != 0) return IE(-1);
#endif
#ifdef JUNK
  if (Test_Set_ss(SgOps) != 0) return IE(-1);
#endif

  return 0;
}

/* the following code gets rid of a few annoying warnings on GCC
   so that we can watch out for genuine problems... */

static void suppress_compiler_warnings(void)
{
   Test_GetRefSetNormAddlG(NULL);
   Test_TidyCBMx(NULL);
   Test_Set_ss(NULL);
   suppress_compiler_warnings();
}
