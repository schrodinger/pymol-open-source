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


#define HmulR(Hm, H, R) \
  (Hm[0]) = (R)[0] * (H[0]) + (R)[3] * (H[1]) + (R)[6] * (H[2]); \
  (Hm[1]) = (R)[1] * (H[0]) + (R)[4] * (H[1]) + (R)[7] * (H[2]); \
  (Hm[2]) = (R)[2] * (H[0]) + (R)[5] * (H[1]) + (R)[8] * (H[2])


static int Is000(const int H[3])
{
  int  i;
  rangei(3) if (H[i]) return 0;
  return 1;
}


static int AreSameMIx(const int H1[3], const int H2[3])
{
  int  i;
  rangei(3) if (H1[i] != H2[i]) return 0;
  return 1;
}


static int AreFriedelMates(const int H1[3], const int H2[3])
{
  int  i;
  rangei(3) if (H1[i] != -H2[i]) return 0;
  return 1;
}


int CmpEqMIx(const int H1[3], const int H2[3])
{
  int        i;
  const int  P[3] = { 2, 0, 1 };

  rangei(3) {
    if (H1[P[i]] >= 0 && H2[P[i]] <  0) return -1;
    if (H1[P[i]] <  0 && H2[P[i]] >= 0) return  1;
  }

  rangei(3) {
    if (abs(H1[P[i]]) < abs(H2[P[i]])) return -1;
    if (abs(H1[P[i]]) > abs(H2[P[i]])) return  1;
  }

  return 0;
}


int IsSysAbsMIx(const T_SgOps *SgOps, const int H[3], int *TH_Restriction)
{
  int        iSMx, iLTr, i;
  int        Hm[3];
  int        TH, InvTmT[3];
  const int  *T, *TS, *TR;

  /* systematically absent reflection: if HR == H and HT != 0 mod 1
     restricted phase: if HR == -H: phi(H) = pi*HT + n*pi
   */

  if (TH_Restriction) *TH_Restriction = -1;

  range1(iSMx, SgOps->nSMx)
  {
    HmulR(Hm, H, SgOps->SMx[iSMx].s.R);
             T = SgOps->SMx[iSMx].s.T;
    TS = NULL;
    TR = NULL;

    if      (AreSameMIx(H, Hm)) {
      TS = T;
      if (TH_Restriction && SgOps->fInv == 2) {
        rangei(3) InvTmT[i] = SgOps->InvT[i] - T[i];
        TR = InvTmT;
      }
    }
    else if (AreFriedelMates(H, Hm)) {
      if (TH_Restriction) TR = T;
      if (SgOps->fInv == 2) {
        rangei(3) InvTmT[i] = SgOps->InvT[i] - T[i];
        TS = InvTmT;
      }
    }

    if (TS) {
      range1(iLTr, SgOps->nLTr) {
        TH = 0;
        rangei(3) TH += (TS[i] + SgOps->LTr[iLTr].v[i]) * H[i];
        TH %= STBF;
        if (TH != 0) return 1;
      }
    }

    if (TR) {
      range1(iLTr, SgOps->nLTr) {
        TH = 0;
        rangei(3) TH += (TR[i] + SgOps->LTr[iLTr].v[i]) * H[i];
        TH %= STBF; if (TH < 0) TH += STBF;
        if      (*TH_Restriction < 0) *TH_Restriction = TH;
        else if (*TH_Restriction != TH) return 1;
      }
    }
  }

  return 0;
}


int IsCentricMIx(const T_SgOps *SgOps, const int H[3])
{
  int  iSMx, Hm[3];

  if (SgOps->fInv == 2) return 1;

  range1(iSMx, SgOps->nSMx) {
    HmulR(Hm, H, SgOps->SMx[iSMx].s.R);
    if (AreFriedelMates(H, Hm)) return 1;
  }

  return 0;
}


int GetPhaseRestriction(const T_SgOps *SgOps, const int H[3])
{
  int        TH, iSMx, Hm[3], i;
  const int  *TR;

  /* restricted phase: if HR == -H: phi(H) = pi*HT + n*pi
   */

  TH = -1;
  TR = NULL;

  if (SgOps->fInv == 2) {
    TR = SgOps->InvT;
  }
  else {
    range1(iSMx, SgOps->nSMx) {
      HmulR(Hm, H, SgOps->SMx[iSMx].s.R);
      if (AreFriedelMates(H, Hm)) {
        TR = SgOps->SMx[iSMx].s.T;
        break;
      }
    }
  }

  if (TR) {
    TH = 0;
    rangei(3) TH += TR[i] * H[i];
    TH %= STBF; if (TH < 0) TH += STBF;
  }

  return TH;
}


int EpsilonMIx(const T_SgOps *SgOps, const int H[3])
{
  int  Epsilon;
  int  iSMx, Hm[3];

  Epsilon = 0;

  range1(iSMx, SgOps->nSMx) {
    HmulR(Hm, H, SgOps->SMx[iSMx].s.R);
    if (AreSameMIx(H, Hm) || (SgOps->fInv == 2 && AreFriedelMates(H, Hm)))
      Epsilon++;
  }

  if (Epsilon == 0 || SgOps->nSMx % Epsilon != 0) return IE(-1);

  return Epsilon;
}


int MultMIx(const T_SgOps *SgOps, int FriedelSym, const int H[3])
{
  int  Centro;
  int  M, R;
  int  iSMx, Hm[3];

  if (Is000(H)) return 1;

  Centro = ((SgOps->fInv == 2) || FriedelSym);

  M = 0;
  R = 0;

  range1(iSMx, SgOps->nSMx) {
    HmulR(Hm, H, SgOps->SMx[iSMx].s.R);
    if      (AreSameMIx(H, Hm)) M++;
    else if (AreFriedelMates(H, Hm)) R++;
  }

  if (M == 0 || SgOps->nSMx % M != 0 || (R != 0 && R != M)) return IE(-1);

  M = SgOps->nSMx / M;

  if (Centro && R == 0)
    M *= 2;

  return M;
}


int BuildEqMIx(const T_SgOps *SgOps, int FriedelSym, const int H[3],
               T_EqMIx *EqMIx)
{
  int      iSMx, Hm[3], iEq, i;
  T_EqMIx  BufEqMIx[1];

  if (EqMIx == NULL)
      EqMIx = BufEqMIx;

  EqMIx->fInv = 1;
  if ((SgOps->fInv == 2 || FriedelSym) && ! Is000(H))
    EqMIx->fInv = 2;

  EqMIx->N = 0;

  range1(iSMx, SgOps->nSMx)
  {
    HmulR(Hm, H, SgOps->SMx[iSMx].s.R);

    range1(iEq, EqMIx->N) {
      if (AreSameMIx(Hm, EqMIx->H[iEq])) break;
      if (EqMIx->fInv != 2) continue;
      if (AreFriedelMates(Hm, EqMIx->H[iEq])) break;
    }

    if (iEq == EqMIx->N)
    {
      if (EqMIx->N >= sizeof EqMIx->TH / sizeof (*EqMIx->TH))
        return IE(0);

      MemCpy(EqMIx->H[iEq], Hm, 3);

          EqMIx->TH[iEq] = 0;
      rangei(3)
          EqMIx->TH[iEq] += SgOps->SMx[iSMx].s.T[i] * H[i];
          EqMIx->TH[iEq] %= STBF;
      if (EqMIx->TH[iEq] < 0)
          EqMIx->TH[iEq] += STBF;

      EqMIx->N++;
    }
  }

  if (SgOps->nSMx % EqMIx->N) return IE(0); /* another error trap */

  return EqMIx->fInv * EqMIx->N;
}


static int OneMxCutPRange(const int R[9])
{
  int  m, ir, ic, s;

  m = 0;

  range1(ic, 3) {
    s = 0;
    range1(ir, 3) s += abs(R[ir * 3 + ic]);
    if (m < s)
        m = s;
  }

  return m + 1;
}


static int SetCheckCutPRange(const T_SgOps *SgOps)
{
  int  Range, iSMx, m;

  Range = 0;

  range1(iSMx, SgOps->nSMx) {
                m = OneMxCutPRange(SgOps->SMx[iSMx].s.R);
    if (Range < m)
        Range = m;
  }

  return Range;
}


static int CheckCutParam(const T_SgOps *SgOps, int FriedelSym, int CutP[3],
                         int Range, int FullBlock)
{
  int      iBV, iEq, i;
  int      AdjRange[3], Step[3], H[3];
  T_EqMIx  EqMIx[1];

  rangei(3) AdjRange[i] = Range;

  range1(iBV, 3)
  {
    rangei(3) Step[i] = 1;
    if (FullBlock == 0) Step[iBV] = 2 * Range;

    for (H[0] = -AdjRange[0]; H[0] <= AdjRange[0]; H[0] += Step[0])
    for (H[1] = -AdjRange[1]; H[1] <= AdjRange[1]; H[1] += Step[1])
    for (H[2] = -AdjRange[2]; H[2] <= AdjRange[2]; H[2] += Step[2])
    {
      if (BuildEqMIx(SgOps, FriedelSym, H, EqMIx) < 1)
        return IE(-1);

      /* search for equivalent hkl in an active octant */
      range1(iEq, EqMIx->N) {
        rangei(3) if (CutP[i] == 0 && EqMIx->H[iEq][i] < 0) break;
        if (i == 3) break;
        if (EqMIx->fInv != 2) continue;
        rangei(3) if (CutP[i] == 0 && EqMIx->H[iEq][i] > 0) break;
        if (i == 3) break;
      }

      if (iEq == EqMIx->N) return 0; /* CutParam does not work */
    }

    if (FullBlock != 0) break;

    AdjRange[iBV]--;
  }

  return 1; /* CutParam works fine */
}


int GetCutParamMIx(const T_SgOps *SgOps, int FriedelSym, int CutP[3])
{
  int  Range, iTrial, status, i;

  static int ListTrialCutP[7][3] =
    {
      {  0,  0,  0 },
      {  0, -1,  0 },
      { -1,  0,  0 },
      {  0,  0, -1 },
      { -1, -1,  0 },
      {  0, -1, -1 },
      { -1,  0, -1 }
    };

  Range = SetCheckCutPRange(SgOps);
#ifdef JUNK
printf("Range=%d\n", Range);
#endif

  range1(iTrial, 7) {
        status = CheckCutParam(SgOps, FriedelSym, ListTrialCutP[iTrial],
                               Range, 0);
    if (status < 0) return IE(-1);
#ifdef JUNK
{
  int stat2;
  stat2 = CheckCutParam(SgOps, FriedelSym, ListTrialCutP[iTrial], Range + 2, 1);
  if (stat2 != status) {
    printf("CutP %d %d %d\n",
      ListTrialCutP[iTrial][0],
      ListTrialCutP[iTrial][1],
      ListTrialCutP[iTrial][2]);
    DumpSgOps(SgOps, stdout);
    return IE(-1);
  }
}
#endif
    if (status > 0) {
      MemCpy(CutP, ListTrialCutP[iTrial], 3);
      return 0;
    }
  }

  rangei(3) CutP[i] = -1;

  return 0;
}


static int IsInActiveArea(const int CutP[3], const int H[3])
{
  int  i;

  rangei(3) if (CutP[i] == 0 && H[i] < 0) return 0;
  return 1;
}


int GetMasterMIx(const T_EqMIx *EqMIx, const int CutP[3], int MasterH[3])
{
  int  iEq, iInv, i;
  int  HaveMaster, EqH[3];

  HaveMaster = 0;

  range1(iEq, EqMIx->N) {
    MemCpy(EqH, EqMIx->H[iEq], 3);
    range1(iInv, EqMIx->fInv) {
      if (iInv) rangei(3) EqH[i] *= -1;
      if (IsInActiveArea(CutP, EqH)) {
        if (HaveMaster == 0 || CmpEqMIx(MasterH, EqH) > 0) {
          MemCpy(MasterH, EqH, 3);
          HaveMaster = 1;
        }
      }
    }
  }

  if (HaveMaster == 0) return IE(-1);

  return 0;
}


int GetMasterMIx_and_MateID(const T_SgOps *SgOps,
                            const int CutP[3], const int MIx[3],
                            int MasterMIx[3], int *MateID)
{
  int      MateMIx[3], MateMasterMIx[3], i;
  T_EqMIx  EqMIx[1];

  if (BuildEqMIx(SgOps, 0, MIx, EqMIx) == 0) return IE(-1);
  if (GetMasterMIx(EqMIx, CutP, MasterMIx) != 0) return IE(-1);
  *MateID = 0;
  if (SgOps->fInv == 1) {
    rangei(3) MateMIx[i] = -MIx[i];
    if (BuildEqMIx(SgOps, 0, MateMIx, EqMIx) == 0) return IE(-1);
    if (GetMasterMIx(EqMIx, CutP, MateMasterMIx) != 0) return IE(-1);
    if (CmpEqMIx(MasterMIx, MateMasterMIx) > 0) {
      rangei(3) MasterMIx[i] = MateMasterMIx[i];
      *MateID = 1;
    }
  }

  return 0;
}
