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
#include "sgconst.h"


typedef struct
  {
    const int  Sym;
    const int  v[3];
  }
  T_HallTr;


static int IsHSymSpace(int c)
{
  if (c == '\0') return 0;
  if (c == '_') return 1;
  return isspace(c);
}


static int IsHSymChar(int c)
{
  if (c == '\0') return 0;
  return ! IsHSymSpace(c);
}


static int GetAbsOrder(int c)
{
  if (c == '1') return 1;
  if (c == '2') return 2;
  if (c == '3') return 3;
  if (c == '4') return 4;
  if (c == '6') return 6;
  return 0;
}


static int GetScrew(int c)
{
  if (c == '1') return 1;
  if (c == '2') return 2;
  if (c == '3') return 3;
  if (c == '4') return 4;
  if (c == '5') return 5;
  return 0;
}


static int GetRefAxis(int c)
{
  c = tolower(c);

  if (   c == 'x'
      || c == 'y'
      || c == 'z') return c;

  return '\0';
}


static int GetDirCode(int c)
{
  if (   c == '\''
      || c ==  '"'
      || c ==  '*') return c;

  if (   c == ','
      || c == '.') return '\'';

  if (   c == ';'
      || c == ':') return  '"';

  return '\0';
}


static const T_HallTr *GetTr(int Sym)
{
#define T(i) ((i) * (STBF / 12))
#define V(i, j, k) { T(i), T(j), T(k) }

  static const T_HallTr  HallTr[] =
  {
    { 'a', V(6, 0, 0) },
    { 'b', V(0, 6, 0) },
    { 'c', V(0, 0, 6) },
    { 'n', V(6, 6, 6) },
    { 'u', V(3, 0, 0) },
    { 'v', V(0, 3, 0) },
    { 'w', V(0, 0, 3) },
    { 'd', V(3, 3, 3) }
  };

#undef  T
#undef  V

  const int  nHallTr = sizeof HallTr / sizeof (*HallTr);

  int             iHTr;
  const T_HallTr  *HTr;


  Sym = tolower(Sym);

  HTr = HallTr;

  for (iHTr = 0; iHTr < nHallTr; iHTr++, HTr++)
    if (HTr->Sym == Sym) return HTr;

  return NULL;
}


static int GetRMx(int Improper, int AbsOrder,
                  int RefAxis, int DirCode,
                  T_RTMx *SMx)
{
  struct T_TabRMx
  {
    const int  Order;
    const int  DirCode;
    const int  *RMx;
  };

  const struct T_TabRMx  TabRMx[] =
    {
      { 1, '\0', R_1_000 },
      { 2, '\0', R_2_001 },
      { 2, '\'', R_2_1b0 },
      { 2,  '"', R_2_110 },
      { 3, '\0', R_3_001 },
      { 3,  '*', R_3_111 },
      { 4, '\0', R_4_001 },
      { 6, '\0', R_6_001 }
    };

  const int  nTabRMx = (sizeof TabRMx / sizeof (*TabRMx));

  int                    iTRMx, i;
  const struct T_TabRMx  *TRMx;


  TRMx = TabRMx;

  for (iTRMx = 0; iTRMx < nTabRMx; iTRMx++, TRMx++)
    if (TRMx->Order == AbsOrder)
      break;

  for (         ; iTRMx < nTabRMx; iTRMx++, TRMx++)
  {
    if (TRMx->Order != AbsOrder)
      break;

    if (TRMx->DirCode == DirCode)
    {
      if (! Improper)
        for (i = 0; i < 9; i++) SMx->s.R[i] =  TRMx->RMx[i];
      else
        for (i = 0; i < 9; i++) SMx->s.R[i] = -TRMx->RMx[i];

      if      (RefAxis == 'x')
        RotateRotMx(SMx->s.R, R_3_111, R_3i111);
      else if (RefAxis == 'y')
        RotateRotMx(SMx->s.R, R_3i111, R_3_111);

      return 0;
    }
  }

  return -1;
}


static int ParseShortCBO(const char *HSym, int StopChar, int *T, int TBF)
{
  int  iHSym, Row, i, n;

#define cHSym HSym[iHSym]
#define ReturnErr return -(++iHSym)


  iHSym = 0;

  for (Row = 0; Row < 3; Row++)
  {
    while (IsHSymSpace(cHSym)) iHSym++;

    if (Row && cHSym == ',') {
      iHSym++;
      while (IsHSymSpace(cHSym)) iHSym++;
    }

    if (cHSym == '\0' || cHSym == StopChar) ReturnErr;

    i = 1;
    n = sscanf(&cHSym, "%d%n", &T[Row], &i);
    iHSym += (i - 1);
    if (n != 1) ReturnErr;
    iHSym++;

    T[Row] *= (TBF / 12);
  }

#undef cHSym
#undef ReturnErr

  return ++iHSym;
}


int ParseHallSymbolCBMx(const char *HSym, T_SgOps *SgOps, int Options,
                        T_RTMx CBMx[2], int *HaveCBMx)
{
  int      Pedantic, NoCType;
  int      iHSym, nAddedMx, iMxSym, i;
  int      Improper, AbsOrder, Screw;
  int      RefAxis, DirCode;
  int      FirstAbsOrder;
  int      FirstRefAxis;
  T_RTMx   SMx[1];

  const T_HallTr  *HTr;

#define cHSym HSym[iHSym]
#define ReturnErr return -(++iHSym)

  rangei(2) InitRTMx(&CBMx[i], CRBF);
  *HaveCBMx = 0;

  Pedantic = NoCType = 0;

  if (Options & PHSymOptPedantic) Pedantic = 1;
  if (Options & PHSymOptNoCType)  NoCType  = 1;

  iHSym = 0;
  nAddedMx = 0;

  if (! NoCType)
  {
    while (IsHSymSpace(cHSym)) iHSym++;

    if (cHSym == '-') {
      if (ExpSgInv(SgOps, NULL) < 0) ReturnErr;
      iHSym++;
      nAddedMx++;
    }

    if (cHSym == '\0') {
      SetSgError("Error: Lattice type not specified");
      ReturnErr;
    }

        i = ExpSgSymCType(SgOps, cHSym);
    if (i < 0) ReturnErr;
    iHSym++;
    nAddedMx += i;
  }

  i = iHSym;

  while (IsHSymSpace(cHSym)) iHSym++;

  if (cHSym == '\0' || cHSym == '(')
  {
    if (Pedantic) {
      SetSgError("Error: Matrix symbol expected");
      ReturnErr;
    }

    if (cHSym == '\0') return nAddedMx;
  }

  if (! NoCType && Pedantic && iHSym == i) {
    SetSgError("Error: Space expected after lattice type symbol");
    ReturnErr;
  }

  iMxSym = 0;
  FirstAbsOrder = 0;
  FirstRefAxis  = '\0';

  while (cHSym != '\0' && cHSym != '(')
  {
    Improper = AbsOrder = Screw = 0;
    RefAxis = DirCode = '\0';
    for (i = 0; i < 3; i++) SMx->s.T[i] = 0;

    if (cHSym == '-')
    {
      Improper = 1;
      iHSym++;

      if (! IsHSymChar(cHSym)) {
        SetSgError("Error: Incomplete matrix symbol");
        ReturnErr;
      }
    }

          AbsOrder = GetAbsOrder(cHSym);
    if (! AbsOrder) {
      SetSgError("Error: Expected a symbol for rotational order");
      ReturnErr;
    }

    iHSym++;

        Screw = GetScrew(cHSym);
    if (Screw)
    {
      if (Screw >= AbsOrder) {
        SetSgError("Error: Improper screw translation");
        ReturnErr;
      }

      iHSym++;
    }

    while (IsHSymChar(cHSym))
    {
      if (  RefAxis == '\0') {
            RefAxis = GetRefAxis(cHSym);
        if (RefAxis != '\0')
        {
          if (    AbsOrder == 1
              || (AbsOrder == 3 && DirCode == '*')) {
            SetSgError("Error: Inconsistent matrix symbol");
            ReturnErr;
          }

          iHSym++;
          continue;
        }
      }
      else if (GetRefAxis(cHSym) != '\0') {
        SetSgError("Error: Multiple axis symbols");
        ReturnErr;
      }

      if (  DirCode == '\0') {
            DirCode = GetDirCode(cHSym);
        if (DirCode != '\0')
        {
          if (   ! (AbsOrder == 2 && (   DirCode ==  '"'
                                      || DirCode == '\''))
              && ! (AbsOrder == 3 && DirCode == '*'))
          {
            SetSgError("Error: Inconsistent matrix symbol");
            ReturnErr;
          }

          if (Screw) {
            SetSgError("Error: Screw translation for non-principal direction");
            ReturnErr;
          }

          iHSym++;
          continue;
        }
      }
      else if (GetDirCode(cHSym) != '\0') {
        SetSgError("Error: Multiple axis symbols");
        ReturnErr;
      }

          HTr = GetTr(cHSym);
      if (HTr)
      {
        for (i = 0; i < 3; i++)
          SMx->s.T[i] = (SMx->s.T[i] + HTr->v[i]) % STBF;

        iHSym++;
        continue;
      }

      if (cHSym == '(')
      {
        if (Pedantic) {
          SetSgError("Error: Space expected before change-of-basis operator");
          ReturnErr;
        }

        break;
      }

      SetSgError("Error: Malformed matrix symbol");
      ReturnErr;
    }

    if (RefAxis == '\0')
    {
      if      (iMxSym == 0)
      {
        if (      AbsOrder != 1
            && ! (AbsOrder == 3 && DirCode == '*'))
          RefAxis = 'z';
      }
      else if (iMxSym == 1)
      {
        if      (AbsOrder == 2)
        {
          if      (FirstAbsOrder == 2 || FirstAbsOrder == 4)
          {
            if (DirCode == '\0')
              RefAxis = 'x';
          }
          else if (FirstAbsOrder == 3 || FirstAbsOrder == 6)
          {
            if (DirCode == '\0')
              DirCode = '\'';

            RefAxis = FirstRefAxis;
          }
        }
        else if (   AbsOrder == 3
                 && (FirstAbsOrder == 2 || FirstAbsOrder == 4)
                 && DirCode == '\0')
          DirCode = '*';
      }
      else if (iMxSym == 2)
      {
        if (AbsOrder == 3 && DirCode == '\0')
          DirCode = '*';
      }
    }

    if (RefAxis == '\0' && (   DirCode ==  '"'
                            || DirCode == '\''))
      RefAxis = 'z';

    if (RefAxis == '\0' && AbsOrder != 1 && DirCode != '*') {
      SetSgError("Error: Need explicit axis symbol");
      ReturnErr;
    }

    if (GetRMx(Improper, AbsOrder, RefAxis, DirCode, SMx) < 0) {
      SetSgError("Internal Error: GetRMx() failed");
      ReturnErr;
    }

    if (Screw)
    {
      switch (RefAxis)
      {
        case 'x': i = 0; break;
        case 'y': i = 1; break;
        default:  i = 2; break;
      }

      SMx->s.T[i] += STBF * Screw / AbsOrder;
    }

    if (ExpSgSMx(SgOps, SMx) < 0)
      ReturnErr;

    if (iMxSym == 0) {
      FirstAbsOrder = AbsOrder;
      FirstRefAxis  = RefAxis;
    }

    iMxSym++;

    if (Improper || AbsOrder != 1)
      nAddedMx++;

    while (IsHSymSpace(cHSym)) iHSym++;
  }

  if (cHSym == '(')
  {
    iHSym++;

        i = ParseShortCBO(&cHSym, ')', CBMx[0].s.T, CTBF);
    if (i <= 0) {
          i = ParseStrXYZ(&cHSym, ')', &CBMx[0], CRBF, CTBF);
      if (i < 0) {
        iHSym += -i - 1;
        SetSgError("Error: Malformed change-of-basis operator");
        ReturnErr;
      }
    }

    iHSym += i - 1;

    while (IsHSymSpace(cHSym)) iHSym++;

    if (cHSym != ')') {
      SetSgError(
        "Error: Closing parenthesis expected after change-of-basis operator");
      ReturnErr;
    }

    if (InverseRTMx(&CBMx[0], &CBMx[1], CRBF) == 0) {
      SetSgError("Error: Change-of-basis operator is not invertible");
      ReturnErr;
    }

    iHSym++;
    *HaveCBMx = -iHSym;
  }

  while (IsHSymSpace(cHSym)) iHSym++;

  if (cHSym != '\0') {
    SetSgError("Error: Unexpected extra character");
    ReturnErr;
  }

#undef cHSym
#undef ReturnErr

  return nAddedMx;
}


int ParseHallSymbol(const char *HSym, T_SgOps *SgOps, int Options)
{
  int      status, HaveCBMx;
  T_RTMx   CBMx[2];
  T_SgOps  LocSgOps[2];


  if (SgOps)
    SgOpsCpy(LocSgOps, SgOps);
  else
    ResetSgOps(LocSgOps);

      status = ParseHallSymbolCBMx(HSym, LocSgOps, Options, CBMx, &HaveCBMx);
  if (status < 0) return status;

  if (HaveCBMx != 0) {
    if (SgOps == NULL) SgOps = &LocSgOps[2];
    ResetSgOps(SgOps);
    SgOps->NoExpand = LocSgOps->NoExpand;
    if (CB_SgOps(LocSgOps, &CBMx[0], &CBMx[1], SgOps) != 0)
      return HaveCBMx;
  }
  else if (SgOps)
    SgOpsCpy(SgOps, LocSgOps);

  return status;
}
