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


int isChiralSpaceGroup(const T_SgOps *SgOps)
{
  int  iSMx;

  if (SgOps->fInv == 2) return 0;

  range2(iSMx, 1, SgOps->nSMx)
    if (GetRtype(SgOps->SMx[iSMx].s.R) < 0) return 0;

  return 1;
}


int isEnantiomorphicSpaceGroup(const T_SgOps *SgOps)
{
  int      T1, T2;
  T_RTMx   FlipCBMx[2];
  T_SgOps  FlipSgOps[1];


  InitRTMx(FlipCBMx, -CRBF);

  ResetSgOps(FlipSgOps);
  if (CB_SgOps(SgOps, FlipCBMx, FlipCBMx, FlipSgOps) != 0) return IE(-1);

      T1 = GetSpaceGroupType(SgOps, NULL, NULL);
  if (T1 < 1) return IE(-1);
      T2 = GetSpaceGroupType(FlipSgOps, NULL, NULL);
  if (T2 < 1) return IE(-1);
  if (T1 != T2) return T2;
  return 0;
}
