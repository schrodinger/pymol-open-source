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


int GetRefSetNormAddlG(int SgNumber, int affine, int UseK2L, int UseL2N,
                       T_RTMx *AddlG)
{
  int         nAddlG, iType, nAddedMx, i;
  const char  *HMxSym;
  T_SgOps     SgOps[1];


  if (SgNumber < 1 || SgNumber > 230) return IE(-1);

  nAddlG = 0;

  range1(iType, 2)
  {
    HMxSym = NULL;

    if      (iType == 0 && UseK2L)
      HMxSym = RefSetNormAddlG[SgNumber].K2L;
    else if (iType == 1 && UseL2N && (SgNumber >= 75 || affine))
      HMxSym = RefSetNormAddlG[SgNumber].L2N;

    if (HMxSym)
    {
      ResetSgOps(SgOps);
      SgOps->NoExpand = 1;

             nAddedMx = ParseHallSymbol(HMxSym, SgOps, PHSymOptNoCType);
      if (   nAddedMx < 1
          || SgOps->nLTr != 1
          || SgOps->fInv - 1 + SgOps->nSMx - 1 + nAddlG > 3)
        return IE(-1);

      if (SgOps->fInv == 2) {
        InitRotMx(AddlG[nAddlG].s.R, -1);
        rangei(3) AddlG[nAddlG].s.T[i] = SgOps->InvT[i];
        nAddlG++;
      }

      if (SgOps->nSMx > 1)
        MemCpy(&AddlG[nAddlG], &SgOps->SMx[1], SgOps->nSMx - 1);

      nAddlG += SgOps->nSMx - 1;
    }
  }

  return nAddlG;
}


int CheckMonoRefSetAffNormRestrictions(int SgNumber, const int M[9], int BF)
{
  /* International Tables Volume A, chapter 15, tables 15.3.3 & 15.3.4.
   */
  switch (SgNumber) {
    case  3:
    case  4:
    case  6:
    case 10:
    case 11: /* M2 */
      break;

    case  5:
    case  8:
    case 12: /* M4 */
    case  9:
    case 15: /* M6 or M12 */
      if (M[0] % (2 * BF) == 0) return -1;
      if (M[6] % (2 * BF) != 0) return -1;
      if (M[8] % (2 * BF) == 0) return -1;
      break;

    case  7:
    case 13:
    case 14: /* M5 */
      if (M[0] % (2 * BF) == 0) return -1;
      if (M[2] % (2 * BF) != 0) return -1;
      if (M[8] % (2 * BF) == 0) return -1;
      break;
  }

  return 0;
}
