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


static int Show_wI_Tr(const T_RTMx *SMx, FILE *fp)
{
  const char  *ff;
  T_RMxI      RI[1];
  int         wI[3], Tr[3], i;


  if (SetRotMxInfo(SMx->s.R, RI) == 0)
    return IE(-1);

  if (Set_wI_Tr(SMx->a, NULL, RI, wI, Tr) != 0)
    return IE(-1);

  fprintf(fp, " [%2d %2d %2d] %d",
    RI->EV[0], RI->EV[1], RI->EV[2],
    RI->Rtype);

  if (RI->SenseOfRotation < 0)
    fputs("^-1", fp);
  else
    fputs("   ", fp);

  fputs(" (", fp);

  rangei(3) {
        ff = FormatFraction(wI[i], STBF, 0, NULL, 0);
    if (ff == NULL) return IE(-1);
    if (i) putc(' ', fp);
    fputs(ff, fp);
  }

  fputs(")", fp);

  rangei(3) {
        ff = FormatFraction(Tr[i], CTBF, 0, NULL, 0);
    if (ff == NULL) return IE(-1);
    putc(' ', fp);
    fputs(ff, fp);
  }

  return 0;
}


int DumpSgOps(const T_SgOps *SgOps, FILE *fp)
{
  int         iLTr, iSMx, i;
  T_RTMx      SMx[1];
  const char  *xyz;


  fprintf(fp, "+ nLTr %d\n", SgOps->nLTr);

  fprintf(fp, "+ fInv %d (%d %d %d)", SgOps->fInv,
    SgOps->InvT[0],
    SgOps->InvT[1],
    SgOps->InvT[2]);
  if (SgOps->fInv == 2) {
    InitRTMx(SMx, -1);
    rangei(3) SMx->s.T[i] = SgOps->InvT[i];
    xyz = RTMx2XYZ(SMx, 1, STBF, 0, 0, 1, NULL, NULL, 0);
    if (! xyz) return IE(-1);
    fprintf(fp, " %s", xyz);
  }
  putc('\n', fp);

  fprintf(fp, "+ nSMx %d\n", SgOps->nSMx);

  range1(iLTr, SgOps->nLTr)
    fprintf(fp, "+ LTr[%d] (%d %d %d)\n", iLTr,
      SgOps->LTr[iLTr].v[0],
      SgOps->LTr[iLTr].v[1],
      SgOps->LTr[iLTr].v[2]);

  range1(iSMx, SgOps->nSMx)
  {
    fprintf(fp, "+ SMx[%02d] ", iSMx);

    xyz = RTMx2XYZ(&SgOps->SMx[iSMx], 1, STBF, 0, 0, 1, NULL, NULL, 0);
    if (! xyz) return IE(-1);
    fprintf(fp, " %-26s", xyz);

    if (Show_wI_Tr(&SgOps->SMx[iSMx], fp) != 0) return -1;

    putc('\n', fp);
  }

  return 0;
}
