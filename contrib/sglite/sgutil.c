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


int *IntSwap(int *a, int *b, int n)
{
  int  *a0, i;

  for (a0 = a; n--; a++, b++) {
    i = *a; *a = *b; *b = i;
  }

  return a0;
}

void IntSetZero(int *a, int n)
{
  int i; rangei(n) a[i] = 0;
}

int IntIsZero(const int *a, int n)
{
  while (n--) if (a[n]) return 0;
  return 1;
}


int iModPositive(int ix, int iy)
{
  if (iy > 0)
  {
    ix %= iy;
    if (ix < 0) ix += iy;
  }

  return ix;
}

int iModShort(int ix, int iy)
{
      ix = iModPositive(ix, iy);
  if (ix > iy / 2)
      ix -= iy;

  return ix;
}

void ViModShort(int *ix, int n, int iy)
{
  int i; for (i = 0; i < n; i++) ix[i] = iModShort(ix[i], iy);
}

void ViModPositive(int *ix, int n, int iy)
{
  int i; for (i = 0; i < n; i++) ix[i] = iModPositive(ix[i], iy);
}


T_RTMx *SetLISMx(const T_SgOps *SgOps, int iLTr, int iInv, int iSMx,
                 T_RTMx *LISMx)
{
  int  i;

  MemCpy(LISMx, &SgOps->SMx[iSMx], 1);
  if (iInv) SMx_t_InvT(LISMx, SgOps->InvT, LISMx);
  rangei(3) LISMx->s.T[i] += SgOps->LTr[iLTr].v[i];

  return LISMx;
}


int Discretize(double fVal, int *iVal, int Fac)
{
  if (Fac == 0) return -1;
      fVal *= Fac;
  if (fVal < 0.) (*iVal) = (int)(fVal - .5);
  else           (*iVal) = (int)(fVal + .5);
      fVal -= (*iVal);
      fVal /= Fac;
  if (fVal < 0.) fVal = -fVal;
  if (fVal > .0001) return -1;
  return 0;
}


int ChangeBaseFactor(const int *Old, int OldBF, int *New, int NewBF, int n)
{
  int  i;

  rangei(n) {
        New[i] = Old[i] * NewBF;
    if (New[i] %  OldBF) return -1;
        New[i] /= OldBF;
  }

  return 0;
}


int SignHemisphere(int h, int k, int l)
{
  if (l >  0) return  1;
  if (l == 0) {
    if (k >  0) return  1;
    if (k == 0) {
      if (h >  0) return  1;
      if (h == 0)
        return 0;
    }
  }

  return -1;
}


void SetRminusI(const int *R, int *RmI, int Inv)
{
  int  i;

  if (Inv == 0)
    for (i = 0; i < 9; i++)
      RmI[i] =  R[i];
  else
    for (i = 0; i < 9; i++)
      RmI[i] = -R[i];

  for (i = 0; i < 9; i += 4)
    RmI[i] -= 1;
}


int NextOf_n_from_m(int m, int n, int *ix)
{
  int  p, l;

  p = l = n - 1;

  for (; p >= 0;)
  {
        ix[p]++;
    if (ix[p] == m - l + p)
      p--;
    else if (p < l)
    {
      ix[p + 1] = ix[p];
         p++;
    }
    else
      return 1;
  }

  return 0;
}


void SgOpsCpy(T_SgOps *t, const T_SgOps *s)
{
  t->NoExpand = s->NoExpand;
  t->nLSL = s->nLSL;
  t->nSSL = s->nSSL;
  t->nLTr = s->nLTr;
  t->fInv = s->fInv;
  t->nSMx = s->nSMx;
  MemCpy(t->LTr, s->LTr, SgOps_mLTr);
  MemCpy(t->InvT, s->InvT, 3);
  MemCpy(t->SMx, s->SMx, SgOps_mSMx);
}


int SgOpsCmp(const T_SgOps *t, const T_SgOps *s)
{
  int  c;

  if (t->NoExpand < s->NoExpand) return -1;
  if (t->NoExpand > s->NoExpand) return  1;

  if (t->nLSL < s->nLSL) return -1;
  if (t->nLSL > s->nLSL) return  1;

  if (t->nSSL < s->nSSL) return -1;
  if (t->nSSL > s->nSSL) return  1;

  if (t->nLTr < s->nLTr) return -1;
  if (t->nLTr > s->nLTr) return  1;

  if (t->fInv < s->fInv) return -1;
  if (t->fInv > s->fInv) return  1;

  if (t->nSMx < s->nSMx) return -1;
  if (t->nSMx > s->nSMx) return  1;

  c = MemCmp(t->LTr, s->LTr, SgOps_mLTr);
  if (c) return c;

  c = MemCmp(t->InvT, s->InvT, 3);
  if (c) return c;

  return MemCmp(t->SMx, s->SMx, SgOps_mSMx);
}
