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


int iGCD(const int a, const int b)
{
  int  ri, rj, rk;


      ri = a;
  if (ri < 0) ri = -ri;

  if ((rj = b) != 0)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    if (ri < 0) ri = -ri;
  }

  return ri;
}


int FindGCD(const int *S, int nS)
{
  int  ri, rj, rk;


  if (nS-- == 0) return 0;

      ri = *S++;
  if (ri < 0) ri = -ri;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      if (ri < 0) ri = -ri;

      if (ri == 1) break;
    }
  }

  return ri;
}


int CancelGCD(int *S, int nS)
{
  int  GCD, i;


      GCD = FindGCD(S, nS);
  if (GCD)
    rangei(3) S[i] /= GCD;

  return GCD;
}


int CancelBFGCD(int *S, int nS, int BF)
{
  int  GCD, i;

      GCD = FindGCD(S, nS);
      GCD = iGCD(GCD, BF);
  if (GCD) {
    rangei(3) S[i] /= GCD;
    return BF / GCD;
  }

  return GCD;
}


int iLCM(const int a, const int b)
{
  int  al, ri, rj, rk;


      al = a;
  if (al == 0) al = 1;

       ri = al;
  if ((rj = b) != 0)
  {
    for (;;)
    {
      rk = ri % rj; if (rk == 0) { ri = rj; break; }
      ri = rj % rk; if (ri == 0) { ri = rk; break; }
      rj = rk % ri; if (rj == 0) {          break; }
    }

    ri = al / ri * b;
  }

  if (ri < 0) return -ri;
              return  ri;
}


int FindLCM(const int *S, int nS)
{
  int  a, b, ri, rj, rk;


  if (nS-- == 0) return 1;

      ri = *S++;
  if (ri == 0) ri = 1;

  while (nS--)
  {
    if ((rj = *S++) != 0)
    {
      a = ri; b = rj;

      for (;;)
      {
        rk = ri % rj; if (rk == 0) { ri = rj; break; }
        ri = rj % rk; if (ri == 0) { ri = rk; break; }
        rj = rk % ri; if (rj == 0) {          break; }
      }

      ri = a / ri * b;
    }
  }

  if (ri < 0) return -ri;
              return  ri;
}


void SimplifyFraction(int nume, int deno, int *o_nume, int *o_deno)
{
  int gcd = iGCD(nume, deno);
  if (gcd)
  {
    *o_nume = nume / gcd;
    *o_deno = deno / gcd;

    if (*o_deno < 0) {
      *o_nume *= -1;
      *o_deno *= -1;
    }
  }
}


#define GmulV(GV, G, V) \
  (GV)[0] = (G)[0] * (V)[0] + (G)[1] * (V)[1] + (G)[2] * (V)[2]; \
  (GV)[1] = (G)[3] * (V)[0] + (G)[4] * (V)[1] + (G)[5] * (V)[2]; \
  (GV)[2] = (G)[6] * (V)[0] + (G)[7] * (V)[1] + (G)[8] * (V)[2]


int iScalProd(const int *u, const int *v, const int *PseudoG)
{
  int  Prod, Gv[3];


  if (PseudoG)
  {
    GmulV(Gv, PseudoG, v);
    v = Gv;
  }

  Prod =   u[0] * v[0]
         + u[1] * v[1]
         + u[2] * v[2];

  return Prod;
}


void iCrossProd(int *rxs, const int *r, const int *s, const int *PseudoG)
{
  int  Gr[3], Gs[3];


  if (PseudoG)
  {
    GmulV(Gr, PseudoG, r);
    GmulV(Gs, PseudoG, s);
    r = Gr;
    s = Gs;
  }

  rxs[0] = r[1] * s[2] - r[2] * s[1];
  rxs[1] = r[2] * s[0] - r[0] * s[2];
  rxs[2] = r[0] * s[1] - r[1] * s[0];
}


#undef GmulV


int AreLinDepV(const int a[3], const int b[3])
{
  int        axb[3], i;
  const int  V000[3] = { 0, 0, 0 };


  iCrossProd(axb, a, b, NULL);
  if (MemCmp(axb, V000, 3) == 0)
  {
    rangei(3) {
      if (a[i]) {
        if (abs(a[i]) > abs(b[i])) return 1;
        return -1;
      }
    }
  }

  return 0;
}


void SMx_t_InvT(const T_RTMx *SMx, const int InvT[3], T_RTMx *ProdSMx)
{
  int   i;

  rangei(9) ProdSMx->s.R[i] = -SMx->s.R[i];
  rangei(3) ProdSMx->s.T[i] = -SMx->s.T[i] + InvT[i];
}


void RotMx_t_Vector(int *R_t_V, const int *RotMx, const int *Vector, int FacTr)
{
  const int  *vec;


  if (FacTr > 0)
  {
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec;
    *R_t_V   %= FacTr; if (*R_t_V < 0) *R_t_V += FacTr;
     R_t_V++;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec;
    *R_t_V   %= FacTr; if (*R_t_V < 0) *R_t_V += FacTr;
     R_t_V++;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx   * *vec;
    *R_t_V   %= FacTr; if (*R_t_V < 0) *R_t_V += FacTr;
  }
  else
  {
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V++ += *RotMx++ * *vec;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V++ += *RotMx++ * *vec;
                            vec = Vector;
    *R_t_V   =  *RotMx++ * *vec++;
    *R_t_V   += *RotMx++ * *vec++;
    *R_t_V   += *RotMx   * *vec;
  }
}


void RotMxMultiply(int *rmxab, const int *rmxa, const int *rmxb)
{
  const int  *a, *b;

  /* no loops to be as fast as posslible */

  a = rmxa;
  b = rmxb;
  *rmxab  = *a++ * *b; b += 3; /* r11 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r12 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r13 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a++ * *b; b = rmxb;
   rmxab++;

  rmxa = a;
  *rmxab  = *a++ * *b; b += 3; /* r21 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r22 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r23 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a++ * *b; b = rmxb;
   rmxab++;

  rmxa = a;
  *rmxab  = *a++ * *b; b += 3; /* r31 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r32 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b; b -= 5;
   rmxab++;

  a = rmxa;
  *rmxab  = *a++ * *b; b += 3; /* r33 */
  *rmxab += *a++ * *b; b += 3;
  *rmxab += *a   * *b;
}


void RotateRotMx(int *RotMx, const int *RMx, const int *InvRMx)
{
  int  BufMx[9];


  RotMxMultiply(BufMx, RotMx, InvRMx);
  RotMxMultiply(RotMx, RMx,   BufMx);
}


void SeitzMxMultiply(T_RTMx *smxab, const T_RTMx *smxa, const T_RTMx *smxb)
{
  const int  *ar, *a, *b, *bt;
  int        *ab;

  /* no loops to be as fast as posslible */

  ar = smxa->a;
  a  = smxa->a;
  b  = smxb->a;
  ab = smxab->a;

  *ab  = *a++ * *b; b += 3; /* r11 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r12 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r13 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = smxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r21 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r22 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r23 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = smxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r31 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r32 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r33 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b++; bt = b;
   ab++;

  ar = smxa->a;
  *ab  = *ar++ * *b++; /* t1 */
  *ab += *ar++ * *b++;
  *ab += *ar++ * *b; b = bt;
  *ab += *a++;
  *ab %= STBF; if (*ab < 0) *ab += STBF;
   ab++;

  *ab  = *ar++ * *b++; /* t2 */
  *ab += *ar++ * *b++;
  *ab += *ar++ * *b; b = bt;
  *ab += *a++;
  *ab %= STBF; if (*ab < 0) *ab += STBF;
   ab++;

  *ab  = *ar++ * *b++; /* t3 */
  *ab += *ar++ * *b++;
  *ab += *ar   * *b;
  *ab += *a;
  *ab %= STBF; if (*ab < 0) *ab += STBF;
}


void RTMxMultiply(T_RTMx *rtmxab, const T_RTMx *rtmxa, const T_RTMx *rtmxb,
                  int FacAug, int FacTr)
{
  const int  *ar, *a, *b, *bt;
  int        *ab;

  /* no loops to be as fast as posslible */

  ar = rtmxa->a;
  a  = rtmxa->a;
  b  = rtmxb->a;
  ab = rtmxab->a;

  *ab  = *a++ * *b; b += 3; /* r11 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r12 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r13 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = rtmxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r21 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r22 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r23 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b; b = rtmxb->a;
   ab++;

  ar = a;
  *ab  = *a++ * *b; b += 3; /* r31 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r32 */
  *ab += *a++ * *b; b += 3;
  *ab += *a   * *b; b -= 5;
   ab++;

  a = ar;
  *ab  = *a++ * *b; b += 3; /* r33 */
  *ab += *a++ * *b; b += 3;
  *ab += *a++ * *b++; bt = b;
   ab++;

  if (FacTr > 0)
  {
    ar = rtmxa->a;
    *ab  = *ar++ * *b++; /* t1 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
    *ab %= FacTr; if (*ab < 0) *ab += FacTr;
     ab++;

    *ab  = *ar++ * *b++; /* t2 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
    *ab %= FacTr; if (*ab < 0) *ab += FacTr;
     ab++;

    *ab  = *ar++ * *b++; /* t3 */
    *ab += *ar++ * *b++;
    *ab += *ar   * *b;
    *ab += *a   * FacAug;
    *ab %= FacTr; if (*ab < 0) *ab += FacTr;
  }
  else
  {
    ar = rtmxa->a;
    *ab  = *ar++ * *b++; /* t1 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
     ab++;

    *ab  = *ar++ * *b++; /* t2 */
    *ab += *ar++ * *b++;
    *ab += *ar++ * *b; b = bt;
    *ab += *a++ * FacAug;
     ab++;

    *ab  = *ar++ * *b++; /* t3 */
    *ab += *ar++ * *b++;
    *ab += *ar   * *b;
    *ab += *a   * FacAug;
  }
}


int CBMxMultiply(T_RTMx *ab, const T_RTMx *a, const T_RTMx *b)
{
  T_RTMx  abf[1];

  RTMxMultiply(abf, a, b, CRBF, CRBF * CTBF);
  if (ChangeBaseFactor(abf->a, CRBF, ab->a, 1, 12) != 0)
    return IE(-1);

  return 0;
}


int CBMx2Multiply(T_RTMx ab[2], const T_RTMx a[2], const T_RTMx b[2])
{
  if (CBMxMultiply(&ab[0], &a[0], &b[0]) != 0) return -1;
  if (CBMxMultiply(&ab[1], &b[1], &a[1]) != 0) return -1;
  return 0;
}


int CBMx2Update(T_RTMx ab[2], const T_RTMx a[2])
{
  if (CBMxMultiply(&ab[0],  &a[0], &ab[0]) != 0) return -1;
  if (CBMxMultiply(&ab[1], &ab[1],  &a[1]) != 0) return -1;
  return 0;
}


int traceRotMx(const int *RotMx)
{
  return RotMx[0] + RotMx[4] + RotMx[8];
}


int deterRotMx(const int *RotMx)
{
  int     det;

  det =  RotMx[0] * (RotMx[4] * RotMx[8] - RotMx[5] * RotMx[7]);
  det -= RotMx[1] * (RotMx[3] * RotMx[8] - RotMx[5] * RotMx[6]);
  det += RotMx[2] * (RotMx[3] * RotMx[7] - RotMx[4] * RotMx[6]);

  return det;
}


void iCoFactorMxTp(const int *Mx, int *CFMxTp)
{
  CFMxTp[0] =  Mx[4] * Mx[8] - Mx[5] * Mx[7];
  CFMxTp[1] = -Mx[1] * Mx[8] + Mx[2] * Mx[7];
  CFMxTp[2] =  Mx[1] * Mx[5] - Mx[2] * Mx[4];
  CFMxTp[3] = -Mx[3] * Mx[8] + Mx[5] * Mx[6];
  CFMxTp[4] =  Mx[0] * Mx[8] - Mx[2] * Mx[6];
  CFMxTp[5] = -Mx[0] * Mx[5] + Mx[2] * Mx[3];
  CFMxTp[6] =  Mx[3] * Mx[7] - Mx[4] * Mx[6];
  CFMxTp[7] = -Mx[0] * Mx[7] + Mx[1] * Mx[6];
  CFMxTp[8] =  Mx[0] * Mx[4] - Mx[1] * Mx[3];
}


int InverseRotMx(const int R[9], int InvR[9], const int RBF)
{
  /*
     InvR = Inverse[R]
     DetF = |R| * RBF^3
   */

  int  DetF, i;

      DetF = deterRotMx(R);
  if (DetF == 0) return 0;

  iCoFactorMxTp(R, InvR);

  rangei(9) InvR[i] *= (RBF * RBF);

  rangei(9) {
    if (InvR[i] %  DetF) return 0;
        InvR[i] /= DetF;
  }

  return DetF;
}


int InverseRTMx(const T_RTMx *RTMx, T_RTMx *InvRTMx, int RBF)
{
  /*
     InvRTMx->s.R =  Inverse[RTMx->s.R]
     InvRTMx->s.T = -Inverse[RTMx->s.R] * RTMx->s.T
     DetF = |RTMx| * RBF^3;
   */

  int  DetF, i;

      DetF = InverseRotMx(RTMx->s.R, InvRTMx->s.R, RBF);
  if (DetF == 0) return 0;

  RotMx_t_Vector(InvRTMx->s.T, InvRTMx->s.R, RTMx->s.T, 0);

  for (i = 0; i < 3; i++) {
    if (InvRTMx->s.T[i] %  RBF) return 0;
        InvRTMx->s.T[i] /= RBF;
        InvRTMx->s.T[i] *= -1;
  }

  return DetF;
}


void iMxMultiply(int *ab, const int *a, const int *b,
                 const int ma, const int na, const int nb)
{
  int        i, j, k;
  const int  *ai, *aij, *bk, *bkj;

  ai = a;

  for (i = 0; i < ma; i++)
  {
    bk = b;

    for (k = 0; k < nb; k++)
    {
      aij = ai;
      bkj = bk;

      *ab = 0;

      for (j = 0; j < na; j++)
      {
        *ab += (*aij) * (*bkj);

        aij++;
        bkj += nb;
      }

      ab++;
      bk++;
    }

    ai += na;
  }
}


int *IdentityMat(int *M, int m)
{
  int  i;

  rangei(m * m) M[i] = 0;
  rangei(m)     M[i * (m + 1)] = 1;

  return M;
}


int *TransposedMat(int *M, int mr, int mc)
{
  int  ir, ic, i;
  int  *Mt;


  if (mc <= 0 || mr <= 0) return NULL;

  nxs_malloc(Mt, mc * mr);
  if (Mt == NULL) {
    (void) SetSgNotEnoughCore(0);
    return NULL;
  }
  i = 0; range1(ir, mr) range1(ic, mc) Mt[ic * mr + ir] = M[i++];
  (void) MemCpy(M, Mt, mc * mr);
  free(Mt);

  return M;
}


int iRowEchelonFormT(int *M, int mr, int mc, int *T, int tc)
{
  /* C version of RowEchelonFormT from the CrystGAP package
         (GAP Version 3.4.4).
     B. Eick, F. Ga"hler and W. Nickel
         Computing Maximal Subgroups and Wyckoff Positions of Space Groups
         Acta Cryst. (1997). A53, 467 - 474
   */

  int  a, i, j, k, ic, Cleared;


#define  M2D(i, j) M[i * mc + j]
#define  T2D(i, j) T[i * tc + j]

  for (i = j = 0; i < mr && j < mc;)
  {
    k = i; while (k < mr && M2D(k, j) == 0) k++;

    if (k == mr)
      j++;
    else
    {
      if (i != k) {
               (void) IntSwap(&M2D(i, 0), &M2D(k, 0), mc);
        if (T) (void) IntSwap(&T2D(i, 0), &T2D(k, 0), tc);
      }

      range2(k, k + 1, mr) {
        a = abs(M2D(k, j));
        if (a != 0 && a < abs(M2D(i, j))) {
                 (void) IntSwap(&M2D(i, 0), &M2D(k, 0), mc);
          if (T) (void) IntSwap(&T2D(i, 0), &T2D(k, 0), tc);
        }
      }

      if (M2D(i, j) < 0) {
               range1(ic, mc) M2D(i, ic) *= -1;
        if (T) range1(ic, tc) T2D(i, ic) *= -1;
      }

      Cleared = 1;
      range2(k, i + 1, mr) {
        a = M2D(k, j) / M2D(i, j);
        if (a != 0) {
                 range1(ic, mc) M2D(k, ic) -= a * M2D(i, ic);
          if (T) range1(ic, tc) T2D(k, ic) -= a * T2D(i, ic);
        }
        if (M2D(k, j) != 0) Cleared = 0;
      }
      if (Cleared) { i++; j++; }
    }
  }

#undef  M2D
#undef  T2D

  return i;
}


static int IsDiagonalMat(const int *M, int mr, int mc)
{
  int  ir, ic;

  if (mr != mc) return 0;
  range1(ir, mr) range1(ic, mc) if (ir != ic && M[ir * mc + ic]) return 0;

  return 1;
}


int iREBacksubst(const int *M, const int *V,
                 const int nr, const int nc,
                 int *Sol, int *FlagIndep)
{
  int  ir, ic, icp, nv;
  int  d, m, f, jc;

#define M2D(i, j) M[i * nc + j]

  if (FlagIndep)
    for (ic = 0; ic < nc; ic++) FlagIndep[ic] = 1;

  d = 1;

  for (ir = nr - 1; ir >= 0; ir--)
  {
    range1(ic, nc) if (M2D(ir, ic)) goto Set_Sol_ic;

    if (V && V[ir] != 0) return 0;
    continue;

    Set_Sol_ic:

    if (FlagIndep) FlagIndep[ic] = 0;

    if (Sol) {
                    icp = ic + 1;
          nv = nc - icp;
      if (nv) {
        iMxMultiply(&Sol[ic], &M2D(ir, icp), &Sol[icp], 1, nv, 1);
                     Sol[ic] *= -1;
      }
      else
        Sol[ic] = 0;

      if (V) Sol[ic] += d * V[ir];

                        m = M2D(ir, ic);
      f = iGCD(Sol[ic], m);
      if (m < 0) f *= -1;
      Sol[ic] /= f;
      f = m / f;
      if (f != 1) {
        range1(jc, nc) if (jc != ic) Sol[jc] *= f;
        d *= f;
      }
    }
  }

#undef M2D

  return d;
}


int iRESetIxIndep(const int *REMx, int nr, int nc, int *IxIndep, int mIndep)
{
  int  nIndep, ic;
  int  FlagIndep[6];

  if (nc > sizeof FlagIndep / sizeof (*FlagIndep))
    return IE(-1);

  if (iREBacksubst(REMx, NULL, nr, nc, NULL, FlagIndep) < 1)
    return IE(-1);

  nIndep = 0;

  range1(ic, nc) {
    if (FlagIndep[ic]) {
      if (nIndep == mIndep) return -1;
      IxIndep[nIndep++] = ic;
    }
  }

  return nIndep;
}


int SolveHomRE2(const int REMx[9], int EV[3])
{
  /* REMx must be in row echelon form with Rank 2.
   */

  int  IxIndep[1], i;

  if (iRESetIxIndep(REMx, 2, 3, IxIndep, 1) != 1)
    return IE(-1);

  rangei(3) EV[i] = 0;
  EV[IxIndep[0]] = 1;

  if (iREBacksubst(REMx, NULL, 2, 3, EV, NULL) < 1)
    return IE(-1);

  if (SignHemisphere(EV[0], EV[1], EV[2]) < 0)
    rangei(3) EV[i] *= -1;

  return 0;
}


int SolveHomRE1(const int REMx[3], const int IxIndep[2], int Sol[4][3])
{
  int        iPV, i;
  const int  TrialV[4][2] =
    {{ 1,  0 },
     { 0,  1 },
     { 1,  1 },
     { 1, -1 },
    };

  range1(iPV, 4)
  {
    rangei(3) Sol[iPV][i] = 0;
    rangei(2) Sol[iPV][IxIndep[i]] = TrialV[iPV][i];

    if (iREBacksubst(REMx, NULL, 2, 3, Sol[iPV], NULL) < 1)
      return IE(-1);
  }

  return 0;
}


int SmithNormalForm(int *M, int mr, int mc, int *P, int *Q)
{
  int  rr, rc;


  rr = mr;
  rc = mc;

  if (P) (void) IdentityMat(P, mr);
  if (Q) (void) IdentityMat(Q, mc);

  for (;;)
  {
    rr = iRowEchelonFormT(M, rr, rc, P, mr);
    if (IsDiagonalMat(M, rr, rc)) break;
    (void) TransposedMat(M, rr, rc);

    rc = iRowEchelonFormT(M, rc, rr, Q, mc);
    if (IsDiagonalMat(M, rc, rr)) break;
    (void) TransposedMat(M, rc, rr);
  }

  return rr;
}
