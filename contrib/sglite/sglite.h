/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifndef SGLITE_H__
#define SGLITE_H__

#ifdef PythonTypes
#include "Python.h"
#else
#define PyObject_HEAD
#endif

#define MemDbg(t, s, n)  printf("MemDbg %ld %ld %ld",\
                         (long)(t), (long)(s), (long)((n) * sizeof (*(t))));\
                         puts("")
#define MemCpy(t, s, n)  memcpy((t), (s), (n) * sizeof (*(t)))
#define MemCmp(t, s, n)  memcmp((t), (s), (n) * sizeof (*(t)))
#define MemSet(s, c, n)  memset((s), (c), (n) * sizeof (*(s)))

#define nxs_malloc(ptr, n) (ptr) = malloc((n) * sizeof (*(ptr)))
#define free_null(ptr) { if (ptr) { free(ptr); (ptr) = NULL; } }

#define ToThe2(i) ((i) * (i))
#define ToThe3(i) ((i) * (i) * (i))

#ifndef rangei
#define rangei(n)          for ((i) =   0; (i) < (n); (i)++)
#define range1(i, n)       for ((i) =   0; (i) < (n); (i)++)
#define range2(i, f, n)    for ((i) = (f); (i) < (n); (i)++)
#define range3(i, f, n, s) for ((i) = (f); (i) < (n); (i) += (s))
#endif

#define CheckPoint printf("@ %s %d\n", __FILE__, __LINE__); fflush(stdout)

#ifndef SG_GLOBAL
extern
const char *SgError;
extern
char        SgErrorBuffer[128];
#else
const char  *SgError = NULL;
char         SgErrorBuffer[128];
#endif

#define  IE(status) \
        SetSg_InternalError((status), __FILE__, __LINE__)

#define pIE(ptr) \
       pSetSg_InternalError((ptr), __FILE__, __LINE__)

#define SetSgNotEnoughCore(status) \
       SetSg_NotEnoughCore((status), __FILE__, __LINE__)

#define pSetSgNotEnoughCore(ptr) \
       pSetSg_NotEnoughCore((ptr), __FILE__, __LINE__)


#define STBF  12 /* Seitz           Matrix Translation Base Factor */

#define CRBF  12 /* Change of Basis Matrix Rotation    Base Factor */
#define CTBF  72 /* Change of Basis Matrix Translation Base Factor */

/* CAUTION: (CTBF / STBF) has to be an INTEGER */
/* CAUTION: CTBF >= 2 * STBF                   */


typedef struct
  {
    int  v[3];
  }
  T_LTr;


typedef union
  {
    struct { int R[9], T[3]; } s;
    int                        a[12];
  }
  T_RTMx;


#define SgOps_mLTr  (108)
#define SgOps_mSMx  (24)

typedef struct
  {
    PyObject_HEAD
    int     NoExpand;
    int     nLSL;
    int     nSSL;
    int     nLTr;
    int     fInv;
    int     nSMx;
    T_LTr    LTr[SgOps_mLTr];
    int     InvT[3];
    T_RTMx   SMx[SgOps_mSMx];
  }
  T_SgOps;

#define SgOpsOrderP(SgOps)                 ((SgOps)->nSMx * (SgOps)->fInv)
#define SgOpsOrderL(SgOps) ((SgOps)->nLTr * (SgOps)->nSMx * (SgOps)->fInv)


typedef struct
  {
    int  Rtype;
    int  EV[3];
    int  SenseOfRotation;
  }
  T_RMxI;


typedef struct
  {
    PyObject_HEAD
    int  fInv;
    int  N;        /* Number of equivalent hkl to follow */
    int  H[24][3]; /* Indices */
    int  TH[24];   /* Phase shift relative to H[0] */
  }
  T_EqMIx;


#define InitRotMx(RotMx, diagonal)\
  {\
    int  private_i_;\
    for (private_i_ = 0; private_i_ <  9; private_i_++)\
        (RotMx)[private_i_] = (private_i_ % 4 ? 0 : diagonal);\
  }

#define InitRTMx(RTMx_, diagonal)\
  {\
    int  private_i_;\
    for (private_i_ = 0; private_i_ < 12; private_i_++)\
      (RTMx_)->a[private_i_] = (private_i_ % 4 ? 0 : diagonal);\
  }


/* sgglobal.c */

void SetSgError(const char *msg);
void ClrSgError(void);
int SetSg_InternalError(int status, const char *file, const int line);
const void *pSetSg_InternalError(const void *ptr,
                                 const char *file, const int line);
int SetSg_NotEnoughCore(int status, const char *file, const int line);
const void *pSetSg_NotEnoughCore(const void *ptr,
                                 const char *file, const int line);


/* sggen.c */

void ResetLLTr(T_LTr *LLTr, int *nLLTr);
int ExpLLTr(int LTBF, int mLLTr,
            T_LTr *LLTr, int *nLLTr, const int *NewLTr);
void ResetSgOps(T_SgOps *SgOps);
int ExpSgLTr(T_SgOps *SgOps, const int *NewLTr);
int ExpSgInv(T_SgOps *SgOps, const int *InvT);
int ExpSgSMx(T_SgOps *SgOps, const T_RTMx *NewSMx);
int ExpSgRMx(T_SgOps *SgOps, const int NewRMx[9]);


/* sgcb.c */

int CB_IT(const int SignI, const int T[3],
          const T_RTMx *CBMx, const T_RTMx *InvCBMx,
          int BC_T[3]);
int CB_RMx(int CRiC[9],
           const int CBMxR[9], const int RMx[9], const int InvCBMxR[9]);
int CB_SMx(T_RTMx *CSiC,
           const T_RTMx *CBMx, const T_RTMx *SMx, const T_RTMx *InvCBMx);
int CB_SgLTr(const T_SgOps *SgOps,
             const T_RTMx *CBMx, const T_RTMx *InvCBMx,
             T_SgOps *BC_SgOps);
int CB_SgOps(const T_SgOps *SgOps,
             const T_RTMx *CBMx, const T_RTMx *InvCBMx,
             T_SgOps *BC_SgOps);


/* sgltr.c */

typedef struct
  {
    int          Sym;
    int          nLTr;
    const T_LTr  *LTr;
  }
  T_ConvCType;

int ExpSgSymCType(T_SgOps *SgOps, int Sym);
int GetSymCType(int nLTr, const T_LTr *LTr);
int GetZ2PCBMx(const T_SgOps *SgOps, T_RTMx CBMx[2]);


/* sgcharmx.c */

int SenseOfRotation(const int *R, int Rtype, const int *EV);
int GetRtype(const int *RotMx);
int SetRotMxInfo(const int *R, T_RMxI *RI);
int OrderOfRtype(int Rtype);
int MakeCumRMx(const int *R, int Rtype, int *CumRMx);
int Set_wI_Tr(const int *R, const int *T, const T_RMxI *RI,
              int wI[3] /* STBF */,
              int Tr[3] /* CTBF */);


/* sgss.c */

typedef struct
  {
    int  V[3];
    int  M;
  }
  T_ssVM;

int Is_ss(const T_ssVM *ssVM, int n_ssVM, int h, int k, int l);
void Set_uvw(const T_ssVM *ssVM, int n_ssVM, int h, int k, int l, int uvw[3]);
int Set_ss(const T_SgOps *SgOps, T_ssVM ssVM[3]);


/* sgstr.c */

int ParseStrXYZ(const char *StrXYZ, int StopChar,
                T_RTMx *RTMx, int FacR, int FacT);
const char *FormatFraction(int nume, int deno, int Decimal,
                           char *Buffer, int SizeBuffer);
const char *RTMx2XYZ(const T_RTMx *RTMx, int RBF, int TBF,
                     int Decimal, int TrFirst, int Low,
                     const char *Separator,
                     char *BufferXYZ, int SizeBufferXYZ);


/* sgfile.c */

int DumpSgOps(const T_SgOps *SgOps, FILE *fp);


/* sgsymbols.c */

typedef struct
  {
    int         SgNumber;
    const char  *Schoenfl;
    const char  *Qualif;
    const char  *HM;
    int         Extension;
    const char  *Hall;
  }
  T_HM_as_Hall;

int SgSymbolLookup(int TableID, const char *Symbol, T_HM_as_Hall *HM_as_Hall);
int MatchTabulatedSettings(const T_SgOps *SgOps, T_HM_as_Hall *HM_as_Hall);


/* sghall.c */

#define PHSymOptPedantic  (1)
#define PHSymOptNoCType   (2)

int ParseHallSymbolCBMx(const char *HSym, T_SgOps *SgOps, int Options,
                        T_RTMx CBMx[2], int *HaveCBMx);
int ParseHallSymbol(const char *HSym, T_SgOps *SgOps, int Options);


/* sghkl.c */

int CmpEqMIx(const int H1[3], const int H2[3]);
int IsSysAbsMIx(const T_SgOps *SgOps, const int H[3], int *TH_Restriction);
int IsCentricMIx(const T_SgOps *SgOps, const int H[3]);
int GetPhaseRestriction(const T_SgOps *SgOps, const int H[3]);
int EpsilonMIx(const T_SgOps *SgOps, const int H[3]);
int MultMIx(const T_SgOps *SgOps, int FriedelSym, const int H[3]);
int BuildEqMIx(const T_SgOps *SgOps, int FriedelSym, const int H[3],
               T_EqMIx *EqMIx);
int GetCutParamMIx(const T_SgOps *SgOps, int FriedelSym, int CutP[3]);
int GetMasterMIx(const T_EqMIx *EqMIx, const int CutP[3], int MasterH[3]);
int GetMasterMIx_and_MateID(const T_SgOps *SgOps,
                            const int CutP[3], const int MIx[3],
                            int MasterMIx[3], int *MateID);


/* sgmath.c */

int iGCD(const int a, const int b);
int FindGCD(const int *S, int nS);
int CancelGCD(int *S, int nS);
int CancelBFGCD(int *S, int nS, int BF);
int iLCM(const int a, const int b);
int FindLCM(const int *S, int nS);
void SimplifyFraction(int nume, int deno, int *o_nume, int *o_deno);
int iScalProd(const int *u, const int *v, const int *PseudoG);
void iCrossProd(int *rxs, const int *r, const int *s, const int *PseudoG);
int AreLinDepV(const int a[3], const int b[3]);
void SMx_t_InvT(const T_RTMx *SMx, const int InvT[3], T_RTMx *ProdSMx);
void RotMx_t_Vector(int *R_t_V, const int *RotMx, const int *Vector, int FacTr);
void RotMxMultiply(int *rmxab, const int *rmxa, const int *rmxb);
void RotateRotMx(int *RotMx, const int *RMx, const int *InvRMx);
void SeitzMxMultiply(T_RTMx *smxab, const T_RTMx *smxa, const T_RTMx *smxb);
void RTMxMultiply(T_RTMx *rtmxab, const T_RTMx *rtmxa, const T_RTMx *rtmxb,
                  int FacAug, int FacTr);
int CBMxMultiply(T_RTMx *ab, const T_RTMx *a, const T_RTMx *b);
int CBMx2Multiply(T_RTMx ab[2], const T_RTMx a[2], const T_RTMx b[2]);
int CBMx2Update(T_RTMx ab[2], const T_RTMx a[2]);
int traceRotMx(const int *RotMx);
int deterRotMx(const int *RotMx);
void iCoFactorMxTp(const int *Mx, int *CFMxTp);
int InverseRotMx(const int R[9], int InvR[9], const int RBF);
int InverseRTMx(const T_RTMx *RTMx, T_RTMx *InvRTMx, int RBF);
void iMxMultiply(int *ab, const int *a, const int *b,
                 const int ma, const int na, const int nb);
int *IdentityMat(int *M, int m);
int *TransposedMat(int *M, int mr, int mc);
int iRowEchelonFormT(int *M, int mr, int mc, int *T, int tc);
int iREBacksubst(const int *M, const int *V,
                 const int nr, const int nc,
                 int *Sol, int *FlagIndep);
int iRESetIxIndep(const int *REMx, int nr, int nc, int *IxIndep, int mIndep);
int SolveHomRE2(const int REMx[9], int EV[3]);
int SolveHomRE1(const int REMx[3], const int IxIndep[2], int Sol[4][3]);
int SmithNormalForm(int *M, int mr, int mc, int *P, int *Q);


/* sgmetric.c */

int CheckMetricalMatrix(const T_SgOps *SgOps, const double *G,
                        double tolerance);


/* sgnorm.c */

int GetRefSetNormAddlG(int SgNumber, int affine, int UseK2L, int UseL2N,
                       T_RTMx *AddlG);
int CheckMonoRefSetAffNormRestrictions(int SgNumber, const int M[9], int BF);


/* sgprop.c */

int isChiralSpaceGroup(const T_SgOps *SgOps);
int isEnantiomorphicSpaceGroup(const T_SgOps *SgOps);


/* sgtidy.c */

int CmpiVect(const int *a, const int *b, int n);
int TidySgOps(T_SgOps *SgOps);


/* sgtype.c */

int GetPG(const T_SgOps *SgOps);
int GetSpaceGroupType(const T_SgOps *SgOps, T_RTMx *CBMx, T_RTMx *InvCBMx);
int TidyCBMx(const T_SgOps *SgOps, int SgNumber, T_RTMx CBMx[2]);
int BuildHallSymbol(const T_SgOps *SgOps, int SgNumber, const T_RTMx CBMx[2],
                    char *HallSymbol, int sizeHallSymbol);


/* sgutil.c */

int *IntSwap(int *a, int *b, int n);
void IntSetZero(int *a, int n);
int IntIsZero(const int *a, int n);
int iModPositive(int ix, int iy);
int iModShort(int ix, int iy);
void ViModShort(int *ix, int n, int iy);
void ViModPositive(int *ix, int n, int iy);
T_RTMx *SetLISMx(const T_SgOps *SgOps, int iLTr, int iInv, int iSMx,
                 T_RTMx *LISMx);
int Discretize(double fVal, int *iVal, int Fac);
int ChangeBaseFactor(const int *Old, int OldBF, int *New, int NewBF, int n);
int SignHemisphere(int h, int k, int l);
void SetRminusI(const int *R, int *RmI, int Inv);
int NextOf_n_from_m(int m, int n, int *ix);
void SgOpsCpy(T_SgOps *t, const T_SgOps *s);
int SgOpsCmp(const T_SgOps *t, const T_SgOps *s);


/* runtests.c */

int RunSgLiteTests(const char *HallSymbol, const char *Mode, int Range);


#endif /* SGLITE_H__ */
