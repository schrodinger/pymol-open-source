/* $Id$ */

/* The source code contained in this file is            */
/* Copyright (C) 1994-2000 by Ralf W. Grosse-Kunstleve. */
/* Please see the LICENSE file for more information.    */

#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif

#include <ctype.h>

#include "ExtensionClass.h"

#undef SG_GLOBAL
#include "sglite.h"
#include "sgconst.h"


staticforward PyExtensionClass SgOpsType;
staticforward PyExtensionClass EqMIxType;

#define is_SgOpsObject(v)  ((v)->ob_type == (void *) &SgOpsType)
#define is_EqMIxObject(v)  ((v)->ob_type == (void *) &EqMIxType)


static PyObject *ErrorObject;

#define pReturnPyError(message) \
  { PyErr_SetString(ErrorObject, message); return NULL; }
#define iReturnPyError(message, i) \
  { PyErr_SetString(ErrorObject, message); return (i); }

#define pReturnPySgError() \
  { PyErr_SetString(ErrorObject, SgError); ClrSgError(); return NULL; }
#define iReturnPySgError(i) \
  { PyErr_SetString(ErrorObject, SgError); ClrSgError(); return (i); }


static PyObject *IntArray_as_PyList(const int *a, int n)
{
  int       i;
  PyObject  *l, *v;

        l = PyList_New(n);
  if (! l) return NULL;
  for(i = 0; i < n; i++) {
         v = PyInt_FromLong((long) a[i]);
    if(! v) {
      Py_DECREF(l);
      return NULL;
    }
    PyList_SET_ITEM(l, i, v);
  }

  return l;
}


static PyObject *IntArray_as_PyTuple(const int *a, int n)
{
  int       i;
  PyObject  *l, *v;

        l = PyTuple_New(n);
  if (! l) return NULL;
  for(i = 0; i < n; i++) {
         v = PyInt_FromLong((long) a[i]);
    if(! v) {
      Py_DECREF(l);
      return NULL;
    }
    PyTuple_SET_ITEM(l, i, v);
  }

  return l;
}


typedef struct {
  int  *a;
  int  m;
  int  n;
}
T_IntArray;

static int PySequence_as_IntArray(PyObject *Seq, T_IntArray *a)
{
  int       n, i;
  PyObject  *v;

  if (! PySequence_Check(Seq)) iReturnPyError("integer sequence expected", 0);
  n = PySequence_Length(Seq);
  if (n > a->m) iReturnPyError("sequence too long", 0);
  if (a->n && a->n != n) iReturnPyError("sequence too short", 0);
  rangei(n) {
    v = PySequence_GetItem(Seq, i);
    if (! v) return 0;
    if (! PyNumber_Check(v)) {
      Py_DECREF(v);
      iReturnPyError("sequence may only contain numbers", 0);
    }
    a->a[i] = PyInt_AsLong(v);
    Py_DECREF(v);
    if (PyErr_Occurred()) return 0;
  }
  a->n = n;

  return 1;
}


static PyObject *SgOps__init__(PyObject *self,
                               PyObject *args, PyObject *keywds)
{
  char         *HallSymbol;
  static char  *kwlist[] = { "HallSymbol", NULL };

  ResetSgOps((T_SgOps *) self);

  HallSymbol = NULL;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "|s", kwlist, &HallSymbol))
    return NULL;

  if (HallSymbol) {
    if (ParseHallSymbol(HallSymbol, (T_SgOps *) self, PHSymOptPedantic) < 0)
      pReturnPySgError();
  }
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *SgOps__getstate__(PyObject *self, PyObject *args)
{
  int             i, l, iLTr, iSMx, offs;
  const T_SgOps   *SgOps;
  int             list[6 + SgOps_mLTr * 3 + 3 + SgOps_mSMx * 12];


  SgOps = (const T_SgOps *) self;

  l = 6 + SgOps->nLTr * 3 + 3 + SgOps->nSMx * 12;

  list[0] = SgOps->NoExpand;
  list[1] = SgOps->nLSL;
  list[2] = SgOps->nSSL;
  list[3] = SgOps->nLTr;
  list[4] = SgOps->fInv;
  list[5] = SgOps->nSMx;

  offs = 6;
  range1(iLTr, SgOps->nLTr)
    rangei(3) list[offs + iLTr * 3 + i] = SgOps->LTr[iLTr].v[i];

  offs = 6 + SgOps->nLTr * 3;
  rangei(3) list[offs + i] = SgOps->InvT[i];

  offs = 6 + SgOps->nLTr * 3 + 3;
  range1(iSMx, SgOps->nSMx)
    rangei(12) list[offs + iSMx * 12 + i] = SgOps->SMx[iSMx].a[i];

  return IntArray_as_PyTuple(list, l);
}


static PyObject * SgOps__setstate__(PyObject *self, PyObject *args)
{
  PyObject  *state, *parent;
  PyObject  *v;
  int       i, l, iLTr, iSMx, offs;
  T_SgOps   *SgOps;
  int       list[6 + SgOps_mLTr * 3 + 3 + SgOps_mSMx * 12];


  SgOps = (T_SgOps *) self;
  ResetSgOps(SgOps);

  state = NULL;
  if(! PyArg_ParseTuple(args, "|OO", &state, &parent)) return NULL;

  if(state)
  {
    if (PyDict_Check(state))
      pReturnPyError("Internal Error");

        l = PyObject_Length(state);
    if (l < 0) return NULL;
    if (l < 6 + 3 || l >= sizeof list / sizeof (*list))
      pReturnPyError("Internal Error");

    for (i = 0; i < l; i++) {
            v = PySequence_GetItem(state, i);
      if (! v) return NULL;
      if (! PyInt_Check(v)) {
        Py_DECREF(v);
        pReturnPyError("Internal Error");
      }
      list[i] = (int) PyInt_AsLong(v);
      Py_DECREF(v);
    }

    SgOps->NoExpand = list[0];
    SgOps->nLSL     = list[1];
    SgOps->nSSL     = list[2];
    SgOps->nLTr     = list[3];
    SgOps->fInv     = list[4];
    SgOps->nSMx     = list[5];

    if (l != 6 + SgOps->nLTr * 3 + 3 + SgOps->nSMx * 12)
      pReturnPyError("Internal Error");

    offs = 6;
    range1(iLTr, SgOps->nLTr)
      rangei(3) SgOps->LTr[iLTr].v[i] = list[offs + iLTr * 3 + i];

    offs = 6 + SgOps->nLTr * 3;
    rangei(3) SgOps->InvT[i] = list[offs + i];

    offs = 6 + SgOps->nLTr * 3 + 3;
    range1(iSMx, SgOps->nSMx)
      rangei(12) SgOps->SMx[iSMx].a[i] = list[offs + iSMx * 12 + i];
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static void DEL_SgOpsObject(T_SgOps *self)
{
  PyObject_Del(self);
}


static int PRI_SgOpsObject(T_SgOps *self, FILE *fp, int flags)
{
  if (DumpSgOps(self, fp) != 0) iReturnPySgError(-1);
  return 0;
}


static int CMP_SgOpsObject(PyObject *a, PyObject *b)
{
  T_SgOps  ta[1], tb[1];

  if (! (is_SgOpsObject(a) && is_SgOpsObject(b))) {
    PyErr_SetString(PyExc_TypeError, "can only compare two SgOps objects");
    return 1;
  }

  SgOpsCpy(ta, (const T_SgOps *) a);
  SgOpsCpy(tb, (const T_SgOps *) b);
  if (TidySgOps(ta) != 0) iReturnPySgError(-1);
  if (TidySgOps(tb) != 0) iReturnPySgError( 1);

  return SgOpsCmp(ta, tb);
}


static char w_ParseHallSymbol__doc__[] = "parse Hall symbol";

static PyObject *w_ParseHallSymbol(PyObject *self,
                                   PyObject *args, PyObject *keywds)
{
  char         *HallSymbol;
  static char  *kwlist[] = { "HallSymbol", NULL };

  HallSymbol = NULL;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "s", kwlist, &HallSymbol))
    return NULL;

  if (ParseHallSymbol(HallSymbol, (T_SgOps *) self, PHSymOptPedantic) < 0)
    pReturnPySgError();

  Py_INCREF(Py_None);
  return Py_None;
}


static char w_ExpSgSMx__doc__[] =
  "expand space group by adding a symmetry operation";

static PyObject *w_ExpSgSMx(PyObject *self, PyObject *args, PyObject *keywds)
{
  T_IntArray   a[1];
  T_RTMx       SMx[1];
  static char  *kwlist[] = { "SMx", NULL };

  a->a = SMx->a;
  a->m = a->n = 12;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "O&", kwlist,
      PySequence_as_IntArray, a)) return NULL;
  if (ExpSgSMx((T_SgOps *) self, SMx) != 0) pReturnPySgError();
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject *get_nLTr(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", ((T_SgOps *) self)->nLTr);
}

static PyObject *get_fInv(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", ((T_SgOps *) self)->fInv);
}

static PyObject *get_nSMx(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", ((T_SgOps *) self)->nSMx);
}


static PyObject *getLISMx(PyObject *self, PyObject *args, PyObject *keywds)
{
  const T_SgOps  *SgOps;
  int            iLTr, iInv, iSMx, modulus;
  T_RTMx         LISMx[1];
  static char    *kwlist[] = { "iLTr", "iInv", "iSMx", "modulus", NULL };

  iLTr = 0;
  iInv = 0;
  iSMx = 0;
  modulus = 0;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "|iiii", kwlist,
      &iLTr, &iInv, &iSMx, &modulus))
    return NULL;

  SgOps = (const T_SgOps *) self;
  if (iLTr < 0 || iLTr >= SgOps->nLTr) pReturnPyError("iLTr out of range");
  if (iInv < 0 || iInv >= SgOps->fInv) pReturnPyError("fInv out of range");
  if (iSMx < 0 || iSMx >= SgOps->nSMx) pReturnPyError("iSMx out of range");
  SetLISMx(SgOps, iLTr, iInv, iSMx, LISMx);

  if      (modulus > 0) ViModPositive(LISMx->s.T, 3, STBF);
  else if (modulus < 0) ViModShort(LISMx->s.T, 3, STBF);

  return IntArray_as_PyList(LISMx->a, 12);
}


static PyObject *isChiral(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", isChiralSpaceGroup((const T_SgOps *) self));
}

static PyObject *isEnantiomorphic(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i",
                       isEnantiomorphicSpaceGroup((const T_SgOps *) self));
}


static char w_getSpaceGroupType__doc__[] = "get space group type";

static PyObject *w_getSpaceGroupType(PyObject *self,
                                     PyObject *args, PyObject *keywds)
{
  int            tidyCBMx, buildHallSymbol, SgNumber;
  T_RTMx         CBMx[2];
  const T_SgOps  *SgOps;
  char           HallSymbol[128];
  static char    *kwlist[] = { "tidyCBMx", "buildHallSymbol", NULL };

  tidyCBMx = 0;
  buildHallSymbol = 0;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "|ii", kwlist,
      &tidyCBMx, &buildHallSymbol)) return NULL;

  SgOps = (const T_SgOps *) self;

      SgNumber = GetSpaceGroupType(SgOps, &CBMx[0], &CBMx[1]);
  if (SgNumber < 0)
    pReturnPySgError();

  if (tidyCBMx) {
    if (TidyCBMx(SgOps, SgNumber, CBMx) != 0) pReturnPySgError();
  }

  if (buildHallSymbol) {
    if (BuildHallSymbol(SgOps, SgNumber, CBMx,
                        HallSymbol, sizeof HallSymbol) != 0)
      pReturnPySgError();
    return Py_BuildValue("{s:i,s:O,s:O,s:s}",
      "SgNumber", SgNumber,
      "CBMx",     IntArray_as_PyList(CBMx[0].a, 12),
      "InvCBMx",  IntArray_as_PyList(CBMx[1].a, 12),
      "Hall",     HallSymbol);
  }

  return Py_BuildValue("{s:i,s:O,s:O}",
    "SgNumber", SgNumber,
    "CBMx",     IntArray_as_PyList(CBMx[0].a, 12),
    "InvCBMx",  IntArray_as_PyList(CBMx[1].a, 12));
}


static PyObject *BuildSymbolDict(const T_HM_as_Hall *HM_as_Hall)
{
  char  Buf[2];

  Buf[0] = HM_as_Hall->Extension;
  Buf[1] = '\0';
  return Py_BuildValue("{s:i,s:s,s:s,s:s,s:s,s:s}",
    "SgNumber",  HM_as_Hall->SgNumber,
    "Schoenfl",  HM_as_Hall->Schoenfl,
    "Qualif",    HM_as_Hall->Qualif,
    "HM",        HM_as_Hall->HM,
    "Extension", Buf,
    "Hall",      HM_as_Hall->Hall);
}


static char w_MatchTabulatedSettings__doc__[] =
"Match symmetry operations against tabulated settings";

static PyObject *w_MatchTabulatedSettings(PyObject *self, PyObject *args)
{
  int           SgNumber;
  T_HM_as_Hall  HM_as_Hall[1];


  if (! PyArg_ParseTuple(args, "")) return NULL;

      SgNumber = MatchTabulatedSettings(((const T_SgOps *) self), HM_as_Hall);
  if (SgNumber < 0) pReturnPySgError();
  if (SgNumber == 0) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  return BuildSymbolDict(HM_as_Hall);
}


static char w_get_ss__doc__[] =
"get structure-seminvariant vectors and moduli";

static PyObject *w_get_ss(PyObject *self, PyObject *args)
{
  int       n_ssVM, i;
  T_ssVM    ssVM[3];
  PyObject  *vm, *vmi, *v;

  if (! PyArg_ParseTuple(args, "")) return NULL;

      n_ssVM = Set_ss((const T_SgOps *) self, ssVM);
  if (n_ssVM < 0) pReturnPySgError();

        vm = PyList_New(n_ssVM);
  if (! vm) return NULL;

  rangei(n_ssVM) {
          v = PyList_New(2);
    if (! v) goto finally;
    PyList_SET_ITEM(vm, i, v);
    vmi = v;
    v = IntArray_as_PyList(ssVM[i].V, 3);
    if (! v) goto finally;
    PyList_SET_ITEM(vmi, 0, v);
    v = PyInt_FromLong((long) ssVM[i].M);
    if (! v) goto finally;
    PyList_SET_ITEM(vmi, 1, v);
  }

  return Py_BuildValue("{s:i,s:O}",
    "N",  n_ssVM,
    "VM", vm);

  finally:
    Py_XDECREF(vm);
    return NULL;
}


static char w_get_AddlGenEuclNorm__doc__[] =
"get additional generators of Euclidean normalizer";

static PyObject *w_get_AddlGenEuclNorm(PyObject *self,
                                       PyObject *args, PyObject *keywds)
{
  int          K2L, L2N;
  int          SgNumber;
  T_RTMx       CBMx[2];
  int          nAddlG, i;
  T_RTMx        AddlG[3], CB_AddlG[3];
  PyObject     *ag, *v;
  static char  *kwlist[] = { "K2L", "L2N", NULL };

  K2L = 0;
  L2N = 0;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "|ii", kwlist, &K2L, &L2N))
    return NULL;

      SgNumber = GetSpaceGroupType((const T_SgOps *) self, &CBMx[0], &CBMx[1]);
  if (SgNumber < 1)
    pReturnPySgError();

  if (TidyCBMx((const T_SgOps *) self, SgNumber, CBMx) != 0)
    pReturnPySgError();

      nAddlG = GetRefSetNormAddlG(SgNumber, 0, K2L, L2N, AddlG);
  if (nAddlG < 0)
    pReturnPySgError();

  rangei(nAddlG) {
    if (CB_SMx(&CB_AddlG[i], &CBMx[1], &AddlG[i], &CBMx[0]) != 0)
      pReturnPySgError();
  }

        ag = PyList_New(nAddlG);
  if (! ag) return NULL;

  rangei(nAddlG) {
    v = IntArray_as_PyList(CB_AddlG[i].a, 12);
    if (! v) {
      Py_XDECREF(ag);
      return NULL;
    }
    PyList_SET_ITEM(ag, i, v);
  }

  return Py_BuildValue("{s:i,s:O}",
    "N",  nAddlG,
    "SMx", ag);
}


static int ParseTuple_Int_3_1(PyObject *args, PyObject *keywds, int H[],
                              char **kwlist_3, char **kwlist_1)
{
  T_IntArray  a[1];

  if (! PyArg_ParseTupleAndKeywords(args, keywds, "iii", kwlist_3,
      &H[0], &H[1], &H[2])) {
    PyErr_Clear();
    a->a = H;
    a->m = a->n = 3;
    if (! PyArg_ParseTupleAndKeywords(args, keywds, "O&", kwlist_1,
        PySequence_as_IntArray, a))
      return -1;
  }

  return 0;
}

static int ParseTuple_hkl(PyObject *args, PyObject *keywds, int H[3])
{
  static char  *kwlist_3[] = { "h", "k", "l", NULL };
  static char  *kwlist_1[] = { "hkl", NULL };
  return ParseTuple_Int_3_1(args, keywds, H, kwlist_3, kwlist_1);
}


static char w_isSysAbsMIx__doc__[] =
  "check condition for systematically absent reflections";

static PyObject *w_isSysAbsMIx(PyObject *self,
                               PyObject *args, PyObject *keywds)
{
  int  H[3];

  if (ParseTuple_hkl(args, keywds, H) != 0) return NULL;
  return Py_BuildValue("i", IsSysAbsMIx((const T_SgOps *) self, H, NULL));
}


static char w_isCentricMIx__doc__[] =
  "check if reflection is centric";

static PyObject *w_isCentricMIx(PyObject *self,
                                PyObject *args, PyObject *keywds)
{
  int  H[3];

  if (ParseTuple_hkl(args, keywds, H) != 0) return NULL;
  return Py_BuildValue("i", IsCentricMIx((const T_SgOps *) self, H));
}


static char w_get_PhaseRestriction__doc__[] =
  "get phase restriction for given reflection (as multiple of STBF/pi)";

static PyObject *w_get_PhaseRestriction(PyObject *self,
                                        PyObject *args, PyObject *keywds)
{
  int  H[3];

  if (ParseTuple_hkl(args, keywds, H) != 0) return NULL;
  return Py_BuildValue("i", GetPhaseRestriction((const T_SgOps *) self, H));
}


static char w_get_EpsilonMIx__doc__[] = "get epsilon for given hkl";

static PyObject *w_get_EpsilonMIx(PyObject *self,
                                  PyObject *args, PyObject *keywds)
{
  int  H[3], e;

  if (ParseTuple_hkl(args, keywds, H) != 0) return NULL;

      e = EpsilonMIx((T_SgOps *) self, H);
  if (e < 1)
    pReturnPySgError();

  return Py_BuildValue("i", e);
}


static int ParseTuple_Int_4_2(PyObject *args, PyObject *keywds,
                              int *F, int H[],
                              char **kwlist_4, char **kwlist_2)
{
  T_IntArray  a[1];

  if (! PyArg_ParseTupleAndKeywords(args, keywds, "iiii", kwlist_4,
      F, &H[0], &H[1], &H[2])) {
    PyErr_Clear();
    a->a = H;
    a->m = a->n = 3;
    if (! PyArg_ParseTupleAndKeywords(args, keywds, "iO&", kwlist_2,
        F, PySequence_as_IntArray, a))
      return -1;
  }

  return 0;
}

static int ParseTuple_Fhkl(PyObject *args, PyObject *keywds, int *F, int H[3])
{
  static char  *kwlist_4[] = { "FriedelSymmetry", "h", "k", "l", NULL };
  static char  *kwlist_2[] = { "FriedelSymmetry", "hkl", NULL };
  return ParseTuple_Int_4_2(args, keywds, F, H, kwlist_4, kwlist_2);
}


static char w_get_MultMIx__doc__[] = "get multiplicity for given hkl";

static PyObject *w_get_MultMIx(PyObject *self,
                               PyObject *args, PyObject *keywds)
{
  int  FriedelSym, H[3], M;

  if (ParseTuple_Fhkl(args, keywds, &FriedelSym, H) != 0) return NULL;

      M = MultMIx((T_SgOps *) self, FriedelSym, H);
  if (M < 1)
    pReturnPySgError();

  return Py_BuildValue("i", M);
}


static char getCutParameters__doc__[] = "get ispace cut parameters";

static PyObject *getCutParameters(PyObject *self,
                                  PyObject *args, PyObject *keywds)
{
  int          FriedelSym, CutP[3];
  static char  *kwlist[] = { "FriedelSymmetry", NULL };

  if (! PyArg_ParseTupleAndKeywords(args, keywds, "i", kwlist,
      &FriedelSym)) return NULL;

  if (GetCutParamMIx((const T_SgOps *) self, FriedelSym, CutP) != 0)
    pReturnPySgError();

  return Py_BuildValue("(iii)", CutP[0], CutP[1], CutP[2]);
}


static char get_CBMx_to_primitive__doc__[] =
"get a change-of-basis matrix that transforms the given space group "
"representation to a primitive setting";

static PyObject *get_CBMx_to_primitive(PyObject *self, PyObject *args)
{
  T_RTMx  Z2PCBMx[2];

  if (! PyArg_ParseTuple(args, "")) return NULL;

  if (GetZ2PCBMx((const T_SgOps *) self, Z2PCBMx) != 0) pReturnPySgError();

  return Py_BuildValue("{s:O,s:O}",
    "CBMx",     IntArray_as_PyList(Z2PCBMx[0].a, 12),
    "InvCBMx",  IntArray_as_PyList(Z2PCBMx[1].a, 12));
}


static char w_SgOps_change_basis__doc__[] =
"apply change-of-basis matrix to given space group, return new instance";

static PyObject *w_SgOps_change_basis(PyObject *self,
                                      PyObject *args, PyObject *keywds)
{
  int            GotCBMx[2], i;
  T_IntArray     a[2];
  T_RTMx         CBMx[2];
  T_SgOps        *BC_SgOps;
  const T_SgOps  *SgOps;
  static char    *kwlist[] = { "CBMx", "InvCBMx", NULL };

  rangei(2) {
    InitRTMx(&CBMx[i], CRBF);
    a[i].a = CBMx[i].a;
    a[i].m = a[i].n = 12;
  }
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "|O&O&", kwlist,
    PySequence_as_IntArray, &a[0],
    PySequence_as_IntArray, &a[1])) return NULL;

  rangei(2) GotCBMx[i] = MemCmp(&CBMx[i], CBMx_1_000, 1);

  if      (GotCBMx[0] && ! GotCBMx[1]) {
    if (InverseRTMx(&CBMx[0], &CBMx[1], CRBF) == 0)
      pReturnPyError("CBMx is not invertible");
  }
  else if (GotCBMx[1] && ! GotCBMx[0]) {
    if (InverseRTMx(&CBMx[1], &CBMx[0], CRBF) == 0)
      pReturnPyError("InvCBMx is not invertible");
  }

      BC_SgOps = PyObject_NEW(T_SgOps, (PyTypeObject *) &SgOpsType);
  if (BC_SgOps == NULL) return NULL;

  SgOps = (const T_SgOps *) self;
  ResetSgOps(BC_SgOps);
  if (CB_SgOps(SgOps, &CBMx[0], &CBMx[1], BC_SgOps) != 0) {
    DEL_SgOpsObject(BC_SgOps);
    pReturnPySgError();
  }

  return (PyObject *) BC_SgOps;
}


typedef struct {
  double  *a;
  int  m;
  int  n;
}
T_DoubleArray;

static int PySequence_as_DoubleArray(PyObject *Seq, T_DoubleArray *a)
{
  int       n, i;
  PyObject  *v;

  if (! PySequence_Check(Seq)) iReturnPyError("float sequence expected", 0);
  n = PySequence_Length(Seq);
  if (n > a->m) iReturnPyError("sequence too long", 0);
  if (a->n && a->n != n) iReturnPyError("sequence too short", 0);
  rangei(n) {
    v = PySequence_GetItem(Seq, i);
    if (! v) return 0;
    if (! PyNumber_Check(v)) {
      Py_DECREF(v);
      iReturnPyError("sequence may only contain numbers", 0);
    }
    a->a[i] = PyFloat_AsDouble(v);
    Py_DECREF(v);
    if (PyErr_Occurred()) return 0;
  }
  a->n = n;

  return 1;
}

static char w_check_MetricalMatrix__doc__[] =
"check if metrical matrix is compatible with symmetry operations";

static PyObject *w_check_MetricalMatrix(PyObject *self,
                                        PyObject *args, PyObject *keywds)
{
  double         G[9], tolerance;
  T_DoubleArray  a[1];
  static char    *kwlist[] = { "MetricalMatrix", "tolerance", NULL };

  a->a = G;
  a->m = a->n = 9;
  tolerance = -1.;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "O&|d", kwlist,
    PySequence_as_DoubleArray, a, &tolerance)) return NULL;

  if (CheckMetricalMatrix((const T_SgOps *) self, G, tolerance) != 0) {
    PyErr_SetString(PyExc_ValueError, SgError);
    ClrSgError();
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}


static T_EqMIx *NEW_EqMIxObject(void)
{
  T_EqMIx  *NewObj;

  NewObj = PyObject_NEW(T_EqMIx, (PyTypeObject *) &EqMIxType);
  if (NewObj == NULL) return NULL;

  return NewObj;
}


static PyObject *EqMIx__init__(PyObject *self, PyObject *args)
{
  int  i;

  if (! PyArg_ParseTuple(args, "")) return NULL;

  ((T_EqMIx *) self)->fInv = 0;
  ((T_EqMIx *) self)->N = 0;
  for (i = 0; i < 3; i++) ((T_EqMIx *) self)->H[0][i] = 0;

  Py_INCREF(Py_None);
  return Py_None;
}


static void DEL_EqMIxObject(T_EqMIx *self)
{
  PyObject_Del(self);
}


static int PRI_EqMIxObject(T_EqMIx *self, FILE *fp, int flags)
{
  fprintf(fp, "fInv = %d, N = %d, H = (%d, %d, %d)",
    self->fInv, self->N,
    self->H[0][0], self->H[0][1], self->H[0][2]);

  return 0;
}


static char w_BuildEqMIx__doc__[] = "build equivalent ispace indices";

static PyObject *w_BuildEqMIx(PyObject *self,
                              PyObject *args, PyObject *keywds)
{
  int          FriedelSym, H[3];
  T_EqMIx      *EqMIx;

  if (ParseTuple_Fhkl(args, keywds, &FriedelSym, H) != 0) return NULL;

      EqMIx = NEW_EqMIxObject();
  if (EqMIx == NULL) return NULL;

  if (BuildEqMIx((T_SgOps *) self, FriedelSym, H, EqMIx) < 1) {
    DEL_EqMIxObject(EqMIx);
    pReturnPySgError();
  }

  return (PyObject *) EqMIx;
}


static char EqMIx_get_fInv__doc__[] = "get EqMIx fInv";

static PyObject *EqMIx_get_fInv(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", ((const T_EqMIx *) self)->fInv);
}


static char EqMIx_get_N__doc__[] = "get EqMIx N";

static PyObject *EqMIx_get_N(PyObject *self, PyObject *args)
{
  if (! PyArg_ParseTuple(args, "")) return NULL;
  return Py_BuildValue("i", ((const T_EqMIx *) self)->N);
}


static char EqMIx_get_H__doc__[] = "get EqMIx H";

static PyObject *EqMIx_get_H(PyObject *self, PyObject *args, PyObject *keywds)
{
  int            iInv, iEq, i;
  int            H[3];
  const T_EqMIx  *EqMIx;
  static char    *kwlist[] = { "iInv", "iEq", NULL };


  iInv = 0;
  iEq = 0;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "|ii", kwlist,
      &iInv, &iEq))
    return NULL;

  EqMIx = (const T_EqMIx *) self;
  if (iEq < 0 || iEq >= EqMIx->N) pReturnPyError("iEq out of range");
  MemCpy(H, EqMIx->H[iEq], 3);
  if (iInv) rangei(3) H[i] *= -1;

  return Py_BuildValue("(iii)", H[0], H[1], H[2]);
}


static char get_MasterMIx__doc__[] = "get master index";

static PyObject *get_MasterMIx(PyObject *self,
                               PyObject *args, PyObject *keywds)
{
  int          CutP[3], MasterH[3];
  static char  *kwlist_3[] = { "CutP_h", "CutP_k", "CutP_l", NULL };
  static char  *kwlist_1[] = { "CutP", NULL };

  if (ParseTuple_Int_3_1(args, keywds, CutP, kwlist_3, kwlist_1) != 0)
    return NULL;

  if (GetMasterMIx((const T_EqMIx *) self, CutP, MasterH) != 0)
    pReturnPySgError();

  return Py_BuildValue("(iii)", MasterH[0], MasterH[1], MasterH[2]);
}


static char get_MasterMIx_and_MateID__doc__[] =
"get master index and MateID";

static PyObject *get_MasterMIx_and_MateID(PyObject *self,
                                          PyObject *args, PyObject *keywds)
{
  static char    *kwlist[] = { "CutP", "hkl", "testSysAbs", NULL };
  int            CutP[3], MIx[3], testSysAbs;
  int            MasterMIx[3], MateID;
  const T_SgOps  *SgOps;
  T_IntArray     aC[1], aM[1];

  aC->a = CutP;
  aC->m = aC->n = 3;
  aM->a = MIx;
  aM->m = aM->n = 3;
  testSysAbs = 1;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "O&O&|i", kwlist,
      PySequence_as_IntArray, aC,
      PySequence_as_IntArray, aM,
      &testSysAbs)) return NULL;

  SgOps = (const T_SgOps *) self;

  if (testSysAbs && IsSysAbsMIx(SgOps, MIx, NULL)) {
    PyErr_SetString(PyExc_ValueError, "systematically absent reflection");
    return NULL;
  }

  if (GetMasterMIx_and_MateID(SgOps, CutP, MIx, MasterMIx, &MateID) != 0)
    pReturnPySgError();

  return Py_BuildValue("(iii)i",
    MasterMIx[0], MasterMIx[1], MasterMIx[2],
    MateID);
}


#define METH_VK (METH_VARARGS | METH_KEYWORDS)

static struct PyMethodDef SgOps_methods[] = {
  { "__init__", (PyCFunction) SgOps__init__, METH_VK,
    NULL
  },
  { "__getstate__", SgOps__getstate__, METH_VARARGS, NULL },
  { "__setstate__", SgOps__setstate__, METH_VARARGS, NULL },
  { "ParseHallSymbol", (PyCFunction) w_ParseHallSymbol, METH_VK,
   w_ParseHallSymbol__doc__
  },
  { "ExpSgSMx", (PyCFunction) w_ExpSgSMx, METH_VK,
   w_ExpSgSMx__doc__
  },
  { "get_nLTr", get_nLTr, METH_VARARGS, NULL },
  { "get_fInv", get_fInv, METH_VARARGS, NULL },
  { "get_nSMx", get_nSMx, METH_VARARGS, NULL },
  { "getLISMx", (PyCFunction) getLISMx, METH_VK, NULL },
  { "isChiral", isChiral, METH_VARARGS, NULL },
  { "isEnantiomorphic", isEnantiomorphic, METH_VARARGS, NULL },
  { "getSpaceGroupType", (PyCFunction) w_getSpaceGroupType, METH_VK,
   w_getSpaceGroupType__doc__ },
  { "MatchTabulatedSettings", w_MatchTabulatedSettings, METH_VARARGS,
   w_MatchTabulatedSettings__doc__ },
  { "get_ss", w_get_ss, METH_VARARGS,
   w_get_ss__doc__
  },
  { "get_AddlGenEuclNorm", (PyCFunction) w_get_AddlGenEuclNorm, METH_VK,
   w_get_AddlGenEuclNorm__doc__
  },
  { "isSysAbsMIx", (PyCFunction) w_isSysAbsMIx, METH_VK,
   w_isSysAbsMIx__doc__
  },
  { "isCentricMIx", (PyCFunction) w_isCentricMIx, METH_VK,
   w_isCentricMIx__doc__
  },
  { "get_PhaseRestriction", (PyCFunction) w_get_PhaseRestriction, METH_VK,
   w_get_PhaseRestriction__doc__
  },
  { "get_EpsilonMIx", (PyCFunction) w_get_EpsilonMIx, METH_VK,
   w_get_EpsilonMIx__doc__
  },
  { "get_MultMIx", (PyCFunction) w_get_MultMIx, METH_VK,
   w_get_MultMIx__doc__
  },
  { "getCutParameters", (PyCFunction) getCutParameters, METH_VK,
     getCutParameters__doc__
  },
  { "get_CBMx_to_primitive", get_CBMx_to_primitive, METH_VARARGS,
     get_CBMx_to_primitive__doc__
  },
  { "change_basis", (PyCFunction) w_SgOps_change_basis, METH_VK,
   w_SgOps_change_basis__doc__
  },
  { "check_MetricalMatrix", (PyCFunction) w_check_MetricalMatrix, METH_VK,
   w_check_MetricalMatrix__doc__
  },
  { "get_MasterMIx_and_MateID", (PyCFunction)get_MasterMIx_and_MateID,METH_VK,
     get_MasterMIx_and_MateID__doc__
  },
  { "BuildEqMIx", (PyCFunction) w_BuildEqMIx, METH_VK,
   w_BuildEqMIx__doc__
  },
  { NULL, NULL }
};

static PyObject *GETATTR_SgOpsObject(T_SgOps *self, char *name)
{
  return Py_FindMethod(SgOps_methods, (PyObject *) self, name);
}


static struct PyMethodDef EqMIx_methods[] = {
  { "__init__", EqMIx__init__, METH_VARARGS, NULL },
  { "get_fInv", EqMIx_get_fInv, METH_VARARGS, EqMIx_get_fInv__doc__ },
  { "get_N",    EqMIx_get_N,    METH_VARARGS, EqMIx_get_N__doc__ },
  { "get_H", (PyCFunction) EqMIx_get_H, METH_VK,
    EqMIx_get_H__doc__
  },
  { "get_MasterMIx", (PyCFunction) get_MasterMIx, METH_VK,
     get_MasterMIx__doc__
  },
  { NULL, NULL }
};

static PyObject *GETATTR_EqMIxObject(T_EqMIx *self, char *name)
{
  return Py_FindMethod(EqMIx_methods, (PyObject *) self, name);
}


static char w_SgSymbolLookup__doc__[] =
"look up space group number, HM symbol or Schoenflies symbol";

static PyObject *w_SgSymbolLookup(PyObject *self,
                                  PyObject *args, PyObject *keywds)
{
  char          *Symbol, *Convention;
  int           TableID, status;
  T_HM_as_Hall  HM_as_Hall[1];
  static char   *kwlist[] = { "Symbol", "Convention", NULL };

  Convention = "";
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "s|s", kwlist,
      &Symbol, &Convention)) return NULL;

  for (; *Convention; Convention++) if (! isspace(*Convention)) break;
  TableID = *Convention;
  if (*Convention) {
    for (Convention++; *Convention; Convention++) {
      if (! isspace(*Convention)) {
        pReturnPyError("TableID not recognized");
      }
    }
  }

      status = SgSymbolLookup(TableID, Symbol, HM_as_Hall);
  if (status < 0) pReturnPySgError();
  if (status == 0) {
    if (HM_as_Hall->Hall == NULL) {
      PyErr_SetString(PyExc_ValueError, "space group symbol not recognized");
      return NULL;
    }
    else {
      return Py_BuildValue("{s:s}", "Hall", HM_as_Hall->Hall);
    }
  }

  return BuildSymbolDict(HM_as_Hall);
}


static char w_ParseStrXYZ__doc__[]="parse xyz and return tuple of 12 integer";

static PyObject *w_ParseStrXYZ(PyObject *self,
                               PyObject *args, PyObject *keywds)
{
  char          *StrXYZ;
  int           RBF, TBF;
  T_RTMx        RTMx[1];
  static char   *kwlist[] = { "xyz", "RBF", "TBF", NULL };

  if (! PyArg_ParseTupleAndKeywords(args, keywds, "sii", kwlist,
      &StrXYZ, &RBF, &TBF)) return NULL;
  if (RBF < 1) pReturnPyError("rotation base factor < 1");
  if (TBF < 1) pReturnPyError("translation base factor < 1");
  if (ParseStrXYZ(StrXYZ, 0, RTMx, RBF, TBF) < 0)
    pReturnPyError("xyz string not recognized");

  return IntArray_as_PyList(RTMx->a, 12);
}


static PyObject *w_RTMx2XYZ(PyObject *self, PyObject *args, PyObject *keywds)
{
  T_IntArray   a[1];
  T_RTMx       Mx[1];
  int          RBF, TBF, Decimal, TrFirst, LowerCase;
  char         *Separator;
  char         xyz[256];
  static char  *kwlist[] = { "RTMx", "RBF", "TBF",
                             "Decimal", "TrFirst", "LowerCase",
                             "Separator", NULL };

  a->a = Mx->a;
  a->m = a->n = 12;
  Decimal = 0;
  TrFirst = 0;
  LowerCase = 1;
  Separator = NULL;
  if (! PyArg_ParseTupleAndKeywords(args, keywds, "O&ii|iiis", kwlist,
    PySequence_as_IntArray, a,
    &RBF, &TBF,
    &Decimal, &TrFirst, &LowerCase,
    &Separator)) return NULL;

  if (! RTMx2XYZ(Mx, RBF, TBF, Decimal, TrFirst, LowerCase, Separator,
                 xyz, sizeof xyz))
    pReturnPySgError();

  return Py_BuildValue("s", xyz);
}


static PyObject *w_CmpEqMIx(PyObject *self, PyObject *args)
{
  int   H1[3], H2[3];

  if (! PyArg_ParseTuple(args, "iiiiii",
    &H1[0], &H1[1], &H1[2],
    &H2[0], &H2[1], &H2[2])) {
    PyErr_Clear();
    if (! PyArg_ParseTuple(args, "(iii)(iii)",
      &H1[0], &H1[1], &H1[2],
      &H2[0], &H2[1], &H2[2])) return NULL;
  }

  return Py_BuildValue("i", CmpEqMIx(H1, H2));
}


static char RunTests__doc__[] = "run tests";

static PyObject *RunTests(PyObject *self, PyObject *args)
{
  char  *HallSymbol;
  char  *Mode;
  int   Range;

  Mode = "";
  Range = 1;
  if (! PyArg_ParseTuple(args, "s|si", &HallSymbol, &Mode, &Range))
    return NULL;

  if (RunSgLiteTests(HallSymbol, Mode, Range) < 0) pReturnPySgError();

  Py_INCREF(Py_None);
  return Py_None;
}


/*****************************************************************************
 * TYPE DESCRIPTORS
 *****************************************************************************/

static PyExtensionClass SgOpsType = {
        PyObject_HEAD_INIT(NULL)
        0,                                /*ob_size*/
        "SgOps",                          /*tp_name*/
        sizeof(T_SgOps),                  /*tp_basicsize*/
        0,                                /*tp_itemsize*/
        /* methods */
        (destructor)DEL_SgOpsObject,      /*tp_dealloc*/
        (printfunc)PRI_SgOpsObject,       /*tp_print*/
        (getattrfunc)GETATTR_SgOpsObject, /*tp_getattr*/
        (setattrfunc)0,                   /*tp_setattr*/
        (cmpfunc)CMP_SgOpsObject,         /*tp_compare*/
        (reprfunc)0,                      /*tp_repr*/
        0,                                /*tp_as_number*/
        0,                                /*tp_as_sequence*/
        0,                                /*tp_as_mapping*/
        (hashfunc)0,                      /*tp_hash*/
        (ternaryfunc)0,                   /*tp_call*/
        (reprfunc)0,                      /*tp_str*/

        /* Space for future expansion */
        0L,0L,0L,0L,
        NULL, /* Documentation string */
        METHOD_CHAIN(SgOps_methods)
};


static PyExtensionClass EqMIxType = {
        PyObject_HEAD_INIT(NULL)
        0,                                 /*ob_size*/
        "EqMIx",                           /*tp_name*/
        sizeof(T_EqMIx),                   /*tp_basicsize*/
        0,                                 /*tp_itemsize*/
        /* methods */
        (destructor)DEL_EqMIxObject,       /*tp_dealloc*/
        (printfunc)PRI_EqMIxObject,        /*tp_print*/
        (getattrfunc)GETATTR_EqMIxObject,  /*tp_getattr*/
        (setattrfunc)0,                    /*tp_setattr*/
        (cmpfunc)0,                        /*tp_compare*/
        (reprfunc)0,                       /*tp_repr*/
        0,                                 /*tp_as_number*/
        0,                                 /*tp_as_sequence*/
        0,                                 /*tp_as_mapping*/
        (hashfunc)0,                       /*tp_hash*/
        (ternaryfunc)0,                    /*tp_call*/
        (reprfunc)0,                       /*tp_str*/

        /* Space for future expansion */
        0L,0L,0L,0L,
        NULL, /* Documentation string */
        METHOD_CHAIN(EqMIx_methods)
};


/*****************************************************************************
 * MODULE LOGIC
 *****************************************************************************/

static struct PyMethodDef module_methods[] = {
  { "RunTests", RunTests, METH_VARARGS, RunTests__doc__ },
  { "SgSymbolLookup", (PyCFunction) w_SgSymbolLookup, METH_VK,
   w_SgSymbolLookup__doc__
  },
  { "ParseStrXYZ", (PyCFunction) w_ParseStrXYZ, METH_VK,
   w_ParseStrXYZ__doc__
  },
  { "RTMx2XYZ", (PyCFunction) w_RTMx2XYZ, METH_VK, NULL },
  { "CmpEqMIx", w_CmpEqMIx, METH_VARARGS, NULL },
  { NULL, NULL}
};

#undef METH_VK

static char *module_documentation = "sglite - space group library";

/* C prototype to suppress GCC compiler warning...*/
DL_EXPORT(void)
     initsglite(void);


DL_EXPORT(void)
initsglite(void)
{
  PyObject  *m, *d, *s;

  char *Revision = "$Revision$";

  m = Py_InitModule4("sglite", module_methods, module_documentation,
                     NULL, PYTHON_API_VERSION);

  d = PyModule_GetDict(m);
  s = PyString_FromStringAndSize(Revision + 11, strlen(Revision + 11) - 2);
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);

  PyExtensionClass_Export(d, "SgOps", SgOpsType);
  PyExtensionClass_Export(d, "EqMIx", EqMIxType);
  ErrorObject = PyString_FromString("sglite.error");
  PyDict_SetItemString(d, "error", ErrorObject);

  PyDict_SetItemString(d, "SRBF", Py_BuildValue("i",    1));
  PyDict_SetItemString(d, "STBF", Py_BuildValue("i", STBF));
  PyDict_SetItemString(d, "CRBF", Py_BuildValue("i", CRBF));
  PyDict_SetItemString(d, "CTBF", Py_BuildValue("i", CTBF));

  if (PyErr_Occurred())
    Py_FatalError("can't initialize module sglite");
}
