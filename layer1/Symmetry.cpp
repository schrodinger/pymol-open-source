
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#include"os_python.h"

#include"os_predef.h"
#include"os_python.h"
#include"os_std.h"

#include"MemoryDebug.h"
#include"Err.h"
#include"Base.h"
#include"OOMac.h"
#include"Symmetry.h"
#include"Setting.h"
#include"Ortho.h"
#include"Matrix.h"
#include"P.h"
#include"PConv.h"
#include"Util.h"
#include"PConv.h"

PyObject *SymmetryAsPyList(CSymmetry * I)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(2);
    PyList_SetItem(result, 0, CrystalAsPyList(I->Crystal));
    PyList_SetItem(result, 1, PyString_FromString(I->SpaceGroup));
  }
  return (PConvAutoNone(result));
}

static int SymmetryFromPyList(CSymmetry * I, PyObject * list)
{
  int ok = true;
  ov_size ll;
  PyObject *secondval;

  if(ok)
    ok = (I != NULL);
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if (ok && ll>=2){
    secondval = PyList_GetItem(list, 1);
    if (PyList_Check(secondval)){
      /* if only the crystal, read it */
      if(ok)
	ok = CrystalFromPyList(I->Crystal, list);    
    } else {
      if(ok)
	ok = CrystalFromPyList(I->Crystal, PyList_GetItem(list, 0));
      if(ok)
	PConvPyStrToStr(PyList_GetItem(list, 1), I->SpaceGroup, sizeof(WordType));
    }
  }
  if(ok) {
      SymmetryUpdate(I);
  }
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  return (ok);
}

CSymmetry *SymmetryNewFromPyList(PyMOLGlobals * G, PyObject * list)
{
  CSymmetry *I = NULL;
  I = SymmetryNew(G);
  if(I) {
    if(!SymmetryFromPyList(I, list)) {
      SymmetryFree(I);
      I = NULL;
    }
  }
  return (I);
}

#ifndef _PYMOL_NOPY
#ifdef _PYMOL_XRAY
static void SymmetryDump44f(PyMOLGlobals * G, const float *m, const char *prefix)
{
  if(prefix) {
    PRINTF "%s %12.5f %12.5f %12.5f %12.5f\n", prefix, m[0], m[1], m[2], m[3] ENDF(G);
    PRINTF "%s %12.5f %12.5f %12.5f %12.5f\n", prefix, m[4], m[5], m[6], m[7] ENDF(G);
    PRINTF "%s %12.5f %12.5f %12.5f %12.5f\n", prefix, m[8], m[9], m[10], m[11] ENDF(G);
    PRINTF "%s %12.5f %12.5f %12.5f %12.5f\n", prefix, m[12], m[13], m[14], m[15] ENDF(G);
  } else {
    PRINTF "%12.5f %12.5f %12.5f %12.5f\n", m[0], m[1], m[2], m[3] ENDF(G);
    PRINTF "%12.5f %12.5f %12.5f %12.5f\n", m[4], m[5], m[6], m[7] ENDF(G);
    PRINTF "%12.5f %12.5f %12.5f %12.5f\n", m[8], m[9], m[10], m[11] ENDF(G);
    PRINTF "%12.5f %12.5f %12.5f %12.5f\n", m[12], m[13], m[14], m[15] ENDF(G);
  }
}
#endif
#endif

/*
 * Lookup the symmetry operations by space group symbol (from Python with
 * pymol.xray) and populate SymMatVLA.
 *
 * Return false if space group unknown.
 */
int SymmetryAttemptGeneration(CSymmetry * I, int quiet)
{
  if (I->SymMatVLA) {
    // don't re-run unless SymmetryUpdate() was called
    return true;
  }

  int ok = false;
#ifndef _PYMOL_NOPY
#ifdef _PYMOL_XRAY
  PyMOLGlobals *G = I->G;
  CrystalUpdate(I->Crystal);
  if(!quiet) {
    if(Feedback(G, FB_Symmetry, FB_Blather)) {
      CrystalDump(I->Crystal);
    }
  }
  /* TAKEN OUT BB 2/2012  SpaceGroup can be blank, 
     sg_sym_to_mat_list has a blank entry with no operations
  if(!I->SpaceGroup[0]) {
    ErrMessage(G, "Symmetry", "Missing space group symbol");
    } else */
  if(P_xray) {
    int blocked = PAutoBlock(G);
    ov_size a, l;
    PyObject *mats;
    mats = PYOBJECT_CALLMETHOD(P_xray, "sg_sym_to_mat_list", "s", I->SpaceGroup);
    if(mats && (mats != Py_None)) {
      l = PyList_Size(mats);
      I->SymMatVLA = VLAlloc(float, 16 * l);
      if(!quiet) {
        PRINTFB(G, FB_Symmetry, FB_Details)
        " Symmetry: Found %d symmetry operators.\n", (int) l ENDFB(G);
      }
      for(a = 0; a < l; a++) {
        PConv44PyListTo44f(PyList_GetItem(mats, a), I->SymMatVLA + (a * 16));
        if(!quiet) {
          if(Feedback(G, FB_Symmetry, FB_Blather)) {
            SymmetryDump44f(G, I->SymMatVLA + (a * 16), " Symmetry:");
          }
        }
      }
      ok = true;
      Py_DECREF(mats);
    } else {
      ErrMessage(G, "Symmetry", "Unable to get matrices.");
    }
    PAutoUnblock(G, blocked);
  }
#endif
#endif
  return (ok);
}

void SymmetryFree(CSymmetry * I)
{
  if (!I)
    return;

  SymmetryClear(I);
  OOFreeP(I);
}

void SymmetryClear(CSymmetry * I)
{
  if(I->Crystal)
    CrystalFree(I->Crystal);
  VLAFreeP(I->SymMatVLA);
}

CSymmetry *SymmetryNew(PyMOLGlobals * G)
{
  OOCalloc(G, CSymmetry);
  I->G = G;
  I->Crystal = CrystalNew(G);
  return (I);
}

CSymmetry *SymmetryCopy(const CSymmetry * other)
{
  if (!other) {
    return NULL;
  }

  OOAlloc(other->G, CSymmetry);
  ok_assert(1, I);

  UtilCopyMem(I, other, sizeof(CSymmetry));
  I->Crystal = CrystalCopy(I->Crystal);
  I->SymMatVLA = NULL;

  ok_assert(2, I->Crystal);

  return (I);
ok_except2:
  SymmetryFree(I);
ok_except1:
  return NULL;
}

void SymmetryUpdate(CSymmetry * I)
{
  if(I->Crystal)
    CrystalUpdate(I->Crystal);
  VLAFreeP(I->SymMatVLA);
}

void SymmetryDump(CSymmetry * I)
{
}

/*
 * Get the number of symmetry matrices
 */
int CSymmetry::getNSymMat() const {
  if (!SymMatVLA)
    return 0;
  return VLAGetSize(SymMatVLA) / 16;
}

/*
 * Register a a space group with symmetry operations (if not already registered)
 *
 * sg: space group symbol, e.g. "P 1"
 * sym_op: list of symmetry operations, e.g. ["x,y,z", "-x,-y,z"]
 */
void SymmetrySpaceGroupRegister(PyMOLGlobals * G, const char* sg, const std::vector<std::string>& sym_op) {
#if !defined(_PYMOL_NOPY) && defined(_PYMOL_XRAY)
  if (!P_xray)
    return;

  int blocked = PAutoBlock(G);
  PYOBJECT_CALLMETHOD(P_xray,
      "sg_register_if_unknown", "sN", sg, PConvToPyObject(sym_op));
  PAutoUnblock(G, blocked);
#endif
}
