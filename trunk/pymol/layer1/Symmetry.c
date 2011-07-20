
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
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  if(I) {
    result = PyList_New(2);
    PyList_SetItem(result, 0, CrystalAsPyList(I->Crystal));
    PyList_SetItem(result, 1, PyString_FromString(I->SpaceGroup));
  }
  return (PConvAutoNone(result));
#endif

}

int SymmetryFromPyList(CSymmetry * I, PyObject * list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  int ok = true;
  ov_size ll;

  if(ok)
    ok = (I != NULL);
  if(ok)
    SymmetryReset(I);
  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if(ok)
    ok = CrystalFromPyList(I->Crystal, PyList_GetItem(list, 0));
  if(ok)
    ok = PConvPyStrToStr(PyList_GetItem(list, 1), I->SpaceGroup, sizeof(WordType));
  if(ok) {
    ok = SymmetryAttemptGeneration(I, true);
  }
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  return (ok);
#endif
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
static void SymmetryDump44f(PyMOLGlobals * G, float *m, char *prefix)
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

int SymmetryAttemptGeneration(CSymmetry * I, int quiet)
{
  int ok = false;
#ifndef _PYMOL_NOPY
#ifdef _PYMOL_XRAY
  PyMOLGlobals *G = I->G;
  PyObject *mats;
  ov_size a, l;
  CrystalUpdate(I->Crystal);
  if(!quiet) {
    if(Feedback(G, FB_Symmetry, FB_Blather)) {
      CrystalDump(I->Crystal);
    }
  }
  if(!I->SpaceGroup[0]) {
    ErrMessage(G, "Symmetry", "Missing space group symbol");
  } else if(P_xray) {
    int blocked = PAutoBlock(G);
    mats = PyObject_CallMethod(P_xray, "sg_sym_to_mat_list", "s", I->SpaceGroup);
    if(mats && (mats != Py_None)) {
      l = PyList_Size(mats);
      VLACheck(I->SymMatVLA, float, 16 * l);
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
      I->NSymMat = l;
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
  if(I->Crystal)
    CrystalFree(I->Crystal);
  VLAFreeP(I->SymMatVLA);
  VLAFreeP(I->SymOpVLA);
  OOFreeP(I);
}

void SymmetryReset(CSymmetry * I)
{
  I->SpaceGroup[0] = 0;
  I->NSymMat = 0;
  I->NSymOp = 0;
}

CSymmetry *SymmetryNew(PyMOLGlobals * G)
{
  OOAlloc(G, CSymmetry);
  I->G = G;
  I->Crystal = CrystalNew(G);
  I->SpaceGroup[0] = 0;
  I->NSymMat = 0;
  I->SymMatVLA = VLAlloc(float, 16);
  I->NSymOp = 0;
  I->SymOpVLA = VLAlloc(WordType, 1);
  return (I);
}

CSymmetry *SymmetryCopy(CSymmetry * other)
{
  OOAlloc(other->G, CSymmetry);
  if(!other) {
    OOFreeP(I);
    return NULL;
  }
  UtilCopyMem(I, other, sizeof(CSymmetry));
  I->Crystal = CrystalCopy(I->Crystal);
  I->SymMatVLA = VLACopy(I->SymMatVLA, float);
  I->SymOpVLA = VLACopy(I->SymOpVLA, WordType);
  return (I);
}

void SymmetryUpdate(CSymmetry * I)
{
}

void SymmetryDump(CSymmetry * I)
{
}
