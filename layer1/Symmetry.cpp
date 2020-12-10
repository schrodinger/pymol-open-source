
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
#include"os_predef.h"
#include"os_python.h"
#include"os_std.h"

#include"Err.h"
#include"Feedback.h"
#include"Symmetry.h"
#include"P.h"
#include"PConv.h"
#ifdef SYM_TO_MAT_LIST_IN_C
#endif

PyObject *SymmetryAsPyList(CSymmetry * I)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(2);
    PyList_SetItem(result, 0, CrystalAsPyList(&I->Crystal));
    PyList_SetItem(result, 1, PyString_FromString(I->spaceGroup()));
  }
  return (PConvAutoNone(result));
}

static int SymmetryFromPyList(CSymmetry * I, PyObject * list)
{
  auto G = I->G;
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
	ok = CrystalFromPyList(&I->Crystal, list);    
    } else {
      if(ok){
	CPythonVal *val = CPythonVal_PyList_GetItem(I->G, list, 0);
	ok = CrystalFromPyList(&I->Crystal, val);
	CPythonVal_Free(val);
      }
      if(ok) {
        std::string sg;
        ok = PConvFromPyListItem(G, list, 1, sg);
        I->setSpaceGroup(sg.c_str());
      }
    }
    CPythonVal_Free(secondval);
  }
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  return (ok);
}

CSymmetry *SymmetryNewFromPyList(PyMOLGlobals * G, PyObject * list)
{
  CSymmetry *I = NULL;
  I = new CSymmetry(G);
  if(I) {
    if(!SymmetryFromPyList(I, list)) {
      delete I;
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

/**
 * Lookup the symmetry operations by space group symbol (from Python with
 * pymol.xray) and populate SymMatVLA.
 *
 * Return false if space group unknown.
 */
bool CSymmetry::updateSymMatVLA() const
{
  auto I = const_cast<CSymmetry*>(this);

  constexpr bool quiet = false;

  if (I->SymMatVLA) {
    // don't re-run unless setSpaceGroup() was called
    return true;
  }

  int ok = false;
#ifdef SYM_TO_MAT_LIST_IN_C
#else
#ifndef _PYMOL_NOPY
#ifdef _PYMOL_XRAY
  /* TAKEN OUT BB 2/2012  SpaceGroup can be blank, 
     sg_sym_to_mat_list has a blank entry with no operations
  if(!I->SpaceGroup[0]) {
    ErrMessage(G, "Symmetry", "Missing space group symbol");
    } else */
  if(P_xray) {
    int blocked = PAutoBlock(G);
    ov_size a, l;
    PyObject *mats;
    mats = PYOBJECT_CALLMETHOD(P_xray, "sg_sym_to_mat_list", "s", spaceGroup());
    if(mats && (mats != Py_None)) {
      l = PyList_Size(mats);
      I->SymMatVLA = pymol::vla<float>(16*l);
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
#endif

  return (ok);
}

void CSymmetry::setSpaceGroup(const char* sg)
{
  strncpy(SpaceGroup, sg, sizeof(SpaceGroup) - 1);

  SymMatVLA.freeP();
}

int CSymmetry::getNSymMat() const {
  if (!updateSymMatVLA())
    return 0;
  if (!SymMatVLA)
    return 0;
  return SymMatVLA.size() / 16;
}

/**
 * Register a a space group with symmetry operations (if not already registered)
 *
 * @param sg Space group symbol, e.g. "P 1"
 * @param sym_op List of symmetry operations, e.g. ["x,y,z", "-x,-y,z"]
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
