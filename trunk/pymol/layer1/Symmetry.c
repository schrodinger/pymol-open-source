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

int SymmetryAttemptGeneration(CSymmetry *I)
{
  int ok = false;

#ifdef _PYMOL_XRAY
  PyObject *mats;
  int a,l;
  CrystalUpdate(I->Crystal);
  if(Feedback(FB_Symmetry,FB_Details)) {
    CrystalDump(I->Crystal);
  }
  if(!I->SpaceGroup[0]) {
    ErrMessage("Symmetry","Missing space group symbol");
  } else {
    PBlock();
    mats = PyObject_CallMethod(P_xray,"sg_sym_to_mat_list","s",I->SpaceGroup);
    if(mats&&(mats!=Py_None)) {
      l = PyList_Size(mats);
      VLACheck(I->SymMatVLA,float,16*l);
      for(a=0;a<l;a++) {
        PConv44PyListTo44f(PyList_GetItem(mats,a),I->SymMatVLA+(a*16));
        MatrixDump44f(I->SymMatVLA+(a*16)," Symmetry:");
      }
      I->NSymMat = l;
      ok = true;
      Py_DECREF(mats);
    } else {
      ErrMessage("Symmetry","Unable to get matrices from sglite");
    }
    PUnblock();
  }
#endif
  return(ok);
}

void SymmetryFree(CSymmetry *I)
{
  if(I->Crystal) CrystalFree(I->Crystal);
  VLAFreeP(I->SymMatVLA);
  VLAFreeP(I->SymOpVLA);
  OOFreeP(I);
}

CSymmetry *SymmetryNew(void)
{
  OOAlloc(CSymmetry);

  I->Crystal=CrystalNew();
  I->SpaceGroup[0]=0;
  I->NSymMat=0;
  I->SymMatVLA=VLAlloc(float,16);
  I->NSymOp=0;
  I->SymOpVLA=VLAlloc(WordType,1);
  return(I);
}

CSymmetry *SymmetryCopy(CSymmetry *other)
{
  OOAlloc(CSymmetry);
  if(!other) {
    OOFreeP(I);
    return NULL;
  }
  UtilCopyMem(I,other,sizeof(CSymmetry));
  I->Crystal=CrystalCopy(I->Crystal);
  I->SymMatVLA=VLACopy(I->SymMatVLA,float);
  I->SymOpVLA=VLACopy(I->SymOpVLA,WordType);
  return(I);
}

void SymmetryUpdate(CSymmetry *I) 
{
}

void SymmetryDump(CSymmetry *I) 
{
  
}


