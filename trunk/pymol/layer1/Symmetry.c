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
#include"PConv.h"



PyObject *SymmetryAsPyList(CSymmetry *I)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(2);
    PyList_SetItem(result,0,CrystalAsPyList(I->Crystal));
    PyList_SetItem(result,1,PyString_FromString(I->SpaceGroup));
  }
  return(PConvAutoNone(result));
  
}

int SymmetryFromPyList(CSymmetry *I,PyObject *list)
{
  int ok=true;
  int ll;

  if(ok) ok = (I!=NULL);
  if(ok) SymmetryReset(I);
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  if(ok) ok = CrystalFromPyList(I->Crystal,PyList_GetItem(list,0));
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,1),I->SpaceGroup,sizeof(WordType));
   if(ok) {
    ok = SymmetryAttemptGeneration(I,true,true);
  }
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
 return(ok);
}

CSymmetry *SymmetryNewFromPyList(PyObject *list)
{
  CSymmetry *I=NULL;
  I=SymmetryNew();
  if(I) {
    if(!SymmetryFromPyList(I,list)) {
      SymmetryFree(I);
      I=NULL;
    }
  }
  return(I);
}

int SymmetryAttemptGeneration(CSymmetry *I,int blocked,int quiet)
{
  int ok = false;

#ifdef _PYMOL_XRAY
  PyObject *mats;
  int a,l;
  CrystalUpdate(I->Crystal);
  if(!quiet) {
    if(Feedback(FB_Symmetry,FB_Details)) {
      CrystalDump(I->Crystal);
    }
  }
  if(!I->SpaceGroup[0]) {
    ErrMessage("Symmetry","Missing space group symbol");
  } else {
    if(!blocked) 
      PBlock();
    mats = PyObject_CallMethod(P_xray,"sg_sym_to_mat_list","s",I->SpaceGroup);
    if(mats&&(mats!=Py_None)) {
      l = PyList_Size(mats);
      VLACheck(I->SymMatVLA,float,16*l);
      for(a=0;a<l;a++) {
        PConv44PyListTo44f(PyList_GetItem(mats,a),I->SymMatVLA+(a*16));
        if(!quiet) {
          if(Feedback(FB_Symmetry,FB_Details)) {
            MatrixDump44f(I->SymMatVLA+(a*16)," Symmetry:");
          }
        }
      }
      I->NSymMat = l;
      ok = true;
      Py_DECREF(mats);
    } else {
      ErrMessage("Symmetry","Unable to get matrices from sglite");
    }
    if(!blocked) 
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

void SymmetryReset(CSymmetry *I)
{
  I->SpaceGroup[0]=0;
  I->NSymMat=0;
  I->NSymOp=0;
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


