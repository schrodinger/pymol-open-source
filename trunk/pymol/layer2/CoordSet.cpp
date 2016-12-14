
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
#include"os_numpy.h"

#include"os_std.h"

#include <algorithm>

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Scene.h"
#include"CoordSet.h"
#include"Color.h"
#include"PConv.h"
#include"P.h"
#include"Matrix.h"
#include"Sphere.h"
#include"Util.h"
#include"Feedback.h"
#include"RepWireBond.h"
#include"RepCylBond.h"
#include"RepDot.h"
#include"RepMesh.h"
#include"RepSphere.h"
#include"RepRibbon.h"
#include"RepCartoon.h"
#include"RepSurface.h"
#include"RepLabel.h"
#include"RepNonbonded.h"
#include"RepNonbondedSphere.h"
#include"RepEllipsoid.h"

#include"PyMOLGlobals.h"
#include"PyMOLObject.h"
#ifdef _PYMOL_IP_EXTRAS
#endif
#include "Executive.h"
#include "Lex.h"


/*========================================================================*/
/*
 * Get coordinate index for given atom index
 */
int CoordSet::atmToIdx(int atm) const {
  if (Obj->DiscreteFlag) {
    if (this == Obj->DiscreteCSet[atm])
      return Obj->DiscreteAtmToIdx[atm];
    return -1;
  }
  return AtmToIdx[atm];
}

/*========================================================================*/
static char sATOM[] = "ATOM  ";
static char sHETATM[] = "HETATM";

int CoordSetValidateRefPos(CoordSet * I)
{
  if(I->RefPos) {
    VLACheck(I->RefPos, RefPosType, I->NIndex);
    return true;
  } else {
    int ok = true && (I->RefPos = VLACalloc(RefPosType, I->NIndex));
    if(ok) {
      int a;
      for(a = 0; a < I->NIndex; a++) {
        float *src = I->Coord + 3 * a;
        copy3f(src, I->RefPos[a].coord);
        I->RefPos[a].specified = true;
      }
    }
    return ok;
  }
}


/*========================================================================*/
int BondCompare(BondType * a, BondType * b)
{
  int ai0 = a->index[0];
  int bi0 = b->index[0];
  if(ai0 == bi0) {
    int ai1 = a->index[1];
    int bi1 = b->index[1];
    if(ai1 == bi1) {
      return 0;
    } else if(ai1 > bi1) {
      return 1;
    } else {
      return -1;
    }
  } else if(ai0 > bi0) {
    return 1;
  } else {
    return -1;
  }
}


/*========================================================================*/
int BondInOrder(BondType * a, int b1, int b2)
{
  return (BondCompare(a + b1, a + b2) <= 0);
}

int CoordSetFromPyList(PyMOLGlobals * G, PyObject * list, CoordSet ** cs)
{
  CoordSet *I = NULL;
  int ok = true;
  int ll = 0;

  if(*cs) {
    (*cs)->fFree();
    *cs = NULL;
  }

  if(list == Py_None) {         /* allow None for CSet */
    *cs = NULL;
  } else {

    if(ok)
      I = CoordSetNew(G);
    if(ok)
      ok = (I != NULL);
    if(ok)
      ok = (list != NULL);
    if(ok)
      ok = PyList_Check(list);
    if(ok)
      ll = PyList_Size(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */
    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->NIndex);
    if(ok)
      ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->NAtIndex);
    if(ok)
      ok = PConvPyListToFloatVLA(PyList_GetItem(list, 2), &I->Coord);
    if(ok)
      ok = PConvPyListToIntVLA(PyList_GetItem(list, 3), &I->IdxToAtm);
    if(ok && (ll > 5))
      ok = PConvPyStrToStr(PyList_GetItem(list, 5), I->Name, sizeof(WordType));
    if(ok && (ll > 6))
      ok = ObjectStateFromPyList(G, PyList_GetItem(list, 6), &I->State);
    if(ok && (ll > 7))
      I->Setting = SettingNewFromPyList(G, PyList_GetItem(list, 7));
    if(ok && (ll > 8))
      ok = PConvPyListToLabPosVLA(PyList_GetItem(list, 8), &I->LabPos);
    if(!ok) {
      if(I)
        I->fFree();
        *cs = NULL;
    } else {
      *cs = I;
    }
  }
  return (ok);
}

/*
 * Coord set as numpy array
 */
PyObject *CoordSetAsNumPyArray(CoordSet * cs, short copy)
{
#ifndef _PYMOL_NUMPY
  PRINTFB(cs->State.G, FB_CoordSet, FB_Errors)
    "No numpy support\n" ENDFB(cs->State.G);
  return NULL;
#else

  PyObject *result = NULL;
  const int base_size = sizeof(float);
  int typenum = -1;
  npy_intp dims[2] = {0, 3};

  import_array1(NULL);

  switch(base_size) {
    case 4: typenum = NPY_FLOAT32; break;
    case 8: typenum = NPY_FLOAT64; break;
  }

  if(typenum == -1) {
    printf("error: no typenum for float size %d\n", base_size);
    return NULL;
  }

  dims[0] = cs->NIndex;

  if(copy) {
    if((result = PyArray_SimpleNew(2, dims, typenum)))
      memcpy(PyArray_DATA((PyArrayObject *)result), cs->Coord, cs->NIndex * 3 * base_size);
  } else {
    result = PyArray_SimpleNewFromData(2, dims, typenum, cs->Coord);
  }

  return result;
#endif
}

PyObject *CoordSetAsPyList(CoordSet * I)
{
  PyObject *result = NULL;

  if(I) {
    int pse_export_version = SettingGetGlobal_f(I->State.G, cSetting_pse_export_version) * 1000;
    bool dump_binary = SettingGetGlobal_b(I->State.G, cSetting_pse_binary_dump) && (!pse_export_version || pse_export_version >= 1765);
    result = PyList_New(9);
    PyList_SetItem(result, 0, PyInt_FromLong(I->NIndex));
    PyList_SetItem(result, 1, PyInt_FromLong(I->NAtIndex));
    PyList_SetItem(result, 2, PConvFloatArrayToPyList(I->Coord, I->NIndex * 3, dump_binary));
    PyList_SetItem(result, 3, PConvIntArrayToPyList(I->IdxToAtm, I->NIndex, dump_binary));
    if(I->AtmToIdx
        && pse_export_version < 1770)
      PyList_SetItem(result, 4, PConvIntArrayToPyList(I->AtmToIdx, I->NAtIndex, dump_binary));
    else
      PyList_SetItem(result, 4, PConvAutoNone(NULL));
    PyList_SetItem(result, 5, PyString_FromString(I->Name));
    PyList_SetItem(result, 6, ObjectStateAsPyList(&I->State));
    PyList_SetItem(result, 7, SettingAsPyList(I->Setting));
    PyList_SetItem(result, 8, PConvLabPosVLAToPyList(I->LabPos, I->NIndex));
    /* TODO symmetry, spheroid, periodic box ... */
  }
  return (PConvAutoNone(result));
}

void CoordSetAdjustAtmIdx(CoordSet * I, int *lookup, int nAtom)

/* performs second half of removal */
{
  /* NOTE: only works in a compressive mode, where lookup[a]<=a  */
  int a;
  int a0;

  PRINTFD(I->State.G, FB_CoordSet)
    " CoordSetAdjustAtmIdx-Debug: entered NAtIndex: %d NIndex %d\n I->AtmToIdx %p\n",
    I->NAtIndex, I->NIndex, (void *) I->AtmToIdx ENDFD;

  if (I->AtmToIdx){
    for(a = 0; a < I->NAtIndex; a++) {
      a0 = lookup[a];
      if(a0 >= 0) {
	I->AtmToIdx[a0] = I->AtmToIdx[a];
      }
    }
  }
  I->NAtIndex = nAtom;
  if (I->AtmToIdx){
    VLASize(I->AtmToIdx, int, nAtom);
  }
  for(a = 0; a < I->NIndex; a++) {
    I->IdxToAtm[a] = lookup[I->IdxToAtm[a]];
  }
  PRINTFD(I->State.G, FB_CoordSet)
    " CoordSetAdjustAtmIdx-Debug: leaving... NAtIndex: %d NIndex %d\n",
    I->NAtIndex, I->NIndex ENDFD;

}


/*========================================================================*/
int CoordSetMerge(ObjectMolecule *OM, CoordSet * I, CoordSet * cs)
{                               /* must be non-overlapping */
  int nIndex;
  int a, i0;
  int ok = true;
  /* calculate new size and make room for new data */
  nIndex = I->NIndex + cs->NIndex;
  VLASize(I->IdxToAtm, int, nIndex);
  CHECKOK(ok, I->IdxToAtm);
  if (ok)
    VLACheck(I->Coord, float, nIndex * 3);
  CHECKOK(ok, I->Coord);
  if (ok){
    for(a = 0; a < cs->NIndex; a++) {
      i0 = a + I->NIndex;
      I->IdxToAtm[i0] = cs->IdxToAtm[a];
      if (OM->DiscreteFlag){
	int idx = cs->IdxToAtm[a];
	OM->DiscreteAtmToIdx[idx] = i0;
	OM->DiscreteCSet[idx] = I;
      } else {
	I->AtmToIdx[cs->IdxToAtm[a]] = i0;
      }
      copy3f(cs->Coord + a * 3, I->Coord + i0 * 3);
    }
  }
  if (ok){
    if(cs->LabPos) {
      if(!I->LabPos)
	I->LabPos = VLACalloc(LabPosType, nIndex);
      else
	VLACheck(I->LabPos, LabPosType, nIndex);
      if(I->LabPos) {
	UtilCopyMem(I->LabPos + I->NIndex, cs->LabPos, sizeof(LabPosType) * cs->NIndex);
      }
    } else if(I->LabPos) {
      VLACheck(I->LabPos, LabPosType, nIndex);
    }
  }
  if (ok){
    if(cs->RefPos) {
      if(!I->RefPos)
	I->RefPos = VLACalloc(RefPosType, nIndex);
      else
	VLACheck(I->RefPos, RefPosType, nIndex);
      if(I->RefPos) {
	UtilCopyMem(I->RefPos + I->NIndex, cs->RefPos, sizeof(RefPosType) * cs->NIndex);
      }
    } else if(I->RefPos) {
      VLACheck(I->RefPos, RefPosType, nIndex);
    }
    I->invalidateRep(cRepAll, cRepInvAll);
  }
  I->NIndex = nIndex;

  return ok;
}


/*========================================================================*/
void CoordSetPurge(CoordSet * I)

/* performs first half of removal  */
{
  int offset = 0;
  int a, a1, ao;
  AtomInfoType *ai;
  ObjectMolecule *obj;
  float *c0, *c1;
  LabPosType *l0, *l1;
  RefPosType *r0, *r1;
  obj = I->Obj;

  PRINTFD(I->State.G, FB_CoordSet)
    " CoordSetPurge-Debug: entering..." ENDFD;

  c0 = c1 = I->Coord;
  r0 = r1 = I->RefPos;
  l0 = l1 = I->LabPos;

  /* This loop slides down the atoms that are not deleted (deleteFlag)
     it moves the Coord, RefPos, and LabPos */
  for(a = 0; a < I->NIndex; a++) {
    a1 = I->IdxToAtm[a];
    ai = obj->AtomInfo + a1;
    if(ai->deleteFlag) {
      offset--;
      c0 += 3;
      if(l0)
        l0++;
      if(r0)
        r0++;
    } else if(offset) {
      ao = a + offset;
      *(c1++) = *(c0++);
      *(c1++) = *(c0++);
      *(c1++) = *(c0++);
      if(r1) {
        *(r1++) = *(r0++);
      }
      if(l0) {
        *(l1++) = *(l0++);
      }
      if (I->AtmToIdx)
	I->AtmToIdx[a1] = ao;
      I->IdxToAtm[ao] = a1;     /* no adjustment of these indexes yet... */
      if (I->Obj->DiscreteFlag){
	I->Obj->DiscreteAtmToIdx[a1] = ao;
	I->Obj->DiscreteCSet[a1] = I;
      }
    } else {
      c0 += 3;
      c1 += 3;
      if(r1) {
        r0++;
        r1++;
      }
      if(l0) {
        l0++;
        l1++;
      }
    }
  }
  if(offset) {
    /* If there were deleted atoms, (offset < 0), then
       re-adjust the array sizes */
    I->NIndex += offset;
    VLASize(I->Coord, float, I->NIndex * 3);
    if(I->LabPos) {
      VLASize(I->LabPos, LabPosType, I->NIndex);
    }
    if(I->RefPos) {
      VLASize(I->RefPos, RefPosType, I->NIndex);
    }
    VLASize(I->IdxToAtm, int, I->NIndex);
    PRINTFD(I->State.G, FB_CoordSet)
      " CoordSetPurge-Debug: I->IdxToAtm shrunk to %d\n", I->NIndex ENDFD;
    I->invalidateRep(cRepAll, cRepInvAtoms);      /* this will free Color */
  }
  PRINTFD(I->State.G, FB_CoordSet)
    " CoordSetPurge-Debug: leaving NAtIndex %d NIndex %d...\n",
    I->NAtIndex, I->NIndex ENDFD;
}


/*========================================================================*/
int CoordSetTransformAtomTTTf(CoordSet * I, int at, const float *TTT)
{
  int a1 = I->atmToIdx(at);
  float *v1;

  if(a1 < 0)
    return false;

  v1 = I->Coord + 3 * a1;
  MatrixTransformTTTfN3f(1, v1, TTT, v1);
  return true;
}


/*========================================================================*/
int CoordSetTransformAtomR44f(CoordSet * I, int at, const float *matrix)
{
  int a1 = I->atmToIdx(at);
  float *v1;

  if(a1 < 0)
    return false;

  v1 = I->Coord + 3 * a1;
  MatrixTransformR44fN3f(1, v1, matrix, v1);
  return true;
}


/*========================================================================*/
void CoordSetRecordTxfApplied(CoordSet * I, const float *matrix, int homogenous)
{
  double temp[16];

  if(!homogenous) {
    convertTTTfR44d(matrix, temp);
  } else {
    convert44f44d(matrix, temp);
  }

  ObjectStateLeftCombineMatrixR44d(&I->State, temp);
}


/*========================================================================*/
int CoordSetMoveAtom(CoordSet * I, int at, const float *v, int mode)
{
  int a1 = I->atmToIdx(at);
  float *v1;

  if(a1 < 0)
    return false;

  v1 = I->Coord + 3 * a1;
  if(mode) {
    add3f(v, v1, v1);
  } else {
    copy3f(v, v1);
  }
  return true;
}


/*========================================================================*/
int CoordSetMoveAtomLabel(CoordSet * I, int at, const float *v, int mode)
{
  ObjectMolecule *obj = I->Obj;
  int a1 = I->atmToIdx(at);
  int result = 0;
  LabPosType *lp;

  /* if label is valid, get the label offset
   * and set the new position relative to that */
  if(a1 >= 0) {
    if(!I->LabPos)
      I->LabPos = VLACalloc(LabPosType, I->NIndex);
    if(I->LabPos) {
      result = 1;
      lp = I->LabPos + a1;
      if(!lp->mode) {
        const float *lab_pos =
          SettingGet_3fv(obj->Obj.G, I->Setting, obj->Obj.Setting,
                         cSetting_label_position);
        copy3f(lab_pos, lp->pos);
      }
      lp->mode = 1;
      if(mode) {
        add3f(v, lp->offset, lp->offset);
      } else {
        copy3f(v, lp->offset);
      }
    }
  }

  return (result);
}


/*========================================================================*/
int CoordSetGetAtomVertex(CoordSet * I, int at, float *v)
{
  int a1 = I->atmToIdx(at);

  if(a1 < 0)
    return false;

  copy3f(I->Coord + 3 * a1, v);
  return true;
}


/*========================================================================*/
int CoordSetGetAtomTxfVertex(CoordSet * I, int at, float *v)
{
  ObjectMolecule *obj = I->Obj;
  int a1 = I->atmToIdx(at);

  if(a1 < 0)
    return false;

  copy3f(I->Coord + 3 * a1, v);

  /* apply state transformation */
  if(I->State.Matrix && (SettingGet_i(I->State.G,
          obj->Obj.Setting, I->Setting,
          cSetting_matrix_mode) > 0)) {
    transform44d3f(I->State.Matrix, v, v);
  }

  /* object transformation */
  if(obj->Obj.TTTFlag) {
    transformTTT44f3f(obj->Obj.TTT, v, v);
  }

  return true;
}


/*========================================================================*/
int CoordSetSetAtomVertex(CoordSet * I, int at, const float *v)
{
  int a1 = I->atmToIdx(at);

  if(a1 < 0)
   return false;

  copy3f(v, I->Coord + 3 * a1);
  return true;
}


/*========================================================================*/
void CoordSetRealToFrac(CoordSet * I, const CCrystal * cryst)
{
  int a;
  float *v;
  v = I->Coord;
  for(a = 0; a < I->NIndex; a++) {
    transform33f3f(cryst->RealToFrac, v, v);
    v += 3;
  }
}


/*========================================================================*/
void CoordSetTransform44f(CoordSet * I, const float *mat)
{
  int a;
  float *v;
  v = I->Coord;
  for(a = 0; a < I->NIndex; a++) {
    transform44f3f(mat, v, v);
    v += 3;
  }
}


/*========================================================================*/

void CoordSetTransform33f(CoordSet * I, const float *mat)
{
  int a;
  float *v;
  v = I->Coord;
  for(a = 0; a < I->NIndex; a++) {
    transform33f3f(mat, v, v);
    v += 3;
  }
}


/*========================================================================*/
void CoordSetGetAverage(CoordSet * I, float *v0)
{
  int a;
  float *v;
  double accum[3];
  if(I->NIndex) {
    v = I->Coord;
    accum[0] = *(v++);
    accum[1] = *(v++);
    accum[2] = *(v++);
    for(a = 1; a < I->NIndex; a++) {
      accum[0] += *(v++);
      accum[1] += *(v++);
      accum[2] += *(v++);
    }
    v0[0] = (float) (accum[0] / I->NIndex);
    v0[1] = (float) (accum[1] / I->NIndex);
    v0[2] = (float) (accum[2] / I->NIndex);
  }
}


/*========================================================================*/
void CoordSetFracToReal(CoordSet * I, const CCrystal * cryst)
{
  int a;
  float *v;
  v = I->Coord;
  for(a = 0; a < I->NIndex; a++) {
    transform33f3f(cryst->FracToReal, v, v);
    v += 3;
  }
}

/*
 * Apply the `sca` (SCALEn) transformation to transform to fractional space,
 * and then the crystals `FracToReal` transformation to transform back to
 * cartesian space.
 *
 * Don't do anything if pdb_insure_orthogonal=off.
 *
 * Don't do anything if SCALEn or CRYST1 look bogus. There is a number of
 * structures in the PDB which have meaningless values for those.
 *
 * Without this, creating symmetry mates might produce wrong results.
 */
bool CoordSetInsureOrthogonal(PyMOLGlobals * G,
    CoordSet * cset,            // coord set to modify
    const float * sca,          // 4x4 SCALE
    const CCrystal *cryst,
    bool quiet)
{
  if (!SettingGetGlobal_b(G, cSetting_pdb_insure_orthogonal))
    return false;

  if (!cryst)
    cryst = cset->Symmetry->Crystal;

  const float * r2f = cryst->RealToFrac;

  // are the matrices sufficiently close to be the same?
  if (!sca[3] && !sca[7] && !sca[11] &&
      is_allclosef(3, r2f, 3, sca, 4, R_SMALL4)) {
    return false;
  }

  // is the cell a orthogonal 1x1x1? If so, then it should probably be ignored...
  // is SCALEn the identity matrix?  If so, then it should probably be ignored...
  if (is_identityf(3, r2f, R_SMALL4) ||
      is_identityf(4, sca, R_SMALL4)) {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " ObjectMolReadPDBStr: ignoring SCALEn (identity matrix).\n" ENDFB(G);
    return false;
  }

  // is SCALEn invalid?  If so, then it should definitely be ignored...
  if (determinant33f(sca, 4) < R_SMALL8 ||
      determinant33f(r2f, 3) < R_SMALL8) {
    PRINTFB(G, FB_ObjectMolecule, FB_Blather)
      " ObjectMolReadPDBStr: ignoring SCALEn (invalid matrix).\n" ENDFB(G);
    return false;
  }

  PRINTFB(G, FB_ObjectMolecule, quiet ? FB_Blather : FB_Actions)
    " ObjectMolecule: using SCALEn to compute orthogonal coordinates.\n"
    ENDFB(G);

  CoordSetTransform44f(cset, sca);
  CoordSetFracToReal(cset, cryst);

  return true;
}

/*========================================================================*/
static char RotateU(const double *matrix, float *anisou)
/* Rotates the ANISOU vector
 *
 * matrix: flat 4x4, but only rotation (upper left 3x3) is considered
 * anisou: has 6 elements (of symmetric 3x3) and will be rotated in-place
 */
{
  int i, j, k;
  float Re[3][3];
  double e_val[3], e_vec[3][3];
  double U[3][3] = {
    { anisou[0], anisou[3], anisou[4] },
    { anisou[3], anisou[1], anisou[5] },
    { anisou[4], anisou[5], anisou[2] },
  };

  // e_val, e_vec = linalg.eigh(U)
  if(!xx_matrix_jacobi_solve(*e_vec, e_val, &i, *U, 3))
    return false;

  // Re = dot(matrix[:3,:3], e_vec)
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      Re[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        Re[i][j] += matrix[i * 4 + k] * e_vec[k][j];
    }

  // U = dot(Re * e_val, Re.T)
  for (i = 0; i < 3; i++)
    for (j = 0; j <= i; j++) {
      U[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        U[i][j] += Re[i][k] * e_val[k] * Re[j][k];
    }

  anisou[0] = U[0][0];
  anisou[1] = U[1][1];
  anisou[2] = U[2][2];
  anisou[3] = U[1][0];
  anisou[4] = U[2][0];
  anisou[5] = U[2][1];

  return true;
}

/*========================================================================*/
void CoordSetAtomToPDBStrVLA(PyMOLGlobals * G, char **charVLA, int *c,
                             const AtomInfoType * ai,
                             const float *v, int cnt,
                             const PDBInfoRec * pdb_info,
                             const double *matrix)
/*
 * v: 3x1 vertex in final output space
 * matrix: 4x4 homogenous transformation matrix from model space to output
 *         space (view matrix * state matrix). Used for ANISOU.
 */
{
  char *aType;
  AtomName name;
  ResName resn;
  ov_word chain;

  char formalCharge[4];
  int literal = SettingGetGlobal_b(G, cSetting_pdb_literal_names);
  int reformat = SettingGetGlobal_i(G, cSetting_pdb_reformat_names_mode);
  int ignore_pdb_segi = SettingGetGlobal_b(G, cSetting_ignore_pdb_segi);
  WordType x, y, z;

  formalCharge[0] = 0;
  sprintf(resn, "%3.4s", LexStr(G, ai->resn));
  if(SettingGetGlobal_b(G, cSetting_pdb_truncate_residue_name)) {
    resn[3] = 0;                /* enforce 3-letter residue name in PDB files */
  }
  if(SettingGetGlobal_b(G, cSetting_pdb_formal_charges)) {
    if((ai->formalCharge > 0) && (ai->formalCharge < 10)) {
      sprintf(formalCharge, "%d+", ai->formalCharge);
    } else if((ai->formalCharge < 0) && (ai->formalCharge > -10)) {
      sprintf(formalCharge, "%d-", -ai->formalCharge);
    }
  }

  if(ai->hetatm)
    aType = sHETATM;
  else
    aType = sATOM;

  char inscode = ai->getInscode(true);

  VLACheck(*charVLA, char, (*c) + 1000);

  const char * ai_name = LexStr(G, ai->name);
  auto ai_name_len = strlen(ai_name);
  bool start_column_1 = false;

  UtilNCopy(name, ai_name, 5);

  if(!ai->name) {
    if(!ai->elem[1])
      sprintf(name, " %s", ai->elem);
    else
      sprintf(name, "%s", ai->elem);
  } else if(!literal) {
    if(ai_name_len < 4) {  /* atom name less than length 4 */
      if(!isdigit(name[0])) {     /* doesn't start with a number */
        if((toupper(ai->elem[0]) == toupper(name[0])) && ((!ai->elem[1]) || /* symbol len = 1 */
              (toupper(ai->elem[1]) == toupper(name[1])))) {        /* matched len 2 */
          /* starts with corrent atomic symbol, so */
          if(!ai->elem[1]) { /* symbol len = 1 */
            switch (reformat) {
            case 1:            /* pdb with internal pdb */
            case 3:            /* pdb with internal iupac */
              if((ai->elem[0] == 'H') && ai_name_len > 2) {
                AtomInfoGetPDB3LetHydroName(G, resn, ai_name, name);
                break;
              }
            default:           /* otherwise, start in column 1 */
              start_column_1 = true;
              break;
            }
          }
        } else {                /* name doesn't start with atomic symbol */
          /* then just place it in column 1 as usual */
          start_column_1 = true;
        }
      } else {                  /* name starts with a number */
        switch (reformat) {
        case 2:                /* make Amber compliant */
          if((ai->elem[0] == name[1]) &&
             ((!ai->elem[1]) || (toupper(ai->elem[1]) == toupper(name[2])))) {
            /* rotate the name to place atom symbol in column 0 to comply with Amber PDB format */
            name[3] = name[0];
            name[0] = ' ';
          }
          break;
        }
      }                         /* just stick it in column 0 and hope for the best */
    } else {                    /* if name is length 4 */
      if((ai->elem[0] == name[0]) && ((!ai->elem[1]) ||     /* symbol len = 1 */
                                          (toupper(ai->elem[1]) == toupper(name[1])))) {    /* matched len 2 */
        /* name starts with the atomic symbol */
        if((!ai->elem[1]) && (ai->elem[0])) {   /* but if element is one letter... */
          switch (reformat) {
          case 1:              /* retaining PDB compliance throughout, or */
          case 3:              /* saving as PDB compliant, but use IUPAC within PyMOL */
            if(isdigit(name[3])) {  /* and last character is a number */
              /* rotate the name to place atom symbol in column 1 to comply with PDB format */
              name[0] = ai_name[3];
              name[1] = ai_name[0];
              name[2] = ai_name[1];
              name[3] = ai_name[2];
              name[4] = 0;
            }
            break;
          }
        }
      } else {                  /* name does not start with the symbol... */
        if(reformat == 2) {     /* AMBER compliance mode */
          if(isdigit(name[0])) {
            if((ai->elem[0] == name[1]) &&
               ((!(ai->elem[1])) || (toupper(ai->elem[1]) == toupper(name[2])))) {
              /* rotate the name to place atom symbol in column 0 to comply with Amber PDB format */
              name[0] = ai_name[1];
              name[1] = ai_name[2];
              name[2] = ai_name[3];
              name[3] = ai_name[0];
              name[4] = 0;
            }
          }
        }
      }
    }
  } else {                      /* LITERAL mode: preserve what was in the original PDB as best PyMOL can 
                                   this should enable people to open and save amber pdb files without issues */
    if (ai_name_len < 4 && !(ai->elem[1] && /* elem len = 2 */
          toupper(ai->elem[0]) == toupper(name[0]) &&
          toupper(ai->elem[1]) == toupper(name[1]))) {
      start_column_1 = true;
    }
  }

  if (start_column_1) {
    name[0] = ' ';
    UtilNCopy(name + 1, ai_name, 4);
  }

  if(SettingGetGlobal_b(G, cSetting_pdb_retain_ids)) {
    cnt = ai->id - 1;
  }
  if(cnt > 99998)
    cnt = 99998;

  name[4] = 0;

  if((!pdb_info) || (!pdb_info->is_pqr_file())) { /* relying upon short-circuit */
    short linelen;
    sprintf(x, "%8.3f", v[0]);
    x[8] = 0;
    sprintf(y, "%8.3f", v[1]);
    y[8] = 0;
    sprintf(z, "%8.3f", v[2]);
    z[8] = 0;
    linelen =
      sprintf((*charVLA) + (*c),
              "%6s%5i %-4s%1s%-4s%1.1s%4i%c   %s%s%s%6.2f%6.2f      %-4.4s%2s%2s\n", aType,
              cnt + 1, name, ai->alt, resn, LexStr(G, ai->chain), ai->resv % 10000, inscode, x, y, z, ai->q, ai->b,
              ignore_pdb_segi ? "" :
              LexStr(G, ai->segi), ai->elem, formalCharge);
    if(ai->anisou) {
      // Columns 7 - 27 and 73 - 80 are identical to the corresponding ATOM/HETATM record.
      char *atomline = (*charVLA) + (*c);
      char *anisoline = atomline + linelen;
      float anisou[6];
      memcpy(anisou, ai->anisou, 6 * sizeof(float));

      if(matrix && !RotateU(matrix, anisou)) {
        PRINTFB(G, FB_CoordSet, FB_Errors) "RotateU failed\n" ENDFB(G);
        return;
      }
      strncpy(anisoline + 6, atomline + 6, 22);
      sprintf(anisoline + 28,
          "%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f",
          anisou[0] * 1e4, anisou[1] * 1e4, anisou[2] * 1e4,
          anisou[3] * 1e4, anisou[4] * 1e4, anisou[5] * 1e4);
      strcpy(anisoline + 70, atomline + 70);
      strncpy(anisoline, "ANISOU", 6);
      (*c) += linelen;
    }
    (*c) += linelen;
  } else {
    Chain alt;
    if(pdb_info->is_pqr_file() && pdb_info->pqr_workarounds) {
      inscode = ' ';
      chain = 0;                /* no chain IDs */
      alt[0] = 0;               /* not alt conf identifiers */
    } else {
      alt[0] = ai->alt[0];
      alt[1] = 0;
      chain = ai->chain;
    }
    sprintf(x, "%8.3f", v[0]);
    if(x[0] != 32)
      sprintf(x, " %7.2f", v[0]);
    x[8] = 0;
    sprintf(y, "%8.3f", v[1]);
    y[8] = 0;
    if(y[0] != 32)
      sprintf(y, " %7.2f", v[1]);
    y[8] = 0;
    sprintf(z, "%8.3f", v[2]);
    if(z[0] != 32)
      sprintf(z, " %7.2f", v[2]);
    z[8] = 0;

    (*c) += sprintf((*charVLA) + (*c), "%6s%5i %-4s%1s%-4s%1.1s%4i%c   %s%s%s %11.8f %7.3f\n",
                    aType, cnt + 1, name, alt, resn,
                    LexStr(G, chain), ai->resv, inscode, x, y, z, ai->partialCharge, ai->elec_radius);
  }

}


/*========================================================================*/
PyObject *CoordSetAtomToChemPyAtom(PyMOLGlobals * G, AtomInfoType * ai, const float *v,
                                   const float *ref, int index, const double *matrix)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *atom = PYOBJECT_CALLMETHOD(P_chempy, "Atom", "");
  if(!atom)
    ErrMessage(G, "CoordSetAtomToChemPyAtom", "can't create atom");
  else {
    float tmp_array[6] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };

    if (ai->anisou) {
      memcpy(tmp_array, ai->anisou, 6 * sizeof(float));
      if (matrix)
      RotateU(matrix, tmp_array);
    }

    char resi[8];
    AtomResiFromResv(resi, sizeof(resi), ai);

    PConvFloat3ToPyObjAttr(atom, "coord", v);
    if(ref)
      PConvFloat3ToPyObjAttr(atom, "ref_coord", ref);
    PConvStringToPyObjAttr(atom, "name", LexStr(G, ai->name));
    PConvStringToPyObjAttr(atom, "symbol", ai->elem);
    PConvStringToPyObjAttr(atom, "resn", LexStr(G, ai->resn));
    PConvStringToPyObjAttr(atom, "resi", resi);
    PConvStringToPyObjAttr(atom, "ss", ai->ssType);
    PConvIntToPyObjAttr(atom, "resi_number", ai->resv);
    PConvIntToPyObjAttr(atom, "stereo", ai->stereo);
    PConvStringToPyObjAttr(atom, "chain", LexStr(G, ai->chain));
    if(ai->alt[0])
      PConvStringToPyObjAttr(atom, "alt", ai->alt);
    PConvStringToPyObjAttr(atom, "segi", LexStr(G, ai->segi));
    PConvFloatToPyObjAttr(atom, "q", ai->q);
    PConvFloatToPyObjAttr(atom, "b", ai->b);
    {
      {
        PyObject *tmp_obj = PConvFloatArrayToPyList(tmp_array, 6);
        if(tmp_obj) {
          PyObject_SetAttrString(atom, "u_aniso", tmp_obj);
          Py_XDECREF(tmp_obj);
        }
      }
    }
    PConvFloatToPyObjAttr(atom, "vdw", ai->vdw);
    PConvFloatToPyObjAttr(atom, "elec_radius", ai->elec_radius);
    PConvFloatToPyObjAttr(atom, "partial_charge", ai->partialCharge);
    PConvIntToPyObjAttr(atom, "formal_charge", ai->formalCharge);
    if(ai->customType != -9999)
      PConvIntToPyObjAttr(atom, "numeric_type", ai->customType);
    PConvStringToPyObjAttr(atom, "text_type", LexStr(G, ai->textType));
    PConvStringToPyObjAttr(atom, "custom", LexStr(G, ai->custom));

    PConvIntToPyObjAttr(atom, "hetatm", ai->hetatm);
    PConvIntToPyObjAttr(atom, "flags", ai->flags);
    PConvIntToPyObjAttr(atom, "id", ai->id);    /* not necc. unique */
    PConvIntToPyObjAttr(atom, "index", index + 1);      /* fragile */
  }
  if(PyErr_Occurred())
    PyErr_Print();
  return (atom);
#endif
}


/*========================================================================*/
void CoordSetAtomToTERStrVLA(PyMOLGlobals * G, char **charVLA, int *c, AtomInfoType * ai,
                             int cnt)
{
  int retain_ids = SettingGetGlobal_b(G, cSetting_pdb_retain_ids);
  int ter_id;

  VLACheck(*charVLA, char, (*c) + 1000);

  if(retain_ids) {
    ter_id = ai->id + 1;
  } else {
    ter_id = cnt + 1;
  }

  (*c) += sprintf((*charVLA) + (*c),
                  "TER   %5i      %3.3s %1.1s%4d%c\n", ter_id, 
                  LexStr(G, ai->resn),
                  LexStr(G, ai->chain), ai->resv,
                  ai->getInscode(true));

}


/*========================================================================*/
void CoordSet::invalidateRep(int type, int level)
{
  CoordSet * I = this;
  int a;
  if(level >= cRepInvVisib) {
    if (I->Obj)
      I->Obj->RepVisCacheValid = false;
  }
  /* graphical representations need redrawing */
  if(level == cRepInvVisib) {
    /* cartoon_side_chain_helper */
    if(SettingGet_b(I->State.G, I->Setting, I->Obj->Obj.Setting,
                    cSetting_cartoon_side_chain_helper)) {
      if((type == cRepCyl) || (type == cRepLine) || (type == cRepSphere))
        invalidateRep(cRepCartoon, cRepInvVisib2);
      else if(type == cRepCartoon) {
        invalidateRep(cRepLine, cRepInvVisib2);
        invalidateRep(cRepCyl, cRepInvVisib2);
        invalidateRep(cRepSphere, cRepInvVisib2);
      }
    }
    /* ribbon_side_chain_helper */
    if(SettingGet_b(I->State.G, I->Setting, I->Obj->Obj.Setting,
                    cSetting_ribbon_side_chain_helper)) {
      if((type == cRepCyl) || (type == cRepLine) || (type == cRepSphere))
        invalidateRep(cRepRibbon, cRepInvVisib2);
      else if(type == cRepRibbon) {
        invalidateRep(cRepLine, cRepInvVisib2);
        invalidateRep(cRepCyl, cRepInvVisib2);
        invalidateRep(cRepSphere, cRepInvVisib2);
      }
    }
    /* line_stick helper  */
    if(SettingGet_b(I->State.G, I->Setting, I->Obj->Obj.Setting,
                    cSetting_line_stick_helper)) {
      if(type == cRepCyl)
        invalidateRep(cRepLine, cRepInvVisib2);
      else if(type == cRepLine) {
        invalidateRep(cRepCyl, cRepInvVisib2);
      }
    }
  }

  if(I->Spheroid)
    if(I->NSpheroid != I->NAtIndex * I->SpheroidSphereSize) {
      FreeP(I->Spheroid);
      FreeP(I->SpheroidNormal);
    }

  /* invalidate basd on one representation, 'type' */
  for (RepIterator iter(I->State.G, type); iter.next(); ){
    int eff_level = level;
    a = iter.rep;
    if(level == cRepInvPick) {
      switch (a) {
      case cRepSurface:
      case cRepMesh:
      case cRepDot:
        /* skip the expensive to recompute, non-pickable
           representations */
        break;
      default:               /* default behavior is to blow away the representation */
        eff_level = cRepInvRep;
        break;
      }
    }
    if(eff_level >= cRepInvVisib)     /* make active if visibility has changed */
      I->Active[a] = true;
    if(I->Rep[a]) {
      if(I->Rep[a]->fInvalidate && (eff_level < cRepInvPurge))
        I->Rep[a]->fInvalidate(I->Rep[a], I, eff_level);
      else if(eff_level >= cRepInvExtColor) {
        I->Rep[a]->fFree(I->Rep[a]);
        I->Rep[a] = NULL;
      }
    }
  }

  if(level >= cRepInvCoord) {   /* if coordinates change, then this map becomes invalid */
    MapFree(I->Coord2Idx);
    I->Coord2Idx = NULL;
    /* invalidate distances */
  }
  SceneChanged(I->State.G);
}


/*========================================================================*/

#define RepUpdateMacro(I,rep,new_fn,state) {\
  if(I->Active[rep]&&(!G->Interrupt)) {\
    if(!I->Rep[rep]) {\
      I->Rep[rep]=new_fn(I,state);\
      if(I->Rep[rep]){ \
         I->Rep[rep]->fNew=(struct Rep *(*)(struct CoordSet *,int state))new_fn;\
      } else {  \
	I->Active[rep] = false;			\
      }         \
    } else {\
      if(I->Rep[rep]->fUpdate)\
         I->Rep[rep] = I->Rep[rep]->fUpdate(I->Rep[rep],I,state,rep);\
    }\
  }\
OrthoBusyFast(I->State.G,rep,cRepCnt);\
}

/*========================================================================*/
void CoordSet::update(int state)
{
  CoordSet * I = this;
  int a;
  PyMOLGlobals *G = I->Obj->Obj.G;

  PRINTFB(G, FB_CoordSet, FB_Blather) " CoordSetUpdate-Entered: object %s state %d cset %p\n",
    I->Obj->Obj.Name, state, (void *) I
    ENDFB(G);

  OrthoBusyFast(G, 0, cRepCnt);
  RepUpdateMacro(I, cRepLine, RepWireBondNew, state);
  RepUpdateMacro(I, cRepCyl, RepCylBondNew, state);
  RepUpdateMacro(I, cRepDot, RepDotNew, state);
  RepUpdateMacro(I, cRepMesh, RepMeshNew, state);
  RepUpdateMacro(I, cRepSphere, RepSphereNew, state);
  RepUpdateMacro(I, cRepRibbon, RepRibbonNew, state);
  RepUpdateMacro(I, cRepCartoon, RepCartoonNew, state);
  RepUpdateMacro(I, cRepSurface, RepSurfaceNew, state);
  RepUpdateMacro(I, cRepLabel, RepLabelNew, state);
  RepUpdateMacro(I, cRepNonbonded, RepNonbondedNew, state);
  RepUpdateMacro(I, cRepNonbondedSphere, RepNonbondedSphereNew, state);
  RepUpdateMacro(I, cRepEllipsoid, RepEllipsoidNew, state);

  for(a = 0; a < cRepCnt; a++)
    if(!I->Rep[a])
      I->Active[a] = false;

  SceneInvalidate(G);
  OrthoBusyFast(G, 1, 1);
  if(Feedback(G, FB_CoordSet, FB_Blather)) {
    printf(" CoordSetUpdate-Leaving: object %s state %d cset %p\n",
           I->Obj->Obj.Name, state, (void *) I);
  }
}


/*========================================================================*/
void CoordSetUpdateCoord2IdxMap(CoordSet * I, float cutoff)
{
  if(cutoff < R_SMALL4)
    cutoff = R_SMALL4;
  if(I->NIndex > 10) {
    if(I->Coord2Idx) {
      if((I->Coord2IdxDiv < cutoff) ||
         (((cutoff - I->Coord2IdxReq) / I->Coord2IdxReq) < -0.5F)) {
        MapFree(I->Coord2Idx);
        I->Coord2Idx = NULL;
      }
    }
    if(I->NIndex && (!I->Coord2Idx)) {  /* NOTE: map based on stored coords */
      I->Coord2IdxReq = cutoff;
      I->Coord2IdxDiv = cutoff * 1.25F;
      I->Coord2Idx = MapNew(I->State.G, I->Coord2IdxDiv, I->Coord, I->NIndex, NULL);
      if(I->Coord2IdxDiv < I->Coord2Idx->Div)
        I->Coord2IdxDiv = I->Coord2Idx->Div;
    }
  }
}


/*========================================================================*/
void CoordSet::render(RenderInfo * info)
{
  CoordSet * I = this;
  PyMOLGlobals *G = I->State.G;
  PRINTFD(G, FB_CoordSet)
    " CoordSetRender: entered (%p).\n", (void *) I ENDFD;

  if(!(info->ray || info->pick) &&
     (SettingGet_i(G, I->Setting, I->Obj->Obj.Setting,
                   cSetting_defer_builds_mode) == 5)) {
    if(!info->pass) {
      ObjectUseColor((CObject *) I->Obj);
      if(I->Active[cRepLine])
        RepWireBondRenderImmediate(I, info);
      if(I->Active[cRepNonbonded])
        RepNonbondedRenderImmediate(I, info);
      if(I->Active[cRepSphere])
        RepSphereRenderImmediate(I, info);
      if(I->Active[cRepCyl])
        RepCylBondRenderImmediate(I, info);
      if(I->Active[cRepRibbon])
        RepRibbonRenderImmediate(I, info);
    }
  } else {
    int pass = info->pass;
    CRay *ray = info->ray;
    Picking **pick = info->pick;
    int a, aa;
    ::Rep *r;
    int float_labels = SettingGet_i(G, I->Setting,
                                    I->Obj->Obj.Setting,
                                    cSetting_float_labels);
    int sculpt_vdw_vis_mode = SettingGet_i(G, I->Setting,
					   I->Obj->Obj.Setting,
					   cSetting_sculpt_vdw_vis_mode);
    if((!pass) && sculpt_vdw_vis_mode && 
       I->SculptCGO && (I->Obj->Obj.visRep & cRepCGOBit)) {
      if(ray) {
        int ok = CGORenderRay(I->SculptCGO, ray,
			      ColorGet(G, I->Obj->Obj.Color), I->Setting, I->Obj->Obj.Setting);
	if (!ok){
	  CGOFree(I->SculptCGO);
	  CGOFree(I->SculptShaderCGO);
	  I->SculptShaderCGO = I->SculptCGO = NULL;
	}
      } else if(G->HaveGUI && G->ValidContext) {
        if(!pick) {
	  int use_shader = SettingGetGlobal_b(G, cSetting_use_shaders);
	  if (use_shader){
	    if (!I->SculptShaderCGO){
	      CGO *convertcgo = NULL;
	      convertcgo = CGOCombineBeginEnd(I->SculptCGO, 0);
	      if (convertcgo){
		I->SculptShaderCGO = CGOOptimizeToVBONotIndexed(convertcgo, 0);
		I->SculptShaderCGO->use_shader = I->SculptShaderCGO->enable_shaders = true;
		CGOFree(convertcgo);
	      }
	    }
	  } else if (I->SculptShaderCGO){
	    CGOFree(I->SculptShaderCGO);
	    I->SculptShaderCGO = NULL;
	  }
	  if (I->SculptShaderCGO){
	    CGORenderGL(I->SculptShaderCGO, ColorGet(G, I->Obj->Obj.Color),
			I->Setting, I->Obj->Obj.Setting, info, NULL);
	  } else {
	    CGORenderGL(I->SculptCGO, ColorGet(G, I->Obj->Obj.Color),
			I->Setting, I->Obj->Obj.Setting, info, NULL);
	  }
        }
      }
    }

    for(aa = 0; aa < cRepCnt; aa++) {
      if(aa == cRepSurface) {   /* reorder */
        a = cRepCell;
      } else if(aa == cRepCell) {
        a = cRepSurface;
      } else {
        a = aa;
      }

      if(I->Active[a] && I->Rep[a]) {
        r = I->Rep[a];
        if(!ray) {
          ObjectUseColor((CObject *) I->Obj);
        } else {
          if(I->Obj)
            ray->wobble(
                         SettingGet_i(G, I->Setting,
                                      I->Obj->Obj.Setting,
                                      cSetting_ray_texture),
                         SettingGet_3fv(G, I->Setting,
                                        I->Obj->Obj.Setting,
                                        cSetting_ray_texture_settings));
          else
            ray->wobble(
                         SettingGet_i(G, I->Setting,
                                      NULL, cSetting_ray_texture),
                         SettingGet_3fv(G, I->Setting, NULL,
                                        cSetting_ray_texture_settings));
          ray->color3fv(ColorGet(G, I->Obj->Obj.Color));
        }

        if(r->fRender) {        /* do OpenGL rendering in three passes */
          if(ray || pick) {

            /* here we need to iterate through and apply coordinate set matrices */

            r->fRender(r, info);
          } else {

            /* here we need to iterate through and apply coordinate set matrices */

            switch (a) {
            case cRepLabel:
              if(float_labels && (pass == -1))
                r->fRender(r, info);
              else if(pass == 1)
                r->fRender(r, info);
              break;
            case cRepNonbondedSphere:
            case cRepRibbon:
            case cRepDot:
            case cRepCGO:
            case cRepCallback:
              if(pass == 1)
                r->fRender(r, info);
              break;
            case cRepLine:
            case cRepMesh:
            case cRepDash:
            case cRepNonbonded:
            case cRepCell:
            case cRepExtent:
              if(!pass)
                r->fRender(r, info);
              break;
            case cRepCyl:      /* render sticks differently depending on transparency */
              if(SettingGet_f(G, r->cs->Setting,
                              r->obj->Setting, cSetting_stick_transparency) > 0.0001) {
                if(pass == -1)
                  r->fRender(r, info);
              } else if(pass == 1){
                r->fRender(r, info);
              }
              break;

            case cRepSurface:
              if(info->alpha_cgo) {
                if(pass == 1)
                  r->fRender(r, info);
              } else {
                if(SettingGet_f(G, r->cs->Setting,
                                r->obj->Setting, cSetting_transparency) > 0.0001) {

                  if(pass == -1)
                    r->fRender(r, info);
                } else if(pass == 1)
                  r->fRender(r, info);
              }
              break;
            case cRepSphere:   /* render spheres differently depending on transparency */
              if(SettingGet_f(G, r->cs->Setting,
                              r->obj->Setting, cSetting_sphere_transparency) > 0.0001) {
                if(pass == -1)
                  r->fRender(r, info);
              } else if(pass == 1)
                r->fRender(r, info);
              break;
            case cRepEllipsoid:        /* render spheres differently depending on transparency */
              if(SettingGet_f(G, r->cs->Setting,
                              r->obj->Setting, cSetting_ellipsoid_transparency) > 0.0001) {
                if(pass == -1)
                  r->fRender(r, info);
              } else if(pass == 1)
                r->fRender(r, info);
              break;
            case cRepCartoon:
              if(info->alpha_cgo) {
                if(pass == 1)
                  r->fRender(r, info);
              } else {
                if(SettingGet_f(G, r->cs->Setting,
                                r->obj->Setting,
                                cSetting_cartoon_transparency) > 0.0001) {
                  if(pass == -1)
                    r->fRender(r, info);
                } else if(pass == 1)
                  r->fRender(r, info);
              }
              break;
            }
          }
        }
        /*          if(ray)
           ray->fWobble(ray,0,NULL); */
      }
    }
  }
  PRINTFD(G, FB_CoordSet)
    " CoordSetRender: leaving...\n" ENDFD;
}


/*========================================================================*/
CoordSet *CoordSetNew(PyMOLGlobals * G)
{
  OOCalloc(G, CoordSet);        /* NULL-initializes all fields */
  ObjectStateInit(G, &I->State);
  I->State.G = G;
  I->PeriodicBoxType = cCSet_NoPeriodicity;

  I->SpheroidSphereSize = I->State.G->Sphere->Sphere[1]->nDot;  /* does this make any sense? */

  return (I);
}

/*========================================================================*/
CoordSet *CoordSetCopy(const CoordSet * cs)
{
  if (!cs)
    return NULL;

  PyMOLGlobals * G = const_cast<PyMOLGlobals*>(cs->State.G);

  /* OOAlloc declares and defines, I:
   * I = ... */
  OOCalloc(G, CoordSet);
  /* shallow copy */
  (*I) = (*cs);                 /* NOTE: must deep-copy all pointers in this struct */
  /* deep copy state struct */
  ObjectStateCopy(&I->State, &cs->State);
  /* deep copy & return ptr to new symmetry */
  I->Symmetry = SymmetryCopy(cs->Symmetry);

  if(I->PeriodicBox)
    I->PeriodicBox = CrystalCopy(I->PeriodicBox);

  // copy VLAs
  I->Coord      = VLACopy2(cs->Coord);
  I->LabPos     = VLACopy2(cs->LabPos);
  I->RefPos     = VLACopy2(cs->RefPos);
  I->AtmToIdx   = VLACopy2(cs->AtmToIdx);
  I->IdxToAtm   = VLACopy2(cs->IdxToAtm);

  UtilZeroMem(I->Rep, sizeof(::Rep *) * cRepCnt);

  I->Setting = NULL;
  I->SculptCGO = NULL;
  I->SculptShaderCGO = NULL;
  I->TmpLinkBond = NULL;
  I->TmpBond = NULL;
  I->Spheroid = NULL;
  I->SpheroidNormal = NULL;
  I->Coord2Idx = NULL;
  return (I);
}


/*========================================================================*/
int CoordSet::extendIndices(int nAtom)
{
  CoordSet * I = this;
  int a, b;
  ObjectMolecule *obj = I->Obj;
  int ok = true;
  if(obj->DiscreteFlag) {
    ok = obj->setNDiscrete(nAtom);

    if(I->AtmToIdx) {           /* convert to discrete if necessary */
      VLAFree(I->AtmToIdx);
      I->AtmToIdx = NULL;
      if (ok){
	for(a = 0; a < I->NIndex; a++) {
	  b = I->IdxToAtm[a];
	  obj->DiscreteAtmToIdx[b] = a;
	  obj->DiscreteCSet[b] = I;
	}
      }
    }
  }
  if(ok && I->NAtIndex < nAtom) {
    if(I->AtmToIdx) {
      VLASize(I->AtmToIdx, int, nAtom);
      CHECKOK(ok, I->AtmToIdx);
      if(ok && nAtom) {
        for(a = I->NAtIndex; a < nAtom; a++)
          I->AtmToIdx[a] = -1;
      }
      I->NAtIndex = nAtom;
    } else if(!obj->DiscreteFlag) {
      I->AtmToIdx = VLACalloc(int, nAtom);
      CHECKOK(ok, I->AtmToIdx);
      if (ok){
	for(a = 0; a < nAtom; a++)
	  I->AtmToIdx[a] = -1;
	I->NAtIndex = nAtom;
      }
    }
  }
  return ok;
}


/*========================================================================*/
void CoordSet::appendIndices(int offset)
{
  CoordSet * I = this;
  int a, b;
  ObjectMolecule *obj = I->Obj;

  I->IdxToAtm = VLACalloc(int, I->NIndex);
  if(I->NIndex) {
    ErrChkPtr(I->State.G, I->IdxToAtm);
    for(a = 0; a < I->NIndex; a++)
      I->IdxToAtm[a] = a + offset;
  }
  if(obj->DiscreteFlag) {
    VLACheck(obj->DiscreteAtmToIdx, int, I->NIndex + offset);
    VLACheck(obj->DiscreteCSet, CoordSet *, I->NIndex + offset);
    for(a = 0; a < I->NIndex; a++) {
      b = a + offset;
      obj->DiscreteAtmToIdx[b] = a;
      obj->DiscreteCSet[b] = I;
    }
  } else {
    I->AtmToIdx = VLACalloc(int, I->NIndex + offset);
    if(I->NIndex + offset) {
      ErrChkPtr(I->State.G, I->AtmToIdx);
      for(a = 0; a < offset; a++)
        I->AtmToIdx[a] = -1;
      for(a = 0; a < I->NIndex; a++)
        I->AtmToIdx[a + offset] = a;
    }
  }
  I->NAtIndex = I->NIndex + offset;
}


/*========================================================================*/
void CoordSet::enumIndices()
{
  CoordSet * I = this;
  /* set up for simple case where 1 = 1, etc. */
  int a;
  I->AtmToIdx = VLACalloc(int, I->NIndex);
  I->IdxToAtm = VLACalloc(int, I->NIndex);
  if(I->NIndex) {
    ErrChkPtr(I->State.G, I->AtmToIdx);
    ErrChkPtr(I->State.G, I->IdxToAtm);
    for(a = 0; a < I->NIndex; a++) {
      I->AtmToIdx[a] = a;
      I->IdxToAtm[a] = a;
    }
  }
  I->NAtIndex = I->NIndex;
}


/*========================================================================*/
void CoordSet::fFree()
{
  CoordSet * I = this;
  int a;
  ObjectMolecule *obj;
  if(I) {
    for(a = 0; a < cRepCnt; a++)
      if(I->Rep[a])
        I->Rep[a]->fFree(I->Rep[a]);
    obj = I->Obj;
    if(obj)
      if(obj->DiscreteFlag)     /* remove references to the atoms in discrete objects */
        for(a = 0; a < I->NIndex; a++) {
          obj->DiscreteAtmToIdx[I->IdxToAtm[a]] = -1;
          obj->DiscreteCSet[I->IdxToAtm[a]] = NULL;
        }
    VLAFreeP(I->AtmToIdx);
    VLAFreeP(I->IdxToAtm);
    MapFree(I->Coord2Idx);
    VLAFreeP(I->Coord);
    VLAFreeP(I->TmpBond);
    if(I->Symmetry)
      SymmetryFree(I->Symmetry);
    if(I->PeriodicBox)
      CrystalFree(I->PeriodicBox);
    FreeP(I->Spheroid);
    FreeP(I->SpheroidNormal);
    SettingFreeP(I->Setting);
    ObjectStatePurge(&I->State);
    CGOFree(I->SculptCGO);
    VLAFreeP(I->LabPos);
    VLAFreeP(I->RefPos);
    /* free and make null */
    OOFreeP(I);
  }
}
void LabPosTypeCopy(const LabPosType * src, LabPosType * dst){
  dst->mode = src->mode;
  copy3f(src->pos, dst->pos);
  copy3f(src->offset, dst->offset);
}
void RefPosTypeCopy(const RefPosType * src, RefPosType * dst){
  copy3f(src->coord, dst->coord);
  dst->specified = src->specified;
}
