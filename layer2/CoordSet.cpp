
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
#include"RepSphereImmediate.h"
#include"RepRibbon.h"
#include"RepCartoon.h"
#include"RepSurface.h"
#include"RepLabel.h"
#include"RepNonbonded.h"
#include"RepNonbondedSphere.h"
#include"RepEllipsoid.h"
#include"Symmetry.h"

#include"PyMOLGlobals.h"
#include"PyMOLObject.h"
#include "Executive.h"
#include "Lex.h"

#ifdef _PYMOL_IP_PROPERTIES
#include "Property.h"
#endif

/**
 * Get atom coordinates, taking symmetry operation into account.
 *
 * @param v_out Buffer for return value, may or may not be used.
 * @param inv If true, apply the inverse of symop
 * @return Coordinate pointer, or NULL if symmetry operation was invalid.
 */
float const* CoordSet::coordPtrSym(
    int idx, pymol::SymOp const& symop, float* v_out, bool inv) const
{
  auto const* v_in = coordPtr(idx);

  // default symop evaluates to false
  if (!symop) {
    return v_in;
  }

  auto* sym = getSymmetry();

  if (!sym || (symop.index && symop.index >= sym->getNSymMat())) {
    return nullptr;
  }

  copy3f(v_in, v_out);

  double const* const statemat = getPremultipliedMatrix();

  if (statemat) {
    transform44d3f(ObjectStateGetInvMatrix(this), v_out, v_out);
  }

  transform33f3f(sym->Crystal.realToFrac(), v_out, v_out);

  if (inv) {
    v_out[0] -= symop.x;
    v_out[1] -= symop.y;
    v_out[2] -= symop.z;
  }

  if (symop.index) {
    auto mat = sym->getSymMat(symop.index);
    if (inv) {
      inverse_transform44f3f(mat, v_out, v_out);
    } else {
      transform44f3f(mat, v_out, v_out);
    }
  }

  if (!inv) {
    v_out[0] += symop.x;
    v_out[1] += symop.y;
    v_out[2] += symop.z;
  }

  transform33f3f(sym->Crystal.fracToReal(), v_out, v_out);

  if (statemat) {
    transform44d3f(statemat, v_out, v_out);
  }

  return v_out;
}

/*========================================================================*/
/**
 * Get coordinate index for given atom index
 */
int CoordSet::atmToIdx(int atm) const {
  if (Obj->DiscreteFlag) {
    if (this == Obj->DiscreteCSet[atm])
      return Obj->DiscreteAtmToIdx[atm];
    return -1;
  }
  assert(atm < AtmToIdx.size());
  return AtmToIdx[atm];
}

void CoordSet::updateNonDiscreteAtmToIdx(unsigned natom)
{
  assert(!Obj || natom == Obj->NAtom);
  AtmToIdx.resize(natom);
  std::fill_n(AtmToIdx.data(), natom, -1);
  for (unsigned idx = 0, idx_end = getNIndex(); idx != idx_end; ++idx) {
    auto const atm = IdxToAtm[idx];
    assert(atm < natom);
    AtmToIdx[atm] = idx;
  }
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
    int ok = true && (I->RefPos = pymol::vla<RefPosType>(I->NIndex));
    if(ok) {
      int a;
      for(a = 0; a < I->NIndex; a++) {
        const float* src = I->coordPtr(a);
        copy3f(src, I->RefPos[a].coord);
        I->RefPos[a].specified = true;
      }
    }
    return ok;
  }
}

namespace
{
/// Return negative if lhs<rhs, zero if lhs==rhs, positive if lhs>rhs.
inline int cmp(unsigned lhs, unsigned rhs)
{
  return lhs < rhs ? -1 : (lhs == rhs ? 0 : 1);
}

inline int cmp(pymol::SymOp const& lhs, pymol::SymOp const& rhs)
{
  return cmp( //
      reinterpret_cast<unsigned const&>(lhs),
      reinterpret_cast<unsigned const&>(rhs));
}
} // namespace

/*========================================================================*/
/**
 * Compare bonds by atom indices and symmetry operations. Swapped indices are
 * not considered as equal. Bond order, ids and settings are not considered.
 */
int BondCompare(BondType const* a, BondType const* b)
{
  int c;
  (c = cmp(a->index[0], b->index[0])) == 0 &&
      (c = cmp(a->index[1], b->index[1])) == 0 &&
      (c = cmp(a->symop_2, b->symop_2));
  return c;
}


/*========================================================================*/
int BondInOrder(BondType const* a, int b1, int b2)
{
  return (BondCompare(a + b1, a + b2) <= 0);
}

int CoordSetFromPyList(PyMOLGlobals * G, PyObject * list, CoordSet ** cs)
{
  CoordSet *I = NULL;
  int ok = true;
  int ll = 0;

  if(*cs) {
    delete *cs;
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
      ok = PConvPyListToFloatVLA(PyList_GetItem(list, 2), &I->Coord);
    if(ok){
      PConvFromPyListItem(G, list, 3, I->IdxToAtm);
    }
    if(ok && (ll > 5))
      ok = CPythonVal_PConvPyStrToStr_From_List(G, list, 5, I->Name, sizeof(WordType));
    if(ok && (ll > 6)){
      CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 6);
      ok = ObjectStateFromPyList(G, val, I);
      CPythonVal_Free(val);
    }
    if(ok && (ll > 7)){
      CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 7);
      I->Setting.reset(SettingNewFromPyList(G, val));
      CPythonVal_Free(val);
    }

#ifdef _PYMOL_IP_PROPERTIES
#endif

    if(ok && (ll > 10)){
      CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 10);
      if (!CPythonVal_IsNone(val))
	I->SculptCGO = CGONewFromPyList(G, val, 0, 1);
      else {
	I->SculptShaderCGO = I->SculptCGO = NULL;
      }
      CPythonVal_Free(val);
    }
    if(ok){
      if(ll > 11){
	// load in atom-state settings
	CPythonVal *val = CPythonVal_PyList_GetItem(G, list, 11);
	if (!CPythonVal_IsNone(val)){
	  int a;
	  I->atom_state_setting_id = pymol::vla<int>(I->NIndex);
	  for (a=0; a<I->NIndex; a++){
	    CPythonVal *val2 = CPythonVal_PyList_GetItem(G, val, a);
	    if (!CPythonVal_IsNone(val2)){
	      CPythonVal_PConvPyIntToInt(val2, &I->atom_state_setting_id[a]);
	      if (I->atom_state_setting_id[a]) {
	        I->atom_state_setting_id[a] = SettingUniqueConvertOldSessionID(G, I->atom_state_setting_id[a]);
	      }
	    }
	    CPythonVal_Free(val2);
	  }
	} else {
	  I->atom_state_setting_id = NULL;
	}
	CPythonVal_Free(val);
      } else {
        // translate legacy LabPos to label_placement_offset
        auto val = CPythonVal_PyList_GetItem(G, list, 8);
        auto res = CPythonVal_PConvPyListToLabPosVec(G, val);
        if (res && !res->empty()) {
          auto const& LabPos = *res;
          for (int a = 0; a < I->NIndex; ++a) {
            if (length3f(LabPos[a].offset) > R_SMALL4) {
              SettingSet(
                  cSetting_label_placement_offset, LabPos[a].offset, I, a);
            }
          }
        }
        CPythonVal_Free(val);
      }
    }

    if (ll > 12) {
      CPythonVal* val = CPythonVal_PyList_GetItem(G, list, 12);
      I->Symmetry.reset(SymmetryNewFromPyList(G, val));
      CPythonVal_Free(val);
    }

    if(!ok) {
      delete I;
      *cs = NULL;
    } else {
      *cs = I;
    }
  }
  return (ok);
}

/**
 * Coord set as numpy array
 */
PyObject *CoordSetAsNumPyArray(CoordSet * cs, short copy)
{
#ifndef _PYMOL_NUMPY
  PRINTFB(cs->G, FB_CoordSet, FB_Errors)
    "No numpy support\n" ENDFB(cs->G);
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
    result = PyArray_SimpleNewFromData(2, dims, typenum, cs->Coord.data());
  }

  return result;
#endif
}

PyObject *CoordSetAsPyList(CoordSet * I)
{
  PyObject *result = NULL;

  if(I) {
    auto G = I->G;
    int pse_export_version = SettingGet<float>(G, cSetting_pse_export_version) * 1000;
    bool dump_binary = SettingGet<bool>(G, cSetting_pse_binary_dump) && (!pse_export_version || pse_export_version >= 1765);
    result = PyList_New(13);
    PyList_SetItem(result, 0, PyInt_FromLong(I->NIndex));
    int const NAtIndex = I->AtmToIdx.size();
    PyList_SetItem(result, 1, PyInt_FromLong(NAtIndex ? NAtIndex : I->Obj->NAtom)); // legacy
    PyList_SetItem(result, 2, PConvFloatArrayToPyList(I->Coord, I->NIndex * 3, dump_binary));
    PyList_SetItem(result, 3, PConvIntArrayToPyList(I->IdxToAtm.data(), I->NIndex, dump_binary));
    if (!I->AtmToIdx.empty()
        && pse_export_version < 1770)
      PyList_SetItem(result, 4, PConvIntArrayToPyList(I->AtmToIdx.data(), NAtIndex, dump_binary));
    else
      PyList_SetItem(result, 4, PConvAutoNone(NULL));
    PyList_SetItem(result, 5, PyString_FromString(I->Name));
    PyList_SetItem(result, 6, ObjectStateAsPyList(I));
    PyList_SetItem(result, 7, SettingAsPyList(I->Setting.get()));
    PyList_SetItem(result, 8, PConvAutoNone(nullptr) /* LabPos */);

    PyList_SetItem(result, 9,
#ifdef _PYMOL_IP_PROPERTIES
#endif
        PConvAutoNone(Py_None));

    if(I->SculptCGO) {
      PyList_SetItem(result, 10, CGOAsPyList(I->SculptCGO));
    } else {
      PyList_SetItem(result, 10, PConvAutoNone(NULL));
    }
    if (I->has_any_atom_state_settings()) {
      int a;
      PyObject *settings_list = NULL;
      settings_list = PyList_New(I->NIndex);
      for (a=0; a<I->NIndex; a++){
        if (I->has_atom_state_settings(a)) {
	  PyList_SetItem(settings_list, a, PyInt_FromLong(I->atom_state_setting_id[a]));
	} else {
	  PyList_SetItem(settings_list, a, PConvAutoNone(NULL));
	}
      }
      PyList_SetItem(result, 11, settings_list);
    } else {
      PyList_SetItem(result, 11, PConvAutoNone(NULL));
    }
    PyList_SetItem(result, 12, SymmetryAsPyList(I->Symmetry.get()));
    /* TODO spheroid, periodic box ... */
  }
  return (PConvAutoNone(result));
}

/**
 * Updates IdxToAtm and adjusts all NIndex sized arrays.
 *
 * @param lookup atm_old to atm_new index mapping.
 * Deleted atoms have `atm_new = -1`.
 *
 * @pre `lookup[a] <= a`
 */
void CoordSetAdjustAtmIdx(CoordSet* I, const int* lookup)
{
  auto G = I->G;
  int offset = 0;

  for (int idx = 0; idx < I->getNIndex(); ++idx) {
    auto const idx_new = idx + offset;
    auto const atm_new = lookup[I->IdxToAtm[idx]];

    assert(I->IdxToAtm[idx] >= atm_new);
    I->IdxToAtm[idx_new] = atm_new;

    if (atm_new == -1) {
      --offset;

      if (I->has_atom_state_settings(idx)) {
        SettingUniqueDetachChain(G, I->atom_state_setting_id[idx]);
        I->atom_state_setting_id[idx] = 0;
      }
    } else if (offset) {
      copy3f(I->coordPtr(idx), I->coordPtr(idx_new));

      if (I->RefPos) {
        I->RefPos[idx_new] = std::move(I->RefPos[idx]);
      }

      if (I->has_atom_state_settings(idx)) {
        I->atom_state_setting_id[idx_new] = std::move(I->atom_state_setting_id[idx]);
        I->atom_state_setting_id[idx] = 0;
      }
    }
  }

  assert(offset <= 0);

  if (offset < 0) {
    I->setNIndex(I->getNIndex() + offset);
    I->invalidateRep(cRepAll, cRepInvAtoms); /* this will free Color */
  }
}

CoordSet* CoordSetCopy(const CoordSet* src)
{
  if(!src) {
    return nullptr;
  }
  return new CoordSet(*src);
}


/*========================================================================*/
/**
 * Append cs to I.
 *
 * @pre I and cs are non-overlapping
 */
int CoordSetMerge(ObjectMolecule *OM, CoordSet * I, const CoordSet * cs)
{
  assert(OM == I->Obj);

  auto const nIndexOld = I->getNIndex();
  I->setNIndex(I->getNIndex() + cs->getNIndex());

  for (int idx_src = 0; idx_src < cs->getNIndex(); ++idx_src) {
    int const idx = idx_src + nIndexOld;
    int const atm = cs->IdxToAtm[idx_src];
    I->IdxToAtm[idx] = cs->IdxToAtm[idx_src];
    if (OM->DiscreteFlag) {
      OM->DiscreteAtmToIdx[atm] = idx;
      OM->DiscreteCSet[atm] = I;
    } else {
      I->AtmToIdx[atm] = idx;
    }
    copy3f(cs->coordPtr(idx_src), I->coordPtr(idx));
  }

  if (cs->RefPos) {
    I->RefPos.resize(I->getNIndex());
    std::copy_n(cs->RefPos.data(), cs->getNIndex(), I->RefPos.data() + nIndexOld);
  }

  // TODO copy or move atom_state_setting_id here?

  I->invalidateRep(cRepAll, cRepInvAll);

  return true;
}


/*========================================================================*/
int CoordSetTransformAtomTTTf(CoordSet * I, int at, const float *TTT)
{
  int a1 = I->atmToIdx(at);
  if(a1 < 0)
    return false;

  auto* v1 = I->coordPtr(a1);
  MatrixTransformTTTfN3f(1, v1, TTT, v1);
  return true;
}


/*========================================================================*/
int CoordSetTransformAtomR44f(CoordSet * I, int at, const float *matrix)
{
  int a1 = I->atmToIdx(at);
  if(a1 < 0)
    return false;

  auto* v1 = I->coordPtr(a1);
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

  ObjectStateLeftCombineMatrixR44d(I, temp);
}


/*========================================================================*/
int CoordSetMoveAtom(CoordSet * I, int at, const float *v, int mode)
{
  int a1 = I->atmToIdx(at);
  if(a1 < 0)
    return false;

  auto* v1 = I->coordPtr(a1);
  if(mode) {
    add3f(v, v1, v1);
  } else {
    copy3f(v, v1);
  }
  return true;
}


/*========================================================================*/

/**
 * Retrieve the label offset for a given atom
 * @param I target Coordset
 * @param atm atom index
 */

pymol::Result<pymol::Vec3> CoordSet::getAtomLabelOffset(int atm) const
{
  int idx = this->atmToIdx(atm);
  if (idx < 0) {
    return pymol::make_error("Invalid atom Idx");
  }
  auto G = this->G;
  auto obj = this->Obj;
  pymol::Vec3 result{};
  const float* at_offset_ptr;
  auto ai = &obj->AtomInfo[atm];
  int at_label_relative_mode = 0;
  AtomStateGetSetting_i(G, obj, this, idx, ai, cSetting_label_relative_mode,
      &at_label_relative_mode);
  switch (at_label_relative_mode) {
  case cLabelRelativeMode::Default:
    AtomStateGetSetting(
        G, obj, this, idx, ai, cSetting_label_placement_offset, &at_offset_ptr);
    break;
  case cLabelRelativeMode::ScreenRelative:
  case cLabelRelativeMode::ScreenPixelSpace:
    AtomStateGetSetting(
        G, obj, this, idx, ai, cSetting_label_screen_point, &at_offset_ptr);
    break;
  }
  std::copy_n(at_offset_ptr, result.size(), result.data());
  return result;
}

/**
 * Retrieve the label offset for a given atom
 * @param I target Coordset
 * @param atm atom index
 * @param offset offset to be set at
 */

pymol::Result<> CoordSet::setAtomLabelOffset(
    int atm, const float* offset)
{
  int idx = this->atmToIdx(atm);
  if (idx < 0) {
    return pymol::make_error("Invalid atom Idx");
  }
  auto G = this->G;
  auto obj = this->Obj;
  auto ai = &obj->AtomInfo[atm];
  int at_label_relative_mode = 0;
  AtomStateGetSetting_i(G, obj, this, idx, ai, cSetting_label_relative_mode,
      &at_label_relative_mode);
  switch (at_label_relative_mode) {
  case cLabelRelativeMode::Default:
    SettingSet(cSetting_label_placement_offset, offset, this, idx);
  case cLabelRelativeMode::ScreenRelative:
  case cLabelRelativeMode::ScreenPixelSpace:
    SettingSet(cSetting_label_screen_point, offset, this, idx);
    break;
  }
  return {};
}

void CoordSet::setTitle(pymol::zstring_view title)
{
  UtilNCopy(Name, title.c_str(), sizeof(WordType));
}

/**
 * Get the transformation matrix which is pre-multiplied to the coordiantes, or
 * NULL if there is either no matrix set, or it's not pre-multiplied
 * (matrix_mode=1).
 */
double const* CoordSet::getPremultipliedMatrix() const {
  return (SettingGet<int>(*this, cSetting_matrix_mode) > 0)
             ? nullptr
             : ObjectStateGetMatrix(this);
}

int CoordSetMoveAtomLabel(CoordSet * I, int at, const float *v, const float *diff)
{
  auto G = I->G;
  ObjectMolecule *obj = I->Obj;
  int a1 = I->atmToIdx(at);
  int result = 0;

  /* if label is valid, get the label offset
   * and set the new position relative to that */
  if(a1 >= 0) {
    float at_offset[3];
    const float * at_offset_ptr;
    int at_label_relative_mode = 0;
    AtomInfoType *ai = obj->AtomInfo + at;

    AtomStateGetSetting_i(G, obj, I, a1, ai, cSetting_label_relative_mode, &at_label_relative_mode);
    switch (at_label_relative_mode){
    case cLabelRelativeMode::Default:
      AtomStateGetSetting(G, obj, I, a1, ai, cSetting_label_placement_offset, &at_offset_ptr);
      add3f(v, at_offset_ptr, at_offset);
      SettingSet(cSetting_label_placement_offset, at_offset, I, a1);
      break;
    case cLabelRelativeMode::ScreenRelative: // screen relative
    case cLabelRelativeMode::ScreenPixelSpace: // screen pixel space
      {
	float voff[3];
	int width, height;
	SceneGetWidthHeight(G, &width, &height);
	if (at_label_relative_mode==1){
	  voff[0] = 2.f * diff[0] / width;
	  voff[1] = 2.f * diff[1] / height;
	} else {
	  voff[0] = diff[0];
	  voff[1] = diff[1];
	}
	voff[2] = 0.f;
	AtomStateGetSetting(G, obj, I, a1, ai, cSetting_label_screen_point, &at_offset_ptr);
	add3f(voff, at_offset_ptr, at_offset);
	SettingSet(cSetting_label_screen_point, at_offset, I, a1);
      }
      break;
    }
  }

  return (result);
}


/*========================================================================*/
int CoordSetGetAtomVertex(const CoordSet * I, int at, float *v)
{
  int a1 = I->atmToIdx(at);

  if(a1 < 0)
    return false;

  copy3f(I->coordPtr(a1), v);
  return true;
}


/*========================================================================*/
int CoordSetGetAtomTxfVertex(const CoordSet * I, int at, float *v)
{
  ObjectMolecule *obj = I->Obj;
  int a1 = I->atmToIdx(at);

  if(a1 < 0)
    return false;

  copy3f(I->coordPtr(a1), v);

  /* apply state transformation */
  if (!I->Matrix.empty() && SettingGet<int>(*I, cSetting_matrix_mode) > 0) {
    transform44d3f(I->Matrix.data(), v, v);
  }

  /* object transformation */
  if(obj->TTTFlag) {
    transformTTT44f3f(obj->TTT, v, v);
  }

  return true;
}


/*========================================================================*/
int CoordSetSetAtomVertex(CoordSet * I, int at, const float *v)
{
  int a1 = I->atmToIdx(at);

  if(a1 < 0)
   return false;

  copy3f(v, I->coordPtr(a1));
  return true;
}


/*========================================================================*/
void CoordSetRealToFrac(CoordSet * I, const CCrystal * cryst)
{
  if (I->getPremultipliedMatrix()) {
    float mat[16];
    copy44d44f(ObjectStateGetInvMatrix(I), mat);
    CoordSetTransform44f(I, mat);
  }

  CoordSetTransform33f(I, cryst->realToFrac());
}


/*========================================================================*/
void CoordSetTransform44f(CoordSet * I, const float *mat)
{
  int a;
  float* v = I->Coord.data();
  for(a = 0; a < I->NIndex; a++) {
    transform44f3f(mat, v, v);
    v += 3;
  }
}


/*========================================================================*/

void CoordSetTransform33f(CoordSet * I, const float *mat)
{
  int a;
  float* v = I->Coord.data();
  for(a = 0; a < I->NIndex; a++) {
    transform33f3f(mat, v, v);
    v += 3;
  }
}


/*========================================================================*/
void CoordSetGetAverage(const CoordSet * I, float *v0)
{
  int a;
  double accum[3];
  if(I->NIndex) {
    const float* v = I->Coord.data();
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
  CoordSetTransform33f(I, cryst->fracToReal());

  if (double const* statemat = I->getPremultipliedMatrix()) {
    float mat[16];
    copy44d44f(statemat, mat);
    CoordSetTransform44f(I, mat);
  }
}

/**
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
  if (!SettingGet<bool>(G, cSetting_pdb_insure_orthogonal))
    return false;

  if (!cryst)
    cryst = &cset->Symmetry->Crystal;

  auto const r2f = cryst->realToFrac();

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
/* Rotates the ANISOU vector
 *
 * matrix: flat 4x4, but only rotation (upper left 3x3) is considered
 * anisou: has 6 elements (of symmetric 3x3) and will be rotated in-place
 */
bool RotateU(const double *matrix, float *anisou)
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
/**
 * @param v 3x1 vertex in final output space
 * @param matrix 4x4 homogenous transformation matrix from model space to output
 *         space (view matrix * state matrix). Used for ANISOU.
 */
void CoordSetAtomToPDBStrVLA(PyMOLGlobals * G, char **charVLA, int *c,
                             const AtomInfoType * ai,
                             const float *v, int cnt,
                             const PDBInfoRec * pdb_info,
                             const double *matrix)
{
  char *aType;
  AtomName name;
  ResName resn;
  lexidx_t chain;

  char formalCharge[4];
  bool const ignore_pdb_segi = SettingGet<bool>(G, cSetting_ignore_pdb_segi);
  WordType x, y, z;

  AtomInfoGetAlignedPDBResidueName(G, ai, resn);
  AtomInfoGetAlignedPDBAtomName(G, ai, resn, name);

  formalCharge[0] = 0;
  if (SettingGet<bool>(G, cSetting_pdb_formal_charges)) {
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

  if (SettingGet<bool>(G, cSetting_pdb_retain_ids)) {
    cnt = ai->id - 1;
  }
  if(cnt > 99998)
    cnt = 99998;

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
      std::copy_n(ai->anisou, 6, anisou);

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
    if (pymol::zstring_view(resn).find_first_not_of(' ') ==
        pymol::zstring_view::npos) {
      // APBS fails with empty residue names
      assert(resn[0] == ' ');
      resn[0] = '.';
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
      std::copy_n(ai->anisou, 6, tmp_array);
      if (matrix)
      RotateU(matrix, tmp_array);
    }

    PConvFloat3ToPyObjAttr(atom, "coord", v);
    if(ref)
      PConvFloat3ToPyObjAttr(atom, "ref_coord", ref);
    if (ai->name)
    PConvStringToPyObjAttr(atom, "name", LexStr(G, ai->name));
    PConvStringToPyObjAttr(atom, "symbol", ai->elem);
// TODO defaults to UNK    if (ai->resn)
    PConvStringToPyObjAttr(atom, "resn", LexStr(G, ai->resn));
    if (ai->inscode) {
      char ins_code[2] = {ai->inscode, 0};
      PConvStringToPyObjAttr(atom, "ins_code", ins_code);
    }
    if (ai->ssType[0])
    PConvStringToPyObjAttr(atom, "ss", ai->ssType);
// TODO defaults to 1    if (ai->resv)
    PConvIntToPyObjAttr(atom, "resi_number", ai->resv);
    if (ai->stereo)
    PConvIntToPyObjAttr(atom, "stereo", ai->stereo);
    if (ai->chain)
    PConvStringToPyObjAttr(atom, "chain", LexStr(G, ai->chain));
    if(ai->alt[0])
      PConvStringToPyObjAttr(atom, "alt", ai->alt);
    if (ai->segi)
    PConvStringToPyObjAttr(atom, "segi", LexStr(G, ai->segi));
    if (ai->q != 1.)
    PConvFloatToPyObjAttr(atom, "q", ai->q);
    if (ai->b)
    PConvFloatToPyObjAttr(atom, "b", ai->b);
    if (ai->anisou)
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
    if (ai->elec_radius)
    PConvFloatToPyObjAttr(atom, "elec_radius", ai->elec_radius);
    if (ai->partialCharge)
    PConvFloatToPyObjAttr(atom, "partial_charge", ai->partialCharge);
    if (ai->formalCharge)
    PConvIntToPyObjAttr(atom, "formal_charge", ai->formalCharge);
// TODO customType=0 from most files
    if(ai->customType != -9999)
      PConvIntToPyObjAttr(atom, "numeric_type", ai->customType);
    if (ai->textType)
    PConvStringToPyObjAttr(atom, "text_type", LexStr(G, ai->textType));
    if (ai->custom)
    PConvStringToPyObjAttr(atom, "custom", LexStr(G, ai->custom));

    PConvIntToPyObjAttr(atom, "hetatm", ai->hetatm);
    PConvIntToPyObjAttr(atom, "flags", ai->flags);
    PConvIntToPyObjAttr(atom, "id", ai->id);    /* not necc. unique */
    PConvIntToPyObjAttr(atom, "index", index + 1);      /* fragile */

#ifdef _PYMOL_IP_PROPERTIES
#endif

  }
  if(PyErr_Occurred())
    PyErr_Print();
  return (atom);
#endif
}


/*========================================================================*/
void CoordSet::invalidateRep(cRep_t type, cRepInv_t level)
{
  if(level >= cRepInvVisib) {
    if (Obj)
      Obj->RepVisCacheValid = false;
  }
  /* graphical representations need redrawing */
  if(level == cRepInvVisib) {
    /* cartoon_side_chain_helper */
    if (SettingGet<bool>(*this, cSetting_cartoon_side_chain_helper)) {
      if((type == cRepCyl) || (type == cRepLine) || (type == cRepSphere))
        invalidateRep(cRepCartoon, cRepInvVisib2);
      else if(type == cRepCartoon) {
        invalidateRep(cRepLine, cRepInvVisib2);
        invalidateRep(cRepCyl, cRepInvVisib2);
        invalidateRep(cRepSphere, cRepInvVisib2);
      }
    }
    /* ribbon_side_chain_helper */
    if (SettingGet<bool>(*this, cSetting_ribbon_side_chain_helper)) {
      if((type == cRepCyl) || (type == cRepLine) || (type == cRepSphere))
        invalidateRep(cRepRibbon, cRepInvVisib2);
      else if(type == cRepRibbon) {
        invalidateRep(cRepLine, cRepInvVisib2);
        invalidateRep(cRepCyl, cRepInvVisib2);
        invalidateRep(cRepSphere, cRepInvVisib2);
      }
    }
    /* line_stick helper  */
    if (SettingGet<bool>(*this, cSetting_line_stick_helper)) {
      if(type == cRepCyl)
        invalidateRep(cRepLine, cRepInvVisib2);
      else if(type == cRepLine) {
        invalidateRep(cRepCyl, cRepInvVisib2);
      }
    }
  }

  /* invalidate basd on one representation, 'type' */
  for (RepIterator iter(G, type); iter.next(); ){
    auto eff_level = level;
    auto const a = iter.rep;
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
      Active[a] = true;
    if (Rep[a]) {
      if (eff_level < cRepInvPurge) {
        Rep[a]->invalidate(eff_level);
      } else {
        delete Rep[a];
        Rep[a] = nullptr;
      }
    }
  }

  if(level >= cRepInvCoord) {   /* if coordinates change, then this map becomes invalid */
    MapFree(Coord2Idx);
    Coord2Idx = nullptr;
    ExecutiveInvalidateSelectionIndicatorsCGO(G);
    SceneInvalidatePicking(G);
    /* invalidate distances */
  }

#ifndef NO_MMLIBS
  if (level >= cRepInvProp) {
    if (validMMStereo == MMPYMOLX_PROP_STATE_AUTO)
      validMMStereo = MMPYMOLX_PROP_STATE_NULL;
    if (validTextType == MMPYMOLX_PROP_STATE_AUTO)
      validTextType = MMPYMOLX_PROP_STATE_NULL;
  }
#endif

  SceneChanged(G);
}


/*========================================================================*/

#define RepUpdateMacro(rep, new_fn, state)                                     \
  {                                                                            \
    if (Active[rep] && !G->Interrupt) {                                        \
      if (Rep[rep]) {                                                          \
        assert(Rep[rep]->cs == this);                                          \
        assert(Rep[rep]->getState() == state);                                 \
        Rep[rep] = Rep[rep]->update();                                         \
      } else {                                                                 \
        Rep[rep] = new_fn(this, state);                                        \
        if (Rep[rep]) {                                                        \
          Rep[rep]->fNew = new_fn;                                             \
          SceneInvalidatePicking(G);                                           \
        } else {                                                               \
          Active[rep] = false;                                                 \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    OrthoBusyFast(G, rep, cRepCnt);                                            \
  }

/*========================================================================*/
void CoordSet::update(int state)
{
  assert(G == Obj->G);

  OrthoBusyFast(G, 0, cRepCnt);
  RepUpdateMacro(cRepLine, RepWireBondNew, state);
  RepUpdateMacro(cRepCyl, RepCylBondNew, state);
  RepUpdateMacro(cRepDot, RepDotNew, state);
  RepUpdateMacro(cRepMesh, RepMeshNew, state);
  RepUpdateMacro(cRepSphere, RepSphereNew, state);
  RepUpdateMacro(cRepRibbon, RepRibbonNew, state);
  RepUpdateMacro(cRepCartoon, RepCartoonNew, state);
  RepUpdateMacro(cRepSurface, RepSurfaceNew, state);
  RepUpdateMacro(cRepLabel, RepLabelNew, state);
  RepUpdateMacro(cRepNonbonded, RepNonbondedNew, state);
  RepUpdateMacro(cRepNonbondedSphere, RepNonbondedSphereNew, state);
  RepUpdateMacro(cRepEllipsoid, RepEllipsoidNew, state);

  for (int a = 0; a < cRepCnt; ++a) {
    if (!Rep[a])
      Active[a] = false;
  }

  // cell
  if ((Obj->visRep & cRepCellBit) && !UnitCellCGO) {
    if (auto const* sym = getSymmetry()) {
      UnitCellCGO.reset(CrystalGetUnitCellCGO(&sym->Crystal));
    }
  }

  SceneInvalidate(G);
  OrthoBusyFast(G, 1, 1);
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
      I->Coord2Idx = MapNew(I->G, I->Coord2IdxDiv, I->Coord, I->NIndex, NULL);
      if(I->Coord2IdxDiv < I->Coord2Idx->Div)
        I->Coord2IdxDiv = I->Coord2Idx->Div;
    }
  }
}

/*========================================================================*/
void CoordSet::render(RenderInfo * info)
{
  const auto pass = info->pass;
  const auto pick = bool(info->pick);
  auto* ray = info->ray;

  if (!(ray || pick) &&
      (SettingGet<int>(*this, cSetting_defer_builds_mode) == 5)) {
    if (pass == RenderPass::Antialias) {
      ObjectUseColor(Obj);
      if (Active[cRepLine])
        RepWireBondRenderImmediate(this, info);
      if (Active[cRepNonbonded])
        RepNonbondedRenderImmediate(this, info);
      if (Active[cRepSphere])
        RepSphereRenderImmediate(this, info);
      if (Active[cRepCyl])
        RepCylBondRenderImmediate(this, info);
      if (Active[cRepRibbon])
        RepRibbonRenderImmediate(this, info);
    }

    return;
  }

  const auto sculpt_vdw_vis_mode =
      SettingGet<int>(*this, cSetting_sculpt_vdw_vis_mode);

  if (pass == RenderPass::Antialias && sculpt_vdw_vis_mode && SculptCGO &&
      (Obj->visRep & cRepCGOBit)) {
    if (ray) {
      CGORenderRay(SculptCGO, ray, info, ColorGet(G, Obj->Color), nullptr,
          Setting.get(), Obj->Setting.get());
    } else if (G->HaveGUI && G->ValidContext && !pick) {
      if (!info->use_shaders) {
        CGOFree(SculptShaderCGO);
      } else if (!SculptShaderCGO) {
        SculptShaderCGO = CGOOptimizeToVBONotIndexed(SculptCGO);
      }

      auto* cgo = SculptShaderCGO ? SculptShaderCGO : SculptCGO;
      CGORenderGL(
          cgo, nullptr, Setting.get(), Obj->Setting.get(), info, nullptr);
    }
  }

  // cell
  if (UnitCellCGO && (Obj->visRep & cRepCellBit)) {
    if (ray) {
      CGORenderRay(UnitCellCGO.get(), ray, info, ColorGet(G, Obj->Color),
          nullptr, Setting.get(), Obj->Setting.get());
    } else if (!pick && pass == RenderPass::Opaque && G->HaveGUI &&
               G->ValidContext) {
      ObjectUseColor(Obj);
      CGORenderGL(UnitCellCGO.get(), ColorGet(G, Obj->Color), Setting.get(),
          Obj->Setting.get(), info, nullptr);
    }
  }

  auto aastart = cRep_t(0), aaend = cRepCnt;

  if (pick) {
    int pick_labels = SettingGet<int>(*this, cSetting_pick_labels);
    if (pick_labels == 2) { // only pick labels
      aastart = cRepLabel;
      aaend = cRep_t(aastart + 1);
    }
  }

  for (auto a = aastart; a != aaend; ++a) {
    auto* const r = Rep[a];

    if (!(r && Active[a])) {
      continue;
    }

    if (!ray) {
      ObjectUseColor(Obj);
    } else {
      ray->wobble(SettingGet<int>(*this, cSetting_ray_texture),
          SettingGet<const float*>(*this, cSetting_ray_texture_settings));
      ray->color3fv(ColorGet(G, Obj->Color));
    }

    if (ray || pick) {
      r->render(info);
      continue;
    }

    // OIT transparency mode
    bool const t_mode_3 = (3 == SettingGet<int>(G, cSetting_transparency_mode));

    // render both opaque and transparent pass
    bool render_both = t_mode_3;

    switch (a) {
    case cRepVolume: {
      if (t_mode_3 && pass == RenderPass::Transparent) {
        r->render(info);
      }
    } break;
    case cRepLabel:
      if (pass == RenderPass::Transparent ||
          (t_mode_3 && pass == RenderPass::Opaque /* TODO_OPENVR 0 */)) {
        r->render(info);
      }
      break;
    case cRepDot:
    case cRepCGO:
    case cRepCallback:
      if (pass == RenderPass::Opaque) {
        r->render(info);
      }
      break;
    case cRepCell:
      // not implemented
      assert(false);
    case cRepLine:
    case cRepMesh:
    case cRepDash:
    case cRepExtent:
      if (pass == RenderPass::Antialias) {
        r->render(info);
      }
      break;
    case cRepNonbonded:
    case cRepRibbon:
      render_both = false; // cartoon and nonbonded do not have atom-level
                           // transparency
    case cRepNonbondedSphere:
    case cRepSurface:
    case cRepEllipsoid:
    case cRepCyl:
    case cRepCartoon:
    case cRepSphere:
      if (render_both) {
        if (pass != RenderPass::Antialias) {
          r->render(info);
        }
      } else if (info->alpha_cgo && a == cRepSurface) {
        if (pass == RenderPass::Opaque) {
          r->render(info);
        }
      } else if (r->hasTransparency()) {
        if (pass == RenderPass::Transparent) {
          r->render(info);
        }
      } else if (pass == RenderPass::Opaque) {
        r->render(info);
      }
      break;
    }
  }
}

/*========================================================================*/
CoordSet::CoordSet(PyMOLGlobals* G)
    : CObjectState(G)
{
}

/*========================================================================*/
CoordSet::CoordSet(const CoordSet& cs)
    : CObjectState(cs)
{
  this->Obj = cs.Obj;
  this->Coord = cs.Coord;
  this->IdxToAtm = cs.IdxToAtm;
  this->NIndex = cs.NIndex;
  std::copy(std::begin(cs.Rep), std::end(cs.Rep), std::begin(this->Rep));
  std::copy(std::begin(cs.Active), std::end(cs.Active), std::begin(this->Active));
  this->NTmpBond = cs.NTmpBond;
  this->NTmpLinkBond = cs.NTmpLinkBond;
  /* deep copy & return ptr to new symmetry */
  if(cs.Symmetry) {
    this->Symmetry = pymol::make_unique<CSymmetry>(*cs.Symmetry);
  }
  std::copy(std::begin(cs.Name), std::end(cs.Name), std::begin(this->Name));
  this->PeriodicBoxType = cs.PeriodicBoxType;
  this->tmp_index = cs.tmp_index;
  this->Coord2IdxReq = cs.Coord2IdxReq;
  this->Coord2IdxDiv = cs.Coord2IdxDiv;
  this->objMolOpInvalidated = cs.objMolOpInvalidated;

  // copy VLAs
  this->RefPos     = cs.RefPos;
  this->AtmToIdx   = cs.AtmToIdx;

  UtilZeroMem(this->Rep, sizeof(::Rep *) * cRepCnt);

#ifdef _PYMOL_IP_PROPERTIES
#endif

}


/*========================================================================*/
/**
 * Sets the number of atoms with coordinates.
 *
 * @post NIndex updated
 * @post All relevant arrays resized
 */
void CoordSet::setNIndex(unsigned nindex)
{
  NIndex = nindex;
  IdxToAtm.resize(nindex);

  if (nindex == 0) {
    return;
  }

  auto const idx = nindex - 1;

  Coord.reserve(nindex * 3);

  if (has_any_atom_state_settings()) {
    atom_state_setting_id.check(idx);
  }

  if (RefPos) {
    RefPos.check(idx);
  }
}

/**
 * Sets the number of atoms.
 *
 * For non-discrete objects:
 *   - Extends AtmToIdx (appends -1s)
 *
 * For discrete objects:
 *   - Converts this coordset to discrete if necessary
 *   - Calls ObjectMolecule::setNDiscrete()
 *
 * @pre NIndex and IdxToAtm are valid
 */
int CoordSet::extendIndices(int nAtom)
{
  bool ok = true;
  if (Obj->DiscreteFlag) {
    ok = Obj->setNDiscrete(nAtom);

    // convert to discrete if necessary
    if (!AtmToIdx.empty()) {
      AtmToIdx.clear();
      if (ok) {
        for (int a = 0; a < NIndex; a++) {
          auto const b = IdxToAtm[a];
          Obj->DiscreteAtmToIdx[b] = a;
          Obj->DiscreteCSet[b] = this;
        }
      }
    }
  } else {
    const auto NAtIndex = AtmToIdx.size();
    assert(NAtIndex <= nAtom);
    if (NAtIndex < nAtom) {
      AtmToIdx.resize(nAtom);
      if (ok && nAtom) {
        for (int a = NAtIndex; a < nAtom; a++)
          AtmToIdx[a] = -1;
      }
    }
  }
  return ok;
}

/*========================================================================*/
/**
 * Set up for simple case where 1 = 1, etc.
 */
void CoordSet::enumIndices()
{
  AtmToIdx.resize(NIndex);
  IdxToAtm.resize(NIndex);
  if (NIndex) {
    for (int a = 0; a < NIndex; ++a) {
      AtmToIdx[a] = a;
      IdxToAtm[a] = a;
    }
  }
}

/*========================================================================*/
CoordSet::~CoordSet()
{
#ifdef _PYMOL_IP_PROPERTIES
#endif

  if (has_any_atom_state_settings()) {
    for (int a = 0; a < NIndex; ++a) {
      if (has_atom_state_settings(a)) {
        SettingUniqueDetachChain(G, atom_state_setting_id[a]);
      }
    }
  }

  for (int a = 0; a < cRepCnt; ++a) {
    delete Rep[a];
  }

  MapFree(Coord2Idx);
  CGOFree(SculptCGO);
  CGOFree(SculptShaderCGO);
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

#ifndef _PYMOL_NOPY
int CoordSetSetSettingFromPyObject(PyMOLGlobals * G, CoordSet *cs, int at, int setting_id, PyObject *val){
  if (val == Py_None)
    val = NULL;

  if (!val && !cs->has_atom_state_settings(at)) {
    return true;
  }

  CoordSetCheckUniqueID(G, cs, at);

  return SettingUniqueSetPyObject(G, cs->atom_state_setting_id[at], setting_id, val);
}
#endif

int CoordSetCheckSetting(PyMOLGlobals * G, CoordSet *cs, int at, int setting_id){
  if (!cs->has_atom_state_settings(at)) {
    return 0;
  }

  return SettingUniqueCheck(G, cs->atom_state_setting_id[at], setting_id);
}

PyObject *SettingGetIfDefinedPyObject(PyMOLGlobals * G, CoordSet *cs, int at, int setting_id){
  if (cs->has_atom_state_settings(at)) {
    return SettingUniqueGetPyObject(G, cs->atom_state_setting_id[at], setting_id);
  }
  return NULL;
}

int CoordSetCheckUniqueID(PyMOLGlobals * G, CoordSet *I, int at){
  if (!I->atom_state_setting_id){
    I->atom_state_setting_id = pymol::vla<int>(I->NIndex);
  }
  if (!I->atom_state_setting_id[at]){
    I->atom_state_setting_id[at] = AtomInfoGetNewUniqueID(G);
  }
  return I->atom_state_setting_id[at];
}

template <typename V>
void AtomStateGetSetting(ATOMSTATEGETSETTINGARGS, V * out) {
  if (cs->has_atom_state_settings(idx) &&
      SettingUniqueGetIfDefined(G, cs->atom_state_setting_id[idx], setting_id, out))
    return;

  if (AtomSettingGetIfDefined(G, ai, setting_id, out))
    return;

  *out = SettingGet<V>(*cs, setting_id);
}

template void AtomStateGetSetting(ATOMSTATEGETSETTINGARGS, int * out);
template void AtomStateGetSetting(ATOMSTATEGETSETTINGARGS, float * out);
template void AtomStateGetSetting(ATOMSTATEGETSETTINGARGS, const float ** out);
