
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
#include"os_std.h"

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"DistSet.h"
#include"Scene.h"
#include"Color.h"
#include"RepDistDash.h"
#include"RepDistLabel.h"
#include"RepAngle.h"
#include"RepDihedral.h"
#include"PConv.h"
#include"ObjectDist.h"
#include"ObjectMolecule.h"
#include"ListMacros.h"
#include"Selector.h"
#include "PyMOL.h"
#include "Executive.h"

int DistSetGetLabelVertex(DistSet * I, int at, float *v)
{
  if (at >= 0 && at < I->LabCoord.size()) {
    const auto& vv = I->LabCoord[at];
    copy3f(vv.data(), v);
    return true;
  }
  return false;
}

/**
 * Retrieves label offset for given distance set
 * @param atm atom index
 */

pymol::Result<pymol::Vec3> DistSet::getLabelOffset(int atm) const
{
  std::array<float, 3> result;
  if (atm < 0 || atm >= this->LabPos.size()) {
    return pymol::make_error("Invalid index");
  }

  auto lp = &this->LabPos[atm];
  if (lp->mode) {
    std::copy_n(lp->offset, result.size(), result.data());
    return result;
  }

  auto obj = this->Obj;
  auto G = obj->G;
  auto lab_pos =
      SettingGet_3fv(G, nullptr, obj->Setting.get(), cSetting_label_position);
  std::copy_n(lab_pos, result.size(), result.data());
  return result;
}

/**
 * Retrieves label offset for given distance set
 * @param atm atom index
 * @param pos offset to set the label at
 */

pymol::Result<> DistSet::setLabelOffset(int atm, const float* pos)
{
  if (atm < 0) {
    return pymol::make_error("Invalid index");
  }
  VecCheck(this->LabPos, atm);
  auto& lp = this->LabPos[atm];
  lp.mode = 1;
  copy3f(pos, lp.offset);
  return {};
}

int DistSetMoveLabel(DistSet * I, int a1, float *v, int mode)
{
  if(a1 < 0) {
    return 0;
  }

  VecCheck(I->LabPos, a1);
  auto& lp = I->LabPos[a1];
  if(!lp.mode) {
    auto obj = I->Obj;
    const float *lab_pos = SettingGet_3fv(obj->G, NULL, obj->Setting.get(),
                                    cSetting_label_position);
    copy3f(lab_pos, lp.pos);
  }
  lp.mode = 1;
  if(mode) {
    add3f(v, lp.offset, lp.offset);
  } else {
    copy3f(v, lp.offset);
  }

  return 1;
}


/**
 * @param I   measurement set, must not be NULL
 * @param obj object molecule, can be NULL so then all items in I will be updated
 * @return  number of updated coordinates
 */
int DistSetMoveWithObject(DistSet * I, struct ObjectMolecule *obj)
{
  PyMOLGlobals * G = I->G;

  int i, N, rVal = 0;
  float * varDst;

  PRINTFD(G, FB_DistSet)
    " DistSet: adjusting distance vertex\n" ENDFD;

  for (auto& listitem : I->MeasureInfo) {
    auto* memb = &listitem;
    varDst = NULL;

    switch(memb->measureType) {
    case cRepDash:
      N = 2;
      if(memb->offset < I->NIndex + 1)
        varDst = I->Coord.data();
      break;
    case cRepAngle:
      N = 3;
      if(memb->offset < I->NAngleIndex + 2)
        varDst = I->AngleCoord.data();
      break;
    case cRepDihedral:
      N = 4;
      if(memb->offset < I->NDihedralIndex + 3)
        varDst = I->DihedralCoord.data();
      break;
    }

    if(!varDst)
      continue;

    varDst += 3 * memb->offset;

    for(i = 0; i < N; i++) {
      auto eoo = ExecutiveUniqueIDAtomDictGet(G, memb->id[i]);

      if(!eoo || (obj && obj != eoo->obj))
        continue;

      if(ObjectMoleculeGetAtomVertex(
            eoo->obj, memb->state[i],
            eoo->atm, varDst + i * 3))
        rVal++;
    }
  }

  if (rVal)
    I->invalidateRep(cRepAll, cRepInvCoord);

  PRINTFD(G, FB_DistSet)
    " DistSet: done updating distance set's vertex\n" ENDFD;

  return rVal;
}

static decltype(DistSet::MeasureInfo) MeasureInfoListFromPyList(
    PyMOLGlobals* G, PyObject* list)
{
  int i, ll, N;
  decltype(DistSet::MeasureInfo) I;
  CPythonVal *val, *tmp;

  ok_assert(1, list && PyList_Check(list));
  ll = PyList_Size(list);

  for (i = 0; i < ll; i++) {
    I.emplace_front();
    auto* item = &I.front();

    val = CPythonVal_PyList_GetItem(G, list, i);
    if(val && PyList_Check(val) &&
              PyList_Size(val) > 2) {

      tmp = CPythonVal_PyList_GetItem(G, val, 1);
      N = PyList_Size(tmp);
      ok_assert(1, N < 5);

      item->measureType = (N == 2) ? cRepDash :
                          (N == 3) ? cRepAngle : cRepDihedral;

      CPythonVal_PConvPyIntToInt_From_List(G, val, 0, &item->offset);
      CPythonVal_PConvPyListToIntArrayInPlace(G, tmp, item->id, N);
      CPythonVal_PConvPyListToIntArrayInPlace_From_List(G, val, 2, item->state, N);
      CPythonVal_Free(tmp);

      for (int j = 0; j < N; ++j) {
        item->id[j] = SettingUniqueConvertOldSessionID(G, item->id[j]);
      }
    }
    CPythonVal_Free(val);
  }

ok_except1:
  return I;
}

static PyObject* MeasureInfoListAsPyList(
    decltype(DistSet::MeasureInfo) const& list)
{
  int N;
  PyObject *item, *result = PyList_New(0);
  ok_assert(1, result);

  for (auto& listitem : list) {
    auto const* I = &listitem;
    switch(I->measureType) {
      case cRepDash: N = 2; break;
      case cRepAngle: N = 3; break;
      default: N = 4;
    }

    ok_assert(1, item = PyList_New(3));

    PyList_SetItem(item, 0, PyInt_FromLong(I->offset));
    PyList_SetItem(item, 1, PConvIntArrayToPyList(I->id, N));
    PyList_SetItem(item, 2, PConvIntArrayToPyList(I->state, N));

    PyList_Append(result, item);
    Py_DECREF(item);
  }

ok_except1:
  return PConvAutoNone(result);
}

DistSet* DistSetFromPyList(PyMOLGlobals * G, PyObject * list)
{
  DistSet *I = NULL;
  int ll = 0;
  CPythonVal *val;

  if(CPythonVal_IsNone(list)) {         /* allow None for CSet */
    return nullptr;
  }

  ok_assert(1, list && PyList_Check(list));
  ok_assert(1, I = DistSetNew(G));

  ll = PyList_Size(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */

  ok_assert(1, CPythonVal_PConvPyIntToInt_From_List(G, list, 0, &I->NIndex));
  ok_assert(1, CPythonVal_PConvPyListToFloatVLANoneOkay_From_List(G, list, 1, &I->Coord));

  ok_assert(2, ll > 2);

  ok_assert(1, CPythonVal_PConvPyIntToInt_From_List(G, list, 3, &I->NAngleIndex));
  ok_assert(1, CPythonVal_PConvPyListToFloatVLANoneOkay_From_List(G, list, 4, &I->AngleCoord));
  ok_assert(1, CPythonVal_PConvPyIntToInt_From_List(G, list, 5, &I->NDihedralIndex));
  ok_assert(1, CPythonVal_PConvPyListToFloatVLANoneOkay_From_List(G, list, 6, &I->DihedralCoord));

  ok_assert(2, ll > 7);

  // DistSet->Setting never gets set (removed BB 11/14), was state settings?
  /*  val = CPythonVal_PyList_GetItem(G, list, 7);
      I->Setting = SettingNewFromPyList(G, val);
      CPythonVal_Free(val); */

  ok_assert(2, ll > 8);

  val = CPythonVal_PyList_GetItem(G, list, 8);
  {
    auto labPosRes = CPythonVal_PConvPyListToLabPosVec(G, val);
    ok_assert(1, (bool)labPosRes);
    I->LabPos = std::move(*labPosRes);
  }
  CPythonVal_Free(val);

  ok_assert(2, ll > 9);

  val = CPythonVal_PyList_GetItem(G, list, 9);
  I->MeasureInfo = MeasureInfoListFromPyList(G, val);
  CPythonVal_Free(val);

ok_except2:
  return I;
ok_except1:
  delete I;
  return nullptr;

}

PyObject *DistSetAsPyList(DistSet * I)
{
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(10);

    PyList_SetItem(result, 0, PyInt_FromLong(I->NIndex));
    PyList_SetItem(result, 1, PConvFloatArrayToPyListNullOkay(I->Coord, I->NIndex * 3));
    PyList_SetItem(result, 2, PConvAutoNone(NULL)); // I->LabCoord recalculated in RepDistLabelNew
    PyList_SetItem(result, 3, PyInt_FromLong(I->NAngleIndex));
    PyList_SetItem(result, 4,
                   PConvFloatArrayToPyListNullOkay(I->AngleCoord, I->NAngleIndex * 3));
    PyList_SetItem(result, 5, PyInt_FromLong(I->NDihedralIndex));
    PyList_SetItem(result, 6,
                   PConvFloatArrayToPyListNullOkay(I->DihedralCoord,
                                                   I->NDihedralIndex * 3));
    // DistSet->Setting never gets set (removed BB 11/14), was state settings?
    PyList_SetItem(result, 7, PConvAutoNone(NULL) /* SettingAsPyList(I->Setting) */);
    if(!I->LabPos.empty()) {
      PyList_SetItem(result, 8, PConvLabPosVecToPyList(I->LabPos));
    } else {
      PyList_SetItem(result, 8, PConvAutoNone(NULL));
    }
    PyList_SetItem(result, 9, MeasureInfoListAsPyList(I->MeasureInfo));
    /* TODO setting ... */
  }
  return (PConvAutoNone(result));
}


/*========================================================================*/
int DistSetGetExtent(DistSet * I, float *mn, float *mx)
{
  int a;
  int c;
  const float* v = I->Coord.data();
  for(a = 0; a < I->NIndex; a++) {
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
  }

  v = I->AngleCoord;
  c = I->NAngleIndex / 5;
  for(a = 0; a < c; a++) {
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 9;
  }

  v = I->DihedralCoord;
  c = I->NDihedralIndex / 6;
  for(a = 0; a < c; a++) {
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 3;
    min3f(v, mn, mn);
    max3f(v, mx, mx);
    v += 9;
  }
  return (I->NIndex + I->NAngleIndex + I->NDihedralIndex);
}


/*========================================================================*/
/**
 * Invalidate reps
 *
 * type: rep enum, e.g. cRepDash
 * level: e.g. cRepInvColor
 */
void DistSet::invalidateRep(cRep_t type, cRepInv_t level)
{
  int a = 0, a_stop = getNRep();
  bool changed = false;

  /* if representation type is specified, adjust it */
  if(type >= 0) {
    if(type >= getNRep())
      return;

    a = type;
    a_stop = a + 1;
  }

  for(; a < a_stop; a++) {
    if(Rep[a]) {
      changed = true;
      Rep[a].reset();
    }
  }

  if (changed)
    SceneChanged(G);
}


/*========================================================================*/
void DistSet::update(int state)
{
  DistSet * I = this;
  /* status bar 0% */
  OrthoBusyFast(G, 0, I->getNRep());
  if(!I->Rep[cRepDash]) {
    /* query the dist set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepDash].reset(RepDistDashNew(I,state));
    SceneInvalidate(G);
  }
  if(!I->Rep[cRepLabel]) {
    /* query the dist set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepLabel].reset(RepDistLabelNew(I, state));
    SceneInvalidate(G);
  }
  if(!I->Rep[cRepAngle]) {
    /* query the angle set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepAngle].reset(RepAngleNew(I, state));
    SceneInvalidate(G);
  }
  if(!I->Rep[cRepDihedral]) {
    /* query the dihedral set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepDihedral].reset(RepDihedralNew(I, state));
    SceneInvalidate(G);
  }
  /* status bar 100% */
  OrthoBusyFast(G, 1, 1);
}


/*========================================================================*/
void DistSet::render(RenderInfo * info)
{
  DistSet * I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  for (int a = 0; a < I->getNRep(); a++) {
    if(!GET_BIT(I->Obj->visRep, a))
      continue;
    if(!I->Rep[a]) {
      switch(a) {
      case cRepDash:
        I->Rep[a].reset(RepDistDashNew(I, -1));
        break;
      case cRepLabel:
        I->Rep[a].reset(RepDistLabelNew(I, -1));
        break;
      case cRepAngle:
        I->Rep[a].reset(RepAngleNew(I, -1));
        break;
      case cRepDihedral:
        I->Rep[a].reset(RepDihedralNew(I, -1));
        break;
      }
    }
    if(I->Rep[a])
      {
        if(ray) {
            ray->color3fv(ColorGet(G, I->Obj->Color));
        } else if (!pick) {
          ObjectUseColor(I->Obj);
        }
        Rep[a]->render(info);
      }
  }
}


/*========================================================================*/
DistSet::DistSet(PyMOLGlobals* G)
    : CObjectState(G)
{
}

