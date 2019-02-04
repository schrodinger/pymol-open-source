
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
#include"ObjectMolecule.h"
#include"ListMacros.h"
#include"Selector.h"
#include "PyMOL.h"
#include "Executive.h"

int DistSetGetLabelVertex(DistSet * I, int at, float *v)
{
  if((at >= 0) && (at < I->NLabel) && I->LabCoord) {
    float *vv = I->LabCoord + 3 * at;
    copy3f(vv, v);
    return true;
  }
  return false;
}

int DistSetMoveLabel(DistSet * I, int at, float *v, int mode)
{
  ObjectDist *obj;
  int a1 = at;
  int result = 0;
  LabPosType *lp;

  obj = I->Obj;

  if(a1 >= 0) {

    if(!I->LabPos)
      I->LabPos = VLACalloc(LabPosType, I->NLabel);
    if(I->LabPos) {
      result = 1;
      lp = I->LabPos + a1;
      if(!lp->mode) {
        const float *lab_pos = SettingGet_3fv(obj->Obj.G, NULL, obj->Obj.Setting,
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


/* -- JV, refactored by TH */
/*
 * PARAMS
 *   I:   measurement set, must not be NULL
 *   obj: object molecule, can be NULL so then all items in I will be updated
 * RETURNS
 *   number of updated coordinates
 */
int DistSetMoveWithObject(DistSet * I, struct ObjectMolecule *obj)
{
  PyMOLGlobals * G = I->State.G;

  int i, N, rVal = 0;
  CMeasureInfo * memb = NULL;
  float * varDst;

  PRINTFD(G, FB_DistSet)
    " DistSet: adjusting distance vertex\n" ENDFD;

  for(memb = I->MeasureInfo; memb; memb = memb->next) {
    varDst = NULL;

    switch(memb->measureType) {
    case cRepDash:
      N = 2;
      if(memb->offset < I->NIndex + 1)
        varDst = I->Coord;
      break;
    case cRepAngle:
      N = 3;
      if(memb->offset < I->NAngleIndex + 2)
        varDst = I->AngleCoord;
      break;
    case cRepDihedral:
      N = 4;
      if(memb->offset < I->NDihedralIndex + 3)
        varDst = I->DihedralCoord;
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
    I->invalidateRep(-1, cRepInvCoord);

  PRINTFD(G, FB_DistSet)
    " DistSet: done updating distance set's vertex\n" ENDFD;

  return rVal;
}

static CMeasureInfo * MeasureInfoListFromPyList(PyMOLGlobals * G, PyObject * list)
{
  int i, ll, N;
  CMeasureInfo *item = NULL, *I= NULL;
  CPythonVal *val, *tmp;

  ok_assert(1, list && PyList_Check(list));
  ll = PyList_Size(list);

  for (i = 0; i < ll; i++) {
    ok_assert(1, item = pymol::malloc<CMeasureInfo>(1));
    ListPrepend(I, item, next);

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

static PyObject *MeasureInfoListAsPyList(CMeasureInfo * I)
{
  int N;
  PyObject *item, *result = PyList_New(0);
  ok_assert(1, result);

  while (I) {
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
    I = I->next;
  }

ok_except1:
  return PConvAutoNone(result);
}

int DistSetFromPyList(PyMOLGlobals * G, PyObject * list, DistSet ** cs)
{
  DistSet *I = NULL;
  int ll = 0;
  CPythonVal *val;

  if(*cs) {
    (*cs)->fFree();
    *cs = NULL;
  }

  if(list == Py_None) {         /* allow None for CSet */
    *cs = NULL;
    return true;
  }

  ok_assert(1, list && PyList_Check(list));
  ok_assert(1, I = DistSetNew(G));

  ll = PyList_Size(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */

  ok_assert(1, CPythonVal_PConvPyIntToInt_From_List(G, list, 0, &I->NIndex));
  ok_assert(1, CPythonVal_PConvPyListToFloatVLANoneOkay_From_List(G, list, 1, &I->Coord));

  ok_assert(2, ll > 2);

  I->LabCoord = NULL; // will be calculated in RepDistLabelNew
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
  ok_assert(1, CPythonVal_PConvPyListToLabPosVLA(G, val, &I->LabPos));
  CPythonVal_Free(val);

  ok_assert(2, ll > 9);

  val = CPythonVal_PyList_GetItem(G, list, 9);
  I->MeasureInfo = MeasureInfoListFromPyList(G, val);
  CPythonVal_Free(val);

ok_except2:
  *cs = I;
  return true;
ok_except1:
  I->fFree();
  return false;
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
    if(I->LabPos) {
      PyList_SetItem(result, 8, PConvLabPosVLAToPyList(I->LabPos, VLAGetSize(I->LabPos)));
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
  float *v;
  int a;
  int c;
  v = I->Coord;
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
/*
 * Invalidate reps
 *
 * type: rep enum, e.g. cRepDash
 * level: e.g. cRepInvColor
 */
void DistSet::invalidateRep(int type, int level)
{
  int a = 0, a_stop = NRep;
  bool changed = false;

  /* if representation type is specified, adjust it */
  if(type >= 0) {
    if(type >= NRep)
      return;

    a = type;
    a_stop = a + 1;
  }

  for(; a < a_stop; a++) {
    if(Rep[a]) {
      changed = true;
      Rep[a]->fFree(Rep[a]);
      Rep[a] = NULL;
    }
  }

  if (changed)
    SceneChanged(State.G);
}


/*========================================================================*/
void DistSet::update(int state)
{
  DistSet * I = this;
  /* status bar 0% */
  OrthoBusyFast(I->State.G, 0, I->NRep);
  if(!I->Rep[cRepDash]) {
    /* query the dist set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepDash] = RepDistDashNew(I,state);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepLabel]) {
    /* query the dist set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepLabel] = RepDistLabelNew(I, state);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepAngle]) {
    /* query the angle set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepAngle] = RepAngleNew(I, state);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepDihedral]) {
    /* query the dihedral set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepDihedral] = RepDihedralNew(I, state);
    SceneInvalidate(I->State.G);
  }
  /* status bar 100% */
  OrthoBusyFast(I->State.G, 1, 1);
}


/*========================================================================*/
void DistSet::render(RenderInfo * info)
{
  DistSet * I = this;
  CRay *ray = info->ray;
  auto pick = info->pick;
  int a;
  ::Rep *r;
  for(a = 0; a < I->NRep; a++)
  {
    if(!GET_BIT(I->Obj->Obj.visRep, a))
      continue;
    if(!I->Rep[a]) {
      switch(a) {
      case cRepDash:
        I->Rep[a] = RepDistDashNew(I, -1);
        break;
      case cRepLabel:
        I->Rep[a] = RepDistLabelNew(I, -1);
        break;
      case cRepAngle:
        I->Rep[a] = RepAngleNew(I, -1);
        break;
      case cRepDihedral:
        I->Rep[a] = RepDihedralNew(I, -1);
        break;
      }
    }
    if(I->Rep[a])
      {
        r = I->Rep[a];
        if(ray || pick) {
          if(ray)
            ray->color3fv(ColorGet(I->State.G, I->Obj->Obj.Color));
          r->fRender(r, info);
        } else {
          ObjectUseColor((CObject *) I->Obj);
          // now all filtering for passes happen in ObjectDistRender
              r->fRender(r, info);
        }
      }
  }
}


/*========================================================================*/
DistSet *DistSetNew(PyMOLGlobals * G)
{
  int a;
  OOAlloc(G, DistSet);
  ObjectStateInit(G, &I->State);
  I->NIndex = 0;
  I->Coord = NULL;
  I->Rep = VLAlloc(Rep *, cRepCnt);
  I->NRep = cRepCnt;
  I->LabPos = NULL;
  I->LabCoord = NULL;
  I->AngleCoord = NULL;
  I->NAngleIndex = 0;
  I->DihedralCoord = NULL;
  I->NDihedralIndex = 0;
  I->NLabel = 0;
  for(a = 0; a < I->NRep; a++)
    I->Rep[a] = NULL;
  I->MeasureInfo = NULL;
  return (I);
}


/*========================================================================*/
#if 0
static void DistSetStrip(DistSet * I)
{
  int a;
  for(a = 0; a < I->NRep; a++)
    if(I->Rep[a])
      I->Rep[a]->fFree(I->Rep[a]);
  I->NRep = 0;
}
#endif

void DistSet::fFree()
{
  DistSet * I = this;
  int a;
  CMeasureInfo * ptr, *target;
  if(I) {
    for(a = 0; a < I->NRep; a++)
      if(I->Rep[a])
	I->Rep[a]->fFree(I->Rep[a]);
    VLAFreeP(I->AngleCoord);
    VLAFreeP(I->DihedralCoord);
    VLAFreeP(I->LabCoord);
    VLAFreeP(I->LabPos);
    VLAFreeP(I->Coord);
    VLAFreeP(I->Rep);

    ptr = I->MeasureInfo;
    while((target = ptr)) {
      ptr = target->next;
      ListElemFree(target);
    }

      /* need to find and decrement the number of dist sets on the objects */
    OOFreeP(I);
  }
}
