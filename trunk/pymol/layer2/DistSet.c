
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
#include"Scene.h"
#include"DistSet.h"
#include"Color.h"
#include"RepDistDash.h"
#include"RepDistLabel.h"
#include"RepAngle.h"
#include"RepDihedral.h"
#include"PConv.h"
#include"ObjectMolecule.h"
#include"ListMacros.h"
#include"Selector.h"

static void DistSetUpdate(DistSet * I, int state);
static void DistSetFree(DistSet * I);
static void DistSetInvalidateRep(DistSet * I, int type, int level);

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
        float *lab_pos = SettingGet_3fv(obj->Obj.G, I->Setting, obj->Obj.Setting,
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


/* -- JV */
int DistSetMoveWithObject(DistSet * I, struct ObjectMolecule * O)
{
  /* TODO:
   *  3.  Make this multi-threaded
   */
  PyMOLGlobals * G = I->State.G;

  int a, idx, rVal=0;
  float *src, *dst;
  CMeasureInfo * memb = NULL;
  CoordSet * cs;

  PRINTFD(G, FB_DistSet)
    " DistSet: adjusting distance vertex\n" ENDFD;

  /* If the user doesn't provide the Object O, we should consider
   * looping over all atoms in Selector->Table instead of bailing */
  if (!I || !O) return 0;

  DListIterate(I->MeasureInfo, memb, next)
    {
      /* if this distance set belongs to this object */
      if (memb && memb->obj==O) {
	for (a=0; a<O->NAtom; a++ ) {
	  /* matched on OBJ above and now ATOM ID within the object */
	  if (O->AtomInfo[a].id == memb->id) {

	    /* if the state we found in our record could be right */
	    if (memb->state < O->NCSet) {
	      /* get the coordinate set for this state */
	      cs = O->CSet[memb->state];
	      
	      /* get proper index based on state and discrete flags */
	      if (O->DiscreteFlag) {
		if (cs==O->DiscreteCSet[a]) {
		  idx = O->DiscreteAtmToIdx[a];
		} else {
		  idx = -1;
		}
	      } else {
		idx = cs->AtmToIdx[a];
	      }

	      /* given a valid index, update the coordinate, finally */
	      if (idx>=0) {
		/* valid offset */
		float * varDst = NULL;
		short update = 0;
		if (memb->measureType==cRepDash && memb->offset < I->NIndex) {
		  varDst = I->Coord;
		  I->fInvalidateRep(I, cRepDash, cRepInvCoord);
		  update = 1;
		} else if (memb->measureType==cRepAngle && memb->offset < I->NAngleIndex) {
		  varDst = I->AngleCoord;
		  I->fInvalidateRep(I, cRepAngle, cRepInvCoord);
		  update = 1;
		} else if (memb->measureType==cRepDihedral && memb->offset < I->NDihedralIndex) {
		  varDst = I->DihedralCoord;
		  I->fInvalidateRep(I, cRepDihedral, cRepInvCoord);
		  update = 1;
		}
		/* update the coordinates */
		if (varDst) {
		  src = cs->Coord + 3 * idx;
		  dst = varDst + 3 * memb->offset;
		  copy3f(src,dst);
		  
		  /* no matter what, update labels */
		  I->fInvalidateRep(I, cRepLabel, cRepInvCoord);
		  rVal |= 1;
		}
		if (update){ 
		  /* We need to update these representations right away since this could be 
		     called during the scene rendering loop, and update might not be called
		     on this object again.  This is ok since update checks to see if the reps
		     exist, and if they already do, then the overhead is minimal, so multiple
		     calls to update won't normally happen, but when they do, it won't hurt. */
		  I->fUpdate(I, -1);
		}
	      }
	    }
	  }
	}
      }
    }
  
  PRINTFD(G, FB_DistSet)
    " DistSet: done updating distance set's vertex\n" ENDFD;

  return rVal;
}

int DistSetFromPyList(PyMOLGlobals * G, PyObject * list, DistSet ** cs)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  DistSet *I = NULL;
  int ok = true;
  int ll = 0;
  if(*cs) {
    DistSetFree(*cs);
    *cs = NULL;
  }

  if(list == Py_None) {         /* allow None for CSet */
    *cs = NULL;
  } else {

    if(ok)
      I = DistSetNew(G);
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
      ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list, 1), &I->Coord);

    if(ok && (ll > 2)) {        /* have angles and torsions too... */
      if(ok)
        ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list, 2), &I->LabCoord);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 3), &I->NAngleIndex);
      if(ok)
        ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list, 4), &I->AngleCoord);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 5), &I->NDihedralIndex);
      if(ok)
        ok = PConvPyListToFloatVLANoneOkay(PyList_GetItem(list, 6), &I->DihedralCoord);
    }
    if(ok && (ll > 7))
      I->Setting = SettingNewFromPyList(G, PyList_GetItem(list, 7));    /* state settings */
    if(ok && (ll > 8))
      ok = PConvPyListToLabPosVLA(PyList_GetItem(list, 8), &I->LabPos);

    if(!ok) {
      if(I)
        DistSetFree(I);
    } else {
      *cs = I;
    }
  }
  return (ok);
#endif
}

PyObject *DistSetAsPyList(DistSet * I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;

  if(I) {
    result = PyList_New(9);

    PyList_SetItem(result, 0, PyInt_FromLong(I->NIndex));
    PyList_SetItem(result, 1, PConvFloatArrayToPyListNullOkay(I->Coord, I->NIndex * 3));
    PyList_SetItem(result, 2, PConvFloatArrayToPyListNullOkay(I->LabCoord, I->NLabel));
    PyList_SetItem(result, 3, PyInt_FromLong(I->NAngleIndex));
    PyList_SetItem(result, 4,
                   PConvFloatArrayToPyListNullOkay(I->AngleCoord, I->NAngleIndex * 3));
    PyList_SetItem(result, 5, PyInt_FromLong(I->NDihedralIndex));
    PyList_SetItem(result, 6,
                   PConvFloatArrayToPyListNullOkay(I->DihedralCoord,
                                                   I->NDihedralIndex * 3));
    PyList_SetItem(result, 7, SettingAsPyList(I->Setting));
    if(I->LabPos) {
      PyList_SetItem(result, 8, PConvLabPosVLAToPyList(I->LabPos, VLAGetSize(I->LabPos)));
    } else {
      PyList_SetItem(result, 8, PConvAutoNone(NULL));
    }
    /* TODO setting ... */
  }
  return (PConvAutoNone(result));
#endif
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
static void DistSetInvalidateRep(DistSet * I, int type, int level)
{
  int a;
  PRINTFD(I->State.G, FB_DistSet)
    " DistSetInvalidateRep: entered.\n" ENDFD;
  /* if representation type is specified, adjust it */
  if(type >= 0) {
    if(type < I->NRep) {
      SceneChanged(I->State.G);
      if(I->Rep[type]) {
        I->Rep[type]->fFree(I->Rep[type]);
        I->Rep[type] = NULL;
      }
    }
  } else {
    /* reset all representation types */
    for(a = 0; a < I->NRep; a++) {
      SceneChanged(I->State.G);
      if(I->Rep[a]) {
        switch (level) {
        case cRepInvColor:
          if(I->Rep[a]->fRecolor) {
            I->Rep[a]->fInvalidate(I->Rep[a], (struct CoordSet *) I, level);
          } else {
            I->Rep[a]->fFree(I->Rep[a]);
            I->Rep[a] = NULL;
          }
          break;
        default:
          I->Rep[a]->fFree(I->Rep[a]);
          I->Rep[a] = NULL;
          break;
        }
      }
    }
  }
}


/*========================================================================*/
static void DistSetUpdate(DistSet * I, int state)
{
  /* status bar 0% */
  OrthoBusyFast(I->State.G, 0, I->NRep);
  if(!I->Rep[cRepDash]) {
    /* query the dist set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepDash] = RepDistDashNew(I);
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
    I->Rep[cRepAngle] = RepAngleNew(I);
    SceneInvalidate(I->State.G);
  }
  if(!I->Rep[cRepDihedral]) {
    /* query the dihedral set looking for the selected atoms for this distance,
     * then update the *coords */
    I->Rep[cRepDihedral] = RepDihedralNew(I);
    SceneInvalidate(I->State.G);
  }
  /* status bar 100% */
  OrthoBusyFast(I->State.G, 1, 1);
}


/*========================================================================*/
static void DistSetRender(DistSet * I, RenderInfo * info)
{
  CRay *ray = info->ray;
  int pass = info->pass;
  Picking **pick = info->pick;
  int float_labels = SettingGet_i(I->State.G,
                                  I->Setting,
                                  I->Obj->Obj.Setting,
                                  cSetting_float_labels);
  int a;
  Rep *r;
  for(a = 0; a < I->NRep; a++)
    if(I->Rep[a])
      if(I->Obj->Obj.RepVis[a]) {
        r = I->Rep[a];
        if(ray || pick) {
          if(ray)
            ray->fColor3fv(ray, ColorGet(I->State.G, I->Obj->Obj.Color));
          r->fRender(r, info);
        } else {
          ObjectUseColor((CObject *) I->Obj);
          switch (a) {
          case cRepLabel:
            if(float_labels) {
              if(pass == -1)
                r->fRender(r, info);
            } else if(pass == 0)
              r->fRender(r, info);
            break;
          default:
            if(pass == 0) {
              r->fRender(r, info);
            }
            break;
          }
        }
      }
}


/*========================================================================*/
DistSet *DistSetNew(PyMOLGlobals * G)
{
  int a;
  OOAlloc(G, DistSet);
  I->State.G = G;
  I->State.Matrix = NULL;
  I->fFree = DistSetFree;
  I->fRender = DistSetRender;
  I->fUpdate = DistSetUpdate;
  I->fInvalidateRep = DistSetInvalidateRep;
  I->NIndex = 0;
  I->Coord = NULL;
  I->Rep = VLAlloc(Rep *, cRepCnt);
  I->NRep = cRepCnt;
  I->Setting = NULL;
  I->LabPos = NULL;
  I->LabCoord = NULL;
  I->AngleCoord = NULL;
  I->NAngleIndex = 0;
  I->DihedralCoord = NULL;
  I->NDihedralIndex = 0;
  I->NLabel = 0;
  for(a = 0; a < I->NRep; a++)
    I->Rep[a] = NULL;
  DListInit(I->MeasureInfo, prev, next, CMeasureInfo);
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

/*========================================================================*/
static void DistSetFree(DistSet * I)
{
  int a;
  CMeasureInfo * ptr, *target;
  for(a = 0; a < I->NRep; a++)
    if(I->Rep[a])
      I->Rep[a]->fFree(I->Rep[a]);
  if(I) {
    VLAFreeP(I->AngleCoord);
    VLAFreeP(I->DihedralCoord);
    VLAFreeP(I->LabCoord);
    VLAFreeP(I->LabPos);
    VLAFreeP(I->Coord);
    VLAFreeP(I->Rep);
    if (I->MeasureInfo) {
      ptr = I->MeasureInfo->next;
      while(ptr != I->MeasureInfo) {
	target = ptr;
	ptr = ptr->next;
	DListRemove(target,prev,next);
	DListElemFree(target);
	target = NULL;
      }	
    }
      /* need to find and decrement the number of dist sets on the objects */
    SettingFreeP(I->Setting);
    OOFreeP(I);
  }
}
