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
#include"os_std.h"
#include"os_gl.h"

#include"Base.h"
#include"Debug.h"
#include"Parse.h"
#include"OOMac.h"
#include"Vector.h"
#include"PConv.h"
#include"ObjectMolecule.h"

/*========================================================================*/
static PyObject *ObjectMoleculeGetCSetPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  int a;
  result = PyList_New(I->NCSet);
  for(a=0;a<I->NCSet;a++) {
    if(I->CSet[a]) {
      PyList_SetItem(result,a,CoordSetGetPyList(I->CSet[a]));
    } else {
      PyList_SetItem(result,a,Py_None);
      Py_INCREF(Py_None);
    }
  }
  return(PConvAutoNone(result));
}

static PyObject *ObjectMoleculeGetDiscreteCSetPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  return(PConvAutoNone(result));
}

static int ObjectMoleculeSetCSetPyList(ObjectMolecule *I,PyObject *list)
{
  int ok=true;
  int a;
  if(ok) ok=PyList_Check(list);
  if(ok) {
    VLACheck(I->CSet,CoordSet*,I->NCSet);
    for(a=0;a<I->NCSet;a++) {
      if(ok) ok = CoordSetSetPyList(PyList_GetItem(list,a),&I->CSet[a]);
      if(ok) I->CSet[a]->Obj = I;
    }
  }
  return(ok);
}

static PyObject *ObjectMoleculeGetBondPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  PyObject *bond_list;
  BondType *bond;
  int a;

  result = PyList_New(I->NBond);  
  bond = I->Bond;
  for(a=0;a<I->NBond;a++) {
    bond_list=PyList_New(5);
    PyList_SetItem(bond_list,0,PyInt_FromLong(bond->index[0]));
    PyList_SetItem(bond_list,1,PyInt_FromLong(bond->index[1]));
    PyList_SetItem(bond_list,2,PyInt_FromLong(bond->order));
    PyList_SetItem(bond_list,3,PyInt_FromLong(bond->id));
    PyList_SetItem(bond_list,4,PyInt_FromLong(bond->stereo));
    PyList_SetItem(result,a,bond_list);
    bond++;
  }

  return(PConvAutoNone(result));
}

static int ObjectMoleculeSetBondPyList(ObjectMolecule *I,PyObject *list) 
{
  int ok=true;
  int a;
  PyObject *bond_list;
  BondType *bond;
  if(ok) ok=PyList_Check(list);  
  if(ok) ok=((I->Bond=VLAlloc(BondType,I->NBond))!=NULL);
  bond = I->Bond;
  for(a=0;a<I->NBond;a++) {
    if(ok) bond_list = PyList_GetItem(list,a);
    if(ok) ok = PyList_Check(bond_list);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,0),&bond->index[0]);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,1),&bond->index[1]);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,2),&bond->order);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,3),&bond->id);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(bond_list,4),&bond->stereo);
    bond++;
  }
  return(ok);
}

static PyObject *ObjectMoleculeGetAtomPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;
  AtomInfoType *ai;
  int a;

  result = PyList_New(I->NAtom);  
  ai = I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    PyList_SetItem(result,a,AtomInfoGetPyList(ai));
    ai++;
  }
  return(PConvAutoNone(result));
}

static int ObjectMoleculeSetAtomPyList(ObjectMolecule *I,PyObject *list) 
{
  int ok=true;
  int a;
  AtomInfoType *ai;
  if(ok) ok=PyList_Check(list);  
  VLACheck(I->AtomInfo,AtomInfoType,I->NAtom+1);
  ai = I->AtomInfo;
  for(a=0;a<I->NAtom;a++) {
    if(ok) ok = AtomInfoSetPyList(ai,PyList_GetItem(list,a));
    ai++;
  }
  return(ok);
}

int ObjectMoleculeNewFromPyList(PyObject *list,ObjectMolecule **result)
{
  int ok = true;
  ObjectMolecule *I=NULL;
  int discrete_flag;
  (*result) = NULL;
  

  if(ok) ok=PyList_Check(list);

  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,8),&discrete_flag);

  I=ObjectMoleculeNew(discrete_flag);
  if(ok) ok = (I!=NULL);

  if(ok) ok = ObjectSetPyList(PyList_GetItem(list,0),&I->Obj);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&I->NCSet);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->NBond);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,3),&I->NAtom);
  if(ok) ok = ObjectMoleculeSetCSetPyList(I,PyList_GetItem(list,4));
  if(ok) ok = CoordSetSetPyList(PyList_GetItem(list,5),&I->CSTmpl);
  if(ok) ok = ObjectMoleculeSetBondPyList(I,PyList_GetItem(list,6));
  if(ok) ok = ObjectMoleculeSetAtomPyList(I,PyList_GetItem(list,7));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,8),&I->DiscreteFlag);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,9),&I->NDiscrete);
  if(ok) I->Symmetry = SymmetryNewFromPyList(PyList_GetItem(list,10));
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,11),&I->CurCSet);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,12),&I->BondCounter);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,13),&I->AtomCounter);
  if(ok&&I->DiscreteFlag) {
    int *dcs = NULL;
    int a,i;
    CoordSet *cs;
    VLACheck(I->DiscreteAtmToIdx,int,I->NDiscrete);
    VLACheck(I->DiscreteCSet,CoordSet*,I->NDiscrete);
    if(ok) ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list,14),I->DiscreteAtmToIdx,I->NDiscrete);
    if(ok) ok = PConvPyListToIntArray(PyList_GetItem(list,15),&dcs);
    if(ok) {

      for(a=0;a<I->NDiscrete;a++) {
        i = dcs[a];
        I->DiscreteCSet[a] = NULL;
        if(i>=0) {
          if(i<I->NCSet) {
            cs = I->CSet[i];
            if(cs) {
              I->DiscreteCSet[a]=cs;
            }
          }
        }
      }
    }
    FreeP(dcs);
  }
  
  ObjectMoleculeInvalidate(I,cRepAll,cRepInvAll);
  if(ok) 
    (*result) = I;
  else {
    /* cleanup? */
  }

  return(ok);
}


/*========================================================================*/
PyObject *ObjectMoleculeGetPyList(ObjectMolecule *I)
{
  PyObject *result = NULL;


  /* first, dump the atoms */

  result = PyList_New(16);
  PyList_SetItem(result,0,ObjectGetPyList(&I->Obj));
  PyList_SetItem(result,1,PyInt_FromLong(I->NCSet));
  PyList_SetItem(result,2,PyInt_FromLong(I->NBond));
  PyList_SetItem(result,3,PyInt_FromLong(I->NAtom));
  PyList_SetItem(result,4,ObjectMoleculeGetCSetPyList(I));
  PyList_SetItem(result,5,CoordSetGetPyList(I->CSTmpl));
  PyList_SetItem(result,6,ObjectMoleculeGetBondPyList(I));
  PyList_SetItem(result,7,ObjectMoleculeGetAtomPyList(I));
  PyList_SetItem(result,8,PyInt_FromLong(I->DiscreteFlag));
  PyList_SetItem(result,9,PyInt_FromLong(I->NDiscrete));
  PyList_SetItem(result,10,SymmetryGetPyList(I->Symmetry));
  PyList_SetItem(result,11,PyInt_FromLong(I->CurCSet));
  PyList_SetItem(result,12,PyInt_FromLong(I->BondCounter));
  PyList_SetItem(result,13,PyInt_FromLong(I->AtomCounter));


  if(I->DiscreteFlag) {
    int *dcs;
    int a;
    CoordSet *cs;
    
    /* get coordinate set indices */
    
    for(a=0;a<I->NCSet;a++) {
      cs = I->CSet[a];
      if(cs) {
        cs->tmp_index = a;
      }
    }
    
    dcs = Alloc(int,I->NDiscrete);
    
    for(a=0;a<I->NDiscrete;a++) {
      cs = I->DiscreteCSet[a];
      if(cs) 
        dcs[a] = cs->tmp_index;
      else
        dcs[a] = -1;
    }

    PyList_SetItem(result,14,PConvIntArrayToPyList(I->DiscreteAtmToIdx,I->NDiscrete));
    PyList_SetItem(result,15,PConvIntArrayToPyList(dcs,I->NDiscrete));
    FreeP(dcs);
  } else {
    PyList_SetItem(result,14,PConvAutoNone(NULL));
    PyList_SetItem(result,15,PConvAutoNone(NULL));
  }

#if 0
  CObject Obj;
  struct CoordSet **CSet;
  int NCSet;
  struct CoordSet *CSTmpl; /* template for trajectories, etc.*/
  BondType *Bond;
  AtomInfoType *AtomInfo;
  int NAtom;
  int NBond;
  int DiscreteFlag,NDiscrete;
  int *DiscreteAtmToIdx;
  struct CoordSet **DiscreteCSet;
  int CurCSet;
  int SeleBase; /* for internal usage by  selector & only valid during selection process */
  CSymmetry *Symmetry;
  int *Neighbor;
  float *UndoCoord[cUndoMask+1];
  int UndoState[cUndoMask+1];
  int UndoNIndex[cUndoMask+1];
  int UndoIter;
  CGO *UnitCellCGO;
  int BondCounter;
  int AtomCounter;
  struct CSculpt *Sculpt;
#endif

  return(PConvAutoNone(result));  
}


