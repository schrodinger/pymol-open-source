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
#ifndef _H_ObjectDist
#define _H_ObjectDist

#include"PyMOLObject.h"
#include"AtomInfo.h"
#include"Vector.h"
#include"Color.h"

typedef struct ObjectDist {
  CObject Obj;
  struct DistSet **DSet;
  int NDSet;
  int CurDSet;
} ObjectDist;

ObjectDist *ObjectDistNewFromSele(PyMOLGlobals *G,ObjectDist *oldObj,int sele1,int sele2,int mode,float cutoff,
                                  int labels,float *result);
ObjectDist *ObjectDistNew(PyMOLGlobals *G);
void ObjectDistInvalidateRep(ObjectDist *I,int rep);
PyObject *ObjectDistAsPyList(ObjectDist *I);
int ObjectDistNewFromPyList(PyMOLGlobals *G,PyObject *list,ObjectDist **result);

struct M4XBondType;
struct ObjectMolecule;

ObjectDist *ObjectDistNewFromM4XBond(PyMOLGlobals *G,ObjectDist *oldObj,                                      
                                      struct ObjectMolecule *objMol,
                                     struct M4XBondType *hbond,int n_hbond,
                                     int nbr_sele);

#endif











