
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

/* NOTE: "Dist" names & symbols should be updated to "Measurement" */

typedef struct ObjectDist {
  /* base class, this IS-A CObject */
  CObject Obj;
  /* Array of pointers to DistSets */
  struct DistSet **DSet;
  /* number of dist sets */
  int NDSet;
  /* current dist set */
  int CurDSet;
  /* ptr to next, if in list */
} ObjectDist;

ObjectDist *ObjectDistNewFromSele(PyMOLGlobals * G, ObjectDist * oldObj,
                                  int sele1, int sele2, int mode, float cutoff,
                                  int labels, int reset, float *result, int state);

ObjectDist *ObjectDistNewFromAngleSele(PyMOLGlobals * G, ObjectDist * oldObj,
                                       int sele1, int sele2, int sele3, int mode,
                                       int labels, float *result, int reset, int state);

ObjectDist *ObjectDistNewFromDihedralSele(PyMOLGlobals * G, ObjectDist * oldObj,
                                          int sele1, int sele2, int sele3, int sele4,
                                          int mode, int labels, float *result,
                                          int reset, int state);

int ObjectDistGetLabelTxfVertex(ObjectDist * I, int state, int index, float *v);
int ObjectDistMoveLabel(ObjectDist * I, int state, int index, float *v, int mode,
                        int log);

ObjectDist *ObjectDistNew(PyMOLGlobals * G);
void ObjectDistInvalidateRep(ObjectDist * I, int rep);
PyObject *ObjectDistAsPyList(ObjectDist * I);
int ObjectDistNewFromPyList(PyMOLGlobals * G, PyObject * list, ObjectDist ** result);

/* -- JV -- dyndist */
/* An ObjectMolecule moved, causing this function to be called */
struct ObjectMolecule;
int ObjectDistMoveWithObject(ObjectDist * I, struct ObjectMolecule * O);

struct M4XBondType;
ObjectDist *ObjectDistNewFromM4XBond(PyMOLGlobals * G, ObjectDist * oldObj,
                                     struct ObjectMolecule *objMol,
                                     struct M4XBondType *hbond, int n_hbond,
                                     int nbr_sele);

#endif
