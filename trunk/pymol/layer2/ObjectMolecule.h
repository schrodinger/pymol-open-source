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
#ifndef _H_ObjectMolecule
#define _H_ObjectMolecule

#include<Python.h>

#include"Object.h"
#include"AtomInfo.h"
#include"Vector.h"
#include"Color.h"
#include"Symmetry.h"

typedef struct ObjectMolecule {
  Object Obj;
  struct CoordSet **CSet;
  int NCSet;
  int *Bond;
  AtomInfoType *AtomInfo;
  int NAtom;
  int NBond;
  int CurCSet;
  char Color[3];
  float FractionExposed;
  int SeleBase;
  CSymmetry *Symmetry;
} ObjectMolecule;

typedef struct ObjectMoleculeOpRec {
  unsigned long code;
  Vector3f v1,v2;
  int cs1;
  int i1,i2,i3,*vc1;
  float f1,*f1VLA;
  double d[3][3];
  float *vv1,*vv2;
  char *charVLA;
  char *s1;
  float ttt[16];
  int nvv1,nvv2;
} ObjectMoleculeOpRec;

/* these four letter code are left over from an 
   earlier multicharacter constant implementation
   and should be replaced with something more verbose */

#define OMOP_PDB1 1
#define OMOP_AVRT 2
#define OMOP_SFIT 3
#define OMOP_COLR 4
#define OMOP_VISI 5
#define OMOP_TTTF 6
#define OMOP_ALTR 7
#define OMOP_CSOC 8
#define OMOP_SUMC 9
#define OMOP_VERT 10
#define OMOP_SVRT 11
#define OMOP_MOME 12
#define OMOP_INVA 13
#define OMOP_MDST 14
#define OMOP_MNMX 15
#define OMOP_AlterState 16
#define OMOP_Flag 17

#include"CoordSet.h"

ObjectMolecule *ObjectMoleculeNew(void);
void ObjectMoleculeSort(ObjectMolecule *I);
ObjectMolecule *ObjectMoleculeCopy(ObjectMolecule *obj);

ObjectMolecule *ObjectMoleculeLoadPDBFile(ObjectMolecule *obj,char *fname,int frame);
ObjectMolecule *ObjectMoleculeLoadMOLFile(ObjectMolecule *obj,char *fname,int frame);
ObjectMolecule *ObjectMoleculeLoadMMDFile(ObjectMolecule *obj,char *fname,
                                          int frame,char *sepPrefix);
ObjectMolecule *ObjectMoleculeLoadChempyModel(ObjectMolecule *I,PyObject *model,int frame);

ObjectMolecule *ObjectMoleculeReadMOLStr(ObjectMolecule *obj,char *molstr,int frame);
ObjectMolecule *ObjectMoleculeReadPDBStr(ObjectMolecule *obj,char *molstr,int frame);
ObjectMolecule *ObjectMoleculeReadMMDStr(ObjectMolecule *I,char *MMDStr,int frame);
void ObjectMoleculeExtendIndices(ObjectMolecule *I);

void ObjectMoleculeInvalidateRep(ObjectMolecule *I,int rep);

void ObjectMoleculeSeleOp(ObjectMolecule *I,int sele,ObjectMoleculeOpRec *op);

CoordSet *ObjectMoleculeGetCoordSet(ObjectMolecule *I,int setIndex);
void ObjectMoleculeBlindSymMovie(ObjectMolecule *I);
void ObjectMoleculeMerge(ObjectMolecule *I,AtomInfoType *ai,CoordSet *cs,int bondSearchFlag);

#endif











