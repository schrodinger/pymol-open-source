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
#ifndef _H_DistSet
#define _H_DistSet

#include"Base.h"
#include"Rep.h"
#include"Setting.h"
#include"PyMOLObject.h"

typedef struct DistSet {
  void (*fUpdate)(struct DistSet *I,int state);
  void (*fRender)(struct DistSet *I,RenderInfo *);
  void (*fFree)(struct DistSet *I);
  void (*fInvalidateRep)(struct DistSet *I,int type,int level);
  CObjectState State;
  struct ObjectDist *Obj;
  float *Coord;
  int NIndex;
  Rep **Rep; /* an array of pointers to representations */
  int NRep;
  CSetting *Setting;
  /* extended for mobile distance labels */
  float *LabCoord;
  LabPosType *LabPos;
  int NLabel;
  /* extended for angles and torsions, with labels embedded with coordinates */
  float *AngleCoord;
  int NAngleIndex;
  float *DihedralCoord;
  int NDihedralIndex;
} DistSet;

#include"ObjectDist.h"

DistSet *DistSetNew(PyMOLGlobals *G);
PyObject *DistSetAsPyList(DistSet *I);
int DistSetFromPyList(PyMOLGlobals *G,PyObject *list,DistSet **cs);
int DistSetGetExtent(DistSet *I,float *mn,float *mx);
int DistSetMoveLabel(DistSet *I,int at,float *v,int mode);
int DistSetGetLabelVertex(DistSet *I,int at, float *v);

#endif

