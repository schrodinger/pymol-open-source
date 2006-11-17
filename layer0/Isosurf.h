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
#ifndef _H_Isosurf
#define _H_Isosurf

#include"Map.h"
#include"MemoryDebug.h"
#include"Crystal.h"
#include"Field.h"
#include"os_python.h"
#include"PyMOLGlobals.h"
#include"Setting.h"

typedef struct {
  int dimensions[3];
  int save_points;
  CField *points;
  CField *data;
  CField *gradients;
} Isofield;

#define F3(field,P1,P2,P3) Ffloat3(field,P1,P2,P3)
#define F3Ptr(field,P1,P2,P3) Ffloat3p(field,P1,P2,P3)

#define F4(field,P1,P2,P3,P4) Ffloat4(field,P1,P2,P3,P4)
#define F4Ptr(field,P1,P2,P3,P4) Ffloat4p(field,P1,P2,P3,P4)

Isofield *IsosurfFieldAlloc(PyMOLGlobals *G,int *dims);
void IsosurfFieldFree(PyMOLGlobals *G,Isofield *field);

int	IsosurfVolume(PyMOLGlobals *G,CSetting *set1,CSetting *set2,
                  Isofield *field,float level,int **num,
                  float **vert,int *range,int mode,int skip,float alt_level);

void IsosurfGetRange(PyMOLGlobals *G,Isofield *field,CCrystal *cryst,float *mn,float *mx,int *range);

int	IsosurfInit(PyMOLGlobals *G);
void  IsosurfFree(PyMOLGlobals *G);

/* isofield operations -- not part of Isosurf */

void IsofieldComputeGradients(PyMOLGlobals *G,Isofield *field);
PyObject *IsosurfAsPyList(Isofield *I);
Isofield *IsosurfNewFromPyList(PyMOLGlobals *G,PyObject *list);
Isofield *IsosurfNewCopy(PyMOLGlobals *G,Isofield *src);

#endif
