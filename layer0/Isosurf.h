

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#include"Symmetry.h"
#include"Field.h"
#include"os_python.h"
#include"PyMOLGlobals.h"
#include"PyMOLEnums.h"
#include"Setting.h"

struct Isofield {
  int dimensions[3]{};
  int save_points = true;
  pymol::copyable_ptr<CField> points;
  pymol::copyable_ptr<CField> data;
  pymol::cache_ptr<CField> gradients;
  Isofield() = default;
  Isofield(PyMOLGlobals * G, const int * const dims);
};

int IsosurfVolume(PyMOLGlobals* G, CSetting* set1, CSetting* set2,
    Isofield* field, float level, pymol::vla<int>& num, pymol::vla<float>& vert,
    int* range, cIsomeshMode, int skip, float alt_level);

int IsosurfGetRange(PyMOLGlobals * G, Isofield * field, CCrystal * cryst,
                    float *mn, float *mx, int *range, int clamp);
int IsosurfExpand(Isofield * field1, Isofield * field2,
                  CCrystal * cryst, CSymmetry * sym, int *range);

int IsosurfInit(PyMOLGlobals * G);
void IsosurfFree(PyMOLGlobals * G);


/* isofield operations -- not part of Isosurf */

void IsofieldComputeGradients(PyMOLGlobals * G, Isofield * field);
PyObject *IsosurfAsPyList(PyMOLGlobals *G, Isofield * I);
Isofield *IsosurfNewFromPyList(PyMOLGlobals * G, PyObject * list);

void IsofieldGetCorners(PyMOLGlobals *, Isofield *, float *);

#endif
