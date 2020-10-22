

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
#ifndef _H_Tetsurf
#define _H_Tetsurf

#include"PyMOLGlobals.h"
#include"PyMOLEnums.h"
#include"vla.h"

struct CCrystal;
struct Isofield;
class CarveHelper;

int TetsurfVolume(PyMOLGlobals* G, Isofield* field, float level,
    pymol::vla<int>& num,    //
    pymol::vla<float>& vert, //
    const int* range,        //
    cIsosurfaceMode mode,    //
    const CarveHelper*,      //
    cIsosurfaceSide side);

void TetsurfGetRange(PyMOLGlobals* G, const Isofield* field,
    const CCrystal* cryst, const float* mn, const float* mx, int* range);

int TetsurfInit(PyMOLGlobals * G);
void TetsurfFree(PyMOLGlobals * G);

#endif
