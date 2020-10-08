

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

#include"Map.h"
#include"MemoryDebug.h"
#include"Crystal.h"
#include"Field.h"
#include"PyMOLGlobals.h"
#include"PyMOLEnums.h"

int TetsurfVolume(PyMOLGlobals * G, Isofield * field, float level, int **num,
                  float **vert, int *range, cIsosurfaceMode,
                  MapType * voxelmap, const float* a_vert, float carvebuffer, cIsosurfaceSide);
void TetsurfGetRange(PyMOLGlobals * G, Isofield * field, CCrystal * cryst, float *mn,
                     float *mx, int *range);

int TetsurfInit(PyMOLGlobals * G);
void TetsurfFree(PyMOLGlobals * G);

#endif
