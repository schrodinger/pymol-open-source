
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

#ifndef _H_SceneRay
#define _H_SceneRay

#include"Base.h"
#include"PyMOLObject.h"
#include"Ortho.h"
#include"View.h"

bool SceneRay(PyMOLGlobals * G,
              int ray_width, int ray_height, int mode,
              char **headerVLA_ptr,
              char **charVLA_ptr, float angle,
              float shift, int quiet, G3dPrimitive ** g3d,
              int show_timing, int antialias);

void SceneRenderRayVolume(PyMOLGlobals * G, CScene *I);

#endif
