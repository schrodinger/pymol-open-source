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

#ifndef _H_ScenePicking
#define _H_ScenePicking

#include"Base.h"
#include"PyMOLObject.h"
#include"Ortho.h"
#include"View.h"

int SceneDoXYPick(PyMOLGlobals * G, int x, int y, int click_side);

void SceneRenderPicking(PyMOLGlobals * G, int stereo_mode, int *click_side, int stereo_double_pump_mono, 
			Picking * pick, int x, int y, Multipick * smp, SceneUnitContext *context,
			GLenum render_buffer);
int SceneMultipick(PyMOLGlobals * G, Multipick * smp);

#endif

