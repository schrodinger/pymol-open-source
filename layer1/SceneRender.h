
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

#ifndef _H_SceneRender
#define _H_SceneRender

void SceneRender(PyMOLGlobals * G, Picking * pick, int x, int y,
                 Multipick * smp, int oversize_width, int oversize_height,
                 int click_side, int force_copy);
void SceneRenderAll(PyMOLGlobals * G, SceneUnitContext * context,
                    float *normal, Picking ** pickVLA,
                    int pass, int fat, float width_scale,
                    GridInfo * grid, int dynamic_pass, short which_objects, bool picking32bit);

void SceneInitializeViewport(PyMOLGlobals * G, int offscreen);

void GridGetGLViewport(PyMOLGlobals * G, GridInfo * I);
void GridSetGLViewport(GridInfo * I, int slot);

#endif

