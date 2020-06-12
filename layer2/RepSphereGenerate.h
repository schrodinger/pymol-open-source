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
#ifndef _H_RepSphereGenerate
#define _H_RepSphereGenerate

void RepSphere_Generate_Triangles(PyMOLGlobals *G, RepSphere *I,
                                  RenderInfo *info);
void RepSphere_Generate_Impostor_Spheres(PyMOLGlobals *G, RepSphere *I,
                                         RenderInfo *info);
void RepSphere_Generate_Point_Sprites(PyMOLGlobals *G, RepSphere *I,
                                      RenderInfo *info, int sphere_mode);

#endif
