
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
#ifndef _H_RepMesh
#define _H_RepMesh

#include"Rep.h"
#include"CoordSet.h"
#include"CGO.h"

Rep *RepMeshNew(CoordSet * cset, int state);

#define cRepMesh_by_flags     0
#define cRepMesh_all          1
#define cRepMesh_heavy_atoms  2

#endif
