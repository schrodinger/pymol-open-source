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
#ifndef _H_Editor
#define _H_Editor

#include"Ortho.h"
#include"ButMode.h"
#include"ObjectMolecule.h"

#define cEditorSele1 "ed1"
#define cEditorSele2 "ed2"
#define cEditorFragPref "fg"
#define cEditorBasePref "_fbase"
#define cEditorComp   "ed"

void EditorInit(void);
int EditorActive(void);
void EditorRender(int state);

void EditorFree(void);
void EditorPrepareDrag(ObjectMolecule *obj,int index,int state);
void EditorDrag(ObjectMolecule *obj,int index,int mode,int state,float *pt,float *mov);

void EditorSetActiveObject(ObjectMolecule *obj,int state);
ObjectMolecule *EditorDragObject(void);
void EditorReplace(char *elem,int geom,int valence);
void EditorAttach(char *elem,int geom,int valence);
void EditorRemove(void);
void EditorRefill(void);
void EditorCycleValence(void);
void EditorInactive(void);
ObjectMolecule *EditorGetActiveObject(void);

#endif
