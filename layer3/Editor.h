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

#define cEditorSele1 "pk1"
#define cEditorSele2 "pk2"
#define cEditorFragPref "pkfrag"
#define cEditorBasePref "_pkbase"
#define cEditorRes    "pkresi"
#define cEditorComp   "pkchain"

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
void EditorRemove(int hydrogen);
void EditorHFill(void);
void EditorCycleValence(void);
void EditorInactive(void);
ObjectMolecule *EditorGetActiveObject(void);
int EditorSelect(char *s0,char *s1,char *s2,char *s3,int pkresi);
int EditorTorsion(float angle);
int EditorInvert(ObjectMolecule *obj,int isele0,int isele1,int mode);

PyObject *EditorAsPyList(void);
int EditorFromPyList(PyObject *list);

#endif
