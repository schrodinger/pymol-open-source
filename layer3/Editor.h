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
#define cEditorSele3 "pk3"
#define cEditorSele4 "pk4"
#define cEditorFragPref "_pkfrag"
#define cEditorBasePref "_pkbase"
#define cEditorSet    "pkset"
#define cEditorRes    "pkresi"
#define cEditorChain  "pkchain"
#define cEditorObject "pkobject"
#define cEditorComp   "pkmol"
#define cEditorLink   "pkfrag"

void EditorInit(void);
int EditorActive(void); 
void EditorRender(int state);
int EditorLogState(int pkresi);

void EditorFree(void);
void EditorPrepareDrag(ObjectMolecule *obj,int index,int state);
void EditorDrag(ObjectMolecule *obj,int index,int mode,int state,
                float *pt,float *mov,float *z_dir);

void EditorActivate(int state,int enable_bond);
ObjectMolecule *EditorDragObject(void);
void EditorReplace(char *elem,int geom,int valence,char *name,int quiet);
void EditorAttach(char *elem,int geom,int valence,char *name,int quiet);
void EditorRemove(int hydrogen,int quiet);
void EditorHFill(int quiet);
void EditorCycleValence(int quiet);
void EditorInactivate(void);
void EditorUpdateState(void);

int EditorIsAnActiveObject(ObjectMolecule *obj);

int EditorSelect(char *s0,char *s1,char *s2,char *s3,int pkresi,int pkbond,int quiet);
int EditorTorsion(float angle);
int EditorInvert(int quiet);

PyObject *EditorAsPyList(void);
int EditorFromPyList(PyObject *list);
void EditorGetNextMultiatom(char *name);
int EditorGetSinglePicked(char *name);
int EditorIsBondMode(void);
int EditorDeselectIfSelected(ObjectMolecule *obj,int index,int update);
void EditorDefineExtraPks(void);
int EditorGetNFrag(void);

#endif
