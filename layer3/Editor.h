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
#define cEditorLink   "pklink"

void EditorInit(void);
int EditorActive(void); 
void EditorRender(int state);

void EditorFree(void);
void EditorPrepareDrag(ObjectMolecule *obj,int index,int state);
void EditorDrag(ObjectMolecule *obj,int index,int mode,int state,
                float *pt,float *mov,float *z_dir);

void EditorSetActiveObject(ObjectMolecule *obj,int state,int enable_bond);
ObjectMolecule *EditorDragObject(void);
void EditorReplace(char *elem,int geom,int valence,char *name);
void EditorAttach(char *elem,int geom,int valence,char *name);
void EditorRemove(int hydrogen);
void EditorHFill(void);
void EditorCycleValence(void);
void EditorInactivate(void);
void EditorUpdateState(void);
ObjectMolecule *EditorGetActiveObject(void);
int EditorSelect(char *s0,char *s1,char *s2,char *s3,int pkresi,int pkbond,int quiet);
int EditorTorsion(float angle);
int EditorInvert(ObjectMolecule *obj,int isele0,int isele1,int mode);

PyObject *EditorAsPyList(void);
int EditorFromPyList(PyObject *list);
void EditorGetNextMultiatom(char *name);
int EditorIsObjectNotCurrent(ObjectMolecule *obj);
int EditorGetSinglePicked(char *name);
int EditorIsBondMode(void);
int EditorDeselectIfSelected(int index,int update);
void EditorDefineExtraPks(void);
int EditorGetNFrag(void);

#endif
