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
#define cEditorBond   "pkbond"
#define cEditorRes    "pkresi"
#define cEditorChain  "pkchain"
#define cEditorObject "pkobject"
#define cEditorComp   "pkmol"
#define cEditorLink   "pkfrag"
#define cEditorDihedral "_pkdihe"
#define cEditorDihe1    "_pkdihe1"
#define cEditorDihe2    "_pkdihe2"
#define cEditorDrag    "_drag"

#define EDITOR_SCHEME_OBJ 1
#define EDITOR_SCHEME_FRAG 2
#define EDITOR_SCHEME_DRAG 3

int EditorInit(PyMOLGlobals *G);
int EditorActive(PyMOLGlobals *G); 
void EditorRender(PyMOLGlobals *G,int state);
int EditorLogState(PyMOLGlobals *G,int pkresi);
void EditorFavorOrigin(PyMOLGlobals *G, float *v1);

void EditorFree(PyMOLGlobals *G);
void EditorSetDrag(PyMOLGlobals *G,ObjectMolecule *obj,int sele, int quiet,int state);
void EditorReadyDrag(PyMOLGlobals *G,int state);
void EditorPrepareDrag(PyMOLGlobals *G,ObjectMolecule *obj,int sele, int index,int state, int mode);
void EditorDrag(PyMOLGlobals *G,ObjectMolecule *obj,int index,int mode,int state,
                float *pt,float *mov,float *z_dir);

void EditorActivate(PyMOLGlobals *G,int state,int enable_bond);
ObjectMolecule *EditorDragObject(PyMOLGlobals *G);
void EditorReplace(PyMOLGlobals *G,char *elem,int geom,int valence,char *name,int quiet);
void EditorAttach(PyMOLGlobals *G,char *elem,int geom,int valence,char *name,int quiet);
void EditorRemove(PyMOLGlobals *G,int hydrogen,int quiet);
void EditorHFill(PyMOLGlobals *G,int quiet);
void EditorHFix(PyMOLGlobals *G,char *sele,int quiet);
void EditorCycleValence(PyMOLGlobals *G,int quiet);
void EditorInactivate(PyMOLGlobals *G);
void EditorUpdateState(PyMOLGlobals *G);

int EditorIsAnActiveObject(PyMOLGlobals *G,ObjectMolecule *obj);

int EditorSelect(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,int pkresi,int pkbond,int quiet);
int EditorTorsion(PyMOLGlobals *G,float angle);
int EditorInvert(PyMOLGlobals *G,int quiet);

PyObject *EditorAsPyList(PyMOLGlobals *G);
int EditorFromPyList(PyMOLGlobals *G,PyObject *list);
void EditorGetNextMultiatom(PyMOLGlobals *G,char *name);
int EditorGetSinglePicked(PyMOLGlobals *G,char *name);
int EditorIsBondMode(PyMOLGlobals *G);
int EditorDeselectIfSelected(PyMOLGlobals *G,ObjectMolecule *obj,int index,int update);
void EditorDefineExtraPks(PyMOLGlobals *G);
int EditorGetNFrag(PyMOLGlobals *G);
void EditorUpdate(PyMOLGlobals *G);
void EditorMouseInvalid(PyMOLGlobals *G);
int EditorGetScheme(PyMOLGlobals *G);
void EditorDihedralInvalid(PyMOLGlobals *G,ObjectMolecule *obj);

#endif
