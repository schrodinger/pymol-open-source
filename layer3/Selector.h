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
#ifndef _H_Selector
#define _H_Selector

#include"os_python.h"

#include"ObjectMolecule.h"
#include"DistSet.h"
#include"ObjectMap.h"

void SelectorInit(void);
int *SelectorSelect(char *sele);
int SelectorCreate(char *name,char *sele,ObjectMolecule *obj,int quiet,Multipick *mp);
void SelectorToggle(int rep,char *name);
void SelectorCylinder(char *sele,char *onoff);
int SelectorUpdateTable(void);
int SelectorIndexByName(char *sele);
int SelectorIsMember(int start,int sele);
void SelectorFree(void);
void SelectorDelete(char *sele);
void SelectorFreeTmp(char *name);
int SelectorGetTmp(char *input,char *store);
int SelectorGetPDB(char **charVLA,int sele,int state,int conectFlag);
PyObject *SelectorGetChemPyModel(int sele,int state);
float SelectorSumVDWOverlap(int sele1,int state1,int sele2,int state2,float adjust);
DistSet *SelectorGetDistSet(int sele1,int state1,int sele2,int state2,int mode,
                            float cutoff,float *result);
int SelectorGetSeleNCSet(int sele);
void SelectorCreateObjectMolecule(int sele,char *name,int target_state,int state);
int SelectorSubdivideObject(char *pref,ObjectMolecule *obj,int sele1,int sele2,
                            char *fragPref,char *compName);
ObjectMolecule *SelectorGetSingleObjectMolecule(int sele);
void SelectorUpdateObjectSele(ObjectMolecule *obj);
void SelectorDeletePrefixSet(char *pref);
void SelectorUpdateCmd(int sele0,int sele1,int sta0,int sta1);
int SelectorGetSingleAtomVertex(int sele,int state,float *v);
int SelectorGetSingleAtomObjectIndex(int sele,ObjectMolecule **in_obj,int *index);
int *SelectorGetResidueVLA(int sele0);
int  SelectorCreateAlignments(int *pair,int sele1,int *vla1,int sele2,
                              int *vla2,char *name1,char *name2,int identical);
int SelectorGetPairIndices(int sele1,int state1,int sele2,int state2,
                           int mode,float cutoff,float h_angle,
                           int **indexVLA, ObjectMolecule ***objVLA);

void SelectorLogSele(char *name);
int SelectorMapMaskVDW(int sele1,ObjectMapState *oMap,float buffer);

int SelectorMapCoulomb(int sele1,ObjectMapState *oMap,float cutoff);
int SelectorMapGaussian(int sele1,ObjectMapState *oMap,float buffer);
PyObject *SelectorAsPyList(int sele1);
int SelectorFromPyList(char *name,PyObject *list);

#endif
