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

#define cSelectionAll 0
#define cSelectionNone 1

void SelectorInit(void);
int *SelectorSelect(char *sele);
int SelectorCreate(char *name,char *sele,ObjectMolecule *obj,int quiet,Multipick *mp);
int SelectorCreateSimple(char *name, char *sele);

int SelectorCreateOrderedFromObjectIndices(char *sname, ObjectMolecule *obj, int *idx, int n_idx); 
/* if n_idx is negative, then looks for negative *idx as the sentinel */

int SelectorCreateEmpty(char *name);
void SelectorToggle(int rep,char *name);
void SelectorCylinder(char *sele,char *onoff);
int SelectorUpdateTable(void);
int SelectorIndexByName(char *sele);
int SelectorIsMember(int start,int sele);
void SelectorFree(void);
void SelectorDelete(char *sele);
void SelectorFreeTmp(char *name);
int SelectorGetTmp(char *input,char *store);
int SelectorGetPDB(char **charVLA,int cLen,int sele,int state,
                   int conectFlag,PDBInfoRec *pdb_info);
PyObject *SelectorGetChemPyModel(int sele,int state);
float SelectorSumVDWOverlap(int sele1,int state1,int sele2,int state2,float adjust);
DistSet *SelectorGetDistSet(int sele1,int state1,int sele2,int state2,int mode,
                            float cutoff,float *result);
int SelectorGetSeleNCSet(int sele);
void SelectorCreateObjectMolecule(int sele,char *name,int target_state,int state,int discrete);
int SelectorSubdivide(char *pref,int sele1,int sele2,
                            int sele3,int sele4,
                            char *fragPref,char *compName,int *bondMode);
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

int SelectorCountAtoms(int sele);
int SelectorCountStates(int sele);
int SelectorClassifyAtoms(int sele, int preserve,ObjectMolecule *only_object);


void SelectorLogSele(char *name);
int SelectorMapMaskVDW(int sele1,ObjectMapState *oMap,float buffer,int state);

int SelectorMapCoulomb(int sele1,ObjectMapState *oMap,float cutoff,int state);
int SelectorMapGaussian(int sele1,ObjectMapState *oMap,float buffer,int state);
PyObject *SelectorAsPyList(int sele1);
int SelectorFromPyList(char *name,PyObject *list);
ObjectMolecule **SelectorGetObjectMoleculeVLA(int sele);

PyObject *SelectorColorectionGet(char *prefix);
int SelectorColorectionApply(PyObject *list,char *prefix);
int SelectorColorectionFree(PyObject *list,char *prefix);
void SelectorReinit(void);
PyObject *SelectorSecretsAsPyList(void);
int SelectorSecretsFromPyList(PyObject *list);
void SelectorMemoryDump(void);
int SelectorAssignSS(int target,int present,int state_value,int preserve,int quiet);

int SelectorPurgeObjectMembers(ObjectMolecule *obj);
void SelectorDefragment(void);
void SelectorSelectByID(char *name,ObjectMolecule *obj,int *id,int n_id);
void SelectorGetUniqueTmpName(char *name_buffer);
int SelectorIsAtomBondedToSele(ObjectMolecule *obj,int sele1atom,int sele2);
void SelectorComputeFragPos(ObjectMolecule *obj,int state,int n_frag, char *prefix,float **vla);

int SelectorSetName(char *new_name, char *old_name);

ObjectMolecule *SelectorGetCachedSingleAtom(int sele,int *theAtom);

ObjectMolecule *SelectorGetFastSingleAtomObjectIndex(int sele,int *index);
ObjectMolecule *SelectorGetFastSingleObjectMolecule(int sele);
MapType *SelectorGetSpacialMapFromSeleCoord(int sele,int state,float cutoff,float **coord_vla);


#endif
