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

int SelectorInit(PyMOLGlobals *G);
int *SelectorSelect(PyMOLGlobals *G,char *sele);
int SelectorCreate(PyMOLGlobals *G,char *name,char *sele,ObjectMolecule *obj,int quiet,Multipick *mp);
int SelectorCreateSimple(PyMOLGlobals *G,char *name, char *sele);
int SelectorCreateOrderedFromObjectIndices(PyMOLGlobals *G,char *sname, ObjectMolecule *obj, int *idx, int n_idx); 
/* if n_idx is negative, then looks for negative *idx as the sentinel */
int SelectorMoveMember(PyMOLGlobals *G,int s,int sele_old,int sele_new);
int SelectorCreateEmpty(PyMOLGlobals *G,char *name);
void SelectorToggle(PyMOLGlobals *G,int rep,char *name);
void SelectorCylinder(PyMOLGlobals *G,char *sele,char *onoff);
int SelectorUpdateTable(PyMOLGlobals *G);
int SelectorIndexByName(PyMOLGlobals *G,char *sele);
void SelectorFree(PyMOLGlobals *G);
void SelectorDelete(PyMOLGlobals *G,char *sele);
void SelectorFreeTmp(PyMOLGlobals *G,char *name);
int SelectorGetTmp(PyMOLGlobals *G,char *input,char *store);
int SelectorGetPDB(PyMOLGlobals *G,char **charVLA,int cLen,int sele,int state,
                   int conectFlag,PDBInfoRec *pdb_info);
PyObject *SelectorGetChemPyModel(PyMOLGlobals *G,int sele,int state);
float SelectorSumVDWOverlap(PyMOLGlobals *G,int sele1,int state1,int sele2,int state2,float adjust);
DistSet *SelectorGetDistSet(PyMOLGlobals *G,int sele1,int state1,int sele2,int state2,int mode,
                            float cutoff,float *result);
int SelectorGetSeleNCSet(PyMOLGlobals *G,int sele);
void SelectorCreateObjectMolecule(PyMOLGlobals *G,int sele,char *name,int target_state,int state,int discrete);
int SelectorSubdivide(PyMOLGlobals *G,char *pref,int sele1,int sele2,
                            int sele3,int sele4,
                            char *fragPref,char *compName,int *bondMode);
ObjectMolecule *SelectorGetSingleObjectMolecule(PyMOLGlobals *G,int sele);
void SelectorUpdateObjectSele(PyMOLGlobals *G,ObjectMolecule *obj);
void SelectorDeletePrefixSet(PyMOLGlobals *G,char *pref);
void SelectorUpdateCmd(PyMOLGlobals *G,int sele0,int sele1,int sta0,int sta1);
int SelectorGetSingleAtomVertex(PyMOLGlobals *G,int sele,int state,float *v);
int SelectorGetSingleAtomObjectIndex(PyMOLGlobals *G,int sele,ObjectMolecule **in_obj,int *index);
int *SelectorGetResidueVLA(PyMOLGlobals *G,int sele0);
int  SelectorCreateAlignments(PyMOLGlobals *G,int *pair,int sele1,int *vla1,int sele2,
                              int *vla2,char *name1,char *name2,int identical);
int SelectorGetPairIndices(PyMOLGlobals *G,int sele1,int state1,int sele2,int state2,
                           int mode,float cutoff,float h_angle,
                           int **indexVLA, ObjectMolecule ***objVLA);

int SelectorCountAtoms(PyMOLGlobals *G,int sele);
int SelectorCountStates(PyMOLGlobals *G,int sele);
int SelectorClassifyAtoms(PyMOLGlobals *G,int sele, int preserve,ObjectMolecule *only_object);


void SelectorLogSele(PyMOLGlobals *G,char *name);
int SelectorMapMaskVDW(PyMOLGlobals *G,int sele1,ObjectMapState *oMap,float buffer,int state);

int SelectorMapCoulomb(PyMOLGlobals *G,int sele1,ObjectMapState *oMap,float cutoff,int state,
                       int neutral,int shift,float shift_power);

int SelectorMapGaussian(PyMOLGlobals *G,int sele1,ObjectMapState *oMap,float buffer,int state);
PyObject *SelectorAsPyList(PyMOLGlobals *G,int sele1);
int SelectorFromPyList(PyMOLGlobals *G,char *name,PyObject *list);
ObjectMolecule **SelectorGetObjectMoleculeVLA(PyMOLGlobals *G,int sele);

PyObject *SelectorColorectionGet(PyMOLGlobals *G,char *prefix);
int SelectorColorectionApply(PyMOLGlobals *G,PyObject *list,char *prefix);
int SelectorColorectionFree(PyMOLGlobals *G,PyObject *list,char *prefix);
void SelectorReinit(PyMOLGlobals *G);
PyObject *SelectorSecretsAsPyList(PyMOLGlobals *G);
int SelectorSecretsFromPyList(PyMOLGlobals *G,PyObject *list);
void SelectorMemoryDump(PyMOLGlobals *G);
int SelectorAssignSS(PyMOLGlobals *G,int target,int present,int state_value,int preserve,int quiet);

int SelectorPurgeObjectMembers(PyMOLGlobals *G,ObjectMolecule *obj);
void SelectorDefragment(PyMOLGlobals *G);
void SelectorSelectByID(PyMOLGlobals *G,char *name,ObjectMolecule *obj,int *id,int n_id);
void SelectorGetUniqueTmpName(PyMOLGlobals *G,char *name_buffer);
int SelectorIsAtomBondedToSele(PyMOLGlobals *G,ObjectMolecule *obj,int sele1atom,int sele2);
void SelectorComputeFragPos(PyMOLGlobals *G,ObjectMolecule *obj,int state,int n_frag, char *prefix,float **vla);

int SelectorSetName(PyMOLGlobals *G,char *new_name, char *old_name);

ObjectMolecule *SelectorGetCachedSingleAtom(PyMOLGlobals *G,int sele,int *theAtom);

ObjectMolecule *SelectorGetFastSingleAtomObjectIndex(PyMOLGlobals *G,int sele,int *index);
ObjectMolecule *SelectorGetFastSingleObjectMolecule(PyMOLGlobals *G,int sele);
MapType *SelectorGetSpacialMapFromSeleCoord(PyMOLGlobals *G,int sele,int state,float cutoff,float **coord_vla);



#ifndef _PYMOL_INLINE

int SelectorIsMemberSlow(PyMOLGlobals *G,int start,int sele);
#define SelectorIsMember SelectorIsMemberSlow

#else

#ifdef _PYMOL_WIN32
#define __inline__ __inline
#endif

int _SelectorIsMemberInlinePartial(PyMOLGlobals *G,int s,int sele);

__inline__ static int SelectorIsMember(PyMOLGlobals *G,int s, int sele) 
{
  if(sele>1)   
    return _SelectorIsMemberInlinePartial(G,s,sele);
  else if(!sele) 
    return true; /* "all" is selection number 0, unordered */
  else 
    return false; /* no atom is a member of none (1), and negative selections don't exist */
}
#endif

#endif
