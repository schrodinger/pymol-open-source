
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
#include"OVOneToAny.h"
#include"Match.h"

#define cSelectionAll 0
#define cSelectionNone 1

int SelectorInit(PyMOLGlobals * G);
int SelectorCreate(PyMOLGlobals * G, char *name, char *sele, ObjectMolecule * obj,
                   int quiet, Multipick * mp);
int SelectorCreateWithStateDomain(PyMOLGlobals * G, char *name, char *sele,
                                  ObjectMolecule * obj, int quiet, Multipick * mp,
                                  int state, char *domain);
int SelectorCreateSimple(PyMOLGlobals * G, char *name, char *sele);
int SelectorCreateFromObjectIndices(PyMOLGlobals * G, char *sname, ObjectMolecule * obj,
                                    int *idx, int n_idx);
int SelectorCreateOrderedFromObjectIndices(PyMOLGlobals * G, char *sname,
                                           ObjectMolecule * obj, int *idx, int n_idx);
int SelectorCreateOrderedFromMultiObjectIdxTag(PyMOLGlobals * G, char *sname,
                                               ObjectMolecule ** obj, int **pri_idx,
                                               int *n_idx, int n_obj);

int SelectorCreateFromTagDict(PyMOLGlobals * G, char *sname, OVOneToAny * id2tag,
                              int exec_managed);


/* if n_idx is negative, then looks for negative *idx as the sentinel */
int SelectorMoveMember(PyMOLGlobals * G, int s, int sele_old, int sele_new);
int SelectorCreateEmpty(PyMOLGlobals * G, char *name, int exec_managed);
void SelectorToggle(PyMOLGlobals * G, int rep, char *name);
void SelectorCylinder(PyMOLGlobals * G, char *sele, char *onoff);

int SelectorUpdateTable(PyMOLGlobals * G, int req_state, int domain);
#define cSelectorUpdateTableAllStates -1
#define cSelectorUpdateTableCurrentState -2
#define cSelectorUpdateTableEffectiveStates -3

int SelectorIndexByName(PyMOLGlobals * G, char *sele);
char *SelectorGetNameFromIndex(PyMOLGlobals * G, int index);
void SelectorFree(PyMOLGlobals * G);
void SelectorDelete(PyMOLGlobals * G, char *sele);
void SelectorFreeTmp(PyMOLGlobals * G, char *name);
int SelectorGetTmp(PyMOLGlobals * G, char *input, char *store);
int SelectorCheckTmp(PyMOLGlobals * G, char *name);
int SelectorGetPDB(PyMOLGlobals * G, char **charVLA, int cLen, int sele, int state,
                   int conectFlag, PDBInfoRec * pdb_info, int *counter, double *ref,
                   ObjectMolecule * single_object);
PyObject *SelectorGetChemPyModel(PyMOLGlobals * G, int sele, int state, double *ref);
float SelectorSumVDWOverlap(PyMOLGlobals * G, int sele1, int state1,
                            int sele2, int state2, float adjust);
int SelectorVdwFit(PyMOLGlobals * G, int sele1, int state1, int sele2, int state2,
                   float buffer, int quiet);

DistSet *SelectorGetDistSet(PyMOLGlobals * G, DistSet * ds,
                            int sele1, int state1, int sele2,
                            int state2, int mode, float cutoff, float *result);
DistSet *SelectorGetAngleSet(PyMOLGlobals * G, DistSet * ds,
                             int sele1, int state1,
                             int sele2, int state2,
                             int sele3, int state3,
                             int mode, float *angle_sum, int *angle_cnt);
DistSet *SelectorGetDihedralSet(PyMOLGlobals * G, DistSet * ds,
                                int sele1, int state1,
                                int sele2, int state2,
                                int sele3, int state3,
                                int sele4, int state4,
                                int mode, float *angle_sum, int *angle_cnt);
int SelectorGetSeleNCSet(PyMOLGlobals * G, int sele);
int SelectorCreateObjectMolecule(PyMOLGlobals * G, int sele, char *name,
                                 int target_state, int state, int discrete,
                                 int zoom, int quiet, int singletons);
int SelectorSubdivide(PyMOLGlobals * G, char *pref, int sele1, int sele2,
                      int sele3, int sele4,
                      char *fragPref, char *compName, int *bondMode);
ObjectMolecule *SelectorGetSingleObjectMolecule(PyMOLGlobals * G, int sele);
ObjectMolecule *SelectorGetFirstObjectMolecule(PyMOLGlobals * G, int sele);
int SelectorRenameObjectAtoms(PyMOLGlobals * G, ObjectMolecule * obj, int sele, int force,
                              int update_table);
void SelectorUpdateObjectSele(PyMOLGlobals * G, ObjectMolecule * obj);
void SelectorDeletePrefixSet(PyMOLGlobals * G, char *pref);
void SelectorUpdateCmd(PyMOLGlobals * G, int sele0, int sele1, int sta0, int sta1,
                       int method, int quiet);
int SelectorGetSingleAtomVertex(PyMOLGlobals * G, int sele, int state, float *v);
int SelectorGetSingleAtomObjectIndex(PyMOLGlobals * G, int sele, ObjectMolecule ** in_obj,
                                     int *index);
int *SelectorGetResidueVLA(PyMOLGlobals * G, int sele0, int ca_only,
                           ObjectMolecule * exclude);
int SelectorCreateAlignments(PyMOLGlobals * G, int *pair, int sele1, int *vla1, int sele2,
                             int *vla2, char *name1, char *name2, int identical,
                             int atomic_input);
int SelectorGetPairIndices(PyMOLGlobals * G, int sele1, int state1, int sele2, int state2,
                           int mode, float cutoff, float h_angle, int **indexVLA,
                           ObjectMolecule *** objVLA);

int SelectorCountAtoms(PyMOLGlobals * G, int sele, int state);
int SelectorCheckIntersection(PyMOLGlobals * G, int sele1, int sele2);
int SelectorCountStates(PyMOLGlobals * G, int sele);
int SelectorClassifyAtoms(PyMOLGlobals * G, int sele, int preserve,
                          ObjectMolecule * only_object);

void SelectorLogSele(PyMOLGlobals * G, char *name);
int SelectorMapMaskVDW(PyMOLGlobals * G, int sele1, ObjectMapState * oMap, float buffer,
                       int state);

int SelectorMapCoulomb(PyMOLGlobals * G, int sele1, ObjectMapState * oMap, float cutoff,
                       int state, int neutral, int shift, float shift_power);

int SelectorMapGaussian(PyMOLGlobals * G, int sele1, ObjectMapState * oMap,
                        float buffer, int state, int normalize, int use_max, int quiet,
                        float resolution);

PyObject *SelectorAsPyList(PyMOLGlobals * G, int sele1);
int SelectorFromPyList(PyMOLGlobals * G, char *name, PyObject * list);
ObjectMolecule **SelectorGetObjectMoleculeVLA(PyMOLGlobals * G, int sele);

PyObject *SelectorColorectionGet(PyMOLGlobals * G, char *prefix);
int SelectorColorectionApply(PyMOLGlobals * G, PyObject * list, char *prefix);
int SelectorColorectionSetName(PyMOLGlobals * G, PyObject * list, char *prefix,
                               char *new_prefix);
int SelectorColorectionFree(PyMOLGlobals * G, PyObject * list, char *prefix);
void SelectorReinit(PyMOLGlobals * G);
PyObject *SelectorSecretsAsPyList(PyMOLGlobals * G);
int SelectorSecretsFromPyList(PyMOLGlobals * G, PyObject * list);
void SelectorMemoryDump(PyMOLGlobals * G);
int SelectorAssignSS(PyMOLGlobals * G, int target, int present, int state_value,
                     int preserve, ObjectMolecule * single_object, int quiet);

int SelectorPurgeObjectMembers(PyMOLGlobals * G, ObjectMolecule * obj);
void SelectorDefragment(PyMOLGlobals * G);
void SelectorSelectByID(PyMOLGlobals * G, char *name, ObjectMolecule * obj, int *id,
                        int n_id);
void SelectorGetUniqueTmpName(PyMOLGlobals * G, char *name_buffer);
int SelectorIsAtomBondedToSele(PyMOLGlobals * G, ObjectMolecule * obj, int sele1atom,
                               int sele2);
void SelectorComputeFragPos(PyMOLGlobals * G, ObjectMolecule * obj, int state, int n_frag,
                            char *prefix, float **vla);

int SelectorSetName(PyMOLGlobals * G, char *new_name, char *old_name);

ObjectMolecule *SelectorGetCachedSingleAtom(PyMOLGlobals * G, int sele, int *theAtom);

ObjectMolecule *SelectorGetFastSingleAtomObjectIndex(PyMOLGlobals * G, int sele,
                                                     int *index);
ObjectMolecule *SelectorGetFastSingleObjectMolecule(PyMOLGlobals * G, int sele);
MapType *SelectorGetSpacialMapFromSeleCoord(PyMOLGlobals * G, int sele, int state,
                                            float cutoff, float **coord_vla);
int SelectorNameIsKeyword(PyMOLGlobals * G, char *name);

int SelectorResidueVLAsTo3DMatchScores(PyMOLGlobals * G, CMatch * match,
                                       int *vla1, int n1, int state1,
                                       int *vla2, int n2, int state2,
                                       float seq_wt,
                                       float radius, float scale,
                                       float base, float coord_wt, float rms_exp);

PyObject *SelectorAssignAtomTypes(PyMOLGlobals * G, int sele, int state, int quiet, int format);


/* reserve special meaning for tags 1-15 and note that 0 is disallowed */

#define SELECTOR_BASE_TAG 0x10

typedef struct {
  int selection;
  int tag;                      /* must not be zero since it is also used as a boolean test for membership */
  int next;
} MemberType;

#ifndef _PYMOL_INLINE

int SelectorIsMemberSlow(PyMOLGlobals * G, int start, int sele);
#define SelectorIsMember SelectorIsMemberSlow

#else


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _PYMOL_WIN32
#define __inline__ __inline
#endif

/* END PROPRIETARY CODE SEGMENT */

__inline__ static int SelectorIsMember(PyMOLGlobals * G, int s, int sele)
{
  /* this is the most heavily called routine in interactive PyMOL */
  register int s_reg;
  if((s_reg = s) && (sele > 1)) {
    /* the first entry of a Selector is the Member (list) pointer, so
     * we access it via a ptr, and cast it to the MemberType */
    register MemberType *member = *((MemberType **) (G->Selector));
    register int sele_reg = sele;
    register MemberType *mem = member + s_reg;
    register int test_sele;
    do {
      test_sele = mem->selection;
      s_reg = mem->next;
      if(test_sele == sele_reg) {
        return mem->tag;
      }
      mem = member + s_reg;
    } while(s_reg);
    return false;
  } else if(!sele)
    return true;                /* "all" is selection number 0, unordered */
  else
    return false;               /* no atom is a member of none (1), and negative selections don't exist */
}
#endif

#endif
