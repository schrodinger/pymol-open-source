
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

#include <unordered_map>

#include"os_python.h"

#include"ObjectMolecule.h"
#include"CoordSet.h"
#include"DistSet.h"
#include"ObjectMap.h"
#include"OVOneToAny.h"
#include"Match.h"

#include "Result.h"

constexpr SelectorID_t cSelectionInvalid = -1;
constexpr SelectorID_t cSelectionAll = 0;
constexpr SelectorID_t cSelectionNone = 1;

int SelectorInit(PyMOLGlobals * G);
int SelectorInitImpl(PyMOLGlobals * G, CSelector **I, short init2);

/// return type of all SelectorCreate related functions
typedef pymol::Result<int> SelectorCreateResult_t;

SelectorCreateResult_t SelectorCreate(PyMOLGlobals * G, const char *name, const char *sele, ObjectMolecule * obj,
                   int quiet, Multipick * mp);
SelectorCreateResult_t SelectorCreateWithStateDomain(PyMOLGlobals * G, const char *name, const char *sele,
                                  ObjectMolecule * obj, int quiet, Multipick * mp,
                                  int state, const char *domain);
SelectorCreateResult_t SelectorCreateSimple(PyMOLGlobals * G, const char *name, const char *sele);
SelectorCreateResult_t SelectorCreateFromObjectIndices(PyMOLGlobals * G, const char *sname, ObjectMolecule * obj,
                                    int *idx, int n_idx);
SelectorCreateResult_t SelectorCreateOrderedFromObjectIndices(PyMOLGlobals * G, const char *sname,
                                           ObjectMolecule * obj, int *idx, int n_idx);

SelectorCreateResult_t SelectorCreateFromTagDict(PyMOLGlobals * G, const char *sname, const std::unordered_map<int, int>& id2tag,
                              int exec_managed);


/* if n_idx is negative, then looks for negative *idx as the sentinel */
bool SelectorMoveMember(PyMOLGlobals * G, SelectorMemberOffset_t s, SelectorID_t sele_old, SelectorID_t sele_new);
SelectorCreateResult_t SelectorCreateEmpty(PyMOLGlobals * G, const char *name, int exec_managed);

int SelectorUpdateTableImpl(PyMOLGlobals * G, CSelector *I, int req_state, SelectorID_t domain);

constexpr StateIndex_t cSelectorUpdateTableAllStates = cStateAll;
constexpr StateIndex_t cSelectorUpdateTableCurrentState = cStateCurrent;
constexpr StateIndex_t cSelectorUpdateTableEffectiveStates =
    -3; ///< deprecated, use cSelectorUpdateTableCurrentState instead

int SelectorUpdateTable(PyMOLGlobals*,
    StateIndex_t req_state = cSelectorUpdateTableAllStates,
    SelectorID_t domain = cSelectionInvalid);

SelectorID_t SelectorIndexByName(PyMOLGlobals * G, const char *sele, int ignore_case=-1);
const char *SelectorGetNameFromIndex(PyMOLGlobals * G, SelectorID_t index);
void SelectorFree(PyMOLGlobals * G);
void SelectorFreeImpl(PyMOLGlobals * G, CSelector *I, short init2);
void SelectorDelete(PyMOLGlobals * G, const char *sele);
void SelectorFreeTmp(PyMOLGlobals * G, const char *name);
int SelectorGetTmp2(PyMOLGlobals * G, const char *input, char *store, bool quiet=false);
int SelectorGetTmp(PyMOLGlobals * G, const char *input, char *store, bool quiet=false);

pymol::Result<int> SelectorGetTmp2Result(PyMOLGlobals * G, const char *input, char *store, bool quiet=false);
pymol::Result<int> SelectorGetTmpResult(PyMOLGlobals * G, const char *input, char *store, bool quiet=false);

pymol::Result<> SelectorLoadCoords(PyMOLGlobals * G, PyObject * coords, int sele, int state);
PyObject *SelectorGetCoordsAsNumPy(PyMOLGlobals * G, int sele, int state);
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
int SelectorGetSeleNCSet(PyMOLGlobals * G, SelectorID_t sele);
int SelectorCreateObjectMolecule(PyMOLGlobals * G, SelectorID_t sele, const char *name,
                                 int target_state, int state, int discrete,
                                 int zoom, int quiet, int singletons, int copy_properties);
int SelectorSubdivide(PyMOLGlobals * G, const char *pref, SelectorID_t sele1, SelectorID_t sele2,
                      SelectorID_t sele3, SelectorID_t sele4,
                      const char *fragPref, const char *compName, int *bondMode);
ObjectMolecule *SelectorGetSingleObjectMolecule(PyMOLGlobals * G, SelectorID_t sele);
ObjectMolecule *SelectorGetFirstObjectMolecule(PyMOLGlobals * G, SelectorID_t sele);
int SelectorRenameObjectAtoms(PyMOLGlobals* G, ObjectMolecule* obj,
    SelectorID_t sele, bool force, bool update_table);
void SelectorUpdateObjectSele(PyMOLGlobals * G, ObjectMolecule * obj);
void SelectorDeletePrefixSet(PyMOLGlobals * G, const char *pref);
pymol::Result<> SelectorUpdateCmd(PyMOLGlobals* G, SelectorID_t sele0, SelectorID_t sele1,
    int sta0, int sta1, int method, int quiet);
pymol::Result<std::array<float, 3>> SelectorGetSingleAtomVertex(PyMOLGlobals * G, int sele, int state);
pymol::Result<std::pair<ObjectMolecule*, int>> SelectorGetSingleAtomObjectIndex(PyMOLGlobals * G, SelectorID_t sele);
int *SelectorGetResidueVLA(PyMOLGlobals * G, SelectorID_t, int ca_only,
                           ObjectMolecule * exclude);
int SelectorCreateAlignments(PyMOLGlobals * G, int *pair, int sele1, int *vla1, int sele2,
                             int *vla2, const char *name1, const char *name2, int identical,
                             int atomic_input);
int SelectorGetPairIndices(PyMOLGlobals * G, int sele1, int state1, int sele2, int state2,
                           int mode, float cutoff, float h_angle, int **indexVLA,
                           ObjectMolecule *** objVLA);

int SelectorCountAtoms(PyMOLGlobals * G, int sele, int state);
void SelectorSetDeleteFlagOnSelectionInObject(PyMOLGlobals * G, int sele, ObjectMolecule *obj, signed char val);
int SelectorCheckIntersection(PyMOLGlobals * G, int sele1, int sele2);
int SelectorCountStates(PyMOLGlobals * G, int sele);
int SelectorClassifyAtoms(PyMOLGlobals * G, int sele, int preserve,
                          ObjectMolecule * only_object);

void SelectorLogSele(PyMOLGlobals * G, const char *name);
int SelectorMapMaskVDW(PyMOLGlobals * G, int sele1, ObjectMapState * oMap, float buffer,
                       int state);

int SelectorMapCoulomb(PyMOLGlobals * G, int sele1, ObjectMapState * oMap, float cutoff,
                       int state, int neutral, int shift, float shift_power);

int SelectorMapGaussian(PyMOLGlobals * G, int sele1, ObjectMapState * oMap,
                        float buffer, int state, int normalize, int use_max, int quiet,
                        float resolution);

PyObject *SelectorAsPyList(PyMOLGlobals * G, SelectorID_t sele1);
int SelectorFromPyList(PyMOLGlobals * G, const char *name, PyObject * list);
ObjectMolecule **SelectorGetObjectMoleculeVLA(PyMOLGlobals * G, SelectorID_t sele);

PyObject *SelectorColorectionGet(PyMOLGlobals * G, const char *prefix);
int SelectorColorectionApply(PyMOLGlobals * G, PyObject * list, const char *prefix);
int SelectorColorectionSetName(PyMOLGlobals * G, PyObject * list, const char *prefix,
                               char *new_prefix);
int SelectorColorectionFree(PyMOLGlobals * G, PyObject * list, const char *prefix);
void SelectorReinit(PyMOLGlobals * G);
PyObject *SelectorSecretsAsPyList(PyMOLGlobals * G);
int SelectorSecretsFromPyList(PyMOLGlobals * G, PyObject * list);
void SelectorMemoryDump(PyMOLGlobals * G);
int SelectorAssignSS(PyMOLGlobals * G, int target, int present, int state_value,
                     int preserve, ObjectMolecule * single_object, int quiet);

int SelectorPurgeObjectMembers(PyMOLGlobals * G, ObjectMolecule * obj);
void SelectorDefragment(PyMOLGlobals * G);

/**
 * Retrives unique selector name
 * @return a unique string to identify the selection
 */
std::string SelectorGetUniqueTmpName(PyMOLGlobals * G);
int SelectorIsAtomBondedToSele(PyMOLGlobals * G, ObjectMolecule * obj, int sele1atom,
                               int sele2);

int SelectorSetName(PyMOLGlobals * G, const char *new_name, const char *old_name);

ObjectMolecule *SelectorGetFastSingleAtomObjectIndex(PyMOLGlobals * G, SelectorID_t sele,
                                                     int *index);
ObjectMolecule *SelectorGetFastSingleObjectMolecule(PyMOLGlobals * G, SelectorID_t sele);
MapType *SelectorGetSpacialMapFromSeleCoord(PyMOLGlobals * G, int sele, int state,
                                            float cutoff, float **coord_vla);

/**
 * Determines whether a string is a reserved keyword
 * @param str string candidate
 * @return true if the string is a reserved keyword
 * @note: case-insensitive
 */
bool SelectorNameIsKeyword(PyMOLGlobals * G, const char *name);

int SelectorResidueVLAsTo3DMatchScores(PyMOLGlobals * G, CMatch * match,
                                       int *vla1, int n1, int state1,
                                       int *vla2, int n2, int state2,
                                       float seq_wt,
                                       float radius, float scale,
                                       float base, float coord_wt, float rms_exp);

int SelectorAssignAtomTypes(PyMOLGlobals * G, int sele, int state, int quiet, int format);


/* reserve special meaning for tags 1-15 and note that 0 is disallowed */

#define SELECTOR_BASE_TAG 0x10

typedef struct {
  SelectorID_t selection;
  int tag;                      /* must not be zero since it is also used as a boolean test for membership */
  SelectorMemberOffset_t next;
} MemberType;

int SelectorIsMember(PyMOLGlobals * G, SelectorMemberOffset_t, SelectorID_t);

/**
 * Wrapper around SelectorGetTmp/SelectorFreeTmp/SelectorIndexByName.
 *
 * Temporary named selection gets deleted when instance gets out of scope.
 */
class SelectorTmp {
protected:
  PyMOLGlobals* m_G = nullptr;
  char m_name[1024]{}; // OrthoLineType
  int m_count = -1;

public:
  SelectorTmp() = default;
  SelectorTmp(PyMOLGlobals * G, const char * sele) : m_G(G) {
    m_count = SelectorGetTmp(m_G, sele, m_name);
  }
  SelectorTmp(SelectorTmp&& other);
  SelectorTmp& operator=(SelectorTmp&& other) {
    std::swap(m_G, other.m_G);
    std::swap(m_count, other.m_count);
    std::swap(m_name, other.m_name);
    return *this;
  }
  ~SelectorTmp() {
    SelectorFreeTmp(m_G, m_name);
  }
  const char * getName() const { return m_name; }
  int getAtomCount() { return m_count; }
  SelectorID_t getIndex() const {
    return m_name[0] ? SelectorIndexByName(m_G, m_name, false) : cSelectionInvalid;
  }
  // Factory which propagages errors
  static pymol::Result<SelectorTmp> make(
      PyMOLGlobals* G, const char* sele, bool empty_is_error = true);
};

struct SelectorTmp2 : SelectorTmp {
  SelectorTmp2() = default;
  SelectorTmp2(PyMOLGlobals* G, const char* sele)
  {
    m_G = G;
    m_count = SelectorGetTmp2(m_G, sele, m_name);
  }
  // Factory which propagages errors
  static pymol::Result<SelectorTmp2> make(
      PyMOLGlobals* G, const char* sele, bool empty_is_error = false);
};

#endif
