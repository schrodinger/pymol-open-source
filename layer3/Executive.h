
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
#ifndef _H_Executive
#define _H_Executive

#include <string>
#include "pymol/zstring_view.h"
#include"os_python.h"

#include"PyMOLGlobals.h"
#include"PyMOLObject.h"
#include"ObjectMolecule.h"
#include"Field.h"

#include"ExecutiveDef.h"
#include"Tracker.h"
#include"Ortho.h"
#include"Word.h"
#include "PyMOL.h"
#include "Executive_pre.h"
#include "Scene.h"
#include "Result.h"
#include "vla.h"
#include "TrackerList.h"
#include "Selector.h"
#include "SpecRecSpecial.h"

enum cLoadType_t : int {
  cLoadTypeUnknown = -1,
  cLoadTypePDB = 0,
  cLoadTypeMOL = 1,
  cLoadTypeSDF1 = 2, /* SDF1 - python-based loader */
  cLoadTypeMOLStr = 3,
  cLoadTypeMMD = 4,
  cLoadTypeMMDSeparate = 5,
  cLoadTypeMMDStr = 6,
  cLoadTypeXPLORMap = 7,
  cLoadTypeChemPyModel = 8,
  cLoadTypePDBStr = 9,
  cLoadTypeChemPyBrick = 10,
  cLoadTypeChemPyMap = 11,
  cLoadTypeCallback = 12,
  cLoadTypeCGO = 13,
  cLoadTypeR3D = 14,
  cLoadTypeXYZ = 15,
  cLoadTypeCCP4Map = 18,
  cLoadTypePMO = 19,
  cLoadTypeTOP = 21,
  cLoadTypeTRJ = 22,
  cLoadTypeCRD = 23,
  cLoadTypeRST = 24,
  cLoadTypePSE = 25,
  cLoadTypeXPLORStr = 26,
  cLoadTypePHIMap = 27,
  cLoadTypeFLDMap = 28,
  cLoadTypeBRIXMap = 29,
  cLoadTypeGRDMap = 30,
  cLoadTypePQR = 31,
  cLoadTypeDXMap = 32,
  cLoadTypeMOL2 = 33,
  cLoadTypeMOL2Str = 34,
  cLoadTypeP1M = 35,
  cLoadTypeCCP4Str = 36,
  cLoadTypeSDF2 = 37,
  cLoadTypeSDF2Str = 38,
  cLoadTypeXTC = 42,
  cLoadTypeTRR = 43,
  cLoadTypeGRO = 44,
  cLoadTypeTRJ2 = 45,
  cLoadTypeG96 = 46,
  cLoadTypeDCD = 47,
  cLoadTypeCUBEMap = 48,
  cLoadTypeXYZStr = 49,

/* 50 is Python-based CIF (cif1) */
  cLoadTypePHIStr = 51,

/* 52 is PIM */

/* 53 is PWG */

/* 54 is ALN */

/* 55 is FASTA */
  cLoadTypeACNTMap = 56,

  cLoadTypeDTR = 57,
/* 58/59 are pze, pzw */
  cLoadTypeCIF = 60,
  cLoadTypeCIFStr = 61,

  cLoadTypeSpider = 62,
  cLoadTypeCMS = 63,
  cLoadTypePlugin = 64,

  cLoadTypeMAE = 65,
  cLoadTypeMAEStr = 66,
  cLoadTypePDBQT = 67,
/* 68 is Python-based MAE  */
  cLoadTypeVDBStr = 69,

  cLoadTypeMMTF = 71,
  cLoadTypeMMTFStr = 72,

  // Distinguish CCP4-like map types:
  // - cLoadTypeCCP4Str/cLoadTypeCCP4Map : CCP4 format, ignore MRC features
  // - cLoadTypeCCP4Unspecified          : CCP4 or MRC format (auto-detect)
  // - cLoadTypeMRC                      : MRC format, ignore CCP4 features
  cLoadTypeCCP4Unspecified = 73,
  cLoadTypeMRC = 74,

  cLoadTypeDXStr = 75,

  cLoadTypeCCP4UnspecifiedStr = 76,
  cLoadTypeMRCStr = 77,
};

/* NOTE: if you add new content/object type above, then be sure to add
   corresponding code in:
   ExecutiveGetExistingCompatible
   ExecutiveLoad
*/

struct ExecutiveObjectOffset{
  ObjectMolecule *obj;
  int atm;
};

class SpecRec;

/**
 * Iterator over objects (uses SpecRec list)
 */
class ObjectIterator {
  PyMOLGlobals * G;
  SpecRec *rec;

public:
  ObjectIterator(PyMOLGlobals * G) : G(G) {
    reset();
  };

  void reset();
  bool next();
  pymol::CObject* getObject();
  SpecRec * getSpecRec() { return rec; }
};

const char * ExecutiveMapGenerate(PyMOLGlobals * G, const char * name, const char * reflection_file, const char * tempFile,
				  char * amplitudes, const char * phases, const char * weights, double reso_low,
				  double reso_high, const char * space_group, double cell[6], int quiet, int zoom);
				  
int ExecutiveReference(PyMOLGlobals * G, int action, const char *sele, int state, int quiet);
int ExecutiveGetExpandedGroupList(PyMOLGlobals * G, const char *name);
int ExecutiveGetExpandedGroupListFromPattern(PyMOLGlobals * G, const char *name);
void ExecutiveFreeGroupList(PyMOLGlobals * G, int list_id);

int ExecutiveCheckGroupMembership(
    PyMOLGlobals* G, int list_id, pymol::CObject* obj); /* 0.5*N for group size */

#define cExecutiveGroupAdd 1
#define cExecutiveGroupRemove 2
#define cExecutiveGroupOpen 3
#define cExecutiveGroupClose 4
#define cExecutiveGroupToggle 5
#define cExecutiveGroupAuto 6
#define cExecutiveGroupUngroup 7
#define cExecutiveGroupEmpty 8
#define cExecutiveGroupPurge 9
#define cExecutiveGroupExcise 10

enum ExecRec_t
{
  cExecObject = 0,
  cExecSelection = 1,
  cExecAll = 2
};

int ExecutiveGroup(PyMOLGlobals * G, pymol::zstring_view name, pymol::zstring_view members, int action, int quiet);
int ExecutiveScrollTo(PyMOLGlobals * G, const char * name, int);

void ExecutiveInvalidateGroups(PyMOLGlobals * G, bool force);
void ExecutiveUpdateGroups(PyMOLGlobals * G, bool force);
void ExecutiveUpdateSceneMembers(PyMOLGlobals*);

int *ExecutiveGetG3d(PyMOLGlobals * G);
pymol::Result<> ExecutiveOrder(PyMOLGlobals * G, pymol::zstring_view s1, int sort, int location);
pymol::Result<> ExecutiveFixChemistry(
    PyMOLGlobals* G, const char* s1, const char* s2, int invalidate, int quiet);
pymol::Result<std::array<float, 3>> ExecutiveGetAtomVertex(PyMOLGlobals * G, const char *s1, int state, int index);
int ExecutiveProcessPDBFile(PyMOLGlobals * G, pymol::CObject * origObj,
                            const char *fname, const char *buffer, const char *oname,
                            int frame, int discrete, int finish, OrthoLineType buf,
                            int variant, int quiet,
                            int multiplex, int zoom);

const ExecutiveObjectOffset * ExecutiveUniqueIDAtomDictGet(PyMOLGlobals * G, int i);
void ExecutiveUniqueIDAtomDictInvalidate(PyMOLGlobals * G);

pymol::Result<> ExecutiveLoad(
                  PyMOLGlobals * G,
                  const char *fname,
                  const char *content, int content_length,
                  cLoadType_t content_format,
                  const char *object_name,
                  int state, int zoom,
                  int discrete, int finish, int multiplex, int quiet, const char *plugin,
                  const char * object_props=NULL, const char * atom_props=NULL,
                  bool mimic=true);

int ExecutiveDebug(PyMOLGlobals * G, const char *name);

typedef struct {
  int n_residues_aligned;
  float raw_alignment_score;
  int initial_n_atom;
  float initial_rms;
  int n_cycles_run;
  int final_n_atom;
  float final_rms;
} ExecutiveRMSInfo;

std::string ExecutivePreparePseudoatomName(PyMOLGlobals* G, pymol::zstring_view object_name);

pymol::Result<> ExecutivePseudoatom(PyMOLGlobals* G, pymol::zstring_view object_name,
    const char* sele, const char* name, const char* resn, const char* resi,
    const char* chain, const char* segi, const char* elem, float vdw,
    int hetatm, float b, float q, const char* label, const float* pos, int color,
    int state, int mode, int quiet);

pymol::Result<> ExecutiveMapSet(PyMOLGlobals* G, const char* name, int,
    const char* operands, int target_state, int source_state, int zoom,
    int quiet);

int ExecutiveAlign(PyMOLGlobals * G, const char *s1, const char *s2, const char *mat_file,
                   float gap, float extend,
                   int max_gap, int max_skip,
                   float cutoff, int cycles, int quiet,
                   const char *oname, int state1, int state2,
                   ExecutiveRMSInfo * rms_info, int transform, int reset,
                   float seq_wt, float radius, float scale, float base,
                   float coord_wt, float expect, int window, float ante);

void ExecutiveUpdateColorDepends(PyMOLGlobals * G, ObjectMolecule * mol);
void ExecutiveUpdateCoordDepends(PyMOLGlobals * G, ObjectMolecule * mol);
pymol::Result<float> ExecutiveDistance(PyMOLGlobals* G, const char* nam,
    const char* s1, const char* s2, int mode, float cutoff, int labels,
    int quiet, int reset, int state, int zoom, int state1 = -4,
    int state2 = -4);
pymol::Result<> ExecutiveBond(PyMOLGlobals* G, const char* s1,
    const char* s2, int order, int mode, int quiet, pymol::zstring_view symop = "");
pymol::Result<> ExecutiveAddBondByIndices(PyMOLGlobals* G, pymol::zstring_view oname,
    unsigned int atm1, unsigned int atm2, int order);
pymol::Result<> ExecutiveRevalence(PyMOLGlobals* G, const char* s1,
    const char* s2, const char* src, int target_state, int source_state,
    int reset, int quiet);
int ExecutiveVdwFit(PyMOLGlobals * G, const char *s1, int state1, const char *s2, int state2,
                    float buffer, int quiet);
#ifdef _WEBGL
#else
pymol::Result<int> ExecutiveIterate(PyMOLGlobals * G, const char *str1, const char *expr, int read_only, int quiet,
                     PyObject * space);
#endif
pymol::Result<int> ExecutiveIterateList(PyMOLGlobals* G, const char* s1,
    PyObject* list, int read_only, int quiet, PyObject* space);

struct SelectArgs
{
  std::string sname;
  std::string sele;
};

SelectArgs ExecutiveSelectPrepareArgs(PyMOLGlobals* G, pymol::zstring_view sname, pymol::zstring_view sele);

pymol::Result<int> ExecutiveSelect(PyMOLGlobals* G, const SelectArgs& sargs,
    int enable, int quiet, int merge, int state,
    const char* domain);

pymol::Result<int> ExecutiveSelectList(PyMOLGlobals* G, const char* sele_name,
    const char* oname, const int* list, size_t list_len, int state, int mode,
    int quiet);

#define cExecutiveLabelEvalOff    0
#define cExecutiveLabelEvalOn     1
#define cExecutiveLabelEvalAlt    2

pymol::Result<> ExecutiveLabel(PyMOLGlobals * G, const char *s1, const char *expr, int quiet, int eval_mode);

int ExecutiveCountMotions(PyMOLGlobals * G);
#ifdef _WEBGL
#else
pymol::Result<int> ExecutiveIterateState(PyMOLGlobals* G, int state,
    const char* str1, const char* expr, int read_only,
    int quiet, PyObject* space);
#endif
pymol::Result<> ExecutiveColor(
    PyMOLGlobals* G, const char* name, const char* color, int flags, int quiet);
pymol::Result<> ExecutiveColorFromSele(
    PyMOLGlobals* G, const char* sele, const char* color, int flags, int quiet);
int ExecutiveInit(PyMOLGlobals * G);
void ExecutiveFree(PyMOLGlobals * G);
int ExecutivePop(PyMOLGlobals * G, const char *target, const char *source, int quiet);
void ExecutiveManageObject(
    PyMOLGlobals* G, pymol::CObject* obj, int allow_zoom, int quiet);
void ExecutiveUpdateObjectSelection(PyMOLGlobals* G, pymol::CObject* obj);
void ExecutiveManageSelection(PyMOLGlobals * G, const char *name);
Block *ExecutiveGetBlock(PyMOLGlobals * G);
pymol::CObject* ExecutiveFindObjectByName(
    PyMOLGlobals* G, pymol::zstring_view name);

/**
 * Typed version of `ExecutiveFindObjectByName`
 * @param name Object name
 * @return NULL if object can't be found or has the wrong type
 */
template <typename T> T* ExecutiveFindObject(PyMOLGlobals* G, pymol::zstring_view name)
{
  return dynamic_cast<T*>(ExecutiveFindObjectByName(G, name));
}

#define ExecutiveFindObjectMoleculeByName ExecutiveFindObject<ObjectMolecule>
#define ExecutiveFindObjectMapByName ExecutiveFindObject<ObjectMap>

pymol::CObject** ExecutiveFindObjectsByType(PyMOLGlobals * G, int objType);
int ExecutiveIterateObject(PyMOLGlobals * G, pymol::CObject ** obj, void **hidden);
pymol::Result<std::vector<DiscardedRec>> ExecutiveDelete(PyMOLGlobals * G, pymol::zstring_view nameView, bool save = false);

/**
 * @brief Unregisters the specification record from PyMOL
 * @param rec specification record to be purged/removed
 * @save if true, the rec's associated object/cgo will not be destroyed, and the caller is responsible for repurging w/o saving.
 */

void ExecutivePurgeSpec(PyMOLGlobals* G, SpecRec* rec, bool save = false);

/**
 * @brief Readds specs to Executive
 * @param G pointer to global pymol singletons
 * @param specs list of specs to be readded
 */

void ExecutiveReAddSpec(PyMOLGlobals* G, std::vector<DiscardedRec>& specs);
void ExecutiveDump(PyMOLGlobals * G, const char *fname, const char *obj, int state, int quiet);
pymol::Result<> ExecutiveSort(PyMOLGlobals * G, const char *name);
PyObject *ExecutiveGetBondSetting(PyMOLGlobals * G, int index, 
				  char *s1, const char *s2, int state, int quiet, int updates);
int ExecutiveSetBondSetting(PyMOLGlobals * G, int index, PyObject * tuple,
                            const char *s1, const char *s2, int state, int quiet, int updates);
int ExecutiveSetBondSettingFromString(PyMOLGlobals * G,
                                      int index, const char *value,
                                      const char *s1, const char *s2, int state,
                                      int quiet, int updates);

int ExecutiveUnsetBondSetting(PyMOLGlobals * G, int index, const char *s1, const char *s2,
                              int state, int quiet, int updates);

pymol::Result<> ExecutiveSetSetting(PyMOLGlobals * G, int index, PyObject * tuple,
                        pymol::zstring_view preSele, int state, int quiet, int updates);
int ExecutiveGetSettingFromString(PyMOLGlobals * G, PyMOLreturn_value *result, 
                                  int index, const char *sele,
                                  int state, int quiet);
int ExecutiveSetSettingFromString(PyMOLGlobals * G, int index, const char *value,
                                  const char *sele, int state, int quiet, int updates);
int ExecutiveSetObjSettingFromString(PyMOLGlobals * G,
                                     int index, const char *value, pymol::CObject * obj,
                                     int state, int quiet, int updates);
int ExecutiveRay(PyMOLGlobals * G, int width, int height, int mode,
                 float angle, float shift, int quiet, int defer, int antialias);
pymol::Result<float> ExecutiveGetDistance(
    PyMOLGlobals* G, const char* s0, const char* s1, int state);
pymol::Result<float> ExecutiveGetAngle(
    PyMOLGlobals* G, const char* s0, const char* s1, const char* s2, int state);
pymol::Result<float> ExecutiveGetDihe(PyMOLGlobals* G, const char* s0,
    const char* s1, const char* s2, const char* s3, int state);
pymol::Result<> ExecutiveSetDihe(PyMOLGlobals* G, const char* s0,
    const char* s1, const char* s2, const char* s3, float value, int state = 0,
    int quiet = 1);
int ExecutiveRMS(PyMOLGlobals * G, const char *sele1, const char *sele2, int mode, float refine,
                 int max_cyc, int quiet, const char *oname, int state1, int state2,
                 int ordered_selections, int matchmaker, ExecutiveRMSInfo * rms_info);

pymol::Result<> ExecutiveUpdateCmd(PyMOLGlobals* G, const char* sele1,
    const char* sele2, int sta1, int sta2, int method, int quiet);
float ExecutiveRMSPairs(PyMOLGlobals* G, const std::vector<SelectorTmp>& sele, int mode, bool quiet);
pymol::Result<pymol::vla<float>> ExecutiveRMSStates(PyMOLGlobals* G,
    const char* s1, int target, int mode, int quiet, int mix, bool pbc = true);
int ExecutiveIndex(PyMOLGlobals * G, const char *s1, int mode, int **indexVLA,
                   ObjectMolecule *** objVLA);
pymol::Result<> ExecutiveReset(PyMOLGlobals*, pymol::zstring_view);
pymol::Result<> ExecutiveResetMatrix(
    PyMOLGlobals* G, const char* name, int mode, int state, int log, int quiet);
void ExecutiveDrawNow(PyMOLGlobals * G);
int ExecutiveDrawCmd(PyMOLGlobals * G, int width, int height, int antialias,
                     int entire_window, int quiet);
pymol::Result<int> ExecutiveCartoon(PyMOLGlobals* G, int type, const char* s1);
pymol::Result<> ExecutiveSetRepVisib(PyMOLGlobals * G, pymol::zstring_view name, int rep, int state);
pymol::Result<> ExecutiveSetRepVisMask(PyMOLGlobals * G, pymol::zstring_view name, int repmask, int state);
pymol::Result<> ExecutiveSetRepVisMaskFromSele(PyMOLGlobals* G, pymol::zstring_view sele, int repmask, int state);
pymol::Result<> ExecutiveToggleRepVisib(PyMOLGlobals * G, const char *name, int rep);

/**
 * Sets the visibility of an object/selection
 *
 * @param name name pattern
 * @param onoff to activate or deactivate recs referenced in name
 * @param parents if true, also activate their parent rec
 * @return true if any visibility has been changed.
 */

pymol::Result<bool> ExecutiveSetObjVisib(PyMOLGlobals * G, pymol::zstring_view name, int onoff, int parents);


/**
 * Sets the rotation center for the scene or object around a selection's center (or provided coords)
 * @param sele (null-safe) selection expression used to calculate rotation center
 * @param preserve preserves the current viewing location
 * @param oname (null-safe) object name
 * @param pos optional explicit coordinates for rotation center
 * @param state state(s) considered to calculate rotation center
 * Note: Either seletion or explicit center position must be provided
 * Note: selection's center is preferred over explicit pos if both are provided
 */
pymol::Result<> ExecutiveOrigin(PyMOLGlobals* G, const char* sele, int preserve,
    const char* oname, const float* pos, int state);
pymol::Result<> ExecutiveCenter(PyMOLGlobals* G, const char* name, int state,
    int inclusive, float animate, float* pos, int quiet);
void ExecutiveDoZoom(PyMOLGlobals * G, pymol::CObject * obj, int is_new, int zoom, int quiet);
pymol::Result<> ExecutiveWindowZoom(PyMOLGlobals* G,
    const char* name, float buffer, int state, int inclusive, float animate,
    int quiet);
int ExecutiveGetMoment(PyMOLGlobals* G, const char* name, double* mi, int state);

pymol::Result<std::vector<const char*>> ExecutiveGetChains(PyMOLGlobals* G, const char* sele, int state);

pymol::Result<> ExecutiveOrient(PyMOLGlobals* G, const char* sele,
    int state, float animate, int complete, float buffer, int quiet);
pymol::Result<> ExecutiveMove(
    PyMOLGlobals* G, pymol::zstring_view axis, float dist);
char *ExecutiveNameToSeqAlignStrVLA(PyMOLGlobals * G, const char *name, int state, int format,
                                    int quiet);

pymol::Result<> ExecutiveCopy(
    PyMOLGlobals* G, const char* src, const char* dst, int zoom);
pymol::Result<> ExecutiveStereo(PyMOLGlobals * G, int flag);
float ExecutiveOverlap(PyMOLGlobals * G, const char *s1, int state1, const char *s2, int state2,
                       float adjust);
int ExecutiveCountStates(PyMOLGlobals * G, const char *s1);
void ExecutiveSymExp(PyMOLGlobals * G, const char *name, const char *obj, const char *sele, float cutoff,
                     int segi, int quiet);
int ExecutiveGetExtent(PyMOLGlobals * G, const char *name, float *mn, float *mx,
                       int transformed, int state, int weighted);
int ExecutiveGetCameraExtent(PyMOLGlobals * G, const char *name, float *mn, float *mx,
                             int transformed, int state);
pymol::Result<> ExecutiveSeleToObject(PyMOLGlobals* G, const char* name,
    const char* s1, int source, int target, int discrete, int zoom, int quiet,
    int singletons, int copy_properties = 0);

PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals * G, const char *s1, int state,
                                     const char *ref_object, int ref_state);
pymol::Result<> ExecutiveInvalidateRep(
    PyMOLGlobals* G, const char* name, cRep_t rep, cRepInv_t level);
pymol::Result<> ExecutiveFlag(PyMOLGlobals * G, int flag, const char *sele, int action, int quiet);
pymol::Result<> ExecutiveRemoveAtoms(PyMOLGlobals * G, const char *s1, int quiet);
pymol::Result<> ExecutiveProtect(PyMOLGlobals * G, const char *s1, int mode, int quiet);
pymol::Result<> ExecutiveMask(
    PyMOLGlobals* G, const char* s1, int mode = 1, int quiet = 1);
void ExecutiveRebuildAll(PyMOLGlobals * G);
pymol::Result<> ExecutiveSpheroid(PyMOLGlobals * G, const char *name, int average);
pymol::Result<> ExecutiveAddHydrogens(PyMOLGlobals* G, const char* s1 = "(all)",
    int quiet = 1, int state = cStateAll, bool legacy = false);
void ExecutiveFixHydrogens(PyMOLGlobals * G, const char *s1, int quiet);
pymol::Result<> ExecutiveFuse(PyMOLGlobals* G, const char* s0 = "(pk1)",
    const char* s1 = "(pk2)", int mode = 0, int recolor = 1, int move_flag = 1);
pymol::Result<> ExecutiveRenameObjectAtoms(
    PyMOLGlobals* G, const char* name, int force, int quiet);

pymol::Result<std::vector<const char*>> ExecutiveGetNames(PyMOLGlobals*, int, int, const char*);
bool ExecutiveIsMoleculeOrSelection(PyMOLGlobals * G, const char *name);
pymol::Result<char const*> ExecutiveGetType(PyMOLGlobals* G, const char* name);

pymol::Result<float> ExecutiveGetArea(
    PyMOLGlobals*, const char* sele, int state, bool load_b);

void ExecutiveInvalidateSceneMembers(PyMOLGlobals * G);
void ExecutiveInvalidateSelectionIndicators(PyMOLGlobals *G);
void ExecutiveInvalidateSelectionIndicatorsCGO(PyMOLGlobals *G);
void ExecutiveRenderSelections(PyMOLGlobals * G, int curState, int slot, GridInfo *grid);
void ExecutiveHideSelections(PyMOLGlobals * G);
pymol::Result<> ExecutiveSetTitle(PyMOLGlobals * G, const char *name, int state, const char *text);
const char *ExecutiveGetTitle(PyMOLGlobals * G, const char *name, int state);
void ExecutiveSetLastObjectEdited(PyMOLGlobals * G, pymol::CObject * o);
pymol::CObject* ExecutiveGetLastObjectEdited(PyMOLGlobals* G);
bool ExecutiveIsFullScreen(PyMOLGlobals * G);
void ExecutiveFullScreen(PyMOLGlobals * G, int flag);

#ifndef _PYMOL_NOPY
PyObject *ExecutiveGetSettingOfType(PyMOLGlobals * G, int index, const char *object, int state,
                                    int type);
#endif

pymol::vla<ObjectMolecule*> ExecutiveGetObjectMoleculeVLA(PyMOLGlobals * G, const char *sele);
pymol::Result<int> ExecutivePairIndices(PyMOLGlobals * G, const char *s1, const char *s2, int state1, int state2,
                         int mode, float cutoff, float h_angle,
                         int **indexVLA, ObjectMolecule *** objVLA);
void ExecutiveRebuildAllObjectDist(PyMOLGlobals * G);
int ExecutivePhiPsi(PyMOLGlobals * G, const char *s1, ObjectMolecule *** objVLA, int **iVLA,
                    float **phiVLA, float **psiVLA, int state);
float *ExecutiveGetVertexVLA(PyMOLGlobals * G, const char *s1, int state);
int ExecutiveValidName(PyMOLGlobals * G, const char *name);
int ExecutiveValidNamePattern(PyMOLGlobals * G, const char *name);
void ExecutiveMakeUnusedName(PyMOLGlobals * G, char * prefix, int length, bool alwaysnumber=true, int start=1, const char * pattern="%02d");
std::string ExecutiveGetUnusedName(PyMOLGlobals * G, const char * prefix="tmp", bool alwaysnumber=true);
int ExecutiveProcessObjectName(PyMOLGlobals * G, const char *proposed, char *actual);

pymol::Result<> ExecutiveIsolevel(PyMOLGlobals* G, const char* name, float level, int state, int quiet);
pymol::Result<float> ExecutiveGetIsolevel(PyMOLGlobals* G, const char* name, int state);
pymol::Result<> ExecutiveTransformObjectSelection(PyMOLGlobals* G,
    const char* name, int state, const char* s1, int log, const float* matrix,
    int homogenous, int global);
pymol::Result<> ExecutiveTransformSelection(PyMOLGlobals* G,
    int state, const char* s1, int log, const float* ttt, int homogenous);
pymol::Result<> ExecutiveTranslateAtom(
    PyMOLGlobals* G, const char* sele, const float* v, int state, int mode, int log);

/**
 * Contains both the number of selected atoms returned by SelectorCreate
 * and the change in selected atoms compared to before the call.
 */

struct NetSelect
{
  int currSelected = 0;
  int netSelected = 0;
};

/**
 * Selects atoms via screen-space 2D rect
 * @param rect rectangle that encloses desired atoms
 * @param mode button mode
 * @return number of atoms selected in active selection and change in that selected
 */

pymol::Result<NetSelect> ExecutiveSelectRect(PyMOLGlobals * G, BlockRect * rect, int mode);
int ExecutiveMapSetBorder(PyMOLGlobals * G, const char *name, float level, int state);
pymol::Result<> ExecutiveMapTrim(PyMOLGlobals* G, const char* name,
    const char* sele, float buffer, int map_state, int sele_state, int quiet);
pymol::Result<> ExecutiveMapDouble(PyMOLGlobals * G, const char *name, int state);
pymol::Result<> ExecutiveMapHalve(PyMOLGlobals * G, const char *name, int state, int smooth);

int ExecutiveIdentifyObjects(PyMOLGlobals * G, const char *s1, int mode, int **indexVLA,
                             ObjectMolecule *** objVLA);
pymol::Result<> ExecutiveTranslateObjectTTT(PyMOLGlobals* G,
    pymol::zstring_view name, const float* trans, int store, int quiet);
pymol::Result<> ExecutiveCombineObjectTTT(PyMOLGlobals* G,
    pymol::zstring_view name, const float* ttt, int reverse_order, int store);
pymol::Result<> ExecutiveSetObjectTTT(PyMOLGlobals* G,
    pymol::zstring_view name, const float* ttt, int state, int quiet, int store);
int ExecutiveGetObjectTTT(PyMOLGlobals * G, const char *name, const float **ttt, int state,
                          int quiet);
int ExecutiveGetObjectMatrix(PyMOLGlobals * G, const char *name, int state, double **matrix,
                             int incl_ttt);
int ExecutiveSetObjectMatrix(PyMOLGlobals * G, const char *name, int state, double *matrix);

pymol::Result<> ExecutiveSetGeometry(
    PyMOLGlobals* G, const char* s1, int geom, int valence);
int ExecutiveSculptIterateAll(PyMOLGlobals * G);
pymol::Result<> ExecutiveSmooth(PyMOLGlobals* G, const char* name, int cycles,
    int window, int first, int last, int ends, int quiet,
    float dist_cutoff = -1, bool pbc = true);
int ExecutiveSculptDeactivate(PyMOLGlobals * G, const char *name);
int ExecutiveSculptActivate(PyMOLGlobals* G, const char* name,
    int state = cStateAll, int match_state = cStateCurrent,
    int match_by_segment = 0);
float ExecutiveSculptIterate(
    PyMOLGlobals* G, const char* name, int state = cStateAll, int n_cycle = 10);
pymol::Result<> ExecutiveMapNew(PyMOLGlobals* G, const char* name, int type,
    float grid_spacing, const char* sele, float buffer, const float* minCorner,
    const float* maxCorner, int state, int have_corners, int quiet, int zoom,
    int normalize, float clamp_floor, float clamp_ceiling, float resolution);

int ***ExecutiveGetBondPrint(PyMOLGlobals * G, const char *name, int max_bond, int max_type,
                             int *dim);

pymol::Result<bool>
ExecutiveGetSymmetry(PyMOLGlobals * G, const char *sele, int state, float *a, float *b, float *c,
                        float *alpha, float *beta, float *gamma, char *sgroup);
pymol::Result<> ExecutiveSetSymmetry(PyMOLGlobals* G, const char* sele,
    int state, float a, float b, float c, float alpha, float beta, float gamma,
    const char* sgroup, int quiet);
pymol::Result<> ExecutiveSymmetryCopy(PyMOLGlobals * G,
			  const char *source_name, const char *target_name,
			  int source_state, int target_state, int quiet);
int ExecutiveGetSession(PyMOLGlobals * G, PyObject * dict, const char *names, int partial,
                        int quiet);
int ExecutiveSetSession(PyMOLGlobals * G, PyObject * session, int partial_restore,
                        int quiet);
int ExecutiveSetSessionNoMLock(PyMOLGlobals* G, PyObject* session);

pymol::Result<> ExecutiveUnsetSetting(PyMOLGlobals * G, int index, pymol::zstring_view preSele,
                          int state, int quiet, int updates);

pymol::Result<> ExecutiveAssignSS(PyMOLGlobals* G,
    const char* target, int state, const char* context, int preserve,
    ObjectMolecule* single_object, int quiet);

pymol::Result<> ExecutiveRampNew(PyMOLGlobals* G, const char* name,
    const char* src_name, pymol::vla<float> range, pymol::vla<float> color,
    int src_state, const char* sele, float beyond, float within, float sigma,
    int zero, int calc_mode, int quiet);

int ExecutiveValidateObjectPtr(
    PyMOLGlobals* G, pymol::CObject* ptr, int object_type);

/**
 * Colors atoms with a spectrum of colors based on atomic property
 * @param G pointer to global pymol singletons
 * @param s1 selection of atoms to color
 * @param expr property expression
 * @param min minimum of color range
 * @param max maximum of color range
 * @param first color idx of first color
 * @param last color idx of last color
 * @param prefix specifies spectrum mode
 * @param digits (unknown)
 * @param byres determines whether coloring is applied per-residue
 * @param quiet determines whether feedback is displayed
 * @return pair of minimum and maximum of spectrum range
 */
pymol::Result<std::pair<float, float>> ExecutiveSpectrum(PyMOLGlobals* G,
    pymol::zstring_view s1, pymol::zstring_view expr, float min, float max, int first, int last,
    pymol::zstring_view prefix, int digits, int byres, int quiet);

pymol::Result<> ExecutiveReinitialize(PyMOLGlobals * G, int what, pymol::zstring_view pattern);
const char *ExecutiveFindBestNameMatch(PyMOLGlobals * G, const char *name);
int ExecutiveSetVisFromPyDict(PyMOLGlobals * G, PyObject * dict);
PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals * G);
CField   *ExecutiveGetVolumeField(PyMOLGlobals * G, const char * objName, int state);
pymol::Result<> ExecutiveSetVolumeRamp(PyMOLGlobals * G, const char * objName, std::vector<float> ramp_list);
PyObject *ExecutiveGetVolumeRamp(PyMOLGlobals * G, const char * objName);

pymol::Result<std::vector<float>>
ExecutiveGetHistogram(PyMOLGlobals * G, const char * objName, int n_points,
        float min_val, float max_val);

int ExecutiveIterateObjectMolecule(PyMOLGlobals * G, ObjectMolecule ** obj,
                                   void **hidden);
pymol::Result<> ExecutiveSetObjectColor(PyMOLGlobals * G, const char *name, const char *color, int quiet);
int ExecutiveGetObjectColorIndex(PyMOLGlobals * G, const char *name);
pymol::Result<> ExecutiveSetOnOffBySele(PyMOLGlobals * G, pymol::zstring_view sname, int onoff);

/**
 * Changes a rec's name to another.
 *
 * @param old_name target rec name
 * @param new_name rec's new name
 * @param save discarded recs will be saved
 * @return discarded recs if saved
 *
 * Note: Caller is responsible for calling ExecutivePurgeSpec on discarded recs
 */

pymol::Result<std::vector<DiscardedRec>>ExecutiveSetName(
        PyMOLGlobals * G, pymol::zstring_view old_name, pymol::zstring_view new_name, bool save = false);
int ExecutiveSetDrag(PyMOLGlobals * G, const char *name, int quiet,int mode);
int ExecutiveGetActiveSeleName(PyMOLGlobals * G, char *name, int create_new, int log);
int ExecutiveGetActiveSeleName(PyMOLGlobals * G, std::string& name, int create_new, int log);
int ExecutiveGetActiveSele(PyMOLGlobals * G);
int ExecutiveGetActiveAlignmentSele(PyMOLGlobals * G);
const char* ExecutiveGetActiveAlignment(PyMOLGlobals* G);
pymol::CObject* ExecutiveGetExistingCompatible(
    PyMOLGlobals* G, const char* oname, cLoadType_t type);
pymol::Result<float> ExecutiveAngle(PyMOLGlobals* G, const char* nam,
    const char* s1, const char* s2, const char* s3, int mode, int labels,
    int reset, int zoom, int quiet, int state, int state1 = -4, int state2 = -4,
    int state3 = -4);

pymol::Result<float> ExecutiveDihedral(PyMOLGlobals* G, const char* nam,
    const char* s1, const char* s2, const char* s3, const char* s4, int mode,
    int labels, int reset, int zoom, int quiet, int state);

int ExecutiveMatrixCopy2(PyMOLGlobals * G, pymol::CObject* source_obj,
    pymol::CObject* target_obj,
                         int source_mode, int target_mode,
                         int source_state, int target_state,
                         int target_undo, int log, int quiet);

int ExecutiveMatrixCopy(PyMOLGlobals* G, const char* source_name,
    const char* target_name, int source_mode, int target_mode, int source_state,
    int target_state, int target_undo, int log, int quiet);

void ExecutiveMemoryDump(PyMOLGlobals * G);
bool ExecutiveObjMolSeleOp(PyMOLGlobals * G, int sele, ObjectMoleculeOpRec * op);

pymol::Result<> ExecutiveIsomeshEtc(PyMOLGlobals * G,
                        const char *mesh_name, const char *map_name, float lvl,
                        const char *sele, float fbuf, int state,
                        float carve, int map_state, int quiet,
                        int mesh_mode, float alt_lvl);

pymol::Result<> ExecutiveIsosurfaceEtc(PyMOLGlobals * G,
                           const char *surf_name, const char *map_name, float lvl,
                           const char *sele, float fbuf, int state,
                           float carve, int map_state, int side,
                           int quiet, int surf_mode);

pymol::Result<> ExecutiveVolume(PyMOLGlobals * G, const char *volume_name, const char *map_name,
		    float lvl,
		    const char *sele, float fbuf, int state,
		    float carve, int map_state, int quiet);

int ExecutiveVolumeColor(PyMOLGlobals * G, const char * volume_name, float * colors, int ncolors );

void ExecutiveMotionDraw(PyMOLGlobals * G, BlockRect *rect, int expected ORTHOCGOARG);

void ExecutiveMotionReinterpolate(PyMOLGlobals * G);

pymol::Result<> ExecutiveMotionViewModify(PyMOLGlobals* G, int action,
    int index, int count, int target, const char* name, int freeze, int quiet);

void ExecutiveMotionMenuActivate(PyMOLGlobals * G, BlockRect *rect, int expected,int passive, 
                                 int x, int y, int same);
int ExecutiveGroupMotionModify(PyMOLGlobals *G, pymol::CObject *group, int action, 
                               int index, int count, int target, int freeze);
int ExecutiveGroupMotion(PyMOLGlobals *G, pymol::CObject *group,int action, int first,
                         int last, float power, float bias,
                         int simple, float linear, int wrap,
                         int hand, int window, int cycles, int state, int quiet);
int ExecutiveGroupCombineTTT(PyMOLGlobals *G, pymol::CObject *group, const float *ttt, int reverse_order, int store);
int ExecutiveGroupTranslateTTT(PyMOLGlobals *G, pymol::CObject *group, const float *v, int store);
void ExecutiveMotionClick(PyMOLGlobals * G, BlockRect *rect, int mode, int expected, int x, int y, int nearest);
void ExecutiveMotionTrim(PyMOLGlobals * G);
void ExecutiveMotionExtend(PyMOLGlobals * G, int freeze);
int ExecutiveMotionView(PyMOLGlobals *G, int action, int first,
                        int last, float power, float bias,
                        int simple, float linear, const char *name, int wrap,
                        int hand, int window, int cycles,
                        const char *scene_name, float scene_cut, int state, int quiet, int autogen);
int ExecutiveAssignAtomTypes(PyMOLGlobals * G, const char *s1, int format, int state, int quiet);

PyObject * ExecutiveCEAlign(PyMOLGlobals * G, PyObject * listA, PyObject * listB, int lenA, int lenB,
			    float d0, float d1, int windowSize, int gapMax);

pymol::Result<> ExecutiveSetFeedbackMask(
    PyMOLGlobals* G, int action, unsigned int sysmod, unsigned char mask);
pymol::Result<> ExecutiveClip(PyMOLGlobals* G, pymol::zstring_view clipStr);

pymol::Result<> ExecutiveSliceNew(PyMOLGlobals* G, const char* slice_name,
    const char* map_name, int state, int map_state);

pymol::Result<ExecutiveRMSInfo> ExecutiveFit(PyMOLGlobals* G,
    pymol::zstring_view str1, pymol::zstring_view str2, int mode, int cutoff,
    int cycles, int quiet, pymol::zstring_view object, int state1, int state2,
    int matchmaker);

#ifdef _PYMOL_LIB
int *ExecutiveGetRepsInSceneForObject(PyMOLGlobals *G, const char *name);
int *ExecutiveGetRepsForObject(PyMOLGlobals *G, const char *name);
#endif

int ExecutiveGetNamesListFromPattern(PyMOLGlobals * G, const char *name,
                                     int allow_partial, int expand_groups);

/**
 * Retrieves a list of candidates provided by a pattern. Similar to
 * ExecutiveGetNamesListFromPattern but returns a managed tracker list.
 *
 * @param str pattern string provided by user
 * @param allow_partial allows for partial name matching
 * @param expand_groups members of the group candidate will be expanded onto the
 * list.
 * @return a managed list of candidate records
 */
pymol::TrackerAdapter<SpecRec> ExecutiveGetSpecRecsFromPattern(PyMOLGlobals* G,
    pymol::zstring_view str, bool allow_partial = true,
    bool expand_groups = true);

char *ExecutiveGetObjectNames(PyMOLGlobals * G, int mode, const char *name, int enabled_only, int *numstrs);

CoordSet * ExecutiveGetCoordSet(PyMOLGlobals * G, const char * name, int state, ObjectMolecule ** omp=NULL);
pymol::Result<> ExecutiveLoadCoordset(
    PyMOLGlobals* G, pymol::zstring_view oname, PyObject* model, int frame);

void ExecutiveUndo(PyMOLGlobals * G, int dir);
int ExecutiveSaveUndo(PyMOLGlobals * G, const char *s1, int state);

pymol::Result<> ExecutiveRebond(
    PyMOLGlobals* G, const char* oname, int state, bool pbc = false);

/**
 * Determines whether the given name is of the Executive type (Selection or Object) provided
 * @param name name of object or selection
 * @param execType type to check against
 * @return true if name is of the type execType
 */

bool ExecutiveIsSpecRecType(PyMOLGlobals* G, pymol::zstring_view name, int execType);

void ExecutiveInvalidateMapDependents(
    PyMOLGlobals* G, const char* map_name, const char* new_name = nullptr);
pymol::Result<> ExecutiveLoadTraj(PyMOLGlobals* G, pymol::zstring_view oname,
    pymol::zstring_view fname, int frame, int type, int interval, int average,
    int start, int stop, int max, pymol::zstring_view str1, int image,
    const float* shift, pymol::zstring_view plugin, int quiet);

pymol::TrackerAdapter<SpecRec> ExecutiveGetSpecRecParents(
    PyMOLGlobals* G, SpecRec& rec);
SpecRec* ExecutiveFindSpec(PyMOLGlobals* G, pymol::zstring_view name_view);

/**
 * Retrives a list rec names a part of a group rec
 * @param groupName group rec name
 * @return space-delimited list of recs
 */

std::string ExecutiveGetGroupMemberNames(PyMOLGlobals* G, pymol::zstring_view groupName);

/**
 * Determines of the rec is of an object type
 * @param rec target rec
 * @param cObjectType object type enum
 */

bool ExecutiveIsObjectType(const SpecRec& rec, int cObjectType);

/**
 * Retrives order of Specification records in list
 * @param nameList expression for recs
 * @return rec names with their position in the Spec Rec list
 */

std::vector<OrderRec> ExecutiveGetOrderOf(PyMOLGlobals* G, pymol::zstring_view nameList);

/**
 * Sets recs in the Specification records in list at certain positions
 * @param recs recs to be reordered in list
 */

void ExecutiveSetOrderOf(PyMOLGlobals* G, const std::vector<OrderRec>& recs);
void ExecutiveInvalidatePanelList(PyMOLGlobals * G);
void ExecutiveSpecSetVisibility(PyMOLGlobals * G, SpecRec * rec,
                                int new_vis, int mod, int parents);

void ExecutiveSpecSetVisibility(PyMOLGlobals * G, SpecRec * rec,
                                int new_vis, int mod, int parents);

/**
 * Sets background color
 * @param color new bg color
 */
pymol::Result<> ExecutiveBackgroundColor(PyMOLGlobals* G, pymol::zstring_view color);
#endif
