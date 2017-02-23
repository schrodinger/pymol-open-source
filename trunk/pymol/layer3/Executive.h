
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

#include"os_python.h"

#include"PyMOLGlobals.h"
#include"PyMOLObject.h"
#include"ObjectMolecule.h"
#include"ObjectMap.h"

#include"Ortho.h"
#include"Word.h"
#include "PyMOL.h"
#include "Executive_pre.h"
#include "Scene.h"

#define cLoadTypeUnknown -1
#define cLoadTypePDB 0
#define cLoadTypeMOL 1
#define cLoadTypeSDF1 2         /* SDF1 - python-based loader */
#define cLoadTypeMOLStr 3
#define cLoadTypeMMD 4
#define cLoadTypeMMDSeparate 5
#define cLoadTypeMMDStr 6
#define cLoadTypeXPLORMap 7
#define cLoadTypeChemPyModel 8
#define cLoadTypePDBStr 9
#define cLoadTypeChemPyBrick 10
#define cLoadTypeChemPyMap 11
#define cLoadTypeCallback 12
#define cLoadTypeCGO 13
#define cLoadTypeR3D 14
#define cLoadTypeXYZ 15
#define cLoadTypeCCP4Map 18
#define cLoadTypePMO  19
#define cLoadTypeTOP  21
#define cLoadTypeTRJ  22
#define cLoadTypeCRD  23
#define cLoadTypeRST  24
#define cLoadTypePSE  25
#define cLoadTypeXPLORStr 26
#define cLoadTypePHIMap 27
#define cLoadTypeFLDMap 28
#define cLoadTypeBRIXMap 29
#define cLoadTypeGRDMap 30
#define cLoadTypePQR 31
#define cLoadTypeDXMap 32
#define cLoadTypeMOL2 33
#define cLoadTypeMOL2Str 34
#define cLoadTypeP1M 35
#define cLoadTypeCCP4Str 36
#define cLoadTypeSDF2 37
#define cLoadTypeSDF2Str 38
#define cLoadTypeXTC 42
#define cLoadTypeTRR 43
#define cLoadTypeGRO 44
#define cLoadTypeTRJ2 45
#define cLoadTypeG96 46
#define cLoadTypeDCD 47
#define cLoadTypeCUBEMap 48
#define cLoadTypeXYZStr 49

/* 50 is Python-based CIF (cif1) */
#define cLoadTypePHIStr 51

/* 52 is PIM */

/* 53 is PWG */

/* 54 is ALN */

/* 55 is FASTA */
#define cLoadTypeACNTMap 56

#define cLoadTypeDTR 57
/* 58/59 are pze, pzw */
#define cLoadTypeCIF 60
#define cLoadTypeCIFStr 61

#define cLoadTypeSpider 62
#define cLoadTypeCMS 63
#define cLoadTypePlugin 64

#define cLoadTypeMAE 65
#define cLoadTypeMAEStr 66
#define cLoadTypePDBQT 67
/* 68 is Python-based MAE  */
#define cLoadTypeVDBStr 69

#define cLoadTypeMMTF 71
#define cLoadTypeMMTFStr 72

/* NOTE: if you add new content/object type above, then be sure to add
   corresponding code in:
   ExecutiveGetExistingCompatible
   ExecutiveLoad
*/

typedef struct {
  ObjectMolecule *obj;
  int atm;
} ExecutiveObjectOffset;

class SpecRec;

/*
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
  CObject * getObject();
};


const char * ExecutiveMapGenerate(PyMOLGlobals * G, const char * name, const char * reflection_file, const char * tempFile,
				  char * amplitudes, const char * phases, const char * weights, double reso_low,
				  double reso_high, const char * space_group, double cell[6], int quiet, int zoom);
				  
int ExecutiveReference(PyMOLGlobals * G, int action, const char *sele, int state, int quiet);
int ExecutiveGetExpandedGroupList(PyMOLGlobals * G, const char *name);
int ExecutiveGetExpandedGroupListFromPattern(PyMOLGlobals * G, const char *name);
void ExecutiveFreeGroupList(PyMOLGlobals * G, int list_id);

int ExecutiveCheckGroupMembership(PyMOLGlobals * G, int list_id, CObject * obj);        /* 0.5*N for group size */

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

int ExecutiveGroup(PyMOLGlobals * G, const char *name, const char *members, int action, int quiet);

void ExecutiveInvalidateGroups(PyMOLGlobals * G, int force);
void ExecutiveUpdateGroups(PyMOLGlobals * G, int force);

int *ExecutiveGetG3d(PyMOLGlobals * G);
int ExecutiveOrder(PyMOLGlobals * G, const char *s1, int sort, int location);
int ExecutiveFixChemistry(PyMOLGlobals * G, const char *s1, const char *s2, int invalidate,
                          int quiet);
int ExecutiveGetAtomVertex(PyMOLGlobals * G, const char *s1, int state, int index, float *v);
int ExecutiveProcessPDBFile(PyMOLGlobals * G, CObject * origObj,
                            const char *fname, const char *buffer, const char *oname,
                            int frame, int discrete, int finish, OrthoLineType buf,
                            int variant, int quiet,
                            int multiplex, int zoom);

const ExecutiveObjectOffset * ExecutiveUniqueIDAtomDictGet(PyMOLGlobals * G, int i);
void ExecutiveUniqueIDAtomDictInvalidate(PyMOLGlobals * G);

int ExecutiveLoad(PyMOLGlobals * G,
                  const char *content, int content_length,
                  int content_format,
                  const char *object_name,
                  int state, int zoom,
                  int discrete, int finish, int multiplex, int quiet, const char *plugin,
                  const char * object_props=NULL, const char * atom_props=NULL);

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

int ExecutivePseudoatom(PyMOLGlobals * G, const char *object_name, const char *sele,
                        const char *name, const char *resn, const char *resi, const char *chain,
                        const char *segi, const char *elem, float vdw, int hetatm,
                        float b, float q, const char *label, float *pos, int color,
                        int state, int mode, int quiet);

int ExecutiveMapSet(PyMOLGlobals * G, const char *name, int, const char *operands,
                    int target_state, int source_state, int zoom, int quiet);

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
int ExecutiveDist(PyMOLGlobals * G, float *result,
                  const char *nam, const char *s1, const char *s2, int mode,
                  float cutoff, int labels, int quiet, int reset, int state, int zoom);
int ExecutiveBond(PyMOLGlobals * G, const char *s1, const char *s2, int order, int mode, int quiet);
int ExecutiveRevalence(PyMOLGlobals * G, const char *s1, const char *s2, const char *src,
                       int target_state, int source_state, int reset, int quiet);
int ExecutiveVdwFit(PyMOLGlobals * G, const char *s1, int state1, const char *s2, int state2,
                    float buffer, int quiet);
#ifdef _WEBGL
#else
int ExecutiveIterate(PyMOLGlobals * G, const char *str1, const char *expr, int read_only, int quiet,
                     PyObject * space);
#endif
int ExecutiveIterateList(PyMOLGlobals * G, const char *s1, PyObject * list, int read_only,
                         int quiet, PyObject * space);
int ExecutiveSelectList(PyMOLGlobals * G, const char *sele_name, const char *s1, int *list,
                        int list_len, int state, int mode, int quiet);

#define cExecutiveLabelEvalOff    0
#define cExecutiveLabelEvalOn     1
#define cExecutiveLabelEvalAlt    2

int ExecutiveLabel(PyMOLGlobals * G, const char *s1, const char *expr, int quiet, int eval_mode);

int ExecutiveCountMotions(PyMOLGlobals * G);
#ifdef _WEBGL
#else
int ExecutiveIterateState(PyMOLGlobals * G, int state, const char *str1, const char *expr, int read_only,
                           int atomic_props, int quiet, PyObject * space);
#endif
int ExecutiveColor(PyMOLGlobals * G, const char *name, const char *color, int flags, int quiet);
int ExecutiveInit(PyMOLGlobals * G);
void ExecutiveFree(PyMOLGlobals * G);
int ExecutivePop(PyMOLGlobals * G, const char *target, const char *source, int quiet);
void ExecutiveManageObject(PyMOLGlobals * G, CObject * obj, int allow_zoom, int quiet);
void ExecutiveUpdateObjectSelection(PyMOLGlobals * G, CObject * obj);
void ExecutiveManageSelection(PyMOLGlobals * G, const char *name);
Block *ExecutiveGetBlock(PyMOLGlobals * G);
CObject *ExecutiveFindObjectByName(PyMOLGlobals * G, const char *name);
ObjectMolecule *ExecutiveFindObjectMoleculeByName(PyMOLGlobals * G, const char *name);
CObject ** ExecutiveFindObjectsByType(PyMOLGlobals * G, int objType);
int ExecutiveIterateObject(PyMOLGlobals * G, CObject ** obj, void **hidden);
void ExecutiveDelete(PyMOLGlobals * G, const char *name);
void ExecutiveDump(PyMOLGlobals * G, const char *fname, const char *obj);
void ExecutiveSort(PyMOLGlobals * G, const char *name);
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

int ExecutiveSetSetting(PyMOLGlobals * G, int index, PyObject * tuple,
                        const char *sele, int state, int quiet, int updates);
int ExecutiveGetSettingFromString(PyMOLGlobals * G, PyMOLreturn_value *result, 
                                  int index, const char *sele,
                                  int state, int quiet);
int ExecutiveSetSettingFromString(PyMOLGlobals * G, int index, const char *value,
                                  const char *sele, int state, int quiet, int updates);
int ExecutiveSetObjSettingFromString(PyMOLGlobals * G,
                                     int index, const char *value, CObject * obj,
                                     int state, int quiet, int updates);
int ExecutiveRay(PyMOLGlobals * G, int width, int height, int mode,
                 float angle, float shift, int quiet, int defer, int antialias);
int ExecutiveGetDistance(PyMOLGlobals * G, const char *s0, const char *s1, float *value, int state);
int ExecutiveGetAngle(PyMOLGlobals * G, const char *s0, const char *s1, const char *s2, float *value,
                      int state);
int ExecutiveGetDihe(PyMOLGlobals * G, const char *s0, const char *s1, const char *s2, const char *s3,
                     float *value, int state);
int ExecutiveSetDihe(PyMOLGlobals * G, const char *s0, const char *s1, const char *s2, const char *s3,
                     float value, int state, int quiet);
int ExecutiveRMS(PyMOLGlobals * G, const char *sele1, const char *sele2, int mode, float refine,
                 int max_cyc, int quiet, const char *oname, int state1, int state2,
                 int ordered_selections, int matchmaker, ExecutiveRMSInfo * rms_info);

void ExecutiveUpdateCmd(PyMOLGlobals * G, const char *sele1, const char *sele2, int sta1, int sta2,
                        int method, int quiet);
float ExecutiveRMSPairs(PyMOLGlobals * G, WordType * sele, int pairs, int mode);
float *ExecutiveRMSStates(PyMOLGlobals * G, const char *s1, int target, int mode, int quiet,
                          int mix);
int *ExecutiveIdentify(PyMOLGlobals * G, const char *s1, int mode);
int ExecutiveIndex(PyMOLGlobals * G, const char *s1, int mode, int **indexVLA,
                   ObjectMolecule *** objVLA);
int ExecutiveReset(PyMOLGlobals * G, int cmd, const char *name);
void ExecutiveResetMatrix(PyMOLGlobals * G,
                          const char *name, int mode, int state, int log, int quiet);
void ExecutiveDrawNow(PyMOLGlobals * G);
int ExecutiveDrawCmd(PyMOLGlobals * G, int width, int height, int antialias,
                     int entire_window, int quiet);
int ExecutiveCartoon(PyMOLGlobals * G, int type, const char *sele);
void ExecutiveSetAllVisib(PyMOLGlobals * G, int state);
void ExecutiveSetRepVisib(PyMOLGlobals * G, const char *name, int rep, int state);
int ExecutiveToggleRepVisib(PyMOLGlobals * G, const char *name, int rep);

int ExecutiveSetObjVisib(PyMOLGlobals * G, const char *name, int onoff, int parents);

int ExecutiveOrigin(PyMOLGlobals * G, const char *name, int preserve, const char *oname, float *pos,
                    int state);
int ExecutiveCenter(PyMOLGlobals * G, const char *name, int state, int inclusive, float animate,
                    float *pos, int quiet);
void ExecutiveDoZoom(PyMOLGlobals * G, CObject * obj, int is_new, int zoom, int quiet);
int ExecutiveWindowZoom(PyMOLGlobals * G, const char *name, float buffer,
                        int state, int inclusive, float animate, int quiet);
int ExecutiveGetMoment(PyMOLGlobals * G, const char *name, double *mi, int state);

const char **ExecutiveGetChains(PyMOLGlobals * G, const char *sele, int state);

void ExecutiveOrient(PyMOLGlobals * G, const char *sele, double *mi,
                     int state, float animate, int complete, float buffer, int quiet);
char *ExecutiveSeleToPDBStr(PyMOLGlobals * G, const char *s1, int state, int conectFlag,
                            int mode, const char *ref_object, int ref_state,
                            ObjectMolecule * single_object, int quiet);
char *ExecutiveNameToSeqAlignStrVLA(PyMOLGlobals * G, const char *name, int state, int format,
                                    int quiet);

int ExecutiveStereo(PyMOLGlobals * G, int flag);
void ExecutiveCopy(PyMOLGlobals * G, const char *src, const char *dst, int zoom);
float ExecutiveOverlap(PyMOLGlobals * G, const char *s1, int state1, const char *s2, int state2,
                       float adjust);
int ExecutiveCountStates(PyMOLGlobals * G, const char *s1);
void ExecutiveSymExp(PyMOLGlobals * G, const char *name, const char *obj, const char *sele, float cutoff,
                     int segi, int quiet);
int ExecutiveGetExtent(PyMOLGlobals * G, const char *name, float *mn, float *mx,
                       int transformed, int state, int weighted);
int ExecutiveGetCameraExtent(PyMOLGlobals * G, const char *name, float *mn, float *mx,
                             int transformed, int state);
int ExecutiveSeleToObject(PyMOLGlobals * G, const char *name, const char *s1, int source, int target,
                          int discrete, int zoom, int quiet, int singletons, int copy_properties=0);
PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals * G, const char *s1, int state,
                                     const char *ref_object, int ref_state);
void ExecutiveInvalidateRep(PyMOLGlobals * G, const char *name, int rep, int level);
void ExecutiveFlag(PyMOLGlobals * G, int flag, const char *s1, int action, int quiet);
void ExecutiveRemoveAtoms(PyMOLGlobals * G, const char *s1, int quiet);
void ExecutiveProtect(PyMOLGlobals * G, const char *s1, int mode, int quiet);
void ExecutiveMask(PyMOLGlobals * G, const char *s1, int mode, int quiet);
void ExecutiveRebuildAll(PyMOLGlobals * G);
void ExecutiveSpheroid(PyMOLGlobals * G, const char *name, int average);
void ExecutiveAddHydrogens(PyMOLGlobals * G, const char *s1, int quiet);
void ExecutiveFixHydrogens(PyMOLGlobals * G, const char *s1, int quiet);
void ExecutiveFuse(PyMOLGlobals * G, const char *s0, const char *s1, int mode, int recolor,
                   int move_flag);
void ExecutiveRenameObjectAtoms(PyMOLGlobals * G, const char *name, int force, int quiet);
int ExecutiveInvert(PyMOLGlobals * G, int quiet);

char *ExecutiveGetNames(PyMOLGlobals * G, int mode, int enabled_only, const char *s0);
bool ExecutiveIsMoleculeOrSelection(PyMOLGlobals * G, const char *name);
int ExecutiveGetType(PyMOLGlobals * G, const char *name, WordType type);
float ExecutiveGetArea(PyMOLGlobals * G, const char *s0, int sta0, int load_b);
void ExecutiveInvalidateSelectionIndicators(PyMOLGlobals *G);
void ExecutiveInvalidateSelectionIndicatorsCGO(PyMOLGlobals *G);
void ExecutiveRenderSelections(PyMOLGlobals * G, int curState, int slot, GridInfo *grid);
void ExecutiveHideSelections(PyMOLGlobals * G);
int ExecutiveSetTitle(PyMOLGlobals * G, const char *name, int state, const char *text);
const char *ExecutiveGetTitle(PyMOLGlobals * G, const char *name, int state);
void ExecutiveSetLastObjectEdited(PyMOLGlobals * G, CObject * o);
CObject *ExecutiveGetLastObjectEdited(PyMOLGlobals * G);
bool ExecutiveIsFullScreen(PyMOLGlobals * G);
void ExecutiveFullScreen(PyMOLGlobals * G, int flag);
PyObject *ExecutiveGetSettingTuple(PyMOLGlobals * G, int index, const char *object, int state);
PyObject *ExecutiveGetSettingText(PyMOLGlobals * G, int index, const char *object, int state);
PyObject *ExecutiveGetSettingOfType(PyMOLGlobals * G, int index, const char *object, int state,
                                    int type);
ObjectMolecule **ExecutiveGetObjectMoleculeVLA(PyMOLGlobals * G, const char *sele);
int ExecutivePairIndices(PyMOLGlobals * G, const char *s1, const char *s2, int state1, int state2,
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

int ExecutiveIsolevel(PyMOLGlobals * G, const char *name, float level, int state, int query,
                      float *result, int quiet);
int ExecutiveTransformObjectSelection(PyMOLGlobals * G, const char *name, int state, const char *s1,
                                      int log, float *matrix, int homogenous, int global);
int ExecutiveTransformSelection(PyMOLGlobals * G, int state, const char *s1, int log,
                                float *ttt, int homogenous);
int ExecutiveTranslateAtom(PyMOLGlobals * G, const char *sele, float *v, int state, int mode,
                           int log);
void ExecutiveSelectRect(PyMOLGlobals * G, BlockRect * rect, int mode);
int ExecutiveMapSetBorder(PyMOLGlobals * G, const char *name, float level, int state);
int ExecutiveMapTrim(PyMOLGlobals * G, const char *name, const char *sele,
                     float buffer, int map_state, int sele_state, int quiet);
int ExecutiveMapDouble(PyMOLGlobals * G, const char *name, int state);
int ExecutiveMapHalve(PyMOLGlobals * G, const char *name, int state, int smooth);

int ExecutiveIdentifyObjects(PyMOLGlobals * G, const char *s1, int mode, int **indexVLA,
                             ObjectMolecule *** objVLA);
int ExecutiveTranslateObjectTTT(PyMOLGlobals * G, const char *name, float *trans, int store, int quiet);
int ExecutiveCombineObjectTTT(PyMOLGlobals * G, const char *name, float *ttt,
                              int reverse_order, int store);
int ExecutiveSetObjectTTT(PyMOLGlobals * G, const char *name, const float *ttt, int state, int quiet, int store);
int ExecutiveGetObjectTTT(PyMOLGlobals * G, const char *name, const float **ttt, int state,
                          int quiet);
int ExecutiveGetObjectMatrix(PyMOLGlobals * G, const char *name, int state, double **matrix,
                             int incl_ttt);
int ExecutiveSetObjectMatrix(PyMOLGlobals * G, const char *name, int state, double *matrix);

int ExecutiveSetGeometry(PyMOLGlobals * G, const char *s1, int geom, int valence);
int ExecutiveSculptIterateAll(PyMOLGlobals * G);
int ExecutiveSmooth(PyMOLGlobals * G, const char *name, int cycles, int window,
                    int first, int last, int ends, int quiet);
int ExecutiveSculptDeactivate(PyMOLGlobals * G, const char *name);
int ExecutiveSculptActivate(PyMOLGlobals * G, const char *name, int state, int match_state,
                            int match_by_segment);
float ExecutiveSculptIterate(PyMOLGlobals * G, const char *name, int state, int n_cycle);
int ExecutiveMapNew(PyMOLGlobals * G, const char *name, int type, float *grid, const char *sele,
                    float buffer, float *minCorner, float *maxCorner,
                    int state, int have_corners, int quiet, int zoom, int normalize,
                    float clamp_floor, float clamp_ceiling, float resolution);

int ***ExecutiveGetBondPrint(PyMOLGlobals * G, const char *name, int max_bond, int max_type,
                             int *dim);
int ExecutiveGetSymmetry(PyMOLGlobals * G, const char *sele, int state, float *a, float *b, float *c,
                        float *alpha, float *beta, float *gamma, char *sgroup,
                        int *defined);
int ExecutiveSetSymmetry(PyMOLGlobals * G, const char *sele, int state, float a, float b, float c,
                        float alpha, float beta, float gamma, const char *sgroup);
int ExecutiveSymmetryCopy(PyMOLGlobals * G,
			  const char *source_name, const char *target_name,
			  int source_mode, int target_mode,
			  int source_state, int target_state,
			  int target_undo, int log, int quiet);
int ExecutiveGetSession(PyMOLGlobals * G, PyObject * dict, const char *names, int partial,
                        int quiet);
int ExecutiveSetSession(PyMOLGlobals * G, PyObject * session, int partial_restore,
                        int quiet);

ObjectMap *ExecutiveFindObjectMapByName(PyMOLGlobals * G, const char *name);

int ExecutiveUnsetSetting(PyMOLGlobals * G, int index, const char *sele,
                          int state, int quiet, int updates);

int ExecutiveAssignSS(PyMOLGlobals * G, const char *target, int state, const char *context,
                      int preserve, ObjectMolecule * single_object, int quiet);

int ExecutiveRampNew(PyMOLGlobals * G, const char *name, const char *src_name, float *range,
                     float *color, int src_state, const char *sele,
                     float beyond, float within, float sigma, int zero, int calc_mode,
                     int quiet);

int ExecutiveValidateObjectPtr(PyMOLGlobals * G, CObject * ptr, int object_type);

int ExecutiveSpectrum(PyMOLGlobals * G, const char *s1, const char *expr, float min, float max,
                      int first, int last, const char *prefix, int digits, int byres, int quiet,
                      float *min_ret, float *max_ret);

int ExecutiveReinitialize(PyMOLGlobals * G, int what, const char *pattern);
const char *ExecutiveFindBestNameMatch(PyMOLGlobals * G, const char *name);
int ExecutiveSetVisFromPyDict(PyMOLGlobals * G, PyObject * dict);
PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals * G);
CField   *ExecutiveGetVolumeField(PyMOLGlobals * G, const char * objName, int state);
int       ExecutiveSetVolumeRamp(PyMOLGlobals * G, const char * objName, float *ramp_list, int list_size);
PyObject *ExecutiveGetVolumeRamp(PyMOLGlobals * G, const char * objName);

float * ExecutiveGetHistogram(PyMOLGlobals * G, const char * objName, int n_points,
        float min_val, float max_val);

int ExecutiveIterateObjectMolecule(PyMOLGlobals * G, ObjectMolecule ** obj,
                                   void **hidden);
int ExecutiveSetObjectColor(PyMOLGlobals * G, const char *name, const char *color, int quiet);
int ExecutiveGetObjectColorIndex(PyMOLGlobals * G, const char *name);
int ExecutiveSetOnOffBySele(PyMOLGlobals * G, const char *name, int onoff);
int ExecutiveSetName(PyMOLGlobals * G, const char *old_name, const char *new_name);
int ExecutiveSetDrag(PyMOLGlobals * G, const char *name, int quiet,int mode);
int ExecutiveGetActiveSeleName(PyMOLGlobals * G, char *name, int create_new, int log);
int ExecutiveGetActiveSele(PyMOLGlobals * G);
int ExecutiveGetActiveAlignmentSele(PyMOLGlobals * G);
CObject *ExecutiveGetExistingCompatible(PyMOLGlobals * G, const char *oname, int type);
int ExecutiveAngle(PyMOLGlobals * G, float *result,
                   const char *nam, const char *s1, const char *s2, const char *s3, int mode,
                   int labels, int reset, int zoom, int quiet, int state);

int ExecutiveDihedral(PyMOLGlobals * G, float *result,
                      const char *nam, const char *s1, const char *s2, const char *s3, const char *s4, int mode,
                      int labels, int reset, int zoom, int quiet, int state);

int ExecutiveMatrixCopy2(PyMOLGlobals * G,
                         CObject * source_obj, CObject * target_obj,
                         int source_mode, int target_mode,
                         int source_state, int target_state,
                         int target_undo, int log, int quiet);

int ExecutiveMatrixCopy(PyMOLGlobals * G,
                        const char *source_name, const char *target_name,
                        int source_mode, int target_mode,
                        int source_state, int target_state,
                        int target_undo, int log, int quiet);

void ExecutiveMemoryDump(PyMOLGlobals * G);
void ExecutiveObjMolSeleOp(PyMOLGlobals * G, int sele, ObjectMoleculeOpRec * op);

int ExecutiveIsomeshEtc(PyMOLGlobals * G,
                        const char *mesh_name, const char *map_name, float lvl,
                        const char *sele, float fbuf, int state,
                        float carve, int map_state, int quiet,
                        int mesh_mode, int box_mode, float alt_lvl);

int ExecutiveIsosurfaceEtc(PyMOLGlobals * G,
                           const char *surf_name, const char *map_name, float lvl,
                           const char *sele, float fbuf, int state,
                           float carve, int map_state, int side,
                           int quiet, int surf_mode, int box_mode);

int ExecutiveVolume(PyMOLGlobals * G, const char *volume_name, const char *map_name,
		    float lvl,
		    const char *sele, float fbuf, int state,
		    float carve, int map_state, int quiet,
		    int mesh_mode, int box_mode, float alt_lvl);

int ExecutiveVolumeColor(PyMOLGlobals * G, const char * volume_name, float * colors, int ncolors );

void ExecutiveMotionDraw(PyMOLGlobals * G, BlockRect *rect, int expected ORTHOCGOARG);

void ExecutiveMotionReinterpolate(PyMOLGlobals * G);

void ExecutiveMotionViewModify(PyMOLGlobals *G, int action, 
                               int index, int count, int target, const char *name, int freeze, int quiet);

void ExecutiveMotionMenuActivate(PyMOLGlobals * G, BlockRect *rect, int expected,int passive, 
                                 int x, int y, int same);
int ExecutiveGroupMotionModify(PyMOLGlobals *G, CObject *group, int action, 
                               int index, int count, int target, int freeze);
int ExecutiveGroupMotion(PyMOLGlobals *G, CObject *group,int action, int first,
                         int last, float power, float bias,
                         int simple, float linear, int wrap,
                         int hand, int window, int cycles, int state, int quiet);
int ExecutiveGroupCombineTTT(PyMOLGlobals *G, CObject *group, float *ttt, int reverse_order, int store);
int ExecutiveGroupTranslateTTT(PyMOLGlobals *G, CObject *group, float *v, int store);
void ExecutiveMotionClick(PyMOLGlobals * G, BlockRect *rect, int mode, int expected, int x, int y, int nearest);
void ExecutiveMotionTrim(PyMOLGlobals * G);
void ExecutiveMotionExtend(PyMOLGlobals * G, int freeze);
int ExecutiveMotionView(PyMOLGlobals *G, int action, int first,
                        int last, float power, float bias,
                        int simple, float linear, const char *name, int wrap,
                        int hand, int window, int cycles,
                        const char *scene_name, float scene_cut, int state, int quiet, int autogen);
int ExecutiveIsomeshEtc(PyMOLGlobals * G,
                        const char *volume_name, const char *map_name, float lvl,
                        const char *sele, float fbuf, int state,
                        float carve, int map_state, int quiet,
                        int mesh_mode, int box_mode, float alt_lvl);
int ExecutiveAssignAtomTypes(PyMOLGlobals * G, const char *s1, int format, int state, int quiet);

PyObject * ExecutiveCEAlign(PyMOLGlobals * G, PyObject * listA, PyObject * listB, int lenA, int lenB,
			    float d0, float d1, int windowSize, int gapMax);

#ifdef _PYMOL_LIB
int *ExecutiveGetRepsInSceneForObject(PyMOLGlobals *G, const char *name);
int *ExecutiveGetRepsForObject(PyMOLGlobals *G, const char *name);
#endif

char *ExecutiveGetObjectNames(PyMOLGlobals * G, int mode, const char *name, int enabled_only, int *numstrs);

CoordSet * ExecutiveGetCoordSet(PyMOLGlobals * G, const char * name, int state, ObjectMolecule ** omp=NULL);

void ExecutiveUndo(PyMOLGlobals * G, int dir);
int ExecutiveSaveUndo(PyMOLGlobals * G, const char *s1, int state);

#endif
