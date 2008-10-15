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

#include"PyMOLObject.h"
#include"ObjectMolecule.h"
#include"ObjectMap.h"

#include"Ortho.h"
#include"Word.h"

#define cLoadTypeUnknown -1
#define cLoadTypePDB 0
#define cLoadTypeMOL 1
#define cLoadTypeSDF1 2 /* SDF1 - python-based loader */
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
/* 50 is CIFStr */
#define cLoadTypePHIStr 51
/* 51 is PIM */
/* 52 is PWG */

/* NOTE: if you add new content/object type above, then be sure to add
   corresponding code in: CmdLoad ExecutiveGetExistingCompatible
   ExecutiveLoad
*/


typedef struct {
  ObjectMolecule *obj;
  int offset;
} ExecutiveObjectOffset;

int ExecutiveGetExpandedGroupList(PyMOLGlobals *G,char *name);
int ExecutiveGetExpandedGroupListFromPattern(PyMOLGlobals *G,char *name);
void ExecutiveFreeGroupList(PyMOLGlobals *G,int list_id);

int ExecutiveCheckGroupMembership(PyMOLGlobals *G,int list_id,CObject *obj); /* 0.5*N for group size */


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

int ExecutiveGroup(PyMOLGlobals *G,char *name,char *members,int action, int quiet);

void ExecutiveInvalidateGroups(PyMOLGlobals *G,int force);
void ExecutiveUpdateGroups(PyMOLGlobals *G,int force);

int *ExecutiveGetG3d(PyMOLGlobals *G);
int ExecutiveOrder(PyMOLGlobals *G, char *s1, int sort, int location);
int ExecutiveFixChemistry(PyMOLGlobals *G,char *s1,char *s2,int invalidate, int quiet);
int ExecutiveGetAtomVertex(PyMOLGlobals *G,char *s1,int state,int index,float *v);
int ExecutiveProcessPDBFile(PyMOLGlobals *G,CObject *origObj,char *fname, char *oname,
                             int frame, int discrete,int finish,OrthoLineType buf,
                             PDBInfoRec *pdb_info,int quiet,int is_string,int multiplex,int zoom);

int ExecutiveGetUniqueIDObjectOffsetVLADict(PyMOLGlobals *G, 
                                            ExecutiveObjectOffset **vla, 
                                            OVOneToOne **dict);

#if 0
void ExecutiveLoadMOL2(PyMOLGlobals *G,CObject *origObj,char *fname,
                       char *oname, int frame, int discrete,int finish,
                       OrthoLineType buf,int multiplex,int quiet,
                       int is_string,int zoom);
#endif

int ExecutiveLoad(PyMOLGlobals *G,CObject *origObj, 
                   char *content, int content_length,
                   int content_format,
                   char *object_name, 
                   int state, int zoom, 
                   int discrete, int finish, 
                  int multiplex, int quiet, char *plugin);

int ExecutiveDebug(PyMOLGlobals *G,char *name);

typedef struct {
  int n_residues_aligned;
  float raw_alignment_score;
  int initial_n_atom;
  float initial_rms;
  int n_cycles_run;
  int final_n_atom;
  float final_rms;
} ExecutiveRMSInfo;

int ExecutivePseudoatom(PyMOLGlobals *G, char *object_name, char *sele,
                        char *name, char *resn, char *resi, char *chain,
                        char *segi, char *elem, float vdw, int hetatm,
                        float b, float q, char *label, float *pos, int color, 
                        int state, int mode, int quiet);

int ExecutiveMapSet(PyMOLGlobals *G,char *name,int operator,char *operands,
                    int target_state,int source_state,int zoom, int quiet);

int ExecutiveAlign(PyMOLGlobals *G,char *s1,char *s2,char *mat_file,
                   float gap,float extend,
                   int max_gap, int max_skip,
                   float cutoff,int cycles,int quiet,
                   char *oname,int state1,int state2,
                   ExecutiveRMSInfo *rms_info,int transform,int reset,
                   float seq_wt,float radius,float scale,float base, 
                   float coord_wt, float expect, int window, float ante);

void ExecutiveUpdateColorDepends(PyMOLGlobals *G,ObjectMolecule *mol);
void ExecutiveUpdateCoordDepends(PyMOLGlobals *G,ObjectMolecule *mol);
float ExecutiveDistance(PyMOLGlobals *G,char *sele1,char *sele2);
int  ExecutiveDist(PyMOLGlobals *G,float *result,
                   char *nam,char *s1,char *s2,int mode,
                   float cutoff,int labels,int quiet,
                   int reset,int state,int zoom);
int ExecutiveBond(PyMOLGlobals *G,char *s1,char *s2,int order,int mode,int quiet);
int ExecutiveRevalence(PyMOLGlobals *G,char *s1,char *s2,char *src,
                       int target_state,int source_state, int reset, int quiet);
int ExecutiveVdwFit(PyMOLGlobals *G,char *s1,int state1,char *s2,int state2,float buffer, int quiet);
int ExecutiveIterate(PyMOLGlobals *G,char *s1,char *expr,int read_only,int quiet,PyObject *space);
int ExecutiveIterateList(PyMOLGlobals *G,char *s1,PyObject *list,int read_only,int quiet,PyObject *space);
int ExecutiveSelectList(PyMOLGlobals *G,char *sele_name,char *s1,
                        int *list,int list_len,int state, int mode, int quiet);
int ExecutiveLabel(PyMOLGlobals *G,char *s1,char *expr,int quiet,int eval);
void ExecutiveIterateState(PyMOLGlobals *G,int i1,char *s1,char *expr,int read_only,
                           int atomic_props,int quiet,PyObject *space);
int ExecutiveColor(PyMOLGlobals *G,char *name,char *color,int flags,int quiet);
int ExecutiveInit(PyMOLGlobals *G);
void ExecutiveFree(PyMOLGlobals *G);
int ExecutivePop(PyMOLGlobals *G,char *target,char *source,int quiet);
void ExecutiveManageObject(PyMOLGlobals *G,CObject *obj,int allow_zoom,int quiet);
void ExecutiveUpdateObjectSelection(PyMOLGlobals *G,CObject *obj);
void ExecutiveManageSelection(PyMOLGlobals *G,char *name);
Block *ExecutiveGetBlock(PyMOLGlobals *G);
CObject *ExecutiveFindObjectByName(PyMOLGlobals *G,char *name);
ObjectMolecule *ExecutiveFindObjectMoleculeByName(PyMOLGlobals *G,char *name);
int ExecutiveIterateObject(PyMOLGlobals *G,CObject **obj,void **hidden);
void ExecutiveDelete(PyMOLGlobals *G,char *name);
void ExecutiveDump(PyMOLGlobals *G,char *fname,char *obj);
void ExecutiveSetControlsOff(PyMOLGlobals *G,char *name);
void ExecutiveSort(PyMOLGlobals *G,char *name);
int  ExecutiveSetBondSetting(PyMOLGlobals *G,int index,PyObject *tuple,
                             char *s1,char *s2,
                             int state,int quiet,int updates);
int  ExecutiveUnsetBondSetting(PyMOLGlobals *G,int index,char *s1,char *s2,
                               int state,int quiet,int updates);

int ExecutiveSetSetting(PyMOLGlobals *G,int index,PyObject *tuple,
                        char *sele,int state,
                         int quiet,int updates);
int  ExecutiveSetSettingFromString(PyMOLGlobals *G,int index,char *value,
                                   char *sele,
                                   int state,int quiet,int updates);
int  ExecutiveSetObjSettingFromString(PyMOLGlobals *G,
                                      int index,char *value,CObject *obj,
                                      int state,int quiet,int updates);

int ExecutiveRay(PyMOLGlobals *G,int width,int height,int mode,
                  float angle,float shift,int quiet,int defer, int antialias);
int ExecutiveGetDistance(PyMOLGlobals *G,char *s0,char *s1,float *value,int state);
int ExecutiveGetAngle(PyMOLGlobals *G,char *s0,char *s1,char *s2,float *value,int state);
int ExecutiveGetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float *value,int state);
int ExecutiveSetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float value,int state,int quiet);
int ExecutiveRMS(PyMOLGlobals *G,char *sele1,char *sele2,int mode,float refine,int max_cyc,
                   int quiet,char *oname,int state1,int state2,
                   int ordered_selections, int matchmaker, 
                   ExecutiveRMSInfo *rms_info);

void ExecutiveUpdateCmd(PyMOLGlobals *G,char *sele1,char *sele2,int sta1,int sta2,int method,int quiet);
float ExecutiveRMSPairs(PyMOLGlobals *G,WordType *sele,int pairs,int mode);
float *ExecutiveRMSStates(PyMOLGlobals *G,char *s1,int target,int mode,int quiet, int mix);
int *ExecutiveIdentify(PyMOLGlobals *G,char *s1,int mode);
int ExecutiveIndex(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA);
int ExecutiveReset(PyMOLGlobals *G,int cmd,char *name);
void ExecutiveResetMatrix(PyMOLGlobals *G,
                          char *name,
                          int   mode,
                          int   state,
                          int   log,  
                          int quiet);
void ExecutiveDrawNow(PyMOLGlobals *G);
int ExecutiveDrawCmd(PyMOLGlobals *G, int width, int height,int antialias, int entire_window, int quiet);
int ExecutiveCartoon(PyMOLGlobals *G,int type,char *sele);
void ExecutiveSetAllVisib(PyMOLGlobals *G,int state);
void ExecutiveSetRepVisib(PyMOLGlobals *G,char *name,int rep,int state);
int ExecutiveToggleRepVisib(PyMOLGlobals *G,char *name,int rep);

int ExecutiveSetObjVisib(PyMOLGlobals *G,char *name,int onoff,int parents);

int ExecutiveOrigin(PyMOLGlobals *G,char *name,int preserve,char *oname,float *pos,int state);
int ExecutiveCenter(PyMOLGlobals *G,char *name,int state,int inclusive, float animate, float *pos,int quiet);
void ExecutiveDoZoom(PyMOLGlobals *G,CObject *obj,int is_new, int zoom,int quiet);
int ExecutiveWindowZoom(PyMOLGlobals *G,char *name,float buffer,
                        int state,int inclusive,float animate,int quiet);
int ExecutiveGetMoment(PyMOLGlobals *G,char *name,double *mi,int state);

char *ExecutiveGetChains(PyMOLGlobals *G,char *sele,int state,int *null_chain);

void ExecutiveOrient(PyMOLGlobals *G,char *sele,double *mi,
                     int state,float animate,int complete,float buffer,int quiet);
char *ExecutiveSeleToPDBStr(PyMOLGlobals *G,char *s1,int state,int conectFlag,
                            int mode,char *ref_object,int ref_state,int quiet);
char *ExecutiveNameToSeqAlignStrVLA(PyMOLGlobals *G,char *name,int state,int format,int quiet);

int ExecutiveStereo(PyMOLGlobals *G,int flag);
void ExecutiveCopy(PyMOLGlobals *G,char *src,char *dst,int zoom);
float ExecutiveOverlap(PyMOLGlobals *G,char *s1,int state1,char *s2,int state2,float adjust);
int ExecutiveCountStates(PyMOLGlobals *G,char *s1);
void ExecutiveSymExp(PyMOLGlobals *G,char *name,char *obj,char *sele,float cutoff,int segi,int quiet);
int ExecutiveGetExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state,int weighted);
int ExecutiveGetCameraExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state);
int ExecutiveSeleToObject(PyMOLGlobals *G,char *name,char *s1,
                          int source,int target,int discrete,int zoom,
                          int quiet,int singletons);
PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals *G,char *s1,int state,char *ref_object,int ref_state);
void ExecutiveInvalidateRep(PyMOLGlobals *G,char *name,int rep,int level);
void ExecutiveFlag(PyMOLGlobals *G,int flag,char *s1,int action,int quiet);
void ExecutiveRemoveAtoms(PyMOLGlobals *G,char *s1,int quiet);
void ExecutiveProtect(PyMOLGlobals *G,char *s1,int mode,int quiet);
void ExecutiveMask(PyMOLGlobals *G,char *s1,int mode,int quiet);
void ExecutiveUndo(PyMOLGlobals *G,int dir);
void ExecutiveRebuildAll(PyMOLGlobals *G);
void ExecutiveSpheroid(PyMOLGlobals *G,char *name,int average);
void ExecutiveAddHydrogens(PyMOLGlobals *G,char *s1,int quiet);
void ExecutiveFixHydrogens(PyMOLGlobals *G,char *s1,int quiet);
void ExecutiveFuse(PyMOLGlobals *G,char *s0,char *s1,int mode,int recolor,int move_flag);
void ExecutiveRenameObjectAtoms(PyMOLGlobals *G,char *name,int force);
int ExecutiveInvert(PyMOLGlobals *G,int quiet);
char *ExecutiveGetNames(PyMOLGlobals *G,int mode,int enabled_only,char *s1);
int ExecutiveGetType(PyMOLGlobals *G,char *name,WordType type);
float ExecutiveGetArea(PyMOLGlobals *G,char *s0,int sta0,int load_b);
void ExecutiveRenderSelections(PyMOLGlobals *G,int curState);
void ExecutiveHideSelections(PyMOLGlobals *G);
int ExecutiveSetTitle(PyMOLGlobals *G,char *name,int state,char *text);
char *ExecutiveGetTitle(PyMOLGlobals *G,char *name,int state);
int ExecutiveSaveUndo(PyMOLGlobals *G,char *s1,int state);
void ExecutiveSetLastObjectEdited(PyMOLGlobals *G,CObject *o);
CObject *ExecutiveGetLastObjectEdited(PyMOLGlobals *G);
void ExecutiveFullScreen(PyMOLGlobals *G,int flag);
PyObject *ExecutiveGetSettingTuple(PyMOLGlobals *G,int index,char *object,int state);
PyObject *ExecutiveGetSettingText(PyMOLGlobals *G,int index,char *object,int state);
PyObject *ExecutiveGetSettingOfType(PyMOLGlobals *G,int index,char *object,int state,int type);
ObjectMolecule **ExecutiveGetObjectMoleculeVLA(PyMOLGlobals *G,char *sele);
int ExecutivePairIndices(PyMOLGlobals *G,char *s1,char *s2,int state1,int state2,
                         int mode,float cutoff,float h_angle,
                         int **indexVLA, ObjectMolecule ***objVLA);
void ExecutiveRebuildAllObjectDist(PyMOLGlobals *G);
int ExecutivePhiPsi(PyMOLGlobals *G,char *s1,ObjectMolecule ***objVLA,int **iVLA,
                    float **phiVLA,float **psiVLA,int state) ;
float *ExecutiveGetVertexVLA(PyMOLGlobals *G,char *s1,int state);
int ExecutiveValidName(PyMOLGlobals *G,char *name);
int ExecutiveValidNamePattern(PyMOLGlobals *G,char *name);
int ExecutiveProcessObjectName(PyMOLGlobals *G,char *proposed,char *actual);

int ExecutiveIsolevel(PyMOLGlobals *G,char *name,float level,int state,int query,float *result,int quiet);
int ExecutiveTransformObjectSelection(PyMOLGlobals *G,char *name,int state,
                                      char *s1,int log,float *matrix,
                                      int homogenous,int global);
int ExecutiveTransformSelection(PyMOLGlobals *G,int state,char *s1,int log,float *ttt,int homogenous);
int ExecutiveTranslateAtom(PyMOLGlobals *G,char *sele,float *v,int state,int mode,int log);
void ExecutiveSelectRect(PyMOLGlobals *G,BlockRect *rect,int mode);
int ExecutiveMapSetBorder(PyMOLGlobals *G,char *name,float level,int state);
int ExecutiveMapTrim(PyMOLGlobals *G,char *name,char *sele,
                         float buffer,
                         int map_state,int sele_state,int quiet);
int ExecutiveMapDouble(PyMOLGlobals *G,char *name,int state);
int ExecutiveMapHalve(PyMOLGlobals *G,char *name,int state,int smooth);

int ExecutiveMultiSave(PyMOLGlobals *G,char *fname,char *name,int state,int append);
int ExecutiveIdentifyObjects(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA);
int ExecutiveCombineObjectTTT(PyMOLGlobals *G,char *name,float *ttt,int reverse_order);
int ExecutiveSetObjectTTT(PyMOLGlobals *G,char *name,float *ttt,int state,int quiet);
int ExecutiveGetObjectTTT(PyMOLGlobals *G,char *name,float **ttt,int state,int quiet);
int ExecutiveGetObjectMatrix(PyMOLGlobals *G,char *name,int state,double **matrix,int incl_ttt);
int ExecutiveSetObjectMatrix(PyMOLGlobals *G,char *name,int state,double *matrix);

int ExecutiveSetGeometry(PyMOLGlobals *G,char *s1,int geom,int valence);
int ExecutiveSculptIterateAll(PyMOLGlobals *G);
int ExecutiveSmooth(PyMOLGlobals *G,char *name,int cycles,int window,
                    int first, int last, int ends,int quiet);
int ExecutiveSculptDeactivate(PyMOLGlobals *G,char *name);
int ExecutiveSculptActivate(PyMOLGlobals *G,char *name,int state,int match_state,int match_by_segment);
float ExecutiveSculptIterate(PyMOLGlobals *G,char *name,int state,int n_cycle);
int ExecutiveMapNew(PyMOLGlobals *G,char *name,int type,float *grid,char *sele,
                    float buffer,float *minCorner,float *maxCorner,
                    int state,int have_corners, int quiet,int zoom,int normalize,
                    float clamp_floor, float clamp_ceiling);

int ***ExecutiveGetBondPrint(PyMOLGlobals *G,char *name,int max_bond,int max_type,int *dim);
int ExecutiveSetCrystal(PyMOLGlobals *G,char *sele,float a,float b,float c,
                         float alpha,float beta,float gamma,char *sgroup);
int ExecutiveGetSession(PyMOLGlobals *G,PyObject *dict,char *names,int partial,int quiet);
int ExecutiveSetSession(PyMOLGlobals *G,PyObject *session,int partial_restore,int quiet);

ObjectMap *ExecutiveFindObjectMapByName(PyMOLGlobals *G,char *name);

int  ExecutiveUnsetSetting(PyMOLGlobals *G,int index,char *sele,
                           int state,int quiet,int updates);


int  ExecutiveAssignSS(PyMOLGlobals *G,char *target,int state,char *context,int preserve,int quiet);

int ExecutiveRampNew(PyMOLGlobals *G,char *name,char *src_name,float *range,
                        float *color,int src_state,char *sele,
                        float beyond,float within,float sigma,int zero,int calc_mode,int quiet);

int ExecutiveValidateObjectPtr(PyMOLGlobals *G,CObject *ptr,int object_type);

int ExecutiveSpectrum(PyMOLGlobals *G,char *s1,char *expr,float min,float max,int first,int last,
                      char *prefix,int digits,int byres,int quiet,float *min_ret,float *max_ret);

int ExecutiveReinitialize(PyMOLGlobals *G,int what, char *pattern);
char *ExecutiveFindBestNameMatch(PyMOLGlobals *G,char *name);
int ExecutiveSetVisFromPyDict(PyMOLGlobals *G,PyObject *dict);
PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals *G);
int ExecutiveGetCrystal(PyMOLGlobals *G,char *sele,float *a,float *b,float *c,
                        float *alpha,float *beta,float *gamma,char *sgroup,int *defined);
int ExecutiveIterateObjectMolecule(PyMOLGlobals *G,ObjectMolecule **obj,void **hidden);
int ExecutiveSetObjectColor(PyMOLGlobals *G,char *name,char *color,int quiet);
int ExecutiveGetObjectColorIndex(PyMOLGlobals *G,char *name);
int ExecutiveSetOnOffBySele(PyMOLGlobals *G,char *name,int onoff);
int ExecutiveSetName(PyMOLGlobals *G,char *old_name, char *new_name);
int ExecutiveSetDrag(PyMOLGlobals *G,char *name, int quiet);
int ExecutiveGetActiveSeleName(PyMOLGlobals *G,char *name, int create_new,int log);
int ExecutiveGetActiveSele(PyMOLGlobals *G);
int ExecutiveGetActiveAlignmentSele(PyMOLGlobals *G);
CObject *ExecutiveGetExistingCompatible(PyMOLGlobals *G,char *oname,int type);
int ExecutiveAngle(PyMOLGlobals *G,float *result,
                     char *nam,char *s1,char *s2,char *s3,int mode,
                     int labels,int reset,int zoom,int quiet,int state);

int ExecutiveDihedral(PyMOLGlobals *G,float *result,
                      char *nam,char *s1,char *s2,char *s3,char *s4,int mode,
                      int labels,int reset,int zoom,int quiet,int state);

int ExecutiveMatrixCopy2(PyMOLGlobals *G,
                         CObject *source_obj, CObject *target_obj,
                         int   source_mode,  int target_mode, 
                         int   source_state, int target_state,
                         int   target_undo,
                         int   log,          int quiet);

int ExecutiveMatrixCopy(PyMOLGlobals *G,
                             char *source_name, char *target_name,
                             int source_mode, int target_mode, 
                             int source_state, int target_state,
                             int target_undo,
                             int log, int quiet);

void ExecutiveMemoryDump(PyMOLGlobals *G);
void ExecutiveObjMolSeleOp(PyMOLGlobals *G,int sele,ObjectMoleculeOpRec *op);

int ExecutiveIsomeshEtc(PyMOLGlobals *G, 
                        char *mesh_name, char *map_name, float lvl, 
                        char *sele, float fbuf, int state, 
                        float carve, int map_state, int quiet,
                        int mesh_mode, int box_mode, float alt_lvl);

int ExecutiveIsosurfaceEtc(PyMOLGlobals *G, 
                           char *surf_name, char *map_name, float lvl, 
                           char *sele, float fbuf, int state, 
                           float carve, int map_state, int side,
                           int quiet, int surf_mode, int box_mode);

#endif



