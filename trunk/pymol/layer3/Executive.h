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

int ExecutiveOrder(PyMOLGlobals *G, char *s1, int sort, int location);
int ExecutiveFixChemistry(PyMOLGlobals *G,char *s1,char *s2,int invalidate, int quiet);
int ExecutiveGetAtomVertex(PyMOLGlobals *G,char *s1,int state,int index,float *v);
void ExecutiveProcessPDBFile(PyMOLGlobals *G,CObject *origObj,char *fname, char *oname,
                             int frame, int discrete,int finish,OrthoLineType buf,
                             PDBInfoRec *pdb_info,int quiet,int is_string);

void ExecutiveLoadMOL2(PyMOLGlobals *G,CObject *origObj,char *fname,
                       char *oname, int frame, int discrete,int finish,
                       OrthoLineType buf,int multiplex,int quiet,
                       int is_string);

int ExecutiveDebug(PyMOLGlobals *G,char *name);
float ExecutiveAlign(PyMOLGlobals *G,char *s1,char *s2,char *mat_file,float gap,float extend,int skip,
                     float cutoff,int cycles,int quiet,char *oname,int state1,int state2);

float ExecutiveDistance(PyMOLGlobals *G,char *sele1,char *sele2);
float ExecutiveDist(PyMOLGlobals *G,char *nam,char *s1,char *s2,int mode,float cutoff,int labels,int quiet);
void ExecutiveBond(PyMOLGlobals *G,char *s1,char *s2,int order,int add);
int ExecutiveIterate(PyMOLGlobals *G,char *s1,char *expr,int read_only,int quiet);
int ExecutiveIterateList(PyMOLGlobals *G,char *s1,PyObject *list,int read_only,int quiet);
int ExecutiveSelectList(PyMOLGlobals *G,char *sele_name,char *s1,
                        PyObject *list,int quiet,int id_type);
void ExecutiveLabel(PyMOLGlobals *G,char *s1,char *expr,int quiet);
void ExecutiveIterateState(PyMOLGlobals *G,int i1,char *s1,char *expr,int read_only,
                           int atomic_props,int quiet);
int ExecutiveColor(PyMOLGlobals *G,char *name,char *color,int flags,int quiet);
int ExecutiveInit(PyMOLGlobals *G);
void ExecutiveFree(PyMOLGlobals *G);
int ExecutivePop(PyMOLGlobals *G,char *target,char *source,int quiet);
void ExecutiveManageObject(PyMOLGlobals *G,struct CObject *obj,int allow_zoom,int quiet);
void ExecutiveUpdateObjectSelection(PyMOLGlobals *G,struct CObject *obj);
void ExecutiveManageSelection(PyMOLGlobals *G,char *name);
Block *ExecutiveGetBlock(PyMOLGlobals *G);
CObject *ExecutiveFindObjectByName(PyMOLGlobals *G,char *name);
ObjectMolecule *ExecutiveFindObjectMoleculeByName(PyMOLGlobals *G,char *name);
int ExecutiveIterateObject(PyMOLGlobals *G,CObject **obj,void **hidden);
void ExecutiveDelete(PyMOLGlobals *G,char *name);
void ExecutiveDump(PyMOLGlobals *G,char *fname,char *obj);
void ExecutiveSetControlsOff(PyMOLGlobals *G,char *name);
void ExecutiveSort(PyMOLGlobals *G,char *name);
int ExecutiveSetSetting(PyMOLGlobals *G,int index,PyObject *tuple,char *sele,int state,
                         int quiet,int updates);
void ExecutiveRay(PyMOLGlobals *G,int width,int height,int mode,float angle,float shift,int quiet);
int ExecutiveGetDistance(PyMOLGlobals *G,char *s0,char *s1,float *value,int state);
int ExecutiveGetAngle(PyMOLGlobals *G,char *s0,char *s1,char *s2,float *value,int state);
int ExecutiveGetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float *value,int state);
int ExecutiveSetDihe(PyMOLGlobals *G,char *s0,char *s1,char *s2,char *s3,float value,int state,int quiet);
float ExecutiveRMS(PyMOLGlobals *G,char *sele1,char *sele2,int mode,float refine,int max_cyc,
                   int quiet,char *oname,int state1,int state2,
                   int ordered_selections);
void ExecutiveUpdateCmd(PyMOLGlobals *G,char *sele1,char *sele2,int sta1,int sta2);
float ExecutiveRMSPairs(PyMOLGlobals *G,WordType *sele,int pairs,int mode);
float *ExecutiveRMSStates(PyMOLGlobals *G,char *s1,int target,int mode,int quiet);
int *ExecutiveIdentify(PyMOLGlobals *G,char *s1,int mode);
int ExecutiveIndex(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA);
int ExecutiveReset(PyMOLGlobals *G,int cmd,char *name);
void ExecutiveDrawNow(PyMOLGlobals *G);
int ExecutiveCartoon(PyMOLGlobals *G,int type,char *sele);
void ExecutiveSetAllVisib(PyMOLGlobals *G,int state);
void ExecutiveSetRepVisib(PyMOLGlobals *G,char *name,int rep,int state);
int ExecutiveToggleRepVisib(PyMOLGlobals *G,char *name,int rep);

void ExecutiveSetAllRepVisib(PyMOLGlobals *G,char *name,int rep,int state);
void ExecutiveSetObjVisib(PyMOLGlobals *G,char *name,int state);

int ExecutiveOrigin(PyMOLGlobals *G,char *name,int preserve,char *oname,float *pos,int state);
int ExecutiveCenter(PyMOLGlobals *G,char *name,int state,int inclusive);
int ExecutiveWindowZoom(PyMOLGlobals *G,char *name,float buffer,int state,int inclusive);
int ExecutiveGetMoment(PyMOLGlobals *G,char *name,Matrix33d mi,int state);

char *ExecutiveGetChains(PyMOLGlobals *G,char *sele,int state,int *null_chain);

void ExecutiveOrient(PyMOLGlobals *G,char *sele,Matrix33d mi,int state);
char *ExecutiveSeleToPDBStr(PyMOLGlobals *G,char *s1,int state,int conectFlag,int mode);
int ExecutiveStereo(PyMOLGlobals *G,int flag);
void ExecutiveCopy(PyMOLGlobals *G,char *src,char *dst);
float ExecutiveOverlap(PyMOLGlobals *G,char *s1,int state1,char *s2,int state2,float adjust);
int ExecutiveCountStates(PyMOLGlobals *G,char *s1);
void ExecutiveSymExp(PyMOLGlobals *G,char *name,char *obj,char *sele,float cutoff);
int ExecutiveGetExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state,int weighted);
int ExecutiveGetCameraExtent(PyMOLGlobals *G,char *name,float *mn,float *mx,int transformed,int state);
void ExecutiveSeleToObject(PyMOLGlobals *G,char *name,char *s1,int source,int target,int discrete);
PyObject *ExecutiveSeleToChemPyModel(PyMOLGlobals *G,char *s1,int state);
void ExecutiveInvalidateRep(PyMOLGlobals *G,char *name,int rep,int level);
void ExecutiveFlag(PyMOLGlobals *G,int flag,char *s1,int action,int quiet);
void ExecutiveRemoveAtoms(PyMOLGlobals *G,char *s1,int quiet);
void ExecutiveProtect(PyMOLGlobals *G,char *s1,int mode,int quiet);
void ExecutiveMask(PyMOLGlobals *G,char *s1,int mode);
void ExecutiveUndo(PyMOLGlobals *G,int dir);
void ExecutiveRebuildAll(PyMOLGlobals *G);
void ExecutiveSpheroid(PyMOLGlobals *G,char *name,int average);
void ExecutiveAddHydrogens(PyMOLGlobals *G,char *s1,int quiet);
void ExecutiveFuse(PyMOLGlobals *G,char *s0,char *s1,int mode);
void ExecutiveRenameObjectAtoms(PyMOLGlobals *G,char *name,int force);
int ExecutiveInvert(PyMOLGlobals *G,int quiet);
char *ExecutiveGetNames(PyMOLGlobals *G,int mode,int enabled_only);
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
void ExecutiveFocus(PyMOLGlobals *G);
PyObject *ExecutiveGetSettingTuple(PyMOLGlobals *G,int index,char *object,int state);
PyObject *ExecutiveGetSettingText(PyMOLGlobals *G,int index,char *object,int state);
ObjectMolecule **ExecutiveGetObjectMoleculeVLA(PyMOLGlobals *G,char *sele);
int ExecutivePairIndices(PyMOLGlobals *G,char *s1,char *s2,int state1,int state2,
                         int mode,float cutoff,float h_angle,
                         int **indexVLA, ObjectMolecule ***objVLA);
void ExecutiveRebuildAllObjectDist(PyMOLGlobals *G);
int ExecutivePhiPsi(PyMOLGlobals *G,char *s1,ObjectMolecule ***objVLA,int **iVLA,
                    float **phiVLA,float **psiVLA,int state) ;
float *ExecutiveGetVertexVLA(PyMOLGlobals *G,char *s1,int state);
int ExecutiveValidName(PyMOLGlobals *G,char *name);
int ExecutiveIsolevel(PyMOLGlobals *G,char *name,float level,int state);
int ExecutiveTransformObjectSelection(PyMOLGlobals *G,char *name,int state,char *s1,int log,float *ttt);
int ExecutiveTransformSelection(PyMOLGlobals *G,int state,char *s1,int log,float *ttt);
int ExecutiveTranslateAtom(PyMOLGlobals *G,char *sele,float *v,int state,int mode,int log);
void ExecutiveSelectRect(PyMOLGlobals *G,BlockRect *rect,int mode);
int ExecutiveMapSetBorder(PyMOLGlobals *G,char *name,float level);
int ExecutiveMapDouble(PyMOLGlobals *G,char *name,int state);

int ExecutiveMultiSave(PyMOLGlobals *G,char *fname,char *name,int state,int append);
int ExecutiveIdentifyObjects(PyMOLGlobals *G,char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA);
int ExecutiveCombineObjectTTT(PyMOLGlobals *G,char *name,float *ttt);
int ExecutiveSetGeometry(PyMOLGlobals *G,char *s1,int geom,int valence);
int ExecutiveSculptIterateAll(PyMOLGlobals *G);
int ExecutiveSmooth(PyMOLGlobals *G,char *name,int cycles,int window,int first, int last, int ends,int quiet);
int ExecutiveSculptDeactivate(PyMOLGlobals *G,char *name);
int ExecutiveSculptActivate(PyMOLGlobals *G,char *name,int state);
float ExecutiveSculptIterate(PyMOLGlobals *G,char *name,int state,int n_cycle);
int ExecutiveMapNew(PyMOLGlobals *G,char *name,int type,float *grid,char *sele,
                    float buffer,float *minCorner,float *maxCorner,
                    int state,int have_corners, int quiet);

int ***ExecutiveGetBondPrint(PyMOLGlobals *G,char *name,int max_bond,int max_type,int *dim);
int ExecutiveSetCrystal(PyMOLGlobals *G,char *sele,float a,float b,float c,
                         float alpha,float beta,float gamma,char *sgroup);
int ExecutiveGetSession(PyMOLGlobals *G,PyObject *dict);
int ExecutiveSetSession(PyMOLGlobals *G,PyObject *session);

ObjectMap *ExecutiveFindObjectMapByName(PyMOLGlobals *G,char *name);

int  ExecutiveUnsetSetting(PyMOLGlobals *G,int index,char *sele,
                           int state,int quiet,int updates);


int  ExecutiveAssignSS(PyMOLGlobals *G,char *target,int state,char *context,int preserve,int quiet);

int ExecutiveRampMapNew(PyMOLGlobals *G,char *name,char *map_name,PyObject *range,
                        PyObject *color,int map_state,char *sele,
                        float beyond,float within,float sigma,int zero);

int ExecutiveValidateObjectPtr(PyMOLGlobals *G,CObject *ptr,int object_type);

int ExecutiveSpectrum(PyMOLGlobals *G,char *s1,char *expr,float min,float max,int first,int last,
                      char *prefix,int digits,int byres,int quiet,float *min_ret,float *max_ret);

int ExecutiveReinitialize(PyMOLGlobals *G);
char *ExecutiveFindBestNameMatch(PyMOLGlobals *G,char *name);
int ExecutiveSetVisFromPyDict(PyMOLGlobals *G,PyObject *dict);
PyObject *ExecutiveGetVisAsPyDict(PyMOLGlobals *G);
int ExecutiveGetCrystal(PyMOLGlobals *G,char *sele,float *a,float *b,float *c,
                        float *alpha,float *beta,float *gamma,char *sgroup,int *defined);
int ExecutiveIterateObjectMolecule(PyMOLGlobals *G,ObjectMolecule **obj,void **hidden);
int ExecutiveGetObjectColorIndex(PyMOLGlobals *G,char *name);
void ExecutiveToggleAllRepVisib(PyMOLGlobals *G,char *name,int rep);
int ExecutiveSetOnOffBySele(PyMOLGlobals *G,char *name,int onoff);
int ExecutiveSetName(PyMOLGlobals *G,char *old_name, char *new_name);
int ExecutiveGetActiveSeleName(PyMOLGlobals *G,char *name, int create_new);
int ExecutiveGetActiveSele(PyMOLGlobals *G);

#endif



