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

int ExecutiveDebug(char *name);
float ExecutiveAlign(char *s1,char *s2,char *mat_file,float gap,float extend,int skip,
                     float cutoff,int cycles,int quiet,char *oname,int state1,int state2);

float ExecutiveDistance(char *sele1,char *sele2);
float ExecutiveDist(char *nam,char *s1,char *s2,int mode,float cutoff,int labels,int quiet);
void ExecutiveBond(char *s1,char *s2,int order,int add);
int ExecutiveIterate(char *s1,char *expr,int read_only);
void ExecutiveLabel(char *s1,char *expr,int quiet);
void ExecutiveIterateState(int i1,char *s1,char *expr,int read_only,int atomic_props);
int ExecutiveColor(char *name,char *color,int flags,int quiet);
void ExecutiveInit(void);
void ExecutiveFree(void);
void ExecutiveManageObject(struct CObject *obj,int allow_zoom,int quiet);
void ExecutiveUpdateObjectSelection(struct CObject *obj);
void ExecutiveManageSelection(char *name);
Block *ExecutiveGetBlock(void);
CObject *ExecutiveFindObjectByName(char *name);
ObjectMolecule *ExecutiveFindObjectMoleculeByName(char *name);
int ExecutiveIterateObject(CObject **obj,void **hidden);
void ExecutiveDelete(char *name);
void ExecutiveDump(char *fname,char *obj);
void ExecutiveSetControlsOff(char *name);
void ExecutiveSort(char *name);
int ExecutiveSetSetting(int index,PyObject *tuple,char *sele,int state,
                         int quiet,int updates);
void ExecutiveRay(int width,int height,int mode,float angle,float shift);
int ExecutiveGetDihe(char *s0,char *s1,char *s2,char *s3,float *value,int state);
int ExecutiveSetDihe(char *s0,char *s1,char *s2,char *s3,float value,int state);
float ExecutiveRMS(char *sele1,char *sele2,int mode,float refine,int max_cyc,
                   int quiet,char *oname,int state1,int state2);
void ExecutiveUpdateCmd(char *sele1,char *sele2,int sta1,int sta2);
float ExecutiveRMSPairs(WordType *sele,int pairs,int mode);
float *ExecutiveRMSStates(char *s1,int target,int mode,int quiet);
int *ExecutiveIdentify(char *s1,int mode);
int ExecutiveIndex(char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA);
int ExecutiveReset(int cmd,char *name);
void ExecutiveDrawNow(void);
int ExecutiveCartoon(int type,char *sele);
void ExecutiveSetAllVisib(int state);
void ExecutiveSetRepVisib(char *name,int rep,int state);
void ExecutiveSetAllRepVisib(char *name,int rep,int state);
void ExecutiveSetObjVisib(char *name,int state);
int ExecutiveOrigin(char *name,int preserve,char *oname,float *pos,int state);
int ExecutiveCenter(char *name,int state,int inclusive);
int ExecutiveWindowZoom(char *name,float buffer,int state,int inclusive);
int ExecutiveGetMoment(char *name,Matrix33d mi,int state);

char *ExecutiveGetChains(char *sele,int state,int *null_chain);

void ExecutiveOrient(char *sele,Matrix33d mi,int state);
char *ExecutiveSeleToPDBStr(char *s1,int state,int conectFlag);
int ExecutiveStereo(int flag);
void ExecutiveCopy(char *src,char *dst);
float ExecutiveOverlap(char *s1,int state1,char *s2,int state2,float adjust);
int ExecutiveCountStates(char *s1);
void ExecutiveSymExp(char *name,char *obj,char *sele,float cutoff);
int ExecutiveGetExtent(char *name,float *mn,float *mx,int transformed,int state,int weighted);
int ExecutiveGetCameraExtent(char *name,float *mn,float *mx,int transformed,int state);
void ExecutiveSeleToObject(char *name,char *s1,int source,int target);
PyObject *ExecutiveSeleToChemPyModel(char *s1,int state);
void ExecutiveInvalidateRep(char *name,int rep,int level);
void ExecutiveFlag(int flag,char *s1,int action,int quiet);
void ExecutiveRemoveAtoms(char *s1);
void ExecutiveProtect(char *s1,int mode);
void ExecutiveMask(char *s1,int mode);
void ExecutiveUndo(int dir);
void ExecutiveRebuildAll(void);
void ExecutiveSpheroid(char *name,int average);
void ExecutiveAddHydrogens(char *s1);
void ExecutiveFuse(char *s0,char *s1,int mode);
void ExecutiveRenameObjectAtoms(char *name,int force);
int ExecutiveInvert(char *s0,char *s1,int mode);
char *ExecutiveGetNames(int mode,int enabled_only);
int ExecutiveGetType(char *name,WordType type);
float ExecutiveGetArea(char *s0,int sta0,int load_b);
void ExecutiveRenderSelections(int curState);
void ExecutiveHideSelections(void);
int ExecutiveSetTitle(char *name,int state,char *text);
char *ExecutiveGetTitle(char *name,int state);
int ExecutiveSaveUndo(char *s1,int state);
void ExecutiveSetLastObjectEdited(CObject *o);
CObject *ExecutiveGetLastObjectEdited(void);
void ExecutiveFullScreen(int flag);
void ExecutiveFocus(void);
PyObject *ExecutiveGetSettingTuple(int index,char *object,int state);
PyObject *ExecutiveGetSettingText(int index,char *object,int state);
int ExecutivePairIndices(char *s1,char *s2,int state1,int state2,
                         int mode,float cutoff,float h_angle,
                         int **indexVLA, ObjectMolecule ***objVLA);
void ExecutiveRebuildAllObjectDist(void);
int ExecutivePhiPsi(char *s1,ObjectMolecule ***objVLA,int **iVLA,
                    float **phiVLA,float **psiVLA,int state) ;
float *ExecutiveGetVertexVLA(char *s1,int state);
int ExecutiveValidName(char *name);
int ExecutiveIsolevel(char *name,float level,int state);
int ExecutiveTransformObjectSelection(char *name,int state,char *s1,int log,float *ttt);
int ExecutiveTransformSelection(int state,char *s1,int log,float *ttt);
int ExecutiveTranslateAtom(char *sele,float *v,int state,int mode,int log);
void ExecutiveSelectRect(BlockRect *rect,int mode);
int ExecutiveMapSetBorder(char *name,float level);
int ExecutiveMultiSave(char *fname,char *name,int state,int append);
int ExecutiveIdentifyObjects(char *s1,int mode,int **indexVLA,ObjectMolecule ***objVLA);
int ExecutiveCombineObjectTTT(char *name,float *ttt);
int ExecutiveSetGeometry(char *s1,int geom,int valence);
int ExecutiveSculptIterateAll(void);
int ExecutiveSmooth(char *name,int cycles,int window,int first, int last, int ends);
int ExecutiveSculptDeactivate(char *name);
int ExecutiveSculptActivate(char *name,int state);
int ExecutiveSculptIterate(char *name,int state,int n_cycle);
int ExecutiveMapNew(char *name,int type,float *grid,char *sele,
                    float buffer,float *minCorner,float *maxCorner,
                    int state);

int ***ExecutiveGetBondPrint(char *name,int max_bond,int max_type,int *dim);
int ExecutiveSetCrystal(char *sele,float a,float b,float c,
                         float alpha,float beta,float gamma,char *sgroup);
int ExecutiveGetSession(PyObject *dict);
int ExecutiveSetSession(PyObject *session);

ObjectMap *ExecutiveFindObjectMapByName(char *name);

int  ExecutiveUnsetSetting(int index,char *sele,
                           int state,int quiet,int updates);

int ExecutiveRampMapNew(char *name,char *map_name,PyObject *range,PyObject *color,int map_state);

int ExecutiveValidateObjectPtr(CObject *ptr,int object_type);

int ExecutiveSpectrum(char *s1,char *expr,float min,float max,int first,int last,
                      char *prefix,int digits,int byres,int quiet,float *min_ret,float *max_ret);

int ExecutiveReinitialize(void);
char *ExecutiveFindBestNameMatch(char *name);
int ExecutiveSetVisFromPyDict(PyObject *dict);
     PyObject *ExecutiveGetVisAsPyDict(void);

#endif



