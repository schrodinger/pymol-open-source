
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
#ifndef _H_PyMOLObject
#define _H_PyMOLObject


/* literally a 3-D object...also an object object */

#include"Base.h"
#include"Ray.h"
#include"Rep.h"
#include"Setting.h"
#include"PyMOLGlobals.h"
#include"View.h"
#include"Word.h"

typedef char ObjectNameType[WordLength];

#define cObjectMolecule     1
#define cObjectMap          2
#define cObjectMesh         3

/* cObjectDist now cObjectMeasurement */
#define cObjectMeasurement  4
#define cObjectCallback     5
#define cObjectCGO          6
#define cObjectSurface      7
#define cObjectGadget       8
#define cObjectCalculator   9
#define cObjectSlice        10
#define cObjectAlignment    11
#define cObjectGroup        12
#define cObjectVolume       13

/* 
   the object base class is in the process of being converted to support
   states explicitly (an unfortunate early omission), which will allow
   for simplified implementation of future multi-state objects.
 */

typedef struct CObjectState {
  PyMOLGlobals *G;
  double *Matrix;
} CObjectState;

#ifndef CObject_DEFINED
typedef struct _CObject CObject;
#define CObject_DEFINED
#endif

struct _CObject {
  PyMOLGlobals *G;
  void (*fUpdate) (CObject * I);        /* update representations */
  void (*fRender) (CObject * I, RenderInfo * info);
  void (*fFree) (CObject * I);
  int (*fGetNFrame) (CObject * I);
  void (*fDescribeElement) (CObject * I, int index, char *buffer);
  void (*fInvalidate) (CObject * I, int rep, int level, int state);
  CSetting **(*fGetSettingHandle) (CObject * I, int state);
  char *(*fGetCaption) (CObject * I, char * ch, int len);
  CObjectState *(*fGetObjectState) (CObject * I, int state);
  int type;
  ObjectNameType Name;
  int Color;
  int RepVis[cRepCnt];
  float ExtentMin[3], ExtentMax[3];
  int ExtentFlag, TTTFlag;
  float TTT[16];                /* translate, transform, translate matrix (to apply when rendering) */
  CSetting *Setting;
  int Enabled;                  /* read-only... maintained by Scene */
  int Context;                  /* 0 = Camera, 1 = Unit Window, 2 = Scaled Window */
  CViewElem *ViewElem;          /* for animating objects via the TTT */

  /* not pickled */
  int grid_slot;
  int Grabbed;
};

void ObjectInit(PyMOLGlobals * G, CObject * I);
int ObjectCopyHeader(CObject * I, CObject * src);
void ObjectPurge(CObject * I);
void ObjectSetName(CObject * I, char *name);
void ObjectMakeValidName(char *name);
void ObjectPurgeSettings(CObject * I);
void ObjectFree(CObject * I);
void ObjectUseColor(CObject * I);
void ObjectSetRepVis(CObject * I, int rep, int state);
void ObjectToggleRepVis(CObject * I, int rep);
void ObjectPrepareContext(CObject * I, CRay * ray);
void ObjectSetTTT(CObject * I, float *ttt, int state,int store);
int ObjectGetTTT(CObject * I, float **ttt, int state);
int ObjectGetTotalMatrix(CObject * I, int state, int history, double *matrix);
void ObjectCombineTTT(CObject * I, float *ttt, int reverse_order, int store);
void ObjectTranslateTTT(CObject * T, float *v,int store);
void ObjectSetTTTOrigin(CObject * I, float *origin);
void ObjectResetTTT(CObject * I,int store);
PyObject *ObjectAsPyList(CObject * I);
int ObjectFromPyList(PyMOLGlobals * G, PyObject * list, CObject * I);
int ObjectGetCurrentState(CObject * I, int ignore_all_states);
void ObjectAdjustStateRebuildRange(CObject * I, int *start, int *stop);
int ObjectMotion(CObject * I, int action, int first,
                 int last, float power, float bias,
                 int simple, float linear, int wrap,
                 int hand, int window, int cycles, int state, int quiet);
int ObjectGetSpecLevel(CObject * I, int frame);
void ObjectMotionTrim(CObject *I, int n_frame);
void ObjectMotionExtend(CObject *I, int n_frame);
void ObjectDrawViewElem(CObject *I, BlockRect *rect, int frames);
void ObjectStateInit(PyMOLGlobals * G, CObjectState * I);
void ObjectStateCopy(CObjectState * dst, CObjectState * src);
void ObjectStatePurge(CObjectState * I);
void ObjectStateSetMatrix(CObjectState * I, double *matrix);
double *ObjectStateGetMatrix(CObjectState * I);
void ObjectStateTransformMatrix(CObjectState * I, double *matrix);
void ObjectStateResetMatrix(CObjectState * I);
PyObject *ObjectStateAsPyList(CObjectState * I);
int ObjectStateFromPyList(PyMOLGlobals * G, PyObject * list, CObjectState * I);
int ObjectStatePushAndApplyMatrix(CObjectState * I, RenderInfo * info);
void ObjectStatePopMatrix(CObjectState * I, RenderInfo * info);
void ObjectStateRightCombineMatrixR44d(CObjectState * I, double *matrix);
void ObjectStateLeftCombineMatrixR44d(CObjectState * I, double *matrix);
void ObjectStateCombineMatrixTTT(CObjectState * I, float *matrix);
int ObjectMotionModify(CObject *I,int action, int index, int count,int target, int freeze, int localize);
void ObjectMotionReinterpolate(CObject *I);
int ObjectMotionGetLength(CObject *I);

typedef struct _CObjectUpdateThreadInfo CObjectUpdateThreadInfo;

#endif
