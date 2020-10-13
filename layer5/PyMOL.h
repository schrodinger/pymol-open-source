
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

#ifndef _H_CPyMOL
#define _H_CPyMOL

/* WARNING: DO NOT DEVELOP APPLICATIONS AGAINST THIS PyMOL_* API.  
 * 
 * This is a internal interface to PyMOL for the exclusive use of
 * Schrodinger, LLC in development of PyMOL itself, for
 * development of "wrapped" PyMOL products, and for use as foundation
 * layer for support of stable public APIs.
 *
 * Schrodinger, LLC changes this interface regularly and without
 * notice.  It may even vanish altogether.  Code you develop against
 * this interface is likely to be fragile, time-consuming, and
 * expensive to maintain.
 *
 * Our recommended public APIs for PyMOL development are Python-based:
 * (#1) the old "from pymol import cmd" module-based interface, or
 * (#2) the new "from pymol2 import PyMOL; cmd=PyMOL().cmd"
 * instance-based interface.  If you need a C, C++, Java, or
 * ActiveX/COM interface or widget for molecular visualization, then
 * please contact sales@pymol.org about obtaining access to our
 * independent developer and/or OEM product APIs (e.g. JyMOL).
 * 
 * If you feel that you absolutely must rely upon this PyMOL_* API,
 * then please be sure to create your own lightweight wrapper layer
 * around it so that the rest of your code will be relatively
 * insensitive to changes made here.  Failure to do so will almost
 * certainly result in code that is impossible to maintain over time.
 *
 * You have been warned!
 */

#include "OVreturns.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PYMOL_BUTTON_DOWN           0
#define PYMOL_BUTTON_UP             1

#define PYMOL_KEY_F1         1
#define PYMOL_KEY_F2         2
#define PYMOL_KEY_F3         3
#define PYMOL_KEY_F4         4
#define PYMOL_KEY_F5         5
#define PYMOL_KEY_F6         6
#define PYMOL_KEY_F7         7
#define PYMOL_KEY_F8         8
#define PYMOL_KEY_F9         9
#define PYMOL_KEY_F10        10
#define PYMOL_KEY_F11        11
#define PYMOL_KEY_F12        12
#define PYMOL_KEY_LEFT       100
#define PYMOL_KEY_UP         101
#define PYMOL_KEY_RIGHT      102
#define PYMOL_KEY_DOWN       103
#define PYMOL_KEY_PAGE_UP    104
#define PYMOL_KEY_PAGE_DOWN  105
#define PYMOL_KEY_HOME       106
#define PYMOL_KEY_END        107
#define PYMOL_KEY_INSERT     108

#define PYMOL_BUTTON_LEFT    0
#define PYMOL_BUTTON_MIDDLE  1
#define PYMOL_BUTTON_RIGHT   2
#define PYMOL_BUTTON_SCROLL_FORWARD 3
#define PYMOL_BUTTON_SCROLL_REVERSE 4

#define PYMOL_MODIFIER_SHIFT   1
#define PYMOL_MODIFIER_CTRL    2
#define PYMOL_MODIFIER_ALT     4

#define PYMOL_FALSE    0
#define PYMOL_TRUE     1

#define PYMOL_DEFAULT  -1

#define PYMOL_PROGRESS_SLOW 0
#define PYMOL_PROGRESS_MED  2
#define PYMOL_PROGRESS_FAST 4

#define PYMOL_PROGRESS_SIZE 6
#define PYMOL_RESHAPE_SIZE 5


/* configuration */

#ifndef CPyMOLOptions_DEFINED
typedef struct _CPyMOLOptions CPyMOLOptions;
#define CPyMOLOptions_DEFINED
#endif

CPyMOLOptions *PyMOLOptions_New(void);
void PyMOLOptions_Free(CPyMOLOptions * option);

/* PyMOL instance type */

#ifndef CPyMOL_DEFINED
typedef struct _CPyMOL CPyMOL;
#define CPyMOL_DEFINED
#endif

#define ATOM_PROP_MODEL 0
#define ATOM_PROP_INDEX 1
#define ATOM_PROP_TYPE  2
#define ATOM_PROP_NAME  3
#define ATOM_PROP_RESN  4
#define ATOM_PROP_RESI  5
#define ATOM_PROP_RESV  6
#define ATOM_PROP_CHAIN 7
#define ATOM_PROP_ALT   8
#define ATOM_PROP_SEGI  9
#define ATOM_PROP_ELEM 10
#define ATOM_PROP_SS        11
#define ATOM_PROP_TEXT_TYPE 12
#define ATOM_PROP_CUSTOM    13
#define ATOM_PROP_LABEL     14
#define ATOM_PROP_NUMERIC_TYPE 15
#define ATOM_PROP_Q 16
#define ATOM_PROP_B 17
#define ATOM_PROP_VDW 18
#define ATOM_PROP_ELEC_RADIUS 19
#define ATOM_PROP_PARTIAL_CHARGE 20
#define ATOM_PROP_FORMAL_CHARGE 21
#define ATOM_PROP_STEREO 22
#define ATOM_PROP_CARTOON 23
#define ATOM_PROP_COLOR 24
#define ATOM_PROP_ID 25
#define ATOM_PROP_RANK 26
#define ATOM_PROP_FLAGS 27
#define ATOM_PROP_GEOM  28
#define ATOM_PROP_VALENCE 29
#define ATOM_PROP_X 30
#define ATOM_PROP_Y 31
#define ATOM_PROP_Z 32
#define ATOM_PROP_SETTINGS 33
#define ATOM_PROP_PROPERTIES 34
#define ATOM_PROP_ONELETTER 40

/* return status values */

#define PyMOLstatus_YES               1
#define PyMOLstatus_NO                0

#define PyMOLstatus_SUCCESS           0
#define PyMOLstatus_FAILURE          -1


/* return types */

typedef int PyMOLstatus;

typedef struct {
  PyMOLstatus status;
} PyMOLreturn_status;

typedef struct {
  PyMOLstatus status;
  int value;
} PyMOLreturn_int;

typedef struct {
  PyMOLstatus status;
  float value;
} PyMOLreturn_float;

typedef struct {
  PyMOLstatus status;
  int size;
  float *array;
} PyMOLreturn_float_array;

typedef struct {
  PyMOLstatus status;
  int size;
  int *array;
} PyMOLreturn_int_array;

typedef struct {
  PyMOLstatus status;
  int size;
  char **array;
} PyMOLreturn_string_array;

#define PYMOL_RETURN_VALUE_IS_STRING       0x01
#define PYMOL_RETURN_VALUE_IS_INT          0x02
#define PYMOL_RETURN_VALUE_IS_FLOAT        0x04
#define PYMOL_RETURN_VALUE_IS_FLOAT_ARRAY  0x08

typedef struct {
  PyMOLstatus status;
  short int type;
  char *string;
  int int_value;
  float float_value;
  int array_length;
  float *float_array;
} PyMOLreturn_value;

typedef void PyMOLModalDrawFn(void *G);


/* creation */

CPyMOL *PyMOL_New(void);
CPyMOL *PyMOL_NewWithOptions(const CPyMOLOptions * option);


/* destruction */

void PyMOL_Free(CPyMOL * I);


/* starting and stopping */

void PyMOL_Start(CPyMOL * I);
void PyMOL_ConfigureShadersGL(CPyMOL * I);
void PyMOL_StartWithPython(CPyMOL * I);
void PyMOL_Stop(CPyMOL * I);


/* upstream invalidation and configuration events */

void PyMOL_NeedFakeDrag(CPyMOL * I);
void PyMOL_NeedRedisplay(CPyMOL * I);
void PyMOL_NeedSwap(CPyMOL * I);
void PyMOL_SetClickReady(CPyMOL * I, const char *name, int index, int button, int mod, int x,
                         int y, const float *pos, int state, int bond = -1);
void PyMOL_SetPassive(CPyMOL * I, int onOff);
void PyMOL_NeedReshape(CPyMOL * I, int mode, int x, int y, int width, int height);


/* valid context management */

void PyMOL_PushValidContext(CPyMOL * I);
void PyMOL_PopValidContext(CPyMOL * I);


/* methods requiring a valid OpenGL context*/

void PyMOL_AdaptToHardware(CPyMOL * I);
void PyMOL_Draw(CPyMOL * I);


/* methods that do not require a valid OpenGL context */

void PyMOL_Key(CPyMOL * I, unsigned char k, int x, int y, int modifiers);
void PyMOL_Special(CPyMOL * I, int k, int x, int y, int modifiers);
void PyMOL_Reshape(CPyMOL * I, int width, int height, int force);
void PyMOL_Drag(CPyMOL * I, int x, int y, int modifiers);
void PyMOL_Button(CPyMOL * I, int button, int state, int x, int y, int modifiers);
int PyMOL_Idle(CPyMOL * I);     /* return true if PyMOL is busy doing real
                                   work (not simply racing the CPU) */

void PyMOL_ExpireIfIdle(CPyMOL * I);    /* auto-termination for command-line mode */

typedef void PyMOLSwapBuffersFn(void);

void PyMOL_SetSwapBuffersFn(CPyMOL * I, PyMOLSwapBuffersFn * fn);
void PyMOL_SwapBuffers(CPyMOL * I);     /* only works if above function has been set */
void PyMOL_SetDefaultMouse(CPyMOL * I);
void PyMOL_SetStereoCapable(CPyMOL * I, int stereoCapable);
void PyMOL_InitializeCMol(CPyMOL * I);


/* host query methods */

int PyMOL_GetRedisplay(CPyMOL * I, int reset);
int PyMOL_GetPassive(CPyMOL * I, int reset);
int PyMOL_GetSwap(CPyMOL * I, int reset);
int PyMOL_GetClickReady(CPyMOL * I, int reset);

int PyMOL_GetImageReady(CPyMOL * I, int reset);
PyMOLreturn_int_array PyMOL_GetImageInfo(CPyMOL * I);
int PyMOL_GetImageData(CPyMOL * I, int width, int height, int row_bytes,
                       void *buffer, int mode, int reset);
PyMOLreturn_int_array PyMOL_GetImageDataReturned(CPyMOL * I, int width, int height, int row_bytes,
						 int mode, int reset);

int PyMOL_GetReshape(CPyMOL * I);
int PyMOL_GetIdleAndReady(CPyMOL * I);


/* int array results */

PyMOLreturn_int_array PyMOL_GetReshapeInfo(CPyMOL * I, int reset);

/* string results */

char *PyMOL_GetClickString(CPyMOL * I, int reset);
int PyMOL_FreeResultString(CPyMOL * I, char *st);


/* asynchronous processing (only useful for Python-based multithreaded builds right now) */

int PyMOL_GetBusy(CPyMOL * I, int reset);
void PyMOL_SetBusy(CPyMOL * I, int value);

void PyMOL_ResetProgress(CPyMOL * I);

void PyMOL_SetProgress(CPyMOL * I, int offset, int current, int range);

int PyMOL_GetProgress(CPyMOL * I, int *progress, int reset);
int PyMOL_GetProgressChanged(CPyMOL * I, int reset);

int PyMOL_GetInterrupt(CPyMOL * I, int reset);
void PyMOL_SetInterrupt(CPyMOL * I, int value);


/* modal updates -- PyMOL is busy with some complex task, but we have
   to return control to the host in order to get a valid draw callback */

int PyMOL_GetModalDraw(CPyMOL * I);
void PyMOL_SetModalDraw(CPyMOL * I, PyMOLModalDrawFn * fn);     /* for internal use only */


/* developer/transient privates */

struct PyMOLGlobals* PyMOL_GetGlobals(CPyMOL * I);
struct PyMOLGlobals** PyMOL_GetGlobalsHandle(CPyMOL * I);

void PyMOL_RunTest(CPyMOL * I, int group, int test);


/* for python integration */

void PyMOL_LockAPIAndUnblock(CPyMOL * I);
void PyMOL_BlockAndUnlockAPI(CPyMOL * I);


/* for Jmol integration */

int PyMOL_NewG3DStream(CPyMOL * I, int **array_ptr);
int PyMOL_DelG3DStream(CPyMOL * I, int *array_ptr);


/* Command API */

PyMOLreturn_status PyMOL_CmdBackgroundColor(CPyMOL * I, const char *value);

PyMOLreturn_status PyMOL_CmdReinitialize(CPyMOL * I,
    const char *what,
    const char *object_name);

PyMOLreturn_status PyMOL_CmdLoad(CPyMOL * I,
                                 const char *content,
                                 const char *content_type,
                                 const char *content_format,
                                 const char *object_name,
                                 int state, int discrete, int finish,
                                 int quiet, int multiplex, int zoom);

PyMOLreturn_status PyMOL_CmdLoadRaw(CPyMOL * I,
                                    const char *content,
                                    int content_length,
                                    const char *content_format,
                                    const char *object_name, int state,
                                    int discrete, int finish,
                                    int quiet, int multiplex, int zoom);

PyMOLreturn_status PyMOL_CmdLoadCGO(CPyMOL * I,
                                    const float *content,
                                    int content_length,
                                    const char *object_name, int state, int quiet, int zoom);

PyMOLreturn_status PyMOL_CmdCreate(CPyMOL * I,
                                   const char *name,
                                   const char *selection, int source_state,
                                   int target_state, int discrete,
                                   int zoom, int quiet, int singletons,
                                   const char *extract, int copy_properties);

PyMOLreturn_status PyMOL_CmdZoom(CPyMOL * I, const char *selection, float buffer,
                                 int state, int complete, float animate, int quiet);

PyMOLreturn_status PyMOL_CmdCenter(CPyMOL * I, const char *selection, int state, int origin,
                                   float animate, int quiet);

PyMOLreturn_status PyMOL_CmdOrient(CPyMOL * I, const char *selection, float buffer, int state,
                                   int complete, float animate, int quiet);

PyMOLreturn_status PyMOL_CmdOrigin(CPyMOL * I, const char *selection, int state, int quiet);

PyMOLreturn_status PyMOL_CmdOriginAt(CPyMOL * I, float x, float y, float z, int quiet);

PyMOLreturn_status PyMOL_CmdClip(CPyMOL * I, const char *mode, float amount,
                                 const char *selection,
                                 int state, int quiet);

PyMOLreturn_status PyMOL_CmdShow(CPyMOL * I,
                                 const char *representation,
                                 const char *selection,
                                 int quiet);

PyMOLreturn_status PyMOL_CmdHide(CPyMOL * I,
                                 const char *representation,
                                 const char *selection,
                                 int quiet);

PyMOLreturn_status PyMOL_CmdEnable(CPyMOL * I, const char *name, int quiet);

PyMOLreturn_status PyMOL_CmdDisable(CPyMOL * I, const char *name, int quiet);

PyMOLreturn_status PyMOL_CmdDelete(CPyMOL * I, const char *name, int quiet);

PyMOLreturn_status PyMOL_CmdSet(CPyMOL * I,
                                const char *setting,
                                const char *value,
                                const char *selection,
                                int state, int quiet, int side_effects);

PyMOLreturn_value PyMOL_CmdGet(CPyMOL * I, const char *setting,
			       const char *selection,
			       int state, int quiet);

PyMOLreturn_status PyMOL_CmdUnset(CPyMOL * I, const char *setting, const char *selection,
                                  int state, int quiet, int side_effects);

PyMOLreturn_status PyMOL_CmdSetBond(CPyMOL * I, const char *setting, const char *value,
                                    const char *selection1, const char *selection2,
                                    int state, int quiet, int side_effects);

PyMOLreturn_status PyMOL_CmdUnsetBond(CPyMOL * I, const char *setting,
                                      const char *selection1, const char *selection2,
                                      int state, int quiet, int side_effects);

PyMOLreturn_status PyMOL_CmdColor(CPyMOL * I, const char *color, const char *selection, int flags,
                                  int quiet);

PyMOLreturn_status PyMOL_CmdLabel(CPyMOL * I, const char *selection, const char *text, int quiet);

PyMOLreturn_status PyMOL_CmdSelect(CPyMOL * I, const char *name, const char *selection, int quiet);

PyMOLreturn_status PyMOL_CmdSelectList(CPyMOL * I, const char *name, const char *object, int *list,
                                       int list_len, int state, const char *mode, int quiet);

PyMOLreturn_int PyMOL_CmdGetMovieLength(CPyMOL * I,int quiet);

PyMOLreturn_float PyMOL_CmdGetDistance(CPyMOL * I,
                                       const char *selection1,
                                       const char *selection2, int state, int quiet);

PyMOLreturn_float PyMOL_CmdDistance(CPyMOL * I,
                                    const char *name,
                                    const char *selection1,
                                    const char *selection2,
                                    int mode,
                                    float cutoff,
                                    int label, int reset, int zoom, int state, int quiet);

PyMOLreturn_float PyMOL_CmdGetAngle(CPyMOL * I,
                                    const char *selection1,
                                    const char *selection2,
                                    const char *selection3, int state, int quiet);

PyMOLreturn_float PyMOL_CmdAngle(CPyMOL * I,
                                 const char *name,
                                 const char *selection1,
                                 const char *selection2,
                                 const char *selection3,
                                 int mode,
                                 int label, int reset, int zoom, int state, int quiet);

PyMOLreturn_float PyMOL_CmdGetDihedral(CPyMOL * I,
                                       const char *selection1,
                                       const char *selection2,
                                       const char *selection3,
                                       const char *selection4, int state, int quiet);

PyMOLreturn_float PyMOL_CmdDihedral(CPyMOL * I,
                                    const char *name,
                                    const char *selection1,
                                    const char *selection2,
                                    const char *selection3,
                                    const char *selection4,
                                    int mode,
                                    int label, int reset, int zoom, int state, int quiet);

PyMOLreturn_float_array PyMOL_CmdAlign(CPyMOL * I, const char *source, const char *target,
                                       float cutoff, int cycles, float gap, float extend,
                                       int max_gap, const char *object, const char *matrix,
                                       int source_state, int target_state, int quiet,
                                       int max_skip, int transform, int reset);

PyMOLreturn_status PyMOL_CmdSetView(CPyMOL * I, float *view, int view_len, float animate,
                                    int quiet);
PyMOLreturn_float_array PyMOL_CmdGetView(CPyMOL * I, int quiet);

PyMOLreturn_status PyMOL_CmdDraw(CPyMOL * I, int width, int height,
                                 int antialias, int quiet);

PyMOLreturn_status PyMOL_CmdCapture(CPyMOL * I, int quiet);

PyMOLreturn_status PyMOL_CmdRay(CPyMOL * I, int width, int height, int antialias,
                                float angle, float shift, int renderer, int defer,
                                int quiet);

PyMOLreturn_status PyMOL_CmdIsodot(CPyMOL * I, const char *name, const char *map_name, float level,
                                   const char *selection, float buffer, int state, float carve,
                                   int source_state, int quiet);

PyMOLreturn_status PyMOL_CmdIsomesh(CPyMOL * I, const char *name, const char *map_name, float level,
                                    const char *selection, float buffer, int state, float carve,
                                    int source_state, int quiet);

PyMOLreturn_status PyMOL_CmdIsosurface(CPyMOL * I, const char *name, const char *map_name,
                                       float level, const char *selection, float buffer,
                                       int state, float carve, int source_state, int side,
                                       int mode, int quiet);

PyMOLreturn_status PyMOL_CmdGradient(CPyMOL * I, const char *name, const char *map_name,
                                     float minimum, float maximum, const char *selection,
                                     float buffer, int state, float carve,
                                     int source_state, int quiet);

PyMOLreturn_float PyMOL_CmdIsolevel(CPyMOL * I, const char *name, float level, int state,
                                    int query, int quiet);

PyMOLreturn_status PyMOL_CmdRampNew(CPyMOL * I, const char *name, const char *map, float *range,
                                    int n_range, const char *color, int state, const char *selection,
                                    float beyond, float within, float sigma,
                                    int zero, int calc_mode, int quiet);

PyMOLreturn_status PyMOL_CmdPseudoatom(CPyMOL * I, const char *object_name, const char *sele,
				       const char *name, const char *resn, const char *resi, const char *chain,
				       const char *segi, const char *elem, float vdw, int hetatm,
				       float b, float q, const char *color, const char *label, 
				       int set_xyz, float x, float y, float z,
				       int state, int mode, int quiet);

PyMOLreturn_status PyMOL_CmdTurn(CPyMOL * I, char axis, float angle);

PyMOLreturn_status PyMOL_CmdMPlay(CPyMOL * I, int cmd);

PyMOLreturn_status PyMOL_CmdSetFeedbackMask(CPyMOL * I, int action, int module, int mask);

/* releasing returned values */

int PyMOL_FreeResultArray(CPyMOL * I, void *array);

PyMOLreturn_status PyMOL_CmdRock(CPyMOL * I, int mode);

PyMOLreturn_string_array PyMOL_CmdGetNames(CPyMOL * I, int mode, const char *s0, int enabled_only);

PyMOLreturn_status PyMOL_CmdMapNew(CPyMOL * I, const char *name, int type, float grid_spacing, 
				   const char *selection, int state, int normalize,
				   int zoom, int quiet);

#ifdef _PYMOL_LIB
PyMOLreturn_string_array PyMOL_GetObjectList(CPyMOL * I, const char *s0);
PyMOLreturn_status PyMOL_SetIsEnabledCallback(CPyMOL * I, void *CallbackObject, void (*enabledCallback)(void *, const char *, int ));
PyMOLreturn_int_array PyMOL_GetRepsInSceneForObject(CPyMOL * I, const char *name);
PyMOLreturn_int_array PyMOL_GetRepsForObject(CPyMOL * I, const char *name);
PyMOLreturn_status PyMOL_SetButton(CPyMOL * I, const char *button, const char *modifier, const char *action);
PyMOLreturn_status PyMOL_SetMouseButtonMode(CPyMOL * I, const char *modename);
PyMOLreturn_float_array PyMOL_Spectrum(CPyMOL * I, const char *expression, const char *palette, const char *selection, float minimum, float maximum, int byres, int quiet);

#endif

PyMOLreturn_status PyMOL_ZoomScene(CPyMOL * I, float scale);
PyMOLreturn_status PyMOL_TranslateScene(CPyMOL * I, float x, float y, float z);

PyMOLreturn_value PyMOL_GetVersion(CPyMOL * I);

typedef struct _P_AtomProperty {
  int id;
  short Ptype;
  int offset;
  int maxlen;
} AtomPropertyInfo;

AtomPropertyInfo *PyMOL_GetAtomPropertyInfo(CPyMOL * I, const char *atompropname);

#ifdef __cplusplus
}
#endif
#endif
