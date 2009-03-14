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

#ifndef _H_CPyMOL
#define _H_CPyMOL

/* WARNING: DO NOT DEVELOP APPLICATIONS AGAINST THIS PyMOL_* API.  
 * 
 * This is a internal interface to PyMOL for the exclusive use of
 * DeLano Scientific LLC in development of PyMOL itself, for
 * development of "wrapped" PyMOL products, and for use as foundation
 * layer for support of stable public APIs.
 *
 * DeLano Scientific LLC changes this interface regularly and without
 * notice.  It may even vanish altogether.  Code you develop against
 * this interface is likely to be fragile, time-consuming, and
 * expensive to maintain.
 *
 * Our recommended public APIs for PyMOL development are Python-based:
 * (#1) the old "from pymol import cmd" module-based interface, or
 * (#2) the new "from pymol2 import PyMOL; cmd=PyMOL().cmd"
 * instance-based interface.  If you need a C, C++, Java, or
 * ActiveX/COM interface or widget for molecular visualization, then
 * please contact sales@delsci.com about obtaining access to our
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
void PyMOLOptions_Free(CPyMOLOptions *option);

#ifndef PYMOL_NO_PY
CPyMOLOptions *PyMOLOptions_NewWithPython(int argc, char *argv[]);
#endif

/* PyMOL instance type */

#ifndef CPyMOL_DEFINED
typedef struct _CPyMOL CPyMOL;
#define CPyMOL_DEFINED
#endif

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

typedef void PyMOLModalDrawFn(void *G);

/* creation */

CPyMOL *PyMOL_New(void);
CPyMOL *PyMOL_NewWithOptions(CPyMOLOptions *option);

/* destruction */

void PyMOL_Free(CPyMOL *I);

/* starting and stopping */

void PyMOL_Start(CPyMOL *I);
#ifndef PYMOL_NO_PY
void PyMOL_StartWithPython(CPyMOL *I);
#endif
void PyMOL_Stop(CPyMOL *I);

/* upstream invalidation and configuration events */

void PyMOL_NeedFakeDrag(CPyMOL *I); 
void PyMOL_NeedRedisplay(CPyMOL *I); 
void PyMOL_NeedSwap(CPyMOL *I); 
void PyMOL_SetClickReady(CPyMOL *I, char *name, int index, int button, int mod,int x,int y);
void PyMOL_SetPassive(CPyMOL *I, int onOff);
void PyMOL_NeedReshape(CPyMOL *I,int mode, int x, int y, int width, int height);

/* valid context management */

void PyMOL_PushValidContext(CPyMOL *I);
void PyMOL_PopValidContext(CPyMOL *I);

/* methods requiring a valid OpenGL context*/

void PyMOL_AdaptToHardware(CPyMOL *I);
void PyMOL_Draw(CPyMOL *I);

/* methods that do not require a valid OpenGL context */

void PyMOL_Key(CPyMOL *I,unsigned char k, int x, int y, int modifiers);
void PyMOL_Special(CPyMOL *I,int k, int x, int y, int modifiers);
void PyMOL_Reshape(CPyMOL *I,int width, int height, int force);
void PyMOL_Drag(CPyMOL *I,int x, int y, int modifiers);
void PyMOL_Button(CPyMOL *I,int button, int state,int x, int y, int modifiers);
int  PyMOL_Idle(CPyMOL *I); /* return true if PyMOL is busy doing real
                               work (not simply racing the CPU) */

void PyMOL_ExpireIfIdle(CPyMOL *I); /* auto-termination for command-line mode */

typedef void PyMOLSwapBuffersFn(void);

void PyMOL_SetSwapBuffersFn(CPyMOL *I, PyMOLSwapBuffersFn *fn);
void PyMOL_SwapBuffers(CPyMOL *I); /* only works if above function has been set */
void PyMOL_SetDefaultMouse(CPyMOL *I); 

/* host query methods */

int PyMOL_GetRedisplay(CPyMOL *I, int reset);
int PyMOL_GetPassive(CPyMOL *I, int reset);
int PyMOL_GetSwap(CPyMOL *I, int reset);
int PyMOL_GetClickReady(CPyMOL *I, int reset);

int PyMOL_GetImageReady(CPyMOL *I, int reset);
PyMOLreturn_int_array PyMOL_GetImageInfo(CPyMOL *I);
int PyMOL_GetImageData(CPyMOL *I, int width, int height, int row_bytes, 
		       void *buffer, int mode, int reset);

int PyMOL_GetReshape(CPyMOL *I);
int PyMOL_GetIdleAndReady(CPyMOL *I);

/* int array results */

PyMOLreturn_int_array PyMOL_GetReshapeInfo(CPyMOL *I,int reset);
/*PyMOLreturn_int_array PyMOL_GetStereoInfo(CPyMOL *I,int reset);  * to come for MacPyMOL -- blue line, etc */

/* string results */

char *PyMOL_GetClickString(CPyMOL *I,int reset);
int PyMOL_FreeResultString(CPyMOL *I,char *st);

/* asynchronous processing (only useful for Python-based multithreaded builds right now) */

int PyMOL_GetBusy(CPyMOL *I, int reset);
void PyMOL_SetBusy(CPyMOL *I, int value);

void PyMOL_ResetProgress(CPyMOL *I);

void PyMOL_SetProgress(CPyMOL *I,int offset, int current, int range);

int PyMOL_GetProgress(CPyMOL *I,int *progress,int reset);
int PyMOL_GetProgressChanged(CPyMOL *I,int reset);
                       
int PyMOL_GetInterrupt(CPyMOL *I, int reset);
void PyMOL_SetInterrupt(CPyMOL *I, int value);

/* modal updates -- PyMOL is busy with some complex task, but we have
   to return control to the host in order to get a valid draw callback */

int PyMOL_GetModalDraw(CPyMOL *I);
void PyMOL_SetModalDraw(CPyMOL *I, PyMOLModalDrawFn *fn); /* for internal use only */

/* developer/transient privates */

struct _PyMOLGlobals *PyMOL_GetGlobals(CPyMOL *I);
struct _PyMOLGlobals **PyMOL_GetGlobalsHandle(CPyMOL *I);

void PyMOL_RunTest(CPyMOL *I, int group, int test);

/* for python integration */

void PyMOL_LockAPIAndUnblock(CPyMOL *I);
void PyMOL_BlockAndUnlockAPI(CPyMOL *I);

/* for Jmol integration */

int PyMOL_NewG3DStream(CPyMOL *I,int **array_ptr);
int PyMOL_DelG3DStream(CPyMOL *I,int *array_ptr);


/* Command API */

PyMOLreturn_status PyMOL_CmdReinitialize(CPyMOL *I,char *what, char *object_name);

PyMOLreturn_status PyMOL_CmdLoad(CPyMOL *I,char *content, 
                                 char *content_type, 
                                 char *content_format, 
                                 char *object_name, 
                                 int state, int discrete, int finish, 
                                 int quiet, int multiplex, int zoom);

PyMOLreturn_status PyMOL_CmdLoadRaw(CPyMOL *I,char *content, 
                                    int content_length,
                                    char *content_format, 
                                    char *object_name, int state, 
                                    int discrete, int finish, 
                                    int quiet, int multiplex, int zoom);

PyMOLreturn_status PyMOL_CmdLoadCGO(CPyMOL *I,float *content, 
                                    int content_length,
                                    char *object_name, int state, 
                                    int quiet, int zoom);

PyMOLreturn_status PyMOL_CmdCreate(CPyMOL *I, char *name, 
                                   char *selection, int source_state,
                                   int target_state, int discrete, 
                                   int zoom, int quiet, 
                                   int singletons, char *extract);

PyMOLreturn_status PyMOL_CmdZoom(CPyMOL *I,char *selection, float buffer,
               int state, int complete, float animate, int quiet);

PyMOLreturn_status PyMOL_CmdCenter(CPyMOL *I,char *selection, int state, int origin, float animate, int quiet);

PyMOLreturn_status PyMOL_CmdOrient(CPyMOL *I,char *selection, float buffer, int state, int complete, float animate, int quiet);

PyMOLreturn_status PyMOL_CmdOrigin(CPyMOL *I,char *selection, int state, int quiet);

PyMOLreturn_status PyMOL_CmdOriginAt(CPyMOL *I,float x, float y, float z, int quiet);

PyMOLreturn_status PyMOL_CmdClip(CPyMOL *I,char *mode, float amount, char *selection, int state, int quiet);

PyMOLreturn_status PyMOL_CmdShow(CPyMOL *I,char *representation, char *selection,int quiet);

PyMOLreturn_status PyMOL_CmdHide(CPyMOL *I,char *representation, char *selection,int quiet);

PyMOLreturn_status PyMOL_CmdEnable(CPyMOL *I,char *name,int quiet);

PyMOLreturn_status PyMOL_CmdDisable(CPyMOL *I,char *name,int quiet);

PyMOLreturn_status PyMOL_CmdDelete(CPyMOL *I,char *name, int quiet);

PyMOLreturn_status PyMOL_CmdSet(CPyMOL *I,char *setting, char *value, char *selection,
                                int state, int quiet, int side_effects);

PyMOLreturn_status PyMOL_CmdColor(CPyMOL *I,char *color, char *selection, int flags, int quiet);

PyMOLreturn_status PyMOL_CmdLabel(CPyMOL *I,char *selection, char *text, int quiet);

PyMOLreturn_status PyMOL_CmdSelect(CPyMOL *I,char *name, char *selection, int quiet);

PyMOLreturn_status PyMOL_CmdSelectList(CPyMOL *I,char *name, char *object, int *list,
                                       int list_len, int state, char *mode, int quiet);

PyMOLreturn_float PyMOL_CmdGetDistance(CPyMOL *I,
                                       char *selection1,
                                       char *selection2, 
                                       int state, int quiet);

PyMOLreturn_float PyMOL_CmdDistance(CPyMOL *I,
                                    char *name,
                                    char *selection1,
                                    char *selection2, 
                                    int mode,
                                    float cutoff,
                                    int label, 
                                    int reset,
                                    int zoom,
                                    int state,
                                    int quiet);


PyMOLreturn_float PyMOL_CmdGetAngle(CPyMOL *I,
                                    char *selection1,
                                    char *selection2,
                                    char *selection3,
                                    int state,int quiet);

PyMOLreturn_float PyMOL_CmdAngle(CPyMOL *I,
                                 char *name,
                                 char *selection1,
                                 char *selection2, 
                                 char *selection3, 
                                 int mode,
                                 int label, 
                                 int reset,
                                 int zoom,
                                 int state,
                                 int quiet);

PyMOLreturn_float PyMOL_CmdGetDihedral(CPyMOL *I,
                                       char *selection1,
                                       char *selection2,
                                       char *selection3,
                                       char *selection4,
                                       int state, int quiet);

PyMOLreturn_float PyMOL_CmdDihedral(CPyMOL *I,
                                    char *name,
                                    char *selection1,
                                    char *selection2, 
                                    char *selection3, 
                                    char *selection4, 
                                    int mode,
                                    int label, 
                                    int reset,
                                    int zoom,
                                    int state,
                                    int quiet);

PyMOLreturn_float_array PyMOL_CmdAlign(CPyMOL *I, char *source, char *target, float cutoff, 
                                       int cycles, float gap, float extend, int max_gap, 
                                       char *object, char *matrix, int source_state,
                                       int target_state, 
                                       int quiet, int max_skip,
                                       int transform, int reset);

PyMOLreturn_status PyMOL_CmdSetView(CPyMOL *I, float *view, int view_len, float animate, int quiet);
PyMOLreturn_float_array PyMOL_CmdGetView(CPyMOL *I,int quiet);

PyMOLreturn_status PyMOL_CmdDraw(CPyMOL *I,int width, int height,
                                 int antialias, int quiet);

PyMOLreturn_status PyMOL_CmdCapture(CPyMOL *I, int quiet);

PyMOLreturn_status PyMOL_CmdRay(CPyMOL *I,int width, int height,int antialias,
                                float angle, float shift,int renderer, int defer, int quiet);

PyMOLreturn_status PyMOL_CmdIsodot(CPyMOL *I, char *name, char *map_name, float level, char *selection,
                                   float buffer, int state, float carve, int source_state, int quiet);

PyMOLreturn_status PyMOL_CmdIsomesh(CPyMOL *I, char *name, char *map_name, float level, char *selection,
                                    float buffer, int state, float carve, int source_state, int quiet);

PyMOLreturn_status PyMOL_CmdIsosurface(CPyMOL *I,char *name, char *map_name, float level, char *selection,
                                       float buffer, int state, float carve, int source_state, int side, 
                                       int mode, int quiet);

PyMOLreturn_status PyMOL_CmdGradient(CPyMOL *I, char *name, char *map_name, float minimum, float maximum,
                                     char *selection, float buffer, int state, float carve, int source_state, 
                                     int quiet);

PyMOLreturn_float PyMOL_CmdIsolevel(CPyMOL *I,char *name, float level, int state, int query, int quiet);

PyMOLreturn_status PyMOL_CmdRampNew(CPyMOL *I,char *name, char *map, float *range, 
                                    int n_range, char *color, int state, char *selection,
                                    float beyond, float within, float sigma,
                                    int zero, int calc_mode, int quiet);

/* releasing returned values */

int PyMOL_FreeResultArray(CPyMOL *I,void *array);

#endif
