/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warrn Lyford Delano of DeLano Scientific. 
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

/* WARNING: This is a private interface to PyMOL for the exclusive
 * internal use of DeLano Scientific LLC in development of wrapped
 * PyMOL applications and as an interface layer for supporting
 * public PyMOL APIs.
 *
 * DeLano Scientific LLC will change this interface regularly and
 * without notice.  It may even vanish altogether.  Any and all code
 * you develop against this interface is guaranteed to be fragile,
 * time-consuming, and expenses to maintain.
 * 
 * For these reasons, DO NOT UNDER ANY CIRCUMSTANCE USE THIS API!
 *
 *                       You have been warned.
 *
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

/* configuration */

#ifndef CPyMOLOptions_DEFINED
typedef struct _CPyMOLOptions CPyMOLOptions;
#define CPyMOLOptionsOption_DEFINED
#endif

CPyMOLOptions *PyMOLOptions_New(void);
void PyMOLOptions_Free(CPyMOLOptions *option);

/* PyMOL instance type */

#ifndef CPyMOL_DEFINED
typedef struct _CPyMOL CPyMOL;
#define CPyMOL_DEFINED
#endif

/* creation and destruction */

CPyMOL *PyMOL_New(void);
CPyMOL *PyMOL_NewWithOptions(CPyMOLOptions *option);
void PyMOL_Free(CPyMOL *I);

/* starting and stopping */

void PyMOL_Start(CPyMOL *I);
void PyMOL_Stop(CPyMOL *I);

/* upstream invalidation and configuration events */

void PyMOL_NeedFakeDrag(CPyMOL *I);
void PyMOL_NeedRedisplay(CPyMOL *I);
void PyMOL_NeedSwap(CPyMOL *I);
void PyMOL_SetPassive(CPyMOL *I, int onOff);

/* valid context management */

void PyMOL_PushValidContext(CPyMOL *I);
void PyMOL_PopValidContext(CPyMOL *I);

/* methods requiring a valid OpenGL context*/

void PyMOL_Draw(CPyMOL *I);

/* methods that do not require a valid OpenGL context */

void PyMOL_Key(CPyMOL *I,unsigned char k, int x, int y, int modifiers);
void PyMOL_Special(CPyMOL *I,int k, int x, int y, int modifiers);
void PyMOL_Reshape(CPyMOL *I,int width, int height, int force);
void PyMOL_Drag(CPyMOL *I,int x, int y, int modifiers);
void PyMOL_Button(CPyMOL *I,int button, int state,int x, int y, int modifiers);
int  PyMOL_Idle(CPyMOL *I); /* return true if PyMOL is busy doing real
                               work (not simply racing the CPU) */

typedef void PyMOLSwapBuffersFn(void);

void PyMOL_SetSwapBuffersFn(CPyMOL *I, PyMOLSwapBuffersFn *fn);
void PyMOL_SwapBuffers(CPyMOL *I); /* only works if above function has been set */
void PyMOL_SetDefaultMouse(CPyMOL *I); 

/* host query methods */

int PyMOL_GetRedisplay(CPyMOL *I, int reset);
int PyMOL_GetPassive(CPyMOL *I, int reset);
int PyMOL_GetSwap(CPyMOL *I, int reset);

/* developer/transient privates */

struct _PyMOLGlobals *PyMOL_GetGlobals(CPyMOL *I);
void PyMOL_RunTest(CPyMOL *I, int group, int test);

/* for Jmol integration */

int PyMOL_NewG3DStream(CPyMOL *I,int **array_ptr);
int PyMOL_DelG3DStream(CPyMOL *I,int *array_ptr);

/* Command API */

int PyMOL_Reinitialize(CPyMOL *I);

int PyMOL_Load(CPyMOL *I,char *content, char *content_type, 
               int content_length, char *content_format, 
               char *object_name, 
               int state, int discrete, int finish, 
               int quiet, int multiplex, int zoom);

int PyMOL_Zoom(CPyMOL *I,char *selection, float buffer,
               int state, int complete, float animate, int quiet);

int PyMOL_Center(CPyMOL *I,char *selection, int state, int origin, float animate, int quiet);

int PyMOL_Orient(CPyMOL *I,char *selection, float buffer, int state, int complete, float animate, int quiet);

int PyMOL_Origin(CPyMOL *I,char *selection, int state, int quiet);

int PyMOL_OriginAt(CPyMOL *I,float x, float y, float z, int quiet);

int PyMOL_Clip(CPyMOL *I,char *mode, float amount, char *selection, int state, int quiet);

int PyMOL_Show(CPyMOL *I,char *representation, char *selection,int quiet);

int PyMOL_Hide(CPyMOL *I,char *representation, char *selection,int quiet);

int PyMOL_Delete(CPyMOL *I,char *name, int quiet);

int PyMOL_Set(CPyMOL *I,char *setting, char *value, char *selection, int state, int quiet, int side_effects);

int PyMOL_Color(CPyMOL *I,char *color, char *selection, int flags, int quiet);

int PyMOL_Select(CPyMOL *I,char *name, char *selection, int quiet);

#endif
