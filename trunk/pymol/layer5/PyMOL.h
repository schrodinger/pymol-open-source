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
 * internal use of DeLano Scientific LLC in the development of wrapped 
 * PyMOL applications.
 *
 * DeLano Scientific LLC can and will change this interface suddenly
 * and without notice.  It may even vanish altogether.  Any and all
 * code you develop against this interface is guaranteed to be
 * fragile, expensive, and time-consuming to maintain.
 * 
 * DO NOT UNDER ANY CIRCUMSTANCE MAKE USE OF THIS API...
 * 
 * YOU HAVE BEEN WARNED!
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
void PyMOL_SwapBuffers(CPyMOL *I); /* only works if above  function has been set */

/* host query methods */

int PyMOL_GetRedisplay(CPyMOL *I, int reset);
int PyMOL_GetPassive(CPyMOL *I, int reset);
int PyMOL_GetSwap(CPyMOL *I, int reset);

/* developer/transient privates */

struct _PyMOLGlobals *PyMOL_GetGlobals(CPyMOL *I);
void PyMOL_RunTest(CPyMOL *I, int group, int test);

#endif
