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

/* methods requiring a valid OpenGL context*/

void PyMOL_Draw(CPyMOL *I);
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
