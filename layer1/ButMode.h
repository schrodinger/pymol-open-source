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
#ifndef _H_ButMode
#define _H_ButMode

#include"Ortho.h"
#include"Word.h"

#define cButModeRotXYZ   0 
#define cButModeTransXY  1
#define cButModeTransZ   2
#define cButModeClipNF   3
#define cButModeRotZ     4
#define cButModeClipN    5
#define cButModeClipF    6
#define cButModeLB       7
#define cButModeMB       8  
#define cButModeRB       9  
#define cButModeAddToLB 10
#define cButModeAddToMB 11
#define cButModeAddToRB 12
#define cButModePickAtom 13
#define cButModePickBond 14
#define cButModeRotFrag  15
#define cButModeTorFrag  16
#define cButModeMovFrag  17
#define cButModeOrigAt   18
#define cButModeRectAdd  19
#define cButModeRectSub  20
#define cButModeRect     21
#define cButModeNone     22
#define cButModeCent     23
#define cButModePkTorBnd 24
#define cButModeScaleSlab 25
#define cButModeMoveSlab  26
#define cButModePickAtom1 27
#define cButModeMoveAtom 28
#define cButModeMenu     29

#define cButModeSeleSet  30
#define cButModeSeleToggle 31
#define cButModeSeleAdd  32
#define cButModeSeleSub  33


#define cButModeCount    34

/* remaineder only used in scene... */

#define cButModeScaleSlabShrink 101
#define cButModeScaleSlabExpand 102
#define cButModeMoveSlabForward 103
#define cButModeMoveSlabBackward 104

typedef char CodeType[10];


int ButModeInit(PyMOLGlobals *G);
void ButModeFree(PyMOLGlobals *G);
Block *ButModeGetBlock(PyMOLGlobals *G);
void ButModeSetRate(PyMOLGlobals *G,float renderTime);
void ButModeResetRate(PyMOLGlobals *G);
void ButModeSet(PyMOLGlobals *G,int button,int action);
void ButModeCaption(PyMOLGlobals *G,char *text);
void ButModeCaptionReset(PyMOLGlobals *G);
int ButModeTranslate(PyMOLGlobals *G,int button,int mod);

#endif
