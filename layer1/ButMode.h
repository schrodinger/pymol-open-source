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
#define cButModeMoveSlabAndZoom 34

/* simple click action for JyMOL -- no selections involved */

#define cButModeSimpleClick 35

/* new drag actions */

#define cButModeRotDrag 36
#define cButModeMovDrag 37
#define cButModeMovDragZ 38

#define cButModeRotObj 39
#define cButModeMovObj 40
#define cButModeMovObjZ 41
#define cButModeMovFragZ 42
#define cButModeMoveAtomZ 43

#define cButModeCount    44

/* remainder only used in scene... */

#define cButModeScaleSlabShrink 101
#define cButModeScaleSlabExpand 102
#define cButModeMoveSlabForward 103
#define cButModeMoveSlabBackward 104
#define cButModeMoveSlabAndZoomForward 105
#define cButModeMoveSlabAndZoomBackward 106
#define cButModeZoomForward 107
#define cButModeZoomBackward 108

/* conversion */

#define cButModeLeftNone   0
#define cButModeMiddleNone 1
#define cButModeRightNone  2
#define cButModeLeftShft   3
#define cButModeMiddleShft 4
#define cButModeRightShft  5
#define cButModeLeftCtrl   6
#define cButModeMiddleCtrl 7
#define cButModeRightCtrl  8
#define cButModeLeftCtSh   9
#define cButModeMiddleCtSh 10
#define cButModeRightCtSh  11
#define cButModeWheelNone  12
#define cButModeWheelShft  13
#define cButModeWheelCtrl  14
#define cButModeWheelCtSh  15
#define cButModeLeftDouble 16
#define cButModeMiddleDouble 17
#define cButModeRightDouble  18
#define cButModeLeftSingle   19
#define cButModeMiddleSingle 20
#define cButModeRightSingle  21

typedef char CodeType[10];


int ButModeInit(PyMOLGlobals *G);
void ButModeFree(PyMOLGlobals *G);
Block *ButModeGetBlock(PyMOLGlobals *G);
void ButModeSetRate(PyMOLGlobals *G,float renderTime);
void ButModeResetRate(PyMOLGlobals *G);
int ButModeGet(PyMOLGlobals *G,int button);
void ButModeSet(PyMOLGlobals *G,int button,int action);
void ButModeCaption(PyMOLGlobals *G,char *text);
void ButModeCaptionReset(PyMOLGlobals *G);
int ButModeTranslate(PyMOLGlobals *G,int button,int mod);

#endif
