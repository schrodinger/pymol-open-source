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
#define cButModeDragMol  44

#define cButModeRotView 45
#define cButModeMovView 46
#define cButModeMovViewZ 47

#define cButModePotentialClick 48

#define cButModeCount    49

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

#define cButModeLeftDouble   16
#define cButModeMiddleDouble 17
#define cButModeRightDouble  18
#define cButModeLeftSingle   19
#define cButModeMiddleSingle 20
#define cButModeRightSingle  21

#define cButModeLeftShftDouble    22
#define cButModeMiddleShftDouble  23
#define cButModeRightShftDouble   24
#define cButModeLeftShftSingle    25
#define cButModeMiddleShftSingle  26
#define cButModeRightShftSingle   27

#define cButModeLeftCtrlDouble    28
#define cButModeMiddleCtrlDouble  29
#define cButModeRightCtrlDouble   30
#define cButModeLeftCtrlSingle    31
#define cButModeMiddleCtrlSingle  32
#define cButModeRightCtrlSingle   33

#define cButModeLeftCtShDouble    34
#define cButModeMiddleCtShDouble  35
#define cButModeRightCtShDouble   36
#define cButModeLeftCtShSingle    37
#define cButModeMiddleCtShSingle  38
#define cButModeRightCtShSingle   39

#define cButModeLeftAltDouble     40
#define cButModeMiddleAltDouble   41
#define cButModeRightAltDouble    42
#define cButModeLeftAltSingle     43
#define cButModeMiddleAltSingle   44
#define cButModeRightAltSingle    45

#define cButModeLeftAltShftDouble    46
#define cButModeMiddleAltShftDouble  47
#define cButModeRightAltShftDouble   48
#define cButModeLeftAltShftSingle    49
#define cButModeMiddleAltShftSingle  50
#define cButModeRightAltShftSingle   51

#define cButModeLeftCtrlAltDouble    52
#define cButModeMiddleCtrlAltDouble  53
#define cButModeRightCtrlAltDouble   54
#define cButModeLeftCtrlAltSingle    55
#define cButModeMiddleCtrlAltSingle  56
#define cButModeRightCtrlAltSingle   57

#define cButModeLeftCtrlAltShftDouble    58
#define cButModeMiddleCtrlAltShftDouble  59
#define cButModeRightCtrlAltShftDouble   60
#define cButModeLeftCtrlAltShftSingle    61
#define cButModeMiddleCtrlAltShftSingle  62
#define cButModeRightCtrlAltShftSingle   63

#define cButModeWheelAlt   64
#define cButModeWheelAltShft 65
#define cButModeWheelCtrlAlt  66
#define cButModeWheelCtrlAltShft  67

#define cButModeLeftAlt    68
#define cButModeMiddleAlt  69
#define cButModeRightAlt   70

#define cButModeLeftAltShft 71
#define cButModeMiddleAltShft 72
#define cButModeRightAltShft 73

#define cButModeLeftCtrlAlt  74
#define cButModeMiddleCtrlAlt 75
#define cButModeRightCtrlAlt 76

#define cButModeLeftCtrlAltShft 77
#define cButModeMiddleCtrlAltShft 78
#define cButModeRightCtrlAltShft 79

#define cButModeInputCount 80

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
