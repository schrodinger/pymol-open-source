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

typedef struct {
  Block *Block;
  CodeType Code[cButModeCount+1];
  int NCode;
  int Mode[20];
  int NBut;
  float Rate;
  float Samples;
  WordType Caption;
  float TextColor1[3];
  float TextColor2[3];
  float TextColor3[3];
}  CButMode;

extern CButMode ButMode;

void ButModeInit(void);
void ButModeFree(void);
Block *ButModeGetBlock(void);
void ButModeSetRate(float renderTime);
void ButModeResetRate(void);
void ButModeSet(int button,int action);
void ButModeCaption(char *text);
void ButModeCaptionReset(void);
int ButModeTranslate(int button,int mod);

#endif
