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
#define cButModePk1      7
#define cButModePk2      8  
#define cButModePk3      9  
#define cButModeAddToPk1 10
#define cButModeAddToPk2 11
#define cButModeAddToPk3 12
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

#define cButModeCount    24

typedef char CodeType[25];

typedef struct {
  Block *Block;
  CodeType Code[cButModeCount+1];
  int NCode;
  int Mode[12];
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
