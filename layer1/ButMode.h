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

#define cButModeRotXYZ   0 
#define cButModeTransXY  1
#define cButModeTransZ   2
#define cButModeClipZZ   3
#define cButModeRotZ     4
#define cButModeClipN    5
#define cButModeClipF    6

typedef char CodeType[25];

typedef struct {
  Block *Block;
  CodeType Code[7];
  int Mode[6];
  int NCode;
  float Rate;
  float Samples;
}  CButMode;

extern CButMode ButMode;

void ButModeInit(void);
void ButModeFree(void);
Block *ButModeGetBlock(void);
void ButModeSetRate(float renderTime);
void ButModeResetRate(void);
void ButModeChange(int mode);

#endif
