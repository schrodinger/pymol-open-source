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
#ifndef _H_Ortho
#define _H_Ortho

#define cOrthoSHIFT 1
#define cOrthoCTRL 2
#define cOrthoALT 4


#define cOrthoRightSceneMargin 220
#define cOrthoBottomSceneMargin 110

#include"Block.h"

#define cOrthoScene 1
#define cOrthoTool 2

void OrthoInit(void);
void OrthoFree(void);

void OrthoAttach(Block *block,int type);
void OrthoDetach(Block *block);

void OrthoSetMargins(Block *block,int t,int l,int b,int r);

Block *OrthoNewBlock(Block *block);
void OrthoFreeBlock(Block *block);

void OrthoReshape(int width,int height);
void OrthoDoDraw(void);

void OrthoPushMatrix(void);
void OrthoPopMatrix(void);

int OrthoButton(int button,int state,int x,int y,int mod);

void OrthoKey(unsigned char k,int x,int y);

void OrthoAddOutput(char *str);
void OrthoNewLine(char *prompt);

int OrthoDrag(int x,int y,int mod);
int OrthoCursor(int x,int y);

void OrthoGrab(Block *block);
void OrthoUngrab(void);

void OrthoRestorePrompt(void);

void OrthoDirty(void);
void OrthoWorking(void);

void OrthoBusyDraw(int force);
void OrthoBusyMessage(char *message);
void OrthoBusySlow(int progress,int total);
void OrthoBusyFast(int progress,int total);
void OrthoBusyPrime(void);

#define OrthoLineLength 1024
typedef char OrthoLineType[OrthoLineLength];

#endif














