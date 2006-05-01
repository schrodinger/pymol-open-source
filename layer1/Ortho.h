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
#define cOrthoBottomSceneMargin 16
#define cOrthoLineHeight 12

#include"os_gl.h"
#include"Block.h"
#include"Feedback.h"
#include"Deferred.h"

#define cOrthoScene 1
#define cOrthoTool 2
#define cOrthoHidden 3

int OrthoInit(PyMOLGlobals *G,int showSplash);
void OrthoFree(PyMOLGlobals *G);

void OrthoAttach(PyMOLGlobals *G,Block *block,int type);
void OrthoDetach(PyMOLGlobals *G,Block *block);

void OrthoSetMargins(Block *block,int t,int l,int b,int r);

Block *OrthoNewBlock(PyMOLGlobals *G,Block *block);
void OrthoFreeBlock(PyMOLGlobals *G,Block *block);

void OrthoReshape(PyMOLGlobals *G,int width,int height,int force);
int OrthoGetWidth(PyMOLGlobals *G);
int OrthoGetHeight(PyMOLGlobals *G);
void OrthoDoDraw(PyMOLGlobals *G,int render_mode);

void OrthoPushMatrix(PyMOLGlobals *G);
void OrthoPopMatrix(PyMOLGlobals *G);
int OrthoGetPushed(PyMOLGlobals *G);

int OrthoButton(PyMOLGlobals *G,int button,int state,int x,int y,int mod);

void OrthoKey(PyMOLGlobals *G,unsigned char k,int x,int y,int mod);

void OrthoAddOutput(PyMOLGlobals *G,char *str);
void OrthoNewLine(PyMOLGlobals *G,char *prompt,int crlf);

int OrthoDrag(PyMOLGlobals *G,int x,int y,int mod);

void OrthoGrab(PyMOLGlobals *G,Block *block);
void OrthoUngrab(PyMOLGlobals *G);
void OrthoSetLoopRect(PyMOLGlobals *G,int flag, BlockRect *rect);

void OrthoRestorePrompt(PyMOLGlobals *G);
void OrthoBusyDraw(PyMOLGlobals *G,int force);

void OrthoDirty(PyMOLGlobals *G);
int OrthoGetDirty(PyMOLGlobals *G);
void OrthoWorking(PyMOLGlobals *G);
void OrthoClear(PyMOLGlobals *G);
void OrthoFakeDrag(PyMOLGlobals *G);
void OrthoBusyMessage(PyMOLGlobals *G,char *message);
void OrthoBusySlow(PyMOLGlobals *G,int progress,int total);
void OrthoBusyFast(PyMOLGlobals *G,int progress,int total);
void OrthoBusyPrime(PyMOLGlobals *G);
void OrthoCommandIn(PyMOLGlobals *G,char *buffer);
int  OrthoCommandOut(PyMOLGlobals *G,char *buffer);
void OrthoFeedbackIn(PyMOLGlobals *G,char *buffer);
int OrthoFeedbackOut(PyMOLGlobals *G,char *buffer);

void OrthoSetWizardPrompt(PyMOLGlobals *G,char *vla);

void OrthoPasteIn(PyMOLGlobals *G,char *buffer);
void OrthoRemoveSplash(PyMOLGlobals *G);
void OrthoSplash(PyMOLGlobals *G);
int OrthoArrowsGrabbed(PyMOLGlobals *G);
void OrthoSpecial(PyMOLGlobals *G,int k,int x,int y,int mod);
int  OrthoCommandWaiting(PyMOLGlobals *G);

int OrthoTextVisible(PyMOLGlobals *G);
void OrthoReshapeWizard(PyMOLGlobals *G,int height);
void OrthoDefer(PyMOLGlobals *G,CDeferred *D);
void OrthoExecDeferred(PyMOLGlobals *G);
int OrthoDeferredWaiting(PyMOLGlobals *G);

void OrthoSetLoop(PyMOLGlobals *G,int flag, int l, int r, int t, int b);
int OrthoGetRenderMode(PyMOLGlobals *G);
void OrthoDrawBuffer(PyMOLGlobals *G,GLenum mode);

#define OrthoLineLength 1024
typedef char OrthoLineType[OrthoLineLength];


#endif














