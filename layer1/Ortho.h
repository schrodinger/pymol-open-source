
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

#define cOrthoRightSceneMargin DIP2PIXEL(220)
#define cOrthoBottomSceneMargin DIP2PIXEL(18)
#define cOrthoLineHeight DIP2PIXEL(12)

#include <string>

#include"os_gl.h"
#include"Block.h"
#include"Feedback.h"
#include"Deferred.h"
#include"Image.h"
#include"pymol/memory.h"
#include"PyMOLEnums.h"

#define cOrthoScene 1
#define cOrthoTool 2
#define cOrthoHidden 3

int OrthoInit(PyMOLGlobals * G, int showSplash);
void OrthoFree(PyMOLGlobals * G);

void OrthoAttach(PyMOLGlobals * G, Block * block, int type);
void OrthoDetach(PyMOLGlobals * G, Block * block);

void OrthoReshape(PyMOLGlobals * G, int width, int height, int force);
int OrthoGetWidth(PyMOLGlobals * G);
int OrthoGetHeight(PyMOLGlobals * G);
void OrthoDoDraw(PyMOLGlobals * G, OrthoRenderMode render_mode);
void OrthoDoViewportWhenReleased(PyMOLGlobals *G);
void OrthoPushMatrix(PyMOLGlobals * G);
void OrthoPopMatrix(PyMOLGlobals * G);
int OrthoGetPushed(PyMOLGlobals * G);

int OrthoButton(PyMOLGlobals * G, int button, int state, int x, int y, int mod);
int OrthoButtonDefer(PyMOLGlobals * G, int button, int state, int x, int y, int mod);

void OrthoKey(PyMOLGlobals * G, unsigned char k, int x, int y, int mod);

void OrthoAddOutput(PyMOLGlobals * G, const char *str);
void OrthoNewLine(PyMOLGlobals * G, const char *prompt, int crlf);

int OrthoDrag(PyMOLGlobals * G, int x, int y, int mod);

void OrthoGrab(PyMOLGlobals * G, Block * block);
int OrthoGrabbedBy(PyMOLGlobals * G, Block * block);
void OrthoUngrab(PyMOLGlobals * G);
void OrthoSetLoopRect(PyMOLGlobals * G, int flag, BlockRect * rect);

void OrthoRestorePrompt(PyMOLGlobals * G);
void OrthoBusyDraw(PyMOLGlobals * G, int force);

void OrthoDirty(PyMOLGlobals * G);
int OrthoGetDirty(PyMOLGlobals * G);
void OrthoClear(PyMOLGlobals * G);
void OrthoFakeDrag(PyMOLGlobals * G);
void OrthoBusyMessage(PyMOLGlobals * G, const char *message);
void OrthoBusySlow(PyMOLGlobals * G, int progress, int total);
void OrthoBusyFast(PyMOLGlobals * G, int progress, int total);
void OrthoBusyPrime(PyMOLGlobals * G);
void OrthoCommandSetBusy(PyMOLGlobals * G, int busy);
void OrthoCommandIn(COrtho&, const char *buffer);
inline void OrthoCommandIn(PyMOLGlobals * G, const char *buffer){
  OrthoCommandIn(*G->Ortho, buffer);
}
std::string OrthoCommandOut(COrtho& ortho);
void OrthoCommandNest(PyMOLGlobals * G, int dir);
bool OrthoCommandIsEmpty(COrtho& ortho);

void OrthoFeedbackIn(COrtho& ortho, std::string str);
std::string OrthoFeedbackOut(PyMOLGlobals* G, COrtho& ortho);

void OrthoSetWizardPrompt(PyMOLGlobals * G, char *vla);

int OrthoGetOverlayStatus(PyMOLGlobals * G);
void OrthoPasteIn(PyMOLGlobals * G, const char *buffer);

void OrthoRemoveSplash(PyMOLGlobals * G);
void OrthoRemoveAutoOverlay(PyMOLGlobals * G);

void OrthoSplash(PyMOLGlobals * G);
int OrthoArrowsGrabbed(PyMOLGlobals * G);
void OrthoSpecial(PyMOLGlobals * G, int k, int x, int y, int mod);
int OrthoCommandWaiting(PyMOLGlobals * G);

int OrthoTextVisible(PyMOLGlobals * G);
void OrthoReshapeWizard(PyMOLGlobals * G, ov_size height);
void OrthoDefer(PyMOLGlobals * G, std::unique_ptr<CDeferred> && D);
void OrthoExecDeferred(PyMOLGlobals * G);
int OrthoDeferredWaiting(PyMOLGlobals * G);

OrthoRenderMode OrthoGetRenderMode(PyMOLGlobals * G);
void OrthoDrawBuffer(PyMOLGlobals * G, GLenum mode);
int OrthoGetWrapClickSide(PyMOLGlobals * G);
float *OrthoGetOverlayColor(PyMOLGlobals * G);
void OrthoDrawWizardPrompt(PyMOLGlobals * G, CGO *orthoCGO);

void bg_grad(PyMOLGlobals * G);
GLuint OrthoGetBackgroundTextureID(PyMOLGlobals * G);
void OrthoInvalidateBackgroundTexture(PyMOLGlobals * G);
void OrthoBackgroundTextureNeedsUpdate(PyMOLGlobals * G);

std::pair<int, int> OrthoGetBackgroundSize(const COrtho& ortho);

void OrthoSetBackgroundImage(PyMOLGlobals * G, const char *image_data, int width, int height);

bool OrthoBackgroundDataIsSet(const COrtho& ortho);
std::shared_ptr<pymol::Image> OrthoBackgroundDataGet(const COrtho& ortho);
std::pair<int, int> OrthoGetSize(const COrtho& ortho);

void OrthoInvalidateDoDraw(PyMOLGlobals * G);
void OrthoRenderCGO(PyMOLGlobals * G);

#define OrthoLineLength 1024
typedef char OrthoLineType[OrthoLineLength];

#endif
