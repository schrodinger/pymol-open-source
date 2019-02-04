

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

#include"os_predef.h"
#include"os_std.h"
#include"MemoryDebug.h"
#include"Feedback.h"
#include"Ortho.h"

int FeedbackInit(PyMOLGlobals * G, int quiet)
{
  int a;

  CFeedback *I;
  I = (G->Feedback = pymol::calloc<CFeedback>(1));

  I->Stack = VLAlloc(char, FB_Total);
  I->Depth = 0;
  G->Feedback->Mask = I->Stack;

  if(quiet) {
    for(a = 0; a < FB_Total; a++) {
      G->Feedback->Mask[a] = 0;
    }
  } else {
    for(a = 0; a < FB_Total; a++) {
      G->Feedback->Mask[a] =
        FB_Output | FB_Results | FB_Errors | FB_Warnings | FB_Actions | FB_Details;
    }

    G->Feedback->Mask[FB_Main] &= ~(FB_Errors); /* suppress opengl errors in main */

  }

  const char *fb_env = getenv("PYMOL_FEEDBACK");
  if(fb_env) {
    int n, sysmod, mask;
    while(sscanf(fb_env, "%i:%i%n", &sysmod, &mask, &n) > 1) {
      FeedbackSetMask(G, sysmod, mask);
      fb_env += n;
    }
  }

  return 1;
}

void FeedbackFree(PyMOLGlobals * G)
{
  CFeedback *I = G->Feedback;

  VLAFreeP(I->Stack);
  FreeP(G->Feedback);

}


/* below we'll presume that any standard feedback on the feedback
module itself will be effected at the Python level, since feedback
levels will be changed as a matter of course inside of PyMOL in order
to quietly perform complex actions.  */

void FeedbackPush(PyMOLGlobals * G)
{
  CFeedback *I = G->Feedback;
  int a;
  I->Depth++;
  VLACheck(I->Stack, char, (I->Depth + 1) * FB_Total);
  G->Feedback->Mask = I->Stack + (I->Depth * FB_Total);
  for(a = 0; a < FB_Total; a++) {
    G->Feedback->Mask[a] = G->Feedback->Mask[a - FB_Total];
  }
  PRINTFD(G, FB_Feedback) " Feedback: push\n" ENDFD;
}

void FeedbackPop(PyMOLGlobals * G)
{
  CFeedback *I = G->Feedback;
  if(I->Depth) {
    I->Depth--;
    G->Feedback->Mask = I->Stack + (I->Depth * FB_Total);
  }
  PRINTFD(G, FB_Feedback) " Feedback: pop\n" ENDFD;
}

void FeedbackSetMask(PyMOLGlobals * G, unsigned int sysmod, unsigned char mask)
{
  int a;
  if((sysmod > 0) && (sysmod < FB_Total)) {
    G->Feedback->Mask[sysmod] = mask;
  } else if(!sysmod) {
    for(a = 0; a < FB_Total; a++) {
      G->Feedback->Mask[a] = mask;
    }
  }
  PRINTFD(G, FB_Feedback)
    " FeedbackSetMask: sysmod %d, mask 0x%02X\n", sysmod, mask ENDFD;
}

void FeedbackDisable(PyMOLGlobals * G, unsigned int sysmod, unsigned char mask)
{
  int a;
  if((sysmod > 0) && (sysmod < FB_Total)) {
    G->Feedback->Mask[sysmod] = G->Feedback->Mask[sysmod] & (0xFF - mask);
  } else if(!sysmod) {
    for(a = 0; a < FB_Total; a++) {
      G->Feedback->Mask[a] = G->Feedback->Mask[a] & (0xFF - mask);
    }
  }
  PRINTFD(G, FB_Feedback)
    " FeedbackDisable: sysmod %d, mask 0x%02X\n", sysmod, mask ENDFD;

}

void FeedbackEnable(PyMOLGlobals * G, unsigned int sysmod, unsigned char mask)
{
  int a;
  if((sysmod > 0) && (sysmod < FB_Total)) {
    G->Feedback->Mask[sysmod] = G->Feedback->Mask[sysmod] | mask;
  } else if(!sysmod) {
    for(a = 0; a < FB_Total; a++) {
      G->Feedback->Mask[a] = G->Feedback->Mask[a] | mask;
    }
  }
  PRINTFD(G, FB_Feedback)
    " FeedbackEnable: sysmod %d, mask 0x%02X\n", sysmod, mask ENDFD;

}

void FeedbackAutoAdd(PyMOLGlobals * G, unsigned int sysmod, unsigned char mask, const char *str)
{
  if(Feedback(G, sysmod, mask))
    FeedbackAddColored(G, str, mask);
}

void FeedbackAdd(PyMOLGlobals * G, const char *str)
{
  OrthoAddOutput(G, str);
}

void FeedbackAddColored(PyMOLGlobals * G, const char *str, unsigned char mask)
{
  FeedbackAdd(G, str);
}
