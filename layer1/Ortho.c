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

#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include"main.h"
#include"Version.h"
#include"MemoryDebug.h"
#include"Err.h"
#include"Util.h"
#include"ListMacros.h"
#include"Ortho.h"
#include"P.h"
#include"Scene.h"
#include"Executive.h"
#include"ButMode.h"
#include "Seq.h"
#include"Control.h"
#include"Setting.h"
#include"Wizard.h"
#include"Queue.h"
#include"Pop.h"
#include"Seq.h"
#include"Text.h"
#include"PyMOLOptions.h"
#include"PyMOL.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#define OrthoSaveLines 0xFF
#define OrthoHistoryLines 0xFF

#define cOrthoCharWidth 8
#define cOrthoLeftMargin 3
#define cOrthoBottomMargin 5

#define WizardMargin1 144
#define WizardMargin2 60

#define ButModeMargin 20
#define ControlMargin 0

#define CMD_QUEUE_MASK 0x3

struct _COrtho {
  Block *Blocks;
  Block *GrabbedBy,*ClickedIn;
  int X,Y,Height,Width;
  int LastX,LastY,LastModifiers;
  int ActiveButton;
  int DrawText;
  int InputFlag; /* whether or not we have active input on the line */

  OrthoLineType Line[OrthoSaveLines+1];
  OrthoLineType History[OrthoHistoryLines+1];
  int HistoryLine,HistoryView;
  int CurLine,CurChar,PromptChar,CursorChar;
  FILE *Pipe;
  char Prompt[255];
  int ShowLines;
  char Saved[OrthoLineLength];
  int SavedPC,SavedCC;
  float TextColor[3],OverlayColor[3],WizardBackColor[3],WizardTextColor[3];
  int DirtyFlag;
  double BusyLast, BusyLastUpdate;
  int BusyStatus[4];
  char BusyMessage[255];
  char *WizardPromptVLA;
  int SplashFlag;
  int HaveSeqViewer;
  BlockRect LoopRect;
  int LoopFlag;
  int cmdNestLevel;
  CQueue *cmdQueue[CMD_QUEUE_MASK+1], *cmdActiveQueue;
  CQueue *feedback;
  int Pushed;
  CDeferred *deferred;
  int RenderMode;
  GLint ViewPort[4];
  int WrapXFlag;
  GLenum ActiveGLBuffer;
  double DrawTime, LastDraw;
  int WrapClickSide; /* ugly kludge for finding click side in geowall stereo mode */
};


void OrthoParseCurrentLine(PyMOLGlobals *G);
static void OrthoDrawWizardPrompt(PyMOLGlobals *G);

Block *OrthoFindBlock(PyMOLGlobals *G,int x,int y);
void OrthoKeyControl(PyMOLGlobals *G,unsigned char k);
void OrthoKeyAlt(PyMOLGlobals *G,unsigned char k);

#define cBusyWidth 240
#define cBusyHeight 60
#define cBusyMargin 10
#define cBusyBar 10
#define cBusySpacing 15

#define cBusyUpdate 0.2

#define cWizardTopMargin 15
#define cWizardLeftMargin 15
#define cWizardBorder 7

static int get_wrap_x(int x, int *last_x, int width, int *click_side)
{
  int width_2 = width/2;
  int width_3 = width/3;
  if(!last_x) {
    if(x>width_2) {
      x-=width_2;
      if(click_side) *click_side = 1;
    } else {
      if(click_side) *click_side = -1;
    }
  } else {
    if((x-(*last_x))>width_3) {
      x-=width_2;
      if(click_side) *click_side = 1;
    } else if(((*last_x)-x)>width_3) {
      x+=width_2;
      if(click_side) *click_side = 1;
    } else {
      if(click_side) *click_side = -1;
    }
  }
  return x;
}

void OrthoDrawBuffer(PyMOLGlobals *G,GLenum mode)
{
  register COrtho *I=G->Ortho;
  if((mode!=I->ActiveGLBuffer) && G->HaveGUI && G->ValidContext) {  
    glDrawBuffer(mode);
    I->ActiveGLBuffer = mode;
  }
}

int OrthoGetDirty(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  return I->DirtyFlag;
}
int OrthoGetRenderMode(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  return I->RenderMode;
}

void OrthoSetLoopRect(PyMOLGlobals *G,int flag, BlockRect *rect)
{
  register COrtho *I=G->Ortho;
  I->LoopRect = (*rect);
  I->LoopFlag=flag;
  OrthoDirty(G);
}

int OrthoDeferredWaiting(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  return (I->deferred!=NULL);
}

void OrthoExecDeferred(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  /* execute all deferred actions that happened to require a
   * valid OpenGL context (such as atom picks, etc.) */

  I->deferred = DeferredExec(I->deferred); 
}

void OrthoDefer(PyMOLGlobals *G,CDeferred *D)
{
  register COrtho *I=G->Ortho;
  register CDeferred *d = I->deferred;
  if(d) {
    while(d->next)
      d = d->next;
    d->next = D; 
  } else {
    I->deferred = D;
  }
  OrthoDirty(G);
}

int OrthoGetWidth(PyMOLGlobals *G)
{
  if(G) {
    register COrtho *I=G->Ortho;
    return(I->Width);
  }
  return 0;
}

int OrthoGetHeight(PyMOLGlobals *G)
{
  if(G) {
    register COrtho *I=G->Ortho;
    return(I->Height);
  }
  return 0;
}
/*========================================================================*/
void OrthoFakeDrag(PyMOLGlobals *G) /* for timing-based events, such as pop-ups */
{
  register COrtho *I=G->Ortho;  
  if(I->GrabbedBy)
    OrthoDrag(G,I->LastX,I->LastY,I->LastModifiers);
}
/*========================================================================*/


void OrthoSetWizardPrompt(PyMOLGlobals *G,char *vla)
{
  register COrtho *I=G->Ortho;
  VLAFreeP(I->WizardPromptVLA);
  I->WizardPromptVLA=vla;
}
/*========================================================================*/
void OrthoSpecial(PyMOLGlobals *G,int k,int x,int y,int mod)
{
  register COrtho *I=G->Ortho;
  int curLine = I->CurLine&OrthoSaveLines;
  switch(k) {
  case P_GLUT_KEY_DOWN:
    if(I->CurChar&&(I->HistoryView==I->HistoryLine)) {
      strcpy(I->History[I->HistoryLine],I->Line[curLine]+I->PromptChar);
    }
    I->HistoryView = (I->HistoryView+1)&OrthoHistoryLines;
    strcpy(I->Line[curLine],I->Prompt);
    I->PromptChar = strlen(I->Prompt);
    if(I->History[I->HistoryView][0]) {
      strcat(I->Line[curLine],I->History[I->HistoryView]);
      I->CurChar = strlen(I->Line[curLine]);
    } else {
      I->CurChar=I->PromptChar;
    }
    I->InputFlag=1;
    I->CursorChar=-1;
    break;
  case P_GLUT_KEY_UP:
    if(I->CurChar&&(I->HistoryView==I->HistoryLine)) {
      strcpy(I->History[I->HistoryLine],I->Line[curLine]+I->PromptChar);
    }
    I->HistoryView = (I->HistoryView-1)&OrthoHistoryLines;
    strcpy(I->Line[curLine],I->Prompt);
    I->PromptChar = strlen(I->Prompt);
    if(I->History[I->HistoryView][0]) {
      strcat(I->Line[curLine],I->History[I->HistoryView]);
      I->CurChar = strlen(I->Line[curLine]);
    } else {
      I->CurChar=I->PromptChar;
    }
    I->CursorChar=-1;
    I->InputFlag=1;
    break;
  case P_GLUT_KEY_LEFT:
    if(I->CursorChar>=0) {
      I->CursorChar--;
    } else {
      I->CursorChar = I->CurChar-1;
    }
    if(I->CursorChar<I->PromptChar)
      I->CursorChar=I->PromptChar;
    break;
  case P_GLUT_KEY_RIGHT:
    if(I->CursorChar>=0) {
      I->CursorChar++;
    } else {
      I->CursorChar = I->CurChar-1;
    }
    if((unsigned)I->CursorChar>strlen(I->Line[curLine]))
      I->CursorChar=strlen(I->Line[curLine]);
    break;
  }
  OrthoDirty(G);
}
/*========================================================================*/
int OrthoTextVisible(PyMOLGlobals *G) {
  return(SettingGet(G,cSetting_internal_feedback)||
         SettingGet(G,cSetting_text)||
         SettingGet(G,cSetting_overlay));
}
/*========================================================================*/

int OrthoArrowsGrabbed(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  return((I->CurChar>I->PromptChar)&&OrthoTextVisible(G)); 
  /* arrows can't be grabbed if text isn't visible */
}
/*========================================================================*/
void  OrthoRemoveSplash(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  I->SplashFlag=false;
}
/*========================================================================*/
void OrthoCommandNest(PyMOLGlobals *G,int dir)
{
  register COrtho *I=G->Ortho;
  I->cmdNestLevel += dir;
  {
    int level = I->cmdNestLevel;
    if(level<0)
      level = 0;
    if(level>CMD_QUEUE_MASK)
      level = CMD_QUEUE_MASK;
    I->cmdActiveQueue = I->cmdQueue[level];
  }
}
/*========================================================================*/
int  OrthoCommandOut(PyMOLGlobals *G,char *buffer)
{
  if(G&&buffer) {
    register COrtho *I=G->Ortho;
    
    if(I && I->cmdActiveQueue) {
      int result;
      result = QueueStrOut(I->cmdActiveQueue,buffer);
      return(result);
    } else
      return 0;
  } else
    return 0;
}
/*========================================================================*/
int  OrthoCommandWaiting(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  return(QueueStrCheck(I->cmdActiveQueue));
}
/*========================================================================*/
void OrthoClear(PyMOLGlobals *G)
{
  int a;
  register COrtho *I=G->Ortho;
  for(a=0;a<=OrthoSaveLines;a++)
    I->Line[a][0]=0;
  OrthoNewLine(G,NULL,true);
  OrthoRestorePrompt(G);
  OrthoDirty(G);
}
/*========================================================================*/
void OrthoFeedbackIn(PyMOLGlobals *G,char *buffer)
{
  register COrtho *I=G->Ortho;
  if(G->HaveGUI) {
    if(I->feedback)
      QueueStrIn(I->feedback,buffer);
  }
}
/*========================================================================*/
int OrthoFeedbackOut(PyMOLGlobals *G,char *buffer)
{
  register COrtho *I=G->Ortho;
  if(I->feedback)
	return(QueueStrOut(I->feedback,buffer));
  else
	return(0);
}
/*========================================================================*/
void OrthoDirty(PyMOLGlobals *G) {
  register COrtho *I=G->Ortho;
  PRINTFD(G,FB_Ortho)
    " OrthoDirty: called.\n"
    ENDFD;
  if(!I->DirtyFlag) {
    I->DirtyFlag = true;
  }

  PyMOL_NeedRedisplay(G->PyMOL);
}
/*========================================================================*/
void OrthoBusyMessage(PyMOLGlobals *G,char *message)
{
  register COrtho *I=G->Ortho;
  if(strlen(message)<255)
    strcpy(I->BusyMessage,message);
}
/*========================================================================*/
void OrthoBusySlow(PyMOLGlobals *G,int progress,int total)
{
  register COrtho *I=G->Ortho;
  double   time_yet = (-I->BusyLastUpdate) + UtilGetSeconds(G);

  PRINTFD(G,FB_Ortho)
    " OrthoBusySlow-DEBUG: progress %d total %d\n",progress,total
    ENDFD;
  I->BusyStatus[0]=progress;
  I->BusyStatus[1]=total;
  if(SettingGetGlobal_b(G,cSetting_show_progress)&&(time_yet>0.15F)) {
    if(PyMOL_GetBusy(G->PyMOL,false)) { /* harmless race condition */
#ifndef _PYMOL_NOPY
      int blocked = PAutoBlock(G);
      if(PLockStatusAttempt(G)) {
#endif
        PyMOL_SetProgress(G->PyMOL,PYMOL_PROGRESS_SLOW,progress,total);
        I->BusyLastUpdate = UtilGetSeconds(G);

#ifndef _PYMOL_NOPY
        PUnlockStatus(G);
      }
      PAutoUnblock(G,blocked);
#endif
    }
    OrthoBusyDraw(G,false);
  }
}
/*========================================================================*/
void OrthoBusyFast(PyMOLGlobals *G,int progress,int total)
{
  register COrtho *I=G->Ortho;
  double   time_yet = (-I->BusyLastUpdate) + UtilGetSeconds(G);
  PRINTFD(G,FB_Ortho)
    " OrthoBusyFast-DEBUG: progress %d total %d\n",progress,total
    ENDFD;
  I->BusyStatus[2]=progress;
  I->BusyStatus[3]=total;
  if(SettingGetGlobal_b(G,cSetting_show_progress)&&(time_yet>0.15F)) {
    if(PyMOL_GetBusy(G->PyMOL,false)) { /* harmless race condition */
#ifndef _PYMOL_NOPY
      int blocked = PAutoBlock(G);
      if(PLockStatusAttempt(G)) {
#endif
        PyMOL_SetProgress(G->PyMOL,PYMOL_PROGRESS_FAST,progress,total);
        I->BusyLastUpdate = UtilGetSeconds(G);
#ifndef _PYMOL_NOPY
        PUnlockStatus(G);
      }
      PAutoUnblock(G,blocked);
#endif
    }
    OrthoBusyDraw(G,false);
  }
}
/*========================================================================*/
void OrthoBusyPrime(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  int a;
  for(a=0;a<4;a++)
    I->BusyStatus[a]=0;
  I->BusyMessage[0]=0;
  I->BusyLast = UtilGetSeconds(G);
  I->BusyLastUpdate = UtilGetSeconds(G);
}

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _MACPYMOL_XCODE
void MacPyMOLSetProgress(float value);
#endif
/* END PROPRIETARY CODE SEGMENT */

/*========================================================================*/
void OrthoBusyDraw(PyMOLGlobals *G,int force)
{
  register COrtho *I=G->Ortho;
  double now;
  double busyTime;
  
  PRINTFD(G,FB_Ortho)
    " OrthoBusyDraw: entered.\n"
    ENDFD;
  now = UtilGetSeconds(G);
  busyTime = (-I->BusyLast) + now;
  if(SettingGet(G,cSetting_show_progress)&&(force||(busyTime>cBusyUpdate))) {
    
    I->BusyLast=now;
    if(PIsGlutThread()) {

#ifdef _MACPYMOL_XCODE
      /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
      float busyValue;
      if(I->BusyStatus[1]) {
        busyValue=(I->BusyStatus[0]*1.0F/I->BusyStatus[1]);
      }
      if(I->BusyStatus[3]) {
        busyValue=(I->BusyStatus[2]*1.0F/I->BusyStatus[3]);
      }
      MacPyMOLSetProgress(busyValue);
      /* END PROPRIETARY CODE SEGMENT */
#else
      if(G->HaveGUI && G->ValidContext) {      
        char *c;
        int x,y;
        float black[3] = {0,0,0};
        float white[3] = {1,1,1};
		int draw_both = SceneMustDrawBoth(G);
        OrthoPushMatrix(G);

        {
          int pass = 0;
          glClear(GL_DEPTH_BUFFER_BIT);
          while(1) {
            if(draw_both) {
              if(!pass) 
                OrthoDrawBuffer(G,GL_FRONT_LEFT); 
              else
                OrthoDrawBuffer(G,GL_FRONT_RIGHT);
            } else {
              OrthoDrawBuffer(G,GL_FRONT); /* draw into the front buffer */
            }
        
        
            glColor3fv(black);
            glBegin(GL_POLYGON);
            glVertex2i(0,I->Height);
            glVertex2i(cBusyWidth,I->Height);
            glVertex2i(cBusyWidth,I->Height-cBusyHeight);
            glVertex2i(0,I->Height-cBusyHeight);
            glVertex2i(0,I->Height); /* needed on old buggy Mesa */
            glEnd();
            
            glColor3fv(white);	 
            
            y=I->Height-cBusyMargin;
            c=I->BusyMessage;
            if(*c) {
              TextSetColor(G,white);
              TextSetPos2i(G,cBusyMargin,y-(cBusySpacing/2));
              TextDrawStr(G,c);
              y-=cBusySpacing;
            }
            
            if(I->BusyStatus[1]) {
              glBegin(GL_LINE_LOOP);
              glVertex2i(cBusyMargin,y);
              glVertex2i(cBusyWidth-cBusyMargin,y);
              glVertex2i(cBusyWidth-cBusyMargin,y-cBusyBar);
              glVertex2i(cBusyMargin,y-cBusyBar);
              glVertex2i(cBusyMargin,y); /* needed on old buggy Mesa */
              glEnd();
              glColor3fv(white);	 
              glBegin(GL_POLYGON);
              glVertex2i(cBusyMargin,y);
              x=(I->BusyStatus[0]*(cBusyWidth-2*cBusyMargin)/I->BusyStatus[1])+cBusyMargin;
              glVertex2i(x,y);
              glVertex2i(x,y-cBusyBar);
              glVertex2i(cBusyMargin,y-cBusyBar);
              glVertex2i(cBusyMargin,y); /* needed on old buggy Mesa */
              glEnd();
              y-=cBusySpacing;
            }
            
            if(I->BusyStatus[3]) {
              glColor3fv(white);	 
              glBegin(GL_LINE_LOOP);
              glVertex2i(cBusyMargin,y);
              glVertex2i(cBusyWidth-cBusyMargin,y);
              glVertex2i(cBusyWidth-cBusyMargin,y-cBusyBar);
              glVertex2i(cBusyMargin,y-cBusyBar);
              glVertex2i(cBusyMargin,y); /* needed on old buggy Mesa */
              glEnd();
              x=(I->BusyStatus[2]*(cBusyWidth-2*cBusyMargin)/I->BusyStatus[3])+cBusyMargin;
              glColor3fv(white);	 
              glBegin(GL_POLYGON);
              glVertex2i(cBusyMargin,y);
              glVertex2i(x,y);
              glVertex2i(x,y-cBusyBar);
              glVertex2i(cBusyMargin,y-cBusyBar);
              glVertex2i(cBusyMargin,y); /* needed on old buggy Mesa */
              glEnd();
              y-=cBusySpacing;
            }
            
            if(!draw_both)
              break;
            if(pass>1)
              break;
            pass++;
          }
        
          glFlush();
          glFinish();
		  
          if(draw_both)
            OrthoDrawBuffer(G,GL_BACK_LEFT);
          else
            OrthoDrawBuffer(G,GL_BACK);      
        }
        OrthoPopMatrix(G);
        OrthoDirty(G);
      }
#endif               

    }
  }

  PRINTFD(G,FB_Ortho)
    " OrthoBusyDraw: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void OrthoRestorePrompt(PyMOLGlobals *G) 
{
  register COrtho *I=G->Ortho;
  int curLine;
  if(!I->InputFlag) 
    {
      if(I->Saved[0]) 
		{
		  if(I->CurChar) {
            OrthoNewLine(G,NULL,true);
		  }
		  curLine = I->CurLine&OrthoSaveLines;
		  strcpy(I->Line[curLine],I->Saved);
		  I->Saved[0]=0;
		  I->CurChar = I->SavedCC;
		  I->PromptChar = I->SavedPC;
		} 
      else 
		{
		  if(I->CurChar) 
            OrthoNewLine(G,I->Prompt,true);
		  else
            {
              curLine = I->CurLine&OrthoSaveLines;
              strcpy(I->Line[curLine],I->Prompt);
              I->CurChar = (I->PromptChar = strlen(I->Prompt));
            }
		}
      I->InputFlag=1;
    }
}
/*========================================================================*/
void OrthoKeyControl(PyMOLGlobals *G,unsigned char k) {
  char buffer[OrthoLineLength];

  /* safer...*/

  sprintf(buffer,"cmd._ctrl(chr(%d))",k);
  /* sprintf(buffer,"_ctrl %c",k); */
  PLog(G,buffer,cPLog_pym);
  PParse(G,buffer);
  PFlush(G);

}
/*========================================================================*/
void OrthoKeyAlt(PyMOLGlobals *G,unsigned char k) {
  char buffer[OrthoLineLength];

  /* safer...*/

  if(k=='@') { 
    /* option G produces '@' on some non-US keyboards, so simply
       ignore the modifier */
    OrthoKey(G,k,0,0,0);
  } else {
    sprintf(buffer,"cmd._alt(chr(%d))",k);
    /* sprintf(buffer,"_alt %c",k); */
    PLog(G,buffer,cPLog_pym);
    PParse(G,buffer);
    PFlush(G);
  }

}
static int add_normal_char(COrtho *I,unsigned char k)
{
  char buffer[OrthoLineLength];
  int curLine=I->CurLine&OrthoSaveLines;
  if(I->CursorChar>=0) {
    strcpy(buffer,I->Line[curLine]+I->CursorChar);
    I->Line[curLine][I->CursorChar]=k;
    I->CursorChar++;
    I->CurChar++;
    strcpy(I->Line[curLine]+I->CursorChar,buffer);
  } else {
    I->Line[curLine][I->CurChar]=k;
    I->CurChar++;
    I->Line[curLine][I->CurChar]=0;
  }
  return curLine;
}
/*========================================================================*/
void OrthoKey(PyMOLGlobals *G,unsigned char k,int x,int y,int mod)
{
  register COrtho *I=G->Ortho;
  char buffer[OrthoLineLength];
  int curLine;

  PRINTFD(G,FB_Ortho)
    " Ortho: %c (%d), x %d y %d, mod %d\n",k,k,x,y,mod
    ENDFD;
    
  if(!I->InputFlag) 
    {
      if(I->Saved[0]) 
		{
		  if(I->CurChar) {
            OrthoNewLine(G,NULL,true);
		  }
		  curLine = I->CurLine&OrthoSaveLines;
		  strcpy(I->Line[curLine],I->Saved);
		  I->Saved[0]=0;
		  I->CurChar = I->SavedCC;
		  I->PromptChar = I->SavedPC;
		} 
      else 
		{
		  if(I->CurChar) 
            OrthoNewLine(G,I->Prompt,true);
		  else
            {
              curLine = I->CurLine&OrthoSaveLines;
              strcpy(I->Line[curLine],I->Prompt);
              I->CurChar = (I->PromptChar = strlen(I->Prompt));
            }
		}
      I->InputFlag=1;
    }
  if(mod==4) { /* alt */
    OrthoKeyAlt(G,k);
  } else if((k>32)&&(k!=127)) {
    curLine=add_normal_char(I,k);
  } else switch(k) {
  case 32: /* spacebar */
    if(SettingGetGlobal_b(G,cSetting_presentation)&&!OrthoTextVisible(G)&&(I->CurChar==I->PromptChar)) {
      PParse(G,"cmd.scene('','next')");
    } else {
      curLine=add_normal_char(I,k);
    }
    break;
  case 127: /* delete */     
#ifndef _PYMOL_OSX
    if((!I->CurChar)||(I->CurChar==I->PromptChar)||!OrthoTextVisible(G)) {
      OrthoKeyControl(G,4+64);
    } else {
      if(I->CursorChar>=0) {
        if(I->CursorChar<I->CurChar)
          I->CursorChar++;
        if(I->CursorChar==I->CurChar)
          I->CursorChar=-1;
      }
      if(I->CurChar>I->PromptChar)
        {
          curLine=I->CurLine&OrthoSaveLines;
          if(I->CursorChar>=0) {
            if(I->CursorChar>I->PromptChar) {
              strcpy(buffer,I->Line[curLine]+I->CursorChar);
              I->CursorChar--;
              I->CurChar--;
              strcpy(I->Line[curLine]+I->CursorChar,buffer);
            }
          } else {
            I->CurChar--;
            I->Line[curLine][I->CurChar]=0;
          }
        }
    } 
    break;
  case 8: /* backspace */
#endif
    if(I->CurChar>I->PromptChar)
      {
        curLine=I->CurLine&OrthoSaveLines;
        if(I->CursorChar>=0) {
          if(I->CursorChar>I->PromptChar) {
            strcpy(buffer,I->Line[curLine]+I->CursorChar);
            I->Line[curLine][I->CursorChar]=k;
            I->CursorChar--;
            I->CurChar--;
            strcpy(I->Line[curLine]+I->CursorChar,buffer);
          }
        } else {
          I->CurChar--;
          I->Line[curLine][I->CurChar]=0;
        }
      }
    break;
  case 5: /* CTRL E -- ending */
    if(OrthoArrowsGrabbed(G)) {
      I->CursorChar=-1;
    } else 
      OrthoKeyControl(G,(unsigned char)(k+64));
    break;
  case 1: /* CTRL A -- beginning */
    if(OrthoArrowsGrabbed(G)) {
      if(I->CurChar)
        I->CursorChar=I->PromptChar;        
    } else 
      OrthoKeyControl(G,(unsigned char)(k+64));
    break;
  case 4: /* CTRL D */
    if((!I->CurChar)||(I->CurChar==I->PromptChar)||!OrthoTextVisible(G)) {
      OrthoKeyControl(G,(unsigned char)(4+64));
    } else if((I->CurChar>I->PromptChar)&&
              (I->CursorChar>=0)&&
              (I->CursorChar<I->CurChar)) { /* deleting */
      curLine=I->CurLine&OrthoSaveLines;
      strcpy(buffer,I->Line[curLine]+I->CursorChar+1);
      I->CurChar--;
      strcpy(I->Line[curLine]+I->CursorChar,buffer);
    } else { /* filename completion query */
      curLine=I->CurLine&OrthoSaveLines;
      if(I->PromptChar) {
        strcpy(buffer,I->Line[curLine]);
        if(PComplete(G,buffer+I->PromptChar,
                     sizeof(OrthoLineType)-I->PromptChar)); /* just print, don't complete */
      }
    }
    break;
  case 9: /* CTRL I -- tab */
    if(mod&cOrthoCTRL) {
      OrthoKeyControl(G,(unsigned char)(k+64)); 
    } else {
      curLine=I->CurLine&OrthoSaveLines;
      if(I->PromptChar) {
        strcpy(buffer,I->Line[curLine]);
          
        if(PComplete(G,buffer+I->PromptChar,
                     sizeof(OrthoLineType)-I->PromptChar))
          {
            OrthoRestorePrompt(G);
            curLine=I->CurLine&OrthoSaveLines;
            strcpy(I->Line[curLine],buffer);
            I->CurChar = strlen(I->Line[curLine]);
          }
      }
    }
    break;
  case 27: /* ESCAPE */
    if(SettingGetGlobal_b(G,cSetting_presentation)&&!(mod&(cOrthoCTRL||cOrthoSHIFT))) {
      PParse(G,"_quit");
    } else {
      if(I->SplashFlag) {
        OrthoRemoveSplash(G);
      } else {
        if(mod&cOrthoSHIFT) 
          SettingSet(G,cSetting_overlay,(float)(!((int)SettingGet(G,cSetting_overlay))));
        else
          SettingSet(G,cSetting_text,(float)(!((int)SettingGet(G,cSetting_text))));
      }
    }
    break;
  case 13: /* CTRL M -- carriage return */
    OrthoParseCurrentLine(G);
    break;
  case 11: /* CTRL K -- truncate */
    if(OrthoArrowsGrabbed(G)) {
      if(I->CursorChar>=0) { 
        I->Line[I->CurLine&OrthoSaveLines][I->CursorChar]=0;
        I->CurChar=I->CursorChar;
        I->CursorChar=-1;
      }
    } else {
      if(mod&cOrthoCTRL) {
        OrthoKeyControl(G,(unsigned char)(k+64));
      }
    }
    break;
  case 22: /* CTRL V -- paste */
#ifndef _PYMOL_NOPY
    PBlockAndUnlockAPI(G);
    PRunStringInstance(G,"cmd.paste()");
    PLockAPIAndUnblock(G);
#endif
    break;
  default:
    OrthoKeyControl(G,(unsigned char)(k+64));
    break;
  }
  PyMOL_NeedRedisplay(G->PyMOL);
}
/*========================================================================*/
void OrthoParseCurrentLine(PyMOLGlobals *G) 
{
  register COrtho *I=G->Ortho;
  char buffer[OrthoLineLength];
  int curLine;
  curLine=I->CurLine&OrthoSaveLines;
  I->Line[curLine][I->CurChar]=0;
  strcpy(buffer,I->Line[curLine]+I->PromptChar);
#ifndef _PYMOL_NOPY
  if(buffer[0]) {
    strcpy(I->History[I->HistoryLine],buffer);
    I->HistoryLine = (I->HistoryLine+1)&OrthoHistoryLines;
    I->History[I->HistoryLine][0]=0;
    I->HistoryView=I->HistoryLine;
    OrthoNewLine(G,NULL,true);
    if(WordMatch(G,buffer,"quit",true)==0) /* don't log quit */
      PLog(G,buffer,cPLog_pml);
    OrthoDirty(G); /* this will force a redraw, if necessary */
    PParse(G,buffer);
    OrthoRestorePrompt(G);
  }
#endif
  I->CursorChar=-1;
}
/*========================================================================*/
void OrthoAddOutput(PyMOLGlobals *G,char *str)
{
  register COrtho *I=G->Ortho;
  int curLine;
  char *p,*q;
  int cc;
  int wrap;
  curLine = I->CurLine&OrthoSaveLines;
  if(I->InputFlag)
	 {
		strcpy(I->Saved,I->Line[curLine]);
		I->SavedPC=I->PromptChar;
		I->SavedCC=I->CurChar;
		I->PromptChar=0;
		I->CurChar=0;
		I->Line[curLine][0]=0;
		I->InputFlag=0;
	 }
  curLine = I->CurLine&OrthoSaveLines;
  p=str;
  q=I->Line[curLine]+I->CurChar;
  cc=I->CurChar;
  while(*p)
	 {
		if(*p>=32)
		  {
			 cc++;
          wrap = (int)SettingGet(G,cSetting_wrap_output);

          if(wrap>0) {
            if(cc>wrap)
              {
                *q=0;
                I->CurChar = cc;
                OrthoNewLine(G,NULL,true);
                cc=0;
                q=I->Line[I->CurLine&OrthoSaveLines];
                curLine = I->CurLine&OrthoSaveLines;
              }
          } 
          if(cc>=OrthoLineLength-6) { /* fail safe */
            *q=0;
            I->CurChar = cc;
            OrthoNewLine(G,NULL,false);
            cc=0;
            q=I->Line[I->CurLine&OrthoSaveLines];
            curLine = I->CurLine&OrthoSaveLines;
          }
			 *q++=*p++;
		  }
		else if((*p==13)||(*p==10))
		  {
			 *q=0;
			 I->CurChar = cc;
			 OrthoNewLine(G,NULL,true);
			 q=I->Line[I->CurLine&OrthoSaveLines];
			 curLine = I->CurLine&OrthoSaveLines;
			 p++;
			 cc=0;
		  }
		else
		  p++;
	 }
  *q=0;
  I->CurChar = strlen(I->Line[curLine]);
  if((SettingGet(G,cSetting_internal_feedback)>1)||SettingGet(G,cSetting_overlay))
    OrthoDirty(G);
}
/*========================================================================*/
void OrthoNewLine(PyMOLGlobals *G,char *prompt,int crlf)
{
  int curLine;
  register COrtho *I=G->Ortho;

  /*  printf("orthoNewLine: CC: %d CL:%d PC: %d IF:L %d\n",I->CurChar,I->CurLine,
		I->PromptChar,I->InputFlag);*/
  /*  if(I->CurChar)
	 {*/
  if(I->CurChar)
    OrthoFeedbackIn(G,I->Line[I->CurLine&OrthoSaveLines]);
  else
    OrthoFeedbackIn(G," ");
  if(Feedback(G,FB_Python,FB_Output)) {
    if(crlf) {
      printf("%s\n",I->Line[I->CurLine&OrthoSaveLines]);
    } else {
      printf("%s",I->Line[I->CurLine&OrthoSaveLines]);
    }
    fflush(stdout);
  }
      /*	 }*/

      /*  if(I->Line[I->CurLine&OrthoSaveLines][0])*/
  I->CurLine++;
  curLine = I->CurLine&OrthoSaveLines;

  if(prompt)
	 {
		strcpy(I->Line[curLine],prompt);
		I->CurChar = (I->PromptChar = strlen(prompt));
		I->InputFlag=1;
	 }
  else
	 {
		I->CurChar = 0;
		I->Line[curLine][0] = 0;
		I->PromptChar = 0;
		I->InputFlag = 0;
	 }
  /*printf("orthoNewLine: CC: %d CL:%d PC: %d IF:L %d\n",I->CurChar,I->CurLine,
	 I->PromptChar,I->InputFlag);*/

}
/*========================================================================*/
void OrthoGrab(PyMOLGlobals *G,Block *block)
{
  register COrtho *I=G->Ortho;
  I->GrabbedBy = block;
}
/*========================================================================*/
void OrthoUngrab(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  I->GrabbedBy = NULL;
}

/*========================================================================*/
Block *OrthoNewBlock(PyMOLGlobals *G,Block *block)
{
  if(!block)
	 ListElemAlloc(G,block,Block);
  UtilZeroMem(block,sizeof(Block));
  BlockInit(G,block);
  return(block);
}
/*========================================================================*/
void OrthoFreeBlock(PyMOLGlobals *G,Block *block)
{
  if(block) 
    ListElemFree(block);
}
/*========================================================================*/
void OrthoAttach(PyMOLGlobals *G,Block *block,int type)
{
  register COrtho *I=G->Ortho;
  ListInsert(I->Blocks,block,NULL,next,Block);
}
/*========================================================================*/
void OrthoDetach(PyMOLGlobals *G,Block *block)
{
  register COrtho *I=G->Ortho;
  if(I->GrabbedBy == block)
    I->GrabbedBy = NULL;
  ListDetach(I->Blocks,block,next,Block);
}
/*========================================================================*/
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef PYMOL_EVAL
#include "OrthoEvalMessage.h"
#endif
#ifdef PYMOL_BETA
#include "OrthoBetaMessage.h"
#endif
#ifdef JYMOL_EVAL
#include "OrthoJyMolEvalMessage.h"
#endif
#ifdef PYMOL_EDU
#include "OrthoEduMessage.h"
#endif
#ifdef PYMOL_COLL
#include "OrthoCollMessage.h"
#endif
#ifdef _PYMOL_ACTIVEX
#include "OrthoAxMessage.h"
#endif
/* END PROPRIETARY CODE SEGMENT */

void OrthoDoDraw(PyMOLGlobals *G,int render_mode)
{
  register COrtho *I=G->Ortho;

  int x,y;
  int l,lcount;
  char *str;
  float *v;
  int showLines;
  int height;
  int overlay,text;
  int rightSceneMargin;
  int internal_feedback;
  int times=1;
  int double_pump=false;
  float *bg_color;
  int skip_prompt = 0;
  int render = false;

  I->RenderMode = render_mode;
  if(SettingGetGlobal_b(G,cSetting_seq_view)) {
    SeqUpdate(G); 
    I->HaveSeqViewer = true;
  } else if(I->HaveSeqViewer) {
    SeqUpdate(G);
    I->HaveSeqViewer = false;
  }

  if(SettingGet_i(G,NULL,NULL,cSetting_internal_prompt))
    skip_prompt = 0;
  else
    skip_prompt = 1;
     
  double_pump=SettingGet_i(G,NULL,NULL,cSetting_stereo_double_pump_mono);
  bg_color=SettingGet_3fv(G,NULL,NULL,cSetting_bg_rgb);

  I->OverlayColor[0]=1.0F-bg_color[0];
  I->OverlayColor[1]=1.0F-bg_color[1];
  I->OverlayColor[2]=1.0F-bg_color[2];
  if(diff3f(I->OverlayColor,bg_color)<0.25)
    zero3f(I->OverlayColor);

  PRINTFD(G,FB_Ortho)
    " OrthoDoDraw: entered.\n"
    ENDFD;
  if(G->HaveGUI && G->ValidContext) {

    if(Feedback(G,FB_OpenGL,FB_Debugging))
      PyMOLCheckOpenGLErr("OrthoDoDraw checkpoint 0");
    
    if(SettingGetGlobal_b(G,cSetting_internal_gui)) {
      switch(SettingGetGlobal_i(G,cSetting_internal_gui_mode)) {
      case 0:
        rightSceneMargin=(int)SettingGet(G,cSetting_internal_gui_width);
        break;
      default:
        rightSceneMargin = 0;
        break;
      }
    } else {
      rightSceneMargin=0;
    }

    internal_feedback=(int)SettingGet(G,cSetting_internal_feedback);

    v=SettingGetfv(G,cSetting_bg_rgb);
    overlay = (int)SettingGet(G,cSetting_overlay);
    if(overlay==1) {
      overlay = (int)SettingGet(G,cSetting_overlay_lines);
    }
    text = (int)SettingGet(G,cSetting_text);

    if(text) overlay=0;
    
    {
      float alpha = (SettingGetGlobal_b(G,cSetting_opaque_background) ? 1.0F : 0.0F);
      glClearColor(v[0],v[1],v[2],alpha);
    }

    if(overlay||(!text)) 
      if(!SceneRenderCached(G))
        render=true;
    
    if(render_mode<2) {
      if(SceneMustDrawBoth(G)) {
        OrthoDrawBuffer(G,GL_BACK_LEFT);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        OrthoDrawBuffer(G,GL_BACK_RIGHT);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        times = 2;
        double_pump = true;
      } else {
        OrthoDrawBuffer(G,GL_BACK);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        times = 1;
        double_pump=false;
      }    
    } else {
      times = 1;
      double_pump=false;
    }

    I->DrawTime = -I->LastDraw;
    I->LastDraw = UtilGetSeconds(G);
    I->DrawTime += I->LastDraw;
    ButModeSetRate(G,(float)I->DrawTime); 
    
    if(render&&(render_mode<2))
      SceneRender(G,NULL,0,0,NULL,0,0,0,SettingGetGlobal_b(G,cSetting_image_copy_always));

    glClearColor(0.0,0.0,0.0,1.0);
    
    while(times--) {

      switch(times) {
      case 1:
        OrthoDrawBuffer(G,GL_BACK_LEFT);

        break;
      case 0:
        if(double_pump) {
          OrthoDrawBuffer(G,GL_BACK_RIGHT);
        } else
          OrthoDrawBuffer(G,GL_BACK);
        break;
      }

      OrthoPushMatrix(G);
      
      x = I->X;
      y = I->Y;
      
      
      if(I->DrawText&&internal_feedback) { /* moved to avoid conflict with menus */
        glColor3f(0.0,0.0,0.0);
        glBegin(GL_POLYGON);
        height=(internal_feedback-1)*cOrthoLineHeight+cOrthoBottomSceneMargin;
        glVertex2i(I->Width-rightSceneMargin,height-1);
        glVertex2i(I->Width-rightSceneMargin,0);
        glVertex2i(0,0);
        glVertex2i(0,cOrthoBottomSceneMargin-1);
        glEnd();
      }

      
      PRINTFD(G,FB_Ortho)
        " OrthoDoDraw: drawing blocks...\n"
        ENDFD;

      if((int)SettingGet(G,cSetting_text)||I->SplashFlag) {
        Block *block;
        int active_tmp;
        block=SeqGetBlock(G);
        active_tmp = block->active;
        block->active = false; 
        BlockRecursiveDraw(I->Blocks);
        block->active = active_tmp;
      } else {
        BlockRecursiveDraw(I->Blocks);
      }
      
      PRINTFD(G,FB_Ortho)
        " OrthoDoDraw: blocks drawn.\n"
        ENDFD;
      
      if(I->LoopFlag) {
        glColor3f(1.0,1.0,1.0);
        glBegin(GL_LINE_LOOP);
        glVertex2i(I->LoopRect.left,I->LoopRect.top);
        glVertex2i(I->LoopRect.right,I->LoopRect.top);
        glVertex2i(I->LoopRect.right,I->LoopRect.bottom);
        glVertex2i(I->LoopRect.left,I->LoopRect.bottom);
        glVertex2i(I->LoopRect.left,I->LoopRect.top);
        glEnd();
      }
      
      OrthoRestorePrompt(G);
      
      if(I->DrawText) {	 
        /* now print the text */
        
        lcount = 0;
        x = cOrthoLeftMargin;
        y = cOrthoBottomMargin;

#ifdef _PYMOL_SHARP3D
        if(SceneGetStereo(G)&&SettingGetGlobal_b(G,cSetting_overlay)) {
          y+=(7*cOrthoLineHeight)/10;
        }
#endif
        if((int)SettingGet(G,cSetting_text)||I->SplashFlag)
          showLines=I->ShowLines;
        else {
          int overlay2;
          overlay2 = (int)SettingGet(G,cSetting_overlay);
          if(overlay2==1) {
            overlay2 = (int)SettingGet(G,cSetting_overlay_lines);
          }
          showLines=internal_feedback+overlay2;
        }

        l=(I->CurLine-(lcount+skip_prompt))&OrthoSaveLines;

        glColor3fv(I->TextColor);
        while(l>=0)
          {
            lcount++;
            if(lcount>showLines)
              break;
            str = I->Line[l&OrthoSaveLines];
            if(strncmp(str,I->Prompt,6)==0)
              TextSetColor(G,I->TextColor);            
            else
              TextSetColor(G,I->OverlayColor);
            TextSetPos2i(G,x,y);
            if(str)
              {
                TextDrawStr(G,str);
                if((lcount==1)&&(I->InputFlag)) 
                  {
                    if(!skip_prompt) {
                      if(I->CursorChar>=0) {
                        TextSetPos2i(G,x+8*I->CursorChar,y);
                      }
                      TextDrawChar(G,'_');
                    }
                  }
              }
            l=(I->CurLine-(lcount+skip_prompt))&OrthoSaveLines;
            y=y+cOrthoLineHeight;
          }
      }
      
      OrthoDrawWizardPrompt(G);
 
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef PYMOL_EVAL
      OrthoDrawEvalMessage(G);
#endif
#ifdef PYMOL_BETA
      OrthoDrawBetaMessage(G);
#endif
#ifdef JYMOL_EVAL
      OrthoDrawEvalMessage(G);
#endif
#ifdef PYMOL_EDU
      OrthoDrawEduMessage(G);
#endif
#ifdef PYMOL_COLL
      OrthoDrawCollMessage(G);
#endif
#ifdef _PYMOL_ACTIVEX
      OrthoDrawAxMessage(G);
#endif

/* END PROPRIETARY CODE SEGMENT */


      OrthoPopMatrix(G);

      if(Feedback(G,FB_OpenGL,FB_Debugging))
        PyMOLCheckOpenGLErr("OrthoDoDraw final checkpoint");

    } /* while */

  }

  I->DirtyFlag =false;
  PRINTFD(G,FB_Ortho)
    " OrthoDoDraw: leaving...\n"
    ENDFD;

}
/*========================================================================*/

static void OrthoDrawWizardPrompt(PyMOLGlobals *G)
{
  /* assumes PMGUI */

  register COrtho *I=G->Ortho;
  
  char *vla,*p;
  int nLine;
  int x,y,xx;
  int nChar,c,ll;
  int maxLen;
  BlockRect rect;
  int prompt_mode = SettingGetGlobal_i(G,cSetting_wizard_prompt_mode);

  if(I->WizardPromptVLA && prompt_mode) {
    vla = I->WizardPromptVLA;
    
    nLine = UtilCountStringVLA(vla);
    if(nLine) {
      nChar = VLAGetSize(I->WizardPromptVLA);
      
      /* count max line length */
      
      maxLen = 0;
      p = vla;
      ll = 0;
      c=nChar;
      while(c>0) {
        if(!*p) {
          if(maxLen<ll)
            maxLen = ll;
          ll=0;
          p++;
          c--;
        } else if(((*p)=='\\')&& /* color encoded */
                  (p[1]>='0')&&(p[1]<='9')&&
                  (p[2]>='0')&&(p[2]<='9')&&
                  (p[3]>='0')&&(p[3]<='9')) /* relying upon short-circuit logic to avoid overrun */ {
          p+=4;
          c-=4;
        } else {
          ll++;
          p++;
          c--;
        }
      }
      
      rect.top = I->Height;
      if(I->HaveSeqViewer)
        if(!SettingGetGlobal_b(G,cSetting_seq_view_location)) {
          rect.top -= SeqGetHeight(G);
        }

      if(prompt_mode!=3) {
        rect.top -= cWizardTopMargin;
        rect.left = cWizardLeftMargin;
      } else {
        rect.top -= 1;
        rect.left = 1;
      }

      rect.bottom = rect.top-(nLine*cOrthoLineHeight+2*cWizardBorder)-2;
      rect.right = rect.left + cOrthoCharWidth*maxLen + 2*cWizardBorder+1;
      
      if(prompt_mode==1) {
        glColor3fv(I->WizardBackColor);
        
        glBegin(GL_POLYGON);
        glVertex2i(rect.right,rect.top);
        glVertex2i(rect.right,rect.bottom);
        glVertex2i(rect.left,rect.bottom);
        glVertex2i(rect.left,rect.top);
        glEnd();
      }

      glColor3fv(I->WizardTextColor);
      
      x = rect.left+cWizardBorder;
      y = rect.top-(cWizardBorder+cOrthoLineHeight);

      vla = I->WizardPromptVLA;
      
      /* count max line length */
      
      TextSetColor(G,I->WizardTextColor);
      TextSetPos2i(G,x,y);
      xx = x;
      p = vla;
      ll = 0;
      c=nChar;
      while(c>0) {
        if(*p) {
          if((*p=='\\')&&(*(p+1))&&(*(p+2))&&(*(p+3))) {
            if(*(p+1)=='-') {
              TextSetColor(G,I->WizardTextColor);
              p+=4;
              c-=4;
            } else {
              TextSetColor3f(G,(*(p+1)-'0')/9.0F,(*(p+2)-'0')/9.0F,(*(p+3)-'0')/9.0F);
              p+=4;
              c-=4;
            }
            TextSetPos2i(G,xx,y);
          }
        }
        if(c--) {
          if(*p) {
            TextDrawChar(G,*p);
            xx = xx + 8;
          }
          if(!*(p++)) {
            y=y-cOrthoLineHeight;
            xx = x;
            TextSetPos2i(G,x,y);
          }
        }
      }
    }
  }
}

/*========================================================================*/
void OrthoReshape(PyMOLGlobals *G,int width, int height,int force)
{
  register COrtho *I=G->Ortho;

  Block *block = NULL;
  int sceneBottom,sceneRight = 0;
  int internal_gui_width;
  int internal_feedback;
  int sceneTop = 0;

  PRINTFD(G,FB_Ortho)
    " OrthoReshape-Debug: %d %d\n",width,height
    ENDFD;

  if((width>0)&&(SettingGetGlobal_i(G,cSetting_stereo_mode)==4)) {
    width = width / 2;
    I->WrapXFlag = true;
  } else {
    I->WrapXFlag = false;
  }

  if((width!=I->Width)||(height!=I->Height)||force) {
  if(width<0) width=I->Width;
  if(height<0) height=I->Height;

  I->Height=height;
  I->Width=width;
  I->ShowLines = height/cOrthoLineHeight;
  
  internal_feedback = (int)SettingGet(G,cSetting_internal_feedback);
  if(internal_feedback)
    sceneBottom = (internal_feedback-1)*cOrthoLineHeight + cOrthoBottomSceneMargin;
  else
    sceneBottom = 0;
    
  internal_gui_width = (int)SettingGet(G,cSetting_internal_gui_width);
  if(!SettingGetGlobal_b(G,cSetting_internal_gui)) {
    internal_gui_width = 0;
    sceneRight = 0;
  } else {
    switch(SettingGetGlobal_i(G,cSetting_internal_gui_mode)) {
    case 1:
      sceneRight = 0;
      break;
    default:
      sceneRight = internal_gui_width;
      break;
    }
  }


  {
    int seqHeight;
    block=SeqGetBlock(G);
    block->active=true;
    
    if(SettingGetGlobal_b(G,cSetting_seq_view_location)) {
      
      BlockSetMargin(block,height-sceneBottom-10,0,sceneBottom,sceneRight);
      if(block->fReshape)
        block->fReshape(block,width,height);			
      seqHeight = SeqGetHeight(G);
      BlockSetMargin(block,height-sceneBottom-seqHeight,0,sceneBottom,sceneRight);
      if(!SettingGetGlobal_b(G,cSetting_seq_view_overlay)) {
        sceneBottom +=seqHeight;
      }
      
    } else {
      
      BlockSetMargin(block,0,0,height-10,sceneRight);
      if(block->fReshape)
        block->fReshape(block,width,height);			
      seqHeight = SeqGetHeight(G);
      BlockSetMargin(block,0,0,height-seqHeight,sceneRight);
      if(!SettingGetGlobal_b(G,cSetting_seq_view_overlay)) {
        sceneTop = seqHeight;
      }
    }
  }

  {
    int WizardMargin = WizardMargin1;
    
    if(!SettingGet(G,cSetting_mouse_grid)) {
      WizardMargin = WizardMargin2;
    }
    if(SettingGet(G,cSetting_internal_gui)) {

#ifndef _PYMOL_NOPY
      block=ExecutiveGetBlock(G);
      block->active=true;
      BlockSetMargin(block,0,width-internal_gui_width,WizardMargin,0);
      block=WizardGetBlock(G);
      BlockSetMargin(block,height-WizardMargin+1,width-internal_gui_width,WizardMargin,0);
      block->active=false;
      block=ButModeGetBlock(G);
      BlockSetMargin(block,height-WizardMargin+1,width-internal_gui_width,ButModeMargin,0);
      block->active=true;
#else
      block=ExecutiveGetBlock(G);
      block->active=true;
      BlockSetMargin(block,0,width-internal_gui_width,ButModeMargin,0);
      block=WizardGetBlock(G);
      BlockSetMargin(block,height-WizardMargin+1,width-internal_gui_width,ButModeMargin,0);
      block->active=false;
      block=ButModeGetBlock(G);
      BlockSetMargin(block,height-WizardMargin+1,width-internal_gui_width,ButModeMargin,0);
      block->active=false;
#endif

      block=ControlGetBlock(G);
      BlockSetMargin(block,height-ButModeMargin+1,width-internal_gui_width,ControlMargin,0);
      block->active=true;
    } else {
      block=ExecutiveGetBlock(G);
      block->active=false;
      BlockSetMargin(block,0,width-internal_gui_width,WizardMargin,0);
      block=WizardGetBlock(G);
      BlockSetMargin(block,height-WizardMargin+1,width-internal_gui_width,WizardMargin,0);
      block->active=false;
      block=ButModeGetBlock(G);
      BlockSetMargin(block,height-WizardMargin+1,width-internal_gui_width,ButModeMargin,0);
      block->active=false;
      block=ControlGetBlock(G);
      BlockSetMargin(block,height-ButModeMargin+1,width-internal_gui_width,ControlMargin,0);
      block->active=false;
    }

  }
  block=SceneGetBlock(G);
  BlockSetMargin(block,sceneTop,0,sceneBottom,sceneRight);

  block=NULL;
  while(ListIterate(I->Blocks,block,next))
	 if(block->fReshape) {
		block->fReshape(block,width,height);			
    }

  WizardRefresh(G); /* safe to call even if no wizard exists */
  }
}

/*========================================================================*/
void OrthoReshapeWizard(PyMOLGlobals *G,ov_size wizHeight)
{
  Block *block;
  register COrtho *I=G->Ortho;
  int height,width;
  int internal_gui_width;

  height=I->Height;
  width=I->Width;

  if(SettingGet(G,cSetting_internal_gui)>0.0) {
    int WizardMargin = WizardMargin1;
    internal_gui_width = (int)SettingGet(G,cSetting_internal_gui_width);
    block=ExecutiveGetBlock(G);

    if(!SettingGet(G,cSetting_mouse_grid)) {
      WizardMargin = WizardMargin2;
    }

    if(height) {
      int wh=wizHeight;
      if(wh) wh++;
      BlockSetMargin(block,0,width-internal_gui_width,WizardMargin+wh,0);
    } else {
      BlockSetMargin(block,0,width-internal_gui_width,WizardMargin,0);
    }
    block->fReshape(block,width,height);

    block=WizardGetBlock(G);

    if(wizHeight) {
      BlockSetMargin(block,height-(WizardMargin+wizHeight),width-internal_gui_width,WizardMargin,0);
      block->active=true;
    } else {
      BlockSetMargin(block,height-WizardMargin,width-internal_gui_width,WizardMargin,0);
      block->active=false;
    }
    block->fReshape(block,width,height);
  }
}

/*========================================================================*/
Block *OrthoFindBlock(PyMOLGlobals *G,int x,int y)
{
  register COrtho *I=G->Ortho;

  return(BlockRecursiveFind(I->Blocks,x,y));
}
/*========================================================================*/
int OrthoGetWrapClickSide(PyMOLGlobals *G)
{
  return G->Ortho->WrapClickSide;
}
/*========================================================================*/
int OrthoButton(PyMOLGlobals *G,int button,int state,int x,int y,int mod)
{
  register COrtho *I=G->Ortho;

  Block *block=NULL;
  int handled = 0; 


  switch(button) {
  case 3:
  case 4:
    block = SceneGetBlock(G);
    break;
  }

  if(I->WrapXFlag) {
    if(state==P_GLUT_DOWN) {
      x = get_wrap_x(x,NULL,G->Option->winX,&I->WrapClickSide);
    } else {
      x = get_wrap_x(x,&I->LastX,G->Option->winX,&I->WrapClickSide);
    }
  } else {
    I->WrapClickSide = 0;
  }

  OrthoRemoveSplash(G);
  I->X=x;
  I->Y=y;
  I->LastX = x;
  I->LastY = y;
  I->LastModifiers = mod;

  if(state==P_GLUT_DOWN) {
    I->ActiveButton = button;
    if(I->GrabbedBy) {
      if(I->GrabbedBy->inside)
        block = BlockRecursiveFind(I->GrabbedBy->inside,x,y);
      else
        block = I->GrabbedBy;
    } else if(!block)
      block = OrthoFindBlock(G,x,y);
    if(block) {
      I->ClickedIn = block;
      if(block->fClick) {
        handled = block->fClick(block,button,x,y,mod);
      }
    }
  } else if(state==P_GLUT_UP) {
    if(I->GrabbedBy) {
      block=I->GrabbedBy;
      if(block->fRelease)
        handled = block->fRelease(block,button,x,y,mod);
      I->ClickedIn = NULL;
    }
    if(I->ClickedIn) {
      block=I->ClickedIn;
      if(block->fRelease)
        handled = block->fRelease(block,button,x,y,mod);
      I->ClickedIn = NULL;
    }
  }
#if 0
  if(block&&!handled) {
    if(SceneGetBlock(G)==block) {
      if(state==P_GLUT_DOWN) {
        I->LoopRect.left=x;
        I->LoopRect.top=y;
        I->LoopRect.right=x;
        I->LoopRect.bottom=y;
        I->LoopFlag=true;
        I->LoopMod = mod;
        I->GrabbedBy=&I->LoopBlock;
        OrthoDirty(G);
      } 
    }
  }
#endif
  return(handled);
}
/*========================================================================*/ 
int OrthoDrag(PyMOLGlobals *G,int x, int y,int mod)
{
  register COrtho *I=G->Ortho;

  Block *block=NULL;
  int handled = 0;

 if(I->WrapXFlag) {
   x = get_wrap_x(x,&I->LastX,G->Option->winX, NULL);
 }

  I->LastX = x;
  I->LastY = y;
  I->LastModifiers = mod;

  I->X=x;
  I->Y=y;
  if(I->GrabbedBy) 
    {
		block = I->GrabbedBy;
		if(block->fDrag)
        handled = block->fDrag(block,x,y,mod);
    }
  else if(I->ClickedIn)
	 {
		block = I->ClickedIn;
		if(block->fDrag)
        handled = block->fDrag(block,x,y,mod);
	 }
  return(handled);
}

/*========================================================================*/
void OrthoSplash(PyMOLGlobals *G) 
{
/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */ 
#ifdef _PYMOL_IP_SPLASH
#include"OrthoIPSplash.h"
#else
  if(G->Option->incentive_product) {
    PRINTF " PyMOL(TM) Incentive Product - Copyright (C) 2008 DeLano Scientific LLC.\n \n" ENDF(G);
    PRINTF " This Executable Build integrates and extends Open-Source PyMOL " ENDF(G);
    PRINTF _PyMOL_VERSION ENDF(G);
    PRINTF ".\n" ENDF(G);
  } else 
/* END PROPRIETARY CODE SEGMENT */
    {
      /* Splash message for unrestricted access open-source versions... */
      PRINTF " PyMOL(TM) Molecular Graphics System, Version " ENDF(G);
      PRINTF _PyMOL_VERSION ENDF(G);
      PRINTF ".\n" ENDF(G);
      PRINTF " Copyright (C) 2008 by DeLano Scientific LLC.\n All Rights Reserved.\n \n" ENDF(G);
      
      PRINTF "    Created by Warren L. DeLano, Ph.D. \n \n" ENDF(G);
      
      PRINTF "    Other Major Authors and Contributors:\n\n" ENDF(G);
      PRINTF "       Ralf W. Grosse-Kunstleve, Ph.D.\n \n" ENDF(G); /* note: sglite no longer used */
      
      PRINTF "    PyMOL is user-supported open-source software.  Although some versions\n" ENDF(G);
      PRINTF "    are freely available, PyMOL is not in the public domain.\n \n" ENDF(G);
      
      PRINTF "    If PyMOL is helpful in your work or study, then please volunteer \n" ENDF(G);
      PRINTF "    support for our ongoing efforts to create open and affordable scientific\n" ENDF(G);
      PRINTF "    software by purchasing a PyMOL Maintenance and/or Support subscription.\n\n" ENDF(G);
      
      PRINTF "    More information can be found at \"http://www.pymol.org\".\n \n" ENDF(G);
      
      PRINTF "    Enter \"help\" for a list of commands.\n" ENDF(G);
      PRINTF 
        "    Enter \"help <command-name>\" for information on a specific command.\n\n"
        ENDF(G);
      
      PRINTF " Hit ESC anytime to toggle between text and graphics.\n\n" ENDF(G);
    }
#endif
}
/*========================================================================*/
int OrthoInit(PyMOLGlobals *G,int showSplash)
{
  register COrtho *I=NULL;

  if( (I=(G->Ortho=Calloc(COrtho,1)))) {


  ListInit(I->Blocks);

  I->Pushed = 0;
  {
    int a;
    for(a=0;a<=CMD_QUEUE_MASK;a++) 
      I->cmdQueue[a] = QueueNew(G,0x7FFF); /* 32K ea. level for commands */
    I->cmdActiveQueue = I->cmdQueue[0];
    I->cmdNestLevel = 0;
  }
  I->feedback = QueueNew(G,0x3FFFF); /* ~256K for output */
  I->deferred = NULL;
  I->RenderMode = 0;
  I->WrapXFlag = false;

  I->WizardBackColor[0]=0.2F;
  I->WizardBackColor[1]=0.2F;
  I->WizardBackColor[2]=0.2F;
  I->WizardTextColor[0]=0.2F;
  I->WizardTextColor[1]=1.0F;
  I->WizardTextColor[2]=0.2F;

  I->GrabbedBy = NULL;
  I->ClickedIn = NULL;
  I->DrawText=1;
  I->HaveSeqViewer = false;
  I->TextColor[0]=0.82F;
  I->TextColor[1]=0.82F;
  I->TextColor[2]=1.0;
  I->OverlayColor[0]=1.0;
  I->OverlayColor[1]=1.0;
  I->OverlayColor[2]=1.0;
  I->CurLine=1000;
  I->PromptChar=0;
  I->CurChar=0;
  I->CurLine=0;
  I->CursorChar=-1;
  I->HistoryLine=0;
  I->HistoryView=0;
  I->Line[I->CurLine&OrthoSaveLines][I->CurChar]=0;
  I->WizardPromptVLA=NULL;
  I->SplashFlag = false;
  I->ShowLines = 1;
  I->Saved[0]=0;
  I->DirtyFlag = true;
  I->ActiveGLBuffer = GL_NONE;
  I->LastDraw = UtilGetSeconds(G);
  I->DrawTime = 0.0;
  if(showSplash) {
	 OrthoSplash(G);
    I->SplashFlag=true;
  }
  /*  OrthoFeedbackIn(G," ");*/
  I->CurLine++;
  strcpy(I->Prompt,"PyMOL>");
  strcpy(I->Line[I->CurLine],I->Prompt);
  I->CurChar = (I->PromptChar = strlen(I->Prompt));
  I->InputFlag=1;

  /*printf("orthoNewLine: CC: %d CL:%d PC: %d IF:L %d\n",I->CurChar,I->CurLine,
	 I->PromptChar,I->InputFlag);*/
  
  PopInit(G);
  {
    int a;
    for(a=0;a<=OrthoHistoryLines;a++)
      I->History[a][0]=0;
  }

  return 1;
  } else {
    return 0;
  }
}

/*========================================================================*/
void OrthoFree(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;

  VLAFreeP(I->WizardPromptVLA);
  PopFree(G);
  {
    int a;
    I->cmdActiveQueue = NULL;
    for(a=0;a<=CMD_QUEUE_MASK;a++) 
      QueueFree(I->cmdQueue[a]);
    I->cmdQueue[a] = NULL;
  }
  QueueFree(I->feedback);
  I->feedback=NULL;
  if(I->deferred) {
    DeferredFree(I->deferred);
    I->deferred = NULL;
  }
  FreeP(G->Ortho);
}
/*========================================================================*/
void OrthoPushMatrix(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;

  if(G->HaveGUI && G->ValidContext) {

    if(!I->Pushed) {
      glGetIntegerv(GL_VIEWPORT,I->ViewPort);
    }
    switch(I->RenderMode) {
    case 1:
      glViewport(I->ViewPort[0],I->ViewPort[1],I->ViewPort[2],I->ViewPort[3]);
      break;
    case 2:
      glViewport(I->ViewPort[0]+I->ViewPort[2],I->ViewPort[1],
                 I->ViewPort[2],I->ViewPort[3]);
      break;
    default:
      glViewport(I->ViewPort[0],I->ViewPort[1],I->ViewPort[2],I->ViewPort[3]);
    }

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0,I->ViewPort[2],0,I->ViewPort[3],-100,100);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0.33F,0.33F,0.0F); /* this generates better 
                                       rasterization on macs */

    if(!SettingGetGlobal_b(G,cSetting_texture_fonts)){
      glDisable(GL_ALPHA_TEST);
    } else {
      glEnable(GL_ALPHA_TEST);
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
    if(G->Option->multisample)    
      glDisable(0x809D); /* GL_MULTISAMPLE_ARB */
    
    I->Pushed++;
  }
  /*  glDisable(GL_ALPHA_TEST);
  glDisable(GL_CULL_FACE);
  glDisable(GL_POINT_SMOOTH);*/
  
}
/*========================================================================*/
void OrthoPopMatrix(PyMOLGlobals *G)
{
  register COrtho *I=G->Ortho;
  if(G->HaveGUI && G->ValidContext) {

    if(I->Pushed>=0) {
      glViewport(I->ViewPort[0],I->ViewPort[1],I->ViewPort[2],I->ViewPort[3]);
      glPopMatrix();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      
      I->Pushed--;
    }
  }
}

int OrthoGetPushed(PyMOLGlobals *G)
{
  return G->Ortho->Pushed;
}
/*========================================================================*/
void OrthoCommandIn(PyMOLGlobals *G,char *buffer)
{
  register COrtho *I=G->Ortho;
  if(I->cmdActiveQueue)
	QueueStrIn(I->cmdActiveQueue,buffer);
}
/*========================================================================*/
void OrthoPasteIn(PyMOLGlobals *G,char *buffer)
{
  register COrtho *I=G->Ortho;
  int curLine = I->CurLine&OrthoSaveLines;
  int execFlag=false;
  OrthoLineType buf2;

  if(I->InputFlag) {
    if(I->CursorChar>=0) {
      strcpy(buf2,I->Line[curLine]+I->CursorChar);
      strcpy(I->Line[curLine]+I->CursorChar,buffer);
      I->CurChar = strlen(I->Line[curLine]);
      I->CursorChar = I->CurChar;
      while((I->Line[curLine][I->CurChar-1]==10)||(I->Line[curLine][I->CurChar-1]==13)) 
        {
          execFlag=true;
          I->CurChar--;
          I->Line[curLine][I->CurChar]=0;
          if(I->CurChar<=I->PromptChar)
            break;
        }
      if(!execFlag) {
        strcpy(I->Line[curLine]+I->CursorChar,buf2);
        I->CurChar=strlen(I->Line[curLine]);
      }
    } else {
      strcat(I->Line[curLine],buffer);
      I->CurChar=strlen(I->Line[curLine]);
      while((I->Line[curLine][I->CurChar-1]==10)||(I->Line[curLine][I->CurChar-1]==13)) 
        {
          execFlag=true;
          I->CurChar--;
          I->Line[curLine][I->CurChar]=0;
          if(I->CurChar<=I->PromptChar)
            break;
        }
    }
  } else {
    OrthoRestorePrompt(G);
    
    while((I->Line[curLine][I->CurChar-1]==10)||(I->Line[curLine][I->CurChar-1]==13)) 
      {
        execFlag=true;
        I->CurChar--;
        I->Line[curLine][I->CurChar]=0;
        if(I->CurChar<=I->PromptChar)
          break;
      }
  }
  if(execFlag) {
    printf("[%s]\n",I->Line[curLine]);
    OrthoParseCurrentLine(G);
    } else
    I->InputFlag=true;
}








