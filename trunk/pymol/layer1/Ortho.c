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
#include"Control.h"
#include"Setting.h"
#include"Wizard.h"
#include"Queue.h"
#include"Pop.h"

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

ListVarDeclare(BlockList,Block);

#define OrthoSaveLines 0xFF
#define OrthoHistoryLines 0xFF

#define cOrthoCharWidth 8
#define cOrthoLeftMargin 8
#define cOrthoBottomMargin 10

#define WizardMargin 119

#define ButModeMargin 26
#define ControlMargin 0

typedef struct {
  Block *Blocks;
  Block *GrabbedBy,*ClickedIn;
  Block LoopBlock;
  GLint ViewPort[4];
  int X,Y,Height,Width;
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
  float BusyLast;
  int BusyStatus[4];
  char BusyMessage[255];
  char *WizardPromptVLA;
  int SplashFlag;
  int LoopFlag;
  int LoopMod;
  BlockRect LoopRect;
  CQueue *cmds;
  CQueue *feedback;
} OrthoObject;

static OrthoObject Ortho;
void OrthoParseCurrentLine(void);
static void OrthoDrawWizardPrompt(void);

Block *OrthoFindBlock(int x,int y);
void OrthoKeyControl(unsigned char k);
void OrthoKeyAlt(unsigned char k);

#define cBusyWidth 240
#define cBusyHeight 60
#define cBusyMargin 10
#define cBusyBar 10
#define cBusySpacing 15

#define cBusyUpdate 0.2

#define cWizardTopMargin 15
#define cWizardLeftMargin 15
#define cWizardBorder 7
int OrthoLoopBlockDrag(Block *block,int x,int y,int mod);
int OrthoLoopBlockRelease(Block *block,int button,int x,int y,int mod);

int OrthoGetWidth(void)
{
  OrthoObject *I=&Ortho;
  return(I->Width);
}
/*========================================================================*/
int OrthoLoopBlockDrag(Block *block,int x,int y,int mod)
{
  OrthoObject *I=&Ortho;  
  I->LoopRect.right=x;
  I->LoopRect.bottom=y;
  OrthoDirty();
  return(1);
}
/*========================================================================*/
int OrthoLoopBlockRelease(Block *block,int button,int x,int y,int mod)
{
  OrthoObject *I=&Ortho;
  int tmp;
  int mode;
  mode = ButModeTranslate(button,I->LoopMod);

  if(I->LoopRect.top<I->LoopRect.bottom) {
    tmp=I->LoopRect.top;
    I->LoopRect.top=I->LoopRect.bottom;
    I->LoopRect.bottom=tmp;
  }
  if(I->LoopRect.right<I->LoopRect.left) {
    tmp=I->LoopRect.right;
    I->LoopRect.right=I->LoopRect.left;
    I->LoopRect.left=tmp;
  }
  ExecutiveSelectRect(&I->LoopRect,mode);
  I->LoopFlag=false;
  I->GrabbedBy=NULL;
  OrthoDirty();
  return(1);
}
/*========================================================================*/
void OrthoSetWizardPrompt(char *vla)
{
  OrthoObject *I=&Ortho;
  VLAFreeP(I->WizardPromptVLA);
  I->WizardPromptVLA=vla;
}
/*========================================================================*/

void OrthoSpecial(int k,int x,int y)
{
  OrthoObject *I=&Ortho;
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
  OrthoDirty();
}
/*========================================================================*/
int OrthoTextVisible(void) {
  return(SettingGet(cSetting_internal_feedback)||
         SettingGet(cSetting_text)||
         SettingGet(cSetting_overlay));
}
/*========================================================================*/

int OrthoArrowsGrabbed(void)
{
  OrthoObject *I=&Ortho;
  return(I->CurChar>I->PromptChar&&OrthoTextVisible()); 
  /* arrows can't be grabbed if text isn't visible */
}
/*========================================================================*/
void  OrthoRemoveSplash(void)
{
  OrthoObject *I=&Ortho;
  I->SplashFlag=false;
}
/*========================================================================*/
int  OrthoCommandOut(char *buffer)
{
  OrthoObject *I=&Ortho;
  if(I->cmds)
	return(QueueStrOut(I->cmds,buffer));
  else
	return(0);
}
/*========================================================================*/
int  OrthoCommandWaiting(void)
{
  OrthoObject *I=&Ortho;
  return(QueueStrCheck(I->cmds));
}
/*========================================================================*/
void OrthoClear(void)
{
  int a;
  OrthoObject *I=&Ortho;
  for(a=0;a<=OrthoSaveLines;a++)
    I->Line[a][0]=0;
  OrthoNewLine(NULL,true);
  OrthoRestorePrompt();
  OrthoDirty();
}
/*========================================================================*/
void OrthoFeedbackIn(char *buffer)
{
  OrthoObject *I=&Ortho;
  if(PMGUI) {
    if(I->feedback)
      QueueStrIn(I->feedback,buffer);
  }
}
/*========================================================================*/
int OrthoFeedbackOut(char *buffer)
{
  OrthoObject *I=&Ortho;
  if(I->feedback)
	return(QueueStrOut(I->feedback,buffer));
  else
	return(0);
}
/*========================================================================*/
void OrthoDirty(void) {
  OrthoObject *I=&Ortho;
  PRINTFD(FB_Ortho)
    " OrthoDirty: called.\n"
    ENDFD;
  if(!I->DirtyFlag) {
	 I->DirtyFlag = true;
  }

  MainDirty();
}
/*========================================================================*/
void OrthoBusyMessage(char *message)
{
  OrthoObject *I=&Ortho;
  if(strlen(message)<255)
	 strcpy(I->BusyMessage,message);
}
/*========================================================================*/
void OrthoBusySlow(int progress,int total)
{
  OrthoObject *I=&Ortho;
  PRINTFD(FB_Ortho)
    " OrthoBusySlow-DEBUG: progress %d total %d\n",progress,total
    ENDFD;
  I->BusyStatus[0]=progress;
  I->BusyStatus[1]=total;
  OrthoBusyDraw(false);
}
/*========================================================================*/
void OrthoBusyFast(int progress,int total)
{
  OrthoObject *I=&Ortho;
  PRINTFD(FB_Ortho)
    " OrthoBusyFast-DEBUG: progress %d total %d\n",progress,total
    ENDFD;
  I->BusyStatus[2]=progress;
  I->BusyStatus[3]=total;
  OrthoBusyDraw(false);
}
/*========================================================================*/
void OrthoBusyPrime(void)
{
  OrthoObject *I=&Ortho;
  int a;
  for(a=0;a<4;a++)
	 I->BusyStatus[a]=0;
  I->BusyMessage[0]=0;
  I->BusyLast = UtilGetSeconds();
}
/*========================================================================*/
void OrthoBusyDraw(int force)
{
  OrthoObject *I=&Ortho;
  char *c;
  int x,y;
  float black[3] = {0,0,0};
  float white[3] = {1,1,1};

  float now;
  float busyTime;

  PRINTFD(FB_Ortho)
    " OrthoBusyDraw: entered.\n"
    ENDFD;
  now = UtilGetSeconds();
  busyTime = (-I->BusyLast) + now;
  if(SettingGet(cSetting_show_progress)&&(force||(busyTime>cBusyUpdate))) {
    
    if(PIsGlutThread()) {
      OrthoPushMatrix();
      
      if(PMGUI) {
        glDrawBuffer(GL_FRONT);
        glClear(GL_DEPTH_BUFFER_BIT);
        
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
          glRasterPos4d(cBusyMargin,y-(cBusySpacing/2),0.0,1.0);
          while(*c)
            p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(c++));
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
        
        glDrawBuffer(GL_BACK);
        glFlush();
      }
      OrthoPopMatrix();
      OrthoDirty();/* switched from SceneDirty */
    I->BusyLast=now;
    }    
  }
  PRINTFD(FB_Ortho)
    " OrthoBusyDraw: leaving...\n"
    ENDFD;

}
/*========================================================================*/
void OrthoRestorePrompt(void) 
{
  OrthoObject *I=&Ortho;
  int curLine;
  if(!I->InputFlag) 
	 {
	 if(I->Saved[0]) 
		{
		  if(I->CurChar) {
			 OrthoNewLine(NULL,true);
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
			 OrthoNewLine(I->Prompt,true);
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
void OrthoKeyControl(unsigned char k) {
  char buffer[OrthoLineLength];

  /* safer...*/

  sprintf(buffer,"cmd._ctrl(chr(%d))",k);
  PLog(buffer,cPLog_pym);
  PParse(buffer);
  PFlush();

  /*  PBlockAndUnlockAPI();
      sprintf(buffer,"cmd._ctrl('%c')",k+64);
      PRunString(buffer);
      PLockAPIAndUnblock(); */

}
/*========================================================================*/
void OrthoKeyAlt(unsigned char k) {
  char buffer[OrthoLineLength];

  /* safer...*/

  sprintf(buffer,"cmd._alt(chr(%d))",k);
  PLog(buffer,cPLog_pym);
  PParse(buffer);
  PFlush();

  /*  PBlockAndUnlockAPI();
      sprintf(buffer,"cmd._ctrl('%c')",k+64);
      PRunString(buffer);
      PLockAPIAndUnblock(); */

}
/*========================================================================*/
void OrthoKey(unsigned char k,int x,int y,int mod)
{
  OrthoObject *I=&Ortho;
  char buffer[OrthoLineLength];
  int curLine;

  PRINTFD(FB_Ortho)
    " Ortho: %c (%d), x %d y %d, mod %d\n",k,k,x,y,mod
    ENDFD;
    
  if(!I->InputFlag) 
	 {
	 if(I->Saved[0]) 
		{
		  if(I->CurChar) {
			 OrthoNewLine(NULL,true);
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
			 OrthoNewLine(I->Prompt,true);
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
    OrthoKeyAlt(k);
  } else if((k>=32)&&(k!=127))
	 {
      curLine=I->CurLine&OrthoSaveLines;
      
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
	 }
  else switch(k)
	 {
    case 127: /* delete */     
#ifndef _PYMOL_OSX
      if((!I->CurChar)||(I->CurChar==I->PromptChar)||!OrthoTextVisible()) {
        OrthoKeyControl(4+64);
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
      if(OrthoArrowsGrabbed()) {
        I->CursorChar=-1;
      } else 
        OrthoKeyControl((unsigned char)(k+64));
      break;
    case 1: /* CTRL A -- beginning */
      if(OrthoArrowsGrabbed()) {
        if(I->CurChar)
          I->CursorChar=I->PromptChar;        
      } else 
        OrthoKeyControl((unsigned char)(k+64));
      break;
    case 4: /* CTRL D */
      if((!I->CurChar)||(I->CurChar==I->PromptChar)||!OrthoTextVisible()) {
        OrthoKeyControl((unsigned char)(4+64));
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
          PComplete(buffer+I->PromptChar,
                    sizeof(OrthoLineType)-I->PromptChar); /* just print, don't complete */
        }
      }
      break;
	 case 9: /* CTRL I -- tab */
      if(mod&cOrthoCTRL) {
        OrthoKeyControl((unsigned char)(k+64)); 
      } else {
        curLine=I->CurLine&OrthoSaveLines;
        if(I->PromptChar) {
          strcpy(buffer,I->Line[curLine]);
          
          if(PComplete(buffer+I->PromptChar,
                       sizeof(OrthoLineType)-I->PromptChar))
            {
              OrthoRestorePrompt();
              curLine=I->CurLine&OrthoSaveLines;
              strcpy(I->Line[curLine],buffer);
              I->CurChar = strlen(I->Line[curLine]);
            }
        }
      }
      break;
    case 27: /* ESCAPE */
      if(I->SplashFlag) {
        OrthoRemoveSplash();
      } else {
        SettingSet(cSetting_text,(float)(!((int)SettingGet(cSetting_text))));
        if(mod&cOrthoSHIFT) 
          SettingSet(cSetting_overlay,(float)(!((int)SettingGet(cSetting_overlay))));
      }
		break;
	 case 13: /* CTRL M -- carriage return */
      OrthoParseCurrentLine();
		break;
	 case 11: /* CTRL K -- truncate */
      if(OrthoArrowsGrabbed()) {
        if(I->CursorChar>=0) { 
          I->Line[I->CurLine&OrthoSaveLines][I->CursorChar]=0;
          I->CurChar=I->CursorChar;
          I->CursorChar=-1;
        }
      } else {
        if(mod&cOrthoCTRL) {
          OrthoKeyControl((unsigned char)(k+64));
        }
      }
      break;
    case 22: /* CTRL V -- paste */
      PParse("cmd.paste()");
      PFlush();
      /* PBlockAndUnlockAPI();
        PRunString("cmd.paste()");
        PLockAPIAndUnblock(); */
      break;
	 default:
      OrthoKeyControl((unsigned char)(k+64));
		break;
	 }
  MainDirty();
}
/*========================================================================*/
void OrthoParseCurrentLine(void) 
{
  OrthoObject *I=&Ortho;
  char buffer[OrthoLineLength];
  int curLine;
  curLine=I->CurLine&OrthoSaveLines;
  I->Line[curLine][I->CurChar]=0;
  strcpy(buffer,I->Line[curLine]+I->PromptChar);
  if(buffer[0])
    {
      strcpy(I->History[I->HistoryLine],buffer);
      I->HistoryLine = (I->HistoryLine+1)&OrthoHistoryLines;
      I->History[I->HistoryLine][0]=0;
      I->HistoryView=I->HistoryLine;
      if(WordMatch(buffer,"quit",true)==0) /* don't log quit */
        PLog(buffer,cPLog_pml);
      OrthoNewLine(NULL,true);
      ExecutiveDrawNow();
      PParse(buffer);
      OrthoRestorePrompt();
    }
  I->CursorChar=-1;
}
/*========================================================================*/
void OrthoAddOutput(char *str)
{
  OrthoObject *I=&Ortho;
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
          wrap = (int)SettingGet(cSetting_wrap_output);

          if(wrap>0) {
            if(cc>wrap)
              {
                *q=0;
                I->CurChar = cc;
                OrthoNewLine(NULL,true);
                cc=0;
                q=I->Line[I->CurLine&OrthoSaveLines];
                curLine = I->CurLine&OrthoSaveLines;
              }
          } 
          if(cc>=OrthoLineLength-6) { /* fail safe */
            *q=0;
            I->CurChar = cc;
            OrthoNewLine(NULL,false);
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
			 OrthoNewLine(NULL,true);
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
  if((SettingGet(cSetting_internal_feedback)>1)||SettingGet(cSetting_overlay))
    OrthoDirty();
}
/*========================================================================*/
void OrthoNewLine(char *prompt,int crlf)
{
  int curLine;
  OrthoObject *I=&Ortho;

  /*  printf("orthoNewLine: CC: %d CL:%d PC: %d IF:L %d\n",I->CurChar,I->CurLine,
		I->PromptChar,I->InputFlag);*/
  /*  if(I->CurChar)
	 {*/
  if(I->CurChar)
    OrthoFeedbackIn(I->Line[I->CurLine&OrthoSaveLines]);
  else
    OrthoFeedbackIn(" ");
  if(Feedback(FB_Python,FB_Output)) {
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
void OrthoGrab(Block *block)
{
  OrthoObject *I=&Ortho;
  I->GrabbedBy = block;
}
/*========================================================================*/
void OrthoUngrab(void)
{
  OrthoObject *I=&Ortho;
  I->GrabbedBy = NULL;
}
/*========================================================================*/
Block *OrthoNewBlock(Block *block)
{
  if(!block)
	 ListElemAlloc(block,Block);
  UtilZeroMem(block,sizeof(Block));
  BlockInit(block);
  return(block);
}
/*========================================================================*/
void OrthoFreeBlock(Block *block)
{
  if(block) 
    ListElemFree(block);
}
/*========================================================================*/
void OrthoAttach(Block *block,int type)
{
  OrthoObject *I=&Ortho;
  ListInsert(I->Blocks,block,NULL,next,BlockList);
}
/*========================================================================*/
void OrthoDetach(Block *block)
{
  OrthoObject *I=&Ortho;
  ListDetach(I->Blocks,block,next,BlockList);
}
/*========================================================================*/
void OrthoDoDraw()
{
  OrthoObject *I=&Ortho;

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
  
  int render = false;
  double_pump=SettingGet_i(NULL,NULL,cSetting_stereo_double_pump_mono);
  bg_color=SettingGet_3fv(NULL,NULL,cSetting_bg_rgb);

  I->OverlayColor[0]=1.0-bg_color[0];
  I->OverlayColor[1]=1.0-bg_color[1];
  I->OverlayColor[2]=1.0-bg_color[2];
  if(diff3f(I->OverlayColor,bg_color)<0.25)
    zero3f(I->OverlayColor);

  PRINTFD(FB_Ortho)
    " OrthoDoDraw: entered.\n"
    ENDFD;
  if(PMGUI) {

    rightSceneMargin=(int)SettingGet(cSetting_internal_gui_width);
    internal_feedback=(int)SettingGet(cSetting_internal_feedback);

    v=SettingGetfv(cSetting_bg_rgb);
    overlay = (int)SettingGet(cSetting_overlay);
    text = (int)SettingGet(cSetting_text);

    if(text) overlay=0;
    
    glClearColor(v[0],v[1],v[2],1.0);

    if(overlay||(!text)) 
      if(!SceneRenderCached())
        render=true;

    if(SceneGetStereo()||double_pump) {
      glDrawBuffer(GL_BACK_LEFT);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glDrawBuffer(GL_BACK_RIGHT);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glClearColor(0.0,0.0,0.0,1.0);
      times = 2;
    } else {
      glDrawBuffer(GL_BACK);
      glClearColor(v[0],v[1],v[2],1.0);
      glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glClearColor(0.0,0.0,0.0,1.0);
      times = 1;
      double_pump=false;
    }    

    if(render)
      SceneRender(NULL,0,0,NULL);

    while(times--) {
      switch(times) {
      case 1:
          glDrawBuffer(GL_BACK_RIGHT);
        break;
      case 0:
        if(double_pump) {
          glDrawBuffer(GL_BACK_LEFT);
        } else
          glDrawBuffer(GL_BACK);
        break;
      }
      OrthoPushMatrix();
      
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
      
      PRINTFD(FB_Ortho)
        " OrthoDoDraw: drawing blocks...\n"
        ENDFD;
      
      BlockRecursiveDraw(I->Blocks);
      
      PRINTFD(FB_Ortho)
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
      
      OrthoRestorePrompt();
      
      if(I->DrawText) {	 
        /* now print the text */
        
        lcount = 0;
        l=(I->CurLine-lcount)&OrthoSaveLines;
        x = cOrthoLeftMargin;
        y = cOrthoBottomMargin;
        
        
        if((int)SettingGet(cSetting_text)||I->SplashFlag)
          showLines=I->ShowLines;
        else
          showLines=internal_feedback+(int)SettingGet(cSetting_overlay);
        
        glColor3fv(I->TextColor);
        while(l>=0)
          {
            lcount++;
            if(lcount>showLines)
              break;
            str = I->Line[l&OrthoSaveLines];
            if(strncmp(str,I->Prompt,6)==0)
              glColor3fv(I->TextColor);            
            else
              glColor3fv(I->OverlayColor);
            glRasterPos4d((double)x,(double)y,0.0,1.0);
            if(str)
              {
                while(*str)
                  p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*(str++));
                if((lcount==1)&&(I->InputFlag)) 
                  {
                    if(I->CursorChar>=0)  
                      glRasterPos4d((double)(x+8*I->CursorChar),(double)y,0.0,1.0);
                    p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,'_');
                  }
              }
            l=(I->CurLine-lcount)&OrthoSaveLines;
            y=y+cOrthoLineHeight;
          }
      }
      
      OrthoDrawWizardPrompt();
      OrthoPopMatrix();
    }
  }

  I->DirtyFlag =false;
  PRINTFD(FB_Ortho)
    " OrthoDoDraw: leaving...\n"
    ENDFD;

}
/*========================================================================*/
static void OrthoDrawWizardPrompt(void)
{
  /* assumes PMGUI */

  OrthoObject *I=&Ortho;
  
  char *vla,*p;
  int nLine;
  int x,y;
  int nChar,c,ll;
  int maxLen;
  BlockRect rect;

  if(I->WizardPromptVLA) {
    vla = I->WizardPromptVLA;
    
    nLine = UtilCountStringVLA(vla);
    if(nLine) {
      nChar = VLAGetSize(I->WizardPromptVLA);
      
      /* count max line length */
      
      maxLen = 0;
      p = vla;
      ll = 0;
      c=nChar;
      while(c--) {
        if(!*(p++)) {
          if(maxLen<ll)
            maxLen = ll;
          ll=0;
        } else 
          ll++;
      }
      
      rect.top = I->Height-cWizardTopMargin;
      rect.bottom = rect.top-(nLine*cOrthoLineHeight+2*cWizardBorder)-2;
      rect.left = cWizardLeftMargin;
      rect.right = rect.left + cOrthoCharWidth*maxLen + 2*cWizardBorder+1;
      
      glColor3fv(I->WizardBackColor);

      glBegin(GL_POLYGON);
      glVertex2i(rect.right,rect.top);
      glVertex2i(rect.right,rect.bottom);
      glVertex2i(rect.left,rect.bottom);
      glVertex2i(rect.left,rect.top);
      glEnd();
      
      glColor3fv(I->WizardTextColor);
      
      x = rect.left+cWizardBorder;
      y = rect.top-(cWizardBorder+cOrthoLineHeight);

      vla = I->WizardPromptVLA;
      
      /* count max line length */
      
      glRasterPos4d((double)x,(double)y,0.0,1.0);
      p = vla;
      ll = 0;
      c=nChar;
      while(c--) {
        if(*p)
          p_glutBitmapCharacter(P_GLUT_BITMAP_8_BY_13,*p);
        if(!*(p++)) {
          y=y-cOrthoLineHeight;
          glRasterPos4d((double)x,(double)y,0.0,1.0);          
        }
      }
    }
  }
}

/*========================================================================*/
void OrthoReshape(int width, int height)
{
  OrthoObject *I=&Ortho;

  Block *block = NULL;
  int sceneBottom,sceneRight;
  int internal_gui_width;
  int internal_feedback;

  PRINTFD(FB_Ortho)
    " OrthoReshape-Debug: %d %d\n",width,height
    ENDFD;

  if((width!=I->Width)||(height!=I->Height)) {
  if(width<0) width=I->Width;
  if(height<0) height=I->Height;

  I->Height=height;
  I->Width=width;
  I->ShowLines = height/cOrthoLineHeight;
  
  internal_feedback = (int)SettingGet(cSetting_internal_feedback);
  if(internal_feedback)
    sceneBottom = (internal_feedback-1)*cOrthoLineHeight + cOrthoBottomSceneMargin;
  else
    sceneBottom = 0;
    
  internal_gui_width = (int)SettingGet(cSetting_internal_gui_width);
  if(SettingGet(cSetting_internal_gui)>0.0) {
    sceneRight = internal_gui_width;
  } else {
    sceneRight = 0;
    internal_gui_width=0;
  }
  block=SceneGetBlock();
  BlockSetMargin(block,0,0,sceneBottom,sceneRight);
  BlockSetMargin(&I->LoopBlock,0,0,sceneBottom,sceneRight);

  if(SettingGet(cSetting_internal_gui)) {
    block=ExecutiveGetBlock();
    block->active=true;
    BlockSetMargin(block,0,width-internal_gui_width,WizardMargin,0);
    block=WizardGetBlock();
    BlockSetMargin(block,height-WizardMargin,width-internal_gui_width,WizardMargin,0);
    block->active=false;
    block=ButModeGetBlock();
    BlockSetMargin(block,height-WizardMargin,width-internal_gui_width,ButModeMargin,0);
    block->active=true;
    block=ControlGetBlock();
    BlockSetMargin(block,height-ButModeMargin,width-internal_gui_width,ControlMargin,0);
    block->active=true;
  } else {
    block=ExecutiveGetBlock();
    block->active=false;
    BlockSetMargin(block,0,width-internal_gui_width,WizardMargin,0);
    block=WizardGetBlock();
    BlockSetMargin(block,height-WizardMargin,width-internal_gui_width,WizardMargin,0);
    block->active=false;
    block=ButModeGetBlock();
    BlockSetMargin(block,height-WizardMargin,width-internal_gui_width,ButModeMargin,0);
    block->active=false;
    block=ControlGetBlock();
    BlockSetMargin(block,height-ButModeMargin,width-internal_gui_width,ControlMargin,0);
    block->active=false;
  }

  if(PMGUI) 
    glGetIntegerv(GL_VIEWPORT,I->ViewPort);

  OrthoPushMatrix();
  block=NULL;
  while(ListIterate(I->Blocks,block,next))
	 if(block->fReshape)
		block->fReshape(block,width,height);			
  OrthoPopMatrix();

  WizardRefresh(); /* safe to call even if no wizard exists */
  }
}

/*========================================================================*/
void OrthoReshapeWizard(int wizHeight)
{
  Block *block;
  OrthoObject *I=&Ortho;
  int height,width;
  int internal_gui_width;

  height=I->Height;
  width=I->Width;

  if(SettingGet(cSetting_internal_gui)>0.0) {
    internal_gui_width = (int)SettingGet(cSetting_internal_gui_width);
    block=ExecutiveGetBlock();
    if(height) {
      BlockSetMargin(block,0,width-internal_gui_width,WizardMargin+wizHeight,0);
    } else {
      BlockSetMargin(block,0,width-internal_gui_width,WizardMargin,0);
    }
    block->fReshape(block,width,height);

    block=WizardGetBlock();

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
Block *OrthoFindBlock(int x,int y)
{
  OrthoObject *I=&Ortho;

  return(BlockRecursiveFind(I->Blocks,x,y));
}
/*========================================================================*/
int OrthoButton(int button,int state,int x,int y,int mod)
{
  OrthoObject *I=&Ortho;

  Block *block=NULL;
  int handled = 0; 

  OrthoRemoveSplash();
  I->X=x;
  I->Y=y;

  if(state==P_GLUT_DOWN)
	 {
		I->ActiveButton = button;
		if(I->GrabbedBy)
		  {
			 if(I->GrabbedBy->inside)
				block = BlockRecursiveFind(I->GrabbedBy->inside,x,y);
			 else
				block = I->GrabbedBy;
		  }
		else
		  block = OrthoFindBlock(x,y);
		if(block)
		  {
			 I->ClickedIn = block;
			 if(block->fClick)
				{
				  handled = block->fClick(block,button,x,y,mod);
				}
		  }
	 }
  else if(state==P_GLUT_UP)
	 {
      if(I->GrabbedBy)
        {
			 block=I->GrabbedBy;
			 if(block->fRelease)
            handled = block->fRelease(block,button,x,y,mod);
			 I->ClickedIn = NULL;
        }
		if(I->ClickedIn)
		  {
			 block=I->ClickedIn;
			 if(block->fRelease)
            handled = block->fRelease(block,button,x,y,mod);
			 I->ClickedIn = NULL;
		  }
	 }
  if(block&&!handled) {
    if(SceneGetBlock()==block) {
      if(state==P_GLUT_DOWN) {
        I->LoopRect.left=x;
        I->LoopRect.top=y;
        I->LoopRect.right=x;
        I->LoopRect.bottom=y;
        I->LoopFlag=true;
        I->LoopMod = mod;
        I->GrabbedBy=&I->LoopBlock;
        OrthoDirty();
      } 
    }
  }
  return(handled);
}
/*========================================================================*/ 
int OrthoDrag(int x, int y,int mod)
{
  OrthoObject *I=&Ortho;

  Block *block=NULL;
  int handled = 0;

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
void OrthoSplash(void) 
{
  OrthoNewLine(NULL,true);
  PRINTF " PyMOL(TM) Molecular Graphics System, Version " ENDF;
  PRINTF _PyMOL_VERSION ENDF;
  PRINTF ".\n" ENDF;
  PRINTF " Copyright (C) 1998-2003 by DeLano Scientific LLC.\n All Rights Reserved.\n \n" ENDF;

  PRINTF "    Created by Warren L. DeLano, Ph.D. \n \n" ENDF;

  PRINTF "    Other Major Authors and Contributors:\n\n" ENDF;
  PRINTF "       Ralf W. Grosse-Kunstleve, Ph.D.\n \n" ENDF;

  PRINTF "    PyMOL is user-supported open-source software.  Although most versions\n" ENDF;
  PRINTF "    are freely available, PyMOL is not in the public domain.\n \n" ENDF;

  PRINTF "    If PyMOL is helpful in your work or study, then please volunteer \n" ENDF;
  PRINTF "    support for our ongoing campaign to create open and affordable software \n" ENDF;
  PRINTF "    for molecular research.\n\n" ENDF;

  PRINTF "    Updates and other information can be found at \"http://www.pymol.org\".\n \n" ENDF;

  PRINTF "    Please cite PyMOL in publications and presentations:\n \n" ENDF;
  PRINTF "       Warren L. DeLano \"The PyMOL Molecular Graphics System.\"\n" ENDF;
  PRINTF "       DeLano Scientific LLC, San Carlos, CA, USA. http://www.pymol.org\n \n" ENDF;


  PRINTF "    Enter \"help\" for a list of commands.\n" ENDF;
  PRINTF 
    "    Enter \"help <command-name>\" for information on a specific command.\n\n"
    ENDF;
  
  PRINTF " Hit ESC anytime to toggle between text and graphics.\n\n" ENDF;
}
/*========================================================================*/
void OrthoInit(int showSplash)
{
  OrthoObject *I=&Ortho;
  int a;

  I->cmds = QueueNew(0xFFFF);
  I->feedback = QueueNew(0xFFFF);

  I->WizardBackColor[0]=0.2F;
  I->WizardBackColor[1]=0.2F;
  I->WizardBackColor[2]=0.2F;
  I->WizardTextColor[0]=0.2F;
  I->WizardTextColor[1]=1.0F;
  I->WizardTextColor[2]=0.2F;
  I->Blocks = NULL;
  I->GrabbedBy = NULL;
  I->ClickedIn = NULL;
  I->DrawText=1;
  I->TextColor[0]=0.7F;
  I->TextColor[1]=0.7F;
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
  I->LoopFlag=0;
  I->LoopMod=0;
  I->SplashFlag = false;
  I->ShowLines = 1;
  I->Saved[0]=0;
  I->DirtyFlag = true;
  BlockInit(&I->LoopBlock);
  I->LoopBlock.fDrag = OrthoLoopBlockDrag;
  I->LoopBlock.fRelease = OrthoLoopBlockRelease;
  
  if(showSplash) {
	 OrthoSplash();
    I->SplashFlag=true;
  }
  OrthoFeedbackIn(" ");
  I->CurLine++;
  strcpy(I->Prompt,"PyMOL>");
  strcpy(I->Line[I->CurLine],I->Prompt);
  I->CurChar = (I->PromptChar = strlen(I->Prompt));
  I->InputFlag=1;

  /*printf("orthoNewLine: CC: %d CL:%d PC: %d IF:L %d\n",I->CurChar,I->CurLine,
	 I->PromptChar,I->InputFlag);*/
  
  ButModeInit();
  ControlInit();
  PopInit();
  for(a=0;a<=OrthoHistoryLines;a++)
    I->History[a][0]=0;
}
/*========================================================================*/
void OrthoFree(void)
{
  OrthoObject *I=&Ortho;

  VLAFreeP(I->WizardPromptVLA);
  PopFree();
  ButModeFree();
  ControlFree();
  QueueFree(I->cmds);
  I->cmds=NULL;
  QueueFree(I->feedback);
  I->feedback=NULL;
}
/*========================================================================*/
void OrthoPushMatrix(void)
{
  OrthoObject *I=&Ortho;

  if(PMGUI) {
    glGetIntegerv(GL_VIEWPORT,I->ViewPort);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0,I->ViewPort[2],0,I->ViewPort[3],-100,100);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
  }
  /*  glDisable(GL_ALPHA_TEST);
  glDisable(GL_CULL_FACE);
  glDisable(GL_POINT_SMOOTH);*/
  
}
/*========================================================================*/
void OrthoPopMatrix(void)
{

  if(PMGUI) {
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
  }
}

/*========================================================================*/
void OrthoCommandIn(char *buffer)
{
  OrthoObject *I=&Ortho;
  if(I->cmds) 
	QueueStrIn(I->cmds,buffer);
}
/*========================================================================*/
void OrthoPasteIn(char *buffer)
{
  OrthoObject *I=&Ortho;
  int curLine = I->CurLine&OrthoSaveLines;
  int execFlag=false;
  OrthoLineType buf2;

  if(I->CurChar) {
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
    strcpy(I->Line[curLine]+I->PromptChar,buffer);
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
  if(execFlag)
    OrthoParseCurrentLine();
}








