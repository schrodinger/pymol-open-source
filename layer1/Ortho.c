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

#define cOrthoLineHeight 12
#define cOrthoLeftMargin 8
#define cOrthoBottomMargin 10

#define ExecutiveMargin 130
#define ButModeMargin 40
#define ControlMargin 0

typedef struct {
  Block *Blocks;
  Block *GrabbedBy,*ClickedIn;
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
  float TextColor[3],OverlayColor[3];
  int DirtyFlag;
  float BusyLast;
  int BusyStatus[4];
  char BusyMessage[255];
  int SplashFlag;
  CQueue *cmds;
  CQueue *feedback;
} OrthoObject;

static OrthoObject Ortho;
void OrthoParseCurrentLine(void);

Block *OrthoFindBlock(int x,int y);

#define cBusyWidth 180
#define cBusyHeight 60
#define cBusyMargin 10
#define cBusyBar 10
#define cBusySpacing 15

#define cBusyUpdate 1.0
/*========================================================================*/

void OrthoSpecial(int k,int x,int y)
{
  OrthoObject *I=&Ortho;
  int curLine = I->CurLine&OrthoSaveLines;
  switch(k) {
  case GLUT_KEY_DOWN:
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
  case GLUT_KEY_UP:
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
  case GLUT_KEY_LEFT:
    if(I->CursorChar>=0) {
      I->CursorChar--;
    } else {
      I->CursorChar = I->CurChar-1;
    }
    if(I->CursorChar<I->PromptChar)
      I->CursorChar=I->PromptChar;
    break;
  case GLUT_KEY_RIGHT:
    if(I->CursorChar>=0) {
      I->CursorChar++;
    } else {
      I->CursorChar = I->CurChar-1;
    }
    if(I->CursorChar>strlen(I->Line[curLine]))
      I->CursorChar=strlen(I->Line[curLine]);
    break;
  }
  OrthoDirty();
}
/*========================================================================*/

int OrthoArrowsGrabbed(void)
{
  OrthoObject *I=&Ortho;
  return(I->CurChar>I->PromptChar);
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
  OrthoNewLine(NULL);
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
  I->BusyStatus[0]=progress;
  I->BusyStatus[1]=total;
  OrthoBusyDraw(false);
}
/*========================================================================*/
void OrthoBusyFast(int progress,int total)
{
  OrthoObject *I=&Ortho;
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

  now = UtilGetSeconds();
  busyTime = (-I->BusyLast) + now;
  if(force||(busyTime>cBusyUpdate)) {

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
          glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
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
    }
	 OrthoPopMatrix();
	 SceneDirty();
	 I->BusyLast=now;
  }
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
			 OrthoNewLine(NULL);
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
			 OrthoNewLine(I->Prompt);
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
void OrthoKey(unsigned char k,int x,int y,int mod)
{
  OrthoObject *I=&Ortho;
  char buffer[OrthoLineLength];
  int curLine;
  
  if(!I->InputFlag) 
	 {
	 if(I->Saved[0]) 
		{
		  if(I->CurChar) {
			 OrthoNewLine(NULL);
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
			 OrthoNewLine(I->Prompt);
		  else
			 {
				curLine = I->CurLine&OrthoSaveLines;
				strcpy(I->Line[curLine],I->Prompt);
				I->CurChar = (I->PromptChar = strlen(I->Prompt));
			 }
		}
	 I->InputFlag=1;
  }
  if((k>=32)&&(k!=127))
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
	 case 4: /* ctrl D */
      if((!I->CurChar)||(I->CurChar==I->PromptChar))
        exit(0); 
      /* otherwise */
    case 127: /* delete */
      if(I->CursorChar>=0) {
        if(I->CursorChar<I->CurChar)
          I->CursorChar++;
        if(I->CursorChar==I->CurChar)
          I->CursorChar=-1;
      }
	 case 8:
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
      if(OrthoArrowsGrabbed())
        I->CursorChar=-1;
      break;
    case 1: /* CTRL A -- beginning */
      if(OrthoArrowsGrabbed())
        if(I->CurChar)
          I->CursorChar=I->PromptChar;        
      break;
	 case 9: /* CTRL I -- tab */
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
      if(OrthoArrowsGrabbed())
        if(I->CursorChar>=0) { 
          I->Line[I->CurLine&OrthoSaveLines][I->CursorChar]=0;
          I->CurChar=I->CursorChar;
          I->CursorChar=-1;
        }
      break;
    case 22: /* CTRL V -- paste */
      PBlockAndUnlockAPI();
      PRunString("cmd.paste()");
      PLockAPIAndUnblock();
      break;
	 default:
      PBlockAndUnlockAPI();
      sprintf(buffer,"cmd._ctrl('%c')",k+65);
      PRunString(buffer);
      PLockAPIAndUnblock();      
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
      OrthoNewLine(NULL);
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
			 if(cc>76)
				{
				  *q=0;
				  I->CurChar = cc;
				  OrthoNewLine(NULL);
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
			 OrthoNewLine(NULL);
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
  OrthoDirty();
}
/*========================================================================*/
void OrthoNewLine(char *prompt)
{
  int curLine;
  OrthoObject *I=&Ortho;

  /*  printf("orthoNewLine: CC: %d CL:%d PC: %d IF:L %d\n",I->CurChar,I->CurLine,
		I->PromptChar,I->InputFlag);*/
  if(I->CurChar)
	 {
	   OrthoFeedbackIn(I->Line[I->CurLine&OrthoSaveLines]);
	   printf("%s\n",I->Line[I->CurLine&OrthoSaveLines]);
	   fflush(stdout);
	 }

  if(I->Line[I->CurLine&OrthoSaveLines][0])
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
  int overlay,text;

  
  if(PMGUI) {

    v=SettingGetfv(cSetting_bg_rgb);
    overlay = SettingGet(cSetting_overlay);
    text = SettingGet(cSetting_text);

    glDrawBuffer(GL_BACK);
    glClearColor(v[0],v[1],v[2],1.0);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glClearColor(0.0,0.0,0.0,1.0);
    
    if(overlay||(!text)) 
      if(!SceneRenderCached())
        SceneRender(NULL,0,0);
    
    OrthoPushMatrix();
    
    x = I->X;
    y = I->Y;
    
    BlockRecursiveDraw(I->Blocks);
    
    OrthoRestorePrompt();

    if(I->DrawText) {	 
      /* now print the text */
	 
      lcount = 0;
      l=(I->CurLine-lcount)&OrthoSaveLines;
      x = cOrthoLeftMargin;
      y = cOrthoBottomMargin;

      glColor3f(0.0,0.0,0.0);
      glBegin(GL_POLYGON);
      glVertex2i(I->Width-cOrthoRightSceneMargin,cOrthoBottomSceneMargin-1);
      glVertex2i(I->Width-cOrthoRightSceneMargin,0);
      glVertex2i(0,0);
      glVertex2i(0,cOrthoBottomSceneMargin-1);
      glEnd();

      if((int)SettingGet(cSetting_text)||I->SplashFlag)
        showLines=I->ShowLines;
      else
        showLines=1;

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
                glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(str++));
              if((lcount==1)&&(I->InputFlag)) 
                {
                  if(I->CursorChar>=0)  
                    glRasterPos4d((double)(x+8*I->CursorChar),(double)y,0.0,1.0);
                  glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'_');
                }
            }
          l=(I->CurLine-lcount)&OrthoSaveLines;
          y=y+cOrthoLineHeight;
        }
    }

    OrthoPopMatrix();
  }

  I->DirtyFlag =false;
}
/*========================================================================*/
void OrthoReshape(int width, int height)
{
  OrthoObject *I=&Ortho;

  Block *block = NULL;
  int sceneBottom,sceneRight;

  I->Height=height;
  I->Width=width;
  I->ShowLines = height/cOrthoLineHeight;

  sceneBottom = cOrthoBottomSceneMargin;
    
  if(!InternalGUI) {
    sceneRight = 0;
  } else {
    sceneRight = cOrthoRightSceneMargin;
  }

  block=SceneGetBlock();
  BlockSetMargin(block,0,0,sceneBottom,sceneRight);

  if(InternalGUI) {
    block=ExecutiveGetBlock();
    BlockSetMargin(block,0,width-cOrthoRightSceneMargin,ExecutiveMargin,0);
    block=ButModeGetBlock();
    BlockSetMargin(block,height-ExecutiveMargin,width-cOrthoRightSceneMargin,ButModeMargin,0);
    block=ControlGetBlock();
    BlockSetMargin(block,height-ButModeMargin,width-cOrthoRightSceneMargin,ControlMargin,0);
  } else {
    block=ExecutiveGetBlock();
    block->active=false;
    block=ButModeGetBlock();
    block->active=false;
    block=ControlGetBlock();
    block->active=false;
  }
  
  glGetIntegerv(GL_VIEWPORT,I->ViewPort);

  OrthoPushMatrix();
  block=NULL;
  while(ListIterate(I->Blocks,block,next,BlockList))
	 if(block->fReshape)
		block->fReshape(block,width,height);			
  OrthoPopMatrix();
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

  Block *block;
  int handled = 0; 

  OrthoRemoveSplash();
  I->X=x;
  I->Y=y;

  if(state==GLUT_DOWN)
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
  else if(state==GLUT_UP)
	 {
      if(I->GrabbedBy)
        {
			 block=I->GrabbedBy;
			 if(block->fRelease)
            handled = block->fRelease(block,x,y,mod);
			 I->ClickedIn = NULL;
        }
		if(I->ClickedIn)
		  {
			 block=I->ClickedIn;
			 if(block->fRelease)
            handled = block->fRelease(block,x,y,mod);
			 I->ClickedIn = NULL;
		  }
	 }
  return(handled);
}
/*========================================================================*/ 
int OrthoDrag(int x, int y,int mod)
{
  OrthoObject *I=&Ortho;

  Block *block;
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
  OrthoNewLine(NULL);
  OrthoAddOutput("PyMOL Molecular Graphics System, Version ");
  OrthoAddOutput(_PyMOL_VERSION);
  OrthoAddOutput(".");
  OrthoNewLine(NULL);
  OrthoAddOutput("Copyright (C) 1998-2000 by DeLano Scientific.\nAll Rights Reserved.\n \n");
  OrthoAddOutput("Principle Author:  Warren L. DeLano, Ph.D.\n \n");
  OrthoAddOutput("Major Authors and Contributors:\n \n");
  OrthoAddOutput("      Ralf W. Grosse-Kunstleve, Ph.D.\n \n");
  OrthoAddOutput("This software is open source and freely available.\n");
  OrthoAddOutput("Updates can be found at \"http://www.pymol.org\".\n \n");
  OrthoAddOutput("Enter \"help commands\" for a list of commands.\n");
  OrthoAddOutput("Enter \"help <command-name>\" for information on a specific command.\n \n");
  OrthoAddOutput("Hit TAB to toggle text or type \"cls\" to clear.\n \n");
}
/*========================================================================*/
void OrthoInit(void)
{
  OrthoObject *I=&Ortho;
  int a;

  I->cmds = QueueNew(0xFFFF);
  I->feedback = QueueNew(0xFFFF);

  I->Blocks = NULL;
  I->GrabbedBy = NULL;
  I->ClickedIn = NULL;
  I->DrawText=1;
  I->TextColor[0]=0.7;
  I->TextColor[1]=0.7;
  I->TextColor[2]=1.0;
  I->OverlayColor[0]=1.0;
  I->OverlayColor[1]=1.0;
  I->OverlayColor[2]=1.0;
  I->CurLine=1000;
  I->PromptChar=0;
  I->CurChar=0;
  I->CursorChar=-1;
  I->HistoryLine=0;
  I->HistoryView=0;
  I->Line[I->CurLine&OrthoSaveLines][I->CurChar]=0;

  I->SplashFlag = true;
  I->ShowLines = 1;
  I->Saved[0]=0;
  I->DirtyFlag = true;

  OrthoSplash();
  strcpy(I->Prompt,"PyMOL>");
  OrthoNewLine(I->Prompt);

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
  /*  glDisable(GL_ALPHA_TEST);
  glDisable(GL_CULL_FACE);
  glDisable(GL_POINT_SMOOTH);*/
  
}
/*========================================================================*/
void OrthoPopMatrix(void)
{

  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
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








