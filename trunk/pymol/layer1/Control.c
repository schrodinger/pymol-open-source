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

#include "main.h"
#include "Base.h"
#include "Control.h"
#include "Scene.h"
#include "Executive.h"
#include "Setting.h"
#include "P.h"

#include"Movie.h"

#define cControlBoxSize 17
#define cControlLeftMargin 4
#define cControlTopMargin 5
#define cControlSpacing 2
#define cControlInnerMargin 4
#define cControlSpread 6

CControl Control;

void ControlDraw(Block *block);
int ControlClick(Block *block,int button,int x,int y,int mod);

Block *ControlGetBlock(void)
{
  CControl *I=&Control;
  {return(I->Block);}
}
/*========================================================================*/
int ControlIdling(void)
{
  CControl *I=&Control;
  return(MoviePlaying()||I->Rocking||SettingGet(cSetting_sculpting));
}
/*========================================================================*/
void ControlInterrupt(void)
{
  /*  CControl *I=&Control;*/
  MoviePlay(cMovieStop);
  ExecutiveDrawNow();
}
/*========================================================================*/
void ControlInit(void)
{
  CControl *I=&Control;

  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick = ControlClick;
  I->Block->fDraw    = ControlDraw;
  I->Block->fReshape = BlockReshape;
  I->Block->active = true;

  I->Block->TextColor[0]=1.0;
  I->Block->TextColor[1]=0.2;
  I->Block->TextColor[2]=0.2;

  OrthoAttach(I->Block,cOrthoTool);

  I->Rocking=false;

}
/*========================================================================*/
void ControlFree(void)
{
  CControl *I=&Control;
  OrthoFreeBlock(I->Block);
}
/*========================================================================*/
void ControlRock(int mode)
{
  CControl *I=&Control;
  switch(mode) {
  case -1:
	I->Rocking=!I->Rocking;
	break;
  case 0:
	I->Rocking=false;
	break;
  case 1:
	I->Rocking=true;
	break;
  }
  SceneRestartTimers();
  OrthoDirty();
}

/*========================================================================*/
int ControlClick(Block *block,int button,int x,int y,int mod)
{
  CControl *I=&Control;
  int sel = 0;
  int flag = false;
	
  x -= I->Block->rect.left+cControlLeftMargin;
  y -= I->Block->rect.top-cControlTopMargin;
  if((y<=0)&&(y>(-cControlBoxSize)))
	 while(x>=0) {
		if(x<cControlBoxSize)
		  {
			 flag=true;
			 break;
		  }
		x-=cControlBoxSize+cControlSpacing;
		sel++;
	 }
  if(flag) {
	 switch(sel) {
	 case 0:
      SceneSetFrame(4,0);
      PLog("cmd.rewind()",cPLog_pym);
		break;
    case 1:
      SceneSetFrame(1,-1);
      PLog("cmd.back()",cPLog_pym);
      break;
	 case 2:
		MoviePlay(cMovieStop);
      if(SettingGet(cSetting_sculpting)) SettingSet(cSetting_sculpting,0);
      if(I->Rocking) I->Rocking=false;
		ExecutiveDrawNow();
      OrthoDirty();
      PLog("cmd.mstop()",cPLog_pym);
		break;
	 case 3:
      if(!MoviePlaying()) {
        if(mod&cOrthoCTRL) {
          PLog("cmd.rewind()",cPLog_pym);
          PLog("cmd.mplay()",cPLog_pym);
          SceneSetFrame(0,0);		
          MoviePlay(cMoviePlay);
        } else {
          PLog("cmd.mplay()",cPLog_pym);
          MoviePlay(cMoviePlay);
        }
      } else {
        MoviePlay(cMovieStop);
        ExecutiveDrawNow();
        OrthoDirty();
        PLog("cmd.mstop()",cPLog_pym);
      }
		break;
    case 4:
      SceneSetFrame(1,1);
      PLog("cmd.forward()",cPLog_pym);
      break;
	 case 5:
		if(mod&cOrthoCTRL) {
		  SceneSetFrame(3,0);
        PLog("cmd.middle()",cPLog_pym);
		} else {
		  SceneSetFrame(2,0);
        PLog("cmd.ending()",cPLog_pym);
		}
		break;
    case 6:
      if(SettingGet(cSetting_sculpting)) {
        SettingSet(cSetting_sculpting,0.0);
        PLog("cmd.set('sculpting',0)",cPLog_pym);
      } else {
        SettingSet(cSetting_sculpting,1.0);
        PLog("cmd.set('sculpting',1)",cPLog_pym);        
      }
      OrthoDirty();
      break;
	 case 7:
		I->Rocking=!I->Rocking;
      if(I->Rocking)
        PLog("cmd.set('rocking',1)",cPLog_pym);
      else
        PLog("cmd.set('rocking',0)",cPLog_pym);
		SceneRestartTimers();
      OrthoDirty();
		break;
	 }
  }
  return(1);
}
/*========================================================================*/
void ControlDraw(Block *block)
{
  CControl *I=&Control;
  int x,y;

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    glColor3fv(I->Block->TextColor);
    
    x = I->Block->rect.left+cControlLeftMargin;
    y = I->Block->rect.top-cControlTopMargin;
    
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(x,y);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize/2));  
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize/2));  
    glEnd();
    x+=cControlBoxSize+cControlSpacing;
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(x,y);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize/2));  
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing;


    glBegin(GL_LINE_LOOP);
    glVertex2i(x,y);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing;

    if(MoviePlaying()) {
      glBegin(GL_TRIANGLE_STRIP);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize,y+1);
      glVertex2i(x+cControlBoxSize,y-(cControlBoxSize-1));

      glEnd();
      glColor3fv(I->Block->BackColor);
      glBegin(GL_LINE_LOOP);
      glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin+1);
      glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin-1);
      glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                 y-(cControlBoxSize/2));  
      glEnd();
      glColor3fv(I->Block->TextColor);
    } else {
      glBegin(GL_LINE_LOOP);
      glVertex2i(x,y);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y);
      glEnd();
      glBegin(GL_LINE_LOOP);
      glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
      glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
                 y-(cControlBoxSize/2));  
      glEnd();
    }
    x+=cControlBoxSize+cControlSpacing;

    glBegin(GL_LINE_LOOP);
    glVertex2i(x,y);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize/2));  
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing;
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(x,y);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize/2));  
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-(cControlBoxSize/2));  
    glEnd();
    x+=cControlBoxSize+cControlSpacing;

    if(SettingGet(cSetting_sculpting)) {
      glBegin(GL_TRIANGLE_STRIP);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize,y+1);
      glVertex2i(x+cControlBoxSize,y-(cControlBoxSize-1));
      
      glEnd();
      glColor3fv(I->Block->BackColor);
      glBegin(GL_LINE_STRIP);
      glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                 y-cControlInnerMargin+1);
      glVertex2i(x+cControlInnerMargin,
                 y-(cControlBoxSize/3));  
      glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                 y-(2*cControlBoxSize/3)+cControlInnerMargin-2);  
      glVertex2i(x+cControlInnerMargin,
                 y-(cControlBoxSize-1)+cControlInnerMargin-1);
      
      glEnd();
      glColor3fv(I->Block->TextColor);    
    } else {

      glBegin(GL_LINE_LOOP);
      glVertex2i(x,y);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y);
      glEnd();
      glBegin(GL_LINE_LOOP);
      glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                 y-cControlInnerMargin+1);
      glVertex2i(x+cControlInnerMargin,
                 y-(cControlBoxSize/3));  
      glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                 y-(2*cControlBoxSize/3)+cControlInnerMargin-2);  
      glVertex2i(x+cControlInnerMargin,
                 y-(cControlBoxSize-1)+cControlInnerMargin-1);
      glEnd();
      

    }
    x+=cControlBoxSize+cControlSpacing; 
    
    if(I->Rocking) {
      
      glBegin(GL_TRIANGLE_STRIP);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y+1);
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      
      glEnd();
      glColor3fv(I->Block->BackColor);
      glBegin(GL_LINE_LOOP);
      glVertex2i(x+(cControlBoxSize/2)+cControlSpread,
                 y-cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2),
                 y-(cControlBoxSize)+cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2)-cControlSpread,
                 y-cControlInnerMargin);
      glEnd();
      glColor3fv(I->Block->TextColor);
    } else {
      glBegin(GL_LINE_LOOP);
      glVertex2i(x,y);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y);
      glEnd();
      glBegin(GL_LINE_LOOP);
      glVertex2i(x+(cControlBoxSize/2)+cControlSpread,
                 y-cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2),
                 y-(cControlBoxSize)+cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2)-cControlSpread,
                 y-cControlInnerMargin);
      glEnd();
      
    }
    x+=cControlBoxSize+cControlSpacing; 



  }
}



