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

#include "main.h"
#include "Base.h"
#include "Control.h"
#include "Scene.h"
#include "Executive.h"
#include "Setting.h"
#include "P.h"
#include "Grap.h"

#include"Movie.h"

#define cControlBoxSize 17
#define cControlLeftMargin 8
#define cControlTopMargin 5
#define cControlSpacing 2
#define cControlInnerMargin 4
#define cControlSpread 6
#define cControlSize 160

#define cControlButtons 7

CControl Control;

void ControlDraw(Block *block);
int ControlClick(Block *block,int button,int x,int y,int mod);

static void ControlReshape(Block *block,int width, int height)
{
  CControl *I=&Control;
  BlockReshape(block,width,height);

  I->ExtraSpace = ((block->rect.right-block->rect.left)-cControlSize);
  if(I->ExtraSpace<0)
    I->ExtraSpace=0;
  
}

static int ControlDrag(Block *block,int x,int y,int mod)
{
  int width;
  int delta;
  int gui_width;
  CControl *I=&Control;
  delta = x-I->LastPos;
  if(I->DragFlag) {
    if(delta) {
      gui_width = SettingGet(cSetting_internal_gui_width)-delta;
      if(gui_width<3)
        gui_width = 3;
      delta = SettingGet(cSetting_internal_gui_width)-gui_width;
      width = OrthoGetWidth()+delta;
    I->LastPos = x;
    SettingSet(cSetting_internal_gui_width,(float)gui_width);
    OrthoReshape(-1,-1);
    }
  }
  return(1);
}

static int ControlRelease(Block *block,int button,int x,int y,int mod)
{
  CControl *I=&Control;  
  OrthoUngrab();
  I->DragFlag=false;
  return(1);
}

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
  I->Block->fDrag = ControlDrag;
  I->Block->fRelease = ControlRelease;
  I->Block->fReshape = ControlReshape;
  I->Block->active = true;
  I->Block->TextColor[0]=1.0;
  I->Block->TextColor[1]=0.75;
  I->Block->TextColor[2]=0.75;
  I->ButtonColor[0]=0.5;
  I->ButtonColor[1]=0.5;
  I->ButtonColor[2]=0.5;
  I->ActiveColor[0]=0.8;
  I->ActiveColor[1]=0.8;
  I->ActiveColor[2]=0.8;
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

static int gap(int c)
{
  CControl *I=&Control;
  return((c*I->ExtraSpace)/cControlButtons-((c-1)*I->ExtraSpace)/cControlButtons);
}
/*========================================================================*/
int ControlClick(Block *block,int button,int x,int y,int mod)
{
  CControl *I=&Control;
  int sel = 0;
  int flag = false;
  int c=0;

  I->LastPos =x;
  x -= I->Block->rect.left+cControlLeftMargin;
  y -= I->Block->rect.top-cControlTopMargin;
  if(x<2) {
    OrthoGrab(block);
    I->DragFlag=true;
  }
  if((y<=0)&&(y>(-cControlBoxSize)))
    c=1;
	 while(x>=0) {
		if(x<cControlBoxSize)
		  {
			 flag=true;
			 break;
		  }
		x-=cControlBoxSize+cControlSpacing+gap(c++);
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
  int c=1;

  if(PMGUI) {

    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    glColor3fv(I->Block->TextColor);

    {
      int top,left,bottom,right;

      left = I->Block->rect.left;
      bottom = I->Block->rect.bottom+3;
      top = bottom+cControlBoxSize+1;
      right=left+4;
      
      
      glColor3f(0.8,0.8,0.8);
      glBegin(GL_POLYGON);
      glVertex2i(right,top);
      glVertex2i(right,bottom+1);
      glVertex2i(left,bottom+1);
      glVertex2i(left,top);
      glEnd();
      
      glColor3f(0.3,0.3,0.3);
      glBegin(GL_POLYGON);
      glVertex2i(right,top-1);
      glVertex2i(right,bottom);
      glVertex2i(left+1,bottom);
      glVertex2i(left+1,top-1);
      glEnd();
      
      glColor3f(0.3,0.3,0.3);
      glBegin(GL_POLYGON);
      glVertex2i(right,bottom+1);
      glVertex2i(right,bottom);
      glVertex2i(left,bottom);
      glVertex2i(left,bottom+1);
      glEnd();

      glColor3fv(I->ButtonColor);
      
      glBegin(GL_POLYGON);
      glVertex2i(right-1,top-1);
      glVertex2i(right-1,bottom+1);
      glVertex2i(left+1,bottom+1);
      glVertex2i(left+1,top-1);
      glEnd();
    }


    
    x = I->Block->rect.left+cControlLeftMargin;
    y = I->Block->rect.top-cControlTopMargin;

    glColor3fv(I->ButtonColor);
    glBegin(GL_POLYGON);
    glVertex2i(x,y+1);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y+1);
    glEnd();
    glColor3fv(I->Block->TextColor);

    glBegin(GL_TRIANGLES);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize/2));  
    glEnd();
    glBegin(GL_LINES);
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing+gap(c++);

    glColor3fv(I->ButtonColor);
    glBegin(GL_POLYGON);
    glVertex2i(x,y+1);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y+1);
    glEnd();
    glColor3fv(I->Block->TextColor);

    glBegin(GL_POLYGON);
    glVertex2i(x+cControlBoxSize/2+2,
               y-(cControlBoxSize/2));  
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,
               y-(cControlBoxSize/2));  
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing+gap(c++);


    glColor3fv(I->ButtonColor);
    glBegin(GL_POLYGON);
    glVertex2i(x,y+1);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y+1);
    glEnd();
    glColor3fv(I->Block->TextColor);

    glBegin(GL_POLYGON);
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing+gap(c++);

    if(MoviePlaying()) {
      glColor3fv(I->ActiveColor);
      glBegin(GL_TRIANGLE_STRIP);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize,y+1);
      glVertex2i(x+cControlBoxSize,y-(cControlBoxSize-1));

      glEnd();
      glColor3fv(I->Block->BackColor);
      glBegin(GL_TRIANGLES);
      glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin+1);
      glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin-1);
      glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                 y-(cControlBoxSize/2));  
      glEnd();
      glColor3fv(I->Block->TextColor);
    } else {

      glColor3fv(I->ButtonColor);
      glBegin(GL_POLYGON);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y+1);
      glEnd();
      glColor3fv(I->Block->TextColor);

      glBegin(GL_TRIANGLES);
      glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
      glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
                 y-(cControlBoxSize/2));  
      glEnd();
    }
    x+=cControlBoxSize+cControlSpacing+gap(c++);

    glColor3fv(I->ButtonColor);
    glBegin(GL_POLYGON);
    glVertex2i(x,y+1);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y+1);
    glEnd();
    glColor3fv(I->Block->TextColor);

    glBegin(GL_POLYGON);
    glVertex2i(x+cControlBoxSize/2-2,
               y-(cControlBoxSize/2));  
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize/2));  
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing+gap(c++);
    
    glColor3fv(I->ButtonColor);
    glBegin(GL_POLYGON);
    glVertex2i(x,y+1);
    glVertex2i(x,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
    glVertex2i(x+cControlBoxSize-1,y+1);
    glEnd();
    glColor3fv(I->Block->TextColor);

    glBegin(GL_TRIANGLES);
    glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize/2));  
    glEnd();
    glBegin(GL_LINES);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
    glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
               y-(cControlBoxSize-1)+cControlInnerMargin);
    glEnd();
    x+=cControlBoxSize+cControlSpacing+gap(c++);

    if(SettingGet(cSetting_sculpting)) {

      glColor3fv(I->ActiveColor);
      glBegin(GL_POLYGON);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y+1);
      glEnd();

      glColor3fv(I->Block->BackColor);

      GrapDrawStr("S",x+cControlInnerMargin,
                  y-cControlBoxSize+cControlInnerMargin+1);
      glColor3fv(I->Block->TextColor);    
    } else {
      
      glColor3fv(I->ButtonColor);
      glBegin(GL_POLYGON);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y+1);
      glEnd();
      glColor3fv(I->Block->TextColor);
      GrapDrawStr("S",x+cControlInnerMargin,
                  y-cControlBoxSize+cControlInnerMargin+1);

    }
    x+=cControlBoxSize+cControlSpacing+gap(c++);
    
    if(I->Rocking) {
      glColor3fv(I->ActiveColor);
      glBegin(GL_TRIANGLE_STRIP);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y+1);
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      
      glEnd();
      glColor3fv(I->Block->BackColor);
      glBegin(GL_POLYGON);
      glVertex2i(x+(cControlBoxSize/2)+cControlSpread,
                 y-cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2),
                 y-(cControlBoxSize)+cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2)-cControlSpread,
                 y-cControlInnerMargin);
      glEnd();
      glColor3fv(I->Block->TextColor);
    } else {
      glColor3fv(I->ActiveColor);
      glColor3fv(I->ButtonColor);
      glBegin(GL_POLYGON);
      glVertex2i(x,y+1);
      glVertex2i(x,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y-(cControlBoxSize-1));
      glVertex2i(x+cControlBoxSize-1,y+1);
      glEnd();
      glColor3fv(I->Block->TextColor);
      glBegin(GL_POLYGON);
      glVertex2i(x+(cControlBoxSize/2)+cControlSpread,
                 y-cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2),
                 y-(cControlBoxSize)+cControlInnerMargin);
      glVertex2i(x+(cControlBoxSize/2)-cControlSpread,
                 y-cControlInnerMargin);
      glEnd();
      
    }
    x+=cControlBoxSize+cControlSpacing+gap(c++);



  }
}



