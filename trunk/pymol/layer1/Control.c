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
#define cControlTopMargin 2
#define cControlSpacing 2
#define cControlInnerMargin 4
#define cControlSpread 6
#define cControlSize 160

#define cControlButtons 7

CControl Control;

int ControlClick(Block *block,int button,int x,int y,int mod);

static void ControlReshape(Block *block,int width, int height)
{
  CControl *I=&Control;
  BlockReshape(block,width,height);

  I->ExtraSpace = ((block->rect.right-block->rect.left)-cControlSize);
  if(I->ExtraSpace<0)
    I->ExtraSpace=0;
  
}

static int which_button(int x,int y)
{
  int result = -1;
  CControl *I=&Control;
  x -= I->Block->rect.left+cControlLeftMargin;
  y -= I->Block->rect.top-cControlTopMargin;
  if(x>=0) 
    if((y<=0)&&(y>(-cControlBoxSize))) {
      int control_width = I->Block->rect.right - (I->Block->rect.left+cControlLeftMargin);
      int nButton = 8;
      result = (nButton*x)/control_width;
    }
  return result;
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
      gui_width = (int)SettingGet(cSetting_internal_gui_width)-delta;
      if(gui_width<3)
        gui_width = 3;
      delta = (int)SettingGet(cSetting_internal_gui_width)-gui_width;
      width = OrthoGetWidth()+delta;
    I->LastPos = x;
    SettingSet(cSetting_internal_gui_width,(float)gui_width);
    OrthoReshape(-1,-1);
    }
  } else {
    I->Active = which_button(x,y);
    if(I->Active!=I->Pressed)
      I->Active = -1;
    OrthoDirty();
  }
  return(1);
}

static int ControlRelease(Block *block,int button,int x,int y,int mod)
{
  CControl *I=&Control;  

  int sel = 0;

  I->LastPos =x;
  sel = which_button(x,y);

  switch(sel) {
  case 0:
    SceneSetFrame(4,0);
    PLog("cmd.rewind()",cPLog_pym);
    break;
  case 1:
    SceneSetFrame(5,-1);
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
        SceneSetFrame(4,0);		
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
    SceneSetFrame(5,1);
    PLog("cmd.forward()",cPLog_pym);
    break;
  case 5:
    if(mod&cOrthoCTRL) {
      SceneSetFrame(3,0);
      PLog("cmd.middle()",cPLog_pym);
    } else {
      SceneSetFrame(6,0);
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
  OrthoDirty();
  OrthoUngrab();
  I->DragFlag=false;
  I->Active = -1;
  I->Pressed = -1;

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
  
  if(x<(I->Block->rect.left + cControlLeftMargin)) {
    y -= I->Block->rect.top-cControlTopMargin;
    if((y<=0)&&(y>(-cControlBoxSize))) {
      I->LastPos = x;
      OrthoGrab(block);
      I->DragFlag=true;
    }
  } else {
    I->Pressed = which_button(x,y);
    I->Active = I->Pressed;
    if(I->Pressed)
      OrthoGrab(block);
    OrthoDirty();
  }
  return(1);
}

static void draw_button(int x2,int y2, int w, int h, float *light, float *dark, float *inside)
{
  glColor3fv(light);
  glBegin(GL_POLYGON);
  glVertex2i(x2,y2);
  glVertex2i(x2,y2+h);
  glVertex2i(x2+w,y2+h);
  glVertex2i(x2+w,y2);
  glEnd();
  
  glColor3fv(dark);
  glBegin(GL_POLYGON);
  glVertex2i(x2+1,y2);
  glVertex2i(x2+1,y2+h-1);
  glVertex2i(x2+w,y2+h-1);
  glVertex2i(x2+w,y2);
  glEnd();
  
  glColor3fv(inside);
  glBegin(GL_POLYGON);
  glVertex2i(x2+1,y2+1);
  glVertex2i(x2+1,y2+h-1);
  glVertex2i(x2+w-1,y2+h-1);
  glVertex2i(x2+w-1,y2+1);
  glEnd();

}
/*========================================================================*/
static void ControlDraw(Block *block)
{
  CControl *I=&Control;
  int x,y;
  int nButton = 8;
  int but_num;
  float lightEdge[3] = { 0.7F, 0.7F, 0.7F};
  float darkEdge[3] = {0.3F, 0.3F, 0.3F};
  float active[3] = {0.8F,0.8F,0.8F};

  if(PMGUI) {
    int control_width = I->Block->rect.right - (I->Block->rect.left+cControlLeftMargin);

    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    glColor3fv(I->Block->TextColor);

    {
      int top,left,bottom,right;
      
      left = I->Block->rect.left+1;
      bottom = I->Block->rect.bottom+1;
      top = I->Block->rect.top - (cControlTopMargin-1);
      right=left+5;
      
      glColor3f(0.8F,0.8F,0.8F);
      glBegin(GL_POLYGON);
      glVertex2i(right,top);
      glVertex2i(right,bottom);
      glVertex2i(left,bottom);
      glVertex2i(left,top);
      glEnd();
      
      glColor3f(0.3F,0.3F,0.3F);
      glBegin(GL_POLYGON);
      glVertex2i(right,top-1);
      glVertex2i(right,bottom);
      glVertex2i(left+1,bottom);
      glVertex2i(left+1,top-1);
      glEnd();

      glColor3fv(I->ButtonColor);
      
      glBegin(GL_POLYGON);
      glVertex2i(right-1,top-1);
      glVertex2i(right-1,bottom+1);
      glVertex2i(left+1,bottom+1);
      glVertex2i(left+1,top-1);
      glEnd();
    }

    y = I->Block->rect.top-cControlTopMargin;
    
    for(but_num=0;but_num<nButton;but_num++) {
      int but_width;
      int but_left;
      int but_bottom;
      int but_height;

      but_left = I->Block->rect.left + cControlLeftMargin + (but_num*control_width)/nButton;
      but_width = (((but_num+1)*control_width/nButton) - 
                   ((but_num)*control_width/nButton)) - 1;
      
      but_bottom = y-(cControlBoxSize-1);
      but_height = cControlBoxSize;
      

      if( ( but_num==I->Active ) ) {
        draw_button(but_left,but_bottom,
                    but_width, but_height, lightEdge,darkEdge,active);
      } else if(((but_num==6)&&((int)SettingGet(cSetting_sculpting))) ||
                ((but_num==3)&&(MoviePlaying())) ||
                ((but_num==7)&&(I->Rocking))) {
        draw_button(but_left,but_bottom,
                    but_width, but_height, lightEdge,darkEdge,I->ActiveColor);
      } else {
        draw_button(but_left,but_bottom,
                    but_width, but_height, lightEdge,darkEdge,I->ButtonColor);
      }

      x = but_left + (but_width-cControlBoxSize)/2;

      glColor3fv(I->Block->TextColor);
      switch(but_num) {
      case 0:
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
        break;
        
      case 1:

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
        break;
      case 2:
        glBegin(GL_POLYGON);
        glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
        glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
        glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
                   y-(cControlBoxSize-1)+cControlInnerMargin);
        glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
        glEnd();
        break;

      case 3:
        glBegin(GL_TRIANGLES);
        glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin+1);
        glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin-1);
        glVertex2i(x+(cControlBoxSize)-cControlInnerMargin,
                   y-(cControlBoxSize/2));  
        glEnd();
        break;
      case 4:
        glBegin(GL_POLYGON);
        glVertex2i(x+cControlBoxSize/2-2,
                   y-(cControlBoxSize/2));  
        glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
        glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
                   y-(cControlBoxSize/2));  
        glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
        glEnd();

        break;
      case 5:
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
        break;
      case 6:
        GrapDrawStr("S",x+cControlInnerMargin,
                    y-cControlBoxSize+cControlInnerMargin+1);
        break;
      case 7:
        glBegin(GL_POLYGON);
        glVertex2i(x+(cControlBoxSize/2)+cControlSpread,
                   y-cControlInnerMargin);
        glVertex2i(x+(cControlBoxSize/2),
                   y-(cControlBoxSize)+cControlInnerMargin);
        glVertex2i(x+(cControlBoxSize/2)-cControlSpread,
                   y-cControlInnerMargin);
        glEnd();
        break;
      }
    }

  }
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
  I->ButtonColor[0]=0.5F;
  I->ButtonColor[1]=0.5F;
  I->ButtonColor[2]=0.5F;
  I->ActiveColor[0]=0.7F;
  I->ActiveColor[1]=0.7F;
  I->ActiveColor[2]=0.7F;
  I->Pressed = -1;
  I->Active = -1;
  OrthoAttach(I->Block,cOrthoTool);

  I->Rocking=false;

}
