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
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <GL/glut.h>

#include "Base.h"
#include "Control.h"
#include "Scene.h"
#include "Executive.h"
#include"Movie.h"

#define cControlBoxSize 25
#define cControlLeftMargin 4
#define cControlTopMargin 5
#define cControlSpacing 5
#define cControlInnerMargin 4
#define cControlSpread 6

CControl Control;

int ControlRelease(Block *block,int x,int y);
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
  return(MoviePlaying()||I->Rocking);
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
		MoviePlay(cMovieStop);
		ExecutiveDrawNow();
		break;
	 case 1:
		if(!(mod&cOrthoCTRL)) {
		  SceneSetFrame(0,0);		
		  if(mod&cOrthoSHIFT) 
			 MoviePlay(cMoviePlayRendered);
		  else 
			 MoviePlay(cMoviePlay);
		} else {
		  MoviePlay(cMoviePlay);
		}
		break;
	 case 2:
		if(mod&cOrthoSHIFT) {
		  SceneSetFrame(1,-1);
		} else {
		  SceneSetFrame(0,0);
		}
		break;
	 case 3:
		if(mod&cOrthoSHIFT) {
		  SceneSetFrame(1,1);
		} else if(mod&cOrthoCTRL) {
		  SceneSetFrame(3,0);
		} else {
		  SceneSetFrame(2,0);
		}
		break;
	 case 4:
		I->Rocking=!I->Rocking;
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
  glVertex2i(x+cControlInnerMargin,y-cControlInnerMargin);
  glVertex2i(x+cControlInnerMargin,y-(cControlBoxSize-1)+cControlInnerMargin);
  glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,
				 y-(cControlBoxSize-1)+cControlInnerMargin);
  glVertex2i(x+(cControlBoxSize-1)-cControlInnerMargin,y-cControlInnerMargin);
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
  glEnd();
  x+=cControlBoxSize+cControlSpacing;

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
  x+=cControlBoxSize+cControlSpacing;  
}


