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

#include "main.h"
#include "Base.h"
#include "ButMode.h"
#include "Scene.h"

#define cButModeLineHeight 12
#define cButModeLeftMargin 2
#define cButModeTopMargin 1

CButMode ButMode;

int ButModeRelease(Block *block,int x,int y,int mod);
void ButModeDraw(Block *block);
int ButModeClick(Block *block,int button,int x,int y,int mod);

Block *ButModeGetBlock(void)
{
  CButMode *I=&ButMode;
  {return(I->Block);}
}
/*========================================================================*/
void ButModeSetRate(float interval)
{
  CButMode *I=&ButMode;
  I->Samples*=0.9;
  I->Rate*=0.9;

  I->Samples++;

  if(interval>=0.01)
	 I->Rate += 1/interval;
  else
	 I->Rate += 99;
  
}
/*========================================================================*/
void ButModeResetRate(void)
{
  CButMode *I=&ButMode;
  I->Samples=0.0;
  I->Rate=0.0;
}
/*========================================================================*/
void ButModeFree(void)
{
  CButMode *I=&ButMode;
  OrthoFreeBlock(I->Block);
}
/*========================================================================*/
void ButModeInit(void)
{

  CButMode *I=&ButMode;

  I->Rate=0.0;
  I->Samples = 0.0;
  I->Mode[0]=0;
  I->Mode[1]=1;
  I->Mode[2]=2;
  I->Mode[3]=3;
  I->Mode[4]=2;
  I->Mode[5]=1;
  I->NMode = 4;

  strcpy(I->Code[0],"R-XYZ");
  strcpy(I->Code[1],"T-XY");
  strcpy(I->Code[2],"T-Z");
  strcpy(I->Code[3],"C-ZZ");

  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick = ButModeClick;
  I->Block->fDraw    = ButModeDraw;
  I->Block->fReshape = BlockReshape;
  I->Block->active = true;

  I->Block->TextColor[0]=0.2;
  I->Block->TextColor[1]=1.0;
  I->Block->TextColor[2]=0.2;

  OrthoAttach(I->Block,cOrthoTool);

}


/*========================================================================*/
int ButModeClick(Block *block,int button,int x,int y,int mod)
{
  CButMode *I=&ButMode;
  int n;
  n=((I->Block->rect.top-(y+2))-cButModeTopMargin)/cButModeLineHeight;
  SceneDirty();
  if((n>=0)&&(n<3)) {
	 if((x>(I->Block->rect.left+I->Block->rect.right)/2)) n+=3;
	 I->Mode[n]=(I->Mode[n]+1)%I->NMode;
  }
  return(1);
}

void ButModeDraw(Block *block)
{
  CButMode *I=&ButMode;
  int x,y,a;
  char *c;
  float rate;

  char rateStr[255];

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    glColor3fv(I->Block->TextColor);

    x = I->Block->rect.left+cButModeLeftMargin;
    y = (I->Block->rect.top-cButModeLineHeight)-cButModeTopMargin;
  
    for(a=0;a<3;a++)
      {
        glRasterPos4d((double)(x),(double)(y),0.0,1.0);
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,a+'1');
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,':');
        c=I->Code[I->Mode[a]];
        while(*c) 
          glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
        y-=cButModeLineHeight;
      }

    if(I->Samples) 
      rate = I->Rate/I->Samples;
    else 
      rate = 0;
    sprintf(rateStr,"R: %1.1f Frame: %2d",rate,SceneGetFrame()+1);
  
    glRasterPos4d((double)(x),(double)y,0.0,1.0);
    c=rateStr;
    while(*c)
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));

    x = (I->Block->rect.left+I->Block->rect.right)/2;
    y = (I->Block->rect.top-cButModeLineHeight)-cButModeTopMargin;
    for(a=3;a<6;a++)
      {
        glRasterPos4d((double)(x),(double)(y),0.0,1.0);
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,a+'1');
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,':');
        c=I->Code[I->Mode[a]];
        while(*c) 
          glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
        y-=cButModeLineHeight;
      }

  }
}


