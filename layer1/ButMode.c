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
void ButModeChange(int mode)
{
  /* THIS STUFF SHOULD DEFINITELY BE MOVED INTO PYTHON */
  CButMode *I=&ButMode;
  switch(mode) {
  case 0: 
    I->Mode[0]=cButModeRotXYZ;
    I->Mode[1]=cButModeTransXY;
    I->Mode[2]=cButModeTransZ;
    I->Mode[3]=cButModeRotZ;
    I->Mode[4]=cButModeTransXY;
    I->Mode[5]=cButModeClipZZ;
    break;
  case 1: 
    I->Mode[0]=cButModeRotXYZ;
    I->Mode[1]=cButModeTransXY;
    I->Mode[2]=cButModeTransZ;
    I->Mode[3]=cButModeTransXY;
    I->Mode[4]=cButModeRotZ;
    I->Mode[5]=cButModeClipZZ;
    break;
  case 2: 
    I->Mode[0]=cButModeRotXYZ;
    I->Mode[1]=cButModeTransXY;
    I->Mode[2]=cButModeTransZ;
    I->Mode[3]=cButModeRotZ;
    I->Mode[4]=cButModeClipN;
    I->Mode[5]=cButModeClipF;
    break;
  case 3: 
    I->Mode[2]=cButModeRotXYZ;
    I->Mode[1]=cButModeTransXY;
    I->Mode[0]=cButModeTransZ;
    I->Mode[5]=cButModeRotZ;
    I->Mode[4]=cButModeTransXY;
    I->Mode[3]=cButModeClipZZ;
    break;
  case 4: 
    I->Mode[2]=cButModeRotXYZ;
    I->Mode[1]=cButModeTransXY;
    I->Mode[0]=cButModeTransZ;
    I->Mode[5]=cButModeTransXY;
    I->Mode[4]=cButModeRotZ;
    I->Mode[3]=cButModeClipZZ;
    break;
  case 5: 
    I->Mode[2]=cButModeRotXYZ;
    I->Mode[1]=cButModeTransXY;
    I->Mode[0]=cButModeTransZ;
    I->Mode[5]=cButModeRotZ;
    I->Mode[4]=cButModeClipN;
    I->Mode[3]=cButModeClipF;
    break;

  }
  OrthoDirty();
}
/*========================================================================*/
void ButModeSetRate(float interval)
{
  CButMode *I=&ButMode;
  I->Samples*=0.9;
  I->Rate*=0.9;

  I->Samples++;

  if(interval>=0.001)
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
  
  I->Mode[0]=cButModeRotXYZ;
  I->Mode[1]=cButModeTransXY;
  I->Mode[2]=cButModeTransZ;
  I->Mode[3]=cButModeRotZ;
  I->Mode[4]=cButModeTransXY;
  I->Mode[5]=cButModeClipZZ;

  I->NCode = 7;

  strcpy(I->Code[cButModeRotXYZ],"R-XYZ");
  strcpy(I->Code[cButModeRotZ],"R-Z");  
  strcpy(I->Code[cButModeTransXY],"T-XY");
  strcpy(I->Code[cButModeTransZ],"T-Z");
  strcpy(I->Code[cButModeClipZZ],"C-NF");
  strcpy(I->Code[cButModeClipN],"C-N");  
  strcpy(I->Code[cButModeClipF],"C-F");

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
	 I->Mode[n]=(I->Mode[n]+1)%I->NCode;
  }
  return(1);
}

void ButModeDraw(Block *block)
{
  CButMode *I=&ButMode;
  int x,y,a;
  char *c;
  float rate;
  char but[3] = { 'L', 'M', 'R' };

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
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,but[a]);
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
    sprintf(rateStr,"Frame:%3d %4.1f FPS",SceneGetFrame()+1,rate);
  
    glRasterPos4d((double)(x),(double)y,0.0,1.0);
    c=rateStr;
    while(*c)
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));

    x = (I->Block->rect.left+80+cButModeLeftMargin);
    y = (I->Block->rect.top-cButModeLineHeight)-cButModeTopMargin;
    for(a=3;a<6;a++)
      {
        glRasterPos4d((double)(x),(double)(y),0.0,1.0);
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,'S');
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,but[a-3]);
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13,':');
        c=I->Code[I->Mode[a]];
        while(*c) 
          glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*(c++));
        y-=cButModeLineHeight;
      }

  }
}


