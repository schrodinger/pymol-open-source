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
#include "ButMode.h"
#include "Scene.h"
#include "Util.h"
#include "Grap.h"
#include "Ortho.h"

#define cButModeLineHeight 12
#define cButModeLeftMargin 2
#define cButModeTopMargin 1

CButMode ButMode;

int ButModeRelease(Block *block,int x,int y,int mod);
void ButModeDraw(Block *block);
int ButModeClick(Block *block,int button,int x,int y,int mod);

/*========================================================================*/
Block *ButModeGetBlock(void)
{
  CButMode *I=&ButMode;
  {return(I->Block);}
}
/*========================================================================*/
void ButModeSet(int button,int action)
{
  CButMode *I=&ButMode;
  if((button>=0)&&(button<I->NBut)&&
     (action>=0)&&(action<I->NCode)) {
    I->Mode[button]=action;
    OrthoDirty();
  }
}
/*========================================================================*/
void ButModeCaption(char *text)
{
  CButMode *I=&ButMode;
  int l;
  l = strlen(I->Caption);
  if((l>0)&&(l<(sizeof(WordType)-1)))
    strcat(I->Caption,",");
  l = (sizeof(WordType)-2)-l;
  UtilNConcat(I->Caption,text,l);
}
/*========================================================================*/
void ButModeCaptionReset(void)
{
  CButMode *I=&ButMode;
  I->Caption[0]=0;
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
  int a;

  I->Rate=0.0;
  I->Samples = 0.0;

  /* Outdated code below:  Mouse configuration now handled completely from Python
     I->Mode[0]=cButModeRotXYZ;
     I->Mode[1]=cButModeTransXY;
     I->Mode[2]=cButModeTransZ;
     I->Mode[3]=cButModeRotZ;
     I->Mode[4]=cButModeTransXY;
     I->Mode[5]=cButModeClipNF;
     I->Mode[6]=cButModePk1;
     I->Mode[7]=cButModePk2;
     I->Mode[8]=cButModePk3;
     I->Mode[9]=cButModeAddToPk1;
     I->Mode[10]=cButModeOrigAt;
     I->Mode[11]=cButModeAddToPk3;
  */

  I->Caption[0] = 0;

  I->NCode = 19;
  I->NBut = 12;

  for(a=0;a<I->NBut;a++) {
    I->Mode[a]=-1;
  }
    
  strcpy(I->Code[cButModeRotXYZ],  "Rota ");
  strcpy(I->Code[cButModeRotZ],    "RotZ ");  
  strcpy(I->Code[cButModeTransXY], "Move ");
  strcpy(I->Code[cButModeTransZ],  "MovZ ");
  strcpy(I->Code[cButModeClipNF],  "Clip ");
  strcpy(I->Code[cButModeClipN],   "ClpN ");  
  strcpy(I->Code[cButModeClipF],   "ClpF ");
  strcpy(I->Code[cButModePickAtom],"PkAt ");
  strcpy(I->Code[cButModePickBond],"PkBd ");
  strcpy(I->Code[cButModeTorFrag], "TorF ");
  strcpy(I->Code[cButModeRotFrag], "RotF ");
  strcpy(I->Code[cButModeMovFrag], "MovF ");
  strcpy(I->Code[cButModePk1],     " Pk1 ");
  strcpy(I->Code[cButModePk2],     " Pk2 ");
  strcpy(I->Code[cButModePk3],     " Pk3 ");
  strcpy(I->Code[cButModeAddToPk1],"+Pk1 ");
  strcpy(I->Code[cButModeAddToPk2],"+Pk2 ");
  strcpy(I->Code[cButModeAddToPk3],"+Pk3 ");
  strcpy(I->Code[cButModeOrigAt],  "Orig ");

  I->Block = OrthoNewBlock(NULL);
  I->Block->fClick = ButModeClick;
  I->Block->fDraw    = ButModeDraw;
  I->Block->fReshape = BlockReshape;
  I->Block->active = true;

  I->Block->TextColor[0]=0.2;
  I->Block->TextColor[1]=1.0;
  I->Block->TextColor[2]=0.2;

  I->TextColor1[0]=0.5;
  I->TextColor1[1]=0.5;
  I->TextColor1[2]=1.0;

  I->TextColor2[0]=0.8;
  I->TextColor2[1]=0.8;
  I->TextColor2[2]=0.8;

  OrthoAttach(I->Block,cOrthoTool);

}
/*========================================================================*/
int ButModeTranslate(int button, int mod)
{
  int mode = 0;
  CButMode *I=&ButMode;
  switch(button) {
  case GLUT_LEFT_BUTTON:
    mode = 0;
    break;
  case GLUT_MIDDLE_BUTTON:
    mode = 1;
    break;
  case GLUT_RIGHT_BUTTON:
    mode = 2;
    break;
  }
  switch(mod) {
  case 0:
    break;
  case cOrthoSHIFT:
    mode+=3;
    break;
  case cOrthoCTRL:
    mode+=6;
    break;
  case (cOrthoCTRL+cOrthoSHIFT):
    mode+=9;
    break;
  }
  return(I->Mode[mode]);
}
/*========================================================================*/
int ButModeClick(Block *block,int button,int x,int y,int mod)
{
  OrthoCommandIn("edit");
  return(1);
}
/*========================================================================*/
void ButModeDraw(Block *block)
{
  CButMode *I=&ButMode;
  int x,y,a;
  float rate;
  char rateStr[255];
  int mode;

  if(PMGUI) {
    glColor3fv(I->Block->BackColor);
    BlockFill(I->Block);
    glColor3fv(I->Block->TextColor);

    x = I->Block->rect.left+cButModeLeftMargin;
    y = (I->Block->rect.top-cButModeLineHeight)-cButModeTopMargin;

    GrapDrawStr("Mouse:",x,y);
    glColor3fv(I->TextColor1);
    GrapDrawStr("  L    M    R",x+40,y);

    y-=cButModeLineHeight;
    GrapDrawStr("None ",x,y);
    glColor3fv(I->TextColor2);
    glRasterPos4d(x+40,y,0,1);
    for(a=0;a<3;a++) {
      mode = I->Mode[a];
      if(mode<0)
        GrapContStr("    ");
      else
        GrapContStr(I->Code[mode]);
    }

    y-=cButModeLineHeight;
    glColor3fv(I->TextColor1);
    GrapDrawStr("Shft ",x,y);
    glColor3fv(I->TextColor2);
    glRasterPos4d(x+40,y,0,1);
    for(a=3;a<6;a++) {
      mode = I->Mode[a];
      if(mode<0)
        GrapContStr("    ");
      else
        GrapContStr(I->Code[mode]);
    }

    y-=cButModeLineHeight;
    glColor3fv(I->TextColor1);
    GrapDrawStr("Ctrl ",x,y);
    glColor3fv(I->TextColor2);
    glRasterPos4d(x+40,y,0,1);
    for(a=6;a<9;a++) {
      mode = I->Mode[a];
      if(mode<0)
        GrapContStr("    ");
      else
        GrapContStr(I->Code[mode]);
    }

    y-=cButModeLineHeight;
    glColor3fv(I->TextColor1);
    GrapDrawStr("CtSh ",x,y);
    glColor3fv(I->TextColor2);
    glRasterPos4d(x+40,y,0,1);
    for(a=9;a<12;a++) {
      mode = I->Mode[a];
      if(mode<0)
        GrapContStr("    ");
      else
        GrapContStr(I->Code[mode]);
    }
    glColor3fv(I->Block->TextColor);
    y-=cButModeLineHeight;

    if(I->Caption[0]) GrapDrawStr(I->Caption,x,y);

    y-=cButModeLineHeight;
    if(I->Samples) 
      rate = I->Rate/I->Samples;
    else 
      rate = 0;
    sprintf(rateStr,"Frame:%3d %4.1f FPS",SceneGetFrame()+1,rate);
    GrapDrawStr(rateStr,x,y);

  }
}


