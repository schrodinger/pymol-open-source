/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2002 by Warren Lyford Delano of DeLano Scientific. 
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
#include"os_gl.h"
#include"Base.h"
#include"Block.h"
#include"OOMac.h"
#include"Ortho.h"
#include"ScrollBar.h"

typedef struct CScrollBar {
  Block *Block;
  int HorV;
  float BackColor[3];
  float BarColor[3];
  int ListSize;
  int DisplaySize;
  int BarSize;
  float ExactBarSize;
  float Value;
  float StartValue;
  float ValueMax;
  int BarRange;
  int BarMin;
  int BarMax;
  int StartPos;
} CScrollBar;

void ScrollBarMaxOut(struct CScrollBar *I)
{
  I->Value = I->ValueMax;
}

int ScrollBarIsMaxed(struct CScrollBar *I)
{
  if(I->ValueMax>0.0F) {
    if(I->Value==I->ValueMax)
      return true;
    else
      return false;
  } else 
    return false;
}

void ScrollBarUpdate(struct CScrollBar *I)
{
  int range;

  if(I->HorV) {
    range = (I->Block->rect.right-I->Block->rect.left);
  } else {
    range = (I->Block->rect.top-I->Block->rect.bottom);
  }
  I->ExactBarSize = (range*I->DisplaySize)/(float)I->ListSize;
  I->BarSize = (int)I->ExactBarSize;
  if(I->BarSize<4)
    I->BarSize=4;
  I->BarRange = range - I->BarSize;
  if(I->BarRange<2)
    I->BarRange=2;
  I->ValueMax = (float)I->ListSize-I->DisplaySize;
  if(I->ValueMax<1)
    I->ValueMax=1;
  if(I->Value>I->ValueMax)
    I->Value=(float)I->ValueMax;

}

static  void ScrollBarDraw(Block *block)
{
  PyMOLGlobals *G=block->G;
  float value;
  int top,left,bottom,right;

  CScrollBar *I = (CScrollBar*)block->reference;
  glColor3fv(I->BackColor);
  BlockFill(I->Block);

  ScrollBarUpdate(I);

  value = I->Value;
  if(value>I->ValueMax)
    value=I->ValueMax;

  if(I->HorV) {
    top = block->rect.top-1;
    bottom = block->rect.bottom+1;
    left = (int)(block->rect.left+(I->BarRange*value)/I->ValueMax);
    right = left+I->BarSize;
    I->BarMin = left;
    I->BarMax = right;
  } else {
    top = (int)(block->rect.top-(I->BarRange*value)/I->ValueMax);
    bottom = top-I->BarSize;
    left = block->rect.left+1;
    right = block->rect.right-1;
    I->BarMin = top;
    I->BarMax = bottom;
  }

  if(G->HaveGUI && G->ValidContext) {

    glColor3f(0.8F,0.8F,0.8F);
    glBegin(GL_POLYGON);
    glVertex2i(right,top);
    glVertex2i(right,bottom+1);
    glVertex2i(left,bottom+1);
    glVertex2i(left,top);
    glEnd();
    
    glColor3f(0.3F,0.3F,0.3F);
    glBegin(GL_POLYGON);
    glVertex2i(right,top-1);
    glVertex2i(right,bottom);
    glVertex2i(left+1,bottom);
    glVertex2i(left+1,top-1);
    glEnd();
    
    glColor3f(0.3F,0.3F,0.3F);
    glBegin(GL_POLYGON);
    glVertex2i(right,bottom+1);
    glVertex2i(right,bottom);
    glVertex2i(left,bottom);
    glVertex2i(left,bottom+1);
    glEnd();
    
    glColor3fv(I->BarColor);
    glBegin(GL_POLYGON);
    glVertex2i(right-1,top-1);
    glVertex2i(right-1,bottom+1);
    glVertex2i(left+1,bottom+1);
    glVertex2i(left+1,top-1);
    glEnd();
  }
}

void ScrollBarDrawHandle(struct CScrollBar *I,float alpha)
{
  float value;
  int top,left,bottom,right;
  Block *block = I->Block;
  PyMOLGlobals *G=block->G;

  value = I->Value;
  if(value>I->ValueMax)
    value=I->ValueMax;

  if(I->HorV) {
    top = block->rect.top-1;
    bottom = block->rect.bottom+1;
    left = (int)(block->rect.left+(I->BarRange*value)/I->ValueMax);
    right = left+I->BarSize;
  } else {
    top = (int)(block->rect.top-(I->BarRange*value)/I->ValueMax);
    bottom = top-I->BarSize;
    left = block->rect.left+1;
    right = block->rect.right-1;
  }

  if(G->HaveGUI && G->ValidContext ) {

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    
    glColor4f(0.8F,0.8F,0.8F,alpha);
    glBegin(GL_POLYGON);
    glVertex2i(right,top);
    glVertex2i(right,bottom+1);
    glVertex2i(left,bottom+1);
    glVertex2i(left,top);
    glEnd();
    
    glColor4f(0.3F,0.3F,0.3F,alpha);
    glBegin(GL_POLYGON);
    glVertex2i(right,top-1);
    glVertex2i(right,bottom);
    glVertex2i(left+1,bottom);
    glVertex2i(left+1,top-1);
    glEnd();
    
    glColor4f(0.3F,0.3F,0.3F,alpha);
    glBegin(GL_POLYGON);
    glVertex2i(right,bottom+1);
    glVertex2i(right,bottom);
    glVertex2i(left,bottom);
    glVertex2i(left,bottom+1);
    glEnd();
    
    glColor4f(I->BarColor[0],I->BarColor[1],I->BarColor[2],alpha);
    glBegin(GL_POLYGON);
    glVertex2i(right-1,top-1);
    glVertex2i(right-1,bottom+1);
    glVertex2i(left+1,bottom+1);
    glVertex2i(left+1,top-1);
    glEnd();

    glDisable(GL_BLEND);
  }
}


void ScrollBarSetValue(struct CScrollBar *I,float value)
{
  I->Value=value;
  ScrollBarUpdate(I);
}
float ScrollBarGetValue(struct CScrollBar *I)
{
  return(I->Value);
}
static void ScrollBarReshape(Block *block,int width,int height)
{
}

static int ScrollBarClick(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  CScrollBar *I = (CScrollBar*)block->reference;

  if(I->HorV) {
    if(x>I->BarMax) {
      switch(button) {
      case P_GLUT_MIDDLE_BUTTON:
        {
          I->Value= (I->ListSize*(x-block->rect.left))/(block->rect.right - block->rect.left) - I->DisplaySize*0.5F;
          if(I->Value > I->ValueMax)
            I->Value = I->ValueMax;
          OrthoGrab(G,I->Block);
          I->StartPos = x;
          I->StartValue=I->Value;
        }
        break;
      default:
        I->Value+=I->DisplaySize;
        if(I->Value > I->ValueMax)
          I->Value = I->ValueMax;
      }

      OrthoDirty(G);
    } else if(x<I->BarMin){
      switch(button) {
      case P_GLUT_MIDDLE_BUTTON:
        {
          I->Value= (I->ListSize*(x-block->rect.left))/(block->rect.right - block->rect.left) - I->DisplaySize*0.5F;
          if(I->Value<0.0)
            I->Value=0.0F;
          OrthoGrab(G,I->Block);
          I->StartPos = x;
          I->StartValue=I->Value;
        }
        break;
      default:
        I->Value-=I->DisplaySize;
        if(I->Value<0.0)
          I->Value=0.0F;
      }
      OrthoDirty(G);
    } else {
      OrthoGrab(G,I->Block);
      I->StartPos = x;
      I->StartValue=I->Value;
      OrthoDirty(G);
    } 
  } else {
    if(y>I->BarMin) {
      switch(button) {
      case P_GLUT_MIDDLE_BUTTON:
        {
          I->Value= (I->ListSize*(y-block->rect.top))/(block->rect.bottom - block->rect.top) - I->DisplaySize*0.5F;
          if(I->Value<0.0)
            I->Value=0.0F;
          OrthoGrab(G,I->Block);
          I->StartPos = y;
          I->StartValue=I->Value;
        }
        break;
      default:
        I->Value-=I->DisplaySize;
        if(I->Value<0.0)
          I->Value=0.0F;
      }
      OrthoDirty(G);
    } else if(y<I->BarMax) {
      switch(button) {
      case P_GLUT_MIDDLE_BUTTON:
        {
          I->Value= (I->ListSize*(y-block->rect.top))/(block->rect.bottom - block->rect.top) - I->DisplaySize*0.5F;
          if(I->Value > I->ValueMax)
            I->Value = I->ValueMax;
          OrthoGrab(G,I->Block);
          I->StartPos = y;
          I->StartValue=I->Value;
        }
        break;
      default:
        I->Value+=I->DisplaySize;
        if(I->Value > I->ValueMax)
          I->Value = I->ValueMax;
      }
      OrthoDirty(G);
    } else {
      OrthoGrab(G,I->Block);
      I->StartPos = y;
      I->StartValue=I->Value;
      OrthoDirty(G);
    }
  }
  return 0;
}

static int ScrollBarDrag(Block *block,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  CScrollBar *I = (CScrollBar*)block->reference;
  int displ;
  if(I->HorV) 
    displ = I->StartPos-x;
  else 
    displ = y-I->StartPos;    
  I->Value = I->StartValue - (I->ValueMax*displ)/I->BarRange;
  /*  if(displ>0.0)
    I->Value-=0.5;
  else
    I->Value+=0.5;
  */

  if(I->Value<0.0) I->Value=0.0;
  if(I->Value>I->ValueMax) I->Value=I->ValueMax;
  OrthoDirty(G);
  return 0;
}
static int ScrollBarRelease(Block *block,int button,int x,int y,int mod)
{
  PyMOLGlobals *G=block->G;
  OrthoUngrab(G);
  return 0;
}

Block *ScrollBarGetBlock(struct CScrollBar *I)
{
  return(I->Block);
}

void ScrollBarSetLimits(struct CScrollBar *I,int list_size,int display_size)
{
  I->ListSize = list_size;
  I->DisplaySize = display_size;
  ScrollBarUpdate(I);
}


void ScrollBarSetBox(struct CScrollBar *I,int top,int left,int bottom, int right)
{
  I->Block->rect.top = top;
  I->Block->rect.left = left;
  I->Block->rect.bottom = bottom;
  I->Block->rect.right = right;
}
void ScrollBarDoDraw(struct CScrollBar *I)
{
  if(I->Block->fDraw)
    I->Block->fDraw(I->Block);
}

void ScrollBarDoRelease(struct CScrollBar *I,int button,int x,int y,int mod)
{
  if(I->Block->fRelease)
    I->Block->fRelease(I->Block,button,x,y,mod);
}

void ScrollBarDoDrag(struct CScrollBar *I,int x,int y,int mod)
{
  if(I->Block->fDrag)
    I->Block->fDrag(I->Block,x,y,mod);
}

void ScrollBarDoClick(struct CScrollBar *I,int button,int x,int y,int mod)
{
  if(I->Block->fClick)
    I->Block->fClick(I->Block,button,x,y,mod);
}

struct CScrollBar *ScrollBarNew(PyMOLGlobals *G,int horizontal)
{
  OOAlloc(G,CScrollBar)

  I->Block = OrthoNewBlock(G,NULL);  
  I->Block->fRelease = ScrollBarRelease;
  I->Block->fClick   = ScrollBarClick;
  I->Block->fDrag    = ScrollBarDrag;
  I->Block->fDraw    = ScrollBarDraw;
  I->Block->fReshape = ScrollBarReshape;
  I->Block->active = false;
  I->Block->reference = (void*)I;
  I->HorV = horizontal;
  I->BackColor[0]=0.1F;
  I->BackColor[1]=0.1F;
  I->BackColor[2]=0.1F;
  I->BarColor[0]=0.5F;
  I->BarColor[1]=0.5F;
  I->BarColor[2]=0.5F;
  I->ListSize = 10;
  I->DisplaySize = 7;
  I->Value = 0.0F;
  I->ValueMax = 0.0F;
  return(I);
}

void ScrollBarFree(struct CScrollBar *I)
{
  PyMOLGlobals *G=I->Block->G;
  OrthoFreeBlock(G,I->Block);
  OOFreeP(I);
}
