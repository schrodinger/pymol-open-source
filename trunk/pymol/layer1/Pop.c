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
#include"os_predef.h"
#include"Base.h"
#include"Ortho.h"
#include "Pop.h"

#define cPopMargin 4

typedef struct {
  Block *Block;
} CPop;

CPop Pop;

void PopReshape(Block *I,int width, int height);

/*========================================================================*/
void PopReshape(Block *I,int width, int height)
{
  I->rect.top=height;
  I->rect.right=width;
}
/*========================================================================*/
Block *PopGetBlock(void)
{
  CPop *I=&Pop;
  {return(I->Block);}
}
/*========================================================================*/
void PopFree(void)
{
  CPop *I=&Pop;
  OrthoFreeBlock(I->Block);
}
/*========================================================================*/
void PopInit(void)
{

  CPop *I=&Pop;

  I->Block = OrthoNewBlock(NULL);
  I->Block->fReshape = PopReshape;
  I->Block->active = false;

  I->Block->rect.top=10;
  I->Block->rect.bottom=0;
  I->Block->rect.left=0;
  I->Block->rect.right=10;

  OrthoAttach(I->Block,cOrthoHidden);

}
/*========================================================================*/
void PopFitBlock(Block *block)
{
  CPop *I=&Pop;
  int delta;

  if((block->rect.bottom - cPopMargin)< I->Block->rect.bottom) {
    delta=(I->Block->rect.bottom - block->rect.bottom) + cPopMargin;
    block->rect.top+=delta;
    block->rect.bottom+=delta;
  }

  if((block->rect.right+cPopMargin) > I->Block->rect.right) {
    delta=(block->rect.right - I->Block->rect.right)  + cPopMargin;
    block->rect.left-=delta;
    block->rect.right-=delta;
  }

  if((block->rect.left-cPopMargin) < I->Block->rect.left) {
    delta=(I->Block->rect.left - block->rect.left) + cPopMargin;
    block->rect.right+=delta;
    block->rect.left+=delta;
  }

  if((block->rect.top+cPopMargin) > I->Block->rect.top) {
    delta=(block->rect.top - I->Block->rect.top)  + cPopMargin;
    block->rect.top-=delta;
    block->rect.bottom-=delta;
  }

}

