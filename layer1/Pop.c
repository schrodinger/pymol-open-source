
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
#include"Pop.h"
#include"MemoryDebug.h"

#define cPopMargin 3

struct _CPop {
  Block *Block;
};

void PopReshape(Block * I, int width, int height);


/*========================================================================*/
void PopReshape(Block * I, int width, int height)
{
  I->rect.top = height;
  I->rect.right = width;
}


/*========================================================================*/
Block *PopGetBlock(PyMOLGlobals * G)
{
  register CPop *I = G->Pop;
  {
    return (I->Block);
  }
}


/*========================================================================*/
void PopFree(PyMOLGlobals * G)
{
  register CPop *I = G->Pop;
  OrthoFreeBlock(G, I->Block);
  FreeP(G->Pop);
}


/*========================================================================*/
int PopInit(PyMOLGlobals * G)
{
  register CPop *I = NULL;
  if((I = (G->Pop = Calloc(CPop, 1)))) {

    I->Block = OrthoNewBlock(G, NULL);
    I->Block->fReshape = PopReshape;
    I->Block->active = false;

    I->Block->rect.top = 10;
    I->Block->rect.bottom = 14;
    I->Block->rect.left = 0;
    I->Block->rect.right = 10;

    OrthoAttach(G, I->Block, cOrthoHidden);
    return 1;
  } else
    return 0;
}


/*========================================================================*/
void PopFitBlock(Block * block)
{
  register CPop *I = block->G->Pop;
  int delta;

  if((block->rect.bottom - cPopMargin) < (I->Block->rect.bottom)) {
    delta = (I->Block->rect.bottom - block->rect.bottom) + cPopMargin;
    block->rect.top += delta;
    block->rect.bottom += delta;
  }

  if((block->rect.right + cPopMargin) > (I->Block->rect.right)) {
    delta = (block->rect.right - (I->Block->rect.right)) + cPopMargin;
    block->rect.left -= delta;
    block->rect.right -= delta;
  }

  if((block->rect.left - cPopMargin) < (I->Block->rect.left)) {
    delta = (I->Block->rect.left - block->rect.left) + cPopMargin;
    block->rect.right += delta;
    block->rect.left += delta;
  }

  if((block->rect.top + cPopMargin) > (I->Block->rect.top)) {
    delta = (block->rect.top - (I->Block->rect.top)) + cPopMargin;
    block->rect.top -= delta;
    block->rect.bottom -= delta;
  }
}


/*========================================================================*/
int PopPlaceChild(Block * block, int left_x, int right_x, int row_y, int affinity)
{

  int width = block->rect.right - block->rect.left;
  int height = block->rect.top - block->rect.bottom;
  int target_x;

  block->rect.top = row_y;
  block->rect.bottom = row_y - height;

  if(affinity >= 0) {
    affinity = 1;
    target_x = right_x - 2;
    block->rect.left = target_x;
    block->rect.right = target_x + width;
  } else {
    affinity = -1;
    target_x = left_x - width + 2;
    block->rect.left = target_x;
    block->rect.right = target_x + width;
  }
  PopFitBlock(block);
  if(affinity >= 0) {
    if(block->rect.left != target_x) {
      affinity = -1;
      target_x = left_x - width + 2;
      block->rect.left = target_x;
      block->rect.right = target_x + width;
      PopFitBlock(block);
    }
  } else {
    if(block->rect.left != target_x) {
      affinity = 1;
      target_x = right_x - 2;
      block->rect.left = target_x;
      block->rect.right = target_x + width;
      PopFitBlock(block);
    }
  }
  return affinity;
}
