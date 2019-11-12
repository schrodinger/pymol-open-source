
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

struct CPop : public Block {
  CPop(PyMOLGlobals * G) : Block(G){}

  void reshape(int width, int height) override;
};


/*========================================================================*/
void CPop::reshape(int width, int height)
{
  rect.top = height;
  rect.right = width;
}


/*========================================================================*/
Block *PopGetBlock(PyMOLGlobals * G)
{
  CPop *I = G->Pop;
  {
    return (I);
  }
}


/*========================================================================*/
void PopFree(PyMOLGlobals * G)
{
  DeleteP(G->Pop);
}


/*========================================================================*/
int PopInit(PyMOLGlobals * G)
{
  CPop *I = NULL;
  if((I = (G->Pop = new CPop(G)))) {

    I->active = false;

    I->rect.top = 10;
    I->rect.bottom = 14;
    I->rect.left = 0;
    I->rect.right = 10;

    OrthoAttach(G, I, cOrthoHidden);
    return 1;
  } else
    return 0;
}


/*========================================================================*/
void PopFitBlock(Block * block)
{
  CPop *I = block->m_G->Pop; // TODO: Three indirections for a 'this' lol
  int delta;

  if((block->rect.bottom - cPopMargin) < (I->rect.bottom)) {
    delta = (I->rect.bottom - block->rect.bottom) + cPopMargin;
    block->rect.top += delta;
    block->rect.bottom += delta;
  }

  if((block->rect.right + cPopMargin) > (I->rect.right)) {
    delta = (block->rect.right - (I->rect.right)) + cPopMargin;
    block->rect.left -= delta;
    block->rect.right -= delta;
  }

  if((block->rect.left - cPopMargin) < (I->rect.left)) {
    delta = (I->rect.left - block->rect.left) + cPopMargin;
    block->rect.right += delta;
    block->rect.left += delta;
  }

  if((block->rect.top + cPopMargin) > (I->rect.top)) {
    delta = (block->rect.top - (I->rect.top)) + cPopMargin;
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
