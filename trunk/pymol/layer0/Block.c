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
#include<GL/gl.h>

#include"Block.h"
#include"main.h"

void BlockInit(Block *I)
{
  I->BackColor[0]=0.2;
  I->BackColor[1]=0.2;
  I->BackColor[2]=0.2;
  I->TextColor[0]=1.0;
  I->TextColor[1]=1.0;
  I->TextColor[2]=1.0;
}
/*========================================================================*/
void BlockFill(Block *I) 
{
  if(PMGUI) {
    glBegin(GL_POLYGON);
    glVertex2i(I->rect.right,I->rect.top);
    glVertex2i(I->rect.right,I->rect.bottom);
    glVertex2i(I->rect.left,I->rect.bottom);
    glVertex2i(I->rect.left,I->rect.top);
    glEnd();
  }
}
/*========================================================================*/
void BlockOutline(Block *I) 
{
  if(PMGUI) {
    glBegin(GL_LINE_LOOP);
    glVertex2i(I->rect.right,I->rect.top);
    glVertex2i(I->rect.right,I->rect.bottom);
    glVertex2i(I->rect.left,I->rect.bottom);
    glVertex2i(I->rect.left,I->rect.top);
    glEnd();
  }
}
/*========================================================================*/
void BlockSetMargin(Block *block,int t,int l,int b,int r)
{
  block->margin.top=t;
  block->margin.left=l;
  block->margin.bottom=b;
  block->margin.right=r;
}
/*========================================================================*/
void BlockReshape(Block *I,int width, int height)
{
  I->rect.top = (height-I->margin.top)-1;
  I->rect.left = I->margin.left;
  I->rect.bottom = I->margin.bottom;
  I->rect.right = (width-I->margin.right)-1;
}
/*========================================================================*/
void BlockRecursiveDraw(Block *block)
{
  if(block)
	 {
		if(block->next)
		  BlockRecursiveDraw(block->next);
		if(block->active)
		  {
			 if(block->fDraw)
				block->fDraw(block);
			 if(block->inside)
				BlockRecursiveDraw(block->inside);
		  }
	 }
}
/*========================================================================*/
Block *BlockRecursiveFind(Block *block,int x,int y)
{
  Block *check;
  if(block){
	 if(!block->active)
		block = BlockRecursiveFind(block->next,x,y);
	 else if ( ! (( block->rect.top    >  y ) &&
					  ( block->rect.bottom <= y ) &&
					  ( block->rect.left   <= x ) &&
					  ( block->rect.right  >  x )))
		block = BlockRecursiveFind(block->next,x,y);
	 else if(block->inside)
		if((check = BlockRecursiveFind(block->inside,x,y)))
		  block = check;
  }
  return(block);
}


