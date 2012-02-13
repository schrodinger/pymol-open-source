

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#include"os_python.h"
#include"os_predef.h"
#include"os_gl.h"

#include"Block.h"
#include"main.h"

void BlockGetSize(Block * I, int *width, int *height)
{
  *width = I->rect.right - I->rect.left;
  *height = I->rect.top - I->rect.bottom;
}


/*========================================================================*/
void BlockInit(PyMOLGlobals * G, Block * I)
{
  I->G = G;
  I->BackColor[0] = 0.2F;
  I->BackColor[1] = 0.2F;
  I->BackColor[2] = 0.2F;
  I->TextColor[0] = 1.0F;
  I->TextColor[1] = 1.0F;
  I->TextColor[2] = 1.0F;
}


/*========================================================================*/
void BlockFill(Block * I)
{
  register PyMOLGlobals *G = I->G;
  if(G->HaveGUI && G->ValidContext) {
#if defined(PURE_OPENGL_ES_2)
    /* TODO */
#else
#ifdef _PYMOL_GL_DRAWARRAYS
    {
      const GLint polyVerts[] = {
	I->rect.right, I->rect.top,
	I->rect.right, I->rect.bottom,
	I->rect.left, I->rect.top,
	I->rect.left, I->rect.bottom
      };
#if defined(OPENGL_ES_1)      
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(2, GL_INT, 0, polyVerts);
      glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
      glDisableClientState(GL_VERTEX_ARRAY);
#elif defined(OPENGL_ES_2)
      glVertexAttribPointer(0, 2, GL_INT, GL_FALSE, 0, polyVerts);
      glEnableVertexAttribArray(0);
      glDrawArrays(GL_TRIANGLES, 0, 4);
      glDisableVertexAttribArray(0);
#endif
    }
#else
    glBegin(GL_POLYGON);
    glVertex2i(I->rect.right, I->rect.top);
    glVertex2i(I->rect.right, I->rect.bottom);
    glVertex2i(I->rect.left, I->rect.bottom);
    glVertex2i(I->rect.left, I->rect.top);
    glEnd();
#endif
#endif
  }
}


/*========================================================================*/
void BlockDrawLeftEdge(Block * I)
{
  register PyMOLGlobals *G = I->G;
#if defined(PURE_OPENGL_ES_2)
    /* TODO */
  {
#else
  if(G->HaveGUI && G->ValidContext) {
    glColor3f(0.3, 0.3, 0.3);
#ifdef _PYMOL_GL_DRAWARRAYS
    {
      const GLint lineVerts[] = {
	I->rect.left, I->rect.bottom,
	I->rect.left, I->rect.top
      };
#if defined(OPENGL_ES_1)
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(2, GL_INT, 0, lineVerts);
      glDrawArrays(GL_LINES, 0, 2);
      glDisableClientState(GL_VERTEX_ARRAY);
#elif defined(OPENGL_ES_2)
      glVertexAttribPointer(0, 2, GL_INT, GL_FALSE, 0, lineVerts);
      glEnableVertexAttribArray(0);
      glDrawArrays(GL_LINES, 0, 2);      
      glDisableVertexAttribArray(0);
#endif
    }
#else
    glBegin(GL_LINES);
    glVertex2i(I->rect.left, I->rect.bottom);
    glVertex2i(I->rect.left, I->rect.top);
    glEnd();
#endif
#endif
  }
}


/*========================================================================*/
void BlockOutline(Block * I)
{
  register PyMOLGlobals *G = I->G;
  if(G->HaveGUI && G->ValidContext) {
#if defined(PURE_OPENGL_ES_2)
    /* TODO */
#else
#ifdef _PYMOL_GL_DRAWARRAYS
    const GLint lineVerts[] = {
      I->rect.right, I->rect.top,
      I->rect.right, I->rect.bottom,
      I->rect.left, I->rect.bottom,
      I->rect.left, I->rect.top
    };
#if defined(OPENGL_ES_1)
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_INT, 0, lineVerts);
    glDrawArrays(GL_LINE_LOOP, 0, 4);
    glDisableClientState(GL_VERTEX_ARRAY);
#elif defined(OPENGL_ES_2)
    glVertexAttribPointer(0, 2, GL_INT, GL_FALSE, 0, lineVerts);
    glEnableVertexAttribArray(0);
    glDrawArrays(GL_LINE_LOOP, 0, 4);
    glDisableVertexAttribArray(0);
#endif

#else
    glBegin(GL_LINE_LOOP);
    glVertex2i(I->rect.right, I->rect.top);
    glVertex2i(I->rect.right, I->rect.bottom);
    glVertex2i(I->rect.left, I->rect.bottom);
    glVertex2i(I->rect.left, I->rect.top);
    glEnd();
#endif
#endif
  }
}


/*========================================================================*/
void BlockDrawTopEdge(Block * I)
{
  register PyMOLGlobals *G = I->G;
#ifdef PURE_OPENGL_ES_2
    /* TODO */
  {
#else
  if(G->HaveGUI && G->ValidContext) {
    glColor3f(0.3, 0.3, 0.3);
#ifdef _PYMOL_GL_DRAWARRAYS
    {
      const GLint lineVerts[] = {
	I->rect.right, I->rect.top,
	I->rect.left, I->rect.top
      };
#if defined(OPENGL_ES_1)
      glEnableClientState(GL_VERTEX_ARRAY);
      glVertexPointer(2, GL_INT, 0, lineVerts);
      glDrawArrays(GL_LINES, 0, 2);
      glDisableClientState(GL_VERTEX_ARRAY);
#elif defined(OPENGL_ES_2)
      glVertexAttribPointer(0, 2, GL_INT, GL_FALSE, 0, lineVerts);
      glEnableVertexAttribArray(0);
      glDrawArrays(GL_LINES, 0, 2);
      glDisableVertexAttribArray(0);
#endif
    }
#else
    glBegin(GL_LINES);
    glVertex2i(I->rect.right, I->rect.top);
    glVertex2i(I->rect.left, I->rect.top);
    glEnd();
#endif
#endif
  }
}


/*========================================================================*/
void BlockSetMargin(Block * block, int t, int l, int b, int r)
{
  block->margin.top = t;
  block->margin.left = l;
  block->margin.bottom = b;
  block->margin.right = r;
}

/*========================================================================*/
void BlockReshape(Block * I, int width, int height)
{
  I->rect.top = (height - I->margin.top);
  I->rect.left = I->margin.left;
  I->rect.bottom = I->margin.bottom;
  I->rect.right = (width - I->margin.right);
}


/*========================================================================*/
void BlockTranslate(Block * I, int dx, int dy)
{
  I->rect.top += dy;
  I->rect.left += dx;
  I->rect.bottom += dy;
  I->rect.right += dx;
}


/*========================================================================*/
void BlockRecursiveDraw(Block * block)
{
  if(block) {
    if(block->next)
      BlockRecursiveDraw(block->next);
    if(block->active) {
      if(block->fDraw)
        block->fDraw(block);
      if(block->inside)
        BlockRecursiveDraw(block->inside);
    }
  }
}


/*========================================================================*/
Block *BlockRecursiveFind(Block * block, int x, int y)
{
  Block *check;
  if(block) {
    if(!block->active)
      block = BlockRecursiveFind(block->next, x, y);
    else if(!((block->rect.top >= y) &&
              (block->rect.bottom <= y) &&
              (block->rect.left <= x) && (block->rect.right >= x)))
      block = BlockRecursiveFind(block->next, x, y);
    else if(block->inside)
      if((check = BlockRecursiveFind(block->inside, x, y)))
        block = check;
  }
  return (block);
}

int BlockRectXYInside(BlockRect *rect, int x, int y)
{
  return ((y <= rect->top) && (y >= rect->bottom) &&
          (x <= rect->right) && (x >= rect->left));
}
