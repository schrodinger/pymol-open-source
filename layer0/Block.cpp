

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
#include"CGO.h"

int Block::getWidth() const {
    return rect.right - rect.left;
}

/*========================================================================*/
int Block::getHeight() const {
    return rect.top - rect.bottom;
}

/*========================================================================*/
void Block::fill(CGO *orthoCGO)
{
  if(m_G->HaveGUI && m_G->ValidContext) {
    if (orthoCGO){
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.right, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.right, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.bottom, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glBegin(GL_POLYGON);
      glVertex2i(rect.right, rect.top);
      glVertex2i(rect.right, rect.bottom);
      glVertex2i(rect.left, rect.bottom);
      glVertex2i(rect.left, rect.top);
      glEnd();
    }
  }
}


/*========================================================================*/
void Block::drawLeftEdge(CGO *orthoCGO)
{
  if(m_G->HaveGUI && m_G->ValidContext) {
    if (orthoCGO){
      CGOColor(orthoCGO, .3f, .3f, .3f);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, rect.left, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left + 1.f, rect.bottom, 0.f);
      CGOVertex(orthoCGO, rect.left, rect.top, 0.f);
      CGOVertex(orthoCGO, rect.left + 1.f, rect.top, 0.f);
      CGOEnd(orthoCGO);
    } else {
      if(m_G->HaveGUI && m_G->ValidContext) {
	glColor3f(0.3, 0.3, 0.3);
	glBegin(GL_LINES);
	glVertex2i(rect.left, rect.bottom);
	glVertex2i(rect.left, rect.top);
	glEnd();
      }
    }
  }
}

/*========================================================================*/
void Block::drawTopEdge()
{
#ifndef PURE_OPENGL_ES_2
  if(m_G->HaveGUI && m_G->ValidContext) {
    glColor3f(0.3, 0.3, 0.3);
    glBegin(GL_LINES);
    glVertex2i(rect.right, rect.top);
    glVertex2i(rect.left, rect.top);
    glEnd();
  }
#endif
}


/*========================================================================*/
void Block::setMargin(int t, int l, int b, int r)
{
  margin.top = t;
  margin.left = l;
  margin.bottom = b;
  margin.right = r;
}

/*========================================================================*/
void Block::reshape(int width, int height)
{
  rect.top = (height - margin.top);
  rect.left = margin.left;
  rect.bottom = margin.bottom;
  rect.right = (width - margin.right);
}

/*========================================================================*/
void Block::translate(int dx, int dy)
{
  rect.top += dy;
  rect.left += dx;
  rect.bottom += dy;
  rect.right += dx;
}


/*========================================================================*/
void Block::recursiveDraw(CGO *orthoCGO)
{
  if (this->next)
    next->recursiveDraw(orthoCGO);
  if (active) {
      draw(orthoCGO);
    if (inside)
      inside->recursiveDraw(orthoCGO);
  }
}

/*========================================================================*/
bool Block::recursiveFastDraw(CGO *orthoCGO)
{
  bool ret = false;
  if (next)
    ret |= next->recursiveFastDraw(orthoCGO);
  if (active) {
      ret |= this->fastDraw(orthoCGO);
    if (inside)
      ret |= inside->recursiveFastDraw(orthoCGO);
  }
  return ret;
}

/*========================================================================*/
Block *Block::recursiveFind(int x, int y)
{
  Block *check;
  Block *block = this;
  if(block) {
    if(!block->active)
      block = block->next->recursiveFind(x, y);
    else if(!rectXYInside(x, y))
      block = block->next->recursiveFind(x, y);
    else if(block->inside)
      if((check = block->inside->recursiveFind(x, y)))
        block = check;
  }
  return (block);
}

bool Block::rectXYInside(int x, int y) const
{
  return ((y <= rect.top) && (y >= rect.bottom) &&
          (x <= rect.right) && (x >= rect.left));
}
