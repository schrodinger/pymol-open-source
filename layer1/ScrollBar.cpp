
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
#include"Ortho.h"
#include"ScrollBar.h"
#include"CGO.h"

bool ScrollBar::isMaxed() const
{
  if(m_ValueMax > 0.0F) {
    return m_Value == m_ValueMax;
  } else
    return false;
}

void ScrollBar::update()
{
  int range;

  if(m_HorV) {
    range = (rect.right - rect.left);
  } else {
    range = (rect.top - rect.bottom);
  }
  m_ExactBarSize = (range * m_DisplaySize) / static_cast<float>(m_ListSize);
  m_BarSize = static_cast<int> (0.499F + m_ExactBarSize);
  if(m_BarSize < 4)
    m_BarSize = DIP2PIXEL(4);
  m_BarRange = range - m_BarSize;
  if(m_BarRange < 2)
    m_BarRange = 2;
  m_ValueMax = static_cast<float>(m_ListSize - m_DisplaySize);
  if(m_ValueMax < 1)
    m_ValueMax = 1;
  m_Value = pymol::clamp(m_Value, 0.0f, m_ValueMax);
}

void ScrollBar::fill(CGO* orthoCGO)
{
  if (orthoCGO)
    CGOColorv(orthoCGO, BackColor);
#ifndef PURE_OPENGL_ES_2
  else
    glColor3fv(BackColor);
#endif
  Block::fill(orthoCGO);
}

void ScrollBar::drawImpl(bool bFill, CGO* orthoCGO)
{
  int top, left, bottom, right;

  if (bFill)
    fill(orthoCGO);

  update();

  float value = std::min(m_Value, m_ValueMax);

  if(m_HorV) {
    top = rect.top - 1;
    bottom = rect.bottom + 1;
    left = (int) (0.499F + rect.left + (m_BarRange * value) / m_ValueMax);
    right = left + m_BarSize;
    m_BarMin = left;
    m_BarMax = right;
  } else {
    top = (int) (0.499F + rect.top - (m_BarRange * value) / m_ValueMax);
    bottom = top - m_BarSize;
    left = rect.left + 1;
    right = rect.right - 1;
    m_BarMin = top;
    m_BarMax = bottom;
  }

  if(m_G->HaveGUI && m_G->ValidContext) {

    if (orthoCGO){
      CGOColor(orthoCGO, 0.8F, 0.8F, 0.8F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right, top, 0.f);
      CGOVertex(orthoCGO, right, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left, top, 0.f);
      CGOVertex(orthoCGO, left, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.8F, 0.8F, 0.8F);
      glBegin(GL_POLYGON);
      glVertex2i(right, top);
      glVertex2i(right, bottom + 1);
      glVertex2i(left, bottom + 1);
      glVertex2i(left, top);
      glEnd();
    }

    if (orthoCGO){
      CGOColor(orthoCGO, 0.3F, 0.3F, 0.3F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right, top - 1, 0.f);
      CGOVertex(orthoCGO, right, bottom, 0.f);
      CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
      CGOVertex(orthoCGO, left + 1, bottom, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.3F, 0.3F, 0.3F);
      glBegin(GL_POLYGON);
      glVertex2i(right, top - 1);
      glVertex2i(right, bottom);
      glVertex2i(left + 1, bottom);
      glVertex2i(left + 1, top - 1);
      glEnd();
    }

    if (orthoCGO){
      CGOColor(orthoCGO, 0.3F, 0.3F, 0.3F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right, bottom + 1, 0.f);
      CGOVertex(orthoCGO, right, bottom, 0.f);
      CGOVertex(orthoCGO, left, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left, bottom, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3f(0.3F, 0.3F, 0.3F);
      glBegin(GL_POLYGON);
      glVertex2i(right, bottom + 1);
      glVertex2i(right, bottom);
      glVertex2i(left, bottom);
      glVertex2i(left, bottom + 1);
      glEnd();
    }

    if (orthoCGO){
      CGOColorv(orthoCGO, m_BarColor);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right - 1, top - 1, 0.f);
      CGOVertex(orthoCGO, right - 1, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
      CGOVertex(orthoCGO, left + 1, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3fv(m_BarColor);
      glBegin(GL_POLYGON);
      glVertex2i(right - 1, top - 1);
      glVertex2i(right - 1, bottom + 1);
      glVertex2i(left + 1, bottom + 1);
      glVertex2i(left + 1, top - 1);
      glEnd();
    }
  }
}

void ScrollBar::drawHandle(float alpha, CGO* orthoCGO)
{
  float value;
  int top, left, bottom, right;

  value = std::min(m_Value, m_ValueMax);

  if(m_HorV) {
    top = rect.top - 1;
    bottom = rect.bottom + 1;
    left = (int) (0.499F + rect.left + (m_BarRange * value) / m_ValueMax);
    right = left + m_BarSize;
  } else {
    top = (int) (0.499F + rect.top - (m_BarRange * value) / m_ValueMax);
    bottom = top - m_BarSize;
    left = rect.left + 1;
    right = rect.right - 1;
  }

  if(m_G->HaveGUI && m_G->ValidContext) {

    glEnable(GL_BLEND);
    if (orthoCGO){
      CGOAlpha(orthoCGO, alpha);
      CGOColor(orthoCGO, 0.8F, 0.8F, 0.8F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right, top, 0.f);
      CGOVertex(orthoCGO, right, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left, top, 0.f);
      CGOVertex(orthoCGO, left, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
      CGOAlpha(orthoCGO, 1.f);
    } else {
      glColor4f(0.8F, 0.8F, 0.8F, alpha);
      glBegin(GL_POLYGON);
      glVertex2i(right, top);
      glVertex2i(right, bottom + 1);
      glVertex2i(left, bottom + 1);
      glVertex2i(left, top);
      glEnd();
    }

    if (orthoCGO){
      CGOAlpha(orthoCGO, alpha);
      CGOColor(orthoCGO, 0.3F, 0.3F, 0.3F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right, top - 1, 0.f);
      CGOVertex(orthoCGO, right, bottom, 0.f);
      CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
      CGOVertex(orthoCGO, left + 1, bottom, 0.f);
      CGOEnd(orthoCGO);
      CGOAlpha(orthoCGO, 1.f);
    } else {
      glColor4f(0.3F, 0.3F, 0.3F, alpha);
      glBegin(GL_POLYGON);
      glVertex2i(right, top - 1);
      glVertex2i(right, bottom);
      glVertex2i(left + 1, bottom);
      glVertex2i(left + 1, top - 1);
      glEnd();
    }

    if (orthoCGO){
      CGOAlpha(orthoCGO, alpha);
      CGOColor(orthoCGO, 0.3F, 0.3F, 0.3F);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right, bottom + 1, 0.f);
      CGOVertex(orthoCGO, right, bottom, 0.f);
      CGOVertex(orthoCGO, left, bottom, 0.f);
      CGOVertex(orthoCGO, left, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
      CGOAlpha(orthoCGO, 1.f);
    } else {
      glColor4f(0.3F, 0.3F, 0.3F, alpha);
      glBegin(GL_POLYGON);
      glVertex2i(right, bottom + 1);
      glVertex2i(right, bottom);
      glVertex2i(left, bottom);
      glVertex2i(left, bottom + 1);
      glEnd();
    }

    if (orthoCGO){
      CGOAlpha(orthoCGO, alpha);
      CGOColor(orthoCGO, m_BarColor[0], m_BarColor[1], m_BarColor[2]);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right - 1, top - 1, 0.f);
      CGOVertex(orthoCGO, right - 1, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
      CGOVertex(orthoCGO, left + 1, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
      CGOAlpha(orthoCGO, 1.f);
    } else {
      glColor4f(m_BarColor[0], m_BarColor[1], m_BarColor[2], alpha);
      glBegin(GL_POLYGON);
      glVertex2i(right - 1, top - 1);
      glVertex2i(right - 1, bottom + 1);
      glVertex2i(left + 1, bottom + 1);
      glVertex2i(left + 1, top - 1);
      glEnd();
    }
    glDisable(GL_BLEND);
  }
}

int ScrollBar::click(int button, int x, int y, int mod)
{
  int grab = 0;

  if(button == P_GLUT_MIDDLE_BUTTON) {
    if(m_HorV) {
      if(x < m_BarMin || x > m_BarMax)
        setValue((m_ListSize * (x - rect.left)) /
            (rect.right - rect.left) - m_DisplaySize * 0.5F);
      grab = x;
    } else {
      if(y > m_BarMin || y < m_BarMax)
        setValue((m_ListSize * (y - rect.top)) /
            (rect.bottom - rect.top) - m_DisplaySize * 0.5F);
      grab = y;
    }
  } else {
    if(m_HorV) {
      if(x > m_BarMax) {
        m_Value += m_DisplaySize;
      } else if(x < m_BarMin) {
        m_Value -= m_DisplaySize;
      } else {
        grab = x;
      }
    } else {
      if(y > m_BarMin) {
        m_Value -= m_DisplaySize;
      } else if(y < m_BarMax) {
        m_Value += m_DisplaySize;
      } else {
        grab = y;
      }
    }
  }

  if(grab) {
    OrthoGrab(m_G, this);
    m_StartPos = grab;
    m_StartValue = m_Value;
  }

  OrthoDirty(m_G);
  return 0;
}

int ScrollBar::drag(int x, int y, int mod)
{
  int displ;
  if(m_HorV)
    displ = m_StartPos - x;
  else
    displ = y - m_StartPos;
  setValue(m_StartValue - (m_ValueMax * displ) / m_BarRange);
  OrthoDirty(m_G);
  return true;
}

int ScrollBar::release(int button, int x, int y, int mod)
{
  OrthoUngrab(m_G);
  OrthoDirty(m_G);
  return 0;
}

void ScrollBar::setLimits(int list_size, int display_size)
{
  m_ListSize = list_size;
  m_DisplaySize = display_size;
  update();
}

void ScrollBar::setBox(int top, int left, int bottom, int right)
{
  rect.top = top;
  rect.left = left;
  rect.bottom = bottom;
  rect.right = right;
}

