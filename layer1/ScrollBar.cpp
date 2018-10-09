
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
#include"CGO.h"

struct CScrollBar : public Block {
  int HorV {};
  float BackColor[3] = { 0.1f, 0.1f, 0.1f };
  float BarColor[3] = { 0.5f, 0.5f, 0.5f };
  int ListSize { 10 };
  int DisplaySize { 7 };
  int BarSize {};
  int StartPos {};
  float ExactBarSize {};
  float Value { 0.0f };
  float StartValue {};
  float ValueMax { 0.0f };
  int BarRange {};
  int BarMin {};
  int BarMax {};
  int Grabbed {};

  CScrollBar(PyMOLGlobals * G) : Block(G){}

  virtual int release(int button, int x, int y, int mod) override;
  virtual int click(int button, int x, int y, int mod) override;
  virtual int drag(int x, int y, int mod) override;
  virtual void draw(CGO *orthoCGO) override;
  virtual void reshape(int width, int height) override;
};

void ScrollBarMaxOut(struct CScrollBar *I)
{
  I->Value = I->ValueMax;
}

int ScrollBarIsMaxed(struct CScrollBar *I)
{
  if(I->ValueMax > 0.0F) {
    if(I->Value == I->ValueMax)
      return true;
    else
      return false;
  } else
    return false;
}

static void ScrollBarUpdate(struct CScrollBar *I)
{
  int range;

  if(I->HorV) {
    range = (I->rect.right - I->rect.left);
  } else {
    range = (I->rect.top - I->rect.bottom);
  }
  I->ExactBarSize = (range * I->DisplaySize) / (float) I->ListSize;
  I->BarSize = (int) (0.499F + I->ExactBarSize);
  if(I->BarSize < 4)
    I->BarSize = DIP2PIXEL(4);
  I->BarRange = range - I->BarSize;
  if(I->BarRange < 2)
    I->BarRange = 2;
  I->ValueMax = (float) I->ListSize - I->DisplaySize;
  if(I->ValueMax < 1)
    I->ValueMax = 1;
  if(I->Value > I->ValueMax)
    I->Value = (float) I->ValueMax;
  else if(I->Value < 0.0)
    I->Value = 0.0F;
}

void ScrollBarFill(struct CScrollBar *I ORTHOCGOARG)
{
  if (orthoCGO)
    CGOColorv(orthoCGO, I->BackColor);
  else
    glColor3fv(I->BackColor);
  I->fill(orthoCGO);
}

void CScrollBar::draw(CGO* orthoCGO)
{
  ScrollBarDrawImpl(this, true, orthoCGO);
}

void ScrollBarDrawImpl(Block * block, short fill  ORTHOCGOARG)
{
  PyMOLGlobals *G = block->G;
  float value;
  int top, left, bottom, right;

  CScrollBar *I = (CScrollBar *) block->reference;

  if (fill)
    ScrollBarFill(I ORTHOCGOARGVAR);

  ScrollBarUpdate(I);

  value = I->Value;
  if(value > I->ValueMax)
    value = I->ValueMax;

  if(I->HorV) {
    top = block->rect.top - 1;
    bottom = block->rect.bottom + 1;
    left = (int) (0.499F + block->rect.left + (I->BarRange * value) / I->ValueMax);
    right = left + I->BarSize;
    I->BarMin = left;
    I->BarMax = right;
  } else {
    top = (int) (0.499F + block->rect.top - (I->BarRange * value) / I->ValueMax);
    bottom = top - I->BarSize;
    left = block->rect.left + 1;
    right = block->rect.right - 1;
    I->BarMin = top;
    I->BarMax = bottom;
  }

  if(G->HaveGUI && G->ValidContext) {
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
      CGOColorv(orthoCGO, I->BarColor);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right - 1, top - 1, 0.f);
      CGOVertex(orthoCGO, right - 1, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
      CGOVertex(orthoCGO, left + 1, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
    } else {
      glColor3fv(I->BarColor);
      glBegin(GL_POLYGON);
      glVertex2i(right - 1, top - 1);
      glVertex2i(right - 1, bottom + 1);
      glVertex2i(left + 1, bottom + 1);
      glVertex2i(left + 1, top - 1);
      glEnd();
    }
  }
}

void ScrollBarDrawHandle(struct CScrollBar *I, float alpha ORTHOCGOARG)
{
  float value;
  int top, left, bottom, right;
  Block *block = I; // TODO: Remove during ScrollBar refactor
  PyMOLGlobals *G = block->G;

  value = I->Value;
  if(value > I->ValueMax)
    value = I->ValueMax;

  if(I->HorV) {
    top = block->rect.top - 1;
    bottom = block->rect.bottom + 1;
    left = (int) (0.499F + block->rect.left + (I->BarRange * value) / I->ValueMax);
    right = left + I->BarSize;
  } else {
    top = (int) (0.499F + block->rect.top - (I->BarRange * value) / I->ValueMax);
    bottom = top - I->BarSize;
    left = block->rect.left + 1;
    right = block->rect.right - 1;
  }

  if(G->HaveGUI && G->ValidContext) {

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
      CGOColor(orthoCGO, I->BarColor[0], I->BarColor[1], I->BarColor[2]);
      CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
      CGOVertex(orthoCGO, right - 1, top - 1, 0.f);
      CGOVertex(orthoCGO, right - 1, bottom + 1, 0.f);
      CGOVertex(orthoCGO, left + 1, top - 1, 0.f);
      CGOVertex(orthoCGO, left + 1, bottom + 1, 0.f);
      CGOEnd(orthoCGO);
      CGOAlpha(orthoCGO, 1.f);
    } else {
      glColor4f(I->BarColor[0], I->BarColor[1], I->BarColor[2], alpha);
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

void ScrollBarSetValueNoCheck(struct CScrollBar *I, float value)
{
  I->Value = value;
}

void ScrollBarSetValue(struct CScrollBar *I, float value)
{
  I->Value = value > I->ValueMax ? I->ValueMax :
             value < 0.0 ? 0.0 : value;
}

void ScrollBarMoveBy(struct CScrollBar *I, float value) {
  ScrollBarSetValue(I, I->Value + value);
}

float ScrollBarGetValue(struct CScrollBar *I)
{
  return (I->Value);
}

void CScrollBar::reshape(int width, int height)
{
}

int ScrollBarGrabbed(struct CScrollBar *I)
{
  return OrthoGrabbedBy(I->G, I);
}

int CScrollBar::click(int button, int x, int y, int mod)
{
  CScrollBar *I = (CScrollBar *) reference;
  int grab = 0;

  if(button == P_GLUT_MIDDLE_BUTTON) {
    if(I->HorV) {
      if(x < I->BarMin || x > I->BarMax)
        ScrollBarSetValue(I, (I->ListSize * (x - rect.left)) /
            (rect.right - rect.left) - I->DisplaySize * 0.5F);
      grab = x;
    } else {
      if(y > I->BarMin || y < I->BarMax)
        ScrollBarSetValue(I, (I->ListSize * (y - rect.top)) /
            (rect.bottom - rect.top) - I->DisplaySize * 0.5F);
      grab = y;
    }
  } else {
    if(I->HorV) {
      if(x > I->BarMax) {
        I->Value += I->DisplaySize;
      } else if(x < I->BarMin) {
        I->Value -= I->DisplaySize;
      } else {
        grab = x;
      }
    } else {
      if(y > I->BarMin) {
        I->Value -= I->DisplaySize;
      } else if(y < I->BarMax) {
        I->Value += I->DisplaySize;
      } else {
        grab = y;
      }
    }
  }

  if(grab) {
    OrthoGrab(G, this);
    I->StartPos = grab;
    I->StartValue = I->Value;
  }

  OrthoDirty(G);
  return 0;
}

int CScrollBar::drag(int x, int y, int mod)
{
  CScrollBar *I = (CScrollBar *) reference;
  int displ;
  if(I->HorV)
    displ = I->StartPos - x;
  else
    displ = y - I->StartPos;
  ScrollBarSetValue(I, I->StartValue - (I->ValueMax * displ) / I->BarRange);
  OrthoDirty(G);
  return true;
}

int CScrollBar::release(int button, int x, int y, int mod)
{
  OrthoUngrab(G);
  OrthoDirty(G);
  return 0;
}

Block *ScrollBarGetBlock(struct CScrollBar * I)
{
  return (I);
}

void ScrollBarSetLimits(struct CScrollBar *I, int list_size, int display_size)
{
  I->ListSize = list_size;
  I->DisplaySize = display_size;
  ScrollBarUpdate(I);
}

void ScrollBarSetBox(struct CScrollBar *I, int top, int left, int bottom, int right)
{
  I->rect.top = top;
  I->rect.left = left;
  I->rect.bottom = bottom;
  I->rect.right = right;
}

void ScrollBarDoDraw(struct CScrollBar *I ORTHOCGOARG)
{
    I->draw(orthoCGO);
}

void ScrollBarDoDrawNoFill(struct CScrollBar *I ORTHOCGOARG)
{
  ScrollBarDrawImpl(I, false, orthoCGO);
}

void ScrollBarDoRelease(struct CScrollBar *I, int button, int x, int y, int mod)
{
    I->release(button, x, y, mod);
}

void ScrollBarDoDrag(struct CScrollBar *I, int x, int y, int mod)
{
    I->drag(x, y, mod);
}

void ScrollBarDoClick(struct CScrollBar *I, int button, int x, int y, int mod)
{
    I->click(button, x, y, mod);
}

struct CScrollBar *ScrollBarNew(PyMOLGlobals * G, int horizontal)
{
  CScrollBar *I = new CScrollBar(G);

  I->active = false;
  I->reference = (void *) I;
  I->HorV = horizontal;
  return (I);
}

void ScrollBarFree(struct CScrollBar *I)
{
  DeleteP(I);
}
