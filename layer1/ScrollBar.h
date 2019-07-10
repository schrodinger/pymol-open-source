
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
#ifndef _H_ScrollBar
#define _H_ScrollBar

#include "Block.h"
#include "Ortho.h"
#include "PyMOLGlobals.h"
#include "pymol/algorithm.h"

class ScrollBar : public Block {
private:
  bool m_HorV {};
  float m_BarColor[3] = { 0.5f, 0.5f, 0.5f };
  int m_ListSize { 10 };
  int m_DisplaySize { 7 };
  int m_BarSize {};
  int m_StartPos {};
  float m_ExactBarSize {};
  float m_StartValue {};
  int m_BarRange {};
  int m_BarMin {};
  int m_BarMax {};
  int Grabbed {};
  float m_Value { 0.0f };
  float m_ValueMax { 0.0f };

  void drawImpl(bool fill, CGO *orthoCGO);
  void update();

public:
  ScrollBar(PyMOLGlobals * G, int horizontal)
  : Block(G)
  , m_HorV(horizontal){
    BackColor[0] = 0.1f;
    BackColor[1] = 0.1f;
    BackColor[2] = 0.1f;
  }

  virtual int release(int button, int x, int y, int mod) override;
  virtual int click(int button, int x, int y, int mod) override;
  virtual int drag(int x, int y, int mod) override;
  virtual void draw(CGO *orthoCGO) override { drawImpl(true, orthoCGO); }
  virtual void reshape(int width, int height) override {};

  void maxOut() { m_Value = m_ValueMax; };
  bool isMaxed() const;
  void fill(CGO *orthoCGO);
  void setValueNoCheck(float _value) { m_Value = _value; }
  void setValue(float _value){ m_Value = pymol::clamp(_value, 0.0f, m_ValueMax); }
  void drawNoFill(CGO *orthoCGO) { drawImpl(false, orthoCGO); }
  int grabbed() { return OrthoGrabbedBy(m_G, this); }
  float getValue() const { return m_Value; }
  void setLimits(int list_size, int display_size);
  void moveBy(float value) { setValue(m_Value + value); }
  void drawHandle(float alpha, CGO *orthoCGO);
  void setBox(int top, int left, int bottom, int right);
};
#endif
