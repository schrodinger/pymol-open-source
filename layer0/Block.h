

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
#ifndef _H_Block
#define _H_Block

#include "PyMOLGlobals.h"

struct BlockRect{
  int top;
  int left;
  int bottom;
  int right;
};

struct Block {
  PyMOLGlobals * G;
  Block *next = nullptr, *inside = nullptr, *parent = nullptr;
  void *reference = nullptr;
  BlockRect rect {}, margin {};
  bool active = false;
  float BackColor[3] = {0.2F, 0.2F, 0.2F };
  float TextColor[3] = {1.0F, 1.0F, 1.0F };

  virtual void draw(CGO *orthoCGO) {};
  virtual bool fastDraw(CGO *orthoCGO) { return false; }
  virtual void reshape(int width, int height);
  virtual int click(int button, int x, int y, int mod) { return 0; }
  virtual int cursor (int x, int y, int mod) { return 0; }
  virtual int drag (int x, int y, int mod) { return 0; }
  virtual int release (int button, int x, int y, int mod) { return 0; }
  Block(PyMOLGlobals * _G) : G(_G){}

  void globalToLocal(int x, int y, int *lx, int *ly);
  void recursiveDraw(CGO *orthoCGO);
  bool recursiveFastDraw(CGO *orthoCGO);
  Block *recursiveFind(int x, int y);
  void setMargin(int t, int l, int b, int r);
  void fill(CGO *orthoCGO);
  int getWidth() const;
  int getHeight() const;
  void translate(int dx, int dy);
  void drawLeftEdge(CGO *orthoCGO);
  void drawTopEdge();
  bool rectXYInside(int x, int y) const;
};

typedef Block **CBlock;

#endif
