

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
  int top = 0;
  int left = 0;
  int bottom = 0;
  int right = 0;
};

struct Block {
  PyMOLGlobals * G;
  Block *next = nullptr, *inside = nullptr, *parent = nullptr;
  void *reference = nullptr;
  BlockRect rect, margin;
  bool active = false;
  float BackColor[3] = {0.2F, 0.2F, 0.2F };
  float TextColor[3] = {1.0F, 1.0F, 1.0F };

  void (*fDraw) (Block * block, CGO *orthoCGO) = nullptr;
  short (*fFastDraw) (Block * block, CGO *orthoCGO) = nullptr;
  void (*fReshape) (Block * block, int width, int height) = nullptr;
  int (*fClick) (Block * block, int button, int x, int y, int mod) = nullptr;
  int (*fCursor) (Block * block, int x, int y, int mod) = nullptr;
  int (*fDrag) (Block * block, int x, int y, int mod) = nullptr;
  int (*fRelease) (Block * block, int button, int x, int y, int mod) = nullptr;
  int (*fTranslate) (Block * block, int dx, int dy) = nullptr;

  Block(PyMOLGlobals* _G) : G(_G){};
  Block(const Block&) = default;
  Block& operator=(const Block&) = default;
  Block(Block&&) = default;
  Block& operator=(Block&&) = default;

  void globalToLocal(int x, int y, int *lx, int *ly);
  void recursiveDraw(CGO *orthoCGO);
  bool recursiveFastDraw(CGO *orthoCGO);
  Block *recursiveFind(int x, int y);
  void setMargin(int t, int l, int b, int r);
  void reshape(int width, int height);
  void fill(CGO *orthoCGO);
  int getWidth() const;
  int getHeight() const;
  void translate(int dx, int dy);
  void drawLeftEdge(CGO *orthoCGO);
  void drawTopEdge();
  bool rectXYInside(int x, int y) const;
};

typedef Block **CBlock;

void BlockReshape(Block * block, int width, int height);
#endif
