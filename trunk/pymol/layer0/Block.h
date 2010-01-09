

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

typedef struct {
  int top, left, bottom, right;
} BlockRect;

typedef struct Block {
  PyMOLGlobals *G;
  struct Block *next, *inside, *parent;
  void *reference;
  BlockRect rect, margin;
  int active;
  float BackColor[3];
  float TextColor[3];
  void (*fDraw) (struct Block * block);
  void (*fReshape) (struct Block * block, int width, int height);
  int (*fClick) (struct Block * block, int button, int x, int y, int mod);
  int (*fCursor) (struct Block * block, int x, int y, int mod);
  int (*fDrag) (struct Block * block, int x, int y, int mod);
  int (*fRelease) (struct Block * block, int button, int x, int y, int mod);
  int (*fTranslate) (struct Block * block, int dx, int dy);
} Block;

typedef Block **CBlock;

void BlockGlobalToLocal(Block * block, int x, int y, int *lx, int *ly);
void BlockRecursiveDraw(Block * block);
Block *BlockRecursiveFind(Block * block, int x, int y);
void BlockSetMargin(Block * block, int t, int l, int b, int r);
void BlockReshape(Block * block, int width, int height);
void BlockFill(Block * I);
void BlockDrawLeftEdge(Block * I);
void BlockGetSize(Block * I, int *width, int *height);
void BlockOutline(Block * I);
void BlockInit(PyMOLGlobals * G, Block * I);
void BlockTranslate(Block * I, int dx, int dy);
void BlockDrawTopEdge(Block * I);
int BlockRectXYInside(BlockRect *rect, int x, int y);
#endif
