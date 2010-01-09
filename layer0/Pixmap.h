

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
#ifndef _H_Pixmap
#define _H_Pixmap

#include"PyMOLGlobals.h"


/* for the sake of simplicity, all pixmaps are 32-bit RGBA */

typedef struct {
  PyMOLGlobals *G;
  int height, width;
  unsigned char *buffer;
} CPixmap;

void PixmapInit(PyMOLGlobals * G, CPixmap * I, int width, int height);

CPixmap *PixmapNew(PyMOLGlobals * G, int width, int height);
void PixmapInitFromBitmap(PyMOLGlobals * G, CPixmap * I,
                          int width,
                          int height,
                          unsigned char *bitmap, unsigned char *rgba, int sampling);
void PixmapInitFromBytemap(PyMOLGlobals * G, CPixmap * I,
                           int width,
                           int height,
                           int pitch,
                           unsigned char *bitmap,
                           unsigned char *rgba, unsigned char *outline_rgb, int flat);
void PixmapPurge(CPixmap * I);
void PixmapFreeP(CPixmap * I);

#endif
