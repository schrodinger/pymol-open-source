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
#ifndef _H_Pixmap
#define _H_Pixmap


/* for the sake of simplicity, all pixmaps are 32-bit RGBA */

typedef struct {
  int height,width;
  unsigned char *buffer;
} CPixmap;

void PixmapInit(CPixmap *I,int width,int height);

CPixmap *PixmapNew(int width,int height);
void PixmapInitFromBitmap(CPixmap *I,int width, int height,
                             unsigned char *bitmap,
                             unsigned char *rgba);
void PixmapPurge(CPixmap *I);
void PixmapFreeP(CPixmap *I);

#endif

