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

#include "Pixmap.h"
#include "OOMac.h"
#include "Util.h"

void PixmapInit(PyMOLGlobals *G,CPixmap *I,int width,int height)
{
  UtilZeroMem(I,sizeof(CPixmap));
  I->G=G;
  I->height = height;
  I->width = width;
  if((height>=0)&&(width>=0)) {
    I->buffer = Alloc(unsigned char,4*height*width);
  }
}

CPixmap *PixmapNew(PyMOLGlobals *G,int width,int height)
{
  OOAlloc(G,CPixmap);
  PixmapInit(G,I,width,height);
  return I;
}

void PixmapInitFromBitmap(PyMOLGlobals *G,CPixmap *I,int width, int height,
                             unsigned char *bitmap,
                             unsigned char *rgba)
{
  if(I) {
 
    int x,y,bit_cnt;
    unsigned char cur=0;
    unsigned char *src;
    unsigned char *dst;
    register unsigned char red,blue,green,alpha;
    PixmapInit(G,I,width,height);
    red = rgba[0];
    green = rgba[1];
    blue = rgba[2];
    alpha = rgba[3];
    UtilZeroMem(I->buffer,4*width*height);
    src = bitmap;
    dst = I->buffer;
    for(y=0;y<height;y++) {
      bit_cnt = 7;
      for(x=0;x<width;x++) {
        bit_cnt++;
        if(bit_cnt>7) {
          cur = *(src++);
          bit_cnt = 0;
        }
        if(cur&0x80) {
          *(dst++)=red;
          *(dst++)=green;
          *(dst++)=blue;
          *(dst++)=alpha;
        } else {
          *(dst++)=0;
          *(dst++)=0;
          *(dst++)=0;
          *(dst++)=0;
        }
        cur <<= 1;
      }
    }
  }
}

void PixmapPurge(CPixmap *I)
{
  if(I) {
    FreeP(I->buffer);
  }
}

void PixmapFreeP(CPixmap *I)
{
  PixmapPurge(I);
  OOFreeP(I);
}

