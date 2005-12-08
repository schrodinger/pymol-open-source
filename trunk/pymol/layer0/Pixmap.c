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
                             unsigned char *rgba,int sampling)
{
  if(I) {
 
    int x,y,bit_cnt;
    unsigned char cur=0;
    unsigned char *src;
    unsigned char *dst;
    register unsigned char red,blue,green,alpha;
    PixmapInit(G,I,width*sampling,height*sampling);
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
    if(sampling>1) {
      unsigned int *p, *pp, *q, *row;
      int row_cnt,col_cnt, width_sampling = width*sampling;

      p = (unsigned int*)(I->buffer + 4*width*height);
      q = (unsigned int*)(I->buffer + 4*width*height*sampling*sampling);
      while(p>(unsigned int*)I->buffer) {
        row_cnt = sampling - 1;
        row = q;
        for(x=0;x<width;x++) { /* first row */
          col_cnt = sampling;
          p--;
          while(col_cnt--) {
            *(--q) = *p;
          }
        }
        if(row_cnt) {
          while(row_cnt--) { /* remaining rows */
            pp = row;
            for(x=0;x<width_sampling;x++) {
              *(--q) = *(--pp);
            }
          }
        }
      }
    }
  }
}
void PixmapInitFromBytemap(PyMOLGlobals *G,CPixmap *I,
                           int width, 
                           int height,
                           int pitch,
                           unsigned char *bytemap,
                           unsigned char *rgba,
                           unsigned char *outline_rgb,
                           int flat
                           )
{
  if(I) {
    
    int x,y;
    unsigned char *src,*sa,alp;
    unsigned char *dst;
    register unsigned char red,blue,green,alpha,no_alpha;
    register unsigned char ored=0,oblue=0,ogreen=0;
    if(!outline_rgb[3])
      outline_rgb = NULL;
    else {
      ored = outline_rgb[0];
      oblue = outline_rgb[1];
      ogreen = outline_rgb[2];
    }
    PixmapInit(G,I,width,height);
    red = rgba[0];
    green = rgba[1];
    blue = rgba[2];
    alpha = rgba[3];
    UtilZeroMem(I->buffer,4*width*height);
    src = bytemap;
    dst = I->buffer;
    no_alpha = flat;
    for(y=0;y<height;y++) {
      sa = src;
      if(no_alpha) {
        for(x=0;x<width;x++) {
          alp = *(sa++);
          if(alp) {
            *(dst++)=red;
            *(dst++)=green;
            *(dst++)=blue;
            *(dst++)=0xFF;
          } else {
            *(dst++)=0;
            *(dst++)=0;
            *(dst++)=0;
            *(dst++)=0;
          }
        }
      } else {
        for(x=0;x<width;x++) {
          if(outline_rgb) {
            unsigned char amax = 0,amin;
            if(y>0) {
              alp = 255 - *(sa - pitch);
            } else {
              alp = 255;
            }
            if(amax<alp) amax = alp;
            if(y<(height-1)) {
              alp = 255 - *(sa + pitch);
            } else {
              alp = 255;
            }
            if(amax<alp) amax = alp;
            if(x>0) {
              alp = 255 - *(sa - 1);
            } else {
              alp = 255;
            }
            if(amax<alp) amax = alp;
            if(x<(width-1)) {
              alp = 255 - *(sa + 1);  
            } else {
              alp = 255;
            }
            if(amax<alp) amax = alp;
            amin = 255 - amax;
            alp = *(sa++);
            if(alp) {
              *(dst++)=(red * amin + ored * amax)/255;
              *(dst++)=(green * amin + ogreen * amax)/255;
              *(dst++)=(blue * amin+ oblue * amax)/255;
              *(dst++)=(alpha * alp)/255;
            } else {
              *(dst++)=0;
              *(dst++)=0;
              *(dst++)=0;
              *(dst++)=0;
            }
          } else {
            alp = *(sa++);
            if(alp) {
              *(dst++)=red;
              *(dst++)=green;
              *(dst++)=blue;
              *(dst++)=(alpha * alp)>>8;
            } else {
              *(dst++)=0;
              *(dst++)=0;
              *(dst++)=0;
              *(dst++)=0;
            }
          }
          
        }
      }
      src+=pitch;
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

