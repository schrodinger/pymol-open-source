
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/

#ifndef _H_FontGLUT
#define _H_FontGLUT

#include"Font.h"


/* GLUT font codes */


/* BEGIN GLUT EXCERPT.  THE FOLLOWING CODE IS:

 * Copyright (c) Mark J. Kilgard, 1994. 
 
 * This program is freely distributable without licensing fees 
 * and is provided without guarantee or warrantee expressed or 
 * implied. This program is -not- in the public domain. 

 * Modifications by Warren L. DeLano, 2004.  

*/

typedef struct {
  const int width;
  const int height;
  const float xorig;
  const float yorig;
  const float advance;
  const unsigned char *bitmap;
} FontGLUTBitmapCharRec, *FontGLUTBitmapCharPtr;

typedef struct {
  const char *name;
  const int num_chars;
  const int first;
  const FontGLUTBitmapCharRec *const *ch;
} FontGLUTBitmapFontRec, *FontGLUTBitmapFontPtr;

typedef void *GLUTbitmapFont;


/* end GLUT Excerpt */

#define cFontGLUT8x13   0
#define cFontGLUT9x15   1
#define cFontGLUTHel10  2
#define cFontGLUTHel12  3
#define cFontGLUTHel18  4

extern FontGLUTBitmapFontRec FontGLUTBitmap8By13;
extern FontGLUTBitmapFontRec FontGLUTBitmap9By15;
extern FontGLUTBitmapFontRec FontGLUTBitmapHelvetica10;
extern FontGLUTBitmapFontRec FontGLUTBitmapHelvetica12;
extern FontGLUTBitmapFontRec FontGLUTBitmapHelvetica18;

typedef struct {
  CFont Font;
  FontGLUTBitmapFontRec *glutFont;

  /* save fields */
  int swapbytes, lsbfirst, rowlength;
  int skiprows, skippixels, alignment;

} CFontGLUT;

CFont *FontGLUTNew(PyMOLGlobals * G, int font_code);
void FontGLUTFree(CFont * I);

#endif
