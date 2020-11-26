
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

extern FontGLUTBitmapFontRec FontGLUTBitmap8By13;
extern FontGLUTBitmapFontRec FontGLUTBitmap9By15;
extern FontGLUTBitmapFontRec FontGLUTBitmapHelvetica10;
extern FontGLUTBitmapFontRec FontGLUTBitmapHelvetica12;
extern FontGLUTBitmapFontRec FontGLUTBitmapHelvetica18;

struct CFontGLUT : public CFont {
  const FontGLUTBitmapFontRec *glutFont;

  /* save fields */
  int swapbytes, lsbfirst, rowlength;
  int skiprows, skippixels, alignment;
  const char* RenderOpenGL(const RenderInfo* info, const char* text,
      float size, const float* rpos, bool needSize, short relativeMode,
      bool shouldRender, CGO* shaderCGO) override;
  const char* RenderOpenGLFlat(const RenderInfo* info, const char* text,
      float size, const float* rpos, bool needSize, short relativeMode,
      bool shouldRender, CGO* shaderCGO) override;
  const char* RenderRay(CRay* ray, const char* text, float size,
      const float* rpos, bool needSize, short relativeMode) override;
  CFontGLUT(PyMOLGlobals* G, const FontGLUTBitmapFontRec*);
};

#endif
