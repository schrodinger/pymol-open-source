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

#ifndef _H_Font
#define _H_Font

#define cFontNULL        0
#define cFontGLUT        1
#define cFontBitmap      2
#define cFontVector      3
#define cFontExtrude     4
#define cFontPixmap      5
#define cFontTexture     6

#define cFontModePixels    0
#define cFontModeSpacial   1

struct CFont;
struct CRay;

typedef char *FontRenderOpenGLFn(struct CFont *,char *);
typedef char *FontRenderRayFn(struct CRay *ray,struct CFont *,char *);

typedef struct CFont {
  int TextID;
  void (*fFree)(struct CFont *);
  FontRenderOpenGLFn *fRenderOpenGL;
  FontRenderRayFn *fRenderRay;
} CFont;

int FontInit(CFont *I);

void FontPurge(CFont *I);

#endif

