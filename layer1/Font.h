
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

#include"PyMOLGlobals.h"
#include"Base.h"

#define cFontNULL        0
#define cFontGLUT        1
#define cFontBitmap      2
#define cFontVector      3
#define cFontExtrude     4
#define cFontPixmap      5
#define cFontTexture     6

#define cFontModePixels    0
#define cFontModeSpacial   1

#include "os_gl.h"

#define SHADERCGOARG , CGO *shaderCGO
#define SHADERCGOARGVAR ,shaderCGO
#define SHADERCGOARGV shaderCGO

struct CFont {
  PyMOLGlobals* G = nullptr;
  int TextID = 0;
  virtual const char* RenderOpenGL(const RenderInfo* info, const char* text,
      float size, const float* rpos, bool needSize, short relativeMode,
      bool shouldRender, CGO* shaderCGO)
  {
    return nullptr;
  }
  virtual const char* RenderOpenGLFlat(const RenderInfo* info, const char* text,
      float size, const float* rpos, bool needSize, short relativeMode,
      bool shouldRender, CGO* shaderCGO)
  {
    return nullptr;
  }
  virtual const char* RenderRay(CRay* ray, const char* text, float size,
      const float* rpos, bool needSize, short relativeMode)
  {
    return nullptr;
  }
  CFont(PyMOLGlobals* G)
      : G(G){};
  virtual ~CFont() = 0;
};

#endif
