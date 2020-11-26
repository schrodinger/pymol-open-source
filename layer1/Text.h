
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
#ifndef _H_Text
#define _H_Text

#include"PyMOLGlobals.h"
#include"Base.h"


/* Here are the issues:

- Many different types of fonts:
   bitmap (glut)
   stroke vector
   3D-vector (extrusion) 
   pixmap (raster)
   texture

- Availability: 
   * we don't want to load fonts which consume resources before they're needed.
   * we don't want to slow down PyMOL launching any more than it already has been

- Performance:
   * really, only vector and texture fonts deliver the performance we seek
   * want zero overhead when using a pre-loaded font (no Python locks allowed)

- Raytracing: 
   * glut font cannot be raytraced, 
   * pixmap and texture mapped fonts need resolution resolution >= PNG file
   * need a geometric stand-in for the font for rendering (filled rect + alpha okay)

- Quality: 
   * vector font looks klunky
   * antialiasing essential if labels are to work with PyMOL imagery

- Flexibility: want to be able to color fonts on the fly, but also support multi-color fonts

- Threading: 
   * font loading will require the Python interpreter lock, since some fonts are defined
   * fonts are global PyMOL data

- Geometric:
   * are fonts rendered in screen space, model space, or both?
   * when we specify size, do we mean screen pixels or model distances?
   * do we support kerning? what other metrics are critical?

- Portability:
   * can we support fonts without C++?  FreeType2 is ANSI C...

*/


/* font_sources */

#define cTextSrcGeneral   0
#define cTextSrcGLUT      1
#define cTextSrcFreeType  2

int TextInit(PyMOLGlobals * G);
void TextFree(PyMOLGlobals * G);

void TextSetColorFromUColor(PyMOLGlobals * G);

void TextSetWorldPos(PyMOLGlobals * G, const float *pos);
float *TextGetWorldPos(PyMOLGlobals * G);
void TextSetLabelPushPos(PyMOLGlobals * G, const float *pos);
float *TextGetLabelPushPos(PyMOLGlobals * G);
void TextSetLabelPos(PyMOLGlobals * G, const float *pos);
float *TextGetLabelPos(PyMOLGlobals * G);
void TextSetLabelPosIsSet(PyMOLGlobals * G, unsigned char isSet);
unsigned char TextGetLabelPosIsSet(PyMOLGlobals * G);
void TextSetScreenWorldOffset(PyMOLGlobals * G, const float *pos);
float *TextGetScreenWorldOffset(PyMOLGlobals * G);
void TextSetTargetPos(PyMOLGlobals * G, const float *pos);
float *TextGetTargetPos(PyMOLGlobals * G);
float TextGetWidth(PyMOLGlobals * G);
float TextGetHeight(PyMOLGlobals * G);
void TextSetWidth(PyMOLGlobals * G, float text_width);
void TextSetHeight(PyMOLGlobals * G, float text_height);
void TextSetIndentFactorX(PyMOLGlobals * G, float factor);
void TextSetIndentFactorY(PyMOLGlobals * G, float factor);
float *TextGetIndentFactor(PyMOLGlobals * G);
void TextSetPos(PyMOLGlobals * G, const float *pos);
void TextSetColor(PyMOLGlobals * G, const float *color);
void TextSetColor3f(PyMOLGlobals * G, float red, float green, float blue);
void TextGetOutlineColor(PyMOLGlobals * G,
                         unsigned char *red,
                         unsigned char *green, unsigned char *blue, unsigned char *alpha);
void TextSetOutlineColor(PyMOLGlobals * G, int color);
void TextSetPosNColor(PyMOLGlobals * G, const float *pos, const float *color);
float *TextGetColor(PyMOLGlobals * G);
float *TextGetPos(PyMOLGlobals * G);
void TextGetColorUChar(PyMOLGlobals * G, unsigned char *red,
                       unsigned char *green, unsigned char *blue, unsigned char *alpha);
unsigned char *TextGetColorUChar4uv(PyMOLGlobals * G);
const char *TextRenderOpenGL(PyMOLGlobals * G, const RenderInfo * info, int text_id, const char *st,
                       float size, const float *rpos, short needSize, short relativeMode, short shouldRender, CGO *shaderCGO);
const char *TextRenderRay(PyMOLGlobals * G, CRay * ray, int text_id, const char *st, float size,
                    const float *rpos, short needSize, short relativeMode);

void TextDrawStrAt(PyMOLGlobals * G, const char *st, int x, int y ORTHOCGOARG);
void TextDrawStr(PyMOLGlobals * G, const char *st ORTHOCGOARG);
void TextIndent(PyMOLGlobals * G, float x, float y);
void TextAdvance(PyMOLGlobals * G, float advance);
void TextSetPos2i(PyMOLGlobals * G, int x, int y);
void TextDrawChar(PyMOLGlobals * G, char ch ORTHOCGOARG);
void TextDrawSubStrFast(PyMOLGlobals * G, const char *c, int x, int y, int start, int n ORTHOCGOARG);
void TextDrawCharRepeat(PyMOLGlobals * G, char c, int x, int y, int start, int n ORTHOCGOARG);

void TextSetLabelBkgrdInfo(PyMOLGlobals * G, float label_spacing, float label_just, const float *buff);

float TextGetSpacing(PyMOLGlobals * G);
float TextGetJustification(PyMOLGlobals * G);
float *TextGetLabelBuffer(PyMOLGlobals * G);

void TextSetIsPicking(PyMOLGlobals * G, bool IsPicking);
bool TextGetIsPicking(PyMOLGlobals * G);

bool TextStartsWithColorCode(const char *);
bool TextSetColorFromCode(PyMOLGlobals *, const char *, const float *);

#endif
