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

int TextInit(PyMOLGlobals *G);
int TextGetFontID(PyMOLGlobals *G,int src, int code, char *name,int size_mode, int size, int style);

void TextFree(PyMOLGlobals *G);

void TextSetPos(PyMOLGlobals *G,float *pos);

void TextSetColor(PyMOLGlobals *G,float *color);

void TextSetPosNColor(PyMOLGlobals *G,float *pos,float *color);
float *TextGetColor(PyMOLGlobals *G);
float *TextGetPos(PyMOLGlobals *G);
void TextGetColorUChar(PyMOLGlobals *G,unsigned char *red,
                       unsigned char *green, 
                       unsigned char *blue,
                       unsigned char *alpha);
char *TextRenderOpenGL(PyMOLGlobals *G,int text_id,char *st);
char *TextRenderRay(PyMOLGlobals *G,CRay *ray,int text_id,char *st);

#endif
