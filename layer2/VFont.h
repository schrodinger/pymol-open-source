
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
#ifndef _H_VFont
#define _H_VFont

struct PyMOLGlobals;
class CGO;

int VFontInit(PyMOLGlobals * G);
void VFontFree(PyMOLGlobals * G);

int VFontLoad(PyMOLGlobals * G, float size, int face, int style, int can_load_new);
int VFontWriteToCGO(PyMOLGlobals * G, int font_id, CGO * cgo, const char *text,
                    float *pos, float *scale, float *matrix, float *color=nullptr);

int VFontIndent(PyMOLGlobals * G, int font_id, const char *text,
                float *pos, float *scale, float *matrix, float dir);

#endif
