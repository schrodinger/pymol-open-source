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
#ifndef _H_Texture
#define _H_Texture

#include "PyMOLGlobals.h"

int TextureInit(PyMOLGlobals *G);
void TextureFree(PyMOLGlobals *G);
int TextureGetFromChar(PyMOLGlobals *G, int char_id,float *extent);


#endif
