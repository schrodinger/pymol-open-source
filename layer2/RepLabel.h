
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
#ifndef _H_RepLabel
#define _H_RepLabel

#include"Rep.h"
#include"CoordSet.h"

struct cLabelRelativeMode
{
  enum  {
    Default = 0,
    ScreenRelative = 1,
    ScreenPixelSpace = 2,
  };
};

Rep *RepLabelNew(CoordSet * cset, int state);

short InvalidateShaderCGOIfTextureNeedsUpdate(PyMOLGlobals *G, float font_size, int texture_font_size, int *sizeArg);

#endif
