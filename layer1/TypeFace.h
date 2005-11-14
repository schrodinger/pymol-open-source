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
#ifndef _H_TypeFace
#define _H_TypeFace

#include "PyMOLGlobals.h"
#include "Character.h"


int TypeInit(PyMOLGlobals *G);
void TypeFree(PyMOLGlobals *G);


typedef struct _CTypeFace CTypeFace;


CTypeFace *TypeFaceLoad(PyMOLGlobals *G,unsigned char *dat, unsigned int len);
void TypeFaceFree(CTypeFace *face);

int TypeFaceCharacterNew(CTypeFace *I,
                         CharFngrprnt *fprnt,
                         float size);

float TypeFaceGetKerning(CTypeFace *I,
                         unsigned int last, 
                         unsigned int current,
                         float size);

#endif
