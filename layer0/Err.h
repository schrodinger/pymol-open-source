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
#ifndef _H_Err
#define _H_Err

#include "PyMOLGlobals.h"

void ErrFatal(PyMOLGlobals *G,const char *where,const char *what);
void ErrPointer(PyMOLGlobals *G,const char *file,int line);
int ErrMessage(PyMOLGlobals *G,const char *where,const char *what);

#define ErrChkPtr(G,p) {if(!p) ErrPointer(G,__FILE__,__LINE__);}
#endif
