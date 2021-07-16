

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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

struct PyMOLGlobals;

void ErrFatal(const PyMOLGlobals * G, const char *where, const char *what);
void ErrPointer(const PyMOLGlobals * G, const char *file, int line);
int ErrMessage(PyMOLGlobals * G, const char *where, const char *what);

#define ErrChkPtr(G,p) {if(!p) ErrPointer(G,__FILE__,__LINE__);}

#define CHECKOK(ok, var) ok &= var ? true : false;

#define ok_raise(x) goto ok_except ## x
#define ok_assert(x, c) {if(!(c)) ok_raise(x);}

#endif
