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

void ErrFatal(const char *where,const char *what);
void ErrPointer(const char *file,int line);
int ErrOk(const char *where,const char *what);
int ErrMessage(const char *where,const char *what);

#define ErrChkPtr(p) {if(!p) ErrPointer(__FILE__,__LINE__);}
#endif
