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
#ifndef _H_Parse
#define _H_Parse

unsigned char *ParseNextLine(unsigned char *p);
unsigned char *ParseWordCopy(unsigned char *dst,unsigned char *src,int n);
unsigned char *ParseNCopy(unsigned char *dst,unsigned char *src,int n);
unsigned char *ParseNSkip(unsigned char *p,int n);

#endif
