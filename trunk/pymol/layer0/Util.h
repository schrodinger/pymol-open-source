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
#ifndef _H_Util
#define _H_Util

void UtilZeroMem(void *ptr,unsigned int howMuch);
void *UtilArrayMalloc(unsigned int *dim,int ndim,unsigned int atom_size);
char *UtilConcat(char *where,char *what);
void UtilCleanStr(char *s);
double UtilGetSeconds(void);
void UtilInit(void);

#endif
