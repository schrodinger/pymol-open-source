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
#ifndef _H_StrBlock
#define _H_StrBlock

/* All data structures which use StrBlock must have this header */

typedef struct { 
  int next_unused;
} StrBlock;

char  *StrBlockNew(int init_size);
int   StrBlockNewStr(char **block,const char *st,int len);
void  StrBlockFree(char *block);
void  StrBlockFreeStr(char *block,int elem);

#endif









