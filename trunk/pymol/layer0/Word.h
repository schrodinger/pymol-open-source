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
#ifndef _H_Word
#define _H_Word

#define WordLength 64

typedef char WordType[WordLength];

int WordMatch(char *p,char *q,int ignCase); 
/* (<0) exact match, (>0) inexact match, =0 no match */

int WordCompare(char *p,char *q,int ignCase);
unsigned int WordChoose(WordType *list, char *word,int minMatch,int ignCase);
int WordIndex(WordType *list,char *word,int minMatch,int ignCase);

#endif
