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

#include "PyMOLGlobals.h"

#define WordLength 64

typedef char WordType[WordLength];

typedef struct {
  WordType word;
  int value;
} WordKeyValue;

typedef struct {
  char *word;
  char **start;
  int n_word;
} CWordList;

CWordList *WordListNew(PyMOLGlobals *G,char *st);
void WordListFreeP(CWordList *I);
void WordListDump(CWordList *I,char *prefix);
int WordListMatch(PyMOLGlobals *G,CWordList *I,char *name, int ignore_case);

int WordInit(PyMOLGlobals *G);
void WordFree(PyMOLGlobals *G);

void WordSetWildcard(PyMOLGlobals *G,char wc);
int WordMatch(PyMOLGlobals *G,char *p,char *q,int ignCase); 
int WordMatchExact(PyMOLGlobals *G,char *p,char *q,int ignCase); 
void WordPrimeCommaMatch(PyMOLGlobals *G,char *p);
int WordMatchComma(PyMOLGlobals *G,char *p,char *q,int ignCase); 
int WordMatchCommaInt(PyMOLGlobals *G,char *p,int number);
int WordMatchCommaExact(PyMOLGlobals *G,char *p,char *q,int ignCase);

/* (<0) exact match, (>0) inexact match, =0 no match */

int WordCompare(PyMOLGlobals *G,char *p,char *q,int ignCase);
int WordIndex(PyMOLGlobals *G,WordType *list,char *word,int minMatch,int ignCase);
int WordKey(PyMOLGlobals *G,WordKeyValue *list,char *word,int minMatch,int ignCase,int *exact);


#endif
