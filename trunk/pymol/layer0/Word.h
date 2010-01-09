

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
#ifndef _H_Word
#define _H_Word

#include "PyMOLGlobals.h"

#define WordLength 256

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

#define cWordMatchOptionNoRanges 0
#define cWordMatchOptionNumericRanges 1
#define cWordMatchOptionAlphaRanges 2

typedef struct {
  int range_mode;               /* 0 = none, 1 = numeric, 2 = alpha */
  int lists;
  int ignore_case;
  int allow_hyphen;
  int allow_plus;
  int space_lists;
  char wildcard;
} CWordMatchOptions;

typedef struct _CWordMatcher CWordMatcher;

void WordMatchOptionsConfigInteger(CWordMatchOptions * I);
void WordMatchOptionsConfigAlpha(CWordMatchOptions * I, char wildcard, int ignore_case);
void WordMatchOptionsConfigAlphaList(CWordMatchOptions * I, char wildcard,
                                     int ignore_case);
void WordMatchOptionsConfigMixed(CWordMatchOptions * I, char wildcard, int ignore_case);
void WordMatchOptionsConfigNameList(CWordMatchOptions * I, char wildcard,
                                    int ignore_case);

CWordMatcher *WordMatcherNew(PyMOLGlobals * G, char *st, CWordMatchOptions * option,
                             int force);
int WordMatcherMatchAlpha(CWordMatcher * I, char *text);
int WordMatcherMatchMixed(CWordMatcher * I, char *text, int value);
int WordMatcherMatchInteger(CWordMatcher * I, int value);
void WordMatcherFree(CWordMatcher * I);
int WordMatchNoWild(PyMOLGlobals * G, char *p, char *q, int ignCase);

CWordList *WordListNew(PyMOLGlobals * G, char *st);
void WordListFreeP(CWordList * I);
void WordListDump(CWordList * I, char *prefix);
int WordListIterate(PyMOLGlobals * G, CWordList * I, char **ptr, int *hidden);
int WordListMatch(PyMOLGlobals * G, CWordList * I, char *name, int ignore_case);

int WordInit(PyMOLGlobals * G);
void WordFree(PyMOLGlobals * G);

void WordSetWildcard(PyMOLGlobals * G, char wc);
int WordMatch(PyMOLGlobals * G, char *p, char *q, int ignCase);
int WordMatchExact(PyMOLGlobals * G, char *p, char *q, int ignCase);
void WordPrimeCommaMatch(PyMOLGlobals * G, char *p);
int WordMatchComma(PyMOLGlobals * G, char *p, char *q, int ignCase);
int WordMatchCommaInt(PyMOLGlobals * G, char *p, int number);
int WordMatchCommaExact(PyMOLGlobals * G, char *p, char *q, int ignCase);


/* (<0) exact match, (>0) inexact match, =0 no match */

int WordIndex(PyMOLGlobals * G, WordType * list, char *word, int minMatch, int ignCase);
int WordKey(PyMOLGlobals * G, WordKeyValue * list, char *word, int minMatch, int ignCase,
            int *exact);

#ifdef _PYMOL_INLINE
__inline__ static int WordCompare(PyMOLGlobals * G, char *p, char *q, int ignCase)


/* all things equal, shorter is smaller */
{
  int result = 0;
  register char cp, cq, tlp, tlq;
  if(ignCase) {
    while((cp = *p) && (cq = *q)) {
      p++;
      q++;
      if(cp != cq) {
        (tlp = tolower(cp));
        (tlq = tolower(cq));
        if(tlp < tlq)
          return -1;
        else if(tlp > tlq) {
          return 1;
        }
      }
    }
  } else {
    while((cp = *p) && (cq = *q)) {
      p++;
      q++;
      if(cp != cq) {
        if(cp < cq) {
          return -1;
        } else if(cp > cq) {
          return 1;
        }
      }
    }
  }
  if((!result) && (!*p) && (*q))
    return -1;
  else if((!result) && (*p) && (!*q))
    return 1;
  return 0;
}
#else
int WordCompare(PyMOLGlobals * G, char *p, char *q, int ignCase);

#endif

#endif
