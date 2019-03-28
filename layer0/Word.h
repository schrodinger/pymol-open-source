

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
#include "Lex.h"

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

CWordMatcher *WordMatcherNew(PyMOLGlobals * G, const char *st, CWordMatchOptions * option,
                             int force);
int WordMatcherMatchAlpha(CWordMatcher * I, const char *text);
int WordMatcherMatchMixed(CWordMatcher * I, const char *text, int value);
int WordMatcherMatchInteger(CWordMatcher * I, int value);
void WordMatcherFree(CWordMatcher * I);
int WordMatchNoWild(PyMOLGlobals * G, const char *p, const char *q, int ignCase);

CWordList *WordListNew(PyMOLGlobals * G, const char *st);
void WordListFreeP(CWordList * I);
void WordListDump(CWordList * I, const char *prefix);
int WordListIterate(PyMOLGlobals * G, CWordList * I, const char **ptr, int *hidden);
int WordListMatch(PyMOLGlobals * G, CWordList * I, const char *name, int ignore_case);

int WordInit(PyMOLGlobals * G);
void WordFree(PyMOLGlobals * G);

int WordMatch(PyMOLGlobals * G, const char *p, const char *q, int ignCase);
int WordMatchExact(PyMOLGlobals * G, const char *p, const char *q, int ignCase);
void WordPrimeCommaMatch(PyMOLGlobals * G, char *p);
int WordMatchComma(PyMOLGlobals * G, const char *p, const char *q, int ignCase);
int WordMatchCommaInt(PyMOLGlobals * G, const char *p, int number);
int WordMatchCommaExact(PyMOLGlobals * G, const char *p, const char *q, int ignCase);


/* (<0) exact match, (>0) inexact match, =0 no match */

int WordIndex(PyMOLGlobals * G, WordType * list, const char *word, int minMatch, int ignCase);
int WordKey(PyMOLGlobals * G, WordKeyValue * list, const char *word, int minMatch, int ignCase,
            int *exact);

int WordCompare(PyMOLGlobals * G, const char *p, const char *q, int ignCase);

inline int WordCompare(PyMOLGlobals * G, const lexidx_t& s1, const lexidx_t& s2, int ignCase) {
  if (s1 == s2)
    return 0;
  return WordCompare(G, LexStr(G, s1), LexStr(G, s2), ignCase);
}

inline int WordMatch(PyMOLGlobals * G, const lexidx_t& s1, const lexidx_t& s2, int ignCase) {
  if (s1 == s2)
    return -1; // negative = perfect match
  return WordMatch(G, LexStr(G, s1), LexStr(G, s2), ignCase);
}

inline int WordMatchNoWild(PyMOLGlobals * G, const lexidx_t& s1, const lexidx_t& s2, int ignCase) {
  if (s1 == s2)
    return -1; // negative = perfect match
  return WordMatchNoWild(G, LexStr(G, s1), LexStr(G, s2), ignCase);
}

inline int WordMatchExact(PyMOLGlobals * G, const lexidx_t& s1, const lexidx_t& s2, int ignCase) {
  if (s1 == s2)
    return 1; // non-zero = perfect match
  if (!ignCase)
    return 0; // 0 = no match
  return WordMatchExact(G, LexStr(G, s1), LexStr(G, s2), ignCase);
}

inline int WordMatchExact(PyMOLGlobals * G, char c1, char c2, int ignCase) {
  if (c1 == c2)
    return 1; // non-zero = perfect match
  if (!ignCase)
    return 0; // 0 = no match
  return c1 && c2 && toupper(c1) == toupper(c2);
}

#endif
