

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
#include"os_python.h"

#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"OOMac.h"
#include"Word.h"
#include"Parse.h"
#include"PyMOLObject.h"
#include"MemoryDebug.h"

struct _CWord {
  int no_state_at_present;
};

typedef struct {
  int match_mode;
  int continued;
  int literal1, literal2;       /* offsets into charVLA */
  int numeric1, numeric2;
  int has1, has2;
} MatchNode;

struct _CWordMatcher {
  PyMOLGlobals *G;
  MatchNode *node;
  int n_node;
  char *charVLA;
  int n_char;
  int ignore_case;
};

#define cMatchLiteral  0
#define cMatchNumericRange  cWordMatchOptionNumericRanges
#define cMatchAlphaRange  cWordMatchOptionAlphaRanges
#define cMatchWildcard 3
#if 0
static void WordMatcherDump(CWordMatcher * I)
{
  MatchNode *cur_node = I->node;
  register int n_node = I->n_node;
  while(n_node--) {
    printf("match_mode %d\n", cur_node->match_mode);
    printf("literal1 '%s' %d literal2 '%s' %d\n",
           I->charVLA + cur_node->literal1,
           cur_node->has1, I->charVLA + cur_node->literal2, cur_node->has2);
    cur_node++;
  }
}
#endif

#ifndef _PYMOL_INLINE
int WordCompare(PyMOLGlobals * G, char *p, char *q, int ignCase)


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
#endif

void WordMatchOptionsConfigInteger(CWordMatchOptions * I)
{
  I->range_mode = cWordMatchOptionNumericRanges;
  I->lists = true;
  I->ignore_case = true;
  I->wildcard = 0;              /* no wildcard for numbers */
  I->allow_hyphen = true;
  I->allow_plus = true;
  I->space_lists = false;
}

void WordMatchOptionsConfigAlpha(CWordMatchOptions * I, char wildcard, int ignore_case)
{
  I->range_mode = cWordMatchOptionAlphaRanges;
  I->lists = true;
  I->ignore_case = ignore_case;
  I->wildcard = wildcard;
  I->allow_hyphen = false;
  I->allow_plus = false;
  I->space_lists = false;
}

void WordMatchOptionsConfigAlphaList(CWordMatchOptions * I, char wildcard,
                                     int ignore_case)
{                               /* here we expect '+' to be used in lists */
  I->range_mode = cWordMatchOptionAlphaRanges;
  I->lists = true;
  I->ignore_case = ignore_case;
  I->wildcard = wildcard;
  I->allow_hyphen = false;
  I->allow_plus = true;
  I->space_lists = false;
}

void WordMatchOptionsConfigMixed(CWordMatchOptions * I, char wildcard, int ignore_case)
{
  I->range_mode = cWordMatchOptionNumericRanges;
  I->lists = true;
  I->ignore_case = ignore_case;
  I->wildcard = wildcard;
  I->allow_hyphen = true;
  I->allow_plus = true;
  I->space_lists = false;
}

void WordMatchOptionsConfigNameList(CWordMatchOptions * I, char wildcard, int ignore_case)
{
  I->range_mode = cWordMatchOptionAlphaRanges;
  I->lists = true;
  I->ignore_case = ignore_case;
  I->wildcard = wildcard;
  I->allow_hyphen = false;
  I->allow_plus = false;
  I->space_lists = true;
}

CWordMatcher *WordMatcherNew(PyMOLGlobals * G, char *st, CWordMatchOptions * option,
                             int force)
{

  CWordMatcher *result = NULL;
  int needed = force;
  char wildcard = option->wildcard;

  if(wildcard == 32)
    wildcard = 0;               /* space as wildcard means no wildcard */

  if(!st)
    return NULL;
  {                             /* first determine if we need to incur the overhead of the matcher */
    int escape = false;
    char *p = st;
    while((*p) && (!needed)) {
      if(!escape) {
        switch (*p) {
        case '\\':
          escape = true;
          needed = true;
          break;
        case '+':
          if((option->lists) && (option->allow_plus))
            needed = true;
          break;
        case ',':              /* list operators */
          if(option->lists)
            needed = true;
          break;
        case '-':              /* range operators */
          if(option->allow_hyphen)
            needed = true;
          break;
        case ':':
          if(option->range_mode)
            needed = true;
          break;
        case ' ':
          if(option->space_lists)
            needed = true;
          break;
        default:
          if(*p == wildcard)
            needed = true;
          break;
        }
      } else
        escape = false;
      p++;
    }
  }

  if(needed) {                  /* if so, then convert the expression into a match tree */
    register int n_char = 0;
    register int n_node = 0;
    OOCalloc(G, CWordMatcher);
    I->charVLA = VLACalloc(char, 10);   /* auto_zeroing... */
    I->node = VLACalloc(MatchNode, 10);
    I->ignore_case = option->ignore_case;
    I->G = G;
    /* build up the matcher structure... */
    {
      char *p = st, c, *q;
      int escape = false;
      int token_active = false;
      int node_active = false;
      int char_handled = false;
      int cur_node = 0;
      int expectation = 1;

      while(1) {
        c = *p;
        char_handled = false;
        if(!escape) {
          switch (c) {
          case '\\':
            escape = true;
            char_handled = true;
            break;
          case 0:
            if(option->lists) {
              char_handled = true;
              node_active = false;
              token_active = false;
            }
            break;
          case '+':            /* list operator */
            if(option->lists && option->allow_plus) {
              if(n_node < expectation) {        /* create empty node */
                VLACheck(I->node, MatchNode, n_node);
                n_node++;
              } else {
                expectation = n_node + 1;
              }
              char_handled = true;
              node_active = false;
              token_active = false;
            }
            break;
          case ',':            /* list operator */
            if(option->lists) {
              if(n_node < expectation) {        /* create empty node */
                VLACheck(I->node, MatchNode, n_node);
                n_node++;
              } else {
                expectation = n_node + 1;
              }
              char_handled = true;
              node_active = false;
              token_active = false;
            }
            break;
          case ' ':            /* space list */
            if(option->space_lists) {
              if(n_node < expectation) {        /* create empty node */
                VLACheck(I->node, MatchNode, n_node);
                n_node++;
              } else {
                expectation = n_node + 1;
              }
              char_handled = true;
              node_active = false;
              token_active = false;
            }
            break;
          case '-':            /* range operators */
            if(option->allow_hyphen && option->range_mode) {
              if(!node_active) {
                cur_node = n_node;
                VLACheck(I->node, MatchNode, cur_node);
                node_active = true;
                n_node++;
              }
              I->node[cur_node].match_mode = option->range_mode;
              token_active = false;
              char_handled = true;
            }
            break;
          case ':':
            if(option->range_mode) {
              if(!node_active) {
                cur_node = n_node;
                VLACheck(I->node, MatchNode, cur_node);
                node_active = true;
                n_node++;
              }
              I->node[cur_node].match_mode = option->range_mode;
              token_active = false;
              char_handled = true;
            }
            break;
          default:
            if(c == wildcard) {
              if(node_active) {
                I->node[cur_node].continued = true;
              }
              VLACheck(I->node, MatchNode, n_node);
              cur_node = n_node;
              I->node[cur_node].match_mode = cMatchWildcard;
              n_node++;
              node_active = true;
              token_active = false;
              char_handled = true;
            }
            break;
          }
        } else
          escape = false;
        if(!char_handled) {
          if(!token_active) {
            n_char++;
            VLACheck(I->charVLA, char, n_char);
            token_active = true;

            if((!node_active) || (I->node[cur_node].match_mode == cMatchWildcard)) {
              if(node_active)   /* must be extending after a wildcard */
                I->node[cur_node].continued = true;
              else
                node_active = true;
              VLACheck(I->node, MatchNode, n_node);
              cur_node = n_node;
              I->node[cur_node].literal1 = n_char;      /* the first literal */
              n_node++;
            } else {
              I->node[cur_node].literal2 = n_char;      /* must be the second literal */
            }
          }

          /* copy character into auto-terminated string */
          VLACheck(I->charVLA, char, n_char + 1);
          q = I->charVLA + n_char;
          (*q++) = c;
          n_char++;
        }
        if(c)
          p++;
        else
          break;
      }
      if(n_node < expectation) {        /* create empty node */
        VLACheck(I->node, MatchNode, n_node);
        n_node++;
      }

    }

    {
      int a;
      int tmp;
      MatchNode *node = I->node;
      for(a = 0; a < n_node; a++) {
        switch (node->match_mode) {
        case cMatchLiteral:
          if(option->range_mode == cWordMatchOptionNumericRanges) {
            if(node->literal1) {
              if(sscanf(I->charVLA + node->literal1, "%d", &tmp) == 1) {
                node->numeric1 = tmp;
                node->has1 = true;
              }
            }
          }
          break;
        case cMatchAlphaRange:
          if(node->literal1)
            node->has1 = true;
          if(node->literal2)
            node->has2 = true;

          break;
        case cMatchNumericRange:
          if(node->literal1) {
            if(sscanf(I->charVLA + node->literal1, "%d", &tmp) == 1) {
              node->numeric1 = tmp;
              node->has1 = true;
            }
          }
          if(node->literal2) {
            if(sscanf(I->charVLA + node->literal2, "%d", &tmp) == 1) {
              node->numeric2 = tmp;
              node->has2 = true;
            }
          }
          break;
        }
        node++;
      }
    }
    I->n_char = n_char;
    I->n_node = n_node;
    /*        WordMatcherDump(I); */
    result = I;
  }
  return result;
}

static int recursive_match(CWordMatcher * I, MatchNode * cur_node, char *text,
                           int *value_ptr)
{
  int ignore_case = I->ignore_case;
  switch (cur_node->match_mode) {
  case cMatchLiteral:
    {
      char *q = I->charVLA + cur_node->literal1;
      char *p = text;
      while((*p) && (*q)) {
        if(*p != *q) {
          if(!ignore_case)
            return false;
          else if(tolower(*p) != tolower(*q))
            return false;
        }
        p++;
        q++;
      }
      if((!(*q)) && (!(*p)))
        return true;

      if((*p) && (!*q) && cur_node->continued)
        return recursive_match(I, cur_node + 1, p, value_ptr);
      if((*p) != (*q))
        return false;
    }
    break;
  case cMatchWildcard:
    {
      char *p;
      p = text;
      if(!cur_node->continued)
        return true;
      else {
        while(*p) {
          if(recursive_match(I, cur_node + 1, p, value_ptr))
            return 1;
          p++;
        }
      }
    }
    break;
  case cMatchAlphaRange:
    {
      char *l1 = I->charVLA + cur_node->literal1;
      char *l2 = I->charVLA + cur_node->literal2;
      if(((!cur_node->has1) ||
          (WordCompare(I->G, l1, text, ignore_case) <= 0)) &&
         ((!cur_node->has2) || (WordCompare(I->G, l2, text, ignore_case) >= 0)))
        return true;
      else
        return false;
    }
    break;
  case cMatchNumericRange:
    if(value_ptr) {
      int value = *value_ptr;

      if(((!cur_node->has1) ||
          (cur_node->numeric1 <= value)) &&
         ((!cur_node->has2) || (cur_node->numeric2 >= value)))
        return true;
    } else {
      int value;
      if(sscanf(text, "%d", &value) == 1)
        if(((!cur_node->has1) ||
            (cur_node->numeric1 <= value)) &&
           ((!cur_node->has2) || (cur_node->numeric2 >= value)))
          return true;
    }
    break;
  }
  return false;
}

int WordMatcherMatchAlpha(CWordMatcher * I, char *text)
{
  register MatchNode *cur_node = I->node;
  register int n_node = I->n_node;

  while((n_node--) > 0) {
    if(recursive_match(I, cur_node, text, NULL))
      return true;
    else {
      while(cur_node->continued) {
        cur_node++;
        n_node--;
      }
      cur_node++;
    }
  }
  return false;
}

int WordMatcherMatchMixed(CWordMatcher * I, char *text, int value)
{
  register MatchNode *cur_node = I->node;
  register int n_node = I->n_node;

  while((n_node--) > 0) {
    if(recursive_match(I, cur_node, text, &value))
      return true;
    else {
      while(cur_node->continued) {
        cur_node++;
        n_node--;
      }
      cur_node++;
    }
  }
  return false;
}

static int integer_match(CWordMatcher * I, MatchNode * cur_node, int value)
{
  switch (cur_node->match_mode) {
  case cMatchLiteral:
    if((cur_node->has1) && (cur_node->numeric1 == value))
      return true;
    break;
  case cMatchNumericRange:
    if(((!cur_node->has1) ||
        (cur_node->numeric1 <= value)) &&
       ((!cur_node->has2) || (cur_node->numeric2 >= value)))
      return true;
    break;
  }
  return false;
}

int WordMatcherMatchInteger(CWordMatcher * I, int value)
{
  register MatchNode *cur_node = I->node;
  register int n_node = I->n_node;

  while((n_node--) > 0) {
    if(integer_match(I, cur_node, value))
      return true;
    else {
      while(cur_node->continued) {
        cur_node++;
        n_node--;
      }
      cur_node++;
    }
  }
  return false;
}

void WordMatcherFree(CWordMatcher * I)
{
  if(I) {
    VLAFreeP(I->node);
    VLAFreeP(I->charVLA);
  }
  OOFreeP(I);
}

CWordList *WordListNew(PyMOLGlobals * G, char *st)
{
  int n_word = 0;
  char *p;
  int len = 0;
  OOCalloc(G, CWordList);

  if(I) {
    p = st;
    /* first, count how many words we have */
    while(*p) {
      if(*p > 32) {
        n_word++;
        while((*p) > 32) {
          len++;
          p++;
        }
        len++;
      } else
        p++;
    }
    /* allocate the storage we'll need to hold the words */
    {
      I->word = Alloc(char, len);
      I->start = Alloc(char *, n_word);

      /* and copy the words */

      if(I->word && I->start) {
        char *q = I->word;
        char **q_ptr = I->start;
        p = st;
        while(*p) {
          if(*p > 32) {
            *(q_ptr++) = q;
            while((*p) > 32) {
              *(q++) = *(p++);
            }
            *(q++) = 0;
            len++;
          } else
            p++;
        }
        I->n_word = n_word;
      }
    }
  }
  return I;
}

void WordListFreeP(CWordList * I)
{
  if(I) {
    FreeP(I->word);
    FreeP(I->start);
    FreeP(I);
  }
}

void WordListDump(CWordList * I, char *prefix)
{
  if(I) {
    int a;
    printf(" %s: n_word %d\n", prefix, I->n_word);
    for(a = 0; a < I->n_word; a++) {
      printf(" %s: word %d=[%s]\n", prefix, a, I->start[a]);
    }
  }
}

int WordListIterate(PyMOLGlobals * G, CWordList * I, char **ptr, int *hidden)
{
  int result = true;
  if(*hidden >= 0) {
    if(*hidden < I->n_word) {
      (*ptr) = I->start[(*hidden)++];
    } else {
      result = false;
    }
  }
  return result;
}

int WordListMatch(PyMOLGlobals * G, CWordList * I, char *name, int ignore_case)
{
  int result = -1;
  if(I) {
    int a;
    for(a = 0; a < I->n_word; a++) {
      if(WordMatch(G, I->start[a], name, ignore_case)) {
        result = a;
        break;
      }
    }
  }
  return result;
}

int WordInit(PyMOLGlobals * G)
{
  register CWord *I = NULL;

  I = (G->Word = Calloc(CWord, 1));
  if(I) {
    return 1;
  } else
    return 0;

}

void WordFree(PyMOLGlobals * G)
{
  FreeP(G->Word);
}

void WordPrimeCommaMatch(PyMOLGlobals * G, char *p)
{                               /* replace '+' with ',' */
  while(*p) {                   /* this should not be done here... */
    if(*p == '+')
      if(!((*(p + 1) == 0) || (*(p + 1) == ',') || (*(p + 1) == '+')))
        *p = ',';
    p++;
  }
}

int WordMatchExact(PyMOLGlobals * G, char *p, char *q, int ignCase)


/* 0 = no match
   non-zero = perfect match  */
{
  while((*p) && (*q)) {
    if(*p != *q) {
      if(!ignCase)
        return 0;
      else if(tolower(*p) != tolower(*q))
        return 0;
    }
    p++;
    q++;
  }
  if((*p) != (*q))
    return 0;
  return 1;
}

int WordMatchNoWild(PyMOLGlobals * G, char *p, char *q, int ignCase)


/* allows for p to match when shorter than q.

Returns:
0 = no match
positive = match out to N characters
negative = perfect match  */
{
  int i = 1;
  while((*p) && (*q)) {
    if(*p != *q) {
      if(ignCase) {
        if(tolower(*p) != tolower(*q)) {
          i = 0;
          break;
        }
      } else {
        i = 0;
        break;
      }
    }
    i++;
    p++;
    q++;
  }
  if((*p) && (!*q))
    i = 0;
  if(i && ((!*p) && (!*q)))     /*exact match gives negative value */
    i = -i;
  return (i);
}

int WordMatch(PyMOLGlobals * G, char *p, char *q, int ignCase)


/* allows for terminal wildcard (*) in p
 * and allows for p to match when shorter than q.

Returns:
0 = no match
positive = match out to N characters
negative = perfect/wildcard match  */
{
  int i = 1;
  register char WILDCARD = '*';
  while((*p) && (*q)) {
    if(*p != *q) {
      if(*p == WILDCARD) {
        i = -i;
        break;
      }
      if(ignCase) {
        if(tolower(*p) != tolower(*q)) {
          i = 0;
          break;
        }
      } else {
        i = 0;
        break;
      }
    }
    i++;
    p++;
    q++;
  }
  if((!*q) && (*p == WILDCARD))
    i = -i;
  if(*p != WILDCARD) {
    if((*p) && (!*q))
      i = 0;
  }
  if(i && ((!*p) && (!*q)))     /*exact match gives negative value */
    i = -i;
  return (i);
}

int WordMatchComma(PyMOLGlobals * G, char *pp, char *qq, int ignCase)


/* allows for comma list in p, also allows wildcards (*) in p */
{
  register char *p = pp, *q = qq;
  register int i = 0;
  register char WILDCARD = '*';
  register char pc, qc;
  register int ic = ignCase;
  int best_i = 0;
  char *q_copy;
  int blank;
  int trailing_comma = 0;

  blank = (!*p);
  q_copy = q;
  while(((*p) || (blank)) && (best_i >= 0)) {
    blank = 0;
    i = 1;
    q = q_copy;
    while((pc = (*p)) && (qc = (*q))) {
      if(pc == ',')
        break;
      if(pc != qc) {
        if(pc == WILDCARD) {
          i = -i;
          break;
        }
        if(ic) {
          if(tolower(pc) != tolower(qc)) {
            i = 0;
            break;
          }
        } else {
          i = 0;
          break;
        }
      }
      p++;
      q++;
      i++;
    }
    if((!*q) && ((*p == WILDCARD) || (*p == ',')))
      i = -i;
    if((*p != WILDCARD) && (*p != ','))
      if((*p) && (!*q))
        i = 0;
    if(i && ((!*p) && (!*q)))   /*exact match */
      i = -i;

    if(i < 0)
      best_i = i;
    else if((best_i >= 0))
      if(i > best_i)
        best_i = i;
    if(best_i >= 0) {
      while(*p) {
        if(*p == ',')
          break;
        p++;
      }
      if(*p == ',') {           /* handle special case, trailing comma */
        if(*(p + 1))
          p++;
        else if(!trailing_comma)
          trailing_comma = 1;
        else
          p++;
      }
    }
  }
  return (best_i);
}

int WordMatchCommaExact(PyMOLGlobals * G, char *p, char *q, int ignCase)


/* allows for comma list in p, no wildcards */
{
  int i = 0;
  int best_i = 0;
  char *q_copy;
  int blank;
  int trailing_comma = 0;

  /*  printf("match? [%s] [%s] ",p,q); */

  blank = (!*p);
  q_copy = q;
  while(((*p) || (blank)) && (best_i >= 0)) {
    blank = 0;
    i = 1;
    q = q_copy;
    while((*p) && (*q)) {
      if(*p == ',')
        break;
      if(*p != *q) {
        if(ignCase) {
          if(tolower(*p) != tolower(*q)) {
            i = 0;
            break;
          }
        } else {
          i = 0;
          break;
        }
      }
      i++;
      p++;
      q++;
    }
    if((!*q) && (*p == ','))
      i = -i;
    if(*p != ',')
      if((*p) && (!*q))
        i = 0;
    if(i && ((!*p) && (!*q)))   /*exact match */
      i = -i;

    if(i < 0)
      best_i = i;
    else if((best_i >= 0))
      if(i > best_i)
        best_i = i;
    if(best_i >= 0) {
      while(*p) {
        if(*p == ',')
          break;
        p++;
      }
      if(*p == ',') {           /* handle special case, trailing comma */
        if(*(p + 1))
          p++;
        else if(!trailing_comma)
          trailing_comma = 1;
        else
          p++;
      }
    }
  }
  /*  printf("result: %d\n",best_i); */

  return (best_i);
}

int WordMatchCommaInt(PyMOLGlobals * G, char *p, int number)
{
  WordType buffer;
  sprintf(buffer, "%d", number);
  return (WordMatchComma(G, p, buffer, 1));
}

int WordIndex(PyMOLGlobals * G, WordType * list, char *word, int minMatch, int ignCase)
{
  int c, i, mi, mc;
  int result = -1;
  c = 0;
  mc = -1;
  mi = -1;
  while(list[c][0]) {
    i = WordMatch(G, word, list[c], ignCase);
    if(i > 0) {
      if(mi < i) {
        mi = i;
        mc = c;
      }
    } else if(i < 0) {
      if((-i) < minMatch)
        mi = minMatch + 1;      /*exact match always matches */
      else
        mi = (-i);
      mc = c;
    }
    c++;
  }
  if((mi > minMatch))
    result = mc;
  return (result);

}

int WordKey(PyMOLGlobals * G, WordKeyValue * list, char *word, int minMatch, int ignCase,
            int *exact)
{
  int c, i, mi, mc;
  int result = 0;
  c = 0;
  mc = -1;
  mi = -1;
  *exact = false;
  while(list[c].word[0]) {
    i = WordMatchNoWild(G, word, list[c].word, ignCase);
    if(i > 0) {
      if(mi < i) {
        mi = i;
        mc = list[c].value;
      }
    } else if(i < 0) {
      *exact = true;
      if((-i) <= minMatch) {
        mi = minMatch + 1;      /*exact match always matches */
      } else
        mi = (-i);
      mc = list[c].value;
    }
    c++;
  }
  if((mi >= minMatch))
    result = mc;
  return (result);
}
