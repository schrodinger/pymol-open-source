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

#include"os_predef.h"
#include"os_std.h"

#include"Base.h"
#include"OOMac.h"
#include"Word.h"
#include"Parse.h"
#include"PyMOLObject.h"
#include"MemoryDebug.h"

struct _CWord {
  char Wildcard;
};

CWordList *WordListNew(PyMOLGlobals *G,char *st)
{
  int n_word = 0;
  char *p;
  int len = 0;
  OOCalloc(G,CWordList);

  if(I) {
    p=st;
    /* first, count how many words we have */
    while(*p) {
      if(*p>32) {
        n_word++;
        while((*p)>32) {
          len++;
          p++;
        }
        len++;
      } else
        p++;
    }
    /* allocate the storage we'll need to hold the words */
    {
      I->word = Alloc(char,len);
      I->start = Alloc(char*,n_word);
      
      /* and copy the words */
      
      if(I->word && I->start) {
        char *q = I->word;
        char **q_ptr = I->start;
        p = st;
        while(*p) {
          if(*p>32) {
            *(q_ptr++) = q;
            while((*p)>32) {
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

void WordListFreeP(CWordList *I)
{
  if(I) {
    FreeP(I->word);
    FreeP(I->start);
    FreeP(I);
  }
}

void WordListDump(CWordList *I,char *prefix)
{
  if(I) {
    int a;
    printf(" %s: n_word %d\n",prefix,I->n_word);
    for(a=0;a<I->n_word;a++) {
      printf(" %s: word %d=[%s]\n",prefix,a,I->start[a]);
    }
  }
}

int WordListMatch(PyMOLGlobals *G,CWordList *I,char *name, int ignore_case)
{
  int result = -1;
  if(I) {
    int a;
    for(a=0;a<I->n_word;a++) {
      if(WordMatch(G,I->start[a],name,ignore_case)) {
        result = a;
        break;
      }
    }
  }
  return result;
}

int WordInit(PyMOLGlobals *G)
{
  register CWord *I= NULL;
  
  I = (G->Word = Calloc(CWord,1));
  if(I) {
    I->Wildcard='*';
    return 1;
  } else 
    return 0;

}

void WordFree(PyMOLGlobals *G)
{
  FreeP(G->Word);
}

void WordSetWildcard(PyMOLGlobals *G,char wc)
{
  register CWord *I = G->Word;
  
  I->Wildcard=wc;
}

void WordPrimeCommaMatch(PyMOLGlobals *G,char *p)
{ /* replace '+' with ',' */
  while(*p) { /* this should not be done here... */
    if(*p=='+')
      if(!((*(p+1)==0)||(*(p+1)==',')||(*(p+1)=='+')))
        *p=',';
    p++;
  }
}

int WordMatchExact(PyMOLGlobals *G,char *p,char *q,int ignCase) 

/* 0 = no match
   non-zero = perfect match  */

{
  while((*p)&&(*q))
	 {
		if(*p!=*q)
		  {
			 if(!ignCase)
            return 0;
          else if(tolower(*p)!=tolower(*q))
            return 0;
		  }
		p++;
		q++;
	 }
  if((*p)!=(*q))
    return 0;
  return 1;
}


int WordMatch(PyMOLGlobals *G,char *p,char *q,int ignCase) 
/* allows for terminal wildcard (*) in p
 * and allows for p to match when shorter than q.

Returns:
0 = no match
positive = match out to N characters
negative = perfect/wildcard match  */

{
  int i=1;
  register char WILDCARD = G->Word->Wildcard;
  while((*p)&&(*q))
	 {
		if(*p!=*q)
		  {
			 if(*p==WILDCARD)
				{
				  i=-i;
				  break;
				}
			 if(ignCase)
				{
				  if(tolower(*p)!=tolower(*q))
					 {
						i=0;
						break;
					 }
				}
			 else
				{
				  i=0;
				  break;
				}
		  }
		i++;
		p++;
		q++;
	 }
  if((!*q)&&(*p==WILDCARD))
	 i=-i;
  if(*p!=WILDCARD) {
	 if((*p)&&(!*q))
		i=0;
  }
  if(i&&((!*p)&&(!*q))) /*exact match gives negative value */
	 i=-i;
  return(i);
}

int WordMatchComma(PyMOLGlobals *G,char *pp,char *qq,int ignCase) 
     /* allows for comma list in p, also allows wildcards (*) in p */
{
  register char *p=pp, *q=qq;
  register i=0;
  register char WILDCARD = G->Word->Wildcard;
  register char pc,qc;
  register int ic = ignCase;
  int best_i=0;
  char *q_copy;
  int blank;
  int trailing_comma=0;

  blank = (!*p);
  q_copy=q;
  while(((*p)||(blank))&&(best_i>=0)) {
    blank=0;
    i=1;
    q=q_copy;
    while((pc=(*p))&&(qc=(*q)))
      {
        if(pc==',')
          break;
        if(pc!=qc)
          {
            if(pc==WILDCARD)
              {
                i=-i;
                break;
              }
            if(ic)
              {
                if(tolower(pc)!=tolower(qc))
                  {
                    i=0;
                    break;
                  }
              }
            else 
              {
                i=0;
                break;
              }
          }
        p++;
        q++;
        i++;
      }
    if((!*q)&&((*p==WILDCARD)||(*p==',')))
      i=-i;
    if((*p!=WILDCARD)&&(*p!=','))
      if((*p)&&(!*q))
        i=0;
    if(i&&((!*p)&&(!*q))) /*exact match*/
      i=-i;

    if(i<0)
      best_i=i;
    else if((best_i>=0))
      if(i>best_i)
        best_i=i;
    if(best_i>=0) {
      while(*p) {
        if(*p==',')
          break;
        p++;
      }
      if(*p==',') { /* handle special case, trailing comma */
        if(*(p+1))
          p++;
        else if(!trailing_comma)
          trailing_comma = 1;
        else
          p++;
      }
    }
  }
  return(best_i);
}

int WordMatchCommaExact(PyMOLGlobals *G,char *p,char *q,int ignCase) 
/* allows for comma list in p, no wildcards */
{
  int i=0;
  int best_i=0;
  char *q_copy;
  int blank;
  int trailing_comma=0;

  /*  printf("match? [%s] [%s] ",p,q);*/

  blank = (!*p);
  q_copy=q;
  while(((*p)||(blank))&&(best_i>=0)) {
    blank=0;
    i=1;
    q=q_copy;
    while((*p)&&(*q))
      {
        if(*p==',')
          break;
        if(*p!=*q)
          {
            if(ignCase)
              {
                if(tolower(*p)!=tolower(*q))
                  {
                    i=0;
                    break;
                  }
              }
            else 
              {
                i=0;
                break;
              }
          }
        i++;
        p++;
        q++;
      }
    if((!*q)&&(*p==','))
      i=-i;
    if(*p!=',')
      if((*p)&&(!*q))
        i=0;
    if(i&&((!*p)&&(!*q))) /*exact match*/
      i=-i;

    if(i<0)
      best_i=i;
    else if((best_i>=0))
      if(i>best_i)
        best_i=i;
    if(best_i>=0) {
      while(*p) {
        if(*p==',')
          break;
        p++;
      }
      if(*p==',') { /* handle special case, trailing comma */
        if(*(p+1))
          p++;
        else if(!trailing_comma)
          trailing_comma = 1;
        else
          p++;
      }
    }
  }
  /*  printf("result: %d\n",best_i);*/

  return(best_i);
}

int WordMatchCommaInt(PyMOLGlobals *G,char *p,int number) 
{
  WordType buffer;
  sprintf(buffer,"%d",number);
  return(WordMatchComma(G,p,buffer,1));
}

int WordCompare(PyMOLGlobals *G,char *p,char *q,int ignCase) 
/* all things equal, shorter is smaller */
{
  int result=0;
  if(ignCase) {
    while((*p)&&(*q))	{
      if(*p!=*q) {
        if(tolower(*p)<tolower(*q)) {
          return -1;
        }
        else if(tolower(*p)>tolower(*q)) {
          return 1;
        }
      }
      p++;
      q++;
    }
  } else {
    while((*p)&&(*q))	{
      if(*p!=*q) {
        if(*p<*q) {
          return -1;
        } else if(*p>*q) {
          return 1;
        }
      }
      p++;
      q++;
    }
  }
  
  if((!result)&&(!*p)&&(*q))
    return -1;
  else if((!result)&&(*p)&&(!*q))
    return 1;
  return 0;
}

int WordIndex(PyMOLGlobals *G,WordType *list,char *word,int minMatch,int ignCase)
{
  int c,i,mi,mc;
  int result = -1;
  c=0;
  mc=-1;
  mi=-1;
  while(list[c][0])
	 {
		i=WordMatch(G,word,list[c],ignCase);
		if(i>0)
		  {
			 if(mi<i)
				{
				  mi=i;
				  mc=c;
				}
		  }
		else if(i<0)
		  {
			 if((-i)<minMatch)
				mi=minMatch+1; /*exact match always matches */
			 else
				mi=(-i);
			 mc=c;
		  }
		c++;
	 }
  if((mi>minMatch))
	 result=mc;
  return(result);  

}

int WordKey(PyMOLGlobals *G,WordKeyValue *list,char *word,int minMatch,int ignCase,int *exact)
{
  int c,i,mi,mc;
  int result = 0;
  c=0;
  mc=-1;
  mi=-1;
  *exact = false;
  while(list[c].word[0])
    {
      i=WordMatch(G,word,list[c].word,ignCase);
		if(i>0)
		  {
			 if(mi<i)
				{
				  mi=i;
				  mc=list[c].value;
				}
        }
      else if(i<0)
		  {
          *exact = true;
			 if((-i)<=minMatch) {
            mi=minMatch+1; /*exact match always matches */
          }
          else
            mi=(-i);
          mc=list[c].value;
		  }
		c++;
	 }
  if((mi>=minMatch))
    result=mc;
  return(result);  
}





