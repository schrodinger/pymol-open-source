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
#include"Word.h"
#include"Parse.h"

void WordPrimeCommaMatch(char *p)
{ /* replace '+' with ',' */
  while(*p) { /* this should not be done here... */
    if(*p=='+')
      if(!((*(p+1)==0)||(*(p+1)==',')||(*(p+1)=='+')))
        *p=',';
    p++;
  }
}


int WordMatch(char *p,char *q,int ignCase) 
/* allows for terminal wildcard (*) in p
 * and allows for p to match when shorter than q.

Returns:
0 = no match
positive = match out to N characters
negative = perfect match  */

{
  int i=1;
  while((*p)&&(*q))
	 {
		if(*p!=*q)
		  {
			 if(*p=='*')
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
  if((!*q)&&(*p=='*'))
	 i=-i;
  if(*p!='*') {
	 if((*p)&&(!*q))
		i=0;
  }
  if(i&&((!*p)&&(!*q))) /*exact match gives negative value */
	 i=-i;
  return(i);
}

int WordMatchComma(char *p,char *q,int ignCase) 
     /* allows for comma list in p, also allows wildcards (*) in p */
{
  int i=0;
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
    while((*p)&&(*q))
      {
        if(*p==',')
          break;
        if(*p!=*q)
          {
            if(*p=='*')
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
    if((!*q)&&((*p=='*')||(*p==',')))
      i=-i;
    if((*p!='*')&&(*p!=','))
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

int WordMatchCommaInt(char *p,int number) 
{
  WordType buffer;
  sprintf(buffer,"%d",number);
  return(WordMatchComma(p,buffer,1));
}

int WordCompare(char *p,char *q,int ignCase) 
/* all things equal, shorter is smaller */
{
  int result=0;
  while((*p)&&(*q))	{
	if(*p!=*q) {
	  if(ignCase) {
		if(tolower(*p)<tolower(*q)) {
		  result=-1;
		  break;
		}
		else if(tolower(*p)>tolower(*q)) {
		  result=1;
		  break;
		}
	  } else {
		if(*p<*q) {
		  result=-1;
		  break;
		} else if(*p>*q) {
		  result=1;
		  break;
		}
	  }
	}
	p++;
	q++;
  }
  if((!result)&&(!*p)&&(*q))
	result=-1;
  else if((!result)&&(*p)&&(!*q))
	result=1;
  return(result);
}

int WordIndex(WordType *list,char *word,int minMatch,int ignCase)
{
  int c,i,mi,mc;
  int result = -1;
  c=0;
  mc=-1;
  mi=-1;
  while(list[c][0])
	 {
		i=WordMatch(word,list[c],ignCase);
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

int WordKey(WordKeyValue *list,char *word,int minMatch,int ignCase,int *exact)
{
  int c,i,mi,mc;
  int result = 0;
  c=0;
  mc=-1;
  mi=-1;
  *exact = false;
  while(list[c].word[0])
    {
      i=WordMatch(word,list[c].word,ignCase);
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





