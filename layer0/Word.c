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
#include<ctype.h>

#include"Word.h"


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
  if(i&&((!*p)&&(!*q))) /*exact match*/
	 i=-i;
  return(i);
}

int WordMatchComma(char *p,char *q,int ignCase) 
     /* allows for comma list in p, also allows wildcards (*) in p */
{
  int i=0;
  int best_i=0;
  char *q_copy;
  q_copy=q;
  while((*p)&&(best_i>=0)) {
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
      if(*p==',') p++;
    }
  }
  return(best_i);
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
  unsigned int result = -1;
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

unsigned int WordChoose(WordType *list, char *word,int minMatch,int ignCase)
{
  int c,i,mi,mc;
  unsigned int result = 0;
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
			 mi=minMatch+1; /*exact match always matches */
			 mc=c;
		  }
		c+=2;
	 }
  if((mi>minMatch))
	 {
		result = (list[mc+1][0]<<24) + 
		  (list[mc+1][1]<<16) + 
		  (list[mc+1][2]<<8) + 
		  (list[mc+1][3]);
	 }
  return(result);  
}



