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
#include"Parse.h"

#ifndef _PYMOL_INLINE

char *ParseNextLine(char *p) {
  register char ch;
  const char mask = -16; /* 0xF0 */
  while((mask & p[0]) &&
	(mask & p[1]) &&
	(mask & p[2]) &&
	(mask & p[3])) /* trusting short-circuit to avoid overrun */
    p+=4;
  while( (ch=*p) ) {
    p++;
    if(ch==0xD) { /* Mac or PC */
      if((*p)==0xA) /* PC */
        return p+1;
      return p;
    } else if(ch==0xA) { /* Unix */
      return p;
    }
  }
  return p;
}
/*========================================================================*/

char *ParseNCopy(char *q,char *p,int n) {  /* n character copy */
  register char ch;
  while( (ch=*p) ) {
    if((ch==0xD)||(ch==0xA)) /* don't copy end of lines */
      break;
    if(!n)
      break;
    n--;
    p++;
    *(q++)=ch;
  }
  *q=0;
  return p;
}

#endif

char *ParseSkipEquals(char *p)
{
  while(*p) {
    if(*p!='=')
      p++;
    else
      break;
  }

  if(*p) {
    p++;
    while(*p) {
      if(*p<33)  /* skip whitespace */
        p++; 
      else
        break;
    }
  }
  return p;
}


/*========================================================================*/
char *ParseIntCopy(char *q,char *p,int n) { /* integer copy */
  while(*p) {
	 if((*p==0xD)||(*p==0xA)) /* don't skip end of lines */
      break;
	 if(*p<=32||!((*p>='0')&&(*p<='9')))
		p++;
	 else
		break;
  }
  while(*p) {
	 if(*p<=32)
		break;
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
    if(!((*p>='0')&&(*p<='9')))
      break;
    *(q++)=*(p++);
    n--;
  }
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseAlphaCopy(char *q,char *p,int n) { /* integer copy */
  while(*p) {
	 if((*p==0xD)||(*p==0xA)) /* don't skip end of lines */
      break;
	 if(*p<=32||!(
                 ((*p>='A')&&(*p<='Z'))||
                 ((*p>='a')&&(*p<='z'))))
		p++;
	 else
		break;
  }
  while(*p) {
	 if(*p<=32)
		break;
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
    if(!( ((*p>='A')&&(*p<='Z'))||
          ((*p>='a')&&(*p<='z'))))
      break;
    *(q++)=*(p++);
    n--;
  }
  *q=0;
  return p;
}

/*========================================================================*/
char *ParseWordCopy(char *q,char *p,int n) { /* word copy */
  while(*p) {
	 if((*p==0xD)||(*p==0xA)) /* don't skip end of lines */
      break;
	 if(*p<=32) 
		p++;
	 else
		break;
  }
  while(*p) {
	 if(*p<=32)
		break;
	 if(!n) {
       while(*p>32) /* finish scanning word, but don't copy into field */
         p++;
		break;
     }
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
  } 
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseWord(char *q,char *p,int n) { /* word copy, across lines */
  while(*p) {
	 if(*p<=32) 
		p++;
	 else
		break;
  }
  while(*p) {
	 if(*p<=32)
		break;
	 if(!n)
		break;
	 *(q++)=*(p++);
	 n--;
  }
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseNTrim(char *q,char *p,int n) {  /* n character trimmed copy */
  char *q_orig = q;
  while(*p) {
	 if((*p==0xD)||(*p==0xA)) /* don't skip end of lines */
      break;
	 if(*p<=32) {
		p++;
      n--;
    } else
		break;
  }
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
  }
  while(q>q_orig) {
    if(*(q-1)<=32) 
      q--;
    else
      break;
  }
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseNTrimRight(char *q,char *p,int n) {  /* n character trimmed copy */
  char *q_orig = q;
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
  }
  while(q>q_orig) {
    if(*(q-1)<=32) 
      q--;
    else
      break;
  }
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseCommaCopy(char *q,char *p,int n) {  /* n character copy up to comma */
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
    if(*p==',')
      break;
	 *(q++)=*(p++);
	 n--;
  }
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseNSkip(char *p,int n) {  /* n character skip */
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* stop at newlines */
		break;
    p++;
	 n--;
  }
  return p;
}

