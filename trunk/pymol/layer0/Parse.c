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

char *ParseNextLine(char *p) {
  while(*p) {
	 if(*p==0xD) { /* Mac or PC */
		if(*(p+1)==0xA) /* PC */
		  p++;
		p++;
		break;
	 }
	 if(*p==0xA) /* Unix */
		{
		  p++;
		  break;
		}
	 p++;
  }
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
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
  }
  *q=0;
  return p;
}
/*========================================================================*/
char *ParseNCopy(char *q,char *p,int n) {  /* n character copy */
  while(*p) {
	 if(!n)
		break;
	 if((*p==0xD)||(*p==0xA)) /* don't copy end of lines */
		break;
	 *(q++)=*(p++);
	 n--;
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

