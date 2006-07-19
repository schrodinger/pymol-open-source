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
#ifndef _H_Parse
#define _H_Parse

char *ParseWordCopy(char *dst,char *src,int n);
char *ParseWord(char *dst,char *src,int n);
char *ParseNTrim(char *q,char *p,int n);
char *ParseNTrimRight(char *q,char *p,int n);
char *ParseNSkip(char *p,int n);
char *ParseCommaCopy(char *q,char *p,int n);
char *ParseSkipEquals(char *p);
char *ParseIntCopy(char *q,char *p,int n);
char *ParseAlphaCopy(char *q,char *p,int n);

#ifdef _PYMOL_INLINE

__inline__ static char *ParseNextLine(char *p)
{
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

__inline__ static char *ParseNCopy(char *q,char *p,int n) {  /* n character copy */
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

#else

char *ParseNextLine(char *p);
char *ParseNCopy(char *dst,char *src,int n);

#endif

#endif
