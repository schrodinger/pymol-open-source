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
#include"MemoryDebug.h"
#include"OOMac.h"
#include"Queue.h"

CQueue *QueueNew(PyMOLGlobals *G,unsigned int mask)
{
  OOAlloc(G,CQueue);
  I->size=mask+1;
  I->ptr=Alloc(char,I->size);
  I->mask=mask;
  I->inp=0;
  I->out=0;
  return(I);
}

void QueueStrIn(CQueue *I,char *c)
{
  int i = I->inp;
  while(*c) {
	 *(I->ptr+i) = *(c++);
	 i=(i+1)&I->mask;
  }
  *(I->ptr+i) = *c;
  i=(i+1)&I->mask;
  I->inp=i; /* important not to do this until null has been written! */
}

int QueueStrCheck(CQueue *I)
{
  return(((I->inp+I->size)-I->out)&I->mask);
}
int QueueStrOut(CQueue *I,char *c)
{
  if(((I->inp+I->size)-I->out)&I->mask) {
	 while(1) {
		*c=*(I->ptr+I->out);
		I->out=(I->out+1)&I->mask;
		if(!*(c++)) {
		  return 1;
      }
	 }
  }
  return 0;
}

void QueueFree(CQueue *I)
{
  FreeP(I->ptr);
  OOFreeP(I);
}


