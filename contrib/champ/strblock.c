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

#include"os_std.h"
#include"os_memory.h"
#include"mac.h"
#include"vla.h"
#include"strblock.h"

/* TEST */

#ifdef _STRLIST_UT

int main(int argc, char *argv[])
{
  char *sb;
  char buffer[256];
  int a,b,c;
  int index;
  
  buffer[0]=0;
  sb = StrBlockNew(10);
  b = rand()&0xFF;
  for(a=0;a<100000;a++) {
    index = StrBlockNewStr(&sb,NULL,b);
    for(c=0;c<b;c++) {
      sb[c+index]='a';
    }
  }
  return 0;
}

#endif

/* end unit test */

/* right now just a braindead, memory-leaky system */

int StrBlockNewStr(char **list_ptr,char *st,int len)
{
  StrBlock *I;
  register int a;
  register char *p,*q;
  char *str;
  int result, new_extent;
  I=*((StrBlock**)list_ptr);
  new_extent = len + 1 + I->next_unused;
  vla_check(I,StrBlock,new_extent);
  (*((StrBlock**)list_ptr))=I;
  result = I->next_unused;
  str = (char*)(((char*)I)+I->next_unused);
  if(st) {
    p = st;
    q = str;
    for(a=0;a<len;a++) 
      *(q++)=*(p++);
  } else
    str[0]=0;
  str[len]=0;
  I->next_unused = new_extent;
  return result;
}

char *StrBlockNew(int init_size)
{
  StrBlock *I;
  I = (StrBlock*)VLAMalloc(init_size+sizeof(StrBlock),sizeof(char),5,0);
  I->next_unused = sizeof(StrBlock);
  return (char*)I;
}

void StrBlockFreeStr(char *list,int index)
{
  /* just leak within the array */
}

void StrBlockFree(char *list)
{
  VLAFree(list);
}

