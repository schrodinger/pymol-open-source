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
#include"Util.h"
#include"MemoryDebug.h"
#include<sys/time.h>

static unsigned int UtilStartSec;

void UtilInit(void) {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  UtilStartSec = tv.tv_sec;
}

double UtilGetSeconds(void)
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return((tv.tv_sec-UtilStartSec)+(tv.tv_usec/((double)1000000.0)));
}

char *UtilConcat(char *where,char *what)
{
  while(*what)
	 *where++=*what++;
  *where=0;
  return(where);
}

void UtilCleanStr(char *s) /*remove flanking white and all unprintables*/
{
  char *p,*q;
  p=s;
  q=s;
  while(*p)
	 if(*p>32)
		break;
	 else 
		p++;
  while(*p)
	 if(*p>=32)
		(*q++)=(*p++);
	 else
		p++;
  *q=0;
  while(q>=s)
	 {
		if(*q>32)
		  break;
		else
		  {
			(*q)=0;
			q--;
		  }
	 }
}
void UtilZeroMem(void *ptr,unsigned int howMuch )
{
  unsigned char *c;
  c=ptr;
  while(howMuch--)
	 *c++=0;
}

void *UtilArrayMalloc(unsigned int *dim,int ndim,unsigned int atom_size)
{
  unsigned int size,sum,product;
  unsigned int chunk;
  int a,b,c;
  void *result,**p;
  char *q;
  
  sum = 0;
  for(a=0;a<(ndim-1);a++)
	 {
	 product = dim[0];
	 for(b=1;b<=a;b++)
		product = product * dim[b];
	 sum = sum + product * sizeof(void*);
	 }
  size = atom_size;
  for(a=0;a<ndim;a++)
	 size = size * dim[a];
  size = size + sum;
  result = (void*)mmalloc(size*2);

  if(result)
	 {
		chunk = 1;
		p = result;
		for(c=0;c<(ndim-1);c++)
		  {
			 if(c<(ndim-2))
            {
				chunk = dim[c+1] * sizeof(void*);
           }
			 else
            {
				chunk = dim[c+1] * atom_size;
            }

			 product = dim[0];
			 for(b=1;b<=c;b++)
				product = product * dim[b];
          q = ((char*)p) + product * sizeof(void*); 
			 for(a=0;a<product;a++)
				{
              *p = q;
				  p++;
				  q+=chunk;
				}
		  }
	 }
  return(result);
}
