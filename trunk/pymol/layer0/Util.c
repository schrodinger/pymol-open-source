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
#include"os_time.h"

#include"Util.h"
#include"MemoryDebug.h"
#include"Err.h"

static unsigned int UtilStartSec;


void UtilInit(void) {
#ifndef WIN32
  struct timeval tv;
  gettimeofday(&tv,NULL);
  UtilStartSec = tv.tv_sec;
#endif
}

double UtilGetSeconds(void)
{
#ifndef WIN32
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return((tv.tv_sec-UtilStartSec)+(tv.tv_usec/((double)1000000.0)));
#else
  return(0.0);
#endif
}

char *UtilConcat(char *where,char *what)
{
  while(*what)
	 *where++=*what++;
  *where=0;
  return(where);
}

void UtilNConcat(char *dst,char *src,int n) { /* copies up to N-1 chars */
  int l;
  l=strlen(dst);
  UtilNCopy(dst+l,src,n-l);
}

void UtilNCopy(char *dst,char *src,int n) { /* copies up to N-1 chars */
  n--;
  while((n--)>=0) {
    if(!*src)
      break;
    else
      *(dst++)=*(src++);
  }
  *dst=0;
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

void UtilSortIndex(int n,void *array,int *x,UtilOrderFn* fOrdered)
{
  int l,a,r,t,i;

  if(n<1) return;
  else if(n==1) { x[0]=0; return; }
  x--;
  for(a=1;a<=n;a++) x[a]=a;
  l=(n>>1)+1;
  r=n;
  while(1) {
	if(l>1)
	  t = x[--l];
	else {
	  t = x[r];
	  x[r] = x[1];
	  if( --r == 1) {
		x[1] = t;
		break;
	  }
	}
	i=l;
	a=l << 1;
	while (a <= r) {
	  if (a < r && (!fOrdered(array,x[a+1]-1,x[a]-1))) a++;
	  if (!fOrdered(array,x[a]-1,t-1)) {
		x[i] = x[a];
		a += (i=a);
	  } else
		a = r + 1;
	}
	x[i] = t;
  }
  x++;
  for(a=0;a<n;a++) x[a]--;
}


void UtilSortInPlace(void *array,int nItem,
					 unsigned int itemSize,
					 UtilOrderFn *fOrdered)

{
  char *tmp;
  int *index;
  int ia;
  int a;
  if(nItem>0)
	 {
	   tmp = Alloc(char,(itemSize*nItem));
	   index = Alloc(int,nItem+1);
	   ErrChkPtr(tmp);
	   ErrChkPtr(index);
	   UtilSortIndex(nItem,array,index,fOrdered);
	   for(a=0;a<nItem;a++) index[a]++; /* ^tricky index adjustment to avoid flag array */
	   for(a=0;a<nItem;a++)
		 {
		   ia = abs(index[a])-1; /* ^ */
		   if(ia!=a)
			 {
			   if(index[a]>0) /* this record not yet copied, so save copy */
				 {
				   memcpy(((char*)tmp  )+(a*itemSize),
						  ((char*)array)+(a*itemSize),
						  itemSize);
				   index[a] = -index[a]; /* set nega-flag */
				 }
			   if(index[ia]<0) /* nega-flag, so record is stored in tmp */
				 memcpy(((char*)array)+(a*itemSize),
						((char*)tmp  )+(ia*itemSize),
						itemSize);
			   else
				 {
				   memcpy(((char*)array)+(a*itemSize),
						  ((char*)array)+(ia*itemSize),
						  itemSize);
				   index[ia] = -index[ia]; 
				   /* nega-flag: record doesn't need to be backed up */
				 }
			 }
		 }
	   mfree(tmp);
	   mfree(index);
	 }
}
