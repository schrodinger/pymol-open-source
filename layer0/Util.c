/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
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
#include"os_time.h"

#include"Util.h"
#include"MemoryDebug.h"
#include"Err.h"

struct _CUtil {
  double StartSec;
};


int UtilInit(PyMOLGlobals *G) 
{
  G->Util = Calloc(CUtil,1);
  G->Util->StartSec = UtilGetSeconds(G);
  return 1;
}

void UtilFree(PyMOLGlobals *G) 
{
  FreeP(G->Util);
}

int UtilShouldWePrintQuantity(int quantity)
{
  if(quantity<10)
    return 1;
  if((quantity>0)&&(quantity<0x07FFFFFF)) /* avoids overflow, just in case */ {
    int factor = 10;
    while((factor*10)<quantity)
      factor *= 10;
    return ((quantity/factor)*factor == quantity);
  }
  return 0;
}
int UtilCountStringVLA(char *vla)
{
  int result=0;
  int cc;
  if (vla) {
    cc=VLAGetSize(vla);
    while(cc--) {
      if(!*vla) 
        result++;
      vla++;
    }
  }
  return(result);
}

double UtilGetSeconds(PyMOLGlobals *G)
{
#ifndef _WIN32
  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return((tv.tv_sec+(tv.tv_usec/((double)1000000.0)))-G->Util->StartSec);
  /* END PROPRIETARY CODE SEGMENT */
#else
   struct _timeb timebuffer;
   _ftime( &timebuffer );
   return((timebuffer.time+(timebuffer.millitm/((double)1000.0)))-G->Util->StartSec);
#endif
}

char *UtilConcat(char *where,char *what)
{
  while(*what)
	 *(where++)=*(what++);
  *where=0;
  return(where);
}

void UtilConcatVLA(char **vla,ov_size *cc,char *str)
{
  char *what;
  char *where;
  ov_size len;

  len=strlen(str);
  VLACheck((*vla),char,len+*cc+1); 
  where = (*cc)+(*vla);
  what = str;
  while(*what)
	 *(where++)=*(what++);
  *where=0;
  *(cc)+=len;
}

void UtilNPadVLA(char **vla,ov_size *cc,char *str,ov_size len)
{
  char *what;
  char *where;
  ov_size n = 0;
  VLACheck((*vla),char,len + *cc +1); 
  where = (*cc)+(*vla);
  what = str;
  while(*what) {
    if(n>=len) 
      break;
    *(where++)=*(what++);
    n++;
  }
  while(n<len) {
    *(where++) = ' ';
    n++;
  }
  *where=0;
  *(cc)+=len;
}

void UtilFillVLA(char **vla,ov_size *cc,char what,ov_size len)
{
  char *where;
  VLACheck((*vla),char,len+(*cc)+1); 
  where = (*cc)+(*vla);
  *(cc)+=len;
  while((len--)>0)
    *(where++)=what;
  *where=0;
}


void UtilNConcat(char *dst,char *src,ov_size n) { /* copies up to N-1 chars */
  ov_size l;
  l=strlen(dst);
  if(n>l) {
    UtilNCopy(dst+l,src,n-l);
  }
}

void UtilNCopy(char *dst,char *src,ov_size n) 
{ /* copies up to N-1 chars */
  if(n--) {
    while(n--) {
      if(!*src)
        break;
      else
        *(dst++)=*(src++);
    }
  }
  *dst=0;
}

void UtilNCopyToLower(char *dst,char *src,ov_size n)
{
  if(n--) {
    while(n--) {
      if(!*src)
        break;
      else
        *(dst++)=tolower(*(src++));
    }
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

void UtilZeroMem(void *ptr,ov_size howMuch)
{
  char *p,*q;
  p=(char*)ptr;
  q=p+howMuch;
  MemoryZero(p,q);
}

void UtilCopyMem(void *dst,void *src,ov_size howMuch) /* optimize! */
{
  /* need to determine the memory is non-overlapping.  If so, then use memcpy. */
  char *c,*d;
  c=(char*)dst;
  d=(char*)src;
  while(howMuch--)
	 *(c++)=*(d++);
}

void UtilExpandArrayElements(void *src,void *dst,int n_entries,int old_rec_size,int new_rec_size)
{
  /* simple but ineffient byte-based copy */
  register char *p,*q,*p_stop,*q_stop;
  int a;
  for(a=0;a<n_entries;a++) {
    p=((char*)src)+(old_rec_size*a); 
    p_stop=p+old_rec_size; 
    q=((char*)dst)+(new_rec_size*a); 
    q_stop=q+new_rec_size; 
    while(p!=p_stop) {
      *(q++)=*(p++);
    }
    while(q!=q_stop) {
      *(q++)=0;
    }
  }
}

void *UtilArrayCalloc(unsigned int *dim,ov_size ndim,ov_size atom_size)
{
  ov_size size;
  ov_size sum,product;
  ov_size chunk;
  ov_size a,b,c;
  void *result,**p;
  char *q;
  
  sum = 0;
  for(a=0;a<(ndim-1);a++) {
    product = dim[0];
    for(b=1;b<=a;b++)
      product = product * dim[b];
    sum = sum + product * sizeof(void*);
  }
  size = atom_size;
  for(a=0;a<ndim;a++)
	 size = size * dim[a];
  size = size + sum;
  result = (void*)mcalloc(size*2,1); /* what is this *2 for ??? */

  if(result) {
    chunk = 1;
    p = result;
    for(c=0;c<(ndim-1);c++) {
      if(c<(ndim-2)) {
        chunk = dim[c+1] * sizeof(void*);
      } else {
        chunk = dim[c+1] * atom_size;
      }
      
      product = dim[0];
      for(b=1;b<=c;b++)
        product = product * dim[b];
      q = ((char*)p) + product * sizeof(void*); 
      for(a=0;a<product;a++) {
        *p = q;
        p++;
        q+=chunk;
      }
    }
  }
  return(result);
}

void UtilApplySortedIndices(int n,int *x, int rec_size, void *src, void *dst)
{
  register int a;
  for(a=0;a<n;a++) {
    memcpy(((char*)dst)+(a*rec_size),
           ((char*)src)+(x[a]*rec_size),
           rec_size);
  }
}


void UtilSortIndex(int n,void *array,int *x,UtilOrderFn* fOrdered)
{
  register int l,a,r,t,i;

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

void UtilSortIndexGlobals(PyMOLGlobals *G,int n,void *array,int *x,UtilOrderFnGlobals* fOrdered)
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
	  if (a < r && (!fOrdered(G,array,x[a+1]-1,x[a]-1))) a++;
	  if (!fOrdered(G,array,x[a]-1,t-1)) {
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

#define MAX_BIN = 100

#ifndef R_SMALL8
#define R_SMALL8 0.00000001F
#endif

int UtilSemiSortFloatIndex(int n,float *array,int *x, int forward)
{
  /* approximate sort, for quick handling of transparency values */
  /* this sort uses 2 arrays start1 and next1 to keep track of */
  /* the indexes.  The values in start1 are set to the index */
  /* relative to the array value within the min/max values.  If */
  /* there is a collision, the value in next1 is set to the value */
  /* that is collided, and start1[idx] is set to the index plus 1 (a+1) */
  /* This makes it easy to go through the 2 arrays and write into the */
  /* x array the approximate order of the floating point values in array */
  /* by indexes. */
  /* Since there are two arrays of n length, this guarentees that there */
  /* will be enough memory to hold all indexes.  If there are many collisions, */
  /* the next1 array will hold a link to most of the indexes, which are traversed */
  /* when the first index is found in start1.  If there are few collisions, then */
  /* the majority of the start1 array is used. The total number of items used in */
  /* both arrays will always be the number of values, i.e., n. */
  int ok = true;
  if(n>0) {
    register float min,max,*f,v;
    register float range, scale;
    register int a;
    register int *start1;
    register int *next1;
    register int idx1;
    register int n_minus_one;

    start1 = Calloc(int,n*2);
    CHECKOK(ok, start1);
    if (!ok){
      return false;
    }

    next1 = start1 + n;
    max = (min = array[0]);
    f = array + 1;
    for(a=1;a<n;a++) {
      v = *(f++);
      if(max<v) max=v;
      if(min>v) min=v;
    }
    range = (max-min)*1.0001F; /* for boundary conditions */
    if(range<R_SMALL8) { 
      for(a=0;a<n;a++)
        x[a] = a;
    } else {
      scale = n/range;
      f = array;
      /* hash by value (actually binning) */
      if(forward) {
        for(a=0;a<n;a++) {
          idx1 = (int)((*(f++)-min)*scale);
          next1[a] = start1[idx1];
          start1[idx1] = a+1;
        }
      } else {
        n_minus_one = n-1;
        for(a=0;a<n;a++) {
          idx1 = n_minus_one-(int)((*(f++)-min)*scale);
          next1[a] = start1[idx1];
          start1[idx1] = a+1;
        }
      }
      /* now read out */
      {
        register int c=0;
        register int cur1;        
        a=0;
        while(a<n) {
          if( (cur1 = start1[a]) ) {
            idx1 = cur1 - 1;
            while(1) {
              x[c++] = idx1;
              if(! (cur1 = next1[idx1]))
                break;
              idx1 = cur1 - 1;
            }
          }
          a++;
        }
      }
    }
    mfree(start1);
  }
  return true;
}

void UtilSortInPlace(PyMOLGlobals *G,void *array,int nItem,
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
	   ErrChkPtr(G,tmp);
	   ErrChkPtr(G,index);
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
