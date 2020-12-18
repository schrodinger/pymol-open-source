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

#include <algorithm>
#include <iterator>

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
  G->Util = pymol::calloc<CUtil>(1);
  G->Util->StartSec = UtilGetSecondsEpoch();
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

/**
 * Get a timestamp in seconds since PyMOL was started
 */
double UtilGetSeconds(PyMOLGlobals *G)
{
  return UtilGetSecondsEpoch() - G->Util->StartSec;
}

/**
 * Get a timestamp in seconds since 1970-01-01
 */
double UtilGetSecondsEpoch()
{
#ifndef _WIN32
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec + (tv.tv_usec / 1e6);
#else
   struct __timeb64 timebuffer;
   _ftime64( &timebuffer );
   return timebuffer.time + (timebuffer.millitm / 1e3);
#endif
}

char *UtilConcat(char *where,const char *what)
{
  while(*what)
	 *(where++)=*(what++);
  *where=0;
  return(where);
}

void UtilConcatVLA(char **vla,ov_size *cc,const char *str)
{
  const char *what;
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

void UtilNPadVLA(char **vla,ov_size *cc,const char *str,ov_size len)
{
  const char *what;
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


void UtilNConcat(char *dst,const char *src,ov_size n) { /* copies up to N-1 chars */
  ov_size l;
  l=strlen(dst);
  if(n>l) {
    UtilNCopy(dst+l,src,n-l);
  }
}

void UtilNCopy(char *dst,const char *src,ov_size n)
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

void UtilNCopyToLower(char *dst,const char *src,ov_size n)
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

/**
 * Removes unprintables and flanking whitespace
 * @param s string to be cleaned
 * @return whitespace-trimmed string with readable characters
 */
std::string UtilCleanStdStr(const std::string& s)
{
  std::string ret;

  auto is_not_whitespace = [](char c) { return c > ' '; };
  auto is_printable = [](char c) { return c >= ' '; };

  auto leading = std::find_if(s.begin(), s.end(), is_not_whitespace);
  auto trailing = std::find_if(s.rbegin(), s.rend(), is_not_whitespace);

  std::copy_if(leading, trailing.base(), std::back_inserter(ret), is_printable);
  return ret;
}

/**
 * Remove ANSI Escape sequences in-place
 */
void UtilStripANSIEscapes(char *s)
{
  for (const char *p = s;; ++p, ++s) {
    while (p[0] == '\033' && p[1] == '[') {
      while (' ' <= p[2] && p[2] < '@') ++p;
      p += 3;
    }
    if (p != s)
      *s = *p;
    if (!p[0])
      break;
  }
}
void UtilStripANSIEscapes(std::string& str)
{
  UtilStripANSIEscapes(&str[0]);
  str.resize(strlen(str.c_str()));
}

void UtilZeroMem(void *ptr,ov_size howMuch)
{
  char *p,*q;
  p=(char*)ptr;
  q=p+howMuch;
  MemoryZero(p,q);
}

void UtilCopyMem(void *dst,const void *src,ov_size howMuch) /* optimize! */
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
  char *p,*q,*p_stop,*q_stop;
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
  void *result;
  char **p;
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
  result = pymol::calloc<char>(size);

  if(result) {
    chunk = 1;
    p = (char**) result;
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
  int a;
  for(a=0;a<n;a++) {
    memcpy(((char*)dst)+(a*rec_size),
           ((char*)src)+(x[a]*rec_size),
           rec_size);
  }
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

void UtilSortIndexGlobals(PyMOLGlobals *G,int n,const void *array,int *x,UtilOrderFnGlobals* fOrdered)
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
  return UtilSemiSortFloatIndexWithNBins(n, n, array, x, forward);
}

int UtilSemiSortFloatIndexWithNBins(int n, int nbins, float *array, int *destx, int forward)
{
  int *start1 = pymol::calloc<int>(n + nbins);
  int ret = UtilSemiSortFloatIndexWithNBinsImpl(start1, n, nbins, array, destx, forward);
  mfree(start1);
  return ret;
}

int UtilSemiSortFloatIndexWithNBinsImpl(int *start1, int n, int nbins, float *array, int *destx, int forward)
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
  /* Since there are two arrays, this guarentees that there */
  /* will be enough memory to hold all indexes.  If there are many collisions, */
  /* the next1 array will hold a link to most of the indexes, which are traversed */
  /* when the first index is found in start1.  If there are few collisions, then */
  /* the majority of the start1 array is used. The total number of items used in */
  /* both arrays will always be the number of values, i.e., n. */
  /* 9/9/14: BB - added start1 and nbins argument
     start1 - pre-allocated memory
     nbins - allows the first array to be controled as the number of bins, 
             to match how CGORenderGLAlpha() sorts its triangles.
  */
  int ok = true;
  if(n>0) {
    float min,max,*f,v;
    float range, scale;
    int a;
    int *next1;
    int idx1;

    CHECKOK(ok, start1);
    if (!ok){
      return false;
    }

    next1 = start1 + nbins;
    max = (min = array[0]);
    f = array + 1;
    for(a=1;a<n;a++) {
      v = *(f++);
      if(max<v) max=v;
      if(min>v) min=v;
    }
    range = (max-min)/.9999F; /* for boundary conditions */
    if(range<R_SMALL8) { 
      for(a=0;a<n;a++)
        destx[a] = a;
    } else {
      scale = nbins/range;
      f = array;
      /* hash by value (actually binning) */
      if(forward) {
        for(a=0;a<n;a++) {
          idx1 = (int)((*(f++)-min)*scale);
          next1[a] = start1[idx1];
          start1[idx1] = a+1;
        }
      } else {
        for(a=0;a<n;a++) {
          idx1 = (nbins-1) - (int)((*(f++)-min)*scale);
          next1[a] = start1[idx1];
          start1[idx1] = a+1;
        }
      }
      /* now read out */
      {
        int c=0;
        int cur1;        
        a=0;
        while(a<nbins) {
          if( (cur1 = start1[a]) ) {
            idx1 = cur1 - 1;
            while(1) {
              destx[c++] = idx1;
              if(! (cur1 = next1[idx1]))
                break;
              idx1 = cur1 - 1;
            }
          }
          a++;
        }
      }
    }
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
	   tmp = pymol::malloc<char>((itemSize*nItem));
	   index = pymol::malloc<int>(nItem+1);
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
