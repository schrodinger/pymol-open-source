#include "ov_port.h"
#include "ov_utility.h"

void ov_utility_zero_range(ov_pointer start, ov_pointer stop)
{
  char *p = (char*)start;
  char *q = (char*)stop;
#if 1
  if(q-p)
    memset(p,0,q-p);
#else
  register unsigned long count;
  register long *a;
  int mask;

  count = q-p;
  mask=sizeof(long)-1;
  /* get us word aligned */
  while(count&&(((int)p)&mask)) {
    count--;
    *p++=0;
  }
  a=(long*)p;
  /* now blank efficiently */
  while(count>(sizeof(long)*16))
	 {
		count-=(sizeof(long)*16);
		*a++=0;
		*a++=0;
		*a++=0;
		*a++=0;

		*a++=0;
		*a++=0;
		*a++=0;
		*a++=0;

		*a++=0;
		*a++=0;
		*a++=0;
		*a++=0;

		*a++=0;
		*a++=0;
		*a++=0;
		*a++=0;
	 }
  p=(char*)a;
  while(count>0)
	 {
		*p++=0;
		count--;
	 }
#endif

}

void ov_utility_zero_bytes(ov_pointer start, ov_size n_bytes)
{
  ov_utility_zero_range(start, ((char*)start)+n_bytes);
}
