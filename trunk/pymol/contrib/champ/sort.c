#include"sort.h"

void  SortIntIndex(int n,int *v,int *i)
{
  int l,a,r,t,d;

  if(n<1) return;
  else if(n==1) { i[0]=0; return; }

  for(a=0;a<n;a++) i[a]=a;
  
  l=(n>>1);
  r=n-1;
  while(1) {
    if(l>0)
      t = i[--l]; 
    else {  
      t = i[r]; 
      i[r] = i[0]; 
      if( --r == 0) { 
        i[0] = t; 
        break;
      }
    }
    d=l; 
    a= ((l+1) << 1)-1;
    while (a <= r) {
      if ((a < r) && (v[i[a+1]]>v[i[a]])) a++; 
      if (v[i[a]]>v[t]) { 
       i[d] = i[a]; 
       a += (d=a)+1; 
	  } else
       a = r + 1; 
	}
	i[d] = t; 
  }
}

#if 0

/* vidation code */

main()
{
  int a,b,c,d;
  int v[1000];
  int i[1000];
   
  for(a=0;a<1000;a++) { /* list length */
    for(b=0;b<(100-(a/100));b++) { /* number of trials */
      for(c=0;c<a;c++) { /* prime list */
        v[c] = (int)(random()&0xFFF);
        /*        printf("%d ",v[c]);*/
      }
      /*      printf("\n");*/
      SortIntIndex(a,v,i);
      for(c=0;c<a;c++) { /* prime list */
        /*        printf("%d ",v[i[c]]);*/
      }
      for(c=0;c<a-1;c++) { /* vidate */
        for(d=c+1;d<a;d++) {
          if(v[i[c]]>v[i[d]])
            printf("broken!\n");
        }
      }
    }
    printf("list len %d\n",a);

  }
}

#endif

