
#include"os_std.h"
#include"sort.h"
#include"chiral.h"

static int chirality_lookup[256];

#define chi(a,b,c,d) (0xFF&(a<<6|b<<4|c<<2|d))

void ChiralInit(void)
{
  int a;
  for(a=0;a<256;a++) {
    chirality_lookup[a]=0;
  }
  
  /* these orderings have input chirality
     (equivalent to 0,1,2,3) */

  chirality_lookup[chi(0,1,2,3)]=1;
  chirality_lookup[chi(0,2,3,1)]=1;
  chirality_lookup[chi(0,3,1,2)]=1;

  chirality_lookup[chi(1,0,3,2)]=1;
  chirality_lookup[chi(1,3,2,0)]=1;
  chirality_lookup[chi(1,2,0,3)]=1;

  chirality_lookup[chi(2,0,1,3)]=1;
  chirality_lookup[chi(2,1,3,0)]=1;
  chirality_lookup[chi(2,3,0,1)]=1;

  chirality_lookup[chi(3,2,1,0)]=1;
  chirality_lookup[chi(3,1,0,2)]=1;
  chirality_lookup[chi(3,0,2,1)]=1;

  /* these orderings have reverse chirality
     (opposite to 0,1,3,2) */

  chirality_lookup[chi(0,1,3,2)]=-1;
  chirality_lookup[chi(0,2,1,3)]=-1;
  chirality_lookup[chi(0,3,2,1)]=-1;

  chirality_lookup[chi(1,0,2,3)]=-1;
  chirality_lookup[chi(1,3,0,2)]=-1;
  chirality_lookup[chi(1,2,3,0)]=-1;

  chirality_lookup[chi(2,0,3,1)]=-1;
  chirality_lookup[chi(2,1,0,3)]=-1;
  chirality_lookup[chi(2,3,1,0)]=-1;

  chirality_lookup[chi(3,2,0,1)]=-1;
  chirality_lookup[chi(3,1,2,0)]=-1;
  chirality_lookup[chi(3,0,1,2)]=-1;
  
}

int ChiralHandedness(int *a)
{
  int idx[4],ord[4];
  int result;
  SortIntIndex(4,a,idx);
  ord[idx[0]]=0;
  ord[idx[1]]=1;
  ord[idx[2]]=2;
  ord[idx[3]]=3;
  result = chirality_lookup[chi(ord[0],ord[1],ord[2],ord[3])];
  /* printf("%d %d %d %d => %d\n",ord[0],ord[1],ord[2],ord[3],result);*/
  return(result);
  
}

#if 0

main()
{
  int a;
  int val[][4] = {
    {1,10,50,60},
    {2,9,6,12},
    {12,1,2,7},
    {2,10,6,9},
    {1,12,6,3},
  };
  ChiralInit();
  for(a=0;a<5;a++) {
    printf("%d %d %d %d => %d\n\n",
           val[a][0],val[a][1],val[a][2],val[a][3],
           ChiralHandedness(val[a]));

  }
}

#endif
