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
#include"list.h"

/* TEST */

#ifdef _LIST_UT

int main(int argc, char *argv[])
{
  ListInt *int_list,*l;
  List *hl;
  int a,b,c,pass,n_cycles;
  int list[256];

  int_list = (ListInt*)ListNew(10,sizeof(ListInt));
  hl = (List*)int_list;

  c = hl->next_avail;
  l = int_list + c;
  while(l->link) {
    printf(" record %p %d link %d\n",int_list,c,l->link);
    c = l->link;
    l = int_list + c;
  }
  printf(" record %p %d link %d\n",int_list,c,l->link);

  for(a=0;a<100;a++) {
    c = ListElemNew(&int_list);
    printf(" ListElemNew %p %d\n",int_list,c);
    fflush(stdout);
  }
  
  for(a=0;a<256;a++)
    list[a]=0;

  for(pass=0;pass<20;pass++) {
    printf("pass %d\n",pass);
    n_cycles = 10000*pass*pass;
    for(a=0;a<n_cycles;a++) {
      b = rand()&0xFF;
      if(rand()&0x1) {
        c = ListElemNew(&int_list);
        int_list[c].link = list[b];
        list[b] = c;
      } else {
        c = list[b];
        if(c) {
          list[b]=int_list[c].link;
          ListElemFree(int_list,c);
        }
      }
    }
    printf(" cleaning up\n");
    hl = (List*)int_list;
    c = ListLen(int_list,hl->next_avail);
    printf(" available %d\n",c);
    for(a=0;a<256;a++) {
      b = ListLen(int_list,list[a]);
      c+=b;
      ListElemFreeChain(int_list,list[a]);
      list[a]=0;
      printf(" list %d len %d exp %d actual %d \n",
             a,b,c,ListLen(int_list,hl->next_avail));
    }
  }
  return 0;
}

#endif

/* end unit test */


int  ListLen(void *list,int start)
{
  List *I;
  ListElem *hle = NULL;
  register int rec_size,len;
  
  if(start) {
    len = 1;
    I=(List*)list;
    rec_size = I->rec_size;
    hle = (ListElem*)(((char*)I)+rec_size*start);
    while(hle->link) {
      hle = (ListElem*)(((char*)I)+rec_size*hle->link);
      len++;
    }
  } else 
    len = 0;
  return len;
}

void ListPrime(void *list,int start,int stop);

int ListGetNAlloc(void *list)
{
  List *I;

  I = ((List*)list);
  return(vla_get_size(I));
}

int ListElemPush(void *list_ptr_ptr,int elem)
{
  List *I;
  ListElem *hle;
  int start,stop;
  int next;

  I=*((List**)list_ptr_ptr);
  if(!I->next_avail) {
    start = vla_get_size(I);
    vla_check(I,List,start+1);
    (*((List**)list_ptr_ptr))=I;
    stop = vla_get_size(I);
    ListPrime(I,start,stop);
  }
  next = I->next_avail;
  hle = (ListElem*)(((char*)I)+I->rec_size*next);
  I->next_avail = hle->link;
  hle->link = elem;
  return next;
}

int   ListElemPushInt(ListInt **list,int elem,int value)
{
  elem = ListElemPush(list,elem);
  (*list)[elem].value=value;
  return(elem);
}

int   ListElemPopInt(ListInt *list,int elem,int *value)
{
  *value = list[elem].value;
  elem = ListElemPop(list,elem);
  return(elem);
}

int   ListElemGetInt(ListInt *list,int elem,int *value)
{
  *value = list[elem].value;
  elem = list[elem].link;
  return(elem);
}

int   ListElemPurgeInt(ListInt *list,int start,int value)
{
  int last = 0;
  int result = start;
  while(start) {
    if(list[start].value==value) {
      if(!last) {
        result = list[start].link;
      } else {
        list[last].link = list[start].link;
      }
      ListElemFree(list,start);
      break;
    }
    start = list[start].link;
  }
  return(result);
}

int ListElemNew(void *list_ptr_ptr)
{
  return(ListElemPush(list_ptr_ptr,0));
}

int ListElemNewZero(void *list_ptr)
{
  List *I;
  ListElem *hle;
  int start,stop;
  int next;

  I=*((List**)list_ptr);
  if(!I->next_avail) {
    start = vla_get_size(I);
    vla_check(I,List,start+1);
    (*((List**)list_ptr))=I;
    stop = vla_get_size(I);
    ListPrime(I,start,stop);
  }
  next = I->next_avail;
  hle = (ListElem*)(((char*)I)+I->rec_size*next);
  I->next_avail = hle->link;
  os_zero((char*)hle,((char*)hle)+I->rec_size);
  return next;
}
  
void *ListNew(int init_size,int rec_size)
{
  List *I;
  I = (List*)VLAMalloc(init_size+1,rec_size,5,0);
  I->rec_size = rec_size;
  I->next_avail = 0;
  ListPrime(I,1,init_size+1);
  return I;
}

void ListPrime(void *list,int start,int stop)
{  /* sets up next free entry */
  List *I;
  ListElem *hle;
  int a,next,rec_size;

  I = ((List*)list);
  rec_size = I->rec_size;
  next = I->next_avail;
  hle = (ListElem*)(((char*)I)+I->rec_size*(stop-1));
  for(a=stop-1;a>=start;a--) {
    hle->link = next;
    next = a;
    hle = (ListElem*)(((char*)hle)-rec_size);
  }
  I->next_avail = next;
}

int ListElemPop(void *list,int elem)
{
  List *I;
  ListElem *hle;
  int result;
  if(elem) {
    I=(List*)list;
    hle = (ListElem*)(((char*)I)+I->rec_size*elem);
    result = hle->link;
    hle->link = I->next_avail;
    I->next_avail = elem;
  } else 
    result=0;
  return(result);
}

void ListElemFree(void *list,int elem)
{
  List *I;
  ListElem *hle;
  if(elem) {
    I=(List*)list;
    hle = (ListElem*)(((char*)I)+I->rec_size*elem);
    hle->link = I->next_avail;
    I->next_avail = elem;
  }
}

void ListElemFreeChain(void *list,int start)
{
  List *I;
  ListElem *hle = NULL;
  register int rec_size;
  
  if(start) {
    I=(List*)list;
    rec_size = I->rec_size;
    hle = (ListElem*)(((char*)I)+rec_size*start);
    while(hle->link) {
      hle = (ListElem*)(((char*)I)+rec_size*hle->link);
    }
    hle->link = I->next_avail;
    I->next_avail = start;
  }
}

void ListFree(void *list)
{
  VLAFree(list);
}











