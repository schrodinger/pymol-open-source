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

#include"os_memory.h"

#ifdef _os_memory_debug_on

#include<stdio.h>
#include<string.h>

#define GDB_ENTRY

#ifndef NULL
#define NULL (void *)0
#endif

#ifndef false
#define false 0
#endif

#ifndef true
#define true 1
#endif

typedef struct DebugRec {
  struct DebugRec *next;
  char file[64], note[64];
  int line;
  unsigned int size;
  int type;
} DebugRec;

#define HASH(x) ((x>>11)&0x3FF)

static DebugRec *HashTable[1024];
static int InitFlag=true;
static int Count;
static int MaxCount;

DebugRec *OSMemoryRemove(void *ptr);

void OSMemoryInit(void);

void OSMemoryZero(char *p,char *q)
{
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
}



void OSMemoryInit(void)
{
  int a;
  for(a=0;a<1024;a++)
    HashTable[a]=NULL;
  InitFlag=false;
  Count=0;
  MaxCount=0;
}

void OSMemoryDump(void)
{
  int a;
  int cnt=0;
  unsigned int tot  = 0;
  DebugRec *rec;
  if(InitFlag) OSMemoryInit();
  for(a=0;a<1024;a++)
    {
      rec=HashTable[a];
      while(rec)
	{
     tot+=rec->size;
     printf(" OSMemory: @%10p:%7x:%i %s:%i     \n",
            (char*)rec+1,rec->size,rec->type,rec->file,rec->line);
	  rec=rec->next;
	  cnt++;
	}
    }
  printf(" Memory: %d blocks expected, %d found, %d maximum allocated.\n",
			Count,cnt,MaxCount);
  printf(" Memory: current memory allocated %x bytes (%0.1f MB).\n",tot,tot/(1024.0*1024));

}

void OSMemoryHashAdd(DebugRec *rec);
DebugRec *OSMemoryHashRemove(void *ptr);

void OSMemoryHashAdd(DebugRec *rec)
{
  int hash;

  hash=(int)rec;
  hash=HASH(hash);
  rec->next=HashTable[hash];
  HashTable[hash]=rec;
}

DebugRec *OSMemoryHashRemove(void *ptr)
{
  DebugRec *rec,*cur,*last;
  int hash;

  rec=(DebugRec*)ptr;
  rec--;
  hash=(int)rec;
  hash=HASH(hash);
  last=NULL;
  cur=HashTable[hash];
  while(cur)
    {
      if(cur==rec)
	{
	  if(last)
	    last->next=cur->next;
	  else
	    HashTable[hash]=cur->next;
	  break;
	}
      last=cur;
      cur=cur->next;
    }  
  return(cur);
}


void *OSMemoryMalloc(unsigned int size,const char *file,int line,int type)
{
  DebugRec *rec;

  if(InitFlag) OSMemoryInit();
  rec=(DebugRec*)malloc(sizeof(DebugRec)+size);
  if(!rec)
    return(NULL);
  strcpy(rec->file,file);
  rec->line=line;
  rec->size=size;
  rec->type=type;
  OSMemoryHashAdd(rec);
  rec++;
  Count++;
  if(MaxCount<Count) MaxCount=Count;
  return((void *) rec);
}

void *OSMemoryCalloc(unsigned int count,unsigned int size,const char *file,int line,int type)
{
  DebugRec *rec;

  if(InitFlag) OSMemoryInit();
  size = size*count;
  rec=(DebugRec*)calloc(1,sizeof(DebugRec)+size);
  if(!rec)
    return(NULL);
  strcpy(rec->file,file);
  rec->line=line;
  rec->size=size;
  rec->type=type;
  OSMemoryHashAdd(rec);
  rec++;
  Count++;
  if(MaxCount<Count) MaxCount=Count;
  return((void *) rec);
}

void *OSMemoryRealloc(void *ptr,unsigned int size,const char *file,
			 int line,int type)
{
  DebugRec *rec;

  if(InitFlag) OSMemoryInit();
  if((!ptr)&&(!size))
    {
      printf(
	     "OSMemory-ERR: realloc given (NULL,zero) (%s:%i)\n",
	     file,line);
#ifdef GDB_ENTRY
  OSMemoryDump();
  printf("hit ctrl/c to enter debugger\n");
  while(true);
#endif
      exit(EXIT_FAILURE);
    }
  if(!ptr)
      return(OSMemoryMalloc(size,file,line,type));
  else if(!size)
    {
      OSMemoryFree(ptr,file,line,type);
      return(NULL);
    }
  else
    {
      rec=OSMemoryHashRemove(ptr);
      if(!rec)
		  {
			 printf(
					  "OSMemory-ERR: realloc() corrupted tree or bad ptr! (%s:%i @%p)\n",
					  file,line,ptr);
#ifdef GDB_ENTRY
  OSMemoryDump();
  printf("hit ctrl/c to enter debugger\n");
			 while(true);
#endif
			 exit(EXIT_FAILURE);
		  }	
      else
		  {
			 if(rec->type!=type)
				{
				  printf("OSMemory-ERR: ptr is of wrong type: %i!=%i (%s:%i)\n",
							rec->type,type,file,line);
#ifdef GDB_ENTRY
              OSMemoryDump();
				  printf("hit ctrl/c to enter debugger\n");
				  while(true);
#endif
              exit(EXIT_FAILURE);
            }
			 rec=(DebugRec*)realloc(rec,size+sizeof(DebugRec));
			 if(!rec)
				{
				  printf("OSMemory-ERR: realloc() failed reallocation! (%s:%i)\n",
							file,line);
#ifdef GDB_ENTRY
              OSMemoryDump();
				  printf("hit ctrl/c to enter debugger\n");
				  while(true);
#endif
				  exit(EXIT_FAILURE);
				}
			 else
				{
				  OSMemoryHashAdd(rec);
				  rec->size=size;
				  rec++;
				  return((void*)rec);
				}
		  }
    }
  return(ptr);
}

void OSMemoryFree(void *ptr,const char*file,int line,int type)
{
  DebugRec *rec;
  
  if(InitFlag) OSMemoryInit();
  if(!ptr)
    {
      printf("OSMemory-ERR: free() called with NULL pointer (%s:%i)\n",
				 file,line);
#ifdef GDB_ENTRY
  OSMemoryDump();
  printf("hit ctrl/c to enter debugger\n");
		while(true);
#endif
      exit(EXIT_FAILURE);
    }
  rec=OSMemoryHashRemove(ptr);
  if(rec)
    {
      if(rec->type!=type)
		  {
			 printf("OSMemory-ERR: ptr is of wrong type: %i!=%i (%s:%i)\n",
					  rec->type,type,file,line);
#ifdef GDB_ENTRY
          OSMemoryDump();
          printf("hit ctrl/c to enter debugger\n");
			 while(true);
#endif
			 exit(EXIT_FAILURE);
		  }
      free(rec);
    }
  else
    {
		
      printf(
				 "OSMemory-ERR: free(): corrupted tree or bad ptr! (%s:%i @%p)\n",
				 file,line,ptr);
#ifdef GDB_ENTRY
      OSMemoryDump();
      printf("hit ctrl/c to enter debugger\n");
		while(true);
#endif
      exit(EXIT_FAILURE);
    }
  Count--;
}


#endif










