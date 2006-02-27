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
/* This file can be compiled under C as a .c file, or under C++ as a .cc file*/

#include"os_predef.h"

#ifdef __cplusplus
#include<new.h>
#endif

#include"MemoryDebug.h"
#include"MemoryCache.h"

#define GDB_ENTRY

void *MemoryReallocForSureSafe(void *ptr, unsigned int newSize, unsigned int oldSize)
{
  if(newSize<oldSize) {
    float *tmp = mmalloc(newSize);
    if(tmp && newSize && oldSize) {
        memcpy(tmp, ptr, newSize);
    }
    FreeP(ptr);
    return tmp;
  } else {
    return mrealloc(ptr,newSize);
  }
}

void *MemoryReallocForSure(void *ptr, unsigned int newSize) /* unsafe -- replace with above */
{
  float *tmp = mmalloc(newSize);
  if(tmp)
    memcpy(tmp,ptr,newSize);
  FreeP(ptr);
  return tmp;
}

static void DieOutOfMemory(void) 
{
  printf("****************************************************************************\n");
  printf("*** EEK!  PyMOL just ran out of memory and crashed.  To get around this, ***\n");
  printf("*** you may need to reduce the quality, size, or complexity of the scene ***\n");
  printf("*** that you are viewing or rendering.    Sorry for the inconvenience... ***\n");
  printf("****************************************************************************\n");
#ifdef GDB_ENTRY
  abort();
#endif

  exit(EXIT_FAILURE);
}

void MemoryZero(char *p,char *q)
{
#if 1
  if(q-p)
    memset(p,0,q-p);
#else
  register unsigned long count;
  register long *a;
  int mask;
  /*  fprintf(stderr,"MemoryZero: start %p stop %p\n",p,q);
      fflush(stderr);
  */

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

void *VLAExpand(void *ptr,unsigned int rec)
{
  VLARec *vla;
  char *start,*stop;
  unsigned int soffset=0;
  vla = &(((VLARec*)ptr)[-1]);
  if(rec>=vla->nAlloc) {
    if(vla->autoZero)
      soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
    vla->nAlloc = ((unsigned int)(rec*vla->growFactor))+1;
    if(vla->nAlloc<=rec) vla->nAlloc = rec+1;
    vla=(void*)mrealloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec));
    if(!vla) {
      printf("VLAExpand-ERR: realloc failed.\n");
      DieOutOfMemory();
    }
    if(vla->autoZero) {
      start = ((char*)vla) + soffset;
      stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
      MemoryZero(start,stop);
    }
  }
  return((void*)&(vla[1]));
}

#ifdef _MemoryCache_ON
void *VLACacheExpand(PyMOLGlobals *G,void *ptr,unsigned int rec,int thread_index,int block_id)
{
  VLARec *vla;
  char *start,*stop;
  unsigned int soffset=0;
  vla = &(((VLARec*)ptr)[-1]);
  if(rec>=vla->nAlloc)
	 {
		if(vla->autoZero)
		  soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		vla->nAlloc = ((unsigned int)(rec*vla->growFactor))+1;
        if(vla->nAlloc<=rec) vla->nAlloc = rec+1;
		vla=(void*)_MemoryCacheRealloc(G,vla,
                                       (vla->recSize*vla->nAlloc)+sizeof(VLARec),
                                       thread_index,block_id MD_FILE_LINE_Call);
		if(!vla)
		  {
			 printf("VLAExpand-ERR: realloc failed.\n");
          DieOutOfMemory();
		  }
		if(vla->autoZero)
		  {
			 start = ((char*)vla) + soffset;
			 stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
			 MemoryZero(start,stop);
		  }
	 }
  return((void*)&(vla[1]));
}
#endif

#ifndef _MemoryDebug_ON
void *VLAMalloc(unsigned int initSize,unsigned int recSize,unsigned int growFactor,int autoZero)
#else
void *_VLAMalloc(const char *file,int line,unsigned int initSize,
                 unsigned int recSize,unsigned int growFactor,int autoZero)
#endif
{
  VLARec *vla;
  char *start,*stop;
#ifndef _MemoryDebug_ON
  vla=(void*)mmalloc((initSize*recSize)+sizeof(VLARec));
#else
  vla=MemoryDebugMalloc((initSize*recSize)+sizeof(VLARec),file,line,_MDPointer);
#endif

  if(!vla) {
    printf("VLAMalloc-ERR: realloc failed\n");
    DieOutOfMemory();
  }
  vla->nAlloc=initSize;
  vla->recSize=recSize;
  vla->growFactor=(1.0F + growFactor*0.1F);
  vla->autoZero=autoZero;
  if(vla->autoZero) {
    start = ((char*)vla)+sizeof(VLARec);
    stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
    MemoryZero(start,stop);
  }
  return((void*)&(vla[1]));
}


#ifdef _MemoryCache_ON
#ifndef _MemoryDebug_ON
void *VLACacheMalloc(PyMOLGlobals *G,unsigned int initSize,unsigned int recSize,
                     unsigned int growFactor,int autoZero,int thread,int id)
#else
void *_VLACacheMalloc(PyMOLGlobals *G,const char *file,int line,
                      unsigned int initSize,unsigned int recSize,unsigned int growFactor,
                      int autoZero,int thread,int id)
#endif
{
  VLARec *vla;
  char *start,*stop;

  vla=(void*)_MemoryCacheMalloc(G,(initSize*recSize)+sizeof(VLARec),
                                thread,id 
                                MD_FILE_LINE_Nest);

  if(!vla)
	 {
		printf("VLAMalloc-ERR: realloc failed\n");
      DieOutOfMemory();
	 }
  vla->nAlloc=initSize;
  vla->recSize=recSize;
  vla->growFactor=(1.0F + growFactor*0.1F);
  vla->autoZero=autoZero;
  if(vla->autoZero)
	 {
		start = ((char*)vla)+sizeof(VLARec);
		stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}

#endif

void VLAFree(void *ptr)
{
  VLARec *vla;
  if(!ptr)
	 {
		printf("VLAFree-ERR: tried to free NULL pointer!\n");
		exit(EXIT_FAILURE);
	 }
  vla = &(((VLARec*)ptr)[-1]);
  mfree(vla);
}

#ifdef _MemoryCache_ON
void VLACacheFree(PyMOLGlobals *G,void *ptr,int thread,int id,int force)
{
  VLARec *vla;
  if(!ptr)
	 {
		printf("VLAFree-ERR: tried to free NULL pointer!\n");
		exit(EXIT_FAILURE);
	 }
  vla = &(((VLARec*)ptr)[-1]);
  _MemoryCacheFree(G,vla,thread,id,force MD_FILE_LINE_Call);
}
#endif

unsigned int VLAGetSize(void *ptr)
{
  VLARec *vla;
  vla = &((VLARec*)ptr)[-1];
  return(vla->nAlloc);
}

void *VLANewCopy(void *ptr)
{
  VLARec *vla,*new_vla;
  unsigned int size;
  vla = &((VLARec*)ptr)[-1];
  size = (vla->recSize*vla->nAlloc)+sizeof(VLARec);
  new_vla=(void*)mmalloc(size);
  if(!new_vla)
	 {
		printf("VLACopy-ERR: mmalloc failed\n");
		exit(EXIT_FAILURE);
	 }
  else
    {
      memcpy(new_vla,vla,size);
    }
  return((void*)&(new_vla[1]));
}
#ifdef _MemoryCache_ON
void *VLACacheSetSize(PyMOLGlobals *G,void *ptr,unsigned int newSize,int group_id,int block_id)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->autoZero) {
	 soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
  }
  vla->nAlloc = newSize;
  vla=(void*)_MemoryCacheRealloc(G,vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec),
                                 group_id,block_id MD_FILE_LINE_Call);
  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed.\n");
      DieOutOfMemory();
	 }
  if(vla->autoZero)
	 {
      start = ((char*)vla)+soffset;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		if(start<stop)
		  MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}
void *VLACacheSetSizeForSure(PyMOLGlobals *G,void *ptr,unsigned int newSize,int group_id,int block_id)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->autoZero) {
	 soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
  }
  if(newSize<vla->nAlloc) {
    vla->nAlloc = newSize;
    vla=(void*)_MemoryCacheShrinkForSure(G,vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec),
                                         group_id,block_id MD_FILE_LINE_Call);
  } else {
    vla->nAlloc = newSize;
    vla=(void*)_MemoryCacheRealloc(G,vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec),
                                   group_id,block_id MD_FILE_LINE_Call);
  }
  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed.\n");
      DieOutOfMemory();
	 }
  if(vla->autoZero)
	 {
      start = ((char*)vla)+soffset;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		if(start<stop)
		  MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}
#endif


void *VLASetSize(void *ptr,unsigned int newSize)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->autoZero) {
	 soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
  }
  vla->nAlloc = newSize;
  vla=(void*)mrealloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec));
  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed.\n");
      DieOutOfMemory();
	 }
  if(vla->autoZero)
	 {
      start = ((char*)vla)+soffset;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
		if(start<stop)
		  MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}

void *VLASetSizeForSure(void *ptr,unsigned int newSize)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->autoZero) {
    soffset = sizeof(VLARec)+(vla->recSize*vla->nAlloc);
  }
  if(newSize<vla->nAlloc) {
    vla=MemoryReallocForSureSafe(vla,
                                 (vla->recSize*newSize)+sizeof(VLARec),
                                 (vla->recSize*vla->nAlloc)+sizeof(VLARec));
    vla->nAlloc = newSize;
  } else {
    vla->nAlloc = newSize;
    vla=(void*)mrealloc(vla,(vla->recSize*vla->nAlloc)+sizeof(VLARec));
  }
  if(!vla) {
    printf("VLASetSize-ERR: realloc failed.\n");
    DieOutOfMemory();
  }
  if(vla->autoZero) {
    start = ((char*)vla)+soffset;
    stop = ((char*)vla)+sizeof(VLARec)+(vla->recSize*vla->nAlloc);
    if(start<stop)
      MemoryZero(start,stop);
  }
  return((void*)&(vla[1]));
}

#ifdef _MemoryDebug_ON

#ifndef NULL
#define NULL (void *)
#endif

#ifndef false
#define false 0
#endif

#ifndef true
#define true 1
#endif

#define _NO_GDB_ENTRY


#ifdef __cplusplus
extern "C" {
#endif


typedef struct DebugRec {
  struct DebugRec *next;
  char file[32], note[64];
  int line;
  size_t size;
  int type;
} DebugRec;

#define HASH(x) ((x>>11)&0x3FF)

static DebugRec *HashTable[1024];
static int InitFlag=true;
static int Count;
static int MaxCount;

DebugRec *MemoryDebugRemove(void *ptr);

void MemoryDebugInit(void);

void MemoryDebugInit(void)
{
  int a;
  for(a=0;a<1024;a++)
    HashTable[a]=NULL;
  InitFlag=false;
  Count=0;
  MaxCount=0;
}
void MemoryDebugRegister(void *addr,const char *note,
			 const char *file,int line)
{
  DebugRec *rec,*str;
  int hash;

  if(InitFlag) MemoryDebugInit();

  rec=(DebugRec*)malloc(sizeof(DebugRec)+strlen(note));
  if(!rec)
    {
      printf("MemoryDebugRegister-ERR: memory allocation failure\n"); 
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
  DieOutOfMemory();
    }
  rec->size=(size_t)addr;
  rec->type=_MDMarker;
  rec->line=line;
  strcpy(rec->file,file);
  str=rec+1;
  strcpy((char*)str,note);
  hash=(int)addr;
  hash=HASH(hash);
  rec->next=HashTable[hash];
  HashTable[hash]=rec;
  Count++;
  if(MaxCount<Count) MaxCount=Count;
}

void MemoryDebugForget(void *addr,const char *file,int line)
{
  DebugRec *rec = NULL;
  DebugRec *cur,*last;
  int hash;

  if(InitFlag) MemoryDebugInit();
  hash=(int)addr;
  hash=HASH(hash);
  last=NULL;
  cur=HashTable[hash];
  while(cur)
    {
      if((cur->size==(size_t)addr)&&(cur->type==_MDMarker))
	{
	  rec=cur;
	  if(last)
	    last->next=cur->next;
	  else
	    HashTable[hash]=cur->next;
	  break;
	}
      last=cur;
      cur=cur->next;
    }  
   if(rec)
    {
      free(rec);
    }
  else
    {
      printf(
   "MemoryDebug-ERR: free(): corrupted tree or bad ptr! (%s:%i @%p)\n",
	     file,line,addr);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
      exit(EXIT_FAILURE);
    }
  Count--;
}

int MemoryDebugUsage(void) 
{
  int a;
  unsigned int tot  = 0;
  DebugRec *rec;
  if(InitFlag) MemoryDebugInit();
  for(a=0;a<1024;a++)
    {
      rec=HashTable[a];
      while(rec)
        {
          if(rec->type!=_MDMarker)
            tot+=rec->size;
          rec=rec->next;
        }
    }
  return(tot);
}

void MemoryDebugDump(void)
{
  int a;
  int cnt=0;
  unsigned int tot  = 0;
  DebugRec *rec,*str;
  if(InitFlag) MemoryDebugInit();
  for(a=0;a<1024;a++)
    {
      rec=HashTable[a];
      while(rec)
	{
	  if(rec->type==_MDMarker)
	    {
	      str=rec+1;
	    printf(" Memory: %s:%i <%s>\n",
		   rec->file,rec->line,(char*)str);
	    }
	  else {
       tot+=rec->size;
	    printf(" Memory: %12p %12p %8x %3.1f %s:%i\n",
				  rec+1,((char*)(rec+1)+rec->size),rec->size,rec->size/1048576.0F,rec->file,rec->line);
     }
	  rec=rec->next;
	  cnt++;
	}
    }
  printf(" Memory: %d blocks expected, %d found, %d maximum allocated.\n",
			Count,cnt,MaxCount);
  printf(" Memory: current memory allocated %x bytes (%0.1f MB).\n",tot,tot/(1024.0*1024));

}

void MemoryDebugHashAdd(DebugRec *rec);
DebugRec *MemoryDebugHashRemove(void *ptr);

void MemoryDebugHashAdd(DebugRec *rec)
{
  int hash;

  hash=(int)rec;
  hash=HASH(hash);
  rec->next=HashTable[hash];
  HashTable[hash]=rec;
}

DebugRec *MemoryDebugHashRemove(void *ptr)
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


void *MemoryDebugMalloc(size_t size,const char *file,int line,int type)
{
  DebugRec *rec;

  if(InitFlag) MemoryDebugInit();
  rec=(DebugRec*)malloc(sizeof(DebugRec)+size);
  if(!rec) {
    if(!size) 
      return(NULL);
    else {
      printf("MemoryDebugMalloc-ERR: alloc failed.\n");
      DieOutOfMemory();
    }
  }
  strcpy(rec->file,file);
  rec->line=line;
  rec->size=size;
  rec->type=type;
  MemoryDebugHashAdd(rec);
  rec++;
  Count++;
  if(MaxCount<Count) MaxCount=Count;
  return((void *) rec);
}

void *MemoryDebugCalloc(size_t num,size_t size,const char *file,int line,int type)
{
  DebugRec *rec;

  if(InitFlag) MemoryDebugInit();
  rec=(DebugRec*)calloc(1,sizeof(DebugRec)+size*num);
  if(!rec)
    return(NULL);
  strcpy(rec->file,file);
  rec->line=line;
  rec->size=size;
  rec->type=type;
  MemoryDebugHashAdd(rec);
  rec++;
  Count++;
  if(MaxCount<Count) MaxCount=Count;
  return((void *) rec);
}

void *MemoryDebugRealloc(void *ptr,size_t size,const char *file,
			 int line,int type)
{
  DebugRec *rec;

  if(InitFlag) MemoryDebugInit();
  if((!ptr)&&(!size))
    {
      printf(
	     "MemoryDebug-ERR: realloc given (NULL,zero) (%s:%i)\n",
	     file,line);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
      exit(EXIT_FAILURE);
    }
  if(!ptr)
      return(MemoryDebugMalloc(size,file,line,type));
  else if(!size)
    {
      MemoryDebugFree(ptr,file,line,type);
      return(NULL);
    }
  else
    {
      rec=MemoryDebugHashRemove(ptr);
      if(!rec)
		  {
			 printf(
					  "MemoryDebug-ERR: realloc() corrupted tree or bad ptr! (%s:%i @%p)\n",
					  file,line,ptr);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
			 exit(EXIT_FAILURE);
		  }	
      else
		  {
			 if(rec->type!=type)
				{
				  printf("MemoryDebug-ERR: ptr %p is of wrong type: %i!=%i (%s:%i)\n",
							ptr,rec->type,type,file,line);
#ifdef GDB_ENTRY
              MemoryDebugDump();
              abort();
#endif
              exit(EXIT_FAILURE);
            }
			 rec=(DebugRec*)realloc(rec,size+sizeof(DebugRec));
			 if(!rec)
				{
				  printf("MemoryDebug-ERR: realloc() failed reallocation! (%s:%i)\n",
							file,line);
#ifdef GDB_ENTRY
              MemoryDebugDump();
              abort();
#endif
              DieOutOfMemory();
				}
			 else
				{
				  MemoryDebugHashAdd(rec);
				  rec->size=size;
				  rec++;
				  return((void*)rec);
				}
		  }
    }
  return(ptr);
}



void *MemoryDebugReallocForSure(void *ptr,size_t size,const char *file,
			 int line,int type)
{
  DebugRec *rec,*new_rec;

  if(InitFlag) MemoryDebugInit();
  if((!ptr)&&(!size))
    {
      printf(
	     "MemoryDebug-ERR: realloc given (NULL,zero) (%s:%i)\n",
	     file,line);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
      exit(EXIT_FAILURE);
    }
  if(!ptr)
      return(MemoryDebugMalloc(size,file,line,type));
  else if(!size)
    {
      MemoryDebugFree(ptr,file,line,type);
      return(NULL);
    }
  else
    {
      rec=MemoryDebugHashRemove(ptr);
      if(!rec)
		  {
			 printf(
					  "MemoryDebug-ERR: realloc() corrupted tree or bad ptr! (%s:%i @%p)\n",
					  file,line,ptr);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
			 exit(EXIT_FAILURE);
		  }	
      else
		  {
			 if(rec->type!=type)
				{
				  printf("MemoryDebug-ERR: ptr %p is of wrong type: %i!=%i (%s:%i)\n",
							ptr,rec->type,type,file,line);
#ifdef GDB_ENTRY
              MemoryDebugDump();
              abort();
#endif
              exit(EXIT_FAILURE);
            }

          new_rec=malloc(size+sizeof(DebugRec));
          if(new_rec)
            memcpy(new_rec,rec,size+sizeof(DebugRec));
          free(rec);
          rec=new_rec;
			 if(!rec)
				{
				  printf("MemoryDebug-ERR: realloc() failed reallocation! (%s:%i)\n",
							file,line);
#ifdef GDB_ENTRY
              MemoryDebugDump();
              abort();
#endif
              DieOutOfMemory();
				}
			 else
				{
				  MemoryDebugHashAdd(rec);
				  rec->size=size;
				  rec++;
				  return((void*)rec);
				}
		  }
    }
  return(ptr);
}

void *MemoryDebugReallocForSureSafe(void *ptr,size_t size,size_t old_size,
                                    const char *file,int line,int type)
{
  DebugRec *rec,*new_rec;

  if(InitFlag) MemoryDebugInit();
  if((!ptr)&&(!size)) {
    printf(
           "MemoryDebug-ERR: realloc given (NULL,zero) (%s:%i)\n",
           file,line);
#ifdef GDB_ENTRY
    MemoryDebugDump();
    abort();
#endif
    exit(EXIT_FAILURE);
  }
  if(!ptr)
    return(MemoryDebugMalloc(size,file,line,type));
  else if(!size) {
    MemoryDebugFree(ptr,file,line,type);
    return(NULL);
  }  else    {
    rec=MemoryDebugHashRemove(ptr);
    if(!rec) {
      printf(
             "MemoryDebug-ERR: realloc() corrupted tree or bad ptr! (%s:%i @%p)\n",
             file,line,ptr);
#ifdef GDB_ENTRY
      MemoryDebugDump();
      abort();
#endif
      exit(EXIT_FAILURE);
    } else {
      if(rec->type!=type) {
        printf("MemoryDebug-ERR: ptr %p is of wrong type: %i!=%i (%s:%i)\n",
               ptr,rec->type,type,file,line);
#ifdef GDB_ENTRY
        MemoryDebugDump();
        abort();
#endif
        exit(EXIT_FAILURE);
      }
      if(old_size>size) {
        new_rec=malloc(size+sizeof(DebugRec));
        if(new_rec)
          memcpy(new_rec,rec,size+sizeof(DebugRec));
        free(rec);
        rec=new_rec;
      } else {
        rec=realloc(rec,size+sizeof(DebugRec));
      }
      if(!rec) {
        printf("MemoryDebug-ERR: realloc() failed reallocation! (%s:%i)\n",
               file,line);
#ifdef GDB_ENTRY
        MemoryDebugDump();
        abort();
#endif
        DieOutOfMemory();
      } else {
        MemoryDebugHashAdd(rec);
        rec->size=size;
        rec++;
        return((void*)rec);
      }
    }
  }
  return(ptr);
}


void MemoryDebugQuietFree(void *ptr,int type)
{
  DebugRec *rec;

  if(InitFlag) MemoryDebugInit();
  if(!ptr)
    {
      printf("MemoryDebug-ERR: MemoryDebugQuietFree() given NULL pointer\n");
    }
  rec=MemoryDebugHashRemove(ptr);
  if(rec)
    {
      if(rec->type!=type)
		  {
			 printf("MemoryDebug-ERR: ptr %p is of wrong type: %i!=%i (allocated %s:%i)\n",
					  ptr,rec->type,type,rec->file,rec->line);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
			 exit(EXIT_FAILURE);
		  }
      free(rec);
      Count--;
    }
  else
    free(ptr);
}

void MemoryDebugFree(void *ptr,const char*file,int line,int type)
{
  DebugRec *rec;
  
  if(InitFlag) MemoryDebugInit();
  if(!ptr)
    {
      printf("MemoryDebug-ERR: free() called with NULL pointer (%s:%i)\n",
				 file,line);
#ifdef GDB_ENTRY
  MemoryDebugDump();
  abort();
#endif
      exit(EXIT_FAILURE);
    }
  rec=MemoryDebugHashRemove(ptr);
  if(rec)
    {
      if(rec->type!=type)
		  {
			 printf("MemoryDebug-ERR: ptr %p is of wrong type: %i!=%i (%s:%i)\n",
					  ptr,rec->type,type,file,line);
#ifdef GDB_ENTRY
          MemoryDebugDump();
          abort();
#endif
			 exit(EXIT_FAILURE);
		  }
      free(rec);
    }
  else
    {
		
      printf(
				 "MemoryDebug-ERR: free(): corrupted tree or bad ptr! (%s:%i @%p)\n",
				 file,line,ptr);
#ifdef GDB_ENTRY
      MemoryDebugDump();
      abort();
#endif
      exit(EXIT_FAILURE);
    }
  Count--;
}

#ifdef __cplusplus
}

#undef new

void *operator new(size_t size, const char *file,int line)
{
  if(InitFlag) MemoryDebugInit();
  return(MemoryDebugMalloc(size,file,line,_MDObject));
}

void operator delete(void *ptr)
{
  if(InitFlag) MemoryDebugInit();
  MemoryDebugQuietFree(ptr,_MDObject);
}
#endif


#endif










