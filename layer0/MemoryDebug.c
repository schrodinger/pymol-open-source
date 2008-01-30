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
#ifndef OV_JENARIX

#include"MemoryDebug.h"
#include"MemoryCache.h"

#define GDB_ENTRY

void *MemoryReallocForSureSafe(void *ptr, unsigned int new_size, unsigned int old_size)
{
  if(new_size<old_size) {
    float *tmp = mmalloc(new_size);
    if(tmp && new_size && old_size) {
        memcpy(tmp, ptr, new_size);
    }
    FreeP(ptr);
    return tmp;
  } else {
    return mrealloc(ptr,new_size);
  }
}

void *MemoryReallocForSure(void *ptr, unsigned int new_size) /* unsafe -- replace with above */
{
  float *tmp = mmalloc(new_size);
  if(tmp)
    memcpy(tmp,ptr,new_size);
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

void *VLAExpand(void *ptr,ov_size rec)
{
  VLARec *vla;
  char *start,*stop;
  unsigned int soffset=0;
  vla = &(((VLARec*)ptr)[-1]);
  if(rec>=vla->size) {
    if(vla->auto_zero)
      soffset = sizeof(VLARec)+(vla->unit_size*vla->size);
    vla->size = ((unsigned int)(rec*vla->grow_factor))+1;
    if(vla->size<=rec) vla->size = rec+1;
    {
      VLARec *old_vla = vla;
      vla=(void*)mrealloc(vla,(vla->unit_size*vla->size)+sizeof(VLARec));
      while(!vla) {  /* back off on the request size until it actually fits */
        vla = old_vla;
        vla->grow_factor = (vla->grow_factor-1.0F)/2.0F + 1.0F;   
        vla->size = ((unsigned int)(rec*vla->grow_factor))+1;
        vla=(void*)mrealloc(vla,(vla->unit_size*vla->size)+sizeof(VLARec));
        if(!vla) {
          if(old_vla->grow_factor<1.001F) {
            printf("VLAExpand-ERR: realloc failed.\n");
            DieOutOfMemory();
          }
        }
      }
    }
    if(vla->auto_zero) {
      start = ((char*)vla) + soffset;
      stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
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
  if(rec>=vla->size)
	 {
		if(vla->auto_zero)
		  soffset = sizeof(VLARec)+(vla->unit_size*vla->size);
		vla->size = ((unsigned int)(rec*vla->grow_factor))+1;
        if(vla->size<=rec) vla->size = rec+1;
		vla=(void*)_MemoryCacheRealloc(G,vla,
                                       (vla->unit_size*vla->size)+sizeof(VLARec),
                                       thread_index,block_id MD_FILE_LINE_Call);
		if(!vla)
		  {
			 printf("VLAExpand-ERR: realloc failed.\n");
          DieOutOfMemory();
		  }
		if(vla->auto_zero)
		  {
			 start = ((char*)vla) + soffset;
			 stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
			 MemoryZero(start,stop);
		  }
	 }
  return((void*)&(vla[1]));
}
#endif

#ifndef _MemoryDebug_ON
void *VLAMalloc(ov_size init_size,ov_size unit_size,unsigned int grow_factor,int auto_zero)
#else
void *_VLAMalloc(const char *file,int line,ov_size init_size,
                 ov_size unit_size,unsigned int grow_factor,int auto_zero)
#endif
{
  VLARec *vla;
  char *start,*stop;
#ifndef _MemoryDebug_ON
  vla=(void*)mmalloc((init_size*unit_size)+sizeof(VLARec));
#else
  vla=MemoryDebugMalloc((init_size*unit_size)+sizeof(VLARec),file,line,_MDPointer);
#endif

  if(!vla) {
    printf("VLAMalloc-ERR: malloc failed\n");
    DieOutOfMemory();
  }
  vla->size=init_size;
  vla->unit_size=unit_size;
  vla->grow_factor=(1.0F + grow_factor*0.1F);
  vla->auto_zero=auto_zero;
  if(vla->auto_zero) {
    start = ((char*)vla)+sizeof(VLARec);
    stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
    MemoryZero(start,stop);
  }
  return((void*)&(vla[1]));
}


#ifdef _MemoryCache_ON
#ifndef _MemoryDebug_ON
void *VLACacheMalloc(PyMOLGlobals *G,unsigned int init_size,unsigned int unit_size,
                     unsigned int grow_factor,int auto_zero,int thread,int id)
#else
void *_VLACacheMalloc(PyMOLGlobals *G,const char *file,int line,
                      unsigned int init_size,unsigned int unit_size,unsigned int grow_factor,
                      int auto_zero,int thread,int id)
#endif
{
  VLARec *vla;
  char *start,*stop;

  vla=(void*)_MemoryCacheMalloc(G,(init_size*unit_size)+sizeof(VLARec),
                                thread,id 
                                MD_FILE_LINE_Nest);

  if(!vla)
	 {
		printf("VLAMalloc-ERR: realloc failed\n");
      DieOutOfMemory();
	 }
  vla->size=init_size;
  vla->unit_size=unit_size;
  vla->grow_factor=(1.0F + grow_factor*0.1F);
  vla->auto_zero=auto_zero;
  if(vla->auto_zero)
	 {
		start = ((char*)vla)+sizeof(VLARec);
		stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
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
  return(vla->size);
}

void *VLANewCopy(void *ptr)
{
  if(ptr) { /* NULL protected */
    VLARec *vla,*new_vla;
    unsigned int size;
    vla = &((VLARec*)ptr)[-1];
    size = (vla->unit_size*vla->size)+sizeof(VLARec);
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
  } else {
    return NULL;
  }
}
#ifdef _MemoryCache_ON
void *VLACacheSetSize(PyMOLGlobals *G,void *ptr,unsigned int new_size,int group_id,int block_id)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->auto_zero) {
	 soffset = sizeof(VLARec)+(vla->unit_size*vla->size);
  }
  vla->size = new_size;
  vla=(void*)_MemoryCacheRealloc(G,vla,(vla->unit_size*vla->size)+sizeof(VLARec),
                                 group_id,block_id MD_FILE_LINE_Call);
  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed.\n");
      DieOutOfMemory();
	 }
  if(vla->auto_zero)
	 {
      start = ((char*)vla)+soffset;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
		if(start<stop)
		  MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}
void *VLACacheSetSizeForSure(PyMOLGlobals *G,void *ptr,unsigned int new_size,int group_id,int block_id)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->auto_zero) {
	 soffset = sizeof(VLARec)+(vla->unit_size*vla->size);
  }
  if(new_size<vla->size) {
    vla->size = new_size;
    vla=(void*)_MemoryCacheShrinkForSure(G,vla,(vla->unit_size*vla->size)+sizeof(VLARec),
                                         group_id,block_id MD_FILE_LINE_Call);
  } else {
    vla->size = new_size;
    vla=(void*)_MemoryCacheRealloc(G,vla,(vla->unit_size*vla->size)+sizeof(VLARec),
                                   group_id,block_id MD_FILE_LINE_Call);
  }
  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed.\n");
      DieOutOfMemory();
	 }
  if(vla->auto_zero)
	 {
      start = ((char*)vla)+soffset;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
		if(start<stop)
		  MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}
#endif


void *VLASetSize(void *ptr,unsigned int new_size)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->auto_zero) {
	 soffset = sizeof(VLARec)+(vla->unit_size*vla->size);
  }
  vla->size = new_size;
  vla=(void*)mrealloc(vla,(vla->unit_size*vla->size)+sizeof(VLARec));
  if(!vla)
	 {
		printf("VLASetSize-ERR: realloc failed.\n");
      DieOutOfMemory();
	 }
  if(vla->auto_zero)
	 {
      start = ((char*)vla)+soffset;
		stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
		if(start<stop)
		  MemoryZero(start,stop);
	 }
  return((void*)&(vla[1]));
}

void *VLASetSizeForSure(void *ptr,unsigned int new_size)
{
  VLARec *vla;
  char *start=NULL;
  char *stop;
  unsigned int soffset=0;
  vla = &((VLARec*)ptr)[-1];
  if(vla->auto_zero) {
    soffset = sizeof(VLARec)+(vla->unit_size*vla->size);
  }
  if(new_size<vla->size) {
    vla=MemoryReallocForSureSafe(vla,
                                 (vla->unit_size*new_size)+sizeof(VLARec),
                                 (vla->unit_size*vla->size)+sizeof(VLARec));
    vla->size = new_size;
  } else {
    vla->size = new_size;
    vla=(void*)mrealloc(vla,(vla->unit_size*vla->size)+sizeof(VLARec));
  }
  if(!vla) {
    printf("VLASetSize-ERR: realloc failed.\n");
    DieOutOfMemory();
  }
  if(vla->auto_zero) {
    start = ((char*)vla)+soffset;
    stop = ((char*)vla)+sizeof(VLARec)+(vla->unit_size*vla->size);
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

int MemoryDebugUsage(void) 
{
  int a;
  unsigned int tot  = 0;
  DebugRec *rec;
  if(InitFlag) MemoryDebugInit();
  for(a=0;a<1024;a++) {
    rec=HashTable[a];
    while(rec) {
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
  DebugRec *rec;
  if(InitFlag) MemoryDebugInit();
  for(a=0;a<1024;a++)  {
    rec=HashTable[a];
    while(rec)	{
	  tot+=rec->size;
      printf("Memory: %12p %12p %8x %3.1f %s:%i\n",
             (void*)(rec+1),
             ((char*)(rec+1)+rec->size),(unsigned int)rec->size,
             rec->size/1048576.0F,rec->file,rec->line);
	  rec=rec->next;
	  cnt++;
	}
  }
  printf("Memory: %d blocks expected, %d found, %d maximum allocated.\n",
         Count,cnt,MaxCount);
  printf("Memory: current memory allocated %x bytes (%0.1f MB).\n",tot,tot/(1024.0*1024));

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
#endif


#endif
#endif









