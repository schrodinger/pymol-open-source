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

#include"Base.h"
#include"Character.h"
#include"Pixmap.h"
#include"Util.h"
#include"MemoryDebug.h"
#include"OOMac.h"
#include"Vector.h"

#define HASH_MASK 0x4FFF

static unsigned int get_hash(CharFngrprnt *fprnt)
{
  register unsigned int result = 0;
  register unsigned short int *data = fprnt->u.d.data;
  result = (data[0]<<1) + data[1] ;
  result =  (result<<4) + data[2] ;
  result = ((result<<7) + data[3]) + (result>>16);
  result = ((result<<10) + data[4]) + (result>>16);
  result = ((result<<13) + data[5]) + (result>>16);
  result = result                  + (result>>16);

  return (HASH_MASK&result);
}

static int equal_fprnt(CharFngrprnt *f1, CharFngrprnt *f2)
{
  register unsigned short int *data1 = f1->u.d.data;
  register unsigned short int *data2 = f2->u.d.data;

  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  return 1;
}

CCharacter Character; /* global singleton */

int CharacterFind(CharFngrprnt *fprnt)
{
  CCharacter *I = &Character;
  unsigned int hash_code = get_hash(fprnt);
  int id = I->Hash[hash_code];
  while(id) {
    if(equal_fprnt(fprnt,&I->Char[id].Fngrprnt)) {

      /* pop character to top of retention list 
         (is this worth the effort?) */
      CharRec *rec = I->Char + id;
      int prev=rec->Prev;
      int next=rec->Next;
      
      if(prev&&next) { /* only act if character is in middle of list */
        I->Char[prev].Next = next;
        I->Char[next].Prev = prev;
        
        prev = I->NewestUsed;
        I->NewestUsed = id;
        I->Char[prev].Next = id;
        rec->Prev = prev;
        rec->Next = 0;
      }
      return id;
    } else
      id = I->Char[id].HashNext;
  }

  return 0;
}

int CharacterNewFromBitmap(int width, int height,
                           unsigned char *bitmap,
                           CharFngrprnt *fprnt)
{
  CCharacter *I = &Character;
  int id=CharacterGetNew();
  if((id>0)&&(id<=I->MaxAlloc)) {
    CharRec *rec = I->Char + id;
    PixmapInitFromBitmap(&rec->Pixmap,width,height,bitmap,
                         fprnt->u.i.color);    
    rec->Width = width;
    rec->Height = height;

    { /* add this character to the hash */
      int hash_code = get_hash(fprnt);
      int cur_entry;
      rec->Fngrprnt = *(fprnt);
      rec->Fngrprnt.hash_code = hash_code;
      cur_entry = I->Hash[hash_code];
      if(cur_entry) {
        I->Char[cur_entry].HashPrev = id;
      }
      I->Char[id].HashNext = I->Hash[hash_code];
      I->Hash[hash_code]=id;
    }
  }
  return id;
}

int CharacterGetWidth(int id)
{
  CCharacter *I = &Character;
  if((id>0)&&(id<=I->MaxAlloc)) {
    return I->Char[id].Width;
  }
  return 0;
}
const float _inv255 = 1.0F/255.0F;

float CharacterInterpolate(int id,float *v)
{
  CCharacter *I = &Character;
  int x = (int)v[0];
  int y = (int)v[1];
  unsigned char *src;

  if((id>0)&&(id<=I->MaxAlloc)) {
    CPixmap *pm = &I->Char[id].Pixmap;
    if(pm) {

      if(x<0) x=0;
      else if(x>=pm->width) x=pm->width-1; /* clamp */
      if(y<0) y=0;
      else if(y>=pm->height) y=pm->height-1;

      src = pm->buffer+((pm->width<<2)*y)+(x<<2);
      v[0] = *(src++) * _inv255;
      v[1] = *(src++) * _inv255;
      v[2] = *(src++) * _inv255;
      return (255-*(src++)) * _inv255;
    } else {
      zero3f(v);
      return 1.0F;
    }
  }
  return 1.0F;
}

int CharacterGetHeight(int id)
{
  CCharacter *I = &Character;
  if((id>0)&&(id<=I->MaxAlloc)) {
    return I->Char[id].Height;
  }
  return 0;
}

void CharacterInit(void)
{
  CCharacter *I = &Character;
  UtilZeroMem(I,sizeof(CCharacter));
  I->MaxAlloc = 10;
  I->Char = VLACalloc(CharRec,I->MaxAlloc+1);
  {
    int a;
    for(a=2;a<=I->MaxAlloc;a++)
      I->Char[a].Prev=a-1;
    I->LastFree = I->MaxAlloc;
  }
  I->Hash = Calloc(int,(HASH_MASK+1));
  I->TargetMaxUsage = 25000; 
}

static void CharacterAllocMore(void)
{
  CCharacter *I = &Character;
  int new_max = I->MaxAlloc * 2;
  VLACheck(I->Char,CharRec,new_max);
  {
    int a;
    I->Char[I->MaxAlloc+1].Prev = I->LastFree;
    for(a=I->MaxAlloc+2;a<=new_max;a++)
      I->Char[a].Prev = a-1;
    I->LastFree = new_max;
    I->MaxAlloc = new_max;
  }
}
void CharacterSetRetention(int retain_all)
{
  CCharacter *I = &Character;
  I->RetainAll = retain_all;
}

static void CharacterPurgeOldest(void)
{
  CCharacter *I = &Character;
  int max_kill = 10;

  while(I->NUsed > I->TargetMaxUsage) {
    if(!(max_kill--)) break; /* if over, only purge a few entries at a time */

    int id = I->OldestUsed;
    
    if(id) {
      int next;
      
      /* trim from end of list */

      if((next = I->Char[id].Next)) {
        I->Char[next].Prev = 0;
        I->OldestUsed = next;
      }
      
      { /* excise character from hash table linked list */
        int hash_code = I->Char[id].Fngrprnt.hash_code;
        int hash_prev = I->Char[id].HashPrev;;
        int hash_next = I->Char[id].HashNext ;
        
        if(hash_prev) {
          I->Char[hash_prev].HashNext = hash_next;
        } else {
          I->Hash[hash_code] = hash_next;
        }
        if(hash_next) {
          I->Char[hash_next].HashPrev = hash_prev;
        }

      }

      /* free and reinitialize */

      PixmapPurge(&I->Char[id].Pixmap);
      UtilZeroMem(I->Char+id,sizeof(CharRec));

      /* add to free chain */

      I->Char[id].Prev = I->LastFree;
      I->LastFree = id;
      I->NUsed--;
    }
  }
}

int CharacterGetNew(void)
{
  CCharacter *I = &Character;
  int result = 0;
  if(!I->LastFree)
    CharacterAllocMore();
  if(I->LastFree) {

    /* remove from free chain */
    result = I->LastFree;
    I->LastFree = I->Char[result].Prev;

    /* backwards-link (for continuous GC) */

    if(I->NewestUsed) {
      I->Char[I->NewestUsed].Next = result; /* double-link list */
    } else {
      I->OldestUsed = result;
    }

    /* forwards-link */

    I->Char[result].Prev = I->NewestUsed;
    I->NewestUsed = result;

    I->NUsed++;
        
    if(!I->RetainAll)
      CharacterPurgeOldest();
  }

  return result;
}

void CharacterFree(void)
{
  CCharacter *I = &Character;
  {
    int a;
    a = I->NewestUsed;
    while(a) {
      PixmapPurge(&I->Char[a].Pixmap);
      a = I->Char[a].Prev;
    }
  }
  FreeP(I->Hash);
  VLAFreeP(I->Char);
}

