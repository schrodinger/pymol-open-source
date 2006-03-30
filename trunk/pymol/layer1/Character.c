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

#include"os_gl.h"

#include"Base.h"
#include"Character.h"
#include"Pixmap.h"
#include"Util.h"
#include"MemoryDebug.h"
#include"OOMac.h"
#include"Vector.h"
#include"Text.h"
#include"Texture.h"

#define HASH_MASK 0x2FFF

static unsigned int get_hash(CharFngrprnt *fprnt)
{
  register unsigned int result = 0;
  register unsigned short int *data = fprnt->u.d.data;
  result = (data[0]<< 1) + data[1] ;
  result = ((result<< 4) + data[2]) ;
  result = ((result<< 7) + data[3]) + (result>>16);
  result = ((result<<10) + data[4]) + (result>>16);
  result = ((result<<13) + data[5]) + (result>>16);
  result = ((result<<15) + data[6]) + (result>>16);
  result = ((result<<15) + data[7]) + (result>>16);
  result = ((result<<15) + data[8]) + (result>>16);
  result = ((result<< 1) + data[9]) + (result>>16);
  return (HASH_MASK&result);
}

static int equal_fprnt(CharFngrprnt *f1, CharFngrprnt *f2)
{
  register unsigned short int *data1 = f1->u.d.data;
  register unsigned short int *data2 = f2->u.d.data;

  /* must compare all fields in fingerprint */

  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  if(*(data1++)!=*(data2++)) return 0;
  return 1;
}

int CharacterFind(PyMOLGlobals *G,CharFngrprnt *fprnt)
{
  register CCharacter *I = G->Character;
  unsigned int hash_code = get_hash(fprnt);
  int id = I->Hash[hash_code];
  /*
  printf("seeking %d %d %d %d %d %d %d\n",
         fprnt->u.i.text_id,
         fprnt->u.i.ch,
         fprnt->u.i.height,
         fprnt->u.i.color[0],
         fprnt->u.i.color[1],
         fprnt->u.i.color[2],
         fprnt->u.i.color[3]);
  */

         
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

unsigned char *CharacterGetPixmapBuffer(PyMOLGlobals *G,int id)
{
  register CCharacter *I = G->Character;
  if(id) {
    CharRec *rec = I->Char + id;
    return rec->Pixmap.buffer;
  }
  return NULL;
}
                                   

int CharacterNewFromBitmap(PyMOLGlobals *G, int width, int height,
                           unsigned char *bitmap,
                           float x_orig, float y_orig, float advance,
                           CharFngrprnt *fprnt,int sampling)
{
  register CCharacter *I = G->Character;
  int id=CharacterGetNew(G);
  if((id>0)&&(id<=I->MaxAlloc)) {
    CharRec *rec = I->Char + id;
    PixmapInitFromBitmap(G,&rec->Pixmap,width,height,bitmap,
                         fprnt->u.i.color,sampling);    
    rec->Width = width * sampling;
    rec->Height = height * sampling;
    rec->XOrig = x_orig * sampling;
    rec->YOrig = y_orig * sampling;
    rec->Advance = advance * sampling;
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

int CharacterNewFromBytemap(PyMOLGlobals *G, int width, int height,
                             int pitch, unsigned char *bytemap,
                            float x_orig, float y_orig, float advance,
                            CharFngrprnt *fprnt)
{
  register CCharacter *I = G->Character;
  int id=CharacterGetNew(G);
  if((id>0)&&(id<=I->MaxAlloc)) {
    CharRec *rec = I->Char + id;
    PixmapInitFromBytemap(G,&rec->Pixmap,width,height,pitch,bytemap,
                         fprnt->u.i.color,fprnt->u.i.outline_color,fprnt->u.i.flat);    
    rec->Width = width;
    rec->Height = height;
    rec->XOrig = x_orig;
    rec->YOrig = y_orig;
    rec->Advance = advance;
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

float CharacterGetAdvance(PyMOLGlobals *G,int sampling, int id)
{
  register CCharacter *I = G->Character;
  CharRec *rec = I->Char + id;
  return rec->Advance/sampling;
}

void CharacterRenderOpenGL(PyMOLGlobals *G,RenderInfo *info,int id)
/* need orientation matrix */
{
  register CCharacter *I = G->Character;
  CharRec *rec = I->Char + id;

  int texture_id = TextureGetFromChar(G,id,rec->extent);
  float sampling = 1.0F;
  if(G->HaveGUI &&  G->ValidContext && texture_id) {
    if(info)
      sampling = (float)info->sampling;
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    if(glIsTexture(texture_id)) {
      float *v,v0[3];
      float v1[3];
      glBindTexture(GL_TEXTURE_2D, texture_id);
      v = TextGetPos(G);
      copy3f(v,v0);
      v0[0]-=rec->XOrig/sampling;
      v0[1]-=rec->YOrig/sampling;
      copy3f(v0,v1);
      v1[0] += rec->Width/sampling;
      v1[1] += rec->Height/sampling;
      /*      glColor4f(0.5F,0.5F,0.5F,1.0F);*/
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0F, 0.0F); glVertex3f(v0[0],v0[1],v0[2]);
      glTexCoord2f(0.0F, rec->extent[1]); glVertex3f(v0[0],v1[1],v0[2]);
      glTexCoord2f(rec->extent[0], rec->extent[1]); glVertex3f(v1[0],v1[1],v0[2]);
      glTexCoord2f(rec->extent[0], 0.0F); glVertex3f(v1[0],v0[1],v0[2]);
      glEnd();
      TextAdvance(G,rec->Advance/sampling);
    }
    glDisable(GL_TEXTURE_2D);
  }
}

int CharacterGetWidth(PyMOLGlobals *G,int id)
{
  register CCharacter *I = G->Character;
  if((id>0)&&(id<=I->MaxAlloc)) {
    return I->Char[id].Width;
  }
  return 0;
}
const float _inv255 = 1.0F/255.0F;

float CharacterInterpolate(PyMOLGlobals *G,int id,float *v)
{
  register CCharacter *I = G->Character;
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

int CharacterGetHeight(PyMOLGlobals *G,int id)
{
  register CCharacter *I = G->Character;
  if((id>0)&&(id<=I->MaxAlloc)) {
    return I->Char[id].Height;
  }
  return 0;
}

int CharacterGetGeometry(PyMOLGlobals *G,int id,
                         int *width, int *height, 
                         float *xorig, float *yorig, float *advance)
{
  register CCharacter *I = G->Character;
  if((id>0)&&(id<=I->MaxAlloc)) {
    CharRec *ch = I->Char + id;
    *width = ch->Width;
    *height = ch->Height;
    *xorig = ch->XOrig;
    *yorig = ch->YOrig;
    *advance = ch->Advance;
  }
  return 0;
}

int CharacterInit(PyMOLGlobals *G)
{
  register CCharacter *I=NULL;
  if( (I=(G->Character=Calloc(CCharacter,1)))) {
    I->MaxAlloc = 5;
    I->Char = VLACalloc(CharRec,I->MaxAlloc+1);
    {
      int a;
      for(a=2;a<=I->MaxAlloc;a++)
        I->Char[a].Prev=a-1;
      I->LastFree = I->MaxAlloc;
    }
    I->Hash = Calloc(int,(HASH_MASK+1));
    I->TargetMaxUsage = 25000; 
    return 1;
  } else 
    return 0;
}

static void CharacterAllocMore(PyMOLGlobals *G)
{
  register CCharacter *I = G->Character;
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
void CharacterSetRetention(PyMOLGlobals *G,int retain_all)
{
  register CCharacter *I = G->Character;
  I->RetainAll = retain_all;
}

static void CharacterPurgeOldest(PyMOLGlobals *G)
{
  register CCharacter *I = G->Character;
  int max_kill = 10;

  while(I->NUsed > I->TargetMaxUsage) {
    if(!(max_kill--)) break; /* if over, only purge a few entries at a time */
    {
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
        int hash_prev = I->Char[id].HashPrev;
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
}

int CharacterGetNew(PyMOLGlobals *G)
{
  register CCharacter *I = G->Character;
  int result = 0;
  if(!I->LastFree)
    CharacterAllocMore(G);
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
      CharacterPurgeOldest(G);
  }

  return result;
}

void CharacterFree(PyMOLGlobals *G)
{
  register CCharacter *I = G->Character;
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
  FreeP(G->Character);
}

