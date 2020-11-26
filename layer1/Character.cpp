
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
#include"CGO.h"
#include <glm/gtc/type_ptr.hpp>

#define HASH_MASK 0x2FFF

static unsigned int get_hash(CharFngrprnt * fprnt)
{
  unsigned int result = 0;
  unsigned short int *data = fprnt->u.d.data;
  result = (data[0] << 1) + data[1];
  result = ((result << 4) + data[2]);
  result = ((result << 7) + data[3]) + (result >> 16);
  result = ((result << 10) + data[4]) + (result >> 16);
  result = ((result << 13) + data[5]) + (result >> 16);
  result = ((result << 15) + data[6]) + (result >> 16);
  result = ((result << 15) + data[7]) + (result >> 16);
  result = ((result << 15) + data[8]) + (result >> 16);
  result = ((result << 1) + data[9]) + (result >> 16);
  return (HASH_MASK & result);
}

static int equal_fprnt(CharFngrprnt * f1, CharFngrprnt * f2)
{
  unsigned short int *data1 = f1->u.d.data;
  unsigned short int *data2 = f2->u.d.data;

  /* must compare all fields in fingerprint */

  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  if(*(data1++) != *(data2++))
    return 0;
  return 1;
}

int CharacterFind(PyMOLGlobals * G, CharFngrprnt * fprnt)
{
  CCharacter *I = G->Character;
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
    if(equal_fprnt(fprnt, &I->Char[id].Fngrprnt)) {

      /* pop character to top of retention list 
         (is this worth the effort?) */
      CharRec *rec = I->Char + id;
      int prev = rec->Prev;
      int next = rec->Next;

      if(prev && next) {        /* only act if character is in middle of list */
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

unsigned char *CharacterGetPixmapBuffer(PyMOLGlobals * G, int id)
{
  CCharacter *I = G->Character;
  if(id) {
    CharRec *rec = I->Char + id;
    return rec->Pixmap.buffer;
  }
  return NULL;
}

int CharacterNewFromBitmap(PyMOLGlobals * G, int width, int height,
                           unsigned char *bitmap,
                           float x_orig, float y_orig, float advance,
                           CharFngrprnt * fprnt, int sampling)
{
  CCharacter *I = G->Character;
  int id = CharacterGetNew(G);
  if((id > 0) && (id <= I->MaxAlloc)) {
    CharRec *rec = I->Char + id;
    PixmapInitFromBitmap(G, &rec->Pixmap, width, height, bitmap,
                         fprnt->u.i.color, sampling);
    rec->Width = width * sampling;
    rec->Height = height * sampling;
    rec->XOrig = x_orig * sampling;
    rec->YOrig = y_orig * sampling;
    rec->Advance = advance * sampling;
    {                           /* add this character to the hash */
      int hash_code = get_hash(fprnt);
      int cur_entry;
      rec->Fngrprnt = *(fprnt);
      rec->Fngrprnt.hash_code = hash_code;
      cur_entry = I->Hash[hash_code];
      if(cur_entry) {
        I->Char[cur_entry].HashPrev = id;
      }
      I->Char[id].HashNext = I->Hash[hash_code];
      I->Hash[hash_code] = id;
    }
  }
  return id;
}

int CharacterNewFromBytemap(PyMOLGlobals * G, int width, int height,
                            int pitch, unsigned char *bytemap,
                            float x_orig, float y_orig, float advance,
                            CharFngrprnt * fprnt)
{
  CCharacter *I = G->Character;
  int id = CharacterGetNew(G);
  if((id > 0) && (id <= I->MaxAlloc)) {
    CharRec *rec = I->Char + id;
    PixmapInitFromBytemap(G, &rec->Pixmap, width, height, pitch, bytemap,
                          fprnt->u.i.color, fprnt->u.i.outline_color, fprnt->u.i.flat);
    rec->Width = width;
    rec->Height = height;
    rec->XOrig = x_orig;
    rec->YOrig = y_orig;
    rec->Advance = advance;
    {                           /* add this character to the hash */
      int hash_code = get_hash(fprnt);
      int cur_entry;
      rec->Fngrprnt = *(fprnt);
      rec->Fngrprnt.hash_code = hash_code;
      cur_entry = I->Hash[hash_code];
      if(cur_entry) {
        I->Char[cur_entry].HashPrev = id;
      }
      I->Char[id].HashNext = I->Hash[hash_code];
      I->Hash[hash_code] = id;
    }
  }
  return id;
}

float CharacterGetAdvance(PyMOLGlobals * G, int sampling, int id)
{
  CCharacter *I = G->Character;
  CharRec *rec = I->Char + id;
  return rec->Advance / sampling;
}

void CharacterRenderOpenGLPrime(PyMOLGlobals * G, const RenderInfo * info)
{
  if(G->HaveGUI && G->ValidContext) {
    if ((info && !info->use_shaders) || (!info && !SettingGetGlobal_b(G, cSetting_use_shaders))){
      glEnable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    }
    /*    glEnable(GL_BLEND);
       glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); */
  }
}

void CharacterRenderOpenGLDone(PyMOLGlobals * G, const RenderInfo * info)
{
  if(G->HaveGUI && G->ValidContext) {
    if ((info && !info->use_shaders) || (!info && !SettingGetGlobal_b(G, cSetting_use_shaders))){
      glDisable(GL_TEXTURE_2D);
    }
    /*    glDisable(GL_BLEND); */
  }
}

short CharacterRenderOpenGL(PyMOLGlobals * G, const RenderInfo * info, int id, short isworldlabel, short relativeMode SHADERCGOARG)

/* need orientation matrix */
{
  CCharacter *I = G->Character;
  CharRec *rec = I->Char + id;
  short success = 1;
  int texture_id = TextureGetFromChar(G, id, rec->extent);
  float sampling = 1.0F;

  if(G->HaveGUI && G->ValidContext && texture_id) {
    if(info)
      sampling = (float) info->sampling;
    if(texture_id) {
    /*    if(glIsTexture(texture_id)) -- BAD -- impacts performance */
      float *v, v0[3];
      float v1[3];
      if (!shaderCGO){
        glBindTexture(GL_TEXTURE_2D, TextGetIsPicking(G) ? 0 : texture_id);
      }
      v = TextGetPos(G);
      copy3f(v, v0);
      v0[0] -= rec->XOrig / sampling;
      v0[1] -= rec->YOrig / sampling;
      copy3f(v0, v1);
      v1[0] += rec->Width / sampling;
      v1[1] += rec->Height / sampling;
      /*      glColor4f(0.5F,0.5F,0.5F,1.0F); */
      if (shaderCGO){
          float *worldPos = TextGetWorldPos(G);
	  if (isworldlabel){
	    float *targetPos = TextGetTargetPos(G);
	    float *screenWorldOffset = TextGetScreenWorldOffset(G);
            shaderCGO->add<cgo::draw::label>(glm::make_vec3(worldPos),
                                             glm::make_vec3(screenWorldOffset),
                                             glm::make_vec3(v0),
                                             glm::make_vec3(v1),
                                             glm::make_vec4(rec->extent),
                                             relativeMode,
                                             glm::make_vec3(targetPos));
	  } else {
	    CGODrawTexture(shaderCGO, texture_id, worldPos, v0, v1, rec->extent);
	  }
      } else {
#ifndef PURE_OPENGL_ES_2
	  glBegin(GL_QUADS);
        if (TextGetIsPicking(G)){
          unsigned char *cptr = TextGetColorUChar4uv(G);
          glColor4ubv(cptr);
          glVertex3f(v0[0], v0[1], v0[2]);
          glVertex3f(v0[0], v1[1], v0[2]);
          glVertex3f(v1[0], v1[1], v0[2]);
          glVertex3f(v1[0], v0[1], v0[2]);
          glEnd();
        } else {
	  glTexCoord2f(rec->extent[0], rec->extent[1]);
	  glVertex3f(v0[0], v0[1], v0[2]);
	  glTexCoord2f(rec->extent[0], rec->extent[3]);
	  glVertex3f(v0[0], v1[1], v0[2]);
	  glTexCoord2f(rec->extent[2], rec->extent[3]);
	  glVertex3f(v1[0], v1[1], v0[2]);
	  glTexCoord2f(rec->extent[2], rec->extent[1]);
	  glVertex3f(v1[0], v0[1], v0[2]);
	  glEnd();
        }
#endif
      }
    }
     TextAdvance(G, rec->Advance / sampling);
  } else {
    if (!texture_id)
      success = 0;
  }
  return success;
}

int CharacterGetWidth(PyMOLGlobals * G, int id)
{
  CCharacter *I = G->Character;
  if((id > 0) && (id <= I->MaxAlloc)) {
    return I->Char[id].Width;
  }
  return 0;
}

const float _inv255 = 1.0F / 255.0F;

const unsigned char zerouc[4] = { 0, 0, 0, 0 };


/* CharacterInterpolate: This function implements bilinear interpolation
   on looking up the pixel value in the texture.
*/
float CharacterInterpolate(PyMOLGlobals * G, int id, float *v)
{
  CCharacter *I = G->Character;
  int x = (int) v[0];
  int y = (int) v[1];

  if((id > 0) && (id <= I->MaxAlloc)) {
    CPixmap *pm = &I->Char[id].Pixmap;
    if(pm) {
      unsigned char *srcx0, *srcx1, *srcy0, *srcy1;
      int x1 = x + 1, y1 = y + 1;
      float xdiff = v[0] - x, ydiff = v[1] - y;
      float xdiff1 = (1.f-xdiff), ydiff1 = (1.f-ydiff);
      float interp0[4], interp1[4]; // interpolated in x for y0 and y1

      if (x < 0 || x > (pm->width - 1))
	srcx0 = (unsigned char *)zerouc;
      else
	srcx0 = pm->buffer + ((pm->width << 2) * y) + (x << 2);
      if (x1 < 0 || x1 > (pm->width - 1))
	srcx1 = (unsigned char *)zerouc;
      else
	srcx1 = pm->buffer + ((pm->width << 2) * y) + (x1 << 2);

      if (y1 < 0 || y1 > (pm->height - 1))
	srcy0 = (unsigned char *)zerouc;
      else
	srcy0 = pm->buffer + ((pm->width << 2) * y1) + (x << 2);

      if (x1 < 0 || x1 > (pm->width - 1) || y1 < 0 || y1 > (pm->height - 1))
	srcy1 = (unsigned char *)zerouc;
      else
	srcy1 = pm->buffer + ((pm->width << 2) * y1) + (x1 << 2);

      interp0[0] = (xdiff1 * (*(srcx0++)) + (xdiff * (*(srcx1++))));
      interp0[1] = (xdiff1 * (*(srcx0++)) + (xdiff * (*(srcx1++))));
      interp0[2] = (xdiff1 * (*(srcx0++)) + (xdiff * (*(srcx1++))));
      interp0[3] = (xdiff1 * (*(srcx0++)) + (xdiff * (*(srcx1++))));

      interp1[0] = (xdiff1 * (*(srcy0++)) + (xdiff * (*(srcy1++))));
      interp1[1] = (xdiff1 * (*(srcy0++)) + (xdiff * (*(srcy1++))));
      interp1[2] = (xdiff1 * (*(srcy0++)) + (xdiff * (*(srcy1++))));
      interp1[3] = (xdiff1 * (*(srcy0++)) + (xdiff * (*(srcy1++))));

      v[0] = ((ydiff1 * interp0[0]) + (ydiff * interp1[0])) * _inv255;
      v[1] = ((ydiff1 * interp0[1]) + (ydiff * interp1[1])) * _inv255;
      v[2] = ((ydiff1 * interp0[2]) + (ydiff * interp1[2])) * _inv255;
      return (255 - ((ydiff1 * interp0[3]) + (ydiff * interp1[3]))) * _inv255;
    } else {
      zero3f(v);
      return 1.0F;
    }
  }
  zero3f(v);
  return 1.0F;
}

int CharacterGetHeight(PyMOLGlobals * G, int id)
{
  CCharacter *I = G->Character;
  if((id > 0) && (id <= I->MaxAlloc)) {
    return I->Char[id].Height;
  }
  return 0;
}

int CharacterGetGeometry(PyMOLGlobals * G, int id,
                         int *width, int *height,
                         float *xorig, float *yorig, float *advance)
{
  CCharacter *I = G->Character;
  if((id > 0) && (id <= I->MaxAlloc)) {
    CharRec *ch = I->Char + id;
    *width = ch->Width;
    *height = ch->Height;
    *xorig = ch->XOrig;
    *yorig = ch->YOrig;
    *advance = ch->Advance;
  }
  return 0;
}

int CharacterInit(PyMOLGlobals * G)
{
  CCharacter *I = NULL;
  if((I = (G->Character = pymol::calloc<CCharacter>(1)))) {
    I->MaxAlloc = 5;
    I->Char = VLACalloc(CharRec, I->MaxAlloc + 1);
    {
      int a;
      for(a = 2; a <= I->MaxAlloc; a++)
        I->Char[a].Prev = a - 1;
      I->LastFree = I->MaxAlloc;
    }
    I->Hash = pymol::calloc<int>((HASH_MASK + 1));
    I->TargetMaxUsage = 25000;
    return 1;
  } else
    return 0;
}

static void CharacterAllocMore(PyMOLGlobals * G)
{
  CCharacter *I = G->Character;
  int new_max = I->MaxAlloc * 2;
  VLACheck(I->Char, CharRec, new_max);
  {
    int a;
    I->Char[I->MaxAlloc + 1].Prev = I->LastFree;
    for(a = I->MaxAlloc + 2; a <= new_max; a++)
      I->Char[a].Prev = a - 1;
    I->LastFree = new_max;
    I->MaxAlloc = new_max;
  }
}

void CharacterSetRetention(PyMOLGlobals * G, int retain_all)
{
  CCharacter *I = G->Character;
  I->RetainAll = retain_all;
}

static void CharacterPurgeOldest(PyMOLGlobals * G)
{
  CCharacter *I = G->Character;
  int max_kill = 10;

  while(I->NUsed > I->TargetMaxUsage) {
    if(!(max_kill--))
      break;                    /* if over, only purge a few entries at a time */
    {
      int id = I->OldestUsed;

      if(id) {
        int next;

        /* trim from end of list */

        if((next = I->Char[id].Next)) {
          I->Char[next].Prev = 0;
          I->OldestUsed = next;
        }

        {                       /* excise character from hash table linked list */
          int hash_code = I->Char[id].Fngrprnt.hash_code;
          int hash_prev = I->Char[id].HashPrev;
          int hash_next = I->Char[id].HashNext;

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
        UtilZeroMem(I->Char + id, sizeof(CharRec));

        /* add to free chain */

        I->Char[id].Prev = I->LastFree;
        I->LastFree = id;
        I->NUsed--;
      }
    }
  }
}

int CharacterGetNew(PyMOLGlobals * G)
{
  CCharacter *I = G->Character;
  int result = 0;
  if(!I->LastFree)
    CharacterAllocMore(G);
  if(I->LastFree) {

    /* remove from free chain */
    result = I->LastFree;
    I->LastFree = I->Char[result].Prev;

    /* backwards-link (for continuous GC) */

    if(I->NewestUsed) {
      I->Char[I->NewestUsed].Next = result;     /* double-link list */
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

void CharacterFree(PyMOLGlobals * G)
{
  CCharacter *I = G->Character;
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
