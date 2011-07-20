

/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* Copyright (c) Schrodinger, LLC. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* --------------------------------------------------\-----------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#include"os_python.h"
#include "os_gl.h"

#include "Base.h"
#include "PyMOLGlobals.h"
#include "Texture.h"
#include "OOMac.h"

#include "OVContext.h"
#include "OVOneToOne.h"
#include "OVHeapArray.h"

#include "Setting.h"
#include "Character.h"
#include "Util.h"

typedef struct  {
  int id,dim;
} texture_info;

struct _CTexture {
  OVOneToOne *ch2tex;
  texture_info *info_list;
  int next_slot;
  int max_active;
};

int TextureInit(PyMOLGlobals * G)
{
  OOAlloc(G, CTexture);

  I->max_active = 2500;

  I->next_slot = 0;
  I->ch2tex = OVOneToOne_New(G->Context->heap);
  I->info_list = OVHeapArray_CALLOC(G->Context->heap, texture_info, I->max_active);

  G->Texture = I;
  return (I && I->ch2tex && I->info_list);
}

int TextureGetFromChar(PyMOLGlobals * G, int char_id, float *extent)
{
  OVreturn_word result;
  CTexture *I = G->Texture;
  int is_new = false;
  int tex_dim = 16;

  if(G->HaveGUI && G->ValidContext) {
    if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->ch2tex, char_id))) {
      if(glIsTexture(result.word))
        return result.word;
      else {
        OVOneToOne_DelReverse(I->ch2tex, result.word);
      }
    }
    {
      unsigned char *buffer = CharacterGetPixmapBuffer(G, char_id);
      if(buffer) {
        int w = CharacterGetWidth(G, char_id);
        int h = CharacterGetHeight(G, char_id);
        GLuint texture_id = 0;
        unsigned char *temp_buffer, scratch[32][32][4];
        while((tex_dim < w) || (tex_dim < h)) {
          tex_dim = (tex_dim << 1);
        }
        if(tex_dim < 33)
          temp_buffer = (unsigned char *) scratch;
        else
          temp_buffer = Alloc(unsigned char, tex_dim * tex_dim * 4);

        {
          int a, b;
          unsigned char *p = buffer, *q;
          UtilZeroMem(temp_buffer, tex_dim * tex_dim * 4);
          for(b = 0; b < h; b++) {
            for(a = 0; a < w; a++) {
              q = temp_buffer + (4 * tex_dim * b) + 4 * a;
              *(q++) = *(p++);
              *(q++) = *(p++);
              *(q++) = *(p++);
              *(q++) = *(p++);
            }
          }
          extent[0] = w / (float) tex_dim;
          extent[1] = h / (float) tex_dim;
        }

        if(!I->info_list[I->next_slot].id) {
          glGenTextures(1, &texture_id);
          is_new = true;
          I->info_list[I->next_slot].id = texture_id;
        } else {
          texture_id = I->info_list[I->next_slot].id;
          OVOneToOne_DelReverse(I->ch2tex, texture_id);
        }
        I->next_slot++;
        if(I->next_slot >= I->max_active)
          I->next_slot = 0;

        if(texture_id && OVreturn_IS_OK(OVOneToOne_Set(I->ch2tex, char_id, texture_id))) {

          glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
          glBindTexture(GL_TEXTURE_2D, texture_id);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
          if(is_new || (tex_dim != I->info_list[I->next_slot].dim)) {
	      I->info_list[I->next_slot].dim = tex_dim;
	      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
			   tex_dim, tex_dim, 0, GL_RGBA, GL_UNSIGNED_BYTE, temp_buffer);
          } else {
	      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0,
			      tex_dim, tex_dim, GL_RGBA, GL_UNSIGNED_BYTE, temp_buffer);
          }
        }
        if(temp_buffer != (unsigned char *) scratch)
          FreeP(temp_buffer);
        return texture_id;
      }
    }
  }
  return 0;
}

void TextureFree(PyMOLGlobals * G)
{
  CTexture *I = G->Texture;
  /* TODO -- free all the resident textures */

  OVOneToOne_DEL_AUTO_NULL(I->ch2tex);
  OVHeapArray_FREE_AUTO_NULL(I->info_list);

  OOFreeP(I);
}
