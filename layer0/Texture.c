/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
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
#include "os_gl.h"

#include "PyMOLGlobals.h"
#include "Texture.h"
#include "OOMac.h"

#include "OVContext.h"
#include "OVOneToOne.h"
#include "OVHeapArray.h"

#include "Setting.h"
#include "Character.h"

struct _CTexture {
  OVOneToOne *ch2tex;
  int *ch_list;
  int next_slot;
  int n_active;
};

int TextureInit(PyMOLGlobals *G)
{
  OOAlloc(G,CTexture);  

  I->n_active = 0;
  I->next_slot = 1;
  I->ch2tex = OVOneToOne_New(G->Context->heap);
  I->ch_list = OVHeapArray_MALLOC(G->Context->heap, int, 1);
  
  G->Texture = I;
  return (I && I->ch2tex && I->ch_list);
}


int TextureGetFromChar(PyMOLGlobals *G, int char_id,float *extent)
{
  OVreturn_word result;
  CTexture *I=G->Texture;

  if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->ch2tex,char_id))) {
    return result.word;
  } else {
    unsigned char *buffer = CharacterGetPixmapBuffer(G,char_id);
    if(buffer) {
      int w = CharacterGetWidth(G,char_id);
      int h = CharacterGetHeight(G,char_id);
      unsigned char temp_buffer[64][64][4];
      {
        int a,b,c;
        unsigned char *p = buffer;
          
        UtilZeroMem(temp_buffer,64*64*4);
        for(b=0;b<h;b++)
          for(a=0;a<w;a++)
            {
              temp_buffer[b][a][0] = *(p++);
              temp_buffer[b][a][1] = *(p++);
              temp_buffer[b][a][2] = *(p++);
              temp_buffer[b][a][3] = *(p++);
              /*              printf("%3d %3d %3d %3d\n",
                     temp_buffer[b][a][0],
                     temp_buffer[b][a][1],
                     temp_buffer[b][a][2],
                     temp_buffer[b][a][3]);*/
            }
        extent[0]=w/64.0F;
        extent[1]=h/64.0F;
      }
      GLuint texture_id;
      glGenTextures(1,&texture_id);
      if(texture_id &&
         OVreturn_IS_OK(OVOneToOne_Set(I->ch2tex,char_id,texture_id))) {
        
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glBindTexture(GL_TEXTURE_2D, texture_id);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
                     64,
                     64,
                     0,
                     GL_RGBA,
                     GL_UNSIGNED_BYTE,
                     temp_buffer);
      }
      return texture_id;
    }
  }
  return 0;
}

void TextureFree(PyMOLGlobals *G)
{
  CTexture *I=G->Texture;
  /* TODO -- free all the resident textures */

  OVOneToOne_DEL_AUTO_NULL(I->ch2tex);
  OVHeapArray_FREE_AUTO_NULL(I->ch_list);

  OOFreeP(I);
}
