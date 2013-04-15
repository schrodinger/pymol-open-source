

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
#include"os_gl.h"
#include"ShaderMgr.h"
#include"Executive.h"
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

#define POS_START  2

typedef struct  {
  int id,dim;
} texture_info;

struct _CTexture {
  OVOneToOne *ch2tex;
  GLuint text_texture_id;
  int xpos, ypos, maxypos;
  int num_chars;
  int text_texture_dim;
};

GLuint TextureGetTextTextureID(PyMOLGlobals * G){
  CTexture *I = G->Texture;
  return I->text_texture_id;
}
int TextureGetTextTextureSize(PyMOLGlobals * G){
  CTexture *I = G->Texture;
  return I->text_texture_dim;
}

#define INIT_TEXTURE_SIZE 512

int TextureInit(PyMOLGlobals * G)
{
  OOAlloc(G, CTexture);

  G->Texture = I;

  I->ch2tex = OVOneToOne_New(G->Context->heap);
  I->text_texture_id = 0;
  I->text_texture_dim = I->ypos = I->maxypos = I->num_chars = 0;
  I->xpos = POS_START;
  return (I ? 1 : 0);
}

void TextureInitTextTexture(PyMOLGlobals *G){
  short is_new = 0;
  CTexture *I = G->Texture;
  if (!I->text_texture_id){
    glGenTextures(1, &I->text_texture_id);
    is_new = 1;
  }
  if(I->text_texture_id){
    if (CShaderMgr_ShadersPresent(G->ShaderMgr)){
#ifdef OPENGL_ES_1
      glActiveTexture(GL_TEXTURE1);
#else
      glActiveTexture(GL_TEXTURE3);
#endif
    }
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_2D, I->text_texture_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    if (is_new){
      int tex_dim = INIT_TEXTURE_SIZE;
      int buff_total = tex_dim * tex_dim;
      unsigned char *temp_buffer = Alloc(unsigned char, buff_total * 4);
      UtilZeroMem(temp_buffer, buff_total * 4);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
		   tex_dim, tex_dim, 0, GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid*)temp_buffer);
      I->text_texture_dim = INIT_TEXTURE_SIZE;
      FreeP(temp_buffer);
      I->xpos = POS_START; I->ypos = 0; I->maxypos = POS_START;
    }
  }
}

int TextureGetFromChar(PyMOLGlobals * G, int char_id, float *extent)
{
  OVreturn_word result;
  CTexture *I = G->Texture;
  int is_new = false;
  int tex_dim = INIT_TEXTURE_SIZE;
  short use_shader = (short) SettingGetGlobal_b(G, cSetting_use_shaders);

  if(G->HaveGUI && G->ValidContext) {
    if(OVreturn_IS_OK(result = OVOneToOne_GetForward(I->ch2tex, char_id))) {
      if(glIsTexture(I->text_texture_id))
        return I->text_texture_id;
      else {
        OVOneToOne_DelReverse(I->ch2tex, result.word);
      }
    }
    {
      unsigned char *buffer = NULL;
      if (!I->text_texture_id){
          is_new = true;
      }
      buffer = CharacterGetPixmapBuffer(G, char_id);
      if(buffer) {
        int w = CharacterGetWidth(G, char_id);
        int h = CharacterGetHeight(G, char_id);
        GLuint texture_id = 0;
        int buff_incr = is_new ? tex_dim : w;
        int buff_total = is_new ? tex_dim * tex_dim : w * h;
        unsigned char *temp_buffer = Alloc(unsigned char, buff_total * 4);

        {
          int a, b;
          unsigned char *p = buffer, *q;
	  int fa = 0, ta = w;
	  if (is_new){
	    fa += I->xpos; ta += I->xpos;
	  }
          UtilZeroMem(temp_buffer, buff_total * 4);
          for(b = 0; b < h; b++) {
	    for(a = fa; a < ta; a++) {
              q = temp_buffer + (4 * buff_incr * b) + 4 * a;
              *(q++) = *(p++);
              *(q++) = *(p++);
              *(q++) = *(p++);
              *(q++) = *(p++);
            }
          }
	  if (I->xpos + w > tex_dim){
	    // if the size of the texture goes off the side, go to next row
	    I->xpos = 0;
	    I->ypos = I->maxypos;
	  }
	  if ((I->xpos + w) >= INIT_TEXTURE_SIZE && (I->ypos + h) >= INIT_TEXTURE_SIZE){
	    I->xpos = POS_START; I->ypos = 0; I->maxypos = POS_START;
	    OVOneToOne_Reset(I->ch2tex);
	    /* Also need to reload the selection markers into the texture, since
	       we are wiping everything out from the texture and starting from the origin */
	    ExecutiveInvalidateSelectionIndicators(G);
	  }
          extent[0] = (I->xpos / (float) tex_dim);// + .0002f;
          extent[1] = (I->ypos / (float) tex_dim);// + .0002f;
          extent[2] = ((I->xpos + w) / (float) tex_dim);// - .0002f;
          extent[3] = ((I->ypos + h) / (float) tex_dim);// - .0002f;
        }

	if (!I->text_texture_id){
          glGenTextures(1, &I->text_texture_id);
	}
	texture_id = I->text_texture_id;
	if(I->text_texture_id && OVreturn_IS_OK(OVOneToOne_Set(I->ch2tex, char_id, I->num_chars++))) {
	  if (use_shader && CShaderMgr_ShadersPresent(G->ShaderMgr)){
#ifdef OPENGL_ES_1
	    glActiveTexture(GL_TEXTURE1);
#else
	    glActiveTexture(GL_TEXTURE3);
#endif
	  }
          glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
          glBindTexture(GL_TEXTURE_2D, texture_id);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
          if(is_new) {
	    I->text_texture_dim = tex_dim;
	      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
			   tex_dim, tex_dim, 0, GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid*)temp_buffer);
          } else {
	    int xoff = 0, yoff = 0;
	    xoff = I->xpos;
	    yoff = I->ypos;
	      glTexSubImage2D(GL_TEXTURE_2D, 0, xoff, yoff,
			      w, h, GL_RGBA, GL_UNSIGNED_BYTE, temp_buffer);
          }
        }
	if (I->ypos + h > I->maxypos){
	  I->maxypos = I->ypos + h + 1;  // added space for running on Ipad/Iphone (weird artifacts)
	}
	if (I->xpos + w > tex_dim){
	  I->xpos = 0;
	  I->ypos = I->maxypos;
	} else {
	  I->xpos += w + 1; // added space for running on Ipad/Iphone (weird artifacts)
	}
          FreeP(temp_buffer);
        return texture_id;
      }
    }
  }
  return 0;
}

void TextureGetPlacementForNewSubtexture(PyMOLGlobals * G, int new_texture_width, int new_texture_height, int *new_texture_posx, int *new_texture_posy){
  CTexture *I = G->Texture;
  if (I->xpos + new_texture_width > INIT_TEXTURE_SIZE){
    I->xpos = 0;
    I->ypos = I->maxypos;
  }
  if (I->ypos + new_texture_height > I->maxypos){
    I->maxypos = I->ypos + new_texture_height + 1;  // added space for running on Ipad/Iphone (weird artifacts)
  }
  *new_texture_posx = I->xpos;
  *new_texture_posy = I->ypos;
  I->xpos += new_texture_width + 1; // added space for running on Ipad/Iphone (weird artifacts)
  
}

void TextureFree(PyMOLGlobals * G)
{
  CTexture *I = G->Texture;
  /* TODO -- free all the resident textures */

  OVOneToOne_DEL_AUTO_NULL(I->ch2tex);

  OOFreeP(I);
}
