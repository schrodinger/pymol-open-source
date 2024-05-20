

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

#include <memory>
#include <unordered_set>

#include"os_python.h"
#include"os_gl.h"

#include "Base.h"
#include "PyMOLGlobals.h"
#include "Texture.h"

#include "Setting.h"
#include "Character.h"
#include "Util.h"
#include "Feedback.h"

#include "Executive.h"
#include "Ortho.h"
#include "ShaderMgr.h"
#include "GenericBuffer.h"

#define POS_START  2

/**
 * Texture unit assigned to the global text texture.
 * See ShaderMgr.cpp for more details.
 */
constexpr std::uint8_t GlobalTextTextureUnit = 3;

struct CTexture {
  std::unordered_set<int> texturedCharIDs; // set of char ids that have been textured
  std::unique_ptr<textureBuffer_t> texture;
  int xpos{};
  int ypos{};
  int maxypos{};
  int text_texture_dim{};
};

void TextureBindTexture(PyMOLGlobals* G) {
  auto I = G->Texture;
  if (I->texture) {
    I->texture->bind();
  }
}

int TextureGetTextTextureSize(PyMOLGlobals * G){
  CTexture *I = G->Texture;
  return I->text_texture_dim;
}

#define INIT_TEXTURE_SIZE 512

int TextureInit(PyMOLGlobals * G)
{
  auto I = new CTexture();
  G->Texture = I;

  I->text_texture_dim = INIT_TEXTURE_SIZE;
  I->ypos = I->maxypos = 0;
  I->xpos = POS_START;
  return (I ? 1 : 0);
}

static
void TextureInitTextTextureImpl(PyMOLGlobals *G, int textureSize);

void TextureInitTextTexture(PyMOLGlobals *G){
  TextureInitTextTextureImpl(G, INIT_TEXTURE_SIZE);
}
void TextureInvalidateTextTexture(PyMOLGlobals * G){
  CTexture *I = G->Texture;
  if (I->texture) {
    I->texturedCharIDs.clear();
    I->texture.reset();
    I->text_texture_dim = INIT_TEXTURE_SIZE;
    I->xpos = POS_START; I->ypos = 0; I->maxypos = POS_START;
  }
}

void TextureInitTextTextureImpl(PyMOLGlobals *G, int textureSizeArg){
  short is_new = 0;
  CTexture *I = G->Texture;
  int textureSize = textureSizeArg;
  if (!textureSize)
    textureSize = INIT_TEXTURE_SIZE;
  if (!I->texture) {
    using namespace tex;
    I->texture =
        std::make_unique<textureBuffer_t>(format::RGBA, data_type::UBYTE,
            filter::NEAREST, filter::NEAREST, wrap::CLAMP, wrap::CLAMP);
    is_new = true;
  }
  if (I->texture) {
    if (is_new) {
      int tex_dim = textureSize;
      int buff_total = tex_dim * tex_dim;
      std::vector<unsigned char> temp_buffer(
          buff_total * GetSizeOfVertexFormat(VertexFormat::UByte4), 0);
      I->texture->bindToTextureUnit(GlobalTextTextureUnit);
      I->texture->texture_data_2D(tex_dim, tex_dim, temp_buffer.data());
      I->text_texture_dim = textureSize;
      I->xpos = POS_START; I->ypos = 0; I->maxypos = POS_START;
    }
  }
}

#include "Rep.h"
bool TextureIsCharTextured(PyMOLGlobals* G, int char_id, float* extent)
{
  CTexture *I = G->Texture;
  int is_new = false;
  int tex_dim = I->text_texture_dim;
  short use_shader = (short) SettingGetGlobal_b(G, cSetting_use_shaders);

  if(G->HaveGUI && G->ValidContext) {
    auto it = I->texturedCharIDs.find(char_id);
    if (it != I->texturedCharIDs.end()) {
      if (I->texture) {
        return true;
      }
      else {
        I->texturedCharIDs.erase(char_id);
      }
    }
    {
      unsigned char *buffer = nullptr;
      if (!I->texture) {
          is_new = true;
      }
      buffer = CharacterGetPixmapBuffer(G, char_id);
      if(buffer) {
        int w = CharacterGetWidth(G, char_id);
        int h = CharacterGetHeight(G, char_id);
        int buff_incr = is_new ? tex_dim : w;
        int buff_total = is_new ? tex_dim * tex_dim : w * h;
        std::vector<unsigned char> temp_buffer(
            buff_total * GetSizeOfVertexFormat(VertexFormat::UByte4), 0);

        {
          int a, b;
          unsigned char *p = buffer, *q;
	  int fa = 0, ta = w;
	  if (is_new){
	    fa += I->xpos; ta += I->xpos;
	  }
          for(b = 0; b < h; b++) {
	    for(a = fa; a < ta; a++) {
              q = temp_buffer.data() + (4 * buff_incr * b) + 4 * a;
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
	  if ((I->ypos + h) >= I->text_texture_dim){ // only need to check y since x gets reset above
	    int nrefreshes;
	    I->xpos = POS_START; I->ypos = 0; I->maxypos = POS_START;
	    I->texturedCharIDs.clear();
	    /* Also need to reload the selection markers into the texture, since
	       we are wiping everything out from the texture and starting from the origin */
	    if ((nrefreshes=SceneIncrementTextureRefreshes(G)) > 1){
	      /* Texture was refreshed more than once for this frame, increase size of texture */
	      int newDim = I->text_texture_dim * 2;
	      I->texture.reset();
	      TextureInitTextTextureImpl(G, newDim);
	      PRINTFB(G, FB_OpenGL, FB_Output)
		" Texture OpenGL: nrefreshes=%d newDim=%d\n", nrefreshes, newDim ENDFB(G);

	      //	      printf("nrefreshes=%d newDim=%d\n", nrefreshes, newDim);
	      I->xpos = POS_START; I->ypos = 0; I->maxypos = POS_START;
	      SceneResetTextureRefreshes(G);
	    }
	    ExecutiveInvalidateRep(G, "all", cRepLabel, cRepInvRep);
	    ExecutiveInvalidateSelectionIndicators(G);
	    OrthoInvalidateDoDraw(G);
	    return false;
	  }
          extent[0] = (I->xpos / (float) tex_dim);
          extent[1] = (I->ypos / (float) tex_dim);
          extent[2] = ((I->xpos + w) / (float) tex_dim);
          extent[3] = ((I->ypos + h) / (float) tex_dim);
        }

        if (!I->texture) {
          using namespace tex;
          I->texture =
              std::make_unique<textureBuffer_t>(format::RGBA, data_type::UBYTE,
                  filter::NEAREST, filter::NEAREST, wrap::CLAMP, wrap::CLAMP);
        }
	if (I->texture) {
	  I->texturedCharIDs.emplace(char_id);
          if (is_new) {
	    I->text_texture_dim = tex_dim;
        I->texture->bindToTextureUnit(GlobalTextTextureUnit);
        I->texture->texture_data_2D(tex_dim, tex_dim, temp_buffer.data());
          } else {
	    int xoff = 0, yoff = 0;
	    xoff = I->xpos;
	    yoff = I->ypos;
        I->texture->texture_subdata_2D(xoff, yoff, w, h, temp_buffer.data());
#ifdef _WEBGL
              static bool error_flag = false;
              if ((glGetError()) != 0 && error_flag) {
                PRINTFB(G, FB_OpenGL, FB_Output)
                  " IE11 texSubImage2D bug detected, supressing...\n" ENDFB(G);
                error_flag = true;
              }
#endif
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
        return I->texture != nullptr;
      }
    }
  }
  return true;
}

void TextureGetPlacementForNewSubtexture(PyMOLGlobals * G, int new_texture_width, int new_texture_height, int *new_texture_posx, int *new_texture_posy){
  CTexture *I = G->Texture;
  if (I->xpos + new_texture_width > I->text_texture_dim){
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

void TextureFillNewSubtexture(PyMOLGlobals* G, int width, int height, int x, int y, const void* data)
{
  CTexture *I = G->Texture;
  if (I->texture) {
    I->texture->texture_subdata_2D(x, y, width, height, data);
  }
}

void TextureFree(PyMOLGlobals * G)
{
  /* TODO -- free all the resident textures */
  DeleteP(G->Texture);
}
