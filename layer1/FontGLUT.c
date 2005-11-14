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
#include "Base.h"
#include "OOMac.h"
#include "FontGLUT.h"
#include "Text.h"
#include "Ray.h"
#include "Character.h"
#include "Scene.h"
#include "Util.h"

static void FontGLUTSave(CFontGLUT *I)
{
  glGetIntegerv(GL_UNPACK_SWAP_BYTES,(GLint*)&I->swapbytes);
  glGetIntegerv(GL_UNPACK_LSB_FIRST, (GLint*)&I->lsbfirst);
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, (GLint*)&I->rowlength);
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, (GLint*)&I->skiprows);
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, (GLint*)&I->skippixels);
  glGetIntegerv(GL_UNPACK_ALIGNMENT, (GLint*)&I->alignment);
  
  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
}

static void FontGLUTRestore(CFontGLUT *I)
{
  glPixelStorei(GL_UNPACK_SWAP_BYTES, I->swapbytes);
  glPixelStorei(GL_UNPACK_LSB_FIRST, I->lsbfirst);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, I->rowlength);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, I->skiprows);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, I->skippixels);
  glPixelStorei(GL_UNPACK_ALIGNMENT, I->alignment);
}

static char *FontGLUTRenderOpenGL(RenderInfo *info,CFontGLUT *I,char *st,float size)
{
  register PyMOLGlobals *G = I->Font.G;
  if(G->ValidContext) {
    int c;
    FontGLUTBitmapFontRec *font_info = I->glutFont;
    int first,last;
    FontGLUTBitmapCharRec const *ch;
    int textured = SettingGetGlobal_b(G,cSetting_texture_fonts);
    int pushed = OrthoGetPushed(G);
    int sampling = 1;
    if(info)
      sampling = info->sampling;
    
    if(sampling>1) 
      textured = true;

    if(st&&(*st)) {
      
      if(!textured) {
        glColor3fv(TextGetColor(G));
        glRasterPos4fv(TextGetPos(G));
        FontGLUTSave(I);
      }
      
      first = font_info->first;
      last = first + font_info->num_chars;
      
      if(textured && !pushed) {
        float *v = TextGetPos(G);
        float zero[3]= {0.0F,0.0F,0.0F};
        ScenePushRasterMatrix(G,v);
        TextSetPos(G,zero);
      }
      
      while((c=*(st++))) {
        if ((c >= first) && (c < last))
          {
            ch = font_info->ch[c - first];
            if (ch) {
              if(!textured) {
                
                glBitmap(ch->width, ch->height, 
                         ch->xorig, ch->yorig,
                         ch->advance, 0, ch->bitmap);
                TextAdvance(G,ch->advance);
              } else {
                CharFngrprnt fprnt;
                unsigned char *rgba;
                UtilZeroMem(&fprnt,sizeof(fprnt));
                fprnt.u.i.text_id = I->Font.TextID;
                fprnt.u.i.size = sampling;
                rgba = fprnt.u.i.color;
                TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);
                fprnt.u.i.ch = c;
                {
                  int id = CharacterFind(G,&fprnt);
                  if(!id) {
                    id = CharacterNewFromBitmap(G,ch->width, 
                                                ch->height, 
                                                (unsigned char*)ch->bitmap,
                                                (float)ch->xorig,
                                                (float)ch->yorig,
                                                (float)ch->advance,
                                                &fprnt,sampling);
                  }
                  if(id) {
                    CharacterRenderOpenGL(G,info,id); /* handles advance */
                  }
                }
              }
            }
          }
      }
                
      if(textured && !pushed) {
        ScenePopRasterMatrix(G);
      }
      if(!textured) {
        FontGLUTRestore(I); 
        glFlush(); /* workaround for screen flashes on late-model nVidia hardware */
      }
    }
  }
  return st;
}

static char *FontGLUTRenderRay(CRay *ray, CFontGLUT *I,char *st,float size)
{
  int c;
  FontGLUTBitmapFontRec *font_info = I->glutFont;
  int first,last;
  FontGLUTBitmapCharRec const *ch;
  CharFngrprnt fprnt;
  unsigned char *rgba;
  int sampling = 1;
  sampling = ray->Sampling;
  
  if(st&&(*st)) {
    UtilZeroMem(&fprnt,sizeof(fprnt));
    first = font_info->first;
    last = first + font_info->num_chars;
    fprnt.u.i.text_id = I->Font.TextID;           
    fprnt.u.i.size = sampling;
    rgba = fprnt.u.i.color;
    TextGetColorUChar(I->Font.G,rgba,rgba+1,rgba+2,rgba+3);

    while((c=*(st++))) {
      if ((c >= first) && (c < last))
        {
          ch = font_info->ch[c - first];
          if (ch) {
            fprnt.u.i.ch = c;
            {
              int id = CharacterFind(I->Font.G,&fprnt);
              if(!id) {
                id = CharacterNewFromBitmap(I->Font.G,ch->width, 
                                            ch->height, 
                                            (unsigned char*)ch->bitmap,
                                            (float)ch->xorig,
                                            (float)ch->yorig,
                                            (float)ch->advance,
                                            &fprnt,sampling);
              }
              if(id) ray->fCharacter(ray,id); /* handles advance */
            }
          }
        }
    }
  }
  return st;
}

void FontGLUTFree(CFont *I)
{
  OOFreeP(I);
}

CFont* FontGLUTNew(PyMOLGlobals *G,int font_code)
{  
  OOAlloc(G,CFontGLUT);
  FontInit(G,&I->Font);
  I->Font.fRenderOpenGL = (FontRenderOpenGLFn*)FontGLUTRenderOpenGL;
  I->Font.fRenderRay = (FontRenderRayFn*)FontGLUTRenderRay;
  I->Font.fFree = FontGLUTFree;
  switch(font_code) {
  case cFontGLUT9x15:
    I->glutFont = &FontGLUTBitmap9By15;
    break;
  case cFontGLUTHel10:
    I->glutFont = &FontGLUTBitmapHelvetica10;
    break;
  case cFontGLUTHel12:
    I->glutFont = &FontGLUTBitmapHelvetica12;
    break;
  case cFontGLUTHel18:
    I->glutFont = &FontGLUTBitmapHelvetica18;
    break;
  case cFontGLUT8x13:
  default:
    I->glutFont = &FontGLUTBitmap8By13;
    break;
  }
  return (CFont*)I;
}


