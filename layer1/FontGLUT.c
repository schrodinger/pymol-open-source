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

static char *FontGLUTRenderOpenGL(CFontGLUT *I,char *st)
{
  int c;
  FontGLUTBitmapFontRec *font_info = I->glutFont;
  int first,last;
  FontGLUTBitmapCharRec *ch;

  if(st&&(*st)) {
    FontGLUTSave(I);
    
    glColor3fv(TextGetColor());
    glRasterPos4fv(TextGetPos());

    first = font_info->first;
    last = first + font_info->num_chars;

    while((c=*(st++))) {
      if ((c >= first) && (c < last))
        {
          ch = font_info->ch[c - first];
          if (ch) {
            glBitmap(ch->width, ch->height, 
                     ch->xorig, ch->yorig,
                     ch->advance, 0, ch->bitmap);

#if 0
            /* testing code for Character engine */
            { 
              CharFngrprnt fprnt;
              unsigned char *rgba;
              fprnt.u.i.text_id = I->Font.TextID;
              rgba = fprnt.u.i.color;
              TextGetColorUChar(rgba,rgba+1,rgba+2,rgba+3);
              fprnt.u.i.ch = c;
              {
                int id = CharacterFind(&fprnt);
                if(!id) {
                  id = CharacterNewFromBitmap(ch->width, 
                                              ch->height, 
                                              (unsigned char*)ch->bitmap,
                                            &fprnt);
                }
              }
            }
            /* end testing code */
#endif

          }
        }
    }
    FontGLUTRestore(I); 
  }
  return st;
}

static char *FontGLUTRenderRay(CRay *ray, CFontGLUT *I,char *st)
{
  int c;
  FontGLUTBitmapFontRec *font_info = I->glutFont;
  int first,last;
  FontGLUTBitmapCharRec *ch;
  CharFngrprnt fprnt;
  unsigned char *rgba;
  
  if(st&&(*st)) {
    first = font_info->first;
    last = first + font_info->num_chars;
    fprnt.u.i.text_id = I->Font.TextID;
    rgba = fprnt.u.i.color;
    TextGetColorUChar(rgba,rgba+1,rgba+2,rgba+3);

    while((c=*(st++))) {
      if ((c >= first) && (c < last))
        {
          ch = font_info->ch[c - first];
          if (ch) {
            fprnt.u.i.ch = c;
            int id = CharacterFind(&fprnt);
            if(!id) {
              id = CharacterNewFromBitmap(ch->width, 
                                          ch->height, 
                                          (unsigned char*)ch->bitmap,
                                          &fprnt);
            }
            ray->fCharacter(ray,id,
                            (float)ch->xorig,
                            (float)ch->yorig,
                            (float)ch->advance);
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

CFont* FontGLUTNew(int font_code)
{  
  OOAlloc(CFontGLUT);
  FontInit(&I->Font);
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


