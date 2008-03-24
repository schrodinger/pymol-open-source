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
#include "Matrix.h"

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

static char *FontGLUTRenderOpenGL(RenderInfo *info,CFontGLUT *I,char *st,float size, float *rpos)
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
    const float _0 = 0.0F, _1 = 1.0F, _m1 = -1.0F;
    float x_indent=0.0F, y_indent=0.0F, z_indent=0.0F;

    if(info)
      sampling = info->sampling;
    
    if(sampling>1) 
      textured = true;

    if(st&&(*st)) {

      float origin[3], v_scale;
      SceneOriginGet(G,origin);
      v_scale = SceneGetScreenVertexScale(G,origin);
      
      first = font_info->first;
      last = first + font_info->num_chars;
 
      if(rpos) {
        if(rpos[0]<_1) { /* we need to measure the string width before starting to draw */
          float factor = rpos[0]/2.0F - 0.5F;
          char *sst = st;
          if(factor<_m1) factor = _m1;
          if(factor>_0) factor = _0;
          while((c=*(sst++))) {
            if ((c >= first) && (c < last)) {
              ch = font_info->ch[c - first];
              if (ch) {
                x_indent -= factor * ch->advance;
              }
            }
          }
        }
        if(rpos[0]<_m1) {
          x_indent -= (rpos[0]+_1)/v_scale;
        } else if(rpos[0]>_1) {
          x_indent -= (rpos[0]-_1)/v_scale;
        }
        if(rpos[1]<_1) {
          float factor = -rpos[1]/2.0F + 0.5F;
          if(factor>_1) factor = _1;
          if(factor<_0) factor = _0;
          y_indent = 0.75*size*factor;
        }
        if(rpos[1]<_m1) {
          y_indent -= (rpos[1]+_1)/v_scale;
        } else if(rpos[1]>_1) {
          y_indent -= (rpos[1]-_1)/v_scale;
        }
        z_indent = rpos[2];
        if(z_indent<_0) { /* leave room for fonts of finite depth */
          z_indent+= _1;
          if(z_indent>_0) z_indent = _0;
        } else if(z_indent>_0) {
          z_indent-= _1;
          if(z_indent<_0) z_indent = _0;
        }
      }
      
      if(textured && !pushed) {
        float *v = TextGetPos(G);
        float loc[3];
        float zero[3]= {0.0F,0.0F,0.0F};
        if(rpos) {
          if(info->ortho) {
            float orig[3];
            SceneOriginGet(G,orig);
            SceneGetEyeNormal(G,orig,loc);
          } else {
            SceneGetEyeNormal(G,v,loc);
          }
          scale3f(loc,z_indent,loc);
          add3f(v,loc,loc);
          v = loc;
        }
        ScenePushRasterMatrix(G,v);
        TextSetPos(G,zero);
      } else if(!textured) {
        if(rpos) {
          float *v = TextGetPos(G);
          float loc[3];
          if(info->ortho) {
            float orig[3];
            SceneOriginGet(G,orig);
            SceneGetEyeNormal(G,orig,loc);
          } else {
            SceneGetEyeNormal(G,v,loc);
          }
          scale3f(loc,z_indent,loc);
          add3f(v,loc,loc);
          TextSetPos(G,loc);
        }
      }

      if(rpos) {
        if(textured) {
          TextIndent(G,x_indent,y_indent);
        } else {
          float *v = TextGetPos(G);
          float indent[3];
          float loc[3];
          float *matrix = SceneGetMatrix(G);
          indent[0] = -v_scale * x_indent;
          indent[1] = -v_scale * y_indent;
          indent[2] = _0;
          MatrixInvTransformC44fAs33f3f(matrix,indent,indent);
          add3f(indent,v,loc);
          TextSetPos(G,loc);
        }
      }
      
      if(!textured) {
        glColor3fv(TextGetColor(G));
        glRasterPos4fv(TextGetPos(G));
        FontGLUTSave(I);
      }

      if(textured) 
        CharacterRenderOpenGLPrime(G,info);
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
                fprnt.u.i.ch = (unsigned int)c;
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
      if(textured) 
        CharacterRenderOpenGLDone(G,info);
                
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

static char *FontGLUTRenderRay(CRay *ray, CFontGLUT *I,char *st,float size, float *rpos)
{
  PyMOLGlobals *G = I->Font.G;
  int c;
  FontGLUTBitmapFontRec *font_info = I->glutFont;
  int first,last;
  FontGLUTBitmapCharRec const *ch;
  CharFngrprnt fprnt;
  unsigned char *rgba;
  int sampling = 1;
  float xn[3], yn[3], x_adj[3], y_adj[3], pos[3], *v;
  const float _0 = 0.0F, _1 = 1.0F, _m1 = -1.0F;
  float x_indent=0.0F, y_indent = 0.0F, z_indent = 0.0F;
  sampling = ray->Sampling;
  
  if(st&&(*st)) {
    float origin[3],v_scale;
    SceneOriginGet(G,origin);
    v_scale = SceneGetScreenVertexScale(G,origin);

    if(rpos) {
      float loc[3];
      v = TextGetPos(G);
      if(ray->Ortho) {
        float orig[3];
        SceneOriginGet(G,orig);
        SceneGetEyeNormal(G,orig,loc);
      } else { 
        SceneGetEyeNormal(G,v,loc);
      }
      scale3f(loc,rpos[2],loc);
      add3f(v,loc,loc);
      TextSetPos(G,loc);
    }

    RayGetScaledAxes(ray,xn,yn);

    UtilZeroMem(&fprnt,sizeof(fprnt));
    first = font_info->first;
    last = first + font_info->num_chars;
    fprnt.u.i.text_id = I->Font.TextID;           
    fprnt.u.i.size = sampling;
    rgba = fprnt.u.i.color;
    TextGetColorUChar(G,rgba,rgba+1,rgba+2,rgba+3);

    if(rpos) {

      if(rpos[0]<_1) { /* we need to measure the string width before starting to draw */
        float factor = rpos[0]/2.0F - 0.5F;
        char *sst = st;
        if(factor<_m1) factor = -_1;
        if(factor>_0) factor = _0;
        
        while((c=*(sst++))) {
          fprnt.u.i.ch = (unsigned int) c;
          ch = font_info->ch[c - first];
          if (ch) {
            x_indent -= 2*factor * ch->advance;
          }
        }
      }
      if(rpos[0]<_m1) {
        x_indent -= 2*(rpos[0]+_1)/v_scale;
      } else if(rpos[0]>_1) {
        x_indent -= 2*(rpos[0]-_1)/v_scale;
      }
      if(rpos[1]<_1) {
        float factor = -rpos[1]/2.0F + 0.5F;
        if(factor>_1) factor = _1;
        if(factor<_0) factor = _0;
        y_indent = 0.75F*sampling*size*factor;
      }
      if(rpos[1]<_m1) {
        y_indent -= 2*(rpos[1]+_1)/v_scale;
      } else if(rpos[1]>_1) {
        y_indent -= 2*(rpos[1]-_1)/v_scale;
      }
      z_indent = rpos[2];
      if(z_indent<_0) { /* leave room for fonts of finite depth */
        z_indent+= _1;
        if(z_indent>_0) z_indent = _0;
      } else if(z_indent>_0) {
        z_indent-= _1;
        if(z_indent<_0) z_indent = _0;
      }
      v = TextGetPos(G);
      scale3f(xn, x_indent, x_adj);
      scale3f(yn, y_indent, y_adj);
      subtract3f(v,x_adj,pos);
      subtract3f(pos,y_adj,pos);
      TextSetPos(G,pos);
    }

    while((c=*(st++))) {
      if ((c >= first) && (c < last))
        {
          ch = font_info->ch[c - first];
          if (ch) {
            fprnt.u.i.ch = (unsigned int)c;
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
  I->Font.fRenderOpenGLFlat = (FontRenderOpenGLFn*)FontGLUTRenderOpenGL;
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


