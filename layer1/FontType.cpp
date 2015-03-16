
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
#include"os_python.h"

#include "MemoryDebug.h"
#include "OOMac.h"
#include "os_gl.h"
#include "FontType.h"
#include "Text.h"
#include "Ortho.h"
#include "Scene.h"
#include "Character.h"
#include "Util.h"
#include "TypeFace.h"

typedef struct {
  CFont Font;                   /* must be first */
  PyMOLGlobals *G;
  CTypeFace *TypeFace;
} CFontType;

#ifdef _PYMOL_INLINE
__inline__
#endif
static const char *_FontTypeRenderOpenGL(RenderInfo * info,
                                   CFontType * I, const char *st,
                                   float size, int flat, float *rpos SHADERCGOARG)
{
  PyMOLGlobals *G = I->Font.G;
  if(G->ValidContext) {
    unsigned int c;
    int pushed = OrthoGetPushed(G);
    int kern_flag = false;
    unsigned int last_c = 0;
    int sampling = 1;
    const float _0 = 0.0F, _1 = 1.0F, _m1 = -1.0F;
    float x_indent = 0.0F, y_indent = 0.0F, z_indent = 0.0F;
    float text_width = 0.f;
    int unicode = 0;
    int unicnt = 0;

    sampling = info->sampling;
    if(st && (*st)) {
      float v_scale;
      float screenWorldOffset[3] = { 0.0F, 0.0F, 0.0F };
      
      v_scale = SceneGetScreenVertexScale(G, NULL);
      if(size < _0) {
        size = (int) (0.5F - size / v_scale);
      }

      if(rpos) {
        if(rpos[0] < _1) {      /* we need to measure the string width before starting to draw */
          const char *sst = st;
          while((c = *(sst++))) {
            if(unicnt) {
              if(!(c & 0x80))   /* corrupt UTF8 */
                unicnt = 0;
              else {
                unicode = (unicode << 6) | (0x3F & c);
                unicnt--;
                c = unicode;
              }
            } else if(c & 0x80) {
              while(c & 0x80) {
                c = (c << 1) & 0xFF;
                unicnt++;
              }
              unicode = (c >> (unicnt--));
            }
            if(!unicnt) {
              CharFngrprnt fprnt;
              unsigned char *rgba;
              UtilZeroMem(&fprnt, sizeof(fprnt));
              fprnt.u.i.text_id = I->Font.TextID;
              fprnt.u.i.size = (int) (size * 64 * sampling);
              rgba = fprnt.u.i.color;
              TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
              rgba = fprnt.u.i.outline_color;
              if(!flat) {
                TextGetOutlineColor(G, rgba, rgba + 1, rgba + 2, rgba + 3);
              } else {
                TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
              }
              fprnt.u.i.ch = c;
              fprnt.u.i.flat = flat;
              {
                int id = CharacterFind(G, &fprnt);
                if(!id) {
                  id = TypeFaceCharacterNew(I->TypeFace, &fprnt, size * sampling);
                }
                if(id) {
                  if(kern_flag) {
                    text_width += (TypeFaceGetKerning(I->TypeFace,
							       last_c, c, size) / sampling);
                  }
                  text_width += CharacterGetAdvance(G, sampling, id);
                }
              }
              kern_flag = true;
              last_c = c;
            }
          }
	  {
	    float factor = rpos[0] / 2.0F - 0.5F;
	    /* if -1. < rpos[0] < 1. , then we need to determine the label's width
	       so that we can justify it appropriately */
	    if(factor < _m1) // if rpos[0] < -1., right justified
	      factor = -_1;
	    if(factor > _0)  // if rpos[0] > 1., left justified, label width not needed
	      factor = _0;
	    x_indent -= factor * text_width;
	  }
        }
	/* if label_position x is -1 to 1, the label is placed in x such that
	   0 - centered 
	   -1 - right justified on projected point
	    1 - left justified on projected point
	*/
	if(rpos[0] < _m1) {
	  screenWorldOffset[0] -= (rpos[0] + _1);// / v_scale;
	} else if(rpos[0] > _1) {
	  screenWorldOffset[0] -= (rpos[0] - _1);// / v_scale;
	}
        if(rpos[1] < _1) {
          float factor = -rpos[1] / 2.0F + 0.5F;
          if(factor > _1)
            factor = _1;
          if(factor < _0)
            factor = _0;
          y_indent = 0.75 * size * factor;
        }
        if(rpos[1] < _m1) {
          screenWorldOffset[1] -= (rpos[1] + _1);
        } else if(rpos[1] > _1) {
          screenWorldOffset[1] -= (rpos[1] - _1);
        }
	/* leave room for fonts of finite depth */
	z_indent = (rpos[2] < _m1) ? (rpos[2]+1.f) : (rpos[2] > 1.f) ? (rpos[2]-1.f) : 0.f;

	if (!shaderCGO){
	  x_indent += screenWorldOffset[0] / v_scale;
	  y_indent += screenWorldOffset[1] / v_scale;
	}
	screenWorldOffset[2] += z_indent; // need to take into account weird -1 to 1, and sub 1 from abs val 
      }
      if(!pushed) {
        float *v = TextGetPos(G);
        float loc[3];
        float zero[3] = { 0.0F, 0.0F, 0.0F };
	TextSetScreenWorldOffset(G, screenWorldOffset);
        TextSetWorldPos(G, v);
        if(rpos) {
          if(info->ortho) {
            float orig[3];
            SceneOriginGet(G, orig);
            SceneGetEyeNormal(G, orig, loc);
          } else {
            SceneGetEyeNormal(G, v, loc);
          }
          scale3f(loc, z_indent, loc);
          add3f(v, loc, loc);
          v = loc;
        }
        ScenePushRasterMatrix(G, v);
        TextSetPos(G, zero);
      }
      if(rpos) {
	//	float ax = x_indent * v_scale, ay = y_indent * v_scale;
	//	printf("x_indent=%f y_indent=%f ax=%f ay=%f v_scale=%f rpos=%f %f %f\n", x_indent, y_indent, ax, ay, v_scale, rpos[0], rpos[1], rpos[2]);
        TextIndent(G, x_indent, y_indent);
      }
      CharacterRenderOpenGLPrime(G, info);
      while((c = *(st++))) {
        if(unicnt) {
          if(!(c & 0x80))       /* corrupt UTF8 */
            unicnt = 0;
          else {
            unicode = (unicode << 6) | (0x3F & c);
            unicnt--;
            c = unicode;
          }
        } else if(c & 0x80) {
          while(c & 0x80) {
            c = (c << 1) & 0xFF;
            unicnt++;
          }
          unicode = (c >> (unicnt--));
        }
        if(!unicnt) {
          CharFngrprnt fprnt;
          unsigned char *rgba;
          UtilZeroMem(&fprnt, sizeof(fprnt));
          fprnt.u.i.text_id = I->Font.TextID;
          fprnt.u.i.size = (int) (size * 64 * sampling);
          rgba = fprnt.u.i.color;
          TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
          rgba = fprnt.u.i.outline_color;
          TextGetOutlineColor(G, rgba, rgba + 1, rgba + 2, rgba + 3);
          fprnt.u.i.ch = c;
          fprnt.u.i.flat = flat;
          {
            int id = CharacterFind(G, &fprnt);
            if(!id) {
              id = TypeFaceCharacterNew(I->TypeFace, &fprnt, size * sampling);
            }
            if(id) {
              if(kern_flag) {
                TextAdvance(G, TypeFaceGetKerning(I->TypeFace,
                                                  last_c, c, size) / sampling);
              }
              CharacterRenderOpenGL(G, info, id, true SHADERCGOARGVAR);       /* handles advance */
            }
          }
          kern_flag = true;
          last_c = c;
        }
      }
      CharacterRenderOpenGLDone(G, info);
      if(!pushed) {
        ScenePopRasterMatrix(G);
      }
    }
  }
  return st;
}

static const char *FontTypeRenderOpenGL(RenderInfo * info, CFontType * I, const char *st, float size,
                                  float *rpos SHADERCGOARG)
{
  return _FontTypeRenderOpenGL(info, I, st, size, false, rpos SHADERCGOARGVAR);
}

static const char *FontTypeRenderOpenGLFlat(RenderInfo * info, CFontType * I, const char *st,
                                      float size, float *rpos SHADERCGOARG)
{
  return _FontTypeRenderOpenGL(info, I, st, size, true, rpos SHADERCGOARGVAR);
}

static const char *FontTypeRenderRay(CRay * ray, CFontType * I, const char *st, float size,
                               float *rpos)
{
  PyMOLGlobals *G = I->Font.G;
  unsigned int c;
  int kern_flag = false;
  unsigned int last_c = 0;
  int sampling = ray->Sampling;
  const float _0 = 0.0F, _1 = 1.0F, _m1 = -1.0F;
  float x_indent = 0.0F, y_indent = 0.0F, z_indent = 0.0F;
  float xn[3], yn[3], x_adj[3], y_adj[3], pos[3], *v;
  int unicode = 0;
  int unicnt = 0;

  if(st && (*st)) {
    float v_scale = SceneGetScreenVertexScale(G, NULL);

    if(rpos) {
      float loc[3];
      /* leave room for fonts of finite depth */
      z_indent = (rpos[2] < _m1) ? (rpos[2]+1.f) : (rpos[2] > 1.f) ? (rpos[2]-1.f) : 0.f;
      v = TextGetPos(I->G);

      if(ray->Ortho) {
        float orig[3];
        SceneOriginGet(G, orig);
        SceneGetEyeNormal(G, orig, loc);
      } else {
        SceneGetEyeNormal(G, v, loc);
      }
      scale3f(loc, z_indent, loc);
      add3f(v, loc, loc);
      TextSetPos(I->G, loc);
    }

    RayGetScaledAxes(ray, xn, yn);

    if(size < _0) {

      size = (int) (0.5F - size / v_scale);
    }

    if(rpos) {

      if(rpos[0] < _1) {        /* we need to measure the string width before starting to draw */
        float factor = rpos[0] / 2.0F - 0.5F;
        const char *sst = st;
	factor = (factor < _m1) ? _m1 : (factor > _0) ? _0 : factor;
        while((c = *(sst++))) {
          if(unicnt) {
            if(!(c & 0x80))     /* corrupt UTF8 */
              unicnt = 0;
            else {
              unicode = (unicode << 6) | (0x3F & c);
              unicnt--;
              c = unicode;
            }
          } else if(c & 0x80) {
            while(c & 0x80) {
              c = (c << 1) & 0xFF;
              unicnt++;
            }
            unicode = (c >> (unicnt--));
          }
          if(!unicnt) {
            CharFngrprnt fprnt;
            unsigned char *rgba;
            UtilZeroMem(&fprnt, sizeof(fprnt));
            fprnt.u.i.text_id = I->Font.TextID;
            fprnt.u.i.size = (int) (size * 64 * sampling);
            rgba = fprnt.u.i.color;
            TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
            rgba = fprnt.u.i.outline_color;
            TextGetOutlineColor(G, rgba, rgba + 1, rgba + 2, rgba + 3);
            fprnt.u.i.ch = c;
            {
              int id = CharacterFind(G, &fprnt);
              if(!id) {
                id = TypeFaceCharacterNew(I->TypeFace, &fprnt, size * sampling);
              }
              if(id) {
                if(kern_flag) {
                  x_indent -= factor * TypeFaceGetKerning(I->TypeFace,
                                                          last_c,
                                                          c, size * sampling) / sampling;
                }
                x_indent -= factor * CharacterGetAdvance(G, 1, id);
                kern_flag = true;
                last_c = c;
              }
            }
          }
        }
      }
      if(rpos[0] < _m1) {
        x_indent -= 2 * (rpos[0] + _1) / v_scale;
      } else if(rpos[0] > _1) {
        x_indent -= 2 * (rpos[0] - _1) / v_scale;
      }
      if(rpos[1] < _1) {
        float factor = -rpos[1] / 2.0F + 0.5F;
        if(factor > _1)
          factor = _1;
        if(factor < _0)
          factor = _0;
        y_indent = 0.75F * sampling * size * factor;
      }
      if(rpos[1] < _m1) {
        y_indent -= 2 * (rpos[1] + _1) / v_scale;
      } else if(rpos[1] > _1) {
        y_indent -= 2 * (rpos[1] - _1) / v_scale;
      }
      v = TextGetPos(I->G);
      scale3f(xn, x_indent, x_adj);
      scale3f(yn, y_indent, y_adj);
      subtract3f(v, x_adj, pos);
      subtract3f(pos, y_adj, pos);
      TextSetPos(I->G, pos);
    }
    kern_flag = false;

    while((c = *(st++))) {
      if(unicnt) {
        if(!(c & 0x80))         /* corrupt UTF8 */
          unicnt = 0;
        else {
          unicode = (unicode << 6) | (0x3F & c);
          unicnt--;
          c = unicode;
        }
      } else if(c & 0x80) {
        while(c & 0x80) {
          c = (c << 1) & 0xFF;
          unicnt++;
        }
        unicode = (c >> (unicnt--));
      }
      if(!unicnt) {
        CharFngrprnt fprnt;
        unsigned char *rgba;
        UtilZeroMem(&fprnt, sizeof(fprnt));
        fprnt.u.i.text_id = I->Font.TextID;
        fprnt.u.i.size = (int) (size * 64 * sampling);
        rgba = fprnt.u.i.color;
        TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
        rgba = fprnt.u.i.outline_color;
        TextGetOutlineColor(G, rgba, rgba + 1, rgba + 2, rgba + 3);
        fprnt.u.i.ch = c;
        {
          int id = CharacterFind(G, &fprnt);
          if(!id) {
            id = TypeFaceCharacterNew(I->TypeFace, &fprnt, size * sampling);
          }
          if(id) {
            if(kern_flag) {
              float kern = TypeFaceGetKerning(I->TypeFace,
                                              last_c,
                                              c,
                                              size * sampling) / sampling;
              v = TextGetPos(I->G);
              scale3f(xn, kern, x_adj);
              add3f(v, x_adj, pos);
              TextSetPos(I->G, pos);
            }
            ray->character(id);   /* handles advance */

            kern_flag = true;
            last_c = c;
          }
        }
      }
    }
  }
  return st;
}

static void FontTypeFree(CFont * font)
{
  CFontType *I = (CFontType *) font;
  TypeFaceFree(I->TypeFace);
  OOFreeP(I);
}

CFont *FontTypeNew(PyMOLGlobals * G, unsigned char *dat, unsigned int len)
{
  OOAlloc(G, CFontType);
  FontInit(G, &I->Font);
  I->G = G;
  I->Font.fRenderOpenGL = (FontRenderOpenGLFn *) FontTypeRenderOpenGL;
  I->Font.fRenderOpenGLFlat = (FontRenderOpenGLFn *) FontTypeRenderOpenGLFlat;
  I->Font.fRenderRay = (FontRenderRayFn *) FontTypeRenderRay;
  I->Font.fFree = FontTypeFree;
  I->TypeFace = TypeFaceLoad(G, dat, len);
  if(!I->TypeFace) {
    OOFreeP(I);
  }
  return (CFont *) I;
}
