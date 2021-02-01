
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

#include <algorithm>

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
#include "FontType.h"

#define max2 std::max

static void CheckUnicode(unsigned int *c, int *unicnt, int *unicode){
  if(*unicnt) {
    if(!(*c & 0x80))   /* corrupt UTF8 */
      *unicnt = 0;
    else {
      *unicode = ((*unicode) << 6) | (0x3F & *c);
      (*unicnt)--;
      *c = *unicode;
    }
  } else if(*c & 0x80) {
    while(*c & 0x80) {
      *c = (*c << 1) & 0xFF;
      (*unicnt)++;
    }
    *unicode = (*c >> ((*unicnt)--));
  }
}

static void GenerateCharFngrprnt(PyMOLGlobals *G, CharFngrprnt *fprnt, unsigned int c, int TextID, float size, int sampling, short no_flat, int flat){
  unsigned char *rgba;
  UtilZeroMem(fprnt, sizeof(CharFngrprnt));
  fprnt->u.i.text_id = TextID;
  fprnt->u.i.size = (int) (size * 64 * sampling);
  rgba = fprnt->u.i.color;
  if (!TextGetIsPicking(G)){
    TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
    rgba = fprnt->u.i.outline_color;
    if (no_flat || !flat){
      TextGetOutlineColor(G, rgba, rgba + 1, rgba + 2, rgba + 3);
    } else if(flat) {
      TextGetColorUChar(G, rgba, rgba + 1, rgba + 2, rgba + 3);
    }
  }
  fprnt->u.i.ch = c;
  fprnt->u.i.flat = flat;
}

static const char* FontTypeRenderOpenGLImpl(const RenderInfo* info, CFontType* I,
    const char* st, float size, int flat, const float* rpos, bool needSize,
    short relativeMode, bool shouldRender, CGO* shaderCGO)
{
  PyMOLGlobals *G = I->G;
  if(G->ValidContext) {
    unsigned int c;
    int pushed = OrthoGetPushed(G);
    int kern_flag = false;
    unsigned int last_c = 0;
    int sampling = 1;
    const float _0 = 0.0F, _1 = 1.0F, _m1 = -1.0F;
    float x_indent = 0.0F, y_indent = 0.0F, z_indent = 0.0F;
    float text_width = 0.f, line_width = 0.f, tot_text_width;
    float text_just = 1.f - TextGetJustification(G);
    float text_spacing = TextGetSpacing(G);
    float text_buffer[2];
    int unicode = 0;
    int unicnt = 0;
    int nlines = countchrs(st, '\n') + 1;
    int linenum = 0;
    short cont = 0;
    float descender = TypeFaceGetDescender(I->TypeFace) / 2.f;
    float v_scale = SceneGetScreenVertexScale(G, NULL);
    copy2f(TextGetLabelBuffer(G), text_buffer);
    sampling = info->sampling;
    if(st && (*st)) {
      float *line_widths = NULL;
      float screenWorldOffset[3] = { 0.0F, 0.0F, 0.0F };
      float tot_height;
      if (nlines>1){
	line_widths = pymol::calloc<float>(nlines);
      }
      if(size < _0) {
        size = (int) (0.5F - size / v_scale);
        if (size <= 0)
          size = 1;
      } else {
        size = DIP2PIXEL(size);
      }

      text_buffer[0] *= size;
      text_buffer[1] *= size;

      if(rpos) {
	TextSetIndentFactorX(G, rpos[0] < _m1 ? _m1 : rpos[0] > 1.f ? 1.f : rpos[0]);
	TextSetIndentFactorY(G, rpos[1] < _m1 ? _m1 : rpos[1] > 1.f ? 1.f : rpos[1]);

        if(needSize || rpos[0] < _1) {      /* we need to measure the string width before starting to draw */
          const char *sst = st;
          while((c = *(sst++))) {
	    if (c == '\n'){
	      line_widths[linenum] = line_width;
	      text_width = max2(text_width, line_width);
	      line_width = 0.f;
	      linenum++;
	      continue;
	    }
	    CheckUnicode(&c, &unicnt, &unicode);

            if(!unicnt) {
              CharFngrprnt fprnt;
	      GenerateCharFngrprnt(G, &fprnt, c, I->TextID, size, sampling, 0, flat);
              {
                int id = CharacterFind(G, &fprnt);
                if(!id) {
                  id = TypeFaceCharacterNew(I->TypeFace, &fprnt, size * sampling);
                }
                if(id) {
                  if(kern_flag) {
                    line_width += (TypeFaceGetKerning(I->TypeFace,
							       last_c, c, size) / sampling);
                  }
                  line_width += CharacterGetAdvance(G, sampling, id);
                }
              }
              kern_flag = true;
              last_c = c;
            }
          }
	  if (line_widths)
	    line_widths[linenum] = line_width;
	  text_width = max2(text_width, line_width) ;
	  tot_text_width = text_width + 2.f * text_buffer[0];
	  TextSetWidth(G, tot_text_width);
	  tot_height = pymol_roundf(size * (nlines + (nlines-1) * (text_spacing-1.f)) + 2.f * text_buffer[1]);
	  TextSetHeight(G, tot_height);
	  {
	    float factor = rpos[0] / 2.0F - 0.5F;
	    /* if -1. < rpos[0] < 1. , then we need to determine the label's width
	       so that we can justify it appropriately */
	    // if rpos[0] < -1., right justified
	    // if rpos[0] > 1., left justified, label width not needed
	    factor = (factor < _m1) ? _m1 : (factor > _0) ? _0 : factor;
	    x_indent -= factor * tot_text_width;
	  }
        }
	/* if label_position x is -1 to 1, the label is placed in x such that
	   0 - centered 
	   -1 - right justified on projected point
	    1 - left justified on projected point
	*/
	if(rpos[0] < _m1) {
	  screenWorldOffset[0] -= (rpos[0] + _1);
	} else if(rpos[0] > _1) {
	  screenWorldOffset[0] -= (rpos[0] - _1);
	}
        if(rpos[1] < _1) {
          float factor = -rpos[1] / 2.0F + 0.5F;
	  factor = (factor > _1) ? _1 : (factor < _0) ? _0 : factor;
	  y_indent = pymol_roundf((size * factor) * (1.f + text_spacing * max2(0.f, nlines-1.f)));
	}
	{
          float factor = rpos[1];
	  factor = (factor > _1) ? _1 : (factor < _m1) ? _m1 : factor;
	  y_indent -= pymol_roundf(factor * text_buffer[1]);
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
	  y_indent += pymol_roundf(screenWorldOffset[1] / v_scale);
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
	if (!shaderCGO){
	  unsigned char posIsSet = TextGetLabelPosIsSet(G);
	  switch (posIsSet){
	  case 1:
	    SceneAdjustZtoScreenZ(G, v, TextGetLabelPos(G)[0]);
	    break;
	  case 2:
	    SceneSetPointToWorldScreenRelative(G, v, TextGetLabelPos(G));
	    break;
	  }
	}
	TextSetLabelPushPos(G, v);
	if (!shouldRender)
	  return st;
        ScenePushRasterMatrix(G, v);
        TextSetPos(G, zero);
#ifndef PURE_OPENGL_ES_2
      } else {
#endif
	if (!shouldRender)
	  return st;
      }
      if(rpos) {
        TextIndent(G, x_indent, y_indent);
      }
      CharacterRenderOpenGLPrime(G, info);
      kern_flag = false;
      TextGetPos(G)[1] += (int)pymol_roundf(size * (text_spacing * (nlines-1) + descender));
      if (line_widths){
	TextGetPos(G)[0] += text_just * (text_width - line_widths[0])/2.f;
      }
      TextGetPos(G)[0] += text_buffer[0];
      linenum = 0;
      cont = 1;
      while(cont && (c = *(st++))) {
	if (c == '\n'){
	  float zero[3] = { 0.0F, 0.0F, 0.0F };
	  TextSetPos(G, zero);
	  if (rpos){
	    TextIndent(G, x_indent, y_indent);
	  }
	  linenum++;
	  TextGetPos(G)[1] += (int)pymol_roundf(size * (text_spacing * (nlines - 1 - linenum) + descender));
	  if (line_widths)
	    TextGetPos(G)[0] += text_just * (text_width - line_widths[linenum])/2.f + text_buffer[0];
	  kern_flag = false;
	  continue;
	}
	CheckUnicode(&c, &unicnt, &unicode);
        if(!unicnt) {
          CharFngrprnt fprnt;
	  GenerateCharFngrprnt(G, &fprnt, c, I->TextID, size, sampling, 1, flat);
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
              cont &= CharacterRenderOpenGL(G, info, id, true, relativeMode SHADERCGOARGVAR);       /* handles advance */
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
      FreeP(line_widths);
    }
    if (!cont)
      return st;
  }
  return st;
}

const char* CFontType::RenderOpenGL(const RenderInfo* info, const char* st,
    float size, const float* rpos, bool needSize, short relativeMode,
    bool shouldRender, CGO* shaderCGO)
{
  return FontTypeRenderOpenGLImpl(info, this, st, size, false, rpos, needSize,
      relativeMode, shouldRender, shaderCGO);
}

const char* CFontType::RenderOpenGLFlat(const RenderInfo* info, const char* st,
    float size, const float* rpos, bool needSize, short relativeMode,
    bool shouldRender, CGO* shaderCGO)
{
  return FontTypeRenderOpenGLImpl(info, this, st, size, true, rpos, needSize,
      relativeMode, shouldRender, shaderCGO);
}

const char* CFontType::RenderRay(CRay* ray, const char* st, float size,
    const float* rpos, bool needSize, short relativeMode)
{
  auto I = this;
  PyMOLGlobals *G = I->G;
  unsigned int c;
  int kern_flag = false;
  unsigned int last_c = 0;
  int sampling = ray->Sampling;
  const float _0 = 0.0F, _1 = 1.0F, _m1 = -1.0F;
  float x_indent = 0.0F, y_indent = 0.0F, z_indent = 0.0F;
  float text_width = 0.f, line_width = 0.f, tot_text_width;
  float text_just = 1.f - TextGetJustification(G);
  float text_spacing = TextGetSpacing(G);
  float text_buffer[2];
  float xn[3], yn[3], x_adj[3], y_adj[3], pos[3], *v;
  int unicode = 0;
  int unicnt = 0;
  int nlines = countchrs(st, '\n') + 1;
  int linenum = 0;
  float descender = TypeFaceGetDescender(I->TypeFace);
  float v_scale = SceneGetScreenVertexScale(G, NULL);
  copy2f(TextGetLabelBuffer(G), text_buffer);
  if(st && (*st)) {
    float origpos[3];
    float *line_widths = NULL;
    float tot_height;
    if (nlines>1){
      line_widths = pymol::calloc<float>(nlines);
    }
    if(size < _0) {
      size = (int) (0.5F - size / v_scale);
    } else {
      size = DIP2PIXEL(size);
    }

    text_buffer[0] *= size;
    text_buffer[1] *= size;

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

      {
	unsigned char posIsSet = TextGetLabelPosIsSet(G);
	switch (posIsSet){
	case 1:
	  RayAdjustZtoScreenZ(ray, loc, TextGetLabelPos(G)[0]);
	  break;
	case 2:
	  RaySetPointToWorldScreenRelative(ray, loc, TextGetLabelPos(G));
	  break;
	}
      }
      TextSetLabelPushPos(G, loc);
      TextSetPos(I->G, loc);
    }

    RayGetScaledAxes(ray, xn, yn);

    if(rpos) {
      if(needSize || rpos[0] < _1) {        /* we need to measure the string width before starting to draw */
        float factor = rpos[0] / 2.0F - 0.5F;
        const char *sst = st;
	float max_x_indent = 0.f;
	// factor -1 to 0 based on justification (i.e., rpos[0])
	factor = (factor < _m1) ? _m1 : (factor > _0) ? _0 : factor;
        while((c = *(sst++))) {
	  if (c == '\n'){
	    text_width = max2(text_width, line_width);
	    line_widths[linenum] = line_width;
	    line_width = 0.f;
	    kern_flag = false;
	    linenum++;
	    max_x_indent = max2(x_indent, max_x_indent);
	    x_indent = 0.f;
	    continue;
	  }
	  CheckUnicode(&c, &unicnt, &unicode);
          if(!unicnt) {
            CharFngrprnt fprnt;
	    GenerateCharFngrprnt(G, &fprnt, c, I->TextID, size, sampling, 1, 0 /* flat is not set */);
            {
              int id = CharacterFind(G, &fprnt);
              if(!id) {
                id = TypeFaceCharacterNew(I->TypeFace, &fprnt, size * sampling);
              }
              if(id) {
		float adv;
                if(kern_flag) {
		  float kern = TypeFaceGetKerning(I->TypeFace,
						  last_c,
						  c, size * sampling) ;
		  line_width += kern;
                  x_indent -= factor * kern;
                }
		adv = CharacterGetAdvance(G, 1, id);
		line_width += adv;
                x_indent -= factor * adv;
                kern_flag = true;
                last_c = c;
              }
            }
          }
        }
	max_x_indent = max2(x_indent, max_x_indent);
	x_indent = max_x_indent;
      }
      text_width = max2(text_width, line_width);
      tot_text_width = text_width + 2.f * text_buffer[0] * sampling;
      if (line_widths)
	line_widths[linenum] = line_width;
      TextSetWidth(G, tot_text_width/(float)sampling);
      tot_height = size * (nlines + (nlines-1) * (text_spacing-1.f)) + 2.f * text_buffer[1];
      TextSetHeight(G, tot_height);
      if(rpos[0] < _m1) {
        x_indent -= 2 * (rpos[0] + _1) / v_scale;
      } else if(rpos[0] > _1) {
        x_indent -= 2 * (rpos[0] - _1) / v_scale;
      }
      if(rpos[1] < _1) {
        float factor = -rpos[1] / 2.0F + 0.5F;
        factor = (factor > _1) ? _1 : (factor < _0) ? _0 : factor;
	y_indent = (sampling * size * factor) * (1.f + text_spacing * max2(0.f, nlines-1.f));
      }
      {
	float factor = rpos[1];
	factor = (factor > _1) ? _1 : (factor < _m1) ? _m1 : factor;
	y_indent -= factor * sampling * text_buffer[1];
      }
      if(rpos[1] < _m1) {
        y_indent -= 2 * (rpos[1] + _1) / v_scale;
      } else if(rpos[1] > _1) {
        y_indent -= 2 * (rpos[1] - _1) / v_scale;
      }
      v = TextGetPos(I->G);
      if (line_widths){
	x_indent -= text_just * (text_width - line_widths[0])/2.f;
      }
      {
	float factor = rpos[0];
	factor = (factor > _1) ? _1 : (factor < _m1) ? _m1 : factor;
	x_indent -= factor * sampling * text_buffer[0];
      }
      scale3f(xn, x_indent, x_adj);
      scale3f(yn, y_indent - size * (sampling * text_spacing * (nlines-1) + descender), y_adj);
      subtract3f(v, x_adj, pos);
      subtract3f(pos, y_adj, pos);
      TextSetPos(I->G, pos);
    }
    kern_flag = false;
    copy3f(TextGetPos(I->G), origpos);
    linenum = 0;
    while((c = *(st++))) {
      if (c == '\n'){
	copy3f(origpos, TextGetPos(I->G));
	kern_flag = false;
	linenum++;
	scale3f(yn, pymol_roundf(text_spacing * size * linenum * sampling), y_adj); // need to round to pixel in y
	subtract3f(TextGetPos(G), y_adj, TextGetPos(G));
	if (line_widths){
	  scale3f(xn, text_just * (line_widths[0] - line_widths[linenum])/2.f, x_adj);
	  add3f(TextGetPos(G), x_adj, TextGetPos(G));
	}
	kern_flag = false;
	continue;
      }
      CheckUnicode(&c, &unicnt, &unicode);
      if(!unicnt) {
        CharFngrprnt fprnt;
	GenerateCharFngrprnt(G, &fprnt, c, I->TextID, size, sampling, 0, 0);
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
                                              size * sampling);
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
    FreeP(line_widths);
  }
  return st;
}

CFontType::~CFontType()
{
  if (TypeFace) {
    TypeFaceFree(TypeFace);
  }
}

CFontType::CFontType(PyMOLGlobals* G, unsigned char* dat, unsigned int len)
    : CFont(G)
{
  TypeFace = TypeFaceLoad(G, dat, len);
}

CFont* FontTypeNew(PyMOLGlobals* G, unsigned char* dat, unsigned int len)
{
  auto fontType = new CFontType(G, dat, len);
  if (!fontType->TypeFace) {
    DeleteP(fontType);
  }
  return fontType;
}
