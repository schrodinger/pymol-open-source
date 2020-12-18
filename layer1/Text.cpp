
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

#include"MemoryDebug.h"
#include"Text.h"
#include"Font.h"

#include"FontGLUT.h"
#include"FontType.h"
#include"Color.h"
#include"Vector.h"
#include"Executive.h"
#include"Picking.h"

#ifdef _PYMOL_FREETYPE
#include "FontTTF.h"
#include "FontTTF2.h"
#endif

#define FONT_NAME_MAX 255

#define TEXT_DEFAULT_SIZE 12.0F
static const float _255 = 255.0F;
static const float _499 = 0.4999F;

#define NFONTS 20

struct CText {
  float Pos[4];
  float WorldPos[4];
  float ScreenWorldOffset[3];
  float TargetPos[3];
  float LabelPushPos[3];
  float LabelPos[3];
  unsigned char LabelPosIsSet; // 1 for just z, 2 for 3f
  float TextIndentFactor[2];
  float Color[4];
  unsigned char UColor[4];
  unsigned char OutlineColor[4];
  int Default_ID = 0;
  float Height, Width;
  float Spacing, Just;
  float LabelBuf[2];

#ifdef _WEBGL
  bool XHRFetched[NFONTS] = {};
  bool XHRFailed[NFONTS] = {};
#endif

  bool Flat = false;
  bool IsPicking = false;

private:
  std::vector<std::unique_ptr<CFont>> m_fonts;

public:
  // TODO make const
  CFont* getFont(unsigned font_id) {
    if (font_id < m_fonts.size()) {
      return m_fonts[font_id].get();
    }
    return nullptr;
  }

  //! Takes ownership of pointer
  void addFont(unsigned font_id, CFont* font) {
    if (!font)
      return;
    VecCheck(m_fonts, font_id);
    m_fonts[font_id].reset(font);
    font->TextID = font_id + 1;
  }
};

static void TextUpdateUColor(CText * I)
{
  I->UColor[0] = (unsigned char) (_255 * I->Color[0] + _499);
  I->UColor[1] = (unsigned char) (_255 * I->Color[1] + _499);
  I->UColor[2] = (unsigned char) (_255 * I->Color[2] + _499);
  I->UColor[3] = (unsigned char) (_255 * I->Color[3] + _499);
}

void TextSetLabelBkgrdInfo(PyMOLGlobals * G, float label_spacing, float label_just, const float *buff){
  CText *I = G->Text;
  I->Spacing = label_spacing;
  I->Just = label_just;
  if (buff){
    I->LabelBuf[0] = buff[0];
    I->LabelBuf[1] = buff[1];
  } else {
    I->LabelBuf[0] = I->LabelBuf[1] = .2f;
  }
}

void TextSetIsPicking(PyMOLGlobals * G, bool IsPicking)
{
  CText *I = G->Text;
  I->IsPicking = IsPicking;
}

bool TextGetIsPicking(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->IsPicking;
}

void TextSetPosNColor(PyMOLGlobals * G, const float *pos, const float *color)
{
  CText *I = G->Text;
  copy3f(pos, I->Pos);
  copy3f(color, I->Color);
  I->Flat = false;
  I->Pos[3] = 1.0F;
  I->Color[3] = 1.0F;
  TextUpdateUColor(I);
}

void TextAdvance(PyMOLGlobals * G, float advance)
{
  G->Text->Pos[0] += advance;
}

void TextIndent(PyMOLGlobals * G, float x, float y)
{
  CText *I = G->Text;
  I->Pos[0] -= x;
  I->Pos[1] -= y;
}

void TextSetPos(PyMOLGlobals * G, const float *pos)
{
  CText *I = G->Text;
  copy3f(pos, I->Pos);
  I->Pos[3] = 1.0F;
}

void TextSetWorldPos(PyMOLGlobals * G, const float *pos)
{
  CText *I = G->Text;
  copy3f(pos, I->WorldPos);
  I->WorldPos[3] = 1.0F;
}
float *TextGetWorldPos(PyMOLGlobals * G){
  CText *I = G->Text;
  return I->WorldPos;
}
void TextSetLabelPos(PyMOLGlobals * G, const float *pos)
{
  CText *I = G->Text;
  copy3f(pos, I->LabelPos);
}

float *TextGetLabelPos(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return (I->LabelPos);
}

void TextSetLabelPosIsSet(PyMOLGlobals * G, unsigned char isSet)
{
  CText *I = G->Text;
  I->LabelPosIsSet = isSet;
}

unsigned char TextGetLabelPosIsSet(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->LabelPosIsSet;
}

void TextSetLabelPushPos(PyMOLGlobals * G, const float *pos)
{
  CText *I = G->Text;
  copy3f(pos, I->LabelPushPos);
}
float *TextGetLabelPushPos(PyMOLGlobals * G){
  CText *I = G->Text;
  return I->LabelPushPos;
}
void TextSetScreenWorldOffset(PyMOLGlobals * G, const float *pos)
{
  CText *I = G->Text;
  I->ScreenWorldOffset[0] = -pos[0];
  I->ScreenWorldOffset[1] = -pos[1];
  I->ScreenWorldOffset[2] = -pos[2];
}
float *TextGetScreenWorldOffset(PyMOLGlobals * G){
  CText *I = G->Text;
  return I->ScreenWorldOffset;
}
void TextSetTargetPos(PyMOLGlobals * G, const float *pos){
  CText *I = G->Text;
  copy3f(pos, I->TargetPos);
}
float *TextGetTargetPos(PyMOLGlobals * G){
  CText *I = G->Text;
  return I->TargetPos;
}

void TextDrawSubStrFast(PyMOLGlobals * G, const char *c, int x, int y, int start, int n ORTHOCGOARG)
{
  c += start;
  TextSetPos2i(G, x, y);
  if(n)
    while(*c) {
      n--;
      TextDrawChar(G, *(c++) ORTHOCGOARGVAR);
      if(n <= 0)
        break;
    }
}

void TextDrawCharRepeat(PyMOLGlobals * G, char c, int x, int y, int start, int n ORTHOCGOARG)
{
  c += start;
  TextSetPos2i(G, x, y);
  while(n) {
    n--;
    TextDrawChar(G, c ORTHOCGOARGVAR);
  }
}

void TextSetPos2i(PyMOLGlobals * G, int x, int y)
{
  CText *I = G->Text;
  I->Pos[0] = (float) x;
  I->Pos[1] = (float) y;
  I->Pos[2] = 0.0F;
  I->Pos[3] = 1.0F;
}

static void TextSetPos3f(PyMOLGlobals * G, float x, float y, float z)
{
  CText *I = G->Text;
  I->Pos[0] = x;
  I->Pos[1] = y;
  I->Pos[2] = z;
  I->Pos[3] = 1.0F;
}

void TextSetColor(PyMOLGlobals * G, const float *color)
{
  CText *I = G->Text;
  copy3f(color, I->Color);
  I->Color[3] = 1.0F;
  I->Flat = false;
  TextUpdateUColor(I);
}

void TextSetColor3f(PyMOLGlobals * G, float red, float green, float blue)
{
  CText *I = G->Text;
  I->Flat = false;
  I->Color[0] = red;
  I->Color[1] = green;
  I->Color[2] = blue;
  I->Color[3] = 1.0F;
  TextUpdateUColor(I);
}

void TextSetOutlineColor(PyMOLGlobals * G, int color)
{
  CText *I = G->Text;
  if(color >= 0) {
    const float *fcolor = ColorGet(G, color);
    I->OutlineColor[0] = (unsigned char) (_255 * fcolor[0]);
    I->OutlineColor[1] = (unsigned char) (_255 * fcolor[1]);
    I->OutlineColor[2] = (unsigned char) (_255 * fcolor[2]);
    I->OutlineColor[3] = 0xFF;
  } else {
    I->OutlineColor[3] = 0;
  }
}

static const float _inv255 = 1.0F / 255.0F;

void TextSetColorFromUColor(PyMOLGlobals * G)
{
  CText *I = G->Text;
  I->Color[0] = I->UColor[0] * _inv255;
  I->Color[1] = I->UColor[1] * _inv255;
  I->Color[2] = I->UColor[2] * _inv255;
  I->Color[3] = 1.0F;
}

float *TextGetPos(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->Pos;
}

float TextGetWidth(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->Width;
}

float TextGetHeight(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->Height;
}

void TextSetIndentFactorX(PyMOLGlobals * G, float factor){
  CText *I = G->Text;
  I->TextIndentFactor[0] = factor;
}

void TextSetIndentFactorY(PyMOLGlobals * G, float factor){
  CText *I = G->Text;
  I->TextIndentFactor[1] = factor;
}
float *TextGetIndentFactor(PyMOLGlobals * G){
  CText *I = G->Text;
  return I->TextIndentFactor;
}

void TextSetWidth(PyMOLGlobals * G, float text_width)
{
  CText *I = G->Text;
  I->Width = text_width;
}

void TextSetHeight(PyMOLGlobals * G, float text_height)
{
  CText *I = G->Text;
  I->Height = text_height;
}

float *TextGetColor(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->Color;
}

unsigned char *TextGetColorUChar4uv(PyMOLGlobals * G){
  CText *I = G->Text;
  return I->UColor;
}

void TextGetColorUChar(PyMOLGlobals * G, unsigned char *red,
                       unsigned char *green, unsigned char *blue, unsigned char *alpha)
{
  CText *I = G->Text;
  *red = I->UColor[0];
  *green = I->UColor[1];
  *blue = I->UColor[2];
  *alpha = I->UColor[3];
}

void TextGetOutlineColor(PyMOLGlobals * G,
                         unsigned char *red,
                         unsigned char *green, unsigned char *blue, unsigned char *alpha)
{
  CText *I = G->Text;
  *red = I->OutlineColor[0];
  *green = I->OutlineColor[1];
  *blue = I->OutlineColor[2];
  *alpha = I->OutlineColor[3];
}

const char *TextRenderOpenGL(PyMOLGlobals * G, const RenderInfo * info, int text_id,
    const char *st, float size, const float *rpos,
    short needSize, short relativeMode, short shouldRender,
    CGO *shaderCGO)
{
  CText *I = G->Text;

  if(st && (*st)) {
    auto font = I->getFont(text_id);

    if (font) {
        if (I->Flat) {
          return font->RenderOpenGLFlat(info, st, size, rpos, needSize,
              relativeMode, shouldRender, shaderCGO);
        } else {
          return font->RenderOpenGL(info, st, size, rpos, needSize,
              relativeMode, shouldRender, shaderCGO);
        }
    }
    /* make sure we got to end of string */
    if(*st)
      while(*(st++));
  }
  return st;
}

void TextDrawStrAt(PyMOLGlobals * G, const char *st, int x, int y ORTHOCGOARG)
{
  CText *I = G->Text;
  TextSetPos3f(G, (float) x, (float) y, 0.0F);
  TextRenderOpenGL(G, NULL, I->Default_ID, st, TEXT_DEFAULT_SIZE, NULL, false, 0, 1 ORTHOCGOARGVAR);
}

void TextDrawStr(PyMOLGlobals * G, const char *st ORTHOCGOARG)
{
  CText *I = G->Text;
  TextRenderOpenGL(G, NULL, I->Default_ID, st, TEXT_DEFAULT_SIZE, NULL, false, 0, 1 ORTHOCGOARGVAR);
}

void TextDrawChar(PyMOLGlobals * G, char ch ORTHOCGOARG)
{
  char st[2] = { 0, 0 };
  CText *I = G->Text;
  st[0] = ch;
  TextRenderOpenGL(G, NULL, I->Default_ID, st, TEXT_DEFAULT_SIZE, NULL, false, 0, 1 ORTHOCGOARGVAR);
}

const char *TextRenderRay(PyMOLGlobals * G, CRay * ray, int text_id,
    const char *st, float size, const float *rpos, short needSize, short relativeMode)
{
  CText *I = G->Text;

  if(st && (*st)) {
    auto font = I->getFont(text_id);
    if (font) {
      if(size >= 0.0F)
        size *= ray->Magnified;

      return font->RenderRay(ray, st, size, rpos, needSize, relativeMode);
    }
    /* make sure we got to end of string */
    if(*st)
      while(*(st++));
  }
  return st;
}

int TextInit(PyMOLGlobals * G)
{
  assert(!G->Text);
  G->Text = new CText();
  auto I = G->Text;

  I->addFont(0, new CFontGLUT(G, &FontGLUTBitmap8By13));
  I->addFont(1, new CFontGLUT(G, &FontGLUTBitmap9By15));
  I->addFont(2, new CFontGLUT(G, &FontGLUTBitmapHelvetica10));
  I->addFont(3, new CFontGLUT(G, &FontGLUTBitmapHelvetica12));
  I->addFont(4, new CFontGLUT(G, &FontGLUTBitmapHelvetica18));

#ifdef _PYMOL_FREETYPE
#if !defined(_WEBGL) || defined(_WEBGL_INCLUDE_DEFAULT_FONT)
  I->addFont(5, FontTypeNew(G, TTF_DejaVuSans_dat, TTF_DejaVuSans_len));
#endif
#ifndef _WEBGL
  I->addFont(6, FontTypeNew(G, TTF_DejaVuSans_Oblique_dat, TTF_DejaVuSans_Oblique_len));
  I->addFont(7, FontTypeNew(G, TTF_DejaVuSans_Bold_dat, TTF_DejaVuSans_Bold_len));
  I->addFont(8, FontTypeNew(G, TTF_DejaVuSans_BoldOblique_dat, TTF_DejaVuSans_BoldOblique_len));
  I->addFont(9, FontTypeNew(G, TTF_DejaVuSerif_dat, TTF_DejaVuSerif_len));
  I->addFont(10, FontTypeNew(G, TTF_DejaVuSerif_Bold_dat, TTF_DejaVuSerif_Bold_len));
  I->addFont(11, FontTypeNew(G, TTF_DejaVuSansMono_dat, TTF_DejaVuSansMono_len));
  I->addFont(12, FontTypeNew(G, TTF_DejaVuSansMono_Oblique_dat, TTF_DejaVuSansMono_Oblique_len));
  I->addFont(13, FontTypeNew(G, TTF_DejaVuSansMono_Bold_dat, TTF_DejaVuSansMono_Bold_len));
  I->addFont(14, FontTypeNew(G, TTF_DejaVuSansMono_BoldOblique_dat, TTF_DejaVuSansMono_BoldOblique_len));
  I->addFont(15, FontTypeNew(G, TTF_GenR102_dat, TTF_GenR102_len));
  I->addFont(16, FontTypeNew(G, TTF_GenI102_dat, TTF_GenI102_len));
  I->addFont(17, FontTypeNew(G, TTF_DejaVuSerif_Oblique_dat, TTF_DejaVuSerif_Oblique_len));
  I->addFont(18, FontTypeNew(G, TTF_DejaVuSerif_BoldOblique_dat, TTF_DejaVuSerif_BoldOblique_len));
#endif
#endif

  return true;
}

void TextFree(PyMOLGlobals * G)
{
  DeleteP(G->Text);
}
float TextGetSpacing(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->Spacing;
}
float TextGetJustification(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->Just;
}
float *TextGetLabelBuffer(PyMOLGlobals * G)
{
  CText *I = G->Text;
  return I->LabelBuf;
}

/**
 * GUI elements like internal menus or the wizard prompt can handle text
 * color markup in the form "\\RGB" where RGB are three digits (0-9) or
 * "---" to reset the color.
 *
 * Return true if `p` starts with "\\RGB" or "\\---".
 */
bool TextStartsWithColorCode(const char *p)
{
  if (p[0] != '\\') {
    return false;
  }

  if (p[1] == '-') {
    return p[2] == '-' && p[3] == '-';
  }

  return (
      ('0' <= p[1] && p[1] <= '9') &&
      ('0' <= p[2] && p[2] <= '9') &&
      ('0' <= p[3] && p[3] <= '9'));
}

/**
 * Set text color from "\\RGB" code.
 *
 * "\\---" -> defaultcolor
 *
 * Return false if `p` does not start with a color code.
 */
bool TextSetColorFromCode(PyMOLGlobals * G,
    const char *p,
    const float *defaultcolor)
{
  if (!TextStartsWithColorCode(p)) {
    return false;
  }

  if (p[1] == '-') {
    TextSetColor(G, defaultcolor);
  } else {
    TextSetColor3f(G,
        (p[1] - '0') / 9.0F,
        (p[2] - '0') / 9.0F,
        (p[3] - '0') / 9.0F);
  }

  return true;
}
