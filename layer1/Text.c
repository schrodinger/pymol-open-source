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

#include"MemoryDebug.h"
#include"Text.h"
#include"Font.h"

#include"FontGLUT.h"
#include"FontType.h"
#include"Color.h"
#include"Vector.h"

#ifdef _PYMOL_FREETYPE
#include "FontTTF.h"
#include "FontTTF2.h"
#endif

#define FONT_NAME_MAX 255

#define TEXT_DEFAULT_SIZE 12.0F
static const float _255 = 255.0F;
static const float _499 = 0.4999F;

typedef struct {
  int Src;
  int Code;
  char Name[FONT_NAME_MAX];
  int Mode;
  int Style;
  CFont *Font;
} ActiveRec;

struct _CText {
  int NActive;
  ActiveRec *Active;
  float Pos[4];
  float Color[4];
  unsigned char UColor[4];
  unsigned char OutlineColor[4];
  int Default_ID;
  int Flat;
};

static void TextUpdateUColor(CText *I)
{
  I->UColor[0] = (unsigned char)(_255*I->Color[0]+_499);
  I->UColor[1] = (unsigned char)(_255*I->Color[1]+_499);
  I->UColor[2] = (unsigned char)(_255*I->Color[2]+_499);
  I->UColor[3] = (unsigned char)(_255*I->Color[3]+_499);
}

void TextSetPosNColor(PyMOLGlobals *G,float *pos,float *color)
{
  register CText *I=G->Text;
  copy3f(pos,I->Pos);
  copy3f(color,I->Color);
  I->Flat = false;
  I->Pos[3]=1.0F;
  I->Color[3]=1.0F;
  TextUpdateUColor(I);
}
void TextAdvance(PyMOLGlobals *G,float advance)
{
  G->Text->Pos[0]+=advance;
}

void TextSetLabPos(PyMOLGlobals *G,float *pos, LabPosType *labpos, char *text)
{
  if((!labpos)||(!labpos->mode))
    TextSetPos(G,pos);
  else {
    register CText *I=G->Text;
    switch(labpos->mode) {
    default:
      copy3f(pos,I->Pos);
      add3f(labpos->offset, I->Pos, I->Pos);
      break;
    }
  }
}

void TextIndent(PyMOLGlobals *G,float x,float y)
{
  register CText *I=G->Text;
  I->Pos[0]-=x;
  I->Pos[1]-=y;
}

void TextSetPos(PyMOLGlobals *G,float *pos)
{
  register CText *I=G->Text;
  copy3f(pos,I->Pos);
  I->Pos[3]=1.0F;
}

void TextDrawSubStrFast(PyMOLGlobals *G,char *c,int x,int y,int start,int n)
{
  c+=start;
  TextSetPos2i(G,x,y);
  if(n)
    while(*c) {
      n--;
      TextDrawChar(G,*(c++));
      if(n<=0) break;
    }
}
void TextDrawCharRepeat(PyMOLGlobals *G,char c,int x,int y,int start,int n)
{
  c+=start;
  TextSetPos2i(G,x,y);
  while(n) {
    n--;
    TextDrawChar(G,c);
  }
}

void TextSetPos2i(PyMOLGlobals *G,int x,int y)
{
  register CText *I=G->Text;
  I->Pos[0]=(float)x;
  I->Pos[1]=(float)y;
  I->Pos[2]=0.0F;
  I->Pos[3]=1.0F;
}

static void TextSetPos3f(PyMOLGlobals *G,float x,float y, float z)
{
  register CText *I=G->Text;
  I->Pos[0]=x;
  I->Pos[1]=y;
  I->Pos[2]=z;
  I->Pos[3]=1.0F;
}


void TextSetColor(PyMOLGlobals *G,float *color)
{
  register CText *I=G->Text;
  copy3f(color,I->Color);
  I->Color[3]=1.0F;
  I->Flat = false;
  TextUpdateUColor(I);
}

void TextSetColor3f(PyMOLGlobals *G,float red, float green, float blue)
{
  register CText *I=G->Text;
  I->Flat = false;
  I->Color[0]=red;
  I->Color[1]=green;
  I->Color[2]=blue;
  I->Color[3]=1.0F;
  TextUpdateUColor(I);
}

void TextSetOutlineColor(PyMOLGlobals *G,int color)
{
  register CText *I=G->Text;
  if(color>=0) {
    float *fcolor = ColorGet(G,color);
    I->OutlineColor[0] = (unsigned char)(_255*fcolor[0]);
    I->OutlineColor[1] = (unsigned char)(_255*fcolor[1]);
    I->OutlineColor[2] = (unsigned char)(_255*fcolor[2]);
    I->OutlineColor[3] = 0xFF;
  } else {
    I->OutlineColor[3] = 0;
  }
}

static const float _inv255 = 1.0F/255.0F;

void TextSetPickColor(PyMOLGlobals *G,int first_pass, int index)
{
  register CText *I=G->Text;
  if(!first_pass)
    index = (index>>12); /* high order bits */

  I->Flat = true;
  I->UColor[0] = ((unsigned char)((index&0xF)<<4));
  I->UColor[1] = ((unsigned char)((index&0xF0)|0x8));
  I->UColor[2] = ((unsigned char)((index&0xF00)>>4));
  I->UColor[3] = 0xFF;

  I->Color[0] = I->UColor[0] * _inv255;
  I->Color[1] = I->UColor[1] * _inv255;
  I->Color[2] = I->UColor[2] * _inv255;
  I->Color[3] = 1.0F;
}

float *TextGetPos(PyMOLGlobals *G)
{
  register CText *I=G->Text;
  return I->Pos;
}

float *TextGetColor(PyMOLGlobals *G)
{
  register CText *I=G->Text;
  return I->Color;
}

void TextGetColorUChar(PyMOLGlobals *G,unsigned char *red,
                       unsigned char *green, 
                       unsigned char *blue,
                       unsigned char *alpha)
{
  register CText *I=G->Text;
  *red = I->UColor[0];
  *green = I->UColor[1];
  *blue = I->UColor[2];
  *alpha = I->UColor[3];
}

void TextGetOutlineColor(PyMOLGlobals *G,
                         unsigned char *red,
                         unsigned char *green, 
                         unsigned char *blue,
                         unsigned char *alpha)
{
  register CText *I=G->Text;
  *red = I->OutlineColor[0];
  *green = I->OutlineColor[1];
  *blue = I->OutlineColor[2];
  *alpha = I->OutlineColor[3];
}

char *TextRenderOpenGL(PyMOLGlobals *G,RenderInfo *info,int text_id,
                       char *st,float size, float *rpos)
{
  register CText *I=G->Text;
  CFont *font;
  FontRenderOpenGLFn *fn;
  if((text_id<0)||(text_id>=I->NActive)) 
    text_id=0;
      
  if(st&&(*st)) {
    if((text_id>=0)&&(text_id<I->NActive)) {
      font = I->Active[text_id].Font;
      if(I->Flat) 
        fn = font->fRenderOpenGLFlat;
      else
        fn = font->fRenderOpenGL;
      if(fn)
        return fn(info,font,st,size,rpos);
    }
    /* make sure we got to end of string */
    if(*st) while(*(st++)); 
  }
  return st;
}

void TextDrawStrAt(PyMOLGlobals *G,char *st, int x, int y)
{
  register CText *I=G->Text;
  TextSetPos3f(G, (float)x, (float)y, 0.0F);
  TextRenderOpenGL(G,NULL,I->Default_ID,st,TEXT_DEFAULT_SIZE,NULL);
}

void TextDrawStr(PyMOLGlobals *G,char *st)
{
  register CText *I=G->Text;
  TextRenderOpenGL(G,NULL,I->Default_ID,st,TEXT_DEFAULT_SIZE,NULL);
}

void TextDrawChar(PyMOLGlobals *G,char ch)
{
  char st[2] = { 0 , 0 };
  register CText *I=G->Text;
  st[0] = ch;
  TextRenderOpenGL(G,NULL,I->Default_ID,st,TEXT_DEFAULT_SIZE,NULL);
}


char *TextRenderRay(PyMOLGlobals *G,CRay *ray,int text_id,
                    char *st,float size, float *rpos)
{
  register CText *I=G->Text;
  CFont *font;
  FontRenderRayFn *fn;

  if((text_id<0)||(text_id>=I->NActive)) 
    text_id=0;

  if(st&&(*st)) {
    if((text_id>=0)&&(text_id<I->NActive)) {
      font = I->Active[text_id].Font;
      if(size>=0.0F)
        size*=ray->Magnified;
        
      fn = font->fRenderRay;
      if(fn)
        return fn(ray,font,st,size,rpos);
    }
    /* make sure we got to end of string */
    if(*st) while(*(st++)); 
  }
  return st;
}


int TextInit(PyMOLGlobals *G)
{
  register CText *I=NULL;
  if( (I=(G->Text=Calloc(CText,1)))) {

    I->NActive = 0;
    I->Active = VLACalloc(ActiveRec,10);
    I->Default_ID = 0;
    I->Flat = false;

    /* font 0 is old reliable GLUT 8x13 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontGLUTNew(G,cFontGLUT8x13);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = cFontGLUT8x13;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* font 1 is GLUT 9x15 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontGLUTNew(G,cFontGLUT9x15);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = cFontGLUT9x15;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* font 2 is GLUT Helvetica10 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontGLUTNew(G,cFontGLUTHel10);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = cFontGLUTHel10;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* font 3 is GLUT Helvetica12 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontGLUTNew(G,cFontGLUTHel12);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = cFontGLUTHel12;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* font 4 is GLUT Helvetica18 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontGLUTNew(G,cFontGLUTHel18);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = cFontGLUTHel18;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

#ifdef _PYMOL_FREETYPE

    /* 5 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSans_dat, TTF_DejaVuSans_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 6 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSans_Oblique_dat, TTF_DejaVuSans_Oblique_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 7 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSans_Bold_dat, TTF_DejaVuSans_Bold_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 8 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSans_BoldOblique_dat, TTF_DejaVuSans_BoldOblique_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 9 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSerif_dat, TTF_DejaVuSerif_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 10 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSerif_Bold_dat, TTF_DejaVuSerif_Bold_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 11 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSansMono_dat, TTF_DejaVuSansMono_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 12 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSansMono_Oblique_dat, TTF_DejaVuSansMono_Oblique_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 13 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSansMono_Bold_dat, TTF_DejaVuSansMono_Bold_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 14 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSansMono_BoldOblique_dat, TTF_DejaVuSansMono_BoldOblique_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* Gentium */

    /* 15 */
    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_GenR102_dat, TTF_GenR102_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 16 */
    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_GenI102_dat, TTF_GenI102_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* back to DejaVu for the last two */

    /* 17 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSerif_Oblique_dat, TTF_DejaVuSerif_Oblique_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    /* 18 */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_DejaVuSerif_BoldOblique_dat, TTF_DejaVuSerif_BoldOblique_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

#if 0

    /* Linux Libertine not inluded due to already excessive bloat of Text.o */

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_LinLibertine_R_2_2_0_dat, TTF_LinLibertine_R_2_2_0_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_LinLibertine_It_2_2_0rc1_dat, TTF_LinLibertine_It_2_2_0rc1_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_LinLibertine_Bd_2_2_0rc10_dat, TTF_LinLibertine_Bd_2_2_0rc10_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontTypeNew(G,TTF_LinLibertine_BdIt_2_1_6v2_dat, TTF_LinLibertine_BdIt_2_1_6v2_len);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcFreeType;
      I->Active[I->NActive].Font->TextID = I->NActive;
      I->NActive++;
    }

#endif

#endif

    return 1;
  } else 
    return 0;

}

int TextGetFontID(PyMOLGlobals *G,int src, int code, char *name,int mode, int style)
{
  /* first, return the font code if it is already active */
  register CText *I=G->Text;
  {
    int a;
    ActiveRec *rec = I->Active;
    for(a=0;I->NActive;a++) {
      if((src == rec->Src) &&
         (code == rec->Code) &&
         (mode == rec->Mode) &&
         (style == rec->Style))
        if(((!name)&&(!rec->Name[0])) ||
           ( name &&( strcmp(name,rec->Name) == 0 ))) {
          return a;
        }
      rec++;
    }
  }

  switch(src) {
  case cTextSrcGLUT:
    VLACheck(I->Active,ActiveRec,I->NActive);
    I->Active[I->NActive].Font = FontGLUTNew(G,code);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = code;
      I->NActive++;
    }
    break;
  case cTextSrcFreeType:
    
    break;
  }
  return -1;
}

void TextFree(PyMOLGlobals *G)
{
  register CText *I=G->Text;
  int a;
  CFont *fp;
  for(a=0;a<I->NActive;a++) {
    fp = I->Active[a].Font;
    if(fp && fp->fFree)
      fp->fFree(fp);
  }
  VLAFreeP(I->Active);
  FreeP(G->Text);
}


