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
#include"Vector.h"

#define FONT_NAME_MAX 255

typedef struct {
  int Src;
  int Code;
  char Name[FONT_NAME_MAX];
  int Size;
  int Mode;
  int Style;
  CFont *Font;
} ActiveRec;

typedef struct {
  int NActive;
  ActiveRec *Active;
  float Pos[4];
  float Color[4];

} CText;

CText Text; /* global singleton */

void TextSetPosNColor(float *pos,float *color)
{
  CText *I=&Text;
  copy3f(pos,I->Pos);
  copy3f(color,I->Color);
  I->Pos[3]=1.0F;
  I->Color[3]=1.0F;
}
void TextSetPos(float *pos)
{
  CText *I=&Text;
  copy3f(pos,I->Pos);
  I->Pos[3]=1.0F;
}
void TextSetColor(float *color)
{
  CText *I=&Text;
  copy3f(color,I->Color);
  I->Color[3]=1.0F;
}

float *TextGetPos(void)
{
  CText *I=&Text;
  return I->Pos;
}

float *TextGetColor(void)
{
  CText *I=&Text;
  return I->Color;
}

const float _255 = 255.0F;
void TextGetColorUChar(unsigned char *red,
                       unsigned char *green, 
                       unsigned char *blue,
                       unsigned char *alpha)
{
  CText *I=&Text;
  *red = (unsigned char)(_255*I->Color[0]);
  *green = (unsigned char)(_255*I->Color[1]);
  *blue = (unsigned char)(_255*I->Color[2]);
  *alpha = (unsigned char)(_255*I->Color[3]);
}

char *TextRenderOpenGL(int text_id,char *st)
{
  CText *I=&Text;
  CFont *font;
  FontRenderOpenGLFn *fn;

  if(st&&(*st)) {
    if((text_id>=0)&&(text_id<I->NActive)) {
      font = I->Active[text_id].Font;
      fn = font->fRenderOpenGL;
      if(fn)
        return fn(font,st);
    }
    /* make sure we got to end of string */
    while(*(st++)); 
  }
  return st;
}

char *TextRenderRay(struct CRay *ray,int text_id,char *st)
{
  CText *I=&Text;
  CFont *font;
  FontRenderRayFn *fn;

  if(st&&(*st)) {
    if((text_id>=0)&&(text_id<I->NActive)) {
      font = I->Active[text_id].Font;
      fn = font->fRenderRay;
      if(fn)
        return fn(ray,font,st);
    }
    /* make sure we got to end of string */
    while(*(st++)); 
  }
  return st;
}


int TextInit(void)
{
  CText *I=&Text;
  I->NActive = 0;
  I->Active = VLACalloc(ActiveRec,10);
  
  /* font 0 is old reliable GLUT 8x13 */

  VLACheck(I->Active,ActiveRec,I->NActive);
  I->Active[I->NActive].Font = FontGLUTNew(cFontGLUT8x13);
  if(I->Active[I->NActive].Font) {
    I->Active[I->NActive].Src = cTextSrcGLUT;
    I->Active[I->NActive].Code = cFontGLUT8x13;
    I->Active[I->NActive].Font->TextID = I->NActive;
    I->NActive++;
  }

  /* font 1 is GLUT 9x15 */

  VLACheck(I->Active,ActiveRec,I->NActive);
  I->Active[I->NActive].Font = FontGLUTNew(cFontGLUT9x15);
  if(I->Active[I->NActive].Font) {
    I->Active[I->NActive].Src = cTextSrcGLUT;
    I->Active[I->NActive].Code = cFontGLUT9x15;
    I->Active[I->NActive].Font->TextID = I->NActive;
    I->NActive++;
  }

  /* font 2 is GLUT Helvetica10 */

  VLACheck(I->Active,ActiveRec,I->NActive);
  I->Active[I->NActive].Font = FontGLUTNew(cFontGLUTHel10);
  if(I->Active[I->NActive].Font) {
    I->Active[I->NActive].Src = cTextSrcGLUT;
    I->Active[I->NActive].Code = cFontGLUTHel10;
    I->Active[I->NActive].Font->TextID = I->NActive;
    I->NActive++;
  }

  /* font 3 is GLUT Helvetica12 */

  VLACheck(I->Active,ActiveRec,I->NActive);
  I->Active[I->NActive].Font = FontGLUTNew(cFontGLUTHel12);
  if(I->Active[I->NActive].Font) {
    I->Active[I->NActive].Src = cTextSrcGLUT;
    I->Active[I->NActive].Code = cFontGLUTHel12;
    I->Active[I->NActive].Font->TextID = I->NActive;
    I->NActive++;
  }

  /* font 4 is GLUT Helvetica18 */

  VLACheck(I->Active,ActiveRec,I->NActive);
  I->Active[I->NActive].Font = FontGLUTNew(cFontGLUTHel18);
  if(I->Active[I->NActive].Font) {
    I->Active[I->NActive].Src = cTextSrcGLUT;
    I->Active[I->NActive].Code = cFontGLUTHel18;
    I->Active[I->NActive].Font->TextID = I->NActive;
    I->NActive++;
  }

  return 1;
}

int TextGetFontID(int src, int code, char *name,int mode, int size, int style)
{
  /* first, return the font code if it is already active */
  CText *I=&Text;
  {
    int a;
    ActiveRec *rec = I->Active;
    for(a=0;I->NActive;a++) {
      if((src == rec->Src) &&
         (code == rec->Code) &&
         (mode == rec->Mode)&&
         (size == rec->Size) &&
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
    I->Active[I->NActive].Font = FontGLUTNew(code);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = code;
      I->NActive++;
    }
    break;
  }
  return -1;
}

void TextFree(void)
{
  CText *I=&Text;
  int a;
  CFont *fp;
  for(a=0;a<I->NActive;a++) {
    fp = I->Active[a].Font;
    if(fp && fp->fFree)
      fp->fFree(fp);
  }
  VLAFreeP(I->Active);
}


