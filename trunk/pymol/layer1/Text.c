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

struct _CText {
  int NActive;
  ActiveRec *Active;
  float Pos[4];
  float Color[4];
  int Default_ID;
};

void TextSetPosNColor(PyMOLGlobals *G,float *pos,float *color)
{
  register CText *I=G->Text;
  copy3f(pos,I->Pos);
  copy3f(color,I->Color);
  I->Pos[3]=1.0F;
  I->Color[3]=1.0F;
}
void TextAdvance(PyMOLGlobals *G,float advance)
{
  G->Text->Pos[0]+=advance;
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
}

void TextSetColor3f(PyMOLGlobals *G,float red, float green, float blue)
{
  register CText *I=G->Text;
  I->Color[0]=red;
  I->Color[1]=green;
  I->Color[2]=blue;
  I->Color[3]=1.0F;
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

const static float _255 = 255.0F;

void TextGetColorUChar(PyMOLGlobals *G,unsigned char *red,
                       unsigned char *green, 
                       unsigned char *blue,
                       unsigned char *alpha)
{
  register CText *I=G->Text;
  *red = (unsigned char)(_255*I->Color[0]);
  *green = (unsigned char)(_255*I->Color[1]);
  *blue = (unsigned char)(_255*I->Color[2]);
  *alpha = (unsigned char)(_255*I->Color[3]);
}

char *TextRenderOpenGL(PyMOLGlobals *G,int text_id,char *st)
{
  register CText *I=G->Text;
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

void TextDrawStrAt(PyMOLGlobals *G,char *st, int x, int y)
{
  register CText *I=G->Text;
  TextSetPos3f(G, (float)x, (float)y, 0.0F);
  TextRenderOpenGL(G,I->Default_ID,st);
}

void TextDrawStr(PyMOLGlobals *G,char *st)
{
  register CText *I=G->Text;
  TextRenderOpenGL(G,I->Default_ID,st);
}

void TextDrawChar(PyMOLGlobals *G,char ch)
{
  char st[2] = { ch , 0 };
  register CText *I=G->Text;
  TextRenderOpenGL(G,I->Default_ID,st);
}


char *TextRenderRay(PyMOLGlobals *G,CRay *ray,int text_id,char *st)
{
  register CText *I=G->Text;
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


int TextInit(PyMOLGlobals *G)
{
  register CText *I=NULL;
  if( (I=(G->Text=Calloc(CText,1)))) {

    I->NActive = 0;
    I->Active = VLACalloc(ActiveRec,10);
    I->Default_ID = 0;

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

    return 1;
  } else 
    return 0;

}

int TextGetFontID(PyMOLGlobals *G,int src, int code, char *name,int mode, int size, int style)
{
  /* first, return the font code if it is already active */
  register CText *I=G->Text;
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
    I->Active[I->NActive].Font = FontGLUTNew(G,code);
    if(I->Active[I->NActive].Font) {
      I->Active[I->NActive].Src = cTextSrcGLUT;
      I->Active[I->NActive].Code = code;
      I->NActive++;
    }
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


