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

#include"os_gl.h"
#include"os_python.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Feedback.h"
#include"P.h"
#include"PConv.h"
#include"VFont.h"

#define VFONT_MASK 0xFF

typedef struct {
  int face;
  float size;
  int style;
  int offset[VFONT_MASK+1];
  float advance[VFONT_MASK+1];
  float *pen;
} VFontRec;

typedef struct  {
  VFontRec **Font;
  int NFont;
} CVFont;

static CVFont VFont;

VFontRec *VFontRecNew(void);
int VFontRecLoad(VFontRec *I,PyObject *dict);
void VFontRecFree(VFontRec *I);

VFontRec *VFontRecNew(void)
{
  int a;
  OOAlloc(VFontRec);
  for(a=0;a<=VFONT_MASK;a++) {
    I->advance[a]=0.0F;
    I->offset[a]=-1;
  }
  I->pen = VLAlloc(float,1000);
  return(I);
}

int VFontWriteToCGO(int font_id,CGO *cgo,char *text,float *pos,float *scale,float *matrix)
{
  CVFont *I=&VFont;
  VFontRec *fr = NULL;
  int ok=true;
  float base[3],pen[3];
  float *pc;
  unsigned char c;
  int drawing,stroke;
  int offset;

  if((font_id>0)&&(font_id<=I->NFont)) {
    fr = I->Font[font_id];
    if(fr) 
      while(1) {
        c = *(text++);
        if(!c) break;
        offset = fr->offset[c];
        if(offset>=0) {
          pc = fr->pen + offset;
          copy3f(pos,base);
          drawing = true; 
          stroke = false;
          while(drawing) {
            switch((int)*(pc++)) {
            case -1: /* sentinel */
              drawing=false;
              break;
            case 0: /* moveto */
              pen[0] = (*(pc++))*scale[0];
              pen[1] = (*(pc++))*scale[1];
              pen[2] = 0.0;
              if(matrix)
                transform33f3f(matrix,pen,pen);
              add3f(base,pen,pen);
              if(stroke) {
                CGOEnd(cgo);
              }
              CGOBegin(cgo,GL_LINE_STRIP);
              CGOVertexv(cgo,pen);
              stroke=true;
              break;
            case 1: /* drawto */
              pen[0] = (*(pc++))*scale[0];
              pen[1] = (*(pc++))*scale[1];
              pen[2] = 0.0;
              if(matrix)
                transform33f3f(matrix,pen,pen);
              add3f(base,pen,pen);
              if(stroke) {
                CGOVertexv(cgo,pen);
              }
              break;
            default:
              drawing=false;
              break;
            }
          }
          pen[0] = fr->advance[c]*scale[0];
          pen[1] = 0.0;
          pen[2] = 0.0;
          if(matrix)
            transform33f3f(matrix,pen,pen);            
          add3f(pen,pos,pos);
          if(stroke)
            CGOEnd(cgo);
        }
      }
  } else {
    PRINTFB(FB_VFont,FB_Errors) 
      "VFontWriteToCGO-Error: invalid font identifier\n"
      ENDFB;
    ok=false;
  }
  return(ok);
}

int VFontIndent(int font_id,char *text,float *pos,float *scale,float *matrix,float dir)
{
  CVFont *I=&VFont;
  VFontRec *fr = NULL;
  int ok=true;
  float base[3],pen[3];
  float *pc;
  unsigned char c;
  int offset;

  if((font_id>0)&&(font_id<=I->NFont)) {
    fr = I->Font[font_id];
    if(fr) 
      while(1) {
        c = *(text++);
        if(!c) break;
        offset = fr->offset[c];
        if(offset>=0) {
          pc = fr->pen + offset;
          copy3f(pos,base);
          pen[0] = fr->advance[c]*scale[0]*dir;
          pen[1] = 0.0;
          pen[2] = 0.0;
          if(matrix)
            transform33f3f(matrix,pen,pen);            
          add3f(pen,pos,pos);
        }
      }
  } else {
    PRINTFB(FB_VFont,FB_Errors) 
      "VFontWriteToCGO-Error: invalid font identifier\n"
      ENDFB;
    ok=false;
  }
  return(ok);
}



int VFontRecLoad(VFontRec *I,PyObject *dict)
{ /* assumes blocked Python interpreter */
  int used=0;
  int ok=true;
  PyObject *key,*char_list;
  PyObject *stroke_list = NULL;
  int pos = 0;
  unsigned char code[2];
  float adv;
  int n_float;
  while (PyDict_Next(dict, &pos, &key, &char_list)) {
    if(!PConvPyStrToStr(key,code,1)) {
      PRINTFB(FB_VFont,FB_Errors) 
        "VFont-Error: Bad character code."
        ENDFB;
      ok=false;
    } else {
      if(ok) ok = (char_list!=NULL);
      if(ok) ok = PyList_Check(char_list);
      if(ok) ok = (PyList_Size(char_list)>=2);
      if(ok) ok = PConvPyObjectToFloat(PyList_GetItem(char_list,0),&adv);
      if(ok) {
        stroke_list = PyList_GetItem(char_list,1);
        if(ok) ok = (stroke_list!=NULL);
        if(ok) ok = PyList_Check(stroke_list);
        if(ok) {
          n_float = PyList_Size(stroke_list);
          VLACheck(I->pen,float,n_float+used+1);
          ok = PConvPyListToFloatArrayInPlace(stroke_list,I->pen+used,n_float);
          I->offset[code[0]] = used;
          I->advance[code[0]] = adv;
          I->pen[used+n_float] = -1.0F; /* sentinel */
          PRINTFD(FB_VFont)
            " VFontRecLoad-Debug: Added '%c' adv: %0.3f n_float: %d\n",code[0],adv,n_float
            ENDFD;
          if(ok) used+=n_float+1;

        }
      }
    }
  }
  return(ok);
}

void VFontRecFree(VFontRec *I)
{
  VLAFreeP(I->pen);
  OOFreeP(I);
}

void VFontInit(void)
{
  CVFont *I=&VFont;
  I->Font=VLAlloc(VFontRec*,10);
  I->NFont = 0;
}

void VFontFree(void)
{
  CVFont *I=&VFont;
  int a;
  for(a=1;a<=I->NFont;a++) {
    VFontRecFree(I->Font[a]);
  }
  VLAFreeP(I->Font);
}

int VFontLoad(float size,int face,int style,int can_load_new)
{ 
  CVFont *I=&VFont;
  VFontRec *fr;
  PyObject *vfont = NULL;
  int a;
  int result = 0;

  PRINTFD(FB_VFont)
    " VFontLoad-Debug: Entered %f %d %d\n",size,face,style
    ENDFD;

  for(a=1;a<=I->NFont;a++) {
    fr = I->Font[a];
    if((fr->size==size)&&
       (fr->face==face)&&
       (fr->style==style)) {
      result=a;
      break;
    }
  }
  if(!result) {
    if(can_load_new) {
      vfont = PGetFontDict(size,face,style);
      if(vfont) {
        if(PyDict_Check(vfont)) {
          VLACheck(I->Font,VFontRec*,I->NFont+1);
          fr = VFontRecNew();
          if(!VFontRecLoad(fr,vfont))
            VFontRecFree(fr);
          else {
            I->NFont++; /* always start at 1 */
            I->Font[I->NFont]=fr;
            result = I->NFont;
            fr->size = size;
            fr->face = face;
            fr->style = style;
          }
        }
        Py_DECREF(vfont);
      }
    }
  }
  PRINTFD(FB_VFont)
    " VFontLoad-Debug: Leaving with result %d  (0 = failure)\n",result
    ENDFD;
  return(result);
}
