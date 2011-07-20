
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

#include"os_gl.h"
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
  ov_diff offset[VFONT_MASK + 1];
  float advance[VFONT_MASK + 1];
  float *pen;
} VFontRec;

struct _CVFont {
  VFontRec **Font;
  int NFont;
};

//#ifndef _PYMOL_NOPY
static VFontRec *VFontRecNew(PyMOLGlobals * G)
{
  int a;
  OOAlloc(G, VFontRec);
  for(a = 0; a <= VFONT_MASK; a++) {
    I->advance[a] = 0.0F;
    I->offset[a] = -1;
  }
  I->pen = VLAlloc(float, 1000);
  return (I);
}
//#endif

int VFontWriteToCGO(PyMOLGlobals * G, int font_id, CGO * cgo,
                    char *text, float *pos, float *scale, float *matrix)
{
  register CVFont *I = G->VFont;
  VFontRec *fr = NULL;
  int ok = true;
  float base[3], pen[3];
  float *pc;
  unsigned char c;
  int drawing, stroke;
  ov_diff offset;
  if((font_id > 0) && (font_id <= I->NFont)) {
    fr = I->Font[font_id];
    if(fr)
      while(1) {
        c = *(text++);
        if(!c)
          break;
        offset = fr->offset[c];
        if(offset >= 0) {
          pc = fr->pen + offset;
          copy3f(pos, base);
          drawing = true;
          stroke = false;
          while(drawing) {
            switch ((int) *(pc++)) {
            case -1:           /* sentinel */
              drawing = false;
              break;
            case 0:            /* moveto */
              pen[0] = (*(pc++)) * scale[0];
              pen[1] = (*(pc++)) * scale[1];
              pen[2] = 0.0;
              if(matrix)
                transform33f3f(matrix, pen, pen);
              add3f(base, pen, pen);
              if(stroke) {
                CGOEnd(cgo);
              }
              CGOBegin(cgo, GL_LINE_STRIP);
              CGOVertexv(cgo, pen);
              stroke = true;
              break;
            case 1:            /* drawto */
              pen[0] = (*(pc++)) * scale[0];
              pen[1] = (*(pc++)) * scale[1];
              pen[2] = 0.0;
              if(matrix)
                transform33f3f(matrix, pen, pen);
              add3f(base, pen, pen);
              if(stroke) {
                CGOVertexv(cgo, pen);
              }
              break;
            default:
              drawing = false;
              break;
            }
          }
          pen[0] = fr->advance[c] * scale[0];
          pen[1] = 0.0;
          pen[2] = 0.0;
          if(matrix)
            transform33f3f(matrix, pen, pen);
          add3f(pen, pos, pos);
          if(stroke)
            CGOEnd(cgo);
        }
      }
  } else {
    PRINTFB(G, FB_VFont, FB_Errors)
      "VFontWriteToCGO-Error: invalid font identifier (%d)\n", font_id ENDFB(G);
    ok = false;
  }
  return (ok);
}

int VFontIndent(PyMOLGlobals * G, int font_id, char *text, float *pos, float *scale,
                float *matrix, float dir)
{
  register CVFont *I = G->VFont;
  VFontRec *fr = NULL;
  int ok = true;
  float base[3], pen[3];
  float *pc;
  unsigned char c;
  ov_diff offset;

  if((font_id > 0) && (font_id <= I->NFont)) {
    fr = I->Font[font_id];
    if(fr)
      while(1) {
        c = *(text++);
        if(!c)
          break;
        offset = fr->offset[c];
        if(offset >= 0) {
          pc = fr->pen + offset;
          copy3f(pos, base);
          pen[0] = fr->advance[c] * scale[0] * dir;
          pen[1] = 0.0;
          pen[2] = 0.0;
          if(matrix)
            transform33f3f(matrix, pen, pen);
          add3f(pen, pos, pos);
        }
      }
  } else {
    PRINTFB(G, FB_VFont, FB_Errors)
      "VFontIndent-Error: invalid font identifier  (%d)\n", font_id ENDFB(G);
    ok = false;
  }
  return (ok);
}

#ifndef _PYMOL_NOPY
static int VFontRecLoad(PyMOLGlobals * G, VFontRec * I, PyObject * dict)
{                               /* assumes blocked Python interpreter */

  ov_diff used = 0;
  int ok = true;
  PyObject *key, *char_list;
  PyObject *stroke_list = NULL;
#if (PY_MAJOR_VERSION>=2)&&(PY_MINOR_VERSION>=5)
  Py_ssize_t pos = 0;
#else
  int pos = 0;
#endif
  unsigned char code[2];
  float adv;
  ov_diff n_float;
  while(PyDict_Next(dict, &pos, &key, &char_list)) {
    if(!PConvPyStrToStr(key, (char *) code, 2)) {
      PRINTFB(G, FB_VFont, FB_Errors)
        "VFont-Error: Bad character code." ENDFB(G);
      ok = false;
    } else {
      if(ok)
        ok = (char_list != NULL);
      if(ok)
        ok = PyList_Check(char_list);
      if(ok)
        ok = (PyList_Size(char_list) >= 2);
      if(ok)
        ok = PConvPyObjectToFloat(PyList_GetItem(char_list, 0), &adv);
      if(ok) {
        stroke_list = PyList_GetItem(char_list, 1);
        if(ok)
          ok = (stroke_list != NULL);
        if(ok)
          ok = PyList_Check(stroke_list);
        if(ok) {
          n_float = PyList_Size(stroke_list);
          VLACheck(I->pen, float, n_float + used + 1);
          ok = PConvPyListToFloatArrayInPlace(stroke_list, I->pen + used, n_float);
          I->offset[code[0]] = used;
          I->advance[code[0]] = adv;
          I->pen[used + n_float] = -1.0F;       /* sentinel */
          PRINTFD(G, FB_VFont)
            " VFontRecLoad-Debug: Added '%c' adv: %0.3f n_float: %d\n", code[0], adv,
            (int)n_float ENDFD;
          if(ok)
            used += n_float + 1;

        }
      }
    }
  }
  return (ok);
}
#else
#include "vfontdata.h"

static int VFontRecLoad(PyMOLGlobals * G, VFontRec * I)
{
  ov_diff used = 0;
  int ok = true;
  int chidx, n_float, i, off;
  float adv;
  
  for (chidx=0;chidx<VFONT_NUMBER_OF_CHARS; chidx++){
    adv = advs[chidx];
    n_float = n_floats[chidx];
    VLACheck(I->pen, float, n_float + used + 1);
    off = stroke_list_place[chidx];
    for (i=0;i<n_float;i++){
      *(I->pen + used + i) = stroke_lists[off+i];
    }
    I->offset[ch[chidx]] = used;
    I->advance[ch[chidx]] = adv;
    I->pen[used + n_float] = -1.0F;       /* sentinel */
    PRINTFD(G, FB_VFont)
      " VFontRecLoad-Debug: Added '%c' adv: %0.3f n_float: %d\n", ch[chidx], adv,
      (int)n_float ENDFD;
    if(ok)
      used += n_float + 1;
  }
  return (ok);
}

#endif

static void VFontRecFree(PyMOLGlobals * G, VFontRec * I)
{
  VLAFreeP(I->pen);
  OOFreeP(I);
}

int VFontInit(PyMOLGlobals * G)
{
  register CVFont *I = NULL;
  if((I = (G->VFont = Calloc(CVFont, 1)))) {

    register CVFont *I = G->VFont;
    I->Font = VLAlloc(VFontRec *, 10);
    I->NFont = 0;
    return 1;
  } else {
    return 0;
  }
}

void VFontFree(PyMOLGlobals * G)
{
  register CVFont *I = G->VFont;
  int a;
  for(a = 1; a <= I->NFont; a++) {
    VFontRecFree(G, I->Font[a]);
  }
  VLAFreeP(I->Font);
  FreeP(G->VFont);
}

int VFontLoad(PyMOLGlobals * G, float size, int face, int style, int can_load_new)
{
  register CVFont *I = G->VFont;
  VFontRec *fr;
  int a;
  int result = 0;
#ifndef _PYMOL_NOPY
  PyObject *vfont = NULL;
#endif

  PRINTFD(G, FB_VFont)
    " VFontLoad-Debug: Entered %f %d %d\n", size, face, style ENDFD;

  for(a = 1; a <= I->NFont; a++) {
    fr = I->Font[a];
    if((fr->size == size) && (fr->face == face) && (fr->style == style)) {
      result = a;
      break;
    }
  }
  if(!result) {
    if(can_load_new) {
#ifndef _PYMOL_NOPY
      vfont = PGetFontDict(G, size, face, style);
      if(vfont) {
        if(PyDict_Check(vfont)) {
#endif
          VLACheck(I->Font, VFontRec *, I->NFont + 1);
          fr = VFontRecNew(G);
#ifndef _PYMOL_NOPY
          if(!VFontRecLoad(G, fr, vfont))
#else
          if(!VFontRecLoad(G, fr))
#endif
            VFontRecFree(G, fr);
          else {
            I->NFont++;         /* always start at 1 */
            I->Font[I->NFont] = fr;
            result = I->NFont;
            fr->size = size;
            fr->face = face;
            fr->style = style;
          }
#ifndef _PYMOL_NOPY
        }
        Py_DECREF(vfont);
      }
#endif
    }
  }
  PRINTFD(G, FB_VFont)
    " VFontLoad-Debug: Leaving with result %d  (0 = failure)\n", result ENDFD;
  return (result);
}
