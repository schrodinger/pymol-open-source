
/* 
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information. 
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-* 
-* 
-*
Z* -------------------------------------------------------------------
*/
#ifndef _H_View
#define _H_View

#include"os_python.h"
#include"os_gl.h"

#include"Ray.h"
#include"Block.h"

typedef struct CViewElem {

  int matrix_flag;              /* homogenous transformation, column-major for OpenGL compatibility */
  double matrix[16];

  int pre_flag;                 /* pre-tranformation translation */
  double pre[3];

  int post_flag;                /* post-transformation translation */
  double post[3];

  int clip_flag;                /* clip planes (only useful for top-level views...
                                   only applied through Scene) */
  float front, back;

  int ortho_flag;
  float ortho;

  int state_flag;               /* only applies to object views */
  int state;

  int view_mode;                /* 0 = relative/subordinate, 1 = absolute/top-level */

  int specification_level;

  int timing_flag;
  double timing;

  int scene_flag;               /* only applies to main movie view */
  int scene_name;               /* lexicon key */

  int power_flag;
  float power;

  int bias_flag;
  float bias;

} CViewElem;

PyObject *ViewElemAsPyList(PyMOLGlobals * G, CViewElem * view);
int ViewElemFromPyList(PyMOLGlobals * G, PyObject * list, CViewElem * view);

int ViewElemVLAFromPyList(PyMOLGlobals * G, PyObject * list, CViewElem ** vla,
                          int nFrame);
PyObject *ViewElemVLAAsPyList(PyMOLGlobals * G, CViewElem * vla, int nFrame);

void ViewElemArrayPurge(PyMOLGlobals * G, CViewElem * view, int nFrame);
void ViewElemCopy(PyMOLGlobals * G, CViewElem * src, CViewElem * dst);

typedef struct CView {
  PyMOLGlobals *G;
  int NView;
  CViewElem *View;

} CView;

typedef int CViewIterator;

CView *ViewNew(PyMOLGlobals * G);
void ViewFree(CView * I);

CViewIterator ViewGetIterator(CView * I);
int ViewIterate(CView * I, CViewIterator * iter, CRay * ray, int at_least_once);
int ViewElemSmooth(CViewElem * first, CViewElem * last, int window, int loop);

int ViewElemInterpolate(PyMOLGlobals * G, CViewElem * first, CViewElem * last,
                        float power, float bias,
                        int simple, float linearity, int hand, float cut);
void ViewElemDraw(PyMOLGlobals *G, CViewElem * src, BlockRect *rect, int frames, char *title);

#define cViewElemModifyInsert 1
#define cViewElemModifyDelete -1
#define cViewElemModifyMove    2
#define cViewElemModifyCopy    3

int ViewElemModify(PyMOLGlobals *G, CViewElem **handle, int action, int index, int count, int target);
int ViewElemXtoFrame(PyMOLGlobals *G, CViewElem * view_elem, BlockRect *rect, int frames, int x, int nearest);
void ViewElemDrawBox(PyMOLGlobals *G, BlockRect *rect,int first, int last,
                     int frames, float *color4, int fill);

#endif
