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

#include"Ray.h"

typedef struct CViewElem {

  int matrix_flag; /* homogenous transformation, column-major for OpenGL compatibility */
  double matrix[16];

  int pre_flag; /* pre-tranformation translation */
  double pre[3];

  int post_flag; /* post-transformation translation */
  double post[3];

  int clip_flag; /* clip planes (only useful for top-level views...only applied through Scene) */
  float front, back; 
  
  int ortho_flag;
  int ortho;

  int state_flag; /* only applies to object views */
  int state;

  int view_mode; /* 0 = relative/subordinate, 1 = absolute/top-level */

  int specification_level;

  int timing_flag;
  double timing;

} CViewElem;

PyObject *ViewElemAsPyList(CViewElem *view);
int ViewElemFromPyList(PyObject *list, CViewElem *view);

int ViewElemVLAFromPyList(PyObject *list, CViewElem **vla, int nFrame);
PyObject *ViewElemVLAAsPyList(CViewElem *vla, int nFrame);

typedef struct CView {
  PyMOLGlobals *G;
  int NView;
  CViewElem *View;

} CView;

typedef int CViewIterator;

CView *ViewNew(PyMOLGlobals *G);
void ViewFree(CView *I);

CViewIterator ViewGetIterator(CView *I);
int ViewIterate(CView *I,CViewIterator *iter,CRay *ray,int at_least_once);

int ViewElemInterpolate(CViewElem *first,CViewElem *last,float power,float bias);

#endif


