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

#include"Ray.h"

typedef struct CViewElem {

  int matrix_flag; /* homogenous transformation (rotation, translation, etc.) */
  double matrix[16];

  int pre_flag; /* pre-translation (position) */
  double pre[3];

  int post_flag; /* post-translation (origin) */
  double post[3];

  int clip_flag; /* clip planes (only useful for top-level views...only applied through Scene) */
  float front, back; 
  
  int ortho_flag;
  int ortho;

  int state_flag; /* only applies to object views */
  int state;

  int view_mode; /* 0 = relative/subordinate, 1 = absolute/top-level */
  int specified;
} CViewElem;

typedef struct CView {
  int NView;
  CViewElem *View;

} CView;

typedef int CViewIterator;

CView *ViewNew(void);
void ViewFree(CView *I);

CViewIterator ViewGetIterator(CView *I);
int ViewIterate(CView *I,CViewIterator *iter,CRay *ray,int at_least_once);

#endif


