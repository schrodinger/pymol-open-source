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

#include"os_gl.h"

#include"Base.h"
#include"OOMac.h"
#include"main.h"
#include"View.h"
#include"Ray.h"
#include"Setting.h"

CView *ViewNew(void)
{
  OOAlloc(CView);
  I->View = NULL;
  return I;
}

void ViewFree(CView *I)
{
  if(I) 
    VLAFreeP(I->View);
}


CViewIterator ViewGetIterator(CView *I)
{
  return 0;
}

int ViewIterate(CView *I,CViewIterator *iter,CRay *ray,int at_least_once)
{
  int result;
  CViewElem *elem = NULL;

  if((!I)||(!I->NView)) { /* trusting short-circuit to avoid segfault */
    if(at_least_once) {
      if(!*iter) { /* do loop at least once if asked to do so */
        *iter = 1;
        result =  true;
      } else 
        result =  false;
    } else {
      result =  false;
    }
  } else {
    if(*iter<I->NView) {
      elem = I->View + (*iter)++;
      result = true;
    } else
      result = false;
  }
  if(elem) { /* are we to apply a transformation? */
    if(ray) {

    } else if(PMGUI) {
      
      if(elem->pre_flag) {
        /* move the camera to the location we are looking at */
        glTranslated(elem->pre[0],elem->pre[1],elem->pre[2]);
      }
      
      if(elem->matrix_flag) {
        /* rotate about the origin (the the center of rotation) */
        glMultMatrixd(elem->matrix);
      }
      
      if(elem->post_flag) {
        /* move the origin to the center of rotation */
        glTranslated(elem->post[0],elem->post[1],elem->post[2]);
      }
      
    }
  }
  return result;
}

