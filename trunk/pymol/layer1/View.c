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
#include"PConv.h"


PyObject *ViewElemAsPyList(CViewElem *view)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;

  result=PyList_New(13);

  if(result) {
    PyList_SetItem(result,0,PyInt_FromLong(view->matrix_flag));
    if(view->matrix_flag) {
      PyList_SetItem(result,1,PConvDoubleArrayToPyList(view->matrix,16));
    } else {
      PyList_SetItem(result,1,PConvAutoNone(NULL));
    }
    
    PyList_SetItem(result,2,PyInt_FromLong(view->pre_flag));
    if(view->pre_flag) {
      PyList_SetItem(result,3,PConvDoubleArrayToPyList(view->pre,3));
    } else {
      PyList_SetItem(result,3,PConvAutoNone(NULL));
    }
    
    PyList_SetItem(result,4,PyInt_FromLong(view->post_flag));
    if(view->post_flag) {
      PyList_SetItem(result,5,PConvDoubleArrayToPyList(view->post,3));
    } else {
      PyList_SetItem(result,5,PConvAutoNone(NULL));
    }
    
    PyList_SetItem(result,6,PyInt_FromLong(view->clip_flag));
    if(view->post_flag) {
      PyList_SetItem(result,7,PyFloat_FromDouble((double)view->front));
      PyList_SetItem(result,8,PyFloat_FromDouble((double)view->back));
    } else {
      PyList_SetItem(result,7,PConvAutoNone(NULL));
      PyList_SetItem(result,8,PConvAutoNone(NULL));
    }
    
    PyList_SetItem(result,9,PyInt_FromLong(view->ortho_flag));
    if(view->ortho_flag) {
      PyList_SetItem(result,10,PyInt_FromLong(view->ortho));
    } else {
      PyList_SetItem(result,10,PConvAutoNone(NULL));
    }
    
    PyList_SetItem(result,11,PyInt_FromLong(view->view_mode));
    
    PyList_SetItem(result,12,PyInt_FromLong(view->specification_level));
  }

  return PConvAutoNone(result);
#endif
}

int ViewElemFromPyList(PyObject *list, CViewElem *view)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  int ll=0;

  if(ok) ok= (list!=NULL);
  if(ok) ok= PyList_Check(list);
  if(ok) ok= ((ll=PyList_Size(list))>11);

  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,0),&view->matrix_flag);
  if(ok&&view->matrix_flag) ok= PConvPyListToDoubleArrayInPlace(PyList_GetItem(list,1),view->matrix,16);
  
  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,2),&view->pre_flag);
  if(ok&&view->pre_flag) ok= PConvPyListToDoubleArrayInPlace(PyList_GetItem(list,3),view->pre,3);

  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,4),&view->post_flag);
  if(ok&&view->post_flag) ok= PConvPyListToDoubleArrayInPlace(PyList_GetItem(list,5),view->post,3);
  
  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,6),&view->clip_flag);
  if(view->post_flag) {
    if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,7),&view->front);
    if(ok) ok = PConvPyFloatToFloat(PyList_GetItem(list,8),&view->back);
  }

  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,9),&view->ortho_flag);
  if(ok&&view->ortho_flag) ok= PConvPyIntToInt(PyList_GetItem(list,10),&view->ortho_flag);

  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,11),&view->view_mode);
  if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,12),&view->specification_level);

  return ok;
#endif
}

int ViewElemVLAFromPyList(PyObject *list, CViewElem **vla_ptr, int nFrame)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;

  CViewElem *vla = NULL;


  if(ok) ok= (list!=NULL);
  if(ok) ok= PyList_Check(list);
  if(ok) ok= (PyList_Size(list)==nFrame);
  if(ok) ok= ((vla=VLACalloc(CViewElem,nFrame))!=NULL);
  if(ok) {
    int a;
    for(a=0;a<nFrame;a++) {
      if(ok) 
        ok=ViewElemFromPyList(PyList_GetItem(list,a),vla+a);
      else
        break;
    }
  }
  if(!ok) {
    VLAFreeP(vla);
  } else
    *vla_ptr = vla;
  return ok;
#endif
}

PyObject *ViewElemVLAAsPyList(CViewElem *vla,int nFrame)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;
  int a;
  result = PyList_New(nFrame);
  for(a=0;a<nFrame;a++) {
    PyList_SetItem(result,a,ViewElemAsPyList(vla+a));
  }
  return(PConvAutoNone(result));
#endif
}

CView *ViewNew(PyMOLGlobals *G)
{
  OOAlloc(G,CView);
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

