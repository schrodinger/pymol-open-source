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
  I->G=G;
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

    } else if(I->G->HaveGUI && I->G->ValidContext) {
      
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

int ViewElemInterpolate(CViewElem *first,CViewElem *last,float power,float bias)
{
  float first3x3[9];
  float last3x3[9];
  float inverse3x3[9];
  float inter3x3[9];
  float rot_axis[3],trans[3]={0.0F,0.0F,0.0F};
  float angle;
  CViewElem *current;
  int n = (last-first)-1;
  Matrix53f rot,imat;
  int a;
  float tVector[3],tCenter[3],tDir[3];
  float tLen;
  float pLen;
  float v1[3],v2[3],sProj[3];
  float tAngle=0.0F;
  int tLinear = true; /* always do linear for now... */
  float pivot[3];
  const float _1 = 1.0F, _p5 = 0.5F;
  int parabolic = true;
  int timing_flag;
  double timing = 0.0F;

  if(power<0.0F) {
    parabolic = false;
    power = -power;
  }
  
  /* I have no clue whether we're column or row major at this
     point...but it hardly matters!  */

  copy44d33f(first->matrix,first3x3);
  copy44d33f(last->matrix,last3x3);
  transpose33f33f(first3x3,inverse3x3);

  rot[3][0] = (float)first->pre[0];
  rot[3][1] = (float)first->pre[1];
  rot[3][2] = (float)first->pre[2];
  rot[4][0] = (float)last->pre[0];
  rot[4][1] = (float)last->pre[1];
  rot[4][2] = (float)last->pre[2];
/*
  copy3f(first->pre,&rot[3][0]);
  copy3f(last->pre,&rot[4][0]);
*/
  multiply33f33f(inverse3x3,last3x3,&rot[0][0]);

  matrix_to_rotation(rot,rot_axis,&angle);

  if(!tLinear) {
    tVector[0] = (float)(last->pre[0] - first->pre[0]);
    tVector[1] = (float)(last->pre[1] - first->pre[1]);
    tVector[2] = (float)(last->pre[2] - first->pre[2]);
/*
    subtract3f(last->pre,first->pre,tVector);
*/
    normalize23f(tVector,tDir);
    tLen = (float)length3f(tVector); 
    
    /* center of translation */

    tCenter[0] = (float)(last->pre[0] + first->pre[0])/2.0F;
    tCenter[1] = (float)(last->pre[1] + first->pre[1])/2.0F;
    tCenter[2] = (float)(last->pre[2] + first->pre[2])/2.0F;
/*
    average3f(last->pre,first->pre,tCenter); 
*/
    
    pLen = dot_product3f(tDir,rot_axis); /* project translation vector onto rotation axis */
    
    trans[0] -= tDir[0] * pLen; /* subtract component of translation from rotation axis*/
    trans[1] -= tDir[1] * pLen; /* why??? */
    trans[2] -= tDir[2] * pLen;
    normalize3f(trans);
    
    /* rotation axis must be perpendicular to the direction of translation */
    
    cross_product3f(tDir,trans,v1);
    
    normalize3f(v1);
    
    transform33f3f(&rot[0][0],v1,v2);

    project3f(v2,trans,sProj);  /* project vector onto plane _|_ to axis */
    subtract3f(v2,sProj,v2); 
    normalize3f(v2);
  
    tAngle = (float)acos(dot_product3f(v1,v2));
    
    if(fabs(tAngle)>fabs(angle))
      tAngle = (float)(fabs(angle)*(tAngle/fabs(angle)));
    
    /* if translation angle > rotation angle then sets translation angle 
     * to same as rotation angle, with proper sign of course */
    
    if (tAngle>0.0001F) 
      {
        pLen = (float)tan(tAngle/2); 
        if(fabs(pLen)>0.0000001)
          pLen = (tLen/2)/pLen;
        else
          {
            pLen = 1000.0;
            tLinear = true;
          }
      }
    else
      {
        pLen = 1000.0;
        tLinear = true;
      }
  
    pivot[0] = tCenter[0] - pLen * v1[0];
    pivot[1] = tCenter[1] - pLen * v1[1];
    pivot[2] = tCenter[2] - pLen * v1[2];
  }

  timing_flag = first->timing_flag && last->timing_flag;

  current = first+1;
  
  for(a=0;a<n;a++) {
    float fxn = ((float)a+1)/(n+1);
    float fxn_1;

    if(timing_flag) {
      timing = (first->timing * (1.0F - fxn) + (last->timing * fxn));
    }

    /*    printf("Fxn: %8.3f ",fxn);*/
    if(bias!=1.0F) {
      fxn = 1-(float)pow(1-pow(fxn,bias),_1/bias);
    }
    /*    printf("%8.3f bias %8.3f\n",fxn,bias);*/

    if(power!=1.0F) {
      if(fxn<0.5F) {
        if(parabolic) 
          fxn = (float)pow(fxn*2.0F,power)*_p5; /* parabolic */
        else
          fxn = (_1-(float)pow(_1-pow((fxn*2.0F),power),_1/power))*_p5; /* circular */
      } else if(fxn>0.5F) {
        fxn = _1 - fxn;
        if(parabolic) 
          fxn = (float)pow(fxn*2.0F,power)*_p5; /* parabolic */
        else
          fxn = (_1-(float)pow(_1-pow((fxn*2.0F),power),_1/power))*_p5; /* circular */
        fxn = _1 - fxn;
      }
    }

    fxn_1 = 1.0F - fxn;

    *current = *first;
    matrix_interpolate(imat,rot,pivot,rot_axis,angle,tAngle,false,tLinear,fxn);
    current->matrix_flag = true;
    multiply33f33f(first3x3,&imat[0][0],inter3x3);

    copy33f44d(inter3x3,current->matrix);

    if(first->pre_flag && last->pre_flag) {
      copy3f(&imat[4][0],current->pre);
      current->pre_flag=true;
    } else {
      current->pre_flag=false;
    }

    if(first->clip_flag && last->clip_flag) {
      current->front = first->front * fxn_1 + last->front * fxn;
      current->back = first->back * fxn_1 + last->back * fxn;
      current->clip_flag = true;
    } else {
      current->clip_flag = false;
    }

    if(first->post_flag && last->post_flag) {
      mix3d(first->post,last->post,(double)fxn,current->post);
      current->post_flag=true;
    } else {
      current->post_flag=false;
    }
    current->specification_level = 1;

    if(timing_flag) {
      current->timing_flag = true;
      current->timing = timing;
    }

    current++;
  }

  return 1;
}
