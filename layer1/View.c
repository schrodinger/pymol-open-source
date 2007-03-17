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

void ViewElemCopy(PyMOLGlobals *G,CViewElem *src,CViewElem *dst)
{
  *dst = *src;
  if(dst->scene_flag && dst->scene_name) {
    OVLexicon_IncRef(G->Lexicon,dst->scene_name);
    dst->scene_name=0;
    dst->scene_flag=false;
  }
}

void ViewElemArrayPurge(PyMOLGlobals *G,CViewElem *view,int nFrame)
{
  int a;
  for(a=0;a<nFrame;a++) {
    if(view->scene_flag && view->scene_name ) {
      OVLexicon_DecRef(G->Lexicon,view->scene_name);
      view->scene_name=0;
      view->scene_flag=false;
    }
    view++;
  }
}

PyObject *ViewElemAsPyList(PyMOLGlobals *G, CViewElem *view)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyObject *result = NULL;

  result=PyList_New(15);

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

    PyList_SetItem(result,13,PyInt_FromLong(view->scene_flag));  
    
    if(view->scene_flag && view->scene_name) {
      char null_st[1] = "";
      char *st = null_st;
      
      st = OVLexicon_FetchCString(G->Lexicon,view->scene_name);
      PyList_SetItem(result, 14,PyString_FromString(st));
    } else {
      PyList_SetItem(result, 14,PyInt_FromLong(0));
    }
  }

  return PConvAutoNone(result);
#endif
}

int ViewElemFromPyList(PyMOLGlobals *G, PyObject *list, CViewElem *view)
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

  if(ok&(ll>14)) {
    if(ok) ok= PConvPyIntToInt(PyList_GetItem(list,13),&view->scene_flag);
    if(ok && view->scene_flag) {
      char *ptr = NULL;
      view->scene_flag = false;
      if(PConvPyStrToStrPtr(PyList_GetItem(list, 14),&ptr)) {
        OVreturn_word result = OVLexicon_GetFromCString(G->Lexicon,ptr);
        if(OVreturn_IS_OK(result)) {
          view->scene_name = result.word;
          view->scene_flag = true;
        }
      }
    }
  }
    
  return ok;
#endif
}

int ViewElemVLAFromPyList(PyMOLGlobals *G, PyObject *list, CViewElem **vla_ptr, int nFrame)
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
        ok=ViewElemFromPyList(G,PyList_GetItem(list,a),vla+a);
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

PyObject *ViewElemVLAAsPyList(PyMOLGlobals *G, CViewElem *vla,int nFrame)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;
  int a;
  result = PyList_New(nFrame);
  for(a=0;a<nFrame;a++) {
    PyList_SetItem(result,a,ViewElemAsPyList(G,vla+a));
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


static void matrix_interpolate(Matrix53f imat,Matrix53f mat,
                               float *pivot_point,
                               float *bisect_dir,
                               float *rot_axis,
                               float rotate_angle,
                               float *trans_axis,
                               float translate_angle, 
                               float fxn, float linearity)
{
  int a;
  float pos[3],adj[3],opp[3],oppdir[3];
  float p0[3],p1[3],center[3];
  float hyplen,adjlen,opplen;
  float tAlpha;

  rotation_to_matrix(imat,rot_axis,fxn*rotate_angle);
  
  
  /*           ______--------______
   *        /____________          \
   *     /   \   opp     |adj         \
   *   |      \          |        trans   | 
   * (CM)---------------------------------->(CM)
   *    \       \        |               /
   *       \      \hyp   |-bisect_dir  /
   *          \p0   \    |        p1/
   *             \    \  |        /
   *                \  \ |     /     
   *                   \\v  /
   *            <--------O pivot
   *                      F-raxis
   */

  subtract3f(&mat[3][0],pivot_point,p0);
  subtract3f(&mat[4][0],pivot_point,p1);
  
  hyplen = (float)length3f(p0);

  average3f(&mat[3][0],&mat[4][0],center);

  cross_product3f(bisect_dir,trans_axis,oppdir);
  normalize3f(oppdir);

  tAlpha = (float)(fabs(0.5-fxn)*translate_angle);
  opplen = (float)fabs(hyplen * sin(tAlpha));
  adjlen = (float)fabs(hyplen * cos(tAlpha));
  
  scale3f(oppdir,opplen,opp);
  scale3f(bisect_dir,-adjlen,adj);
  
  add3f(pivot_point,adj,pos);
  
  if(fxn<=0.5) {
    add3f(pos,opp,pos);
  } else {
    subtract3f(pos,opp,pos);
  }
  
  /* straight linear for now... */

  for(a=0;a<3;a++) {
    imat[4][a] = (((1.0-fxn)*mat[3][a] + fxn*mat[4][a]) * linearity) +
      (1.0-linearity) * pos[a];
  }
}

int ViewElemSmooth(CViewElem *first,CViewElem *last,int window,int loop)
{
  int n = (last-first)+1;
  int delta;
  loop = 0; /* FOR NOW */
  if(window>n)
    window=n;
  delta = (window-1)/2;
  if(n&&delta) {
    CViewElem *cpy = Alloc(CViewElem,(n+2*delta));
    CViewElem *src,*dst;
    int a,b,c,cnt;
    memcpy(cpy+delta,first,sizeof(CViewElem)*n);
    if(loop) {
      for(a=0;a<delta;a++) {
        memcpy(cpy+a,last-delta+a,sizeof(CViewElem));
        memcpy(cpy+(delta+n)+a,first+delta-a,sizeof(CViewElem));
      }
    } else {
      for(a=0;a<delta;a++) {
        memcpy(cpy+a,first,sizeof(CViewElem));
        memcpy(cpy+(delta+n)+a,last,sizeof(CViewElem));
      }
    }
    for(a=0;a<n;a++) {
      int above,below;
      dst = first+a;
      
      above = delta;
      below = delta;
      if(above>a) above = a;
      if(below>((n-1)-a)) below = (n-1)-a;

      if(dst->specification_level) { /* has to be specified */

        if(dst->matrix_flag) {
          cnt=1;
          for(b=-below;b<=above;b++) {
            if(b) {
              src = cpy+delta+a+b;
              if(src->matrix_flag) {
                cnt++;
                for(c=0;c<16;c++) {
                  dst->matrix[c]+= src->matrix[c];
                }
              }
            }
          }
          for(c=0;c<16;c++) {
            dst->matrix[c] /= cnt;
          }
          reorient44d(dst->matrix); /* convert those averages into a valid matrix */
        }
        
        if(dst->pre_flag) {
          cnt=1;
          for(b=-below;b<=above;b++) {
            if(b) {
              src = cpy+delta+a+b;
              if(src->pre_flag) {
                cnt++;
                for(c=0;c<3;c++) {
                  dst->pre[c]+= src->pre[c];
                }
              }
            }
          }
          for(c=0;c<3;c++) {
            dst->pre[c] /= cnt;
          }
        }

        if(dst->post_flag) {
          cnt=1;
          for(b=-below;b<=above;b++) {
            if(b) {
              src = cpy+delta+a+b;
              if(src->post_flag) {
                cnt++;
                for(c=0;c<3;c++) {
                  dst->post[c]+= src->post[c];
                }
              }
            }
          }
          for(c=0;c<3;c++) {
            dst->post[c] /= cnt;
          }
        }

        if(dst->clip_flag) {
          cnt=1;
          for(b=-below;b<=above;b++) {
            if(b) {
              src = cpy+delta+a+b;
              if(src->clip_flag) {
                cnt++;
                dst->front += src->front;
                dst->back += src->back;
              }
            }
          }
          dst->front /= cnt;
          dst->back /= cnt;
        }
        
      }
    }
  }
  return 1;
}

int ViewElemInterpolate(PyMOLGlobals *G,CViewElem *first,CViewElem *last,
                        float power,float bias,
                        int simple, float linearity, int hand, float cut)
{
  float first3x3[9];
  float last3x3[9];
  float inverse3x3[9];
  float inter3x3[9];
  float rot_axis[3],trans_axis[3]={0.0F,0.0F,0.0F};
  float angle;
  CViewElem *current;
  int n = (last-first)-1;
  Matrix53f rot,imat;
  int a;
  float tVector[3],tCenter[3],tDir[3];
  float tLen=0.0F;
  float bisect[3],v2[3];
  float translate_angle=0.0F;
  float pivot[3] = {0.0F,0.0F,0.0F};
  const float _1 = 1.0F, _p5 = 0.5F;
  int parabolic = true;
  int timing_flag;
  double timing = 0.0F;
  float pre[3];
  float firstC44f[16], firstRTTT[16], firstR44f[16];
  float lastC44f[16], lastRTTT[16], lastR44f[16];
  int linear = false;
  int debug = false;

  if(power<0.0F) {
    parabolic = false;
    power = -power;
  }

  /* WARNING: this routine is operating on column-major matrices!!! */

  copy44d33f(first->matrix,first3x3);
  copy44d33f(last->matrix,last3x3);

  transpose33f33f(first3x3,inverse3x3);
  
  multiply33f33f(inverse3x3,last3x3,&rot[0][0]);  /* [rot] = [first]^-1 [last] */
  matrix_to_rotation(rot,rot_axis,&angle);

  if(debug) dump3f(rot_axis,"rot_axis");
  
  if(hand) {
    if((cPI-fabs(angle))<0.01F) { /* this a complete 180 degree motion */
      if(((rot_axis[0]*0.7F + rot_axis[1]*0.8F + rot_axis[2]*0.9F) * hand * angle)>0.0F) {
        invert3f(rot_axis);
        if(angle>0) {
          angle = (2*cPI) - angle;
        } else {
          angle = -(2*cPI) - angle;
        }
      }
    }
  }

  if(!simple) {
    /* switch back into row major to promote developer sanity */

    copy33f44f(first3x3,firstC44f);
    copy33f44f(last3x3,lastC44f);
    
    transpose44f44f(firstC44f,firstRTTT);
    transpose44f44f(lastC44f,lastRTTT);
    
    /* form TTTs */

    firstRTTT[12] = -first->pre[0];
    firstRTTT[13] = -first->pre[1];
    firstRTTT[14] = -first->pre[2];

    firstRTTT[ 3] = first->post[0];
    firstRTTT[ 7] = first->post[1];
    firstRTTT[11] = first->post[2];

    lastRTTT[12] = -last->pre[0];
    lastRTTT[13] = -last->pre[1];
    lastRTTT[14] = -last->pre[2];

    lastRTTT[ 3] = last->post[0];
    lastRTTT[ 7] = last->post[1];
    lastRTTT[11] = last->post[2];

    if(debug) dump44f(firstRTTT,"firstRTTT");
    if(debug) dump44f(lastRTTT,"lastRTTT");

    /* convert to homogenous */

    convertTTTfR44f(firstRTTT, firstR44f);
    convertTTTfR44f(lastRTTT, lastR44f);
    
    /* reset both matrices to a common origin */

    {
      float first_pre[3],last_pre[3];
      float post[4], *dst;
      
      copy3d3f(first->pre,first_pre);
      copy3d3f(last->pre,last_pre);
      average3f(first_pre, last_pre, pre);
    
      transform44f3fas33f3f(firstR44f, pre, post);
      copy44f(firstR44f,firstRTTT);
      firstRTTT[ 3] += post[0];
      firstRTTT[ 7] += post[1];
      firstRTTT[11] += post[2];
      dst = firstRTTT+12;
      invert3f3f(pre,dst);
      
      transform44f3fas33f3f(lastR44f, pre, post);
      copy44f(lastR44f,lastRTTT);
      lastRTTT[ 3] += post[0];
      lastRTTT[ 7] += post[1];
      lastRTTT[11] += post[2];
      dst = lastRTTT+12;
      invert3f3f(pre,dst);
    }
    
    if(debug) dump44f(firstRTTT,"firstRTTT");
    if(debug) dump44f(lastRTTT,"lastRTTT");

    /*    convertTTTfR44f(firstRTTT, firstR44f);
          convertTTTfR44f(lastRTTT, lastR44f);*/

    /* now populate the translation fields */

    rot[3][0] = firstRTTT[ 3];
    rot[3][1] = firstRTTT[ 7];
    rot[3][2] = firstRTTT[11];

    rot[4][0] = lastRTTT[ 3];
    rot[4][1] = lastRTTT[ 7];
    rot[4][2] = lastRTTT[11];

    /* now set up the interpolation */

    subtract3f(&rot[4][0],&rot[3][0],tVector);
    tLen = (float)length3f(tVector); 
    average3f(&rot[4][0],&rot[3][0],tCenter);

    if(tLen<0.0001F) {
      if(debug) printf("translation too short %8.3f\n",tLen);
      simple=true;
    }
  }

  if(!simple) {
    
    normalize23f(tVector,tDir);
    if(debug) dump3f(tDir,"tDir");
    cross_product3f(tDir,rot_axis,bisect);
    /* bisect is a vector in the translation arc */
    if(length3f(bisect)<0.0001F) {
      if(debug) printf("rotation coincident with translation\n");
      linear = true;
    }
  }

  if(!(simple||linear)) {
    normalize3f(bisect);

    /* this section needs work... */

    cross_product3f(bisect, tDir, trans_axis);
    normalize3f(trans_axis);

    transform33Tf3f(&rot[0][0],bisect,v2); /* column major */

    remove_component3f(v2,trans_axis,v2);
    normalize3f(v2);  /* project vector onto plane _|_ to axis */

    if(debug) {
      dump3f(rot_axis,"rot_axis");
      dump3f(tDir,"tDir");
      dump3f(bisect,"bisect");
      dump3f(trans_axis,"trans_axis");
      dump3f(v2,"v2");
    }

    {
      double dot = dot_product3f(bisect,v2);
      if(dot<-1.0F) dot=-1.0F;
      if(dot>1.0F) dot=1.0F;
      translate_angle = (float)acos(dot);
      
      /* if translation angle > rotation angle then sets translation angle 
       * to same as rotation angle, with proper sign of course */
      
      if(fabs(translate_angle)>fabs(angle))
        translate_angle = (float)(fabs(angle)*(translate_angle/fabs(angle)));

      if(fabs(translate_angle)<0.0001F) {
        linear = true;
        if(debug) printf("no significant rotation\n");
      }
      
      if((translate_angle * angle)<0.0F) {
        /* if motions are in opposing directions, then flip translation axis and location */
        invert3f(bisect);
        invert3f(trans_axis);
      }

    }
  }

  if(!(simple||linear)) {
    float pLen = (float)tan(translate_angle/2); 
    if(fabs(pLen)>0.0000001)
      pLen = (tLen/2)/pLen;
    else {
      if(debug) printf("pLen too short %8.3f\n",pLen);
      simple = true;
    }
    
    if(!simple) {
      pivot[0] = tCenter[0] + pLen * bisect[0];
      pivot[1] = tCenter[1] + pLen * bisect[1];
      pivot[2] = tCenter[2] + pLen * bisect[2];
    }

    if(debug && !simple) {
      dump3f(tCenter,"center");
      dump3f(pivot,"pivot");
      printf("pLen %8.3f angle %8.3f translate_angle %8.3f\n",
             pLen,angle,translate_angle);
    }
  }

  /* now interpolate */

  timing_flag = first->timing_flag && last->timing_flag;

  current = first+1;
  
  if(debug) dump44f(firstR44f,"first");
  
  for(a=0;a<n;a++) {
    float fxn = ((float)a+1)/(n+1);
    float fxn_1;

    if(timing_flag) {
      timing = (first->timing * (1.0F - fxn) + (last->timing * fxn));
    }

    if(bias!=1.0F) {
      fxn = 1-(float)pow(1-pow(fxn,bias),_1/bias);
    }

    if((power!=1.0F)||(!parabolic)) {
      if(fxn<0.5F) {
        if(!parabolic) 
          fxn = (_1-cos(cPI*fxn))*_p5; /* circular */
        fxn = (float)pow(fxn*2.0F,power)*_p5; /* parabolic*/
      } else if(fxn>0.5F) {
        fxn = _1 - fxn;
        if(!parabolic) 
          fxn = (_1-cos(cPI*fxn))*_p5;
        fxn = (float)pow(fxn*2.0F,power)*_p5; /* parabolic */
        fxn = _1 - fxn;
      }
    }

    fxn_1 = 1.0F - fxn;

    ViewElemCopy(G,first,current);

    if(simple) {
      rotation_matrix3f(fxn*angle,rot_axis[0],rot_axis[1],rot_axis[2],&imat[0][0]);

/* [cur] = [first] [partial-rot], so....
   at start: [cur] = [first] [identity] = [first]
   at end: [cur] = [first] [first]^-1 [last] = [last] 
*/
      current->matrix_flag = true;
      multiply33f33f(first3x3,&imat[0][0],inter3x3); 
      copy33f44d(inter3x3,current->matrix);
      
      if(first->pre_flag && last->pre_flag) {
        mix3d(first->pre,last->pre,(double)fxn,current->pre);
        current->pre_flag=true;
      } else {
        current->pre_flag=false;
      }
      if(first->post_flag && last->post_flag) {
        mix3d(first->post,last->post,(double)fxn,current->post);
        current->post_flag=true;
      } else {
        current->post_flag=false;
      }

    } else if(linear) {
      int b;
      rotation_matrix3f(fxn*angle,rot_axis[0],rot_axis[1],rot_axis[2],&imat[0][0]);
      current->matrix_flag = true;
      multiply33f33f(first3x3,&imat[0][0],inter3x3); 
      copy33f44d(inter3x3,current->matrix);
      
      current->pre_flag = true;
      copy3f3d(pre,current->pre);
      
      current->post_flag = true;
      for(b=0;b<3;b++) {
        imat[4][b] = ((1.0-fxn)*rot[3][b] + fxn*rot[4][b]);
      }
      copy3f3d(&imat[4][0],current->post);
       
    } else {

      matrix_interpolate(imat,rot,
                         pivot,bisect,
                         rot_axis,angle,
                         trans_axis,translate_angle,
                         fxn,linearity);

      current->matrix_flag = true;
      multiply33f33f(first3x3,&imat[0][0],inter3x3); 
      copy33f44d(inter3x3,current->matrix);

      current->pre_flag = true;
      copy3f3d(pre,current->pre);

      current->post_flag = true;
      copy3f3d(&imat[4][0],current->post);

    }

    if(debug) {
      if((a==0)||(a==n-1)) {
        float curC44f[16],curRTTT[16],curR44f[16];
        
        copy33f44f(inter3x3,curC44f);
        
        transpose44f44f(curC44f,curRTTT);
        
        /* form TTTs */
        
        curRTTT[12] = -current->pre[0];
        curRTTT[13] = -current->pre[1];
        curRTTT[14] = -current->pre[2];
        
        curRTTT[ 3] = current->post[0];
        curRTTT[ 7] = current->post[1];
        curRTTT[11] = current->post[2];
        
        convertTTTfR44f(curRTTT, curR44f);
        dump44f(curR44f,"cur");
      }
    }
    
    if(first->clip_flag && last->clip_flag) {
      current->front = first->front * fxn_1 + last->front * fxn;
      current->back = first->back * fxn_1 + last->back * fxn;
      current->clip_flag = true;
    } else {
      current->clip_flag = false;
    }
    
    if(first->ortho_flag && last->ortho_flag) {
      current->ortho = ((first->ortho * fxn_1 + last->ortho * fxn)>0.5F);
    }
    current->specification_level = 1;
    
    if(timing_flag) {
      current->timing_flag = true;
      current->timing = timing;
    }

    if(first->scene_flag && last->scene_flag) {
      if(current->scene_name) {
        OVLexicon_DecRef(G->Lexicon,current->scene_name);
      }
      current->scene_flag = true;
      if(fxn>=cut) {
        current->scene_name = last->scene_name;
      } else {
        current->scene_name = first->scene_name;
      }
      OVLexicon_IncRef(G->Lexicon,current->scene_name);
    }
    current++;
  }
  if(debug) dump44f(lastR44f,"last");

  return 1;
}
