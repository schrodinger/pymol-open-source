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

#include"os_predef.h"
#include"os_gl.h"
#include"os_std.h"

#include"main.h"
#include"PyMOLObject.h"
#include"Color.h"
#include"Ortho.h"
#include"Scene.h"
#include"Util.h"
#include"Ray.h"
#include"PConv.h"
#include"MemoryDebug.h"
#include"Movie.h"

#include"Executive.h" 

int ObjectGetNFrames(CObject *I);

void ObjectDescribeElement(struct CObject *I,int index,char *buffer);
CSetting **ObjectGetSettingHandle(struct CObject *I,int state);

void ObjectPurgeSettings(CObject *I)
{
  SettingFreeP(I->Setting);
  I->Setting = NULL;
}

static void TTTToViewElem(float *TTT,CViewElem *elem)
{
  register float *fp = TTT;
  register double *dp;

  /* convert row-major TTT to column-major ViewElem */

  elem->matrix_flag = true;
  dp = elem->matrix;

  dp[ 0] = (double) fp[ 0];
  dp[ 1] = (double) fp[ 4];
  dp[ 2] = (double) fp[ 8];
  dp[ 3] = 0.0;

  dp[ 4] = (double) fp[ 1];
  dp[ 5] = (double) fp[ 5];
  dp[ 6] = (double) fp[ 9];
  dp[ 7] = 0.0;

  dp[ 8] = (double) fp[ 2];
  dp[ 9] = (double) fp[ 6];
  dp[10] = (double) fp[10];
  dp[11] = 0.0;

  dp[12] = 0.0;
  dp[13] = 0.0;
  dp[14] = 0.0;
  dp[15] = 1.0;

  /* copy inverse pre */

  elem->pre_flag = true;
  dp = elem->pre;
  *(dp++) = (double) -TTT[12];
  *(dp++) = (double) -TTT[13];
  *(dp++) = (double) -TTT[14];

  /* copy post */

  elem->post_flag = true;
  dp = elem->post;
  *(dp++) = (double) TTT[ 3];
  *(dp++) = (double) TTT[ 7];
  *(dp++) = (double) TTT[11];
  
}

static void TTTFromViewElem(float *TTT,CViewElem *elem)
{
  register float *fp = TTT;
  register double *dp;

  if(elem->matrix_flag) {
    dp = elem->matrix;

    fp[ 0] = (float) dp[ 0];
    fp[ 1] = (float) dp[ 4];
    fp[ 2] = (float) dp[ 8];
    fp[ 3] = 0.0;
    
    fp[ 4] = (float) dp[ 1];
    fp[ 5] = (float) dp[ 5];
    fp[ 6] = (float) dp[ 9];
    fp[ 7] = 0.0;
    
    fp[ 8] = (float) dp[ 2];
    fp[ 9] = (float) dp[ 6];
    fp[10] = (float) dp[10];
    fp[11] = 0.0;
    
    fp[12] = 0.0;
    fp[13] = 0.0;
    fp[14] = 0.0;
    fp[15] = 1.0;
  }
  
  if(elem->pre_flag) {
    dp = elem->pre;
    fp[12] = (float) (-*(dp++));
    fp[13] = (float) (-*(dp++));
    fp[14] = (float) (-*(dp++));
  }
  
  if(elem->post_flag) {
    dp = elem->post;
    fp[ 3] = (float) *(dp++);
    fp[ 7] = (float) *(dp++);
    fp[11] = (float) *(dp++);
  }   
  fp[15] = 1.0F;
}

int ObjectView(CObject *I,int action,int first,
               int last,float power,float bias,
               int simple, float linear,int wrap,
               int hand,int window,int cycles,int quiet)
{
  register PyMOLGlobals *G = I->G;
  int frame;
  int nFrame = MovieGetLength(I->G);

  if(nFrame<0)
    nFrame = -nFrame;

  if(!I->ViewElem) {
    I->ViewElem = VLACalloc(CViewElem, 0);    
  }

  switch(action) {
  case 0: /* set */
    if(!I->TTTFlag) {
      float mn[3],mx[3],orig[3];
      if(ExecutiveGetExtent(G,I->Name,mn,mx,true,-1,true)) {
        average3f(mn,mx,orig);
        ObjectSetTTTOrigin(I,orig);
      } else {
        initializeTTT44f(I->TTT);
        I->TTTFlag = true;
      }
    }
    if(I->ViewElem && I->TTTFlag) {
      if(first<0)
        first = SceneGetFrame(G);
      if(last<0)
        last = first;
      for(frame=first;frame<=last;frame++) {
        if((frame>=0)&&(frame<nFrame)) {
          VLACheck(I->ViewElem,CViewElem,frame);
          if(!quiet) {
            PRINTFB(G,FB_Object,FB_Details)
              " ObjectView: Setting frame %d.\n",frame+1
              ENDFB(G);
          }
          TTTToViewElem(I->TTT,I->ViewElem+frame);          
          I->ViewElem[frame].specification_level = 2;
        }
      }
    }
    break;
  case 1: /* clear */
    if(I->ViewElem) {
      if(first<0)
        first = SceneGetFrame(G);
      if(last<0)
        last = first;
      for(frame=first;frame<=last;frame++) {
        if((frame>=0)&&(frame<nFrame)) {
          VLACheck(I->ViewElem,CViewElem,frame);
          ViewElemArrayPurge(G,I->ViewElem+frame,1);
          UtilZeroMem((void*)(I->ViewElem+frame),sizeof(CViewElem));
        }
      }
    }
    break;
  case 2: /* interpolate & reinterpolate */
  case 3:
    {
      CViewElem *first_view=NULL,*last_view=NULL;
      int zero_flag = -1;

      if(first<0)
        first = 0;
      if(last<0)
        last = nFrame;
      /* note that we're leaving a blank frame at the end... */

      if(last<0) {
        last = nFrame;
        if(last && !wrap)
          last--;
      }
      if(last>=nFrame) {
        last = nFrame;
        if(last && !wrap)
          last--;
      }

      VLACheck(I->ViewElem,CViewElem,last);
      
      if(wrap && (last == nFrame)) { /* if we're interpolating beyond the
                                 last frame, then wrap by copying
                                 first to last */
        ViewElemCopy(G,I->ViewElem, I->ViewElem+last);
        zero_flag = last;
      }

      VLACheck(I->ViewElem,CViewElem,last);
      if(!quiet) {
        if(action==2) {
          if(last == nFrame) {
            PRINTFB(G,FB_Object,FB_Details)
              " ObjectView: interpolating unspecified frames %d to %d (wrapping).\n",first+1,last
              ENDFB(G);
          } else {
            PRINTFB(G,FB_Object,FB_Details)
              " ObjectView: interpolating unspecified frames %d to %d.\n",first+1,last+1
              ENDFB(G);
          }
        } else {
          if(last == nFrame) {
            PRINTFB(G,FB_Object,FB_Details)
              " ObjectView: reinterpolating all frames %d to %d (wrapping).\n",first+1,last
              ENDFB(G);
          } else {
            PRINTFB(G,FB_Object,FB_Details)
              " ObjectView: reinterpolating all frames %d to %d.\n",first+1,last+1
              ENDFB(G);
          }
        }
      }
      for(frame=first;frame<=last;frame++) {
        if(!first_view) {
          if(I->ViewElem[frame].specification_level==2) { /* specified */
            first_view = I->ViewElem + frame;
          }
        } else {
          CViewElem *view;
          int interpolate_flag = false;
          if(I->ViewElem[frame].specification_level==2) { /* specified */
            last_view = I->ViewElem + frame;
            if(action==2) {/* interpolate */
              for(view=first_view+1;view<last_view;view++) {
                if(!view->specification_level)
                  interpolate_flag = true;
              }
            } else {
              interpolate_flag=true;
            }
            if(interpolate_flag) {
              ViewElemInterpolate(G,first_view,last_view,
                                  power,bias,
                                  simple,linear,hand,0.0F);
            }
            first_view = last_view;
            last_view = NULL;
          }
        }
      }
      if(zero_flag>=0) { /* erase temporary view */
        UtilZeroMem((void*)(I->ViewElem + last), sizeof(CViewElem));
      }
    }
    break;
   case 4: /* smooth */
   {
      if(first<0)
        first = 0;

      if(last<0) {
        last = nFrame;
      }
      if(last>=nFrame) {
        last = nFrame-1;
      }
      if(first<=last) {
        int a;
        VLACheck(I->ViewElem,CViewElem,last);
        for(a=0;a<cycles;a++) {
          ViewElemSmooth(I->ViewElem+first, I->ViewElem + last, window,wrap);
        }
      }
      break;
   }
  case 5: /* reset */
    if(I->ViewElem) {
      int size = VLAGetSize(I->ViewElem);
      VLAFreeP(I->ViewElem);
      I->ViewElem = VLACalloc(CViewElem, size);
    }
    break;

 }
  return 1;
}

void ObjectAdjustStateRebuildRange(CObject *I,int *start, int *stop)
{
  /* on entry, start and stop should hold the valid range for the object */
  int defer_builds_mode = SettingGet_i(I->G,NULL,I->Setting,cSetting_defer_builds_mode);
  int async_builds = SettingGet_b(I->G,NULL,I->Setting,cSetting_async_builds);
  int max_threads =  SettingGet_i(I->G,NULL,I->Setting,cSetting_max_threads);
  if(defer_builds_mode==3) {
    if(SceneObjectIsActive(I->G,I))
      defer_builds_mode=2;
  }
  switch(defer_builds_mode) {
  case 1: /* defer geometry builds until needed */
  case 2: /* defer and destroy continuously for increase memory conservation */
    {
      int min = *start;
      int max = *stop;

      *start = ObjectGetCurrentState(I,false);
      if((!async_builds)||(max_threads<1)) {
        *stop = *start + 1;
      } else {
        int base = (*start / max_threads);
        *start = (base) * max_threads;
        *stop = (base+1) * max_threads;
      }
      if(*start<min) *start = min;
      if(*start>max) *start = max;
      if(*stop <min)  *stop  = min;
      if(*stop >max)  *stop  = max;
    }
    break;
  case 3: /* object not active, so do not rebuild anything */
    *stop = *start;
    break;
  }
}

void ObjectMakeValidName(char *name)
{
  char *p=name,*q;
  if(p) {

    /* currently legal are A to Z, a to z, 0 to 9, -, _, + */
    while(*p) {
      if((*p<43)||(*p>122)||
         ((*p>57)&&(*p<65))||
         ((*p>90)&&(*p<94))||
         (*p==44)||(*p==47)||(*p==60))
        /* must be an ASCII-visible character */
        *p=1; /* placeholder for non-printable */
      p++;
    }
    /* eliminate sequential and terminal nonprintables */
    p=name;
    q=name;
    while(*p) {
      if(q==name)
        while(*p==1)
          p++;
      while((*p==1)&&(p[1]==1))
        p++;
      *q++=*p++;
      if(!p[-1])
        break;
    }
    *q=0;
    while(q>name) {
      if(q[-1]==1) {
        q[-1]=0;
        q--;
      } else
        break;
    }
    /* convert invalides to underscore */
    p=name;
    while(*p) {
      if(*p==1)
        *p='_';
      p++;
    }
  }
}

int ObjectGetCurrentState(CObject *I,int ignore_all_states)
{
  int state=-2;
  int objState;

  if(SettingGetIfDefined_i(I->G,I->Setting,cSetting_state,&objState)) {
    if(objState>0) { /* a specific state */
      state=objState-1;
    } if(objState<0) {
      state=-1; /* all states */
    }
  }
  if(state==-2) { /* default -- use global state */
    state = SettingGetGlobal_i(I->G,cSetting_state)-1;
  }
  if(!(ignore_all_states)&&(state>=0) )
    if(SettingGet_i(I->G,I->Setting,NULL,cSetting_all_states))
      state=-1;
  if(state<-1)
    state=-1;
  return(state);
}

PyObject *ObjectAsPyList(CObject *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;
  result = PyList_New(14);
  PyList_SetItem(result,0,PyInt_FromLong(I->type));
  PyList_SetItem(result,1,PyString_FromString(I->Name));
  PyList_SetItem(result,2,PyInt_FromLong(I->Color));
  PyList_SetItem(result,3,PConvIntArrayToPyList(I->RepVis,cRepCnt));
  PyList_SetItem(result,4,PConvFloatArrayToPyList(I->ExtentMin,3));
  PyList_SetItem(result,5,PConvFloatArrayToPyList(I->ExtentMax,3));
  PyList_SetItem(result,6,PyInt_FromLong(I->ExtentFlag));
  PyList_SetItem(result,7,PyInt_FromLong(I->TTTFlag));
  PyList_SetItem(result,8,SettingAsPyList(I->Setting));
  
  PyList_SetItem(result,9,PyInt_FromLong(I->Enabled));
  PyList_SetItem(result,10,PyInt_FromLong(I->Context));
  PyList_SetItem(result,11,PConvFloatArrayToPyList(I->TTT,16));
  PyList_SetItem(result,11,PConvFloatArrayToPyList(I->TTT,16));
  if(I->ViewElem) {
    int nFrame = VLAGetSize(I->ViewElem);
    PyList_SetItem(result,12,PyInt_FromLong(nFrame));
    PyList_SetItem(result,13,ViewElemVLAAsPyList(I->G,I->ViewElem,nFrame));
  } else {
    PyList_SetItem(result,12,PyInt_FromLong(0));
    PyList_SetItem(result,13,PConvAutoNone(NULL));
  }
  return(PConvAutoNone(result));
#endif
}

int ObjectFromPyList(PyMOLGlobals *G,PyObject *list,CObject *I)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok=true;
  int ll=0;
  I->G = G;
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) ll=PyList_Size(list);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&I->type);
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,1),I->Name,WordLength);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->Color);
  if(ok) I->Color = ColorConvertOldSessionIndex(G,I->Color);
  if(ok) ok = PConvPyListToIntArrayInPlaceAutoZero(PyList_GetItem(list,3),I->RepVis,cRepCnt);
  if(ok) ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list,4),I->ExtentMin,3);
  if(ok) ok = PConvPyListToFloatArrayInPlaceAutoZero(PyList_GetItem(list,5),I->ExtentMax,3);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,6),&I->ExtentFlag);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,7),&I->TTTFlag);
  if(ok) I->Setting=SettingNewFromPyList(G,PyList_GetItem(list,8));
  if(ok&&(ll>9)) ok = PConvPyIntToInt(PyList_GetItem(list,9),&I->Enabled);
  if(ok&&(ll>10)) ok = PConvPyIntToInt(PyList_GetItem(list,10),&I->Context);
  if(ok&&(ll>11)) ok = 
      PConvPyListToFloatArrayInPlaceAutoZero(
            PyList_GetItem(list,11),I->TTT,16);
  if(ok&&(ll>13)) {
    PyObject *tmp;
    int nFrame;
    VLAFreeP(I->ViewElem);
    I->ViewElem = NULL;
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,12),&nFrame);
    if(ok && nFrame) {
      tmp = PyList_GetItem(list,13);
      if(tmp && !(tmp == Py_None))
        ok = ViewElemVLAFromPyList(G,tmp,&I->ViewElem,nFrame);
    }
  }
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  
  return(ok);
#endif
}

int ObjectCopyHeader(CObject *I,CObject *src)
{
  int ok=true;

  I->G = src->G;
  I->type = src->type;
  UtilNCopy(I->Name,src->Name,WordLength);
  I->Color = src->Color;
  {
    int a;
    for(a=0;a<cRepCnt;a++)
      I->RepVis[a] = src->RepVis[a];
  }
  copy3f(src->ExtentMin, I->ExtentMin);
  copy3f(src->ExtentMax, I->ExtentMax);

  I->ExtentFlag = src->ExtentFlag;
  I->TTTFlag = src->TTTFlag;
  I->Setting = NULL; /* to do */
  I->Enabled = src->Enabled;
  I->Context = src->Context;
  {
    int a;
    for(a=0;a<16;a++)
      I->TTT[a] = src->TTT[a];
  }
  I->ViewElem = NULL; /* to do */
  
  return(ok);
}

/*========================================================================*/
void ObjectCombineTTT(CObject *I,float *ttt,int reverse_order)
{
  float cpy[16];
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(cpy);
  } else {
    UtilCopyMem(cpy,I->TTT,sizeof(float)*16);
  }
  if(reverse_order) {
    combineTTT44f44f(cpy,ttt,I->TTT);
  } else {
    combineTTT44f44f(ttt,cpy,I->TTT);
  }
}
/*========================================================================*/
void ObjectTranslateTTT(CObject *I,float *v)
{
#if 1
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(I->TTT);
  }
  I->TTT[3]+=v[0];
  I->TTT[7]+=v[1];
  I->TTT[11]+=v[2];
#else
  float cpy[16];
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(cpy);
  } else {
    UtilCopyMem(cpy,I->TTT,sizeof(float)*16);
  }
 combineTTT44f44f(ttt,cpy,I->TTT);
#endif

}
/*========================================================================*/
void ObjectSetTTT(CObject *I, float *ttt, int state)
{
  if(state<0) {
    if(ttt) {
      UtilCopyMem(I->TTT,ttt,sizeof(float)*16);
      I->TTTFlag=true;
    } else {
      I->TTTFlag=false;
    }
  } else {
    /* to do */
  }
}
/*========================================================================*/
int ObjectGetTTT(CObject *I, float **ttt, int state)
{
  if(state<0) {
    if(I->TTTFlag) {
      *ttt = I->TTT;
      return 1;
    } else {
      *ttt = NULL;
    }
  } else {
  }
  return 0;
}
/*========================================================================*/
void ObjectResetTTT(CObject *I)
{
  I->TTTFlag=false;
}
/*========================================================================*/
int ObjectGetTotalMatrix(CObject *I, int state, int history, double *matrix)
{
  int result = 0;
  if(I->TTTFlag) {
    convertTTTfR44d(I->TTT,matrix);
    result = true;
  }
  
  {
    int use_matrices = SettingGet_b(I->G,I->Setting,NULL,cSetting_matrix_mode);
    if(use_matrices||history) {
      if(I->fGetObjectState) {
        CObjectState *obj_state = I->fGetObjectState(I,state);
        if(obj_state) {
          double *state_matrix = obj_state->Matrix;
          if(state_matrix) {
            if(result) {
              right_multiply44d44d(matrix,state_matrix);
            } else {
              copy44d(state_matrix, matrix);
            }
            result = true;
          }
        }
      }
    }
  }
  return result;
}

/*========================================================================*/
void ObjectPrepareContext(CObject *I,CRay *ray)
{
  if(I->ViewElem) {
    int frame = SceneGetFrame(I->G);
    if(frame>=0) {
      VLACheck(I->ViewElem,CViewElem,frame);


      if(I->ViewElem[frame].specification_level) {
        TTTFromViewElem(I->TTT,I->ViewElem + frame);
        I->TTTFlag=true;
      }
    }
  }
  if(ray) {
    RaySetTTT(ray,I->TTTFlag,I->TTT);
  } else {
    PyMOLGlobals *G = I->G;
    if(G->HaveGUI && G->ValidContext ) {
      if(I->TTTFlag) {
        /* convert the row-major TTT matrix to a column-major OpenGL matrix */
        float gl[16],*ttt;
  
        ttt=I->TTT;
        gl[ 0] = ttt[ 0];
        gl[ 4] = ttt[ 1];
        gl[ 8] = ttt[ 2];
        gl[12] = ttt[ 3];
        gl[ 1] = ttt[ 4];
        gl[ 5] = ttt[ 5];
        gl[ 9] = ttt[ 6];
        gl[13] = ttt[ 7];
        gl[ 2] = ttt[ 8];
        gl[ 6] = ttt[ 9];
        gl[10] = ttt[10];
        gl[14] = ttt[11];
        gl[ 3] = 0.0;
        gl[ 7] = 0.0;
        gl[11] = 0.0;
        gl[15] = 1.0;

        glMultMatrixf(gl);

        /* include the pre-translation */
        glTranslatef(ttt[12],ttt[13],ttt[14]);
      }
    }
  }
}
/*========================================================================*/
void ObjectSetTTTOrigin(CObject *I,float *origin)
{
#if 1
  float homo[16];
  float *dst;
  float post[3];

  if(!I->TTTFlag) {
    I->TTTFlag = true;
    initializeTTT44f(I->TTT);
  }

  /* convert the existing TTT into a homogenous transformation matrix */
  
  convertTTTfR44f(I->TTT, homo);
  
  /* now reset to the passed-in origin */
  
  transform44f3fas33f3f(homo, origin, post);

  homo[ 3] += post[0];
  homo[ 7] += post[1];
  homo[11] += post[2];

  dst = homo+12;
  
  invert3f3f(origin, dst);

  copy44f(homo,I->TTT);
  
#else
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(I->TTT);
  }

  I->TTT[3]+=I->TTT[12]; /* remove existing origin from overall translation */
  I->TTT[7]+=I->TTT[13];
  I->TTT[11]+=I->TTT[14];

  scale3f(origin,-1.0F,I->TTT+12); /* set new origin */

  I->TTT[3]+=origin[0]; /* add new origin into overall translation */
  I->TTT[7]+=origin[1];
  I->TTT[11]+=origin[2];
#endif

}
/*========================================================================*/
CSetting **ObjectGetSettingHandle(struct CObject *I,int state)
{
  return(&I->Setting);
}
/*========================================================================*/
void ObjectDescribeElement(struct CObject *I,int index,char *buffer)
{
  buffer[0]=0;
}
/*========================================================================*/
void ObjectToggleRepVis(CObject *I,int rep)
{
  if((rep>=0)&&(rep<cRepCnt))
    I->RepVis[rep]=!I->RepVis[rep];
}
/*========================================================================*/
void ObjectSetRepVis(CObject *I,int rep,int state)
{
  if((rep>=0)&&(rep<cRepCnt))
    I->RepVis[rep]=state;
}
/*========================================================================*/
void ObjectSetName(CObject *I,char *name)
{
  UtilNCopy(I->Name,name,WordLength);
  if(SettingGetGlobal_b(I->G,cSetting_validate_object_names))
    ObjectMakeValidName(I->Name);
}
/*========================================================================*/
void ObjectUpdate(struct CObject *I);

/*========================================================================*/
void ObjectUpdate(struct CObject *I)
{
  
}
/*========================================================================*/
void ObjectPurge(CObject *I)
{
  if(I) {
    SettingFreeP(I->Setting);
    VLAFreeP(I->ViewElem);
  }
}
/*========================================================================*/
void ObjectFree(CObject *I)
{
  if(I)
    ObjectPurge(I);
}
/*========================================================================*/
int ObjectGetNFrames(CObject *I)
{
  return 1;
}
/*========================================================================*/
void ObjectUseColor(CObject *I)
{
  register PyMOLGlobals *G = I->G;    
  if(G->HaveGUI && G->ValidContext) {
    glColor3fv(ColorGet(I->G,I->Color));
  }
}
/*========================================================================*/
static void ObjectInvalidate(CObject *this,int rep,int level,int state)
{
  
}
/*========================================================================*/
static void ObjectRenderUnitBox(CObject *this,RenderInfo *info)
{
  register PyMOLGlobals *G = this->G;    
  if(G->HaveGUI && G->ValidContext) {

    glBegin(GL_LINE_LOOP);
    glVertex3i(-1,-1,-1);
    glVertex3i(-1,-1, 1);
    glVertex3i(-1, 1, 1);
    glVertex3i(-1, 1,-1);
    
    glVertex3i( 1, 1,-1);
    glVertex3i( 1, 1, 1);
    glVertex3i( 1,-1, 1);
    glVertex3i( 1,-1,-1);
    glEnd();
    
    glBegin(GL_LINES);
    glVertex3i(0,0,0);
    glVertex3i(1,0,0);
    
    glVertex3i(0,0,0);
    glVertex3i(0,3,0);
    
    glVertex3i(0,0,0);
    glVertex3i(0,0,9);

    glEnd();
  }
}
/*========================================================================*/
void ObjectInit(PyMOLGlobals *G,CObject *I)
{
  int a;
  I->G = G;
  I->fFree = ObjectFree;
  I->fRender = ObjectRenderUnitBox;
  I->fUpdate = ObjectUpdate;
  I->fGetNFrame = ObjectGetNFrames;
  I->fDescribeElement = ObjectDescribeElement;
  I->fGetSettingHandle = ObjectGetSettingHandle;
  I->fInvalidate = ObjectInvalidate;
  I->fGetCaption = NULL;
  I->fGetObjectState = NULL;
  I->Name[0]=0;
  I->Color=0; /* white */
  I->ExtentFlag=false;
  I->Setting=NULL;
  I->TTTFlag=false;
  I->Enabled=false;
  zero3f(I->ExtentMin);
  zero3f(I->ExtentMax);
  OrthoRemoveSplash(G);
  for(a=0;a<cRepCnt;a++) I->RepVis[a]=true;
  for(a=0;a<16;a++) I->TTT[a]=0.0F;
  I->RepVis[cRepCell]=false;
  I->RepVis[cRepExtent]=false;
  I->Context=0;
  I->ViewElem = NULL;
}

void ObjectStateInit(PyMOLGlobals *G,CObjectState *I)
{
  I->G = G;
  I->Matrix = NULL;
}

void ObjectStateCopy(CObjectState *dst, CObjectState *src)
{
  *dst = *src;
  if(src->Matrix) {
    dst->Matrix = Alloc(double,16);
    if(dst->Matrix) {
      copy44d(src->Matrix,dst->Matrix);
    }
  }
}

void ObjectStatePurge(CObjectState *I)
{
  FreeP(I->Matrix);
}

void ObjectStateSetMatrix(CObjectState *I, double *matrix)
{
  if(matrix) {
    if(!I->Matrix)
      I->Matrix = Alloc(double,16);
    if(I->Matrix) {
      copy44d(matrix,I->Matrix);
    }
  } else if(I->Matrix) {
    FreeP(I->Matrix);
  }
}

void ObjectStateRightCombineMatrixR44d(CObjectState *I, double *matrix)
{
  if(matrix) {
    if(!I->Matrix) {
      I->Matrix = Alloc(double,16);
      copy44d(matrix,I->Matrix);
    } else {
      right_multiply44d44d(I->Matrix,matrix);
      recondition44d(I->Matrix);
    }
  }
}

void ObjectStateLeftCombineMatrixR44d(CObjectState *I, double *matrix)
{
  if(matrix) {
    if(!I->Matrix) {
      I->Matrix = Alloc(double,16);
      copy44d(matrix,I->Matrix);
    } else {
      left_multiply44d44d(matrix,I->Matrix);
      recondition44d(I->Matrix);
    }
  }
}

void ObjectStateCombineMatrixTTT(CObjectState *I, float *matrix)
{

  if(matrix) {
    if(!I->Matrix) {
      I->Matrix = Alloc(double,16);
      convertTTTfR44d(matrix,I->Matrix);
    } else {
      double tmp[16];
      convertTTTfR44d(matrix,tmp);
      right_multiply44d44d(I->Matrix,tmp);
      recondition44d(I->Matrix);
    }
  }
}

double *ObjectStateGetMatrix(CObjectState *I)
{
  return I->Matrix;
}

void ObjectStateTransformMatrix(CObjectState *I, double *matrix)
{
  if(!I->Matrix) {
    I->Matrix = Alloc(double,16);
    if(I->Matrix) {
      copy44d(matrix,I->Matrix);
    }
  } else {
    right_multiply44d44d(I->Matrix,matrix);
  }
}

int ObjectStatePushAndApplyMatrix(CObjectState *I,RenderInfo *info) 
{
  register PyMOLGlobals *G = I->G;
  float matrix[16];
  register double *i_matrix = I->Matrix;
  int result = false;
  if(i_matrix) {
    if(info->ray) {
      float ttt[16],matrix[16],i_matrixf[16];
      RayPushTTT(info->ray);
      RayGetTTT(info->ray,ttt);
      convertTTTfR44f(ttt,matrix);
      copy44d44f(i_matrix,i_matrixf);
      right_multiply44f44f(matrix,i_matrixf);
      RaySetTTT(info->ray,true,matrix);
      result = true;
    } else if( G->HaveGUI && G->ValidContext) {
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      matrix[ 0] = i_matrix[ 0];
      matrix[ 1] = i_matrix[ 4];
      matrix[ 2] = i_matrix[ 8];
      matrix[ 3] = i_matrix[12];
      matrix[ 4] = i_matrix[ 1];
      matrix[ 5] = i_matrix[ 5];
      matrix[ 6] = i_matrix[ 9];
      matrix[ 7] = i_matrix[13];
      matrix[ 8] = i_matrix[ 2];
      matrix[ 9] = i_matrix[ 6];
      matrix[10] = i_matrix[10];
      matrix[11] = i_matrix[14];
      matrix[12] = i_matrix[ 3];
      matrix[13] = i_matrix[ 7];
      matrix[14] = i_matrix[11];
      matrix[15] = i_matrix[15];
      glMultMatrixf(matrix);
      result = true;
    }
  }
  return result;
}

void ObjectStatePopMatrix(CObjectState *I,RenderInfo *info)
{
  register PyMOLGlobals *G = I->G;
  if(info->ray) {
    RayPopTTT(info->ray);
  } else if(G->HaveGUI && G->ValidContext) {
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
  }
}

void ObjectStateResetMatrix(CObjectState *I)
{
  FreeP(I->Matrix);
}

PyObject *ObjectStateAsPyList(CObjectState *I)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  PyObject *result = NULL;

  if(I) {
    result = PyList_New(1);
    
    if(I->Matrix) {
      PyList_SetItem(result,0,PConvDoubleArrayToPyList(I->Matrix,16));
    } else {
      PyList_SetItem(result,0,Py_None);
    }
  }
  return(PConvAutoNone(result));
#endif  
}

int ObjectStateFromPyList(PyMOLGlobals *G,PyObject *list,CObjectState *I)
{
#ifdef _PYMOL_NOPY
  return 0;
#else
  PyObject *tmp;
  int ok = true;
  int ll = 0;

  ObjectStateInit(G,I);

  if(list&&(list!=Py_None)) { /* allow None */
    if(ok) ok = (list!=NULL);
    if(ok) ok = PyList_Check(list);
    if(ok) ll = PyList_Size(list);
    /* TO SUPPORT BACKWARDS COMPATIBILITY...
       Always check ll when adding new PyList_GetItem's */
    if(ok) {
      tmp = PyList_GetItem(list,0);
      if(tmp!=Py_None) 
        ok = PConvPyListToDoubleArray(tmp,&I->Matrix);
    }
  }
  return(ok);
#endif
}

