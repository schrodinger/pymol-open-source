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

int ObjectGetNFrames(CObject *I);

void ObjectDescribeElement(struct CObject *I,int index,char *buffer);
CSetting **ObjectGetSettingHandle(struct CObject *I,int state);

void ObjectMakeValidName(char *name)
{
  char *p=name,*q;
  if(p) {
    while(*p) {
      if((*p<45)||(*p>122)||
         ((*p>57)&&(*p<65))||
         ((*p>90)&&(*p<94))||
         (*p==47)||(*p==60))
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
  result = PyList_New(12);
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
  if(ok) ok = PConvPyStrToStr(PyList_GetItem(list,1),I->Name,ObjNameMax);
  if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,2),&I->Color);
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

  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  
  return(ok);
#endif
}
/*========================================================================*/
void ObjectCombineTTT(CObject *I,float *ttt)
{
  float cpy[16];
  if(!I->TTTFlag) {
    I->TTTFlag=true;
    initializeTTT44f(cpy);
  } else {
    UtilCopyMem(cpy,I->TTT,sizeof(float)*16);
  }
  combineTTT44f44f(ttt,cpy,I->TTT);
}
/*========================================================================*/
void ObjectResetTTT(CObject *I)
{
  I->TTTFlag=false;
  SceneDirty(I->G);
}
/*========================================================================*/
void ObjectPrepareContext(CObject *I,CRay *ray)
{
  float gl[16],*ttt;
  if(ray) {
    RaySetTTT(ray,I->TTTFlag,I->TTT);
  } else {
    if(I->G->HaveGUI) {
      if(I->TTTFlag) {
        /* form standard 4x4 GL matrix with TTT rotation and 2nd translation */
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
        /*        dump44f(gl,"ttt");
                  dump44f(gl,"gl");*/
        glMultMatrixf(gl);
        /* now add in the first translation */
        glTranslatef(ttt[12],ttt[13],ttt[14]);
      }
    }
  }
}

/*========================================================================*/
void ObjectSetTTTOrigin(CObject *I,float *origin)
{
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

  SceneDirty(I->G);

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
  UtilNCopy(I->Name,name,ObjNameMax);
  ObjectMakeValidName(I->Name);
}
/*========================================================================*/
void ObjectRenderUnitBox(struct CObject *this,int frame,CRay *ray,Pickable **pick,int pass);
void ObjectUpdate(struct CObject *I);

/*========================================================================*/
void ObjectUpdate(struct CObject *I)
{
  
}
/*========================================================================*/
void ObjectPurge(CObject *I)
{
  if(I) 
    SettingFreeP(I->Setting);
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
  if(I->G->HaveGUI) glColor3fv(ColorGet(I->G,I->Color));
}
/*========================================================================*/
static void ObjectInvalidate(CObject *this,int rep,int level,int state)
{
  
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
}
/*========================================================================*/
void ObjectRenderUnitBox(CObject *this,int frame,
                         CRay *ray,Pickable **pick,int pass)
{
  if(this->G->HaveGUI) {
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




