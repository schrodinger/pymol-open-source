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
#include"os_std.h"
#include"os_gl.h"

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Executive.h"
#include"Ortho.h"
#include"Movie.h"
#include"Scene.h"
#include"MyPNG.h"
#include"P.h"
#include"Setting.h"
#include"main.h"
#include"PConv.h"
#include"Util.h"
#include"Parse.h"

struct _CMovie {
  ImageType *Image;
  int *Sequence;
  MovieCmdType *Cmd;
  int NImage,NFrame;
  unsigned Width,Height;
  int MatrixFlag;
  SceneViewType Matrix;
  int Playing;
  int Locked;
  int CacheSave;
  CViewElem *ViewElem;
};

void MovieCopyPrepare(PyMOLGlobals *G,int *width,int *height,int *length)
{
  /* assumed locked api, blocked threads, and master thread on entry */

  /* writes the current movie sequence to a set of PNG files 
  * this routine can take a LOT of time -- 
  * TODO: develop an interrupt mechanism */

  int start,stop;
  register CMovie *I=G->Movie;
  int nFrame;

  I->CacheSave = (int)SettingGet(G,cSetting_cache_frames); 
  if(!I->CacheSave)
    MovieClearImages(G);
  SettingSet(G,cSetting_cache_frames,1.0);
  nFrame = I->NFrame;
  if(!nFrame) {
	 nFrame=SceneGetNFrame(G);
  }
  start=0;
  stop=nFrame;
  if((start!=0)||(stop!=(nFrame+1)))
    SceneSetFrame(G,0,0);
  MoviePlay(G,cMoviePlay);
  VLACheck(I->Image,ImageType,nFrame);
  SceneGetWidthHeight(G,width,height);
  *length = nFrame;
}

int MovieCopyFrame(PyMOLGlobals *G,int frame,int width,int height,int rowbytes,void *ptr)
{
  register CMovie *I=G->Movie;
  int result=false;
  int nFrame;
  
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame=SceneGetNFrame(G);
  }
  
  if((width==I->Width)&&(height==I->Height)&&(frame<nFrame)&&(ptr)) {
    int a = frame;
    int i;
    SceneSetFrame(G,0,a);
    MovieDoFrameCommand(G,a);
    PFlush();
    i=MovieFrameToImage(G,a);
    VLACheck(I->Image,ImageType,i);
    if(!I->Image[i]) {
      SceneMakeMovieImage(G);
    }
    if(!I->Image[i]) {
      PRINTFB(G,FB_Movie,FB_Errors) 
        "MoviePNG-Error: Missing rendered image.\n"
        ENDFB(G);
    } else {
      {
        unsigned char *srcImage = (unsigned char*)I->Image[i];
        int i,j;
        for (i=0; i< height; i++)
          {
            unsigned char *dst = ((unsigned char*)ptr) + i * rowbytes;
            unsigned char *src = srcImage + ((height-1)-i) * width*4;
            for (j = 0; j < width; j++)
              {
                *dst++ = src[3];
                *dst++ = src[0];
                *dst++ = src[1];
                *dst++ = src[2];
                src+=4;
              }
          }
        result = true;
      }
      ExecutiveDrawNow(G);
      if(PMGUI) p_glutSwapBuffers();
    }
    if(!I->CacheSave) {
      if(I->Image[i])
        mfree(I->Image[i]);
      I->Image[i]=NULL;
    }
  }
  return result;
}

void MovieCopyFinish(PyMOLGlobals *G) 
{
  register CMovie *I=G->Movie;
  SceneDirty(G); /* important */
  SettingSet(G,cSetting_cache_frames,(float)I->CacheSave);
  MoviePlay(G,cMovieStop);
  if(!I->CacheSave) {
    MovieClearImages(G); 
    SceneSuppressMovieFrame(G);
  }
}

int MovieLocked(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  return(I->Locked);
}

void MovieSetLock(PyMOLGlobals *G,int lock)
{
  register CMovie *I=G->Movie;
  I->Locked=lock;
}
/*========================================================================*/
void MovieDump(PyMOLGlobals *G)
{
  int a;
  register CMovie *I=G->Movie;
  int flag=false;
  char buffer[OrthoLineLength+100];

  for(a=0;a<I->NFrame;a++) {
    if(I->Cmd[a][0]) {
      flag=true;
      break;
    }
  }
  if(flag&&I->NFrame) {
    PRINTFB(G,FB_Movie,FB_Results)
      " Movie: General Purpose Commands:\n"
      ENDFB(G);
    for(a=0;a<I->NFrame;a++) {
      if(I->Cmd[a][0]) {
        sprintf(buffer,"%5d: %s\n",a+1,I->Cmd[a]);
        OrthoAddOutput(G,buffer);
      }
    }
  } else {
    PRINTFB(G,FB_Movie,FB_Results)
      " Movie: No movie commands are defined.\n"
      ENDFB(G);
  }
}
/*========================================================================*/
static int MovieCmdFromPyList(PyMOLGlobals *G,PyObject *list,int *warning)
{
  register CMovie *I=G->Movie;
  int ok=true;
  int a;
  int warn=false;

  if(ok) ok=(list!=NULL);
  if(ok) ok=PyList_Check(list);
  
  for(a=0;a<I->NFrame;a++) {
    if(ok) ok=PConvPyStrToStr(PyList_GetItem(list,a),I->Cmd[a],OrthoLineLength);
    if(ok) warn= (warn||I->Cmd[a][0]);
  }
  *warning=warn;
  return(ok);
}
/*========================================================================*/
int MovieFromPyList(PyMOLGlobals *G,PyObject *list,int *warning)
{
  int ok=true;
  register CMovie *I=G->Movie;
  int ll = 0;

  MovieReset(G);
  if(ok) ok=PyList_Check(list);
  if(ok) ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
   Always check ll when adding new PyList_GetItem's */
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,0),&I->NFrame);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,1),&I->MatrixFlag);
  if(ok&&I->MatrixFlag) 
    ok=PConvPyListToFloatArrayInPlace(PyList_GetItem(list,2),I->Matrix,cSceneViewSize);
  if(ok) ok=PConvPyIntToInt(PyList_GetItem(list,3),&I->Playing);
  if(ok&&I->NFrame) {
    I->Sequence=VLACalloc(int,I->NFrame);
    I->Cmd=VLACalloc(MovieCmdType,I->NFrame);
    if(ok) ok=PConvPyListToIntArrayInPlace(PyList_GetItem(list,4),I->Sequence,I->NFrame);
    if(ok) ok=MovieCmdFromPyList(G,PyList_GetItem(list,5),warning);
    if((*warning)&&Security) {
      MovieSetLock(G,true);
    }
  }
  if(ok&&(ll>6)) {
    PyObject *tmp;
    VLAFreeP(I->ViewElem);
    I->ViewElem = NULL;
    tmp = PyList_GetItem(list,6);
    if(tmp && !(tmp == Py_None))
      ok = ViewElemVLAFromPyList(tmp,&I->ViewElem,I->NFrame);
  }
  if(!ok) {
    MovieReset(G);
  }
  return(ok);
}
/*========================================================================*/
static PyObject *MovieCmdAsPyList(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  PyObject *result=NULL;
  int a;

  result = PyList_New(I->NFrame);
  if(result) 
    for(a=0;a<I->NFrame;a++) {
      PyList_SetItem(result,a,PyString_FromString(I->Cmd[a]));
    }
  return(PConvAutoNone(result));
}

/*========================================================================*/
PyObject *MovieAsPyList(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  PyObject *result = NULL;

  result = PyList_New(7);
  PyList_SetItem(result,0,PyInt_FromLong(I->NFrame));
  PyList_SetItem(result,1,PyInt_FromLong(I->MatrixFlag));
  PyList_SetItem(result,2,PConvFloatArrayToPyList(I->Matrix,cSceneViewSize));
  PyList_SetItem(result,3,PyInt_FromLong(I->Playing));
  if(I->Sequence) {
    PyList_SetItem(result,4,PConvIntArrayToPyList(I->Sequence,I->NFrame));
  } else {
    PyList_SetItem(result,4,PConvAutoNone(NULL));
  }
  if(I->Cmd) {
    PyList_SetItem(result,5,MovieCmdAsPyList(G));
  } else {
    PyList_SetItem(result,5,PConvAutoNone(NULL));
  }
  if(I->ViewElem) {
    PyList_SetItem(result,6,ViewElemVLAAsPyList(I->ViewElem,I->NFrame));
  } else {
    PyList_SetItem(result,6,PConvAutoNone(NULL));
  }

/*   ImageType *Image;
  int *Sequence;
  MovieCmdType *Cmd;
  int NImage,NFrame;
  unsigned Width,Height;
  int MatrixFlag;
  float Matrix[16];
  int Playing;
*/
  return(PConvAutoNone(result));
}
/*========================================================================*/
int MoviePlaying(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  if(I->Locked)
    return false;
  return(I->Playing);
}
/*========================================================================*/
void MoviePlay(PyMOLGlobals *G,int cmd)
{
  register CMovie *I=G->Movie;
  switch(cmd) {
  case cMovieStop:
	 I->Playing=false;
	 break;
  case cMoviePlay:
    if(!(int)SettingGet(G,cSetting_movie_loop)) { 
      /* if not looping, and at end of movie, then automatically rewind
       and force execution of the first movie command */
      if((SettingGetGlobal_i(G,cSetting_frame))==(SceneGetNFrame(G))) {
        SceneSetFrame(G,7,0);
      }
    }
	 I->Playing=true;
	 break;
  }
  OrthoDirty(G);
  SceneRestartTimers(G);
}
/*========================================================================*/
int  MovieMatrix(PyMOLGlobals *G,int action)
{
  register CMovie *I=G->Movie;
  int result = false;
  switch(action) {
  case cMovieMatrixClear:
	 I->MatrixFlag=false;
    result = 1;
	 break;
  case cMovieMatrixStore:
    SceneGetView(G,I->Matrix);
	 I->MatrixFlag=true;
    result = 1;
	 break;
  case cMovieMatrixRecall:
	 if(I->MatrixFlag) 
		SceneSetView(G,I->Matrix,true);
    else
      result = 0;
	 break;
  case cMovieMatrixCheck:
    result = I->MatrixFlag;
    break;
  }
  return result;
}
/*========================================================================*/
void MovieSetSize(PyMOLGlobals *G,unsigned int width,unsigned int height)
{  
  register CMovie *I=G->Movie;
  I->Width=width;
  I->Height=height;
}
/*========================================================================*/
int MoviePNG(PyMOLGlobals *G,char *prefix,int save,int start,int stop)
{
  /* assumed locked api, blocked threads, and master thread on entry */

  /* writes the current movie sequence to a set of PNG files 
  * this routine can take a LOT of time -- 
  * TODO: develop an interrupt mechanism */

  register CMovie *I=G->Movie;
  int a;
  int i;
  char fname[255];
  char buffer[255];
  int nFrame;

  save = (int)SettingGet(G,cSetting_cache_frames); 
  if(!save)
    MovieClearImages(G);
  SettingSet(G,cSetting_cache_frames,1.0);
  OrthoBusyPrime(G);
  nFrame = I->NFrame;
  if(!nFrame) {
	 nFrame=SceneGetNFrame(G);
  }
  if(start<0) start=0;
  if(start>nFrame) start=nFrame;
  if(stop<0) stop=nFrame;
  if(stop>nFrame) stop=nFrame;
  sprintf(buffer,"Creating movie (%d frames)...",nFrame);
  OrthoBusyMessage(G,buffer);
  if((start!=0)||(stop!=(nFrame+1)))
    
  SceneSetFrame(G,0,0);
  MoviePlay(G,cMoviePlay);
  VLACheck(I->Image,ImageType,nFrame);

  OrthoBusySlow(G,0,nFrame);
  for(a=0;a<nFrame;a++)
	 {
      PRINTFB(G,FB_Movie,FB_Debugging)
        " MoviePNG-DEBUG: Cycle %d...\n",a
        ENDFB(G);
		sprintf(fname,"%s%04d.png",prefix,a+1);
		SceneSetFrame(G,0,a);
		MovieDoFrameCommand(G,a);
		PFlush();
		i=MovieFrameToImage(G,a);
      VLACheck(I->Image,ImageType,i);
      if((a>=start)&&(a<=stop)) { /* only render frames in the specified interval */
        if(!I->Image[i]) {
          SceneMakeMovieImage(G);
        }
        if(!I->Image[i]) {
          PRINTFB(G,FB_Movie,FB_Errors) 
            "MoviePNG-Error: Missing rendered image.\n"
            ENDFB(G);
        } else {
          MyPNGWrite(G,fname,I->Image[i],I->Width,I->Height);		
          ExecutiveDrawNow(G);
          OrthoBusySlow(G,a,nFrame);
          if(PMGUI) p_glutSwapBuffers();
          PRINTFB(G,FB_Movie,FB_Debugging)
            " MoviePNG-DEBUG: i = %d, I->Image[i] = %p\n",i,I->Image[i]
            ENDFB(G);
          if(Feedback(G,FB_Movie,FB_Actions)) {
            printf(" MoviePNG: wrote %s\n",fname);
          }
        }
      }
      if(I->Image[i])
        mfree(I->Image[i]);
		I->Image[i]=NULL;
	 }
  SceneDirty(G); /* important */
  PRINTFB(G,FB_Movie,FB_Debugging)
    " MoviePNG-DEBUG: done.\n"
    ENDFB(G);
  SettingSet(G,cSetting_cache_frames,(float)save);
  MoviePlay(G,cMovieStop);
  MovieClearImages(G);
  SceneSuppressMovieFrame(G);
  return(true);
}
/*========================================================================*/
void MovieAppendSequence(PyMOLGlobals *G,char *str,int start_from)
{
  register CMovie *I=G->Movie;
  int c=0;
  int i;
  char *s,number[20];
 
  if(start_from<0)
    start_from = I->NFrame;

  c=start_from;

  PRINTFB(G,FB_Movie,FB_Debugging)
    " MovieSequence: entered. str:%s\n",str
    ENDFB(G);

  s=str;
  while(*s) {
    s=ParseWord(number,s,20);
	 if(sscanf(number,"%i",&i)) { /* slow */
		c++;
	 }
  }

  if(!c) {
    VLAFreeP(I->Sequence); 
    VLAFreeP(I->Cmd);
    VLAFreeP(I->ViewElem);
    I->NFrame=0;
  } else {
    if(!I->Sequence) {
      I->Sequence = VLACalloc(int,c);
    } else {
      VLASize(I->Sequence, int, start_from); /* to clear */
      VLASize(I->Sequence, int, c);
    }
    if(!I->Cmd) {
      I->Cmd = VLACalloc(MovieCmdType, c);
    } else {
      VLASize(I->Cmd, MovieCmdType, start_from);      
      VLASize(I->Cmd, MovieCmdType, c);
    }
    if(!I->ViewElem) {
      I->ViewElem = VLACalloc(CViewElem, c);    
    } else {
      VLASize(I->ViewElem, CViewElem, start_from);
      VLASize(I->ViewElem, CViewElem, c);
    }
  }

  if(c&&str[0]) { /* not just a reset */
    for(i=start_from;i<c;i++)
      I->Cmd[i][0]=0;
    c=start_from;
    s=str;
    while(*s) {
      s=ParseWord(number,s,20);
      if(sscanf(number,"%i",&I->Sequence[c])) {
        c++;
      }
    }
    I->NFrame=c;
  } else if(!str[0]) {
    I->NFrame=start_from;
  }
           
  VLACheck(I->Image,ImageType,I->NFrame);
  PRINTFB(G,FB_Movie,FB_Debugging)
    " MovieSequence: leaving... I->NFrame%d\n",I->NFrame
    ENDFB(G);

}
/*========================================================================*/
int MovieFrameToImage(PyMOLGlobals *G,int frame)
{
  int result = 0;
  int single_image = (int)SettingGet(G,cSetting_single_image);
  if(single_image)
	 result = MovieFrameToIndex(G,frame);
  else
    result = frame;
  PRINTFB(G,FB_Movie,FB_Debugging)
    " MovieFrameToImage-DEBUG: result %d\n",result
    ENDFB(G);
  return(result);
}
/*========================================================================*/
int MovieFrameToIndex(PyMOLGlobals *G,int frame)
{
  register CMovie *I=G->Movie;
  if(I->Sequence&&I->NFrame) {
	 if(frame<I->NFrame)
		return(I->Sequence[frame]);
	 else 
		return(I->Sequence[I->NFrame-1]);
  } else return(frame);
}
/*========================================================================*/
void MovieSetImage(PyMOLGlobals *G,int index,ImageType image)
{
  register CMovie *I=G->Movie;

  PRINTFB(G,FB_Movie,FB_Blather)
    " MovieSetImage: setting movie image %d\n",index+1
    ENDFB(G);

  VLACheck(I->Image,ImageType,index);
  if(I->Image[index]) FreeP(I->Image[index]);
  I->Image[index]=image;
  if(I->NImage<(index+1))
	 I->NImage=index+1;
}
/*========================================================================*/
void MovieDoFrameCommand(PyMOLGlobals *G,int frame)
{
  register CMovie *I=G->Movie;
  if(frame==0)
	 MovieMatrix(G,cMovieMatrixRecall);
  if(!I->Locked) {
    if((frame>=0)&&(frame<I->NFrame))
      {
        if(I->Cmd[frame][0]) 
          PParse(I->Cmd[frame]);
        if(I->ViewElem) 
          SceneFromViewElem(G,I->ViewElem+frame);
      }
  }
}
/*========================================================================*/
void MovieSetCommand(PyMOLGlobals *G,int frame,char *command)
{
  register CMovie *I=G->Movie;
  int a,len;
  if((frame>=0)&&(frame<I->NFrame))
	 {
		len=strlen(command);
		if(len>(sizeof(MovieCmdType)-1))
		  len=sizeof(MovieCmdType)-1;
		for(a=0;a<len;a++)
		  I->Cmd[frame][a]=command[a];
		I->Cmd[frame][len]=0;
	 }
  else {
    PRINTFB(G,FB_Movie,FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",frame+1
      ENDFB(G);
  }
}

static int interpolate_view(CViewElem *first,CViewElem *last,float power,float bias)
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
  float tLinear = true; /* always do linear for now... */
  float pivot[3];
  const float _1 = 1.0F, _p5 = 0.5F;
  int parabolic = true;

  if(power<0.0F) {
    parabolic = false;
    power = -power;
  }
  
  /* I have no clue whether we're column or row major at this
     point...but it hardly matters!  */

  copy44d33f(first->matrix,first3x3);
  copy44d33f(last->matrix,last3x3);
  transpose33f33f(first3x3,inverse3x3);

  copy3f(first->pre,&rot[3][0]);
  copy3f(last->pre,&rot[4][0]);
  multiply33f33f(inverse3x3,last3x3,&rot[0][0]);

  matrix_to_rotation(rot,rot_axis,&angle);

  if(!tLinear) {
    subtract3f(last->pre,first->pre,tVector);
    normalize23f(tVector,tDir);
    tLen = length3f(tVector); 
    
    average3f(last->pre,first->pre,tCenter); /* center of translation */
    
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
  
    tAngle = acos(dot_product3f(v1,v2));
    
    if(fabs(tAngle)>fabs(angle))
      tAngle = fabs(angle)*(tAngle/fabs(angle));
    
    /* if translation angle > rotation angle then sets translation angle 
     * to same as rotation angle, with proper sign of course */
    
    if (tAngle>0.0001) 
      {
        pLen = tan(tAngle/2); 
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

  current = first+1;
  
  for(a=0;a<n;a++) {
    float fxn = ((float)a+1)/(n+1);
    float fxn_1;

    /*    printf("Fxn: %8.3f ",fxn);*/
    if(bias!=1.0F) {
      fxn = 1-pow(1-pow(fxn,bias),_1/bias);
    }
    /*    printf("%8.3f bias %8.3f\n",fxn,bias);*/

    if(power!=1.0F) {
      if(fxn<0.5F) {
        if(parabolic) 
          fxn = (float)pow(fxn*2.0F,power)*_p5; /* parabolic */
        else
          fxn = (_1-pow(_1-pow((fxn*2.0F),power),_1/power))*_p5; /* circular */
      } else if(fxn>0.5F) {
        fxn = _1 - fxn;
        if(parabolic) 
          fxn = (float)pow(fxn*2.0F,power)*_p5; /* parabolic */
        else
          fxn = (_1-pow(_1-pow((fxn*2.0F),power),_1/power))*_p5; /* circular */
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
    current++;
  }

  return 1;
}
/*========================================================================*/
int MovieView(PyMOLGlobals *G,int action,int first,int last,float power,float bias)
{
  register CMovie *I=G->Movie;
  int frame;
  switch(action) {
  case 0: /* set */
    if(I->ViewElem) {
      if(first<0)
        first = SceneGetFrame(G);
      if(last<0)
        last = first;
      for(frame=first;frame<=last;frame++) {
        if((frame>=0)&&(frame<I->NFrame)) {
          VLACheck(I->ViewElem,CViewElem,frame);
          PRINTFB(G,FB_Movie,FB_Details)
            " MovieView: Setting frame %d.\n",frame+1
            ENDFB(G);
          SceneToViewElem(G,I->ViewElem+frame);          
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
        if((frame>=0)&&(frame<I->NFrame)) {
          VLACheck(I->ViewElem,CViewElem,frame);
          UtilZeroMem((void*)(I->ViewElem+frame),sizeof(CViewElem));
        }
      }
    }
    break;
  case 2: /* interpolate & reinterpolate */
  case 3:
    {
      CViewElem *first_view=NULL,*last_view=NULL;
      if(first<0)
        first = 0;
      if(last<0)
        last = SceneGetNFrame(G)-1;

      VLACheck(I->ViewElem,CViewElem,last);
      if(action==2) {
        PRINTFB(G,FB_Movie,FB_Details)
          " MovieView: interpolating unspecified frames %d to %d.\n",first+1,last+1
          ENDFB(G);
      } else {
        PRINTFB(G,FB_Movie,FB_Details)
          " MovieView: reinterpolating all frames %d to %d.\n",first+1,last+1
          ENDFB(G);
      }
      for(frame=first;frame<=last;frame++) {
        if((frame>=0)&&(frame<I->NFrame)) {
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
                interpolate_view(first_view,last_view,power,bias);
              }
              first_view = last_view;
              last_view = NULL;
            }
          }
        }
      }
    }
    break;
  }
  return 1;
}
/*========================================================================*/
void MovieAppendCommand(PyMOLGlobals *G,int frame,char *command)
{
  register CMovie *I=G->Movie;
  int a,len,cur_len;
  if((frame>=0)&&(frame<I->NFrame))
	 {
		len=strlen(command);
      cur_len=strlen(I->Cmd[frame]);
		if((unsigned)len>(sizeof(MovieCmdType)+cur_len-1))
        len=sizeof(MovieCmdType)+cur_len-1;
      for(a=0;a<len;a++)
        I->Cmd[frame][cur_len+a]=command[a];
		I->Cmd[frame][cur_len+len]=0;
	 }
  else {
    PRINTFB(G,FB_Movie,FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",frame+1
      ENDFB(G);
  }
}
/*========================================================================*/
ImageType MovieGetImage(PyMOLGlobals *G,int index)
{
  register CMovie *I=G->Movie;
  if((index>=0)&&(index<I->NImage))
	 return(I->Image[index]);
  else
	 return(NULL);
}
/*========================================================================*/
int MovieDefined(PyMOLGlobals *G) 
{
  register CMovie *I=G->Movie;
  return(I->NFrame>0);
}
/*========================================================================*/
int MovieGetLength(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  int len;
  if(!I->NFrame)
	 len = -I->NImage;
  else
	 len = I->NFrame;
  PRINTFD(G,FB_Movie)
    " MovieGetLength: leaving...result %d\n",len
    ENDFD;
  return(len);
}
/*========================================================================*/
void MovieClearImages(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  int a;

  PRINTFB(G,FB_Movie,FB_Blather)
    " MovieClearImages: clearing...\n"
    ENDFB(G);
  for(a=0;a<I->NImage;a++)
	 {
		if(I->Image[a]) {
		  FreeP(I->Image[a]);
		  I->Image[a]=NULL;
		}
	 }
  I->NImage=0;
  SceneDirty(G);
}
/*========================================================================*/
void MovieReset(PyMOLGlobals *G) {
  register CMovie *I=G->Movie;
  MovieClearImages(G);

  VLAFreeP(I->Cmd);  
  VLAFreeP(I->Sequence);  
  VLAFreeP(I->ViewElem);  

  I->NFrame=0;
  I->MatrixFlag=false;
  I->Locked=false;
  I->Playing=false;
}
/*========================================================================*/
void MovieFree(PyMOLGlobals *G)
{
  register CMovie *I=G->Movie;
  MovieClearImages(G);
  VLAFree(I->Image);
  VLAFreeP(I->ViewElem);
  VLAFreeP(I->Cmd);
  VLAFreeP(I->Sequence);
  FreeP(G->Movie);
}
/*========================================================================*/
int MovieInit(PyMOLGlobals *G)
{
  register CMovie *I=NULL;

  if( (I=(G->Movie=Calloc(CMovie,1)))) {
    int a;
    I->Playing=false;
    I->Image=VLAMalloc(10,sizeof(ImageType),5,true); /* auto-zero */
    I->Sequence=NULL;
    I->Cmd=NULL;
    I->ViewElem = NULL;
    I->NImage=0;
    I->NFrame=0;
    for(a=0;a<16;a++)
      I->Matrix[a]=0.0F;
    I->MatrixFlag=false;
    return 1;
  } else {
    return 0;
  }
}




