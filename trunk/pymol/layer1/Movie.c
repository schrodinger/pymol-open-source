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


CMovie Movie;

void MovieCopyPrepare(int *width,int *height,int *length)
{
  /* assumed locked api, blocked threads, and master thread on entry */

  /* writes the current movie sequence to a set of PNG files 
  * this routine can take a LOT of time -- 
  * TODO: develop an interrupt mechanism */

  int start,stop;
  CMovie *I=&Movie;
  int nFrame;

  I->CacheSave = (int)SettingGet(cSetting_cache_frames); 
  if(!I->CacheSave)
    MovieClearImages();
  SettingSet(cSetting_cache_frames,1.0);
  nFrame = I->NFrame;
  if(!nFrame) {
	 nFrame=SceneGetNFrame();
  }
  start=0;
  stop=nFrame;
  if((start!=0)||(stop!=(nFrame+1)))
    SceneSetFrame(0,0);
  MoviePlay(cMoviePlay);
  VLACheck(I->Image,ImageType,nFrame);
  SceneGetWidthHeight(width,height);
  *length = nFrame;
}

int MovieCopyFrame(int frame,int width,int height,int rowbytes,void *ptr)
{
  CMovie *I=&Movie;
  int result=false;
  int nFrame;
  
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame=SceneGetNFrame();
  }
  
  if((width==I->Width)&&(height==I->Height)&&(frame<nFrame)&&(ptr)) {
    int a = frame;
    int i;
    SceneSetFrame(0,a);
    MovieDoFrameCommand(a);
    PFlush();
    i=MovieFrameToImage(a);
    VLACheck(I->Image,ImageType,i);
    if(!I->Image[i]) {
      SceneMakeMovieImage();
    }
    if(!I->Image[i]) {
      PRINTFB(FB_Movie,FB_Errors) 
        "MoviePNG-Error: Missing rendered image.\n"
        ENDFB;
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
      ExecutiveDrawNow();
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

void MovieCopyFinish(void) 
{
  CMovie *I=&Movie;
  SceneDirty(); /* important */
  SettingSet(cSetting_cache_frames,(float)I->CacheSave);
  MoviePlay(cMovieStop);
  if(!I->CacheSave) {
    MovieClearImages(); 
    SceneSuppressMovieFrame();
  }
}

int MovieLocked(void)
{
  CMovie *I=&Movie;
  return(I->Locked);
}

void MovieSetLock(int lock)
{
  CMovie *I=&Movie;
  I->Locked=lock;
}
/*========================================================================*/
void MovieDump(void)
{
  int a;
  CMovie *I=&Movie;
  int flag=false;
  char buffer[OrthoLineLength+100];

  for(a=0;a<I->NFrame;a++) {
    if(I->Cmd[a][0]) {
      flag=true;
      break;
    }
  }
  if(flag&&I->NFrame) {
    PRINTFB(FB_Movie,FB_Results)
      " Movie: General Purpose Commands:\n"
      ENDFB;
    for(a=0;a<I->NFrame;a++) {
      if(I->Cmd[a][0]) {
        sprintf(buffer,"%5d: %s\n",a+1,I->Cmd[a]);
        OrthoAddOutput(buffer);
      }
    }
  } else {
    PRINTFB(FB_Movie,FB_Results)
      " Movie: No movie commands are defined.\n"
      ENDFB;
  }
}
/*========================================================================*/
static int MovieCmdFromPyList(PyObject *list,int *warning)
{
  CMovie *I=&Movie;
  int ok=true;
  int a;
  int warn=false;

  if(ok) ok=PyList_Check(list);
  
  for(a=0;a<I->NFrame;a++) {
    if(ok) ok=PConvPyStrToStr(PyList_GetItem(list,a),I->Cmd[a],OrthoLineLength);
    if(ok) warn= (warn||I->Cmd[a][0]);
  }
  *warning=warn;
  return(ok);
}
/*========================================================================*/
int MovieFromPyList(PyObject *list,int *warning)
{
  int ok=true;
  CMovie *I=&Movie;
  int ll;

  MovieReset();
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
    I->Sequence=Alloc(int,I->NFrame+1);
    I->Cmd=Alloc(MovieCmdType,I->NFrame+1);
    if(ok) ok=PConvPyListToIntArrayInPlace(PyList_GetItem(list,4),I->Sequence,I->NFrame);
    if(ok) ok=MovieCmdFromPyList(PyList_GetItem(list,5),warning);
    if((*warning)&&Security) {
      MovieSetLock(true);
    }
  }
  if(!ok) {
    MovieReset();
  }
  return(ok);
}
/*========================================================================*/
static PyObject *MovieCmdAsPyList(void)
{
  CMovie *I=&Movie;
  PyObject *result=NULL;
  int a;

  result = PyList_New(I->NFrame);
  for(a=0;a<I->NFrame;a++) {
    PyList_SetItem(result,a,PyString_FromString(I->Cmd[a]));
  }
  return(PConvAutoNone(result));
}

/*========================================================================*/
PyObject *MovieAsPyList(void)
{
  CMovie *I=&Movie;
  PyObject *result = NULL;

  result = PyList_New(6);
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
    PyList_SetItem(result,5,MovieCmdAsPyList());
  } else {
    PyList_SetItem(result,5,PConvAutoNone(NULL));
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
int MoviePlaying(void)
{
  CMovie *I=&Movie;
  if(I->Locked)
    return false;
  return(I->Playing);
}
/*========================================================================*/
void MoviePlay(int cmd)
{
  CMovie *I=&Movie;
  switch(cmd) {
  case cMovieStop:
	 I->Playing=false;
	 break;
  case cMoviePlay:
    if(!(int)SettingGet(cSetting_movie_loop)) { 
      /* if not looping, and at end of movie, then automatically rewind
       and force execution of the first movie command */
      if((SettingGetGlobal_i(cSetting_frame))==(SceneGetNFrame())) {
        SceneSetFrame(7,0);
      }
    }
	 I->Playing=true;
	 break;
  }
  OrthoDirty();
  SceneRestartTimers();
}
/*========================================================================*/
int  MovieMatrix(int action)
{
  CMovie *I=&Movie;
  int result = false;
  switch(action) {
  case cMovieMatrixClear:
	 I->MatrixFlag=false;
    result = 1;
	 break;
  case cMovieMatrixStore:
    SceneGetView(I->Matrix);
	 I->MatrixFlag=true;
    result = 1;
	 break;
  case cMovieMatrixRecall:
	 if(I->MatrixFlag) 
		SceneSetView(I->Matrix,true);
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
void MovieSetSize(unsigned int width,unsigned int height)
{  
  CMovie *I=&Movie;
  I->Width=width;
  I->Height=height;
}
/*========================================================================*/
int MoviePNG(char *prefix,int save,int start,int stop)
{
  /* assumed locked api, blocked threads, and master thread on entry */

  /* writes the current movie sequence to a set of PNG files 
  * this routine can take a LOT of time -- 
  * TODO: develop an interrupt mechanism */

  CMovie *I=&Movie;
  int a;
  int i;
  char fname[255];
  char buffer[255];
  int nFrame;

  save = (int)SettingGet(cSetting_cache_frames); 
  if(!save)
    MovieClearImages();
  SettingSet(cSetting_cache_frames,1.0);
  OrthoBusyPrime();
  nFrame = I->NFrame;
  if(!nFrame) {
	 nFrame=SceneGetNFrame();
  }
  if(start<0) start=0;
  if(start>nFrame) start=nFrame;
  if(stop<0) stop=nFrame;
  if(stop>nFrame) stop=nFrame;
  sprintf(buffer,"Creating movie (%d frames)...",nFrame);
  OrthoBusyMessage(buffer);
  if((start!=0)||(stop!=(nFrame+1)))
    
  SceneSetFrame(0,0);
  MoviePlay(cMoviePlay);
  VLACheck(I->Image,ImageType,nFrame);

  OrthoBusySlow(0,nFrame);
  for(a=0;a<nFrame;a++)
	 {
      PRINTFB(FB_Movie,FB_Debugging)
        " MoviePNG-DEBUG: Cycle %d...\n",a
        ENDFB;
		sprintf(fname,"%s%04d.png",prefix,a+1);
		SceneSetFrame(0,a);
		MovieDoFrameCommand(a);
		PFlush();
		i=MovieFrameToImage(a);
      VLACheck(I->Image,ImageType,i);
      if((a>=start)&&(a<=stop)) { /* only render frames in the specified interval */
        if(!I->Image[i]) {
          SceneMakeMovieImage();
        }
        if(!I->Image[i]) {
          PRINTFB(FB_Movie,FB_Errors) 
            "MoviePNG-Error: Missing rendered image.\n"
            ENDFB;
        } else {
          MyPNGWrite(fname,I->Image[i],I->Width,I->Height);		
          ExecutiveDrawNow();
          OrthoBusySlow(a,nFrame);
          if(PMGUI) p_glutSwapBuffers();
          PRINTFB(FB_Movie,FB_Debugging)
            " MoviePNG-DEBUG: i = %d, I->Image[i] = %p\n",i,I->Image[i]
            ENDFB;
          if(Feedback(FB_Movie,FB_Actions)) {
            printf(" MoviePNG: wrote %s\n",fname);
          }
        }
      }
      if(I->Image[i])
        mfree(I->Image[i]);
		I->Image[i]=NULL;
	 }
  SceneDirty(); /* important */
  PRINTFB(FB_Movie,FB_Debugging)
    " MoviePNG-DEBUG: done.\n"
    ENDFB;
  SettingSet(cSetting_cache_frames,(float)save);
  MoviePlay(cMovieStop);
  MovieClearImages();
  SceneSuppressMovieFrame();
  return(true);
}
/*========================================================================*/
void MovieSequence(char *str)
{
  CMovie *I=&Movie;
  int c=0;
  int i;
  char *s;
  c=0;

  PRINTFB(FB_Movie,FB_Debugging)
    " MovieSequence: entered. str:%s\n",str
    ENDFB;

  s=str;
  while(*s) {
	 if(sscanf(s,"%i",&i)) {
		c++;
	 }
	 s++;
	 while(*s) {
		if(*s==' ') break;
		s++;
	 }
  }
  FreeP(I->Sequence);
  FreeP(I->Cmd);
  VLAFreeP(I->ViewElem);
  I->NFrame=0;
  if(str[0]) { /* not just a reset */
    I->Sequence=Alloc(int,c+1);
    I->Cmd=Alloc(MovieCmdType,c+1);
    I->ViewElem=VLACalloc(CViewElem,c+1);
    for(i=0;i<c;i++)
      I->Cmd[i][0]=0;
    c=0;
    s=str;
    while(*s) {
      if(sscanf(s,"%i",&I->Sequence[c])) {
        c++;
      }
      s++;
      while(*s) {
        if(*s==' ') break;
        s++;
      }
    }
    I->Sequence[c]=-1;
    I->NFrame=c;
  }
  VLACheck(I->Image,ImageType,I->NFrame);
  PRINTFB(FB_Movie,FB_Debugging)
    " MovieSequence: leaving... I->NFrame%d\n",I->NFrame
    ENDFB;

}
/*========================================================================*/
int MovieFrameToImage(int frame)
{
  int result = 0;
  int single_image = (int)SettingGet(cSetting_single_image);
  if(single_image)
	 result = MovieFrameToIndex(frame);
  else
    result = frame;
  PRINTFB(FB_Movie,FB_Debugging)
    " MovieFrameToImage-DEBUG: result %d\n",result
    ENDFB;
  return(result);
}
/*========================================================================*/
int MovieFrameToIndex(int frame)
{
  CMovie *I=&Movie;
  if(I->Sequence&&I->NFrame) {
	 if(frame<I->NFrame)
		return(I->Sequence[frame]);
	 else 
		return(I->Sequence[I->NFrame-1]);
  } else return(frame);
}
/*========================================================================*/
void MovieSetImage(int index,ImageType image)
{
  CMovie *I=&Movie;

  PRINTFB(FB_Movie,FB_Blather)
    " MovieSetImage: setting movie image %d\n",index+1
    ENDFB;

  VLACheck(I->Image,ImageType,index);
  if(I->Image[index]) FreeP(I->Image[index]);
  I->Image[index]=image;
  if(I->NImage<(index+1))
	 I->NImage=index+1;
}
/*========================================================================*/
void MovieDoFrameCommand(int frame)
{
  CMovie *I=&Movie;
  if(frame==0)
	 MovieMatrix(cMovieMatrixRecall);
  if(!I->Locked) {
    if((frame>=0)&&(frame<I->NFrame))
      {
        if(I->Cmd[frame][0]) 
          PParse(I->Cmd[frame]);
        if(I->ViewElem) 
          SceneFromViewElem(I->ViewElem+frame);
      }
  }
}
/*========================================================================*/
void MovieSetCommand(int frame,char *command)
{
  CMovie *I=&Movie;
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
    PRINTFB(FB_Movie,FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",frame+1
      ENDFB;
  }
}

static int interpolate_view(CViewElem *first,CViewElem *last,float power)
{
  float first3x3[9];
  float last3x3[9];
  float inverse3x3[9];
  float inter3x3[9];
  float rot_axis[3],trans[3];
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
    
    if(power!=1.0F) {
      if(fxn<0.5F) {
        fxn = (float)pow(fxn*2.0F,power)*0.5F;
      } else if(fxn>0.5F) {
        fxn = 1.0F - fxn;
        fxn = (float)pow(fxn*2.0F,power)*0.5F;
        fxn = 1.0F - fxn;
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
    current->specified = false;
    current++;
  }

  return 1;
}
/*========================================================================*/
int MovieView(int action,int first,int last,float power)
{
  CMovie *I=&Movie;
  int frame;
  switch(action) {
  case 0: /* set */
    if(I->ViewElem) {
      if(first<0)
        first = SceneGetFrame();
      if(last<0)
        last = first;
      for(frame=first;frame<=last;frame++) {
        if((frame>=0)&&(frame<I->NFrame)) {
          VLACheck(I->ViewElem,CViewElem,frame);
          PRINTFB(FB_Movie,FB_Details)
            " MovieView: Setting frame %d.\n",frame+1
            ENDFB;
          SceneToViewElem(I->ViewElem+frame);          
          I->ViewElem[frame].specified = true;
        }
      }
    }
    break;
  case 1: /* clear */
    if(I->ViewElem) {
      if(first<0)
        first = SceneGetFrame();
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
  case 2: /* interpolate */
    {
      CViewElem *first_view=NULL,*last_view=NULL;
      if(first<0)
        first = 0;
      if(last<0)
        last = SceneGetNFrame()-1;

      VLACheck(I->ViewElem,CViewElem,last);
      PRINTFB(FB_Movie,FB_Details)
        " MovieView: interpolating frames %d to %d.\n",first+1,last+1
        ENDFB;
      for(frame=first;frame<=last;frame++) {
        if((frame>=0)&&(frame<I->NFrame)) {
          if(!first_view) {
            if(I->ViewElem[frame].specified) {
              first_view = I->ViewElem + frame;
            }
          } else {
            if(I->ViewElem[frame].specified) {
              last_view = I->ViewElem + frame;
              interpolate_view(first_view,last_view,power);
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
void MovieAppendCommand(int frame,char *command)
{
  CMovie *I=&Movie;
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
    PRINTFB(FB_Movie,FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",frame+1
      ENDFB;
  }
}
/*========================================================================*/
ImageType MovieGetImage(int index)
{
  CMovie *I=&Movie;
  if((index>=0)&&(index<I->NImage))
	 return(I->Image[index]);
  else
	 return(NULL);
}
/*========================================================================*/
int MovieDefined(void) 
{
  CMovie *I=&Movie;
  return(I->NFrame>0);
}
/*========================================================================*/
int MovieGetLength(void)
{
  CMovie *I=&Movie;
  int len;
  if(!I->NFrame)
	 len = -I->NImage;
  else
	 len = I->NFrame;
  PRINTFD(FB_Movie)
    " MovieGetLength: leaving...result %d\n",len
    ENDFD;
  return(len);
}
/*========================================================================*/
void MovieClearImages(void)
{
  CMovie *I=&Movie;
  int a;

  PRINTFB(FB_Movie,FB_Blather)
    " MovieClearImages: clearing...\n"
    ENDFB;
  for(a=0;a<I->NImage;a++)
	 {
		if(I->Image[a]) {
		  FreeP(I->Image[a]);
		  I->Image[a]=NULL;
		}
	 }
  I->NImage=0;
  SceneDirty();
}
/*========================================================================*/
void MovieReset(void) {
  CMovie *I=&Movie;
  MovieClearImages();
  FreeP(I->Cmd);
  FreeP(I->Sequence);
  I->NFrame=0;
  I->MatrixFlag=false;
  I->Locked=false;
  I->Playing=false;
}
/*========================================================================*/
void MovieFree(void)
{
  CMovie *I=&Movie;
  MovieClearImages();
  VLAFree(I->Image);
  VLAFreeP(I->ViewElem);
  FreeP(I->Cmd);
  FreeP(I->Sequence);
}
/*========================================================================*/
void MovieInit(void)
{
  CMovie *I=&Movie;
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
}




