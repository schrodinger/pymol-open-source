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

CMovie Movie;

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
  for(a=0;a<I->NFrame;a++) {
    if(I->Cmd[a][0]) {
      flag=true;
      break;
    }
  }
  if(flag&&I->NFrame) {
    PRINTFB(FB_Movie,FB_Results)
      " Movie commands:\n"
      ENDFB;
    for(a=0;a<I->NFrame;a++) {
      if(I->Cmd[a][0]) {
        PRINTFB(FB_Movie,FB_Results)
          "%5d: %s\n",a+1,I->Cmd[a]
          ENDFB;
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
    if(ok) warn=warn||I->Cmd[a][0];
  }
  *warning=warn;
  return(ok);
}
/*========================================================================*/
int MovieFromPyList(PyObject *list,int *warning)
{
  int ok=true;
  CMovie *I=&Movie;

  MovieReset();
  if(ok) ok=PyList_Check(list);
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
    if(warning&&Security) {
      MovieSetLock(true);
      PParse("wizard security");
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
	 I->Playing=true;
	 break;
  }
  SceneRestartTimers();
}
/*========================================================================*/
void MovieMatrix(int action)
{
  CMovie *I=&Movie;

  switch(action) {
  case cMovieMatrixClear:
	 I->MatrixFlag=false;
	 break;
  case cMovieMatrixStore:
    SceneGetView(I->Matrix);
	 I->MatrixFlag=true;
	 break;
  case cMovieMatrixRecall:
	 if(I->MatrixFlag) 
		SceneSetView(I->Matrix);
	 break;
  }
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

  save = SettingGet(cSetting_cache_frames); 
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
		sprintf(fname,"%s_%04d.png",prefix,a+1);
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
  SettingSet(cSetting_cache_frames,save);
  MoviePlay(cMovieStop);
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
  I->NFrame=0;
  if(str[0]) { /* not just a reset */
    I->Sequence=Alloc(int,c+1);
    I->Cmd=Alloc(MovieCmdType,c+1);
    
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
		if(len>(sizeof(MovieCmdType)+cur_len-1))
        len=sizeof(MovieCmdType)+cur_len-1;
      for(a=0;a<len;a++)
        I->Cmd[frame][cur_len+a]=command[a];
		I->Cmd[frame][cur_len+len]=0;
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
  I->NImage=0;
  I->NFrame=0;
  for(a=0;a<16;a++)
    I->Matrix[a]=0.0F;
  I->MatrixFlag=false;
}


