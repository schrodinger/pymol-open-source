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

#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<GL/glut.h>

#include"Base.h"
#include"OOMac.h"
#include"MemoryDebug.h"
#include"Executive.h"
#include"Ortho.h"
#include"Movie.h"
#include"Scene.h"
#include"MyPNG.h"
#include"PUtils.h"
#include"Setting.h"
#include"main.h"

CMovie Movie;

/*========================================================================*/
int MoviePlaying(void)
{
  CMovie *I=&Movie;
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
  int a;
  float *m;

  switch(action) {
  case cMovieMatrixClear:
	 I->MatrixFlag=false;
	 break;
  case cMovieMatrixStore:
	 I->MatrixFlag=true;
	 m=SceneGetMatrix();
	 for(a=0;a<16;a++)
		I->Matrix[a]=m[a];
	 break;
  case cMovieMatrixRecall:
	 if(I->MatrixFlag) 
		SceneSetMatrix(I->Matrix);
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
void MoviePNG(char *prefix,int save)
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

  fflush(stdout);
  save = SettingGet(cSetting_cache_frames);
  SettingSet(cSetting_cache_frames,1.0);
  OrthoBusyPrime();
  sprintf(buffer,"Creating movie (%d frames)...",I->NFrame);
  OrthoBusyMessage(buffer);
  SceneSetFrame(0,0);
  MoviePlay(cMoviePlay);
  for(a=0;a<I->NFrame;a++)
	 {
		OrthoBusySlow(a,I->NFrame);
		sprintf(fname,"%s_%04d.png",prefix,a+1);
		SceneSetFrame(0,a);
		MovieDoFrameCommand(a);
		PFlush();
		i=MovieFrameToImage(a);
		if(!I->Image[i]) {
		  SceneMakeMovieImage();
		}
		if(!I->Image[i])
		  ErrFatal("MoviePNG","Missing rendering movie image!");
		fflush(stdout);		  
		MyPNGWrite(fname,I->Image[i],I->Width,I->Height);		
		ExecutiveDrawNow();
		if(PMGUI) glutSwapBuffers();
		printf(" MoviePNG: wrote %s\n",fname);
		if(!save) {
		  SceneDirty();
		  mfree(I->Image[i]);
		  I->Image[i]=NULL;
		}
	 }
  SettingSet(cSetting_cache_frames,save);
  MoviePlay(cMovieStop);
}
/*========================================================================*/
void MovieSequence(char *str)
{
  CMovie *I=&Movie;
  int c;
  int i;
  char *s;
  c=0;

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
  if(I->Sequence) FreeP(I->Sequence);
  if(I->Cmd) FreeP(I->Cmd);
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
  VLACheck(I->Image,ImageType,I->NFrame);
}
/*========================================================================*/
int MovieFrameToImage(int frame)
{
  int single_image = (int)SettingGet(cSetting_single_image);
  if(single_image)
	 return(MovieFrameToIndex(frame));
  else
	 return(frame);
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
  if((frame>=0)&&(frame<I->NFrame))
	 {
		if(I->Cmd[frame][0]) 
		  PParse(I->Cmd[frame]);
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
  if(!I->NFrame)
	 return(I->NImage);
  else
	 return(I->NFrame);
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
void MovieFree(void)
{
  CMovie *I=&Movie;
  MovieClearImages();
  VLAFree(I->Image);
  FreeP(I->Cmd);
  FreeP(I->Sequence);
  I->Sequence=NULL;
}
/*========================================================================*/
void MovieInit(void)
{
  CMovie *I=&Movie;

  I->Playing=false;
  I->Image=VLAMalloc(10,sizeof(ImageType),5,true); /* auto-zero */
  I->Sequence=NULL;
  I->Cmd=NULL;
  I->NImage=0;
  I->NFrame=0;
  I->MatrixFlag=false;

}


