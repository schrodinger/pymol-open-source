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
#ifndef _H_Movie
#define _H_Movie

#include"os_python.h"

#include"Ortho.h"
#include"Scene.h"
#include"View.h"

typedef unsigned char *ImageType;
typedef char MovieCmdType[OrthoLineLength];

typedef struct  {
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
} CMovie;

int MovieFromPyList(PyObject *list,int *warning);
PyObject *MovieAsPyList(void);

void MovieInit(void);
void MovieFree(void);
void MovieReset(void);
void MovieDump(void);
void MovieSequence(char *seq);
int MoviePNG(char *prefix,int save,int start,int stop);
void MovieSetCommand(int frame,char *command);
void MovieAppendCommand(int frame,char *command);

void MovieDoFrameCommand(int frame);

void MovieCopyPrepare(int *width,int *height,int *length);
int MovieCopyFrame(int frame,int width,int height,int rowbytes,void *ptr);
void MovieCopyFinish(void);

#define cMovieStop 0
#define cMoviePlay 1

void MoviePlay(int cmd);
int MoviePlaying(void);
void MovieSetSize(unsigned int width,unsigned int height);

void MovieClearImages(void);
ImageType MovieGetImage(int image);
void MovieSetImage(int index,ImageType image);

int MovieGetLength(void);
int MovieFrameToImage(int frame);
int MovieFrameToIndex(int frame);
int MovieLocked(void);
void MovieSetLock(int);
int MovieDefined(void);
int MovieView(int action,int first,int last,float power);

#define cMovieMatrixClear  0
#define cMovieMatrixStore  1
#define cMovieMatrixRecall 2
#define cMovieMatrixCheck  3

int MovieMatrix(int action); /* 0 clear, 1 remember, 2 recall */

/*void MovieSave(char *fname);
  void MovieLoad(char *fname);*/

#endif
