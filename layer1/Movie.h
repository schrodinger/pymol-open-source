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


typedef unsigned char *ImageType;
typedef char MovieCmdType[255];

typedef struct  {
  ImageType *Image;
  int *Sequence;
  MovieCmdType *Cmd;
  int NImage,NFrame;
  unsigned Width,Height;
  int MatrixFlag;
  float Matrix[16];
  int Playing;
} CMovie;

void MovieInit(void);
void MovieFree(void);
void MovieSequence(char *seq);
int MoviePNG(char *prefix,int save,int start,int stop);
void MovieSetCommand(int frame,char *command);
void MovieAppendCommand(int frame,char *command);

void MovieDoFrameCommand(int frame);

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

#define cMovieMatrixClear  0
#define cMovieMatrixStore  1
#define cMovieMatrixRecall 2

void MovieMatrix(int action); /* 0 clear, 1 remember, 2 recall */

/*void MovieSave(char *fname);
  void MovieLoad(char *fname);*/

#endif
