
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

typedef char MovieCmdType[OrthoLineLength];

int MovieFromPyList(PyMOLGlobals * G, PyObject * list, int *warning);
PyObject *MovieAsPyList(PyMOLGlobals * G);
int MovieGetSpecLevel(PyMOLGlobals *G,int frame);
void MovieDrawViewElem(PyMOLGlobals *G, BlockRect *rect,int frames);

int MovieInit(PyMOLGlobals * G);
Block *MovieGetBlock(PyMOLGlobals * G);
void MovieFree(PyMOLGlobals * G);
void MovieReset(PyMOLGlobals * G);
void MovieDump(PyMOLGlobals * G);
void MovieAppendSequence(PyMOLGlobals * G, char *seq, int start_from,int freeze);
int MovieSeekScene(PyMOLGlobals * G, int loop);
int MoviePNG(PyMOLGlobals * G, char *prefix, int save, int start, int stop,
             int missing_only, int modal, int format, int mode, int quiet);
void MovieSetScrollBarFrame(PyMOLGlobals * G, int frame);
void MovieSetCommand(PyMOLGlobals * G, int frame, char *command);
void MovieAppendCommand(PyMOLGlobals * G, int frame, char *command);

void MovieDoFrameCommand(PyMOLGlobals * G, int frame);

void MovieCopyPrepare(PyMOLGlobals * G, int *width, int *height, int *length);
int MovieCopyFrame(PyMOLGlobals * G, int frame, int width, int height, int rowbytes,
                   void *ptr);
int MoviePurgeFrame(PyMOLGlobals * G, int frame);
void MovieCopyFinish(PyMOLGlobals * G);

#define cMovieStop 0
#define cMoviePlay 1
#define cMovieToggle -1

void MoviePlay(PyMOLGlobals * G, int cmd);
int MoviePlaying(PyMOLGlobals * G);
void MovieSetSize(PyMOLGlobals * G, unsigned int width, unsigned int height);

void MovieClearImages(PyMOLGlobals * G);
ImageType *MovieGetImage(PyMOLGlobals * G, int image);
void MovieSetImage(PyMOLGlobals * G, int index, ImageType * image);

int MovieGetLength(PyMOLGlobals * G);
int MovieGetPanelHeight(PyMOLGlobals * G);
int MovieFrameToImage(PyMOLGlobals * G, int frame);
int MovieFrameToIndex(PyMOLGlobals * G, int frame);
int MovieLocked(PyMOLGlobals * G);
void MovieSetLock(PyMOLGlobals * G, int);
int MovieDefined(PyMOLGlobals * G);
int MovieView(PyMOLGlobals * G, int action, int first,
              int last, float power, float bias,
              int simple, float linear, int wrap,
              int hand, int window, int cycles,
              char *scene_name, float scene_cut, int state, int quiet);
void MovieFlushCommands(PyMOLGlobals * G);
void MovieSetRealtime(PyMOLGlobals * G, int realtime);
int MovieGetRealtime(PyMOLGlobals * G);

#define cMovieMatrixClear  0
#define cMovieMatrixStore  1
#define cMovieMatrixRecall 2
#define cMovieMatrixCheck  3

int MovieMatrix(PyMOLGlobals * G, int action);  /* 0 clear, 1 remember, 2 recall */

int MovieViewModify(PyMOLGlobals *G,int action, int index, int count, int target, int freeze, int localize);
void MovieViewReinterpolate(PyMOLGlobals *G);
void MovieViewTrim(PyMOLGlobals *G,int n_frame);

void MoviePrepareDrag(PyMOLGlobals *G, BlockRect * rect, CObject * obj, int mode, int x, int y, int nearest);
int MovieXtoFrame(PyMOLGlobals *G, BlockRect *rect, int frames, int x, int nearest);

/*void MovieSave(char *fname);
  void MovieLoad(char *fname);*/

#endif
