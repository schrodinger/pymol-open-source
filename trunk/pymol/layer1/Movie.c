
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
#include"os_python.h"

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
#include"PyMOL.h"
#include"ScrollBar.h"
#include"Menu.h"

#define cMovieDragModeMoveKey   1
#define cMovieDragModeInsDel    2
#define cMovieDragModeCopyKey   3
#define cMovieDragModeOblate    4

typedef struct {
  int stage;

  /* input parameters */
  OrthoLineType prefix;
  int save, start, stop, missing_only;
  int modal, mode;

  /* job / local parameters */
  int frame;
  int image;
  int nFrame;
  double accumTiming;
  double timing;
  int complete;
  int file_missing;
  int format;
  int quiet;
  OrthoLineType fname;

} CMovieModal;

struct _CMovie {
  Block *Block;
  ImageType **Image;
  int *Sequence;
  MovieCmdType *Cmd;
  int NImage, NFrame;
  int MatrixFlag;
  SceneViewType Matrix;
  int Playing;
  int Locked;
  int CacheSave;
  int OverlaySave;
  CViewElem *ViewElem;
  int RecursionFlag;
  int RealtimeFlag;
  CMovieModal Modal;
  int Width, Height;
  struct CScrollBar *ScrollBar;
  int DragMode;
  int Dragging;
  CObject *DragObj; /* if not dragging all */
  BlockRect DragRect;
  int DragX, DragY, DragMenu;
  int DragStartFrame, DragCurFrame, DragNearest, DragDraw;
  int DragColumn;
  int LabelIndent;
  int PanelActive;
};

void MovieViewReinterpolate(PyMOLGlobals *G)
{
  float power  = SettingGetGlobal_f(G, cSetting_motion_power);
  float bias   = SettingGetGlobal_f(G, cSetting_motion_bias);
  float linear = SettingGetGlobal_f(G, cSetting_motion_linear);
  int hand     = SettingGetGlobal_i(G, cSetting_motion_hand);

  MovieView(G, 3, -1, -1, power, bias, 1, /* note simple always = 1 for camera motion...*/
            linear, 
            SettingGetGlobal_b(G,cSetting_movie_loop) ? 1 : 0 ,
            hand, 5, 1, NULL, 0.5, -1, 1); 
}

int MovieXtoFrame(PyMOLGlobals *G, BlockRect *rect, int frames, int x, int nearest)
{
  register CMovie *I = G->Movie;
  return ViewElemXtoFrame(G,I->ViewElem,rect,frames,x,nearest);
}

void MovieViewTrim(PyMOLGlobals *G,int n_frame)
{
  register CMovie *I = G->Movie;
  if(n_frame>=0) {
    if(!I->Sequence) {
      I->Sequence = VLACalloc(int, n_frame);
    } else {
      VLASize(I->Sequence, int, n_frame);
    }
    if(!I->Cmd) {
      I->Cmd = VLACalloc(MovieCmdType, n_frame);
    } else {
      VLASize(I->Cmd, MovieCmdType, n_frame);
    }
    if(!I->ViewElem) {
      I->ViewElem = VLACalloc(CViewElem, n_frame);
    } else {
      VLASize(I->ViewElem, CViewElem, n_frame);
    }
    I->NFrame = n_frame;
  }
}


int MovieViewModify(PyMOLGlobals *G,int action, int index, int count,int target, int freeze, int localize)
{
  register CMovie *I = G->Movie;
  int ok = true;
  MovieClearImages(G);
  if( (ok = ViewElemModify(G,&I->ViewElem, action, index, count, target)) ) {
    switch(action) {
    case cViewElemModifyInsert: 
      VLAInsert(I->Sequence,int,index,count);
      VLAInsert(I->Cmd,MovieCmdType,index,count);
      I->NFrame = VLAGetSize(I->Sequence);
      {
        int frame = SceneGetFrame(G);
        if(frame >= index) {
          SceneSetFrame(G,0,frame + count);
        }
      }
      break;
    case cViewElemModifyDelete:
      VLADelete(I->Sequence,int,index,count);
      VLADelete(I->Cmd,MovieCmdType,index,count);
      I->NFrame = VLAGetSize(I->Sequence);
      break;
    case cViewElemModifyMove:
      if((index>=0) && (target>=0) && (index<I->NFrame) && (target<I->NFrame)) {
        int i;
        for(i=0;i<count;i++) {
          if( ((i+index)<I->NFrame) && ((i+target)<I->NFrame)) {
            int src,dst;
            if(index>target) {
              src = index+i;
              dst = target+i;
            } else {
              src = index+(count-1)-i;
              dst = target+(count-1)-i;
            }
            I->Sequence[dst] = I->Sequence[src];
            memcpy(I->Cmd + dst, I->Cmd + src, sizeof(MovieCmdType));
            I->Cmd[src][0]=0;
          }
        }
      }
      break;
    case cViewElemModifyCopy:
      if((index>=0) && (target>=0) && (index<I->NFrame) && (target<I->NFrame)) {
        int i;
        for(i=0;i<count;i++) {
          if( ((i+index)<I->NFrame) && ((i+target)<I->NFrame)) {
            int src,dst;
            if(index>target) {
              src = index+i;
              dst = target+i;
            } else {
              src = index+(count-1)-i;
              dst = target+(count-1)-i;
            }
            memcpy(I->Cmd + dst, I->Cmd + src, sizeof(MovieCmdType));
          }
        }
      }
      break;
    }
  }
  if(ok && ((!freeze)&&(!localize))) {
    ExecutiveMotionExtend(G,freeze);
  }
  return ok;
}

int MovieGetSpecLevel(PyMOLGlobals *G,int frame)
{
  register CMovie *I = G->Movie;
  if(I->ViewElem) {
    int size = VLAGetSize(I->ViewElem);
    if(frame<0) {
      int max_level = 0;
      int i;
      for(i=0;i<size;i++) {
        if(max_level < I->ViewElem[i].specification_level)
          max_level = I->ViewElem[i].specification_level;
      }
      return max_level;
    }
    if((frame>=0) && (frame<size))
      return I->ViewElem[frame].specification_level;
    return 0;
  }
  return -1;
}

void MovieSetRealtime(PyMOLGlobals * G, int realtime)
{
  register CMovie *I = G->Movie;
  I->RealtimeFlag = realtime;
}

int MovieGetRealtime(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  return I->RealtimeFlag;
}

void MovieCopyPrepare(PyMOLGlobals * G, int *width, int *height, int *length)
{
  /* assumed locked api, blocked threads, and master thread on entry */

  /* writes the current movie sequence to a set of PNG files 
   * this routine can take a LOT of time -- 
   * TODO: develop an interrupt mechanism */

  int start, stop;
  register CMovie *I = G->Movie;
  int nFrame;

  I->CacheSave = (int) SettingGet(G, cSetting_cache_frames);
  I->OverlaySave = (int) SettingGet(G, cSetting_overlay);
  if(!I->CacheSave)
    MovieClearImages(G);
  SettingSet(G, cSetting_cache_frames, 1.0);
  SettingSet(G, cSetting_overlay, 5);
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame = SceneGetNFrame(G, NULL);
  }
  start = 0;
  stop = nFrame;
  if((start != 0) || (stop != (nFrame + 1)))
    SceneSetFrame(G, 0, 0);
  MoviePlay(G, cMoviePlay);
  VLACheck(I->Image, ImageType *, nFrame);
  SceneGetWidthHeight(G, width, height);
  {
    int uniform_height = -1;
    int uniform_width = -1;
    int uniform_flag = false;
    int scene_match = true;
    int a;
    ImageType *image;
    /* make sure all the movie frames match the screen size or are pre-rendered and are already the same size */
    for(a = 0; a < nFrame; a++) {
      image = I->Image[a];
      if(image) {
        if((image->height != *height) || (image->width != *width)) {
          scene_match = false;
          if(uniform_height < 0) {
            uniform_height = image->height;
            uniform_width = image->width;
          } else {
            if((image->height != uniform_height) || (image->width != uniform_width))
              uniform_flag = false;
          }
        }
      } else
        uniform_flag = false;   /* missing at least one image, so not uniform */
    }
    if(!scene_match) {
      if(uniform_flag) {
        *height = uniform_height;
        *width = uniform_width;
      } else {
        MovieClearImages(G);
      }
    }

  }
  *length = nFrame;
}

void MovieFlushCommands(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  I->RecursionFlag = true;
  PFlush(G);
  I->RecursionFlag = false;
}

int MovieCopyFrame(PyMOLGlobals * G, int frame, int width, int height, int rowbytes,
                   void *ptr)
{
  register CMovie *I = G->Movie;
  int result = false;
  int nFrame;
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame = SceneGetNFrame(G, NULL);
  }

  if((frame < nFrame) && (ptr)) {
    int a = frame;
    int i;
    SceneSetFrame(G, 0, a);
    MovieDoFrameCommand(G, a);
    MovieFlushCommands(G);
    i = MovieFrameToImage(G, a);
    VLACheck(I->Image, ImageType *, i);
    if(!I->Image[i]) {
      SceneUpdate(G, false);
      SceneMakeMovieImage(G, false, false, cSceneImage_Default);
    }
    if(!I->Image[i]) {
      PRINTFB(G, FB_Movie, FB_Errors)
        "MoviePNG-Error: Missing rendered image.\n" ENDFB(G);
    } else {
      if((I->Image[i]->height == height) && (I->Image[i]->width == width)) {
        unsigned char *srcImage = (unsigned char *) I->Image[i]->data;
        int i, j;
        for(i = 0; i < height; i++) {
          unsigned char *dst = ((unsigned char *) ptr) + i * rowbytes;
          unsigned char *src = srcImage + ((height - 1) - i) * width * 4;
          for(j = 0; j < width; j++) {
            *dst++ = src[3];
            *dst++ = src[0];
            *dst++ = src[1];
            *dst++ = src[2];
            src += 4;
          }
        }
        result = true;
      } else {
        /* mismatched dimensions, so furnish a white image */
        memset(ptr, 0xFF, 4 * height * width);
      }
      ExecutiveDrawNow(G);
      if(G->HaveGUI)
        PyMOL_SwapBuffers(G->PyMOL);
    }
    if(!I->CacheSave) {
      if(I->Image[i]) {
        FreeP(I->Image[i]->data);
        FreeP(I->Image[i]);
      }
    }
  }
  return result;
}

int MoviePurgeFrame(PyMOLGlobals * G, int frame)
{
  register CMovie *I = G->Movie;
  int result = false;
  int nFrame;
  int i;
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame = SceneGetNFrame(G, NULL);
  }
  if(!I->CacheSave) {
    if(frame < nFrame) {
      int a = frame;
      i = MovieFrameToImage(G, a);
      VLACheck(I->Image, ImageType *, i);
      if(I->Image[i]) {
        FreeP(I->Image[i]->data);
        FreeP(I->Image[i]);
        I->Image[i] = NULL;
        result = true;
      }
    }
  }
  return result;
}

void MovieCopyFinish(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  SceneInvalidate(G);           /* important */
  SettingSet(G, cSetting_cache_frames, (float) I->CacheSave);
  SettingSet(G, cSetting_overlay, (float) I->OverlaySave);
  MoviePlay(G, cMovieStop);
  if(!I->CacheSave) {
    MovieClearImages(G);
  }
}

int MovieLocked(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  return (I->Locked);
}

void MovieSetLock(PyMOLGlobals * G, int lock)
{
  register CMovie *I = G->Movie;
  I->Locked = lock;
}


/*========================================================================*/
void MovieDump(PyMOLGlobals * G)
{
  int a;
  register CMovie *I = G->Movie;
  int flag = false;
  char buffer[OrthoLineLength + 100];

  for(a = 0; a < I->NFrame; a++) {
    if(I->Cmd[a][0]) {
      flag = true;
      break;
    }
  }
  if(flag && I->NFrame) {
    PRINTFB(G, FB_Movie, FB_Results)
      " Movie: General Purpose Commands:\n" ENDFB(G);
    for(a = 0; a < I->NFrame; a++) {
      if(I->Cmd[a][0]) {
        sprintf(buffer, "%5d: %s\n", a + 1, I->Cmd[a]);
        OrthoAddOutput(G, buffer);
      }
    }
  } else {
    PRINTFB(G, FB_Movie, FB_Results)
      " Movie: No movie commands are defined.\n" ENDFB(G);
  }
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
static int MovieCmdFromPyList(PyMOLGlobals * G, PyObject * list, int *warning)
{

  register CMovie *I = G->Movie;
  int ok = true;
  int a;
  int warn = false;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);

  for(a = 0; a < I->NFrame; a++) {
    if(ok)
      ok = PConvPyStrToStr(PyList_GetItem(list, a), I->Cmd[a], OrthoLineLength);
    if(ok)
      warn = (warn || I->Cmd[a][0]);
  }
  *warning = warn;
  return (ok);

}
#endif

/*========================================================================*/
int MovieFromPyList(PyMOLGlobals * G, PyObject * list, int *warning)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  register CMovie *I = G->Movie;
  int ll = 0;

  MovieReset(G);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  /* TO SUPPORT BACKWARDS COMPATIBILITY...
     Always check ll when adding new PyList_GetItem's */
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 0), &I->NFrame);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 1), &I->MatrixFlag);
  if(ok && I->MatrixFlag)
    ok =
      PConvPyListToFloatArrayInPlace(PyList_GetItem(list, 2), I->Matrix, cSceneViewSize);
  if(ok)
    ok = PConvPyIntToInt(PyList_GetItem(list, 3), &I->Playing);
  if(ok && I->NFrame) {
    I->Sequence = VLACalloc(int, I->NFrame);
    I->Cmd = VLACalloc(MovieCmdType, I->NFrame);
    if(ok)
      ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 4), I->Sequence, I->NFrame);
    if(ok)
      ok = MovieCmdFromPyList(G, PyList_GetItem(list, 5), warning);
    if((*warning) && G->Security) {
      MovieSetLock(G, true);
    }
  }
  if(ok && (ll > 6)) {
    PyObject *tmp;
    VLAFreeP(I->ViewElem);
    I->ViewElem = NULL;
    tmp = PyList_GetItem(list, 6);
    if(tmp && !(tmp == Py_None))
      ok = ViewElemVLAFromPyList(G, tmp, &I->ViewElem, I->NFrame);
  }
  if(!ok) {
    MovieReset(G);
  } else if(MovieDefined(G)) {
    OrthoReshape(G,-1,-1,true);
  }
  return (ok);
#endif
}


/*========================================================================*/
#ifndef _PYMOL_NOPY
static PyObject *MovieCmdAsPyList(PyMOLGlobals * G)
{

  register CMovie *I = G->Movie;
  PyObject *result = NULL;
  int a;

  result = PyList_New(I->NFrame);
  if(result)
    for(a = 0; a < I->NFrame; a++) {
      PyList_SetItem(result, a, PyString_FromString(I->Cmd[a]));
    }
  return (PConvAutoNone(result));

}
#endif

/*========================================================================*/
PyObject *MovieAsPyList(PyMOLGlobals * G)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else

  register CMovie *I = G->Movie;
  PyObject *result = NULL;

  result = PyList_New(7);
  PyList_SetItem(result, 0, PyInt_FromLong(I->NFrame));
  PyList_SetItem(result, 1, PyInt_FromLong(I->MatrixFlag));
  PyList_SetItem(result, 2, PConvFloatArrayToPyList(I->Matrix, cSceneViewSize));
  PyList_SetItem(result, 3, PyInt_FromLong(I->Playing));
  if(I->Sequence) {
    PyList_SetItem(result, 4, PConvIntArrayToPyList(I->Sequence, I->NFrame));
  } else {
    PyList_SetItem(result, 4, PConvAutoNone(NULL));
  }
  if(I->Cmd) {
    PyList_SetItem(result, 5, MovieCmdAsPyList(G));
  } else {
    PyList_SetItem(result, 5, PConvAutoNone(NULL));
  }
  if(I->ViewElem) {
    PyList_SetItem(result, 6, ViewElemVLAAsPyList(G, I->ViewElem, I->NFrame));
  } else {
    PyList_SetItem(result, 6, PConvAutoNone(NULL));
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
  return (PConvAutoNone(result));
#endif
}


/*========================================================================*/
int MoviePlaying(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  if(I->Locked)
    return false;
  if(I->Playing && G->Interrupt) {
    I->Playing = false;
  }
  return (I->Playing || I->RecursionFlag);
  /* returns true if movie is playing OR if we're in the process of evaluating
     movie commands */
}


/*========================================================================*/
void MoviePlay(PyMOLGlobals * G, int cmd)
{
  register CMovie *I = G->Movie;
  switch (cmd) {
  case cMovieToggle:
    I->Playing = !I->Playing;
    if(I->Playing && !(int) SettingGet(G, cSetting_movie_loop)) {
      /* if not looping, and at end of movie, then automatically rewind
         and force execution of the first movie command */
      if((SettingGetGlobal_i(G, cSetting_frame)) == (SceneGetNFrame(G, NULL))) {
        SceneSetFrame(G, 7, 0);
      }
    }
    break;
  case cMovieStop:
    I->Playing = false;
    break;
  case cMoviePlay:
    if(!(int) SettingGet(G, cSetting_movie_loop)) {
      /* if not looping, and at end of movie, then automatically rewind
         and force execution of the first movie command */
      if((SettingGetGlobal_i(G, cSetting_frame)) == (SceneGetNFrame(G, NULL))) {
        SceneSetFrame(G, 7, 0);
      }
    }
    I->Playing = true;
    break;
  }
  OrthoDirty(G);
  SceneRestartFrameTimer(G);
}


/*========================================================================*/
int MovieMatrix(PyMOLGlobals * G, int action)
{
  register CMovie *I = G->Movie;
  int result = false;
  switch (action) {
  case cMovieMatrixClear:
    I->MatrixFlag = false;
    result = 1;
    break;
  case cMovieMatrixStore:
    SceneGetView(G, I->Matrix);
    I->MatrixFlag = true;
    result = 1;
    break;
  case cMovieMatrixRecall:
    if(I->MatrixFlag) {
      SceneSetView(G, I->Matrix, true, 0, 0);
      result = 1;
    } else
      result = 0;
    break;
  case cMovieMatrixCheck:
    result = I->MatrixFlag;
    break;
  }
  return result;
}


/*========================================================================*/
void MovieSetSize(PyMOLGlobals * G, unsigned int width, unsigned int height)
{
  return;
}


/*========================================================================*/
static void MovieModalPNG(PyMOLGlobals * G, CMovie * I, CMovieModal * M)
{
  switch (M->stage) {
  case 0:                      /* setup */
    MovieSetRealtime(G, false);
    M->save = (int) SettingGet(G, cSetting_cache_frames);
    if(!M->save)
      MovieClearImages(G);
    SettingSet(G, cSetting_cache_frames, 1.0);
    OrthoBusyPrime(G);
    M->nFrame = I->NFrame;
    if(!M->nFrame) {
      M->nFrame = SceneGetNFrame(G, NULL);
      if(M->nFrame < 1) {
        M->nFrame = 1;
      }
    }
    if(M->start < 0)
      M->start = 0;
    if(M->start > M->nFrame)
      M->start = M->nFrame;
    if(M->stop < 0)
      M->stop = M->nFrame;
    if(M->stop > M->nFrame)
      M->stop = M->nFrame;
    {
      OrthoLineType buffer;
      sprintf(buffer, "Creating movie (%d frames)...", M->nFrame);
      OrthoBusyMessage(G, buffer);
    }
    if((M->start != 0) || (M->stop != (M->nFrame + 1)))
      SceneSetFrame(G, 0, 0);
    MoviePlay(G, cMoviePlay);
    VLACheck(I->Image, ImageType *, M->nFrame);
    M->frame = 0;
    M->stage = 1;
    if(G->Interrupt) {
      M->stage = 5;             /* abort */
    }
    break;
  case 1:                      /* RENDER LOOP: advance a frame */
    if(M->frame < M->nFrame) {

      M->file_missing = true;
      M->timing = UtilGetSeconds(G);    /* start timing the process */

      PRINTFB(G, FB_Movie, FB_Debugging)
        " MoviePNG-DEBUG: Cycle %d...\n", M->frame ENDFB(G);
      switch (M->format) {
      case cMyPNG_FormatPPM:
        sprintf(M->fname, "%s%04d.ppm", M->prefix, M->frame + 1);
        break;
      case cMyPNG_FormatPNG:
      default:
        sprintf(M->fname, "%s%04d.png", M->prefix, M->frame + 1);
        break;
      }

      if(M->missing_only) {
        FILE *tmp = fopen(M->fname, "rb");
        if(tmp) {
          fclose(tmp);
          M->file_missing = false;
        } else {
          M->file_missing = true;
        }
      }
      SceneSetFrame(G, 0, M->frame);
      MovieDoFrameCommand(G, M->frame);
      MovieFlushCommands(G);

      M->image = MovieFrameToImage(G, M->frame);
      M->stage = 2;
      if(G->Interrupt) {
        M->stage = 5;           /* abort */
      }
    }
    break;
  }

  switch (M->stage) {
  case 2:                      /* IN RENDER LOOP: create the image */
    VLACheck(I->Image, ImageType *, M->image);
    if((M->frame >= M->start) &&        /* only render frames in the specified interval... */
       (M->frame <= M->stop) && (M->file_missing)) {    /* ...that don't already exist */
      if(!I->Image[M->image]) {
        SceneUpdate(G, false);
        if(SceneMakeMovieImage(G, false, M->modal, M->mode) || (!M->modal)) {
          M->stage = 3;
        } else {
          /* didn't get an image... */
          PRINTFB(G, FB_Movie, FB_Errors)
            " MoviePNG-Error: unable to obtain a valid OpenGL image.  Trying again...\n"
            ENDFB(G);
        }
      } else {                  /* frame already rendered */
        M->stage = 3;
      }
    } else {                    /* we don't need to render this frame */
      M->stage = 4;
    }
    if(G->Interrupt) {
      M->stage = 5;             /* abort */
    }
    break;
  }

  switch (M->stage) {
  case 3:                      /* IN RENDER LOOP: have image, so write to file */
    if(!I->Image[M->image]) {
      PRINTFB(G, FB_Movie, FB_Errors)
        "MoviePNG-Error: Missing rendered image.\n" ENDFB(G);
    } else {
      if(!MyPNGWrite(G, M->fname, I->Image[M->image]->data,
                     I->Image[M->image]->width,
                     I->Image[M->image]->height,
                     SettingGetGlobal_f(G, cSetting_image_dots_per_inch),
                     M->format, M->quiet)) {
        PRINTFB(G, FB_Movie, FB_Errors)
          " MoviePNG-Error: unable to write '%s'\n", M->fname ENDFB(G);
      }
      ExecutiveDrawNow(G);
      OrthoBusySlow(G, M->frame, M->nFrame);
      if(G->HaveGUI)
        PyMOL_SwapBuffers(G->PyMOL);
      PRINTFB(G, FB_Movie, FB_Debugging)
        " MoviePNG-DEBUG: i = %d, I->Image[image] = %p\n", M->image,
        I->Image[M->image]->data ENDFB(G);
    }
    if(I->Image[M->image]) {
      FreeP(I->Image[M->image]->data);
      FreeP(I->Image[M->image]);
    }
    M->timing = UtilGetSeconds(G) - M->timing;
    M->accumTiming += M->timing;
    {
      double est1 = (M->nFrame - M->frame) * M->timing;
      double est2 = ((M->nFrame - M->frame) / (float) (M->frame + 1)) * M->accumTiming;

      PRINTFB(G, FB_Movie, FB_Details)
        " Movie: frame %4d of %4d, %4.2f sec. (%d:%02d:%02d - %d:%02d:%02d to go).\n",
        M->frame + 1, M->nFrame,
        M->timing,
        (int) (est1 / 3600),
        ((int) (est1 / 60)) % 60,
        ((int) est1) % 60,
        (int) (est2 / 3600), ((int) (est2 / 60)) % 60, ((int) est2) % 60 ENDFB(G);
    }
    M->stage = 4;
    if(G->Interrupt) {
      M->stage = 5;             /* abort */
    }
    break;
  }

  switch (M->stage) {
  case 4:                      /* IN RENDER LOOP: advance to next frame or kick out when done */
    M->frame++;
    if(M->frame >= M->nFrame) {
      M->stage = 5;
    } else {
      M->stage = 1;             /* RESTART LOOP */
    }
    if(G->Interrupt) {
      M->stage = 5;             /* abort */
    }
    break;
  }

  switch (M->stage) {
  case 5:                      /* finish up */

    SceneInvalidate(G);         /* important */
    PRINTFB(G, FB_Movie, FB_Debugging)
      " MoviePNG-DEBUG: done.\n" ENDFB(G);
    SettingSet(G, cSetting_cache_frames, (float) M->save);
    MoviePlay(G, cMovieStop);
    MovieClearImages(G);
    MovieSetRealtime(G, true);
    M->complete = true;
    M->stage = 6;
    break;
  }
}

static void MovieModalDraw(PyMOLGlobals * G);

static void MovieModalDraw(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  MovieModalPNG(G, I, &I->Modal);
  if(!I->Modal.complete)        /* force modalic return until job is done */
    PyMOL_SetModalDraw(G->PyMOL, (PyMOLModalDrawFn *) MovieModalDraw);
}

int MoviePNG(PyMOLGlobals * G, char *prefix, int save, int start,
             int stop, int missing_only, int modal, int format, int mode, int quiet)
{
  /* assumes locked api, blocked threads, and master thread on entry */
  register CMovie *I = G->Movie;

  /* new routine allows PyMOL to workaround XP and Vista Windowing
     issues */

  CMovieModal *M = &I->Modal;

  UtilZeroMem(M, sizeof(CMovieModal));

  UtilNCopy(M->prefix, prefix, sizeof(OrthoLineType));
  M->save = save;
  M->start = start;
  M->stop = stop;
  M->missing_only = missing_only;
  M->stage = 0;
  M->format = format;
  M->mode = mode;
  M->quiet = quiet;

  if(modal < 0) {
    /* default behavior is to go modal unless we're ray tracing */
    if((mode < cSceneImage_Ray) || (!SettingGetGlobal_b(G, cSetting_ray_trace_frames))) {
      modal = 1;
    }
  }
  M->modal = modal;

  if(modal) {
    PyMOL_SetModalDraw(G->PyMOL, (PyMOLModalDrawFn *) MovieModalDraw);
  } else {
    while(!M->complete) {
      MovieModalPNG(G, I, &I->Modal);
    }
  }
  return true;
}


/*========================================================================*/
void MovieAppendSequence(PyMOLGlobals * G, char *str, int start_from,int freeze)
{
  register CMovie *I = G->Movie;
  int c = 0;
  int i;
  /*  int old_NFrame = I->NFrame; */
  char *s, number[20];

  if(start_from < 0)
    start_from = I->NFrame;

  c = start_from;

  PRINTFB(G, FB_Movie, FB_Debugging)
    " MovieSequence: entered. str:%s\n", str ENDFB(G);

  s = str;
  while(*s) {
    s = ParseWord(number, s, 20);
    if(sscanf(number, "%i", &i)) {      /* slow */
      c++;
    }
  }

  if(!c) {
    VLAFreeP(I->Sequence);
    VLAFreeP(I->Cmd);
    VLAFreeP(I->ViewElem);
    I->NFrame = 0;
  } else {
    if(!I->Sequence) {
      I->Sequence = VLACalloc(int, c);
    } else {
      VLASize(I->Sequence, int, start_from);    /* to clear */
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

  if(c && str[0]) {             /* not just a reset */
    for(i = start_from; i < c; i++)
      I->Cmd[i][0] = 0;
    c = start_from;
    s = str;
    while(*s) {
      s = ParseWord(number, s, 20);
      if(sscanf(number, "%i", &I->Sequence[c])) {
        c++;
      }
    }
    I->NFrame = c;
  } else if(!str[0]) {
    I->NFrame = start_from;
  }

  VLASize(I->Image, ImageType *, I->NFrame);
  PRINTFB(G, FB_Movie, FB_Debugging)
    " MovieSequence: leaving... I->NFrame%d\n", I->NFrame ENDFB(G);

  if(!freeze) {
    if(SettingGetGlobal_b(G,cSetting_movie_auto_interpolate)) 
      ExecutiveMotionReinterpolate(G);
  }
  ExecutiveCountMotions(G);
}


/*========================================================================*/
int MovieFrameToImage(PyMOLGlobals * G, int frame)
{
  int result = 0;
  int single_image = (int) SettingGet(G, cSetting_single_image);
  if(single_image)
    result = MovieFrameToIndex(G, frame);
  else
    result = frame;
  PRINTFB(G, FB_Movie, FB_Debugging)
    " MovieFrameToImage-DEBUG: result %d\n", result ENDFB(G);
  return (result);
}


/*========================================================================*/
int MovieFrameToIndex(PyMOLGlobals * G, int frame)
{
  register CMovie *I = G->Movie;
  if(I->Sequence && I->NFrame) {
    if(frame >= I->NFrame) {
      frame = I->NFrame - 1;
    }
    if(I->ViewElem && I->ViewElem[frame].state_flag) {
      return I->ViewElem[frame].state;
    }
    return (I->Sequence[frame]);
  } else {
    return (frame);
  }
}
/*========================================================================*/
void MovieSetImage(PyMOLGlobals * G, int index, ImageType * image)
{
  register CMovie *I = G->Movie;

  PRINTFB(G, FB_Movie, FB_Blather)
    " MovieSetImage: setting movie image %d\n", index + 1 ENDFB(G);

  VLACheck(I->Image, ImageType *, index);
  if(I->Image[index])
    FreeP(I->Image[index]);
  I->Image[index] = image;
  if(I->NImage < (index + 1))
    I->NImage = index + 1;
}

int MovieSeekScene(PyMOLGlobals * G, int loop)
{
  register CMovie *I = G->Movie;
  int result = -1;
  OVreturn_word ret;
  char *scene_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
  if(OVreturn_IS_OK
     ((ret = OVLexicon_BorrowFromCString
       (G->Lexicon, scene_name)))) {
    if(I->ViewElem) {
      int i,len = MovieGetLength(G);
      for(i = SceneGetFrame(G); i < len; i++) {
	if(I->ViewElem[i].scene_flag) {
	  if(I->ViewElem[i].scene_name == ret.word) {
	    result = i;
	    break;
	  }
	}
      }
      if(loop) {
	len = SceneGetFrame(G);
	for(i = 0; i < len; i++ ) {
	  if(I->ViewElem[i].scene_flag) {
	    if(I->ViewElem[i].scene_name == ret.word) {
	      result = i;
	      break;
	    }
	  }
	}
      }
    }
  }

  return result;
}

/*========================================================================*/
void MovieDoFrameCommand(PyMOLGlobals * G, int frame)
{
  register CMovie *I = G->Movie;
  if(frame == 0)
    MovieMatrix(G, cMovieMatrixRecall);
  if(!I->Locked) {
    if((frame >= 0) && (frame < I->NFrame)) {
      if(I->Cmd[frame][0]) {
        if(!I->RecursionFlag) {
          PParse(G, I->Cmd[frame]);
        }
      }
      if(I->ViewElem) {
        if(I->ViewElem[frame].scene_flag) {
          char *st = OVLexicon_FetchCString(G->Lexicon, I->ViewElem[frame].scene_name);
          if(strcmp(st, SettingGetGlobal_s(G, cSetting_scene_current_name))) {
#ifndef _PYMOL_NOPY
            PBlock(G);
            PXDecRef(PyObject_CallMethod(G->P_inst->cmd, "scene",
                                         "sssiiiii", st, "recall", NULL, 0, 1, 1, 1, 0, 0));
            if(PyErr_Occurred()) {
              PyErr_Clear();
            }

            PUnblock(G);
#endif
          }
        }
        SceneFromViewElem(G, I->ViewElem + frame, true);
      }
    }
  }
}


/*========================================================================*/
void MovieSetCommand(PyMOLGlobals * G, int frame, char *command)
{
  register CMovie *I = G->Movie;
  int a, len;
  if((frame >= 0) && (frame < I->NFrame)) {
    len = strlen(command);
    if(len > (sizeof(MovieCmdType) - 1))
      len = sizeof(MovieCmdType) - 1;
    for(a = 0; a < len; a++)
      I->Cmd[frame][a] = command[a];
    I->Cmd[frame][len] = 0;
  } else {
    PRINTFB(G, FB_Movie, FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",
      frame + 1 ENDFB(G);
  }
}


/*========================================================================*/
int MovieView(PyMOLGlobals * G, int action, int first,
              int last, float power, float bias,
              int simple, float linear, int wrap,
              int hand, int window, int cycles,
              char *scene_name, float scene_cut, int state, int quiet)
{
  register CMovie *I = G->Movie;
  int frame;
  
  if(wrap<0) {
    wrap = SettingGetGlobal_b(G,cSetting_movie_loop);
  }
  if((action == 7) || (action == 8)) { /* toggle */
    frame = first;
    if(first < 0)
      frame = SceneGetFrame(G);
    VLACheck(I->ViewElem, CViewElem, frame);
    if(action == 7) {
      if(I->ViewElem[frame].specification_level>1) {
        action = 1;
      } else {
        action = 0;
      }
    } else if(action == 8) {
      if(I->ViewElem[frame].specification_level>1) {
        int frame;
        action = 3;
        for(frame=0;frame<I->NFrame;frame++) {
          if(I->ViewElem[frame].specification_level==1) {
            action = 6;
            break;
          }
        }
      }
      else if(I->ViewElem[frame].specification_level>0) {
        action = 6;
      } else {
        action = 3;
      }
    }
  }
  if(action == 4) {
    if(I->ViewElem) {
      int save_last = last;
      if(first < 0)
        first = 0;
      
      if(last < 0) {
        last = SceneGetNFrame(G, NULL) - 1;
      }
      if(last >= I->NFrame) {
        last = I->NFrame - 1;
      }
      if(first <= last) {
        int a;
        VLACheck(I->ViewElem, CViewElem, last);
        for(a = 0; a < cycles; a++) {
          ViewElemSmooth(I->ViewElem + first, I->ViewElem + last, window, wrap);
        }
      }
      if(SettingGetGlobal_b(G, cSetting_movie_auto_interpolate)) {
        action = 3; /* reinterpolate */
        last = save_last;
      }
    }
  }
  switch (action) {
  case 0:                      /* store */
    if(I->ViewElem) {
      if(first < 0)
        first = SceneGetFrame(G);
      if(last < 0)
        last = first;

      VLACheck(I->ViewElem, CViewElem, last);

      for(frame = first; frame <= last; frame++) {
        if((frame >= 0) && (frame < I->NFrame)) {
          VLACheck(I->ViewElem, CViewElem, frame);
          if(!quiet) {
            PRINTFB(G, FB_Movie, FB_Details)
              " MovieView: Setting frame %d.\n", frame + 1 ENDFB(G);
          }
          if(scene_name && (!scene_name[0]))
            scene_name = NULL;
          SceneToViewElem(G, I->ViewElem + frame, scene_name);
          if(state>=0) {
            I->ViewElem[frame].state = state;
            I->ViewElem[frame].state_flag = true;
          } else if(state==-2) {
            I->ViewElem[frame].state = SceneGetState(G);
            I->ViewElem[frame].state_flag = true;
          }

          if(power!=0.0F) {
            I->ViewElem[frame].power_flag = true;
            I->ViewElem[frame].power = power;
          }

          if(bias > 0.0F) {
            I->ViewElem[frame].bias_flag = true;
            I->ViewElem[frame].bias = bias;
          }

          I->ViewElem[frame].specification_level = 2;
        }
      }
    }
    break;
  case 1:                      /* clear */
    if(I->ViewElem) {
      if(first < 0)
        first = SceneGetFrame(G);
      if(last < 0)
        last = first;
      for(frame = first; frame <= last; frame++) {
        if((frame >= 0) && (frame < I->NFrame)) {
          VLACheck(I->ViewElem, CViewElem, frame);
          ViewElemArrayPurge(G, I->ViewElem + frame, 1);
          UtilZeroMem((void *) (I->ViewElem + frame), sizeof(CViewElem));
        }
      }
    }
    break;
  case 2:                      /* interpolate & reinterpolate */
  case 3:
    if(I->ViewElem) {
      int view_found = false;
      CViewElem *first_view = NULL, *last_view = NULL;
      if(first < 0)
        first = 0;

      if(first > I->NFrame) {
        first = I->NFrame - 1;
      }
      /* note that we're leaving a blank frame at the end... */
      if(last < 0) {
        last = MovieGetLength(G);

        if(last) {
          if(!wrap)
            last--;
          else {
            int frame = 0;
            VLACheck(I->ViewElem, CViewElem, last);
            for(frame = 0; frame < last; frame++) {
              if(I->ViewElem[frame].specification_level > 1) {
                last += frame;
                break;
              }
            }
          }
        }
      } else if(last >= I->NFrame) {
        last = I->NFrame;
        if(last && (!wrap))
          last--;
      }
      
      VLACheck(I->ViewElem, CViewElem, last);

      if(wrap && (last >= I->NFrame)) {
        /* if we're interpolating beyond the last frame, then wrap by
           copying early frames to last frames */
        int a;
        for(a = I->NFrame; a <= last; a++) {
          ViewElemCopy(G, I->ViewElem + a - I->NFrame, I->ViewElem + a);
        }
      } else if(!wrap) { 
        /* if we're not wrapping, then make sure we nuke any stray / old
           interpolated frames */
        frame = I->NFrame - 1;
        while(frame>=0) {
          if(I->ViewElem[frame].specification_level > 1) 
            break;
          else
            UtilZeroMem((void *) (I->ViewElem + frame), sizeof(CViewElem));
          frame--;
        }
      }

      if(!quiet) {
        if(action == 2) {
          if(last == I->NFrame) {
            PRINTFB(G, FB_Movie, FB_Details)
              " MovieView: interpolating unspecified frames %d to %d (wrapping)\n",
              first + 1, last ENDFB(G);
          } else {
            PRINTFB(G, FB_Movie, FB_Details)
              " MovieView: interpolating unspecified frames %d to %d.\n", first + 1,
              last + 1 ENDFB(G);
          }

        } else {
          if(last == I->NFrame) {
            PRINTFB(G, FB_Movie, FB_Details)
              " MovieView: reinterpolating all frames %d to %d (wrapping).\n", first + 1,
              last ENDFB(G);
          } else {
            PRINTFB(G, FB_Movie, FB_Details)
              " MovieView: reinterpolating all frames %d to %d.\n", first + 1, last + 1
              ENDFB(G);
          }
        }
      }

      for(frame = first; frame <= last; frame++) {
        if(!first_view) {
          if(I->ViewElem[frame].specification_level == 2) {     /* specified */
            first_view = I->ViewElem + frame;
            view_found = true;
          }
        } else {
          CViewElem *view;
          int interpolate_flag = false;
          if(I->ViewElem[frame].specification_level == 2) {     /* specified */
            last_view = I->ViewElem + frame;
            if(action == 2) {   /* interpolate */
              for(view = first_view + 1; view < last_view; view++) {
                if(!view->specification_level)
                  interpolate_flag = true;
              }
            } else {
              interpolate_flag = true;
            }
            if(interpolate_flag) {
              ViewElemInterpolate(G, first_view, last_view,
                                  power, bias, simple, linear, hand, scene_cut);
            }
            first_view = last_view;
            last_view = NULL;
          }
        }
      }

      if(first_view) {
        if(wrap && (last >= I->NFrame)) {
          /* if we're interpolating beyond the last frame, then wrap by
             copying the last frames back over the early frames */
          int a;
          for(a = I->NFrame; a <= last; a++) {
            ViewElemCopy(G, I->ViewElem + a, I->ViewElem + a - I->NFrame);
          }
        }
      }

      if((!view_found) && (last>=first) && (first>=0) && (last<=I->NFrame)) {
        UtilZeroMem(I->ViewElem + first, sizeof(CViewElem) * (1 + (last-first)));
      }

      if(last >= I->NFrame) {   /* now erase temporary views */
        ViewElemArrayPurge(G, I->ViewElem + I->NFrame, (1 + last - I->NFrame));
        UtilZeroMem((void *) (I->ViewElem + I->NFrame),
                    sizeof(CViewElem) * (1 + last - I->NFrame));
      }
    }
    break;
  case 5:                      /* reset */
    if(I->ViewElem) {
      int size = VLAGetSize(I->ViewElem);
      VLAFreeP(I->ViewElem);
      I->ViewElem = VLACalloc(CViewElem, size);
    }
    break;
  case 6:                      /* uninterpolate */
    if(I->ViewElem) {
      if(first < 0)
        first = 0;
      if(last < 0) {
        last = SceneGetNFrame(G, NULL) - 1;
      }
      for(frame = first; frame <= last; frame++) {
        if((frame >= 0) && (frame < I->NFrame)) {
          VLACheck(I->ViewElem, CViewElem, frame);
          if(I->ViewElem[frame].specification_level < 2) {
            ViewElemArrayPurge(G, I->ViewElem + frame, 1);
            UtilZeroMem((void *) (I->ViewElem + frame), sizeof(CViewElem));
          }
        }
      }
    }
    break;
  }
  /* adjust VLA to current movie length */
  if(I->ViewElem) {
    VLASize(I->ViewElem,CViewElem,I->NFrame);
  }
  SceneSetFrame(G,1,0); /* force frame update */
  return 1;
}


/*========================================================================*/
void MovieAppendCommand(PyMOLGlobals * G, int frame, char *command)
{
  register CMovie *I = G->Movie;
  int a, len, cur_len;
  if((frame >= 0) && (frame < I->NFrame)) {
    len = strlen(command);
    cur_len = strlen(I->Cmd[frame]);
    if((unsigned) len > (sizeof(MovieCmdType) + cur_len - 1))
      len = sizeof(MovieCmdType) + cur_len - 1;
    for(a = 0; a < len; a++)
      I->Cmd[frame][cur_len + a] = command[a];
    I->Cmd[frame][cur_len + len] = 0;
  } else {
    PRINTFB(G, FB_Movie, FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",
      frame + 1 ENDFB(G);
  }
}


/*========================================================================*/
ImageType *MovieGetImage(PyMOLGlobals * G, int index)
{
  register CMovie *I = G->Movie;
  if((index >= 0) && (index < I->NImage))
    return (I->Image[index]);
  else
    return (NULL);
}


/*========================================================================*/
int MovieDefined(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  return (I->NFrame > 0);
}


/*========================================================================*/
int MovieGetLength(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  int len;
  if(!I->NFrame)
    len = -I->NImage;
  else
    len = I->NFrame;
  return (len);
}


/*========================================================================*/
void MovieClearImages(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  int a;

  PRINTFB(G, FB_Movie, FB_Blather)
    " MovieClearImages: clearing...\n" ENDFB(G);
  if(I->Image) {
    for(a = 0; a < I->NImage; a++) {
      if(I->Image[a]) {
        FreeP(I->Image[a]->data);
        FreeP(I->Image[a]);
        I->Image[a] = NULL;
      }
    }
  }
  I->NImage = 0;
  SceneInvalidate(G);
  SceneSuppressMovieFrame(G);
}


/*========================================================================*/
void MovieReset(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  MovieClearImages(G);

  VLAFreeP(I->Cmd);
  VLAFreeP(I->Sequence);
  VLAFreeP(I->ViewElem);

  I->NFrame = 0;
  I->MatrixFlag = false;
  I->Locked = false;
  I->Playing = false;
}


/*========================================================================*/
void MovieFree(PyMOLGlobals * G)
{
  register CMovie *I = G->Movie;
  MovieClearImages(G);
  VLAFree(I->Image);
  VLAFreeP(I->ViewElem);
  VLAFreeP(I->Cmd);
  VLAFreeP(I->Sequence);
  ScrollBarFree(I->ScrollBar);
  OrthoFreeBlock(G, I->Block);
  FreeP(G->Movie);
}

Block *MovieGetBlock(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  return (I->Block);
}


void MoviePrepareDrag(PyMOLGlobals *G, BlockRect * rect, 
                      CObject * obj, int mode, int x, int y, int nearest)
{
  CMovie *I = G->Movie;
  I->DragMode = mode;
  I->DragObj = obj;
  I->DragX = x;
  I->DragY = y;
  I->DragRect = *rect;
  if(I->DragColumn) {
    I->DragRect.top = I->Block->rect.top - 1;
    I->DragRect.bottom = I->Block->rect.bottom + 1;
  }
  I->DragStartFrame = ViewElemXtoFrame(G,I->ViewElem,rect,MovieGetLength(G),x,nearest);
  if(I->DragStartFrame > MovieGetLength(G))
    I->DragStartFrame = MovieGetLength(G);
  I->DragCurFrame = ViewElemXtoFrame(G,I->ViewElem,rect,MovieGetLength(G),x,nearest);
  I->DragNearest = nearest;
}

static int MovieClick(Block * block, int button, int x, int y, int mod)
{
  PyMOLGlobals *G = block->G;
  CMovie *I = G->Movie;
  int count = ExecutiveCountMotions(G);
  BlockRect rect = block->rect;
  rect.right -= I->LabelIndent;

  switch(button) {
  case P_GLUT_RIGHT_BUTTON:
    {
      int n_frame = MovieGetLength(G);
      if(mod == (cOrthoCTRL | cOrthoSHIFT))
        I->DragColumn = true;
      if(mod == (cOrthoSHIFT)) 
        ExecutiveMotionClick(G,&rect,cMovieDragModeCopyKey,count,x,y,false);
      else
        ExecutiveMotionClick(G,&rect,cMovieDragModeMoveKey,count,x,y,false);
      if(I->DragStartFrame<n_frame) {
        I->DragDraw = true;
        I->DragMenu = true;
        OrthoDirty(G);
      } else {
        ExecutiveMotionMenuActivate(G,&rect,count,false,x,y,I->DragColumn);
      }
    }
    break;
  case P_GLUT_LEFT_BUTTON:
    {
      switch(mod) {
      case cOrthoSHIFT: /* TEMPORAL SELECTIONS -- TO COME in PYMOL 1.3+ */
        break;
      case (cOrthoSHIFT | cOrthoCTRL): 
        I->DragColumn = true;
        ExecutiveMotionClick(G,&rect,cMovieDragModeInsDel,count,x,y, true);
        I->DragDraw = true;
        OrthoDirty(G);
        break;
      case cOrthoCTRL:
        ExecutiveMotionClick(G,&rect,cMovieDragModeInsDel,count,x,y, true);
        I->DragDraw = true;
        OrthoDirty(G);
        break;
      default:
        ScrollBarDoClick(I->ScrollBar, button, x, y, mod);
        break;
      }
    }
    break;
  case P_GLUT_MIDDLE_BUTTON:
    {
      switch(mod) {
      case (cOrthoCTRL | cOrthoSHIFT):
        I->DragColumn = true;
        /* intentional fall-through */
      case cOrthoCTRL: 
        I->DragDraw = true;
        ExecutiveMotionClick(G,&rect,cMovieDragModeOblate,count,x,y,false);
        break;
      default:
        ScrollBarDoClick(I->ScrollBar, button, x, y, mod);
        break;
      }
    }
    break;
  }
  return 1;
}

static int MovieDrag(Block * block, int x, int y, int mod)
{
  PyMOLGlobals *G = block->G;

  CMovie *I = G->Movie;
  if(I->DragMode) {
    I->DragDraw = ((y < (block->rect.top + 50)) && (y > (block->rect.bottom - 50)));
    switch(I->DragMode) {
    case cMovieDragModeMoveKey:
    case cMovieDragModeCopyKey:
      {
        int n_frame = MovieGetLength(G);
        I->DragCurFrame = ViewElemXtoFrame(G,I->ViewElem,&I->DragRect,n_frame,x,false);
        if(I->DragStartFrame<n_frame) {
          if((abs(x-I->DragX)>3) ||
             (abs(y-I->DragY)>5)) {
            I->DragMenu = false;
          }
          OrthoDirty(G);
        }
      }
      break;
    case cMovieDragModeOblate:
      I->DragCurFrame = ViewElemXtoFrame(G,I->ViewElem,&I->DragRect,MovieGetLength(G),x,false);
      OrthoDirty(G);
      break;
    case cMovieDragModeInsDel:
      I->DragCurFrame = ViewElemXtoFrame(G,I->ViewElem,&I->DragRect,MovieGetLength(G),x,true);
      OrthoDirty(G);
      break;
    }
  }
  return 1;
}

static int MovieRelease(Block * block, int button, int x, int y, int mod)
{
  PyMOLGlobals *G = block->G;
  CMovie *I = G->Movie;
  ScrollBarDoRelease(I->ScrollBar, button, x, y, mod);
  if(I->DragMode) {
    OrthoLineType buffer = "";
    OrthoLineType extra  = "";
    
    int n_frame = MovieGetLength(G);
    
    if(I->DragColumn) {
      sprintf(extra,",object=''");
    } else if(I->DragObj && ExecutiveValidateObjectPtr(G,I->DragObj,0)) {
      sprintf(extra,",object='%s'",I->DragObj->Name);
    } else {
      sprintf(extra,",object='none'");
    }
    switch(I->DragMode) {
    case cMovieDragModeMoveKey:
      if((I->DragCurFrame == I->DragStartFrame) && (I->DragMenu)) {
        int count = ExecutiveCountMotions(G);
        BlockRect rect = block->rect;
        rect.right -= I->LabelIndent;
        ExecutiveMotionMenuActivate(G,&rect,count,true,x,y,I->DragColumn);
        I->DragMenu = false;
      } else if(I->DragDraw &&
                (I->DragCurFrame!=I->DragStartFrame) && 
                (I->DragCurFrame >= 0) && 
                (I->DragCurFrame < n_frame)) {
        sprintf(buffer,"cmd.mmove(%d,%d,%d%s)", 1+I->DragCurFrame, 1+I->DragStartFrame, 1, extra);
      }
      break;
    case cMovieDragModeCopyKey:
      if((I->DragCurFrame == I->DragStartFrame) && (I->DragMenu)) {
        int count = ExecutiveCountMotions(G);
        BlockRect rect = block->rect;
        rect.right -= I->LabelIndent;
        ExecutiveMotionMenuActivate(G,&rect,count,true,x,y,I->DragColumn);
        I->DragMenu = false;
      } else if(I->DragDraw &&
                (I->DragCurFrame!=I->DragStartFrame) && 
                (I->DragCurFrame >= 0) && 
                (I->DragCurFrame < n_frame)) {
        sprintf(buffer,"cmd.mcopy(%d,%d,%d%s)", 1+I->DragCurFrame, 1+I->DragStartFrame, 1, extra);
      }
      break;
    case cMovieDragModeOblate:
      if(I->DragDraw) {
        int min_frame = (I->DragStartFrame < I->DragCurFrame) ? I->DragStartFrame : I->DragCurFrame;
        int max_frame = (I->DragStartFrame > I->DragCurFrame) ? I->DragStartFrame : I->DragCurFrame;
        if(min_frame<0) min_frame = 0;
        if(max_frame<0) max_frame = 0;
        if(min_frame>=n_frame) min_frame = n_frame - 1;
        if(max_frame>=n_frame) max_frame = n_frame - 1;
        if(I->DragColumn) {
          sprintf(extra,",object='same'");
        }
        sprintf(buffer,"cmd.mview('clear',first=%d,last=%d%s)", 
                1 + min_frame, 1 + max_frame, extra);
      }
      break;
    case cMovieDragModeInsDel:
      if(I->DragDraw) {
        if(I->DragCurFrame<0) 
          I->DragCurFrame = 0;
        if(I->DragCurFrame>I->DragStartFrame) {
          int first = I->DragStartFrame + 1;
          if(first<0) first = 0;
          sprintf(buffer,"cmd.minsert(%d,%d%s)", I->DragCurFrame - I->DragStartFrame, first, extra);
        } else {
          int first = I->DragCurFrame;
          if(first<0) first = 0;
          sprintf(buffer,"cmd.mdelete(%d,%d%s)", I->DragStartFrame - I->DragCurFrame, first+1, extra);
        }
      }
      break;
    }
    if(buffer[0]) {
      PParse(G, buffer);
      PFlush(G);
      PLog(G, buffer, cPLog_pym);
    }
  }
  I->DragMode = 0;
  I->DragDraw = false;
  I->DragMenu = false;
  I->DragColumn = false;

  return 1;
}

int MovieGetPanelHeight(PyMOLGlobals * G)
{
  int movie_panel = SettingGetGlobal_i(G, cSetting_movie_panel);
  CMovie *I = G->Movie;
  if(movie_panel != 0) {
    if(MovieGetLength(G)) {
      movie_panel = 1;
    } else {
      movie_panel = 0;
    }
  }
  
  if(movie_panel) {
    int row_height = SettingGetGlobal_i(G,cSetting_movie_panel_row_height);
    I->PanelActive = true;
    if(SettingGetGlobal_b(G, cSetting_presentation)) { 
      /* show camera line only when in presentation mode */
      return row_height;
    } else {
      return row_height * ExecutiveCountMotions(G); 
    }
  } else {
    I->PanelActive = false;
    return 0;
  }
}

void MovieDrawViewElem(PyMOLGlobals *G, BlockRect *rect,int frames)
{
  CMovie *I = G->Movie;
  if(I->ViewElem) {
    ViewElemDraw(G,I->ViewElem,rect,frames,"camera");
  }
}

static void MovieDraw(Block * block)
{
  PyMOLGlobals *G = block->G;
  CMovie *I = G->Movie;
  if(I->PanelActive) {
    int n_frame = MovieGetLength(G);
    int frame = SceneGetFrame(G);
    int count = ExecutiveCountMotions(G);
    BlockRect rect = block->rect;
    if(count) {
      rect.right -= I->LabelIndent;

      if(G->HaveGUI && G->ValidContext) {
        float black[3] = {0.0F,0.0F,0.0F};
        glColor3fv(black);
        glBegin(GL_POLYGON);
        glVertex2f(rect.right, rect.bottom);
        glVertex2f(rect.right, rect.top);
        glVertex2f(block->rect.right, rect.top);
        glVertex2f(block->rect.right, rect.bottom);
        glEnd();
      }

      if(!n_frame) {
        ScrollBarSetLimits(I->ScrollBar, 1, 1);
        ScrollBarSetValue(I->ScrollBar, 0);
      } else {
        float scroll_value = ScrollBarGetValue(I->ScrollBar);
        int new_frame = (int) (scroll_value + 0.5F);
        if(new_frame != frame) {
          frame = new_frame;
          SceneSetFrame(G, 7, frame);
        }
        if(!ScrollBarGrabbed(I->ScrollBar)) {
          ScrollBarSetValue(I->ScrollBar, frame);
        }
        ScrollBarSetLimits(I->ScrollBar, n_frame, 1);
      }
      ScrollBarSetBox(I->ScrollBar, rect.top,
                      rect.left, rect.bottom, rect.right);
      ScrollBarDoDraw(I->ScrollBar);
      ExecutiveMotionDraw(G,&rect,count);
      ScrollBarDrawHandle(I->ScrollBar, 0.35F);

      /* drag selection box */
      if(I->DragDraw) {

        float white[4] = {1.0F, 1.0F, 1.0F,0.5F};

        switch(I->DragMode) {
        case cMovieDragModeMoveKey:
        case cMovieDragModeCopyKey:
          {
            float grey[4] = {0.75F,0.75F,0.75f,0.5};
            if(I->DragStartFrame<n_frame) 
              ViewElemDrawBox(G,&I->DragRect, I->DragStartFrame, I->DragStartFrame+1, n_frame, white, false);        
            if((I->DragCurFrame>=0) && (I->DragCurFrame<n_frame)) {
              ViewElemDrawBox(G,&I->DragRect, I->DragCurFrame, I->DragCurFrame+1, n_frame, grey, true);
            }
          }
          break;
        case cMovieDragModeOblate:
          {
            float grey[4] = {0.75F,0.75F,0.75f,0.5};

            int min_frame = (I->DragStartFrame < I->DragCurFrame) ? I->DragStartFrame : I->DragCurFrame;
            int max_frame = (I->DragStartFrame > I->DragCurFrame) ? I->DragStartFrame : I->DragCurFrame;
            if(min_frame<0) min_frame = 0;
            if(max_frame<0) max_frame = 0;
            if(min_frame>=n_frame) min_frame = n_frame - 1;
            if(max_frame>=n_frame) max_frame = n_frame - 1;
            ViewElemDrawBox(G,&I->DragRect, min_frame, max_frame+1, n_frame, white, false);        
            ViewElemDrawBox(G,&I->DragRect, min_frame, max_frame+1, n_frame, grey, true);
          }
          break;
        case cMovieDragModeInsDel:
          if(I->DragCurFrame==I->DragStartFrame) {
            ViewElemDrawBox(G,&I->DragRect, I->DragStartFrame, I->DragStartFrame, n_frame, white, true);        
          } else if(I->DragCurFrame>=I->DragStartFrame) {
            float green[4] = {0.5F, 1.0F, 0.5F,0.5F};
            ViewElemDrawBox(G,&I->DragRect, I->DragStartFrame, I->DragCurFrame, n_frame, green, true);        
          } else {
            float red[4] = {1.0F, 0.5F, 0.5F,0.5F};          
            ViewElemDrawBox(G,&I->DragRect, I->DragCurFrame, I->DragStartFrame, n_frame, red, true);        
          }
          break;
        }

      }
    }
  }
}

void MovieSetScrollBarFrame(PyMOLGlobals * G, int frame)
{
  CMovie *I = G->Movie;
  if(!ScrollBarGrabbed(I->ScrollBar)) {
    ScrollBarSetValue(I->ScrollBar, frame);
  }
}

static void MovieReshape(Block * block, int width, int height)
{
  PyMOLGlobals *G = block->G;
  CMovie *I = G->Movie;
  BlockReshape(block, width, height);
  I->Width = block->rect.right - block->rect.left + 1;
  I->Height = block->rect.top - block->rect.bottom + 1;
  if(SettingGetGlobal_b(G, cSetting_presentation)) { 
    I->LabelIndent = 0;
  } else {
    I->LabelIndent = 8 * 8;
  }
}


/*========================================================================*/
int MovieInit(PyMOLGlobals * G)
{
  register CMovie *I = NULL;

  if((I = (G->Movie = Calloc(CMovie, 1)))) {
    int a;
    I->Block = OrthoNewBlock(G, NULL);
    I->Block->fRelease = MovieRelease;
    I->Block->fClick = MovieClick;
    I->Block->fDrag = MovieDrag;
    I->Block->fDraw = MovieDraw;
    I->Block->fReshape = MovieReshape;
    I->Block->active = true;
    I->ScrollBar = ScrollBarNew(G, true);
    OrthoAttach(G, I->Block, cOrthoTool);

    I->Playing = false;
    I->Image = VLAMalloc(10, sizeof(ImageType), 5, true);       /* auto-zero */
    I->Sequence = NULL;
    I->Cmd = NULL;
    I->ViewElem = NULL;
    I->NImage = 0;
    I->NFrame = 0;
    I->RecursionFlag = false;
    I->RealtimeFlag = true;
    for(a = 0; a < 16; a++)
      I->Matrix[a] = 0.0F;
    I->MatrixFlag = false;
    return 1;
  } else {
    return 0;
  }
}
