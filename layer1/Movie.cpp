
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
#include"Util2.h"
#include"Parse.h"
#include"PyMOL.h"
#include"ScrollBar.h"
#include"Menu.h"
#include"View.h"
#include"Seq.h"
#include"CGO.h"
#include"MovieScene.h"
#include"Feedback.h"

#define cMovieDragModeMoveKey   1
#define cMovieDragModeInsDel    2
#define cMovieDragModeCopyKey   3
#define cMovieDragModeOblate    4

CMovie::CMovie(PyMOLGlobals* G)
    : Block(G)
    , m_ScrollBar(G, true)
{
  active = true;
  OrthoAttach(G, this, cOrthoTool);
}

void MovieViewReinterpolate(PyMOLGlobals *G)
{
  float power  = SettingGetGlobal_f(G, cSetting_motion_power);
  float bias   = SettingGetGlobal_f(G, cSetting_motion_bias);
  float linear = SettingGetGlobal_f(G, cSetting_motion_linear);
  int hand     = SettingGetGlobal_i(G, cSetting_motion_hand);

  MovieView(G, 3, -1, -1, power, bias, 1, /* note simple always = 1 for camera motion...*/
            linear, 
            SettingGetGlobal_b(G,cSetting_movie_loop) ? 1 : 0 ,
            hand, 5, 1, nullptr, 0.5, -1, 1); 
}

int MovieXtoFrame(PyMOLGlobals *G, BlockRect *rect, int frames, int x, int nearest)
{
  return ViewElemXtoFrame(rect,frames,x,nearest);
}

void MovieViewTrim(PyMOLGlobals *G,int n_frame)
{
  CMovie *I = G->Movie;
  if(n_frame>=0) {
    I->Sequence.resize(n_frame);
    I->Cmd.resize(n_frame);
    I->ViewElem.resize(n_frame);
    I->NFrame = n_frame;
  }
}


int MovieViewModify(PyMOLGlobals *G,int action, int index, int count,int target, int freeze, int localize)
{
  CMovie *I = G->Movie;
  int ok = true;
  MovieClearImages(G);
  if( (ok = ViewElemModify(G,&I->ViewElem, action, index, count, target)) ) {
    switch(action) {
    case cViewElemModifyInsert:
      if (index >= 0 && index < I->NFrame) {
        I->Sequence.insert(index, count);
        I->Cmd.insert(I->Cmd.begin() + index, count, "");
        I->NFrame = VLAGetSize(I->Sequence);
        {
          int frame = SceneGetFrame(G);
          if (frame >= index) {
            SceneSetFrame(G, 0, frame + count);
          }
        }
      }
      break;
    case cViewElemModifyDelete:
      if (index >= 0 && index < I->NFrame) {
        I->Sequence.erase(index, count);
        int end_pos = std::min(index + count, static_cast<int>(I->Cmd.size()));
        I->Cmd.erase(I->Cmd.begin() + index, I->Cmd.begin() + end_pos);
        I->NFrame = VLAGetSize(I->Sequence);
      }
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
            I->Cmd[dst] = std::move(I->Cmd[src]);
            I->Cmd[src].clear();
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
            I->Cmd[dst] = I->Cmd[src];
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
  CMovie *I = G->Movie;
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
  CMovie *I = G->Movie;
  I->RealtimeFlag = realtime;
}

int MovieGetRealtime(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  return I->RealtimeFlag;
}

void MovieCopyPrepare(PyMOLGlobals * G, int *width, int *height, int *length)
{
  /* assumed locked api, blocked threads, and master thread on entry */

  /* writes the current movie sequence to a set of PNG files 
   * this routine can take a LOT of time -- 
   * TODO: develop an interrupt mechanism */

  int start, stop;
  CMovie *I = G->Movie;
  int nFrame;

  I->CacheSave = SettingGetGlobal_b(G, cSetting_cache_frames);
  I->OverlaySave = SettingGetGlobal_i(G, cSetting_overlay);
  if(!I->CacheSave)
    MovieClearImages(G);
  SettingSetGlobal_b(G, cSetting_cache_frames, 1);
  SettingSetGlobal_i(G, cSetting_overlay, 5);
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame = SceneGetNFrame(G, nullptr);
  }
  start = 0;
  stop = nFrame;
  if((start != 0) || (stop != (nFrame + 1)))
    SceneSetFrame(G, 0, 0);
  MoviePlay(G, cMoviePlay);
  VecCheck(I->Image, nFrame);
  SceneGetWidthHeight(G, width, height);
  {
    int uniform_height = -1;
    int uniform_width = -1;
    int uniform_flag = false;
    int scene_match = true;
    int a;
    /* make sure all the movie frames match the screen size or are pre-rendered and are already the same size */
    for(a = 0; a < nFrame; a++) {
      const pymol::Image* image = I->Image[a].get();
      if(image) {
        if((image->getHeight() != *height) || (image->getWidth() != *width)) {
          scene_match = false;
          if(uniform_height < 0) {
            uniform_height = image->getHeight();
            uniform_width = image->getWidth();
          } else {
            if((image->getHeight() != uniform_height) || (image->getWidth() != uniform_width))
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
  CMovie *I = G->Movie;
  I->RecursionFlag = true;
  PFlush(G);
  I->RecursionFlag = false;
}

int MovieCopyFrame(PyMOLGlobals * G, int frame, int width, int height, int rowbytes,
                   void *ptr)
{
  CMovie *I = G->Movie;
  int result = false;
  int nFrame;
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame = SceneGetNFrame(G, nullptr);
  }

  if((frame < nFrame) && (ptr)) {
    int a = frame;
    int i;
    SceneSetFrame(G, 0, a);
    MovieDoFrameCommand(G, a);
    MovieFlushCommands(G);
    i = MovieFrameToImage(G, a);
    VecCheck(I->Image, i);
    if(!I->Image[i]) {
      SceneUpdate(G, false);
      SceneMakeMovieImage(G, false, false, cSceneImage_Default);
    }
    if(!I->Image[i]) {
      PRINTFB(G, FB_Movie, FB_Errors)
        "MoviePNG-Error: Missing rendered image.\n" ENDFB(G);
    } else {
      if((I->Image[i]->getHeight() == height) && (I->Image[i]->getWidth() == width)) {
        unsigned char *srcImage = I->Image[i]->bits();
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
        I->Image[i] = nullptr;
      }
    }
  }
  return result;
}

int MoviePurgeFrame(PyMOLGlobals * G, int frame)
{
  CMovie *I = G->Movie;
  int result = false;
  int nFrame;
  int i;
  nFrame = I->NFrame;
  if(!nFrame) {
    nFrame = SceneGetNFrame(G, nullptr);
  }
  if(!I->CacheSave) {
    if(frame < nFrame) {
      int a = frame;
      i = MovieFrameToImage(G, a);
      VecCheck(I->Image, i);
      if(I->Image[i]) {
        I->Image[i] = nullptr;
        result = true;
      }
    }
  }
  return result;
}

void MovieCopyFinish(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  SceneInvalidate(G);           /* important */
  SettingSetGlobal_b(G, cSetting_cache_frames, I->CacheSave);
  SettingSetGlobal_i(G, cSetting_overlay, I->OverlaySave);
  MoviePlay(G, cMovieStop);
  if(!I->CacheSave) {
    MovieClearImages(G);
  }
}

int MovieLocked(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  return (I->Locked);
}

void MovieSetLock(PyMOLGlobals * G, int lock)
{
  CMovie *I = G->Movie;
  I->Locked = lock;
}


/*========================================================================*/
void MovieDump(PyMOLGlobals * G)
{
  int a;
  CMovie *I = G->Movie;
  int flag = false;

  for(a = 0; a < I->NFrame; a++) {
    if(!I->Cmd[a].empty()) {
      flag = true;
      break;
    }
  }
  if(flag && I->NFrame) {
    PRINTFB(G, FB_Movie, FB_Results)
      " Movie: General Purpose Commands:\n" ENDFB(G);
    for(a = 0; a < I->NFrame; a++) {
      if(!I->Cmd[a].empty()) {
        auto buffer = pymol::string_format("%5d: %s\n", a + 1, I->Cmd[a]);
        OrthoAddOutput(G, buffer.c_str());
      }
    }
  } else {
    PRINTFB(G, FB_Movie, FB_Results)
      " Movie: No movie commands are defined.\n" ENDFB(G);
  }
}

static int MovieCmdFromPyList(PyMOLGlobals * G, PyObject * list, int *warning)
{

  CMovie *I = G->Movie;
  int ok = true;
  int a;
  int warn = false;

  if(ok)
    ok = (list != nullptr);
  if(ok)
    ok = PyList_Check(list);

  for(a = 0; a < I->NFrame; a++) {
    if(ok)
      ok = PConvFromPyListItem(G, list, a, I->Cmd[a]);
    if(ok)
      warn = (warn || !I->Cmd[a].empty());
  }
  *warning = warn;
  return (ok);

}

/*========================================================================*/
int MovieFromPyList(PyMOLGlobals * G, PyObject * list, int *warning)
{
  int ok = true;
  CMovie *I = G->Movie;
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
    I->Sequence = pymol::vla<int>(I->NFrame);
    I->Cmd = std::vector<std::string>(I->NFrame);
    if(ok)
      ok = PConvPyListToIntArrayInPlace(PyList_GetItem(list, 4), I->Sequence.data(), I->NFrame);
    if(ok)
      ok = MovieCmdFromPyList(G, PyList_GetItem(list, 5), warning);
    if((*warning) && G->Security) {
      MovieSetLock(G, true);
    }
  }
  if(ok && (ll > 6)) {
    PyObject *tmp;
    VLAFreeP(I->ViewElem);
    I->ViewElem = nullptr;
    tmp = PyList_GetItem(list, 6);
    if(tmp && !(tmp == Py_None))
      ok = ViewElemVLAFromPyList(G, tmp, &I->ViewElem, I->NFrame);
  }
  if(!ok) {
    MovieReset(G);
  } else if(MovieDefined(G)) {
    OrthoReshape(G,-1,-1,true);
    SceneCountFrames(G);
  }
  return (ok);
}


/*========================================================================*/
static PyObject *MovieCmdAsPyList(PyMOLGlobals * G)
{

  CMovie *I = G->Movie;
  PyObject *result = nullptr;
  int a;

  result = PyList_New(I->NFrame);
  if(result)
    for(a = 0; a < I->NFrame; a++) {
      PyList_SetItem(result, a, PyString_FromString(I->Cmd[a].data()));
    }
  return (PConvAutoNone(result));

}

/*========================================================================*/
PyObject *MovieAsPyList(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  PyObject *result = nullptr;

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
  if(!I->Cmd.empty()) {
    PyList_SetItem(result, 5, MovieCmdAsPyList(G));
  } else {
    PyList_SetItem(result, 5, PConvAutoNone(NULL));
  }
  if(I->ViewElem) {
    PyList_SetItem(result, 6, ViewElemVLAAsPyList(G, I->ViewElem, I->NFrame));
  } else {
    PyList_SetItem(result, 6, PConvAutoNone(NULL));
  }


  /*   pymol::Image *Image;
       int *Sequence;
       MovieCmdType *Cmd;
       int NImage,NFrame;
       unsigned Width,Height;
       int MatrixFlag;
       float Matrix[16];
       int Playing;
  */
  return (PConvAutoNone(result));
}


/*========================================================================*/
int MoviePlaying(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
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
  CMovie *I = G->Movie;
  switch (cmd) {
  case cMovieToggle:
    I->Playing = !I->Playing;
    if(I->Playing && !SettingGetGlobal_b(G, cSetting_movie_loop)) {
      /* if not looping, and at end of movie, then automatically rewind
         and force execution of the first movie command */
      if((SettingGetGlobal_i(G, cSetting_frame)) == (SceneGetNFrame(G, nullptr))) {
        SceneSetFrame(G, 7, 0);
      }
    }
    break;
  case cMovieStop:
    I->Playing = false;
    break;
  case cMoviePlay:
    if(!SettingGetGlobal_b(G, cSetting_movie_loop)) {
      /* if not looping, and at end of movie, then automatically rewind
         and force execution of the first movie command */
      if((SettingGetGlobal_i(G, cSetting_frame)) == (SceneGetNFrame(G, nullptr))) {
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
  CMovie *I = G->Movie;
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
    M->save = SettingGetGlobal_b(G, cSetting_cache_frames);
    if(!M->save)
      MovieClearImages(G);
    SettingSetGlobal_b(G, cSetting_cache_frames, 1);
    OrthoBusyPrime(G);
    M->nFrame = I->NFrame;
    if(!M->nFrame) {
      M->nFrame = SceneGetNFrame(G, nullptr);
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
      auto buffer = pymol::string_format("Creating movie (%d frames)...", M->nFrame);
      OrthoBusyMessage(G, buffer.c_str());
    }
    if((M->start != 0) || (M->stop != (M->nFrame + 1)))
      SceneSetFrame(G, 0, 0);
    MoviePlay(G, cMoviePlay);
    VecCheck(I->Image, M->nFrame);
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
        M->fname = pymol::string_format("%s%04d.ppm", M->prefix.c_str(), M->frame + 1);
        break;
      case cMyPNG_FormatPNG:
      default:
        M->fname = pymol::string_format("%s%04d.png", M->prefix.c_str(), M->frame + 1);
        break;
      }

      if(M->missing_only) {
        FILE *tmp = fopen(M->fname.c_str(), "rb");
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
    VecCheck(I->Image, M->image);
    if((M->frame >= M->start) &&        /* only render frames in the specified interval... */
       (M->frame <= M->stop) && (M->file_missing)) {    /* ...that don't already exist */
      if(!I->Image[M->image]) {
        SceneUpdate(G, false);
        if(SceneMakeMovieImage(G, false, M->modal, M->mode, M->width, M->height)
            || !M->modal) {
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
      if (!MyPNGWrite(M->fname.c_str(), *I->Image[M->image],
              SettingGetGlobal_f(G, cSetting_image_dots_per_inch), M->format,
              M->quiet, SettingGetGlobal_f(G, cSetting_png_screen_gamma),
              SettingGetGlobal_f(G, cSetting_png_file_gamma))) {
        PRINTFB(G, FB_Movie, FB_Errors)
          " MoviePNG-Error: unable to write '%s'\n", M->fname.c_str() ENDFB(G);
      }
      ExecutiveDrawNow(G);
      OrthoBusySlow(G, M->frame, M->nFrame);
      if(G->HaveGUI)
        PyMOL_SwapBuffers(G->PyMOL);
      PRINTFB(G, FB_Movie, FB_Debugging)
        " MoviePNG-DEBUG: i = %d, I->Image[image] = %p\n", M->image,
        I->Image[M->image]->bits() ENDFB(G);
    }
    if(I->Image[M->image]) {
      I->Image[M->image] = nullptr;
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
    SettingSetGlobal_b(G, cSetting_cache_frames, M->save);
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
  CMovie *I = G->Movie;
  MovieModalPNG(G, I, &I->Modal);
  if(!I->Modal.complete)        /* force modalic return until job is done */
    PyMOL_SetModalDraw(G->PyMOL, (PyMOLModalDrawFn *) MovieModalDraw);
}

int MoviePNG(PyMOLGlobals * G, const char* prefix, int save, int start,
             int stop, int missing_only, int modal, int format, int mode, int quiet,
             int width, int height)
{
  /* assumes locked api, blocked threads, and master thread on entry */
  CMovie *I = G->Movie;

  /* new routine allows PyMOL to workaround XP and Vista Windowing
     issues */

  CMovieModal *M = &I->Modal;

  *M = CMovieModal();

  mode = SceneValidateImageMode(G, mode, width || height);

  /* default behavior is to go modal unless we're ray tracing */
  if(modal < 0 && mode == cSceneImage_Ray) {
    modal = 0;
  }

  M->prefix = prefix;
  M->save = save;
  M->start = start;
  M->stop = stop;
  M->missing_only = missing_only;
  M->stage = 0;
  M->format = format;
  M->mode = mode;
  M->quiet = quiet;
  M->width = width;
  M->height = height;

  if(SettingGetGlobal_b(G, cSetting_seq_view)) {
    PRINTFB(G, FB_Movie, FB_Warnings)
      " MoviePNG-Warning: disabling seq_view, may conflict with movie export\n" ENDFB(G);
    SettingSetGlobal_b(G, cSetting_seq_view, 0);
    // force viewport update
    SeqChanged(G);
    OrthoDoDraw(G, OrthoRenderMode::Main);
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
void MovieSet(PyMOLGlobals* G, pymol::zstring_view specification,
    int start_from, bool freeze)
{
  MovieAppendSequence(G, specification.c_str(), start_from, freeze);
  SceneCountFrames(G);

  // fix for PYMOL-1465
  // force GUI update for movie panel
  if(G->HaveGUI)
  OrthoReshape(G, -1, -1, false);
}

/*========================================================================*/
void MovieAppendSequence(PyMOLGlobals* G, const char* str, int start_from, bool freeze)
{
  CMovie *I = G->Movie;
  int c = 0;
  int i;
  /*  int old_NFrame = I->NFrame; */
  const char *s;
  char number[20];

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
    I->Cmd.clear();
    VLAFreeP(I->ViewElem);
    I->NFrame = 0;
  } else {
    // to clear
    I->Sequence.resize(start_from);
    I->Cmd.resize(start_from);
    I->ViewElem.resize(start_from);

    // append (c - start_from) zero-initialized elements
    I->Sequence.resize(c);
    I->Cmd.resize(c);
    I->ViewElem.resize(c);
  }

  if(c && str[0]) {             /* not just a reset */
    for(i = start_from; i < c; i++)
      I->Cmd[i].clear();
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

  // fixes PYMOL-2710
  MovieClearImages(G);

  I->Image.resize(I->NFrame);
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
  int single_image = SettingGetGlobal_b(G, cSetting_single_image);
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
  CMovie *I = G->Movie;
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
void MovieSetImage(PyMOLGlobals * G, int index, std::shared_ptr<pymol::Image> image)
{
  CMovie *I = G->Movie;

  PRINTFB(G, FB_Movie, FB_Blather)
    " MovieSetImage: setting movie image %d\n", index + 1 ENDFB(G);

  VecCheck(I->Image, index);
  I->Image[index] = image;
  if(I->NImage < (index + 1))
    I->NImage = index + 1;
}

int MovieSeekScene(PyMOLGlobals * G, int loop)
{
  CMovie *I = G->Movie;
  int result = -1;
  OVreturn_word ret;
  const char *scene_name = SettingGetGlobal_s(G,cSetting_scene_current_name);
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
  CMovie *I = G->Movie;
  if(frame == 0)
    MovieMatrix(G, cMovieMatrixRecall);
  if(!I->Locked) {
    if((frame >= 0) && (frame < I->NFrame)) {
      if(!I->Cmd[frame].empty()) {
        if(!I->RecursionFlag) {
          PParse(G, I->Cmd[frame].data());
        }
      }
      if(I->ViewElem) {
        if(I->ViewElem[frame].scene_flag) {
          char *st = OVLexicon_FetchCString(G->Lexicon, I->ViewElem[frame].scene_name);
          if(strcmp(st, SettingGetGlobal_s(G, cSetting_scene_current_name))) {
            MovieSceneRecall(G, st, 0.0,
                /* view */ false, true, true, true,
                /* frame */ false);
          }
        }
        SceneFromViewElem(G, I->ViewElem + frame, true);
      }
    }
  }
}


/*========================================================================*/
void MovieSetCommand(PyMOLGlobals* G, int frame, const char* command)
{
  CMovie *I = G->Movie;
  if((frame >= 0) && (frame < I->NFrame)) {
    I->Cmd[frame] = command;
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
              const char *scene_name, float scene_cut, int state, int quiet)
{
  CMovie *I = G->Movie;
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
        last = SceneGetNFrame(G, nullptr) - 1;
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
            scene_name = nullptr;
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
      CViewElem *first_view = nullptr, *last_view = nullptr;
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
            last_view = nullptr;
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
      I->ViewElem = pymol::vla<CViewElem>(size);
    }
    break;
  case 6:                      /* uninterpolate */
    if(I->ViewElem) {
      if(first < 0)
        first = 0;
      if(last < 0) {
        last = SceneGetNFrame(G, nullptr) - 1;
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
void MovieAppendCommand(PyMOLGlobals * G, int frame, const char* command)
{
  CMovie *I = G->Movie;
  if((frame >= 0) && (frame < I->NFrame)) {
    I->Cmd[frame] += command;
  } else {
    PRINTFB(G, FB_Movie, FB_Errors)
      " Movie-Error: frame %d does not exist.  Use 'mset' to define movie first.\n",
      frame + 1 ENDFB(G);
  }
}


/*========================================================================*/
std::shared_ptr<pymol::Image> MovieGetImage(PyMOLGlobals * G, int index)
{
  CMovie *I = G->Movie;
  if((index >= 0) && (index < I->NImage))
    return I->Image[index];
  else
    return nullptr;
}


/*========================================================================*/
int MovieDefined(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  return (I->NFrame > 0);
}


/*========================================================================*/
int MovieGetLength(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  int len;
  if(!I->NFrame)
    len = -I->NImage;
  else
    len = I->NFrame;
  return (len);
}


/*========================================================================*/
void MovieClearImages(PyMOLGlobals * G, CMovie* I)
{
  I->Image.clear();
  I->NImage = 0;
  SceneInvalidate(G);
  SceneSuppressMovieFrame(G);
}

void MovieClearImages(PyMOLGlobals * G)
{
  PRINTFB(G, FB_Movie, FB_Blather)
    " MovieClearImages: clearing...\n" ENDFB(G);
  MovieClearImages(G, G->Movie);
}

/*========================================================================*/
void MovieReset(PyMOLGlobals * G)
{
  CMovie *I = G->Movie;
  MovieClearImages(G);

  I->Cmd.clear();
  VLAFreeP(I->Sequence);
  VLAFreeP(I->ViewElem);

  I->NFrame = 0;
  I->MatrixFlag = false;
  I->Locked = false;
  I->Playing = false;
}


/*========================================================================*/
CMovie::~CMovie()
{
  MovieClearImages(m_G, this);
}

Block *MovieGetBlock(PyMOLGlobals * G)
{
  return G->Movie;
}


void MoviePrepareDrag(PyMOLGlobals *G, BlockRect * rect, 
                      pymol::CObject * obj, int mode, int x, int y, int nearest)
{
  CMovie *I = G->Movie;
  I->DragMode = mode;
  I->DragObj = obj;
  I->DragX = x;
  I->DragY = y;
  I->DragRect = *rect;
  if(I->DragColumn) {
    I->DragRect.top = I->rect.top - 1;
    I->DragRect.bottom = I->rect.bottom + 1;
  }
  I->DragStartFrame = ViewElemXtoFrame(rect,MovieGetLength(G),x,nearest);
  if(I->DragStartFrame > MovieGetLength(G))
    I->DragStartFrame = MovieGetLength(G);
  I->DragCurFrame = ViewElemXtoFrame(rect,MovieGetLength(G),x,nearest);
  I->DragNearest = nearest;
}

int CMovie::click(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CMovie *I = G->Movie;
  int count = ExecutiveCountMotions(G);
  short scrolldir = 1;
  BlockRect tmpRect = rect;
  tmpRect.right -= I->LabelIndent;

  switch(button) {
  case P_GLUT_RIGHT_BUTTON:
    {
      int n_frame = MovieGetLength(G);
      if(mod == (cOrthoCTRL | cOrthoSHIFT))
        I->DragColumn = true;
      if(mod == (cOrthoSHIFT)) 
        ExecutiveMotionClick(G,&tmpRect,cMovieDragModeCopyKey,count,x,y,false);
      else
        ExecutiveMotionClick(G,&tmpRect,cMovieDragModeMoveKey,count,x,y,false);
      if(I->DragStartFrame<n_frame) {
        I->DragDraw = true;
        I->DragMenu = true;
        OrthoDirty(G);
      } else {
        ExecutiveMotionMenuActivate(G,&tmpRect,count,false,x,y,I->DragColumn);
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
        ExecutiveMotionClick(G,&tmpRect,cMovieDragModeInsDel,count,x,y, true);
        I->DragDraw = true;
        OrthoDirty(G);
        break;
      case cOrthoCTRL:
        ExecutiveMotionClick(G,&tmpRect,cMovieDragModeInsDel,count,x,y, true);
        I->DragDraw = true;
        OrthoDirty(G);
        break;
      default:
        I->m_ScrollBar.click(button, x, y, mod);
        {
	  SceneSetFrame(G, 7, I->m_ScrollBar.getValue());
	}
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
        ExecutiveMotionClick(G,&tmpRect,cMovieDragModeOblate,count,x,y,false);
        break;
      default:
        I->m_ScrollBar.click(button, x, y, mod);
        break;
      }
    }
    break;
  case P_GLUT_BUTTON_SCROLL_FORWARD:
    scrolldir = -1;
  case P_GLUT_BUTTON_SCROLL_BACKWARD:
    switch(mod) {
      case (cOrthoCTRL | cOrthoSHIFT):
        SettingSetGlobal_i(G, cSetting_movie_panel_row_height,
            SettingGetGlobal_i(G, cSetting_movie_panel_row_height) - scrolldir);
        OrthoReshape(G,-1,-1,true);
        break;
      default:
        SceneSetFrame(G, 5, scrolldir);
    }
    break;
  }
  return 1;
}

int CMovie::drag(int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;

  CMovie *I = this; // TODO: Remove all I's in Movie refactor
  if(I->DragMode) {
    I->DragDraw = ((y < (rect.top + 50)) && (y > (rect.bottom - 50)));
    switch(I->DragMode) {
    case cMovieDragModeMoveKey:
    case cMovieDragModeCopyKey:
      {
        int n_frame = MovieGetLength(G);
        I->DragCurFrame = ViewElemXtoFrame(&I->DragRect,n_frame,x,false);
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
      I->DragCurFrame = ViewElemXtoFrame(&I->DragRect,MovieGetLength(G),x,false);
      OrthoDirty(G);
      break;
    case cMovieDragModeInsDel:
      I->DragCurFrame = ViewElemXtoFrame(&I->DragRect,MovieGetLength(G),x,true);
      OrthoDirty(G);
      break;
    }
  }
  return 1;
}

int CMovie::release(int button, int x, int y, int mod)
{
  PyMOLGlobals *G = m_G;
  CMovie *I = G->Movie;
  I->m_ScrollBar.release(button, x, y, mod);
  if(I->DragMode) {
    std::string buffer;
    std::string extra;
    
    int n_frame = MovieGetLength(G);
    
    if(I->DragColumn) {
      extra = ",object=''";
    } else if(I->DragObj && ExecutiveValidateObjectPtr(G,I->DragObj,0)) {
      extra = pymol::string_format(",object='%s'",I->DragObj->Name);
    } else {
      extra = ",object='none'";
    }
    switch(I->DragMode) {
    case cMovieDragModeMoveKey:
      if((I->DragCurFrame == I->DragStartFrame) && (I->DragMenu)) {
        int count = ExecutiveCountMotions(G);
        BlockRect tmpRect = rect;
        tmpRect.right -= I->LabelIndent;
        ExecutiveMotionMenuActivate(G,&tmpRect,count,true,x,y,I->DragColumn);
        I->DragMenu = false;
      } else if(I->DragDraw &&
                (I->DragCurFrame!=I->DragStartFrame) && 
                (I->DragCurFrame >= 0) && 
                (I->DragCurFrame < n_frame)) {
        buffer = pymol::string_format("cmd.mmove(%d,%d,%d%s)", 1+I->DragCurFrame, 1+I->DragStartFrame, 1, extra);
      }
      break;
    case cMovieDragModeCopyKey:
      if((I->DragCurFrame == I->DragStartFrame) && (I->DragMenu)) {
        int count = ExecutiveCountMotions(G);
        BlockRect tmpRect = rect;
        tmpRect.right -= I->LabelIndent;
        ExecutiveMotionMenuActivate(G,&tmpRect,count,true,x,y,I->DragColumn);
        I->DragMenu = false;
      } else if(I->DragDraw &&
                (I->DragCurFrame!=I->DragStartFrame) && 
                (I->DragCurFrame >= 0) && 
                (I->DragCurFrame < n_frame)) {
        buffer = pymol::string_format("cmd.mcopy(%d,%d,%d%s)", 1+I->DragCurFrame, 1+I->DragStartFrame, 1, extra);
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
          extra = ",object='same'";
        }
        buffer = pymol::string_format("cmd.mview('clear',first=%d,last=%d%s)",
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
          buffer = pymol::string_format("cmd.minsert(%d,%d%s)", I->DragCurFrame - I->DragStartFrame, first, extra);
        } else {
          int first = I->DragCurFrame;
          if(first<0) first = 0;
          buffer = pymol::string_format("cmd.mdelete(%d,%d%s)", I->DragStartFrame - I->DragCurFrame, first+1, extra);
        }
      }
      break;
    }
    if(!buffer.empty()) {
      PParse(G, buffer.c_str());
      PFlush(G);
      PLog(G, buffer.c_str(), cPLog_pym);
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
      movie_panel = SceneGetNFrame(G) > 1;
    }
  }
  
  if(movie_panel) {
    int row_height = DIP2PIXEL(SettingGetGlobal_i(G,cSetting_movie_panel_row_height));
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

void MovieDrawViewElem(PyMOLGlobals *G, BlockRect *rect,int frames , CGO *orthoCGO)
{
  CMovie *I = G->Movie;
  if(I->ViewElem) {
    ViewElemDraw(G,I->ViewElem,rect,frames,"camera", orthoCGO);
  }
}

bool CMovie::fastDraw(CGO* orthoCGO)
{
  return true;
}

void CMovie::draw(CGO* orthoCGO)
{
  PyMOLGlobals *G = m_G;
  CMovie *I = G->Movie;
  if(I->PanelActive) {
    int n_frame = SceneGetNFrame(G);
    int frame = SceneGetFrame(G);
    int count = ExecutiveCountMotions(G);
    BlockRect tmpRect = rect;
    if(count) {
      tmpRect.right -= I->LabelIndent;

      if(G->HaveGUI && G->ValidContext) {
        float black[3] = {0.0F,0.0F,0.0F};
	if (orthoCGO){
	  CGOColorv(orthoCGO, black);
	  CGOBegin(orthoCGO, GL_TRIANGLE_STRIP);
	  CGOVertex(orthoCGO, tmpRect.right, tmpRect.bottom, 0.f);
	  CGOVertex(orthoCGO, tmpRect.right, tmpRect.top, 0.f);
	  CGOVertex(orthoCGO, rect.right, tmpRect.bottom, 0.f);
	  CGOVertex(orthoCGO, rect.right, tmpRect.top, 0.f);
	  CGOEnd(orthoCGO);
	} else {
	  glColor3fv(black);
	  glBegin(GL_POLYGON);
	  glVertex2f(tmpRect.right, tmpRect.bottom);
	  glVertex2f(tmpRect.right, tmpRect.top);
	  glVertex2f(rect.right, tmpRect.top);
	  glVertex2f(rect.right, tmpRect.bottom);
	  glEnd();
	}
      }

      if(!n_frame) {
        I->m_ScrollBar.setLimits(1, 1);
        I->m_ScrollBar.setValue(0);
      } else {
        float scroll_value = I->m_ScrollBar.getValue();
        int new_frame = (int) (scroll_value + 0.5F);
        if(I->m_ScrollBar.grabbed()) {
	  if(new_frame != frame) {
	    frame = new_frame;
	    SceneSetFrame(G, 7, frame);
	  }
        }
        I->m_ScrollBar.setLimits(n_frame, 1);
      }
      I->m_ScrollBar.setBox(tmpRect.top, tmpRect.left,
                             tmpRect.bottom, tmpRect.right);
      {
	I->m_ScrollBar.draw(orthoCGO);
	ExecutiveMotionDraw(G,&tmpRect,count, orthoCGO);
	I->m_ScrollBar.drawHandle(0.35F, orthoCGO);
      }

      /* drag selection box */
      if(I->DragDraw) {

        float white[4] = {1.0F, 1.0F, 1.0F,0.5F};

        switch(I->DragMode) {
        case cMovieDragModeMoveKey:
        case cMovieDragModeCopyKey:
          {
            float grey[4] = {0.75F,0.75F,0.75f,0.5};
            if(I->DragStartFrame<n_frame) 
              ViewElemDrawBox(G,&I->DragRect, I->DragStartFrame, I->DragStartFrame+1, n_frame, white, false, orthoCGO);
            if((I->DragCurFrame>=0) && (I->DragCurFrame<n_frame)) {
              ViewElemDrawBox(G,&I->DragRect, I->DragCurFrame, I->DragCurFrame+1, n_frame, grey, true, orthoCGO);
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
            ViewElemDrawBox(G,&I->DragRect, min_frame, max_frame+1, n_frame, white, false, orthoCGO);
            ViewElemDrawBox(G,&I->DragRect, min_frame, max_frame+1, n_frame, grey, true, orthoCGO);
          }
          break;
        case cMovieDragModeInsDel:
          if(I->DragCurFrame==I->DragStartFrame) {
            ViewElemDrawBox(G,&I->DragRect, I->DragStartFrame, I->DragStartFrame, n_frame, white, true, orthoCGO);
          } else if(I->DragCurFrame>=I->DragStartFrame) {
            float green[4] = {0.5F, 1.0F, 0.5F,0.5F};
            ViewElemDrawBox(G,&I->DragRect, I->DragStartFrame, I->DragCurFrame, n_frame, green, true, orthoCGO);
          } else {
            float red[4] = {1.0F, 0.5F, 0.5F,0.5F};          
            ViewElemDrawBox(G,&I->DragRect, I->DragCurFrame, I->DragStartFrame, n_frame, red, true, orthoCGO);
          }
          break;
        }

      }

      if (!ViewElem) {
        ViewElemDrawLabel(G, "states", &tmpRect, orthoCGO);
      }
    }
  }
}

void MovieSetScrollBarFrame(PyMOLGlobals * G, int frame)
{
  CMovie *I = G->Movie;
  if(!I->m_ScrollBar.grabbed()) {
    I->m_ScrollBar.setValue(frame);
  }
}

void CMovie::reshape(int width, int height)
{
  PyMOLGlobals *G = m_G;
  CMovie *I = G->Movie;
  Block::reshape(width, height);
  I->Width = rect.right - rect.left + 1;
  I->Height = rect.top - rect.bottom + 1;
  if(SettingGetGlobal_b(G, cSetting_presentation)) { 
    I->LabelIndent = 0;
  } else {
    I->LabelIndent = DIP2PIXEL(8 * 8);
  }
}

