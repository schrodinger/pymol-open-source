
/* 
   A* -------------------------------------------------------------------
   B* This file contains source code for the PyMOL computer program
   C* Copyright (c) Schrodinger, LLC. 
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
#include "os_python.h"

#include "os_predef.h"


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#include <signal.h>
#include <tchar.h>
#include <stdafx.h>
#endif

/* END PROPRIETARY CODE SEGMENT */

#include "os_std.h"
#include "os_gl.h"

#ifdef _PYMOL_MODULE
#ifdef _DRI_WORKAROUND
#include <dlfcn.h>
#endif
#endif

#include "PyMOLGlobals.h"
#include "PyMOL.h"
#include "PyMOLOptions.h"

#include "MemoryDebug.h"
#include "Base.h"
#include "main.h"

#include "P.h"
#include "PConv.h"
#include "Ortho.h"
#include "Scene.h"
#include "Seq.h"
#include "Util.h"
#include "Control.h"
#include "Movie.h"
#include "Shader.h"
#ifdef _PYMOL_NO_MAIN

int MainSavingUnderWhileIdle(void)
{
  return 0;
}

#ifdef _MACPYMOL_XCODE

void MainBlock(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PBlock(G);
}

void MainUnblock(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PUnblock(G);
}

int MainLockAPIAsGlut(int a)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  return PLockAPIAsGlut(G, a);
}

void MainUnlockAPIAsGlut()
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PUnlockAPIAsGlut(G);
}

int MainFeedbackOut(char *st)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  return OrthoFeedbackOut(G, st);
}

void MainRunCommand(char *str1)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  if(PLockAPIAsGlut(G, true)) {

    if(str1[0] != '_') {        /* suppress internal call-backs */
      if(strncmp(str1, "cmd._", 5)) {
        OrthoAddOutput(G, "PyMOL>");
        OrthoAddOutput(G, str1);
        OrthoNewLine(G, NULL, true);
        if(WordMatch(G, str1, "quit", true) == 0)       /* don't log quit */
          PLog(G, str1, cPLog_pml);
      }
      PParse(G, str1);
    } else if(str1[1] == ' ') { /* "_ command" suppresses echoing of command, but it is still logged */
      if(WordMatch(G, str1 + 2, "quit", true) >= 0)     /* don't log quit */
        PLog(G, str1 + 2, cPLog_pml);
      PParse(G, str1 + 2);
    } else {
      PParse(G, str1);
    }
    PUnlockAPIAsGlut(G);
  }
}

PyObject *MainGetStringResult(char *str)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PyObject *result;
  result = PyRun_String(str, Py_eval_input, G->P_inst->dict, G->P_inst->dict);
  return (result);
}

void MainRunString(char *str)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PBlock(G);
  PLockStatus(G);
  PRunStringModule(G, str);
  PUnlockStatus(G);
  PUnblock(G);
}

void MainMovieCopyPrepare(int *width, int *height, int *length)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    MovieCopyPrepare(G, width, height, length);
    MovieSetRealtime(G, false);
    PUnlockAPIAsGlut(G);
  }
}

int MainMovieCopyFrame(int frame, int width, int height, int rowbytes, void *ptr)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  int result = false;
  if(PLockAPIAsGlut(G, true)) {
    result = MovieCopyFrame(G, frame, width, height, rowbytes, ptr);
    PUnlockAPIAsGlut(G);
  }
  return result;
}

int MainMoviePurgeFrame(int frame)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  int result = false;
  if(PLockAPIAsGlut(G, true)) {
    result = MoviePurgeFrame(G, frame);
    PUnlockAPIAsGlut(G);
  }
  return result;
}

void MainMovieCopyFinish(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    MovieCopyFinish(G);
    MovieSetRealtime(G, true);
    PUnlockAPIAsGlut(G);
  }
}

void MainFlushAsync(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    PFlush(G);
    PUnlockAPIAsGlut(G);
  }
}

int MainCheckRedundantOpen(char *file)
{
  int result = false;
#ifndef _PYMOL_NOPY
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PBlock(G);
  result = PTruthCallStr(G->P_inst->cmd, "check_redundant_open", file);
  PUnblock(G);
#endif
  return result;
}

void MainSceneGetSize(int *width, int *height)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    SceneGetImageSize(G, width, height);
    PUnlockAPIAsGlut(G);
  }
}

int MainSceneCopy(int width, int height, int rowbytes, void *ptr)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  int result = false;
  if(PLockAPIAsGlut(G, true)) {
    result = SceneCopyExternal(G, width, height, rowbytes, (unsigned char *) ptr, 0x1);
    PUnlockAPIAsGlut(G);
  }
  return result;
}

PyObject *MainComplete(char *str)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PyObject *result = NULL;
  result = PyObject_CallFunction(G->P_inst->complete, "s", str);
  return (result);
}

#endif

#else


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#if 0
I think this code is dead ...
#ifdef _PYMOL_OSX
int *MacPyMOLReady = NULL;
CPyMOLOptions *MacPyMOLOption = NULL;
#endif
#endif

/* END PROPRIETARY CODE SEGMENT */

void MainFree(void);
void MainReshape(int width, int height);
static void MainDrawLocked(void);
static void MainDrag(int x, int y);

static CPyMOL *PyMOLInstance = NULL;    /* eliminate */

static char **myArgv, *myArgvv[2], myArgvvv[1024];
static int myArgc;

struct _CMain {
  int IdleMode;
  double IdleTime;
  int IdleCount;
  int Modifiers;
  int FinalInitCounter, FinalInitTrigger, FinalInitDone;
  int TheWindow;
  int WindowIsVisible;
  double ReshapeTime;
  double DrawAfter, DrawDelay;
  int DrawGovernorActive, DrawDeferred, DrawSignalled;
  int DrawnFlag, DeferReshapeDeferral;
  int MaximizeCheck;
  CPyMOLOptions *OwnedOptions;
  /* if Main owns a reference to a copy of the startup options that
     needs to be freed upon shutdown to fully scrub the heap */
};


/* global options */

void MainOnExit(void);

static void MainPushValidContext(PyMOLGlobals * G)
{
  PLockStatus(G);
  PyMOL_PushValidContext(G->PyMOL);
  PUnlockStatus(G);
}

static void MainPopValidContext(PyMOLGlobals * G)
{
  PLockStatus(G);
  PyMOL_PopValidContext(G->PyMOL);
  PUnlockStatus(G);
}

static void DrawBlueLine(PyMOLGlobals * G)
{
  if(G->Option->blue_line) {
    GLint i;
    unsigned long buffer;
    GLint window_width, window_height;

    window_width = G->Option->winX;
    window_height = G->Option->winY;

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDisable(GL_ALPHA_TEST);
    glDisable(GL_BLEND);
    for(i = 0; i < 6; i++)
      glDisable(GL_CLIP_PLANE0 + i);
    glDisable(GL_COLOR_LOGIC_OP);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DITHER);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LINE_STIPPLE);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_STENCIL_TEST);


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _PYMOL_OSX
    /*glDisable(GL_SHARED_TEXTURE_PALETTE_EXT); */
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_3D);
    glDisable(GL_TEXTURE_CUBE_MAP);
    glDisable(GL_TEXTURE_RECTANGLE_EXT);
    glDisable(GL_VERTEX_PROGRAM_ARB);
#endif

/* END PROPRIETARY CODE SEGMENT */

    for(buffer = GL_BACK_LEFT; buffer <= GL_BACK_RIGHT; buffer++) {
      GLint matrixMode;
      GLint vp[4];

      OrthoDrawBuffer(G, buffer);

      glGetIntegerv(GL_VIEWPORT, vp);
      glViewport(0, 0, window_width, window_height);

      glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();
      glScalef(2.0f / window_width, -2.0f / window_height, 1.0f);
      glTranslatef(-window_width / 2.0f, -window_height / 2.0f, 0.0f);

      /* draw sync lines */
      glColor3d(0.0f, 0.0f, 0.0f);
      glBegin(GL_LINES);        /* Draw a background line */
      glVertex3f(0.0F, window_height - 0.5F, 0.0F);
      glVertex3f((float) window_width, window_height - 0.5F, 0.0F);
      glEnd();
      glColor3d(0.0f, 0.0f, 1.0f);
      glBegin(GL_LINES);        /* Draw a line of the correct length (the cross
                                   over is about 40% across the screen from the left */
      glVertex3f(0.0f, window_height - 0.5f, 0.0f);
      if(buffer == GL_BACK_LEFT)
        glVertex3f(window_width * 0.30f, window_height - 0.5f, 0.0f);
      else
        glVertex3f(window_width * 0.80f, window_height - 0.5f, 0.0f);
      glEnd();

      glPopMatrix();
      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(matrixMode);

      glViewport(vp[0], vp[1], vp[2], vp[3]);
    }
    glPopAttrib();
  }
}


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#if 0
I think this code is dead ...
#ifdef _PYMOL_OSX

/* SPECIAL HOOKS FOR MacPyMOL */
int MainCheckRedundantOpen(char *file)
{
  int result = false;
#ifndef _PYMOL_NOPY
  PBlock(G);
  result = PTruthCallStr(G->P_inst->cmd_do, "check_redundant_open", file);
  PUnblock(G);
#endif
  return result;
}

void MainMovieCopyPrepare(int *width, int *height, int *length)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(true)) {
    MovieCopyPrepare(G, width, height, length);
    PUnlockAPIAsGlut(G);
  }
}

int MainMovieCopyFrame(int frame, int width, int height, int rowbytes, void *ptr)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  int result = false;
  if(PLockAPIAsGlut(G, true)) {
    result = MovieCopyFrame(G, frame, width, height, rowbytes, ptr);
    PUnlockAPIAsGlut(G);
  }
  return result;
}

void MainMovieCopyFinish(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    MovieCopyFinish(G);
    PUnlockAPIAsGlut(G);
  }
}

void MainSceneGetSize(int *width, int *height)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    SceneGetWidthHeight(G, width, height);
    PUnlockAPIAsGlut(G);
  }
}

int MainSceneCopy(int width, int height, int rowbytes, void *ptr)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  int result = false;
  if(PLockAPIAsGlut(G, true)) {
    result = SceneCopyExternal(G, width, height, rowbytes, (unsigned char *) ptr, 0);
    PUnlockAPIAsGlut(G);
  }
  return result;
}

void MainDoCommand(char *str1)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  if(PLockAPIAsGlut(G, true)) {
    if(str1[0] != '_') {        /* suppress internal call-backs */
      if(strncmp(str1, "cmd._", 5)) {
        OrthoAddOutput(G, "PyMOL>");
        OrthoAddOutput(G, str1);
        OrthoNewLine(G, NULL, true);
      }
      PDo(str1);
    } else if(str1[1] == ' ') { /* "_ command" suppresses echoing of command, but it is still logged */
      PDo(str1 + 2);
    } else {
      PDo(str1);
    }
    PUnlockAPIAsGlut(G);
  }
}

void MainRunCommand(char *str1)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  if(PLockAPIAsGlut(G, true)) {

    if(str1[0] != '_') {        /* suppress internal call-backs */
      if(strncmp(str1, "cmd._", 5)) {
        OrthoAddOutput(G, "PyMOL>");
        OrthoAddOutput(G, str1);
        OrthoNewLine(G, NULL, true);
        if(WordMatch(G, str1, "quit", true) == 0)       /* don't log quit */
          PLog(G, str1, cPLog_pml);
      }
      PParse(G, str1);
    } else if(str1[1] == ' ') { /* "_ command" suppresses echoing of command, but it is still logged */
      if(WordMatch(G, str1 + 2, "quit", true) >= 0)     /* don't log quit */
        PLog(G, str1 + 2, cPLog_pml);
      PParse(G, str1 + 2);
    } else {
      PParse(G, str1);
    }
    PUnlockAPIAsGlut(G);
  }
}

void MainFlushAsync(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  if(PLockAPIAsGlut(G, true)) {
    PFlush(G);
    PUnlockAPIAsGlut(G);
  }
}

void MainFlush(void)
{                               /* assumes GIL held */
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  MainPushValidContext(G);

  PFlush(G);

  MainPopValidContext(G);

}

void MainRunString(char *str)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PBlock(G);
  PLockStatus(G);
  MainPushValidContext(G);
  PRunStringModule(G, str);
  MainPopValidContext(G);
  PUnlockStatus(G);
  PUnblock(G);
}

PyObject *MainGetStringResult(char *str)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  PyObject *result;
  MainPushValidContext(G);
  result = PyRun_String(str, Py_eval_input, G->P_inst->dict, G->P_inst->dict);
  MainPopValidContext(G);
  return (result);
}
#endif

#endif

/* END PROPRIETARY CODE SEGMENT */


/*========================================================================*/

void MainOnExit(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  /* 
     here we enter not knowing anything about the current state...
     so, no graceful exit is possible -- in fact under Window's we'll
     crash unless we take the drastic way out 
   */
  if(!G->Terminating) {
    G->Terminating = true;
    printf(" PyMOL: abrupt program termination.\n");

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    TerminateProcess(GetCurrentProcess(), 0);   /* only way to avoid a crash */
#endif

/* END PROPRIETARY CODE SEGMENT */
    exit(EXIT_SUCCESS);
  }
}


/*========================================================================*/
int MainSavingUnderWhileIdle(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  CMain *I = G->Main;
  return (I && (I->IdleMode > 3));
}


/*========================================================================*/
void MainSetWindowVisibility(int mode)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  G->Option->window_visible = mode;

}


/*========================================================================*/
static void MainDrag(int x, int y)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  CMain *I = G->Main;

  if(PLockAPIAsGlut(G, false)) {

    y = G->Option->winY - y;

    PyMOL_Drag(PyMOLInstance, x, y, I->Modifiers);

    if(PyMOL_GetRedisplay(PyMOLInstance, true)) {
      if(G->HaveGUI) {
        p_glutPostRedisplay();
      }
    }
    I->IdleMode = 0;

    PUnlockAPIAsGlut(G);
  }
}


/*========================================================================*/
static void MainButton(int button, int state, int x, int y)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  int glMod;

  CMain *I = G->Main;

  glMod = p_glutGetModifiers();

  if(PLockAPIAsGlut(G, false)) {

    I->IdleMode = 0;            /* restore responsiveness */

    if(PyMOL_GetPassive(PyMOLInstance, (button < 3) ? true : false)) {
      MainDrag(x, y);
    } else {
      /* stay blocked here because Clicks->SexFrame->PParse */

      y = G->Option->winY - y;

      I->Modifiers = ((glMod & P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
        ((glMod & P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
        ((glMod & P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

      switch (button) {
      case P_GLUT_BUTTON_SCROLL_FORWARD:
      case P_GLUT_BUTTON_SCROLL_BACKWARD:
        x = G->Option->winX / 2;
        y = G->Option->winY / 2;        /* force into scene */
        break;
      }
      PyMOL_Button(PyMOLInstance, button, state, x, y, I->Modifiers);
    }

    PUnlockAPIAsGlut(G);
  }
}


/*========================================================================*/
static void MainPassive(int x, int y)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  CMain *I = G->Main;

#define PASSIVE_EDGE 20

  if(PyMOL_GetPassive(G->PyMOL, false)) {       /* a harmless race condition -- we don't want
                                                   to slow Python down buy locking on passive
                                                   mouse motion */

    if(PLockAPIAsGlut(G, false)) {
      if((y < -PASSIVE_EDGE) || (x < -PASSIVE_EDGE) ||
         (x > (G->Option->winX + PASSIVE_EDGE)) ||
         (y > (G->Option->winY + PASSIVE_EDGE))) {

        /* release passive drag if mouse leaves window... 
           [IF we continue to receive mouse events...] */

        y = G->Option->winY - y;

        PyMOL_Button(PyMOLInstance, P_GLUT_LEFT_BUTTON, P_GLUT_UP, x, y, I->Modifiers);

        PyMOL_GetPassive(G->PyMOL, true);       /* reset the flag */

      } else {

        y = G->Option->winY - y;

        PyMOL_Drag(PyMOLInstance, x, y, I->Modifiers);

      }

      if(PyMOL_GetRedisplay(PyMOLInstance, true)) {
        if(G->HaveGUI) {
          p_glutPostRedisplay();
        }
        I->IdleMode = 0;
      }

      PUnlockAPIAsGlut(G);
    }
  }

}


/*========================================================================*/
static void MainDrawLocked(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  CMain *I = G->Main;
  if(I->FinalInitTrigger) {

    I->FinalInitTrigger = false;

#ifndef _PYMOL_NOPY

    /* okay, on the way in, the API is locked but the interpreter is
       running async, so we need to block it */

    PBlock(G);

    if(PyErr_Occurred())
      PyErr_Print();

    /* next, we need to let PyMOL know we have a valid context,
       because some initializations, involve GL calls (such as going
       into full-screen or stereo mode) */

    if(G->HaveGUI)
      MainPushValidContext(G);

    /* restore working directory if asked to */
    PRunStringModule(G,
                     "if os.environ.has_key('PYMOL_WD'): os.chdir(os.environ['PYMOL_WD'])");

    if(PyErr_Occurred())
      PyErr_Print();


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _PYMOL_OSX
    PRunStringModule(G,
                     "if os.getcwd()[-23:]=='.app/Contents/Resources': os.chdir('../../..')");
#endif

/* END PROPRIETARY CODE SEGMENT */

    if(PyErr_Occurred())
      PyErr_Print();

    PXDecRef(PyObject_CallMethod(G->P_inst->obj, "launch_gui", "O", G->P_inst->obj));

    if(PyErr_Occurred())
      PyErr_Print();

    PXDecRef(PyObject_CallMethod
             (G->P_inst->obj, "adapt_to_hardware", "O", G->P_inst->obj));

    if(PyErr_Occurred())
      PyErr_Print();

    if(G->StereoCapable) {
      OrthoAddOutput(G,
                     " OpenGL quad-buffer stereo 3D detected and enabled.\n");;
      } else {
        if(G->LaunchStatus & cPyMOLGlobals_LaunchStatus_StereoFailed) {
          OrthoAddOutput(G,
                         "Error: The requested stereo 3D visualization mode is not available.\n");
        }
      }
    
    if(G->LaunchStatus & cPyMOLGlobals_LaunchStatus_MultisampleFailed) {
      OrthoAddOutput(G,
                     "Error: The requested multisampling mode is not available.\n");
    }

    if(G->Option->incentive_product) {  /* perform incentive product initialization (if any) */
      PyRun_SimpleString("try:\n   import ipymol\nexcept:\n   pass\n");
      if(PyErr_Occurred())
        PyErr_Print();
    }

    PXDecRef(PyObject_CallMethod(G->P_inst->obj, "exec_deferred", "O", G->P_inst->obj));

    if(PyErr_Occurred())
      PyErr_Print();

#ifdef _PYMOL_SHARP3D
    /*PParse("load $TUT/1hpv.pdb;hide;show sticks;show surface;set surface_color,white;set transparency,0.5;stereo on"); */
    /*PParse("stereo on");
       PParse("wizard demo,cartoon"); */
#endif

    if(G->HaveGUI)
      MainPopValidContext(G);

    PUnblock(G);

#endif
    I->FinalInitDone = true;
  }

  PyMOL_Draw(PyMOLInstance);

  if(G->HaveGUI) {
    if(Feedback(G, FB_OpenGL, FB_Debugging)) {
      PyMOLCheckOpenGLErr("During Rendering");
    }
  }

  if(PyMOL_GetSwap(G->PyMOL, true)) {
    if(!(int) SettingGet(G, cSetting_suspend_updates))
      if(G->HaveGUI) {
        DrawBlueLine(G);
        p_glutSwapBuffers();
      }
  }
}

static void MainDrawProgress(PyMOLGlobals * G)
{
  int progress[PYMOL_PROGRESS_SIZE];
  int update = false;
  PBlock(G);
  PLockStatus(G);
  update = PyMOL_GetProgress(G->PyMOL, progress, true);
  PUnlockStatus(G);
  PUnblock(G);

  /*
     printf("show progress %d %d %d %d %d %d\n",
     progress[0],progress[1],progress[2],
     progress[3],progress[4],progress[5]); */

  if(update &&
     (progress[PYMOL_PROGRESS_SLOW] ||
      progress[PYMOL_PROGRESS_MED] || progress[PYMOL_PROGRESS_FAST])) {

    int offset;
    int x = 0, y;
    float black[3] = { 0, 0, 0 };
    float white[3] = { 1, 1, 1 };
    GLint ViewPort[4];

#define cBusyWidth 240
#define cBusyHeight 60
#define cBusyMargin 10
#define cBusyBar 10
#define cBusySpacing 15

    glGetIntegerv(GL_VIEWPORT, ViewPort);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, ViewPort[2], 0, ViewPort[3], -100, 100);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0.33F, 0.33F, 0.0F);   /* this generates better 
                                           rasterization on macs */

    glDisable(GL_ALPHA_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);

    {
      int pass = 0;

      int draw_both = SceneMustDrawBoth(G);

      glClear(GL_DEPTH_BUFFER_BIT);
      while(1) {
        if(draw_both) {
          if(!pass)
            OrthoDrawBuffer(G, GL_FRONT_LEFT);
          else
            OrthoDrawBuffer(G, GL_FRONT_RIGHT);
        } else {
          OrthoDrawBuffer(G, GL_FRONT); /* draw into the front buffer */
        }

        glColor3fv(black);
        glBegin(GL_POLYGON);
        glVertex2i(0, ViewPort[3]);
        glVertex2i(cBusyWidth, ViewPort[3]);
        glVertex2i(cBusyWidth, ViewPort[3] - cBusyHeight);
        glVertex2i(0, ViewPort[3] - cBusyHeight);
        glVertex2i(0, ViewPort[3]);     /* needed on old buggy Mesa */
        glEnd();

        y = ViewPort[3] - cBusyMargin;

        glColor3fv(white);

        /* 
           c=I->BusyMessage;
           if(*c) {
           TextSetColor(G,white);
           TextSetPos2i(G,cBusyMargin,y-(cBusySpacing/2));
           TextDrawStr(G,c);
           y-=cBusySpacing;
           }
         */

        for(offset = 0; offset < PYMOL_PROGRESS_SIZE; offset += 2) {

          if(progress[offset + 1]) {

            glBegin(GL_LINE_LOOP);
            glVertex2i(cBusyMargin, y);
            glVertex2i(cBusyWidth - cBusyMargin, y);
            glVertex2i(cBusyWidth - cBusyMargin, y - cBusyBar);
            glVertex2i(cBusyMargin, y - cBusyBar);
            glVertex2i(cBusyMargin, y); /* needed on old buggy Mesa */
            glEnd();
            glColor3fv(white);
            glBegin(GL_POLYGON);
            glVertex2i(cBusyMargin, y);
            x =
              (progress[offset] * (cBusyWidth - 2 * cBusyMargin) / progress[offset + 1]) +
              cBusyMargin;
            glVertex2i(x, y);
            glVertex2i(x, y - cBusyBar);
            glVertex2i(cBusyMargin, y - cBusyBar);
            glVertex2i(cBusyMargin, y); /* needed on old buggy Mesa */
            glEnd();
            y -= cBusySpacing;
          }
        }
        if(!draw_both)
          break;
        if(pass > 1)
          break;
        pass++;
      }

      glFlush();
      glFinish();
      if(draw_both)
        OrthoDrawBuffer(G, GL_BACK_LEFT);
      else
        OrthoDrawBuffer(G, GL_BACK);
    }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

  }
}


/*========================================================================*/
static void MainDraw(void)
{
  int usec_sleep_after_draw = 0;
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  PRINTFD(G, FB_Main)
    " MainDraw: called.\n" ENDFD;
  if(PLockAPIAsGlut(G, false)) {
    CMain *I = G->Main;
    if(I->DrawGovernorActive &&
       (!I->DrawSignalled) && (UtilGetSeconds(G) < I->DrawAfter) && I->FinalInitDone) {
      I->DrawDeferred = true;
    } else {
      int skip = false;
      if(I->MaximizeCheck) {
        {
          /* is the window manager screwing us over??? */

          int height = p_glutGet(P_GLUT_SCREEN_HEIGHT);
          int width = p_glutGet(P_GLUT_SCREEN_WIDTH);
          int actual_x = p_glutGet(P_GLUT_WINDOW_X);
          int actual_y = p_glutGet(P_GLUT_WINDOW_Y);

          I->MaximizeCheck = false;

          if(actual_x != 0) {
            width -= 2 * actual_x;
            height -= actual_x;
          }
          if(actual_y != 0) {
            height -= actual_y;
          }
          p_glutPositionWindow(0, 0);
          p_glutReshapeWindow(width, height);
          skip = true;
        }
      }
      if((!skip) && (!I->DrawnFlag) && I->FinalInitDone) {
        if(I->DeferReshapeDeferral > 0)
          I->DeferReshapeDeferral--;
        else {
          double time_since_reshape = UtilGetSeconds(G) - I->ReshapeTime;
          if(time_since_reshape < 0.05) {
            /* defer screen updates while it's being actively resized */
            skip = true;
          }
        }
      }
      if(!skip) {
        MainDrawLocked();
        I->DrawnFlag = true;
        if(PyMOL_GetModalDraw(PyMOLInstance)) {
          usec_sleep_after_draw = 10000;        /* give other threads a chance to run */
        }
      } else {
        PyMOL_NeedRedisplay(PyMOLInstance);
      }

      I->DrawAfter = UtilGetSeconds(G) + I->DrawDelay;
      I->DrawSignalled = false;
      I->DrawDeferred = false;
    }
    PUnlockAPIAsGlut(G);
  } else {                      /* we're busy -- so try to display a progress indicator */
    MainDrawProgress(G);
  }
  PRINTFD(G, FB_Main)
    " MainDraw: completed.\n" ENDFD;
  if(usec_sleep_after_draw && G) {
    PSleepUnlocked(G, usec_sleep_after_draw);
  }
}


/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  CMain *I = G->Main;
  int glMod;

  glMod = p_glutGetModifiers();

  PRINTFB(G, FB_Main, FB_Blather)
    " MainKey: code:%d modifiers:0x%02x x:%d y:%d\n", k, glMod, x, y ENDFB(G);
  if(PLockAPIAsGlut(G, false)) {

    I->IdleMode = 0;            /* restore responsiveness */

    I->Modifiers = ((glMod & P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
      ((glMod & P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
      ((glMod & P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

    PyMOL_Key(PyMOLInstance, k, x, y, I->Modifiers);

    PUnlockAPIAsGlut(G);
  } else {
    if((k == 8) || (k == 127)) {        /* interrupt busy state (if possibele) */
      PBlock(G);
      PLockStatus(G);
      PyMOL_SetInterrupt(G->PyMOL, true);
      PUnlockStatus(G);
      PUnblock(G);
    }
  }
}


/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  CMain *I = G->Main;
  int glMod;

  PRINTFB(G, FB_Main, FB_Blather)
    " MainSpecial: %d %d %d\n", k, x, y ENDFB(G);
  glMod = p_glutGetModifiers();
  if(PLockAPIAsGlut(G, false)) {

    I->Modifiers = ((glMod & P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
      ((glMod & P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
      ((glMod & P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

    PyMOL_Special(PyMOLInstance, k, x, y, I->Modifiers);

    PUnlockAPIAsGlut(G);
  }
}


/* new window size or exposure */

/*========================================================================*/
void MainReshape(int width, int height)
{                               /* called by Glut */
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  if(G) {
    CMain *I = G->Main;

    I->ReshapeTime = (double) UtilGetSeconds(G);
    I->DrawnFlag = false;

    if(PLockAPIAsGlut(G, true)) {
      if(G->HaveGUI) {
        glViewport(0, 0, (GLint) width, (GLint) height);
        if((!PyMOLInstance) ||
           (width != OrthoGetWidth(G)) || (height != OrthoGetHeight(G))) {
          /* wipe the screen ASAP to prevent display of garbage... */

          int draw_both = G->StereoCapable &&
            ((SceneGetStereo(G) == 1) ||
             SettingGetGlobal_b(G, cSetting_stereo_double_pump_mono));

          glClearColor(0.0, 0.0, 0.0, 1.0);
          if(draw_both) {
            OrthoDrawBuffer(G, GL_FRONT_LEFT);
            glClear(GL_COLOR_BUFFER_BIT);
            OrthoDrawBuffer(G, GL_FRONT_RIGHT);
            glClear(GL_COLOR_BUFFER_BIT);
            OrthoDrawBuffer(G, GL_BACK_LEFT);
            glClear(GL_COLOR_BUFFER_BIT);
            OrthoDrawBuffer(G, GL_BACK_RIGHT);
            glClear(GL_COLOR_BUFFER_BIT);
          } else {
            OrthoDrawBuffer(G, GL_FRONT);
            glClear(GL_COLOR_BUFFER_BIT);
            OrthoDrawBuffer(G, GL_BACK);
            glClear(GL_COLOR_BUFFER_BIT);
          }
        }
        PyMOL_SwapBuffers(PyMOLInstance);
      }
    }
    if(PyMOLInstance)
      PyMOL_Reshape(PyMOLInstance, width, height, false);
    PUnlockAPIAsGlut(G);
  }
}


/*========================================================================*/
void MainDoReshape(int width, int height)
{                               /* called internally */
  int h, w;
  int internal_feedback;
  int force = false;
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  if(G) {
    /* if width is negative, force a reshape based on the current width */

    if(width < 0) {
      BlockGetSize(SceneGetBlock(G), &width, &h);
      if(SettingGetGlobal_b(G, cSetting_internal_gui))
        width += SettingGetGlobal_i(G, cSetting_internal_gui_width);
      force = true;
    } 

    /* if height is negative, force a reshape based on the current height */

    if(height < 0) {
      BlockGetSize(SceneGetBlock(G), &w, &height);
      internal_feedback = (int) SettingGet(G, cSetting_internal_feedback);
      if(internal_feedback)
        height += (internal_feedback - 1) * cOrthoLineHeight + cOrthoBottomSceneMargin;
      if(SettingGetGlobal_b(G, cSetting_seq_view)
         && !SettingGetGlobal_b(G, cSetting_seq_view_overlay))
        height += SeqGetHeight(G);
      height += MovieGetPanelHeight(G);
      force = true;
    } 

    /* if we have a GUI, for a reshape event */

    if(G->HaveGUI && G->ValidContext && width && height) {
      p_glutReshapeWindow(width, height);
      glViewport(0, 0, (GLint) width, (GLint) height);
    }

    if((!width)||(!height)) {
      int actual_width = width ? width : G->Option->winX;
      int actual_height = height ? height : G->Option->winY;

      PyMOL_Reshape(G->PyMOL, actual_width, actual_height, true);

    } else {
      PyMOL_Reshape(G->PyMOL, width, height, force);
      if(G->Main) {
        G->Main->DeferReshapeDeferral = 1;
      }

      /* do we need to become full-screen? */

      if(SettingGet(G, cSetting_full_screen) && G->HaveGUI && G->ValidContext) {
#ifndef __APPLE__
        p_glutFullScreen();
#else
        int height = p_glutGet(P_GLUT_SCREEN_HEIGHT);
        int width = p_glutGet(P_GLUT_SCREEN_WIDTH);
        height = height - 44;
        p_glutInitWindowPosition(0, 0);
        p_glutInitWindowSize(width, height);
#endif
      }
    }
  }

}


/*========================================================================*/
static void MainInit(PyMOLGlobals * G)
{

  CMain *I = (G->Main = Calloc(CMain, 1));
  /* Data structure is zeroed on start...no need for explicit zero inits */

  I->DeferReshapeDeferral = 1;

  PyMOL_Start(PyMOLInstance);
  PyMOL_ConfigureShadersGL(PyMOLInstance);

  PyMOL_SetSwapBuffersFn(PyMOLInstance, (PyMOLSwapBuffersFn *) p_glutSwapBuffers);
  I->ReshapeTime = (I->IdleTime = UtilGetSeconds(G));

  I->DrawGovernorActive = false;
  I->DrawSignalled = true;
  I->DrawDelay = 0.01;          /* 100 FPS max */
  I->DrawAfter = 0.0;
  I->DrawDeferred = false;
}


/*========================================================================*/
void MainFree(void)
{
  PyMOLGlobals *G = PyMOL_GetGlobals(PyMOLInstance);    /* temporary -- will change */

  CPyMOLOptions *owned_options = G->Main->OwnedOptions;

  int show_message = G->Option->show_splash && !G->Option->quiet;


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  int haveGUI = G->HaveGUI;
  int theWindow = G->Main->TheWindow;
#endif
#ifdef _PYMOL_OSX
  int game_mode = G->Option->game_mode;
  int haveGUI = G->HaveGUI;
  int theWindow = G->Main->TheWindow;
#endif

/* END PROPRIETARY CODE SEGMENT */

  PyMOL_PushValidContext(PyMOLInstance);
  PyMOL_Stop(PyMOLInstance);
  PyMOL_PopValidContext(PyMOLInstance);

  FreeP(G->Main);
  PyMOL_Free(PyMOLInstance);

  if(owned_options)
    PyMOLOptions_Free(owned_options);   /* clean up launch options if we're supposed to */

#ifndef OV_JX
  MemoryDebugDump();            /* this is a no-op unless memory debugging is enabled */
#endif

  if(show_message) {
    printf(" PyMOL: normal program termination.\n");
  }

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
  if(haveGUI)
    p_glutDestroyWindow(theWindow);
  TerminateProcess(GetCurrentProcess(), 0);     /* only way to avoid a crash */
#endif
#ifdef _PYMOL_OSX
  if(haveGUI) {
    if(game_mode) {
      p_glutLeaveGameMode();
      /* force a full-screen refresh to eliminate garbage on screen 
       * NOTE that we currently have to patch Apple's GLUT to make this work */
      p_glutInitWindowPosition(0, 0);
      p_glutInitWindowSize(640, 480);
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE);
      if(p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
        theWindow = p_glutCreateWindow("PyMOL Viewer");
        p_glutFullScreen();
        p_glutDestroyWindow(theWindow);
      }
    } else
      p_glutDestroyWindow(theWindow);
  }
#endif

/* END PROPRIETARY CODE SEGMENT */

}


/*========================================================================*/
void MainRefreshNow(void)
{                               /* should only be called by the master thread, with a locked API */
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  CMain *I = G->Main;
  if(PyMOL_GetSwap(G->PyMOL, true)) {
    if(G->HaveGUI) {
      DrawBlueLine(G);
      p_glutSwapBuffers();
    }
  }
  if(PyMOL_GetRedisplay(PyMOLInstance, true)) {
    if(G->HaveGUI)
      p_glutPostRedisplay();
    else
      MainDrawLocked();
    I->IdleMode = 0;
  }
}


/*========================================================================*/

#ifdef _PYMOL_SHARP3D
static int Sharp3DLastWindowX = -1000;
#endif

static void MainBusyIdle(void)
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  /* This is one of the few places in the program where we can be sure 
   * that we have the "glut" thread...glut doesn't seem to be completely
   * thread safe or rather thread consistent
   */

  CMain *I = G->Main;

  PRINTFD(G, FB_Main)
    " MainBusyIdle: called.\n" ENDFD;

#if 0
#ifdef  _PYMOL_SHARP3D
  /* keep the window on even coordinates to preserve L/R stereo... */
  {
    int x, y;
    x = glutGet(P_GLUT_WINDOW_X);
    if(x != Sharp3DLastWindowX) {
      Sharp3DLastWindowX = x;
      if(x & 0x1) {
        y = glutGet(P_GLUT_WINDOW_Y);
        glutPositionWindow(x - 1, y);
      }
    }
  }
#endif
#endif

  /* flush command and output queues */

  /*  PRINTFD(G,FB_Main)
     " MainBusyIdle: entered, IdleMode %d\n",
     I->IdleMode
     ENDFD; */

  if(PLockAPIAsGlut(G, false)) {

    PRINTFD(G, FB_Main)
      " MainBusyIdle: got lock.\n" ENDFD;

    /* change window visibility & refresh, if necessary */

    if(G->HaveGUI) {
      if(I->WindowIsVisible != G->Option->window_visible) {
        I->WindowIsVisible = G->Option->window_visible;
        if(I->WindowIsVisible) {
          p_glutShowWindow();
          OrthoDirty(G);
        } else {
          p_glutHideWindow();
        }
      }
    }

    PRINTFD(G, FB_Main)
      " MainBusyIdle: calling idle function.\n" ENDFD;

    if(PyMOL_Idle(PyMOLInstance)) {
      /* we did some work, so keep PyMOL responsive */
      I->IdleMode = 1;
    } else if(!I->IdleMode) {
      I->IdleMode = 1;
    } else if(I->IdleMode == 1) {
      I->IdleMode = 2;
      I->IdleTime = UtilGetSeconds(G);
    }
    PRINTFD(G, FB_Main)
      " MainBusyIdle: swap check.\n" ENDFD;

    if(PyMOL_GetSwap(G->PyMOL, true)) {
      if(G->HaveGUI) {
        DrawBlueLine(G);
        p_glutSwapBuffers();
      }
    }

    /* if the screen has become dirty, post a redisplay event, or if
       we're running without a GUI, then call the draw routine (if we */

    PRINTFD(G, FB_Main)
      " MainBusyIdle: redisplay.\n" ENDFD;

    if(PyMOL_GetRedisplay(PyMOLInstance, true)) {
      if(G->HaveGUI)
        p_glutPostRedisplay();
      else
        MainDrawLocked();
      if(I->IdleMode > 1) {
        I->IdleMode = 1;
      }
    }

    PRINTFD(G, FB_Main)
      " MainBusyIdle: redisplay.\n" ENDFD;

    /* the following code enables PyMOL to avoid busy-idling 
     * even though we're using GLUT! */

#ifndef _PYMOL_WX_GLUT
    /* however, don't spend any extra time sleeping in PYMOL if we're
       running under wxPython though... */

    switch (I->IdleMode) {
    case 2:                    /* avoid racing the CPU */
      if((UtilGetSeconds(G) - I->IdleTime) > (SettingGet(G, cSetting_idle_delay) / 5.0)) {
        I->IdleMode = 3;
        I->IdleTime = UtilGetSeconds(G);
      }
      break;
    case 3:
      if((UtilGetSeconds(G) - I->IdleTime) > (SettingGet(G, cSetting_idle_delay))) {
        I->IdleMode = 4;
        if(G->HaveGUI)
          if(SettingGet(G, cSetting_cache_display)) {
            p_glutPostRedisplay();      /* trigger caching of the current scene */
          }
        break;
      }
    }

#endif

    PRINTFD(G, FB_Main)
      " MainBusyIdle: unlocking.\n" ENDFD;

    {
      int control_idling = false;
      if(I->IdleMode == 1) {
        control_idling = ControlIdling(G);
      }
      PUnlockAPIAsGlut(G);

      switch (I->IdleMode) {
      case 4:
        PSleepUnlocked(G, SettingGetGlobal_i(G, cSetting_slow_idle));   /* slow idle - save CPU cycles */
        break;
      case 3:
        PSleepUnlocked(G, SettingGetGlobal_i(G, cSetting_fast_idle));   /* fast idle - more responsive */
        break;
      case 2:
        PSleepUnlocked(G, SettingGetGlobal_i(G, cSetting_no_idle));     /* give Tcl/Tk a chance to run */
        break;
      case 1:
        if(control_idling) {
          PSleepUnlocked(G, SettingGetGlobal_i(G, cSetting_no_idle));   /* give Tcl/Tk a chance to run */
        }
        break;
      default:
        break;
      }
    }

    /* run final initilization code for Python-based PyMOL implementations. */

#define FINAL_INIT_AT 10
    if(I->FinalInitCounter < FINAL_INIT_AT) {
      I->FinalInitCounter = I->FinalInitCounter + 1;
      if(I->FinalInitCounter == FINAL_INIT_AT) {
        I->FinalInitTrigger = true;
        PyMOL_NeedRedisplay(PyMOLInstance);
      }
    }

    /* when running in command-line mode, if we're not reading from
     * standard input and if we're not keeping the thread alive, then 
     * we can have no further input.  Therefore die. */

    if(!G->HaveGUI) {
      if(!(OrthoCommandWaiting(G) ||
           PyMOL_GetModalDraw(G->PyMOL) ||
           OrthoDeferredWaiting(G) || SettingGetGlobal_b(G, cSetting_keep_alive))) {
        if((!G->Option->keep_thread_alive) &&
           (!G->Option->read_stdin) && (I->FinalInitCounter >= FINAL_INIT_AT)) {
          I->IdleCount++;
          if(I->IdleCount == 10) {
            if(PLockAPIAsGlut(G, true)) {
              PParse(G, "_quit");
              PFlush(G);
              PUnlockAPIAsGlut(G);
            }
          }
        }
      } else {
        I->IdleCount = 0;
      }
    }

    {
      int max_ups = SettingGetGlobal_i(G, cSetting_max_ups);
      if(max_ups < 1) {
        I->DrawGovernorActive = false;
        if(I->DrawDeferred) {
          p_glutPostRedisplay();
        }
      } else {
        I->DrawDelay = 1.0 / max_ups;
        I->DrawGovernorActive = true;
        if(I->DrawDeferred) {
          if(UtilGetSeconds(G) > I->DrawAfter) {
            I->DrawSignalled = true;
          }
          if(I->DrawSignalled) {
            I->DrawDeferred = false;
            p_glutPostRedisplay();
          }
        }
      }
    }
  } else {
    PRINTFD(G, FB_Main)
      " MainBusyIdle: lock not obtained...\n" ENDFD;

    PSleepWhileBusy(G, 100000); /* 10 per second */
    if(G->HaveGUI) {
      PBlock(G);
      PLockStatus(G);
      if(PyMOL_GetProgressChanged(G->PyMOL, false))
        p_glutPostRedisplay();
      PUnlockStatus(G);
      PUnblock(G);
    }
  }

  PRINTFD(G, FB_Main)
    " MainBusyIdle: leaving... IdleMode %d\n", I->IdleMode ENDFD;

}

void MainRepositionWindowDefault(PyMOLGlobals * G)
{
#ifdef _PYMOL_FINK
  p_glutPositionWindow(G->Option->winPX, G->Option->winPY - 22);
  /* something is wrong with FinkGlut's window positioning... */
#else
  p_glutPositionWindow(G->Option->winPX, G->Option->winPY);
#endif

  p_glutReshapeWindow(G->Option->winX, G->Option->winY);
}

void MainSetWindowPosition(PyMOLGlobals * G, int x, int y)
{
#ifdef _PYMOL_FINK
  y -= 22;                      /* something is wrong with FinkGlut's window positioning... */
#endif
  p_glutPositionWindow(x, y);
}

void MainSetWindowSize(PyMOLGlobals * G, int w, int h)
{
  G->Main->DeferReshapeDeferral = 1;
  p_glutReshapeWindow(w, h);
}

void MainMaximizeWindow(PyMOLGlobals * G)
{
  int height = p_glutGet(P_GLUT_SCREEN_HEIGHT);
  int width = p_glutGet(P_GLUT_SCREEN_WIDTH);
  G->Main->DeferReshapeDeferral = 1;
  G->Main->MaximizeCheck = true;
  p_glutPositionWindow(0, 0);
  p_glutReshapeWindow(width, height);

}

void MainCheckWindowFit(PyMOLGlobals * G)
{
  CMain *I = G->Main;
  if(G && G->Main) {
    int height = p_glutGet(P_GLUT_SCREEN_HEIGHT);
    int width = p_glutGet(P_GLUT_SCREEN_WIDTH);
    int actual_x = p_glutGet(P_GLUT_WINDOW_X);
    int actual_y = p_glutGet(P_GLUT_WINDOW_Y);
    int actual_width = p_glutGet(P_GLUT_WINDOW_WIDTH);
    int actual_height = p_glutGet(P_GLUT_WINDOW_HEIGHT);
    int new_width = -1;
    int new_height = -1;

    I->DeferReshapeDeferral = 1;

    if((actual_x + actual_width) > width)
      new_width = (width - actual_x) - 5;       /* allow room for decoration */
    if((actual_y + actual_height) > height)
      new_height = (height - actual_y) - 5;     /* allow room for decoration */
    if((new_width > 0) || (new_height > 0)) {
      if(new_width < 0)
        new_width = actual_width;
      if(new_height < 0)
        new_height = actual_height;
      MainSetWindowSize(G, new_width, new_height);
    }
  }

}


/*========================================================================*/


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
BOOL WINAPI HandlerRoutine(DWORD dwCtrlType     /*  control signal type */
  )
{
  switch (dwCtrlType) {
  case CTRL_CLOSE_EVENT:
  case CTRL_BREAK_EVENT:
  case CTRL_C_EVENT:
    TerminateProcess(GetCurrentProcess(), 0);   /* only way to avoid a crash */

    break;
  }
  return 1;
}
#endif

#ifdef _PYMOL_SHARP3D
void sharp3d_prepare_context(void);
#endif

static void launch(CPyMOLOptions * options, int own_the_options)
{

  int multisample_mask = 0;
  int theWindow = 0;
  PyMOLGlobals *G = NULL;

  PyMOLInstance = PyMOL_NewWithOptions(options);
  G = PyMOL_GetGlobals(PyMOLInstance);

  /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _MACPYMOL_XCODE
  MacPyMOLOption = G->Option;
  MacPyMOLReady = &G->Ready;
#endif
  /* END PROPRIETARY CODE SEGMENT */

  if(G->Option->multisample)
    multisample_mask = P_GLUT_MULTISAMPLE;

  /* if were running GLUT, then we have the ability to increase the
   * size of the window in order to accomodate the context */

  if(G->Option->internal_gui && (!G->Option->game_mode))
    G->Option->winX += cOrthoRightSceneMargin;

  if(G->Option->internal_feedback && (!G->Option->game_mode))
    G->Option->winY +=
      (G->Option->internal_feedback - 1) * cOrthoLineHeight + cOrthoBottomSceneMargin;

  if(G->HaveGUI) {
#ifndef _PYMOL_OSX
    atexit(MainOnExit);         /* register callback to help prevent crashes
                                   when GLUT spontaneously kills us */
#endif

    /* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
    SetConsoleCtrlHandler(HandlerRoutine,       // address of handler function
                          true  // handler to add or remove
      );
    SetConsoleTitle(_T("PyMOL"));
#endif

#ifdef _PYMOL_SHARP3D
    sharp3d_prepare_context();
#endif

    p_glutInit(&myArgc, myArgv);

    {
      int display_mode_possible = false;
      if(G->Option->stereo_mode > 1)
        G->Option->force_stereo = 0;

      switch (G->Option->force_stereo) {
      case -1:                 /* force mono */
        G->StereoCapable = 0;
        break;
      case 0:                  /* default/autodetect */
#ifdef _PYMOL_SHARP3D
        G->Option->stereo_mode = cStereo_stencil_custom;
#endif
        switch (G->Option->stereo_mode) {
        case cStereo_default:
        case cStereo_quadbuffer:
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask |
                                P_GLUT_DOUBLE | P_GLUT_STEREO);
          display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
          if(multisample_mask && (!display_mode_possible)) {
            G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_MultisampleFailed;
            p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE |
                                  P_GLUT_STEREO);
            display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
          }
          if(display_mode_possible) {
            G->StereoCapable = 1;
          } else if(G->Option->stereo_mode == cStereo_quadbuffer) {
            G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_StereoFailed;
          }
          break;
        case cStereo_clone_dynamic:
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE | P_GLUT_ACCUM
                                | P_GLUT_STEREO);
          display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
          if(!display_mode_possible) {
            G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_StereoFailed;
            G->Option->stereo_mode = 0;
          } else {
            G->StereoCapable = 1;
          }
          break;
        case cStereo_dynamic:
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE |
                                P_GLUT_ACCUM);
          display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
          if(!display_mode_possible) {
            G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_StereoFailed;
            G->Option->stereo_mode = 0;
          }
          break;
        case cStereo_stencil_by_row:
        case cStereo_stencil_by_column:
        case cStereo_stencil_checkerboard:
        case cStereo_stencil_custom:
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE |
                                P_GLUT_STENCIL);
          display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
          if(!display_mode_possible) {
            G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_StereoFailed;
            G->Option->stereo_mode = 0;
          }
          break;
        case cStereo_anaglyph: /* nothing special required for anaglyph */
          G->StereoCapable = 1;
          break;
        default:               /* fall through */
          break;
        }
        break;
      case 1:                  /* force quad buffer stereo (if possible) */
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask |
                              P_GLUT_DOUBLE | P_GLUT_STEREO);
        display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
        if(multisample_mask && (!display_mode_possible)) {
          G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_MultisampleFailed;
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE |
                                P_GLUT_STEREO);
          display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
        }
        if(display_mode_possible) {
          G->StereoCapable = 1;
        } else {
          G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_StereoFailed;
        }
        break;
      }
      /* fallback behavior */

      if(!display_mode_possible) {
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask |
                              P_GLUT_DOUBLE);
        display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
        G->StereoCapable = 0;
      }

      if(multisample_mask && (!display_mode_possible)) {
        G->LaunchStatus |= cPyMOLGlobals_LaunchStatus_MultisampleFailed;
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE);
        display_mode_possible = p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE);
        G->StereoCapable = 0;
      }
    }

    if(!G->Option->game_mode) {
      if((G->Option->winPX > -10000) && (G->Option->winPY > -10000)) {
#ifndef _PYMOL_FINK
        p_glutInitWindowPosition(G->Option->winPX, G->Option->winPY);
#else
        p_glutInitWindowPosition(G->Option->winPX, G->Option->winPY - 22);      /* somethings wrong with FinkGlut's window positioning... */
#endif
      }

      p_glutInitWindowSize(G->Option->winX, G->Option->winY);

      if(G->Option->full_screen) {
        int height = p_glutGet(P_GLUT_SCREEN_HEIGHT);
        int width = p_glutGet(P_GLUT_SCREEN_WIDTH);
#ifdef __APPLE__
        height = height - 44;
#endif
        p_glutInitWindowPosition(0, 0);
        p_glutInitWindowSize(width, height);
      }

      theWindow = p_glutCreateWindow("PyMOL Viewer");

      if(G->Option->window_visible) {
        p_glutShowWindow();
      } else {
        p_glutHideWindow();
      }

    } else {
      char str[255];
      sprintf(str, "%dx%d:32@120", G->Option->winX, G->Option->winY);
      p_glutGameModeString(str);
      p_glutEnterGameMode();
    }
  }
#ifdef _PYMOL_SHARP3D
  /* Setup OpenGL */
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
#endif
  MainInit(G);
  if(own_the_options)
    G->Main->OwnedOptions = options;
  /* make sure we can clean up these options later if we've been asked to do so */
  {
    CMain *I = G->Main;

    I->TheWindow = theWindow;

    PInit(G, true);

#ifdef _PYMOL_SHARP3D
    SettingSetGlobal_f(G, cSetting_stereo_shift, 2.5);  /* increase strength */
#endif

    if(G->HaveGUI) {
      p_glutDisplayFunc(MainDraw);
      p_glutReshapeFunc(MainReshape);
      p_glutKeyboardFunc(MainKey);
      p_glutMouseFunc(MainButton);
      p_glutMotionFunc(MainDrag);
      p_glutPassiveMotionFunc(MainPassive);
      p_glutSpecialFunc(MainSpecial);
      p_glutIdleFunc(MainBusyIdle);
    }

    PUnblock(G);

    if(G->HaveGUI) {
      if(!I->WindowIsVisible)
        MainReshape(G->Option->winX, G->Option->winY);
      I->IdleMode = 3;
      p_glutMainLoop();         /* never returns with traditional GLUT implementation */
      PBlock(G);                /* if we've gotten here, then we're heading back to Python... */
    } else {
      SceneSetCardInfo(G, "none", "ray trace only", "none");
      if(G->Option->show_splash && !G->Option->quiet)
        printf(" Command mode. No graphics front end.\n");
      MainReshape(G->Option->winX, G->Option->winY);
      MainDraw();               /* for command line processing */
      while(1) {
        MainBusyIdle();
        MainDraw();
      }
    }
  }
}


/*========================================================================*/

#ifndef _PYMOL_MODULE
int main(int argc, char *argv[])
{
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  myArgc = argc;
  myArgv = argv;

  fflush(stdout);
  PSetupEmbedded(G, argc, argv);

#else
int was_main(void)
{
  myArgc = 1;
  strcpy(myArgvvv, "pymol");
  myArgvv[0] = myArgvvv;
  myArgvv[1] = NULL;
  myArgv = myArgvv;
#ifdef _DRI_WORKAROUND
  dlopen("libGL.so.1", RTLD_LAZY | RTLD_GLOBAL);
#endif

#endif

  {                             /* no matter how PyMOL was built, we always come through here... */

    CPyMOLOptions *options = PyMOLOptions_New();

    if(options) {

      PGetOptions(options);

      launch(options, true);
      /* this only returns when PyMOL is not running under GLUT */

    }
  }
  return 0;
}

#endif


/*========================================================================*/

PyObject *MainAsPyList(void);

PyObject *MainAsPyList(void)
{
#ifdef _PYMOL_NOPY
  return NULL;
#else
  PyMOLGlobals *G = SingletonPyMOLGlobals;

  PyObject *result = NULL;
  int width, height;
  result = PyList_New(2);
  BlockGetSize(SceneGetBlock(G), &width, &height);
  if(SettingGetGlobal_b(G, cSetting_seq_view)
     && !SettingGetGlobal_b(G, cSetting_seq_view_overlay))
    height += SeqGetHeight(G);
  PyList_SetItem(result, 0, PyInt_FromLong(width));
  PyList_SetItem(result, 1, PyInt_FromLong(height));
  return (PConvAutoNone(result));
#endif
}

int MainFromPyList(PyObject * list);
int MainFromPyList(PyObject * list)
{
#ifdef _PYMOL_NOPY
  return 0;
#else

  int ok = true;
  int win_x, win_y;
  int ll = 0;
  PyMOLGlobals *G = SingletonPyMOLGlobals;
  OrthoLineType buffer;

  if(ok)
    ok = (list != NULL);
  if(ok)
    ok = PyList_Check(list);
  if(ok)
    ll = PyList_Size(list);
  if(ok && (ll >= 2)) {
    if(!G->Option->presentation) {
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 0), &win_x);
      if(ok)
        ok = PConvPyIntToInt(PyList_GetItem(list, 1), &win_y);
      /* BlockGetSize(SceneGetBlock(G),&win_x,&win_y); * so how did this get into 0.98beta29?  */
      if(ok) {

        sprintf(buffer, "viewport %d, %d", win_x, win_y);
        PParse(G, buffer);
      }
    }
  }
  return (ok);
#endif
}
