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
#include"os_python.h"

#ifdef WIN32
#include<signal.h>
#endif

#include"os_std.h"
#include"os_gl.h"

#ifdef _PYMOL_MODULE
#ifdef _DRI_WORKAROUND
#include <dlfcn.h>
#endif
#endif

#include"MemoryDebug.h"
#include"MemoryCache.h"
#include"Err.h"
#include"Util.h"
#include"Selector.h"
#include"Color.h"
#include"Ortho.h"
#include"Scene.h"
#include"PyMOLObject.h"
#include"P.h"
#include"Executive.h"
#include"Word.h"
#include"RepMesh.h"
#include"ObjectMolecule.h"
#include"Control.h"
#include"Sphere.h"
#include"Setting.h"
#include"Ray.h"
#include"Util.h"
#include"Movie.h"
#include"main.h"
#include"Editor.h"
#include"SculptCache.h"
#include"Isosurf.h"
#include"Tetsurf.h"
#include"PConv.h"
#include"VFont.h"
#include"Wizard.h"
#include"Text.h"
#include"Character.h"
#include"Seq.h"
#include"Seeker.h"

#include"PyMOLGlobals.h"
#include"ClassPyMOL.h"

void MainFree(void);
void MainTest(void);
void MainBusyIdle(void);
static void MainInit(PyMOLGlobals *G);
void MainReshape(int width, int height);
static void MainDrawLocked(void);
static void MainDrag(int x,int y);

static  ClassPyMOL *PyMOLInstance; /* eliminate */

GLuint obj;

/* GLUT window size */

static int WinX = 640; /* these values are overridden! */
static int WinY = 480;

/* GLUT window initial position */

static int WinPX = 0;    /* these values are overridden! */
static int WinPY = 175;

static GLint Modifiers = 0;

static char **myArgv,*myArgvv[2],myArgvvv[1024];
static int myArgc;

static int FinalInitFlag=1;

int TheWindow;

typedef struct {
  int DirtyFlag;
  int IdleMode;
  int SwapFlag;
  double IdleTime;
  int IdleCount;
  int ReshapeFlag;
  int DragDirtyFlag;
  int DragPassive;
} CMain;

static CMain Main;
int PyMOLReady = false;
int PyMOLTerminating = false;

/* global options */

int ExternalGUI=1;
int PMGUI = true;
int StereoCapable=false;
int Security = true;
int ForceStereo = 0; 
  /* 1 = force stereo (if possible); -1 = force mono; 0 = autodetect */
int GameMode = false;
int BlueLine = false;

static int InternalGUI = true;
static int InternalFeedback = true;
static int WindowIsVisible = false;

int ShowSplash=true;
int PyMOLRegisterSigIntHandler = true;
int PyMOLMultisample = 0;

static PyMOLOptionRec PyMOLOptionGlobal = {
  true, /* pmgui */
  true, /* internal_gui*/
  true, /* show_splash */
  1,   /* internal_feedback */
  true, /* security */
  false, /* game mode */
  0, /* force_stereo */
  640, /* winX */
  480, /* winY */
  false, /* blue_line */
  0, /* winPX */
  175, /* winPY */
  true, /* external_gui */
  true, /* siginthand */
  false, /* reuse helper */
  false, /* auto reinitialize */
  false, /* keep thread alive */
  false, /* quiet */
  false, /* incentive product */
  "", /* after_load_script */
  0, /* multisample */
  1, /* window_visible */
  0, /* read_stdin */
};

struct PyMOLOptionRec *PyMOLOption = &PyMOLOptionGlobal;

void launch(void);

void MainOnExit(void);

void MainDragDirty(void)
{
  CMain *I = &Main;
  I->DragDirtyFlag = 1;
}

static void DrawBlueLine(void)
{
  GLint i;
  unsigned long buffer;
  GLint window_width, window_height;
  window_width=WinX;
  window_height=WinY;

  if(BlueLine) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
 
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_BLEND);
    for(i = 0; i < 6; i++) glDisable(GL_CLIP_PLANE0 + i);
    glDisable(GL_COLOR_LOGIC_OP);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_DITHER);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LINE_STIPPLE);
    glDisable(GL_SCISSOR_TEST);
#ifdef _PYMOL_OSX
    glDisable(GL_SHARED_TEXTURE_PALETTE_EXT);
#endif
    glDisable(GL_STENCIL_TEST);
#ifdef _PYMOL_OSX
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_3D);
    glDisable(GL_TEXTURE_CUBE_MAP);
    glDisable(GL_TEXTURE_RECTANGLE_EXT);
    glDisable(GL_VERTEX_PROGRAM_ARB);
#endif
  
    for(buffer = GL_BACK_LEFT; buffer <= GL_BACK_RIGHT; buffer++) {
      GLint matrixMode;
      GLint vp[4];
  
      glDrawBuffer(buffer);
  
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
 
      /* draw sync lines*/
      glColor3d(0.0f, 0.0f, 0.0f);
      glBegin(GL_LINES); /* Draw a background line*/
      glVertex3f(0.0F, window_height - 0.5F, 0.0F);
      glVertex3f((float)window_width, window_height - 0.5F, 0.0F);
      glEnd();
      glColor3d(0.0f, 0.0f, 1.0f);
      glBegin(GL_LINES); /* Draw a line of the correct length (the cross
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

#ifdef _PYMOL_OSX

void MainMovieCopyPrepare(int *width,int *height,int *length)
{
  PLockAPIAsGlut();
  MovieCopyPrepare(TempPyMOLGlobals,width,height,length);
  PUnlockAPIAsGlut();
}
int MainMovieCopyFrame(int frame,int width,int height,int rowbytes,void *ptr)
{
  int result;
  PLockAPIAsGlut();
  result = MovieCopyFrame(TempPyMOLGlobals,frame,width,height,rowbytes,ptr);
  PUnlockAPIAsGlut();
  return result;
}
void MainMovieCopyFinish(void)
{
  PLockAPIAsGlut();
  MovieCopyFinish(TempPyMOLGlobals);
  PUnlockAPIAsGlut();
}
void MainSceneGetSize(int *width,int *height)
{
  PLockAPIAsGlut();
  SceneGetWidthHeight(TempPyMOLGlobals,width,height);
  PUnlockAPIAsGlut();
}
int MainSceneCopy(int width,int height,int rowbytes,void *ptr)
{ 
  int result;
  PLockAPIAsGlut();
  result = SceneCopyExternal(TempPyMOLGlobals,width, height,rowbytes,(unsigned char *)ptr);
  PUnlockAPIAsGlut();
  return result;
}
/*========================================================================*/
void MainDoCommand(char *str1)
{
  PLockAPIAsGlut();
  if(str1[0]!='_') { /* suppress internal call-backs */
    if(strncmp(str1,"cmd._",5)) {
      OrthoAddOutput(TempPyMOLGlobals,"PyMOL>");
      OrthoAddOutput(TempPyMOLGlobals,str1);
      OrthoNewLine(TempPyMOLGlobals,NULL,true);
    }
    PDo(str1);
  } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
    PDo(str1+2);    
  } else { 
    PDo(str1);
  }
  PUnlockAPIAsGlut();
}
/*========================================================================*/
void MainRunCommand(char *str1)
{
  PLockAPIAsGlut();
  if(str1[0]!='_') { /* suppress internal call-backs */
    if(strncmp(str1,"cmd._",5)) {
      OrthoAddOutput(TempPyMOLGlobals,"PyMOL>");
      OrthoAddOutput(TempPyMOLGlobals,str1);
      OrthoNewLine(TempPyMOLGlobals,NULL,true);
      if(WordMatch(TempPyMOLGlobals,str1,"quit",true)==0) /* don't log quit */
        PLog(str1,cPLog_pml);
    }
    PParse(str1);
  } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
    if(WordMatch(TempPyMOLGlobals,str1+2,"quit",true)>=0) /* don't log quit */
      PLog(str1+2,cPLog_pml);
    PParse(str1+2);    
  } else { 
    PParse(str1);
  }
  PUnlockAPIAsGlut();
}

/*========================================================================*/
void MainFlushAsync(void)
{
  PLockAPIAsGlut();
  PFlush();
  PUnlockAPIAsGlut();
}
/*========================================================================*/
void MainFlush(void)
{
  PFlush();
}
/*========================================================================*/
void MainRunString(char *str)
{
  PBlock();
  PRunString(str);
  PUnblock();
}
/*========================================================================*/
PyObject *MainGetStringResult(char *str)
{
  PyObject *result;
  result = PyRun_String(str,Py_eval_input,P_globals,P_globals);
  return(result);
}

#endif

/*========================================================================*/

void MainOnExit(void)
{ /* 
     here we enter not knowing anything about the current state...
     so, no graceful exit is possible -- in fact under Window's we'll
     crash unless we take the drastic way out 
  */
  if(!PyMOLTerminating) {
    PyMOLTerminating=true;
	printf(" PyMOL: abrupt program termination.\n");
#ifdef WIN32
	TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
#endif
    exit(EXIT_SUCCESS);
  }
}
/*========================================================================*/
int MainSavingUnderWhileIdle(void)
{
  CMain *I = &Main;
  return(I->IdleMode==2);
}
/*========================================================================*/
void MainResetIdle(void)
{
  CMain *I = &Main;
  I->IdleMode = 0;
}
/*========================================================================*/
void MainDirty(void)
{
  CMain *I = &Main;
  /*  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainDirty: called.\n"
    ENDFD; */
  I->DirtyFlag=true;
  I->IdleMode = 0;
}
/*========================================================================*/
void MainSwapBuffers(void)
{
  CMain *I = &Main;
  I->SwapFlag=true;
}
/*========================================================================*/
void MainSetPassiveDrag(int onOrOff)
{
  CMain *I = &Main;
  I->DragPassive = onOrOff;
}
/*========================================================================*/
void MainSetWindowVisibility(int mode)
{
  PyMOLOption->window_visible = mode;
}
/*========================================================================*/
void MainTest(void)
{
}
/*========================================================================*/
static void MainButton(int button,int state,int x,int y)
{
  static int glMod;  
  CMain *I = &Main;

  glMod = p_glutGetModifiers();

  PLockAPIAsGlut();

  if(I->DragPassive) {
    I->DragPassive = false;
    MainDrag(x,y);
  } else {
    /* stay blocked here because Clicks->SexFrame->PParse */
    
    y=WinY-y;
    
    Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
      ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
      ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);
    
    switch(button) {
    case P_GLUT_BUTTON_SCROLL_FORWARD:
    case P_GLUT_BUTTON_SCROLL_BACKWARD:
      x=WinX/2;y=WinY/2; /* force into scene */
      break;
    }
    if(!OrthoButton(TempPyMOLGlobals,button,state,x,y,Modifiers))
      {
      }
  }

  PUnlockAPIAsGlut();

}
/*========================================================================*/
static void MainDrag(int x,int y)
{
  CMain *I = &Main;
  
  PLockAPIAsGlut();
  
  y=WinY-y;
  if(!OrthoDrag(TempPyMOLGlobals,x,y,Modifiers))
    {
    }
  
  if(I->DirtyFlag)
    if(PMGUI) {
      p_glutPostRedisplay();
      I->DirtyFlag=false;
    }
  
  PUnlockAPIAsGlut();

}
/*========================================================================*/
static void MainPassive(int x,int y)
{
  CMain *I = &Main;

  if(I->DragPassive) { /* a harmless race condition -- we don't want
                          to slow Python down buy locking on passive
                          mouse motion */

    PLockAPIAsGlut();

    if((y<0)||(x<0)||(x>WinX)||(y>WinY)) {       /* release passive drag if mouse leaves window... */

      y=WinY-y;
     

      OrthoButton(TempPyMOLGlobals,P_GLUT_LEFT_BUTTON, P_GLUT_UP,x,y,Modifiers);
      I->DragPassive = false;

    } else {
      
      y=WinY-y;
      
      if(!OrthoDrag(TempPyMOLGlobals,x,y,Modifiers))
        {
        }
      
      if(I->DirtyFlag)
        if(PMGUI) {
          p_glutPostRedisplay();
          I->DirtyFlag=false;
        }
    }
    
    PUnlockAPIAsGlut();
  }
  
}

/*========================================================================*/
static void MainDrawLocked(void)
{
  CMain *I = &Main;
  
  
  if(I->DirtyFlag) {
    I->DirtyFlag=false;
  }
  
  OrthoBusyPrime(TempPyMOLGlobals);
  ExecutiveDrawNow(TempPyMOLGlobals);

  if(PMGUI) {
    if(Feedback(TempPyMOLGlobals,FB_OpenGL,FB_Debugging)) {
      PyMOLCheckOpenGLErr("During Rendering");
    }
  }

  if(I->SwapFlag)
    {
      if(!(int)SettingGet(TempPyMOLGlobals,cSetting_suspend_updates))
        if(PMGUI) {
          DrawBlueLine();
          p_glutSwapBuffers();
        }
      I->SwapFlag=false;
    }
  /* don't claim to be ready until we've drawn at least once */
  PyMOLReady = true; 
}
/*========================================================================*/
static void MainDraw(void)
{
  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainDraw: called.\n"
    ENDFD;
  PLockAPIAsGlut();

  MainDrawLocked();
  PUnlockAPIAsGlut();
  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainDraw: completed.\n"
    ENDFD;
}
/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  /*  CMain *I = &Main;*/
  int glMod;

  glMod = p_glutGetModifiers();

  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainKey: %d %d %d\n",k,x,y
    ENDFD;
  PLockAPIAsGlut();

  Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!WizardDoKey(TempPyMOLGlobals,k,x,y,Modifiers))
    OrthoKey(TempPyMOLGlobals,k,x,y,Modifiers);

  PUnlockAPIAsGlut();
  
}

/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  char buffer[255];
  int grabbed = false;
  static int glMod;  

  glMod = p_glutGetModifiers();
  PLockAPIAsGlut();

  Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!grabbed)
    grabbed = WizardDoKey(TempPyMOLGlobals,(unsigned char)k,x,y,Modifiers);

  switch(k) {
    case P_GLUT_KEY_UP:
    case P_GLUT_KEY_DOWN:
      grabbed=1;
      OrthoSpecial(TempPyMOLGlobals,k,x,y,Modifiers);
      break;
    case P_GLUT_KEY_LEFT:
    case P_GLUT_KEY_RIGHT:      
      if(OrthoArrowsGrabbed(TempPyMOLGlobals)) {
        grabbed=1;
        OrthoSpecial(TempPyMOLGlobals,k,x,y,Modifiers);
      }
      break;
  }

  if(!grabbed) {
    sprintf(buffer,"_special %d,%d,%d,%d",k,x,y,Modifiers);
    PLog(buffer,cPLog_pml);
    PParse(buffer);
    PFlush();
  }
  PUnlockAPIAsGlut();
}

/* new window size or exposure */
/*========================================================================*/
void MainReshape(int width, int height) /* called by Glut */
{
  float h;

  WinX = width;
  WinY = height;

  h = ((float)height)/width;
  if(PMGUI) glViewport(0, 0, (GLint) width, (GLint) height);
  
  OrthoReshape(TempPyMOLGlobals,width,height,false);

}
/*========================================================================*/

PyObject *MainAsPyList(void) 
{
  PyObject *result=NULL;
  int width,height;
  result = PyList_New(2);
  BlockGetSize(SceneGetBlock(TempPyMOLGlobals),&width,&height);
  if(SettingGetGlobal_b(TempPyMOLGlobals,cSetting_seq_view)&&!SettingGetGlobal_b(TempPyMOLGlobals,cSetting_seq_view_overlay))
    height+=SeqGetHeight(TempPyMOLGlobals);
  PyList_SetItem(result,0,PyInt_FromLong(width));
  PyList_SetItem(result,1,PyInt_FromLong(height));
  return(PConvAutoNone(result));
}

int MainFromPyList(PyObject *list)
{
  int ok=true;
  int win_x,win_y;
  int ll=0;
  OrthoLineType buffer;
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) ll=PyList_Size(list);
  if(ok&&(ll>=2)) {
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,0),&win_x);
    if(ok) ok = PConvPyIntToInt(PyList_GetItem(list,1),&win_y);
    if(ok) {
      sprintf(buffer,"viewport %d, %d",win_x,win_y);
      PParse(buffer);
    }
  }
  return(ok);
}
/*========================================================================*/
void MainDoReshape(int width, int height) /* called internally */
{
  int h,w;
  int internal_feedback;
  int force =false;

  if(width<0) {
    BlockGetSize(SceneGetBlock(TempPyMOLGlobals),&width,&h);
    if(SettingGetGlobal_b(TempPyMOLGlobals,cSetting_internal_gui))
      width+=SettingGetGlobal_i(TempPyMOLGlobals,cSetting_internal_gui_width);
    force = true;
  }

  if(height<0) { 
    BlockGetSize(SceneGetBlock(TempPyMOLGlobals),&w,&height);
    internal_feedback = (int)SettingGet(TempPyMOLGlobals,cSetting_internal_feedback);
    if(internal_feedback)
      height+=(internal_feedback-1)*cOrthoLineHeight+cOrthoBottomSceneMargin;
    if(SettingGetGlobal_b(TempPyMOLGlobals,cSetting_seq_view)&&!SettingGetGlobal_b(TempPyMOLGlobals,cSetting_seq_view_overlay))
      height+=SeqGetHeight(TempPyMOLGlobals);
    force = true;
  }
  if(PMGUI) {
    p_glutReshapeWindow(width,height);
    glViewport(0, 0, (GLint) width, (GLint) height);
  }
  OrthoReshape(TempPyMOLGlobals,width,height,force);
  if(SettingGet(TempPyMOLGlobals,cSetting_full_screen))
    p_glutFullScreen();

}
/*========================================================================*/
static void MainInit(PyMOLGlobals *G) 
{
  /*  GLfloat one[4] = { 1,1,1,1 }; 
      GLfloat low[4] = { 0.20,0.20,0.20,1 };*/

  CMain *I = &Main;

  I->DirtyFlag=true;
  I->IdleMode=2;
  I->IdleCount = 0;
  I->ReshapeFlag = false;
  I->DragDirtyFlag=0;
  I->DragPassive = false;

  if(PMGUI) {

    /* get us into a well defined GL state */
    
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
  }

  MemoryCacheInit(G);
  FeedbackInit(G,PyMOLOption->quiet);
  WordInit(G);
  UtilInit(G);

  I->IdleTime=(float)UtilGetSeconds(G);



  ColorInit(G);
  CGORendererInit(G);
  SettingInitGlobal(G,true,true);  
  SettingSet(G,cSetting_internal_gui,(float)InternalGUI);
  SettingSet(G,cSetting_internal_feedback,(float)InternalFeedback);
  TextInit(G);
  CharacterInit(G);
  SphereInit(G);
  OrthoInit(G,ShowSplash);
  WizardInit(G); /* must come after ortho */
  MovieInit(G);
  SceneInit(G);
  SelectorInit(G);
  SeqInit(G);
  SeekerInit(G);
  ButModeInit(G);
  ControlInit(G);
  AtomInfoInit(G);
  SculptCacheInit(G);
  VFontInit(G);
  ExecutiveInit(G);
  IsosurfInit(G);
  TetsurfInit(G);
  EditorInit(G);
}




/*========================================================================*/
void MainFree(void)
{
  PyMOLGlobals *G = ClassPyMOLGetGlobals(PyMOLInstance); /* temporary -- will change */

  PyMOLTerminating=true;
  TetsurfFree(G);
  IsosurfFree(G);
  WizardFree(G);
  SceneCleanupStereo(G);
  EditorFree(G);
  ExecutiveFree(G);
  VFontFree(G);
  SculptCacheFree(G);
  AtomInfoFree(G);
  ButModeFree(G);
  ControlFree(G);
  SeekerFree(G);
  SeqFree(G);
  SelectorFree(G);
  SceneFree(G);
  MovieFree(G);
  OrthoFree(G);
  SettingFreeGlobal(G);
  CharacterFree(G);
  TextFree(G);
  SphereFree(G);
  PFree();
  CGORendererFree(G);
  ColorFree(G);
  UtilFree(G);
  WordFree(G);
  FeedbackFree(G);
  MemoryCacheDone(G);
                     
  ClassPyMOLFree(PyMOLInstance);

  if(ShowSplash) {
    MemoryDebugDump();
    printf(" PyMOL: normal program termination.\n");
  }
  

#ifdef WIN32
  if(PMGUI) p_glutDestroyWindow(TheWindow);
  TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
#endif
#ifdef _PYMOL_OSX
  if(PMGUI) {
    if(GameMode) {
      p_glutLeaveGameMode();
      /* force a full-screen refresh to eliminate garbage on screen */
      p_glutInitWindowPosition(0,0);
      p_glutInitWindowSize(640,480);
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );            
      if(p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
        TheWindow = p_glutCreateWindow("PyMOL Viewer");
        p_glutFullScreen();
        p_glutDestroyWindow(TheWindow);
      }
    } else 
      p_glutDestroyWindow(TheWindow);
  }
#endif

}
/*========================================================================*/
void MainRefreshNow(void) 
{ /* should only be called by the master thread, with a locked API */

  CMain *I = &Main;
  if(I->SwapFlag)
    {
      if(PMGUI) {
        DrawBlueLine();
        p_glutSwapBuffers();
      }
      I->SwapFlag=false;
    }
  if(I->DirtyFlag)
    {
      if(PMGUI) 
        p_glutPostRedisplay();
      else
        MainDrawLocked();
      I->DirtyFlag=false;
    }
}
/*========================================================================*/

#ifdef _PYMOL_SHARP3D
static int Sharp3DLastWindowX = -1000;
#endif

void MainBusyIdle(void) 
{
  /* This is one of the few places in the program where we can be sure 
	* that we have the "glut" thread...glut doesn't seem to be completely
	* thread safe or rather thread consistent
   */

  CMain *I = &Main;


  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainBusyIdle: called.\n"
    ENDFD;

#ifdef _PYMOL_SHARP3D
// keep the window on even coordinates to preserve L/R stereo...
 {
   int x,y;
   x = glutGet(GLUT_WINDOW_X);
   if(x!=Sharp3DLastWindowX) {
     Sharp3DLastWindowX=x;
     if(x&0x1) {
       y = glutGet(GLUT_WINDOW_Y);
       glutPositionWindow(x-1,y);
     }
   }
 }
#endif


  /* flush command and output queues */
  
  /*  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainBusyIdle: entered, IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;*/
  PLockAPIAsGlut();

  if(PMGUI) {
    if(WindowIsVisible!=PyMOLOption->window_visible) {
      WindowIsVisible = PyMOLOption->window_visible;
      if(WindowIsVisible) {
        p_glutShowWindow();
        OrthoDirty(TempPyMOLGlobals);
      } else {
        p_glutHideWindow();
      }
    }
  }

  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainBusyIdle: got lock.\n"
    ENDFD;

  if(I->DragDirtyFlag==1) {
    I->DragDirtyFlag = 0;
    OrthoFakeDrag(TempPyMOLGlobals);
  }

  if(ControlIdling(TempPyMOLGlobals)) {
    ExecutiveSculptIterateAll(TempPyMOLGlobals);
	 SceneIdle(TempPyMOLGlobals); 
	 I->IdleMode=0;
  } else {
    if(!I->IdleMode) {
      I->IdleTime=UtilGetSeconds(TempPyMOLGlobals);
      I->IdleMode=1;
    }
  }

  if(SceneRovingCheckDirty(TempPyMOLGlobals)) {
    SceneRovingUpdate(TempPyMOLGlobals);
    I->IdleMode = 0;
  }
  PFlush();

  if(I->SwapFlag) {
    if(PMGUI) {
      DrawBlueLine();
      p_glutSwapBuffers();
    }
    I->SwapFlag=false;
  }
  if(I->DirtyFlag&&
     ((!SettingGet_b(TempPyMOLGlobals,NULL,NULL,cSetting_defer_updates)||
       ControlIdling(TempPyMOLGlobals)))
     ) {
    if(PMGUI) 
      p_glutPostRedisplay();
    else
      MainDrawLocked();
    I->DirtyFlag=false;
  }

  if(I->IdleMode) { /* avoid racing the CPU */
    if(I->IdleMode==1) {
      if(UtilGetSeconds(TempPyMOLGlobals)-I->IdleTime>SettingGet(TempPyMOLGlobals,cSetting_idle_delay)) { 
        I->IdleMode=2;
        if(PMGUI)
          if(SettingGet(TempPyMOLGlobals,cSetting_cache_display))
             p_glutPostRedisplay(); /* trigger caching of the current scene */
      }
    }
    if(I->IdleMode==1)
      PSleep((int)SettingGet(TempPyMOLGlobals,cSetting_fast_idle)); /* fast idle - more responsive */
    else
      PSleep((int)SettingGet(TempPyMOLGlobals,cSetting_slow_idle)); /* slow idle - save CPU cycles */
  } else {
    PSleep((int)SettingGet(TempPyMOLGlobals,cSetting_no_idle)); /* give Tcl/Tk a chance to run */
  }
  PUnlockAPIAsGlut();

  if(FinalInitFlag)
	 {
      FinalInitFlag=FinalInitFlag+1;
      if(FinalInitFlag>=10) {
        FinalInitFlag=0;
        PBlock();
        /* restore working directory if asked to */
        PRunString("if os.environ.has_key('PYMOL_WD'): os.chdir(os.environ['PYMOL_WD'])");

        #ifdef _PYMOL_OSX
        PRunString("if os.getcwd()[-23:]=='.app/Contents/Resources': os.chdir('../../..')");
        #endif

        PRunString("launch_gui()");

        /*#ifndef _PYMOL_WX_GLUT
         PRunString("launch_gui()");
        #endif*/

        PRunString("adapt_to_hardware()");

        if(PyMOLOption->incentive_product) { /* perform incentive product initialization (if any) */
          PyRun_SimpleString("try:\n   import ipymol\nexcept:\n   pass\n");
        }

        PRunString("exec_deferred()");
#ifdef _PYMOL_SHARP3D
        //PParse("load $TUT/1hpv.pdb;hide;show sticks;show surface;set surface_color,white;set transparency,0.5;stereo on");
        PParse("stereo on");
        //PParse("wizard demo,cartoon");
#endif
        PUnblock();
      }
    }
  if(I->ReshapeFlag) {
    MainDoReshape(WinX,WinY);
    I->ReshapeFlag=false;
  }
  if(!PMGUI) {
    if(!OrthoCommandWaiting(TempPyMOLGlobals)) {
      if((!PyMOLOption->keep_thread_alive)&&(!PyMOLOption->read_stdin)) {
        I->IdleCount++;
        if(I->IdleCount==10) {
          PLockAPIAsGlut();
          PParse("_quit");
          PFlush();
          PUnlockAPIAsGlut();
        }
      }
    }
  }
  PRINTFD(TempPyMOLGlobals,FB_Main)
    " MainBusyIdle: leaving... IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;

}
void MainRepositionWindowDefault(void)
{
  p_glutPositionWindow(WinPX,WinPY);
  p_glutReshapeWindow(WinX,WinY);
}
/*========================================================================*/

#ifdef WIN32
BOOL WINAPI HandlerRoutine(
						     DWORD dwCtrlType   //  control signal type
)
{
	switch(dwCtrlType) {
	case CTRL_CLOSE_EVENT:
	case CTRL_BREAK_EVENT:
	case CTRL_C_EVENT:
		TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
	
		break;
	}
	return 1;
}
#endif

#ifdef _PYMOL_SHARP3D
void sharp3d_prepare_context(void);
#endif

void launch(void)
{
  int multisample_mask = 0;
  PyMOLInstance = ClassPyMOLNew();

  if(PyMOLOption->multisample)
    multisample_mask = P_GLUT_MULTISAMPLE;
  
  if(InternalGUI&&(!GameMode))
    WinX+=cOrthoRightSceneMargin;
  if(InternalFeedback&&(!GameMode))
    WinY+= (InternalFeedback-1)*cOrthoLineHeight + cOrthoBottomSceneMargin;

  if(PMGUI) {
    #ifndef _PYMOL_OSX
    atexit(MainOnExit); /* register callback to help prevent crashes
                                 when GLUT spontaneously kills us */
    #endif

#ifdef WIN32
SetConsoleCtrlHandler(
  HandlerRoutine,  // address of handler function
  true                          // handler to add or remove
);
 
#endif
#ifdef _PYMOL_SHARP3D
      sharp3d_prepare_context();
#endif

      p_glutInit(&myArgc, myArgv);

    switch(ForceStereo) {

    case -1: /* force mono */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask | P_GLUT_DOUBLE );
      StereoCapable = 0;
      break;

    case 0: /* default/autodetect (stereo on win/unix; mono on macs) */
#ifdef _PYMOL_SHARP3D
      WinX = 794+220;
      WinY = 547;
      // WinPX = 0;
      //WinPY = 0;

      glutInitDisplayMode( P_GLUT_RGBA| P_GLUT_DOUBLE| P_GLUT_STENCIL );
      StereoCapable = 1;
#else
#ifndef _PYMOL_OSX
      /* don't try stereo on OS X unless asked to do so */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask | P_GLUT_DOUBLE | P_GLUT_STEREO );
      if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
#endif
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | multisample_mask | P_GLUT_DOUBLE );
        if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
          p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );
        }
        StereoCapable = 0;
#ifndef _PYMOL_OSX
      } else {
        StereoCapable = 1;
      }
#endif
#endif
      break;

    case 1: /* force stereo (if possible) */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE | P_GLUT_STEREO );
      if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {

        StereoCapable = 0;
      } else {
        StereoCapable = 1;
      }      
      break;
    }

    if(!GameMode) {
      if((WinPX>-1000)&&(WinPY>-1000)) {
        #ifndef _PYMOL_FINK
        p_glutInitWindowPosition(WinPX,WinPY);
        #else
        p_glutInitWindowPosition(WinPX,WinPY-22); /* somethings wrong with FinkGlut's window positioning...*/
        #endif
      }
      p_glutInitWindowSize(WinX, WinY);
 
      TheWindow = p_glutCreateWindow("PyMOL Viewer");
      if(WindowIsVisible) {
        p_glutShowWindow();
      } else {
        p_glutHideWindow();
      }
        
    } else {
      char str[255];
      sprintf(str,"%dx%d:32@120",WinX,WinY);
      p_glutGameModeString(str);
      p_glutEnterGameMode(); 
    }
  } 

#ifdef _PYMOL_SHARP3D
  	// Setup OpenGL
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
#endif

   MainInit(ClassPyMOLGetGlobals(PyMOLInstance));
   PInit();
  
#ifdef _PYMOL_SHARP3D
  SettingSetGlobal_b(TempPyMOLGlobals,cSetting_overlay,1);
  SettingSetGlobal_b(TempPyMOLGlobals,cSetting_overlay_lines,1);
#endif

  if(PMGUI) {
    p_glutDisplayFunc(         MainDraw );
    p_glutReshapeFunc(         MainReshape );
    p_glutKeyboardFunc(        MainKey );
    p_glutMouseFunc(           MainButton );
    p_glutMotionFunc(          MainDrag );
    p_glutPassiveMotionFunc(   MainPassive );
    p_glutSpecialFunc(         MainSpecial );
    p_glutIdleFunc(            MainBusyIdle );

    p_glutPostRedisplay();
  }

  PUnblock();

  if(PMGUI) {
    SceneSetCardInfo(TempPyMOLGlobals,(char*)glGetString(GL_VENDOR),
                     (char*)glGetString(GL_RENDERER),
                     (char*)glGetString(GL_VERSION));
    if(ShowSplash) {
      
      printf(" OpenGL-based graphics engine:\n");
      printf("  GL_VENDOR: %s\n",(char*)glGetString(GL_VENDOR));
      printf("  GL_RENDERER: %s\n",(char*)glGetString(GL_RENDERER));
      printf("  GL_VERSION: %s\n",(char*)glGetString(GL_VERSION));
      /*        printf("  GL_EXTENSIONS: %s\n",(char*)glGetString(GL_EXTENSIONS));*/
      if(StereoCapable) {
        printf(" Hardware stereo capability detected.\n");
      } else if((ForceStereo==1)&&(!StereoCapable)) {
        printf(" Hardware stereo not present (unable to force).\n");
      }
    } 
    if(!WindowIsVisible)
      MainReshape(WinX,WinY);
    p_glutMainLoop();
    PBlock(); /* if we've gotten here, then we're heading back to Python... */
  } else {
    SceneSetCardInfo(TempPyMOLGlobals,"none","ray trace only","none");
    if(ShowSplash) 
      printf(" Command mode. No graphics front end.\n");
    MainReshape(WinX,WinY);
    MainDraw(); /* for command line processing */
    while(1) {
      MainBusyIdle();
      MainDraw();
    }
  }
}

int MainCheckRedundantOpen(char *file)
{
  int result = false;
  PBlock();
  result = PTruthCallStr(P_cmd,"check_redundant_open",file);
  PUnblock();
  return result;
}

/*========================================================================*/
#ifndef _PYMOL_MODULE
int main(int argc, char *argv[])
{
  myArgc=argc;
  myArgv=argv;

  PInitEmbedded(argc,argv);

#else
int was_main(void)
{
  myArgc=1;
  strcpy(myArgvvv,"pymol");
  myArgvv[0]=myArgvvv;
  myArgvv[1]=NULL;
  myArgv=myArgvv;

#ifdef _DRI_WORKAROUND
  dlopen("libGL.so.1",RTLD_LAZY|RTLD_GLOBAL);
#endif

#endif  

  PGetOptions(PyMOLOption);

  /* below need to be phased out by modifying code to use
     PyMOLOption global */

  PMGUI = PyMOLOption->pmgui;
  InternalGUI=PyMOLOption->internal_gui;
  ShowSplash = PyMOLOption->show_splash;
  InternalFeedback = PyMOLOption->internal_feedback;
  Security = PyMOLOption->security;
  GameMode = PyMOLOption->game_mode;
  ForceStereo = PyMOLOption->force_stereo;
  WinX = PyMOLOption->winX;
  WinY = PyMOLOption->winY;
  BlueLine = PyMOLOption->blue_line;
  WinPX = PyMOLOption->winPX;
  WinPY = PyMOLOption->winPY;
  ExternalGUI = PyMOLOption->external_gui;
  PyMOLMultisample = PyMOLOption->multisample; /* currently only used by MacOSX */
  WindowIsVisible = PyMOLOption->window_visible;

#ifdef _PYMOL_SHARP3D
  //  InternalGUI = 0;
  InternalFeedback = 0;
  ShowSplash = 0;
#endif

  launch();

  return 0;

}






